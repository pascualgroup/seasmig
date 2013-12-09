package seasmig.treelikelihood.trees;

import java.util.ArrayList;
import java.util.HashMap;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.MigrationBaseModel.Transition;
import seasmig.util.Util;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


@SuppressWarnings("serial")
public class TreeWithLocationsAndNucleotides implements LikelihoodTree {

	public static final int UNKNOWN_TAXA = -1;
	public static final int UNKNOWN_LOCATION = -1;
	public static final int ERR_LOCATION = -2;
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	// Tree & Model
	TreeWithLocationsAndNucleotidesNode root = null;		
	private MigrationBaseModel likelihoodModel = null;

	int numLocations = 0;
	private int numIdentifiedLocations;

	// Taxa
	HashMap<String, Integer> taxaIndices = new HashMap<String,Integer>();

	double[] ZERO_LOG_PROBS;
	private double logLike = 0;

	// Tree generate parameters for test purpose
	static final private double testBranchLengthMean = 0.5;
	static final private double testBranchLengthVariance = 1.0;

	// Load a tree from a basic jebl tree
	// locations are loaded from nexsus tree trait location_attribute name
	public TreeWithLocationsAndNucleotides(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String,Integer> taxaIndices_, String locationAttributeName, int num_locations_, double lastTipTime) {
		taxaIndices = taxaIndices_;
		numLocations=num_locations_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex = UNKNOWN_TAXA;
		if (tree.getRootNode().getAttribute(locationAttributeName)!=null) {
			root = new TreeWithLocationsAndNucleotidesNode((Integer) tree.getRootNode().getAttribute(locationAttributeName),rootTaxonIndex,0,null);
			numIdentifiedLocations+=1;
		}
		else {
			root = new TreeWithLocationsAndNucleotidesNode(UNKNOWN_TAXA,rootTaxonIndex,0,null);
		}

		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
		recalibrateTimes(root, lastTipTime);
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocationsAndNucleotides(jebl.evolution.trees.SimpleRootedTree tree,HashMap<String,Integer> taxaIndices_, HashMap<String, Integer> locationMap, int num_locations_, double lastTipTime/*, HashMap<String, Double> stateMap*/) {
		taxaIndices = taxaIndices_;
		numLocations=num_locations_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;		

		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex= UNKNOWN_TAXA;
		root = new TreeWithLocationsAndNucleotidesNode(location,rootTaxonIndex,0,null);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
		recalibrateTimes(root, lastTipTime);		
	}

	private void recalibrateTimes(TreeWithLocationsAndNucleotidesNode root, double lastTipTime) {		
		double maxTime=Double.NEGATIVE_INFINITY;
		for (TreeWithLocationsAndNucleotidesNode node : root) {
			if (node.time>maxTime) {
				maxTime=node.time;
			}
		}
		for (TreeWithLocationsAndNucleotidesNode node : root) {
			node.time = node.time - maxTime + lastTipTime; 
		}
	}

	public TreeWithLocationsAndNucleotides(TreeWithLocationsAndNucleotidesNode root_, MigrationBaseModel likelihoodModel_) {
		root = root_;
		likelihoodModel=likelihoodModel_;
		numLocations=likelihoodModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
	}

	@Override
	public double logLikelihood() {
		for (TreeWithLocationsAndNucleotidesNode node : root) { // Postorder 

			if (node.loc==TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION) {
				node.logProbs =new double[numLocations]; // Internal node initialization
			}
			else {
				node.logProbs= ZERO_LOG_PROBS.clone();
				node.logProbs[node.loc]=0; // Tip node
			}

			if (node.children.size()!=0) { // this is an internal node			
				for (int from = 0; from < numLocations; from++) {
					for (TreeWithLocationsAndNucleotidesNode child : node.children ) {
						// for now caching is done inside likelihood model...
						DoubleMatrix2D p;						
						// TODO: check if clause (here for numerics issues) 
						if (node.time!=child.time && node.loc==child.loc) 												
							p = likelihoodModel.transitionMatrix(node.time, child.time);
						else
							p = likelihoodModel.transitionMatrix(node.time, child.time+Util.minValue);
						double[] alphas = new double[numLocations];						
						for (int to = 0; to < numLocations; to++) { // Integrate over all possible locations
							alphas[to]=(Math.log(p.get(from,to)) + child.logProbs[to]);							
						}						
						node.logProbs[from] += logSumExp(alphas); // Probability of internal node state based on children 
					}								
				}

			}


		}

		// Calculate root base frequency contribution... 
		DoubleMatrix1D rootFreq = likelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbs[i] + Math.log(rootFreq.get(i));
		}	
		logLike = logSumExp(alphas);
		return logLike;		
	}

	private void removeInternalLocations(TreeWithLocationsAndNucleotidesNode node) {
		if (node.children.size()!=0) {
			node.loc=TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
			for (TreeWithLocationsAndNucleotidesNode child : node.children) {
				removeInternalLocations(child);				
			}
		}				
	}

	@Override
	public LikelihoodTree copy() {
		// TODO: test this... and or remove ...
		TreeWithLocationsAndNucleotides copyTree = new TreeWithLocationsAndNucleotides();
		copyTree.likelihoodModel=this.likelihoodModel;
		copyTree.numIdentifiedLocations=this.numIdentifiedLocations;
		copyTree.numLocations=this.numLocations;		
		copyTree.ZERO_LOG_PROBS=this.ZERO_LOG_PROBS;
		copyTree.root = new TreeWithLocationsAndNucleotidesNode(root.loc,root.taxonIndex,root.time,null);
		copyTree.taxaIndices = taxaIndices;
		treeCopy(this.root, copyTree.root);  
		return copyTree;
	}

	private void treeCopy(TreeWithLocationsAndNucleotidesNode from, TreeWithLocationsAndNucleotidesNode to) {
		for (TreeWithLocationsAndNucleotidesNode child : from.children) {
			TreeWithLocationsAndNucleotidesNode newChild = new TreeWithLocationsAndNucleotidesNode(child.loc,child.taxonIndex,child.time, to);
			to.children.add(newChild);			
			treeCopy(child, newChild);
		}		
	}


	@Override
	public String print() {
		return print(root);
	}

	public String print(TreeWithLocationsAndNucleotidesNode treePart) {
		String returnValue = Integer.toString(treePart.loc);
		if (treePart.children.size()>0) {
			returnValue+=" (";
			returnValue+=print(treePart.children.get(0));
			if (treePart.children.size()>1) {
				for (int i=1;i<treePart.children.size();i++) {
					returnValue+=",";
					returnValue+=treePart.children.get(i);					
				}				
			}
			returnValue+=")";
		}
		return returnValue;
	}


	public void removeInternalLocations() {
		removeInternalLocations(root);

	}
	protected TreeWithLocationsAndNucleotides() {
	}

	private void makeSubTree(SimpleRootedTree inputTree,
			HashMap<String, Integer> locationMap, TreeWithLocationsAndNucleotidesNode root,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
			Integer taxonIndex = TreeWithLocationsAndNucleotides.UNKNOWN_TAXA;			
			if (taxon!=null) {
				location = locationMap.get(taxon.toString());				
				if (location==null) 
					location=TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;
				taxonIndex = taxaIndices.get(taxon.toString());
				if (taxonIndex==null) 
					taxonIndex = TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;

			}			
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;
			root.children.add(new TreeWithLocationsAndNucleotidesNode(location,taxonIndex,root.time+inputTree.getLength(node),root));			
			makeSubTree(inputTree,locationMap, root.children.get(root.children.size()-1), node);			
		}
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, String locationAttributeName, TreeWithLocationsAndNucleotidesNode outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {
			Taxon taxon = inputTree.getTaxon(node);
			Integer taxonIndex = UNKNOWN_TAXA;
			if (taxon!=null)
				taxonIndex = taxaIndices.get(taxon.getName());			
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;
			if (locationAttributeName!=null) {
				if (node.getAttribute(locationAttributeName)!=null) {
					outputSubTree.children.add(new TreeWithLocationsAndNucleotidesNode((Integer)node.getAttribute(locationAttributeName),taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
					numIdentifiedLocations+=1;
				}
				else {
					outputSubTree.children.add(new TreeWithLocationsAndNucleotidesNode(TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
				}
			}
			else 
				outputSubTree.children.add(new TreeWithLocationsAndNucleotidesNode(TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, TreeWithLocationsAndNucleotidesNode root, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.loc, location, root.time, to_time));
					if (cern.jet.random.Uniform.staticNextDouble()<=p) {
						root.children.add(new TreeWithLocationsAndNucleotidesNode(location,TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,to_time,root));
						break;
					}
				}			
			}

			for (TreeWithLocationsAndNucleotidesNode child : root.children) {
				makeRandomTree(m,child,(int) Math.round(nNodes/2.0));
			}
		}

	}

	@Override
	public int getNumLocations() {		
		return numLocations;
	}

	public int getNumIdentifiedLocations() {
		// TODO: this for trees with states...
		return numIdentifiedLocations;
	}

	@Override 
	public void setLikelihoodModel(Object likelihoodModel_) {
		likelihoodModel = (MigrationBaseModel) likelihoodModel_;
	}

	public String newickProbs() {	 
		return newickProbs(root,likelihoodModel.rootfreq(root.time).toArray()) + "\n";
	}

	private String newickProbs(TreeWithLocationsAndNucleotidesNode treePart, double[] rootFreq) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.parent.time);
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+treePart.parseProbs(rootFreq)+":"+branchLength);
		}
		else if (treePart.children.size()>0) {
			returnValue+="(";
			returnValue+=newickProbs(treePart.children.get(0),rootFreq);
			for (int i = 1; i < treePart.children.size(); i++){
				returnValue+=",";
				returnValue+=newickProbs(treePart.children.get(i),rootFreq);	
			}
			returnValue+=")";
			double parentTime=0;
			if (treePart.parent!=null) {
				parentTime=treePart.parent.time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+=treePart.parseProbs(rootFreq)+":"+branchLength;
		}		
		return returnValue;
	}

	@Override
	public double cachedLogLikelihood() {		
		return logLike;
	}

	public final double logSumExp(double[] alphas) {
		// TODO: improve this
		double sumExp = 0;
		double minWithoutNegInf = 0;
		for (int i=0;i<alphas.length;i++) {
			if (!Double.isInfinite(alphas[i])) {
				if (alphas[i]<minWithoutNegInf) {
					minWithoutNegInf = alphas[i];
				}
			}
		}
		for (int i=0;i<alphas.length;i++) {			
			sumExp=sumExp+cern.jet.math.Functions.exp.apply(alphas[i]-minWithoutNegInf);
		}
		double returnValue=minWithoutNegInf+cern.jet.math.Functions.log.apply(sumExp);
		if (Double.isNaN(returnValue) ){
			return minNegative;
		}
		return returnValue;
	}

	@Override
	public String newickStochasticMapping() {
		asr(); // Ancestral state reconstruction
		stochsticMapping(root);

		return newickSM(root);
	}

	private String newickSM(TreeWithLocationsAndNucleotidesNode treePart) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.parent.time);
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+ "[&states="+Integer.toString(treePart.loc)+"]:"+treePart.parseMap()+branchLength);
		}
		else {
			returnValue+="(";
			returnValue+=newickSM(treePart.children.get(0));
			for (int i = 1; i < treePart.children.size(); i++){
				returnValue+=",";
				returnValue+=newickSM(treePart.children.get(i));	
			}
			returnValue+=")";
			double parentTime=0;
			if (treePart.parent!=null) {
				parentTime=treePart.parent.time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+="[&states="+Integer.toString(treePart.loc)+"]:"+treePart.parseMap()+branchLength;
		}		
		return returnValue;
	}

	private void stochsticMapping(TreeWithLocationsAndNucleotidesNode root) {
		// TODO: test
		// TODO: cite
		// TODO: preorder iterator			
		for (TreeWithLocationsAndNucleotidesNode child : root.children) { 
			int currentLoc = root.loc;
			double currentTime = root.time;
			boolean doneWithBranch = false;
			Transition event = null;
			do {
				event = likelihoodModel.nextEvent(currentTime, currentLoc);
				if (event.time < child.time) {
					if (child.transitions==null) child.transitions = new ArrayList<Transition>();					
					child.transitions.add(event);
					currentLoc = event.loc;
					currentTime = event.time;
				}
				else if (currentLoc!=child.loc) {
					// If there is a mismatch between stochastic mapping and child ASR than we 
					// restart the entire branch
					if (child.transitions!=null) {
						child.transitions.clear();
					}
					currentLoc = root.loc;
					currentTime = root.time;
				} else {
					doneWithBranch = true;								
				}
			} while (!doneWithBranch);

			stochsticMapping(child);
		}		
	}

	private int normalizeAndGetRandomSampleFromLogProbs(double[] logProbs) {
		double min=logProbs[0];		
		for (int i=1;i<logProbs.length;i++) {
			if (Double.isInfinite(min)) 
				min = logProbs[i];			
			if (min>logProbs[i] &&  !Double.isInfinite(logProbs[i]))
				min=logProbs[i];
		}
		double sum=0;			
		for (int i=0;i<logProbs.length;i++) {
			sum+=cern.jet.math.Functions.exp.apply(logProbs[i]-min);
		}

		double p=0;	
		for (int i=0;i<logProbs.length;i++) {
			p=p+cern.jet.math.Functions.exp.apply(logProbs[i]-min);
			if (cern.jet.random.Uniform.staticNextDouble()<=(p/sum)) {
				return i;				
			}			
		}
		return TreeWithLocationsAndNucleotides.ERR_LOCATION;
	}

	private void asr() {

		// Calculate root state
		// TODO: check this
		DoubleMatrix1D rootFreq = likelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbs[i] + Math.log(rootFreq.get(i));
		}		
		root.loc = normalizeAndGetRandomSampleFromLogProbs(alphas);		
		for (TreeWithLocationsAndNucleotidesNode node : root.children) {
			asr(node);
		}		
	}

	private void asr(TreeWithLocationsAndNucleotidesNode node) {
		// TODO: check this
		TreeWithLocationsAndNucleotidesNode parent = node.parent;		
		double[] alphas = new double[numLocations];	
		// TODO: check if clause (here for numerics issues)
		DoubleMatrix1D p = likelihoodModel.probability(parent.loc, parent.time, node.time);
		for (int i=0; i < numLocations; i++) {								
			alphas[i] = cern.jet.math.Functions.log.apply(p.get(i)) + node.logProbs[i];
		}		
		node.loc = normalizeAndGetRandomSampleFromLogProbs(alphas);

		for (TreeWithLocationsAndNucleotidesNode child : node.children) {
			asr(child);
		}
	}	

	@Override
	public String newickAncestralStateReconstruction() {
		// TODO: Check this
		asr();
		String returnValue = newickStates(root) + "\n";
		return returnValue;

	}

	private String newickStates(TreeWithLocationsAndNucleotidesNode treePart) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.parent.time);			
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+ "[&states="+Integer.toString(treePart.loc)+"]:"+branchLength);
		}
		else {
			returnValue+="(";
			returnValue+=newickStates(treePart.children.get(0));
			for (int i = 1; i < treePart.children.size(); i++){
				returnValue+=",";
				returnValue+=newickStates(treePart.children.get(i));	
			}
			returnValue+=")";
			double parentTime=0;
			if (treePart.parent!=null) {
				parentTime=treePart.parent.time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+="[&states="+Integer.toString(treePart.loc)+"]:"+branchLength;
		}		
		return returnValue;
	}

	@Override
	public String smTransitions() {

		String returnValue = "{";
		String[][] transitionTimes = new String[numLocations][numLocations];		

		for (TreeWithLocationsAndNucleotidesNode node : root) {
			if (node==root) continue;
			if (node.transitions!=null) {
				int fromLocation = node.parent.loc;
				for (Transition transition : node.transitions) {
					if (transitionTimes[fromLocation][transition.loc]==null) {
						transitionTimes[fromLocation][transition.loc]=new String();
					}
					if (fromLocation!=transition.loc) 
						transitionTimes[fromLocation][transition.loc]+=String.format("%.3f,",transition.time);
					fromLocation=transition.loc;
				}				
			}
		}

		for (int i=0;i<numLocations;i++) {
			returnValue+="{";
			for (int j=0;j<numLocations;j++) {						
				if (transitionTimes[i][j]!=null) {
					returnValue+="{"+transitionTimes[i][j].substring(0, transitionTimes[i][j].length()-1)+"}";
				}
				else {
					returnValue+="{}";
				}
				if (j!=(numLocations-1)) {
					returnValue+=",";
				}				
			}
			returnValue+="}";
			if (i!=(numLocations-1)) {
				returnValue+=",";
			}
		}
		returnValue+="}";
		return returnValue;
	}

	@Override
	public String smTipDwellings() {
		String returnValue = "{";
		String[] tipDwellings = new String[numLocations];	

		for (int i=0;i<numLocations;i++) {	
			tipDwellings[i]=new String();			
		}

		for (TreeWithLocationsAndNucleotidesNode node : root) {			
			if (node==root) continue;
			if (!node.isTip()) continue;
			int fromLocation = node.loc;
			double fromTime = tracebackWithNoChanges(node);
			tipDwellings[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
		}

		for (int i=0;i<numLocations;i++) {
			if (tipDwellings[i].length()>0) {
				returnValue+="{"+tipDwellings[i].substring(0, tipDwellings[i].length()-1)+"}";
			}
			else {
				returnValue+="{}";
			}
			if (i!=(numLocations-1)) {
				returnValue+=",";
			}			
		}
		returnValue+="}";
		return returnValue;
	}

	private double tracebackWithNoChanges(TreeWithLocationsAndNucleotidesNode node) {
		// TODO: test this
		if (node.parent==null) {
			return node.time;
		}
		if (node.transitions!=null) {
			if (node.transitions.size()>0) {
				return (node.transitions.get(node.transitions.size()-1).time);
			}
			else {
				return (tracebackWithNoChanges(node.parent));
			}
		}
		else {
			return (tracebackWithNoChanges(node.parent));
		}
	}

	@Override
	public String smLineages() {
		// TODO: test this
		String returnValue = "{";
		String[] lineages = new String[numLocations];	
		
		for (int i=0;i<numLocations;i++) {	
			lineages[i]=new String();			
		}

		for (TreeWithLocationsAndNucleotidesNode node : root) {
			// Only record lineages from node.parent.time to node.time  
			if (node==root) continue;
			int fromLocation = node.loc;
			double fromTime = node.parent.time;
			if (node.transitions==null) {
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}
			else if (node.transitions.size()==0) {
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}
			else {	
				for (Transition transition : node.transitions) { 														
					lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,transition.time);
					fromLocation=transition.loc;
					fromTime=transition.time;
				}
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}			
		}

		for (int i=0;i<numLocations;i++) {	
			if (lineages[i]!=null) {
				returnValue+="{"+lineages[i].substring(0, lineages[i].length()-1)+"}";
			} else {
				returnValue+="{}";
			}				
			if (i!=(numLocations-1)) {
				returnValue+=",";
			}			
		}
		returnValue+="}";
		return returnValue;

	}


}



