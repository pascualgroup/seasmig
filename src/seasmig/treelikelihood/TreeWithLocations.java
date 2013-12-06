package seasmig.treelikelihood;

import java.util.ArrayList;
import java.util.HashMap;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.treelikelihood.MigrationBaseModel.Event;
import seasmig.util.Util;


@SuppressWarnings("serial")
public class TreeWithLocations implements LikelihoodTree {

	public static final int UNKNOWN_TAXA = -1;
	public static final int UNKNOWN_LOCATION = -1;
	public static final int ERR_LOCATION = -2;
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	// Tree & Model
	TreeWithLocationsNode root = null;		
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

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocations(MigrationBaseModel createTreeModel, int numNodes) {		
		numLocations=createTreeModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		int rootLocation = getRandomSampleFrom(createTreeModel.rootfreq(0));
		root = new TreeWithLocationsNode(rootLocation,TreeWithLocations.UNKNOWN_TAXA,0,null);
		makeRandomTree(createTreeModel, root, numNodes);		
	}

	// Generate random tree states based on input tree topology and model .... 
	public TreeWithLocations(MigrationBaseModel createTreeModel, jebl.evolution.trees.SimpleRootedTree tree) {
		numLocations=createTreeModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		likelihoodModel=createTreeModel;
		int rootLocation = getRandomSampleFrom(likelihoodModel.rootfreq(0));
		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex = UNKNOWN_TAXA;
		root = new TreeWithLocationsNode(rootLocation,rootTaxonIndex,0,null);
		makeSubTree(tree,(String)null, root,tree.getRootNode());		
	}

	private void fillRandomTraits(TreeWithLocationsNode root) {
		if (root.children!=null) {
			for (TreeWithLocationsNode child : root.children) {	
				double p=0;
				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(likelihoodModel.logprobability(root.loc, location, root.time, child.time));
					if (cern.jet.random.Uniform.staticNextDouble()<=p) {
						child.loc=location;
						break;
					}
				}			
			}
			for (TreeWithLocationsNode child : root.children) {
				fillRandomTraits(child);
			}
		}
	}

	public void fillRandomTraits() {
		fillRandomTraits(root);
	}


	// Load a tree from a basic jebl tree
	// locations are loaded from nexsus tree trait location_attribute name
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String,Integer> taxaIndices_, String locationAttributeName, int num_locations_, double lastTipTime) {
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
			root = new TreeWithLocationsNode((Integer) tree.getRootNode().getAttribute(locationAttributeName),rootTaxonIndex,0,null);
			numIdentifiedLocations+=1;
		}
		else {
			root = new TreeWithLocationsNode(UNKNOWN_TAXA,rootTaxonIndex,0,null);
		}

		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
		recalibrateTimes(root, lastTipTime);
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree,HashMap<String,Integer> taxaIndices_, HashMap<String, Integer> locationMap, int num_locations_, double lastTipTime/*, HashMap<String, Double> stateMap*/) {
		taxaIndices = taxaIndices_;
		numLocations=num_locations_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=TreeWithLocations.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;		

		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex= UNKNOWN_TAXA;
		root = new TreeWithLocationsNode(location,rootTaxonIndex,0,null);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
		recalibrateTimes(root, lastTipTime);		
	}

	private void recalibrateTimes(TreeWithLocationsNode root, double lastTipTime) {		
		double maxTime=Double.NEGATIVE_INFINITY;
		for (TreeWithLocationsNode node : root) {
			if (node.time>maxTime) {
				maxTime=node.time;
			}
		}
		for (TreeWithLocationsNode node : root) {
			node.time = node.time - maxTime + lastTipTime; 
		}
	}

	public TreeWithLocations(TreeWithLocationsNode root_, MigrationBaseModel likelihoodModel_) {
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
		for (TreeWithLocationsNode node : root) { // Postorder 

			if (node.loc==TreeWithLocations.UNKNOWN_LOCATION) {
				node.logProbs =new double[numLocations]; // Internal node initialization
			}
			else {
				node.logProbs= ZERO_LOG_PROBS.clone();
				node.logProbs[node.loc]=0; // Tip node
			}

			if (node.children.size()!=0) { // this is an internal node			
				for (int from = 0; from < numLocations; from++) {
					for (TreeWithLocationsNode child : node.children ) {
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

	private void removeInternalLocations(TreeWithLocationsNode node) {
		if (node.children.size()!=0) {
			node.loc=TreeWithLocations.UNKNOWN_LOCATION;
			for (TreeWithLocationsNode child : node.children) {
				removeInternalLocations(child);				
			}
		}				
	}

	@Override
	public LikelihoodTree copy() {
		// TOOD: test this...
		TreeWithLocations copyTree = new TreeWithLocations();
		copyTree.likelihoodModel=this.likelihoodModel;
		copyTree.numIdentifiedLocations=this.numIdentifiedLocations;
		copyTree.numLocations=this.numLocations;		
		copyTree.ZERO_LOG_PROBS=this.ZERO_LOG_PROBS;
		copyTree.root = new TreeWithLocationsNode(root.loc,root.taxonIndex,root.time,null);
		copyTree.taxaIndices = taxaIndices;
		treeCopy(this.root, copyTree.root);  
		return copyTree;
	}

	private void treeCopy(TreeWithLocationsNode from, TreeWithLocationsNode to) {
		for (TreeWithLocationsNode child : from.children) {
			TreeWithLocationsNode newChild = new TreeWithLocationsNode(child.loc,child.taxonIndex,child.time, to);
			to.children.add(newChild);			
			treeCopy(child, newChild);
		}		
	}


	@Override
	public String print() {
		return print(root);
	}

	public String print(TreeWithLocationsNode treePart) {
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
	protected TreeWithLocations() {
	}

	private void makeSubTree(SimpleRootedTree inputTree,
			HashMap<String, Integer> locationMap, TreeWithLocationsNode root,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = TreeWithLocations.UNKNOWN_LOCATION;
			Integer taxonIndex = TreeWithLocations.UNKNOWN_TAXA;			
			if (taxon!=null) {
				location = locationMap.get(taxon.toString());				
				if (location==null) 
					location=TreeWithLocations.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;
				taxonIndex = taxaIndices.get(taxon.toString());
				if (taxonIndex==null) 
					taxonIndex = TreeWithLocations.UNKNOWN_LOCATION;

			}			
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;
			root.children.add(new TreeWithLocationsNode(location,taxonIndex,root.time+inputTree.getLength(node),root));			
			makeSubTree(inputTree,locationMap, root.children.get(root.children.size()-1), node);			
		}
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, String locationAttributeName, TreeWithLocationsNode outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {
			Taxon taxon = inputTree.getTaxon(node);
			Integer taxonIndex = UNKNOWN_TAXA;
			if (taxon!=null)
				taxonIndex = taxaIndices.get(taxon.getName());			
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;
			if (locationAttributeName!=null) {
				if (node.getAttribute(locationAttributeName)!=null) {
					outputSubTree.children.add(new TreeWithLocationsNode((Integer)node.getAttribute(locationAttributeName),taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
					numIdentifiedLocations+=1;
				}
				else {
					outputSubTree.children.add(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
				}
			}
			else 
				outputSubTree.children.add(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, TreeWithLocationsNode root, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.loc, location, root.time, to_time));
					if (cern.jet.random.Uniform.staticNextDouble()<=p) {
						root.children.add(new TreeWithLocationsNode(location,TreeWithLocations.UNKNOWN_LOCATION,to_time,root));
						break;
					}
				}			
			}

			for (TreeWithLocationsNode child : root.children) {
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

	// TODO: (348[&antigenic={-6.00510611736,5.84199000915},rate=1.1478703001047978,states="japan_korea"]:2.44, ....
	private String newickProbs(TreeWithLocationsNode treePart, double[] rootFreq) {
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

	
	
	private String newickSM(TreeWithLocationsNode treePart) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.parent.time);
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+"[&states="+Integer.toString(treePart.loc)+"]:"+treePart.parseMap()+branchLength);
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
			returnValue+=treePart.loc+":"+branchLength;
		}		
		return returnValue;
	}

	private void stochsticMapping(TreeWithLocationsNode root) {
		// TODO Auto-generated method stub
		// TODO: this
		// TODO: cite
        // TODO: preorder iterator			
		for (TreeWithLocationsNode child : root.children) { 
			int currentLoc = root.loc;
			double currentTime = root.time;
			boolean doneWithBranch = false;
			Event event = null;
			do {
				event = likelihoodModel.nextEvent(currentTime, currentLoc);
				if (event.time < child.time) {
					if (child.changes==null) child.changes = new ArrayList<Event>();					
					child.changes.add(event);
					currentLoc = event.loc;
					currentTime = event.time;
				}
				else if (currentLoc!=child.loc) {
				// If there is a mismatch between stochastic mapping and child ASR than we 
				// restart the entire branch
					if (child.changes!=null) {
						child.changes.clear();
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

	private int getRandomSampleFrom(DoubleMatrix1D doubleMatrix1D) {
		double p=0;		
		for (int i=0;i<doubleMatrix1D.size();i++) {
			p=p+doubleMatrix1D.get(i);
			if (cern.jet.random.Uniform.staticNextDouble()<=p) {
				return i;				
			}			
		}
		return -1;
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
		return TreeWithLocations.ERR_LOCATION;
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
		for (TreeWithLocationsNode node : root.children) {
			asr(node);
		}		
	}

	private void asr(TreeWithLocationsNode node) {
		// TODO: check this
		TreeWithLocationsNode parent = node.parent;		
		double[] alphas = new double[numLocations];	
		// TODO: check if clause (here for numerics issues)
		DoubleMatrix1D p = likelihoodModel.probability(parent.loc, parent.time, node.time);
		for (int i=0; i < numLocations; i++) {								
			alphas[i] = cern.jet.math.Functions.log.apply(p.get(i)) + node.logProbs[i];
		}		
		node.loc = normalizeAndGetRandomSampleFromLogProbs(alphas);

		for (TreeWithLocationsNode child : node.children) {
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

	private String newickStates(TreeWithLocationsNode treePart) {
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


}



