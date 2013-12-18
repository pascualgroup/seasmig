package seasmig.treelikelihood.trees;

import java.util.ArrayList;
import java.util.HashMap;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.TransitionModel.Transition;
import seasmig.util.Util;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


@SuppressWarnings("serial")
public class TreeWithLocationsAndNucleotides implements LikelihoodTree {

	// ENCOUDING FOR UNKNOWN TAXA NUMBER (i.e. internal nodes)
	public static final int UNKNOWN_TAXA = -1;

	// ENCOUDING FOR UNKNOWN GEOGRAPHIC LOCATION (i.e. internal nodes)
	public static final int UNKNOWN_LOCATION = -1;

	// ENCODING FOR LOCATION ERROR 
	public static final int ERR_LOCATION = -2;

	// CONST 
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	// Tree & Models
	TreeWithLocationsAndNucleotidesNode root = null;		
	private TransitionModel migrationLikelihoodModel = null;
	private TransitionModel[] codonLikelihoodModel = new TransitionModel[3]; // CP1, CP2, CP3

	// For Location Reading
	int numLocations = 0;
	private int numIdentifiedLocations;
	private int numIdentifiedSequences;

	// Taxa
	HashMap<String, Integer> taxaIndices = new HashMap<String,Integer>();

	double[] ZERO_LOG_PROBS;	

	int seqLength = 0;
	private double locationLogLike = 0;
	private double seqLogLike = 0;
	private double logLike = 0;


	// for test purpose
	public TreeWithLocationsAndNucleotides(TreeWithLocationsAndNucleotidesNode root, int numLocations, int seqLength) {
		taxaIndices = null;
		this.numLocations=numLocations;				

		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}

		// PARSE SEQUENCE
		this.seqLength = seqLength;

		// Create new node
		this.root = root;
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	// sequences are loaded from a hashmap
	public TreeWithLocationsAndNucleotides(jebl.evolution.trees.SimpleRootedTree tree,HashMap<String,Integer> taxaIndices_, HashMap<String, Integer> locationMap, HashMap<String, String> seqStrMap, int num_locations_, double lastTipTime) {		
		taxaIndices = taxaIndices_;
		numLocations=num_locations_;				

		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}

		// PARSE LOCATION
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;		

		// PARSE SEQUENCE
		seqLength = ((String) seqStrMap.values().toArray()[0]).length();

		String seqStr = seqStrMap.get(tree.getTaxon(tree.getRootNode()));
		Sequence sequence;
		if (seqStr==null) 
			sequence=new Sequence(seqLength);
		else {
			sequence=new Sequence("ABCD123",seqStr);
			numIdentifiedSequences+=1;
		}

		// PARSE TAXA
		Integer rootTaxonIndex = UNKNOWN_TAXA;		
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex= UNKNOWN_TAXA;

		// Create new node
		root = new TreeWithLocationsAndNucleotidesNode(sequence, location,rootTaxonIndex,0,null);
		makeSubTree(tree,locationMap,seqStrMap, root,tree.getRootNode());
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


	public double locLogLikelihood() {
		for (TreeWithLocationsAndNucleotidesNode node : root) { // Postorder 

			if (node.loc==TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION) {
				node.logProbsLOC =new double[numLocations]; // Internal node initialization
			}
			else {
				node.logProbsLOC= ZERO_LOG_PROBS.clone();
				node.logProbsLOC[node.loc]=0; // Tip node
			}
			
			if (node.children.size()!=0) { // this is an internal node		
				for (TreeWithLocationsAndNucleotidesNode child : node.children ) {
					// for now caching is done inside likelihood model...
					DoubleMatrix2D p;						
					// TODO: check if clause (here for numerics issues) 
					if (node.time!=child.time && node.loc==child.loc) 												
						p = migrationLikelihoodModel.transitionMatrix(node.time, child.time);
					else
						p = migrationLikelihoodModel.transitionMatrix(node.time, child.time+Util.minValue);

					for (int from = 0; from < numLocations; from++) {
						double[] alphas = new double[numLocations];						
						for (int to = 0; to < numLocations; to++) { // Integrate over all possible locations
							alphas[to]=(Math.log(p.get(from,to)) + child.logProbsLOC[to]);							
						}						
						node.logProbsLOC[from] += logSumExp(alphas); // Probability of internal node state based on children 
					}								
				}

			}
		}

		// Calculate root base frequency contribution... 
		DoubleMatrix1D rootFreq = migrationLikelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbsLOC[i] + Math.log(rootFreq.get(i));
		}	
		locationLogLike = logSumExp(alphas);
		return locationLogLike;		
	}

	public double seqLogLikelihood() {
		seqLogLike=0;
		for (TreeWithLocationsAndNucleotidesNode node : root) { // Postorder
			if (node.logProbsCP==null) {
				node.logProbsCP = new double[3][(seqLength-1)/3+1][4];

				for (int codonPos = 0; codonPos <3 ; codonPos++) {
					for (int loc=0;loc<(seqLength-1)/3+1;loc++) {
						if ((loc*3+codonPos) >= seqLength) continue;
						node.logProbsCP[codonPos][loc]=node.seq.get(loc*3+codonPos).clone();									
					}
				}
			}

			if (node.children.size()!=0) { // this is an internal node			
				for (TreeWithLocationsAndNucleotidesNode child : node.children ) {					
					for (int codonPos = 0; codonPos <3 ; codonPos++) {
						DoubleMatrix2D p = null;						
						// TODO: Util.minValue here for numerics issues, check this... 
						p = codonLikelihoodModel[codonPos].transitionMatrix(node.time, child.time+Util.minValue);

						for (int loc=0;loc<((seqLength-1)/3+1);loc++) {
							for (int from = 0; from < 4; from++) {													
								if ((loc*3+codonPos) >= seqLength) continue;							
								double[] alphas = new double[4];						
								for (int to = 0; to < 4; to++) { // Integrate over all possible nucleotides									
									alphas[to]=(Math.log(p.get(from,to)) + child.logProbsCP[codonPos][loc][to]);							
								}
								node.logProbsCP[codonPos][loc][from] += logSumExp(alphas);							
							}
						}
					}								
				}
			}
		}

		// Add root frequency		
		DoubleMatrix1D pi[] = new DoubleMatrix1D[3];

		for (int i=0;i<3;i++) {
			pi[i]=codonLikelihoodModel[i].rootfreq(root.time);
		}


		for (int codonPos = 0; codonPos <3 ; codonPos++) {
			for (int loc=0;loc<((seqLength-1)/3+1);loc++) {				
				if ((loc*3+codonPos) >= seqLength) continue;			
				double[] alphas = new double[4];	
				for (int i = 0; i < 4; i++) {				
					alphas[i]=root.logProbsCP[codonPos][loc][i] + Math.log(pi[codonPos].get(i));
				}
				seqLogLike += logSumExp(alphas);
			}			
		}			
		return seqLogLike;		
	}

	@Override
	public double logLikelihood() {
		locLogLikelihood();		
		seqLogLikelihood();
		logLike = locationLogLike+seqLogLike;
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
		copyTree.migrationLikelihoodModel=this.migrationLikelihoodModel;
		copyTree.numIdentifiedLocations=this.numIdentifiedLocations;
		copyTree.numLocations=this.numLocations;		
		copyTree.ZERO_LOG_PROBS=this.ZERO_LOG_PROBS;
		copyTree.codonLikelihoodModel=this.codonLikelihoodModel;
		copyTree.numIdentifiedSequences=this.numIdentifiedSequences;
		copyTree.root = new TreeWithLocationsAndNucleotidesNode(root.seq, root.loc,root.taxonIndex,root.time,null);
		copyTree.taxaIndices = taxaIndices;
		copyTree.seqLength = seqLength;
		treeCopy(this.root, copyTree.root);  
		return copyTree;
	}

	private void treeCopy(TreeWithLocationsAndNucleotidesNode from, TreeWithLocationsAndNucleotidesNode to) {
		for (TreeWithLocationsAndNucleotidesNode child : from.children) {
			TreeWithLocationsAndNucleotidesNode newChild = new TreeWithLocationsAndNucleotidesNode(child.seq, child.loc,child.taxonIndex,child.time, to);
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
			HashMap<String, Integer> locationMap, HashMap<String, String> seqStrMap, TreeWithLocationsAndNucleotidesNode root,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
			Integer taxonIndex = TreeWithLocationsAndNucleotides.UNKNOWN_TAXA;	
			Sequence sequence = null;
			if (taxon!=null) {

				// Parse location
				location = locationMap.get(taxon.toString());				
				if (location==null) 
					location=TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;

				// Parse sequence
				String seqStr = seqStrMap.get(taxon.toString());	
				if (seqStr==null) 
					sequence=new Sequence(seqLength);
				else {
					sequence=new Sequence("ABCD123",seqStr);
					numIdentifiedSequences+=1;
				}

				// Parse taxon
				taxonIndex = taxaIndices.get(taxon.toString());
				if (taxonIndex==null) 
					taxonIndex = TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION;

			}												
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;



			root.children.add(new TreeWithLocationsAndNucleotidesNode(sequence, location,taxonIndex,root.time+inputTree.getLength(node),root));			
			makeSubTree(inputTree,locationMap, seqStrMap, root.children.get(root.children.size()-1), node);			
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
	public void setMigrationModel(Object likelihoodModel_) {
		migrationLikelihoodModel = (TransitionModel) likelihoodModel_;
	}

	public String newickProbs() {	 
		return newickProbs(root,migrationLikelihoodModel.rootfreq(root.time).toArray()) + "\n";
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
		return locationLogLike;
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
	public String newickStochasticMapping(int maxBranchRetries) {
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
				event = migrationLikelihoodModel.nextEvent(currentTime, currentLoc);
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
		DoubleMatrix1D rootFreq = migrationLikelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbsLOC[i] + Math.log(rootFreq.get(i));
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
		DoubleMatrix1D p = migrationLikelihoodModel.probability(parent.loc, parent.time, node.time);
		for (int i=0; i < numLocations; i++) {								
			alphas[i] = cern.jet.math.Functions.log.apply(p.get(i)) + node.logProbsLOC[i];
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


	@Override
	public void setCodonModel(Object codonModel) {
		this.codonLikelihoodModel = (TransitionModel[]) codonModel;		
	}

	@Override
	public String smDescendants() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String smTrunkStats(double presentDayTipInterval,
			double timeToDesignateTrunk) {
		// TODO Auto-generated method stub
		return null;
	}	


}



