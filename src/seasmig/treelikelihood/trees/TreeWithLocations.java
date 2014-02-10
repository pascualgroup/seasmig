package seasmig.treelikelihood.trees;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Stack;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.migrationmain.Config;
import seasmig.migrationmain.Config.SeqModelType;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.TransitionModel.Transition;
import seasmig.util.AltTreeOutput;
import seasmig.util.Util;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


@SuppressWarnings("serial")
public class TreeWithLocations implements LikelihoodTree {
	// TODO: check -Infinity likelihood source 

	// ENCOUDING FOR UNKNOWN TAXA NUMBER (i.e. internal nodes)
	public static final int UNKNOWN_TAXA = -1;

	// ENCOUDING FOR UNKNOWN GEOGRAPHIC LOCATION (i.e. internal nodes)
	public static final int UNKNOWN_LOCATION = -1;

	// ENCODING FOR LOCATION ERROR // TODO: replace with error 
	public static final int ERR_LOCATION = -2;

	// CONST 
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	// Tree & Models
	TreeWithLocationsNode root = null;		
	private TransitionModel migrationModel = null;
	private TransitionModel[] codonLikelihoodModel = new TransitionModel[3]; // CP1, CP2, CP3

	int numLocations = 0; // number of location categories

	// For Location Parsing
	private int numIdentifiedLocations;

	// For Location Parsing
	private int numIdentifiedSeqs;		

	// Taxa
	HashMap<String, Integer> taxaIndices = new HashMap<String,Integer>();

	double[] ZERO_LOG_PROBS ;	
	private Double logLike = null;

	// for sequences
	private int seqLength;
	private HashMap<String, Sequence> seqMap;

	// for likelihood
	private double seqLogLike;
	private double locationLogLike;

	private Config config = null;

	private boolean stochasticallyMapped = false;

	private boolean asrDone =false;
	private boolean asrSeqDone =false;
	private boolean sortingDone = false;

	static final Comparator<TreeWithLocationsNode> descendantOrder = new Comparator<TreeWithLocationsNode>() {
		public int compare(TreeWithLocationsNode v1, TreeWithLocationsNode v2) {
			Integer descendantsV1 = new Integer(getNumberOfDescendants(v1));
			Integer descendantsV2 = new Integer(getNumberOfDescendants(v2));
			return descendantsV1.compareTo(descendantsV2);
		}
	};	

	// for test purpose
	public TreeWithLocations(TreeWithLocationsNode root, int numLocations, int seqLength, Config config) {
		this.config = config;
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

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocations(TransitionModel createTreeModel, int numNodes, Config config) {
		this.config = config;
		numLocations=createTreeModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		int rootLocation = getRandomSampleFrom(createTreeModel.rootfreq(0));
		if (config.seqModelType==SeqModelType.NONE)
			root = new TreeWithLocationsNode(null, rootLocation,TreeWithLocations.UNKNOWN_TAXA,0,null,false);
		else
			root = new TreeWithLocationsNode(null, rootLocation,TreeWithLocations.UNKNOWN_TAXA,0,null,true);
		makeRandomTree(createTreeModel, root, numNodes);	
	}

	// Generate random tree states based on input tree topology and model .... 
	public TreeWithLocations(TransitionModel createTreeModel, jebl.evolution.trees.SimpleRootedTree tree, Config config) {
		this.config = config;
		numLocations=createTreeModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations]; // TODO: ???
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		migrationModel=createTreeModel;
		int rootLocation = getRandomSampleFrom(migrationModel.rootfreq(0));
		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex = UNKNOWN_TAXA;
		if (config.seqModelType==SeqModelType.NONE)
			root = new TreeWithLocationsNode(null, rootLocation,rootTaxonIndex,0,null,false);
		else
			root = new TreeWithLocationsNode(null, rootLocation,rootTaxonIndex,0,null,true);
		makeSubTree(tree,(String)null, root,tree.getRootNode());		
	}

	private void fillRandomTraits(TreeWithLocationsNode root) {
		if (root.children!=null) {
			for (TreeWithLocationsNode child : root.children) {	
				double p=0;
				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(migrationModel.logprobability(root.getLoc(), location, root.time, child.time));
					if (cern.jet.random.Uniform.staticNextDouble()<=p) {
						child.setLoc(location);
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
	// locations are loaded from a hashmap	
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree,HashMap<String,Integer> taxaIndices_, HashMap<String, Integer> locationMap, int num_locations_, double lastTipTime, HashMap<String, Sequence> seqMap_, int seqLength_, Config config) {
		this.config = config;
		taxaIndices = taxaIndices_;
		numLocations=num_locations_;
		seqMap=seqMap_;
		seqLength=seqLength_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=TreeWithLocations.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;		

		Sequence seq = null;
		if (seqMap!=null) {
			seq = seqMap.get(tree.getTaxon(tree.getRootNode()));
			if (seq==null) 
				seq=new Sequence(seqLength);
			else
				numIdentifiedSeqs+=1;
		}

		Integer rootTaxonIndex = UNKNOWN_TAXA;
		Taxon rootTaxon = tree.getTaxon(tree.getRootNode());
		if (rootTaxon!=null) {
			rootTaxonIndex = taxaIndices.get(rootTaxon.getName());
		}
		if (rootTaxonIndex==null)
			rootTaxonIndex= UNKNOWN_TAXA;

		if (config.seqModelType==SeqModelType.NONE)
			root = new TreeWithLocationsNode(seq,location,rootTaxonIndex,0,null,false);
		else
			root = new TreeWithLocationsNode(seq,location,rootTaxonIndex,0,null,true);

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

	public TreeWithLocations(TreeWithLocationsNode root_, TransitionModel likelihoodModel_, Config config) {
		this.config = config;
		root = root_;
		migrationModel=likelihoodModel_;
		numLocations=migrationModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
	}

	public double locLogLikelihood() {
		for (TreeWithLocationsNode node : root) { // Postorder 

			if (node.getLoc()==TreeWithLocations.UNKNOWN_LOCATION) {
				node.logProbsLOC =new double[numLocations]; // Internal node initialization
			}
			else {
				node.logProbsLOC= ZERO_LOG_PROBS.clone();
				node.logProbsLOC[node.getLoc()]=0; // Tip node
			}

			if (node.children.size()!=0) { // this is an internal node		
				for (TreeWithLocationsNode child : node.children ) {
					// for now caching is done inside likelihood model...
					DoubleMatrix2D p;						
					//					// TODO: check if clause (here for numerics issues) 
					p = migrationModel.transitionMatrix(node.time, child.time);

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
		DoubleMatrix1D rootFreq = migrationModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbsLOC[i] + Math.log(rootFreq.get(i));
		}	
		locationLogLike = logSumExp(alphas);
		return locationLogLike;		
	}

	public double seqLogLikelihood() {
		seqLogLike=0;
		if (seqLength>0) {
			for (TreeWithLocationsNode node : root) { // Postorder

				if (node.children.size()!=0) { // this is an internal node		
					if (node.logProbsCP==null) {
						node.logProbsCP = new double[3][(seqLength-1)/3+1][4];
					}					
				}
				else { // this is a tip
					node.logProbsCP = new double[3][(seqLength-1)/3+1][];
					for (int codonPos = 0; codonPos<3 ; codonPos++) {
						for (int loc=0;loc<(seqLength-1)/3+1;loc++) {
							if ((loc*3+codonPos) >= seqLength) continue;					
							node.logProbsCP[codonPos][loc]=node.seq.get(loc*3+codonPos);
						}
					}		
				}

				if (node.children.size()!=0) { // this is an internal node			
					for (TreeWithLocationsNode child : node.children ) {					
						for (int codonPos = 0; codonPos <3 ; codonPos++) {
							DoubleMatrix2D p = null;						
							// TODO: Util.minValue here for numerics issues, check this... 
							p = codonLikelihoodModel[codonPos].transitionMatrix(node.time, child.time);

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
		}
		return seqLogLike;		
	}

	public Iterable<TreeWithLocationsNode> eachPreorder() {
		ArrayList<TreeWithLocationsNode> returnValue = new ArrayList<TreeWithLocationsNode>();
		if(root==null) return returnValue;
		TreeWithLocationsNode current = root;
		Stack<TreeWithLocationsNode> stack=new Stack<TreeWithLocationsNode>();

		while(true){
			while(current!=null){
				//process current Node
				returnValue.add(current);
				stack.push(current);
				current=current.left();
			}
			if(stack.isEmpty()) break;
			current=stack.pop();
			current=current.right();
		}
		return returnValue;
	}

	@Override
	public double logLikelihood() {
		if (logLike==null) {
			logLike = locLogLikelihood()+seqLogLikelihood();
		}
		return logLike;
	}	

	@Override
	public LikelihoodTree copy() {
		TreeWithLocations copyTree = new TreeWithLocations();
		copyTree.migrationModel=this.migrationModel;
		copyTree.codonLikelihoodModel=this.codonLikelihoodModel;
		copyTree.numIdentifiedLocations=this.numIdentifiedLocations;
		copyTree.seqLength = this.seqLength;		
		copyTree.numIdentifiedSeqs=this.numIdentifiedSeqs;
		copyTree.numLocations=this.numLocations;		
		copyTree.ZERO_LOG_PROBS=this.ZERO_LOG_PROBS;
		copyTree.logLike=null;
		if (config.seqModelType==SeqModelType.NONE)
			copyTree.root = new TreeWithLocationsNode(root.seq, root.getLoc(),root.taxonIndex,root.time,null,false);
		else
			copyTree.root = new TreeWithLocationsNode(root.seq, root.getLoc(),root.taxonIndex,root.time,null,true);
		copyTree.taxaIndices = taxaIndices;
		copyTree.config=config;
		copyTree.stochasticallyMapped=false;
		copyTree.asrDone=false;
		copyTree.asrSeqDone=false;
		copyTree.sortingDone=false;
		treeCopy(this.root, copyTree.root);  
		return copyTree;			
	}

	private void treeCopy(TreeWithLocationsNode from, TreeWithLocationsNode to) {
		for (TreeWithLocationsNode child : from.children) {
			TreeWithLocationsNode newChild = null;
			if (config.seqModelType==SeqModelType.NONE)
				newChild = new TreeWithLocationsNode(child.seq, child.getLoc(),child.taxonIndex,child.time, to,false);
			else
				newChild = new TreeWithLocationsNode(child.seq, child.getLoc(),child.taxonIndex,child.time, to,true);
			to.children.add(newChild);			
			treeCopy(child, newChild);
		}		
	}

	@Override
	public void clearInternalNodes() {
		for (TreeWithLocationsNode node : root) {
			node.mutations=null;
			node.migrations=null;
			if (node.children!=null) {
				if (node.children.size()>0) {
					node.setLoc(TreeWithLocations.UNKNOWN_LOCATION);
					node.seq=null;
				}				
			}
		}

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
			Sequence seq = null;
			if (taxon!=null) {
				location = locationMap.get(taxon.toString());				
				if (location==null) 
					location=TreeWithLocations.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;

				if (seqMap!=null) {
					seq = seqMap.get(taxon.toString());
					if (seq==null) 
						seq=new Sequence(seqLength);
					else
						numIdentifiedSeqs+=1;
				}

				taxonIndex = taxaIndices.get(taxon.toString());
				if (taxonIndex==null) 
					taxonIndex = TreeWithLocations.UNKNOWN_LOCATION;

			}			
			else {
				seq = new Sequence(seqLength);
			}
			if (taxonIndex==null) taxonIndex = UNKNOWN_TAXA;
			if (config.seqModelType==SeqModelType.NONE)
				root.children.add(new TreeWithLocationsNode(seq, location,taxonIndex,root.time+inputTree.getLength(node),root,false));
			else
				root.children.add(new TreeWithLocationsNode(seq, location,taxonIndex,root.time+inputTree.getLength(node),root,true));
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
			if (config.seqModelType==SeqModelType.NONE)
				outputSubTree.children.add(new TreeWithLocationsNode(new Sequence(seqLength), TreeWithLocations.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree,false));
			else
				outputSubTree.children.add(new TreeWithLocationsNode(new Sequence(seqLength), TreeWithLocations.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree,true));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(TransitionModel m, TreeWithLocationsNode root, int nNodes) {
		// for test purpose

		// tree generate parameters for test purpose
		final double testBranchLengthMean = 0.5;
		final double testBranchLengthVariance = 1.0;

		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.getLoc(), location, root.time, to_time));
					if (cern.jet.random.Uniform.staticNextDouble()<=p) {
						if (config.seqModelType==SeqModelType.NONE)
							root.children.add(new TreeWithLocationsNode(new Sequence(seqLength), location,TreeWithLocations.UNKNOWN_LOCATION,to_time,root,false));
						else
							root.children.add(new TreeWithLocationsNode(new Sequence(seqLength), location,TreeWithLocations.UNKNOWN_LOCATION,to_time,root,true));
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
	public void setMigrationModel(Object likelihoodModel_) {
		migrationModel = (TransitionModel) likelihoodModel_;
	}

	public String newickProbs() {	 
		return newickProbs(root,migrationModel.rootfreq(root.time).toArray());
	}

	private String newickProbs(TreeWithLocationsNode treePart, double[] rootFreq) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.getParent().time);
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
			if (treePart.getParent()!=null) {
				parentTime=treePart.getParent().time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+=treePart.parseProbs(rootFreq)+":"+branchLength;
		}		
		return returnValue;
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
			if (!Double.isInfinite(alphas[i]) && !Double.isNaN(alphas[i])) {
				sumExp=sumExp+cern.jet.math.Functions.exp.apply(alphas[i]-minWithoutNegInf);
			}
		}
		double returnValue=minWithoutNegInf+cern.jet.math.Functions.log.apply(sumExp);
		if (Double.isNaN(returnValue) ){
			return minNegative;
		}
		return returnValue;
	}

	@Override
	public String newickSM(int maxBranchRetries) {
		sortChildrenByDescendants();
		asr(); // Ancestral state reconstruction
		stochasticMapping(maxBranchRetries);

		return newickSM(root);
	}

	private String newickSM(TreeWithLocationsNode treePart) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.getParent().time);
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+ "[&states="+Integer.toString(treePart.getLoc())+"]:"+treePart.parseMap()+branchLength);
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
			if (treePart.getParent()!=null) {
				parentTime=treePart.getParent().time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+="[&states="+Integer.toString(treePart.getLoc())+"]:"+treePart.parseMap()+branchLength;
		}		
		return returnValue;
	}

	private void stochasticMapping(int maxBranchRetries) {
		// TODO: test
		// TODO: cite		
		//System.err.println(" Q="+migrationModel.parse());
		if (stochasticallyMapped) return;		
		for (TreeWithLocationsNode node : eachPreorder()) { 
			if (node.getParent()==null) continue;
			node.migrations=new ArrayList<Transition>();
			int currentLoc = node.getParent().getLoc();
			double currentTime = node.getParent().time;
			boolean doneWithBranch = false;
			boolean failedMapping = false;
			Transition event = null;
			int repeats = 0;
			do {
				repeats+=1;
				event = migrationModel.nextEvent(currentTime, currentLoc);				
				if (event.time < node.time) {				
					node.migrations.add(event);					
					currentLoc = event.toTrait;
					currentTime = event.time;
				}
				else {
					if (currentLoc!=node.getLoc()) {
						// If there is a mismatch between stochastic mapping and child ASR than we 
						// restart the entire branch
						node.migrations.clear();
						currentLoc = node.getParent().getLoc();
						currentTime = node.getParent().time;
					} else {						
						doneWithBranch = true;								
					}
				}
				if (repeats>maxBranchRetries) {
					// TODO: 
					System.err.println("Failed to stochasticaly map branch locations after "+Integer.toString(maxBranchRetries)+" iterations\nReasons include: 1. Low retry limit. 2. Unconverged matrices 3. Long branch 4. Aggresive rounding of tip time, consider adding jitter to tip times.\n{("+Integer.toString(node.getParent().taxonIndex)+","+Double.toString(node.getParent().time)+","+Double.toString(node.getParent().getLoc())+"),("+Integer.toString(node.taxonIndex)+","+Double.toString(node.time)+","+Double.toString(node.getLoc())+")}");					
					failedMapping=true;
					node.migrations.clear();
				}
			} while (!doneWithBranch && !failedMapping);

			// TODO: organize this
			if (node.migrations!=null) {
				if (node.migrations.size()>0) {
					if (node.migrations.get(node.migrations.size()-1).toTrait!=node.getLoc()) {
						System.err.println("failed branch location stochastic mapping (last migration)!");
						System.exit(-1);						
					}
					if (node.migrations.get(0).toTrait==node.getParent().getLoc()) {
						System.err.println("failed branch location stochastic mapping (first migration)!");
						System.exit(-1);						
					}
				}
			}			
		}			
		stochasticallyMapped=true;
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
			if (min<logProbs[i] &&  !Double.isInfinite(logProbs[i]))
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
		System.err.println("Error sampling random location from log probs!");
		System.exit(-1);
		return TreeWithLocations.ERR_LOCATION;
	}

	private void asr() {
		if (asrDone) return;		
		// Calculate root state
		// TODO: check this
		DoubleMatrix1D rootFreq = migrationModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbsLOC[i] + Math.log(rootFreq.get(i));
		}		
		root.setLoc(normalizeAndGetRandomSampleFromLogProbs(alphas));		
		for (TreeWithLocationsNode node : eachPreorder()) {
			if (node.getParent()==null) continue;
			TreeWithLocationsNode parent = node.getParent();					
			// TODO: check if clause (here for numerics issues)
			DoubleMatrix1D p = migrationModel.probability(parent.getLoc(), parent.time, node.time);
			for (int i=0; i < numLocations; i++) {								
				alphas[i] = cern.jet.math.Functions.log.apply(p.get(i)) + node.logProbsLOC[i];
			}		
			node.setLoc(normalizeAndGetRandomSampleFromLogProbs(alphas));

		}	
		asrDone=true;
	}	

	@Override
	public String newickASR() {
		// TODO: Check this
		sortChildrenByDescendants();
		asr();
		String returnValue = newickStates(root);
		return returnValue;

	}

	private String newickStates(TreeWithLocationsNode treePart) {
		sortChildrenByDescendants();
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.getParent().time);			
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+ "[&states="+Integer.toString(treePart.getLoc())+"]:"+branchLength);
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
			if (treePart.getParent()!=null) {
				parentTime=treePart.getParent().time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+="[&states="+Integer.toString(treePart.getLoc())+"]:"+branchLength;
		}		
		return returnValue;
	}

	@Override
	public String smTransitions() {


		String returnValue = "{";
		String[][] transitionTimes = new String[numLocations][numLocations];		

		int postOrderIndex=0;
		for (TreeWithLocationsNode node : root) {
			postOrderIndex++;
			if (node==root) continue;
			if (node.migrations!=null) {
				int fromLocation = node.getParent().getLoc();
				for (Transition transition : node.migrations) {
					if (transitionTimes[fromLocation][transition.toTrait]==null) {
						transitionTimes[fromLocation][transition.toTrait]=new String();
					}
					assert(fromLocation!=transition.toTrait);
					if (config.smMigrationNodeNumTipAndSequenceData) { 
						String nodeString = "";
						// is tip
						nodeString+="{"+node.isTip()+",";
						// post order node number
						nodeString+=Integer.toString(postOrderIndex)+",";
						// migraration time
						nodeString+=String.format("%.3f",transition.time)+"},";
						//
						transitionTimes[fromLocation][transition.toTrait]+=nodeString;
					}
					else  {
						transitionTimes[fromLocation][transition.toTrait]+=String.format("%.3f,",transition.time);
					}
					fromLocation=transition.toTrait;
				}				
			}
		}

		for (int i=0;i<numLocations;i++) {
			returnValue+="{";
			for (int j=0;j<numLocations;j++) {						
				if (transitionTimes[i][j]!=null) {
					if (transitionTimes[i][j].length()>=2) { 
						returnValue+="{"+transitionTimes[i][j].substring(0, transitionTimes[i][j].length()-1)+"}";
					}
					else {
						returnValue+="{}";
					}
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
		sortChildrenByDescendants();
		String returnValue = "{";
		String[] tipDwellings = new String[numLocations];	

		for (int i=0;i<numLocations;i++) {	
			tipDwellings[i]=new String();			
		}

		for (TreeWithLocationsNode node : root) {			
			if (node==root) continue;
			if (!node.isTip()) continue;
			int fromLocation = node.getLoc();
			double fromTime = tracebackWithNoChanges(node);
			tipDwellings[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
		}

		for (int i=0;i<numLocations;i++) {
			if (tipDwellings[i].length()>=2) {
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

	private double tracebackWithNoChanges(TreeWithLocationsNode node) {
		// TODO: test this
		if (node.getParent()==null) {
			return node.time;
		}
		if (node.migrations!=null) {
			if (node.migrations.size()>0) {
				return (node.migrations.get(node.migrations.size()-1).time);
			}
			else {
				return (tracebackWithNoChanges(node.getParent()));
			}
		}
		else {
			return (tracebackWithNoChanges(node.getParent()));
		}
	}

	@Override
	public String smLineages() {
		// TODO: test this
		sortChildrenByDescendants();
		String returnValue = "{";
		String[] lineages = new String[numLocations];	

		for (int i=0;i<numLocations;i++) {	
			lineages[i]=new String();			
		}

		for (TreeWithLocationsNode node : root) {
			// Only record lineages from node.parent.time to node.time  
			if (node==root) continue;
			int fromLocation = node.getLoc();
			double fromTime = node.getParent().time;
			if (node.migrations==null) {
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}
			else if (node.migrations.size()==0) {
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}
			else {	
				for (Transition transition : node.migrations) { 														
					lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,transition.time);
					fromLocation=transition.toTrait;
					fromTime=transition.time;
				}
				lineages[fromLocation]+=String.format("{%.3f,%.3f},",fromTime,node.time);
			}			
		}

		for (int i=0;i<numLocations;i++) {	
			if (lineages[i]!=null) {
				if (lineages[i].length()>=2) { 
					returnValue+="{"+lineages[i].substring(0, lineages[i].length()-1)+"}";
				}
				else {
					returnValue+="{}";
				}

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
	public String smTrunkStats(double presentDayTipInterval, double timeToDesignateTrunk) {
		// TODO: test this
		sortChildrenByDescendants();
		List<TreeWithLocationsNode> trunkNodes = markTrunk(presentDayTipInterval, timeToDesignateTrunk);

		String returnValue = "{";


		for (TreeWithLocationsNode node : trunkNodes) {
			if (node.getParent()!=null) {
				if (node.migrations!=null) {
					if (node.migrations.size()==0) {
						if (returnValue.charAt(returnValue.length()-1)=='}') returnValue+=",";
						returnValue+=("{"+node.getLoc()+","+node.getParent().time+","+node.time+"}");
					}
					else {
						double fromTime = node.getParent().time;
						int fromLoc = node.getLoc();
						for (Transition transition : node.migrations) {
							if (returnValue.charAt(returnValue.length()-1)=='}') returnValue+=",";
							returnValue+=("{"+fromLoc+","+fromTime+","+transition.time+"}"); // TODO:
							fromTime = transition.time;
							fromLoc = transition.toTrait;
						}
					}
				}
				else {
					if (returnValue.charAt(returnValue.length()-1)=='}') returnValue+=",";
					returnValue+=("{"+node.getLoc()+","+node.getParent().time+","+node.time+"}");
				}
			}
		}		
		returnValue+="}";
		return returnValue;

	}

	private List<TreeWithLocationsNode> markTrunk(double presentDayTipInterval, double timeToDesignateTrunk) {
		double maxTime = root.time;
		List<TreeWithLocationsNode> returnValue = new ArrayList<TreeWithLocationsNode>();
		for (TreeWithLocationsNode node : root) {
			if (node.time>maxTime) {
				maxTime = node.time;
			}
		}

		List<TreeWithLocationsNode> presentDayTips = new ArrayList<TreeWithLocationsNode>();
		for (TreeWithLocationsNode node : root) {
			if ((maxTime-node.time)<presentDayTipInterval) {
				presentDayTips.add(node);
			}
		}

		for (TreeWithLocationsNode node : presentDayTips) {
			returnValue.addAll(markTrunkAnc(node, maxTime, timeToDesignateTrunk));
		}
		return returnValue;

	}

	private List<TreeWithLocationsNode> markTrunkAnc(TreeWithLocationsNode node, double maxTime, double timeToDesignateTrunk) {
		List<TreeWithLocationsNode> returnValue = new ArrayList<TreeWithLocationsNode>();
		if ((maxTime-node.time)>timeToDesignateTrunk) {
			node.setTrunk();
			returnValue.add(node);
		}
		if (node.getParent()!=null && !node.getParent().isTrunk()) {
			returnValue.addAll(markTrunkAnc(node.getParent(), maxTime, timeToDesignateTrunk));
		}
		return returnValue;
	}

	@Override
	public String smDescendants() {
		// TODO: this
		sortChildrenByDescendants();
		String returnValue = "{";
		String[][] decendants = new String[numLocations][numLocations];		

		for (TreeWithLocationsNode node : root) {
			if (node==root) continue;
			if (node.migrations!=null) {
				int fromLocation = node.getParent().getLoc();
				for (Transition transition : node.migrations) {
					if (decendants[fromLocation][transition.toTrait]==null) {
						decendants[fromLocation][transition.toTrait]=new String();
					}
					if (fromLocation!=transition.toTrait) 
						decendants[fromLocation][transition.toTrait]+=String.format("{%.3f,%d,%.3f,%d,%.3f},",transition.time,getTotalNumTips(node,transition),getTotalBranchLength(node,transition),getLocalNumTips(node,transition),getLocalBranchLength(node,transition));
					fromLocation=transition.toTrait;
				}				
			}
		}

		for (int i=0;i<numLocations;i++) {
			returnValue+="{";
			for (int j=0;j<numLocations;j++) {						
				if (decendants[i][j]!=null) {
					if (decendants[i][j].length()>=2) { 
						returnValue+="{"+decendants[i][j].substring(0,decendants[i][j].length()-1)+"}";
					}
					else {
						returnValue+="{}";
					}
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


	private double getLocalBranchLength(TreeWithLocationsNode node,
			Transition transition) {
		// TODO Auto-generated method stub
		// TODO test this		
		// List with at least one element...
		int index = node.migrations.indexOf(transition);
		if (index!=(node.migrations.size()-1)) {
			return (node.migrations.get(index+1).time - transition.time);			
		}
		else {
			return (node.time-transition.time) + getLocalBranchLength(node);
		}
	}

	private double getLocalBranchLength(TreeWithLocationsNode node) {
		// TODO test this		
		double returnValue = 0;
		for (TreeWithLocationsNode child : node.children) {		
			if (child.migrations!=null) {
				if (child.migrations.size()>0) {
					returnValue+=child.migrations.get(0).time - node.time;
				}
				else 
					returnValue+=(child.time - node.time + getLocalBranchLength(child));
			}		
		}
		return returnValue;
	}

	private int getLocalNumTips(TreeWithLocationsNode node,
			Transition transition) {
		// TODO test this
		int index = node.migrations.indexOf(transition);
		if (index!=(node.migrations.size()-1)) {
			return 0;			
		}
		else {			
			return (node.isTip() ? 1 : 0)+getNumLocalTips(node);
		}
	}

	private int getNumLocalTips(TreeWithLocationsNode node) {
		// TODO test this 		
		int returnValue = 0;
		for (TreeWithLocationsNode child : node.children) {		
			if (child.migrations!=null) {
				if (child.migrations.size()==0) {
					if (child.isTip()) {
						returnValue+=1;
					}
					else {
						returnValue+=getNumLocalTips(child);
					}
				}
			}		
		}
		return returnValue;
	}

	private double getTotalBranchLength(TreeWithLocationsNode node,
			Transition transition) { 
		// TODO test this		
		return (node.time-transition.time) + getTotalBranchLength(node);		
	}

	private double getTotalBranchLength(TreeWithLocationsNode node) {
		// TODO test this
		if (node.isTip()) { 
			return 0;
		}
		else {
			int returnValue=0;
			for (TreeWithLocationsNode child : node.children) {
				returnValue+=(child.time - node.time) + getTotalBranchLength(child);
			}
			return returnValue;
		}
	}

	private int getTotalNumTips(TreeWithLocationsNode node,	Transition transition) {	
		// TODO test this
		return getTotalNumTips(node);
	}

	private int getTotalNumTips(TreeWithLocationsNode node) {
		// TODO test this
		if (node.isTip()) 
			return 1;
		else {
			int returnValue=0;
			for (TreeWithLocationsNode child : node.children) {
				returnValue+=getTotalNumTips(child);
			}
			return returnValue;
		}
	}

	@Override
	public void setCodonModel(Object codonModel) {
		codonLikelihoodModel = (TransitionModel[]) codonModel;
	}

	public double getNumIdentifiedSeqs() {
		return numIdentifiedSeqs;
	}

	@Override
	public String seqMutationStats(int maxBranchRetries) throws Exception {

		// TODO: check SM time at change in node...
		System.out.print("Sorting children by the number of descendants for consistant node ordering...");	
		sortChildrenByDescendants();
		System.out.print("done! ");
		System.out.println("Stochastically mapping mutations:");
		System.out.print("Sequence ASR....");
		asrSeq();
		System.out.print("done! ");
		System.out.print("Sequence SM....");
		stochsticMappingSeq(maxBranchRetries);
		System.out.print("done! ");

		System.out.print("Generating output");
		String returnValue = "{";	
		int i=0;		
		for (TreeWithLocationsNode node : root) {
			i++;
			if (i%10==0) System.out.print(".");
			if (i%1000==0) System.out.println();
			for (int codonPosition = 0; codonPosition<3;codonPosition++) {				
				for (int loc=0; loc<((seqLength-1)/3+1);loc++) {
					if ((loc*3+codonPosition) >= seqLength) continue;
					if (node==root) continue;
					if (node.mutations.get(codonPosition).get(loc)!=null) {						
						int fromNuc = 0;
						fromNuc = node.getParent().seq.getNuc(loc*3+codonPosition);						
						for (Transition transition : node.mutations.get(codonPosition).get(loc)) {
							if (transition.time>config.seqStochasticMappingStartTime) {
								if (returnValue.charAt(returnValue.length()-1)=='}') returnValue+=",";
								returnValue+=String.format("{%.3f,",transition.time);
								returnValue+=Integer.toString(loc*3+codonPosition)+",";
								returnValue+=Sequence.toChar(fromNuc)+",";								
								returnValue+=Sequence.toChar(transition.toTrait)+",";
								/// codon
								if (config.seqMutationsStatsCodonOutput) {
									String fromCodon=mapMutationsToCodon(node.getParent().seq, loc, node.mutations, transition.time);								
									returnValue+=fromCodon+",";
									returnValue+=(codonPosition==0 ? Sequence.toChar(transition.toTrait) : fromCodon.charAt(0));
									returnValue+=(codonPosition==1 ? Sequence.toChar(transition.toTrait) : fromCodon.charAt(1));
									returnValue+=(codonPosition==2 ? Sequence.toChar(transition.toTrait) : fromCodon.charAt(2));
									returnValue+=",";
								}
								// source seq						
								if (config.seqMutationsStatsSeqOutput) {
									String fromSeq=mapMutationsToSequence(node.getParent().seq, node.mutations, transition.time);
									returnValue+=fromSeq+",";
								}							
								// location
								returnValue+=getSequenceTransitionLocation(node,transition.time) + ",";
								// is tip
								returnValue+=node.isTip()+",";
								// post order node number
								returnValue+=Integer.toString(i);
								//
								returnValue+="}";								
							}
							fromNuc=transition.toTrait;
						}									
					}
				}
			}
		}
		returnValue+="}";
		System.out.println("done!");
		return returnValue;
	}

	private String mapMutationsToSequence(Sequence seq,
			ArrayList<ArrayList<ArrayList<Transition>>> mutations, double upToTime) {
		if (mutations==null) 
			return seq.toString();
		else {
			Sequence seqCopy = seq.copy();		
			for (int loc=0; loc<((seqLength-1)/3+1);loc++) {
				for (int cp = 0; cp<3; cp++) {
					if ((loc*3+cp) >= seqLength) continue;	

					for (Transition mutation : mutations.get(cp).get(loc)) {
						if (mutation.time<upToTime) {
							seqCopy.set(loc*3+cp, mutation.toTrait);
						}
					}
				}
			}		
			return seqCopy.toString();
		}
	}

	private String mapMutationsToCodon(Sequence seq, int loc,
			ArrayList<ArrayList<ArrayList<Transition>>> mutations, double upToTime) {
		Sequence seqCopy = seq.copy();
		for (int cp = 0; cp<3; cp++) {
			for (Transition mutation : mutations.get(cp).get(loc)) {
				if (mutation.time<upToTime) {
					seqCopy.set(loc*3+cp, mutation.toTrait);
				}
			}
		}
		String returnValue=new String();
		for (int cp = 0; cp<3; cp++) {
			try {
				returnValue+=Sequence.toChar(seqCopy.getNuc(loc*3+cp));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return returnValue;
	}

	private int getSequenceTransitionLocation(TreeWithLocationsNode node,
			double time) {
		if (node.migrations==null) {
			return node.getLoc();
		}
		else if (node.migrations.size()==0) {
			return node.getLoc();
		}
		else if (node.migrations.get(0).time>time) 	
			return node.getParent().getLoc();
		else
			return binarySearchMigrations(node.migrations,time).toTrait;			
	}



	private Transition binarySearchMigrations(List<Transition> migrations, double time) {
		int start=0;
		int end=migrations.size()-1;
		int current=(start+end)/2;
		while ((end-start)>1) {
			if (migrations.get(current).time<time) {
				start=current;
			}
			else {
				end=current;
			}
			current = (start+end)/2;			
		}
		if (migrations.get(end).time<time) {
			return migrations.get(end);
		}
		else {
			return migrations.get(start);
		}
	}

	@Override
	public String pies() {

		String returnValue = "{";	
		for (int i=0; i<3; i++) {
			returnValue+=Integer.toString(i)+": {";			
			for (int j=0; j<3; j++) {
				returnValue+=Integer.toString(j)+": "+codonLikelihoodModel[i].rootfreq(0).get(j);
				if (j!=2) returnValue+=",";
				returnValue+=" ";
			}			
			returnValue+="}";
			if (i!=2) returnValue+=",";
		}
		returnValue+="}";			
		return returnValue;
	}




	private void stochsticMappingSeq(int maxBranchRetries) throws Exception {
		// TODO: test
		// TODO: cite	
		
		sortChildrenByDescendants();
		for (TreeWithLocationsNode node : eachPreorder()) {
			if (node.getParent()==null) continue;
			if (node.mutations==null) {
				node.mutations = new ArrayList<ArrayList<ArrayList<TransitionModel.Transition>>>();
				for (int i=0;i<3;i++) {
					node.mutations.add(new ArrayList<ArrayList<TransitionModel.Transition>>());
					for (int j=0;j<((seqLength-1)/3+1);j++) {
						node.mutations.get(i).add(new ArrayList<TransitionModel.Transition>());
					}
				}
			}

			for (int codonPosition = 0; codonPosition<3;codonPosition++) {
				for (int loc=0; loc<((seqLength-1)/3+1);loc++) { 
					if ((loc*3+codonPosition) >= seqLength) continue;
					node.mutations.get(codonPosition).get(loc).clear();
					int currentNuc =  node.getParent().seq.getNuc(loc*3+codonPosition);
					double currentTime = node.getParent().time;
					boolean doneWithBranch = false;
					boolean failedMapping = false;
					Transition event = null;
					int repeats = 0;
					do {
						repeats+=1;
						event = codonLikelihoodModel[codonPosition].nextEvent(currentTime, currentNuc);						
						if (event.time < node.time) {					
							node.mutations.get(codonPosition).get(loc).add(event);
							currentNuc = event.toTrait;
							currentTime = event.time;
						} 
						else if (currentNuc!=node.seq.getNuc(loc*3+codonPosition)) {
							// If there is a mismatch between stochastic mapping and child ASR than we 
							// restart the entire branch
							node.mutations.get(codonPosition).get(loc).clear();							
							currentNuc = node.getParent().seq.getNuc(loc*3+codonPosition);
							currentTime = node.getParent().time;
						} 
						else {							
							doneWithBranch = true;								
						}
						if (repeats>maxBranchRetries) {
							System.err.println("Failed to stochasticaly map branch sequence after "+Integer.toString(maxBranchRetries)+" iterations at loci: "+(loc*3+codonPosition)+"\nusing last attempt for output. sing last attempt for output.\nReasons include: 1. Low retry limit. 2. Unconverged matrices 3. Long branch 4. Aggresive rounding of tip time, consider adding jitter to tip times\n{("+Integer.toString(node.getParent().taxonIndex)+","+Double.toString(node.getParent().time)+","+Double.toString(node.getParent().seq.getNuc(loc*3+codonPosition))+"),("+Integer.toString(node.taxonIndex)+","+Double.toString(node.time)+","+Double.toString(node.seq.getNuc(loc*3+codonPosition))+")}");											
							failedMapping=true;
						}
					} while (!doneWithBranch && !failedMapping);
				}				
			}			
		}		
	}

	private void asrSeq() throws Exception {

		if (asrSeqDone) return;	

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
				root.seq.set(loc*3+codonPos,normalizeAndGetRandomSampleFromLogProbs(alphas));
			}			
		}			

		for (TreeWithLocationsNode node : eachPreorder()) {
			if (node.getParent()==null) continue;
			for (int codonPosition = 0; codonPosition<3;codonPosition++) {
				for (int loc=0; loc<((seqLength-1)/3+1);loc++) { 
					if ((loc*3+codonPosition) >= seqLength) continue;

					// TODO: check this
					TreeWithLocationsNode parent = node.getParent();

					double[] alphas = new double[4];	
					// TODO: check need for Util.minValue
					DoubleMatrix1D p = codonLikelihoodModel[codonPosition].probability(parent.seq.getNuc(loc*3+codonPosition), parent.time, node.time+Util.minValue);

					for (int i=0; i < 4; i++) {								
						alphas[i] = cern.jet.math.Functions.log.apply(p.get(i)) + node.logProbsCP[codonPosition][loc][i];
					}		
					node.seq.set(loc*3+codonPosition,normalizeAndGetRandomSampleFromLogProbs(alphas));
				}
			}
		}		

		asrSeqDone=true;
	}

	@Override
	public Double seqLikelihood() {		
		return this.seqLogLike;
	}

	@Override
	public Double locLikelihood() {
		return this.locationLogLike;
	}

	@Override
	public String seqMigrationsSeqOutput() throws Exception {
		
		sortChildrenByDescendants();
		String returnValue = "{";
		String[][] transitionTimes = new String[numLocations][numLocations];		

		int postOrderIndex=0;
		for (TreeWithLocationsNode node : root) {
			postOrderIndex++;
			if (node==root) continue;
			if (node.migrations!=null) {
				int fromLocation = node.getParent().getLoc();
				for (Transition transition : node.migrations) {
					if (transitionTimes[fromLocation][transition.toTrait]==null) {
						transitionTimes[fromLocation][transition.toTrait]=new String();
					}
					assert(fromLocation!=transition.toTrait);
					if (config.smMigrationNodeNumTipAndSequenceData) { 
						String nodeString = "";
						// is tip
						nodeString+="{"+node.isTip()+",";
						// post order node number
						nodeString+=Integer.toString(postOrderIndex)+",";
						// migraration time
						nodeString+=String.format("%.3f",transition.time)+",";
						// sequence at migration
						nodeString+=mapMutationsToSequence(node.getParent().seq, node.mutations, transition.time)+"},";
						// add to list
						transitionTimes[fromLocation][transition.toTrait]+=nodeString;
					}
					else  {
						transitionTimes[fromLocation][transition.toTrait]+=String.format("%.3f,",transition.time);
					}
					fromLocation=transition.toTrait;
				}				
			}
		}

		for (int i=0;i<numLocations;i++) {
			returnValue+="{";
			for (int j=0;j<numLocations;j++) {						
				if (transitionTimes[i][j]!=null) {
					if (transitionTimes[i][j].length()>=2) { 
						returnValue+="{"+transitionTimes[i][j].substring(0, transitionTimes[i][j].length()-1)+"}";
					}
					else {
						returnValue+="{}";
					}
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
	public AltTreeOutput smAlternativeTreeOutput() {

		setLayoutByDescendants();

		AltTreeOutput altTreeOutput = new AltTreeOutput();

		int postOrderIndex=0;
		for (TreeWithLocationsNode node : root) {
			postOrderIndex++;	
			altTreeOutput.addNode(postOrderIndex, node.time, node.getLoc(), node.getLayout(), mapMutationsToSequence(node.seq, node.mutations, node.time),node.isTip());
			node.postOrderIndex=postOrderIndex;
		}

		for (TreeWithLocationsNode node : root) {
			for (TreeWithLocationsNode child : node.children) {
				altTreeOutput.addBranch(node.postOrderIndex, child.postOrderIndex);
			}
		}


		return altTreeOutput;

	}

	// sets node layout based on a postorder traversal & number of descendants 
	public void setLayoutByDescendants() {

		// sort children so that the child with the more descendants is left most
		sortChildrenByDescendants();

		// set layout of tips based on post order traversal order of tips 
		float y = 0;
		for (TreeWithLocationsNode s : root) {
			if (s.isTip()) {
				s.setLayout(y);
				y++;
			}
		}

		// update layout of internal nodes based on average of children layout
		for (TreeWithLocationsNode s : root) {
			if (s.children.size() > 0) {
				float mean = 0;
				for (TreeWithLocationsNode child : s.children) {
					mean += child.getLayout();
				}
				mean /= s.children.size();
				s.setLayout(mean);
			}
		}

	}	


	// sorts children lists so that first member is child with more descendents than second member
	public static void sortChildrenByDescendants(TreeWithLocationsNode r) {

		Collections.sort(r.children, descendantOrder);

		Stack<TreeWithLocationsNode> S = new Stack<TreeWithLocationsNode>();
		TreeWithLocationsNode u;

		S.push(r);
		while (!S.isEmpty()) {
			u = S.pop();
			for (TreeWithLocationsNode s : u.children) {
				Collections.sort(s.children, descendantOrder);
				S.push(s);                
			}
		}        
	}	

	public void sortChildrenByDescendants() {		
		if (!sortingDone) 
			sortChildrenByDescendants(root);
		sortingDone=true;
	}


	public static int getNumberOfDescendants(TreeWithLocationsNode r) {

		int numberOfDescendants = 0;

		Stack<TreeWithLocationsNode> S = new Stack<TreeWithLocationsNode>();
		TreeWithLocationsNode u;

		S.push(r);
		while (!S.isEmpty()) {
			u = S.pop();
			for (TreeWithLocationsNode s : u.children) {
				numberOfDescendants+=1;
				S.push(s);                
			}
		}        

		return numberOfDescendants;
	}


}



