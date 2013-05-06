package seasmig.treelikelihood;

import java.util.HashMap;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.util.Util;


@SuppressWarnings("serial")
public class TreeWithLocations implements LikelihoodTree {

	public static final int UNKNOWN_TAXA = -1;
	public static final int UNKNOWN_LOCATION = -1;

	// Tree & Model
	TreeWithLocationsNode root = null;		
	private MigrationBaseModel likelihoodModel = null;

	int numLocations = 0;
	private int numIdentifiedLocations;
	
	// Taxa
	HashMap<String, Integer> taxaIndices = new HashMap<String,Integer>();
	
	// TODO: this..
	double[] ZERO_LOG_PROBS;
	private double logLike = 0;
	
	// Tree generate parameters for test purpose
	static final private double testBranchLengthMean = 0.1;
	static final private double testBranchLengthVariance = 3.0;

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocations(MigrationBaseModel createTreeModel, int numNodes) {		
		numLocations=createTreeModel.getNumLocations();
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		double p=createTreeModel.rootfreq(0)[0];
		int rootLocation =0;
		for (int i=0;i<numLocations;i++) {
			if (cern.jet.random.Uniform.staticNextDouble()<=p) {
				rootLocation=i;
				break;
			}
			else {
				p=p+createTreeModel.rootfreq(0)[i];
			}
		}
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
		double p=likelihoodModel.rootfreq(0)[0];
		int rootLocation =0;
		for (int i=0;i<numLocations;i++) {
			if (cern.jet.random.Uniform.staticNextDouble()<=p) {
				rootLocation=i;
				break;
			}
			else {
				p=p+likelihoodModel.rootfreq(0)[i];
			}
		}
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
				double d = cern.jet.random.Uniform.staticNextDouble();
				double p=0;
				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(likelihoodModel.logprobability(root.loc, location, root.time, child.time));
					if (d<=p) {
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
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String,Integer> taxaIndices_, String locationAttributeName, int num_locations_) {
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
		root = new TreeWithLocationsNode(Integer.parseInt((String)tree.getRootNode().getAttribute(locationAttributeName))-1,rootTaxonIndex,0,null);
		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree,HashMap<String,Integer> taxaIndices_, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
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
		logLike  = 0;
		for (TreeWithLocationsNode node : root) {
			if (node.loc==TreeWithLocations.UNKNOWN_LOCATION) {
				node.logProbs =new double[numLocations];
			}
			else {
				node.logProbs= ZERO_LOG_PROBS.clone();
				node.logProbs[node.loc]=0;
			}

			if (node.children.size()!=0) { // this is an internal node			
				for (int from = 0; from < numLocations; from++) {
					for (TreeWithLocationsNode child : node.children ) {
						// for now caching is done inside likelihood model...
						double[][] p = likelihoodModel.transitionMatrix(node.time, child.time); 
						double[] alphas = new double[numLocations];						
						for (int to = 0; to < numLocations; to++) {
							alphas[to]=(Math.log(p[from][to]) + child.logProbs[to]);
						}						
						node.logProbs[from] += Util.logSumExp(alphas);
					}								
				}
			}		
		
		}

		// Calculate root base frequency contribution... 
		double[] rootFreq = likelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logProbs[i] + Math.log(rootFreq[i]);
		}			
		logLike += Util.logSumExp(alphas);

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
			if (locationAttributeName!=null)
				outputSubTree.children.add(new TreeWithLocationsNode(Integer.parseInt((String)node.getAttribute(locationAttributeName))-1,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			else 
				outputSubTree.children.add(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,taxonIndex,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, TreeWithLocationsNode root, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.loc, location, root.time, to_time));
					if (d<=p) {
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
	
	public String newick() {	 
		return newick(root,likelihoodModel.rootfreq(root.time)) + "\n";
	}

	// TODO: (348[&antigenic={-6.00510611736,5.84199000915},rate=1.1478703001047978,states="japan_korea"]:2.44, ....
	private String newick(TreeWithLocationsNode treePart, double[] rootFreq) {
		String returnValue = new String();

		if (treePart.isTip()) {
			String branchLength = String.format("%.3f", treePart.time-treePart.parent.time);
			returnValue+=(Integer.toString(treePart.getTaxonIndex())+treePart.parseAncestralStates(rootFreq)+":"+branchLength);
		}
		else if (treePart.children.size()>0) {
			returnValue+="(";
			returnValue+=newick(treePart.children.get(0),rootFreq);
			for (int i = 1; i < treePart.children.size(); i++){
				returnValue+=",";
				returnValue+=newick(treePart.children.get(i),rootFreq);	
			}
			returnValue+=")";
			double parentTime=0;
			if (treePart.parent!=null) {
				parentTime=treePart.parent.time;
			}
			String branchLength = String.format("%.3f", treePart.time-parentTime);
			returnValue+=treePart.parseAncestralStates(rootFreq)+":"+branchLength;
		}		
		return returnValue;
	}

	@Override
	public double cachedLogLikelihood() {		
		return logLike;
	}

}



