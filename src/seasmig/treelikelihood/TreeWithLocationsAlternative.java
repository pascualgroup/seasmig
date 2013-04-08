package seasmig.treelikelihood;

import java.util.HashMap;

import seasmig.util.Util;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;


@SuppressWarnings("serial")
public class TreeWithLocationsAlternative implements LikelihoodTree {

	// TODO: this..
	double[] ZERO_LOG_PROBS;

	// Tree generate parameters for test purpose
	static final private double testBranchLengthMean = 0.1;
	static final private double testBranchLengthVariance = 3.0;

	// Tree & Model
	TreeWithLocationsAlternativeNode root = null;		
	private MigrationBaseModel likelihoodModel = null;

	int numLocations = 0;
	private int numIdentifiedLocations;

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocationsAlternative(MigrationBaseModel createTreeModel, int numNodes) {		
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
		root = new TreeWithLocationsAlternativeNode(rootLocation,0,null);
		makeRandomTree(createTreeModel, root, numNodes);		
	}

	// Generate random tree states based on input tree topology and model .... 
	public TreeWithLocationsAlternative(MigrationBaseModel createTreeModel, jebl.evolution.trees.SimpleRootedTree tree) {
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
		root = new TreeWithLocationsAlternativeNode(rootLocation,0,null);
		makeSubTree(tree,(String)null, root,tree.getRootNode());
	}

	private void fillRandomTraits(TreeWithLocationsAlternativeNode root) {
		if (root.children!=null) {
			for (TreeWithLocationsAlternativeNode child : root.children) {	
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
			for (TreeWithLocationsAlternativeNode child : root.children) {
				fillRandomTraits(child);
			}
		}
	}

	public void fillRandomTraits() {
		fillRandomTraits(root);
	}


	// Load a tree from a basic jebl tree
	// locations are loaded from nexsus tree trait location_attribute name
	public TreeWithLocationsAlternative(jebl.evolution.trees.SimpleRootedTree tree, String locationAttributeName, int num_locations_) {
		numLocations=num_locations_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		root = new TreeWithLocationsAlternativeNode(Integer.parseInt((String)tree.getRootNode().getAttribute(locationAttributeName))-1,0,null);
		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocationsAlternative(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
		numLocations=num_locations_;
		ZERO_LOG_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			ZERO_LOG_PROBS[i]=Double.NEGATIVE_INFINITY;
		}
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=TreeWithLocationsAlternativeNode.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;
		root = new TreeWithLocationsAlternativeNode(location,0,null);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
	}

	public TreeWithLocationsAlternative(TreeWithLocationsAlternativeNode root_, MigrationBaseModel likelihoodModel_) {
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
		// TODO: check from/to nomenclature...

		double logLike = 0;
		
		for (TreeWithLocationsAlternativeNode node : root) {
			if (node.loc==TreeWithLocationsAlternativeNode.UNKNOWN_LOCATION) {
				node.logprobs =new double[numLocations];
			}
			else {
				node.logprobs= ZERO_LOG_PROBS.clone();
				node.logprobs[node.loc]=0;
			}

			if (node.children.size()!=0) { // this is an internal node			
				for (int from = 0; from < numLocations; from++) {
					for (TreeWithLocationsAlternativeNode child : node.children ) {
						// for now caching is done inside likelihood model...
						double[][] p = likelihoodModel.transitionMatrix(node.time, child.time); 
						double[] alphas = new double[numLocations];						
						for (int to = 0; to < numLocations; to++) {
							alphas[to]=(Math.log(p[from][to]) + child.logprobs[to]);
						}						
						node.logprobs[from] += Util.logSumExp(alphas);
					}								
				}
			}		

		}

		// Calculate root base frequency contribution... 
		double[] rootFreq = likelihoodModel.rootfreq(root.time);
		double[] alphas = new double[numLocations];	
		for (int i = 0; i < numLocations; i++) {
			alphas[i]=root.logprobs[i] + Math.log(rootFreq[i]);
		}			
		logLike += Util.logSumExp(alphas);

		return logLike;		
	}

	private void removeInternalLocations(TreeWithLocationsAlternativeNode node) {
		if (node.children.size()!=0) {
			node.loc=TreeWithLocationsAlternativeNode.UNKNOWN_LOCATION;
			for (TreeWithLocationsAlternativeNode child : node.children) {
				removeInternalLocations(child);				
			}
		}				
	}

	@Override
	public LikelihoodTree copy() {
		// TOOD: test this...
		TreeWithLocationsAlternative copyTree = new TreeWithLocationsAlternative();
		copyTree.likelihoodModel=this.likelihoodModel;
		copyTree.numIdentifiedLocations=this.numIdentifiedLocations;
		copyTree.numLocations=this.numLocations;		
		copyTree.ZERO_LOG_PROBS=this.ZERO_LOG_PROBS;
		copyTree.root = new TreeWithLocationsAlternativeNode(root.loc,root.time,null);		
		treeCopy(this.root, copyTree.root);  
		return copyTree;
	}

	private void treeCopy(TreeWithLocationsAlternativeNode from, TreeWithLocationsAlternativeNode to) {
		for (TreeWithLocationsAlternativeNode child : from.children) {
			TreeWithLocationsAlternativeNode newChild = new TreeWithLocationsAlternativeNode(child.loc,child.time, to);
			to.children.add(newChild);			
			treeCopy(child, newChild);
		}		
	}


	@Override
	public String print() {
		return print(root);
	}

	public String print(TreeWithLocationsAlternativeNode treePart) {
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
	protected TreeWithLocationsAlternative() {
	}

	private void makeSubTree(SimpleRootedTree inputTree,
			HashMap<String, Integer> locationMap, TreeWithLocationsAlternativeNode root,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = TreeWithLocationsAlternativeNode.UNKNOWN_LOCATION;
			if (taxon!=null) {
				location = locationMap.get(inputTree.getTaxon(node).toString());
				if (location==null) 
					location=TreeWithLocationsAlternativeNode.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;
			}			
			root.children.add(new TreeWithLocationsAlternativeNode(location,root.time+inputTree.getLength(node)/*+jitter*(cern.jet.random.Uniform.staticNextDouble()-0.5)*/,root));
			makeSubTree(inputTree,locationMap, root.children.get(root.children.size()-1), node);			
		}
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, String locationAttributeName, TreeWithLocationsAlternativeNode outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		// Null attribute means 0 location tree
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {
			if (locationAttributeName!=null)
				outputSubTree.children.add(new TreeWithLocationsAlternativeNode(Integer.parseInt((String)node.getAttribute(locationAttributeName))-1,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			else 
				outputSubTree.children.add(new TreeWithLocationsAlternativeNode(0,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, TreeWithLocationsAlternativeNode root, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.loc, location, root.time, to_time));
					if (d<=p) {
						root.children.add(new TreeWithLocationsAlternativeNode(location,to_time,root));
						break;
					}
				}			
			}

			for (TreeWithLocationsAlternativeNode child : root.children) {
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

}



