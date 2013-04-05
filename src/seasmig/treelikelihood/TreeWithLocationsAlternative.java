package seasmig.treelikelihood;

import java.util.HashMap;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;


@SuppressWarnings("serial")
public class TreeWithLocationsAlternative implements LikelihoodTree {

	// TODO: this..
	double[] UNKNOWN_LOCATION_PROBS;

	// Tree generate parameters for test purpose
	static final private double testBranchLengthMean = 0.1;
	static final private double testBranchLengthVariance = 3.0;

	// Tree & Model
	LocationTreeNode root = null;		
	private MigrationBaseModel likelihoodModel = null;

	int numLocations = 0;
	private int numIdentifiedLocations;

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocationsAlternative(MigrationBaseModel createTreeModel, int numNodes) {		
		numLocations=createTreeModel.getNumLocations();
		UNKNOWN_LOCATION_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			UNKNOWN_LOCATION_PROBS[i]=1.0;
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
		root = new LocationTreeNode(rootLocation,0,null);
		makeRandomTree(createTreeModel, root, numNodes);		
	}

	// Generate random tree states based on input tree topology and model .... 
	public TreeWithLocationsAlternative(MigrationBaseModel createTreeModel, jebl.evolution.trees.SimpleRootedTree tree) {
		numLocations=createTreeModel.getNumLocations();
		UNKNOWN_LOCATION_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			UNKNOWN_LOCATION_PROBS[i]=1.0;
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
		root = new LocationTreeNode(rootLocation,0,null);
		makeSubTree(tree,(String)null, root,tree.getRootNode());
	}

	private void fillRandomTraits(LocationTreeNode root) {
		if (root.children!=null) {
			for (LocationTreeNode child : root.children) {	
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
			for (LocationTreeNode child : root.children) {
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
		UNKNOWN_LOCATION_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			UNKNOWN_LOCATION_PROBS[i]=1.0;
		}
		root = new LocationTreeNode(Integer.parseInt((String)tree.getRootNode().getAttribute(locationAttributeName))-1,0,null);
		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocationsAlternative(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
		numLocations=num_locations_;
		UNKNOWN_LOCATION_PROBS = new double[numLocations];
		for (int i=0;i<numLocations;i++){
			UNKNOWN_LOCATION_PROBS[i]=1.0;
		}
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=LocationTreeNode.UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;
		root = new LocationTreeNode(location,0,null);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
	}

	public TreeWithLocationsAlternative(LocationTreeNode root_, MigrationBaseModel likelihoodModel_) {
		root = root_;
		likelihoodModel=likelihoodModel_;
	}

	@Override
	public double logLikelihood() {
		// TODO: check from/to nomenclature...

		double logLike = 0;

		// Calculate tree likelihood for site
		double tempRetLike = 0;
		for (LocationTreeNode node : root) {
			if (node.children.size()!=0) { // this is an internal node\
				// this doesn't assume no information in internal nodes...
				if (node.loc==LocationTreeNode.UNKNOWN_LOCATION) {
					node.probs = UNKNOWN_LOCATION_PROBS.clone();
				}
				else {
					node.probs=new double[numLocations];
					node.probs[node.loc]=1.0;
				}
				
				for (int from = 0; from < numLocations; from++) {
					for (LocationTreeNode child : node.children ) {
						// for now caching is done inside likelihood model...
						double[][] p = likelihoodModel.transitionMatrix(node.time, child.time);
						double tempLike = 0; 
						for (int to = 0; to < numLocations; to++) {
							tempLike += (p[from][to] * child.probs[to]);
						}
						node.probs[from] *= tempLike;
					}
				}
			}
			else { // this is a tip
				node.probs = new double[numLocations];
				node.probs[node.loc]=1.0;
			}
		}

		// Calculate root base frequency contribution... 
		tempRetLike = 0;
		double[] rootFreq = likelihoodModel.rootfreq(root.time).clone();
		for(int i = 0; i < numLocations; i++) {
			tempRetLike += root.probs[i]*rootFreq[i];
		}
		logLike += Math.log(tempRetLike);


		return logLike;		
	}

	private void removeInternalLocations(LocationTreeNode node) {
		if (node.children.size()!=0) {
			node.loc=LocationTreeNode.UNKNOWN_LOCATION;
			for (LocationTreeNode child : node.children) {
				removeInternalLocations(child);				
			}
		}				
	}

	@Override
	public LikelihoodTree copy() {
		// TODO: check this...
		return this;
	}

	@Override
	public String print() {
		return print(root);
	}

	public String print(LocationTreeNode treePart) {
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
			HashMap<String, Integer> locationMap, LocationTreeNode root,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = LocationTreeNode.UNKNOWN_LOCATION;
			if (taxon!=null) {
				location = locationMap.get(inputTree.getTaxon(node).toString());
				if (location==null) 
					location=LocationTreeNode.UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;
			}			
			root.children.add(new LocationTreeNode(location,root.time+inputTree.getLength(node)/*+jitter*(cern.jet.random.Uniform.staticNextDouble()-0.5)*/,root));
			makeSubTree(inputTree,locationMap, root.children.get(root.children.size()-1), node);			
		}
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, String locationAttributeName, LocationTreeNode outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		// Null attribute means 0 location tree
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {
			if (locationAttributeName!=null)
				outputSubTree.children.add(new LocationTreeNode(Integer.parseInt((String)node.getAttribute(locationAttributeName))-1,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			else 
				outputSubTree.children.add(new LocationTreeNode(0,outputSubTree.time+inputTree.getLength(node),outputSubTree));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, LocationTreeNode root, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				// Decide on branch length
				double to_time = root.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<numLocations;location++) {
					p=p+Math.exp(m.logprobability(root.loc, location, root.time, to_time));
					if (d<=p) {
						root.children.add(new LocationTreeNode(location,to_time,root));
						break;
					}
				}			
			}

			for (LocationTreeNode child : root.children) {
				makeRandomTree(m,child,(int) Math.round(nNodes/2.0));
			}
		}

	}

//	private void removeInternalLocations(Node node) {
//		if (node.children.size()!=0) {
//			node.location=UNKNOWN_LOCATION;
//			for (Node child : node.children) {
//				removeInternalLocations(child);				
//			}
//		}				
//	}

//	private void clearCachedLikelihood(Node node) {	
//		node.cachedConditionalLogLikelihood=new double[numLocations];
//		if (node.children.size()!=0) {			
//			for (Node child : node.children) {
//				clearCachedLikelihood(child);				
//			}
//		}				
//	}

	@Override
	public int getNumLocations() {		
		return numLocations;
	}

	public int getNumIdentifiedLocations() {
		// TODO: this for trees with states...
		return numIdentifiedLocations;
	}

	@Override
	public void setLikelihoodModel(Object likelihoodModel) {
		// TODO Auto-generated method stub
		
	}

}



