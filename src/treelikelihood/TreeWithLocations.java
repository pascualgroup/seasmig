package treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;


@SuppressWarnings("serial")
public class TreeWithLocations implements LikelihoodTree {

	// Tree generate parameters for test purpose
	static final private double testBranchLengthMean = 0.1;
	static final private double testBranchLengthVariance = 3.0;

	// Tree & Model
	Node root = null;		
	int num_locations = 0;
	private MigrationBaseModel likelihoodModel = null;
	private int UNKNOWN_LOCATION;
	private int numIdentifiedLocations; 

	// Generate a random tree based on createTreeModel .... 
	public TreeWithLocations(MigrationBaseModel createTreeModel, int numNodes) {		
		num_locations=createTreeModel.getNumLocations();
		UNKNOWN_LOCATION=num_locations;
		root = new Node(UNKNOWN_LOCATION,0,num_locations);
		makeRandomTree(createTreeModel, root, numNodes);		
	}

	// Generate random tree states based on input tree topology and model .... 
	public TreeWithLocations(MigrationBaseModel createTreeModel, jebl.evolution.trees.SimpleRootedTree tree) {
		num_locations=createTreeModel.getNumLocations();
		UNKNOWN_LOCATION=num_locations;
		likelihoodModel=createTreeModel;
		root = new Node(0,0,num_locations);
		makeSubTree(tree,(String)null, root,tree.getRootNode());
	}

	private void fillRandomTraits(Node parent) {
		if (parent.children!=null) {
			for (Node child : parent.children) {	
				double d = cern.jet.random.Uniform.staticNextDouble();
				double p=0;
				for (int location=0;location<num_locations;location++) {
					p=p+Math.exp(likelihoodModel.logprobability(parent.location, location, parent.time, child.time));
					if (d<=p) {
						child.location=location;
						break;
					}
				}			
			}
			for (Node child : parent.children) {
				fillRandomTraits(child);
			}
		}
	}
	
	public void fillRandomTraits() {
		fillRandomTraits(root);
	}


	// Load a tree from a basic jebl tree
	// locations are loaded from nexsus tree trait location_attribute name
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree, String locationAttributeName, int num_locations_) {
		num_locations=num_locations_;
		UNKNOWN_LOCATION=num_locations;
		root = new Node(Integer.parseInt((String)tree.getRootNode().getAttribute(locationAttributeName))-1,0,num_locations);
		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocations(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
		num_locations=num_locations_;
		UNKNOWN_LOCATION=num_locations;
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=UNKNOWN_LOCATION;
		else
			numIdentifiedLocations+=1;
		root = new Node(location,0,num_locations);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
	}

	@Override 
	public void setLikelihoodModel(Object likelihoodModel_) {
		likelihoodModel = (MigrationBaseModel) likelihoodModel_;
	}

	@Override
	public double logLikelihood() {
		double[] alphas=new double[num_locations];
		double min = Double.POSITIVE_INFINITY;
		if (root.location==UNKNOWN_LOCATION) {
			for (int rootLocation=0;rootLocation<num_locations;rootLocation++) {				
				double alpha=conditionalLogLikelihood(root, rootLocation);
				if (Double.isNaN(alpha))  
					alpha=Util.minNegative;
				alphas[rootLocation]=alpha;	
				if (alpha<min) min=alpha;				
			}
			return Util.logSumExp(alphas,min);
		}
		else {
			return conditionalLogLikelihood(root,root.location);					
		}		
	}

	@Override
	public LikelihoodTree copy() {
		TreeWithLocations newTree = new TreeWithLocations();
		newTree.num_locations=num_locations;
		newTree.UNKNOWN_LOCATION=UNKNOWN_LOCATION;
		newTree.root=new Node(root.location,root.time,num_locations);
		treeCopyNoCache(root,newTree.root);		
		return newTree;
	}

	private double conditionalLogLikelihood(Node node, int nodeLocation) {

		if (node.cachedConditionalLogLikelihood[nodeLocation]!=0) {			
			return node.cachedConditionalLogLikelihood[nodeLocation];
		}
		else {
			double loglikelihood=0;
			for (Node child : node.children) {
				if (child.location!=UNKNOWN_LOCATION) {
					assert child.time>node.time;
					loglikelihood=loglikelihood+conditionalLogLikelihood(child,child.location)+likelihoodModel.logprobability(nodeLocation, child.location, node.time, child.time);
				}
				else {
					double[] alphas=new double[num_locations];
					double min = Double.POSITIVE_INFINITY;
					for (int childLocation=0;childLocation<num_locations;childLocation++) {
						double alpha = likelihoodModel.logprobability(nodeLocation, childLocation, node.time, child.time)+conditionalLogLikelihood(child, childLocation);
						if (Double.isNaN(alpha))  
							alpha=Util.minNegative;
						alphas[childLocation]=alpha;
						if (alpha<min) min=alpha;
					}
					loglikelihood=loglikelihood+Util.logSumExp(alphas,min);					
				}			 
			}
			node.cachedConditionalLogLikelihood[nodeLocation]=loglikelihood;
			return loglikelihood;		
		}
	}


	@Override
	public String print() {
		return print(root);
	}

	public String print(Node treePart) {
		String returnValue = Integer.toString(treePart.location);
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

	public void clearCachedLikelihood() {
		clearCachedLikelihood(root);	
	}


	protected TreeWithLocations() {
	}

	private void treeCopyNoCache(Node from, Node to) {
		for (Node child : from.children) {
			Node newChild = new Node(child.location,child.time,num_locations);
			newChild.parent=from;
			to.children.add(newChild);			
			treeCopyNoCache(child, newChild);
		}		
	}

	private void makeSubTree(SimpleRootedTree inputTree,
			HashMap<String, Integer> locationMap, Node outputSubTree,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer location = UNKNOWN_LOCATION;
			if (taxon!=null) {
				location = locationMap.get(inputTree.getTaxon(node).toString());
				if (location==null) 
					location=UNKNOWN_LOCATION;
				else
					numIdentifiedLocations+=1;
			}			
			outputSubTree.children.add(new Node(location,outputSubTree.time+inputTree.getLength(node)/*+jitter*(cern.jet.random.Uniform.staticNextDouble()-0.5)*/,num_locations));
			makeSubTree(inputTree,locationMap, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}


	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, String locationAttributeName, Node outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		// Null attribute means 0 location tree
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {
			if (locationAttributeName!=null)
				outputSubTree.children.add(new Node(Integer.parseInt((String)node.getAttribute(locationAttributeName))-1,outputSubTree.time+inputTree.getLength(node),num_locations));
			else 
				outputSubTree.children.add(new Node(0,outputSubTree.time+inputTree.getLength(node),num_locations));
			makeSubTree(inputTree, locationAttributeName, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}

	}

	public void makeRandomTree(MigrationBaseModel m, Node from, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				// Decide on branch length
				double to_time = from.time+cern.jet.random.Gamma.staticNextDouble(testBranchLengthMean*testBranchLengthMean/testBranchLengthVariance,1.0/(testBranchLengthVariance/testBranchLengthMean));
				double p=0;		

				for (int location=0;location<num_locations;location++) {
					p=p+Math.exp(m.logprobability(from.location, location, from.time, to_time));
					if (d<=p) {
						from.children.add(new Node(location,to_time,num_locations));
						break;
					}
				}			
			}

			for (Node child : from.children) {
				makeRandomTree(m,child,(int) Math.round(nNodes/2.0));
			}
		}

	}

	private void removeInternalLocations(Node node) {
		if (node.children.size()!=0) {
			node.location=UNKNOWN_LOCATION;
			for (Node child : node.children) {
				removeInternalLocations(child);				
			}
		}				
	}

	private void clearCachedLikelihood(Node node) {	
		node.cachedConditionalLogLikelihood=new double[num_locations];
		if (node.children.size()!=0) {			
			for (Node child : node.children) {
				clearCachedLikelihood(child);				
			}
		}				
	}

	// Nodes in this tree...
	public class Node implements Serializable {	
		Node parent = null;
		List<Node> children = new ArrayList<Node>();
		int location = 0; // Locations are translated to a discrete integer index 0...num_locations
		double time = 0;
		double[] cachedConditionalLogLikelihood = null;

		protected Node() {};
		
		public Node(int location_, double time_, int num_locations_) {
			location=location_;
			time=time_;
			// TODO: add caching of unknown states...
			// TODO: change to a throw phrase...
			if (location>=num_locations_ && location!=UNKNOWN_LOCATION) {
				System.err.println("error creating node with location: "+location+" out of range: 0<location<"+(num_locations-1));		
			}
			cachedConditionalLogLikelihood = new double[num_locations_+1];
		}
	}

	@Override
	public int getNumLocations() {		
		return num_locations;
	}

	public int getNumIdentifiedLocations() {
		// TODO: this for trees with states...
		return numIdentifiedLocations;
	}

}



