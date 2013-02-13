package treelikelihood;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;


public class TreeWithLocationsAlternative implements LikelihoodTree {

	// Tree & Model
	Node root = null;		
	int num_locations = 0;
	private MigrationBaseModel likelihoodModel = null;


	// Load a tree from a basic jebl tree
	// locations are loaded from nexsus tree trait location_attribute name
	public TreeWithLocationsAlternative(jebl.evolution.trees.SimpleRootedTree tree, String locationAttributeName, int num_locations_) {
		num_locations=num_locations_;
		root = new Node(Integer.parseInt((String)tree.getRootNode().getAttribute(locationAttributeName))-1,0,num_locations);
		makeSubTree(tree,locationAttributeName, root,tree.getRootNode());
	}

	// Load a tree from a basic jebl tree
	// locations are loaded from a hashmap	
	public TreeWithLocationsAlternative(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
		num_locations=num_locations_;
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
			location=-1;
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
		double min = Util.minValue;
		if (root.location==MigrationBaseModel.UNKNOWN_LOCATION) {
			for (int rootLocation=0;rootLocation<num_locations;rootLocation++) {				
				double alpha=conditionalLogLikelihood(root, rootLocation);
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
		TreeWithLocationsAlternative newTree = new TreeWithLocationsAlternative();
		newTree.num_locations=num_locations;
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
				if (child.location!=MigrationBaseModel.UNKNOWN_LOCATION) {
					assert child.time>node.time;
					loglikelihood=loglikelihood+conditionalLogLikelihood(child,child.location)+likelihoodModel.logprobability(nodeLocation, child.location, node.time, child.time);
				}
				else {
					double[] alphas=new double[num_locations];
					double min = Double.POSITIVE_INFINITY;
					for (int childLocation=0;childLocation<num_locations;childLocation++) {
						double alpha = likelihoodModel.logprobability(nodeLocation, childLocation, node.time, child.time)+conditionalLogLikelihood(child, childLocation);						
						alphas[childLocation]=alpha;
						if (alpha<min) min=alpha;
					}
					loglikelihood=loglikelihood+ Util.logSumExp(alphas,min);
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

	private TreeWithLocationsAlternative() {
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
			Integer location = MigrationBaseModel.UNKNOWN_LOCATION;
			if (taxon!=null) {
				location = locationMap.get(inputTree.getTaxon(node).toString());
				if (location==null) 
					location=-1;
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

	// Nodes in this tree...
	public class Node {	
		Node parent = null;
		List<Node> children = new ArrayList<Node>();
		int location = 0; // Locations are translated to a discrete integer index 0...num_locations
		double time = 0;
		double[] cachedConditionalLogLikelihood = null;

		public Node(int location_, double time_, int num_locations_) {
			location=location_;
			time=time_;
			// TODO: add caching of unknown states...
			// TODO: change to a throw phrase...
			if (location>=num_locations_) {
				System.err.println("error creating node with location: "+location+" out of range: 0<location<"+(num_locations-1));		
			}
			cachedConditionalLogLikelihood = new double[num_locations_];
		}
	}


}



