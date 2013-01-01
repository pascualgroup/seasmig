package treelikelihood;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;


public class Tree {
	
	public class Node {	
		Node parent = null;
		List<Node> children = new ArrayList<Node>();
		int location = 0;	
		//double state = 0;
		double time = 0;
		double[] cachedConditionalLogLikelihood = null;
		
		public Node(int location_, double time_, int num_locations_) {
			location=location_;
			time=time_;
			cachedConditionalLogLikelihood = new double[num_locations_];
		}
	}

	Node root = null;		
	int num_locations = 0;
	
	// TODO: organize location and state names and move to config...
	String location_attribute_name = "states";
	String state_attribute_name = null;
	
	public Tree(MigrationBaseModel createTreeModel, int numNodes) {
		num_locations=createTreeModel.getNumLocations();
		root = new Node(0,0,num_locations);
		makeRandomTree(createTreeModel, root, numNodes);		
	}
	
	private Tree() {
	}
	
	public Tree(jebl.evolution.trees.SimpleRootedTree tree, int num_locations_) {
		num_locations=num_locations_;
		root = new Node(Integer.parseInt((String)tree.getRootNode().getAttribute(location_attribute_name))-1,0,num_locations);
		makeSubTree(tree,root,tree.getRootNode());
	}

	// TODO: add load states from tree and not just file....
	public Tree(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> locationMap, int num_locations_/*, HashMap<String, Double> stateMap*/) {
		num_locations=num_locations_;
		Integer location = locationMap.get(tree.getTaxon(tree.getRootNode()));
		if (location==null) 
        	location=-1;
		//Double state = stateMap.get(tree.getTaxon(tree.getRootNode()));
        //if (state==null) 
        	
		root = new Node(location,0,num_locations);
		makeSubTree(tree,locationMap,root,tree.getRootNode());
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
			outputSubTree.children.add(new Node(location,outputSubTree.time+inputTree.getLength(node),num_locations));
			makeSubTree(inputTree,locationMap, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}
		
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, Node outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {						
			outputSubTree.children.add(new Node(Integer.parseInt((String)node.getAttribute(location_attribute_name))-1,outputSubTree.time+inputTree.getLength(node),num_locations));
			makeSubTree(inputTree, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}
		
	}
	
	public void makeRandomTree(MigrationBaseModel m, Node from, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				double to_time = from.time+cern.jet.random.Gamma.staticNextDouble(2, 2);
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
			node.location=MigrationBaseModel.UNKNOWN_LOCATION;
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

	public void print(Node treePart) {

		System.out.print(treePart.location);
		if (treePart.children.size()>0) {
			System.out.print(" (");
			print(treePart.children.get(0));
			if (treePart.children.size()>1) {
				for (int i=1;i<treePart.children.size();i++) {
					System.out.print(",");
					print(treePart.children.get(i));					
				}				
			}
			System.out.print(")");
		}
	}


	public double logLikelihood(MigrationBaseModel likelihoodModel) {
		double[] alphas=new double[num_locations];
		double min = Double.MIN_VALUE;
		if (root.location==MigrationBaseModel.UNKNOWN_LOCATION) {
			for (int rootLocation=0;rootLocation<num_locations;rootLocation++) {				
				double alpha=conditionalLogLikelihood(likelihoodModel,root, rootLocation);
				alphas[rootLocation]=alpha;
				if (alpha<min) min=alpha;				
			}
			return logSumExp(alphas,min);
		}
		else {
			return conditionalLogLikelihood(likelihoodModel,root,root.location);
		}		
	}

	private double conditionalLogLikelihood(MigrationBaseModel likelihoodModel, Node node, int nodeLocation) {

		if (node.cachedConditionalLogLikelihood[nodeLocation]!=0) {			
			return node.cachedConditionalLogLikelihood[nodeLocation];
		}
		else {
			double loglikelihood=0;
			for (Node child : node.children) {
				if (child.location!=MigrationBaseModel.UNKNOWN_LOCATION) {
					assert child.time>node.time;
					loglikelihood=loglikelihood+conditionalLogLikelihood(likelihoodModel,child,child.location)+likelihoodModel.logprobability(nodeLocation, child.location, node.time, child.time);
				}
				else {
					double[] alphas=new double[num_locations];
					double min = Double.MIN_VALUE;
					for (int childLocation=0;childLocation<num_locations;childLocation++) {
						double alpha = likelihoodModel.logprobability(nodeLocation, childLocation, node.time, child.time)+conditionalLogLikelihood(likelihoodModel,child, childLocation);						
						alphas[childLocation]=alpha;
						if (alpha<min) min=alpha;
					}
					loglikelihood=loglikelihood+logSumExp(alphas,min);
				}			 
			}
			node.cachedConditionalLogLikelihood[nodeLocation]=loglikelihood;
			return loglikelihood;		
		}
	}

	static double logSumExp(double[] alphas,double min) {

		if (min>Double.NEGATIVE_INFINITY) {
			double sumExp = 0;
			for (int i=0;i<alphas.length;i++) {			
				sumExp=sumExp+Math.exp(alphas[i]-min);
			}
			return min+Math.log(sumExp);
		}
		else {
			double sumExp = 0;
			for (int i=0;i<alphas.length;i++) {			
				sumExp=sumExp+Math.exp(alphas[i]);
			}
			return Math.log(sumExp);
		}
	}

	public void print() {
		print(root);
		
	}

	public void removeInternalLocations() {
		removeInternalLocations(root);
		
	}

	public void clearCachedLikelihood() {
		clearCachedLikelihood(root);	
	}
	
	public Tree copyWithNoCache() {
		Tree newTree = new Tree();
		newTree.num_locations=num_locations;
		newTree.root=new Node(root.location,root.time,num_locations);
		treeCopyNoCache(root,newTree.root);		
		return newTree;
	}

	private void treeCopyNoCache(Node from, Node to) {
		for (Node child : from.children) {
			Node newChild = new Node(child.location,child.time,num_locations);
			newChild.parent=from;
			to.children.add(newChild);			
			treeCopyNoCache(child, newChild);
		}		
	}
}



