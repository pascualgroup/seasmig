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
		int state = 0;	
		double time = 0;
		Double[] cachedConditionalLogLikelihood = null;
		
		public Node(int trait_, double time_, int num_states_) {
			state=trait_;
			time=time_;
			cachedConditionalLogLikelihood = new Double[num_states_];
		}
	}

	Node root = null;		
	int num_states = 0;
	
	public Tree(MigrationBaseModel createTreeModel, int numNodes) {
		num_states=createTreeModel.getNumStates();
		root = new Node(0,0,num_states);
		makeRandomTree(createTreeModel, root, numNodes);		
	}
	
	private Tree() {
	}
	
	public Tree(jebl.evolution.trees.SimpleRootedTree tree, int num_states_) {
		num_states=num_states_;
		root = new Node(Integer.parseInt((String)tree.getRootNode().getAttribute("states"))-1,0,num_states);
		makeSubTree(tree,root,tree.getRootNode());
	}
	
	public Tree(jebl.evolution.trees.SimpleRootedTree tree, HashMap<String, Integer> traitMap, int num_states_) {
		num_states=num_states_;
		Integer trait = traitMap.get(tree.getTaxon(tree.getRootNode()));
        if (trait==null) 
        	trait=-1;
		root = new Node(trait,0,num_states);
		makeSubTree(tree,traitMap,root,tree.getRootNode());
	}
	
	private void makeSubTree(SimpleRootedTree inputTree,
			HashMap<String, Integer> traitMap, Node outputSubTree,
			jebl.evolution.graphs.Node inputSubTree) {
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {	
			Taxon taxon = inputTree.getTaxon(node);
			Integer trait = MigrationBaseModel.UNKNOWN_STATE;
			if (taxon!=null) {
				trait = traitMap.get(inputTree.getTaxon(node).toString());
				if (trait==null) 
					trait=-1;
			}			
			outputSubTree.children.add(new Node(trait,outputSubTree.time+inputTree.getLength(node),num_states));
			makeSubTree(inputTree,traitMap, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}
		
	}

	public void makeSubTree(jebl.evolution.trees.SimpleRootedTree inputTree, Node outputSubTree, jebl.evolution.graphs.Node inputSubTree) {
		
		for (jebl.evolution.graphs.Node node : inputTree.getChildren(inputSubTree)) {						
			outputSubTree.children.add(new Node(Integer.parseInt((String)node.getAttribute("states"))-1,outputSubTree.time+inputTree.getLength(node),num_states));
			makeSubTree(inputTree, outputSubTree.children.get(outputSubTree.children.size()-1), node);			
		}
		
	}
	
	public void makeRandomTree(MigrationBaseModel m, Node from, int nNodes) {		
		if (nNodes>1) {
			for (int child=0;child<2;child++) {
				double d = cern.jet.random.Uniform.staticNextDouble();
				double to_time = from.time+cern.jet.random.Gamma.staticNextDouble(2, 2);
				double p=0;		
				
				for (int state=0;state<num_states;state++) {
					p=p+Math.exp(m.logprobability(from.state, state, from.time, to_time));
					if (d<=p) {
						from.children.add(new Node(state,to_time,num_states));
						break;
					}
				}			
			}

			for (Node child : from.children) {
				makeRandomTree(m,child,(int) Math.round(nNodes/2.0));
			}
		}

	}
		

	private void removeInternalStates(Node node) {
		if (node.children.size()!=0) {
			node.state=MigrationBaseModel.UNKNOWN_STATE;
			for (Node child : node.children) {
				removeInternalStates(child);				
			}
		}				
	}
	
	private void clearCachedLikelihood(Node node) {	
		node.cachedConditionalLogLikelihood=new Double[num_states];
		if (node.children.size()!=0) {			
			for (Node child : node.children) {
				clearCachedLikelihood(child);				
			}
		}				
	}

	public void print(Node treePart) {

		System.out.print(treePart.state);
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
		double[] alphas=new double[num_states];
		double min = Double.MIN_VALUE;
		if (root.state==MigrationBaseModel.UNKNOWN_STATE) {
			for (int rootState=0;rootState<num_states;rootState++) {				
				double alpha=conditionalLogLikelihood(likelihoodModel,root, rootState);
				alphas[rootState]=alpha;
				if (alpha<min) min=alpha;				
			}
			return logSumExp(alphas,min);
		}
		else {
			return conditionalLogLikelihood(likelihoodModel,root,root.state);
		}		
	}

	private double conditionalLogLikelihood(MigrationBaseModel likelihoodModel, Node node, int nodeState) {

		if (node.cachedConditionalLogLikelihood[nodeState]!=null) {			
			return node.cachedConditionalLogLikelihood[nodeState];
		}
		else {
			double loglikelihood=0;
			for (Node child : node.children) {
				if (child.state!=MigrationBaseModel.UNKNOWN_STATE) {
					assert child.time>node.time;
					loglikelihood=loglikelihood+conditionalLogLikelihood(likelihoodModel,child,child.state)+likelihoodModel.logprobability(nodeState, child.state, node.time, child.time);
				}
				else {
					double[] alphas=new double[num_states];
					double min = Double.MIN_VALUE;
					for (int childState=0;childState<num_states;childState++) {
						double alpha = likelihoodModel.logprobability(nodeState, childState, node.time, child.time)+conditionalLogLikelihood(likelihoodModel,child, childState);						
						alphas[childState]=alpha;
						if (alpha<min) min=alpha;
					}
					loglikelihood=loglikelihood+logSumExp(alphas,min);
				}			 
			}
			node.cachedConditionalLogLikelihood[nodeState]=loglikelihood;
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

	public void removeInternalStates() {
		removeInternalStates(root);
		
	}

	public void clearCachedLikelihood() {
		clearCachedLikelihood(root);	
	}
	
	public Tree copyWithNoCache() {
		Tree newTree = new Tree();
		newTree.num_states=num_states;
		newTree.root=new Node(root.state,root.time,num_states);
		treeCopyNoCache(root,newTree.root);		
		return newTree;
	}

	private void treeCopyNoCache(Node from, Node to) {
		for (Node child : from.children) {
			Node newChild = new Node(child.state,child.time,num_states);
			newChild.parent=from;
			to.children.add(newChild);			
			treeCopyNoCache(child, newChild);
		}		
	}
}



