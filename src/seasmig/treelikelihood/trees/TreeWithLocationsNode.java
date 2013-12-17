package seasmig.treelikelihood.trees;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import seasmig.treelikelihood.TransitionModel.Transition;

//Nodes in this tree...
@SuppressWarnings("serial")
public class TreeWithLocationsNode implements Serializable, Iterable<TreeWithLocationsNode> {	

	public double[] logProbs;
	TreeWithLocationsNode parent = null;
	List<TreeWithLocationsNode> children = new ArrayList<TreeWithLocationsNode>();
	double time = 0;
	int taxonIndex = 0;

	// For tips and for ASR
	int loc = TreeWithLocations.UNKNOWN_LOCATION;

	// For stochastic mapping
	public List<Transition> transitions = null;
	private boolean isTrunk = false;

	public static final double minNegative = Double.NEGATIVE_INFINITY;

	protected TreeWithLocationsNode() {};

	// Internal Node constructor
	public TreeWithLocationsNode(int loc_, int taxonIndex_, double time_, TreeWithLocationsNode parent_) {
		loc=loc_;
		time=time_;
		parent=parent_;	
		taxonIndex=taxonIndex_;
	}

	@Override
	public Iterator<TreeWithLocationsNode> iterator() {
		return new postOrderIter(this);
	}

	public boolean isTip() {
		return this.children.size()==0;
	}
	
	public boolean isTrunk() {
		return isTrunk;
	}
	
	public void setTrunk() {
		isTrunk=true;;
	}

	public class postOrderIter implements Iterator<TreeWithLocationsNode>,  Serializable {

		TreeWithLocationsNode next;
		int nextIndex = 0;

		public postOrderIter(TreeWithLocationsNode root) {
			// Start at left most leaf
			if (root.parent!=null) {
				nextIndex=(root.parent.children.get(0)==root ? 0 : 1);
			}
			next=root;
			while (next.children.size()>0) {
				next=next.children.get(0);
			}
		}

		@Override
		public boolean hasNext() { 
			return (next!=null); 
		}

		@Override
		public TreeWithLocationsNode next() {
			TreeWithLocationsNode returnValue=next;
			// If more children after this one reach leftmost child of next child
			if (next.parent!=null) {
				if ((nextIndex+1)<next.parent.children.size()){
					nextIndex+=1;
					next=next.parent.children.get(nextIndex);
					while (next.children.size()>0) {
						next=next.children.get(0);
						nextIndex=0;
					}
				}
				else {
					next=next.parent;
					if (next!=null) {
						if (next.parent!=null)
							nextIndex = (next.parent.children.get(0)==next ? 0 : 1);
						else
							nextIndex=0;
					}
				}
			}
			else {
				next=null;
			}
			return returnValue;
		}

		@Override
		public void remove() {
			if (next!=null) {		
				if (next.parent!=null) {
					next.parent.children.remove((next.parent.children.get(0)==next ? 0 : 1)/*positionBinaryTree(next,next.parent.children)*/);
					nextIndex-=1;
				}
			}				
		}
	}


	public void addChild(TreeWithLocationsNode locationTreeNode) {
		this.children.add(locationTreeNode);
		locationTreeNode.parent=this;

	}

	public int getTaxonIndex() {		
		return taxonIndex;
	}

	// TODO: (348[&antigenic={-6.00510611736,5.84199000915},rate=1.1478703001047978,states="japan_korea"]:2.44, ....
	public String parseProbs(double[] rootfreq) {			
		String returnValue = "[&prob={";		

		double[] logAncStateProbs = new double[logProbs.length];
		double[] logAncProbs = new double[logProbs.length];

		for (int i=0;i<logProbs.length;i++) {	
			if (parent!=null)
				logAncStateProbs[i] = logProbs[i]+parent.logProbs[i];
			else
				logAncStateProbs[i] = logProbs[i]+rootfreq[i];					
		}

		double logTotalProb = logSumExp(logAncStateProbs);
		for (int i=0;i<logProbs.length;i++) {
			logAncStateProbs[i]=logAncStateProbs[i]-logTotalProb;
		}

		// Deal with numeric issues of not summing up to 1.0
		double totalProb = 0;
		for (int i=0;i<logProbs.length;i++) {			
			logAncProbs[i]=Math.exp(logAncStateProbs[i]);
			totalProb+=logAncProbs[i];
		}

		for (int i=0;i<logProbs.length;i++) {
			String prob = String.format("%.3f",logAncProbs[i]/totalProb);
			returnValue+=prob;
			if (i!=(logProbs.length-1)) 
				returnValue+=",";
		}
		returnValue+="}]";
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
			sumExp=sumExp+cern.jet.math.Functions.exp.apply(alphas[i]-minWithoutNegInf);
		}
		double returnValue=minWithoutNegInf+cern.jet.math.Functions.log.apply(sumExp);
		if (Double.isNaN(returnValue) ){
			return minNegative;
		}
		return returnValue;
	}

	public String parseMap() {
		String returnValue = new String();

		if (transitions!=null) {
			if (transitions.size()>0) {	
				returnValue+="[&map={";
				double timeFrom = parent.time;

				for (int i=0;i<transitions.size();i++) {
					returnValue+=String.format("%.3f", transitions.get(i).time-timeFrom);						
					if (i!=(transitions.size()-1)) {
						returnValue+=","+Integer.toString(transitions.get(i).loc)+",";
					}
					timeFrom = transitions.get(i).time;
				}
				returnValue+="}]";
			}	
		}
		return returnValue;
	}


}