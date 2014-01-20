package seasmig.treelikelihood.trees;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import seasmig.treelikelihood.TransitionModel.Transition;

//Nodes in this tree...
@SuppressWarnings("serial")
public class TreeWithLocationsNode implements Serializable, Iterable<TreeWithLocationsNode> {	

	private TreeWithLocationsNode parent = null;
	public List<TreeWithLocationsNode> children = new ArrayList<TreeWithLocationsNode>();
	double time = 0;
	int taxonIndex = 0;

	private int loc = TreeWithLocations.UNKNOWN_LOCATION;
	public Sequence seq;
	
	public double[] logProbsLOC;
	public double[][][] logProbsCP; // codon position 0,1,2; sequence position TODO: 

	// For stochastic mapping
	public List<Transition> migrations = null;
	public ArrayList<ArrayList<ArrayList<Transition>>> mutations = null;
	
	private boolean isTrunk = false;
	private Sequence noSequence = new Sequence(0);

	public static final double minNegative = Double.NEGATIVE_INFINITY;

	protected TreeWithLocationsNode() {};

	// Internal Node constructor
	public TreeWithLocationsNode(Sequence seq_, int loc_, int taxonIndex_, double time_, TreeWithLocationsNode parent_, boolean useSequenceData) {
		if (seq==null) {
			seq = noSequence;
		} 	
		if (useSequenceData)
			if (!seq.isTip()) seq=seq_.copy();		
		setLoc(loc_);
		time=time_;
		setParent(parent_);	
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
			if (root.getParent()!=null) {
				nextIndex=(root.getParent().children.get(0)==root ? 0 : 1);
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
			if (next.getParent()!=null) {
				if ((nextIndex+1)<next.getParent().children.size()){
					nextIndex+=1;
					next=next.getParent().children.get(nextIndex);
					while (next.children.size()>0) {
						next=next.children.get(0);
						nextIndex=0;
					}
				}
				else {
					next=next.getParent();
					if (next!=null) {
						if (next.getParent()!=null)
							nextIndex = (next.getParent().children.get(0)==next ? 0 : 1);
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
				if (next.getParent()!=null) {
					next.getParent().children.remove((next.getParent().children.get(0)==next ? 0 : 1)/*positionBinaryTree(next,next.parent.children)*/);
					nextIndex-=1;
				}
			}				
		}
	}


	public void addChild(TreeWithLocationsNode locationTreeNode) {
		this.children.add(locationTreeNode);
		locationTreeNode.setParent(this);

	}

	public int getTaxonIndex() {		
		return taxonIndex;
	}

	// TODO: (348[&antigenic={-6.00510611736,5.84199000915},rate=1.1478703001047978,states="japan_korea"]:2.44, ....
	public String parseProbs(double[] rootfreq) {			
		String returnValue = "[&prob={";		

		double[] logAncStateProbs = new double[logProbsLOC.length];
		double[] logAncProbs = new double[logProbsLOC.length];

		for (int i=0;i<logProbsLOC.length;i++) {	
			if (getParent()!=null)
				logAncStateProbs[i] = logProbsLOC[i]+getParent().logProbsLOC[i];
			else
				logAncStateProbs[i] = logProbsLOC[i]+rootfreq[i];					
		}

		double logTotalProb = logSumExp(logAncStateProbs);
		for (int i=0;i<logProbsLOC.length;i++) {
			logAncStateProbs[i]=logAncStateProbs[i]-logTotalProb;
		}

		// Deal with numeric issues of not summing up to 1.0
		double totalProb = 0;
		for (int i=0;i<logProbsLOC.length;i++) {			
			logAncProbs[i]=Math.exp(logAncStateProbs[i]);
			totalProb+=logAncProbs[i];
		}

		for (int i=0;i<logProbsLOC.length;i++) {
			String prob = String.format("%.3f",logAncProbs[i]/totalProb);
			returnValue+=prob;
			if (i!=(logProbsLOC.length-1)) 
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

		if (migrations!=null) {
			if (migrations.size()>0) {	
				returnValue+="[&map={";
				double timeFrom = getParent().time;

				for (int i=0;i<migrations.size();i++) {
					returnValue+=String.format("%.3f", migrations.get(i).time-timeFrom);						
					if (i!=(migrations.size()-1)) {
						returnValue+=","+Integer.toString(migrations.get(i).toTrait)+",";
					}
					timeFrom = migrations.get(i).time;
				}
				returnValue+="}]";
			}	
		}
		return returnValue;
	}

	public TreeWithLocationsNode left() {
		if (this.children!=null) 
			if (this.children.size()>0) 
				return this.children.get(0);
		return null;			
	}
	
	public TreeWithLocationsNode right() {
		if (this.children!=null) 
			if (this.children.size()>1) 
				return this.children.get(1);
		return null;			
	}

	public int getLoc() {
		return loc;
	}

	public void setLoc(int loc) {
		this.loc = loc;
	}

	public TreeWithLocationsNode getParent() {
		return parent;
	}

	public void setParent(TreeWithLocationsNode parent) {
		this.parent = parent;
	}
}