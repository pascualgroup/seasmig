package seasmig.treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

//Nodes in this tree...
@SuppressWarnings("serial")
public class StochasticMappingNode implements Serializable, Iterable<StochasticMappingNode> {	
	
	StochasticMappingNode parent = null;
	List<StochasticMappingNode> children = new ArrayList<StochasticMappingNode>();
	int loc = TreeWithLocations.UNKNOWN_LOCATION; 
	double time = 0;
	
	protected StochasticMappingNode() {};

	// Internal Node constructor
	public StochasticMappingNode(int loc_, double time_, StochasticMappingNode parent_) {
		loc=loc_;
		time=time_;
		parent=parent_;	
	}

	@Override
	public Iterator<StochasticMappingNode> iterator() {
		return new postOrderIter(this);
	}

	public boolean isTip() {
		return this.children.size()==0;
	}

	public class postOrderIter implements Iterator<StochasticMappingNode>,  Serializable {

		StochasticMappingNode next;
		int nextIndex = 0;

		public postOrderIter(StochasticMappingNode root) {
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
		public StochasticMappingNode next() {
			StochasticMappingNode returnValue=next;
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

	public void addChild(StochasticMappingNode locationTreeNode) {
		this.children.add(locationTreeNode);
		locationTreeNode.parent=this;

	}

	// TODO: (348[&loc=2,rate=1.1478703001047978,states="japan_korea"]:2.44, ....
	public String parseState(double[] rootfreq) {			
		return "[&loc="+loc+"]";		
	}
	
	


}