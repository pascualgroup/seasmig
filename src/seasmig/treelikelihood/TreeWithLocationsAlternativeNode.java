package seasmig.treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;

//Nodes in this tree...
@SuppressWarnings("serial")
public class TreeWithLocationsAlternativeNode implements Serializable, Iterable<TreeWithLocationsAlternativeNode> {	

	public static final int UNKNOWN_LOCATION = -1;
	public double[] logprobs;
	TreeWithLocationsAlternativeNode parent = null;
	List<TreeWithLocationsAlternativeNode> children = new ArrayList<TreeWithLocationsAlternativeNode>();
	int loc = UNKNOWN_LOCATION; 
	double time = 0;

	protected TreeWithLocationsAlternativeNode() {};

	// Internal Node constructor
	public TreeWithLocationsAlternativeNode(int loc_, double time_, TreeWithLocationsAlternativeNode parent_) {
		loc=loc_;
		time=time_;
		parent=parent_;	
	}

	@Override
	public Iterator<TreeWithLocationsAlternativeNode> iterator() {
		return new postOrderIter(this);
	}
	
	public boolean isTip() {
		return this.children.size()==0;
	}

	public class postOrderIter implements Iterator<TreeWithLocationsAlternativeNode>,  Serializable {

		TreeWithLocationsAlternativeNode next;
		int nextIndex = 0;

		public postOrderIter(TreeWithLocationsAlternativeNode root) {
			// Start at left most leaf
			if (root.parent!=null)
				nextIndex=position(root, root.parent.children);
			next=root;
			while (next.children.size()>0) {
				next=next.children.get(0);
			}
		}

		@Override
		public boolean hasNext() { 
			return (next!=null); 
		}

		private int position(TreeWithLocationsAlternativeNode node, List<TreeWithLocationsAlternativeNode> list) {
			for (int i=0;i<list.size();i++) {
				if (list.get(i)==node) {
					return i;
				}
			}
			return list.size();
		}

		@Override
		public TreeWithLocationsAlternativeNode next() {
			TreeWithLocationsAlternativeNode returnValue=next;
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
							nextIndex=position(next, next.parent.children);
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
					next.parent.children.remove(position(next,next.parent.children));
					nextIndex-=1;
				}
			}				
		}
	}

	public void addChild(TreeWithLocationsAlternativeNode locationTreeNode) {
		this.children.add(locationTreeNode);
		locationTreeNode.parent=this;
		
	}
	
}