package seasmig.treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;

//Nodes in this tree...
@SuppressWarnings("serial")
public class LocationTreeNode implements Serializable, Iterable<LocationTreeNode> {	

	public static final int UNKNOWN_LOCATION = -1;
	public double[] logprobs;
	LocationTreeNode parent = null;
	List<LocationTreeNode> children = new ArrayList<LocationTreeNode>();
	BitSet cladeTaxaIndices = new BitSet(); // Stores the taxa indices of all children  
	int loc = UNKNOWN_LOCATION; 
	double time = 0;

	protected LocationTreeNode() {};

	// Internal Node constructor
	public LocationTreeNode(int loc_, double time_, LocationTreeNode parent_) {
		loc=loc_;
		time=time_;
		parent=parent_;	
	}

	// Tip Constructor
	public LocationTreeNode(int loc_, int taxaIndex, double time_, LocationTreeNode parent_) {
		loc=loc_;
		time=time_;
		parent=parent_;	
		cladeTaxaIndices.set(taxaIndex);
	}
	
	public BitSet fillCladeTaxaIndices() {
		for (LocationTreeNode child : children) {
			cladeTaxaIndices.or(child.fillCladeTaxaIndices());
		}
		return cladeTaxaIndices;
	}

	@Override
	public Iterator<LocationTreeNode> iterator() {
		return new postOrderIter(this);
	}
	
	public boolean isTip() {
		return this.children.size()==0;
	}

	public class postOrderIter implements Iterator<LocationTreeNode>,  Serializable {

		LocationTreeNode next;
		int nextIndex = 0;

		public postOrderIter(LocationTreeNode root) {
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

		private int position(LocationTreeNode node, List<LocationTreeNode> list) {
			for (int i=0;i<list.size();i++) {
				if (list.get(i)==node) {
					return i;
				}
			}
			return list.size();
		}

		@Override
		public LocationTreeNode next() {
			LocationTreeNode returnValue=next;
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
	
}