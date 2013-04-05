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

	public class postOrderIter implements Iterator<LocationTreeNode> {

		LocationTreeNode current;
		int currentIndex = 0;

		public postOrderIter(LocationTreeNode root) {
			// Start at left most leaf
			if (root.parent!=null)
				currentIndex=position(root, root.parent.children);
			current=root;
			while (current.children.size()>0) {
				current=current.children.get(0);
			}
		}

		@Override
		public boolean hasNext() { 
			return (current!=null); 
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
			LocationTreeNode returnValue=current;
			// If more children after this one reach leftmost child of next child
			if (current.parent!=null) {
				if ((currentIndex+1)<current.parent.children.size()){
					currentIndex+=1;
					current=current.parent.children.get(currentIndex);
					while (current.children.size()>0) {
						current=current.children.get(0);
						currentIndex=0;
					}
				}
				else {
					current=current.parent;
					if (current!=null) {
						if (current.parent!=null)
							currentIndex=position(current, current.parent.children);
						else
							currentIndex=0;
					}
				}
			}
			else {
				current=null;
			}
			return returnValue;
		}

		@Override
		public void remove() {
			// TODO Test This!!!!
			LocationTreeNode toRemove=current;
			// If more children after this one reach leftmost child of next child
			if (current.parent!=null) {				
				if ((currentIndex+1)<current.parent.children.size()){
					currentIndex+=1;
					current=current.parent.children.get(currentIndex);
					while (current.children.size()>0) {
						current=current.children.get(0);
						currentIndex=0;
					}
				}
				else {
					current=current.parent;
					if (current!=null) {
						if (current.parent!=null)
							currentIndex=position(current, current.parent.children);
						else
							currentIndex=0;						
					}

				}
				toRemove.parent.children.remove(position(toRemove,toRemove.parent.children));
			}
			else {
				current=null;
			}		
		}
	}
}