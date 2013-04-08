package seasmig.treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

//Nodes in this tree...
	public class TreeWithLocationsNode implements Serializable {	
		//Node parent = null;
		List<TreeWithLocationsNode> children = new ArrayList<TreeWithLocationsNode>();
		int location = 0; // Locations are translated to a discrete integer index 0...num_locations
		double time = 0;
		double[] cachedConditionalLogLikelihood = null;

		protected TreeWithLocationsNode() {};

		public TreeWithLocationsNode(int location_, double time_, int num_locations_) {
			location=location_;
			time=time_;
			cachedConditionalLogLikelihood = new double[num_locations_+1];
		}
	}