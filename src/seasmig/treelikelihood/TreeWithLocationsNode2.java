package seasmig.treelikelihood;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

//Nodes in this tree...
	@SuppressWarnings("serial")
	public class TreeWithLocationsNode2 implements Serializable {	
		//Node parent = null;
		List<TreeWithLocationsNode2> children = new ArrayList<TreeWithLocationsNode2>();
		int location = 0; // Locations are translated to a discrete integer index 0...num_locations
		double time = 0;
		double[] cachedConditionalLogLikelihood = null;

		protected TreeWithLocationsNode2() {};

		public TreeWithLocationsNode2(int location_, double time_, int num_locations_) {
			location=location_;
			time=time_;
			cachedConditionalLogLikelihood = new double[num_locations_+1];
		}
	}