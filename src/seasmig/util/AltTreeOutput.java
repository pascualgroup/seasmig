package seasmig.util;

import seasmig.treelikelihood.trees.TreeWithLocationsNode.NodeType;

public class AltTreeOutput {
//	private String branches="";
	private String nodes="";

	public void addNode(int postOrderIndex, double time, int loc, float layout, String seq, NodeType nodeType, int taxonIndex, int parent) {
		if (nodes.length()>0)
			nodes+=",";		
		nodes+="{"+"\"index\":"+Integer.toString(postOrderIndex)+","+"\"time\":"+String.format("%.3f", time)+","+"\"location\":"+Integer.toString(loc)+","+"\"layout\":"+String.format("%.3f", layout)+","+"\"seq\":"+seq+","+"\"type\":"+nodeType+","+"\"taxon\":"+Integer.toString(taxonIndex)+","+"\"parent\":"+Integer.toString(parent)+"}";			
	}
//	public void addBranch(int from, int to) {
//		if (branches.length()>0)
//			branches+=",";
//		branches+="{"+Integer.toString(from)+","+Integer.toString(to)+"}";	
//		
//	}
	public String getNodes() {
		return nodes;
	}

//	public String getBranches() {
//		return branches;
//	}
	
}