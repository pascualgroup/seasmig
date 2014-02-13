package seasmig.util;

public class AltTreeOutput {
	private String branches="";
	private String nodes="";
	public void addNode(int postOrderIndex, double time, int loc, float layout, String seq, boolean isTip, int taxonIndex) {
		if (nodes.length()>0)
			nodes+=",";
		nodes+="{"+Integer.toString(postOrderIndex)+","+String.format("%.3f", time)+","+Integer.toString(loc)+","+String.format("%.3f", layout)+","+seq+","+isTip+","+Integer.toString(taxonIndex)+"}";			
	}
	public void addBranch(int from, int to) {
		if (branches.length()>0)
			branches+=",";
		branches+="{"+Integer.toString(from)+","+Integer.toString(to)+"}";	
		
	}
	public String getNodes() {
		return nodes;
	}

	public String getBranches() {
		return branches;
	}
	
}