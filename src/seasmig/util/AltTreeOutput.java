package seasmig.util;

public class AltTreeOutput {
	String branches="";
	String nodes="";
	public void addNode(int postOrderIndex, double time, int loc, float layout, String seq) {
		if (nodes.length()>0)
			nodes+=",";
		nodes+="{"+Integer.toString(postOrderIndex)+","+String.format("%.3f", time)+","+Integer.toString(loc)+","+String.format("%.3f", layout)+","+seq+"}";			
	}
	public void addBranch(int from, int to) {
		if (nodes.length()>0)
			nodes+=",";
		nodes+="{"+Integer.toString(from)+","+Integer.toString(to)+"}";	
		
	}
}