package seasmig.util;

public class AltTreeOutput {
	private String branches="";
	private String nodes="";
	public void addNode(int postOrderIndex, double time, int loc, float layout, String seq) {
		if (getNodes().length()>0)
			setNodes(getNodes() + ",");
		setNodes(getNodes() + "{"+Integer.toString(postOrderIndex)+","+String.format("%.3f", time)+","+Integer.toString(loc)+","+String.format("%.3f", layout)+","+seq+"}");			
	}
	public void addBranch(int from, int to) {
		if (getNodes().length()>0)
			setNodes(getNodes() + ",");
		setNodes(getNodes() + "{"+Integer.toString(from)+","+Integer.toString(to)+"}");	
		
	}
	public String getNodes() {
		return nodes;
	}
	public void setNodes(String nodes) {
		this.nodes = nodes;
	}
	public String getBranches() {
		return branches;
	}
	public void setBranches(String branches) {
		this.branches = branches;
	}
}