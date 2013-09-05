package seasmig.treelikelihood;

import java.io.Serializable;


public interface MigrationBaseModel extends Serializable {		
	public double logprobability(int from_location, int to_location, double from_time, double to_time, boolean reverseTime);
	public double[][] transitionMatrix(double from_time, double to_time, boolean reverseTime);
	public int getNumLocations();
	public String print();
	public String parse();
	public String getModelName();
	public double[] rootfreq(double when);
}
