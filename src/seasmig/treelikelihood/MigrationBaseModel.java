package seasmig.treelikelihood;

import java.io.Serializable;


public interface MigrationBaseModel extends Serializable {		
	public double logprobability(int from_location, int to_location, double from_time, double to_time);
	public double[][] transitionMatrix(double from_time, double to_time);
	public int getNumLocations();
	public String print();
	public String parse();
	public String getModelName();
	public double[] probability(int from_location, double from_time, double to_time);
	public double[] rootfreq(double when);
	public Event getNextRandomEvent(double time, int loc);
	
	class Event {
		public double time;
		int loc; 
	}
}
