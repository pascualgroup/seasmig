package seasmig.treelikelihood;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


public interface MigrationBaseModel extends Serializable {
	
	class Event {
		Event(double time, int loc) {
			this.time = time;
			this.loc = loc;
		}
		public double time;
		public int loc; 
	}
	
	
	public double logprobability(int from_location, int to_location, double from_time, double to_time);
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time);
	public int getNumLocations();
	public String print();
	public String parse();
	public String getModelName();
	public DoubleMatrix1D probability(int from_location, double from_time, double to_time);
	public DoubleMatrix1D rootfreq(double when);
	public Event nextEvent(double time, int loc);
	
}
