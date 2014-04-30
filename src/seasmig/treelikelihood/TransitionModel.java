package seasmig.treelikelihood;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


public interface TransitionModel extends Serializable {
	
	@SuppressWarnings("serial")
	public class Transition implements Serializable {
		
		protected Transition() {};
		public Transition(double time, int toTrait) {
			this.time = time;
			this.toTrait = toTrait;
		}
		public double time;
		public int toTrait; 
	}
	
	
	public double logprobability(int from_trait, int to_trait, double from_time, double to_time);
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time);
	public int getNumLocations();
	public String print();
	public String parse();
	public String getModelName();
	public DoubleMatrix1D probability(int from_trait, double from_time, double to_time);
	public DoubleMatrix1D rootfreq(double when);
	public Transition nextEvent(double time, int trait);
	
}
