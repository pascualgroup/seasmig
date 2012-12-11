package treelikelihood;
import cern.colt.matrix.tdouble.DoubleMatrix2D;


public interface MigrationBaseModel  {		
	static final int UNKNOWN_STATE = -1;
	public double logprobability(int from_state, int to_state, double from_time, double to_time);
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time);
	public int getNumStates();
	public String print();
}
