package treelikelihood;
import cern.colt.matrix.tdouble.DoubleMatrix2D;


public interface MigrationBaseModel  {		
	static final int UNKNOWN_LOCATION = -1;
	public double logprobability(int from_location, int to_location, double from_time, double to_time);
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time);
	public int getNumLocations();
	public String print();
	public String parse();
}
