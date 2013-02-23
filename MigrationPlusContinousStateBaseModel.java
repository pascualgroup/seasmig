package treelikelihood;


public interface MigrationPlusContinousStateBaseModel  {				
	public String print();
	double logprobability(int from_location, int to_location, double from_time, double to_time, double from_state, double to_state);
	// TODO: maybe add get kernel function
}
