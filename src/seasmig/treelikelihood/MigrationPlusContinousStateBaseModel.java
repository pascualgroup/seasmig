package seasmig.treelikelihood;

public interface MigrationPlusContinousStateBaseModel {

	double logprobability(int from_location, int to_location, double from_time,
			double to_time, double from_state, double to_state);

	String print();

}
