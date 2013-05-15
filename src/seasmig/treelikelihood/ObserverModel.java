package seasmig.treelikelihood;

public interface ObserverModel {
	double logObservationProbability(int location, double time);
}
