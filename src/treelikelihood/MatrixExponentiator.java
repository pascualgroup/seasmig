package treelikelihood;

public interface MatrixExponentiator {
	
	// Should return exp(Q*t)
	double[][] expm(double t);

}
