package treelikelihood;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public interface MatrixExponentiator {
	
	// Should return exp(Q*t)
	DoubleMatrix2D expm(double t);

}
