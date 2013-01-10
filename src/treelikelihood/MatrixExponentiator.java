package treelikelihood;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public interface MatrixExponentiator {
	
	// TODO: add closed form implementation
	// Compute Taylor expansion to calculate P(t)=Exp(Qt) 
	
	// Should return exp(Q*t)
	DoubleMatrix2D expm(double t);

}
