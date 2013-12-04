package seasmig.treelikelihood;

import java.io.Serializable;

import cern.colt.matrix.DoubleMatrix2D;

public interface MatrixExponentiator extends Serializable {
	
	DoubleMatrix2D expm(double t); 	// return exp(Q*t)
	boolean checkMethod(); // check if method is applicable to a specific Q matrix
	String getMethodName();

}
