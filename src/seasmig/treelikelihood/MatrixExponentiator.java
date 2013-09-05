package seasmig.treelikelihood;

import java.io.Serializable;

public interface MatrixExponentiator extends Serializable {
	
	double[][] expm(double t); 	// return exp(Q*t)
	double[][] expm(double t, boolean transpose); // return exp(Q*t) or Transpose(exp(Q*t))
	boolean checkMethod(); // check if method is applicable to a specific Q matrix
	String getMethodName();
}
