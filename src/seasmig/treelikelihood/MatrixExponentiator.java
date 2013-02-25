package seasmig.treelikelihood;

import java.io.Serializable;

public interface MatrixExponentiator extends Serializable {
	
	// Should return exp(Q*t)
	double[][] expm(double t);

}
