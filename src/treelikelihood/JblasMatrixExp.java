package treelikelihood;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;


public class JblasMatrixExp implements MatrixExponentiator {
	
	private DoubleMatrix Q;
	
	public JblasMatrixExp(double[][] Q_) {
		Q=new DoubleMatrix(Q_);
	}
	@Override
	public double[][] expm(double t) {
		DoubleMatrix A = Q.mul(t);
		
		return MatrixFunctions.expm(A).toArray2();
	}

}

	