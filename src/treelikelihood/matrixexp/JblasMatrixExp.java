package treelikelihood.matrixexp;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import treelikelihood.MatrixExponentiator;


public class JblasMatrixExp implements MatrixExponentiator {
	
	private DoubleMatrix Q;
	
	protected JblasMatrixExp() {};
	
	public JblasMatrixExp(double[][] Q_) {
		Q=new DoubleMatrix(Q_);
	}
	@Override
	public double[][] expm(double t) {
		return MatrixFunctions.expm(Q.mul(t)).toArray2();
	}

}

	