package seasmig.treelikelihood.matrixexp;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
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
	
	@Override
	public boolean checkMethod() {
		// TODO: this....
		return true;
	}

}

	