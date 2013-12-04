package seasmig.treelikelihood.matrixexp;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class JblasMatrixExp implements MatrixExponentiator {
	
	private DoubleMatrix Q;
	
	protected JblasMatrixExp() {};
	
	public JblasMatrixExp(double[][] Q_) {
		Q=new DoubleMatrix(Q_);
	}
	@Override
	public DoubleMatrix2D expm(double t) {
		return DoubleFactory2D.dense.make(MatrixFunctions.expm(Q.mul(t)).toArray2());
	}
	
	@Override
	public boolean checkMethod() {
		// TODO: this....
		return true;
	}

	@Override
	public String getMethodName() {		
		Class<?> enclosingClass = getClass().getEnclosingClass();
		if (enclosingClass != null) 
		    return enclosingClass.getName();
		else 
		    return getClass().getName();
		
	}
}

	