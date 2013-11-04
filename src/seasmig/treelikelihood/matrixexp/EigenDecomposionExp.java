package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.util.Util;

@SuppressWarnings("serial")
public class EigenDecomposionExp implements MatrixExponentiator {

	static final Algebra algebra = new Algebra(Util.minValue);
	static final DoubleFactory2D F = DoubleFactory2D.dense;
	
	private DoubleMatrix2D V;
	private DoubleMatrix2D invV;
	private DoubleMatrix1D D;

	// Q = V*D*invV
	// Exp(Q*t) = V*Exp(D*t)*invV
	
	protected EigenDecomposionExp() {};
	
	public EigenDecomposionExp(double[][] Q_) {
		DoubleMatrix2D Q = F.make(Q_);	
		EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(Q);
		V = eigenDecomposition.getV();
		invV = algebra.inverse(V);
		D = eigenDecomposition.getRealEigenvalues();
	}
	
	@Override
	public double[][] expm(double tt) {
		/*
		 * DoubleMatrix1D 	zMult(DoubleMatrix1D y, DoubleMatrix1D z, double alpha, double beta, boolean transposeA)
          Linear algebraic matrix-vector multiplication; z = alpha * A * y + beta*z.
 			
 			DoubleMatrix2D 	zMult(DoubleMatrix2D B, DoubleMatrix2D C, double alpha, double beta, boolean transposeA, boolean transposeB)
          Linear algebraic matrix-matrix multiplication; C = alpha * A x B + beta*C.
		 */
		DoubleMatrix1D Dt = D.copy().assign(cern.jet.math.Functions.mult(tt)).assign(cern.jet.math.Functions.exp);
		if (D.size()==3) {
			System.err.println(Util.print(F.diagonal(D).toArray()));
			System.err.println(Util.print(F.diagonal(Dt).toArray()));
		}
		DoubleMatrix2D matDt = F.diagonal(Dt); 
		return V.zMult(matDt.zMult(invV,null,1,0,false,false),null,1,0,false,false).toArray();
	}
	
	@Override
	public boolean checkMethod() {
		// TODO: this....
		// return Q == U*D*invU
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
