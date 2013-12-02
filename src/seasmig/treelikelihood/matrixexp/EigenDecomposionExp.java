package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.util.Util;

@SuppressWarnings("serial")
public class EigenDecomposionExp implements MatrixExponentiator {

	static final Algebra algebra = new Algebra(Util.minValue);
	static final DoubleFactory2D F = DoubleFactory2D.dense;
	static final DoubleFactory1D F1 = DoubleFactory1D.dense;
	
	private DoubleMatrix2D V;
	private DoubleMatrix2D invV;
	private DoubleMatrix1D diagD;
	private boolean checkMethod;


	// Q = V*D*invV
	// Exp(Q*t) = V*Exp(D*t)*invV
	
	protected EigenDecomposionExp() {};
	
	public EigenDecomposionExp(double[][] Q_) {
		DoubleMatrix2D Q = F.make(Q_);	
		EigenvalueDecomposition eigenDecomposition = new EigenvalueDecomposition(Q);
		V = eigenDecomposition.getV();
		invV = algebra.inverse(V);		
		diagD=F.diagonal(eigenDecomposition.getD());

		// Method works when eigenvalues are not repeating (need check why this is)
		checkMethod=true;
		for (int i=1;i<diagD.size();i++) {
			checkMethod=checkMethod && (diagD.get(i-1)!=diagD.get(i));
		}
		
	}

		
	@Override
	public double[][] expm(double tt) {
		/*
		 * DoubleMatrix1D 	zMult(DoubleMatrix1D y, DoubleMatrix1D z, double alpha, double beta, boolean transposeA)
          Linear algebraic matrix-vector multiplication; z = alpha * A * y + beta*z.
 			
 			DoubleMatrix2D 	zMult(DoubleMatrix2D B, DoubleMatrix2D C, double alpha, double beta, boolean transposeA, boolean transposeB)
          Linear algebraic matrix-matrix multiplication; C = alpha * A x B + beta*C.
		 */
		// TODO: remove for
		DoubleMatrix2D expDt = F.diagonal(diagD.copy().assign(cern.jet.math.Functions.mult(tt)).assign(cern.jet.math.Functions.exp));	
		return V.zMult(expDt.zMult(invV,null,1,0,false,false),null,1,0,false,false).toArray();
	}
	
	@Override
	public boolean checkMethod() {
		// TODO: this....	
		return checkMethod;
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
