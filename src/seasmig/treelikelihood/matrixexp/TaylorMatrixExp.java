package seasmig.treelikelihood.matrixexp;

import java.util.Vector; // 

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.util.Util;

@SuppressWarnings("serial")
public class TaylorMatrixExp implements MatrixExponentiator {
	// TODO: this has many numerical issues...

	// Precision Parameters...	
	static final double precision = Util.minValue;
	static int nTaylor = 500; // Maximum number of Taylor series terms

	// Cache 
	Vector<DoubleMatrix2D> cachedQnDivFactorialN = new Vector<DoubleMatrix2D>();
	Vector<Double> cachedScale = new Vector<Double>();
	DoubleMatrix2D Q;
	DoubleMatrix2D zeroMatrix;
	boolean checkMethod = true;

	DoubleFactory2D F = DoubleFactory2D.dense; 

	protected TaylorMatrixExp() {};

	public TaylorMatrixExp(double[][] testMatrix) {		

		// TODO: this
		if (testMatrix.length>40) 
			checkMethod=false;
		else {
			Q = DoubleFactory2D.dense.make(testMatrix);

			zeroMatrix = F.make(Q.rows(),Q.rows(),0);
			DoubleMatrix2D next = F.identity(Q.rows()); // Q^0/0! = I
			double scale = getScale(next);		

			for (int i=0;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*Q/N
				cachedQnDivFactorialN.add(i,next);
				cachedScale.add(i,scale);
				next = next.zMult(Q,null,scale/(i+1.0),0,false,false);
				scale = getScale(next);
				next.assign(cern.jet.math.Functions.div(scale));	
				if (scale<precision) {
					nTaylor = i-1;
				}
			}
		}
	}

	private double getScale(DoubleMatrix2D Q) {	
		double norm1=0;
		for (int i=0; i<Q.columns(); i++) {
			double absRowSum=0;
			for (int j=0; j<Q.rows();j++) {
				absRowSum+=Q.get(j,i);
			}
			if (absRowSum>norm1)
				norm1=absRowSum;
		}
		return norm1;
	}

	@Override
	public double[][] expm(double t, boolean transpose) {
		// TODO: remove extra construction steps...
		// TODO: deal with different scales of t
		if (t>0.1) {
			DoubleMatrix2D halfProb = DoubleFactory2D.dense.make(expm(t/2.0));
			return halfProb.zMult(halfProb, null).toArray();
		}

		DoubleMatrix2D result = zeroMatrix.copy();

		for (int i=nTaylor-1;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.Functions.mult(cachedScale.get(i)*Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.Functions.plus);
		}
		for (int i=nTaylor-2;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.Functions.mult(cachedScale.get(i)*Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.Functions.plus);
		}

		if (transpose)
			return result.viewDice().toArray();
		else
			return result.toArray();

	}

	@Override
	public boolean checkMethod() {
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

	@Override
	public double[][] expm(double t) {
		return expm(t, false);
	}

}
