package treelikelihood;

import java.util.Vector;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class TaylorMatrixExp implements MatrixExponentiator {

	Vector<DoubleMatrix2D> cachedQnDivFactorialN = new Vector<DoubleMatrix2D>();
	Vector<Double> cachedScale = new Vector<Double>();
	DoubleMatrix2D Q;
	// Precision Parameters...
	int nTaylor = 500;	
	
	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;
	
	DoubleFactory2D F = DoubleFactory2D.dense; // TODO: is parallelizable? 	
	

	public TaylorMatrixExp(DoubleMatrix2D Q_) {
		Q = Q_;
		DoubleMatrix2D next = F.identity(Q.rows()); // Q^0/0! = I
		double scale = getScale(next);		
			
		for (int i=0;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*Q/N
			cachedQnDivFactorialN.add(i,next);
			cachedScale.add(i,scale);
			next = next.zMult(Q,null,scale/(i+1.0),0,false,false);
			scale = getScale(next);
			next.assign(cern.jet.math.tdouble.DoubleFunctions.mult(1.0/scale));
		}
	}

	private double getScale(DoubleMatrix2D Q) {
		double scale = Math.sqrt(Math.abs(Q.getMinLocation()[0])*Math.abs(Q.getMaxLocation()[0])); 
		if (scale!=0) 
			return scale;
		else
			return 1.0;
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		DoubleMatrix2D result = cachedQnDivFactorialN.get(0).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(cachedScale.get(0)*Math.pow(t, 0)));

		for (int i=1;i<nTaylor;i++) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(cachedScale.get(i)*Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);				
		}
		return result; 	

	}
	

}
