package treelikelihood;

import java.util.Vector; // 

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class TaylorMatrixExp implements MatrixExponentiator {

	// Precision Parameters...	
	static final double precision = 1E-15;
	
	// Cache 
	static final int maxCachedTransitionMatrices = 16000;		
	Vector<DoubleMatrix2D> cachedQnDivFactorialN = new Vector<DoubleMatrix2D>();
	DoubleMatrix2D Q;
	DoubleMatrix2D zeroMatrix;
	private int nTaylor = 100;
	
	DoubleFactory2D F = DoubleFactory2D.dense; 
	
	public TaylorMatrixExp(DoubleMatrix2D Q_) {
		Q = Q_;
		zeroMatrix = F.make(Q.rows(),Q.rows(),0);
		DoubleMatrix2D next = F.identity(Q.rows()); // Q^0/0! = I
		int nTaylor = Integer.MAX_VALUE;		
		for (int i=0;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*Q/N
			cachedQnDivFactorialN.add(i,next);
			next = next.zMult(Q,null,1.0/(i+1.0),0,false,false);
		}
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		
		// TODO: deal with different scales of t
		if (t>0.025) {
			DoubleMatrix2D halfProb = expm(t/2.0);
			return halfProb.zMult(halfProb, null);
		}
		
		DoubleMatrix2D result = zeroMatrix.copy();

		for (int i=nTaylor-1;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(/*cachedScale.get(i)**/Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);
		}
		for (int i=nTaylor-2;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(/*cachedScale.get(i)**/Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);
		}
		
		// TODO: remove debug, at least partially
//		if (Math.random()<0.0001) {
//			System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());	
//		}
				
		for (int i=0;i<result.rows();i++) {
			for (int j=0;j<result.rows();j++) {
				if (result.get(i, j)<0) {
					result.set(i, j, precision);
					System.err.println("result.get(i, j)<0");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (result.get(i,j)>1) {
					result.set(i, j, 1-precision);
					System.err.println("result.get(i, j)>1");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (Double.isNaN(result.get(i, j))) {
					System.err.println("result.get(i, j)==Double.NaN");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
			}
		}
		
		return result; 	
		
	}
	

}
