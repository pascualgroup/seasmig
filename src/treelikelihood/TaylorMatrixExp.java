package treelikelihood;

import java.util.Vector; // 

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class TaylorMatrixExp implements MatrixExponentiator {

	Vector<DoubleMatrix2D> cachedQnDivFactorialN = new Vector<DoubleMatrix2D>();
	Vector<Double> cachedScale = new Vector<Double>();
	DoubleMatrix2D Q;
	DoubleMatrix2D zeroMatrix;
	// Precision Parameters...
	int nTaylor = 500;	
	double minValue = 1E-150;
	double maxValue = 1E150;
	
	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;
	
	DoubleFactory2D F = DoubleFactory2D.dense; // TODO: is parallelizable? 	
	

	public TaylorMatrixExp(DoubleMatrix2D Q_) {
		Q = Q_;
		zeroMatrix = F.make(Q.rows(),Q.rows(),0);
		DoubleMatrix2D next = F.identity(Q.rows()); // Q^0/0! = I
		double scale = getScale(next);		
			
		for (int i=0;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*Q/N
			cachedQnDivFactorialN.add(i,next);
			cachedScale.add(i,scale);
			next = next.zMult(Q,null,scale/(i+1.0),0,false,false);
			scale = getScale(next);
			next.assign(cern.jet.math.tdouble.DoubleFunctions.div(scale));	
			if (scale<minValue) {
				nTaylor = i-1;
			}
		}
	}

	private double getScale(DoubleMatrix2D Q) {
		// TODO: redo this...
		double minScale=Double.MAX_VALUE;;
		double maxScale=0;
		for (int i=0; i<Q.rows(); i++) {
			for (int j=0; j<Q.rows(); j++) {
				if (Math.abs(Q.get(i, j))<minScale) {
					minScale = Math.abs(Q.get(i, j));
				}
				if (Math.abs(Q.get(i, j))>maxScale) {
					maxScale = Math.abs(Q.get(i, j));
				}
			}
		}
		double scale = Math.sqrt(minScale*maxScale);
		if (scale!=0) 
			return scale;
		else
			return 1.0;
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		
		if (t>0.025) {
			DoubleMatrix2D halfProb = expm(t/2.0);
			return halfProb.zMult(halfProb, null);
		}
		
		DoubleMatrix2D result = zeroMatrix.copy();

		for (int i=nTaylor-1;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(cachedScale.get(i)*Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);
		}
		for (int i=nTaylor-2;i>=0;i=i-2) {
			DoubleMatrix2D taylorn = cachedQnDivFactorialN.get(i).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(cachedScale.get(i)*Math.pow(t, i)));
			result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);
		}
		
//		if (Math.random()<0.0001) {
//			System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());	
//		}
		
		// TODO: remove debug, at least partially
		for (int i=0;i<result.rows();i++) {
			for (int j=0;j<result.rows();j++) {
				if (result.get(i, j)<0) {
					result.set(i, j, minValue);
					System.err.println("result.get(i, j)<0");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (result.get(i,j)>1) {
					result.set(i, j, 1-minValue);
					System.err.println("result.get(i, j)>1");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (result.get(i, j)==Double.NaN) {
					System.err.println("result.get(i, j)==Double.NaN");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (result.get(i, j)==Double.MAX_VALUE) {
					System.err.println("result.get(i, j)==Double.MAX_VALUE");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}
				if (result.get(i, j)==Double.MIN_VALUE) {
					System.err.println("result.get(i, j)==Double.MIN_VALUE");
					System.err.println("\nQ:"+Q.toString()+"\nt: "+t+"\nexpm:"+result.toString());
				}				
			}
		}
		//
		return result; 	
		
	}
	

}
