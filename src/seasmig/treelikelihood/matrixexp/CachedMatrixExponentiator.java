package seasmig.treelikelihood.matrixexp;

import java.util.HashMap;

import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class CachedMatrixExponentiator implements MatrixExponentiator {

	boolean methodOK;
	private MatrixExponentiator matrixExp;
	
	// cache
	// TODO: improve caching method
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();
	final static int maxCachedTransitionMatrices=16000;

	protected CachedMatrixExponentiator() {};	

	public CachedMatrixExponentiator(MatrixExponentiator matrixExp_) {			
		matrixExp = matrixExp_;
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		DoubleMatrix2D cached = cachedTransitionMatrices.get(t);
		if (cached!=null) {
			return cached;
		}
		else {
			DoubleMatrix2D result = matrixExp.expm(t);
			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				for (int i=0;i<cachedTransitionMatrices.size()/3;i++) {	
					cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
				}				
			}			
			cachedTransitionMatrices.put(t, result);
			return result;
		}		
	}

	@Override
	public boolean checkMethod() {		 
		return matrixExp.checkMethod();
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


