package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class JC69MatrixExp implements MatrixExponentiator {

	boolean methodOK;
	private double mu;

	protected JC69MatrixExp() {};	

	public JC69MatrixExp(double mu) {
		methodOK = mu>0;
		this.mu = mu;
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		DoubleMatrix2D returnValue = DoubleFactory2D.dense.make(4,4);

		double exptmu = cern.jet.math.Functions.exp.apply(-1*t*mu);
		double offdiag = 1.0/4.0 - 1.0/4.0*exptmu;
		double diag = 1.0/4.0 + 3.0/4.0*exptmu;
		 
		returnValue.setQuick(0, 1, offdiag);
		returnValue.setQuick(0, 2, offdiag);
		returnValue.setQuick(0, 3, offdiag);
		returnValue.setQuick(1, 0, offdiag);
		returnValue.setQuick(1, 2, offdiag);
		returnValue.setQuick(1, 3, offdiag);
		returnValue.setQuick(2, 0, offdiag);
		returnValue.setQuick(2, 1, offdiag);
		returnValue.setQuick(2, 3, offdiag);
		returnValue.setQuick(3, 0, offdiag);
		returnValue.setQuick(3, 1, offdiag);
		returnValue.setQuick(3, 2, offdiag);
		returnValue.setQuick(0, 0, diag);
		returnValue.setQuick(1, 1, diag);
		returnValue.setQuick(2, 2, diag);
		returnValue.setQuick(3, 3, diag);
						
		return returnValue;
	}


	@Override
	public boolean checkMethod() {		
		return methodOK;
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
