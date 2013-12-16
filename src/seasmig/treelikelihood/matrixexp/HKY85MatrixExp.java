package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class HKY85MatrixExp implements MatrixExponentiator {

	boolean methodOK;
	private double k;
	private double piC;
	private double piA;
	private double piG;
	private double piT;
	private double mu;

	protected HKY85MatrixExp() {};	

	public HKY85MatrixExp(double mu, double kappa, double piC, double piA, double piG) {			
		// kappa ratio of transitions to transversions
		// mu mutaiton rate
		// pi_j background or equilibrium frequency of base j

		this.mu = mu;
		this.k = kappa;
		this.piC = piC;
		this.piA = piA;
		this.piG = piG;
		this.piT = 1 - piC - piA - piG; 	
		
		methodOK = (mu>0) && (kappa>0) && (piC>=0) && (piC<=1.0) && (piA>=0) && (piA<=1.0) && (piG>=0) && (piG<=1.0) && ((piA+piC+piG)<1);

	}

	@Override
	public DoubleMatrix2D expm(double t) {
		// see HKY.nb
		// 

		// TCAG is position 0,1,2,3
		DoubleMatrix2D returnValue = DoubleFactory2D.dense.make(4,4);	

		double exp_mu_t = cern.jet.math.Functions.exp.apply(mu*t);
		
		returnValue.setQuick(0, 0, (piC/cern.jet.math.Functions.exp.apply((-1 + k)*mu*(piC + piT)*t) - piT*(-1 + piC + piT) + exp_mu_t*piT*(piC + piT))/(exp_mu_t*(piC + piT)));
		returnValue.setQuick(0, 1, (piC*(1 - cern.jet.math.Functions.exp.apply(-((-1 + k)*mu*(piC + piT)*t)) - piC - piT + exp_mu_t*(piC + piT)))/(exp_mu_t*(piC+piT)));
		returnValue.setQuick(0, 2, -(((-1 + exp_mu_t)*(-1 + piC + piG + piT))/exp_mu_t));
		returnValue.setQuick(0, 3, piG - piG/exp_mu_t);
		returnValue.setQuick(1, 0, (piT*(1 - cern.jet.math.Functions.exp.apply(-((-1 + k)*mu*(piC + piT)*t)) - piC - piT + exp_mu_t*(piC + piT)))/(exp_mu_t*(piC + piT)));
		returnValue.setQuick(1, 1, (piT/cern.jet.math.Functions.exp.apply((-1 + k)*mu*(piC + piT)*t) - piC*(-1 + piC + piT) + exp_mu_t*piC*(piC + piT))/(exp_mu_t*(piC + piT)));
		returnValue.setQuick(1, 2, returnValue.getQuick(0,2));
		returnValue.setQuick(1, 3, returnValue.getQuick(0,3));
		returnValue.setQuick(2, 0, piT - piT/exp_mu_t);
		returnValue.setQuick(2, 1, piC - piC/exp_mu_t);
		returnValue.setQuick(2, 2, (cern.jet.math.Functions.exp.apply((-1 + k)*mu*(-1 + piC + piT)*t)*piG + exp_mu_t*(-1 + piC + piT)*(-1 + piC + piG + piT) - (piC + piT)*(-1 + piC + piG + piT))/
				 (exp_mu_t*(piA + piG)));
		returnValue.setQuick(2, 3, (piG*(-cern.jet.math.Functions.exp.apply((-1 + k)*mu*(-1 + piC + piT)*t) + piC + piT - exp_mu_t*(-1 + piC + piT)))/(exp_mu_t*(piA + piG)));
		returnValue.setQuick(3, 0, returnValue.getQuick(2,0));
		returnValue.setQuick(3, 1, returnValue.getQuick(2,1));
		returnValue.setQuick(3, 2, ((-1 + piC + piG + piT)*(cern.jet.math.Functions.exp.apply((-1 + k)*mu*(-1 + piC + piT)*t) - piC - piT + exp_mu_t*(-1 + piC + piT)))/(exp_mu_t*(piA + piG)));
		returnValue.setQuick(3, 3, (-(exp_mu_t*piG*(-1 + piC + piT)) + piG*(piC + piT) - cern.jet.math.Functions.exp.apply((-1 + k)*mu*(-1 + piC + piT)*t)*(-1 + piC + piG + piT))/(exp_mu_t*(piA + piG)));

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


