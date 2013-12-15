package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class HKY85MatrixExp implements MatrixExponentiator {

	boolean methodOK;
	private double mu;
	private double kappa;
	private double piC;
	private double piA;
	private double piG;
	private double piT;
	private double beta;

	protected HKY85MatrixExp() {};	

	public HKY85MatrixExp(double mu, double kappa, double piC, double piA, double piG) {		
		this.mu = mu;
		this.kappa = kappa;
		this.piC = piC;
		this.piA = piA;
		this.piG = piG;
		this.piT = 1 - piC - piA - piG; 
		methodOK = (mu>0) && (kappa>0) && (piC>=0) && (piC<=1.0) && (piA>=0) && (piA<=1.0) && (piG>=0) && (piG<=1.0) && (Math.abs(1-piC-piA-piG-piT)<0.000001);
		this.beta = 1.0/( 2.0*(piA+piG)*(piC+piT)+2*kappa*((piA*piG)+(piC*piT)) );
		
	}

	@Override
	public DoubleMatrix2D expm(double t) {
		// TODO: check transpose
		// TCAG is position 0,1,2,3
		DoubleMatrix2D returnValue = DoubleFactory2D.dense.make(4,4);

		double expbetamut = cern.jet.math.Functions.exp.apply(-1*t*mu*beta);
	
		returnValue.setQuick(0, 1, ( piC*(piT+piC+(piG+piA)*expbetamut)-piC*expbetamut*cern.jet.math.Functions.exp.apply(1+(piT+piC)/(kappa-1)) )/(piT+piC)); // P_TC 
		returnValue.setQuick(0, 2, piA*(1-expbetamut)); // P_TA
		returnValue.setQuick(0, 3, piG*(1-expbetamut)); // P_TG
		returnValue.setQuick(1, 0, ( piT*(piC+piT+(piG+piA)*expbetamut)-piT*expbetamut*cern.jet.math.Functions.exp.apply(1+(piC+piT)/(kappa-1)) )/(piC+piT)); // P_CT
		returnValue.setQuick(1, 2, piA*(1-expbetamut)); // P_CA
		returnValue.setQuick(1, 3, piG*(1-expbetamut)); // P_CG
		returnValue.setQuick(2, 0, piT*(1-expbetamut)); // P_AT
		returnValue.setQuick(2, 1, piC*(1-expbetamut)); // P_AC
		returnValue.setQuick(2, 3, ( piG*(piA+piG+(piC+piT)*expbetamut)-piG*expbetamut*cern.jet.math.Functions.exp.apply(1+(piA+piG)/(kappa-1)) )/(piA+piG)); // P_AG
		returnValue.setQuick(3, 0, piT*(1-expbetamut)); // P_GT
		returnValue.setQuick(3, 1, piC*(1-expbetamut)); // P_GC
		
		returnValue.setQuick(3, 2, ( piA*(piG+piA+(piC+piT)*expbetamut)-piA*expbetamut*cern.jet.math.Functions.exp.apply(1+(piG+piA)/(kappa-1)) )/(piG+piA)); // P_GA;		
		returnValue.setQuick(0, 0, ( piT*(piT+piG+(piC+piA)*expbetamut)+piG*expbetamut*cern.jet.math.Functions.exp.apply(1+(piT+piG)*(kappa-1)) )/(piT+piG)); // P_TT
		returnValue.setQuick(1, 1, ( piC*(piC+piG+(piA+piT)*expbetamut)+piG*expbetamut*cern.jet.math.Functions.exp.apply(1+(piC+piG)*(kappa-1)) )/(piC+piG)); // P_CC
		returnValue.setQuick(2, 2, ( piA*(piA+piG+(piC+piT)*expbetamut)+piG*expbetamut*cern.jet.math.Functions.exp.apply(1+(piA+piG)*(kappa-1)) )/(piA+piG)); // P_AA
		returnValue.setQuick(3, 3, ( piG*(piG+piA+(piC+piT)*expbetamut)+piA*expbetamut*cern.jet.math.Functions.exp.apply(1+(piG+piA)*(kappa-1)) )/(piG+piA)); // P_GG
						
		return returnValue;
	}


	@Override
	public boolean checkMethod() {		
		// TODO: 
		//return methodOK;
		return false;
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
