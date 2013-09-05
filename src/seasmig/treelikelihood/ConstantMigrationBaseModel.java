package seasmig.treelikelihood;
import java.util.HashMap;

import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;

@SuppressWarnings("serial")
public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameter
	static final double timePrecision = 1.0E-5;
	static final double veryLongTime = 1000; // TODO: move to config...

	// Rate Matrix  
	double[][] Q;
	private int num_locations = 0;		

	// cache
	HashMap<Double, double[][]> cachedTransitionMatrices = new HashMap<Double, double[][]>();
	final static int maxCachedTransitionMatrices=16000;

	// Matrix Exponentiation
	MatrixExponentiator matrixExponentiator;

	private double[] basefreq;	

	protected ConstantMigrationBaseModel() {};
	
	// Constructor	
	public ConstantMigrationBaseModel(double[][] Q_) {	
		Q = Q_;
		num_locations=Q_.length;
		switch (num_locations) {
		case 2:
			matrixExponentiator=new AnalyticMatrixExp2(Q);
			if (!matrixExponentiator.checkMethod()) {
				matrixExponentiator=new Matlab7MatrixExp(Q);
			}
			break;
		case 3:
			matrixExponentiator=new AnalyticMatrixExp3(Q);
			if (!matrixExponentiator.checkMethod()) {
				matrixExponentiator=new Matlab7MatrixExp(Q);
			}
			break;			
		default:
			matrixExponentiator=new Matlab7MatrixExp(Q);	
		}
		// TODO: Maybe use other method to get s.s. freq 
		basefreq=matrixExponentiator.expm(veryLongTime)[1];
		// TODO: Figure out if calculating base freq helps or interferes. 
//		basefreq=new double[num_locations];
//		for (int i=0;i<num_locations;i++) {
//			basefreq[i]=1.0/num_locations;
//		}
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time, boolean reverseTime) {
		// TODO: note: can't mix reverseTime and no reverseTime in same instance because of current cache implementation
		double dt = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 
		double[][] cached = cachedTransitionMatrices.get(dt);

		if (cached!=null)  
			return Math.log(cached[from_location][to_location]);		
		else 		
			return Math.log(transitionMatrix(from_time, to_time,reverseTime)[from_location][to_location]);		
	}

	@Override
	public double[][] transitionMatrix(double from_time, double to_time, boolean reverseTime) {	
		// TODO: note: can't mix reverseTime and no reverseTime in same instance because of current cache implementation
		double dt = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 
		
		double[][] cached = cachedTransitionMatrices.get(dt);
		if (cached!=null) {
			return cached;
		}
		else {
			double[][] result = matrixExponentiator.expm(dt, reverseTime);
			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			cachedTransitionMatrices.put(dt, result);
			return result;
		}
	}

	@Override
	public String print() {		
		String returnValue = "{";
		for (int i=0; i<Q.length;i++) {
			if (i!=0) 
				returnValue+=" ";
			returnValue+="{";
			for (int j=0; j<Q.length;j++) {
				returnValue+=String.format("%6.4f",Q[i][j]);
				if (j!=Q.length-1) {
					returnValue+=",";
				}
			}			
			returnValue+="}";
			if (i!=Q.length-1) {
				returnValue+=",\n";
			}			
		}
		returnValue+="}\n";
		return returnValue;
	}
	
	@Override
	public String parse() {		
		String returnValue = "{";
		for (int i=0; i<Q.length;i++) {
			if (i!=0) 
				returnValue+=" ";
			returnValue+="{";
			for (int j=0; j<Q.length;j++) {
				returnValue+=String.format("%6.4f",Q[i][j]);
				if (j!=Q.length-1) {
					returnValue+=",";
				}
			}			
			returnValue+="}";
			if (i!=Q.length-1) {
				returnValue+=",";
			}			
		}
		returnValue+="}";
		return returnValue;
	}

	@Override
	public int getNumLocations() {
		return num_locations ;
	}

	@Override
	public String getModelName() {		
		return "Constant";
	}

	@Override
	public double[] rootfreq(double when) {
		return basefreq;
	}

}
