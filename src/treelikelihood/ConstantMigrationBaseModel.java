package treelikelihood;
import java.util.HashMap;
import java.util.Vector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoublePlusMultSecond;

public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameters...
	static final double precisionGoal = 1E-15; // Individual probability calculation precision goal
	static final int maxN = 10000; // Maximum number of Taylor series values

	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;

	// Rate Matrix  
	DoubleMatrix2D Q;
	private int num_locations = 0;		
	
	// Caching 
	DoubleFactory2D F = DoubleFactory2D.dense;	
	Vector<DoubleMatrix2D> cachedMatrixPower = new Vector<DoubleMatrix2D>();	
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();

	// Constructor	
	public ConstantMigrationBaseModel(double[][] Q_) {	
		Q = F.make(Q_);	
		num_locations=Q.rows();
		cachedMatrixPower.add(0,F.identity(Q.rows())); // Q^0 = I
		cachedMatrixPower.add(1,Q.copy());  // Q^1 = Q	
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {

		if (to_location==MigrationBaseModel.UNKNOWN_LOCATION) 
			return 0;

		double result=0;

		DoubleMatrix2D cached = cachedTransitionMatrices.get(to_time-from_time);
		
		if (cached!=null)  
			result=cached.get(from_location, to_location);		
		else 		
			result=transitionMatrix(from_time, to_time).get(from_location, to_location);		

		if (result<0) 
			result=precisionGoal;
		if (result>1) 
			result=1-precisionGoal;
		if (result==0) 
			return Double.MIN_VALUE;
		if (result==1)
			return Double.MAX_VALUE;

		return Math.log(result);
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {		

		DoubleMatrix2D cached = cachedTransitionMatrices.get(to_time-from_time);
		if (cached!=null) {
			return cached;
		}
		else {
			// Compute Taylor expansion to calculate P(t)=Exp(Qt) 
			// TODO: Replace with closed form expression for Exp(Qt) at least for specific n's 
			int n = 0;
			DoubleMatrix2D result = matrixPowerQ(n).copy(); 
			double t = (to_time - from_time);
			double logt = Math.log(t);		
			double precision=0;

			do {
				n=n+1;		
				DoubleMatrix2D Qn = matrixPowerQ(n);
				precision=Math.max(Math.abs(Qn.getMinLocation()[0]),Math.abs(Qn.getMaxLocation()[0]))*Math.exp(n*logt-cern.jet.math.tdouble.DoubleArithmetic.logFactorial(n)); 				
				result.assign(Qn,DoublePlusMultSecond.plusMult(Math.exp(n*logt-cern.jet.math.tdouble.DoubleArithmetic.logFactorial(n))));
			} while (precision>precisionGoal && (n<maxN));	

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			cachedTransitionMatrices.put(to_time-from_time, result);
			return result;
		}
	}
	
	@Override
	public String print() {		
		return Q.toString();
	}
	
	@Override
	public int getNumLocations() {
		return num_locations ;
	}

	private DoubleMatrix2D matrixPowerQ(int n) {
		for (int i=cachedMatrixPower.size();i<=n;i++) {
			cachedMatrixPower.add(i,cachedMatrixPower.get(i-1).zMult(Q,null));						
		}		
		return cachedMatrixPower.get(n);
	}
	
	

}
