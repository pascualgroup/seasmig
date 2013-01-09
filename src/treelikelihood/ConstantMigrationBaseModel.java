package treelikelihood;
import java.util.HashMap;
import java.util.Vector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.PlusMult;
import cern.jet.math.tdouble.DoublePlusMultSecond;

public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameters...
	public static final double precisionGoal = 1E-25; // Individual probability calculation precision goal
	static final int maxIter = 10000; // Maximum number of Taylor series values
	static final int minIter = 30; // Minimum number of Taylor series values
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

//		if (result==0) 
//			return Math.log(precisionGoal);
//		if (result==Double.NaN)
//			return Math.log(precisionGoal);
//		if (result==1)
//			return Math.log(1-precisionGoal);

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
			boolean probabilityOK=false;
			do {
				n=n+1;		
				DoubleMatrix2D Qn = matrixPowerQ(n);
				// TODO: Fix precision to handle both 0 and 1 side... 
				// TODO: this...
				DoubleMatrix2D diff = F.dense.make(Q.rows(),Q.rows(),0);
				double taylorCoeff = Math.exp(n*logt-cern.jet.math.tdouble.DoubleArithmetic.logFactorial(n));
				diff.assign(Qn,DoublePlusMultSecond.plusMult(taylorCoeff));
				result.assign(diff, cern.jet.math.tdouble.DoubleFunctions.plus);
				precision=Math.max(Math.abs(diff.getMinLocation()[0]),Math.abs(diff.getMaxLocation()[0]));
				probabilityOK = (result.getMaxLocation()[0]<=1.0) && (result.getMinLocation()[0]>=0); 
			} while ((precision>precisionGoal || !probabilityOK || n<minIter) && n<maxIter);	

			if (n==maxIter) {
				for (int i=0; i<Q.rows();i++) {
					for (int j=0; j<Q.rows();j++) {
						if (result.get(i,j)<=0) {
							result.set(i,j, precisionGoal);
						}
						if (result.get(i, j)>=1) {
							result.set(i, j, 1-precisionGoal);
						}
					}
				}
			}
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
