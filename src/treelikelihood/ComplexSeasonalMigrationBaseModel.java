package treelikelihood;
//TODO: This code has major flaws... :)

import java.util.HashMap;
import java.util.Vector;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdcomplex.DComplexFactory2D;
import cern.colt.matrix.tdcomplex.DComplexMatrix2D;
import cern.colt.matrix.tdcomplex.impl.DenseDComplexMatrix2D;
import org.javatuples.Pair;

public class ComplexSeasonalMigrationBaseModel implements MigrationBaseModel {
	
	// Precision Parameters...
	static final double precisionGoal = 1E-15; // Individual probability calculation precision goal
	static final int maxN = 500; // Maximum number of taylor series values

	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;

	// Rate Matrix  
	DenseDComplexMatrix2D Q;
	private int num_locations = 0;	


	// Caching 
	DoubleFactory2D F = DoubleFactory2D.dense;	
	Vector<DComplexMatrix2D> cachedMatrixPower = new Vector<DComplexMatrix2D>();	
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, DoubleMatrix2D>();
	double[] logFactorial = new double[maxN+1];

	// Constructor	
	public ComplexSeasonalMigrationBaseModel(double[][]c1) {	
		Q = new DenseDComplexMatrix2D(c1);		
		num_locations =Q.rows();
		cachedMatrixPower.add(0, DComplexFactory2D.dense.identity(Q.rows())); // Q^0 = I
		cachedMatrixPower.add(1,Q);  // Q^1 = Q
		for (int i=0;i<logFactorial.length;i++) {
			logFactorial[i]=cern.jet.math.Arithmetic.logFactorial(i);
		}				
		
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

		if (result<0) // TODO: deal with negative probability....
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

		DoubleMatrix2D cached = cachedTransitionMatrices.get(new Pair<Double,Double>(to_time,from_time));
		if (cached!=null) {
			return cached;
		}
		else {
			// Compute Taylor expansion to calculate P(t)=Exp(Qt) 
			int n = 0;
			DComplexMatrix2D result = matrixPowerQ(n).copy(); 
			double t = (to_time - from_time);
			double logt = Math.log(t);		
			double precision=0;

			do {
				n=n+1;		
				DComplexMatrix2D Qn = matrixPowerQ(n);
				precision=Math.max(Math.abs(Qn.getRealPart().getMinLocation()[0]),Math.abs(Qn.getRealPart().getMaxLocation()[0]))*Math.exp(n*logt-logFactorial[n]); 	
				double[] alpha = {Math.exp(n*logt-logFactorial[n]),0.0};				
				result.assign(Qn,cern.jet.math.tdcomplex.DComplexPlusMultSecond.plusMult(alpha));
			} while (precision>precisionGoal && (n<maxN));	

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			// TODO: for yearly cycles, change cache 
			cachedTransitionMatrices.put(new Pair<Double,Double>(to_time,from_time), result.getRealPart());
			return result.getRealPart();
		}
	}
	
	@Override
	public String print() {		
		return Q.toString();
	}
	
	@Override
	public int getNumLocations() {
		return num_locations;
	}
	
	private DComplexMatrix2D matrixPowerQ(int n) {
		for (int i=cachedMatrixPower.size();i<=n;i++) {
			cachedMatrixPower.add(i,cachedMatrixPower.get(i-1).zMult(Q,null));						
		}		
		return cachedMatrixPower.get(n);
	}

}