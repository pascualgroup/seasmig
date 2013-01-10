package treelikelihood;
import java.util.HashMap;
import java.util.Vector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.PlusMult;
import cern.jet.math.tdouble.DoublePlusMultSecond;

public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameters...
	int nTaylor = 1000;	
	static final double maxValue = Double.MAX_VALUE/1024;
	static final double minValue = Double.MIN_VALUE*1024;

	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;

	// Rate Matrix  
	DoubleMatrix2D Q;
	private int num_locations = 0;		

	// Caching 
	DoubleFactory2D F = DoubleFactory2D.dense;	
	Vector<DoubleMatrix2D> cachedQnDivFactorialN = new Vector<DoubleMatrix2D>();
	Vector<Double> cachedMatrixScale = new Vector<Double>();
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();

	// Constructor	
	public ConstantMigrationBaseModel(double[][] Q_) {	
		num_locations=Q_.length;

		Q = F.make(Q_);
// TODO:scale cache
		addToCache(0,F.identity(Q.rows()));	// Q^0/0! = I
		DoubleMatrix2D next = getFromCache(0);
		for (int i=1;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*scale*Q/N
			next = getFromCache(i-1).zMult(Q,null,getScale(i-1)/(i+1.0),0,false,false);						
			cachedQnDivFactorialN.add(i,next);
			// TODO: check for overflow or underflow	
			for (int j=0;j<next.rows();j++) {
				for (int k=0;k<next.rows();k++) {
					if (Math.abs(next.get(j, k))>maxValue) {
						nTaylor=i-1;
					}
					if (Math.abs(next.get(j, k))<minValue) {
						nTaylor=i-1;
					}
				}
			}
		}

	}

	private DoubleMatrix2D getFromCache(int i) {
		// TODO Auto-generated method stub
		return null;
	}

	private double getScale(int i) {
		// TODO Auto-generated method stub
		return cachedMatrixScale.get(i);
	}

	private void addToCache(int i, DoubleMatrix2D identity) {
		// TODO Auto-generated method stub
		
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
			double t = (to_time - from_time);
			DoubleMatrix2D result = taylorSeriesQ(0,t);

			for (int i=1;i<nTaylor;i++) {
				DoubleMatrix2D taylorn = taylorSeriesQ(i,t);
				result.assign(taylorn, cern.jet.math.tdouble.DoubleFunctions.plus);				
			} 	

			// TODO: this
			// fix negative and greater than 1 values
			for (int i=0;i<result.rows();i++) {
				for (int j=0;j<result.rows();j++) {
					if (result.get(i, j)<0) {
						result.set(i, j,minValue);
					}
					if (result.get(i, j)>1) {
						result.set(i, j,maxValue);
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

	private DoubleMatrix2D taylorSeriesQ(int n, double t) {
		// TODO: fix and check this...
		if (n<nTaylor) {
			return cachedQnDivFactorialN.get(n).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(Math.pow(t, n)*getScale(n)));
		}
		else {
			// TODO: deal with this...
			System.err.println("shouldn't be here debug!");
			return null;

		}

	}



}
