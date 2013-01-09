package treelikelihood;
import java.util.HashMap;
import java.util.Vector;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.PlusMult;
import cern.jet.math.tdouble.DoublePlusMultSecond;

public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameters...
	int nTaylor = 100;	
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
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();

	// Constructor	
	public ConstantMigrationBaseModel(double[][] Q_) {	
		num_locations=Q_.length;

		Q = F.make(Q_);
		cachedQnDivFactorialN.add(0,F.identity(Q.rows())); // Q^0/0! = I
		cachedQnDivFactorialN.add(1,Q.copy());  // Q^1/1! = Q
		DoubleMatrix2D next = Q.copy();
		for (int i=1;i<nTaylor;i++) { // Q^N/N! = Q^(N-1)/(N-1)!*Q/N
			next = next.zMult(Q,null,1.0/(i+1.0),0,false,false);						
			cachedQnDivFactorialN.add(i+1,next);
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
			return cachedQnDivFactorialN.get(n).copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(Math.pow(t, n)));
		}
		else {
			// TODO: deal with this...
			System.err.println("shouldn't be here debug!");
			return null;

		}

	}



}
