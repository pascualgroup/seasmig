package treelikelihood;
import java.util.HashMap;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class ConstantMigrationBaseModel implements MigrationBaseModel {

	// Rate Matrix  
	DoubleMatrix2D Q;
	private int num_locations = 0;		

	// cache
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();
	final static int maxCachedTransitionMatrices=6000;
	
	// Matrix Exponentiation
	MatrixExponentiator matrixExponentiator;
	DoubleFactory2D F = DoubleFactory2D.dense;	
	
	// Constructor	
	public ConstantMigrationBaseModel(double[][] Q_) {	
		Q = F.make(Q_);
		num_locations=Q_.length;
		matrixExponentiator=new Matlab7MatrixExp(Q);
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {

		if (to_location==MigrationBaseModel.UNKNOWN_LOCATION) 
			return 0;

		DoubleMatrix2D cached = cachedTransitionMatrices.get(to_time-from_time);
		
		if (cached!=null)  
			return Math.log(cached.get(from_location, to_location));		
		else 		
			return Math.log(transitionMatrix(from_time, to_time).get(from_location, to_location));		
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {		

		DoubleMatrix2D cached = cachedTransitionMatrices.get(to_time-from_time);
		if (cached!=null) {
			return cached;
		}
		else {
			DoubleMatrix2D result = matrixExponentiator.expm(to_time-from_time);
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





}
