package seasmig.treelikelihood.models;
import java.util.HashMap;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.EigenDecomposionExp;
import seasmig.treelikelihood.matrixexp.HKY85MatrixExp;
import seasmig.treelikelihood.matrixexp.JC69MatrixExp;
//import seasmig.treelikelihood.matrixexp.EigenDecomposionExp;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;

@SuppressWarnings("serial")
public class ConstantTransitionBaseModel implements TransitionModel {

	// Precision Parameter
	static final double timePrecision = 1.0E-5;
	static final double veryLongTime = 1000; // TODO: move to config...

	// Rate Matrix  
	double[][] Q;
	private int num_locations = 0;		

	// cache
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();
	final static int maxCachedTransitionMatrices=16000;

	// Matrix Exponentiation
	MatrixExponentiator matrixExponentiator;



	private DoubleMatrix1D basefreq;	

	protected ConstantTransitionBaseModel() {};

	// Constructor	
	public ConstantTransitionBaseModel(double[][] Q_) {	
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
			// TODO: check if this helps
			matrixExponentiator=new EigenDecomposionExp(Q);
			if (!matrixExponentiator.checkMethod()) {
				matrixExponentiator=new Matlab7MatrixExp(Q);
			}
		}
		// TODO: Maybe use other method to get s.s. freq 
		basefreq=matrixExponentiator.expm(veryLongTime).viewRow(1);
		// TODO: Figure out if calculating base freq helps or interferes. 
		//		basefreq=new double[num_locations];
		//		for (int i=0;i<num_locations;i++) {
		//			basefreq[i]=1.0/num_locations;
		//		}
	}

	// Constructor	
	public ConstantTransitionBaseModel(double mu) {	// JC69
		matrixExponentiator=new JC69MatrixExp(mu);
	}
	
	// Constructor	
	public ConstantTransitionBaseModel(double mu, double kappa, double piC, double piA, double piG) {	// HKY85
		matrixExponentiator=new HKY85MatrixExp(mu, kappa, piC, piA, piG);
	}


	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {

		double dt = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 
		DoubleMatrix2D cached = cachedTransitionMatrices.get(dt);

		if (cached!=null)  
			return Math.log(cached.get(from_location,to_location));		
		else 		
			return Math.log(transitionMatrix(from_time, to_time).get(from_location,to_location));		
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {		
		double dt = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 

		DoubleMatrix2D cached = cachedTransitionMatrices.get(dt);
		if (cached!=null) {
			return cached;
		}
		else {
			DoubleMatrix2D result = matrixExponentiator.expm(dt);
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
	public DoubleMatrix1D probability(int from_location, double from_time,	double to_time) {
		double dt = Math.max(timePrecision,cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 
		DoubleMatrix2D cached = cachedTransitionMatrices.get(dt);

		if (cached!=null)  
			return cached.viewRow(from_location);		
		else 		
			return transitionMatrix(from_time, to_time).viewRow(from_location);	
	}

	@Override
	public DoubleMatrix1D rootfreq(double when) {
		return basefreq;
	}

	@Override
	public Transition nextEvent(double from_time, int from) {
		double lambda = -Q[from][from];

		double time = cern.jet.random.Exponential.staticNextDouble(lambda); 
		// mean of this exponential is 1/lambda, higher the rate, the shorter the time --> this is the correct direction.
		
		int first = (from+1)%num_locations;
		double p=-Q[from][first]/Q[from][from];
		double rnd = cern.jet.random.Uniform.staticNextDouble();
		int loc = first;
		while (rnd>p) {
			loc=(loc+1)%num_locations;
			p = p-Q[from][loc]/Q[from][from];
		}
		// TODO: check this

		return new Transition(time+from_time, loc);
	}
	

}
