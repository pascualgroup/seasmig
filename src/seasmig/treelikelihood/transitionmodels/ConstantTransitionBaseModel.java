package seasmig.treelikelihood.transitionmodels;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.EigenDecomposionExp;
import seasmig.treelikelihood.matrixexp.HKY85MatrixExp;
import seasmig.treelikelihood.matrixexp.JC69MatrixExp;
import seasmig.treelikelihood.matrixexp.CachedMatrixExponentiator;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;

@SuppressWarnings("serial")
public class ConstantTransitionBaseModel implements TransitionModel {

	// Precision Parameter
	static final double timePrecision = 1.0E-5;
	static final double veryLongTime = 1000; // TODO: move to config...

	// Rate Matrix  
	public double[][] Q = null;
	private int dimension = 0;		

	// Matrix Exponentiation
	MatrixExponentiator matrixExponentiator;

	private DoubleMatrix1D basefreq;	

	protected ConstantTransitionBaseModel() {};

	// Constructor	
	public ConstantTransitionBaseModel(double[][] Q_) {	
		Q = Q_;
		dimension=Q_.length;
		switch (dimension) {
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
			matrixExponentiator=new CachedMatrixExponentiator(new EigenDecomposionExp(Q));
			if (!matrixExponentiator.checkMethod()) {
				matrixExponentiator=new CachedMatrixExponentiator(new Matlab7MatrixExp(Q));
			}
		}
		// TODO: Maybe use other method to get s.s. freq		
		DoubleMatrix2D stationaryDist = matrixExponentiator.expm(veryLongTime);
		double[] baseFreqdoubleform = new double[dimension];
		for (int i=0;i<dimension;i++) {
			baseFreqdoubleform[i]=stationaryDist.get(0, i);
		}
		basefreq=new DenseDoubleMatrix1D(baseFreqdoubleform);	

		// TODO: Figure out if calculating base freq helps or interferes. 
		//		basefreq=new double[num_locations];
		//		for (int i=0;i<num_locations;i++) {
		//			basefreq[i]=1.0/num_locations;
		//		}
	}

	// Constructor	
	public ConstantTransitionBaseModel(double mu) {	// JC69
		// TODO: check if cache helps
		matrixExponentiator=new JC69MatrixExp(mu);		
		double[][] Q = {
				{0,mu*0.25, mu*0.25, mu*0.25},
				{mu*0.25, 0, mu*0.25, mu*0.25},
				{mu*0.25, mu*0.25, 0, mu*0.25},
				{mu*0.25, mu*0.25, mu*0.25, 0}};		

		for (int i=0;i<4;i++) {
			double rowsum=0;
			for (int j=0;j<4;j++) {
				rowsum=rowsum+Q[i][j];
			}
			Q[i][i]=-rowsum;
		}
		dimension=4;
		basefreq=new DenseDoubleMatrix1D(new double[]{0.25,0.25,0.25,0.25});
		this.Q=Q;
	}


	// Constructor	
	public ConstantTransitionBaseModel(double mu, double kappa, double piC, double piA, double piG) {	// HKY85
		// TODO: check if cache helps
		matrixExponentiator=new CachedMatrixExponentiator(new HKY85MatrixExp(mu, kappa, piC, piA, piG));
		basefreq=new DenseDoubleMatrix1D(new double[]{1.0-piC-piA-piG,piC,piA,piG});
		double piT = 1.0 - piA - piG - piC;
		double[][] Q = {
				{0,mu*kappa*piC, mu*piA, mu*piG},
				{mu*kappa*piT, 0, mu*piA, mu*piG},
				{mu*piT, mu*piC, 0, mu*kappa*piG},
				{mu*piT, mu*piC, mu*kappa*piA, 0}};		
		dimension=4;
		for (int i=0;i<4;i++) {
			double rowsum=0;
			for (int j=0;j<4;j++) {
				rowsum=rowsum+Q[i][j];
			}
			Q[i][i]=-rowsum;
		}
		this.Q=Q;
	}


	public ConstantTransitionBaseModel(String readLine) {
		// TODO Auto-generated constructor stub
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {
		return Math.log(transitionMatrix(from_time, to_time).get(from_location,to_location));		
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {		
		//double dt = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time-from_time)); 
		double dt = to_time-from_time;
		return matrixExponentiator.expm(dt);		
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
		return dimension ;
	}

	@Override
	public String getModelName() {		
		return "Constant";
	}

	@Override
	public DoubleMatrix1D probability(int from_location, double from_time,	double to_time) {
		return transitionMatrix(from_time, to_time).viewRow(from_location);	
	}

	@Override
	public DoubleMatrix1D rootfreq(double when) {
		return basefreq;
	}

	@Override
	public Transition nextEvent(double from_time, int from) {
		double lambda = -Q[from][from];
		double time = cern.jet.random.Exponential.staticNextDouble(lambda); // TODO: ??????????
		// mean of this exponential is 1/lambda, higher the rate, the shorter the time interval --> nextDouble(lambda) is the correct direction.
		//System.err.println("lambda: "+String.format("%.3f", lambda)+" interval: "+String.format("%.3f", time)+" Q:"+Util.parse(Q));
		
		int first = (from+1)%dimension;
		double p=-Q[from][first]/Q[from][from];
		double rnd = cern.jet.random.Uniform.staticNextDouble();
		int loc = first;
		while (rnd>p) {
			loc=(loc+1)%dimension;
			p = p-Q[from][loc]/Q[from][from];
			if (loc==first) {
				loc=(loc-1)%dimension;
				break;
			}			
		}
		// TODO: check this
		if (loc>=dimension) {
			System.err.println("error in stochastic mapping nextEvent, dimension="+dimension+" loc="+loc+" from="+from);
			System.exit(-1);
		}
		if (loc==from) {
			System.err.println("error in stochastic mapping nextEvent, dimension="+dimension+" loc="+loc+" from="+from);
			System.exit(-1);
		}
		return new Transition(time+from_time, loc);
	}


}
