package seasmig.treelikelihood;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;


@SuppressWarnings("serial")
public class SinusoidialSeasonalMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...

	private int num_states = 0;
	GeneralSeasonalMigrationBaseModel baseModel;

	protected SinusoidialSeasonalMigrationBaseModel() {};
	
	public class SeasonalRatePlusSineFunction implements DoubleFunction {

		private double rate;
		private double amp;
		private double phase;


		public SeasonalRatePlusSineFunction(double rate_, double amp_, double phase_) {
			rate=rate_;
			phase=phase_;
			amp=amp_;
		}

		@Override
		public double apply(double t) {
			return rate*(1+amp*Math.sin(2*Math.PI*t+2*Math.PI*phase));
		}

		@Override
		public String toString() {
			return String.format("%.4f",rate)+"(1+"+String.format("%.4f",amp)+"*sin(2Pi*t+2Pi*"+String.format("%.4f",phase)+"))";			
		}

	}


	// Constructor	
	public SinusoidialSeasonalMigrationBaseModel(double[][] rates, double[][] amp, double[][] phase){
		num_states=rates.length;
		DoubleFunction seasonalRates[][] = new DoubleFunction[rates.length][rates.length];
		// TODO: calculate rootFreq ?
		DoubleFunction[] rootFreq = new DoubleFunction[rates.length];
		
		for (int i=0;i<rates.length;i++) {
			rootFreq[i] = cern.jet.math.Functions.constant(1.0/rates.length);
			for (int j=0;j<rates.length;j++) {
				if (i!=j) {
					seasonalRates[i][j]=new SeasonalRatePlusSineFunction(rates[i][j],amp[i][j],phase[i][j]);
				}
			}
		}
					
		baseModel = new GeneralSeasonalMigrationBaseModel(seasonalRates, rootFreq,24);
	}

	// Methods
	@Override
	public double logprobability(int from_state, int to_state, double from_time, double to_time) {		
		return Math.log(transitionMatrix(from_time, to_time).get(from_state,to_state));
	}

	// Methods
	@Override
	public DoubleMatrix1D probability(int from_state,  double from_time, double to_time) {		
		return transitionMatrix(from_time, to_time).viewRow(from_state);
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {
		return baseModel.transitionMatrix(from_time, to_time);
	}

	@Override
	public String print() {
		return baseModel.print();	
	}


	@Override
	public int getNumLocations() {
		return num_states ;
	}

	@Override
	public String parse() {
		// TODO Auto-generated method stub
		return print();
	}

	@Override
	public String getModelName() {		
		return "Sinusoidial Seasonal";
	}

	@Override
	public DoubleMatrix1D rootfreq(double when) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Event nextEvent(double time, int from) {
		// TODO Auto-generated method stub
		return null;
	}

}
