package treelikelihood;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.tdouble.*;


public class SinusoidialSeasonalMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...
	
	private int num_states = 0;
	GeneralSeasonalMigrationBaseModel baseModel;
	
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
		for (int i=0;i<rates.length;i++) {
			for (int j=0;j<rates.length;j++) {
				if (i!=j) {
					seasonalRates[i][j]=new SeasonalRatePlusSineFunction(rates[i][j],amp[i][j],phase[i][j]);
				}
			}
		}
		baseModel = new GeneralSeasonalMigrationBaseModel(seasonalRates);
	}

	// Methods
	@Override
	public double logprobability(int from_state, int to_state, double from_time, double to_time) {		
		return Math.log(transitionMatrix(from_time, to_time).get(from_state, to_state));
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



}
