package treelikelihood;

public class MigrationPlusBrowninanBaseModel implements MigrationPlusContinousStateBaseModel {
	// TODO: go over this...
	// P_location(from_location,to_location)*P_state(from_location,to_location)(from_state,to_state)
	
	
	BrownianMotionPlusSeasonalityFunction logSeasonalStatesTransitionProbabilities[][] = null;
	private MigrationBaseModel migrationBaseModel;	
	
	public class BrownianMotionPlusSeasonalityFunction {

		private double alpha;
		private double amp_from;
		private double phase_from;
		private double amp_to;
		private double phase_to;
		

		public BrownianMotionPlusSeasonalityFunction(double alpha_, double amp_from_, double phase_from_, double amp_to_, double phase_to) {
			alpha=alpha_;
			phase_from=phase_from_;
			amp_from=amp_from_;
			amp_to=amp_to_;
			amp_from=amp_from_;
		}

		@Override
		public String toString() {
			return "Seasonality: "+amp_from+"*sin(2Pi*t+2Pi*"+phase_from+"),"+amp_to+"*sin(2Pi*t+2Pi*"+phase_to+"), Wiener: N(0,"+alpha+"*t)";			
		}
		
		public double apply(double from_time, double to_time, double from_state, double to_state) {		
			// TODO: go over this...
			// When removing the seasonality the adjusted state is assumed to be a Wiener process (Brownian Motion)
			// Wt-Ws ~ N(0, alpha*(t-s)) (for 0<s<t)
			double xs = from_state - amp_from*Math.sin(2*Math.PI*from_time+2*Math.PI*phase_from);
			double xt = to_state - amp_to*Math.sin(2*Math.PI*to_time+2*Math.PI*phase_to);
			double var = Math.abs(to_time-from_time)*alpha;			
			return -Math.log(2*Math.PI*var)/2.0+(xt-xs)*(xt-xs)/(2*var); 
		}

	}

	// Constructor	
	public MigrationPlusBrowninanBaseModel(MigrationBaseModel migrationBaseModel_, double[][] alphas, double[] amp, double[] phase) {
		// TODO: Lookup i and j from to....
		migrationBaseModel=migrationBaseModel_;
		logSeasonalStatesTransitionProbabilities = new BrownianMotionPlusSeasonalityFunction[alphas.length][alphas.length];		
		for (int i=0;i<alphas.length;i++) {
			for (int j=0;j<alphas.length;j++) {
				logSeasonalStatesTransitionProbabilities[i][j]=new BrownianMotionPlusSeasonalityFunction(alphas[i][j],amp[i],phase[i],amp[j],phase[j]);
			}
		}
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time, double from_state, double to_state) {
		return migrationBaseModel.logprobability(from_location, to_location, from_time, to_time)+logSeasonalStatesTransitionProbabilities[from_location][to_location].apply(from_time, to_time, from_state, to_state);
	}

	@Override
	public String print() {
		String returnValue = "Continuous Seasonal Migration Model:\n";
		returnValue=returnValue+migrationBaseModel.print();
		returnValue+="\n[";
		for (int i=0;i<logSeasonalStatesTransitionProbabilities.length;i++) {
			if (i!=0) returnValue+=" ";
			returnValue+="[";
			for (int j=0;j<logSeasonalStatesTransitionProbabilities[i].length;j++) {
				returnValue=returnValue+String.format("%30s",logSeasonalStatesTransitionProbabilities[i][j].toString());				
				if (j!=logSeasonalStatesTransitionProbabilities[i].length-1) returnValue+=",\t";
			}
			returnValue+="]";
			if (i!=logSeasonalStatesTransitionProbabilities.length-1) returnValue+="\n";
		}
		returnValue+="]\n";
		return returnValue;	
	}

}
