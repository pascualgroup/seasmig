package seasmig.treelikelihood.models;

import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.MigrationPlusContinousStateBaseModel;
import seasmig.util.Util;

public class MigrationPlusSeasonalBrowninanBaseModel implements MigrationPlusContinousStateBaseModel {
	// TODO: go over this...
	// P_location(from_location,to_location)*P_state(from_location,to_location)(from_state,to_state)
	
	
	BrownianMotionPlusSeasonalityFunction logSeasonalStatesTransitionProbabilities[] = null;
	private MigrationBaseModel migrationBaseModel;	
	
	protected MigrationPlusSeasonalBrowninanBaseModel() {};
	
	public class BrownianMotionPlusSeasonalityFunction {

		private double alpha; 				
		private double ampAnnual;
		private double phaseAnnual;
		private double ampBiannual;
		private double phaseBiannual;
		

		public BrownianMotionPlusSeasonalityFunction(double alpha_, double ampAnnual_, double phaseAnnual_, double ampBiannual_, double phaseBiannual_) {
			alpha=alpha_;
			ampAnnual=phaseAnnual_;
			phaseAnnual=phaseAnnual_;
			ampBiannual=ampBiannual_;
			phaseBiannual=phaseBiannual_;
		}

		@Override
		public String toString() {
			return "Seasonality: "+ampAnnual+"*sin(2Pi*t+2Pi*"+phaseAnnual+")+"+ampBiannual+"*sin(4Pi*t+2Pi*"+phaseBiannual+"), Wiener: N(0,"+alpha+"*t)";			
		}
		
		public double apply(double from_time, double to_time, double from_state, double to_state) {		
			// TODO: go over this...
			// When removing the seasonality the adjusted state is assumed to be a Wiener process (Brownian Motion)
			// Wt-Ws ~ N(0, alpha*(t-s)) (for 0<s<t)
			// When locations are different we assume that the baseline migration model estimates the probability regardless of the actual state values 			
			double xs = from_state - ampAnnual*Math.sin(2*Math.PI*from_time+2*Math.PI*phaseAnnual)+ampBiannual*Math.sin(4*Math.PI*from_time+2*Math.PI*phaseBiannual);
			double xt = to_state -ampAnnual*Math.sin(2*Math.PI*to_time+2*Math.PI*phaseAnnual)+ampBiannual*Math.sin(4*Math.PI*to_time+2*Math.PI*phaseBiannual);
			double var = Math.abs(to_time-from_time+Util.minValue)*alpha;			
			return -Math.log(2*Math.PI*var)/2.0+(xt-xs)*(xt-xs)/(2*var); 
		}

	}

	// Constructor	
	public MigrationPlusSeasonalBrowninanBaseModel(MigrationBaseModel migrationBaseModel_, double[] alphas, double[] ampAnnual, double[] phaseAnnual, double[] ampBiannual, double[] phaseBiannual) {
		migrationBaseModel=migrationBaseModel_;
		logSeasonalStatesTransitionProbabilities = new BrownianMotionPlusSeasonalityFunction[alphas.length];		
		for (int i=0;i<alphas.length;i++) {
			logSeasonalStatesTransitionProbabilities[i]=new BrownianMotionPlusSeasonalityFunction(alphas[i],ampAnnual[i],phaseAnnual[i],ampBiannual[i],phaseBiannual[i]);
		}
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time, double from_state, double to_state) {
		double returnValue = migrationBaseModel.logprobability(from_location, to_location, from_time, to_time);
		if (from_location==to_location) {
			returnValue+=logSeasonalStatesTransitionProbabilities[from_location].apply(from_time, to_time, from_state, to_state); 
		}
		return returnValue; 
	}

	@Override
	public String print() {
		String returnValue = "Continuous Seasonal Migration Model:\n";
		returnValue=returnValue+migrationBaseModel.print();
		returnValue+="\n[";
		for (int i=0;i<logSeasonalStatesTransitionProbabilities.length;i++) {
			returnValue=returnValue+String.format("%s\n",logSeasonalStatesTransitionProbabilities[i].toString());				
			if (i!=logSeasonalStatesTransitionProbabilities.length-1) returnValue+=",\t";
		}			
		returnValue+="]\n";
		return returnValue;	
	}

}
