package seasmig.treelikelihood.models;

import java.util.HashMap;

import org.javatuples.Pair;

import seasmig.treelikelihood.MigrationBaseModel;
import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

@SuppressWarnings("serial")
public class GeneralSeasonalMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...
	// TODO: Check zMult order... sould be ok for Q where rows sum to 1...
	
	// Precision Parameter
	static final double infinitesimalTime = 1E-5;
	
	// Cache Parameters 
	static final int maxCachedTransitionMatrices = 1600;

	// Precision Parameters
	int nYearParts;

	// Origin Model
	DoubleFunction[][] seasonalRates;	
	DoubleFunction[] rootFreq;
	
	// Constant Migration Models
	MigrationBaseModel constantModels[];

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, DoubleMatrix2D>();

	private int num_locations = 0;
	private double dt = 1.0/(double)nYearParts; 

	protected GeneralSeasonalMigrationBaseModel() {};
	
	// Constructor	
	public GeneralSeasonalMigrationBaseModel(DoubleFunction[][] seasonalRates_, DoubleFunction[] rootFreq_, int nYearParts_) {	
		// TODO: Check this...
		// diagonal rates functions are calculated through row sums and are ignored...
		num_locations=seasonalRates_.length;	
		seasonalRates=seasonalRates_;
		nYearParts = nYearParts_;
		dt = 1.0/(double)nYearParts;
		constantModels = new ConstantMigrationBaseModel[nYearParts];
		double t=dt/2.0;
		for (int i=0;i<nYearParts;i++) {
			double[][] migrationMatrix = new double[num_locations][num_locations];
			for (int j=0; j<num_locations; j++) {
				double row_sum = 0;
				for (int k=0; k<num_locations; k++) {
					if (j!=k) {
						migrationMatrix[j][k]=seasonalRates[j][k].apply(t);
						row_sum+=migrationMatrix[j][k];
					}
				}
				migrationMatrix[j][j]=-row_sum;
			}
			constantModels[i]=new ConstantMigrationBaseModel(migrationMatrix);
			t+=dt;
		}		
		rootFreq = rootFreq_;
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {		
		return Math.log(transitionMatrix(from_time, to_time).get(from_location,to_location));
	}
	
	// Methods
	@Override
	public DoubleMatrix1D probability(int from_state,  double from_time, double to_time) {		
		return transitionMatrix(from_time, to_time).viewRow(from_state);
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {
		// TODO: organize this...
		double from_time_reminder = from_time % 1.0;
		double from_time_div = from_time - from_time_reminder;		
		double to_time_reminder = to_time - from_time_div;
		DoubleMatrix2D cached = cachedTransitionMatrices.get(new Pair<Double,Double>(from_time_reminder,to_time_reminder));
		if (cached!=null) {
			return cached;
		}
		else {			
			// first step: 
			double step_start_time = from_time_reminder;
			double step_end_time = Math.min(to_time_reminder, Math.floor(step_start_time/dt)*dt+dt);
			DoubleMatrix2D result = F.identity(num_locations);	 
			
			while (step_start_time<to_time_reminder) {
				int yearPartIndex = (int) Math.floor(step_start_time%1.0/dt);
				// TODO: replace with other matrix mult
				result = result.zMult(constantModels[yearPartIndex].transitionMatrix(step_start_time, step_end_time),null);	
				step_start_time = step_end_time;
				step_end_time = Math.min(to_time_reminder, Math.floor((step_start_time+infinitesimalTime)/dt)*dt+dt);
			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time_reminder, to_time_reminder),result);

			// TODO: replace with no conversion step
			return result;
		}
	}

	@Override
	public String print() {
		String returnValue = "General Seasonal Migration Model:\n";
		returnValue+="[";
		for (int i=0;i<seasonalRates.length;i++) {
			if (i!=0) returnValue+=" ";
			returnValue+="[";
			for (int j=0;j<seasonalRates[i].length;j++) {
				if (i!=j) 
					returnValue=returnValue+String.format("%40s",seasonalRates[i][j].toString());
				else 
					returnValue+=String.format("%40s", "NA");
				if (j!=seasonalRates[i].length-1) returnValue+=",";
			}
			returnValue+="]";
			if (i!=seasonalRates.length-1) returnValue+="\n";
		}
		returnValue+="]\n";
		return returnValue;	
	}


	@Override
	public int getNumLocations() {
		return num_locations ;
	}

	@Override
	public String parse() {
		// TODO Auto-generated method stub
		return print();
	}

	@Override
	public String getModelName() {		
		return "General Seasonal";
	}

	@Override
	public DoubleMatrix1D rootfreq(double when) {
		double[] returnValue = new double[num_locations];
		for (int i=0;i<num_locations;i++) {
			returnValue[i]=rootFreq[i].apply(when);
		}
		return DoubleFactory1D.dense.make(returnValue);
	}

	@Override
	public Transition nextEvent(double time, int from) {
		// TODO Auto-generated method stub
		return null;
	}


}
