package seasmig.treelikelihood;

import java.util.HashMap;
import java.util.Vector;

import org.javatuples.Pair;
import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

@SuppressWarnings("serial")
public class GeneralSeasonalMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...
	// TODO: Fails degeneracy test....
	
	// Precision Parameter
	static final double timePrecision = 1E-5;
	
	// Cache Parameters 
	static final int maxCachedTransitionMatrices = 16000;

	// Precision Parameters
	static final int nYearParts = 36;
	// Figure out reason why changing nYearParts can generate an error...
	
	// Origin Model
	DoubleFunction[][] seasonalRates;
	
	// Constant Migration Models
	MigrationBaseModel constantModels[] = new ConstantMigrationBaseModel[nYearParts];

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	Vector<double[][]> cachedMatrixPower = new Vector<double[][]>();	
	HashMap<Pair<Double,Double>, double[][]> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, double[][]>();

	private int num_locations = 0;
	private double dt = 1.0/(double)nYearParts; 

	protected GeneralSeasonalMigrationBaseModel() {};
	
	// Constructor	
	public GeneralSeasonalMigrationBaseModel(DoubleFunction[][] seasonalRates_) {	
		// TODO: Check this...
		// diagonal rates are calculated on row sums and are ignored...
		num_locations=seasonalRates_.length;	
		seasonalRates=seasonalRates_;
		dt = 1.0/(double)nYearParts;
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
	}

	// Methods
	@Override
	public double logprobability(int from_location, int to_location, double from_time, double to_time) {		
		return Math.log(transitionMatrix(from_time, to_time)[from_location][to_location]);
	}
	
	// Methods
	@Override
	public double[] probability(int from_state,  double from_time, double to_time) {		
		return transitionMatrix(from_time, to_time)[from_state];
	}

	@Override
	public double[][] transitionMatrix(double from_time, double to_time) {
		// TODO: organize this...
		double from_time_reminder = from_time % 1.0;
		double from_time_div = from_time - from_time_reminder;
		double to_time_reminder = to_time - from_time_div;
		double from_time_reminder_round = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(from_time_reminder));
		double to_time_reminder_round = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time_reminder));
		double[][] cached = cachedTransitionMatrices.get(new Pair<Double,Double>(from_time_reminder_round,to_time_reminder_round));
		if (cached!=null) {
			return cached;
		}
		else {			
			// first step: 
			double step_start_time = from_time_reminder_round;
			double to_time_round = Math.max(timePrecision, cern.jet.math.Functions.round(timePrecision).apply(to_time));
			double step_end_time = Math.min(to_time_round, Math.max(0,Math.ceil(step_start_time/dt)*dt));
			DoubleMatrix2D result = F.identity(num_locations);	 
			
			while (step_end_time<to_time_round) {
				int yearPartIndex = (int) Math.floor(step_start_time%1.0/dt);
				// TODO: replace with other matrix mult
				result = result.zMult(DoubleFactory2D.dense.make(constantModels[yearPartIndex].transitionMatrix(step_start_time, step_end_time)),null);	
				step_start_time = step_end_time;
				step_end_time = Math.min(to_time_round, step_start_time+dt);
			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			double[][] returnValue=result.toArray();
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time_reminder_round, to_time_reminder_round),returnValue);

			// TODO: replace with no conversion step
			return returnValue;
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
	public double[] rootfreq(double when) {
		// TODO Auto-generated method stub
		return null;
	}



}
