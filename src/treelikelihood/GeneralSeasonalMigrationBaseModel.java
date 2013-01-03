package treelikelihood;

import java.util.HashMap;
import java.util.Vector;

import org.javatuples.Pair;

import corejava.Format;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.tdouble.*;



public class GeneralSeasonalMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...
	
	// Cache Parameters 
	static final int maxCachedTransitionMatrices = 16000;

	// Precision Parameters
	static final int nYearParts = 12;
	
	// Origin Model
	DoubleFunction[][] seasonalRates;
	
	// Constant Migration Models
	MigrationBaseModel constantModels[] = new ConstantMigrationBaseModel[nYearParts];

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	Vector<DoubleMatrix2D> cachedMatrixPower = new Vector<DoubleMatrix2D>();	
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, DoubleMatrix2D>();

	private int num_locations = 0;
	private double dt = 1.0/(double)nYearParts; 

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
		return Math.log(transitionMatrix(from_time, to_time).get(from_location, to_location));
	}

	@Override
	public DoubleMatrix2D transitionMatrix(double from_time, double to_time) {
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
			double step_end_time = Math.min(to_time, Math.max(0,Math.ceil(step_start_time/dt)*dt));
			DoubleMatrix2D result = F.identity(num_locations);	 
			
			while (step_end_time<to_time) {
				int yearPartIndex = (int) Math.floor(step_start_time%1.0/dt);
				result = result.zMult(constantModels[yearPartIndex].transitionMatrix(step_start_time, step_end_time),null);	
				step_start_time = step_end_time;
				step_end_time = Math.min(to_time, step_start_time+dt);
			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time_reminder, to_time_reminder),result);

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



}
