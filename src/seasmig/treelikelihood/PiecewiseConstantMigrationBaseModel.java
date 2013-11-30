package seasmig.treelikelihood;

import java.util.HashMap;

import org.javatuples.Pair;

import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

@SuppressWarnings("serial")
public class PiecewiseConstantMigrationBaseModel implements MigrationBaseModel {
	// TODO: Check this...
	// TODO: Check zMult order... sould be ok for Q where rows sum to 1...

	// Precision Parameter
	static final double infinitesimalTime = 1E-5;

	// Cache Parameters 
	static final int maxCachedTransitionMatrices = 1600;

	// Precision Parameters
	int nYearParts;

	// Origin Model
	double[][][] seasonalRates;	
	DoubleFunction[] rootFreq;

	// Constant Migration Models
	MigrationBaseModel constantModels[];

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	HashMap<Pair<Double,Double>, double[][]> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, double[][]>();

	private int num_locations = 0;
	private double dt; 

	protected PiecewiseConstantMigrationBaseModel() {};

	// Constructor	
	public PiecewiseConstantMigrationBaseModel(double[][][] seasonalRates_, DoubleFunction[] rootFreq_, int nYearParts_) {	
		// TODO: Check this...
		// diagonal rates functions are calculated through row sums and are ignored...
		num_locations=seasonalRates_[0].length;	
		seasonalRates=seasonalRates_;
		nYearParts = nYearParts_;
		dt = 1.0/(double)nYearParts;
		constantModels = new ConstantMigrationBaseModel[nYearParts];

		for (int i=0;i<nYearParts;i++) {
			double[][] migrationMatrix = new double[num_locations][num_locations];
			for (int j=0; j<num_locations; j++) {
				double row_sum = 0;
				for (int k=0; k<num_locations; k++) {
					if (j!=k) {
						migrationMatrix[j][k]=seasonalRates[i][j][k];
						row_sum+=migrationMatrix[j][k];
					}
				}
				migrationMatrix[j][j]=-row_sum;
			}
			constantModels[i]=new ConstantMigrationBaseModel(migrationMatrix);
		}		
		rootFreq = rootFreq_;
	}

	// Constructor	
	public PiecewiseConstantMigrationBaseModel(double[][][] seasonalRates_, int nYearParts_) {	
		this(seasonalRates_,null,nYearParts_);
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
		double[][] cached = cachedTransitionMatrices.get(new Pair<Double,Double>(from_time_reminder,to_time_reminder));
		if (cached!=null) {
			return cached;
		}
		else {			
			// first step: 
			double step_start_time = from_time_reminder;
			double step_end_time = Math.min(to_time_reminder, Math.floor(step_start_time/dt)*dt+dt);
			DoubleMatrix2D result = F.identity(num_locations);	 

			while (step_start_time<to_time_reminder) {
				int yearPartIndex = Math.max(Math.min(0,(int) Math.floor(step_start_time%1.0/dt)),nYearParts-1);
				//assert yearPartIndex<nYearParts;
				// TODO: replace with other matrix mult
				result = result.zMult(DoubleFactory2D.dense.make(constantModels[yearPartIndex].transitionMatrix(step_start_time, step_end_time)),null);	
				step_start_time = step_end_time;
				step_end_time = Math.min(to_time_reminder, Math.floor((step_start_time+infinitesimalTime)/dt)*dt+dt);
			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			double[][] returnValue=result.toArray();
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time_reminder, to_time_reminder),returnValue);

			// TODO: replace with no conversion step
			return returnValue;
		}
	}

	@Override
	public String print() {
		String returnValue = "Picewise Constant Migration Model:\n";
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
		double[] returnValue = new double[num_locations];
		for (int i=0;i<num_locations;i++) {
			if (rootFreq!=null) 
				returnValue[i]=rootFreq[i].apply(when);		
			else 
				returnValue[i]=1.0/num_locations;
		}
		return returnValue;
	}

	@Override
	public Pair<Double, Integer> getNextRandomEvent(double time, int loc) {
		// TODO Auto-generated method stub
		return null;
	}



}

