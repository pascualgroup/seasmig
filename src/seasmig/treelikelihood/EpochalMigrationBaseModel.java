package seasmig.treelikelihood;

import java.util.HashMap;

import mc3kit.DoubleVariable;

import org.javatuples.Pair;

import seasmig.treelikelihood.MigrationBaseModel.Event;
import cern.colt.function.DoubleFunction;
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

@SuppressWarnings("serial")
public class EpochalMigrationBaseModel implements MigrationBaseModel {

	// Precision Parameter
	static final double infinitesimalTime = 1E-5;

	// Cache Parameters 
	static final int maxCachedTransitionMatrices = 1600;

	// Precision Parameters
	int nParts;
	double[] epochs;
	DoubleVariable[] epochTimes;

	// Origin Model
	double[][][] epochRates;	
	DoubleFunction[] rootFreq;

	// Constant Migration Models
	MigrationBaseModel constantModels[];

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>,DoubleMatrix2D>();

	private int num_locations = 0;

	protected EpochalMigrationBaseModel() {};

	// Constructor	
	public EpochalMigrationBaseModel(double[][][] seasonalRates_, DoubleFunction[] rootFreq_, double[] epochs_) {	
		// TODO: Check this...
		// diagonal rates functions are calculated through row sums and are ignored...
		num_locations=seasonalRates_[0].length;	
		epochRates=seasonalRates_;
		epochs = epochs_;
		nParts = epochs.length;
		constantModels = new ConstantMigrationBaseModel[nParts];

		for (int i=0;i<nParts;i++) {
			double[][] migrationMatrix = new double[num_locations][num_locations];
			for (int j=0; j<num_locations; j++) {
				double row_sum = 0;
				for (int k=0; k<num_locations; k++) {
					if (j!=k) {
						migrationMatrix[j][k]=epochRates[i][j][k];
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
	public EpochalMigrationBaseModel(double[][][] seasonalRates_, double[] epochs_) {	
		this(seasonalRates_, null, epochs_);
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
		DoubleMatrix2D cached = cachedTransitionMatrices.get(new Pair<Double,Double>(from_time,to_time));
		if (cached!=null) {
			return cached;
		}
		else {			
			// first step: 
			double step_start_time = from_time;
			int epochIndex = epochIndex(from_time);
			double step_end_time = Math.min(to_time, epochEndTime(epochIndex));
			DoubleMatrix2D result = F.identity(num_locations);	 

			while (step_start_time<to_time) {
				// TODO: replace with other matrix mult
				result = result.zMult(constantModels[epochIndex].transitionMatrix(step_start_time, step_end_time),null);	
				step_start_time = step_end_time;
				epochIndex=epochIndex+1;
				step_end_time = Math.min(to_time, epochEndTime(epochIndex));
			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time, to_time),result);

			// TODO: replace with no conversion step
			return result;
		}
	}

	private int epochIndex(double from_time) {
		for (int i=0;i<nParts-1;i++) {
			if (epochs[i]>from_time) {
				return i;
			}
		}
		return nParts-1;
	}

	private double epochEndTime(int index) {
		if (index<(nParts-1))
			return epochs[index];
		else
			return Double.MAX_VALUE;
	}

	@Override
	public String print() {
		String returnValue = "Epochal Migration Model:\n";
		returnValue+="[";
		for (int i=0;i<epochRates.length;i++) {
			if (i!=0) returnValue+=" ";
			returnValue+="[";
			for (int j=0;j<epochRates[i].length;j++) {
				if (i!=j) 
					returnValue=returnValue+String.format("%40s",epochRates[i][j].toString());
				else 
					returnValue+=String.format("%40s", "NA");
				if (j!=epochRates[i].length-1) returnValue+=",";
			}
			returnValue+="]";
			if (i!=epochRates.length-1) returnValue+="\n";
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
			if (rootFreq!=null) 
				returnValue[i]=rootFreq[i].apply(when);		
			else 
				returnValue[i]=1.0/num_locations;
		}
		return DoubleFactory1D.dense.make(returnValue);
	}

	@Override
	public Event nextEvent(double time, int from) {
		// TODO: check this...
		Event nextEvent = null;
		boolean done = false;
		double currentTime = time;
		int currentLoc = from;
		do {
			int currentEpochIndex = epochIndex(currentTime);
			nextEvent = constantModels[currentEpochIndex].nextEvent(currentTime, currentLoc);
			int nextEpochIndex = epochIndex(nextEvent.time);
			done = (nextEpochIndex==currentEpochIndex);
			if (!done) {
				currentTime = epochs[nextEpochIndex];				
			}									
		} while (!done);
		return nextEvent;
	}

}

