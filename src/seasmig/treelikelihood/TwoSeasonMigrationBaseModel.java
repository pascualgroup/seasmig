package seasmig.treelikelihood;

import java.util.HashMap;
import java.util.Vector;

import org.javatuples.Pair;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

//TODO: test multiplication order... ok for Q where rows sum to 1


@SuppressWarnings("serial")
public class TwoSeasonMigrationBaseModel implements MigrationBaseModel {

	// Cache Parameters
	static final int maxCachedTransitionMatrices = 16000;

	// Seasonal Migration Models
	MigrationBaseModel season1MigrationModel = null;
	MigrationBaseModel season2MigrationModel = null;

	// Seasonal Range [...season2....|...season1......|...season2..]
	//                0             S1S              S1E        1 year
	// season1 CAN NOT cross year borders...

	double season1Start = 0;
	double season1Length = 0;
	double season2Length = 0;
	
	final static double infitesimalTimeInterval = 1E-6;

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	Vector<DoubleMatrix2D> cachedMatrixPower = new Vector<DoubleMatrix2D>();	
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, DoubleMatrix2D>();

	private int num_states = 0;

	
	
	protected TwoSeasonMigrationBaseModel() {};

	// Constructor	
	public TwoSeasonMigrationBaseModel(double[][] Q1_,double[][] Q2_, double season1Start_, double season1End_) {		
		season1Start=season1Start_;		
		season1Length=season1End_-season1Start_;
		season2Length=1-season1Length;
		season1MigrationModel=new ConstantMigrationBaseModel(Q1_);
		season2MigrationModel=new ConstantMigrationBaseModel(Q2_);
		num_states=Q1_.length;
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
		// TODO: remove make
		double from_time_reminder = from_time % 1.0;
		double from_time_div = from_time - from_time_reminder;
		double to_time_reminder = to_time - from_time_div;
		DoubleMatrix2D cached = cachedTransitionMatrices.get(new Pair<Double,Double>(from_time_reminder,to_time_reminder));
		if (cached!=null) {
			return cached;
		}
		else {
			double step_start_time = from_time;
			double step_end_time = step_start_time;
			DoubleMatrix2D result = F.identity(num_states);

			while (step_start_time<to_time) {
				double step_start_time_reminder=step_start_time%1.0;
				double step_start_time_div=step_start_time - step_start_time_reminder;
				if (isInSeason1(step_start_time)) {					
					step_end_time = Math.min(to_time,step_start_time_div+season1Start+season1Length+ infitesimalTimeInterval);
					result = result.zMult(season1MigrationModel.transitionMatrix(step_start_time, step_end_time),null);									
				} 
				else { // In Season 2 		
					if (step_start_time_reminder<season1Start) {
						step_end_time = Math.min(to_time,step_start_time_div+season1Start+ infitesimalTimeInterval);
					}
					else {
						step_end_time = Math.min(to_time,step_start_time_div+1.0+season1Start+ infitesimalTimeInterval);
					}
					result = result.zMult(season2MigrationModel.transitionMatrix(step_start_time, step_end_time),null);								
				}
				step_start_time=step_end_time;	

			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				for (int i=0;i<cachedTransitionMatrices.size()/2;i++) {
					cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
				}
			}			
			cachedTransitionMatrices.put(new Pair<Double,Double>(from_time_reminder, to_time_reminder),result);
			
			return result;
		}
	}

	@Override
	public String print() {	
		String returnValue="";
		returnValue+="{phase,length,rates1,rates2},\n";
		returnValue+="{"+season1Start+","+season1Length+",\n"+season1MigrationModel.print()+","+season2MigrationModel.print()+"}";
		return returnValue;
	}
	
	@Override
	public String parse() {	
		String returnValue="";
		returnValue+="{\"phase\",\"length\",\"rates1\",\"rates2\"}\n";
		returnValue+="{"+season1Start+","+season1Length+"}\n"+season1MigrationModel.parse()+"\n"+season2MigrationModel.parse()+"\n";
		return returnValue;
	}

	private boolean isInSeason1(double time) {
		return (time%1.0>=season1Start) && ((time-season1Start)%1.0<season1Length);			
	}

	@Override
	public int getNumLocations() {
		return num_states ;
	}

	@Override
	public String getModelName() {		
		return "Two Seasons";
	}

	@Override
	public DoubleMatrix1D rootfreq(double when) {
		if (isInSeason1(when))
			return season1MigrationModel.rootfreq(when);
		else
			return season2MigrationModel.rootfreq(when);
	}

}
