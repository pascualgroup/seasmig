package treelikelihood;

import java.util.HashMap;
import java.util.Vector;

import org.javatuples.Pair;

import cern.colt.matrix.tdouble.*;



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

	// Caching
	DoubleFactory2D F = DoubleFactory2D.dense;
	Vector<DoubleMatrix2D> cachedMatrixPower = new Vector<DoubleMatrix2D>();	
	HashMap<Pair<Double,Double>, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Pair<Double,Double>, DoubleMatrix2D>();

	private int num_states = 0;

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
		return Math.log(transitionMatrix(from_time, to_time).get(from_state, to_state));
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
			double step_start_time = from_time;
			double step_end_time = step_start_time;
			DoubleMatrix2D result = F.identity(num_states);

			int n=0;
			while (step_end_time<to_time) {
				n=n+1;
				if (isInSeason1(step_start_time)) {
					step_end_time = Math.min(to_time,step_start_time+season1Length);
					result = result.zMult(season1MigrationModel.transitionMatrix(step_start_time, step_end_time),null);	
					step_start_time=step_end_time;				
				} else {
					step_end_time = Math.min(to_time,step_start_time+season2Length);
					result = result.zMult(season2MigrationModel.transitionMatrix(step_start_time, step_end_time),null);	
					step_start_time=step_end_time;			
				}

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
		return "\nphase: "+season1Start+" length: "+season1Length+"\n"+season1MigrationModel.print()+"\n"+season2MigrationModel.print();
	}

	private boolean isInSeason1(double time) {
		return (time%1.0>=season1Start) && ((time-season1Start)%1.0<season1Length);			
	}

	@Override
	public int getNumStates() {
		return num_states ;
	}



}
