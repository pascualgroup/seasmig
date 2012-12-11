package treelikelihood;

import java.util.HashMap;
import java.util.Vector;

import cern.colt.matrix.tdouble.*;
import cern.colt.matrix.tdouble.DoubleMatrix2D;



public class TwoSeasonalMigrationBaseModel implements MigrationBaseModel {

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
	HashMap<Double, DoubleMatrix2D> cachedTransitionMatrices = new HashMap<Double, DoubleMatrix2D>();

	private int num_states = 0;

	// Constructor	
	public TwoSeasonalMigrationBaseModel(double[][] Q1_,double[][] Q2_, double season1Start_, double season1End_) {		
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

		DoubleMatrix2D cached = cachedTransitionMatrices.get(to_time-from_time);
		if (cached!=null) {
			return cached;
		}
		else {
			double start_time = from_time;
			double end_time = start_time;
			DoubleMatrix2D result = F.identity(num_states);

			int n=0;
			while (end_time<to_time) {
				n=n+1;
				if (isInSeason1(start_time)) {
					end_time = Math.min(to_time,start_time+season1Length);
					result = result.zMult(season1MigrationModel.transitionMatrix(start_time, end_time),null);	
					start_time=end_time;				
				} else {
					end_time = Math.min(to_time,start_time+season2Length);
					result = result.zMult(season2MigrationModel.transitionMatrix(start_time, end_time),null);	
					start_time=end_time;			
				}

			}

			// cache result
			if (cachedTransitionMatrices.size()>=maxCachedTransitionMatrices) {
				cachedTransitionMatrices.remove(cachedTransitionMatrices.keySet().iterator().next());
			}			
			cachedTransitionMatrices.put(to_time-from_time, result);

			return result;
		}
	}

	@Override
	public String print() {		
		return season1MigrationModel.print()+","+season2MigrationModel.print();
	}

	private boolean isInSeason1(double time) {
		return (time%1.0>=season1Start) && ((time-season1Start)%1.0<season1Length);			
	}

	@Override
	public int getNumStates() {
		return num_states ;
	}



}
