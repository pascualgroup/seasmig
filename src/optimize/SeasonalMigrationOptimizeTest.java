package optimize;

import static org.junit.Assert.assertEquals;

import java.util.Vector;

import org.javatuples.Pair;
import org.junit.Test;


public class SeasonalMigrationOptimizeTest {
	// TODO: this...

	@Test
	public void test() {
		NelderMead myOptimize = new NelderMead();

		double f_goal = 0;		
		Pair<Vector<Double>, Double> result = myOptimize.optimize();
		System.err.println("goal: [1 1 1...1]\nf(goal):"+f_goal);
		System.err.println("optimization result: "+result.getValue0()+"\nf(result): "+result.getValue1());
		
		assertEquals(result.getValue1(),f_goal,0.0001);
	}


	
}
