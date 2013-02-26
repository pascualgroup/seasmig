package optimize;
// THIS IS A PRESKETCH!!!

import static org.junit.Assert.assertEquals;

import java.util.Vector;

import optimize.NelderMead.SimpleDoubleFunction;

import org.javatuples.Pair;
import org.junit.Test;


public class SeasonalMigrationOptimizeTest {
	// TODO: this...

	public class MyFunc implements NelderMead.SimpleDoubleFunction {
		public double eval(Vector<Double> x) {
			double returnValue = 0; 
			for (int i=0;i<x.size();i++) 
				returnValue+=Math.pow(Math.abs(x.get(i)-i),i);
			return returnValue;
		}
	}
	
	@Test
	public void test() {
		
		SimpleDoubleFunction myFunc = new MyFunc();
		Vector<Vector<Double>> x0s = new Vector<Vector<Double>>(10);
		for (int i=0;i<10;i++) {
			x0s.add(new Vector<Double>(5));
			for (int j=0;j<5;j++) {
				x0s.get(i).add(Math.random());
			}
		}
		
		NelderMead myOptimize = new NelderMead(myFunc, x0s);

		double f_goal = 0;		
		Pair<Vector<Double>, Double> result = myOptimize.optimize();
		System.err.println("goal: [1 1 1...1]\nf(goal):"+f_goal);
		System.err.println("optimization result: "+result.getValue0()+"\nf(result): "+result.getValue1());
		
		assertEquals(result.getValue1(),f_goal,0.0001);
	}


	
}
