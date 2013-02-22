package seasmig;

import org.junit.Test;


public class SeasonalMigrationOptimizeTest {
	// TODO: this...

	@Test
	public void test() {
		NelderMead myOptimize = new NelderMead();

		myOptimize.NM(myOptimize.makeInitPoints(15, 15));				

	}


	
}
