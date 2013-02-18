package seasmig;


public class SeasonalMigrationOptimizeMain {

	public static void main(String[] args) {
		NelderMead myOptimize = new NelderMead();

		myOptimize.NM(myOptimize.makeInitPoints(15, 15));				

	}


	
}
