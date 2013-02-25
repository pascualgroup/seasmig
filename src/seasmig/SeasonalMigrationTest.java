package seasmig;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.GregorianCalendar;
import jebl.evolution.io.ImportException;

import org.junit.Test;

import seasmig.TestData.TestType;
import seasmig.treelikelihood.*;
import seasmig.treelikelihood.matrixexp.JamaMolerMatrixExp;
import seasmig.treelikelihood.matrixexp.JblasMatrixExp;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;
import seasmig.treelikelihood.matrixexp.MolerMatrixExp;
import seasmig.treelikelihood.matrixexp.TaylorMatrixExp;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class SeasonalMigrationTest {
// TODO: Implement this...	
	
	private static int numTestRepeats = 100;
	private static int numLocations = 3;
	
	@Test
	public void testMatrixExponentiation() {
		System.out.println("Testing matrix exponentiation:");

		for (int i=1;i<100;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);
			MatrixExponentiator test4 = new JamaMolerMatrixExp(testMatrix);
			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);

			for (double t=0;t<500;t=(t+0.1)*2) {
				String res1=seasmig.util.Util.parse(test1.expm(t));
				String res2=seasmig.util.Util.parse(test2.expm(t));
				String res3=seasmig.util.Util.parse(test3.expm(t));
				String res4=seasmig.util.Util.parse(test4.expm(t));
				String res5=seasmig.util.Util.parse(test5.expm(t));

				if (res1.equalsIgnoreCase(res2) && res2.equalsIgnoreCase(res3) && res3.equalsIgnoreCase(res4)  && res4.equalsIgnoreCase(res5)) {
					System.out.print(".");
				}
				else {
					System.out.println("\nMolerMatrixExp: "+seasmig.util.Util.print(test1.expm(t)));					
					System.out.println("Matlab7MatrixExp: "+seasmig.util.Util.print(test2.expm(t)));
					System.out.println("TaylorMatrixExp: "+seasmig.util.Util.print(test3.expm(t)));
					System.out.println("JamaMolerMatrixExp: "+seasmig.util.Util.print(test4.expm(t)));
					System.out.println("JblasMatrixExp: "+seasmig.util.Util.print(test5.expm(t)));

				}
			}
			if (i%10==0) System.out.println();
		}

		long startTime1= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res1=test1.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res1);
			}
		}
		long time1= System.currentTimeMillis()-startTime1;

		long startTime2= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res2=test2.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res2);
			}
		}
		long time2= System.currentTimeMillis()-startTime2;

		long startTime3= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res3=test3.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res3);
			}
		}
		long time3= System.currentTimeMillis()-startTime3;

		long startTime4= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test4 = new JamaMolerMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res4=test4.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res4);
			}
		}
		long time4= System.currentTimeMillis()-startTime4;

		long startTime5= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = TestData.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res5=test5.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res5);
			}
		}
		long time5= System.currentTimeMillis()-startTime5;

		System.out.println("\nMolerMatrixExp: "+time1+"ms");
		System.out.println("Matlab7MatrixExp: "+time2+"ms");
		System.out.println("TaylorMatrixExp: "+time3+"ms");
		System.out.println("JamaMolerMatrixExp: "+time4+"ms");
		System.out.println("JblasMatrixExp: "+time5+"ms");

		System.out.println("\nCompleted matrix exponentiation test");

	}

//	@Test
//	public void testMain() throws Throwable
//	{
//		// Load config
//		System.out.print("Loading config file... ");
//		Gson gson = new GsonBuilder().setPrettyPrinting().create();
//		Config config = null;
//		try {
//			config = gson.fromJson(new FileReader("config.json"), Config.class);
//			System.out.println(" done");
//		}
//		catch(Throwable e)	{
//			config=new Config();
//			System.out.println("config.json file not found, using default config. See out.config.json for details");			
//		}			
//
//		System.out.print("Writing full config options to out.config...");
//		config.outputToFile("out.config",gson);
//		System.out.println(" done");
//
//		// Load data files and prepare data....			
//		Data data = new TestData(config,)
//
//		// Setup MCMC
//		System.out.print("Setting up MCMC....");
//		MCMC mcmc = new MCMC();
//		mcmc.setRandomSeed(config.randomSeed); // TODO: check that this changes...
//
//		MCMC.setLogLevel(config.logLevel);
//
//		ModelFactory mf = new SeasonalMigrationFactory(config, data) {
//			@Override
//			public Model createModel(Chain initialChain) throws MC3KitException {
//				Model m = new Model(initialChain);
//
//				m.beginConstruction();
//				new IntVariable(m, "v", new UniformIntDistribution(m, 1,4));
//				m.endConstruction();
//
//				return m;
//			}
//		};
//
//		mcmc.setModelFactory(mf);
//
//		UnivariateProposalStep proposalStep = new UnivariateProposalStep(0.25, 100, config.burnIn);
//		mcmc.addStep(proposalStep);
//
//		System.out.println(" done!");
//		// Run, collect statistics, and check moments against expected distribution
//
//		System.out.println("Running MCMC...");
//
//		System.out.println("state\tlikelihood\thours/million states");
//		double sum = 0;
//		for(long i = 0; i < config.iterCount; i++) {
//			mcmc.step();
//			mcmc.getModel().recalculate();
//
//			assertEquals(i + 1, mcmc.getIterationCount());
//
//			if(i >= config.burnIn) {
//				int val = mcmc.getModel().getIntVariable("v").getValue();
//				sum += val;
//			}
//		}
//
//		double N = config.iterCount - config.burnIn;
//		System.out.println(sum/N);		
//		System.out.println("done!");	
//	}

	@Test
	public void testModelDegeneracy() throws IOException, ImportException {
		// Load config
		System.out.print("Loading config file... ");
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		Config config = null;
		try {
			config = gson.fromJson(new FileReader("config.json"), Config.class);
			System.out.println(" done");
		}
		catch(Throwable e)	{
			config=new Config();
			System.out.println("config.json file not found, using default config. See out.config.json for details");			
		}			

		System.out.print("Writing full config options to out.config...");
		config.outputToFile("out.config",gson);
		System.out.println(" done");

		// Load data files and prepare data....			
		Data data = new TestData(config,TestType.TEST_MODEL_DEGENERACY,3,3,10);

		// Creating test file 
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		System.out.println("Calculating tree likelihood using degenerate models:");				


		for (int i=0;i<((TestData) data).testModels.size();i++) {
			System.out.println("SEASONALITY "+((TestData) data).testModels.get(i).getModelName());						
			long startTime= System.currentTimeMillis();	
			double testLikelihood = 0;
			for (LikelihoodTree tree : data.getTrees()) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setLikelihoodModel(((TestData) data).testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.getTrees().size();
			System.out.println(testLikelihood);
			long duration= System.currentTimeMillis()-startTime;
			System.out.println("duration: "+duration+"ms");
		}

		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();
		System.exit(0);

	}

	@Test
	public void testLikelihood() throws IOException, ImportException {

		// Load config
		System.out.print("Loading config file... ");
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		Config config = null;
		try {
			config = gson.fromJson(new FileReader("config.json"), Config.class);
			System.out.println(" done");
		}
		catch(Throwable e)	{
			config=new Config();
			System.out.println("config.json file not found, using default config. See out.config.json for details");			
		}			

		System.out.print("Writing full config options to out.config...");
		config.outputToFile("out.config",gson);
		System.out.println(" done");

		// Load data files and prepare data....			
		Data data = new TestData(config,TestType.TEST_USING_GENERATED_TREES,10,3,5);

		// Creating test file 
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		System.out.println("Calculating tree likelihood using the same model used to create the tree: SEASONALITY "+config.migrationSeasonality+",");				
		testStream.print("{\""+config.migrationSeasonality+"\",");
		testStream.print(((TestData)data).createModel.parse());
		System.out.println(((TestData)data).createModel.print());
		double createLikelihood = 0;
		for (LikelihoodTree tree : data.getTrees()) {
			System.out.print(".");
			// TODO: maybe get likelihood to not require copy...
			LikelihoodTree workingCopy = tree.copy();
			workingCopy.setLikelihoodModel(((TestData)data).createModel);
			createLikelihood+=workingCopy.logLikelihood();
		}
		createLikelihood=createLikelihood/data.getTrees().size();
		System.out.println(createLikelihood);

		System.out.println("\nCalculating tree likelihood using test models with increasing noise:");
		for (int i=0;i<((TestData)data).testModels.size();i++) {
			if (i%numTestRepeats ==0) {						
				System.out.println("SEASONALITY "+((TestData)data).testModels.get(i).getModelName());						
			}

			double testLikelihood = 0;
			for (LikelihoodTree tree : data.getTrees()) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setLikelihoodModel(((TestData)data).testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.getTrees().size();
			System.out.println(testLikelihood);
		}
		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();

	}

	
}
