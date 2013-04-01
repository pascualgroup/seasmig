package seasmig.treelikelihood.matrixexp;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.GregorianCalendar;
import jebl.evolution.io.ImportException;

import org.junit.Test;

import seasmig.Config;
import seasmig.Data;
import seasmig.data.DataForTests;
import seasmig.data.DataForTests.TestType;
import seasmig.treelikelihood.*;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.JamaMolerMatrixExp;
import seasmig.treelikelihood.matrixexp.JblasMatrixExp;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;
import seasmig.treelikelihood.matrixexp.MolerMatrixExp;
import seasmig.treelikelihood.matrixexp.TaylorMatrixExp;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class MatrixExpTest {
// TODO: Implement this...	
	
	private static int numTestRepeats = 100;
	private static int numLocations = 3;
	
	@Test
	public void testMatrixExponentiation() {
		System.out.println("Testing matrix exponentiation:");

		for (int i=1;i<100;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);
			MatrixExponentiator test4 = new JamaMolerMatrixExp(testMatrix);
			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);
			MatrixExponentiator test6 = null;
			if (numLocations==3) {				
				test6 = new AnalyticMatrixExp3(testMatrix);
				if (!test6.checkMethod())
					test6=null;
			}
			if (numLocations==2) {
				test6 = new AnalyticMatrixExp2(testMatrix);
				if (!test6.checkMethod())
					test6=null;
			}

			for (double t=0;t<500;t=(t+0.1)*2) {
				String res1=seasmig.util.Util.parse(test1.expm(t));
				String res2=seasmig.util.Util.parse(test2.expm(t));
				String res3=seasmig.util.Util.parse(test3.expm(t));
				String res4=seasmig.util.Util.parse(test4.expm(t));
				String res5=seasmig.util.Util.parse(test5.expm(t));
				String res6 = null;
				if (test6!=null) res6 = seasmig.util.Util.parse(test6.expm(t));

				boolean resultsMatch=(res1.equalsIgnoreCase(res2) && res2.equalsIgnoreCase(res3) && res3.equalsIgnoreCase(res4)  && res4.equalsIgnoreCase(res5));
				if (test6!=null) {
					resultsMatch=resultsMatch&&(res6.equalsIgnoreCase(res1));
				}
				if (resultsMatch) {
					System.out.print(".");
				}
				else {
					System.out.println("\nMatrix exponentiation results fail to match!\n t="+t+" \nQ=\n"+seasmig.util.Util.print(testMatrix));
					System.out.println("MolerMatrixExp - checkMethod: "+test1.checkMethod()+"\n"+seasmig.util.Util.print(test1.expm(t)));					
					System.out.println("Matlab7MatrixExp - checkMethod: "+test2.checkMethod()+"\n"+seasmig.util.Util.print(test2.expm(t)));
					System.out.println("TaylorMatrixExp - checkMethod: "+test3.checkMethod()+"\n"+seasmig.util.Util.print(test3.expm(t)));
					System.out.println("JamaMolerMatrixExp: - checkMethod: "+test4.checkMethod()+"\n"+seasmig.util.Util.print(test4.expm(t)));
					System.out.println("JblasMatrixExp: - checkMethod: "+test5.checkMethod()+"\n"+seasmig.util.Util.print(test5.expm(t)));
					if (numLocations==3) System.out.println("AnalyticMatrixExp3: - checkMethod: "+test6.checkMethod()+"\n"+seasmig.util.Util.print(test6.expm(t)));;
					if (numLocations==2) System.out.println("AnalyticMatrixExp2: - checkMethod: "+test6.checkMethod()+"\n"+seasmig.util.Util.print(test6.expm(t)));;

				}
			}
			if (i%10==0) System.out.println();
		}

		long startTime1= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res1=test1.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res1);
			}
		}
		long time1= System.currentTimeMillis()-startTime1;

		long startTime2= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res2=test2.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res2);
			}
		}
		long time2= System.currentTimeMillis()-startTime2;

		long startTime3= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res3=test3.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res3);
			}
		}
		long time3= System.currentTimeMillis()-startTime3;

		long startTime4= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test4 = new JamaMolerMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res4=test4.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res4);
			}
		}
		long time4= System.currentTimeMillis()-startTime4;

		long startTime5= System.currentTimeMillis();		
		for (int i=1;i<50;i++) {
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res5=test5.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res5);
			}
		}
		long time5= System.currentTimeMillis()-startTime5;
		
		long time6 = 0;
		if (numLocations==3) {
			long startTime6= System.currentTimeMillis();		
			for (int i=1;i<50;i++) {
				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
				MatrixExponentiator test6 = new AnalyticMatrixExp3(testMatrix);;
				for (int j=0;j<1200;j++) {
					double[][] res6=test6.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
					if (Math.random()<0.0000000000001) System.out.println(res6);
				}
			}
			time6= System.currentTimeMillis()-startTime6;
		}
		if (numLocations==2) {
			long startTime6= System.currentTimeMillis();		
			for (int i=1;i<50;i++) {
				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
				MatrixExponentiator test6 = new AnalyticMatrixExp2(testMatrix);;
				for (int j=0;j<1200;j++) {
					double[][] res6=test6.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
					if (Math.random()<0.0000000000001) System.out.println(res6);
				}
			}
			time6= System.currentTimeMillis()-startTime6;
		}
	

		System.out.println("\nMolerMatrixExp: "+time1+"ms");
		System.out.println("Matlab7MatrixExp: "+time2+"ms");
		System.out.println("TaylorMatrixExp: "+time3+"ms");
		System.out.println("JamaMolerMatrixExp: "+time4+"ms");
		System.out.println("JblasMatrixExp: "+time5+"ms");
		if (numLocations==3) System.out.println("AnalyticMatrixExp3: "+time6+"ms");
		if (numLocations==2) System.out.println("AnalyticMatrixExp2: "+time6+"ms");

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
//		Data data = new DataForTests(config,)
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
		Data data = new DataForTests(config,TestType.TEST_MODEL_DEGENERACY,3,3,10);

		// Creating test file 
		System.out.println("Calculating tree likelihood using degenerate models:");				
        double[] results = new double[((DataForTests) data).testModels.size()];

		for (int i=0;i<((DataForTests) data).testModels.size();i++) {
			System.out.println("SEASONALITY "+((DataForTests) data).testModels.get(i).getModelName());						
			long startTime= System.currentTimeMillis();	
			double testLikelihood = 0;
			for (LikelihoodTree tree : data.getTrees().get(0)) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setLikelihoodModel(((DataForTests) data).testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.getTrees().size();
			System.out.println(testLikelihood);
			results[i]=testLikelihood;
			long duration= System.currentTimeMillis()-startTime;
			System.out.println("duration: "+duration+"ms");
		}
				
		for (int i=1;i<results.length;i++) {
			assertEquals(results[i],results[i-1], 1E-3);			
		}
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
		Data data = new DataForTests(config,TestType.TEST_USING_GENERATED_TREES,10,3,5);

		// Creating test file 
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		System.out.println("Calculating tree likelihood using the same model used to create the tree: SEASONALITY "+config.migrationSeasonality+",");				
		testStream.print("{\""+config.migrationSeasonality+"\",");
		testStream.print(((DataForTests)data).createModel.parse());
		System.out.println(((DataForTests)data).createModel.print());
		double createLikelihood = 0;
		for (LikelihoodTree tree : data.getTrees().get(0)) {
			System.out.print(".");
			// TODO: maybe get likelihood to not require copy...
			LikelihoodTree workingCopy = tree.copy();
			workingCopy.setLikelihoodModel(((DataForTests)data).createModel);
			createLikelihood+=workingCopy.logLikelihood();
		}
		createLikelihood=createLikelihood/data.getTrees().size();
		System.out.println(createLikelihood);

		System.out.println("\nCalculating tree likelihood using test models with increasing noise:");
		for (int i=0;i<((DataForTests)data).testModels.size();i++) {
			if (i%numTestRepeats ==0) {						
				System.out.println("SEASONALITY "+((DataForTests)data).testModels.get(i).getModelName());						
			}

			double testLikelihood = 0;
			for (LikelihoodTree tree : data.getTrees().get(0)) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setLikelihoodModel(((DataForTests)data).testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.getTrees().size();
			System.out.println(testLikelihood);
		}
		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();

	}

	
}
