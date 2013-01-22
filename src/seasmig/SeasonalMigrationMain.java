package seasmig;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.GregorianCalendar;

import mc3kit.FormattingLogger;
import mc3kit.MCMC;
import mc3kit.PowerHeatFunction;
import mc3kit.graphical.step.PriorLikelihoodOutputStep;
import mc3kit.graphical.step.SampleOutputStep;
import mc3kit.graphical.step.SwapStep;
import mc3kit.graphical.step.UnivariateStep;

import org.apache.log4j.ConsoleAppender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Layout;
import org.apache.log4j.Logger;
import org.apache.log4j.TTCCLayout;

import seasmig.Config.RunMode;
import treelikelihood.*;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.sun.xml.internal.ws.message.Util;

public class SeasonalMigrationMain
{
	public static void main(String[] args)
	{
		try
		{
			// Load config
			System.out.print("Loading config file...");
			Gson gson = new Gson();
			Config config = gson.fromJson(new FileReader("config.json"), Config.class);
			System.out.println(" done");
			System.out.print("Writing full config options to out.config...");
			config.outputToFile("out.config",gson);
			System.out.println(" done");

			if (config.runMode==RunMode.TEST1 || config.runMode==RunMode.TEST2) {
				// Roughly comparing results of several different exponentiation algorithms 
				testMatrixExponentiation(config.numLocations);
			}

			// Load data files and prepare data....			
			Data data = new Data(config);

			// Tests...			
			if (config.runMode==RunMode.TEST1 || config.runMode==RunMode.TEST2) {
				System.out.print("Running likelihood test...\n");
				testLikelihood(config, data);				
				System.out.print("Completed likelihood test!\n\n");

			}

			System.out.print("Initializing MCMC STEP 1...");
			MCMC mcmc = new MCMC();
			mcmc.setRandomSeed(config.randomSeed);
			mcmc.setChainCount(config.chainCount);
			System.out.println(" done");

			// Initialize logger
			System.out.print("Initializing logger...");
			Logger logger = Logger.getRootLogger();
			logger.setLevel(config.logLevel.getLog4jLevel());
			Layout layout = new TTCCLayout("ISO8601");
			if(config.logFilename.equals("-")) 	{
				logger.addAppender(new ConsoleAppender(layout));
			}
			else {
				logger.addAppender(new FileAppender(layout, config.logFilename));
			}
			FormattingLogger fmtLogger = new FormattingLogger(logger);

			String configJson = new GsonBuilder().setPrettyPrinting().create().toJson(config, Config.class);
			fmtLogger.info("Parsed config: %s", configJson);

			mcmc.setLogger(fmtLogger);

			System.out.println(" done");

			System.out.print("Initializing MCMC STEP 2...");

			// Set heating function for chains, interpolating from 1.0 to 0.0 
			// with an accelerating heating schedule
			PowerHeatFunction heatFunc = new PowerHeatFunction();
			heatFunc.setHeatPower(config.heatPower);
			heatFunc.setMinHeatExponent(0.0);
			mcmc.setHeatFunction(heatFunc);

			// Set model factory, which generates model instances to run on multiple chains
			SeasonalMigrationFactory modelFactory = new SeasonalMigrationFactory(config, data);
			mcmc.setModelFactory(modelFactory);

			// Step representing each chain going through each variable one at a time in random order
			UnivariateStep univarStep = new UnivariateStep();
			univarStep.setTuneFor(config.tuneFor);
			univarStep.setTuneEvery(config.tuneEvery);
			univarStep.setStatsFilename(config.varStatsFilename);
			univarStep.setRecordHeatedStats(config.recordHeatedStats);
			univarStep.setRecordStatsAfterTuning(true);

			// Differential evolution step
//			BlockDifferentialEvolutionStep deStep = new BlockDifferentialEvolutionStep();
//			deStep.setStatsFilename(config.demcStatsFilename);
//			deStep.setRecordHistoryEvery(config.thin);
//			deStep.setRecordHistoryAfter(config.recordHistoryAfter);
//			deStep.setTuneEvery(config.tuneEvery);
//			deStep.setTuneFor(config.tuneFor);
//			deStep.setInitialHistoryCount(config.initialHistoryCount);
//			deStep.setRecordHeatedStats(false);
//			deStep.setRecordStatsAfterTuning(true);

			// Swap steps: even (0,1), (2,3), etc. Odd (1,2), (3,4), etc.
			// Set up to allow parallelization of all swaps.
			SwapStep evenSwapStep = new SwapStep(SwapStep.SwapParity.EVEN);
			evenSwapStep.setStatsFilename(config.swapStatsFilename);
			evenSwapStep.setStatsEvery(config.tuneEvery * config.chainCount);
			SwapStep oddSwapStep = new SwapStep(SwapStep.SwapParity.ODD);
			oddSwapStep.setStatsFilename(config.swapStatsFilename);
			oddSwapStep.setStatsEvery(config.tuneEvery * config.chainCount);

			// Sample output step: all model parameters just for cold chain
			SampleOutputStep sampOutStep = new SampleOutputStep();
			sampOutStep.setChainId(0);
			sampOutStep.setFilename(config.sampleFilename);
			sampOutStep.setThin(config.thin);

			// Prior-likelihood output step: log-prior/log-likelihood for all chains
			PriorLikelihoodOutputStep plOutStep = new PriorLikelihoodOutputStep();
			plOutStep.setFilename(config.priorLikelihoodFilename);
			plOutStep.setThin(config.thin);

			// Set up execution order for the steps
			mcmc.addStep(univarStep);
//			mcmc.addStep(deStep);
			for(int i = 0; i < config.chainCount; i++)	{
				mcmc.addStep(evenSwapStep);
				mcmc.addStep(oddSwapStep);				
			}
			mcmc.addStep(sampOutStep);
			mcmc.addStep(plOutStep);

			System.out.println(" done");
			System.out.println("Running MCMC...");
			System.out.println("state\tlikelihood\thours/million states");
			mcmc.run();			
			System.out.println("done");
		}
		catch(Throwable e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static void testLikelihood(Config config, Data data) throws IOException {
		// Creating test file 
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		System.out.println("Calculating tree likelihood using the same model used to create the tree: SEASONALITY "+config.migrationSeasonality+",");				
		testStream.print("{\""+config.migrationSeasonality+"\",");
		testStream.print(data.createModel.parse());
		System.out.println(data.createModel.print());
		double createLikelihood = 0;
		for (LikelihoodTree tree : data.trees) {
			System.out.print(".");
			// TODO: maybe get likelihood to not require copy...
			LikelihoodTree workingCopy = tree.copy();
			workingCopy.setLikelihoodModel(data.createModel);
			createLikelihood+=workingCopy.logLikelihood();
		}
		createLikelihood=createLikelihood/data.trees.size();
		System.out.println(createLikelihood);

		System.out.println("\nCalculating tree likelihood using test models with increasing noise:");
		for (int i=0;i<data.testModels.size();i++) {
			if (i%config.numTestRepeats==0) {						
				System.out.println("SEASONALITY "+Config.Seasonality.values()[i/config.numTestRepeats]);						
			}

			double testLikelihood = 0;
			for (LikelihoodTree tree : data.trees) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setLikelihoodModel(data.testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.trees.size();
			System.out.println(testLikelihood);
		}
		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();
		
	}

	private static void testMatrixExponentiation(int n) {
		System.out.println("Testing matrix exponentiation:");

		for (int i=1;i<100;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);
			MatrixExponentiator test4 = new ReMolerMatrixExp(testMatrix);

			for (double t=0;t<500;t=(t+0.000001)*2) {
				String res1=treelikelihood.Util.parse(test1.expm(t));
				String res2=treelikelihood.Util.parse(test2.expm(t));
				String res3=treelikelihood.Util.parse(test3.expm(t));
				String res4=treelikelihood.Util.parse(test4.expm(t));
	
				if (res1.equalsIgnoreCase(res2) && res2.equalsIgnoreCase(res3) && res3.equalsIgnoreCase(res4)) {
					System.out.print(".");
				}
				else {
					System.out.println("----------------------------");
					System.out.println("MolerMatrixExp:"+treelikelihood.Util.print(test1.expm(t)));					
					System.out.println("Matlab7MatrixExp:"+treelikelihood.Util.print(test2.expm(t)));
					System.out.println("TaylorMatrixExp:"+treelikelihood.Util.print(test3.expm(t)));
					System.out.println("ReMolerMatrixExp:"+treelikelihood.Util.print(test4.expm(t)));
	
				}
			}
			if (i%10==0) System.out.println();
		}
		
		
		
		long startTime1= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res1=test1.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res1);
			}
		}
		long time1= System.currentTimeMillis()-startTime1;
		
		long startTime2= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
			for (int j=0;j<1200;j++) {
				double[][] res2=test2.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res2);
			}
		}
		long time2= System.currentTimeMillis()-startTime2;
		
		long startTime3= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res3=test3.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res3);
			}
		}
		long time3= System.currentTimeMillis()-startTime3;
		
		long startTime4= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test4 = new ReMolerMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res4=test4.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res4);
			}
		}
		long time4= System.currentTimeMillis()-startTime4;
		
		long startTime5= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res5=test5.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res5);
			}
		}
		long time5= System.currentTimeMillis()-startTime5;
		
		long startTime6= System.currentTimeMillis();		
		for (int i=1;i<200;i++) {
			double[][] testMatrix = Data.makeRandomMigrationMatrix(n,(double) i/100.0);
			MatrixExponentiator test6 = new JblasMatrixExp(testMatrix);;
			for (int j=0;j<1200;j++) {
				double[][] res6=test6.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
				if (Math.random()<0.0000000000001) System.out.println(res6);
			}
		}
		long time6= System.currentTimeMillis()-startTime6;
		
		System.out.println("----------------------------");
		System.out.println("MolerMatrixExp: "+time1+"ms");
		System.out.println("Matlab7MatrixExp: "+time2+"ms");
		System.out.println("TaylorMatrixExp: "+time3+"ms");
		System.out.println("ReMolerMatrixExp: "+time4+"ms");
		System.out.println("JblasMatrixExp: "+time5+"ms");
		System.out.println("JeblMatrixExp: "+time6+"ms");
		
		System.out.println("\nCompleted matrix exponentiation test");

	}
}
