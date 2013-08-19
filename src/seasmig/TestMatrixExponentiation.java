package seasmig;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.GregorianCalendar;
import java.util.Vector;

import jebl.evolution.io.ImportException;
import mc3kit.ChainParity;
import mc3kit.MCMC;
import mc3kit.Step;
import mc3kit.SwapStep;
import mc3kit.VerificationStep;
import mc3kit.monitoring.MarginalLikelihoodStep;
import mc3kit.output.PriorLikelihoodOutputStep;
import mc3kit.output.SampleOutputStep;
import mc3kit.proposal.DEMCProposalStep;
import mc3kit.proposal.UnivariateProposalStep;

import org.junit.Test;

import seasmig.data.DataForTests;
import seasmig.data.DataForTests.TestType;
import seasmig.models.SeasonalMigrationModelFactory;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.JamaMolerMatrixExp;
import seasmig.treelikelihood.matrixexp.JblasMatrixExp;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;
import seasmig.treelikelihood.matrixexp.MolerMatrixExp;
import seasmig.treelikelihood.matrixexp.TaylorMatrixExp;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class TestMatrixExponentiation {
	// TODO: Organize this....

	@Test
	public void testMatrixExponentiation() {
		final int numLocations = 3;
		final int numScaleSteps = 100;
		final double maxMatrixScale = 10.0;
		final double minMatrixScale = 0.001;
		final double minTime = 0.001;
		final double maxTime = 100.0;
		final double numTimeSteps = 100.0;
		final double tol = 0.001;
		System.out.println("Testing matrix exponentiation: nLocations="+numLocations+ "tol="+tol);
		System.out.println("minMatrixScale: "+minMatrixScale+" maxMatrixScale: "+maxMatrixScale+" steps: "+numScaleSteps);
		System.out.println("minTime: "+minTime+" maxTime: "+maxTime+" numTimesteps: "+numTimeSteps);
		System.out.println("Comparing expm(Q*t)");
		int dotIter =0;
		for (int scaleIter=0;scaleIter<100;scaleIter++) {
			dotIter+=1;
			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) scaleIter*(maxMatrixScale-minMatrixScale)/(double)numScaleSteps);
			Vector<MatrixExponentiator> tests = new Vector<MatrixExponentiator>();
			tests.add(new MolerMatrixExp(testMatrix));

			tests.add(new Matlab7MatrixExp(testMatrix));
			tests.add(new TaylorMatrixExp(testMatrix));
			tests.add(new JamaMolerMatrixExp(testMatrix));
			tests.add(new JblasMatrixExp(testMatrix));
			if (numLocations==3) 				
				tests.add(new AnalyticMatrixExp3(testMatrix));
			if (numLocations==2) 
				tests.add(new AnalyticMatrixExp2(testMatrix));

			for (double t=minTime;t<maxTime;t+=(maxTime-minTime)/numTimeSteps) {
				Vector<double[][]> results = new Vector<double[][]>();
			
				for (MatrixExponentiator expMethod : tests) {
					if (expMethod.checkMethod()) {
						results.add(expMethod.expm(t));
					}
					else {
						results.add(null);
					}										
				}
				
				for (int i=0;i<results.size()-1;i++) {
					for (int j=(i+1);j<results.size();j++) {
						assertEqualMatrices(results.get(i),results.get(j),tol);
					}
				}								
				System.out.print(".");
				if (dotIter%10==0) {
					System.out.println();
				}
				
			}
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

	private boolean assertEqualMatrixExp(double[][] mat1, double[][] mat2, double tol) {
		assertEquals(mat1.length,mat2.length,0);
		
		for (int i=0;i<mat1.length;i++) {
			for (int j=0;j<mat1.length;j++) {
				A
			}
		}
		
	}

	@Test
	public void testMain() {
		final int numLocations = 3;
		final int numTestTrees=1;
		TestType testType = TestType.TEST_USING_INPUT_TREES;
		

		System.out.println("Running testMain!");
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

		try {
			MCMC mcmc;

			// If a checkpoint file exists, pick up where we left off
			File checkpointFile = new File(config.checkpointFilename);
			if(checkpointFile.exists()) {
				System.out.println("Checking out from: "+config.checkpointFilename);
				mcmc = MCMC.loadFromFile(checkpointFile.getPath());
			}
			// Otherwise, construct a new MCMC
			else {

				// Load data files and prepare data....										
			
				Data data = new DataForTests(config,testType,1,numLocations,numTestTrees);
				System.out.println("Running MCMC...");
				mcmc = new MCMC(); 

				if (config.randomSeed!=null)
					mcmc.setRandomSeed(config.randomSeed);

				// Object that will be asked to create model objects for each chain
				mcmc.setModelFactory(new SeasonalMigrationModelFactory(config,data));

				// Number of chains
				mcmc.setChainCount(config.chainCount);

				// Chains will explore Prior * [Likelihood^tau]
				// where tau == x^heatPower
				// x == 0.0 for the hottest chain (chainId: chainCount - 1);
				// x == 1.0 for the coldest chain (chainId: 0)
				mcmc.setHeatFunction(config.heatPower);

				// Simple default Metropolis-Hastings step for each variable;
				// these will tune to try to reach the target acceptance rate during the burn-in
				// period.
				// 
				// Proposals are minimally intelligent:
				// - normal proposals for distributions with real support (e.g., normal)
				// - multiplier proposals for distributions with positive real support (e.g., gamma)
				// - uniform proposals restricted to min, max for distributions with finite support
				//   (e.g., uniform, beta)
				// - Gibbs sample for binary-valued variables
				// A natural extension would be to automatically choose Gibbs samplers
				// intelligently, but that's not in the current version.
				Step univarStep = new UnivariateProposalStep(config.targetAcceptanceRate, config.burnIn, config.tuneEvery);

				// "Differential evolution MCMC" step, which proposes changes to multiple
				// variables simultaneously, taking into account their correlations in a clever way.
				// Proposes at multiple scales: 8, 16, 32, ... variables at a time.
				// (Since this model only has 5 parameters, it'll just do all 5.)
				// For details of method and parameters, see source for DEMCProposalStep and methods paper
				Step demcStep = new DEMCProposalStep(
						config.targetAcceptanceRate,
						config.burnIn, // Only tune during the burn-in period
						config.tuneEvery,
						config.thin, // Thin historical samples for DEMC in memory
						config.initialHistoryCount, // Collect this many samples before starting DEMC proposals
						8, // Minimum number of variables to propose at a time
						128, // Maximum number of variables to propose at a time
						true, // Use standard "parallel" DEMC proposals (no projection)
						true, // Also use double-size parallel DEMC proposals
						true // Also use snooker proposals
						);

				// Swap steps: even (0/1, 2/3, 4/5) and odd (1/2, 3/4, 5/6);
				// alternating these sets of pairs of chains ensures up to chainCount/2
				// parallelization while swapping, where
				// No tuning, but tuneEvery used to print swap statistics to log file
				// every so often
				Step evenSwapStep = new SwapStep(ChainParity.EVEN, config.tuneEvery);
				Step oddSwapStep = new SwapStep(ChainParity.ODD, config.tuneEvery);

				// Verification step: just asks all models to recalculate
				// log prior, likelihood from scratch and compares to existing value;
				// throws an exception if too much error has accumulated.
				Step verificationStep = new VerificationStep(config.thin, 1E-4);

				// Sample output step
				Step sampOutStep = new SampleOutputStep(config.sampleFilename, config.thin);

				// Marginal-likelihood calculation during run
				Step mlOutStep = new MarginalLikelihoodStep(config.mlFilename, config.burnIn, config.mlthin);

				// Prior-likelihood output step for marginal likelihood calculation
				Step plOutStep = new PriorLikelihoodOutputStep(config.priorLikelihoodFilename, config.thin);

				// Assemble all steps into a sequence; repeat swaps chainCount times
				// since they're so cheap and beneficial for mixing.
				// Each iteration thus includes many little steps:
				// - Proposes changes to all individual parameters
				// - Proposes changes using parallel DEMC, double-size parallel DEMC,
				//   and snooker DEMC at multiple scales
				// - Proposes chainCount * 2 swaps
				// - Every thin iterations, writes samples to a file
				mcmc.addStep(univarStep);
				mcmc.addStep(demcStep);
				for(int i = 0; i < config.chainCount; i++) {
					mcmc.addStep(evenSwapStep);
					mcmc.addStep(oddSwapStep);
				}
				mcmc.addStep(verificationStep);
				mcmc.addStep(sampOutStep);
				mcmc.addStep(mlOutStep);
				mcmc.addStep(plOutStep);
			}

			// Run the thing until checkpointEvery steps at a time;
			// write the checkpoint file in between.
			// The runFor call automatically parallelizes chains.
			while(mcmc.getIterationCount() < config.iterationCount) {
				mcmc.runFor(config.checkpointEvery);
				mcmc.writeToFile(config.checkpointFilename);
			}

			// Tells the MCMC to stop the thread pool so this program will exit
			mcmc.shutdown();
		}
		catch(Throwable e) {
			e.printStackTrace();
			System.exit(1);
		}
	}


}
