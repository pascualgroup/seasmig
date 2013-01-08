package seasmig;

import java.io.FileReader;

import mc3kit.FormattingLogger;
import mc3kit.MCMC;
import mc3kit.PowerHeatFunction;
import mc3kit.graphical.step.BlockDifferentialEvolutionStep;
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
import treelikelihood.LikelihoodTree;
import treelikelihood.TreeWithLocations;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

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
			
			// Load data files and prepare data....			
			Data data = new Data(config);
			
			// Tests...			
			if (config.runMode==RunMode.TEST) {
				System.out.print("Running likelihood test...\n");
				System.out.println("Calculating tree likelihood using the same model used to create the tree: SEASONALITY "+config.migrationSeasonality);
				System.out.println(data.createModel.print());
				double createLikelihood = 0;
				for (LikelihoodTree tree : data.trees) {
					System.out.print(".");
					// TODO: maybe get likelihood to not require copy...
					LikelihoodTree workingCopy = tree.workingCopy();
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
						LikelihoodTree workingCopy = tree.workingCopy();
						workingCopy.setLikelihoodModel(data.testModels.get(i));
						testLikelihood+=workingCopy.logLikelihood();
					}
					testLikelihood=testLikelihood/data.trees.size();
					System.out.println(testLikelihood);
				}
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
			BlockDifferentialEvolutionStep deStep = new BlockDifferentialEvolutionStep();
			deStep.setStatsFilename(config.demcStatsFilename);
			deStep.setRecordHistoryEvery(config.thin);
			deStep.setRecordHistoryAfter(config.recordHistoryAfter);
			deStep.setTuneEvery(config.tuneEvery);
			deStep.setTuneFor(config.tuneFor);
			deStep.setInitialHistoryCount(config.initialHistoryCount);
			deStep.setRecordHeatedStats(false);
			deStep.setRecordStatsAfterTuning(true);
			
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
			mcmc.addStep(deStep);
			for(int i = 0; i < config.chainCount; i++)	{
				mcmc.addStep(evenSwapStep);
				mcmc.addStep(oddSwapStep);				
			}
			mcmc.addStep(sampOutStep);
			mcmc.addStep(plOutStep);
			
			System.out.println(" done");
			System.out.print("Running MCMC...");
			mcmc.run();			
			System.out.println("done");
		}
		catch(Throwable e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
}
