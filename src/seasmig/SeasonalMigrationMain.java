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

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class SeasonalMigrationMain
{
	public static void main(String[] args)
	{
		try
		{
			// Load config
			Gson gson = new Gson();
			Config config = gson.fromJson(new FileReader("config.json"), Config.class);
			
			// Load data files
			Data data = new Data(config);
				
			MCMC mcmc = new MCMC();
			mcmc.setRandomSeed(config.randomSeed);
			mcmc.setChainCount(config.chainCount);
			
			// Initialize debugging logger
			Logger logger = Logger.getRootLogger();
			logger.setLevel(config.logLevel.getLog4jLevel());
			Layout layout = new TTCCLayout("ISO8601");
			if(config.logFilename.equals("-"))
			{
				logger.addAppender(new ConsoleAppender(layout));
			}
			else
			{
				logger.addAppender(new FileAppender(layout, config.logFilename));
			}
			FormattingLogger fmtLogger = new FormattingLogger(logger);
			
			String configJson = new GsonBuilder().setPrettyPrinting().create().toJson(config, Config.class);
			fmtLogger.info("Parsed config: %s", configJson);
			
			mcmc.setLogger(fmtLogger);
			
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
			for(int i = 0; i < config.chainCount; i++)
			{
				mcmc.addStep(evenSwapStep);
				mcmc.addStep(oddSwapStep);
			}
			mcmc.addStep(sampOutStep);
			mcmc.addStep(plOutStep);
			
			mcmc.run();
		}
		catch(Throwable e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
}
