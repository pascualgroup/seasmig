package seasmig.migrationmain;

import java.io.File;
import java.io.FileReader;

import seasmig.data.Data;
import seasmig.data.DataFromFiles;
import seasmig.models.migrationandsequencemodels.SequenceAndMigrationModelFactory;
import seasmig.models.migrationmodels.MigrationModelFactory;
import mc3kit.ChainParity;
import mc3kit.MCMC;
import mc3kit.Step;
import mc3kit.VerificationStep;
import mc3kit.monitoring.MarginalLikelihoodStep;
import mc3kit.output.PriorLikelihoodOutputStep;
import mc3kit.output.SampleOutputStep;
import mc3kit.proposal.DEMCProposalStep;
import mc3kit.proposal.UnivariateProposalStep;
import mc3kit.step.swap.IntervalSwapStep;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class SeasonalMigrationMain
{
	public static void main(String[] args) throws Throwable
	{
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
				
				System.out.println("Loading data from files...");
				// Load data files and prepare data....			
				Data data = new DataFromFiles(config);

				System.out.println("Running MCMC");
				mcmc = new MCMC(); 

				if (config.randomSeed!=null)
					mcmc.setRandomSeed(config.randomSeed);

				// Object that will be asked to create model objects for each chain
				switch (config.seqModelType) {
				case NONE:
					mcmc.setModelFactory(new MigrationModelFactory(config,data));
					break;
				case HKY_3CP: case HKY_3CP_AS_INPUT: 
					mcmc.setModelFactory(new SequenceAndMigrationModelFactory(config,data));
					break;
				}
					

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
						32, // Maximum number of variables to propose at a time
						true, // Use standard "parallel" DEMC proposals (no projection)
						true, // Also use double-size parallel DEMC proposals
						true // Also use snooker proposals
						);

				// Swap steps: even (0/1, 2/3, 4/5) and odd (1/2, 3/4, 5/6);
				// alternating these sets of pairs of chains ensures up to chainCount/2
				// parallelization while swapping, where
				// No tuning, but tuneEvery used to print swap statistics to log file
				// every so often
				Step evenSwapStep = new IntervalSwapStep(ChainParity.EVEN, config.tuneEvery, config.swapInterval);
				Step oddSwapStep = new IntervalSwapStep(ChainParity.ODD, config.tuneEvery, config.swapInterval);

				// Verification step: just asks all models to recalculate
				// log prior, likelihood from scratch and compares to existing value;
				// throws an exception if too much error has accumulated.
				Step verificationStep = new VerificationStep(config.thin, config.verificationTolerance);
				
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



