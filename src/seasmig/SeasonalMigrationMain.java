package seasmig;

import static org.junit.Assert.assertEquals;

import java.io.FileReader;
import java.io.IOException;

import jebl.evolution.io.ImportException;

import mc3kit.Chain;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.Model;
import mc3kit.ModelFactory;
import mc3kit.distributions.UniformIntDistribution;
import mc3kit.proposal.UnivariateProposalStep;


import treelikelihood.*;

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

		// Load data files and prepare data....			
		Data data = new Data(config);
		
		// Setup MCMC
		System.out.print("Setting up MCMC....");
		MCMC mcmc = new MCMC();
		mcmc.setRandomSeed(config.randomSeed); // TODO: check that this changes...

		MCMC.setLogLevel(config.logLevel);

		ModelFactory mf = new SeasonalMigrationFactory(config, data) {
			@Override
			public Model createModel(Chain initialChain) throws MC3KitException {
				Model m = new Model(initialChain);
				
				m.beginConstruction();
				new IntVariable(m, "v", new UniformIntDistribution(m, 1,4));
				m.endConstruction();

				return m;
			}
		};

		mcmc.setModelFactory(mf);

		UnivariateProposalStep proposalStep = new UnivariateProposalStep(0.25, 100, config.burnIn);
		mcmc.addStep(proposalStep);

		System.out.println(" done!");
		// Run, collect statistics, and check moments against expected distribution
		
		System.out.println("Running MCMC...");
		
		System.out.println("state\tlikelihood\thours/million states");
		double sum = 0;
		for(long i = 0; i < config.iterCount; i++) {
			mcmc.step();
			mcmc.getModel().recalculate();

			assertEquals(i + 1, mcmc.getIterationCount());

			if(i >= config.burnIn) {
				int val = mcmc.getModel().getIntVariable("v").getValue();
				sum += val;
			}
		}

		double N = config.iterCount - config.burnIn;
		System.out.println(sum/N);		
		System.out.println("done!");	
	}
	
}
