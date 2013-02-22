package seasmig;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

import java.util.*;

import seasmig.Config.Seasonality;
import treelikelihood.LikelihoodTree;

import mc3kit.Chain;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.Model;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformIntDistribution;

public class SeasonalMigrationModelNoSeasonality extends Model
// TODO: Go over this...
{
	Config config;
	Collection<LikelihoodTree> trees;
	RandomEngine rng;

	DoubleVariable ratePriorRate;
	DoubleDistribution ratePrior;

	SeasonalMigrationLikelihood likeVar;

	public SeasonalMigrationModelNoSeasonality(Chain initialChain, Config config, Collection<LikelihoodTree> trees) throws MC3KitException
	{
		this.config = config;
		this.trees = trees;	
		
		 Model m = new Model(initialChain);
		 rng = m.getRng();
	        
	     m.beginConstruction();
	     new IntVariable(m, "v", new UniformIntDistribution(m, 1,4));
	   
	     new DoubleVariable(m, "ratePriorRate", new ExponentialDistribution(m, 1.0));
	     new DoubleVariable(m, "ratePrior", new ExponentialDistribution(m,"ratePriorRate"));
		
		 for(int i = 0; i < config.numLocations; i++) {
			for(int j = 0; j < config.numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				new DoubleVariable(m, "rateParams"+Integer.toString(i)+Integer.toString(j), new ExponentialDistribution(m,"ratePrior"));
			}
		 }

		likeVar = new SeasonalMigrationLikelihood(m,rng);
		
		m.endConstruction();

	} 

}
