package seasmig;

import java.util.Collection;

import treelikelihood.LikelihoodTree;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.Model;
import mc3kit.ModelFactory;

public class SeasonalMigrationFactory implements ModelFactory
{
	Config config;
	Collection<LikelihoodTree> trees;

	public SeasonalMigrationFactory(Config config, Collection<LikelihoodTree> trees)
	{
		this.config = config;
		this.trees = trees;
	}

	@Override
	public Model createModel(Chain initialChain) throws MC3KitException {
		switch (config.migrationSeasonality) {
		case NONE:			
			return new SeasonalMigrationModelNoSeasonality(initialChain, config, trees);			

		case TWO_CONSTANT_SEASONS:			
			return new SeasonalMigrationModelTwoConstantSeasons(initialChain, config, trees);			

		case TWO_CONSTANT_SEASONS_FIXED_PHASE:			
			return new SeasonalMigrationModelTwoConstantSeasonsFixedPhase(initialChain, config, trees);			

		case SINUSOIDAL:			
			return new SeasonalMigrationModelTwoConstantSeasonsSinusoidal(initialChain, config, trees);			
		}
	}
}
