package seasmig;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.ModelFactory;

@SuppressWarnings("serial")
public class SeasonalMigrationFactory implements ModelFactory
{
	Config config;
	Data data;

	public SeasonalMigrationFactory(Config config, Data data)
	{
		this.config = config;
		this.data = data;
	}

	@Override
	public Model createModel(Chain initialChain) throws MC3KitException {
		switch (config.migrationSeasonality) {
		case NONE:			
			return new SeasonalMigrationModelNoSeasonality(initialChain, config, data);			

//		case TWO_CONSTANT_SEASONS:			
//			return new SeasonalMigrationModelTwoConstantSeasons(initialChain, config, trees);			
//
//		case TWO_CONSTANT_SEASONS_FIXED_PHASE:			
//			return new SeasonalMigrationModelTwoConstantSeasonsFixedPhase(initialChain, config, trees);			
//
//		case SINUSOIDAL:			
//			return new SeasonalMigrationModelTwoConstantSeasonsSinusoidal(initialChain, config, trees);			
		}
		return null;
	}
}
