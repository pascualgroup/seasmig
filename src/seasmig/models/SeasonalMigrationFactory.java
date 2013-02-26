package seasmig.models;

import seasmig.Config;
import seasmig.Config.TwoConstantSeasonsParameterization;
import seasmig.Data;
import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.ModelFactory;

@SuppressWarnings("serial")
public class SeasonalMigrationFactory implements ModelFactory
{
	Config config;
	Data data;

	protected SeasonalMigrationFactory() {};

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

		case TWO_CONSTANT_SEASONS: 
			boolean fixPhase = false;
			switch (config.twoSeasonPhase) {
			case FIXED:
				fixPhase=true;
				break;
			case FREE:
				fixPhase=false;
				break;
			default: 
				System.err.println(config.twoSeasonPhase+" two season phase not implemented!");
			}
			boolean fixTo = config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_TO_DIFF;
			boolean fixFrom = config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_DIFF;

			if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.ORIGNIAL_PARAMETERIZATION_ALL_FREE) {
				return new SeasonalMigrationModelTwoConstantSeasonsOrigParametarization(initialChain, config, data,fixPhase);
			}
			else if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.DIFF_PARAMETERIZATION || config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_DIFF || config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_TO_DIFF) {
				return new SeasonalMigrationModelTwoConstantSeasons(initialChain, config, data,fixPhase, fixFrom, fixTo);
			}
			else {
				System.err.println("twiSeasonParameterization"+config.twoSeasonParameterization+" not implemented!");
				System.exit(1);
			}
			break;
			default:
				
			//
			//		case SINUSOIDAL:			
			//			return new SeasonalMigrationModelTwoConstantSeasonsSinusoidal(initialChain, config, trees);			

		}
		return null;
	}
}
