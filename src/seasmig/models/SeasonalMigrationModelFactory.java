package seasmig.models;

import seasmig.Config;
import seasmig.Config.TwoConstantSeasonsParameterization;
import seasmig.Data;
import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.ModelFactory;

@SuppressWarnings("serial")
public class SeasonalMigrationModelFactory implements ModelFactory
{
	Config config;
	Data data;

	protected SeasonalMigrationModelFactory() {};

	public SeasonalMigrationModelFactory(Config config, Data data)
	{
		this.config = config;
		this.data = data;
	}

	@Override
	public Model createModel(Chain initialChain) throws MC3KitException {
		switch (config.migrationSeasonality) {
		case NONE:	
			switch (config.noSeasonalityParameterization) {
			case ALL:
				return new SeasonalMigrationModelNoSeasonality(initialChain, config, data);
			case VARIABLE_SELECTION:
				return new SeasonalMigrationModelNoSeasonalityVarSelect(initialChain, config, data);
			}
		case N_CONSTANT_SEASONS: 	
			return new SeasonalMigrationModelNConstantSeasons(initialChain,config,data,config.nSeasonalParts);
		case N_CONSTANT_SEASONS_VAR_SELECT: 	
			return new SeasonalMigrationModelNConstantSeasonsVarSelect(initialChain,config,data,config.nSeasonalParts);
		case TWO_CONSTANT_SEASONS: 
			boolean fixPhase = false;
			boolean fixLength = false;
			switch (config.twoSeasonPhase) {
			case FIXED_PHASE_FIXED_LENGTH:
				fixPhase=true;
				fixLength=true;
				break;
			case FREE_PHASE_FREE_LENGTH:
				fixPhase=false;
				fixLength=false;
				break;
			case FIXED_PHASE_FREE_LENGTH:
				fixPhase=true;
				fixLength=false;
				break;
			case FREE_PHASE_FIXED_LENGTH:
				fixPhase=false;
				fixLength=true;
				break;
			default: 
				System.err.println(config.twoSeasonPhase+" two season phase not implemented!");
			}
			if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.VARIABLE_SELECTION)
				return new SeasonalMigrationModelTwoConstantSeasonsFullVariableSelection(initialChain,config,data,fixPhase,fixLength,config.fixRate);
			
			if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.VARIABLE_SELECTION_DIFF)
				return new SeasonalMigrationModelTwoConstantSeasonsVariableSelection(initialChain,config,data,fixPhase,fixLength,config.fixRate);
			
			if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.VARIABLE_SELECTION_GTR)
				return new SeasonalMigrationModelTwoConstantSeasonsVariableSelectionGTR(initialChain,config,data,fixPhase,fixLength);
			
			boolean fixTo = config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_TO_DIFF || config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_TO_DIFF;
			boolean fixFrom = config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_DIFF || config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_TO_DIFF;
		
			if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.RATES12_PARAMETERZIATION) {
				return new SeasonalMigrationModelTwoConstantSeasonsOrigParametarization(initialChain, config, data,fixPhase);
			}
			else if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.RATES12_VARIABLE_SELECTION) {
				return new SeasonalMigrationModelTwoConstantSeasonsOrigParametarizationVarSelection(initialChain, config, data,fixPhase);
			}
			else if (config.twoSeasonParameterization==TwoConstantSeasonsParameterization.DIFF_PARAMETERIZATION ||
					 config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_DIFF || 
					 config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_TO_DIFF ||
					 config.twoSeasonParameterization==TwoConstantSeasonsParameterization.FIX_FROM_TO_DIFF) {
				return new SeasonalMigrationModelTwoConstantSeasons(initialChain, config, data,fixPhase, fixFrom, fixTo);
			}
			else {
				System.err.println("twoSeasonParameterization: "+config.twoSeasonParameterization+" not implemented!");
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
