package seasmig.models.migrationandsequencemodels;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.ModelFactory;
import seasmig.data.Data;
import seasmig.migrationmain.Config;
import seasmig.migrationmain.Config.MigrationModelType;
import seasmig.migrationmain.Config.SeqModelType;

@SuppressWarnings("serial")
public class SequenceAndMigrationModelFactory implements ModelFactory
{
	Config config;
	Data data;

	protected SequenceAndMigrationModelFactory() {};

	public SequenceAndMigrationModelFactory(Config config, Data data)
	{
		this.config = config;
		this.data = data;
	}

	@Override
	public Model createModel(Chain initialChain) throws MC3KitException {
		switch (config.migrationModelType) {
		case CONSTANT: case CONSTANT_AS_INPUT:	
			boolean inputMigrationModel=(config.migrationModelType==MigrationModelType.CONSTANT_AS_INPUT);
			switch (config.seqModelType) {
			case HKY_3CP: case HKY_3CP_AS_INPUT:
				boolean inputSeqModel=(config.seqModelType==SeqModelType.HKY_3CP_AS_INPUT);
				return new HKY_3CP_NoMigrationSeasonality(initialChain, config, data,inputMigrationModel, inputSeqModel);			
			case NONE:
				System.err.println("config error... seqModelType==NONE and ???");
				System.exit(-1);
			}
				
		default:
			System.err.println("only constant migration model implemented with sequence stochastic mapping...");
			System.exit(-1);
		}
		return null;
	}
}
