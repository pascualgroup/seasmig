package seasmig.models.migrationandsequencemodels;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.ModelFactory;
import seasmig.data.Data;
import seasmig.migrationmain.Config;

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
		case CONSTANT:	
			switch (config.seqModelType) {
			case HKY_3CP:
				return new HKY3CPConstSeasonalMigrationModelNoSeasonality(initialChain, config, data);
			}
				
		default:
			System.err.println("only constant migration model implemented with seasonality...");
		}
		return null;
	}
}
