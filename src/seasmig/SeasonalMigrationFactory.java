package seasmig;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.Model;
import mc3kit.ModelFactory;

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
		// TODO Auto-generated method stub
		return new SeasonalMigrationModel(mcmc, config, data);
	}
}
