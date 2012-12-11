package seasmig;

import java.util.Set;

import seasmig.Config.Seasonality;

import mc3kit.MC3KitException;
import mc3kit.graphical.GraphicalModel;
import mc3kit.graphical.Jack;
import mc3kit.graphical.NoDistribution;
import mc3kit.graphical.RandomVariable;
import mc3kit.graphical.types.DoubleValued;

public class SeasonalMigrationLikelihood extends RandomVariable<NoDistribution>
{
	SeasonalMigrationModel model;
	Config config;
	Data data;

	public SeasonalMigrationLikelihood(SeasonalMigrationModel model) throws MC3KitException
	{
		super(model, "likelihood");
		this.model = model;
		this.config = model.config;
		this.data = model.data;

		// This is the weird way dependencies are declared. It will become unweird someday.

		if(config.seasonality == Seasonality.TWO_MATRICES)
		{
			new Jack<DoubleValued>(model, this, model.twoMatrixPhase);
		}

		for(int i = 0; i < config.stateCount; i++)
		{
			for(int j = 0; j < config.stateCount; j++)
			{
				switch (config.seasonality) {
				case NONE:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);										
					break;
					// TODO: Add SINUSOIDIAL SEASONALITY
					//				case SINUSOIDAL:
					//					new Jack<DoubleValued>(model, this, model.rateParams[i][j].amplitude);
					//					new Jack<DoubleValued>(model, this, model.rateParams[i][j].phase);
					//					break;
				case TWO_MATRICES:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate2);
					break;
				default:
				}
			}
		}
	}

	@Override
	public boolean update(Set<Jack<?>> changedJacks) throws MC3KitException
	{
		recalculate();
		return true;
	}

	void recalculate() throws MC3KitException
	{
		double logLikelihood = 0.0;
		
		switch (config.seasonality) {
		case NONE:
			double[][] rates = new double[config.stateCount][config.stateCount];
			for (int i=0;i<config.stateCount;i++) {
				for (int j=0;j<config.stateCount;j++) {
					rates[i][j]=model.getRateParams(i, j).getRate();
				}
			}
			
			MigrationModel 
			
			break;
			
		case TWO_MATRICES:		
			rates = new double[config.stateCount][config.stateCount];
			double[][] rates2 = new double[config.stateCount][config.stateCount];
			for (int i=0;i<config.stateCount;i++) {
				for (int j=0;j<config.stateCount;j++) {
					rates[i][j]=model.getRateParams(i, j).getRate();
					rates2[i][j]=model.getRateParams(i, j).getRate2();
				}
			}
			double twoMatrixPhase = model.getTwoMatrixPhase();
			double twoMatrixLength = 0.5; // TODO: Add parameter for first season length
			
			break;
			
		//TODO: case SINUSIODIAL:
			
		}

		// TODO: calculate log-likelihood here.

		// Get parameters this way (all return double):
		// model.getRateParams(i,j).getRate() // rate #1
		// model.getRateParams(i,j).getRate2() // rate #2 for two-matrix model
		// model.getRateParams(i,j).getAmplitude() // seasonal amplitude for sinusoidal model
		// model.getRateParams(i,j).getPhase() // seasonal phase (between 0 and 1) for sinusoidal model
		// model.getTwoMatrixPhase() // change point in year for two matrices (between 0 and 1)

		setLogP(logLikelihood);
	}

	@Override
	public Object makeOutputObject()
	{
		return null;
	}
}
