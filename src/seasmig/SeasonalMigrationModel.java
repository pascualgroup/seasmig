package seasmig;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

import java.util.*;

import seasmig.Config.Seasonality;

import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.MC3KitException;
import mc3kit.MCMC;

@SuppressWarnings("serial")
public class SeasonalMigrationModel extends GraphicalModel
// TODO: Go over this...
{
	Config config;
	Data data;
	RandomEngine rng;

	DoubleVariable ratePriorRate;
	DoubleDistribution ratePrior;

	DoubleVariable amplitudePriorAlpha;
	DoubleVariable amplitudePriorBeta;
	DoubleDistribution amplitudePrior;

	DoubleVariable twoMatrixPhase;
	DoubleDistribution phasePrior;

	SeasonalMigrationLikelihood likeVar;

	RateParams[][] rateParams;

	public SeasonalMigrationModel(MCMC mcmc, Config config, Data data) throws MC3KitException
	{
		super(mcmc);

		this.config = config;
		this.data = data;	
	} 

	@Override
	protected void buildModel(RandomEngine rng) throws MC3KitException
	{
		// TODO: when Ed changes things, eliminate rng member variable
		// and exclusively use RNG passed into constructor and update methods.		
		this.rng=new MersenneTwister(rng.nextInt()); //TOOD: might be weirdly correlated 

		ratePriorRate = new ExponentialVariable(this, "ratePriorRate", 1.0);
		ratePrior = new ExponentialDistribution(this, "ratePrior", ratePriorRate);

		switch (config.migrationSeasonality) { 
		case SINUSOIDAL :		
			amplitudePriorAlpha = new ExponentialVariable(this, "amplitudePriorAlpha", 1.0);
			amplitudePriorBeta = new ExponentialVariable(this, "amplitudePriorBeta", 1.0);
			amplitudePrior = new BetaDistribution(this, "amplitudePrior", amplitudePriorAlpha, amplitudePriorBeta);

			phasePrior = new UniformDoubleDistribution(this, "phasePrior", 0.0, 1.0);
			break;
		case TWO_CONSTANT_SEASONS:
			twoMatrixPhase = new UniformDoubleVariable(this, "twoMatrixPhase", 0.0, 0.49999);
			break;
		}

		rateParams = new RateParams[config.numLocations][config.numLocations];
		for(int i = 0; i < config.numLocations; i++) {
			for(int j = 0; j < config.numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rateParams[i][j] = new RateParams(i, j);
			}
		}

		likeVar = new SeasonalMigrationLikelihood(this, rng);
	}

	@Override
	public Map<String, Object> makeOutputObject() throws MC3KitException
	{
		Map<String, Object> obj = new LinkedHashMap<String, Object>();

		obj.put("ratePriorRate", ratePriorRate.getValue());

		switch (config.migrationSeasonality) { 
		case SINUSOIDAL :	
			obj.put("amplitudePriorAlpha", amplitudePriorAlpha.getValue());
			obj.put("amplitudePriorBeta", amplitudePriorBeta.getValue());
			break;
		case TWO_CONSTANT_SEASONS:
			obj.put("twoMatrixPhase", twoMatrixPhase.getValue());
			break;
		}

		double[][] rates = new double[config.numLocations][config.numLocations];
		for(int i = 0; i < config.numLocations; i++)
		{
			for(int j = 0; j < config.numLocations; j++)
			{
				if(i == j) continue;
				rates[i][j] = rateParams[i][j].getRate();
			}
		}
		obj.put("rates", rates);

		switch (config.migrationSeasonality) { 
		case SINUSOIDAL :	
			double[][] amplitudes = new double[config.numLocations][config.numLocations];
			double[][] phases = new double[config.numLocations][config.numLocations];

			for(int i = 0; i < config.numLocations; i++)
			{
				for(int j = 0; j < config.numLocations; j++)
				{
					if(i == j) continue;

					amplitudes[i][j] = rateParams[i][j].getAmplitude();
					phases[i][j] = rateParams[i][j].getPhase();
				}
			}
			obj.put("amplitudes", amplitudes);
			obj.put("phases", phases);
			break;
			
		case TWO_CONSTANT_SEASONS: case TWO_CONSTANT_SEASONS_FIXED_PHASE: 
			double[][] rates2 = new double[config.numLocations][config.numLocations];
			for(int i = 0; i < config.numLocations; i++) {
				for(int j = 0; j < config.numLocations; j++) {
					if(i == j) continue;
					rates2[i][j] = rateParams[i][j].getRate2();
				}
			}
			obj.put("rates2", rates2);
			break;
		}

		obj.put("likelihood", this.getLogLikelihood());
		obj.put("prior", this.getLogPrior());
		obj.put("posterior",this.getLogPrior()+this.getLogLikelihood());

		return obj;
	}

	class RateParams
	{
		DoubleVariable rate;
		DoubleVariable rate2;
		DoubleVariable amplitude;
		DoubleVariable phase;

		RateParams(int i, int j) throws MC3KitException
		{
			rate = new DoubleVariable(
					SeasonalMigrationModel.this, "rate_" + i + "_" + j, ratePrior
					);

			switch (config.migrationSeasonality) { 
			case SINUSOIDAL :	
				amplitude = new DoubleVariable(
						SeasonalMigrationModel.this, "amplitude_" + i + "_" + j, amplitudePrior
						);
				phase = new DoubleVariable(
						SeasonalMigrationModel.this, "phase_" + i + "_" + j, phasePrior
						);
				break;
			case TWO_CONSTANT_SEASONS: case TWO_CONSTANT_SEASONS_FIXED_PHASE: 			
				rate2 = new DoubleVariable(
						SeasonalMigrationModel.this, "rate2_" + i + "_" + j, ratePrior
						);
				break;
			}
		}

		public double getRate()
		{
			return rate.getValue();
		}

		public double getRate2()
		{
			return rate2.getValue();
		}

		public double getAmplitude()
		{
			return amplitude.getValue();
		}

		public double getPhase()
		{
			return phase.getValue();
		}
	}

	double getTwoMatrixPhase()
	{
		return twoMatrixPhase.getValue();
	}

	RateParams getRateParams(int i, int j)
	{
		return rateParams[i][j];
	}
}
