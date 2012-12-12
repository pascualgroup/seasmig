package seasmig;

import cern.jet.random.engine.RandomEngine;

import java.util.*;

import seasmig.Config.Seasonality;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.graphical.*;
import mc3kit.graphical.distributions.*;
import mc3kit.graphical.types.*;

@SuppressWarnings("serial")
public class SeasonalMigrationModel extends GraphicalModel
{
	Config config;
	Data data;

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
		ratePriorRate = new ExponentialVariable(this, "ratePriorRate", 1.0);
		ratePrior = new ExponentialDistribution(this, "ratePrior", ratePriorRate);

		switch (config.seasonality) {

		case SINUSOIDAL :
			amplitudePriorAlpha = new ExponentialVariable(this, "amplitudePriorAlpha", 1.0);
			amplitudePriorBeta = new ExponentialVariable(this, "amplitudePriorBeta", 1.0);
			amplitudePrior = new BetaDistribution(this, "amplitudePrior", amplitudePriorAlpha, amplitudePriorBeta);
			phasePrior = new UniformDoubleDistribution(this, "phasePrior", 0.0, 1.0);
			break;
		case TWO_MATRICES:
			twoMatrixPhase = new UniformDoubleVariable(this, "twoMatrixPhase", 0.0, 1.0);
			break;
		}

		rateParams = new RateParams[config.stateCount][config.stateCount];
		for(int i = 0; i < config.stateCount; i++)
		{
			for(int j = 0; j < config.stateCount; j++)
			{
				if(i == j) continue; // rateParams[i,i] remains null
				rateParams[i][j] = new RateParams(i, j);
				// TODO: Ask Ed about prior for TWO_MATRICES...
			}
		}
	}

	@Override
	public Map<String, Object> makeOutputObject() throws MC3KitException
	{
		Map<String, Object> obj = new LinkedHashMap<String, Object>();

		switch (config.seasonality) {
		case  NONE: 
			double[][] rates = new double[config.stateCount][config.stateCount];

			for(int i = 0; i < config.stateCount; i++)
			{
				for(int j = 0; j < config.stateCount; j++)
				{
					if(i == j) continue;
					rates[i][j] = rateParams[i][j].getRate();
				}
			}
			obj.put("ratePriorRate", ratePriorRate.getValue());
			obj.put("rates", rates);
			break;
		case TWO_MATRICES :
			rates = new double[config.stateCount][config.stateCount];
			double[][] rates2 = new double[config.stateCount][config.stateCount];

			for(int i = 0; i < config.stateCount; i++)
			{
				for(int j = 0; j < config.stateCount; j++)
				{
					if(i == j) continue;
					rates[i][j] = rateParams[i][j].getRate();
					rates[i][j] = rateParams[i][j].getRate2();
				}
			}
			obj.put("ratePriorRate", ratePriorRate.getValue());
			obj.put("twoMatrixPhase", twoMatrixPhase.getValue());
			obj.put("rates", rates);
			obj.put("rates2", rates2);
			break;

		case SINUSOIDAL :
			rates = new double[config.stateCount][config.stateCount];
			double[][] amplitudes = new double[config.stateCount][config.stateCount];
			double[][] phases = new double[config.stateCount][config.stateCount];

			for(int i = 0; i < config.stateCount; i++)
			{
				for(int j = 0; j < config.stateCount; j++)
				{
					if(i == j) continue;
					rates[i][j] = rateParams[i][j].getRate();
					amplitudes[i][j] = rateParams[i][j].getAmplitude();
					phases[i][j] = rateParams[i][j].getPhase();
				}
			}
			obj.put("ratePriorRate", ratePriorRate.getValue());
			obj.put("amplitudePriorAlpha", amplitudePriorAlpha.getValue());
			obj.put("amplitudePriorBeta", amplitudePriorBeta.getValue());
			obj.put("rates", rates);
			obj.put("amplitudes", amplitudes);
			obj.put("phases", phases);
		}

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
			switch (config.seasonality) {
			case NONE:
				rate = new DoubleVariable(
						SeasonalMigrationModel.this, "rate_" + i + "_" + j, ratePrior
						);
				break;
			case TWO_MATRICES:
				rate = new DoubleVariable(
						SeasonalMigrationModel.this, "rate_" + i + "_" + j, ratePrior
						);
				rate2 = new DoubleVariable(
						SeasonalMigrationModel.this, "rate2_" + i + "_" + j, ratePrior
						);			
				break;

			case SINUSOIDAL :
				rate = new DoubleVariable(
						SeasonalMigrationModel.this, "rate_" + i + "_" + j, ratePrior
						);
				amplitude = new DoubleVariable(
						SeasonalMigrationModel.this, "amplitude_" + i + "_" + j, amplitudePrior
						);
				phase = new DoubleVariable(
						SeasonalMigrationModel.this, "phase_" + i + "_" + j, phasePrior
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
