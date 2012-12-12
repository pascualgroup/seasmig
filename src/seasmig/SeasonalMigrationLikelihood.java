package seasmig;

import java.util.Set;

import cern.jet.random.engine.RandomEngine;

import mc3kit.MC3KitException;
import mc3kit.graphical.Jack;
import mc3kit.graphical.NoDistribution;
import mc3kit.graphical.RandomVariable;
import mc3kit.graphical.types.DoubleValued;

import treelikelihood.*;

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
		for(int i = 0; i < config.stateCount; i++)
		{
			for(int j = 0; j < config.stateCount-1; j++)
			{
				switch (config.seasonality) {
				case NONE:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);										
					break;									
				case SINUSOIDAL:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate2);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].amplitude);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].phase);
					break;
				case TWO_MATRICES:
					new Jack<DoubleValued>(model, this, model.twoMatrixPhase);
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
		
		RandomEngine rng = model.getRNG();

		switch (config.seasonality) {
		case NONE:
			double[][] rates = new double[config.stateCount][config.stateCount];
			for (int i=0;i<config.stateCount;i++) {
				double rowsum=0;
				for (int j=0;j<config.stateCount;j++) {
					if (i!=j) {
						rates[i][j]=model.getIndexAdjustedRateParams(i, j).getRate();
						rowsum+=rates[i][j];
					}
				}
				rates[i][i]=rowsum;
			}

			MigrationBaseModel likelihoodModel = new ConstantMigrationBaseModel(rates);

			// TODO: Replace with random int from..to			
			data.trees.get(rng.nextInt()%data.trees.size()).copyWithNoCache().logLikelihood(likelihoodModel);

			break;

		case TWO_MATRICES:		
			rates = new double[config.stateCount][config.stateCount];
			double[][] rates2 = new double[config.stateCount][config.stateCount];
			for (int i=0;i<config.stateCount;i++) {
				double row1sum=0;
				double row2sum=0;
				for (int j=0;j<config.stateCount;j++) {		
					if (i!=j) {
						rates[i][j]=model.getIndexAdjustedRateParams(i, j).getRate();
						row1sum+=rates[i][j];
						rates2[i][j]=model.getIndexAdjustedRateParams(i, j).getRate2();
						row2sum+=rates2[i][j];
					}
				}
				rates[i][i]=row1sum;
				rates2[i][i]=row2sum;
			}
			double twoMatrixPhase = model.getTwoMatrixPhase();
			// TODO: Add parameter for first season length
		
			double season1Start=0;			
			if (twoMatrixPhase>0.5) 
				season1Start=(0.5+twoMatrixPhase)%1;
			else 
				season1Start=twoMatrixPhase;			
			double season1End=0.5+season1Start;

			likelihoodModel = new TwoMatrixMigrationBaseModel(rates,rates2,season1Start,season1End);

			// TODO: Replace with random int from..to			
			logLikelihood = data.trees.get(rng.nextInt()%data.trees.size()).copyWithNoCache().logLikelihood(likelihoodModel);

			break;

			//TODO: case SINUSIODIAL:
			// model.getRateParams(i,j).getRate() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getAmplitude() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getPhase() // seasonal phase (between 0 and 1) for sinusoidal model
		}

		setLogP(logLikelihood);
	}

	@Override
	public Object makeOutputObject()
	{
		return null;
	}
}
