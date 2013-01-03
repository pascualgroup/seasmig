package seasmig;

import java.util.Set;

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
		
		// This is the weird way dependencies are declared. It will become unweird someday.
		
		if(config.seasonality == Config.Seasonality.TWO_CONSTANT_SEASONS)
		{
			new Jack<DoubleValued>(model, this, model.twoMatrixPhase);
		}
		
		for(int i = 0; i < config.numLocations; i++)
		{
			for(int j = 0; j < config.numLocations; j++)
			{
				new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
				if(config.seasonality == Config.Seasonality.SINUSOIDAL)
				{
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].amplitude);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].phase);
				}
				else if(config.seasonality ==Config.Seasonality.TWO_CONSTANT_SEASONS)
				{
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate2);
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

		MigrationBaseModel migrationBaseModel = null;

		switch (config.seasonality) {
		case NONE:
			double[][] rates = new double[config.numLocations][config.numLocations];
			for (int i=0;i<config.numLocations;i++) {
				double rowsum=0;
				for (int j=0;j<config.numLocations;j++) {
					if (i!=j) {
						rates[i][j]=model.getRateParams(i, j).getRate();
						rowsum-=rates[i][j];
					}
				}
				rates[i][i]=rowsum;
			}

			migrationBaseModel = new ConstantMigrationBaseModel(rates);
			break;

		case TWO_CONSTANT_SEASONS:		
			rates = new double[config.numLocations][config.numLocations];
			double[][] rates2 = new double[config.numLocations][config.numLocations];
			for (int i=0;i<config.numLocations;i++) {
				double row1sum=0;
				double row2sum=0;
				for (int j=0;j<config.numLocations;j++) {		
					if (i!=j) {
						rates[i][j]=model.getRateParams(i, j).getRate();
						row1sum-=rates[i][j];
						rates2[i][j]=model.getRateParams(i, j).getRate2();
						row2sum-=rates2[i][j];
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

			migrationBaseModel = new TwoSeasonMigrationBaseModel(rates,rates2,season1Start,season1End);	
			break;

		case SINUSOIDAL: 
			// model.getRateParams(i,j).getRate() // overall rate for sinusoidal model
			// model.getRateParams(i,j).getAmplitude() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getPhase() // seasonal phase (between 0 and 1) for sinusoidal model
			// ??? rate*(1+amp*sin(2*pi*t+2*pi*phase)) ???

			rates = new double[config.numLocations][config.numLocations];
			double[][] amp = new double[config.numLocations][config.numLocations];
			double[][] phase = new double[config.numLocations][config.numLocations];
			for (int i=0;i<config.numLocations;i++) {
				for (int j=0;j<config.numLocations;j++) {		
					if (i!=j) {
						rates[i][j]=model.getRateParams(i, j).getRate();
						amp[i][j]=model.getRateParams(i, j).getAmplitude();
						amp[i][j]=model.getRateParams(i, j).getPhase();
					}
				}
			}
		
			migrationBaseModel = new SinusoidialSeasonalMigrationBaseModel(rates, amp, phase);
		}
		
		// TODO: maybe get likelihood to work without copy...
		LikelihoodTree workingCopy = data.trees.get(model.rng.nextInt()%data.trees.size()).workingCopy(); 
		workingCopy.setLikelihoodModel(migrationBaseModel);
		logLikelihood=workingCopy.logLikelihood();
		
		//TODO: Figure out zero log likelihood in files
		setLogP(logLikelihood);
	}

	@Override
	public Object makeOutputObject() {
		// TODO: Ask Ed what goes here :)
		return null;
	}
}
