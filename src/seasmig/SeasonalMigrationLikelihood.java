package seasmig;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

import jebl.math.Random;

import org.apache.log4j.Priority;

import cern.jet.random.engine.RandomEngine;

import mc3kit.FormattingLogger;
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
		for(int i = 0; i < config.locationCount; i++)
		{
			for(int j = 0; j < config.locationCount; j++)
			{
				switch (config.seasonality) {
				case NONE:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);										
					break;									
				case SINUSOIDAL:
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].amplitude);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].phase);
					break;
				case TWO_CONSTANT_SEASONS:
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

		//RandomEngine rng = model.getRNG();

		switch (config.seasonality) {
		case NONE:
			double[][] rates = new double[config.locationCount][config.locationCount];
			for (int i=0;i<config.locationCount;i++) {
				double rowsum=0;
				for (int j=0;j<config.locationCount;j++) {
					if (i!=j) {
						rates[i][j]=model.getRateParams(i, j).getRate();
						rowsum-=rates[i][j];
					}
				}
				rates[i][i]=rowsum;
			}

			MigrationBaseModel likelihoodModel = new ConstantMigrationBaseModel(rates);

			data.trees.get(model.rng.nextInt()%data.trees.size()).copyWithNoCache().logLikelihood(likelihoodModel);

			break;

		case TWO_CONSTANT_SEASONS:		
			rates = new double[config.locationCount][config.locationCount];
			double[][] rates2 = new double[config.locationCount][config.locationCount];
			for (int i=0;i<config.locationCount;i++) {
				double row1sum=0;
				double row2sum=0;
				for (int j=0;j<config.locationCount;j++) {		
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

			likelihoodModel = new TwoSeasonMigrationBaseModel(rates,rates2,season1Start,season1End);

			logLikelihood = data.trees.get(model.rng.nextInt()%data.trees.size()).copyWithNoCache().logLikelihood(likelihoodModel);

			break;

		case SINUSOIDAL: 
			// model.getRateParams(i,j).getRate() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getAmplitude() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getPhase() // seasonal phase (between 0 and 1) for sinusoidal model
			// ??? rate*(1+amp*sin(2*pi*t+2*pi*phase)) ???

			rates = new double[config.locationCount][config.locationCount];
			double[][] amp = new double[config.locationCount][config.locationCount];
			double[][] phase = new double[config.locationCount][config.locationCount];
			for (int i=0;i<config.locationCount;i++) {
				for (int j=0;j<config.locationCount;j++) {		
					if (i!=j) {
						rates[i][j]=model.getRateParams(i, j).getRate();
						amp[i][j]=model.getRateParams(i, j).getAmplitude();
						amp[i][j]=model.getRateParams(i, j).getPhase();
					}
				}
			}
		
			likelihoodModel = new SinusoidialSeasonalMigrationBaseModel(rates, amp, phase);

			logLikelihood = data.trees.get(model.rng.nextInt()%data.trees.size()).copyWithNoCache().logLikelihood(likelihoodModel);
		}
		
		
		//TODO: Figure out zero log likelihood in files
		setLogP(logLikelihood);
		
	}

	@Override
	public Object makeOutputObject() {
//		// TODO: Ask Ed if this makes sense (probably doesn't :)
//		
//		Map<String, Object> obj = new LinkedHashMap<String, Object>();
//
//		obj.put("model", model);
//		
//		return obj;
		return null;
	}


}
