package seasmig;

import java.util.Set;

import cern.jet.random.Uniform;
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

	double oldLogLikelihood;
	long stateCount=1;
	private long time;
	private long oldTime;

	public SeasonalMigrationLikelihood(SeasonalMigrationModel model, RandomEngine rng) throws MC3KitException
	{
		super(model, "likelihood");

		oldTime=System.currentTimeMillis();

		setObserved(true);

		this.model = model;
		this.config = model.config;		
		this.data = model.data;


		switch (config.migrationSeasonality) {

		case TWO_CONSTANT_SEASONS:
			// This is the weird way dependencies are declared. It will become unweird someday.
			new Jack<DoubleValued>(model, this, model.twoMatrixPhase);
			for(int i = 0; i < config.numLocations; i++) 		{
				for(int j = 0; j < config.numLocations; j++) 			{
					if(i == j) continue;
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate2);					
				}
			}
			break;
		case TWO_CONSTANT_SEASONS_FIXED_PHASE:
			// This is the weird way dependencies are declared. It will become unweird someday.
			for(int i = 0; i < config.numLocations; i++) 		{
				for(int j = 0; j < config.numLocations; j++) 			{
					if(i == j) continue;
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate2);					
				}
			}
			break;
		case NONE:
			for(int i = 0; i < config.numLocations; i++) 		{
				for(int j = 0; j < config.numLocations; j++) 			{
					if(i == j) continue;					
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);					
				}
			}
			break;
		case SINUSOIDAL:

			for(int i = 0; i < config.numLocations; i++) 		{
				for(int j = 0; j < config.numLocations; j++) 			{
					if(i == j) continue;
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].rate);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].amplitude);
					new Jack<DoubleValued>(model, this, model.rateParams[i][j].phase);			
				}
				break;
			}
		}
	}

	@Override
	public boolean update(Set<Jack<?>> changedJacks) throws MC3KitException 	{
		recalculate();
		return true;
	}

	void recalculate() throws MC3KitException 	{
		oldLogLikelihood = getLogP();

		double logLikelihood = 0.0;

		MigrationBaseModel migrationBaseModel = null;

		switch (config.migrationSeasonality) {
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
		case TWO_CONSTANT_SEASONS_FIXED_PHASE:		
			rates = new double[config.numLocations][config.numLocations];
			rates2 = new double[config.numLocations][config.numLocations];
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
			twoMatrixPhase = config.fixedPhase;
			// TODO: Add parameter for first season length

			season1Start=0;			
			season1Start=twoMatrixPhase;
			season1End=0.5+season1Start;

			migrationBaseModel = new TwoSeasonMigrationBaseModel(rates,rates2,season1Start,season1End);	
			break;

		case SINUSOIDAL: 
			// model.getRateParams(i,j).getRate() // overall rate for sinusoidal model
			// model.getRateParams(i,j).getAmplitude() // seasonal amplitude for sinusoidal model
			// model.getRateParams(i,j).getPhase() // seasonal phase (between 0 and 1) for sinusoidal model
			// rate*(1+amp*sin(2*pi*t+2*pi*phase))

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

		Uniform uniform = new Uniform(model.rng);
		LikelihoodTree workingCopy = data.trees.get(uniform.nextIntFromTo(0,data.trees.size()-1)).copy(); 
		workingCopy.setLikelihoodModel(migrationBaseModel);
		logLikelihood=workingCopy.logLikelihood();

		setLogP(logLikelihood);

		// TODO: organize this...
		// Display state per hour....
		stateCount+=1;
		if (stateCount==25) {
			time = System.currentTimeMillis();
			System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/100L*1000000L/(60L*60L*1000L));
		}
		if (stateCount==250) {
			time = System.currentTimeMillis();
			System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/1000L*1000000L/(60L*60L*1000L));
		}
		if (stateCount%(config.printEveryNStates/config.chainCount)==0) {
			time = System.currentTimeMillis();
			System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/10000L*1000000L/(60L*60L*1000L));
			oldTime=time;
		}
		//

	}



	@Override
	public boolean updateAfterRejection(Set<Jack<?>> changedParents)
			throws MC3KitException
			{
		setLogP(oldLogLikelihood);
		return true;
			}

	@Override
	public Object makeOutputObject() {
		// TODO: Ask Ed what goes here :)
		return null;
	}
}
