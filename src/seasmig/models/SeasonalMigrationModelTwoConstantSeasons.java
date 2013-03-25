package seasmig.models;

import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;
import mc3kit.*;
import mc3kit.distributions.*;

@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasons extends Model {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	DoubleVariable[][] diffMultipliers;
	DoubleVariable[] diffMultipliersRowCol;
	DoubleVariable seasonalPhase;
	double seasonalPhaseRealization;
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;	

	boolean fixedPhase;
	boolean fixTo;
	boolean fixFrom;
	private ExponentialDistribution ratePriorDist;

	protected SeasonalMigrationModelTwoConstantSeasons() { }

	public SeasonalMigrationModelTwoConstantSeasons(Chain initialChain, Config config, Data data, boolean fixedPhase, boolean fixFrom, boolean fixTo) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.config = config;
		this.data = data;
		this.fixedPhase=fixedPhase;
		this.fixTo=fixTo;
		this.fixFrom=fixFrom;
		numLocations=data.getNumLocations();
		nTrees=data.getTrees().size();		
		rates = new DoubleVariable[numLocations][numLocations];

		diffMultipliersRowCol = new DoubleVariable[numLocations];
		diffMultipliers = new DoubleVariable[numLocations][numLocations];

		beginConstruction();

		if (nTrees>1) {
			treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));
		}
		
		if (!fixedPhase) {
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));
		}
		else {
			seasonalPhaseRealization=config.fixedPhase;
		}

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);
		DoubleDistribution diffMultiplierPriorDist = new UniformDistribution(this,-1.0,1.0);

		for (int i=0; i< numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "meanRates."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				if (!fixFrom && !fixTo) 
					diffMultipliers[i][j] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i)+"."+Integer.toString(j), diffMultiplierPriorDist);
			}		
		}

		if (fixFrom || fixTo) {
			for(int i = 0; i < numLocations; i++) {
				diffMultipliersRowCol[i] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i), diffMultiplierPriorDist);
			}		
		}
	
		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogLikelihood;
		private double logMaxLikelihood = Double.NEGATIVE_INFINITY;

		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasons m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			if (nTrees>1)
				m.addEdge(this, m.treeIndex);
			if (!fixedPhase)
				m.addEdge(this, m.seasonalPhase);

			for(int i = 0; i < numLocations; i++) {
				if (diffMultipliersRowCol[i]!=null ) 
					m.addEdge(this, diffMultipliersRowCol[i]);				
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;
					m.addEdge(this,rates[i][j]);
					if (diffMultipliers[i][j]!=null) 
						m.addEdge(this,diffMultipliers[i][j]);		
				}
			}
		}


		/*
		 * The update method is called whenever a node this variable
		 * depends on has changed; random variables should typically
		 * just recalculate their log-probability and call setLogP here.
		 * 
		 * The method returns a boolean indicating whether this node actually
		 * changed during the update process. Typically you just return true,
		 * but this can be used to optimize updating: e.g., a maximum
		 * function may get updated very frequently without its value
		 * actually changing.
		 * 
		 * This example assumes that the data is a mixture between two normal
		 * distributions.
		 */
		@Override
		public boolean update() {

			double logLikelihood = 0.0;
			//double logPrior = 0.0;

			double[][] rates1doubleForm = new double[numLocations][numLocations];
			double[][] rates2doubleForm = new double[numLocations][numLocations];
			for (int i=0;i<numLocations;i++) {
				double rowsum1=0;
				double rowsum2=0;
				for (int j=0;j<numLocations;j++) {
					if (i!=j) {
						if (diffMultipliers[i][j]!=null) {
							rates1doubleForm[i][j]=rates[i][j].getValue()*(1-diffMultipliers[i][j].getValue());
							rates2doubleForm[i][j]=rates[i][j].getValue()*(1+diffMultipliers[i][j].getValue());
						}
						else if (diffMultipliersRowCol[j]!=null && (fixTo)) {
							rates1doubleForm[i][j]=rates[i][j].getValue()*(1-diffMultipliersRowCol[j].getValue());
							rates2doubleForm[i][j]=rates[i][j].getValue()*(1+diffMultipliersRowCol[j].getValue());
						}
						else if (diffMultipliersRowCol[i]!=null && (fixFrom)) {
							rates1doubleForm[i][j]=rates[i][j].getValue()*(1-diffMultipliersRowCol[i].getValue());
							rates2doubleForm[i][j]=rates[i][j].getValue()*(1+diffMultipliersRowCol[i].getValue());
						}
						else {
							rates1doubleForm[i][j]=rates[i][j].getValue();
							rates2doubleForm[i][j]=rates[i][j].getValue();{
							}
						}
						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;		
			}

			// TODO: add update to migration model instead of reconstructing...
			if (!fixedPhase)
				seasonalPhaseRealization=seasonalPhase.getValue();
			MigrationBaseModel migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhaseRealization,seasonalPhaseRealization+0.5);
			LikelihoodTree workingCopy;
			if (nTrees>1)
				workingCopy = data.getTrees().get((int)treeIndex.getValue()).copy(); 
			else
				workingCopy = data.getTrees().get(0).copy();
			
			workingCopy.setLikelihoodModel(migrationBaseModel);
			logLikelihood=workingCopy.logLikelihood();								

			setLogP(logLikelihood);			
			oldLogLikelihood=logLikelihood;
			if (logLikelihood>logMaxLikelihood) {
				logMaxLikelihood=logLikelihood;
			}
			return true;
		}

		/*
		 * If you want to avoid calculating the log-probability again
		 * when a proposal is rejected, override these methods to restore
		 * a cached value of the log-probability.
		 */
		@Override
		public boolean updateAfterRejection() {
			setLogP(oldLogLikelihood);
			return true;
		}

	}
}
