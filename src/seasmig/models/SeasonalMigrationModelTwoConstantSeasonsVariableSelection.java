package seasmig.models;

import mc3kit.BinaryDistribution;
import mc3kit.BinaryVariable;
import mc3kit.Chain;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.Variable;
import mc3kit.distributions.BernoulliDistribution;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.Config;
import seasmig.Data;
import seasmig.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;

@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsVariableSelection extends Model {

	Config config;
	Data data;
	int numLocations;
	final static double minSeasonalWindowLength = 0.083333333333333*2.0; // 2 month 	

	DoubleVariable[][] rates;	
	DoubleVariable[][] diffMultipliers;
	BinaryVariable[][] diffIndicators;
	
	DoubleVariable seasonalPhase;
	DoubleVariable seasonalLength;
	double seasonStart;
	double seasonEnd;
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;	

	boolean fixedPhase;
	boolean fixedPhaseLength;
	private ExponentialDistribution ratePriorDist;

	protected SeasonalMigrationModelTwoConstantSeasonsVariableSelection() { }

	public SeasonalMigrationModelTwoConstantSeasonsVariableSelection(Chain initialChain, Config config, Data data, boolean fixedPhase, boolean fixedPhaseLength) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.config = config;
		this.data = data;
		this.fixedPhase=fixedPhase;
		this.fixedPhaseLength=fixedPhaseLength;
		numLocations=data.getNumLocations();
		nTrees=data.getTrees().size();		
		rates = new DoubleVariable[numLocations][numLocations];

		diffMultipliers = new DoubleVariable[numLocations][numLocations];
		diffIndicators = new BinaryVariable[numLocations][numLocations];
		
		beginConstruction();

		if (nTrees>1)
			treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);
		
		if (fixedPhase && fixedPhaseLength) {
			seasonStart=config.fixedPhase;			
		}
		else if (!fixedPhase && fixedPhaseLength) {
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));			
		}
		else if (!fixedPhase && !fixedPhaseLength) {	
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,1));
			seasonalLength = new DoubleVariable(this,"seasonalLength", new UniformDistribution(this,minSeasonalWindowLength,0.5));			
		}
		else /* fixedPhase && !fixedLength */ {
			seasonStart=config.fixedPhase;
			seasonalLength = new DoubleVariable(this,"seasonalLength", new UniformDistribution(this,minSeasonalWindowLength,0.5));			
		}
		
		DoubleDistribution diffMultiplierPriorDist = new UniformDistribution(this,-1.0,1.0);
		BinaryDistribution diffIndicatorPriorDist = new BernoulliDistribution(this, 0.5);
		
		for (int i=0; i< numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "meanRates."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				diffMultipliers[i][j] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i)+"."+Integer.toString(j), diffMultiplierPriorDist);
				diffIndicators[i][j] = new BinaryVariable(this, "diffIndicators."+Integer.toString(i)+"."+Integer.toString(j), diffIndicatorPriorDist);
			}		
		}


		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogLikelihood;
		private double logMaxLikelihood = Double.NEGATIVE_INFINITY;

		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsVariableSelection m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			if (nTrees>1)
				m.addEdge(this, m.treeIndex);
			if (!fixedPhase)
				m.addEdge(this, m.seasonalPhase);
			if (!fixedPhaseLength)
				m.addEdge(this, m.seasonalLength);

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;
					m.addEdge(this,rates[i][j]);
					m.addEdge(this,diffMultipliers[i][j]);
					m.addEdge(this,diffIndicators[i][j]);
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
						if (diffIndicators[i][j].getValue()==true) {
							rates1doubleForm[i][j]=rates[i][j].getValue()*(1-diffMultipliers[i][j].getValue());
							rates2doubleForm[i][j]=rates[i][j].getValue()*(1+diffMultipliers[i][j].getValue());
						}
						else {
							rates1doubleForm[i][j]=rates[i][j].getValue();
							rates2doubleForm[i][j]=rates[i][j].getValue();
						}
						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}
			
			// TODO: add update to migration model instead of reconstructing...
			if (fixedPhase && fixedPhaseLength) {
				seasonStart=config.fixedPhase;	
				seasonEnd = seasonStart+0.5;
			}
			else if (!fixedPhase && fixedPhaseLength) {
				seasonStart=seasonalPhase.getValue();
				seasonEnd=seasonalPhase.getValue()+0.5;			
			}
			else if (!fixedPhase && !fixedPhaseLength) {
				if (seasonalPhase.getValue()+seasonalLength.getValue()<1) {
					seasonStart=seasonalPhase.getValue();
					seasonEnd=seasonalPhase.getValue()+seasonalLength.getValue();
				} 
				else {
					seasonStart=seasonalPhase.getValue()+seasonalLength.getValue()-1.0;
					seasonEnd=seasonalPhase.getValue();
				}						
			}
			else /* fixedPhase && !fixedLength */ {					
				if (config.fixedPhase+seasonalLength.getValue()<1) {
					seasonStart=config.fixedPhase;
					seasonEnd=seasonStart+seasonalLength.getValue();
				} 
				else {
					seasonStart=config.fixedPhase+seasonalLength.getValue()-1.0;
					seasonEnd=config.fixedPhase;
				}			
			}			
			
			MigrationBaseModel migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonStart,seasonEnd);
			LikelihoodTree workingCopy;
			if (nTrees>1)
				workingCopy = data.getTrees().get((int)treeIndex.getValue()).copy(); 
			else
				workingCopy = data.getTrees().get(0).copy();

			workingCopy.setLikelihoodModel(migrationBaseModel);
			logLikelihood=workingCopy.logLikelihood();								

			// TODO: think about how to integrate over more than one tree
			
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
