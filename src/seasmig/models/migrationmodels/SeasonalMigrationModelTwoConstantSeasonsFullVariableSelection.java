package seasmig.models.migrationmodels;

import java.util.ArrayList;
import java.util.List;
import seasmig.models.*;
import mc3kit.BinaryDistribution;
import mc3kit.BinaryVariable;
import mc3kit.Chain;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.distributions.BernoulliDistribution;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.migrationmain.Config;
import seasmig.models.TreesLikelihoodVariable;
import seasmig.data.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;
import seasmig.treelikelihood.transitionmodels.TwoSeasonMigrationBaseModel;

@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsFullVariableSelection extends MigrationModel {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	DoubleVariable[][] diffMultipliers;
	BinaryVariable[][] diffIndicators;

	DoubleVariable seasonalPhase;
	DoubleVariable seasonalLength;
	double seasonStart;
	double seasonEnd;
	IntVariable treeIndices[];
	LikelihoodVariable likeVar;
	private int nTrees[];	

	boolean fixedPhase;
	boolean fixedPhaseLength;
	boolean fixRate;
	private ExponentialDistribution ratePriorDist;
	private BinaryVariable[][] rateIndicators;
	private DoubleVariable rateHyperPrior;

	protected SeasonalMigrationModelTwoConstantSeasonsFullVariableSelection() { }

	public SeasonalMigrationModelTwoConstantSeasonsFullVariableSelection(Chain initialChain, Config config, Data data, boolean fixedPhase, boolean fixedPhaseLength, boolean fixRate) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.config = config;
		this.data = data;
		this.fixedPhase=fixedPhase;
		this.fixedPhaseLength=fixedPhaseLength;
		this.fixRate=fixRate;
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}	
		rates = new DoubleVariable[numLocations][numLocations];

		diffMultipliers = new DoubleVariable[numLocations][numLocations];
		diffIndicators = new BinaryVariable[numLocations][numLocations];
		rateIndicators = new BinaryVariable[numLocations][numLocations];

		beginConstruction();

		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}
				
		rateHyperPrior= new DoubleVariable(this, "rateHyperPrior",  new ExponentialDistribution(this));
		ratePriorDist = new ExponentialDistribution(this,"ratePriorDist");
		ratePriorDist.setRate(rateHyperPrior);		

		if (fixedPhase && fixedPhaseLength) {
			seasonStart=config.fixedPhase;			
		}
		else if (!fixedPhase && fixedPhaseLength) {
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));			
		}
		else if (!fixedPhase && !fixedPhaseLength) {	
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,1));
			seasonalLength = new DoubleVariable(this,"seasonalLength", new UniformDistribution(this,config.minSeasonLength,0.5));			
		}
		else /* fixedPhase && !fixedLength */ {
			seasonStart=config.fixedPhase;
			seasonalLength = new DoubleVariable(this,"seasonalLength", new UniformDistribution(this,config.minSeasonLength,0.5));			
		}

		DoubleDistribution diffMultiplierPriorDist = new UniformDistribution(this,-1.0,1.0);
		BinaryDistribution diffIndicatorPriorDist = new BernoulliDistribution(this, 0.5);
		BinaryDistribution rateIndicatorPriorDist = new BernoulliDistribution(this, 0.5);

		for (int i=0; i< numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "meanRates."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				diffMultipliers[i][j] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i)+"."+Integer.toString(j), diffMultiplierPriorDist);
				diffIndicators[i][j] = new BinaryVariable(this, "diffIndicators."+Integer.toString(i)+"."+Integer.toString(j), diffIndicatorPriorDist);
				rateIndicators[i][j] = new BinaryVariable(this, "rateIndicators."+Integer.toString(i)+"."+Integer.toString(j), rateIndicatorPriorDist);
			}		
		}


		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {
	
		private double infinitesimalRate = 1E-3;

		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsFullVariableSelection m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length,config);

			// TODO: check if this is required
			m.addEdge(ratePriorDist,rateHyperPrior);
			
			// Add dependencies between likelihood variable and parameters
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1) {
					m.addEdge(this, m.treeIndices[i]);
				}
			}
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
					m.addEdge(this,rateIndicators[i][j]);
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

			double[][] rates1doubleForm = new double[numLocations][numLocations];
			double[][] rates2doubleForm = new double[numLocations][numLocations];			

			for (int i=0;i<numLocations;i++) {
				double rowsum1=0;
				double rowsum2=0;
				for (int j=0;j<numLocations;j++) {
					if (i!=j) {
						if (rateIndicators[i][j].getValue()==true) {
							if (diffIndicators[i][j].getValue()==true) {
								rates1doubleForm[i][j]=rates[i][j].getValue()*(1-diffMultipliers[i][j].getValue());
								rates2doubleForm[i][j]=rates[i][j].getValue()*(1+diffMultipliers[i][j].getValue());
							}
							else {
								rates1doubleForm[i][j]=rates[i][j].getValue();
								rates2doubleForm[i][j]=rates[i][j].getValue();
							}
						}
						else {
							rates1doubleForm[i][j]=infinitesimalRate*(rates[i][j].getValue()%1.0);
							rates2doubleForm[i][j]=infinitesimalRate*(rates[i][j].getValue()%1.0);
						}
						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}

			if (fixRate) {
				double[][] ratesdoubleForm = new double[numLocations][numLocations];
				
				for (int i=0;i<numLocations;i++) {
					double rowsum=0;
					for (int j=0;j<numLocations;j++) {
						if (i!=j) 
							ratesdoubleForm[i][j]=rates[i][j].getValue();						
						rowsum-=ratesdoubleForm[i][j];
					}
					ratesdoubleForm[i][i]=rowsum;
				}
				
				double rate=getRate(ratesdoubleForm);	
				adjustRate(rates1doubleForm, rate);
				adjustRate(rates2doubleForm, rate);
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

			TransitionModel migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonStart,seasonEnd);
			LikelihoodTree workingCopy;		
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1)
					workingCopy = data.getTrees().get(i).get((int)treeIndices[i].getValue()).copy(); 
				else
					workingCopy = data.getTrees().get(i).get(0).copy();
				workingCopy.setMigrationModel(migrationBaseModel);
				logLikelihood+=config.treeWeights[i]*workingCopy.logLikelihood();
				trees[i]=workingCopy;
			}						

			setLogP(logLikelihood);				
			oldLogP=logLikelihood;
			return true;
		}

		private void adjustRate(double[][] Q, double adjustedRate) {
			double inputRate = getRate(Q);
			for (int i=0;i<Q.length;i++) {
				for (int j=0;j<Q.length;j++) {
					Q[i][j]=Q[i][j]/inputRate*adjustedRate;
				}
			}
		}

		private double getRate(double[][] ratesdoubleForm) {
			MatrixExponentiator matrixExp = new Matlab7MatrixExp(ratesdoubleForm);
			double[] pi = matrixExp.expm(config.veryLongTime).viewRow(1).toArray();	
			double rate=0;
			for (int i=0;i<numLocations;i++) {
				rate-=pi[i]*ratesdoubleForm[i][i];
			}
			return rate;
		}

		/*
		 * If you want to avoid calculating the log-probability again
		 * when a proposal is rejected, override these methods to restore
		 * a cached value of the log-probability.
		 */
		@Override
		public boolean updateAfterRejection() {
			if (!firstCall)
				setLogP(oldLogP);
			else {
				super.updateAfterRejection();
				firstCall=false;
			}
			return true;
		}

	}
}
