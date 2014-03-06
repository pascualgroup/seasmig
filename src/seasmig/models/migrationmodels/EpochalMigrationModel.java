package seasmig.models.migrationmodels;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import mc3kit.BinaryVariable;
import mc3kit.Chain;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.distributions.BernoulliDistribution;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.data.Data;
import seasmig.migrationmain.Config;
import seasmig.models.*;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.EpochalMigrationBaseModel;

@SuppressWarnings("serial")
public class EpochalMigrationModel extends MigrationModel {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][][] rates;
	BinaryVariable[][][] rateIndicators;
	DoubleVariable[] epochs;
	IntVariable treeIndices[];
	LikelihoodVariable likeVar;
	private int nTrees[];	

	private ExponentialDistribution ratePriorDist;
	private int nParts;
	private boolean freeTimes;
	private boolean vs;
	private BernoulliDistribution rateIndicatorPriorDist;

	protected EpochalMigrationModel() { }

	public EpochalMigrationModel(Chain initialChain, Config config, Data data, boolean freeTimes_, boolean vs_) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.freeTimes=freeTimes_;
		this.vs = vs_;
		this.config = config;
		this.data = data;
		this.nParts = config.nEpochs;
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}

		rates = new DoubleVariable[nParts][numLocations][numLocations];
		
		if (vs) {
			rateIndicators = new BinaryVariable[nParts][numLocations][numLocations];
		}

		if (freeTimes) {
			epochs= new DoubleVariable[nParts-1];
		}
		
		beginConstruction();
		
		if (freeTimes) {
			// TODO: add sort and Dirichlet distribution
			for (int i=0;i<epochs.length;i++) {
				epochs[i]=new DoubleVariable(this, "epochTime."+i, new UniformDistribution(this, config.minEpochTime, config.maxEpochTime));
			}
		}
		
		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}		

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);
		if (vs) rateIndicatorPriorDist = new BernoulliDistribution(this, "rateIndicatorPrior", 0.5);

		for (int i=0; i<nParts; i++) {
			for (int j=0; j<numLocations; j++) {
				for(int k=0; k<numLocations; k++) {
					if(j == k) continue; // rateParams[i,i] remains null			
					rates[i][j][k] = new DoubleVariable(this, "rates."+Integer.toString(i)+"."+Integer.toString(j)+"."+Integer.toString(k), ratePriorDist);
					if (vs) rateIndicators[i][j][k] = new BinaryVariable(this, "rateIndicators."+Integer.toString(i)+"."+Integer.toString(j)+"."+Integer.toString(k), rateIndicatorPriorDist);					
				}
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {
		
		private double infinitesimalRate = 1E-4;

		LikelihoodVariable(EpochalMigrationModel m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length,config);

			// Add dependencies between likelihood variable and parameters
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1) {
					m.addEdge(this, m.treeIndices[i]);
				}
			}
			
			if (freeTimes) {
				// TODO: add sort and Dirichlet distribution
				for (int i=0;i<(nParts-1);i++) {
					m.addEdge(this,m.epochs[i]);
				}
			}
			
			// Should priors be included here? guess not...
			// if (vs) m.addEdge(this,rateIndicatorPriorDist);
			// m.addEdge(this,ratePriorDist);
			
			for (int i=0; i<nParts; i++) {
				for (int j=0; j<numLocations; j++) {
					for(int k=0; k<numLocations; k++) {
						if(j == k) continue; // rateParams[i,i] remains null			
						m.addEdge(this,rates[i][j][k]);
						if (vs) m.addEdge(this, rateIndicators[i][j][k]);
					}
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

			double[][][] ratesDoubleForm = new double[nParts][numLocations][numLocations];

			for (int i=0;i<nParts;i++) {
				for (int j=0;j<numLocations;j++) {	
					double rowsum=0;
					for (int k=0;k<numLocations;k++) {
						if (j!=k) {
							if (vs) {
								ratesDoubleForm[i][j][k]=(rateIndicators[i][j][k].getValue() ? rates[i][j][k].getValue() : infinitesimalRate*(rates[i][j][k].getValue()%1.0));
							}								
							else {
								ratesDoubleForm[i][j][k]=rates[i][j][k].getValue();
							}
						}		
						rowsum-=ratesDoubleForm[i][j][k];						
					}
					ratesDoubleForm[i][j][j]=rowsum;
				}
			}

			double[] epochsDoubleForm;
			
			if (!freeTimes) {
				epochsDoubleForm = config.epochTimes;
			}
			else {
				// TODO: add Dirichlet ... 
				Vector<Double> epochsDoubleFormVector = new Vector<Double>(nParts-1);
				for (int i=0;i<epochs.length;i++) {
					epochsDoubleFormVector.add(epochs[i].getValue());
				}
				Collections.sort(epochsDoubleFormVector);
				epochsDoubleForm = new double[nParts-1];
				for (int i=0;i<epochsDoubleForm.length;i++)
					epochsDoubleForm[i] = epochsDoubleFormVector.get(i);
			}
	
			TransitionModel migrationBaseModel = new EpochalMigrationBaseModel(ratesDoubleForm,epochsDoubleForm);
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
