package seasmig.models;

import java.util.ArrayList;
import java.util.List;

import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.PiecewiseConstantMigrationBaseModel;
import mc3kit.*;
import mc3kit.distributions.*;

@SuppressWarnings("serial")
public class SeasonalMigrationModelNConstantSeasonsVarSelect extends SeasonalMigrationModel {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][][] rates;			
	IntVariable treeIndices[];
	BinaryVariable[][][] sameIndicators;
	LikelihoodVariable likeVar;
	private int nTrees[];	

	private ExponentialDistribution ratePriorDist;
	private int nParts;
	private BernoulliDistribution indicatorPriorDist;
	//private DoubleVariable indicatorHyperprior;
	
	protected SeasonalMigrationModelNConstantSeasonsVarSelect() { }

	public SeasonalMigrationModelNConstantSeasonsVarSelect(Chain initialChain, Config config, Data data, int nParts) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.config = config;
		this.data = data;
		this.nParts = nParts;
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}

		rates = new DoubleVariable[nParts][numLocations][numLocations];
		sameIndicators = new BinaryVariable[nParts][numLocations][numLocations];

		beginConstruction();

		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}		

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);
		//indicatorHyperprior = new DoubleVariable(this, "indicatorHyperprior", new UniformDistribution(this));
		indicatorPriorDist = new BernoulliDistribution(this, "indicatorPrior",0.5);

		for (int i=0; i<nParts; i++) {
			for (int j=0; j<numLocations; j++) {
				for(int k=0; k<numLocations; k++) {
					if(j == k) continue; // rateParams[i,i] remains null			
					rates[i][j][k] = new DoubleVariable(this, "rates."+Integer.toString(i)+"."+Integer.toString(j)+"."+Integer.toString(k), ratePriorDist);
					if (i>=1) sameIndicators[i][j][k] = new BinaryVariable(this, "indicators."+Integer.toString(i)+"."+Integer.toString(j)+"."+Integer.toString(k), indicatorPriorDist);
				}
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {


		LikelihoodVariable(SeasonalMigrationModelNConstantSeasonsVarSelect m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length);

			// Add dependencies between likelihood variable and parameters
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1) {
					m.addEdge(this, m.treeIndices[i]);
				}
			}
			
			for (int i=0; i<nParts; i++) {
				for (int j=0; j<numLocations; j++) {
					for(int k=0; k<numLocations; k++) {
						if(j == k) continue; // rateParams[i,i] remains null			
						m.addEdge(this,rates[i][j][k]);
						if (i>=1) m.addEdge(this,sameIndicators[i][j][k]);
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
					for (int k=0;k<numLocations;k++) {
						if (j!=k) {
							if (i>=1) {
								ratesDoubleForm[i][j][k]=(sameIndicators[i][j][k].getValue() ? ratesDoubleForm[i-1][j][k] : rates[i][j][k].getValue());
							}
							else {
								ratesDoubleForm[i][j][k]=rates[i][j][k].getValue();
							}
						}
					}
				}
			}

			MigrationBaseModel migrationBaseModel = new PiecewiseConstantMigrationBaseModel(ratesDoubleForm,nParts);
			LikelihoodTree workingCopy;		
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1)
					workingCopy = data.getTrees().get(i).get((int)treeIndices[i].getValue()).copy(); 
				else
					workingCopy = data.getTrees().get(i).get(0).copy();
				workingCopy.setLikelihoodModel(migrationBaseModel);
				logLikelihood+=workingCopy.logLikelihood();
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
