package seasmig.models.migrationmodels;

import java.util.ArrayList;
import java.util.List;

import mc3kit.Chain;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.migrationmain.Config;
import seasmig.models.*;
import seasmig.data.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.PiecewiseConstantMigrationBaseModel;

@SuppressWarnings("serial")
public class SeasonalMigrationModelNConstantSeasons extends MigrationModel {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][][] rates;			
	IntVariable treeIndices[];
	LikelihoodVariable likeVar;
	private int nTrees[];	

	private ExponentialDistribution ratePriorDist;
	private int nParts;

	protected SeasonalMigrationModelNConstantSeasons() { }

	public SeasonalMigrationModelNConstantSeasons(Chain initialChain, Config config, Data data, int nParts) throws MC3KitException
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

		beginConstruction();

		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}		

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);

		for (int i=0; i<nParts; i++) {
			for (int j=0; j<numLocations; j++) {
				for(int k=0; k<numLocations; k++) {
					if(j == k) continue; // rateParams[i,i] remains null			
					rates[i][j][k] = new DoubleVariable(this, "rates."+Integer.toString(i)+"."+Integer.toString(j)+"."+Integer.toString(k), ratePriorDist);
				}
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {


		LikelihoodVariable(SeasonalMigrationModelNConstantSeasons m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length,config);

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
							ratesDoubleForm[i][j][k]=rates[i][j][k].getValue();				
						}
					}
				}
			}

			TransitionModel migrationBaseModel = new PiecewiseConstantMigrationBaseModel(ratesDoubleForm,nParts);
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
