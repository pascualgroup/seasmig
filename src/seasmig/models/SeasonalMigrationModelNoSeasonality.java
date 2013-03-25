package seasmig.models;

import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.ConstantMigrationBaseModel;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import mc3kit.*;
import mc3kit.distributions.*;


@SuppressWarnings("serial")
public class SeasonalMigrationModelNoSeasonality extends Model {


	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;
	private ExponentialDistribution ratePriorDist;	

	protected SeasonalMigrationModelNoSeasonality() { }

	public SeasonalMigrationModelNoSeasonality(Chain initialChain, Config config, Data data) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.data = data;		
		numLocations=data.getNumLocations();
		nTrees=data.getTrees().size();		
		rates = new DoubleVariable[numLocations][numLocations];
		
		beginConstruction();
		
		if (nTrees>1)
			treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));

		ratePriorDist = new ExponentialDistribution(this,"ratePrior");
		
		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j),ratePriorDist);
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogP;


		LikelihoodVariable(SeasonalMigrationModelNoSeasonality m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			if (nTrees>1)
				m.addEdge(this, m.treeIndex);

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;				
					m.addEdge(this,rates[i][j]);
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

			double logP = 0.0;

			double[][] ratesdoubleForm = new double[numLocations][numLocations];
			for (int i=0;i<numLocations;i++) {
				double rowsum=0;
				for (int j=0;j<numLocations;j++) {
					if (i!=j) {
						ratesdoubleForm[i][j]=rates[i][j].getValue();
						rowsum-=rates[i][j].getValue();
					}
				}
				ratesdoubleForm[i][i]=rowsum;
			}

			// TODO: add update to migration model instead of reconstructing...
			MigrationBaseModel migrationBaseModel = new ConstantMigrationBaseModel(ratesdoubleForm);
			LikelihoodTree workingCopy;
			if (nTrees>1)
				workingCopy = data.getTrees().get((int)treeIndex.getValue()).copy(); 
			else
				workingCopy = data.getTrees().get(0).copy();
			workingCopy.setLikelihoodModel(migrationBaseModel);
			logP=workingCopy.logLikelihood();
			
			setLogP(logP);
			oldLogP=logP;
			return true;
		}


		/*
		 * If you want to avoid calculating the log-probability again
		 * when a proposal is rejected, override these methods to restore
		 * a cached value of the log-probability.
		 */
		@Override
		public boolean updateAfterRejection() {
			setLogP(oldLogP);
			return true;
		}

	}
}
