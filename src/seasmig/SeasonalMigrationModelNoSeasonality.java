package seasmig;

import java.util.Collection;

import cern.jet.random.Uniform;

import treelikelihood.ConstantMigrationBaseModel;
import treelikelihood.LikelihoodTree;
import treelikelihood.MigrationBaseModel;
import treelikelihood.SinusoidialSeasonalMigrationBaseModel;
import treelikelihood.TwoSeasonMigrationBaseModel;

import mc3kit.*;
import mc3kit.distributions.*;
import mc3kit.example.ExampleModel;
import static java.lang.String.format;
import static mc3kit.util.Math.*;


@SuppressWarnings("serial")
public class SeasonalMigrationModelNoSeasonality extends Model {


	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;	

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
		
		treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));

		DoubleDistribution ratePrior = new ExponentialDistribution(this,1.0);

		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j), ratePrior);
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

			LikelihoodTree workingCopy = data.getTrees().get(treeIndex.getValue()).copy(); 
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
