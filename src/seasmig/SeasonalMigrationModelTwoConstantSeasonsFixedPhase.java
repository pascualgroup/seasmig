package seasmig;

import treelikelihood.LikelihoodTree;
import treelikelihood.MigrationBaseModel;
import treelikelihood.TwoSeasonMigrationBaseModel;

import mc3kit.*;
import mc3kit.distributions.*;


@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsFixedPhase extends Model {


	Config config;
	Data data;
	int numLocations;
	double seasonalPhase;

	DoubleVariable[][] rates;	
	DoubleVariable[][] logRatios;
	
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;	

	protected SeasonalMigrationModelTwoConstantSeasonsFixedPhase() { }

	public SeasonalMigrationModelTwoConstantSeasonsFixedPhase(Chain initialChain, Config config, Data data) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.data = data;		
		numLocations=data.getNumLocations();
		seasonalPhase = config.fixedPhase;
		
		nTrees=data.getTrees().size();		
		rates = new DoubleVariable[numLocations][numLocations];
		logRatios = new DoubleVariable[numLocations][numLocations];
		
		beginConstruction();
		
		treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));
		DoubleDistribution logRatioPrior = new UniformDistribution(this,-10.0,10.0);
		
		DoubleDistribution ratePrior = new ExponentialDistribution(this,1.0);

		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j), ratePrior);
				logRatios[i][j] = new DoubleVariable(this, "ratioParams."+Integer.toString(i)+"."+Integer.toString(j), logRatioPrior);
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogP;


		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsFixedPhase m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			m.addEdge(this, m.treeIndex);
		
			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;				
					m.addEdge(this,rates[i][j]);
					m.addEdge(this,logRatios[i][j]);
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

			double[][] rates1doubleForm = new double[numLocations][numLocations];
			double[][] rates2doubleForm = new double[numLocations][numLocations];
			for (int i=0;i<numLocations;i++) {
				double rowsum1=0;
				double rowsum2=0;
				for (int j=0;j<numLocations;j++) {
					if (i!=j) {
						rates1doubleForm[i][j]=rates[i][j].getValue()*Math.exp(-logRatios[i][j].getValue()/2);
						rates2doubleForm[i][j]=rates[i][j].getValue()*Math.exp(logRatios[i][j].getValue()/2);
						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}

			// TODO: add update to migration model instead of reconstructing...
			MigrationBaseModel migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhase,seasonalPhase+0.5);

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