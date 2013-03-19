package seasmig.models;

import seasmig.Config;
import seasmig.Data;
import seasmig.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;
import seasmig.util.Util;
import mc3kit.*;
import mc3kit.distributions.*;


@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsOrigParametarization extends Model {


	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates1;	
	DoubleVariable[][] rates2;
	DoubleVariable seasonalPhase;
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;
	private boolean fixPhase;
	private double seasonalPhaseRealization;	

	protected SeasonalMigrationModelTwoConstantSeasonsOrigParametarization() { }

	public SeasonalMigrationModelTwoConstantSeasonsOrigParametarization(Chain initialChain, Config config, Data data, boolean fixPhase) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.data = data;		
		this.fixPhase = fixPhase;
		numLocations=data.getNumLocations();
		nTrees=data.getTrees().size();		
		rates1 = new DoubleVariable[numLocations][numLocations];
		rates2 = new DoubleVariable[numLocations][numLocations];
		
		beginConstruction();
		
		treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));
		
		if (!fixPhase)
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));
		else
			seasonalPhaseRealization = config.fixedPhase;
		DoubleDistribution ratePrior = new ExponentialDistribution(this,1.0);

		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates1[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j), ratePrior);
				rates2[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j), ratePrior);
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogP;


		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsOrigParametarization m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			m.addEdge(this, m.treeIndex);
			
			if (!fixPhase) {
				m.addEdge(this, m.seasonalPhase);
			}

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;				
					m.addEdge(this,rates1[i][j]);
					m.addEdge(this,rates2[i][j]);
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
						rates1doubleForm[i][j]=rates1[i][j].getValue();
						rates2doubleForm[i][j]=rates2[i][j].getValue();
						rowsum1-=rates1[i][j].getValue();
						rowsum2-=rates2[i][j].getValue();
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}

			// TODO: check this...
			Util.transposeSquareMatrix(rates1doubleForm);
			Util.transposeSquareMatrix(rates2doubleForm);
			
			// TODO: add update to migration model instead of reconstructing...
			MigrationBaseModel migrationBaseModel;
			if (!fixPhase)
				migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhase.getValue(),seasonalPhase.getValue()+0.5);
			else 
				migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhaseRealization,seasonalPhaseRealization+0.5);
			
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
