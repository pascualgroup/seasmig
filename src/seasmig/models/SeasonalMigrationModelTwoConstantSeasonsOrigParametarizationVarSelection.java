package seasmig.models;

import java.util.ArrayList;
import java.util.List;

import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;
import mc3kit.*;
import mc3kit.distributions.*;


@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsOrigParametarizationVarSelection extends SeasonalMigrationModel {


	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates1;	
	DoubleVariable[][] rates2;
	DoubleVariable seasonalPhase;
	IntVariable treeIndices[];
	LikelihoodVariable likeVar;
	private int nTrees[];
	private boolean fixPhase;
	private double seasonalPhaseRealization;
	private ExponentialDistribution ratePriorDist;	
	private BinaryVariable[][] rateIndicators1;
	private BinaryVariable[][] rateIndicators2;
	private BernoulliDistribution rateIndicatorPriorDist;

	protected SeasonalMigrationModelTwoConstantSeasonsOrigParametarizationVarSelection() { }

	public SeasonalMigrationModelTwoConstantSeasonsOrigParametarizationVarSelection(Chain initialChain, Config config, Data data, boolean fixPhase) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.data = data;		
		this.fixPhase = fixPhase;
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}
		rates1 = new DoubleVariable[numLocations][numLocations];
		rates2 = new DoubleVariable[numLocations][numLocations];
		rateIndicators1 = new BinaryVariable[numLocations][numLocations];
		rateIndicators2 = new BinaryVariable[numLocations][numLocations];
		
		beginConstruction();
		
		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}
		if (!fixPhase)
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));
		else
			seasonalPhaseRealization = config.fixedPhase;
		
		ratePriorDist = new ExponentialDistribution(this,"ratePrior");
		rateIndicatorPriorDist = new BernoulliDistribution(this, "rateIndicatorPriorDist",0.5);
		
		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates1[i][j] = new DoubleVariable(this, "rateParams1."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				rates2[i][j] = new DoubleVariable(this, "rateParams2."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				rateIndicators1[i][j] = new BinaryVariable(this, "rateIndicators1."+Integer.toString(i)+"."+Integer.toString(j), rateIndicatorPriorDist);
				rateIndicators2[i][j] = new BinaryVariable(this, "rateIndicators2."+Integer.toString(i)+"."+Integer.toString(j), rateIndicatorPriorDist);					
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {
	
		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsOrigParametarizationVarSelection m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length);

			// Add dependencies between likelihood variable and parameters
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1) {
					m.addEdge(this, m.treeIndices[i]);
				}
			}
			if (!fixPhase) {
				m.addEdge(this, m.seasonalPhase);
			}

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;				
					m.addEdge(this,rates1[i][j]);
					m.addEdge(this,rates2[i][j]);
					m.addEdge(this,rateIndicators1[i][j]);
					m.addEdge(this,rateIndicators2[i][j]);
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
						rates1doubleForm[i][j]=rates1[i][j].getValue()*(rateIndicators1[i][j].getValue() ? 1 : 0);						
						rates2doubleForm[i][j]=rates2[i][j].getValue()*(rateIndicators2[i][j].getValue() ? 1 : 0);
						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}
			
			// TODO: add update to migration model instead of reconstructing...
			MigrationBaseModel migrationBaseModel;
			if (!fixPhase)
				migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhase.getValue(),seasonalPhase.getValue()+0.5);
			else 
				migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhaseRealization,seasonalPhaseRealization+0.5);
			
			LikelihoodTree workingCopy;		
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1)
					workingCopy = data.getTrees().get(i).get((int)treeIndices[i].getValue()).copy(); 
				else
					workingCopy = data.getTrees().get(i).get(0).copy();
				workingCopy.setLikelihoodModel(migrationBaseModel);
				logP+=workingCopy.logLikelihood();
				trees[i]=workingCopy;
			}						

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
