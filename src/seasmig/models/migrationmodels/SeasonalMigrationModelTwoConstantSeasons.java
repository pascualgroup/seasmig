package seasmig.models.migrationmodels;

import java.util.ArrayList;
import java.util.List;
import seasmig.models.*;
import mc3kit.Chain;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.UniformDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.migrationmain.Config;
import seasmig.models.TreesLikelihoodVariable;
import seasmig.data.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.TwoSeasonMigrationBaseModel;

@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasons extends MigrationModel {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	DoubleVariable[][] diffMultipliers;
	DoubleVariable[] diffMultipliersFrom;
	DoubleVariable[] diffMultipliersTo;
	DoubleVariable seasonalPhase;
	IntVariable treeIndices[];
	LikelihoodVariable likeVar;
	private int nTrees[];	

	boolean fixedPhase;
	boolean fixTo;
	boolean fixFrom;
	private ExponentialDistribution ratePriorDist;
	private double seasonStart;
	private DoubleVariable seasonalLength;
	private boolean fixedPhaseLength;

	protected SeasonalMigrationModelTwoConstantSeasons() { }

	public SeasonalMigrationModelTwoConstantSeasons(Chain initialChain, Config config, Data data, boolean fixedPhase, boolean fixedPhaseLength, boolean fixFrom, boolean fixTo) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		this.config = config;
		this.data = data;
		this.fixedPhase=fixedPhase;
		this.fixTo=fixTo;
		this.fixFrom=fixFrom;
		this.fixedPhaseLength=fixedPhaseLength;
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}

		rates = new DoubleVariable[numLocations][numLocations];

		diffMultipliersFrom = new DoubleVariable[numLocations];
		diffMultipliersTo = new DoubleVariable[numLocations];
		diffMultipliers = new DoubleVariable[numLocations][numLocations];

		beginConstruction();

		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}		

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

		ratePriorDist = new ExponentialDistribution(this,"ratePrior",1.0);
		DoubleDistribution diffMultiplierPriorDist = new UniformDistribution(this,-1.0,1.0);

		for (int i=0; i< numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "meanRates."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				if (!fixFrom && !fixTo) 
					diffMultipliers[i][j] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i)+"."+Integer.toString(j), diffMultiplierPriorDist);
			}		
		}

		if (fixFrom) {
			for(int i = 0; i < numLocations; i++) {
				diffMultipliersFrom[i] = new DoubleVariable(this, "diffMultipliersFrom."+Integer.toString(i), diffMultiplierPriorDist);
			}		
		}
		if (fixTo) {
			for(int i = 0; i < numLocations; i++) {
				diffMultipliersTo[i] = new DoubleVariable(this, "diffMultipliersTo."+Integer.toString(i), diffMultiplierPriorDist);
			}		
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends TreesLikelihoodVariable {


		private double seasonEnd;

		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasons m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true,nTrees.length,config);

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
				if (diffMultipliersFrom[i]!=null ) 
					m.addEdge(this, diffMultipliersFrom[i]);				
				if (diffMultipliersTo[i]!=null ) 
					m.addEdge(this, diffMultipliersTo[i]);				
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;
					m.addEdge(this,rates[i][j]);
					if (diffMultipliers[i][j]!=null) 
						m.addEdge(this,diffMultipliers[i][j]);		
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
						rates1doubleForm[i][j]=rates[i][j].getValue();
						rates2doubleForm[i][j]=rates[i][j].getValue();
						if (diffMultipliers[i][j]!=null) {
							rates1doubleForm[i][j]*=(1-diffMultipliers[i][j].getValue());
							rates2doubleForm[i][j]*=(1+diffMultipliers[i][j].getValue());
						}
						else if (diffMultipliersTo[j]!=null && (fixTo)) {
							rates1doubleForm[i][j]*=(1-diffMultipliersTo[j].getValue());
							rates2doubleForm[i][j]*=(1+diffMultipliersTo[j].getValue());
						}
						else if (diffMultipliersFrom[i]!=null && (fixFrom)) {
							rates1doubleForm[i][j]*=(1-diffMultipliersFrom[i].getValue());
							rates2doubleForm[i][j]*=(1+diffMultipliersFrom[i].getValue());
						}
						rowsum1+=rates1doubleForm[i][j];
						rowsum2+=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=-rowsum1;
				rates2doubleForm[i][i]=-rowsum2;		
			}

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
