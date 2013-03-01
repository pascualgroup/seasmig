package seasmig.models;

import java.util.Collections;
import java.util.Vector;

import seasmig.Config;
import seasmig.Data;
import seasmig.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;
import mc3kit.*;
import mc3kit.distributions.*;

@SuppressWarnings("serial")
public class SeasonalMigrationModelTwoConstantSeasonsVariableSelectionToFromGTR extends Model {

	Config config;
	Data data;
	int numLocations;

	DoubleVariable[][] rates;	
	DoubleVariable[] diffFromMultipliers;
	BinaryVariable[] diffFromIndicators;
	DoubleVariable[] diffToMultipliers;
	BinaryVariable[] diffToIndicators;
	DoubleVariable[][] diffMultipliers;
	BinaryVariable[][] diffIndicators;
	DoubleVariable[] baseFreq;

	DoubleVariable seasonalPhase;
	double seasonalPhaseRealization;
	IntVariable treeIndex;
	LikelihoodVariable likeVar;
	private int nTrees;	

	private boolean fixedPhase;
	private boolean diffFrom;
	private boolean diffTo;
	private boolean diffIndividual;
	private boolean GTR;

	protected SeasonalMigrationModelTwoConstantSeasonsVariableSelectionToFromGTR() { }

	public SeasonalMigrationModelTwoConstantSeasonsVariableSelectionToFromGTR(Chain initialChain, Config config, Data data, boolean fixedPhase, boolean diffIndividual, boolean diffFrom, boolean diffTo, boolean GTR) throws MC3KitException
	{
		// Either rows or columns or none of them can be set to have the same differential rates for season one vs. season two....
		super(initialChain);

		// TODO: there is no from to in gtr... 
		this.config = config;
		this.data = data;
		this.fixedPhase=fixedPhase;
		this.diffTo=diffTo;
		this.diffFrom=diffFrom;
		this.diffIndividual=diffIndividual;
		this.GTR=GTR;
		numLocations=data.getNumLocations();
		nTrees=data.getTrees().size();		
		rates = new DoubleVariable[numLocations][numLocations];

		diffFromMultipliers = new DoubleVariable[numLocations];
		diffFromIndicators = new BinaryVariable[numLocations];
		diffToMultipliers = new DoubleVariable[numLocations];
		diffToIndicators = new BinaryVariable[numLocations];
		diffMultipliers = new DoubleVariable[numLocations][numLocations];
		diffIndicators = new BinaryVariable[numLocations][numLocations];
		baseFreq= new DoubleVariable[numLocations-1]; // TODO: setup a Dirichlet prior 

		beginConstruction();

		treeIndex = new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0, nTrees-1));

		if (!fixedPhase) {
			seasonalPhase = new DoubleVariable(this,"seasonalPhase", new UniformDistribution(this,0,0.5));
		}
		else {
			seasonalPhaseRealization=config.fixedPhase;
		}

		DoubleDistribution ratePriorDist = new ExponentialDistribution(this,1.0);
		DoubleDistribution diffMultiplierPriorDist = new UniformDistribution(this,-1.0,1.0);
		BinaryDistribution diffIndicatorPriorDist = new BernoulliDistribution(this, 0.5);
		DoubleDistribution baseFreqPriorDist = new UniformDistribution(this,0.0,1.0);

		for (int i=0; i< numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null		
				if(GTR && i>j) continue; // half matrix for GTR
				rates[i][j] = new DoubleVariable(this, "meanRates."+Integer.toString(i)+"."+Integer.toString(j), ratePriorDist);
				if (diffIndividual) {
					diffMultipliers[i][j] = new DoubleVariable(this, "diffMultipliers."+Integer.toString(i)+"."+Integer.toString(j), diffMultiplierPriorDist);
					diffIndicators[i][j] = new BinaryVariable(this, "diffIndicators."+Integer.toString(i)+"."+Integer.toString(j), diffIndicatorPriorDist);
				}
			}
			if (diffFrom) {
				diffFromMultipliers[i] = new DoubleVariable(this, "diffFromMultipliers."+Integer.toString(i), diffMultiplierPriorDist);
				diffFromIndicators[i] = new BinaryVariable(this, "diffFromIndicators."+Integer.toString(i), diffIndicatorPriorDist);
			}
			if (diffTo) {
				diffToMultipliers[i] = new DoubleVariable(this, "diffToMultipliers."+Integer.toString(i), diffMultiplierPriorDist);
				diffToIndicators[i] = new BinaryVariable(this, "diffToIndicators."+Integer.toString(i), diffIndicatorPriorDist);
			}			
			if (GTR && i!=(numLocations-1)) {
				baseFreq[i] = new DoubleVariable(this, "baseFreq."+Integer.toString(i),baseFreqPriorDist);
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogLikelihood;
		private double logMaxLikelihood = Double.NEGATIVE_INFINITY;

		LikelihoodVariable(SeasonalMigrationModelTwoConstantSeasonsVariableSelectionToFromGTR m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			m.addEdge(this, m.treeIndex);
			if (!fixedPhase)
				m.addEdge(this, m.seasonalPhase);

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;
					if (GTR && i>j) continue;
					m.addEdge(this,rates[i][j]);
					if (diffIndividual) {
						m.addEdge(this,diffMultipliers[i][j]);
						m.addEdge(this,diffIndicators[i][j]);
					}
				}
				if (diffFrom) {
					m.addEdge(this,diffFromMultipliers[i]);
					m.addEdge(this,diffFromIndicators[i]);
				}
				if (diffTo) {
					m.addEdge(this,diffToMultipliers[i]);
					m.addEdge(this,diffToIndicators[i]);
				}
				if (GTR && i!=numLocations-1) {
					m.addEdge(this,baseFreq[i]);
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
						if (GTR && i>j) {
							rates1doubleForm[i][j]=rates[j][i].getValue();
							rates2doubleForm[i][j]=rates[j][i].getValue();
							if (diffFrom==true) {
								if (diffFromIndicators[i].getValue()==true) {						
									rates1doubleForm[i][j]*=(1-diffFromMultipliers[i].getValue());
									rates2doubleForm[i][j]*=(1+diffFromMultipliers[i].getValue());
								}
							}
							if (diffTo==true) {
								if (diffToIndicators[j].getValue()==true) {
									rates1doubleForm[i][j]*=(1-diffToMultipliers[j].getValue());
									rates2doubleForm[i][j]*=(1+diffToMultipliers[j].getValue());
								}
							}
							if (diffIndividual==true) {
								if (diffIndicators[j][i].getValue()==true) {
									rates1doubleForm[i][j]*=(1-diffMultipliers[j][i].getValue());
									rates2doubleForm[i][j]*=(1+diffMultipliers[j][i].getValue());
								}
							}
						}
						else {
							rates1doubleForm[i][j]=rates[i][j].getValue();
							rates2doubleForm[i][j]=rates[i][j].getValue();
							if (diffFrom==true) {
								if (diffFromIndicators[i].getValue()==true) {						
									rates1doubleForm[i][j]*=(1-diffFromMultipliers[i].getValue());
									rates2doubleForm[i][j]*=(1+diffFromMultipliers[i].getValue());
								}
							}
							if (diffTo==true) {
								if (diffToIndicators[j].getValue()==true) {
									rates1doubleForm[i][j]*=(1-diffToMultipliers[j].getValue());
									rates2doubleForm[i][j]*=(1+diffToMultipliers[j].getValue());
								}
							}
							if (diffIndividual==true) {
								if (diffIndicators[i][j].getValue()==true) {
									rates1doubleForm[i][j]*=(1-diffMultipliers[i][j].getValue());
									rates2doubleForm[i][j]*=(1+diffMultipliers[i][j].getValue());
								}
							}
						}

						rowsum1-=rates1doubleForm[i][j];
						rowsum2-=rates2doubleForm[i][j];
					}
				}
				rates1doubleForm[i][i]=rowsum1;
				rates2doubleForm[i][i]=rowsum2;
			}


			if (GTR) { 
				double[] baseFreqdoubleForm = convertUniformToDirichlet(baseFreq);		

				// Multiply by diags Q=R*Diag(Pi)
				double[][] newRates1doubleForm = new double[numLocations][numLocations];
				double[][] newRates2doubleForm = new double[numLocations][numLocations];

				for(int i = 0; i < numLocations; i++) { // aRow
					for(int j = 0; j < numLocations; j++) { // bColumn
						newRates1doubleForm[i][j] += rates1doubleForm[i][j] * baseFreqdoubleForm[j];					
						newRates2doubleForm[i][j] += rates2doubleForm[i][j] * baseFreqdoubleForm[j];
					} 
				}

				rates1doubleForm=newRates1doubleForm;
				rates2doubleForm=newRates2doubleForm;	
			}

			// TODO: add update to migration model instead of reconstructing...
			if (!fixedPhase)
				seasonalPhaseRealization=seasonalPhase.getValue();
			MigrationBaseModel migrationBaseModel = new TwoSeasonMigrationBaseModel(rates1doubleForm,rates2doubleForm,seasonalPhaseRealization,seasonalPhaseRealization+0.5);
			LikelihoodTree workingCopy = data.getTrees().get((int)treeIndex.getValue()).copy(); 
			workingCopy.setLikelihoodModel(migrationBaseModel);
			logLikelihood=workingCopy.logLikelihood();								

			setLogP(logLikelihood);			
			oldLogLikelihood=logLikelihood;
			if (logLikelihood>logMaxLikelihood) {
				logMaxLikelihood=logLikelihood;
			}
			return true;
		}


		private double[] convertUniformToDirichlet(DoubleVariable[] baseFreq) {
			// TODO: simplify this... and or move to mc3kit
			Vector<Double> baseFreqVector = new Vector<Double>(numLocations+1);
			baseFreqVector.add(0.0);
			for (int i=0;i<numLocations-1;i++) {
				baseFreqVector.add(baseFreq[i].getValue());									
			}

			Collections.sort(baseFreqVector);
			baseFreqVector.add(1.0);
			double[] baseFreqdoubleForm = new double[numLocations];				
			for (int i=1;i<(numLocations+1);i++) {
				baseFreqdoubleForm[i-1] = baseFreqVector.get(i)-baseFreqVector.get(i-1);		
			}
			return baseFreqdoubleForm;
		}

		/*
		 * If you want to avoid calculating the log-probability again
		 * when a proposal is rejected, override these methods to restore
		 * a cached value of the log-probability.
		 */
		@Override
		public boolean updateAfterRejection() {
			setLogP(oldLogLikelihood);
			return true;
		}

	}
}
