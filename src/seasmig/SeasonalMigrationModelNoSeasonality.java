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
	Collection<LikelihoodTree> trees;
	int numLocations;

	DoubleVariable[][] rates;	
	IntVariable treeIndex;
	LikelihoodVariable likeVar;	

	protected SeasonalMigrationModelNoSeasonality() { }

	public SeasonalMigrationModelNoSeasonality(Chain initialChain, Config config, Collection<LikelihoodTree> trees) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.trees = trees;		
		numLocations=trees.iterator().next().getNumLocations();

		beginConstruction();

		new IntVariable(this, "treeIndex", new UniformIntDistribution(this, 0,trees.size()-1));

		new DoubleVariable(this, "ratePriorRate", new ExponentialDistribution(this, 1.0));
		new DoubleVariable(this, "ratePrior", new ExponentialDistribution(this,"ratePriorRate"));


		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j), new ExponentialDistribution(this,"ratePrior"));
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 

	private class LikelihoodVariable extends Variable {
		private double oldLogP;


		LikelihoodVariable(SeasonalMigrationModelNoSeasonality m) throws MC3KitException
		{
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true);

			// Add dependencies between likelihood variable and parameters
			m.addEdge(this, m.treeIndex);

			for(int i = 0; i < config.numLocations; i++) {
				for(int j = 0; j < config.numLocations; j++) {
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
			
			oldLogLikelihood = getLogP();

			double logLikelihood = 0.0;

			MigrationBaseModel migrationBaseModel = null;

			double[][] rates = new double[config.numLocations][config.numLocations];
				for (int i=0;i<config.numLocations;i++) {
					double rowsum=0;
					for (int j=0;j<config.numLocations;j++) {
						if (i!=j) {
							rates[i][j]=model.getRateParams(i, j).getRate();
							rowsum-=rates[i][j];
						}
					}
					rates[i][i]=rowsum;
				}

				migrationBaseModel = new ConstantMigrationBaseModel(rates);
				break;

			
				
		
			LikelihoodTree workingCopy = trees.get(treeIndex).copy(); 
			workingCopy.setLikelihoodModel(migrationBaseModel);
			logLikelihood=workingCopy.logLikelihood();
			

			setLogP(logLikelihood);

			// TODO: organize this...
			// Display state per hour....
			stateCount+=1;
			if (stateCount==25) {
				time = System.currentTimeMillis();
				System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/100L*1000000L/(60L*60L*1000L));
			}
			if (stateCount==250) {
				time = System.currentTimeMillis();
				System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/1000L*1000000L/(60L*60L*1000L));
			}
			if (stateCount%(config.printEveryNStates/config.chainCount)==0) {
				time = System.currentTimeMillis();
				System.out.printf("%d\t%.2f\t%.2f hours/million states\n",stateCount*config.chainCount,logLikelihood,(time-oldTime)/((double)config.chainCount)/10000L*1000000L/(60L*60L*1000L));
				oldTime=time;
			}
			double logP = 0.0;

			for(double x : data) {
				double w = weight.getValue();
				double[] logPs = new double[2];
				for(int i = 0; i < 2; i++) {
					double mean = means[i].getValue();
					double prec = precs[i].getValue();
					logPs[i] = NormalDistribution.getLogPPrecision(mean, prec, x);
				}
				double[] coeffs = new double[] { w, 1.0 - w };

				logP += logSumExp(logPs, coeffs);
			}
			setLogP(logP);
			oldLogP = logP;

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
