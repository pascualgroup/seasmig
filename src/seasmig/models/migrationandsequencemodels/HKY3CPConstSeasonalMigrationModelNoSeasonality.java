package seasmig.models.migrationandsequencemodels;

import java.util.ArrayList;
import java.util.List;

import mc3kit.Chain;
import mc3kit.DoubleVariable;
import mc3kit.IntVariable;
import mc3kit.MC3KitException;
import mc3kit.distributions.ExponentialDistribution;
import mc3kit.distributions.NormalDistribution;
import mc3kit.distributions.UniformDistribution;
import mc3kit.distributions.UniformIntDistribution;
import seasmig.data.Data;
import seasmig.migrationmain.Config;
import seasmig.models.MigrationModel;
import seasmig.models.TreesLikelihoodVariable;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.ConstantTransitionBaseModel;
import seasmig.util.Util;


@SuppressWarnings("serial")
public class HKY3CPConstSeasonalMigrationModelNoSeasonality extends MigrationModel {

	// TODO: seperate into two variables for seqence and location likelihood...

	Config config;
	Data data;
	int numLocations;

	// HKY sequence model parameters
	DoubleVariable[] logk;	
	DoubleVariable[][] forPis; // will be converted to Dirichlet pis
	DoubleVariable[] mu;

	// Constant migration model parameters
	DoubleVariable[][] rates;
	
	// Trees
	IntVariable treeIndices[];
	
	// Likelihood Variable
	LikelihoodVariable likeVar;
	
	// Realization trees
	private int nTrees[];
	
	// Priors
	private ExponentialDistribution ratePriorDist;
	private ExponentialDistribution muPriorDist;
	private NormalDistribution logkPriorDist;
	private ExponentialDistribution rateHyperPriorDist;
	private DoubleVariable rateHyperPrior;
	private ExponentialDistribution muHyperPriorDist;
	private DoubleVariable muHyperPrior;
	protected HKY3CPConstSeasonalMigrationModelNoSeasonality() { }

	public HKY3CPConstSeasonalMigrationModelNoSeasonality(Chain initialChain, Config config, Data data) throws MC3KitException
	{
		super(initialChain);
		this.config = config;
		this.data = data;		
		numLocations=data.getNumLocations();
		List<ArrayList<LikelihoodTree>> trees = data.getTrees();		
		nTrees = new int[trees.size()];
		for (int i=0;i<trees.size();i++) {
			nTrees[i]=trees.get(i).size();
		}
		rates = new DoubleVariable[numLocations][numLocations];
		forPis = new DoubleVariable[3][3]; // Codon position, rest will be converted to piC, piA, piG
		logk = new DoubleVariable[3]; // Codon position
		mu = new DoubleVariable[3]; // Codon position

		beginConstruction();
		treeIndices = new IntVariable[trees.size()];
		for (int i=0;i<trees.size();i++) {
			if (nTrees[i]>1) {
				treeIndices[i] = new IntVariable(this, "treeIndex."+i, new UniformIntDistribution(this, 0, nTrees[i]-1));
			}
		}
		
		rateHyperPriorDist = new ExponentialDistribution(this,"rateHyperPriorDist");
		rateHyperPrior= new DoubleVariable(this, "rateHyperPrior", rateHyperPriorDist);
		ratePriorDist = new ExponentialDistribution(this,"ratePriorDist");
		ratePriorDist.setRate(rateHyperPrior);		
		// TODO: add 1/x distribution...
		muHyperPriorDist = new ExponentialDistribution(this,"muHyperPriorDist");
		muHyperPrior = new DoubleVariable(this, "muHyperPrior", muHyperPriorDist);
		muPriorDist = new ExponentialDistribution(this,"muPrior");
		muPriorDist.setRate(muHyperPrior);
		// TODO: add Log Normal distribution...
		logkPriorDist = new NormalDistribution(this,"logkPrior",0,50);

		for(int i = 0; i < numLocations; i++) {
			for(int j = 0; j < numLocations; j++) {
				if(i == j) continue; // rateParams[i,i] remains null			
				rates[i][j] = new DoubleVariable(this, "rateParams."+Integer.toString(i)+"."+Integer.toString(j),ratePriorDist);
			}
		}
		
		for(int i = 0; i < 3; i++) {
			mu[i]= new DoubleVariable(this, "mu."+Integer.toString(i),muPriorDist);
		}
		
		for(int i = 0; i < 3; i++) {
			logk[i]= new DoubleVariable(this, "logk."+Integer.toString(i),logkPriorDist);
		}
		
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				forPis[i][j] = new DoubleVariable(this, "forPis."+Integer.toString(i)+"."+Integer.toString(j),new UniformDistribution(this));
			}
		}

		// Custom likelihood variable
		likeVar = new LikelihoodVariable(this);

		endConstruction();

	} 	

	private class LikelihoodVariable extends TreesLikelihoodVariable {
	

		LikelihoodVariable(HKY3CPConstSeasonalMigrationModelNoSeasonality m) throws MC3KitException {
			// Call superclass constructor specifying that this is an
			// OBSERVED random variable (true for last parameter).
			super(m, "likeVar", true, nTrees.length,config);

			// TODO: make sure you don't need to add m.addEdge(m.muHyperPrior,m.muPriorDist); etc... ??
			
			// Add dependencies between likelihood variable and parameters
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1) {
					m.addEdge(this, m.treeIndices[i]);
				}
			}

			for(int i = 0; i < numLocations; i++) {
				for(int j = 0; j < numLocations; j++) {
					if (i==j) continue;				
					m.addEdge(this,rates[i][j]);
				}
			}
			
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {			
					m.addEdge(this,forPis[i][j]);
				}
			}
			
			for(int i = 0; i < 3; i++) {			
				m.addEdge(this,logk[i]);
			}
			
			for(int i = 0; i < 3; i++) {			
				m.addEdge(this,mu[i]);
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
			
			double[] mudoubleForm = new double[3];
			for (int i=0;i<3;i++) {				
				mudoubleForm[i]=mu[i].getValue();
			}
			
			double[] kdoubleForm = new double[3];
			for (int i=0;i<3;i++) {				
				kdoubleForm[i]=cern.jet.math.Functions.exp.apply(logk[i].getValue());
			}
			
			// TODO: add Dirichlet distribution...
			
			double[][] forPisdoubleform= new double[3][3];
			for (int i=0;i<3;i++) {
				for (int j=0;j<3;j++) {
					forPisdoubleform[i][j] = forPis[i][j].getValue();
				}
			}
			
			double[][] pisdoubleform = new double[3][3];
			for (int i=0;i<3;i++) {
				pisdoubleform[i] = Util.toDirichletNonDegenerate(forPisdoubleform[i]);
			}			

			// TODO: add update to migration model instead of reconstructing...
			TransitionModel migrationBaseModel = new ConstantTransitionBaseModel(ratesdoubleForm);
			TransitionModel[] codonModel = new TransitionModel[3];
			for (int i=0; i<3; i++) {
				codonModel[i]=new ConstantTransitionBaseModel(mudoubleForm[i],kdoubleForm[i],pisdoubleform[i][0],pisdoubleform[i][1],pisdoubleform[i][2]);
			}
			
			LikelihoodTree workingCopy;
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1)
					workingCopy = data.getTrees().get(i).get((int)treeIndices[i].getValue()).copy(); 
				else
					workingCopy = data.getTrees().get(i).get(0).copy();
				workingCopy.setMigrationModel(migrationBaseModel);
				workingCopy.setCodonModel(codonModel);
				logP+=config.treeWeights[i]*workingCopy.logLikelihood();
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
