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
public class HKY_3CP_FromInput_NoMigrationSeasonality extends MigrationModel {

	// TODO: seperate into two variables for seqence and location likelihood...

	Config config;
	Data data;
	int numLocations;

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
	private ExponentialDistribution rateHyperPriorDist;
	private DoubleVariable rateHyperPrior;
	private double[] mudoubleForm;
	private double[] kdoubleForm;
	private double[][] pisdoubleForm = new double[3][];
	
	protected HKY_3CP_FromInput_NoMigrationSeasonality() { }

	public HKY_3CP_FromInput_NoMigrationSeasonality(Chain initialChain, Config config, Data data) throws MC3KitException
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
		
		mudoubleForm = config.HKY_3CP_Input_mu;
		kdoubleForm = config.HKY_3CP_Input_k;
		pisdoubleForm[0]= config.HKY_3CP_Input_pis0;
		pisdoubleForm[1]= config.HKY_3CP_Input_pis1;
		pisdoubleForm[2]= config.HKY_3CP_Input_pis2;
	
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

	private class LikelihoodVariable extends TreesLikelihoodVariable {
	
		LikelihoodVariable(HKY_3CP_FromInput_NoMigrationSeasonality m) throws MC3KitException {
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
			TransitionModel migrationBaseModel = new ConstantTransitionBaseModel(ratesdoubleForm);
			TransitionModel[] codonModel = new TransitionModel[3];
			for (int i=0; i<3; i++) {
				codonModel[i]=new ConstantTransitionBaseModel(mudoubleForm[i],kdoubleForm[i],pisdoubleForm[i][0],pisdoubleForm[i][1],pisdoubleForm[i][2]);
			}

			LikelihoodTree workingCopy;
			for (int i=0;i<nTrees.length;i++) {
				if (nTrees[i]>1)
					workingCopy = data.getTrees().get(i).get((int)treeIndices[i].getValue()).copy(); 
				else
					workingCopy = data.getTrees().get(i).get(0).copy();
				workingCopy.setMigrationModel(migrationBaseModel);
				workingCopy.setCodonModel(codonModel);
				logP+=config.treeWeights[i]*workingCopy.locLikelihood();
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
