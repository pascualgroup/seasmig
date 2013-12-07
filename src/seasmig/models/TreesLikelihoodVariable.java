package seasmig.models;

import seasmig.Config;
import seasmig.treelikelihood.LikelihoodTree;
import mc3kit.Model;
import mc3kit.Variable;

@SuppressWarnings("serial")
public class TreesLikelihoodVariable extends Variable {
	protected double oldLogP;
	protected boolean firstCall;
	protected LikelihoodTree[] trees;
	private Config config;

	public TreesLikelihoodVariable(Model m,
			String name, boolean observable, int numTrees, Config config_) {
		super(m,name,observable);
		trees = new LikelihoodTree[numTrees];
		firstCall = true;
		config = config_;
	}

	@Override
	public Object makeOutputObject() {		
		switch (config.stateReconstructionAndTreeOutput) {
		case NONE : {
			return null;		
		}
		case PROBS: {
			String[] returnValue=new String[trees.length];
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickProbs()));
			}
			return returnValue;
		}
		case ASR: {
			String[] returnValue=new String[trees.length];
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickAncestralStateReconstruction()));			
			}
			return returnValue;
		}
		case STOCHASTIC_MAPPING: {				
			class OutputObject  {
				@SuppressWarnings("unused")
				String[] smTrees = null;
				@SuppressWarnings("unused")
				String[] smTransitions = null;
				@SuppressWarnings("unused")
				String[] smDwelings = null;
				@SuppressWarnings("unused")
				String[] smLineages = null;
			}
			OutputObject outputObject = new OutputObject();
			if (config.smTrees) {				
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					Model model = this.getModel();		 	
					// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
					String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
					returnValue[i]=(header + (trees[i].newickStochasticMapping()));
				}				
				outputObject.smTrees=returnValue;
			}		
			if (config.smTransitions) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smTransitions();
				}				
				outputObject.smTransitions=returnValue;				
			}
			if (config.smDwellings) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smDwellings();
				}				
				outputObject.smDwelings=returnValue;				
			}
			if (config.smLineages) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smLineages();
				}				
				outputObject.smLineages=returnValue;				
			}
			return outputObject;
		}
		default: 
			break;	
		}	

		return null;
	}
}
