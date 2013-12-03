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
		String[] returnValue=new String[trees.length];
		switch (config.stateReconstructionAndTreeOutput) {
		case NONE : {
			for (int i=0;i<trees.length;i++) {
				returnValue[i] = "";
			}					
		}
		break;
		case PROBS: {
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickProbs()+String.format("\n")));
			}
		}
		break;
		case ASR: {
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickAncestralStateReconstruction()+String.format("\n")));
			}
		}
		break;	
		case STOCHASTIC_MAPPING: {
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickStochasticMapping()+String.format("\n")));
			}
		}
		default: 
			break;	
		}	

		return returnValue;
	}
}
