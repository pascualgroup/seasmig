package seasmig.models;

import java.io.Serializable;

import seasmig.migrationmain.Config;
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

	public class TreesLikelihoodVariableOutputObject implements Serializable {		
		protected TreesLikelihoodVariableOutputObject() {};
		String[] probTrees = null;
		public String[] asrTrees = null;
		public String[] smTrees = null;
		public String[] smTransitions = null;
		public String[] smTipDwellings = null;
		public String[] smLineages = null;	
		public String[] smDescendants = null;
		public String[] smTrunkStats = null;
		public String[] smSeqMigStats = null;
		public String[] pis = null;
		public String[] seqMutationStats = null;
	}
	

	@Override
	public Object makeOutputObject() {		
		TreesLikelihoodVariableOutputObject outputObject = new TreesLikelihoodVariableOutputObject();
		switch (config.stateReconstructionAndTreeOutput) {
		case NONE : {
			return outputObject;		
		}
		case PROBS: {
			String[] returnValue=new String[trees.length];
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickProbs()));
			}
			outputObject.probTrees=returnValue;
			return outputObject;
		}
		case ASR: {
			String[] returnValue=new String[trees.length];
			for (int i=0;i<trees.length;i++) {
				Model model = this.getModel();		 	
				// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
				String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
				returnValue[i]=(header + (trees[i].newickAncestralStateReconstruction()));			
			}
			outputObject.asrTrees=returnValue;
			return outputObject;
		}
		case STOCHASTIC_MAPPING: {					
			if (config.smTrees) {				
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					Model model = this.getModel();		 	
					// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
					String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
					returnValue[i]=(header + (trees[i].newickStochasticMapping(config.maxSMBranchRetries)));
				}				
				outputObject.smTrees=returnValue;
			}		
			if (config.asrTrees) {				
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					Model model = this.getModel();		 	
					// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
					String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
					returnValue[i]=(header + (trees[i].newickAncestralStateReconstruction()));
				}				
				outputObject.asrTrees=returnValue;
			}
			if (config.smTransitions) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smTransitions();
				}				
				outputObject.smTransitions=returnValue;				
			}
			if (config.smTipDwellings) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smTipDwellings();
				}				
				outputObject.smTipDwellings=returnValue;				
			}
			if (config.smLineages) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smLineages();
				}				
				outputObject.smLineages=returnValue;				
			}
			if (config.smDescendants) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smDescendants();
				}				
				outputObject.smDescendants=returnValue;				
			}
			if (config.smTrunkStats) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smTrunkStats(config.presentDayTipInterval, config.timeToDesignateTrunk);
				}				
				outputObject.smTrunkStats=returnValue;				
			}
			return outputObject;
		}
		case SEQ_STOCHASTIC_MAPPING: {	
			if (config.smTrees) {				
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					Model model = this.getModel();		 	
					// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
					String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
					returnValue[i]=(header + (trees[i].newickStochasticMapping(config.maxSMBranchRetries)));
				}				
				outputObject.smTrees=returnValue;
			}		
			if (config.asrTrees) {				
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					Model model = this.getModel();		 	
					// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
					String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
					returnValue[i]=(header + (trees[i].newickAncestralStateReconstruction()));
				}				
				outputObject.asrTrees=returnValue;
			}
			if (config.smTransitions) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smTransitions();
				}				
				outputObject.smTransitions=returnValue;				
			}
			if (config.smTipDwellings) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {
					returnValue[i]=trees[i].smTipDwellings();
				}				
				outputObject.smTipDwellings=returnValue;				
			}
			if (config.smLineages) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smLineages();
				}				
				outputObject.smLineages=returnValue;				
			}
			if (config.smDescendants) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smDescendants();
				}				
				outputObject.smDescendants=returnValue;				
			}
			if (config.smTrunkStats) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					returnValue[i]=trees[i].smTrunkStats(config.presentDayTipInterval, config.timeToDesignateTrunk);
				}				
				outputObject.smTrunkStats=returnValue;				
			}
			
			/////////////////////////////////////
			if (config.seqMutationStats) {
				String[] returnValue=new String[trees.length];
				for (int i=0;i<trees.length;i++) {	 			
					try {
						returnValue[i]=trees[i].seqMutationStats(config.maxSMBranchRetries);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}				
				outputObject.seqMutationStats=returnValue;				
			}
			
			String[] returnValue=new String[trees.length];
			for (int i=0;i<trees.length;i++) {	 			
				returnValue[i]=trees[i].pis();
			}				
			outputObject.pis=returnValue;											

			return outputObject;
		}
		default: 
			break;	
		}	

		return null;
	}
}
