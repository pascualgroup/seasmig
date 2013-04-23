package seasmig.models;

import seasmig.treelikelihood.LikelihoodTree;
import mc3kit.Model;
import mc3kit.Variable;

@SuppressWarnings("serial")
public class TreesLikelihoodVariable extends Variable {
	protected double oldLogP;
	protected LikelihoodTree[] trees;
	
	public TreesLikelihoodVariable(Model m,
			String name, boolean observable, int numTrees) {
		super(m,name,observable);
		trees = new LikelihoodTree[numTrees];
	}

	@Override
	public Object makeOutputObject() {
		String[] returnValue=new String[trees.length];
		for (int i=0;i<trees.length;i++) {
			Model model = this.getModel();			
			// "tree STATE_50850000 [&lnP=-34291.617355973016,posterior=-34291.617355973016] = [&R]"
			String header = "tree STATE_" + getChain().getIterationCount() + " [&lnP=" + (model.getLogPrior()+trees[i].cachedLogLikelihood()) + "]" +  " = [&R] ";			
			returnValue[i]=(header + (trees[i].newick()+String.format("\n")));
		}			
		return returnValue;
	}
}
