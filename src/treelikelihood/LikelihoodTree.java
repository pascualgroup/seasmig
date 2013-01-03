package treelikelihood;

public interface LikelihoodTree {
	// TODO: think about likelihood model typing enforcement
	public void setLikelihoodModel(Object likelihoodModel);
	public double logLikelihood();
	String print();
	public LikelihoodTree workingCopy(); 
}
