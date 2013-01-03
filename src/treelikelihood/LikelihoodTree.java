package treelikelihood;

public interface LikelihoodTree {
	public void setLikelihoodModel(Object likelihoodModel);
	public double logLikelihood();
	String print(); 
}
