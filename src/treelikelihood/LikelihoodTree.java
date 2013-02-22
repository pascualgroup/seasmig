package treelikelihood;

public interface LikelihoodTree {
	public void setLikelihoodModel(Object likelihoodModel);
	public double logLikelihood();
	String print();
	public LikelihoodTree copy(); 
	public int getNumLocations();
}
