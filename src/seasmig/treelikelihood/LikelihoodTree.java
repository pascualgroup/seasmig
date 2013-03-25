package seasmig.treelikelihood;

import java.io.Serializable;

public interface LikelihoodTree extends Serializable {
	public void setLikelihoodModel(Object likelihoodModel);
	public double logLikelihood();
	String print();
	public LikelihoodTree copy(); 
	public int getNumLocations();
}
