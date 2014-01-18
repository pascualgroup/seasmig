package seasmig.treelikelihood;

import java.io.Serializable;

public interface LikelihoodTree extends Serializable {
	public void setMigrationModel(Object migrationModel);
	public void setCodonModel(Object codonModel);
	public double logLikelihood();
	String print();
	public LikelihoodTree copy(); 
	public int getNumLocations();
	public String newickProbs();
	public double cachedLogLikelihood();
	
	// Stochastic Mapping & State Reconstruction of Migration
	public String newickAncestralStateReconstruction();
	public String newickStochasticMapping(int maxBranchRetries);
	public String smTransitions();	
	public String smTipDwellings();
	public String smLineages();
	public String smDescendants();
	String smTrunkStats(double presentDayTipInterval,
			double timeToDesignateTrunk);
	
	// Stochastic Mapping & State Reconstruction of Sequence and Location
	public String seqMutationStats(int maxBranchRetries) throws Exception;
	String pis();
	public Double seqLikelihood();
	public Double locLikelihood();
	
	// Cleanup after reconstruction
	public void clearInternalNodes();
	public Object getRoot(); // TODO: remove this
}
