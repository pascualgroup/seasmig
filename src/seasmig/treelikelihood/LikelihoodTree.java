package seasmig.treelikelihood;

import java.io.Serializable;

import seasmig.util.AltTreeOutput;

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
	public String smTrunkStats(double presentDayTipInterval, double timeToDesignateTrunk);
	
	// Alternative tree output nodes and branches following stochastic mapping
	public AltTreeOutput smAlternativeTreeOutput();
	
	// Stochastic Mapping & State Reconstruction of Sequence and Location
	public String seqMutationStats(int maxBranchRetries) throws Exception;
	public String seqMigrationsSeqOutput() throws Exception;
	String pies();
	public Double seqLikelihood();
	public Double locLikelihood();
	
	// Cleanup after reconstruction
	public void clearInternalNodes();
}
