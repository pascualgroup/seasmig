package seasmig.treelikelihood;

import java.io.Serializable;

import seasmig.util.AltTreeOutput;

public interface LikelihoodTree extends Serializable {
	
	// SET MODELS
	public void setMigrationModel(Object migrationModel);
	public void setCodonModel(Object codonModel);
	
	// GET LIKELIHOODS
	public double logLikelihood(); // doesn't recalculate likelihood if already calculated after copy
	public Double seqLikelihood();
	public Double locLikelihood();

	// MAKE CLEAN WORKING COPY
	public LikelihoodTree copy();

	// GET OUTPUT
	public int getNumLocations();
	public String newickProbs();

	
	// Stochastic Mapping & State Reconstruction of Migration
	public String newickASR();
	public String newickSM(int maxBranchRetries);
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
	
	// Cleanup after reconstruction
	public void clearInternalNodes();
	
}
