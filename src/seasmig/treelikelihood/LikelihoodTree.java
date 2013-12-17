package seasmig.treelikelihood;

import java.io.Serializable;

public interface LikelihoodTree extends Serializable {
	public void setMigrationModel(Object migrationModel);
	public void setCodonModel(Object condonModel);
	public double logLikelihood();
	String print();
	public LikelihoodTree copy(); 
	public int getNumLocations();
	public String newickProbs();
	public double cachedLogLikelihood();
	public String newickAncestralStateReconstruction();
	public String newickStochasticMapping();
	public String smTransitions();	
	public String smTipDwellings();
	public String smLineages();
	public String smDescendants();
	String smTrunkStats(double presentDayTipInterval,
			double timeToDesignateTrunk);
}
