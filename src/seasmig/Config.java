package seasmig;


import mc3kit.LogLevel;

public class Config
{
	enum Seasonality {	NONE, TWO_CONSTANT_SEASONS,	SINUSOIDAL }; //TODO: DEBUG SINUSODIAL
																  //TODO: ADD CONTINUOUS SEASONAL MODEL
	enum RunMode {	NORMAL,	TEST };
	
	public Long randomSeed;
	
	public RunMode runMode = RunMode.TEST;
	
	// LOG RELATED PARAMETERS
	public String sampleFilename = "samples.jsons";
	public String priorLikelihoodFilename = "prior_likelihood.txt";
	public String varStatsFilename = "var_stats.jsons";
	public String demcStatsFilename = "demc_stats.jsons";
	public boolean recordHeatedStats = false;
	public String swapStatsFilename = "swap_stats.txt";
	public LogLevel logLevel = LogLevel.INFO;
	public String logFilename = "debug.log";

	// MCMC & LOG RELATED PARAMETERS
	public long thin = 10;

	public long tuneEvery = 100;
	public long tuneFor = 10000;
	
	public long initialHistoryCount = 100;
	public long recordHistoryAfter = 1000;
	
	public int chainCount = 1;
	public double heatPower = 3.0;
	

	// MODEL RELATED PARAMETERS
	public Seasonality seasonality = Seasonality.NONE;
	public int locationCount = 4;  // TODO: add support of one location....
								   // TODO: add as an attribute loaded with attribute loader...
	
	// MODEL DATA RELATED PARAMETERS
	// TODO: add states & combine files 
	public String stateFilename =null; // null if states are loaded from tree or non-existent 
	public String locationFilename ="locations.txt"; // null if locations are loaded from tree
	public String treeFilename = "beastInput.trees"; // null for test generated data 
	public int numTreesFromTail = 50; // at most number of trees to read from tree file's tail
	
	// TEST RELATED PARAMETERS
	public Seasonality testTreesCreateSeasonality = Seasonality.SINUSOIDAL;
	public int numTestTrees = 10;
	public int numTestTips = 400; 
	
}
