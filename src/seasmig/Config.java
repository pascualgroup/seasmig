package seasmig;


import mc3kit.LogLevel;

public class Config
{
	enum Seasonality {	NONE, TWO_CONSTANT_SEASONS,	SINUSOIDAL }; //TODO: DEBUG SINUSODIAL
																  //TODO: ADD CONTINUOUS SEASONAL MODEL
	enum StateModel { NONE, BROWNIAN, BROWNIAN_SEASONAL };   // TODO: IMPLEMENT THIS... 
	
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
	public long thin = 50;

	public long tuneEvery = 10;
	public long tuneFor = 100;
	
	public long initialHistoryCount = 100;
	public long recordHistoryAfter = 1000;
	
	public int chainCount = 4;
	public double heatPower = 3.0;

	// MODEL RELATED PARAMETERS
	public Seasonality migrationSeasonality = Seasonality.NONE;
	public StateModel stateModel = StateModel.NONE; // TODO: this...
	public int numLocations = 3;  // TODO: add support of one location....
							      // TODO: add as an attribute loaded with attribute loader...
	
	// MODEL DATA RELATED PARAMETERS
	// TODO: add statese & combine files 
	public String stateFilename =null; // null if states are loaded from tree or non-existent 
	public String locationFilename ="locations.txt"; // null if locations are loaded from tree
	public String treeFilename = "beastInput.trees"; // null for test generated data 
	public String locationAttributeNameInTree = "states"; // location attribute in jebl tree
	public String stateAttributeNameInTree = null; // location attribute in jebl tree
	public int numTreesFromTail = 99; // at most number of trees to read from tree file's tail
	
	// TEST RELATED PARAMETERS
	public int numTestTrees = 50;
	public int numTestTips = 500;
	public int numTestRepeats = 5; 
	public double disturbanceScale = 0.3;
	
}
