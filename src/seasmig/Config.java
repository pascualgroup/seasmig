package seasmig;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.google.gson.Gson;

import mc3kit.LogLevel;

public class Config
{
	enum Seasonality {	NONE, TWO_CONSTANT_SEASONS,TWO_CONSTANT_SEASONS_FIXED_PHASE,SINUSOIDAL }; //TODO: DEBUG SINUSODIAL //TODO: ADD CONTINUOUS SEASONAL MODEL
	enum StateModel { NONE, BROWNIAN, BROWNIAN_SEASONAL };   // TODO: IMPLEMENT THIS... 
	
	enum RunMode {	NORMAL,	TEST_USING_GENERATED_TREES, TEST_USING_INPUT_TREES, TEST_MODEL_DEGENERACY}

	public boolean useAllTreesForSingleLikelihoodCalculation = true;
	
	public Long randomSeed;
	
	public RunMode runMode = RunMode.NORMAL;
	
	// LOG RELATED PARAMETERS
	public String sampleFilename = "samples.jsons";
	public String priorLikelihoodFilename = "prior_likelihood.txt";
	public String varStatsFilename = "var_stats.jsons";
	public String demcStatsFilename = "demc_stats.jsons";
	public boolean recordHeatedStats = false;
	public String swapStatsFilename = "swap_stats.txt";
	public LogLevel logLevel = LogLevel.ERROR;
	public String logFilename = "debug.log";

	// MCMC & LOG RELATED PARAMETERS
	// in iterations
	public long thin = 50;

	public long tuneEvery = 50; 
	public long tuneFor = 100000; 
	
	public long initialHistoryCount = 20000;
	public long recordHistoryAfter = 40000;
	
	public int chainCount = 16;
	public double heatPower = 3.0;
	
	// DISPLAY RELATED PARAMTERS
	public int printEveryNStates = 100000;

	// MODEL RELATED PARAMETERS
	public Seasonality migrationSeasonality = Seasonality.TWO_CONSTANT_SEASONS_FIXED_PHASE;
	public StateModel stateModel = StateModel.NONE; // TODO: this...
	public int numLocations = 3;  // TODO: add support of one location....
							      // TODO: add as an attribute loaded with attribute loader...
	public double fixedPhase = 0.1;
	
	// MODEL DATA RELATED PARAMETERS
	// TODO: add statese & combine files 
	public String stateFilename =null; // null if states are loaded from tree or non-existent 
	public String locationFilename ="regions.txt"; // null if locations are loaded from tree
	public String treeFilename = "beastInput.trees"; // null for test generated data 
	public String locationAttributeNameInTree = "states"; // location attribute in jebl tree
	public String stateAttributeNameInTree = null; // state attribute in jebl tree
	public int numTreesFromTail = 2; // at most number of trees to read from tree file's tail
	
	// TEST RELATED PARAMETERS	
	// TEST1+TEST2
	public int numTestTrees = 2;

	// TEST1
	public int numTestTips = 1200;
	public int numTestRepeats = 5; 
	public double disturbanceScale = 0.3;

	
		
	// OUTPUT CONFIG TO FILE
	public void outputToFile(String outfilename, Gson gson) throws IOException {
		File configOutputFile = new File("out.config.json");
		configOutputFile.delete();
		configOutputFile.createNewFile();
		PrintStream configOutputStream = new PrintStream(configOutputFile);
		configOutputStream.print(gson.toJson(this).toString());
		configOutputStream.close();
	}

	
	
}
