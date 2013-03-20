package seasmig;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.logging.Level;
import com.google.gson.Gson;

@SuppressWarnings("serial")
public class Config implements Serializable
{
	public static enum Seasonality { NONE, TWO_CONSTANT_SEASONS, SINUSOIDAL};  //TODO: IMPLEMENT ADD CONTINUOUS SEASONAL MODEL //TODO: IMPLEMENT SINUSOIDAL 
	public static enum TwoConstantSeasonsParameterization 
	// TODO: check from to corresponds correctly to from to...
	{  DIRECT_ALL_FREE, 	  // Separate rates for two matrices
	   DIFF_PARAMETERIZATION, // (1-DiffMult[i][j])*rate[i][j] is used for season1 
	   						  // (1+DiffMult[i][j])*rate[i][j] is used for season2
	   FIX_SOME,	   	  	  // Same as above, but only a partial set (define by fixSome)  
	   						  // of DiffMults is used   
	   FIX_FROM_DIFF,  		  // Same as above, but each rows uses the same DiffMult:
	   						  // (1+-DiffMult[i])*rate[i][j] is used for season1/2 
	   FIX_TO_DIFF,    	      // Same as above, but each col uses the same DiffMult
	   						  // (1+-DiffMult[j])*rate[i][j] is used for season1/2
	   FIX_SOME_FROM,		  // Same as above, but only a partial set (defined by fixSomeFromTo
	   						  // of DiffMults is used 
	   FIX_SOME_TO,
	   VARIABLE_SELECTION,     // (1-DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season1 
		  					  // (1+DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season2};
	   VARIABLE_SELECTION_TO,
	   VARIABLE_SELECTION_FROM	   
	}
	public static enum TwoConstantSeasonsPhase { FIXED_PHASE_FIXED_LENGTH, FREE_PHASE_FREE_LENGTH, FIXED_PHASE_FREE_LENGTH, FREE_PHASE_FIXED_LENGTH };
	public static enum StateModel { NONE, BROWNIAN, BROWNIAN_SEASONAL };   // TODO: IMPLEMENT THIS... 
	
	public Long randomSeed;
	
	// IO RELATED PARAMETERS
	public String sampleFilename = "samples.jsons";
	public String swapStatsFilename = "swap_stats.txt";
	public Level logLevel = Level.INFO;
	public String checkpointFilename = "checkpoint.bin";
	public String priorLikelihoodFilename = "prior_likelihood.txt";
	public String mlFilename = "ml.txt";
	public long checkpointEvery = 150;
	public long thin = 50;
	
	// MCMC RELATED PARAMETERS
	public long burnIn = 200; 	// in iterations	
	public long iterationCount = 100000000L;
	
	public long tuneEvery = 50; 
	public long tuneFor = 200;
	public long mlthin = 10;
	
	public long initialHistoryCount = 50;

	public int chainCount = 4;
	public double heatPower = 3.0;
	
	public double targetAcceptanceRate = 0.25;

	// MODEL RELATED PARAMETERS
	public Seasonality migrationSeasonality = Seasonality.TWO_CONSTANT_SEASONS;
	public TwoConstantSeasonsParameterization twoSeasonParameterization = TwoConstantSeasonsParameterization.VARIABLE_SELECTION; 
	public TwoConstantSeasonsPhase twoSeasonPhase = TwoConstantSeasonsPhase.FREE_PHASE_FREE_LENGTH; // FREE LENGTH ONLY IMPLEMENTED FOR VARIABLE SELECTION...
	public double fixedPhase = 0.12;
	public int[][] fixSome = {{0,1},{1,0},{2,0},{2,1}};
	public int[] fixSomeFromTo = {0,1};
	
	public StateModel stateModel = StateModel.NONE; // TODO: IMPLEMENT THIS...
	
	// MODEL DATA RELATED PARAMETERS
	// TODO: add statese & combine files 
	public String stateFilename =null; // null if states are loaded from tree or non-existent 
	public String locationFilename ="regions.txt"; // null if locations are loaded from tree
	public String treeFilename = "beastInput.trees"; // null for test generated data 
	public String locationAttributeNameInTree = "states"; // location attribute in jebl tree
	public String stateAttributeNameInTree = null; // state attribute in jebl tree
	public int numTreesFromTail = 10; // at most number of trees to read from tree file's tail
	public Integer numLocations = null; // needs to be specified if locations are loaded from trees....
		
	protected Config() {};
	
	// OUTPUT CONFIG TO FILE
	public void outputToFile(String outfilename, Gson gson) {
		try {
			File configOutputFile = new File("out.config.json");
			configOutputFile.delete();
			configOutputFile.createNewFile();
			PrintStream configOutputStream = new PrintStream(configOutputFile);
			configOutputStream.print(gson.toJson(this).toString());
			configOutputStream.close();
		} catch (IOException e) {
			System.err.println("Failed to generate output file out.config.json!");
			e.printStackTrace();
		}
	}	
	
}
