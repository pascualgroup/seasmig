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
	public static enum ModelType { CONSTANT, TWO_CONSTANT_SEASONS, SINUSOIDAL, N_CONSTANT_SEASONS, EPOCHAL};
	
	//TODO: IMPLEMENT ADD CONTINUOUS SEASONAL MODEL //TODO: IMPLEMENT SINUSOIDAL
	//TODO: GET N_CONSTANT_SEASONS and N_CONSTANT_SEASONS_VAR_SELECT COMBINED FOR THIS ENUM
	
	public static enum TwoConstantSeasonsParameterization 
	{  RATES12_PARAMETERZIATION,    // Separate rates for the two parts of the years rates1[i][j], rates2[i][j]
	   RATES12_VARIABLE_SELECTION,  // rate1[i][j]*indicators1[i][j], rate2[i][j]*indicators2[i][j]
	   DIFF_PARAMETERIZATION,       // (1-DiffMult[i][j])*rate[i][j] is used for season1 
	   						   		// (1+DiffMult[i][j])*rate[i][j] is used for season2
	   FIX_FROM_DIFF,  		   		// Same as above, but each rows uses the same DiffMult:
	   						   		// (1+-DiffMult[i])*rate[i][j] is used for season1/2 
	   FIX_TO_DIFF,    	       		// Same as above, but each col uses the same DiffMult
	   						   		// (1+-DiffMult[j])*rate[i][j] is used for season1/2
	   FIX_FROM_TO_DIFF,       		// Both rows and columns...
	   VARIABLE_SELECTION,     		// (1-DiffMult[i][j]*Indicator[i][j])*rate[i][j]*rateIndicator[i][j] is used for season1 
		  					   		// (1+DiffMult[i][j]*Indicator[i][j])*rate[i][j]*rateIndicator[i][j] is used for season2};
	   VARIABLE_SELECTION_DIFF,		// (1-DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season1 
	   								// (1+DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season2};
	   VARIABLE_SELECTION_GTR	   
	}
	
	public static enum NoSeasonalityParameterization 
	{  VARIABLE_SELECTION,		
	   ALL,
	}
	
	public static enum EpochParameterization 
	{  EPOCHAL,
	   EPOCHAL_VS, 
	   EPOCHAL_FREE_TIMES, 	   
	   EPOCHAL_FREE_TIMES_VS
	}
	
	public static enum NConstantSeasonsParameterization 
	{  ALL,
	   VARIABLE_SELECTION
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
	public long thin = 5;
	
	// MCMC RELATED PARAMETERS
	public long burnIn = 2000; 	// in iterations	
	public long iterationCount = 100000000L;	
	public long tuneEvery = 500; 
	public long tuneFor = 5000;
	public long mlthin = 5;	
	public long initialHistoryCount = 50;
	public int chainCount = 4;
	public double heatPower = 3.0;
	public long swapInterval = 10;	
	public double targetAcceptanceRate = 0.25;
	
	// RESTORE FROM CHECKPOINT 
	public boolean verificationStep = false; // When false allows for "hot swaps" of tree files, when resuming from checkpoint.bin
	public long checkpointEvery = 10000;
	public boolean restoreFromDatabase = false;

	// SEASONAL MODEL RELATED PARAMETERS
	public ModelType modelType = ModelType.TWO_CONSTANT_SEASONS;
	public TwoConstantSeasonsParameterization twoSeasonParameterization = TwoConstantSeasonsParameterization.VARIABLE_SELECTION_DIFF;
	public TwoConstantSeasonsPhase twoSeasonPhase = TwoConstantSeasonsPhase.FREE_PHASE_FIXED_LENGTH; // FREE LENGTH ONLY IMPLEMENTED FOR VARIABLE SELECTION TWO SEASONS...
	public double fixedPhase = 0.3;
	public double minSeasonLength = 0.3333333;
	public boolean fixRate = false;
	
	// N_CONSTANT_SEASONS 
	public int nSeasonalParts=4; 
	public NConstantSeasonsParameterization nConstantSeasonsParameterization = NConstantSeasonsParameterization.ALL;
	
	// NON-SEASONAL 
	public NoSeasonalityParameterization noSeasonalityParameterization = NoSeasonalityParameterization.ALL;
	
	// EPOCHAL 
	public EpochParameterization epochParameterization = EpochParameterization.EPOCHAL;
	public int nEpochs = 2;
	public double[] epochTimes = {2006.0};
	public double minEpochTime = 1985;
	public double maxEpochTime = 2011;
	
	// GENERAL PARAMETERS
	public double veryLongTime = 1000;		
	public StateModel stateModel = StateModel.NONE; // TODO: IMPLEMENT THIS...
	
	// MODEL DATA RELATED PARAMETERS
	public String stateFilename =null; // TODO: null if states are loaded from tree or non-existent 
	public String[] locationFilenames ={"regionsHA.txt","regionsNA.txt"}; // null if locations are loaded from tree
	public String[] treeFilenames = {"beastInputHA.trees","beastInputNA.trees"}; // null for test generated data 
	public double[] treeWeights = {0.5,0.5};
	public String locationAttributeNameInTree = "states"; // location attribute in jebl tree
	public String stateAttributeNameInTree = null; // state attribute in jebl tree
	public int numTreesFromTail = 10; // at most number of trees to read from tree file's tail
	public int numLocations = 9; // needs to be specified if locations are loaded from trees....
	
	// TIME CALIBRATION (THIS PARAMETER IS ABSOLUTLY CRUCIAL) 
	public double lastTipTime = 2012.74; // time of most recent tip on tree, used to calibrate all tree nodes 
			
	public Config() {};
	
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
