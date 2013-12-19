package seasmig.migrationmain;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Serializable;
import java.util.logging.Level;

import com.google.gson.Gson;

@SuppressWarnings("serial")
public class Config implements Serializable
{
	// OPTION DEFINITIONS
	
	public static enum MigrationModelType { 
		CONSTANT, 
		TWO_CONSTANT_SEASONS, 
		SINUSOIDAL, // TODO:
		N_CONSTANT_SEASONS,  
		EPOCHAL
	}
	
	public static enum SeqModelType { 
		HKY_3CP // HKY model for 3 codon positions
	}
	
	// CONSTANT
	public static enum NoSeasonalityParameterization {  
		VARIABLE_SELECTION,		
		ALL,
	}
	
	// TWO CONSTANT SEASONS
	public static enum TwoConstantSeasonsParameterization 	{  
	   RATES12_PARAMETERZIATION,    // Separate rates for the two parts of the years rates1[i][j], rates2[i][j]
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
	   VARIABLE_SELECTION_DIFF		// (1-DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season1 
	   								// (1+DiffMult[i][j]*Indicator[i][j])*rate[i][j] is used for season2};	   
	}
	
	public static enum TwoConstantSeasonsPhase	{ 
		FIXED_PHASE_FIXED_LENGTH, 
		FREE_PHASE_FREE_LENGTH, 
		FIXED_PHASE_FREE_LENGTH, 
		FREE_PHASE_FIXED_LENGTH 
	}
	
	// EPOCHAL
	public static enum EpochParameterization {  
		EPOCHAL,
		EPOCHAL_VS, 
		EPOCHAL_FREE_TIMES, 	   
		EPOCHAL_FREE_TIMES_VS
	}

	// N_CONSTANT_SEASONS
	public static enum NConstantSeasonsParameterization {  
		ALL,
		VARIABLE_SELECTION
	}
	
	// STATE MODEL
	public static enum StateModel { NONE, BROWNIAN, BROWNIAN_SEASONAL };   // TODO: 

	// RECONSTRUCTION
	public static enum StateReconstructionAndTreeOutput { 
		NONE, // don't output trees
		PROBS, // output conditional probability of internal node states based on children
		ASR, // output realization of internal node states // TODO: test 
		STOCHASTIC_MAPPING, // output realization of branches and internal node states // TODO: test  
		EXACT_MAPPING_EXPECTATION  // output realization of branches and internal node states // TODO: 
	};

	// SETTINGS
	
	public Long randomSeed;
	
	// IO RELATED PARAMETERS
	public String sampleFilename = "samples.jsons";
	public String swapStatsFilename = "swap_stats.txt";
	public Level logLevel = Level.INFO;
	public String checkpointFilename = "checkpoint.bin";
	public String priorLikelihoodFilename = "prior_likelihood.txt";
	public String mlFilename = "ml.txt";
	public long thin = 1;
	
	// MCMC RELATED PARAMETERS
	public long burnIn = 10; 	// in iterations	
	public long iterationCount = 100000000L;	
	public long tuneEvery = 500000; 
	public long tuneFor = 1;
	public long mlthin = 10;	
	public long initialHistoryCount = 5;
	public int chainCount = 4;
	public double heatPower = 3.0;
	public long swapInterval = 1;	
	public double targetAcceptanceRate = 0.25;
	
	// RESTORE FROM CHECKPOINT 
	public long checkpointEvery = 1000;
	public boolean restoreFromDatabase = false;

	// SEQ EVOLUTION MODEL RELATED PARAMETERS
	
	
	// MIGRATION MODEL RELATED PARAMETERS
	public MigrationModelType migrationModelType = MigrationModelType.TWO_CONSTANT_SEASONS;
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
	
	// STOCHASTIC MAPPING OUTPUT
	public boolean asrTrees = true;  
	public boolean smTrees = true;  // TODO: add alternative with single child branches instead of &map 
	public boolean smTransitions = true; 
	public boolean smTipDwellings = true; 
	public boolean smLineages = true;
	public boolean smDescendants = true;
	public boolean smTrunkStats = true; 
	
	// MODEL DATA RELATED PARAMETERS
	public String[] locationFilenames ={"regionsHA.txt","regionsNA.txt"}; // null if locations are loaded from tree
	public String[] treeFilenames = {"beastInputHA.trees","beastInputNA.trees"}; // null for test generated data 
	public double[] treeWeights = {0.5,0.5};
	public int numTreesFromTail = 10; // at most number of trees to read from tree file's tail
	public int numLocations = 8; // needs to be specified if locations are loaded from trees....
	
	// VARIABLE SELECTION
	public double rateIndicatorPrior = 0.5; // Prior for including any migration between two locations 
	
	// ANCESTRAL STATE RECONSTRUCTIOn
	public StateReconstructionAndTreeOutput stateReconstructionAndTreeOutput = StateReconstructionAndTreeOutput.STOCHASTIC_MAPPING;
	
	// TIME CALIBRATION (THIS PARAMETER IS ABSOLUTLY CRUCIAL) 
	public double[] lastTipTime = {2012.74, 2011.7}; // time of most recent tip on tree, used to calibrate all tree nodes 

	// STOCHASTIC MAPPING OF TRUNK
	public double presentDayTipInterval = 0.25; // width of time interval of recent considered to have "survived" for trunk designation purpose  
	public double timeToDesignateTrunk = 2.0; // time back from present day tip ancestry designated as trunk 
	public int maxSMBranchRetries = 20000; // maximum number of retries for stochastically mapping a single branch

	public String[] alignmentFilenames = null;//{"ha.fasta","na.fasta"};

	public SeqModelType seqModelType = SeqModelType.HKY_3CP;
	
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
