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
		EPOCHAL,
		CONSTANT_AS_INPUT
	}
	
	public static enum SeqModelType { 
		NONE, // Sequences are not included in likelihood
		HKY_3CP, // HKY model for 3 codon positions
		HKY_3CP_AS_INPUT 
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
		ASR, // output realization of internal node states // TODO: better test
		STOCHASTIC_MAPPING, // output realization of branches and internal node states // TODO: better test  
		EXACT_MAPPING_EXPECTATION,  // output realization of branches and internal node states // TODO: 
		SEQ_STOCHASTIC_MAPPING 
	};

	// SETTINGS
	
	public Long randomSeed;
	
	// IO RELATED PARAMETERS
	public String sampleFilename = "F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\samples.jsons";
	public String swapStatsFilename = "F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\swap_stats.txt";
	public Level logLevel = Level.INFO;
	public String checkpointFilename = "F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\checkpoint.bin";
	public String priorLikelihoodFilename = "F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\prior_likelihood.txt";
	public String mlFilename = "F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\ml.txt";
	public long thin = 20;
	
	// MCMC RELATED PARAMETERS
	public long burnIn = 200; 	// in iterations	
	public long iterationCount = 100000000L;	
	public long tuneEvery = 50; 
	public long tuneFor = 200;
	public long mlthin = 5;	
	public long initialHistoryCount = 5;
	public int chainCount = 1;
	public double heatPower = 1.0;
	public long swapInterval = 1;	
	public double targetAcceptanceRate = 0.25;
	
	// RESTORE FROM CHECKPOINT 
	public long checkpointEvery = 1000;
	public boolean restoreFromDatabase = false;

	// SEQ EVOLUTION MODEL RELATED PARAMETERS
	
	
	// MIGRATION MODEL RELATED PARAMETERS
	public MigrationModelType migrationModelType = MigrationModelType.CONSTANT_AS_INPUT;
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
	
	// trees
	public boolean asrTrees = true;  
	public boolean smTrees = true;  
	public boolean smAlternativeTreeOutput = true; // nodes and branches
	
	// events
	public boolean smTransitions = true; 
	public boolean smTipDwellings = true; 
	public boolean smLineages = false;
	public boolean smDescendants = false;
	public boolean smTrunkStats = false;
	public boolean smMigrationNodeNumTipAndSequenceData = true; // output sequence and node data for each migration event
	public boolean seqMutationStats = true;
	public boolean seqMutationsStatsCodonOutput = true;
	public boolean seqMutationsStatsSeqOutput = true;
	public double seqStochasticMappingStartTime = 1970; 
	
	// MODEL DATA RELATED PARAMETERS
	public String[] locationFilenames ={"F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\accession_locid.txt"}; // null if locations are loaded from tree
	public String[] treeFilenames = {"F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\newTrait_tree_with_trait.trees"}; // null for test generated data 
	public double[] treeWeights = {1.0};
	public int numTreesFromTail = 50; // at most number of trees to read from tree file's tail
	public int numLocations = 7; // needs to be specified if locations are loaded from trees....
	public boolean sampleTreesSequentially = true;
	
	// ANCESTRAL STATE RECONSTRUCTIOn
	public StateReconstructionAndTreeOutput stateReconstructionAndTreeOutput = StateReconstructionAndTreeOutput.SEQ_STOCHASTIC_MAPPING ;
	
	// TIME CALIBRATION (THIS PARAMETER IS ABSOLUTLY CRUCIAL) 
	public double[] lastTipTime = {2012.0}; // time of most recent tip on tree, used to calibrate all tree nodes 

	// STOCHASTIC MAPPING OF TRUNK
	public double presentDayTipInterval = 2.25; // width of time interval of recent tips considered to have "survived" for trunk designation purpose  
	public double timeToDesignateTrunk = 4.0; // time back from present day tip ancestry designated as trunk 
	public int maxSMBranchRetries = 200000; // maximum number of retries for stochastically mapping a single branch

	public String[] alignmentFilenames = {"F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\all_epitopes.txt"};

	public SeqModelType seqModelType = SeqModelType.HKY_3CP_AS_INPUT;

	public double verificationTolerance = 3; // TODO: check why verification fails ?ed? with seq data, and is it just a convergence things.

	public String[] migrationModelFilenames = {"F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\migration_models.txt"};
	public String[] codonModelFilenames = {"F:\\Daniel\\Dropbox\\SharedFolderBobDanielRV\\who_and_us\\for_testing_seasmig_mapping_without_mcmc\\codon_models.txt"};


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
