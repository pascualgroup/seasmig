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
	enum Seasonality {	NONE, TWO_CONSTANT_SEASONS,TWO_CONSTANT_SEASONS_FIXED_PHASE, SINUSOIDAL};  //TODO: IMPLEMENT ADD CONTINUOUS SEASONAL MODEL //TODO: IMPLEMENT SINUSOIDAL 
	enum StateModel { NONE, BROWNIAN, BROWNIAN_SEASONAL };   // TODO: IMPLEMENT THIS... 
	
	public Long randomSeed;
	
	// IO RELATED PARAMETERS
	public String sampleFilename = "samples.jsons";
	public String swapStatsFilename = "swap_stats.txt";
	public Level logLevel = Level.INFO;
	public String checkpointFilename = "checkpoint.bin";
	public String priorLikelihoodFilename = "prior_likelihood.txt";	

	// MCMC RELATED PARAMETERS
	public long burnIn = 2000; 	// in iterations	
	public long iterationCount = 100000000L;
	
	public long tuneEvery = 1000; 
	public long tuneFor = 10000;
	public long thin = 50;
	
	public long initialHistoryCount = 20000;
	public long checkpointEvery = 100;

	public int chainCount = 16;
	public double heatPower = 2.0;
	
	public double targetAcceptanceRate = 0.25;

	// MODEL RELATED PARAMETERS
	public Seasonality migrationSeasonality = Seasonality.TWO_CONSTANT_SEASONS;
	public StateModel stateModel = StateModel.NONE; // TODO: IMPLEMENT THIS...
	public double fixedPhase = 0.1;
	
	// MODEL DATA RELATED PARAMETERS
	// TODO: add statese & combine files 
	public String stateFilename =null; // null if states are loaded from tree or non-existent 
	public String locationFilename ="regions.txt"; // null if locations are loaded from tree
	public String treeFilename = "beastInput.trees"; // null for test generated data 
	public String locationAttributeNameInTree = "states"; // location attribute in jebl tree
	public String stateAttributeNameInTree = null; // state attribute in jebl tree
	public int numTreesFromTail = 100; // at most number of trees to read from tree file's tail
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
