package seasmig;


import mc3kit.LogLevel;

public class Config
{
	enum Seasonality
	{
		NONE,
		TWO_MATRICES,
		SINUSOIDAL
	}
	
	public String sampleFilename = "samples.jsons";
	public String priorLikelihoodFilename = "prior_likelihood.txt";
	public String varStatsFilename = "var_stats.jsons";
	public String demcStatsFilename = "demc_stats.jsons";
	public boolean recordHeatedStats = false;
	public String swapStatsFilename = "swap_stats.txt";
	
	public Long randomSeed;
	
	public long thin = 10;
	
	public long tuneEvery = 100;
	public long tuneFor = 10000;
	
	public long initialHistoryCount = 100;
	public long recordHistoryAfter = 1000;
	
	public int chainCount = 1;
	public double heatPower = 3.0;
	
	public LogLevel logLevel = LogLevel.INFO;
	public String logFilename = "debug.log";
	
	public String dataFilename = "data.csv";
	
	public Seasonality seasonality = Seasonality.NONE;
	public int stateCount = 2;
}
