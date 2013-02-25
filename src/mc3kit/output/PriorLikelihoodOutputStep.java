package mc3kit.output;

import java.util.*;
import java.io.*;

import mc3kit.*;
import mc3kit.util.*;

@SuppressWarnings("serial")
public class PriorLikelihoodOutputStep implements Step, Serializable
{
	private String filename;
	private long thin;
	private int chainCount;
	
	transient PrintWriter writer;
	transient Collector<LogPriorLikelihoodValue> collector;
//	Map<Long, ValueCollector> collectedValues;
	
	/*** METHODS ***/
	
	protected PriorLikelihoodOutputStep() { };
	
	public PriorLikelihoodOutputStep(String filename, long thin) throws FileNotFoundException {
	  this.filename = filename;
	  this.thin = thin;
	  
    writer = new PrintWriter(new FileOutputStream(filename, false));
    writer.println("iteration\tchainId\tlogPrior\tlogLikelihood");
    writer.flush();
	}
  
	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
	  in.defaultReadObject();
	  writer = new PrintWriter(new FileOutputStream(filename, true));
	}
	
	private synchronized List<LogPriorLikelihoodValue> takeValue(long iteration, int index, double logPrior, double logLikelihood) {
	  if(collector == null) {
	    collector = new Collector<LogPriorLikelihoodValue>(chainCount);
	  }
	  return collector.takeValue(iteration, index, new LogPriorLikelihoodValue(logPrior, logLikelihood));
	}
	
	@Override
	public List<Task> makeTasks(int chainCount) throws MC3KitException
	{
		this.chainCount = chainCount;
		List<Task> tasks = new ArrayList<Task>(chainCount);
		for(int i = 0; i < chainCount; i++)
		{
			tasks.add(new PLOutputTask(i));
		}
		return tasks;
	}
	
	/*** TASK CLASS ***/
	
	private class PLOutputTask implements Task
	{
		int chainId;
		
		private long iterationCount;
		
		@Override
		public int[] getChainIds()
		{
			return new int[] { chainId };
		}
		
		PLOutputTask(int chainId) throws MC3KitException
		{
			this.chainId = chainId;
		}
		
		@Override
		public void step(Chain[] chains) throws MC3KitException
		{
			iterationCount++;
			
			if(iterationCount % thin == 0)
			{
				assert(chains.length == 1);
				Chain chain = chains[0];
				Model model = chain.getModel();
				
				List<LogPriorLikelihoodValue> plValues = takeValue(iterationCount, chainId, model.getLogPrior(), model.getLogLikelihood());
				
				if(plValues != null)
				{
					for(int i = 0; i < chainCount; i++)
					{
						writer.printf("%d\t%d\t%.3f\t%.3f\n", iterationCount, i,
							plValues.get(i).logPrior,
							plValues.get(i).logLikelihood);
					}
          writer.flush();
				}
			}
		}
	}
	
	private class LogPriorLikelihoodValue implements Serializable
	{
		double logPrior;
		double logLikelihood;
		
		LogPriorLikelihoodValue(double logPrior, double logLikelihood)
		{
			this.logPrior = logPrior;
			this.logLikelihood = logLikelihood;
		}
	}
}
