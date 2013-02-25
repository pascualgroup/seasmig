package mc3kit.monitoring;

import java.io.*;
import java.util.*;

import mc3kit.*;
import mc3kit.util.*;
import static mc3kit.util.Math.*;

@SuppressWarnings("serial")
public class MarginalLikelihoodStep implements Step {
  private String filename;
  private long burnIn;
  private long thin;
  
  transient PrintWriter writer;
  transient Collector<MargLikeValue> collector;
  
  protected MarginalLikelihoodStep() { }
  
  public MarginalLikelihoodStep(String filename, long burnIn, long thin) throws FileNotFoundException {
    this.filename = filename;
    this.burnIn = burnIn;
    this.thin = thin;
    
    writer = new PrintWriter(new FileOutputStream(filename, false));
    writer.println("iteration\tmarginalLikelihood\tmarginalLikelihoodAuto");
    writer.flush();
  }
  
  private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
    in.defaultReadObject();
    writer = new PrintWriter(new FileOutputStream(filename, true));
  }

  @Override
  public List<Task> makeTasks(int chainCount) throws MC3KitException {
    List<Task> tasks = new ArrayList<Task>(chainCount);
    for(int i = 0; i < chainCount; i++) {
      tasks.add(new MarginalLikelihoodTask(i));
    }
    return tasks;
  }
  
  private synchronized void takeValue(long iteration, Chain chain, double meanLogLikes, double meanAutoLogLikes) {
    int chainCount = chain.getChainCount();
    
    if(collector == null) {
      collector = new Collector<MargLikeValue>(chainCount);
    }
    
    List<MargLikeValue> values = collector.takeValue(
      iteration, chain.getChainId(), new MargLikeValue(chain.getLikelihoodHeatExponent(), meanLogLikes, meanAutoLogLikes)
    );
    
    if(values != null) {
      double[] taus = new double[chainCount];
      double[] means = new double[chainCount];
      double[] autoMeans = new double[chainCount];
      
      // Go hottest to coldest
      for(int i = 0; i < chainCount; i++) {
        MargLikeValue vals = values.get(chainCount - 1 - i);
        taus[i] = vals.tau;
        means[i] = vals.meanLogLikes;
        autoMeans[i] = vals.meanAutoLogLikes;
      }
      
      double margLike = integrateTrapezoid(taus, means);
      double autoMargLike = integrateTrapezoid(taus, autoMeans);
      
      writer.printf("%d\t%.3f\t%.3f", iteration, margLike, autoMargLike);
      writer.println();
      writer.flush();
    }
  }

  
  private class MarginalLikelihoodTask implements Task {
    int chainId;
    
    long iterationCount;
    long logLikeCount;
    double sumLogLikes;
    List<Double> autoLogLikes;
    double sumAutoLogLikes;
    
    @SuppressWarnings("unused")
    protected MarginalLikelihoodTask() { }
    
    MarginalLikelihoodTask(int chainId) {
      this.chainId = chainId;
      
      logLikeCount = 0;
      sumLogLikes = 0.0;
      
      autoLogLikes = new LinkedList<Double>();
      sumAutoLogLikes = 0;
    }
    
    @Override
    public int[] getChainIds() {
      return new int[] { chainId };
    }

    @Override
    public void step(Chain[] chains) throws MC3KitException {
      iterationCount++;
      
      if(iterationCount % thin == 0) {
        double logLike = chains[0].getModel().getLogLikelihood();
        
        // Marginal likelihood from all samples after burn-in period
        if(iterationCount >= burnIn) {
          logLikeCount++;
          sumLogLikes += logLike;
        }
        
        // Marginal likelihood from last half of samples
        if(autoLogLikes.size() % 2 == 1) {
          double oldLogLike = autoLogLikes.remove(0);
          sumAutoLogLikes -= oldLogLike;
        }
        autoLogLikes.add(logLike);
        sumAutoLogLikes += logLike;
        
        takeValue(iterationCount, chains[0], sumLogLikes / logLikeCount, sumAutoLogLikes / autoLogLikes.size());
      }
    }
  }
  
  private class MargLikeValue implements Serializable {
    double tau;
    double meanLogLikes;
    double meanAutoLogLikes;
    
    MargLikeValue(double tau, double meanLogLikes, double meanAutoLogLikes)
    {
      this.tau = tau;
      this.meanLogLikes = meanLogLikes;
      this.meanAutoLogLikes = meanAutoLogLikes;
    }
  }
}
