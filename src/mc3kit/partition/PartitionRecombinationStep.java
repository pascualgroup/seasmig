package mc3kit.partition;

import static java.lang.String.format;
import static mc3kit.util.Utils.makeMap;

import java.util.*;
import java.util.logging.*;

import static java.lang.Math.*;

import cern.jet.random.Uniform;
import cern.jet.random.engine.RandomEngine;

import mc3kit.*;
import mc3kit.util.*;

@SuppressWarnings("serial")
public class PartitionRecombinationStep implements Step {
  ChainParity parity;
  long statsEvery;
  
  private Collector<Stats> collector;
  private int chainCount;
  
  protected PartitionRecombinationStep() { }
  
  public PartitionRecombinationStep(ChainParity parity, long statsEvery) {
    this.parity = parity;
    this.statsEvery = statsEvery;
  }
  
  @Override
  public List<Task> makeTasks(int chainCount) {
    this.chainCount = chainCount;
    List<Task> tasks = new ArrayList<Task>(chainCount);
    for(int i = parity == ChainParity.EVEN ? 0 : 1; i + 1 < chainCount; i += 2) {
      tasks.add(new PartitionRecombiner(i, i+1));
    }
    collector = new Collector<Stats>(tasks.size());
    return tasks;
  }
  
  private void takeStats(MCMC mcmc, long iterationCount, int[] chainIds, Map<String, Long> acceptanceCounts) {
    List<Stats> statsList = collector.takeValue(iterationCount, chainIds[0]/2, new Stats(chainIds, acceptanceCounts));
    
    // If all the stats are in, write them to the log.
    Logger logger = mcmc.getLogger("mc3kit.partition.PartitionRecombinationStep");
    if(statsList != null && logger.isLoggable(Level.INFO))
    {
      Map<String, Object> allStats = new LinkedHashMap<String, Object>();
      for(int i = (parity == ChainParity.EVEN ? 0 : 1); i + 1 < chainCount; i += 2)
      {
        Stats stats = statsList.get(i/2);
        allStats.put(
            format("(%d,%d)", stats.chainIds[0], stats.chainIds[1]),
            stats.acceptanceRates
        );
      }
      
      logger.log(Level.INFO,
        "PartitionRecombinationStep acceptance rates",
        makeMap("parity", parity.toString(), "stats", allStats)
      );
    }
  }
  
  
  /*** RECOMBINER TASK CLASS ***/
  
  private class PartitionRecombiner implements Task {
    int[] chainIds;
    long iterationCount;
    
    List<String> varNames;
    boolean initialized;
    
    private Map<String, Long> acceptanceCounts;
    
    PartitionRecombiner(int... chainIds) {
      this.chainIds = chainIds;
      this.iterationCount = 0;
    }
    
    @Override
    public int[] getChainIds() {
      return chainIds;
    }

    @Override
    public void step(Chain[] chains) throws MC3KitException {
      initialize(chains[0].getModel());
      iterationCount++;
      
      for(String varName : varNames) {
        step(chains, varName);
      }
      
      if(iterationCount % statsEvery == 0) {
        takeStats(chains[0].getMCMC(), iterationCount, chainIds, acceptanceCounts);
        acceptanceCounts = new HashMap<String, Long>();
      }
    }
    
    private void initialize(Model model) {
      if(initialized) {
        return;
      }
      initialized = true;
      
      varNames = new ArrayList<String>();
      for(Variable var : model.getUnobservedVariables()) {
        if(var instanceof PartitionVariable) {
          varNames.add(var.getName());
        }
      }
      acceptanceCounts = new HashMap<String, Long>(varNames.size());
    }
    
    private void step(Chain[] chains, String varName) throws MC3KitException {
      RandomEngine rng = chains[0].getRng();
      
      Model model0 = chains[0].getModel();
      PartitionVariable var0 = (PartitionVariable)model0.getVariable(varName);
      
      if(var0.getGroupCount() == 1) {
        return;
      }
      
      double oldLP0 = model0.getLogPrior();
      double oldLL0 = model0.getLogLikelihood();
      
      Model model1 = chains[1].getModel();
      PartitionVariable var1 = (PartitionVariable)model1.getVariable(varName);
      double oldLP1 = model1.getLogPrior();
      double oldLL1 = model1.getLogLikelihood();
      
      // Get new assignments
      Uniform unif = new Uniform(rng);
      int[] oldAssign0 = var0.assignment.clone();
      int[] oldAssign1 = var1.assignment.clone();
      int[][] newAssign = generateNewAssignments(rng, unif, var0.getGroupCount(), oldAssign0, oldAssign1, var0.allowsEmptyGroups);
      
      
      // Make proposals
      model0.beginProposal();
      assignGroups(var0, newAssign[0]);
      model0.endProposal();
      double newLP0 = model0.getLogPrior();
      double newLL0 = model0.getLogLikelihood();
      
      model1.beginProposal();
      assignGroups(var1, newAssign[1]);
      model1.endProposal();
      double newLP1 = model1.getLogPrior();
      double newLL1 = model1.getLogLikelihood();
      
      // Calculate acceptance probability and determine acceptance
      double tp0 = chains[0].getPriorHeatExponent();
      double tl0 = chains[0].getLikelihoodHeatExponent();
      double tp1 = chains[1].getPriorHeatExponent();
      double tl1 = chains[1].getLikelihoodHeatExponent();
      double logR = tp0 * newLP0 + tl0 * newLL0 + tp1 * newLP1 + tl1 * newLL1
                  - tp0 * oldLP0 - tl0 * oldLL0 - tp1 * oldLP1 - tl1 * oldLL1;
      boolean accepted = (logR >= 0.0) || (log(rng.nextDouble()) < logR);
      
      if(accepted) {
        model0.acceptProposal();
        model1.acceptProposal();
        
        if(acceptanceCounts.containsKey(varName)) {
          acceptanceCounts.put(varName, acceptanceCounts.get(varName) + 1);
        }
        else {
          acceptanceCounts.put(varName, 1L);
        }
      }
      else {
        model0.beginRejection();
        assignGroups(var0, oldAssign0);
        model0.endRejection();
        
        model1.beginRejection();
        assignGroups(var1, oldAssign1);
        model1.endRejection();
      }
    }
    
    private int[][] generateNewAssignments(RandomEngine rng, Uniform unif, int groupCount,
        int[] oldAssign0, int[] oldAssign1, boolean allowEmptyGroups) {
      int n = oldAssign0.length;
      
      int[] newAssign0 = new int[n];
      int[] newAssign1 = new int[n];
      boolean hasEmptyGroups;
      do {
        int nRecomb = 2;
        IterableBitSet recombPtBS = mc3kit.util.Random.uniformRandomSubset(unif, n, nRecomb);
        int[] recombPts = recombPtBS.getSetBits();
        
        int[] newCount0 = new int[groupCount];
        int[] newCount1 = new int[groupCount];
        
        int recombIndex = 0;
        boolean swapping = true;
        for(int i = 0; i < n; i++) {
          int j = (recombPts[0] + i) % n;
          
          newAssign0[j] = swapping ? oldAssign1[j] : oldAssign0[j];
          newAssign1[j] = swapping ? oldAssign0[j] : oldAssign1[j];
          
          newCount0[newAssign0[j]]++;
          newCount1[newAssign1[j]]++;
          
          if(recombIndex < nRecomb - 1 && j == recombPts[recombIndex + 1]) {
            recombIndex++;
            swapping = !swapping;
          }
        }
        
        hasEmptyGroups = false;
        if(!allowEmptyGroups) {
          for(int i = 0; i < groupCount; i++) {
            if(newCount0[i] == 0 || newCount1[i] == 0) {
              hasEmptyGroups = true;
              break;
            }
          }
        }
        
      } while(!allowEmptyGroups && hasEmptyGroups);
      
      return new int[][] {newAssign0, newAssign1};
    }
    
    private void assignGroups(PartitionVariable var, int[] assignments) throws MC3KitException {
      for(int i = 0; i < assignments.length; i++) {
        int oldGroup = var.getGroupId(i);
        if(oldGroup != assignments[i]) {
          var.setGroup(i, assignments[i]);
        }
      }
    }
  }
  
  private class Stats
  {
    int[] chainIds;
    Map<String, Double> acceptanceRates;
    
    Stats(int[] chainIds, Map<String, Long> acceptanceCounts)
    {
      this.chainIds = chainIds;
      acceptanceRates = new HashMap<String, Double>();
      for(String key : acceptanceCounts.keySet()) {
        acceptanceRates.put(key, acceptanceCounts.get(key) / (double)statsEvery);
      }
    }
  }
}
