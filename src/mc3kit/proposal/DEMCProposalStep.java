/***
  This file is part of mc3kit.
  
  Copyright (C) 2013 Edward B. Baskerville

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
***/

package mc3kit.proposal;

import static java.lang.Math.*;
import static mc3kit.util.Math.*;
import static java.lang.String.format;
import static mc3kit.util.Utils.makeMap;
import mc3kit.*;
import mc3kit.util.*;

import java.io.Serializable;
import java.util.*;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import com.google.gson.*;

import cern.colt.*;
import cern.colt.function.*;
import cern.colt.matrix.*;
import cern.colt.matrix.impl.*;
import cern.jet.random.*;
import cern.jet.random.engine.*;

@SuppressWarnings("serial")
public class DEMCProposalStep implements Step {
  private static enum CounterType 
  {
    ACCEPTANCE,
    REJECTION,
    IMPOSSIBLE
  }
  
  double targetAcceptanceRate;
  long tuneFor;
  long tuneEvery;
  long historyThin;
  long initialHistoryCount;
  int minBlockSize;
  int maxBlockSize;
  boolean useParallel;
  boolean useLarge;
  boolean useSnooker;
  
  protected DEMCProposalStep() { }
  
  /**
   * Constructor for a differential evolution MCMC (DEMC) step.
   * 
   * Cajo J.F. ter Braak and Jasper A. Vrugt
   * Differential Evolution Markov Chain with snooker updater and fewer chains 
   * Stat Comput (2008) 18: 435--446
   * DOI 10.1007/s11222-008-9104-9
   * 
   * The paper uses multiple chains; this implementation exists on a single chain and uses
   * historical samples exclusively to generate proposals.
   * 
   * The step will do the following:
   * (1) If useParallel is true, then do standard DEMC proposals (difference of historical samples).
   * (2) If useParallel AND useLarge is true, then also do double-size proposals.
   * (3) If useSnooker is on, also perform snooker proposals (project differences into a sensible direction).
   * 
   * For each of these, proposals will be done at multiple block sizes:
   * minBlockSize, 2 * minBlockSize, 4 * minBlockSize, 8 * minBlockSize, ..., N
   * where N = minBlockSize * 2^n and N <= maxBlockSize.
   * 
   * If all three proposal types are on, then a single step will make 3m proposals to each
   * variable, where m = # of different block sizes.
   * 
   * @param targetAcceptanceRate What fraction of proposals should be accepted?
   * @param tuneFor How many iterations the tuning period lasts.
   * @param tuneEvery How often to tune within the tuning period.
   * @param historyThin How many samples go by for every historical sample recorded for the DEMC.
   * @param initialHistoryCount The number of historical samples to accumulate before doing DEMC.
   * @param minBlockSize The smallest number of variables to propose at a time.
   * @param maxBlockSize The
   * @param useParallel
   * @param useLarge
   * @param useSnooker
   */
  public DEMCProposalStep(
      double targetAcceptanceRate,
      long tuneFor,
      long tuneEvery,
      long historyThin,
      long initialHistoryCount,
      int minBlockSize,
      int maxBlockSize,
      boolean useParallel,
      boolean useLarge,
      boolean useSnooker
  ) {
    super();
    this.targetAcceptanceRate = targetAcceptanceRate;
    this.tuneFor = tuneFor;
    this.tuneEvery = tuneEvery;
    this.historyThin = historyThin;
    this.initialHistoryCount = initialHistoryCount;
    this.minBlockSize = minBlockSize;
    this.maxBlockSize = maxBlockSize;
    this.useParallel = useParallel;
    this.useLarge = useLarge;
    this.useSnooker = useSnooker;
  }

  @Override
  public List<Task> makeTasks(int chainCount) throws MC3KitException {
    List<Task> tasks = new ArrayList<Task>();
    for(int i = 0; i < chainCount; i++) {
      tasks.add(new DEMCProposalTask(i));
    }
    return tasks;
  }
  
  /*** TASK INTERFACE IMPLEMENTATION ***/
  
  class DEMCProposalTask implements Task {
    boolean initialized;
    List<String> varNames;
    List<DoubleMatrix1D> history;
    double[] historySums;
    double[] historySumSqs;
    double[] historyMeans;
    double[] historyStdDevs;
    
    int chainId;
    private long iterationCount;
    
    List<BlockSizeManager> blockSizeManagers;
    
    /*** CONSTRUCTOR ***/
    
    public DEMCProposalTask(int chainId) {
      this.chainId = chainId;
    }
    
    @Override
    public int[] getChainIds() {
      return new int[] { chainId };
    }

    @Override
    public void step(Chain[] chains) throws MC3KitException {
      iterationCount++;
      
      assert (chains.length == 1);

      Chain chain = chains[0];
      Logger logger = chain.getLogger();
      if(logger.isLoggable(Level.FINE)) {
        logger.fine(format("DEMCProposalStep stepping %d", chainId));
      }
      Model model = chain.getModel();
      
      initialize(model);
      
      assert(chains.length == 1);
      
      // Only iterate if we have enough initial values
      if(history.size() >= initialHistoryCount) {
        for(BlockSizeManager bsm : blockSizeManagers) {
          bsm.step(chain, model);
        }
      }
      
      // Record history
      if(iterationCount % historyThin == 0) {
        recordHistory(model);
      }
    }
    
    /*** PRIVATE METHODS ***/
    
    private void recordHistory(Model model)
    {
      DoubleMatrix1D vec = makeVector(model);
      history.add(vec);
      
      for(int i = 0; i < vec.size(); i++)
      {
        double xi = vec.getQuick(i);
        historySums[i] += xi;
        historySumSqs[i] += xi * xi;
        historyMeans[i] = historySums[i] / history.size();
        historyStdDevs[i] = historySumSqs[i] / history.size() - historyMeans[i] * historyMeans[i];
      }
    }
    
    DoubleMatrix1D makeVector(Model model) {
      DoubleMatrix1D vec = new DenseDoubleMatrix1D(varNames.size());
      for(int i = 0; i < varNames.size(); i++) {
        vec.setQuick(i, model.getDoubleVariable(varNames.get(i)).getValue());
      }
      return vec;
    }
    
    DoubleMatrix1D makeVector(Model model, int[] block) {
      DoubleMatrix1D vec = new DenseDoubleMatrix1D(block.length);
      for(int i = 0; i < block.length; i++) {
        vec.setQuick(i, model.getDoubleVariable(varNames.get(block[i])).getValue());
      }
      return vec;
    }
    
    boolean vectorIsValid(Model model, int[] block, DoubleMatrix1D xNew) throws MC3KitException {
      boolean valid = true;
      for(int i = 0; i < block.length; i++) {
        if(!model.getDoubleVariable(varNames.get(block[i])).valueIsValid(xNew.get(i))) {
          valid = false;
          break;
        }
      }
      return valid;
    }
    
    void setVector(Model model, int[] block, DoubleMatrix1D xNew) {
      for(int i = 0; i < block.length; i++) {
        model.getDoubleVariable(varNames.get(block[i])).setValue(xNew.get(i));
      }
    }
    
    private void initialize(Model model) throws MC3KitException {
      if(initialized) return;
      initialized = true;
      
      varNames = new ArrayList<String>();
      for(String varName : model.getUnobservedVariableNames()) {
        if(model.getVariable(varName) instanceof DoubleVariable) {
          varNames.add(varName);
        }
      }
      
      history = new ArrayList<DoubleMatrix1D>();
      historySums = new double[varNames.size()];
      historySumSqs = new double[varNames.size()];
      historyMeans = new double[varNames.size()];
      historyStdDevs = new double[varNames.size()];
      
      // Create managers for each block size.
      // Block sizes are minBlockSize, 2*minBlockSize, 4*minBlockSize, ...
      int blockSize = minBlockSize;
      blockSizeManagers = new ArrayList<BlockSizeManager>();
      while(blockSize <= maxBlockSize && blockSize <= varNames.size())
      {
        if(model.getLogger().isLoggable(Level.FINE)) {
          model.getLogger().fine(format("Creating block size manager for block size = %d", blockSize));
        }
        blockSizeManagers.add(new BlockSizeManager(blockSize));
        blockSize *= 2;
      }
    }
    
    private DoubleMatrix1D[] getRandomSamples(RandomEngine rng, int count)
    {
      Uniform unif = new Uniform(rng);
      
      int[] indexes = new int[count];
      DoubleMatrix1D[] samples = new DoubleMatrix1D[count];
      
      int historyCount = history.size();
      assert(historyCount >= count);
      for(int i = 0; i < count; i++)
      {
        boolean done = false;
        while(!done)
        {
          indexes[i] = unif.nextIntFromTo(0, historyCount - 1);
          done = true;
          for(int j = 0; j < i; j++)
          {
            if(indexes[i] == indexes[j])
            {
              done = false;
              break;
            }
          }
        }
        samples[i] = history.get(indexes[i]);
      }
      
      return samples;
    }
    
    /*** MANAGER FOR PROPOSALS FOR A SINGLE BLOCK SIZE ***/
    
    private class BlockSizeManager implements Serializable
    {
      int blockSize;
      
      double snookerGammaFactor;
      double parallelGammaFactor;
      
      MultiCounter<CounterType> parallelSmallCounter;
      MultiCounter<CounterType> parallelLargeCounter;
      MultiCounter<CounterType> snookerCounter;
      
      BlockSizeManager(int scale)
      {
        this.blockSize = scale;
        
        snookerGammaFactor = 1.0;
        parallelGammaFactor = 1.0;
        
        parallelSmallCounter = new MultiCounter<CounterType>();
        parallelLargeCounter = new MultiCounter<CounterType>();
        snookerCounter = new MultiCounter<CounterType>();
      }
      
      void step(Chain chain, Model model) throws MC3KitException
      {
        if(chain.getLogger().isLoggable(Level.FINER)) {
          chain.getLogger().finer(format("Stepping block size %d", blockSize));
        }
        
        proposeDEMC(chain, model);
        recordStats(chain, chainId);
        tune(chain, chainId);
      }
      
      void recordStats(Chain chain, int chainId) throws MC3KitException
      {
        if(history.size() > initialHistoryCount && iterationCount % tuneEvery == 0) {
          Map<String, Object> outputObj = makeMap(
              "iteration", iterationCount,
              "chainId", chainId,
              "blockSize", blockSize,
              "snookerGammaFactor", snookerGammaFactor,
              "snookerRates", getCounterObject(snookerCounter),
              "parallelGammaFactor", parallelGammaFactor,
              "parallelSmallRates", getCounterObject(parallelSmallCounter),
              "parallelLargeRates", getCounterObject(parallelLargeCounter)
            );
          
          chain.getLogger().info(format("Recording stats for %d...", blockSize));
          chain.getLogger().info(new GsonBuilder().setPrettyPrinting().create().toJson(outputObj));
        }
      }
      
      void tune(Chain chain, int chainId) throws MC3KitException
      {
        if(history.size() > initialHistoryCount && iterationCount % tuneEvery == 0)
        {
          Logger logger = chain.getLogger();
          
          if(iterationCount <= tuneFor)
          {
            if(logger.isLoggable(Level.FINE)) {
              logger.fine(format("Tuning for %d...", blockSize));
            }
            
            if(useParallel) {
              double parallelRate = parallelSmallCounter.getRate(CounterType.ACCEPTANCE);
              if(logger.isLoggable(Level.FINE)) {
                logger.fine(format("Old parallelGammaFactor: %f", parallelGammaFactor));
              }
              parallelGammaFactor = adjustTuningParameter(
                parallelGammaFactor, parallelRate, targetAcceptanceRate
              );
              if(logger.isLoggable(Level.FINE)) {
                logger.fine(format("New parallelGammaFactor: %f", parallelGammaFactor));
              }
            }
            
            if(useSnooker) {
              double snookerRate = snookerCounter.getRate(CounterType.ACCEPTANCE);
              if(logger.isLoggable(Level.FINE)) {
                logger.fine(format("Old snookerGammaFactor: %f", snookerGammaFactor));
              }
              snookerGammaFactor = adjustTuningParameter(
                snookerGammaFactor, snookerRate, targetAcceptanceRate
              );
              if(logger.isLoggable(Level.FINE)) {
                logger.fine(format("New snookerGammaFactor: %f", snookerGammaFactor));
              }
            }
          }
          parallelSmallCounter.reset();
          parallelLargeCounter.reset();
          snookerCounter.reset();
          
        }
      }
      
      Map<String, Object> getCounterObject(MultiCounter<CounterType> counter) throws MC3KitException
      {
        return makeMap(
          "count", counter.getCount(),
          "acceptance", counter.getRate(CounterType.ACCEPTANCE),
          "rejection", counter.getRate(CounterType.REJECTION),
          "impossible", counter.getRate(CounterType.IMPOSSIBLE)
        );
      }
      
      void proposeDEMC(Chain chain, Model xModel) throws MC3KitException
      {
        // Implementation of DEMC-Z and DEMC-ZS algorithms from
        //
        // Cajo J.F. ter Braak and Jasper A. Vrugt
        // Differential Evolution Markov Chain with snooker updater and fewer chains 
        // Stat Comput (2008) 18: 435--446
        // DOI 10.1007/s11222-008-9104-9
        // 
        // Modifications: one chain. All non-proposal vectors are sampled
        // from past. No noise: noise needs to be introduced in another
        // step (i.e. one-at-a-time proposals to variables).
        // 
        // Variable names chosen to match those in paper.
        
        Logger logger = xModel.getLogger();
        
        if(logger.isLoggable(Level.FINE)) {
          logger.finer(format("Proposing %d...", blockSize));
        }
        
        double priorHeatExp = chain.getPriorHeatExponent();
        double likeHeatExp = chain.getLikelihoodHeatExponent();
        RandomEngine rng = chain.getRng();
        
        // Random sample from past to determine variable order and alignment
        DoubleMatrix1D refVec = getRandomSamples(rng, 1)[0];
        
        // Get order of entries in a way that makes covarying/anti-covarying
        // entries tend to get lumped together
        int[] entryOrder = getEntryOrder(refVec);
        
        if(logger.isLoggable(Level.FINER)) {
          logger.finer(format("Entry order %d: %s", blockSize, Arrays.toString(entryOrder)));
        }
        
        // Do a snooker, small parallel, and large parallel proposal for each block
        for(int i = 0; i < entryOrder.length; i += blockSize)
        {
          if(logger.isLoggable(Level.FINER)) {
            logger.finer(format("Proposing blockSize %d, blockStart %d", blockSize, i));
          }
          int blockEnd = i + blockSize;
          if(blockEnd > entryOrder.length) blockEnd = entryOrder.length;
          int[] block = Arrays.copyOfRange(entryOrder, i, blockEnd);
          
          if(useParallel) {
            proposeBlockDEMCParallel(priorHeatExp, likeHeatExp, false, block, xModel, rng);
            if(useLarge) {
              proposeBlockDEMCParallel(priorHeatExp, likeHeatExp, true, block, xModel, rng);
            }
          }
          if(useSnooker) {
            proposeBlockDEMCSnooker(priorHeatExp, likeHeatExp, block, xModel, rng);
          }
        }
      }
      
      void proposeBlockDEMCSnooker(
        double priorHeatExp, double likeHeatExp,
        int[] block,
        Model xModel,
        RandomEngine rng) throws MC3KitException
      {
        // Scale factor of 1.7 is optimal for normal/Student posteriors;
        // adjusted by factor to target this distribution
        double gamma = snookerGammaFactor * 1.7;
        
        proposeBlockDEMC(priorHeatExp, likeHeatExp, gamma, false, true, block, xModel, rng);
      }
      
      void proposeBlockDEMCParallel(
          double priorHeatExp, double likeHeatExp,
          boolean isLarge,
          int[] block,
          Model xModel,
          RandomEngine rng) throws MC3KitException 
      {
        // Scale factor of 2.38/sqrt(2d) is optimal for normal/Student
        // posteriors.
        // For "large" proposals, to help avoid local minima,
        // gamma = 2 * base scale factor
        double gamma = parallelGammaFactor * 2.38 / sqrt(2 * block.length);
        if(isLarge)
          gamma *= 2.0;
        
        proposeBlockDEMC(priorHeatExp, likeHeatExp, gamma, isLarge, false, block, xModel, rng);
      }
      
      void proposeBlockDEMC(
          double priorHeatExp, double likeHeatExp,
          double gamma,
          boolean isLarge,
          boolean isSnooker,
          int[] block,
          Model xModel,
          RandomEngine rng) throws MC3KitException 
      {
        Logger logger = xModel.getLogger();
        int d = block.length;
        
        DoubleMatrix1D xOld = makeVector(xModel, block);
        
        DoubleMatrix1D z = null;
        double xMinusZNormOld = 0.0;
        DoubleMatrix1D z1;
        DoubleMatrix1D z2;
        if(isSnooker) {
          DoubleMatrix1D[] samps;
          DoubleMatrix1D xMinusZOld;
          
          do
          {
            // Get three random samples:
            // * z: to define projection vector (x - z)
            // * zR1, zR2: to project onto (x - z) to define difference
            samps = getRandomSamples(rng, 3);
            
            z = samps[0].viewSelection(block);
            xMinusZOld = subtract(xOld, z);
            xMinusZNormOld = norm2(xMinusZOld);
          }
          while(xMinusZNormOld == 0.0);
          DoubleMatrix1D zR1 = samps[1].viewSelection(block);
          DoubleMatrix1D zR2 = samps[2].viewSelection(block);
          
          // Project zR1, zR2 onto x - z
          divideInPlace(xMinusZOld, xMinusZNormOld);
          z1 = project(zR1, xMinusZOld);
          z2 = project(zR2, xMinusZOld);
        }
        else { 
          // Get two random samples to define difference
          DoubleMatrix1D[] samps = getRandomSamples(rng, 2);
          z1 = samps[0].viewSelection(block);
          z2 = samps[1].viewSelection(block);
        }
        
        // Update x
        DoubleMatrix1D xNew = new DenseDoubleMatrix1D(d);
        for(int i = 0; i < d; i++)
        {
          xNew.setQuick(i,
            xOld.getQuick(i)
            + gamma * (z1.getQuick(i) - z2.getQuick(i))
          );
        }
        
        // Perform update
        if(logger.isLoggable(Level.FINE)) {
          logger.fine(format("Setting entries %s to %s", Arrays.toString(block), xNew));
        }
        double oldLogPrior = xModel.getLogPrior();
        double oldLogLike = xModel.getLogLikelihood();
        
        // Calculate final norm of difference, make sure it's not zero
        boolean impossible = false;
        double xMinusZNormNew = 0.0;
        if(isSnooker) { 
          DoubleMatrix1D xMinusZNew = subtract(xNew, z);
          xMinusZNormNew = norm2(xMinusZNew);
          
          if(xMinusZNormNew == 0.0) {
            impossible = true;
          }
        }
        if(!impossible) {
          impossible = !vectorIsValid(xModel, block, xNew);
        }
        
        double newLogPrior = Double.NEGATIVE_INFINITY;
        double newLogLike = Double.NEGATIVE_INFINITY;
        if(!impossible)
        {
          xModel.beginProposal();
          setVector(xModel, block, xNew);
          xModel.endProposal();
          newLogPrior = xModel.getLogPrior();
          newLogLike = xModel.getLogLikelihood();
          
          assert(!Double.isInfinite(newLogPrior));
          assert(!Double.isInfinite(newLogLike));
        }
        
        // Acceptance/rejection
        boolean accepted;
        if(impossible)
        {
          accepted = false;
        }
        else
        {
          double logProposalRatio;
          if(isSnooker) {
            logProposalRatio = (d - 1) * (log(xMinusZNormNew) - log(xMinusZNormOld));
          }
          else {
            logProposalRatio = 0.0;
          }
          
          accepted = shouldAcceptMetropolisHastings(
            rng,
            priorHeatExp, likeHeatExp,
            oldLogPrior, oldLogLike,
            newLogPrior, newLogLike,
            logProposalRatio
          );
        }
        
        MultiCounter<CounterType> counter = isSnooker ? snookerCounter :
          (isLarge ? parallelLargeCounter : parallelSmallCounter);
        if(accepted)
        {
          logger.fine("Accepted");
          xModel.acceptProposal();
          counter.record(CounterType.ACCEPTANCE);
        }
        else if(impossible)
        {
          logger.fine("Impossible");
          counter.record(CounterType.REJECTION, CounterType.IMPOSSIBLE);
        }
        else
        {
          logger.fine("Rejected");
          xModel.beginRejection();
          setVector(xModel, block, xOld);
          xModel.endRejection();
          
          counter.record(CounterType.REJECTION);
        }
      }
      
      // Sort entries by abs(std dev-normalized distance from mean)
      // so that covarying or anti-covarying quantities will tend to
      // cluster together
      int[] getEntryOrder(DoubleMatrix1D x)
      {
        final double[] xRel = new double[x.size()];
        final int[] order = new int[x.size()];
        
        // Generate original order and x values standardized by mean/stddev estimates
        for(int i = 0; i < xRel.length; i++)
        {
          order[i] = i;
          if(historyStdDevs[i] == 0)
            xRel[i] = 0;
          else
            xRel[i] = abs((x.getQuick(i) - historyMeans[i]) / historyStdDevs[i]);
          assert(!Double.isInfinite(xRel[i]));
          assert(!Double.isNaN(xRel[i]));
        }
        
        IntComparator comparator = new IntComparator()
        {
          @Override
          public int compare(int a, int b)
          {
            return xRel[a] == xRel[b] ? 0 : (xRel[a] < xRel[b] ? -1 : 1);
          }
        };
        
        Swapper swapper = new Swapper()
        {
          @Override
          public void swap(int a, int b)
          {
            int tmpOrder;
            double tmpXRel;
            
            tmpOrder = order[a];
            tmpXRel = xRel[a];
            order[a] = order[b];
            xRel[a] = xRel[b];
            order[b] = tmpOrder;
            xRel[b] = tmpXRel;
          }
        };
        
        GenericSorting.quickSort(0, xRel.length, comparator, swapper);
        
        return order;
      }
    }
  }
}
