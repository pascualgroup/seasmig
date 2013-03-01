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

import static mc3kit.util.Math.getRandomPermutation;
import static mc3kit.util.Utils.makeMap;
import static java.lang.String.format;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import mc3kit.Chain;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.Step;
import mc3kit.Task;
import mc3kit.Variable;
import mc3kit.VariableProposer;

import cern.jet.random.Uniform;
import cern.jet.random.engine.RandomEngine;

@SuppressWarnings("serial")
public class UnivariateProposalStep implements Step {
  
  double targetAcceptanceRate = 0.25;
  long tuneFor = 1000;
  long tuneEvery = 100;

  protected UnivariateProposalStep() { }
  
  public UnivariateProposalStep(double targetAcceptanceRate, long tuneFor, long tuneEvery) {
    this.targetAcceptanceRate = targetAcceptanceRate;
    this.tuneFor = tuneFor;
    this.tuneEvery = tuneEvery;
  }

  @Override
  public List<Task> makeTasks(int chainCount) throws MC3KitException {
    List<Task> tasks = new ArrayList<Task>();
    for(int i = 0; i < chainCount; i++) {
      tasks.add(new UnivariateProposalTask(i));
    }
    return tasks;
  }
  
  /*** TASK INTERFACE IMPLEMENTATION ***/
  
  class UnivariateProposalTask implements Task {
    boolean initialized;
    
    int chainId;
    
    private long iterationCount;
    
    VariableProposer[] proposers;
    
    /*** CONSTRUCTOR ***/
    
    public UnivariateProposalTask(int chainId) {
      this.chainId = chainId;
    }
    
    
    /*** TASK INTERFACE IMPLEMENTATION ***/
    
    @Override
    public int[] getChainIds() {
      return new int[] { chainId };
    }

    @Override
    public void step(Chain[] chains) throws MC3KitException {
      assert (chains.length == 1);

      Chain chain = chains[0];
      Logger logger = chain.getLogger();
      
      Model model = chain.getModel();
      initialize(model);
      
      RandomEngine rng = chain.getRng();

      Uniform unif = new Uniform(rng);

      // Run all proposers in random order
      for(int i : getRandomPermutation(proposers.length, unif)) {
        proposers[i].step(model);
      }

      iterationCount++;
      
      // Write out acceptance rates
      if(iterationCount % tuneEvery == 0 && logger.isLoggable(Level.INFO)) {
        Map<String, Double> acceptanceRates = new LinkedHashMap<String, Double>();
        for(VariableProposer proposer : proposers) {
          acceptanceRates.put(proposer.getName(), proposer.getAcceptanceRate());
        }
        Map<String, Object> infoObj = makeMap(
          "iteration", iterationCount,
          "chainId", chainId,
          "acceptanceRates", acceptanceRates
        );
        logger.log(Level.INFO, "UnivariateProposalStep acceptance rates", infoObj);
      }
      
      // If we're still in the tuning period, tune
      if((iterationCount <= tuneFor) && iterationCount % tuneEvery == 0) {
        for(VariableProposer proposer : proposers) {
          proposer.tune(targetAcceptanceRate);
          proposer.resetTuningPeriod();
        }
      }
    }
    
    /*** PRIVATE METHODS ***/
    
    private void initialize(Model model) throws MC3KitException {
      if(initialized) return;
      initialized = true;

      String[] varNames = model.getUnobservedVariableNames();
      proposers = new VariableProposer[varNames.length];
      for(int i = 0; i < varNames.length; i++) {
        proposers[i] = makeVariableProposer(model, varNames[i]);
      }
    }
    
    private VariableProposer makeVariableProposer(Model model, String varName) throws MC3KitException {
      Variable var = model.getVariable(varName);
      
      if(var == null) {
        throw new MC3KitException(format("No variable named %s", varName));
      }
      
      return var.makeProposer();
    }
  }
}
