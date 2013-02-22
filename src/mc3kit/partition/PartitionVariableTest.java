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

package mc3kit.partition;

import static org.junit.Assert.*;

import mc3kit.Chain;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.MC3KitException;
import mc3kit.MCMC;
import mc3kit.Model;
import mc3kit.ModelFactory;
import mc3kit.distributions.NormalDistribution;
import mc3kit.proposal.UnivariateProposalStep;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PartitionVariableTest {

  @Before
  public void setUp() throws Exception {
  }

  @After
  public void tearDown() throws Exception {
  }

  @SuppressWarnings("serial")
  @Test
  public void normalMixture() throws Throwable {
    long burnIn = 5000;
    long iterCount = 40000;
    
    MCMC mcmc = new MCMC();
    mcmc.setRandomSeed(100L);
    
    mcmc.setModelFactory(new ModelFactory() {
      
      @Override
      public Model createModel(Chain initialChain) throws MC3KitException {
        Model m = new Model(initialChain);
        
        m.beginConstruction();
        
        DoubleDistribution[] d = new DoubleDistribution[2];
        d[0] = new NormalDistribution(m, -1.0, 1.0);
        d[1] = new NormalDistribution(m, 3.0, 1.0);
        
        DoubleVariable[] v = new DoubleVariable[50];
        for(int i = 0; i < 50; i++) {
          v[i] = new DoubleVariable(m, "v" + i);
        }
        
        PartitionVariable pv = new PartitionVariable(m, "part", v.length, d.length);
        pv.associateVariablesWithDistributions(v, d);
        pv.sample();
        
        m.endConstruction();
        
        return m;
      }
    });

    UnivariateProposalStep proposalStep = new UnivariateProposalStep(0.25, 100, burnIn);
    mcmc.addStep(proposalStep);
    
    // Run, collect statistics, and check moments against expected distribution
    double sum = 0;
    double sumSq = 0;
    for(long i = 0; i < iterCount; i++) {
      System.err.println("stepping " + i);
      mcmc.step();
      mcmc.getModel().recalculate();
      
      assertEquals(i + 1, mcmc.getIterationCount());
      
      if(i >= burnIn) {
        for(int j = 0; j < 50; j++) {
          double val = mcmc.getModel().getDoubleVariable("v" + j).getValue();
          sum += val;
          sumSq += val * val;
        }
      }
    }
    
    double N = 50 * (iterCount - burnIn);
    
    double mean = sum / N;
    System.err.printf("mean = %f\n", mean);
    assertEquals(1.0, mean, 0.02);
    
    double var = N / (N - 1) * (sumSq/N - mean * mean);
    System.err.printf("var = %f\n", var);
    assertEquals(5.0, var, 0.02);
  }
}
