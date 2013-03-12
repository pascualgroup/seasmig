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

package mc3kit.distributions;

import static org.junit.Assert.*;

import java.util.logging.Level;

import mc3kit.*;
import mc3kit.proposal.UnivariateProposalStep;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

@SuppressWarnings("serial")
public class BernoulliDistributionTest {

  @Before
  public void setUp() throws Exception {
  }

  @After
  public void tearDown() throws Exception {
  }

  @Test
  public void test() throws Throwable {
    long burnIn = 5000;
    long iterCount = 10000;
    
    MCMC mcmc = new MCMC();
    mcmc.setRandomSeed(100L);
    
    mcmc.setLogLevel(Level.INFO);
    
    ModelFactory mf = new ModelFactory() {
      @Override
      public Model createModel(Chain initialChain) throws MC3KitException {
        Model m = new Model(initialChain);
        
        m.beginConstruction();
        new BinaryVariable(m, "v", new BernoulliDistribution(m, 0.3));
        m.endConstruction();
        
        return m;
      }
    };
    
    mcmc.setModelFactory(mf);
    
    UnivariateProposalStep proposalStep = new UnivariateProposalStep(0.25, 100, burnIn);
    mcmc.addStep(proposalStep);
    
    // Run, collect statistics, and check moments against expected distribution
    int sum = 0;
    for(long i = 0; i < iterCount; i++) {
      mcmc.step();
      mcmc.getModel().recalculate(1e-8);
      
      assertEquals(i + 1, mcmc.getIterationCount());
      
      if(i >= burnIn) {
        boolean val = mcmc.getModel().getBinaryVariable("v").getValue();
        sum += val ? 1 : 0;
      }
    }
    
    double N = iterCount - burnIn;
    
    double mean = sum / (double)N;
    System.err.printf("mean = %f\n", mean);
    assertEquals(0.3, mean, 0.02);
  }
}
