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

package mc3kit.example;

import mc3kit.*;
import mc3kit.distributions.*;
import static java.lang.String.format;
import static mc3kit.util.Math.*;

@SuppressWarnings("serial")
public class ExampleModel extends Model {
  double[] data;
  
  DoubleVariable[] means;
  DoubleVariable[] precs;
  DoubleVariable weight;
  LikelihoodVariable likeVar;
  
  double oldLogP; // to be restored upon rejection
  
  protected ExampleModel() { }

  public ExampleModel(Chain initialChain, double[] data) throws MC3KitException {
    super(initialChain);
    this.data = data;
    
    // Initial construction of the model must be bounded by
    // beginConstruction() and endConstruction() and should include:
    // (1) Adding random variables to model, both parameters (unobserved)
    //     and data (observed)
    // (2) Connection of random variables to other nodes (distributions
    //     and functions) and to each other
    // (3) Setting any parameters to specific desired initial values.
    // 
    // The endConstruction() call causes the model to first sample
    // initial values for any parameters that have not already been set,
    // and then to calculate the log-prior and log-likelihood by traversing
    // the model graph.
    beginConstruction();
    
    // Normal(mean=0, stdDev=8) prior on distribution means
    DoubleDistribution meanPrior = new NormalDistribution(this, 0.0, 8.0);
    
    // Gamma(shape=2, rate=2) prior on distribution precision
    DoubleDistribution precPrior = new GammaDistribution(this, 2.0, 2.0);
    
    // Uniform prior on the weight of the first distribution
    weight = new DoubleVariable(this, "weight", new UniformDistribution(this, 0.0, 1.0));
    
    // Two means & precisions. Naming scheme means that sample output
    // will look like this:
    //  {
    //    "iteration" : 100,
    //    "logPrior" : -444.3,
    //    "logLikelihood" : -1.5222
    //    "d0" : { "mean" : 2.42, "prec" : 1.11 },
    //    "d1" : { "mean" : -1.4, "prec" : 0.44 }
    //  }
    means = new DoubleVariable[2];
    precs = new DoubleVariable[2];
    for(int i = 0; i < 2; i++) {
      means[i] = new DoubleVariable(this, format("d%d.mean", i), meanPrior);
      precs[i] = new DoubleVariable(this, format("d%d.prec", i), precPrior);
    }
    
    // Custom likelihood variable
    likeVar = new LikelihoodVariable(this);
    
    endConstruction();
  }
  
  private class LikelihoodVariable extends Variable {
    LikelihoodVariable(ExampleModel m) throws MC3KitException
    {
      // Call superclass constructor specifying that this is an
      // OBSERVED random variable (true for last parameter).
      super(m, "likeVar", true);
      
      // Add dependencies between likelihood variable and parameters
      m.addEdge(this, m.weight);
      for(int i = 0; i < 2; i++) {
        m.addEdge(this, m.means[i]);
        m.addEdge(this, m.precs[i]);
      }
    }

    /*
     * The update method is called whenever a node this variable
     * depends on has changed; random variables should typically
     * just recalculate their log-probability and call setLogP here.
     * 
     * The method returns a boolean indicating whether this node actually
     * changed during the update process. Typically you just return true,
     * but this can be used to optimize updating: e.g., a maximum
     * function may get updated very frequently without its value
     * actually changing.
     * 
     * This example assumes that the data is a mixture between two normal
     * distributions.
     */
    @Override
    public boolean update() {
      double logP = 0.0;
      
      for(double x : data) {
        double w = weight.getValue();
        double[] logPs = new double[2];
        for(int i = 0; i < 2; i++) {
          double mean = means[i].getValue();
          double prec = precs[i].getValue();
          logPs[i] = NormalDistribution.getLogPPrecision(mean, prec, x);
        }
        double[] coeffs = new double[] { w, 1.0 - w };
        
        logP += logSumExp(logPs, coeffs);
      }
      setLogP(logP);
      oldLogP = logP;
      
      return true;
    }
    

    /*
     * If you want to avoid calculating the log-probability again
     * when a proposal is rejected, override these methods to restore
     * a cached value of the log-probability.
     */
    @Override
    public boolean updateAfterRejection() {
      setLogP(oldLogP);
      return true;
    }
    
    /*
     * You can override these methods to optimize updating of your node
     * based on which edges changed. By default, the optimized update
     * simply calls the standard update, which should recalculate from
     * scratch.
     */
//    @Override
//    public boolean update(Set<ModelEdge> changedEdges) {
//      return update();
//    }
//    
//    @Override
//    public boolean updateAfterRejection(Set<ModelEdge> changedEdges) {
//      return update(changedEdges);
//    }  
    
  }
}
