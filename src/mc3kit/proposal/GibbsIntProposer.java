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

import static java.lang.String.format;
import cern.jet.random.engine.RandomEngine;
import mc3kit.*;
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class GibbsIntProposer extends VariableProposer {
	
  double min; // TODO: int...
  double max; // TODO: int...

  protected GibbsIntProposer() { }
  
  public GibbsIntProposer(String name, double min, double max) {
    super(name);
    this.min=min;
    this.max=max;
  }

  @Override
  public void step(Model model) throws MC3KitException {
    Chain chain = model.getChain();
    RandomEngine rng = chain.getRng();
    
    chain.getLogger().finest("GibbsIntProposer stepping");

    double oldLogPrior = model.getLogPrior();
    double oldLogLike = model.getLogLikelihood();
    
    chain.getLogger().finest(format("oldLP, oldLL: %f, %f", oldLogPrior, oldLogLike));
    
    IntVariable v = model.getIntVariable(getName());
    int oldValue = v.getValue();
    
    // TODO: check this...
    int newValue =  (int)oldValue+1 + (int)(rng.nextDouble() * ((max-1 - min) + 1));
    if (newValue>max) {
    	newValue=newValue-(int)max+(int)min-1;
    }
    
    // Propose from uniform distribution except old value....
    model.beginProposal();
    v.setValue(newValue);
    model.endProposal();

    double newLogPrior = model.getLogPrior();
    double newLogLike = model.getLogLikelihood();
    
    chain.getLogger().finest(format("newLP, newLL: %f, %f", newLogPrior, newLogLike));
    
    // TODO: ask Ed....
    // For this to be a Gibbs step, the ratio of the probability of acceptance (change the value)
    // to rejection (keep the value what it was) needs to be the ratio of the conditional densities
    // of those two states:
    // p / (1 - p) = f(!oldValue) / f(oldValue) = r
    // i.e.
    // p = r / (1 + r)
    // where f = prior^priorHeatExponent x likelihood^heatExponent
    double priorHeatExp = chain.getPriorHeatExponent();
    double likeHeatExp = chain.getLikelihoodHeatExponent();
    double logR = priorHeatExp * (newLogPrior - oldLogPrior) + likeHeatExp * (newLogLike - oldLogLike);
    double logP = logR - log1p(exp(logR));
    
    boolean accepted = log(rng.nextDouble()) < logP;
    if(accepted) {
      model.acceptProposal();
      recordAcceptance();
    }
    else {
      model.beginRejection();
      v.setValue(oldValue);
      model.endRejection();
      recordRejection();
    }
  }

}
