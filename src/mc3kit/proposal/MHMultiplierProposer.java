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

import cern.jet.random.engine.RandomEngine;
import mc3kit.*;
import static java.lang.Math.exp;
import static mc3kit.util.Math.*;
import static java.lang.String.format;

@SuppressWarnings("serial")
public class MHMultiplierProposer extends VariableProposer {

  double lambda;
  
  protected MHMultiplierProposer() { }
  
  public MHMultiplierProposer(String name) {
    super(name);
    this.lambda = 0.5;
  } 

  @Override
  public void step(Model model) throws MC3KitException {
    Chain chain = model.getChain();
    RandomEngine rng = chain.getRng();
    
    chain.getLogger().finest("MHMultiplierProposer stepping");
    
    double oldLogLikelihood = model.getLogLikelihood();
    double oldLogPrior = model.getLogPrior();
    
    chain.getLogger().finest(format("oldLP, oldLL: %f, %f", oldLogPrior, oldLogLikelihood));
    
    DoubleVariable rv = model.getDoubleVariable(getName());
    
    double oldValue = rv.getValue();
    double logMultiplier = lambda * (rng.nextDouble() - 0.5);
    double multiplier = exp(logMultiplier);
    
    double newValue = multiplier * oldValue;
    
    chain.getLogger().finest(format("oldVal, newVal = %f, %f", oldValue, newValue));
    
    if(!rv.valueIsValid(newValue))
    {
      recordRejection();
      return;
    }
    
    model.beginProposal();
    rv.setValue(newValue);
    model.endProposal();

    double newLogPrior = model.getLogPrior();
    double newLogLikelihood = model.getLogLikelihood();
    
    chain.getLogger().finest(format("newLP, newLL: %f, %f", newLogPrior, newLogLikelihood));
    
    boolean accepted = shouldAcceptMetropolisHastings(rng,
      chain.getPriorHeatExponent(), chain.getLikelihoodHeatExponent(),
      oldLogPrior, oldLogLikelihood,
      newLogPrior, newLogLikelihood,
      logMultiplier
    );
    
    if(accepted)
    {
      model.acceptProposal();
      recordAcceptance();
    }
    else
    {
      model.beginRejection();
      rv.setValue(oldValue);
      model.endRejection();
      
      recordRejection();
    }
  }

  @Override
  public void tune(double targetRate) throws MC3KitException {
    lambda = adjustTuningParameter(lambda, getAcceptanceRate(), targetRate);
  }
}
