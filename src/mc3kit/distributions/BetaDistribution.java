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

import mc3kit.*;
import mc3kit.proposal.MHUniformProposer;
import cern.jet.random.Beta;
import static mc3kit.util.Math.*;
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class BetaDistribution extends DoubleDistribution {
  ModelEdge alphaEdge;
  double alpha;
  
  ModelEdge betaEdge;
  double beta;
  
  double logBetaAB;
  
  protected BetaDistribution() { }
  
  public BetaDistribution(Model model) {
    this(model, null);
  }

  public BetaDistribution(Model model, String name) {
    this(model, name, 1.0, 1.0);
  }
  
  public BetaDistribution(Model model, double alpha, double beta) throws IllegalArgumentException {
    this(model, null, alpha, beta);
  }
  
  public BetaDistribution(Model model, String name, double alpha, double beta) throws IllegalArgumentException {
    super(model, name);
    
    if(Double.isNaN(alpha) || Double.isNaN(beta)) {
      throw new IllegalArgumentException("alpha or beta is NaN");
    }
    
    if(Double.isInfinite(alpha) || Double.isInfinite(beta)) {
      throw new IllegalArgumentException("alpha or beta is infinite");
    }
    
    if(alpha <= 0.0 || beta <= 0.0) {
      throw new IllegalArgumentException("alpha and beta must be larger than zero");
    }
    
    this.alpha = alpha;
    this.beta = beta;
    updateConstant();
  }
  
  public <T extends ModelNode & DoubleValued> BetaDistribution setAlpha(T alphaNode) throws MC3KitException {
    alphaEdge = updateEdge(alphaEdge, alphaNode);
    return this;
  }
  
  public <T extends ModelNode & DoubleValued> BetaDistribution setBeta(T betaNode) throws MC3KitException {
    betaEdge = updateEdge(betaEdge, betaNode);
    return this;
  }
  
  private void updateConstant() {
    logBetaAB = logBeta(alpha, beta);
  }

  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHUniformProposer(varName, 0.0, 1.0);
  }

  @Override
  public boolean valueIsValid(double value) {
    return value > 0.0 && value < 1.0;
  }

  @Override
  public double getLogP(Variable var) {
    double x = ((DoubleVariable)var).getValue();
    double logP = (alpha - 1.0) * log(x) + (beta - 1.0) * log1p(-x) - logBetaAB;
    return logP;
  }

  @Override
  public void sample(Variable var) {
    assert alpha > 0.0;
    assert beta > 0.0;
    
    Beta betaGen = new Beta(alpha, beta, getRng());
    double newVal;
    do {
      newVal = betaGen.nextDouble();
    } while(newVal == 0.0 || newVal == 1.0);
    
    assert !Double.isNaN(newVal);
    assert newVal >= 0.0;
    assert newVal > 0.0;
    assert newVal <= 1.0;
    assert newVal < 1.0;
    
    ((DoubleVariable)var).setValue(newVal);
  }

  @Override
  public boolean update() {
    boolean changed = false;
    if(alphaEdge != null) {
      double oldAlpha = alpha;
      alpha = getDoubleValue(alphaEdge);
      if(alpha != oldAlpha) {
        changed = true;
      }
    }
    if(betaEdge != null) {
      double oldBeta = beta;
      beta = getDoubleValue(betaEdge);
      if(beta != oldBeta) {
        changed = true;
      }
    }
    if(changed) {
      updateConstant();
    }
    return changed;
  }
}
