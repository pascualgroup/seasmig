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

import static java.lang.Math.log;
import cern.jet.random.Exponential;
import mc3kit.*;
import mc3kit.proposal.*;

@SuppressWarnings("serial")
public class ExponentialDistribution extends DoubleDistribution {

  double rate;
  
  ModelEdge rateEdge;
  ModelEdge scaleEdge;
  
  protected ExponentialDistribution() { }
  
  public ExponentialDistribution(Model model) {
    this(model, null);
  }
  
  public ExponentialDistribution(Model model, String name) {
    this(model, name, 1.0);
  }
  
  public ExponentialDistribution(Model model, double rate) {
    this(model, null, rate);
  }
  
  public ExponentialDistribution(Model model, String name, double rate) {
    super(model, name);
    this.rate = rate;
  }

  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHMultiplierProposer(varName);
  }
  
  public <T extends ModelNode & DoubleValued> void setRate(T rateNode) throws MC3KitException {
    scaleEdge = updateEdge(scaleEdge, null);
    rateEdge = updateEdge(rateEdge, rateNode);
  }
  
  public <T extends ModelNode & DoubleValued> void setScale(T scaleNode) throws MC3KitException {
    rateEdge = updateEdge(rateEdge, null);
    scaleEdge = updateEdge(scaleEdge, scaleNode);
  }

  @Override
  public double getLogP(Variable v) {
    double x = ((DoubleVariable)v).getValue();

    if(scaleEdge != null) {
      assert rateEdge == null;
      return getLogPScale(x, getDoubleValue(scaleEdge));
    }
    else {
      assert scaleEdge == null;
      return getLogPRate(x, rateEdge == null ? rate : getDoubleValue(rateEdge));
    }
  }
  
  public static double getLogPRate(double x, double rate) {
    assert x > 0;
    assert rate > 0;
    
    return log(rate) - rate * x;
  }
  
  public static double getLogPScale(double x, double scale) {
    assert x > 0;
    assert scale > 0;
    
    return -log(scale) - x / scale;
  }

  @Override
  public boolean valueIsValid(double val) {
    if(Double.isInfinite(val) || Double.isNaN(val) || val <= 0.0) {
      return false;
    }
    return true;
  }

  @Override
  public void sample(Variable var) {
    double rate;
    if(rateEdge != null) {
      rate = getDoubleValue(rateEdge);
    }
    else if(scaleEdge != null) {
      rate = 1.0 / getDoubleValue(scaleEdge);
    }
    else {
      rate = this.rate;
    }
    double newVal = new Exponential(rate, getRng()).nextDouble();
    
    assert !Double.isNaN(newVal);
    assert !Double.isInfinite(newVal);
    assert newVal > 0.0;
    
    ((DoubleVariable)var).setValue(newVal);
  }
}
