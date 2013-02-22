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

import static cern.jet.stat.Gamma.logGamma;
import static java.lang.Math.log;
import cern.jet.random.Exponential;
import mc3kit.*;
import mc3kit.proposal.*;

@SuppressWarnings("serial")
public class GammaDistribution extends DoubleDistribution {
  
  double shape;
  ModelEdge shapeEdge;
  
  double rate;
  ModelEdge rateEdge;
  
  ModelEdge scaleEdge;
  
  protected GammaDistribution() { }
  
  public GammaDistribution(Model model) {
    this(model, null);
  }
  
  public GammaDistribution(Model model, String name) {
    this(model, name, 1.0, 1.0);
  }
  
  public GammaDistribution(Model model, double shape, double rate) {
    this(model, null, shape, rate);
  }
  
  public GammaDistribution(Model model, String name, double shape, double rate) {
    super(model, name);
    this.shape = shape;
    this.rate = rate;
  }
  
  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHMultiplierProposer(varName);
  }
  
  public <T extends ModelNode & DoubleValued> GammaDistribution setShape(T shapeNode) throws MC3KitException {
    shapeEdge = updateEdge(shapeEdge, shapeNode);
    return this;
  }
  
  public <T extends ModelNode & DoubleValued> GammaDistribution setRate(T rateNode) throws MC3KitException {
    scaleEdge = updateEdge(scaleEdge, null);
    rateEdge = updateEdge(rateEdge, rateNode);
    return this;
  }
  
  public <T extends ModelNode & DoubleValued> GammaDistribution setScale(T scaleNode) throws MC3KitException {
    rateEdge = updateEdge(rateEdge, null);
    scaleEdge = updateEdge(scaleEdge, scaleNode);
    return this;
  }

  @Override
  public double getLogP(Variable v) {
    double x = ((DoubleVariable)v).getValue();
    
    double shape = shapeEdge == null ? this.shape : getDoubleValue(shapeEdge);
    
    if(scaleEdge != null) {
      assert rateEdge == null;
      getLogPScale(x, shape, getDoubleValue(scaleEdge));
    }
    
    double rate = rateEdge == null ? this.rate : getDoubleValue(rateEdge);
    return getLogPRate(x, shape, rate);
  }
  
  public static double getLogPRate(double x, double shape, double rate) {
    assert x > 0;
    assert shape > 0;
    assert rate > 0;
    
    return shape * log(rate) - logGamma(shape) + (shape - 1.0) * log(x) - rate * x;
  }
  
  public static double getLogPScale(double x, double shape, double scale) {
    assert x > 0;
    assert shape > 0;
    assert scale > 0;
    
    return - shape * log(scale) - logGamma(shape) + (shape - 1.0) * log(x) - x / scale;
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
      rate = 1.0;
    }
    
    double newVal = new Exponential(rate, getRng()).nextDouble();
    assert !Double.isNaN(newVal);
    assert !Double.isInfinite(newVal);
    assert newVal > 0.0;
    ((DoubleVariable)var).setValue(newVal);
  }
}
