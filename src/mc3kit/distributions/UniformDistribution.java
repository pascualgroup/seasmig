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
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class UniformDistribution extends DoubleDistribution {
  double min;
  double max;
  double logP;
  
  protected UniformDistribution() { }
  
  public UniformDistribution(Model model) {
    this(model, null);
  }

  public UniformDistribution(Model model, String name) {
    this(model, name, 0.0, 1.0);
  }
  
  public UniformDistribution(Model model, double min, double max) throws IllegalArgumentException {
    this(model, null, min, max);
  }
  
  public UniformDistribution(Model model, String name, double min, double max) throws IllegalArgumentException {
    super(model, name);
    
    if(Double.isNaN(min) || Double.isNaN(max)) {
      throw new IllegalArgumentException("min or max is NaN");
    }
    
    if(min >= max) {
      throw new IllegalArgumentException("min must be less than max.");
    }
    
    if(Double.isInfinite(min) || Double.isInfinite(max)) {
      throw new IllegalArgumentException("min or max is infinite");
    }
    
    this.min = min;
    this.max = max;
    logP = getLogP(min, max);
  }

  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHUniformProposer(varName, min, max);
  }

  @Override
  public boolean valueIsValid(double value) {
    return value > min && value < max;
  }

  @Override
  public double getLogP(Variable var) {
    return logP;
  }

  @Override
  public void sample(Variable var) {
    double newVal = min + (max - min) * getRng().nextDouble();
    
    assert !Double.isNaN(newVal);
    assert !Double.isInfinite(newVal);
    assert newVal > min;
    assert newVal < max;
    
    ((DoubleVariable)var).setValue(newVal);
  }

  @Override
  public boolean update() {
    return false;
  }
  
  public static double getLogP(double min, double max) {
    return -log(max - min);
  }
}
