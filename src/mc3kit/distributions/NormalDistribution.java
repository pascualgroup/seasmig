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

import cern.jet.random.Normal;
import mc3kit.*;

import static java.lang.Math.*;
import static mc3kit.util.Math.*;
import mc3kit.DoubleDistribution;
import mc3kit.DoubleVariable;
import mc3kit.MC3KitException;
import mc3kit.VariableProposer;

import mc3kit.proposal.*;

@SuppressWarnings("serial")
public class NormalDistribution extends DoubleDistribution {

  double mean;
  double stdDev;
  
  ModelEdge meanEdge;
  ModelEdge stdDevEdge;
  ModelEdge varEdge;
  ModelEdge precEdge;
  
  protected NormalDistribution() { }
  
  public NormalDistribution(Model model) {
    this(model, null);
  }
  
  public NormalDistribution(Model model, String name) {
    this(model, name, 0.0, 1.0);
  }

  public NormalDistribution(Model model, double mean, double stdDev) {
    this(model, null, mean, stdDev);
  }
  
  public NormalDistribution(Model model, String name, double mean, double stdDev) {
    super(model, name);
    this.mean = mean;
    this.stdDev = stdDev;
  }
  
  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHNormalProposer(varName);
  }
  
  public NormalDistribution setMean(DoubleValued meanNode) throws MC3KitException {
    meanEdge = updateEdge(meanEdge, (ModelNode)meanNode);
    return this;
  }
  
  public NormalDistribution setVariance(DoubleValued varNode) throws MC3KitException {
    precEdge = updateEdge(precEdge, null);
    stdDevEdge = updateEdge(stdDevEdge, null);
    varEdge = updateEdge(varEdge, (ModelNode)varNode);
    return this;
  }
  
  public NormalDistribution setStdDev(DoubleValued stdDevNode) throws MC3KitException {
    precEdge = updateEdge(precEdge, null);
    varEdge = updateEdge(varEdge, null);
    stdDevEdge = updateEdge(stdDevEdge, (ModelNode)stdDevNode);
    return this;
  }
  
  public NormalDistribution setPrecision(DoubleValued precNode) throws MC3KitException {
    varEdge = updateEdge(varEdge, null);
    stdDevEdge = updateEdge(stdDevEdge, null);
    precEdge = updateEdge(precEdge, (ModelNode)precNode);
    return this;
  }

  @Override
  public double getLogP(Variable v) {
    double x = ((DoubleVariable)v).getValue();
    double mean = meanEdge == null ? this.mean : getDoubleValue(meanEdge);
    
    if(precEdge != null) {
      assert varEdge == null;
      assert stdDevEdge == null;
      return getLogPPrecision(mean, getDoubleValue(precEdge), x);
    }
    else if(varEdge != null) {
      assert precEdge == null;
      assert stdDevEdge == null;
      return getLogPVar(mean, getDoubleValue(varEdge), x);
    }
    else  {
      assert precEdge == null;
      assert varEdge == null;
      return getLogPStdDev(mean, stdDevEdge == null ? this.stdDev : getDoubleValue(stdDevEdge), x);
    }
  }
  
  public static double getLogPPrecision(double mean, double prec, double x) {
    assert prec > 0.0;
    double d = x - mean;
    return 0.5 * (log(prec) - LOG_TWO_PI) - d * d * prec / 2.0;
  }
  
  public static double getLogPStdDev(double mean, double stdDev, double x) {
    assert stdDev > 0.0;
    double d = x - mean;
    return -(log(stdDev) + 0.5 * LOG_TWO_PI) - d * d / (2.0 * stdDev * stdDev);
  }
  
  public static double getLogPVar(double mean, double var, double x) {
    assert var > 0.0;
    double d = x - mean;
    return -0.5 * (log(var) + LOG_TWO_PI) - d * d / (2.0 * var);
  }

  @Override
  public void sample(Variable var) {
    double mean = this.mean;
    if(meanEdge != null) {
      mean = getDoubleValue(meanEdge);
    }
    
    double sd = this.stdDev;
    if(stdDevEdge != null) {
      assert precEdge == null;
      assert varEdge == null;
      sd = getDoubleValue(stdDevEdge);
    }
    else if(varEdge != null) {
      assert precEdge == null;
      assert stdDevEdge == null;
      sd = sqrt(getDoubleValue(varEdge));
    }
    else if(precEdge != null) {
      assert stdDevEdge == null;
      assert varEdge == null;
      sd = sqrt(1.0/getDoubleValue(precEdge));
    }
    
    assert !Double.isNaN(mean);
    assert !Double.isInfinite(mean);
    
    assert sd != 0;
    assert !Double.isNaN(sd);
    assert !Double.isInfinite(sd);
    
    double newVal = new Normal(mean, sd, getRng()).nextDouble();
    assert !Double.isNaN(newVal);
    assert !Double.isInfinite(newVal);
    
    ((DoubleVariable)var).setValue(newVal);
  }

  @Override
  public boolean valueIsValid(double val) {
    if(Double.isInfinite(val) || Double.isNaN(val)) {
      return false;
    }
    return true;
  }
}
