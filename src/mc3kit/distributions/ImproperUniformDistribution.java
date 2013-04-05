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
import mc3kit.proposal.MHNormalProposer;

@SuppressWarnings("serial")
public class ImproperUniformDistribution extends DoubleDistribution {
  protected ImproperUniformDistribution() { }
  
  public ImproperUniformDistribution(Model model) {
    this(model, null);
  }

  public ImproperUniformDistribution(Model model, String name) {
    super(model, name);
  }

  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new MHNormalProposer(varName);
  }

  @Override
  public boolean valueIsValid(double value) {
    return !Double.isNaN(value) && !Double.isInfinite(value);
  }

  @Override
  public double getLogP(Variable var) {
    return 0.0;
  }

  @Override
  public void sample(Variable var) throws MC3KitException {
    throw new MC3KitException("Cannot sample from improper prior. Initialize value during construction.");
  }

  @Override
  public boolean update() {
    return false;
  }
}
