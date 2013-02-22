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
import mc3kit.proposal.GibbsBinaryProposer;
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class BernoulliDistribution extends BinaryDistribution {

  double p;
  ModelEdge pEdge;
  
  protected BernoulliDistribution() { }
  
  public BernoulliDistribution(Model model) {
    this(model, null, 0.5);
  }

  public BernoulliDistribution(Model model, String name) {
    this(model, name, 0.5);
  }
  
  public BernoulliDistribution(Model model, double p) {
    this(model, null, p);
  }
  
  public BernoulliDistribution(Model model, String name, double p) {
    super(model, name);
    this.p = p;
  }
  
  
  public <T extends ModelNode & DoubleValued> BernoulliDistribution setP(T node) throws MC3KitException {
    pEdge = updateEdge(pEdge, node);
    return this;
  }
  
  @Override
  public double getLogP(Variable var) {
    boolean value = ((BinaryVariable)var).getValue();
    double p = pEdge == null ? this.p : getDoubleValue(pEdge);
    
    return value ? log(p) : log1p(-p);
  }

  @Override
  public VariableProposer makeVariableProposer(String varName) {
    return new GibbsBinaryProposer(varName);
  }

  @Override
  public void sample(Variable var) {
    double p = pEdge == null ? this.p : getDoubleValue(pEdge);
    ((BinaryVariable)var).setValue(getRng().nextDouble() < p);
  }
}
