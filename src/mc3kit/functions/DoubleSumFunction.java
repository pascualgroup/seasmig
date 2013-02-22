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

package mc3kit.functions;

import mc3kit.*;

import java.util.*;

@SuppressWarnings("serial")
public class DoubleSumFunction extends DoubleFunction {
  Map<DoubleValued, Summand> summandMap;
  Map<ModelEdge, Summand> edgeMap;
  
  protected DoubleSumFunction() { }
  
  public DoubleSumFunction(Model model) {
    this(model, null);
  }

  public DoubleSumFunction(Model model, String name) {
    super(model, name);
    summandMap = new LinkedHashMap<DoubleValued, Summand>(2);
  }
  
  public <T extends ModelNode & DoubleValued> DoubleSumFunction add(T summandNode) throws ModelException {
    return add(summandNode, 1.0);
  }
  
  public <T extends ModelNode & DoubleValued> DoubleSumFunction add(T summandNode, double coeff) throws ModelException {
    if(summandMap.containsKey(summandNode)) {
      throw new IllegalArgumentException("Summand already present.");
    }
    
    ModelEdge edge = new ModelEdge(getModel(), this, summandNode);
    Summand summand = new Summand(edge, coeff);
    summandMap.put(summandNode, summand);
    
    return this;
  }
  
  public <T extends ModelNode & DoubleValued> DoubleSumFunction remove(T summandNode) throws ModelException {
    Summand summand = summandMap.remove(summandNode);
    if(summand == null) {
      throw new IllegalArgumentException("Summand not present.");
    }
    getModel().removeEdge(summand.edge);
    return this;
  }
  
  @Override
  public boolean update() {
    double oldVal = getValue();
    double newVal = 0.0;
    for(Summand summand : summandMap.values()) {
      newVal += summand.coeff * getDoubleValue(summand.edge);
    }
    if(newVal != oldVal) {
      setValue(newVal);
      return true;
    }
    return false;
  }
  
  private class Summand {
    double coeff;
    ModelEdge edge;
    
    Summand(ModelEdge edge, double coeff) {
      this.coeff = coeff;
      this.edge = edge;
    }
  }
}
