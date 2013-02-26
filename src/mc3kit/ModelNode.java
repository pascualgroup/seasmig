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

package mc3kit;

import java.util.*;
import java.util.logging.Logger;

import cern.jet.random.engine.RandomEngine;
import mc3kit.graph.*;

@SuppressWarnings("serial")
public abstract class ModelNode extends Node {
  Model model;
  
  protected ModelNode() { }

  protected ModelNode(String name) {
    super(name);
  }
  
  public boolean update() {
    return true;
  }
  
  public boolean update(Set<ModelEdge> fromEdges) {
    return update();
  }
  
  public boolean updateAfterRejection() {
    return update();
  }
  
  public boolean updateAfterRejection(Set<ModelEdge> fromEdges) {
    return update(fromEdges);
  }
  
  public Model getModel() {
    return model;
  }
  
  protected ModelEdge updateEdge(ModelEdge edge, ModelNode headNode) throws MC3KitException {
    if(edge != null) {
      if(edge.getHead() == headNode) {
        return edge;
      }
      getModel().removeEdge(edge);
    }
    
    if(headNode == null) {
      return null;
    }
    
    edge = new ModelEdge(getModel(), this, headNode);
    return edge;
  }
  
  protected double getDoubleValue(ModelEdge edge) {
    return ((DoubleValued)edge.getHead()).getValue();
  }
  
  public Chain getChain() {
    return model.getChain();
  }
  
  public RandomEngine getRng() {
    return getChain().getRng();
  }
  
  public Logger getLogger() {
    return model.getLogger();
  }
}
