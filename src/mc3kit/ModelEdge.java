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

import mc3kit.graph.Edge;

@SuppressWarnings("serial")
public class ModelEdge extends Edge {
  Model model;
  
  protected ModelEdge() { }
  
  public ModelEdge(Model model, ModelNode tail, ModelNode head) throws ModelException {
    super(tail, head);
    model.addEdge(this);
  }

  public ModelEdge(Model model, String name, ModelNode tail, ModelNode head) throws ModelException {
    super(name, tail, head);
    model.addEdge(this);
  }
  
  public ModelNode getTail() {
    return (ModelNode)super.getTail();
  }
  
  public ModelNode getHead() {
    return (ModelNode)super.getHead();
  }
}
