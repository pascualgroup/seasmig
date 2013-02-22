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

package mc3kit.graph;

import java.io.Serializable;
import java.util.Observable;

@SuppressWarnings("serial")
public class Edge extends Observable implements Serializable {
  Node tail;
  Node head;
  String name;
  Graph graph;
  
  protected Edge() { }
  
  public Edge(Node tail, Node head) {
    this.tail = tail;
    this.head = head;
  }
  
  public Edge(String name, Node tail, Node head) {
    this(tail, head);
    this.name = name;
  }
  
  public Node getTail() {
    return tail;
  }
  
  public Node getHead() {
    return head;
  }
}
