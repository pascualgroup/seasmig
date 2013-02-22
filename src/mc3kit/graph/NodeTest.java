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

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class NodeTest {
  
  Node node;
  
  @Before
  public void setUp() throws Exception {
    node = new Node();
  }

  @After
  public void tearDown() throws Exception {
  }

  @Test
  public void notInGraph() {
    Node node = new Node();
    
    assertTrue(node.graph == null);
  }
  
  @Test
  public void addToGraph() throws Exception {
    Graph graph = new Graph();
    
    graph.addNode(node);
    assertTrue(node.getGraph() == graph);
  }
  
  @Test
  public void removeFromGraph() throws Exception {
    Graph graph = new Graph();
    
    graph.addNode(node);
    graph.removeNode(node);
    
    assertTrue(node.graph == null);
  }
  
  @Test
  public void alreadyInGraph() throws Exception {
    Graph graph = new Graph();
    graph.addNode(node);
    
    try {
      Graph graph2 = new Graph();
      graph2.addNode(node);
      fail("Should not have been able to add to second graph.");
    }
    catch(Exception e) {
    }
  }
}
