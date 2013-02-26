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
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Represents a directed graph.
 * @author Ed Baskerville
 *
 */
@SuppressWarnings("serial")
public class Graph extends Observable implements Serializable {
  Set<Node> nodes;
  Map<String, Node> nodeNameMap;
  
  Map<Node, Integer> nodeOrderMap;
  NavigableMap<Integer, Node> orderNodeMap;
  
  Set<Edge> edges;
  Map<String, Edge> edgeNameMap;

  Map<Node, Set<Edge>> tailNodeMap;
  Map<Node, Set<Edge>> headNodeMap;
  
  transient Logger _logger;
  
  public Graph() {
    nodes = new LinkedHashSet<Node>();
    nodeNameMap = new HashMap<String, Node>();
    
    nodeOrderMap = new HashMap<Node, Integer>();
    orderNodeMap = new TreeMap<Integer, Node>();
    
    edges = new LinkedHashSet<Edge>();
    edgeNameMap = new HashMap<String, Edge>();

    tailNodeMap = new HashMap<Node, Set<Edge>>();
    headNodeMap = new HashMap<Node, Set<Edge>>();
  }
  
  /**
   * Adds a node to the graph.
   * @param node The node to add.
   */
  public Graph addNode(Node node) {
    if(nodes.contains(node)) {
      throw new IllegalArgumentException("Node already in this graph.");
    }
    
    if(node.graph != null) {
      throw new IllegalArgumentException("Node already in another graph.");
    }
    
    if(node.name != null && nodeNameMap.containsKey(node.name)) {
      throw new IllegalArgumentException("Node with this name already in graph.");
    }
    
    node.graph = this;
    nodes.add(node);
    if(node.name != null) {
      nodeNameMap.put(node.name, node);
    }
    tailNodeMap.put(node, new HashSet<Edge>());
    headNodeMap.put(node, new HashSet<Edge>());
    
    setOrder(node, getNextOrder());
    
    return this;
  }
  
  /**
   * Removes a node from the graph.
   * @param node The node to remove.
   */
  public void removeNode(Node node) {
    if(!nodes.contains(node)) {
      throw new IllegalArgumentException("Node not in graph.");
    }
    assert(node.getGraph() == this);
    
    if(!(tailNodeMap.get(node).isEmpty() && headNodeMap.get(node).isEmpty())) {
      throw new IllegalArgumentException("Node has edges in graph.");
    }
    
    removeOrder(node);
    
    nodes.remove(node);
    node.graph = null;
    if(node.name != null) {
      nodeNameMap.remove(node.name);
    }
    tailNodeMap.remove(node);
    headNodeMap.remove(node);
  }
  
  /**
   * Adds an edge to the graph.
   * @param edge The edge to remove. 
   */
  public Graph addEdge(Edge edge) {
    if(edges.contains(edge)) {
      throw new IllegalArgumentException("Edge already in graph.");
    }
    
    if(edge.graph != null) {
      throw new IllegalArgumentException("Edge already in another graph.");
    }
    
    if(edge.tail == null) {
      throw new IllegalArgumentException("Edge tail is null.");
    }
    
    if(edge.head == null) {
      throw new IllegalArgumentException("Edge head is null.");
    }
    
    if(!nodes.contains(edge.tail)) {
      throw new IllegalArgumentException("Edge tail isn't in this graph.");
    }
    
    if(!nodes.contains(edge.head)) {
      throw new IllegalArgumentException("Edge head isn't in this graph.");
    }
    
    if(edge.name != null && edgeNameMap.containsKey(edge.name)) {
      throw new IllegalArgumentException("Edge with this name already in graph.");
    }
    
    // Update order before manipulating other state. If edge introduces a cycle,
    // an exception will be thrown but the graph will be unchanged from its previous state.
    updateOrder(edge);
    
    edge.graph = this;
    edges.add(edge);
    if(edge.name != null) {
      edgeNameMap.put(edge.name, edge);
    }
    
    tailNodeMap.get(edge.tail).add(edge);
    headNodeMap.get(edge.head).add(edge);
    
    return this;
  }
  
  public void removeEdge(Edge edge) {
    if(!edges.contains(edge)) {
      throw new IllegalArgumentException("Edge not in graph.");
    }
    assert(edge.graph == this);
    
    edges.remove(edge);
    edge.graph = null;
    if(edge.name != null) {
      edgeNameMap.remove(edge.name);
    }
    tailNodeMap.get(edge.tail).remove(edge);
    headNodeMap.get(edge.head).remove(edge);
  }
  
  public Set<Edge> getHeadEdges(Node node) {
    return headNodeMap.get(node);
  }
  
  public Set<Edge> getTailEdges(Node node) {
    return tailNodeMap.get(node);
  }
  
  /**
   * Gets a node by name.
   * @param name
   * @return The corresponding node, or null if it does not exist.
   */
  public Node getNode(String name) {
    return nodeNameMap.get(name);
  }
  
  /**
   * Gets an edge by name
   * @param name
   * @return The corresponding edge, or null if it does not exist.
   */
  public Edge getEdge(String name) {
    return edgeNameMap.get(name);
  }
  
  public int nodeCount() {
    return nodes.size();
  }
  
  public int edgeCount() {
    return edges.size();
  }
  
  public Collection<Node> orderedNodesHeadToTail() {
    return orderNodeMap.values();
  }
  
  public int getOrder(Node node) {
    return nodeOrderMap.get(node);
  }
  
  private int getNextOrder() {
    if(nodeOrderMap.isEmpty()) {
      return 0;
    }
    else {
      return orderNodeMap.firstKey() - 1;
    }
  }
  
  private void removeOrder(Node node) {
    assert node != null;
    
    Integer order = nodeOrderMap.remove(node);
    assert order != null;
    
    node = orderNodeMap.remove(order);
    assert node != null;
  }
  
  private void setOrder(Node node, Integer order) {
    assert node != null;
    assert order != null;
    
    nodeOrderMap.put(node, order);
    orderNodeMap.put(order, node);
  }
  
  /*** DYNAMIC MAINTENANCE OF TOPOLOGICAL ORDERING 
   * @throws IllegalArgumentException ***/
  
  private void updateOrder(Edge edge) throws IllegalArgumentException {
    Logger logger = getLogger();
    
    // Implementation of the PK algorithm, published in
    // Pearce, David J. and Paul H. J. Kelly. 2007.
    // A dynamic topological sort algorithm for directed acyclic graphs.
    // ACM Journal of Experimental Algorithmics 11.
    // http://doi.acm.org/10.1145/1187436.1210590
    // http://homepages.ecs.vuw.ac.nz/~djp/dts.html
    
    
    Node tail = edge.getTail();
    Node head = edge.getHead();
    int tailOrder = tail.getOrder();
    int headOrder = head.getOrder();
    assert tailOrder != headOrder;
    if(logger.isLoggable(Level.FINEST)) {
      logger.finest(format("old order: tail %d, head %d", tailOrder, headOrder));
    }
    if(headOrder < tailOrder) {
      // If the head already has a lower order, then do nothing
      return;
    }
    
    // Find affected forward & backward nodes
    // (delta_xy^F, delta_xy^B in paper)
    Set<Node> fwNodes = findAffectedForwardNodes(edge);
    Set<Node> bwNodes = findAffectedBackwardNodes(edge);
    if(logger.isLoggable(Level.FINEST)) {
      logger.finest(format("affected fw: %s", fwNodes));
      logger.finest(format("affected bw: %s", bwNodes));
    }

    // Get sorted list of all indexes
    applyNewIndexes(fwNodes, bwNodes);
    if(logger.isLoggable(Level.FINEST)) {
      logger.finest(format("new order: tail %d, head %d", tail.getOrder(), head.getOrder()));
    }
  }

  private Set<Node> findAffectedForwardNodes(Edge edge) throws IllegalArgumentException {
    Stack<Node> stack = new Stack<Node>();
    Set<Node> visited = new HashSet<Node>();
    Set<Node> outOfOrder = new HashSet<Node>();
    stack.push(edge.tail);
    while(!stack.isEmpty()) {
      Node top = stack.pop();
      
      if(top == edge.head) {
        throw new IllegalArgumentException("Adding edge generates directed cycle in graph."); 
      }
      
      if(!visited.contains(top)) {
        visited.add(top);
        if(top.getOrder() < edge.head.getOrder()) {
          outOfOrder.add(top);
          for(Edge topEdge : getHeadEdges(top)) {
            stack.push(topEdge.getTail());
          }
        }
      }
    }

    return outOfOrder;
  }

  private Set<Node> findAffectedBackwardNodes(Edge edge) {
    Stack<Node> stack = new Stack<Node>();
    Set<Node> visited = new HashSet<Node>();
    Set<Node> outOfOrder = new HashSet<Node>();
    stack.push(edge.head);
    while(!stack.isEmpty()) {
      Node top = stack.pop();
      
      if(top == edge.tail) {
        throw new IllegalArgumentException("Adding edge generates directed cycle in graph."); 
      }
      
      if(!visited.contains(top)) {
        visited.add(top);
        if(top.getOrder() > edge.tail.getOrder()) {
          outOfOrder.add(top);
          for(Edge topEdge : getTailEdges(top)) {
            stack.push(topEdge.getHead());
          }
        }
      }
    }

    return outOfOrder;
  }
  
  private void applyNewIndexes(Set<Node> fwNodes, Set<Node> bwNodes) {
    Logger logger = getLogger();
    
    Set<Node> intersection = new HashSet<Node>(fwNodes);
    intersection.retainAll(bwNodes);
    assert intersection.isEmpty();
    
    // Get sorted list of all indexes
    int[] indexes = new int[fwNodes.size() + bwNodes.size()];
    int i = 0;
    for (Node node : bwNodes)
      indexes[i++] = node.getOrder();
    for (Node node : fwNodes)
      indexes[i++] = node.getOrder();
    Arrays.sort(indexes);
    if(logger.isLoggable(Level.FINEST)) {
      logger.finest(format("indexes: %s", Arrays.toString(indexes)));
    }
    
    Comparator<Node> comparator = new Comparator<Node>() {
      @Override
      public int compare(Node o1, Node o2) {
        return getOrder(o1) - getOrder(o2);
      }
    };
    
    // Get nodes sorted by desired order
    Node[] bwNodesSorted = bwNodes.toArray(new Node[0]);
    Arrays.sort(bwNodesSorted, comparator);
    Node[] fwNodesSorted = fwNodes.toArray(new Node[0]);
    Arrays.sort(fwNodesSorted, comparator);
    
    // Apply indexes to backward nodes, then forward nodes
    i = 0;
    for (Node node : bwNodesSorted)
      setOrder(node, indexes[i++]);
    for (Node node : fwNodesSorted)
      setOrder(node, indexes[i++]);
  }
  
  public Set<Edge> getEdges() {
    return edges;
  }
  
  public boolean verifyOrder() {
    boolean valid = true;
    Set<Node> visited = new HashSet<Node>();
    for(Node node : orderNodeMap.values()) {
      for(Edge edge : tailNodeMap.get(node)) {
        Node head = edge.getHead();
        if(!visited.contains(head)) {
          getLogger().severe(format("node %s: dependency %s not yet visited\n", node, head));
          valid = false;
        }
      }
      
      visited.add(node);
    }
    return valid;
  }
  
  private Logger getLogger() {
    if(_logger == null) {
      _logger = Logger.getLogger("mc3kit.graph.Graph");
    }
    return _logger;
  }
}
