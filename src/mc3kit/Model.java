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
import java.util.logging.*;
import java.io.*;

import static java.lang.String.format;
import static java.lang.Math.*;
import static mc3kit.util.Utils.*;

import cern.jet.random.engine.RandomEngine;

import mc3kit.graph.*;

/**
 * Represents a directed probabilistic graphical model consisting of random
 * variables, distributions, and (deterministic) functions in a directed
 * acyclic graph.
 * @author Ed Baskerville
 *
 */
@SuppressWarnings("serial")
public class Model implements Observer, Serializable {
  
  Chain chain;
  Graph graph;
  
  List<Variable> metaVariables;
  List<Variable> unobservedVariables;
  
  private double logPrior;
  private double logLikelihood;
  
  private double oldLogPrior;
  private double oldLogLikelihood;
  
  State state;
  Set<Variable> changedValueVars;
  Set<ModelNode> newEdgeHeads;
  
  protected Model() { }
  
  public Model(Chain initialChain) {
    this.chain = initialChain;
    graph = new Graph();
    metaVariables = new ArrayList<Variable>();
    unobservedVariables = new ArrayList<Variable>();
    state = State.UNINITIALIZED;
    changedValueVars = new HashSet<Variable>();
    newEdgeHeads = new HashSet<ModelNode>();
  }
  
  private void readObject(java.io.ObjectInputStream stream) throws IOException, ClassNotFoundException {
    stream.defaultReadObject();
    for(Variable var : unobservedVariables) {
      var.addObserver(this);
    }
  }
  
  public List<Variable> getUnobservedVariables() {
    return unobservedVariables;
  }
  
  public String[] getUnobservedVariableNames() {
    String[] varNames = new String[unobservedVariables.size()];
    for(int i = 0; i < varNames.length; i++) {
      varNames[i] = unobservedVariables.get(i).getName();
    }
    return varNames;
  }
  
  /*** CALCULATIONS ***/
  
  public void beginConstruction() throws MC3KitException {
    if(state != State.UNINITIALIZED) {
      throw new MC3KitException("beginConstruction called with wrong state");
    }
    
    state = State.IN_CONSTRUCTION;
  }
  
  public void endConstruction() throws MC3KitException {
    Logger logger = getLogger();
    
    logPrior = 0.0;
    logLikelihood = 0.0;
    
    if(state != State.IN_CONSTRUCTION) {
      throw new MC3KitException("endConstruction called with wrong state");
    }
    
    if(logger.isLoggable(Level.FINE)) {
      logger.fine("NODE ORDER:");
      for(Node node : graph.orderedNodesHeadToTail()) {
        if(node instanceof Variable) {
          logger.fine(format("  %s", node.getName()));
        }
      }
      
      logger.fine("EDGES:");
      for(Edge edge : graph.getEdges()) {
        logger.fine(format("  %s -> %s", edge.getTail(), edge.getHead()));
      }
    }
    
    assert graph.verifyOrder();
    
    // Sample from meta-graph: variables that can manipulate the graph and variables they depend on
    logger.fine("Sampling from meta-graph...");
    Collection<Node> metaGraphNodes = graph.orderedNodesHeadToTails(metaVariables);
    Set<Node> metaGraphNodeSet = new HashSet<Node>(metaGraphNodes);
    for(Node node : metaGraphNodes) {
      sampleNode(node);
    }
    
    // Sample from full graph, ignoring nodes in meta-graph
    logger.fine("Sampling from rest of graph...");
    for(Node node : graph.orderedNodesHeadToTail()) {
      if(!metaGraphNodeSet.contains(node)) {
        sampleNode(node);
      }
    }
    
    changedValueVars.clear();
    newEdgeHeads.clear();
    
    state = State.READY;
  }
  
  private void sampleNode(Node node) throws MC3KitException {
    Logger logger = getLogger();
    
    if(node instanceof Variable) {
      Variable var = (Variable)node;
      if(!var.isObserved() && !changedValueVars.contains(var)) {
        var.sample();
        if(logger.isLoggable(Level.FINE)) {
          logger.fine(format("Sampling %s: %s", var, var.makeOutputString()));
        }
      }
      else if(!var.isObserved()) {
        if(logger.isLoggable(Level.FINE)) {
          logger.fine(format("Not sampling %s", var));
        }
      }
    }
    
    ((ModelNode)node).update();
    if(logger.isLoggable(Level.FINE)) {
      logger.fine(format("Updating %s", node));
    }
    
    if(node instanceof Variable) {
      Variable var = (Variable)node;
      if(var.isObserved()) {
        logLikelihood += var.getLogP();
      }
      else {
        logPrior += var.getLogP();
      }
    }
  }
  
  public void recalculate(double tol) throws MC3KitException {
    double preLogPrior = logPrior;
    double preLogLike = logLikelihood;
    
    logPrior = 0.0;
    logLikelihood = 0.0;
    
    for(Node node : graph.orderedNodesHeadToTail()) {
      ((ModelNode)node).update();
      
      if(node instanceof Variable) {
        Variable var = (Variable)node;
        if(var.isObserved()) {
          logLikelihood += var.getLogP();
        }
        else {
          logPrior += var.getLogP();
        }
      }
    }
    
    double logPriorDiff = abs(preLogPrior - logPrior);
    double logLikeDiff = abs(preLogLike - logLikelihood);
    
    if(logPriorDiff > tol || logLikeDiff > tol) {
      throw new MC3KitException(format("Too much error in prior (%f, should be %f, diff %f) or likelihood (%f, should be %f, diff %f)", 
          preLogPrior, logPrior, logPriorDiff, preLogLike, logLikelihood, logLikeDiff)
      );
    }
  }
  
  public void beginProposal() throws MC3KitException {
    if(state != State.READY) {
      throw new MC3KitException("beginProposal called with wrong state");
    }

    oldLogPrior = logPrior;
    oldLogLikelihood = logLikelihood;
    state = State.IN_PROPOSAL;
  }
  
  public void endProposal() throws MC3KitException {
    if(state != State.IN_PROPOSAL) {
      throw new MC3KitException("endProposal called with wrong state");
    }
    
    propagateChanges(true);
    
    state = State.PROPOSAL_COMPLETE;
  }
  
  public void acceptProposal() throws MC3KitException {
    if(state != State.PROPOSAL_COMPLETE) {
      throw new MC3KitException("acceptProposal called with wrong state");
    }
    state = State.READY;
  }
  
  public void beginRejection() throws MC3KitException {
    if(state != State.PROPOSAL_COMPLETE) {
      throw new MC3KitException("beginRejection called with wrong state");
    }
    
    state = State.IN_REJECTION;
  }
  
  public void endRejection() throws MC3KitException {
    if(state != State.IN_REJECTION) {
      throw new MC3KitException("endRejection called with wrong state");
    }
    
    propagateChanges(false);
    logPrior = oldLogPrior;
    logLikelihood = oldLogLikelihood;
    
    state = State.READY;
  }
  
  private void propagateChanges(boolean isProposal) throws MC3KitException {
    Set<ModelEdge> emptyEdgeSet = new HashSet<ModelEdge>(0);
    
    // Map of nodes to parent edges that have changed, for the sake of efficient updating
    Map<ModelNode, Set<ModelEdge>> visitedEdges = new HashMap<ModelNode, Set<ModelEdge>>();
    
    // Queue of nodes to update in topological order, starting with variables
    // whose values have changed and heads of new edges
    SortedMap<Integer, ModelNode> updateQueue = new TreeMap<Integer, ModelNode>();
    for(Variable var : changedValueVars) {
      updateQueue.put(var.getOrder(), var);
    }
    for(ModelNode node : newEdgeHeads) {
      updateQueue.put(node.getOrder(), node);
    }
    
    Set<ModelNode> visitedNodes = new HashSet<ModelNode>();
    
    // Traverse graph in topological order
    int lastOrder = Integer.MIN_VALUE;
    while(!updateQueue.isEmpty()) {
      int order = updateQueue.firstKey();
      assert order > lastOrder;
      lastOrder = order;
      
      ModelNode node = updateQueue.get(order);
      visitedNodes.add(node);
      updateQueue.remove(order);
      
      Set<ModelEdge> fromEdges = visitedEdges.get(node);
      if(fromEdges == null) fromEdges = emptyEdgeSet;
      
      // Get old log-probability for random variables
      double oldLogP = 0.0;
      if(node instanceof Variable) {
        oldLogP = ((Variable)node).getLogP();
      }
      
      // Update the node
      boolean changed;
      if(isProposal) {
        changed = node.update(fromEdges);
      }
      else {
        changed = node.updateAfterRejection(fromEdges);
      }
      
      // Update log-prior or log-likelihood for random variables
      if(node instanceof Variable) {
        Variable var = (Variable)node;
        double newLogP = var.getLogP();
        if(var.isObserved()) {
          logLikelihood += (newLogP - oldLogP);
        }
        else {
          logPrior += (newLogP - oldLogP);
        }
      }
      
      // If the node changed, add all the edges representing dependencies
      // on it (the head) to the visited-edge map for the dependent nodes
      // (the tails)
      if(changed || changedValueVars.contains(node) || newEdgeHeads.contains(node)) {
        for(Edge edge : node.getHeadEdges()) {
          ModelNode tail = (ModelNode)edge.getTail();
          if(!visitedEdges.containsKey(tail)) {
            visitedEdges.put(tail, new HashSet<ModelEdge>());
          }
          visitedEdges.get(tail).add((ModelEdge)edge);
          assert !visitedNodes.contains(tail);
          updateQueue.put(tail.getOrder(), tail);
        }
      }
    }
    
    changedValueVars.clear();
    newEdgeHeads.clear();
  }
  
  /*** GRAPH CONSTRUCTION/MANIPULATION ***/
  
  public <T extends ModelNode> T addNode(T node) {
    if(node instanceof Variable) {
      addVariable((Variable)node);
    }
    else if(node instanceof Function) {
      addFunction((Function)node);
    }
    else if(node instanceof Distribution) {
      addDistribution((Distribution)node);
    }
    else {
      throw new IllegalArgumentException("Unknown node type.");
    }
    return node;
  }
  
  public <V extends Variable> V addVariable(V var) {
    if(var.model != null) {
      throw new IllegalArgumentException("Variable already in model");
    }
    
    graph.addNode(var);
    var.model = this;
    
    if(!var.isObserved()) {
      if(var.getName() == null) {
        throw new IllegalArgumentException("Unobserved random variables must have a name.");
      }
      unobservedVariables.add(var);
      var.addObserver(this);
      
      if(var.canManipulateGraph()) {
        metaVariables.add(var);
      }
    }
    
    return var;
  }
  
  public <F extends Function> F addFunction(F func) {
    if(func.model != null) {
      throw new IllegalArgumentException("Function already in model.");
    }
    
    graph.addNode(func);
    func.model = this;
    
    return func;
  }
  
  public <D extends Distribution> D addDistribution(D dist) {
    if(dist.model != null) {
      throw new IllegalArgumentException("Distribution already in model.");
    }
    
    graph.addNode(dist);
    dist.model = this;
    
    return dist;
  }
  
  public ModelEdge addEdge(ModelNode tail, ModelNode head) throws MC3KitException {
    if(!(state == State.IN_CONSTRUCTION || state == State.IN_PROPOSAL || state == State.IN_REJECTION)) {
      throw new MC3KitException("Adding edge in wrong state");
    }
    
    ModelEdge edge = new ModelEdge(this, tail, head);
    return edge;
  }
  
  public void addEdge(ModelEdge edge) throws MC3KitException {
    if(!(state == State.IN_CONSTRUCTION || state == State.IN_PROPOSAL || state == State.IN_REJECTION)) {
      throw new MC3KitException("Adding edge in wrong state");
    }
    
    ModelNode head = edge.getHead();
    if(head instanceof Variable) {
      assert !(((Variable)head).canManipulateGraph());
    }
    
    graph.addEdge(edge);
    newEdgeHeads.add(edge.getHead());
  }
  
  public void removeEdge(ModelEdge edge) throws MC3KitException {
    if(!(state == State.IN_CONSTRUCTION || state == State.IN_PROPOSAL || state == State.IN_REJECTION)) {
      throw new MC3KitException("Removing edge in wrong state");
    }
    
    graph.removeEdge(edge);
  }
  
  public Map<String, Object> makeHierarchicalSample() {
    Map<String, Object> flatMap = new LinkedHashMap<String, Object>();
    
    flatMap.put("iterationCount", getChain().getIterationCount() + 1);
    flatMap.put("logPrior", logPrior);
    flatMap.put("logLikelihood", logLikelihood);
    
    for(Variable var : unobservedVariables) {
      flatMap.put(var.getName(), var.makeOutputObject());
    }
    
    return makeHierarchicalMap(flatMap);
  }
  
  public Map<String, String> makeFlatSample() {
    Map<String, String> samp = new LinkedHashMap<String, String>();
    
    samp.put("iterationCount", Long.toString(getChain().getIterationCount()));
    samp.put("logPrior", Double.toString(logPrior));
    samp.put("logLikelihood", Double.toString(logLikelihood));
    
    for(Variable var : unobservedVariables) {
      samp.put(var.getName(), var.makeOutputString());
    }
    
    return samp;
  }
  
  /*** GETTERS ***/
  
  public double getLogPrior() {
    return logPrior;
  }
  
  public double getLogLikelihood() {
    return logLikelihood;
  }
  
  public ModelNode get(String name) {
    return (ModelNode)graph.getNode(name);
  }
  
  public Variable getVariable(String name) {
    return (Variable)graph.getNode(name);
  }
  
  public BinaryVariable getBinaryVariable(String name) {
    return (BinaryVariable)graph.getNode(name);
  }
  
  public DoubleVariable getDoubleVariable(String name) {
    return (DoubleVariable)graph.getNode(name);
  }
  
  public DoubleFunction getDoubleFunction(String name) {
    return (DoubleFunction)graph.getNode(name);
  }
  
  public Distribution getDistribution(String name) {
    return (Distribution)graph.getNode(name);
  }
  
  private enum State {
    UNINITIALIZED,
    IN_CONSTRUCTION,
    READY,
    IN_PROPOSAL,
    PROPOSAL_COMPLETE,
    IN_REJECTION
  }

  @Override
  public void update(Observable obj, Object arg1) {
    if(obj instanceof Variable) {
      Variable var = (Variable)obj;
      assert state == State.IN_CONSTRUCTION || state == State.IN_PROPOSAL || state == State.IN_REJECTION;
      if(state == State.IN_PROPOSAL || state == State.IN_REJECTION) {
        assert !var.isObserved();
      }
      
      if(!var.isObserved()) {
        changedValueVars.add(var);
      }
    }
  }
  
  public void setChain(Chain chain) {
    this.chain = chain;
  }
  
  public Chain getChain() {
    return chain;
  }
  
  public RandomEngine getRng() {
    return chain.getRng();
  }
  
  public Logger getLogger() {
    return chain.getLogger();
  }
}
