package mc3kit.partition;

import cern.jet.random.Empirical;
import cern.jet.random.EmpiricalWalker;
import cern.jet.random.Gamma;
import cern.jet.random.engine.RandomEngine;
import mc3kit.*;
import static mc3kit.util.Math.*;

@SuppressWarnings("serial")
public class DirichletCategoricalDistribution extends PartitionDistribution {
  
  double alpha;
  ModelEdge alphaEdge;
  
  double[] alphaArray = null;
  ModelEdge alphaArrayEdge;
  
  protected DirichletCategoricalDistribution() { }
  
  public DirichletCategoricalDistribution(Model model) {
    this(model, null, 1.0);
  }
  
  public DirichletCategoricalDistribution(Model model, String name) {
    this(model, name, 1.0);
  }
  
  public DirichletCategoricalDistribution(Model model, double alpha) {
    this(model, null, alpha);
  }
  
  public DirichletCategoricalDistribution(Model model, String name, double alpha) {
    super(model, name);
    this.alpha = alpha;
  }
  
  public DirichletCategoricalDistribution(Model model, String name, double[] alphaArray) {
    super(model, name);
    this.alphaArray = alphaArray;
  }
  
  public DirichletCategoricalDistribution setAlpha(DoubleValued alphaNode) throws MC3KitException {
    alphaEdge = updateEdge(alphaEdge, (ModelNode)alphaNode);
    alphaArrayEdge = updateEdge(alphaArrayEdge, null);
    
    return this;
  }
  
  public DirichletCategoricalDistribution setAlphaArray(DoubleArrayValued alphaArrayNode) throws MC3KitException {
    alphaArrayEdge = updateEdge(alphaArrayEdge, (ModelNode)alphaArrayNode);
    alphaEdge = updateEdge(alphaEdge, null);
    
    return this;
  }
  
  private double getAlpha(int index) {
    if(alphaEdge != null) {
      assert alphaArrayEdge == null;
      
      return ((DoubleValued)alphaEdge.getHead()).getValue();
    }
    else if(alphaArrayEdge != null) {
      assert alphaEdge == null;
      
      return ((DoubleArrayValued)alphaArrayEdge.getHead()).getValue(index);
    }
    else if(alphaArray != null) {
      assert alphaEdge == null;
      assert alphaArrayEdge == null;
      
      return alphaArray[index];
    }
    else {
      return alpha;
    }
  }
  
  @Override
  public double getLogP(Variable var) {
    PartitionVariable partVar = (PartitionVariable)var;
    
    double logP = 0.0;
    
    double sumAlpha = 0.0;
    for(int g = 0; g < partVar.getGroupCount(); g++) {
      int ng = partVar.getGroupSize(g);
      double ag = getAlpha(g);
      sumAlpha += ag;
      logP += logGamma(ng + ag) - logGamma(ag);
    }
    int N = partVar.getElementCount();
    logP += logGamma(sumAlpha) - logGamma(N + sumAlpha);
    
    return logP;
  }
  
  @Override
  public void sample(Variable var) throws MC3KitException {
    PartitionVariable partVar = (PartitionVariable)var;
    
    assert partVar.allowsEmptyGroups == true;
    
    int n = partVar.n;
    int k = partVar.k;
    
    RandomEngine rng = getRng();
    
    // Draw group weights from Dirichlet distribution
    double[] weights = new double[k];
    Gamma gamma = new Gamma(1.0, 1.0, rng);
    for(int g = 0; g < k; g++) {
      weights[g] = gamma.nextDouble(getAlpha(g), 1.0);
    }
    
    // Draw group membership from categorical distribution with Dirichlet-drawn weights
    EmpiricalWalker assigner = new EmpiricalWalker(weights, Empirical.NO_INTERPOLATION, rng);
    int[] gs = new int[n];
    for(int i = 0; i < n; i++) {
      int gi = assigner.nextInt();
      gs[i] = gi;
    }
    
    partVar.setGroups(gs);
  }
}
