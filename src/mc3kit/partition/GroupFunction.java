package mc3kit.partition;

import mc3kit.Function;
import mc3kit.MC3KitException;
import mc3kit.Model;
import mc3kit.util.IterableBitSet;

@SuppressWarnings("serial")
public class GroupFunction extends Function {
  private PartitionVariable partVar;
  private int groupNum;
  private IterableBitSet group;
  
  protected GroupFunction() { }
  
  public GroupFunction(Model model, PartitionVariable partVar, int groupNum) throws MC3KitException {
    super(model);
    
    this.partVar = partVar;
    this.groupNum = groupNum;
    model.addEdge(this, partVar);
  }
  
  @Override
  public boolean update() {
    IterableBitSet oldGroup = group;
    group = (IterableBitSet)partVar.getGroup(groupNum).clone();
    if(oldGroup == null || !group.equals(oldGroup)) {
      return true;
    }
    return false;
  }
  
  public IterableBitSet getGroup() {
    return group;
  }
}
