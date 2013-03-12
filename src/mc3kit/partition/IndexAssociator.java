package mc3kit.partition;

import java.io.Serializable;

import mc3kit.*;

public interface IndexAssociator extends Serializable {
  void associate(int itemIndex, int groupIndex) throws MC3KitException;
}
