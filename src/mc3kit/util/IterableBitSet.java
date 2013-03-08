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

package mc3kit.util;

import java.util.*;

@SuppressWarnings("serial")
public class IterableBitSet extends BitSet implements Iterable<Integer> {
  
  public IterableBitSet() { }

  public IterableBitSet(int n) {
    super(n);
  }
  
  public int[] getSetBits() {
    int[] setBits = new int[cardinality()];
    int j = 0;
    for(int i = nextSetBit(0); i >= 0; i = nextSetBit(i+1)) {
      setBits[j++] = i;
    }
    return setBits;
  }
  
  @Override
  public Iterator<Integer> iterator() {
    return new BitSetIterator(this);
  }
  
  @Override
  public Object clone() {
    IterableBitSet bs = new IterableBitSet();
    for(int i : this) {
      bs.set(i);
    }
    return bs;
  }
}
