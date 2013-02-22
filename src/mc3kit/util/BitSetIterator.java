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

import java.io.Serializable;
import java.util.BitSet;
import java.util.Iterator;
import java.util.NoSuchElementException;

@SuppressWarnings("serial")
public class BitSetIterator implements Iterator<Integer>, Serializable {
  private BitSet bitSet;
  private int index;
  
  protected BitSetIterator() { }
  
  public BitSetIterator(BitSet bitSet) {
   this.bitSet = bitSet;
   index = -1;
  }

  @Override
  public boolean hasNext() {
    return bitSet.nextSetBit(index + 1) != -1;
  }

  @Override
  public Integer next() {
    index = bitSet.nextSetBit(index + 1);
    if(index == -1) throw new NoSuchElementException("No more elements.");
    
    return index;
  }

  @Override
  public void remove() {
    bitSet.clear(index);
  }
}
