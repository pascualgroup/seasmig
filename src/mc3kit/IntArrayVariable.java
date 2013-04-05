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

import java.util.Arrays;

@SuppressWarnings("serial")
public class IntArrayVariable extends Variable implements IntArrayValued {
  
  private int[] value;
  
  protected IntArrayVariable() { }
  
  /**
   * Constructor for an <i>observed</i> int-array-valued variable.
   * @param model
   * @param value
   */
  public IntArrayVariable(Model model, int[] value) {
    this(model, (String)null, value);
  }
  
  /**
   * Constructor for an <i>observed</i> int-array-valued variable with a name.
   * @param model
   * @param name
   * @param value
   */
  public IntArrayVariable(Model model, String name, int[] value) {
    super(model, name, true);
    this.value = value;
  }
  
  /**
   * Constructor for an <i>observed</i> int-array-valued variable with a distribution.
   * @param model
   * @param name
   * @param value
   * @throws MC3KitException 
   */
  public IntArrayVariable(Model model, int[] value, IntArrayDistribution dist) throws MC3KitException {
    this(model, null, value, dist);
  }
  
  /**
   * Constructor for an <i>observed</i> int-array-valued variable with a name and distribution.
   * @param model
   * @param name
   * @param value
   * @param dist
   * @throws MC3KitException
   */
  public IntArrayVariable(Model model, String name, int[] value, IntArrayDistribution dist) throws MC3KitException {
    super(model, name, true, dist);
    this.value = value;
  }

  @Override
  public int[] getValue() {
    return value;
  }
  
  public IntArrayVariable setDistribution(IntArrayDistribution dist) throws MC3KitException {
    super.setDistribution(dist);
    return this;
  }
  
  @Override
  public Object makeOutputObject() {
    return value;
  }
  
  @Override
  public String makeOutputString() {
    return Arrays.toString(value);
  }

  @Override
  public int getLength() {
    return value.length;
  }

  @Override
  public void setValue(int[] value) {
    if(value.length != this.value.length) {
      throw new IllegalArgumentException("int array variables have fixed length");
    }
    
    for(int i = 0; i < value.length; i++) {
      this.value[i] = value[i];
    }
  }

  @Override
  public int getValue(int index) {
    return value[index];
  }

  @Override
  public void setValue(int index, int value) {
    this.value[index] = value;
  }
}
