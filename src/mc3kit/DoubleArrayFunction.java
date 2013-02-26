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

@SuppressWarnings("serial")
public abstract class DoubleArrayFunction extends Function implements DoubleArrayValued {
  private double[] value;
  
  protected DoubleArrayFunction() { }

  protected DoubleArrayFunction(Model model, int length) {
    super(model);
    this.value = new double[length];
  }

  @Override
  public double getValue(int index) {
    return value[index];
  }

  @Override
  public void setValue(int index, double value) {
    this.value[index] = value;
  }
  
  @Override
  public double[] getValue() {
    return value.clone();
  }
  
  @Override
  public void setValue(double[] value) {
    if(value.length != this.value.length) {
      throw new IllegalArgumentException("Double array functions have fixed length.");
    }
    
    for(int i = 0; i < value.length; i++) {
      this.value[i] = value[i];
    }
  }
  
  @Override
  public int getLength() {
    return value.length;
  }
}
