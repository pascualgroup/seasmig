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
public class IntVariable extends Variable implements IntValued {
  
  private int value;
  
  protected IntVariable() { }
  
  public IntVariable(Model model, int value) {
    this(model, (String)null, value);
  }
  
  public IntVariable(Model model, String name) {
    super(model, name, false);
  }
  
  public IntVariable(Model model, String name, int value) {
    super(model, name, true);
    this.value = value;
  }
  
  public IntVariable(Model model, int value, IntDistribution dist) throws ModelException {
    this(model, null, value, dist);
  }
  
  public IntVariable(Model model, String name, int value, IntDistribution dist) throws ModelException {
    super(model, name, true, dist);
    this.value = value;
  }
  
  public IntVariable(Model model, String name, IntDistribution dist) throws ModelException {
    super(model, name, false, dist);
  }

  @Override
  public int getValue() {
    return value;
  }

  @Override
  public void setValue(int value) {
    if(isObserved()) {
      throw new UnsupportedOperationException("Can't set value on an observed variable.");
    }
    
    this.value = value;
    setChanged();
    notifyObservers();
  }
  
  @Override
  public IntVariable setDistribution(Distribution dist) throws ModelException {
    super.setDistribution(dist);
    return this;
  }

  @Override
  public boolean update() {
    IntDistribution dist = (IntDistribution)getDistribution();
    setLogP(dist.getLogP(this));
    return false;
  }

  @Override
  public Object makeOutputObject() {
    return value;
  }
  
  @Override
  public String makeOutputString() {
    return Integer.toString(value);
  }
}
