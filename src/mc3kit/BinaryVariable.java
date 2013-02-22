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
public class BinaryVariable extends Variable implements BinaryValued {
  
  private boolean value;
  
  protected BinaryVariable() { }
  
  public BinaryVariable(Model model, boolean value) {
    this(model, (String)null, value);
  }
  
  public BinaryVariable(Model model, String name) {
    super(model, name, false);
  }
  
  public BinaryVariable(Model model, String name, boolean value) {
    super(model, name, true);
    this.value = value;
  }
  
  public BinaryVariable(Model model, boolean value, BinaryDistribution dist) throws MC3KitException {
    this(model, null, value, dist);
  }
  
  public BinaryVariable(Model model, String name, boolean value, BinaryDistribution dist) throws MC3KitException {
    super(model, name, true, dist);
    this.value = value;
  }
  
  public BinaryVariable(Model model, String name, BinaryDistribution dist) throws MC3KitException {
    super(model, name, false, dist);
  }

  @Override
  public boolean getValue() {
    return value;
  }

  @Override
  public void setValue(boolean value) {
    if(isObserved()) {
      throw new UnsupportedOperationException("Can't set value on an observed variable.");
    }
    
    this.value = value;
    setChanged();
    notifyObservers();
  }
  
  @Override
  public BinaryVariable setDistribution(Distribution dist) throws MC3KitException {
    super.setDistribution(dist);
    return this;
  }

  @Override
  public boolean update() {
    BinaryDistribution dist = (BinaryDistribution)getDistribution();
    setLogP(dist.getLogP(this));
    return false;
  }

  @Override
  public Object makeOutputObject() {
    return value ? 1 : 0;
  }
  
  @Override
  public String makeOutputString() {
    return value ? "1" : "0";
  }
}
