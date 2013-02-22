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

package mc3kit.example;

import mc3kit.*;

@SuppressWarnings("serial")
public class ExampleModelFactory implements ModelFactory {
  double[] data;
  
  public ExampleModelFactory(double[] data)
  {
    this.data = data;
  }

  @Override
  public Model createModel(Chain initialChain) throws MC3KitException {
    // Models can also be created without subclassing directly in this method,
    // although subclassing makes it more convenient to keep track of
    // objects within the model. E.g., you could just do
    // 
    // Model m = new Model(initialChain);
    // m.beginConstruction();
    // DoubleDistribution d = new NormalDistribution(m);
    // DoubleVariable v = new DoubleVariable(m, "v", d);
    // v.setValue(-10 + 20 * m.getRng().nextDouble());
    // m.endConstruction();
    // return m;
    //
    // Construction must be bounded by beginConstruction() and endConstruction().
    
    return new ExampleModel(initialChain, data);
  }
}
