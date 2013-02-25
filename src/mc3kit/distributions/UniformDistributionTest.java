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

package mc3kit.distributions;

import static org.junit.Assert.*;

import mc3kit.DoubleVariable;

import org.junit.Test;

public class UniformDistributionTest {
  
  @Test
  public void minMaxEqual() {
    try {
      new UniformDistribution(null, "ud", 0.0, 0.0);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }

  @Test
  public void minLTMax() {
    try {
      new UniformDistribution(null, "ud", 1.0, 0.0);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }
  
  @Test
  public void minNaN() {
    try {
      new UniformDistribution(null, "ud", Double.NaN, 0.0);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }
  
  @Test
  public void maxNaN() {
    try {
      new UniformDistribution(null, "ud", 0.0, Double.NaN);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }
  
  @Test
  public void minInf() {
    try {
      new UniformDistribution(null, "ud", Double.NEGATIVE_INFINITY, 0.0);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }
  
  @Test
  public void maxInf() {
    try {
      new UniformDistribution(null, "ud", 0.0, Double.POSITIVE_INFINITY);
      fail();
    }
    catch(IllegalArgumentException e) {
    }
  }
  
  @Test
  public void standard() {
    UniformDistribution ud = new UniformDistribution(null);
    DoubleVariable var = new DoubleVariable(null, 0.5);
    double logP = ud.getLogP(var);
    assertEquals(0.0, logP, 1e-12);
  }
  
  @Test
  public void standardOutOfRange() {
    UniformDistribution ud = new UniformDistribution(null);
    assertFalse(ud.valueIsValid(0.0));
    assertFalse(ud.valueIsValid(1.0));
    assertTrue(ud.valueIsValid(1e-12));
    assertTrue(ud.valueIsValid(1.0 - 1e-12));
  }
  
  @Test
  public void severalRandom() {
    double[] mins = { 6.30683091934770, 9.26587566733360, -8.11127731576562 };
    double[] maxs = { 13.9267403981648, 14.6373395924456, -0.1198287634179 };
    double[] xs = { 8.33015870757078, 10.96825034973461, -1.20515368875235 };
    double[] logPs = { -2.03076449020767, -1.68110048312165, -2.07837203900915 };
    
    for(int i = 0; i < mins.length; i++) {
      UniformDistribution ud = new UniformDistribution(null, mins[i], maxs[i]);
      DoubleVariable var = new DoubleVariable(null, xs[i]);
      double logP = ud.getLogP(var);
      assertEquals(logPs[i], logP, 1e-12);
    }
  }
}
