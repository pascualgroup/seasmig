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

import org.junit.*;
import static org.junit.Assert.*;

public class ExponentialDistributionTest {

  @Test
  public void standardAtOne() {
    assertEquals(-1.0, ExponentialDistribution.getLogPRate(1.0, 1.0), 1e-10);
    assertEquals(-1.0, ExponentialDistribution.getLogPScale(1.0, 1.0), 1e-10);
  }
  
  @Test
  public void standardScaledAtOne() {
    double logP = -1.30685281944005;
    assertEquals(logP, ExponentialDistribution.getLogPRate(1.0, 2.0), 1e-10);
    assertEquals(logP, ExponentialDistribution.getLogPScale(1.0, 0.5), 1e-10);
  }
  
  @Test
  public void severalRandom() {
    double[] xs = {
         1.312683843659727, 0.670858537743178, 0.194262114674794
    };
    
    double[] rates = {
        0.3129079407081008, 0.0783979706466198, 0.7181518440646054
    };
    
    double[] logPs = {
        -1.57259544916731,  -2.59855118452025, -0.470583946194538
    };
    
    for(int i = 0; i < xs.length; i++) {
      double logP = logPs[i];
      double x = xs[i];
      double rate = rates[i];
      double scale = 1.0 / rate;
      
      assertEquals(logP, ExponentialDistribution.getLogPRate(x, rate), 1e-10);
      assertEquals(logP, ExponentialDistribution.getLogPScale(x, scale), 1e-10);
    }
  }
}
