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

public class GammaDistributionTest {

  @Test
  public void standardExpAtOne() {
    assertEquals(-1.0, GammaDistribution.getLogPRate(1.0, 1.0, 1.0), 1e-10);
    assertEquals(-1.0, GammaDistribution.getLogPScale(1.0, 1.0, 1.0), 1e-10);
  }
  
  @Test
  public void standardExpScaledAtOne() {
    double logP = -1.30685281944005;
    assertEquals(logP, GammaDistribution.getLogPRate(1.0, 1.0, 2.0), 1e-10);
    assertEquals(logP, GammaDistribution.getLogPScale(1.0, 1.0, 0.5), 1e-10);
  }
  
  @Test
  public void gamma22AtOne() {
    double logP = -0.613705638880109;
    assertEquals(logP, GammaDistribution.getLogPRate(1.0, 2.0, 2.0), 1e-10);
    assertEquals(logP, GammaDistribution.getLogPScale(1.0, 2.0, 0.5), 1e-10);
  }
  
  @Test
  public void gammaHalfHalfAtOne() {
    double logP = -1.41893853320467;
    assertEquals(logP, GammaDistribution.getLogPRate(1.0, 0.5, 0.5), 1e-10);
    assertEquals(logP, GammaDistribution.getLogPScale(1.0, 0.5, 2.0), 1e-10);
  }
  
  @Test
  public void gamma2HalfAtOne() {
    double logP = -1.88629436111989;
    assertEquals(logP, GammaDistribution.getLogPRate(1.0, 2.0, 0.5), 1e-10);
    assertEquals(logP, GammaDistribution.getLogPScale(1.0, 2.0, 2.0), 1e-10);
  }
  
  @Test
  public void gammaHalfTwoAtOne() {
    double logP = -2.22579135264473;
    assertEquals(logP, GammaDistribution.getLogPRate(1.0, 0.5, 2.0), 1e-10);
    assertEquals(logP, GammaDistribution.getLogPScale(1.0, 0.5, 0.5), 1e-10);
  }
  
  @Test
  public void severalExpRandom() {
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
      
      assertEquals(logP, GammaDistribution.getLogPRate(x, 1.0, rate), 1e-10);
      assertEquals(logP, GammaDistribution.getLogPScale(x, 1.0, scale), 1e-10);
    }
  }
  
  @Test
  public void severalGammaRandom() {
    double[] xs = {
        0.165520979559969, 1.744962668356796, 3.419642100591502
    };
    
    double[] shapes = {
        4.3757828851382348, 0.0159208724896113, 2.1933361329138279
    };
    
    double[] rates = {
        3.1516953440264173, 0.1952994993719014, 0.0416127828806077
    };
    
    double[] logPs = {
        -3.85356377761384, -5.04580351596200, -5.74176289694918
    };
    
    for(int i = 0; i < xs.length; i++) {
      double logP = logPs[i];
      double x = xs[i];
      double shape = shapes[i];
      double rate = rates[i];
      double scale = 1.0 / rate;
      
      assertEquals(logP, GammaDistribution.getLogPRate(x, shape, rate), 1e-10);
      assertEquals(logP, GammaDistribution.getLogPScale(x, shape, scale), 1e-10);
    }
  }
}
