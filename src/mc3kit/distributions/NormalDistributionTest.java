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
import static mc3kit.util.Math.*;

public class NormalDistributionTest {

  @Test
  public void standardAtZero() {
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPStdDev(0.0, 1.0, 0.0), 1e-10);
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPVar(0.0, 1.0, 0.0), 1e-10);
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPPrecision(0.0, 1.0, 0.0), 1e-10);
  }
  
  @Test
  public void standardShiftedAtShift() {
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPStdDev(0.34, 1.0, 0.34), 1e-10);
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPVar(0.34, 1.0, 0.34), 1e-10);
    assertEquals(-0.5 * LOG_TWO_PI, NormalDistribution.getLogPPrecision(0.34, 1.0, 0.34), 1e-10);
  }
  
  @Test
  public void standardScaledAtZero() {
    double logP = -1.61208571376462;
    assertEquals(logP, NormalDistribution.getLogPStdDev(0, 2.0, 0), 1e-10);
    assertEquals(logP, NormalDistribution.getLogPVar(0, 4.0, 0), 1e-10);
    assertEquals(logP, NormalDistribution.getLogPPrecision(0, 0.25, 0), 1e-10);
  }
  
  @Test
  public void standardAtStdDev() {
    double logP = -1.41893853320467;
    assertEquals(logP, NormalDistribution.getLogPStdDev(1.0, 1.0, 0.0), 1e-10);
    assertEquals(logP, NormalDistribution.getLogPVar(1.0, 1.0, 0.0), 1e-10);
    assertEquals(logP, NormalDistribution.getLogPPrecision(1.0, 1.0, 0.0), 1e-10);
  }
  
  @Test
  public void severalRandom() {
    double[] xs = {
         1.312683843659727, 0.670858537743178, 0.194262114674794, -0.899823881697790,
         0.534662494413241, 0.386875294383253, 0.382038536116346,  0.238970471627243,
        -0.744874795402239, -0.214736224175856
    };
    
    double[] means = {
         0.9671728708049787,  0.0346330241961732, -0.8940425209987085,
        -0.6548770527143790, -0.9173025666306984, -0.1425695917268636,
        -0.2158604813894655,  0.8463053340741612,  1.4645341093977491,
         0.3884219976312813
    };
    
    double[] stdDevs = {
        0.3129079407081008, 0.0783979706466198, 0.7181518440646054, 0.3216339857317507,
        0.0680680340155959, 0.5779900937341154, 0.0830698376521468, 1.3251915406015637,
        1.3522130120125535, 0.5493702352978289
    };
    
    double[] logPs = {
        -0.3667141319607241,  -31.3022832310961014, -1.7361194343763757,
        -0.0745919770383709, -225.7392976579421884, -0.7902774297499704,
        -24.3332025845342947,  -1.3055148406771628, -2.5555295676328753,
        -0.9226573285574932
    };
    
    for(int i = 0; i < xs.length; i++) {
      double logP = logPs[i];
      double x = xs[i];
      double mean = means[i];
      double stdDev = stdDevs[i];
      double var = stdDev * stdDev;
      double prec = 1.0 / var;
      
      assertEquals(logP, NormalDistribution.getLogPStdDev(mean, stdDev, x), 1e-10);
      assertEquals(logP, NormalDistribution.getLogPVar(mean, var, x), 1e-10);
      assertEquals(logP, NormalDistribution.getLogPPrecision(mean, prec, x), 1e-10);
    }
  }
}
