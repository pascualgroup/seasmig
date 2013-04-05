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

import java.io.Serializable;

@SuppressWarnings("serial")
public class ConstantHeatFunction implements HeatFunction, Serializable
{
  double priorHeatExponent;
  double likelihoodHeatExponent;
  
  public ConstantHeatFunction() {
    this(1.0, 1.0);
  }
  
  public ConstantHeatFunction(double priorHeatExponent, double likelihoodHeatExponent) {
    this.priorHeatExponent = priorHeatExponent;
    this.likelihoodHeatExponent = likelihoodHeatExponent;
  }
  
	@Override
	public double[] getPriorHeatExponents(int chainCount)
	{
		double[] heatExponents = new double[chainCount];
		for(int i = 0; i < chainCount; i++)
			heatExponents[i] = priorHeatExponent;
		return heatExponents;
	}

	@Override
	public double[] getLikelihoodHeatExponents(int chainCount)
	{
	  System.err.printf("likelihoodheatexpnoent: %f\n", likelihoodHeatExponent);
	  
    double[] heatExponents = new double[chainCount];
    for(int i = 0; i < chainCount; i++)
      heatExponents[i] = likelihoodHeatExponent;
    return heatExponents;
	}

}
