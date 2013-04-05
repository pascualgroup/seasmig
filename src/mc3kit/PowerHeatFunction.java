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

import static java.lang.Math.pow;

import java.io.Serializable;

@SuppressWarnings("serial")
public class PowerHeatFunction implements HeatFunction, Serializable
{
	double likeHeatPower;
  double minLikeHeatExp;
  double maxLikeHeatExp;
  
  double priorHeatPower;
  double minPriorHeatExp;
  double maxPriorHeatExp;
  
  protected PowerHeatFunction() { }
  
  public PowerHeatFunction(double heatPower, double minHeatExponent) {
    this(heatPower, minHeatExponent, 1.0, 1.0);
  }
  
  public PowerHeatFunction(double likeHeatPow, double minLikeHeatExp, double priorHeatPow, double minPriorHeatExp) {
    this(likeHeatPow, minLikeHeatExp, 1.0, priorHeatPow, minPriorHeatExp, 1.0);
  }
  
  public PowerHeatFunction(double likeHeatPow, double minLikeHeatExp, double maxLikeHeatExp,
    double priorHeatPow, double minPriorHeatExp, double maxPriorHeatExp) {
    
    this.likeHeatPower = likeHeatPow;
    this.minLikeHeatExp = minLikeHeatExp;
    this.maxLikeHeatExp = maxLikeHeatExp;
    
    this.priorHeatPower = priorHeatPow;
    this.minPriorHeatExp = minPriorHeatExp;
    this.maxPriorHeatExp = maxPriorHeatExp;
  }
	
	@Override
	public double[] getPriorHeatExponents(int chainCount)
	{
    double[] heatExponents = new double[chainCount];
    heatExponents[0] = 1.0;
    if(chainCount > 1)
    {
      for(int i = 1; i < chainCount; i++)
      {
        double linearValue = 1.0 - i / ((double)chainCount-1);
        heatExponents[i] = minPriorHeatExp + pow(linearValue, priorHeatPower) * (maxPriorHeatExp - minPriorHeatExp);
      }
    }
      
    return heatExponents;
	}

	@Override
	public double[] getLikelihoodHeatExponents(int chainCount)
	{
		double[] heatExponents = new double[chainCount];
		heatExponents[0] = 1.0;
		if(chainCount > 1)
		{
			for(int i = 1; i < chainCount; i++)
			{
				double linearValue = 1.0 - i / ((double)chainCount-1);
				heatExponents[i] = minLikeHeatExp + pow(linearValue, likeHeatPower) * (maxLikeHeatExp - minLikeHeatExp);
			}
		}
			
		return heatExponents;
	}

}
