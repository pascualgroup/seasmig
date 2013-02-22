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
	double likeHeatPower = 1.0;
  double minLikeHeatExponent = 0.0;
  
  double priorHeatPower = 1.0;
  double minPriorHeatExponent = 1.0;
  
  protected PowerHeatFunction() { }
  
  public PowerHeatFunction(double likeHeatPower, double minLikeHeatExponent, double priorHeatPower, double minPriorHeatExponent) {
    this.likeHeatPower = likeHeatPower;
    this.minLikeHeatExponent = minLikeHeatExponent;
    this.priorHeatPower = priorHeatPower;
    this.minPriorHeatExponent = minPriorHeatExponent;
  }
  
  public PowerHeatFunction(double heatPower, double minHeatExponent) {
    this.likeHeatPower = heatPower;
    this.minLikeHeatExponent = minHeatExponent;
  }
  
	public double getHeatPower()
	{
		return likeHeatPower;
	}

	public void setHeatPower(double heatPower)
	{
		this.likeHeatPower = heatPower;
	}

	public double getMinHeatExponent()
	{
		return minLikeHeatExponent;
	}
	
	public void setMinHeatExponent(double minHeatExponent)
	{
		this.minLikeHeatExponent = minHeatExponent;
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
        heatExponents[i] = minPriorHeatExponent + pow(linearValue, priorHeatPower) * (1.0 - minPriorHeatExponent);
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
				heatExponents[i] = minLikeHeatExponent + pow(linearValue, likeHeatPower) * (1.0 - minLikeHeatExponent);
			}
		}
			
		return heatExponents;
	}

}
