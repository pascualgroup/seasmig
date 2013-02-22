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
import java.util.logging.Logger;

import cern.jet.random.engine.RandomEngine;

@SuppressWarnings("serial")
public class Chain implements Serializable
{
	MCMC mcmc;
	int chainId;
	int chainCount;
	double priorHeatExponent;
	double likelihoodHeatExponent;
	RandomEngine rng;
	
	Model model;
	
	long iterationCount;
	
	protected Chain() { }
	
	Chain(MCMC mcmc, int chainId, int chainCount, double priorHeatExponent, double likelihoodHeatExponent, RandomEngine rng)
	{
		this.mcmc = mcmc;
		this.chainId = chainId;
		this.chainCount = chainCount;
		this.priorHeatExponent = priorHeatExponent;
		this.likelihoodHeatExponent = likelihoodHeatExponent;
		this.rng = rng;
		
		this.iterationCount = 0;
	}
	
	public MCMC getMCMC()
	{
		return mcmc;
	}
  
  void setModel(Model model)
  {
    this.model = model;
  }
	
	public Model getModel()
	{
		return model;
	}
	
	public RandomEngine getRng()
	{
		return rng;
	}
	
	public int getChainId()
	{
		return chainId;
	}
	
	public int getChainCount()
	{
		return chainCount;
	}
	
	public double getPriorHeatExponent()
	{
		return priorHeatExponent;
	}
	
	public double getLikelihoodHeatExponent()
	{
		return likelihoodHeatExponent;
	}
	
	public Chain getChainAbove()
	{
		if(chainId < chainCount - 1)
		{
			return mcmc.getChain(chainId + 1);
		}
		return null;
	}
	
	public Logger getLogger() {
	  return Logger.getLogger("mc3kit.Chain." + chainId);
	}
	
	public long getIterationCount() {
	  return iterationCount;
	}
}
