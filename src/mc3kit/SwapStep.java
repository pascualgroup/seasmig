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

import static java.lang.Math.log;

import java.io.*;
import java.util.*;
import static java.lang.String.*;

import cern.jet.random.engine.RandomEngine;
import mc3kit.util.Collector;

@SuppressWarnings("serial")
public class SwapStep implements Step
{
  SwapParity swapParity;
  long statsEvery = 100;
  
	/*** ENUMS ***/
	
	public static enum SwapParity
	{
		EVEN,
		ODD
	}

	/*** STATE ***/
	
	private Collector<SwapStats> collector;
	int chainCount;
	
	/*** METHODS 
	 * @throws IOException ***/
	
	protected SwapStep() { }
	
	public SwapStep(SwapParity swapParity, long statsEvery) throws IOException
	{
		this.swapParity = swapParity;
		this.statsEvery = statsEvery;
	}

	@Override
	public List<Task> makeTasks(int chainCount)
	{
		this.chainCount = chainCount;
		List<Task> Tasks = new ArrayList<Task>(chainCount);
		for(int i = swapParity == SwapParity.EVEN ? 0 : 1; i + 1 < chainCount; i += 2)
		{
			Tasks.add(new Swapper(i, i+1));
		}
		
		collector = new Collector<SwapStats>(Tasks.size());
		
		return Tasks;
	}
	
	/*** Task CLASS ***/
	
	private class Swapper implements Task
	{
		int[] chainIds;
		
		private long iterationCount;
		private long acceptanceCount;
		
		Swapper(int... chainIds)
		{
			this.chainIds = chainIds;
		}
		
		@Override
		public int[] getChainIds()
		{
			return chainIds;
		}
		
		@Override
		public void step(Chain[] chains) throws MC3KitException
		{
			assert(chains.length == 2);

			Model model0 = chains[0].getModel();
			Model model1 = chains[1].getModel();

			double priorHeatExponent0 = chains[0].getPriorHeatExponent();
			double priorHeatExponent1 = chains[1].getPriorHeatExponent();
			double likelihoodHeatExponent0 = chains[0].getLikelihoodHeatExponent();
			double likelihoodHeatExponent1 = chains[1].getLikelihoodHeatExponent();
			
			double logPriorDiff = model1.getLogPrior() - model0.getLogPrior();
			double logLikelihoodDiff = model1.getLogLikelihood() - model0.getLogLikelihood();
			
			RandomEngine rng = chains[0].getRng();
			
			boolean accepted = false;
			assert(priorHeatExponent1 != Double.POSITIVE_INFINITY);
			assert(likelihoodHeatExponent1 != Double.POSITIVE_INFINITY);
			if(likelihoodHeatExponent0 == Double.POSITIVE_INFINITY
				&& priorHeatExponent0 == Double.POSITIVE_INFINITY)
			{
				if(logPriorDiff + logLikelihoodDiff > 0)
				{
					accepted = true;
				}
			}
			else if(likelihoodHeatExponent0 == Double.POSITIVE_INFINITY)
			{
				if(logLikelihoodDiff > 0)
				{
					accepted = true;
				}
			}
			else if(priorHeatExponent0 == Double.POSITIVE_INFINITY)
			{
				if(logPriorDiff > 0)
				{
					accepted = true;
				}
			}
			else
			{
				double logSwapProb = (priorHeatExponent0 - priorHeatExponent1) * logPriorDiff
					+ (likelihoodHeatExponent0 - likelihoodHeatExponent1) * logLikelihoodDiff;
				
				if(logSwapProb >= 0.0 || log(rng.nextDouble()) < logSwapProb)
				{
					accepted = true;
				}
			}
			
			if(accepted)
			{
				chains[0].setModel(model1);
				model1.setChain(chains[0]);
				chains[1].setModel(model0);
				model0.setChain(chains[1]);
				acceptanceCount++;
			}
			
			iterationCount++;
			
			if(iterationCount % statsEvery == 0)
			{
				List<SwapStats> statsList = collector.takeValue(iterationCount, chainIds[0] / 2,
					new SwapStats(acceptanceCount));
				
				if(statsList != null)
				{
					for(int i = swapParity == SwapParity.EVEN ? 0 : 1; i + 1 < chainCount; i += 2)
					{
					  int lower = chains[0].chainId;
					  int upper = chains[1].chainId;
					  chains[0].getLogger().info(format(
					    "swap rate (%d,%d) for iteration %d: %f", lower, upper, iterationCount, statsList.get(i/2).acceptanceRate)
					  );
					}
				}
				
				acceptanceCount = 0;
			}
		}
	}
	
  class SwapStats
  {
    double acceptanceRate;
    
    SwapStats(long acceptanceCount)
    {
      acceptanceRate = acceptanceCount / (double)statsEvery;
    }
  }
}
