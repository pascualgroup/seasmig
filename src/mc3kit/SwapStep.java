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

import static mc3kit.util.Utils.*;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import static java.lang.String.*;

import cern.jet.random.engine.RandomEngine;
import mc3kit.util.Collector;

@SuppressWarnings("serial")
public class SwapStep implements Step
{
  ChainParity swapParity;
  long statsEvery = 100;


	/*** STATE ***/
	
	private Collector<SwapStats> collector;
	int chainCount;
	
	/*** METHODS 
	 * @throws IOException ***/
	
	protected SwapStep() { }
	
	public SwapStep(ChainParity swapParity, long statsEvery) throws IOException
	{
		this.swapParity = swapParity;
		this.statsEvery = statsEvery;
	}
	
	@Override
	public List<Task> makeTasks(int chainCount)
	{
		this.chainCount = chainCount;
		List<Task> tasks = new ArrayList<Task>(chainCount);
		for(int i = swapParity == ChainParity.EVEN ? 0 : 1; i + 1 < chainCount; i += 2)
		{
			tasks.add(new Swapper(i, i+1));
		}
		
		collector = new Collector<SwapStats>(tasks.size());
		
		return tasks;
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
					new SwapStats(chainIds, acceptanceCount));
				
				Logger logger = chains[0].getMCMC().getLogger("mc3kit.SwapStep");
				if(statsList != null && logger.isLoggable(Level.INFO))
				{
				  Map<String, Object> swapStats = new LinkedHashMap<String, Object>();
					for(int i = swapParity == ChainParity.EVEN ? 0 : 1; i + 1 < chainCount; i += 2)
					{
					  SwapStats stats = statsList.get(i/2);
	          swapStats.put(
	              format("(%d,%d)", stats.chainIds[0], stats.chainIds[1]),
	              stats.acceptanceRate
	          );
					}
					
					logger.log(Level.INFO,
					  "SwapStep acceptance rates",
					  makeMap("swapParity", swapParity.toString(), "swapStats", swapStats)
					);
				}
				
				acceptanceCount = 0;
			}
		}
	}
	
  class SwapStats
  {
    int[] chainIds;
    double acceptanceRate;
    
    SwapStats(int[] chainIds, long acceptanceCount)
    {
      this.chainIds = chainIds;
      acceptanceRate = acceptanceCount / (double)statsEvery;
    }
  }
}
