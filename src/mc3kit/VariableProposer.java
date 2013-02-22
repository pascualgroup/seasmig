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
public abstract class VariableProposer implements Serializable
{
	private String name;
	
	private int proposalCount;
	private int acceptanceCount;
	
	protected VariableProposer() { }
	
	protected VariableProposer(String name) {
	  this.name = name;
	}
	
	public String getName() {
	  return name;
	}
	
	public void resetTuningPeriod()
	{
		proposalCount = 0;
		acceptanceCount = 0;
	}
	
	public double getAcceptanceRate()
	{
		return (proposalCount == 0)
			? 0.0
			: acceptanceCount / (double)proposalCount;
	}
	
	protected void recordRejection()
	{
		proposalCount++;
	}
	
	protected void recordAcceptance()
	{
		proposalCount++;
		acceptanceCount++;
	}
	
	public abstract void step(Model model)
	    throws MC3KitException;
	public void tune(double targetRate) throws MC3KitException {}
}
