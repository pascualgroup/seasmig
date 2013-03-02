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

package mc3kit.output;

import java.io.FileNotFoundException;
import static java.lang.String.*;
import java.util.*;

import mc3kit.*;
import mc3kit.output.SampleWriterFactory;

@SuppressWarnings("serial")
public class SampleOutputStep implements Step
{
  String filename;
  String format;
  boolean useQuotes;
  long thin;
  int chainId;
	
  public SampleOutputStep(String filename, long thin) {
    this(filename, null, false, thin, 0);
  }
  
	public SampleOutputStep(String filename, String format, boolean useQuotes, long thin, int chainId) {
	  this.filename = filename;
	  this.format = format;
	  this.useQuotes = useQuotes;
	  this.thin = thin;
	  this.chainId = chainId;
	}

	/*** METHODS ***/

	@Override
	public List<Task> makeTasks(int chainCount) throws MC3KitException
	{
	  List<Task> Tasks = new ArrayList<Task>();
		Tasks.add(new SampleOutputTask());
		return Tasks;
	}
	
	/*** Task CLASS ***/
	
	private class SampleOutputTask implements Task
	{
		SampleWriter writer;
		
		private long iterationCount;
		
		@Override
		public int[] getChainIds()
		{
			return new int[] { chainId };
		}
		
		SampleOutputTask() throws MC3KitException
		{
			try {
        writer = SampleWriterFactory.getFactory().createSampleWriter(filename, format, useQuotes);
      }
      catch(FileNotFoundException e) {
        throw new MC3KitException("File not found", e);
      }
		}
		
		@Override
		public void step(Chain[] chains) throws MC3KitException
		{
      iterationCount++;
      
      Chain chain = chains[0];
			if(iterationCount % thin == 0)
			{
	      chain.getLogger().fine(format("Writing sample %d", iterationCount));
				Model model = chain.getModel();
				
				writer.writeSample(model);
			}
		}
	}
}
