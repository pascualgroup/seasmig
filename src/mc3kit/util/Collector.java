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

package mc3kit.util;

import java.io.Serializable;
import java.util.*;

/**
 * Class that accumulates values that may arrive out of order
 * from multiple threads.
 * 
 * @param <V> The value type. 
 */
@SuppressWarnings("serial")
public class Collector<V> implements Serializable
{
	protected int count;
	protected int counter;
	protected Map<Long, IterationData> iterationMap;
	
	protected Collector() { }
	
	public Collector(int count)
	{
		this.count = count;
		iterationMap = new HashMap<Long, IterationData>();
	}
	
	public synchronized List<V> takeValue(long iteration, int index, V value)
	{
		counter++;
		
		IterationData iterData = iterationMap.get(iteration);
		if(iterData == null)
		{
			iterData = new IterationData();
			iterationMap.put(iteration, iterData);
		}
		iterData.set(index, value);
		
		List<V> returnValues = null;
		if(iterData.isDone())
		{
			returnValues = iterData.getValues();
			iterationMap.remove(iteration);
		}
		return returnValues;
	}
	
	private class IterationData
	{
		int counter;
		ArrayList<V> values;
		
		IterationData()
		{
			counter = 0;
			values = new ArrayList<V>(count);
			for(int i = 0; i < count; i++)
				values.add(null);
		}
		
		void set(int index, V value)
		{
			assert(values.get(index) == null);
			values.set(index,  value);
			counter++;
		}
		
		ArrayList<V> getValues()
		{
			return values;
		}
		
		boolean isDone()
		{
			return counter == count;
		}
	}
}
