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

@SuppressWarnings("serial")
public class MultiCounter<T extends Enum<?>> implements Serializable
{
	long count;
	Map<T, Long> counts;
	
	public MultiCounter()
	{
		count = 0;
		counts = new HashMap<T, Long>();
	}
	
	public long getCount()
	{
		return count;
	}
	
	public long getCount(T type)
	{
		return counts.containsKey(type) ? counts.get(type) : 0;
	}
	
	public double getRate(T type)
	{
		return count == 0 ? 0.0 : getCount(type) / (double)count;
	}
	
	public void reset()
	{
		count = 0;
		counts.clear();
	}
	
	public void record(T... types)
	{
		count++;
		for(T type : types)
		{
			increment(type);
		}
	}
	
	private void increment(T type)
	{
		counts.put(type, getCount(type) + 1);
	}
}
