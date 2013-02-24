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

import java.io.*;
import java.util.*;

import mc3kit.*;

@SuppressWarnings("serial")
public class CsvSampleWriter implements SampleWriter
{
  private String filename;
	private transient PrintWriter writer;
	private String delimiter;
	private boolean useQuotes;
	
	private Set<String> keys;
	
	protected CsvSampleWriter() { }
	
	public CsvSampleWriter(String filename, String delimiter, boolean useQuotes) throws FileNotFoundException
	{
	  this.filename = filename;
		this.writer = new PrintWriter(new FileOutputStream(filename, false));
		this.delimiter = delimiter;
		this.useQuotes = useQuotes;
	}
	
	private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
	  in.defaultReadObject();
	  this.writer = new PrintWriter(new FileOutputStream(filename, true));
	}
	
	@Override
	public synchronized void writeSample(Model model) throws MC3KitException
	{
	  Map<String, String> samp = model.makeFlatSample();
	  
		// Establish keys to use and write headers if unestablished
		if(keys == null)
		{
			keys = new LinkedHashSet<String>(samp.keySet());
			
			int i = 0;
			for(String key : keys)
			{
				writer.write(formattedString(key));
				if(i < keys.size() - 1)
				{
					writer.write(delimiter);
				}
				i++;
			}
			writer.println();
		}
		
		// Write quoted, quote-escaped data as string
		int i = 0;
		for(String key : keys)
		{
			String value = samp.get(key);
			if(value != null)
			{
				String str = formattedString(value);
				writer.write(str);
			}
			if(i < keys.size() - 1)
			{
				writer.write(delimiter);
			}
			i++;
		}
		writer.println();
	}
	
	private String formattedString(String str) throws MC3KitException
	{
		if(useQuotes)
			return String.format("\"%s\"", str.replaceAll("\"", "\"\""));
		else
		{
			if(str.contains(delimiter))
				throw new MC3KitException(
					String.format("Delimiter found in string in non-quote mode: %s.", str)
				);
			return str;
		}
	}
}
