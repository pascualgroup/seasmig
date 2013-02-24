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
import java.util.Map;

import mc3kit.*;

import com.google.gson.*;

@SuppressWarnings("serial")
public class JsonsSampleWriter implements SampleWriter
{
  private String filename;
	private transient PrintWriter writer;
	private transient Gson gson;
	
	protected JsonsSampleWriter() { }
	
	public JsonsSampleWriter(String filename) throws FileNotFoundException
	{
	  this.filename = filename;
	  initializeTransientObjects(false);
	}
  
	private void initializeTransientObjects(boolean append) throws FileNotFoundException {
    this.gson = new GsonBuilder()
    .serializeSpecialFloatingPointValues()
    .setPrettyPrinting()
    .create();
    this.writer = new PrintWriter(new FileOutputStream(filename, append));
	}
	
  private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException {
    in.defaultReadObject();
    initializeTransientObjects(true);
  }
	
	@Override
	public synchronized void writeSample(Model model) throws MC3KitException
	{
	  writer.println("---");
		gson.toJson(model.makeHierarchicalSample(), writer);
		writer.println();
		writer.flush();
	}

	@Override
	public void writeFlatData(Map<String, String> flatData)
			throws MC3KitException {
		// TODO: check this!!!!
		writer.println("---");
		gson.toJson(flatData, writer);
		writer.println();
		writer.flush();
	}
}
