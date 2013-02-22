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

import com.google.gson.*;
import java.util.*;
import java.util.Map.Entry;
import java.io.*;

public class JsonUtils
{
	public static <T> T parseObject(Class<T> cls, String filename) throws FileNotFoundException
	{
		return parseObject(cls, new FileReader(filename));
	}
	
	public static  <T> T parseObject(Class<T> cls, Reader reader)
	{
		Gson gson = new Gson();
		JsonParser parser = new JsonParser();
		JsonObject raw = parser.parse(reader).getAsJsonObject();
		return gson.fromJson(raw, cls);
	}
	
	public static <T> List<T> parseList(Class<T> cls, String filename) throws FileNotFoundException
	{
		return parseList(cls, new FileReader(filename));
	}
	
	public static <T> List<T> parseList(Class<T> cls, Reader reader)
	{
		Gson gson = new Gson();
		JsonParser parser = new JsonParser();
		
		JsonObject raw = parser.parse(reader).getAsJsonObject();
		
		List<T> list = new ArrayList<T>();
		for(Entry<String, JsonElement> entry : raw.entrySet())
		{
			T obj = gson.fromJson(entry.getValue(), cls);
			list.add(obj);
		}
		
		return list;
	}
	
	public static <T> Map<String, T> parseMap(Class<T> cls, String filename) throws FileNotFoundException
	{
		return parseMap(cls, new FileReader(filename));
	}
	
	public static <T> Map<String, T> parseMap(Class<T> cls, Reader reader)
	{
		Gson gson = new Gson();
		JsonParser parser = new JsonParser();
		
		JsonObject raw = parser.parse(reader).getAsJsonObject();
		
		Map<String, T> map = new LinkedHashMap<String, T>();
		for(Entry<String, JsonElement> entry : raw.entrySet())
		{
			T obj = gson.fromJson(entry.getValue(), cls);
			map.put(entry.getKey(), obj);
		}
		
		return map;
	}
}
