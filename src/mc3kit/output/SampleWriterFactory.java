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

import mc3kit.MC3KitException;

public class SampleWriterFactory {
  private static SampleWriterFactory factory;

  public static synchronized SampleWriterFactory getFactory() {
    if(factory == null) {
      factory = new SampleWriterFactory();
    }
    return factory;
  }

  private SampleWriterFactory() {
  }

  public SampleWriter createSampleWriter(String filename) throws MC3KitException, FileNotFoundException {
    return createSampleWriter(filename, null, false);
  }

  public SampleWriter createSampleWriter(String filename, String format)
      throws MC3KitException, FileNotFoundException {
    return createSampleWriter(filename, format, false);
  }

  public SampleWriter createSampleWriter(String filename, String format,
      boolean useQuotes) throws FileNotFoundException, MC3KitException {
    if(format == null) {
      String[] filenamePieces = filename.split("\\.");
      if(filenamePieces.length > 1)
        format = filenamePieces[filenamePieces.length - 1];
    }

    if(format.equalsIgnoreCase("jsons")) {
      return new JsonsSampleWriter(filename);
    }
    else if(format.equalsIgnoreCase("csv")) {
      return new CsvSampleWriter(filename, ",", useQuotes);
    }
    else if(format.equalsIgnoreCase("txt")) {
      return new CsvSampleWriter(filename, "\t", useQuotes);
    }
    else {
      throw new MC3KitException(String.format("Unknown format %s.", format));
    }
  }
}
