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

import java.util.*;
import static org.junit.Assert.*;
import static mc3kit.util.Utils.*;

import mc3kit.MC3KitException;

import org.junit.*;

public class UtilsTest {

  @Before
  public void setUp() throws Exception {
  }

  @After
  public void tearDown() throws Exception {
  }

  @Test
  public void flatMap() throws MC3KitException {
    Map<String, Object> flatMap = makeMap(
      "key1", 2, "key2", new double[] {5, 4, 3}, "key3", new LinkedHashMap<String, Object>()
    );
    assertEquals(2, flatMap.get("key1"));
    assertArrayEquals(new double[] {5, 4, 3}, (double[])flatMap.get("key2"), 0.0);
    assertEquals(new LinkedHashMap<String, Object>(), flatMap.get("key3"));
  }
  
  @SuppressWarnings("unchecked")
  @Test
  public void oneLevelHMap() throws Exception {
    Map<String, Object> flatMap = makeMap(
      "0.mean", 4.5,
      "0.prec", 4.5,
      "0.name", "zero"
    );
    Map<String, Object> hMap = makeHierarchicalMap(flatMap);
    Map<String, Object> obj = (Map<String, Object>)hMap.get("0");
    System.err.printf("%s\n", hMap);
    assertEquals(4.5, obj.get("mean"));
    assertEquals(4.5, obj.get("prec"));
    assertEquals("zero", obj.get("name"));
  }
}
