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

@SuppressWarnings("serial")
public class MC3KitException extends Exception {
  public MC3KitException() {
    super();
  }
  
  public MC3KitException(Throwable cause) {
    super(cause);
  }
  
  public MC3KitException(String msg) {
    super(msg);
  }

  public MC3KitException(String msg, Throwable cause) {
    super(msg, cause);
  }
}
