package mc3kit.util;

import java.util.*;
import static mc3kit.util.Utils.*;
import java.util.logging.*;
import com.google.gson.*;

public class JsonsLogFormatter extends java.util.logging.Formatter {
  private Gson gson;
  
  public JsonsLogFormatter()
  {
    gson = new GsonBuilder().setPrettyPrinting().create();
  }
  
  @Override
  public String format(LogRecord r) {
    Map<String, Object> map = makeMap(
      "date", new Date(r.getMillis()).toString(),
      "logger", r.getLoggerName(),
      "class", r.getSourceClassName(),
      "level", r.getLevel().toString(),
      "message", r.getMessage()
    );
    
    // Passed-in objects, assumed JSON-compatible
    Object[] params = r.getParameters();
    if(params != null) {
      if(params.length == 1) {
        map.put("info", params[0]);
      }
      else {
        map.put("info", params);
      }
    }
    
    // Stack trace as array of strings
    if(r.getThrown() != null) {
      StackTraceElement[] stackTrace = r.getThrown().getStackTrace();
      String[] stArray = new String[stackTrace.length];
      for(int i = 0; i < stackTrace.length; i++) {
        stArray[i] = stackTrace[i].toString();
      }
      map.put("thrown", stArray);
    }
    
    return String.format("---\n%s\n", gson.toJson(map));
  }
  
}
