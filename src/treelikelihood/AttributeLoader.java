package treelikelihood;
import java.util.HashMap;
import java.util.List;

public interface AttributeLoader {
	
	HashMap<String, Object> getAttributes();
	List<String> getAttributeNames();

}
