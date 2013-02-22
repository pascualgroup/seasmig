package treelikelihood;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public interface AttributeLoader {
	HashMap<String, Object> getAttributes();
	Set<String> getAttributeNames();
}
