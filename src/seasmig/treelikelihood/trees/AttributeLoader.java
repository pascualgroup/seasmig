package seasmig.treelikelihood.trees;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Set;

public interface AttributeLoader extends Serializable {
	HashMap<String, Object> getAttributes();
	Set<String> getAttributeNames();
}
