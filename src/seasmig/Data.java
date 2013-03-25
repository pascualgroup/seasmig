package seasmig;
import java.io.Serializable;
import java.util.List;

import seasmig.treelikelihood.LikelihoodTree;

public interface Data extends Serializable{
	List<LikelihoodTree> getTrees();

	int getNumLocations();
}
