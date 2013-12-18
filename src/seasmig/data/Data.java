package seasmig.data;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import seasmig.treelikelihood.LikelihoodTree;

public interface Data extends Serializable{
	List<ArrayList<LikelihoodTree>> getTrees();

	int getNumLocations();
}
