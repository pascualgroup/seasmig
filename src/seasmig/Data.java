package seasmig;
import java.util.List;

import treelikelihood.*;

public interface Data {
	List<LikelihoodTree> getTrees();

	int getNumLocations();
}
