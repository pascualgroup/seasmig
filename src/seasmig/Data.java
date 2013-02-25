package seasmig;
import java.io.Serializable;
import java.util.List;

public interface Data extends Serializable{
	List<LikelihoodTree> getTrees();

	int getNumLocations();
}
