package seasmig.models;

import static mc3kit.util.Utils.makeHierarchicalMap;

import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

import mc3kit.Chain;
import mc3kit.Model;
import mc3kit.Variable;

@SuppressWarnings("serial")
public abstract class SeasonalMigrationModel extends Model implements Serializable {
	
	protected SeasonalMigrationModel() { };
	
	public SeasonalMigrationModel(Chain initialChain) {
		super(initialChain);
	}

	@Override
	public Map<String, Object> makeHierarchicalSample() {		
		Map<String, Object> flatMap = new LinkedHashMap<String, Object>();

		flatMap.put("iterationCount", getChain().getIterationCount() + 1);
		flatMap.put("logPrior", getLogPrior());
		flatMap.put("logLikelihood", getLogLikelihood());

		for(Variable var : getUnobservedVariables()) {
			flatMap.put(var.getName(), var.makeOutputObject());
		}
		
		String[] trees = (String[]) getVariable("likeVar").makeOutputObject();
		for (int i=0;i<trees.length;i++) {
			flatMap.put("trees."+i,String.format("%s",trees[i]));
		}
		
		return makeHierarchicalMap(flatMap);
	}

}
