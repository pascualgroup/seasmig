package seasmig.models;

import static mc3kit.util.Utils.makeHierarchicalMap;

import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

import seasmig.models.TreesLikelihoodVariable.TreesLikelihoodVariableOutputObject;
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

		TreesLikelihoodVariableOutputObject outputObject = (TreesLikelihoodVariable.TreesLikelihoodVariableOutputObject) getVariable("likeVar").makeOutputObject();
		if (outputObject.trees!=null) {
			for (int i=0;i<outputObject.trees.length;i++) {
				flatMap.put("trees."+i,String.format("%s",outputObject.trees[i]));
			}
		}
		if (outputObject.smTransitions!=null) {
			for (int i=0;i<outputObject.trees.length;i++) {
				flatMap.put("smTransitions."+i,String.format("%s",outputObject.smTransitions[i]));
			}
		}
		if (outputObject.smTipDwelings!=null) {
			for (int i=0;i<outputObject.trees.length;i++) {
				flatMap.put("smDwelings."+i,String.format("%s",outputObject.smTipDwelings[i]));
			}
		}
		if (outputObject.smLineages!=null) {
			for (int i=0;i<outputObject.trees.length;i++) {
				flatMap.put("smLineages."+i,String.format("%s",outputObject.smLineages[i]));
			}
		}
		
		return makeHierarchicalMap(flatMap);
	}

}
