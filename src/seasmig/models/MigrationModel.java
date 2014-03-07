package seasmig.models;

import static mc3kit.util.Utils.makeHierarchicalMap;

import java.io.Serializable;
import java.util.LinkedHashMap;
import java.util.Map;

import mc3kit.Chain;
import mc3kit.Model;
import mc3kit.Variable;
import seasmig.models.TreesLikelihoodVariable.TreesLikelihoodVariableOutputObject;

@SuppressWarnings("serial")
public abstract class MigrationModel extends Model implements Serializable {
	
	protected MigrationModel() { };
	
	public MigrationModel(Chain initialChain) {
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
		
		if (outputObject.seqLikelihood!=null) {
			for (int i=0;i<outputObject.seqLikelihood.length;i++) {
				flatMap.put("seqLikelihood."+i, outputObject.seqLikelihood[i]);
			}
		}
		
		if (outputObject.locLikelihood!=null) {
			for (int i=0;i<outputObject.locLikelihood.length;i++) {
				flatMap.put("locLikelihood."+i, outputObject.locLikelihood[i]);
			}
		}
		
		if (outputObject.asrTrees!=null) {
			for (int i=0;i<outputObject.asrTrees.length;i++) {
				flatMap.put("trees."+i,String.format("%s",outputObject.asrTrees[i]));
			}
		}
		if (outputObject.smTrees!=null) {
			for (int i=0;i<outputObject.smTrees.length;i++) {
				flatMap.put("smTrees."+i,String.format("%s",outputObject.smTrees[i]));
			}
		}
		if (outputObject.smTransitions!=null) {
			for (int i=0;i<outputObject.smTransitions.length;i++) {
				flatMap.put("smTransitions."+i,String.format("%s",outputObject.smTransitions[i]));
			}
		}
		if (outputObject.smTipDwellings!=null) {
			for (int i=0;i<outputObject.smTipDwellings.length;i++) {
				flatMap.put("smTipDwellings."+i,String.format("%s",outputObject.smTipDwellings[i]));
			}
		}
		if (outputObject.smLineages!=null) {
			for (int i=0;i<outputObject.smLineages.length;i++) {
				flatMap.put("smLineages."+i,String.format("%s",outputObject.smLineages[i]));
			}
		}
		if (outputObject.smDescendants!=null) {
			for (int i=0;i<outputObject.smDescendants.length;i++) {
				flatMap.put("smDescendants."+i,String.format("%s",outputObject.smDescendants[i]));
			}
		}
		if (outputObject.smTrunkStats!=null) {
			for (int i=0;i<outputObject.smTrunkStats.length;i++) {
				flatMap.put("smTrunkStats."+i,String.format("%s",outputObject.smTrunkStats[i]));
			}
		}	

		///////////////////////////////////////
		if (outputObject.pies!=null) {
			for (int i=0;i<outputObject.pies.length;i++) {
				flatMap.put("pi."+i,String.format("%s",outputObject.pies[i]));
			}
		}
		
		///////////////////////////////////////
		if (outputObject.seqMutationStats!=null) {
			for (int i=0;i<outputObject.seqMutationStats.length;i++) {
				flatMap.put("seqMutationStats."+i,String.format("%s",outputObject.seqMutationStats[i]));
			}
		}
		
		// Alternative Tree Output
		if (outputObject.altBranches!=null && outputObject.altNodes!=null) {
			for (int i=0;i<outputObject.altBranches.length;i++) {				
				flatMap.put("altTreeOutput."+i+"."+"nodes",outputObject.altNodes[i]);
			}
		}
				
		return makeHierarchicalMap(flatMap);
	}

}
