/***
  DZ
 ***/

package mc3kit.distributions;

import mc3kit.*;
import mc3kit.proposal.MHUniformIntRandomProposer;
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class UniformIntDistribution extends IntDistribution {

	// Random uniform distribution from min to max inclusive .....

	static int MaxRange = 1000000; // TODO: find precision error free working range
	static int MinRange = -1000000; // TODO: find precision error working range
	int min; 
	int max; 
	double logP; 
	ModelEdge minEdge;
	ModelEdge maxEdge;

	protected UniformIntDistribution() { }

	public UniformIntDistribution(Model model) {
		this(model, null, MinRange, MaxRange); 
	}

	public UniformIntDistribution(Model model, String name) {
		this(model, name, MinRange,MaxRange);
	}

	public UniformIntDistribution(Model model,int min, int max){
		this(model, null, min, max);
	}

	public UniformIntDistribution(Model model, String name,int min, int max){
		super(model, name);

		if(min >= max) {
			throw new IllegalArgumentException("min: "+min+" must be less than max: "+max);
		}

		if(min >= MaxRange) {
			throw new IllegalArgumentException("min: "+min+" must be less than MaxRange: "+MaxRange);
		}

		if(max <= MinRange) {
			throw new IllegalArgumentException("max: "+max+" must be greater than MinRange: "+MinRange);
		}

		this.min = min;
		this.max = max;
		logP=-log(max-min+1); // TODO: check this
	}

	public <T extends ModelNode & IntValued> UniformIntDistribution setP(T node) throws MC3KitException {
		minEdge = updateEdge(minEdge, node);
		maxEdge = updateEdge(maxEdge, node);
		return this;
	}

	@Override
	public VariableProposer makeVariableProposer(String varName) {
		return new MHUniformIntRandomProposer(varName,min,max);
	}

	@Override
	public double getLogP(Variable var) {  
		return logP; 
	}

	@Override
	public void sample(Variable var) {
		// TODO: check this...
		double min = minEdge == null ? this.min : getDoubleValue(minEdge);
		double max = maxEdge == null ? this.max : getDoubleValue(maxEdge);
		double doubleValue = this.getRng().nextDouble();
		((IntVariable)var).setValue( (int)min + (int)(doubleValue * ((max - min) + 1))); 
	}
}
