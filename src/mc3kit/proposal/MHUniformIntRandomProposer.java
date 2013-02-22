package mc3kit.proposal;

import static mc3kit.util.Math.*;
import static java.lang.String.format;
import cern.jet.random.engine.RandomEngine;
import mc3kit.*;

@SuppressWarnings("serial")
public class MHUniformIntRandomProposer extends VariableProposer {

	double min;
	double max;

	protected MHUniformIntRandomProposer() { }

	public MHUniformIntRandomProposer(String name, double min, double max) {
		super(name);
		this.min = min;
		this.max = max;
	} 

	@Override
	public void step(Model model) throws MC3KitException {
		Chain chain = model.getChain();
		RandomEngine rng = chain.getRng();

		chain.getLogger().finest("MHUniformProposer stepping");

		double oldLogLikelihood = model.getLogLikelihood();
		double oldLogPrior = model.getLogPrior();

		chain.getLogger().finest(format("oldLP, oldLL: %f, %f", oldLogPrior, oldLogLikelihood));

		IntVariable rv = model.getIntVariable(getName());

		int oldValue = rv.getValue();

		// TODO: check this...
		int newValue =  oldValue+1 + (int)(rng.nextDouble() * ((max-1 - min) + 1));
		if (newValue>max) {
			newValue=newValue-(int)max+(int)min-1;
		}

		chain.getLogger().finest(format("oldVal, newVal = %d, %d", oldValue, newValue));

		model.beginProposal();
		rv.setValue(newValue);
		model.endProposal();

		double newLogPrior = model.getLogPrior();
		double newLogLikelihood = model.getLogLikelihood();

		chain.getLogger().finest(format("newLP, newLL: %f, %f", newLogPrior, newLogLikelihood));

		boolean accepted = shouldAcceptMetropolisHastings(rng,
				chain.getPriorHeatExponent(), chain.getLikelihoodHeatExponent(),
				oldLogPrior, oldLogLikelihood,
				newLogPrior, newLogLikelihood,
				0.0
				);

		if(accepted)
		{
			model.acceptProposal();
			recordAcceptance();
		}
		else
		{
			model.beginRejection();
			rv.setValue(oldValue);
			model.endRejection();

			recordRejection();
		}
	}

	// TODO: tune...
}

