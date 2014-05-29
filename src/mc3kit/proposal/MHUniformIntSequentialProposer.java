package mc3kit.proposal;

import static mc3kit.util.Math.*;
import static java.lang.String.format;
import cern.jet.random.engine.RandomEngine;
import mc3kit.*;



@SuppressWarnings("serial")
public class MHUniformIntSequentialProposer extends VariableProposer {

	int min;
	int max;

	protected MHUniformIntSequentialProposer() { }

	public MHUniformIntSequentialProposer(String name, int min, int max) {
		super(name);
		this.min = min;
		this.max = max;
	} 

	@Override
	public void step(Model model) throws MC3KitException {
		Chain chain = model.getChain();
		
		chain.getLogger().finest("MHUniformIntProposer stepping");

		double oldLogLikelihood = model.getLogLikelihood();
		double oldLogPrior = model.getLogPrior();

		chain.getLogger().finest(format("oldLP, oldLL: %f, %f", oldLogPrior, oldLogLikelihood));

		IntVariable rv = (IntVariable) model.getVariable(getName());

		int oldValue = rv.getValue();

		RandomEngine rng = model.getRng();
		
		int newValue = ((oldValue+1)>max?1:(oldValue+1));;

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

