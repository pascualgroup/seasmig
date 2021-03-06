package mc3kit.distributions;

import static org.junit.Assert.*;

import mc3kit.*;
import mc3kit.proposal.UnivariateProposalStep;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

@SuppressWarnings("serial")
public class UniformIntIntegratorTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}


	@Test
	public void test() throws Throwable {

		long burnIn = 5000;
		long iterCount = 100000;

		MCMC mcmc = new MCMC();
		mcmc.setRandomSeed(101L);

		ModelFactory mf = new ModelFactory() {
			@Override
			public Model createModel(Chain initialChain) throws MC3KitException {
				Model m = new Model(initialChain);

				m.beginConstruction();
				new IntVariable(m, "v", new UniformIntIntegrator(m, 1,4));
				m.endConstruction();

				return m;
			}
		};

		mcmc.setModelFactory(mf);

		UnivariateProposalStep proposalStep = new UnivariateProposalStep(0.25, 100, burnIn);
		mcmc.addStep(proposalStep);

		// Run, collect statistics, and check moments against expected distribution
		double sum = 0;
		for(long i = 0; i < iterCount; i++) {
			mcmc.step();
			mcmc.getModel().recalculate(1E-4);

			assertEquals(i + 1, mcmc.getIterationCount());
			if(i >= burnIn) {				
				int val = ((IntVariable)mcmc.getModel().getVariable("v")).getValue();
				System.err.println(val);
				sum += val;								 
			}
		}

		double N = iterCount - burnIn;

		double mean = sum / (double)N;
		System.err.printf("mean = %f\n", mean);
		assertEquals(2.5, mean, 0.02);
	}
}


