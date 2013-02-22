/***
  This file is part of mc3kit.

  Copyright (C) 2013 Edward B. Baskerville

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***/

package mc3kit.distributions;

import mc3kit.*;
import mc3kit.proposal.GibbsIntProposer;
import static java.lang.Math.*;

@SuppressWarnings("serial")
public class UniformIntDistribution extends IntDistribution {

	// Random uniform distribution from min to max inclusive .....

	static int MaxRange = 1000000; // TODO: find precision error free working range
	static int MinRange = -1000000; // TODO: find precision error working range
	double min; // TODO: int?
	double max; // TODO: int?
	double logP; 
	ModelEdge minEdge;
	ModelEdge maxEdge;

	protected UniformIntDistribution() { }

	public UniformIntDistribution(Model model) {
		this(model, null, MinRange, MaxRange); // TODO: will work for large ints correctly ????
	}

	public UniformIntDistribution(Model model, String name) {
		this(model, name, MinRange,MaxRange); // TODO: will work for large ints correctly ????
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

	public <T extends ModelNode & IntValued> UniformIntDistribution setP(T node) throws ModelException {
		minEdge = updateEdge(minEdge, node);
		maxEdge = updateEdge(maxEdge, node);
		return this;
	}

	@Override
	public double getLogP(Variable var) {  
		return logP; 
	}


	@Override
	public VariableProposer makeVariableProposer(String varName) {
		return new GibbsIntProposer(varName,min,max);
	}

	@Override
	public void sample(Variable var) {
		double min = minEdge == null ? this.min : getDoubleValue(minEdge);
		double max = maxEdge == null ? this.max : getDoubleValue(maxEdge);
		((IntVariable)var).setValue( (int)min + (int)(getRng().nextDouble() * ((max - min) + 1))); // TODO: check this
	}
}
