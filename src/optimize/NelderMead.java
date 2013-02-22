package optimize;

import java.util.Collections;
import java.util.Comparator;
import java.util.Vector;

import jebl.math.Random;

import org.javatuples.Pair;

public class NelderMead {

	double alpha = 1;
	double gamma = 2;
	double rho = -0.5;
	double sigma = 0.5;
	double nu = 0.00001;
	private double terminationCondition = 1e-4;
	private Vector<Vector<Double>> x0s;
	
	
	public NelderMead() {
		this(makeInitPoints(15, 5));
	}
	
	public NelderMead(Vector<Vector<Double>> x0s) {
		this(x0s,1e-4);
	}
	
	public NelderMead(Vector<Vector<Double>> x0s, double terminationCondition) {
		this.x0s=x0s;
		this.terminationCondition=terminationCondition;
	}
	

	// Keeping pair of <x, f(x)>		
	public class XsFxsComparator implements Comparator<Pair<Vector<Double>,Double>> {
		@Override
		public int compare(Pair<Vector<Double>,Double> o1,Pair<Vector<Double>,Double> o2) {
			return o1.getValue1().compareTo(o2.getValue1());
		}
	}

	void addMultSecond(Vector<Double> first, Vector<Double> second, Double alpha) {						
		for (int i=0;i<first.size();i++) {
			first.set(i, first.get(i)+second.get(i)*alpha);
		}				
	}

	void mult(Vector<Double> x, Double alpha) {						
		for (int i=0;i<x.size();i++) {
			x.set(i, x.get(i)*alpha);
		}				
	}

	void order(Vector<Pair<Vector<Double>,Double>> simplex) {			
		Collections.sort(simplex, new XsFxsComparator());
	}

	Vector<Double> centeroid_of_all_expect_last(Vector<Pair<Vector<Double>,Double>> simplex) {
		Vector<Double> xo =zeros(simplex.get(0).getValue0().size());
 
		for (int i=0;i<simplex.size()-1;i++) {
			addMultSecond(xo, simplex.get(i).getValue0(), 1.0);			
		}

		mult(xo, 1.0/(simplex.size()-1));			

		return xo;
	}

	private Vector<Double> zeros(int size) {
		Vector<Double> returnValue = new Vector<Double>();
		for (int i=0;i<size;i++) {
			returnValue.add(0.0);
		}
		
		return returnValue;
	}

	Vector<Double> adjust(Vector<Double> xo, Vector<Double> xn, Double delta) {
		Vector<Double> xr = zeros(xo.size());

		for (int i=0;i<xo.size();i++) {
			xr.set(i, xo.get(i)+delta*(xo.get(i)-xn.get(i)));
		}		

		return xr;
	}

	void reduce(Vector<Pair<Vector<Double>, Double>> simplex) {	
		for (int i=1; i<(simplex.size()-1);i++) {
			sigma = -Math.log(Random.nextDouble())*(0.5);
			Vector<Double> xi_new = adjust(simplex.get(i).getValue0(),simplex.get(0).getValue0(),sigma);
			simplex.set(i, new Pair<Vector<Double>,Double>(xi_new,feval(xi_new))); 
		}			
	}
	
	void noise(Vector<Pair<Vector<Double>, Double>> simplex) {
		for (int i=0; i<(simplex.size()-1);i++) {
			nu = -Math.log(Random.nextDouble())*0.0001;
			Vector<Double> xi_new = simplex.get(i).getValue0();
			for (int j=0; j<simplex.get(i).getValue0().size();j++) {
				xi_new.set(j, xi_new.get(j)+Random.nextDouble()*nu);
			}
			
			simplex.set(i, new Pair<Vector<Double>,Double>(xi_new,feval(xi_new))); 
		}			
	}


	static Vector<Vector<Double>> makeInitPoints(int n, int m) {			
		Vector<Vector<Double>> x0s = new Vector<Vector<Double>>(n);
		for (int i=0;i<n;i++) {
			x0s.add(new Vector<Double>(m));
			for (int j=0;j<m;j++) {
				x0s.get(i).add(Random.nextDouble()*100);
			}
		}
		
		return x0s;
	}

	
	public Pair<Vector<Double>,Double> optimize() {
		int n = x0s.size();

		// Evaluate						
		Vector<Pair<Vector<Double>,Double>> simplex = evaluate(x0s);
		
		while (!termination_condition(simplex)) {
			alpha = -Math.log(Random.nextDouble())*1.0;
			gamma = -Math.log(Random.nextDouble())*2.0;
			rho = -Math.log(Random.nextDouble())*(-0.5);	
			
			// Order
			order(simplex);
			Double fxn=simplex.get(simplex.size()-1).getValue1();
			Double fx0=simplex.get(0).getValue1();
			Vector<Double> xn=simplex.get(simplex.size()-1).getValue0();
			// Center of gravitiy
			Vector<Double> xo = centeroid_of_all_expect_last(simplex);
			// Reflect
			Vector<Double> xr = adjust(xo,xn,alpha);
			Double fxr = feval(xr);
			if (fxr>=fx0 && fxr<simplex.get(n-2).getValue1()) {
				simplex.set(n-1, new Pair<Vector<Double>,Double>(xr,fxr));
				continue;
			}
			// Expand
			// if the reflected point is the best point so far f(xr)<f(x0)
			if (fxr<simplex.get(0).getValue1()) {
				Vector<Double> xe = adjust(xo,xn,gamma);
				Double fxe = feval(xe);
				if (fxe<fxr) {
					simplex.set(n-1,new Pair<Vector<Double>,Double>(xe,fxe));
					continue;
				}
				simplex.set(n-1, new Pair<Vector<Double>,Double>(xr,fxr));
				continue;
			}
			// Contract
			assert fxr>=simplex.get(n-1).getValue1();
			Vector<Double> xc = adjust(xo,xn,rho);
			Double fxc = feval(xc);
			if (fxc<fxn) {
				simplex.set(n-1, new Pair<Vector<Double>,Double>(xc,fxc));
				continue;
			}				
			// Reduce
			reduce(simplex);
			// Noise
			noise(simplex);
		}
		
		order(simplex);
		return (simplex.get(0));

	}

	private boolean termination_condition(
			Vector<Pair<Vector<Double>, Double>> simplex) {			
		return (Math.abs(simplex.get(0).getValue1())<terminationCondition );
	}

	private Vector<Pair<Vector<Double>, Double>> evaluate(Vector<Vector<Double>> x0s) {
		Vector<Pair<Vector<Double>, Double>> simplex = new Vector<Pair<Vector<Double>, Double>>(x0s.size());
		for (int i=0; i<x0s.size();i++) {
			simplex.add(new Pair<Vector<Double>, Double>((Vector<Double>) x0s.get(i).clone(),feval(x0s.get(i))));
		}
		return simplex;
	}

	private Double feval(Vector<Double> vector) {
		Double returnValue=0.0;
		for (int i=0;i<vector.size();i++) {
			returnValue+=Math.abs((1-Math.pow(vector.get(i), i+1)))+Math.abs((vector.get(i)-1)*Random.nextDouble());
		}
		return returnValue;
	}

}