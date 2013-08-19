package seasmig.treelikelihood.matrixexp;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.util.Util;
import seasmig.util.Util.FRexpResult;

// TODO: Failes Matrix Exp Test !!!
//Matrix exponentiation results fail to match!
//t=51.00000000000001 
//Q=
//{{ -0.01,  0.01,  0.00},
//{  0.01, -0.02,  0.01},
//{  0.00,  0.00, -0.01}}
//
//MolerMatrixExp - checkMethod: true
//{{  0.37,  0.28,  0.35},
//{  0.34,  0.27,  0.38},
//{  0.30,  0.25,  0.46}}
//
//Matlab7MatrixExp - checkMethod: true
//{{  0.63,  0.23,  0.14},
//{  0.27,  0.49,  0.24},
//{  0.13,  0.14,  0.73}}
//
//TaylorMatrixExp - checkMethod: true
//{{  0.63,  0.23,  0.14},
//{  0.27,  0.49,  0.24},
//{  0.13,  0.14,  0.73}}
//
//JamaMolerMatrixExp: - checkMethod: true
//{{  0.63,  0.23,  0.14},
//{  0.27,  0.49,  0.24},
//{  0.13,  0.14,  0.73}}
//
//JblasMatrixExp: - checkMethod: true
//{{  0.63,  0.23,  0.14},
//{  0.27,  0.49,  0.24},
//{  0.13,  0.14,  0.73}}
//
//AnalyticMatrixExp3: - checkMethod: true
//{{  0.63,  0.23,  0.14},
//{  0.27,  0.49,  0.24},
//{  0.13,  0.14,  0.73}}

@SuppressWarnings("serial")
public class MolerMatrixExp implements MatrixExponentiator {

	double[][] Q;
	Algebra algebra = new Algebra();
	DoubleMatrix2D eye;
	private double normInfQ; 
	static final double q = 6.0;
	
	protected MolerMatrixExp() {};
	
	public MolerMatrixExp(double[][] testMatrix) {
		Q = testMatrix;	
		eye = DoubleFactory2D.dense.identity(Q.length);
		normInfQ=norm_inf();		
	}
	
	public DoubleMatrix2D eye() {
		return DoubleFactory2D.dense.identity(Q.length); // TODO: change to dynamic call...
	}
	
	public double norm_inf() {
		double returnValue = Double.NEGATIVE_INFINITY;
		for (int i=0; i<Q.length;i++) {
			double rowSum=0;
			for (int j=0; j<Q.length;j++) {
				rowSum=rowSum+Math.abs(Q[i][j]);
			}
			if (rowSum>returnValue) {
				returnValue=rowSum;
			}
		}
		return returnValue;
	}
	
	/*
	function E = expm1 ( A )

			%*****************************************************************************
			% EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
			%
			%  Licensing:    This code is distributed under the GNU LGPL license.
			%  Modified:    18 March 2011
			%  Author:    Cleve Moler, Charles Van Loan
			%  Reference: 
			%    Cleve Moler, Charles VanLoan,
			%    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
			%    Twenty-Five Years Later, SIAM Review, Volume 45, Number 1, March 2003, pages 3-49.
			%  Parameters:
			%    Input, real A(N,N), the matrix.
			%    Output, real E(N,N), the estimate for exp(A).
			%
			  [ f, e ] = log2 ( norm ( A, 'inf' ) );
			  s = max ( 0, e + 1 );
			  A = A / 2^s;

			  X = A; 
			  c = 1 / 2;
			  E = eye ( size ( A ) ) + c * A;
			  D = eye ( size ( A ) ) - c * A;
			  q = 6;
			  p = 1;

			  for k = 2 :  q
			    c = c * ( q - k + 1 ) / ( k * ( 2 * q - k + 1 ) );
			    X = A * X;
			    cX = c * X;
			    E = E + cX;
			    if ( p )
			      D = D + cX;
			    else
			      D = D - cX;
			    end
			    p = ~p;
			  end

			  E = D \ E;

			  for k = 1 : s
			    E = E * E;
			  end

			  return
			end
	 */
	@Override
	public double[][] expm(double t) {
		
		DoubleMatrix2D A = DoubleFactory2D.dense.make(Q).assign(cern.jet.math.Functions.mult(t));
		FRexpResult fe = Util.log2(normInfQ*t); // [ f, e ] = log2 ( norm ( A, 'inf' ) );
		long s = Math.max(0, fe.e+1); //  s = max ( 0, e + 1 );
		A.assign(cern.jet.math.Functions.div(1L<<s)); // A = A / 2^s;

		DoubleMatrix2D X = A.copy(); // X = A
		double c = 0.5;
		DoubleMatrix2D E = A.copy().assign(eye,cern.jet.math.Functions.plusMult(c)); // I + c * A; 
		DoubleMatrix2D D = A.copy().assign(eye,cern.jet.math.Functions.minusMult(c)); // I - c * A; 
		boolean p = true;
		for (int k = 2;k<=q;k++) {
			c = c*(q-k+1)/(k*(2*q-k+1));
			X=X.zMult(A, null); // X = A * X;			    
			E.assign(X,cern.jet.math.Functions.plusMult(c)); //E = E + cX;
			if ( p ) 
				D.assign(X,cern.jet.math.Functions.plusMult(c)); // D = D + cX;
			else 
				D.assign(X,cern.jet.math.Functions.minusMult(c)); // D = D - cX;	    
			p = !p;
		}	
		E=algebra.inverse(D).zMult(E, null); // E = D \ E;
		for (int k=1;k<=s;k++) 
			E=E.zMult(E, null); // E = E*E 
		System.err.println("Method Failes MatrixExp Test... Exiting\n");
		System.exit(-1); 
		return E.toArray();
	}

	@Override
	public boolean checkMethod() {
		// TODO: this....
		return false;
	}

}
