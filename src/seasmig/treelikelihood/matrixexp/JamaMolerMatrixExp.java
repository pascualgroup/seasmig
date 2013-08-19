package seasmig.treelikelihood.matrixexp;

import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.util.Util;
import seasmig.util.Util.FRexpResult;
import Jama.Matrix;


@SuppressWarnings("serial")
public class JamaMolerMatrixExp implements MatrixExponentiator {

	double[][] Q;	
	Matrix eye; 
	static final double q = 6.0;
	
	protected JamaMolerMatrixExp() {};
	
	public JamaMolerMatrixExp(double[][] Q_) {
		Q = Q_;	
		eye = Matrix.identity(Q.length,Q.length);
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
		
		Matrix A = new Matrix(Q).times(t);
		FRexpResult fe = Util.log2(A.normInf()); // [ f, e ] = log2 ( norm ( A, 'inf' ) );
		long s = Math.max(0, fe.e+1); //  s = max ( 0, e + 1 );
		A=A.times(1.0/((long)1L<<s)); // A = A / 2^s;

		Matrix X = A.copy(); // X = A
		double c = 0.5;
		Matrix cA = A.times(c);
		Matrix E = eye.plus(cA); // I + c * A; 
		Matrix D = eye.minus(cA); // I - c * A; 
		boolean p = true;
		for (int k = 2;k<=q;k++) {
			c = c*(q-k+1)/(k*(2*q-k+1));
			X=A.times(X); // X = A * X;
			Matrix cX=X.times(c);
			E=E.plus(cX); //E = E + cX;
			if ( p ) 
				D=D.plus(cX); // D = D + cX;
			else 
				D=D.minus(cX); // D = D - cX;	    
			p = !p;
		}	
		E=D.inverse().times(E); // E = D \ E;
		for (int k=1;k<=s;k++) 
			E=E.times(E); // E = E*E 
		return E.getArray();
	}

	@Override
	public boolean checkMethod() {
		// TODO: this....
		return true;
	}

	@Override
	public String getMethodName() {		
		Class<?> enclosingClass = getClass().getEnclosingClass();
		if (enclosingClass != null) 
		    return enclosingClass.getName();
		else 
		    return getClass().getName();
		
	}

}
