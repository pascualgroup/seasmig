package treelikelihood;

import java.util.Vector;

import org.javatuples.Pair;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;


public class MatlabMatrixExp implements MatrixExponentiator {

	// Cache 
	DoubleMatrix2D Q;

	DoubleFactory2D F = DoubleFactory2D.dense; 

	public MatlabMatrixExp(DoubleMatrix2D Q_) {
		Q = Q_;		
	}
	
	public DoubleMatrix2D eye() {
		return F.dense.identity(Q.rows());
	}
	
	public double norm_inf() {
		double returnValue = Double.NEGATIVE_INFINITY;
		for (int i=0; i<Q.rows();i++) {
			double rowSum=0;
			for (int j=0; j<Q.rows();j++) {
				rowSum=rowSum+Math.abs(Q.get(i, j));
			}
			if (rowSum>returnValue) {
				returnValue=rowSum;
			}
		}
		return returnValue;
	}
	
	
	
	public static FRexpResult log2(double value)
	{
	   final FRexpResult result = new FRexpResult();
	   long bits = Double.doubleToLongBits(value);
	   double realMant = 1.;

	   // Test for NaN, infinity, and zero.
	   if (Double.isNaN(value) || value + value == value || Double.isInfinite(value))   {
	      result.e = 0;
	      result.f = value;
	   }
	   else  {

	      boolean neg = (bits < 0);
	      int exponent = (int)((bits >> 52) & 0x7ffL);
	      long mantissa = bits & 0xfffffffffffffL;

	      if(exponent == 0) {
	         exponent++;
	      }
	      else  {
	         mantissa = mantissa | (1L<<52);
	      }

	      // bias the exponent - actually biased by 1023.
	      // we are treating the mantissa as m.0 instead of 0.m
	      //  so subtract another 52.
	      exponent -= 1075;
	      realMant = mantissa;

	      // normalize
	      while(realMant > 1.0) {
	         mantissa >>= 1;
	         realMant /= 2.;
	         exponent++;
	      }

	      if(neg) {
	         realMant = realMant * -1;
	      }

	      result.e = exponent;
	      result.f = realMant;
	   }
	   return result;
	}
	
	public static class FRexpResult
	{
	   public int e = 0;
	   public double f = 0.;
	}

	@Override
	public DoubleMatrix2D expm(double t) {
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
		DoubleMatrix2D A = Q.copy();
		MatlabMatrixExp.FRexpResult fe = log2(norm_inf()); // [ f, e ] = log2 ( norm ( A, 'inf' ) );
		long s = Math.max(0, fe.e+1); //  s = max ( 0, e + 1 );
		A.assign(cern.jet.math.tdouble.DoubleFunctions.div(1L<<s)); // A = A / 2^s;
		
		DoubleMatrix2D X = A.copy();
		double c = 0.5;
		DoubleMatrix2D E = A.copy().assign(eye(),cern.jet.math.tdouble.DoublePlusMultFirst.plusMult(c)); // I + c * A;
		DoubleMatrix2D D = A.copy().assign(eye(),cern.jet.math.tdouble.DoublePlusMultFirst.minusMult(c)); // I - c * A;
		double q = 6.0;
		boolean p = true;
		for (int k = 2;k<=q;k++) {
			    c = c*(q-k+1)/(k*(2*q-k+1));
			    X.zMult(A, null); // X = A * X;
			    DoubleMatrix2D cX=X.copy().assign(cern.jet.math.tdouble.DoubleFunctions.mult(c)); // cX = c * X; 
			    E.assign(X,cern.jet.math.tdouble.DoublePlusMultSecond.plusMult(c)); //E = E + cX;
			    if ( p ) {
			      D.assign(X,cern.jet.math.tdouble.DoublePlusMultSecond.plusMult(c)); // D = D + cX;
			    }
			    else {
			      D.assign(X,cern.jet.math.tdouble.DoublePlusMultSecond.minusMult(c)); // D = D - cX;
			    }
			    
			    p = !p;
		}
		
		// TODO: check if can use E.zMult(E,E);
		DenseDoubleAlgebra myAlgebra = new DenseDoubleAlgebra();
		E=myAlgebra.inverse(D).zMult(E, null); // E = D \ E;

		for (int k=1;k<=s;k++) {
			E=E.zMult(E, null); // E = E*E 
		} 
		 		
		return E;
	}

	

}
