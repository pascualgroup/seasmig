package treelikelihood;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.algo.DenseDoubleAlgebra;
import cern.jet.math.tdouble.DoublePlusMultSecond;
import cern.jet.math.tdouble.DoubleFunctions;

/*
function F = expm(A)
		%EXPM   Matrix exponential.
		%   EXPM(X) is the matrix exponential of X.  EXPM is computed using
		%   a scaling and squaring algorithm with a Pade approximation.
		%
		%   Although it is not computed this way, if X has a full set
		%   of eigenvectors V with corresponding eigenvalues D, then
		%   [V,D] = EIG(X) and EXPM(X) = V*diag(exp(diag(D)))/V.
		%
		%   EXP(X) computes the exponential of X element-by-element.
		%
		%   See also LOGM, SQRTM, FUNM.

		%   Reference:
		%   N. J. Higham, The scaling and squaring method for the matrix
		%   exponential revisited. SIAM J. Matrix Anal. Appl.,
		%   26(4) (2005), pp. 1179-1193.
		%
		%   Nicholas J. Higham
		%   Copyright 1984-2005 The MathWorks, Inc.
		%   $Revision: 5.10.4.6 $  $Date: 2005/11/18 14:15:53 $

		[m_vals, theta, classA] = expmchk; % Initialization
		normA = norm(A,1);

		if normA <= theta(end)
		    % no scaling and squaring is required.
		    for i = 1:length(m_vals)
		        if normA <= theta(i)
		            F = PadeApproximantOfDegree(m_vals(i));
		            break;
		        end
		    end
		else

		    [t s] = log2(normA/theta(end));
		    s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
		    A = A/2^s;    % Scaling
		    F = PadeApproximantOfDegree(m_vals(end));
		    for i = 1:s
		        F = F*F;  % Squaring
		    end
		end
		% End of expm


/*
 * function c = getPadeCoefficients
		            % GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
		            %    C = GETPADECOEFFICIENTS returns coefficients of numerator
		            %    of [M/M] Pade approximant, where M = 3,5,7,9,13.
		            switch m
		                case 3
		                    c = [120, 60, 12, 1];
		                case 5
		                    c = [30240, 15120, 3360, 420, 30, 1];
		                case 7
		                    c = [17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1];
		                case 9
		                    c = [17643225600, 8821612800, 2075673600, 302702400, 30270240, ...
		                         2162160, 110880, 3960, 90, 1];
		                case 13
		                    c = [64764752532480000, 32382376266240000, 7771770303897600, ...
		                         1187353796428800,  129060195264000,   10559470521600, ...
		                         670442572800,      33522128640,       1323241920,...
		                         40840800,          960960,            16380,  182,  1];
		            end
		        end
 */


/*
function F = PadeApproximantOfDegree(m)
        %PADEAPPROXIMANTOFDEGREE  Pade approximant to exponential.
        %   F = PADEAPPROXIMANTOFDEGREE(M) is the degree M diagonal
        %   Pade approximant to EXP(A), where M = 3, 5, 7, 9 or 13.
        %   Series are evaluated in decreasing order of powers, which is
        %   in approx. increasing order of maximum norms of the terms.

        n = length(A);
        c = getPadeCoefficients;

        % Evaluate Pade approximant.
        switch m

            case {3, 5, 7, 9}

                Apowers = cell(ceil((m+1)/2),1);
                Apowers{1} = eye(n,classA);
                Apowers{2} = A*A;
                for j = 3:ceil((m+1)/2)
                    Apowers{j} = Apowers{j-1}*Apowers{2};
                end
                U = zeros(n,classA); V = zeros(n,classA);

                for j = m+1:-2:2
                    U = U + c(j)*Apowers{j/2};
                end
                U = A*U;
                for j = m:-2:1
                    V = V + c(j)*Apowers{(j+1)/2};
                end
                F = (-U+V)\(U+V);

            case 13

                % For optimal evaluation need different formula for m >= 12.
                A2 = A*A; A4 = A2*A2; A6 = A2*A4;
                U = A * (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) ...
                    + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA) );
                V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) ...
                    + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA);
                F = (-U+V)\(U+V);
        end
    end
 */


public class Matlab7MatrixExp implements MatrixExponentiator {

	// Pade Coefficients m based ("not zero based")
	static final long[][] c = {null,null,null,
		{120L, 60L, 12L, 1L} /*c3*/,
		null,
		{30240L, 15120L, 3360L, 420L, 30L, 1L} /*c5*/,
		null,
		{17297280L, 8648640L, 1995840L, 277200L, 25200L, 1512L, 56L, 1L} /*c7*/,
		null,
		{17643225600L, 8821612800L, 2075673600L, 302702400L, 30270240L, 
			2162160L, 110880L, 3960L, 90L, 1L} /*c9*/,
			null, null, null,
			{64764752532480000L, 32382376266240000L, 7771770303897600L, 
				1187353796428800L,  129060195264000L,   10559470521600L, 
				670442572800L,      33522128640L,       1323241920L,
				40840800L,          960960L,            16380L,  182L,  1L} /*c13*/}; 

	// Assumes double precision
	int[] m_vals = new int[]{3,5,7,9,13};
	double[] theta = new double[]{1.495585217958292e-002,  // m_vals = 3
			2.539398330063230e-001,  // m_vals = 5
			9.504178996162932e-001,  // m_vals = 7
			2.097847961257068e+000,  // m_vals = 9
			5.371920351148152e+000}; // m_vals = 13


	DoubleMatrix2D Q;
	DenseDoubleAlgebra algebra = new DenseDoubleAlgebra(Double.MIN_VALUE);

	private DoubleMatrix2D eye;

	public Matlab7MatrixExp(DoubleMatrix2D Q_) {
		Q = Q_;		
		eye = DoubleFactory2D.dense.identity(Q.rows());
	}

	
	DoubleMatrix2D padeApproximantOfDegree(int m, DoubleMatrix2D A) {
		// TODO: use Q powers to quickly calculate A powers when precision allows!	
		
		//	c = getPadeCoefficients (in constructor)
		DoubleMatrix2D[] Apowers; 
		DoubleMatrix2D U = DoubleFactory2D.dense.make(A.rows(),A.rows());
		DoubleMatrix2D V = DoubleFactory2D.dense.make(A.rows(),A.rows());
		DoubleMatrix2D F;

		// Evaluate Pade approximant.
		switch (m) {

		case 3: case 5: case 7: case 9:
			//Apowers = cell(ceil((m+1)/2),1);
			Apowers = new DoubleMatrix2D[(m+1)/2]; 
			// Apowers{1} = eye(n,classA);
			Apowers[0]=eye;
			// Apowers{2} = A*A;
			Apowers[1]=A.zMult(A,null);
			for (int j=2;j<(m+1)/2;j++) { // for j = 3:ceil((m+1)/2)
				Apowers[j]=Apowers[j-1].zMult(Apowers[1], null); // Apowers{j} = Apowers{j-1}*Apowers{2};
			}
			// for j = m+1:-2:2
			for (int j=m;j>=1;j-=2) { // convert j to zero based index
				U.assign(Apowers[(j-1)/2],DoublePlusMultSecond.plusMult(c[m][j])); //  U = U + c(j)*Apowers{j/2};
			}
			U=A.zMult(U, null); // U = A*U;
			
			// for j = m:-2:1
			for (int j=(m-1);j>=0;j-=2) {
				//				V = V + c(j)*Apowers{(j+1)/2};
				V.assign(Apowers[j/2],DoublePlusMultSecond.plusMult(c[m][j]));
			}

			break;
		case 13: 
			// % For optimal evaluation need different formula for m >= 12.
			DoubleMatrix2D A2 = A.zMult(A, null); 	// A2 = A*A; 
			DoubleMatrix2D A4 = A2.zMult(A2, null);
			DoubleMatrix2D A6 = A4.zMult(A2, null);
			//  U = A * 
			// (A6*(c(14)*A6 + c(12)*A4 + c(10)*A2) 
			//  + c(8)*A6 + c(6)*A4 + c(4)*A2 + c(2)*eye(n,classA)  
			// );
			U = A.zMult(A6.zMult(zSum3(A6,A4,A2,c[m][13],c[m][11],c[m][9]),null).assign(zSum4(A6,A4,A2,eye,c[m][7],c[m][5],c[m][3],c[m][1]),DoubleFunctions.plus),null);

			//	V = A6*(c(13)*A6 + c(11)*A4 + c(9)*A2) 
			//  + c(7)*A6 + c(5)*A4 + c(3)*A2 + c(1)*eye(n,classA);
			V = A6.zMult(zSum3(A6,A4,A2,c[m][12],c[m][10],c[m][8]),null).assign(zSum4(A6,A4,A2,eye,c[m][6],c[m][4],c[m][2],c[m][0]),DoubleFunctions.plus);	                    
		}
		// TODO: optimize this
		DoubleMatrix2D VplusU = V.copy().assign(U,DoubleFunctions.plus); 
		DoubleMatrix2D VminusU = V.assign(U,DoubleFunctions.minus);
		F = algebra.inverse(VminusU).zMult(VplusU,null); // F = (-U+V)\(U+V);
		return F;
	}
	
	DoubleMatrix2D zSum3(DoubleMatrix2D A, DoubleMatrix2D B,DoubleMatrix2D C,double alpha,double beta, double gamma) {
		DoubleMatrix2D returnValue = A.copy().assign(DoubleFunctions.mult(alpha));
		returnValue.assign(B,DoublePlusMultSecond.plusMult(beta));
		returnValue.assign(C,DoublePlusMultSecond.plusMult(gamma));
		return returnValue;
	}
	
	DoubleMatrix2D zSum4(DoubleMatrix2D A, DoubleMatrix2D B,DoubleMatrix2D C,DoubleMatrix2D D, double alpha,double beta, double gamma, double delta) {
		DoubleMatrix2D returnValue = A.copy().assign(DoubleFunctions.mult(alpha));
		returnValue.assign(B,DoublePlusMultSecond.plusMult(beta));
		returnValue.assign(C,DoublePlusMultSecond.plusMult(gamma));
		returnValue.assign(D,DoublePlusMultSecond.plusMult(delta));
		return returnValue;
	}

	@Override
	public DoubleMatrix2D expm(double tt) {
	
		//	Initialization is in constructor [m_vals, theta, classA=='double'] = expmchk;
		DoubleMatrix2D A = Q.copy().assign(DoubleFunctions.mult(tt));
		double normA = algebra.norm1(A);
		DoubleMatrix2D F=null;

		if (normA <= theta[theta.length-1]) { //		if normA <= theta(end)
			//	no scaling and squaring is required.
			for (int i=0;i<m_vals.length;i++) {
				if (normA <= theta[i]) {
					F = padeApproximantOfDegree(m_vals[i],A);
					break;
				}
			}
		}
		else {
			FRexpResult ts = log2(normA/theta[theta.length-1]); 	// [t s] = log2(normA/theta(end));
			int s = ts.e;
			double t = ts.f;
			if (t==0.5) s = s - 1; // s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
			A.assign(DoubleFunctions.div(1L<<s)); // A = A/2^s;    % Scaling
			F = padeApproximantOfDegree(m_vals[m_vals.length-1],A);
			for (int i=0;i<s;i++) {
				F=F.zMult(F, null); // F = F*F;  % Squaring
			}
		}
		return F;
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

}
