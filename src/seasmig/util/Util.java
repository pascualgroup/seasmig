package seasmig.util;

import cern.jet.random.engine.RandomEngine;

public class Util {

	public static final double minValue = Double.MIN_VALUE;
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	protected Util() {};

	public static double logSumExp(double[] alphas, double min) {
		double sumExp = 0;		
		for (int i=0;i<alphas.length;i++) {			
			sumExp=sumExp+Math.exp(alphas[i]-min);
		}
		double returnValue=min+Math.log(sumExp);
		if (Double.isNaN(returnValue) ){
			return minNegative;
		}
		return returnValue;

	}

	static public String print(double[][] Q) {		
		String returnValue = "{";
		for (int i=0; i<Q.length;i++) {
			if (i!=0) 
				returnValue+=" ";
			returnValue+="{";
			for (int j=0; j<Q.length;j++) {
				returnValue+=String.format("%6.2f",Q[i][j]);
				if (j!=Q.length-1) {
					returnValue+=",";
				}
			}			
			returnValue+="}";
			if (i!=Q.length-1) {
				returnValue+=",\n";
			}			
		}
		returnValue+="}\n";
		return returnValue;
	}

	static public String parse(double[][] Q) {		
		String returnValue = "{";
		for (int i=0; i<Q.length;i++) {
			if (i!=0) 
				returnValue+=" ";
			returnValue+="{";
			for (int j=0; j<Q.length;j++) {
				returnValue+=String.format("%6.2f",Q[i][j]);
				if (j!=Q.length-1) {
					returnValue+=",";
				}
			}			
			returnValue+="}";
			if (i!=Q.length-1) {
				returnValue+=",";
			}			
		}
		returnValue+="}";
		return returnValue;
	}

//	def calc_nuc_q_matrix(rmatrix,basefreq):
//	bigpi = ones((4,4))
//	count = 0
//	for i in bigpi:
//		count2 = 0
//		for j in i:
//			if count != count2:
//				bigpi[count][count2] = basefreq[count] * basefreq[count2]
//			else:
//				bigpi[count][count2] = basefreq[count]
//			count2 += 1
//		count += 1
//	t = rmatrix * bigpi
//	count = 0
//	fill_diagonal(t,0)
//	tscale = sum(t)
//	t= t/tscale
//	#make it so the diags make the rows sum to 0
//	for i in t:
//		t[count][count] = 0-sum(i)
//		count += 1
//	t = t/basefreq
//	t = transpose(t)
//	return t
	
	static public double[][] calcQMatrix(double[][] rMatrix, double[] basefreq) {
		// TODO: check this
		int n = rMatrix.length;
		
		double[][] piMatrix = ones(n,n);				
		for (int i=0;i<n;i++) {
			for (int j=0;j<n;j++) {
				if (i!=j) {
					piMatrix[i][j]=basefreq[i]*basefreq[j];
				}
				else {
					piMatrix[i][j]=basefreq[i];	
				}
			}
		}
		
		double[][] t = mul(rMatrix, piMatrix);
	
		//	fillDiagonal of t with zeros;
		for (int i=0;i<n;i++)
			t[i][i]=0;
					
		double tscale = sum(t);
		
		// t=t/tscale
		for (int i=0;i<n;i++) {
			for (int j=0;j<n;j++) {
				t[i][j]=t[i][j]/tscale;
			}
		}
		
		// make it so the diags make the rows sum to 0
		for (int i=0;i<n;i++) {
			t[i][i]=0.0-sum(t[i]);
		}
		
		// t = t / basefreq
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				t[i][j]=t[i][j]/basefreq[j];
			}
		}
	
		transposeSquareMatrix(t);
		return t; 		
	}
	

	private static double sum(double[] vector) {
		double returnValue=0;
		for (int i=0;i<vector.length;i++){
			returnValue+=vector[i];			
		}
		return returnValue;
	}

	private static double sum(double[][] t) {
		int n = t.length;
		int m = t[0].length;
		double returnValue=0;
		for (int i=0;i<n;i++){
			for (int j=0;j<m;j++){
				returnValue+=t[i][j];
			}
		}
		return returnValue;
	}

	private static double[][] mul(double[][] a, double[][] b) {
		int n = a.length;
		int m = a[0].length;
		double[][] returnValue = new double[n][m];
		for (int i=0;i<n;i++){
			for (int j=0;j<m;j++){
				returnValue[i][j]=a[i][j]*b[i][j];
			}
		}
		return returnValue;
	}


	static public double[][] ones(int n, int m) {
		double[][] returnValue = new double[n][m];
		for (int i=0;i<n;i++){
			for (int j=0;j<m;j++){
				returnValue[i][j]=1.0;
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

	public static int nextIntFromToExcept(RandomEngine rng, int min, int max, int except) {
		int newValue =  min + (int)(rng.nextDouble() * ((max-1 - min) + 1));
		if (newValue==except) {
			newValue=max;
		}
		return newValue;
	}

	public static void transposeSquareMatrix(double[][] mat) {	
		int n = mat.length;		
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				if (i!=j) {
					double temp=mat[i][j];
					mat[i][j]=mat[j][i];
					mat[j][i]=temp;
				}
			}
		}			
	}


}
