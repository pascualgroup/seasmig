package seasmig.util;

import java.util.Arrays;
import java.util.Comparator;

import cern.jet.random.engine.RandomEngine;

public class Util {

	public static final double minValue = Double.MIN_VALUE;
	public static final double minNegative = Double.NEGATIVE_INFINITY;

	protected Util() {};

	public static double logSumExp(double[] alphas, double min) {
		// TODO: improve this
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
				returnValue+=String.format("%f",Q[i][j]);
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
	
//	private static double sum(double[] vector) {
//		double returnValue=0;
//		for (int i=0;i<vector.length;i++){
//			returnValue+=vector[i];			
//		}
//		return returnValue;
//	}
//
//	private static double sum(double[][] t) {
//		int n = t.length;
//		int m = t[0].length;
//		double returnValue=0;
//		for (int i=0;i<n;i++){
//			for (int j=0;j<m;j++){
//				returnValue+=t[i][j];
//			}
//		}
//		return returnValue;
//	}

//	private static double[][] mul(double[][] a, double[][] b) {
//		int n = a.length;
//		int m = a[0].length;
//		double[][] returnValue = new double[n][m];
//		for (int i=0;i<n;i++){
//			for (int j=0;j<m;j++){
//				returnValue[i][j]=a[i][j]*b[i][j];
//			}
//		}
//		return returnValue;
//	}
	
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

	public static class ArrayIndexComparator implements Comparator<Integer>
	{
	    private final double[] array;

	    public ArrayIndexComparator(double[] array)
	    {
	        this.array = array;
	    }

	    public Integer[] createIndexArray()
	    {
	        Integer[] indexes = new Integer[array.length];
	        for (int i = 0; i < array.length; i++)
	        {
	            indexes[i] = i; // Autoboxing
	        }
	        return indexes;
	    }

	    @Override
	    public int compare(Integer index1, Integer index2)
	    {
	         // Autounbox from Integer to int to use as array indexes
	        return Double.compare(array[index1],array[index2]);
	    }
	}
	
	public static double[] toDirichletDegenerate(double[] unif) {
		// TODO: test this
		ArrayIndexComparator comparator = new ArrayIndexComparator(unif);
		Integer[] indexes = comparator.createIndexArray();	
		Arrays.sort(indexes, comparator);
		double[] returnValue = new double[unif.length];
		returnValue[0]=unif[indexes[0]]-0;		
		for (int i=1;i<indexes.length;i++) {
			returnValue[i]=unif[indexes[i]]-unif[indexes[i-1]];		
		}			
		
		return returnValue; 
	}
	
	public static double[] toDirichletNonDegenerate(double[] unif) {  // 3 states only
		// TODO: test this
		// converts a vector of n i.i.d uniformly distributed U(0,1) variables to n Dirichlet distributed variables the last (1-x0-x1-x2...) is not included
		// The order of the original variables remains intact in order to prevent multiple mappings to the same realization
		// (I haven't found an online example of this (the order part), so should be careful about this part of the code)
		if (unif.length!=3) {
			System.err.println("non degenerate conversion to Dirichlet only supported for 3 states\n");
			System.exit(0);
		}
		ArrayIndexComparator comparator = new ArrayIndexComparator(unif);
		Integer[] indexes = comparator.createIndexArray();	
		Arrays.sort(indexes, comparator);
		double[] sorted = new double[unif.length+1];
		sorted[0]=unif[indexes[0]]-0;
		double sum=sorted[0];
		for (int i=1;i<indexes.length;i++) {
			sorted[i]=unif[indexes[i]]-unif[indexes[i-1]];
			sum=sum+sorted[i];
		}
		sorted[indexes.length]=1-sum;
		// Now sorted is n+1 variables from a Dirichlet distribution x0+x1+x2+..+xn = 1 x0>0, x1>0 ... xn>0
		// I would return sorted now but I also want a 1-1 mapping
		
		// Which order to plug them into the result, to keep a one-on-one mapping
		// this is the dodgy part:				
		
		// lets say n = 3 so we have x0,x1,x2
		// which in order are xa,xb,xc
		// and r0 = xa-0, r1 = xb-xa, r2 = xc-xb, r3 = 1 -xc		
		// there are 3! orderings of r
		// We still have 4 (n) options of which of the r0,r1,r2,r3 (...rn) to exclude
		// there is a degeneracy xc could be closer to 1 or to xb with the same result
		// there is a degeneracy xa could be closer to 0 or to xb with the same result 

		int itemToDrop = ((1-unif[indexes[2]]) > (unif[indexes[2]]-unif[indexes[1]])?2:0) +
				         ((unif[indexes[0]]-0) > (unif[indexes[1]]-unif[indexes[0]])?1:0);
		sorted[itemToDrop]=2; // will have > 1 value so will be max
		Arrays.sort(sorted);
		double[] returnValue = new double[unif.length];
		for (int i=0;i<indexes.length;i++) {
			returnValue[indexes[i]]=sorted[i];
		}
		return returnValue;
	}		

}
