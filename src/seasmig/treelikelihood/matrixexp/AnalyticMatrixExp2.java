package seasmig.treelikelihood.matrixexp;

import org.junit.Test;

import seasmig.treelikelihood.*;
import seasmig.util.Util;

@SuppressWarnings("serial")
public class AnalyticMatrixExp2 implements MatrixExponentiator {

	// Q has to be a proper rate matrix with Q[i][i] = -Sum of rest of row i

	double[][] Q = new double[2][2]; 
	double denum;
	boolean methodOK;

	protected AnalyticMatrixExp2() {};

	public AnalyticMatrixExp2(double[][] Q) {
		methodOK = (Q.length==2);
		if (methodOK) 
			methodOK = (Q[0].length==2);
		if (methodOK) {  
			for (int i=0;i<Q.length;i++) {
				for (int j=0;j<Q.length;j++) {
					this.Q[i][j]=Q[i][j];
				}
			}
			denum=Q[0][1]+Q[1][0]+Util.minValue;
			methodOK=(denum>0);
		}
	}

	@Override
	public double[][] expm(double t) {

		double[][] returnValue = new double[2][2];

		if (denum!=0) {
			double exp = Math.exp(t*(-Q[0][1]-Q[1][0]));			
			returnValue[0][0]=(exp*Q[0][1]+Q[1][0])/denum;
			returnValue[0][1]=(1-exp)*Q[0][1]/denum;
			returnValue[1][0]=(1-exp)*Q[1][0]/denum;
			returnValue[1][1]=(exp*Q[1][0]+Q[0][1])/denum;
			return returnValue;
		} 

		return returnValue;		
	}

	@Test
	public void test()
	{
		// TODO: This...
		MatrixExponentiator matrixExponentiator1 = new AnalyticMatrixExp2(new double[][]{{-3,3},{0.5,-0.5}});	
		MatrixExponentiator matrixExponentiator2 = new Matlab7MatrixExp(new double[][]{{-3,3},{0.5,-0.5}});
		double[][] res1=matrixExponentiator1.expm(0.1);
		double[][] res2=matrixExponentiator2.expm(0.1);

		System.out.println("res1:");
		for (int i=0;i<res1.length;i++) {
			for (int j=0;j<res1[0].length;j++) {
				System.out.print(res1[i][j]+"\t");				
			}
			System.out.println();
		}
		System.out.println("res2:");
		for (int i=0;i<res2.length;i++) {
			for (int j=0;j<res2[0].length;j++) {
				System.out.print(res2[i][j]+"\t");								
			}
			System.out.println();
		}

		System.out.println("timing:");
		long startTime1= System.currentTimeMillis();	
		for (int rep=0;rep<1000000;rep++) {
			res1=matrixExponentiator1.expm(rep/10000);			
			if (Math.random()<0.00000001) {
				System.out.println(res1[0][0]);
			}
		}
		long time1= System.currentTimeMillis()-startTime1;

		long startTime2= System.currentTimeMillis();
		for (int rep=0;rep<1000000;rep++) {
			res2=matrixExponentiator2.expm(rep/10000);			
			if (Math.random()<0.00000001) {
				System.out.println(res2[0][0]);
			}
		}
		long time2= System.currentTimeMillis()-startTime2;
		System.out.println("time1: "+time1+"[ms] time2: "+time2+"[ms]");	
	}

	@Override
	public boolean checkMethod() {		
		return methodOK;
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
