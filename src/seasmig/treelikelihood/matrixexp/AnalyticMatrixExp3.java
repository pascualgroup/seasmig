package seasmig.treelikelihood.matrixexp;

import org.junit.Test;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import seasmig.treelikelihood.MatrixExponentiator;

@SuppressWarnings("serial")
public class AnalyticMatrixExp3 implements MatrixExponentiator {

	// Q has to be a proper rate matrix with Q[i][i] = -Sum of rest of row i

	double[][] Q = new double[3][3];
	double Lambda;
	double Delta; 
	double Epsilon;

	boolean methodOK;
	private double denum1;
	private double denum3;

	protected AnalyticMatrixExp3() {};

	public AnalyticMatrixExp3(double[][] Q) {

		methodOK = (Q.length==3);
		if (methodOK) 
			methodOK = (Q[0].length==3);

		if (methodOK) {
			for (int i=0;i<Q.length;i++) {
				for (int j=0;j<Q.length;j++) {
					this.Q[i][j]=Q[i][j];
				}
			}

			Lambda = Q[0][1]+Q[0][2]+Q[1][0]+Q[1][2]+Q[2][0]+Q[2][1];

			Delta = Q[0][2]*Q[1][0]+Q[0][1]*Q[1][2]+Q[0][2]*Q[1][2]+
					Q[0][1]*Q[2][0]+Q[1][0]*Q[2][0]+Q[1][2]*Q[2][0]+
					Q[0][1]*Q[2][1]+Q[0][2]*Q[2][1]+Q[1][0]*Q[2][1];

			Epsilon = Math.sqrt(-4*Delta+Lambda*Lambda);

			denum1=4*Delta+(Epsilon-Lambda)*Lambda;
			denum3=-4*Delta+(Epsilon+Lambda)*Lambda;

			methodOK = (denum1>0) && (denum3>0) && !Double.isNaN(Epsilon);
		}

	}

	@Override
	public DoubleMatrix2D expm(double t) {
		double[][] returnValue = new double[3][3];


		if (denum1!=0 && denum3!=0 && Delta!=0) {
			double exp1 = Math.exp(0.5*t*(Epsilon-Lambda));			
			double exp3 = Math.exp(-0.5*t*(Epsilon+Lambda));

			for (int i=0;i<3;i++) {	
				returnValue[i][i]=
						exp1*(-2*Q[i][(i+1)%3]*Q[i][(i+1)%3]+Q[i][(i+1)%3]*(Epsilon+Lambda-4*Q[i][(i+2)%3]-2*Q[(i+1)%3][i])+Q[i][(i+2)%3]*(Epsilon+Lambda-2*Q[i][(i+2)%3]-2*Q[(i+2)%3][i]))/denum1+
						(Delta+Q[i][(i+1)%3]*Q[i][(i+1)%3]+Q[i][(i+2)%3]*Q[i][(i+2)%3]+Q[i][(i+1)%3]*(-Lambda+2*Q[i][(i+2)%3]+Q[(i+1)%3][i])+Q[i][(i+2)%3]*(-Lambda+Q[(i+2)%3][i]))/Delta+
						exp3*(2*Q[i][(i+1)%3]*Q[i][(i+1)%3]+Q[i][(i+1)%3]*(Epsilon-Lambda+4*Q[i][(i+2)%3]+2*Q[(i+1)%3][i])+Q[i][(i+2)%3]*(Epsilon-Lambda+2*Q[i][(i+2)%3]+2*Q[(i+2)%3][i]))/denum3;
			}

			for (int i=0;i<3;i++) {	
				for (int j=0;j<3;j++) {
					if (i==j) continue;
					int notij=3-(i|j);					
					returnValue[i][j]=
							-exp1*(Q[i][j]*(Epsilon+Lambda-2*Q[i][j]-2*Q[j][i]-2*Q[j][notij]-2*Q[i][notij])+2*Q[i][notij]*Q[notij][j])/denum1+
							(Q[i][notij]*Q[notij][j]+Q[i][j]*(Lambda-Q[i][notij]-Q[i][j]-Q[j][i]-Q[j][notij]))/Delta+
							-exp3*(Q[i][j]*(Epsilon-Lambda+2*Q[i][j]+2*Q[j][i]+2*Q[j][notij]+2*Q[i][notij])-2*Q[i][notij]*Q[notij][j])/denum3;

				}
			}

			return DoubleFactory2D.dense.make(returnValue);
		} 
		return null;
	}

	@Test
	public void test() {
		// TODO: This...
		MatrixExponentiator matrixExponentiator1 = new AnalyticMatrixExp3(new double[][]{{-3,1,2},{0.5,-4,3.5},{4,6,-10}});
		MatrixExponentiator matrixExponentiator2 = new Matlab7MatrixExp(new double[][]{{-3,1,2},{0.5,-4,3.5},{4,6,-10}});
		double[][] res1=matrixExponentiator1.expm(0.1).toArray();
		double[][] res2=matrixExponentiator2.expm(0.1).toArray();

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
			res1=matrixExponentiator1.expm(rep/10000).toArray();			
			if (Math.random()<0.00000001) {
				System.out.println(res1[0][0]);
			}
		}
		long time1= System.currentTimeMillis()-startTime1;

		long startTime2= System.currentTimeMillis();
		for (int rep=0;rep<1000000;rep++) {
			res2=matrixExponentiator2.expm(rep/10000).toArray();			
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
