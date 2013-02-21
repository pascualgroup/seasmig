package treelikelihood;

public class AnalyticMatrixExp2 implements MatrixExponentiator {
	
	double[][] Q = new double[2][2]; 
	double denum;
	
	public AnalyticMatrixExp2(double[][] Q) {
		for (int i=0;i<Q.length;i++) {
			for (int j=0;j<Q.length;j++) {
				this.Q[i][j]=Q[i][j];
			}
		}
		denum=Q[0][1]+Q[1][0];
	}
	
	@Override
	public double[][] expm(double t) {
		
		double[][] returnValue = new double[2][2];
		
		if (denum!=0) {
			double exp = Math.exp(t*(-Q[0][1]-Q[1][0]));			
			returnValue[0][0]=(exp*Q[0][1]+Q[1][0])/denum;
			returnValue[0][1]=(exp-1)*Q[0][1]/denum;
			returnValue[1][0]=(exp-1)*Q[1][0]/denum;
			returnValue[1][1]=(exp*Q[1][0]+Q[0][1])/denum;
			return returnValue;
		} 
		
		return returnValue;		
	}

}
