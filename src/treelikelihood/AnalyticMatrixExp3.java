package treelikelihood;

public class AnalyticMatrixExp3 implements MatrixExponentiator {

	double[][] Q = new double[3][3];
	double Lambda;
	double Delta; 
	double Disc;
	
	public AnalyticMatrixExp3(double[][] Q) {
		for (int i=0;i<Q.length;i++) {
			for (int j=0;j<Q.length;j++) {
				this.Q[i][j]=Q[i][j];
			}
		}
		
		Lambda = Q[0][1]+Q[0][2]+Q[1][0]+Q[1][2]+Q[2][0]+Q[2][1];
		
		Delta = Q[0][2]*Q[1][0]+Q[0][1]*Q[1][2]+Q[0][2]*Q[1][2]+
				Q[0][1]*Q[2][0]+Q[1][0]*Q[2][0]+Q[1][2]*Q[2][0]+
				Q[0][1]*Q[2][1]+Q[0][2]*Q[2][1]+Q[1][0]*Q[2][1];
		
		Disc = (Math.sqrt(Delta));
	}
	
	@Override
	public double[][] expm(double t) {
		
		double[][] returnValue = new double[2][2];
		
		double denum=Q[0][1]+Q[1][0];
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
