package seasmig;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.util.Vector;

import org.junit.Test;

import seasmig.data.DataForTests;
import seasmig.treelikelihood.MatrixExponentiator;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp2;
import seasmig.treelikelihood.matrixexp.AnalyticMatrixExp3;
import seasmig.treelikelihood.matrixexp.EigenDecomposionExp;
import seasmig.treelikelihood.matrixexp.HKY85MatrixExp;
import seasmig.treelikelihood.matrixexp.JC69MatrixExp;
import seasmig.treelikelihood.matrixexp.JamaMolerMatrixExp;
import seasmig.treelikelihood.matrixexp.JblasMatrixExp;
import seasmig.treelikelihood.matrixexp.Matlab7MatrixExp;
import seasmig.treelikelihood.matrixexp.TaylorMatrixExp;
import seasmig.util.Util;

public class TestMatrixExponentiation {
	// TODO: Organize this....

	final int numScaleSteps = 100;
	final double maxMatrixScale = 10.0;
	final double minMatrixScale = 0.001;
	final double minTime = 0.000;
	final double maxTime = 2;
	final double numTimeSteps = 30.0;
	final double tol = 0.001;
	final int[] testDimensions = {2,3,4,5,6,7,8,9,10,20,40};

	@Test
	public void testJC69() {
		MatrixExponentiator matrixExponentiator1 = new JC69MatrixExp(0.1);
		MatrixExponentiator matrixExponentiator2 = new Matlab7MatrixExp(
				new double[][]{
						{-0.1*3.0/4.0,0.1/4.0,0.1/4.0,0.1/4.0},
						{0.1/4.0,-0.1*3.0/4.0,0.1/4.0,0.1/4.0},
						{0.1/4.0,0.1/4.0,-0.1*3.0/4.0,0.1/4.0},
						{0.1/4.0,0.1/4.0,0.1/4.0,-0.1*3.0/4.0}});
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
		
		for (int i=0; i<4; i++) 
			assertArrayEquals(res1[i],res2[i],0.0000001);
		
	

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

	@Test
	public void testHKY85() {
		double mu = 0.2;
		double kappa = 0.5;
		double piC = 0.3;
		double piA = 0.2;
		double piG = 0.26;
		double piT = 1.0 - piC - piA - piG;
		
		double[][] Q = {
				{0,mu*kappa*piC, mu*piA, mu*piG},
				{mu*kappa*piT, 0, mu*piA, mu*piG},
				{mu*piT, mu*piC, 0, mu*kappa*piG},
				{mu*piT, mu*piC, mu*kappa*piA, 0}};		
		
		for (int i=0;i<4;i++) {
			double rowsum=0;
			for (int j=0;j<4;j++) {
				rowsum=rowsum+Q[i][j];
			}
			Q[i][i]=-rowsum;
		}
		
		MatrixExponentiator matrixExponentiator1 = new HKY85MatrixExp(mu,kappa,piC,piA,piG);
		MatrixExponentiator matrixExponentiator2 = new Matlab7MatrixExp(Q);

		double[][] res1=matrixExponentiator1.expm(0.5).toArray();
		double[][] res2=matrixExponentiator2.expm(0.5).toArray();
	
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
		
		for (int i=0; i<4; i++) 
			for (int j=0; j<4; j++) {
				System.out.println(res1[i][j]);
				System.out.println(res2[i][j]);
				assertEquals(res1[i][j],res2[i][j], 0.000001);
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

	@Test
	public void testMatrixExponentiation() {	
		int dotIter=0;

		System.out.println("Testing matrix exponentiation tol="+tol);
		System.out.println("minMatrixScale: "+minMatrixScale+" maxMatrixScale: "+maxMatrixScale+" steps: "+numScaleSteps);
		System.out.println("minTime: "+minTime+" maxTime: "+maxTime+" numTimesteps: "+numTimeSteps);
		System.out.println("Comparing expm(Q*t)");
		for (int numLocations : testDimensions) {	
			System.out.println("\nnumDimensions: "+numLocations);
			for (int scaleIter=0;scaleIter<numScaleSteps;scaleIter++) {				
				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) scaleIter*(maxMatrixScale-minMatrixScale)/(double)numScaleSteps+minMatrixScale);
				Vector<MatrixExponentiator> tests = new Vector<MatrixExponentiator>();
				tests.add(new Matlab7MatrixExp(testMatrix));
				tests.add(new TaylorMatrixExp(testMatrix));
				tests.add(new JamaMolerMatrixExp(testMatrix));
				tests.add(new JblasMatrixExp(testMatrix));
				tests.add(new AnalyticMatrixExp3(testMatrix));
				tests.add(new AnalyticMatrixExp2(testMatrix));
				tests.add(new EigenDecomposionExp(testMatrix));

				for (double t=minTime;t<maxTime;t+=(maxTime-minTime)/numTimeSteps) {
					dotIter++;
					Vector<double[][]> results = new Vector<double[][]>();

					for (MatrixExponentiator expMethod : tests) {
						if (expMethod.checkMethod()) {
							results.add(expMethod.expm(t).toArray());
						}
						else {
							results.add(null);
						}										
					}

					for (int i=0;i<results.size()-1;i++) {
						if (results.get(i)==null) continue;
						for (int j=(i+1);j<results.size();j++) {
							if (results.get(j)==null) continue;
							assertEqualMatrixExp(tests.get(i),tests.get(j),results.get(i),results.get(j),testMatrix,t,tol);							
						}
					}	

					System.out.print(".");
					if (dotIter%200==0) {
						System.out.println();
					}
					if (dotIter%1000==0) {
						System.out.println();
						System.out.println("Example MatrixExp Test:");
						System.out.println("t: "+t+" Q:");
						System.out.println(Util.print(testMatrix));
						for (MatrixExponentiator expMethod : tests) {
							System.out.println("method: "+expMethod.getMethodName());
							if (expMethod.checkMethod()) {
								System.out.println(Util.print(expMethod.expm(t).toArray()));
							}
							else {
								System.out.println("Incompatible Method!!!");
							}
						}
					}
				}				
			}
		}

		System.out.println("\nCompleted matrix exponentiation test");

	}	

	private void assertEqualMatrixExp(MatrixExponentiator matrixExp1, MatrixExponentiator matrixExp2, 
			double[][] mat1, double[][] mat2, double[][] Q, double t, double tol) {
		assertEquals(mat1.length,mat2.length,0);
		assertEquals(mat1[0].length,mat2[0].length,0);
		for (int i=0;i<mat1.length;i++) {
			for (int j=0;j<mat1.length;j++) {
				if (Math.abs(mat1[i][j]-mat2[i][j])>tol) {
					System.err.println("Failed MatrixExp Test:");
					System.err.println("t: "+t+" Q:");
					System.err.println(Util.print(Q));
					System.err.println("method: "+matrixExp1.getMethodName());
					System.err.println(Util.print(mat1));
					System.err.println("method: "+matrixExp2.getMethodName());
					System.err.println(Util.print(mat2));
					assertEquals(true,false);
				}
			}
		}

	}

	//	@Test 
	//	public void testMatrixExpTimingWithCache() {
	//		int dotIter=0;
	//		
	//		System.out.println("Testing matrix exponentiation tol="+tol);
	//		System.out.println("minMatrixScale: "+minMatrixScale+" maxMatrixScale: "+maxMatrixScale+" steps: "+numScaleSteps);
	//		System.out.println("minTime: "+minTime+" maxTime: "+maxTime+" numTimesteps: "+numTimeSteps);
	//		System.out.println("Comparing expm(Q*t)");
	//		for (int numLocations : testDimensions) {	
	//			System.out.println("\nnumDimensions: "+numLocations);
	//			
	//			for (MatrixExponentiator expMethod : tests) {
	//			
	//				for (int scaleIter=0;scaleIter<numScaleSteps;scaleIter++) {				
	//					double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) scaleIter*(maxMatrixScale-minMatrixScale)/(double)numScaleSteps+minMatrixScale);
	//					Vector<MatrixExponentiator> tests = new Vector<MatrixExponentiator>();
	//					tests.add(new Matlab7MatrixExp(testMatrix));
	//					tests.add(new TaylorMatrixExp(testMatrix));
	//					tests.add(new JamaMolerMatrixExp(testMatrix));
	//					tests.add(new JblasMatrixExp(testMatrix));
	//					tests.add(new AnalyticMatrixExp3(testMatrix));
	//					tests.add(new AnalyticMatrixExp2(testMatrix));
	//
	//					for (double t=minTime;t<maxTime;t+=(maxTime-minTime)/numTimeSteps) {
	//						dotIter++;
	//						Vector<double[][]> results = new Vector<double[][]>();
	//
	//				
	//						if (expMethod.checkMethod()) {
	//							results.add(expMethod.expm(t));
	//						}
	//						else {
	//							results.add(null);
	//						}										
	//					}
	//
	//					for (int i=0;i<results.size()-1;i++) {
	//						if (results.get(i)==null) continue;
	//						for (int j=(i+1);j<results.size();j++) {
	//							if (results.get(j)==null) continue;
	//							assertEqualMatrixExp(tests.get(i),tests.get(j),results.get(i),results.get(j),testMatrix,t,tol);							
	//						}
	//					}	
	//
	//					System.out.print(".");
	//					if (dotIter%200==0) {
	//						System.out.println();
	//					}
	//					if (dotIter%1000==0) {
	//						System.out.println();
	//						System.out.println("Example MatrixExp Test:");
	//						System.out.println("t: "+t+" Q:");
	//						System.out.println(Util.print(testMatrix));
	//						for (MatrixExponentiator expMethod : tests) {
	//							System.out.println("method: "+expMethod.getMethodName());
	//							if (expMethod.checkMethod()) {
	//								System.out.println(Util.print(expMethod.expm(t)));
	//							}
	//							else {
	//								System.out.println("Incompatible Method!!!");
	//							}
	//						}
	//					}
	//				}				
	//			}
	//		}
	//
	//		System.out.println("\nCompleted matrix exponentiation timing test");
	//
	//	}	


	//	@Test
	//	public void testMatrixExpTiming() {
	//		
	//		final int[] testDimensions = {2,3,4,5,6,7,8,9,10,20,40};
	//		
	//		for (int numLocations : testDimensions) {	
	//			System.out.println("\nnumDimensions: "+numLocations);
	//			for (int scaleIter=0;scaleIter<numScaleSteps;scaleIter++) {				
	//				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) scaleIter*(maxMatrixScale-minMatrixScale)/(double)numScaleSteps+minMatrixScale);
	//		
	//		Vector<MatrixExponentiator> tests = new Vector<MatrixExponentiator>();
	//		tests.add(new Matlab7MatrixExp(testMatrix));
	//		tests.add(new TaylorMatrixExp(testMatrix));
	//		tests.add(new JamaMolerMatrixExp(testMatrix));
	//		tests.add(new JblasMatrixExp(testMatrix));
	//		tests.add(new AnalyticMatrixExp3(testMatrix));
	//		tests.add(new AnalyticMatrixExp2(testMatrix));
	//		tests.add(new EigenDecomposionExp(testMatrix));
	//		
	//		// TODO: this....
	//		long startTime1= System.currentTimeMillis();		
	//		for (int i=1;i<50;i++) {
	//			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//			MatrixExponentiator test1 = new MolerMatrixExp(testMatrix);
	//			for (int j=0;j<1200;j++) {
	//				double[][] res1=test1.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//				if (Math.random()<0.0000000000001) System.out.println(res1);
	//			}
	//		}
	//		long time1= System.currentTimeMillis()-startTime1;
	//
	//		long startTime2= System.currentTimeMillis();		
	//		for (int i=1;i<50;i++) {
	//			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//			MatrixExponentiator test2 = new Matlab7MatrixExp(testMatrix);
	//			for (int j=0;j<1200;j++) {
	//				double[][] res2=test2.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//				if (Math.random()<0.0000000000001) System.out.println(res2);
	//			}
	//		}
	//		long time2= System.currentTimeMillis()-startTime2;
	//
	//		long startTime3= System.currentTimeMillis();		
	//		for (int i=1;i<50;i++) {
	//			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//			MatrixExponentiator test3 = new TaylorMatrixExp(testMatrix);;
	//			for (int j=0;j<1200;j++) {
	//				double[][] res3=test3.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//				if (Math.random()<0.0000000000001) System.out.println(res3);
	//			}
	//		}
	//		long time3= System.currentTimeMillis()-startTime3;
	//
	//		long startTime4= System.currentTimeMillis();		
	//		for (int i=1;i<50;i++) {
	//			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//			MatrixExponentiator test4 = new JamaMolerMatrixExp(testMatrix);;
	//			for (int j=0;j<1200;j++) {
	//				double[][] res4=test4.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//				if (Math.random()<0.0000000000001) System.out.println(res4);
	//			}
	//		}
	//		long time4= System.currentTimeMillis()-startTime4;
	//
	//		long startTime5= System.currentTimeMillis();		
	//		for (int i=1;i<50;i++) {
	//			double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//			MatrixExponentiator test5 = new JblasMatrixExp(testMatrix);;
	//			for (int j=0;j<1200;j++) {
	//				double[][] res5=test5.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//				if (Math.random()<0.0000000000001) System.out.println(res5);
	//			}
	//		}
	//		long time5= System.currentTimeMillis()-startTime5;
	//
	//		long time6 = 0;
	//		if (numLocations==3) {
	//			long startTime6= System.currentTimeMillis();		
	//			for (int i=1;i<50;i++) {
	//				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//				MatrixExponentiator test6 = new AnalyticMatrixExp3(testMatrix);;
	//				for (int j=0;j<1200;j++) {
	//					double[][] res6=test6.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//					if (Math.random()<0.0000000000001) System.out.println(res6);
	//				}
	//			}
	//			time6= System.currentTimeMillis()-startTime6;
	//		}
	//		if (numLocations==2) {
	//			long startTime6= System.currentTimeMillis();		
	//			for (int i=1;i<50;i++) {
	//				double[][] testMatrix = DataForTests.makeRandomMigrationMatrix(numLocations,(double) i/100.0);
	//				MatrixExponentiator test6 = new AnalyticMatrixExp2(testMatrix);;
	//				for (int j=0;j<1200;j++) {
	//					double[][] res6=test6.expm(cern.jet.random.Uniform.staticNextDoubleFromTo(0, 5));
	//					if (Math.random()<0.0000000000001) System.out.println(res6);
	//				}
	//			}
	//			time6= System.currentTimeMillis()-startTime6;
	//		}
	//
	//
	//		System.out.println("\nMolerMatrixExp: "+time1+"ms");
	//		System.out.println("Matlab7MatrixExp: "+time2+"ms");
	//		System.out.println("TaylorMatrixExp: "+time3+"ms");
	//		System.out.println("JamaMolerMatrixExp: "+time4+"ms");
	//		System.out.println("JblasMatrixExp: "+time5+"ms");
	//		if (numLocations==3) System.out.println("AnalyticMatrixExp3: "+time6+"ms");
	//		if (numLocations==2) System.out.println("AnalyticMatrixExp2: "+time6+"ms");
	//
	//
	//	}

}
