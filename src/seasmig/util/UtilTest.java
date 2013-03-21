package seasmig.util;

import static org.junit.Assert.*;

import org.junit.*;

public class UtilTest {

  @Before
  public void setUp() throws Exception {
  }

  @After
  public void tearDown() throws Exception {
  }

  @Test
  public void testCalcQMatrix() {
	  double[][] R = new double[][]{{0,1,1.33,1},{1,0,1,1.33},{1.33,1,0,1},{1,1.33,1,0}};
	  double[] b = new double[]{0.1,0.4,0.2,0.3};
	  double[][] Q = Util.calcQMatrix(R, b);
	  double[][] expectedQ = {{-1.218,0.504,0.336,0.378},{0.126, -0.882, 0.252, 0.504},{0.168, 0.504,-1.05, 0.378},{0.126, 0.672, 0.252,-1.05}};
			 
	  for (int i=0;i<Q.length;i++)
	    assertArrayEquals(expectedQ[i],Q[i],0.05);
	  
	  // double b half Q for the same R...
	  R = new double[][]{{0,1,1.33,1},{1,0,1,1.33},{1.33,1,0,1},{1,1.33,1,0}};
	  b = new double[]{0.2,0.8,0.4,0.6};
	  Q = Util.calcQMatrix(R, b);
	  expectedQ = new double[][]{{-1.218/2.0,0.504/2.0,0.336/2.0,0.378/2.0},{0.126/2.0, -0.882/2.0, 0.252/2.0, 0.504/2.0},{0.168/2.0, 0.504/2.0,-1.05/2.0, 0.378/2.0},{0.126/2.0, 0.672/2.0, 0.252/2.0,-1.05/2.0}};
			 
	  for (int i=0;i<Q.length;i++)
	    assertArrayEquals(expectedQ[i],Q[i],0.05);
  }
  
  // TODO: test logSumExp...
  
}
