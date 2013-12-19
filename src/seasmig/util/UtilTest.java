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
  public void testDirichlet() {
	  double[] result =Util.toDirichletNonDegenerate(new double[]{0.1,0.3,0.6});
	  assertArrayEquals(new double[]{0.1,0.2,0.4},result,0.00001);
	  System.out.println(result[0]);
	  System.out.println(result[1]);
	  System.out.println(result[2]);
	  System.out.println();
	  result =Util.toDirichletNonDegenerate(new double[]{0.6,0.1,0.3});
	  assertArrayEquals(new double[]{0.4,0.1,0.2},result,0.00001);
	  System.out.println(result[0]);
	  System.out.println(result[1]);
	  System.out.println(result[2]);
	  System.out.println();
	  result =Util.toDirichletNonDegenerate(new double[]{0.6,0.2,0.3});
	  assertArrayEquals(new double[]{0.3,0.1,0.2},result,0.00001);
	  System.out.println(result[0]);
	  System.out.println(result[1]);
	  System.out.println(result[2]);
	  System.out.println();
	  result =Util.toDirichletNonDegenerate(new double[]{0.7,0.3,0.1});
	  assertArrayEquals(new double[]{0.4,0.3,0.2},result,0.00001);
	  System.out.println(result[0]);
	  System.out.println(result[1]);
	  System.out.println(result[2]);
	  System.out.println();
	  result =Util.toDirichletNonDegenerate(new double[]{0.1,0.3,0.7});
	  assertArrayEquals(new double[]{0.2,0.3,0.4},result,0.00001);
	  System.out.println(result[0]);
	  System.out.println(result[1]);
	  System.out.println(result[2]);
	  System.out.println();
  }
  
  // TODO: test logSumExp...
  
}
