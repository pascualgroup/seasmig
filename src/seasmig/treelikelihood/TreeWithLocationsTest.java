package seasmig.treelikelihood;

import static org.junit.Assert.assertEquals;

import org.junit.Test;


public class TreeWithLocationsTest {

	@Test
	public void testLikelihood() throws Exception {
		/*
		 
		   --
		  /  \
		 0   --
		    /  \
		   1    0
		 
		 For the JC model 
		 for no-change we have p=1/4 E^(-4 t/3) (3 + E^(4 t/3)), 
		 for a change we have p=1/4 E^(-4 t/3) (-1 + E^(4 t/3))
		 
		 branch length is 1....
		 First site:
		         changes, no changes
		 2,2 --> 3, 1
		 2,0 --> 3, 1
		 2,3 --> 4
		 2,1 --> 3, 1
		 ---
         0,2 --> 3, 1
         0,0 --> 1, 3
         0,3 --> 3, 1
         0,1 --> 2, 2
         ---
         3,2 --> 4
         3,0 --> 3, 1
         3,3 --> 3, 1
         3,1 --> 3, 1
         ---
         1,2 --> 4
         1,0 --> 3, 1
         1,3 --> 4
         1,1 --> 2, 2
         
         Total: Log[9*(c^3*nc^1)+4*(c^4)+1*(c^1*nc^3)+2*(c^2*nc^2)]=-2.81603
         plus Log[1/4] for root freq =  
         2*(-2.81603+Log[1/4])=-4.2023210731
		 */
		
		double expectedResult = -4.2023210731;
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		locTree.logLikelihood();
		System.out.println(locTree.newick());
		
		assertEquals(expectedResult, locTree.logLikelihood(),0.01);
			
	}
	
	public static double[][] makeRandomMigrationMatrix(int size, double scale) {
		// For test purposes...
		double[][] returnValue = new double[size][size];
		for (int i=0;i<size;i++) {
			double rowSum=0;
			for (int j=0;j<size;j++) {
				if (i!=j) {
					returnValue[i][j]=Math.random()*scale;
					rowSum+=returnValue[i][j];
				}
			}
			returnValue[i][i]=-rowSum;
		}
		return returnValue;
	}
	
	@Test
	public void testEqualLikelihoodForTwoDifferentImplementations() throws Exception {
		/*
		 
		   --
		  /  \
		 0   --
		    /  \
		   1    0
		 
	
		 */
		
		// Model
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(makeRandomMigrationMatrix(cern.jet.random.Uniform.staticNextIntFromTo(2, 10), 10*cern.jet.random.Uniform.staticNextDouble()));
		
		// Tree Likelihood Method 1
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,root.children.get(1)));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,root.children.get(1)));
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		// Tree Likelihood Method 2

		int UNKNOWN_LOCATION = equalModel.getNumLocations();
		int NUM_LOCATIONS= equalModel.getNumLocations();
		TreeWithLocationsNode2 root2 = new TreeWithLocationsNode2(UNKNOWN_LOCATION,0.0, NUM_LOCATIONS);		
		root2.children.add(new TreeWithLocationsNode2(0,1.0,NUM_LOCATIONS));
		root2.children.add(new TreeWithLocationsNode2(UNKNOWN_LOCATION,1.0, NUM_LOCATIONS));
		root2.children.get(1).children.add(new TreeWithLocationsNode2(1,2.0, NUM_LOCATIONS));
		root2.children.get(1).children.add(new  TreeWithLocationsNode2(0,2.0, NUM_LOCATIONS));

		TreeWithLocations2 locTree2 = new TreeWithLocations2(equalModel,root2);
		
		
		// Test
		locTree.logLikelihood();
		System.out.println(locTree.newick());
		assertEquals(locTree.logLikelihood(), locTree2.logLikelihood(),0.01);
	
			
	}
	
	@Test
	public void testLikelihoodAsymetric() throws Exception {
		
		// TODO: calculate likelihood for this...
		/*
		 
		   --
		  /  \
		 0   --
		    /  \
		   1    0
		 
		
		 branch length is 1....
		 First site:
		         changes, no changes
		 2,2 --> 3, 1
		 2,0 --> 3, 1
		 2,3 --> 4
		 2,1 --> 3, 1
		 ---
         0,2 --> 3, 1
         0,0 --> 1, 3
         0,3 --> 3, 1
         0,1 --> 2, 2
         ---
         3,2 --> 4
         3,0 --> 3, 1
         3,3 --> 3, 1
         3,1 --> 3, 1
         ---
         1,2 --> 4
         1,0 --> 3, 1
         1,3 --> 4
         1,1 --> 2, 2
         
         Total: ...... Log[9*(c^3*nc^1)+4*(c^4)+1*(c^1*nc^3)+2*(c^2*nc^2)]=.....
         plus ..... Log[1/4] for root freq =  
         2*(-2.81603+Log[1/4])=....
		 */
		
		double expectedResult = -1234567890;
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocationsNode.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		assertEquals(expectedResult, locTree.logLikelihood(),0.01);
			
	}
}
