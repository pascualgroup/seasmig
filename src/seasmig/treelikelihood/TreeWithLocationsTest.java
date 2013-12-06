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
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		
		assertEquals(expectedResult, locTree.logLikelihood(),0.01);			
	}
	
	@Test
	public void testLikelihoodTwoSeasons() throws Exception {
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
		 2,2 --> 3, 1 ;;; 		 2,0 --> 3, 1 ;;; 		 2,3 --> 4 ;;; 		 2,1 --> 3, 1
		 ---
         0,2 --> 3, 1 ;;;        0,0 --> 1, 3 ;;;        0,3 --> 3, 1 ;;;    0,1 --> 2, 2
         ---
         3,2 --> 4 ;;;           3,0 --> 3, 1 ;;;        3,3 --> 3, 1 ;;;    3,1 --> 3, 1
         ---
         1,2 --> 4 ;;;           1,0 --> 3, 1 ;;;        1,3 --> 4 ;;;       1,1 --> 2, 2
         
         Total: Log[9*(c^3*nc^1)+4*(c^4)+1*(c^1*nc^3)+2*(c^2*nc^2)]=-2.81603
         plus Log[1/4] for root freq =  
         2*(-2.81603+Log[1/4])=-4.2023210731
		 */
		
		double expectedResult = -4.2023210731;
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new TwoSeasonMigrationBaseModel(
				new double[][]{{-1*2.0/3.0,0.333333*2.0/3.0,0.333333*2.0/3.0,0.333333*2.0/3.0},
						 	   { 0.333333*2.0/3.0,-1*2.0/3.0,0.333333*2.0/3.0,0.333333*2.0/3.0},
						       { 0.333333*2.0/3.0,0.333333*2.0/3.0,-1*2.0/3.0,0.333333*2.0/3.0},
						       { 0.333333*2.0/3.0,0.333333*2.0/3.0,0.333333*2.0/3.0,-1*2.0/3.0}},
				new double[][]{{-1*4.0/3.0,0.333333*4.0/3.0,0.333333*4.0/3.0,0.333333*4.0/3.0},
				  			   { 0.333333*4.0/3.0,-1*4.0/3.0,0.333333*4.0/3.0,0.333333*4.0/3.0},
							   { 0.333333*4.0/3.0,0.333333*4.0/3.0,-1*4.0/3.0,0.333333*4.0/3.0},
							   { 0.333333*4.0/3.0,0.333333*4.0/3.0,0.333333*4.0/3.0,-1*4.0/3.0}},
				  	 		   0.3,0.8);
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		
		assertEquals(expectedResult, locTree.logLikelihood(),0.01);			
	}
	
	@Test
	public void testLikelihoodTwoSeasons2() throws Exception {		
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));

		double[][] singleMatrix = makeRandomMigrationMatrix(cern.jet.random.Uniform.staticNextIntFromTo(2, 10),				
				                                            10*cern.jet.random.Uniform.staticNextDouble());
		double[][] matrixSeas1 = singleMatrix.clone();
		double[][] matrixSeas2 = singleMatrix.clone();
		for (int i=0;i<matrixSeas1.length;i++) {
			for (int j=0;j<matrixSeas1[0].length;j++) {
				matrixSeas1[i][j]=matrixSeas1[i][j]*2.0/3.0;
				matrixSeas2[i][j]=matrixSeas2[i][j]*4.0/3.0;
			}
		}
		double seasonalPhase = cern.jet.random.Uniform.staticNextDouble()*0.5;
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(singleMatrix);		
		MigrationBaseModel twoSeasonModel = new TwoSeasonMigrationBaseModel(matrixSeas1 ,matrixSeas2 , seasonalPhase, seasonalPhase+0.5);
		TreeWithLocations locTreeEqual = new TreeWithLocations(root,equalModel);
		TreeWithLocations locTreeTwoSeasons = new TreeWithLocations(root,twoSeasonModel);
		
		assertEquals(locTreeEqual.logLikelihood(), locTreeTwoSeasons.logLikelihood(),1E-10);			
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
		MigrationBaseModel equalModel = 
				new 
				ConstantMigrationBaseModel(makeRandomMigrationMatrix(cern.jet.random.Uniform.staticNextIntFromTo(2, 10), 
				10*cern.jet.random.Uniform.staticNextDouble()));
		
		// Tree Likelihood Method 1
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
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
		System.out.println(locTree.newickProbs());
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
		
		 */
		
		double expectedResult = -3.5715396865307345; // This wasn't calculated manually but is the output of the run
		                                            // i.e. good for regression testing
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.2,0.3,0.5},
																 { 0.333333,-1,0.333333,0.333333},
																 { 1.333333,0.333333,-4,2.333333},
																 { 0.333333,0.666666,1.0,-2}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);
		
		assertEquals(expectedResult, locTree.logLikelihood(),0.01);
			
	}
	
	@Test
	public void testASR() throws Exception {
		/* 
		   --
		  /  \
		 0   --
		    /  \
		   1    0
		 		
         
		 */
		
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		
		boolean asrOK = false;
		int countParsimonious = 0;
		for (int repeat=1; repeat<1000; repeat++) {		

			root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
			root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
			root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
			root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
			root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
			equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																	 { 0.333333,-1,0.333333,0.333333},
																	 { 0.333333,0.333333,-1,0.333333},
																	 { 0.333333,0.333333,0.333333,-1}});
			locTree = new TreeWithLocations(root,equalModel);		
			locTree.logLikelihood();
			String asrNewick = locTree.newickAncestralStateReconstruction();			
			System.out.println(asrNewick);
			for (int i=0;i<4;i++) {
				for (int j=0; j<4; j++) {									     
					if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states="+Integer.toString(i)+"]:1.000)[&states="+Integer.toString(j)+"]:0.000\n")) {
						asrOK = true;
					}					
				}
			}			
			if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states=0]:1.000)[&states=0]:0.000\n")) {
				countParsimonious+=1;
			}
			assertEquals(asrOK, true);
		}
		assertEquals(countParsimonious > 190 ,true);
	}
	
	
	@Test
	public void testSM() throws Exception {
		// TODO: better test...
		/* 
		   --
		  /  \
		 0   --
		    /  \
		   1    0
		 		
         
		 */
		
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(0,1,1.0,root));
		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(1,2,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(0,3,2.0,null));
		MigrationBaseModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		System.out.println(locTree.newickAncestralStateReconstruction());	
		System.out.println(locTree.newickStochasticMapping());
		assertEquals(false, true);
	}
}
