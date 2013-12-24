package seasmig.treelikelihood.trees;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.TransitionModel.Transition;
import seasmig.treelikelihood.transitionmodels.ConstantTransitionBaseModel;
import seasmig.treelikelihood.transitionmodels.TwoSeasonMigrationBaseModel;


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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null, TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(null, TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null, 1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new TwoSeasonMigrationBaseModel(
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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));

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
		TransitionModel equalModel = new ConstantTransitionBaseModel(singleMatrix);		
		TransitionModel twoSeasonModel = new TwoSeasonMigrationBaseModel(matrixSeas1 ,matrixSeas2 , seasonalPhase, seasonalPhase+0.5);
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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.2,0.3,0.5},
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


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
				{ 0.333333,-1,0.333333,0.333333},
				{ 0.333333,0.333333,-1,0.333333},
				{ 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());

		boolean asrOK = false;
		int countParsimonious = 0;
		for (int repeat=1; repeat<1000; repeat++) {		

			root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
			root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
			root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
			root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
			root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
			equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
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
		// TODO:test...
		/* 
		   --
		  /  \
		 0   --
		    /  \
		   1    0


		 */


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,1,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,2,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,3,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
				{ 0.333333,-1,0.333333,0.333333},
				{ 0.333333,0.333333,-1,0.333333},
				{ 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		System.out.println(locTree.newickAncestralStateReconstruction());	
		System.out.println(locTree.newickStochasticMapping(10000));
		assertEquals(false, true);
	}

	@Test
	public void testDescendantsStats() throws Exception {
		// TODO:test...
		/* 
		   -- 
		  /  \
		 0   --
		    /  \
		   1    0


		 */		
		/* 
		   -- (1)
		  /  \
		 0   -- (1)
		    /  \
		   1    0


		 */		


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(null,0,1,1.0,root));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,2,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,3,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1.0/5.0,0.333333/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,-1.0/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,-1.0/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,0.333333/5,-1.0/5.0}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		System.out.println(locTree.newickAncestralStateReconstruction());	
		System.out.println(locTree.newickStochasticMapping(1000));
		System.out.println(locTree.smDescendants());
		assertEquals(false, true);
	}

	@Test
	public void testNextEvent() throws Exception {

		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1.0/5.0,0.333333/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,-1.0/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,-1.0/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,0.333333/5,-1.0/5.0}});
		
		double[] a = new double[4];
		double[] interval = new double[4];
		Transition[] transitions = new Transition[4];
		for (int i=0;i<5000000;i++) {
			transitions[0]=equalModel.nextEvent(0.1, 0);
			transitions[1]=equalModel.nextEvent(0.2, 0);
			transitions[2]=equalModel.nextEvent(0.3, 0);
			transitions[3]=equalModel.nextEvent(0.4, 0);
			a[transitions[0].toTrait]+=1.0/20000000.0;
			a[transitions[1].toTrait]+=1.0/20000000.0;
			a[transitions[2].toTrait]+=1.0/20000000.0;
			a[transitions[3].toTrait]+=1.0/20000000.0;			
			interval[0]+=(transitions[0].time-0.1)/5000000;
			interval[1]+=(transitions[1].time-0.2)/5000000;
			interval[2]+=(transitions[2].time-0.3)/5000000;
			interval[3]+=(transitions[3].time-0.4)/5000000;
		}
		System.out.println(a[0]+","+a[1]+","+a[2]+","+a[3]);		
		assertEquals(0,a[0], 0.00001);
		assertEquals(0.33333333333333,a[1], 0.01 );
		assertEquals(0.33333333333333,a[2], 0.01 );
		assertEquals(0.33333333333333,a[3], 0.01 );
		assertEquals(5,interval[0], 0.01);
		assertEquals(5,interval[1], 0.01 );
		assertEquals(5,interval[2], 0.01 );
		assertEquals(5,interval[3], 0.01 );
		
		// TODO: better test
	}
	
	@Test
	public void testLikelihoodWithSeq() throws Exception {
		/* 
		   --
		  /   \
		 AG0   --
		      /  \
		    CT1   AG0
		 
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
         2*(-2.81603+Log[1/4])=-4.2023210731*3
		 */
		
		double expectedResult = -4.2023210731*3;
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(new Sequence(2), TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsNode(new Sequence("ABCD123","AG"),0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsNode(new Sequence(2),TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","CT"), 1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","AG"),0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TransitionModel[] codonModel = new TransitionModel[]{equalModel,equalModel,equalModel};
		
		TreeWithLocations tree = new TreeWithLocations(root, 4, 2);
		tree.setMigrationModel(equalModel);
		tree.setCodonModel(codonModel);
		
		double logLikelihood = tree.logLikelihood();
		System.out.println(tree.newickProbs());
		
		assertEquals(expectedResult, logLikelihood,0.01);			
	}
}
