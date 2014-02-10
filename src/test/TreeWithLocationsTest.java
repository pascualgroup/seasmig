package test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.TransitionModel.Transition;
import seasmig.treelikelihood.transitionmodels.ConstantTransitionBaseModel;
import seasmig.treelikelihood.transitionmodels.TwoSeasonMigrationBaseModel;
import seasmig.treelikelihood.trees.Sequence;
import seasmig.treelikelihood.trees.TreeWithLocations;
import seasmig.treelikelihood.trees.TreeWithLocationsNode;


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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null, TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null, TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null, 1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
				{ 0.333333,-1,0.333333,0.333333},
				{ 0.333333,0.333333,-1,0.333333},
				{ 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);

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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
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
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);

		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());

		assertEquals(expectedResult, locTree.logLikelihood(),0.01);			
	}

	@Test
	public void testLikelihoodTwoSeasons2() throws Exception {		

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));

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
		TreeWithLocations locTreeEqual = new TreeWithLocations(root,equalModel,null);
		TreeWithLocations locTreeTwoSeasons = new TreeWithLocations(root,twoSeasonModel,null);

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

		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.2,0.3,0.5},
				{ 0.333333,-1,0.333333,0.333333},
				{ 1.333333,0.333333,-4,2.333333},
				{ 0.333333,0.666666,1.0,-2}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);

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


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
				{ 0.333333,-1,0.333333,0.333333},
				{ 0.333333,0.333333,-1,0.333333},
				{ 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());

		boolean asrOK = false;
		int countParsimonious = 0;
		for (int repeat=1; repeat<1000; repeat++) {		

			root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
			root.addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,false));
			root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
			root.children.get(1).addChild(new TreeWithLocationsNode(null,1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
			root.children.get(1).addChild(new TreeWithLocationsNode(null,0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,false));
			equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
					{ 0.333333,-1,0.333333,0.333333},
					{ 0.333333,0.333333,-1,0.333333},
					{ 0.333333,0.333333,0.333333,-1}});
			locTree = new TreeWithLocations(root,equalModel,null);		
			locTree.logLikelihood();
			String asrNewick = locTree.newickASR();
			for (int i=0;i<4;i++) {
				for (int j=0; j<4; j++) {	
					System.out.println(asrNewick);
					if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states="+Integer.toString(i)+"]:1.000)[&states="+Integer.toString(j)+"]:0.000")) {
						asrOK = true;
					}					
				}
			}			
			if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states=0]:1.000)[&states=0]:0.000")) {
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
		   --0
		  /  \
		 0   --1
		    /  \
		   1    0


		 */


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,1,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,2,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,3,2.0,null,false));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
				{ 0.333333,-1,0.333333,0.333333},
				{ 0.333333,0.333333,-1,0.333333},
				{ 0.333333,0.333333,0.333333,-1}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		System.out.println(locTree.newickASR());	
		System.out.println(locTree.newickSM(10000));
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


		TreeWithLocationsNode root = new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,false);		
		root.addChild(new TreeWithLocationsNode(null,0,1,1.0,root,false));
		root.addChild(new TreeWithLocationsNode(null,TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,1,2,2.0,null,false));
		root.children.get(1).addChild(new TreeWithLocationsNode(null,0,3,2.0,null,false));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1.0/5.0,0.333333/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,-1.0/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,-1.0/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,0.333333/5,-1.0/5.0}});
		TreeWithLocations locTree = new TreeWithLocations(root,equalModel,null);		
		locTree.logLikelihood();
		System.out.println(locTree.newickProbs());
		System.out.println(locTree.newickASR());	
		System.out.println(locTree.newickSM(1000));
		System.out.println(locTree.smDescendants());
		assertEquals(false, true);
	}

	@Test
	public void testNextEvent() throws Exception {

		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]
				{{-1.0*5,0.333333*5,0.333333*5,0.333333*5},
				{ 0.333333/5.0,-1.0/5.0,0.333333/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,-1.0/5.0,0.333333/5.0},
				{ 0.333333/5.0,0.333333/5.0,0.333333/5,-1.0/5.0}});
		
		double[] a = new double[4];
		double[] interval = new double[4];
		Transition[] transitions = new Transition[4];
		for (int i=0;i<500000;i++) {
			transitions[0]=equalModel.nextEvent(0.1, 0);
			transitions[1]=equalModel.nextEvent(0.2, 0);
			transitions[2]=equalModel.nextEvent(0.3, 0);
			transitions[3]=equalModel.nextEvent(0.4, 0);
			a[transitions[0].toTrait]+=1.0/2000000.0;
			a[transitions[1].toTrait]+=1.0/2000000.0;
			a[transitions[2].toTrait]+=1.0/2000000.0;
			a[transitions[3].toTrait]+=1.0/2000000.0;			
			interval[0]+=(transitions[0].time-0.1)/500000;
			interval[1]+=(transitions[1].time-0.2)/500000;
			interval[2]+=(transitions[2].time-0.3)/500000;
			interval[3]+=(transitions[3].time-0.4)/500000;
		}
		System.out.println(a[0]+","+a[1]+","+a[2]+","+a[3]);		
		assertEquals(0,a[0], 0.00001);
		assertEquals(0.33333333333333,a[1], 0.01 );
		assertEquals(0.33333333333333,a[2], 0.01 );
		assertEquals(0.33333333333333,a[3], 0.01 );
		assertEquals(0.2,interval[0], 0.01);
		assertEquals(0.2,interval[1], 0.01 );
		assertEquals(0.2,interval[2], 0.01 );
		assertEquals(0.2,interval[3], 0.01 );
		
		TransitionModel jkModel = new ConstantTransitionBaseModel(0.1);
		
		a = new double[4];
		interval = new double[4];
		transitions = new Transition[4];
		for (int i=0;i<500000;i++) {
			transitions[0]=jkModel.nextEvent(0.1, 0);
			transitions[1]=jkModel.nextEvent(0.2, 0);
			transitions[2]=jkModel.nextEvent(0.3, 0);
			transitions[3]=jkModel.nextEvent(0.4, 0);
			a[transitions[0].toTrait]+=1.0/2000000.0;
			a[transitions[1].toTrait]+=1.0/2000000.0;
			a[transitions[2].toTrait]+=1.0/2000000.0;
			a[transitions[3].toTrait]+=1.0/2000000.0;			
			interval[0]+=(transitions[0].time-0.1)/500000;
			interval[1]+=(transitions[1].time-0.2)/500000;
			interval[2]+=(transitions[2].time-0.3)/500000;
			interval[3]+=(transitions[3].time-0.4)/500000;
		}
		System.out.println(a[0]+","+a[1]+","+a[2]+","+a[3]);		
		assertEquals(0,a[0], 0.00001);
		assertEquals(0.33333333333333,a[1], 0.01 );
		assertEquals(0.33333333333333,a[2], 0.01 );
		assertEquals(0.33333333333333,a[3], 0.01 );
		assertEquals(13.33,interval[0], 0.1);
		assertEquals(13.33,interval[1], 0.1 );
		assertEquals(13.33,interval[2], 0.1 );
		assertEquals(13.33,interval[3], 0.1 );
		
		TransitionModel hkyModel = new ConstantTransitionBaseModel(0.1,1.0,0.125,0.25,0.25);
		
		a = new double[4];
		interval = new double[4];
		transitions = new Transition[4];
		for (int i=0;i<500000;i++) {
			transitions[0]=hkyModel.nextEvent(0.1, 0);
			transitions[1]=hkyModel.nextEvent(0.2, 1);
			transitions[2]=hkyModel.nextEvent(0.3, 2);
			transitions[3]=hkyModel.nextEvent(0.4, 3);
			a[transitions[0].toTrait]+=1.0/2000000.0;
			a[transitions[1].toTrait]+=1.0/2000000.0;
			a[transitions[2].toTrait]+=1.0/200000.0;
			a[transitions[3].toTrait]+=1.0/20000000.0;			
			interval[0]+=(transitions[0].time-0.1)/500000;
			interval[1]+=(transitions[1].time-0.2)/500000;
			interval[2]+=(transitions[2].time-0.3)/500000;
			interval[3]+=(transitions[3].time-0.4)/500000;
		}
		System.out.println(a[0]+","+a[1]+","+a[2]+","+a[3]);		
//		assertEquals(0,a[0], 0.00001); // TODO:
//		assertEquals(0.33333333333333,a[1], 0.01 ); // TODO:
//		assertEquals(0.33333333333333,a[2], 0.01 ); // TODO:
//		assertEquals(0.33333333333333,a[3], 0.01 ); // TODO:
		assertEquals(13.333,interval[0], 0.1);
		assertEquals(13.333,interval[1], 0.1 );
		assertEquals(13.333,interval[2], 0.1 );
		assertEquals(13.333,interval[3], 0.1 );
		
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
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(new Sequence(2), TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,true);		
		root.addChild(new TreeWithLocationsNode(new Sequence("ABCD123","AG"),0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,true));
		root.addChild(new TreeWithLocationsNode(new Sequence(2),TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,true));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","CT"), 1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,true));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","AG"),0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,true));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TransitionModel[] codonModel = new TransitionModel[]{equalModel,equalModel,equalModel};
		
		TreeWithLocations tree = new TreeWithLocations(root, 4, 2,null);
		tree.setMigrationModel(equalModel);
		tree.setCodonModel(codonModel);
		
		double logLikelihood = tree.logLikelihood();
		System.out.println(tree.newickProbs());
		
		assertEquals(expectedResult, logLikelihood,0.01);			
	}
	
	@Test
	public void testPreorderIter() throws Exception {
		/* 
		   --
		  /   \
		 AG0   --
		      /  \
		    CT1   CC0
		 
		 */
		
		
		TreeWithLocationsNode root = new TreeWithLocationsNode(new Sequence(2), TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null,true);		
		root.addChild(new TreeWithLocationsNode(new Sequence("ABCD123","AG"),0,TreeWithLocations.UNKNOWN_TAXA,1.0,root,true));
		root.addChild(new TreeWithLocationsNode(new Sequence(2),TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null,true));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","CT"), 1,TreeWithLocations.UNKNOWN_TAXA,2.0,null,true));
		root.children.get(1).addChild(new TreeWithLocationsNode(new Sequence("ABCD123","CC"),0,TreeWithLocations.UNKNOWN_TAXA,2.0,null,true));
		TransitionModel equalModel = new ConstantTransitionBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TransitionModel[] codonModel = new TransitionModel[]{equalModel,equalModel,equalModel};
		
		TreeWithLocations tree = new TreeWithLocations(root, 4, 2,null);
		tree.setMigrationModel(equalModel);
		tree.setCodonModel(codonModel);
			
		String expectedResult = "??AG??CTCC";
		String returnValue="";
		for (TreeWithLocationsNode node : tree.eachPreorder()) {
			returnValue+=node.seq;			
		}
		assertEquals(returnValue,expectedResult);
	}
}
