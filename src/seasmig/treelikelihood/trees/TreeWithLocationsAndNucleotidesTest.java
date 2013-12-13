package seasmig.treelikelihood.trees;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.models.ConstantMigrationBaseModel;
import seasmig.treelikelihood.models.TwoSeasonMigrationBaseModel;


public class TreeWithLocationsAndNucleotidesTest {

	@Test
	public void testLikelihood() throws Exception {
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
		
		TreeWithLocationsAndNucleotidesNode root = new TreeWithLocationsAndNucleotidesNode(new Sequence(1), TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("A"),0,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence(1),TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("C"), 1,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("A"),0,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TransitionModel[] codonModel = new TransitionModel[]{equalModel,equalModel,equalModel};
		
		TreeWithLocationsAndNucleotides tree = new TreeWithLocationsAndNucleotides(root, 4, 1);
		tree.setMigrationModel(equalModel);
		tree.setCodonModel(codonModel);
		
		tree.logLikelihood();
		System.out.println(tree.newickProbs());
		
		assertEquals(expectedResult, tree.logLikelihood(),0.01);			
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
		
		
//		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
//		root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
//		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
//		root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
//		root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
//		TransitionModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
//																 { 0.333333,-1,0.333333,0.333333},
//																 { 0.333333,0.333333,-1,0.333333},
//																 { 0.333333,0.333333,0.333333,-1}});
//		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
//		locTree.logLikelihood();
//		System.out.println(locTree.newickProbs());
//		
//		boolean asrOK = false;
//		int countParsimonious = 0;
//		for (int repeat=1; repeat<1000; repeat++) {		
//
//			root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
//			root.addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,1.0,root));
//			root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
//			root.children.get(1).addChild(new TreeWithLocationsNode(1,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
//			root.children.get(1).addChild(new TreeWithLocationsNode(0,TreeWithLocations.UNKNOWN_TAXA,2.0,null));
//			equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
//																	 { 0.333333,-1,0.333333,0.333333},
//																	 { 0.333333,0.333333,-1,0.333333},
//																	 { 0.333333,0.333333,0.333333,-1}});
//			locTree = new TreeWithLocations(root,equalModel);		
//			locTree.logLikelihood();
//			String asrNewick = locTree.newickAncestralStateReconstruction();			
//			System.out.println(asrNewick);
//			for (int i=0;i<4;i++) {
//				for (int j=0; j<4; j++) {									     
//					if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states="+Integer.toString(i)+"]:1.000)[&states="+Integer.toString(j)+"]:0.000\n")) {
//						asrOK = true;
//					}					
//				}
//			}			
//			if (asrNewick.equals("(-1[&states=0]:1.000,(-1[&states=1]:1.000,-1[&states=0]:1.000)[&states=0]:1.000)[&states=0]:0.000\n")) {
//				countParsimonious+=1;
//			}
//			assertEquals(asrOK, true);
//		}
//		assertEquals(countParsimonious > 190 ,true);
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
		
		
//		TreeWithLocationsNode root = new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,0.0,null);		
//		root.addChild(new TreeWithLocationsNode(0,1,1.0,root));
//		root.addChild(new TreeWithLocationsNode(TreeWithLocations.UNKNOWN_LOCATION,TreeWithLocations.UNKNOWN_TAXA,1.0,null));
//		root.children.get(1).addChild(new TreeWithLocationsNode(1,2,2.0,null));
//		root.children.get(1).addChild(new TreeWithLocationsNode(0,3,2.0,null));
//		TransitionModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
//																 { 0.333333,-1,0.333333,0.333333},
//																 { 0.333333,0.333333,-1,0.333333},
//																 { 0.333333,0.333333,0.333333,-1}});
//		TreeWithLocations locTree = new TreeWithLocations(root,equalModel);		
//		locTree.logLikelihood();
//		System.out.println(locTree.newickProbs());
//		System.out.println(locTree.newickAncestralStateReconstruction());	
//		System.out.println(locTree.newickStochasticMapping());
//		assertEquals(false, true);
	}
}
