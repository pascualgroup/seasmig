package seasmig.treelikelihood.trees;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.models.ConstantMigrationBaseModel;


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
		
		TreeWithLocationsAndNucleotidesNode root = new TreeWithLocationsAndNucleotidesNode(new Sequence(2), TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,0.0,null);		
		root.addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("AG"),0,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,1.0,root));
		root.addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence(2),TreeWithLocationsAndNucleotides.UNKNOWN_LOCATION,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,1.0,null));
		root.children.get(1).addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("CT"), 1,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,2.0,null));
		root.children.get(1).addChild(new TreeWithLocationsAndNucleotidesNode(new Sequence("AG"),0,TreeWithLocationsAndNucleotides.UNKNOWN_TAXA,2.0,null));
		TransitionModel equalModel = new ConstantMigrationBaseModel(new double[][]{{-1,0.333333,0.333333,0.333333},
																 { 0.333333,-1,0.333333,0.333333},
																 { 0.333333,0.333333,-1,0.333333},
																 { 0.333333,0.333333,0.333333,-1}});
		TransitionModel[] codonModel = new TransitionModel[]{equalModel,equalModel,equalModel};
		
		TreeWithLocationsAndNucleotides tree = new TreeWithLocationsAndNucleotides(root, 4, 2);
		tree.setMigrationModel(equalModel);
		tree.setCodonModel(codonModel);
		
		double logLikelihood = tree.logLikelihood();
		System.out.println(tree.newickProbs());
		
		assertEquals(expectedResult, logLikelihood,0.01);			
	}
}
