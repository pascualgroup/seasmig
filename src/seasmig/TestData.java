package seasmig;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import treelikelihood.*;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.math.Random;

public class TestData implements Data {
	// TODO: This...
	public Vector<LikelihoodTree> trees = new Vector<LikelihoodTree>();
	Config config = null;
	enum TestType {NORMAL, TEST_USING_GENERATED_TREES, TEST_USING_INPUT_TREES, TEST_MODEL_DEGENERACY};

	// TEST MODELS
	MigrationBaseModel createModel = null;
	List<MigrationBaseModel> testModels = null;
	private SimpleRootedTree numTestTips;
	private int numLocations;
	private int disturbanceScale;

	public TestData(Config config_, TestType testType, int numTestRepeats, int numTestLocations, int numTestTrees) throws IOException, ImportException 	{

		config = config_;		

		switch (testType) {
		case NORMAL:
			// Load trees
			System.out.print("Loading trees... ");			
			File treeFile = new File(config.treeFilename);
			FileReader reader = new FileReader(treeFile);
			NexusImporter nexusImporter = new NexusImporter(reader);
			List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();
			System.out.println("loaded "+nexsusTrees.size()+" trees");

			System.out.print("Keeping tail... ");		
			List<jebl.evolution.trees.Tree> nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
				nexsusTreeTail.add(nexsusTrees.get(i));
			}
			System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");			

			// TODO: add states....
			// Convert trees to internal tree representation
			if (config.locationFilename!=null) {
				System.out.print("Loading traits... ");
				AttributeLoader attributeLoader= new SimpleAttributeLoader(config.locationFilename, config.stateFilename);	
				HashMap<String,Object> attributes = attributeLoader.getAttributes();
				HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributes.get("locations");
				HashMap<String,Double> stateMap = (HashMap<String,Double>) attributes.get("states");
				int numLocations = (Integer) attributes.get("numLocations");
				System.out.println("loaded "+locationMap.size()+" taxon traits");

				System.out.print("Reparsing trees... ");
				if (stateMap==null) {
					for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
						trees.add(new TreeWithLocations((SimpleRootedTree) tree,locationMap,numLocations));
					}
				}
				else {
					// TODO: this...
				}
				System.out.println(" reparsed "+trees.size()+" trees");
			}
			else {
				// TODO: add load states from trees...
				numLocations=config.numLocations;
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new TreeWithLocations((SimpleRootedTree) tree,config.locationAttributeNameInTree, numLocations));
				}		
				System.out.println(" reparsed "+trees.size()+" trees");
			}	
			break;

		case TEST_USING_GENERATED_TREES:
			
			numLocations=numTestLocations;
			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Generating test trees... ");

			// Generate test data and trees

			// For constant model...
			double[][] Q = makeRandomMigrationMatrix(numTestLocations,2); 

			// For two seasonal model...
			double[][] QW = myMatrixCopy(Q);
			double[][] QS = makeRandomMigrationMatrix(numTestLocations,3); 

			// For sinusoidal model...
			double[][] rates = myMatrixCopy(Q);
			double[][] amps = makeRandomMigrationMatrix(numTestLocations,1);
			double[][] phases = makeRandomMigrationMatrix(numTestLocations,1);

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: case TWO_CONSTANT_SEASONS_FIXED_PHASE:
				double phase = 0.3; double length = 0.5;
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL:
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			for (int i=0;i<numTestTrees;i++) {
				TreeWithLocations testTree = new TreeWithLocations(createModel,numTestTips);
				testTree.removeInternalLocations();
				trees.add(testTree);
			}

			System.out.println(" generated "+trees.size()+" trees");
			System.out.print("Generating test models... ");

			testModels = new ArrayList<MigrationBaseModel>();

			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new ConstantMigrationBaseModel(disturbMigrationMatrix(Q, disturbanceScale*i,99999)));
			}
			for (int i=0; i<numTestRepeats; i++) {
				double phase = Math.max(0,Math.min(1,0.3+i/3*(Random.nextDouble()-0.5))); double length = 0.5;
				testModels.add(new TwoSeasonMigrationBaseModel(disturbMigrationMatrix(QW,disturbanceScale*i/3,99999),disturbMigrationMatrix(QS,disturbanceScale*i/3,99999),phase, phase+length));
			}
			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new SinusoidialSeasonalMigrationBaseModel(disturbMigrationMatrix(rates,disturbanceScale*i/3,999999),disturbMigrationMatrix(amps,disturbanceScale*i/3,1),disturbMigrationMatrix(phases,disturbanceScale*i/3,1)));
			}
			System.out.println(" generated "+testModels.size()+" test models");
			break;

		case TEST_USING_INPUT_TREES:{
			
			numLocations=numTestLocations;
			// TODO: get number of locations from trees...
			
			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Generating test trees based on input tree topology ... ");

			// Generate test data and trees

			// For constant model...
			Q = makeRandomMigrationMatrix(numTestLocations,2); 

			// For two seasonal model...
			QW = myMatrixCopy(Q);
			QS = makeRandomMigrationMatrix(numTestLocations,3); 

			// For sinusoidal model...
			rates = myMatrixCopy(Q);
			amps = makeRandomMigrationMatrix(numTestLocations,1);
			phases = makeRandomMigrationMatrix(numTestLocations,1);

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: case TWO_CONSTANT_SEASONS_FIXED_PHASE:
				double phase = 0.3; double length = 0.5;
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL: 
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			System.out.print("Loading trees... ");			
			treeFile = new File(config.treeFilename);
			reader = new FileReader(treeFile);
			nexusImporter = new NexusImporter(reader);
			nexsusTrees = nexusImporter.importTrees();
			System.out.println("loaded "+nexsusTrees.size()+" trees");

			System.out.print("Keeping tail... ");		
			nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
				nexsusTreeTail.add(nexsusTrees.get(i));
			}
			System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");			
			for (int i=0; i<numTestTrees;i++) {				
				TreeWithLocations testTree = new TreeWithLocations(createModel,(SimpleRootedTree) nexsusTreeTail.get(Random.nextInt(nexsusTreeTail.size())));
				testTree.fillRandomTraits();
				testTree.removeInternalLocations();
				trees.add(testTree);				
			}

			System.out.println("Generated "+trees.size()+" model based random tip annotations with input tree topology");
			System.out.print("Generating test models... ");

			testModels = new ArrayList<MigrationBaseModel>();

			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new ConstantMigrationBaseModel(disturbMigrationMatrix(Q, disturbanceScale*i,99999)));
			}
			for (int i=0; i<numTestRepeats; i++) {
				double phase = Math.max(0,Math.min(1,0.3+i/3*(Random.nextDouble()-0.5))); double length = 0.5;
				testModels.add(new TwoSeasonMigrationBaseModel(disturbMigrationMatrix(QW,disturbanceScale*i/3,99999),disturbMigrationMatrix(QS,disturbanceScale*i/3,99999),phase, phase+length));
			}
			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new SinusoidialSeasonalMigrationBaseModel(disturbMigrationMatrix(rates,disturbanceScale*i/3,999999),disturbMigrationMatrix(amps,disturbanceScale*i/3,1),disturbMigrationMatrix(phases,disturbanceScale*i/3,1)));
			}
			System.out.println(" generated "+testModels.size()+" test models");

		}
		break;

		case TEST_MODEL_DEGENERACY:{
			
			numLocations=numTestLocations;
			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Building degenerate test models... ");

			// Generate test data and trees

			// For constant model...
			Q = makeRandomMigrationMatrix(numLocations,2); 

			// For two seasonal model...
			QW = myMatrixCopy(Q);
			QS = myMatrixCopy(Q);

			// For sinusoidal model...
			rates = myMatrixCopy(Q);
			amps = makeRandomMigrationMatrix(numLocations,0);
			phases = makeRandomMigrationMatrix(numLocations,1);

			// For two constant seasons model...
			double phase = 0.3; double length = 0.5;

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: case TWO_CONSTANT_SEASONS_FIXED_PHASE:
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL:
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			System.out.print(" done!\nLoading trees... ");			
			treeFile = new File(config.treeFilename);
			reader = new FileReader(treeFile);
			nexusImporter = new NexusImporter(reader);
			nexsusTrees = nexusImporter.importTrees();
			System.out.println("loaded "+nexsusTrees.size()+" trees");

			System.out.print("Keeping tail... ");		
			nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
				nexsusTreeTail.add(nexsusTrees.get(i));
			}
			System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");			

			// Convert trees to internal tree representation
			if (config.locationFilename!=null) {
				System.out.print("Loading traits... ");
				AttributeLoader attributeLoader= new SimpleAttributeLoader(config.locationFilename, config.stateFilename);
				// TODO: think about this...
				HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributeLoader.getAttributes().get("locations");
				HashMap<String,Double> stateMap = (HashMap<String,Double>) attributeLoader.getAttributes().get("states");
				System.out.println("loaded "+locationMap.size()+" taxon traits");

				System.out.print("Reparsing trees... ");
				if (stateMap==null) {
					for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
						trees.add(new TreeWithLocations((SimpleRootedTree) tree,locationMap,numLocations));
					}
				}
				else {
					// TODO: this...
				}
				System.out.println(" reparsed "+trees.size()+" trees");
			}
			else {
				// TODO: add load states from trees...
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new TreeWithLocations((SimpleRootedTree) tree,config.locationAttributeNameInTree, numLocations));
				}		
				System.out.println(" reparsed "+trees.size()+" trees");
			}				

			testModels = new ArrayList<MigrationBaseModel>();

			testModels.add(new ConstantMigrationBaseModel(Q));
			testModels.add(new TwoSeasonMigrationBaseModel(QW,QS,phase, phase+length));
			testModels.add(new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases));

			System.out.println(" generated "+testModels.size()+" test models");

		}
		break;

		}

	}

	private double[][] myMatrixCopy(double[][] q) {
		double [][] returnValue = new double[q.length][];
		for(int i = 0; i < q.length; i++) {
			returnValue[i] = q[i].clone();
		}				
		return returnValue;
	}

	static double[][] makeRandomMigrationMatrix(int size, double scale) {
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

	private double[][] disturbMigrationMatrix(double[][] migrationMatrix, double disturbanceMagnitude, double max) {
		// For test purposes...
		double[][] returnValue = myMatrixCopy(migrationMatrix);
		for (int i=0;i<migrationMatrix.length;i++) {
			double rowSum=0;
			for (int j=0;j<migrationMatrix.length;j++) {
				if (i!=j) {
					do {
						returnValue[i][j]+=(Math.random()-0.5)*disturbanceMagnitude;
					} while ((returnValue[i][j]<0) || (returnValue[i][j]>max));
					rowSum+=returnValue[i][j];
				}
			}
			returnValue[i][i]=-rowSum;
		}
		return returnValue;
	}

	@Override
	public List<LikelihoodTree> getTrees() {
		return trees;
	}

	@Override
	public int getNumLocations() {
		return numLocations;
	}






}
