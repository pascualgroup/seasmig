package seasmig;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import treelikelihood.*;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;

public class Data
{


	public Vector<Tree> trees = new Vector<Tree>();
	Config config = null;

	// TEST MODELS
	MigrationBaseModel createModel = null;
	MigrationBaseModel testModel = null;

	public Data(Config config_) throws IOException, ImportException 	{

		config = config_;		

		switch (config.runMode) {
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
				// TODO: think about this...
				HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributeLoader.getAttributes().get("locations");
				System.out.println("loaded "+locationMap.size()+" taxon traits");

				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new Tree((SimpleRootedTree) tree,locationMap,config.locationCount));
				}
				System.out.println(" reparsed "+trees.size()+" trees");
			}
			else {
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new Tree((SimpleRootedTree) tree,config.locationCount));
				}		
				System.out.println(" reparsed "+trees.size()+" trees");
			}	
			break;
		case TEST:
			// TODO: add test files instead of hardcoding matrices ....
			System.out.print("Generating test trees... ");

			// Generate test data and trees

			// For constant model...
			double[][] Q = {{-0.9,0.4,0.3,0.2},
					{0.3,-0.6,0.2,0.1},
					{0.2,0.1,-0.3,0.0},
					{0.1,0.0,0.0,-0.1}};

			// For two seasonal model...
			double[][] QW = {{-0.9,0.8,0.1,0.0},
					{0.2,-0.5,0.2,0.1},
					{0.2,0.0,-0.2,0.0},
					{0.3,0.0,0.0,-0.3}};

			double[][] QS = {{-0.9,0.4,0.3,0.2},
					{0.3,-0.6,0.2,0.1},
					{0.2,0.1,-0.3,0.0},
					{0.1,0.0,0.0,-0.1}};

			// For sinusoidal model...
			double[][] rates = {{0,0.4,0.3,0.2},
					{0.3,0,0.2,0.1},
					{0.2,0.1,0,0.0},
					{0.1,0.0,0.0,0}};
			double[][] amps = {{0,0.8,0.1,0.0},
					{0.2,0,0.2,0.1},
					{0.2,0.0,0,0.0},
					{0.3,0.0,1.9,0}};
			double[][] phases = {{0,1.0,1.0,0.0},
					{0.5,0,0.5,0.1},
					{0.25,0.0,0,1.0},
					{0.25,0.5,1.0,0}};

			switch (config.testTreesCreateSeasonality) {
			case NONE:

				createModel = new ConstantMigrationBaseModel(Q);

				for (int i=0;i<config.numTestTrees;i++) {
					Tree testTree = new Tree(createModel,config.numTestTips);
					testTree.removeInternalLocations();
					trees.add(testTree);
				}

				/////////////

				double phase = 0.3;
				double length = 0.5;
				testModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,length);

				break;

			case TWO_CONSTANT_SEASONS:

				phase = 0.3;
				length = 0.5;
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);

				for (int i=0;i<config.numTestTrees;i++) {
					Tree testTree = new Tree(createModel,config.numTestTips);
					testTree.removeInternalLocations();
					trees.add(testTree);
				}

				/////////////

				testModel = new ConstantMigrationBaseModel(Q);
				break;

			case SINUSOIDAL:
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);

				for (int i=0;i<config.numTestTrees;i++) {
					Tree testTree = new Tree(createModel,config.numTestTips);
					testTree.removeInternalLocations();
					trees.add(testTree);
				}

				/////////////

				testModel = new ConstantMigrationBaseModel(Q);
				break;

			}
			System.out.println(" generated "+trees.size()+" trees");
			break;
		}
	}

	


}
