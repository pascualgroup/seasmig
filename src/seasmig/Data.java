package seasmig;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import cern.colt.function.DoubleFunction;

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

			// Convert trees to internal tree representation
			if (config.locationFilename!=null) {
				System.out.print("Loading traits... ");
				HashMap<String,Integer> locationMap = readLocations(config.locationFilename);
				System.out.println("loaded "+locationMap.size()+" taxon locations");

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

	HashMap<String, Integer> readLocations(String fileName) throws NumberFormatException, IOException {

		FileInputStream locationFIStream = new FileInputStream(fileName);
		DataInputStream locationDIStream = new DataInputStream(locationFIStream);
		BufferedReader locationReader = new BufferedReader(new InputStreamReader(locationDIStream));

		HashMap<String,Integer> locationMap = new HashMap<String,Integer>();

		String strLine;
		//Read File Line By Line
		while ((strLine = locationReader.readLine()) != null)   {
			String taxa = strLine;
			Integer state = Integer.parseInt(locationReader.readLine());
			locationMap.put(taxa, state);
		}

		return locationMap;
	}
	
	HashMap<String, Double> readStates(String fileName) throws NumberFormatException, IOException {

		FileInputStream stateFIStream = new FileInputStream(fileName);
		DataInputStream stateDIStream = new DataInputStream(stateFIStream);
		BufferedReader stateReader = new BufferedReader(new InputStreamReader(stateDIStream));

		HashMap<String,Double> stateMap = new HashMap<String,Double>();

		String strLine;
		//Read File Line By Line
		while ((strLine = stateReader.readLine()) != null)   {
			String taxa = strLine;
			Double state = Double.parseDouble(stateReader.readLine());
			stateMap.put(taxa, state);
		}

		return stateMap;
	}


}
