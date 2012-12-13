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

import treelikelihood.*;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusExporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;

public class Data
{


	public Vector<Tree> trees = new Vector<Tree>();
	Config config = null;

	public Data(Config config_) throws IOException, ImportException 	{

		config = config_;		

		if (config.) {

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
			if (config.traitFilename!=null) {
				System.out.print("Loading traits... ");
				HashMap<String,Integer> traitMap = readTraits("traits.txt");
				System.out.println("loaded "+traitMap.size()+" unique traits");

				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new Tree((SimpleRootedTree) tree,traitMap,config.stateCount));
				}
				System.out.println(" reparsed "+trees.size()+" trees");
			}
			else {
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.add(new Tree((SimpleRootedTree) tree,config.stateCount));
				}		
				System.out.println(" reparsed "+trees.size()+" trees");
			}	
		}
		else {

			System.out.print("Generating test trees... ");

			switch (config.treeCreateSeasonality) {
			case NONE:

				// Generate test data and trees 
				double[][] Q1 = {{-0.9,0.4,0.3,0.2},
								{0.3,-0.6,0.2,0.1},
								{0.2,0.1,-0.3,0.0},
								{0.1,0.0,0.0,-0.1}};

				MigrationBaseModel createModel = new ConstantMigrationBaseModel(Q1);

				for (int i=0;i<config.numTestTrees;i++) {
					Tree testTree = new Tree(createModel,config.numTestTips);
					testTree.removeInternalStates();
					trees.add(testTree);
				}
				break;

			case TWO_MATRICES:

				// Generate test data and trees 
				double[][] QW = {{-0.9,0.4,0.3,0.2},
								{0.3,-0.6,0.2,0.1},
								{0.2,0.1,-0.3,0.0},
								{0.1,0.0,0.0,-0.1}};

				double[][] QS = {{-0.9,0.4,0.3,0.2},
								{0.3,-0.6,0.2,0.1},
								{0.2,0.1,-0.3,0.0},
								{0.1,0.0,0.0,-0.1}};

				double phase = 0.3;
				double length = 0.5;
				createModel = new TwoMatrixMigrationBaseModel(QW,QS,phase,length);

				for (int i=0;i<config.numTestTrees;i++) {
					Tree testTree = new Tree(createModel,config.numTestTips);
					testTree.removeInternalStates();
					trees.add(testTree);
				}
				break;
				
			case SINUSOIDAL:
				// TODO:
				System.err.println("TODO: SEASONAL MODEL\n");
				break;

			}
			System.out.println(" generated "+trees.size()+" trees");
		}
	}

	HashMap<String, Integer> readTraits(String fileName) throws NumberFormatException, IOException {

		FileInputStream traitFIStream = new FileInputStream(fileName);
		DataInputStream traitDIStream = new DataInputStream(traitFIStream);
		BufferedReader traitReader = new BufferedReader(new InputStreamReader(traitDIStream));

		HashMap<String,Integer> traitMap = new HashMap<String,Integer>();

		String strLine;
		//Read File Line By Line
		while ((strLine = traitReader.readLine()) != null)   {
			String taxa = strLine;
			Integer state = Integer.parseInt(traitReader.readLine());
			traitMap.put(taxa, state);
		}

		return traitMap;
	}


}
