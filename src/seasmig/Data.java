package seasmig;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
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

	public Data(Config config_) throws IOException, ImportException 	{

		config = config_;		

		if (config.treeFilename!=null) {

			// Load trees
			File treeFile = new File(config.treeFilename);
			FileReader reader = new FileReader(treeFile);
			NexusImporter nexusImporter = new NexusImporter(reader);		
			List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();

			// Convert trees to internal tree representation
			if (config.traitFilename!=null) {
				HashMap<String,Integer> traitMap = readTraits("traits.txt");

				int treeLoadStartIndex = Math.max(0,nexsusTrees.size()-config.numTreesFromTail); 
				for (int i=treeLoadStartIndex;i<nexsusTrees.size();i++); {
					trees.add(new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-1),traitMap,config.stateCount));
				}
			}
			else {

				int treeLoadStartIndex = Math.max(0,nexsusTrees.size()-config.numTreesFromTail);
				for (int i=treeLoadStartIndex;i<nexsusTrees.size();i++); {
					trees.add(new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-1),config.stateCount));
				}		
			}	
		}
		else {
			// Generate test data and trees (random 20 trees with 4 states)
			double[][] Q1 = {{-0.75,0.25,0.25,0.25},{0.25,-0.75,0.25,0.25},{0.25,0.25,-0.75,0.25},{0.25,0.25,0.25,-0.75}};

			MigrationBaseModel createModel = new ConstantMigrationBaseModel(Q1);

			for (int i=0;i<20;i++) {
				Tree myTree = new Tree(createModel,1000);
				myTree.removeInternalStates();
			}	
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
