package seasmig;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import treelikelihood.ComplexSeasonalMigrationBaseModel;
import treelikelihood.ConstantMigrationBaseModel;
import treelikelihood.MigrationBaseModel;
import treelikelihood.Tree;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;

public class Data
{

	public int stateCount = 3;
	public Vector<Tree> trees = new Vector<Tree>();

	public Data(String treeFilename,String traitFilename) throws IOException, ImportException 	{			

		File treeFile = new File("beastInput.trees");
		FileReader reader = new FileReader(treeFile);

		NexusImporter nexusImporter = new NexusImporter(reader);		
		List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();

		HashMap<String,Integer> traitMap = readTraits("traits.txt");

		int treeLoadStartIndex = Math.max(0,trees.size()-20); //TODO: add config param,
		for (int i=treeLoadStartIndex;i<nexsusTrees.size();i++); {
			trees.add(new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-1),traitMap,stateCount));
		}

	}

	public Data(String treeFilename) throws IOException, ImportException 	{
		File treeFile = new File(treeFilename);
		FileReader reader = new FileReader(treeFile);

		NexusImporter nexusImporter = new NexusImporter(reader);		
		List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();

		int treeLoadStartIndex = Math.max(0,trees.size()-20); //TODO: add config param,
		for (int i=treeLoadStartIndex;i<nexsusTrees.size();i++); {
			trees.add(new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-1),stateCount));
		}		
	}	

	public Data() {
		
		// TODO:		
		// generate data...		

	}

	static HashMap<String, Integer> readTraits(String fileName) throws NumberFormatException, IOException {

		FileInputStream traitFIStream = new FileInputStream(fileName);
		DataInputStream traitDIStream = new DataInputStream(traitFIStream);
		BufferedReader traitReader = new BufferedReader(new InputStreamReader(traitDIStream));

		HashMap<String,Integer> traitMap = new HashMap<String,Integer>();

		String strLine;
		//Read File Line By Line
		while ((strLine = traitReader.readLine()) != null)   {
			// Print the content on the console
			//System.out.println (strLine);
			String taxa = strLine;
			Integer state = Integer.parseInt(traitReader.readLine());
			traitMap.put(taxa, state);
		}

		return traitMap;
	}


}
