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

import treelikelihood.ComplexSeasonalMigrationBaseModel;
import treelikelihood.ConstantMigrationBaseModel;
import treelikelihood.MigrationBaseModel;
import treelikelihood.Tree;

import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;

public class Data
{

	public int stateCount = 3;
	
		
	public Data(String treeFilename,String traitFilename) 	{			

		treeFile = new File("beastInput.trees");
		reader = new FileReader(treeFile);

		nexusImporter = new NexusImporter(reader);		
		nexsusTrees = nexusImporter.importTrees();

		HashMap<String,Integer> traitMap = readTraits("traits.txt");

		MigrationBaseModel constatnStateModel3 = new ConstantMigrationBaseModel(Q3S);

		System.out.println(nexsusTrees.get(0).getTaxa().iterator().next().getName());
		myTree = new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-1),traitMap,3);
		//myTree.removeInternalStates();		
		myTree.print();

		time1 = System.currentTimeMillis();
		double ll_Q3 = myTree.logLikelihood(constatnStateModel3);
		time1 = System.currentTimeMillis() - time1;
		System.out.println("\n\nLog-Likelihood using State Model: "+ll_Q3+" time: "+time1+"ms");	
	}

	public Data(String treeFilename) 	{
		File treeFile = new File("h3n2.trees");
		FileReader reader = new FileReader(treeFile);

		NexusImporter nexusImporter = new NexusImporter(reader);		
		List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();

		myTree = new Tree((SimpleRootedTree) nexsusTrees.get(nexsusTrees.size()-5),14);
		//myTree.removeInternalStates();		
		
	}	

	public Data() {
	
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
