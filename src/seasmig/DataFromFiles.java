package seasmig;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;

import treelikelihood.*;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.math.Random;

@SuppressWarnings("serial")
public class DataFromFiles implements Data
{
	public Vector<LikelihoodTree> trees = new Vector<LikelihoodTree>();
	Config config = null;
	int numLocations = 0;

	protected DataFromFiles() {};
	
	public DataFromFiles(Config config_) throws IOException, ImportException 	{

		config = config_;		

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
			numLocations = (Integer) attributes.get("numLocations");
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
			numLocations=config.numLocations; // TODO: get this to be automatically loaded from trees
			System.out.print("Reparsing trees... ");
			for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
				trees.add(new TreeWithLocations((SimpleRootedTree) tree,config.locationAttributeNameInTree, numLocations));
			}		
			System.out.println(" reparsed "+trees.size()+" trees");
		}	


	}

	@Override
	public int getNumLocations() {		
		return numLocations;
	}

	@Override
	public List<LikelihoodTree> getTrees() {
		return trees;
	}
	
}
