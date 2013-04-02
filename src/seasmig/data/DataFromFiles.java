package seasmig.data;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectStreamException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.JsonSyntaxException;

import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.AttributeLoader;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.SimpleAttributeLoader;
import seasmig.treelikelihood.TreeWithLocations;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.trees.SimpleRootedTree;

@SuppressWarnings("serial")
public class DataFromFiles implements Data
{
	public List<ArrayList<LikelihoodTree>> trees;
	Config config = null;
	int numLocations;

	protected DataFromFiles() {};

	private void writeObject(java.io.ObjectOutputStream out) throws IOException {		 
	}

	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException {
		Gson gson = new GsonBuilder().setPrettyPrinting().create();	
		config = gson.fromJson(new FileReader("config.json"), Config.class);	
		try {
			loadDataFromFiles(config);
		} catch (ImportException e) {
			System.err.println("Failed to unpack config.json");
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unused")
	private void readObjectNoData() throws ObjectStreamException {
		Gson gson = new GsonBuilder().setPrettyPrinting().create();	
		try {
			config = gson.fromJson(new FileReader("config.json"), Config.class);
		} catch (JsonSyntaxException e) {
			System.err.println("Failed to unpack config.json JsonSyntaxException ");			
			e.printStackTrace();
		} catch (JsonIOException e) {
			System.err.println("Failed to unpack config.json JsonIOException");
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			System.err.println("Failed to unpack config.json FileNotFoundException");
			e.printStackTrace();
		}	
		try {
			loadDataFromFiles(config);			
		} catch (IOException e) {
			System.err.println("Failed to unpack config.json IOException");
			e.printStackTrace();
			e.printStackTrace();
		} catch (ImportException e) {
			System.err.println("Failed to unpack config.json ImportExceptoin");
			e.printStackTrace();
			e.printStackTrace();		}
	}


	public DataFromFiles(Config config_) throws IOException, ImportException 	{
		loadDataFromFiles(config_);
	}

	private void loadDataFromFiles(Config config_) throws IOException, ImportException {
		config = config_;
		trees = new ArrayList<ArrayList<LikelihoodTree>>();
		numLocations=0;

		// Load trees
		System.out.println("Loading trees... ");
		for (int h=0;h<config.treeFilenames.length;h++) {
			trees.add(new ArrayList<LikelihoodTree>());
			String treeFilename = config.treeFilenames[h];				
			System.out.println("Loading tree file: "+treeFilename);

			File treeFile = new File(treeFilename);
			FileReader reader = new FileReader(treeFile);
			NexusImporter nexusImporter = new NexusImporter(reader);
			List<jebl.evolution.trees.Tree> nexsusTrees = nexusImporter.importTrees();
			System.out.println("loaded "+nexsusTrees.size()+" trees");

			System.out.print("Keeping tail... ");
			double meanNumTaxa=0;
			List<jebl.evolution.trees.Tree> nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
				nexsusTreeTail.add(nexsusTrees.get(i));
				meanNumTaxa+=nexsusTrees.get(i).getTaxa().size();
			}
			meanNumTaxa/=nexsusTreeTail.size();
			System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");
			System.out.println(meanNumTaxa+" taxa on average per tree");

			// TODO: add states....
			// Convert trees to internal tree representation
			if (config.locationFilenames[h]!=null) {
				System.out.print("Loading traits... ");
				String locationFilename = config.locationFilenames[h];
				AttributeLoader attributeLoader= new SimpleAttributeLoader(locationFilename, config.stateFilename);	
				HashMap<String,Object> attributes = attributeLoader.getAttributes();
				HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributes.get("locations");
				HashMap<String,Double> stateMap = (HashMap<String,Double>) attributes.get("states");
				numLocations = (Integer) attributes.get("numLocations");
				System.out.println("loaded "+locationMap.size()+" taxon traits");

				System.out.print("Reparsing trees... ");
				if (stateMap==null) {
					double numIdentifiedLocations=0;
					for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
						trees.get(h).add(new TreeWithLocations((SimpleRootedTree) tree,locationMap,numLocations));
						numIdentifiedLocations+=((TreeWithLocations)trees.get(h).get(trees.get(h).size()-1)).getNumIdentifiedLocations();
					}
					numIdentifiedLocations=numIdentifiedLocations/trees.get(h).size();
					System.out.print("identified "+numIdentifiedLocations+" tip locations on average per tree");
				}
				else {
					// TODO: this...
				}
				System.out.println(" reparsed "+trees.get(h).size()+" trees");
			}
			else {
				// TODO: add load states from trees...
				numLocations=config.numLocations; // TODO: get this to be automatically loaded from trees
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.get(h).add(new TreeWithLocations((SimpleRootedTree) tree,config.locationAttributeNameInTree, numLocations));
				}		
				System.out.println(" reparsed "+trees.get(h).size()+" trees");
			}	
		}

	}

	@Override
	public int getNumLocations() {		
		return numLocations;
	}

	@Override
	public List<ArrayList<LikelihoodTree>> getTrees() {
		return trees;
	}

}
