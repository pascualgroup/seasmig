package seasmig.data;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectStreamException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import seasmig.migrationmain.Config;
import seasmig.data.Data;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.trees.AttributeLoader;
import seasmig.treelikelihood.trees.Sequence;
import seasmig.treelikelihood.trees.SimpleAttributeLoader;
import seasmig.treelikelihood.trees.TreeWithLocations;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonIOException;
import com.google.gson.JsonSyntaxException;

@SuppressWarnings("serial")
public class DataFromFiles implements Data
{
	public List<ArrayList<LikelihoodTree>> trees;
	Config config = null;
	int numLocations;
	long iteration = 0;

	protected DataFromFiles() {};

	private void writeObject(java.io.ObjectOutputStream out) throws IOException {
		// Data from files isn't serialized it is reloaded...

		// TODO: move report of iteration to somewhere else...
		iteration += config.checkpointEvery;
		System.out.print("\riteration: "+iteration);
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
			List<Taxon> taxa = nexusImporter.parseTaxaBlock();		
			HashMap<String,Integer> taxaIndices = new HashMap<String,Integer>();			
			for (int i=0;i<taxa.size();i++) {
				taxaIndices.put(taxa.get(i).getName(), i);
			}

			List<jebl.evolution.trees.Tree> nexusTrees = nexusImporter.importTrees();			
			System.out.println("loaded "+nexusTrees.size()+" trees");
			//System.err.println(nexusTrees.get(0).getAttribute("posterior"));

			System.out.print("Keeping tail... ");
			double meanNumTaxa=0;
			List<jebl.evolution.trees.Tree> nexusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexusTrees.size()-config.numTreesFromTail);i<nexusTrees.size();i++) {
				nexusTreeTail.add(nexusTrees.get(i));
				meanNumTaxa+=nexusTrees.get(i).getTaxa().size();
			}
			meanNumTaxa/=nexusTreeTail.size();
			System.out.println(" keeping last "+nexusTreeTail.size()+ " trees");
			System.out.println(meanNumTaxa+" taxa on average per tree");

			// Convert trees to internal tree representation
			String locationFilename = null;
			if (config.locationFilenames[h]!=null) {
				if (config.locationFilenames[h].length()>0) {
					locationFilename = config.locationFilenames[h];
				}
			}

			String alignmentFilename = null;
			if (config.alignmentFilenames!=null) {
				if (config.alignmentFilenames.length>0) {
					if (config.alignmentFilenames[h]!=null) {
						if (config.alignmentFilenames[h].length()>0) {
							alignmentFilename = config.alignmentFilenames[h];
						}
					}
				}
			}

			// Convert trees to internal tree representation
			String migrationModelFilename = null;
			if (config.migrationModelFilenames!=null) {
				if (config.migrationModelFilenames[h]!=null) {
					if (config.migrationModelFilenames[h].length()>0) {
						migrationModelFilename = config.migrationModelFilenames[h];
					}
				}
			}

			// Convert trees to internal tree representation
			String codonModelFilename = null;
			if (config.codonModelFilenames!=null) {
				if (config.codonModelFilenames[h]!=null) {
					if (config.codonModelFilenames[h].length()>0) {
						codonModelFilename = config.codonModelFilenames[h];
					}
				}
			}

			System.out.print("Loading traits, alignments, migration and codon models... ");

			AttributeLoader attributeLoader= new SimpleAttributeLoader(config, locationFilename, null, alignmentFilename, migrationModelFilename, codonModelFilename);	
			HashMap<String,Object> attributes = attributeLoader.getAttributes();
			@SuppressWarnings("unchecked")
			HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributes.get("locations");
			@SuppressWarnings("unchecked") 
			HashMap<String,Sequence> seqMap = (HashMap<String,Sequence>) attributes.get("alignments");
			@SuppressWarnings("unchecked")
			HashMap<String,TransitionModel> migrationModelMap = (HashMap<String,TransitionModel>) attributes.get("migrationModels");
			@SuppressWarnings("unchecked")
			HashMap<String,TransitionModel> codonModelMap = (HashMap<String,TransitionModel>) attributes.get("codonModels");		

			numLocations = (Integer) attributes.get("numLocations");
			System.out.println("loaded "+locationMap.size()+" taxon traits"+" from "+locationFilename);	
			Integer seqLength = (Integer) attributes.get("seqLength");
			if (seqLength==null) 
				seqLength=0;
			if (seqMap!=null)
				System.out.println("loaded "+seqMap.size()+" sequences"+" from "+alignmentFilename);

			if (codonModelMap!=null)
				System.out.println("loaded "+codonModelMap.size()/3+" codon models"+" from "+codonModelFilename);
			
			if (migrationModelMap!=null)
				System.out.println("loaded "+migrationModelMap.size()+" migration models"+" from "+migrationModelFilename);
			
			System.out.print("Reparsing trees... ");

			if (seqMap!=null) {
				if (seqMap.size()>0)
					seqLength=seqMap.get(seqMap.keySet().iterator().next()).length();
			}

			double numIdentifiedLocations=0;
			double numIdentifiedSeqs=0;
			int numid = 0;
			for (jebl.evolution.trees.Tree tree : nexusTreeTail) {
				TransitionModel migrationModel = null;
				if (migrationModelMap!=null) 
					migrationModel = migrationModelMap.get(Integer.toString(numid));
				TransitionModel[] codonModel = null;
				if (codonModelMap!=null) {
					codonModel = new TransitionModel[3];
					codonModel[0]=codonModelMap.get(Integer.toString(numid)+"."+"0");
					codonModel[1]=codonModelMap.get(Integer.toString(numid)+"."+"1");
					codonModel[2]=codonModelMap.get(Integer.toString(numid)+"."+"2");					
				}
				trees.get(h).add(new TreeWithLocations((SimpleRootedTree) tree, taxaIndices, locationMap,numLocations,config.lastTipTime[h],seqMap,seqLength,config, migrationModel, codonModel));
				numid++;								
				if (numid%10==0) System.out.print(".");
				numIdentifiedLocations+=((TreeWithLocations)trees.get(h).get(trees.get(h).size()-1)).getNumIdentifiedLocations();
				numIdentifiedSeqs+=((TreeWithLocations)trees.get(h).get(trees.get(h).size()-1)).getNumIdentifiedSeqs();
			}
			System.out.println("done!");
			numIdentifiedLocations=numIdentifiedLocations/trees.get(h).size();
			numIdentifiedSeqs=numIdentifiedSeqs/trees.get(h).size();
			System.out.println("identified "+numIdentifiedLocations+" locations on average per tree");
			System.out.println("identified "+numIdentifiedSeqs+" sequences on average per tree");

			System.out.println(" reparsed "+trees.get(h).size()+" trees");			
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
