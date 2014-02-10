package seasmig.treelikelihood.trees;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.ConstantTransitionBaseModel;

@SuppressWarnings("serial")
public class SimpleAttributeLoader implements AttributeLoader{
	HashMap<String, Object> attributes = new HashMap<String,Object>();

	protected SimpleAttributeLoader() {};
	
	public SimpleAttributeLoader(String locationFilename,String stateFilename, String alignmentFilename, String migrationModelFilename, String codonModelFilename) throws NumberFormatException, IOException {
		if (locationFilename!=null) {
			processLocations(locationFilename);							
		}
		if (stateFilename!=null) {
			processStates(stateFilename);
		}
		if (alignmentFilename!=null) {
			processAlignments(alignmentFilename);
		}
		if (migrationModelFilename!=null) {
			processMigrationModels(migrationModelFilename);
		}
		
		if (codonModelFilename!=null) {
			processCodonModels(codonModelFilename);
		}
	}

	
	private void processCodonModels(String filename) throws IOException {
		FileInputStream codonModelFIStream = new FileInputStream(filename);
		DataInputStream codonModelDIStream = new DataInputStream(codonModelFIStream);
		BufferedReader codonModelReader = new BufferedReader(new InputStreamReader(codonModelDIStream));

		HashMap<String, TransitionModel> codonModelMap = new HashMap<String, TransitionModel>();
	
		String strLine;
		//Read File Line By Line
		int i=0;
		while ((strLine = codonModelReader.readLine()) != null)   {
			try {
				TransitionModel model = new ConstantTransitionBaseModel(codonModelReader.readLine());
				codonModelMap.put(Integer.toString(i), model);				
			}
			catch (NumberFormatException e) {
				System.err.println("Failed to parse model on line: "+i+" of "+filename);				
			}
			i=i+1;					
		}
		
		codonModelReader.close();

		attributes.put("codonModels", codonModelMap);
	}

	private void processMigrationModels(String filename) throws IOException {

		FileInputStream migrationModelFIStream = new FileInputStream(filename);
		DataInputStream migrationModelDIStream = new DataInputStream(migrationModelFIStream);
		BufferedReader migrationModelReader = new BufferedReader(new InputStreamReader(migrationModelDIStream));

		HashMap<String, TransitionModel> migrationModelMap = new HashMap<String, TransitionModel>();
	
		String strLine;
		//Read File Line By Line
		int i=0;
		while ((strLine = migrationModelReader.readLine()) != null)   {
			try {
				TransitionModel model = new ConstantTransitionBaseModel(migrationModelReader.readLine());
				migrationModelMap.put(Integer.toString(i), model);				
			}
			catch (NumberFormatException e) {
				System.err.println("Failed to parse model on line: "+i+" of "+filename);				
			}
			i=i+1;					
		}
		
		migrationModelReader.close();

		attributes.put("migrationModels", migrationModelMap);
		
	}

	@Override
	public HashMap<String, Object> getAttributes() {
		return attributes;
	}
	
	void processLocations(String fileName) throws IOException {

		FileInputStream locationFIStream = new FileInputStream(fileName);
		DataInputStream locationDIStream = new DataInputStream(locationFIStream);
		BufferedReader locationReader = new BufferedReader(new InputStreamReader(locationDIStream));

		HashMap<String,Integer> locationMap = new HashMap<String,Integer>();
		Set<Integer> locations = new HashSet<Integer>();

		String strLine;
		//Read File Line By Line
		int i=0;
		while ((strLine = locationReader.readLine()) != null)   {
			i=i+2;
			String taxa = strLine;
			try {
				Integer state = Integer.parseInt(locationReader.readLine());
				locationMap.put(taxa, state);
				locations.add(state);
			}
			catch (NumberFormatException e) {
				System.err.println("Failed to parse trait on line: "+i+" of "+fileName);				
			}
						
		}
		
		locationReader.close();

		attributes.put("locations",locationMap);
		attributes.put("numLocations",locations.size());
		
	}
	
	void processAlignments(String fileName) throws IOException {
		// Read fasta file:
		// >header1
		// ACCCCCCAAAGTTAATATATRRYYSSSWWWW
		// AAAAAAAAAAAAA---------AACCCTGTG
		// ACG

		FileInputStream seqFIStream = new FileInputStream(fileName);
		DataInputStream seqDIStream = new DataInputStream(seqFIStream);
		BufferedReader seqReader = new BufferedReader(new InputStreamReader(seqDIStream));

		HashMap<String,Sequence> seqMap = new HashMap<String,Sequence>();
		
		String strLine;
		//Read File Line By Line
		int i=0;
		String taxa = null;
		String seq = "";
		while ((strLine = seqReader.readLine()) != null)   {
			i=i+1;
			if (strLine.substring(0,1).equals(">")) {
				if (seq.length()>0) {
					seq.replaceAll("\n", "");
					seq.replaceAll("\r", "");
					seqMap.put(taxa, new Sequence(taxa,seq));
				}				
				taxa = strLine.substring(1,strLine.length());
				seq = "";
			}
			else {
				seq = seq+strLine;
			}
		}
		
		if (seq.length()>0) {
			seq.replaceAll("\n", "");
			seq.replaceAll("\r", "");
			seqMap.put(taxa, new Sequence(taxa,seq));
		}	
		
		seqReader.close();

		attributes.put("alignments",seqMap);
		if (!seqMap.isEmpty()) {
			attributes.put("seqLength",((Sequence) seqMap.values().toArray()[0]).length());
		}
		else {
			attributes.put("seqLength",0);
		}
	}
	
	void processStates(String fileName) throws NumberFormatException, IOException {

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
		
		stateReader.close();
		attributes.put("states",stateMap);
	}

	@Override
	public Set<String> getAttributeNames() {
		return attributes.keySet();
	}

}
