package seasmig.treelikelihood;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

@SuppressWarnings("serial")
public class SimpleAttributeLoader implements AttributeLoader{
	// TODO: change this to a different and single file type for attributes....
	// TODO: add grouping of locations together .....
	
	HashMap<String, Object> attributes = new HashMap<String,Object>();

	protected SimpleAttributeLoader() {};
	
	public SimpleAttributeLoader(String locationFilename,String stateFilename) throws NumberFormatException, IOException {
		if (locationFilename!=null) {
			processLocations(locationFilename);							
		}
		if (stateFilename!=null) {
			processStates(stateFilename);
		}
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

		attributes.put("locations",locationMap);
		attributes.put("numLocations",locations.size());
		
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

		attributes.put("states",stateMap);
	}

	@Override
	public Set<String> getAttributeNames() {
		return attributes.keySet();
	}

}
