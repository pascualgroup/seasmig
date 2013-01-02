package treelikelihood;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SimpleAttributeLoader implements AttributeLoader{
	// TODO: change this to a different and single file type for attributes....
	// TODO: add grouping of locations together .....
	
	HashMap<String, Object> attributes = new HashMap<String,Object>();
	List<String> attributeNames = new ArrayList<String>();

	public SimpleAttributeLoader(String locationFilename,String stateFilename) throws NumberFormatException, IOException {
		if (locationFilename!=null) {
			attributes.put("locations",readLocations(locationFilename));
			attributeNames.add("locations");
		}
		if (stateFilename!=null) {
			attributes.put("states",readStates(stateFilename));
			attributeNames.add("states");
		}
	}

	@Override
	public HashMap<String, Object> getAttributes() {
		return attributes;
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

	@Override
	public List<String> getAttributeNames() {
		return attributeNames;
	}

}
