package seasmig.treelikelihood.trees;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import seasmig.migrationmain.Config;
import seasmig.treelikelihood.TransitionModel;
import seasmig.treelikelihood.transitionmodels.ConstantTransitionBaseModel;

@SuppressWarnings("serial")
public class SimpleAttributeLoader implements AttributeLoader{
	HashMap<String, Object> attributes = new HashMap<String,Object>();
	private Config config;

	protected SimpleAttributeLoader() {};
	
	public SimpleAttributeLoader(Config config, String locationFilename,String stateFilename, String alignmentFilename, String migrationModelFilename, String codonModelFilename) throws NumberFormatException, IOException {
		this.config=config;
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
	
		ArrayList<String> matrixLines = new ArrayList<String>(); 
		
		//Read File Line By Line
		String strLine = null;
		int i=0;
		while ((strLine = codonModelReader.readLine()) != null)   {
			i++;
			try {
				matrixLines.add(strLine);				
			}
			catch (NumberFormatException e) {
				System.err.println("Failed to parse model on line: "+i+" of "+filename);				
			}
		}
		
		codonModelReader.close();
			
		int treeNum=0;
		for (i=Math.max(0,matrixLines.size()-config.numTreesFromTail);i<matrixLines.size();i++) {
			String[] splitString=matrixLines.get(i).split("\t");
			double k0=Double.parseDouble(splitString[1]);										
			double k1=Double.parseDouble(splitString[2]);
			double k2=Double.parseDouble(splitString[3]);
			double mu0=Double.parseDouble(splitString[4]);										
			double mu1=Double.parseDouble(splitString[5]);
			double mu2=Double.parseDouble(splitString[6]);
			double piA0; double piA1; double piA2;
			double piC0; double piC1; double piC2;
			double piT0; double piT1; double piT2;
			double piG0; double piG1; double piG2;
			piA0=(Double) attributes.get("piA0");
			piC0=(Double) attributes.get("piC0");
			piT0=(Double) attributes.get("piT0");
			piG0=(Double) attributes.get("piG0");
			
			piA1=(Double) attributes.get("piA1");
			piC1=(Double) attributes.get("piC1");
			piT1=(Double) attributes.get("piT1");
			piG1=(Double) attributes.get("piG1");
			
			piA2=(Double) attributes.get("piA2");
			piC2=(Double) attributes.get("piC2");
			piT2=(Double) attributes.get("piT2");
			piG2=(Double) attributes.get("piG2");
			
			TransitionModel[] models = new TransitionModel[3];
			models[0]=new ConstantTransitionBaseModel(mu0, k0, piC0, piA0, piG0);
			models[1]=new ConstantTransitionBaseModel(mu1, k1, piC1, piA1, piG1);
			models[2]=new ConstantTransitionBaseModel(mu2, k2, piC2, piA2, piG2);
			codonModelMap.put(Integer.toString(treeNum)+".0", models[0]);
			codonModelMap.put(Integer.toString(treeNum)+".1", models[1]);
			codonModelMap.put(Integer.toString(treeNum)+".2", models[2]);
			treeNum=treeNum+1;
		}

		if (i>0) attributes.put("codonModels", codonModelMap);	
	}

	private void processMigrationModels(String filename) throws IOException {

		FileInputStream migrationModelFIStream = new FileInputStream(filename);
		DataInputStream migrationModelDIStream = new DataInputStream(migrationModelFIStream);
		BufferedReader migrationModelReader = new BufferedReader(new InputStreamReader(migrationModelDIStream));

		HashMap<String, TransitionModel> migrationModelMap = new HashMap<String, TransitionModel>();
	
		ArrayList<String> matrixLines = new ArrayList<String>(); 
		
		//Read File Line By Line
		String strLine = null;
		int i=0;
		while ((strLine = migrationModelReader.readLine()) != null)   {
			i++;
			try {
				matrixLines.add(strLine);				
			}
			catch (NumberFormatException e) {
				System.err.println("Failed to parse model on line: "+i+" of "+filename);				
			}
		}
		
		migrationModelReader.close();
		
		//System.err.println();		
		int treeNum=0;
		for (i=Math.max(0,matrixLines.size()-config.numTreesFromTail);i<matrixLines.size();i++) {
			System.out.print(".");
			String[] splitString=matrixLines.get(i).split("\t");
			//System.err.println(splitString.length);
			double[][] Q = new double[config.numLocations][];
			for (int j=0;j<config.numLocations;j++) {
				Q[j]=new double[config.numLocations];
				for (int k=0;k<config.numLocations;k++) {
					Q[j][k]=Double.parseDouble(splitString[j*config.numLocations+k+1]);										
				}
			}
			//System.err.println("parsing migration model from iteration "+splitString[0]);
			//System.err.println(Util.parse(Q));
			TransitionModel model = new ConstantTransitionBaseModel(Q);
			migrationModelMap.put(Integer.toString(treeNum), model);
			treeNum=treeNum+1;
		}

		if (i>0) attributes.put("migrationModels", migrationModelMap);		
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
	
	static String stringTake(String str, int from, int to, int step) {
		String result=new String();
		for (int i=from;i<to;i+=step){
			result+=str.charAt(i);
		}
		return result;
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
		double countA = 0;	double countC = 0;	double countT = 0;	double countG = 0;
		double countA0 = 0;	double countC0 = 0;	double countT0 = 0;	double countG0 = 0;
		double countA1 = 0;	double countC1 = 0;	double countT1 = 0;	double countG1 = 0;
		double countA2 = 0;	double countC2 = 0;	double countT2 = 0;	double countG2 = 0;
		while ((strLine = seqReader.readLine()) != null)   {
			i=i+1;
			if (strLine.substring(0,1).equals(">")) {
				if (seq.length()>0) {
					seq.replaceAll("\n", "");
					seq.replaceAll("\r", "");
					seqMap.put(taxa, new Sequence(taxa,seq));
					
					// COUNT ACTG FREQ
					String tempA=seq; String tempC=seq; String tempG=seq; String tempT=seq;
					tempA=tempA.replaceAll("t", ""); tempA=tempA.replaceAll("T", ""); tempA=tempA.replaceAll("C", ""); tempA=tempA.replaceAll("c", "");	tempA=tempA.replaceAll("G", "");;	tempA=tempA.replaceAll("g", "");
					tempC=tempC.replaceAll("t", ""); tempC=tempC.replaceAll("T", ""); tempC=tempC.replaceAll("A", ""); tempC=tempC.replaceAll("a", "");	tempC=tempC.replaceAll("G", "");;	tempC=tempC.replaceAll("g", "");
					tempG=tempG.replaceAll("t", ""); tempG=tempG.replaceAll("T", ""); tempG=tempG.replaceAll("C", ""); tempG=tempG.replaceAll("c", "");	tempG=tempG.replaceAll("A", "");;	tempG=tempG.replaceAll("a", "");
					tempT=tempT.replaceAll("a", ""); tempT=tempT.replaceAll("A", ""); tempT=tempT.replaceAll("C", ""); tempT=tempT.replaceAll("c", "");	tempT=tempT.replaceAll("G", "");;	tempT=tempT.replaceAll("g", "");
					countA+=tempA.length();	countC+=tempC.length();	countT+=tempT.length();	countG+=tempG.length();
					//System.err.println("seq length: "+seq.length()+" ACTG length: "+(tempA.length()+tempC.length()+tempT.length()+tempG.length()));
					// for specific codon position 0
					String cp0 = stringTake(seq,0,seq.length(),3);
					String tempA0=cp0; String tempC0=cp0; String tempG0=cp0; String tempT0=cp0;
					tempA0=tempA0.replaceAll("t", ""); tempA0=tempA0.replaceAll("T", ""); tempA0=tempA0.replaceAll("C", ""); tempA0=tempA0.replaceAll("c", "");	tempA0=tempA0.replaceAll("G", "");;	tempA0=tempA0.replaceAll("g", "");
					tempC0=tempC0.replaceAll("t", ""); tempC0=tempC0.replaceAll("T", ""); tempC0=tempC0.replaceAll("A", ""); tempC0=tempC0.replaceAll("a", "");	tempC0=tempC0.replaceAll("G", "");;	tempC0=tempC0.replaceAll("g", "");
					tempG0=tempG0.replaceAll("t", ""); tempG0=tempG0.replaceAll("T", ""); tempG0=tempG0.replaceAll("C", ""); tempG0=tempG0.replaceAll("c", "");	tempG0=tempG0.replaceAll("A", "");;	tempG0=tempG0.replaceAll("a", "");
					tempT0=tempT0.replaceAll("a", ""); tempT0=tempT0.replaceAll("A", ""); tempT0=tempT0.replaceAll("C", ""); tempT0=tempT0.replaceAll("c", "");	tempT0=tempT0.replaceAll("G", "");;	tempT0=tempT0.replaceAll("g", "");
					countA0+=tempA0.length();	countC0+=tempC0.length();	countT0+=tempT0.length();	countG0+=tempG0.length();
					//System.err.println("seq length: "+seq.length()+" ACTG length: "+(tempA.length()+tempC.length()+tempT.length()+tempG.length()));
					// for specific codon position 1
					String cp1 = stringTake(seq,1,seq.length(),3);
					String tempA1=cp1; String tempC1=cp1; String tempG1=cp1; String tempT1=cp1;
					tempA1=tempA1.replaceAll("t", ""); tempA1=tempA1.replaceAll("T", ""); tempA1=tempA1.replaceAll("C", ""); tempA1=tempA1.replaceAll("c", "");	tempA1=tempA1.replaceAll("G", "");;	tempA1=tempA1.replaceAll("g", "");
					tempC1=tempC1.replaceAll("t", ""); tempC1=tempC1.replaceAll("T", ""); tempC1=tempC1.replaceAll("A", ""); tempC1=tempC1.replaceAll("a", "");	tempC1=tempC1.replaceAll("G", "");;	tempC1=tempC1.replaceAll("g", "");
					tempG1=tempG1.replaceAll("t", ""); tempG1=tempG1.replaceAll("T", ""); tempG1=tempG1.replaceAll("C", ""); tempG1=tempG1.replaceAll("c", "");	tempG1=tempG1.replaceAll("A", "");;	tempG1=tempG1.replaceAll("a", "");
					tempT1=tempT1.replaceAll("a", ""); tempT1=tempT1.replaceAll("A", ""); tempT1=tempT1.replaceAll("C", ""); tempT1=tempT1.replaceAll("c", "");	tempT1=tempT1.replaceAll("G", "");;	tempT1=tempT1.replaceAll("g", "");
					countA1+=tempA1.length();	countC1+=tempC1.length();	countT1+=tempT1.length();	countG1+=tempG1.length();
					// for specific codon position 2
					String cp2 = stringTake(seq,2,seq.length(),3);
					String tempA2=cp2; String tempC2=cp2; String tempG2=cp2; String tempT2=cp2;
					tempA2=tempA2.replaceAll("t", ""); tempA2=tempA2.replaceAll("T", ""); tempA2=tempA2.replaceAll("C", ""); tempA2=tempA2.replaceAll("c", "");	tempA2=tempA2.replaceAll("G", "");;	tempA2=tempA2.replaceAll("g", "");
					tempC2=tempC2.replaceAll("t", ""); tempC2=tempC2.replaceAll("T", ""); tempC2=tempC2.replaceAll("A", ""); tempC2=tempC2.replaceAll("a", "");	tempC2=tempC2.replaceAll("G", "");;	tempC2=tempC2.replaceAll("g", "");
					tempG2=tempG2.replaceAll("t", ""); tempG2=tempG2.replaceAll("T", ""); tempG2=tempG2.replaceAll("C", ""); tempG2=tempG2.replaceAll("c", "");	tempG2=tempG2.replaceAll("A", "");;	tempG2=tempG2.replaceAll("a", "");
					tempT2=tempT2.replaceAll("a", ""); tempT2=tempT2.replaceAll("A", ""); tempT2=tempT2.replaceAll("C", ""); tempT2=tempT2.replaceAll("c", "");	tempT2=tempT2.replaceAll("G", "");;	tempT2=tempT2.replaceAll("g", "");
					countA2+=tempA2.length();	countC2+=tempC2.length();	countT2+=tempT2.length();	countG2+=tempG2.length();
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
			
			// COUNT ACTG FREQ
			String tempA=seq; String tempC=seq; String tempG=seq; String tempT=seq;
			tempA=tempA.replaceAll("t", ""); tempA=tempA.replaceAll("T", ""); tempA=tempA.replaceAll("C", ""); tempA=tempA.replaceAll("c", "");	tempA=tempA.replaceAll("G", "");;	tempA=tempA.replaceAll("g", "");
			tempC=tempC.replaceAll("t", ""); tempC=tempC.replaceAll("T", ""); tempC=tempC.replaceAll("A", ""); tempC=tempC.replaceAll("a", "");	tempC=tempC.replaceAll("G", "");;	tempC=tempC.replaceAll("g", "");
			tempG=tempG.replaceAll("t", ""); tempG=tempG.replaceAll("T", ""); tempG=tempG.replaceAll("C", ""); tempG=tempG.replaceAll("c", "");	tempG=tempG.replaceAll("A", "");;	tempG=tempG.replaceAll("a", "");
			tempT=tempT.replaceAll("a", ""); tempT=tempT.replaceAll("A", ""); tempT=tempT.replaceAll("C", ""); tempT=tempT.replaceAll("c", "");	tempT=tempT.replaceAll("G", "");;	tempT=tempT.replaceAll("g", "");
			countA+=tempA.length();	countC+=tempC.length();	countT+=tempT.length();	countG+=tempG.length();
			// for specific codon position 0
			String cp0 = stringTake(seq,0,seq.length(),3);
			String tempA0=cp0; String tempC0=cp0; String tempG0=cp0; String tempT0=cp0;
			tempA0=tempA0.replaceAll("t", ""); tempA0=tempA0.replaceAll("T", ""); tempA0=tempA0.replaceAll("C", ""); tempA0=tempA0.replaceAll("c", "");	tempA0=tempA0.replaceAll("G", "");;	tempA0=tempA0.replaceAll("g", "");
			tempC0=tempC0.replaceAll("t", ""); tempC0=tempC0.replaceAll("T", ""); tempC0=tempC0.replaceAll("A", ""); tempC0=tempC0.replaceAll("a", "");	tempC0=tempC0.replaceAll("G", "");;	tempC0=tempC0.replaceAll("g", "");
			tempG0=tempG0.replaceAll("t", ""); tempG0=tempG0.replaceAll("T", ""); tempG0=tempG0.replaceAll("C", ""); tempG0=tempG0.replaceAll("c", "");	tempG0=tempG0.replaceAll("A", "");;	tempG0=tempG0.replaceAll("a", "");
			tempT0=tempT0.replaceAll("a", ""); tempT0=tempT0.replaceAll("A", ""); tempT0=tempT0.replaceAll("C", ""); tempT0=tempT0.replaceAll("c", "");	tempT0=tempT0.replaceAll("G", "");;	tempT0=tempT0.replaceAll("g", "");
			countA0+=tempA0.length();	countC0+=tempC0.length();	countT0+=tempT0.length();	countG0+=tempG0.length();
			//System.err.println("seq length: "+seq.length()+" ACTG length: "+(tempA.length()+tempC.length()+tempT.length()+tempG.length()));
			// for specific codon position 1
			String cp1 = stringTake(seq,1,seq.length(),3);
			String tempA1=cp1; String tempC1=cp1; String tempG1=cp1; String tempT1=cp1;
			tempA1=tempA1.replaceAll("t", ""); tempA1=tempA1.replaceAll("T", ""); tempA1=tempA1.replaceAll("C", ""); tempA1=tempA1.replaceAll("c", "");	tempA1=tempA1.replaceAll("G", "");;	tempA1=tempA1.replaceAll("g", "");
			tempC1=tempC1.replaceAll("t", ""); tempC1=tempC1.replaceAll("T", ""); tempC1=tempC1.replaceAll("A", ""); tempC1=tempC1.replaceAll("a", "");	tempC1=tempC1.replaceAll("G", "");;	tempC1=tempC1.replaceAll("g", "");
			tempG1=tempG1.replaceAll("t", ""); tempG1=tempG1.replaceAll("T", ""); tempG1=tempG1.replaceAll("C", ""); tempG1=tempG1.replaceAll("c", "");	tempG1=tempG1.replaceAll("A", "");;	tempG1=tempG1.replaceAll("a", "");
			tempT1=tempT1.replaceAll("a", ""); tempT1=tempT1.replaceAll("A", ""); tempT1=tempT1.replaceAll("C", ""); tempT1=tempT1.replaceAll("c", "");	tempT1=tempT1.replaceAll("G", "");;	tempT1=tempT1.replaceAll("g", "");
			countA1+=tempA1.length();	countC1+=tempC1.length();	countT1+=tempT1.length();	countG1+=tempG1.length();
			// for specific codon position 2
			String cp2 = stringTake(seq,2,seq.length(),3);
			String tempA2=cp2; String tempC2=cp2; String tempG2=cp2; String tempT2=cp2;
			tempA2=tempA2.replaceAll("t", ""); tempA2=tempA2.replaceAll("T", ""); tempA2=tempA2.replaceAll("C", ""); tempA2=tempA2.replaceAll("c", "");	tempA2=tempA2.replaceAll("G", "");;	tempA2=tempA2.replaceAll("g", "");
			tempC2=tempC2.replaceAll("t", ""); tempC2=tempC2.replaceAll("T", ""); tempC2=tempC2.replaceAll("A", ""); tempC2=tempC2.replaceAll("a", "");	tempC2=tempC2.replaceAll("G", "");;	tempC2=tempC2.replaceAll("g", "");
			tempG2=tempG2.replaceAll("t", ""); tempG2=tempG2.replaceAll("T", ""); tempG2=tempG2.replaceAll("C", ""); tempG2=tempG2.replaceAll("c", "");	tempG2=tempG2.replaceAll("A", "");;	tempG2=tempG2.replaceAll("a", "");
			tempT2=tempT2.replaceAll("a", ""); tempT2=tempT2.replaceAll("A", ""); tempT2=tempT2.replaceAll("C", ""); tempT2=tempT2.replaceAll("c", "");	tempT2=tempT2.replaceAll("G", "");;	tempT2=tempT2.replaceAll("g", "");
			countA2+=tempA2.length();	countC2+=tempC2.length();	countT2+=tempT2.length();	countG2+=tempG2.length();
		}	
		
		seqReader.close();

		attributes.put("alignments",seqMap);
		if (!seqMap.isEmpty()) {
			attributes.put("seqLength",((Sequence) seqMap.values().toArray()[0]).length());
			attributes.put("piA", countA/(countA+countC+countT+countG));
			attributes.put("piC", countC/(countA+countC+countT+countG));
			attributes.put("piT", countT/(countA+countC+countT+countG));
			attributes.put("piG", countG/(countA+countC+countT+countG));
			// for specific codon position
			attributes.put("piA0", countA0/(countA0+countC0+countT0+countG0));
			attributes.put("piC0", countC0/(countA0+countC0+countT0+countG0));
			attributes.put("piT0", countT0/(countA0+countC0+countT0+countG0));
			attributes.put("piG0", countG0/(countA0+countC0+countT0+countG0));
			// for specific codon position
			attributes.put("piA1", countA1/(countA1+countC1+countT1+countG1));
			attributes.put("piC1", countC1/(countA1+countC1+countT1+countG1));
			attributes.put("piT1", countT1/(countA1+countC1+countT1+countG1));
			attributes.put("piG1", countG1/(countA1+countC1+countT1+countG1));
			// for specific codon position
			attributes.put("piA2", countA2/(countA2+countC2+countT2+countG2));
			attributes.put("piC2", countC2/(countA2+countC2+countT2+countG2));
			attributes.put("piT2", countT2/(countA2+countC2+countT2+countG2));
			attributes.put("piG2", countG2/(countA2+countC2+countT2+countG2));
			
			System.out.println("empirical base freq (ACTG)");
			System.out.println("cp0=("+(Double) attributes.get("piA0")+","+(Double) attributes.get("piC0")+","+(Double) attributes.get("piT0")+","+(Double) attributes.get("piG0")+")");
			System.out.println("cp1=("+(Double) attributes.get("piA1")+","+(Double) attributes.get("piC1")+","+(Double) attributes.get("piT1")+","+(Double) attributes.get("piG1")+")");
			System.out.println("cp2=("+(Double) attributes.get("piA2")+","+(Double) attributes.get("piC2")+","+(Double) attributes.get("piT2")+","+(Double) attributes.get("piG2")+")");
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
