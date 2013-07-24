package seasmig.data;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.List;

import cern.colt.function.DoubleFunction;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.trees.SimpleRootedTree;
import jebl.evolution.trees.Tree;
import jebl.math.Random;
import seasmig.Config;
import seasmig.Data;
import seasmig.treelikelihood.AttributeLoader;
import seasmig.treelikelihood.ConstantMigrationBaseModel;
import seasmig.treelikelihood.GeneralSeasonalMigrationBaseModel;
import seasmig.treelikelihood.LikelihoodTree;
import seasmig.treelikelihood.MigrationBaseModel;
import seasmig.treelikelihood.SimpleAttributeLoader;
import seasmig.treelikelihood.SinusoidialSeasonalMigrationBaseModel;
import seasmig.treelikelihood.TreeWithLocations;
import seasmig.treelikelihood.TwoSeasonMigrationBaseModel;

@SuppressWarnings("serial")
public class DataForTests implements Data {
	// TODO: This...
	public List<ArrayList<LikelihoodTree>> trees = new ArrayList<ArrayList<LikelihoodTree>>();
	Config config = null;
	long iteration = 0;
	public enum TestType {TEST_USING_GENERATED_TREES, TEST_USING_INPUT_TREES, TEST_MODEL_DEGENERACY};

	// TEST MODELS
	public MigrationBaseModel createModel = null;
	public List<MigrationBaseModel> testModels = null;
	private int numTestTips=3000;
	private int numLocations;
	private int disturbanceScale;

	protected DataForTests() {};
	
	private void writeObject(java.io.ObjectOutputStream out) throws IOException {	
		// TODO: Support Serialization of Test Data...
		// TODO: move report of iteration to somewhere else...
		iteration += config.checkpointEvery;
		System.out.print("\riteration: "+iteration);
	}

	public DataForTests(Config config_, TestType testType, int numTestRepeats, int numTestLocations, int numTestTrees) throws IOException, ImportException 	{				

		config = config_;		

		File treeFile;
		FileReader reader;
		NexusImporter nexusImporter;
		List<Tree> nexsusTrees;
		ArrayList<Tree> nexsusTreeTail;
		switch (testType) {

		case TEST_USING_GENERATED_TREES:

			numLocations=numTestLocations;
			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Generating test trees... ");

			// Generate test data and trees

			// For constant model...
			double[][] Q = makeRandomMigrationMatrix(numTestLocations,5); 
//			Q = new double[][]{{-1,0.333333,0.333333,0.333333},
//					 { 0.333333,-1,0.333333,0.333333},
//					 { 0.333333,0.333333,-1,0.333333},
//					 { 0.333333,0.333333,0.333333,-1}};

			// For two seasonal model...
			double[][] QW = myMatrixCopy(Q);
			double[][] QS = makeRandomMigrationMatrix(numTestLocations,3); 

			// For sinusoidal model...
			double[][] rates = myMatrixCopy(Q);
			double[][] amps = makeRandomMigrationMatrix(numTestLocations,1);
			double[][] phases = makeRandomMigrationMatrix(numTestLocations,1);

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: 
				double phase = 0.3; double length = 0.5;
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL:
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			trees.add(new ArrayList<LikelihoodTree>());
			for (int i=0;i<=numTestTrees;i++) {
				TreeWithLocations testTree = new TreeWithLocations(createModel,numTestTips);
				testTree.removeInternalLocations();
				trees.get(0).add(testTree);
			}

			System.out.println(" generated "+trees.get(0).size()+" trees");
			System.out.print("Generating test models... ");

			testModels = new ArrayList<MigrationBaseModel>();

			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new ConstantMigrationBaseModel(disturbMigrationMatrix(Q, disturbanceScale*i,99999)));
			}
			for (int i=0; i<numTestRepeats; i++) {
				double phase = Math.max(0,Math.min(1,0.3+i/3*(Random.nextDouble()-0.5))); double length = 0.5;
				testModels.add(new TwoSeasonMigrationBaseModel(disturbMigrationMatrix(QW,disturbanceScale*i/3,99999),disturbMigrationMatrix(QS,disturbanceScale*i/3,99999),phase, phase+length));
			}
			
			DoubleFunction[] rootFreqFunction = new DoubleFunction[Q.length]; 
			DoubleFunction[][] generalMigrationFunction = new DoubleFunction[Q.length][Q[0].length]; 
			for (int i=0;i<Q.length;i++) {
				rootFreqFunction[i]=cern.jet.math.Functions.constant(1.0/Q.length);
				for (int j=0;j<Q[0].length;j++) {
					generalMigrationFunction[i][j]=cern.jet.math.Functions.constant(Q[i][j]);
				}
			}
								
			for (int i=0; i<numTestRepeats; i++) {
				testModels.add(new GeneralSeasonalMigrationBaseModel(generalMigrationFunction,rootFreqFunction,2));
			}
			break;

		case TEST_USING_INPUT_TREES:{

			numLocations=numTestLocations;
			// TODO: get number of locations from trees...

			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Generating test trees based on input tree topology ... ");

			// Generate test data and trees

			// For constant model...
			Q = makeRandomMigrationMatrix(numTestLocations,2); 
			
			// For two seasonal model...
			QW = myMatrixCopy(Q);
			QS = makeRandomMigrationMatrix(numTestLocations,3); 

			// For sinusoidal model...
			rates = myMatrixCopy(Q);
			amps = makeRandomMigrationMatrix(numTestLocations,1);
			phases = makeRandomMigrationMatrix(numTestLocations,1);

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: 
				double phase = 0.3; double length = 0.5;
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL: 
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			System.out.print("Loading trees... ");

			for (int h=0;h<config.treeFilenames.length;h++) {
				trees.add(new ArrayList<LikelihoodTree>());
				treeFile = new File(config.treeFilenames[h]);
				reader = new FileReader(treeFile);
				nexusImporter = new NexusImporter(reader);
				nexsusTrees = nexusImporter.importTrees();
				System.out.println("loaded "+nexsusTrees.size()+" trees");

				System.out.print("Keeping tail... ");		
				nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
				for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
					nexsusTreeTail.add(nexsusTrees.get(i));
				}
				System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");			
				for (int i=0; i<numTestTrees;i++) {				
					TreeWithLocations testTree = new TreeWithLocations(createModel,(SimpleRootedTree) nexsusTreeTail.get(Random.nextInt(nexsusTreeTail.size())));
					testTree.fillRandomTraits();
					testTree.removeInternalLocations();
					trees.get(h).add(testTree);				
				}
			}

			System.out.println("Generated "+trees.get(0).size()*trees.size()+" trees with model generated random tip annotations and input tree topology");

		}
		break;

		case TEST_MODEL_DEGENERACY:{

			numLocations=numTestLocations;
			// TODO: add more tests + test files ....
			// TODO: add tests for states...
			System.out.print("Building degenerate test models... ");

			// Generate test data and trees

			// For constant model...
			Q = makeRandomMigrationMatrix(numLocations,2); 

			// For two seasonal model...
			QW = myMatrixCopy(Q);
			QS = myMatrixCopy(Q);

			// For sinusoidal model...
			rates = myMatrixCopy(Q);
			amps = makeRandomMigrationMatrix(numLocations,0);
			phases = makeRandomMigrationMatrix(numLocations,1);

			// For two constant seasons model...
			double phase = 0.3; double length = 0.5;

			switch (config.migrationSeasonality) {
			case NONE:	
				createModel = new ConstantMigrationBaseModel(Q); 
				break;
			case TWO_CONSTANT_SEASONS: 
				createModel = new TwoSeasonMigrationBaseModel(QW,QS,phase,phase+length);
				break;
			case SINUSOIDAL:
				createModel = new SinusoidialSeasonalMigrationBaseModel(rates,amps,phases);
				break;
			default: 
				System.err.println("Migration Seasonality: "+config.migrationSeasonality+" not implemented for this configuration!!!");
				System.exit(-1);
			}

			System.out.print(" done!\nLoading trees... ");	
			
			trees.add(new ArrayList<LikelihoodTree>());
			
			treeFile = new File(config.treeFilenames[0]);
			reader = new FileReader(treeFile);
			nexusImporter = new NexusImporter(reader);
			List<Taxon> taxa = nexusImporter.parseTaxaBlock();
			HashMap<String,Integer> taxaIndices = new HashMap<String,Integer>();			
			for (int i=0;i<taxa.size();i++) {
				taxaIndices.put(taxa.get(i).getName(), i);
			}
			nexsusTrees = nexusImporter.importTrees();
			System.out.println("loaded "+nexsusTrees.size()+" trees");

			System.out.print("Keeping tail... ");		
			nexsusTreeTail = new ArrayList<jebl.evolution.trees.Tree>();
			for (int i=Math.max(0,nexsusTrees.size()-config.numTreesFromTail);i<nexsusTrees.size();i++) {
				nexsusTreeTail.add(nexsusTrees.get(i));
			}
			System.out.println(" keeping last "+nexsusTreeTail.size()+ " trees");			

			// Convert trees to internal tree representation
			if (config.locationFilenames[0]!=null) {
				System.out.print("Loading traits... ");
				AttributeLoader attributeLoader= new SimpleAttributeLoader(config.locationFilenames[0], config.stateFilename);
				// TODO: think about this...
				HashMap<String,Integer> locationMap = (HashMap<String,Integer>) attributeLoader.getAttributes().get("locations");
				HashMap<String,Double> stateMap = (HashMap<String,Double>) attributeLoader.getAttributes().get("states");
				System.out.println("loaded "+locationMap.size()+" taxon traits");

				System.out.print("Reparsing trees... ");
				if (stateMap==null) {
					for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
						trees.get(0).add(new TreeWithLocations((SimpleRootedTree) tree,taxaIndices,locationMap,numLocations));
					}
				}
				else {
					// TODO: this...
				}
				System.out.println(" reparsed "+trees.size()+" trees");
			}
			else {
				// TODO: add load states from trees...
				System.out.print("Reparsing trees... ");
				for (jebl.evolution.trees.Tree tree : nexsusTreeTail) {
					trees.get(0).add(new TreeWithLocations((SimpleRootedTree) tree,taxaIndices,config.locationAttributeNameInTree, numLocations));
				}		
				System.out.println(" reparsed "+trees.size()+" trees");
			}				

			testModels = new ArrayList<MigrationBaseModel>();

			testModels.add(new ConstantMigrationBaseModel(Q));
			testModels.add(new TwoSeasonMigrationBaseModel(QW,QS,phase, phase+length));
			
			generalMigrationFunction = new DoubleFunction[Q.length][Q[0].length]; 
			for (int i=0;i<Q.length;i++) {
				for (int j=0;j<Q[0].length;j++) {
					generalMigrationFunction[i][j]=cern.jet.math.Functions.constant(Q[i][j]);
				}
			}
			
			rootFreqFunction = new DoubleFunction[Q.length]; 
			for (int i=0;i<Q.length;i++) {
				rootFreqFunction[i]=cern.jet.math.Functions.constant(1.0/Q.length);
				for (int j=0;j<Q[0].length;j++) {
					generalMigrationFunction[i][j]=cern.jet.math.Functions.constant(Q[i][j]);
				}
			}
			
			testModels.add(new GeneralSeasonalMigrationBaseModel(generalMigrationFunction, rootFreqFunction,25));
		
			System.out.println(" generated "+testModels.size()+" test models");

		}
		break;

		}

		// Creating test file 
		System.out.print("Writing test model to out.test ...");
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		testStream.print(createModel.parse());
		testStream.println();
		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();
		System.out.println("done");

	}

	private double[][] myMatrixCopy(double[][] q) {
		double [][] returnValue = new double[q.length][];
		for(int i = 0; i < q.length; i++) {
			returnValue[i] = q[i].clone();
		}				
		return returnValue;
	}

	public static double[][] makeRandomMigrationMatrix(int size, double scale) {
		// For test purposes...
		double[][] returnValue = new double[size][size];
		for (int i=0;i<size;i++) {
			double rowSum=0;
			for (int j=0;j<size;j++) {
				if (i!=j) {
					returnValue[i][j]=Math.random()*scale;
					rowSum+=returnValue[i][j];
				}
			}
			returnValue[i][i]=-rowSum;
		}
		return returnValue;
	}

	private double[][] disturbMigrationMatrix(double[][] migrationMatrix, double disturbanceMagnitude, double max) {
		// For test purposes...
		double[][] returnValue = myMatrixCopy(migrationMatrix);
		for (int i=0;i<migrationMatrix.length;i++) {
			double rowSum=0;
			for (int j=0;j<migrationMatrix.length;j++) {
				if (i!=j) {
					do {
						returnValue[i][j]+=(Math.random()-0.5)*disturbanceMagnitude;
					} while ((returnValue[i][j]<0) || (returnValue[i][j]>max));
					rowSum+=returnValue[i][j];
				}
			}
			returnValue[i][i]=-rowSum;
		}
		return returnValue;
	}

	@Override
	public List<ArrayList<LikelihoodTree>> getTrees() {
		return trees;
	}

	@Override
	public int getNumLocations() {
		return numLocations;
	}






}
