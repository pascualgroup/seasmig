package seasmig;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.GregorianCalendar;
import jebl.evolution.io.ImportException;
import org.junit.Test;

import seasmig.data.DataForTests;
import seasmig.data.DataForTests.TestType;
import seasmig.treelikelihood.LikelihoodTree;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class TestOther {
	
	// TODO: Organize this....
	
	@Test
	public void testModelDegeneracy() throws IOException, ImportException {
		// Load config
		System.out.print("Loading config file... ");
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
		Config config = null;
		try {
			config = gson.fromJson(new FileReader("config.json"), Config.class);
			System.out.println(" done");
		}
		catch(Throwable e)	{
			config=new Config();
			System.out.println("config.json file not found, using default config. See out.config.json for details");			
		}			

		System.out.print("Writing full config options to out.config...");
		config.outputToFile("out.config",gson);
		System.out.println(" done");

		// Load data files and prepare data....			
		Data data = new DataForTests(config,TestType.TEST_MODEL_DEGENERACY,3,3,10);

		// Creating test file 
		File testFile = new File("out.test");
		testFile.delete();
		testFile.createNewFile();
		PrintStream testStream = new PrintStream(testFile);
		System.out.println("Calculating tree likelihood using degenerate models:");				
		double[] results = new double[((DataForTests) data).testModels.size()];

		for (int i=0;i<((DataForTests) data).testModels.size();i++) {
			System.out.println("SEASONALITY "+((DataForTests) data).testModels.get(i).getModelName());						
			long startTime= System.currentTimeMillis();	
			double testLikelihood = 0;
			for (LikelihoodTree tree : data.getTrees().get(0)) {
				System.out.print(".");
				LikelihoodTree workingCopy = tree.copy();
				workingCopy.setMigrationModel(((DataForTests) data).testModels.get(i));
				testLikelihood+=workingCopy.logLikelihood();
			}
			testLikelihood=testLikelihood/data.getTrees().size();
			System.out.println(testLikelihood);
			results[i]=testLikelihood;
			long duration= System.currentTimeMillis()-startTime;
			System.out.println("duration: "+duration+"ms");
		}

		testStream.print(",\""+(new GregorianCalendar()).getTime()+"\"}");
		testStream.close();

		for (int i=1;i<results.length;i++) {
			assertEquals(results[i],results[i-1], 1E-3);			
		}
		System.exit(0);
	}	
}
