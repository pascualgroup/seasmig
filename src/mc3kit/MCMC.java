/***
  This file is part of mc3kit.
  
  Copyright (C) 2013 Edward B. Baskerville

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
***/

package mc3kit;

import java.util.*;
import java.util.concurrent.*;
import java.util.logging.*;
import java.io.*;

import mc3kit.util.JsonsLogFormatter;

import static java.lang.String.format;

import cern.jet.random.engine.*;

@SuppressWarnings("serial")
public class MCMC implements Serializable {
  
  /*** SETTABLE FIELDS ***/
  
  ModelFactory modelFactory;
  int chainCount = 1;

  private HeatFunction heatFunction;
  private double[] priorHeatExponents;
  private double[] likelihoodHeatExponents;

  private Long randomSeed = null;
  private List<Step> steps;
  
  /*** STATE ***/
  
  boolean initialized;
  long iterationCount;

  RandomSeedGenerator seedGen;
  Chain[] chains;

  TaskManager[][] taskManagers;

  Level logLevel;
  String logFilename;
  boolean logAllChains;
  
  transient Logger _logger;
  transient Map<String, Logger> _loggers;
  
  /*** THREAD POOL MANAGEMENT ***/

  transient ThreadPoolExecutor threadPool;
  transient ExecutorCompletionService<Object> completionService;
  transient TerminationManager terminationManager;
  transient long terminationCount;
  
  public MCMC() {
    steps = new ArrayList<Step>();
    
    initializeLoggers();
  }
  
  private void initializeLoggers() {
    _logger = Logger.getLogger("");
    for(Handler handler : _logger.getHandlers()) {
      _logger.removeHandler(handler);
    }
    _loggers = new HashMap<String, Logger>();
  }
  
  private Logger getLogger() {
    return getLogger("mc3kit.MCMC");
  }
  
  public Logger getLogger(String name) {
    synchronized(_loggers) {
      String[] pieces = name.split("\\.");
      if(!pieces[0].equals("mc3kit")) {
        throw new IllegalArgumentException("name must begin with mc3kit.");
      }
      String subName = "mc3kit";
      for(int i = 1; i < pieces.length; i++) {
        subName = subName + "." + pieces[i];
        if(!_loggers.containsKey(subName)) {
          Logger logger = Logger.getLogger(subName);
          logger.setLevel(logLevel);
          _loggers.put(subName, logger);
        }
      }
      
      assert _loggers.containsKey(name);
      
      return _loggers.get(name);
    }
  }
  
  private void initialize() throws Throwable {
    if(initialized) {
      return;
    }
    
    initializeThreadPool();
    
    if(heatFunction == null) {
      heatFunction = new ConstantHeatFunction();
    }
     
    priorHeatExponents = heatFunction.getPriorHeatExponents(chainCount);
    likelihoodHeatExponents = heatFunction.getLikelihoodHeatExponents(chainCount);

    for (int i = 0; i < chainCount; i++) {
      getLogger().info(format("chain %d: prior heat exp. %f, like heat exp %f", i,
          priorHeatExponents[i], likelihoodHeatExponents[i]));
    }

    // Initialize a top-level RNG to generate start seed-table
    // start locations for each replicate
    Random random = new Random();
    if (randomSeed == null)
      randomSeed = random.nextLong();
    random.setSeed(randomSeed);
    int randomSeedStartLocation = random.nextInt(Integer.MAX_VALUE);
    seedGen = new RandomSeedGenerator(randomSeedStartLocation, 0);

    initializeChains();
    initializeSteps();
    
    initialized = true;
  }
  
  private void readObject(ObjectInputStream in) throws IOException, ClassNotFoundException, SecurityException, MC3KitException {
    in.defaultReadObject();
    
    initializeLoggers();
    
    if(logLevel != null) {
      setLogLevelPrivate(logLevel);
    }
    
    if(logFilename != null) {
      setLogFilenamePrivate(logFilename);
    }
    
    if(initialized) {
      initializeThreadPool();
    }
  }
  
  private void initializeThreadPool() {
    threadPool = new ThreadPoolExecutor(chainCount, 2 * chainCount, 60,
        TimeUnit.SECONDS, new SynchronousQueue<Runnable>());
    completionService = new ExecutorCompletionService<Object>(threadPool);
  }

  private void initializeChains() throws Throwable {
    if(chainCount < 1) {
      throw new MC3KitException("There must be at least one chain.");
    }

    if(modelFactory == null) {
      throw new MC3KitException("modelFactory must be set."); 
    }
    
    // Make chains
    chains = new Chain[chainCount];
    for(int i = 0; i < chainCount; i++) {
      RandomEngine rng = makeRandomEngine();
      chains[i] = new Chain(this, i, chainCount, priorHeatExponents[i],
          likelihoodHeatExponents[i], rng);
    }
    
    // Use thread pool to construct models for different chains
    for(int i = 0; i < chainCount; i++) {
      final int chainId = i;
      safeSubmit(new Callable<Object>() {
        @Override
        public Object call() throws Exception {
          
          Model model = modelFactory.createModel(chains[chainId]);
          assert model != null;
          chains[chainId].setModel(model);
          model.setChain(chains[chainId]);
          return this;
        }
      });
    }
    
    // Wait for completion
    for(int i = 0; i < chainCount; i++) {
      Future<Object> completedTask = completionService.take();
      completedTask.get();
    }
    
    getLogger().info("Chains created.");
  }

  private void initializeSteps() throws MC3KitException {
    Map<Step, List<Task>> taskMap = new HashMap<Step, List<Task>>();

    taskManagers = new TaskManager[steps.size()][];
    TaskManager[][] taskManagersByChain = new TaskManager[steps.size()][chainCount];
    for (int i = 0; i < steps.size(); i++) {
      // Get existing tasks for the Step object, or create new ones
      Step step = steps.get(i);
      List<Task> tasks;
      if (taskMap.containsKey(step)) {
        tasks = taskMap.get(step);
      }
      else {
        tasks = step.makeTasks(chainCount);
        taskMap.put(step, tasks);
      }

      // Create TaskManager objects for each task at each step.
      // There will be multiple TaskManagers for a task if a step appears
      // multiple times.
      taskManagers[i] = new TaskManager[tasks.size()];
      for (int j = 0; j < tasks.size(); j++) {
        Task task = tasks.get(j);
        int[] chainIds = task.getChainIds();

        taskManagers[i][j] = new TaskManager(task, i == 0);

        for (int chainId : chainIds) {
          if (taskManagersByChain[i][chainId] != null) {
            throw new MC3KitException(
              String.format("Two tasks handling chain %d in step %d.", chainId, i)
            );
          }
          taskManagersByChain[i][chainId] = taskManagers[i][j];
        }
      }
    }

    // Make sure that every chain is covered in the first step
    for (int chainId = 0; chainId < chainCount; chainId++) {
      if (taskManagersByChain[0][chainId] == null)
        throw new MC3KitException(
          String.format("Chain %d is not covered by first step.", chainId)
        );
    }

    // Identify forward dependencies for each step method
    for (int i = 0; i < steps.size(); i++) {
      for (TaskManager taskManager : taskManagers[i]) {
        Chain[] handledChains = taskManager.handledChains;
        TaskManager[] nextTaskManagers = new TaskManager[handledChains.length];
        for (int j = 0; j < handledChains.length; j++) {
          int chainId = handledChains[j].getChainId();
          
          // Find the next step method sequentially that uses this chain,
          // looping back if necessary.
          for (int k = i + 1; k <= steps.size() + i; k++) {
            nextTaskManagers[j] = taskManagersByChain[k % steps.size()][chainId];
            if (nextTaskManagers[j] != null)
              break;
          }
        }

        taskManager.nextTaskManagers = nextTaskManagers;
      }
    }
    getLogger().info("Steps set up.");
  }
  
  public static MCMC loadFromFile(String filename) throws ClassNotFoundException, IOException {
    ObjectInput input = new ObjectInputStream(new BufferedInputStream(new FileInputStream(filename)));
    MCMC mcmc = (MCMC)input.readObject();
    input.close();
    
    return mcmc;
  }
  
  /**
   * @throws IOException 
   * @throws FileNotFoundException 
   * 
   */
  public synchronized void writeToFile(String filename) throws FileNotFoundException, IOException {
    File file = new File(filename);
    File tmpFile = new File(filename + ".partial");
    
    ObjectOutputStream stream = new ObjectOutputStream(new FileOutputStream(tmpFile));
    stream.writeObject(this);
    stream.close();
    
    if(file.exists()) {
      file.delete();
    }
    if(!tmpFile.renameTo(file)) {
      getLogger().severe("COULD NOT MOVE CHECKPOINT FILE");
    }
    assert tmpFile.getName().endsWith(".partial");
    if(tmpFile.exists()) {
      tmpFile.delete();
    }
  }
  
  /**
   * Run MCMC for a single step.
   * @throws MC3KitException 
   */
  public synchronized void step() throws Throwable {
    runFor(1);
  }
  
  /**
   * Run MCMC for multiple steps.
   * @param iterationCount
   * @throws MC3KitException 
   */
  public synchronized void runFor(long runForCount) throws Throwable {
    if(runForCount == 0)
      return;
    
    terminationCount = iterationCount + runForCount;
    run();
  }
  
  /**
   * Shut down the thread pool so the program can exit.
   */
  public synchronized void shutdown() throws Throwable {
    if(initialized) {
      threadPool.shutdown();
    }
  }
  
  public synchronized long getIterationCount() {
    return iterationCount;
  }
  
  private synchronized void run() throws Throwable {
    assert(terminationCount > iterationCount);
    initialize();
    
    getLogger().fine(format("Running until %d", terminationCount));
    
    terminationManager = new TerminationManager();
    
    getLogger().fine("Starting run.");

    try {
      for (TaskManager taskManager : taskManagers[0]) {
        assert(taskManager.iterationCount == iterationCount);
        safeSubmit(taskManager);
      }

      boolean done = false;
      while(!done) {
        Future<Object> completedTask = completionService.take();
        Object result = completedTask.get();
        if(result == terminationManager) {
          getLogger().fine("Got termination task.");
          done = true;
          assert(iterationCount == terminationCount);
        }
        else {
          getLogger().finer("Got non-termination task.");
        }
      }
    }
    catch (ExecutionException e) {
      throw e.getCause();
    }
    
    terminationManager = null;
  }
  
  Chain getChain(int chainId) {
    if (chainId < 0 || chainId >= chainCount)
      return null;

    return chains[chainId];
  }
  
  public synchronized ModelFactory getModelFactory() { 
    return modelFactory;
  }
  
  public synchronized Model getModel() {
    return getChain(0).getModel();
  }
  
  public synchronized Model getModel(int chainId) {
    return getChain(chainId).getModel();
  }

  synchronized RandomEngine makeRandomEngine() {
    return new MersenneTwister(seedGen.nextSeed());
  }
  
  class TerminationManager implements Callable<Object>, Serializable {
    BitSet completedChains;
    
    public TerminationManager() {
      completedChains = new BitSet();
    }
    
    public void tallyComplete(int[] chainIds) {
      boolean complete = false;
      
      synchronized(completedChains) {
        for(int chainId : chainIds) {
          assert(!completedChains.get(chainId));
          
          completedChains.set(chainId);
        }
        
        if(completedChains.cardinality() == chainCount) {
          complete = true;
        }
      }
      
      if(complete) {
        safeSubmit(this);
      }
    }

    @Override
    public Object call() throws Exception {
      iterationCount = terminationCount;
      return this;
    }
  }
  
  class TaskManager implements Callable<Object>, Serializable {
    long iterationCount;
    int completedChainCount;
    Task task;
    boolean isTerminationHandler;
    Chain[] handledChains;

    TaskManager[] nextTaskManagers;

    public TaskManager(Task task, boolean isTerminationHandler) {
      this.task = task;
      this.isTerminationHandler = isTerminationHandler;

      int[] chainIds = task.getChainIds();
      handledChains = new Chain[chainIds.length];
      for (int i = 0; i < chainIds.length; i++) {
        handledChains[i] = chains[chainIds[i]];
      }
    }

    public void tallyComplete() {
      boolean allComplete = false;
      synchronized (this) {
        completedChainCount++;

        if (completedChainCount == handledChains.length) {
          allComplete = true;
          completedChainCount = 0;
        }
      }
      if (allComplete) {
        iterationCount++;
        
        if(isTerminationHandler) {
          for(Chain chain : handledChains) {
            chain.iterationCount++;
          }
        }
        
        if(isTerminationHandler && iterationCount == terminationCount) {
          terminationManager.tallyComplete(task.getChainIds());
        }
        else {
          safeSubmit(this);
        }
      }
    }

    @Override
    public TaskManager call() throws Exception {
      try {
        task.step(handledChains);
      }
      catch (Exception e) {
        throw new MC3KitException("Exception thrown from task", e);
      }

      // Let the next step for this chain know we're complete;
      // it will start itself up once all its chains are ready
      for (TaskManager manager : nextTaskManagers)
        manager.tallyComplete();

      return this;
    }

    @Override
    public String toString() {
      return task.toString();
    }
  }

  private void throwIfInitialized() throws MC3KitException {
    if (initialized)
      throw new MC3KitException("Already initialized");
  }
  
  private Future<Object> safeSubmit(Callable<Object> task) {

    boolean submitted = false;
    Future<Object> result = null;
    int sleepMillis = 1;
    while (!submitted) {
      try {
        result = completionService.submit(task);
        submitted = true;
      }
      catch (RejectedExecutionException e) {
        try {
          Thread.sleep(sleepMillis);
        }
        catch (InterruptedException e1) {
        }
        sleepMillis *= 2;
      }
    }
    return result;
  }
  
  /*** ACCESSORS ***/

  public synchronized void setModelFactory(ModelFactory modelFactory) throws MC3KitException {
    throwIfInitialized();
    this.modelFactory = modelFactory;
  }

  public synchronized void setChainCount(int chainCount) throws MC3KitException {
    throwIfInitialized();
    this.chainCount = chainCount;
  }
  
  public synchronized void addStep(Step step) throws MC3KitException {
    throwIfInitialized();
    steps.add(step);
  }

  public synchronized void setHeatFunction(double heatPower) throws MC3KitException {
    throwIfInitialized();
    setHeatFunction(new PowerHeatFunction(heatPower, 0.0));
  }
  
  public synchronized void setHeatFunction(double heatPower, double minHeatExponent) throws MC3KitException {
    throwIfInitialized();
    setHeatFunction(new PowerHeatFunction(heatPower, minHeatExponent));
  }
  
  public synchronized void setHeatFunction(double likeHeatPower, double minLikeHeatExponent, double priorHeatPower, double minPriorHeatExponent)
      throws MC3KitException {
    throwIfInitialized();
    setHeatFunction(new PowerHeatFunction(likeHeatPower, minLikeHeatExponent, priorHeatPower, minPriorHeatExponent));
  }
  
  public synchronized void setHeatFunction(HeatFunction heatFunction) throws MC3KitException {
    throwIfInitialized();
    this.heatFunction = heatFunction;
  }

  public synchronized void setRandomSeed(Long randomSeed) throws MC3KitException {
    throwIfInitialized();
    this.randomSeed = randomSeed;
  }
  
  public synchronized void setLogLevel(Level level) throws MC3KitException {
    throwIfInitialized();
    setLogLevelPrivate(level);
  }
  
  private void setLogLevelPrivate(Level level) throws MC3KitException {
    this.logLevel = level;
    _logger.setLevel(level);
    for(Handler handler : _logger.getHandlers()) {
      handler.setLevel(level);
    }
  }
  
  public synchronized void setLogFilename(String filename) throws SecurityException, IOException, MC3KitException {
    throwIfInitialized();
    setLogFilenamePrivate(filename);
  }
  
  private void setLogFilenamePrivate(String filename) throws SecurityException, IOException, MC3KitException {
    this.logFilename = filename;
    
    for(Handler handler : _logger.getHandlers()) {
      _logger.removeHandler(handler);
    }
    
    JsonsLogFormatter formatter = new JsonsLogFormatter();
    
    if(filename == null || filename.equals("-") || filename.equals("")) {
      ConsoleHandler handler = new ConsoleHandler();
      handler.setFormatter(formatter);
      handler.setLevel(logLevel);
      _logger.addHandler(handler);
    }
    else {
      FileHandler handler = new FileHandler(filename);
      handler.setFormatter(formatter);
      handler.setLevel(logLevel);
      _logger.addHandler(handler);
    }
  }
  
  public synchronized void setLogAllChains(boolean logAllChains) throws MC3KitException {
    throwIfInitialized();
    this.logAllChains = logAllChains;
  }
}
