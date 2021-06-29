import java.util.*;
import java.lang.Math;

public class SearchAlgorithm {

  // Your search algorithm should return a solution in the form of a valid
  // schedule before the deadline given (deadline is given by system time in ms)
  public Schedule solve(SchedulingProblem problem, long deadline) {

    // get an empty solution to start from
    Schedule solution = problem.getEmptySchedule();

    // YOUR CODE HERE

    return solution;
  }

  // This is a very naive baseline scheduling strategy
  // It should be easily beaten by any reasonable strategy
  public Schedule naiveBaseline(SchedulingProblem problem, long deadline) {

    // get an empty solution to start from
    Schedule solution = problem.getEmptySchedule();

    for (int i = 0; i < problem.courses.size(); i++) {
      Course c = problem.courses.get(i);
      boolean scheduled = false;
      for (int j = 0; j < c.timeSlotValues.length; j++) {
        if (scheduled) break;
        if (c.timeSlotValues[j] > 0) {
          for (int k = 0; k < problem.rooms.size(); k++) {
            if (solution.schedule[k][j] < 0) {
              solution.schedule[k][j] = i;
              scheduled = true;
              break;
            }
          }
        }
      }
    }

    return solution;
  }
// needed to create new population for genetic algo
  public Schedule naiveBaseline(SchedulingProblem problem) {

    // get an empty solution to start from
    Schedule solution = problem.getEmptySchedule();

    for (int i = 0; i < problem.courses.size(); i++) {
      Course c = problem.courses.get(i);
      boolean scheduled = false;
      for (int j = 0; j < c.timeSlotValues.length; j++) {
        if (scheduled) break;
        if (c.timeSlotValues[j] > 0) {
          for (int k = 0; k < problem.rooms.size(); k++) {
            if (solution.schedule[k][j] < 0) {
              solution.schedule[k][j] = i;
              scheduled = true;
              break;
            }
          }
        }
      }
    }

    return solution;
  }
  ///////////////////////
  // BackTracking///////
  /////////////////////
  // Calls recursive backtracking
  public Schedule BackTracking(SchedulingProblem problem, long deadline)
  {
    Schedule solution = problem.getEmptySchedule();
    //Solution, Problem and size
    return recursiveBackTracking(solution,problem,0);
  }
  public Schedule recursiveBackTracking(Schedule sol, SchedulingProblem problem, int size) {
    //Base case
    if(problem.courses.size() == size){
      return sol;
    }
    //Temp course created
    Course tmp = problem.courses.get(size);
    //loop till we have gone through the slotvalues
    for (int row = 0; row < tmp.timeSlotValues.length; row++){
      if (tmp.timeSlotValues[row] > -1) {
        //Loop through all the rooms
        for (int col = 0; col < problem.rooms.size(); col++){
          //Check if the schedule is open and place a class
          if (sol.schedule[col][row] < 0) {
            sol.schedule[col][row] = size;

            //Recursively look at the rest of the problem
            Schedule tmp_solution = recursiveBackTracking(sol,problem, size = size + 1);
            if(tmp_solution !=  null){
              return sol;
            }
            sol.schedule[col][row] = -1;
          }

        }
      }
    }
    return null;
  }
  /////////////////////
  ////Genetic Algo////
  ////////////////////
  public Schedule solveGenetic(SchedulingProblem problem, long deadline){
    // Create population space of Scheudles to be used for our individuals
    Schedule[] population = new Schedule[8];
    //loop for fittest individual and store it
    for(int i = 0; i <population.length; i++){
      //Generations here
      SchedulingProblem newProb = problem;
      newProb.createRandomInstance(problem.buildings.size(),problem.rooms.size(),problem.courses.size());
      population[i] = (naiveBaseline(newProb));
    }
    // After populating our population, we select the fittest schedule and return it
    return getFittest(evolvePopulation(population,problem),problem);

  }
  // Generate child schedules from parent
  public Schedule crossover(Schedule first, Schedule second,SchedulingProblem problem) {

    // Iterate and get the smallest schedule that is not out of range
    Schedule inboundSchedule;
    if(first.schedule.length > second.schedule.length){
      inboundSchedule = second;
    }
    else{
      inboundSchedule = first;
    }
    // Create new solution with the dimensions of the small schedule
    Schedule newSol = new Schedule(inboundSchedule.schedule.length,inboundSchedule.schedule[0].length);
    //Swap genes based on random value threshold crossover and swap gene from parent if within the threshold
    for (int row = 0; row < inboundSchedule.schedule.length; row++) {
      for(int col = 0; col < newSol.schedule[row].length ;col++){
        if (Math.random() <= 0.32) {
          //Return index of gene/time slot from the individual
          newSol.schedule[row][col] = first.schedule[row][col];
        } else {
          //Return index of gene/time slot from the individual
          newSol.schedule[row][col] = second.schedule[row][col];
        }
      }
    }
    return newSol;
  }

  public void mutateSchedule(Schedule individual) {
    //Create new random value to be used for solution
    Random random = new Random();
    //Iterate row col of individual and assign new value
    for (int row = 0; row < individual.schedule.length; row++) {
      for(int col = 0; col < 10; col++){
        if (Math.random() <= 0.14) {
          int newGene = random.nextInt(individual.schedule.length);
          // Sets a gene in a given individual
          individual.schedule[row][col] = newGene;
        }
      }
    }
  }

  // This method randomly generates a fit individual from a population
  public Schedule randomSelection(Schedule[] population, SchedulingProblem problem) {
    // Generate the population
    Schedule[] newPopulation = new Schedule[population.length];
    // Populate the population with random individuals from another population
    for (int index = 0; index < newPopulation.length; index++) {
      int randomIndex = (int) (Math.random() * newPopulation.length);
      newPopulation[index] = population[randomIndex];
    }
    // We then return the fittest of the population that we generated
    Schedule fittestIndividual = getFittest(newPopulation,problem);
    return fittestIndividual;
  }
  //Take population and problem and find fittest individual
  public Schedule getFittest(Schedule[] population, SchedulingProblem problem){
    //Create fit individual schedule to be compared too
    Schedule fittestIndividual = population[0];
    //Iterate over entire population
    for(int i = 0; i < population.length; i++){
      //Compare to get fittest individual
      //When evaluating we will get errors about schedule dimensions being invalid can comment errow printing out in schedulingProblem
      if(problem.evaluateSchedule(population[i]) > problem.evaluateSchedule(fittestIndividual)){
        fittestIndividual = population[i];
      }
    }
    return fittestIndividual;
  }


  //Populations created randomly and used them to crossover genes into a child
  public Schedule[] evolvePopulation(Schedule[] population, SchedulingProblem problem) {
    // Generate a new population to store the size of the current population
    Schedule[] newPop = new Schedule[population.length];
    // Randomly select two of the fittest parents and crossover their features into a child
    for (int index = 0; index < population.length; index++) {
      //Iterations here
      Schedule firstInd= randomSelection(population, problem);
      Schedule secondInd = randomSelection(population, problem);
      Schedule newIndividual = crossover(firstInd, secondInd,problem);
      //Put child into new generation
      newPop[index] = newIndividual;
    }
    for (int index = 0; index < newPop.length; index++) {
      // Mutate random features in a schedule
      mutateSchedule(population[index]);
    }

    return newPop;
  }





}


