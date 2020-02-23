/*
 * maindesign.cpp
 *
 *  Created on: 24-nov-2018
 *      Author: M. El-Kebir
 */

#include <lemon/arg_parser.h>
#include "utils.h"
#include "inputinstance.h"
#include "setcover.h"
#include "setcoverilp.h"
#include <vector>
#include <fstream>

void getSetCoverInput(const InputInstance& input,
                      int treeIndex,
                      int sampleIndex,
                      CloneTree& tree,
                      CloneTreeVector& otherTrees,
                      std::map<std::string, double>& freqMap)
{
  tree = input.getTrees()[treeIndex];
  for (int i = 0; i < input.getNrTrees(); ++i)
  {
    if (i != treeIndex)
      otherTrees.push_back(input.getTrees()[i]);
  }
  
  for (NodeIt v(tree.tree()); v != lemon::INVALID; ++v)
  {
    int idx = input.getFrequencies().characterToIndex(tree.label(v));
    freqMap[tree.label(v)] = input.getFrequencies().min(sampleIndex, idx);
  }
}

int main(int argc, char** argv)
{
  int treeIndex = 0;
  int sampleIndex = 0;
  
  lemon::ArgParser ap(argc, argv);
  ap.other("input", "Input");
  ap.refOption("i", "Tree index", treeIndex);
  ap.refOption("p", "Sample index", sampleIndex);
  ap.run();
  
  if (ap.files().empty())
  {
    std::cerr << "Error: no input specified" << std::endl;
    return 1;
  }
  
  std::string inputFilename(ap.files()[0]);
  
  InputInstance input;
  try
  {
    std::ifstream in(inputFilename.c_str());
    if (!in.good())
    {
      std::cerr << "Error: could not open '" << inputFilename << "' for reading" << std::endl;
      return 1;
    }
    
    in >> input;
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  
  // let's print mutations of first tree
  if (!(0 <= treeIndex && treeIndex < input.getNrTrees()))
  {
    std::cerr << "Error: invalid tree index" << std::endl;
    return 1;
  }
  
  if (!(0 <= sampleIndex && sampleIndex < input.getFrequencies().getNrSamples()))
  {
    std::cerr << "Error: invalid sample index" << std::endl;
    return 1;
  }
  
  //  const FrequencyMatrix& F = input.getFrequencies();
  //  const CloneTree& T = input.getTrees()[treeIndex];
  //  for (NodeIt v(T.tree()); v != lemon::INVALID; ++v)
  //  {
  //    for (const std::string& mutation : T.nodeToMutations(v))
  //    {
  //      std::cout << mutation << " ";
  //    }
  //    const int i = F.characterToIndex(T.label(v));
  //    std::cout << F.min(sampleIndex, i) << std::endl;
  //  }
  
  CloneTree T;
  CloneTreeVector otherTrees;
  
  std::map<std::string, double> freqMap;
  getSetCoverInput(input, treeIndex, sampleIndex,
                   T, otherTrees, freqMap);
  
  SetCover sc(T, otherTrees, freqMap);
  sc.writeDOT(std::cout);
  
  
  std::vector <std::vector<int>> solutions;
  int first_val;
  
  do{
    std::cout << "solving ip" << std::endl;
    SetCoverIlp sci(T, otherTrees, freqMap, solutions);
    std::vector<int> new_solution = sci.solve();
    first_val = new_solution[0];
    if(new_solution[0] >= 0){
      solutions.push_back(new_solution);
    }
    
  }while(first_val > 0);

  SetCoverIlp allSol(T, otherTrees, freqMap, solutions);
  allSol.printSolutions(inputFilename, "distFeats", treeIndex,solutions );
  
  
  
   return 0;
}
