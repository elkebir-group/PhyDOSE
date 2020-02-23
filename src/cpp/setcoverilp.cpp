/*
 * setcoverilp.cpp
 *
 *  Created on: 5-dec-2019
 *      Author: M. El-Kebir
 */

#include "setcoverilp.h"
#include <vector>
#include "utils.h"

SetCoverIlp::SetCoverIlp(const CloneTree& T,
                         const CloneTreeVector& scriptT,
                         const std::map<std::string, double>& freqMap,
                         const std::vector<std::vector<int>> solution_set)
                         
  : SetCover(T, scriptT, freqMap)
  , _env()
  , _model(_env)
  , _cplex(_model)
{
  init(solution_set);
}

void SetCoverIlp::init(std::vector<std::vector<int>> solution_set)
{
  char buf[1024];
  
  const int nrFeatures = _featureVector.size();

  
  
  
  // Initialize variables
  _x = IloBoolVarArray(_env, nrFeatures);
  for (int i = 0; i < nrFeatures; ++i)
  {
    snprintf(buf, 1024, "x:%d", i);
    _x[i] = IloBoolVar(_env, buf);
  }
  
  //_z = IloNumVar(_env, -IloInfinity, 0, "z");
  
  // Initialize constraints
  IloExpr sum(_env);

  IloExpr obj(_env);
  
  // // there are _k features in the cover
  // for (int i = 0; i < nrFeatures; ++i)
  // {
  //   sum += _x[i];
  // }
  // _model.add(sum == _k);
  // sum.clear();
  
  lemon::DynArcLookUp<Graph> lookUp(_G);
  for (int j = 0; j < _nodesTree.size(); ++j)
  {
    Graph::Node v_j = _nodesTree[j];
    for (int i = 0; i < nrFeatures; ++i)
    {
      Graph::Node v_i = _nodesFeatures[i];
      if (lookUp(v_i, v_j) != lemon::INVALID)
      {
        sum += _x[i];
      }
    }
    // every tree is covered by at least one feature
    _model.add(sum >= 1);
    sum.clear();
  }

  for( std::vector<int> sol: solution_set){
    int card = sol.size();
    for( int s: sol){
      sum += _x[s];
    }
    _model.add(sum < card);
    sum.clear();
  }
 
  //create the objective function
  for (int i = 0; i < nrFeatures; ++i)
  {
    //assert(_freqMap.count(_mutationVector[i]) == 1);
    //double freq = log(_freqMap.find(_mutationVector[i])->second);
    obj += _x[i];
  }
  
  _model.add(IloMinimize(_env, obj));
}

std::vector<int> SetCoverIlp::solve()
{
  std::vector<int> new_sol;
  _cplex.setParam(IloCplex::Threads, 1);
  if (!_cplex.solve())
  {
    new_sol.push_back(-1);
    return new_sol;
  }
  
  _objValue = _cplex.getBestObjValue();
  std::cout << "Obj value: " << _objValue << std::endl;
 
  const int nrFeatures = _featureVector.size();
  for (int i = 0; i < nrFeatures; ++i)
  {
    if (_cplex.getValue(_x[i]) > 0.4)
    {
      
      new_sol.push_back(i);
      std::cout << "Featurette:";
      for(std::string s: _featureVector[i]){
         std::cout  << s << "->";
      }

        

        
        
        
    }

      std::cout << "" << std::endl;
      
      
      //std::cout << _mutationVector[i] << " " << _freqMap.find(_mutationVector[i])->second << std::endl;
    
  }
  
  return new_sol;
}

void SetCoverIlp::printSolutions(std::string origFile, std::string subdir, 
                                  int treeNum,
                                  const std::vector<std::vector<int>> solutions)
{
  std::string fname =  subdir + "/" + origFile + "_T" + 
                      std::to_string(treeNum + 1) + "_DF.csv";
  std::ofstream myfile;
  myfile.open (fname);
   if(!myfile) 
   { 
       std::cout<<"Error in creating file!!!"; 
       
   } 

  std::string dfFam = "";
  //int nr = _featureVector.size();
  std::vector<StringSet> features = _featureVector;
  for(const std::vector<int> df: solutions){
    for(auto it = df.begin(); it != df.end(); ++it){
      dfFam += "'";
      StringSet featurettes = features[(*it)];
      for( auto feat = featurettes.begin(); feat != featurettes.end(); ++feat ){
        dfFam += (*feat);
        if(std::next(feat) == featurettes.end()){
          dfFam += "'";
        }else{
          dfFam += " ";
        }
      }
      if(std::next(it) == df.end()){
        dfFam += "\n";

      }else{
        dfFam += ",";
      }
    }
  }
    myfile << dfFam;
    myfile.close();
}


