/*
 * setcoverilp.h
 *
 *  Created on: 5-dec-2019
 *      Author: M. El-Kebir
 */

#ifndef SETCOVERILP_H
#define SETCOVERILP_H

#include "setcover.h"
#include <ilcplex/ilocplex.h>
#include <fstream>

class SetCoverIlp : public SetCover
{
public:
  SetCoverIlp(const CloneTree& T,
              const CloneTreeVector& scriptT,
              const std::map<std::string, double>& freqMap,
              const std::vector<std::vector<int>> solution_set);
  
  std::vector<int> solve();
  void printSolutions(std::string, std::string, int ,const std::vector<std::vector<int>> );


  
private:
  void init(std::vector<std::vector<int>> solution_set);
  
private:
  /// Number of features in distinguishing feature set
  //const int _k;
  /// Environment
  IloEnv _env;
  /// CPlex model
  IloModel _model;
  /// Solver
  IloCplex _cplex;
  /// Cover variables
  IloBoolVarArray _x;
  /// Minimum weight variable
  //IloNumVar _z;
  /// Objective value
  double _objValue;
};

#endif // SETCOVERILP_H
