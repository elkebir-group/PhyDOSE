/*
 * setcover.h
 *
 *  Created on: 24-nov-2019
 *      Author: M. El-Kebir
 */

#ifndef SETCOVER_H
#define SETCOVER_H

#include "utils.h"
#include "clonetree.h"

class SetCover
{
public:
  SetCover(const CloneTree& T,
           const CloneTreeVector& scriptT,
           const std::map<std::string, double>& freqMap);
  
  void writeDOT(std::ostream& out) const;
  
private:
  void init();
  
protected:
  const CloneTree& _T;
  
  const CloneTreeVector& _scriptT;
  
  std::vector<StringSet> _featureVector;
  
  StringVector _mutationVector;
  
  const std::map<std::string, double>& _freqMap;
  
  Graph _G;
  
  std::vector<Graph::Node> _nodesTree;
  
  std::vector<Graph::Node> _nodesFeatures;
};

#endif // SETCOVER_H
