/*
 * setcover.cpp
 *
 *  Created on: 24-nov-2019
 *      Author: M. El-Kebir
 */

#include "setcover.h"

SetCover::SetCover(const CloneTree& T,
                   const CloneTreeVector& scriptT,
                   const std::map<std::string, double>& freqMap)
  : _T(T)
  , _scriptT(scriptT)
  , _featureVector()
  , _mutationVector()
  , _freqMap(freqMap)
  , _G()
  , _nodesTree()
  , _nodesFeatures()
{
  init();
}

void SetCover::writeDOT(std::ostream& out) const
{
  out << "graph G {" << std::endl;
  
  for (int i = 0; i < _nodesFeatures.size(); ++i)
  {
    std::vector<std::pair<int, Node>> levels;
    const StringSet& feature_i = _featureVector[i];
    for (const std::string& mutation : feature_i)
    {
      Node v = _T.getNodeByLabel(mutation);
      levels.push_back(std::make_pair(_T.level(v), v));
    }
    std::sort(levels.begin(), levels.end());
    
    out << "\t" << i << " [label=\"";
    for (const auto& pair : levels)
    {
      const std::string& mutation = _T.label(pair.second);
      if (mutation == _mutationVector[i])
        std::cout << "*";
      
      std::cout << mutation << "\\n";
    }
    out << _freqMap.find(_mutationVector[i])->second;
    out << "\"]" << std::endl;
  }
  
  for (int j = 0; j < _nodesTree.size(); ++j)
  {
    out << "\tT" << j << std::endl;
  }
  
  lemon::DynArcLookUp<Graph> lookUp(_G);
  for (int i = 0; i < _nodesFeatures.size(); ++i)
  {
    Graph::Node v_i = _nodesFeatures[i];
    for (int j = 0; j < _nodesTree.size(); ++j)
    {
      Graph::Node v_j = _nodesTree[j];
      if (lookUp(v_i, v_j) != lemon::INVALID)
      {
        out << "\t" << i << " -- T" << j << std::endl;
      }
    }
  }
  
  out << "}" << std::endl;
}

void SetCover::init()
{
  // 0. Set feature vector
  std::vector<std::pair<int, Node>> levels;
  for (NodeIt v(_T.tree()); v != lemon::INVALID; ++v)
  {
    levels.push_back(std::make_pair(_T.level(v), v));
  }
  std::sort(levels.begin(), levels.end());
  
  for (const auto& pair : levels)
  {
    Node v = pair.second;
    _mutationVector.push_back(_T.label(v));
    _featureVector.push_back(_T.nodeToMutations(v));
  }
  
  // 1. One node for every tree in _scriptT
  for (const CloneTree& otherT : _scriptT)
  {
    _nodesTree.push_back(_G.addNode());
  }
  
  // 2. Another node for every feature in _T
  for (const StringSet& clone : _featureVector)
  {
    _nodesFeatures.push_back(_G.addNode());
  }
  
  // 3. Add edges
  for (int i = 0; i < _nodesFeatures.size(); ++i)
  {
    const StringSet& feature_i = _featureVector[i];
    Graph::Node v_i = _nodesFeatures[i];
    for (int j = 0; j < _nodesTree.size(); ++j)
    {
      Graph::Node v_j = _nodesTree[j];
      if (_scriptT[j].getCloneSet().count(feature_i) == 0)
      {
        _G.addEdge(v_i, v_j);
      }
    }
  }
}
