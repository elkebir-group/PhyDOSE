/*
 * inputinstance.h
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#ifndef INPUTINSTANCE_H
#define INPUTINSTANCE_H

#include "utils.h"
#include "clonetree.h"
#include "frequencymatrix.h"

class InputInstance
{
public:
  /// Constructor
  InputInstance();
  
  /// Return clone trees
  const CloneTreeVector& getTrees() const
  {
    return _scriptT;
  }
  
  /// Return frequency matrix
  const FrequencyMatrix& getFrequencies() const
  {
    return _F;
  }
  
  /// Return total number of trees across all patients
  int getNrTrees() const
  {
    return _scriptT.size();
  }
  
private:
  CloneTreeVector _scriptT;
  
  FrequencyMatrix _F;
  
  friend std::ostream& operator<<(std::ostream& out, const InputInstance& input);
  friend std::istream& operator>>(std::istream& in, InputInstance& input);
};

std::ostream& operator<<(std::ostream& out, const InputInstance& input);

std::istream& operator>>(std::istream& in, InputInstance& input);

#endif // INPUTINSTANCE_H
