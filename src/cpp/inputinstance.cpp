/*
 * inputinstance.cpp
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#include "inputinstance.h"

InputInstance::InputInstance()
  : _scriptT()
{
}

std::ostream& operator<<(std::ostream& out, const InputInstance& input)
{
  out << std::endl;
  out << input._scriptT;
  out << input._F;
  return out;
}

std::istream& operator>>(std::istream& in, InputInstance& input)
{
  std::string line;
  getline(in, line);
  in >> input._scriptT;
  in >> input._F;

  return in;
}
