/* inputTools.h
 * 
 * Library to deal with input file operations
 * 
 */

#ifndef _INPUT_TOOLS_H
#define _INPUT_TOOLS_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

struct split
{
  enum empties_t { empties_ok, no_empties };
};

/// Split something into a Container of something using a certain delimiter
/// from <http:\\www.cplusplus.com/faq/sequences/strings/split/>
/// @param result the container where output should be placed
/// @param s the object to be searched into
/// @param delimiter
/// @param empties flag to account for empty spaces
template <typename Container>
Container& split( Container& result, const typename Container::value_type& s,
		typename Container::value_type::value_type delimiter,
		split::empties_t empties = split::empties_ok )
{
  result.clear();
  std::istringstream ss( s );
  while (!ss.eof())
  {
    typename Container::value_type field;
    getline( ss, field, delimiter );
    // TODO next line is not compiling due to the following error message:
    // error: ‘split’ is not a class, namespace, or enumeration
    //
//    if ((empties == split::no_empties) && field.empty()) continue;
    result.push_back( field );
  }
  return result;
}

/// Given a text file, searches for a specific string and returns its associated values
/// @param inFile input file as its location char string
/// @param searchText string to be searched
/// @returns the output as a strings
std::string searchParameter(const std::string &fileLocation,const std::string &searchText);

#endif
