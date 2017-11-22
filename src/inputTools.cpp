#include "inputTools.h"

/**
 * searchParameter searches a file for a certain flag string and returns a
 * string after the '=' sign
 * @param fileLocation
 * @param searchText
 * @return
 */
std::string searchParameter(const std::string &fileLocation,const std::string &searchText){
  
  std::string output = "No string found";
  std::string currentLine;
  std::vector<std::string> splittedCurrLine;
  std::ifstream file;
  bool found = false;

  file.open(fileLocation.c_str());
  // open error handling
  if(!file){ std::cout << "File is not available! Open it first." << std::endl; }
  // if no error occurs...
  else{
    // get the first line
    // cycles through the file until the end or if the found flag is true
    while( getline(file,currentLine) && !found ){
      // if no occurrence of the key is found, continues while loop
      if( currentLine.find(searchText) == std::string::npos ) continue;
      // if something is found, set flag to true
      else{
	found = true;
	split(splittedCurrLine,currentLine,'=');
	output = splittedCurrLine[1];
      }
    }
  }
  file.close();
  return output;
}
