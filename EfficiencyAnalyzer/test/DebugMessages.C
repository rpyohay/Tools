#include <string>

//error codes
enum ErrorCode {SUCCESS = 0, CANNOT_RETRIEVE_OBJECT};

//error: file cannot be opened
string errorCannotOpenFile(const string& fnName, const string& fileName)
{
  return ("Error in " + fnName + ":\nCannot open file \"" + fileName + "\".\n");
}

//error: object could not be retrieved from file
string errorCannotRetrieveObject(const string& fnName, const string& fileName, 
				 const string& objName)
{
  return ("Error in " + fnName + ":\nCannot retrieve object \"" + objName + "\" from file \"" + 
	  fileName + "\".\n");
}
