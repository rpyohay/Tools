#include <string>

//add trailing / to save path if missing
void formatSavePath(string& savePath)
{
  if (savePath.find('/') != (savePath.length() - 1)) savePath+="/";
}
