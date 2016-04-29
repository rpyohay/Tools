#include <string>
#include "TFile.h"
#include "TCanvas.h"

//get object from canvas
template<typename T>
T* getObjectFromCanvas(TFile& file, const string& objName, const string& canvasName, 
		       const unsigned int pad)
{
  TCanvas* canvas = NULL;
  file.GetObject(canvasName.c_str(), canvas);
  T* obj = NULL;
  if (canvas != NULL) {
    canvas->cd(pad);
    obj = (T*)canvas->GetPrimitive(objName.c_str())/*->Clone()*/;
  }
  return obj;
}
