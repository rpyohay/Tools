{
  gROOT->Reset();

  //get CMSSW path
  const char* CMSSWPathCString = gSystem->Getenv("CMSSW_BASE");
  if (!CMSSWPathCString) {
    CMSSWPathCString = "";
    cout << "Error: environment variable CMSSW_BASE is not set.  ";
    cout << "Please run cmsenv from within your CMSSW project area.\n";
  }
  string CMSSWPathCPPString(CMSSWPathCString);

  //load
  string macroPath(CMSSWPathCPPString + "/src/Tools/EfficiencyAnalyzer/test/");
  // gSystem->Load((macroPath + "STLDictionary.so").c_str());

  //compile
  // gROOT->LoadMacro((macroPath + "PlotStyle.C++").c_str());
  // gROOT->LoadMacro((macroPath + "DebugMessages.C++").c_str());
  // gROOT->LoadMacro((macroPath + "StringManipulation.C++").c_str());
  // gROOT->LoadMacro((macroPath + "Plot.C++").c_str());
  gROOT->LoadMacro((macroPath + "tdrstyle.C").c_str());
  gROOT->LoadMacro((macroPath + "CMS_lumi.C").c_str());
}
