{
  //initial
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
  gROOT->ProcessLine("#include <utility>");
  gSystem->Load((macroPath + "STLDictionary.so").c_str());
  gSystem->Load((macroPath + "Plot_C.so").c_str());

  //unit strings
  string unitPTTau("#tau_{X} p_{T} (GeV)");
  string unitEtaTau("#tau_{X} #eta");

  //map of bin labels for certain efficiency plots
  map<string, vector<string> > binLabelMap;

  //map of inputs to 1D efficiency histograms
  map<string, pair<string, string> > effHistMap1DTau;
  effHistMap1DTau["visiblePTnumerator"] = make_pair(string("visiblePTdenominator"), unitPTTau);
  effHistMap1DTau["visibleEtanumerator"] = make_pair(string("visibleEtadenominator"), unitEtaTau);

  //map of inputs to 2D efficiency histograms
  map<pair<string, string>, pair<string, string> > effHistMap2DTau;

  //map of inputs to 1D histograms
  map<string, string> hist1DMapTau;
  hist1DMapTau["visiblePTnumerator"] = unitPTTau;
  hist1DMapTau["visiblePTdenominator"] = unitPTTau;
  hist1DMapTau["visibleEtanumerator"] = unitEtaTau;
  hist1DMapTau["visibleEtadenominator"] = unitEtaTau;

  //make efficiency plots for WH and Z-->tautau
  vector<string> effInputFiles;
  effInputFiles.push_back("gentau_mi_Z.root");
  effInputFiles.push_back("recomatched_gen_tau_mi.root");
  for (vector<string>::const_iterator iFile = effInputFiles.begin(); iFile != effInputFiles.end(); 
       ++iFile) {
    const unsigned int strLen = iFile->find(".root");
    string outputFileName(iFile->substr(0, strLen) + "_final.root");
    float weight = 2.26606084150487; //SM Z M > 50 GeV cross section
    if ((iFile - effInputFiles.begin()) == 1) weight = 0.138806; /*SM WH cross section * 
								   BR(W-->leptons)*/
    plotNice(*iFile, effHistMap1DTau, effHistMap2DTau, binLabelMap, hist1DMapTau, outputFileName, 
	     "noPDF", weight);
  }

  //map of canvases to files to draw options
  map<pair<string, string>, vector<pair<pair<string, Option_t*>, 
    pair<Color_t, Style_t> > > > canvasMap;
  vector<pair<pair<string, Option_t*>, pair<Color_t, Style_t> > > graphMap;
  graphMap.push_back(make_pair(make_pair("gentau_mi_Z_final.root", "AP"), 
			       make_pair(kRed, kFullTriangleDown)));
  graphMap.push_back(make_pair(make_pair("recomatched_gen_tau_mi_final.root", "PSAME"), 
  			       make_pair(kBlue, kFullTriangleUp)));
  canvasMap[make_pair("eff_visiblePTnumerator_over_visiblePTdenominator", unitPTTau)] = graphMap;
  canvasMap[make_pair("eff_visibleEtanumerator_over_visibleEtadenominator", unitEtaTau)] = graphMap;

  //overlay tau efficiency graphs
  plotNiceDifferentFiles(canvasMap, binLabelMap, "WH_vs_ZTauTau_tauX_eff.root", "noPDF");
}
