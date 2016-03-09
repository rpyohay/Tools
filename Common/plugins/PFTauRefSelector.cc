#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"

typedef SingleObjectSelector<
  reco::PFTauCollection, 
  StringCutObjectSelector<reco::PFTau>,
  reco::PFTauRefVector
  > PFTauRefSelector;

DEFINE_FWK_MODULE(PFTauRefSelector);
