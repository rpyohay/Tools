#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/NMSSM_Analysis/CMSSW_7_4_12_patch4/src/BoostedTauAnalysis/CleanJets/BSUB/DIRNAME/PRODUCER.py . 
cmsRun PRODUCER.py
cmsStage -f DIRNAME_NUM.root /store/user/ktos/DIRNAME
rm PRODUCER.py 
exit 0
