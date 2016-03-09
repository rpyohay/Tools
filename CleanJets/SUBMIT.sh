#!/bin/bash

export SCRAM_ARCH=slc6_amd64_gcc491
cd /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/ 
eval `scramv1 runtime -sh`
cd -
source /afs/cern.ch/cms/caf/setup.sh
cp /afs/cern.ch/user/k/ktos/GroupDir/CMSSW_7_6_3/src/Tools/CleanJets/BSUB/DIRNAME/PRODUCER.py . 
cmsRun PRODUCER.py
cmsStage -f DIRNAME_NUM.root /store/user/ktos/DIRNAME
rm PRODUCER.py 
exit 0
