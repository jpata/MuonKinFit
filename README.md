# MuonKinFit
Dimuon kinematic fit as a CMSSW module. Based on the [UFLX2MuMu/Ntupliser](https://github.com/UFLX2MuMu/Ntupliser/blob/49e4fd57ffbdc18a98bcea64db8c736090d42eaf/DiMuons/src/MuonHelper.cc#L29]).

~~~
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
mkdir MyAnalysis
cd MyAnalysis
git clone https://github.com/jpata/MuonKinFit.git
scram b
cmsRun MuonKinFit/test/testMuonKinFit_cfg.py
~~~

