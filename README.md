Forked from charis's original [repo@GitLab ](https://gitlab.cern.ch/ckoraka/EGTagAndProbe/)


# EGTagAndProbe
Set of tools to evaluate L1EG trigger performance on T&amp;P

Based on TauTagAndProbe package developed by L. Cadamuro & O. Davignon

### Install instructions
```
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
cmsenv
git clone git@github.com:ats2008/EGTagAndProbe.git
cd EGTagAndProbe 
git checkout raw_reco_TandP # to use developements targeting the RAW-RECO dataset
scram b -j 4
```

### Producing TagAndProbe ntuples with emulated L1EG


1. Compile and cd to relevant directory:
```
cd CMSSW_X_Y_Z/src
scram b -j4
${CMSSW_BASE}/src/EGTagAndProbe/EGTagAndProbe/test
```

2. Set flag isMC and isMINIAOD according to sample in `TandPRobeNtuplizerWithReEmulation_data.py`

```
cmsRun TandPRobeNtuplizerWithReEmulation_data.py
```



### Submit job on the Grid
Modify crab3_config.py: change requestName, inputDataSet, outLFNDirBase, outputDatasetTag, storageSite
```
cd ${CMSSW_BASE}/src/EGTagAndProbe/EGTagAndProbe/test
source /cvmfs/cms.cern.ch/crab3/crab.sh
voms-proxy-init -voms cms
crab submit -c crab3_config.py
```

### Producing turn-on plots
``
cd EGTagAndProbe/test/fitter/
``
Use the readme `EGTagAndProbe/test/fitter/Readme.md`


### Producing TagAndProbe ntuples with unpacked L1EG (no re-emulation)
This workflow is depricated in this baranch , the re-emulated ntuples also has branches for the unpacked studies.

cd to relevant directory
```
${CMSSW_BASE}/src/EGTagAndProbe/EGTagAndProbe
```
2. Set flag isMC and isMINIAOD according to sample in test/test.py
3. HLT path used specified in python/MCAnalysis_cff.py (MC) or python/tagAndProbe_cff.py (data)
4. Launch test/test.py

```
cd test
cmsRun test.py
```


