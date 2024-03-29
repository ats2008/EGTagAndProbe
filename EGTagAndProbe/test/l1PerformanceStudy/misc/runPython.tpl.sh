#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh 
set -x
export HOME=@@HOME
export X509_USER_PROXY=@@proxy_path
cd @@DIRNAME
set +x
eval `scramv1 runtime -sh`
set -x
TMPDIR=`mktemp -d`
cd $TMPDIR
cp @@DIRNAME/@@CFGFILENAME .
mv @@RUNSCRIPT @@RUNSCRIPT.busy 
python3 @@EXECUTABLE @@CFGFILENAME
if [ $? -eq 0 ]; then 
    cp *.root @@DESTINATION
    if [ $? -eq 0 ] ; then
        mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess 
        echo OK
    else
        mv @@RUNSCRIPT.busy @@RUNSCRIPT
        echo FAIL
    fi
else
    mv @@RUNSCRIPT.busy @@RUNSCRIPT 
    echo FAIL
fi
