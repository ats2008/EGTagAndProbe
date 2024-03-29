#!/usr/bin/env python3
import os
import sys
import argparse
import json

version='v1'

useStr="\
Usage\n\
    ./jobMakerForCmsRuns.py <executable> <InputFileListFname> <destination> <NJOBS> <RUNS_PER_JOB> <jobPrefix>\n\
\n\
"

executable='recoStepUL2018.py'

NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
RUNS_PER_JOB=1
maxMeterialize=100
offsetStep=10000

pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
xrdRedirector="root://cms-xrd-global.cern.ch/"

FileSource ="bmm5FileList.txt"
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/runLumiList/'
tag=""



parser = argparse.ArgumentParser()
parser.add_argument("--exe", help="Executable")
parser.add_argument("--runFileMap", help="{ run:[ Files] } Map ")
parser.add_argument("--runs", help="runs To Process",default="*")
parser.add_argument("--runScript", help="template runscript")
parser.add_argument("--cfg", help="Configuration template file")
parser.add_argument("--dest", help="destination",default='./')
parser.add_argument("--jn", help="Max number of Jobs",default='10')
parser.add_argument("--rn", help="File per Job",default=1)
parser.add_argument("--fn", help="File per Job",default=1)
parser.add_argument("--maxEvt", help="Mac Events per Job , -1 for all",default=100)
parser.add_argument("--tag", help="Job Tag",default="")
parser.add_argument("--offsetStep", help="offsetStep for event number re-defenition",default="5000")
parser.add_argument("--maxMeterialize", help="Executable",default=None)

args = parser.parse_args()

executable=args.exe
FileSource=args.runFileMap
runScriptTemplate=args.runScript
cfgScriptTemplate=args.cfg
destination=args.dest
NJOBS=int(args.jn)
RUNS_PER_JOB=int(args.rn)
FILES_PER_JOB=int(args.fn)
NEVENTS_PER_JOB=int(args.maxEvt)
tag=args.tag
maxMeterialize=int(args.maxMeterialize)
offsetStep=int(args.offsetStep)

runs=args.runs.split(",")

if(not os.path.exists(destination)):
    os.system("mkdir -p "+destination)

destination=os.path.abspath(destination)


print(" Executable ",executable)
print("Source file list ",FileSource)
print("Run Script Template ",runScriptTemplate)
print("Configuration ",cfgScriptTemplate )
print("destination : ",destination)
print("NJOBS : ",NJOBS)
print("RUNS_PER_JOB : ",RUNS_PER_JOB)
print("NEVENTS_PER_JOB : ",NEVENTS_PER_JOB)
print("maxMeterialize : ",maxMeterialize)
print("tag : ",tag)
print()
print("Runs to process : ",",".join(runs))
print()

FileMap=open(FileSource,'r')
runFileMap=json.load(FileMap)
if '*' in runs:
    runs=list(runFileMap.keys())
FileMap.close()
print("Number avilable runs = ",len(runFileMap))

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"longlunch\"\n\
"
if maxMeterialize >0 :
    condorScriptString+="max_materialize="+str(maxMeterialize)+"\n"

f=open(runScriptTemplate,'r')
runScriptTxt=f.readlines()
f.close()
runScriptTxt=''.join(runScriptTxt)
runScriptTxt=runScriptTxt.replace("@@proxy_path",proxy_path)
runScriptTxt=runScriptTxt.replace("@@HOME",HOME)
runScriptTxt=runScriptTxt.replace("@@DESTINATION",destination)
#print(runScriptTxt)

f=open(cfgScriptTemplate,'r')
configurationTxt=f.readlines()
f.close()
configurationTxt=''.join(configurationTxt)
#print(configurationTxt)


head='Condor/Jobs'+tag
if not os.path.exists(head ):
    os.system('mkdir -p '+head)

condorScriptName=head+'/job'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)

n=int(len(runs)/RUNS_PER_JOB) + 1
if n < NJOBS:
    NJOBS=n
print("Making Run ",NJOBS," Jobs ")

njobs=0
for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    dirNameM= pwd+'/'+head+'/Job_'+str(i)
    
    if(ii%10==0) : print("\nJob Made : ",end = " " )
    print(ii,end =" ")
    
    k=0
    rIdx=''
    fnamesToLoad=[]
    while k <RUNS_PER_JOB:
        if len(runs) < 1 :
            break
        r=runs.pop(0)
        rIdx=r+'_'
        for fname in runFileMap[r]:
          fnamesToLoad.append(fname)
        k+=1
    rIdxM=rIdx[:-1]
    if len(rIdx) > 6:
        rIdx=str(i)
    
    filesLeft = len(fnamesToLoad)> 0
    nJobInRun=0
    while filesLeft:
        nJobInRun+=1
        tmp=''
        for j in range(FILES_PER_JOB):
            fname=fnamesToLoad.pop(0)
            tmp+=fname+"\n"
            filesLeft = len(fnamesToLoad) > 0
            if not filesLeft:
                break
        dirName= dirNameM+f"_{nJobInRun}"
        rIdx=rIdxM+f"_{nJobInRun}"
        if os.path.exists(dirName):
            k=True
        else:
            os.system('mkdir '+dirName)
        
        cfgFileName='customize_'+str(i)+'.cfg'
        cfgFile=open(dirName+'/'+cfgFileName,'w')
        
        tmp=configurationTxt.replace("@@FNAMES",tmp[:-1])
        tmp=tmp.replace("@@TAG",tag)
        tmp=tmp.replace("@@IDX",str(rIdx))
        tmp=tmp.replace("@@OFFSET_IDX",str(i*offsetStep))
        tmp=tmp.replace("@@MAXEVENTS",str(NEVENTS_PER_JOB))
        cfgFile.write(tmp)
        cfgFile.close()   
        
        runScriptName=dirName+'/'+tag+'run'+str(i)+f'_{rIdx}'+'.sh'
        if os.path.exists(runScriptName+'.sucess'):
           os.system('rm '+runScriptName+'.sucess')
        runScript=open(runScriptName,'w')
        tmp=runScriptTxt.replace("@@DIRNAME",dirName)
        tmp=tmp.replace("@@IDX",str(i))
        tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
        tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
        tmp=tmp.replace("@@EXECUTABLE",executable)
        runScript.write(tmp)
        runScript.close()
        os.system('chmod +x '+runScriptName)
        if maxMeterialize < 0:
            condorScript.write("queue filename matching ("+runScriptName+")\n")
        njobs+=1
if maxMeterialize >=0:
    condorScript.write("queue filename matching ("+head+"/*/*.sh)\n")
    
print()
print(" Number of jobs made : ", njobs)
print(" Condor submit file  : ", condorScriptName)
condorScript.close()
