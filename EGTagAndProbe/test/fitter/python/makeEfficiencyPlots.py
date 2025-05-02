import json,os
from Util import *
from plotter import *
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputConfigJson", help="Input Configuration File")
parser.add_argument('-o',"--outputFolder", help="Destination path for plots",default='results/plots/efficiencies/')
args = parser.parse_args()

with open(args.inputConfigJson) as f:
    config=json.load(f)

prefix=args.outputFolder

inputFileNames=config['efficiency_files']
tagsToPlot=config['tagsToPlot']
colours_=config['colours']
colours ={c:getattr(ROOT,colours_[c]) for c in colours_}
legend=config['legend']
era=None
if 'era' in config:
    era=config['era']
lumi=None
if 'lumi' in config:
    lumi=config['lumi']

if not os.path.exists(prefix):
    os.system('mkdir -p '+prefix)

histStore={}
histFileStore={}
for tag in inputFileNames:
    if tag in histFileStore:
        histFileStore[tag].Close()
    print("loading  tag ",tag," with filename ",inputFileNames[tag])
    histFileStore[tag]=ROOT.TFile(inputFileNames[tag],'READ')
    histStore[tag]=getTheObjectsFromFile(histFileStore[tag])['0']


xlabel_base=config['xlabel']
plotEffVsPu(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEffBarrelECapSplit(histStore,[5,15,25,32,36,40],tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEffBarrelECapSplitTightIso(histStore,[18,32,34,40],tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEffBarrelECapSplitLooseIso(histStore,[18,32,34,40],tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEffComparisionLooseVsTightVsInclusive(histStore,[18,32,34,40],tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
       
