import json,os
from Util import *
from plotter import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputConfigJson", help="Input Configuration File")
parser.add_argument('-o',"--outputFolder", help="Destination path for plots",default='results/plots/responses/')
args = parser.parse_args()

with open(args.inputConfigJson) as f:
    config=json.load(f)

prefix=args.outputFolder

inputFileNames=config['response_files']
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
    histFileStore[tag]=ROOT.TFile(inputFileNames[tag],'READ')
    print("loading  tag ",tag)
    histStore[tag]=getTheObjectsFromFile(histFileStore[tag])['0']

xlabel_base=config['xlabel']
plotEtVsdEt(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEtResolutions(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEtResolutionsGreaterThan20(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEtResolutionInPt(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEtResolutionInEta(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
plotEtaPhiresolutions(histStore,tagsToPlot,colours,legend,era,lumi,prefix,xlabel_base)
       
