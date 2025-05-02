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

inputFileNames=config['rate_files']
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
    histStore[tag]=getTheObjectsFromFile(histFileStore[tag])


plotRateGeneric(histStore,plotName="SingleEG_rate_inclusive",hname="SingleEG_rate_inclusive",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG Rate"],xlims=[5,80],ylims=[1e0,1e5],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)


plotRateGeneric(histStore,plotName="SingleEG_rate_inclusive_ROI",hname="SingleEG_rate_inclusive",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG Rate"],xlims=[25,50],ylims=[0,50],doLogY=False,
                    era=era,lumi=lumi,prefix=prefix)


plotRateGeneric(histStore,plotName="SingleEG_rate_TightIso",hname="SingleEG_rate_TightIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG TightIso Rate"],xlims=[5,80],ylims=[1e0,1e5],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)


plotRateGeneric(histStore,plotName="SingleEG_rate_TightIso_ROI",hname="SingleEG_rate_TightIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG TightIso Rate"],xlims=[25,50],ylims=[0,50],doLogY=False,
                    era=era,lumi=lumi,prefix=prefix)

plotRateGeneric(histStore,plotName="SingleEG_rate_LooseIso",hname="SingleEG_rate_LooseIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG LooseIso Rate"],xlims=[5,80],ylims=[1e0,1e5],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)


plotRateGeneric(histStore,plotName="SingleEG_rate_LooseIso_ROI",hname="SingleEG_rate_LooseIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Single EG LooseIso Rate"],xlims=[8,40],ylims=[1e1,1e4],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)
  
plotRateGeneric(histStore,plotName="DoubleEG_rate_inclusive",hname="DoubleEG_rate_inclusive",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Double EG Rate"],xlims=[5.0,80.0],ylims=[1e-1,5e4],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)

plotRateGeneric(histStore,plotName="DoubleEG_rate_inclusive_ROI",hname="DoubleEG_rate_inclusive",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Double EG Rate"],xlims=[5,30],ylims=[0,50],doLogY=False,
                    era=era,lumi=lumi,prefix=prefix)


plotRateGeneric(histStore,plotName="DoubleEG_rate_LooseIso",hname="DoubleEG_rate_LooseIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Double EG LooseIso Rate"],xlims=[5,80],ylims=[1e-1,1e3],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)

plotRateGeneric(histStore,plotName="DoubleEG_rate_LooseIso_ROI",hname="DoubleEG_rate_LooseIso",
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias","Double EG LooseIso Rate"],xlims=[5,30],ylims=[0,20],doLogY=False,
                    era=era,lumi=lumi,prefix=prefix)
  


hlist=[
    'SingleEG_rate_inclusive',
    'SingleEG_rate_TightIso',
    'DoubleEG_rate_inclusive',
    'DoubleEG_rate_LooseIso',
]
hist_legend={
   'SingleEG_rate_inclusive': 'SingleEG' , 
   'SingleEG_rate_TightIso' : 'SingleEG,Tight' ,
   'DoubleEG_rate_inclusive': 'DoubleEG' ,
   'DoubleEG_rate_LooseIso' : 'DoubleEG,Loose' ,
}

for tag in tagsToPlot:
    plotRatesAll(histStore,plotName=f"{tag}_rates",hlist=hlist,hist_legend=hist_legend,addTagLegend=False,
                        tagsToPlot=[tag],colors=colours,legend=legend,
                        description=["L1 EG","Zero-Bias"],xlims=[5,80],ylims=[1e0,1e5],doLogY=True,
                        era=era,lumi=lumi,prefix=prefix)

plotRatesAll(histStore,plotName="all_rates",hlist=hlist,hist_legend=hist_legend,addTagLegend=True,
                    tagsToPlot=tagsToPlot,colors=colours,legend=legend,
                    description=["L1 EG","Zero-Bias"],xlims=[5,80],ylims=[1e0,1e5],doLogY=True,
                    era=era,lumi=lumi,prefix=prefix)
  




