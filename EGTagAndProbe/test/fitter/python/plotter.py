from Util import *

baseCanvas='Run3'

##  EFFICIENCIES
def plotEffVsPu(histStoreTurnOns,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
   
    plot=getDefaultPlot(prefix=prefix,name='effVsPU',cPars=cPars) 
    plot.legendPosition = (0.5,0.41,0.85,0.52)
    plot.descPosition   = (0.55,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","|#eta| <2.4","E_{T}^{"+varToChose+"} > 40 , E_{T}^{L1} > 32"]
    plot.yTitle = "#epsilon"  
    plot.xTitle = "PU" 
    plot.xRange = (5.0,100.0) 
    plot.logx = False
    plot.yRange = (0.0,1.01) 
    plots.append(plot) 

    i=2
    hname='Graph_from_hist_L1_32vsPUoffline40GeVEfficiency'
    
    for tag in tagsToPlot: #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        print(tag)
        hist=histStoreTurnOns[tag][hname].Clone()
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
        pparams['Options']='pe'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=False ; aplot.drawFit=True
        aplot.verbose= False
        aplot.drawLegend = True
        plots[-1].addPlot(aplot)
    
    canvas = []
    for plot in plots:
        canvas.append(plot.plot()) 
    


def plotEffBarrelECapSplit(histStoreTurnOns,ETsToPlot,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    for eT in ETsToPlot:
        plot=getDefaultPlot(prefix=prefix,name='effComponentBarrel_'+str(eT)+varToChose,cPars=cPars)
        plot.legendPosition = (0.60,0.35,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_{T} > "+str(eT),"Barrel"]
        plot.yTitle = "#epsilon"
        plot.xTitle = "E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,500.0)
        plot.logx = True
        plot.yRange = (0.0,1.01)
        plots.append(plot)

        i=2
        hname='Graph_from_hist_L1Et'+str(eT)+'_default_BarrelEfficiency'

        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


        plot=getDefaultPlot(prefix=prefix,name='effComponentECap_'+str(eT)+varToChose,cPars=cPars)
        plot.legendPosition = (0.60,0.35,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_T > "+str(eT),"End Cap"]
        plot.yTitle = "#epsilon"
        plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,500.0)
        plot.logx = True
        plot.yRange = (0.0,1.01)
        plots.append(plot)
        hname='Graph_from_hist_L1Et'+str(eT)+'_default_ECapEfficiency'

        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


        plot=getDefaultPlot(prefix=prefix,name='effComponentInclusive_'+str(eT)+varToChose,cPars=cPars)
        plot.legendPosition = (0.60,0.35,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_{T} > "+str(eT),"Inclusive"]
        plot.yTitle = "#epsilon"
        plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,500.0)
        plot.logx = True
        plot.yRange = (0.0,1.01)
        plots.append(plot)
        hname='Graph_from_hist_L1Et'+str(eT)+'_defaultEfficiency'

        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


    canvas = []
    for plot in plots:
        canvas.append(plot.plot())

def plotEffBarrelECapSplitTightIso(histStoreTurnOns,ETsToPlot,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams('Run3Data',lumi)
    rBin=1
    for eT in ETsToPlot:
        plot=getDefaultPlot(prefix=prefix,name='effComponentBarrel_tightIso_'+str(eT)+varToChose,cPars=cPars) 
        plot.legendPosition = (0.60,0.48,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_T > "+str(eT),"Tight Isolation","Barrel"]
        plot.yTitle = "#epsilon"  
        plot.xTitle = "E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,80.0) 
        plot.logx = False
        plot.yRange = (0.0,1.01) 
        plots.append(plot) 

        i=2
        hname='Graph_from_hist_L1Et'+str(eT)+'_tightiso_BarrelEfficiency'
    
        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)
        
        
        plot=getDefaultPlot(prefix=prefix,name='effComponentECap_tight_'+str(eT)+varToChose,cPars=cPars) 
        plot.legendPosition = (0.60,0.35,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_T > "+str(eT),"Tight Isolation","Endcap"]
        plot.yTitle = "#epsilon"  
        plot.xTitle = "E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,100.0) 
        plot.logx = False
        plot.yRange = (0.0,1.01) 
        plots.append(plot) 
        hname='Graph_from_hist_L1Et'+str(eT)+'_tightiso_ECapEfficiency'
    
        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot()) 
    
def plotEffBarrelECapSplitLooseIso(histStoreTurnOns,ETsToPlot,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams('Run3Data',lumi)
    rBin=1
    for eT in [15,22,25]:
        plot=getDefaultPlot(prefix=prefix,name='effComponentBarrel_looseIso_'+str(eT)+varToChose,cPars=cPars) 
        plot.legendPosition = (0.60,0.45,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_T_DoubleEG_"+str(eT-10)+"_"+str(eT),"Loose Isolation","Barrel"]
        plot.desc =["E_{T} > "+str(eT),"Loose Isolation","Barrel"]
        plot.yTitle = "#epsilon"  
        plot.xTitle = "E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,80.0) 
        plot.logx = False
        plot.yRange = (0.0,1.01) 
        plots.append(plot) 

        i=2
        hname='Graph_from_hist_L1Et'+str(eT)+'_looseiso_BarrelEfficiency'
    
        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)
        
        
        plot=getDefaultPlot(prefix=prefix,name='effComponentECap_loose_'+str(eT)+varToChose,cPars=cPars) 
        plot.legendPosition = (0.60,0.45,0.85,0.60)
        plot.descPosition   = (0.65,0.75)
        plot.desc =["E_T > "+str(eT),"Loose Isolation","Endcap"]
        plot.yTitle = "#epsilon"  
        plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
        plot.xRange = (5.0,100.0) 
        plot.logx = False
        plot.yRange = (0.0,1.01) 
        plots.append(plot) 
        hname='Graph_from_hist_L1Et'+str(eT)+'_looseiso_ECapEfficiency'
    
        for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


    canvas = []
    for plot in plots:
        canvas.append(plot.plot()) 
    

def plotEffComparisionLooseVsTightVsInclusive(histStoreTurnOns,ETsToPlot,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    col=[ROOT.kGreen+2,ROOT.kBlue, ROOT.kRed]
    plots=[]
    cPars=getCanvasParams('Run3Data',lumi)
    rBin=1
    for tag in  tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        for eT in ETsToPlot:
            plot=getDefaultPlot(prefix=prefix,name='tightVsLoose_'+str(eT)+'_'+tag+'_'+varToChose,cPars=cPars) 
            plot.legendPosition = (0.65,0.51,0.90,0.60)
            plot.descPosition   = (0.65,0.75)
            plot.desc =["E_T > "+str(eT),legend[tag]]
            plot.yTitle = "#epsilon"  
            plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
            plot.xRange = (5.0,500.0) 
            plot.logx = True
            plot.yRange = (0.0,1.01) 
            plots.append(plot) 

            i=0
            hname='Graph_from_hist_L1Et'+str(eT)+'_defaultEfficiency'
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']='Inclusive'
            pparams['MarkerColor']= col[i]; pparams['LineColor'] = col[i]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


            i=1
            hname='Graph_from_hist_L1Et'+str(eT)+'_looseiso_defEfficiency'
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']='Loose Isolation'
            pparams['MarkerColor']= col[i]; pparams['LineColor'] = col[i]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)

            i=2
            hname='Graph_from_hist_L1Et'+str(eT)+'_tightiso_defEfficiency'
            hist=histStoreTurnOns[tag][hname].Clone()
    #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']='Tight Isolation'
            pparams['MarkerColor']= col[i]; pparams['LineColor'] = col[i]
            pparams['Options']='pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=True
            aplot.verbose= False
            plots[-1].addPlot(aplot)


        canvas = []
        for plot in plots:
            canvas.append(plot.plot())




##  RESPONSES 
def plotEtVsdEt(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    xAxisLbl=varToChose
    if True:  # Inclusive Et barrel Ecap vs Gen Et
        plots=[]
        cPars=getCanvasParams(baseCanvas,era,lumi)
        rBin=1
        
        plot=getDefaultPlot(prefix=prefix,name='EtVsdEtScaleBarrel_'+xAxisLbl,cPars=cPars) 
        plot.legendPosition = (0.45,0.75,0.89,0.90)
        plot.descPosition   = (0.55,0.70)
        plot.desc =["Barrel"]
        plot.yTitle = "E_{T}^{L1} / E_{T}^{"+xAxisLbl+"}"  
        plot.xTitle = "E_{T}^{"+xAxisLbl+"}"
        plot.xRange = (0.,150.0) 
        plot.logx = False
        plot.yRange = (0.7,1.2) 
        plots.append(plot) 
        
    
        i=2
        for tag in tagsToPlot:
            hname='EtVsdEtScaleBarrel'
            hist=histStore[tag][hname].Clone()
        #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
            pparams['MarkerStyle']=8
            pparams['Options']='HIST pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.drawLegend = True
            aplot.doFit=False ; aplot.drawFit=False
            aplot.verbose= False
            plots[-1].addPlot(aplot)
    
    
        plot=getDefaultPlot(prefix=prefix,name='EtVsdEtScaleECap_'+xAxisLbl,cPars=cPars)   
        plot.legendPosition = (0.45,0.75,0.89,0.90)
        plot.descPosition   = (0.55,0.70)
        plot.desc =["Endcap"]
        plot.yTitle = "E_{T}^{L1} / E_{T}^{"+xAxisLbl+"}"  
        plot.xTitle = "E_{T}^{"+xAxisLbl+"}"
        plot.xRange = (0.0,150) 
        plot.logx = False
        plot.yRange = (0.7,1.2) 
        plots.append(plot) 
        
        i=4
        hname='EtVsdEtScaleECap'
        for tag in tagsToPlot:
            hist=histStore[tag][hname].Clone()
        #         hist.Rebin(rBin)
            pparams=getDefaultPlotParams(col=i,marker=22)
            pparams['Legend']=legend[tag]
            pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
            pparams['MarkerStyle']=8
            pparams['Options']='HIST pe'
            aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
            aplot.normalize = False ;aplot.scaleFactor =None
            aplot.doFit=False ; aplot.drawFit=False
            aplot.drawLegend = True
            aplot.verbose= False
            plots[-1].addPlot(aplot)
    
        for tag in tagsToPlot:
            plot=getDefaultPlot(prefix=prefix,name='EtVsdEtScale_'+tag+'_'+xAxisLbl,cPars=cPars)   
            plot.legendPosition = (0.45,0.75,0.89,0.90)
            plot.descPosition   = (0.55,0.70)
            plot.desc =[legend[tag]]
            plot.yTitle = "E_{T}^{L1} / E_{T}^{"+xAxisLbl+"}"  
            plot.xTitle = "E_{T}^{"+xAxisLbl+"}"
            plot.xRange = (0.0,150) 
            plot.logx = False
            plot.yRange = (0.7,1.2) 
            plots.append(plot) 
            cc={'EtVsdEtScaleBarrel':ROOT.kRed,'EtVsdEtScaleECap':ROOT.kBlue}
            ll={'EtVsdEtScaleBarrel':'Barrel','EtVsdEtScaleECap':'Endcap'}
            for hname in ['EtVsdEtScaleBarrel','EtVsdEtScaleECap']:
                hist=histStore[tag][hname].Clone()
                pparams=getDefaultPlotParams(col=i,marker=22)
                pparams['Legend']=ll[hname]
                pparams['MarkerColor']= cc[hname]; pparams['LineColor'] = cc[hname] ;i+=1;
                pparams['MarkerStyle']=8
                pparams['Options']='HIST pe'
                aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
                aplot.normalize = False ;aplot.scaleFactor =None
                aplot.doFit=False ; aplot.drawFit=False
                aplot.drawLegend = True
                aplot.verbose= False
                plots[-1].addPlot(aplot)
            
        canvas = []
        for plot in plots:
            canvas.append(plot.plot())
     
def plotEtResolutions(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    
    plot=getDefaultPlot(prefix=prefix,name='resolutionET_Barrel'+varToChose,cPars=cPars) 
    plot.legendPosition = (0.16,0.70,0.66,0.88)
    plot.descPosition   = (0.25,0.65)
    plot.desc =[" "]
    plot.yTitle = "a.u"  
    plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
    plot.xRange = (0.7,1.20) 
    plot.logx = False
    plot.yRange = (0.0,0.14) 
    plots.append(plot) 

    i=2
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hname='EtVsdEtBareResBarrel'
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']='Barrel , '+legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pl'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)


    plot=getDefaultPlot(prefix=prefix,name='resolutionET_ECap',cPars=cPars) 
    plot.legendPosition = (0.2,0.72,0.72,0.90)
    plot.descPosition   = (0.65,0.75)
    plot.desc =[" "]
    plot.yTitle = "a.u"  
    plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
    plot.xRange = (0.62,1.60) 
    plot.logx = False
    plot.yRange = (0.0,0.14) 
    plots.append(plot) 
    
    i=4
    hname='EtVsdEtBareResECap'
    for tag in tagsToPlot: #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']='Endcap , '+legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pl'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())

def plotEtResolutionsGreaterThan20(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    
    plot=getDefaultPlot(prefix=prefix,name='resolutionETGThan20_Barrel',cPars=cPars) 
    plot.legendPosition = (0.16,0.65,0.60,0.80)
    plot.descPosition   = (0.25,0.55)
    plot.desc =["L1 EG "," E_{T}^{"+varToChose+"} > 20 GeV"]
    plot.yTitle = "a.u"  
    plot.xTitle = "E_{T}^{L1} / E_{T}^{"+varToChose+"}"
    plot.xRange = (0.7,1.20) 
    plot.logx = False
    plot.yRange = (0.0,0.24) 
    plots.append(plot) 

    i=2
    for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        print(tag)
        hname='PtGThan20EtVsdEtBareResBarrel'
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']='Barrel , '+legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pl'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)


    plot=getDefaultPlot(prefix=prefix,name='resolutionETGthan20_ECap',cPars=cPars) 
    plot.legendPosition = (0.2,0.72,0.72,0.90)
    plot.descPosition   = (0.65,0.65)
    plot.desc =["L1 EG "," E_{T}^{offline} > 20 GeV"]
    plot.yTitle = "a.u"  
    plot.xTitle = "E_{T}^{L1} / E_{T}^{offline}"
    plot.xRange = (0.62,1.60) 
    plot.logx = False
    plot.yRange = (0.0,0.14) 
    plots.append(plot) 
    
    i=4
    hname='PtGThan20EtVsdEtBareResECap'
    for tag in tagsToPlot:#['run3MC','dataR3Unpaked_postCalib','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']='Endcap , '+legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pl'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())

def plotEtResolutionInPt(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):

    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    plot=getDefaultPlot(prefix=prefix,name='EtVsdEt_highResoRes_Barrel',cPars=cPars) 
    plot.legendPosition = (0.40,0.40,0.90,0.60)
    plot.descPosition   = (0.55,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","Barrel Region"]
    plot.yTitle = "FWHM <E_{T}^{L1} / E_{T}^{"+varToChose+"}>"  
    plot.xTitle = "E_{T}^{"+varToChose+"}"
    plot.xRange = (5.0,100.0) 
    plot.logx = False
    plot.yRange = (0.0,0.35) 
    plots.append(plot) 

    i=2
    hname='EtVsdEt_highResoResBarrel'
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pe'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=False ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)
    
    plot=getDefaultPlot(prefix=prefix,name='EtVsdEt_highResoRes_ECap',cPars=cPars) 
    plot.legendPosition = (0.40,0.40,0.90,0.60)
    plot.descPosition   = (0.55,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","Endcap Region"]
    plot.yTitle = "FWHM <E_{T}^{L1} / E_{T}^{"+varToChose+"}>"  
    plot.xTitle = "E_{T}^{"+varToChose+"}"
    plot.xRange = (5.0,100.0) 
    plot.logx = False
    plot.yRange = (0.0,0.7) 
    plots.append(plot) 

    i=4
    hname='EtVsdEt_highResoResECap'
    for tag in tagsToPlot:# ['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pe'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=False ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())
    


def plotEtResolutionInEta(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    plot=getDefaultPlot(prefix=prefix,name='EtVsdEt_highResoRes_Barrel',cPars=cPars) 
    plot.legendPosition = (0.40,0.40,0.90,0.60)
    plot.descPosition   = (0.55,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","Barrel Region"]
    plot.yTitle = "FWHM <E_{T}^{L1} / E_{T}^{"+varToChose+"}>"  
    plot.xTitle = "E_{T}^{"+varToChose+"}"
    plot.xRange = (5.0,100.0) 
    plot.logx = False
    plot.yRange = (0.0,0.35) 
    plots.append(plot) 

    i=2
    hname='EtVsdEt_highResoResBarrel'
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pe'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=False ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)
    
    plot=getDefaultPlot(prefix=prefix,name='EtVsdEt_highResoRes_ECap',cPars=cPars) 
    plot.legendPosition = (0.40,0.40,0.90,0.60)
    plot.descPosition   = (0.55,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","Endcap Region"]
    plot.yTitle = "FWHM <E_{T}^{L1} / E_{T}^{"+varToChose+"}>"  
    plot.xTitle = "E_{T}^{"+varToChose+"}"
    plot.xRange = (5.0,100.0) 
    plot.logx = False
    plot.yRange = (0.0,0.7) 
    plots.append(plot) 

    i=4
    hname='EtVsdEt_highResoResECap'
    for tag in tagsToPlot:# ['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag] ;i+=1;
        pparams['Options']='HIST pe'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = False ;aplot.scaleFactor =None
        aplot.doFit=False ; aplot.drawFit=False
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())
    


def plotEtaPhiresolutions(histStore,tagsToPlot,colours,legend,era,lumi,prefix,varToChose):
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    plot=getDefaultPlot(prefix=prefix,name='EtaResolution',cPars=cPars) 
    plot.legendPosition = (0.55,0.25,0.90,0.39)
    plot.descPosition   = (0.65,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","|#eta| <2.4"]
    plot.yTitle = "a.u"  
    plot.xTitle = "#eta_{L1}-#eta_{"+varToChose+"}"
    plot.xRange = (-0.1,0.2) 
    plot.logx = False
    plot.yRange = (0.0,0.15) 
    plots.append(plot) 

    i=2
    hname='EtaVsdEtaBareResInclusive'
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerStyle']=8
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
        pparams['Options']='HIST pec'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = True ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=True
        aplot.drawLegend=True
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())

    
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    plot=getDefaultPlot(prefix=prefix,name='EtaResolution',cPars=cPars) 
    plot.legendPosition = (0.55,0.25,0.90,0.39)
    plot.descPosition   = (0.65,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","|#eta| <2.4"]
    plot.yTitle = "a.u"  
    plot.xTitle = "#eta_{L1}-#eta_{"+varToChose+"}"
    plot.xRange = (-0.1,0.2) 
    plot.logx = False
    plot.yRange = (0.0,0.15) 
    plots.append(plot) 

    i=2
    hname='EtaVsdEtaBareResInclusive'
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerStyle']=8
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
        pparams['Options']='HIST pec'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = True ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=True
        aplot.drawLegend=True
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())

        
    plots=[]
    cPars=getCanvasParams(baseCanvas,era,lumi)
    rBin=1
    plot=getDefaultPlot(prefix=prefix,name='PhiResolution',cPars=cPars) 
    plot.legendPosition = (0.55,0.25,0.90,0.39)
    plot.descPosition   = (0.65,0.75)
    plot.desc =["L1 EG","Loose Offline Electron","|#eta| <2.4"]
    plot.yTitle = "a.u"  
    plot.xTitle = "#phi_{L1}-#phi_{offline}"
    plot.xRange = (-0.15,0.25) 
    plot.logx = False
    plot.yRange = (0.0,0.15) 
    plots.append(plot) 

    i=2
    hname='PhiVsdPhiBareResInclusive'
    for tag in tagsToPlot : #['run3MC','dataR3Unpaked_postCalib','data2018','run3MCWithGenEt','run3MCWithGenEtVsGenEt']:
        hist=histStore[tag][hname].Clone()
    #         hist.Rebin(rBin)
        pparams=getDefaultPlotParams(col=i,marker=22)
        pparams['Legend']=legend[tag]
        pparams['MarkerStyle']=8
        pparams['MarkerColor']= colours[tag]; pparams['LineColor'] = colours[tag]
        pparams['Options']='HIST pec'
        aplot = Plot(Name=hist.GetName(), Histo=hist,**pparams)
        aplot.normalize = True ;aplot.scaleFactor =None
        aplot.doFit=True ; aplot.drawFit=True
        aplot.drawLegend=True
        aplot.verbose= False
        plots[-1].addPlot(aplot)

    canvas = []
    for plot in plots:
        canvas.append(plot.plot())
                  
