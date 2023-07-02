#!/usr/bin/env python
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os

def GetWeights(pt, is2016 = False):

    params = [-2.02274e-01, 1.09734e-04, -1.30088e-07, 5.83494e+01, 1.96252e+02]
    params2016 = [1.04554e+00, 5.19012e-02, -1.72927e+00, 2.57113e-03]

    if is2016:
        weight = params2016[0] + params2016[1] * math.tanh(params2016[2] + params2016[3]*pt)
    else:
        weight = math.exp(params[0] + params[1]*pt + params[2]*pt*pt + (params[3]/(pt + params[4])))

    return weight

class topPtReweightProduce(Module):
    def __init__(self, year):
        self.year = year
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch('TopPtWeight','F', lenVar='1')
        self.out.branch('TopPtWeight_up','F', lenVar='1')
        self.out.branch('TopPtWeight_down','F', lenVar='1')

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        # Using HH->bbWW example: https://github.com/FlorianBury/HHbbWWAnalysis/blob/aaad8763c2fcadfff45d6c008de6df3940ace1ae/BaseHHtobbWW.py#L1037
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m -> do not use when ttbar is background
        # Correct : https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
        #       -> Weight formula           : slide 2
        #       -> top/antitop SF           : slide 12 bottom left 
        #       -> 2016 CUET extra addition : slide 13 bottom left 

        topPt = 0
        topPtWeight = 1.0
        antiopPt = 0
        antitopPtWeight = 1.0
        eventWeight = 1.0

        genParticles = Collection(event, 'GenPart')
        for igen in range(0, event.nGenPart):
            if genParticles[igen].pdgId == 6 and genParticles[igen].statusFlags == 13:
                topPt = genParticles[igen].pt

            if genParticles[igen].pdgId == -6 and genParticles[igen].statusFlags == 13:
                antitopPt = genParticles[igen].pt

            if topPt > 0 and antitopPt > 0:
                break

        if topPt > 0 and antitopPt > 0:
            topPtWeight = GetWeights(topPt, is2016=False)
            antitopPtWeight = GetWeights(antitopPt, is2016=False)
        eventWeight = math.sqrt(topPtWeight*antitopPtWeight)

        if '2016' in self.year:
            topPtWeight = GetWeights(topPt, is2016=True)
            antitopPtWeight = GetWeights(antitopPt, is2016=True)
            eventWeight *= math.sqrt(topPtWeight*antitopPtWeight)

        eventWeight_up = eventWeight*eventWeight
        eventWeight_down = 1.0

        self.out.fillBranch('TopPtWeight', eventWeight)
        self.out.fillBranch('TopPtWeight_up', eventWeight_up)
        self.out.fillBranch('TopPtWeight_down', eventWeight_down)

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

topPtReweightModule2016 = lambda : topPtReweightProduce(year='2016')
topPtReweightModule2017 = lambda : topPtReweightProduce(year='2017')
topPtReweightModule2018 = lambda : topPtReweightProduce(year='2018')