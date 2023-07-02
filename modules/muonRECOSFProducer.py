import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import json

class muonRECOSFProducer(Module):
  def __init__( self , year ):
    self.year = year
    self.reco_input = "muon_RECO_SF.json"
    self.SF_location_path = "%s/src/PhysicsTools/NanoAODTools/python/postprocessing/analysis/data/year%s/" %(os.environ['CMSSW_BASE'], self.year)
#    print 'SF location:', self.SF_location_path

  def beginJob(self):
    print 'begin to set Muon RECO SF --->>>'
    print 'start to open SF root file --->>>'
    self.file_json = open(self.SF_location_path+self.reco_input)
    self.reco_dict = json.load(self.file_json)
    print 'open SF files successfully --->>>'

  def endJob(self):
    print 'close SF root file --->>>'
    self.file_json.close()
    print 'finish setting Muon RECO and SF --->>>'
    
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch('Muon_Reco_LooseRelTkIsoandHighPtIDandIPCut_SF','F', lenVar='nMuon')
    self.out.branch('Muon_Reco_LooseRelTkIsoandHighPtIDandIPCut_SFerr','F', lenVar='nMuon')  
  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    muons = Collection(event, "Muon")
    if not (len(muons)>0): pass
    Muon_RECO_SF = []
    Muon_RECO_SFerr = []
    
    for imu in range(0, len(event.nMuons)):

      # Muon momentum should be TuneP for High-pT ID muons
      # Not necessary to use Rochester + resolution corrected momentum
      # From https://indico.cern.ch/event/835143/contributions/3499965/attachments/1880865/3098909/MuonPOG_July2019.pdf:
      # "No bias from the momentum resolutions is seen in the reconstruction efficiency"

      muonVec = TLorentzVector()
      muonVec.SetPtEtaPhiM(muons[imu].pt * muons[imu].tunepRelPt, muons[imu].eta, muons[imu].phi, muons[imu].mass)
      muonAbsEta = abs(muonVec.Eta())
      muonP = muonVec.P()

      for abseta in self.reco_dict:
        etaBin = abseta.split(':')[-1]
        etaBinLow = float(etaBin.strip('[]').split(',')[0])
        etaBinHigh = float(etaBin.strip('[]').split(',')[-1])
        if 'abseta' in abseta and etaBinLow < muonAbsEta < etaBinHigh: 
          for p in self.reco_dict[abseta]:
            pBin = p.split(':')[-1]
            pBinLow = float(pBin.strip('[]').split(',')[0])
            pBinHigh = float(pBin.strip('[]').split(',')[-1])
            if 'p' in p and pBinLow < muonP < pBinHigh:
              Muon_RECO_SF.append(self.reco_dict[abseta][p]["SF"])
              Muon_RECO_SFerr.append(self.reco_dict[abseta][p]["Error"])

    self.out.fillBranch('Muon_Reco_LooseRelTkIsoandHighPtIDandIPCut_SF', Muon_RECO_SF)
    self.out.fillBranch('Muon_Reco_LooseRelTkIsoandHighPtIDandIPCut_SFerr', Muon_RECO_SFerr)

    return True

muonRECOSF2016apv = lambda: muonRECOSFProducer("2016apv")
muonRECOSF2016 = lambda: muonRECOSFProducer("2016")
muonRECOSF2017 = lambda: muonRECOSFProducer("2017")
muonRECOSF2018 = lambda: muonRECOSFProducer("2018")
