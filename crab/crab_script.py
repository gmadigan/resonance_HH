import os
import sys
import optparse
import ROOT
import re

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.modules.common.countHistogramsModule import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.eleIDSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.highPtMuonScaleResProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.topPtReweightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonIDISOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.muonRECOSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.modules.LQProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.others.for_jme.jetmetHelperRun2 import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.others.for_btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.others.for_pileup.puWeightProducer import *
from PhysicsTools.NanoAODTools.postprocessing.analysis.others.for_prefiring.PrefireCorr import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.puWeightProducer import *
#from PhysicsTools.NanoAODTools.postprocessing.modules.common.PrefireCorr import *
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles, runsAndLumis
### main python file to run ###

def main():

  usage = 'usage: %prog [options]'
  parser = optparse.OptionParser(usage)
  parser.add_option('--year', dest='year', help='which year sample', default='2018', type='string')
  parser.add_option('-m', dest='ismc', help='to apply sf correction or not', default=True, action='store_true')
  parser.add_option('-d', dest='ismc', help='to apply sf correction or not', action='store_false')
  (opt, args) = parser.parse_args()

  # full preselection for mumujj channel
  # do not apply here; will lose some events needed for systematic up/down variations
  # preselection = "(lead_muon_pt > 53)*(trail_muon_pt > 53)*(-2.4 < lead_muon_eta < 2.4)*(-2.4 < trail_muon_eta < 2.4)*(lead_jet_pt > 50)*(trail_jet_pt > 50)*(-2.4 < lead_jet_eta < 2.4)*(-2.4 < trail_jet_eta < 2.4)*(DR_mumu < 0.3)*(M_mumu < 50)*(ST < 300)"
  
  # light skim on the leading-pT muon
  preselection = "lead_muon_pt > 20"

  # Modules add branches to NanoAOD ntuples under "events" branch
  # For MC, mostly includes corrections to match MC to data; entirely year-dependant
  #
  # countHistogramModule()
  # Adds a histogram with the original (generated) number of events; used for normalizing MC (on same level as "events" branch; not under "events" branch)
  #
  # puWeight_[year]()
  # Provides weights to correct pileup mismodeling
  #
  # PrefCorr[year]() 
  # Provides weights that correct the L1 prefire rates for the ECAL and Muon system in 2016 and 2017 only
  #
  # muonIDISOSF[year]()
  # Provides scale factors to match the muon identification and isolation efficiencies in data
  #
  # highPtMuonScaleRes[year]()
  # Provides corrected muon pT, eta, and phi measurements to match Muon Energy Scale (MES) and Muon Energy Resolution (MER) seen in data
  # Follows reccomendations for high-pT muons:
  # Uses muon pT measurements with the TuneP algorithm, "corrected" from Particle Flow (PF);
  # Muon pT- and eta-dependent Rochester scale/resolution corrections/smearing for muons with pT < 200 GeV;
  # An "extra" muon p- (3-momentum) and eta-dependant resolution smearing is applied to all muons; not pT-dependant, but changes in p are translated to pT, eta, and phi measurements
  #
  # eleRECOSF[year]()
  # Provides scale factors to match the electron reconstruction efficieny in data 
  #
  # eleIDSF[year]() 
  # Provides scale factors to match the electron identification efficieny in data 
  #
  # jmeCorrections_UL[year]MC()
  # Provides systematic variations for Jet Energy Corrections (JECs) on AK4PFCHS jets;
  # Jet Energy Resolution (JER) smearing also applied to AK4PFCHS jets and propogated to jet kinematics
  # JER systematic variations provided
  # Systematic uncertainties of JECs for jets in HEM failure region fixed
  #
  # btagSF[year]()
  # Provides scale factors to match the b-tag efficiency (using the DeepFlavorB--aka, DeepJet---algorithm) measured in data
  #
  # topPtReweighting[year]()
  # Provides weights that reweight top pT distributions in tt-bar MC only to correct mismodeling
  #
  # LQ[year]() 
  # Performs the object selection and computes all the kinematic quantities needed for the LQ analysis in mumujj, mumuj, munujj, and munuj channels;
  # Muons are high-pT ID muons; leading and sub-leading (trailing) pT muons kept for analyis;
  # Electrons are HEEP ID electrons; leading and sub-leading (trailing) pT electrons kept for analyis;
  # Jets are AK4PFCHS jets; leading and sub-leading (trailing) pT jets kept for analysis;
  # Missing Transverse Energy (MET) uses PF algorithm with JECs propogated;


  if opt.ismc:
    if opt.year == "2016a":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[countHistogramsModule(),puWeight_2016_preAPV(),PrefCorr2016(),muonIDISOSF2016apv(),muonRECOSF2016apv(),highPtMuonScaleRes2016a(),eleRECOSF2016apv(),eleIDSF2016apv(),jmeCorrections_UL2016APVMC(),btagSF2016ULapv(),topPtReweighting2016(),LQ2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016b":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[countHistogramsModule(),puWeight_2016_postAPV(),PrefCorr2016(),muonIDISOSF2016(),muonRECOSF2016(),highPtMuonScaleRes2016b(),eleRECOSF2016(),eleIDSF2016(),jmeCorrections_UL2016MC(),btagSF2016UL(),topPtReweighting2016(),LQ2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[countHistogramsModule(),puWeight_2017(),PrefCorr2017(),muonIDISOSF2017(),muonRECOSF2017(),highPtMuonScaleRes2017(),eleRECOSF2017(),eleIDSF2017(), jmeCorrections_UL2017MC(),btagSF2017UL(),topPtReweighting2017(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[countHistogramsModule(),puWeight_2018(),muonIDISOSF2018(),muonRECOSF2018(),highPtMuonScaleRes2018(),eleRECOSF2018(),eleIDSF2018(),jmeCorrections_UL2018MC(), btagSF2018UL(),topPtReweighting2018(),LQ2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")


# Sequence for data
  if not (opt.ismc):
    if opt.year == "2016b":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016a(),jmeCorrections_UL2016B(),LQ2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016c":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016a(),jmeCorrections_UL2016C(),LQ2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016d":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016a(),jmeCorrections_UL2016D(),LQ2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016e":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016a(),jmeCorrections_UL2016E(),LQ2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016f_apv":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016a(),jmeCorrections_UL2016APVF(),LQ2016apv()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016f":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016b(),jmeCorrections_UL2016F(),LQ2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016g":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016b(),jmeCorrections_UL2016G(),LQ2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2016h":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2016b(),jmeCorrections_UL2016H(),LQ2016()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017b":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2017(),jmeCorrections_UL2017B(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017c":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2017(),jmeCorrections_UL2017C(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017d":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2017(),jmeCorrections_UL2017D(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017e":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2017(),jmeCorrections_UL2017E(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2017f":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2017(),jmeCorrections_UL2017F(),LQ2017()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018a":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2018(),jmeCorrections_UL2018A(),LQ2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018b":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2018(),jmeCorrections_UL2018B(),LQ2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018c":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2018(),jmeCorrections_UL2018C(),LQ2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
    if opt.year == "2018d":
      p = PostProcessor(".", inputFiles(), cut=preselection, modules=[highPtMuonScaleRes2018(),jmeCorrections_UL2018D(),LQ2018()], provenance=True,fwkJobReport=True, jsonInput=runsAndLumis(),outputbranchsel="keep_and_drop.txt")
  p.run()

if __name__ == "__main__":
    sys.exit(main())
