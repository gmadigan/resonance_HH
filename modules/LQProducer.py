import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign

def SortCollections(listOfCollections):
  # Takes a list of collections (e.g., pT, eta, phi) and sorts the first collection by value (high to low).
  # Simultaniously, the other collections are reordered identically to the first. 
  # Returns the reordered collections as a list.
  # Make sure Pt cocktails the first collection in ListOfCollections!

  zippedCollections = zip(*listOfCollections) # zip unpacked (*) collections together into a list of tuples
  sortedCollections = sorted(zippedCollections, reverse=True) # reorder all collections simultaniously by pT
  unzippedCollections = zip(*sortedCollections) # unzip collections; still a list of tuples
  returnListOfCollections = map(list, unzippedCollections) # Map list of tuples back into separate lists and unpack as collections

  return returnListOfCollections

class LQProducer(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
  
    self.out.branch("lead_muon_pt", "F")
    self.out.branch("lead_muon_eta", "F")
    self.out.branch("lead_muon_phi", "F")
    self.out.branch("lead_muon_mass", "F")
    self.out.branch("lead_muon_ind", "I")

    self.out.branch("lead_muon_RecoSF", "F")
    self.out.branch("lead_muon_RecoSF_up", "F")
    self.out.branch("lead_muon_RecoSF_down", "F")
    self.out.branch("lead_muon_IsoSF", "F")
    self.out.branch("lead_muon_IsoSF_up", "F")
    self.out.branch("lead_muon_IsoSF_down", "F")
    self.out.branch("lead_muon_IDSF", "F")
    self.out.branch("lead_muon_IDSF_up", "F")
    self.out.branch("lead_muon_IDSF_down", "F")
    self.out.branch("lead_muon_HLTSF", "F")
    self.out.branch("lead_muon_HLTSF_up", "F")
    self.out.branch("lead_muon_HLTSF_down", "F")

    self.out.branch("trail_muon_pt", "F")
    self.out.branch("trail_muon_eta", "F")
    self.out.branch("trail_muon_phi", "F")
    self.out.branch("trail_muon_mass", "F")
    self.out.branch("trail_muon_ind", "I")
    
    self.out.branch("trail_muon_RecoSF", "F")
    self.out.branch("trail_muon_RecoSF_up", "F")
    self.out.branch("trail_muon_RecoSF_down", "F")
    self.out.branch("trail_muon_IsoSF", "F")
    self.out.branch("trail_muon_IsoSF_up", "F")
    self.out.branch("trail_muon_IsoSF_down", "F")
    self.out.branch("trail_muon_IDSF", "F")
    self.out.branch("trail_muon_IDSF_up", "F")
    self.out.branch("trail_muon_IDSF_down", "F")
    self.out.branch("trail_muon_HLTSF", "F")
    self.out.branch("trail_muon_HLTSF_up", "F")
    self.out.branch("trail_muon_HLTSF_down", "F")

    self.out.branch("ngood_muons", "I")

    self.out.branch("lead_ele_pt", "F")
    self.out.branch("lead_ele_eta", "F")
    self.out.branch("lead_ele_phi", "F")
    self.out.branch("lead_ele_mass", "F")
    self.out.branch("lead_ele_ind", "I")

    self.out.branch("trail_ele_pt", "F")
    self.out.branch("trail_ele_eta", "F")
    self.out.branch("trail_ele_phi", "F")
    self.out.branch("trail_ele_mass", "F")
    self.out.branch("trail_ele_ind", "I")

    self.out.branch("ngood_electrons", "I")

    self.out.branch("lead_jet_pt", "F")
    self.out.branch("lead_jet_eta", "F")
    self.out.branch("lead_jet_phi", "F")
    self.out.branch("lead_jet_mass", "F")
    self.out.branch("lead_jet_ind", "I")
    self.out.branch("lead_jet_btag", "F")
    self.out.branch("lead_jet_btagSF", "F")
    self.out.branch("lead_jet_btagSF_up", "F")
    self.out.branch("lead_jet_btagSF_down", "F")

    self.out.branch("trail_jet_pt", "F")
    self.out.branch("trail_jet_eta", "F")
    self.out.branch("trail_jet_phi", "F")
    self.out.branch("trail_jet_mass", "F")
    self.out.branch("trail_jet_ind", "I")
    self.out.branch("trail_jet_btag", "F")
    self.out.branch("trail_jet_btagSF", "F")
    self.out.branch("trail_jet_btagSF_up", "F")
    self.out.branch("trail_jet_btagSF_down", "F")

    self.out.branch("ngood_jets", "I")
    self.out.branch("ngood_bjets", "I")

    self.out.branch("M_mumu", "F")
    self.out.branch("M_mumujj", "F")
    self.out.branch("ST_mumujj", "F")
    self.out.branch("M_mumuj", "F")
    self.out.branch("ST_mumuj", "F")

    self.out.branch("M_ee", "F")
    self.out.branch("M_eejj", "F")
    self.out.branch("ST_eejj", "F")
    self.out.branch("M_eej", "F")
    self.out.branch("ST_eej", "F")

    self.out.branch("DeltaR_mumu", "F")
    self.out.branch("DeltaR_jj", "F")
    self.out.branch("DeltaR_leadMuon_leadJet", "F")
    self.out.branch("DeltaR_leadMuon_trailJet", "F")
    self.out.branch("DeltaR_trailMuon_leadJet", "F")
    self.out.branch("DeltaR_trailMuon_trailJet", "F")
    self.out.branch("DeltaR_mumu_trailJet", "F")

    self.out.branch("DeltaPhi_mumu", "F")
    self.out.branch("DeltaPhi_jj", "F")
    self.out.branch("DeltaPhi_leadMuon_leadJet", "F")
    self.out.branch("DeltaPhi_leadMuon_trailJet", "F")
    self.out.branch("DeltaPhi_trailMuon_leadJet", "F")
    self.out.branch("DeltaPhi_trailMuon_trailJet", "F")
    self.out.branch("DeltaPhi_mumu_trailJet", "F")

    self.out.branch("DeltaR_ee", "F")
    self.out.branch("DeltaR_leadEle_leadJet", "F")
    self.out.branch("DeltaR_leadEle_trailJet", "F")
    self.out.branch("DeltaR_trailEle_leadJet", "F")
    self.out.branch("DeltaR_trailEle_trailJet", "F")
    self.out.branch("DeltaR_ee_trailJet", "F")

    self.out.branch("DeltaPhi_ee", "F")
    self.out.branch("DeltaPhi_leadEle_leadJet", "F")
    self.out.branch("DeltaPhi_leadEle_trailJet", "F")
    self.out.branch("DeltaPhi_trailEle_leadJet", "F")
    self.out.branch("DeltaPhi_trailEle_trailJet", "F")
    self.out.branch("DeltaPhi_ee_trailJet", "F")

    self.out.branch("M_mumujj_lead", "F")
    self.out.branch("M_mumujj_trail", "F")
    self.out.branch("M_mumuj_lead", "F")
    self.out.branch("M_mumuj_trail", "F")

    self.out.branch("M_eejj_lead", "F")
    self.out.branch("M_eejj_trail", "F")
    self.out.branch("M_eej_lead", "F")
    self.out.branch("M_eej_trail", "F")

    self.out.branch("met_pt", "F")
    self.out.branch("met_phi", "F")
    self.out.branch("met_filter", "I")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):

    ################################################
    ##### single LQ analysis requirements here #####
    ################################################

    emptyVector = TLorentzVector()
    emptyVector.SetPtEtaPhiM(0,0,0,0)

    ################################################
    #####            MET calculation           #####
    ################################################

    # Define MET object
    # Will propogate muon, electron, and jet corrections to MET later in function
    met = TLorentzVector()

    # Apply correction to EE noise issue in 2017 data and MC
    # Type-1 JER smearing applied and propogated to MET only in MC
    if self.isMC:
      if self.year == "2017": 
        met.SetPtEtaPhiM(T.METFixEE2017_T1Smear_pt,0,T.METFixEE2017_T1Smear_phi,0)
      else: 
        met.SetPtEtaPhiM(T.MET_T1Smear_pt,0,T.MET_T1Smear_phi,0)
    else:
      if self.year == "2017": 
        met.SetPtEtaPhiM(T.METFixEE2017_T1_pt,0,T.METFixEE2017_T1_phi,0)
      else: 
        met.SetPtEtaPhiM(T.MET_T1_pt,0,T.MET_T1_phi,0)


    ################################################
    #####            Muon selection            #####
    ################################################

    # Require muons to pass high-pT ID and loose track-based relative isolation
    # Use TuneP algorithm for momentum with Rochester scale/resolution corrections/smearing applied
    # For high-pT muons, additional resolution smearing is applied
    # Select the two highest pT (leading- and trailing-pT) muons for analysis
    # Create branches for kinematics, scale factors, and up and down systematic variations for these two muons
    # If less than two muons pass the ID and Iso requirements, create zero-valued placeholders
    # Events with placeholders will be rejected with preselection requirements on muons
    #
    # No kinematic cuts are applied in this module; cuts can kill events needed in the estimation of the systematic uncertainties,
    # e.g., the muon resolution systematic variations will shift the pT and eta measurements up and down, 
    # applying a pT or eta cut could remove the down-varied events leading to an underestimation of the uncertainty (done at the histogram level)


    muons = Collection(event, 'Muon')

    highPtMuons = []
    highPtMuonPt = []
    highPtMuonEta = []
    highPtMuonPhi = []
    highPtMuonMass = []
    highPtMuons_ind = []
    nMuons = 0

    for imu in range(0, event.nMuon):
      # Select muons passing the high-pT ID working point; requires TuneP momentum (automatically used in the "highPtResSmeared" corrected quantities)
      # Track-based relative isolation; loose working point is 0.05 and corresponds to a 95% efficiency; measured within DeltaR < 0.3
      if muons[imu].tkRelIso > 0.05 or not muons[imu].highPtId > 1: 
        continue
      else: 
        tmpMuon = TLorentzVector()
        tmpMuon.SetPtEtaPhiM(muons[imu].highPtResSmeared_pt, muons[imu].highPtResSmeared_eta, muons[imu].highPtResSmeared_phi, muons[imu].mass)
        highPtMuons.append(tmpMuon)
        highPtMuonPt.append(tmpMuon.Pt())
        highPtMuonEta.append(tmpMuon.Eta())
        highPtMuonPhi.append(tmpMuon.Phi())
        highPtMuonMass.append(tmpMuon.M())
        highPtMuons_ind.append(imu)

        if self.isMC:
          # Get muon efficiency scale factors and systematic variations
          muonRecoSFs.append(muons[imu].Reco_LooseRelTkIsoandHighPtIDandIPCut_SF)
          muonRecoSFs_up.append(muons[imu].Reco_LooseRelTkIsoandHighPtIDandIPCut_SF + muons[imu].Reco_LooseRelTkIsoandHighPtIDandIPCut_SFerr)
          muonRecoSFs_down.append(muons[imu].Reco_LooseRelTkIsoandHighPtIDandIPCut_SF - muons[imu].Reco_LooseRelTkIsoandHighPtIDandIPCut_SFerr)
          muonIsoSFs.append(muons[imu].TightRelTkIso_HighPtIDandIPCut_SF)
          muonIsoSFs_up.append(muons[imu].TightRelTkIso_HighPtIDandIPCut_SF + muons[imu].TightRelTkIso_HighPtIDandIPCut_SFerr)
          muonIsoSFs_down.append(muons[imu].TightRelTkIso_HighPtIDandIPCut_SF - muons[imu].TightRelTkIso_HighPtIDandIPCut_SFerr)
          muonIDSFs.append(muons[imu].CutBased_HighPtID_SF)
          muonIDSFs_up.append(muons[imu].CutBased_HighPtID_SF + muons[imu].CutBased_HighPtID_SFerr)
          muonIDSFs_down.append(muons[imu].CutBased_HighPtID_SF - muons[imu].CutBased_HighPtID_SFerr)
          muonHLTSFs.append(1.0) # Temporary while trigger module in dev
          muonHLTSFs_up.append(1.0) # Temporary while trigger module in dev
          muonHLTSFs_down.append(1.0) # Temporary while trigger module in dev
        else:
          muonRecoSFs.append(1.0)
          muonRecoSFs_up.append(1.0)
          muonRecoSFs_down.append(1.0)
          muonIsoSFs.append(1.0)
          muonIsoSFs_up.append(1.0)
          muonIsoSFs_down.append(1.0)
          muonIDSFs.append(1.0)
          muonIDSFs_up.append(1.0)
          muonIDSFs_down.append(1.0)
          muonHLTSFs.append(1.0)
          muonHLTSFs_up.append(1.0)
          muonHLTSFs_down.append(1.0)

        # propogate changes to MET
        met += tmpMuon.Pt()
        met -= muons[imu].pt
        met_phi += tmpMuon.Phi()
        met_phi -= muons[imu].phi

        nMuons += 1

    while len(highPtMuons) < 2: 
      highPtMuons.append(emptyVector)
      highPtMuonPt.append(-99.9)
      highPtMuonEta.append(-99.9)
      highPtMuonPhi.append(-99.9)
      highPtMuonMass.append(-99.9)
      highPtMuons_ind.append(-99)
      muonRecoSFs.append(0.0)
      muonRecoSFs_up.append(0.0)
      muonRecoSFs_down.append(0.0)
      muonIsoSFs.append(0.0)
      muonIsoSFs_up.append(0.0)
      muonIsoSFs_down.append(0.0)
      muonIDSFs.append(0.0)
      muonIDSFs_up.append(0.0)
      muonIDSFs_down.append(0.0)
      muonHLTSFs.append(0.0)
      muonHLTSFs_up.append(0.0)
      muonHLTSFs_down.append(0.0)

    # Sort the muon collection by pT to get the two highest-pT muons
    highPtMuonPt, highPtMuonEta, highPtMuonPhi, highPtMuonMass, highPtMuons_ind, muonRecoSFs, muonRecoSFs_up, muonRecoSFs_down, muonIsoSFs, muonIsoSFs_up, muonIsoSFs_down, muonIDSFs, muonIDSFs_up, muonIDSFs_down = SortCollectionsByPt([highPtMuonPt, highPtMuonEta, highPtMuonPhi, highPtMuonMass, highPtMuons_ind, muonRecoSFs, muonRecoSFs_up, muonRecoSFs_down, muonIsoSFs, muonIsoSFs_up, muonIsoSFs_down, muonIDSFs, muonIDSFs_up, muonIDSFs_down])

    mu1_pt = highPtMuonPt[0]
    mu1_eta = highPtMuonEta[0]
    mu1_phi = highPtMuonPhi[0]
    mu1_mass = highPtMuonMass[0]
    mu1_ind = highPtMuons_id[0]
    mu1_recoSF = muonRecoSFs[0]
    mu1_recoSF_up = muonRecoSFs_up[0]
    mu1_recoSF_down = muonRecoSFs_down[0]
    mu1_IsoSF = muonIsoSFs[0]
    mu1_IsoSF_up = muonIsoSFs_up[0]
    mu1_IsoSF_down = muonIsoSFs_down[0]
    mu1_IDSF = muonIDSFs[0]
    mu1_IDSF_up = muonIDSFs_up[0]
    mu1_IDSF_down = muonIDSFs_down[0]
    mu1_HLTSF = muonHLTSFs[0]
    mu1_HLTSF_up = muonHLTSFs_up[0]
    mu1_HLTSF_down = muonHLTSFs_down[0]

    mu2_pt = highPtMuonPt[1]
    mu2_eta = highPtMuonEta[1]
    mu2_phi = highPtMuonPhi[1]
    mu2_mass = highPtMuonMass[1]
    mu2_ind = highPtMuons_id[1]
    mu2_recoSF = muonRecoSFs[1]
    mu2_recoSF_up = muonRecoSFs_up[1]
    mu2_recoSF_down = muonRecoSFs_down[1]
    mu2_IsoSF = muonIsoSFs[1]
    mu2_IsoSF_up = muonIsoSFs_up[1]
    mu2_IsoSF_down = muonIsoSFs_down[1]
    mu2_IDSF = muonIDSFs[1]
    mu2_IDSF_up = muonIDSFs_up[1]
    mu2_IDSF_down = muonIDSFs_down[1]
    mu2_HLTSF = muonHLTSFs[1]
    mu2_HLTSF_up = muonHLTSFs_up[1]
    mu2_HLTSF_down = muonHLTSFs_down[1]

    self.out.fillBranch("lead_muon_pt",mu1_pt)
    self.out.fillBranch("lead_muon_eta",mu1_eta)
    self.out.fillBranch("lead_muon_phi",mu1_phi)
    self.out.fillBranch("lead_muon_mass",mu1_mass)
    self.out.fillBranch("lead_muon_ind",mu1_ind)
    self.out.fillBranch("lead_muon_RecoSF",mu1_recoSF)
    self.out.fillBranch("lead_muon_RecoSF_up",mu1_recoSF_up)
    self.out.fillBranch("lead_muon_RecoSF_down",mu1_recoSF_down)
    self.out.fillBranch("lead_muon_IsoSF",mu1_IsoSF)
    self.out.fillBranch("lead_muon_IsoSF_up",mu1_IsoSF_up)
    self.out.fillBranch("lead_muon_IsoSF_down",mu1_IsoSF_down)
    self.out.fillBranch("lead_muon_IDSF",mu1_IDSF)
    self.out.fillBranch("lead_muon_IDSF_up",mu1_IDSF_up)
    self.out.fillBranch("lead_muon_IDSF_down",mu1_IDSF_down)
    self.out.fillBranch("lead_muon_HLTSF",mu1_HLTSF)
    self.out.fillBranch("lead_muon_HLTSF_up",mu1_HLTSF_up)
    self.out.fillBranch("lead_muon_HLTSF_down",mu1_HLTSF_down)

    self.out.fillBranch("trail_muon_pt",mu2_pt)
    self.out.fillBranch("trail_muon_eta",mu2_eta)
    self.out.fillBranch("trail_muon_phi",mu2_phi)
    self.out.fillBranch("trail_muon_mass",mu2_mass)
    self.out.fillBranch("trail_muon_ind",mu2_ind)
    self.out.fillBranch("trail_muon_RecoSF",mu2_recoSF)
    self.out.fillBranch("trail_muon_RecoSF_up",mu2_recoSF_up)
    self.out.fillBranch("trail_muon_RecoSF_down",mu2_recoSF_down)
    self.out.fillBranch("trail_muon_IsoSF",mu2_IsoSF)
    self.out.fillBranch("trail_muon_IsoSF_up",mu2_IsoSF_up)
    self.out.fillBranch("trail_muon_IsoSF_down",mu2_IsoSF_down)
    self.out.fillBranch("trail_muon_IDSF",mu2_IDSF)
    self.out.fillBranch("trail_muon_IDSF_up",mu2_IDSF_up)
    self.out.fillBranch("trail_muon_IDSF_down",mu2_IDSF_down)
    self.out.fillBranch("trail_muon_HLTSF",mu2_HLTSF)
    self.out.fillBranch("trail_muon_HLTSF_up",mu2_HLTSF_up)
    self.out.fillBranch("trail_muon_HLTSF_down",mu2_HLTSF_down)

    self.out.fillBranch("ngood_muons", nMuons)

    ################################################
    #####          Electron selection          #####
    ################################################

    # Require electrons to pass the HEEP ID
    electrons = Collection(event, 'Electron')

    heepElectrons = []
    heepElectrons_ind = []
    nEle = 0

    for iele in range(0, event.nElectron):
      if not electrons[iele].cutBased_HEEP: 
        continue
      else: 
        tmpElectron = TLorentzVector()
        tmpElectron.SetPtEtaPhiM(electrons[iele].pt, electrons[iele].eta, electrons[iele].phi, electrons[iele].mass)
        heepElectrons.append(tmpElectron)
        heepElectrons_ind.append(iele)

        # propogate changes to MET
        met += tmpElectron.Pt()
        met -= electrons[iele].pt
        met_phi += tmpElectron.Phi()
        met_phi -= electrons[iele].phi

        nEle += 1

    while len(heepElectrons) < 2: 
      heepElectrons.append(emptyVector)
      heepElectrons_ind.append(-99)

    ele1_pt = heepElectrons[0].Pt()
    ele1_eta = heepElectrons[0].Eta()
    ele1_phi = heepElectrons[0].Phi()
    ele1_mass = heepElectrons[0].M()
    ele1_ind = heepElectrons_ind[0]

    ele2_pt = heepElectrons[1].Pt()
    ele2_eta = heepElectrons[1].Eta()
    ele2_phi = heepElectrons[1].Phi()
    ele2_mass = heepElectrons[1].M()
    ele2_ind = heepElectrons_ind[1]

    m_ee = (heepElectrons[0]+heepElectrons[1]).M()

    self.out.fillBranch("lead_ele_pt",ele1_pt)
    self.out.fillBranch("lead_ele_eta",ele1_eta)
    self.out.fillBranch("lead_ele_phi",ele1_phi)
    self.out.fillBranch("lead_ele_mass",ele1_mass)
    self.out.fillBranch("lead_ele_ind",ele1_ind)

    self.out.fillBranch("trail_ele_pt",ele2_pt)
    self.out.fillBranch("trail_ele_eta",ele2_eta)
    self.out.fillBranch("trail_ele_phi",ele2_phi)
    self.out.fillBranch("trail_ele_mass",ele2_mass)
    self.out.fillBranch("trail_ele_ind",ele2_ind)

    self.out.fillBranch("ngood_electrons", nEle)

    ################################################
    #####            Jet selection             #####
    ################################################

    # Require jets to pass the tight and tightLepVeto ID
    jets = Collection(event, "Jet")

    tightJets = []
    tightJets_ind = []
    tightJets_mediumDeepFlavB_id = []
    tightJets_btag = []
    tightJets_btagSF = []
    tightJets_btagSF_up = []
    tightJets_btagSF_down = []

    medium_Bcut = 0.2598
    if self.year=="2016":
      medium_Bcut = 0.2489
    if self.year=="2017":
      medium_Bcut = 0.3040
    if self.year=="2018":
      medium_Bcut = 0.2783

    for ijet in range(0, event.nJet):
      if not (jets[ijet].jetId == 6):
        continue
      else:
        tmpJet = TLorentzVector()
        tmpJet.SetPtEtaPhiM(jets[ijet].pt, jets[ijet].eta, jets[ijet].phi, jets[ijet].mass)
        tightJets.append(tmpJet)
        tightJets_ind.append(ijet)

        # Get DeepFlavorB b-tag scores and scale factors for each tight ID jet    
        tightJets_btag.append(jets[ijet].btagDeepFlavB)
        if jets[ijet].btagDeepFlavB > medium_Bcut:
          tightJets_mediumDeepFlavB_ind.append(ijet)

        if self.isMC:
          tightJets_btagSF.append(jets[ijet].btagSF_deepjet_M)
          tightJets_btagSF_up.append(jets[ijet].btagSF_deepjet_M_up)
          tightJets_btagSF_down.append(jets[ijet].btagSF_deepjet_M_down)
        else:
          tightJets_btagSF.append(1.0)
          tightJets_btagSF_up.append(1.0)
          tightJets_btagSF_down.append(1.0)

        # propogate changes to MET
        met += tmpJet.Pt()
        met -= jets[ijet].pt
        met_phi += tmpJet.Phi()
        met_phi -= jets[ijet].phi

    nJets = len(tightJets)
      
    while len(tightJets) < 2: 
      tightJets.append(emptyVector)
      tightJets_ind.append(-99)
      tightJets_btag.append(-99.0)
      tightJets_btagSF.append(1.0)
      tightJets_btagSF_up.append(1.0)
      tightJets_btagSF_down.append(1.0)

    jet1_pt = tightJets[0].Pt()
    jet1_eta = tightJets[0].Eta()
    jet1_phi = tightJets[0].Phi()
    jet1_mass = tightJets[0].M()
    jet1_ind = tightJets_ind[0]
    jet1_btag = tightJets_btag[0]
    jet1_btagSF = tightJets_btagSF[0]
    jet1_btagSF_up = tightJets_btagSF_up[0]
    jet1_btagSF_down = tightJets_btagSF_down[0]

    jet2_pt = tightJets[1].Pt()
    jet2_eta = tightJets[1].Eta()
    jet2_phi = tightJets[1].Phi()
    jet2_mass = tightJets[1].M()
    jet2_ind = tightJets_ind[1]
    jet2_btag = tightJets_btag[1]
    jet1_btagSF = tightJets_btagSF[0]
    jet1_btagSF_up = tightJets_btagSF_up[0]
    jet1_btagSF_down = tightJets_btagSF_down[0]

    nbjets = 0
    if jet1_ind in tightJets_mediumDeepFlavB_ind:
      nbjets += 1
    if jet2_ind in tightJets_mediumDeepFlavB_ind:
      nbjets += 1

    self.out.fillBranch("lead_jet_pt", jet1_pt)
    self.out.fillBranch("lead_jet_eta", jet1_eta)
    self.out.fillBranch("lead_jet_phi", jet1_phi)
    self.out.fillBranch("lead_jet_mass", jet1_mass)
    self.out.fillBranch("lead_jet_ind", jet1_ind)
    self.out.fillBranch("lead_jet_btag", jet1_btag)
    self.out.fillBranch("lead_jet_btagSF", jet1_btagSF)
    self.out.fillBranch("lead_jet_btagSF_up", jet1_btagSF_up)
    self.out.fillBranch("lead_jet_btagSF_down", jet1_btagSF_down)

    self.out.fillBranch("trail_jet_pt", jet2_pt)
    self.out.fillBranch("trail_jet_eta", jet2_eta)
    self.out.fillBranch("trail_jet_phi", jet2_phi)
    self.out.fillBranch("trail_jet_mass", jet2_mass)
    self.out.fillBranch("trail_jet_ind", jet2_ind)
    self.out.fillBranch("trail_jet_btag", jet2_btag)
    self.out.fillBranch("trail_jet_btagSF", jet2_btagSF)
    self.out.fillBranch("trail_jet_btagSF_up", jet2_btagSF_up)
    self.out.fillBranch("trail_jet_btagSF_down", jet2_btagSF_down)

    self.out.fillBranch("ngood_jets", nJets)
    self.out.fillBranch("ngood_bjets", nbjets)

    ################################################
    #####    Composite object kinematics       #####
    ################################################

    # Muon channel

    # Pair-production variables
    m_mumu = (highPtMuons[0] + highPtMuons[1]).M()
    m_mumujj = (highPtMuons[0] + highPtMuons[1] + tightJets[0] + tightJets[1]).M()
    st_mumujj = highPtMuons[0].Pt() + highPtMuons[1].Pt() + tightJets[0].Pt() + tightJets[1].Pt()

    # Single-production variables
    m_mumuj = (highPtMuons[0] + highPtMuons[1] + tightJets[0]).M()
    st_mumuj = highPtMuons[0].Pt() + highPtMuons[1].Pt() + tightJets[0].Pt()

    # DeltaR
    dr_mumu = abs(highPtMuons[0].DeltaR(highPtMuons[1]))
    dr_jj = abs(tightJets[0].DeltaR(tightJets[1]))
    dr_mu1j1 = abs(highPtMuons[0].DeltaR(tightJets[0]))
    dr_mu1j2 = abs(highPtMuons[0].DeltaR(tightJets[1]))
    dr_mu2j1 = abs(highPtMuons[1].DeltaR(tightJets[0]))
    dr_mu2j2 = abs(highPtMuons[1].DeltaR(tightJets[1]))
    dr_mumuj1 = abs(m_mumu.DeltaR(tightJets[1]))

    # DeltaPhi
    dphi_mumu = abs(highPtMuons[0].DeltaPhi(highPtMuons[1]))
    dphi_jj = abs(tightJets[0].DeltaPhi(tightJets[1]))
    dphi_mu1j1 = abs(highPtMuons[0].DeltaPhi(tightJets[0]))
    dphi_mu1j2 = abs(highPtMuons[0].DeltaPhi(tightJets[1]))
    dphi_mu2j1 = abs(highPtMuons[1].DeltaPhi(tightJets[0]))
    dphi_mu2j2 = abs(highPtMuons[1].DeltaPhi(tightJets[1]))
    dphi_mumuj1 = abs(m_mumu.DeltaPhi(tightJets[1]))

    # Electron channel

    # Pair-production variables
    m_ee = (heepElectrons[0] + heepElectrons[1]).M()
    m_eejj = (heepElectrons[0] + heepElectrons[1] + tightJets[0] + tightJets[1]).M()
    st_eejj = heepElectrons[0].Pt() + heepElectrons[1].Pt() + tightJets[0].Pt() + tightJets[1].Pt()

    # Single-production variables
    m_eej = (heepElectrons[0] + heepElectrons[1] + tightJets[0]).M()
    st_eej = heepElectrons[0].Pt() + heepElectrons[1].Pt() + tightJets[0].Pt()

    # DeltaR
    dr_ee = abs(heepElectrons[0].DeltaR(heepElectrons[1]))
    dr_jj = abs(tightJets[0].DeltaR(tightJets[1]))
    dr_e1j1 = abs(heepElectrons[0].DeltaR(tightJets[0]))
    dr_e1j2 = abs(heepElectrons[0].DeltaR(tightJets[1]))
    dr_e2j1 = abs(heepElectrons[1].DeltaR(tightJets[0]))
    dr_e2j2 = abs(heepElectrons[1].DeltaR(tightJets[1]))
    dr_eej1 = abs(m_ee.DeltaR(tightJets[1]))

    # DeltaPhi
    dphi_ee = abs(heepElectrons[0].DeltaPhi(heepElectrons[1]))
    dphi_jj = abs(tightJets[0].DeltaPhi(tightJets[1]))
    dphi_e1j1 = abs(heepElectrons[0].DeltaPhi(tightJets[0]))
    dphi_e1j2 = abs(heepElectrons[0].DeltaPhi(tightJets[1]))
    dphi_e2j1 = abs(heepElectrons[1].DeltaPhi(tightJets[0]))
    dphi_e2j2 = abs(heepElectrons[1].DeltaPhi(tightJets[1]))
    dphi_eej1 = abs(m_ee.DeltaPhi(tightJets[1]))


    # Fill branches with variables

    self.out.fillBranch("M_mumu",m_mumu)
    self.out.fillBranch("M_mumujj",m_mumujj)
    self.out.fillBranch("ST_mumujj",st_mumujj)
    self.out.fillBranch("M_mumuj",m_mumuj)
    self.out.fillBranch("ST_mumuj",st_mumuj)

    self.out.fillBranch("M_ee",m_ee)
    self.out.fillBranch("M_eejj",m_eejj)
    self.out.fillBranch("ST_eejj",st_eejj)
    self.out.fillBranch("M_eej",m_eej)
    self.out.fillBranch("ST_eej",st_eej)

    self.out.fillBranch("DeltaR_mumu",dr_mumu)
    self.out.fillBranch("DeltaR_jj",dr_jj)
    self.out.fillBranch("DeltaR_leadMuon_leadJet",dr_mu1j1)
    self.out.fillBranch("DeltaR_leadMuon_trailJet",dr_mu1j2)
    self.out.fillBranch("DeltaR_trailMuon_leadJet",dr_mu2j1)
    self.out.fillBranch("DeltaR_trailMuon_trailJet",dr_mu2j2)
    self.out.fillBranch("DeltaR_mumu_trailJet",dr_mumuj1)

    self.out.fillBranch("DeltaPhi_mumu",dphi_mumu)
    self.out.fillBranch("DeltaPhi_jj",dphi_jj)
    self.out.fillBranch("DeltaPhi_leadMuon_leadJet",dphi_mu1j1)
    self.out.fillBranch("DeltaPhi_leadMuon_trailJet",dphi_mu1j2)
    self.out.fillBranch("DeltaPhi_trailMuon_leadJet",dphi_mu2j1)
    self.out.fillBranch("DeltaPhi_trailMuon_trailJet",dphi_mu2j2)
    self.out.fillBranch("DeltaPhi_mumu_trailJet",dphi_mumuj1 )

    self.out.fillBranch("DeltaR_ee",dr_ee)
    self.out.fillBranch("DeltaR_leadEle_leadJet",dr_e1j1)
    self.out.fillBranch("DeltaR_leadEle_trailJet",dr_e1j2)
    self.out.fillBranch("DeltaR_trailEle_leadJet",dr_e2j1)
    self.out.fillBranch("DeltaR_trailEle_trailJet",dr_e2j2)
    self.out.fillBranch("DeltaR_ee_trailJet",dr_eej1)

    self.out.fillBranch("DeltaPhi_ee",dphi_ee)
    self.out.fillBranch("DeltaPhi_leadEle_leadJet",dphi_e1j1)
    self.out.fillBranch("DeltaPhi_leadEle_trailJet",dphi_e1j2)
    self.out.fillBranch("DeltaPhi_trailEle_leadJet",dphi_e2j1)
    self.out.fillBranch("DeltaPhi_trailEle_trailJet",dphi_e2j2)
    self.out.fillBranch("DeltaPhi_ee_trailJet",dphi_eej1 )

    ################################################
    #####        Leptoquark candidates         #####
    ################################################

    lqCandidates = {}

    # Muon channel

    # Pair-production variables

    lqCandidates['mumujj'] = {
      'lq1a': max( (highPtMuons[0] + tightJets[0]).M(), (highPtMuons[1] + tightJets[1]).M() ),
      'lq2a': min( (highPtMuons[0] + tightJets[0]).M(), (highPtMuons[1] + tightJets[1]).M() ),
      'lq1b': max( (highPtMuons[0] + tightJets[1]).M(), (highPtMuons[1] + tightJets[0]).M() ),
      'lq2b': min( (highPtMuons[0] + tightJets[1]).M(), (highPtMuons[1] + tightJets[0]).M() )
      }

    lqDiffs_mumujj = {
      'a': abs(lqCandidates['mumujj']['lq1a'] - lqCandidates['mumujj']['lq2a']),
      'b': abs(lqCandidates['mumujj']['lq1b'] - lqCandidates['mumujj']['lq2b']),
    }
    
    minDiff_mumujj = min(lqDiffs_mumujj, key=lqDiffs_mumujj.get)
    m_mumujj1 = lqCandidates['lq1'+minDiff_mumujj]
    m_mumujj2 = lqCandidates['lq2'+minDiff_mumujj]

    # Single-production variables
    
    lqCandidates['mumuj'] = {
      'lqa': max( (highPtMuons[0] + tightJets[0]).M(), (highPtMuons[1] + tightJets[0]).M() ),
      'lqb': min( (highPtMuons[0] + tightJets[0]).M(), (highPtMuons[1] + tightJets[0]).M() )
      }

    m_mumuj1 = lqCandidates['lqa']
    m_mumuj2 = lqCandidates['lqb']

    # Electron channel

    # Pair-production variables

    lqCandidates['eejj'] = {
      'lq1a': max( (heepElectrons[0] + tightJets[0]).M(), (heepElectrons[1] + tightJets[1]).M() ),
      'lq2a': min( (heepElectrons[0] + tightJets[0]).M(), (heepElectrons[1] + tightJets[1]).M() ),
      'lq1b': max( (heepElectrons[0] + tightJets[1]).M(), (heepElectrons[1] + tightJets[0]).M() ),
      'lq2b': min( (heepElectrons[0] + tightJets[1]).M(), (heepElectrons[1] + tightJets[0]).M() )
      }

    lqDiffs_eejj = {
      'a': abs(lqCandidates['eejj']['lq1a'] - lqCandidates['eejj']['lq2a']),
      'b': abs(lqCandidates['eejj']['lq1b'] - lqCandidates['eejj']['lq2b']),
    }
    
    minDiff_eejj = min(lqDiffs_eejj, key=lqDiffs_eejj.get)
    m_eejj1 = lqCandidates['lq1'+minDiff_eejj]
    m_eejj2 = lqCandidates['lq2'+minDiff_eejj]

    # Single-production variables

    lqCandidates['eej'] = {
      'lqa': max( (heepElectrons[0] + tightJets[0]).M(), (heepElectrons[1] + tightJets[0]).M() ),
      'lqb': min( (heepElectrons[0] + tightJets[0]).M(), (heepElectrons[1] + tightJets[0]).M() )
      }

    m_eej1 = lqCandidates['lqa']
    m_eej2 = lqCandidates['lqb']

    # Fill branches with variables
    
    self.out.fillBranch("M_mumujj_lead",m_mumujj1)
    self.out.fillBranch("M_mumujj_trail",m_mumujj2)
    self.out.fillBranch("M_mumuj_lead",m_mumuj1)
    self.out.fillBranch("M_mumuj_trail",m_mumuj2)

    self.out.fillBranch("M_eejj_lead",m_eejj1)
    self.out.fillBranch("M_eejj_trail",m_eejj2)
    self.out.fillBranch("M_eej_lead",m_eej1)
    self.out.fillBranch("M_eej_trail",m_eej2)


    ################################################
    #####            Muon Triggers             #####
    ################################################

    # To do

    ################################################
    #####             MET Filters              #####
    ################################################

    pass_metFilters = True
    pass_metFilters *= bool(event.Flag_goodVertices)
    pass_metFilters *= bool(event.Flag_HBHENoiseFilter)
    pass_metFilters *= bool(event.Flag_HBHENoiseIsoFilter)
    pass_metFilters *= bool(event.Flag_eeBadScFilter)
    pass_metFilters *= bool(event.Flag_EcalDeadCellTriggerPrimitiveFilter)
    pass_metFilters *= bool(event.Flag_globalSuperTightHalo2016Filter)
    pass_metFilters *= bool(event.Flag_BadPFMuonFilter)

    self.out.fillBranch("met_pt",met)
    self.out.fillBranch("met_phi",met_phi)
    self.out.fillBranch("met_filter",pass_metFilters)

    return True

LQ2016apv = lambda: LQProducer("2016apv")
LQ2016 = lambda: LQProducer("2016")
LQ2017 = lambda: LQProducer("2017")
LQ2018 = lambda: LQProducer("2018")
