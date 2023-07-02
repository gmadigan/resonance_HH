from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import ROOT
import os
import random
ROOT.PyConfig.IgnoreCommandLineOptions = True


def mk_safe(fct, *args):
    try:
        return fct(*args)
    except Exception as e:
        if any('Error in function boost::math::erf_inv' in arg
               for arg in e.args):
            print(
                'WARNING: catching exception and returning -1. Exception arguments: %s'
                % e.args)
            return -1.
        else:
            raise e

def MERParametrization(p, eta, year):
    #
    # Definition of muon momentum resolution parameterization function (sigma)
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016#Momentum_Resolution
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017#Momentum_Resolution
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018#Momentum_Resolution
    #
    # Note: p is muon 3-momentum--not transverse momentum

    coeffsAllEta = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    coeffs = coeffsAllEta[0]

    if year == '2016a': coeffsAllEta = [[0.00112, 6.87e-5, -3.88e-8, 9.03e-12],[0.013, 6.93e-5, -3.46e-8, 7.72e-12]]
    elif year == '2016b': coeffsAllEta = [[0.0102, 6.77e-5, -3.72e-8, 8.53e-12],[0.0129, 6.48e-5, -3.04e-8, 6.63e-12]]
    elif year == '2017': coeffsAllEta = [[0.0104, 6.11e-5, -3.31e-8, 6.73e-12],[0.0121, 5.92e-05, -2.61e-8, 5.11e-12]]
    elif year == '2018': coeffsAllEta = [[0.0108, 5.93e-5,  -3.08e-8, 6.04e-12],[0.0136, 5.47e-5, -2.3e-8, 4.66e-12]]

    if abs(eta) < 1.2 : coeffs = coeffsAllEta[0]
    elif 1.2 <= abs(eta) < 2.4 : coeffs = coeffsAllEta[1]

    sigmaToReturn = coeffs[0] + coeffs[1]*p + coeffs[2]*p*p + coeffs[3]*p*p*p

    return sigmaToReturn

def SmearMuonCollections(ptCollection, etaCollection, phiCollection, year, isSystematic):
    #
    # Following perscription for extra muon momentum resolution smearing of 5 to 10% from resolution uncertainty and 10% for systematic uncertainty
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2016#Momentum_Resolution
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2017#Momentum_Resolution
    # https://twiki.cern.ch/twiki/bin/view/CMS/MuonUL2018#Momentum_Resolution

    smearedPtCollection = ptCollection
    smearedEtaCollection = etaCollection
    smearedPhiCollection = phiCollection

    smearConst = [0.0, 0.0]
    if year == '2016a': smearConst = [0.0, 0.46]
    elif year == '2016b': smearConst = [0.0, 0.46]
    elif year == '2017': smearConst = [0.3202, 0.46]
    elif year == '2018': smearConst = [0.46, 0.46]

    # loop through muons
    Plist = []
    smearedPlist = []
    smearingList = []

    if not isSystematic and (year == '2016a' or year == '2016b' or abs(etaCollection[n]) < 1.2): pass
    else:
        for n in range(len(ptCollection)):
            smearedLorentz = TLorentzVector()
            origLorentz = TLorentzVector()
            origLorentz.SetPtEtaPhiM(ptCollection[n], etaCollection[n], phiCollection[n], 0)

            # Smearing is performed on P, convert from Pt, Eta, Phi to cartesian 3-momentum and back
            origPx = origLorentz.Px()
            origPy = origLorentz.Py()
            origPz = origLorentz.Pz()
            origP = origLorentz.P()
            # Smear momenta here
            if isSystematic: smearing = 1 + random.gauss(0.0, MERParametrization(origP, etaCollection[n])*smearConst[1], year)
            else: smearing = 1 + random.gauss(0.0, MERParametrization(origP, etaCollection[n])*smearConst[0], year)
            smearedPx = origPx*smearing
            smearedPy = origPy*smearing
            smearedPz = origPz*smearing
            smearedE = math.sqrt(smearedPx*smearedPx + smearedPy*smearedPy + smearedPz*smearedPz)
            smearedLorentz.SetPxPyPzE(smearedPx, smearedPy, smearedPz, smearedE)

            # Return the smeared Pt, Eta, Phi collections
            smearedPtCollection[n] = smearedLorentz.Pt()
            smearedEtaCollection[n] = smearedLorentz.Eta()
            smearedPhiCollection[n] = smearedLorentz.Phi()

    return [smearedPtCollection, smearedEtaCollection, smearedPhiCollection]

class highPtMuonScaleResProducer(Module):
    def __init__(self, rc_dir, rc_corrections, dataYear):
        p_postproc = '%s/src/PhysicsTools/NanoAODTools/python/postprocessing' % os.environ[
            'CMSSW_BASE']
        p_roccor = p_postproc + '/analysis/data/' + rc_dir
        if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries():
            p_helper = '%s/RoccoR.cc' % p_roccor
            print('Loading C++ helper from ' + p_helper)
            ROOT.gROOT.ProcessLine('.L ' + p_helper)
        self._roccor = ROOT.RoccoR(p_roccor + '/' + rc_corrections)
        self.dataYear =  dataYear

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Muon_highPtScaleRoccor_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtScaleRoccorUp_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtScaleRoccorDown_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmeared_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedUp_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedDown_pt", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmeared_eta", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedUp_eta", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedDown_eta", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmeared_phi", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedUp_phi", "F", lenVar="nMuon")
        self.out.branch("Muon_highPtResSmearedDown_phi", "F", lenVar="nMuon")
        self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        muons = Collection(event, "Muon")
        if self.is_mc:
            genparticles = Collection(event, "GenPart")
        roccor = self._roccor
        if self.is_mc:
            pt_corr = []
            pt_err = []
            eta = []
            phi = []
            for mu in muons:
                tunepPt = mu.pt * mu.tunepRelPt
                genIdx = mu.genPartIdx
                if genIdx >= 0 and genIdx < len(genparticles):
                    genMu = genparticles[genIdx]
                    if tunepPt < 200.0:
                        pt_corr.append(tunepPt *
                                       mk_safe(roccor.kSpreadMC, mu.charge, tunepPt,
                                               mu.eta, mu.phi, genMu.pt))
                        pt_err.append(tunepPt *
                                      mk_safe(roccor.kSpreadMCerror, mu.charge,
                                              tunepPt, mu.eta, mu.phi, genMu.pt))
                    else:
                        pt_corr.append(tunepPt)
                        pt_err.append(0.0)
                else:
                    if tunepPt < 200.0:
                        u1 = random.uniform(0.0, 1.0)
                        pt_corr.append(
                            tunepPt * mk_safe(roccor.kSmearMC, mu.charge, tunepPt,
                                            mu.eta, mu.phi, mu.nTrackerLayers, u1))
                        pt_err.append(
                            tunepPt * mk_safe(roccor.kSmearMCerror, mu.charge, tunepPt,
                                            mu.eta, mu.phi, mu.nTrackerLayers, u1))
                    else:
                        pt_corr.append(muPt)
                        pt_err.append(0.0)
                eta.append(mu.eta)
                phi.append(mu.phi)

        else:
            pt_corr = list( mu.pt * mu.tunepRelPt * mk_safe(roccor.kScaleDT, mu.charge, mu.pt * mu.tunepRelPt, mu.eta, mu.phi) if mu.pt* mu.tunepRelPt < 200.0 else mu.pt * mu.tunepRelPt for mu in muons)
            pt_err = list( mu.pt * mu.tunepRelPt * mk_safe(roccor.kScaleDTerror, mu.charge, mu.pt * mu.tunepRelPt, mu.eta, mu.phi) if mu.pt * mu.tunepRelPt < 200.0 else 0.0 for mu in muons)
            eta = list(mu.eta for mu in muons)
            phi = list(mu.phi for mu in muons)
        
        pt_corr_up = list(
            max(pt_corr[imu] + pt_err[imu], 0.0)
            for imu, mu in enumerate(muons))
        pt_corr_down = list(
            max(pt_corr[imu] - pt_err[imu], 0.0)
            for imu, mu in enumerate(muons))

        if self.isMC:
            [pt_smear, eta_smear, phi_smear] = SmearMuonCollections(pt_corr, eta, phi, self.dataYear, False)
            [pt_smear_up, eta_smear_up, phi_smear_up] = SmearMuonCollections(pt_corr_up, eta, phi, self.dataYear, False)
            [pt_smear_down, eta_smear_down, phi_smear_down] = SmearMuonCollections(pt_corr_down, eta, phi, self.dataYear, False)
        else:
            [pt_smear, eta_smear, phi_smear] = [pt_corr, eta, phi]
            [pt_smear_up, eta_smear_up, phi_smear_up] = [pt_corr_up, eta, phi]
            [pt_smear_down, eta_smear_down, phi_smear_down] = [pt_corr_down, eta, phi]

        self.out.fillBranch("Muon_highPtScaleRoccor_pt", pt_corr)
        self.out.fillBranch("Muon_highPtScaleRoccorUp_pt", pt_corr_up)
        self.out.fillBranch("Muon_highPtScaleRoccorDown_pt", pt_corr_down)
        self.out.fillBranch("Muon_highPtResSmeared_pt", pt_smear)
        self.out.fillBranch("Muon_highPtResSmearedUp_pt", pt_smear_up)
        self.out.fillBranch("Muon_highPtResSmearedDown_pt", pt_smear_down)
        self.out.fillBranch("Muon_highPtResSmeared_eta", eta_smear)
        self.out.fillBranch("Muon_highPtResSmearedUp_eta", eta_smear_up)
        self.out.fillBranch("Muon_highPtResSmearedDown_eta", eta_smear_down)
        self.out.fillBranch("Muon_highPtResSmeared_phi", phi_smear)
        self.out.fillBranch("Muon_highPtResSmearedUp_phi", phi_smear_up)
        self.out.fillBranch("Muon_highPtResSmearedDown_phi", phi_smear_down)

        return True

highPtMuonScaleRes2016a = lambda: highPtMuonScaleResProducer('roccor.Run2UL.v5','RoccoR2016aUL.txt', '2016a')
highPtMuonScaleRes2016b = lambda: highPtMuonScaleResProducer('roccor.Run2UL.v5','RoccoR2016bUL.txt', '2016')
highPtMuonScaleRes2017 = lambda: highPtMuonScaleResProducer('roccor.Run2UL.v5','RoccoR2017UL.txt', '2017')
highPtMuonScaleRes2018 = lambda: highPtMuonScaleResProducer('roccor.Run2UL.v5','RoccoR2018UL.txt', '2018')