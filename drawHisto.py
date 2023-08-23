# conda activate py36
import ROOT
import numpy as np
import os

# Criteria
DATA = "_Semi"
OS = "OS" + DATA
Tree = "Delphes"

# Input
indir = "./samples1/" 
tthh   = ROOT.RDataFrame(Tree, indir + OS + "_tthh.root")
ttbbbb = ROOT.RDataFrame(Tree, indir + OS + "_ttbbbb.root")
ttbbcc = ROOT.RDataFrame(Tree, indir + OS + "_ttbbcc.root")
ttbb   = ROOT.RDataFrame(Tree, indir + OS + "_ttbb.root")
dfs = {"_tthh": tthh, "_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb }

# Output
outdir = "./plots/"+ OS + "/"
try : os.mkdir(outdir)
except : pass

# Definition
def drawHistoSingle(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}

        for df_name, df in dfs.items():
            nbin, xmin, xmax = 40, 0, 400
            _xmax = df.Max(hist_name).GetValue()
            if _xmax < 100: xmax = 100
            if _xmax < 20: nbin, xmax = 12, 6 #int(_xmax+2), int(_xmax+2)
            if _xmax < 0: nbin, xmin, xmax = 20, -4, 4
            hist_title = hist_name.replace("_size", " multiplicity")
            hist_title = hist_title.replace("_", " ")
            hist_title = hist_title + df_name
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            h.GetXaxis().SetTitle(hist_title)
            h.GetYaxis().SetTitle("Normalized Entries")
            h.GetYaxis().SetTitleOffset(1.5)                           
            h.SetLineWidth(2)
            hist_dict[hist_name + df_name] = h

        for _tmp, h in hist_dict.items():
            h.DrawNormalized("hist")
            canvas.Print(outdir + OS + "_" + _tmp + flag + ".pdf")
            canvas.Clear()

# Histogram Config #
hists_S0 = [

    # nGen Quark 
    "nTop", "nHiggs", "nHiggs1", "nHiggs2", "nW", "nGluon", "nGenbQuark", "nGencQuark",
    "nWFromTop", "nbFromTop", "nbFromHiggs", "nbFromHiggs1", "nbFromHiggs2",
    "nbFromGluon", "nHiggsFromTop", "nElectronFromW", "nMuonFromW", "nLepFromW",
    "nGenAddQuark", "nGenAddbQuark", "nGenAddcQuark",

    # nGen Jet 
    "nGenbJet", "nGencJet",
    "nGenbJetFromTop", "nGenbJetFromHiggs", "nGenbJetFromHiggs1", "nGenbJetFromHiggs2",
    "nGenAddJet", "nGenAddbJet", "nGenAddcJet", "nGenAddlfJet",
    "category", "ttbbbb_val", "ttbbcc_val",

    # Gen Quark 4-vector 
    "GenTop_mass", "GenHiggs_mass",
    "GenbQuark_pt", "GenbQuark_eta", "GenbQuark_phi",
    "GenElectronFromW_pt", "GenElectronFromW_eta", "GenElectronFromW_phi",
    "GenMuonFromW_pt", "GenMuonFromW_eta", "GenMuonFromW_phi",
    "GenbQuarkFromTop_pt", "GenbQuarkFromTop_eta", "GenbQuarkFromTop_phi",
    "GenbQuarkFromHiggs_pt", "GenbQuarkFromHiggs_eta", "GenbQuarkFromHiggs_phi",
    "GenAddbQuark_pt", "GenAddbQuark_eta", "GenAddbQuark_phi",

    # Gen Jet 4-vector
    "GenJet_pt", "GenbJet_pt", "GenbJet_eta", "GenbJet_phi",
    "GenbJetFromTop_pt", "GenbJetFromTop_eta", "GenbJetFromTop_phi", "GenbJetFromTop_mass",
    "GenbJetFromHiggs_pt", "GenbJetFromHiggs_eta", "GenbJetFromHiggs_phi", "GenbJetFromHiggs_mass",
    "GenbJetFromHiggs1_pt", "GenbJetFromHiggs1_eta", "GenbJetFromHiggs1_phi", "GenbJetFromHiggs1_mass",
    "GenbJetFromHiggs2_pt", "GenbJetFromHiggs2_eta", "GenbJetFromHiggs2_phi", "GenbJetFromHiggs2_mass",
    "GenAddJet_pt", "GenAddJet_eta", "GenAddJet_phi", "GenAddbJet_pt", "GenAddbJet_eta", "GenAddbJet_phi",

    # Gen dR 
    "GenbJetFromTop_dr", "GenbJetFromHiggs1_dr", "GenbJetFromHiggs2_dr", "GenbJetFromHiggs_dr",
    "HiggsFromWhere",

    # Reco

    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_btag", "Jet_size",
    "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_size",
    "Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "Muon_e", "Muon_size",
    "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "Electron_e", "Electron_size",
    "Lepton_size", "Lepton_pt",
    "Higgs_vars", "Higgs_pt", "Higgs_eta", "Higgs_phi", "Higgs_mass"                      
                      
]

hists_S1 = hists_S0 + [
    "Muon1_pt"
]

hists_S2 = hists_S1 + [
    "Jet1_pt", "Jet2_pt", "Jet3_pt", "Jet4_pt", "Jet5_pt"
]

hists_S3 = hists_S2 + [
    "bJet1_pt", "bJet2_pt", "bJet3_pt" 
]


# Draw #
drawHistoSingle(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton 
for dfname, df in dfs.items():
    df = df.Filter("Muon_size == 1")\
           .Define("Muon1_pt", "Muon_pt[0]") 
    dfs[dfname] = df
drawHistoSingle(hists_S1, dfs, "_S1")

# ES2 : Jet
for dfname, df in dfs.items():
    df = df.Filter("Jet_size >= 5")\
           .Define("Jet1_pt", "Jet_pt[0]")\
           .Define("Jet2_pt", "Jet_pt[1]")\
           .Define("Jet3_pt", "Jet_pt[2]")\
           .Define("Jet4_pt", "Jet_pt[3]")\
           .Define("Jet5_pt", "Jet_pt[4]")
    dfs[dfname] = df
drawHistoSingle(hists_S2, dfs, "_S2")

# ES3 : bJet
for dfname, df in dfs.items():
    df = df.Filter("bJet_size >= 3")\
           .Define("bJet1_pt", "bJet_pt[0]")\
           .Define("bJet2_pt", "bJet_pt[1]")\
           .Define("bJet3_pt", "bJet_pt[2]")
    dfs[dfname] = df
drawHistoSingle(hists_S3, dfs, "_S3")

