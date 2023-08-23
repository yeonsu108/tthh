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
outdir = "./samples2/"
try : os.mkdir(outdir)
except : pass

# Extand variables for DNN
for dfname, df in dfs.items():
    df = df.Filter("Muon_size==1")\
           .Filter("Jet_size>=5")\
           .Filter("bJet_size>=3")\
           .Define("ngoodJets", "Jet_size")\
           .Define("ngoodbJets", "bJet_size")\
           .Define("ngoodElectrons", "Electron_size")\
           .Define("ngoodMuons", "Muon_size")\
           .Define("Jet1_pt", "Jet_pt[0]").Define("Jet1_eta", "Jet_eta[0]").Define("Jet1_mass", "Jet_mass[0]")\
           .Define("Jet2_pt", "Jet_pt[1]").Define("Jet2_eta", "Jet_eta[1]").Define("Jet2_mass", "Jet_mass[1]")\
           .Define("Jet3_pt", "Jet_pt[2]").Define("Jet3_eta", "Jet_eta[2]").Define("Jet3_mass", "Jet_mass[2]")\
           .Define("Jet4_pt", "Jet_pt[3]").Define("Jet4_eta", "Jet_eta[3]").Define("Jet4_mass", "Jet_mass[3]")\
           .Define("Jet5_pt", "Jet_pt[4]").Define("Jet5_eta", "Jet_eta[4]").Define("Jet5_mass", "Jet_mass[4]")\
           .Define("bJet1_pt", "bJet_pt[0]").Define("bJet1_eta", "bJet_eta[0]").Define("bJet1_mass", "bJet_mass[0]")\
           .Define("bJet2_pt", "bJet_pt[1]").Define("bJet2_eta", "bJet_eta[1]").Define("bJet2_mass", "bJet_mass[1]")\
           .Define("bJet3_pt", "bJet_pt[2]").Define("bJet3_eta", "bJet_eta[2]").Define("bJet3_mass", "bJet_mass[2]")\
           .Define("Muon1_pt", "Muon_pt[0]").Define("Muon1_eta", "Muon_eta[0]").Define("Muon1_e", "Muon_e[0]")\
           .Snapshot("Delphes", outdir + OS + dfname + ".root")

