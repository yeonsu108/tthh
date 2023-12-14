# conda activate py36
import ROOT
import numpy as np
import os

# Criteria
OS = "B_1202_14TeV_Matchable"
Tree = "Delphes"

# Input
indir = "./samples2/" 
tthh   = ROOT.RDataFrame(Tree, indir + OS + "_tthh_di.root")
#ttbbbb = ROOT.RDataFrame(Tree, indir + OS + "_ttbbbb_di.root")
#ttbbcc = ROOT.RDataFrame(Tree, indir + OS + "_ttbbcc_di.root")
#ttbb   = ROOT.RDataFrame(Tree, indir + OS + "_ttbb_di.root")
dfs = {"_tthh": tthh} #"_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb }

# Output
outdir = "./plots/"+ OS + "/"
try : os.makedirs(outdir)
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
            if _xmax < 20: nbin, xmin, xmax = 9, -2, 7 #int(_xmax+2), int(_xmax+2)
            if _xmax < 1: nbin, xmax = 10 ,1
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
 "Correct_Chi"   
]

hists_S1 = hists_S0 + [
]

hists_S2 = hists_S1 + [
]

hists_S3 = hists_S2 + [
]


# Draw #
drawHistoSingle(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton 
for dfname, df in dfs.items():
    df = df.Filter("Correct_Chi==1")
    dfs[dfname] = df
drawHistoSingle(hists_S1, dfs, "_S1")

# ES2 : Jet
#for dfname, df in dfs.items():
#    df = df.Filter("Jet_size >= 5")\
#           .Define("Jet1_pt", "Jet_pt[0]")\
#           .Define("Jet2_pt", "Jet_pt[1]")\
#           .Define("Jet3_pt", "Jet_pt[2]")\
#           .Define("Jet4_pt", "Jet_pt[3]")\
#           .Define("Jet5_pt", "Jet_pt[4]")
#    dfs[dfname] = df
#drawHistoSingle(hists_S2, dfs, "_S2")

# ES3 : bJet
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHistoSingle(hists_S3, dfs, "_S3")
