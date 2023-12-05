import ROOT
import numpy as np
import os
ROOT.gStyle.SetOptStat(0)

# RDF
PRE = "FULL_1202_14TeV"
Tree = "Delphes"

# Input
indir = "./samples2/" 
tthh   = ROOT.RDataFrame(Tree, indir + PRE + "_tthh_di.root")
tthbb  = ROOT.RDataFrame(Tree, indir + PRE + "_tthbb_di.root")
ttbbbb = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbbb_di.root")
ttbbcc = ROOT.RDataFrame(Tree, indir + PRE + "_ttbbcc_di.root")
ttbb   = ROOT.RDataFrame(Tree, indir + PRE + "_ttbb_di.root")
dfs = {"_tthh": tthh, "_ttbbbb": ttbbbb, "_ttbbcc": ttbbcc, "_ttbb": ttbb, "_tthbb": tthbb}

# OutDirectory. Modify!!
outdir = "./plots/" + PRE + "/Same/"
try : os.makedirs(outdir)
except : pass

# Definition
def drawHisto(hists, dfs, flag="_S0"):
    canvas = ROOT.TCanvas("c", "c", 400, 400)
    for hist_name in hists:
        hist_dict = {}
        legend = ROOT.TLegend(0.60, 0.9, 0.9, 0.70)
        legend.SetX1(0.65)
        legend.SetX2(0.9)
        legend.SetY1(0.75)
        legend.SetY2(0.9)
        hist_title = hist_name.replace("_size", " multiplicity")
        hist_title = hist_title.replace("_", " ")

        ymax, color = 0, 1
        for df_name, df in dfs.items():
            nbin, xmin, xmax = 60, 0, 300
            _xmax = df.Max(hist_name).GetValue()
            if _xmax < 100: xmax = 100
            if _xmax < 20: nbin, xmin, xmax = 12, -6, 6 #int(_xmax+2), int(_xmax+2)
            if _xmax < 0: nbin, xmin, xmax = 20, -4, 4
            h = df.Histo1D(ROOT.RDF.TH1DModel(hist_name, hist_title, nbin, xmin, xmax), hist_name)
            if ymax < h.GetMaximum(): ymax = h.GetMaximum()
            h.GetXaxis().SetTitle(hist_title)
            h.GetYaxis().SetTitle("Normalized Entries")
            h.GetYaxis().SetTitleOffset(1.6)           
            h.SetLineColor(color)
            color+=1
            if color in [5] :
                color += 1
            h.SetLineWidth(2)
            legend.AddEntry(h.GetValue(), df_name, "l")
            hist_dict[hist_name + "_" + df_name] = h

        first = True
        for _tmp, h in hist_dict.items():
            h.SetMaximum(ymax * 1.2)
            if first:
                h.DrawNormalized("histE")
                first = False
            else: h.DrawNormalized("sameE")
        legend.Draw()
        canvas.Print(outdir + PRE + "_" + hist_name + flag + ".pdf")
        canvas.Clear()


# Histogram Config #
hists_S0 = [
#    "close_Higgs_mass", "Chi_min"

    "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_size",
     "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_size",
     "Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass",
     "Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass",
     "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass",
     "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass",
     "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass",
     "j1j2_dr", "j1j3_dr", "j1j4_dr", "j1j5_dr", "j2j3_dr", "j2j4_dr", "j2j5_dr",
     "j3j4_dr", "j3j5_dr", "j4j5_dr",

     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass",
     "bJet5_pt", "bJet5_eta", "bJet5_phi", "bJet5_mass",
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b1b5_dr", "b2b3_dr", "b2b4_dr", "b2b5_dr", "b3b4_dr", "b3b5_dr", "b4b5_dr",

     "Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "nMuon", "Lep_size",
     "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "nElectron",
     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",

     "MET_E", "MET_Eta", "MET_Phi",

     "bb_dr", "b_Vars", "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
    "jj_dr", "j_Vars", "jj_avg_dr", "jj_max_dr", "jj_min_dr", "j_ht", "jj_dEta_WhenMaxdR", "j_cent", "jj_max_deta", "jj_max_mass", "jj_twist",
    "close_Higgs_pt", "close_Higgs_eta", "close_Higgs_phi", "close_Higgs_mass",

     "isMatchable", "Matched_idx1", "Matched_idx2", "Correct_Chi", "Chi_min",
]

hists_S1 = hists_S0 + [
]

hists_S2 = hists_S1 + [
]

hists_S3 = hists_S2 + [
#    "bJet1_pt", "bJet2_pt", "bJet3_pt"
]


# Draw #
drawHisto(hists_S0, dfs, "_S0")

# Event Selection and Draw
# ES1 : Lepton
for dfname, df in dfs.items():
    df = df.Filter("Lep_size >= 2")
#           .Define("GenbJet1_pt", "GenbJet_pt[0]")\
#           .Define("GenbJet2_pt", "GenbJet_pt[1]")\
#           .Define("GenbJet3_pt", "GenbJet_pt[2]")\
#           .Define("GenbJet4_pt", "GenbJet_pt[3]")
    dfs[dfname] = df
#drawHisto(hists_S1, dfs, "_S1")

# ES2 : Jet
#for dfname, df in dfs.items():
#    df = df.Filter("Jet_size >= 5")\
#           .Define("Jet1_pt", "Jet_pt[0]")\
#           .Define("Jet2_pt", "Jet_pt[1]")\
#           .Define("Jet3_pt", "Jet_pt[2]")\
#           .Define("Jet4_pt", "Jet_pt[3]")\
#           .Define("Jet5_pt", "Jet_pt[4]")
#    dfs[dfname] = df
#drawHisto(hists_S2, dfs, "_S2")

# ES3 : bJet
#for dfname, df in dfs.items():
#    df = df.Filter("bJet_size >= 3")\
#           .Define("bJet1_pt", "bJet_pt[0]")\
#           .Define("bJet2_pt", "bJet_pt[1]")\
#           .Define("bJet3_pt", "bJet_pt[2]")
#    dfs[dfname] = df
#drawHisto(hists_S3, dfs, "_S3")

