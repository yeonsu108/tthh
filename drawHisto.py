import ROOT
import numpy as np
ROOT.gStyle.SetOptStat(0)

inFile1 = "/home/stiger97/TTHH/TTHH_Histo/PROCESSED/tthh.root"
inFile2 = "/home/stiger97/TTHH/TTHH_Histo/PROCESSED/ttbbbb.root"
inFile3 = "/home/stiger97/TTHH/TTHH_Histo/PROCESSED/ttbbcc.root"
inFile4 = "/home/stiger97/TTHH/TTHH_Histo/PROCESSED/ttbb.root"
TreeName = "Delphes"

# RDF
tthh = ROOT.RDataFrame(TreeName, inFile1) # tthh
ttbbbb = ROOT.RDataFrame(TreeName, inFile2) # ttbbbb
ttbbcc = ROOT.RDataFrame(TreeName, inFile3) # ttbbcc
ttbb = ROOT.RDataFrame(TreeName, inFile4) # ttbb

##**##  Additional Branches after ES
def Define_Branch(df):
    df = df.Define("Lepton_size", "Muon_size + Electron_size")\
           .Redefine("bJet1_pt", "bJet_pt[0]")\
           .Redefine("bJet2_pt", "bJet_pt[1]")\
           .Redefine("bJet3_pt", "bJet_pt[2]")\
           .Redefine("bJet4_pt", "bJet_pt[3]")\
           .Redefine("Lepton1_pt", "Lepton_pt[0]")\
           .Redefine("Lepton2_pt", "Lepton_pt[1]")
    return df

# ES
tthh_1S1 = tthh.Filter("Electron_size + Muon_size == 1")
tthh_1S1 = Define_Branch(tthh_1S1)
tthh_1S2 = tthh.Filter("Electron_size + Muon_size == 1 && Jet_size>=4")
tthh_1S2 = Define_Branch(tthh_1S2)
tthh_1S3 = tthh.Filter("Electron_size + Muon_size == 1 && Jet_size>=4 && bJet_size>=4")
tthh_1S3 = Define_Branch(tthh_1S3)
tthh_2S1 = tthh.Filter("Electron_size + Muon_size >= 2")
tthh_2S1 = Define_Branch(tthh_2S1)
tthh_2S2 = tthh.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4")
tthh_2S2 = Define_Branch(tthh_2S2)
tthh_2S3 = tthh.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4 && bJet_size>=4")
tthh_2S3 = Define_Branch(tthh_2S3)

ttbbbb_1S1 = ttbbbb.Filter("Electron_size + Muon_size == 1")
ttbbbb_1S1 = Define_Branch(ttbbbb_1S1)
ttbbbb_1S2 = ttbbbb.Filter("Electron_size + Muon_size == 1 && Jet_size>=4")
ttbbbb_1S2 = Define_Branch(ttbbbb_1S2)
ttbbbb_1S3 = ttbbbb.Filter("Electron_size + Muon_size == 1 && Jet_size>=4 && bJet_size>=4")
ttbbbb_1S3 = Define_Branch(ttbbbb_1S3)
ttbbbb_2S1 = ttbbbb.Filter("Electron_size + Muon_size >= 2")
ttbbbb_2S1 = Define_Branch(ttbbbb_2S1)
ttbbbb_2S2 = ttbbbb.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4")
ttbbbb_2S2 = Define_Branch(ttbbbb_2S2)
ttbbbb_2S3 = ttbbbb.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4 && bJet_size>=4")
ttbbbb_2S3 = Define_Branch(ttbbbb_2S3)

ttbbcc_1S1 = ttbbcc.Filter("Electron_size + Muon_size == 1")
ttbbcc_1S1 = Define_Branch(ttbbcc_1S1)
ttbbcc_1S2 = ttbbcc.Filter("Electron_size + Muon_size == 1 && Jet_size>=4")
ttbbcc_1S2 = Define_Branch(ttbbcc_1S2)
ttbbcc_1S3 = ttbbcc.Filter("Electron_size + Muon_size == 1 && Jet_size>=4 && bJet_size>=4")
ttbbcc_1S3 = Define_Branch(ttbbcc_1S3)
ttbbcc_2S1 = ttbbcc.Filter("Electron_size + Muon_size >= 2")
ttbbcc_2S1 = Define_Branch(ttbbcc_2S1)
ttbbcc_2S2 = ttbbcc.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4")
ttbbcc_2S2 = Define_Branch(ttbbcc_2S2)
ttbbcc_2S3 = ttbbcc.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4 && bJet_size>=4")
ttbbcc_2S3 = Define_Branch(ttbbcc_2S3)

ttbb_1S1 = ttbb.Filter("Electron_size + Muon_size == 1")
ttbb_1S1 = Define_Branch(ttbb_1S1)
ttbb_1S2 = ttbb.Filter("Electron_size + Muon_size == 1 && Jet_size>=4")
ttbb_1S2 = Define_Branch(ttbb_1S2)
ttbb_1S3 = ttbb.Filter("Electron_size + Muon_size == 1 && Jet_size>=4 && bJet_size>=4")
ttbb_1S3 = Define_Branch(ttbb_1S3)
ttbb_2S1 = ttbb.Filter("Electron_size + Muon_size >= 2")
ttbb_2S1 = Define_Branch(ttbb_2S1)
ttbb_2S2 = ttbb.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4")
ttbb_2S2 = Define_Branch(ttbb_2S2)
ttbb_2S3 = ttbb.Filter("Electron_size + Muon_size >= 2 && Jet_size>=4 && bJet_size>=4")
ttbb_2S3 = Define_Branch(ttbb_2S3)

##**##
dfs = {
#"tthh_noCut" : [tthh, 1], "ttbbbb_noCut" : [ttbbbb, 38], "ttbbcc_noCut" : [ttbbcc, 30], "ttbb_noCut" : [ttbb, 46],
#"tthh_1S1" : [tthh_1S1, 1], "ttbbbb_1S1" : [ttbbbb_1S1, 38], "ttbbcc_1S1" : [ttbbcc_1S1, 30], "ttbb_1S1" : [ttbb_1S1, 46],
#"tthh_1S2" : [tthh_1S2, 1], "ttbbbb_1S2" : [ttbbbb_1S2, 38], "ttbbcc_1S2" : [ttbbcc_1S2, 30], "ttbb_1S2" : [ttbb_1S2, 46],
#"tthh_1S3" : [tthh_1S3, 1], "ttbbbb_1S3" : [ttbbbb_1S3, 38], "ttbbcc_1S3" : [ttbbcc_1S3, 30], "ttbb_1S3" : [ttbb_1S3, 46],
#"tthh_2S1" : [tthh_2S1, 1], "ttbbbb_2S1" : [ttbbbb_2S1, 38], "ttbbcc_2S1" : [ttbbcc_2S1, 30], "ttbb_2S1" : [ttbb_2S1, 46],
#"tthh_2S2" : [tthh_2S2, 1], "ttbbbb_2S2" : [ttbbbb_2S2, 38], "ttbbcc_2S2" : [ttbbcc_2S2, 30], "ttbb_2S2" : [ttbb_2S2, 46],
"tthh_2S3" : [tthh_2S3, 1], "ttbbbb_2S3" : [ttbbbb_2S3, 38], "ttbbcc_2S3" : [ttbbcc_2S3, 30], "ttbb_2S3" : [ttbb_2S3, 46],
      }


##**##    Histogram Features  
hists = {
   #     "nGenQuark" : ["nGenQuark", "nGenQuark", 80, 0, 80, "nGenAddQuark", "N"],    
   #     "nGenbQuark" : ["nGenbQuark", "nGenbQuark", 10, 0, 10, "nGenAddbQuark", "N"],   
   #     "nGencQuark" : ["nGencQuark", "nGencQuark", 10, 0, 10, "nGenAddcQuark", "N"],    
   #     "nGenJet" : ["nGenJet", "nGenJet", 12, 0, 12, "nGenAddJet", "N"],
   #     "nGenbJet" : ["nGenbJet", "nGenbJet", 10, 0, 10, "nGenAddbJet", "N"],
   #     "nGencJet" : ["nGencJet", "nGencJet", 10, 0, 10, "nGenAddcJet", "N"],
   #     "nGenlfJet" : ["nGenlfJet", "nGenlfJet", 10, 0, 10, "nGenAddlfJet", "N"],
   #     "nLepformT" : ["nLepFromW", "nLepFromW", 12, 0, 12, "nLepFromW", "N"],
    
   #     "nJet" : ["Jet_size", "Jet_multiplicity", 15, 0, 15, "Jet_size", "N"],
   #     "nbJet" : ["bJet_size", "bJet_multiplicity", 15, 0, 15, "bJet_size", "N"],
   #     "nMuon" : ["Muon_size", "Muon_multiplicity", 15, 0, 15, "Muon_size", "N"],
   #     "nElectron" : ["Electron_size", "Electron_multiplicity", 15, 0, 15, "Electron_size", "N"],
   #     "nLepton" : ["Lepton_size", "Lepton_multiplicity", 15, 0, 15, "Lepton_size", "N"]

   #     "Jet_pt1" : ["Jet1_PT", "Jet_first_leading_order_pt", 40, 0, 400, "Jet1_PT", "Electron p_{T} (GeV)"],
   #     "h_jet_eta1" : ["Jet1_Eta", "Jet_first_leading_order_eta", 8,-4, 4, "Jet1_Eta", "#eta"],
         "bjet1_pt" : ["bJet1_pt", "bJet1_pt", 40, 0, 400, "bJet1_pt", "bJet p_{T} (GeV)"],
         "bjet2_pt" : ["bJet2_pt", "bJet2_pt", 40, 0, 400, "bJet2_pt", "bJet p_{T} (GeV)"],
         "bjet3_pt" : ["bJet3_pt", "bJet3_pt", 40, 0, 400, "bJet3_pt", "bJet p_{T} (GeV)"],
         "bjet4_pt" : ["bJet4_pt", "bJet4_pt", 40, 0, 400, "bJet4_pt", "bJet p_{T} (GeV)"]
   #     "h_bjet_eta1" : ["bJet1_Eta", "bJet_first_leading_order_eta", 16, -4, 4, "bJet1_Eta", "#eta"],
   #     "h_muon_pt1" : ["Muon1_PT", "Muon_first_pt", 40, 0, 400, "Muon1_PT", "Muon p_{T} (GeV)"]
   #     "h_muon_pt2" : ["Muon2_PT", "Muon_second_pt", 40, 0, 400, "Muon2_PT", "Muon p_{T} (GeV)"],
   #     "h_muon_eta1" : ["Muon1_Eta", "Muon_first_leading_order_eta", 16, -4, 4, "Muon1_Eta", "eta"],
   #     "h_electron_pt1" : ["Electron1_PT", "Electron_first_leading_order_pt", 40, 0, 400, "Electron1_PT", "Electron p_{T} (GeV)"]
   #     "h_electron_eta1" : ["Electron1_Eta", "Electron_first_leading_order_eta", 16, -4, 4, "Electron1_Eta", "#eta"]
   #     "h_Lepton_pt1" : ["Lep1_PT", "Lepton_first_pt", 40, 0, 400, "Lep1_PT", "Lep p_{T} (GeV)"]
         }

######################################################################################################

canvas = ROOT.TCanvas("c", "c", 400, 400)

def DrawPrintHisto(hists, dfs):
    for hist_name, hist in hists.items():
        hist_dict = {}
        legend = ROOT.TLegend(0.55, 0.9, 0.75, 0.8)
    
        for df_name, df in dfs.items() :
            h = df[0].Histo1D(ROOT.RDF.TH1DModel(hist[0], hist[0], hist[2], hist[3], hist[4]), hist[5])
            h.SetMaximum(h.GetMaximum()* 2.0)  ##**##
            h.GetXaxis().SetTitle(hist[6])
            h.GetYaxis().SetTitle("Entries")
            h.SetLineColor(df[1])
            h.SetLineWidth(2)
            legend.AddEntry(h.GetValue(), df_name, "l")
            hist_dict[hist_name + "_" + df_name] = h
        
        first = True
        for histName , h in hist_dict.items():  
            if first :
                h.DrawNormalized("hist")
                first = False
            else :
                h.DrawNormalized("same")
        legend.Draw()     
        canvas.Print(histName + ".pdf") ##**##
        canvas.Clear()

DrawPrintHisto(hists, dfs)
