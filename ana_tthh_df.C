
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>
#include <algorithm>

#include "TTree.h"
#include "TFile.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TMath.h>

#include "utility2.h"

// MODIFY!!
void ana_tthh_df(std::string channel, std::string outdir="./samples1/"){
    gSystem->Load("libDelphes");

//    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/test/*.root";
//    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/tag_1*.root";
    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/*.root";
    std::cout << infile << std::endl;
    std::cout << outdir << std::endl;
    auto treename = "Delphes";

    auto _df = ROOT::RDataFrame(treename, infile);
    
    // Constants, Basic Branches //
    auto df0 = _df.Define("T_pid", "int(6)")
                  .Define("H_pid", "int(25)")
                  .Define("g_pid", "int(21)")
                  .Define("W_pid", "int(24)")
                  .Define("b_pid", "int(5)")
                  .Define("c_pid", "int(4)")
                  .Define("e_pid", "int(11)")
                  .Define("mu_pid", "int(13)")
                  .Define("int0", "int(0)").Define("int1", "int(1)").Define("int2", "int(2)")
                  .Define("float0", "float(0)")
                  .Define("drmax1", "float(0.15)").Define("drmax2", "float(0.4)")

                  .Define("ParticlePID", {"Particle.PID"})
                  .Define("ParticlePT", {"Particle.PT"})
                  .Define("D1", {"Particle.D1"})
                  .Define("D2", {"Particle.D2"})
                  .Define("GenJetBTag", {"GenJet.BTag"}).Define("nGenbJet", "Sum(GenJetBTag)")
                  .Define("JetBTag", {"Jet.BTag"}).Define("nbJet", "Sum(JetBTag)")
                  .Define("GenMissingET_met", "GenMissingET.MET")
                  .Define("GenMissingET_eta", "GenMissingET.Eta")
                  .Define("GenMissingET_phi", "GenMissingET.Phi");
                  

    // Gen and Matching //
    auto df1 = df0.Define("isLast", ::isLast, {"Particle.PID", "Particle.D1", "Particle.D2"})
                  .Define("Top", "abs(Particle.PID) == 6 && isLast").Define("nTop", "Sum(Top)")
                  .Define("Higgs", "abs(Particle.PID) == 25 && isLast").Define("nHiggs", "Sum(Higgs)")
                  .Define("W", "abs(Particle.PID) == 24 && isLast").Define("nW", "Sum(W)")
                  .Define("GenbQuark", "abs(Particle.PID) == 5 && isLast")
                  .Define("GencQuark", "abs(Particle.PID) == 4 && isLast")

                  // Find Last Particles
                  .Define("FinalGenPart_idx", ::FirstParticle_idx, {"Particle.PID", "Particle.PT", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Top", "Higgs"})
                  .Define("Top1_idx", "FinalGenPart_idx[0]")
                  .Define("GenbFromTop1_idx", "FinalGenPart_idx[1]")
                  .Define("Top2_idx", "FinalGenPart_idx[2]")
                  .Define("GenbFromTop2_idx", "FinalGenPart_idx[3]")
                  .Define("Higgs1_idx", "FinalGenPart_idx[4]")
                  .Define("Genb1FromHiggs1_idx", "FinalGenPart_idx[5]")
                  .Define("Genb2FromHiggs1_idx", "FinalGenPart_idx[6]")
                  .Define("Higgs2_idx", "FinalGenPart_idx[7]")
                  .Define("Genb1FromHiggs2_idx", "FinalGenPart_idx[8]")
                  .Define("Genb2FromHiggs2_idx", "FinalGenPart_idx[9]")
                  .Define("bQuarkFromTop1_pt", "Particle.PT[GenbFromTop1_idx]")
                  .Define("bQuarkFromTop2_pt", "Particle.PT[GenbFromTop2_idx]")
                  .Define("b1QuarkFromHiggs1_pt", "Particle.PT[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_pt", "Particle.PT[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_pt", "Particle.PT[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_pt", "Particle.PT[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_eta", "Particle.Eta[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_eta", "Particle.Eta[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_eta", "Particle.Eta[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_eta", "Particle.Eta[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_phi", "Particle.Phi[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_phi", "Particle.Phi[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_phi", "Particle.Phi[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_phi", "Particle.Phi[Genb2FromHiggs2_idx]")
                  .Define("b1QuarkFromHiggs1_mass", "Particle.Mass[Genb1FromHiggs1_idx]")
                  .Define("b2QuarkFromHiggs1_mass", "Particle.Mass[Genb2FromHiggs1_idx]")
                  .Define("b1QuarkFromHiggs2_mass", "Particle.Mass[Genb1FromHiggs2_idx]")
                  .Define("b2QuarkFromHiggs2_mass", "Particle.Mass[Genb2FromHiggs2_idx]")
                  .Define("Q_Higgs1_var", ::GenHiggsReco, {"b1QuarkFromHiggs1_pt", "b1QuarkFromHiggs1_eta", "b1QuarkFromHiggs1_phi", "b1QuarkFromHiggs1_mass", "b2QuarkFromHiggs1_pt", "b2QuarkFromHiggs1_eta", "b2QuarkFromHiggs1_phi", "b2QuarkFromHiggs1_mass"})
                  .Define("Q_Higgs2_var", ::GenHiggsReco, {"b1QuarkFromHiggs2_pt", "b1QuarkFromHiggs2_eta", "b1QuarkFromHiggs2_phi", "b1QuarkFromHiggs2_mass", "b2QuarkFromHiggs2_pt", "b2QuarkFromHiggs2_eta", "b2QuarkFromHiggs2_phi", "b2QuarkFromHiggs2_mass"})
                  .Define("Q_Higgs1_mass", "Q_Higgs1_var[3]")
                  .Define("Q_Higgs2_mass", "Q_Higgs2_var[3]")
                  .Define("Q_Higgs_mass", ::ConcatFloat, {"Q_Higgs1_mass", "Q_Higgs2_mass"})


                  
                  // Gen bJet Matching
                  .Define("GenbJetFromTop1_idx", ::dRMatching_idx, {"GenbFromTop1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("GenbJetFromTop2_idx", ::dRMatching_idx, {"GenbFromTop2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb1JetFromHiggs1_idx", ::dRMatching_idx, {"Genb1FromHiggs1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb2JetFromHiggs1_idx", ::dRMatching_idx, {"Genb2FromHiggs1_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb1JetFromHiggs2_idx", ::dRMatching_idx, {"Genb1FromHiggs2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("Genb2JetFromHiggs2_idx", ::dRMatching_idx, {"Genb2FromHiggs2_idx", "drmax2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.Mass", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass"})
                  .Define("nGenOverlap", ::nOverlap, {"GenbJetFromTop1_idx", "GenbJetFromTop2_idx", "Genb1JetFromHiggs1_idx", "Genb2JetFromHiggs1_idx", "Genb1JetFromHiggs2_idx", "Genb2JetFromHiggs2_idx"})
                  .Define("GenOverlap_bt1", "nGenOverlap[0]")
                  .Define("GenOverlap_bt2", "nGenOverlap[1]")
                  .Define("GenOverlap_b1h1", "nGenOverlap[2]")
                  .Define("GenOverlap_b2h1", "nGenOverlap[3]")
                  .Define("GenOverlap_b1h2", "nGenOverlap[4]")
                  .Define("GenOverlap_b2h2", "nGenOverlap[5]")
                  .Define("GenMuddiness", "nGenOverlap[6]")

                  // Reco bJet Matching
                  .Define("bJetFromTop1_idx", ::dRMatching_idx, {"GenbJetFromTop1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("bJetFromTop2_idx", ::dRMatching_idx, {"GenbJetFromTop2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("b1JetFromHiggs1_idx", ::dRMatching_idx, {"Genb1JetFromHiggs1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("b2JetFromHiggs1_idx", ::dRMatching_idx, {"Genb2JetFromHiggs1_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("b1JetFromHiggs2_idx", ::dRMatching_idx, {"Genb1JetFromHiggs2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("b2JetFromHiggs2_idx", ::dRMatching_idx, {"Genb2JetFromHiggs2_idx", "drmax2", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "Jet.PT", "Jet.Eta", "Jet.Phi", "Jet.Mass"})
                  .Define("nOverlap", ::nOverlap, {"bJetFromTop1_idx", "bJetFromTop2_idx", "b1JetFromHiggs1_idx", "b2JetFromHiggs1_idx", "b1JetFromHiggs2_idx", "b2JetFromHiggs2_idx"})
                  .Define("Overlap_bt1", "nOverlap[0]")
                  .Define("Overlap_bt2", "nOverlap[1]")
                  .Define("Overlap_b1h1", "nOverlap[2]")
                  .Define("Overlap_b2h1", "nOverlap[3]")
                  .Define("Overlap_b1h2", "nOverlap[4]")
                  .Define("Overlap_b2h2", "nOverlap[5]")
                  .Define("Muddiness", "nOverlap[6]")
                  
                  .Define("nMatchedbJet", ::NumberOf, {"b1JetFromHiggs1_idx", "b2JetFromHiggs1_idx", "b1JetFromHiggs2_idx", "b2JetFromHiggs2_idx", "bJetFromTop1_idx", "bJetFromTop2_idx"})
                  .Define("nMatchedbJet_FromHiggs1", "nMatchedbJet[0]")
                  .Define("nMatchedbJet_FromHiggs2", "nMatchedbJet[1]")
                  .Define("nMatchedbJet_FromTop1", "nMatchedbJet[2]")
                  .Define("nMatchedbJet_FromTop2", "nMatchedbJet[3]")
                  .Define("nMatchedbJet_all", "nMatchedbJet[4]");

    // 4 Vector of Gen //
    auto df2 = df1.Define("Genb1JetFromHiggs1_pt", ::idx_var, {"GenJet.PT", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_pt", ::idx_var, {"GenJet.PT", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_pt", ::idx_var, {"GenJet.PT", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_pt", ::idx_var, {"GenJet.PT", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_pt", ::idx_var, {"GenJet.PT", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_pt", ::idx_var, {"GenJet.PT", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_eta", ::idx_var, {"GenJet.Eta", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_eta", ::idx_var, {"GenJet.Eta", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_eta", ::idx_var, {"GenJet.Eta", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_eta", ::idx_var, {"GenJet.Eta", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_eta", ::idx_var, {"GenJet.Eta", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_eta", ::idx_var, {"GenJet.Eta", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_phi", ::idx_var, {"GenJet.Phi", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_phi", ::idx_var, {"GenJet.Phi", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_phi", ::idx_var, {"GenJet.Phi", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_phi", ::idx_var, {"GenJet.Phi", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_phi", ::idx_var, {"GenJet.Phi", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_phi", ::idx_var, {"GenJet.Phi", "GenbJetFromTop2_idx"})

                  .Define("Genb1JetFromHiggs1_mass", ::idx_var, {"GenJet.Mass", "Genb1JetFromHiggs1_idx"})
                  .Define("Genb2JetFromHiggs1_mass", ::idx_var, {"GenJet.Mass", "Genb2JetFromHiggs1_idx"})
                  .Define("Genb1JetFromHiggs2_mass", ::idx_var, {"GenJet.Mass", "Genb1JetFromHiggs2_idx"})
                  .Define("Genb2JetFromHiggs2_mass", ::idx_var, {"GenJet.Mass", "Genb2JetFromHiggs2_idx"})
                  .Define("GenbJetFromTop1_mass", ::idx_var, {"GenJet.Mass", "GenbJetFromTop1_idx"})
                  .Define("GenbJetFromTop2_mass", ::idx_var, {"GenJet.Mass", "GenbJetFromTop2_idx"})
                  .Define("GenbJet_pt_scheme", ::pt_scheme, {"Genb1JetFromHiggs1_pt", "Genb2JetFromHiggs1_pt", "Genb1JetFromHiggs2_pt", "Genb2JetFromHiggs2_pt", "GenbJetFromTop1_pt", "GenbJetFromTop2_pt"})

                  // Higgs Reco From GenJets 
                  .Define("GenHiggs1_var", ::GenHiggsReco, {"Genb1JetFromHiggs1_pt", "Genb1JetFromHiggs1_eta", "Genb1JetFromHiggs1_phi", "Genb1JetFromHiggs1_mass", "Genb2JetFromHiggs1_pt", "Genb2JetFromHiggs1_eta", "Genb2JetFromHiggs1_phi", "Genb2JetFromHiggs1_mass"})
                  .Define("GenHiggs2_var", ::GenHiggsReco, {"Genb1JetFromHiggs2_pt", "Genb1JetFromHiggs2_eta", "Genb1JetFromHiggs2_phi", "Genb1JetFromHiggs2_mass", "Genb2JetFromHiggs2_pt", "Genb2JetFromHiggs2_eta", "Genb2JetFromHiggs2_phi", "Genb2JetFromHiggs2_mass"})
                  .Define("GenHiggs1_pt", "GenHiggs1_var[0]")
                  .Define("GenHiggs1_eta", "GenHiggs1_var[1]")
                  .Define("GenHiggs1_phi", "GenHiggs1_var[2]")
                  .Define("GenHiggs1_mass", "GenHiggs1_var[3]")
                  .Define("GenHiggs2_pt", "GenHiggs2_var[0]")
                  .Define("GenHiggs2_eta", "GenHiggs2_var[1]")
                  .Define("GenHiggs2_phi", "GenHiggs2_var[2]")
                  .Define("GenHiggs2_mass", "GenHiggs2_var[3]")
                  .Define("GenHiggs_mass", ::ConcatFloat, {"GenHiggs1_mass", "GenHiggs2_mass"});

    // Reco Matching //
    auto df3 = df2.Define("b1JetFromHiggs1_pt", ::idx_var, {"Jet.PT", "b1JetFromHiggs1_idx"})
                  .Define("b2JetFromHiggs1_pt", ::idx_var, {"Jet.PT", "b2JetFromHiggs1_idx"})
                  .Define("b1JetFromHiggs2_pt", ::idx_var, {"Jet.PT", "b1JetFromHiggs2_idx"})
                  .Define("b2JetFromHiggs2_pt", ::idx_var, {"Jet.PT", "b2JetFromHiggs2_idx"})
                  .Define("bJetFromTop1_pt", ::idx_var, {"Jet.PT", "bJetFromTop1_idx"})
                  .Define("bJetFromTop2_pt", ::idx_var, {"Jet.PT", "bJetFromTop2_idx"});
                  

    // Reco //
    auto df4 = df3.Define("goodJet", "Jet.PT>=30 && abs(Jet.Eta)<2.4")
                  .Define("goodElectron", "Electron.PT>=20 && abs(Electron.Eta)<2.4")
                  .Define("goodMuon", "Muon.PT>=20 && abs(Muon.Eta)<2.4")

                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_mass", "Jet.Mass[goodJet]")
                  .Define("Jet_E", ::GetE, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
                  .Define("Jet_btag", "Jet.BTag[goodJet]")
                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_mass", "Jet_mass[Jet_btag]")
                  .Define("bJet_E", ::GetE, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
                  .Define("bJet_size", "Sum(Jet_btag)")

                  .Define("Muon_pt", "Muon.PT[goodMuon]")
                  .Define("Muon_eta", "Muon.Eta[goodMuon]")
                  .Define("Muon_phi", "Muon.Phi[goodMuon]")
                  .Define("Muon_t", "Muon.T[goodMuon]")
                  .Define("nMuon", "Sum(goodMuon)")
                  .Define("Electron_pt", "Electron.PT[goodElectron]")
                  .Define("Electron_eta", "Electron.Eta[goodElectron]")
                  .Define("Electron_phi", "Electron.Phi[goodElectron]")
                  .Define("Electron_t", "Electron.T[goodElectron]")
                  .Define("nElectron", "Sum(goodElectron)")
                  .Define("Lep_size", "nMuon + nElectron")
                  .Define("Lep_4vec", ::TwoLeptons, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t"})
                  .Define("Lep1_pt", "Lep_4vec[0]")
                  .Define("Lep1_eta", "Lep_4vec[1]")
                  .Define("Lep1_phi", "Lep_4vec[2]")
                  .Define("Lep1_t", "Lep_4vec[3]")
                  .Define("Lep2_pt", "Lep_4vec[4]")
                  .Define("Lep2_eta", "Lep_4vec[5]")
                  .Define("Lep2_phi", "Lep_4vec[6]")
                  .Define("Lep2_t", "Lep_4vec[7]")
                  // You Must Use [bi] Before Filter // 

                  .Filter("Jet_size >= 5")
                  .Define("Jet1_pt", "Jet_pt[0]").Define("Jet1_eta", "Jet_eta[0]").Define("Jet1_phi", "Jet_phi[0]")
                  .Define("Jet1_mass", "bJet_mass[0]")
                  .Define("Jet2_pt", "Jet_pt[1]").Define("Jet2_eta", "Jet_eta[1]").Define("Jet2_phi", "Jet_phi[1]")
                  .Define("Jet2_mass", "Jet_mass[1]")
                  .Define("Jet3_pt", "Jet_pt[2]").Define("Jet3_eta", "Jet_eta[2]").Define("Jet3_phi", "Jet_phi[2]")
                  .Define("Jet3_mass", "Jet_mass[2]")
                  .Define("Jet4_pt", "Jet_pt[3]").Define("Jet4_eta", "Jet_eta[3]").Define("Jet4_phi", "Jet_phi[3]")
                  .Define("Jet4_mass", "Jet_mass[3]")
                  .Define("Jet5_pt", "Jet_pt[4]").Define("Jet5_eta", "Jet_eta[4]").Define("Jet5_phi", "Jet_phi[4]")
                  .Define("Jet5_mass", "Jet_mass[4]")

                  .Filter("bJet_size >= 4")
                  .Define("bJet1_pt", "bJet_pt[0]").Define("bJet1_eta", "bJet_eta[0]").Define("bJet1_phi", "bJet_phi[0]")
                  .Define("bJet1_mass", "bJet_mass[0]")
                  .Define("bJet2_pt", "bJet_pt[1]").Define("bJet2_eta", "bJet_eta[1]").Define("bJet2_phi", "bJet_phi[1]")
                  .Define("bJet2_mass", "bJet_mass[1]")
                  .Define("bJet3_pt", "bJet_pt[2]").Define("bJet3_eta", "bJet_eta[2]").Define("bJet3_phi", "bJet_phi[2]")
                  .Define("bJet3_mass", "bJet_mass[2]")
                  .Define("bJet4_pt", "bJet_pt[3]").Define("bJet4_eta", "bJet_eta[3]").Define("bJet4_phi", "bJet_phi[3]")
                  .Define("bJet4_mass", "bJet_mass[3]")

                  .Define("j1j2_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass"})
                  .Define("j1j3_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass"})
                  .Define("j1j4_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j1j5_dr", ::dR2, {"Jet1_pt", "Jet1_eta", "Jet1_phi", "Jet1_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j2j3_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass"})
                  .Define("j2j4_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j2j5_dr", ::dR2, {"Jet2_pt", "Jet2_eta", "Jet2_phi", "Jet2_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j3j4_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass", "Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass"})
                  .Define("j3j5_dr", ::dR2, {"Jet3_pt", "Jet3_eta", "Jet3_phi", "Jet3_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("j4j5_dr", ::dR2, {"Jet4_pt", "Jet4_eta", "Jet4_phi", "Jet4_mass", "Jet5_pt", "Jet5_eta", "Jet5_phi", "Jet5_mass"})
                  .Define("jj_dr", ::ConcatFloat_withoutSort_10, {"j1j2_dr", "j1j3_dr", "j1j4_dr", "j1j5_dr", "j2j3_dr", "j2j4_dr", "j2j5_dr", "j3j4_dr", "j3j5_dr", "j4j5_dr"}) //Be Careful
                  .Define("j_Vars", ::Vars, {"jj_dr", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_E"})
                  .Define("jj_avg_dr", "j_Vars[0]")
                  .Define("jj_max_dr", "j_Vars[1]")
                  .Define("jj_min_dr", "j_Vars[2]")
                  .Define("jj_dEta_WhenMaxdR", "j_Vars[3]")
                  .Define("j_ht", "j_Vars[4]")
                  .Define("j_cent", "j_Vars[5]")
                  .Define("jj_max_deta", "j_Vars[6]")
                  .Define("jj_max_mass", "j_Vars[7]")
                  .Define("jj_twist", "j_Vars[8]")

                  .Define("b1b2_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass"})
                  .Define("b1b3_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass"})
                  .Define("b1b4_dr", ::dR2, {"bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("b2b3_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass", "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass"})
                  .Define("b2b4_dr", ::dR2, {"bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("b3b4_dr", ::dR2, {"bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass", "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass"})
                  .Define("bb_dr", ::ConcatFloat_withoutSort_6, {"b1b2_dr", "b1b3_dr", "b1b4_dr", "b2b3_dr", "b2b4_dr", "b3b4_dr"})

                  .Define("b_Vars", ::Vars, {"bb_dr", "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_E"})
                  .Define("bb_avg_dr", "b_Vars[0]")
                  .Define("bb_max_dr", "b_Vars[1]")
                  .Define("bb_min_dr", "b_Vars[2]")
                  .Define("bb_dEta_WhenMaxdR", "b_Vars[3]")
                  .Define("b_ht", "b_Vars[4]")
                  .Define("b_cent", "b_Vars[5]")
                  .Define("bb_max_deta", "b_Vars[6]")
                  .Define("bb_max_mass", "b_Vars[7]")
                  .Define("bb_twist", "b_Vars[8]")

                   // bJet Back Tracing //  
                  .Define("bJetFrom", ::bJetFrom, {"bJet_pt", "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", "bJetFromTop1_pt", "bJetFromTop2_pt"})
                  .Define("b1", "bJetFrom[0]")
                  .Define("b2", "bJetFrom[1]")
                  .Define("b3", "bJetFrom[2]")
                  .Define("b4", "bJetFrom[3]")
                  .Define("nMatched_bFromTop", "bJetFrom[4]")
                  .Define("nMatched_bFromHiggs", "bJetFrom[5]")
                  .Define("nMatched_bJet", "bJetFrom[6]")

                  // Answer Categorizations for DNN // 
                  .Define("bCat_top_1", ::bCat_top_1, {"b1", "b2", "b3", "b4"})
                  .Define("bCat_top_2", ::bCat_top_2, {"b1", "b2", "b3", "b4"})
                  .Define("bCat_higgs_1", ::bCat_higgs_1, {"b1", "b2", "b3", "b4"})
                  .Define("bCat_higgs_2", ::bCat_higgs_2, {"b1", "b2", "b3", "b4"});
//                  .Define("bCat_3higgs", ::bCat_3higgs, {"b1", "b2", "b3", "b4"});
                  

    // Reconstruction of Higgs //
    auto df5 = df4.Define("RecoHiggs", ::RecoHiggs, {"bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass"})
                  .Define("close_Higgs_pt", "RecoHiggs[0]")
                  .Define("close_Higgs_eta", "RecoHiggs[1]")
                  .Define("close_Higgs_phi", "RecoHiggs[2]")
                  .Define("close_Higgs_mass", "RecoHiggs[3]");

    std::initializer_list<std::string> variables = {
    "ParticlePID", "ParticlePT", "D1", "D2", "JetBTag", 
    "Top1_idx", "GenbFromTop1_idx", "Top2_idx", "GenbFromTop2_idx", 
    "Higgs1_idx", "Genb1FromHiggs1_idx", "Genb2FromHiggs1_idx",
    "Higgs2_idx", "Genb1FromHiggs2_idx", "Genb2FromHiggs2_idx",
    "b1JetFromHiggs1_idx", "b2JetFromHiggs1_idx", "b1JetFromHiggs2_idx", "b2JetFromHiggs2_idx",

    "nGenbJet", "nbJet",

    "bQuarkFromTop1_pt", "bQuarkFromTop2_pt", "b1QuarkFromHiggs1_pt", "b2QuarkFromHiggs1_pt", "b1QuarkFromHiggs2_pt", "b2QuarkFromHiggs2_pt",
    "GenHiggs_mass", "GenHiggs1_mass", "GenHiggs2_mass",
    "Q_Higgs_mass", "Q_Higgs1_mass", "Q_Higgs2_mass",
    "Genb1JetFromHiggs1_pt", "Genb2JetFromHiggs1_pt", "Genb1JetFromHiggs2_pt", "Genb2JetFromHiggs2_pt", 
    "GenbJetFromTop1_pt", "GenbJetFromTop2_pt",
    "b1JetFromHiggs1_pt", "b2JetFromHiggs1_pt", "b1JetFromHiggs2_pt", "b2JetFromHiggs2_pt", 
    "bJetFromTop1_pt", "bJetFromTop2_pt",
    "GenOverlap_bt1", "GenOverlap_bt2", "GenOverlap_b1h1", "GenOverlap_b2h1", "GenOverlap_b1h2", "GenOverlap_b2h2",
    "GenMuddiness",
    "Overlap_bt1", "Overlap_bt2", "Overlap_b1h1", "Overlap_b2h1", "Overlap_b1h2", "Overlap_b2h2",
    "Muddiness",
                     
    "GenbJet_pt_scheme",
    "b1", "b2", "b3", "b4", "nMatched_bFromTop", "nMatched_bFromHiggs", "nMatched_bJet",
    "bCat_top_1", "bCat_top_2", "bCat_higgs_1", "bCat_higgs_2",

    //----------------Reco---------------------//

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
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b2b3_dr", "b2b4_dr", "b3b4_dr",

     "Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "nMuon", "Lep_size",
     "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "nElectron",
     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",

    //----------------New-----------------------//
     "bb_dr", "b_Vars", "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
    "jj_dr", "j_Vars", "jj_avg_dr", "jj_max_dr", "jj_min_dr", "j_ht", "jj_dEta_WhenMaxdR", "j_cent", "jj_max_deta", "jj_max_mass", "jj_twist", 
    "close_Higgs_pt", "close_Higgs_eta", "close_Higgs_phi", "close_Higgs_mass"
    
    };

    // MODIFY!!
    df5.Snapshot(treename, outdir+ "Full1106_" + channel + ".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
