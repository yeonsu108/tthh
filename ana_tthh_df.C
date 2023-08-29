
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

#include "utility.h"

// MODIFY!!
void ana_tthh_df(std::string channel, std::string outdir="./samples1/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/*.root";
    std::cout << infile << std::endl;
    std::cout << outdir << std::endl;
    auto treename = "Delphes";

    auto _df = ROOT::RDataFrame(treename, infile);
    
    // Constants, Basic Branches //
    auto df0 = _df.Define("TopPID", "int(6)")
                  .Define("HiggsPID", "int(25)")
                  .Define("W_pid", "int(24)")
                  .Define("b_pid", "int(5)")
                  .Define("Electron_pid", "int(11)")
                  .Define("Muon_pid", "int(13)")
                  .Define("GluonPID", "int(21)")
                  .Define("int0", "int(0)")
                  .Define("int1", "int(1)")
                  .Define("int2", "int(2)")
                  .Define("isLast", ::isLast, {"Particle.PID", "Particle.D1", "Particle.D2"})
                  .Define("ParticlePID", {"Particle.PID"})
                  .Define("ParticlePT", {"Particle.PT"})
                  .Define("D1", {"Particle.D1"})
                  .Define("D2", {"Particle.D2"})
                  .Define("JetBTag", {"Jet.BTag"});

    // Gen_bi //   // Quark, FromQuark, AddQuark
    auto df1 = df0.Define("Top", "abs(Particle.PID) == 6 && isLast")
                  .Define("Higgs", ::MakeLastTag, {"Particle.PID","Particle.D1", "Particle.D2", "HiggsPID"})
                  .Define("Higgs1", ::Order, {"Higgs", "Particle.PT", "int1"})
                  .Define("Higgs2", ::Order, {"Higgs", "Particle.PT", "int2"})
                  .Define("W", "abs(Particle.PID) == 24 && isLast")
                  .Define("Gluon", "abs(Particle.PID) == 21 && isLast")
                  .Define("GenbQuark_bi", "abs(Particle.PID) == 5 && isLast")
                  .Define("GencQuark_bi", "abs(Particle.PID) == 4 && isLast")

                  .Define("WFromTop", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "W_pid", "TopPID"})
                  .Define("bFromTop", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "b_pid", "TopPID"})
                  .Define("bFromHiggs", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "b_pid", "HiggsPID"})
                  .Define("bFromHiggs1", ::FromMotherExact, {"Particle.PID", "Particle.D1", "Particle.D2", "Higgs1"})
                  .Define("bFromHiggs2", ::FromMotherExact, {"Particle.PID", "Particle.D1", "Particle.D2", "Higgs2"})
                  .Define("bFromGluon", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "b_pid", "GluonPID"})
                  .Define("HiggsFromTop", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "HiggsPID", "TopPID"})
                  .Define("ElectronFromW", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Electron_pid", "W_pid"})
                  .Define("MuonFromW", ::FromMother, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "Muon_pid", "W_pid"})

                  .Define("GenAddQuark_bi", ::SelectAddQuark, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "TopPID"})
                  .Define("GenAddbQuark_bi", "abs(Particle.PID) == 5 && GenAddQuark_bi == 1")
                  .Define("GenAddcQuark_bi", "abs(Particle.PID) == 4 && GenAddQuark_bi == 1")

                  // Jet, FromJet, AddJet <- Quark
                  .Define("GenbJet_bi", ::dRMatching, {"GenbQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GencJet_bi", ::dRMatching, {"GencQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenbJetFromTop_bi", ::dRMatching, {"bFromTop", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenbJetFromHiggs_bi", ::dRMatching, {"bFromHiggs", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenbJetFromHiggs1_bi", ::dRMatching, {"bFromHiggs1", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenbJetFromHiggs2_bi", ::dRMatching, {"bFromHiggs2", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenAddJet_bi", ::dRMatching, {"GenAddQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenAddbJet_bi", ::dRMatching, {"GenAddbQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenAddcJet_bi", ::dRMatching, {"GenAddcQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Redefine("GenAddcJet_bi", "GenAddcJet_bi == 1 && GenAddbJet_bi != 1")
                  .Define("GenAddlfJet_bi", "GenAddJet_bi == 1 && GenAddbJet_bi != 1 && GenAddcJet_bi != 1");                                  

    // nGen // 
    auto df2 = df1.Define("nTop", "Sum(Top)")
                  .Define("nHiggs", "Sum(Higgs)")
                  .Define("nHiggs1", "Sum(Higgs1)")
                  .Define("nHiggs2", "Sum(Higgs2)")
                  .Define("nW", "Sum(W)")
                  .Define("nGluon", "Sum(Gluon)")
                  .Define("nGenbQuark", "Sum(GenbQuark_bi)")
                  .Define("nGencQuark", "Sum(GencQuark_bi)")

                  .Define("nWFromTop", "Sum(WFromTop)")
                  .Define("nbFromTop", "Sum(bFromTop)")
                  .Define("nbFromHiggs", "Sum(bFromHiggs)")
                  .Define("nbFromHiggs1", "Sum(bFromHiggs1)")
                  .Define("nbFromHiggs2", "Sum(bFromHiggs2)")
                  .Define("nbFromGluon", "Sum(bFromGluon)")
                  .Define("nHiggsFromTop", "Sum(HiggsFromTop)")
                  .Define("nElectronFromW", "Sum(ElectronFromW)")
                  .Define("nMuonFromW", "Sum(MuonFromW)")
                  .Define("nLepFromW", "nElectronFromW+nMuonFromW")

                  .Define("nGenAddQuark", "Sum(GenAddQuark_bi)")
                  .Define("nGenAddbQuark", "Sum(GenAddbQuark_bi)")
                  .Define("nGenAddcQuark", "Sum(GenAddcQuark_bi)")
                    
                  // Jet
                  .Define("nGenbJet", "Sum(GenbJet_bi)")
                  .Define("nGencJet", "Sum(GencJet_bi)")
                  .Define("nGenbJetFromTop", "Sum(GenbJetFromTop_bi)")
                  .Define("nGenbJetFromHiggs", "Sum(GenbJetFromHiggs_bi)")
                  .Define("nGenbJetFromHiggs1", "Sum(GenbJetFromHiggs1_bi)")
                  .Define("nGenbJetFromHiggs2", "Sum(GenbJetFromHiggs2_bi)")
                  .Define("nGenAddJet", "Sum(GenAddJet_bi)")
                  .Define("nGenAddbJet", "Sum(GenAddbJet_bi)")
                  .Define("nGenAddcJet", "Sum(GenAddcJet_bi)")
                  .Define("nGenAddlfJet", "Sum(GenAddcJet_bi)")

                  .Define("category", ::defineCategory, {"nGenAddJet", "nGenAddbJet", "nGenAddcJet", "nGenAddlfJet"})
                  .Define("ttbbbb_val", "nTop ==2 && nGenAddbQuark == 4")
                  .Define("ttbbcc_val", "nTop ==2 && nGenAddcQuark == 2 && nGenAddbQuark == 2");

    // Gen 4-vector //
    auto df3 = df2.Define("GenTop_mass", "Particle.Mass[Top]")
                  .Define("GenHiggs_mass", "Particle.Mass[Higgs]")

                  .Define("GenbQuark_pt", "Particle.PT[GenbQuark_bi]")
                  .Define("GenbQuark_eta", "Particle.Eta[GenbQuark_bi]")
                  .Define("GenbQuark_phi", "Particle.Phi[GenbQuark_bi]")
                  .Define("GenElectronFromW_pt", "Particle.PT[ElectronFromW]")
                  .Define("GenElectronFromW_eta", "Particle.Eta[ElectronFromW]")
                  .Define("GenElectronFromW_phi", "Particle.Phi[ElectronFromW]")
                  .Define("GenMuonFromW_pt", "Particle.PT[MuonFromW]")
                  .Define("GenMuonFromW_eta", "Particle.Eta[MuonFromW]")
                  .Define("GenMuonFromW_phi", "Particle.Phi[MuonFromW]")

                  .Define("GenbQuarkFromTop_pt", "Particle.PT[bFromTop]")
                  .Define("GenbQuarkFromTop_eta", "Particle.Eta[bFromTop]")
                  .Define("GenbQuarkFromTop_phi", "Particle.Phi[bFromTop]")
                  .Define("GenbQuarkFromHiggs_pt", "Particle.PT[bFromHiggs]")
                  .Define("GenbQuarkFromHiggs_eta", "Particle.Eta[bFromHiggs]")
                  .Define("GenbQuarkFromHiggs_phi", "Particle.Phi[bFromHiggs]")

                  .Define("GenAddbQuark_pt", "Particle.PT[GenAddbQuark_bi]")
                  .Define("GenAddbQuark_eta", "Particle.Eta[GenAddbQuark_bi]")
                  .Define("GenAddbQuark_phi", "Particle.Phi[GenAddbQuark_bi]")
                    
                  // bJet
                  .Define("GenJet_pt", "GenJet.PT")
                  .Define("GenbJet_pt", "GenJet.PT[GenbJet_bi]")
                  .Define("GenbJet_eta", "GenJet.Eta[GenbJet_bi]")
                  .Define("GenbJet_phi", "GenJet.Phi[GenbJet_bi]")

                  .Define("GenbJetFromTop_pt", "GenJet.PT[GenbJetFromTop_bi]")
                  .Define("GenbJetFromTop_eta", "GenJet.Eta[GenbJetFromTop_bi]")
                  .Define("GenbJetFromTop_phi", "GenJet.Phi[GenbJetFromTop_bi]")
                  .Define("GenbJetFromTop_mass", "GenJet.Mass[GenbJetFromTop_bi]")
                  .Define("GenbJetFromHiggs_pt", "GenJet.PT[GenbJetFromHiggs_bi]")
                  .Define("GenbJetFromHiggs_eta", "GenJet.Eta[GenbJetFromHiggs_bi]")
                  .Define("GenbJetFromHiggs_phi", "GenJet.Phi[GenbJetFromHiggs_bi]")
                  .Define("GenbJetFromHiggs_mass", "GenJet.Mass[GenbJetFromHiggs_bi]")
                  .Define("GenbJetFromHiggs1_pt", "GenJet.PT[GenbJetFromHiggs1_bi]")//
                  .Define("GenbJetFromHiggs1_eta", "GenJet.Eta[GenbJetFromHiggs1_bi]")
                  .Define("GenbJetFromHiggs1_phi", "GenJet.Phi[GenbJetFromHiggs1_bi]")
                  .Define("GenbJetFromHiggs1_mass", "GenJet.Mass[GenbJetFromHiggs1_bi]")
                  .Define("GenbJetFromHiggs2_pt", "GenJet.PT[GenbJetFromHiggs2_bi]")
                  .Define("GenbJetFromHiggs2_eta", "GenJet.Eta[GenbJetFromHiggs2_bi]")
                  .Define("GenbJetFromHiggs2_phi", "GenJet.Phi[GenbJetFromHiggs2_bi]")
                  .Define("GenbJetFromHiggs2_mass", "GenJet.Mass[GenbJetFromHiggs2_bi]")//

                  .Define("GenAddJet_pt", "GenJet.PT[GenAddJet_bi]")
                  .Define("GenAddJet_eta", "GenJet.Eta[GenAddJet_bi]")
                  .Define("GenAddJet_phi", "GenJet.Phi[GenAddJet_bi]")
                  .Define("GenAddbJet_pt", "GenJet.PT[GenAddbJet_bi]")
                  .Define("GenAddbJet_eta", "GenJet.Eta[GenAddbJet_bi]")
                  .Define("GenAddbJet_phi", "GenJet.Phi[GenAddbJet_bi]");
                     
    // Gen dR //
    auto df4 = df3.Define("GenbJetFromTop_dr", ::dR, {"GenbJetFromTop_pt", "GenbJetFromTop_eta", "GenbJetFromTop_phi", "GenbJetFromTop_mass"})
                  .Define("GenbJetFromHiggs1_dr", ::dR, {"GenbJetFromHiggs1_pt", "GenbJetFromHiggs1_eta", "GenbJetFromHiggs1_phi", "GenbJetFromHiggs1_mass"})
                  .Define("GenbJetFromHiggs2_dr", ::dR, {"GenbJetFromHiggs2_pt", "GenbJetFromHiggs2_eta", "GenbJetFromHiggs2_phi", "GenbJetFromHiggs2_mass"})
                  .Define("GenbJetFromHiggs_dr", ::ConcatFloat, {"GenbJetFromHiggs1_dr", "GenbJetFromHiggs2_dr"})
                  .Define("HiggsFromWhere", ::FromWhere, {"Higgs", "Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2", "GenbQuark_pt", "int0"});

    //Reco //      // OS //
    auto df5 = df4.Define("goodJet", "Jet.PT>=40 && abs(Jet.Eta)<2.4")
                  .Define("goodElectron", "Electron.PT>=25 && abs(Electron.Eta)<2.4")
                  .Define("goodMuon", "Muon.PT>=25 && abs(Muon.Eta)<2.4")

                  .Define("Jet_pt", "Jet.PT[goodJet]")
                  .Define("Jet_eta", "Jet.Eta[goodJet]")
                  .Define("Jet_phi", "Jet.Phi[goodJet]")
                  .Define("Jet_mass", "Jet.Mass[goodJet]")
                  .Define("Jet_btag", "Jet.BTag[goodJet]")
                  .Redefine("Jet_size", "Sum(goodJet)")

                  .Define("bJet_pt", "Jet_pt[Jet_btag]")
                  .Define("bJet_eta", "Jet_eta[Jet_btag]")
                  .Define("bJet_phi", "Jet_phi[Jet_btag]")
                  .Define("bJet_mass", "Jet_mass[Jet_btag]")
                  .Define("bJet_size", "Sum(Jet_btag)")

                  .Define("Muon_pt", "Muon.PT[goodMuon]")
                  .Define("Muon_eta", "Muon.Eta[goodMuon]")
                  .Define("Muon_phi", "Muon.Phi[goodMuon]")
                  .Define("Muon_t", "Muon.T[goodMuon]")
                  .Define("Muon_e", ::GetE, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_t"})
                  .Redefine("Muon_size", "Sum(goodMuon)")

                  .Define("Electron_pt", "Electron.PT[goodElectron]")
                  .Define("Electron_eta", "Electron.Eta[goodElectron]")
                  .Define("Electron_phi", "Electron.Phi[goodElectron]")
                  .Define("Electron_t", "Electron.T[goodElectron]")
                  .Define("Electron_e", ::GetE, {"Electron_pt", "Electron_eta", "Electron_phi", "Electron_t"})
                  .Redefine("Electron_size", "Sum(goodElectron)")
                  .Define("Lepton_size", "Muon_size+Electron_size")
                  .Define("Lepton_pt", ::ConcatVector, {"Muon_pt", "Electron_pt"})
                  
                  .Define("Higgs_vars", ::RecoHiggs, {"Jet_btag", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
                  .Define("Higgs_pt", "Higgs_vars[0]")
                  .Define("Higgs_eta", "Higgs_vars[1]")
                  .Define("Higgs_phi", "Higgs_vars[2]")
                  .Define("Higgs_mass", "Higgs_vars[3]");
                  
// ~df4 is for GenStudy.
    auto df = df5;   

    std::initializer_list<std::string> variables = {
                      // Basic Branches // 
                      "ParticlePID", "ParticlePT", "D1", "D2", "JetBTag",

                      // Gen Quark_bi //
                      "Top", "Higgs", "Higgs1", "Higgs2", "W", "Gluon", "GenbQuark_bi", "GencQuark_bi",

                      "WFromTop", "bFromTop", "bFromHiggs", "bFromHiggs1", "bFromHiggs2", "bFromGluon", 
                      "HiggsFromTop", "ElectronFromW", "MuonFromW",

                      "GenAddbQuark_bi", "GenAddcQuark_bi",

                      // Gen Jet_bi //
                      "GenbJet_bi", "GencJet_bi", 

                      "GenbJetFromTop_bi", "GenbJetFromHiggs_bi", "GenbJetFromHiggs1_bi", "GenbJetFromHiggs2_bi",

                      "GenAddJet_bi", "GenAddbJet_bi", "GenAddcJet_bi", "GenAddlfJet_bi",

                      // nGen Quark//
                      "nTop", "nHiggs", "nHiggs1", "nHiggs2", "nW", "nGluon", "nGenbQuark", "nGencQuark",

                      "nWFromTop", "nbFromTop", "nbFromHiggs", "nbFromHiggs1", "nbFromHiggs2", 
                      "nbFromGluon", "nHiggsFromTop", "nElectronFromW", "nMuonFromW", "nLepFromW",

                      "nGenAddQuark", "nGenAddbQuark", "nGenAddcQuark",

                      // nGen Jet //
                      "nGenbJet", "nGencJet",
                      
                      "nGenbJetFromTop", "nGenbJetFromHiggs", "nGenbJetFromHiggs1", "nGenbJetFromHiggs2", 

                      "nGenAddJet", "nGenAddbJet", "nGenAddcJet", "nGenAddlfJet",

                      "category", "ttbbbb_val", "ttbbcc_val",

                      // Gen Quark 4-vector //
                      "GenTop_mass", "GenHiggs_mass",
                      "GenbQuark_pt", "GenbQuark_eta", "GenbQuark_phi",
                      "GenElectronFromW_pt", "GenElectronFromW_eta", "GenElectronFromW_phi",
                      "GenMuonFromW_pt", "GenMuonFromW_eta", "GenMuonFromW_phi",

                      "GenbQuarkFromTop_pt", "GenbQuarkFromTop_eta", "GenbQuarkFromTop_phi",
                      "GenbQuarkFromHiggs_pt", "GenbQuarkFromHiggs_eta", "GenbQuarkFromHiggs_phi",

                      "GenAddbQuark_pt", "GenAddbQuark_eta", "GenAddbQuark_phi",

                      // Gen Jet 4-vector //
                      "GenJet_pt", "GenbJet_pt", "GenbJet_eta", "GenbJet_phi",

                      "GenbJetFromTop_pt", "GenbJetFromTop_eta", "GenbJetFromTop_phi", "GenbJetFromTop_mass",
                      "GenbJetFromHiggs_pt", "GenbJetFromHiggs_eta", "GenbJetFromHiggs_phi", "GenbJetFromHiggs_mass",
                      "GenbJetFromHiggs1_pt", "GenbJetFromHiggs1_eta", "GenbJetFromHiggs1_phi", "GenbJetFromHiggs1_mass",
                      "GenbJetFromHiggs2_pt", "GenbJetFromHiggs2_eta", "GenbJetFromHiggs2_phi", "GenbJetFromHiggs2_mass",
                      
                      "GenAddJet_pt", "GenAddJet_eta", "GenAddJet_phi", "GenAddbJet_pt", "GenAddbJet_eta", "GenAddbJet_phi",
                        
                      // Gen dR //
                      "GenbJetFromTop_dr", "GenbJetFromHiggs1_dr", "GenbJetFromHiggs2_dr", "GenbJetFromHiggs_dr",
                      "HiggsFromWhere",



                      //-------------------------------------- Reco----------------------------------------//
                      
                      // Object Selection //
                      "goodJet", "goodElectron", "goodMuon",

                      // Objects //
                      "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_btag", "Jet_size",

                      "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_size",

                      "Muon_pt", "Muon_eta", "Muon_phi", "Muon_t", "Muon_e", "Muon_size",

                      "Electron_pt", "Electron_eta", "Electron_phi", "Electron_t", "Electron_e", "Electron_size",

                      "Lepton_size", "Lepton_pt",

                      "Higgs_vars", "Higgs_pt", "Higgs_eta", "Higgs_phi", "Higgs_mass"


    };

    // MODIFY!!
    df.Snapshot(treename, outdir+ "OS" + "_" +channel+".root", variables); 
    std::cout << "done" << std::endl; 

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
