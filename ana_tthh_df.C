
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

void ana_tthh_df(std::string channel, std::string outdir="./output/"){
    gSystem->Load("libDelphes");

    auto infile = "/data1/users/itseyes/tthh/13.6TeV/"+channel+"/Events/*.root";
    std::cout << infile << std::endl;
    std::cout << outdir << std::endl;
    auto treename = "Delphes";

    auto _df = ROOT::RDataFrame(treename, infile);

    //GenParticle Selection
    auto df1 = _df.Define("GenAddQuark_bi", ::SelectAddQuark, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2"})
                  .Define("GenAddbQuark_bi", "abs(Particle.PID == 5) && GenAddQuark_bi == 1")
                  .Define("GenAddcQuark_bi", "abs(Particle.PID == 4) && GenAddQuark_bi == 1")
                  
                  .Define("int0", "int(0)")
                  .Define("int1", "int(1)")
                  .Define("GenAddJet_bi", ::dRMatching, {"GenAddQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenAddbJet_bi", ::dRMatching, {"GenAddbQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Define("GenAddcJet_bi", ::dRMatching, {"GenAddcQuark_bi", "Particle.PT", "Particle.Eta", "Particle.Phi", "Particle.E", "GenJet.PT", "GenJet.Eta", "GenJet.Phi", "GenJet.Mass", "int1"})
                  .Redefine("GenAddcJet_bi", "GenAddcJet_bi == 1 && GenAddbJet_bi != 1")
                  .Define("GenAddlfJet_bi", "GenAddJet_bi == 1 && GenAddbJet_bi != 1 && GenAddcJet_bi != 1")

                  .Define("nGenAddQuark", "Sum(GenAddQuark_bi)")
                  .Define("nGenAddbQuark", "Sum(GenAddbQuark_bi)")
                  .Define("nGenAddcQuark", "Sum(GenAddcQuark_bi)")
                  .Define("nGenAddJet", "Sum(GenAddJet_bi)")
                  .Define("nGenAddbJet", "Sum(GenAddbJet_bi)")
                  .Define("nGenAddcJet", "Sum(GenAddcJet_bi)")
                  .Define("nGenAddlfJet", "Sum(GenAddcJet_bi)")

                  .Define("category", ::defineCategory, {"nGenAddJet", "nGenAddbJet", "nGenAddcJet", "nGenAddlfJet"})

                  .Define("LepFromTop", ::SelectWLep, {"Particle.PID", "Particle.M1", "Particle.M2", "Particle.D1", "Particle.D2"})
                  .Define("nLepFromTop", "Sum(LepFromTop)");
                  

    //object selection
    auto df2 = df1.Define("goodJet", "Jet.PT>30 && abs(Jet.Eta)<2.4")
                  .Define("goodElectron", "Electron.PT>30 && abs(Electron.Eta)<2.4")
                  .Define("goodMuon", "Muon.PT>30 && abs(Muon.Eta)<2.4")

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
                  .Define("Lepton_size", "Muon_size+Electron_size");

    auto df = df2;

    std::initializer_list<std::string> variables = {//"Event", "Jet", "Muon", "Electron",
                      "nGenAddQuark", "nGenAddbQuark", "nGenAddcQuark", // "nGenAddlfQuark",
                      "nGenAddJet", "nGenAddbJet", "nGenAddcJet", "nGenAddlfJet",
                      "GenAddQuark_bi", "GenAddbQuark_bi", "GenAddcQuark_bi",
                      "GenAddJet_bi", "GenAddbJet_bi", "GenAddcJet_bi",
                      "category", "LepFromTop", "nLepFromTop", 
                      "goodJet", "goodElectron", "goodMuon",
                      "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_btag", "Jet_size",
                      "bJet_pt", "bJet_eta", "bJet_phi", "bJet_mass", "bJet_size",
                      "Muon_pt", "Muon_eta", "Muon_phi", "Muon_e", "Muon_size",
                      "Electron_pt", "Electron_eta", "Electron_phi", "Electron_e", "Electron_size",
                      "Lepton_size",
    };
    df.Snapshot(treename, outdir+channel+".root", variables);

    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", variables, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
    //df.Snapshot<TClonesArray, TClonesArray, TClonesArray, TClonesArray>("outputTree", "out.root", {"Event", "Electron", "Muon", "Jet"}, ROOT::RDF::RSnapshotOptions("RECreate", ROOT::kZLIB, 1, 0, 99, false));
}
