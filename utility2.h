
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif
#include <vector>
#include <algorithm>
#include <map>

#include "TTree.h"
#include "TFile.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <TMath.h>

using floats =  ROOT::VecOps::RVec<float>;
using floatsVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>;
using doubles =  ROOT::VecOps::RVec<double>;
using doublesVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>;
using ints =  ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using uchars = ROOT::VecOps::RVec<unsigned char>;
using strings = ROOT::VecOps::RVec<std::string>;

using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;

bool isFromTop (ints pid, ints m1, ints m2, int idx, int motherPID=6){
    int mother = -1;
    if ( m1[idx] < 0 && m2[idx] < 0 ) return false;
    if ( (m1[idx] >= 0 && m2[idx] < 0) || (m1[idx] == m2[idx]) ) mother = m1[idx];
    else if ( m1[idx] < 0 && m2[idx] >= 0 ) mother = m2[idx];
    else{
        if ( abs(pid[m1[idx]]) == motherPID || abs(pid[m2[idx]]) == motherPID ) return true;
        else return (isFromTop(pid, m1, m2, m1[idx], motherPID) || isFromTop(pid, m1, m2, m2[idx], motherPID));
    }
    if ( abs(pid[mother]) == motherPID ) return true;
    else return isFromTop(pid, m1, m2, mother, motherPID);
} 

ints SelectAddQuark(ints pid, ints m1, ints m2, ints d1, ints d2){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (abs(pid[i]) < 1 || abs(pid[i]) > 6) out.push_back(0);
        else if (pid[i] == pid[d1[i]] || pid[i] == pid[d2[i]]) out.push_back(0);
        else if (isFromTop(pid, m1, m2, i)) out.push_back(0);
        else out.push_back(1);
    }
    return out;   
}

ints isLast (ints pid, ints d1, ints d2){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if ((d1[i]>=0 && pid[i]==pid[d1[i]]) || (d2[i]>=0 && pid[i]==pid[d2[i]])) out.push_back(0);
        else out.push_back(1);
    }
    return out;
}

int findLastIdx(int idx, ints pid, ints d1, ints d2){
    while(true){
        if (d1[idx] < 0 && d2[idx] < 0) return idx;
        if (d1[idx] >= 0 && pid[d1[idx]] == pid[idx]) idx = d1[idx];
        else if(d2[idx] >= 0 && pid[d2[idx]] == pid[idx]) idx = d2[idx];
        else return idx;
    }

}


ints FinalGenPart_idx(ints pid, ints m1, ints m2, ints d1, ints d2, ints top, ints higgs){
    ints out;
    int top_idx=-1; int bFromTop=-1; int lepFromTop=-1;
    int b_idx, w_idx;
    //std::cout << "FinalPart Functions " << std::endl;

    for (int i=0; i<int(pid.size()); i++){
        if (top[i] != 0){
            if (abs(pid[d1[i]]) == 5 && abs(pid[d2[i]]) == 24){
               b_idx = d1[i]; w_idx = findLastIdx(d2[i], pid, d1, d2);
            }
            else if (abs(pid[d1[i]]) == 24 && abs(pid[d2[i]]) == 5){
               b_idx = d2[i]; w_idx = findLastIdx(d1[i], pid, d1, d2);
            }
            top_idx = i;
            bFromTop = findLastIdx(b_idx, pid, d1, d2);
            if (abs(pid[d1[w_idx]]) % 2 == 1) {
                lepFromTop = findLastIdx(d1[w_idx], pid, d1, d2);
                }
                else {
                    lepFromTop = findLastIdx(d2[w_idx], pid, d1, d2);
                    //std::cout << " lep: " << d2[w_idx] << " pid: " << pid[d2[w_idx]] << std::endl;
                }
            out.push_back(top_idx);        // 0, 3
            out.push_back(bFromTop);   // 1, 4 
            out.push_back(lepFromTop); // 2, 5
            }
        }
    return out;
}

// Final Particle
ints FinalParticle_idx(ints pid, floats pt, ints m1, ints m2, ints d1, ints d2, ints top, ints higgs){
    ints out; 
    ints top_idx; int top1_idx=-1; int top2_idx=-1; ints bft_idx; floats bft_pt; int bft1=-1; int bft2=-1; 
    ints h_idx; ints bfh; floats bfh_pt; int h1=-1; int h2=-1; int b1fh1=-1; int b2fh1=-1; int b1fh2=-1; int b2fh2=-1;

    for (int i=0; i<int(pid.size()); i++){
        if (top[i] ==1){
            top_idx.push_back(i);
            if (abs(pid[d1[i]]) == 5 && abs(pid[d2[i]]) == 24) {bft_pt.push_back(pt[d1[i]]); bft_idx.push_back(d1[i]);} 
            else {bft_pt.push_back(pt[d2[i]]); bft_idx.push_back(d2[i]);}
        }
        if (higgs[i] ==1){
            h_idx.push_back(i);
            //std::cout << i << endl;
            bfh.push_back(d1[i]); bfh.push_back(d2[i]);
            bfh_pt.push_back(pt[d1[i]]); bfh_pt.push_back(pt[d2[i]]);
        } 
    }
    //Top
    if (bft_pt[0] < bft_pt[1]) {reverse(top_idx.begin(), top_idx.end()); reverse(bft_idx.begin(), bft_idx.end());}
    top1_idx = top_idx[0]; top2_idx = top_idx[1]; bft1 = bft_idx[0]; bft2 = bft_idx[1];
    //Higgs
    int max_idx = std::max_element(bfh_pt.begin(), bfh_pt.end()) - bfh_pt.begin();
//    std::cout << "bpt: " << bfh_pt[0] << " " << bfh_pt[1] << " " << bfh_pt[2] << " " << bfh_pt[3] << endl; 
    if (max_idx == 0){
        if (bfh_pt[2] >= bfh_pt[3]){
            h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[0]; b2fh1 = bfh[1]; b1fh2 = bfh[2]; b2fh2 = bfh[3];}
        else {h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[0]; b2fh1 = bfh[1]; b1fh2 = bfh[3]; b2fh2 = bfh[2];}}
    else if (max_idx == 1){
        if (bfh_pt[2] >= bfh_pt[3]){
            h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[1]; b2fh1 = bfh[0]; b1fh2 = bfh[2]; b2fh2 = bfh[3];}
        else {h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[1]; b2fh1 = bfh[0]; b1fh2 = bfh[3]; b2fh2 = bfh[2];}}
    else if (max_idx == 2){
        if (bfh_pt[0] >= bfh_pt[1]){
            h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[0]; b2fh2 = bfh[1];}
        else {h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[1]; b2fh2 = bfh[0];}}
    else if (max_idx == 3){
        if (bfh_pt[0] >= bfh_pt[1]) {
            h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[3]; b2fh1 = bfh[2]; b1fh2 = bfh[0]; b2fh2 = bfh[1];}
        else {h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[1]; b2fh2 = bfh[0];}}        
    else {h1 = -1; h2 = -1; b1fh1 = -1; b2fh1 = -1; b1fh2 = -1; b2fh2 = -1; std::cout << "Noooo!" << endl;}
    //std::cout << h1 << " " << h2 << " " << b1fh1 << endl << endl;
    out.push_back(top1_idx); //0
    out.push_back(findLastIdx(bft1, pid, d1, d2));
    out.push_back(top2_idx); //2
    out.push_back(findLastIdx(bft2, pid, d1, d2));
    out.push_back(h1);       //4
    out.push_back(findLastIdx(b1fh1, pid, d1, d2));
    out.push_back(findLastIdx(b2fh1, pid, d1, d2));
    out.push_back(h2);       //7 
    out.push_back(findLastIdx(b1fh2, pid, d1, d2));
    out.push_back(findLastIdx(b2fh2, pid, d1, d2));

    return out; 
}

// First Particle for Gen Particle Higgs reconstructino
ints FirstParticle_idx(ints pid, floats pt, ints m1, ints m2, ints d1, ints d2, ints top, ints higgs){
    ints out; 
    ints top_idx; int top1_idx=-1; int top2_idx=-1; ints bft_idx; floats bft_pt; int bft1=-1; int bft2=-1; 
    ints h_idx; ints bfh; floats bfh_pt; int h1=-1; int h2=-1; int b1fh1=-1; int b2fh1=-1; int b1fh2=-1; int b2fh2=-1;

    for (int i=0; i<int(pid.size()); i++){
        if (top[i] ==1){
            top_idx.push_back(i);
            if (abs(pid[d1[i]]) == 5 && abs(pid[d2[i]]) == 24) {bft_pt.push_back(pt[d1[i]]); bft_idx.push_back(d1[i]);} 
            else {bft_pt.push_back(pt[d2[i]]); bft_idx.push_back(d2[i]);}
        }
        if (higgs[i] ==1){
            h_idx.push_back(i);
            //std::cout << i << endl;
            bfh.push_back(d1[i]); bfh.push_back(d2[i]);
            bfh_pt.push_back(pt[d1[i]]); bfh_pt.push_back(pt[d2[i]]);
        } 
    }
    //Top
    if (bft_pt[0] < bft_pt[1]) {reverse(top_idx.begin(), top_idx.end()); reverse(bft_idx.begin(), bft_idx.end());}
    top1_idx = top_idx[0]; top2_idx = top_idx[1]; bft1 = bft_idx[0]; bft2 = bft_idx[1];
    //Higgs
    int max_idx = std::max_element(bfh_pt.begin(), bfh_pt.end()) - bfh_pt.begin();
//    std::cout << "bpt: " << bfh_pt[0] << " " << bfh_pt[1] << " " << bfh_pt[2] << " " << bfh_pt[3] << endl; 
    if (max_idx == 0){
        if (bfh_pt[2] >= bfh_pt[3]){
            h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[0]; b2fh1 = bfh[1]; b1fh2 = bfh[2]; b2fh2 = bfh[3];}
        else {h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[0]; b2fh1 = bfh[1]; b1fh2 = bfh[3]; b2fh2 = bfh[2];}}
    else if (max_idx == 1){
        if (bfh_pt[2] >= bfh_pt[3]){
            h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[1]; b2fh1 = bfh[0]; b1fh2 = bfh[2]; b2fh2 = bfh[3];}
        else {h1 = h_idx[0]; h2 = h_idx[1]; b1fh1 = bfh[1]; b2fh1 = bfh[0]; b1fh2 = bfh[3]; b2fh2 = bfh[2];}}
    else if (max_idx == 2){
        if (bfh_pt[0] >= bfh_pt[1]){
            h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[0]; b2fh2 = bfh[1];}
        else {h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[1]; b2fh2 = bfh[0];}}
    else if (max_idx == 3){
        if (bfh_pt[0] >= bfh_pt[1]) {
            h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[3]; b2fh1 = bfh[2]; b1fh2 = bfh[0]; b2fh2 = bfh[1];}
        else {h1 = h_idx[1]; h2 = h_idx[0]; b1fh1 = bfh[2]; b2fh1 = bfh[3]; b1fh2 = bfh[1]; b2fh2 = bfh[0];}}        
    else {h1 = -1; h2 = -1; b1fh1 = -1; b2fh1 = -1; b1fh2 = -1; b2fh2 = -1; std::cout << "Noooo!" << endl;}
    //std::cout << h1 << " " << h2 << " " << b1fh1 << endl << endl;
    out.push_back(top1_idx); //0
    out.push_back(bft1);
    out.push_back(top2_idx); //2
    out.push_back(bft2);
    out.push_back(h1);       //4
    out.push_back(b1fh1);
    out.push_back(b2fh1);
    out.push_back(h2);       //7 
    out.push_back(b1fh2);
    out.push_back(b2fh2);

    return out; 
}

float GenTopMass(float b_pt, float b_eta, float b_phi, float b_mass, float lep_pt, float lep_eta, float lep_phi, float lep_mass){
    float out;
//    auto tmp0 = TLorentzVector();
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    auto tmp = TLorentzVector();
    float m;
//    float MET_px = met * TMath::Cos(met_phi);
//    float MET_py = met * TMath::Sin(met_phi);
//    tmp0.SetPxPyPzE(MET_px, MET_py, 0, met);
    tmp1.SetPtEtaPhiM(b_pt, b_eta, b_phi, b_mass);
    tmp2.SetPtEtaPhiM(lep_pt, lep_eta, lep_phi, lep_mass);
    tmp = tmp1 + tmp2;
    out = tmp.M();
    return out; 
}

ints make_binary(ints idx, int size){
    ints out;
    for (int i=0; i<size; i++){
        int tag = 0;
        for (int j=0; j<idx.size(); j++){
            if (idx[j] == i) tag=1;
        }
        out.push_back(tag);
    }
    return out;
}

int dRMatching_idx(int idx, float drmax, floats pt1, floats eta1, floats phi1, floats m1, floats pt2, floats eta2, floats phi2, floats m2){
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    if (idx < 0) return -1;
    tmp1.SetPtEtaPhiM(pt1[idx], eta1[idx], phi1[idx], m1[idx]);
    int matched_idx = -1; float mindR = 9999;
    for (int j=0; j<int(pt2.size()); j++){
        if (pt2[j] < 20 || abs(eta2[j]) > 2.5) continue;
        if (abs((pt1[idx]-pt2[j])/pt1[idx])>0.2) continue;
//        if (pt2[j] == m2[j]) m2[j]=0;
        tmp2.SetPtEtaPhiM(pt2[j], eta2[j], phi2[j], m2[j]);
        if (tmp1.DeltaR(tmp2) < mindR) {
            matched_idx = j;
            mindR = tmp1.DeltaR(tmp2);
            // std::cout << "mindR : " << mindR << endl;
        }
    }
    if (mindR > drmax) return -1;
    return matched_idx;
}

ints dRMatching(ints idx, floats pt1, floats eta1, floats phi1, floats m1, floats pt2, floats eta2, floats phi2, floats m2){
    ints out;
    for (int i=0; i<int(pt1.size()); i++){
        if (idx[i] == 0) continue;
        int matched_idx = dRMatching_idx(i, 0.4, pt1, eta1, phi1, m1, pt2, eta2, phi2, m2);
        if (matched_idx < 0) continue;
        out.push_back(matched_idx);
    }
    return make_binary(out, int(pt2.size()));
} 

floats GetE(floats pt, floats eta, floats phi, floats m){
    floats out;
    for (int i=0; i<int(pt.size()); i++){
        auto tmp = TLorentzVector();
        tmp.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]);
        out.push_back(tmp.E());
    }
    return out;
}

floats GenHiggsReco(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2){
    floats out;
    if (pt1 <= 0 || pt2 <= 0) out = {-999,-999,-999,-999,-999};
//        std::cout << "pt: " << pt1 << " " << pt2 << "GenHiggsReco : " << out << endl ; return out;}
    auto tmp1 = TLorentzVector(); tmp1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    auto tmp2 = TLorentzVector(); tmp2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
    out.push_back((tmp1+tmp2).Pt());
    out.push_back((tmp1+tmp2).Eta());
    out.push_back((tmp1+tmp2).Phi());
    out.push_back((tmp1+tmp2).M());
    out.push_back(tmp1.DeltaR(tmp2));
    return out;
}

float dR(floats pt, floats eta, floats phi, floats m){
    float out;
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    //float deltaR;
    if (pt.size() != 2) out = -1;  // Sweep off Null events from ES at DrawHisto.
    else{
        tmp1.SetPtEtaPhiM(pt[0], eta[0], phi[0], m[0]);
        tmp2.SetPtEtaPhiM(pt[1], eta[1], phi[1], m[1]);
        out = tmp1.DeltaR(tmp2);
    }
    return out;
}

float dR2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2){
    float out;
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    tmp1.SetPtEtaPhiM(pt1, eta1, phi1, m1);
    tmp2.SetPtEtaPhiM(pt2, eta2, phi2, m2);
    out = tmp1.DeltaR(tmp2);
    return out;
}

floats RecoHiggs(floats pt, floats eta, floats phi, floats m){
    floats out;
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    auto tmp = TLorentzVector();
    float X_Higgs;
    float X_Higgs_min = 99999;
    float HiggsMass = 125.0;
    float WHiggs = 10.0;
    float dR;
    float Higgs_pt; float Higgs_eta; float Higgs_phi; float Higgs_mass;
    for (int i=0; i<int(pt.size()); i++){
        tmp1.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]);
        for (int j=i+1; j<int(pt.size()); j++){
            tmp2.SetPtEtaPhiM(pt[j], eta[j], phi[j], m[j]);
            tmp = tmp1 + tmp2; dR = tmp1.DeltaR(tmp2);
            Higgs_mass = tmp.M();
            X_Higgs = std::pow((Higgs_mass-HiggsMass)/WHiggs, 2); // + dR;
            if (X_Higgs < X_Higgs_min){
                //std::cout << dR << endl;
                X_Higgs_min = X_Higgs;
                //std::cout << i << " " << X_Higgs_min << endl; 
                Higgs_pt = tmp.Pt();
                Higgs_eta = tmp.Eta();
                Higgs_phi = tmp.Phi();
                Higgs_mass = tmp.M();
                //std::cout << Higgs_mass << endl;
            }
        }
    }
    out.push_back(Higgs_pt);
    out.push_back(Higgs_eta);
    out.push_back(Higgs_phi);
    out.push_back(Higgs_mass);
    return out;
}


floats HadTopReco(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3){
    floats out;
    auto tmp1 = TLorentzVector(); tmp1.SetPtEtaPhiM(pt1, eta1, phi1, m1); //b
    auto tmp2 = TLorentzVector(); tmp2.SetPtEtaPhiM(pt2, eta2, phi2, m2); //q1 from W
    auto tmp3 = TLorentzVector(); tmp3.SetPtEtaPhiM(pt3, eta3, phi3, m3); //q2 from W

    // W boson
    out.push_back((tmp2+tmp3).Pt());
    out.push_back((tmp2+tmp3).Eta());
    out.push_back((tmp2+tmp3).Phi());
    out.push_back((tmp2+tmp3).M());
    out.push_back(tmp2.DeltaR(tmp3));

    // Hadronic Top
    out.push_back((tmp1+(tmp2+tmp3)).Pt());
    out.push_back((tmp1+(tmp2+tmp3)).Eta());
    out.push_back((tmp1+(tmp2+tmp3)).Phi());
    out.push_back((tmp1+(tmp2+tmp3)).M());
    out.push_back(tmp1.DeltaR(tmp2+tmp3));
    return out;
}


floats ConcatFloat(float f1, float f2){
    floats out;
    out.push_back(f1); out.push_back(f2);
    std::sort(out.begin(), out.end(), std::greater<float>());

    return out;
}

floats ConcatFloat_withoutSort_6(float f1, float f2, float f3, float f4, float f5, float f6){
    floats out;
    out.push_back(f1); out.push_back(f2); out.push_back(f3); 
    out.push_back(f4); out.push_back(f5); out.push_back(f6); 
    return out;
}
floats ConcatFloat_withoutSort_10(float f1, float f2, float f3, float f4, float f5, float f6, float f7, float f8, float f9, float f10){
    floats out;
    out.push_back(f1); out.push_back(f2); out.push_back(f3); 
    out.push_back(f4); out.push_back(f5); out.push_back(f6); 
    out.push_back(f7); out.push_back(f8); out.push_back(f9); 
    out.push_back(f10);
    return out;
}

floats ConcatVector(floats v1, floats v2){
    for (const auto& element : v2) {
        v1.push_back(element);
    }
    std::sort(v1.begin(), v1.end(), std::greater<float>());

    return v1;
}

floats Concat(floats v1, floats v2){
    for (const auto& element : v2) {
        v1.push_back(element);
    }
    return v1;
}

ints NumberOf(int b1h1_idx, int b2h1_idx, int b1h2_idx, int b2h2_idx, int bt1_idx, int bt2_idx){
    ints out; int b_h1=0; int b_h2=0; int b_t1=0; int b_t2=0; int b=0;
    if (b1h1_idx>=0) b_h1 += 1; 
    if (b2h1_idx>=0) b_h1 += 1; out.push_back(b_h1);
    if (b1h2_idx>=0) b_h2 += 1; 
    if (b2h2_idx>=0) b_h2 += 1; out.push_back(b_h2);
    if (bt1_idx>=0) b_t1 += 1; out.push_back(b_t1);
    if (bt2_idx>=0) b_t2 += 1; out.push_back(b_t2);
    b = b_h1 + b_h2 + b_t1 + b_t2; out.push_back(b);
    return out;
}

float Var_Diff(float pt1, float pt2){
    float out;
    out = pt1 - pt2;
    return out;
}

floats TwoLeptons(floats pt1, floats eta1, floats phi1, floats t1, floats pt2, floats eta2, floats phi2, floats t2){
    floats out;
    floats pt = Concat(pt1, pt2);
    floats eta = Concat(eta1, eta2);
    floats phi = Concat(phi1, phi2);
    floats t = Concat(t1, t2);
    if (int(pt.size()) != 2) out = {-999, -999, -999, -999, -999, -999, -999, -999};
    else {
        if (pt[0]>pt[1]) {
            out.push_back(pt[0]); out.push_back(eta[0]); out.push_back(phi[0]); out.push_back(t[0]); 
            out.push_back(pt[1]); out.push_back(eta[1]); out.push_back(phi[1]); out.push_back(t[1]); 
        }
        else {
            out.push_back(pt[1]); out.push_back(eta[1]); out.push_back(phi[1]); out.push_back(t[1]);
            out.push_back(pt[0]); out.push_back(eta[0]); out.push_back(phi[0]); out.push_back(t[0]);
        }
    }    
    return out; 
}

ints nOverlap(int bt1, int bt2, int b1h1, int b2h1, int b1h2, int b2h2){
    ints out;
    ints idx; if (bt1!=-1) idx.push_back(bt1); if (bt2!=-1) idx.push_back(bt2); 
    if (b1h1!=-1) idx.push_back(b1h1); if (b2h1!=-1) idx.push_back(b2h1); 
    if (b1h2!=-1) idx.push_back(b1h2); if (b2h2!=-1) idx.push_back(b2h2);
    if (idx.size()==0) idx.push_back(-999);
//    std::cout << "bt1: " << bt1 << ", bt2: " << bt2 << ", b1h1: " << b1h1 << ", b2h1: " << b2h1 << ", b1h1: " << b1h2 << ", b2h2: " << b2h2 << endl;
//    std::cout << "Overlap idx vector : " << idx << endl;
    int count1 = std::count(idx.begin(), idx.end(), bt1);
    int count2 = std::count(idx.begin(), idx.end(), bt2);
    int count3 = std::count(idx.begin(), idx.end(), b1h1);
    int count4 = std::count(idx.begin(), idx.end(), b2h1);
    int count5 = std::count(idx.begin(), idx.end(), b1h2);
    int count6 = std::count(idx.begin(), idx.end(), b2h2);
    int muddiness = count1+count2+count3+count4+count5+count6;
    out.push_back(count1);
    out.push_back(count2);
    out.push_back(count3);
    out.push_back(count4);
    out.push_back(count5);
    out.push_back(count6);
    out.push_back(muddiness);
    return out;   
}

float idx_var(floats var, int idx){
    float out;
    if (idx == -1) out = -999;
    else out = var[idx];
    return out; 
}

ints pt_scheme(float b1h1, float b2h1, float b1h2, float b2h2, float bt1, float bt2){
    ints out; ints err = {-3, -3, -3, -3, -3, -3, -3};
    if (b1h1<0 || b2h1<0 || b1h2<0 || b2h2<0 || bt1<0 || bt2<0) return err; 
    floats tmp1 = {b1h1, b2h1, b1h2, b2h2, bt1, bt2} ;
    floats tmp2 = {b1h1, b2h1, b1h2, b2h2, bt1, bt2} ;
    std::cout << "pT how's it? : " << tmp2 << endl; 
    std::sort(tmp2.begin(), tmp2.end(), std::greater<float>());    
    for (int i=0; i<int(tmp1.size()); i++){
        float target = tmp1[i];
        int order = -1;
        for (int j=0; j<int(tmp2.size()); j++){
            if (tmp2[j] == target) order = j;
        }
        out.push_back(order);
    }    
    if ((out[0]<=3 && out[1]<=3) || (out[2]<=3 && out[3]<=3)) out.push_back(1); // Higgs in the scope.
    else out.push_back(-1);
    return out;
}

//////////////// bJet Categorization /////////////////////////// 
ints bJetFrom(floats bJet, float b1h1, float b2h1, float b1h2, float b2h2, float bt1, float bt2){
    ints out; int nMatched=0; int nfromtop=0; int nfromhiggs=0;
    if (int(bJet.size())<4) {out = {-2, -2, -2, -2, -2, -2, -2, -2}; return out;}  
    for (int i=0; i<4; i++){
        if (bJet[i] == b1h1 || bJet[i] == b2h1) out.push_back(1);
        else if (bJet[i] == b1h2 || bJet[i] == b2h2) out.push_back(2);
        else if (bJet[i] == bt1 || bJet[i] == bt2) out.push_back(0); // Priority for higgs, considering bkg vs sig.
        else out.push_back(-1);

        if (out.back() == 0) nfromtop += 1;
        if (out.back() == 1 || out.back() == 2) nfromhiggs += 1;
        if (out.back() != -1) nMatched += 1;
    }
    out.push_back(nfromtop);
    out.push_back(nfromhiggs);
    out.push_back(nMatched);
//    std::cout << out << endl << endl;
    return out;
}

int bCat_top_1(int b1, int b2, int b3, int b4){
    int out;
    if (b1 == 0 && b2 == 0) out = 0;
    else if (b1 == 0 && b3 == 0) out=1;
    else if (b1 == 0 && b4 == 0) out=2;
    else if (b2 == 0 && b3 == 0) out=3;
    else if (b2 == 0 && b4 == 0) out=4;
    else if (b3 == 0 && b4 == 0) out=5;
    else if (b1 == 0) out=6;
    else if (b2 == 0) out=7;
    else if (b3 == 0) out=8;
    else if (b4 == 0) out=9;
    else out=10;
    return out;
}
int bCat_top_2(int b1, int b2, int b3, int b4){
    int out;
    if (b1 == 0 && b2 == 0) out = 0;
    else if (b1 == 0 && b3 == 0) out=1;
    else if (b1 == 0 && b4 == 0) out=2;
    else if (b2 == 0 && b3 == 0) out=3;
    else if (b2 == 0 && b4 == 0) out=4;
    else if (b3 == 0 && b4 == 0) out=5;
    else out=6;
    return out;
}

int bCat_higgs_1(int b1, int b2, int b3, int b4){
    int out;
    if (b1 == 2) b1 -= 1;
    if (b2 == 2) b2 -= 1;
    if (b3 == 2) b3 -= 1;
    if (b4 == 2) b4 -= 1;

    if (b2 == 1 && b3 == 1) out = 0;
    else if (b2 == 1 && b4 == 1) out=1;
    else if (b3 == 1 && b4 == 1) out=2;
    else if (b2 == 1) out=3;
    else if (b3 == 1) out=4;
    else if (b4 == 1) out=5;
    else out=6;
    return out;
}
int bCat_higgs_2(int b1, int b2, int b3, int b4){
    int out;
    if (b1 == 2) b1 -= 1;
    if (b2 == 2) b2 -= 1;
    if (b3 == 2) b3 -= 1;
    if (b4 == 2) b4 -= 1;

    if (b2 == 1 && b3 == 1) out = 0;
    else if (b2 == 1 && b4 == 1) out=1;
    else if (b3 == 1 && b4 == 1) out=2;
    else out=3;
    return out;
}

//////////////// bJet Categorization /////////////////////////// 
//////////////// bJet Categorization /////////////////////////// 

floats Vars(floats dr, floats pt, floats eta, floats phi, floats mass, floats E){
    floats out; 
    std::map<int, ints> arrow = {{0, {0,1}}, {1, {0,2}}, {2, {0,3}}, {3, {1,2}}, {4, {1,3}}, {5, {2,3}}};
    int max_idx = std::max_element(dr.begin(), dr.end()) - dr.begin(); float max_dr = dr[max_idx];
    int min_idx = std::min_element(dr.begin(), dr.end()) - dr.begin(); float min_dr = dr[min_idx];
 
    float sum_dr = 0;
    int num_dr = 0;
    for (float tmp : dr){
        if (abs(tmp) < 10) {sum_dr += tmp; num_dr += 1;}
    }
    float avg_dr = sum_dr/num_dr;
    out.push_back(avg_dr); // [0] avg_dr 
    out.push_back(max_dr); // [1] max_dr
    out.push_back(min_dr); // [2] min_dr

    int eta_idx1 = arrow[max_idx][0]; int eta_idx2 = arrow[max_idx][1]; // * VALID AT 4JET CASE ONLY *
    out.push_back(abs(eta[eta_idx1] - eta[eta_idx2])); // [3] dEta_WhenMaxdR * VALID AT 4JET CASE ONLY *

    float ht;
    float sum_pt = 0;
    int num_pt = 0;
    for (float tmp : pt){
        if (tmp > 0) {sum_pt += tmp; num_pt += 1;} 
    }
    ht = sum_pt; 
    out.push_back(ht); // [4] ht, scalar sum of pt. 

    float sum_E = 0; 
    for (float tmp : E){
        if (tmp > 0) sum_E += tmp;
    }
    float cent = ht/sum_E;
    out.push_back(cent); // [5] Centrality : ht/sum(E)

    float max_deta = 999;
    for (int i=0; i<int(eta.size()); i++){
        float eta1 = eta[i];
        for (int j=i+1; j<int(eta.size()); j++){
            float eta2 = eta[j];
            float deta = abs(eta1 - eta2);
            if (deta < max_deta) max_deta = deta;
        }
    }
    if (max_deta < 10) out.push_back(max_deta); // [6] Max dEta;
    else out.push_back(-999);

    float twist;
    float tmp_mass; float max_mass=0;
    auto tmp_v1 = TLorentzVector();
    auto tmp_v2 = TLorentzVector();
    auto tmp_v = TLorentzVector();
    for (int i=0; i<int(pt.size()); i++){
        tmp_v1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
        for (int j=i+1; j<int(pt.size()); j++){
            tmp_v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
            tmp_v = tmp_v1 + tmp_v2; tmp_mass = tmp_v.M();
            if (tmp_mass > max_mass){ 
                max_mass = tmp_mass;
                float tdphi = phi[i] - phi[j];
                float tdeta = eta[i] - eta[j];
                twist = std::atan2(tdphi, tdeta);
            }
        }
    }
    if (max_mass > 0) out.push_back(max_mass); // [7] Biggest Mass Combination
    else out.push_back(0);
    out.push_back(abs(twist)); // [8] twist angle @ Biggest Mass Combination
    
    return out;
}
