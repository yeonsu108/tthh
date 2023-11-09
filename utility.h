
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

using floats =  ROOT::VecOps::RVec<float>;
using floatsVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<float>>;
using doubles =  ROOT::VecOps::RVec<double>;
using doublesVec =  ROOT::VecOps::RVec<ROOT::VecOps::RVec<double>>;
using ints =  ROOT::VecOps::RVec<int>;
using unints = ROOT::VecOps::RVec<unsigned int>;
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

ints SelectAddQuark(ints pid, ints m1, ints m2, ints d1, ints d2, int motherPID){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (abs(pid[i]) < 1 || abs(pid[i]) > 6) out.push_back(0);
        else if (pid[i] == pid[d1[i]] || pid[i] == pid[d2[i]]) out.push_back(0);
        else if (isFromTop(pid, m1, m2, i, motherPID)) out.push_back(0);
        else out.push_back(1);
    }
    return out;   
}

ints isLast(ints pid, ints d1, ints d2){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if ((d1[i]>=0 && pid[i]==pid[d1[i]]) || (d2[i]>=0 && pid[i]==pid[d2[i]])) out.push_back(0);
        else out.push_back(1);
    }
    return out;
}

ints make_binary(ints idx, int size){
    ints out;
    for (int i=0; i<size; i++){
        int tag = 0;
        for (int j=0; j<idx.size(); j++){
            if (idx[j] == i) tag= 1;
        }
        out.push_back(tag);
    }
    return out;
}

ints MakeLastTag(ints pid, ints d1, ints d2, int targetPID = 25){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (abs(pid[i]) != targetPID) out.push_back(0);  
        else if ((d1[i]>=0 && pid[i]==pid[d1[i]]) || (d2[i]>=0 && pid[i]==pid[d2[i]])) out.push_back(0);
        else out.push_back(1); 
        }    
    return out;
}

ints FromMother(ints pid, ints m1, ints m2, ints d1, ints d2, int targetPID, int motherPID){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (abs(pid[i]) != targetPID) out.push_back(0);
        else if ((m1[i]>=0 && abs(pid[m1[i]])==motherPID) || (m2[i]>=0 && abs(pid[m2[i]])==motherPID)) out.push_back(1);
        else out.push_back(0);
    }
    return out;   
}

ints FromMotherExact(ints pid, ints d1, ints d2, ints mother_bi){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (mother_bi[i]==0) continue;
        else {
            out.push_back(d1[i]);
            out.push_back(d2[i]);
        }
    }
    return make_binary(out, int(pid.size())); 
}

ints FromWhere(ints target, ints pid, ints m1, ints m2, ints d1, ints d2, floats pt, int order){
    ints out;
    for (int i=0; i<int(pid.size()); i++){
        if (target[i] != 1) continue;
        if (m1[i] < 0 && m2[i] >= 0) out.push_back(abs(pid[m2[i]]));
        else if (m1[i] >= 0 && m2[i] < 0) out.push_back(abs(pid[m1[i]]));
        else out.push_back(-1);
    }
    return out;
}

ints dRMatching(ints idx, floats pt1, floats eta1, floats phi1, floats m1, floats pt2, floats eta2, floats phi2, floats m2, int flag){
    ints out;
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    for (int i=0; i<int(pt1.size()); i++){
        if (idx[i] == 0) continue;
        if (flag) tmp1.SetPtEtaPhiE(pt1[i], eta1[i], phi1[i], m1[i]);
        else tmp1.SetPtEtaPhiM(pt1[i], eta1[i], phi1[i], m1[i]);
        int matched_idx = -1; float mindR = 9999;
        for (int j=0; j<int(pt2.size()); j++){
            if (pt2[j] < 20 || abs(eta2[j]) > 2.5) continue;
            tmp2.SetPtEtaPhiM(pt2[j], eta2[j], phi2[j], m2[j]);
            if (tmp1.DeltaR(tmp2) < mindR) {
                matched_idx = j;
                mindR = tmp1.DeltaR(tmp2);
            }
        }
        if (mindR > 0.4) continue;
        out.push_back(matched_idx);
    }
    return make_binary(out, int(pt2.size()));
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

int defineCategory(int naddjet, int naddbjet, int naddcjet, int naddlfjet){
    if (naddjet < 2) return -1;
    if (naddbjet > 1) return 1; //ttbb
    if (naddbjet > 0) return 2; //ttb
    if (naddcjet > 1) return 3; //ttcc
    return 0;                   //ttjj
}

floats GetE(floats pt, floats eta, floats phi, floats t){
    floats out; 
    if (int(pt.size())<=0) return out;
    auto tmp = TLorentzVector();
    for (int i=0; i<int(pt.size()); i++){
        tmp.SetPtEtaPhiE(pt[i], eta[i], phi[i], 0);
        tmp.SetT(t[i]);
        out.push_back(-tmp.M());
    }
    return out;
}

floats ConcatVector(floats v1, floats v2){
    for (const auto& element : v2) {
        v1.push_back(element);
    }
    std::sort(v1.begin(), v1.end(), std::greater<float>());

    return v1;
}

floats ConcatFloat(float f1, float f2){
    floats out;
    out.push_back(f1); out.push_back(f2);
    std::sort(out.begin(), out.end(), std::greater<float>());

    return out;
}

floats RecoHiggs(unints btag, floats pt, floats eta, floats phi, floats m){
    floats out;
    auto tmp1 = TLorentzVector();
    auto tmp2 = TLorentzVector();
    auto reco = TLorentzVector();
    for (int i=0; i<int(pt.size()); i++){
        if (btag[i]==0) continue; // select b-tagged jet1
        tmp1.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]);
        float mindR = 9999;
        for (int j=i+1; j<int(pt.size()); j++){
            if (btag[j] == 0) continue; // select b-tagged jet2
            tmp2.SetPtEtaPhiM(pt[j], eta[j], phi[j], m[j]);
            if (tmp1.DeltaR(tmp2) < mindR){  // && abs((tmp1+tmp2).M()-125) < 10){
                mindR = tmp1.DeltaR(tmp2);
                reco = tmp1 + tmp2;
            }
            if (mindR < 0.4) continue;    
        }
    }    
    out.push_back(reco.Pt());
    out.push_back(reco.Eta());
    out.push_back(reco.Phi());
    out.push_back(reco.M());
    return out;
}

ints Order(ints bi, floats pt, int order){
    ints out;
    floats tmp_pt;
    floats tmp2_pt;
    ints tmp_idx;
    for (int i=0; i<int(bi.size()); i++){
        if (bi[i] == 0) continue;
        else {
        tmp_pt.push_back(pt[i]);
        tmp2_pt.push_back(pt[i]);
        tmp_idx.push_back(i); 
        }
    }
    sort(tmp2_pt.rbegin(), tmp2_pt.rend());
    float num = tmp2_pt[order-1];
    auto it = find(tmp_pt.begin(), tmp_pt.end(), num) - tmp_pt.begin();
    out.push_back(tmp_idx[it]);

    return make_binary(out, int(bi.size()));
}
