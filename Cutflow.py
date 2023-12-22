#ng histograms and save in a pdf.
# e.g. tthh
import ROOT
import numpy as np

# Criteria
PRE = "CF_1216_14TeV"
Tree = "Delphes"
print(PRE)

# Input
tthh = "samples1/" + PRE + "_tthh_di.root"
ttbbbb = "samples1/" + PRE + "_ttbbbb_di.root"
ttbbcc = "samples1/" + PRE + "_ttbbcc_di.root"
ttbb = "samples1/" + PRE + "_ttbb_di.root"
tthbb = "samples1/" + PRE + "_tthbb_di.root"
TreeName = "Delphes"

# Luminosity [fb^-1]
L = 3000

# [HLLHC : DiLeptonic]
x_tthh = 4.113e-05  *L * 1000
x_ttbbbb =  4.405e-04 *L * 1000
x_ttbbcc = 6.929e-04 *L * 1000
x_ttbb =  1.303e-01 *L * 1000  # 1.315e-01 was original. - x_ttbbbb, bbcc.
x_tthbb = 2.261e-02 *L * 1000

# RDF
tthh = ROOT.RDataFrame(Tree, tthh)
ttbbbb = ROOT.RDataFrame(Tree, ttbbbb)
ttbbcc = ROOT.RDataFrame(Tree, ttbbcc)
ttbb = ROOT.RDataFrame(Tree, ttbb)
tthbb = ROOT.RDataFrame(Tree, tthbb)

print("Calculating Acceptance and Cutflow")

# MODIFY!! E.S.
def Acceptance(df, df_name):
    Accept = []
    S0 = float(df.Count().GetValue())
    df = df.Filter("Lep_size == 2 && Lep1_pt > 17 && Lep2_pt > 10"); S1 = float(df.Count().GetValue())
    df = df.Filter("Jet_size >= 5"); S2 = float(df.Count().GetValue())
    df = df.Filter("bJet_size >= 5"); S3 = float(df.Count().GetValue())
    DNN = float(input(f"{df_name} DNN Efficiency = "))
    S4 = S3 * DNN
    Accept.extend([S0,S1,S2,S3, S4])
    print(Accept)
    return Accept

print("________ACCEPTANCE________")
tthh = Acceptance(tthh, "tthh")
ttbbbb = Acceptance(ttbbbb, "ttbbbb")
ttbbcc = Acceptance(ttbbcc, "ttbbcc")
ttbb = Acceptance(ttbb, "ttbb")
tthbb = Acceptance(tthbb, "tthbb")

Acc = {
    "tthh" : [tthh, x_tthh/tthh[0]], 
    "ttbbbb" : [ttbbbb, x_ttbbbb/ttbbbb[0]],
    "ttbbcc" : [ttbbcc, x_ttbbcc/ttbbcc[0]],
    "ttbb" : [ttbb, x_ttbb/ttbb[0]],
    "tthbb" : [tthbb, x_tthbb/tthbb[0]] 
}

# Cutflow
def Cutflow(Acc):
    for key, value in Acc.items():
        value[0] = [element * value[1] for element in value[0]]
        print(key, value[0])
    return Acc

print("__________CUTFLOW__________")        
CF = Cutflow(Acc)

print(" ")
print("________SIGNIFICANCE________")

# Significance # Modify! # 
for i in range(0,5):
    print("Significance of ES :", i)
    print(CF["tthh"][0][i]/np.sqrt(CF["tthh"][0][i]+CF["ttbbbb"][0][i]+CF["ttbbcc"][0][i]+CF["ttbb"][0][i]+CF["tthbb"][0][i]))
    print("----Done----")

