#ng histograms and save in a pdf.
# e.g. tthh

import ROOT
import numpy as np

# Criteria
DATA = "_Semi"
OS = "OS" + DATA
Tree = "Delphes"
print(OS)

# Input
tthh = "./samples1/" + OS + "_tthh.root"
ttbbbb = "./samples1/" + OS + "_ttbbbb.root"
ttbbcc = "./samples1/" + OS + "_ttbbcc.root"
ttbb = "./samples1/" + OS + "_ttbb.root"
TreeName = "Delphes"

# Luminosity [fb^-1]
L = 300

# MODIFY!!
# Cross section * Luminosity = Expected Events
# [DiLeptonic]
#x_tthh = 3.791e-5 *L * 1000
#x_ttbbbb = 0.0168 *L * 1000
#x_ttbbcc = 0.004701 *L * 1000
#x_ttbb = 0.762 *L * 1000

# [SemiLeptonic] ; *2 -> Generated sampeles partially(W+, W-) and merged.
x_tthh = 1.139e-4 *L * 1000 * 2 
x_ttbbbb = 0.049 *L * 1000 * 2
x_ttbbcc = 0.0137 *L * 1000 * 2
x_ttbb = 2.287 *L * 1000 * 2 


# RDF
tthh = ROOT.RDataFrame(Tree, tthh)
ttbbbb = ROOT.RDataFrame(Tree, ttbbbb)
ttbbcc = ROOT.RDataFrame(Tree, ttbbcc)
ttbb = ROOT.RDataFrame(Tree, ttbb)


# MODIFY!! E.S.
def Acceptance(df):
    Accept = []
    S0 = float(df.Count().GetValue())
    df = df.Filter("Muon_size == 1")
    S1 = float(df.Count().GetValue())
    df = df.Filter("Jet_size >= 5")
    S2 = float(df.Count().GetValue())
    df = df.Filter("bJet_size >= 3")
    S3 = float(df.Count().GetValue())
    Accept.extend([S0,S1,S2,S3])
    print(Accept)
    return Accept
    '''Yeild
    print("noCut : ", S0 , (S0/S0)*100, "%")
    print("S1 : ", S1, (S1/S0)*100, "%")
    print("S2 : ", S2, (S2/S0)*100, "%")
    print("S3 : ", S3, (S3/S0)*100, "%")
    '''

print("________ACCEPTANCE________")
print("tthh")    
tthh = Acceptance(tthh)
print("ttbbbb")
ttbbbb = Acceptance(ttbbbb)
print("ttbbcc")
ttbbcc = Acceptance(ttbbcc)
print("ttbb")
ttbb = Acceptance(ttbb)

Acc = {
    "tthh" : [tthh, x_tthh/tthh[0]], 
    "ttbbbb" : [ttbbbb, x_ttbbbb/ttbbbb[0]],
    "ttbbcc" : [ttbbcc, x_ttbbcc/ttbbcc[0]],
    "ttbb" : [ttbb, x_ttbb/ttbb[0]]
}

# Cutflow
def Cutflow(Acc):
    for key, value in Acc.items():
        value[0] = [element * value[1] for element in value[0]]
        print(key, value[0])
    return Acc

print("__________CUTFLOW__________")        
Acc_re = Cutflow(Acc)

print(" ")
print("________SIGNIFICANCE________")

# Significance
for i in range(0,4):
    print("Siignificance of ES :", i)
    print(Acc_re["tthh"][0][i]/np.sqrt(Acc_re["ttbbbb"][0][i]+Acc_re["ttbbcc"][0][i]+Acc_re["ttbb"][0][i]))

        




