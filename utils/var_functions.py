# Build functions make additional variables with DNN result, "pred_bcat"
import numpy as np

def higgs5_3(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass", "bJet5_mass"]
    arrow = {0:[0,1,2], 1:[0,1,3], 2:[0,1,4], 3:[0,2,3], 4:[0,2,4], 5:[0,3,4], 6:[1,2,3], 7:[1,2,4], 8:[1,3,4], 9:[2,3,4]}
    if cat == 10: return 0;

    else: 
        pt1, pt2, pt3 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]], row[pts[arrow[cat][2]]]
        eta1, eta2, eta3 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]], row[etas[arrow[cat][2]]]
        phi1, phi2, phi3 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]], row[phis[arrow[cat][2]]]
        m1, m2, m3 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]], row[ms[arrow[cat][2]]]
   
        px1, px2, px3 = pt1 * np.cos(phi1), pt2 * np.cos(phi2), pt3 * np.cos(phi3)
        py1, py2, py3 = pt1 * np.sin(phi1), pt2 * np.sin(phi2), pt3 * np.sin(phi3)
        pz1, pz2, pz3 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2), pt3 * np.sinh(eta3)
        E1,E2,E3=np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2),np.sqrt(pt3**2+m3**2)*np.cosh(eta3)
        pxs = [px1, px2, px3]; pys = [py1, py2, py3]; pzs = [pz1, pz2, pz3]; Es = [E1, E2, E3]

    close_mass = 999
    for i in range(len(pxs)):
        px_1, py_1, pz_1, E_1 = pxs[i], pys[i], pzs[i], Es[i]
        for j in range(i+1, len(pxs)):
            px_2, py_2, pz_2, E_2 = pxs[j], pys[j], pzs[j], Es[j]

            px = px_1 + px_2 
            py = py_1 + py_2 
            pz = pz_1 + pz_2 
            E = E_1 + E_2

            mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
            if abs(mass-125) < abs(close_mass-125): 
                close_mass = mass

    if close_mass != 999:
        return close_mass
    else: return 0

def higgs_5_2(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass", "bJet5_mass"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[0,4], 4:[1,2], 5:[1,3], 6:[1,4], 7:[2,3], 8:[2,4], 9:[3,4]}
    if cat == 10: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1, E2 =np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
    pxs = [px1, px2]; pys = [py1, py2]; pzs = [pz1, pz2]; Es = [E1, E2]

    close_mass = 999
    for i in range(len(pxs)):
        px_1, py_1, pz_1, E_1 = pxs[i], pys[i], pzs[i], Es[i]
        for j in range(i+1, len(pxs)):
            px_2, py_2, pz_2, E_2 = pxs[j], pys[j], pzs[j], Es[j]

            px = px_1 + px_2 
            py = py_1 + py_2 
            pz = pz_1 + pz_2 
            E = E_1 + E_2

            mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
            if abs(mass-125) < abs(close_mass-125): 
                close_mass = mass

    if close_mass != 999:
        return close_mass
    else: return 0

'''
def top_1(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass", "bJet5_mass"]
    arrow = {0:[1,2,3,4], 1:[0,2,3,4], 2:[0,1,3,4], 3:[0,1,2,4], 4:[0,1,2,3]}
    if cat == 5: return 0;

    else: 
        pt1, pt2, pt3, pt4 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]], row[pts[arrow[cat][2]]], row[pts[arrow[cat][3]]]
        eta1, eta2, eta3, eta4 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]], row[etas[arrow[cat][2]]], row[etas[arrow[cat][3]]]
        phi1, phi2, phi3, phi4 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]], row[phis[arrow[cat][2]]], row[phis[arrow[cat][3]]]
        m1, m2, m3, m4 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]], row[ms[arrow[cat][2]]], row[ms[arrow[cat][3]]]
   
        px1, px2, px3, px4 = pt1 * np.cos(phi1), pt2 * np.cos(phi2), pt3 * np.cos(phi3)
        py1, py2, py3, py4 = pt1 * np.sin(phi1), pt2 * np.sin(phi2), pt3 * np.sin(phi3)
        pz1, pz2, pz3, pz4 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2), pt3 * np.sinh(eta3)
        E1,E2,E3,E4=np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2),np.sqrt(pt3**2+m3**2)*np.cosh(eta3), np.sqrt(pt4**2+m4**2)*np.cosh(eta4)
        pxs = [px1, px2, px3, px4]; pys = [py1, py2, py3, py4]; pzs = [pz1, pz2, pz3, pz4]; Es = [E1, E2, E3, E4]

    h1_mass = 999
    for i in range(len(pxs)):
        px_1, py_1, pz_1, E_1 = pxs[i], pys[i], pzs[i], Es[i]
        for j in range(i+1, len(pxs)):
            px_2, py_2, pz_2, E_2 = pxs[j], pys[j], pzs[j], Es[j]

            px = px_1 + px_2 
            py = py_1 + py_2 
            pz = pz_1 + pz_2 
            E = E_1 + E_2

            mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
            if abs(mass-125) < abs(h1_mass-125): 
                h1_mass = mass
                idx1 = i; idx2 = j

    if h1_mass != 999: out = [h1_mass]
    else out = [0] 

    my_list = np.array([0, 1, 2, 3, 4])
    excluded_indices = [cat, idx1, idx2]
    result = np.delete(my_list, excluded_indices); k = result[0]; l = result[1]
    px_1, py_1, pz_1, E_1 = pxs[k], pys[k], pzs[k], Es[k]
    px_2, py_2, pz_2, E_2 = pxs[l], pys[l], pzs[l], Es[l]
    px = px_1 + px_2 
    py = py_1 + py_2 
    pz = pz_1 + pz_2 
    E = E_1 + E_2
    h2_mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    out.append(h2_mass)
    return out;
'''

def X_higgs(row):
    HiggsMass = 125.0; HiggsWidth = 4.07;
    mass = row["higgs_mass"]
    X_Higgs = ((HiggsMass - mass)/HiggsWidth)**2
    return X_Higgs

def bfh_dr(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt", "bJet5_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta", "bJet5_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi", "bJet5_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass", "bJet5_mass"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[1,2], 4:[1,3], 5:[2,3], 6:[]}
    if cat == 6: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
    
    out = np.sqrt( (eta1-eta2)**2 + (phi1-phi2)**2 )
    return out

def higgs_4(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass"]
    arrow = {0:[0,1], 1:[0,2], 2:[0,3], 3:[1,2], 4:[1,3], 5:[2,3], 6:[]}
    if cat == 6: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1,E2 = np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
    px = px1+px2; py = py1+py2; pz = pz1+pz2; E = E1+E2
    mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
    return mass
