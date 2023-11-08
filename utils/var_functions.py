# Build functions make additional variables with DNN result, "pred_bcat"
import numpy as np

def higgs_top_1(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass"]
    arrow = {0:[2,3], 1:[1,3], 2:[1,2], 3:[0,3], 4:[0,2], 5:[0,1], 6:[1,2,3], 7:[0,2,3], 8:[0,1,3], 9:[0,1,2], 10:[]}
    if cat == 10: return 0;

    if cat <= 5:
        pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
        eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
        phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
        m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
        px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
        py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
        pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
        E1,E2 = np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
        pxs = [px1, px2]; pys = [py1, py2]; pzs = [pz1, pz2]; Es = [E1, E2]

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
                
def higgs_top_2(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass"]
    arrow = {0:[2,3], 1:[1,3], 2:[1,2], 3:[0,3], 4:[0,2], 5:[0,1], 6:[]}
    if cat == 6: return 0;

    pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
    eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
    phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
    m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
    px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
    py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
    pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
    E1,E2 = np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
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

def higgs_higgs_1(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass"]
    arrow = {0:[0,1,2], 1:[0,1,3], 2:[0,2,3], 3:[0,1], 4:[0,2], 5:[0,3], 6:[]}
    if cat == 6: return 0;

    if cat >= 3:
        pt1, pt2 = row[pts[arrow[cat][0]]], row[pts[arrow[cat][1]]]
        eta1, eta2 = row[etas[arrow[cat][0]]], row[etas[arrow[cat][1]]]
        phi1, phi2 = row[phis[arrow[cat][0]]], row[phis[arrow[cat][1]]]
        m1, m2 = row[ms[arrow[cat][0]]], row[ms[arrow[cat][1]]]
   
        px1, px2 = pt1 * np.cos(phi1), pt2 * np.cos(phi2)
        py1, py2 = pt1 * np.sin(phi1), pt2 * np.sin(phi2)
        pz1, pz2 = pt1 * np.sinh(eta1), pt2 * np.sinh(eta2)
        E1,E2 = np.sqrt(pt1**2+m1**2)*np.cosh(eta1),np.sqrt(pt2**2+m2**2)*np.cosh(eta2)
        pxs = [px1, px2]; pys = [py1, py2]; pzs = [pz1, pz2]; Es = [E1, E2]

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

def higgs_higgs_2(row):
    cat = row["pred_bcat"]
    pts = ["bJet1_pt", "bJet2_pt", "bJet3_pt", "bJet4_pt"]
    etas = ["bJet1_eta", "bJet2_eta", "bJet3_eta", "bJet4_eta"]
    phis = ["bJet1_phi", "bJet2_phi", "bJet3_phi", "bJet4_phi"]
    ms = ["bJet1_mass", "bJet2_mass", "bJet3_mass", "bJet4_mass"]
    arrow = {0:[0,1,2], 1:[0,1,3], 2:[0,2,3], 3:[]}
    if cat == 3: return 0;

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

