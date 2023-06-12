import re
import numpy as np
import argparse as ar
import pandas as pd
import glob


class Ag:
    def __init__(self, path):
        """The object is initizialized with the agpviz file"""
        self.path = path

    def _nom(self):
        """Gets the name of the file and atom from the atomic files extension"""
        self.atom = re.split(r'\d', re.split(r'[/.\\]', self.path)[-2])[0]
        self.nom = self.path.split("_atomicfiles")[0]

    def _rdagpviz(self):
        """Reads the agpviz file, outpout: CPs (object list), CPsU (Union CPs)"""
        CPs, CPsU = [], []
        with open(self.path, mode="r") as file:
            for line in file:
                a, b = line.find("<CP of DelSqRho>"), line.find("Starting CP#")
                if a >= 0:
                    ob_temp = CP(0, 0, 0, 0, 0, 0)
                    ob_temp._rawdata(file)
                    CPs.append(ob_temp)
                if b >= 0:
                    ini, fin = int(line.split()[3]), int(file.readline().split()[3])
                    CPsU.append([ini, fin])
        return CPs, CPsU

    def _identify(self):
        """ Empirical Standars for distances between atom and CP (Various Atoms)"""
        standar = [0, 1]
        rec_dist = {"sc": [0.77, 0.9], "ti": [0.70, 0.85], "v": [0.67, 0.78], "cr": [0.65, 0.75], "mn": [
            0.60, 0.70], "fe": [0.50, 0.65], "co": [0.50, 0.65], "ni": [0.5, 0.65], "cu": [0.5, 0.65], "zn": [0.5, 0.65], "n": [0.7, 0.85], "o": [0.60, 0.74]}
        for atom, dist in rec_dist.items():
            if atom == self.atom:
                standar = dist
        return standar

    def _filter(self, CPs, dist_nucl):
        """filters CPs not in the range of distance between atom and CP"""
        CPs_n = []
        for CP in CPs:
            if dist_nucl[0] <= CP.DistFromNuc <= dist_nucl[1]:
                CPs_n.append(CP)
        return CPs_n

    def _euler(self, CPs):
        """Classifies CPs in v,e,f,checks Euler"""
        self.v, self.f, self.e = [], [], []
        for CP in CPs:
            if CP.Type == "(3,+3)":
                self.v.append(CP)
            elif CP.Type == "(3,+1)":
                self.e.append(CP)
            elif CP.Type == "(3,-1)":
                self.f.append(CP)
        # self.v.sort(key=lambda x: x.DelSqRho), self.f.sort(key=lambda x: x.DelSqRho)
        if len(self.v) + len(self.f) - len(self.e) == 2:
            print("Euler: True")
        else:
            print("Euler: False")

    def _filterU(self):
        """Filters the union list"""
        internal_num = [CP.internal_num for CP in self.e]
        for CP in self.v:
            CP.interaction = list(set(CP.interaction).intersection(internal_num))
        self.C = sum([len(CP.interaction) for CP in self.v]) / len(self.v)
        self.atomic_graph = f"[{len(self.v)}({self.C}),{len(self.e)},{len(self.f)}]"
        print(self.atomic_graph)

    def _matrix_conect(self, rows="v", columns="e"):
        """creates the matrix conectivity between rows and columns"""
        conectivity, self.matrix_conect = [0 for CP in eval("self." + str(rows))], []
        for i, CP in enumerate(eval("self." + str(rows))):
            conectivity = [0 for CP in eval("self." + str(rows))]
            for j, CP2 in enumerate(eval("self." + str(rows))):
                if CP == CP2:
                    conectivity[i] = CP.internal_num
                else:
                    conectivity[j] = list(
                        set(CP.interaction).intersection(set(CP2.interaction)))
            self.matrix_conect.append(conectivity)
        self.matrix_conect = np.array(self.matrix_conect, dtype=list)

    def _matrix_property(self, rows="v", columns="e"):
        """creates the matrix of rho"""
        self.matrix_prop = np.copy(self.matrix_conect)
        np.fill_diagonal(self.matrix_prop, np.array(
            [CP.Rho for CP in eval("self." + str(rows))]))
        for idx, x in np.ndenumerate(self.matrix_prop):
            if idx[0] != idx[1]:
                if len(x) > 0:
                    self.matrix_prop[idx] = sum(
                        [CP.Rho for CP in eval("self." + str(columns)) if any(CP.internal_num == d for d in x)])
                else:
                    self.matrix_prop[idx] = 0
        self.matrix_prop = np.array(self.matrix_prop, dtype=float)
        self.detminant_matrix = np.linalg.det(self.matrix_prop)

    def _agd(self):
        """calculates agd"""
        self.agd = sum([CP.DelSqRho for CP in self.v]) - \
            sum([CPn.DelSqRho for CPn in self.f])

    def _completCPs(self, csv_file, CPs, CPsU):
        """completes CPs for one Ag object (missing CPs)"""
        csv_df = pd.read_csv(csv_file, sep=";")
        for number, row in csv_df.iterrows():
            cp = CP(row["internal_num"], row["type"], eval(row["coords"]),
                    float(row["distfromnuc"]), float(row["rho"]), float(row["delsqrho"]))
            for number in eval(row["connectivity"]):
                CPsU.append([row["internal_num"], number])
            CPs.append(cp)
        return CPs, CPsU

    def _complete_interaction(self, csv_file, CPsU):
        """Completes gradient information"""
        csv_df = pd.read_csv(csv_file, sep=";")
        for number, row in csv_df.iterrows():
            CPsU.append([row["startCP"], row["endCP"]])
        return CPsU

    def _scale_matrix(self):
        """Scales data from the matrix"""
        ver, ed = [CP.DelSqRho for CP in self.v], [CP.DelSqRho for CP in self.e]
        # prom_v, std_v, prom_e, std_e = np.mean(ver), np.std(ver), np.mean(ed), np.std(ed)
        prom_v, prom_e = np.mean(ver), np.mean(ed)
        for idx, x in np.ndenumerate(self.matrix_prop):
            if idx[0] == idx[1]:
                self.matrix_prop[idx] = x / prom_v
            else:
                self.matrix_prop[idx] = x / prom_e
        self.detminant_matrix = np.linalg.det(self.matrix_prop)


class CP:
    def __init__(self, internal_num, Type, Coords, DistFromNuc, Rho, DelSqRho):
        self.Type = Type
        self.internal_num = internal_num
        self.Coords = Coords
        self.DistFromNuc = DistFromNuc
        self.Rho = Rho
        self.DelSqRho = DelSqRho

    def _rawdata(self, file):
        """Intended to be used with rdagpviz(Ag)"""
        self.internal_num = int(file.readline().split()[2])
        self.Type = file.readline().split()[2]
        self.Coords = file.readline().split("=")[1].split()
        self.DistFromNuc = float(file.readline().split("=")[1])
        self.Rho = float(file.readline().split("=")[1])
        self.DelSqRho = float(file.readline().split("=")[1])
        self.Coords = [float(element) for element in self.Coords]

    def _interaction(self, CPsU):
        """Interaction of all the critical points"""
        self.interaction = []
        for union in CPsU:
            if union[0] == self.internal_num:
                self.interaction.append(union[1])
            if union[1] == self.internal_num:
                self.interaction.append(union[0])


def parse_args():
    ap = ar.ArgumentParser()
    ap.add_argument("FILE", help="agpviz File")
    ap.add_argument("-i", "--idist", dest="IN_D", type=float,
                    default=None, required=False, help="Search criteria for CPs (First point)")
    ap.add_argument("-f", "--fdist", dest="FIN_D", type=float, default=None,
                    required=False, help="Search criteria for CPs (Second point)")
    ap.add_argument("-c", "--complete", dest="CCPs", type=ar.FileType(), default=None,
                    required=False, help="CVS file to complete CPs")
    ap.add_argument("-g", "--gradient", dest="C_gradient", type=ar.FileType(), default=None,
                    required=False, help="CVS file to complete gradient connectivity")
    return ap.parse_args()


def main(agpviz, filter=[None, None], complete_cp=None, complete_gradient=None, rows="v", columns="e"):
    ag_m = Ag(agpviz)
    CPs, CPsU = ag_m._rdagpviz()
    ag_m._nom()
    if filter == [None, None]:
        filtro = ag_m._identify()
    else:
        filtro = filter
    print("CP criteria:", filtro)
    if complete_cp is not None:
        CPs, CPsU = ag_m._completCPs(complete_cp, CPs, CPsU)
    CPs = ag_m._filter(CPs, filtro)
    if complete_gradient is not None:
        CPsU = ag_m._complete_interaction(complete_gradient, CPsU)
    for PC in CPs:
        PC._interaction(CPsU)
    ag_m._euler(CPs)
    ag_m._filterU()
    ag_m._matrix_conect(rows, columns)
    ag_m._matrix_property(rows, columns)
    ag_m._agd()
    return ag_m


def multiple_files(path, out_name="data"):
    """Creates the matrixes for a folder with all agpviz"""
    data = []
    for file in glob.glob(path + "/*.agpviz"):
        ag = main(file)
        print(file)
        data.append([ag.nom, ag.atom, ag.atomic_graph, ag.detminant_matrix])
    data = pd.DataFrame(data, index=None)
    data.to_csv(out_name + ".csv")


if __name__ == "__main__":
    args = parse_args()
    ag_m = main(agpviz=args.FILE, filter=[
                args.IN_D, args.FIN_D], complete_cp=args.CCPs, complete_gradient=args.C_gradient)
    print("Atom:", ag_m.atom)
    with np.printoptions(precision=5, linewidth=15 * len(ag_m.v)):
        print(ag_m.matrix_prop)
    print("Determinant:", ag_m.detminant_matrix)
    print("AGD:", ag_m.agd)
    if ag_m.C != int(ag_m.C):
        print(ag_m.matrix_conect)
