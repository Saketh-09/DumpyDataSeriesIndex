import os
import threading
import numpy as np

# Constants
fuzzy_number = 0
mask = None

class Node:
    def __init__(self, layer=0, partition_id=-1, size=0):
        self.layer = layer
        self.partition_id = partition_id
        self.size = size
        self.ch = []
        self.offsets = []

    def generateSaxAndCardIn1stLayer(self, nav_id):
        pass

    def generateSaxAndCardinality(self, child, index):
        pass

    def generateSaxAndCardIn1stLayer4LeafNode(self, nav_id):
        pass

    def generateSaxAndCardinality4LeafNode(self, child, index):
        pass

    def route(self, saxes):
        pass

class NODE_RECORDER:
    def __init__(self, actual_size=0, series_index_list=None):
        self.actual_size = actual_size
        self.series_index_list = series_index_list if series_index_list is not None else []

class Const:
    segmentNum = 0
    th = 0
    delta = 0
    fbl_series_num = 0
    tsLength = 0
    vertexNum = 0
    fuzzyidxfn = ""

class partUnit:
    def __init__(self, id=0, s=0, process_id=-1):
        self.id = id
        self.s = s
        self.process_id = process_id

class CAND:
    def __init__(self, id=0, score=0):
        self.id = id
        self.score = score

    @staticmethod
    def order(a, b):
        return a.score < b.score

def materialize1stLayerFuzzy(datafn, root, navids, index_dir, navigating_tbl):
    global fuzzy_number
    Const.logPrint("Start move data to disk file in 1st layer.")
    with open(datafn, "rb") as f:
        rest = root.size
        total = root.size
        current = 0

        while rest > 0:
            fbl = {}
            number = min(Const.fbl_series_num, rest)
            tss = np.fromfile(f, dtype=np.float32, count=number * Const.tsLength)
            tss = tss.reshape((number, Const.tsLength))

            bound = current + number

            for i in range(current, bound):
                fbl[root.ch[navids[i]]] = {"buffer": [tss[i - current]]}

            for node, recorder in navigating_tbl.items():
                if node.partition_id != -1:
                    cur_pos = recorder.actual_size
                    while cur_pos < len(recorder.series_index_list) and recorder.series_index_list[cur_pos] < bound:
                        index = recorder.series_index_list[cur_pos]
                        fbl[node] = {"buffer": [tss[index - current]]}
                        cur_pos += 1

            for node, data in fbl.items():
                outfile = os.path.join(index_dir, f"U_{node.id}" if node.partition_id == -1 else f"{node.layer}_{node.partition_id}")
                with open(outfile, "ab") as outf:
                    for buf in data["buffer"]:
                        buf.tofile(outf)

            rest -= number
            current += number
            Const.logPrint(f"Now in 1st layer {current / total * 100}% series have been written to disk.(FIRST STAGE)")

        f.seek(0)
        rest = root.size
        current = 0

        while rest > 0:
            fbl = {}
            number = min(Const.fbl_series_num, rest)
            tss = np.fromfile(f, dtype=np.float32, count=number * Const.tsLength)
            tss = tss.reshape((number, Const.tsLength))

            bound = current + number

            for node, recorder in navigating_tbl.items():
                if node.partition_id == -1:
                    cur_pos = recorder.actual_size
                    while cur_pos < len(recorder.series_index_list) and recorder.series_index_list[cur_pos] < bound:
                        index = recorder.series_index_list[cur_pos]
                        fbl[node] = {"buffer": [tss[index - current]]}
                        cur_pos += 1

            for node, data in fbl.items():
                out_file = os.path.join(index_dir, f"U_{node.id}")
                with open(out_file, "ab") as outf:
                    for buf in data["buffer"]:
                        buf.tofile(outf)

            rest -= number
            current += number
            Const.logPrint(f"Now in 1st layer {current / total * 100}% series have been written to disk.(SECOND STAGE)")

def materializeInterNodeFuzzy(node, saxes, actual_size, navigating_tbl):
    fuzzyidxfn = Const.fuzzyidxfn
    file_path = os.path.join(fuzzyidxfn, f"U_{node.id}")
    with open(file_path, "rb") as f:
        rest = node.size
        current = 0

        while rest > 0:
            lbl = {}
            number = min(Const.fbl_series_num, rest)
            tss = np.fromfile(f, dtype=np.float32, count=number * Const.tsLength)
            tss = tss.reshape((number, Const.tsLength))

            bound = current + number

            for i in range(current, bound):
                target = node.route(saxes[node.offsets[i] * Const.segmentNum:(node.offsets[i] + 1) * Const.segmentNum])
                lbl[target] = {"buffer": [tss[i - current]]}

            for n, recorder in navigating_tbl.items():
                cur_pos = recorder.actual_size
                while cur_pos < len(recorder.series_index_list) and recorder.series_index_list[cur_pos] < bound:
                    index = recorder.series_index_list[cur_pos]
                    lbl[n] = {"buffer": [tss[index - current]]}
                    cur_pos += 1

            for n, data in lbl.items():
                out_file = os.path.join(fuzzyidxfn, n.getFileName())
                with open(out_file, "ab") as outf:
                    for buf in data["buffer"]:
                        buf.tofile(outf)

            rest -= number
            current += number

class ConversionUtil:
    @staticmethod
    def extendSax(sax, bits_cardinality, chosenSegments):
        pass

    @staticmethod
    def invSaxHeadFromSax(sax, bitsCardinality, segmentNum):
        pass

    @staticmethod
    def findFirstGE(power_2, start, end, value):
        pass

    @staticmethod
    def getValueRange(sax_symbol, bits_cardinality, lb, ub):
        pass

class FileUtil:
    @staticmethod
    def checkDirClean(dir_path):
        pass

    @staticmethod
    def fileRemove(file_path):
        pass

class MathUtil:
    @staticmethod
    def generateMask(segmentNum):
        pass

def partition(nodeIn1stLayer, vertexNum):
    pass

def BuildIndexFuzzy(datafn, saxfn, paixfn, fuzzyidxfn, vertexNum, s_len, segmentNum, th, fbl_series_num, delta, logPrint):
    Const.logPrint = logPrint
    Const.segmentNum = segmentNum
    Const.tsLength = s_len
    Const.th = th
    Const.delta = delta
    Const.fbl_series_num = fbl_series_num
    Const.vertexNum = vertexNum
    Const.fuzzyidxfn = fuzzyidxfn

    nodeIn1stLayer = Node()
    nodeIn1stLayer.generateSaxAndCardIn1stLayer(range(vertexNum))

    partition(nodeIn1stLayer, vertexNum)

    sax = None
    index = None
    sa = None

    # Construct the sax and paix files
    Const.logPrint("Start construct the sax and paix files.")
    with open(saxfn, "rb") as f:
        sax = np.fromfile(f, dtype=np.int32, count=Const.vertexNum * Const.tsLength)
        sax = sax.reshape((Const.vertexNum, Const.tsLength))

    Const.logPrint("Start build index for 1st layer.")
    navigating_tbl = {}
    root = None
    navids = []

    root = Node()
    root.generateSaxAndCardIn1stLayer(range(Const.vertexNum))
    materialize1stLayerFuzzy(datafn, root, navids, fuzzyidxfn, navigating_tbl)

    Const.logPrint("Start build index for inter nodes.")
    saxes = None
    actual_size = None

    # Construct index for inter nodes
    Const.logPrint("Start construct index for inter nodes.")
    with open(saxfn, "rb") as f:
        saxes = np.fromfile(f, dtype=np.int32, count=Const.vertexNum * Const.tsLength)
        saxes = saxes.reshape((Const.vertexNum, Const.tsLength))

    actual_size = [0] * Const.vertexNum
    for node, recorder in navigating_tbl.items():
        actual_size[node.id] = recorder.actual_size

    materializeInterNodeFuzzy(node, saxes, actual_size, navigating_tbl)

    Const.logPrint("Finish build index for all nodes.")

def main():
    # Specify the parameters
    datafn = "your_data_file_name"
    saxfn = "your_sax_file_name"
    paixfn = "your_paix_file_name"
    fuzzyidxfn = "your_fuzzy_index_file_name"
    vertexNum = 10  # Specify your vertexNum
    s_len = 100  # Specify your s_len
    segmentNum = 7  # Specify your segmentNum
    th = 0.2  # Specify your th
    fbl_series_num = 2  # Specify your fbl_series_num
    delta = 1  # Specify your delta

    # Call the function to build the index
    BuildIndexFuzzy(datafn, saxfn, paixfn, fuzzyidxfn, vertexNum, s_len, segmentNum, th, fbl_series_num, delta, print)

if __name__ == "__main__":
    main()
