import os
import numpy as np

class Const:
    dataset = "deep"
    index = 1
    ops = 1
    query_num = 200
    k = 3
    visit_node_num = 25
    nprobes = 100
    dtw_window_percent = 0.05
    series_num = -1
    visited_node_num = -1
    dtw_window_size = -1

    th = 1000
    segmentNum = 16
    bitsCardinality = 8
    fbl_size = 12480
    max_diff = 2
    fuzzy_f_1st = 0.3
    fuzzy_f = 0.3
    delta = 2
    small_perc = 0.2
    max_mask_bit_percentage = 0.8
    f_low = 0.5
    f_high = 1.5
    alpha = 0.5
    tardis_sample_percent = 1

    breakpointsfn = "../breakpoints.txt"
    graphfn = "../RawGraph_16_3.bin"
    bitsReserve = 3


    pre_read = 1
    bits_reserve = -1
    paafn = "../data/paa/deep-96-100m-16.bin"
    saxfn = "../data/sax/deep-96-100m-16.bin"
    idxfn = "../data/non-mat/deep-96-100m-16.bin_le"
    memoryidxfn = "/mnt/c/Series4Similarity_Search/deep/memory/"
    fuzzyidxfn = "/mnt/c/Series4Similarity_Search/deep/fuzzy-index/"
    datafn = "/mnt/c/Series4Similarity_Search/deep/deep1b-96-100m.bin"
    queryfn = "/mnt/c/Series4Similarity_Search/deep/deep1b-96-1k.bin"
    tsLength = 96
    maxK = 500
    tsLengthPerSegment = -1
    cardinality = -1
    tsLengthBytes = -1
    vertex_num = -1
    offset = -1

    @staticmethod
    def readConfig():
        print("Reading configuration...")
        # Example configuration
        Const.dataset = "ExampleDataset"
        Const.index = 1
        Const.query_num = 10
        Const.k = 5
        Const.graphfn = "example_graph.bin"
        Const.saxfn = "example_sax.txt"
        print("Configuration read successfully.")

    @staticmethod
    def logPrint(msg):
        print(msg)  # Logging method

def main():
    Const.readConfig()

    if not os.path.exists(Const.graphfn):
        print("File not found:", Const.graphfn)
        exit(-1)

    print("Dataset:", Const.dataset)
    print("Index:", Const.index)
    print("Query Num:", Const.query_num)
    print("K:", Const.k)
    print("Graph File:", Const.graphfn)

    ts = np.array([1, 2, 3, 4, 5])
    print("Example Time Series:", ts)

    Const.logPrint("This is a log message.")

