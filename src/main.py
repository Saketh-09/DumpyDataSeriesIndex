import sys
import os
import threading
import time
from typing import List

from src.Const import Const
from src.IndexConstruction.Node import Node
from src.IndexConstruction.graph_construction import GraphConstruction
from src.SearchEngine import Searcher
from Utils.FileUtil import FileUtil
from Utils.MathUtil import MathUtil
from Utils.TimeSeriesUtil import TimeSeriesUtil


def construct_graph():
    GraphConstruction.build_and_save_to_disk()


def load_graph_skeleton() -> List[List[int]]:
    vd = 0
    i = 1
    while i <= Const.bits_reserve:
        vd += MathUtil.n_choose_k(Const.segment_num, i)
        i += 1

    nn_list = [[-1] * vd for _ in range(Const.vertex_num)]

    if not os.path.exists(Const.graphfn):
        print("File not exists!", Const.graphfn)
        sys.exit(-1)

    with open(Const.graphfn, "rb") as f:
        i = 1
        while i < Const.vertex_num:
            nn_list[i] = list(map(int, f.read(vd * sizeof(int))))
            i += 1

    return nn_list


def build_dumpy():
    g = load_graph_skeleton()
    root = Node.building_index(Const.datafn, Const.saxfn)
    root.save_to_disk(Const.idxfn + "root.idx")


def approx_search_one_node():
    root = Node.load_from_disk(Const.saxfn, Const.idxfn + "root.idx", False)
    g = load_graph_skeleton()
    queries = FileUtil.read_queries()
    for i in range(Const.query_num):
        print("Query " + str(i) + ":")
        approx_knn = Searcher.approx_search(root, queries[i * Const.ts_length:(i + 1) * Const.ts_length], Const.k, g,
                                            Const.idxfn)
        print("Results:")
        for j, item in enumerate(approx_knn):
            print(j + 1, ": ", TimeSeriesUtil.time_series_to_line(item.ts))


def approx_search_more_node():
    root = Node.load_from_disk(Const.saxfn, Const.idxfn + "root.idx", False)
    queries = FileUtil.read_queries()
    for i in range(Const.query_num):
        print("Query " + str(i) + ":")
        approx_knn = Searcher.approx_inc_search(root, queries[i * Const.ts_length:(i + 1) * Const.ts_length], Const.k,
                                                Const.idxfn, Const.visited_node_num)
        print("Results:")
        for j, item in enumerate(approx_knn):
            print(j + 1, ": ", TimeSeriesUtil.time_series_to_line(item.ts))


# Define other functions in a similar way


def main():
    Const.read_config()

    if Const.index == 0:
        construct_graph()
    elif Const.index == 1:
        ops = Const.ops
        if ops == 0:
            build_dumpy()
        elif ops == 1:
            approx_search_one_node()
        elif ops == 2:
            # Call corresponding function for exact search
            pass
        elif ops == 3:
            # Call corresponding function for index stats
            pass
        elif ops == 4:
            approx_search_more_node()
        elif ops == 5:
            # Call corresponding function for DTW
            pass
        elif ops == 6:
            # Call corresponding function for DTW
            pass
        elif ops == 7:
            # Call corresponding function for ngSearch
            pass
        elif ops == 8:
            # Call corresponding function for DTW
            pass
    elif Const.index == 2:
        ops = Const.ops
        if ops == 0:
            # Call corresponding function for fuzzy
            pass
        elif ops == 1:
            # Call corresponding function for fuzzy
            pass
        elif ops == 3:
            # Call corresponding function for index stats
            pass
        elif ops == 4:
            # Call corresponding function for fuzzy
            pass
        elif ops == 7:
            # Call corresponding function for ngSearch
            pass


if __name__ == "__main__":
    main()
