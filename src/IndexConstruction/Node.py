import os
from typing import List

import numpy as np
import math

from src.Const import Const
from src.Utils.MathUtil import MathUtil
from collections import defaultdict
from heapq import heappop, heappush
import random

class Node:
    def __init__(self):
        self.size = 0
        self.saxes = None
        self.paas = None
        self.sax = [0] * Const.segmentNum
        self.bits_cardinality = [0] * Const.segmentNum
        self.chosenSegments = []
        self.ch = [None] * Const.vertexNum
        self.offsets = []
        self.layer = 0
        self.leaf_num = 0

    saxes = None
    paas = None
    a = MathUtil.n_choose_k(Const.segmentNum, 1)
    b = MathUtil.n_choose_k(Const.segmentNum, 2)
    c = MathUtil.n_choose_k(Const.segmentNum, 3)
    power_2 = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    combined_number = None
    combines = None

    def choose_one_segment(self, node_paa_info):
        minimum = math.inf
        min_index = -1
        for i in range(Const.segmentNum):
            large = max(node_paa_info.paa_up_size[i], node_paa_info.paa_below_size[i])
            if large < minimum:
                minimum = large
                min_index = i
        return min_index

    def choose_segment(self, node_paa_info, chosen_number):
        self.chosenSegments.clear()
        if chosen_number == 1:
            self.chosenSegments.append(self.choose_one_segment(node_paa_info))
            return

        scores = [None] * Const.segmentNum
        for i in range(Const.segmentNum):
            if self.bits_cardinality[i] >= Const.bitsCardinality:
                scores[i] = Temporary(i, -1)
            else:
                scores[i] = Temporary(i, node_paa_info.paa_variance[i])

        scores.sort(reverse=True)

        for i in range(chosen_number):
            self.chosenSegments.append(scores[i].i)
        self.chosenSegments.sort()

    def generate_sax_and_card_in_1st_layer(self, new_id):
        for i in range(Const.segmentNum - 1, -1, -1):
            t = new_id % 2
            new_id >>= 1
            self.sax[i] = t

    def generate_sax_and_cardinality(self, node, new_id):
        node.sax = self.sax[:]
        node.bits_cardinality = self.bits_cardinality[:]
        for i in range(len(self.chosenSegments) - 1, -1, -1):
            seg = self.chosenSegments[i]
            node.bits_cardinality[seg] += 1
            t = new_id % 2
            new_id >>= 1
            node.sax[seg] = (node.sax[seg] << 1) + t

    def generate_sax_and_card_in_1st_layer_4_leaf_node(self, new_id):
        if self.bits_cardinality[0] == -1:
            self.bits_cardinality = [1] * Const.segmentNum
            self.generate_sax_and_card_in_1st_layer(new_id)
            return
        for i in range(Const.segmentNum - 1, -1, -1):
            t = new_id % 2
            new_id >>= 1
            if self.bits_cardinality[i] == 1 and self.sax[i] != t:
                self.bits_cardinality[i] = 0

    def generate_sax_and_cardinality_4_leaf_node(self, node, new_id):
        if node.bits_cardinality[0] == -1:
            self.generate_sax_and_cardinality(node, new_id)
            return
        for i in range(len(self.chosenSegments) - 1, -1, -1):
            seg = self.chosenSegments[i]
            t = new_id % 2
            new_id >>= 1
            if (node.bits_cardinality[seg] == self.bits_cardinality[seg] + 1
                    and node.sax[seg] % 2 != t):
                node.bits_cardinality[seg] -= 1
                node.sax[seg] >>= 1

    @staticmethod
    def partition(nodes_map: List[partUnit], chosen_segment_number: int) -> int:
        nodes = []
        input_node_number = 1 << chosen_segment_number
        total_size = 0
        for i in range(input_node_number):
            if Const.th * Const.small_perc < nodes_map[i].s < Const.th * Const.small_perc and nodes_map[i].s > 0:
                nodes.append(nodes_map[i])
                total_size += nodes_map[i].s

        node_number = len(nodes)
        if node_number <= 3:
            return 0

        pid = 0
        first_round_pack_num = int(total_size / Const.th)
        if node_number <= first_round_pack_num:
            return 0

        max_mask_num = int(chosen_segment_number * Const.max_mask_bit_percentage)
        packs = [pack(nodes[i], chosen_segment_number, i) for i in range(first_round_pack_num)]
        packs.sort(key=attrgetter('size'))

        for cur_node in nodes:
            if cur_node.process_id != -1:
                continue
            cur_id = cur_node.id
            minimum_cost = chosen_segment_number
            pack_id = -1
            for p in packs:
                if p.tot_size + cur_node.s > Const.th or p.masked_bits_num >= max_mask_num:
                    continue
                cost = p.calc_cost(cur_id, chosen_segment_number)
                if cost + p.masked_bits_num >= max_mask_num:
                    continue
                if cost < minimum_cost:
                    minimum_cost = cost
                    pack_id = p.process_id

            if pack_id == -1:
                packs.append(pack(cur_node, chosen_segment_number, len(packs)))
            else:
                packs[pack_id].insert(cur_node, chosen_segment_number)

        # merge packs
        process_id_map = {i: i for i in range(len(packs))}
        packs.sort(key=attrgetter('size'))
        for i in range(len(packs)):
            cur_pack = packs[i]
            if cur_pack.process_id != process_id_map[cur_pack.process_id]:
                continue
            minimum_cost = chosen_segment_number
            minimum_size = math.inf
            min_pack_id = -1
            for j in range(len(packs)):
                target_pack = packs[j]
                if (target_pack.disabled or cur_pack.process_id == target_pack.process_id
                        or cur_pack.tot_size + target_pack.tot_size > Const.th
                        or cur_pack.masked_bits_num >= max_mask_num
                        or target_pack.masked_bits_num >= max_mask_num):
                    continue
                cur_cost, tar_cost = cur_pack.calc_pack_merge_cost(target_pack, chosen_segment_number)
                if (cur_cost + cur_pack.masked_bits_num >= max_mask_num
                        or tar_cost + target_pack.masked_bits_num >= max_mask_num):
                    continue
                cost = cur_cost + tar_cost
                if cost < minimum_cost or (cost == minimum_cost and cur_pack.tot_size < minimum_size):
                    minimum_cost = cost
                    minimum_size = target_pack.tot_size
                    min_pack_id = j

            if minimum_size < math.inf:
                cur_pack.merge_pack(packs[min_pack_id], chosen_segment_number)

        # re-assign the process ids to the nodes
        maximum_pid = 0
        for node in nodes:
            node.process_id = process_id_map[node.process_id]
            maximum_pid = max(maximum_pid, node.process_id)

        return maximum_pid + 1

    def get_index_stats(self):
        total_leaf_node_num = self.get_leaf_node_num()
        total_size = self.get_total_size()
        print("Total size =", total_size)
        print("Total nodes number =", self.get_node_num())
        print("Leaf node number =", total_leaf_node_num)
        print("1st layer node number =", self.get_1st_layer_nodes_no())
        print("1st layer internal node number =", self.get_1st_layer_inter_nodes_no())
        print("1st layer internal series number =", self.get_1st_layer_inter_node_series_no())
        print("Max. height =", self.get_max_height() - 1)
        print("Avg. Height =", self.get_sum_height() / total_leaf_node_num)
        print("Avg. Filling Factor =", total_size / total_leaf_node_num / Const.th)
        print("Bias leaf node ratio =", self.get_bias_leaf_node_num() / total_leaf_node_num)

    def get_1st_layer_inter_nodes_no(self):
        nodes = set()
        for child in self.ch:
            if child and child.size > Const.th and child not in nodes:
                nodes.add(child)
        return len(nodes)

    def get_1st_layer_nodes_no(self):
        nodes = set()
        for child in self.ch:
            if child and child not in nodes:
                nodes.add(child)
        return len(nodes)

    def get_1st_layer_inter_node_series_no(self):
        nodes = set()
        ret = 0
        for child in self.ch:
            if (child and child.size > Const.th and child not in nodes):
                nodes.add(child)
                ret += child.size
        return ret

    def get_max_height(self):
        if self.is_leaf_node():
            return 1
        max_height = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                max_height = max(child.get_max_height(), max_height)
                hash_map.add(child)
        return max_height + 1

    def get_leaf_node_num(self):
        if self.is_leaf_node():
            return 1
        sum = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                sum += child.get_leaf_node_num()
                hash_map.add(child)
        return sum

    def get_bias_leaf_node_num(self):
        if self.is_leaf_node():
            max_b = max(self.bits_cardinality)
            min_b = min(self.bits_cardinality)
            return max_b - min_b >= 4
        sum = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                sum += child.get_bias_leaf_node_num()
                hash_map.add(child)
        return sum

    def get_total_size(self):
        if self.is_leaf_node():
            return self.size
        sum = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                sum += child.get_total_size()
                hash_map.add(child)
        return sum

    def get_node_num(self):
        if self.is_leaf_node():
            return 1
        sum = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                sum += child.get_node_num()
                hash_map.add(child)
        return sum + 1

    def get_sum_height(self):
        if self.is_leaf_node():
            return self.layer
        sum_height = 0
        hash_map = set()
        for child in self.ch:
            if child and child not in hash_map:
                sum_height += child.get_sum_height()
                hash_map.add(child)
        return sum_height

    def load_sax(self, saxfn):
        f_size = os.path.getsize(saxfn)
        series_num = f_size // (2 * Const.segmentNum)
        self.saxes = np.fromfile(saxfn, dtype=np.uint16, count=series_num * Const.segmentNum)
        print("Finish loading sax")
        return series_num

    def load_paa(self, paafn):
        f_size = os.path.getsize(paafn)
        self.paas = np.fromfile(paafn, dtype=np.float32)
        print("Finish loading paa")

    def generate_sax_and_paa_tbl(self):
        fn = Const.datafn
        fs = os.path.getsize(fn)
        series_num = fs // Const.tsLengthBytes
        print("Total Series Number is", series_num)
        self.saxes = np.zeros(series_num * Const.segmentNum, dtype=np.uint16)
        self.paas = np.zeros(series_num * Const.segmentNum, dtype=np.float32)
        rest = series_num
        cur = 0
        with open(fn, "rb") as f:
            while rest > 0:
                num = min(rest, 4000000)
                tss = np.fromfile(f, dtype=np.float32, count=num * Const.tsLength)
                for i in range(num):
                    if np.isnan(tss[i * Const.tsLength]):
                        self.saxes[i * Const.segmentNum: (i + 1) * Const.segmentNum] = 0
                        print("Dirty data:", i)
                    else:
                        paas_tmp, saxes_tmp = ConversionUtil.paa_and_sax_from_ts(
                            tss[i * Const.tsLength: (i + 1) * Const.tsLength],
                            Const.tsLengthPerSegment, Const.segmentNum, Const.cardinality
                        )
                        self.paas[cur + i] = paas_tmp
                        self.saxes[cur + i] = saxes_tmp
                rest -= num
                cur += num
        return series_num

    def generate_sax_tbl(self):
        fn = Const.datafn
        if Const.series_num == -1:
            fs = os.path.getsize(fn)
            series_num = fs // Const.tsLengthBytes
        else:
            series_num = Const.series_num
        print("Total Series Number is", series_num)
        self.saxes = np.zeros(series_num * Const.segmentNum, dtype=np.uint16)
        rest = series_num
        cur = 0
        with open(fn, "rb") as f:
            while rest > 0:
                num = min(rest, 2000000)
                tss = np.fromfile(f, dtype=np.float32, count=num * Const.tsLength)
                for i in range(num):
                    if np.isnan(tss[i * Const.tsLength]):
                        self.saxes[i * Const.segmentNum: (i + 1) * Const.segmentNum] = 0
                        print("Dirty data:", i)
                    else:
                        self.saxes[cur + i] = ConversionUtil.sax_from_ts(
                            tss[i * Const.tsLength: (i + 1) * Const.tsLength],
                            Const.tsLengthPerSegment, Const.segmentNum, Const.cardinality
                        )
                rest -= num
                cur += num
        return series_num

    @staticmethod
    def load_from_disk(saxfn, idxfn, need_sax):
        node = Node()
        if need_sax:
            node.load_sax(saxfn)
        with open(idxfn, "rb") as f:
            node = pickle.load(f)
        return node

    def assign_leaf_num(self):
        if not self.is_internal_node():
            self.leaf_num = 1
            return 1

        visited = set()
        for child in self.ch:
            if child and child not in visited:
                visited.add(child)
                self.leaf_num += child.assign_leaf_num()

        return self.leaf_num
    @staticmethod
    def loadCombines():
        base = "../combines/" + str(Const.segmentNum) + "-"
        ret = [None] * (Const.segmentNum + 1)
        combined_number = [0] * Const.segmentNum
        with open("../combines/cnum-" + str(Const.segmentNum) + ".txt", "r") as ff:
            for i in range(Const.segmentNum):
                combined_number[i] = int(ff.readline())

        for i in range(1, Const.segmentNum):
            ret[i] = []
            with open(base + str(i) + ".txt", "r") as f:
                for j in range(combined_number[i]):
                    ret[i].append([int(x) for x in f.readline().split()])

        ret[Const.segmentNum] = [[i for i in range(Const.segmentNum)]]
        Node.combined_number = combined_number
        Node.combines = ret

    def route(self, _sax):
        if self.isLeafNode():
            return self
        nav_id = ConversionUtil.extendSax(_sax, self.bits_cardinality, self.chosenSegments)
        if self.ch[nav_id] is None:
            return self
        return self.ch[nav_id].route(_sax)

    def search(self, k, queryTs, heap, index_dir):
        bsf = np.inf if len(heap) < k else heap[0].dist
        fn = index_dir + self.getFileName()

        with open(fn, "rb") as f:
            ts = np.fromfile(f, dtype=np.float32, count=self.size * Const.tsLength)
            ts = ts.reshape(self.size, Const.tsLength)

            for i in range(self.size):
                dist = TimeSeriesUtil.euclideanDist(queryTs.ts, ts[i], Const.tsLength, bsf)

                if len(heap) < k:
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])
                elif dist < bsf:
                    heapq.heappop(heap)
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])

                if len(heap) >= k:
                    bsf = heap[0].dist

        for s in heap:
            if s.needDeepCopy:
                s.copyData()

    def searchWithHashSet(self, k, queryTs, heap, index_dir, hash_set):
        assert self.isLeafNode()
        bsf = np.inf if len(heap) < k else heap[0].dist
        fn = index_dir + self.getFileName()

        with open(fn, "rb") as f:
            ts = np.fromfile(f, dtype=np.float32, count=self.size * Const.tsLength)
            ts = ts.reshape(self.size, Const.tsLength)

            for i in range(self.size):
                if ts[i].data in hash_set:
                    continue
                dist = TimeSeriesUtil.euclideanDist(queryTs.ts, ts[i], Const.tsLength, bsf)

                if len(heap) < k:
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])
                    hash_set.add(ts[i].data)
                elif dist < bsf:
                    heapq.heappop(heap)
                    hash_set.remove(heap[-1].ts.data)
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])
                    hash_set.add(ts[i].data)

                if len(heap) >= k:
                    bsf = heap[0].dist

        for s in heap:
            if s.needDeepCopy:
                s.copyData()

    def searchDTW(self, k, queryTs, heap, index_dir):
        assert not self.isInternalNode()
        bsf = np.inf if len(heap) < k else heap[0].dist
        fn = index_dir + self.getFileName()

        with open(fn, "rb") as f:
            ts = np.fromfile(f, dtype=np.float32, count=self.size * Const.tsLength)
            ts = ts.reshape(self.size, Const.tsLength)

            for i in range(self.size):
                dist = TimeSeriesUtil.dtw(queryTs.ts, ts[i], Const.tsLength, Const.dtw_window_size, bsf)

                if len(heap) < k:
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])
                elif dist < bsf:
                    heapq.heappop(heap)
                    heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                    heapq.heappush(heap, heap[-1])

                if len(heap) >= k:
                    bsf = heap[0].dist

        for s in heap:
            if s.needDeepCopy:
                s.copyData()

    def searchDTWSIMD(self, k, queryTs, heap, index_dir, upperLemire, lowerLemire):
        assert not self.isInternalNode()
        bsf = np.inf if len(heap) < k else heap[0].dist
        fn = index_dir + self.getFileName()

        with open(fn, "rb") as f:
            ts = np.fromfile(f, dtype=np.float32, count=self.size * Const.tsLength)
            ts = ts.reshape(self.size, Const.tsLength)

            cb = np.zeros(Const.tsLength, dtype=np.float32)
            cb1 = np.zeros(Const.tsLength, dtype=np.float32)
            length = 2 * Const.dtw_window_size + 1
            tSum = np.zeros(length, dtype=np.float32)
            pCost = np.zeros(length, dtype=np.float32)
            rDist = np.zeros(length, dtype=np.float32)

            for i in range(self.size):
                dist2 = TimeSeriesUtil.lb_keogh_data_bound(ts[i], upperLemire, lowerLemire, cb1, Const.tsLength, bsf)
                if dist2 < bsf:
                    cb[Const.tsLength - 1] = cb1[Const.tsLength - 1]
                    for ii in range(Const.tsLength - 2, -1, -1):
                        cb[ii] = cb[ii + 1] + cb1[ii]

                    dist = TimeSeriesUtil.dtwsimd(queryTs.ts, ts[i], cb, Const.tsLength, Const.dtw_window_size, bsf, tSum, pCost, rDist)

                    if len(heap) < k:
                        heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                        heapq.heappush(heap, heap[-1])
                    elif dist < bsf:
                        heapq.heappop(heap)
                        heap.append(PqItemSeries(np.copy(ts[i]), dist, False, True))
                        heapq.heappush(heap, heap[-1])

                    if len(heap) >= k:
                        bsf = heap[0].dist

        for s in heap:
            if s.needDeepCopy:
                s.copyData()

    def BuildingIndex(self, datafn, saxfn):
        print("Start building index.")
        # FileUtil.checkDirClean(Const.idxfn.c_str()) # FileUtil implementation is missing
        self.loadCombines()
        series_number = self.generateSaxTbl()

        root = Node()
        root.size = series_number
        root.bits_cardinality = [0] * Const.vertexNum
        nodeIn1stLayer = [partUnit() for _ in range(Const.vertexNum)]
        navids = [0] * series_number

        # Obtain initial layer node size
        for index in range(series_number):
            asax = self.saxes[index * Const.segmentNum: (index + 1) * Const.segmentNum]
            nav_id = ConversionUtil.invSaxHeadFromSax(asax, Const.bitsCardinality, Const.segmentNum)
            navids[index] = nav_id
            nodeIn1stLayer[nav_id].s += 1

        print("Finish statistic size of nodes in the 1st layer.")

        # Partition 1st layer
        partNum = self.partition(nodeIn1stLayer, Const.segmentNum)
        print("Finish partition")
        childrenList = [Node(1, index) for index in range(partNum)]
        root.ch = [None] * Const.vertexNum
        for i in range(Const.vertexNum):
            if nodeIn1stLayer[i].s <= 0:
                continue
            if nodeIn1stLayer[i].s > Const.th:
                root.ch[i] = Node(1, nodeIn1stLayer[i].s, i)
                root.ch[i].generateSaxAndCardIn1stLayer(i)
            elif nodeIn1stLayer[i].process_id == -1:
                root.ch[i] = Node(1, nodeIn1stLayer[i].s, i)
                root.ch[i].generateSaxAndCardIn1stLayer(i)
            else:
                pid = nodeIn1stLayer[i].process_id
                root.ch[i] = childrenList[pid]
                childrenList[pid].size += nodeIn1stLayer[i].s
                childrenList[pid].generateSaxAndCardIn1stLayer4LeafNode(i)

        print("Finish build index structure 1st layer.")

        # Add data offsets to internal nodes in first layer
        for i in range(Const.vertexNum):
            if nodeIn1stLayer[i].s > Const.th:
                root.ch[i].offsets = []
        for i in range(series_number):
            nav_id = navids[i]
            root.ch[nav_id].offsets.append(i)

        print("Data offsets have been put into nodes in the 1st layer.")

        j = 0
        milestone = 0.1 * Const.vertexNum
        print("Start grow the index structure")
        for i in range(Const.vertexNum):
            if nodeIn1stLayer[i].s > Const.th:
                root.ch[i].growIndex()
            if (j + 1) % milestone == 0:
                print(f"{j + 1} nodes in the 1st layer has been processed.")
            j += 1

        print("Build index skeleton finished.")

        print("Start materialize leaves")
        # materializeAllLeavesWithSax(datafn, root, navids, Const.idxfn, saxes) # implementation is missing
        print("Build index successfully!")
        del self.saxes

        return root

    def growIndex(self):
        if self.size <= Const.th:
            return
        self.determineSegments()
        chosen_num = len(self.chosenSegments)

        # Statistic children information in order to partition
        nodes = [partUnit() for _ in range(1 << chosen_num)]
        node_offsets = [[] for _ in range(1 << chosen_num)]

        for i in range(self.size):
            new_id = ConversionUtil.extendSax(self.saxes[i * Const.segmentNum: (i + 1) * Const.segmentNum],
                                              self.bits_cardinality, self.chosenSegments)
            nodes[new_id].s += 1
            node_offsets[new_id].append(self.offsets[i])

        if self.layer > 1:
            self.offsets = []

        partNum = self.partition(nodes, chosen_num)

        childrenList = [Node(self, i) for i in range(partNum)]
        self.ch = [None] * (1 << chosen_num)
        for i in range(1 << chosen_num):
            if nodes[i].s <= 0:
                continue
            elif nodes[i].s > Const.th:
                self.ch[i] = Node(self, nodes[i].s, i)
                self.generateSaxAndCardinality(self.ch[i], i)
                self.ch[i].offsets = node_offsets[i].copy()
            elif self.partition_id == -1:
                self.ch[i] = Node(self, nodes[i].s, i)
                self.generateSaxAndCardinality(self.ch[i], i)
            else:
                _pid = nodes[i].process_id
                self.ch[i] = childrenList[_pid]
                childrenList[_pid].size += nodes[i].s
                self.generateSaxAndCardinality4LeafNode(self.ch[i], i)

        for child in self.ch:
            if child is not None and child.size > Const.th:
                child.growIndex()

    def determineFanout(self, lambda_minimum, lambda_max):
        if self.size < 2 * Const.th:
            lambda_minimum[0] = 1
            lambda_max[0] = 1
            return
        lambda_minimum[0] = -1
        lambda_max[0] = -1
        _min = self.size / (Const.th * Const.f_high)
        _max = self.size / (Const.th * Const.f_low)
        for i in range(1, Const.segmentNum + 1):
            if lambda_minimum[0] == -1:
                if (1 << i) >= _min:
                    lambda_minimum[0] = i
            else:
                if (1 << i) == _max:
                    lambda_max[0] = i
                    break
                elif (1 << i) > _max:
                    lambda_max[0] = max(i - 1, lambda_minimum[0])
                    break

    def determineSegments(self):
        lambda_min = [0]
        lambda_max = [0]
        self.determineFanout(lambda_min, lambda_max)

        unit_size = [0] * Const.vertexNum

        data_seg_symbols = [defaultdict(int) for _ in range(Const.segmentNum)]
        for offset in self.offsets:
            cur_sax = self.saxes[offset * Const.segmentNum: (offset + 1) * Const.segmentNum]
            for i in range(Const.segmentNum):
                data_seg_symbols[i][cur_sax[i]] += 1
            head = ConversionUtil.extendSax(cur_sax, self.bits_cardinality)
            unit_size[head] += 1

        data_seg_mean = [0] * Const.segmentNum
        data_seg_stdev = [0] * Const.segmentNum
        for i in range(Const.segmentNum):
            map = data_seg_symbols[i]
            for symbol, count in map.items():
                mid_value = ConversionUtil.getMidLineFromSaxSymbolbc8(symbol)
                data_seg_mean[i] += mid_value * count
            data_seg_mean[i] /= self.size
            for symbol, count in map.items():
                mid_value = ConversionUtil.getMidLineFromSaxSymbolbc8(symbol)
                data_seg_stdev[i] += (count * ((mid_value - data_seg_mean[i]) ** 2))
            data_seg_stdev[i] /= self.size

        plan_num = 0
        if lambda_max[0] < Const.segmentNum:
            plan_num = self.combined_number[lambda_max[0]]
        else:
            plan_num = 1
        visited = set()
        max_score = 0
        best_plan = []
        for i in range(plan_num):
            plan = self.combines[lambda_max[0]][i]
            plan_node_sizes = [0] * (1 << lambda_max[0])
            mask_code = MathUtil.generateMaskSettingKbits(plan, lambda_max[0], Const.segmentNum)
            max_node_size = defaultdict(int)
            for j in range(Const.vertexNum):
                max_node_size[mask_code & j] += unit_size[j]

            for k, v in max_node_size.items():
                plan_node_sizes[k] = v

            score = self.compute_score(plan_node_sizes, plan, lambda_max[0], data_seg_stdev)
            if score > max_score:
                max_score = score
                best_plan = plan[:]

            if lambda_min[0] <= lambda_max[0] - 1:
                self.visitPlanFromBaseTable(visited, lambda_max[0] - 1, plan, plan_node_sizes, max_score, best_plan,
                                            lambda_min[0], mask_code, data_seg_stdev, score)

    def compute_score(self, node_sizes, plan, lambda_val, data_seg_stdev):
        if self.size < 2 * Const.th:
            if node_sizes[0] > Const.th or node_sizes[1] > Const.th:
                return min(node_sizes[0], node_sizes[1]) / Const.th
            return data_seg_stdev[plan[0]] * 100
        over_th_nodes_no = sum(1 for _ in node_sizes if _ > Const.th)
        w = over_th_nodes_no / len(node_sizes)
        sum_seg = sum(data_seg_stdev[plan[i]] for i in range(lambda_val)) / lambda_val
        sum_seg = math.sqrt(sum_seg)
        sum_seg = math.exp(1 + sum_seg)

        tmp = [node_size / Const.th for node_size in node_sizes]
        stdev_fill_factor = np.std(tmp)

        balance = math.exp(-(1 + w) * stdev_fill_factor)
        ret = sum_seg + Const.alpha * balance
        return ret

    def visitPlanFromBaseTable(self, visited, cur_lambda, plan, base_tbl, max_score, best_plan, lambda_min, mask_code,
                               data_seg_stdev, base_score):
        base_mask = (1 << cur_lambda + 1) - 1

        for i in range(cur_lambda + 1):
            reset_pos = plan[i]
            cur_whole_mask = mask_code - (1 << (Const.segmentNum - 1 - reset_pos))
            if cur_whole_mask in visited:
                continue
            visited.add(cur_whole_mask)

            new_plan = [p for j, p in enumerate(plan) if j != i]
            cur_base_mask = base_mask - (1 << (cur_lambda - i))
            node_size_map = defaultdict(int)
            for j in range(len(base_tbl)):
                node_size_map[cur_base_mask & j] += base_tbl[j]

            new_tbl = [0] * (1 << cur_lambda)
            for k, v in node_size_map.items():
                new_tbl[k] = v

            score = self.compute_score(new_tbl, new_plan, cur_lambda, data_seg_stdev)
            if score > max_score[0]:
                max_score[0] = score
                best_plan.clear()
                best_plan.extend(new_plan)

            if cur_lambda > lambda_min:
                self.visitPlanFromBaseTable(visited, cur_lambda - 1, new_plan, new_tbl, max_score, best_plan,
                                            lambda_min, cur_whole_mask, data_seg_stdev, score)

    def statPaa(self):
        r = PAA_INFO()
        split_line = [0] * Const.segmentNum
        paa_max = [-float('inf')] * Const.segmentNum
        paa_min = [float('inf')] * Const.segmentNum
        paa_mu = [0] * Const.segmentNum

        for i in range(Const.segmentNum):
            lb, split_line[i] = ConversionUtil.getValueRange(self.sax[i] << 1, self.bits_cardinality[i] + 1)
        paa_up_size = [0] * Const.segmentNum
        paa_below_size = [0] * Const.segmentNum
        paa_variance = [0] * Const.segmentNum

        for offset in self.offsets:
            start = self.paas[offset * Const.segmentNum: (offset + 1) * Const.segmentNum]
            for i in range(Const.segmentNum):
                value = start[i]
                paa_mu[i] += value
                paa_min[i] = min(paa_min[i], value)
                paa_max[i] = max(paa_max[i], value)
                if value > split_line[i]:
                    paa_up_size[i] += 1
                else:
                    paa_below_size[i] += 1

        for index in range(len(paa_mu)):
            paa_mu[index] /= self.size

        for offset in self.offsets:
            start = self.paas[offset * Const.segmentNum: (offset + 1) * Const.segmentNum]
            for i in range(Const.segmentNum):
                value = start[i]
                paa_variance[i] += (value - paa_mu[i]) ** 2

        r.paa_up_size = paa_up_size
        r.paa_below_size = paa_below_size
        r.paa_variance = paa_variance
        return r

class pack:
    def __init__(self):
        self.current_bits = []
        self.current_mask = []
        self.tot_size = 0
        self.process_id = 0
        self.masked_bits_num = 0
        self.disabled = False

    def __init__(self, node, chosen_segment_number, _pid):
        self.tot_size = node.s
        self.masked_bits_num = 0
        self.current_bits = [False] * chosen_segment_number
        self.current_mask = [False] * chosen_segment_number
        _id = node.id
        for i in range(chosen_segment_number):
            self.current_bits[chosen_segment_number - 1 - i] = _id % 2
            _id >>= 1
        self.process_id = _pid
        node.process_id = self.process_id
        self.disabled = False

    def calc_cost(self, _id, chosen_seg_num):
        cost = 0
        for index in range(chosen_seg_num):
            if not self.current_mask[chosen_seg_num - 1 - index] and self.current_bits[chosen_seg_num - 1 - index] != (_id % 2):
                cost += 1
            _id >>= 1
        return cost

    def calc_pack_merge_cost(self, p, chosen_seg_num):
        cur_cost = 0
        tar_cost = 0
        cost = 0
        for i in range(chosen_seg_num):
            if self.current_mask[i] and p.current_mask[i]:
                continue
            if self.current_mask[i] and not p.current_mask[i]:
                tar_cost += 1
                cost += 1
            elif not self.current_mask[i] and p.current_mask[i]:
                cur_cost += 1
                cost += 1
            elif self.current_bits[i] != p.current_bits[i]:
                cur_cost += 1
                tar_cost += 1
                cost += 1
        return cost

    def merge_pack(self, p, chosen_seg_num):
        dis_one, res_one = (p, self) if self.process_id < p.process_id else (self, p)
        dis_one.disabled = True
        res_one.tot_size += p.tot_size
        for i in range(chosen_seg_num):
            if dis_one.current_mask[i] and res_one.current_mask[i]:
                continue
            if dis_one.current_mask[i] and not res_one.current_mask[i]:
                res_one.current_mask[i] = True
                res_one.masked_bits_num += 1
            elif not dis_one.current_mask[i] and res_one.current_mask[i]:
                continue
            elif dis_one.current_bits[i] != res_one.current_bits[i]:
                res_one.masked_bits_num += 1
                res_one.current_mask[i] = True

    def insert(self, node, chosen_seg_num):
        node.process_id = self.process_id
        _id = node.id
        self.tot_size += node.s
        for i in range(chosen_seg_num):
            if not self.current_mask[chosen_seg_num - 1 - i] and self.current_bits[chosen_seg_num - 1 - i] != (_id % 2):
                self.current_mask[chosen_seg_num - 1 - i] = True
                self.masked_bits_num += 1
            _id >>= 1

    @staticmethod
    def compare_size(x, y):
        return x.tot_size < y.tot_size

def materializeAllLeavesWithSax(datafn, root, navids, index_dir, sax_tbl):
    start_sax = time.time()
    print("Start move sax to disk file in 1st layer.")

    sax_buffer = defaultdict(list)
    for i in range(root.size):
        sax = sax_tbl[i * Const.segmentNum: (i + 1) * Const.segmentNum]
        node = root.route(sax)
        sax_buffer[node].append(sax)

    for node, buffer in sax_buffer.items():
        outfile = Const.idxfn + node.getFileName() + "_sax"
        if node.partition_id == -1:
            outfile += "_L"
        with open(outfile, "a") as outf:
            for sax in buffer:
                outf.write(" ".join(map(str, sax)) + "\n")

    sax_buffer.clear()

    start_t = time.time()
    print("Start move data to disk file in 1st layer.")
    with open(datafn, "rb") as f:
        rest = root.size
        total = root.size
        current = 0
        lbl = defaultdict(LBL_UNIT)

        while rest > 0:
            lbl.clear()
            n = min(Const.fbl_series_num, rest)
            tss = np.fromfile(f, dtype=np.float32, count=n * Const.tsLength)
            tss = tss.reshape(n, Const.tsLength)

            for index in range(current, current + n):
                sax = sax_tbl[index * Const.segmentNum: (index + 1) * Const.segmentNum]
                node = root.route(sax)
                lbl[node].buffer.append(tss[index - current])

            for node, lbl_unit in lbl.items():
                out_file = Const.idxfn + node.getFileName()
                if node.partition_id == -1:
                    out_file += "_L"
                with open(out_file, "ab") as outf:
                    lbl_unit.buffer.tofile(outf)

            rest -= n
            current += n
            print("Now Materialize all leaves. Progress:", (current / total) * 100, "%")

    print("Done.")


class PqItemSeries:
    def __init__(self, ts, dist, needDeepCopy, isLeafNode):
        self.ts = ts
        self.dist = dist
        self.needDeepCopy = needDeepCopy
        self.isLeafNode = isLeafNode

    def copyData(self):
        # Your implementation to copy data goes here
        pass


class TimeSeries:
    pass


class Const:
    tsLength = None
    dtw_window_size = None
    fbl_series_num = None

    @staticmethod
    def timerStart(io):
        pass


class TimeSeriesUtil:
    @staticmethod
    def euclideanDist(ts1, ts2, length, bsf):
        pass

    @staticmethod
    def dtw(ts1, ts2, length, dtw_window_size, bsf):
        pass

    @staticmethod
    def dtwsimd(ts1, ts2, cb, length, dtw_window_size, bsf, tSum, pCost, rDist):
        pass

    @staticmethod
    def lb_keogh_data_bound(ts, upperLemire, lowerLemire, cb1, length, bsf):
        pass

class Const:
    vertexNum = 100
    th = 50
    segmentNum = 10
    idxfn = ""  # Add path for idxfn
    f_high = 0.3
    f_low = 0.03


class partUnit:
    def __init__(self):
        self.id = 0
        self.s = 0
        self.process_id = -1


class ConversionUtil:
    @staticmethod
    def invSaxHeadFromSax(asax, bitsCardinality, segmentNum):
        pass

    @staticmethod
    def extendSax(sax, bits_cardinality, chosenSegments):
        pass

    @staticmethod
    def getMidLineFromSaxSymbolbc8(symbol):
        pass

    @staticmethod
    def getValueRange(value, cardinality):
        pass


class MathUtil:
    @staticmethod
    def generateMaskSettingKbits(plan, lambda_max, segmentNum):
        pass

    @staticmethod
    def deviation(tmp, size):
        pass


class PAA_INFO:
    def __init__(self):
        self.paa_up_size = []
        self.paa_below_size = []
        self.paa_variance = []


class Const:
    segmentNum = 64
    bitsCardinality = 4
    vertexNum = 2
    th = 128
    small_perc = 0.25
    max_mask_bit_percentage = 0.1
    tsLengthBytes = 1024
    tsLength = 256
    tsLengthPerSegment = 4
    cardinality = 8
    datafn = "data.bin"
    series_num = -1


class partUnit:
    def __init__(self):
        self.s = 0
        self.id = 0
        self.process_id = -1


class Temporary:
    def __init__(self, i, v):
        self.i = i
        self.v = v

    def __lt__(self, other):
        return self.v < other.v


class pack:
    def __init__(self, node, chosen_segment_number, process_id):
        self.tot_size = node.s
        self.masked_bits_num = 0
        self.size = node.s
        self.disabled = False
        self.process_id = process_id
        self.members = [node]
        self.offsets = []
        self.merge_count = defaultdict(int)
        self.calc_offsets(chosen_segment_number)

    def insert(self, node, chosen_segment_number):
        self.tot_size += node.s
        self.size += node.s
        self.members.append(node)
        self.calc_offsets(chosen_segment_number)

    def calc_offsets(self, chosen_segment_number):
        self.offsets.clear()
        for member in self.members:
            self.offsets.append(member.id)
            for seg in range(Const.segmentNum):
                if member.sax[seg] == 0:
                    self.merge_count[(seg, 0)] += 1
                else:
                    self.merge_count[(seg, 1)] += 1

        self.masked_bits_num = sum(v for k, v in self.merge_count.items() if v > 1)
        if self.masked_bits_num >= int(chosen_segment_number * Const.max_mask_bit_percentage):
            self.disabled = True
            return

        while self.size > Const.th:
            seg, bit = max(self.merge_count.items(), key=lambda x: x[1])
            self.size -= 1
            self.merge_count[(seg[0], 0)] -= 1
            self.merge_count[(seg[0], 1)] -= 1
            self.masked_bits_num -= 1
            self.offsets.pop()
            if self.masked_bits_num < int(chosen_segment_number * Const.max_mask_bit_percentage):
                self.disabled = False
                return

    def calc_cost(self, cur_id, chosen_segment_number):
        cost = 0
        for seg in range(Const.segmentNum):
            if self.merge_count[(seg, 0)] + self.merge_count[(seg, 1)] >= 1:
                cost += 1
        return cost

    def calc_pack_merge_cost(self, target_pack, chosen_segment_number):
        cur_cost = self.calc_cost(0, chosen_segment_number)
        tar_cost = target_pack.calc_cost(0, chosen_segment_number)
        return cur_cost, tar_cost

    def merge_pack(self, target_pack, chosen_segment_number):
        for member in target_pack.members:
            self.insert(member, chosen_segment_number)


class ConversionUtil:
    @staticmethod
    def paa_and_sax_from_ts(ts, paas, saxes, ts_length_per_segment, segment_num, cardinality):
        # Your PAA and SAX conversion logic here
        pass

    @staticmethod
    def sax_from_ts(ts, saxes, ts_length_per_segment, segment_num, cardinality):
        # Your SAX conversion logic here
        pass


def main():
    # Example usage
    node = Node.load_from_disk("sax.bin", "index.bin", True)
    node.get_index_stats()


if __name__ == "__main__":
    main()
# Implement FileUtil, PAA_INFO, and other missing classes and methods
