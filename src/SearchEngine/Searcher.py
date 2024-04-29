class Searcher:
    @staticmethod
    def approxSearch(root: Node, query: np.ndarray, k: int, g: List[List[int]], index_dir: str) -> List[PqItemSeries]:
        queryTs = TimeSeries(query)
        heap = []

        sax = np.array(queryTs.sax)
        head = invSaxHeadFromSax(sax, Const.bitsCardinality, Const.segmentNum)
        current = root.ch[head]

        if current is None:
            node = None
            for i in range(Node.a + Node.b + Node.c):
                if root.ch[g[head][i]] is not None:
                    node = root.ch[g[head][i]]
                    break

            if node.isInternalNode():
                Searcher.approxSearchInterNode(node, queryTs, sax, k, heap, index_dir)
            else:
                node.search(k, queryTs, heap, index_dir)
        elif not current.isInternalNode():
            current.search(k, queryTs, heap, index_dir)
        else:
            Searcher.approxSearchInterNode(current, queryTs, sax, k, heap, index_dir)

        queryTs = None
        heap.sort(key=lambda x: x.dist, reverse=True)
        return heap

    @staticmethod
    def approxSearchInterNode(root: Node, queryTs: TimeSeries, sax: np.ndarray, k: int,
                              heap: List[PqItemSeries], index_dir: str):
        current = root.route(sax)
        if not current.isInternalNode():
            current.search(k, queryTs, heap, index_dir)
        else:
            minimum_distance = inf
            maximum_size = 0
            node = None

            for index, child in enumerate(current.ch):
                if child is None:
                    continue
                distance = LowerBound_Paa_iSax(queryTs.paa, child.sax, child.bits_cardinality)
                if distance < minimum_distance:
                    minimum_distance = distance
                    maximum_size = child.size
                    node = child
                elif distance == minimum_distance and child.size > maximum_size:
                    maximum_size = child.size
                    node = child

            if node.isInternalNode():
                Searcher.approxSearchInterNode(node, queryTs, sax, k, heap, index_dir)
            else:
                node.search(k, queryTs, heap, index_dir)

    @staticmethod
    def ngSearch(root: Node, query: np.ndarray, k: int, g: List[List[int]], index_dir: str) -> List[PqItemSeries]:
        queryTs = TimeSeries(query)
        heap = []

        sax = np.array(queryTs.sax)
        head = invSaxHeadFromSax(sax, Const.bitsCardinality, Const.segmentNum)
        current = root.ch[head]

        if current is None:
            node = None
            for i in range(Node.a + Node.b + Node.c):
                if root.ch[g[head][i]] is not None:
                    node = root.ch[g[head][i]]
                    break

            if node.isInternalNode():
                Searcher.approxSearchInterNode(node, queryTs, sax, k, heap, index_dir)
            else:
                node.searchDTW(k, queryTs, heap, index_dir)
        elif not current.isInternalNode():
            current.searchDTW(k, queryTs, heap, index_dir)
        else:
            Searcher.approxSearchInterNode(current, queryTs, sax, k, heap, index_dir)

        queryTs = None
        heap.sort(key=lambda x: x.dist, reverse=True)
        return heap

    @staticmethod
    def ngSearchSIMD(root: Node, query: np.ndarray, k: int, g: List[List[int]], index_dir: str,
                     upperLemire: List[float], lowerLemire: List[float]) -> List[PqItemSeries]:
        queryTs = TimeSeries(query)
        heap = []

        sax = np.array(queryTs.sax)
        head = invSaxHeadFromSax(sax, Const.bitsCardinality, Const.segmentNum)
        current = root.ch[head]

        if current is None:
            node = None
            for i in range(Node.a + Node.b + Node.c):
                if root.ch[g[head][i]] is not None:
                    node = root.ch[g[head][i]]
                    break

            if node.isInternalNode():
                Searcher.approxSearchInterNode(node, queryTs, sax, k, heap, index_dir)
            else:
                node.searchDTWSIMD(k, queryTs, heap, index_dir, upperLemire, lowerLemire)
        elif not current.isInternalNode():
            current.searchDTWSIMD(k, queryTs, heap, index_dir, upperLemire, lowerLemire)
        else:
            Searcher.approxSearchInterNode(current, queryTs, sax, k, heap, index_dir)

        queryTs = None
        heap.sort(key=lambda x: x.dist, reverse=True)
        return heap