import threading
from Utils.MathUtil import MathUtil
from Utils.FileUtil import FileUtil

def func(s, e, neighbor_number, neighbor_bits, segment_number, graph_file_name):
    graph_skeleton = [0] * neighbor_number
    with open(f"{graph_file_name}{s}", "wb") as f:
        a = MathUtil.n_choose_k(segment_number, 1)
        b = MathUtil.n_choose_k(segment_number, 2)
        c = MathUtil.n_choose_k(segment_number, 3)
        for i in range(s, e):
            if (i - s) % 10000 == 0:
                print(f"{s}: {i - s}")
            s = 0
            if neighbor_bits >= 1:
                MathUtil.get_1_bit_changed_nums(i, segment_number, graph_skeleton, s)
                s += a
            if neighbor_bits >= 2:
                MathUtil.get_2_bits_changed_nums(i, segment_number, graph_skeleton, s)
                s += b
            if neighbor_bits >= 3:
                MathUtil.get_3_bits_changed_nums(i, segment_number, graph_skeleton, s)
                s += c
            if neighbor_bits >= 4:
                MathUtil.get_4_bits_changed_nums(i, segment_number, graph_skeleton, s)
            f.write(bytes(graph_skeleton))

class GraphConstruction:
    def __init__(self, segment_number, bits_reserve):
        self.segment_number = segment_number
        self.bits_reserve = bits_reserve
        self.graph_file_name = f"../RawGraph_{self.segment_number}_{self.bits_reserve}.bin"

    def build_and_save_to_disk(self):
        arr_length = 1 << self.segment_number
        neighbor_bits = self.bits_reserve
        neighbor_num = sum(MathUtil.n_choose_k(self.segment_number, index) for index in range(1, neighbor_bits + 1))
        print("neighbor number =", neighbor_num)
        thread_number = 20
        chunk_size = arr_length // thread_number
        threads = []
        for i in range(thread_number - 1):
            threads.append(threading.Thread(target=func, args=(i * chunk_size, (i + 1) * chunk_size,
                                                               neighbor_num, neighbor_bits, self.segment_number,
                                                               self.graph_file_name)))
        threads.append(threading.Thread(target=func, args=((thread_number - 1) * chunk_size, arr_length,
                                                           neighbor_num, neighbor_bits, self.segment_number,
                                                           self.graph_file_name)))

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        sources = [f"{self.graph_file_name}{i * chunk_size}" for i in range(thread_number)]
        FileUtil.merge_files(sources, self.graph_file_name, thread_number)

if __name__ == "__main__":
    graph_construction = GraphConstruction(16, 3)
    graph_construction.build_and_save_to_disk()
