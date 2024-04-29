import math

class TimeSeriesUtil:
    @staticmethod
    def deviation(ts, start, end, dAve):
        length = end - start
        dVar = sum((ts[i] - dAve) ** 2 for i in range(start, end))
        return math.sqrt(dVar / length)

    @staticmethod
    def avg(time_series, start, end):
        return sum(time_series[start:end]) / (end - start)

    @staticmethod
    def bit_count(input_val):
        res = 0
        while input_val != 0:
            res += input_val % 2
            input_val //= 2
        return res

    @staticmethod
    def n_choose_k(n, k):
        k = min(k, (n - k))
        if k <= 1:
            return 1 if k == 0 else n

        limit = 2 ** (n - 1)
        cnk = 0
        for i in range(3, limit):
            if TimeSeriesUtil.bit_count(i) == k:
                cnk += 1
        return cnk

    @staticmethod
    def split(points, min_length, length):
        new_points = [0] * (length + length)
        c = 0
        for i in range(length):
            segment_length = points[i] if i == 0 else points[i] - points[i - 1]

            if segment_length >= min_length * 2:
                start = 0 if i == 0 else points[i - 1]
                new_points[c] = start + segment_length // 2
                c += 1
            new_points[c] = points[i]
            c += 1

        return new_points[:c]

    @staticmethod
    def get_1_bit_changed_nums(x, n):
        return [x ^ (1 << i) for i in range(n)]

    @staticmethod
    def get_2_bits_changed_nums(x, n):
        res = []
        for i in range(1, n):
            mask = 1 << i
            temp = 1 << (i - 1)
            for j in range(i + 1, n + 1):
                real_mask = mask + temp
                res.append(x ^ real_mask)
                mask <<= 1
        return res

    @staticmethod
    def get_3_bits_changed_nums(x, n):
        res = []
        for i in range(1, n - 1):
            for j in range(i + 1, n):
                temp = (1 << (i - 1)) + (1 << (j - 1))
                mask = 1 << j
                for k in range(j + 1, n + 1):
                    real_mask = mask + temp
                    res.append(x ^ real_mask)
                    mask <<= 1
        return res

    @staticmethod
    def get_4_bits_changed_nums(x, n):
        res = []
        for i in range(1, n - 2):
            for j in range(i + 1, n - 1):
                for k in range(j + 1, n):
                    temp = (1 << (i - 1)) + (1 << (j - 1)) + (1 << (k - 1))
                    mask = 1 << k
                    for m in range(k + 1, n + 1):
                        real_mask = mask + temp
                        res.append(x ^ real_mask)
                        mask <<= 1
        return res

    @staticmethod
    def bit_diff_num(i, j, n):
        different_bits = i ^ j
        res = 0
        for k in range(n):
            if different_bits == 0:
                break
            tmp = different_bits >> 1
            if tmp + tmp + 1 == different_bits:
                res += 1
            different_bits = tmp
        return res
