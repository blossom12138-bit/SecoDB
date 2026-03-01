from typing import List
import math
import pandas as pd

class Item:
    def __init__(self, idx: int, w: float, r: float):
        self.idx = idx
        self.w = w      # 写次数（float）
        self.r = r      # 读次数（float）

class Bucket:
    def __init__(self, items: List[Item]):
        self.items = items
        self.rho = 0.0  # 混合更新密度

    def size(self):
        return len(self.items)

    def total_w(self):
        return sum(item.w for item in self.items)

    def total_r(self):
        return sum(item.r for item in self.items)


def compute_rho(bucket: Bucket, epsilon: float = 1e-6):
    w, r, s = bucket.total_w(), bucket.total_r(), bucket.size()
    rho = (w / max(1, math.sqrt(s))) * (1.0 / math.sqrt(r + epsilon))
    return rho


def dynamic_bucket_partition(buckets: List[Bucket],
                             tau_split: float,
                             tau_merge: float,
                             L_max: int,
                             L_min: int):
    """动态桶划分算法：先 split 再 merge（带最小桶大小约束）"""

    # Step 1: 更新 rho
    for b in buckets:
        b.rho = compute_rho(b)

    # Step 2: Splitting（加入 L_min 约束）
    new_buckets: List[Bucket] = []
    for b in buckets:
        # 只有当拆分后两个桶都 >= L_min 才允许拆分
        if (b.rho > tau_split and 
            b.size() >= 2 * L_min):
            mid = b.size() // 2
            b1 = Bucket(b.items[:mid])
            b2 = Bucket(b.items[mid:])
            b1.rho = compute_rho(b1)
            b2.rho = compute_rho(b2)
            new_buckets.extend([b1, b2])
        else:
            new_buckets.append(b)

    # Step 3: Merging（只检查相邻桶）
    merged: List[Bucket] = []
    i = 0
    n = len(new_buckets)

    while i < n:
        if i < n - 1:
            bi, bj = new_buckets[i], new_buckets[i + 1]
            if (bi.rho < tau_merge and
                bj.rho < tau_merge and
                bi.size() + bj.size() <= L_max):
                # merge
                new_bucket = Bucket(bi.items + bj.items)
                new_bucket.rho = compute_rho(new_bucket)
                merged.append(new_bucket)
                i += 2
                continue

        merged.append(new_buckets[i])
        i += 1

    return merged


def initialize_buckets(n_elements: int, p_buckets: int, scenario: str, read_write_ratio):
    items_per_bucket = n_elements // p_buckets
    buckets: List[Bucket] = []

    for i in range(p_buckets):
        items: List[Item] = []
        for j in range(items_per_bucket):
            idx = i * items_per_bucket + j

            if scenario == 'A':
                total_ops = 2.0  # 改为 float
                r = total_ops * read_write_ratio[0]
                w = total_ops * read_write_ratio[1]

            elif scenario == 'B':
                hotspot_ratio = 0.1
                total_ops = 1
                if i < int(p_buckets * hotspot_ratio):
                    r = 10.0 if read_write_ratio[0] > read_write_ratio[1] else 1.0
                    w = 10.0 if read_write_ratio[1] > read_write_ratio[0] else 1.0
                else:
                    r = 0.5
                    w = 0.5
            else:
                r, w = 1.0, 0.0

            items.append(Item(idx, w, r))

        buckets.append(Bucket(items))

    return buckets


def simulate_scenario_A_extended(n_elements, p_buckets,
                                 tau_split_list,
                                 tau_merge_list,
                                 L_max,
                                 L_min,
                                 cycles):

    ratios = [(0.9, 0.1), (0.5, 0.5), (0.1, 0.9)]

    for ratio in ratios:
        results = []
        read_r, write_r = ratio

        for tau_split in tau_split_list:
            for tau_merge in tau_merge_list:

                buckets = initialize_buckets(n_elements, p_buckets, 'A', ratio)

                row = {
                    'ratio': f'Read{read_r}_Write{write_r}',
                    'tau_split': tau_split,
                    'tau_merge': tau_merge
                }

                for t in range(1, cycles + 1):
                    buckets = dynamic_bucket_partition(buckets, tau_split, tau_merge, L_max,L_min)
                    row[f'cycle_{t}'] = len(buckets)

                results.append(row)

        df = pd.DataFrame(results)
        filename = f"scenario_A_Read{read_r}_Write{write_r}.csv"
        df.to_csv(filename, index=False)
        print(f"Saved: {filename}")


if __name__ == '__main__':
    n = 2**20
    p = 2**10
    L_max = 2**12
    L_min = 2**5
    cycles = 10

    tau_split_list = [2]
    tau_merge_list = [0.35]

    simulate_scenario_A_extended(
        n, p,
        tau_split_list,
        tau_merge_list,
        L_max,
        L_min,
        cycles
    )
