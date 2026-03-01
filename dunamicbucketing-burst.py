from typing import List
import math
import random
import pandas as pd
import os
import csv


# ----------------------------------------------------
# Item + Bucket
# ----------------------------------------------------
class Item:
    def __init__(self, idx: int):
        self.idx = idx
        self.w = 1.0    # float baseline
        self.r = 1.0


class Bucket:
    def __init__(self, items: List[Item]):
        self.items = items
        self.rho = 0.0

    def size(self):
        return len(self.items)

    def total_w(self):
        return sum(item.w for item in self.items)

    def total_r(self):
        return sum(item.r for item in self.items)


# ----------------------------------------------------
# rho + split/merge
# ----------------------------------------------------
def compute_rho(bucket: Bucket, epsilon: float = 1e-6):
    w, r, s = bucket.total_w(), bucket.total_r(), bucket.size()
    rho = (w / max(1, math.sqrt(s))) * (1.0 / math.sqrt(r + epsilon))
    return rho


def dynamic_bucket_partition(buckets: List[Bucket],
                             tau_split: float,
                             tau_merge: float,
                             L_max: int,
                             L_min: int):

    for b in buckets:
        b.rho = compute_rho(b)

    new_buckets: List[Bucket] = []
    for b in buckets:
        if (b.rho > tau_split and b.size() >= 2 * L_min):
            mid = b.size() // 2
            b1 = Bucket(b.items[:mid])
            b2 = Bucket(b.items[mid:])
            b1.rho = compute_rho(b1)
            b2.rho = compute_rho(b2)
            new_buckets.extend([b1, b2])
        else:
            new_buckets.append(b)

    merged: List[Bucket] = []
    i = 0
    n = len(new_buckets)

    while i < n:
        if i < n - 1:
            bi, bj = new_buckets[i], new_buckets[i + 1]
            if (bi.rho < tau_merge and
                bj.rho < tau_merge and
                bi.size() + bj.size() <= L_max):
                new_bucket = Bucket(bi.items + bj.items)
                new_bucket.rho = compute_rho(new_bucket)
                merged.append(new_bucket)
                i += 2
                continue

        merged.append(new_buckets[i])
        i += 1

    return merged


# ----------------------------------------------------
# Workloads
# ----------------------------------------------------
def reset_all_items_to_baseline(buckets):
    for b in buckets:
        for item in b.items:
            item.w = 1.0
            item.r = 1.0


def apply_burst_write(buckets, burst_hot_ratio, burst_multiplier):
    reset_all_items_to_baseline(buckets)
    num_hot_buckets = max(1, int(len(buckets) * burst_hot_ratio))
    for b in random.sample(buckets, num_hot_buckets):
        for item in b.items:
            item.w = burst_multiplier


def apply_burst_read(buckets, burst_hot_ratio, burst_multiplier):
    reset_all_items_to_baseline(buckets)
    num_hot_buckets = max(1, int(len(buckets) * burst_hot_ratio))
    for b in random.sample(buckets, num_hot_buckets):
        for item in b.items:
            item.r = burst_multiplier


# ----------------------------------------------------
# Init
# ----------------------------------------------------
def init_buckets(n, p):
    items = [Item(i) for i in range(n)]
    L = n // p
    return [Bucket(items[i*L:(i+1)*L]) for i in range(p)]


def print_bucket_stats(epoch, buckets, csv_path="bucket_stats.csv"):
    """
    将指定 epoch 的每个桶统计信息写入 CSV
    CSV 列: epoch, bucket_id, size, total_r, total_w
    """

    file_exists = os.path.isfile(csv_path)

    with open(csv_path, mode="a", newline="") as f:
        writer = csv.writer(f)

        # 只在文件不存在时写表头
        if not file_exists:
            writer.writerow(["epoch", "bucket_id", "size", "total_r", "total_w"])

        for i, b in enumerate(buckets):
            writer.writerow([
                epoch,
                i,
                b.size(),
                round(b.total_r(), 6),
                round(b.total_w(), 6)
            ])


# ----------------------------------------------------
# Single simulation
# ----------------------------------------------------
def simulate_num_buckets(T, n, p,
                         tau_split, tau_merge, L_max, L_min,
                         write_burst_window, read_burst_window,
                         burst_hot_ratio, burst_multiplier):

    buckets = init_buckets(n, p)
    nums = []

    snapshot_epochs = {5, 10, 25, 30}

    for t in range(T):

        if write_burst_window[0] <= t <= write_burst_window[1]:
            apply_burst_write(buckets, burst_hot_ratio, burst_multiplier)
        elif read_burst_window[0] <= t <= read_burst_window[1]:
            apply_burst_read(buckets, burst_hot_ratio, burst_multiplier)
        else:
            reset_all_items_to_baseline(buckets)

        buckets = dynamic_bucket_partition(
            buckets, tau_split, tau_merge, L_max, L_min
        )

        nums.append(len(buckets))

        if t in snapshot_epochs:
            print_bucket_stats(t, buckets)

    return nums


# ----------------------------------------------------
# Run experiments
# ----------------------------------------------------
def run_and_save_compare():
    random.seed(42)

    T = 50
    n = 53248
    p = 52
    L_max = 2**12
    L_min = 32

    write_burst_window = (0, 0)
    read_burst_window = (25, 30)

    burst_hot_ratio = 1
    burst_multiplier = 100.0

    param_list = [
        (2, 0.8),
        (4, 0.8),
        (2, 0.5),
        (4, 0.5)
      

    ]

    results = {}
    for idx, (tau_split, tau_merge) in enumerate(param_list, start=1):
        print(f"\n######## Experiment {idx} ########")
        random.seed(42 + idx)
        nums = simulate_num_buckets(
            T, n, p,
            tau_split, tau_merge, L_max, L_min,
            write_burst_window, read_burst_window,
            burst_hot_ratio, burst_multiplier
        )
        results[f"num_buckets_exp{idx}"] = nums

    df = pd.DataFrame({"time": list(range(T))})
    for col, series in results.items():
        df[col] = series

    df.to_csv("experiment_compare.csv", index=False)
    print("\nSaved → experiment_compare.csv")


if __name__ == "__main__":
    run_and_save_compare()
