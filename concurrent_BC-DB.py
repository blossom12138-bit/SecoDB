import random
import statistics
import math
import pandas as pd


# =========================
# 全局参数
# =========================
NUM_TX = 10000
BLOCK_CAPACITY = 200

# DB 执行时间
READ_TIME = 0.001
WRITE_TIME = 0.032+0.0012  #VC更新证明+承诺 （向量大小L=2^18）

# 验证时间
VERIFY_READ = 0.0005
VERIFY_WRITE = 0.0012

# 并行控制开销
TSO_OVERHEAD = 0.0008
COMMIT_SERIAL_TIME = 0.002


# =========================
# 严格 slot 驱动的区块链
# =========================
class Blockchain:
    def __init__(self, slot_time):
        self.slot = slot_time
        self.queue = []           # (tx_id, ready_time)
        self.confirm_time = {}

    def submit(self, tx_id, ready_time):
        self.queue.append((tx_id, ready_time))

    def run(self):
        self.queue.sort(key=lambda x: x[1])

        block_id = 0
        i = 0
        n = len(self.queue)

        while i < n:
            block_time = (block_id + 1) * self.slot

            cnt = 0
            while cnt < BLOCK_CAPACITY and i < n:
                tx_id, ready = self.queue[i]
                if ready <= block_time:
                    self.confirm_time[tx_id] = block_time
                    i += 1
                    cnt += 1
                else:
                    break

            block_id += 1

        self.total_time = block_id * self.slot

    def throughput(self):
        return len(self.confirm_time) / self.total_time


# =========================
# 串行模型
# =========================
def simulate_serial(read_ratio, slot_time):
    bc = Blockchain(slot_time)

    t = 0.0
    submit_time = 0.0
    latencies = []

    for tx_id in range(NUM_TX):
        is_read = random.random() < read_ratio

        exec_time = READ_TIME if is_read else WRITE_TIME
        verify = VERIFY_READ if is_read else VERIFY_WRITE

        t += exec_time
        ready = t + verify
        bc.submit(tx_id, ready)

    bc.run()

    for tx_id in range(NUM_TX):
        latency = bc.confirm_time[tx_id] - submit_time
        latencies.append(latency)

    return bc.throughput(), statistics.mean(latencies)


# =========================
# 并行模型（MVCC + TSO）
# =========================
def simulate_parallel(read_ratio, slot_time):
    bc = Blockchain(slot_time)

    submit_time = 0.0
    exec_finish = []
    write_commits = []

    # ---- 并行执行阶段 ----
    for tx_id in range(NUM_TX):
        is_read = random.random() < read_ratio

        exec_time = READ_TIME if is_read else WRITE_TIME
        verify = VERIFY_READ if is_read else VERIFY_WRITE

        finish = exec_time + verify + TSO_OVERHEAD
        exec_finish.append((tx_id, is_read, finish))

        if not is_read:
            write_commits.append((tx_id, finish))

    # ---- 写事务 commit----
    write_commits.sort(key=lambda x: x[1])
    commit_ts = {}

    commit_time = 0.0
    for tx_id, ready in write_commits:
        commit_time = max(commit_time, ready)
        commit_time += COMMIT_SERIAL_TIME
        commit_ts[tx_id] = commit_time

    # ---- 提交区块链 ----
    for tx_id, is_read, finish in exec_finish:
        ready = finish if is_read else commit_ts[tx_id]
        bc.submit(tx_id, ready)

    bc.run()

    latencies = [
        bc.confirm_time[tx_id] - submit_time
        for tx_id in range(NUM_TX)
    ]

    return bc.throughput(), statistics.mean(latencies)


# =========================
# 实验主函数
# =========================
def experiment():
    slot_times = [0.4, 1, 2, 12]
    #read_ratios = [0,0.1, 0.2,0.3, 0.4,0.5,0.6, 0.7,0.8, 0.9,1]
    read_ratios = [0.1,0.3,0.5, 0.7, 0.9]
    print("slot  read  serial_tps  parallel_tps  serial_lat  parallel_lat")

    results = []  # 👉 用于存 Excel

    for slot in slot_times:
        for r in read_ratios:
            st, sl = simulate_serial(r, slot)
            pt, pl = simulate_parallel(r, slot)
            gain = (pt - st) / st *100 if st > 0 else float("nan")

            print(f"{slot:<4}  {r:<4.1f}  "
                  f"{st:<10.2f}  {pt:<12.2f}  "
                  f"{gain:<6.2f}  "
                  f"{sl:<10.4f}  {pl:<10.4f}")

            results.append({
                "slot_time(s)": slot,
                "read_ratio": r,
                "serial_tps": st,
                "parallel_tps": pt,
                "throughput_gain(%)": gain,
                "serial_latency": sl,
                "parallel_latency": pl
            })

    # 👉 写入 Excel
    df = pd.DataFrame(results)
    df.to_excel("concunrrent-BC-DB.xlsx", index=False)

    print("\nResults saved to concunrrent-BC-DB.xlsx")

if __name__ == "__main__":
    experiment()
