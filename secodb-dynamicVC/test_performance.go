package main

import (
	"fmt"
	"math"
	"math/rand"
	"time"

	"encoding/csv"
	"os"

	"github.com/alinush/go-mcl"
	"github.com/wangnick2017/balanceproofs-go/asvc"
	"github.com/wangnick2017/balanceproofs-go/fft"
	"github.com/wangnick2017/balanceproofs-go/vcs"
)

func main() {
	// 初始化 FFT（必须）
	fft.InitGlobals()
	rand.Seed(time.Now().Unix())

	//testStatic()
	testBucket()

	//test_bucket_split_by_size()
	//test_bucket_merge_by_size()

}

func randomVector(n uint64) []mcl.Fr {
	v := make([]mcl.Fr, n)
	for i := uint64(0); i < n; i++ {
		v[i].SetInt64(int64(rand.Intn(1000000)))
	}
	return v
}

func proofSize(p vcs.BasicBucketProof) int {
	b1 := p.BucketProof.Serialize()
	b2 := p.IndividualProof.Serialize()
	return len(b1) + len(b2)
}

func testStatic() {
	file, err := os.Create("static_test.csv")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// CSV header（新增 SingleProofSize 与 AggProofSize）
	header := []string{
		"L", "N",
		"Commit(ms)", "OpenAll(ms)", "Aggregate(ms)",
		"UpdateAll(ms)", "AvgUpdate(ms)",
		"UpdateCommitment(ms)", "AvgUpdateCommit(ms)",
		"AggProofSize(bytes)",
		"VerifySingle(ms)", "VerifyAggregate(ms)",
	}
	writer.Write(header)

	Ls := []uint8{16, 18, 20, 22}

	for _, L := range Ls {
		N := uint64(1) << L
		p := uint64(math.Sqrt(float64(N)))

		fmt.Println("Running static test L =", L)
		fmt.Println("p =", p)

		vb := &vcs.VBAS_B1{}
		vb.Init(L, p)
		fmt.Println("Complete Init")

		vector := randomVector(N)

		// ----------------------
		// Commit
		// ----------------------
		t := time.Now()
		digest := vb.Commit(vector)
		tCommit := time.Since(t).Milliseconds()
		fmt.Println("Complete Commit, time =", tCommit, "ms")

		// ----------------------
		// OpenAll
		// ----------------------
		// t = time.Now()
		// proofs := vb.OpenAll(vector)
		// tOpenAll := time.Since(t).Milliseconds()
		// fmt.Println("Complete OpenAll, time =", tOpenAll/(1000*60), "min")

		// ----------------------
		// Skip OpenAll —— 跳过整个 OpenAll 计算
		// ----------------------
		tOpenAll := int64(0)
		fmt.Println("Skip OpenAll (set time = 0)")

		// 构造一个最小可用的 proofs 数组
		proofs := make([]vcs.BasicBucketProof, N)
		for i := uint64(0); i < N; i++ {
			proofs[i] = vcs.BasicBucketProof{
				Index:           i,
				BucketProof:     mcl.G1{}, // empty
				IndividualProof: mcl.G1{}, // empty
			}
		}

		// ----------------------
		// Aggregate 1000 proofs
		// ----------------------
		randomIndices := make([]uint64, 0, 1000)
		for len(randomIndices) < 1000 && len(randomIndices) < int(N) {
			idx := uint64(rand.Intn(int(N)))
			randomIndices = append(randomIndices, idx)
		}

		randomProofs := make([]vcs.BasicBucketProof, len(randomIndices))
		subsetVector := make([]mcl.Fr, len(randomIndices))
		for i, idx := range randomIndices {
			randomProofs[i] = proofs[idx]
			subsetVector[i] = vector[idx]
		}

		t = time.Now()
		batchProofs := vb.Aggregate(randomProofs)
		tAgg := time.Since(t).Milliseconds()
		fmt.Println("Complete Aggregate, time =", tAgg, "ms")

		// -----------------------------
		// 计算单个证明和聚合证明大小
		// -----------------------------
		aggProofBytes := 0
		if len(batchProofs) > 0 {
			aggProofBytes = ProofSize(batchProofs[0])
		}
		fmt.Println("Batch  Proof Size =", aggProofBytes, "bytes")

		// ===============================
		// 1024-Random-Updates Benchmark（每次随机桶内更新）
		// ===============================
		updateCount := 1024
		updates := make([]asvc.UpdateReq, updateCount)

		pDetail := vb.GetPDetail()

		for i := 0; i < updateCount; i++ {
			// 随机选择桶
			bid := rand.Intn(int(vb.GetP()))
			bucketSize := int(p)
			if bucketSize == 0 {
				i-- // 空桶，重试
				continue
			}

			bucketStart := pDetail[bid][0] // 桶起始位置

			// 桶内随机位置
			posInBucket := rand.Intn(bucketSize)
			idx := bucketStart + uint64(posInBucket)

			var delta mcl.Fr
			delta.Random()
			updates[i] = asvc.UpdateReq{Index: idx, Delta: delta}
		}

		// 按桶整理更新请求
		aux := make([][]asvc.UpdateReq, vb.GetP())
		for _, u := range updates {
			bid := int(u.Index % uint64(vb.GetP()))
			aux[bid] = append(aux[bid], u)
		}

		// Batch UpdateAll
		t = time.Now()
		proofs, vector, aux = vb.UpdateAll(proofs, vector, updates[0], aux)
		tUpdateAll := time.Since(t).Milliseconds()
		avgUpdate := float64(tUpdateAll) / float64(updateCount)
		fmt.Printf("Complete UpdateAll: total=%d ms, avg=%.4f ms/update\n", tUpdateAll, avgUpdate)

		// 1024 次 UpdateCommitment()
		t = time.Now()
		newDigest := digest
		for _, u := range updates {
			newDigest = vb.UpdateCommitment(newDigest, u)
		}
		tUpdateCommitment := time.Since(t).Milliseconds()
		avgUC := float64(tUpdateCommitment) / float64(updateCount)
		fmt.Printf("Complete UpdateCommitment: total=%d ms, avg=%.4f ms/update\n", tUpdateCommitment, avgUC)

		// ----------------------
		// VerifySingle
		// ----------------------
		t = time.Now()
		_ = vb.VerifySingle(digest, proofs[0], vector[1], 1)
		tVerifySingle := time.Since(t).Milliseconds()
		fmt.Println("Complete VerifySingle, time =", tVerifySingle, "ms")

		// ----------------------
		// VerifyAggregate
		// ----------------------
		t = time.Now()
		_ = vb.VerifyAggregation(newDigest, batchProofs, subsetVector)
		tVerifyAgg := time.Since(t).Milliseconds()
		fmt.Println("Complete VerifyAggregate, time =", tVerifyAgg, "ms")

		// ===============================
		// 写入 CSV
		// ===============================
		row := []string{
			fmt.Sprint(L),
			fmt.Sprint(N),
			fmt.Sprint(tCommit),
			fmt.Sprint(tOpenAll),
			fmt.Sprint(tAgg),
			fmt.Sprint(tUpdateAll),
			fmt.Sprintf("%.4f", avgUpdate),
			fmt.Sprint(tUpdateCommitment),
			fmt.Sprintf("%.4f", avgUC),
			fmt.Sprint(aggProofBytes),
			fmt.Sprint(tVerifySingle),
			fmt.Sprint(tVerifyAgg),
		}
		writer.Write(row)
	}

	fmt.Println("✔ static_test.csv generated")
}

func ProofSize(bp vcs.BasicBatchProof) int {
	size := 0

	// BucketProof 直接调用 Serialize()
	size += len(bp.BucketProof.Serialize())

	// 聚合证明也直接调用 Serialize()
	agg := bp.AggIndividualProof() // 返回值 mcl.G1
	size += len(agg.Serialize())

	// 索引数组，每个 uint64 8 字节
	size += len(bp.Indexs) * 8

	return size
}

func testBucket() {
	file, err := os.Create("bucket_test.csv")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// 表头
	header := []string{"p", "AvgUpdate(ms)", "AvgVerify(ms)", "ProofSize(KB)"}
	_ = writer.Write(header)

	L := uint8(20)
	N := uint64(1) << L

	base := uint64(1) << (L / 2)

	pList := []uint64{
		base * 32,
		base / 8,
		base / 4,
		base / 2,
		base,
		base * 2,
		base * 4,
		base * 8,
	}

	fmt.Println("===== Bucket Test Start =====")
	fmt.Println("L =", L, ", N =", N)
	fmt.Println("Candidate p =", pList)
	fmt.Println("=============================")

	for _, p := range pList {
		if p == 0 || p > N {
			continue
		}

		fmt.Println("Running p =", p)

		vb := &vcs.VBAS_B1{}
		vb.Init(L, p)

		vector := randomVector(N)
		digest := vb.Commit(vector)

		// 初始化 proofs
		proofs := make([]vcs.BasicBucketProof, N)
		for i := uint64(0); i < N; i++ {
			proofs[i] = vcs.BasicBucketProof{
				Index:           i,
				BucketProof:     mcl.G1{},
				IndividualProof: mcl.G1{},
			}
		}

		// ==================================================
		// ✅ 1024-Random-Updates + UpdateAll Benchmark
		// ==================================================
		updateCount := 1024
		updates := make([]asvc.UpdateReq, updateCount)

		pDetail := vb.GetPDetail()

		for i := 0; i < updateCount; i++ {
			bid := rand.Intn(int(vb.GetP()))
			bucketSize := int(p)
			if bucketSize == 0 {
				i--
				continue
			}

			bucketStart := pDetail[bid][0]
			posInBucket := rand.Intn(bucketSize)
			idx := bucketStart + uint64(posInBucket)

			var delta mcl.Fr
			delta.Random()

			updates[i] = asvc.UpdateReq{
				Index: idx,
				Delta: delta,
			}
		}

		// 按桶分组
		aux := make([][]asvc.UpdateReq, vb.GetP())
		for _, u := range updates {
			bid := int(u.Index % uint64(vb.GetP()))
			aux[bid] = append(aux[bid], u)
		}

		// 批量 UpdateAll
		t := time.Now()
		proofs, vector, aux = vb.UpdateAll(proofs, vector, updates[0], aux)
		tUpdateAll := time.Since(t).Milliseconds()

		avgUpdate := float64(tUpdateAll) / float64(updateCount)

		// ==================================================
		// ✅ VerifyAggregation（保持不变）
		// ==================================================
		batchSize := 1024
		if batchSize > int(N) {
			batchSize = int(N)
		}

		perm := rand.Perm(int(N))[:batchSize]

		randomProofs := make([]vcs.BasicBucketProof, batchSize)
		subsetVector := make([]mcl.Fr, batchSize)

		for i, idx := range perm {
			u := uint64(idx)
			randomProofs[i] = proofs[u]
			subsetVector[i] = vector[u]
		}

		batchProof := vb.Aggregate(randomProofs)

		t = time.Now()
		_ = vb.VerifyAggregation(digest, batchProof, subsetVector)
		avgVerify := float64(time.Since(t).Microseconds()) / 1000.0

		// ProofSize
		proofSizeKB := float64(48*p) / 1024.0

		row := []string{
			fmt.Sprint(p),
			fmt.Sprintf("%.4f", avgUpdate),
			fmt.Sprintf("%.4f", avgVerify),
			fmt.Sprintf("%.4f", proofSizeKB),
		}
		_ = writer.Write(row)

		fmt.Printf(
			"[p=%d] AvgUpdateAll = %.4f ms/update | VerifyAggregate(1024) = %.4f ms | ProofSize = %.4f KB\n",
			p, avgUpdate, avgVerify, proofSizeKB,
		)
	}

	fmt.Println("✔ bucket_test.csv generated")
}

func test_split() {
	fmt.Println("===== Bucket Split Test =====")

	L := uint8(16)
	N := uint64(1) << L
	p := uint64(1) << 8

	vb := &vcs.VBAS_B1{}
	vb.Init(L, p)

	// 初始化向量 & 承诺
	vector := randomVector(N)
	_ = vb.Commit(vector)

	// 初始化 proofs
	proofs := make([]vcs.BasicBucketProof, N)
	for i := uint64(0); i < N; i++ {
		proofs[i] = vcs.BasicBucketProof{
			Index:           i,
			BucketProof:     mcl.G1{},
			IndividualProof: mcl.G1{},
		}
	}

	// 随机选一个桶
	targetBucket := uint64(rand.Intn(int(vb.GetP())))

	fmt.Printf("Split bucket %d (size=%d) into 2 buckets\n",
		targetBucket,
		vb.GetPDetail()[targetBucket][1]-vb.GetPDetail()[targetBucket][0]+1,
	)

	t := time.Now()
	proofs = vb.BucketSplit(targetBucket, 2, proofs, vector)
	elapsed := time.Since(t).Milliseconds()

	fmt.Printf("BucketSplit finished: %d ms\n", elapsed)
	fmt.Printf("New p = %d\n", vb.GetP())
	fmt.Println("================================\n")
}

func test_merge() {
	fmt.Println("===== Bucket Merge Test =====")

	L := uint8(20)
	N := uint64(1) << L
	p := uint64(1) << 10 // 2^10 buckets

	vb := &vcs.VBAS_B1{}
	vb.Init(L, p)

	// 初始化向量 & 承诺
	vector := randomVector(N)
	_ = vb.Commit(vector)

	// 初始化 proofs
	proofs := make([]vcs.BasicBucketProof, N)
	for i := uint64(0); i < N; i++ {
		proofs[i] = vcs.BasicBucketProof{
			Index:           i,
			BucketProof:     mcl.G1{},
			IndividualProof: mcl.G1{},
		}
	}

	// 随机选两个相邻桶
	startBucket := uint64(rand.Intn(int(vb.GetP() - 1)))

	fmt.Printf("Merge bucket %d and %d\n", startBucket, startBucket+1)

	t := time.Now()
	proofs = vb.BucketMerge(startBucket, 2, proofs, vector)
	elapsed := time.Since(t).Milliseconds()

	fmt.Printf("BucketMerge finished: %d ms\n", elapsed)
	fmt.Printf("New p = %d\n", vb.GetP())
	fmt.Println("================================\n")
}

func test_bucket_merge_by_size() {
	fmt.Println("===== BucketMerge vs Bucket Size =====")

	L := uint8(20)
	N := uint64(1) << L

	fmt.Printf("Global vector size: N = 2^%d = %d\n\n", L, N)

	bucketSizes := []uint64{
		1 << 6,
		1 << 8,
		1 << 10,
		1 << 12,
		1 << 14,
	}

	// ===== CSV 初始化 =====
	file, err := os.Create("bucket_merge_by_size.csv")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// CSV header
	writer.Write([]string{
		"bucket_exp",
		"bucket_size",
		"p",
		"time_s",
	})

	for _, bucketSize := range bucketSizes {
		p := N / bucketSize
		if p < 2 {
			continue
		}

		vb := &vcs.VBAS_B1{}
		vb.Init(L, p)

		vector := randomVector(N)
		_ = vb.Commit(vector)
		vb.InitPhiPoly(vector)

		proofs := make([]vcs.BasicBucketProof, N)
		for i := uint64(0); i < N; i++ {
			proofs[i] = vcs.BasicBucketProof{Index: i}
		}

		startBucket := uint64(rand.Intn(int(p - 1)))

		t := time.Now()
		proofs = vb.BucketMerge(startBucket, 2, proofs, vector)
		elapsed := time.Since(t).Seconds()

		exp := uint(math.Log2(float64(bucketSize)))

		// 控制台输出
		fmt.Printf(
			"BucketMerge | bucketSize = 2^%d | p = %d | time = %f s\n",
			exp, p, elapsed,
		)

		// CSV 写入
		writer.Write([]string{
			fmt.Sprintf("%d", exp),
			fmt.Sprintf("%d", bucketSize),
			fmt.Sprintf("%d", p),
			fmt.Sprintf("%f", elapsed),
		})
	}

	fmt.Println("======================================\n")
}

func test_bucket_split_by_size() {
	fmt.Println("===== BucketSplit vs Bucket Size =====")

	L := uint8(20)
	N := uint64(1) << L

	fmt.Printf("Global vector size: N = 2^%d = %d\n\n", L, N)

	bucketSizes := []uint64{
		1 << 6,
		1 << 8,
		1 << 10,
		1 << 12,
		1 << 14,
	}

	// ===== CSV 初始化 =====
	file, err := os.Create("bucket_split_by_size.csv")
	if err != nil {
		panic(err)
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	writer.Write([]string{
		"bucket_exp",
		"bucket_size",
		"p",
		"time_s",
	})

	for _, bucketSize := range bucketSizes {
		p := N / bucketSize

		vb := &vcs.VBAS_B1{}
		vb.Init(L, p)

		vector := randomVector(N)
		_ = vb.Commit(vector)
		vb.InitPhiPoly(vector)

		proofs := make([]vcs.BasicBucketProof, N)
		for i := uint64(0); i < N; i++ {
			proofs[i] = vcs.BasicBucketProof{Index: i}
		}

		targetBucket := uint64(rand.Intn(int(p)))

		t := time.Now()
		proofs = vb.BucketSplit(targetBucket, 2, proofs, vector)
		elapsed := time.Since(t).Seconds()

		exp := uint(math.Log2(float64(bucketSize)))

		fmt.Printf(
			"BucketSplit | bucketSize = 2^%d | p = %d | time = %f s\n",
			exp, p, elapsed,
		)

		writer.Write([]string{
			fmt.Sprintf("%d", exp),
			fmt.Sprintf("%d", bucketSize),
			fmt.Sprintf("%d", p),
			fmt.Sprintf("%f", elapsed),
		})
	}

	fmt.Println("======================================\n")
}
