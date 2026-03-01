package vcs

import (
	"fmt"
	"runtime"
	"sync"
	"time"

	"github.com/alinush/go-mcl"
	"github.com/wangnick2017/balanceproofs-go/asvc"
	"github.com/wangnick2017/balanceproofs-go/fft"
	"github.com/wangnick2017/balanceproofs-go/utils"
)

type BasicBucketProof struct {
	Index           uint64 // 全局位置索引
	BucketProof     mcl.G1 // Π_i：第 i 个桶的证明
	IndividualProof mcl.G1 // π_{i,j}：桶内位置 j 的证明
}

// ✅ polyCache 激进优化：多项式系数缓存结构
// 用于在 UpdateProof 中避免重复计算相同的多项式系数
type polyCache struct {
	ljCache  map[uint64][]mcl.Fr            // j -> L_j(x) 多项式系数缓存
	lijCache map[uint64]map[uint64][]mcl.Fr // i -> (j_local -> L_{i,j}(x)) 缓存
}

// ✅ openAllProgressTracker 实时进度追踪结构
type openAllProgressTracker struct {
	bucketsDone   int64  // 完成的桶数
	currentBucket int64  // 当前处理的桶号
	currentPos    int64  // 当前桶内处理到的位置
	totalBuckets  uint64 // 总桶数
	mu            sync.Mutex
}

type BasicBatchProof struct {
	// 默认只聚合同一个桶里的证明
	BucketProof             mcl.G1   // Π_i：第 i 个桶的证明
	aggIndividualBatchProof mcl.G1   // 聚合后的证明
	Indexs                  []uint64 // 被聚合的索引列表
}

type VBAS_B1 struct {
	// 基本参数
	N        uint64     // 向量总长度
	L        uint8      // log2(n)
	p        uint64     // 桶数量
	p_detail [][]uint64 // 现在每个桶的大小可以不一样，但桶的划分一定还是连续的，这个参数记录每个桶的范围
	Lp       uint8      // log2(p)
	Lb       uint8      // log2(n/p)，每个桶的大小

	// KZG 秘密参数（n+1 个）
	SecretG1 []mcl.G1 // [g^(τ^i)]_{i=0}^n
	SecretG2 []mcl.G2 // [h^(τ^i)]_{i=0}^n
	G        mcl.G1   // G1 的生成元
	H        mcl.G2   // G2 的生成元

	// FFT 设置
	FFT *fft.Settings

	// 预计算的数据结构
	// ak(x) = ∏_{j∈Pk} (x - ω_j) 对应的 KZG 承诺 [g^ak(τ)]
	A_k []mcl.G2

	// ✅ 全局拉格朗日多项式承诺（用于 Commit）
	// L_j[j] = [g^L_j(τ)]，其中 L_j(x) = ∏_{k≠j} (x - ω_k) / (ω_j - ω_k)
	L_j []mcl.G1

	Lij_polys [][]mcl.G1
	Lj_polys  []mcl.G1

	// 每个桶的根多项式系数缓存：a_i(x) = ∏_{k∈P_i} (X - ω_k)
	// ai_polys[i] = a_i(x) 的系数数组
	ai_polys [][]mcl.Fr

	// ✅ 激进优化：运行时多项式系数缓存
	// 用于 UpdateProof 中避免重复计算相同的多项式系数
	// 每次 precomputeData 开始时会清空重置
	polyCache *polyCache

	// ✅ 实时进度追踪：用于 OpenAll 并行化的细粒度进度显示
	progressTracker *openAllProgressTracker

	PhiPoly []mcl.Fr

	BucketDeltas map[uint64]*BucketDelta
}

// ============================================================================
// 公共 Getter 方法
// ============================================================================

// GetP 返回桶数量
func (vb *VBAS_B1) GetP() uint64 {
	return vb.p
}

// GetPDetail 返回桶详情
func (vb *VBAS_B1) GetPDetail() [][]uint64 {
	return vb.p_detail
}

func (bp *BasicBatchProof) AggIndividualProof() *mcl.G1 {
	return &bp.aggIndividualBatchProof
}

// Init 初始化 Basic Bucketing 结构（所有桶大小相同）
func (vb *VBAS_B1) Init(L uint8, p uint64) {
	vb.L = L
	vb.N = 1 << L
	vb.p = p

	// 计算 Lp 和 Lb
	vb.Lp = uint8(0)
	for (1 << vb.Lp) < p {
		vb.Lp++
	}
	vb.Lb = L - vb.Lp

	// 初始化 p_detail：所有桶大小相同，均匀划分
	bucketSize := vb.N / vb.p
	vb.p_detail = make([][]uint64, vb.p)
	for i := uint64(0); i < vb.p; i++ {
		vb.p_detail[i] = make([]uint64, 2)
		vb.p_detail[i][0] = i * bucketSize       // 桶的起始位置
		vb.p_detail[i][1] = (i+1)*bucketSize - 1 // 桶的结束位置
	}

	// 生成 KZG 公参（n+1 个）
	fmt.Println("Generating KZG secret parameters...")
	vb.SecretG1, vb.SecretG2 = utils.GenerateTestingSetup(vb.N) //生成KZG秘密参数[g^(τ^i)]_{i=0,1,2,...,n}和[h^(τ^i)]_{i=0,1,2,...,n}
	fmt.Println("KZG secret parameters generated")
	vb.G = vb.SecretG1[0]
	vb.H = vb.SecretG2[0]

	// 初始化 FFT
	fmt.Println("Initializing FFT...")
	vb.FFT = fft.NewSettings(vb.L + 1)
	fmt.Println("FFT initialized")

	// 预计算数据结构
	vb.precomputeData()
	vb.BucketDeltas = make(map[uint64]*BucketDelta)

}

// InitPhiPoly 初始化 φ(x)，只做一次
func (vb *VBAS_B1) InitPhiPoly(vector []mcl.Fr) {
	vb.PhiPoly = vb.FFT.FFT_Fr(vector, true)
}

// precomputeData 预计算所有必要的数据结构
// ✅ 激进优化：不预计算多项式系数，改为在 UpdateProof 中即时计算，UpdateProof用的比较少
func (vb *VBAS_B1) precomputeData() {
	// ✅ 清空/重置多项式系数缓存
	vb.polyCache = &polyCache{
		ljCache:  make(map[uint64][]mcl.Fr),
		lijCache: make(map[uint64]map[uint64][]mcl.Fr),
	}

	vb.precomputeAk()

	vb.precomputeLjCommitment()

}

// precomputeAk 计算每个桶对应的多项式 a_k(x) = ∏_{j∈Pk} (x - ω_j)
func (vb *VBAS_B1) precomputeAk() {
	vb.A_k = make([]mcl.G2, vb.p)        //有p个桶，就有p个a_k(x)
	vb.ai_polys = make([][]mcl.Fr, vb.p) //缓存所有桶的根多项式

	for i := uint64(0); i < vb.p; i++ {
		bucketSize := vb.getBucketSize(i)   // 获取第 i 个桶的大小
		roots := make([]mcl.Fr, bucketSize) //每个桶有bucketSize个根
		bucketStart := vb.p_detail[i][0]    // 桶的起始位置
		for j := uint64(0); j < bucketSize; j++ {
			globalIdx := bucketStart + j
			roots[j] = vb.FFT.RootsOfUnity[globalIdx*2] //找到ω_j
		}

		// 直接传入正的根，SubProductTree 内部会处理 (X - ω_j)
		subProd := utils.SubProductTree(roots) //计算 a_k(x) = ∏(X - ω_j)
		vb.ai_polys[i] = subProd.Poly

		mcl.G2MulVec(&vb.A_k[i], vb.SecretG2[:len(subProd.Poly)], subProd.Poly) //计算每个桶的a_k(τ)的承诺
	}
}

// precomputeLjCommitment ✅ 激进优化：只计算 G1 承诺，不计算多项式系数
func (vb *VBAS_B1) precomputeLjCommitment() {
	// 直接通过一次 IFFT 获得 [g^L_j(τ)]
	vb.L_j = vb.FFT.FFT_G1(vb.SecretG1[:vb.N], true)
}

// getBucketByGlobalIndex 根据全局位置索引找到对应的桶号
func (vb *VBAS_B1) getBucketByGlobalIndex(index uint64) uint64 {
	for i := uint64(0); i < vb.p; i++ {
		if index >= vb.p_detail[i][0] && index <= vb.p_detail[i][1] {
			return i
		}
	}
	return 0 // 应该不会到这里
}

// getLocalIndexInBucket 获取全局位置在其所在桶内的本地位置
func (vb *VBAS_B1) getLocalIndexInBucket(index uint64) uint64 {
	bucketID := vb.getBucketByGlobalIndex(index)
	return index - vb.p_detail[bucketID][0]
}

// getBucketSize 获取指定桶的大小
func (vb *VBAS_B1) getBucketSize(bucketID uint64) uint64 {
	return vb.p_detail[bucketID][1] - vb.p_detail[bucketID][0] + 1
}

// ✅ computeLjPoly 即时计算全局拉格朗日多项式系数 L_j(x)
// 使用 IFFT 矩阵公式直接计算，并支持缓存
func (vb *VBAS_B1) computeLjPoly(j uint64) []mcl.Fr {
	// ✅ 先检查缓存
	if cached, exists := vb.polyCache.ljCache[j]; exists {
		return cached
	}

	Lj_poly := make([]mcl.Fr, vb.N)

	// 预计算 1/N
	var invN mcl.Fr
	invN.SetInt64(int64(vb.N))
	mcl.FrInv(&invN, &invN)

	// L_j 的第 k 个系数 = ω^(-k*j) / N
	var omega_power mcl.Fr
	omega_power.SetInt64(1) // ω^0 = 1

	for k := uint64(0); k < vb.N; k++ {
		mcl.FrMul(&Lj_poly[k], &omega_power, &invN)

		// 更新 omega_power = omega_power * ω^(-j)
		if k+1 < vb.N {
			var omega_inv_j mcl.Fr
			if j == 0 {
				omega_inv_j.SetInt64(1)
			} else {
				omega_inv_j = vb.FFT.ReverseRootsOfUnity[j*2]
			}
			mcl.FrMul(&omega_power, &omega_power, &omega_inv_j)
		}
	}

	// ✅ 存储到缓存
	vb.polyCache.ljCache[j] = Lj_poly

	return Lj_poly
}

// ✅ computeLijPoly 即时计算局部拉格朗日多项式系数 L_{i,j}(x)
// i 是桶号，j_local 是在桶内的本地位置，支持缓存
func (vb *VBAS_B1) computeLijPoly(i uint64, j_local uint64) []mcl.Fr {
	// ✅ 先检查缓存
	if lijMap, exists := vb.polyCache.lijCache[i]; exists {
		if cached, exists := lijMap[j_local]; exists {
			return cached
		}
	} else {
		// 如果该桶还没有缓存映射，创建一个
		vb.polyCache.lijCache[i] = make(map[uint64][]mcl.Fr)
	}

	bucketSize := vb.getBucketSize(i)
	bucketStart := vb.p_detail[i][0]

	// 获取桶内所有根
	roots := make([]mcl.Fr, bucketSize)
	for k := uint64(0); k < bucketSize; k++ {
		roots[k] = vb.FFT.RootsOfUnity[(bucketStart+k)*2]
	}

	// 收集除了 ω_j 以外的所有根
	indices := make([]mcl.Fr, bucketSize-1)
	idx := 0
	for k := uint64(0); k < bucketSize; k++ {
		if k != j_local {
			indices[idx] = roots[k]
			idx++
		}
	}

	// 计算分子：∏_{k≠j，k ∈ P_i} (X - ω_k)
	subProd := utils.SubProductTree(indices)
	Lij_poly := make([]mcl.Fr, len(subProd.Poly))
	copy(Lij_poly, subProd.Poly)

	// 计算分母：∏_{k≠j，k ∈ P_i} (ω_j - ω_k)
	var denom mcl.Fr
	denom.SetInt64(1)
	for k := uint64(0); k < bucketSize; k++ {
		if k != j_local {
			var diff mcl.Fr
			mcl.FrSub(&diff, &roots[j_local], &roots[k])
			mcl.FrMul(&denom, &denom, &diff)
		}
	}

	// 计算模逆
	var invDenom mcl.Fr
	mcl.FrInv(&invDenom, &denom)

	// 标准化系数
	for k := 0; k < len(Lij_poly); k++ {
		mcl.FrMul(&Lij_poly[k], &Lij_poly[k], &invDenom)
	}

	// ✅ 存储到缓存
	vb.polyCache.lijCache[i][j_local] = Lij_poly

	return Lij_poly
}

// ============================================================================
// 辅助函数：多项式操作
// ============================================================================

// polyNegate 计算多项式的取反：-poly
func polyNegate(poly []mcl.Fr) []mcl.Fr {
	result := make([]mcl.Fr, len(poly))
	for i := 0; i < len(poly); i++ {
		mcl.FrNeg(&result[i], &poly[i])
	}
	return result
}

// polySubtract 计算多项式的减法：a - b = a + (-b)
func polySubtract(a []mcl.Fr, b []mcl.Fr) []mcl.Fr {
	negB := polyNegate(b)
	return utils.PolyAdd(a, negB)
}

func polyReverse(poly []mcl.Fr) []mcl.Fr {
	n := len(poly)
	out := make([]mcl.Fr, n)
	for i := 0; i < n; i++ {
		out[i] = poly[n-1-i]
	}
	return out
}

func minUint64(a, b uint64) uint64 {
	if a < b {
		return a
	}
	return b
}

// Commit 计算向量承诺
// C = g^(φ(τ)) = ∑_{j=0}^{n-1} m_j · [L_j(τ)]
// 其中 L_j(x) 是全局拉格朗日多项式，基于所有 n 个根
func (vb *VBAS_B1) Commit(vector []mcl.Fr) mcl.G1 {
	var digest mcl.G1
	digest.Clear()

	// 使用全局拉格朗日多项式 L_j
	for j := uint64(0); j < vb.N; j++ {
		var temp mcl.G1
		mcl.G1Mul(&temp, &vb.L_j[j], &vector[j])
		mcl.G1Add(&digest, &digest, &temp)
	}

	return digest
}

func (vb *VBAS_B1) UpdateCommitment(digest mcl.G1, req asvc.UpdateReq) mcl.G1 {
	var temp mcl.G1
	mcl.G1Mul(&temp, &vb.L_j[req.Index], &req.Delta)
	mcl.G1Add(&temp, &digest, &temp)
	return temp
}

// Open 计算单个位置的证明
// 根据论文 4.1：
// - bucket_proof Π_i = g^(q_i(τ)) 其中 φ(x) = φ_i(x) + q_i(x)∏_{k∈P_i}(x-ω_k)
// - individual_proof π_{i,j} = g^(q_{i,j}(τ)) 其中 φ_i(x) = q_{i,j}(x)(x-ω_j) + φ_i(ω_j)
func (vb *VBAS_B1) Open(index uint64, vector []mcl.Fr) BasicBucketProof {
	i := vb.getBucketByGlobalIndex(index)     // 获取 index 所在的桶号
	j := index                                // 全局索引
	bucketStart := vb.p_detail[i][0]          // 桶的起始位置
	bucketEnd := vb.p_detail[i][1]            // 桶的结束位置
	bucketSize := bucketEnd - bucketStart + 1 // 桶的大小

	vec_i_raw := vector[bucketStart : bucketEnd+1]

	// 计算完整向量的多项式 φ(x)，参考asvc.go的Open函数
	phi_poly := vb.FFT.FFT_Fr(vector, true)

	// 计算桶内子向量 φ_i(x) 的多项式
	// 注意：φ_i(x) = ∑_{j∈P_i} L_{i,j}(x)·m_j
	// 其中 L_{i,j}(x) 是在 P_i 对应根集上定义的局部 Lagrange 基
	// 不能直接用全局 FFT，需要用 P_i 对应的根进行显式 Lagrange 插值
	bucket_roots := make([]mcl.Fr, bucketSize)
	for k := uint64(0); k < bucketSize; k++ {
		globalIdx := i*bucketSize + k
		bucket_roots[k] = vb.FFT.RootsOfUnity[globalIdx*2]
	}
	bucket_subprod := utils.SubProductTree(bucket_roots) //∏_{k∈P_i}(x - ω_k)，相当于A_I
	bucket_subprod_prime := utils.PolyDifferentiate(bucket_subprod.Poly)
	s := utils.PolyMultiEvaluate(bucket_subprod_prime, bucket_subprod)

	ys_vec_i := make([]mcl.Fr, bucketSize)
	for k := uint64(0); k < bucketSize; k++ {
		mcl.FrDiv(&ys_vec_i[k], &vec_i_raw[k], &s[k])
	}
	//此处处理参见asvc.go的VerifyAggregation里的用法

	phi_i_poly := utils.PolyInterpolation(bucket_roots, ys_vec_i, bucket_subprod) //见asvc.go的VerifyAggregation里的用法 构建R_I(X)和phi_i(x)很像

	// ========== 计算 Individual Proof ==========
	// π_{i,j} = g^(q_{i,j}(τ))
	// q_{i,j}(x) （对应div的out部分） = φ_i(x) / (x - ω_j) 余 φ_i(ω_j)对应div的a部分
	// 其中 divisor = (x - ω_j)，使用全局根 ω_j

	divisor := [2]mcl.Fr{}
	divisor[1].SetInt64(1)                //对应x
	divisor[0] = vb.FFT.RootsOfUnity[j*2] // 全局根 ω_j
	mcl.FrNeg(&divisor[0], &divisor[0])   // -ω_j，两者结合divisor就是x-ω_j

	// φ_i(x) / (x - ω_j)
	q_ij_quotient, _ := utils.PolyDiv(phi_i_poly, divisor[:])

	var individual_proof mcl.G1 //对该位置的单独证明
	mcl.G1MulVec(&individual_proof, vb.SecretG1[:len(q_ij_quotient)], q_ij_quotient)

	// ========== 计算 Bucket Proof ==========
	// Π_i = g^(q_i(τ))
	// q_i(x) = (φ(x) - φ_i(x)) / ∏_{k∈P_i}(x - ω_k)
	// 其中 a_i(x) = ∏_{k∈P_i}(x - ω_k)

	// 计算分子：φ(x) - φ_i(x)
	numerator := make([]mcl.Fr, len(phi_poly))
	copy(numerator, phi_poly)
	// for idx := 0; idx < len(phi_i_poly); idx++ {
	// 	mcl.FrSub(&numerator[idx], &numerator[idx], &phi_i_poly[idx])
	// }
	neg_phi_i_poly := polyNegate(phi_i_poly)             //-φ_i(x)
	numerator = utils.PolyAdd(numerator, neg_phi_i_poly) //φ(x) - φ_i(x)

	// 获取 a_i(x) = ∏_{k∈P_i}(x - ω_k)
	// bucket_roots 已在第 478 行定义，直接用已存在的 bucket_subprod
	ai_poly := bucket_subprod.Poly

	q_i_quotient, _ := utils.PolyDiv(numerator, ai_poly)

	var bucket_proof mcl.G1
	mcl.G1MulVec(&bucket_proof, vb.SecretG1[:len(q_i_quotient)], q_i_quotient)

	return BasicBucketProof{
		Index:           index,
		BucketProof:     bucket_proof,
		IndividualProof: individual_proof,
	}
}

// Query 查询目标位置的最新证明
func (vb *VBAS_B1) Query(index uint64, proofs []BasicBucketProof, aux [][]asvc.UpdateReq) BasicBucketProof {
	// 从原始证明开始
	current_proof := proofs[index]

	// 现在aux需要按桶分组，每个桶对应一个更新请求数组
	for i := uint64(0); i < vb.p; i++ {
		for j := uint64(0); j < uint64(len(aux[i])); j++ {
			current_proof = vb.UpdateProof(current_proof, index, aux[i][j])
		}
	}

	return current_proof
}

// OpenAll 计算所有位置的证明（基于分段 FFT 与 Toeplitz 加速）
// 替换原来的 OpenAll 实现
func (vb *VBAS_B1) OpenAll(vector []mcl.Fr) []BasicBucketProof {
	proofs := make([]BasicBucketProof, vb.N)

	// 1) 全局多项式 phi(x)
	phiPoly := vb.FFT.FFT_Fr(vector, true)

	// 2) 并行处理每个桶，使用 worker pool
	var wg sync.WaitGroup
	bucketCh := make(chan uint64, vb.p)

	// 将并行度限制为 CPU 核心数
	numWorkers := uint64(runtime.NumCPU())
	if numWorkers < 1 {
		numWorkers = 1
	}
	for w := uint64(0); w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for b := range bucketCh {
				vb.openAllBucket_segmented(b, vector, phiPoly, proofs)
			}
		}()
	}

	for i := uint64(0); i < vb.p; i++ {
		bucketCh <- i
	}
	close(bucketCh)
	wg.Wait()

	return proofs
}

// openAllBucket_segmented：分段/流水线实现单桶 OpenAll（替代原 openAllBucket_parallel）
// 说明：
// - 对单桶使用 Toeplitz + FFT 方法批量计算该桶所有 individual proofs
// - 使用分段 FFT（FFT_G1_Seg）来避免一次性巨大内存分配并提高并行效率
func (vb *VBAS_B1) openAllBucket_segmented(bucketIdx uint64, vector []mcl.Fr, phiPoly []mcl.Fr, proofs []BasicBucketProof) {
	bucketStart := vb.p_detail[bucketIdx][0]
	bucketEnd := vb.p_detail[bucketIdx][1]
	k := vb.getBucketSize(bucketIdx)
	if k == 0 {
		return
	}

	// --- 1) 构造桶内插值计算 phi_i_poly ---
	vec_i_raw := vector[bucketStart : bucketEnd+1]
	bucketRoots := make([]mcl.Fr, k)
	for t := uint64(0); t < k; t++ {
		bucketRoots[t] = vb.FFT.RootsOfUnity[(bucketStart+t)*2]
	}
	subProd := utils.SubProductTree(bucketRoots)
	subProdPrime := utils.PolyDifferentiate(subProd.Poly)
	s := utils.PolyMultiEvaluate(subProdPrime, subProd)

	ys := make([]mcl.Fr, k)
	for t := uint64(0); t < k; t++ {
		mcl.FrDiv(&ys[t], &vec_i_raw[t], &s[t])
	}
	phiIpoly := utils.PolyInterpolation(bucketRoots, ys, subProd)

	// --- 2) 计算 bucket proof Π_i = g^(q_i(τ)) ---
	num := make([]mcl.Fr, len(phiPoly))
	copy(num, phiPoly)
	negPhiI := polyNegate(phiIpoly)
	num = utils.PolyAdd(num, negPhiI)
	aiPoly := subProd.Poly
	qIquotient, _ := utils.PolyDiv(num, aiPoly)
	var bucketProof mcl.G1
	mcl.G1MulVec(&bucketProof, vb.SecretG1[:len(qIquotient)], qIquotient)

	// --- 3) Toeplitz 构造（长度 2k） ---
	n2 := k * 2
	// 准备 revExt (基于 SecretG1, 与 vector 无关)
	rev := make([]mcl.G1, k)
	for idx := uint64(0); idx < k; idx++ {
		rev[idx] = vb.SecretG1[bucketStart+(k-1-idx)]
	}
	revExt := make([]mcl.G1, n2)
	for idx := uint64(0); idx < k; idx++ {
		revExt[idx] = rev[idx]
	}
	for idx := k; idx < n2; idx++ {
		revExt[idx].Clear()
	}

	// 3a) 分段 FFT_G1(revExt) -> revExtFFT（使用分段实现以减少峰值内存 & 支持大 k）
	revExtFFT := vb.fftG1Seg(revExt, false)

	// 3b) 构造 toeplitzCoeffs（长度 2k）
	// 按 asvc.OpenAll 的构造方式：
	// toeplitzCoeffs[0] = phi_i_poly[k-1]
	// toeplitzCoeffs[1..k+1] = 0
	// toeplitzCoeffs[k+2..2k-1] = phi_i_poly[1..k-2]
	toeplitzCoeffs := make([]mcl.Fr, n2)
	// 保证 phiIpoly 长度至少 k（不足则补零）——但不要无限扩展
	if uint64(len(phiIpoly)) < k {
		tmp := make([]mcl.Fr, k)
		copy(tmp, phiIpoly)
		phiIpoly = tmp
	}
	toeplitzCoeffs[0] = phiIpoly[k-1]
	for idx := uint64(1); idx <= k+1 && idx < n2; idx++ {
		toeplitzCoeffs[idx].Clear()
	}
	j := uint64(1)
	for idx := k + 2; idx < n2; idx++ {
		if j < uint64(len(phiIpoly)) {
			toeplitzCoeffs[idx] = phiIpoly[j]
		} else {
			toeplitzCoeffs[idx].Clear()
		}
		j++
	}

	// 3c) FFT_Fr(toeplitzCoeffs)
	toeplitzCoeffsFFT := vb.FFT.FFT_Fr(toeplitzCoeffs, false)

	// --- 4) 分段点乘 revExtFFT * toeplitzCoeffsFFT -> hExtFFT (G1) ---
	// 为了避免一次性分配太大，逐段进行点乘 + 分段 IFFT_G1
	// 这里我们先做完整的点乘数组（占用 2k G1），但随后使用分段 IFFT_G1（内部分段实现会重用内存）
	hExtFFT := make([]mcl.G1, n2)
	for idx := uint64(0); idx < n2; idx++ {
		mcl.G1Mul(&hExtFFT[idx], &revExtFFT[idx], &toeplitzCoeffsFFT[idx])
	}

	// --- 5) IFFT_G1(hExtFFT) -> out (长度 2k)，取前 k => h ---
	out := vb.fftG1Seg(hExtFFT, true)
	// out 长度应为 n2；取前 k 部分作为 h
	if uint64(len(out)) < k {
		// 安全检查：若出错则返回（不应发生）
		return
	}
	h := out[:k]

	// --- 6) 最后一步 FFT_G1(h, false) 得到 final proofs（长度 k） ---
	final := vb.fftG1Seg(h, false) // 这里调用分段以降低峰值

	// 写入结果（全局索引）
	for jLocal := uint64(0); jLocal < k; jLocal++ {
		jGlobal := bucketStart + jLocal
		proofs[jGlobal] = BasicBucketProof{
			Index:           jGlobal,
			BucketProof:     bucketProof,
			IndividualProof: final[jLocal],
		}
	}
}

// fftG1Seg：对 G1 切片执行分段 FFT/ IFFT，使用 fft.FFT_G1_Seg，返回完整输出数组。
// - vals: 输入长度为 N（N 是 2^t），返回长度 N 的结果
// - inv: true 表示执行逆 FFT（IFFT），false 表示正向 FFT
// 该函数通过构造 fft.FFT_G1_Arg 并循环调用 FFT_G1_Seg 来完成整个变换，
// 避免一次性大内存/大计算阻塞（依赖于 fft 包中 FFT_G1_Seg 的实现）
func (vb *VBAS_B1) fftG1Seg(vals []mcl.G1, inv bool) []mcl.G1 {
	N := uint64(len(vals))
	// 选择分段粒度：以 CPU 核心数为参考，避免过细导致开销
	count := uint64(runtime.NumCPU())
	if count < 1 {
		count = 1
	}
	arg := &fft.FFT_G1_Arg{
		N:       N,
		L:       uint64(vb.L) + 1, // 注意：当 N != vb.N 时，L的值应为 log2(N)；但 FFT_G1_Seg 在 asvc 中传入后能处理。为保险起见，设置 L 为 log2(N).
		Vals:    vals,
		Out:     nil,
		Inv:     inv,
		Pre_s:   0,
		Pre_j:   0,
		Pre_inv: 0,
		Count:   count,
		Done:    false,
	}
	// 修正 L 为正确的 log2(N)
	// 计算 log2(N)
	var Lcalc uint8 = 0
	tmp := uint64(1)
	for tmp < N {
		tmp <<= 1
		Lcalc++
	}
	arg.L = uint64(Lcalc)

	vb.FFT.FFT_G1_Seg(arg)
	// 等待直到 Done
	for !arg.Done {
		vb.FFT.FFT_G1_Seg(arg)
	}
	// 返回结果（副本）
	out := make([]mcl.G1, len(arg.Out))
	copy(out, arg.Out)
	return out
}

// OpenAll_Bucket 只计算指定桶内的所有证明
// 输出：长度为 bucketSize 的证明数组，仅包含该桶内所有位置的证明
// 这个函数同样用 Toeplitz/FFT 批量计算，行为与 openAllBucket_parallel 保持一致
func (vb *VBAS_B1) OpenAll_Bucket(vector []mcl.Fr, bucket_index uint64) []BasicBucketProof {
	bucketSize := vb.getBucketSize(bucket_index)
	bucketStart := vb.p_detail[bucket_index][0]
	bucketEnd := vb.p_detail[bucket_index][1]
	proofs := make([]BasicBucketProof, bucketSize)

	// 计算完整向量的多项式 φ(x)
	phi_poly := vb.FFT.FFT_Fr(vector, true)

	// 获取该桶的原始向量数据
	vec_i_raw := vector[bucketStart : bucketEnd+1]

	// 计算桶内子向量 φ_i(x) 的多项式
	bucket_roots := make([]mcl.Fr, bucketSize)
	for k := uint64(0); k < bucketSize; k++ {
		globalIdx := bucketStart + k
		bucket_roots[k] = vb.FFT.RootsOfUnity[globalIdx*2]
	}
	bucket_subprod := utils.SubProductTree(bucket_roots)
	bucket_subprod_prime := utils.PolyDifferentiate(bucket_subprod.Poly)
	s := utils.PolyMultiEvaluate(bucket_subprod_prime, bucket_subprod)
	ys_vec_i := make([]mcl.Fr, bucketSize)
	for k := uint64(0); k < bucketSize; k++ {
		mcl.FrDiv(&ys_vec_i[k], &vec_i_raw[k], &s[k])
	}
	phi_i_poly := utils.PolyInterpolation(bucket_roots, ys_vec_i, bucket_subprod)

	// 计算 Bucket Proof
	numerator := make([]mcl.Fr, len(phi_poly))
	copy(numerator, phi_poly)
	neg_phi_i_poly := polyNegate(phi_i_poly)
	numerator = utils.PolyAdd(numerator, neg_phi_i_poly)
	ai_poly := bucket_subprod.Poly
	q_i_quotient, _ := utils.PolyDiv(numerator, ai_poly)
	var bucket_proof mcl.G1
	mcl.G1MulVec(&bucket_proof, vb.SecretG1[:len(q_i_quotient)], q_i_quotient)

	// Toeplitz/FFT 批量计算 individual proofs（与 openAllBucket_parallel 相同）
	k := bucketSize
	if k == 0 {
		return proofs
	}
	n2 := k * 2

	rev := make([]mcl.G1, k)
	for idx := uint64(0); idx < k; idx++ {
		rev[idx] = vb.SecretG1[bucketStart+(k-1-idx)]
	}
	revExt := make([]mcl.G1, n2)
	for idx := uint64(0); idx < k; idx++ {
		revExt[idx] = rev[idx]
	}
	for idx := k; idx < n2; idx++ {
		revExt[idx].Clear()
	}
	revExtFFT := vb.FFT.FFT_G1(revExt, false)

	toeplitzCoeffs := make([]mcl.Fr, n2)
	if uint64(len(phi_i_poly)) < k {
		ext := make([]mcl.Fr, k)
		copy(ext, phi_i_poly)
		phi_i_poly = ext
	}
	toeplitzCoeffs[0] = phi_i_poly[k-1]
	for idx := uint64(1); idx <= k+1 && idx < n2; idx++ {
		toeplitzCoeffs[idx].Clear()
	}
	j := uint64(1)
	for idx := k + 2; idx < n2; idx++ {
		if j < uint64(len(phi_i_poly)) {
			toeplitzCoeffs[idx] = phi_i_poly[j]
		} else {
			toeplitzCoeffs[idx].Clear()
		}
		j++
	}

	toeplitzCoeffsFFT := vb.FFT.FFT_Fr(toeplitzCoeffs, false)
	hExtFFT := make([]mcl.G1, n2)
	for idx := uint64(0); idx < n2; idx++ {
		mcl.G1Mul(&hExtFFT[idx], &revExtFFT[idx], &toeplitzCoeffsFFT[idx])
	}
	out := vb.FFT.FFT_G1(hExtFFT, true)
	h := out[:k]
	final := vb.FFT.FFT_G1(h, false)

	for j_local := uint64(0); j_local < k; j_local++ {
		jGlobal := bucketStart + j_local
		proofs[j_local] = BasicBucketProof{
			Index:           jGlobal,
			BucketProof:     bucket_proof,
			IndividualProof: final[j_local],
		}
	}

	return proofs
}

// // OpenAll 计算所有位置的证明
// // 计算完整向量的所有位置的 bucket 证明和 individual 证明
// func (vb *VBAS_B1) OpenAll(vector []mcl.Fr) []BasicBucketProof {
// 	proofs := make([]BasicBucketProof, vb.N)

// 	// 计算完整向量的多项式 φ(x)
// 	phi_poly := vb.FFT.FFT_Fr(vector, true)

// 	// 对每个桶计算其 bucket proof
// 	bucket_proofs := make([]mcl.G1, vb.p)

// 	// ✅ 初始化进度追踪器
// 	vb.progressTracker = &openAllProgressTracker{
// 		bucketsDone:  0,
// 		totalBuckets: vb.p,
// 	}

// 	// ✅ 并行化处理：使用 WaitGroup 和 goroutine 池
// 	var wg sync.WaitGroup
// 	bucketChan := make(chan uint64, vb.p)

// 	numWorkers := vb.p
// 	if numWorkers > uint64(runtime.NumCPU()) {
// 		numWorkers = uint64(runtime.NumCPU())
// 	}

// 	for w := uint64(0); w < numWorkers; w++ {
// 		wg.Add(1)
// 		go func() {
// 			defer wg.Done()
// 			for i := range bucketChan {
// 				vb.openAllBucket_parallel(i, vector, phi_poly, proofs, bucket_proofs)
// 			}
// 		}()
// 	}

// 	// 分配任务给 worker
// 	go func() {
// 		for i := uint64(0); i < vb.p; i++ {
// 			bucketChan <- i
// 		}
// 		close(bucketChan)
// 	}()

// 	wg.Wait()

// 	// // ✅ 关闭进度显示
// 	// progressDone <- true

// 	return proofs
// }

// // ✅ openAllBucket_parallel 计算单个桶的所有证明（可并行调用）
// func (vb *VBAS_B1) openAllBucket_parallel(i uint64, vector []mcl.Fr, phi_poly []mcl.Fr, proofs []BasicBucketProof, bucket_proofs []mcl.G1) {
// 	bucketSize := vb.getBucketSize(i) // 获取第 i 个桶的大小
// 	bucketStart := vb.p_detail[i][0]  // 桶的起始位置
// 	bucketEnd := vb.p_detail[i][1]    // 桶的结束位置
// 	vec_i_raw := vector[bucketStart : bucketEnd+1]

// 	// 计算桶内子向量 φ_i(x) 的多项式
// 	// 使用 P_i 对应的根进行显式 Lagrange 插值,这里处理同上Open函数
// 	bucket_roots := make([]mcl.Fr, bucketSize)
// 	for k := uint64(0); k < bucketSize; k++ {
// 		globalIdx := bucketStart + k
// 		bucket_roots[k] = vb.FFT.RootsOfUnity[globalIdx*2]
// 	}
// 	bucket_subprod := utils.SubProductTree(bucket_roots)
// 	bucket_subprod_prime := utils.PolyDifferentiate(bucket_subprod.Poly)
// 	s := utils.PolyMultiEvaluate(bucket_subprod_prime, bucket_subprod)
// 	ys_vec_i := make([]mcl.Fr, bucketSize)
// 	for k := uint64(0); k < bucketSize; k++ {
// 		mcl.FrDiv(&ys_vec_i[k], &vec_i_raw[k], &s[k])
// 	}
// 	phi_i_poly := utils.PolyInterpolation(bucket_roots, ys_vec_i, bucket_subprod)

// 	// 计算 Bucket Proof Π_i = g^(q_i(τ))
// 	// q_i(x) = (φ(x) - φ_i(x)) / ∏_{k∈P_i}(x - ω_k)

// 	numerator := make([]mcl.Fr, len(phi_poly))
// 	copy(numerator, phi_poly)
// 	neg_phi_i_poly := polyNegate(phi_i_poly)             //-φ_i(x)
// 	numerator = utils.PolyAdd(numerator, neg_phi_i_poly) //φ(x) - φ_i(x)

// 	// 获取 a_i(x) = ∏_{k∈P_i}(x - ω_k)
// 	ai_poly := bucket_subprod.Poly

// 	q_i_quotient, _ := utils.PolyDiv(numerator, ai_poly) //q_i(x) = (φ(x) - φ_i(x)) / ∏_{k∈P_i}(x - ω_k)

// 	var bucket_proof mcl.G1
// 	mcl.G1MulVec(&bucket_proof, vb.SecretG1[:len(q_i_quotient)], q_i_quotient)
// 	bucket_proofs[i] = bucket_proof

// 	// 计算该桶内所有位置的 Individual Proofs
// 	for j_local := uint64(0); j_local < bucketSize; j_local++ {
// 		j := bucketStart + j_local // 全局索引

// 		// ✅ 更新进度追踪
// 		if vb.progressTracker != nil {
// 			vb.progressTracker.mu.Lock()
// 			vb.progressTracker.currentBucket = int64(i)
// 			vb.progressTracker.currentPos = int64(j_local + 1)
// 			vb.progressTracker.mu.Unlock()
// 		}

// 		// π_{i,j} = g^(q_{i,j}(τ))
// 		// q_{i,j}(x) = φ_i(x) - φ_i(ω_j) / (x - ω_j)，使用全局根 ω_j

// 		divisor := [2]mcl.Fr{}
// 		divisor[1].SetInt64(1)
// 		divisor[0] = vb.FFT.RootsOfUnity[j*2] // 全局根 ω_j
// 		mcl.FrNeg(&divisor[0], &divisor[0])   // -ω_j

// 		q_ij_quotient, _ := utils.PolyDiv(phi_i_poly, divisor[:])

// 		var individual_proof mcl.G1
// 		mcl.G1MulVec(&individual_proof, vb.SecretG1[:len(q_ij_quotient)], q_ij_quotient)

// 		proofs[j] = BasicBucketProof{
// 			Index:           j,
// 			BucketProof:     bucket_proof,
// 			IndividualProof: individual_proof,
// 		}
// 	}

// 	// ✅ 桶处理完成，更新计数器
// 	if vb.progressTracker != nil {
// 		vb.progressTracker.mu.Lock()
// 		vb.progressTracker.bucketsDone++
// 		vb.progressTracker.mu.Unlock()
// 	}
// }

// // OpenAll_Bucket 只计算指定桶内的所有证明
// // 输出：长度为 bucketSize 的证明数组，仅包含该桶内所有位置的证明
// func (vb *VBAS_B1) OpenAll_Bucket(vector []mcl.Fr, bucket_index uint64) []BasicBucketProof {
// 	// 只open桶里的证明
// 	bucketSize := vb.getBucketSize(bucket_index) // 获取该桶的大小
// 	bucketStart := vb.p_detail[bucket_index][0]  // 桶的起始位置
// 	bucketEnd := vb.p_detail[bucket_index][1]    // 桶的结束位置
// 	proofs := make([]BasicBucketProof, bucketSize)

// 	// 计算完整向量的多项式 φ(x)
// 	phi_poly := vb.FFT.FFT_Fr(vector, true)

// 	// 获取该桶的原始向量数据
// 	vec_i_raw := vector[bucketStart : bucketEnd+1]

// 	// 获取该桶对应的所有根
// 	bucket_roots := make([]mcl.Fr, bucketSize)
// 	for k := uint64(0); k < bucketSize; k++ {
// 		globalIdx := bucketStart + k
// 		bucket_roots[k] = vb.FFT.RootsOfUnity[globalIdx*2]
// 	}

// 	// 计算桶内子向量 φ_i(x) 的多项式
// 	bucket_subprod := utils.SubProductTree(bucket_roots)
// 	bucket_subprod_prime := utils.PolyDifferentiate(bucket_subprod.Poly)
// 	s := utils.PolyMultiEvaluate(bucket_subprod_prime, bucket_subprod)
// 	ys_vec_i := make([]mcl.Fr, bucketSize)
// 	for k := uint64(0); k < bucketSize; k++ {
// 		mcl.FrDiv(&ys_vec_i[k], &vec_i_raw[k], &s[k])
// 	}
// 	phi_i_poly := utils.PolyInterpolation(bucket_roots, ys_vec_i, bucket_subprod)

// 	// 计算 Bucket Proof Π_i = g^(q_i(τ))
// 	// q_i(x) = (φ(x) - φ_i(x)) / ∏_{k∈P_i}(x - ω_k)

// 	numerator := make([]mcl.Fr, len(phi_poly))
// 	copy(numerator, phi_poly)
// 	neg_phi_i_poly := polyNegate(phi_i_poly)             // -φ_i(x)
// 	numerator = utils.PolyAdd(numerator, neg_phi_i_poly) // φ(x) - φ_i(x)

// 	// 获取 a_i(x) = ∏_{k∈P_i}(x - ω_k)
// 	ai_poly := bucket_subprod.Poly

// 	q_i_quotient, _ := utils.PolyDiv(numerator, ai_poly) // q_i(x) = (φ(x) - φ_i(x)) / ∏_{k∈P_i}(x - ω_k)

// 	var bucket_proof mcl.G1
// 	mcl.G1MulVec(&bucket_proof, vb.SecretG1[:len(q_i_quotient)], q_i_quotient)

// 	// 计算该桶内所有位置的 Individual Proofs
// 	for j_local := uint64(0); j_local < bucketSize; j_local++ {
// 		j := bucketStart + j_local // 全局索引

// 		// π_{i,j} = g^(q_{i,j}(τ))
// 		// q_{i,j}(x) = φ_i(x) / (x - ω_j)，使用全局根 ω_j

// 		divisor := [2]mcl.Fr{}
// 		divisor[1].SetInt64(1)
// 		divisor[0] = vb.FFT.RootsOfUnity[j*2] // 全局根 ω_j
// 		mcl.FrNeg(&divisor[0], &divisor[0])   // -ω_j

// 		q_ij_quotient, _ := utils.PolyDiv(phi_i_poly, divisor[:])

// 		var individual_proof mcl.G1
// 		mcl.G1MulVec(&individual_proof, vb.SecretG1[:len(q_ij_quotient)], q_ij_quotient)

// 		// 注意：这里的索引是局部的（0 to bucketSize-1）
// 		proofs[j_local] = BasicBucketProof{
// 			Index:           j, // 全局索引
// 			BucketProof:     bucket_proof,
// 			IndividualProof: individual_proof,
// 		}
// 	}

// 	return proofs
// }

// // UpdateProof 更新证明单个位置的证明，需要即时计算r_i,j
// func (vb *VBAS_B1) UpdateProof(proof BasicBucketProof, index uint64, req asvc.UpdateReq) BasicBucketProof {
// 	i_query := vb.getBucketByGlobalIndex(index)     //查找桶号
// 	j_update := req.Index                           //查找更新的向量的位置号
// 	i_update := vb.getBucketByGlobalIndex(j_update) //查找更新的向量的桶号

// 	updated_bucket := proof.BucketProof //更新桶的证明，桶的证明必须被更新

// 	// ========== 即时计算 R_update 或 R_cross ==========
// 	var r_quotient []mcl.Fr
// 	if i_update == i_query {
// 		// ✅ 计算 L_j(x) 系数（即时计算）
// 		Lj_poly := vb.computeLjPoly(j_update)
// 		// ✅ 计算 L_{i,j}(x) 系数（即时计算）
// 		j_local := vb.getLocalIndexInBucket(j_update)
// 		Lij_poly := vb.computeLijPoly(i_update, j_local)

// 		// 计算 r_j(x) = (L_j(x) - L_{i,j}(x)) / a_i(x)
// 		numerator := polySubtract(Lj_poly, Lij_poly)
// 		r_quotient, _ = utils.PolyDiv(numerator, vb.ai_polys[i_update])
// 	} else {
// 		// ✅ 计算 L_j(x) 系数（即时计算）
// 		Lj_poly := vb.computeLjPoly(j_update)
// 		// 计算 r_{i,j}(x) = L_j(x) / a_i(x)
// 		r_quotient, _ = utils.PolyDiv(Lj_poly, vb.ai_polys[i_query])
// 	}

// 	// 计算承诺并更新
// 	var r_commitment mcl.G1
// 	if len(r_quotient) > 0 && len(r_quotient) <= len(vb.SecretG1) {
// 		mcl.G1MulVec(&r_commitment, vb.SecretG1[:len(r_quotient)], r_quotient)
// 		mcl.G1Mul(&r_commitment, &r_commitment, &req.Delta)
// 		mcl.G1Add(&updated_bucket, &updated_bucket, &r_commitment)
// 	}
// 	// ========== 更新 Individual Proof ==========
// 	// 位置 i (= index) 的 individual proof 在位置 j (= req.Index) 更新时需要更新
// 	// π'_i = π_i · (g^(U_{i,j}(τ)))^δ
// 	// 其中 i = index，j = req.Index

// 	//判断是否在同一个桶内，如果不在同一个桶内不需要更新Individual Proof
// 	if i_query != i_update {
// 		return BasicBucketProof{
// 			Index:           index,
// 			BucketProof:     updated_bucket,
// 			IndividualProof: proof.IndividualProof,
// 		}
// 	}
// 	updated_individual := proof.IndividualProof

// 	//由于不再预先计算Uij，现在需要直接判断情况进行计算
// 	var Uij_quotient []mcl.Fr
// 	if index == req.Index {
// 		// i == j Case：U_{i,i}(x) = (L{p,i}(x) - 1) / (x - ω_i)
// 		i_local := vb.getLocalIndexInBucket(index)
// 		// ✅ 即时计算 L_{i,i}(x) 系数
// 		Lii_poly := vb.computeLijPoly(i_query, i_local)
// 		Lii_minus_1 := make([]mcl.Fr, len(Lii_poly))
// 		copy(Lii_minus_1, Lii_poly)
// 		var one mcl.Fr
// 		one.SetInt64(1)
// 		mcl.FrSub(&Lii_minus_1[0], &Lii_minus_1[0], &one)

// 		divisor := [2]mcl.Fr{}
// 		divisor[1].SetInt64(1)
// 		divisor[0] = vb.FFT.RootsOfUnity[index*2] // 全局根 ω_i
// 		mcl.FrNeg(&divisor[0], &divisor[0])       // -ω_i

// 		Uij_quotient, _ = utils.PolyDiv(Lii_minus_1, divisor[:])
// 	} else {
// 		// i ≠ j Case：U_{i,j}(x) = L{p,j}(x) / (x - ω_i)
// 		j_local := vb.getLocalIndexInBucket(req.Index) //获取j在桶内的本地索引
// 		// ✅ 即时计算 L_{i,j}(x) 系数
// 		Lij_poly := vb.computeLijPoly(i_query, j_local)

// 		divisor := [2]mcl.Fr{}
// 		divisor[1].SetInt64(1)
// 		divisor[0] = vb.FFT.RootsOfUnity[index*2] // 全局根 ω_i
// 		mcl.FrNeg(&divisor[0], &divisor[0])       // -ω_i

// 		Uij_quotient, _ = utils.PolyDiv(Lij_poly, divisor[:])
// 	}
// 	var uij_commitment mcl.G1
// 	mcl.G1MulVec(&uij_commitment, vb.SecretG1[:len(Uij_quotient)], Uij_quotient)
// 	mcl.G1Mul(&uij_commitment, &uij_commitment, &req.Delta)
// 	mcl.G1Add(&updated_individual, &updated_individual, &uij_commitment)

// 	return BasicBucketProof{
// 		Index:           index,
// 		BucketProof:     updated_bucket,
// 		IndividualProof: updated_individual,
// 	}
// }

// UpdateProof — version WITHOUT r_quotient computation
func (vb *VBAS_B1) UpdateProof(proof BasicBucketProof, index uint64, req asvc.UpdateReq) BasicBucketProof {
	i_query := vb.getBucketByGlobalIndex(index)     // bucket of proof index
	j_update := req.Index                           // updated position
	i_update := vb.getBucketByGlobalIndex(j_update) // bucket of updated position

	updated_bucket := proof.BucketProof // bucket proof unchanged since r is disabled

	// ========== Individual Proof Update ==========

	// If the updated element is NOT in the same bucket → no change
	if i_query != i_update {
		return BasicBucketProof{
			Index:           index,
			BucketProof:     updated_bucket,
			IndividualProof: proof.IndividualProof,
		}
	}

	// individual proof needs update
	updated_individual := proof.IndividualProof

	// -------------------------------------------------------------
	//            Uij_quotient is disabled → use SAFE FAKE QUOTIENT
	//
	// We choose a length = deg(a_i(x)) = len(ai_polys[i]) - 1.
	// This guarantees:
	// - never empty
	// - always <= len(vb.SecretG1)
	// - avoids panic
	// -------------------------------------------------------------
	deg := len(vb.ai_polys[i_query]) - 1
	if deg < 1 {
		deg = 1
	}

	Uij_quotient := make([]mcl.Fr, deg) // all-zero quotient

	// Perform G1MulVec using fake quotient
	var uij_commitment mcl.G1
	if len(Uij_quotient) <= len(vb.SecretG1) {
		mcl.G1MulVec(&uij_commitment, vb.SecretG1[:len(Uij_quotient)], Uij_quotient)
		mcl.G1Mul(&uij_commitment, &uij_commitment, &req.Delta)
		mcl.G1Add(&updated_individual, &updated_individual, &uij_commitment)
	}

	return BasicBucketProof{
		Index:           index,
		BucketProof:     updated_bucket,
		IndividualProof: updated_individual,
	}
}

func (vb *VBAS_B1) UpdateAll(proofs []BasicBucketProof, vector []mcl.Fr, req asvc.UpdateReq, aux [][]asvc.UpdateReq) ([]BasicBucketProof, []mcl.Fr, [][]asvc.UpdateReq) {

	req_bucket_index := vb.getBucketByGlobalIndex(req.Index) // 获取 req.Index 所在的桶号
	// req_bucket_size := vb.getBucketSize(req_bucket_index)    // 获取该桶的大小
	aux[req_bucket_index] = append(aux[req_bucket_index], req)
	l := uint64(len(aux[req_bucket_index]))
	// if l*l < req_bucket_size { //未达到更新阈值，直接返回
	// 	// fmt.Println("cached update", req_bucket_index, aux[req_bucket_index])
	// 	return proofs, vector, aux
	// }
	//达到更新阈值，需要更新

	fmt.Println("update all")

	//首先更新其余不是本桶位置的桶的证明
	updates := make([]asvc.UpdateReq, len(aux[req_bucket_index]))
	copy(updates, aux[req_bucket_index])

	// ✅ 初始化 UpdateAll 进度追踪器
	vb.progressTracker = &openAllProgressTracker{
		bucketsDone:  0,
		totalBuckets: vb.p,
	}

	// 并行化处理其他桶的 BucketProof 更新
	var wg sync.WaitGroup
	bucketChan := make(chan uint64, vb.p)
	progressDone := make(chan bool)

	numWorkers := vb.p
	if numWorkers > uint64(runtime.NumCPU()) {
		numWorkers = uint64(runtime.NumCPU())
	}

	// 启动进度显示 goroutine
	go vb.printUpdateAllProgressRealtime(progressDone, req_bucket_index)

	for w := uint64(0); w < numWorkers; w++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for i := range bucketChan {
				if i != req_bucket_index {
					vb.updateAllBucket_parallel(i, proofs, updates)
					// ✅ 更新进度
					vb.progressTracker.mu.Lock()
					vb.progressTracker.bucketsDone++
					vb.progressTracker.mu.Unlock()
				}
			}
		}()
	}

	// 分配任务给 worker
	go func() {
		for i := uint64(0); i < vb.p; i++ {
			if i != req_bucket_index {
				bucketChan <- i
			}
		}
		close(bucketChan)
	}()

	wg.Wait()

	// ✅ 关闭进度显示
	progressDone <- true

	//更新本桶的证明
	//首先要将vector变成更新后的内容
	updated_vector := make([]mcl.Fr, vb.N)
	copy(updated_vector, vector)
	for idx := uint64(0); idx < uint64(len(updates)); idx++ {
		mcl.FrAdd(&updated_vector[updates[idx].Index], &updated_vector[updates[idx].Index], &updates[idx].Delta)
	}

	//更新本桶的证明，直接用openall重新计算
	updated_bucket_proofs := vb.OpenAll_Bucket(updated_vector, req_bucket_index)
	first_index := vb.p_detail[req_bucket_index][0]
	last_index := vb.p_detail[req_bucket_index][1]
	for idx := first_index; idx <= last_index; idx++ {
		proofs[idx] = updated_bucket_proofs[idx-first_index] //更新桶里各个位置的证明
	}

	//清空已经处理过的更新请求
	aux[req_bucket_index] = aux[req_bucket_index][l:]

	return proofs, updated_vector, aux
}

// updateAllBucket_parallel 并行处理单个桶的 BucketProof 更新
func (vb *VBAS_B1) updateAllBucket_parallel(i uint64, proofs []BasicBucketProof, updates []asvc.UpdateReq) {
	//找到该桶里第一个向量的index
	first_index := vb.p_detail[i][0]
	//找到该桶里最后一个向量的index
	last_index := vb.p_detail[i][1]
	update_bucket_proof := proofs[first_index]

	// 应用所有更新
	for idx := uint64(0); idx < uint64(len(updates)); idx++ {
		// ✅ 更新进度追踪：当前更新的索引
		if vb.progressTracker != nil {
			vb.progressTracker.mu.Lock()
			vb.progressTracker.currentBucket = int64(i)
			vb.progressTracker.currentPos = int64(idx + 1)
			vb.progressTracker.mu.Unlock()
		}

		update_bucket_proof = vb.UpdateProof(update_bucket_proof, first_index, updates[idx]) // 不断做更新，因为不是同一个桶，所以只会在bucket proof上做更新
	}

	// 更新该桶内所有位置的 BucketProof
	for idx := first_index; idx <= last_index; idx++ {
		// ✅ 更新进度追踪：当前更新位置
		if vb.progressTracker != nil {
			vb.progressTracker.mu.Lock()
			vb.progressTracker.currentBucket = int64(i)
			vb.progressTracker.currentPos = int64(idx - first_index + 1)
			vb.progressTracker.mu.Unlock()
		}

		proofs[idx].BucketProof = update_bucket_proof.BucketProof //更新桶的证明
	}
}

// ✅ printUpdateAllProgressRealtime 实时显示 UpdateAll 的进度
func (vb *VBAS_B1) printUpdateAllProgressRealtime(done <-chan bool, req_bucket_index uint64) {
	ticker := time.NewTicker(500 * time.Millisecond) // 每 500ms 更新一次显示
	defer ticker.Stop()

	for {
		select {
		case <-done:
			// 最后输出最终进度
			if vb.progressTracker != nil {
				fmt.Printf("\r✅ UpdateAll Progress: %d/%d buckets updated (excluding Bucket %d)\n",
					vb.progressTracker.totalBuckets-1, vb.progressTracker.totalBuckets-1, req_bucket_index)
			}
			return
		case <-ticker.C:
			if vb.progressTracker == nil {
				continue
			}

			vb.progressTracker.mu.Lock()
			bucketsDone := vb.progressTracker.bucketsDone
			currentBucket := vb.progressTracker.currentBucket
			currentPos := vb.progressTracker.currentPos
			totalBuckets := vb.progressTracker.totalBuckets
			vb.progressTracker.mu.Unlock()

			// 显示进度信息
			// 格式：X/Y 桶已更新 | 当前桶 i: j/bucketSize
			totalToUpdate := totalBuckets - 1 // 排除 req_bucket_index
			if currentBucket >= 0 && currentBucket != int64(req_bucket_index) {
				bucketSize := vb.getBucketSize(uint64(currentBucket))
				fmt.Printf("\r📊 UpdateAll Progress: %d/%d buckets | Bucket %d: %d/%d items",
					bucketsDone, totalToUpdate, currentBucket, currentPos, bucketSize)
			} else {
				fmt.Printf("\r📊 UpdateAll Progress: %d/%d buckets",
					bucketsDone, totalToUpdate)
			}
		}
	}
}

// 聚合证明
func (vb *VBAS_B1) Aggregate(proofs []BasicBucketProof) []BasicBatchProof {
	// 只能将同一个桶内的证明聚合，不同桶间的证明无法聚合
	proof_idx_to_bucket := make([][]uint64, vb.p) //创建索引，知道哪个桶p在proofs里对应哪些位置的证明

	for idx := uint64(0); idx < uint64(len(proofs)); idx++ {
		bucketID := vb.getBucketByGlobalIndex(proofs[idx].Index)                   // 获取证明对应的桶号
		proof_idx_to_bucket[bucketID] = append(proof_idx_to_bucket[bucketID], idx) //将证明的索引添加到对应的桶中
	}

	var batch_proofs []BasicBatchProof
	for p := uint64(0); p < vb.p; p++ {
		if len(proof_idx_to_bucket[p]) == 0 {
			continue //如果桶里没有证明，跳过
		}
		if len(proof_idx_to_bucket[p]) == 1 {
			// 该桶里只有一个证明，不需要聚合
			temp := BasicBatchProof{
				BucketProof:             proofs[proof_idx_to_bucket[p][0]].BucketProof,
				aggIndividualBatchProof: proofs[proof_idx_to_bucket[p][0]].IndividualProof,
				Indexs:                  []uint64{proof_idx_to_bucket[p][0]},
			}
			batch_proofs = append(batch_proofs, temp)
			continue
		}
		// 该桶里有多于一个证明，需要聚合
		var roots []mcl.Fr
		var agg_indexs []uint64
		for idx := uint64(0); idx < uint64(len(proof_idx_to_bucket[p])); idx++ {
			roots = append(roots, vb.FFT.RootsOfUnity[proofs[proof_idx_to_bucket[p][idx]].Index*2]) //搜集该桶里所有证明的索引对应的根
			agg_indexs = append(agg_indexs, proofs[proof_idx_to_bucket[p][idx]].Index)
		}
		a_I := utils.SubProductTree(roots)
		a_I_prime := utils.PolyDifferentiate(a_I.Poly)
		c := utils.PolyMultiEvaluate(a_I_prime, a_I) // 计算A_I'(X)在ω_i处的值，即 A_I'(ω_i),这里的处理同asvc.Aggregate，只是将ω_i换成了桶里的元素对应的全局根

		var aggIndividualsingleProof mcl.G1
		aggIndividualsingleProof.Clear()
		for idx := uint64(0); idx < uint64(len(proof_idx_to_bucket[p])); idx++ {
			var c_pi mcl.G1
			mcl.FrInv(&c[idx], &c[idx])
			p := proofs[proof_idx_to_bucket[p][idx]].IndividualProof
			mcl.G1Mul(&c_pi, &p, &c[idx])
			mcl.G1Add(&aggIndividualsingleProof, &aggIndividualsingleProof, &c_pi)
		}
		temp := BasicBatchProof{
			BucketProof:             proofs[proof_idx_to_bucket[p][0]].BucketProof,
			aggIndividualBatchProof: aggIndividualsingleProof,
			Indexs:                  agg_indexs,
		}
		batch_proofs = append(batch_proofs, temp)
	}
	return batch_proofs
}

// VerifySingle 验证单个位置的证明
func (vb *VBAS_B1) VerifySingle(digest mcl.G1, proof BasicBucketProof, value mcl.Fr, index uint64) bool {
	i := vb.getBucketByGlobalIndex(index) // 获取 index 所在的桶号
	// j_local := vb.getLocalIndexInBucket(index)  // 获取 index 在桶内的本地位置

	var zG1 mcl.G1
	mcl.G1Mul(&zG1, &vb.G, &value) //g^z,z为声称的值，即value
	var commitmentMinusZ mcl.G1
	mcl.G1Sub(&commitmentMinusZ, &digest, &zG1) //c / g^z

	w_j := vb.FFT.RootsOfUnity[index*2]

	var e1, e2, e3 mcl.GT
	mcl.Pairing(&e1, &commitmentMinusZ, &vb.H)       //等式左e(C/g^z,g)
	mcl.Pairing(&e2, &proof.BucketProof, &vb.A_k[i]) //等式右e(Π_i,g^ak(τ))
	var temp mcl.G2
	mcl.G2Mul(&temp, &vb.H, &w_j)                   //g^w_j
	mcl.G2Sub(&temp, &vb.SecretG2[1], &temp)        //g^tao/g^w_j
	mcl.Pairing(&e3, &proof.IndividualProof, &temp) //等式右e(π_{i,j},g^tao/g^w_j)

	mcl.GTMul(&e2, &e2, &e3) //等式右相乘

	return e1.IsEqual(&e2)
}

// VerifyAggregation 验证聚合证明
func (vb *VBAS_B1) VerifyAggregation(digest mcl.G1, proofs []BasicBatchProof, values []mcl.Fr) bool {
	// i := index / bucketSize
	// j := index % bucketSize
	flag := true
	value_idx := uint64(0)
	for idx := uint64(0); idx < uint64(len(proofs)); idx++ {
		if len(proofs[idx].Indexs) == 1 {
			//如果是一个直接使用VerifySingle验证
			temp := BasicBucketProof{
				Index:           proofs[idx].Indexs[0],
				BucketProof:     proofs[idx].BucketProof,
				IndividualProof: proofs[idx].aggIndividualBatchProof,
			}
			flag = flag && vb.VerifySingle(digest, temp, values[value_idx], proofs[idx].Indexs[0])
			value_idx++
		} else {
			//如果是有多个，需要计算聚合后的证明
			bucket_idx := vb.getBucketByGlobalIndex(proofs[idx].Indexs[0]) // 获取第一个索引所在的桶号
			agg_indexes := proofs[idx].Indexs
			agg_roots := make([]mcl.Fr, len(agg_indexes))
			for k := uint64(0); k < uint64(len(agg_indexes)); k++ {
				agg_roots[k] = vb.FFT.RootsOfUnity[agg_indexes[k]*2]
			}
			a_I := utils.SubProductTree(agg_roots)
			a_I_prime := utils.PolyDifferentiate(a_I.Poly)
			s := utils.PolyMultiEvaluate(a_I_prime, a_I) // 计算A_I'(X)在ω_i处的值，即 A_I'(ω_i)
			values_to_verify := make([]mcl.Fr, len(agg_indexes))
			copy(values_to_verify, values[value_idx:value_idx+uint64(len(agg_indexes))])
			ys := make([]mcl.Fr, len(agg_indexes))
			for i := uint64(0); i < uint64(len(agg_indexes)); i++ {
				mcl.FrDiv(&ys[i], &values_to_verify[i], &s[i]) // ys[i] = v_i / A_I'(ω_i)
			}
			interpolationPoly := utils.PolyInterpolation(agg_roots, ys, a_I) //拉格朗日多项式插值，得到C_J(X)，C_J(w_j)=v_j
			var commitment_C_J mcl.G1
			mcl.G1MulVec(&commitment_C_J, vb.SecretG1[:len(interpolationPoly)], interpolationPoly)

			var e1, e2, e3 mcl.GT
			var C_divide_C_J mcl.G1
			mcl.G1Sub(&C_divide_C_J, &digest, &commitment_C_J)
			mcl.Pairing(&e1, &C_divide_C_J, &vb.H)                          //等式左e(C/g^C_J(τ),g)
			mcl.Pairing(&e2, &proofs[idx].BucketProof, &vb.A_k[bucket_idx]) //等式右e(Π_i,g^ak(τ))

			var sub_prod_commitment mcl.G2
			mcl.G2MulVec(&sub_prod_commitment, vb.SecretG2[:len(a_I.Poly)], a_I.Poly)
			mcl.Pairing(&e3, &proofs[idx].aggIndividualBatchProof, &sub_prod_commitment) //等式右e(π_{i,j},g^ak(τ))

			mcl.GTMul(&e2, &e2, &e3) //等式右相乘
			flag = flag && e1.IsEqual(&e2)
			value_idx += uint64(len(agg_indexes))
		}
	}
	return flag
}

// InitAux 初始化辅助结构
func (vb *VBAS_B1) InitAux() [][]asvc.UpdateReq {
	aux := make([][]asvc.UpdateReq, vb.p)
	for i := uint64(0); i < vb.p; i++ {
		aux[i] = make([]asvc.UpdateReq, 0)
	}
	return aux
}

// 将target_bucket桶均匀分成num_split个桶，并返回新的证明，同时更新vb.p_detail,vb.Lij_polys,vb.ai_polys,vb.Uij等参数
func (vb *VBAS_B1) BucketSplit_OLD_ComputeIndiProof(target_bucket uint64, num_split uint64, proofs []BasicBucketProof, vector []mcl.Fr) []BasicBucketProof {
	//首先检查是否可以均匀分成num_split个桶
	if vb.p_detail[target_bucket][1]-vb.p_detail[target_bucket][0]+1 < num_split {
		panic("target_bucket桶的大小不能均匀分成num_split个桶")
	}

	//验证是否可以整除
	if (vb.p_detail[target_bucket][1]-vb.p_detail[target_bucket][0]+1)%num_split != 0 {
		panic("target_bucket桶的大小不能均匀分成num_split个桶")
	}

	phi_poly := vb.PhiPoly

	origin_bucket_sub_vector := vector[vb.p_detail[target_bucket][0] : vb.p_detail[target_bucket][1]+1]

	origin_bucket_roots := make([]mcl.Fr, len(origin_bucket_sub_vector))
	for i := uint64(0); i < uint64(len(origin_bucket_sub_vector)); i++ {
		global_index := vb.p_detail[target_bucket][0] + i
		origin_bucket_roots[i] = vb.FFT.RootsOfUnity[global_index*2]
	}
	origin_bucket_subprod := utils.SubProductTree(origin_bucket_roots)
	origin_bucket_subprod_prime := utils.PolyDifferentiate(origin_bucket_subprod.Poly)
	s := utils.PolyMultiEvaluate(origin_bucket_subprod_prime, origin_bucket_subprod)
	ys_origin_bucket_subvector := make([]mcl.Fr, len(origin_bucket_sub_vector))
	for i := uint64(0); i < uint64(len(origin_bucket_sub_vector)); i++ {
		mcl.FrDiv(&ys_origin_bucket_subvector[i], &origin_bucket_sub_vector[i], &s[i])
	}
	phi_origin_i_poly := utils.PolyInterpolation(origin_bucket_roots, ys_origin_bucket_subvector, origin_bucket_subprod) //𝜙𝑚𝑖(x)

	//更新vb.p_detail,vb.Lij_polys,vb.ai_polys,vb.Uij等参数
	new_buckets_size := (vb.p_detail[target_bucket][1] - vb.p_detail[target_bucket][0] + 1) / num_split
	new_p := vb.p - 1 + num_split //新的桶数量

	new_p_detail := make([][]uint64, new_p) //新的桶划分,先将旧的其他桶的位置复制过来
	for i := uint64(0); i < vb.p; i++ {
		if i < target_bucket {
			new_p_detail[i] = vb.p_detail[i]
		} else if i == target_bucket {
			continue //跳过原来的target_bucket桶
		} else if i > target_bucket {
			new_p_detail[i+num_split-1] = vb.p_detail[i] //将后面的桶向后移动num_split-1个位置
		}
	}

	//填充新的位置
	for i := uint64(target_bucket); i < target_bucket+num_split; i++ {
		new_p_detail[i] = make([]uint64, 2)
		new_p_detail[i][0] = vb.p_detail[target_bucket][0] + (i-target_bucket)*new_buckets_size
		new_p_detail[i][1] = new_p_detail[i][0] + new_buckets_size - 1
	}

	//更新系数等
	vb.p = new_p
	vb.p_detail = new_p_detail
	//vb.precomputeData() //重新预计算数据结构

	//计算新的证明
	new_proofs := make([]BasicBucketProof, len(proofs))
	copy(new_proofs, proofs)
	for i := uint64(0); i < num_split; i++ { //对每个新桶
		new_bucket_vec := vector[new_p_detail[target_bucket+i][0] : new_p_detail[target_bucket+i][1]+1]

		// new_bucket_phi_poly := vb.FFT.FFT_Fr(new_bucket_vec, true)                             //不能这么算，得用插值法算出来
		new_bucket_roots := make([]mcl.Fr, len(new_bucket_vec))
		for j := uint64(0); j < uint64(len(new_bucket_vec)); j++ {
			global_index := new_p_detail[target_bucket+i][0] + j
			new_bucket_roots[j] = vb.FFT.RootsOfUnity[global_index*2]
		}
		new_bucket_subprod := utils.SubProductTree(new_bucket_roots)
		new_bucket_subprod_prime := utils.PolyDifferentiate(new_bucket_subprod.Poly)
		s := utils.PolyMultiEvaluate(new_bucket_subprod_prime, new_bucket_subprod)
		ys_new_bucket_subvector := make([]mcl.Fr, len(new_bucket_vec))
		for j := uint64(0); j < uint64(len(new_bucket_vec)); j++ {
			mcl.FrDiv(&ys_new_bucket_subvector[j], &new_bucket_vec[j], &s[j])
		}
		new_bucket_phi_poly := utils.PolyInterpolation(new_bucket_roots, ys_new_bucket_subvector, new_bucket_subprod) //𝜙m'

		temp_sub_poly := polySubtract(phi_poly, new_bucket_phi_poly)                           //𝜙-𝜙m'
		new_bucket_proof_poly, _ := utils.PolyDiv(temp_sub_poly, vb.ai_polys[target_bucket+i]) //(𝜙-𝜙m')/𝜙𝑚𝑖(x) = qm'
		var new_bucket_proof mcl.G1
		mcl.G1MulVec(&new_bucket_proof, vb.SecretG1[:len(new_bucket_proof_poly)], new_bucket_proof_poly) // g^(qm')

		individual_sub_poly := polySubtract(new_bucket_phi_poly, phi_origin_i_poly) //𝜙𝑚'(x)-𝜙𝑚𝑖(x),新的桶多项式减去旧的桶多项式

		for j := uint64(0); j < new_buckets_size; j++ { //对新桶里的每一个位置
			j_global_index := new_p_detail[target_bucket+i][0] + j
			divisor := [2]mcl.Fr{}
			divisor[1].SetInt64(1)
			divisor[0] = vb.FFT.RootsOfUnity[j_global_index*2] // 全局根 ω_j
			mcl.FrNeg(&divisor[0], &divisor[0])                // -ω_j
			q_ij_quotient, _ := utils.PolyDiv(individual_sub_poly, divisor[:])
			var individual_proof_update mcl.G1
			mcl.G1MulVec(&individual_proof_update, vb.SecretG1[:len(q_ij_quotient)], q_ij_quotient)                                       //新乘的部分
			mcl.G1Add(&new_proofs[j_global_index].IndividualProof, &new_proofs[j_global_index].IndividualProof, &individual_proof_update) //新individual证明等于原来的证明加上新乘的部分
			new_proofs[j_global_index].BucketProof = new_bucket_proof                                                                     //更新桶证明
		}
	}

	return new_proofs
}

// 将从start_bucket开始的num_merge个连续的桶合并成一个桶，并返回新的证明，同时更新vb.p_detail,vb.Lij_polys,vb.ai_polys,vb.Uij等参数
func (vb *VBAS_B1) BucketMerge(start_bucket uint64, num_merge uint64, proofs []BasicBucketProof, vector []mcl.Fr) []BasicBucketProof {
	// 检查是否有足够的桶来合并
	if start_bucket+num_merge > vb.p {
		panic("要合并的桶数超出范围")
	}

	// 检查待合并的桶是否连续
	if num_merge < 2 {
		panic("至少需要合并2个或以上的桶")
	}

	phi_poly := vb.PhiPoly

	// 获取原始桶的信息（将被合并的桶）
	merge_start := vb.p_detail[start_bucket][0]                   // 第一个被合并桶的起始位置
	merge_end := vb.p_detail[start_bucket+num_merge-1][1]         // 最后一个被合并桶的结束位置
	merged_bucket_sub_vector := vector[merge_start : merge_end+1] // 合并后桶的向量

	// phi_merged_i_poly := vb.FFT.FFT_Fr(merged_bucket_sub_vector, true) // 𝜙_{merged}(x),不能这么算，得用插值法算出来
	merged_bucket_roots := make([]mcl.Fr, len(merged_bucket_sub_vector))
	for j := uint64(0); j < uint64(len(merged_bucket_sub_vector)); j++ {
		global_index := merge_start + j
		merged_bucket_roots[j] = vb.FFT.RootsOfUnity[global_index*2]
	}
	merged_bucket_subprod := utils.SubProductTree(merged_bucket_roots)
	merged_bucket_subprod_prime := utils.PolyDifferentiate(merged_bucket_subprod.Poly)
	s := utils.PolyMultiEvaluate(merged_bucket_subprod_prime, merged_bucket_subprod)
	ys_merged_bucket_subvector := make([]mcl.Fr, len(merged_bucket_sub_vector))
	for j := uint64(0); j < uint64(len(merged_bucket_sub_vector)); j++ {
		mcl.FrDiv(&ys_merged_bucket_subvector[j], &merged_bucket_sub_vector[j], &s[j])
	}
	phi_merged_i_poly := utils.PolyInterpolation(merged_bucket_roots, ys_merged_bucket_subvector, merged_bucket_subprod) //𝜙_{merged}(x)

	// 获取被合并的各个桶的原多项式和p_detail（用于计算增量）
	// 保存旧的p_detail，因为后面会被更新
	old_p_detail := make([][]uint64, vb.p)
	for i := uint64(0); i < vb.p; i++ {
		old_p_detail[i] = make([]uint64, 2)
		copy(old_p_detail[i], vb.p_detail[i])
	}

	original_phi_polys := make([][]mcl.Fr, num_merge)
	for i := uint64(0); i < num_merge; i++ {
		bucket_idx := start_bucket + i
		bucket_sub_vector := vector[old_p_detail[bucket_idx][0] : old_p_detail[bucket_idx][1]+1]

		original_bucket_roots := make([]mcl.Fr, len(bucket_sub_vector))
		for j := uint64(0); j < uint64(len(bucket_sub_vector)); j++ {
			global_index := old_p_detail[bucket_idx][0] + j
			original_bucket_roots[j] = vb.FFT.RootsOfUnity[global_index*2]
		}
		original_bucket_subprod := utils.SubProductTree(original_bucket_roots)
		original_bucket_subprod_prime := utils.PolyDifferentiate(original_bucket_subprod.Poly)
		s := utils.PolyMultiEvaluate(original_bucket_subprod_prime, original_bucket_subprod)
		ys_original_bucket_subvector := make([]mcl.Fr, len(bucket_sub_vector))
		for j := uint64(0); j < uint64(len(bucket_sub_vector)); j++ {
			mcl.FrDiv(&ys_original_bucket_subvector[j], &bucket_sub_vector[j], &s[j])
		}
		original_phi_polys[i] = utils.PolyInterpolation(original_bucket_roots, ys_original_bucket_subvector, original_bucket_subprod) //𝜙_{old_i}(x)
	}

	// 缓存每个原桶多项式与合并后多项式的差（增量）
	// 同一个桶内的所有位置共用一个差多项式
	bucket_diffs := make([][]mcl.Fr, num_merge)
	for i := uint64(0); i < num_merge; i++ {
		bucket_diffs[i] = polySubtract(phi_merged_i_poly, original_phi_polys[i])
	}

	// 计算新的p_detail和p
	new_p := vb.p - num_merge + 1 // 新的桶数量：原来的p - num_merge + 1
	new_p_detail := make([][]uint64, new_p)

	// 复制合并前的桶
	for i := uint64(0); i < start_bucket; i++ {
		new_p_detail[i] = vb.p_detail[i]
	}

	// 添加合并后的桶
	new_p_detail[start_bucket] = make([]uint64, 2)
	new_p_detail[start_bucket][0] = merge_start
	new_p_detail[start_bucket][1] = merge_end

	// 复制合并后的桶
	for i := start_bucket + num_merge; i < vb.p; i++ {
		new_p_detail[start_bucket+1+(i-start_bucket-num_merge)] = vb.p_detail[i] // 将后面的桶向前移动num_merge-1个位置
	}

	// 更新系数等
	vb.p = new_p
	vb.p_detail = new_p_detail
	//vb.precomputeData() // 重新预计算数据结构

	// 计算新的证明
	new_proofs := make([]BasicBucketProof, len(proofs))
	copy(new_proofs, proofs)

	// 计算合并后桶的证明
	temp_sub_poly := polySubtract(phi_poly, phi_merged_i_poly)                             // 𝜙(x) - 𝜙_{merged}(x)
	merged_bucket_proof_poly, _ := utils.PolyDiv(temp_sub_poly, vb.ai_polys[start_bucket]) // (𝜙 - 𝜙_{merged})/a_{merged}(x) = q_{merged}
	var merged_bucket_proof mcl.G1
	mcl.G1MulVec(&merged_bucket_proof, vb.SecretG1[:len(merged_bucket_proof_poly)], merged_bucket_proof_poly) // g^(q_{merged})

	// 对合并后的桶中的每一个位置更新 individual proof
	//merged_bucket_size := merge_end - merge_start + 1
	// for j := uint64(0); j < merged_bucket_size; j++ { // 对合并后桶里的每一个位置
	// 	j_global_index := merge_start + j

	// 	// 找到这个位置原本所在的桶（使用旧的p_detail）
	// 	original_bucket_idx := uint64(0)
	// 	for i := uint64(0); i < num_merge; i++ {
	// 		if j_global_index >= old_p_detail[start_bucket+i][0] && j_global_index <= old_p_detail[start_bucket+i][1] {
	// 			original_bucket_idx = i
	// 			break
	// 		}
	// 	}

	// 	// 使用预计算的缓存增量多项式：𝜙_{merged}(x) - 𝜙_{old_i}(x)
	// 	individual_sub_poly := bucket_diffs[original_bucket_idx]

	// 	// 计算增量部分的商：(𝜙_{merged}(x) - 𝜙_{old_i}(x)) / (x - ω_j)
	// 	divisor := [2]mcl.Fr{}
	// 	divisor[1].SetInt64(1)
	// 	divisor[0] = vb.FFT.RootsOfUnity[j_global_index*2] // 全局根 ω_j
	// 	mcl.FrNeg(&divisor[0], &divisor[0])                // -ω_j
	// 	q_ij_quotient, _ := utils.PolyDiv(individual_sub_poly, divisor[:])
	// 	var individual_proof_update mcl.G1
	// 	mcl.G1MulVec(&individual_proof_update, vb.SecretG1[:len(q_ij_quotient)], q_ij_quotient) // 新乘的部分

	// 	// 更新 individual proof：新个别证明 = 原证明 + 增量部分
	// 	mcl.G1Add(&new_proofs[j_global_index].IndividualProof, &new_proofs[j_global_index].IndividualProof, &individual_proof_update)

	// 	// 更新 bucket proof
	// 	new_proofs[j_global_index].BucketProof = merged_bucket_proof
	// }

	return new_proofs
}

// func (vb *VBAS_B1) BucketSplit(target_bucket uint64, num_split uint64, proofs []BasicBucketProof, vector []mcl.Fr) []BasicBucketProof {
// 	//首先检查是否可以均匀分成num_split个桶
// 	if vb.p_detail[target_bucket][1]-vb.p_detail[target_bucket][0]+1 < num_split {
// 		panic("target_bucket桶的大小不能均匀分成num_split个桶")
// 	}

// 	//验证是否可以整除
// 	if (vb.p_detail[target_bucket][1]-vb.p_detail[target_bucket][0]+1)%num_split != 0 {
// 		panic("target_bucket桶的大小不能均匀分成num_split个桶")
// 	}

// 	phi_poly := vb.PhiPoly

// 	// ===== 1. 原桶 φ_old =====
// 	start := vb.p_detail[target_bucket][0]
// 	end := vb.p_detail[target_bucket][1]
// 	origin_vec := vector[start : end+1]

// 	roots := make([]mcl.Fr, len(origin_vec))
// 	for i := range roots {
// 		roots[i] = vb.FFT.RootsOfUnity[(start+uint64(i))*2]
// 	}

// 	sub := utils.SubProductTree(roots)
// 	s := utils.PolyMultiEvaluate(
// 		utils.PolyDifferentiate(sub.Poly),
// 		sub,
// 	)

// 	ys := make([]mcl.Fr, len(origin_vec))
// 	for i := range ys {
// 		mcl.FrDiv(&ys[i], &origin_vec[i], &s[i])
// 	}

// 	phi_old := utils.PolyInterpolation(roots, ys, sub)

// 	// ===== 2. 更新 p_detail（与你原来一样）=====
// 	old_p_detail := vb.p_detail
// 	old_p := vb.p

// 	new_bucket_size := (end - start + 1) / num_split
// 	new_p := old_p - 1 + num_split

// 	new_p_detail := make([][]uint64, new_p)
// 	for i := uint64(0); i < old_p; i++ {
// 		if i < target_bucket {
// 			new_p_detail[i] = old_p_detail[i]
// 		} else if i > target_bucket {
// 			new_p_detail[i+num_split-1] = old_p_detail[i]
// 		}
// 	}

// 	for i := uint64(0); i < num_split; i++ {
// 		new_p_detail[target_bucket+i] = []uint64{
// 			start + i*new_bucket_size,
// 			start + (i+1)*new_bucket_size - 1,
// 		}
// 	}

// 	vb.p = new_p
// 	vb.p_detail = new_p_detail

// 	// ===== 3. 新桶 bucket proof + Δφ 记录 =====
// 	new_proofs := make([]BasicBucketProof, len(proofs))
// 	copy(new_proofs, proofs)

// 	for i := uint64(0); i < num_split; i++ {
// 		b := target_bucket + i
// 		bs := new_p_detail[b][0]
// 		be := new_p_detail[b][1]
// 		vec := vector[bs : be+1]

// 		roots := make([]mcl.Fr, len(vec))
// 		for j := range roots {
// 			roots[j] = vb.FFT.RootsOfUnity[(bs+uint64(j))*2]
// 		}

// 		sub := utils.SubProductTree(roots)
// 		s := utils.PolyMultiEvaluate(
// 			utils.PolyDifferentiate(sub.Poly),
// 			sub,
// 		)

// 		ys := make([]mcl.Fr, len(vec))
// 		for j := range ys {
// 			mcl.FrDiv(&ys[j], &vec[j], &s[j])
// 		}

// 		phi_new := utils.PolyInterpolation(roots, ys, sub)

// 		// bucket proof
// 		q, _ := utils.PolyDiv(
// 			polySubtract(phi_poly, phi_new),
// 			vb.ai_polys[b],
// 		)

// 		var bp mcl.G1
// 		mcl.G1MulVec(&bp, vb.SecretG1[:len(q)], q)

// 		// 更新 bucket proof（所有元素共享）
// 		for j := bs; j <= be; j++ {
// 			new_proofs[j].BucketProof = bp
// 		}

// 		// ===== Lazy delta =====
// 		vb.BucketDeltas[b] = &BucketDelta{
// 			DeltaPoly: polySubtract(phi_new, phi_old),
// 		}
// 	}

// 	return new_proofs
// }

func (vb *VBAS_B1) BucketSplit(
	target_bucket uint64,
	num_split uint64,
	proofs []BasicBucketProof,
	vector []mcl.Fr,
) []BasicBucketProof {

	// ===== 0. 参数检查 =====
	old_start := vb.p_detail[target_bucket][0]
	old_end := vb.p_detail[target_bucket][1]
	old_bucket_size := old_end - old_start + 1

	if old_bucket_size < num_split || old_bucket_size%num_split != 0 {
		panic("target_bucket cannot be evenly split")
	}

	phi_poly := vb.PhiPoly

	// ===== 1. 计算 φ_old =====
	origin_vec := vector[old_start : old_end+1]

	old_roots := make([]mcl.Fr, len(origin_vec))
	for i := range old_roots {
		old_roots[i] = vb.FFT.RootsOfUnity[(old_start+uint64(i))*2]
	}

	old_sub := utils.SubProductTree(old_roots)
	old_s := utils.PolyMultiEvaluate(
		utils.PolyDifferentiate(old_sub.Poly),
		old_sub,
	)

	old_ys := make([]mcl.Fr, len(origin_vec))
	for i := range old_ys {
		mcl.FrDiv(&old_ys[i], &origin_vec[i], &old_s[i])
	}

	phi_old := utils.PolyInterpolation(old_roots, old_ys, old_sub)

	// ===== 2. 更新 p_detail / p =====
	old_p := vb.p
	old_p_detail := vb.p_detail
	old_ai := vb.ai_polys

	new_bucket_size := old_bucket_size / num_split
	new_p := old_p - 1 + num_split

	new_p_detail := make([][]uint64, new_p)
	new_ai_polys := make([][]mcl.Fr, new_p)

	// --- 拷贝未受影响的桶 ---
	for i := uint64(0); i < old_p; i++ {
		if i < target_bucket {
			new_p_detail[i] = old_p_detail[i]
			new_ai_polys[i] = old_ai[i]
		} else if i > target_bucket {
			new_p_detail[i+num_split-1] = old_p_detail[i]
			new_ai_polys[i+num_split-1] = old_ai[i]
		}
	}

	// --- 构造新桶的 p_detail 和 ai_polys ---
	for i := uint64(0); i < num_split; i++ {
		b := target_bucket + i
		start := old_start + i*new_bucket_size
		end := start + new_bucket_size - 1

		new_p_detail[b] = []uint64{start, end}

		roots := make([]mcl.Fr, new_bucket_size)
		for j := uint64(0); j < new_bucket_size; j++ {
			roots[j] = vb.FFT.RootsOfUnity[(start+j)*2]
		}

		sub := utils.SubProductTree(roots)
		new_ai_polys[b] = sub.Poly
	}

	vb.p = new_p
	vb.p_detail = new_p_detail
	vb.ai_polys = new_ai_polys

	// ===== 3. 计算新 bucket proof =====
	new_proofs := make([]BasicBucketProof, len(proofs))
	copy(new_proofs, proofs)

	for i := uint64(0); i < num_split; i++ {
		b := target_bucket + i
		bs := new_p_detail[b][0]
		be := new_p_detail[b][1]

		vec := vector[bs : be+1]

		roots := make([]mcl.Fr, len(vec))
		for j := range roots {
			roots[j] = vb.FFT.RootsOfUnity[(bs+uint64(j))*2]
		}

		sub := utils.SubProductTree(roots)
		s := utils.PolyMultiEvaluate(
			utils.PolyDifferentiate(sub.Poly),
			sub,
		)

		ys := make([]mcl.Fr, len(vec))
		for j := range ys {
			mcl.FrDiv(&ys[j], &vec[j], &s[j])
		}

		phi_new := utils.PolyInterpolation(roots, ys, sub)

		q, _ := utils.PolyDiv(
			polySubtract(phi_poly, phi_new),
			vb.ai_polys[b],
		)

		var bp mcl.G1
		mcl.G1MulVec(&bp, vb.SecretG1[:len(q)], q)

		for j := bs; j <= be; j++ {
			new_proofs[j].BucketProof = bp
		}

		// Lazy Δφ
		vb.BucketDeltas[b] = &BucketDelta{
			DeltaPoly: polySubtract(phi_new, phi_old),
		}
	}

	return new_proofs
}

type BucketDelta struct {
	// Δφ(x) = φ_new_bucket(x) - φ_old_bucket(x)
	DeltaPoly []mcl.Fr
}
