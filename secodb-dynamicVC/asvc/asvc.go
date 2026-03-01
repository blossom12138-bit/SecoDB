package asvc

import (
    "github.com/alinush/go-mcl"
    "github.com/wangnick2017/balanceproofs-go/fft"
    "github.com/wangnick2017/balanceproofs-go/utils"
)

type UpdateReq struct {
    Index uint64
    Delta mcl.Fr
}

type Inp struct {
    Index uint64
    Proof mcl.G1
}

type Val struct {
    Index uint64
    Y     mcl.Fr
}

type ASVC struct {
    N        uint64        //n=2^L 向量大小
    L        uint8         //log2(n)，向量的对数大小
    SecretG1 []mcl.G1      //秘密密钥在G1上
    SecretG2 []mcl.G2      //秘密密钥在G2上
    G        mcl.G1        //G1生成元
    H        mcl.G2        //G2生成元
    FFT      *fft.Settings //FFT设置
    L_i      []mcl.G1      //拉格朗日承诺值
    U_i      []mcl.G1      //更新密钥用的数值
    A_i      []mcl.G1      //openones的结果
    a        mcl.G1        //A(tao)的值
    xExtFFT  []mcl.G1      //FFT加速用的扩展值
}

func (asvc *ASVC) Init(L uint8, s1 []mcl.G1, s2 []mcl.G2) {
    asvc.L = L
    asvc.N = 1 << L
    if len(s1) == 0 {
        asvc.SecretG1, asvc.SecretG2 = utils.GenerateTestingSetup(asvc.N) // 输出SecretG1=[g^(τ^i)]和SecretG2=[h^(τ^i)]
    } else {
        asvc.SecretG1 = s1
        asvc.SecretG2 = s2
    }

    asvc.G = asvc.SecretG1[0]
    asvc.H = asvc.SecretG2[0]
    asvc.FFT = fft.NewSettings(L + 1) //FFT计算器创建

    //ToeplitzPart1
    x := make([]mcl.G1, asvc.N, asvc.N) //x=[g^(τ^(N-1)),g^(τ^(N-2)),...,g^(τ^0)]，构造反向的秘密参数列表为x
    for i, j := uint64(0), asvc.N-2; i < asvc.N-1; i, j = i+1, j-1 {
        x[i] = asvc.SecretG1[j]
    }
    x[asvc.N-1].Clear()
    n2 := asvc.N * 2
    xExt := make([]mcl.G1, n2, n2) //xExt=[g^(τ^(N-1)),g^(τ^(N-2)),...,g^(τ^0),0,...,0]，构造双倍的反向的秘密参数列表为xExt，其中后半部分为0
    for i := uint64(0); i < asvc.N; i++ {
        xExt[i] = x[i]
    }
    for i := asvc.N; i < n2; i++ {
        xExt[i].Clear()
    }
    asvc.xExtFFT = asvc.FFT.FFT_G1(xExt, false)

    //compute asvc.L_i, a_i, u_i, a
    mcl.G1Sub(&asvc.a, &asvc.SecretG1[asvc.N], &asvc.SecretG1[0]) //a=g^(τ^(N-1))/g^(τ^0)=g^(τ^(N-1))=g^A(τ)
    asvc.L_i = asvc.FFT.FFT_G1(asvc.SecretG1[:asvc.N], true)      //公共参数g^L_i(τ)
    asvc.A_i = asvc.OpenOnes()                                    //计算公共参数a_i=g^(A(τ)/x-w^i)

    h := make([]mcl.G1, asvc.N)
    for i := uint64(0); i <= asvc.N-2; i++ {
        var j mcl.Fr
        j.SetInt64(int64(i + 1))
        mcl.G1Mul(&h[i], &asvc.SecretG1[asvc.N-2-i], &j)
    }
    h[asvc.N-1].Clear()

    asvc.U_i = asvc.FFT.FFT_G1(h, false) //计算公共参数u_i=g^(L_i(τ)-1/τ-w^i)
    for i := uint64(0); i < asvc.N; i++ {
        var nFr mcl.Fr
        nFr.SetInt64(int64(asvc.N))
        mcl.FrDiv(&nFr, &asvc.FFT.RootsOfUnity[i*2], &nFr)
        mcl.G1Mul(&asvc.U_i[i], &asvc.U_i[i], &nFr)
    }
}

// 计算对整个向量vector的承诺值
func (asvc *ASVC) Commit(vector []mcl.Fr) mcl.G1 {
    var digest mcl.G1
    mcl.G1MulVec(&digest, asvc.L_i, vector)

    return digest
}

// Open could also be implemented in O(n) time through n calls of UpdateProof()
// 生成vector[index]的证明
func (asvc *ASVC) Open(index uint64, vector []mcl.Fr) mcl.G1 {
    poly := asvc.FFT.FFT_Fr(vector, true) //执行FFT逆运算得到整个向量的多项式

    // divisor = [-index, 1]
    divisor := [2]mcl.Fr{}                      //长度为2的Fr数组，{}表示初始化一个长度为2的Fr数组，默认值为0
    divisor[1].SetInt64(1)                      // divisor[1] = 1
    divisor[0] = asvc.FFT.RootsOfUnity[index*2] //divisor[0] = ω_i
    mcl.FrNeg(&divisor[0], &divisor[0])         //divisor[0] = -ω_i,FrNeg是将第二个输入取反，放入第一个输入中
    // quot = poly / divisor
    quotientPolynomial, _ := utils.PolyDiv(poly, divisor[:]) //计算多项式除法，整条向量的多项式除以，除数为(X - ω_i)
    //注意这里的divisor被多项式 (X - ω_i) 的系数表示为 [-ω_i, 1]即：(X - ω_i) = -ω_i + 1·X = divisor[0] + divisor[1]·X

    // evaluate quotient poly at shared secret, in G1
    var proof mcl.G1
    mcl.G1MulVec(&proof, asvc.SecretG1[:len(quotientPolynomial)], quotientPolynomial)
    return proof
}

// 快速计算a_i=g^(A(τ)/x-w^i)
func (asvc *ASVC) OpenOnes() []mcl.G1 {
    //ToeplitzPart1
    x := make([]mcl.G1, asvc.N, asvc.N)
    for i, j := uint64(0), asvc.N-1; i < asvc.N; i, j = i+1, j-1 {
        x[i] = asvc.SecretG1[j]
    }
    n2 := asvc.N * 2
    xExt := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < asvc.N; i++ {
        xExt[i] = x[i]
    }
    for i := asvc.N; i < n2; i++ {
        xExt[i].Clear()
    }
    xExtFFT := asvc.FFT.FFT_G1(xExt, false)

    //ToeplitzcoeffsStep
    // [last poly item] + [0]*(n+1) + [poly items except first and last]
    toeplitzCoeffs := make([]mcl.Fr, n2, n2)
    toeplitzCoeffs[0].SetInt64(1)
    for i := uint64(1); i < n2; i++ {
        toeplitzCoeffs[i].Clear()
    }

    //ToeplitzPart2
    toeplitzCoeffsFFT := asvc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := asvc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return asvc.FFT.FFT_G1(h, false)
}

// 快速计算对整个向量vector的证明
func (asvc *ASVC) OpenAll(vector []mcl.Fr) []mcl.G1 {
    poly := asvc.FFT.FFT_Fr(vector, true)

    //ToeplitzcoeffsStep
    n := asvc.N
    n2 := n * 2
    // [last poly item] + [0]*(n+1) + [poly items except first and last]
    toeplitzCoeffs := make([]mcl.Fr, n2, n2)
    toeplitzCoeffs[0] = poly[n-1]
    for i := uint64(1); i <= n+1; i++ {
        toeplitzCoeffs[i].Clear()
    }
    for i, j := n+2, 1; i < n2; i, j = i+1, j+1 {
        toeplitzCoeffs[i] = poly[j]
    }

    //ToeplitzPart2
    toeplitzCoeffsFFT := asvc.FFT.FFT_Fr(toeplitzCoeffs, false)
    hExtFFT := make([]mcl.G1, n2, n2)
    for i := uint64(0); i < n2; i++ {
        mcl.G1Mul(&hExtFFT[i], &asvc.xExtFFT[i], &toeplitzCoeffsFFT[i])
    }

    //h := fk.ToeplitzPart3(hExtFFT)
    out := asvc.FFT.FFT_G1(hExtFFT, true)
    // Only the top half is the Toeplitz product, the rest is padding
    h := out[:len(out)/2]

    return asvc.FFT.FFT_G1(h, false)
}

type OpenAllStepArg struct {
    Vector []mcl.Fr
    Step   uint64
    Proofs []mcl.G1
    Done   bool
    Count  uint64

    toeplitzCoeffs    []mcl.Fr
    toeplitzCoeffsFFT []mcl.Fr
    hExtFFT           []mcl.G1
    h                 []mcl.G1
    firstFFTDone      bool
    fftArg            *fft.FFT_G1_Arg
}

// 分多步完成OpenAll
func (asvc *ASVC) OpenAllStep(openArg *OpenAllStepArg) {
    n := asvc.N
    n2 := n * 2

    switch {
    case openArg.Step == 0:
        openArg.firstFFTDone = false
        poly := asvc.FFT.FFT_Fr(openArg.Vector, true)

        //ToeplitzcoeffsStep

        // [last poly item] + [0]*(n+1) + [poly items except first and last]
        openArg.toeplitzCoeffs = make([]mcl.Fr, n2, n2)
        openArg.toeplitzCoeffs[0] = poly[n-1]
        for i := uint64(1); i <= n+1; i++ {
            openArg.toeplitzCoeffs[i].Clear()
        }
        for i, j := n+2, 1; i < n2; i, j = i+1, j+1 {
            openArg.toeplitzCoeffs[i] = poly[j]
        }

    case openArg.Step == 1:
        //ToeplitzPart2
        openArg.toeplitzCoeffsFFT = asvc.FFT.FFT_Fr(openArg.toeplitzCoeffs, false)

    case openArg.Step > 1 && openArg.Step <= 1+n2/openArg.Count:
        if openArg.Step == 2 {
            openArg.hExtFFT = make([]mcl.G1, n2, n2)
        }
        for i, j := (openArg.Step-2)*openArg.Count, uint64(0); j < openArg.Count; i, j = i+1, j+1 {
            mcl.G1Mul(&openArg.hExtFFT[i], &asvc.xExtFFT[i], &openArg.toeplitzCoeffsFFT[i])
        }

    default:
        if openArg.Step == 2+n2/openArg.Count {
            openArg.fftArg = &fft.FFT_G1_Arg{
                N:       n2,
                L:       uint64(asvc.L) + 1,
                Vals:    openArg.hExtFFT,
                Out:     nil,
                Inv:     true,
                Pre_s:   0,
                Pre_j:   0,
                Pre_inv: 0,
                Count:   openArg.Count,
                Done:    false,
            }
        }
        if !openArg.firstFFTDone {
            asvc.FFT.FFT_G1_Seg(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.h = openArg.fftArg.Out[:n]
                openArg.firstFFTDone = true
                openArg.fftArg = nil
            }
        } else {
            if openArg.fftArg == nil {
                openArg.fftArg = &fft.FFT_G1_Arg{
                    N:       n,
                    L:       uint64(asvc.L),
                    Vals:    openArg.h,
                    Out:     nil,
                    Inv:     false,
                    Pre_s:   0,
                    Pre_j:   0,
                    Pre_inv: 0,
                    Count:   openArg.Count,
                    Done:    false,
                }
            }
            asvc.FFT.FFT_G1_Seg(openArg.fftArg)
            if openArg.fftArg.Done {
                openArg.Proofs = openArg.fftArg.Out
                openArg.Done = true
            }
        }
    }
}

// 更新整个向量的承诺
func (asvc *ASVC) UpdateCommitment(digest mcl.G1, req UpdateReq) mcl.G1 {
    var temp mcl.G1
    mcl.G1Mul(&temp, &asvc.L_i[req.Index], &req.Delta)
    mcl.G1Add(&temp, &digest, &temp)
    return temp
}

// 更新单个元素的证明
func (asvc *ASVC) UpdateProof(proof mcl.G1, index uint64, req UpdateReq) mcl.G1 {
    if index == req.Index {
        var temp mcl.G1
        mcl.G1Mul(&temp, &asvc.U_i[index], &req.Delta)
        mcl.G1Add(&temp, &proof, &temp)
        return temp
    }

    var omega_i mcl.Fr
    mcl.FrSub(&omega_i, &asvc.FFT.RootsOfUnity[req.Index*2], &asvc.FFT.RootsOfUnity[index*2]) //ω_j-ω_i

    var w_i_j mcl.G1
    mcl.G1Sub(&w_i_j, &asvc.A_i[req.Index], &asvc.A_i[index]) //a_j/a_i 注意这里G1上的sub实际上在论文里相当于除
    //两个都反过来了，等于文章里没反过来的写法
    var n, ao mcl.Fr
    n.SetInt64(int64(asvc.N))
    mcl.FrDiv(&ao, &asvc.FFT.RootsOfUnity[req.Index*2], &n) //w_j/n
    mcl.FrMul(&ao, &ao, &req.Delta)                         //ao = (ω_j / n) * δ
    mcl.FrDiv(&omega_i, &ao, &omega_i)                      //omega_i = ao / (ω_j - ω_i)

    mcl.G1Mul(&w_i_j, &w_i_j, &omega_i) //w_i_j = w_i_j ^ omega_i, 这里的w_i_j是a_j/a_i，omega_i已经是次方上的集成版本

    var temp mcl.G1
    mcl.G1Add(&temp, &proof, &w_i_j) //G1Add对应论文里的点乘操作
    return temp
}

// 聚合多个证明
func (asvc *ASVC) Aggregate(aggs []Inp) mcl.G1 {
    l := len(aggs) //聚合的证明数量
    if l == 1 {
        return aggs[0].Proof
    }
    coef := make([]mcl.Fr, l) //系数数组，长度为l，每个元素为ω_i
    for i := 0; i < l; i++ {
        coef[i] = asvc.FFT.RootsOfUnity[aggs[i].Index*2]
    }
    a_I := utils.SubProductTree(coef)              //A_I(X)=∏(X-ω_i)
    a_I_prime := utils.PolyDifferentiate(a_I.Poly) //A_I'(X)
    c := utils.PolyMultiEvaluate(a_I_prime, a_I)   //A_I'(X)在ω_i处的值，即 A_I'(ω_i)

    var res mcl.G1
    res.Clear()
    for i := 0; i < l; i++ {
        var c_pi mcl.G1
        mcl.FrInv(&c[i], &c[i]) //根据c_i=1/A_I'(ω_i)，计算上面c里每个位置的逆元A_I'(X)^-1
        p := aggs[i].Proof
        mcl.G1Mul(&c_pi, &p, &c[i]) //c_pi = p_i ^ c_i
        mcl.G1Add(&res, &res, &c_pi)
    }
    return res
}

// 验证单个位置的证明
func (asvc *ASVC) VerifySingle(digest mcl.G1, proof mcl.G1, v Val) bool {
    var i mcl.Fr
    i = asvc.FFT.RootsOfUnity[v.Index*2] //ω_i
    var xG2 mcl.G2
    mcl.G2Mul(&xG2, &asvc.H, &i) //xG2=g^w_i
    var sMinuxX mcl.G2
    mcl.G2Sub(&sMinuxX, &asvc.SecretG2[1], &xG2) //g^tao/g^w_i = A_I(tao)，当I={v.Index}时

    var yG1 mcl.G1
    mcl.G1Mul(&yG1, &asvc.G, &v.Y) //yG1=g^v_i，这里v_i是要验证位置的值，其实是g^R_I(X),但单个变量时为R_i(X)=v_i
    var commitmentMinusY mcl.G1
    mcl.G1Sub(&commitmentMinusY, &digest, &yG1) //commitment / yG1 = c / g^v_i

    // e([commitment - y], [1]) = e([proof],  [s - x])
    var e1, e2 mcl.GT
    mcl.Pairing(&e1, &commitmentMinusY, &asvc.H)
    mcl.Pairing(&e2, &proof, &sMinuxX)
    return e1.IsEqual(&e2)
}

// 验证聚合证明
func (asvc *ASVC) VerifyAggregation(digest mcl.G1, proof mcl.G1, aggvs []Val) bool {
    l := len(aggvs)
    x := make([]mcl.Fr, l)

    for i := 0; i < l; i++ {
        x[i] = asvc.FFT.RootsOfUnity[aggvs[i].Index*2] //ω_i，i \in I
    }

    a_I := utils.SubProductTree(x) //A_I(X)=∏(X-ω_i)，i \in I
    a_I_prime := utils.PolyDifferentiate(a_I.Poly)
    s := utils.PolyMultiEvaluate(a_I_prime, a_I) //A_I'(X)在ω_i处的值，即 A_I'(ω_i)

    ys := make([]mcl.Fr, l)
    for i := 0; i < l; i++ {
        mcl.FrDiv(&ys[i], &aggvs[i].Y, &s[i]) // ys[i] = v_i / A_I'(ω_i)
    }
    interpolationPoly := utils.PolyInterpolation(x, ys, a_I) //拉格朗日多项式插值，得到R_I(X)，R_I(w_i)=v_i

    subProductPoly := a_I.Poly
    var sp2 mcl.G2
    mcl.G2MulVec(&sp2, asvc.SecretG2[:len(subProductPoly)], subProductPoly) //sp2=g^(A_I(tao))

    // [interpolation_polynomial(s)]_1
    var is1 mcl.G1
    mcl.G1MulVec(&is1, asvc.SecretG1[:len(interpolationPoly)], interpolationPoly) //is1=g^(R_I(tao))
    // [commitment - interpolation_polynomial(s)]_1 = [commit]_1 - [interpolation_polynomial(s)]_1
    var commitMinusInterpolation mcl.G1
    mcl.G1Sub(&commitMinusInterpolation, &digest, &is1) //commitment / is1 = c / g^(R_I(tao))

    // Verify the pairing equation
    //
    // e([commitment - interpolation_polynomial(s)], [1]) = e([proof],  [s^n - x^n]),原作者写的，不知道这个s^n - x^n是什么意思
    // 等价于 e([c / g^(R_I(tao))], [g]) = e([proof],  [g^(A_I(tao)])
    var e1, e2 mcl.GT
    mcl.Pairing(&e1, &commitMinusInterpolation, &asvc.H)
    mcl.Pairing(&e2, &proof, &sp2)
    return e1.IsEqual(&e2)
}
