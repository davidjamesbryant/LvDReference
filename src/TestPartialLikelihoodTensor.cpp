// Unit tests for PartialLikelihoodTensor.
//
// Strategy: each of the 12 dispatch paths in product() is tested by
// comparing against a simple reference loop implementation (matching
// the reference version commented out at the bottom of the header).
// Multisite tests use S=3 sites with per-site varying data so that
// wrong site-offset calculations produce detectable failures.

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include "PartialLikelihoodTensor.h"

static int nPass = 0, nFail = 0;

static bool approxEq(const std::vector<Scalar>& got, const std::vector<Scalar>& exp, double tol = 1e-10)
{
    if (got.size() != exp.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i)
        if (std::fabs(got[i] - exp[i]) > tol * (1.0 + std::fabs(exp[i]))) return false;
    return true;
}

static void check(const char* name, const std::vector<Scalar>& got, const std::vector<Scalar>& expected)
{
    bool ok = approxEq(got, expected);
    std::cout << (ok ? "PASS" : "FAIL") << ": " << name << "\n";
    if (!ok) {
        std::cout << "  sizes: got=" << got.size() << " expected=" << expected.size() << "\n";
        for (std::size_t i = 0; i < std::min(got.size(), expected.size()); ++i)
            if (std::fabs(got[i] - expected[i]) > 1e-10 * (1.0 + std::fabs(expected[i])))
                std::cout << "  [" << i << "] got=" << got[i] << " expected=" << expected[i] << "\n";
    }
    ok ? ++nPass : ++nFail;
}

// ============================================================
// Reference implementations — simple loops matching the reference
// version at the bottom of PartialLikelihoodTensor.h
// ============================================================

static std::vector<Scalar> ref_vec_vec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    std::vector<Scalar> C(4);
    for (int i = 0; i < 4; i++) C[i] = A[i] * B[i];
    return C;
}

static std::vector<Scalar> ref_vec_multivec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)B.size() / 4;
    std::vector<Scalar> C(4*s);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            C[4*p+i] = A[i] * B[4*p+i];
    return C;
}

static std::vector<Scalar> ref_multivec_multivec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)A.size() / 4;
    std::vector<Scalar> C(4*s);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            C[4*p+i] = A[4*p+i] * B[4*p+i];
    return C;
}

static std::vector<Scalar> ref_mat_vec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    std::vector<Scalar> C(4, 0.0);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            C[i] += A[4*i+j] * B[j];
    return C;
}

static std::vector<Scalar> ref_mat_multivec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)B.size() / 4;
    std::vector<Scalar> C(4*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                C[4*p+i] += A[4*i+j] * B[4*p+j];
    return C;
}

static std::vector<Scalar> ref_mat_mat(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    std::vector<Scalar> C(16, 0.0);
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                C[4*i+j] += A[4*i+k] * B[4*k+j];
    return C;
}

static std::vector<Scalar> ref_mat_multimat(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)B.size() / 16;
    std::vector<Scalar> C(16*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    C[16*p+4*i+j] += A[4*i+k] * B[16*p+4*k+j];
    return C;
}

static std::vector<Scalar> ref_multimat_vec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)A.size() / 16;
    std::vector<Scalar> C(4*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                C[4*p+i] += A[16*p+4*i+j] * B[j];
    return C;
}

static std::vector<Scalar> ref_multimat_multivec(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)A.size() / 16;
    std::vector<Scalar> C(4*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int k = 0; k < 4; k++)
                C[4*p+i] += A[16*p+4*i+k] * B[4*p+k];
    return C;
}

static std::vector<Scalar> ref_multimat_mat(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)A.size() / 16;
    std::vector<Scalar> C(16*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    C[16*p+4*i+j] += A[16*p+4*i+k] * B[4*k+j];
    return C;
}

static std::vector<Scalar> ref_multimat_multimat(const std::vector<Scalar>& A, const std::vector<Scalar>& B)
{
    int s = (int)A.size() / 16;
    std::vector<Scalar> C(16*s, 0.0);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                for (int k = 0; k < 4; k++)
                    C[16*p+4*i+j] += A[16*p+4*i+k] * B[16*p+4*k+j];
    return C;
}

// ============================================================
// Data generators
//
// Storage layout matches PartialLikelihoodTensor:
//   vector site p:  data[4*p + i]       (i in 0..3)
//   matrix site p:  data[16*p + 4*i+j]  (row-major, i=row, j=col)
//
// Per-site values differ so wrong site-offset bugs are detectable.
// ============================================================

// Single vec: {1, 2, 3, 4} * scale
static std::vector<Scalar> singleVec(Scalar scale = 1.0)
{
    return {scale*1, scale*2, scale*3, scale*4};
}

// Multisite vec: site p has components (4p + i + 1) * scale
static std::vector<Scalar> multiVec(int s, Scalar scale = 0.1)
{
    std::vector<Scalar> v(4*s);
    for (int p = 0; p < s; p++)
        for (int i = 0; i < 4; i++)
            v[4*p+i] = scale * (4*p + i + 1);
    return v;
}

// Single mat: element k (row-major) = (k + 1) * scale
static std::vector<Scalar> singleMat(Scalar scale = 0.1)
{
    std::vector<Scalar> m(16);
    for (int k = 0; k < 16; k++) m[k] = scale * (k + 1);
    return m;
}

// Multisite mat: site p element k = (16p + k + 1) * scale
static std::vector<Scalar> multiMat(int s, Scalar scale = 0.01)
{
    std::vector<Scalar> m(16*s);
    for (int p = 0; p < s; p++)
        for (int k = 0; k < 16; k++)
            m[16*p+k] = scale * (16*p + k + 1);
    return m;
}

// ============================================================
// Tests — one per dispatch path in product()
// ============================================================

int main()
{
    const int S = 3;

    // Case 1: vec_vec  (clade !ms × clade !ms → clade !ms)
    {
        auto aD = singleVec(0.5), bD = singleVec(0.3);
        PartialLikelihoodTensor a(true, false, 1), b(true, false, 1), c(true, false, 1);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("vec_vec", c.getData(), ref_vec_vec(aD, bD));
    }

    // Case 2: vec_multivec  (clade !ms × clade ms → clade ms)
    {
        auto aD = singleVec(0.5);
        auto bD = multiVec(S, 0.1);
        PartialLikelihoodTensor a(true, false, 1), b(true, true, S), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("vec_multivec", c.getData(), ref_vec_multivec(aD, bD));
    }

    // Case 3: multivec_vec  (clade ms × clade !ms → clade ms)
    //   dispatch calls vec_multivec(b.data, a.data, c.data) since multiply commutes
    {
        auto aD = multiVec(S, 0.1);
        auto bD = singleVec(0.7);
        PartialLikelihoodTensor a(true, true, S), b(true, false, 1), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multivec_vec", c.getData(), ref_vec_multivec(bD, aD));
    }

    // Case 4: multivec_multivec  (clade ms × clade ms → clade ms)
    {
        auto aD = multiVec(S, 0.1);
        auto bD = multiVec(S, 0.2);
        PartialLikelihoodTensor a(true, true, S), b(true, true, S), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multivec_multivec", c.getData(), ref_multivec_multivec(aD, bD));
    }

    // Case 5: mat_vec  (!clade !ms × clade !ms → clade !ms)
    {
        auto aD = singleMat(0.1);
        auto bD = singleVec(1.0);
        PartialLikelihoodTensor a(false, false, 1), b(true, false, 1), c(true, false, 1);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("mat_vec", c.getData(), ref_mat_vec(aD, bD));
    }

    // Case 6: mat_multivec  (!clade !ms × clade ms → clade ms)
    {
        auto aD = singleMat(0.1);
        auto bD = multiVec(S, 0.1);
        PartialLikelihoodTensor a(false, false, 1), b(true, true, S), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("mat_multivec", c.getData(), ref_mat_multivec(aD, bD));
    }

    // Case 7: mat_mat  (!clade !ms × !clade !ms → !clade !ms)
    {
        auto aD = singleMat(0.1);
        auto bD = singleMat(0.2);
        PartialLikelihoodTensor a(false, false, 1), b(false, false, 1), c(false, false, 1);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("mat_mat", c.getData(), ref_mat_mat(aD, bD));
    }

    // Case 8: mat_multimat  (!clade !ms × !clade ms → !clade ms)
    {
        auto aD = singleMat(0.1);
        auto bD = multiMat(S, 0.01);
        PartialLikelihoodTensor a(false, false, 1), b(false, true, S), c(false, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("mat_multimat", c.getData(), ref_mat_multimat(aD, bD));
    }

    // Case 9: multimat_vec  (!clade ms × clade !ms → clade ms)
    {
        auto aD = multiMat(S, 0.01);
        auto bD = singleVec(1.0);
        PartialLikelihoodTensor a(false, true, S), b(true, false, 1), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multimat_vec", c.getData(), ref_multimat_vec(aD, bD));
    }

    // Case 10: multimat_multivec  (!clade ms × clade ms → clade ms)
    {
        auto aD = multiMat(S, 0.01);
        auto bD = multiVec(S, 0.1);
        PartialLikelihoodTensor a(false, true, S), b(true, true, S), c(true, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multimat_multivec", c.getData(), ref_multimat_multivec(aD, bD));
    }

    // Case 11: multimat_mat  (!clade ms × !clade !ms → !clade ms)
    {
        auto aD = multiMat(S, 0.01);
        auto bD = singleMat(0.2);
        PartialLikelihoodTensor a(false, true, S), b(false, false, 1), c(false, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multimat_mat", c.getData(), ref_multimat_mat(aD, bD));
    }

    // Case 12: multimat_multimat  (!clade ms × !clade ms → !clade ms)
    {
        auto aD = multiMat(S, 0.01);
        auto bD = multiMat(S, 0.02);
        PartialLikelihoodTensor a(false, true, S), b(false, true, S), c(false, true, S);
        a.setData(aD); b.setData(bD);
        product(a, b, c);
        check("multimat_multimat", c.getData(), ref_multimat_multimat(aD, bD));
    }

    std::cout << "\n" << nPass << "/" << (nPass + nFail) << " passed\n";
    return nFail == 0 ? 0 : 1;
}
