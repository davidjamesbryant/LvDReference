//
// Created by David Bryant on 14/05/2026.
// Optimized site-contiguous kernels by ChatGPT.
//

#ifndef LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H
#define LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H

#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <Eigen/Dense>

#include "phylib.h"

#if defined(__GNUC__) || defined(__clang__)
    #define PLT_RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define PLT_RESTRICT __restrict
#else
    #define PLT_RESTRICT
#endif

class PartialLikelihoodTensor {
private:
    std::vector<Scalar> data;
    bool isClade;
    bool multisite;
    int nSites;

    static inline std::size_t siteCountFromN(int n) noexcept
    {
        assert(n > 0);
        return static_cast<std::size_t>(n);
    }

    static inline void assertNoAlias(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& C)
    {
        assert(A.data() != C.data());
    }

    static inline void assertNoAlias(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        const std::vector<Scalar>& C)
    {
        assert(A.data() != C.data());
        assert(B.data() != C.data());
    }

    static inline void assertProductShape(
        const PartialLikelihoodTensor& a,
        const PartialLikelihoodTensor& b,
        const PartialLikelihoodTensor& c)
    {
        assert(a.nSites > 0);
        assert(b.nSites > 0);
        assert(c.nSites > 0);

        // All tensors should belong to the same site universe.
        assert(a.nSites == b.nSites);
        assert(c.nSites == a.nSites);

        // Output must be a separate tensor. The matrix-matrix kernels are not
        // safe for in-place use, and the optimized kernels use non-aliasing
        // local pointers.
        assert(&c != &a);
        assert(&c != &b);

        // A clade can only be multiplied entrywise by another clade.
        assert(!a.isClade || b.isClade);

        const bool expectedIsClade = a.isClade ? true : b.isClade;
        const bool expectedMultisite = a.multisite || b.multisite;

        const std::size_t n = siteCountFromN(a.nSites);
        const std::size_t expectedSize =
            expectedIsClade
                ? (expectedMultisite ? 4 * n : 4)
                : (expectedMultisite ? 16 * n : 16);

        assert(c.isClade == expectedIsClade);
        assert(c.multisite == expectedMultisite);
        assert(c.data.size() == expectedSize);
    }

public:
    PartialLikelihoodTensor(bool _isClade, bool _multisite, int _nSites)
        : isClade(_isClade), multisite(_multisite), nSites(_nSites)
    {
        assert(nSites > 0);
        const std::size_t n = siteCountFromN(nSites);

        if (isClade) {
            data.resize(multisite ? 4 * n : 4);
        } else {
            data.resize(multisite ? 16 * n : 16);
        }
    }

    void fill(const Eigen::Matrix<Scalar, 4, 1>& mat)
    {
        assert(isClade);

        const Scalar x0 = mat(0, 0);
        const Scalar x1 = mat(1, 0);
        const Scalar x2 = mat(2, 0);
        const Scalar x3 = mat(3, 0);

        if (multisite) {
            const std::size_t n = siteCountFromN(nSites);
            assert(data.size() == 4 * n);

            Scalar* PLT_RESTRICT out = data.data();
            for (std::size_t p = 0; p < n; ++p) {
                Scalar* PLT_RESTRICT v = out + 4 * p;
                v[0] = x0;
                v[1] = x1;
                v[2] = x2;
                v[3] = x3;
            }
        } else {
            assert(data.size() == 4);
            data[0] = x0;
            data[1] = x1;
            data[2] = x2;
            data[3] = x3;
        }
    }

    void fill(const Eigen::Matrix<Scalar, 4, 4>& mat)
    {
        assert(!isClade);

        const Scalar m00 = mat(0, 0), m01 = mat(0, 1), m02 = mat(0, 2), m03 = mat(0, 3);
        const Scalar m10 = mat(1, 0), m11 = mat(1, 1), m12 = mat(1, 2), m13 = mat(1, 3);
        const Scalar m20 = mat(2, 0), m21 = mat(2, 1), m22 = mat(2, 2), m23 = mat(2, 3);
        const Scalar m30 = mat(3, 0), m31 = mat(3, 1), m32 = mat(3, 2), m33 = mat(3, 3);

        if (multisite) {
            const std::size_t n = siteCountFromN(nSites);
            assert(data.size() == 16 * n);

            Scalar* PLT_RESTRICT out = data.data();
            for (std::size_t p = 0; p < n; ++p) {
                Scalar* PLT_RESTRICT M = out + 16 * p;

                M[0]  = m00; M[1]  = m01; M[2]  = m02; M[3]  = m03;
                M[4]  = m10; M[5]  = m11; M[6]  = m12; M[7]  = m13;
                M[8]  = m20; M[9]  = m21; M[10] = m22; M[11] = m23;
                M[12] = m30; M[13] = m31; M[14] = m32; M[15] = m33;
            }
        } else {
            assert(data.size() == 16);

            data[0]  = m00; data[1]  = m01; data[2]  = m02; data[3]  = m03;
            data[4]  = m10; data[5]  = m11; data[6]  = m12; data[7]  = m13;
            data[8]  = m20; data[9]  = m21; data[10] = m22; data[11] = m23;
            data[12] = m30; data[13] = m31; data[14] = m32; data[15] = m33;
        }
    }

    friend void product(
        const PartialLikelihoodTensor& a,
        const PartialLikelihoodTensor& b,
        PartialLikelihoodTensor& c)
    {
        assertProductShape(a, b, c);

        if (a.isClade) {
            assert(b.isClade);

            if (!a.multisite) {
                if (!b.multisite) {
                    vec_vec(a.data, b.data, c.data);
                } else {
                    vec_multivec(a.data, b.data, c.data);
                }
            } else {
                if (!b.multisite) {
                    // Entrywise multiplication commutes.
                    vec_multivec(b.data, a.data, c.data);
                } else {
                    multivec_multivec(a.data, b.data, c.data);
                }
            }
            return;
        }

        // From here on, a is a matrix/segment.
        if (!a.multisite) {
            if (b.isClade) {
                if (!b.multisite) {
                    mat_vec(a.data, b.data, c.data);
                } else {
                    mat_multivec(a.data, b.data, c.data);
                }
            } else {
                if (!b.multisite) {
                    mat_mat(a.data, b.data, c.data);
                } else {
                    mat_multimat(a.data, b.data, c.data);
                }
            }
            return;
        }

        // From here on, a is a multisite matrix/segment.
        if (b.isClade) {
            if (!b.multisite) {
                multimat_vec(a.data, b.data, c.data);
            } else {
                multimat_multivec(a.data, b.data, c.data);
            }
        } else {
            if (!b.multisite) {
                multimat_mat(a.data, b.data, c.data);
            } else {
                multimat_multimat(a.data, b.data, c.data);
            }
        }
    }

    // ============================================================
    //
    // Site-contiguous layouts:
    //
    //   single vector:      A[i]                  size 4
    //   multisite vector:  A[4*p + i]            size 4*s
    //
    //   single matrix:      A[4*i + j]            size 16
    //   multisite matrix:  A[16*p + 4*i + j]     size 16*s
    //
    // The kernels below work site-by-site. For multisite objects, each site
    // occupies one contiguous block of 4 or 16 Scalars.
    //
    // ============================================================

    static inline void vec_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 4 && B.size() == 4 && C.size() == 4);
        assertNoAlias(A, B, C);

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar* PLT_RESTRICT Bp = B.data();
        Scalar* PLT_RESTRICT Cp = C.data();

        Cp[0] = Ap[0] * Bp[0];
        Cp[1] = Ap[1] * Bp[1];
        Cp[2] = Ap[2] * Bp[2];
        Cp[3] = Ap[3] * Bp[3];
    }

    static inline void vec_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 4);
        assert(B.size() % 4 == 0);
        assert(C.size() == B.size());
        assertNoAlias(A, B, C);

        const std::size_t s = B.size() / 4;

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar a0 = Ap[0];
        const Scalar a1 = Ap[1];
        const Scalar a2 = Ap[2];
        const Scalar a3 = Ap[3];

        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Bp = Bbase + 4 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 4 * p;

            Cp[0] = a0 * Bp[0];
            Cp[1] = a1 * Bp[1];
            Cp[2] = a2 * Bp[2];
            Cp[3] = a3 * Bp[3];
        }
    }

    static inline void multivec_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 4 == 0);
        assert(B.size() == A.size());
        assert(C.size() == A.size());
        assertNoAlias(A, B, C);

        const std::size_t s = A.size() / 4;

        const Scalar* PLT_RESTRICT Abase = A.data();
        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Ap = Abase + 4 * p;
            const Scalar* PLT_RESTRICT Bp = Bbase + 4 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 4 * p;

            Cp[0] = Ap[0] * Bp[0];
            Cp[1] = Ap[1] * Bp[1];
            Cp[2] = Ap[2] * Bp[2];
            Cp[3] = Ap[3] * Bp[3];
        }
    }

    static inline void mat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() == 4 && C.size() == 4);
        assertNoAlias(A, B, C);

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar* PLT_RESTRICT Bp = B.data();
        Scalar* PLT_RESTRICT Cp = C.data();

        const Scalar b0 = Bp[0];
        const Scalar b1 = Bp[1];
        const Scalar b2 = Bp[2];
        const Scalar b3 = Bp[3];

        Cp[0] = Ap[0]  * b0 + Ap[1]  * b1 + Ap[2]  * b2 + Ap[3]  * b3;
        Cp[1] = Ap[4]  * b0 + Ap[5]  * b1 + Ap[6]  * b2 + Ap[7]  * b3;
        Cp[2] = Ap[8]  * b0 + Ap[9]  * b1 + Ap[10] * b2 + Ap[11] * b3;
        Cp[3] = Ap[12] * b0 + Ap[13] * b1 + Ap[14] * b2 + Ap[15] * b3;
    }

    static inline void mat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16);
        assert(B.size() % 4 == 0);
        assert(C.size() == B.size());
        assertNoAlias(A, B, C);

        const std::size_t s = B.size() / 4;

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Bp = Bbase + 4 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 4 * p;

            const Scalar b0 = Bp[0];
            const Scalar b1 = Bp[1];
            const Scalar b2 = Bp[2];
            const Scalar b3 = Bp[3];

            Cp[0] = Ap[0]  * b0 + Ap[1]  * b1 + Ap[2]  * b2 + Ap[3]  * b3;
            Cp[1] = Ap[4]  * b0 + Ap[5]  * b1 + Ap[6]  * b2 + Ap[7]  * b3;
            Cp[2] = Ap[8]  * b0 + Ap[9]  * b1 + Ap[10] * b2 + Ap[11] * b3;
            Cp[3] = Ap[12] * b0 + Ap[13] * b1 + Ap[14] * b2 + Ap[15] * b3;
        }
    }

    static inline void mat_mat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() == 16 && C.size() == 16);
        assertNoAlias(A, B, C);

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar* PLT_RESTRICT Bp = B.data();
        Scalar* PLT_RESTRICT Cp = C.data();

        for (int i = 0; i < 4; ++i) {
            const Scalar a0 = Ap[4*i + 0];
            const Scalar a1 = Ap[4*i + 1];
            const Scalar a2 = Ap[4*i + 2];
            const Scalar a3 = Ap[4*i + 3];

            Cp[4*i + 0] = a0 * Bp[0]  + a1 * Bp[4]  + a2 * Bp[8]  + a3 * Bp[12];
            Cp[4*i + 1] = a0 * Bp[1]  + a1 * Bp[5]  + a2 * Bp[9]  + a3 * Bp[13];
            Cp[4*i + 2] = a0 * Bp[2]  + a1 * Bp[6]  + a2 * Bp[10] + a3 * Bp[14];
            Cp[4*i + 3] = a0 * Bp[3]  + a1 * Bp[7]  + a2 * Bp[11] + a3 * Bp[15];
        }
    }

    static inline void multimat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == 4);
        assert(C.size() == A.size() / 4);
        assertNoAlias(A, B, C);

        const std::size_t s = A.size() / 16;

        const Scalar* PLT_RESTRICT Abase = A.data();
        const Scalar* PLT_RESTRICT Bp = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        const Scalar b0 = Bp[0];
        const Scalar b1 = Bp[1];
        const Scalar b2 = Bp[2];
        const Scalar b3 = Bp[3];

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Ap = Abase + 16 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 4 * p;

            Cp[0] = Ap[0]  * b0 + Ap[1]  * b1 + Ap[2]  * b2 + Ap[3]  * b3;
            Cp[1] = Ap[4]  * b0 + Ap[5]  * b1 + Ap[6]  * b2 + Ap[7]  * b3;
            Cp[2] = Ap[8]  * b0 + Ap[9]  * b1 + Ap[10] * b2 + Ap[11] * b3;
            Cp[3] = Ap[12] * b0 + Ap[13] * b1 + Ap[14] * b2 + Ap[15] * b3;
        }
    }

    static inline void multimat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);

        const std::size_t s = A.size() / 16;

        assert(B.size() == 4 * s);
        assert(C.size() == 4 * s);
        assertNoAlias(A, B, C);

        const Scalar* PLT_RESTRICT Abase = A.data();
        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Ap = Abase + 16 * p;
            const Scalar* PLT_RESTRICT Bp = Bbase + 4 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 4 * p;

            const Scalar b0 = Bp[0];
            const Scalar b1 = Bp[1];
            const Scalar b2 = Bp[2];
            const Scalar b3 = Bp[3];

            Cp[0] = Ap[0]  * b0 + Ap[1]  * b1 + Ap[2]  * b2 + Ap[3]  * b3;
            Cp[1] = Ap[4]  * b0 + Ap[5]  * b1 + Ap[6]  * b2 + Ap[7]  * b3;
            Cp[2] = Ap[8]  * b0 + Ap[9]  * b1 + Ap[10] * b2 + Ap[11] * b3;
            Cp[3] = Ap[12] * b0 + Ap[13] * b1 + Ap[14] * b2 + Ap[15] * b3;
        }
    }

    static inline void mat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16);
        assert(B.size() % 16 == 0);
        assert(C.size() == B.size());
        assertNoAlias(A, B, C);

        const std::size_t s = B.size() / 16;

        const Scalar* PLT_RESTRICT Ap = A.data();
        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Bp = Bbase + 16 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 16 * p;

            for (int i = 0; i < 4; ++i) {
                const Scalar a0 = Ap[4*i + 0];
                const Scalar a1 = Ap[4*i + 1];
                const Scalar a2 = Ap[4*i + 2];
                const Scalar a3 = Ap[4*i + 3];

                Cp[4*i + 0] = a0 * Bp[0]  + a1 * Bp[4]  + a2 * Bp[8]  + a3 * Bp[12];
                Cp[4*i + 1] = a0 * Bp[1]  + a1 * Bp[5]  + a2 * Bp[9]  + a3 * Bp[13];
                Cp[4*i + 2] = a0 * Bp[2]  + a1 * Bp[6]  + a2 * Bp[10] + a3 * Bp[14];
                Cp[4*i + 3] = a0 * Bp[3]  + a1 * Bp[7]  + a2 * Bp[11] + a3 * Bp[15];
            }
        }
    }

    static inline void multimat_mat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == 16);
        assert(C.size() == A.size());
        assertNoAlias(A, B, C);

        const std::size_t s = A.size() / 16;

        const Scalar* PLT_RESTRICT Abase = A.data();
        const Scalar* PLT_RESTRICT Bp = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Ap = Abase + 16 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 16 * p;

            for (int i = 0; i < 4; ++i) {
                const Scalar a0 = Ap[4*i + 0];
                const Scalar a1 = Ap[4*i + 1];
                const Scalar a2 = Ap[4*i + 2];
                const Scalar a3 = Ap[4*i + 3];

                Cp[4*i + 0] = a0 * Bp[0]  + a1 * Bp[4]  + a2 * Bp[8]  + a3 * Bp[12];
                Cp[4*i + 1] = a0 * Bp[1]  + a1 * Bp[5]  + a2 * Bp[9]  + a3 * Bp[13];
                Cp[4*i + 2] = a0 * Bp[2]  + a1 * Bp[6]  + a2 * Bp[10] + a3 * Bp[14];
                Cp[4*i + 3] = a0 * Bp[3]  + a1 * Bp[7]  + a2 * Bp[11] + a3 * Bp[15];
            }
        }
    }

    static inline void multimat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == A.size());
        assert(C.size() == A.size());
        assertNoAlias(A, B, C);

        const std::size_t s = A.size() / 16;

        const Scalar* PLT_RESTRICT Abase = A.data();
        const Scalar* PLT_RESTRICT Bbase = B.data();
        Scalar* PLT_RESTRICT Cbase = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* PLT_RESTRICT Ap = Abase + 16 * p;
            const Scalar* PLT_RESTRICT Bp = Bbase + 16 * p;
            Scalar* PLT_RESTRICT Cp = Cbase + 16 * p;

            for (int i = 0; i < 4; ++i) {
                const Scalar a0 = Ap[4*i + 0];
                const Scalar a1 = Ap[4*i + 1];
                const Scalar a2 = Ap[4*i + 2];
                const Scalar a3 = Ap[4*i + 3];

                Cp[4*i + 0] = a0 * Bp[0]  + a1 * Bp[4]  + a2 * Bp[8]  + a3 * Bp[12];
                Cp[4*i + 1] = a0 * Bp[1]  + a1 * Bp[5]  + a2 * Bp[9]  + a3 * Bp[13];
                Cp[4*i + 2] = a0 * Bp[2]  + a1 * Bp[6]  + a2 * Bp[10] + a3 * Bp[14];
                Cp[4*i + 3] = a0 * Bp[3]  + a1 * Bp[7]  + a2 * Bp[11] + a3 * Bp[15];
            }
        }
    }
};

#undef PLT_RESTRICT

#endif // LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H
