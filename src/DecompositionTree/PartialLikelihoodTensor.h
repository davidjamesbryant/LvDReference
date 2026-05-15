//
// Created by David Bryant on 14/05/2026.
// Optimized pass: raw pointers, site-local offsets, and unrolled 4-entry kernels.
//

#ifndef LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H
#define LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H


// tiny4_tensor.hpp
#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <Eigen>

#include "phylib.h"


class PartialLikelihoodTensor
{
private:
    std::vector<Scalar> data;
    bool isClade;
    bool multisite;
    int nSites;
public:
    PartialLikelihoodTensor(bool _isClade, bool _multisite, int _nSites): isClade(_isClade), multisite(_multisite), nSites(_nSites)
    {
        if (isClade)
        {
            if (multisite)
                data.resize(4 * nSites);
            else
                data.resize(4);
        } else
        {
            if (multisite)
                data.resize(16*nSites);
            else
                data.resize(16);
        }
    }

    void fill(Eigen::Matrix<Scalar, 4, 1> mat)
    {
        assert(isClade);
        if (multisite)
        {
            for (std::size_t s=0;s<nSites; s++)
            {
                for (std::size_t i=0;i<4;i++)
                    data[s * 4 + i] = mat(i,0);
            }
        }
        else
        {
            for (std::size_t i=0;i<4;i++)
                data[i] = mat(i,0);
        }
    }

    void fill(Eigen::Matrix<Scalar, 4, 4> mat)
    {
        assert(!isClade);
        if (multisite)
        {
            for (std::size_t s=0;s<nSites; s++)
            {
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*s + 4*i + j] = mat(i,j);
            }
        } else
        {
            for (std::size_t i=0;i<4;i++)
                for (std::size_t j=0;j<4;j++)
                    data[4*i + j] = mat(i,j);
        }
    }

    friend void product(const PartialLikelihoodTensor& a, const PartialLikelihoodTensor& b, PartialLikelihoodTensor& c)
    {
        //If a is a clade then b has to be too. But either can be multisite.
        if (a.isClade && !a.multisite) //a is a clade, same for all sites
        {
            assert(b.isClade);
            if (!b.multisite) //b is a clade, same for all sites
                vec_vec(a.data, b.data, c.data);
            else  //b is a clade, different for all sites
                vec_multivec(a.data, b.data, c.data);
        }
        if (a.isClade && !!a.multisite) //a is a clade, different for all sites
        {
            assert(b.isClade);
            if (!b.multisite) //b is a clade, same for all sites
                vec_multivec(b.data, a.data, c.data); //Note - entrywise multiplication commutes
            else //b is a clade, different for all sites
                multivec_multivec(a.data, b.data, c.data);
        }

        //If a is a matrix, then b can be a clade or a matrix. Either can be !multisite
        if (!a.isClade && !a.multisite)
        {
            //a is a segment, same for all sites
            if (b.isClade && !b.multisite) //b is a clade, same for all sites
                mat_vec(a.data, b.data, c.data);
            else if (b.isClade && !!b.multisite) //b is a clade, different for all sites
                mat_multivec(a.data, b.data, c.data);
            else if (!b.isClade && !b.multisite) //b is a segment, same for all sites
                mat_mat(a.data, b.data, c.data);
            else if (!b.isClade && !!b.multisite) //b is a segment, different for all sites
                mat_multimat(a.data, b.data, c.data);
        }

        if (!a.isClade && !!a.multisite)
        {
            //a is a segment, different for all sites
            if (b.isClade && !b.multisite) //b is a clade, same for all sites
                multimat_vec(a.data, b.data, c.data);
            else if (b.isClade && !!b.multisite) //b is a clade, different for all sites
                multimat_multivec(a.data, b.data, c.data);
            else if (!b.isClade && !b.multisite) //b is a segment, same for all sites
                multimat_mat(a.data, b.data, c.data);
            else if (!b.isClade && !!b.multisite) //b is a segment, different for all sites
                multimat_multimat(a.data, b.data, c.data);
        }
    }

private:
    // ============================================================
    //
    // In the following:
    //      A multivec is a 4xs matrix, each column being a vector
    //      A multiMat is a 4x(4s) matrix, each 4 columns making up a 4x4 matrix
    //
    // ============================================================



    // --------------------------------------------------------
    // Vector-vector element-wise multiply
    //
    // Computes:
    //
    //   C[i]
    //      = A[i] B[i]
    //
    // A size: 4
    // B size: 4
    // C size: 4
    // --------------------------------------------------------
    static inline void vec_vec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size()==4 && B.size()==4 && C.size()==4);

        const Scalar* Ap = A.data();
        const Scalar* Bp = B.data();
        Scalar* Cp = C.data();

        Cp[0] = Ap[0] * Bp[0];
        Cp[1] = Ap[1] * Bp[1];
        Cp[2] = Ap[2] * Bp[2];
        Cp[3] = Ap[3] * Bp[3];
    }

    // --------------------------------------------------------
    // Vector-multivector element-wise multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = A[i] B[p][i]
    //
    // A size: 4
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void vec_multivec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 4);
        assert(B.size() % 4 == 0);
        assert(C.size() == B.size());

        const std::size_t s = B.size() / 4;
        const Scalar* Ap = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Bp = Bdata + 4*p;
            Scalar* Cp = Cdata + 4*p;

            Cp[0] = Ap[0] * Bp[0];
            Cp[1] = Ap[1] * Bp[1];
            Cp[2] = Ap[2] * Bp[2];
            Cp[3] = Ap[3] * Bp[3];
        }
    }

    // --------------------------------------------------------
    // Multivector-multivector element-wise multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = A[i] B[p][i]
    //
    // A size: 4
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void multivec_multivec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 4 ==  0);
        assert(B.size() == A.size());
        assert(C.size() == B.size());

        const std::size_t s = A.size() / 4;
        const Scalar* Adata = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Ap = Adata + 4*p;
            const Scalar* Bp = Bdata + 4*p;
            Scalar* Cp = Cdata + 4*p;

            Cp[0] = Ap[0] * Bp[0];
            Cp[1] = Ap[1] * Bp[1];
            Cp[2] = Ap[2] * Bp[2];
            Cp[3] = Ap[3] * Bp[3];
        }
    }


    // --------------------------------------------------------
    // Matrix-vector multiply
    //
    // Computes:
    //
    //   C[i]
    //      = sum_k A[i][k] B[k]
    //
    // A size: 16
    // B size: 4
    // C size: 4
    // --------------------------------------------------------

    static inline void mat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C
    )
    {
        assert(A.size() == 16 && B.size() == 4 && C.size() == 4);

        const Scalar* Ap = A.data();
        const Scalar* Bp = B.data();
        Scalar* Cp = C.data();

        Cp[0] = Ap[0]*Bp[0] + Ap[1]*Bp[1] + Ap[2]*Bp[2] + Ap[3]*Bp[3];
        Cp[1] = Ap[4]*Bp[0] + Ap[5]*Bp[1] + Ap[6]*Bp[2] + Ap[7]*Bp[3];
        Cp[2] = Ap[8]*Bp[0] + Ap[9]*Bp[1] + Ap[10]*Bp[2] + Ap[11]*Bp[3];
        Cp[3] = Ap[12]*Bp[0] + Ap[13]*Bp[1] + Ap[14]*Bp[2] + Ap[15]*Bp[3];
    }

    // --------------------------------------------------------
    // Matrix-multivector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[i][k] B[p][k]
    //
    // A size: 16
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------


    static inline void mat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() % 4 == 0 && C.size() == B.size());

        const std::size_t s = B.size() / 4;
        const Scalar* Ap = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Bp = Bdata + 4*p;
            Scalar* Cp = Cdata + 4*p;

            Cp[0] = Ap[0]*Bp[0] + Ap[1]*Bp[1] + Ap[2]*Bp[2] + Ap[3]*Bp[3];
            Cp[1] = Ap[4]*Bp[0] + Ap[5]*Bp[1] + Ap[6]*Bp[2] + Ap[7]*Bp[3];
            Cp[2] = Ap[8]*Bp[0] + Ap[9]*Bp[1] + Ap[10]*Bp[2] + Ap[11]*Bp[3];
            Cp[3] = Ap[12]*Bp[0] + Ap[13]*Bp[1] + Ap[14]*Bp[2] + Ap[15]*Bp[3];
        }
    }


    // --------------------------------------------------------
    // Matrix-matrix multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[i][k] B[p][k]
    //
    // A size: 16
    // B size: 16
    // C size: 16
    // --------------------------------------------------------

    static inline void mat_mat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() == 16 && C.size() == 16);

        const Scalar* Ap = A.data();
        const Scalar* Bp = B.data();
        Scalar* Cp = C.data();

        Cp[0]  = Ap[0]*Bp[0]  + Ap[1]*Bp[4]  + Ap[2]*Bp[8]  + Ap[3]*Bp[12];
        Cp[1]  = Ap[0]*Bp[1]  + Ap[1]*Bp[5]  + Ap[2]*Bp[9]  + Ap[3]*Bp[13];
        Cp[2]  = Ap[0]*Bp[2]  + Ap[1]*Bp[6]  + Ap[2]*Bp[10] + Ap[3]*Bp[14];
        Cp[3]  = Ap[0]*Bp[3]  + Ap[1]*Bp[7]  + Ap[2]*Bp[11] + Ap[3]*Bp[15];

        Cp[4]  = Ap[4]*Bp[0]  + Ap[5]*Bp[4]  + Ap[6]*Bp[8]  + Ap[7]*Bp[12];
        Cp[5]  = Ap[4]*Bp[1]  + Ap[5]*Bp[5]  + Ap[6]*Bp[9]  + Ap[7]*Bp[13];
        Cp[6]  = Ap[4]*Bp[2]  + Ap[5]*Bp[6]  + Ap[6]*Bp[10] + Ap[7]*Bp[14];
        Cp[7]  = Ap[4]*Bp[3]  + Ap[5]*Bp[7]  + Ap[6]*Bp[11] + Ap[7]*Bp[15];

        Cp[8]  = Ap[8]*Bp[0]  + Ap[9]*Bp[4]  + Ap[10]*Bp[8]  + Ap[11]*Bp[12];
        Cp[9]  = Ap[8]*Bp[1]  + Ap[9]*Bp[5]  + Ap[10]*Bp[9]  + Ap[11]*Bp[13];
        Cp[10] = Ap[8]*Bp[2]  + Ap[9]*Bp[6]  + Ap[10]*Bp[10] + Ap[11]*Bp[14];
        Cp[11] = Ap[8]*Bp[3]  + Ap[9]*Bp[7]  + Ap[10]*Bp[11] + Ap[11]*Bp[15];

        Cp[12] = Ap[12]*Bp[0] + Ap[13]*Bp[4] + Ap[14]*Bp[8]  + Ap[15]*Bp[12];
        Cp[13] = Ap[12]*Bp[1] + Ap[13]*Bp[5] + Ap[14]*Bp[9]  + Ap[15]*Bp[13];
        Cp[14] = Ap[12]*Bp[2] + Ap[13]*Bp[6] + Ap[14]*Bp[10] + Ap[15]*Bp[14];
        Cp[15] = Ap[12]*Bp[3] + Ap[13]*Bp[7] + Ap[14]*Bp[11] + Ap[15]*Bp[15];
    }


    // --------------------------------------------------------
    // Multimatrix-vector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[p][i][k] B[k]
    //
    // A size: 16*s
    // B size: 4
    // C size: 4*s
    // --------------------------------------------------------
    static inline void multimat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0 && B.size() == 4 && C.size() == A.size()/4);

        const std::size_t s = A.size() / 16;
        const Scalar* Adata = A.data();
        const Scalar* Bp = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Ap = Adata + 16*p;
            Scalar* Cp = Cdata + 4*p;

            Cp[0] = Ap[0]*Bp[0] + Ap[1]*Bp[1] + Ap[2]*Bp[2] + Ap[3]*Bp[3];
            Cp[1] = Ap[4]*Bp[0] + Ap[5]*Bp[1] + Ap[6]*Bp[2] + Ap[7]*Bp[3];
            Cp[2] = Ap[8]*Bp[0] + Ap[9]*Bp[1] + Ap[10]*Bp[2] + Ap[11]*Bp[3];
            Cp[3] = Ap[12]*Bp[0] + Ap[13]*Bp[1] + Ap[14]*Bp[2] + Ap[15]*Bp[3];
        }
    }



    // --------------------------------------------------------
    // Multimatrix-multivector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[p][i][k] B[p][k]
    //
    // A size: 16*s
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void multimat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);

        const std::size_t nSites = A.size() / 16;

        assert(B.size() == 4 * nSites);
        assert(C.size() == 4 * nSites);

        const Scalar* Adata = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < nSites; ++p)
        {
            const Scalar* Ap = Adata + 16*p;
            const Scalar* Bp = Bdata + 4*p;
            Scalar* Cp = Cdata + 4*p;

            Cp[0] = Ap[0]*Bp[0] + Ap[1]*Bp[1] + Ap[2]*Bp[2] + Ap[3]*Bp[3];
            Cp[1] = Ap[4]*Bp[0] + Ap[5]*Bp[1] + Ap[6]*Bp[2] + Ap[7]*Bp[3];
            Cp[2] = Ap[8]*Bp[0] + Ap[9]*Bp[1] + Ap[10]*Bp[2] + Ap[11]*Bp[3];
            Cp[3] = Ap[12]*Bp[0] + Ap[13]*Bp[1] + Ap[14]*Bp[2] + Ap[15]*Bp[3];
        }
    }


    // --------------------------------------------------------
    // Matrix-multimatrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[i][k] B[p][k][j]
    //
    // A size: 16
    // B size: 16*s
    // C size: 16*s
    // --------------------------------------------------------

    static inline void mat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16);
        assert(B.size() % 16 == 0);
        std::size_t s = B.size() / 16;
        assert(C.size() == B.size());

        const Scalar* Ap = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Bp = Bdata + 16*p;
            Scalar* Cp = Cdata + 16*p;

            Cp[0]  = Ap[0]*Bp[0]  + Ap[1]*Bp[4]  + Ap[2]*Bp[8]  + Ap[3]*Bp[12];
            Cp[1]  = Ap[0]*Bp[1]  + Ap[1]*Bp[5]  + Ap[2]*Bp[9]  + Ap[3]*Bp[13];
            Cp[2]  = Ap[0]*Bp[2]  + Ap[1]*Bp[6]  + Ap[2]*Bp[10] + Ap[3]*Bp[14];
            Cp[3]  = Ap[0]*Bp[3]  + Ap[1]*Bp[7]  + Ap[2]*Bp[11] + Ap[3]*Bp[15];

            Cp[4]  = Ap[4]*Bp[0]  + Ap[5]*Bp[4]  + Ap[6]*Bp[8]  + Ap[7]*Bp[12];
            Cp[5]  = Ap[4]*Bp[1]  + Ap[5]*Bp[5]  + Ap[6]*Bp[9]  + Ap[7]*Bp[13];
            Cp[6]  = Ap[4]*Bp[2]  + Ap[5]*Bp[6]  + Ap[6]*Bp[10] + Ap[7]*Bp[14];
            Cp[7]  = Ap[4]*Bp[3]  + Ap[5]*Bp[7]  + Ap[6]*Bp[11] + Ap[7]*Bp[15];

            Cp[8]  = Ap[8]*Bp[0]  + Ap[9]*Bp[4]  + Ap[10]*Bp[8]  + Ap[11]*Bp[12];
            Cp[9]  = Ap[8]*Bp[1]  + Ap[9]*Bp[5]  + Ap[10]*Bp[9]  + Ap[11]*Bp[13];
            Cp[10] = Ap[8]*Bp[2]  + Ap[9]*Bp[6]  + Ap[10]*Bp[10] + Ap[11]*Bp[14];
            Cp[11] = Ap[8]*Bp[3]  + Ap[9]*Bp[7]  + Ap[10]*Bp[11] + Ap[11]*Bp[15];

            Cp[12] = Ap[12]*Bp[0] + Ap[13]*Bp[4] + Ap[14]*Bp[8]  + Ap[15]*Bp[12];
            Cp[13] = Ap[12]*Bp[1] + Ap[13]*Bp[5] + Ap[14]*Bp[9]  + Ap[15]*Bp[13];
            Cp[14] = Ap[12]*Bp[2] + Ap[13]*Bp[6] + Ap[14]*Bp[10] + Ap[15]*Bp[14];
            Cp[15] = Ap[12]*Bp[3] + Ap[13]*Bp[7] + Ap[14]*Bp[11] + Ap[15]*Bp[15];
        }
    }


    // --------------------------------------------------------
    // Multimatrix-matrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[p][i][k] B[k][j]
    //
    // A size: 16*s
    // B size: 16
    // C size: 16*s
    // --------------------------------------------------------
    static inline void multimat_mat(
    const std::vector<Scalar>& A,
    const std::vector<Scalar>& B,
    std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == 16);
        assert(C.size() == A.size());

        const std::size_t s = A.size() / 16;
        const Scalar* Adata = A.data();
        const Scalar* Bp = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* Ap = Adata + 16*p;
            Scalar* Cp = Cdata + 16*p;

            Cp[0]  = Ap[0]*Bp[0]  + Ap[1]*Bp[4]  + Ap[2]*Bp[8]  + Ap[3]*Bp[12];
            Cp[1]  = Ap[0]*Bp[1]  + Ap[1]*Bp[5]  + Ap[2]*Bp[9]  + Ap[3]*Bp[13];
            Cp[2]  = Ap[0]*Bp[2]  + Ap[1]*Bp[6]  + Ap[2]*Bp[10] + Ap[3]*Bp[14];
            Cp[3]  = Ap[0]*Bp[3]  + Ap[1]*Bp[7]  + Ap[2]*Bp[11] + Ap[3]*Bp[15];

            Cp[4]  = Ap[4]*Bp[0]  + Ap[5]*Bp[4]  + Ap[6]*Bp[8]  + Ap[7]*Bp[12];
            Cp[5]  = Ap[4]*Bp[1]  + Ap[5]*Bp[5]  + Ap[6]*Bp[9]  + Ap[7]*Bp[13];
            Cp[6]  = Ap[4]*Bp[2]  + Ap[5]*Bp[6]  + Ap[6]*Bp[10] + Ap[7]*Bp[14];
            Cp[7]  = Ap[4]*Bp[3]  + Ap[5]*Bp[7]  + Ap[6]*Bp[11] + Ap[7]*Bp[15];

            Cp[8]  = Ap[8]*Bp[0]  + Ap[9]*Bp[4]  + Ap[10]*Bp[8]  + Ap[11]*Bp[12];
            Cp[9]  = Ap[8]*Bp[1]  + Ap[9]*Bp[5]  + Ap[10]*Bp[9]  + Ap[11]*Bp[13];
            Cp[10] = Ap[8]*Bp[2]  + Ap[9]*Bp[6]  + Ap[10]*Bp[10] + Ap[11]*Bp[14];
            Cp[11] = Ap[8]*Bp[3]  + Ap[9]*Bp[7]  + Ap[10]*Bp[11] + Ap[11]*Bp[15];

            Cp[12] = Ap[12]*Bp[0] + Ap[13]*Bp[4] + Ap[14]*Bp[8]  + Ap[15]*Bp[12];
            Cp[13] = Ap[12]*Bp[1] + Ap[13]*Bp[5] + Ap[14]*Bp[9]  + Ap[15]*Bp[13];
            Cp[14] = Ap[12]*Bp[2] + Ap[13]*Bp[6] + Ap[14]*Bp[10] + Ap[15]*Bp[14];
            Cp[15] = Ap[12]*Bp[3] + Ap[13]*Bp[7] + Ap[14]*Bp[11] + Ap[15]*Bp[15];
        }
    }




    // --------------------------------------------------------
    // Multimatrix-multimatrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[p][i][k] B[p][k][j]
    //
    // A size: 16*s
    // B size: 16*s
    // C size: 16*s
    // --------------------------------------------------------

    static inline void multimat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == A.size() && C.size() == A.size());

        const std::size_t s = A.size() / 16;
        const Scalar* Adata = A.data();
        const Scalar* Bdata = B.data();
        Scalar* Cdata = C.data();

        for (std::size_t p = 0; p < s; ++p)
        {
            const Scalar* Ap = Adata + 16*p;
            const Scalar* Bp = Bdata + 16*p;
            Scalar* Cp = Cdata + 16*p;

            Cp[0]  = Ap[0]*Bp[0]  + Ap[1]*Bp[4]  + Ap[2]*Bp[8]  + Ap[3]*Bp[12];
            Cp[1]  = Ap[0]*Bp[1]  + Ap[1]*Bp[5]  + Ap[2]*Bp[9]  + Ap[3]*Bp[13];
            Cp[2]  = Ap[0]*Bp[2]  + Ap[1]*Bp[6]  + Ap[2]*Bp[10] + Ap[3]*Bp[14];
            Cp[3]  = Ap[0]*Bp[3]  + Ap[1]*Bp[7]  + Ap[2]*Bp[11] + Ap[3]*Bp[15];

            Cp[4]  = Ap[4]*Bp[0]  + Ap[5]*Bp[4]  + Ap[6]*Bp[8]  + Ap[7]*Bp[12];
            Cp[5]  = Ap[4]*Bp[1]  + Ap[5]*Bp[5]  + Ap[6]*Bp[9]  + Ap[7]*Bp[13];
            Cp[6]  = Ap[4]*Bp[2]  + Ap[5]*Bp[6]  + Ap[6]*Bp[10] + Ap[7]*Bp[14];
            Cp[7]  = Ap[4]*Bp[3]  + Ap[5]*Bp[7]  + Ap[6]*Bp[11] + Ap[7]*Bp[15];

            Cp[8]  = Ap[8]*Bp[0]  + Ap[9]*Bp[4]  + Ap[10]*Bp[8]  + Ap[11]*Bp[12];
            Cp[9]  = Ap[8]*Bp[1]  + Ap[9]*Bp[5]  + Ap[10]*Bp[9]  + Ap[11]*Bp[13];
            Cp[10] = Ap[8]*Bp[2]  + Ap[9]*Bp[6]  + Ap[10]*Bp[10] + Ap[11]*Bp[14];
            Cp[11] = Ap[8]*Bp[3]  + Ap[9]*Bp[7]  + Ap[10]*Bp[11] + Ap[11]*Bp[15];

            Cp[12] = Ap[12]*Bp[0] + Ap[13]*Bp[4] + Ap[14]*Bp[8]  + Ap[15]*Bp[12];
            Cp[13] = Ap[12]*Bp[1] + Ap[13]*Bp[5] + Ap[14]*Bp[9]  + Ap[15]*Bp[13];
            Cp[14] = Ap[12]*Bp[2] + Ap[13]*Bp[6] + Ap[14]*Bp[10] + Ap[15]*Bp[14];
            Cp[15] = Ap[12]*Bp[3] + Ap[13]*Bp[7] + Ap[14]*Bp[11] + Ap[15]*Bp[15];
        }
    }
};
#endif //LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H

/*****************
 * Here is the reference version of the code, before it was unrooled by the AI
 **************/

//
// Created by David Bryant on 14/05/2026.
//

#ifdef LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_REFERENCE_H


// tiny4_tensor.hpp
#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <Eigen>

#include "phylib.h"


class PartialLikelihoodTensor
{
private:
    std::vector<Scalar> data;
    bool isClade;
    bool multisite;
    int nSites;
public:
    PartialLikelihoodTensor(bool _isClade, bool _multisite, int _nSites): isClade(_isClade), multisite(_multisite), nSites(_nSites)
    {
        if (isClade)
        {
            if (multisite)
                data.resize(4 * nSites);
            else
                data.resize(4);
        } else
        {
            if (multisite)
                data.resize(16*nSites);
            else
                data.resize(16);
        }
    }

    void fill(Eigen::Matrix<Scalar, 4, 1> mat)
    {
        assert(isClade);
        if (multisite)
        {
            for (std::size_t s=0;s<nSites; s++)
            {
                for (std::size_t i=0;i<4;i++)
                    data[s * 4 + i] = mat(i,0);
            }
        }
        else
        {
            for (std::size_t i=0;i<4;i++)
                data[i] = mat(i,0);
        }
    }

    void fill(Eigen::Matrix<Scalar, 4, 4> mat)
    {
        assert(!isClade);
        if (multisite)
        {
            for (std::size_t s=0;s<nSites; s++)
            {
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*s + 4*i + j] = mat(i,j);
            }
        } else
        {
            for (std::size_t i=0;i<4;i++)
                for (std::size_t j=0;j<4;j++)
                    data[4*i + j] = mat(i,j);
        }
    }

    friend void product(const PartialLikelihoodTensor& a, const PartialLikelihoodTensor& b, PartialLikelihoodTensor& c)
    {
        //If a is a clade then b has to be too. But either can be multisite.
        if (a.isClade && !a.multisite) //a is a clade, same for all sites
        {
            assert(b.isClade);
            if (!b.multisite) //b is a clade, same for all sites
                vec_vec(a.data, b.data, c.data);
            else  //b is a clade, different for all sites
                vec_multivec(a.data, b.data, c.data);
        }
        if (a.isClade && !!a.multisite) //a is a clade, different for all sites
        {
            assert(b.isClade);
            if (!b.multisite) //b is a clade, same for all sites
                vec_multivec(b.data, a.data, c.data); //Note - entrywise multiplication commutes
            else //b is a clade, different for all sites
                multivec_multivec(a.data, b.data, c.data);
        }

        //If a is a matrix, then b can be a clade or a matrix. Either can be !multisite
        if (!a.isClade && !a.multisite)
        {
            //a is a segment, same for all sites
            if (b.isClade && !b.multisite) //b is a clade, same for all sites
                mat_vec(a.data, b.data, c.data);
            else if (b.isClade && !!b.multisite) //b is a clade, different for all sites
                mat_multivec(a.data, b.data, c.data);
            else if (!b.isClade && !b.multisite) //b is a segment, same for all sites
                mat_mat(a.data, b.data, c.data);
            else if (!b.isClade && !!b.multisite) //b is a segment, different for all sites
                mat_multimat(a.data, b.data, c.data);
        }

        if (!a.isClade && !!a.multisite)
        {
            //a is a segment, different for all sites
            if (b.isClade && !b.multisite) //b is a clade, same for all sites
                multimat_vec(a.data, b.data, c.data);
            else if (b.isClade && !!b.multisite) //b is a clade, different for all sites
                multimat_multivec(a.data, b.data, c.data);
            else if (!b.isClade && !b.multisite) //b is a segment, same for all sites
                multimat_mat(a.data, b.data, c.data);
            else if (!b.isClade && !!b.multisite) //b is a segment, different for all sites
                multimat_multimat(a.data, b.data, c.data);
        }
    }

private:
    // ============================================================
    //
    // In the following:
    //      A multivec is a 4xs matrix, each column being a vector
    //      A multiMat is a 4x(4s) matrix, each 4 columns making up a 4x4 matrix
    //
    // ============================================================



    // --------------------------------------------------------
    // Vector-vector element-wise multiply
    //
    // Computes:
    //
    //   C[i]
    //      = A[i] B[i]
    //
    // A size: 4
    // B size: 4
    // C size: 4
    // --------------------------------------------------------
    static inline void vec_vec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size()==4 && B.size()==4 && C.size()==4);
        for (std::size_t i=0;i<4;i++)
            C[i] = A[i]*B[i];
    }

    // --------------------------------------------------------
    // Vector-multivector element-wise multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = A[i] B[p][i]
    //
    // A size: 4
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void vec_multivec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 4);
        assert(B.size() % 4 == 0);
        assert(C.size() == B.size());

        const std::size_t s = B.size() / 4;

        for (std::size_t p = 0; p < s; ++p)
        {
            for (std::size_t i=0;i<4;i++)
                C[4*p + i] = A[i] * B[4*p + i];
        }
    }

    // --------------------------------------------------------
    // Multivector-multivector element-wise multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = A[i] B[p][i]
    //
    // A size: 4
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void multivec_multivec(const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 4 ==  0);
        assert(B.size() == A.size());
        assert(C.size() == B.size());

        const std::size_t s = A.size() / 4;

        for (std::size_t p = 0; p < s; ++p)
        {
            for (std::size_t i=0;i<4;i++)
                C[4*p + i] = A[4*p+i] * B[4*p + i];
        }
    }


    // --------------------------------------------------------
    // Matrix-vector multiply
    //
    // Computes:
    //
    //   C[i]
    //      = sum_k A[i][k] B[k]
    //
    // A size: 16
    // B size: 4
    // C size: 4
    // --------------------------------------------------------

    static inline void mat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C
    )
    {
        assert(A.size() == 16 && B.size() == 4 && C.size() == 4);
        for (std::size_t i=0;i<4;i++)
        {
            Scalar sum = 0.0;
            for (std::size_t j=0;j<4;j++)
                sum += A[4*i + j] * B[j];
            C[i] = sum;
        }
    }

    // --------------------------------------------------------
    // Matrix-multivector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[i][k] B[p][k]
    //
    // A size: 16
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------


    static inline void mat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() % 4 == 0 && C.size() == B.size());
        const std::size_t s = B.size() / 4;
        for (std::size_t p = 0; p < s; ++p)
        {
            for (std::size_t i=0;i<4;i++)
            {
                Scalar sum = 0.0;
                for (std::size_t j=0;j<4;j++)
                    sum += A[4*i + j] * B[4*p + j];
                C[4*p + i] = sum;
            }
        }
    }


    // --------------------------------------------------------
    // Matrix-matrix multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[i][k] B[p][k]
    //
    // A size: 16
    // B size: 16
    // C size: 16
    // --------------------------------------------------------

    static inline void mat_mat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16 && B.size() == 16 && C.size() == 16);
        for (std::size_t i=0;i<4;i++) {
            for (std::size_t j=0;j<4;j++) {
                Scalar sum = 0.0;
                for (std::size_t k=0;k<4;k++)
                    sum += A[4*i + k] * B[4*k + j];
                C[4*i + j] = sum;
            }
        }
    }


    // --------------------------------------------------------
    // Multimatrix-vector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[p][i][k] B[k]
    //
    // A size: 16*s
    // B size: 4
    // C size: 4*s
    // --------------------------------------------------------
    static inline void multimat_vec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0 && B.size() == 4 && C.size() == A.size()/4);
        const std::size_t s = A.size() / 16;
        for (std::size_t p = 0; p < s; ++p)
            for (std::size_t i=0;i<4;i++)
            {
                Scalar sum = 0.0;
                for (std::size_t j=0;j<4;j++)
                    sum += A[16*p + 4*i + j] * B[j];
                C[4*p + i] = sum;
            }
    }



    // --------------------------------------------------------
    // Multimatrix-multivector multiply
    //
    // Computes:
    //
    //   C[p][i]
    //      = sum_k A[p][i][k] B[p][k]
    //
    // A size: 16*s
    // B size: 4*s
    // C size: 4*s
    // --------------------------------------------------------

    static inline void multimat_multivec(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);

        const std::size_t nSites = A.size() / 16;

        assert(B.size() == 4 * nSites);
        assert(C.size() == 4 * nSites);


        for (std::size_t p = 0; p < nSites; ++p)
        {
            for (std::size_t i=0;i<4;i++)
            {
                Scalar sum = 0.0;
                for (std::size_t k=0;k<4;k++)
                    sum += A[16*p + 4*i + k] * B[4*p + k];
                C[4*p + i] = sum;
            }
        }
    }


    // --------------------------------------------------------
    // Matrix-multimatrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[i][k] B[p][k][j]
    //
    // A size: 16
    // B size: 16*s
    // C size: 16*s
    // --------------------------------------------------------

    static inline void mat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() == 16);
        assert(B.size() % 16 == 0);
        std::size_t s = B.size() / 16;
        assert(C.size() == B.size());

        for (std::size_t p = 0; p < s; ++p)
        {
            for (std::size_t i=0;i<4;i++)
                for (std::size_t j=0;j<4;j++)
                {
                    Scalar sum = 0.0;
                    for (std::size_t k=0;k<4;k++)
                        sum += A[4*i + k] * B[16*p + 4*k + j];
                    C[16*p + 4*i + j] = sum;
                }
        }
    }


    // --------------------------------------------------------
    // Multimatrix-matrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[p][i][k] B[k][j]
    //
    // A size: 16*s
    // B size: 16
    // C size: 16*s
    // --------------------------------------------------------
    static inline void multimat_mat(
    const std::vector<Scalar>& A,
    const std::vector<Scalar>& B,
    std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == 16);
        assert(C.size() == A.size());

        const std::size_t s = A.size() / 16;

        for (std::size_t p = 0; p < s; ++p) {
            const Scalar* Ap = A.data() + 16*p;
            Scalar* Cp = C.data() + 16*p;

            for (std::size_t i = 0; i < 4; ++i) {
                for (std::size_t j = 0; j < 4; ++j) {
                    Scalar sum = 0.0;
                    for (std::size_t k = 0; k < 4; ++k) {
                        sum += Ap[4*i + k] * B[4*k + j];
                    }
                    Cp[4*i + j] = sum;
                }
            }
        }
    }




    // --------------------------------------------------------
    // Multimatrix-multimatrix multiply
    //
    // Computes:
    //
    //   C[p][i][j]
    //      = sum_k A[p][i][k] B[p][k][j]
    //
    // A size: 16*s
    // B size: 16*s
    // C size: 16*s
    // --------------------------------------------------------

    static inline void multimat_multimat(
        const std::vector<Scalar>& A,
        const std::vector<Scalar>& B,
        std::vector<Scalar>& C)
    {
        assert(A.size() % 16 == 0);
        assert(B.size() == A.size() && C.size() == A.size());

        const std::size_t s = A.size() / 16;


        for (std::size_t p = 0; p < s; ++p)
        {
            for (std::size_t i=0;i<4;i++)
            {
                for (std::size_t j=0;j<4;j++)
                {
                    Scalar sum = 0.0;
                    for (std::size_t k=0;k<4;k++)
                        sum += A[16*p + 4*i + k] * B[16*p + 4*k + j];
                    C[16*p + 4*i + j] = sum;
                }
            }
        }
    }
};
#endif
