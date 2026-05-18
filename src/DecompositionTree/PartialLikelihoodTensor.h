// DO NOT EDIT — auto-generated from src/Reference/PartialLikelihoodTensorReference.h
// Edit the reference file and regenerate instead.
//
// Define USE_REFERENCE_TENSOR_CALCULATIONS (or set the CMake option
// USE_REFERENCE_TENSOR) to use the plain reference implementation instead.
//
// Optimizations applied in the default (non-reference) version:
//   (1) Raw pointer access: data.data() pointers obtained once per function call,
//       eliminating repeated vector-bounds indirection.
//   (2) Site-local offset pointers: inside the outer p-loop, pointers Ap/Bp/Cp are
//       computed once as data + 16*p (or 4*p), avoiding repeated offset arithmetic
//       in every inner-loop iteration.
//   (3) Inner loops of size 4 are unrolled in the computational kernels.
//       Simple setters/accessors are left as loops — the compiler handles those.

#ifndef LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H
#define LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H

#ifdef USE_REFERENCE_TENSOR_CALCULATIONS
#include "PartialLikelihoodTensorReference.h"
using PartialLikelihoodTensor = PartialLikelihoodTensorReference;
#else

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
    bool isVectors;  //True if this is 4xnsites, false if this is 4x4xnSites
    std::size_t size; //1 if this is a single vector or matrix, more if array
public:
    PartialLikelihoodTensor() : isVectors(false), size(0) {}

    PartialLikelihoodTensor(bool isVectors_, std::size_t size_) : isVectors(isVectors_), size(size_)
    {
        if (isVectors)
            data.resize(4*size);
        else
            data.resize(16*size);
    }

    /**
     * Change the type of storage and reallocate as necessary.
     */
    void resize(bool isVectors_, std::size_t size_)
    {
        isVectors = isVectors_;
        size      = size_;
        if (isVectors)
            data.resize(4*size);
        else
            data.resize(16*size);
    }

    /**
     * Set entry i in position p of this 4×size tensor.
     */
    void set(std::size_t p, std::size_t i, Scalar val)
    {
        assert(isVectors && p < size && i < 4);
        data[4*p + i] = val;
    }

    /**
     * Copy col into position p of this 4×size tensor.
     */
    void set_column(std::size_t p, const Eigen::Matrix<Scalar,4,1>& col)
    {
        assert(isVectors && p < size && col.size() == 4);
        Scalar* Cp = data.data() + 4*p;
        Cp[0] = col(0);  Cp[1] = col(1);  Cp[2] = col(2);  Cp[3] = col(3);
    }

    /**
     * Fill position p with repeated val.
     */
    void fill_column(std::size_t p, Scalar val)
    {
        assert(isVectors && p < size);
        Scalar* Cp = data.data() + 4*p;
        Cp[0] = val;  Cp[1] = val;  Cp[2] = val;  Cp[3] = val;
    }

    /**
     * Copy mat into this single 4×4 tensor (size must be 1).
     */
    void set(const Eigen::Matrix<Scalar,4,4>& mat)
    {
        assert(!isVectors && size == 1);
        Scalar* Dp = data.data();
        for (std::size_t i = 0; i < 4; i++)
            for (std::size_t j = 0; j < 4; j++)
                Dp[4*i + j] = mat(i,j);
    }

    // -----------------------------------------------------------------------
    // Computational kernels — all three optimisations applied.
    // -----------------------------------------------------------------------

    /**
     * this[p][i] = vecsA[p][i] * vecsB[p][i]  (elementwise product)
     */
    void dot_times(const PartialLikelihoodTensor& vecsA, const PartialLikelihoodTensor& vecsB)
    {
        assert(size == vecsA.size && size == vecsB.size);
        assert(isVectors && vecsA.isVectors && vecsB.isVectors);

              Scalar* Cdata = data.data();
        const Scalar* Adata = vecsA.data.data();
        const Scalar* Bdata = vecsB.data.data();

        for (std::size_t p = 0; p < size; p++)
        {
                  Scalar* Cp = Cdata + 4*p;
            const Scalar* Ap = Adata + 4*p;
            const Scalar* Bp = Bdata + 4*p;
            Cp[0] = Ap[0] * Bp[0];
            Cp[1] = Ap[1] * Bp[1];
            Cp[2] = Ap[2] * Bp[2];
            Cp[3] = Ap[3] * Bp[3];
        }
    }

    /**
     * this[p][i][j] = vecsA[p][i] * matsB[p][i][j]
     * If matsB.size == 1, the same matrix is used for all p.
     */
    void scale_rows(const PartialLikelihoodTensor& vecsA, const PartialLikelihoodTensor& matsB)
    {
        assert(!isVectors && size == vecsA.size && vecsA.isVectors && !matsB.isVectors);

              Scalar* Cdata = data.data();
        const Scalar* Adata = vecsA.data.data();
        const Scalar* Bdata = matsB.data.data();

        if (matsB.size == 1)
        {
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 4*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Ap[0] * Bdata[0];   Cp[1]  = Ap[0] * Bdata[1];   Cp[2]  = Ap[0] * Bdata[2];   Cp[3]  = Ap[0] * Bdata[3];
                Cp[4]  = Ap[1] * Bdata[4];   Cp[5]  = Ap[1] * Bdata[5];   Cp[6]  = Ap[1] * Bdata[6];   Cp[7]  = Ap[1] * Bdata[7];
                Cp[8]  = Ap[2] * Bdata[8];   Cp[9]  = Ap[2] * Bdata[9];   Cp[10] = Ap[2] * Bdata[10];  Cp[11] = Ap[2] * Bdata[11];
                Cp[12] = Ap[3] * Bdata[12];  Cp[13] = Ap[3] * Bdata[13];  Cp[14] = Ap[3] * Bdata[14];  Cp[15] = Ap[3] * Bdata[15];
            }
        }
        else
        {
            assert(matsB.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 4*p;
                const Scalar* Bp = Bdata + 16*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Ap[0] * Bp[0];   Cp[1]  = Ap[0] * Bp[1];   Cp[2]  = Ap[0] * Bp[2];   Cp[3]  = Ap[0] * Bp[3];
                Cp[4]  = Ap[1] * Bp[4];   Cp[5]  = Ap[1] * Bp[5];   Cp[6]  = Ap[1] * Bp[6];   Cp[7]  = Ap[1] * Bp[7];
                Cp[8]  = Ap[2] * Bp[8];   Cp[9]  = Ap[2] * Bp[9];   Cp[10] = Ap[2] * Bp[10];  Cp[11] = Ap[2] * Bp[11];
                Cp[12] = Ap[3] * Bp[12];  Cp[13] = Ap[3] * Bp[13];  Cp[14] = Ap[3] * Bp[14];  Cp[15] = Ap[3] * Bp[15];
            }
        }
    }

    /**
     * this[p][i][j] = matsA[p][i][j] * vecsB[p][j]
     * If matsA.size == 1, the same matrix is used for all p.
     */
    void scale_columns(const PartialLikelihoodTensor& matsA, const PartialLikelihoodTensor& vecsB)
    {
        assert(!isVectors && size == vecsB.size && !matsA.isVectors && vecsB.isVectors);

              Scalar* Cdata = data.data();
        const Scalar* Adata = matsA.data.data();
        const Scalar* Bdata = vecsB.data.data();

        if (matsA.size == 1)
        {
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Bp = Bdata + 4*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Adata[0]  * Bp[0];  Cp[1]  = Adata[1]  * Bp[1];  Cp[2]  = Adata[2]  * Bp[2];  Cp[3]  = Adata[3]  * Bp[3];
                Cp[4]  = Adata[4]  * Bp[0];  Cp[5]  = Adata[5]  * Bp[1];  Cp[6]  = Adata[6]  * Bp[2];  Cp[7]  = Adata[7]  * Bp[3];
                Cp[8]  = Adata[8]  * Bp[0];  Cp[9]  = Adata[9]  * Bp[1];  Cp[10] = Adata[10] * Bp[2];  Cp[11] = Adata[11] * Bp[3];
                Cp[12] = Adata[12] * Bp[0];  Cp[13] = Adata[13] * Bp[1];  Cp[14] = Adata[14] * Bp[2];  Cp[15] = Adata[15] * Bp[3];
            }
        }
        else
        {
            assert(matsA.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 16*p;
                const Scalar* Bp = Bdata + 4*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Ap[0]  * Bp[0];  Cp[1]  = Ap[1]  * Bp[1];  Cp[2]  = Ap[2]  * Bp[2];  Cp[3]  = Ap[3]  * Bp[3];
                Cp[4]  = Ap[4]  * Bp[0];  Cp[5]  = Ap[5]  * Bp[1];  Cp[6]  = Ap[6]  * Bp[2];  Cp[7]  = Ap[7]  * Bp[3];
                Cp[8]  = Ap[8]  * Bp[0];  Cp[9]  = Ap[9]  * Bp[1];  Cp[10] = Ap[10] * Bp[2];  Cp[11] = Ap[11] * Bp[3];
                Cp[12] = Ap[12] * Bp[0];  Cp[13] = Ap[13] * Bp[1];  Cp[14] = Ap[14] * Bp[2];  Cp[15] = Ap[15] * Bp[3];
            }
        }
    }

    /**
     * this[p] = matsA[p] * vecsB[p]  (matrix-vector product per site)
     * If matsA.size == 1, the same matrix is applied to every site.
     */
    void matrix_vector_product(const PartialLikelihoodTensor& matsA, const PartialLikelihoodTensor& vecsB)
    {
        assert(size == vecsB.size && isVectors && !matsA.isVectors && vecsB.isVectors);

              Scalar* Cdata = data.data();
        const Scalar* Adata = matsA.data.data();
        const Scalar* Bdata = vecsB.data.data();

        if (matsA.size == 1)
        {
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Bp = Bdata + 4*p;
                      Scalar* Cp = Cdata + 4*p;
                Cp[0] = Adata[0]*Bp[0] + Adata[1]*Bp[1] + Adata[2]*Bp[2]  + Adata[3]*Bp[3];
                Cp[1] = Adata[4]*Bp[0] + Adata[5]*Bp[1] + Adata[6]*Bp[2]  + Adata[7]*Bp[3];
                Cp[2] = Adata[8]*Bp[0] + Adata[9]*Bp[1] + Adata[10]*Bp[2] + Adata[11]*Bp[3];
                Cp[3] = Adata[12]*Bp[0]+ Adata[13]*Bp[1]+ Adata[14]*Bp[2] + Adata[15]*Bp[3];
            }
        }
        else
        {
            assert(matsA.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 16*p;
                const Scalar* Bp = Bdata + 4*p;
                      Scalar* Cp = Cdata + 4*p;
                Cp[0] = Ap[0]*Bp[0] + Ap[1]*Bp[1] + Ap[2]*Bp[2]  + Ap[3]*Bp[3];
                Cp[1] = Ap[4]*Bp[0] + Ap[5]*Bp[1] + Ap[6]*Bp[2]  + Ap[7]*Bp[3];
                Cp[2] = Ap[8]*Bp[0] + Ap[9]*Bp[1] + Ap[10]*Bp[2] + Ap[11]*Bp[3];
                Cp[3] = Ap[12]*Bp[0]+ Ap[13]*Bp[1]+ Ap[14]*Bp[2] + Ap[15]*Bp[3];
            }
        }
    }

    /**
     * this[p] = matsA[p] * matsB[p]  (matrix-matrix product per site)
     * Either input may have size 1, in which case that matrix is broadcast.
     */
    void matrix_matrix_product(const PartialLikelihoodTensor& matsA, const PartialLikelihoodTensor& matsB)
    {
              Scalar* Cdata = data.data();
        const Scalar* Adata = matsA.data.data();
        const Scalar* Bdata = matsB.data.data();

        if (matsA.size == 1 && matsB.size > 1)
        {
            assert(matsB.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Bp = Bdata + 16*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Adata[0]*Bp[0]  + Adata[1]*Bp[4]  + Adata[2]*Bp[8]   + Adata[3]*Bp[12];
                Cp[1]  = Adata[0]*Bp[1]  + Adata[1]*Bp[5]  + Adata[2]*Bp[9]   + Adata[3]*Bp[13];
                Cp[2]  = Adata[0]*Bp[2]  + Adata[1]*Bp[6]  + Adata[2]*Bp[10]  + Adata[3]*Bp[14];
                Cp[3]  = Adata[0]*Bp[3]  + Adata[1]*Bp[7]  + Adata[2]*Bp[11]  + Adata[3]*Bp[15];
                Cp[4]  = Adata[4]*Bp[0]  + Adata[5]*Bp[4]  + Adata[6]*Bp[8]   + Adata[7]*Bp[12];
                Cp[5]  = Adata[4]*Bp[1]  + Adata[5]*Bp[5]  + Adata[6]*Bp[9]   + Adata[7]*Bp[13];
                Cp[6]  = Adata[4]*Bp[2]  + Adata[5]*Bp[6]  + Adata[6]*Bp[10]  + Adata[7]*Bp[14];
                Cp[7]  = Adata[4]*Bp[3]  + Adata[5]*Bp[7]  + Adata[6]*Bp[11]  + Adata[7]*Bp[15];
                Cp[8]  = Adata[8]*Bp[0]  + Adata[9]*Bp[4]  + Adata[10]*Bp[8]  + Adata[11]*Bp[12];
                Cp[9]  = Adata[8]*Bp[1]  + Adata[9]*Bp[5]  + Adata[10]*Bp[9]  + Adata[11]*Bp[13];
                Cp[10] = Adata[8]*Bp[2]  + Adata[9]*Bp[6]  + Adata[10]*Bp[10] + Adata[11]*Bp[14];
                Cp[11] = Adata[8]*Bp[3]  + Adata[9]*Bp[7]  + Adata[10]*Bp[11] + Adata[11]*Bp[15];
                Cp[12] = Adata[12]*Bp[0] + Adata[13]*Bp[4] + Adata[14]*Bp[8]  + Adata[15]*Bp[12];
                Cp[13] = Adata[12]*Bp[1] + Adata[13]*Bp[5] + Adata[14]*Bp[9]  + Adata[15]*Bp[13];
                Cp[14] = Adata[12]*Bp[2] + Adata[13]*Bp[6] + Adata[14]*Bp[10] + Adata[15]*Bp[14];
                Cp[15] = Adata[12]*Bp[3] + Adata[13]*Bp[7] + Adata[14]*Bp[11] + Adata[15]*Bp[15];
            }
        }
        else if (matsA.size > 1 && matsB.size == 1)
        {
            assert(matsA.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 16*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Ap[0]*Bdata[0]  + Ap[1]*Bdata[4]  + Ap[2]*Bdata[8]   + Ap[3]*Bdata[12];
                Cp[1]  = Ap[0]*Bdata[1]  + Ap[1]*Bdata[5]  + Ap[2]*Bdata[9]   + Ap[3]*Bdata[13];
                Cp[2]  = Ap[0]*Bdata[2]  + Ap[1]*Bdata[6]  + Ap[2]*Bdata[10]  + Ap[3]*Bdata[14];
                Cp[3]  = Ap[0]*Bdata[3]  + Ap[1]*Bdata[7]  + Ap[2]*Bdata[11]  + Ap[3]*Bdata[15];
                Cp[4]  = Ap[4]*Bdata[0]  + Ap[5]*Bdata[4]  + Ap[6]*Bdata[8]   + Ap[7]*Bdata[12];
                Cp[5]  = Ap[4]*Bdata[1]  + Ap[5]*Bdata[5]  + Ap[6]*Bdata[9]   + Ap[7]*Bdata[13];
                Cp[6]  = Ap[4]*Bdata[2]  + Ap[5]*Bdata[6]  + Ap[6]*Bdata[10]  + Ap[7]*Bdata[14];
                Cp[7]  = Ap[4]*Bdata[3]  + Ap[5]*Bdata[7]  + Ap[6]*Bdata[11]  + Ap[7]*Bdata[15];
                Cp[8]  = Ap[8]*Bdata[0]  + Ap[9]*Bdata[4]  + Ap[10]*Bdata[8]  + Ap[11]*Bdata[12];
                Cp[9]  = Ap[8]*Bdata[1]  + Ap[9]*Bdata[5]  + Ap[10]*Bdata[9]  + Ap[11]*Bdata[13];
                Cp[10] = Ap[8]*Bdata[2]  + Ap[9]*Bdata[6]  + Ap[10]*Bdata[10] + Ap[11]*Bdata[14];
                Cp[11] = Ap[8]*Bdata[3]  + Ap[9]*Bdata[7]  + Ap[10]*Bdata[11] + Ap[11]*Bdata[15];
                Cp[12] = Ap[12]*Bdata[0] + Ap[13]*Bdata[4] + Ap[14]*Bdata[8]  + Ap[15]*Bdata[12];
                Cp[13] = Ap[12]*Bdata[1] + Ap[13]*Bdata[5] + Ap[14]*Bdata[9]  + Ap[15]*Bdata[13];
                Cp[14] = Ap[12]*Bdata[2] + Ap[13]*Bdata[6] + Ap[14]*Bdata[10] + Ap[15]*Bdata[14];
                Cp[15] = Ap[12]*Bdata[3] + Ap[13]*Bdata[7] + Ap[14]*Bdata[11] + Ap[15]*Bdata[15];
            }
        }
        else
        {
            assert(matsA.size == size && matsB.size == size);
            for (std::size_t p = 0; p < size; p++)
            {
                const Scalar* Ap = Adata + 16*p;
                const Scalar* Bp = Bdata + 16*p;
                      Scalar* Cp = Cdata + 16*p;
                Cp[0]  = Ap[0]*Bp[0]  + Ap[1]*Bp[4]  + Ap[2]*Bp[8]   + Ap[3]*Bp[12];
                Cp[1]  = Ap[0]*Bp[1]  + Ap[1]*Bp[5]  + Ap[2]*Bp[9]   + Ap[3]*Bp[13];
                Cp[2]  = Ap[0]*Bp[2]  + Ap[1]*Bp[6]  + Ap[2]*Bp[10]  + Ap[3]*Bp[14];
                Cp[3]  = Ap[0]*Bp[3]  + Ap[1]*Bp[7]  + Ap[2]*Bp[11]  + Ap[3]*Bp[15];
                Cp[4]  = Ap[4]*Bp[0]  + Ap[5]*Bp[4]  + Ap[6]*Bp[8]   + Ap[7]*Bp[12];
                Cp[5]  = Ap[4]*Bp[1]  + Ap[5]*Bp[5]  + Ap[6]*Bp[9]   + Ap[7]*Bp[13];
                Cp[6]  = Ap[4]*Bp[2]  + Ap[5]*Bp[6]  + Ap[6]*Bp[10]  + Ap[7]*Bp[14];
                Cp[7]  = Ap[4]*Bp[3]  + Ap[5]*Bp[7]  + Ap[6]*Bp[11]  + Ap[7]*Bp[15];
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
    }

    /**
     * Returns a vector of length size where entry p = pi^T * column[p].
     */
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> dot_product(const Eigen::Matrix<Scalar,4,1>& pi)
    {
        assert(isVectors);
        Eigen::Matrix<Scalar, Eigen::Dynamic, 1> lvec(size, 1);
        const Scalar* Ddata = data.data();
        const Scalar pi0 = pi(0), pi1 = pi(1), pi2 = pi(2), pi3 = pi(3);
        for (std::size_t p = 0; p < size; p++)
        {
            const Scalar* Dp = Ddata + 4*p;
            lvec(p) = pi0*Dp[0] + pi1*Dp[1] + pi2*Dp[2] + pi3*Dp[3];
        }
        return lvec;
    }

    // -----------------------------------------------------------------------
    // Utility — raw pointer only; inner loops are trivially handled by the
    // compiler and are not hot-path.
    // -----------------------------------------------------------------------

    /**
     * Returns the maximum entry in position p (used for overflow handling).
     */
    Scalar max_coefficient(std::size_t p) const
    {
        assert(p < size);
        Scalar max_val = 0.0;
        if (isVectors)
        {
            const Scalar* Dp = data.data() + 4*p;
            for (std::size_t i = 0; i < 4; i++)
                max_val = std::max(max_val, Dp[i]);
        }
        else
        {
            const Scalar* Dp = data.data() + 16*p;
            for (std::size_t k = 0; k < 16; k++)
                max_val = std::max(max_val, Dp[k]);
        }
        return max_val;
    }

    /**
     * Multiplies every entry in position p by multiplier (used for overflow handling).
     */
    void rescale(std::size_t p, Scalar multiplier)
    {
        assert(p < size);
        if (isVectors)
        {
            Scalar* Dp = data.data() + 4*p;
            Dp[0] *= multiplier;  Dp[1] *= multiplier;  Dp[2] *= multiplier;  Dp[3] *= multiplier;
        }
        else
        {
            Scalar* Dp = data.data() + 16*p;
            for (std::size_t k = 0; k < 16; k++)
                Dp[k] *= multiplier;
        }
    }
};

#endif // USE_REFERENCE_TENSOR_CALCULATIONS
#endif // LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_H
