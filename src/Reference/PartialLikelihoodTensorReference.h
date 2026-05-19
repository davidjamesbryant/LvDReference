//
// Created by David Bryant on 14/05/2026.
//
// This is the REFERENCE implementation of PartialLikelihoodTensor.
// It is intentionally written for clarity: clean loops, no manual unrolling,
// no platform-specific attributes.
//
// The optimised counterpart lives in src/DecompositionTree/PartialLikelihoodTensor.h.
// When this file changes, apply the following transformations to regenerate it.
// Each optimisation is also flagged with an "OPT-N" marker in the method where it applies.
//
// ═══════════════════════════════════════════════════════════════════════════════
// OPTIMISATION GUIDE — how to derive PartialLikelihoodTensor.h from this file
// ═══════════════════════════════════════════════════════════════════════════════
//
// OPT-1  Raw pointer access
//   Replace every indexed access   data[expr]      with a raw pointer obtained
//   once per function:
//       const Scalar* Adata = matsA.data.data();
//       Scalar*       Cdata = data.data();
//   then use Adata[expr], Cdata[expr] throughout.
//   Rationale: std::vector::operator[] is bounds-checked in Debug builds and
//   carries a pointer-dereference indirection even in Release.  Raw pointers
//   eliminate the indirection and let the compiler keep the base address in a
//   register for the duration of the loop.
//
// OPT-2  Site-local offset pointers
//   Inside the outer p-loop, compute pointers to the start of each site's
//   block once, before the inner i/j/k loops:
//       const Scalar* Ap = Adata + 16*p;   // or 4*p for vectors
//       const Scalar* Bp = Bdata +  4*p;
//             Scalar* Cp = Cdata + 16*p;
//   then write  Ap[4*i+j]  instead of  Adata[16*p + 4*i + j].
//   Rationale: the 16*p or 4*p offset would otherwise be recomputed for every
//   access inside the inner loop.  The compiler often hoists it anyway, but
//   making it explicit guarantees the saving and aids readability.
//
// OPT-3  Inner loop unrolling (size-4 loops only)
//   Replace every inner loop of trip-count 4 with 4 explicit statements.
//   Example for the i-loop in matrix_vector_product (broadcast case):
//       Cp[0] = Adata[0]*Bp[0] + Adata[1]*Bp[1] + Adata[2]*Bp[2] + Adata[3]*Bp[3];
//       Cp[1] = Adata[4]*Bp[0] + Adata[5]*Bp[1] + Adata[6]*Bp[2] + Adata[7]*Bp[3];
//       Cp[2] = Adata[8]*Bp[0] + ...;
//       Cp[3] = Adata[12]*Bp[0] + ...;
//   The outermost p-loop (trip-count = nSites, unknown at compile time) is
//   left as a loop — the compiler vectorises it automatically once the inner
//   loops are gone.  Do NOT unroll the p-loop.
//   Rationale: the compiler cannot unroll a loop whose bounds are runtime
//   values.  Explicit unrolling of the known-size-4 loops removes the loop
//   overhead, exposes instruction-level parallelism, and allows the register
//   allocator to keep all 4 result accumulators live simultaneously.
//   Leave set_column, fill_column, set(mat), max_coefficient, rescale as
//   plain loops — they are not on the hot path.
//
// OPT-4  Force inlining on hot kernels
//   Add  __attribute__((always_inline))  before the return type of every
//   computational kernel method (dot_times, scale_rows, scale_columns,
//   matrix_vector_product, matrix_matrix_product, dot_product).
//   Do NOT add it to max_coefficient, rescale, or the setters.
//   Rationale: sampling profiling showed matrix_vector_product appearing as a
//   separate call frame consuming ~29% of pruning MCMC likelihood time, even
//   at -O3.  The function bodies exceed the compiler's default inlining
//   threshold.  Forcing inlining lets the compiler fuse the arithmetic with
//   the surrounding ancestor-walk loop in decompositionLikelihood.cpp and
//   apply cross-call optimisations (register reuse, better scheduling).
//
// OPT-5  Restrict-qualified pointer locals
//   Add  __restrict__  to every raw pointer local declared in the hot kernels:
//       Scalar* __restrict__       Cdata = data.data();
//       const Scalar* __restrict__ Adata = matsA.data.data();
//       const Scalar* __restrict__ Bdata = vecsB.data.data();
//   Rationale: the three data arrays belong to distinct PartialLikelihoodTensor
//   objects and never alias.  Without __restrict__ the compiler must conservatively
//   assume any write through Cdata could affect what is read through Adata or Bdata,
//   blocking auto-vectorisation of the outer p-loop.  With __restrict__ the compiler
//   can generate SIMD (SSE/AVX) code that processes multiple sites per cycle.
//
// OPT-6  Constructor parameter shadowing fix
//   The constructor  PartialLikelihoodTensorReference(bool isVectors_, std::size_t size)
//   has a parameter named  size  that shadows the member of the same name.
//   In the optimised version rename it to  size_  for clarity.
//
// ═══════════════════════════════════════════════════════════════════════════════

#ifndef LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_REFERENCE_H
#define LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_REFERENCE_H


// tiny4_tensor.hpp
#pragma once

#include <vector>
#include <cstddef>
#include <cassert>
#include <Eigen>

#include "phylib.h"


class PartialLikelihoodTensorReference
{
private:
    std::vector<Scalar> data;
    bool isVectors;  //True if this is 4xnsites, false if this is 4x4xnSites
    std::size_t size; //1 if this is a single vector or matrix, more if this is an array of vectors or matrices
public:
    PartialLikelihoodTensorReference() : isVectors(false), size(0) {}  //Default constructor - no allocation

    PartialLikelihoodTensorReference(const PartialLikelihoodTensorReference& x) : isVectors(x.isVectors), size(x.size)
    {
        resize(isVectors,size);
        copy(x.data.begin(), x.data.end(), data.begin());
    }

    PartialLikelihoodTensorReference(bool isVectors_, std::size_t size): isVectors(isVectors_), size(size) // OPT-6: rename parameter to size_
    {
        if (isVectors)
            data.resize(4*size);
        else
            data.resize(16*size);
    }

    /**
     * Change the type of storage and reallocate as necessary
     * @param isVectors_
     * @param size_
     */
    void resize(bool isVectors_, std::size_t size_)
    {
        isVectors = isVectors_;
        size = size_;
        if (isVectors)
            data.resize(4*size);
        else
            data.resize(16*size);
    }

    /**
     * Set entry i in position p of this 4xs tensor
     * @param s
     * @param i
     * @param val
     */
    void set(std::size_t p, std::size_t  i, Scalar val)
    {
        assert(isVectors && p < size && i < 4);
        data[4*p + i] = val;
    }

    /**
     * Copy col into position p of this  4xs tensor
     * @param p position number
     * @param col column vector
     */
    void set_column(std::size_t p, const Eigen::Matrix<Scalar,4,1>& col)
    {
        assert(isVectors && p < size && col.size() == 4);
        for (std::size_t i=0;i<4;i++)
            data[4*p + i] = col(i);
    }

    /**
     * Fill column p with repeated val values
     * @param p
     * @param val
     */
    void fill_column(std::size_t p, Scalar val)
    {
        assert(isVectors && p<size);
        for (std::size_t i=0;i<4;i++)
        {
            data[p*4 + i] = val;
        }
    }


    /**
     * Copy mat into this 4x4 tensor
     * @param mat
     */
    void set(const Eigen::Matrix<Scalar,4,4>& mat)
    {
        assert(!isVectors&&size==1);
        for (std::size_t i=0;i<4;i++)
            for (std::size_t j=0;j<4;j++)
                data[4*i + j] = mat(i,j);
    }

    /**
     * Make this tensor equal to A.*B
     * @param vecsA
     * @param vecsB
     */
    // OPT-4: add __attribute__((always_inline))
    void dot_times(const PartialLikelihoodTensorReference& vecsA, const PartialLikelihoodTensorReference& vecsB)
    {
        assert(size == vecsA.size && size == vecsB.size);
        assert(isVectors && vecsA.isVectors && vecsB.isVectors);

        // OPT-1/5: raw __restrict__ pointers; OPT-2: hoist Ap/Bp/Cp inside p-loop; OPT-3: unroll i-loop
        for (std::size_t p=0;p<size;p++)
            for (std::size_t i=0;i<4;i++)
                data[4*p + i] = vecsA.data[4*p + i] * vecsB.data[4*p + i];
    }

    /**
     * Scale rows of matrix in position p of matsB with vector from position p of vecsA. Fill this tensor
     * with the result. If matsB has size 1 then pretend like matsB is copies of the same matrix.
     * then
     * @param vecsA
     * @param matsB
     */
    // OPT-4: add __attribute__((always_inline))
    void scale_rows(const PartialLikelihoodTensorReference& vecsA, const PartialLikelihoodTensorReference& matsB)
    {
        assert(!isVectors && size == vecsA.size && vecsA.isVectors && !matsB.isVectors);
        // OPT-1/5: raw __restrict__ pointers; OPT-2: hoist Ap/Bp/Cp inside p-loop; OPT-3: unroll i and j loops
        if (matsB.size==1)
        {
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*p + 4*i + j] = vecsA.data[4*p + i] * matsB.data[4*i + j];
        } else
        {
            assert(matsB.size==size);
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*p + 4*i + j] = vecsA.data[4*p + i] * matsB.data[16*p + 4*i + j];
        }
    }

    /**
     * Scale columns of matrix in position p of matsA with vector from position p of vecsB. Fill this tensor
     * with the result. If matsA has size 1 then pretend like matsA is copies of the same matrix.
     * then
     * @param matsA
     * @param vecsB
     */
    // OPT-4: add __attribute__((always_inline))
    // OPT-1/5: raw __restrict__ pointers; OPT-2: hoist Ap/Bp/Cp inside p-loop; OPT-3: unroll i and j loops
    void scale_columns(const PartialLikelihoodTensorReference& matsA, const PartialLikelihoodTensorReference& vecsB)
    {
        assert(!isVectors && size == vecsB.size && !matsA.isVectors && vecsB.isVectors);
        if (matsA.size==1)
        {
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*p + 4*i + j] = matsA.data[4*i + j] * vecsB.data[4*p + j];
        } else
        {
            assert(matsA.size == size);
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                        data[16*p + 4*i + j] = matsA.data[16*p + 4*i + j] * vecsB.data[4*p + j];
        }
    }

    /**
     * Multiply vector in each position of vecsB by the matrix in the same position of matsA. If matsA has size 1
     * then we multiply each vector of vecsB by the same matrix
     * @param matsA
     * @param vecsB
     */
    // OPT-4: add __attribute__((always_inline))  ← confirmed hot by sample profiler (~29% of pruning MCMC time)
    // OPT-1/5: raw __restrict__ pointers; OPT-2: hoist Ap/Bp/Cp inside p-loop; OPT-3: unroll i and j loops
    void matrix_vector_product(const PartialLikelihoodTensorReference& matsA, const PartialLikelihoodTensorReference& vecsB)
    {
        assert(size==vecsB.size && isVectors && !matsA.isVectors && vecsB.isVectors);
        if (matsA.size==1)
        {
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                {
                    Scalar sum = 0.0;
                    for (std::size_t j=0;j<4;j++)
                        sum+=matsA.data[4*i + j] * vecsB.data[4*p + j];
                    data[4*p + i] = sum;
                }
        } else
        {
            assert(matsA.size == size);
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                {
                    Scalar sum = 0.0;
                    for (std::size_t j=0;j<4;j++)
                        sum+=matsA.data[16*p + 4*i + j] * vecsB.data[4*p + j];
                    data[4*p + i] = sum;
                }
        }
    }

    /**
     * Multiply the matrix in position p of matsA with the matrix in position p of matsB. If either of matsA or matsB
     * have size one then we pretend its an array full of identical copies of that matrix.
     * @param matsA
     * @param matsB
     */
    // OPT-4: add __attribute__((always_inline))
    // OPT-1/5: raw __restrict__ pointers; OPT-2: hoist Ap/Bp/Cp inside p-loop; OPT-3: unroll i, j, k loops
    void matrix_matrix_product(const PartialLikelihoodTensorReference& matsA, const PartialLikelihoodTensorReference& matsB)
    {

        if (matsA.size == 1 && matsB.size > 1)
        {
            assert(matsB.size == size);
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                    {
                        Scalar sum = 0.0;
                        for (std::size_t k=0;k<4;k++)
                            sum += matsA.data[4*i + k] * matsB.data[16*p+4*k + j];
                        data[16*p + 4*i + j] = sum;
                    }

        } else if (matsA.size > 1 && matsB.size == 1)
        {
            assert(matsA.size == size);
            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                    {
                        Scalar sum = 0.0;
                        for (std::size_t k=0;k<4;k++)
                            sum += matsA.data[16*p + 4*i + k] * matsB.data[4*k + j];
                        data[16*p + 4*i + j] = sum;
                    }
        } else //Arrays of matrices same length - multiple matrix in position p of matsA with matrix in position p of matsB
        {
            assert(matsA.size == size && matsB.size == size);

            for (std::size_t p=0;p<size;p++)
                for (std::size_t i=0;i<4;i++)
                    for (std::size_t j=0;j<4;j++)
                    {
                        Scalar sum = 0.0;
                        for (std::size_t k=0;k<4;k++)
                            sum += matsA.data[16*p + 4*i + k] * matsB.data[16*p + 4*k + j];
                        data[16*p + 4*i + j] = sum;
                    }
        }
    }

    /**
     * Compute  pi^Tv for the vector in every position of this array, returning the dot products
     * as an Eigen vector with 'size' entries.
     * @param pi
     * @return
     */
    // OPT-4: add __attribute__((always_inline))
    // OPT-1/5: raw __restrict__ pointer; OPT-3: unroll i-loop; hoist pi components to scalars
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> dot_product(const Eigen::Matrix<Scalar,4,1>& pi)
    {
        assert(isVectors);
        Eigen::Matrix<Scalar, Eigen::Dynamic , 1> lvec(size,1);
        for (std::size_t p=0;p<size;p++)
        {
            Scalar sum = 0.0;
            for (std::size_t i=0;i<4;i++)
                sum += pi(i) * data[4*p + i];
            lvec(p) = sum;
        }
        return lvec;
    }

    /**
     * Computes the maximum value for the vector or matrix in position p. Used for handling overflow.
     * @param p
     * @return
     */
    Scalar max_coefficient(std::size_t p) const
    {
        assert(p < size);
        Scalar max_val = 0.0; //Tensors non-negative
        if (isVectors)
        {
            for (std::size_t i=0;i<4;i++)
                max_val = std::max(max_val,data[4*p + i]);
        } else
        {
            for (std::size_t i=0;i<4;i++)
                for (std::size_t j=0;j<4;j++)
                    max_val = std::max(max_val,data[16*p + 4*i + j]);

        }
        return max_val;
    }

    /**
     * Multiples every entry of the vector or matrix in position p by the multiplier.
     * Used for handling overflow.
     *
     * @param p
     * @param multiplier
     */
    void rescale(std::size_t p, Scalar multiplier)
    {
        assert(p < size);
        if (isVectors)
        {
            for (std::size_t i=0;i<4;i++)
                data[4*p + i] *= multiplier;
        } else
        {
            for (std::size_t i=0;i<4;i++)
                for (std::size_t j=0;j<4;j++)
                    data[16*p + 4*i + j] *= multiplier;
        }
    }
};

#endif // LVDREFERENCE_PARTIAL_LIKELIHOOD_TENSOR_REFERENCE_H