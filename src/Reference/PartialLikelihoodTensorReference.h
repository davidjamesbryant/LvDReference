//
// Created by David Bryant on 14/05/2026.
// Optimized pass: raw pointers, site-local offsets, and unrolled 4-entry kernels.
//

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

    PartialLikelihoodTensorReference(bool isVectors_, std::size_t size): isVectors(isVectors_), size(size)
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
    void dot_times(const PartialLikelihoodTensorReference& vecsA, const PartialLikelihoodTensorReference& vecsB)
    {
        assert(size == vecsA.size && size == vecsB.size);
        assert(isVectors && vecsA.isVectors && vecsB.isVectors);

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
    void scale_rows(const PartialLikelihoodTensorReference& vecsA, const PartialLikelihoodTensorReference& matsB)
    {
        assert(!isVectors && size == vecsA.size && vecsA.isVectors && !matsB.isVectors);
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