#pragma once

// MIT License

// Copyright (c) 2020 xaedes

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#define VEC_USE_NATIVE 1

#ifdef VEC_USE_OPENCV
#include <opencv2/opencv2.hpp>
#endif
#ifdef VEC_USE_NATIVE
#include <cmath>
#endif

#ifdef VEC_USE_OPENCV
    using cv::Vec;
    using cv::norm;
#endif
#ifdef VEC_USE_NATIVE
#define ARITHMETIC_OPERATOR(OP) \
    Vec<TValueType,TNumDimensions>& operator OP ## =(const Vec<TValueType,TNumDimensions>& rhs)                                      \
    {                                                                                                                                \
        for(unsigned int i=0; i<TNumDimensions; ++i)                                                                                 \
            *this[i] OP ## = rhs[i];                                                                                                 \
        return *this;                                                                                                                \
    }                                                                                                                                \
                                                                                                                                     \
    friend Vec<TValueType,TNumDimensions> operator OP(Vec<TValueType,TNumDimensions> lhs, const Vec<TValueType,TNumDimensions>& rhs) \
    {                                                                                                                                \
        for(unsigned int i=0; i<TNumDimensions; ++i)                                                                                 \
            lhs[i] OP ## = rhs[i];                                                                                                   \
        return lhs;                                                                                                                  \
    }                                                                                                                                \
                                                                                                                                     \
    Vec<TValueType,TNumDimensions>& operator OP ## =(const TValueType& rhs)                                                          \
    {                                                                                                                                \
        for(unsigned int i=0; i<TNumDimensions; ++i)                                                                                 \
            *this[i] OP ## = rhs;                                                                                                    \
        return *this;                                                                                                                \
    }                                                                                                                                \
                                                                                                                                     \
    friend Vec<TValueType,TNumDimensions> operator OP(Vec<TValueType,TNumDimensions> lhs, const TValueType& rhs)                     \
    {                                                                                                                                \
        for(unsigned int i=0; i<TNumDimensions; ++i)                                                                                 \
            lhs[i] OP ## = rhs;                                                                                                      \
        return lhs;                                                                                                                  \
    }                                                                                                                                \
    \

    template<typename TValueType, unsigned int TNumDimensions>
    struct Vec
    {
        TValueType operator[](int i) const { return data[i]; }
        TValueType& operator[](int i) { return data[i]; }
        TValueType data[TNumDimensions];
        Vec() : Vec(0) {}
        Vec(TValueType value)
        {
            for (unsigned int i = 0; i < TNumDimensions; ++i)
                data[i] = value;
        }
        template<typename ... TArgs> Vec(TValueType first, TArgs... args) : data{ first, (TValueType)args... } {}
        ARITHMETIC_OPERATOR(+)
        ARITHMETIC_OPERATOR(-)
        ARITHMETIC_OPERATOR(*)
        ARITHMETIC_OPERATOR(/ )
    };

    template<typename TValueType, unsigned int TNumDimensions>
    TValueType norm(Vec<TValueType, TNumDimensions> vec)
    {
        TValueType sum = 0;
        for (unsigned int i = 0; i < TNumDimensions; ++i)
            sum += vec[i] * vec[i];
        return sqrt(sum);
    }
#endif

