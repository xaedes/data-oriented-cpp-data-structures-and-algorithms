#pragma once

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

