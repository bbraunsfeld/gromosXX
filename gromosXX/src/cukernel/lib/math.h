/**
 * @file math.h
 * some basic math operations
 */
#ifndef CUKERNEL_MATH_H
#define CUKERNEL_MATH_H

// float3 operations
#define HOSTDEVICE __device__ inline
#include "float3.h"
#include "double3.h"
#undef HOSTDEVICE

// additional types
#include "types.h"

// For the precision
#include "../macros.h"

/**
 * a very small number
 */
#define EPSILON 0.000001f
#define EPSILOND PREC(0.000001)

#endif
