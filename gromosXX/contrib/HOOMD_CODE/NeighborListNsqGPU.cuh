/*
Highly Optimized Object-Oriented Molecular Dynamics (HOOMD) Open
Source Software License
Copyright (c) 2008 Ames Laboratory Iowa State University
All rights reserved.

Redistribution and use of HOOMD, in source and binary forms, with or
without modification, are permitted, provided that the following
conditions are met:

* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names HOOMD's
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

Disclaimer

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND
CONTRIBUTORS ``AS IS''  AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 

IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS  BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.
*/

// $Id: NeighborListNsqGPU.cuh 1826 2009-04-27 22:18:31Z joaander $
// $URL: https://codeblue.umich.edu/hoomd-blue/svn/tags/hoomd-0.8.2/src/cuda/NeighborListNsqGPU.cuh $
// Maintainer: joaander

#ifndef _NEIGHBORLISTNSQGPU_CUH_
#define _NEIGHBORLISTNSQGPU_CUH_

#include <stdio.h>
#include <cuda_runtime.h>

#include "NeighborList.cuh"

/*! \file NeighborListNsqGPU.cuh
	\brief Declares data structures and methods used by NeighborListNsqGPU
*/

//! Kernel driver to run GPU code called by NeighborListNsqGPU
cudaError_t gpu_compute_nlist_nsq(const gpu_nlist_array &nlist, const gpu_pdata_arrays &pdata, const gpu_boxsize &box, float r_maxsq);

#endif
