#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <time.h>

#include "Aster/impl/kernels/sorting.cl.h"
#include "Aster/building-api/logging.h"

namespace Aster{
namespace GPU{


extern cl_platform_id platform;
extern cl_device_id device;
extern cl_context context;
extern cl_command_queue queue;

constexpr size_t WORKGROUP_SIZE  = 256;
constexpr size_t NUM_BUCKETS     = 256;
constexpr size_t NUM_PASSES      = 4;
constexpr float MAX_MEM_FRACTION = 0.80f;

#define CHECK_CL(err, msg)                                      \
    if ((err) != CL_SUCCESS) {                                  \
        fprintf(stderr, "[ERRORE] %s → %d\n", msg, err);        \
        exit(EXIT_FAILURE);                                     \
    }

/**
* @brief: gets_the maximum amount of data we can put on the GPU (assuming everything is a uint32_t)
* @param device: device to scan
* @returns the amount of uint32_t to load
*/
static size_t compute_max_chunk_size(cl_device_id device) {
    cl_int err;
    cl_ulong global_mem = 0, max_alloc = 0;

    CHECK_CL(clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE,
                              sizeof(cl_ulong), &global_mem, NULL),
             "clGetDeviceInfo GLOBAL_MEM_SIZE");
    CHECK_CL(clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE,
                              sizeof(cl_ulong), &max_alloc, NULL),
             "clGetDeviceInfo MAX_MEM_ALLOC_SIZE");

    double usable_mem = (double)global_mem * MAX_MEM_FRACTION;

    size_t chunk = (size_t)((usable_mem) / (8.0 + 2.0 * (256.0 * 4.0 / WORKGROUP_SIZE)));
    chunk = (chunk / WORKGROUP_SIZE) * WORKGROUP_SIZE;
    if (chunk == 0) chunk = WORKGROUP_SIZE;

    auto mem_needed = [&](size_t M) {
        size_t num_groups = (M + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE;
        double hist_bytes = (double)num_groups * NUM_BUCKETS * sizeof(uint32_t);
        double total =  // input + output
                       (double)M * sizeof(uint32_t) * 2.0
                       // istogrammi + offsets
                       + 2.0 * hist_bytes;
        return total;
    };

    while (chunk >= WORKGROUP_SIZE) {
        double needed = mem_needed(chunk);
        if (needed <= usable_mem) {
            size_t bytes_buffer = chunk * sizeof(uint32_t);
            cl_ulong d_hist_size = (( (chunk + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE )
                                     * NUM_BUCKETS * sizeof(uint32_t));
            if (bytes_buffer <= max_alloc &&
                d_hist_size <= max_alloc) {
                break;
            }
        }
        chunk -= WORKGROUP_SIZE;
    }
    return chunk;
}

/**
* @brief merges two already sorted arrays
* @param left first array
* @param L first array size
* @param right second array
* @param R second array size
* @param merged ptr to the merged array
*/
static void merge_arrays(const uint32_t* left, size_t L,
                         const uint32_t* right, size_t R,
                         uint32_t* merged) {
    size_t i = 0, j = 0, k = 0;
    while (i < L && j < R) {
        if (left[i] <= right[j]) {
            merged[k++] = left[i++];
        } else {
            merged[k++] = right[j++];
        }
    }
    while (i < L) merged[k++] = left[i++];
    while (j < R) merged[k++] = right[j++];
}

inline void print_kernel_log(cl_int& err, cl_program& program){
    size_t log_size;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
    char* log = reinterpret_cast<char*>(malloc(log_size + 1));
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,log_size, log, NULL);
    log[log_size] = '\0';
    fprintf(stderr, "Error encountered while building sorting kernel:\n%s\n", log);
    free(log);
    exit(-1);
}

cl_program compile_sorting_kernels(){
    static bool loaded = false;
    static cl_program program;
    
    if (!loaded) loaded = true;
    else return program;

    const char* source_str = sorting_kernels_cl.c_str();
    size_t src_size = sorting_kernels_cl.size();
    cl_int err;

    program = clCreateProgramWithSource(context, 1, (const char**)&source_str, &src_size, &err);
    CHECK_CL(err, "clCreateProgramWithSource");

    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    return program;
}

void sort(uint32_t* input, uint32_t* output, size_t N) {
    cl_int err;
    static cl_program program = compile_sorting_kernels();
    static size_t max_chunk = compute_max_chunk_size(device);

    // compile the kernels once
    static cl_kernel histogram_kernel = clCreateKernel(program, "histogram_kernel", &err);
    CHECK_CL(err, "clCreateKernel histogram_kernel");
    static cl_kernel scatter_kernel = clCreateKernel(program, "scatter_kernel", &err);
    CHECK_CL(err, "clCreateKernel scatter_kernel");

    // allocate host buffers for one chunk
    uint32_t* h_chunk_in  = (uint32_t*)malloc(sizeof(uint32_t) * max_chunk);
    uint32_t* h_chunk_out = (uint32_t*)malloc(sizeof(uint32_t) * max_chunk);
    if (critical_if(!h_chunk_in || !h_chunk_out, "Out of memory")) exit(-1);

    // compute maximum histogram size for a chunk
    size_t hist_elems_max = ((max_chunk + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE) * NUM_BUCKETS;
    size_t hist_bytes_max = sizeof(uint32_t) * hist_elems_max;

    // allocate host buffers for histogram and offsets
    uint32_t* h_hist       = (uint32_t*)malloc(hist_bytes_max);
    uint32_t* h_offsets_gb = (uint32_t*)malloc(hist_bytes_max);
    uint32_t* total_hist   = (uint32_t*)malloc(sizeof(uint32_t) * NUM_BUCKETS);
    uint32_t* bucket_base  = (uint32_t*)malloc(sizeof(uint32_t) * NUM_BUCKETS);
    if (critical_if(!h_hist || !h_offsets_gb || !total_hist || !bucket_base, "Out of memory")) exit(-1);

    // allocate host buffers for incremental merge
    uint32_t* h_sorted_total = (uint32_t*)malloc(sizeof(uint32_t) * N);
    uint32_t* h_merged       = (uint32_t*)malloc(sizeof(uint32_t) * N);
    if (critical_if(!h_sorted_total || !h_merged, "Out of memory")) exit(-1);

    size_t total_filled = 0;

    // create device buffers for one chunk
    cl_mem d_buf_A = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(uint32_t) * max_chunk, NULL, &err);
    CHECK_CL(err, "clCreateBuffer d_buf_A");
    cl_mem d_buf_B = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(uint32_t) * max_chunk, NULL, &err);
    CHECK_CL(err, "clCreateBuffer d_buf_B");
    cl_mem d_hist = clCreateBuffer(context, CL_MEM_READ_WRITE, hist_bytes_max, NULL, &err);
    CHECK_CL(err, "clCreateBuffer d_hist");
    cl_mem d_offsets_gb = clCreateBuffer(context, CL_MEM_READ_WRITE, hist_bytes_max, NULL, &err);
    CHECK_CL(err, "clCreateBuffer d_offsets_gb");

    // ping-pong memory addresses for device buffers
    cl_mem input_buf  = d_buf_A;
    cl_mem output_buf = d_buf_B;

    // compute number of chunks
    size_t num_chunks = (N + max_chunk - 1) / max_chunk;
    for (size_t c = 0; c < num_chunks; ++c) {
        size_t offset = c * max_chunk;
        size_t M = (offset + max_chunk > N) ? (N - offset) : max_chunk;

        // copy chunk from input array to host buffer
        memcpy(h_chunk_in, &input[offset], M * sizeof(uint32_t));

        // write chunk into device input buffer
        CHECK_CL(clEnqueueWriteBuffer(queue, input_buf, CL_TRUE, 0, M * sizeof(uint32_t), h_chunk_in, 0, NULL, NULL),
                 "clEnqueueWriteBuffer input_buf chunk");

        // compute global and group sizes for this chunk
        size_t global_size = ((M + WORKGROUP_SIZE - 1) / WORKGROUP_SIZE) * WORKGROUP_SIZE;
        size_t num_groups  = global_size / WORKGROUP_SIZE;
        size_t hist_elems  = num_groups * NUM_BUCKETS;
        size_t hist_bytes  = hist_elems * sizeof(uint32_t);

        // perform radix passes for this chunk
        for (uint32_t pass = 0; pass < NUM_PASSES; ++pass) {
            uint32_t shift = pass * 8;

            // zero out histogram buffer on device
            {
                void* zero_array = calloc(hist_bytes, 1);
                CHECK_CL(clEnqueueWriteBuffer(queue, d_hist, CL_TRUE, 0, hist_bytes, zero_array, 0, NULL, NULL),
                         "clEnqueueWriteBuffer d_hist zero");
                free(zero_array);
            }

            // set kernel arguments for histogram pass
            CHECK_CL(clSetKernelArg(histogram_kernel, 0, sizeof(cl_mem), &input_buf),  "Arg 0 histogram");
            CHECK_CL(clSetKernelArg(histogram_kernel, 1, sizeof(cl_mem), &d_hist),      "Arg 1 histogram");
            CHECK_CL(clSetKernelArg(histogram_kernel, 2, sizeof(uint32_t), &shift),    "Arg 2 histogram");
            CHECK_CL(clSetKernelArg(histogram_kernel, 3, sizeof(uint32_t), &M),        "Arg 3 histogram");

            // enqueue histogram kernel
            CHECK_CL(clEnqueueNDRangeKernel(queue, histogram_kernel, 1, NULL, &global_size, &WORKGROUP_SIZE, 0, NULL, NULL),
                     "clEnqueueNDRangeKernel histogram_kernel");
            CHECK_CL(clFinish(queue), "clFinish after histogram_kernel");

            // read histogram buffer back to host
            CHECK_CL(clEnqueueReadBuffer(queue, d_hist, CL_TRUE, 0, hist_bytes, h_hist, 0, NULL, NULL),
                     "clEnqueueReadBuffer h_hist");

            // compute total_hist for each bucket
            for (uint32_t b = 0; b < NUM_BUCKETS; ++b) {
                uint32_t sum_b = 0;
                for (size_t g = 0; g < num_groups; ++g) {
                    sum_b += h_hist[g * NUM_BUCKETS + b];
                }
                total_hist[b] = sum_b;
            }
            // compute bucket_base as prefix sums of total_hist
            bucket_base[0] = 0;
            for (uint32_t b = 1; b < NUM_BUCKETS; ++b) {
                bucket_base[b] = bucket_base[b - 1] + total_hist[b - 1];
            }
            // compute offsets_gb for each group and bucket
            for (uint32_t b = 0; b < NUM_BUCKETS; ++b) {
                uint32_t running = bucket_base[b];
                for (size_t g = 0; g < num_groups; ++g) {
                    h_offsets_gb[g * NUM_BUCKETS + b] = running;
                    running += h_hist[g * NUM_BUCKETS + b];
                }
            }
            // write offsets back to device
            CHECK_CL(clEnqueueWriteBuffer(queue, d_offsets_gb, CL_TRUE, 0, hist_bytes, h_offsets_gb, 0, NULL, NULL),
                     "clEnqueueWriteBuffer d_offsets_gb");

            // set kernel arguments for scatter pass
            CHECK_CL(clSetKernelArg(scatter_kernel, 0, sizeof(cl_mem), &input_buf),     "Arg 0 scatter");
            CHECK_CL(clSetKernelArg(scatter_kernel, 1, sizeof(cl_mem), &output_buf),    "Arg 1 scatter");
            CHECK_CL(clSetKernelArg(scatter_kernel, 2, sizeof(cl_mem), &d_offsets_gb),  "Arg 2 scatter");
            CHECK_CL(clSetKernelArg(scatter_kernel, 3, sizeof(uint32_t), &shift),       "Arg 3 scatter");
            CHECK_CL(clSetKernelArg(scatter_kernel, 4, sizeof(uint32_t), &M),           "Arg 4 scatter");

            // enqueue scatter kernel
            CHECK_CL(clEnqueueNDRangeKernel(queue, scatter_kernel, 1, NULL, &global_size, &WORKGROUP_SIZE, 0, NULL, NULL),
                     "clEnqueueNDRangeKernel scatter_kernel");
            CHECK_CL(clFinish(queue), "clFinish after scatter_kernel");

            // swap device buffers for next pass
            {
                cl_mem tmp = input_buf;
                input_buf  = output_buf;
                output_buf = tmp;
            }
        }

        // read sorted chunk from device into host buffer
        CHECK_CL(clEnqueueReadBuffer(queue, input_buf, CL_TRUE, 0, M * sizeof(uint32_t), h_chunk_out, 0, NULL, NULL),
                 "clEnqueueReadBuffer h_chunk_out");

        // incremental merge on CPU: merge into h_sorted_total
        if (c == 0) {
            // first chunk: copy directly to h_sorted_total
            memcpy(h_sorted_total, h_chunk_out, M * sizeof(uint32_t));
            total_filled = M;
        } else {
            // merge previously sorted data with new chunk
            size_t L = total_filled;
            size_t R = M;
            size_t i = 0, j = 0, k = 0;
            while (i < L && j < R) {
                if (h_sorted_total[i] <= h_chunk_out[j]) {
                    h_merged[k++] = h_sorted_total[i++];
                } else {
                    h_merged[k++] = h_chunk_out[j++];
                }
            }
            while (i < L) h_merged[k++] = h_sorted_total[i++];
            while (j < R) h_merged[k++] = h_chunk_out[j++];
            // copy merged result back to h_sorted_total
            memcpy(h_sorted_total, h_merged, (L + R) * sizeof(uint32_t));
            total_filled += R;
        }

        // reset device buffer pointers for next chunk
        input_buf  = d_buf_A;
        output_buf = d_buf_B;
    }

    // copy final sorted data from host buffer to output array
    memcpy(output, h_sorted_total, N * sizeof(uint32_t));

    // release device resources
    clReleaseMemObject(d_buf_A);
    clReleaseMemObject(d_buf_B);
    clReleaseMemObject(d_hist);
    clReleaseMemObject(d_offsets_gb);

    // free host buffers
    free(h_chunk_in);
    free(h_chunk_out);
    free(h_hist);
    free(h_offsets_gb);
    free(total_hist);
    free(bucket_base);
    free(h_sorted_total);
    free(h_merged);
}

}
}