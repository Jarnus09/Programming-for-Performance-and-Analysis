#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <x86intrin.h>
#include <immintrin.h>
#include <assert.h>
using std::cout;
using std::endl;

using std::chrono::duration_cast;
using HR = std::chrono::high_resolution_clock;
using HRTimer = HR::time_point;
using std::chrono::microseconds;
using std::chrono::milliseconds;

const uint32_t NX = 128;
const uint32_t NY = 128;
const uint32_t NZ = 128;
const uint64_t TOTAL_SIZE = (NX * NY * NZ);
constexpr size_t ALIGN = 32;
#define SSE_WIDTH_BITS (128)
#define AVX2_WIDTH_BITS (256)

const uint32_t N_ITERATIONS = 100;
const uint64_t INITIAL_VAL = 1000000;

long compute_checksum(const uint64_t *grid)
{
  uint64_t sum = 0;
  for (int i = 1; i < (NX - 1); i++)
  {
    for (int j = 0; j < NY; j++)
    {
      for (int k = 0; k < NZ; k++)
      {
        sum += grid[i * NY * NZ + j * NZ + k];
      }
    }
  }
  return sum;
}

std::tuple<std::chrono::milliseconds::rep, uint64_t> profile_gradient(const auto &compute_gradient)
{

  auto *i_grid = static_cast<uint64_t *>(aligned_alloc(ALIGN, TOTAL_SIZE * sizeof(uint64_t)));
  auto *o_grid1 = static_cast<uint64_t *>(aligned_alloc(ALIGN, TOTAL_SIZE * sizeof(uint64_t)));

  for (int i = 0; i < NX; i++)
  {
    for (int j = 0; j < NY; j++)
    {
      for (int k = 0; k < NZ; k++)
      {
        i_grid[i * NY * NZ + j * NZ + k] = INITIAL_VAL + i + 2 * j + 3 * k;
      }
    }
  }

  std::fill_n(o_grid1, TOTAL_SIZE, 0);

  auto start = HR::now();
  for (int iter = 0; iter < N_ITERATIONS; ++iter)
  {
    compute_gradient(i_grid, o_grid1);
  }
  auto end = HR::now();
  auto duration = duration_cast<milliseconds>(end - start).count();

  uint64_t checksum = compute_checksum(o_grid1);

  free(i_grid);
  free(o_grid1);

  return std::make_pair(duration, checksum);
}

void scalar_3d_gradient(const uint64_t *A, uint64_t *B)
{
  const uint64_t stride_i = (NY * NZ);
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; ++k)
      {
        uint64_t base_idx = (i * NY * NZ) + j * NZ + k;
        // A[i+1, j, k]
        int A_right = A[base_idx + stride_i];
        // A[i-1, j, k]
        int A_left = A[base_idx - stride_i];
        B[base_idx] = A_right - A_left;
      }
    }
  }
}

void sse4_3d_gradient(const uint64_t *A, uint64_t *B)
{
  const uint64_t stride_i = (NY * NZ);
  const int stride_vec = SSE_WIDTH_BITS / (sizeof(long long) * CHAR_BIT);
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; k += stride_vec)
      {
        uint64_t base_idx = (i * NY * NZ) + j * NZ + k;
        // A[i+1, j, k]
        __m128i tmp0 = _mm_load_si128((__m128i *)&A[base_idx + stride_i]);
        // A[i-1, j, k]
        __m128i tmp1 = _mm_load_si128((__m128i *)&A[base_idx - stride_i]);

        __m128i tmp2 = _mm_sub_epi64(tmp0, tmp1);
        _mm_storeu_si128((__m128i *)&B[base_idx], tmp2);
      }
    }
  }
}

void avx2_3d_gradient(const uint64_t *A, uint64_t *B)
{
  const uint64_t stride_i = (NY * NZ);
  const int stride_vec = AVX2_WIDTH_BITS / (sizeof(long long) * CHAR_BIT);
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 0; j < NY; ++j)
    {
      for (int k = 0; k < NZ; k += stride_vec)
      {
        uint64_t base_idx = (i * NY * NZ) + j * NZ + k;
        // A[i+1, j, k]
        __m256i tmp0 = _mm256_load_si256((__m256i *)&A[base_idx + stride_i]);
        // A[i-1, j, k]
        __m256i tmp1 = _mm256_load_si256((__m256i *)&A[base_idx - stride_i]);

        __m256i tmp2 = _mm256_sub_epi64(tmp0, tmp1);
        _mm256_storeu_si256((__m256i *)&B[base_idx], tmp2);
      }
    }
  }
}

int main()
{
  auto *i_grid = new uint64_t[TOTAL_SIZE];
  for (int i = 0; i < NX; i++)
  {
    for (int j = 0; j < NY; j++)
    {
      for (int k = 0; k < NZ; k++)
      {
        i_grid[i * NY * NZ + j * NZ + k] = (INITIAL_VAL + i +
                                            2 * j + 3 * k);
      }
    }
  }

  auto *o_grid1 = new uint64_t[TOTAL_SIZE];
  std::fill_n(o_grid1, TOTAL_SIZE, 0);

  auto start = HR::now();
  for (int iter = 0; iter < N_ITERATIONS; ++iter)
  {
    scalar_3d_gradient(i_grid, o_grid1);
  }
  auto end = HR::now();
  auto duration = duration_cast<milliseconds>(end - start).count();
  cout << "Scalar kernel time (ms): " << duration << "\n";

  // Compare checksum with vector versions
  uint64_t scalar_checksum = compute_checksum(o_grid1);
  cout << "Checksum: " << scalar_checksum << "\n";
  cout << "\n";

  auto [duration_sse4, checksum_sse4] = profile_gradient(sse4_3d_gradient);
  cout << "SSE4 kernel time (ms): " << duration_sse4 << "\n";
  cout << "Checksum: " << checksum_sse4 << "\n";
  assert(checksum_sse4 == scalar_checksum);
  cout << "\n";

  auto [duration_avx2, checksum_avx2] = profile_gradient(avx2_3d_gradient);
  cout << "AVX2 kernel time (ms): " << duration_avx2 << "\n";
  cout << "Checksum: " << checksum_avx2 << "\n";
  assert(checksum_avx2 == scalar_checksum);

  // Assert the checksum for vectors variants

  delete[] i_grid;
  delete[] o_grid1;

  return EXIT_SUCCESS;
}

