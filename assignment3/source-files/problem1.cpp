#include <chrono>
#include <cmath>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits.h>
#include <omp.h>
#include <tuple>
using std::cout;
using std::endl;

using std::chrono::duration_cast;
using HR = std::chrono::high_resolution_clock;
using HRTimer = HR::time_point;
using std::chrono::microseconds;
using std::chrono::milliseconds;

const uint64_t TIMESTEPS = 100;

const double W_OWN = (1.0 / 7.0);
const double W_NEIGHBORS = (1.0 / 7.0);

const uint64_t NX = 258; // 256 interior points + 2 boundary points
const uint64_t NY = 258;
const uint64_t NZ = 258;
const uint64_t TOTAL_SIZE = NX * NY * NZ;

const static double EPSILON = 1e-9;

#define ASSERT_APPROX_EQ(a, b) assert(std::abs((a) - (b)) < EPSILON);

std::tuple<std::chrono::milliseconds::rep, double, double> profile_kernel(const auto &compute_kernel)
{
  auto *grid1 = new double[TOTAL_SIZE];
  std::fill_n(grid1, TOTAL_SIZE, 0.0);
  grid1[(NX / 2) * NY * NZ + (NY / 2) * NZ + (NZ / 2)] = 100.0;

  auto *grid2 = new double[TOTAL_SIZE];
  std::fill_n(grid2, TOTAL_SIZE, 0.0);
  grid2[(NX / 2) * NY * NZ + (NY / 2) * NZ + (NZ / 2)] = 100.0;

  double *current_grid = grid1;
  double *next_grid = grid2;

  auto start = HR::now();
  for (int t = 0; t < TIMESTEPS; t++)
  {
    compute_kernel(current_grid, next_grid);
    std::swap(current_grid, next_grid);
  }
  auto end = HR::now();
  auto duration = duration_cast<milliseconds>(end - start).count();

  double final = current_grid[(NX / 2) * NY * NZ + (NY / 2) * NZ + (NZ / 2)];

  double total_sum = 0.0;
  for (size_t i = 0; i < TOTAL_SIZE; i++)
  {
    total_sum += current_grid[i];
  }

  delete[] grid1;
  delete[] grid2;

  return {duration, final, total_sum};
}

void stencil_3d_7pt(const double *curr, double *next)
{
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
      for (int k = 1; k < NZ - 1; ++k)
      {
        double neighbors_sum = 0.0;
        neighbors_sum += curr[(i + 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[(i - 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j + 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j - 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k + 1)];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_vectorised(const double *curr, double *next)
{
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
#pragma omp simd
      for (int k = 1; k < NZ - 1; ++k)
      {
        double neighbors_sum = 0.0;
        neighbors_sum += curr[(i + 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[(i - 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j + 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j - 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k + 1)];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_single_loop(const double *curr, double *next)
{
#pragma omp simd // does not vectorise
  for (int cnt = 0; cnt < (NX) * (NY) * (NZ); ++cnt)
  {
    int k = cnt % NZ;
    int j = (cnt / NZ) % NY;
    int i = ((cnt / NZ) / NY) % NX;
    if (i == 0 || i == NX - 1)
      continue;
    if (j == 0 || j == NY - 1)
      continue;
    if (k == 0 || k == NZ - 1)
      continue;

    double neighbors_sum = 0.0;
    neighbors_sum += curr[(i + 1) * NY * NZ + j * NZ + k];
    neighbors_sum += curr[(i - 1) * NY * NZ + j * NZ + k];
    neighbors_sum += curr[i * NY * NZ + (j + 1) * NZ + k];
    neighbors_sum += curr[i * NY * NZ + (j - 1) * NZ + k];
    neighbors_sum += curr[i * NY * NZ + j * NZ + (k + 1)];
    neighbors_sum += curr[i * NY * NZ + j * NZ + (k - 1)];

    next[i * NY * NZ + j * NZ + k] =
        W_OWN * curr[i * NY * NZ + j * NZ + k] +
        W_NEIGHBORS * neighbors_sum;
  }
}

void stencil_3d_7pt_unrolled4(const double *curr, double *next)
{
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
      int k;

      for (k = 1; k < NZ - 5; k += 4)
      {

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k)] +
              curr[i * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k - 1)];
          next[i * NY * NZ + j * NZ + (k)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k)] +
              W_NEIGHBORS * neighbors_sum;
        }

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k)];
          next[i * NY * NZ + j * NZ + (k + 1)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 1)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 1)];
          next[i * NY * NZ + j * NZ + (k + 2)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 2)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 4)] +
              curr[i * NY * NZ + j * NZ + (k + 2)];
          next[i * NY * NZ + j * NZ + (k + 3)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 3)] +
              W_NEIGHBORS * neighbors_sum;
        }
      }

      for (; k < NZ - 1; ++k)
      {
        double neighbors_sum =
            curr[(i + 1) * NY * NZ + j * NZ + k] +
            curr[(i - 1) * NY * NZ + j * NZ + k] +
            curr[i * NY * NZ + (j + 1) * NZ + k] +
            curr[i * NY * NZ + (j - 1) * NZ + k] +
            curr[i * NY * NZ + j * NZ + (k + 1)] +
            curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_loop_split_vectorised(const double *curr, double *next)
{
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
#pragma omp simd
      for (int k = 1; k < NZ - 1; ++k)
      {
        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * (curr[i * NY * NZ + j * NZ + (k + 1)] + curr[i * NY * NZ + j * NZ + (k - 1)]);
      }
#pragma omp simd
      for (int k = 1; k < NZ - 1; ++k)
      {
        next[i * NY * NZ + j * NZ + k] +=
            W_NEIGHBORS * (curr[i * NY * NZ + (j + 1) * NZ + k] + curr[i * NY * NZ + (j - 1) * NZ + k]);
      }

#pragma omp simd
      for (int k = 1; k < NZ - 1; ++k)
      {
        next[i * NY * NZ + j * NZ + k] +=
            W_NEIGHBORS * (curr[(i + 1) * NY * NZ + j * NZ + k] + curr[(i - 1) * NY * NZ + j * NZ + k]);
      }
    }
  }
}

void stencil_3d_7pt_parallelised(const double *curr, double *next)
{
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (int i = 1; i < NX - 1; ++i)
  {

    for (int j = 1; j < NY - 1; ++j)
    {

      for (int k = 1; k < NZ - 1; ++k)
      {

        double neighbors_sum = 0.0;
        neighbors_sum += curr[(i + 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[(i - 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j + 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j - 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k + 1)];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_parallelised_unrolled4(const double *curr, double *next)
{
#pragma omp parallel for collapse(2) schedule(dynamic)
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
      int k;

      for (k = 1; k < NZ - 5; k += 4)
      {

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k)] +
              curr[i * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k - 1)];
          next[i * NY * NZ + j * NZ + (k)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k)] +
              W_NEIGHBORS * neighbors_sum;
        }

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k)];
          next[i * NY * NZ + j * NZ + (k + 1)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 1)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 1)];
          next[i * NY * NZ + j * NZ + (k + 2)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 2)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 4)] +
              curr[i * NY * NZ + j * NZ + (k + 2)];
          next[i * NY * NZ + j * NZ + (k + 3)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 3)] +
              W_NEIGHBORS * neighbors_sum;
        }
      }

      for (; k < NZ - 1; ++k)
      {
        double neighbors_sum =
            curr[(i + 1) * NY * NZ + j * NZ + k] +
            curr[(i - 1) * NY * NZ + j * NZ + k] +
            curr[i * NY * NZ + (j + 1) * NZ + k] +
            curr[i * NY * NZ + (j - 1) * NZ + k] +
            curr[i * NY * NZ + j * NZ + (k + 1)] +
            curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_vectorised_parallelised(const double *curr, double *next)
{
#pragma omp parallel for collapse(1) schedule(dynamic)
  for (int i = 1; i < NX - 1; ++i)
  {

    for (int j = 1; j < NY - 1; ++j)
    {
#pragma omp simd
      for (int k = 1; k < NZ - 1; ++k)
      {

        double neighbors_sum = 0.0;
        neighbors_sum += curr[(i + 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[(i - 1) * NY * NZ + j * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j + 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + (j - 1) * NZ + k];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k + 1)];
        neighbors_sum += curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

void stencil_3d_7pt_vectorised_parallelised_unrolled4(const double *curr, double *next)
{
#pragma omp parallel for collapse(1) schedule(dynamic)
  for (int i = 1; i < NX - 1; ++i)
  {
    for (int j = 1; j < NY - 1; ++j)
    {
      int k;
#pragma omp simd
      for (k = 1; k < NZ - 5; k += 4)
      {

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k)] +
              curr[i * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k - 1)];
          next[i * NY * NZ + j * NZ + (k)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k)] +
              W_NEIGHBORS * neighbors_sum;
        }

        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 1)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 1)] +
              curr[i * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k)];
          next[i * NY * NZ + j * NZ + (k + 1)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 1)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 2)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 2)] +
              curr[i * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 1)];
          next[i * NY * NZ + j * NZ + (k + 2)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 2)] +
              W_NEIGHBORS * neighbors_sum;
        }
        {
          double neighbors_sum =
              curr[(i + 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[(i - 1) * NY * NZ + j * NZ + (k + 3)] +
              curr[i * NY * NZ + (j + 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + (j - 1) * NZ + (k + 3)] +
              curr[i * NY * NZ + j * NZ + (k + 4)] +
              curr[i * NY * NZ + j * NZ + (k + 2)];
          next[i * NY * NZ + j * NZ + (k + 3)] =
              W_OWN * curr[i * NY * NZ + j * NZ + (k + 3)] +
              W_NEIGHBORS * neighbors_sum;
        }
      }
      int next_itr = k;
#pragma omp simd
      for (int k = next_itr; k < NZ - 1; ++k)
      {
        double neighbors_sum =
            curr[(i + 1) * NY * NZ + j * NZ + k] +
            curr[(i - 1) * NY * NZ + j * NZ + k] +
            curr[i * NY * NZ + (j + 1) * NZ + k] +
            curr[i * NY * NZ + (j - 1) * NZ + k] +
            curr[i * NY * NZ + j * NZ + (k + 1)] +
            curr[i * NY * NZ + j * NZ + (k - 1)];

        next[i * NY * NZ + j * NZ + k] =
            W_OWN * curr[i * NY * NZ + j * NZ + k] +
            W_NEIGHBORS * neighbors_sum;
      }
    }
  }
}

int main()
{

  auto [duration_scalar, value_scalar, total_sum_scalar] = profile_kernel(stencil_3d_7pt);

  cout << "Scalar kernel time: " << duration_scalar << " ms" << endl;
  cout << "Final value at center: " << value_scalar << "\n";
  cout << "Total sum : " << total_sum_scalar << "\n";
  cout << "\n";

  auto [duration_vector, value_vector, total_sum_vector] = profile_kernel(stencil_3d_7pt_vectorised);

  cout << "Vector kernel time: " << duration_vector << " ms" << endl;
  cout << "Final value at center: " << value_vector << "\n";
  cout << "Total sum : " << total_sum_vector << "\n";
  cout << "\n";

  auto [duration_single_loop, value_single_loop, total_sum_single_loop] = profile_kernel(stencil_3d_7pt_single_loop);

  cout << "Single Loop kernel time: " << duration_single_loop << " ms" << endl;
  cout << "Final value at center: " << value_single_loop << "\n";
  cout << "Total sum : " << total_sum_single_loop << "\n";
  cout << "\n";

  auto [duration_loop_split_vectorised, value_loop_split_vectorised, total_sum_loop_split_vectorised] = profile_kernel(stencil_3d_7pt_loop_split_vectorised);

  cout << "Loop split vectorised kernel time: " << duration_loop_split_vectorised << " ms" << endl;
  cout << "Final value at center: " << value_loop_split_vectorised << "\n";
  cout << "Total sum : " << total_sum_loop_split_vectorised << "\n";
  cout << "\n";

  auto [duration_loop_unroll, value_loop_unroll, total_sum_loop_unroll] = profile_kernel(stencil_3d_7pt_unrolled4);

  cout << "Loop unroll kernel time: " << duration_loop_unroll << " ms" << endl;
  cout << "Final value at center: " << value_loop_unroll << "\n";
  cout << "Total sum : " << total_sum_loop_unroll << "\n";
  cout << "\n";

  auto [duration_parallel, value_parallel, total_sum_parallel] = profile_kernel(stencil_3d_7pt_parallelised);

  cout << "Parallel kernel time: " << duration_parallel << " ms" << endl;
  cout << "Final value at center: " << value_parallel << "\n";
  cout << "Total sum : " << total_sum_parallel << "\n";
  cout << "\n";

  auto [duration_parallel_unrolled, value_parallel_unrolled, total_sum_parallel_unrolled] = profile_kernel(stencil_3d_7pt_parallelised_unrolled4);

  cout << "Parallel unrolled kernel time: " << duration_parallel_unrolled << " ms" << endl;
  cout << "Final value at center: " << value_parallel_unrolled << "\n";
  cout << "Total sum : " << total_sum_parallel_unrolled << "\n";
  cout << "\n";

  auto [duration_vector_parallel, value_vector_parallel, total_sum_vector_parallel] = profile_kernel(stencil_3d_7pt_vectorised_parallelised);

  cout << "Vector and Parallel kernel time: " << duration_vector_parallel << " ms" << endl;
  cout << "Final value at center: " << value_vector_parallel << "\n";
  cout << "Total sum : " << total_sum_vector_parallel << "\n";
  cout << "\n";

  auto [duration_vector_parallel_unrolled, value_vector_parallel_unrolled, total_sum_vector_parallel_unrolled] = profile_kernel(stencil_3d_7pt_vectorised_parallelised_unrolled4);

  cout << "Vector Parallel and Unroll kernel time: " << duration_vector_parallel_unrolled << " ms" << endl;
  cout << "Final value at center: " << value_vector_parallel_unrolled << "\n";
  cout << "Total sum : " << total_sum_vector_parallel_unrolled << "\n";
  cout << "\n";

  // Assert on final value at the center and the total sum for the other variants
  // that you come up with. We are dealing with doubles.

  auto value_ref = value_scalar;
  auto total_ref = total_sum_scalar;

  // Vectorized kernel
  ASSERT_APPROX_EQ(value_vector, value_ref);
  ASSERT_APPROX_EQ(total_sum_vector, total_ref);

  // Single loop kernel
  ASSERT_APPROX_EQ(value_single_loop, value_ref);
  ASSERT_APPROX_EQ(total_sum_single_loop, total_ref);

  // Loop split vectorized
  ASSERT_APPROX_EQ(value_loop_split_vectorised, value_ref);
  ASSERT_APPROX_EQ(total_sum_loop_split_vectorised, total_ref);

  // Loop unrolled
  ASSERT_APPROX_EQ(value_loop_unroll, value_ref);
  ASSERT_APPROX_EQ(total_sum_loop_unroll, total_ref);

  // Parallel kernel
  ASSERT_APPROX_EQ(value_parallel, value_ref);
  ASSERT_APPROX_EQ(total_sum_parallel, total_ref);

  // Vector + Parallel kernel
  ASSERT_APPROX_EQ(value_vector_parallel, value_ref);
  ASSERT_APPROX_EQ(total_sum_vector_parallel, total_ref);

  // Vector + Parallel + unrolled kernel
  ASSERT_APPROX_EQ(value_vector_parallel_unrolled, value_ref);
  ASSERT_APPROX_EQ(total_sum_vector_parallel_unrolled, total_ref);

  return EXIT_SUCCESS;
}

