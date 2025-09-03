#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <papi.h>

using std::cerr;
using std::cout;
using std::endl;
using std::uint8_t, std::uint16_t, std::uint32_t, std::uint64_t;

#define INP_H (1 << 9)
#define INP_W (1 << 9)
#define INP_D (1 << 9)
#define FIL_H (3)
#define FIL_W (3)
#define FIL_D (3)

/** Cross-correlation without padding */
void cc_3d_no_padding(const uint64_t *input,
                      const uint64_t (*kernel)[FIL_W][FIL_D], uint64_t *result,
                      const uint64_t outputHeight, const uint64_t outputWidth,
                      const uint64_t outputDepth)
{
  for (uint64_t i = 0; i < outputHeight; i++)
  {
    for (uint64_t j = 0; j < outputWidth; j++)
    {
      for (uint64_t k = 0; k < outputDepth; k++)
      {
        uint64_t sum = 0;
        for (uint64_t ki = 0; ki < FIL_H; ki++)
        {
          for (uint64_t kj = 0; kj < FIL_W; kj++)
          {
            for (uint64_t kk = 0; kk < FIL_D; kk++)
            {
              sum += input[(i + ki) * INP_W * INP_D + (j + kj) * INP_D +
                           (k + kk)] *
                     kernel[ki][kj][kk];
            }
          }
        }
        result[i * outputWidth * outputDepth + j * outputDepth + k] += sum;
      }
    }
  }
}

void cc_3d_with_blocking_input(const uint64_t *input,
                               const uint64_t (*kernel)[FIL_W][FIL_D], uint64_t *result,
                               const uint64_t outputHeight, const uint64_t outputWidth,
                               const uint64_t outputDepth, const uint64_t blockHeight,
                               const uint64_t blockWidth, const uint64_t blockDepth)
{
  uint64_t numBlockHeight = (INP_H + blockHeight - 1) / blockHeight;
  uint64_t numBlockWidth = (INP_W + blockWidth - 1) / blockWidth;
  uint64_t numBlockDepth = (INP_D + blockDepth - 1) / blockDepth;

  for (uint64_t bi = 0; bi < numBlockHeight; bi++)
  {
    for (uint64_t bj = 0; bj < numBlockWidth; bj++)
    {
      for (uint64_t bk = 0; bk < numBlockDepth; bk++)
      {
        for (uint64_t i = bi * blockHeight; i < std::min((bi + 1) * blockHeight, (uint64_t)INP_H); i++)
        {
          for (uint64_t j = bj * blockWidth; j < std::min((bj + 1) * blockWidth, (uint64_t)INP_W); j++)
          {
            for (uint64_t k = bk * blockDepth; k < std::min((bk + 1) * blockDepth, (uint64_t)INP_D); k++)
            {
              {
                uint64_t inp_i = i;
                uint64_t inp_j = j;
                uint64_t inp_k = k;

                uint64_t inp = input[inp_i * INP_W * INP_D + inp_j * INP_D + inp_k];

                for (uint64_t ki = 0; ki < FIL_H; ki++)
                {
                  for (uint64_t kj = 0; kj < FIL_W; kj++)
                  {
                    for (uint64_t kk = 0; kk < FIL_D; kk++)
                    {
                      uint64_t out_i = inp_i - ki;
                      uint64_t out_j = inp_j - kj;
                      uint64_t out_k = inp_k - kk;

                      if (out_i < outputHeight && out_j < outputWidth && out_k < outputDepth && out_i >= 0 && out_j >= 0 && out_k >= 0)
                      {
                        result[out_i * outputWidth * outputDepth + out_j * outputDepth + out_k] +=
                            inp * kernel[ki][kj][kk];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

void cc_3d_with_blocking_output(const uint64_t *input,
                                const uint64_t (*kernel)[FIL_W][FIL_D], uint64_t *result,
                                const uint64_t outputHeight, const uint64_t outputWidth,
                                const uint64_t outputDepth, const uint64_t blockHeight,
                                const uint64_t blockWidth, const uint64_t blockDepth)
{
  uint64_t numBlockHeight = (INP_H + blockHeight - 1) / blockHeight;
  uint64_t numBlockWidth = (INP_W + blockWidth - 1) / blockWidth;
  uint64_t numBlockDepth = (INP_D + blockDepth - 1) / blockDepth;

  for (uint64_t bi = 0; bi < numBlockHeight; bi++)
  {
    for (uint64_t bj = 0; bj < numBlockWidth; bj++)
    {
      for (uint64_t bk = 0; bk < numBlockDepth; bk++)
      {
        for (uint64_t i = bi * blockHeight; i < std::min((bi + 1) * blockHeight, outputHeight); i++)
        {
          for (uint64_t j = bj * blockWidth; j < std::min((bj + 1) * blockWidth, outputWidth); j++)
          {
            for (uint64_t k = bk * blockDepth; k < std::min((bk + 1) * blockDepth, outputDepth); k++)
            {
              {
                uint64_t sum = 0;
                uint64_t inp_i = i;
                uint64_t inp_j = j;
                uint64_t inp_k = k;

                for (uint64_t ki = 0; ki < FIL_H; ki++)
                {
                  for (uint64_t kj = 0; kj < FIL_W; kj++)
                  {
                    for (uint64_t kk = 0; kk < FIL_D; kk++)
                    {
                      sum += input[(inp_i + ki) * INP_W * INP_D + (inp_j + kj) * INP_D +
                                   (inp_k + kk)] *
                             kernel[ki][kj][kk];
                    }
                  }
                }
                result[inp_i * outputWidth * outputDepth + inp_j * outputDepth + inp_k] += sum;
              }
            }
          }
        }
      }
    }
  }
}

int main(int argc, char *argv[])
{

  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0] << " <blockHeight> <blockWidth> <blockDepth>" << std::endl;
    return 1;
  }

  uint64_t blockHeight = std::stoull(argv[1]);
  uint64_t blockWidth = std::stoull(argv[2]);
  uint64_t blockDepth = std::stoull(argv[3]);

  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    std::cerr << "PAPI initialization error!" << std::endl;
    return 1;
  }
  long long start_us, end_us;
  uint64_t *input = new uint64_t[INP_H * INP_W * INP_D];
  std::fill_n(input, INP_H * INP_W * INP_D, 1);

  uint64_t filter[FIL_H][FIL_W][FIL_D] = {{{2, 1, 3}, {2, 1, 3}, {2, 1, 3}},
                                          {{2, 1, 3}, {2, 1, 3}, {2, 1, 3}},
                                          {{2, 1, 3}, {2, 1, 3}, {2, 1, 3}}};

  uint64_t outputHeight = INP_H - FIL_H + 1;
  uint64_t outputWidth = INP_W - FIL_W + 1;
  uint64_t outputDepth = INP_D - FIL_D + 1;

  int eventSet = PAPI_NULL;
  if (PAPI_create_eventset(&eventSet) != PAPI_OK)
  {
    std::cerr << "PAPI create event set error!" << std::endl;
    return 1;
  }

  PAPI_add_event(eventSet, PAPI_L1_DCM);
  PAPI_add_event(eventSet, PAPI_L2_DCM);
  PAPI_add_event(eventSet, PAPI_TOT_CYC);

  // Naive
  start_us = PAPI_get_real_usec();
  PAPI_start(eventSet);
  auto *result = new uint64_t[outputHeight * outputWidth * outputDepth]{0};
  cc_3d_no_padding(input, filter, result, outputHeight, outputWidth,
                   outputDepth);

  long long perf_metrics_naive[4];
  PAPI_stop(eventSet, perf_metrics_naive);
  end_us = PAPI_get_real_usec();
  perf_metrics_naive[3] = end_us - start_us;

  // Blocked input

  start_us = PAPI_get_real_usec();
  PAPI_start(eventSet);
  auto *result_blocked_input = new uint64_t[outputHeight * outputWidth * outputDepth]{0};
  cc_3d_with_blocking_input(input, filter, result_blocked_input, outputHeight, outputWidth,
                            outputDepth, blockHeight, blockWidth, blockDepth);

  long long perf_metrics_blocked_input[4];
  PAPI_stop(eventSet, perf_metrics_blocked_input);
  end_us = PAPI_get_real_usec();
  perf_metrics_blocked_input[3] = end_us - start_us;
  // Blocked output
  start_us = PAPI_get_real_usec();
  PAPI_start(eventSet);
  auto *result_blocked_output = new uint64_t[outputHeight * outputWidth * outputDepth]{0};
  cc_3d_with_blocking_output(input, filter, result_blocked_output, outputHeight, outputWidth,
                             outputDepth, blockHeight, blockWidth, blockDepth);

  long long perf_metrics_blocked_output[4];
  PAPI_stop(eventSet, perf_metrics_blocked_output);
  end_us = PAPI_get_real_usec();
  perf_metrics_blocked_output[3] = end_us - start_us;

  cout << "Blocked height: " << blockHeight << ", Blocked width: " << blockWidth << ", Blocked depth: " << blockDepth << endl;
  cout << "Performance metrics (naive): L1 DCM = " << perf_metrics_naive[0]
       << ", L2 DCM = " << perf_metrics_naive[1] << ", TOT_CYC = " << perf_metrics_naive[2] << ", TIME = " << perf_metrics_naive[3] << " us" << endl;

  cout << "Performance metrics (blocked input): L1 DCM = " << perf_metrics_blocked_input[0]
       << ", L2 DCM = " << perf_metrics_blocked_input[1] << ", TOT_CYC = " << perf_metrics_blocked_input[2] << ", TIME = " << perf_metrics_blocked_input[3] << " us" << endl;

  cout << "Performance metrics (blocked output): L1 DCM = " << perf_metrics_blocked_output[0]
       << ", L2 DCM = " << perf_metrics_blocked_output[1] << ", TOT_CYC = " << perf_metrics_blocked_output[2] << ", TIME = " << perf_metrics_blocked_output[3] << " us" << endl;

  for (uint64_t i = 0; i < outputHeight * outputWidth * outputDepth; i++)
  {
    assert(result[i] == result_blocked_input[i]);
    assert(result[i] == result_blocked_output[i]);
  }

  delete[] result;
  delete[] result_blocked_input;
  delete[] result_blocked_output;

  return EXIT_SUCCESS;
}

