#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <unordered_set>
#include <iomanip>
#include "cache.cpp"

void process_entry_I(uint64_t addr, L2Cache &l2, L3Cache &l3, Prefetcher &prefetcher)
{
    uint64_t l2_result = l2.access(addr, false);
    if (l2_result != HIT)
    {
        l2.misses++;
        uint64_t evicted_addr_l3 = INVALID;
        uint64_t l3_result = l3.access(addr);

        if (l3_result == INVALID_MISS)
        {
            l3.misses++;
        }
        else if (l3_result == FULL)
        {
            evicted_addr_l3 = l3.run_replacement(addr);
            l2.evict(evicted_addr_l3);
            l3.misses++;
        }

        if (l2_result != HIT)
        {
            l2_result = l2.access(addr);
            if (l2_result == FULL)
            {
                l2.run_replacement(addr);
            }
        }
    }
    uint64_t prefetch_addr = INVALID;
    uint64_t result = prefetcher.check_for_prefetch(addr, prefetch_addr);
    if (result == PREFETCH)
    {

        addr = prefetch_addr;
        uint64_t l2_result = l2.access(addr, false, true);
        if (l2_result != HIT)
        {
            prefetcher.prefetches_l2++;
            uint64_t evicted_addr_l3 = INVALID;
            uint64_t l3_result = l3.access(addr, true, true);
            if (l3_result == INVALID_MISS)
            {
                prefetcher.prefetches_l3++;
            }
            if (l3_result == FULL)
            {
                prefetcher.prefetches_l3++;
                evicted_addr_l3 = l3.run_replacement(addr, true);
                l2.evict(evicted_addr_l3);
            }

            if (l2_result != HIT)
            {
                l2_result = l2.access(addr, true, true);
                if (l2_result == FULL)
                {
                    l2.run_replacement(addr, true);
                }
            }
        }
    }
}

void process_entry_II(uint64_t addr, uint64_t pc, L2Cache &l2, L3Cache &l3, Prefetcher &prefetcher)
{
    uint64_t l2_result = l2.access(addr, false);
    if (l2_result != HIT)
    {
        l2.misses++;
        uint64_t evicted_addr_l3 = INVALID;
        uint64_t l3_result = l3.access(addr);

        if (l3_result == INVALID_MISS)
        {
            l3.misses++;
        }
        else if (l3_result == FULL)
        {
            evicted_addr_l3 = l3.run_replacement(addr);
            l2.evict(evicted_addr_l3);
            l3.misses++;
        }

        if (l2_result != HIT)
        {
            l2_result = l2.access(addr);
            if (l2_result == FULL)
            {
                l2.run_replacement(addr);
            }
        }
    }
    uint64_t prefetch_addr = INVALID;
    uint64_t result = prefetcher.check_for_prefetch(addr, prefetch_addr, pc, true);
    if (result == PREFETCH)
    {

        addr = prefetch_addr;
        uint64_t l2_result = l2.access(addr, false, true);
        if (l2_result != HIT)
        {
            prefetcher.prefetches_l2++;
            uint64_t evicted_addr_l3 = INVALID;
            uint64_t l3_result = l3.access(addr, true, true);
            if (l3_result == INVALID_MISS)
            {
                prefetcher.prefetches_l3++;
            }
            if (l3_result == FULL)
            {
                prefetcher.prefetches_l3++;
                evicted_addr_l3 = l3.run_replacement(addr, true);
                l2.evict(evicted_addr_l3);
            }

            if (l2_result != HIT)
            {
                l2_result = l2.access(addr, true, true);
                if (l2_result == FULL)
                {
                    l2.run_replacement(addr, true);
                }
            }
        }
    }
}

void process_entry(uint64_t addr, uint64_t pc, std::vector<L2Cache> &l2, std::vector<L3Cache> &l3, std::vector<Prefetcher> &prefetcher)
{
    process_entry_I(addr, l2[0], l3[0], prefetcher[0]);
    process_entry_II(addr, pc, l2[1], l3[1], prefetcher[1]);
}

int main(int argc, char *argv[])
{
    std::vector<L2Cache> l2(2, L2Cache(8, 32, 512 * 1024));   // 512KB, 32B block, 8-way
    std::vector<L3Cache> l3(2, L3Cache(24, 32, 1536 * 1024)); // 1.5MB, 32B block, 24-way
    std::vector<Prefetcher> prefetcher(2, Prefetcher(32));
    int total_accesses = 1;

    FILE *fp;
    char input_name[256];
    int numtraces = atoi(argv[2]);
    char i_or_d;
    char type;
    unsigned long long addr;
    unsigned pc;
    for (int k = 0; k < numtraces; k++)
    {
        sprintf(input_name, "%s_%d", argv[1], k);
        fp = fopen(input_name, "rb");
        assert(fp != NULL);

        while (!feof(fp))
        {
            size_t a = fread(&i_or_d, sizeof(char), 1, fp);
            if (a != 1)
            {
                if (feof(fp))
                {
                    // Reached end of file
                    break;
                }
                else
                {
                    perror("Error reading i_or_d");
                    exit(EXIT_FAILURE);
                }
            }

            size_t b = fread(&type, sizeof(char), 1, fp);
            if (b != 1)
            {
                perror("Error reading type");
                exit(EXIT_FAILURE);
            }

            size_t c = fread(&addr, sizeof(unsigned long long), 1, fp);
            if (c != 1)
            {
                perror("Error reading addr");
                exit(EXIT_FAILURE);
            }

            size_t d = fread(&pc, sizeof(unsigned), 1, fp);
            if (d != 1)
            {
                perror("Error reading pc");
                exit(EXIT_FAILURE);
            }

            // Process the entry
            if (type != 0)
            {
                process_entry(addr, pc, l2, l3, prefetcher);
                total_accesses++;
            }
        }
        fclose(fp);
    }
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "L2 misses Prefetch I: " << l2[0].misses << ", ";
    std::cout << "L2 hits Prefetches I: " << l2[0].prefetch_hits << ", ";
    std::cout << "L2 Prefetches I: " << prefetcher[0].prefetches_l2 << ", ";
    std::cout << "L2 Accuracy I: " << (l2[0].prefetch_hits) / (double)prefetcher[0].prefetches_l2 << ", ";
    std::cout << "L2 Coverage I: " << l2[0].prefetch_hits / (double)(l2[0].prefetch_hits + l2[0].misses) << "\n";

    std::cout << "L3 misses Prefetch I: " << l3[0].misses << ", ";
    std::cout << "L3 hits Prefetches I: " << l3[0].prefetch_hits << ", ";
    std::cout << "L3 Prefetches I: " << prefetcher[0].prefetches_l3 << ", ";
    std::cout << "L3 Coverage I: " << (l3[0].prefetch_hits / (double)(l3[0].prefetch_hits + l3[0].misses)) << "\n"
              << std::endl;

    std::cout << "L2 misses Prefetch II: " << l2[1].misses << ", ";
    std::cout << "L2 hits Prefetches II: " << l2[1].prefetch_hits << ", ";
    std::cout << "L2 Prefetches II: " << prefetcher[1].prefetches_l2 << ", ";
    std::cout << "L2 Accuracy II: " << (l2[1].prefetch_hits) / (double)prefetcher[1].prefetches_l2 << ", ";
    std::cout << "L2 Coverage II: " << l2[1].prefetch_hits / (double)(l2[1].prefetch_hits + l2[1].misses) << "\n";

    std::cout << "L3 misses Prefetch II: " << l3[1].misses << ", ";
    std::cout << "L3 hits Prefetches II: " << l3[1].prefetch_hits << ", ";
    std::cout << "L3 Prefetches II: " << prefetcher[1].prefetches_l3 << ", ";
    std::cout << "L3 Coverage I: " << (l3[1].prefetch_hits / (double)(l3[1].prefetch_hits + l3[1].misses))
              << std::endl;
}
