#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <unordered_set>
#include "cache.cpp"

void process_entry_unopt(unsigned long long addr, L2Cache &l2, L3Cache &l3)
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
}
void process_entry_opt(unsigned long long addr, L2Cache &l2, L3CacheOPT &l3)
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
    l3.remove_future_use(addr);
}

void process_entry(unsigned long long addr, std::vector<L2Cache> &l2, std::vector<L3Cache> &l3, std::vector<L3CacheOPT> &l3_opt)
{
    process_entry_unopt(addr, l2[0], l3[0]);
    process_entry_unopt(addr, l2[1], l3[1]);

    process_entry_opt(addr, l2[2], l3_opt[0]);
    process_entry_opt(addr, l2[3], l3_opt[1]);
}

struct TraceEntry
{
    char i_or_d;
    char type;
    unsigned long long addr;
    unsigned pc;
};

int main(int argc, char *argv[])
{
    std::vector<L2Cache> l2(4, L2Cache(8, 32, 512 * 1024));   // 512KB, 32B block, 8-way
    std::vector<L3Cache> l3(2, L3Cache(24, 32, 1536 * 1024)); // 1.5MB, 32B block, 24-way
    std::vector<L3CacheOPT> l3_opt(2, L3CacheOPT(24, 32, 1536 * 1024));
    l3[1] = L3Cache(49152, 32, 1536 * 1024);
    l3_opt[1] = L3CacheOPT(49152, 32, 1536 * 1024);
    uint64_t total_accesses = 0;
    std::vector<TraceEntry> entries;
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
                entries.push_back({i_or_d, type, addr, pc});

                l3_opt[0].add_future_use(addr, total_accesses);
                l3_opt[1].add_future_use(addr, total_accesses);
                total_accesses++;
            }
        }
        fclose(fp);
    }

    for (auto it : entries)
    {
        process_entry(it.addr, l2, l3, l3_opt);
    }

    std::cout << "L3 cold misses Original: " << l3[0].cold_misses << std::endl;
    std::cout << "L3 conflict misses Original (Quad - Age): " << l3[0].misses - l3_opt[1].misses << std::endl;
    std::cout << "L3 conflict misses Original (OPT): " << l3_opt[0].misses - l3_opt[1].misses << std::endl;
    std::cout << "L3 capacity misses Original: " << l3_opt[1].misses - l3[0].cold_misses << std::endl;
    std::cout << std::endl;

    std::cout << "L3 cold misses Quad-Age: " << l3[0].cold_misses << std::endl;
    std::cout << "L3 conflict misses Quad-Age (Quad - Age): " << l3[0].misses - l3[1].misses << std::endl;
    std::cout << "L3 conflict misses Quad-Age (OPT): " << l3_opt[0].misses - l3[1].misses << std::endl;
    std::cout << "L3 capacity misses Quad-Age: " << l3[1].misses - l3[0].cold_misses << std::endl;

    std::cout << std::endl;
}
