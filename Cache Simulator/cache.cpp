#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <set>

#define INVALID (UINT64_MAX)

#define FULL (2)
#define HIT (1)
#define INVALID_MISS (0)
#define VALID (1)

#define HIT_AGE (0)
#define NEW_AGE (2)
#define EVICT_AGE (3)
#define INVALID_AGE (4)

#define PAGE_SIZE (1LL << (12))
#define PREFETCH (1)
#define NO_PREFETCH (0)

constexpr uint64_t stream_table_size = 32;

class BaseCache
{
protected:
    uint64_t associativity;
    uint64_t block_size;
    uint64_t cache_size;
    uint64_t num_sets;
    uint64_t get_set_index(uint64_t addr)
    {
        return (addr / block_size) % num_sets;
    }

    uint64_t get_tag(uint64_t addr)
    {
        return addr / (block_size * num_sets);
    }
    uint64_t get_addr(uint64_t set_index, uint64_t tag)
    {
        return ((tag * num_sets) + set_index) * block_size;
    }

public:
    int64_t misses = 0;
    int64_t cold_misses = 0;
    int64_t prefetch_hits = 0;

    BaseCache(uint64_t a, uint64_t b, uint64_t c) : associativity(a), block_size(b), cache_size(c)
    {
        num_sets = cache_size / (block_size * associativity);
    }

    virtual ~BaseCache() {}
    virtual uint64_t access(uint64_t addr, bool put_block = true, bool prefetch_bit = true) = 0;
    virtual void evict(uint64_t addr) = 0;
};

class L2Cache : public BaseCache
{
private:
    std::vector<std::vector<std::array<uint64_t, 2>>> cache;
    std::vector<uint64_t> max_lru_counter;
    std::vector<std::vector<bool>> prefetch_bits;

    uint64_t find_invalid_or_hit_line(uint64_t set_index, uint64_t tag, bool put_block, bool prefetch_bit)
    {
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
            {
                if (prefetch_bits[set_index][i] == true && prefetch_bit == false)
                {
                    prefetch_hits++;
                    prefetch_bits[set_index][i] = false;
                }
                cache[set_index][i][1] = max_lru_counter[set_index]++;
                return HIT;
            }
        }
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][1] == INVALID)
            {
                if (put_block == false)
                    return INVALID_MISS;
                prefetch_bits[set_index][i] = prefetch_bit;
                cache[set_index][i][0] = tag;
                cache[set_index][i][1] = max_lru_counter[set_index]++;
                return INVALID_MISS;
            }
        }
        return FULL;
    }

public:
    L2Cache(uint64_t a, uint64_t b, uint64_t c) : BaseCache(a, b, c)
    {
        cache.resize(num_sets, std::vector<std::array<uint64_t, 2>>(associativity, {0, INVALID}));
        max_lru_counter.resize(num_sets, 0);
        prefetch_bits.resize(num_sets, std::vector<bool>(associativity, false));
    }

    ~L2Cache() override = default;

    uint64_t access(uint64_t addr, bool put_block = true, bool prefetch_bit = false) override
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        uint64_t line_status = find_invalid_or_hit_line(set_index, tag, put_block, prefetch_bit);
        return line_status;
    }

    uint64_t run_replacement(uint64_t addr, bool prefetch_bit = false)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);

        uint64_t lru_index = 0;
        uint64_t lru_counter = UINT64_MAX;
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][1] != INVALID && cache[set_index][i][1] < lru_counter)
            {
                lru_counter = cache[set_index][i][1];
                lru_index = i;
            }
        }
        uint64_t evicted_tag = cache[set_index][lru_index][0];
        uint64_t evicted_addr = ((evicted_tag * num_sets) + set_index) * block_size;
        prefetch_bits[set_index][lru_index] = prefetch_bit;
        cache[set_index][lru_index][0] = tag;
        cache[set_index][lru_index][1] = max_lru_counter[set_index]++;
        return evicted_addr;
    }

    void evict(uint64_t addr) override
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag)
            {
                cache[set_index][i][1] = INVALID;
                return;
            }
        }
    }
};

class L3Cache : public BaseCache
{
private:
    std::vector<std::vector<std::array<uint64_t, 2>>> cache;
    std::vector<std::vector<bool>> l2_notify;
    std::vector<std::unordered_set<uint64_t>> accessed_tags;
    std::vector<std::vector<bool>> prefetch_bits;

    std::vector<std::set<uint64_t>> cold_miss_tracker;

    uint64_t find_invalid_or_hit_line(uint64_t set_index, uint64_t tag, bool put_block, bool prefetch_bit)
    {
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
            {
                if (prefetch_bits[set_index][i] == true && prefetch_bit == false)
                {
                    prefetch_hits++;
                    prefetch_bits[set_index][i] = false;
                }
                l2_notify[set_index][i] = true;
                cache[set_index][i][1] = HIT_AGE;

                return HIT;
            }
        }
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][1] == INVALID)
            {
                if (put_block == false)
                    return INVALID_MISS;
                check_cold_miss(set_index, tag);
                prefetch_bits[set_index][i] = prefetch_bit;
                l2_notify[set_index][i] = true;
                cache[set_index][i][0] = tag;
                cache[set_index][i][1] = NEW_AGE;
                return INVALID_MISS;
            }
        }
        return FULL;
    }
    void check_cold_miss(uint64_t set_index, uint64_t tag)
    {
        if (cold_miss_tracker[set_index].find(tag) == cold_miss_tracker[set_index].end())
        {
            cold_misses++;
        }
        cold_miss_tracker[set_index].insert(tag);
    }
    void increment_ages(uint64_t set_index)
    {
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][1] != INVALID)
                cache[set_index][i][1]++;
        }
    }

public:
    L3Cache(uint64_t a, uint64_t b, uint64_t c) : BaseCache(a, b, c)
    {
        cache.resize(num_sets, std::vector<std::array<uint64_t, 2>>(associativity, {0, INVALID}));
        l2_notify.resize(num_sets, std::vector<bool>(associativity, false));
        accessed_tags.resize(num_sets);
        cold_miss_tracker.resize(num_sets);
        prefetch_bits.resize(num_sets, std::vector<bool>(associativity, false));
    }

    ~L3Cache() override = default;
    uint64_t access(uint64_t addr, bool put_block = true, bool prefetch_bit = false) override
    {

        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        uint64_t line_status = find_invalid_or_hit_line(set_index, tag, put_block, prefetch_bit);
        return line_status;
    }
    uint64_t run_replacement(uint64_t addr, bool prefetch_bit = false)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        while (1)
        {

            for (uint64_t i = 0; i < associativity; i++)
            {
                if (cache[set_index][i][1] == EVICT_AGE)
                {
                    check_cold_miss(set_index, tag);
                    uint64_t evicted_tag = cache[set_index][i][0];
                    cache[set_index][i][0] = tag;
                    cache[set_index][i][1] = NEW_AGE;
                    prefetch_bits[set_index][i] = prefetch_bit;
                    uint64_t evicted_addr = ((evicted_tag * num_sets) + set_index) * block_size;
                    return evicted_addr;
                }
            }
            increment_ages(set_index);
        }
    }

    void notify(uint64_t addr)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
            {
                l2_notify[set_index][i] = false;
                return;
            }
        }
    }

    uint64_t run_replacement_D(uint64_t addr)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        while (1)
        {
            for (uint64_t i = 0; i < associativity; i++)
            {
                if (cache[set_index][i][1] == EVICT_AGE && (l2_notify[set_index][i] == false))
                {
                    uint64_t evicted_tag = cache[set_index][i][0];
                    cache[set_index][i][0] = tag;
                    cache[set_index][i][1] = NEW_AGE;
                    l2_notify[set_index][i] = true;
                    uint64_t evicted_addr = ((evicted_tag * num_sets) + set_index) * block_size;
                    return evicted_addr;
                }
            }

            for (uint64_t i = 0; i < associativity; i++)
            {
                if (cache[set_index][i][1] == EVICT_AGE)
                {
                    uint64_t evicted_tag = cache[set_index][i][0];
                    cache[set_index][i][0] = tag;
                    cache[set_index][i][1] = NEW_AGE;
                    l2_notify[set_index][i] = true;
                    uint64_t evicted_addr = ((evicted_tag * num_sets) + set_index) * block_size;
                    return evicted_addr;
                }
            }
            increment_ages(set_index);
        }
    }

    void update_age(uint64_t addr)
    {
        uint64_t tag = get_tag(addr);
        uint64_t set_index = get_set_index(addr);
        if (accessed_tags[set_index].find(tag) != accessed_tags[set_index].end())
        {
            for (uint64_t i = 0; i < associativity; i++)
            {
                if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
                {
                    cache[set_index][i][1] = HIT_AGE;
                    return;
                }
            }
        }
    }
    void add_unordered_set(uint64_t addr)
    {
        uint64_t tag = get_tag(addr);
        uint64_t set_index = get_set_index(addr);
        accessed_tags[set_index].insert(tag);
    }
    void remove_unordered_set(uint64_t addr)
    {
        uint64_t tag = get_tag(addr);
        uint64_t set_index = get_set_index(addr);
        accessed_tags[set_index].erase(tag);
    }
    void evict(uint64_t addr) override
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag)
            {
                l2_notify[set_index][i] = false;
                cache[set_index][i][1] = INVALID;
                return;
            }
        }
    }
};

class L3CacheOPT : public BaseCache
{
private:
    std::vector<std::set<uint64_t>> cache;
    std::vector<std::unordered_map<uint64_t, std::queue<uint64_t>>> future_uses;
    std::vector<std::set<std::array<uint64_t, 2>>> max_time_tracker;
    std::vector<std::set<uint64_t>> cold_miss_tracker;

    uint64_t find_invalid_or_hit_line(uint64_t set_index, uint64_t tag, bool put_block = true)
    {

        if (cache[set_index].count(tag))
        {

            return HIT;
        }

        if (cache[set_index].size() < associativity)
        {
            cache[set_index].insert(tag);

            return INVALID_MISS;
        }
        return FULL;
    }
    void check_cold_miss(uint64_t set_index, uint64_t tag)
    {
        if (cold_miss_tracker[set_index].find(tag) == cold_miss_tracker[set_index].end())
        {
            cold_misses++;
        }
        cold_miss_tracker[set_index].insert(tag);
    }

public:
    L3CacheOPT(uint64_t a, uint64_t b, uint64_t c) : BaseCache(a, b, c)
    {
        cache.resize(num_sets);
        future_uses.resize(num_sets);
        cold_miss_tracker.resize(num_sets);
        max_time_tracker.resize(num_sets);
    }
    ~L3CacheOPT() override = default;
    uint64_t access(uint64_t addr, bool put_block = true, bool prefetch_bit = true) override
    {

        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        uint64_t line_status = find_invalid_or_hit_line(set_index, tag, put_block);
        return line_status;
    }
    void add_future_use(uint64_t addr, uint64_t use_index)
    {
        uint64_t tag = get_tag(addr);
        uint64_t set_index = get_set_index(addr);
        future_uses[set_index][tag].push(use_index);
    }
    uint64_t run_replacement(uint64_t addr)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);

        auto top = *max_time_tracker[set_index].rbegin();
        max_time_tracker[set_index].erase(top);
        cache[set_index].erase(top[1]);

        uint64_t evicted_tag = top[1];
        uint64_t evicted_addr = ((evicted_tag * num_sets) + set_index) * block_size;
        cache[set_index].insert(tag);

        return evicted_addr;
    }
    void evict(uint64_t addr) override
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);

        if (cache[set_index].count(tag))
        {
            cache[set_index].erase(tag);
            return;
        }
    }

    void remove_future_use(uint64_t addr)
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        if (future_uses[set_index].find(tag) != future_uses[set_index].end())
        {
            uint64_t time = INVALID;
            auto curr = future_uses[set_index][tag].front();
            if (max_time_tracker[set_index].count({curr, tag}))
                max_time_tracker[set_index].erase({curr, tag});
            future_uses[set_index][tag].pop();
            if (future_uses[set_index][tag].empty())
            {
                future_uses[set_index].erase(tag);
            }
            else
            {
                time = future_uses[set_index][tag].front();
            }
            max_time_tracker[set_index].insert({time, tag});
        }
        else
            assert(false);
    }
};

struct StreamEntry
{
    uint64_t page_id;
    uint64_t last_addr;
    uint64_t stride = 0;
    uint64_t pc = 0;
    uint64_t lru_counter;
    bool check_valid = false;
};

class Prefetcher
{
private:
    std::vector<StreamEntry> stream_table;
    uint64_t block_size;
    uint64_t max_lru_counter;

    uint64_t get_block_aligned_addr(uint64_t addr)
    {
        return (addr / block_size) * block_size;
    }
    uint64_t get_page_id(uint64_t addr)
    {
        return (addr / PAGE_SIZE);
    }

public:
    uint64_t prefetches_l2 = 0;
    uint64_t prefetches_l3 = 0;
    Prefetcher(uint64_t b) : stream_table(stream_table_size), block_size(b), max_lru_counter(1) {};

    uint64_t check_for_prefetch(uint64_t addr, uint64_t &prefetch_addr, uint64_t pc = 0, bool algorithmII = false)
    {
        uint64_t page_id = get_page_id(addr);
        uint64_t block_addr = get_block_aligned_addr(addr);

        for (uint64_t i = 0; i < stream_table_size; i++)
        {
            if (algorithmII == false)
            {
                if (page_id == stream_table[i].page_id && stream_table[i].check_valid == true)
                {
                    if (stream_table[i].last_addr != block_addr)
                    {
                        if (stream_table[i].stride == 0)
                        {

                            stream_table[i].stride = block_addr - stream_table[i].last_addr;
                            stream_table[i].lru_counter = max_lru_counter++;
                            stream_table[i].last_addr = block_addr;
                            return NO_PREFETCH;
                        }
                        else
                        {
                            if (block_addr == stream_table[i].stride + stream_table[i].last_addr)
                            {
                                prefetch_addr = block_addr + stream_table[i].stride;
                                stream_table[i].lru_counter = max_lru_counter++;
                                stream_table[i].last_addr = block_addr;
                                if (get_page_id(prefetch_addr) == page_id)
                                {
                                    return PREFETCH;
                                }
                                else
                                    return NO_PREFETCH;
                            }
                            else
                            {
                                stream_table[i].check_valid = true;
                                stream_table[i].stride = block_addr - stream_table[i].last_addr;
                                stream_table[i].last_addr = block_addr;

                                stream_table[i].page_id = page_id;
                                stream_table[i].lru_counter = max_lru_counter++;
                                return NO_PREFETCH;
                            }
                        }
                    }
                    else
                    {
                        stream_table[i].lru_counter = max_lru_counter++;
                        return NO_PREFETCH;
                    }
                }
            }
            else
            {
                if (pc == stream_table[i].pc && stream_table[i].check_valid == true)
                {
                    if (stream_table[i].last_addr != block_addr)
                    {
                        if (stream_table[i].stride == 0)
                        {
                            stream_table[i].stride = block_addr - stream_table[i].last_addr;
                            stream_table[i].lru_counter = max_lru_counter++;
                            stream_table[i].last_addr = block_addr;
                            return NO_PREFETCH;
                        }
                        else
                        {
                            if (block_addr == stream_table[i].stride + stream_table[i].last_addr)
                            {
                                prefetch_addr = block_addr + stream_table[i].stride;
                                stream_table[i].lru_counter = max_lru_counter++;
                                stream_table[i].last_addr = block_addr;
                                if (get_page_id(prefetch_addr) == page_id)
                                {
                                    return PREFETCH;
                                }
                                else
                                    return NO_PREFETCH;
                            }
                            else
                            {
                                stream_table[i].check_valid = true;
                                stream_table[i].stride = block_addr - stream_table[i].last_addr;
                                stream_table[i].last_addr = block_addr;

                                stream_table[i].page_id = page_id;
                                stream_table[i].lru_counter = max_lru_counter++;
                                return NO_PREFETCH;
                            }
                        }
                    }
                    else
                    {
                        stream_table[i].lru_counter = max_lru_counter++;
                        return NO_PREFETCH;
                    }
                }
            }
        }
        uint64_t min_lru = UINT64_MAX;
        uint64_t replace_index = -1;
        for (uint64_t i = 0; i < stream_table_size; i++)
        {
            if (stream_table[i].check_valid == false)
            {
                stream_table[i].check_valid = true;
                stream_table[i].last_addr = block_addr;
                stream_table[i].stride = 0;
                stream_table[i].page_id = page_id;
                stream_table[i].pc = pc;
                stream_table[i].lru_counter = max_lru_counter++;
                return NO_PREFETCH;
            }
            else if (stream_table[i].lru_counter < min_lru)
            {
                min_lru = stream_table[i].lru_counter;
                replace_index = i;
            }
        }
        stream_table[replace_index].pc = pc;
        stream_table[replace_index].check_valid = true;
        stream_table[replace_index].last_addr = block_addr;
        stream_table[replace_index].stride = 0;
        stream_table[replace_index].page_id = page_id;
        stream_table[replace_index].lru_counter = max_lru_counter++;
        return NO_PREFETCH;
    }
};

