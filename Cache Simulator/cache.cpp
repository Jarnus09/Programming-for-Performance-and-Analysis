#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include <unordered_set>
#include <map>
#include <set>

#define INVALID (UINT64_MAX)

#define FULL (2)
#define HIT (1)
#define INVALID_MISS (0)

#define HIT_AGE (0)
#define NEW_AGE (2)
#define EVICT_AGE (3)
#define INVALID_AGE (4)

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

public:
    uint64_t misses = 0;
    uint64_t cold_misses = 0;
    BaseCache(uint64_t a, uint64_t b, uint64_t c) : associativity(a), block_size(b), cache_size(c)
    {
        num_sets = cache_size / (block_size * associativity);
    }
    virtual ~BaseCache() {}
    virtual uint64_t access(uint64_t addr, bool put_block = true) = 0;
    virtual void evict(uint64_t addr) = 0;
};

class L2Cache : public BaseCache
{
private:
    std::vector<std::vector<std::array<uint64_t, 2>>> cache;
    std::vector<uint64_t> max_lru_counter;

    uint64_t find_invalid_or_hit_line(uint64_t set_index, uint64_t tag, bool put_block = true)
    {
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
            {
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
    }

    ~L2Cache() override = default;

    uint64_t access(uint64_t addr, bool put_block = true) override
    {
        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        uint64_t line_status = find_invalid_or_hit_line(set_index, tag, put_block);
        return line_status;
    }

    uint64_t run_replacement(uint64_t addr)
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

    std::vector<std::set<uint64_t>> cold_miss_tracker;

    uint64_t find_invalid_or_hit_line(uint64_t set_index, uint64_t tag, bool put_block = true)
    {
        for (uint64_t i = 0; i < associativity; i++)
        {
            if (cache[set_index][i][0] == tag && cache[set_index][i][1] != INVALID)
            {
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
    }

    ~L3Cache() override = default;
    uint64_t access(uint64_t addr, bool put_block = true) override
    {

        uint64_t set_index = get_set_index(addr);
        uint64_t tag = get_tag(addr);
        uint64_t line_status = find_invalid_or_hit_line(set_index, tag, put_block);
        return line_status;
    }
    uint64_t run_replacement(uint64_t addr)
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
