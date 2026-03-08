#pragma once
#include <cstdint>
#include <algorithm>

struct RingI16 {
    int16_t* buf = nullptr;
    uint32_t cap = 0;   // number of samples in backing buffer
    uint32_t r = 0;     // read index
    uint32_t w = 0;     // write index

    // Note: cap must be >= 2. We leave one slot empty to distinguish full vs empty.
    void init(int16_t* backing, uint32_t capacity) {
        buf = backing;
        cap = capacity;
        r = w = 0;
    }

    uint32_t size() const {              // samples available
        return (w + cap - r) % cap;
    }

    uint32_t free() const {              // free slots
        return (cap - 1) - size();
    }

    uint32_t capacity() const {
        return cap - 1;   // total capacity
    }

    // Push up to n samples. Returns number pushed.
    uint32_t push(const int16_t* in, uint32_t n) {
        uint32_t can = std::min(n, free());
        for (uint32_t i = 0; i < can; ++i) {
            buf[w] = in[i];
            w = (w + 1) % cap;
        }
        return can;
    }

    // Copy up to n samples to out WITHOUT consuming them. Returns number copied.
    uint32_t peek(int16_t* out, uint32_t n) const {
        uint32_t can = std::min(n, size());
        uint32_t idx = r;
        for (uint32_t i = 0; i < can; ++i) {
            out[i] = buf[idx];
            idx = (idx + 1) % cap;
        }
        return can;
    }

    // Consume (drop) up to n samples. Returns number consumed.
    uint32_t drop(uint32_t n) {
        uint32_t can = std::min(n, size());
        r = (r + can) % cap;
        return can;
    }
};