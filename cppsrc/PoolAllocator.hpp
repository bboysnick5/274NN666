//
//  Allocator.hpp
//  274F16NearestSB
//
//  Created by nick on 1/12/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef PoolAllocator_hpp
#define PoolAllocator_hpp

#include <stdio.h>
#include <cstdlib>
#include <utility>


class PooledAllocator {
    /* We maintain memory alignment to word boundaries by requiring that all
     allocations be in multiples of the machine wordsize.  */
    /* Size of machine word in bytes.  Must be power of 2. */
    /* Minimum number of bytes requested at a time from the system.  Must be
     * multiple of WORDSIZE. */
    
    constexpr static std::size_t WORDSIZE = 8;
    constexpr static std::size_t BLOCKSIZE = 8192;
    
    std::size_t remaining; /* Number of bytes left in current block of storage. */
    void *base;       /* PointNDer to base of current block of storage. */
    void *loc;        /* Current location in block to next allocate memory. */
    
    void internal_init() {
        remaining = 0;
        base = nullptr;
    }
    
    /**
     * Returns a pointer to a piece of new memory of the given size in bytes
     * allocated from the pool.
     */
    void *malloc(const std::size_t req_size, bool exact = true) {
        /* Round size up to a multiple of wordsize.  The following expression
         only works for WORDSIZE that is a power of 2, by masking last bits of
         incremented size to zero.
         */
        const std::size_t size = (req_size + (WORDSIZE - 1)) & ~(WORDSIZE - 1);
        
        /* Check whether a new block must be allocated.  Note that the first word
         of a block is reserved for a pointer to the previous block.
         */
        if (size > remaining) {
            /* Allocate new storage. */
            const std::size_t blockSizeReq = size + sizeof(void *) + (WORDSIZE - 1);
            const std::size_t blocksize = exact ? blockSizeReq :
                ((blockSizeReq > BLOCKSIZE) ? blockSizeReq : BLOCKSIZE);
            
            // use the standard C malloc to allocate memory
            void *m = ::malloc(blocksize);
            if (!m) {
                fprintf(stderr, "Failed to allocate memory.\n");
                return nullptr;
            }
            
            /* Fill first word of new block with pointer to previous block. */
            static_cast<void **>(m)[0] = base;
            base = m;
            
            constexpr static std::size_t SIZE_OF_VOID_STAR = sizeof(void*);
            remaining = blocksize - SIZE_OF_VOID_STAR;
            loc = static_cast<char *>(m) + SIZE_OF_VOID_STAR;
        }
        void *rloc = loc;
        loc = static_cast<char *>(loc) + size;
        remaining -= size;
        
        return rloc;
    }
    
public:
    
    /**
     Default constructor. Initializes a new pool.
     */
    PooledAllocator() { internal_init(); }
    
    PooledAllocator(const PooledAllocator &other) = delete;
    PooledAllocator& operator = (const PooledAllocator &other) & = delete;
    
    PooledAllocator(PooledAllocator &&other) noexcept
    : remaining(other.remaining), base(other.base), loc(other.loc) {
        other.base = nullptr;
    }
    
    PooledAllocator& operator = (PooledAllocator &&other) noexcept {
        if (this != &other) {
            free_all();
            remaining = other.remaining;
            base = other.base;
            loc = other.loc;
            other.base = nullptr;
        }
        return *this;
    }
    
    /**
     * Destructor. Frees all the memory allocated in this pool.
     */
    ~PooledAllocator() { free_all(); }
    
    /** Frees all allocated memory chunks */
    void free_all() {
        while (base) {
            void *prev =
            *(static_cast<void **>(base)); /* Get pointer to prev block. */
            ::free(base);
            base = prev;
        }
        internal_init();
    }
    
    template <typename T, class ...Args>
    void static construct(T* p, Args&&... args) {
        ::new (static_cast<void*>(p)) T(std::forward<Args>(args)...);
    }
    
    template <typename T>
    void static destroy(T *p) {
        p->~T();
    }
    
    template <typename T>
    void destroy_and_free_all() {
        while (base) {
            std::destroy_n(reinterpret_cast<T*>(static_cast<char*>(base)+sizeof(void*)),
                           (static_cast<char*>(loc)- static_cast<char*>(base)+remaining)/sizeof(T));
            void *prev = *(static_cast<void **>(base));
            ::free(base);
            base = prev;
        }
        internal_init();
    }
    
    /**
     * Allocates (using this pool) a generic type T.
     *
     * Params:
     *     count = number of instances to allocate.
     * Returns: pointer (of type T*) to memory buffer
     */
    template <typename T>
    T *allocateExact(const std::size_t count = 1) {
        return static_cast<T *>(this->malloc(sizeof(T) * count));
    }
    
    template <typename T>
    T *allocatePool(const std::size_t count = 1) {
        return static_cast<T *>(this->malloc(sizeof(T) * count, false));
    }
};


#endif /* PoolAllocator_hpp */
