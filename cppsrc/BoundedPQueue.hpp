/**
 * File: BoundedPQueue.h
 * Author: Keith Schwarz (htiek@cs.stanford.edu)
 *
 * An implementation of the bounded priority queue abstraction.
 * A bounded priority queue is in many ways like a regular priority
 * queue.  It stores a collection of elements tagged with a real-
 * valued priority, and allows for access to the element whose
 * priority is the smallest.  However, unlike a regular priority
 * queue, the number of elements in a bounded priority queue has
 * a hard limit that is specified in the constructor.  Whenever an
 * element is added to the bounded priority queue such that the
 * size exceeds the maximum, the element with the highest priority
 * value will be ejected from the bounded priority queue.  In this
 * sense, a bounded priority queue is like a high score table for
 * a video game that stores a fixed number of elements and deletes
 * the least-important entry whenever a new value is inserted.
 *
 * When creating a bounded priority queue, you must specify the
 * maximum number of elements to store in the queue as an argument
 * to the constructor.  For example:
 *
 * BoundedPQueue<int> bpq(15); // Holds up to fifteen values.
 *
 * The maximum size of the bounded priority queue can be obtained
 * using the maxSize() function, as in
 *
 * size_t k = bpq.maxSize();
 *
 * Beyond these restrictions, the bounded priority queue behaves
 * similarly to other containers.  You can query its size using
 * size() and check whether it is empty using empty().  You
 * can enqueue an element into the bounded priority queue by
 * writing
 *
 * bpq.enqueue(elem, priority);
 *
 * Note that after enqueuing the element, there is no guarantee
 * that the value will actually be in the queue.  If the queue
 * is full and the new element's priority exceeds the largest
 * priority in the container, it will not be added.
 *
 * You can dequeue elements from a bounded priority queue using
 * the dequeueMin() function, as in
 *
 * int val = bpq.dequeueMin();
 *
 * The bounded priority queue also allows you to query the min
 * and max priorities of the values in the queue.  These values
 * can be queried using the best() and worst() functions, which
 * return the smallest and largest priorities in the queue,
 * respectively.
 */

#ifndef BOUNDED_PQUEUE_INCLUDED
#define BOUNDED_PQUEUE_INCLUDED

#include <map>
#include <algorithm>
#include <limits>
#include <stdlib.h>



template <typename T, typename dist_type>
class BoundedPQueue {
public:
    // Constructor: BoundedPQueue(size_t maxSize);
    // Usage: BoundedPQueue<int> bpq(15);
    // --------------------------------------------------
    // Constructs a new, empty BoundedPQueue with
    // maximum size equal to the constructor argument.
    ///
    explicit BoundedPQueue(size_t maxSize);
    
    // void enqueue(const T& value, dist_type priority);
    // Usage: bpq.enqueue("Hi!", 2.71828);
    // --------------------------------------------------
    // Enqueues a new element into the BoundedPQueue with
    // the specified priority. If this overflows the maximum
    // size of the queue, the element with the highest
    // priority will be deleted from the queue. Note that
    // this might be the element that was just added.
    void enqueue(const T& value, dist_type priority);
    
    // T dequeueMin();
    // Usage: int val = bpq.dequeueMin();
    // --------------------------------------------------
    // Returns the element from the BoundedPQueue with the
    // smallest priority value, then removes that element
    // from the queue.
    T dequeueMin();
    
    // size_t size() const;
    // bool empty() const;
    // Usage: while (!bpq.empty()) { ... }
    // --------------------------------------------------
    // Returns the number of elements in the queue and whether
    // the queue is empty, respectively.
    size_t size() const;
    bool empty() const;
    
    // size_t maxSize() const;
    // Usage: size_t queueSize = bpq.maxSize();
    // --------------------------------------------------
    // Returns the maximum number of elements that can be
    // stored in the queue.
    size_t maxSize() const;
    
    // dist_type best() const;
    // dist_type worst() const;
    // Usage: dist_type highestPriority = bpq.worst();
    // --------------------------------------------------
    // best() returns the smallest priority of an element
    // stored in the container (i.e. the priority of the
    // element that will be dequeued first using dequeueMin).
    // worst() returns the largest priority of an element
    // stored in the container.  If an element is enqueued
    // with a priority above this value, it will automatically
    // be deleted from the queue.  Both functions return
    // numeric_limits<dist_type>::infinity() if the queue is
    // empty.
    dist_type best()  const;
    dist_type worst() const;
    
private:
    // This class is layered on top of a multimap mapping from priorities
    // to elements with those priorities.
    std::multimap<dist_type, T> elems;
    size_t maximumSize;
};

/** BoundedPQueue class implementation details */

template <typename T, typename dist_type>
BoundedPQueue<T, dist_type>::BoundedPQueue(size_t maxSize) {
    maximumSize = maxSize;
}

// enqueue adds the element to the map, then deletes the last element of the
// map if there size exceeds the maximum size.
template <typename T, typename dist_type>
void BoundedPQueue<T, dist_type>::enqueue(const T& value, dist_type priority) {
    // Add the element to the collection.
    if (size() == maxSize() && priority > worst())
        return;
    elems.emplace(priority, value);
    
    // If there are too many elements in the queue, drop off the last one.
    if (size() > maxSize()) {
        typename std::multimap<dist_type, T>::iterator last = elems.end();
        --last; // Now points to highest-priority element
        elems.erase(last);
    }
}

// dequeueMin copies the lowest element of the map (the one pointed at by
// begin()) and then removes it.
template <typename T, typename dist_type>
T BoundedPQueue<T, dist_type>::dequeueMin() {
    // Copy the best value.
    T result = elems.begin()->second;
    
    // Remove it from the map.
    elems.erase(elems.begin());
    
    return result;
}

// size() and empty() call directly down to the underlying map.
template <typename T, typename dist_type>
size_t BoundedPQueue<T, dist_type>::size() const {
    return elems.size();
}

template <typename T, typename dist_type>
bool BoundedPQueue<T, dist_type>::empty() const {
    return elems.empty();
}

// maxSize just returns the appropriate data member.
template <typename T, typename dist_type>
size_t BoundedPQueue<T, dist_type>::maxSize() const {
    return maximumSize;
}

// The best() and worst() functions check if the queue is empty,
// and if so return infinity.
template <typename T, typename dist_type>
dist_type BoundedPQueue<T, dist_type>::best() const {
    return empty()? std::numeric_limits<dist_type>::infinity() : elems.begin()->first;
}

template <typename T, typename dist_type>
dist_type BoundedPQueue<T, dist_type>::worst() const {
    return empty()? std::numeric_limits<dist_type>::infinity() : elems.rbegin()->first;
}

#endif // BOUNDED_PQUEUE_INCLUDED
