///////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2012-2013 DreamWorks Animation LLC
//
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
//
// Redistributions of source code must retain the above copyright
// and license notice and the following restrictions and disclaimer.
//
// *     Neither the name of DreamWorks Animation nor the names of
// its contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// IN NO EVENT SHALL THE COPYRIGHT HOLDERS' AND CONTRIBUTORS' AGGREGATE
// LIABILITY FOR ALL CLAIMS REGARDLESS OF THEIR BASIS EXCEED US$250.00.
//
///////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////
///
/// @author Ken Museth
///
/// @file LeafManager.h
///
/// @brief Given a tree this class defines a linear array of its leafs
/// and optional auxiliary buffers. This is useful for multi-threading
/// computations over leaf values in a static tree, i.e. voxel values,
/// vs topology, is dynamic. The auxiliary buffers can conventiently
/// be used for temporal integration. Efficient methods are provided for
/// multi-threaded swapping and sync'ing (i.e. copying) of these buffers.

#ifndef OPENVDB_TREE_LEAF_MANAGER_HAS_BEEN_INCLUDED
#define OPENVDB_TREE_LEAF_MANAGER_HAS_BEEN_INCLUDED

#include <iostream>
#include <algorithm> // for std::swap
#include <boost/shared_ptr.hpp>
#include <boost/static_assert.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <openvdb/Types.h>

namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tree {
  
/// @brief Given a tree this class defines a linear array of its leafs
/// and optional auxiliary buffers. This is useful for multi-threading
/// computations over leaf values in a static tree, i.e. voxel values,
/// vs topology, is dynamic. The auxiliary buffers can conventiently
/// be used for temporal integration. Efficient methods are provided for
/// multi-threaded swapping and sync'ing (i.e. copying) of these
/// buffers.
///    
/// @note By definition the buffers residing inside LeafNodes
/// correspond to bufferId = 0, and any (optional) auxiliary buffers are
/// indexed starting from one. 
template<typename TreeT>
class LeafManager
{
  public:
    typedef TreeT                               TreeType;
    typedef typename TreeType::LeafNodeType     LeafType;
    typedef typename LeafType::Buffer           BufferType;
    typedef typename tbb::blocked_range<size_t> RangeType;//index range to leafs
    
    class LeafRange {    
    public:

        class Iterator
        {
        public:
            Iterator(const LeafRange& range, size_t pos) : mRange(range), mPos(pos) {assert(this->isValid());}
            Iterator& operator++() { ++mPos; return *this; }
            /// Return a reference to the leaf to which this iterator is pointing.
            LeafType& operator*() const { return mRange.mLeafManager.leaf(mPos); }
            /// Return a pointer to the leaf to which this iterator is pointing.
            LeafType* operator->() const { return &(this->operator*()); }
            BufferType& buffer(size_t bufferId) { return mRange.mLeafManager.getBuffer(mPos, bufferId); }
            size_t pos() const { return mPos; }
            bool isValid() const { return mPos>=mRange.mBegin && mPos<=mRange.mEnd; }
            /// Return @c true if this iterator is not yet exhausted.
            operator bool() const { return mPos < mRange.mEnd; }
            //bool operator<( const Iterator& other ) const { return mPos < other.mPos; }
            void operator=( const Iterator& other) { mRange = other.mRange; mPos = other.mPos; }
            bool operator!=(const Iterator& other) const
            {
                return  (mPos != other.mPos) || (&mRange != &other.mRange);
            }
            bool operator==(const Iterator& other) const { return !(*this != other); }
            const LeafRange& leafRange() const { return mRange; }
        private:
            const LeafRange& mRange;
            size_t mPos;
        };// end Iterator
        
        LeafRange( size_t begin, size_t end, const LeafManager& leafManager, size_t grainSize=1 ) : 
            mEnd(end), mBegin(begin), mGrainSize(grainSize), mLeafManager(leafManager) {}
        
        Iterator begin() const {return Iterator(*this, mBegin);}
 
        Iterator end() const {return Iterator(*this, mEnd);}
        
        size_t size() const { return mEnd - mBegin; }
        
        size_t grainsize() const { return mGrainSize; }

        const LeafManager& leafManager() const { return mLeafManager; }
        
        bool empty() const {return !(mBegin < mEnd);}
        
        bool is_divisible() const {return mGrainSize < this->size();} 
 
        LeafRange( LeafRange& r, tbb::split )
            : mEnd(r.mEnd), mBegin(do_split(r)), mGrainSize(r.mGrainSize),
              mLeafManager(r.mLeafManager) {} 
    private:
        size_t mEnd, mBegin, mGrainSize;
        const LeafManager& mLeafManager;
 
        static size_t do_split( LeafRange& r ) {
            assert(r.is_divisible());
            size_t middle = r.mBegin + (r.mEnd - r.mBegin)/2u;
            r.mEnd = middle;
            return middle;
        }
    };// end of LeafRange
    
    /// @brief Constructor from a tree reference and auxiliary buffers
    /// count (default is no auxiliary buffers).
    LeafManager(TreeType &tree, size_t auxBuffersPerLeaf=0, bool serial=false):
        mTree(&tree),
        mLeafCount(0),
        mAuxBufferCount(0),
        mAuxBuffersPerLeaf(auxBuffersPerLeaf),
        mLeafs(NULL),
        mAuxBuffers(NULL),
        mTask(0),
        mIsMaster(true)
      {
        this->rebuild(serial);
      }

    /// Shallow copy constructor called by tbb::parallel_for() threads
    ///
    /// @note This should never get called directly
    LeafManager(const LeafManager& other):
        mTree(other.mTree),
        mLeafCount(other.mLeafCount),
        mAuxBufferCount(other.mAuxBufferCount),
        mAuxBuffersPerLeaf(other.mAuxBuffersPerLeaf),
        mLeafs(other.mLeafs),
        mAuxBuffers(other.mAuxBuffers),
        mTask(other.mTask),
        mIsMaster(false)
      {
      }

    virtual ~LeafManager()
      {
          if (mIsMaster) {
              delete [] mLeafs;
              delete [] mAuxBuffers;
          }
      }

    /// @brief Use this method for (re)-initialization. It clears all
    /// existing leafs and auxiliary buffers and allocates new
    /// ones. New auxiliary buffers are initialized to the
    /// corresponding leaf node buffer.
    void rebuild(bool serial=false)
    {
        this->initLeafs();
        this->initAuxBuffers(serial);
    }
    /// @brief Clears all existing leafs and auxiliary buffers and
    /// allocates new ones.
    void rebuild(size_t auxBuffersPerLeaf, bool serial=false)
      {
          mAuxBuffersPerLeaf = auxBuffersPerLeaf;
          this->rebuild(serial);
      }
    /// @brief Clears all existing leafs and auxiliary buffers and
    /// allocates new ones.
    void rebuild(TreeType& tree, bool serial=false)
      {
          mTree = &tree;
          this->rebuild(serial);
      }
    /// @brief Clears all existing leafs and auxiliary buffers and
    /// allocates new ones.
    void rebuild(TreeType& tree, size_t auxBuffersPerLeaf, bool serial=false)
      {
          mTree = &tree;
          mAuxBuffersPerLeaf = auxBuffersPerLeaf;
          this->rebuild(serial);
      }
    /// Use this method to change to number of auxiliary buffers. If
    /// auxBuffersPerLeaf=0 all existing auxiliary buffers are deleted,
    /// but the leaf nodes are always unchanged. Also note that new
    /// auxiliary buffers are initialized as a copy of the
    /// corresponding leaf node buffer.
    void rebuildAuxBuffers(size_t auxBuffersPerLeaf, bool serial=false)
      {
          mAuxBuffersPerLeaf = auxBuffersPerLeaf;
          this->initAuxBuffers(serial);
      }
    /// @brief Remove the auxiliary buffers (the leafs are unchanged)
    void removeAuxBuffers() { this->rebuildAuxBuffers(0); }
    
    /// @brief Remove the auxiliary buffers and rebuild the array of leafs
    void rebuildLeafs()
    {
        this->removeAuxBuffers();
        this->initLeafs();
    }
    
    /// The total number of allocated auxiliary buffers
    size_t auxBufferCount()  const { return mAuxBufferCount; }
    /// The number of auxiliary buffers per LeafNode
    size_t auxBuffersPerLeaf() const { return mAuxBuffersPerLeaf; }

    /// @return size of the of the LeafManager, corresponding to the
    /// number of leaf nodes.
    size_t leafCount() const { return mLeafCount; }

    TreeType& tree() { return *mTree; }
    
    /// @return the leaf specified by leafId.
    ///
    /// @note For performance reasons no range check is performed - other than assert!
    LeafType& leaf(size_t leafId) const { assert(leafId<mLeafCount); return *mLeafs[leafId]; }

    /// @return The leaf or auxiliary buffer specified by leafId. If
    /// bufferId is zero the LeafNode buffer is returned, else the
    /// n'th auxiliary buffer where n = bufferId -1.
    ///
    /// @note For performance reasons no range checks are performed on
    /// the inputs - other than asserts! Since auxiliary buffers,
    /// unlike leaf buffers, are likely to not exist, be especially
    /// caseful when specifying the bufferId index.
    BufferType& getBuffer(size_t leafId, size_t bufferId) const
      {
          assert(leafId < mLeafCount);
          assert(bufferId==0 || bufferId -1 < mAuxBuffersPerLeaf);
          return bufferId == 0 ? mLeafs[leafId]->buffer() :
                 mAuxBuffers[ leafId*mAuxBuffersPerLeaf + bufferId -1 ];
      }

    /// Return a tbb::blocked_range of indices to the leafs. Useful for tbb multithreading.
    ///
    /// @note Consider using leafRange() instead which provides access methods to leafs and buffers
    RangeType getRange(size_t grainsize=1)
      {
          return RangeType(0, mLeafCount, grainsize);
      }
    
    /// Return a LeafRange that can conveniently be use for multithreading with tbb.
    LeafRange leafRange(size_t grainsize=1)
      {
          return LeafRange(0, mLeafCount, *this, grainsize);
      }

    /// @brief Swaps the leaf node with the specified auxiliary buffer.
    /// @return true if swap was successfuly
    /// @param bufferId Index of the buffer that will be swapped with
    ///                 the corresponding leaf node buffer. 
    /// @param serial  if false, swap buffers in parallel using multiple threads.
    /// @note Recall that the indexing of auxiliary buffers is 1-based,
    /// since bufferId=0 correponds to the leaf node buffer. So
    /// bufferId=1 corresponds to the first auxiliary buffer.
    bool swapLeafBuffer(size_t bufferId, bool serial = false)
      {
          if (bufferId == 0 || bufferId > mAuxBuffersPerLeaf) return false;
          mTask = boost::bind(&LeafManager::doSwapLeafBuffer, _1, _2, bufferId-1);
          this->cook(serial ? 0 : 512);
          return true;//success
      }
    bool swapBuffer(size_t bufferId1, size_t bufferId2, bool serial = false)
      {
          const size_t b1 = std::min(bufferId1,bufferId2);
          const size_t b2 = std::max(bufferId1,bufferId2);
          if (b1==b2 || b2 > mAuxBuffersPerLeaf) return false;
          if (b1==0) {
              mTask = boost::bind(&LeafManager::doSwapLeafBuffer, _1, _2, b2-1);
          } else {
              mTask = boost::bind(&LeafManager::doSwapAuxBuffer, _1, _2, b1-1, b2-1);
          }
          this->cook(serial ? 0 : 512);
          return true;//success
      }

    /// @brief Syncs up the specified auxiliary buffer with the corresponding leaf node buffer.
    /// @return true if sync was successfuly
    /// @param bufferId Index of the buffer that will contain a
    ///                 copy of the corresponding leaf node buffer. 
    /// @param serial  if false, sync buffers in parallel using multiple threads.
    /// @note Recall that the indexing of auxiliary buffers is 1-based,
    /// since bufferId=0 correponds to the leaf node buffer. So
    /// bufferId=1 corresponds to the first auxiliary buffer.
    bool syncAuxBuffer(size_t bufferId, bool serial = false)
      {
          if (bufferId==0 || bufferId > mAuxBuffersPerLeaf) return false;
          mTask = boost::bind(&LeafManager::doSyncAuxBuffer, _1, _2, bufferId - 1);
          this->cook(serial ? 0 : 64);
          return true;//success
      }
    
    /// @brief Syncs up all auxiliary buffer with the corresponding leaf node buffer.
    /// @return true if sync was successfuly 
    /// @param serial  if false, sync buffers in parallel using multiple threads.
    bool syncAllBuffers(bool serial = false)
      {
          switch (mAuxBuffersPerLeaf) {
            case 0:
                return false;//nothing to do
            case 1:
                mTask = boost::bind(&LeafManager::doSyncAllBuffers1, _1, _2); break;
            case 2:
                mTask = boost::bind(&LeafManager::doSyncAllBuffers2, _1, _2); break;
            default:
                mTask = boost::bind(&LeafManager::doSyncAllBuffersN, _1, _2);
          }
          this->cook(serial ? 0 : 64);
          return true;//success
      }

    ////////////////////////////////////////////////////////////////////////////////////
    /// All methods below are for internal use only and should never be called directly

    /// Used internally by tbb::parallel_for() - never call it directly!
    void operator()(const RangeType& r) const
      {
          if (mTask) mTask(const_cast<LeafManager*>(this), r);
          else OPENVDB_THROW(ValueError, "task is undefined");
      }

  private: 
   
    void initLeafs()
      {
          const size_t leafCount = mTree->leafCount();
          if (leafCount != mLeafCount) {
              delete [] mLeafs;
              mLeafs = (leafCount == 0) ? NULL : new LeafType*[leafCount];
              mLeafCount = leafCount;
          }
          typename TreeType::LeafIter iter = mTree->beginLeaf();
          for (size_t n = 0; n != leafCount; ++n, ++iter) mLeafs[n] = iter.getLeaf();
      }
    void initAuxBuffers(bool serial)
      {
          const size_t auxBufferCount = mLeafCount * mAuxBuffersPerLeaf;
          if (auxBufferCount != mAuxBufferCount) {
              delete [] mAuxBuffers;
              mAuxBuffers = (auxBufferCount == 0) ? NULL : new BufferType[auxBufferCount];
              mAuxBufferCount = auxBufferCount;
          }
          this->syncAllBuffers(serial);
      }
    void cook(size_t grainsize)
      {
          if (grainsize>0) {
              tbb::parallel_for(this->getRange(grainsize), *this);
          } else {
              (*this)(this->getRange());
          }
      }
    void doSwapLeafBuffer(const RangeType& r, size_t auxBufferId)
      {
          for (size_t n = r.begin(), m = r.end(), N = mAuxBuffersPerLeaf; n != m; ++n) {
              mLeafs[n]->swap(mAuxBuffers[ n*N + auxBufferId ]);
          }
      }
    void doSwapAuxBuffer(const RangeType& r, size_t auxBufferId1, size_t auxBufferId2)
       {
           for (size_t N = mAuxBuffersPerLeaf, n = N*r.begin(), m = N*r.end(); n != m; n+=N) {
               mAuxBuffers[n + auxBufferId1].swap(mAuxBuffers[n + auxBufferId2]);
           }
      }
    void doSyncAuxBuffer(const RangeType& r, size_t auxBufferId)
      {
          for (size_t n = r.begin(), m = r.end(), N = mAuxBuffersPerLeaf; n != m; ++n) {
              mAuxBuffers[ n*N + auxBufferId ] = mLeafs[n]->buffer();
          }
      }
    void doSyncAllBuffers1(const RangeType& r)
      {
          for (size_t n = r.begin(), m = r.end(); n != m; ++n) mAuxBuffers[n] = mLeafs[n]->buffer();
      }
     void doSyncAllBuffers2(const RangeType& r)
      {
          for (size_t n = r.begin(), m = r.end(); n != m; ++n) {
              const BufferType& leafBuffer = mLeafs[n]->buffer();
              mAuxBuffers[2*n  ] = leafBuffer;
              mAuxBuffers[2*n+1] = leafBuffer;
          }
      }
     void doSyncAllBuffersN(const RangeType& r)
      {
          for (size_t n = r.begin(), m = r.end(), N = mAuxBuffersPerLeaf; n != m; ++n) {
              const BufferType& leafBuffer = mLeafs[n]->buffer();
              for (size_t i=n*N, j=i+N; i!=j; ++i) mAuxBuffers[i] = leafBuffer;
          }
      }
    typedef typename boost::function<void (LeafManager*, const RangeType&)> FuncType;
    
    TreeType*   mTree;
    size_t      mLeafCount, mAuxBufferCount, mAuxBuffersPerLeaf;
    LeafType**  mLeafs;//array of LeafNode pointers
    BufferType* mAuxBuffers;//array of auxiliary buffers
    FuncType    mTask;
    const bool  mIsMaster;
};//end of LeafManager class

} // namespace tree
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TREE_LEAF_MANAGER_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
