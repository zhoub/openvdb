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

#ifndef OPENVDB_TOOLS_MORPHOLOGY_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_MORPHOLOGY_HAS_BEEN_INCLUDED

#include <openvdb/Types.h>
#include <openvdb/tree/TreeIterator.h>
#include <openvdb/tree/ValueAccessor.h>
#include <openvdb/tree/LeafManager.h>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tools {

/// Topologically dilate all leaf-level active voxels in the given tree,
/// i.e., expand the set of active voxels by one voxel in the +x, -x,
/// +y, -y, +z and -z directions, but don't change the values of any voxels,
/// only their active states.
/// @todo Currently operates only on leaf voxels; need to extend to tiles.
template<typename TreeType> OPENVDB_STATIC_SPECIALIZATION
inline void dilateVoxels(TreeType& tree);

/// Topologically erode all leaf-level active voxels in the given tree,
/// i.e., shrink the set of active voxels by one voxel in the +x, -x,
/// +y, -y, +z and -z directions, but don't change the values of any voxels,
/// only their active states.
/// @todo template<typename TreeType> inline void erodeVoxels(TreeType& tree);


////////////////////////////////////////


/// Mapping from a Log2Dim to a data type of size 2^Log2Dim bits
template<Index Log2Dim> struct DimToWord { typedef uint8_t Type[0]; };
template<> struct DimToWord<3> { typedef uint8_t Type; };
template<> struct DimToWord<4> { typedef uint16_t Type; };
template<> struct DimToWord<5> { typedef uint32_t Type; };
template<> struct DimToWord<6> { typedef uint64_t Type; };


////////////////////////////////////////


#undef OPENVDB_TOOLS_MORPHOLOGY_USE_WORKAROUND_FOR_BROKEN_ICC12_COMPILER
#ifdef __ICC
#if __ICC >= 1210 && __ICC < 1300
// ICC 12.1 generates incorrect optimized code for the following dilateVoxels() implementation.
// A less-efficient alternate implementation is provided instead for use with that compiler.
#define OPENVDB_TOOLS_MORPHOLOGY_USE_WORKAROUND_FOR_BROKEN_ICC12_COMPILER
#endif
#endif


#ifndef OPENVDB_TOOLS_MORPHOLOGY_USE_WORKAROUND_FOR_BROKEN_ICC12_COMPILER

template<typename TreeType>
OPENVDB_STATIC_SPECIALIZATION inline void
dilateVoxels(tree::LeafManager<TreeType>& manager)
{
    typedef typename TreeType::LeafNodeType LeafType;
    typedef typename LeafType::NodeMaskType MaskType;
    typedef tree::ValueAccessor<TreeType>   Accessor;

    static const Index LEAF_DIM     = LeafType::DIM;
    static const Index LEAF_LOG2DIM = LeafType::LOG2DIM;

    typedef typename DimToWord<LEAF_LOG2DIM>::Type Word;

    Accessor acc(manager.tree());

    const Index leafCount = manager.leafCount();

    // Save the value masks of all leaf nodes.
    std::vector<MaskType> savedMasks;
    savedMasks.resize(leafCount);
    for (Index i = 0; i < leafCount; ++i) {
        savedMasks[i] = manager.leaf(i).getValueMask();
    }

    struct Neighbor {
        LeafType* leaf;//null if a tile
        bool      isOn;//true if active tile
        Coord     orig;
        Neighbor(Accessor& acc, const Coord& xyz):
            leaf(acc.probeLeaf(xyz)), isOn(false), orig(xyz)
        {
            if (leaf == NULL) isOn = acc.isValueOn(xyz);
        }
        void createLeaf(Accessor& acc) { if (leaf == NULL) leaf = acc.touchLeaf(orig); }
    };

    Coord origin;
    for (Index i = 0; i < leafCount; ++i) {

        const MaskType& oldMask = savedMasks[i];
        LeafType& leaf = manager.leaf(i);//current leaf node
        leaf.getOrigin(origin);

        // Cache the six neighbor nodes of the current leaf node. Note
        // they can be a leaf node or an active or inactive tile node.
        Neighbor MX(acc, origin.offsetBy(-LEAF_DIM,         0,        0));
        Neighbor PX(acc, origin.offsetBy( LEAF_DIM,         0,        0));
        Neighbor MY(acc, origin.offsetBy(        0, -LEAF_DIM,        0));
        Neighbor PY(acc, origin.offsetBy(        0,  LEAF_DIM,        0));
        Neighbor MZ(acc, origin.offsetBy(        0,         0,-LEAF_DIM));
        Neighbor PZ(acc, origin.offsetBy(        0,         0, LEAF_DIM));

        for (Index x = 0; x < LEAF_DIM; ++x) {
            for (Index y = 0, n = (x << LEAF_LOG2DIM); y < LEAF_DIM; ++y, ++n) {
                // Extract the portion of the original mask that corresponds to a row in z.
                const Word oldWord = oldMask.template getWord<Word>(n);
                if (oldWord == 0) continue; // no active voxels

                // Dilate the current leaf node in the z direction by ORing its mask
                // with itself shifted first left and then right by one bit.
                leaf.getValueMask().template getWord<Word>(n) |= (oldWord >> 1) | (oldWord << 1);

                if (!MZ.isOn && Word(oldWord<<(LEAF_DIM-1))) {
                    MZ.createLeaf(acc);
                    MZ.leaf->getValueMask().template getWord<Word>(n)
                        |= Word(oldWord<<(LEAF_DIM-1));
                }
                if (!PZ.isOn && Word(oldWord>>(LEAF_DIM-1))) {
                    PZ.createLeaf(acc);
                    PZ.leaf->getValueMask().template getWord<Word>(n)
                        |= Word(oldWord>>(LEAF_DIM-1));
                }

                if (x > 0) {
                    // Dilate row (x-1, y, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    leaf.getValueMask().template getWord<Word>(n-LEAF_DIM) |= oldWord;
                } else if (!MX.isOn) {
                    // If the neighbor in the -x direction is an active tile, there's no need
                    // to propagate set bits.  Otherwise, dilate into the neighbor's first row.
                    MX.createLeaf(acc);
                    // Dilate into the neighbor's first row by ORing its mask with that of
                    // the current leaf node's original row (0, y, z).
                    MX.leaf->getValueMask().template getWord<Word>(n+LEAF_DIM*(LEAF_DIM-1))
                        |= oldWord;
                }

                if (x < LEAF_DIM - 1) {
                    // Dilate row (x+1, y, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    leaf.getValueMask().template getWord<Word>(n+LEAF_DIM) |= oldWord;
                } else if (!PX.isOn) {
                    // Dilate into the first row of the neighbor in the +x direction.
                    PX.createLeaf(acc);
                    PX.leaf->getValueMask().template getWord<Word>(n-LEAF_DIM*(LEAF_DIM-1))
                        |= oldWord;
                }

                if (y > 0) {
                    // Dilate row (x, y-1, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    leaf.getValueMask().template getWord<Word>(n-1) |= oldWord;
                } else if (!MY.isOn) {
                    // Dilate into the first row of the neighbor in the -y direction.
                    MY.createLeaf(acc);
                    MY.leaf->getValueMask().template getWord<Word>(n+LEAF_DIM-1) |= oldWord;
                }

                if (y < LEAF_DIM - 1) {
                    // Dilate row (x, y+1, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    leaf.getValueMask().template getWord<Word>(n+1) |= oldWord;
                } else if (!PY.isOn) {
                    // Dilate into the first row of the neighbor in the +y direction.
                    PY.createLeaf(acc);
                    PY.leaf->getValueMask().template getWord<Word>(n-LEAF_DIM+1) |= oldWord;
                }
            }
        }
    }
}

#else // ifdef OPENVDB_TOOLS_MORPHOLOGY_USE_WORKAROUND_FOR_BROKEN_ICC12_COMPILER

//#warning using alternate tools::dilateVoxels() implementation due to broken ICC 12.1 compiler

// ICC 12.1 generates incorrect optimized code for the above implementation.
// This less-efficient implementation is used instead with that compiler.

/// Helper class to find or create the neighbors of a leaf node
template<typename TreeType>
struct LeafLookup
{
    typedef typename TreeType::LeafNodeType LeafType;

    tree::ValueAccessor<TreeType> cache;
    // These dummy leaf nodes can be inserted into the cache, where they
    // represent neighbors that are either active or inactive tiles.
    const boost::shared_ptr<LeafType> onTile, offTile;

    LeafLookup(TreeType& tree):
        cache(tree), onTile(new LeafType), offTile(new LeafType) {}

    bool isActiveTile(const LeafType* leaf) const { return leaf == onTile.get(); }
    bool isInactiveTile(const LeafType* leaf) const { return leaf == offTile.get(); }

    // Return the leaf node that contains voxel (x, y, z), or, if the voxel lies
    // inside a tile, return either the onTile or offTile dummy leaf node.
    LeafType* operator()(const Coord& xyz)
    {
        // Remove any existing leaf node from the cache.
        cache.template eraseNode<LeafType>();
        // Trigger a cache update, which may or may not replace the preloaded node.
        bool on = cache.isValueOn(xyz);
        LeafType* leaf = cache.template getNode<LeafType>();
        // If the "active tile" node is still in the cache, but voxel (x, y, z)
        // is actually inactive, then it must lie inside an inactive tile.
        return leaf != NULL ? leaf : (on ? onTile.get() : offTile.get());
    }

    // Create and return the leaf node that contains voxel (x, y, z).
    LeafType* createLeaf(const Coord& xyz)
    {
        cache.template eraseNode<LeafType>();
        // Mark the voxel as active to force creation of the node.
        cache.setValueOn(xyz);
        LeafType* leaf = cache.template getNode<LeafType>();
        // Mark all voxels as inactive.
        leaf->setValuesOff();
        return leaf;
    }
}; // struct LeafLookup


template<typename TreeType>
inline void
dilateVoxels(tree::LeafManager<TreeType>& manager)
{
    typedef typename TreeType::LeafNodeType LeafType;
    typedef typename LeafType::NodeMaskType MaskType;

    static const Index LEAF_DIM     = LeafType::DIM;
    static const Index LEAF_LOG2DIM = LeafType::LOG2DIM;

    typedef typename DimToWord<LEAF_LOG2DIM>::Type Word;

    /// @todo Currently operates only on leaf voxels; need to extend to tiles.
    const Index leafCount = manager.leafCount();

    // Save the value masks of all leaf nodes.
    std::vector<MaskType> savedMasks;
    savedMasks.resize(leafCount);
    for (Index i = 0; i < leafCount; ++i) {
        savedMasks[i] = manager.leaf(i).getValueMask();
    }

    LeafLookup<TreeType> leafLookup(manager.tree());

    Coord origin;
    for (Index i = 0; i < leafCount; ++i) {
        const MaskType& oldMask = savedMasks[i];

        // Get the origin of the current leaf node.
        manager.leaf(i).getOrigin(origin);

        // This array caches pointers to the current leaf node (nbhd[0]) and four of its
        // neighbors, some of which might not yet exist and might need to be created
        LeafType* nbhd[5] = {
            &manager.leaf(i),
            leafLookup(origin.offsetBy(-LEAF_DIM,         0, 0)),
            leafLookup(origin.offsetBy( LEAF_DIM,         0, 0)),
            leafLookup(origin.offsetBy(        0, -LEAF_DIM, 0)),
            leafLookup(origin.offsetBy(        0,  LEAF_DIM, 0))
        };

        for (Index x = 0; x < LEAF_DIM; ++x) {
            for (Index y = 0, n = (x << LEAF_LOG2DIM); y < LEAF_DIM; ++y, ++n) {
                // Extract the portion of the original mask that corresponds to a row in z.
                const Word oldWord = oldMask.template getWord<Word>(n);
                if (oldWord == 0) continue; // no active voxels

                // Dilate the current leaf node in the z direction by ORing its mask
                // with itself shifted first left and then right by one bit.
                nbhd[0]->getValueMask().template getWord<Word>(n) |=
                    (oldWord >> 1) | (oldWord << 1);

                if (oldWord & 1) {
                    // If the low bit of the original mask was set, it needs to propagate
                    // into the mask of the neighboring leaf node in the -z direction.
                    leafLookup.cache.setValueOn(origin.offsetBy(x, y, -1));
                }
                if (Word(oldWord >> (LEAF_DIM - 1))) {
                    // If the high bit of the original mask was set, it needs to propagate
                    // into the mask of the neighboring leaf node in the +z direction.
                    leafLookup.cache.setValueOn(origin.offsetBy(x, y, LEAF_DIM));
                }

                if (x > 0) {
                    // Dilate row (x-1, y, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    nbhd[0]->getValueMask().template getWord<Word>(n - LEAF_DIM) |= oldWord;
                } else if (!leafLookup.isActiveTile(nbhd[1])) {
                    // If the neighbor in the -x direction is an active tile, there's no need
                    // to propagate set bits.  Otherwise, dilate into the neighbor's first row.
                    if (leafLookup.isInactiveTile(nbhd[1])) {
                        // If the neighbor is an inactive tile, create a neighbor leaf node.
                        nbhd[1] = leafLookup.createLeaf(origin.offsetBy(-LEAF_DIM, 0, 0));
                    }
                    // Dilate into the neighbor's first row by ORing its mask with that of
                    // the current leaf node's original row (0, y, z).
                    nbhd[1]->getValueMask().template getWord<Word>(
                        n + LEAF_DIM * (LEAF_DIM - 1)) |= oldWord;
                }

                if (x < LEAF_DIM - 1) {
                    // Dilate row (x+1, y, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    nbhd[0]->getValueMask().template getWord<Word>(n + LEAF_DIM) |= oldWord;
                } else if (!leafLookup.isActiveTile(nbhd[2])) {
                    // Dilate into the first row of the neighbor in the +x direction.
                    if (leafLookup.isInactiveTile(nbhd[2])) {
                        // If the neighbor is an inactive tile, create a neighbor leaf node.
                        nbhd[2] = leafLookup.createLeaf(origin.offsetBy(LEAF_DIM, 0, 0));
                    }
                    nbhd[2]->getValueMask().template getWord<Word>(
                        n - LEAF_DIM * (LEAF_DIM - 1)) |= oldWord;
                }

                if (y > 0) {
                    // Dilate row (x, y-1, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    nbhd[0]->getValueMask().template getWord<Word>(n - 1) |= oldWord;
                } else if (!leafLookup.isActiveTile(nbhd[3])) {
                    // Dilate into the first row of the neighbor in the -y direction.
                    if (leafLookup.isInactiveTile(nbhd[3])) {
                        // If the neighbor is an inactive tile, create a neighbor leaf node.
                        nbhd[3] = leafLookup.createLeaf(origin.offsetBy(0, -LEAF_DIM, 0));
                    }
                    nbhd[3]->getValueMask().template getWord<Word>(n+LEAF_DIM-1) |= oldWord;
                }

                if (y < LEAF_DIM - 1) {
                    // Dilate row (x, y+1, z) of the current leaf node by ORing its mask
                    // with that of the original row (x, y, z).
                    nbhd[0]->getValueMask().template getWord<Word>(n + 1) |= oldWord;
                } else if (!leafLookup.isActiveTile(nbhd[4])) {
                    // Dilate into the first row of the neighbor in the +y direction.
                    if (leafLookup.isInactiveTile(nbhd[4])) {
                        // If the neighbor is an inactive tile, create a neighbor leaf node.
                        nbhd[4] = leafLookup.createLeaf(origin.offsetBy(0, LEAF_DIM, 0));
                    }
                    nbhd[4]->getValueMask().template getWord<Word>(n-LEAF_DIM+1) |= oldWord;
                }
            }
        }
    }
}

#endif // defined(OPENVDB_TOOLS_MORPHOLOGY_USE_WORKAROUND_FOR_BROKEN_ICC12_COMPILER)


template<typename TreeType>
OPENVDB_STATIC_SPECIALIZATION inline void
dilateVoxels(TreeType& tree)
{
    tree::LeafManager<TreeType> manager(tree);
    dilateVoxels<TreeType>(manager);
}

} // namespace tools
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_MORPHOLOGY_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
