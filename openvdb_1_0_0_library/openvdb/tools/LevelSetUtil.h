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

#ifndef OPENVDB_TOOLS_LEVELSETUTIL_HAS_BEEN_INCLUDED
#define OPENVDB_TOOLS_LEVELSETUTIL_HAS_BEEN_INCLUDED

#include <openvdb/Grid.h>
#include <openvdb/math/Quat.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tree/LeafManager.h>
#include <openvdb/util/NullInterrupter.h>
#include <openvdb/util/Util.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <limits>
#include <list>
#include <deque> // used by segment()


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace tools {

// MS Visual C++ requires this extra level of indirection in order to compile
// THIS MUST EXIST IN AN UNNAMED NAMESPACE IN ORDER TO COMPILE ON WINDOWS
namespace {

template<typename GridType>
inline typename GridType::ValueType lsutilGridMax()
{
    return std::numeric_limits<typename GridType::ValueType>::max();
}

template<typename GridType>
inline typename GridType::ValueType lsutilGridZero()
{
    return zeroVal<typename GridType::ValueType>();
}

} // unnamed namespace


/// @brief  Threaded method to convert a sparse level-set/SDF into a sparse fog volume.
/// @note   The active negative/interior values of the narrow-band are transformed into
/// a smooth [0 to 1] gradient, and the interior inactive values are transformed to active
/// constanct values of 1. All exterior values, including the background, are set to 0 with
/// inactive states. The interior is still a sparse representations but the values arenow
/// active.
///
/// @param grid            Level set / SDF grid to transform.
/// @param cutOffDistance  Optional world space cutoff distance for the gradient
///                        (automatically clamped if greater than the interior
///                        narrow band width)
template<class GridType>
inline void
sdfToFogVolume(
    GridType& grid,
    typename GridType::ValueType cutOffDistance = lsutilGridMax<GridType>());


//////////


/// @brief Threaded method to extract the interior region mask from a level-set/SDF grid.
///
/// @return Shared ptr to a new boolean grid with the same tree configuration and
///         transform as the incoming @c grid. The interior region is made active.
///
/// @param grid     A Level set / SDF grid.
/// @param iso      Values below this threshold are considered to be part of the
///                 interior region.
///
template<class GridType>
inline typename Grid<typename GridType::TreeType::template ValueConverter<bool>::Type>::Ptr
sdfInteriorMask(
    const GridType& grid,
    typename GridType::ValueType iso = lsutilGridZero<GridType>());


//////////


/// @brief Threaded operation to find the min and max active voxel values.
template<class TreeType>
class MinMaxVoxel
{
public:
    typedef tree::LeafManager<TreeType> LeafArray;
    typedef typename TreeType::ValueType ValueType;

    /// LeafArray = openvdb::tree::LeafManager<TreeType> leafs(myTree)
    MinMaxVoxel(LeafArray&);

    void runParallel();
    void runSerial();

    const ValueType& minVoxel() const { return mMin; }
    const ValueType& maxVoxel() const { return mMax; }


    inline MinMaxVoxel(const MinMaxVoxel<TreeType>&, tbb::split);
    inline void operator()(const tbb::blocked_range<size_t>&);
    inline void join(const MinMaxVoxel<TreeType>&);

private:
    LeafArray& mLeafArray;
    ValueType mMin, mMax;
};


//////////


/// @brief Threaded leaf-node transformation scheme.
template<class TreeType, class LeafOp>
class LeafTransformer
{
public:
    typedef tree::LeafManager<TreeType> LeafArray;
    typedef typename TreeType::ValueType ValueType;

    /// LeafArray = openvdb::tree::LeafManager<TreeType> leafs(myTree)
    LeafTransformer(LeafArray&, LeafOp&);

    void runParallel();
    void runSerial();


    inline void operator()(const tbb::blocked_range<size_t>&) const;
    inline LeafTransformer(const LeafTransformer<TreeType, LeafOp>&);

private:
    LeafArray& mLeafArray;
    LeafOp& mLeafOp;
};


//////////


/// @brief Level set fracturing
template<class GridType, class InterruptType = util::NullInterrupter>
class LevelSetFracture
{
public:
    typedef std::vector<Vec3s> Vec3sList;
    typedef std::vector<math::Quats> QuatsList;
    typedef std::list<typename GridType::Ptr> GridPtrList;
    typedef typename GridPtrList::iterator GridPtrListIter;


    /// @brief Default constructor
    ///
    /// @param interrupter  Optional interrupter object.
    LevelSetFracture(InterruptType* interrupter = NULL);


    /// @brief  Fracture level set grids
    ///
    /// @note   The incomming @a grids and the @a cutter are required to have matching
    ///         transforms and narrow band widths!
    ///
    /// @param  grids           List of grids to fracture. The residuals of the
    ///                         fractured grids will remain in this list.
    /// @param  cutter          A level set grid to use as the cutter object.
    /// @param  segment         Toggle to split disjoint fragments into their own girds.
    /// @param  points          Optional list of world space points to instance the cutter
    ///                         object onto. The cutter object is used in place if no points
    ///                         are provided.
    /// @param  rotations       Optional list of custom rotations for each cutter instance.
    /// @param  cutterOverlap   Toggle to allow consecutive cutter instances to fracture
    ///                         previously generated fragments.
    void fracture(GridPtrList& grids, const GridType& cutter, bool segment = false,
        const Vec3sList* points = NULL, const QuatsList* rotations = NULL, bool cutterOverlap = true);


    /// @return List of new fragments (does not include the residuals from the input grids).
    GridPtrList& fragments() { return mFragments; }


    /// @brief Removes all elements from the fragment list.
    void clear() { mFragments.clear(); }


private:
    // disallow copy by assignment
    void operator=(const LevelSetFracture<GridType, InterruptType>&) {}

    bool wasInterrupted(int percent = -1) const {
        return mInterrupter && mInterrupter->wasInterrupted(percent);
    }

    bool isValidFragment(GridType&) const;
    void segmentFragments(GridPtrList&) const;
    void process(GridPtrList&, const GridType& cutter);

    InterruptType* mInterrupter;
    GridPtrList mFragments;
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Internal utility objects and implementation details


namespace internal {


template<typename ValueType>
struct FogVolumeOp
{
    FogVolumeOp(ValueType cutOffDistance)
        : mWeight(ValueType(1.0) / cutOffDistance)
    {
    }

    // cutOff has to be < 0.0
    template <typename LeafNodeType>
    void operator()(LeafNodeType &leaf, size_t/*leafIndex*/) const
    {
        const ValueType zero = zeroVal<ValueType>();

        for (typename LeafNodeType::ValueAllIter iter = leaf.beginValueAll(); iter; ++iter) {
            ValueType& value = const_cast<ValueType&>(iter.getValue());
            if (value > zero) {
                value = zero;
                iter.setValueOff();
            } else {
                value = std::min(ValueType(1.0), value * mWeight);
                iter.setValueOn(value > zero);
            }
        }
    }

private:
    ValueType mWeight;
};


template<typename TreeType>
struct InteriorMaskOp
{
    InteriorMaskOp(const TreeType& tree, typename TreeType::ValueType iso)
        : mTree(tree)
        , mIso(iso)
    {
    }

    template <typename LeafNodeType>
    void operator()(LeafNodeType &leaf, size_t/*leafIndex*/) const
    {
        const Coord origin = leaf.getOrigin();
        const typename TreeType::LeafNodeType* refLeafPt = mTree.probeConstLeaf(origin);

        if (refLeafPt != NULL) {

            const typename TreeType::LeafNodeType& refLeaf = *refLeafPt;
            typename LeafNodeType::ValueAllIter iter = leaf.beginValueAll();

            for (; iter; ++iter) {
                if (refLeaf.getValue(iter.pos()) < mIso) {
                    iter.setValueOn();
                } else {
                    iter.setValueOff();
                }
            }
        }
    }

private:
    const TreeType& mTree;
    typename TreeType::ValueType mIso;
};


/// @brief Segmentation scheme, splits disjoint fragments into separate grids.
/// @note This is a temporary solution and it will be replaced soon.
template<typename GridType, typename InterruptType>
inline
std::vector<typename GridType::Ptr> segment(GridType& grid, InterruptType* interrupter = NULL)
{
    typedef typename GridType::Ptr GridPtr;
    typedef typename GridType::TreeType TreeType;
    typedef typename TreeType::Ptr TreePtr;
    typedef typename TreeType::ValueType ValueType;

    std::vector<GridPtr> segments;

    while (grid.activeVoxelCount() > 0) {

        if (interrupter && interrupter->wasInterrupted()) break;

        // Deep copy the grid's metadata (tree and transform are shared)
        GridPtr segment(new GridType(grid, ShallowCopy()));
        // Make the transform unique and insert an empty tree
        segment->setTransform(grid.transform().copy());
        TreePtr tree(new TreeType(grid.background()));
        segment->setTree(tree);

        std::deque<Coord> coordList;
        coordList.push_back(grid.tree().beginLeaf()->beginValueOn().getCoord());

        Coord ijk, n_ijk;
        ValueType value;

        typename tree::ValueAccessor<TreeType> sourceAcc(grid.tree());
        typename tree::ValueAccessor<TreeType> targetAcc(segment->tree());

        while (!coordList.empty()) {

            if (interrupter && interrupter->wasInterrupted()) break;

            ijk = coordList.back();
            coordList.pop_back();

            if (!sourceAcc.probeValue(ijk, value)) continue;
            if (targetAcc.isValueOn(ijk)) continue;

            targetAcc.setValue(ijk, value);
            sourceAcc.setValueOff(ijk);

            for (int n = 0; n < 6; n++) {
                n_ijk = ijk + util::COORD_OFFSETS[n];
                if (!targetAcc.isValueOn(n_ijk) && sourceAcc.isValueOn(n_ijk)) {
                    coordList.push_back(n_ijk);
                }
            }
        }

        grid.tree().pruneInactive();
        segment->tree().signedFloodFill();
        segments.push_back(segment);
    }



    return segments;
}

} // namespace internal


////////////////////////////////////////


template <class TreeType>
MinMaxVoxel<TreeType>::MinMaxVoxel(LeafArray& leafs)
    : mLeafArray(leafs)
    , mMin(std::numeric_limits<ValueType>::max())
    , mMax(-mMin)
{
}

template <class TreeType>
inline
MinMaxVoxel<TreeType>::MinMaxVoxel(
    const MinMaxVoxel<TreeType>& rhs, tbb::split)
    : mLeafArray(rhs.mLeafArray)
    , mMin(rhs.mMin)
    , mMax(rhs.mMax)
{
}

template <class TreeType>
void
MinMaxVoxel<TreeType>::runParallel()
{
    tbb::parallel_reduce(mLeafArray.getRange(), *this);
}


template <class TreeType>
void
MinMaxVoxel<TreeType>::runSerial()
{
    (*this)(mLeafArray.getRange());
}


template <class TreeType>
inline void
MinMaxVoxel<TreeType>::operator()(const tbb::blocked_range<size_t>& range)
{
    typename TreeType::LeafNodeType::ValueOnCIter iter;

    for (size_t n = range.begin(); n < range.end(); ++n) {
        iter = mLeafArray.leaf(n).cbeginValueOn();
        for (; iter; ++iter) {
            const ValueType value = iter.getValue();
            mMin = std::min(mMin, value);
            mMax = std::max(mMax, value);
        }
    }
}

template <class TreeType>
inline void
MinMaxVoxel<TreeType>::join(const MinMaxVoxel<TreeType>& rhs)
{
    mMin = std::min(mMin, rhs.mMin);
    mMax = std::max(mMax, rhs.mMax);
}


////////////////////////////////////////


template <class TreeType, class LeafOp>
LeafTransformer<TreeType, LeafOp>::
    LeafTransformer(LeafArray& leafs, LeafOp& leafOp)
    : mLeafArray(leafs)
    , mLeafOp(leafOp)
{
}

template <class TreeType, class LeafOp>
inline
LeafTransformer<TreeType, LeafOp>::LeafTransformer(
    const LeafTransformer<TreeType, LeafOp>& rhs)
    : mLeafArray(rhs.mLeafArray)
    , mLeafOp(rhs.mLeafOp)
{
}

template <class TreeType, class LeafOp>
void
LeafTransformer<TreeType, LeafOp>::runParallel()
{
    tbb::parallel_for(mLeafArray.getRange(), *this);
}

template <class TreeType, class LeafOp>
void
LeafTransformer<TreeType, LeafOp>::runSerial()
{
    (*this)(mLeafArray.getRange());
}

template <class TreeType, class LeafOp>
inline void
LeafTransformer<TreeType, LeafOp>::operator()(const tbb::blocked_range<size_t>& range) const
{
    for (size_t n = range.begin(); n < range.end(); ++n) {
        mLeafOp(mLeafArray.leaf(n), n);
    }
}


////////////////////////////////////////


template <class GridType>
inline void
sdfToFogVolume(GridType& grid, typename GridType::ValueType cutOffDistance)
{
    typedef typename GridType::TreeType TreeType;
    typedef typename GridType::ValueType ValueType;

    cutOffDistance = -std::abs(cutOffDistance);

    TreeType& tree = const_cast<TreeType&>(grid.tree());

    { // Transform all voxels (parallel, over leaf nodes)
        tree::LeafManager<TreeType> leafs(tree);

        MinMaxVoxel<TreeType> minmax(leafs);
        minmax.runParallel();

        // Clamp to the interior band width.
        if (minmax.minVoxel() > cutOffDistance) {
            cutOffDistance = minmax.minVoxel();
        }

        internal::FogVolumeOp<ValueType> op(cutOffDistance);
        LeafTransformer<TreeType, internal::FogVolumeOp<ValueType> > transform(leafs, op);
        transform.runParallel();
    }

    // Transform all tile values (serial, but the iteration
    // is constrained from descending into leaf nodes)
    const ValueType zero = zeroVal<ValueType>();
    typename TreeType::ValueAllIter iter(tree);
    iter.setMaxDepth(TreeType::ValueAllIter::LEAF_DEPTH - 1);

    for ( ; iter; ++iter) {
        ValueType& value = const_cast<ValueType&>(iter.getValue());

        if (value > zero) {
            value = zero;
            iter.setValueOff();
        } else {
            value = ValueType(1.0);
            iter.setActiveState(true);
        }
    }

    // Update the tree background value.

    typename TreeType::Ptr newTree(new TreeType(/*background=*/zero));
    newTree->merge(tree);
    // This is faster than calling Tree::setBackground, since we only need
    // to update the value that is returned for coordinates that don't fall
    // inside an allocated node. All inactive tiles and voxels have already
    // been updated in the previous step so the Tree::setBackground method
    // will in this case do a redundant traversal of the tree to update the
    // inactive values once more.

    //newTree->pruneInactive();
    grid.setTree(newTree);

    grid.setGridClass(GRID_FOG_VOLUME);
}


////////////////////////////////////////


template <class GridType>
inline typename Grid<typename GridType::TreeType::template ValueConverter<bool>::Type>::Ptr
sdfInteriorMask(const GridType& grid, typename GridType::ValueType iso)
{
    typedef typename GridType::TreeType::template ValueConverter<bool>::Type BoolTreeType;
    typedef Grid<BoolTreeType> BoolGridType;

    typename BoolGridType::Ptr maskGrid(BoolGridType::create(false));
    maskGrid->setTransform(grid.transform().copy());
    BoolTreeType& maskTree = maskGrid->tree();

    maskTree.topologyUnion(grid.tree());

    { // Evaluate voxels (parallel, over leaf nodes)

        tree::LeafManager<BoolTreeType> leafs(maskTree);

        typedef internal::InteriorMaskOp<typename GridType::TreeType> InteriorMaskOp;
        InteriorMaskOp op(grid.tree(), iso);

        LeafTransformer<BoolTreeType, InteriorMaskOp> transform(leafs, op);
        transform.runParallel();
    }

    // Evaluate tile values (serial, but the iteration
    // is constrained from descending into leaf nodes)

    tree::ValueAccessor<const typename GridType::TreeType> acc(grid.tree());
    typename BoolTreeType::ValueAllIter iter(maskTree);
    iter.setMaxDepth(BoolTreeType::ValueAllIter::LEAF_DEPTH - 1);

    for ( ; iter; ++iter) {
        iter.setActiveState(acc.getValue(iter.getCoord()) < iso);
    }

    maskTree.pruneInactive();

    return maskGrid;
}


////////////////////////////////////////


template<class GridType, class InterruptType>
LevelSetFracture<GridType, InterruptType>::LevelSetFracture(InterruptType* interrupter)
    : mInterrupter(interrupter)
    , mFragments()
{
}


template<class GridType, class InterruptType>
void
LevelSetFracture<GridType, InterruptType>::fracture(GridPtrList& grids, const GridType& cutter,
    bool segmentation, const Vec3sList* points, const QuatsList* rotations, bool cutterOverlap)
{
    // We can process all incoming grids with the same cutter instance,
    // this optimization is enabled by the requirement of having matching
    // transforms between all incoming grids and the cutter object. 
    if (points && points->size() != 0) {

        math::Transform::Ptr originalCutterTransform = cutter.transform().copy();
        GridType cutterGrid(cutter, ShallowCopy());

        const bool hasInstanceRotations = points && rotations && points->size() == rotations->size();

        // for each instance point..
        for (size_t p = 0, P = points->size(); p < P; ++p) {

            int percent = int((float(p) / float(P)) * 100.0);
            if (wasInterrupted(percent)) break;

            GridType instCutterGrid;
            instCutterGrid.setTransform(originalCutterTransform->copy());
            math::Transform::Ptr xform = originalCutterTransform->copy();

            if (hasInstanceRotations) {
                const Vec3s& rot = (*rotations)[p].eulerAngles(math::XYZ_ROTATION);
                xform->preRotate(rot[0], openvdb::math::X_AXIS);
                xform->preRotate(rot[1], openvdb::math::Y_AXIS);
                xform->preRotate(rot[2], openvdb::math::Z_AXIS);
                xform->postTranslate((*points)[p]);
            } else {
                xform->postTranslate((*points)[p]);
            }

            cutterGrid.setTransform(xform);

            if (wasInterrupted()) break;

            resampleToMatch<openvdb::tools::BoxSampler>(cutterGrid, instCutterGrid);

            if (cutterOverlap) process(mFragments, instCutterGrid);
            process(grids, instCutterGrid);
        }

    } else {
        // use cutter in place
        process(grids, cutter);
    }

    if (segmentation) {

        segmentFragments(mFragments);
        segmentFragments(grids);
    }
}


template<class GridType, class InterruptType>
bool
LevelSetFracture<GridType, InterruptType>::isValidFragment(GridType& grid) const
{
    typedef typename GridType::TreeType TreeType;
    if (grid.activeVoxelCount() < 27) return false;

    // Check if valid level-set
    {
        tree::LeafManager<TreeType> leafs(grid.tree());
        MinMaxVoxel<TreeType> minmax(leafs);
        minmax.runParallel();

        if ((minmax.minVoxel() < 0) == (minmax.maxVoxel() < 0)) return false;
    }

    return true;
}


template<class GridType, class InterruptType>
void
LevelSetFracture<GridType, InterruptType>::segmentFragments(GridPtrList& grids) const
{
    GridPtrList newFragments;

    for (GridPtrListIter it = grids.begin(); it != grids.end(); ++it) {

        if (wasInterrupted()) break;

        std::vector<typename GridType::Ptr> segments = internal::segment(*(*it), mInterrupter);
        for (size_t n = 0, N = segments.size(); n < N; ++n) {

            if (wasInterrupted()) break;

            if (isValidFragment(*segments[n])) {
                newFragments.push_back(segments[n]);
            }
        }
    }

    grids.swap(newFragments);
}


template<class GridType, class InterruptType>
void
LevelSetFracture<GridType, InterruptType>::process(
    GridPtrList& grids, const GridType& cutter)
{
    typedef typename GridType::Ptr GridPtr;

    GridPtrList newFragments;

    for (GridPtrListIter it = grids.begin(); it != grids.end(); ++it) {

        if (wasInterrupted()) break;

        GridPtr grid = *it;

        // gen new fragment
        GridPtr fragment = grid->deepCopy();
        openvdb::tools::csgIntersection(*fragment, *cutter.deepCopy());

        if (wasInterrupted()) break;

        if (!isValidFragment(*fragment)) continue;

        // update residual
        GridPtr residual = grid->deepCopy();
        openvdb::tools::csgDifference(*residual, *cutter.deepCopy());

        if (wasInterrupted()) break;

        if (!isValidFragment(*residual)) continue;

        newFragments.push_back(fragment);

        grid->tree().clear();
        grid->tree().merge(residual->tree());
    }

    if(!newFragments.empty()) {
        mFragments.splice(mFragments.end(), newFragments);
    }
}

} // namespace tools
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_LEVELSETUTIL_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
