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

#include <vector>
#include <cppunit/extensions/HelperMacros.h>

#include <openvdb/openvdb.h>
#include <openvdb/Exceptions.h>

#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/util/Util.h>


class TestVolumeToMesh: public CppUnit::TestCase
{
public:
    CPPUNIT_TEST_SUITE(TestVolumeToMesh);
    CPPUNIT_TEST(testAuxData);
    CPPUNIT_TEST_SUITE_END();

    void testAuxData();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestVolumeToMesh);


////////////////////////////////////////


void
TestVolumeToMesh::testAuxData()
{
    openvdb::FloatTree::Ptr tree(new openvdb::FloatTree(0));

    // create one voxel with 3 upwind edges (that have a sign change)
    tree->setValue(openvdb::Coord(0,0,0), -1);
    tree->setValue(openvdb::Coord(1,0,0),  1);
    tree->setValue(openvdb::Coord(0,1,0),  1);
    tree->setValue(openvdb::Coord(0,0,1),  1);


    openvdb::tools::internal::LeafCPtrList<openvdb::FloatTree> leafs(*tree);

    openvdb::tools::internal::AuxiliaryData<openvdb::FloatTree> init(*tree, leafs, 0.0);

    init.runParallel();

    CPPUNIT_ASSERT(init.auxTree()->activeVoxelCount() == 7);
    CPPUNIT_ASSERT(init.edgeTree()->activeVoxelCount() == 1);

    int flags = int(init.edgeTree()->getValue(openvdb::Coord(0,0,0)));

    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::INSIDE));
    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::XEDGE));
    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::YEDGE));
    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::ZEDGE));


    tree->setValueOff(openvdb::Coord(0,0,1), -1);

    init.runParallel();

    CPPUNIT_ASSERT(init.auxTree()->activeVoxelCount() == 7);
    CPPUNIT_ASSERT(init.edgeTree()->activeVoxelCount() == 1);

    flags = int(init.edgeTree()->getValue(openvdb::Coord(0,0,0)));

    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::INSIDE));
    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::XEDGE));
    CPPUNIT_ASSERT(bool(flags & openvdb::tools::internal::YEDGE));
    CPPUNIT_ASSERT(!bool(flags & openvdb::tools::internal::ZEDGE));
}

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
