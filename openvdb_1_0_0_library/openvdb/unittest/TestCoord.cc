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

#include <cppunit/extensions/HelperMacros.h>
#include <openvdb/Types.h>
#include <sstream>


class TestCoord: public CppUnit::TestCase
{
public:
    CPPUNIT_TEST_SUITE(TestCoord);
    CPPUNIT_TEST(testCoord);
    CPPUNIT_TEST(testConversion);
    CPPUNIT_TEST(testIO);
    CPPUNIT_TEST_SUITE_END();

    void testCoord();
    void testConversion();
    void testIO();
};

CPPUNIT_TEST_SUITE_REGISTRATION(TestCoord);


void
TestCoord::testCoord()
{
    using openvdb::Coord;

    Coord xyz(-1, 2, 4);
    Coord xyz2 = -xyz;
    CPPUNIT_ASSERT_EQUAL(Coord(1, -2, -4), xyz2);

    xyz2 = -xyz2;
    CPPUNIT_ASSERT_EQUAL(xyz, xyz2);

    xyz.setX(-xyz.x());
    CPPUNIT_ASSERT_EQUAL(Coord(1, 2, 4), xyz);

    xyz2 = xyz >> 1;
    CPPUNIT_ASSERT_EQUAL(Coord(0, 1, 2), xyz2);

    xyz2 |= 1;
    CPPUNIT_ASSERT_EQUAL(Coord(1, 1, 3), xyz2);

    CPPUNIT_ASSERT(xyz2 != xyz);
    CPPUNIT_ASSERT(xyz2 < xyz);
    CPPUNIT_ASSERT(xyz2 <= xyz);

    xyz2 -= xyz2;
    CPPUNIT_ASSERT_EQUAL(Coord(), xyz2);

    xyz2.reset(0, 4, 4);
    xyz2.offset(-1);
    CPPUNIT_ASSERT_EQUAL(Coord(-1, 3, 3), xyz2);

    // xyz = (1, 2, 4), xyz2 = (-1, 3, 3)
    CPPUNIT_ASSERT_EQUAL(Coord(-1, 2, 3), Coord::minComponent(xyz, xyz2));
    CPPUNIT_ASSERT_EQUAL(Coord(1, 3, 4), Coord::maxComponent(xyz, xyz2));
}


void
TestCoord::testConversion()
{
    using openvdb::Coord;

    openvdb::Vec3I iv(1, 2, 4);
    Coord xyz(iv);
    CPPUNIT_ASSERT_EQUAL(Coord(1, 2, 4), xyz);
    CPPUNIT_ASSERT_EQUAL(iv, xyz.asVec3I());
    CPPUNIT_ASSERT_EQUAL(openvdb::Vec3i(1, 2, 4), xyz.asVec3i());

    iv = (xyz + iv) + xyz;
    CPPUNIT_ASSERT_EQUAL(openvdb::Vec3I(3, 6, 12), iv);
    iv = iv - xyz;
    CPPUNIT_ASSERT_EQUAL(openvdb::Vec3I(2, 4, 8), iv);

    openvdb::Vec3s fv = xyz.asVec3s();
    CPPUNIT_ASSERT(openvdb::math::isExactlyEqual(openvdb::Vec3s(1, 2, 4), fv));
}


void
TestCoord::testIO()
{
    using openvdb::Coord;

    Coord xyz(-1, 2, 4), xyz2;

    std::ostringstream os(std::ios_base::binary);
    CPPUNIT_ASSERT_NO_THROW(xyz.write(os));

    std::istringstream is(os.str(), std::ios_base::binary);
    CPPUNIT_ASSERT_NO_THROW(xyz2.read(is));

    CPPUNIT_ASSERT_EQUAL(xyz, xyz2);

    os.str("");
    os << xyz;
    CPPUNIT_ASSERT_EQUAL(std::string("[-1, 2, 4]"), os.str());
}

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
