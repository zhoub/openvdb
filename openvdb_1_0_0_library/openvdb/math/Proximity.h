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

#ifndef OPENVDB_MATH_PROXIMITY_HAS_BEEN_INCLUDED
#define OPENVDB_MATH_PROXIMITY_HAS_BEEN_INCLUDED

#include <openvdb/Types.h>
//#include <openvdb/openvdb.h>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace math {


/// @brief Squared distance of a line segment p(t) = (1-t)*p0 + t*p1 to point. 
/// @return the closest point on the line segment as a function of t
OPENVDB_API double 
sLineSeg3ToPointDistSqr(const Vec3d &p0, 
                        const Vec3d &p1, 
                        const Vec3d &point, 
                        double &t, 
                        double  epsilon = 1e-10);


////////////////////////////////////////


/// @brief Slightly modified version of the algorithm described in "Geometric Tools for
/// Computer Graphics" pg 376 to 382 by Schneider and Eberly. Extended to handle
/// the case of a degenerate triangle. Also returns barycentric rather than
/// (s,t) coordinates. 
/// 
/// Basic Idea (See book for details): 
///
/// Write the equation of the line as 
///
///     T(s,t) = v0 + s*(v1-v0) + t*(v2-v0)
///
/// Minimize the quadratic function 
///
///     || T(s,t) - point || ^2 
///
/// by solving for when the gradient is 0. This can be done without any 
/// square roots. 
///
/// If the resulting solution satisfies 0 <= s + t <= 1, then the solution lies
/// on the interior of the triangle, and we are done (region 0). If it does 
/// not then the closest solution lies on a boundary and we have to solve for 
/// it by solving a 1D problem where we use one variable as free say "s" and
/// set the other variable t = (1-s) 
///
/// @return the closest point on the triangle and barycentric coordinates.
OPENVDB_API double
sTri3ToPointDistSqr(const Vec3d &v0,
                    const Vec3d &v1,
                    const Vec3d &v2,
                    const Vec3d &point,
                          Vec2d &uv,
                          double epsilon = 1e-10);
                          

////////////////////////////////////////


/// @return the closest point on the triangle.
static inline double
triToPtnDistSqr(const Vec3d &v0,
                const Vec3d &v1,
                const Vec3d &v2,
                const Vec3d &point)
{
    Vec2d uv;
    return sTri3ToPointDistSqr(v0,v1,v2, point, uv);
}


} // namespace math
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_TOOLS_MESH_TO_VOLUME_UTIL_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
