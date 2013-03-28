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

#include "Proximity.h"

namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace math {


////////////////////////////////////////


double 
sLineSeg3ToPointDistSqr(const Vec3d &p0, 
                        const Vec3d &p1, 
                        const Vec3d &point, 
                        double &t, 
                        double  epsilon)
{
    Vec3d  pDelta;
    Vec3d  tDelta;
    double pDeltaDot;
    
    pDelta.sub(p1, p0); 
    tDelta.sub(point, p0); 

    //
    // Line is nearly a point check end points
    //
    pDeltaDot = pDelta.dot(pDelta);
    if (pDeltaDot < epsilon) {
        pDelta.sub(p1, point);
        if (pDelta.dot(pDelta) < tDelta.dot(tDelta)) {
            t = 1;
            return pDelta.dot(pDelta);
        } else {
            t = 0;
            return tDelta.dot(tDelta);
        }
    } 
    t = tDelta.dot(pDelta) / pDeltaDot;
    if (t < 0) {
        t = 0;
    } else if (t > 1) {
        t = 1;
        tDelta.sub(point, p1); 
    } else {
        tDelta -= t * pDelta;
    }
    return tDelta.dot(tDelta);    
}


////////////////////////////////////////


double
sTri3ToPointDistSqr(const Vec3d &v0,
                    const Vec3d &v1,
                    const Vec3d &v2,
                    const Vec3d &point,
                          Vec2d &uv,
                          double epsilon)
{
    Vec3d e0, e1;
    double distSqr;
    
    e0.sub(v1, v0);
    e1.sub(v2, v0);
    
    Vec3d  delta = v0 - point;
    double a00   = e0.dot(e0);
    double a01   = e0.dot(e1);
    double a11   = e1.dot(e1);
    double b0    = delta.dot(e0);
    double b1    = delta.dot(e1);
    double c     = delta.dot(delta);
    double det   = fabs(a00*a11-a01*a01);
    double aMax  = (a00 > a11) ? a00 : a11;
    double epsilon2 = epsilon * epsilon;

    //
    // Triangle is degenerate. Use an absolute test for the length
    // of the edges and a relative test for area squared
    //
    if ((a00 <= epsilon2 && a11 <= epsilon2) || det <= epsilon * aMax * aMax) {

        double t;
        double minDistSqr;
        
        minDistSqr = sLineSeg3ToPointDistSqr(v0, v1, point, t, epsilon);
        uv[0] = 1.0 - t;
        uv[1] = t;

        distSqr = sLineSeg3ToPointDistSqr(v0, v2, point, t, epsilon);
        if (distSqr < minDistSqr) {
            minDistSqr = distSqr;
            uv[0] = 1.0 - t;
            uv[1] = 0;
        }

        distSqr = sLineSeg3ToPointDistSqr(v1, v2, point, t, epsilon);
        if (distSqr < minDistSqr) {
            minDistSqr = distSqr;
            uv[0] = 0;
            uv[1] = 1.0 - t;
        }

        return minDistSqr;
    } 

    double s = a01*b1-a11*b0;
    double t = a01*b0-a00*b1;

    if (s + t <= det ) {
        if (s < 0.0) {
            if (t < 0.0) { 
                // region 4
                if (b0 < 0.0) {
                    t = 0.0;
                    if (-b0 >= a00) {
                        s = 1.0;
                        distSqr = a00+2.0*b0+c;
                    } else {
                        s = -b0/a00;
                        distSqr = b0*s+c;
                    }
                } else {
                    s = 0.0;
                    if (b1 >= 0.0) {
                        t = 0.0;
                        distSqr = c;
                    } else if (-b1 >= a11) {
                        t = 1.0;
                        distSqr = a11+2.0*b1+c;
                    } else {
                        t = -b1/a11;
                        distSqr = b1*t+c;
                    }
                }
            } else  { 
                // region 3  
                s = 0.0;
                if (b1 >= 0.0) {
                    t = 0.0;
                    distSqr = c;
                }
                else if (-b1 >= a11) {
                    t = 1.0;
                    distSqr = a11+2.0*b1+c;
                }
                else {
                    t = -b1/a11;
                    distSqr = b1*t+c;
                }
            }
        } else if (t < 0.0)  {
            // region 5        

            t = 0.0;
            if (b0 >= 0.0) {
                s = 0.0;
                distSqr = c;
            } else if (-b0 >= a00) {
                s = 1.0;
                distSqr = a00+2.0*b0+c;
            } else {
                s = -b0/a00;
                distSqr = b0*s+c;
            }
        } else { 
            // region 0

            // minimum at interior point
            double fInvDet = 1.0/det;
            s *= fInvDet;
            t *= fInvDet;
            distSqr = s*(a00*s+a01*t+2.0*b0) +
                      t*(a01*s+a11*t+2.0*b1)+c;
        }
    } else {
        double tmp0, tmp1, numer, denom;

        if (s < 0.0)  { 
            // region 2 

            tmp0 = a01 + b0;
            tmp1 = a11 + b1;
            if (tmp1 > tmp0) {
                numer = tmp1 - tmp0;
                denom = a00-2.0*a01+a11;
                if (numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                    distSqr = a00+2.0*b0+c;
                } else {
                    s = numer/denom;
                    t = 1.0 - s;
                    distSqr = s*(a00*s+a01*t+2.0*b0) +
                              t*(a01*s+a11*t+2.0*b1)+c;
                }
            } else {
                s = 0.0;
                if (tmp1 <= 0.0) {
                    t = 1.0;
                    distSqr = a11+2.0*b1+c;
                } else if (b1 >= 0.0) {
                    t = 0.0;
                    distSqr = c;
                } else {
                    t = -b1/a11;
                    distSqr = b1*t+c;
                }
            }
        } else if (t < 0.0)  { 
            // region 6

            tmp0 = a01 + b1;
            tmp1 = a00 + b0;
            if (tmp1 > tmp0 ) {
                numer = tmp1 - tmp0;
                denom = a00-2.0*a01+a11;
                if (numer >= denom ) {
                    t = 1.0;
                    s = 0.0;
                    distSqr = a11+2.0*b1+c;
                } else {
                    t = numer/denom;
                    s = 1.0 - t;
                    distSqr = s*(a00*s+a01*t+2.0*b0) +
                              t*(a01*s+a11*t+2.0*b1)+c;
                }
            } else {
                t = 0.0;
                if (tmp1 <= 0.0) {
                    s = 1.0;
                    distSqr = a00+2.0*b0+c;
                } else if (b0 >= 0.0) {
                    s = 0.0;
                    distSqr = c;
                } else {
                    s = -b0/a00;
                    distSqr = b0*s+c;
                }
            }
        } else { 
            // region 1
            numer = a11 + b1 - a01 - b0;
            if (numer <= 0.0) {
                s = 0.0;
                t = 1.0;
                distSqr = a11+2.0*b1+c;
            } else {
                denom = a00-2.0*a01+a11;
                if (numer >= denom ) {
                    s = 1.0;
                    t = 0.0;
                    distSqr = a00+2.0*b0+c;
                } else {
                    s = numer/denom;
                    t = 1.0 - s;
                    distSqr = s*(a00*s+a01*t+2.0*b0) +
                              t*(a01*s+a11*t+2.0*b1)+c;
                }
            }
        }
    }

    // Convert s,t into barycentric coordinates
    uv[0] = 1.0 - s - t;
    uv[1] = s;

    return (distSqr < 0) ? 0.0 : distSqr;
}


} // namespace math
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
