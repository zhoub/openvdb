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

#ifndef OPENVDB_VIEWER_VIEWER_HAS_BEEN_INCLUDED
#define OPENVDB_VIEWER_VIEWER_HAS_BEEN_INCLUDED

#include "RenderModules.h"
#include <openvdb/openvdb.h>
#include <tbb/mutex.h>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <vector>

#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <GL/glfw.h>


class Viewer
{
public:
    static void init(bool verbose = false);
    static void view(const std::vector<openvdb::GridBase::ConstPtr>&,
        int width = 900, int height = 800);
    
    //@{
    /// input callback functions
    static void keyCallback(int key, int action);
    static void mouseButtonCallback(int button, int action);
    static void mousePosCallback(int x, int y);
    static void mouseWheelCallback(int pos);
    //@}
    
private:    
    struct Camera;
    typedef boost::shared_ptr<Camera> CameraPtr;
    typedef boost::shared_ptr<RenderModule> RenderModulePtr;
    
    Viewer();
    static Viewer* sViewer;
    
    void render();
    
    CameraPtr mCamera;
    size_t mFirstRenderModule, mSecondRenderModule;
    std::vector<RenderModulePtr> mRenderModules;
};

#endif // OPENVDB_VIEWER_VIEWER_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
