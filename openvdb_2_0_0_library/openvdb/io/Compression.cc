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

#include "Compression.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/shared_array.hpp>
#include <zlib.h>
#include <openvdb/Exceptions.h>
#include <openvdb/util/logging.h>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace io {

std::string
compressionToString(uint32_t flags)
{
    if (flags == COMPRESS_NONE) return "none";

    std::vector<std::string> words;
    if (flags & COMPRESS_ZIP) words.push_back("zipped");
    if (flags & COMPRESS_ACTIVE_MASK) {
        words.push_back("active values");
    }
    return boost::join(words, " ");
}


////////////////////////////////////////


namespace {
const int ZIP_COMPRESSION_LEVEL = Z_DEFAULT_COMPRESSION; ///< @todo use Z_BEST_SPEED?
}

void
zipToStream(std::ostream& os, const char* data, size_t numBytes)
{
    // Get an upper bound on the size of the compressed data.
    uLongf numZippedBytes = compressBound(numBytes);
    // Compress the data.
    boost::shared_array<Bytef> zippedData(new Bytef[numZippedBytes]);
    int status = compress2(
        /*dest=*/zippedData.get(), &numZippedBytes,
        /*src=*/reinterpret_cast<const Bytef*>(data), numBytes,
        /*level=*/ZIP_COMPRESSION_LEVEL);
    if (status != Z_OK) {
        std::string errDescr;
        if (const char* s = zError(status)) errDescr = s;
        if (!errDescr.empty()) errDescr = " (" + errDescr + ")";
        OPENVDB_LOG_DEBUG("zlib compress2() returned error code " << status << errDescr);
    }
    if (status == Z_OK && numZippedBytes < numBytes) {
        // Write the size of the compressed data.
        Int64 outZippedBytes = numZippedBytes;
        os.write(reinterpret_cast<char*>(&outZippedBytes), 8);
        // Write the compressed data.
        os.write(reinterpret_cast<char*>(zippedData.get()), outZippedBytes);
    } else {
        // Write the size of the uncompressed data.
        Int64 negBytes = -numBytes;
        os.write(reinterpret_cast<char*>(&negBytes), 8);
        // Write the uncompressed data.
        os.write(reinterpret_cast<const char*>(data), numBytes);
    }
}


void
unzipFromStream(std::istream& is, char* data, size_t numBytes)
{
    // Read the size of the compressed data.
    // A negative size indicates uncompressed data.
    Int64 numZippedBytes;
    is.read(reinterpret_cast<char*>(&numZippedBytes), 8);

    if (numZippedBytes <= 0) {
        // Read the uncompressed data.
        is.read(data, -numZippedBytes);
        if (size_t(-numZippedBytes) != numBytes) {
            OPENVDB_THROW(RuntimeError, "Expected to read a " << numBytes
                << "-byte chunk, got a " << -numZippedBytes << "-byte chunk");
        }
    } else {
        // Read the compressed data.
        boost::shared_array<Bytef> zippedData(new Bytef[numZippedBytes]);
        is.read(reinterpret_cast<char*>(zippedData.get()), numZippedBytes);
        // Uncompress the data.
        uLongf numUnzippedBytes = numBytes;
        int status = uncompress(
            /*dest=*/reinterpret_cast<Bytef*>(data), &numUnzippedBytes,
            /*src=*/zippedData.get(), static_cast<uLongf>(numZippedBytes));
        if (status != Z_OK) {
            std::string errDescr;
            if (const char* s = zError(status)) errDescr = s;
            if (!errDescr.empty()) errDescr = " (" + errDescr + ")";
            OPENVDB_LOG_DEBUG("zlib uncompress() returned error code " << status << errDescr);
        }
        if (numUnzippedBytes != numBytes) {
            OPENVDB_THROW(RuntimeError, "Expected to decompress " << numBytes
                << " byte" << (numBytes == 1 ? "" : "s") << ", got "
                << numZippedBytes << " byte" << (numZippedBytes == 1 ? "" : "s"));
        }
    }
}

} // namespace io
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
