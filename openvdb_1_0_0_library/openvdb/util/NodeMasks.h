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
/// @author Ken Museth
///
/// @file NodeMasks.h

#ifndef OPENVDB_UTIL_NODEMASK_HAS_BEEN_INCLUDED
#define OPENVDB_UTIL_NODEMASK_HAS_BEEN_INCLUDED

#include <cassert>
#include <cstring>
#include <iostream>// for cout
#include <openvdb/Types.h>


namespace openvdb {
OPENVDB_USE_VERSION_NAMESPACE
namespace OPENVDB_VERSION_NAME {
namespace util {

inline Index32
CountOn(Index32 v)
{
    v = v - ((v >> 1) & 0x55555555);
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
}

inline Index32 CountOff(Index32 v) { return 32-CountOn(v); }


/// @return least significant bit of v
inline Index32
FindLowestOn(Index32 v)
{
    static const Index32 MultiplyDeBruijnBitPosition[32] = {
        0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
        31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
    };
    return MultiplyDeBruijnBitPosition[Index32((v & -v) * 0x077CB531U) >> 27];
}

/// @return most significant bit of v
inline Index32
FindHighestOn(Index32 v)
{
    static const Index32 MultiplyDeBruijnBitPosition[32] = {
        0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
        8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
    };
    v |= v >> 1; // first round down to one less than a power of 2
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return MultiplyDeBruijnBitPosition[Index32(v * 0x07C4ACDDU) >> 27];
}


// This class is 32-bit specefic, hence the use if Index32 vs Index!
template<Index Log2Dim>
class NodeMask
{
public:
    static const Index32 LOG2DIM  = Log2Dim;
    static const Index32 DIM      = 1<<Log2Dim;
    static const Index32 BIT_SIZE = 1<<3*Log2Dim;
    static const Index32 INT_SIZE = (1<<(3*Log2Dim-3))>>2;
    BOOST_STATIC_ASSERT( Log2Dim>1 );

protected:
    Index32  mBits[INT_SIZE];// only member data!!

public:
    // Default is always off!
    NodeMask() { this->setOff(); }

    NodeMask(bool on) { this->set(on); }

    NodeMask(const NodeMask &B) { *this = B; }

    ~NodeMask() {}

    void operator = (const NodeMask &B)
    {
        for (Index32 i=0; i<INT_SIZE; ++i) mBits[i]=B.mBits[i];
    }

    class BaseIterator
    {
    protected:
        Index32          mPos;//bit position
        const NodeMask*  mParent;//this iterator can't change the parent_mask!
    public:
        BaseIterator() : mPos(BIT_SIZE), mParent(NULL) {}
        BaseIterator(Index32 pos,const NodeMask *parent) :
            mPos(pos), mParent(parent) {
            assert( (parent==NULL && pos==0 ) ||  (parent!=NULL && pos<=BIT_SIZE) );
        }
        bool operator==(const BaseIterator &iter) const {return mPos == iter.mPos;}
        bool operator!=(const BaseIterator &iter) const {return mPos != iter.mPos;}
        bool operator< (const BaseIterator &iter) const {return mPos <  iter.mPos;}
        void operator=(const BaseIterator  &iter) {
            mPos    = iter.mPos;
            mParent = iter.mParent;
        }
        Index32 offset() const {return mPos;}

        Index32 pos() const {return mPos;}

        bool test() const {
            assert(mPos  <= BIT_SIZE);
            return (mPos != BIT_SIZE);
        }

        operator bool() const {return this->test();}
    }; // class BaseIterator

    /// @note This happens to be a const-iterator!
    class OnIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!
    public:
        OnIterator() : BaseIterator() {}
        OnIterator(Index32 pos,const NodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos=mParent->findNextOn(mPos+1);
            assert(mPos <= BIT_SIZE);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return true;}
        OnIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class OnIterator

    class OffIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!
    public:
        OffIterator() : BaseIterator()  {}
        OffIterator(Index32 pos,const NodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos=mParent->findNextOff(mPos+1);
            assert(mPos <= BIT_SIZE);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return false;}
        OffIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class OffIterator

    class DenseIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!

    public:
        DenseIterator() : BaseIterator() {}
        DenseIterator(Index32 pos,const NodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos += 1;//careful - the increment might go beyond the end
            assert(mPos<= BIT_SIZE);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return mParent->isOn(mPos);}
        DenseIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class DenseIterator

    OnIterator beginOn() const       { return OnIterator(this->findFirstOn(),this); }
    OnIterator endOn() const         { return OnIterator(BIT_SIZE,this); }
    OffIterator beginOff() const     { return OffIterator(this->findFirstOff(),this); }
    OffIterator endOff() const       { return OffIterator(BIT_SIZE,this); }
    DenseIterator beginDense() const { return DenseIterator(0,this); }
    DenseIterator endDense() const   { return DenseIterator(BIT_SIZE,this); }

    bool operator == (const NodeMask &B) const {
        for (Index32 i=0; i<INT_SIZE; ++i) if (mBits[i] !=  B.mBits[i]) return false;
        return true;
    }

    bool operator != (const NodeMask &B) const {
        for (Index32 i=0; i<INT_SIZE; ++i) if (mBits[i] !=  B.mBits[i]) return true;
        return false;
    }

    //
    // Bitwise logical operations
    //
    NodeMask operator!() const { NodeMask m(*this); m.toggle(); return m; }
    const NodeMask& operator&=(const NodeMask& other) {
        for (Index32 i = 0; i < INT_SIZE; ++i) mBits[i] &= other.mBits[i];
        return *this;
    }
    const NodeMask& operator|=(const NodeMask& other) {
        for (Index32 i = 0; i < INT_SIZE; ++i) mBits[i] |= other.mBits[i];
        return *this;
    }
    const NodeMask& operator^=(const NodeMask& other) {
        for (Index32 i = 0; i < INT_SIZE; ++i) mBits[i] ^= other.mBits[i];
        return *this;
    }
    NodeMask operator&(const NodeMask& other) const { NodeMask m(*this); m &= other; return m; }
    NodeMask operator|(const NodeMask& other) const { NodeMask m(*this); m |= other; return m; }
    NodeMask operator^(const NodeMask& other) const { NodeMask m(*this); m ^= other; return m; }

    Index32 getMemUsage() const {return sizeof(this);}

    Index32 countOn() const {
        Index32 n=0;
        for (Index32 i=0; i< INT_SIZE; ++i) n += CountOn(mBits[i]);
        return n;
    }

    Index32 countOff() const { return BIT_SIZE-this->countOn(); }

    void setOn(Index32 i) {
        assert( (i>>5) < INT_SIZE);
        mBits[i>>5] |=  1<<(i&31);
    }

    void setOff(Index32 i) {
        assert( (i>>5) < INT_SIZE);
        mBits[i>>5] &=  ~(1<<(i&31));
    }

    void set(Index32 i, bool On) { On ? this->setOn(i) : this->setOff(i); }

    void set(bool On) {
        const Index32 val = On ? 0xFFFFFFFF : 0x00000000;
        for (Index32 i=0; i<INT_SIZE; ++i) mBits[i] = val;
    }

    void setOn()  { for (Index32 i=0; i<INT_SIZE; ++i) mBits[i]=0xFFFFFFFF; }
    void setOff() { for (Index32 i=0; i<INT_SIZE; ++i) mBits[i]=0x00000000; }
    void toggle(Index32 i) {
        assert( (i>>5) < INT_SIZE);
        mBits[i>>5] ^= 1<<(i&31);
    }
    void toggle() { for (Index32 i=0; i<INT_SIZE; ++i) mBits[i]=~mBits[i]; }
    void setFirstOn()  { this->setOn(0); }
    void setLastOn()   { this->setOn(BIT_SIZE-1); }
    void setFirstOff() { this->setOff(0); }
    void setLastOff()  { this->setOff(BIT_SIZE-1); }
    bool isOn(Index32 i) const {
        assert( (i>>5) < INT_SIZE);
        return ( mBits[i >> 5] & (1<<(i&31)) );
    }
    bool isOff(Index32 i) const {
        assert( (i>>5) < INT_SIZE);
        return ( ~mBits[i >> 5] & (1<<(i&31)) );
    }
    bool isOn() const { //all on
        for (Index32 i=0; i<INT_SIZE; ++i) if (mBits[i] != 0xFFFFFFFF) return false;
        return true;
    }
    bool isOff() const { //all off
        for (Index32 i=0; i<INT_SIZE; ++i) if (mBits[i] != 0x00000000) return false;
        return true;
    }
    Index32 findFirstOn() const {
        Index32 i=0;
        while(!mBits[i]) if (++i == INT_SIZE) return BIT_SIZE;//reached end
        return 32*i + FindLowestOn(mBits[i]);
    }
    Index32 findFirstOff() const {
        Index32 i=0;
        while(!(~mBits[i])) if (++i == INT_SIZE) return BIT_SIZE;//reached end
        return 32*i + FindLowestOn(~mBits[i]);
    }

    //@{
    /// @return the <i>n</i>th word of the bit mask, for a word of arbitrary size.
    template<typename WordT>
    WordT getWord(Index n) const {
        assert(n*8*sizeof(WordT) < BIT_SIZE);
        return reinterpret_cast<const WordT*>(mBits)[n];
    }
    template<typename WordT>
    WordT& getWord(Index n) {
        assert(n*8*sizeof(WordT) < BIT_SIZE);
        return reinterpret_cast<WordT*>(mBits)[n];
    }
    //@}

    void save(std::ostream& os) const {
        os.write((const char *)mBits,INT_SIZE*sizeof(Index32));
    }
    void load(std::istream& is) {
        is.read((char *)mBits,INT_SIZE*sizeof(Index32));
    }
    /// @brief simple print method for debugging
    void printInfo(std::ostream& os=std::cout) const {
        os << "NodeMask: Log2Dim="<<DIM<<" Bit-size="<<BIT_SIZE<<" Int-size="<<INT_SIZE<<std::endl;
    }
    void printBits(std::ostream& os=std::cout, Index32 max_out=80u) const {
        const Index32 n=(BIT_SIZE>max_out?max_out:BIT_SIZE);
        for (Index32 i=0; i < n; ++i) {
            if ( !(i&31) )
                os << "||";
            else if ( !(i%8) )
                os << "|";
            os << this->isOn(i);
        }
        os << "|" << std::endl;
    }
    void printAll(std::ostream& os=std::cout, Index32 max_out=80u) const {
        this->printInfo(os);
        this->printBits(os,max_out);
    }

    Index32 findNextOn(Index32 start) const {
        Index32 n=start>>5;//initiate
        if (n>=INT_SIZE) return BIT_SIZE; // check for out of bounds
        Index32 m=start&31, b=mBits[n];
        if (b & (1<<m)) return start;//simpel case: next is on
        b &= 0xFFFFFFFF << m;// mask out lower bits
        while(!b && ++n<INT_SIZE) b = mBits[n];// find next none-zero int
        return (!b ? BIT_SIZE : 32*n + FindLowestOn(b));//catch last-int=0
    }

    Index32 findNextOff(Index32 start) const {
        Index32 n=start>>5;// isolate the relevant int
        if (n>=INT_SIZE) return BIT_SIZE; // check for out of bounds
        Index32 m=start&31, b=~mBits[n];
        if (b & (1<<m)) return start;//simpel case: next is off
        b &= 0xFFFFFFFF<<m;// mask out lower bits
        while(!b && ++n<INT_SIZE) b = ~mBits[n];// find next none-zero int
        return (!b ? BIT_SIZE : 32*n + FindLowestOn(b));//catch last-int=0
    }

    Index32 memUsage() const {
        return INT_SIZE*sizeof(Index32);//in bytes
    }
};


// This class is 32-bit specefic, hence the use if Index32 vs Index!
class RootNodeMask
{
protected:
    Index32   mBitSize, mIntSize;
    Index32  *mBits;

public:
    RootNodeMask(): mBitSize(0), mIntSize(0), mBits(NULL) {}
    RootNodeMask(Index32 bit_size):
        mBitSize(bit_size), mIntSize(((bit_size-1)>>5)+1), mBits(new Index32[mIntSize])
    {
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=0x00000000;
    }
    RootNodeMask(const RootNodeMask& B):
        mBitSize(B.mBitSize), mIntSize(B.mIntSize), mBits(new Index32[mIntSize])
    {
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=B.mBits[i];
    }
    ~RootNodeMask() {delete [] mBits;}

    void init(Index32 bit_size) {
        mBitSize = bit_size;
        mIntSize =((bit_size-1)>>5)+1;
        delete [] mBits;
        mBits = new Index32[mIntSize];
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=0x00000000;
    }

    Index getBitSize() const {return mBitSize;}

    Index getIntSize() const {return mIntSize;}

    void operator = (const RootNodeMask &B) {
        if (mBitSize!=B.mBitSize) {
            mBitSize=B.mBitSize;
            mIntSize=B.mIntSize;
            delete [] mBits;
            mBits = new Index32[mIntSize];
        }
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=B.mBits[i];
    }

    class BaseIterator
    {
    protected:
        Index32             mPos;//bit position
        Index32             mBitSize;
        const RootNodeMask* mParent;//this iterator can't change the parent_mask!
    public:
        BaseIterator() : mPos(0), mBitSize(0), mParent(NULL) {}
        BaseIterator(Index32 pos,const RootNodeMask *parent)
            : mPos(pos), mBitSize(parent->getBitSize()), mParent(parent) {
            assert( pos<=mBitSize );
        }
        bool operator==(const BaseIterator &iter) const {return mPos == iter.mPos;}
        bool operator!=(const BaseIterator &iter) const {return mPos != iter.mPos;}
        bool operator< (const BaseIterator &iter) const {return mPos <  iter.mPos;}
        void operator=(const BaseIterator  &iter) {
            mPos      = iter.mPos;
            mBitSize = iter.mBitSize;
            mParent   = iter.mParent;
        }

        Index32 offset() const {return mPos;}

        Index32 pos() const {return mPos;}

        bool test() const {
            assert(mPos  <= mBitSize);
            return (mPos != mBitSize);
        }

        operator bool() const {return this->test();}
    }; // class BaseIterator

    /// @note This happens to be a const-iterator!
    class OnIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mBitSize;//bit size;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!
    public:
        OnIterator() : BaseIterator() {}
        OnIterator(Index32 pos,const RootNodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos=mParent->findNextOn(mPos+1);
            assert(mPos <= mBitSize);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return true;}
        OnIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class OnIterator

    class OffIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mBitSize;//bit size;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!
    public:
        OffIterator() : BaseIterator()  {}
        OffIterator(Index32 pos,const RootNodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos=mParent->findNextOff(mPos+1);
            assert(mPos <= mBitSize);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return true;}
        OffIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class OffIterator

    class DenseIterator: public BaseIterator
    {
    protected:
        using BaseIterator::mPos;//bit position;
        using BaseIterator::mBitSize;//bit size;
        using BaseIterator::mParent;//this iterator can't change the parent_mask!
    public:
        DenseIterator() : BaseIterator() {}
        DenseIterator(Index32 pos,const RootNodeMask *parent) : BaseIterator(pos,parent) {}
        void increment() {
            assert(mParent!=NULL);
            mPos += 1;//carefull - the increament might go beyond the end
            assert(mPos<= mBitSize);
        }
        void increment(Index n) {
            for (Index i=0; i<n && this->next(); ++i) {}
        }
        bool next() {
            this->increment();
            return this->test();
        }
        bool operator*() const {return mParent->isOn(mPos);}
        DenseIterator& operator++() {
            this->increment();
            return *this;
        }
    }; // class DenseIterator

    OnIterator beginOn() const       { return OnIterator(this->findFirstOn(),this); }
    OnIterator endOn() const         { return OnIterator(mBitSize,this); }
    OffIterator beginOff() const     { return OffIterator(this->findFirstOff(),this); }
    OffIterator endOff() const       { return OffIterator(mBitSize,this); }
    DenseIterator beginDense() const { return DenseIterator(0,this); }
    DenseIterator endDense() const   { return DenseIterator(mBitSize,this); }

    bool operator == (const RootNodeMask &B) const {
        if (mBitSize != B.mBitSize) return false;
        for (Index32 i=0; i<mIntSize; ++i) if (mBits[i] !=  B.mBits[i]) return false;
        return true;
    }

    bool operator != (const RootNodeMask &B) const {
        if (mBitSize != B.mBitSize) return true;
        for (Index32 i=0; i<mIntSize; ++i) if (mBits[i] !=  B.mBits[i]) return true;
        return false;
    }

    //
    // Bitwise logical operations
    //
    RootNodeMask operator!() const { RootNodeMask m = *this; m.toggle(); return m; }
    const RootNodeMask& operator&=(const RootNodeMask& other) {
        assert(mIntSize == other.mIntSize);
        for (Index32 i = 0, N = std::min(mIntSize, other.mIntSize); i < N; ++i) {
            mBits[i] &= other.mBits[i];
        }
        for (Index32 i = other.mIntSize; i < mIntSize; ++i) mBits[i] = 0x00000000;
        return *this;
    }
    const RootNodeMask& operator|=(const RootNodeMask& other) {
        assert(mIntSize == other.mIntSize);
        for (Index32 i = 0, N = std::min(mIntSize, other.mIntSize); i < N; ++i) {
            mBits[i] |= other.mBits[i];
        }
        return *this;
    }
    const RootNodeMask& operator^=(const RootNodeMask& other) {
        assert(mIntSize == other.mIntSize);
        for (Index32 i = 0, N = std::min(mIntSize, other.mIntSize); i < N; ++i) {
            mBits[i] ^= other.mBits[i];
        }
        return *this;
    }
    RootNodeMask operator&(const RootNodeMask& other) const {
        RootNodeMask m(*this); m &= other; return m;
    }
    RootNodeMask operator|(const RootNodeMask& other) const {
        RootNodeMask m(*this); m |= other; return m;
    }
    RootNodeMask operator^(const RootNodeMask& other) const {
        RootNodeMask m(*this); m ^= other; return m;
    }


    Index32 getMemUsage() const {return mIntSize*sizeof(Index32) + sizeof(this);}

    Index32 countOn() const {
        assert(mBits);
        Index32 n=0;
        for (Index32 i=0; i< mIntSize; ++i) n += CountOn(mBits[i]);
        return n;
    }

    Index32 countOff() const { return mBitSize-this->countOn(); }

    void setOn(Index32 i) {
        assert(mBits);
        assert( (i>>5) < mIntSize);
        mBits[i>>5] |=  1<<(i&31);
    }

    void setOff(Index32 i) {
        assert(mBits);
        assert( (i>>5) < mIntSize);
        mBits[i>>5] &=  ~(1<<(i&31));
    }

    void set(Index32 i, bool On) { On ? this->setOn(i) : this->setOff(i); }

    void setOn() {
        assert(mBits);
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=0xFFFFFFFF;
    }
    void setOff() {
        assert(mBits);
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=0x00000000;
    }
    void toggle(Index32 i) {
        assert(mBits);
        assert( (i>>5) < mIntSize);
        mBits[i>>5] ^= 1<<(i&31);
    }
    void toggle() {
        assert(mBits);
        for (Index32 i=0; i<mIntSize; ++i) mBits[i]=~mBits[i];
    }
    void setFirstOn()  { this->setOn(0); }
    void setLastOn()   { this->setOn(mBitSize-1); }
    void setFirstOff() { this->setOff(0); }
    void setLastOff()  { this->setOff(mBitSize-1); }
    bool isOn(Index32 i) const {
        assert(mBits);
        assert( (i>>5) < mIntSize);
        return ( mBits[i >> 5] & (1<<(i&31)) );
    }
    bool isOff(Index32 i) const {
        assert(mBits);
        assert( (i>>5) < mIntSize);
        return ( ~mBits[i >> 5] & (1<<(i&31)) );
    }

    bool isOn() const {
        if (!mBits) return false;//undefined is off
        for (Index32 i=0; i<mIntSize; ++i) if (mBits[i] != 0xFFFFFFFF) return false;
        return true;
    }

    bool isOff() const {
        if (!mBits) return true;//undefined is off
        for (Index32 i=0; i<mIntSize; ++i) if (mBits[i] != 0) return false;
        return true;
    }

    Index32 findFirstOn() const {
        assert(mBits);
        Index32 i=0;
        while(!mBits[i]) if (++i == mIntSize) return mBitSize;//reached end
        return 32*i + FindLowestOn(mBits[i]);
    }

    Index32 findFirstOff() const {
        assert(mBits);
        Index32 i=0;
        while(!(~mBits[i])) if (++i == mIntSize) return mBitSize;//reached end
        return 32*i + FindLowestOn(~mBits[i]);
    }

    void save(std::ostream& os) const {
        assert(mBits);
        os.write((const char *)mBits,mIntSize*sizeof(Index32));
    }
    void load(std::istream& is) {
        assert(mBits);
        is.read((char *)mBits,mIntSize*sizeof(Index32));
    }
    /// @brief simple print method for debugging
    void printInfo(std::ostream& os=std::cout) const {
        os << "RootNodeMask: Bit-size="<<mBitSize<<" Int-size="<<mIntSize<<std::endl;
    }

    void printBits(std::ostream& os=std::cout, Index32 max_out=80u) const {
        const Index32 n=(mBitSize>max_out?max_out:mBitSize);
        for (Index32 i=0; i < n; ++i) {
            if ( !(i&31) )
                os << "||";
            else if ( !(i%8) )
                os << "|";
            os << this->isOn(i);
        }
        os << "|" << std::endl;
    }

    void printAll(std::ostream& os=std::cout, Index32 max_out=80u) const {
        this->printInfo(os);
        this->printBits(os,max_out);
    }

    Index32 findNextOn(Index32 start) const {
        assert(mBits);
        Index32 n = start >> 5, m = start & 31;//initiate
        if (n>=mIntSize) return mBitSize; // check for out of bounds
        Index32 b = mBits[n];
        if (b & (1<<m)) return start;//simple case
        b &= 0xFFFFFFFF << m;// mask lower bits
        while(!b && ++n<mIntSize) b = mBits[n];// find next nonzero int
        return (!b ? mBitSize : 32*n + FindLowestOn(b));//catch last-int=0
    }

    Index32 findNextOff(Index32 start) const {
        assert(mBits);
        Index32 n = start >> 5, m = start & 31;//initiate
        if (n>=mIntSize) return mBitSize; // check for out of bounds
        Index32 b = ~mBits[n];
        if (b & (1<<m)) return start;//simple case
        b &= 0xFFFFFFFF<<m;// mask lower bits
        while(!b && ++n<mIntSize) b = ~mBits[n];// find next nonzero int
        return (!b ? mBitSize : 32*n + FindLowestOn(b));//catch last-int=0
    }

    Index32 memUsage() const {
        assert(mBits);
        return sizeof(Index32*)+(2+mIntSize)*sizeof(Index32);//in bytes
    }
}; // class RootNodeMask

} // namespace util
} // namespace OPENVDB_VERSION_NAME
} // namespace openvdb

#endif // OPENVDB_UTIL_NODEMASK_HAS_BEEN_INCLUDED

// Copyright (c) 2012-2013 DreamWorks Animation LLC
// All rights reserved. This software is distributed under the
// Mozilla Public License 2.0 ( http://www.mozilla.org/MPL/2.0/ )
