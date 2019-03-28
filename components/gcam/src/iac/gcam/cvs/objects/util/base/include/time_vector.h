#ifndef _TIME_VECTOR_H_
#define _TIME_VECTOR_H_
#if defined(_MSC_VER)
#pragma once
#endif

/*
 * LEGAL NOTICE
 * This computer software was prepared by Battelle Memorial Institute,
 * hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830
 * with the Department of Energy (DOE). NEITHER THE GOVERNMENT NOR THE
 * CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LIABILITY FOR THE USE OF THIS SOFTWARE. This notice including this
 * sentence must appear on any copies of this computer software.
 * 
 * EXPORT CONTROL
 * User agrees that the Software will not be shipped, transferred or
 * exported into any country or used in any manner prohibited by the
 * United States Export Administration Act or any other applicable
 * export laws, restrictions or regulations (collectively the "Export Laws").
 * Export of the Software may require some form of license or other
 * authority from the U.S. Government, and failure to obtain such
 * export control license may result in criminal liability under
 * U.S. laws. In addition, if the Software is identified as export controlled
 * items under the Export Laws, User represents and warrants that User
 * is not a citizen, or otherwise located within, an embargoed nation
 * (including without limitation Iran, Syria, Sudan, Cuba, and North Korea)
 *     and that User is not otherwise prohibited
 * under the Export Laws from receiving the Software.
 * 
 * Copyright 2011 Battelle Memorial Institute.  All Rights Reserved.
 * Distributed as open-source under the terms of the Educational Community 
 * License version 2.0 (ECL 2.0). http://www.opensource.org/licenses/ecl2.php
 * 
 * For further details, see: http://www.globalchange.umd.edu/models/gcam/
 * 
 */

#include <cassert>
#include <algorithm>
#include <vector>
#include <iterator>

// TODO: Reduce these includes
#include "util/base/include/model_time.h"
#include "containers/include/scenario.h"
#include "util/base/include/util.h"

extern Scenario* scenario;

/*! 
* \file time_vector.h  
* \ingroup util
* \brief Header file for the TimeVector class.
* \author Josh Lurz
*/
namespace objects {
    // Need to inject the util namespace into the objects namespace.
    
    /*!
     * \brief Base class of vectors indexed by year or period.
     * \details Provides common code for year and period vectors.
     */
    template<class T>
    class TimeVectorBase {
    public:
    /*!
    * \brief Constant input iterator.
    */
    class const_iterator: public std::iterator<std::random_access_iterator_tag,
                                               T,
                                               size_t,
                                               const T&,
                                               const T*>
    {
        public:

            const_iterator();

            virtual ~const_iterator();

            const T& operator*() const;

            const T* operator->() const;

            bool operator==( const const_iterator& aOther ) const;

            bool operator!=( const const_iterator& aOther ) const;
            
            bool operator<( const const_iterator& aOther ) const;

            bool operator>( const const_iterator& aOther ) const;

            const_iterator& operator++();

            const_iterator operator++(int);

            const_iterator& operator--();

            const_iterator operator--(int);

            const_iterator operator+( const size_t aIncrement ) const;

            const_iterator operator-( const size_t aDecrement ) const;
    
            const_iterator& operator+=( const size_t aIncrement );
            
            const_iterator& operator-=( const size_t aDecrement );
            
            size_t operator-( const const_iterator& aOther ) const;

            // TODO: This constructor should be protected, but that
            // is very difficult to get right due to friend declarations
            // with derived classes.
            const_iterator( const unsigned int aPos, 
                            const TimeVectorBase* aParent );
        protected:


            //! Current index into the array.
            size_t mPos;

            //! Parent container.
            const TimeVectorBase* mParent;
        };

        /*!
        * \brief Mutable input iterator.
        */
        class iterator: public const_iterator {
        public:
            iterator();

            T& operator*();

            T* operator->();

            bool operator==( const iterator& aOther ) const;

            bool operator!=( const iterator& aOther ) const;
            
            bool operator<( const iterator& aOther ) const;

            bool operator>( const iterator& aOther ) const;

            iterator& operator++();

            iterator operator++(int);

            iterator& operator--();

            iterator operator--(int);

            iterator operator+( const size_t aIncrement ) const;

            iterator operator-( const size_t aDecrement ) const;
            
            iterator& operator+=( const size_t aIncrement );
            
            iterator& operator-=( const size_t aDecrement );
            
            size_t operator-( const iterator& aOther ) const;

            // TODO: This constructor should be protected, but that
            // is very difficult to get right due to friend declarations
            // with derived classes.
            iterator( const unsigned int aPos,
                      TimeVectorBase* aParent );
        protected:
            // Need to inform the compiler that this derived class will be
            // using base class members.
            using const_iterator::mParent;
            using const_iterator::mPos;
        };

        TimeVectorBase( const unsigned int aSize, const T aDefaultValue );
        virtual ~TimeVectorBase();
        TimeVectorBase( const TimeVectorBase& aOther );
        TimeVectorBase& operator=( const TimeVectorBase& aOther );

        bool operator==( TimeVectorBase& aOther ) const;
        
        bool operator!=( TimeVectorBase& aOther ) const;

        virtual T& operator[]( const size_t aIndex ) = 0;
        virtual const T& operator[]( const size_t aIndex ) const = 0;

        size_t size() const;
        void assign( const size_t aPositions, const T& aValue );
        const_iterator begin() const;
        const_iterator end() const;
        const_iterator last() const;
        iterator begin();
        iterator end();
        iterator last();
    protected:
        //! Dynamic array containing the data.
        T* mData;

        //! Size of the array.
        size_t mSize;
    private:
        void init( const unsigned int aSize,
                   const T aDefaultValue );

        void clear();
    };
    
    /*!
     * \brief Required default constructor.
     */
    template<class T>
    TimeVectorBase<T>::const_iterator::const_iterator()
    : mPos( 0 ), mParent( 0 ) {
    }

    /*!
     * \brief Constructor.
     * \param aPos Position of the iterator.
     * \param aParent TimeVectorBase parent.
     */
    template<class T>
    TimeVectorBase<T>::const_iterator::const_iterator( const unsigned int aPos, 
                                                       const TimeVectorBase<T>* aParent )
    : mPos( aPos ), mParent( aParent ){
    }
    
    /*!
     * \brief Destructor.
     * \note Only necessary to avoid warnings about
     * base classes without virtual destructors.
     */
    template<class T>
    TimeVectorBase<T>::const_iterator::~const_iterator(){
    }
    
    //! Get the contents of the iterator.
    template<class T>
    const T& TimeVectorBase<T>::const_iterator::operator*() const {
        assert( mParent );
        assert( mPos <= mParent->size() );
        return mParent->mData[ mPos ];
    }

    //! Get the contents of the iterator so a function can be called on it.
    template<class T>
    const T* TimeVectorBase<T>::const_iterator::operator->() const {
        assert( mParent );
        assert( mPos <= mParent->size() );
        return &mParent->mData[ mPos ];
    }

    //! Equality operator.
    template<class T>
    bool TimeVectorBase<T>::const_iterator::operator==( const typename TimeVectorBase<T>::const_iterator& aOther ) const {
        return mPos == aOther.mPos;
    }

    //! Inequality operator.
    template<class T>
    bool TimeVectorBase<T>::const_iterator::operator!=( const typename TimeVectorBase<T>::const_iterator& aOther ) const {
        return !( *this == aOther );
    }

    //! Less than operator.
    template<class T>
    bool TimeVectorBase<T>::const_iterator::operator<( const typename TimeVectorBase<T>::const_iterator& aOther ) const {
        return mPos < aOther.mPos;
    }

    //! Greater than operator.
    template<class T>
    bool TimeVectorBase<T>::const_iterator::operator>( const typename TimeVectorBase<T>::const_iterator& aOther ) const {
        return mPos > aOther.mPos;
    }

    //! Prefix increment.
    template<class T>
    typename TimeVectorBase<T>::const_iterator& TimeVectorBase<T>::const_iterator::operator++(){
        ++mPos;
        // If the end iterator is incremented it should remain at the end.
        mPos = std::min( mPos, mParent->size() );
        return *this;
    }

    //! Postfix increment.
    template<class T>
    typename TimeVectorBase<T>::const_iterator TimeVectorBase<T>::const_iterator::operator++(int){
        typename TimeVectorBase<T>::const_iterator prev = *this;
        operator++();
        return *this;
    }

    //! Prefix decrement.
    template<class T>
    typename TimeVectorBase<T>::const_iterator& TimeVectorBase<T>::const_iterator::operator--(){
        // If the decrement causes the iterator to move before the first element
        // of the array set the iterator to the end element.
        if( --mPos < 0 ){
            mPos = mParent->size();
        }
        return *this;
    }

    //! Postfix decrement.
    template<class T>
    typename TimeVectorBase<T>::const_iterator TimeVectorBase<T>::const_iterator::operator--(int){
        typename TimeVectorBase<T>::const_iterator prev = *this;
        operator--();
        return prev;
    }
    
    /*!
     * \brief Return an iterator incremented a specific number of positions from
     *        the current iterator.
     * \param aIncrement Amount by which to increment the iterator.
     */
    template<class T>
        typename TimeVectorBase<T>::const_iterator
        TimeVectorBase<T>::const_iterator::operator+( const size_t aIncrement ) const {
            // Use the increment operator to do the work.
            typename TimeVectorBase<T>::const_iterator temp = *this;
            return temp += aIncrement;
        }

    /*!
     * \brief Return an iterator decremented a specific number of positions from
     *        the current iterator.
     * \param aDecrement Amount by which to increment the iterator.
     */
    template<class T>
        typename TimeVectorBase<T>::const_iterator
        TimeVectorBase<T>::const_iterator::operator-( const size_t aDecrement ) const {
            // Use operator plus to do the work.
            return operator+( -1 * aDecrement );
        }

    /*!
     * \brief Increment the iterator by a specified number of positions.
     * \param aIncrement Amount by which to increment the iterator.
     * \return A reference to the iterator after incrementing.
     */
    template<class T>
        typename TimeVectorBase<T>::const_iterator&
        TimeVectorBase<T>::const_iterator::operator+=( const size_t aIncrement ) {
            // Check that this does not exceed the bounds of the iterator.
            mPos += aIncrement;
            if( mPos < 0 || mPos >= mParent->mSize ){
                mPos = mParent->mSize;
            }
            return *this;
        }

    /*!
     * \brief Decrement the iterator by a specified number of positions.
     * \param aDecrement Amount by which to decrement the iterator.
     * \return A reference to the iterator after decrementing.
     */
    template<class T>
        typename TimeVectorBase<T>::const_iterator&
        TimeVectorBase<T>::const_iterator::operator-=( const size_t aDecrement ) {
            // Use the increment operator to do the work.
            return operator+=( -1 * aDecrement );
        }
    
    /*! 
     * \brief Return the difference in position between two iterators.
     * \param aOther Iterator to return the difference between this iterator and,
     *                  such that i = j + n where n is the difference.
     * \return Difference in position between the two iterators.
     */
    template<class T>
    size_t TimeVectorBase<T>::const_iterator::operator-(const typename TimeVectorBase<T>::const_iterator& aOther) const {
        // Check that they have the same parent.
        // assert( mParent == aOther.mParent ); Perhaps undefined in copy operation?
        return mPos - aOther.mPos;
    }

    /*!
     * \brief Required default constructor.
     */
    template<class T>
    TimeVectorBase<T>::iterator::iterator(){
    }

    /*! \brief Constructor.
     * \param aPos Position of the iterator.
     */
    template<class T>
    TimeVectorBase<T>::iterator::iterator( const unsigned int aPos,
                                           TimeVectorBase<T>* aParent )
    :const_iterator( aPos, aParent )
    {
    }

    //! Get the contents of the iterator.
    template<class T>
    T& TimeVectorBase<T>::iterator::operator*() {
        assert( mParent );
        assert( mPos <= mParent->size() );
        return mParent->mData[ mPos ];
    }
    
    //! Get the contents of the iterator so a function can be called on it.
    template<class T>
    T* TimeVectorBase<T>::iterator::operator->() {
        assert( mParent );
        assert( mPos <= mParent->size() );
        return &mParent->mData[ mPos ];
    }

    //! Equality operator.
    template<class T>
    bool TimeVectorBase<T>::iterator::operator==( const typename TimeVectorBase<T>::iterator& aOther ) const {
        return mPos == aOther.mPos;
    }

    //! Inequality operator.
    template<class T>
    bool TimeVectorBase<T>::iterator::operator!=( const typename TimeVectorBase<T>::iterator& aOther ) const {
        return !( *this == aOther );
    }

    //! Less than operator.
    template<class T>
    bool TimeVectorBase<T>::iterator::operator<( const typename TimeVectorBase<T>::iterator& aOther ) const {
        return mPos < aOther.mPos;
    }

    //! Greater than operator.
    template<class T>
    bool TimeVectorBase<T>::iterator::operator>( const typename TimeVectorBase<T>::iterator& aOther ) const {
        return mPos > aOther.mPos;
    }

    //! Prefix increment.
    template<class T>
    typename TimeVectorBase<T>::iterator& TimeVectorBase<T>::iterator::operator++(){
        ++mPos;
        // If the end iterator is incremented it should remain at the end.
        mPos = std::min( mPos, mParent->size() );
        return *this;
    }

    //! Postfix increment.
    template<class T>
    typename TimeVectorBase<T>::iterator TimeVectorBase<T>::iterator::operator++(int){
        typename TimeVectorBase<T>::iterator prev = *this;
        ++mPos;
        // If the end iterator is incremented it should remain at the end.
        mPos = std::min( mPos, mParent->size() );
        return prev;
    }

    //! Prefix decrement.
    template<class T>
    typename TimeVectorBase<T>::iterator& TimeVectorBase<T>::iterator::operator--(){
        // If the decrement causes the iterator to move before the first element
        // of the array set the iterator to the end element.
        if( --mPos < 0 ){
            mPos = mParent->size();
        }
        return *this;
    }

    //! Postfix decrement.
    template<class T>
    typename TimeVectorBase<T>::iterator TimeVectorBase<T>::iterator::operator--(int){
        typename TimeVectorBase<T>::iterator prev = *this;
        // If the decrement causes the iterator to move before the first element
        // of the array set the iterator to the end element.
        if( --mPos < 0 ){
            mPos = mParent->size();
        }
         return prev;
    }

    /*!
     * \brief Return an iterator incremented a specific number of positions from
     *        the current iterator.
     * \param aIncrement Amount by which to increment the iterator.
     */
    template<class T>
        typename TimeVectorBase<T>::iterator
        TimeVectorBase<T>::iterator::operator+( const size_t aIncrement ) const {
            // Check that this does not exceed the bounds of the iterator.
            if( mPos + aIncrement >= mParent->mSize ){
                return const_cast<TimeVectorBase<T>*>( mParent )->end();
            }
            return typename TimeVectorBase<T>::iterator( mPos + aIncrement,
                                                         const_cast<TimeVectorBase<T>*>( mParent ) );
        }

    /*!
     * \brief Return an iterator decremented a specific number of positions from
     *        the current iterator.
     * \param aDecrement Amount by which to increment the iterator.
     */
    template<class T>
        typename TimeVectorBase<T>::iterator
        TimeVectorBase<T>::iterator::operator-( const size_t aDecrement ) const {
            // Check that this does not exceed the bounds of the iterator.
            if( mPos - aDecrement < 0 ){
                return const_cast<TimeVectorBase<T>*>( mParent )->end();
            }
            return typename TimeVectorBase<T>::iterator( mPos - aDecrement,
                                                         const_cast<TimeVectorBase<T>*>( mParent ) );
        }

	/*!
     * \brief Increment the iterator by a specified number of positions.
     * \param aIncrement Amount by which to increment the iterator.
     * \return A reference to the iterator after incrementing.
     */
    template<class T>
	typename TimeVectorBase<T>::iterator&
	TimeVectorBase<T>::iterator::operator+=( const size_t aIncrement ) {
		// Check that this does not exceed the bounds of the iterator.
		mPos += aIncrement;
		if( mPos < 0 || mPos >= mParent->mSize ){
			mPos = mParent->mSize;
		}
		return *this;
	}
	
    /*!
     * \brief Decrement the iterator by a specified number of positions.
     * \param aDecrement Amount by which to decrement the iterator.
     * \return A reference to the iterator after decrementing.
     */
    template<class T>
	typename TimeVectorBase<T>::iterator&
	TimeVectorBase<T>::iterator::operator-=( const size_t aDecrement ) {
		// Use the increment operator to do the work.
		return operator+=( -1 * aDecrement );
	}
    
    /*! 
     * \brief Return the difference in position between two iterators.
     * \param aOther Iterator to return the difference between this iterator and,
     *                  such that i = j + n where n is the difference.
     * \return Difference in position between the two iterators.
     */
    template<class T>
    size_t TimeVectorBase<T>::iterator::operator-(const typename TimeVectorBase<T>::iterator& aOther) const {
        // Check that they have the same parent.
        // assert( mParent == aOther.mParent ); Perhaps undefined in copy operation?
        return mPos - aOther.mPos;
    }

    /*!
     * \brief Constructor.
     * \param aSize Size of the TimeVectorBase. The size is immutable once
     *              constructed.
     * \param aDefaultValue Default for all values.
     */
    template<class T>
        TimeVectorBase<T>::TimeVectorBase( const unsigned int aSize,
                                           const T aDefaultValue )
    {
            init( aSize, aDefaultValue );
    }

    /*!
     * \brief Destructor which deallocates the array.
     */
    template<class T>
        TimeVectorBase<T>::~TimeVectorBase(){
            clear();
        }

    /*!
     * \brief Private member function to deallocate the memory.
     */
   template<class T>
       void TimeVectorBase<T>::clear(){
            delete[] mData;
       }

   /*! 
    * \brief Initialize the TimeVectorBase.
     * \param aSize Size of the TimeVectorBase. The size is immutable once
     *              constructed.
     * \param aDefaultValue Default for all values.
    */
   template<class T>
       void TimeVectorBase<T>::init( const unsigned int aSize,
                                     const T aDefaultValue )
   {
           mSize = aSize;
           mData = new T[ mSize ];

           // Initialize the data to the default value.
           std::uninitialized_fill( &mData[ 0 ], &mData[ 0 ] + mSize, aDefaultValue );
    }

    /*!
     * \brief Copy constructor.
     * \param aOther TimeVectorBase to copy.
     */
    template<class T>
        TimeVectorBase<T>::TimeVectorBase( const TimeVectorBase<T>& aOther ){
            init( aOther.mSize, T() );
            copy( aOther.begin(), aOther.end(), begin() );
        }

    /*!
     * \brief Assignment operator.
     * \param aOther TimeVectorBase to copy.
     * \return The newly constructed TimeVectorBase by reference(for chaining
     *         assignment).
     */
    template<class T>
        TimeVectorBase<T>& TimeVectorBase<T>::operator=( const TimeVectorBase<T>& aOther ){
            // Check for self-assignment.
            if( this != &aOther ){
                clear();
                init( aOther.size(), T() );
                copy( aOther.begin(), aOther.end(), begin() );
            }
            return *this;
        }

    /*!
     * \brief Equals operator.
     * \param aOther TimeVectorBase to check for equivalence.
     * \details TimeVectorBases are equivalent if they have the same size and
     *          all elements at respective positions are equal.
     * \return Whether the two vectors are equal.
     */
    template<class T>
        bool  TimeVectorBase<T>::operator==( TimeVectorBase& aOther ) const {
            return size() == aOther.size() && equal( begin(), end(), aOther.begin() );
        }

    /*!
     * \brief Not-equals operator.
     * \param aOther TimeVectorBase to check for dis-equivalence.
     * \details TimeVectorBases are not equivalent if they have the different
     *          sizes or all elements at respective positions are not equal.
     * \return Whether the two vectors are equal.
     */      
    template<class T>
        bool  TimeVectorBase<T>::operator!=( TimeVectorBase& aOther ) const {
            return !( *this == aOther );
        }

   /*!
    * \brief Return the size of the vector.
    * \return Size of the vector.
    */
    template<class T>
        size_t TimeVectorBase<T>::size() const {
            return mSize;
        }

    /*!
     * \brief Assign a single value to a number of positions in the vector
     *        start at position zero.
     * \param aPositions Number of positions to assign the value to.
     * \param aValue Value to assign to each position.
     */
    template<class T>
        void TimeVectorBase<T>::assign( const size_t aPositions, const T& aValue ){
            assert( aPositions <= size() );
            for( unsigned int i = 0; i < aPositions; ++i ){
                mData[ i ] = aValue;
            }
        }

   /*!
    * \brief Return the constant iterator to the first position in the vector.
    * \return Constant iterator to the first position.
    */
   template<class T>
       typename TimeVectorBase<T>::const_iterator TimeVectorBase<T>::begin() const {
           return typename TimeVectorBase<T>::const_iterator( 0, this );
       }

   /*!
    * \brief Return the constant iterator to the position past the last position
    *        in the vector.
    * \return Constant iterator to one position past the last position.
    */
   template<class T>
       typename TimeVectorBase<T>::const_iterator TimeVectorBase<T>::end() const {
           return typename TimeVectorBase<T>::const_iterator( size(), this );
       }

   /*!
    * \brief Return the constant iterator to the last position in the vector.
    * \return Constant iterator to the last position
    */
   template<class T>
       typename TimeVectorBase<T>::const_iterator TimeVectorBase<T>::last() const {
           // Return the last position. If the vector is empty make sure this is
           // the zeroth position.
           return typename TimeVectorBase<T>::const_iterator( std::max( int( size() ) - 1, 0 ), this );
       }

  /*!
   * \brief Return a mutable iterator to the first position in the vector.
   * \return Mutable iterator to the first position.
   */
   template<class T>
       typename TimeVectorBase<T>::iterator TimeVectorBase<T>::begin() {
           return typename TimeVectorBase<T>::iterator( 0, this );
       }

   /*!
    * \brief Return a mutable iterator to the position past the last position in
    *        the vector.
    * \return Mutable iterator to one position past the last position.
    */
   template<class T>
       typename TimeVectorBase<T>::iterator TimeVectorBase<T>::end() {
           return typename TimeVectorBase<T>::iterator( size(), this );
       }
    
   /*!
    * \brief Return a mutable iterator to the last position in the vector.
    * \return Mutable iterator to the last position
    */
   template<class T>
       typename TimeVectorBase<T>::iterator TimeVectorBase<T>::last() {
           // Return the last position. If the vector is empty make sure this is
           // the zeroth position.
           return typename TimeVectorBase<T>::iterator( std::max( static_cast<int>( size() ) - 1, 0 ), this );
       }

   /*
    * \brief YearVector is a fixed size vector which contains a value for each
    *        year within a range, and is indexed by year.
    * \details A YearVector is initialized with a start and end year and is
    *          constructed so that it has a value for each year between and
    *          including the start and end years. The size and end points of the
    *          YearVector cannot be changed once it is initialized. The
    *          YearVector is indexed by year, not vector position. Indexing
    *          outside the start and end years is invalid.
    */
   template<class T>
   class YearVector: public TimeVectorBase<T> {
   public:
        using TimeVectorBase<T>::begin;
        using TimeVectorBase<T>::end;
        using TimeVectorBase<T>::last;
        using TimeVectorBase<T>::size;
        using TimeVectorBase<T>::assign;
        using typename TimeVectorBase<T>::const_iterator;
        using typename TimeVectorBase<T>::iterator;

        YearVector( const unsigned int aStartYear,
                    const unsigned int aEndYear,
                    const T aDefaultValue = T() );

        const YearVector& operator=( const YearVector& aOther );

        virtual T& operator[]( const size_t aIndex );
        virtual const T& operator[]( const size_t aIndex ) const;
        typename TimeVectorBase<T>::const_iterator find( const unsigned int aIndex ) const;
        typename TimeVectorBase<T>::iterator find( const unsigned int aIndex );

   protected:
        //! First year of the array. This year is included in the array.
        unsigned int mStartYear;

        //! End year of the array. This year is included in the array.
        unsigned int mEndYear;
        
        // Declare that this class is using several base class members.
        using TimeVectorBase<T>::mData;
        using TimeVectorBase<T>::mSize;
   };

    /*!
     * \brief Constructor which sizes the vector to the specied number of years.
     * \param aStartYear First year of the array.
     * \param aEndYear End year of the array.
     * \param aDefaultValue Default value for each year. This argument is
     *        optional and defaults to the value created by the default
     *        constructor of the type.
     */
    template<class T>
    YearVector<T>::YearVector( const unsigned int aStartYear,
                               const unsigned int aEndYear,
                               const T aDefaultValue )
                               : TimeVectorBase<T>( aEndYear - aStartYear + 1,
                                                    aDefaultValue ),
                                 mStartYear( aStartYear ),
                                 mEndYear( aEndYear )
    {
    }

    /*!
     * \brief Assignment operator.
     * \param aOther YearVector to copy.
     * \return The newly constructed YearVector by reference(for chaining
     *         assignment).
     */
    template<class T>
    const YearVector<T>& YearVector<T>::operator=( const YearVector<T>& aOther ){
            // Check for self-assignment.
            if( this != &aOther ){
                mStartYear = aOther.mStartYear;
                mEndYear = aOther.mEndYear;
                TimeVectorBase<T>::operator=( aOther );
            }
            return *this;
        }

    /*!
     * \brief Operator which references data in the array.
     * \param aYear Year of the value to return.
     * \return Mutable value at the year by reference.
     */
    template<class T>
        T& YearVector<T>::operator[]( const size_t aYear ){
            /*! \pre The index must be between the start year and end year
            *        inclusive. 
            */
            assert( aYear >= mStartYear && aYear <= mEndYear );
            assert( isValidNumber( mData[ aYear - mStartYear ] ) );
            return mData[ aYear - mStartYear ];
        }
    
    /*!
     * \brief Operator which references data in the array.
     * \param aYear Year of the value to return.
     * \return Constant value at the year by reference.
     */
    template<class T>
        const T& YearVector<T>::operator[]( const size_t aYear ) const {
            /*! \pre The index must be between the start year and end year
            *        inclusive. 
            */
            assert( aYear >= mStartYear && aYear <= mEndYear );
            assert( isValidNumber( mData[ aYear - mStartYear ] ) );
            return mData[ aYear - mStartYear ];
        }

    /*!
     * \brief Find a specific year in the vector and return a constant iterator
     *        to it.
     * \param aYear Year for which to return a constant iterator.
     * \return Constant iterator for the year, the end iterator if it is not
     *         found.
     */
    template<class T>
        typename TimeVectorBase<T>::const_iterator YearVector<T>::find( const unsigned int aYear ) const {
            // Check if the year is valid.
            if( aYear < mStartYear || aYear > mEndYear ){
                return end();
            }
            // Return an iterator to the year.
            return typename TimeVectorBase<T>::const_iterator( aYear - mStartYear );
        }

    /*!
     * \brief Find a specific year in the vector and return a mutable iterator
     *        to it.
     * \param aYear Year for which to return a mutable iterator.
     * \return Mutable iterator for the year, the end iterator if it is not
     *         found.
     */
    template<class T>
        typename TimeVectorBase<T>::iterator YearVector<T>::find( const unsigned int aYear ) {
            // Check if the year is valid.
            if( aYear < mStartYear || aYear > mEndYear ){
                return end();
            }
            // Return an iterator to the year.
            return typename TimeVectorBase<T>::iterator( aYear - mStartYear, this );
        }

    /*!
     * \brief Array which when constructed automatically sizes to the maximum
     *          number of periods.
     * \details Array which when created will automatically have one position
     *          per period as defined by the Modeltime object.
     * \note This class is especially useful when using arrays by period within
     *       maps, since maps automatically construct their elements.
     */
    template<class T>
    class PeriodVector: public TimeVectorBase<T> {
    public:
        using typename TimeVectorBase<T>::const_iterator;
        using typename TimeVectorBase<T>::iterator;
        using TimeVectorBase<T>::begin;
        using TimeVectorBase<T>::end;
        using TimeVectorBase<T>::last;
        using TimeVectorBase<T>::size;
        using TimeVectorBase<T>::assign;

        PeriodVector( const T aDefaultValue = T() );
        virtual T& operator[]( const size_t aIndex );
        virtual const T& operator[]( const size_t aIndex ) const;
    protected:
        // Declare that this class is using the base class data and size.
        using TimeVectorBase<T>::mData;
        using TimeVectorBase<T>::mSize;
    };

    /*!
     * \brief Constructor which sizes the vector to the number of periods in the
     *        model.
     * \param aDefaultValue Default value for each year. This argument is
     *        optional and defaults to the value created by the default
     *        constructor of the type.
     */
    template<class T>
        PeriodVector<T>::PeriodVector( const T aDefaultValue )
        :TimeVectorBase<T>( scenario->getModeltime()->getmaxper(),
                            aDefaultValue )
    {
    }

    /*!
     * \brief Operator which references data in the array.
     * \param aIndex Index of the value to return.
     * \return Mutable value at the index by reference.
     */
    template<class T>
        T& PeriodVector<T>::operator[]( const size_t aIndex ){
            assert( aIndex < size() );
            assert( isValidNumber( mData[ aIndex ] ) );
            return mData[ aIndex ];
        }
    
    /*!
     * \brief Operator which references data in the array.
     * \param aIndex Index of the value to return.
     * \return Constant value at the index by reference.
     */
    template<class T>
        const T& PeriodVector<T>::operator[]( const size_t aIndex ) const {
            assert( aIndex < size() );
            assert( isValidNumber( mData[ aIndex ] ) );
            return mData[ aIndex ];
        }
    }
    /*!
     * \brief Convert a PeriodVector into a std::vector.
     * \return Vector equivalent to the period vector.
     * \todo Remove this after the database code is removed.
     */
    template<class T>
    static std::vector<double> convertToVector( const objects::PeriodVector<T>& aTimeVector ) {
        std::vector<double> convVector( aTimeVector.size() );
        copy( aTimeVector.begin(), aTimeVector.end(), convVector.begin() );
        return convVector;
   }

#endif // _TIME_VECTOR_H_
