/*
 * bastructs.h
 *
 * Copyright (C) 2019 Louis-Francois Handfield
 * e-mail: lfhandfield@gmail.com
 *
 * This program is free software; upon notification by email to the licensor
 * of the licencee identity and nature of use, the licencee can redistribute
 * this program and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2
 * of the License, or (at the licencee option) any later version. As such,
 * no further notifications are required.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 *
 * This file contain class template for basic container types, such as:
 *   Vector
 *   Heaptree
 *   Hashmap
 *   Dictionary
 */

// optionnal external dependency, which uses a modified file "GSLfunc.hpp"

#ifndef _defined_Bastruct
#define _defined_Bastruct

#define GNU_SCIENTIFIC_LIBRARY
#ifndef LFH_HAS_RUNNINGTIME_STATISTICS
#define LFH_HAS_RUNNINGTIME_STATISTICS true
#endif


#include "core.h"
#include "polymorphic.h"

using namespace std;
namespace LFHPrimitive{


template<class C = void> class Dictionary;

template<class C>
class Complex{
public:
		Tuple<C,2> data;
        typedef Complex<C> SAFETYPE;
		typedef C REAL_TYPE;
		typedef Complex<C> COMPLEX_TYPE;
		typedef Quaternion<C> QUATERNION_TYPE;
        typedef Complex<C> VECTOR_TYPE;
        typedef std::integral_constant<bool, ExCo<C>::IS_COMMUTATIVE::value > IS_COMMUTATIVE;
        typedef Tuple<C,2> IMPLICIT_TYPE;
            operator Tuple<C,2>& (){return(data);}
            operator const Tuple<C,2>& ()const{return(data);}

        Complex<C> mkVectorization()const{return(*this);}


		// ExOp section
		typedef std::integral_constant<bool, ExCo<C>::IsPOD::value > IsPOD;
		// ExOp section
		static const bool NeedsAddLink = false; // containers needs to update addresses in link registers


		ExCoMeMdEcLaRe( Complex<C> );

		Complex(){}
		Complex(C const & value) {data[0] = value; data[1] = ExCo<C>::mkZero();}
		Complex(C const & real,C const & im) {data[0] = real; data[1] = im;}
		C& operator[](unsigned int);
		const C& operator[](unsigned int) const;

		Complex<C> inverse() const;
        Complex<C>& toOne(){ExOp::toOne(data[0]);ExOp::toZero(data[1]); return(*this);}
        Complex<C>& toZero(){ExOp::toZero(data[0]);ExOp::toZero(data[1]);return(*this);}
        Complex<C>& toRand(){ExOp::toRand(data[0]);ExOp::toRand(data[1]);return(*this);}
        Complex<C> mkInverse()const{return (Complex<C>(data[0],-data[1]) /= (data[0]*data[0] + data[1]*data[1]));}
        Complex<C>& toInverse(){double norm = ExOp::pnorm(data[0])+ ExOp::pnorm(data[1]);  data[0] /= norm; data[1] /= -norm; return *this;}
        Complex<C>& toSquare(){C tmp = data[0] * data[1]; ExOp::toSquare(data[0]); data[0] -= ExOp::mkSquare(data[1]); data[1] = tmp * 2.0f; return *this; }

        Complex<C>& toTrJu(){ExOp::toNegative(data[1]); return *this;}
        Complex<C> mkTrJu() const {return Complex<C>(data[0],ExOp::mknegative(data[1]));}
		double pnorm() const;
		double norm() const {return sqrt(pnorm());}

		double sign() const;

        Complex<C>& operator+=(const Complex<C>& other){ ExOp::toAdd(data[0], other.data[0]); ExOp::toAdd(data[1], other.data[1]); return *this; }
        Complex<C>& operator-=(const Complex<C>& other){ ExOp::toSubt(data[0], other.data[0]); ExOp::toSubt(data[1], other.data[1]); return *this; }
        Complex<C>& operator*=(const Complex<C>&);
        Complex<C>& operator/=(const Complex<C>&);
        Complex<C> operator+(const Complex<C>& other) const{return Complex<C>(data[0] + other.data[0], data[1] + other.data[1]);}
        Complex<C> operator-(const Complex<C>& other) const{return Complex<C>(data[0] - other.data[0], data[1] - other.data[1]);}
        Complex<C> operator*(const Complex<C>& other) const{return Complex<C>( ((*this)[0] * other[0]) - ((*this)[1] * other[1]), ((*this)[1] * other[0])+ ((*this)[0] * other[1]) );}
        Complex<C> operator/(const Complex<C>& other) const{return Complex<C>( ((*this)[0] * other[0]) + ((*this)[1] * other[1]), ((*this)[1] * other[0])- ((*this)[0] * other[1]) ) / ExOp::pnorm(other);}

		template<class A> Complex<C>& operator+=(const Complex<A>&);
		template<class A> Complex<C>& operator-=(const Complex<A>&);
		template<class A> Complex<C>& operator*=(const Complex<A>&);
		template<class A> Complex<C>& operator/=(const Complex<A>& other){return((*this) *= (other.inverse()));}

		template<class A> Complex<C>& operator+=(const A&);
		template<class A> Complex<C>& operator-=(const A&);
		template<class A> Complex<C>& operator*=(const A&);
		template<class A> Complex<C>& operator/=(const A&);


		template<class A> Complex< typename STDRETTYPE2<C,A>::PLUS_TYPE > operator+(const Complex<A>&) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::MINU_TYPE > operator-(const Complex<A>&) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::PROD_TYPE > operator*(Complex<A> const & other) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > operator/(Complex<A> const & other) const;

		template<class A> Complex< typename STDRETTYPE2<C,A>::PROD_TYPE > operator*(A const & other) const;
		template<class A> Complex< typename STDRETTYPE2<C,A>::DIVI_TYPE > operator/(A const & other) const;

		Complex< typename ExCo<C>::NEG_TYPE > operator-() const;

		SETCMP_enum setcmp(const Complex<C> &) const;

		bool operator>(const Complex<C> &)const;
		bool operator>=(const Complex<C> &)const;
		bool operator<(const Complex<C> &)const;
		bool operator<=(const Complex<C> &)const;
		bool operator==(const Complex<C> &)const;
		bool operator!=(const Complex<C> &)const;
		inline C mkrealproj()const{return data[0];}
		inline C mkimmaproj()const{return data[1];}
		inline C mkjmmaproj()const{return data[0];}
		inline C mkkmmaproj()const{return data[1];}
};
typedef Complex<double> mycomplex;
template<class C>
class Quaternion{
public:
    Tuple<C,4> data;
    Quaternion();
    Quaternion(const C _w,const C _x,const C _y,const C _z){data[0] = _w;data[1] = _x;data[2] = _y;data[3] = _z;}

    typedef C REAL_TYPE;
    typedef Quaternion<C> SAFETYPE;
    typedef Quaternion<C> COMPLEX_TYPE;
    typedef Quaternion<C> QUATERNION_TYPE;
    typedef Quaternion<C> VECTOR_TYPE;
    typedef YESNO<false> IS_COMMUTATIVE;
    typedef Tuple<C,4> IMPLICIT_TYPE;
        operator Tuple<C,4>& (){return(data);}
        operator const Tuple<C,4>& ()const{return(data);}
    Quaternion<C> mkVectorization()const{return(*this);}
    ExCoMeMdEcLaRe( Quaternion  <C> );
    const Quaternion<C>& to_normal(const C& _i,const C& _j,const C& _k);
    const Quaternion<C>& to_normal_and_scale(const C& _i,const C& _j,const C& _k,const C& _s);
    void mk_proj_matrix(TMatrix<C,3,3> &)const;
    void mk_proj_matrix(TMatrix<C,4,4> &)const;
    const Quaternion<C>& rotateX(double);
    const Quaternion<C>& rotateY(double);
    const Quaternion<C>& rotateZ(double);

    Tuple<C,3u> mkXvector() const;
    Tuple<C,3u> mkYvector() const;
    Tuple<C,3u> mkZvector() const;
    template<class O> void wrXvector(O*) const;
    template<class O> void wrYvector(O*) const;
    template<class O> void wrZvector(O*) const;

    double pnorm() const;
    double norm() const {return sqrt(pnorm());}

    const Quaternion<C>& inverse() const;
    Quaternion<C>& toInverse(){double norm = ExOp::pnorm(data[0])+ ExOp::pnorm(data[1])+ ExOp::pnorm(data[2])+ ExOp::pnorm(data[3]);  data[0] /= norm; data[1] /= -norm; data[2] /= -norm; data[3] /= -norm; return *this;}
    Quaternion<C>& toOne(){ExOp::toOne(data[0]);ExOp::toZero(data[1]);ExOp::toZero(data[2]); ExOp::toZero(data[3]); return(*this);}
    Quaternion<C>& toZero(){ExOp::toZero(data[0]);ExOp::toZero(data[1]);ExOp::toZero(data[2]);ExOp::toZero(data[3]);return(*this);}
    Quaternion<C>& toRand(){ExOp::toRand(data[0]);ExOp::toRand(data[1]);ExOp::toRand(data[2]);ExOp::toRand(data[3]);return(*this);}
    Quaternion<C>& toUnitary();

    Quaternion<C> mkInterpolated(const Quaternion<C>& origin, float frac) const;
    Quaternion<C> mkCloseInterpolated(const Quaternion<C>& origin, float frac) const;

    template<Tuple_flag CF> Tuple<C, 3, CF> operator*(const Tuple<C, 3, CF> &val)const;
    template<class OC, Tuple_flag CF> Tuple<OC, 3, CF> operator*(const Tuple<OC, 3, CF> &val)const;

    Quaternion<C> operator-()const;

    Quaternion<C> operator*(const Quaternion<C> &val)const;
    template<class OC> Quaternion<C> operator*(const Quaternion<OC> &val)const;

    Quaternion<C> operator/(const Quaternion<C> &val)const;
    template<class OC> Quaternion<C> operator/(const Quaternion<OC> &val)const;

    Quaternion<C>& operator*=(const Quaternion<C> &val);
    Quaternion<C>& operator/=(const Quaternion<C> &val);
    Quaternion<C>& toBackMult(const Quaternion<C> &val);
    Quaternion<C>& toBackDivi(const Quaternion<C> &val);

    C dotProduct(const Quaternion<C>& other)const;


    template<class OC, class OB> void wrMatrix(TMatrix<OC, 4,4>& fout, const OB& scale, bool transpose = false, bool is_final = true) const;
    template<class OC, class OB> void wrMatrix(TMatrix<OC, 3,3>& fout, const OB& scale, bool transpose = false) const;
    const C& operator[](unsigned int) const;
    C& operator[](unsigned int);
    template<class A> Quaternion<C>& operator+=(const A& other){data[0] += other;return(*this);}
    template<class A> Quaternion<C>& operator-=(const A& other){data[0] -= other;return(*this);}
    template<class A> Quaternion<C>& operator*=(const A& other){data[0] *= other;data[1] *= other;data[2] *= other;data[3] *= other; return(*this);}
    template<class A> Quaternion<C>& operator/=(const A& other){data[0] /= other;data[1] /= other;data[2] /= other;data[3] /= other; return(*this);}
    template<class A> Quaternion<C>& operator+=(const Quaternion<A>& other);
    template<class A> Quaternion<C>& operator-=(const Quaternion<A>& other);

//	template<class A> const Quaternion<C>& operator*=(const Quaternion<A>& other);
//	template<class A> const Quaternion<C>& operator/=(const Quaternion<A>& other){return((*this) *= (other.inverse()));}

    // has NO multiplication

    inline C mkrealproj()const{return data[0];}
    inline C mkimmaproj()const{return data[1];}
    inline C mkjmmaproj()const{return data[2];}
    inline C mkkmmaproj()const{return data[3];}
    const Quaternion<C>& toUnitQuaternion(const Tuple<double, 3> value);
};




template<class KEY, class DATA, bool isMax>
class ExtremumScope{
public:
    KEY best_key;
    DATA best;
    inline ExtremumScope<KEY,DATA,isMax>& toZero();
    inline void init(const KEY& key, const DATA& data);
    inline void regist(const KEY& key, const DATA& data);
};
template<class KEY, class DATA>
class ExtremumScope<KEY,DATA,false>{
public:
    KEY best_key;
    DATA best;
    inline ExtremumScope<KEY,DATA,false>&  toZero();
    inline void init(const KEY& key, const DATA& data);
    inline void regist(const KEY& key, const DATA& data);
};
template<class C, class B>
class KeyElem{
public:
	typedef std::integral_constant<bool, ExCo<C>::IsPOD::value && ExCo<B>::IsPOD::value > IsPOD;
    typedef KeyElem<C,B> SAFETYPE;
    typedef C INDEX_TYPE;
	C k;
	B d;

	KeyElem(){}
	KeyElem(const C& _k, const B& _n);

    KeyElem<C,B>& toRand(){ExOp::toRand(k);ExOp::toRand(d); return(*this);}
    KeyElem<C,B>& toZero(){ExOp::toZero(k);ExOp::toZero(d); return(*this);}
    KeyElem<C,B>& toMemmove(KeyElem<C,B>& o){ExOp::toMemmove(k,o.k);ExOp::toMemmove(d,o.d); return(*this);}
    KeyElem<C,B>& toMemswap(KeyElem<C,B>& o){ExOp::toMemswap(k,o.k);ExOp::toMemswap(d,o.d); return(*this);}
    KeyElem<C,B>& toMemfree(){ExOp::toMemfree(k);ExOp::toMemfree(d); return(*this);}

	KeyElem<C,B>& operator=(const KeyElem<C,B> &other){k = other.k; d = other.d; return(*this);}
	template<class C2,class B2> KeyElem<C,B>& operator=(const KeyElem<C2,B2> &other){k = other.k; d = other.d; return(*this);}

	C getIndex()const{return k;}
	C& getIndex(){return k;}


	bool isValid() const;

	bool operator>(const KeyElem<C,B> &other) const {return( ExOp::isGT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isGT(d, other.d)));}
	bool operator<(const KeyElem<C,B> &other) const {return( ExOp::isLT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isLT(d, other.d)));}
	bool operator>=(const KeyElem<C,B> &other) const {return( ExOp::isGT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isGE(d, other.d)));}
	bool operator<=(const KeyElem<C,B> &other) const {return( ExOp::isLT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isLE(d, other.d)));}
	bool operator==(const KeyElem<C,B> &other) const {return( ExOp::isEQ(k, other.k) && ExOp::isEQ(d, other.d));}
	bool operator!=(const KeyElem<C,B> &other) const {return( ExOp::isNQ(k,other.k) || ExOp::isNQ(d, other.d));}
	template<class A> bool operator>(const KeyElem<C,A> &other) const {return( ExOp::isGT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isGT(d, other.d)));}
	template<class A> bool operator<(const KeyElem<C,A> &other) const {return( ExOp::isLT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isLT(d, other.d)));}
	template<class A> bool operator>=(const KeyElem<C,A> &other) const {return( ExOp::isGT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isGE(d, other.d)));}
	template<class A> bool operator<=(const KeyElem<C,A> &other) const {return( ExOp::isLT(k, other.k) || (ExOp::isEQ(k, other.k) && ExOp::isLE(d, other.d)));}
	template<class A> bool operator==(const KeyElem<C,A> &other) const {return( ExOp::isEQ(k, other.k) && ExOp::isEQ(d, other.d));}
	template<class A> bool operator!=(const KeyElem<C,A> &other) const {return( ExOp::isNQ(k,other.k) || ExOp::isNQ(d, other.d));}
    bool operator>(const C &other) const {return( ExOp::isGT(k, other));}
	bool operator<(const C &other) const {return( ExOp::isLT(k, other));}
	bool operator>=(const C &other) const {return( ExOp::isGE(k, other));}
	bool operator<=(const C &other) const {return( ExOp::isLE(k, other));}
	bool operator==(const C &other) const {return( ExOp::isEQ(k, other));}
	bool operator!=(const C &other) const {return( ExOp::isNQ(k, other));}
    template<class A> bool operator>(const A &other) const {return( ExOp::isGT(k, other));}
	template<class A> bool operator<(const A &other) const {return( ExOp::isLT(k, other));}
	template<class A> bool operator>=(const A &other) const {return( ExOp::isGE(k, other));}
	template<class A> bool operator<=(const A &other) const {return( ExOp::isLE(k, other));}
	template<class A> bool operator==(const A &other) const {return( ExOp::isEQ(k, other));}
	template<class A> bool operator!=(const A &other) const {return( ExOp::isNQ(k, other));}

    template<class A,class D> bool operator>(const KeyElem<D,A> &other) const {return( ExOp::isGT(k, other.k));}
	template<class A,class D> bool operator<(const KeyElem<D,A> &other) const {return( ExOp::isLT(k, other.k));}
	template<class A,class D> bool operator>=(const KeyElem<D,A> &other) const {return( ExOp::isGE(k, other.k));}
	template<class A,class D> bool operator<=(const KeyElem<D,A> &other) const {return( ExOp::isLE(k, other.k));}
	template<class A,class D> bool operator==(const KeyElem<D,A> &other) const {return( ExOp::isEQ(k, other.k));}
	template<class A,class D> bool operator!=(const KeyElem<D,A> &other) const {return( ExOp::isNQ(k, other.k));}

    template<class A> bool operator>(const KeyElem<A,B> &other) const {return( ExOp::isGT(k, other.k));}
	template<class A> bool operator<(const KeyElem<A,B> &other) const {return( ExOp::isLT(k, other.k));}
	template<class A> bool operator>=(const KeyElem<A,B> &other) const {return( ExOp::isGE(k, other.k));}
	template<class A> bool operator<=(const KeyElem<A,B> &other) const {return( ExOp::isLE(k, other.k));}
	template<class A> bool operator==(const KeyElem<A,B> &other) const {return( ExOp::isEQ(k, other.k));}
	template<class A> bool operator!=(const KeyElem<A,B> &other) const {return( ExOp::isNQ(k, other.k));}



	void show(FILE* out = stdout, int level=0) const;

	ERRCODE load(const uint8_t * &chunk, const uint8_t* const endpos){if (ExOp::load(k,chunk,endpos) != 0) return 1; return ExOp::load(d,chunk,endpos);}
    ERRCODE save(FILE*f) const{ERRCODE fout = ExOp::save(k,f); return fout | ExOp::save(d,f);}
    ERRCODE load(FILE*f) {return ExOp::load(k,f) | ExOp::load(d,f);}
    string type_tostring()const{return string("KeyElem<") + ExOp::type_tostring(k) + string(",") + ExOp::type_tostring(d) + string(">");}
};

// base struct for recofering mean variance and higher momments of arbritrary weighted objects
template<class C,unsigned int order =1>
class WeightElem{
public:
    typedef WeightElem< C, order> SAFETYPE;
	typedef WeightElem< typename ExCo<C>::COMPLEX_TYPE, order> COMPLEX_TYPE;
	typedef WeightElem< typename ExCo<C>::REAL_TYPE, order> REAL_TYPE;
	typedef DBL_Weight WEIGHT_TYPE;

	// ExOp Section:
	typedef std::integral_constant<bool, ExCo<C>::IsPOD::value  > IsPOD;

	static const bool NeedsAddLink = false; // containers needs to update addresses in link registers

	ExCoMeMdEcLaRe( LFHCONCAT2(WeightElem<C, order>) )

	// ExOp Section end

	Tuple<double, order> w;
	Tuple<C, order> e;
	WeightElem(){}
	WeightElem(const C & ob, double weight =1.0f){
		w[0] =weight;
		for(unsigned int i=1;i<order;i++) w[i] = w[i-1] * weight;
		e[0] = ob * weight;
		for(unsigned int i=1;i<order;i++) e[i] = e[i-1] * ob;

	}
	WeightElem(const WeightElem<C,order>& clonefrom) : w(clonefrom.w), e(clonefrom.e) {}
	template<class O> WeightElem(const WeightElem<O,order>& clonefrom) : w(clonefrom.w), e(clonefrom.e) {}
	WeightElem(const Tuple<C, order>& _e,const Tuple<double, order>& _w) : w(_w), e(_e) {}

	operator C() const{return( (e[0] / w[0]) );}

	WeightElem<C,order>& toZero(){ExOp::toZero(w); ExOp::toZero(e); return(*this);}
	WeightElem<C,order>& toUndefined(){ExOp::toZero(w); ExOp::toZero(e); return(*this);}
	WeightElem<C,order>& toRand(){ExOp::toRand(e[0]); ExOp::toOne(w); for(unsigned int i=1;i<order;i++) e[i] = e[i-1] * e[0];return(*this);}
	WeightElem<C,order>& toOne(){ExOp::toOne(e); ExOp::toOne(w); return(*this);}

	WeightElem<C,order> operator+(const C& s) const {WeightElem<C,order> _out = *this; return(_out += s);}
	WeightElem<C,order> operator-(const C& s) const {WeightElem<C,order> _out = *this; return(_out -= s);}
	WeightElem<C,order> operator*(const C& s) const {WeightElem<C,order> _out = *this; return(_out *= s);}
	WeightElem<C,order> operator/(const C& s) const {WeightElem<C,order> _out = *this; return(_out /= s);}


	void replaceNan(C mean, C var, bool force_var_minimum = false);

inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkrealproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkrealproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkimmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkimmaproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkjmmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkjmmaproj(e),w);}
inline WeightElem< typename ExCo<C>::REAL_TYPE, order> mkkmmaproj()const{return WeightElem< typename ExCo<C>::REAL_TYPE, order>(ExOp::mkkmmaproj(e),w);}
WeightElem<C,order>& operator=(const WeightElem<C,order>& o){w = o.w;e = o.e;return(*this);}
WeightElem<C,order>& operator+=(const C& s){
    Tuple<C,order> buf; buf.toIntPows(s);
    if (order >2) e[1] += (buf[0] * e[1] + buf[1] * e[0]) * 3.0f + buf[2] * w[0];
    if (order >1) e[1] += buf[0] * e[0] * 2.0f + buf[1] * w[0];
    e[0] += s * w[0];
return(*this);}
WeightElem<C,order>& operator-=(const C& s){
    Tuple<C,order> buf; buf.toIntPows(s);
    if (order >2) e[1] += (buf[0] * e[1] - buf[1] * e[0]) * -3.0f - buf[2] * w[0];
    if (order >1) e[1] += (s * e[0]) * -2.0f + ExCo<C>::intPow(s,2) * w[0];
    e[0] -= s * w[0];
return(*this);}
WeightElem<C,order>& operator*=(const C& s){
    Tuple<C,order> buf; buf.toIntPows(s);
    if (order >2) e[2] *= buf[2];
    if (order >1) e[1] *= buf[1];
    e[0] *= s;
return(*this);}
template<class O> WeightElem<C,order>& operator*=(const O& s){
    Tuple<O,order> buf; buf.toIntPows(s);
    if (order >2) e[2] *= buf[2];
    if (order >1) e[1] *= buf[1];
    e[0] *= s;
return(*this);}
WeightElem<C,order>& operator+=(const WeightElem<C,order>& o ){w += o.w;  e += o.e;return(*this);}
WeightElem<C,order>& operator-=(const WeightElem<C,order>& o ){w += o.w;  e -= o.e;return(*this);}
WeightElem<C,order>& operator*=(const DBL_Weight& _w){
    double o = _w;
    double f = o;
    e *= o;
    w[0] *=o;
    for(unsigned int i=1;i<order;i++) {f *=o; w[i] *= f;}
return(*this);}
WeightElem<C,order>& operator/=(const DBL_Weight& _w ){
    double o = 1.0f / _w;
    double f = o;
    e *= o;
    w[0] *=o;
    for(unsigned int i=1;i<order;i++) {f *=o; w[i] *= f;}
    return(*(WeightElem<C,order>*)this);
    }
	WeightElem<C,order> operator+(const WeightElem<C,order>& o )const {return(WeightElem<C,order>(e +o.e,w +o.w)); }
	WeightElem<C,order> operator-(const WeightElem<C,order>& o )const;
	template<class O> WeightElem<C,order> operator+(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out += o);}
	template<class O> WeightElem<C,order> operator-(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out -= o);}
	template<class O> WeightElem<C,order> operator*(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out *= o);}
	template<class O> WeightElem<C,order> operator/(const O&o)const{WeightElem<C,order> f_out = *this; return (f_out /= o);}

void setWeight(double n_w){e *= (n_w / w[0]); if (order >2) w[1] *= (n_w / w[0]); w[0] = n_w; }

	C getMean() const {	return( (e[0] / w[0]) );}
	C getVar() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1])) );}
	C getVar_biaised() const {return( (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0])) );}

	// adds a weight of this->w[1] / this->w[0]*this->w[0] for a given variance
//	C getVar_underPrior(C prior) const {return(  (this->w[1] >= this->w[0]*this->w[0]) ? prior : (this->e[1]*this->w[0] - this->e[0] * this->e[0] + prior * this->w[1] ) * (1.0f / (this->w[0]*this->w[0]))    );}
	// adds a weight of (this->w[1] / this->w[0]*this->w[0])^2 for a given variance
    C getVar_underPrior(C prior) const {double ratio = this->w[1] / (this->w[0]*this->w[0]); return(  (ratio > 0.5f)?  prior :  (this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (1.0f / (this->w[0]*this->w[0] - this->w[1]))) ;}

	bool operator<(const WeightElem<C,order> &other){
		return( w[0] < other.w[0] );
	}
	C getSkew() const {
//		LinkAssert< (order>2) > ass;
		double w2 = this->w[0]*this->w[0];
		/*
		 printf("%f\t%f\t%f\n",this->e[3],this->w[0],this->e[3]/this->w[0]);
		 printf("%f\t%f\t%f\n",this->e[4],this->w[1],this->e[4]/this->w[1]);
		 printf("%f\t%f\t%f\n",this->e[5],this->w[2],this->e[5]/this->w[2]);

		 printf("%f\t%f\t%f\n",- (this->e[0]*this->e[2] - this->e[5])*3.0f,this->w[0]*this->w[1] - this->w[2],- (this->e[0]*this->e[2] - this->e[5])*3.0f/(this->w[0]*this->w[1] - this->w[2]));
		 printf("%f\t%f\t%f\n",- (this->e[0]*this->e[1] - this->e[4])*3.0f,w2 - this->w[1],- (this->e[0]*this->e[1] - this->e[4])*3.0f/(w2 - this->w[1]));


		 printf("%f\t%f\t%f\n",(this->e[0]*(this->e[0]*this->e[0] - this->e[2]*3.0f) + this->e[5] * 2.0f )*2.0f , (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]),(this->e[0]*(this->e[0]*this->e[0] - this->e[2]*3.0f) + this->e[5] * 2.0f )*2.0f / (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]));

		 printf("%f\t%f\t%f\n",this->w[0]*this->w[0]*this->e[3] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0], (this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]), (this->w[0]*this->w[0]*this->e[3] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0])/ ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2])));
		 */
		/*
		 return( this->e[3] * (1 / this->w[0])
		 - (this->e[0]*this->e[1] - this->e[4])*(3.0f / ( w2 - this->w[1]))
		 + (this->e[0]*(this->e[0]*this->e[0] - this->e[2]*2.0f) + this->e[5] * 1.0f )*(2.0f/(this->w[0]*(w2 -2.0f*this->w[1]) + 1.0f* this->w[2]) )
		 );

		 */
		//	 return((w2*this->e[2] -3*this->w[0]*this->e[0]*this->e[1] +2*this->e[0]*this->e[0]*this->e[0])
		//	 / ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2])));
		//		return(((this->e[2]* w2) +this->e[0]*this->e[1]*(-3.0f*this->w[0]) +this->e[0]*this->e[0]*this->e[0]*2.0f) *
		//			   (1.0f / ((this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]))));

				return((this->e[2]* (1.0f/this->w[0])) +this->e[0]*this->e[1]*(-3.0f/w2) +this->e[0]*this->e[0]*this->e[0]*(2.0f / (w2 * this->w[0]))) ;
		//		return( (this->e[2] +this->e[0] * (
		//				((this->e[1]*this->w[0] - this->e[0] * this->e[0]) * (-3.0f / (this->w[0]*this->w[0] - this->w[1]))) + this->e[0]*this->e[0]*(2.0f /(this->w[0]*this->w[0]))
		//			   ) )  * (1.0f/this->w[0])) ;
	}
	C getVar_scaleinv() const {	double v = getVar();if (v == 0.0f) return(0.0f);return sqrt(v) / getMean();}
	C getSkew_scaleinv() const {double v = getVar();if (v == 0.0f) return(0.0f);return getSkew() * pow(v, -1.5f);}
	C getKurt_scaleinv() const {double v = getVar();if (v == 0.0f) return(0.0f); return getKurt() * pow(v, -2.0f);}

	C getKurt() const {
//		LinkAssert< (order>3) > ass;
		//$K = \frac{V_4}{C_1}
		//- \frac{4V_1V_3 + 3V_2^2 - 7W_4}{C_1^2 - C_2}
		// + 12*\frac{V_1^2V_2 - 2V_1W_3 - V_2W_2 + 2W_4}{C_1^3- 3C_1C_2+ 2C_3}
		// - 6*\frac{V_1^4 - 6V_1^2W_2 +8V_1W_3 + 3W_2^2 -6W_4}{C_1^4 - 6C_1^2C_2 + 8C_1C_3 +3C_2^2 -6C_4}$

		/*
		 return( this->e[6] * (1 / this->w[0])
		 - (this->e[0]*this->e[3]*4.0f + this->e[1]*this->e[1]*3.0f + this->e[7]* -7.0f)*(1.0f / ( w2 - this->w[1]))
		 + ( this->e[1] * (e2 - this->e[2]) + this->e[0]* this->e[4] * -2.0f + this->e[8] * 2.0f )
		 *(12.0f/(this->w[0]*(w2 -3.0f*this->w[1]) + 2.0f* this->w[2]) )
		 - ( e2*(e2 + this->e[2]*-6.0f) + this->e[0] * this->e[5] *8.0f + this->e[2] *this->e[2] *3.0f + this->e[9] * -6.0f )
		 *(6.0f/(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]))
		 );

		 */
		double w2 = this->w[0]*this->w[0];

		C e2 = this->e[0]*this->e[0] * (1.0f / w2) ;
		/*
		double den = (1.0f /(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]) );
		return( ((this->e[3]*(this->w[0]) + this->e[2]*this->e[0]*(-4.0f) + this->e[1]*this->e[1]*(-3.0f))*w2 +this->e[1]*e2*(12.0f*this->w[0]) +e2*e2*(-6.0f) ) *den);*/

	//	return( this->e[3]*(1.0f / this->w[0]) + this->e[2]*this->e[0]*(-4.0f / w2) + this->e[1]*this->e[1]*(-3.0f / w2) +this->e[1]*e2*(12.0f /this->w[0]) +e2*e2*(-6.0f) );
		return( this->e[3]*(1.0f / this->w[0]) + this->e[2]*this->e[0]*(-4.0f / w2) +e2* (this->e[1]*(6.0f /this->w[0]) +e2*(-3.0f)) );
		//		 return( (w2*this->w[0]*this->e[3] -4*w2*this->e[2]*this->e[0] +6*this->w[0]*this->e[1]*e2 -3*e2*e2 )
		//		 /(w2*w2 + 8.0f*this->w[0]*this->w[2] - 3.0f*this->w[1]*(2.0f*w2 - this->w[1]) - 6*this->w[3]) );
	}

	C getSecondMomment(){ // E[(X - mu)^2] , biased variance
		return (this->e[1]  - this->e[0] * this->e[0] * (1.0f / this->w[0])) * (1.0f / this->w[0]);
	}

	C getFouthMomment(){ // E[(X - mu)^4]
		C e2 = this->e[0] * (1.0f / this->w[0]);
		return (this->e[3] - e2 * (this->e[2] * 4.0f - e2 * (this->e[1] * 6.0f - this->e[0] * e2 * 3.0f))) * (1.0f / this->w[0]);
	}


	C operator[](int i){
		if (i>order) return(ExCo<C>::zero());
		switch(i){
			case 1: return(this->getMean());
			case 2: return(this->getVar());
			case 3: return(this->getSkew());
			case 4: return(this->getKurt());
			default: return(this->getMean());
		}
	}

	class OpMean : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, order>));
	class OpMom : LFHDECL_OPER2(LFHCONCAT2(Tuple<C, order>),LFHCONCAT2(WeightElem<C, order>));
	class OpVar : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, 2>));
	class OpStd : LFHDECL_OPER2(C,LFHCONCAT2(WeightElem<C, 2>));


	double getHDist(const WeightElem<C,order> &other){return sqrt(1.0f - exp(-getBattDist(other)));}
	double getBattDist(const WeightElem<C,order> &other) const {
	LinkAssert< (order>1) > ass;
		C diff = (other.e[0] * this->w[0] - this->e[0] * other.w[0]);
		double d = ExOp::pnorm(diff) / (w[0] * other.w[0]);
		double v1 = ExOp::norm(this->getVar());
		double v2 = ExOp::norm(other.getVar());

		d /= 0.5f * this->w[0] * other.w[0] * (v1+v2);
		return((d - log (v1) -log(v2)) * 0.5f + log(v1+v2) - log(2));
	}


	void show(FILE* out = stdout, int level=0){
		switch(level){
			case 0:
				fprintf(out,"Weight: %f\n", w[0]);
				fprintf(out,"Mean: ", w[0]); ExOp::show(getMean(),out, 1);
				if (order>1) {fprintf(out,"\nVar: ", w[0]); ExOp::show(getVar(),out, 1);}
				fprintf(out,"\n");
			break;
			case 1:

				break;
		}
	}
	double getWeight() const{return w[0];}
//	Weight getWeight() const{return Weight(w[0]);}
};
// Anything can be int double char*, and optionally have an integer index
// The type itself holds the size of an unit within the 0xFF byte
template<class C, int ISCONST>
class LinearMemoryIterator{
public:
	typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
	typename MetaType<C,ISCONST>::IS_CONST_PTR range_start;
	typename MetaType<C,ISCONST>::IS_CONST_PTR range_end;
	typename MetaType<C,ISCONST>::IS_CONST_PTR cur;
	LinearMemoryIterator(typename MetaType<C,ISCONST>::IS_CONST_PTR _start, typename MetaType<C,ISCONST>::IS_CONST_PTR _end): range_start(_start), range_end(_end){}
	operator bool (){if (range_start == NULL) return false; cur = range_start; return true;}
	bool operator++(int){return ((++cur) != range_end);}
	uint32_t operator()()const{return (uint32_t) (cur - range_start);}
	typename MetaType<C,ISCONST>::IS_CONST_PTR operator->(){return cur;}
	typename MetaType<C,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return *cur;}
	LinearMemoryIterator<C,ISCONST> mkIterator(){return *this;}
};
// for matrix operation, these are treated as a diagonal matrix (not a colunm/row vector, use TMatrix for that)
template<class C,unsigned int TSIZE, Tuple_flag Cflag>
class Tuple{
public:
	typedef std::integral_constant<bool, ExCo<C>::IsPOD::value > IsPOD;
    static const uint32_t DIMENTIONALITY = 1u;
    typedef Tuple<C,TSIZE> TUPLE_TYPE;
    typedef uint32_t INDEX_TYPE;

    static C accessTupleType(const Tuple<C,TSIZE>& val, uint32_t offset) {return val[offset];}
    static C& setTupleType(Tuple<C,TSIZE>& val,uint32_t offset){return val[offset];}

	C data[TSIZE];

	unsigned int getDims() const{return 1;}
    unsigned int getSize() const{return TSIZE;}
    typedef Tuple< C, TSIZE> SAFETYPE;
	typedef Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> REAL_TYPE;
	typedef Tuple< typename ExCo<C>::COMPLEX_TYPE, TSIZE> COMPLEX_TYPE;
    typedef GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::value, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans >,TMatrix<C, TSIZE ,TSIZE> >::TYPE > GAUS_TYPE;
	typedef std::integral_constant<bool, ExCo<C>::IS_COMMUTATIVE::value > IS_COMMUTATIVE;

	// ExOp Section:


	typedef unsigned int ITERATOR_TYPE;
    template <class S> class SUBS_INNER{public: typedef Tuple<S, TSIZE, Cflag> TYPE;};

	ExCoMeMdEcLaRe( LFHCONCAT3(Tuple<C, TSIZE, Cflag>) )

	// ExOp Section end

	Tuple(){}
//	Tuple(C const & val); // update all using single val
	Tuple(const Tuple<C,TSIZE> &clonefrom){for(unsigned int i=0;i<TSIZE;i++) data[i] =clonefrom.data[i];}
	template<class O> Tuple(const Tuple<O,TSIZE> & clonefrom){unsigned int i; for(i=0;i<TSIZE;i++) data[i] = C(clonefrom.data[i]);}

	operator const C*()const{return data;}
	operator C*(){return data;}


	LinearMemoryIterator<C, 0> mkIterator(){return LinearMemoryIterator<C, 0>(data, data + TSIZE);}
	LinearMemoryIterator<C, 1> mkIterator()const{return LinearMemoryIterator<C, 1>(data, data + TSIZE);}

	Tuple(C const * const clonefrom);


	template<class O, unsigned int oTSIZE> Tuple(Tuple<O,oTSIZE,Cflag> const & other);
	int getTSIZE() const{ return(TSIZE);}

	void fourierTransform_routine();
	void invfourierTransform_routine();

    void fourierTransform_routine2();
	void invfourierTransform_routine2();
	void pow2_bitswap_permutation();
    C* operator()(){return data;}
    const C* operator()()const{return data;}

    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int length);

   static void fourierTransform_routine(mycomplex*);
	static void invfourierTransform_routine(mycomplex*);

	static Tuple<mycomplex, TSIZE,Cflag> bluesteinWindow(int wTSIZE);
    static void bluesteinWindow(mycomplex*&,int wTSIZE);

	template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> fourierTransform(const Tuple<mycomplex, superTSIZE,Cflag>& bluewindow) const;
	template<unsigned int superTSIZE> Tuple<C, TSIZE,Cflag> invfourierTransform(const Tuple<mycomplex, superTSIZE,Cflag>& bluewindow) const;
	template<unsigned int superTSIZE> Tuple<Complex<C>, TSIZE,Cflag> fourierTransformReal() const; // C cannot be multiplied by complex

	Tuple<C,TSIZE,Cflag>& toZero(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toZero(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toOne(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toOne(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toRand(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toRand(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toMin(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toMin(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toMax(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toMax(data[i]);return(*this);}
    Tuple<C,TSIZE,Cflag>& toMemfree(){for(unsigned int i=0;i<TSIZE;i++) ExOp::toMemfree(data[i]);return(*this);}

    TMatrix<C,TSIZE,TSIZE> mkOuterProd(const Tuple<C, TSIZE> &)const;

	Tuple<C, TSIZE, Cflag> normalize() const;
	Tuple<C, TSIZE, Cflag> normalize(C const & norm) const;

	Tuple<C, TSIZE, Cflag>& toIntPows(const C& value);

	C max() const;
	C min() const;

	bool isValid() const;

	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkrealproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkrealproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkimmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkimmaproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkjmmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkjmmaproj(data[i]); return fout;}
	inline Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> mkkmmaproj()const{Tuple< typename ExCo<C>::REAL_TYPE, TSIZE> fout; for(unsigned int i=0;i<TSIZE;i++) fout[i] = ExOp::mkkmmaproj(data[i]); return fout;}

	Tuple<C,TSIZE,Cflag>& operator=(Tuple<C,TSIZE,Cflag> const & other);
	Tuple<C,TSIZE,Cflag>& operator=(const C (&other)[TSIZE]);
	template<class O> Tuple<C,TSIZE,Cflag>& operator=(O const & other);
	template<class O> Tuple<C,TSIZE,Cflag>& operator=(const O (&other)[TSIZE]);
    Tuple<C,TSIZE,Cflag>& toMemmove(Tuple<C,TSIZE,Cflag> & other);

//	template<class O> const Tuple<C,TSIZE,Cflag>& operator=(Tuple<O,TSIZE,Cflag> const & other);
	template<class O, unsigned int  OSIZE, Tuple_flag OFLAG> Tuple<C,TSIZE,Cflag>& operator=(const Tuple<O,OSIZE,OFLAG> & other);

#undef LFHTEMP
#define LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag>
	LFHTEMP char compare(const Tuple<O,oTSIZE, Oflag>& other) const;
	LFHTEMP bool operator>(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 1);}
	LFHTEMP bool operator<(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 2);}
	LFHTEMP bool operator>=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 2);}
	LFHTEMP bool operator<=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 1);}
	LFHTEMP bool operator==(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 0);}
	LFHTEMP bool operator!=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 0);}



	template<class A_1> void operator() (Oper1<A_1> const & op); // not a match
	void operator() (Oper1<C> const & op); // match
	const Tuple<C,TSIZE,Cflag>& operator-() const;

	template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & ); // not a match
	template<class C_2, unsigned int TSIZE_2> void operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> const & ); // match
	template<class A_1, class A_2, class C_2, unsigned int TSIZE_2> void operator() (Oper2<A_1,A_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & ); // not a match
	template<class C_2, unsigned int TSIZE_2> void operator() (Oper2<C,C_2> const & op, Tuple<C_2, TSIZE_2,Cflag> & ); // match
#undef LFHTEMP
#define LFHTEMP template<class A_1, class A_2, class A_3, class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3>
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const &, Tuple<C_3, TSIZE_3,Cflag> const & ); // not a match
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> const & ); // not a match
	LFHTEMP	void operator() (Oper3<A_1,A_2,A_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> & ); // not a match
#undef LFHTEMP
#define LFHTEMP template<class C_2, unsigned int TSIZE_2, class C_3, unsigned int TSIZE_3>
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> const &, Tuple<C_3, TSIZE_3,Cflag> const & ); // match
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> const & ); // match
	LFHTEMP void operator() (Oper3<C,C_2,C_3> const & op, Tuple<C_2, TSIZE_2,Cflag> &, Tuple<C_3, TSIZE_3,Cflag> & ); // match




#undef LFHTEMP
#define LFHTEMP template<class O>

	LFHTEMP void addition(Tuple<O,TSIZE,Cflag> const & other){data[0] += other.data[0];}


	LFHTEMP Tuple<C,TSIZE,Cflag>& operator+=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator-=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(Tuple<O,TSIZE,Cflag> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(Tuple<O,TSIZE,Cflag> const & other);

	LFHTEMP auto operator+(Tuple<O,TSIZE,Cflag> const & other) const -> Tuple< decltype( ExOp::mkAdd(data[0], other.data[0]) ),TSIZE,Cflag>;
	LFHTEMP auto operator-(Tuple<O,TSIZE,Cflag> const & other) const -> Tuple< decltype( ExOp::mkSubt(data[0], other.data[0]) ),TSIZE,Cflag>;
	LFHTEMP auto operator*(Tuple<O,TSIZE,Cflag> const & other) const -> Tuple< decltype( ExOp::mkMult(data[0], other.data[0]) ),TSIZE,Cflag>;
	LFHTEMP auto operator/(Tuple<O,TSIZE,Cflag> const & other) const -> Tuple< decltype( ExOp::mkDivi(data[0], other.data[0]) ),TSIZE,Cflag>;

	// expending to all values
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator+=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator-=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(O const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(O const & other);
	LFHTEMP	Tuple< typename STDRETTYPE2<C,O>::PLUS_TYPE ,TSIZE,Cflag> operator+(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::MINU_TYPE ,TSIZE,Cflag> operator-(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::PROD_TYPE ,TSIZE,Cflag> operator*(O const & other) const; //
	LFHTEMP Tuple< typename STDRETTYPE2<C,O>::DIVI_TYPE ,TSIZE,Cflag> operator/(O const & other) const; //

	LFHTEMP Tuple<C,TSIZE,Cflag>& operator*=(TMatrix<O,TSIZE,TSIZE> const & other);
	LFHTEMP Tuple<C,TSIZE,Cflag>& operator/=(TMatrix<O,TSIZE,TSIZE> const & other);

	LFHTEMP	Tuple<C,TSIZE,Cflag> operator+(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator-(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator*(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag> operator/(KeyElem<unsigned int, O> const & other) const;
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator+=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator-=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator*=(KeyElem<unsigned int, O> const & other);
	LFHTEMP	Tuple<C,TSIZE,Cflag>& operator/=(KeyElem<unsigned int, O> const & other);

	//template<unsigned int PSIZE> Tuple<PolyThing<C,PSIZE>,TSIZE,Cflag> operator*(const PolyThing<C,PSIZE>&) const;


	inline C& operator[](int const pos);
	inline const C& operator[](int const pos) const;

	template<class D> operator Tuple<D,TSIZE> () const;
	double weight();

	inline static const Tuple<C,TSIZE,Cflag>& genConvolution_Gaussian(double std);



	// operators

	template<unsigned int pos> class Selector : LFHDECL_OPER2(C,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	template<unsigned int firsTSIZE> class Concatenate : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C,firsTSIZE,Cflag>), LFHCONCAT3(Tuple<C, TSIZE - firsTSIZE,Cflag>));
	template<unsigned int firsTSIZE> class Deconcatenate : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C,firsTSIZE,Cflag>), LFHCONCAT3(Tuple<C, TSIZE + firsTSIZE,Cflag>));


	template<unsigned int pos> class SelectorWrite : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),C);
	template<unsigned int TSIZEin, unsigned int pos_in, unsigned int pos_out> class SelectSelectorWrite : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZEin,Cflag>));


//	class FiniteDifference: LFHDECL_OPER2(LFHCONCAT(Tuple<C, TSIZE>),LFHCONCAT(Tuple<C, TSIZE>));
	class FiniteDifference:  public Oper2< Tuple<C, TSIZE,Cflag> , Tuple<C, TSIZE,Cflag> > { public: void operator()(Tuple<C, TSIZE,Cflag>  &, Tuple<C, TSIZE,Cflag> &) const;};
	class FiniteDifferenceAssign : LFHDECL_OPER1(LFHCONCAT3(Tuple<C, TSIZE,Cflag>));
	class L2Norm : LFHDECL_OPER2(C,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	class ArgMax : LFHDECL_OPER2(int,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));

	class MassCenter : LFHDECL_OPER2(double,LFHCONCAT3(Tuple<C, TSIZE,Cflag>));
	class MassCenter2 : LFHDECL_OPER2(LFHCONCAT3(Tuple<double, 3,Cflag>),LFHCONCAT3(Tuple<C, TSIZE,Cflag>));


	template<unsigned int winTSIZE>	class Convolution : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, winTSIZE,Cflag>));

	class MakeTuple : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>), C);
	class PopWeight : LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE+1,Cflag>));
	class PushWeight : LFHDECL_OPER3(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE-1,Cflag>), C);

	class VarNormalFiniteDifference:  LFHDECL_OPER2(LFHCONCAT3(Tuple<C, TSIZE,Cflag>),LFHCONCAT3(Tuple<C, TSIZE*2,Cflag>));

	C dotProduct(const Tuple<C,TSIZE, Cflag>& other)const{
		C _out = data[0] * other[0];
		int i;
		for(i=1;i<TSIZE;i++) _out += data[i] * other[i];
		return(_out);
	}

template<unsigned int order> Tuple<C,  TEMPLATE_TRIANGLE_NUMBER<TSIZE, order+1>::ans  > genProducts(){
    C partial[order-1];
    Tuple<C,  TEMPLATE_TRIANGLE_NUMBER<TSIZE,order+1>::ans  > _out;
    int coor[order];
    int cor;
    int j;
    for(j=0;j<TSIZE;j++) _out[j] = data[j];
    int k = TSIZE;
    for(cor=2;cor<=order;cor++){
        for(j=0;j<cor-1;j++) {
            coor[j] = j;
            partial[j] = (j ==0) ?  data[0] : data[j] * partial[j-1];
        }
        coor[j] = j;
        do{
            _out[k] = data[k] * partial[cor-2];k++;

            j = cor-1;
            if (coor[j] == TSIZE-1) {
                for(j--; j>=0 ;j--) if (coor[j] != coor[j+1] -1) break;
                if (j < 0) break;
                coor[j]++;
                partial[j] = (j ==0) ?  data[0] : data[j] * partial[j-1];
                for(j++;j<cor-1;j++) {
                    coor[j] = coor[j-1]+1;
                    partial[j] = data[ coor[j] ] * partial[j-1];
                }
            }
        }while(true);
    }
return(_out);}

    template<class I> Tuple<I> mkOrdering()const;
    Tuple<C,TSIZE,Cflag> mkInverse()const;
    Tuple<C,TSIZE,Cflag> mkPowInvInt(int a)const;
    Tuple<C,TSIZE,Cflag> mkPowInt(int a)const;
    Tuple<C,TSIZE,Cflag>& toInverse();


    void wrDeterminant(C&) const;
    template <class O> void wrDeterminant(O&) const;
    void wrTrace(C&) const;
    template <class O> void wrTrace(O&) const;

	double pnorm() const;
	double norm() const;
    template<class O> std::vector<O>& wrStdVector(std::vector<O>& fout)const;
    template<class O> operator std::vector<O> ()const;

    GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::value, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE > mkgaussstat(double &w) const;

    static double overlap(const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::value, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&a, const GaussScope< Tuple<C,TSIZE> , typename MT_IFTYPE<ExCo<C>::IS_COMMUTATIVE::value, Tuple<C, TEMPLATE_TRIANGLE_NUMBER<TSIZE, 2>::ans> , TMatrix<C,TSIZE,TSIZE> >::TYPE >&b);
};
template<class C, Tuple_flag Cflag>
class Tuple<C,0u, Cflag>{
    void setSize_init(uint32_t);
public:
    static const uint32_t DIMENTIONALITY = 1u;
    typedef Tuple<C,0u> TUPLE_TYPE;
    typedef uint32_t INDEX_TYPE;
    static C accessTupleType(const Tuple<C,0u>& val, uint32_t offset) {return val[offset];}
    static C& setTupleType(Tuple<C,0u>& val,uint32_t offset){return val[offset];}

    C* data;
    uint32_t tup_size;
    typedef Tuple< C, 0u> SAFETYPE;
    Tuple(): tup_size(0){}
    Tuple(const Tuple<C,0u, Cflag> &other);
    Tuple(std::initializer_list<C>);
    template<class O> Tuple(std::initializer_list<O>);


    Tuple(Tuple<C,0u, Cflag> &&other): data(other.data), tup_size(other.tup_size) {other.tup_size = 0;}
    Tuple<C,0u, Cflag>& operator=(const Tuple<C,0u, Cflag> &other){this->setSize(other.tup_size); for(unsigned int i=0;i<tup_size;i++) data[i] = other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator=(Tuple<C,0u, Cflag> && other){if (this == &other) return(*this); if (tup_size) delete[](data); data = other.data; tup_size = other.tup_size; other.tup_size=0; return(*this);}
    ~Tuple(){if (tup_size!= 0) delete[](data);}

    template<class O> Tuple(const Vector<O> &other);
    template<class O> Tuple(const std::vector<O> &other);
    template<unsigned int OSIZE> Tuple(const Tuple<C,OSIZE, Cflag> &other);

    operator RemoteMemory<const C,1> ()const{return(RemoteMemory<const C,1>(data, tup_size));}


    template<class B, unsigned int OSIZE, Tuple_flag OFLAG> Tuple<C,0u, Cflag>& operator=(const Tuple<B,OSIZE, OFLAG> &other){this->setSize(other.getSize()); for(unsigned int i=0;i<tup_size;i++) data[i] = other.data[i]; return(*this);}

	ExCoMeMdEcLaRe( LFHCONCAT3(Tuple<C, 0u, Cflag>) )

	class Iterator{
    public:
        Tuple<C,0u, Cflag>& target;
        uint32_t cur, maxval;
        typedef uint32_t KEYITERATORMKR_TYPE;
        Iterator(Tuple<C,0u, Cflag>& trg): target(trg), cur(0),maxval(trg.getSize()){}
        Iterator(Tuple<C,0u, Cflag>& trg, uint32_t _start, uint32_t _endoffset): target(trg), cur(_start),maxval(_endoffset){}
        ~Iterator(){}

        uint32_t getSize(){return target.tup_size;}
        operator LFHPrimitive::Iterator<uint32_t, C>(){printf("Iterator is to be created! itar is %p\n", (cur != maxval) ? &target[cur] :NULL); fflush(stdout); using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, C>(std::bind( &Iterator::next , this, _1,_2), &cur, (cur != maxval) ? &target[cur] :NULL );}
        operator bool (){cur = 0; return (0 != maxval);}
        bool operator++(int){return (++cur < maxval);}
        bool next(uint32_t*& keyt, C*& itar){
            printf("next was called! %i %i\n", cur, maxval);
            if (++cur == maxval) return false; printf("next was called! %i\n", cur); fflush(stdout); itar = &(target[cur]); return true;
        }

        uint32_t operator()()const{return cur;}
        C* operator->(){return target.data + cur;}
        C& operator*(){return target[cur];}
        Iterator& mkIterator(){return *this;}
	};
    class ConstIterator{
    public:
        uint32_t cur, maxval;
        const Tuple<C,0u, Cflag>& target;
        typedef uint32_t KEYITERATORMKR_TYPE;
        ConstIterator(const Tuple<C,0u, Cflag>& trg): target(trg), cur(0),maxval(trg.getSize()){}

        operator LFHPrimitive::Iterator<uint32_t, const C>(){using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, const C>(std::bind( &ConstIterator::next , this, _1,_2), &cur, (cur != maxval) ? &target[cur] :NULL );}
        bool next(uint32_t*& keyt, const C*& itar){if (++cur == maxval) return false; itar = &(target[cur]); return true;}

        operator bool (){cur = 0; return (0 != maxval);}
        bool operator++(int){return (++cur < maxval);}
        uint32_t operator()()const{return cur;}
        const C* operator->(){return target.data + cur;}
        const C& operator*(){return target[cur];}
        ConstIterator& mkIterator(){return *this;}
	};

	template<int ISCONST>
	class ModuloRangeIterator{
	public:
		typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
		uint32_t start, limit, modul, cur;
		typename MetaType<Tuple<C,0u, Cflag> ,ISCONST>::IS_CONST_REF target;
        ModuloRangeIterator(typename MetaType<Tuple<C,0u, Cflag> ,ISCONST>::IS_CONST_REF trg, uint32_t _start, uint32_t _limit, uint32_t _modul): target(trg),start(_start),limit(_limit-1),modul(_modul){}
        operator bool (){if ((target.tup_size <= start)||(start == limit + 1)||(start >= modul)) return false; cur = start; return true;}
        bool operator++(int){if ((cur % modul) == limit) cur += modul + start - (cur % modul); else cur++; return cur < target.tup_size;}
        uint32_t operator()()const{return cur;}

		typename MetaType<C,ISCONST>::IS_CONST_PTR operator->(){return target.data + cur;}
        typename MetaType<C,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target[cur];}
        ModuloRangeIterator<ISCONST>& mkIterator(){return *this;}
	};

    Tuple<C,0u, Cflag>::Iterator getIterator(){return Iterator(*this);}
    Tuple<C,0u, Cflag>::ConstIterator getIterator()const{return ConstIterator(*this);}
    Tuple<C,0u, Cflag>::Iterator mkIterator(){return Iterator(*this);}
    Tuple<C,0u, Cflag>::ConstIterator mkIterator()const{return ConstIterator(*this);}
    Tuple<C,0u, Cflag>::Iterator getPartitionIterator(int which, int nbparts){return Iterator(*this, (which * getSize()) / nbparts, ((which+1) * getSize()) / nbparts );}
    Tuple<C,0u, Cflag>::ConstIterator getPartitionIterator(int which, int nbparts)const{return ConstIterator(*this, (which * getSize()) / nbparts, ((which+1) * getSize()) / nbparts );}
	Tuple<C,0u, Cflag>::ModuloRangeIterator<1> mkModuloRangeIterator(int modulo, int start, int lenght)const {return ModuloRangeIterator<1>(*this, start,start+lenght,modulo);}
	Tuple<C,0u, Cflag>::ModuloRangeIterator<0> mkModuloRangeIterator(int modulo, int start, int lenght){return ModuloRangeIterator<0>(*this, start,start+lenght,modulo);}

    LFHPrimitive::Iterator<uint32_t, C> getPartitionPolyIterator(int which, int nbparts){return (LFHPrimitive::Iterator<uint32_t, C>) Iterator(*this, (which * getSize()) / nbparts, ((which+1) * getSize()) / nbparts );}

    operator LFHPrimitive::Iterator<uint32_t,C> (){return (LFHPrimitive::Iterator<uint32_t, C>) Iterator(*this);}
    operator LFHPrimitive::Iterator<uint32_t, const C> ()const {return (LFHPrimitive::Iterator<uint32_t, const C>) ConstIterator(*this);}

    operator IteratorMaker<uint32_t, C> (){using std::placeholders::_1; using std::placeholders::_2; return IteratorMaker<uint32_t, C>(std::bind( &getPartitionPolyIterator, this, _1,_2));}

    operator Accessor<unsigned int, C> (){using std::placeholders::_1; using std::placeholders::_2;  return LFHPrimitive::Accessor<unsigned int, C>(std::bind( &access, this, _1,_2));}
    bool access(const uint32_t &index, C*& itar){itar = data + index; return(index < tup_size);}

    inline const C& operator[](unsigned int index)const;
    inline C& operator[](unsigned int index);
    inline const C& operator[](int index)const;
    inline C& operator[](int index);


	unsigned int getDims() const{return 1;}
    unsigned int getSize() const{return tup_size;}
    Tuple<C,0u, Cflag>& setSize(unsigned int s); // discards data
    Tuple<C,0u, Cflag>& toResize(unsigned int s); // keeps overlapping data

	bool isValid() const{for(unsigned int i =0; i < tup_size;i++) if (!ExOp::isValid(data[i])) return false; return true;}

    Tuple<C,0u, Cflag>& toZero(){for(unsigned int i=0;i<tup_size;i++) ExOp::toZero(data[i]); return(*this);}
    Tuple<C,0u, Cflag>& toOne(){for(unsigned int i=0;i<tup_size;i++) ExOp::toOne(data[i]); return(*this);}
    Tuple<C,0u, Cflag>& toRand(){for(unsigned int i=0;i<tup_size;i++) ExOp::toRand(data[i]); return(*this);}
	Tuple<C,0u, Cflag>& toMemfree(){if (tup_size != 0) delete[](data); tup_size=0; return(*this);}

    Tuple<C,0u, Cflag>& operator+=(const Tuple<C,0u, Cflag> &other){unsigned int j = (tup_size < other.tup_size) ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] += other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator-=(const Tuple<C,0u, Cflag> &other){unsigned int j = (tup_size < other.tup_size) ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] -= other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator*=(const Tuple<C,0u, Cflag> &other){unsigned int j = (tup_size < other.tup_size) ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] *= other.data[i]; return(*this);}
    Tuple<C,0u, Cflag>& operator/=(const Tuple<C,0u, Cflag> &other){unsigned int j = (tup_size < other.tup_size) ? tup_size : other.tup_size; for(unsigned int i=0;i<j;i++) data[i] /= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator+=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = (tup_size < OSIZE) ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] += other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator-=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = (tup_size < OSIZE) ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] -= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator*=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = (tup_size < OSIZE) ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] *= other.data[i]; return(*this);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag>& operator/=(const Tuple<C,OSIZE, Cflag> &other){unsigned int j = (tup_size < OSIZE) ? tup_size : OSIZE; for(unsigned int i=0;i<j;i++) data[i] /= other.data[i]; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator+=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] += other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator-=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] -= other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator*=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] *= other; return(*this);}
    template<class O> Tuple<C,0u, Cflag>& operator/=(const O &other){for(unsigned int i=0;i<tup_size;i++) data[i] /= other; return(*this);}

    Tuple<C,0u, Cflag> operator+(const Tuple<C,0u, Cflag> &other)const;
    Tuple<C,0u, Cflag> operator-(const Tuple<C,0u, Cflag> &other)const;
    Tuple<C,0u, Cflag> operator*(const Tuple<C,0u, Cflag> &other)const;
    Tuple<C,0u, Cflag> operator/(const Tuple<C,0u, Cflag> &other)const;

    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator+(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) += other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator-(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) -= other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator*(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) *= other);}
    template<unsigned int OSIZE> Tuple<C,0u, Cflag> operator/(const Tuple<C,OSIZE, Cflag> &other)const{return((Tuple<C,0u, Cflag>(*this)) /= other);}

    template<class O> Tuple<C,0u, Cflag> operator+(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) += other);}
    template<class O> Tuple<C,0u, Cflag> operator-(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) -= other);}
    template<class O> Tuple<C,0u, Cflag> operator*(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) *= other);}
    template<class O> Tuple<C,0u, Cflag> operator/(const O &other)const{return((Tuple<C,0u, Cflag>(*this)) /= other);}




    template<class O, uint32_t OSIZE> auto operator+(const Tuple<O,OSIZE, Cflag> &other)const -> Tuple< decltype( ExOp::mkAdd(data[0u], other.data[0u]) ),0u,Cflag>;
    template<class O, uint32_t OSIZE> auto operator-(const Tuple<O,OSIZE, Cflag> &other)const -> Tuple< decltype( ExOp::mkSubt(data[0u], other.data[0u]) ),0u,Cflag>;
    template<class O, uint32_t OSIZE> auto operator*(const Tuple<O,OSIZE, Cflag> &other)const -> Tuple< decltype( ExOp::mkMult(data[0u], other.data[0u]) ),0u,Cflag>;
    template<class O, uint32_t OSIZE> auto operator/(const Tuple<O,OSIZE, Cflag> &other)const -> Tuple< decltype( ExOp::mkDivi(data[0u], other.data[0u]) ),0u,Cflag>;


    template<class O, uint32_t Tuple_flag, class D, uint32_t L> Tuple<C,0u, Cflag>& toMult(const TMatrix<O,0u,0u> & a, const Tuple<D,L> & b) const;
    template<class O, uint32_t Tuple_flag, class D, uint32_t L> Tuple<C,0u, Cflag>& toMult(const Tuple<D,L> & a, const TMatrix<O,0u,0u> & b) const;
    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& toMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const;
    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& toMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const;

    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& addMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const;
    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& addMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const;
    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& subtMult(const TMatrix<O,0u,0u> & a, const SparseTuple<D> & b) const;
    template<class O, uint32_t Tuple_flag, class D> Tuple<C,0u, Cflag>& subtMult(const SparseTuple<D> & a, const TMatrix<O,0u,0u> & b) const;

	//template<class I, ENABLEIF_NOT_ITERATOR(I, uint32_t) =true> void mkInnerProd(I x) const{printf("this is generic\n");}
	template<class I, ENABLEIF_ITERATOR(I, uint32_t) =true> auto mkInnerProd(I x) const -> decltype(data[0] * (*(x.mkIterator())));

	//template<class I, ENABLEIF_NOT_ITERATORMKR(I, uint32_t) =true> void mkInnerProd(const I& x) const{printf("this is generic\n");}
	//template<class I, ENABLEIF_ITERATORMKR(I, uint32_t) =true> auto mkInnerProd(const I& x) const -> decltype(data[0] * (*(x.mkIterator())));
	//template<class I, ENABLEIF_ITERATOR(I, uint32_t) =true> auto mkInnerProd(const I& x) const -> decltype(data[0] * (*(x.getIterator())));

	template<class I, ENABLEIF_WRITE_ITERATOR(I, uint32_t) =true> void mkInnerProdWrite(const I& x){printf("thiw is specific for writing\n");}
	template<class I, ENABLEIF_NOT_ITERATOR(I, uint32_t) =true> void mkInnerProdWrite(const I& x){printf("thiw is generic\n");}

    template<class I> auto mkInnerProdOLD(const I& other, ITERABLE_DECL(I)) const -> decltype(data[0] * other[0]);
    //template<class I, int ISCONST> auto mkInnerProd(KeyIterator<uint32_t, I, ISCONST>& other) const -> decltype(data[0] * (*other));
    //template<class I> void std::enable_if<ExCo<I>:: >

    template<class O> auto mkInnerProd(const SparseTuple<O>& other) const -> decltype(data[0] * other[0]);

    template<class O, uint32_t S, Tuple_flag CFLAG2> auto mkInnerProd(const Tuple<O,S,CFLAG2>& other)const -> decltype(data[0] * other.data[0]);
    template<class O, uint32_t S, Tuple_flag CFLAG2> TMatrix<C,0u,0u> mkOuterProd(const Tuple<O,S,CFLAG2>& other)const;

    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int length );


    C& last(int offset =1){return data[tup_size-offset];}
    const C& last(int offset =1)const{return data[tup_size-offset];}

	Tuple<C,0u, Cflag>& toMemmove(Tuple<C,0u, Cflag>& source){if (&source == this) return *this; if (tup_size) delete[](data); data = source.data; tup_size = source.tup_size; source.tup_size =0; return *this;}
	Tuple<C,0u, Cflag>& toMemmove(Vector<C>& source);
	Tuple<C,0u, Cflag>& operator=(const Vector<C>& source);

	Tuple<C,0u, Cflag>& permute(const Tuple<uint32_t>& permutation);

    C* operator()(){return data;}
    const C* operator()()const{return data;}
	// operator C*()const{return data;} causes ambiguity, use (*this)()

    C& push_back();
    C& push_back(const C& value){return this->push_back() = value;} // warning: the function below if the pointer is from this vector!
    C& push_back_copy(unsigned int offset){return this->push_back() =(*this)[offset];}

    void pop_back(); // O(n)!
    Tuple<C,0u,Cflag>& removeAt(uint32_t position); // O(n)!

    C mkSum()const{C fout; if (tup_size == 0) ExOp::toZero(fout); else{fout = data[0]; for(int i=1;i<tup_size;i++) fout +=  data[i];} return fout;}
    C mkProduct()const{C fout; if (tup_size == 0) ExOp::toOne(fout); else{fout = data[0]; for(int i=1;i<tup_size;i++) fout *=  data[i];} return fout;}

	//void push_back(const C&);
    //C pop_back();
	double pnorm() const;
	double norm() const;

    Tuple<INDEX_TYPE> mkOrdering()const;
    Tuple<C,0u,Cflag> mkInverse()const;
    Tuple<C,0u,Cflag> mkPowInvInt(int a)const;
    Tuple<C,0u,Cflag> mkPowInt(int a)const;
    Tuple<C,0u,Cflag>& toInverse();
    void wrDeterminant(C&)const;
    template <class O> void wrDeterminant(O&) const;
    void wrTrace(C&)const;
    template <class O> void wrTrace(O&) const;

#undef LFHTEMP
#define LFHTEMP template<class O,unsigned int oTSIZE, Tuple_flag Oflag>
	LFHTEMP char compare(const Tuple<O,oTSIZE, Oflag>& other) const;
	LFHTEMP bool operator>(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 1);}
	LFHTEMP bool operator<(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) == 2);}
	LFHTEMP bool operator>=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 2);}
	LFHTEMP bool operator<=(const Tuple<O,oTSIZE, Oflag>& other) const{return(compare(other) != 1);}
	LFHTEMP bool operator==(const Tuple<O,oTSIZE, Oflag>& other) const{return (tup_size == other.getSize()) ? (compare(other) == 0) : false;}
	LFHTEMP bool operator!=(const Tuple<O,oTSIZE, Oflag>& other) const{return (tup_size == other.getSize()) ? (compare(other) != 0) : true;}

    template<class O> std::vector<O>& wrStdVector(std::vector<O>& fout)const;
    void show(char* &ptr, int level) const;
//     string type_tostring()const;
//     string type_tostring()const;
};
template<class C>
class Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>{
public:
    unsigned int tup_size;
    C* data;
    unsigned int getSize() const{return tup_size;}
    void setSize(unsigned int s){tup_size = s;}
    class Iterator{
    public:
        uint32_t cur, maxval;
        Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& target;
        Iterator(Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& trg): target(trg), cur(0),maxval(trg.getSize()){}
        Iterator(Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& trg, uint32_t _start, uint32_t _endoffset): target(trg), cur(_start),maxval(_endoffset){}
        operator LFHPrimitive::Iterator<uint32_t, C>(){using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, C>(std::bind( &Iterator::next , this, _1,_2), &cur, (cur != maxval) ? &target[cur] :NULL );}
        operator bool (){cur =0; maxval = target.getSize(); return (cur != maxval);}
        bool operator++(int){return (++cur != maxval);}
        uint32_t operator()()const{return cur;}
        uint32_t getOffset()const{return cur;}
        C* operator->(){return target.data + cur;}
        C& operator*(){return target[cur];}
        bool next(uint32_t*& keyt, const C*& itar){if (++cur == maxval) return false; itar = &(target[cur]); return true;}
	};
    class ConstIterator{
    public:
        uint32_t cur, maxval;
        const Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& target;
        ConstIterator(const Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& trg): target(trg){}
        operator bool (){cur =0; maxval = target.getSize(); return (cur != maxval);}
        bool operator++(int){return (++cur != maxval);}
        uint32_t operator()()const{return cur;}
        uint32_t getOffset()const{return cur;}
        const C* operator->(){return target.data + cur;}
        const C& operator*(){return target[cur];}

        operator LFHPrimitive::Iterator<uint32_t, const C>(){using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, const C>(std::bind( &ConstIterator::next , this, _1,_2), &cur, (cur != maxval) ? &target[cur] :NULL );}
        bool next(uint32_t*& keyt, const C*& itar){if (++cur == maxval) return false; itar = &(target[cur]); return true;}

	};
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::Iterator getIterator(){return Iterator(*this);}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::ConstIterator getIterator()const{return ConstIterator(*this);}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::Iterator getPartitionIterator(int which, int nbparts){return Iterator(*this, (which * getSize()) /nbparts , ((which+1) * getSize()) / nbparts );}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>::ConstIterator getPartitionIterator(int which, int nbparts)const{return ConstIterator(*this, (which * getSize()) /nbparts , ((which+1) * getSize()) / nbparts  );}
    operator LFHPrimitive::Iterator<uint32_t, const C> ()const {return (LFHPrimitive::Iterator<uint32_t, const C>) ConstIterator(*this);}


    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator=(C& target) {data = &target;return *this;}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator=(C* target) {data = target; return *this;}

    template<unsigned int OSIZE, Tuple_flag OFLAG> Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator=(const Tuple<C,OSIZE,OFLAG> &target){for(unsigned int i=0;i<tup_size;i++) data[i] = target[i];  return(*this);}

    bool isValid()const{for(unsigned int i =0; i < tup_size;i++) if (!ExOp::isValid(data[i])) return false; return true;}

    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& toZero(){for(unsigned int i=0;i<tup_size;i++) ExOp::toZero(data[i]); return(*this);}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& toOne(){for(unsigned int i=0;i<tup_size;i++) ExOp::toOne(data[i]); return(*this);}
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& toRand(){for(unsigned int i=0;i<tup_size;i++) ExOp::toRand(data[i]); return(*this);}
	Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& toMemfree(){if (tup_size != 0) delete[](data); tup_size=0; return(*this);}

	Tuple<C,0u> operator-(LFHPrimitive::Iterator<uint32_t, const C> ite)const{Tuple<C,0u> fout; fout.setSize(tup_size); if (ite()) do{fout[ite()] = data[ite()] - *ite; }while(ite++); return fout;}

    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator+=(LFHPrimitive::Iterator<uint32_t, const C> other);
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator-=(LFHPrimitive::Iterator<uint32_t, const C> other);
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator*=(LFHPrimitive::Iterator<uint32_t, const C> other);
    Tuple<C,0u, TUPLE_FLAG_REMOTE_MEMORY>& operator/=(LFHPrimitive::Iterator<uint32_t, const C> other);



    const C& operator[](unsigned int index)const{return data[index];}
    C& operator[](unsigned int index){return data[index];}
    C& last(){return data[tup_size-1];}
    const C& last()const{return data[tup_size-1];}


    void show(FILE* f = stdout, int level=0)const;
};
template<Tuple_flag Cflag>
class Tuple<Anything,0u, Cflag>{
public:
	uint32_t tup_size;
	uint32_t anytype;
	uint8_t* data;

	Tuple<Anything,0u, Cflag>& setSize(uint32_t nsize);

    void show(FILE* f = stdout, int level= 0)const;
	ERRCODE load(FILE *f);
	ERRCODE save(FILE *f) const;
};

template<class C>
class Tensor{
public:
    Tuple<C, 0u> data;
    Tuple<uint32_t, 0u> shape;

	~Tensor(){}
	Tensor& setShape(const Tuple<uint32_t> &_shape){shape = _shape; data.setSize(shape.mkProduct()); return *this;}

	Tensor& toZero(){data.toZero();return(*this);}
	Tensor& toOne(){data.toOne();return(*this);}


	ERRCODE load(FILE *f){ERRCODE fout = shape.load(f); fout |= data.load(f); return fout;}
	ERRCODE save(FILE *f) const{ERRCODE fout = shape.save(f); fout |= data.save(f); return fout;}
};

template<class C, unsigned int NBDIMS>
class RemoteMemory{
public:
	Tuple<uint32_t, NBDIMS> dims;
	Tuple<int32_t, NBDIMS> offsets; // allows flipped direction

	C* remote_data;

	class Iterator{
		public:
		typedef typename MetaType<uint32_t ,std::is_const<C>::value >::IS_CONST KEYITERATOR_TYPE;
		Tuple<uint32_t, NBDIMS> coor;


		RemoteMemory<C, NBDIMS> &target;
		uint32_t getSize(){uint32_t fout =1; for(int i=0;i< NBDIMS;i++) fout *= target.dims[i]; return fout;}
		C* cur;
        Iterator(RemoteMemory<C, NBDIMS> &trg): target(trg){}
        operator bool (){for(unsigned int i=0;i<target.dims.getSize();i++) if (target.dims[i] == 0) return false; coor.toZero(); cur = target.remote_data; return true;}
        bool operator++(int){uint32_t i; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i]) {cur += target.offsets[i]; break;} cur -= target.offsets[i] * (coor[i] -1); coor[i]=0;} return (i < coor.getSize());}

        const Tuple<uint32_t, NBDIMS>& operator()()const{return coor;}
		C* operator->(){return cur;}
        C& operator*(){return *cur;}
		Iterator& mkIterator(){return *this;}
	};

	/*class ArrayIterator{
		public:
		typedef typename MetaType<uint32_t ,std::is_const<C>::value >::IS_CONST KEYITERATOR_TYPE;
		Tuple<uint32_t, NBDIMS> coor;

		RemoteMemory<C, NBDIMS> &target;
		C* cur;
        Iterator(RemoteMemory<C, NBDIMS> &trg): target(trg){}
        operator bool (){for(unsigned int i=0;i<target.dims.getSize();i++) if (target.dims[i] == 0) return false; coor.toZero(); cur = target.remote_data; return true;}
        bool operator++(int){uint32_t i; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i]) {cur += target.offsets[i]; break;} cur -= target.offsets[i] * (coor[i] -1); coor[i]=0;} return (i < coor.getSize());}

        const Tuple<uint32_t, NBDIMS>& operator()()const{return coor;}
		C* operator->(){return cur;}
        C& operator*(){return *cur;}
		Iterator& mkIterator(){return *this;}
	};*/


	class PixelIterator{
		public:
		Tuple<uint32_t, NBDIMS - 1> coor;
		RemoteMemory<C, NBDIMS> &target;
		C* cur;
        PixelIterator(RemoteMemory<C, NBDIMS> &trg): target(trg){}
        uint32_t getSize(){uint32_t fout =1; for(int i=0;i< NBDIMS - 1;i++) fout *= target.dims[i+1]; return fout;}
        operator bool (){if (target.offsets[0] != 1) return false; for(unsigned int i=0;i<NBDIMS;i++) if (target.dims[i] == 0) return false; coor.toZero(); cur = target.remote_data; return true;}
        bool operator++(int){uint32_t i; for(i=0;i<coor.getSize();i++) {if ((++coor[i]) < target.dims[i+1]) {cur += target.offsets[i+1]; break;} cur -= target.offsets[i+1] * (coor[i] -1); coor[i]=0;} return (i < coor.getSize());}
        Tuple<uint32_t, NBDIMS - 1> operator()()const{return coor;}
		uint32_t getShape(){return target.dims[0];}

		C** operator->(){return &cur;}
        C* operator*(){return cur;}
		PixelIterator& mkIterator(){return *this;}
	};

	C& operator[](const Tuple<uint32_t, NBDIMS>& coor){
		uint32_t pos = coor[0] * offsets[0];
		for(unsigned int i=1;i < dims.getSize();i++) pos += coor[i] * offsets[i];
	return*(remote_data + pos);}

	C& operator[](const uint32_t* coor){
		uint32_t pos = coor[0] * offsets[0];
		for(unsigned int i=1;i < dims.getSize();i++) pos += coor[i] * offsets[i];
	return*(remote_data + pos);}


	Tuple<uint32_t, NBDIMS> getShape()const {return dims;}
	uint32_t getShape(int i)const{return dims[i];}

	typename RemoteMemory<C, NBDIMS>::Iterator mkIterator() const {return RemoteMemory<C, NBDIMS>::Iterator(*const_cast<RemoteMemory<C, NBDIMS> *> (this));}
	RemoteMemory<C, NBDIMS>::Iterator operator()(){return RemoteMemory<C, NBDIMS>::Iterator(*this);}
	RemoteMemory<C, NBDIMS>::PixelIterator mkVectorIterator(){if (offsets[0] != 1) printf("warning, not contiguous on last level...\n"); return RemoteMemory<C, NBDIMS>::PixelIterator(*this);}

	operator RemoteMemory<const C, NBDIMS>(){RemoteMemory<const C, NBDIMS> fout; fout.dims = dims; fout.offsets = offsets; fout.remote_data = remote_data; return fout;}

	RemoteMemory<C, NBDIMS>& toSwapDims(uint32_t x, uint32_t y){if ((x >= NBDIMS)||(y >= NBDIMS)) {printf("warning, swaping unexisting dimention!\n");return *this;} uint32_t tmp = dims[x]; dims[x] = dims[y]; dims[y] = tmp; tmp = offsets[x]; offsets[x] = offsets[y]; offsets[y] = tmp; return *this;}

	RemoteMemory<C, NBDIMS-1> operator[](int i);
	RemoteMemory<C, NBDIMS-1> operator[](unsigned int i);

	RemoteMemory<C, NBDIMS-1> selectSlice(int index, uint32_t which_dim);
	RemoteMemory<C, NBDIMS-1> selectSlice(uint32_t index, uint32_t which_dim);

	RemoteMemory<C, NBDIMS-1> vectorizeOntoDimension(uint32_t dim);
	RemoteMemory<C, NBDIMS+1> splitDimension(uint32_t dim, uint32_t length);



	TMatrix<C> eigenVectorOfInnerProduct(int direction, int nbeigen, Tuple<double> *eigenvalues)const;


	const RemoteMemory<C, NBDIMS>& show(FILE* f= stdout, int level =0) const{
		if (remote_data == NULL) {fprintf(f, "Empty/invalid remote %iD memory\n", NBDIMS); return *this;}
		unsigned int i,j,k;
		switch(NBDIMS){
		case 2:
			fprintf(f, "Remote 2D memory %i x %i\n", dims[0], dims[1]);
			for(i=0;i<dims[0];i++){
				for(j=0;j<dims[1];j++){
					ExOp::show(remote_data[i*offsets[0] + j*offsets[1]],f,level+2);
					fprintf(f,"%c", (j == (dims[1]-1)) ? '\n' : '\t');
				}
			}
		break; case 3:
			fprintf(f, "Remote 3D memory %i x %i x %i\n", dims[0], dims[1], dims[2]);
			for(k=0;k<dims[2];k++){
				fprintf(f,"Slice (%i/%i):\n", k, dims[2]);
				for(i=0;i<dims[0];i++){
					for(j=0;j<dims[1];j++){
						ExOp::show(remote_data[i*offsets[0] + j*offsets[1]+ k*offsets[2]],f,level+2);
						fprintf(f,"%c", (j == (dims[1]-1)) ? '\n' : '\t');
					}
				}
			}
			fprintf(f,"End of Remote 3D memory %i x %i x %i\n", dims[0], dims[1], dims[2]);
		break; default: fprintf(f, "Remote memory %i dims\n", NBDIMS);
		}
	return *this;}
};


/*
	uint32_t i, j;
	i = NBDIMS -1;
	j = offsets.getSize();
	uint32_t pos = coor[i];
	for(i--; j > 0;i--) pos = pos * offsets[--j] + coor[i];
	for(; i !=0xFFFFFFFF;i--) pos = pos * dims[i] + coor[i];
*/
template<class C>
class RemoteMemory<C, 1u>{
public:
	typedef uint32_t KEYITERATORMKR_TYPE;
	Tuple<uint32_t, 1u> dims;
	Tuple<uint32_t, 1u> offsets;
	C* remote_data;
	RemoteMemory()=default;
	RemoteMemory(C* _daptr, uint32_t _dims, uint32_t _offsets = 1u) : remote_data(_daptr){dims[0] = _dims; offsets[0] = _offsets;}
	C& operator[](uint32_t coor) const{return *(remote_data + coor * offsets[0]);}

	class Iterator{
		public:
		typedef typename MetaType<uint32_t , std::is_const<C>::value>::IS_CONST KEYITERATOR_TYPE;
		uint32_t ite;
		RemoteMemory<C, 1u> &target;
		C* cur;
        Iterator(RemoteMemory<C, 1u> &trg): target(trg){}
        operator bool (){if (target.dims[0] == 0) return false; cur = target.remote_data; ite =0; return true;}
        bool operator++(int){cur += target.offsets[0]; return ((++ite) < target.dims[0]);}
        uint32_t operator()()const{return ite;}

		C* operator->(){return cur;}
        C& operator*(){return *cur;}
		Iterator& mkIterator(){return *this;}
	};
	RemoteMemory<C, 1u>::Iterator mkIterator(){return RemoteMemory<C, 1u>::Iterator(*this);}

	C& operator[](int i){return *(remote_data + i * offsets[0]);}
	const C& operator[](int i) const {return *(remote_data + i * offsets[0]);}

	Tuple<uint32_t> getShape()const {Tuple<uint32_t> shape; shape.setSize(1); shape[0] = dims[0]; return shape;}
	Tuple<uint32_t, 1u> getShape(int i)const {Tuple<uint32_t, 1u> fout; fout[0] = dims[0]; return fout;} // i should be 0 right :D
	template<class I, ENABLEIF_ITERATOR(I, uint32_t) =true> auto mkInnerProd(I x) const -> decltype(remote_data[0] * (*(x.mkIterator())));
	template<class I, ENABLEIF_ITERATOR(I, uint32_t) =true> RemoteMemory<C, 1u>& operator+=(I x);
	template<class I, ENABLEIF_ITERATOR(I, uint32_t) =true> RemoteMemory<C, 1u>& operator-=(I x);

	const RemoteMemory<C, 1u>& show(FILE*f =stdout, int level =0) const;
	RemoteMemory<C, 1u>& show(FILE*f =stdout, int level =0) {return const_cast<RemoteMemory<C, 1u>&>(const_cast<const RemoteMemory<C, 1u>*>(this)->show(f,level));}
};

template<class C>
class RemoteMemory<C, 0u>{
public:
	Tuple<uint32_t, 0u> dims;
	Tuple<uint32_t, 0u> offsets;
	C* remote_data;
};




/*
template<class C, int size> struct IteratorScope< C, Tuple<C,size> > : public AbstractIterator< C, Tuple<C,size> > {
	enum {valid = true};
	unsigned int i;
	IteratorScope(): i(0){}
	void init(Tuple<C,size> & ob) {i=0;}
	C* next(Tuple<C,size> & ob){ i++; return((i <= size)? ob.data + i - 1 : NULL);}
};

template<class C, class D, int size> struct IteratorScope<D, Tuple<C,size> > : public AbstractIterator<D, Tuple<C,size> >{
	enum {valid = true};
	unsigned int i;
	IteratorScope<D, C> j;
	IteratorScope(Tuple<C,size> & ob){}
	void init(Tuple<C,size> & ob) {i=0; j.init();}

	C* next(Tuple<C,size> & ob){
		D* _out = j.next(ob.data[i]);
		while (_out == NULL){
			i++;
		if (i == size) return(NULL);
		j.init();
		_out = j.next(ob.data[i]);
		}
		return(_out);
		}
};*/







template<class C, int size,Tuple_flag Cflag> struct isTypeEquivalent< C*, Tuple<C,size, Cflag> > {enum {ans = true }; };
template<class C, int size,Tuple_flag Cflag> struct isTypeEquivalent< Tuple<C,size, Cflag> , C*> {enum {ans = true }; };

// data order is row major order! for now
template<class C, unsigned int sizex, unsigned int sizey, Tuple_flag TF>
class TMatrix{
	void leftHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint); // bottom of matrix
	void rightHouseHolderMultiply(const double * const vec, const int lenght, const double& sqrt_den, bool hint);
public:
	const TMatrix<C,sizex,sizey,TF>& LeftHouseHolderMultiply(const double * const vec, const double& sqrt_den, const int lenght, C* inner, bool hint); // top of matrix
    typedef YESNO<false> IS_COMMUTATIVE;
	// ExOp section
	typedef std::integral_constant<bool, ExCo<C>::IsPOD::value > IsPOD;
	// ExOp section
	Tuple<C, sizex*sizey, TF> data;
	bool fake_inverse;

	TMatrix(){}
	TMatrix(C* idata){memcpy(data,idata,sizeof(C)*sizex*sizey);}
	TMatrix(TMatrix<C, sizex, sizey,TF> const & idata);
	TMatrix(DataGrid<C, 2> const & idata);

	void zero(){unsigned int i = sizex*sizey;for(i--;i != ExCo<unsigned int>::mkMaximum() ; i--) ExOp::toZero(data[i]);}
	TMatrix<C,sizex,sizey,TF>& toZero();
	TMatrix<C,sizex,sizey,TF>& toOne();
	TMatrix<C,sizex,sizey,TF>& toRand();
    TMatrix<C,sizex,sizey,TF>& toRandSymetric();

	Trianglix<C,sizey> operator*(const Trianglix<C,sizex> &input) const;
    uint32_t nbRows()const{return sizex;}
    uint32_t nbCols()const{return sizey;}

	C operator()(const Tuple<unsigned int, 2> &coor) const{return data[coor[0] + coor[1] * sizex];}
	C& operator()(const Tuple<unsigned int, 2> &coor){return data[coor[0] + coor[1] * sizex];}
    C operator()(unsigned int x, unsigned int y) const{return data[x + y * sizex];}
    C& operator()(unsigned int x, unsigned int y){return data[x + y * sizex];}

	void getDims(Tuple<unsigned int, 2> &o_dims) const{o_dims[0] = sizex;o_dims[1] = sizey; }

	void bidiag(); // test

	bool isValid() const{
		int i;
		for(i=0;i<sizex*sizey;i++) if (!ExOp::isValid(data[i])) break;
		return (i == sizex*sizey);
	}

//	template<class O> Tuple<C,size>& operator+=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator-=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator*=(Tuple<O,size> const & other);
//	template<class O> Tuple<C,size>& operator/=(Tuple<O,size> const & other);

	TMatrix<C,sizex,sizey,TF>& operator=(TMatrix<C, sizex, sizey, TF> const & idata);

	void operator()(Tuple<C, sizey> &, Tuple<C, sizex> &) const;
	void derivative(Tuple<C, sizey> &, Tuple<C, sizex > &, int in_direct) const;
	void derivativeMatrix(TMatrix<C, sizex, sizey, TF> &, Tuple<C, sizex > &) const;

	template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& operator+=(TMatrix<O,sizex,sizey,OTF> const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] += other.data[i];return(*this);}
	template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& operator-=(TMatrix<O,sizex,sizey,OTF> const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] -= other.data[i];return(*this);}
	template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& operator*=(TMatrix<O,sizex,sizex,OTF> const & other);
	template<class O, Tuple_flag OTF> TMatrix<C,sizex,sizey,TF>& operator/=(TMatrix<O,sizex,sizex,OTF> const & other);

	template<class O, unsigned int osize, Tuple_flag OTF> TMatrix<C,osize,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator*(TMatrix<O,osize,sizex,OTF> const & other) const;
	Tuple<C,sizey> operator*(Tuple<C,sizex> const & other) const;
	template<class O> auto operator*(Tuple<O,sizex> const & other)const -> Tuple<decltype(data.data[0] *other.data.data[0]),sizey>;
	template<class O> auto mkBackMult(Tuple<O,sizey> const & other)const -> Tuple<decltype(data.data[0] *other.data.data[0]),sizex>;


	// template<class O, class A, Tuple_flag FA, Tuple_flag FO> Tuple<O, sizex,FO> operator*(Tuple<A,sizey, FA> const & in);


template<class O> TMatrix<C,sizex,sizey,TF>& operator+=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] += other;return(*this);}
template<class O> TMatrix<C,sizex,sizey,TF>& operator-=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] -= other;return(*this);}
template<class O> TMatrix<C,sizex,sizey,TF>& operator*=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] *= other;return(*this);}
template<class O> TMatrix<C,sizex,sizey,TF>& operator/=(O const & other){for(unsigned int i=sizex*sizey-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] /= other;return(*this);}

template<class A> TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator+(A const & other) const {return(TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) += other);}
template<class A> TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator-(A const & other) const {return(TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) -= other);}
template<class A> TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator*(A const & other) const {return(TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) *= other);}
template<class A> TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator/(A const & other) const {return(TMatrix<C,sizex,sizey,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) /= other);}
/*
template<class O, class A> TMatrix<O,sizex,sizey> operator+(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator-(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator*(A const & other) const;
template<class O, class A> TMatrix<O,sizex,sizey> operator/(A const & other) const;*/

	template<class O,Tuple_flag FO> TMatrix<C,sizex,sizey,TF> scale_rows(Tuple<O, sizey, FO> const &) const;
	template<class O,Tuple_flag FO> TMatrix<C,sizex,sizey,TF> scale_cols(Tuple<O, sizex, FO> const &) const;


//template<class O, unsigned int osize> TMatrix<C,sizex,sizey> operator-(O const & other) const{return( (TMatrix(*this)) -= other );}

//	template<class O> TMatrix<C,sizex,sizey> operator*(O const & other) const{return( (TMatrix(*this)) *= other );}
//	template<class O> TMatrix<C,sizex,sizey> operator/(O const & other) const{return( (TMatrix(*this)) /= other );}


	TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkInverse() const;
	TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkTranspose() const;
    TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkTrju() const;



//	PolyThing<C,0> makeCharacteristic() const;
	Continuous<C,sizex,C,sizey>* derivative(Tuple<C,sizex> const & where) const; // a TMatrix is it's own derivative!
	void show(FILE* o = stdout, int level=0) const;

    TMatrix<C,sizey,sizex,TUPLE_FLAG_CONSTRUCT_MASK(TF)> inverse_2() const;

	C determinant() const;

};
// data order is column major order!
template<class C, Tuple_flag TF>
class TMatrix<C, 0u, 0u, TF> { //  : public ConstGrid<C,2u>
public:
    uint32_t sizes[2];
    Tuple<C, 0u, TF> data;

    class Iterator{
    public:
    	typedef Tuple<uint32_t, 2u> KEYITERATOR_TYPE;
        TMatrix<C, 0u, 0u,TF>& target;
        Tuple<uint32_t, 2u> coor;
        C* cur;
        C* maxval;
        Iterator(TMatrix<C,0u,0u,TF>& trg): target(trg){}
        uint32_t getSize(){return target.data.tup_size;}
        operator bool (){cur = target.data(); maxval = cur + target.sizes[0]*target.sizes[1]; coor.toZero(); return (maxval!=cur);}
        bool operator++(int){uint32_t i; cur++; for(i=0;i<2u;i++) {if ((++coor[i]) < target.sizes[i]) break; coor[i]=0;} return (i < 2u);}
        uint32_t operator[](int dim)const{return coor[dim];}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        C* operator->(){return cur;}
        C& operator*(){return *cur;}
        Iterator& mkIterator(){return *this;}
	};
    class ConstIterator{
    public:
    	typedef const Tuple<uint32_t, 2u> KEYITERATOR_TYPE;
        Tuple<uint32_t, 2u> coor;
        const C* cur;
        const C* maxval;
        const TMatrix<C, 0u, 0u,TF>& target;
        ConstIterator(const TMatrix<C,0u,0u,TF>& trg): target(trg){}

        uint32_t operator[](int dim)const{return coor[dim];}
        operator bool (){cur = target.data(); maxval = cur + target.sizes[0]*target.sizes[1]; coor.toZero(); return (maxval!=cur);}
        bool operator++(int){uint32_t i; cur++; for(i=0;i<2u;i++) {if ((++coor[i]) < target.sizes[i]) break; coor[i]=0;} return (i < 2u);}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        const C* operator->(){return cur;}
        C operator*(){return *cur;}
        ConstIterator& mkIterator(){return *this;}
	};
    TMatrix<C, 0u, 0u,TF>::Iterator getIterator(){return TMatrix<C, 0u, 0u,TF>::Iterator(*this);};
    TMatrix<C, 0u, 0u,TF>::ConstIterator getIterator()const{return TMatrix<C, 0u, 0u,TF>::ConstIterator(*this);};
    TMatrix<C, 0u, 0u,TF>::Iterator mkIterator(){return TMatrix<C, 0u, 0u,TF>::Iterator(*this);};
    TMatrix<C, 0u, 0u,TF>::ConstIterator mkIterator()const{return TMatrix<C, 0u, 0u,TF>::ConstIterator(*this);};
    class ColumnIteratorOLD{
    public:
        uint32_t col, endc;
        Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> cur;
        TMatrix<C, 0u, 0u,TF>& target;
        ColumnIteratorOLD(TMatrix<C,0u,0u,TF>& trg): target(trg),col(0), endc(trg.sizes[1]){cur.data = trg.data();}
        ColumnIteratorOLD(TMatrix<C,0u,0u,TF>& trg, uint32_t _col, uint32_t _endc): target(trg),col(_col), endc(_endc){cur.data = trg.data() + _col * trg.sizes[0];}

        operator LFHPrimitive::Iterator<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> > (){using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> >(std::bind( &ColumnIteratorOLD::next , this, _1,_2), &col, endc != col ? &cur :NULL);} // ,
        bool next(uint32_t*& keyt, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>*& itar){cur.data += target.sizes[0]; return(endc!=(++col));}
        operator bool (){return (target.sizes[1]!=col);}
        bool operator++(int){cur.data += target.sizes[0];return(endc!=(++col));}
        uint32_t operator()()const{return col;}
        Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>* operator->(){return &cur;}
        Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>& operator*(){return cur;}
	};
    class ColumnConstIteratorOLD{
    public:
        uint32_t col, endc;
        Tuple<const C, 0u, TUPLE_FLAG_REMOTE_MEMORY> cur;
        const TMatrix<C, 0u, 0u,TF>& target;
        ColumnConstIteratorOLD(const TMatrix<C,0u,0u,TF>& trg): target(trg),col(0), endc(trg.sizes[1]){cur.data = trg.data();}
        ColumnConstIteratorOLD(TMatrix<C,0u,0u,TF>& trg, uint32_t _col, uint32_t _endc): target(trg),col(_col), endc(_endc){cur.data = trg.data() + _col * trg.sizes[0];}
        operator bool (){return (target.sizes[1]!= (col));}
        bool operator++(int){cur.data += target.sizes[0];return (target.sizes[1]!= (++col));}
        uint32_t operator()()const{return col;}
        const Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>* operator->(){return &cur;}
        const Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>& operator*(){return cur;}
	};

	template<int ISCONST>
	class ColumnSelector{
	public:
		typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
		uint32_t col, row;
		typename MetaType<TMatrix<C, 0u, 0u,TF> ,ISCONST>::IS_CONST_REF target;
		typename MetaType<C,ISCONST>::IS_CONST_PTR cur;
        ColumnSelector(typename MetaType<TMatrix<C, 0u, 0u,TF> ,ISCONST>::IS_CONST_REF trg, uint32_t _col): target(trg),col(_col){}
        operator bool (){if (target.data.tup_size == 0) return false; cur = target.data() + col * target.sizes[0]; row =0; return true;}
        bool operator++(int){cur++; return ((++row) < target.sizes[0]);}
        uint32_t operator()()const{return row;}

		typename MetaType<C,ISCONST>::IS_CONST_PTR operator->(){return cur;}
        typename MetaType<C,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return *cur;}
        ColumnSelector<ISCONST>& mkIterator(){return *this;}
	};
	template<int ISCONST>
	class RowSelector{
	public:
		typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
		uint32_t col, row;
		typename MetaType<TMatrix<C, 0u, 0u,TF> ,ISCONST>::IS_CONST_REF target;
		typename MetaType<C,ISCONST>::IS_CONST_PTR cur;
        RowSelector(typename MetaType<TMatrix<C, 0u, 0u,TF> ,ISCONST>::IS_CONST_REF trg, uint32_t _row): target(trg),row(_row){}
        operator bool (){if (target.data.tup_size == 0) return false; cur = target.data() + row; col =0; return true;}
        bool operator++(int){cur += target.sizes[0]; return ((++col) < target.sizes[1]);}
        uint32_t operator()()const{return col;}

		typename MetaType<C,ISCONST>::IS_CONST_PTR operator->(){return cur;}
        typename MetaType<C,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return *cur;}
        RowSelector<ISCONST>& mkIterator(){return *this;}
	};


	TMatrix<C, 0u, 0u,TF> cmpEigenVectorOfInnerProduct(Tuple<C> &eigenvalues, int nbeigen =0)const;
	TMatrix<C, 0u, 0u,TF>& leftHouseHolderMultiply(RemoteMemory<const C,1> vec, const double sqrt_den, uint32_t start = 0u, bool applytoall = false);
	TMatrix<C, 0u, 0u,TF>& rightHouseHolderMultiply(RemoteMemory<const C,1> vec, const double sqrt_den, uint32_t start = 0u, bool applytoall = false);
	TMatrix<C, 0u, 0u,TF>& rightHouseHolderMultiply_inplace(const double sqrt_den, uint32_t start = 0u);

	TMatrix<C, 0u, 0u,TF>& rotColPair(uint32_t i, double cos, double sin);
	TMatrix<C, 0u, 0u,TF>& negateColumn(uint32_t i);
	TMatrix<C, 0u, 0u,TF>& rotRowPair(uint32_t i, double cos, double sin);
	TMatrix<C, 0u, 0u,TF>& negateRow(uint32_t i);
	TMatrix<C, 0u, 0u,TF>& toFlipColumns();
	TMatrix<C, 0u, 0u,TF>& toFlipRows();

	TMatrix<C, 0u, 0u,TF>& toUpperTriangular(bool make_squarre = true);
	TMatrix<C, 0u, 0u,TF>& toUpperTriangulariser(RemoteMemory<const C,2> matrix); // Find Holder transformation such that X^t * (matrix) is upper triangular


	// if nbrows <= nbcols, returns U matrix and orthogonal is V^t, otherwise return V^t and orthogonal is U
	Tuple<C> mkSingualDecomposition(TMatrix<C> &U, TMatrix<C> &Vt)const;
	TMatrix<C, 0u, 0u,TF>& toSingularDecomposition(RemoteMemory<const C,2> matrix, Tuple<C> &eigenvalues, TMatrix<C> &orthogonal); // PARTLY WORKING ONLY TODO Find Holder transformation such that X * this is upper triangular

	//TMatrix<C, 0u, 0u,TF>& toUpperMatrixOf()(RemoteMemory<C, 2u> fout);


	TMatrix<C, 0u, 0u,TF>::ColumnSelector<0> selectorColumn(uint32_t col){if (col >= sizes[1]) {printf("Selected an illegal column %i >= %i\n", col, sizes[0]); exit(1);} return TMatrix<C, 0u, 0u,TF>::ColumnSelector<0>(*this, col);};
    TMatrix<C, 0u, 0u,TF>::ColumnSelector<1> selectorColumn(uint32_t col)const{if (col >= sizes[1]) {printf("Selected an illegal column %i >= %i\n", col, sizes[0]); exit(1);} return TMatrix<C, 0u, 0u,TF>::ColumnSelector<1>(*this, col);};
	TMatrix<C, 0u, 0u,TF>::RowSelector<0> selectorRow(uint32_t row){if (row >= sizes[0]) {printf("Selected an illegal row %i >= %i\n", row, sizes[1]); exit(1);} return TMatrix<C, 0u, 0u,TF>::RowSelector<0>(*this, row);};
    TMatrix<C, 0u, 0u,TF>::RowSelector<1> selectorRow(uint32_t row)const{if (row >= sizes[0]) {printf("Selected an illegal row %i >= %i\n", row, sizes[1]); exit(1);} return TMatrix<C, 0u, 0u,TF>::RowSelector<1>(*this, row);};

	operator RemoteMemory<C, 2u>(){RemoteMemory<C, 2u> fout; fout.remote_data = data(); fout.dims[0] = sizes[0];fout.dims[1] = sizes[1]; fout.offsets[0]=1; fout.offsets[1] = sizes[0]; return fout;}
	operator RemoteMemory<const C, 2u>()const{RemoteMemory<const C, 2u> fout; fout.remote_data = data(); fout.dims[0] = sizes[0];fout.dims[1] = sizes[1]; fout.offsets[0]=1; fout.offsets[1] = sizes[0]; return fout;}

	RemoteMemory<C, 2u> operator*(){RemoteMemory<C, 2u> fout; fout.remote_data = data(); fout.dims[0] = sizes[0];fout.dims[1] = sizes[1]; fout.offsets[0]=1; fout.offsets[1] = sizes[0]; return fout;}
	RemoteMemory<const C, 2u> operator*()const{RemoteMemory<const C, 2u> fout; fout.remote_data = data(); fout.dims[0] = sizes[0];fout.dims[1] = sizes[1]; fout.offsets[0]=1; fout.offsets[1] = sizes[0]; return fout;}


	RemoteMemory<C, 1u> selectColumn(uint32_t col){if (col >= sizes[1]) {printf("Selected an illegal column %i >= %i\n", col, sizes[0]); exit(1);}return RemoteMemory<C, 1u>(data() + col * sizes[0], (data.tup_size == 0) ? 0 : sizes[0],1);}
    RemoteMemory<const C, 1u> selectColumn(uint32_t col)const{if (col >= sizes[1]) {printf("Selected an illegal column %i >= %i\n", col, sizes[0]); exit(1);}return RemoteMemory<const C, 1u>(data() + col * sizes[0], (data.tup_size == 0) ? 0 : sizes[0],1);}
	RemoteMemory<C, 1u> selectRow(uint32_t row){if (row >= sizes[0]) {printf("Selected an illegal row %i >= %i\n", row, sizes[1]); exit(1);}return RemoteMemory<C, 1u>(data() + row, (data.tup_size == 0) ? 0 : sizes[1],sizes[0]);}
    RemoteMemory<const C, 1u> selectRow(uint32_t row)const{if (row >= sizes[0]) {printf("Selected an illegal row %i >= %i\n", row, sizes[1]); exit(1);}return RemoteMemory<const C, 1u>(data() + row, (data.tup_size == 0) ? 0 : sizes[1],sizes[0]);}
	RemoteMemory<C, 1u> vectorize(){return RemoteMemory<C, 1u>(data(), (data.tup_size == 0) ? 0 : sizes[1] * sizes[0],1);}
    RemoteMemory<const C, 1u> vectorize()const{return RemoteMemory<const C, 1u>(data(), (data.tup_size == 0) ? 0 : sizes[1] * sizes[0],1);}



    TMatrix<C, 0u, 0u,TF>::ColumnIteratorOLD getColumnIterator(){return TMatrix<C, 0u, 0u,TF>::ColumnIteratorOLD(*this);};
    TMatrix<C, 0u, 0u,TF>::ColumnConstIteratorOLD getColumnIterator()const{return TMatrix<C, 0u, 0u,TF>::ColumnConstIteratorOLD(*this);};
    TMatrix<C, 0u, 0u,TF>::ColumnIteratorOLD getColumnParititionIterator(uint32_t i, uint32_t nb){return TMatrix<C, 0u, 0u,TF>::ColumnIteratorOLD(*this, (i * sizes[1]) / nb, ((i+1) * sizes[1]) / nb);};
    TMatrix<C, 0u, 0u,TF>::ColumnConstIteratorOLD getColumnParititionIterator(uint32_t i, uint32_t nb)const{return TMatrix<C, 0u, 0u,TF>::ColumnConstIteratorOLD(*this, (i * sizes[1]) / nb, ((i+1) * sizes[1]) / nb);};
    LFHPrimitive::Iterator<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> > getColumnPartitionPolyIterator(uint32_t i, uint32_t nb){return (LFHPrimitive::Iterator<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> >)this->getColumnParititionIterator(i,nb);}
//    LFHPrimitive::Iterator<uint32_t, Tuple<const C, 0u, TUPLE_FLAG_REMOTE_MEMORY> > getColumnPartitionPolyIterator(uint32_t i, uint32_t nb)const{return (LFHPrimitive::Iterator<uint32_t, Tuple<const C, 0u, TUPLE_FLAG_REMOTE_MEMORY> >)this->getColumnParititionIterator(i,nb);}


    //explicit operator Accessor<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> > (){using std::placeholders::_1; using std::placeholders::_2;  return LFHPrimitive::Accessor<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> >(std::bind( &accessCol, this, _1,_2));}
    //bool accessCol(const uint32_t& col, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY>*& target){}

    explicit operator IteratorMaker<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> > (){using std::placeholders::_1; using std::placeholders::_2; return IteratorMaker<uint32_t, Tuple<C, 0u, TUPLE_FLAG_REMOTE_MEMORY> >(std::bind( &getColumnPartitionPolyIterator, this, _1,_2));}


    TMatrix();
    ~TMatrix();
    TMatrix(const TMatrix<C, 0u, 0u,TF>& _clone);
    TMatrix(TMatrix<C, 0u, 0u, TF>&& );
    TMatrix(RemoteMemory<const C, 2u> _clone);
    template<class ARR, typename RequireArray<ARR>::type = true> TMatrix(ARR arr);

    // TMatrix(const SparseMatrix<C> & idata);
    TMatrix<C,0u,0u,TF>& operator=(const TMatrix<C, 0u, 0u,TF> &other);
    TMatrix<C,0u,0u,TF>& operator=(TMatrix<C, 0u, 0u,TF> && other);
	template<class ARR, typename RequireArray<ARR>::type = true> TMatrix<C, 0u, 0u,TF>& operator=(ARR arr);

    uint32_t nbRows()const{return sizes[0];}
    uint32_t nbCols()const{return sizes[1];}
	Tuple<C> sumRows() const;
    Tuple<C> sumCols() const;

	template<class D> TMatrix<C, 0u, 0u,TF>& scaleRows(const RemoteMemory<const D, 1u> o){
		if (o.dims[0] != sizes[0]) exit(1);
		C* cur = data(); for(int y=0;y<sizes[1];y++)
		for(int x=0;x<sizes[0];x++) *(cur++) *= o[x];
		return *this;}
	template<class D> TMatrix<C, 0u, 0u,TF>& scaleCols(const RemoteMemory<const D, 1u> o){if (o.dims[0] != sizes[1]) exit(1); C* cur = data(); for(int y=0;y<sizes[1];y++) for(int x=0;x<sizes[0];x++) *(cur++) *= o[y]; return *this;}

	C mkSum()const{return data.mkSum();}

	bool isValid() const{return data.isValid();}

    TMatrix<C, 0u, 0u,TF>& setSizes(unsigned int x, unsigned int y);
    TMatrix<C, 0u, 0u,TF>& setSizes(Tuple<uint32_t,2u> coor){return this->setSizes(coor[0],coor[1]);}

	TMatrix<C, 0u, 0u,TF>& toZero();
	TMatrix<C, 0u, 0u,TF>& toOne();
	TMatrix<C, 0u, 0u,TF>& toRand();




    C operator()(const Tuple<unsigned int, 2u> &coor) const;
    C& operator()(const Tuple<unsigned int, 2u> &coor);
    C operator()(unsigned int x, unsigned int y) const;
    C& operator()(unsigned int x, unsigned int y);

    TMatrix<C,0u,0u,TF>& permuteRows(const Tuple<uint32_t> &permutation);
    TMatrix<C,0u,0u,TF>& permuteCols(const Tuple<uint32_t> &permutation);

    Trianglix<C, 0u> mkInnerProduct() const;
    Trianglix<C, 0u> mkOuterProduct() const;

    template<class O> auto mkInnerMult(const Tuple<O> & left, const Tuple<O> & right)const -> decltype((this->data[0] * right[0]) * left[0]);
    template<class O> auto mkInnerMult(const SparseTuple<O> & left, const SparseTuple<O> & right)const -> decltype((this->data[0] * right[0]) * left[0]);


    void getDims(Tuple<unsigned int, 2u> &o_dims) const;
    TMatrix<C, 0u, 0u,TF>& toMemmove(TMatrix<C,0u,0u,TF>& other);
    Tuple<C,0u> getColumn(uint32_t colID)const;

    template<class O, Tuple_flag OTF> TMatrix<C, 0u, 0u,TF>& operator=(TMatrix<O, 0u, 0u,OTF> const & other);


   	template<class O> TMatrix< decltype(declval<C>() * declval<O>()) ,0u,0u,TF> mkLeftMult(RemoteMemory<const O,2u> & other) const;
	template<class O> TMatrix< decltype(declval<C>() * declval<O>()) ,0u,0u,TF> operator*(RemoteMemory<const O,2u> & other) const;
    template<class O, Tuple_flag OTF> TMatrix<C, 0u, 0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator*(TMatrix<O, 0u, 0u,OTF> const & other)const;
    template<class O> TMatrix<C, 0u, 0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator*(O const & other)const;

    void show(FILE* f = stdout, int level= 0)const;
	ERRCODE load(FILE *f);
	ERRCODE save(FILE *f) const;


    template<class O> TMatrix<C,0u,0u,TF>& operator+=(O const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] += other;return(*this);}
    template<class O> TMatrix<C,0u,0u,TF>& operator-=(O const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] -= other;return(*this);}
    template<class O, Tuple_flag OTF> TMatrix<C,0u,0u,TF>& operator+=(TMatrix<O,0u,0u,OTF> const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] += other.data[i];return(*this);}
    template<class O, Tuple_flag OTF> TMatrix<C,0u,0u,TF>& operator-=(TMatrix<O,0u,0u,OTF> const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] -= other.data[i];return(*this);}
    template<class O> TMatrix<C,0u,0u,TF>& operator*=(O const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] *= other;return(*this);}
    template<class O> TMatrix<C,0u,0u,TF>& operator/=(O const & other){for(unsigned int i=sizes[0]*sizes[1]-1;i!=ExCo<unsigned int>::mkMaximum();i--) data[i] /= other;return(*this);}
   // template<class O> TMatrix<C,0u,0u>& operator*=(TMatrix<O,0u,0u> const & other);
   // template<class O> TMatrix<C,0u,0u>& operator/=(TMatrix<O,0u,0u> const & other);

    template<class A> TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator+(A const & other) const {return(TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) += other);}
    template<class A> TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator-(A const & other) const {return(TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) -= other);}
    template<class A> TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> operator/(A const & other) const {return(TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)>(*this) *= other);}

    template<class O> TMatrix<C,0u,0u,TF>& operator+=(const SparseMatrix<O> & other);
    template<class O> TMatrix<C,0u,0u,TF>& operator-=(const SparseMatrix<O> & other);

    template<class O,class P> TMatrix<C,0u,0u,TF>& toAddOuterProd(const SparseTuple<O> & a, const SparseTuple<P> & b);
    template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& toAddOuterProd(const Tuple<O, S> & a, const SparseTuple<P> & b);
    template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& toAddOuterProd(const SparseTuple<O> & a, const Tuple<P, S> & b);
    template<class O,class P> TMatrix<C,0u,0u,TF>& toAddBackOuterProd(const SparseTuple<O> & a, const SparseTuple<P> & b);
    template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& toAddBackOuterProd(const Tuple<O, S> & a, const SparseTuple<P> & b);
    template<class O,class P, uint32_t S> TMatrix<C,0u,0u,TF>& toAddBackOuterProd(const SparseTuple<O> & a, const Tuple<P, S> & b);

    template<class O> TMatrix<C,0u,0u,TF>& toAddOuterProd(const SparseMatrix<O> & a, const SparseMatrix<O> & b);
    template<class O> TMatrix<C,0u,0u,TF>& toSubOuterProd(const SparseMatrix<O> & a, const SparseMatrix<O> & b);
    template<class O> TMatrix<C,0u,0u,TF>& toAddInnerProd(const SparseMatrix<O> & a, const SparseMatrix<O> & b);
    template<class O> TMatrix<C,0u,0u,TF>& toSubInnerProd(const SparseMatrix<O> & a, const SparseMatrix<O> & b);

    template<class O, Tuple_flag OFLAG> auto operator*(const Tuple<O,0u,OFLAG> & a)const -> Tuple<decltype(data.data[0]* a[0] ),0u,TUPLE_FLAG_NULL>;
    template<class O, Tuple_flag OFLAG> auto mkBackMult(const Tuple<O,0u,OFLAG> & a)const -> Tuple<decltype(data.data[0]* a[0] ),0u,TUPLE_FLAG_NULL>;
    template<class O> auto operator*(const SparseTuple<O> & a)const -> Tuple<decltype(data.data[0]* a[0] ),0u,TUPLE_FLAG_NULL>;
    template<class O> auto mkBackMult(const SparseTuple<O> & a)const -> Tuple<decltype(data.data[0]* a[0] ),0u,TUPLE_FLAG_NULL>;

    TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkTranspose() const;
    TMatrix<C,0u,0u,TF>& toTranspose();
    TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkTrju() const;
    TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> mkInverse() const;
    TMatrix<C,0u,0u,TUPLE_FLAG_CONSTRUCT_MASK(TF)> toInverse();
    #ifdef Rcpp_hpp
        void rdMatrix(const Rcpp::NumericMatrix &object);
        void wrMatrix(Rcpp::NumericMatrix &object)const;
    #endif
};

template<class Key>
class defaultHashFnc{
    public:
       static inline unsigned int makeSeed(const Key &inv) {  return (unsigned int) inv;  }
};
template<class Key>
class defaultHashFnc<EnumBox<Key> >{
    public:
       static inline unsigned int makeSeed(const EnumBox<Key> &inv) {  return (unsigned int) inv.data;  }
};
template<class C, unsigned int SIZE>
class defaultHashFnc< Tuple< C, SIZE>  >{
    public:
       static inline unsigned int makeSeed(const Tuple< C, SIZE> &inv) {
        unsigned int fout = defaultHashFnc<C>::makeSeed(inv[0]);
        for(unsigned int i=1;i< SIZE;i++){
            fout = defaultHashFnc<C>::makeSeed(inv[i]) ^ (fout << 1);
            }
        return (unsigned int) fout;}
};
template<class C>
class defaultHashFnc< Tuple< C, 0u>  >{
    public:
       static inline unsigned int makeSeed(const Tuple< C, 0u> &inv) {
        unsigned int fout = defaultHashFnc<C>::makeSeed(inv[0]);
        for(unsigned int i=1;i< inv.tup_size;i++){
            fout = defaultHashFnc<C>::makeSeed(inv[i]) ^ ((fout << 3) | (fout >> 29));
            }
        return (unsigned int) fout;}
};
template<class Key, class Key2>
class defaultHashFnc< pair< Key, Key2>  >{
    public:
       static inline unsigned int makeSeed(const pair< Key, Key2> &in){ return defaultHashFnc<Key>::makeSeed(in.first) ^ defaultHashFnc<Key2>::makeSeed(in.second) ;  }
};
template<class Key, class Key2>
class defaultHashFnc< KeyElem< Key, Key2>  >{
    public:
       static inline unsigned int makeSeed(const KeyElem< Key, Key2> &in){ return defaultHashFnc<Key>::makeSeed(in.first) ^ defaultHashFnc<Key2>::makeSeed(in.second) ;  }
};
template< > class defaultHashFnc<string>{
public: static inline unsigned int makeSeed(const string __s){   unsigned long __h = 0;
const char* tmp = __s.c_str();
for ( ; *tmp; ++tmp)
__h = 5*__h + *tmp;
return size_t(__h); }
};
template<unsigned int SIZE > class defaultHashFnc<Tuple<char,SIZE> >{
public: static inline unsigned int makeSeed(const Tuple<char,4> __s) {   unsigned long __h = 0;
	const char* tmp = & __s[0];
	for ( ; (unsigned int)(tmp - &__s[0]) < SIZE ; ++tmp)
		__h = 5*__h + *tmp;
return size_t(__h); }
};
// WARNING! this hash function assumes a '/0' exists at all
template< > class defaultHashFnc<char*>{
  public: static inline unsigned int makeSeed(const char* __s) {   unsigned long __h = 0;
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  return size_t(__h); }
};
// WARNING! this hash function assumes a '/0' exists at all
template< > class defaultHashFnc<const char*>{
  public: static inline unsigned int makeSeed(const char* __s) {   unsigned long __h = 0;
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  return size_t(__h); }
};
template< > class defaultHashFnc<char> { public: static inline unsigned int  makeSeed(char __x) { return __x; }};
template< > class defaultHashFnc<unsigned char> { public: static inline unsigned int makeSeed(unsigned char __x) { return __x; }};
template< > class defaultHashFnc<signed char> { public: static inline unsigned int makeSeed(unsigned char __x)  { return __x; }};
template< > class defaultHashFnc<short> { public: static inline unsigned int makeSeed(short __x)  { return __x; }};
template< > class defaultHashFnc<unsigned short> { public: static inline unsigned int makeSeed(unsigned short __x)  { return __x; }};
template< > class defaultHashFnc<int> { public: static inline unsigned int makeSeed(int __x)  { return __x; }};
template< > class defaultHashFnc<unsigned int> { public:  static inline unsigned int makeSeed(unsigned int __x)  { return __x; }};
template< > class defaultHashFnc<long> { public: static inline unsigned int makeSeed(long __x)  { return __x; }};
template< > class defaultHashFnc<unsigned long> { public: static inline unsigned int makeSeed(unsigned long __x)  { return __x; }};
template<class C> class defaultHashFnc<C*> { public: static inline unsigned int makeSeed(C* inv) {
	unsigned long __h = 0;
	const char* tmp = (const char*)  &inv;
	for ( ; (unsigned int)(tmp - ((const char*)&inv)) < sizeof(C*) ; ++tmp)
		__h = 5*__h + *tmp;
	return size_t(__h);
}};  // from address
template <class C>
class Vector{
	template<int supersize> Vector<C> fourierTransform_routine(const Vector<mycomplex>& bluewindow) const;
	template<int supersize> Vector<C> invfourierTransform_routine(const Vector<mycomplex>& bluewindow) const;
    void fourierTransform_routine(Vector<C> & fout, const mycomplex * bluewindow, unsigned char mag) const;
    void invfourierTransform_routine(Vector<C> & fout,const mycomplex * bluewindow, unsigned char mag) const;
    void setSize_init(uint32_t nsize);
public:
    static const uint32_t DIMENTIONALITY = 1u;
	typedef uint32_t KEYITERATORMKR_TYPE;
	uint32_t asize;
	C* darray;

	typedef Vector<C> SAFETYPE;
	typedef DataGrid<C,2> TRJU_TYPE;
	typedef DataGrid<C,2> LMUL_TYPE;
	typedef unsigned int ITERATOR_TYPE;
	typedef C INNER_TYPE;
	typedef TMatrix<C> OUTER_TYPE;

	ExCoMeMdEcLaRe( Vector<C> )

	class Iterator{
    public:
    	typedef uint32_t KEYITERATOR_TYPE;
        uint32_t i,maxval;
        Vector<C>& target;
        Iterator(Vector<C>& trg): target(trg), i(0),maxval(trg.getSize()){}
        Iterator(Vector<C>& trg, uint32_t _start, uint32_t _endoffset): target(trg), i(_start),maxval(_endoffset){}
        operator LFHPrimitive::Iterator<uint32_t, C>(){using std::placeholders::_1; using std::placeholders::_2; return LFHPrimitive::Iterator<uint32_t, C>(std::bind( &Iterator::next , this, _1,_2), &i, (i != maxval) ? &target[i] :NULL);}
        operator bool (){return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        bool next(uint32_t*& keyt, C*& itar){if (++i == maxval) return false; itar = &(target[i]); return true;}
        uint32_t operator()()const{return i;}
        uint32_t getOffset()const{return i;}
        C* operator->(){return target.darray + i;}
        C& operator*(){return target.darray[i];}
        Iterator& mkIterator(){return *this;}
	};
    class ConstIterator{
    public:
    	typedef const uint32_t KEYITERATOR_TYPE;
        uint32_t i,maxval;
        const Vector<C>& target;
        ConstIterator(const Vector<C>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t operator()()const{return i;}
        uint32_t getOffset()const{return i;}
        const C* operator->(){return target.darray + i;}
        const C& operator*(){return target.darray[i];}
        ConstIterator& mkIterator(){return *this;}
	};

	template<int ISCONST>
	class MetaIterator{
	public:
		typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
        uint32_t i,maxval;
        decltype(C().mkIterator()) *subite;
        typedef decltype( *(const_cast<const C*>(new C())->mkIterator()) ) INNER_TYPE;
        typename MetaType<Vector<C>,ISCONST>::IS_CONST_REF target;
        MetaIterator(typename MetaType<Vector<C>,ISCONST>::IS_CONST_REF trg): target(trg), subite(NULL) {}
        MetaIterator(const MetaIterator<ISCONST>& trg);
        ~MetaIterator(){if (subite) delete(subite);}

        operator bool ();
        bool operator++(int);
        uint32_t operator()()const{return i;}

        MetaIterator<ISCONST>& operator=(const MetaIterator<ISCONST>& ite);
        // const INNER_TYPE* operator->(){return (*subite)->;}
        typename MetaType<INNER_TYPE,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return **subite;}
        MetaIterator<ISCONST>& mkIterator(){return *this;}
	};
/*
	template<int ISCONST>
	class MetaIteratorIsaPtr{
	public:
		typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
        uint32_t i,maxval;
        decltype(C()->mkIterator()) *subite;
        typedef decltype( **(const_cast<const C*>(new C())->mkIterator()) ) INNER_TYPE;
        typename MetaType<Vector<C>,ISCONST>::IS_CONST_REF target;
        MetaIteratorIsaPtr(typename MetaType<Vector<C>,ISCONST>::IS_CONST_REF trg): target(trg){}
        operator bool ();
        bool operator++(int);
        uint32_t operator()()const{return i;}
        // const INNER_TYPE* operator->(){return (*subite)->;}
        typename MetaType<INNER_TYPE,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return **subite;}
        MetaIteratorIsaPtr<ISCONST>& mkIterator(){return *this;}
	};*/


	Vector(): asize(0){}
	Vector(const Vector<C> & v){setSize_init(v.getSize()); for(unsigned int i=0;i<v.getSize();i++) darray[i] = v[i];}
	Vector(Vector<C>&& v):asize(v.asize),darray(v.darray){v.asize =0;}
	Vector<C>& operator=(const Vector<C> & v);
	Vector<C>& operator=(Vector<C>&& v){this->toMemfree(); asize = v.asize; darray = v.darray; v.asize=0; return *this;}
	~Vector(){if (asize != 0) delete[](darray);}
    Vector(std::initializer_list<C>);
    template<class O> Vector(std::initializer_list<O>);

	Vector(const C* data, int size);
	template<unsigned int S, Tuple_flag Cflag> Vector(const Tuple<C,S,Cflag> &data){setSize_init(data.getSize());for(unsigned int i=0;i<data.getSize();i++) darray[i] = data[i];}
	template<unsigned int S, Tuple_flag Cflag> Vector(Tuple<C,S,Cflag> &&data){setSize_init(data.getSize());for(unsigned int i=0;i<data.getSize();i++) darray[i] = std::move(data[i]);}
	template<class O, unsigned int Osize, Tuple_flag Cflag> Vector(const Tuple<O,Osize,Cflag> &data){setSize_init(data.getSize());for(unsigned int i=0;i<data.getSize();i++) darray[i] = (C)data[i];}
    template<class D> Vector(const Vector<D>& data){setSize_init(data.getSize()); for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] = (C) data.darray[i];}
    template<class O> Vector(const std::vector<O> &other){this->setSize_init(other.size());for(uint32_t i=0;i<getSize();i++) darray[i] = other[i];}

    Vector<C>::Iterator getIterator(){return Vector<C>::Iterator(*this);}
    Vector<C>::ConstIterator getIterator()const{return Vector<C>::ConstIterator(*this);}
    Vector<C>::Iterator mkIterator(){return Vector<C>::Iterator(*this);}
    Vector<C>::ConstIterator mkIterator()const{return Vector<C>::ConstIterator(*this);}
    Vector<C>::Iterator getPartitionIterator(int which, int nbparts){return Iterator(*this, (which * getSize()) / nbparts, ((which+1) * getSize()) / nbparts );}
    Vector<C>::ConstIterator getPartitionIterator(int which, int nbparts)const{return ConstIterator(*this, (which * getSize()) / nbparts, ((which+1) * getSize()) / nbparts );}

	Vector<C>::MetaIterator<0> mkPartitionedMetaIterator(){return Vector<C>::MetaIterator<0>(*this);}
    Vector<C>::MetaIterator<1> mkPartitionedMetaIterator()const{return Vector<C>::MetaIterator<1>(*this);}
//	Vector<C>::MetaIterator<0> mkPartitionedMetaIteratorIsaPtr(){return Vector<C>::MetaIteratorIsaPtr<0>(*this);}
//  Vector<C>::MetaIterator<1> mkPartitionedMetaIteratorIsaPtr()const{return Vector<C>::MetaIteratorIsaPtr<1>(*this);}


    operator Accessor<unsigned int, C> (){using std::placeholders::_1; using std::placeholders::_2;  return LFHPrimitive::Accessor<unsigned int, C>(std::bind( &access, this, _1,_2));}
    operator std::function<C(uint32_t)> ()const {using std::placeholders::_1; return  std::bind(&operator[], this, _1);}
    bool access(const uint32_t &index, C*& itar){itar = darray + index; return(index < getSize());}

	void setLinkMemory(void* _newmem);
	void* getLinkMemory() const;

	Vector<C>& toZero(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toZero(darray[i]); return(*this);}
	Vector<C>& toOne(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toOne(darray[i]); return(*this);}
	Vector<C>& toRand(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toRand(darray[i]); return(*this);}

	Vector<C>& toMin(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toMin(darray[i]); return(*this);}
	Vector<C>& toMax(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toMax(darray[i]); return(*this);}
	Vector<C>& toMemfree();

	unsigned int allocatedSize() const;


	template<class OC> Vector<C>& operator=(const Vector<OC> & v);
	template<class OC> Vector<C>& toConvert(const Vector<OC> & v);

	Vector<C>& toMemmove(Vector<C>& source);
	Vector<C>& toMemswap(Vector<C>& source);

	template<class A_1> void operator() (Oper1<A_1> const & op); // not a match
	void operator() (Oper1< C> const & op); // match

	template<class A_1, class A_2, class C_I> void operator() (Oper2<A_1,A_2> const & op, Vector<C_I> const & _in ); // not a match
	template<class C_I> void operator() (Oper2<C,C_I> const & op, Vector<C_I> const & _in); // match
	template<class A_1, class A_2, class C_I> void operator() (Oper2<A_1,A_2> const & op, Vector<C_I> & _in ); // not a match
	template<class C_I> void operator() (Oper2<C,C_I> const & op, Vector<C_I>  & _in); // match

	template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> const &, Vector<C_3> const &); // not a match
	template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> const &, Vector<C_3> const &); // match
	template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2>  &, Vector<C_3> const &); // not a match
	template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> &, Vector<C_3> const &); // match
	template<class A_1, class A_2, class A_3, class C_2, class C_3> void operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2>  &, Vector<C_3>  &); // not a match
	template<class C_2, class C_3> void operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2>  &, Vector<C_3> &); // match


	template<class O> Vector<C>& toOrderOf(const O*); // init C using ints, size needs to be predetermined


	C& mempush_back(C& entry) {return ExOp::toMemmove(push_back(),entry);}
	C& push_back(const C& entry) {return (push_back() = entry);}
	C& push_back(C&& entry);
	C& push_back();
    C& push_back_within_copy(unsigned int offset){C& target = push_back();return target = (*this)[offset];}
    C& push_back_within_move(unsigned int offset){C& target = push_back();return ExOp::toMemmove(target,(*this)[offset]);}
	void pop_back();
	void pop_swap(unsigned int a);
	void insert_at_position(unsigned int pos, const C&); // may be slow, shift most entries at any beyong "pos" by 1

	//	void up_alloc();
	//	void down_alloc();

	unsigned int getSize() const;
	unsigned int size() const;
	Vector<C>& setSize(unsigned int nsize);
	Vector<C>& DownSize(unsigned int nsize); // conserve the first elems;
	Vector<C>& upSize(unsigned int nsize); // conserve the first elems; assumes nsize > asize
    Vector<C>& toAppend(const Vector<C>& other);
    Vector<C>& toMemappend(Vector<C>& other);
    Vector<C>& toMemappendInverted(Vector<C>& other);


	const C& last(int i = 1) const{return(darray[(asize & 0x7FFFFFFF)-i]);}
	C& last(int i = 1){return(darray[(asize & 0x7FFFFFFF)-i]);}
	const C* begin() const{return(darray);}
	const C* end() const{return(darray + (asize & 0x7FFFFFFF));}
	C* begin(){return(darray);}
	C* end(){return(darray + (asize & 0x7FFFFFFF));}


	C & operator[](int const which){return(darray[which]);}
	const C & operator[](int const which)const{/*printf("index access... %i/%i\n",which, this->getSize()); fflush(stdout);*/ return(darray[which]);}

	//template<int flag2> operator PolyThing<C2> () const;
	template<int size> operator Tuple<C,size> () const;

	uint32_t searchFor(const C& what, bool _is_sorted = false) const;

	Vector<C> operator-()const {Vector<C> fout; fout.setSize(this->getSize()); for(int i =0;i< this->getSize();i++) fout[i] = -darray[i]; return fout;}
	Vector<C> operator!()const {Vector<C> fout; fout.setSize(this->getSize()); for(int i =0;i< this->getSize();i++) fout[i] = !darray[i]; return fout;}
	template<class O> Vector<C> operator+(const O &other)const {return( ((Vector<C>(*this)) += other) );}
	template<class O> Vector<C> operator-(const O &other)const {return( ((Vector<C>(*this)) -= other) );}
	template<class O> Vector<C> operator*(const O &other)const {return( ((Vector<C>(*this)) *= other) );}
	template<class O> Vector<C> operator/(const O &other)const {return( ((Vector<C>(*this)) /= other) );}


	template<class D> Vector<C>& operator+=(const Vector<D> & v);
	template<class D> Vector<C>& operator+=(const D & v);
	template<class D> Vector<C>& operator-=(const Vector<D> & v);
	template<class D> Vector<C>& operator-=(const D & v);
	template<class D> Vector<C>& operator*=(const Vector<D> & v);
	template<class D> Vector<C>& operator*=(const D & v);
	template<class D> Vector<C>& operator/=(const Vector<D> & v);
	template<class D> Vector<C>& operator/=(const D & v);

	LFH_FAULTY Vector<unsigned int> mksortindex() const; // makes a vector of indexes which access elements in a increasing order


    Vector<C>& sortAndCountPermutation(uint32_t &nbperm);
    Vector<C>& sort();
    Vector<C>& toRandPermute();
    Vector<C>& sortWithinRange(int first, int limit);
    uint32_t sortWithinRangeAndSumPerm(int first, int limit);
	Vector<C>& sort_unique();
    Vector<C>& sortFirstN(uint32_t nbtosort);

  	Vector<C>& sort_decr(); // Decreasing Order
  	Vector<C>& sortWithinRange_decreasing(int first, int limit); // Decreasing Order


    /* These method uses binary search, and assumes the array is sorted*/
    uint32_t findLE(const C& what)const;
    template<class D> uint32_t findLE(const D& what)const;


	template<class COMP> Vector<C>& sort_comp(const COMP& comparator); // COMP needs int operator(const &C, const &C);
	bool issorted() const;
	template<class COMP> bool issorted_comp(const COMP& comparator) const; // COMP needs int operator(const &C, const &C);
	bool issorted_decr() const; // Decreasing Order
	Vector<C>& reverse();
	Vector<C>& reverseLast(uint32_t nbpermuted);
	void random_permute();
	template<class D> void getIntersection(const IntervalSet<D> & query , Vector<C> &out) const;


	static Vector<mycomplex> bluesteinWindow(unsigned int size);
	static Vector<mycomplex> bluesteinWindow_new(unsigned int size);
//	static void bluesteinWindow(complex *&, unsigned int size);
	Vector<C> fourierTransform(const Vector<mycomplex>& bluewindow) const;
	Vector<C> invfourierTransform(const Vector<mycomplex>& bluewindow) const;

	void fourierTransform_semiroutine(const Vector<mycomplex>& bluewindow);
	void invfourierTransform_semiroutine(const Vector<mycomplex>& bluewindow);

	Matrix<C> mkouterprod() const;
/*
	void GP_Covariance(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), double noise_variance =0.0f) const;
	void GP_Covariance_cross(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), Vector<C>& query) const;
	void GP_fromCorrel(DataGrid<double, 2> &f_out, double (*correl)(const C&, const C&)) const;
	void GP_fromCorrel_cross(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const;
	void GP_fromCorrel_complete(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const;
*/
	template<class O> std::vector<O>& wrStdVector(std::vector<O>& fout)const;


#ifdef Rcpp_hpp
    template<int D> void wrRCPP(Rcpp::Vector<D> &fout) const;
    template<int D> void rdRCPP(const Rcpp::Vector<D> &fin);
#endif
    void showElem(uint32_t index, FILE* f, int lvl) const{ ExOp::show(darray[index], f, lvl);}
};
// Objects are to be extracted or inserted elsewhere, but target might be locked by mutex
// this structure temporary hold elements needing to be inserted, so that "writing" threads
// do not get locked, "reading" thread would need to insert elements before any read is allowed
template < > class Vector<Anything>{
    uint8_t* push_back_routine();
    uint8_t* assess_routine(int pos){return(darray + (anytype & 0xFF) * pos);}
    const uint8_t* assess_routine(int pos)const {return(darray + (anytype & 0xFF) * pos);}
    Vector<Anything>& setSize_init(uint32_t nsize);
public:
    static const uint32_t DIMENTIONALITY = 1u;
	uint32_t asize;
	uint32_t anytype;
	uint8_t* darray;
	typedef Vector<Anything> SAFETYPE;
	typedef Anything INNER_TYPE;
	typedef unsigned int ITERATOR_TYPE;
	ExCoMeMdEcLaRe( Vector<Anything> )
	class Iterator{
    public:
        uint32_t i,maxval;
        Vector<Anything>& target;
        Iterator(Vector<Anything>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t operator()()const{return i;}
    //  C* operator->(){return target.darray + i;}
    //  C& operator*(){return target.darray[i];}
	};
    class ConstIterator{
    public:
        uint32_t i,maxval;
        const Vector<Anything>& target;
        ConstIterator(const Vector<Anything>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t operator()()const{return i;}
    //  const C* operator->(){return target.darray + i;}
    //  const C& operator*(){return target.darray[i];}
	};

	Vector(): asize(0){}
	Vector(const Vector<Anything> & v):anytype(v.anytype) {setSize_init(v.getSize()); memcpy(darray,v.darray, v.getSize() * (v.anytype & 0xFF));}
	Vector(Vector<Anything>&& v):asize(v.asize),anytype(v.anytype),darray(v.darray){v.asize =0;}
	Vector<Anything>& operator=(const Vector<Anything> & v){if (asize != 0) delete[](darray); anytype = v.anytype; setSize_init(v.getSize()); memcpy(darray,v.darray, v.getSize() * (v.anytype & 0xFF)); return *this;}
	Vector<Anything>& operator=(Vector<Anything>&& v){this->toMemfree(); asize = v.asize; anytype = v.anytype; darray = v.darray; v.asize=0; return *this;}
	~Vector(){if (asize != 0) delete[](darray);}

	Vector<Anything>::Iterator getIterator(){return Vector<Anything>::Iterator(*this);}
	Vector<Anything>::ConstIterator getIterator()const{return Vector<Anything>::ConstIterator(*this);}
	Vector<Anything>& toZero(){unsigned int s = getSize();for(unsigned int i=0;i<s;i++) ExOp::toZero(darray[i]); return(*this);}
	Vector<Anything>& toMemfree();

	unsigned int allocatedSize() const;

	Vector<Anything>& toMemmove(Vector<Anything>& source);
	Vector<Anything>& toMemswap(Vector<Anything>& source);

//	C& mempush_back(C& entry) {return ExOp::toMemmove(push_back(),entry);}
//	C& push_back(const C& entry) {return (push_back() = entry);}
//	C& push_back(C&& entry);

	Anything operator[](int pos){return(Anything(assess_routine(pos),anytype));}
	AnythingConst operator[](int pos)const{return(AnythingConst(assess_routine(pos),anytype));}

    	Anything push_back(){return(Anything(this->push_back_routine(),anytype));}
	template<class C> void push_back(const C& val){Anything(this->push_back_routine(),anytype) = val;}

	void pop_back();
	void pop_swap(unsigned int a);

	unsigned int getSize() const;
	unsigned int size() const;

	Vector<Anything>& setType(uint8_t _anytype){anytype = _anytype; return(*this);}
	Vector<Anything>& setSize(uint32_t nsize);
	void DownSize(unsigned int nsize); // conserve the first elems;
	void upSize(unsigned int nsize); // conserve the first elems; assumes nsize > asize
	Vector<Anything>& toAppend(const Vector<Anything>& other);
	Vector<Anything>& toMemappend(Vector<Anything>& other);

	void reverse();

#ifdef Rcpp_hpp
    template<int D> void wrRCPP(Rcpp::Vector<D> &fout) const;
    template<int D> void rdRCPP(const Rcpp::Vector<D> &fin);
#endif
    void showElem(uint32_t index, FILE* f, int lvl) const{return;}
};
template < > class Vector<void>{// This class is used by other structures where no data is held, the index itself is returned instead

    public:
    inline uint32_t operator[](uint32_t val)const{return val;}
    inline void push_back(){return;}
    inline Vector<void>& toMemfree(){return *this;} // nothing toDo
    inline Vector<void>& toZero(){return *this;} // nothing toDo
    inline Vector<void>& toOne(){return *this;} // nothing toDo
    inline void showElem(uint32_t index, FILE* f, int lvl) const{return;}
    inline ERRCODE save(FILE*f) const{return 0;}
    inline ERRCODE load(FILE*f){return 0;}
};
template<class C, int SIZE_MAG>
class AsyncBuffer{
    atomic<uint32_t> sema; // MK2
    std::mutex mtx; // MK2
    std::mutex overmtx;  // MK2
public:
    C buffer[StConstList<SIZE_MAG>::twopower];
    std::stack<C> overflow;
    // tries to read,
    class ReadAccess{
        AsyncBuffer& target;
        public:
        std::lock_guard<std::mutex> mlock;
        ReadAccess(AsyncBuffer& _target);
        bool needsUpdate();
        C update();
    };
    class WriteAccess{
        AsyncBuffer& target;
        public:
        std::lock_guard<std::mutex> mlock;
        WriteAccess(AsyncBuffer& _target);

    };

    std::atomic<int> sema_wr;
    std::atomic<int> sema_rd;
    uint32_t pos_wr;
    uint32_t pos_rd;

    AsyncBuffer();

    void insert(const C &what);

    ERRCODE rdSync(C& where); // only use if no other thread can read;
    ERRCODE rdAsync(C& where);
    ERRCODE wrSync(const C& what); // only use if no other thread can write;
    ERRCODE wrAsync(const C& what);
    ERRCODE mvSync(C& what); // only use if no other thread can write;
    ERRCODE mvAsync(C& what);
    void show(FILE*f= stdout, int level=0) const;
};

template <class D>
class defaultHeapBackPointer{
public:
	unsigned int backptr;
	D data;

	defaultHeapBackPointer<D>& operator=(const defaultHeapBackPointer<D> &other){data = other.data; return(*this);}
	bool operator>(const defaultHeapBackPointer<D> &other) const {return( data > other.data);}
	bool operator<(const defaultHeapBackPointer<D> &other) const {return( data < other.data);}
	bool operator>=(const defaultHeapBackPointer<D> &other) const {return( data >= other.data);}
	bool operator<=(const defaultHeapBackPointer<D> &other) const {return( data <= other.data);}
	bool operator==(const defaultHeapBackPointer<D> &other) const {return( data == other.data);}
	bool operator!=(const defaultHeapBackPointer<D> &other) const {return( data != other.data);}
};
template<class C> class HeapTreeBackPointerOffset{
public:
	enum{ans = 0}; // put non-zero if non trivial
	// static int& getBackPtr(C&)
};
template<class D>
class HeapTreeBackPointerOffset< defaultHeapBackPointer<D> >{
public:
	enum{ans = 1};
	static unsigned int& getBackOffset(defaultHeapBackPointer<D>& targ){return(targ.backptr);}
};
// heapTree, may contain an unsorted node at position 0!
template <class C, unsigned int ASYNC_MAG>
class HeapTree{
    typedef C ELEM_TYPE;
    Vector<C> data;
	bool hasunsorted;
	HeapTree<C,0> heaptree;
    AsyncBuffer<C, ASYNC_MAG> asyncbuffer;
    public:
	HeapTree() : hasunsorted(false){data.setSize(1);}
	HeapTree& operator=(const HeapTree& other)=delete;
	bool isEmpty();
	bool top(C& fout);
	bool pop(C& fout);
	void insert(const C&);
	void show(FILE*f = stdout, int level=0)const;
};
template <class C>
class HeapTree<C, 0u>{ // Synced heaptree
public:
    typedef C ELEM_TYPE;
    Vector<C> data;
	bool hasunsorted;
	HeapTree() : hasunsorted(false){data.setSize(1);}

	bool isEmpty() const{return((!hasunsorted) && (data.getSize() == 1));}
    unsigned int getSize()const {return data.getSize() - ( ((hasunsorted) || (data.getSize() == 0)) ? 0: 1);}
    HeapTree& toMemfree(){data.setSize(1);hasunsorted = false;return *this;}
    HeapTree& toMemmove(HeapTree& other){data.toMemmove(other.data); other.data.setSize(1);hasunsorted = other.hasunsorted; other.hasunsorted = false;return *this;}

	void insert(const C &what);
	HeapTree& operator<<=(C &what); // memory insert!

    // {*assumes non-empty*
	C top() const{return( data[ (  (hasunsorted)&&((data.getSize() == 1)||(data[0] < data[1])) ? 0 : 1) ]);}
	C& top(){return( data[ (  (hasunsorted)&&((data.getSize() == 1)||(data[0] < data[1])) ? 0 : 1) ]);}
	void updateTop(){this->update( (hasunsorted)&&(((*this).getSize() == 1)||((*this)[0] < (*this)[1])) ? 0 : 1);}
    void update(unsigned int offset); // if top was altered, move it within HeapTree
	C pop();
	void pop_and_insert(const C &what);
    // }*assumes non-empty*
	void show(FILE*f = stdout, int level=0)const;
};
template <class C>
class HeapTree<C*, 0u> : public Vector<C*>{
	public:
	bool hasunsorted;
	HeapTree();

	//	Node& operator()(const Key & where);
	bool isEmpty();
	void insert(C* what);

	// {*assumes non-empty*
	C& pop();
	void updateTop(){this->update( (hasunsorted)&&(((*this).getSize() == 1)||((*this)[0] < (*this)[1])) ? 0 : 1);}
	void update(unsigned int offset); // Elem was changed! update its position!
	C& top()const{return( *(*this)[ (  (hasunsorted)&&(((*this).getSize() == 1)||((*this)[0] < (*this)[1])) ? 0 : 1) ]);}

	//	KeyElem<Key, Node> pop();
//	KeyElem<Key, Node> ipop(Node what, Key where); // insert and pop!

};
class FoundIndex{
public:
    uint32_t value;
    FoundIndex(uint32_t _value): value(_value){}
    operator uint32_t() const{return value;}
    operator bool ()const{return ((value + 1) != 0);}
};

template<class Key, class Data, class HashFnc>
class myHashmap{
    typedef typename ExCo<Key>::TYPE KeyType;
    unsigned int hashpos(unsigned int seed) const;
    void swap_and_pop(unsigned int ite); // moves entry in vector, into an unlinked location
    void showLinks(FILE* f = stdout) const;
    void remake_links_routine();
    public:

    Vector< KeyElem< KeyElem<Key, Data> , unsigned int > > heap;
    typedef unsigned int ITERATOR_TYPE;
	typedef Key KEYITERATORMKR_TYPE;



    typedef myHashmap<Key, Data, HashFnc> SAFETYPE;
    unsigned int* hash;

    unsigned char hash_mag;
    class ConstIterator{
        uint32_t i,maxval;
    public:
    	typedef const Key KEYITERATOR_TYPE;
        const myHashmap<Key, Data, HashFnc>& target;
        ConstIterator(const myHashmap<Key, Data, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const {return i;}
        Key operator()()const{return target.deref_key(i);}
        const Data* operator->(){return &(target.deref(i));}
        Data operator*(){return target.deref(i);}
        ConstIterator& mkIterator(){return *this;} // no need to make, it is itself in this case
	};

    class Iterator{
    public:
    	typedef Key KEYITERATOR_TYPE;
        uint32_t i,maxval;
        myHashmap<Key, Data, HashFnc>& target;
        Iterator(myHashmap<Key, Data, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const {return i;}
        Key operator()()const{return target.deref_key(i);}
        Data* operator->(){return &(target.deref(i));}
        Data& operator*(){return target.deref(i);}
        Iterator& mkIterator(){return *this;} // no need to make, it is itself in this case

	};

	template<int ISCONST> class SlotIterator{
    public:
        const Key key;
        uint32_t ite;
        typename MetaType<myHashmap<Key, Data, HashFnc>,ISCONST>::IS_CONST_REF target;
        SlotIterator(typename MetaType<myHashmap<Key, Data, HashFnc>,ISCONST>::IS_CONST_REF trg, const Key& _k): target(trg), key(_k){}
        operator bool (){if ((target.hash_mag == 0)||((ite = target.hash[target.hashpos( HashFnc::makeSeed(key) )]) == 0xFFFFFFFF)) return false; return (target.heap[ite].k.k != key) ? (*this)++ : true;}
        bool operator++(int){ do{ite = target.heap[ite].d; if (ite == 0xFFFFFFFF) return false;}while(target.heap[ite].k.k != key); return true;}
        Key operator()()const{return key;}
        typename MetaType<Data*,ISCONST>::IS_CONST operator->(){return &(target.deref(ite));}
        typename MetaType<Data,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.deref(ite);}
		void remove(){if (ite == 0xFFFFFFFF) printf("trying to erase with invalid!\n"); else target.erase_from_iterator(ite);}
        void show(FILE* f =stdout, int level=0)const{printf("Valid Slot iterator with key: "); ExOp::show(key,f, 0);}
	};

    myHashmap(): hash(nullptr),hash_mag(0){}
    myHashmap(const myHashmap<Key, Data, HashFnc>& other);
    myHashmap(myHashmap<Key, Data, HashFnc>&& other);
    myHashmap<Key, Data, HashFnc>& operator=(const myHashmap<Key, Data, HashFnc>& other);
    myHashmap<Key, Data, HashFnc>& operator=(myHashmap<Key, Data, HashFnc>&& other);
    ~myHashmap();

    Data operator[](const Key &) const; // assumes entry exists, use "find" otherwise
    Data& operator[](const Key &); // creates entry if it does not exist

	KeyElem<Key, Data>& getKeyElem(const Key &);

    unsigned int find(const Key &) const;

    myHashmap<Key, Data, HashFnc>::Iterator mkIterator(){return myHashmap<Key, Data, HashFnc>::Iterator(*this);}
    myHashmap<Key, Data, HashFnc>::ConstIterator mkIterator()const{return myHashmap<Key, Data, HashFnc>::ConstIterator(*this);}

    myHashmap<Key, Data, HashFnc>::Iterator operator()(){return myHashmap<Key, Data, HashFnc>::Iterator(*this);}
    myHashmap<Key, Data, HashFnc>::ConstIterator operator()()const{return myHashmap<Key, Data, HashFnc>::ConstIterator(*this);}
    myHashmap<Key, Data, HashFnc>::SlotIterator<0> operator()(const Key &key){return myHashmap<Key, Data, HashFnc>::SlotIterator<0>(*this, key);}
    myHashmap<Key, Data, HashFnc>::SlotIterator<1> operator()(const Key &key)const{return myHashmap<Key, Data, HashFnc>::SlotIterator<1>(*this, key);}

    Data& addEntry(const Key &); // add entry *can have more than 1 data with same key*
    template<class Data_COMP> bool remEntry(const Key &, const  Data_COMP& ); // remove entry if matching key and == other

    myHashmap<Key, Data, HashFnc>& toMemfree(){heap.toMemfree(); if (hash_mag != 0) delete[](hash); hash_mag = 0; return(*this); }
    myHashmap<Key, Data, HashFnc>& toMemmove(myHashmap<Key, Data, HashFnc>& other);
    myHashmap<Key, Data, HashFnc>& toMemswap(myHashmap<Key, Data, HashFnc>& other);

    template<class D2> myHashmap<Key, Data, HashFnc>& operator=(const myHashmap<Key, D2, HashFnc>& other);


    template<class D2> myHashmap<Key, Data, HashFnc>& toConvert(const myHashmap<Key, D2, HashFnc>& other);

    myHashmap<Key, Data, HashFnc>& toZero(){return *this;}
    myHashmap<Key, Data, HashFnc>& toOne(){return *this;}


    Data& deref(unsigned int ite){return heap[ite].k.d;}
    const Data& deref(unsigned int ite) const {return heap[ite].k.d;}

    Key& deref_key(unsigned int ite){return heap[ite].k.k;}
    const Key& deref_key(unsigned int ite) const {return heap[ite].k.k;}


    void erase(const Key &);
    void erase_from_iterator(unsigned int ite);
    void changeKey(unsigned int ite, const Key&);
    inline unsigned int getSize() const{return heap.getSize();}

    void rehash(unsigned char _new_mag);
    void sort();
    void sort_decr();

    void keysort();
    void keysort_decr();
    void keysort_byteswap();
	template<class COMP> void sortKeyElem(const COMP &c);

    void insert_at_position(unsigned int pos, const Key &, const Data&); // may be slow, shift most entries at any beyong "pos" by 1

    ERRCODE load(const uint8_t * &chunk, const uint8_t* const endpos);
    ERRCODE save(FILE*f)const;
    ERRCODE load(FILE*f);
    void show(FILE*f = stdout, int level=0) const;
    class Compareclass{public:
        int operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const{
            if (ExOp::isEQ(a.k.d, b.k.d)) return 0;
            else return (ExOp::isGT(a.k.d, b.k.d)) ? 1 : -1;
        }
    };
    class Compareclass_decr{public:
        int operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const{
            if (ExOp::isEQ(a.k.d, b.k.d)) return 0;
            else return (ExOp::isGT(a.k.d, b.k.d)) ? -1 : 1;
        }
    };
    class Compareclass_key{public:
        int operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const{
            if (ExOp::isEQ(a.k.k, b.k.k)) return 0;
            else return (ExOp::isGT(a.k.k, b.k.k)) ? 1 : -1;
        }
    };
    class Compareclass_keydecr{public:
        int operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const{
            if (ExOp::isEQ(a.k.k, b.k.k)) return 0;
            else return (ExOp::isGT(a.k.k, b.k.k)) ? -1 : 1;
        }
    };
/*  class Compareclass_keybyteswap{public:
        int operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const;
    };*/
    void swapEntries(uint32_t iteA, uint32_t iteB);
};
template<class Key, class HashFnc> // lists keys, no data
class myHashmap<Key, void, HashFnc>{
    typedef typename ExCo<Key>::TYPE KeyType;
    unsigned int hashpos(unsigned int seed) const;
    void showLinks(FILE* f = stdout) const;
    void swap_and_pop(unsigned int ite);  // moves entry in vector, into an unlinked location
    public:
    typedef unsigned int ITERATOR_TYPE;
    Vector< pair< Key , unsigned int > > heap;
    typedef myHashmap<Key, void, HashFnc> SAFETYPE;
    unsigned int* hash;
    unsigned char hash_mag;

    myHashmap(): hash(nullptr),hash_mag(0){}
    myHashmap(const myHashmap<Key, void, HashFnc>&);
    myHashmap(myHashmap<Key, void, HashFnc>&&);
    myHashmap<Key, void, HashFnc>& operator=(const myHashmap<Key, void, HashFnc>& other);
    myHashmap<Key, void, HashFnc>& operator=(myHashmap<Key, void, HashFnc>&& other);
    ~myHashmap();
    // unsigned int find(const Key &) const;

    class ConstIterator{
        uint32_t i,maxval;
    public:
        const myHashmap<Key, void, HashFnc>& target;
        ConstIterator(const myHashmap<Key, void, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const { return i;}
        typename ExCo<Key>::INDEX_TYPE operator()()const{return target.deref_key(i);}
        const Key* operator->(){return &(target.deref(i));}
        Key operator*(){return target.deref(i);}
        void show(FILE *f =stdout, int level=0)const{fprintf(f, "const iterator at %i out of %i\n", i,maxval);}
	};
    class Iterator{
    uint32_t i,maxval;
    public:
        myHashmap<Key, void, HashFnc>& target;
        Iterator(myHashmap<Key, void, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const { return i;}
        typename ExCo<Key>::INDEX_TYPE operator()()const{return target.deref_key(i);}
        Key* operator->(){return &(target.keypreserving_deref(i));}
        Key& operator*(){return target.keypreserving_deref(i);}
        void show(FILE *f =stdout, int level=0)const{fprintf(f, "iterator at %i out of %i\n", i,maxval);}
	};
	template<int ISCONST> class SlotIterator{
    public:
        const typename ExCo<Key>::INDEX_TYPE & key;
        uint32_t ite;
        typename MetaType<myHashmap<Key, void, HashFnc>,ISCONST>::IS_CONST_REF target;
        SlotIterator(typename MetaType<myHashmap<Key, void, HashFnc>,ISCONST>::IS_CONST_REF& trg, const typename ExCo<Key>::INDEX_TYPE & _k): target(trg), key(_k){}
        operator bool (){if ((target.hash_mag == 0)||((ite = target.hash[target.hashpos( HashFnc::makeSeed(key) )] ) == 0xFFFFFFFF)) return false; return ExOp::getIndex(target.heap[ite].first) != key ? (*this)++ : true;}
        bool operator++(int){do{ite = target.heap[ite].second;if (ite == 0xFFFFFFFF) return false;} while(ExOp::getIndex(target.heap[ite].first) != key); return true;}
        typename ExCo<Key>::INDEX_TYPE operator()()const{return key;}
        typename MetaType<Key*,ISCONST>::IS_CONST operator->(){return &(target.deref(ite));}
        typename MetaType<Key,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.deref(ite);}
        void remove(){if (ite == 0xFFFFFFFF) printf("trying to erase with invalid!\n"); else target.erase_from_iterator(ite);}
	};

    myHashmap<Key, void, HashFnc>::Iterator mkIterator(bool modifications_should_preserve_index_reminder = false){return myHashmap<Key, void, HashFnc>::Iterator(*this);}
    myHashmap<Key, void, HashFnc>::ConstIterator mkIterator()const{return myHashmap<Key, void, HashFnc>::ConstIterator(*this);}


    myHashmap<Key, void, HashFnc>::SlotIterator<0> mkIterator(const typename ExCo<Key>::INDEX_TYPE &key){return myHashmap<Key, void, HashFnc>::SlotIterator<0>(*this, key);}
    myHashmap<Key, void, HashFnc>::SlotIterator<1> mkIterator(const typename ExCo<Key>::INDEX_TYPE &key)const{return myHashmap<Key, void, HashFnc>::SlotIterator<1>(*this, key);}

	// SYNONYMOUS TO
    myHashmap<Key, void, HashFnc>::Iterator operator()(bool modifications_should_preserve_index_reminder = false){return myHashmap<Key, void, HashFnc>::Iterator(*this);}
    myHashmap<Key, void, HashFnc>::ConstIterator operator()()const{return myHashmap<Key, void, HashFnc>::ConstIterator(*this);}
    myHashmap<Key, void, HashFnc>::SlotIterator<0> operator()(const typename ExCo<Key>::INDEX_TYPE &key){return myHashmap<Key, void, HashFnc>::SlotIterator<0>(*this, key);}
    myHashmap<Key, void, HashFnc>::SlotIterator<1> operator()(const typename ExCo<Key>::INDEX_TYPE &key)const{return myHashmap<Key, void, HashFnc>::SlotIterator<1>(*this, key);}
	// SYNONYMOUS TO

    unsigned int find(const typename ExCo<Key>::INDEX_TYPE &) const;

    myHashmap<Key, void, HashFnc>& toMemfree(){heap.toMemfree(); if (hash_mag != 0) delete[](hash); hash_mag = 0; return(*this); }
    myHashmap<Key, void, HashFnc>& toMemmove(myHashmap<Key, void, HashFnc>& other);
    myHashmap<Key, void, HashFnc>& toMemswap(myHashmap<Key, void, HashFnc>& other);

    Key operator[](const typename ExCo<Key>::INDEX_TYPE &) const;
    Key& operator[](const typename ExCo<Key>::INDEX_TYPE &);

	const Key& deref(unsigned int ite) const {return heap[ite].first;}
	Key& keypreserving_deref(unsigned int ite){return heap[ite].first;} // warning, assumes the key remains unchanged

	typename ExCo<Key>::INDEX_TYPE deref_key(unsigned int ite) const;

    template<class F> myHashmap<Key, void, HashFnc>& toConvert(F& , ITERABLE_DECL(F) );


    void addEntry(const Key &);
    Key& createEntry(const typename ExCo<Key>::INDEX_TYPE &);

    void erase(const Key &);
    void removeEntry(const typename ExCo<Key>::INDEX_TYPE &);
    void erase_from_iterator(unsigned int ite);
    inline unsigned int getSize() const{return heap.getSize();}

    void rehash(unsigned char _new_mag);
    ERRCODE save(FILE*f)const;
    ERRCODE load(FILE*f);
    void show(FILE*f = stdout, int level=0) const;
};
// key is an integer index of varying size
template<class C>
class myHashmap<void,C, defaultHashFnc<void> >{
public:
    Vector<C> data;
    Vector< Anything > heap;
    uint8_t* hash;
    uint32_t size;
    C& operator[](uint32_t index);
};
template<class Key, class HashFnc> // lists keys, no data
class myHashmap<Key, Anything, HashFnc>{
    typedef Key KeyType;
    unsigned int hashpos(unsigned int seed) const;
    void swap_and_pop(unsigned int ite);  // moves entry in vector, into an unlinked location
    public:
    typedef unsigned int ITERATOR_TYPE;
    Vector< pair< Key , unsigned int > > heap;
    Vector<Anything> data;

    typedef myHashmap<Key, Anything, HashFnc> SAFETYPE;
    unsigned int* hash;
    unsigned char hash_mag;

    myHashmap(): hash_mag(0),hash(nullptr){}
    ~myHashmap(){if (hash_mag) delete[](hash);}
    // unsigned int find(const Key &) const;

    class ConstIterator{
        uint32_t i,maxval;
    public:
        const myHashmap<Key, Anything, HashFnc>& target;
        ConstIterator(const myHashmap<Key, Anything, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const { return i;}
        typename ExCo<Key>::INDEX_TYPE operator()()const{return target.deref_key(i);}
        const Key* operator->(){return &(target.deref(i));}
        Key operator*(){return target.deref(i);}
	};
	/*	Hi Luz, I hope your fight with cell2location is going well, I just want to share the primary tissue data integrated with the Wang et al. dataset.
		The h5ad object can be browsed via cellxgene at

		and the raw file is at the folloing location on lustre:
		https://jupyter.cellgeni.sanger.ac.uk/user/lfhandfield/proxy/5008/?token=fd772e00fd854f6ab2db0880e1b5f66f
		/lustre/scratch117/cellgen/team292/lh20/1gpunobooks/primary/N5-integrated_donors.h5ad
		There is "general_celltypes" and "subcluster_epitelial_balanced" levels of general annotations for the epithelial cells, the later is defined for a subset of the epithelial cells (see "broad_celltypes"), since some samples with low quality have outstanding number of cells, so this was controlled by sammplingtheses out.

		I can make a epithelial subset h5ad if that is helpful, just that I thought that the idea is to get LR pair information for all celltypes to see communication, not exactly sure how one would interpret that in a "eprithelial only" subset.
		If limma is out of commission for some reason (as I heard) I was just wondering what format the DE method required to output for that cellphoneDB to use
		so I can see if I can do something about it, if that helps.
*/
    class Iterator{
        uint32_t i,maxval;
    public:
        myHashmap<Key, Anything, HashFnc>& target;
        Iterator(myHashmap<Key, Anything, HashFnc>& trg): target(trg){}
        operator bool (){i =0; maxval = target.getSize(); return (i != maxval);}
        bool operator++(int){return (++i != maxval);}
        uint32_t getOffset()const { return i;}
        typename ExCo<Key>::INDEX_TYPE operator()()const{return target.deref_key(i);}
        Key* operator->(){return &(target.deref(i));}
        Key& operator*(){return target.deref(i);}
	};
    myHashmap<Key, Anything, HashFnc>::Iterator getIterator(){return myHashmap<Key, Anything, HashFnc>::Iterator(*this);}
    myHashmap<Key, Anything, HashFnc>::ConstIterator getIterator()const{return myHashmap<Key, Anything, HashFnc>::ConstIterator(*this);}

    uint32_t find(const Key &) const;

    myHashmap<Key, Anything, HashFnc>& toMemfree(){heap.toMemfree(); if (hash_mag != 0) delete[](hash); hash_mag = 0; return(*this); }
    myHashmap<Key, Anything, HashFnc>& toMemmove(myHashmap<Key, Anything, HashFnc>& other);
    myHashmap<Key, Anything, HashFnc>& toMemswap(myHashmap<Key, Anything, HashFnc>& other);
    myHashmap<Key, Anything, HashFnc>& operator=(const myHashmap<Key, Anything, HashFnc>& other);

    AnythingConst operator[](const Key &) const;
    Anything operator[](const Key &);
	AnythingConst deref(unsigned int ite) const {return data[ite];}
	Anything keypreserving_deref(unsigned int ite){return data[ite];} // warning, assumes the key remains unchanged

	Key deref_key(unsigned int ite) const{return heap[ite].first;}

    Anything addEntry(const Key &);

    void erase(const Key &);
    void rmEntry(const Key &);
    void erase_from_iterator(unsigned int ite);
    inline unsigned int getSize() const{return heap.getSize();}

    void rehash(unsigned char _new_mag);
    ERRCODE save(FILE*f)const;
    ERRCODE load(FILE*f);
    void show(FILE*f = stdout, int level=0) const;
};


// hopefully dense int -> int mapping, used to fill near-dense data
// if density > 50%, store missing entries instead of existing
// designed for compression, but is mutable
class Sparsity{
    uint32_t asize; // number of items
    uint32_t dsize; // number of hash_items
    union{
        uint16_t* data16;
        uint32_t* data; // tries to store everything in 1 uint32_t, but uses 2 if impossible
    };
    // the hash maps Keys to Continuous range
    // *however* if the key is smaller than the size then the trivial 1:1 map is assumed unless the hashmap contains an entry
    // the "unused" slots are use point to the following unused slot

    // list needs data to make/store the full Key, Data offset and 'next' offset
    // hsize
    // the key needs maxmag - hashdepth (or hashmag if match) bits, and the Data needs sizemag bits, next need hashmag bits
    // hence sizemag + hashmag + maxmag - hashmag bits (+1 bit if we reuse original range)
    // we have that maxmag >= sizemag >= hashmag >= hashdepth , and we only get to choose hashdepth
    // the hash length is at least twice bigger than the number of item contained
    // three ranges of bits: hash-occupency bit, key-next, data

    // unused slots stores the position of next unused

    uint8_t hash_magn;  //
    uint8_t hash_encoding;
    uint8_t max_magn; //
    //uint32_t high_mask; // has hashmag and sizemag, and bit no 1 tells if 2 fields are used

    // 1 bit: data field uses 2nd field
    // 2 bit: data16 is used
    // 4 bit: data32 is used

    // the hashmap has 2 type of entries, one for key->socket, socket->key
    // each have their own 1 hash range
public:
    class ConstIterator{
        public:
        uint32_t cur_value;
        uint32_t nvac;
        uint32_t offset;
        unsigned long mask;
        unsigned long high_mask;
        const Sparsity& target;
        ConstIterator(const Sparsity& trg): target(trg){}
        operator bool ();
        bool operator++(int);
        uint32_t operator()()const{return offset;}
        const uint32_t* operator->(){return &cur_value;}
        uint32_t operator*(){return cur_value;}
    };

    Sparsity();
    ~Sparsity();
    ConstIterator getIterator()const{return ConstIterator(*this);}

    Sparsity& toMemfree();
    Sparsity& toZero(){return this->toMemfree();}
    Sparsity& toOne(){ if (asize != 0) setDense(asize); return *this;}
    Sparsity& setDense(uint32_t length){this->toMemfree(); dsize = asize = length; hash_encoding =0; return *this;}

    uint32_t operator[](uint32_t index)const; // gives index of

    uint32_t genAlias();
    uint32_t getSize() const {return asize;}

    uint32_t find(uint32_t alias)const;
    uint32_t addEntry(uint32_t alias);

    void rehash(uint32_t maxalias, uint32_t maxsize,uint8_t hashsize_mag);
    void show(FILE*f =stdout, int level = 0) const;
};

template<class D>
class AsyncReadHandle{
public:
	std::atomic<uint32_t>* gate;
	uint32_t value;
	D* target;
	AsyncReadHandle(const AsyncReadHandle&)=delete; AsyncReadHandle& operator=(const AsyncReadHandle&)=delete;
    AsyncReadHandle(AsyncReadHandle&&o): gate(o.gate), value(o.value), target(o.target){o.value =0;}
    AsyncReadHandle& operator=(AsyncReadHandle&&o){gate = o.gate; value = o.value; target = o.target; o.value =0;  return *this;}
	const D* operator->(){if (target == NULL){printf("booo!\n");} return target;}
    D operator*(){return *target;}
};

// individual elements in hashmap can be locked
// such elements should be "big enough" to justify locking





template<class K, class D, class HF>
class AsyncHashmap{
	std::atomic<uint32_t> chunkIO;
    void rehash(uint32_t _hsize);
    uint32_t find_routine(const K &key) const;

    uint32_t try_create_routine(const K &key, uint32_t &exitcode);
    uint32_t block_create_routine(const K &key, uint32_t &exitcode);
    //uint32_t find_create_routine(const K &key); // (blocking)
    void moveEntry_routine(uint32_t from, uint32_t to); // finalmutex & permissions obtained
    void delayed_remove_routine(const K &key, std::unique_lock<std::mutex> &lck);
public:
    uint32_t asize;
    uint32_t hsize;
	KeyElem<K, D>* darray;
    uint32_t* dahash;


    AsyncHashmap():asize(0){}
    ~AsyncHashmap(){if (asize != 0) {delete[](darray); delete[](dahash);}}
    AsyncHashmap(const AsyncHashmap<K,D,HF>&)=delete;
    AsyncHashmap(AsyncHashmap<K,D,HF>&&)=delete;
    AsyncHashmap& operator=(const AsyncHashmap<K,D,HF>&)=delete;
    AsyncHashmap& operator=(AsyncHashmap<K,D,HF>&&)=delete;


    // insert/remove   lock for read/write
    class WriteAccess{
        uint32_t try_write_routine(const K &key);
        uint32_t block_routine(const K &key); // returns immediately if object is not found
        friend class AsyncHashmap<K,D,HF>;
    public:
        AsyncHashmap<K,D,HF> &target;
        D* elem;
        uint32_t exitcode;
        WriteAccess(AsyncHashmap<K,D,HF> &_target): target(_target){}
        ~WriteAccess();
        WriteAccess(const WriteAccess&)=delete;
        WriteAccess(WriteAccess&&o): target(o.target), elem(o.elem), exitcode(o.exitcode){o.exitcode =0;}
        WriteAccess& operator=(const WriteAccess&)= delete;
        WriteAccess& operator=(WriteAccess&& o){exitcode = o.exitcode; o.exitcode =0; elem = o.elem; target = o.target; return *this;}
        void unlock(); // preemptively unlock target prior to destruction

        operator bool()const {return exitcode != 0;}
        D* operator->(){return elem;}
        D& operator*(){return *elem;}
        void remove();
    };

    class ReadAccess{
        bool find_only_routine(const K &key);
        uint32_t find_read_syncroutine(const K &key);
        uint32_t block_routine(const K &key); // returns immediately if object is not found
        friend class AsyncHashmap<K,D,HF>;
    public:
        AsyncHashmap<K,D,HF> &target;
        const D* elem;
        uint32_t exitcode;
        ReadAccess(const AsyncHashmap<K,D,HF> &_target): target(const_cast<AsyncHashmap<K,D,HF>&>(_target)){}
        ~ReadAccess();
        ReadAccess(const ReadAccess&)=delete; ReadAccess& operator=(const ReadAccess&)= delete;
        ReadAccess(ReadAccess&&o): target(o.target), elem(o.elem), exitcode(o.exitcode){o.exitcode =0;}
        ReadAccess& operator=(ReadAccess&&o){exitcode = o.exitcode; o.exitcode =0; elem = o.elem; target = o.target; return *this;}
        void unlock(); // preemptively unlock target prior to destruction

        operator bool()const {return exitcode != 0;}
        const D* operator->(){
            if (elem == NULL){
                printf("booo!\n");
            }
            return elem;
        }
        D operator*(){return *elem;}
        void remove(); // this bypass "const", use with care
    };

    class ConstIterator{
        friend class AsyncHashmap<K,D,HF>;
    public:
        uint32_t i;
        AsyncHashmap<K,D,HF> &target;
        uint32_t exitcode;
        ConstIterator(const AsyncHashmap<K,D,HF>& trg): target(const_cast<AsyncHashmap<K,D,HF>&>(trg)), exitcode(0){}
        ConstIterator(const ConstIterator& other)=delete; ConstIterator& operator=(const ConstIterator& other)=delete;
        ConstIterator(ConstIterator&&o): target(o.target), exitcode(o.exitcode){o.exitcode =0;}
        ConstIterator& operator=(ConstIterator&&o){target = o.targetl;exitcode = o.exitcode; o.exitcode =0;}
        ~ConstIterator();
        operator bool ();
        bool operator++(int);
		uint32_t getOffset()const { return i;}
        K operator()()const{return target.darray[i].k;}
        const D* operator->() const{return &(target.darray[i].d);}
        D operator*()const{return target.darray[i].d;}
	};

    ConstIterator getIterator() const {return ConstIterator(*this);}
    //ReadAccess operator()(const K &key); // find (blocking: can only be used by thread owning the structure)
    uint32_t queryChunkIO(){return chunkIO.fetch_add(0);}

    D& createEntry_sync(const K key);
    uint32_t find_sync(const K key)const;
inline D& syncAccess(const K key){return darray[this->find_sync(key)].d;}
inline const D& syncAccess(const K key) const{return darray[this->find_sync(key)].d;}



    // DONE
    void insert(const K kay, D && newelem);
    WriteAccess operator[](const K &key); // find Or create, (blocking: can only be used by thread owning the structure)
inline ReadAccess operator[](const K &key) const{return this->read(key);}
    ReadAccess read(const K &key) const; // if element does not exist, return invalid access
    ReadAccess await(const K &key, atomic<int32_t>* exit_signal) const; // if element does not exist, BLOCKS until created by other thread!

    uint32_t getSize()const{return asize;}
    bool doesContain(const K &key) const;
    ReadAccess tryRead(const K &key) const;
    WriteAccess tryWrite(const K &key);
    WriteAccess tryCreate(const K &key);

    void remove(const K &key);

    // DEBUG
    uint32_t getChunk_please()const{return chunkIO;}
    void show(FILE *f= stdout, int level =0)const;


    void testThisDataStructure(uint32_t nbthreads=32, uint32_t nbitems=17000, uint32_t noisythread = 0xFFFFFFFF); // does some test, insert/removes

};

template<class C, class HF>
class AsyncHashmap<C, void, HF>{
    std::atomic<uint32_t> chunkIO;
    std::atomic<uint32_t> angryIO;
    void rehash(uint32_t _hsize);
    uint32_t find_routine(const typename ExCo<C>::INDEX_TYPE &key) const;

    uint32_t try_create_routine(const typename ExCo<C>::INDEX_TYPE &key, uint32_t &exitcode);
    uint32_t block_create_routine(const typename ExCo<C>::INDEX_TYPE &key, uint32_t &exitcode);
    //uint32_t find_create_routine(const typename ExCo<C>::INDEX_TYPE &key); // (blocking)
    void moveEntry_routine(uint32_t from, uint32_t to); // finalmutex & permissions obtained
    void delayed_remove_routine(const typename ExCo<C>::INDEX_TYPE &key, std::unique_lock<std::mutex> &lck);
public:
    uint32_t asize;
    uint32_t hsize;
	C* darray;
    uint32_t* dahash;


    AsyncHashmap():asize(0){}
    ~AsyncHashmap(){if (asize != 0) {delete[](darray); delete[](dahash);}}
    AsyncHashmap(const AsyncHashmap<C, void, HF>&)=delete;
    AsyncHashmap(AsyncHashmap<C, void, HF>&&)=delete;
    AsyncHashmap& operator=(const AsyncHashmap<C, void, HF>&)=delete;
    AsyncHashmap& operator=(AsyncHashmap<C, void, HF>&&)=delete;


    // insert/remove   lock for read/write
    class WriteAccess{
        uint32_t try_write_routine(const typename ExCo<C>::INDEX_TYPE &key);
        uint32_t block_routine(const typename ExCo<C>::INDEX_TYPE &key); // returns immediately if object is not found
        friend class AsyncHashmap<C, void, HF>;
    public:
        AsyncHashmap<C, void, HF> &target;
        C* elem;
        uint32_t exitcode;
        WriteAccess(AsyncHashmap<C, void, HF> &_target): target(_target){}
        ~WriteAccess();
        WriteAccess(const WriteAccess&)=delete;
        WriteAccess(WriteAccess&&o): target(o.target), elem(o.elem), exitcode(o.exitcode){o.exitcode =0;}
        WriteAccess& operator=(const WriteAccess&)= delete;
        WriteAccess& operator=(WriteAccess&& o){exitcode = o.exitcode; o.exitcode =0; elem = o.elem; target = o.target; return *this;}
        void unlock(); // preemptively unlock target prior to destruction

        operator bool()const {return exitcode != 0;}
        C* operator->(){return elem;}
        C& operator*(){return *elem;}
        void remove();
    };

    class ReadAccess{
        bool find_only_routine(const typename ExCo<C>::INDEX_TYPE &key);
        uint32_t find_read_syncroutine(const typename ExCo<C>::INDEX_TYPE &key);
        uint32_t block_routine(const typename ExCo<C>::INDEX_TYPE &key); // returns immediately if object is not found
        friend class AsyncHashmap<C, void, HF>;
    public:
        AsyncHashmap<C, void, HF> &target;
        const C* elem;
        uint32_t exitcode;
        ReadAccess(const AsyncHashmap<C, void, HF> &_target): target(const_cast<AsyncHashmap<C, void, HF>&>(_target)){}
        ~ReadAccess();
        ReadAccess(const ReadAccess&)=delete; ReadAccess& operator=(const ReadAccess&)= delete;
        ReadAccess(ReadAccess&&o): target(o.target), elem(o.elem), exitcode(o.exitcode){o.exitcode =0;}
        ReadAccess& operator=(ReadAccess&&o){exitcode = o.exitcode; o.exitcode =0; elem = o.elem; target = o.target; return *this;}
        void unlock(); // preemptively unlock target prior to destruction

        operator bool()const {return exitcode != 0;}
        const C* operator->(){
            if (elem == NULL){
                printf("booo!\n");
            }
            return elem;
        }
        C operator*(){return *elem;}
        void remove(); // this bypass "const", use with care
    };

    class ConstIterator{
        friend class AsyncHashmap<C, void, HF>;
    public:
        uint32_t i;
        AsyncHashmap<C, void, HF> &target;
        uint32_t exitcode;
        ConstIterator(const AsyncHashmap<C, void, HF>& trg): target(const_cast<AsyncHashmap<C, void, HF>&>(trg)), exitcode(0){}
        ConstIterator(const ConstIterator& other)=delete; ConstIterator& operator=(const ConstIterator& other)=delete;
        ConstIterator(ConstIterator&&o): target(o.target), exitcode(o.exitcode){o.exitcode =0;}
        ConstIterator& operator=(ConstIterator&&o){target = o.targetl;exitcode = o.exitcode; o.exitcode =0;}
        ~ConstIterator();
        operator bool ();
        bool operator++(int);
		uint32_t getOffset()const { return i;}
        typename ExCo<C>::INDEX_TYPE operator()()const{return ExOp::getIndex(target.darray[i]);}
        const C* operator->() const{return &(target.darray[i]);}
        C operator*()const{return target.darray[i];}
	};

    ConstIterator getIterator() const {return ConstIterator(*this);}
    //ReadAccess operator()(const typename ExCo<C>::INDEX_TYPE &key); // find (blocking: can only be used by thread owning the structure)
    uint32_t queryChunkIO(){return chunkIO.fetch_add(0);}


    // DONE
    void insert(C && newelem);

    C& createEntry_sync(const typename ExCo<C>::INDEX_TYPE key);
    uint32_t find_sync(const typename ExCo<C>::INDEX_TYPE key)const;
inline C& syncAccess(const typename ExCo<C>::INDEX_TYPE key){return darray[this->find_sync(key)];}
inline const C& syncAccess(const typename ExCo<C>::INDEX_TYPE key) const{return darray[this->find_sync(key)];}

    WriteAccess operator[](const typename ExCo<C>::INDEX_TYPE &key); // find Or create, (blocking: can only be used by thread owning the structure)
inline ReadAccess operator[](const typename ExCo<C>::INDEX_TYPE &key) const{return this->read(key);}
    ReadAccess read(const typename ExCo<C>::INDEX_TYPE &key) const; // if element does not exist, return invalid access
    ReadAccess await(const typename ExCo<C>::INDEX_TYPE &key, atomic<int32_t>* exit_signal) const; // if element does not exist, BLOCKS until created by other thread!

    uint32_t getSize()const{return asize;}
    bool doesContain(const typename ExCo<C>::INDEX_TYPE &key) const;
    ReadAccess tryRead(const typename ExCo<C>::INDEX_TYPE &key) const;
    WriteAccess tryWrite(const typename ExCo<C>::INDEX_TYPE &key);
    WriteAccess tryCreate(const typename ExCo<C>::INDEX_TYPE &key);

    void remove(const typename ExCo<C>::INDEX_TYPE &key);



    //uint32_t find_routine(typename ExCo<C>::INDEX_TYPE &key);

    //WriteAccess addEntry(typename ExCo<C>::INDEX_TYPE &key);
    //void removeEntry(typename ExCo<C>::INDEX_TYPE &key);

    // DEBUG
    uint32_t getChunk_please()const{return chunkIO;}
    void show(FILE *f= stdout, int level =0)const;


    void testThisDataStructure(uint32_t nbthreads=32, uint32_t nbitems=17000, uint32_t noisythread = 0xFFFFFFFF); // does some test, insert/removes

};
template<class Key, class Data, class Category, class HashFnc, class Category_HashFnc >
class CategoryHashmap{
// entries are reached either by their key, as a hashmap, or sequentially within sets of entries grouped in categories

    typedef typename ExCo<Key>::TYPE KeyType;
    unsigned int hashpos(unsigned int seed) const;
    void remove_links_routine(uint32_t ite, uint32_t catite);

	void move_up_routine(unsigned int from, unsigned int to);
    void swap_and_pop(unsigned int ite);

    unsigned int enlargeback_routine(unsigned int catite);
    void memmove_routine(uint32_t target, uint32_t source);

    unsigned int getNewSlotForCategory_routine(const Category&);

    public:
    typedef unsigned int ITERATOR_TYPE;

    uint32_t asize;
    pair< KeyElem<Key, Data> , unsigned int >* heap;


    typedef myHashmap<Key, void, HashFnc> SAFETYPE;
    unsigned int* hash;

    unsigned char hash_mag;

    myHashmap<Category, pair<unsigned int, unsigned int >, Category_HashFnc > categories;

    uint32_t nb_categoryless;

    inline int getSize()const{return asize & 0x7FFFFFF;}

    CategoryHashmap(): asize(0){}
    ~CategoryHashmap();
    unsigned int find(const Key &) const;
    Data& deref(unsigned int ite){return heap[ite].first.d;}
    const Data& deref(unsigned int ite) const {return heap[ite].first.d;}

    Key& keypreserving_deref(unsigned int ite){return heap[ite].first.k;}
    const Key& deref_key(unsigned int ite) const {return heap[ite].first.k;}


    //  void insert(const Key &, const Data&);
    void insert(const Key &, const Data&, const Category& );

    uint32_t findCategory(uint32_t ite)const;
    Category deref_category(uint32_t ite)const;
    Category getCategory(const Key &)const;


    Data operator[](const Key &) const; // assumes entry exists, use "find" otherwise
    Data& operator[](const Key &); // creates entry if it does not exist

    CategoryHashmap<Key,Data,Category,HashFnc,Category_HashFnc>& toMemfree(){categories.toMemfree(); if (hash_mag != 0) {delete[](hash); delete[](heap);} hash_mag = 0; return(*this); }
    CategoryHashmap<Key,Data,Category,HashFnc,Category_HashFnc>& toMemmove(CategoryHashmap<Key,Data,Category,HashFnc,Category_HashFnc>& other);
    CategoryHashmap<Key,Data,Category,HashFnc,Category_HashFnc>& operator=(const CategoryHashmap<Key,Data,Category,HashFnc,Category_HashFnc>& other);

    void getCategoryIndexes(const Category&, unsigned int& first, unsigned int& size) const;
    unsigned int getSizeWithinCategory(const Category& what) const{unsigned int fout[2]; this->getCategoryIndexes(what,fout[0],fout[1]); return fout[1];}
    unsigned int getFirstWithinCategory(const Category& what) const{unsigned int fout[2]; this->getCategoryIndexes(what,fout[0],fout[1]); return fout[0];}




    void erase(const Key &);
    void erase_from_iterator(unsigned int ite);
    void setCategory(const Key& , const Category& );

    void rehash(unsigned char _new_mag);

    void show(FILE*f = stdout, int level=0) const;
};
template<class Key, class Category, class HashFnc, class Category_HashFnc >
class CategoryHashmap<Key, void, Category, HashFnc, Category_HashFnc>{
 typedef typename ExCo<Key>::TYPE KeyType;
    unsigned int hashpos(unsigned int seed) const;

    void remove_links_routine(uint32_t ite, uint32_t catite);
	void move_up_routine(unsigned int from, unsigned int to);
    void swap_and_pop(unsigned int ite);

    void memmove_routine(uint32_t target, uint32_t source);
    unsigned int enlargeback_routine(unsigned int catite);

    public:
    typedef unsigned int ITERATOR_TYPE;

    uint32_t nb_categoryless;
    uint32_t asize;
    pair< Key , unsigned int >* heap;

    typedef myHashmap<Key, void, HashFnc> SAFETYPE;
    unsigned int* hash;

    unsigned char hash_mag;

    myHashmap<Category, pair<unsigned int, unsigned int >, Category_HashFnc > categories;


    class QueryIterator{
		public:
		const CategoryHashmap<Key, void, Category, HashFnc, Category_HashFnc>& target;
		unsigned int ite;
		unsigned int nb;
		QueryIterator(const CategoryHashmap<Key, void, Category, HashFnc, Category_HashFnc>& _target):target(_target){}
		const Key* operator->()const;
        bool findFirst( const Category &cat);
        bool findNext();
	};


	QueryIterator mkQueryIterator() const{return CategoryHashmap<Key, void, Category, HashFnc, Category_HashFnc>::QueryIterator(*this);}

    inline int getSize()const{return asize & 0x7FFFFFF;}

    CategoryHashmap(): asize(0){}
    ~CategoryHashmap();
    unsigned int find(const typename ExCo<Key>::INDEX_TYPE &) const;

    Key& keypreserving_deref(unsigned int ite){return heap[ite].first;}
    const Key& deref(unsigned int ite) const {return heap[ite].first;}


    Key& keypreserving_deref(QueryIterator &query){/*if (query.target != this) exit(1); */return heap[query.ite].first;}


    unsigned int insert(const Key &, const Category&); // return iterator to entry

    uint32_t findCategory(uint32_t ite)const;
    Category deref_category(uint32_t ite)const;
    Category getCategory(const typename ExCo<Key>::INDEX_TYPE &)const;

    Key operator[](const typename ExCo<Key>::INDEX_TYPE  &) const; // assumes entry exists, use "find" otherwise
    Key& operator[](const typename ExCo<Key>::INDEX_TYPE  &); // creates entry with default category if it does not exist

    CategoryHashmap<Key,void,Category,HashFnc,Category_HashFnc>& toMemfree(){categories.toMemfree(); if (hash_mag != 0) {delete[](hash); delete[](heap);} hash_mag = 0; return(*this); }
    CategoryHashmap<Key,void,Category,HashFnc,Category_HashFnc>& toMemmove(CategoryHashmap<Key,void,Category,HashFnc,Category_HashFnc>& other);
    CategoryHashmap<Key,void,Category,HashFnc,Category_HashFnc>& operator=(const CategoryHashmap<Key,void,Category,HashFnc,Category_HashFnc>& other);


    void getCategoryIndexes(const Category&, unsigned int& first, unsigned int& size) const;
    unsigned int getSizeWithinCategory(const Category& what) const{unsigned int fout[2]; this->getCategoryIndexes(what,fout[0],fout[1]); return fout[1];}
    unsigned int getFirstWithinCategory(const Category& what) const{unsigned int fout[2]; this->getCategoryIndexes(what,fout[0],fout[1]); return fout[0];}




    void erase(const typename ExCo<Key>::INDEX_TYPE &);
    void erase_from_iterator(unsigned int ite);

    void rehash(unsigned char _new_mag);

    void show(FILE*f = stdout, int level=0) const;
};

// associations are given an uint32_t alias, which can alternatively be used to access the data, and is required to remove entries (unless Key is also uniquely determining)
// preserve order of  elements inserted with identical keys (insert/(query first) in front by default, as a stack)
template<class Key, class Data>
class AliasedHashmap{
    void rehash_routine(uint32_t _new_mag);

    uint32_t prepare_slot_routine(uint32_t& alias, uint8_t protected_bits);
	uint32_t asize;
    public:
	uint32_t hash_mask;
    typedef unsigned int ITERATOR_TYPE;
	pair< KeyElem< Key, Data>, Tuple<uint32_t, 2u> >* heap; //key, data, next, alias-next
    uint32_t* main_hash; // 2^n (uint32_t) main hash, 2^n (uint32_t) secondaty hash, 1 (uint32_t) emptypos, 2^n bits stream
	uint32_t emptypos;

    template<int ISCONST> class Iterator{
    public:
        typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF target;
        uint32_t cur;
        Iterator(typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF trg): target(trg){}
        operator bool (){if (target.asize == 0) return false; cur = 0; if ((target.heap[0].second[0]+1) & 0x80000000) return (*this)++; return true;}
        bool operator++(int){ for(cur++;cur <= target.hash_mask; cur++) if (((target.heap[cur].second[0]+1)  & 0x80000000) == 0) return true; return false;}
		Key operator()() const{return target.heap[cur].first.k;}
        U32_Alias getAlias() const{return target.getAliasAt_routine(cur);}
        typename MetaType<Data*,ISCONST>::IS_CONST operator->(){return &(target.heap[cur].first.d);}
        typename MetaType<Data,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.heap[cur].first.d;}
	};

    template<int ISCONST> class SlotIterator{
    public:
    	typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF target;
        const Key key;
        uint32_t cur;
        SlotIterator(typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF trg, const Key& _k): target(trg), key(_k){}
        operator bool (){if ((target.asize == 0)||((cur = target.main_hash[ExOp::mkHashValue(key) & target.hash_mask]) == 0xFFFFFFFF)) return false; return (target.heap[cur].first.k != key) ? (*this)++ : true;}
        bool operator++(int){ do{cur = target.heap[cur].second[0]; if (cur == 0xFFFFFFFF) return false;}while(target.heap[cur].first.k != key); return true;}
        Key operator()() const{return key;}
        U32_Alias getAlias() const{return target.getAliasAt_routine(cur);}
        typename MetaType<Data*,ISCONST>::IS_CONST operator->(){return &(target.heap[cur].first.d);}
        typename MetaType<Data,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.heap[cur].first.d;}
        void show(FILE* f =stdout, int level=0)const{ printf("Valid Slot iterator with key: "); ExOp::show(key,f, 0);}
	};

	template<int ISCONST> class FindIterator{ public:
        typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF target;
        uint32_t cur;
        FindIterator(typename MetaType<AliasedHashmap<Key, Data>,ISCONST>::IS_CONST_REF trg, uint32_t _cur): target(trg), cur(_cur){}
        operator bool (){return (cur != 0xFFFFFFFF);}
		Key operator()() const{return target.heap[cur].first.k;}
        U32_Alias getAlias() const{return target.getAliasAt_routine(cur);}

        FindIterator& operator=(const FindIterator& other){cur = other.cur;return *this;} // Warning, assumes the target is the same, use with care

        typename MetaType<Data*,ISCONST>::IS_CONST operator->(){return &(target.heap[cur].first.d);}
        typename MetaType<Data,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.heap[cur].first.d;}
        void remove(){ target.erase_from_iterator(cur);} // Iterator is invalid afterwards
	};

    AliasedHashmap(): asize(0){}
    AliasedHashmap(const AliasedHashmap<Key, Data>& other);
	AliasedHashmap(AliasedHashmap<Key, Data>&& o) : asize(o.asize), hash_mask(o.hash_mask),heap(o.heap), main_hash(o.main_hash), emptypos(o.emptypos){o.asize=0;}
    ~AliasedHashmap();

    AliasedHashmap& operator=(const AliasedHashmap<Key, Data>& other);
    AliasedHashmap& operator=(AliasedHashmap<Key, Data>&& other);

    Data& createEntry(const Key &, uint32_t &alias, bool insert_at_end = false, uint8_t protected_bits = 0);
	Data& createEntryWithAlias(const Key &key, uint32_t alias, bool insert_at_end = false){return createEntry(key,alias,insert_at_end,32);}
	Data& createEntryAfter(const U32_Alias &, uint32_t &alias, bool insert_before_instead = false, uint8_t protected_bits = 0);
	Data& createEntryAfterWithAlias(const U32_Alias &prev, uint32_t alias, bool insert_at_end = false){return createEntryAfter(prev,alias,insert_at_end,32);}

	Iterator<0> mkIterator() {return Iterator<0>(*this);}
	Iterator<1> mkIterator()const {return Iterator<1>(*this);}
	SlotIterator<0> mkIterator(const Key &k){return SlotIterator<0>(*this, k);}
	SlotIterator<1> mkIterator(const Key &k)const{return SlotIterator<1>(*this, k);}
	// SYNONYMOUS TO
	Iterator<0> operator()() {return Iterator<0>(*this);}
	Iterator<1> operator()()const {return Iterator<1>(*this);}
	SlotIterator<0> operator()(const Key &k){return SlotIterator<0>(*this, k);}
	SlotIterator<1> operator()(const Key &k)const{return SlotIterator<1>(*this, k);}
	// SYNONYMOUS TO

	SlotIterator<0> mkIteratorFromFirstAlias(uint32_t k){return SlotIterator<0>(*this, heap[this->find(U32_Alias(k))].first.k);}
	SlotIterator<1> mkIteratorFromFirstAlias(uint32_t k)const{return SlotIterator<1>(*this, heap[this->find(U32_Alias(k))].first.k);}


	uint32_t getSize()const{return asize;}
    uint32_t find(const Key &) const;
    uint32_t find(const U32_Alias &) const;
    FindIterator<0> operator()(U32_Alias alias){return FindIterator<0>(*this, const_cast<const AliasedHashmap<Key,Data>* >(this)->find(alias));}
    FindIterator<1> operator()(U32_Alias alias) const{return FindIterator<1>(*this, this->find(alias));}

    AliasedHashmap& toMemfree();
    bool hasAlias(uint32_t) const;

    U32_Alias getAlias(const Key &key) const {uint32_t r = this->find(key); return U32_Alias(r == 0xFFFFFFFF ? 0 : getAliasAt_routine(r));}
    Data& operator[](const Key &);
    Data& operator[](const U32_Alias &);

    Data& deref(unsigned int ite){return heap[ite].first.d;}
	const Data& deref(unsigned int ite)const{return heap[ite].first.d;}
    const Key& deref_key(unsigned int ite) const {return heap[ite].first.k;}
    const Key& deref_key(const U32_Alias &alias) const {uint32_t ite = this->find(alias); if (ite == 0xFFFFFFFF) throw std::domain_error("Alias-assessed element does not exist");	return heap[ite].first.k;}


    uint32_t deref_alias(unsigned int ite) const {return ite == 0xFFFFFFFF ? 0 : getAliasAt_routine(ite);}
	uint32_t getAliasAt_routine(uint32_t pos) const;
    void erase(const Key &k){erase_from_iterator(find(k));}
    void eraseUsingAlias(uint32_t k){erase_from_iterator(find(k));}


    void erase_from_iterator(uint32_t ite, uint32_t * ite_to_maintain = NULL); // maintain todo...

    Key frontKey(const U32_Alias& b) const {return heap[find(b)].first.k;}

    ERRCODE save(FILE*f)const;
    ERRCODE load(FILE*f);

    void show(FILE*f = stdout, int level=0)const ;
    void test(FILE*f = stdout);
};

template<class Key>
class AliasedHashmap<Key, void>{
    void rehash_routine(uint32_t _new_mag);

	uint32_t asize;
	uint32_t prepare_slot_routine(uint32_t& alias, uint8_t protected_bits);
	template<class K2, class D2> friend class StructuredHash;
    public:
	uint32_t hash_mask;

    typedef unsigned int ITERATOR_TYPE;
	pair< Key, Tuple<uint32_t, 2u> >* heap; //key, data, next, alias-next
    uint32_t* main_hash; // 2^n (uint32_t) main hash, 2^n (uint32_t) secondaty hash, 1 (uint32_t) emptypos, 2^n bits stream
	uint32_t emptypos;


	template<int ISCONST> class SlotIterator{
    public:
    	typename MetaType<AliasedHashmap<Key, void>,ISCONST>::IS_CONST_REF target;
        const typename ExCo<Key>::INDEX_TYPE key;
        uint32_t cur;
        SlotIterator(typename MetaType<AliasedHashmap<Key, void>,ISCONST>::IS_CONST_REF trg, const typename ExCo<Key>::INDEX_TYPE & _k): target(trg), key(_k){}
        operator bool (){if ((target.asize == 0)||((cur = target.main_hash[ExOp::mkHashValue(key) & target.hash_mask]) == 0xFFFFFFFF)) return false; return (ExOp::getIndex(target.heap[cur].first) != key) ? (*this)++ : true;}
        bool operator++(int){ do{cur = target.heap[cur].second[0]; if (cur == 0xFFFFFFFF) return false;}while(ExOp::getIndex(target.heap[cur].first) != key); return true;}
        uint32_t operator()()const{return target.getAliasAt_routine(cur);}
        typename MetaType<Key*,ISCONST>::IS_CONST operator->(){return &(target.heap[cur].first);}
        typename MetaType<Key,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.heap[cur].first;}

        //typename std::enable_if<ISCONST == 0, SlotIterator&>::type changeKeyForNext(typename ExCo<Key>::INDEX_TYPE newkey, uint32_t nbnexts=1);
        void show(FILE* f =stdout, int level=0)const{ printf("Valid Slot iterator with key: "); ExOp::show(key,f, 0);}
	};

    AliasedHashmap(): asize(0){}
    AliasedHashmap(const AliasedHashmap<Key, void>& other);
    AliasedHashmap(AliasedHashmap<Key, void>&& o) : asize(o.asize), hash_mask(o.hash_mask),heap(o.heap), main_hash(o.main_hash), emptypos(o.emptypos){o.asize=0;}
    ~AliasedHashmap();
    AliasedHashmap& operator=(const AliasedHashmap<Key, void>& other);
    AliasedHashmap& operator=(AliasedHashmap<Key, void>&& other);


	SlotIterator<0> mkIterator(const typename ExCo<Key>::INDEX_TYPE &k){return SlotIterator<0>(*this, k);}
	SlotIterator<1> mkIterator(const typename ExCo<Key>::INDEX_TYPE &k)const{return SlotIterator<1>(*this, k);}



    Key& createEntry(const Key &, uint32_t &alias, bool insert_at_end = false, uint8_t protected_bits = 0);
	Key& createEntryWithAlias(const Key &key, uint32_t alias, bool insert_at_end = false){return createEntry(key,alias,insert_at_end,32);}
	Key& createEntryAfter(const U32_Alias &, uint32_t &alias, bool insert_before_instead = false, uint8_t protected_bits = 0);
	Key& createEntryAfterWithAlias(const U32_Alias &prev, uint32_t alias, bool insert_at_end = false){return createEntryAfter(prev,alias,insert_at_end,32);}

	AliasedHashmap<Key,void>& moveNextEntry(uint32_t ite, const Key& new_key, int pos = 0);

	AliasedHashmap<Key,void>& changeKeyForNextAfter(uint32_t ite, const Key &new_key, uint32_t nbitem);
	//AliasedHashmap<Key,void>& changeKeyForFirst(const Key &key, const Key &new_key, uint32_t nbitem);

	uint32_t getAliasAt_routine(uint32_t pos) const;

    unsigned int find(const typename ExCo<Key>::INDEX_TYPE &) const;
    unsigned int find(const U32_Alias &) const;

	U32_Alias getAlias(const typename ExCo<Key>::INDEX_TYPE &key) {uint32_t r = this->find(key); return U32_Alias(r == 0xFFFFFFFF ? 0 : getAliasAt_routine(r));}

    Key& operator[](const typename ExCo<Key>::INDEX_TYPE &);
    Key& operator[](const U32_Alias &);

	Key& deref(unsigned int ite){return heap[ite].first;}
	uint32_t deref_alias(unsigned int ite) const {return ite == 0xFFFFFFFF ? 0 : getAliasAt_routine(ite);}

	bool hasAlias(uint32_t alias)const;
	void erase_from_iterator(uint32_t ite);
	AliasedHashmap<Key,void>& toMemfree();

    ERRCODE save(FILE*f) const;
    ERRCODE load(FILE*f);

    void show(FILE*f = stdout, int level=0)const ;
    void test(FILE*f = stdout);
};

// Structure maintaining tree-like structure of nodes
template<class Key, class Data>
class StructuredHash{
public: // relations are either "child" or "next", contiguous "next" are on tiles, contiguous in memory
	AliasedHashmap<Key,Data> data;
	AliasedHashmap<uint32_t> links;
	myHashmap<uint32_t> voidparents; // allows unstructured forests with no parent

	StructuredHash();
    StructuredHash& insertElem(const Key& key, const Data& data);
    StructuredHash& insertElemUnder(const Key& key, const Data& data, const Key& parent);
    StructuredHash& insertElemAfter(const Key& key, const Data& data, const Key& prev);

    // takes X contiguous entries after prec, and make next child of a newly created node;
    StructuredHash& insertParentNodeAt(const Key& key, const Data& data, U32_Alias first, uint32_t nbsub);

    StructuredHash& setChild(const Key& child_key, const Key& parent_key);
    StructuredHash& setNext(const Key& left_key, const Key& right_key);
    StructuredHash& addChild(const Key& child_key, const Key& parent_key);
    StructuredHash& addNext(const Key& left_key, const Key& right_key);

	template<int ISCONST> class TreeIterator{
    public:
    	typename MetaType<StructuredHash<Key, Data>,ISCONST>::IS_CONST_REF target;
        uint32_t cur;
        uint32_t dcur;
        TreeIterator(typename MetaType<StructuredHash<Key, Data>,ISCONST>::IS_CONST_REF trg): target(trg){}
        TreeIterator(typename MetaType<StructuredHash<Key, Data>,ISCONST>::IS_CONST_REF trg, const Key & _k): target(trg){if ((dcur = target.data.main_hash[ExOp::mkHashValue(_k) & target.data.hash_mask]) == 0xFFFFFFFF) return; while(ExOp::getIndex(target.data.heap[dcur].first) != _k){dcur = target.data.heap[dcur].second[0]; if (dcur == 0xFFFFFFFF) return;} U32_Alias daalias(target.data.deref_alias(dcur)); cur = target.links.find(daalias);}
        TreeIterator(typename MetaType<StructuredHash<Key, Data>,ISCONST>::IS_CONST_REF trg, const U32_Alias & _alias): target(trg){dcur = target.data.find(_alias);cur = target.links.find(_alias);}
        operator bool (){return (dcur != 0xFFFFFFFF);}
        bool hasParent()const{return(target.links.heap[cur].first != 0);}
        uint32_t getParentAlias()const{return(target.links.heap[cur].first);}
        bool next(){ uint32_t dapar = target.links.heap[cur].first;   do{cur = target.links.heap[cur].second[0]; if (cur == 0xFFFFFFFF) {dcur = 0xFFFFFFFF; return false;}}while(ExOp::getIndex(target.links.heap[cur].first) != dapar); U32_Alias daalias(target.links.deref_alias(cur)); dcur = target.data.find(daalias); return true;}
		bool preordernext() {
			uint32_t ncur;
			if ((ncur = target.links.find(target.data.getAliasAt_routine(dcur))) != 0xFFFFFFFF) cur = ncur;
			else{
				ncur = target.links.heap[cur].second[0];
				while(ncur != 0xFFFFFFFF) {
					//printf("%X\t%X\t%X\n", target.links.heap[ncur].first, target.links.heap[cur].first, target.links.getAliasAt_routine(ncur));
					if (target.links.heap[ncur].first == target.links.heap[cur].first) break;
					ncur = target.links.heap[ncur].second[0];
				}
				if (ncur != 0xFFFFFFFF) cur = ncur;
				else return false; // TODO
			}
			dcur = target.data.find(U32_Alias(target.links.getAliasAt_routine(cur)));
		return true;}
		TreeIterator<ISCONST>  mkPreordernext() const { TreeIterator<ISCONST> fout(target);
			if ((fout.cur = target.links.find(target.data.getAliasAt_routine(dcur))) == 0xFFFFFFFF){
				fout.cur = target.links.heap[cur].second[0];
				while(fout.cur != 0xFFFFFFFF) {
					//printf("%X\t%X\t%X\n", target.links.heap[ncur].first, target.links.heap[cur].first, target.links.getAliasAt_routine(ncur));
					if (target.links.heap[fout.cur].first == target.links.heap[cur].first) break;
					fout.cur = target.links.heap[fout.cur].second[0];
				}
				if (fout.cur == 0xFFFFFFFF) {fout.dcur = 0xFFFFFFFF; return fout;}// TODO
			}
			fout.dcur = target.data.find(U32_Alias(target.links.getAliasAt_routine(fout.cur)));
		return fout;}
		TreeIterator<ISCONST> mkLastChild() const { TreeIterator<ISCONST> fout(target);
			uint32_t ite;
			if ((ite = target.links.find(target.links.getAliasAt_routine(cur))) == 0xFFFFFFFF){fout.cur = cur; fout.dcur = dcur; return fout;} // no child, return self
			do{ fout.cur = ite;
			}while((ite = target.links.find(target.links.getAliasAt_routine(fout.cur))) != 0xFFFFFFFF);
			fout.dcur = target.data.find(U32_Alias(target.links.getAliasAt_routine(fout.cur)));
		return fout;}

		bool postordernext() {
			//printf("got %X, alias %X %i:%i\n", cur, target.links.getAliasAt_routine(cur), target.data.heap[dcur].first.k, target.data.heap[dcur].first.d);
			//printf("next is...\n");
			uint32_t ncur = target.links.heap[cur].second[0];
			uint32_t tal;
			while(ncur != 0xFFFFFFFF) {
				//printf("%X\t%X\t%X\n", target.links.heap[ncur].first, target.links.heap[cur].first, target.links.getAliasAt_routine(ncur));
				if (target.links.heap[ncur].first == target.links.heap[cur].first) break;
				ncur = target.links.heap[ncur].second[0];
			}
			if (ncur != 0xFFFFFFFF){
				//printf("found next at %i, alias %X\n", ncur, target.links.getAliasAt_routine(ncur));
				while((cur = target.links.find( tal = target.links.getAliasAt_routine(ncur))) != 0xFFFFFFFF) ncur = cur;
				cur = ncur;
			}else {
				//printf("using parent alias %X, alias %X was last\n", target.links.heap[cur].first, target.links.getAliasAt_routine(cur));
				tal = target.links.heap[cur].first;
				if (tal == 0) return false;
				cur = target.links.find(U32_Alias(tal));
				//printf("found eh? %X == %X\n",tal, target.links.getAliasAt_routine(cur));
			}
			dcur = target.data.find(U32_Alias(tal));
			//printf("got %X, alias %X %i:%i\n", cur, target.links.getAliasAt_routine(cur), target.data.heap[dcur].first.k, target.data.heap[dcur].first.d);
		return true;}

        bool operator++(int){ uint32_t dapar = target.links.heap[cur].first;   do{cur = target.links.heap[cur].second[0]; if (cur == 0xFFFFFFFF) return false;}while(ExOp::getIndex(target.links.heap[cur].first) != dapar); U32_Alias daalias(target.links.deref_alias(cur)); dcur = target.data.find(daalias); return true;}

		TreeIterator operator+(int value) const{ TreeIterator fout(target);
			if (value >= 0){
				fout.cur = cur;
				while (value != 0){
					fout.cur = target.links.heap[fout.cur].second[0];
					if (fout.cur == 0xFFFFFFFF) {fout.dcur = 0xFFFFFFFF; return fout;}
					if (target.links.heap[fout.cur].first == target.links.heap[cur].first) value--;
				}
			}else{ // TODO

			}
			fout.dcur = target.data.find(U32_Alias(target.links.getAliasAt_routine(cur)));
		return fout;}

        U32_Alias getAlias()const{ return target.data.getAliasAt_routine(dcur);}
        Key operator()()const{return target.data.heap[dcur].first.k;}
		typename std::enable_if<ISCONST == 0, uint32_t>::type moveNextForward(uint32_t nbpos = 1);
		typename std::enable_if<ISCONST == 0, uint32_t>::type removeNext();
		typename std::enable_if<ISCONST == 0, uint32_t>::type moveNextUnder(int position = 0);
        typename MetaType<Data*,ISCONST>::IS_CONST operator->(){return &(target.data.heap[dcur].first.d);}
        typename MetaType<Data,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return target.data.heap[dcur].first.d;}
        void show(FILE* f =stdout, int level=0)const{ fprintf(f,"Tree iterator with [key="); ExOp::show(target.data.heap[dcur].first.k ,f, 1); fprintf(f,", data="); ExOp::show(target.data.heap[dcur].first.d ,f, 1); fprintf(f,", alias=%X]\n", target.data.getAliasAt_routine(dcur)); }
	};

//	1->2->3->4
//	1->4->0 3->2->3
//	1->4->3->2->0

//	1->2->3
//	1->3->0 2->2->2
//	1->3->2->0

	TreeIterator<0> mkIterator(const Key &k){return TreeIterator<0>(*this, k);}
	TreeIterator<1> mkIterator(const Key &k)const{return TreeIterator<1>(*this, k);}
	TreeIterator<0> mkIterator(const U32_Alias &alias){return TreeIterator<0>(*this, alias);}
	TreeIterator<1> mkIterator(const U32_Alias &alias)const{return TreeIterator<1>(*this, alias);}



    ERRCODE save(FILE*f)const;
    ERRCODE load(FILE*f);
    void show(FILE*f = stdout, int level=0)const;
};

// hashmap where every entries are tagged with an alias, which is an altenative hashable key to retrieve the key-data pairing

class OptionStruct{
    public:
    myHashmap<string, uint32_t> data;
    OptionStruct();
    OptionStruct(string value);
    OptionStruct(const char* value);
    OptionStruct(const std::nullptr_t value);
    template<int L> OptionStruct(const string values[L]);
    template<typename... VARTIC> OptionStruct(string value, VARTIC... vartic);
    bool has(const char*);
    void populateFlags(const char* str, char sep = ',');

};

template<class IC, class D>
class MaskMap{
public:
    myHashmap<Tuple<IC,0u>, D > hashmap;
    void insert(const Tuple<IC>& key, const D& value);
    void getIntersection(Vector<D> &res, const Tuple<IC>& value, const Tuple<IC>& mask) const; // return all entries with  ((value ^ key) & mask) == 0
};
template<class C>
class SparseTuple : public myHashmap<uint32_t,C>{
public:
    SparseTuple(): myHashmap<uint32_t,C>(){}
    SparseTuple(const SparseTuple<C>& other): myHashmap<uint32_t,C>(dynamic_cast<const myHashmap<uint32_t,C>&>(other)){}
    SparseTuple(SparseTuple<C>&& other): myHashmap<uint32_t,C>(std::move(dynamic_cast<myHashmap<uint32_t,C>&>(other))){}
    SparseTuple<C>& operator=(SparseTuple<C>&& other){*dynamic_cast<myHashmap<uint32_t,C>*>(this) = std::move(dynamic_cast<myHashmap<uint32_t,C>&>(other)); return *this;}
    SparseTuple<C>& operator=(const SparseTuple<C>& other){*dynamic_cast<myHashmap<uint32_t,C>*>(this) = dynamic_cast<const myHashmap<uint32_t,C>&>(other); return *this;}
    template<class O> SparseTuple(const std::vector<O> &other);

    SparseTuple(const myHashmap<uint32_t,C>& other): myHashmap<uint32_t,C>(other){}

    template<class D> C mkInnerMult(const SparseTuple<D>& other)const;
    template<class D> SparseMatrix<C> mkOuterMult(const SparseTuple<D>& other)const;

    SparseTuple<C>& toMemmove(SparseTuple<C>& other){
        dynamic_cast<myHashmap<uint32_t,C>*>(this)->toMemmove(dynamic_cast<myHashmap<uint32_t,C>&>(other)); return *this;
    }
    template<class O> std::vector<O>& wrStdVector(std::vector<O>& fout)const;
    SparseTuple<C>& toMult(const Trianglix<C,0u>& a, const SparseTuple<C>&b); // only fills needed entries
    bool makeMaximumAsFirst(); // return (did change)
    bool hasSameSparsity(const SparseTuple<C> &other)const; // assumes both are sorted

    template<class O> bool operator==(const SparseTuple<O>& other)const;
    template<class O> SparseTuple<C>& operator+=(const SparseTuple<O>& other){if (auto ite = other.mkIterator()) do{(*this)[ite()] += *ite;} while(ite++); return *this;}
    template<class O> SparseTuple<C>& operator-=(const SparseTuple<O>& other){if (auto ite = other.mkIterator()) do{(*this)[ite()] -= *ite;} while(ite++); return *this;}
    template<class O> SparseTuple<C> operator+(const SparseTuple<O>& other)const{SparseTuple<C> fout(*this); if (auto ite = other.mkIterator()) do{fout[ite()] += *ite;} while(ite++); return fout;}
    template<class O> SparseTuple<C> operator-(const SparseTuple<O>& other)const{SparseTuple<C> fout(*this); if (auto ite = other.mkIterator()) do{fout[ite()] -= *ite;} while(ite++); return fout;}


    template<class O> SparseTuple<C> operator+(const O& value)const{SparseTuple<C> fout(*this); if (auto ite = fout.mkIterator()) do{ *ite += value;} while(ite++); return fout;}
    template<class O> SparseTuple<C> operator-(const O& value)const{SparseTuple<C> fout(*this); if (auto ite = fout.mkIterator()) do{ *ite -= value;} while(ite++); return fout;}
    template<class O> SparseTuple<C> operator*(const O& value)const{SparseTuple<C> fout(*this); if (auto ite = fout.mkIterator()) do{ *ite *= value;} while(ite++); return fout;}
    template<class O> SparseTuple<C> operator/(const O& value)const{SparseTuple<C> fout(*this); if (auto ite = fout.mkIterator()) do{ *ite /= value;} while(ite++); return fout;}
};
template<class C>
class SparsityComparator{ // assumes row are sorted ('sortByRowNumber' was called )
    public:
        const myHashmap<uint32_t,C>* target;
        bool operator<(const SparsityComparator& other) const;
        bool operator>(const SparsityComparator& other) const;
        bool operator>=(const SparsityComparator& other) const;
        bool operator<=(const SparsityComparator& other) const;
        bool operator==(const SparsityComparator& other) const;
        bool operator!=(const SparsityComparator& other) const;
    };
// make obsolete! std::function for the win
template<class ARG>
class Event{
public:
    typedef ARG TEMPLATE_ARG;

	static const bool IsPOD = false;
	static const bool NeedsAddLink = false;
    virtual uint32_t operator()(ARG)=0; // output 0 to be deleted 0xFFFFFFFF to be ignored, and anything else to be re-inserted (with a increment)
    virtual ~Event(){}

    virtual void show(FILE*f =stdout, int level=0)const{fprintf(f,"An event%c", level == 0 ? '\n' : ' ');}
    };
template< >
class Event<void>{
public:
    typedef uint32_t TEMPLATE_ARG; // to prevent compiler error, should not matter
    virtual unsigned int operator()(){return 0;} // output 0 to be deleted 0xFFFFFFFF to be ignored, and anything else to be re-inserted (with a increment)
    virtual ~Event(){}
	virtual double getProgress()const{return 0.0f;}
	virtual void wrName(char* buffer256)const{strcpy(buffer256,"Unnamed Task");}
    virtual void show(FILE*f =stdout, int level=0)const{fprintf(f,"An event%c", level == 0 ? '\n' : ' ');}
    };

// is anything including pointer, but will not no free data and has default constructor/move/assign/destuctor

/*
template<class T>
class AsyncDeque{
    std::deque<T> data;
    std::mutex mut;
public:
    std::unique_lock<std::mutex> try_lock(){return std::unique_lock<std::mutex>(mut, std::try_to_lock);}
    std::deque<T>& unlocked(){return data;} // beware! must have lock!

    void push_back(const T& val);
    void push_front(const T& val);
    void pop_back();
    void pop_front();
};*/

template<class C>
class AsyncStructure{
    C& target;
    atomic<int> sema;
public:
    class ReadAccess{
        AsyncStructure& target;
        int state;
        public:
        ReadAccess(AsyncStructure& _target):target(_target){}
        ~ReadAccess(){target.sema.fetch_add(-state);}
        operator bool ();
        const C& operator()(){return target.target;}
    };
    class WriteAccess{
        AsyncStructure& target;
        int state;
        public:
        WriteAccess(AsyncStructure& _target):target(_target){}
        ~WriteAccess(){target.sema.fetch_add(-state);}
        operator bool ();
        C& operator()(){return target.target;}
    };
    AsyncStructure():sema(0){}
    ReadAccess tryRead(){return ReadAccess(*this);}
    WriteAccess tryWrite(){return WriteAccess(*this);}
};

class ProgressBarPrint{
    public:
    unsigned int state;
    unsigned int length;
    unsigned int lasttime;
    ProgressBarPrint(uint32_t _length);
    void start(const char*);
    void update(double fraction);
    void update(int nbdone, int total){return this->update(((double)nbdone) / total);}
    void finish();
};

// a task is a routine that needs to be executed on K elements, this can be done in parallel
// a control function is called before and after:
// first with the number of instantiated tasks that groups the K elements (argument is non-zero)
// and last when all instantiated tasks have completed (argument is zero)
class TaskNode{
public:
	typedef uint32_t INDEX_TYPE;
	uint32_t alias;
	string taskname;
	std::function< uint32_t(uint32_t)> ctrl_fnct;
	std::function< uint32_t(uint32_t)> fnct;
	Tuple<uint32_t> inputs;
	uint32_t& getIndex(){return alias;}
	uint32_t getIndex() const {return alias;}
};

class ThreadBase;
class Workflow{
	friend class ThreadBase;
	Workflow(ThreadBase& _tb):tb(_tb){}
public:
	ThreadBase& tb;
	Vector<Tuple<uint32_t,2u> > task_order;
	ERRCODE operator()();
	void show(FILE* F, int level) const;
};

class RangeIterator{
public:
	uint32_t length;
	uint32_t cur;
	RangeIterator(uint32_t _length): length(_length){}
	operator bool(){cur = 0; return length != 0;}
	uint32_t operator ()()const{return cur;}
	bool initLast(){cur = length -1; return cur != 0xFFFFFFFF;}
	bool operator++(int){if (++cur < length) return true; cur--; return false;}
	bool operator--(int){if (--cur != 0xFFFFFFFF) return true; cur++; return false;}
	const uint32_t& operator*(){return cur;}
	const uint32_t* operator->(){return &cur;}
};

class AtomicRangeIterator{
	atomic<uint32_t> range_iterator;
	uint32_t maxval;
public:
	AtomicRangeIterator()= default;
	AtomicRangeIterator(const AtomicRangeIterator&)= delete;
	AtomicRangeIterator(uint32_t _maxval): range_iterator(0), maxval(_maxval){}
	uint32_t operator()(bool& gotitem){uint32_t fout = range_iterator.fetch_add(1); gotitem = fout < maxval; return fout;}
	AtomicRangeIterator& setSize(uint32_t length){range_iterator = 0; maxval = length; return *this;}
};


template<class ITE>
class AtomicIterator{
	ITE ite;
public:
	typedef typename std::remove_reference<decltype(*ite)>::type ELEMENT;
	typedef typename std::remove_reference<decltype(ite())>::type KEYELEMENT;
	typedef typename std::remove_reference<decltype(*ite)>::type ELEMENT_INNERTYPE;

	char flag;
	uint32_t prog; // just dont...
	uint32_t prog_max; // just dont...


	class Iterator{
	public:
		AtomicIterator& producer;
		typename std::conditional<std::is_pointer<decltype(*(producer.ite))>::value, decltype(*(producer.ite)), typename std::remove_reference<decltype(*(producer.ite))>::type* >::type target;
		typename std::remove_reference<decltype(producer.ite())>::type key;
		Iterator(AtomicIterator& _producer): producer(_producer){}
		operator bool(){return producer.produce(target,key);}
		bool operator++(int){return producer.produce(target,key);}
		typename std::remove_reference<decltype(producer.ite())>::type& operator()(){return key;}
		template<typename std::enable_if< std::is_pointer<decltype(*(producer.ite))>::value , bool>::type = true> typename std::remove_reference<decltype(*(producer.ite))>::type* operator->(){return &target;}
		template<typename std::enable_if< !std::is_pointer<decltype(*(producer.ite))>::value , bool>::type = true> typename std::remove_reference<decltype(*(producer.ite))>::type* operator->(){return target;}
		template<typename std::enable_if< std::is_pointer<decltype(*(producer.ite))>::value , bool>::type = true> typename std::remove_reference<decltype(*(producer.ite))>::type& operator*(){return target;}
		template<typename std::enable_if< !std::is_pointer<decltype(*(producer.ite))>::value , bool>::type = true> typename std::remove_reference<decltype(*(producer.ite))>::type& operator*(){return *target;}
	};


	AtomicIterator(ITE _ite): ite(_ite), flag(0){}
	AtomicIterator(ITE&& _ite): ite(_ite), flag(0){}
	AtomicIterator(AtomicIterator&& o): ite(std::move(o.ite)), flag(o.flag),prog(o.prog), prog_max(o.prog_max){}
	AtomicIterator(const AtomicIterator& other)=delete;
	AtomicIterator& operator=(const AtomicIterator& other)=delete;

	Iterator operator()(){return Iterator(*this);}
	bool produce(ELEMENT& writeto, KEYELEMENT &key);
	bool produce(ELEMENT*& writeto, KEYELEMENT &key);
	// , typename std::enable_if< Exlisten_functor<TASK, argument_type<void(TyPe)>::type& (TASK::*)() >::ans , bool>::type =true
	template<class TASK> auto run()->decltype(TASK()());
//	template<class TASK, class A> auto run(A&)->decltype(TASK()());
//	template<class TASK, class A, class B> auto run(A&, B&)->decltype(TASK()());
	template<class TASK, class A> auto run(A&)->decltype(TASK()(declval<A&>()));
	template<class TASK, class A, class B> auto run(A&, B&)->decltype(TASK()(declval<A&>(), declval<B&>()));

};

// Structure that controls threads and holds optional argument information available to tasks
class ThreadBase{
    std::mutex mut;
    std::condition_variable condvar;
    std::condition_variable main_condvar;
    std::mutex endmut;
    std::deque<std::function<void()> > todolist;
    void operator()(); // slave threads loop
    bool running;
    int nbactivethr;

    std::mutex final_mutex; // for foreign light-weight structures
    std::condition_variable final_condvar; // for foreign light-weight structures
    uint32_t nbthreads; // for thread array
	template<class I> friend class AtomicIterator;
public:
    std::thread** thrds;
    uint32_t nbactive;
    std::vector< std::thread  > futures;
    myHashmap<uint32_t, std::thread* > dedicated;
    //std::Vector< KeyElem<std::thread::id, string> > thrname;
    std::thread::id mainthreadid;

    myHashmap<TaskNode,void> tasknodes;
    Workflow makeWorkflow(uint32_t target);


    uint32_t async_progress_maintain;
    atomic<uint32_t> async_progress;
    uint32_t async_progress_max;

    ProgressBarPrint progb;
    myHashmap<string,bool> flags; // wannabe anything!
    HeapTree<string, 3> msgs;
    FILE* log;

	// if a task has N independent things to compute in any order
	// let it have a function that pull asynchoneously jobIDs
    class ThreadArrayScope{
    public:
        ThreadBase &tb;
        uint32_t nbthreads;
        atomic<uint32_t> range_iterator;
        uint32_t range_limit;
        bool hasProgress;
        ThreadArrayScope(ThreadBase &_tb, uint32_t _nbthreads): tb(_tb), nbthreads(_nbthreads),hasProgress(false){if (_nbthreads != 0) tb.startThreadArray(_nbthreads);}
        ~ThreadArrayScope(){if (nbthreads != 0) tb.stopThreadArray();}
        ThreadArrayScope(const ThreadArrayScope&) = delete;
        ThreadArrayScope(ThreadArrayScope&& o): tb(o.tb),nbthreads(o.nbthreads), hasProgress(o.hasProgress){o.nbthreads =0;}
        ThreadArrayScope& operator=(const ThreadArrayScope&) = delete;
        ThreadArrayScope& operator=(ThreadArrayScope&& o)= delete;
        template<class FUNCTOR> void submit(FUNCTOR & func);
        void initRange(uint32_t maxval, const char* progstring =NULL){range_iterator =0; range_limit = maxval; if (progstring) {hasProgress = true; tb.startProgress(progstring, maxval);}}
        bool getRangeVal(uint32_t& val);
    };

    class RangeIterator{
        ThreadBase* tb;
        int cur;
        int limit;
        public:
        int32_t thrID;
        RangeIterator(ThreadBase* _tb,int _start, int _limit, uint32_t _thrID): tb(_tb),cur(_start), limit(_limit), thrID(_thrID){}
        operator bool (){return cur < limit;}
        bool operator++(int){if (tb) tb->updateProgress(thrID); return ((++cur) < limit);}
        uint32_t operator()()const{return cur;}
	};

	RangeIterator getRangeIterator(int thrID, int nbval, bool doprogress=false){return RangeIterator((doprogress) ? this : NULL,(thrID *  nbval) / nbthreads, ((thrID+1) *  nbval) / nbthreads, thrID);}

    ThreadBase();
    ThreadBase(const ThreadBase&)=delete;
    ThreadBase(ThreadBase&&)=delete;
    ThreadBase& operator=(const ThreadBase&)=delete;
    ThreadBase& operator=(ThreadBase&&)=delete;
    ~ThreadBase();

    ThreadBase& toMemfree();
    ThreadBase& toSize(uint32_t s){this->toMemfree(); nbthreads = s; thrds = new std::thread*[nbthreads]; return *this;}
    void setSize(uint32_t s){this->toSize(s);}

    uint32_t startThread(uint32_t thr_input, std::function<int (uint32_t)> fnc); // return thread Alias, 0 if fail
    uint32_t getSize()const{return nbthreads;}
    uint32_t getThreadArraySize() const{return nbthreads;}

    ThreadArrayScope mkThreadArray(uint32_t _nbthreads = std::thread::hardware_concurrency()){return ThreadArrayScope(*this, (nbthreads == 0) ? _nbthreads : 0);}

    void startThreadArray(uint32_t nbthreads = std::thread::hardware_concurrency());
    void stopThreadArray();

    void joinThread(uint32_t thdID); // assumes thread is about to exit
    void joinThreads(); // assumes all threads are about to exit

    std::mutex& accessFinalMutex(){return final_mutex;} // *must the guarantied to not lock anything else or block*
    std::condition_variable& accessFinalCondvar(){return final_condvar;} // *must the guarantied to not lock anything else or block*



	template<class C> void show(const C& obj);
    int print(const char* str, ...);
    int printf(const char* str, ...); // concurrent protected printf
    int printf_l(const char* str, ...); // ASSUMES final mutex is Already locked, if f == NULL uses log!
    int fprintf(FILE* f,const char* str, ...); // concurrent protected printf
    int fprintf_l(FILE* f, const char* str, ...); // ASSUMES final mutex is Already locked, if f == NULL uses log!
    int printLog(const char* str, ...); // ASSUMES final mutex is unlocked
    int printLogF(const char* str, ...); // ASSUMES final mutex is Already locked

    void terminate(const char* str, ...);


    void flushMsgs(FILE* f = stdout);
    void joinArray();
    void joinAll();

    FILE* exitlog(const char* str){running = false; FILE* f =fopen("error.log", "w+"); std::fprintf(f, "Error %s\n", str); msgs.insert(string("Created error log for: ") + string(str)) ; return f;}
    bool isRunning()const{return running;}


    //template<class F, typename... A> ERRCODE startFunctor(const F& ev, A... args);
    //template<class F, typename... A> ERRCODE startFunctor_ThenWait(const F& ev, A... args);

    template<class SCOPE> ERRCODE startEvent(Event<SCOPE>* ev, typename Event<SCOPE>::TEMPLATE_ARG scp);
    template<class SCOPE> void startEvent_ThenWait(Event<SCOPE>* ev, typename Event<SCOPE>::TEMPLATE_ARG scp); // run event, then wait of *all* other threads to end
    ERRCODE startEvent(Event<void>* ev);
    void waitForAllThreads();
    void startEvent_ThenWait(Event<void>* ev); // run event, then wait of *all* other threads to end
    static uint32_t callThatEventVoid(Event<void>* ev);
    template<class S> static  uint32_t callThatEvent(KeyElem<Event<S>* , typename Event<S>::TEMPLATE_ARG> arg);

    template<class F, typename std::enable_if< Exlisten_functor< F, void (F::*)()>::ans , bool>::type  = true> void insert(F* functor);



    void submit(std::function<void()>);
    void submit_ThenWait(std::function<void()>);
    template<class FUNCTOR, typename... ARGS> void submit(FUNCTOR & func, ARGS...);
    template<class FUNCTOR, typename... ARGS> void submit_ThenWait(FUNCTOR & func, ARGS...); // uses current thread, then joins
    template<class FUNCTOR> void submit_Array(FUNCTOR & func, uint32_t);
    template<class OBJECT,class MEMBERFUNCTIONADDRESS, typename... ARGS> void submitFunction(OBJECT & func, MEMBERFUNCTIONADDRESS mf, ARGS...);
    template<class OBJECT,class MEMBERFUNCTIONADDRESS, typename... ARGS> void submitFunction_ThenWait(OBJECT & func, MEMBERFUNCTIONADDRESS mf, ARGS...);
    template<class OBJECT,class MEMBERFUNCTIONADDRESS> void submitFunction_Array(OBJECT & func, MEMBERFUNCTIONADDRESS mf, uint32_t);
    std::function<void()> getFunc();

	// well, the vadiadic way is not really working so... (requires c++x20)


    template<class TASK, class A> int execute(A&);
    template<class ITE> AtomicIterator<ITE> operator[](ITE &ite){return AtomicIterator<ITE>(ite);}
    template<class ITE> AtomicIterator<ITE> operator[](ITE &&ite){return AtomicIterator<ITE>(ite);}
 //   template<class TASK, class A, class B> int execute(A&,B);
 //   template<class TASK, class A, class B, class C> int execute(A&,B,C);
 //   template<class TASK, class A, class B, class C, class D> int execute(A&,B,C,D);

    ERRCODE runTask(Event<uint32_t>* task){uint32_t i; for(i=nbthreads-1;i>0;i--) startEvent(task, i); startEvent_ThenWait(task,i); return 0;}

    void initEqualRanges(Tuple<uint32_t> &fout, uint32_t nbelems, uint32_t nbthreads=0) const;

    void startProgress(const char* text, uint32_t totalsteps); // address of progress int, iti
    void updateProgress(uint32_t threadID);
    void finishProgress(uint32_t threadID);

//    template<class C, class I, ENABLEIF_ARRAY(C, I) = true> void tempmetacheck_array(C& daobj){}

};

extern LFHPrimitive::ThreadBase foreach;


class TaskGraph{
public:
    myHashmap<uint32_t, std::function<void()> > func;
    myHashmap<uint32_t, Tuple<uint32_t> > deps;
    void run();
};


/*
    ThreadBase tb;
    class Localtask{
        public:
        void operator()(int val){printf("%i\n", val % 10);}
        void custom(int val){printf("%c\n", (char)('A' + val));}
    };
    Localtask myy;
    tb.startThreads();

    for(int hhh=0;hhh<100;hhh++) {tb.submit(myy, hhh); tb.submitFunction(myy, Localtask::custom, hhh % 26);}
    tb.joinAll();
    printf("all done!\n");


ERRCODE exampleFunction(ThreadBase &tb, const Something& something) const{
   class Task : public Event<uint32_t>{
    public:
    ThreadBase &tb;
    const Something& something;
    Task(ThreadBase &_tb,const Something& _some):tb(_tb),something(_some){ // reference initializations for shared scope
    }
    operator bool(){return (*this)() != 0;}
    operator ERRCODE (){return (*this)();}
    ERRCODE operator()(){
    return 0;} // main task
    uint32_t operator()(uint32_t threadID){ // parallel task
    return 0;}
    };
return (ERRCODE) Task(tb,something);}
    */
template<class C>
class SparseMatrix{
    static const bool IsPOD = false;
public:
    uint32_t computeNBrows() const;
    Vector<uint32_t> compute_rowsizes(uint32_t nbrow_guess =0u) const;
    Vector<myHashmap<uint32_t,C> > data;
    //Vector<uint32_t> row_sizes; // can be 0, if it is, row sizes are not maintained...
    //uint32_t nbRows;

    /*class KeyIterator{   seems unreliable...
        uint32_t colID;
        public:
        Tuple<uint32_t, 2u> curkey;
        const SparseMatrix<C>* target;
        uint32_t rowIterator;
        KeyIterator(const SparseMatrix<C> &_tar):target(&_tar){}
        uint32_t getCol()const{return colID;}
        uint32_t getRow()const{return this->target->data[colID].deref_key(rowIterator);}
        Tuple<uint32_t, 2u> operator()()const;
        const C* operator->()const{return &(this->target->data[colID].deref(rowIterator));}
        const C& operator*()const{return this->target->data[colID].deref(rowIterator);}
        bool first(){if (this->target->data.getSize() == 0) return false; colID=0; rowIterator=0xFFFFFFFF; return this->next();}
        bool next(){if ((++rowIterator) < this->target->data[colID].getSize()) return true; rowIterator =0; while ((++colID) < this->target->data.getSize()){ if (this->target->data[colID].getSize() != 0) return true; } return false;}
        bool last(){colID = this->target->data.getSize()-1; while(colID != 0xFFFFFFFF){if ((rowIterator =this->target->data[colID].getSize()) != 0) break;} if (colID == 0xFFFFFFFF) return false; rowIterator--; return true;}
        bool prev(){return false;} // TODO?
    };*/

	class Iterator{
        Tuple<uint32_t, 2u> coor;
        uint32_t rowIterator;
        bool next(){if ((++rowIterator) < this->target.data[coor[1]].getSize()) {coor[0] = this->target.data[coor[1]].deref_key(rowIterator); return true;} rowIterator =0; while ((++coor[1]) < this->target.data.getSize()){ if (this->target.data[coor[1]].getSize() != 0) {coor[0] = this->target.data[coor[1]].deref_key(rowIterator); return true;} } return false;}
        public:
        SparseMatrix<C>& target;
        Iterator(SparseMatrix<C>& trg): target(trg){}
        operator bool (){if (this->target.data.getSize() == 0) return false; coor[1]=0; rowIterator=0xFFFFFFFF; return this->next();}
        bool operator++(int){return next();}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        C* operator->(){return &(this->target.data[coor[1]].deref(rowIterator));}
        C& operator*(){return this->target.data[coor[1]].deref(rowIterator);}
        uint32_t getRow()const{return coor[0];}
        uint32_t getCol()const{return coor[1];}
        Iterator& mkIterator(){return *this;}
	};
	class ConstIterator{
        Tuple<uint32_t, 2u> coor;
        uint32_t rowIterator;
        bool next(){if ((++rowIterator) < this->target.data[coor[1]].getSize()) {coor[0] = this->target.data[coor[1]].deref_key(rowIterator); return true;} rowIterator =0; while ((++coor[1]) < this->target.data.getSize()){ if (this->target.data[coor[1]].getSize() != 0) {coor[0] = this->target.data[coor[1]].deref_key(rowIterator); return true;} } return false;}
        public:
        const SparseMatrix<C>& target;
        ConstIterator(const SparseMatrix<C>& trg): target(trg){}
        operator bool (){if (this->target.data.getSize() == 0) return false; coor[1]=0; rowIterator=0xFFFFFFFF; return this->next();}
        bool operator++(int){return next();}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        const C* operator->(){return &(this->target.data[coor[1]].deref(rowIterator));}
        const C& operator*(){return this->target.data[coor[1]].deref(rowIterator);}
        uint32_t getRow()const{return coor[0];}
        uint32_t getCol()const{return coor[1];}
        ConstIterator& mkIterator(){return *this;}
	};



    Tuple<uint32_t, 2u> getDim() const{Tuple<uint32_t, 2u> fout;  fout[0] = computeNBrows(); fout[1] = data.getSize(); return fout;}
    uint32_t getNBcols()const{return data.getSize();}
    uint32_t getColSize(int col_offset)const{return data[col_offset].getSize();}

    explicit operator TMatrix<C>()const;

    operator Accessor<uint32_t, myHashmap<uint32_t,C> > (){return (Accessor<uint32_t, myHashmap<uint32_t,C> >)data;}

    Tuple<C> sumRows() const;
    Tuple<C> sumCols() const;


    SparseMatrix<C>::Iterator getIterator(){return SparseMatrix<C>::Iterator(*this);}
    SparseMatrix<C>::ConstIterator getIterator()const{return SparseMatrix<C>::ConstIterator(*this);}

    SparseMatrix<C>::Iterator mkIterator(){return SparseMatrix<C>::Iterator(*this);}
    SparseMatrix<C>::ConstIterator mkIterator()const{return SparseMatrix<C>::ConstIterator(*this);}

    void setNBcols(uint32_t _col){data.setSize(_col);}
    SparseMatrix<C>& toDims(Tuple<uint32_t,2u> coor){this->setNBcols(coor[1]); return *this;}

    C& addEntry(uint32_t row, uint32_t col);
    bool hasEntry(uint32_t row, uint32_t col)const{return data[col].find(row) != 0xFFFFFFFF;}
    C& addEntry(Tuple<uint32_t, 2u> coor);
    bool hasEntry(Tuple<uint32_t, 2u> coor)const{return data[coor[1]].find(coor[0]) != 0xFFFFFFFF;}
    SparseMatrix<C>& toSparseFromNonEqual(const DataGrid<C, 2u>& data, C nan_value);

    DataGrid<C,2u> mkDataGrid() const;
    template<class F> SparseMatrix<C> subsetRows(F&, ITERABLE_DECL(F) ) const;


    SparseMatrix<C>& purgeValue(C);

    const SparseTuple<C>& getColumn(uint32_t colID) const{return (const SparseTuple<C>&)  data[colID];}; // needs to be constant...

	C operator()(const Tuple<unsigned int, 2> &coor) const{return (*this)(coor[0],coor[1]);}
	C& operator()(const Tuple<unsigned int, 2> &coor){return (*this)(coor[0],coor[1]);}
    C operator()(unsigned int x, unsigned int y) const;
    C& operator()(unsigned int x, unsigned int y);

    SparseMatrix<C>& operator=(const SparseMatrix<C>& other);
    template<class D> SparseMatrix<C>& operator=(const SparseMatrix<D>& other);

    SparseMatrix<C>& toMemmove(SparseMatrix<C>& other){/*nbRows = other.nbRows; other.nbRows =0u; row_sizes.toMemmove(other.row_sizes);*/data.toMemmove(other.data);return *this;}
    SparseMatrix<C>& toMemfree(){data.toMemfree(); /*nbRows =0u; row_sizes.toMemfree();*/ return *this; }

    template<class D> SparseMatrix<C>& operator+=(SparseMatrix<D> const & other);
    template<class D> SparseMatrix<C>& operator-=(SparseMatrix<D> const & other);

    template<class D> auto operator*(SparseTuple<D> const & other) const -> SparseTuple<decltype( data.darray[0].heap[0].k.d * other.heap.k.d[0u] )>;

    void memmoveColumn(uint32_t col, SparseTuple<C>& newcol);

  //  myHashmap<uint32_t, uint32_t> makeSparseClassPartitions()const;

    ERRCODE readDenseTable(const char* path_input, Vector<string> &rownames, Vector<string> &colnames, string &tname, OptionStruct = nullptr, bool is_floatingpoint=false, char separator= '\t');
    ERRCODE readMTXTable(const char* path_prefix , Vector<string> &rownames, Vector<string> &colnames, string &tname, bool is_floatingpoint, bool is_min_addr_one = false);
    ERRCODE writeMTXTable(const char* path_prefix , const Vector<string> &rownames, const Vector<string> &colnames, const string &tname, OptionStruct optstr=nullptr, bool is_floatingpoint =false, bool is_min_addr_one = false)const;


    Tuple<uint32_t, 2u> getSizes()const{Tuple<uint32_t, 2u> fout; fout[0] = 0u; fout[1] = data.getSize();return fout;}
    SparseMatrix<C>& toNormalizedRows(); // mu =0, sigma = 1 !
    SparseMatrix<C>& permuteRows(const Tuple<uint32_t> &permutation);

    void wrDotProductKernel(Trianglix<C>& fout)const; // matrix of partial dot products
    ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f);
    const SparseMatrix<C>& show(FILE* f= stdout, int level=0)const;
    PartialGaussElem<C> mkSingleColCovar(unsigned int col, uint32_t* offsets, int nbrows)const;
    SparseMatrix& toFilter(double fraq_in_row, double fraq_in_col);
    void wrRowDistance(Trianglix<double>& fout)const;
    #ifdef Rcpp_hpp
    void rdRcpp(const Rcpp::NumericMatrix &object);
    void wrRcpp(Rcpp::NumericMatrix &object, uint32_t nbrowsuggest=0u)const;
    void rdRcppdgCMatrix(const Rcpp::S4 &object, bool transpose=false);
    void wrRcppdgCMatrix(Rcpp::S4 &object, uint32_t nbrowsuggest=0u)const;
    void wrRcppRow(uint32_t row, SEXP &object)const;
    void wrRcppCol(uint32_t col, SEXP &object, uint32_t nbrowsuggest=0u)const;
    #endif

    void sortByRowNumber();
    void sortByRowNumber(uint32_t colID);
    bool makeMaximumAsFirst(uint32_t colID); // return (did change)
    SparseMatrix<C>& toTranspose();
    SparseMatrix<C> mkTranspose() const;
    SparseMatrix<C>& toColAppend(const SparseMatrix<C>& other);
    SparseMatrix<C>& toColMemappend(SparseMatrix<C>& other);
    Vector<SparsityComparator<C> > getSparseOrdered()const;


    class UTestPaired_Task;
    TMatrix<double> UTestPaired(bool cz) const; // do a wilcoxon test for the *difference* of all paired rows

    Tuple<double> UTestZScore(const Tuple<uint32_t> &listX, const Tuple<uint32_t> &listY, Tuple<double> *opt_logit_auroc=NULL, Tuple<double> *opt_logit_seen_auroc=NULL, bool print_prog =true, bool return_logitpval=false, bool hypergeo_correction=true)const; // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported
    Tuple<double> UTestSamplingPvalue(const Tuple<double> &logitPval, uint32_t nb_P, uint32_t nb_N, const Tuple<uint32_t> &list_P, const Tuple<uint32_t> &list_N, uint32_t nbsubsamples= 1000)const;
    TMatrix<Zscore> UTestZScoreSplit(const Tuple< Tuple<uint32_t> > &lists, const Tuple<uint32_t> &ordering, uint32_t nbthreads =4u, uint32_t splitgroups = 4u, const Tuple<uint32_t>& partit = Tuple<uint32_t>(), TMatrix<double> *opt_logit_auroc=NULL, TMatrix<double> *opt_mean=NULL, TMatrix<double> *opt_drop_enrich=NULL, bool print_prog =true, bool hypergeo_correction=true, bool downsample=false)const; // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported
    Tuple< TMatrix<double> > UTestZScoreMulti(ThreadBase &tb, const TMatrix<Tuple<uint32_t> > &list_grid, Tuple<TMatrix<double> > *opt_logit_auroc=NULL,Tuple<TMatrix<double> > *opt_mean=NULL, Tuple<TMatrix<double> > *opt_drop_enrich=NULL, const TMatrix<double> *downsample = NULL)const; // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported
//    Tuple<double> wilcoxTestZScore(const Tuple< Tuple<uint32_t> > &lists, Tuple<double> *opt_logit_auroc=NULL, bool print_prog =true, bool return_logitpval=false)const; // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported

    TMatrix<double> computeTFIDF(const Tuple< Tuple<uint32_t> > &lists) const;


    TMatrix<double> KendallColumns()const; // compute rank test for each column for specific rows only statistics , Zscore and optionnaly logit Auroc reported

    class Task_getPrincipalComponents;
    TMatrix<double> getPrincipalComponents(uint32_t nbcomponents) const; // finds components in the variance spawned by collumns
};




template<class C, unsigned int SIZE>
class Trianglix{
    void HouseHolderMultiply(const C * const vec, double denum2, unsigned int  lenght, C * buffer, bool hint); // * (I + vt^H)
    void offdiagelimination_down(const C &fact,unsigned int col, C * buffer);
    void offdiagelimination_up(const C &fact,unsigned int col, C * buffer);
public:
	typedef std::integral_constant<bool, ExCo<C>::IsPOD::value > IsPOD;
    static const unsigned int totsize = TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans;
    typedef Trianglix< C, SIZE> SAFETYPE;
    C data[TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans];

    Trianglix(){}
	Trianglix(const Tuple<C, SIZE>& u); // uu'
	Trianglix(const Tuple<C, SIZE>& u, const Tuple<C, SIZE>& v); // uv' + vu'
    Trianglix(const Trianglix<C, SIZE> &other){for(unsigned int i=0;i<totsize;i++) data[i] = other.data[i];}
    unsigned int getSize()const{return SIZE;}

	void CholeskyStep_up(const C * const vec, unsigned int length, C* buf); // L * T *L'

	bool isValid()const{for(unsigned int i =0; i < totsize;i++) if (!ExOp::isValid(data[i])) return false; return true;}

	C maxEigenValue() const; // finds the maximum eigen value, uses the "Power iteration" method

explicit operator TMatrix<C, SIZE, SIZE> () const{ TMatrix<C, SIZE, SIZE> fout; unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) { fout.data[(k-2) * (SIZE +1)] = data[i]; j += k++;}else {fout.data[(k-2) + (k-2-j+i) * SIZE] = data[i]; fout.data[(k-2-j+i) + (k-2) * SIZE] = ExOp::mkTrju(data[i]);} return fout; }

//	Trianglix<C, SIZE>& makeNonSingular(double factor); // assumes all eigen values are positive, then makes the eigen value converge to the maximum eigne value, depending on factor (0 does nothing, 1 all eigen values matches (diagonal matrix)


	double WishartLogDensity(const Trianglix<C, SIZE>&, unsigned int nbsamples) const;


	//void makeQLambdaDecomposition(TMatrix<double,SIZE> &out_q, Tuple<double,SIZE> out_l) const; // solve (*this) = QDQ^t (lambda diagonal)




    template<class O> Trianglix<C, SIZE>& operator=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] = other.data[i]; return(*this);}
    template<class O> Trianglix<C, SIZE>& operator+=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] += other.data[i]; return(*this);}
    template<class O> Trianglix<C, SIZE>& operator-=( const Trianglix<O,SIZE>& other){ for(unsigned int i=0;i<totsize;i++) data[i] -= other.data[i]; return(*this);}
	Trianglix<C, SIZE>& operator*=( const double& other) {for(unsigned int i=0;i<totsize;i++) data[i] *= other; return(*this);}
	Trianglix<C, SIZE>& operator/=( const double& other) {for(unsigned int i=0;i<totsize;i++) data[i] /= other; return(*this);}


	template<class O> Trianglix<C, SIZE> operator/( const O& other)const{return (Trianglix<C, SIZE>(*this) /= other);}
    template<class O> Trianglix<C, SIZE> operator*( const O& other)const{return (Trianglix<C, SIZE>(*this) *= other);}
    template<class O> Trianglix<C, SIZE> operator+( const Trianglix<O, SIZE>& other)const{ return (Trianglix<C, SIZE>(*this) += other);}
    template<class O> Trianglix<C, SIZE> operator-( const Trianglix<O, SIZE>& other)const{ return (Trianglix<C, SIZE>(*this) -= other);}

    Trianglix<C, SIZE>& toZero(){ unsigned int i; for(i=0;i<totsize;i++) ExOp::toZero(data[i]); return(*this);}
    Trianglix<C, SIZE>& toOne(){ unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) {ExOp::toOne(data[i]); j += k++;}else ExOp::toZero(data[i]);}
    Trianglix<C, SIZE>& toRand(){ unsigned int i,j,k; for(i=0,j=0,k=2;i<totsize;i++) if (i == j) {ExOp::toRand(data[i]); data[i] = ExOp::mkrealproj(data[i]); j += k++;}else ExOp::toRand(data[i]);}
    Trianglix<C, SIZE>& toUndefined(){ unsigned int i,j,k; for(i=0;i<totsize;i++) ExOp::toUndefined(data[i]); return(*this);}

    Trianglix<C,SIZE>& operator=(const Tuple<C, SIZE>& other);
    Trianglix<C, SIZE> inverse() const;
	double pnorm() const;

	C trace_of_division(const Trianglix<C,SIZE> &divisor) const;
	template<class B, class A> Trianglix<C,SIZE>& toAddMult(const Trianglix<B,SIZE>& b, const A& c ) {unsigned int i; for(i=0;i<totsize;i++) data[i] += b.data[i] * c; return *this;}

	C& cell(unsigned int x, unsigned int y){return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}
	C cell(unsigned int x, unsigned int y) const {return data[ (x>= y) ? y + ((x * (x+1)) /2) : x + ((y * (y+1)) /2)];}

	Tuple<C,SIZE> getEigenValues()const;
	Tuple<C,SIZE> getDiagonal()const{Tuple<C,SIZE> fout;  unsigned int i,j; for(i=0,j=0;j<SIZE;i+=j+1) fout[j++] = data[i]; return(fout);}
    C determinant() const;
    C log_determinant() const;

  //  C determinant_singularguard() const; // cannot be zero but if all eigneval are 0, forces eigenval to the geometric mean of not-zero eignevals
  //  C log_determinant_singularguard() const;
    C Xformed_inner_product_of_inverse(const Tuple<C, SIZE>& other) const;
    C inv_Xformed_inner_product_singularguard( const Tuple<C, SIZE>& other, double guard_fraction = 1.0f, double *log_det = NULL) const;


    Trianglix<C, SIZE>& toTrJu(){for(unsigned int i =0; i < totsize;i++) ExOp::toTrju(data[i]); return *this;}
    Trianglix<C, SIZE> mkTrJu()const{Trianglix<C, SIZE> fout; for(unsigned int i =0; i < totsize;i++) fout.data[i] = ExOp::mkTrju(data[i]); return fout;}


    TMatrix<C,SIZE,SIZE> diagonalizer_of_inverse() const;

    double bhattacharryya(const Tuple<C,SIZE>&dev ,const Trianglix<C,SIZE>& other)const{return bhattacharryya_partial(dev,other) - 0.25f * (log(determinant())+ log(other.determinant())); }
	double bhattacharryya_partial(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&)const;

	template< unsigned int SIZE2> C Xformed_inner_product( const Tuple<C, SIZE2>& other) const;

	C Xformed_inner_product( const Tuple<C, 0u>& other) const;


	double weighted_bhattacharryya(const Tuple<C,SIZE>&dev ,const Trianglix<C,SIZE>& other, double this_weight, double other_weight)const{return weighted_battacharya_partial(dev,other,this_weight,other_weight) - 0.5f * (this_weight*log(determinant())+ other_weight*log(other.determinant())); }
	double weighted_bhattacharryya_partial(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&, double this_weight, double other_weight)const;
	double weighted_bhattacharryya_partial_MK2(const Tuple<C,SIZE>& ,const Trianglix<C,SIZE>&, double this_weight, double other_weight)const;


    template<class O> auto mkInvDivi(const Tuple<O,SIZE>& dev) -> Tuple<decltype(data[0]*dev[0]), SIZE>; //TODO
    template<class O> auto mkPosInvDivi(const Tuple<O,SIZE>& dev) -> Tuple<decltype(data[0]*dev[0]), SIZE>;//TODO

    Trianglix<C, SIZE> mkEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc) const;
    Trianglix<C, SIZE>& toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc);




    template<class F> Trianglix<C, SIZE> mkEigenTransform(F fnc) const;
    template<class F> Trianglix<C, SIZE>& toEigenTransform(F fnc);

    Trianglix<C, SIZE>& addOuterProduct(const Tuple<C, SIZE>& other, double scalar = 1.0);

    template<class O, unsigned int OSIZE, Tuple_flag Cflag> auto operator*(const Tuple<O,OSIZE,Cflag> & a)const -> Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL>;


    const Trianglix<C, SIZE>& show(FILE* f = stdout, int level= 0)const;
    string type_tostring()const;

ERRCODE save(FILE* f)const{return (fwrite(data,sizeof(C),TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans,f) == TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans) ? 0 : 1;}
ERRCODE load(FILE* f) {return (fread(data,sizeof(C),TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans,f) == TEMPLATE_TRIANGLE_NUMBER<SIZE,2>::ans) ? 0 : 1;}


    #ifdef Rcpp_hpp
    template<class RCL> void rdMatrix(const arma::Mat<RCL> &where);
    template<class RCL> Trianglix<C, SIZE>& toMatrixRecomb(const arma::Mat<RCL> &ortho, const arma::Col<RCL> &eigen); // VtDV
    template<class RCL> void wrMatrix(arma::Mat<RCL> &where)const;
    void rdMatrix(const Rcpp::NumericMatrix &where);
    void wrMatrix(Rcpp::NumericMatrix &where)const;
    #endif

};



template<class C>
class Trianglix<C, 0u>{ // square, symetric matrix (anti-symmetric in imaginary), has special mulitiplation to ensure symmetry (multiply other to the left and right at the same time)
    void offdiagelimination_up(const C &fact,unsigned int col, C * buffer);
    void offdiagelimination_up_backroutine(const C &fact,unsigned int col, C * buffer); // assumes that the row of interest is triagonal, so some values are assumed to be zero
    Trianglix<C, 0u> mkOuterInverseProduct(const Tuple<C*> outer) const;
    void given_rotation_routine(double& c, double& s, double& c2,uint32_t offset, uint32_t maxsize);
    void toInmemInverse_routine(int start, int lenght);
    void combineBlocks(const Trianglix<C, 0u> &xformed_D_inverse); // assumes top for matrix is A^-1, and B is untouched
public:
    class Iterator{
    public:
        Tuple<uint32_t, 2u> coor;
        C* cur;
        Trianglix<C, 0u>& target;
        Iterator(Trianglix<C, 0u>& trg): target(trg){}
        operator bool (){cur = target.data; coor.toZero(); return (target.t_size != 0u);}
        bool operator++(int){cur++; if (coor[0] == coor[1]) {coor[0] = 0u; return ((++coor[1]) < target.t_size);} coor[0]++; return true;}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        C* operator->()const{return cur;}
        C& operator*()const{return *cur;}
	};
    class ConstIterator{
    public:
        Tuple<uint32_t, 2u> coor;
        C* cur;
        const Trianglix<C, 0u>& target;
        ConstIterator(const Trianglix<C, 0u>& trg): target(trg){}
        operator bool (){cur = target.data; coor.toZero(); return (target.t_size != 0u);}
        bool operator++(int){cur++; if (coor[0] == coor[1]) {coor[0] = 0u; return ((++coor[1]) < target.t_size);} coor[0]++; return true;}
        Tuple<uint32_t, 2u> operator()()const{return coor;}
        const C* operator->()const{return cur;}
        const C& operator*()const{return *cur;}
	};
	Trianglix<C, 0u>::Iterator getIterator(){return Trianglix<C, 0u>::Iterator(*this);};
	Trianglix<C, 0u>::ConstIterator getIterator()const{return Trianglix<C, 0u>::Iterator(*this);};
    // scope for calculations involcing the inverse of a matrix;

    typedef Trianglix< C, 0u> SAFETYPE;

    void QR_back(const C &factC,const C &factS,unsigned int col);

    typedef std::integral_constant<bool, ExCo<C>::IS_COMMUTATIVE::value > IS_COMMUTATIVE;
    C* data;
    unsigned int t_size;

    class HHscope{
    public:
        C* tribuf;
        C* hv;
        double* normbuf;
        HHscope(uint32_t _siz);
        HHscope(const Trianglix<C, 0u> &target);
        ~HHscope();

        double cptInvProjection(C* _out, const C* input, bool do_add_instead = false)const; // Output norm

        //void setSize(uint32_t _siz);

        void runFindInverse(const Trianglix<C, 0u> &target);
        void runFindEigen(const Trianglix<C, 0u> &target);
        void show(int _size);
    };
    Trianglix(): t_size(0){}

	C maxEigenValue() const; // finds the maximum eigen value, uses the "Power iteration" method
    unsigned int getSize()const{return t_size;}

    Trianglix(const Tuple<C, 0u>& u); // = uu'
    Trianglix(const Tuple<C, 0u>& u, const Tuple<C, 0u>& v); // = uv' + vu'
    template<unsigned int SIZE, Tuple_flag TF> Trianglix(const Tuple<C, SIZE, TF>& other);
    template<unsigned int SIZE, Tuple_flag TF> Trianglix(const Tuple<C, SIZE, TF>& u, const Tuple<C, SIZE>& v);

    Trianglix(const Trianglix<C, 0u>& other);
    template<unsigned int SIZE> Trianglix(const Trianglix<C, SIZE>& other);

    ~Trianglix();
	Trianglix<C, 0u>& operator=(const Trianglix<C, 0u>&other){setSize(other.t_size); unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i]; return(*this); }
	template<class B, unsigned int SIZE> Trianglix<C, 0u>& operator=(const Trianglix<B, SIZE>&other){setSize(other.getSize()); unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i]; return(*this); }

	template<class B, class A> Trianglix<C,0u>& toAddMult(const Trianglix<B,0u>& b, const A& c){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] += b.data[i] * c; return(*this);}

	bool isValid()const{if (t_size == 0) return(true); if (data == NULL) return(false); for(unsigned int i =0; i < totsize();i++) if (!ExOp::isValid(data[i])) return false; return true;  }

	template<unsigned int SIZE> Trianglix<C, 0u>& operator=(const Trianglix<C,SIZE>&other){setSize(SIZE); unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] = other.data[i]; return(*this); }
    Trianglix<C, 0u>& toMemmove(Trianglix<C, 0u>& other);

    unsigned int totsize()const{return (t_size & 1) ? t_size * ((t_size>>1)+1) : (t_size>>1) * (t_size+1);}
    Trianglix<C, 0u>& setSize(unsigned int s);

	void offdiagelimination_down(const C &fact,unsigned int col, C * buffer);
	template<class O> void offdiagelimination_down(const O &fact,unsigned int col, C * buffer);
	void HouseHolderMultiply(const C * const vec, double denum2, unsigned int  lenght, C * buffer, bool hint); // * (I + vt^H)
    template<class O> void HouseHolderMultiply(const O* const vec, double denum2, unsigned int  lenght, C * buffer, bool hint);
    template<class O, unsigned int SIZE, Tuple_flag Cflag> auto operator*(const Tuple<O,SIZE,Cflag> & a)const -> Tuple<decltype(this->data[0]* a[0] ),SIZE,TUPLE_FLAG_NULL>;



	template<class O> Trianglix<C, 0u> operator+( const O& other)const{ return (Trianglix<C, 0u>(*this) += other);}
    template<class O> Trianglix<C, 0u> operator-( const O& other)const{ return (Trianglix<C, 0u>(*this) -= other);}
	template<class O> Trianglix<C, 0u> operator*( const O& other)const{ return (Trianglix<C, 0u>(*this) *= other);}
    template<class O> Trianglix<C, 0u> operator/( const O& other)const{ return (Trianglix<C, 0u>(*this) /= other);}

    template<class O> Trianglix<C, 0u> operator+( const Trianglix<O, 0u>& other)const{ return (Trianglix<C, 0u>(*this) += other);}
    template<class O> Trianglix<C, 0u> operator-( const Trianglix<O, 0u>& other)const{ return (Trianglix<C, 0u>(*this) -= other);}
    template<class O> Trianglix<C, 0u>& operator+=( const Trianglix<O, 0u>& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] += other.data[i]; return(*this);}
    template<class O> Trianglix<C, 0u>& operator-=( const Trianglix<O, 0u>& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] -= other.data[i]; return(*this);}

    template<class O> auto operator*(const SparseTuple<O> & a)const -> Tuple<decltype(this->data[0]* a[0] ),0u,TUPLE_FLAG_NULL>;
    template<class O> auto mkInnerMult(const Tuple<O> & a)const -> decltype(this->data[0]* a[0] );
    template<class O> auto mkInnerMult(const SparseTuple<O> & a)const -> decltype(this->data[0]* a[0] );

    template<class O> auto mkInnerMult(const Tuple<O> & a, const Tuple<O> & b)const -> decltype(this->data[0]* a[0] );
    template<class O> auto mkInnerMult(const SparseTuple<O> & a, const SparseTuple<O> & b)const -> decltype(this->data[0]* a[0] );
    template<class O> Trianglix<C, 0u>& toInnerProduct(RemoteMemory<const O,2u> & other);
    template<class O> Trianglix<C, 0u>& toOuterProduct(RemoteMemory<const O,2u> & other);

    Trianglix<C, 0u>& operator*=( const double& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] *= other; return(*this);}
    Trianglix<C, 0u>& operator/=( const double& other){unsigned int i,ts; for(i=0,ts=totsize();i<ts;i++) data[i] /= other; return(*this);}

    Trianglix<C, 0u>& toZero(){ unsigned int i,ts; ts = totsize(); for(i=0;i<ts;i++) ExOp::toZero(data[i]); return(*this);}
    Trianglix<C, 0u>& toOne(){ unsigned int i,j,k,ts; ts = totsize(); for(i=0,j=0,k=2;i<ts;i++) if (i == j) {ExOp::toOne(data[i]); j += k++;}else ExOp::toZero(data[i]); return(*this);}
    Trianglix<C, 0u>& toRand(){ unsigned int i,j,k,ts; ts = totsize(); for(i=0,j=0,k=2;i<ts;i++) if (i == j) {ExOp::toRand(data[i]); data[i] = ExOp::mkrealproj(data[i]);j += k++;}else ExOp::toRand(data[i]);return(*this);}
    Trianglix<C, 0u>& toUndefined(){unsigned int i,ts; ts = totsize(); for(i=0;i<ts;i++) ExOp::toUndefined(data[i]); return(*this);}

    C operator()(unsigned int x, unsigned int y)const{return( (x < y) ? data[((y * (y+1)) / 2) + x]: data[((x * (x+1)) / 2) + y]);}
    C& operator()(unsigned int x, unsigned int y){return( (x < y) ? data[((y * (y+1)) / 2) + x]: data[((x * (x+1)) / 2) + y]);}
    const Trianglix<C, 0u>& show(FILE* f = stdout, int level= 0)const;
    string type_tostring()const;
    Trianglix<C, 0u>& toRotate(int firstrow, double sine, double cosine); // K <- QT*K*Q , where Q rotates row f into f+1

    // to access diagonal of trianglix:
    C& operator[](unsigned int i){ if (i>= t_size) {printf("Accessing beyond diagonal!\n"); LFH_exit(1);} return data[((i * (i+3))>> 1)];}
    const C& operator[](unsigned int i)const{if (i>= t_size) {printf("Accessing beyond diagonal!\n"); LFH_exit(1);} return data[((i * (i+3))>> 1)];}

	template<class O> Trianglix<C, 0u>& tosymmetricMult(const Trianglix<O, 0u> & other); // (other^(1/2)-T) * this * other^(1/2)
    template<class O> Trianglix<C, 0u>& tosymmetricDivi(const Trianglix<O, 0u> & other); // (other^(-1/2)-T) * this * other^(-1/2)
	template<class O> Trianglix<C, 0u> mksymmetricMult(const Trianglix<O, 0u> & other) const{Trianglix<C, 0u> fout = *this; return fout.tosymmetricMult(other);}
    template<class O> Trianglix<C, 0u> mksymmetricDivi(const Trianglix<O, 0u> & other) const{Trianglix<C, 0u> fout = *this; return fout.tosymmetricDivi(other);}

    //Trianglix<C, 0u> mkCholesky

    template<unsigned int SIZE, Tuple_flag TF> C Xformed_inner_product( const Tuple<C, SIZE, TF>& other) const;
    C Xformed_inner_product( const Tuple<C, 0u>& other) const; // x^tTx
    template<class O> C Xformed_inner_product( const Tuple<O, 0u>& other) const; // x^tTx
	template<class O> C Xformed_inner_product(RemoteMemory<const O, 1u> other) const;
	template<class O> C Xformed_inner_product(RemoteMemory<O, 1u> other) const;
//	template<class O> Tuple<d> operator*( const O& other)const{ return (Trianglix<C, 0u>(*this) *= other);}


    template<class O, unsigned int SIZE, Tuple_flag OF> auto Xformed_inner_product_of_squarre(const Tuple<O,SIZE,OF>& other) const -> decltype(ExOp::mkTrjuProd(this->data[0]* other[0])); // x^tT^2x
    C Xformed_inner_product_of_inverse(const Tuple<C, 0u>& other) const; // x^tT^{-1}x
    C gaussLL(const Tuple<C, 0u>& other) const;

	template< unsigned int SIZE> Tuple<C, SIZE,TUPLE_FLAG_NULL> operator*( const Tuple<C, SIZE>& other) const; // Ky
    Tuple<C, 0u,TUPLE_FLAG_NULL> operator*( const Tuple<C, 0u>& other) const; // Ky
    template<class O> Tuple<C, 0u,TUPLE_FLAG_NULL> operator*(RemoteMemory<const O, 1u> other) const; // Ky
    template<class O> Tuple<C, 0u,TUPLE_FLAG_NULL> operator*(RemoteMemory<O, 1u> other) const; // Ky

	template< unsigned int SIZE> Tuple<C, SIZE,TUPLE_FLAG_NULL> leftDivision( const Tuple<C, SIZE>& other) const; // K^(-1)y
    Tuple<C, 0u,TUPLE_FLAG_NULL> leftDivision( const Tuple<C, 0u>& other) const; // K^(-1)y

    Trianglix<C, 0u> mkAbs() const; // makes eigenvalues real positive
    Trianglix<C, 0u> mkAbsInverse() const; // makes eigenvalues real positive, then invert
    Trianglix<C, 0u> mkAbsoft(const double &minimum) const; // makes eigenvalues >= minimum
    Trianglix<C, 0u> mkAbsoftInverse(const double &minimum) const; // makes eigenvalues >= minimum
    Trianglix<C, 0u> mkAbhard(const double &minimum) const; // makes eigenvalues >= minimum
    Trianglix<C, 0u> mkAbhardInverse(const double &minimum) const; // makes eigenvalues >= minimum

    Trianglix<C, 0u> mkEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc) const;
    Trianglix<C, 0u>& toEigenTransform(typename ExCo<C>::FUNCTION_TYPE fnc);

    template<class F> Trianglix<C, 0u> mkEigenTransform(F fnc) const;
    template<class F> Trianglix<C, 0u>& toEigenTransform(F fnc);
    template<class F> Trianglix<C, 0u>& toEigenTransformVec(F fnc);

//    template<class D, int SX, int SY> Trianglix<C, 0u> dualMultiplication(const TMatrix<D,SX,SY> &other)const; // M^t*T*M does not work... yet
    template<class D> Trianglix<C, 0u> dualDivision(const Tuple<D*> &other, uint32_t lenght)const; // M^t*T^-1*M does not work... yet
    template<class D, unsigned int SX, unsigned int SY> auto mkInnerMult(const TMatrix<D,SX,SY> &other)const -> Trianglix<decltype(other.data[0] * this->data[0]* other.data[0]), 0u>; // M^t*T*M does not work... yet
    template<class D> auto mkInnerMult(const SparseMatrix<D> &other)const -> Trianglix<decltype(other.data[0] * this->data[0]* other.data[0]), 0u>; // M^t*T*M does not work... yet


    template<class D, unsigned int SX, unsigned int SY> Trianglix<C, 0u> mkOuterMult(const TMatrix<D,SX,SY> &other)const; // M*T*M^t does not work... yet

    template<class D, unsigned int SX, unsigned int SY> Trianglix<C, 0u> dualDivision(const TMatrix<D,SX,SY> &other)const; // M^t*T^-1*M does not work... yet


    template<class D, Tuple_flag Cflag> Tuple<C, 0u> mkBackMult(const Tuple<D, 0u, Cflag> &input)const;
    template<class D, Tuple_flag Cflag> Tuple<C, 0u> mkBackDivi(const Tuple<D, 0u, Cflag> &input)const;
    template<class D, Tuple_flag Cflag> Tuple<C, 0u> mkInvDivi(const Tuple<D, 0u, Cflag> &input)const;
    template<class D> SparseTuple<C> mkInvDivi(const SparseTuple<D> &input) const;

    Trianglix<C, 0u> Xformed_outer_product( const Tuple<C, 0u>& o) const{ return(Trianglix( (*this) * o ) ); } // Kyy^TK^T
    Trianglix<C, 0u> Xformed_outer_product_of_inverse( const Tuple<C, 0u>& o) const{ return(Trianglix( this->divisionof(o) ) );} // K^{-1}yy^T(K^{-1})^T

    //Trianglix<C, 0u> operator*(const Matrix<C>& other) const; // does not work... yet
    //const Trianglix<C, 0u>& operator*=(const Matrix<C>& other); // does not work... yet
    Trianglix<C, 0u>& operator*=(const Trianglix<C, 0u>& other);

    template<unsigned int SIZE> Trianglix<C>& operator=(const Tuple<C, SIZE>& other);
    Trianglix<C>& addOuterProduct(const Tuple<C>& other, double scalar = 1.0);

	C trace_of_product(const Trianglix<C,0u> &other) const;
	C trace_of_division(const Trianglix<C,0u> &divisor) const;
	template<class O> C trace_of_division(const Trianglix<O,0u> &divisor) const;

    Matrix<C> makeMatrix()const;
    operator Matrix<C>()const;
    operator TMatrix<C,0u,0u>()const;

	double pnorm() const;
    Tuple<C,0u> getEigenValues() const;
	Tuple<C,0u> getDiagonal()const{Tuple<C,0u> fout; fout.setSize(t_size); unsigned int i,j; for(i=0,j=0;j<t_size;i+=j+1) fout[j++] = data[i]; return(fout);}

    void wrDeterminant( C& fout)const;
	template<class O> void wrDeterminant( O& fout)const;
    void wrTrace( C& fout)const;
	template<class O> void wrTrace( O& fout)const;


	C trace()const{if (t_size == 0) LFH_exit("determinant of a size 0 matrix, check for that!\n"); C fout = data[0]; unsigned int i,j; for(i=2,j=1;j<t_size;i+=j+1) {fout += data[i];j++;} return(fout); }
	C trace_of_inverse()const;

    double bhattacharryya(const Tuple<C, 0u>&dev ,const Trianglix<C, 0u>& other)const{return battacharya_partial(dev,other) - 0.25f * (log(determinant())+ log(other.determinant())); }
    double bhattacharryya_partial(const Tuple<C, 0u>& ,const Trianglix<C, 0u>&)const;


    Trianglix<C, 0u> inverse_OLD() const;
	Trianglix<C, 0u>& toInverse() {Trianglix<C, 0u> tmp = this->mkInverse(); return ExOp::toMemmove(*this,tmp); }

	void wrInverse_block(Trianglix<C, 0u>& where, uint32_t _size =0) const;

    Trianglix<C, 0u> inverse_MK2() const;

    Trianglix<C, 0u> mkInverse_m2block() const;

    /** \brief Make Inverse of Matrix, with optional additional outputs
     *
     * \param Magscaled<C> det optionally outputs determinant of this matrix
     * \param short* is_definite optionally outputs the matrix class among (0:undefinite, 1: pos-definite, 2:neg_definite)
     * \return
     *
     */
    Trianglix<C, 0u> mkInverse(Magscaled<C>* det= NULL, short* is_definite = NULL) const; // is_definite: 0 undefinite matrix, 1 positive-definite, 2 negative-definite
// log sum a_0 * a_1 * a_2 * a_3

    Magscaled<C> determinant(short* is_definite = NULL, bool verbose = false) const;
    C log_determinant(short* is_definite = NULL, bool verbose = false) const;
    TMatrix<C> getEigenVectors(Tuple<double> &eigenvalue)const;


    // requires RCPP!!!
    template<class F> TMatrix<C, 0u> makeQDecomposition(F fnc, bool transpose = false) const; // solve (*this)^power = QLQ^t then returns Qfnc[L]
    TMatrix<C, 0u> makeQLambdaDecomposition(Tuple<C, 0u> &out_l, bool transpose = false) const; // solve (*this) = QDQ^t (lambda diagonal)
	// requires RCPP!!!


    Trianglix<C, 0u> mkInverse_older() const;
	C& cell(unsigned int x, unsigned int y){return data[ (x>= y) ? y + ((x * (x+1)) >>1) : x + ((y * (y+1)) >>1)];}
	C cell(unsigned int x, unsigned int y) const {return data[ (x>= y) ? y + ((x * (x+1)) >>1) : x + ((y * (y+1)) >>1)];}

	ERRCODE save(FILE* f) const{ERRCODE fout = (fwrite(&t_size, sizeof(unsigned int), 1,f) == 1) ? 0 : 1; return fout | (fwrite(data,sizeof(C),totsize(),f) == totsize()) ? 0 : 1;}
	ERRCODE load(FILE* f) {unsigned int i; if (fread(&i, sizeof(unsigned int), 1,f) != 1) return 1; this->setSize(i); return (fread(data,sizeof(C),totsize(),f) == totsize())? 0 : 1;}

    void wrSubTrianglix(Trianglix<C> &where, const Tuple<uint32_t> &dims)const;

    #ifdef Rcpp_hpp
    template<class RCL> void rdMatrix(const arma::Mat<RCL> &where);
    template<class RCL> Trianglix<C, 0u>& toMatrixRecomb(const arma::Mat<RCL> &ortho, const arma::Col<RCL> &eigen); // VtDV
    template<class RCL> void wrMatrix(arma::Mat<RCL> &where)const;
    template<class RCL> void wrSubMatrix(arma::Mat<RCL> &where, Tuple<uint32_t> dims)const;
    void rdMatrix(const Rcpp::NumericMatrix &where);
    void wrMatrix(Rcpp::NumericMatrix &where)const;
    #endif
};

// a Task has the following structure:
// (an object-centric contructor with mandatory references and fields):       auto <- ob.tskPermuteThatThing(target, parameters)
// (function to set optional arguments): ob.setoptionflag(true)
// (1+ main function that require thread ressource) result <- ob()  *such functions takes no arguments! part of the constructor/scope*
// optionaly, that function can be triggered by an implicit conversion to boolean/int! (check if sucessful)
// (a thread function that is called if threadID) errorcode <- ob(1)

class HiddenClasses{
    public:
    template<class C> friend class Dictionary;
    class DictionaryBase{
        template<class C> friend class Dictionary;
        uint32_t findWord(const char* word, int word_lenght) const;
        uint32_t insertWord(const char* word, int word_lenght);
        Tuple<uint32_t, 4u> findPrefix(const char* prefix, int word_lenght) const;
        uint32_t flags;
        Vector< Vector<uint32_t> > word_links;
        Vector< myHashmap<uint32_t, uint32_t> > hashes;
        DictionaryBase();
        ~DictionaryBase();
        DictionaryBase(const DictionaryBase&);
        DictionaryBase(DictionaryBase&& o) : flags(o.flags), word_links(std::move(o.word_links)), hashes(std::move(o.hashes)),entries(std::move(o.entries)){}
        DictionaryBase& operator=(const DictionaryBase&);
        DictionaryBase& operator=(DictionaryBase&& o) {flags = o.flags; entries = std::move(o.entries); word_links = std::move(o.word_links); hashes = std::move(o.hashes); return *this;}
        DictionaryBase& toMemfree();
        ERRCODE save(FILE*f) const;
        ERRCODE load(FILE*f);
        public:
        Vector<char*> entries;
        void setIgnoreCaps(bool ignore){ flags = (flags & 0xFFFFFFFE) | (ignore? 1: 0);}
        void addEntry_routine(const char* str, uint32_t strlen); // , bool allow_duplicate = true
        uint32_t findEntry(const char* query) const;

        void findIndexes(Vector<uint32_t> &fout, const char* query, uint32_t max_requested = 16u) const;
       	void findKeys(Vector<const char*> &fout, const char* query, uint32_t max_requested = 16u) const;
        Vector<uint32_t> findEntries(const vector<string> &queries) const;
    };
};

template<class C>
class Dictionary : public HiddenClasses::DictionaryBase{
	public:
    Vector<C> data;
	Dictionary()=default;
	~Dictionary(){HiddenClasses::DictionaryBase::toMemfree();}

	auto addEntry(const std::string str) -> decltype(data.push_back());
	auto addEntry(const char* str) -> decltype(data.push_back());

    auto operator[](uint32_t index)const -> decltype(data[0u]){return data[index];}
	auto operator[](uint32_t index) -> decltype(data[0u]){return data[index];}

	//void registerSubstrings(int maxlen); // maxlenght a multiple of 4 is expected, rounds up otherwise

    int getSize() const{return entries.getSize();}

    uint32_t find(const char* what) const;

    Dictionary<C>& toMemfree();

    void show(FILE*f =stdout, int level=0) const;
    ERRCODE save(FILE*f) const{ERRCODE fout = HiddenClasses::DictionaryBase::save(f); fout |= data.save(f); return fout;}
    ERRCODE load(FILE*f){ERRCODE fout = HiddenClasses::DictionaryBase::load(f); fout |= data.load(f); return fout;}


};
/*
template< >
class Dictionary<void> : HiddenClasses::DictionaryBase{
public:

	auto addEntry(const std::string str)
	void addEntry(const char* str);
    ERRCODE save(FILE*f) const{ERRCODE fout = HiddenClasses::DictionaryBase::save(f); fout |= data.save(f); return fout;}
    ERRCODE load(FILE*f){ERRCODE fout = HiddenClasses::DictionaryBase::load(f); fout |= data.load(f); return fout;}
};*/
/*
template<class C>
class DicoElem{
	public:
	Dictionary base;
	Vector<C> data;
	DicoElem();
//	~dictionary();

	void setIgnoreCaps(bool ignore){base.setIgnoreCaps(ignore);}
	void addEntry(const char* str, const C &data);
	void findSubs(Vector<C> &fout, const char* query, int max_requested = 16);
    void show(FILE*f =stdout, int level=0) const;
};*/


#include "bastructs.hpp"
#include "HashMap.hpp"

} // end of namespace LFHPrimitive



#endif


