/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */

#ifndef _defined_LFHAdvanced
#define _defined_LFHAdvanced

 #include "primitive.h"

 #include<atomic>
namespace LFHPrimitive{

#define ENGINE_MY_VISION_FADE ((double)500000000000.0f)
#define ENGINE_MY_VISION_RANGE ((double)1000000000000.0f)
#define ENGINE_MY_VISION_RESOLUTION ((double)0.0005f)

#define ENGINE_MAP_SIZE ((int)65)
#define MYASSERT(A,B) if (A) exit(MessageBox(NULL,B,"ERROR",MB_OK|MB_ICONSTOP));

	class HyperBox;


	template<class C, int nbdim> class HyperPositionCube;
	template<class C, int nbdim> class HyperPositionProjectCube;

	template<unsigned int mag, bool ISSIGNED = true, unsigned int floatpt =0> class vlongint;


    template<class PARAM, class TIMECLASS = int32_t> class ParamQueue;

	template<class ARG = void> class FIFOQueue;
	template<class KEY = uint32_t> class StochasticQueue;
	template<class KEY = uint32_t, class SCP = void, int ASYNC_MAG = 0> class PriorityQueue;
    template<class TAR, class OBJ = typename TAR::ELEM_TYPE> class AsyncInserter;

    template <class MSG> class Listener;



    template<class C> class Spline;
	class Data;
	class variint;
	template<int size> class Permutation;
	template<class key, class type>	class LinkHeapFind;
	template<class key, class type>	class LinkHeap;
	template<class key, class type, class find>	class Link;
	class RKable;
	class RKenvironnement;
	class image;
	template<int nbdim> class Parobject;
	template<int nbdim> class Sphere;

//	class PathFinderQueryFunction;
	class abstractMesh;
	class abstractMeshFace;
//	class movie;

	class CustomQueue;
	class HeapStore;
	class RemoteLink;

	template<int length,  const double &relativeFrequency> class LowPassFilterConvulutionWindow; // class meant for const initialisation

	template<class intclass, unsigned int nbdim, unsigned int nblead = 0u> class HyperPosition;
    template<class intclass, unsigned int nbdim, unsigned int nblead = 0u> class HyperCursor;
	template<class intclass, unsigned int nbdim, unsigned int nblead, class NODE, unsigned int  compression = 1u> class SpacePartition;

    class AliasHost;
	template<class Host, class Resident, int Hoffoff, int Roffoff> class dxorptr;
	template<class Host, class Resident, int Hoffoff, int Roffoff, int nbptr> class dxorptrarray;
	template<class Host, class Resident, int Hoffoff, int Roffoff> class dxorlist;
	template<class Host, class Resident, int Hoffoff, int Roffoff> class xptr;
	template<class Host, class Resident, int Hoffoff, int Roffoff> class xptrref;


	template<int value> class templateInt;
	template<class C, unsigned int SIZE> class PolyThing;
	template <class Key> class defaultcomparator;
	template <class Key,class Comparator> class RBTreeNode;
	template <class Key,class Comparator> class RBTree;
	template <class Key, class Comparator = void> class RBTofDoom; // , unsigned int ALLOCMAG
   // template <class Key, class Data> class RBTofApocalypse; // , unsigned int ALLOCMAG


	class BaseManager;
	template <class T,class E> class DataManager;
	template <class T> class Ressource;
	template<class Key> class RessourceType;
	template< class Key> class RessourcePtr;
	template< class Key, unsigned int ressource_appends =0> class ResManager; // loads a const object, which can be accessed for read-only purposes
	template< class T > class EntryPtr;
	template< class KEY> class EntryManager;

Tuple<double, 3> RGBfromHCY(double h, double c, double y, bool integer_hue = false); // 0< c,y< 1; if (integer_hue) 0<h <6, otherwise 0 < h < 2pi,
Tuple<double, 3> RGBfromHCY_gamma(double h, double c, double y,double g = 2.2, bool opt_hue = true); // 0< c,y< 1; 0<h <1
Tuple<double, 3> RGBfromHCY_gamma2(double h, double c, double y,bool opt_hue = true); // 0< c,y< 1; 0<h <1
Tuple<double, 3> HCYfromRGB(Tuple<double, 3> rgb, bool integer_hue = false); // 0<h <6, 0< c,y< 1
double gammaCorrHue(double hue, double gamma);

enum FACTOR_OPERATOR_enum{
    FACTOR_OPERATOR_NULL=0,
    FACTOR_OPERATOR_INDEX_LOOKUP=1
};

class FactorStruct{
public:
    uint32_t argA;
    uint32_t argB;
    FACTOR_OPERATOR_enum oper;
};

class FactorFrame{
public:
    myHashmap<string, FactorStruct> collumns;
    Vector< Vector<Anything> > data;

    ERRCODE load(FILE* f);
    ERRCODE save(FILE* f) const;
};


/*
template<class Key, unsigned int allocsize = 1> // allocsize needs to be a power of 2
class RBTofDoom{
public:
	unsigned int size;
	unsigned char alloc_mag;
	void* tree; // key_size*3 + keys units
	void* data; // (allocsize-1) keys unit
    unsigned int min,max;

    class Iterator{
    public:
        RBTofDoom<Key,allocsize> *target;
        void* path;
        unsigned char mag;
        const Key& operator*() const;
        const Key* operator->() const;
        void operator++();
        void operator--();
        bool isValid() const;
        Iterator& toInvalid();
    };
    typename RBTofDoom<Key,allocsize>::Iterator first() const;
    typename RBTofDoom<Key,allocsize>::Iterator last() const;

	void show(FILE* f = stdout, int level =0) const;
};*/


class OrthoRef3{
public:
    char r;
    OrthoRef3(){}
    OrthoRef3(char val):r(val){}
    OrthoRef3(ORTHOREF3_enum val):r((char)val){}


    OrthoRef3& operator/=(const OrthoRef3& other);
    OrthoRef3 operator/(const OrthoRef3& other) const;

    OrthoRef3& operator*=(const OrthoRef3& other);
    OrthoRef3 operator*(const OrthoRef3& other) const;
    template<class I> I& applyTo(I&)const; // assume I is a int class, modifies last 6bit only

	static const char* oristring;
	template<class C> operator TMatrix<C,4,4>() const;

    template<class C> operator TMatrix<C,3,3>() const;
//    template<class B, class DC, class Q> operator Reference3D<B,DC,Q> () const;


    bool isPositive(){int t = (r >> 2) ^ r; return (((t ^ (t >> 1)) & 1) == 0);}

	const char* getOri()const;
	const char* getOriInv()const;


	template<class C> Tuple<C,3u> multAbs(const Tuple<C,3u> &) const; // no reflexion, rotation only (preserve sign
	//template<class C> Tuple<C,3u> multAbs(const Tuple<C,3u> &a ) const; // no reflexion, rotation only (preserve sign


    OrthoRef3 flipX()const{return (*this) * OrthoRef3(1);}
    OrthoRef3 flipY()const{return (*this) * OrthoRef3(2);}
    OrthoRef3 flipZ()const{return (*this) * OrthoRef3(4);}

    OrthoRef3 flipXdir()const{return OrthoRef3(1) * (*this);}
    OrthoRef3 flipYdir()const{return OrthoRef3(2) * (*this);}
    OrthoRef3 flipZdir()const{return OrthoRef3(4) * (*this);}

    OrthoRef3 rotX()const{return OrthoRef3(r ^ 6);} OrthoRef3& toRotX(){r ^= 6; return (*this);}
    OrthoRef3 rotY()const{return OrthoRef3(r ^ 5);} OrthoRef3& toRotY(){r ^= 5; return (*this);}
    OrthoRef3 rotZ()const{return OrthoRef3(r ^ 3);} OrthoRef3& toRotZ(){r ^= 3; return (*this);}

    OrthoRef3 rotposX()const{return (*this) * OrthoRef3(10);}
    OrthoRef3 rotposY()const{return (*this) * OrthoRef3(44);}
    OrthoRef3 rotposZ()const{return (*this) * OrthoRef3(25);}
    OrthoRef3& toRotZP(){return (*this) *= OrthoRef3(25);}
    OrthoRef3& toRotYP(){return (*this) *= OrthoRef3(44);}
    OrthoRef3& toRotXP(){return (*this) *= OrthoRef3(10);}

    OrthoRef3 rotnegX()const{return (*this) * OrthoRef3(12);}
    OrthoRef3 rotnegY()const{return (*this) * OrthoRef3(41);}
    OrthoRef3 rotnegZ()const{return (*this) * OrthoRef3(26);}

    OrthoRef3& toZero(){r = 0;return *this;}
    OrthoRef3& toOne(){r = 0;return *this;}
    OrthoRef3& toRand(){r = rand() % 48;return *this;}

    OrthoRef3& toRotZN(){return (*this) *= OrthoRef3(26);}
    OrthoRef3& toRotYN(){return (*this) *= OrthoRef3(41);}
    OrthoRef3& toRotXN(){return (*this) *= OrthoRef3(12);}



    OrthoRef3 rotXdir()const{return OrthoRef3(6) * (*this);}
    OrthoRef3 rotYdir()const{return OrthoRef3(5) * (*this);}
    OrthoRef3 rotZdir()const{return OrthoRef3(3) * (*this);}

    OrthoRef3 rotposXdir()const{return OrthoRef3(10) * (*this);}
    OrthoRef3 rotposYdir()const{return OrthoRef3(44) * (*this);}
    OrthoRef3 rotposZdir()const{return OrthoRef3(25) * (*this);}

    OrthoRef3 rotnegXdir()const{return OrthoRef3(12) * (*this);}
    OrthoRef3 rotnegYdir()const{return OrthoRef3(41) * (*this);}
    OrthoRef3 rotnegZdir()const{return OrthoRef3(26) * (*this);}

    OrthoRef3 mkInverse() const;
    OrthoRef3& toinverse();

    const OrthoRef3& show(FILE*f = stdout, int level=0) const;
    ERRCODE save(FILE*f) const{return (1 == fwrite(&r,sizeof(char),1, f)) ? 0 : 1;}
    ERRCODE load(FILE*f) {return (1 == fread(&r,sizeof(char),1, f)) ? 0 : 1;}

};

template<class Key, class Comparator>
class RBTofDoom{
	public: class Iterator;
	private:
	void swap_with_next();
	void memswap_childrens();
	void memswap_torightchild();
	void memswap_toleftchild();

	void pop_out();

	bool reach(const Key &newf); // binarysearch check equals and stops
	void reach_par(const Key &newf); // binarysearch

	bool reach_continue(const Key &newf); // found a node, but we want to loop to find all equal nodes
	template<class A> bool reach(const A &newf); // binarysearch check equals and stops

	void removeSubTree();
	void premoveSubTree();

    void remove_maintain_routine(const Key &todelete);
	// binarysearch check equals and stops


//	void compress_toppair();

void batchInit_routine(Vector<Key> &batch, bool needsort = true); // clears current tree!!!

template<class OKEY> void range_move_routine(const OKEY& min,const OKEY& max, Vector<Key> &moveout);

void resolve_doubleblack( typename RBTofDoom<Key,Comparator>::Iterator &ite);
void remove_static_routine(typename RBTofDoom<Key,Comparator>::Iterator &todel,  pair<Vector<unsigned int>, unsigned int> &list_and_lone);

bool areLonersHealty()const;

public:

    void remove_subtree_routine( typename RBTofDoom<Key,Comparator>::Iterator &ite, pair<Vector<unsigned int>, unsigned int> &list_and_lone); //  deletes parent of subtree as well, updates pointer to brother subtree
    void delete_nodes_routine(pair<Vector<unsigned int>, unsigned int> &list_and_lone); // listed nodes are assumed to be all the unreachable ones, *MAKES ANY ITERATOR INVALID*



	//	Key* top; // chunk-o-data, 2 structures
	//	Interval<Key> range;

	Comparator comp;
	unsigned int size;
	unsigned char alloc_mag;
	unsigned int* path;

//	unsigned int toins_double;
//	unsigned int nb_holes;

	int depth;
	void makeEmptyTree(unsigned int size); // initialize a most balanced tree, that it to be filled!
	Key min; //_offset;
	Key max; //_offset;
	// allocated pointer is shifted, dont access position 0!
	pair<Key, unsigned int> *tree; // the int points to the left or right child, depending if it is red or black
	// the tail of the array contains pointers to parents!

	unsigned int firstlone; // min range of lone nodes!, all lone nodes are right childs (>= than parent)

	// call this to create a perfect balanced tree, with ((2^n)-1) nodes and 2^(n+1) space, all nodes are black, but the parents of the leaves are red
	RBTofDoom(): size(0){}
	RBTofDoom(const RBTofDoom<Key,Comparator>&);
	~RBTofDoom(){if (size > 0){delete[](tree+1);delete[](path);}}

    class Iterator{

	//	void toMaxWithinBrotherSubtree_routine(int height);
	//	void toMinWithinBrotherSubtree_routine(int height);

    public:
        const RBTofDoom<Key,Comparator>& target;
        unsigned int path[63];
        Iterator(const RBTofDoom<Key,Comparator>& );

        Iterator(const RBTofDoom<Key,Comparator>& target, const Key& key); // init iterator and findGE
        template <class O> Iterator(const RBTofDoom<Key,Comparator>& target, const O& key); // init iterator and findGE
        Iterator(const RBTofDoom<Key,Comparator>& target, const Key& key, SETCMP_enum cmp); // init iterator and find based on "cmp" {SETCMP_GE,SETCMP_LE,SETCMP_LT,SETCMP_GT}
        template <class O> Iterator(const RBTofDoom<Key,Comparator>& target, const O& key, SETCMP_enum cmp); // init iterator and  find based on "cmp"


		Iterator& operator=(const Iterator& other);
        const typename RBTofDoom<Key,Comparator>::Iterator& operator++(); // prefix
        const typename RBTofDoom<Key,Comparator>::Iterator& operator--(); // (TODO) prefix
        const Key* operator->() const;
        bool isValid() const;
        Iterator& toInvalid();
        const Key& operator*() const; // MUST BE VALID!
        const Key& Prev() const;
        const Key& Next() const;
        unsigned int get_index() const;

        unsigned int Prev_index() const;
        unsigned int Next_index() const;
        bool operator==(const Iterator& other)const;
        bool operator!=(const Iterator& other)const{return !((*this) == other);}
        Key* ordering_preserving_change() const; // allows modification, beware!
        Key& ordering_preserving_reference() const; // allows modification, beware!

        bool hsPrev() const;
        bool hsNext() const;

        bool hasLeftChild() const{return (target->tree[get_index()].second - 2) < target->firstlone-1;}
        bool hasRightChild() const{return target->tree[get_index()].second > 1;}

		bool findFirst();
		bool findLast();

		bool findGE(const Key&);
		bool findLE(const Key&);
		bool findGT(const Key&);
		bool findLT(const Key&);
		bool findEQ(const Key&);

		template<class O> bool findGE(const O&);
		template<class O> bool findLE(const O&);
		template<class O> bool findGT(const O&);
		template<class O> bool findLT(const O&);
		template<class O> bool findEQ(const O&);


    };

    class QueryIterator{
    public:
        unsigned int path[62];
        unsigned int iaindex[2];
        Key iaval[64];
    };

	//void minmax_update();

	Key& orderpreserve_deref(const typename RBTofDoom<Key,Comparator>::Iterator& ite); // warning, any modification must not perturb ordering!


    template<class F> bool first(QueryIterator&, const F&)const;
    template<class F> bool next(QueryIterator&, const F&)const;

    typename RBTofDoom<Key,Comparator>::Iterator find(const Key& tval) const;

    //typename RBTofDoom<Key>::Iterator find_first(const Key& tval) const; // find the smallest which is >= key
    //typename RBTofDoom<Key>::Iterator find_last(const Key& tval) const; // find the largest which is <= key

    //template<class A> typename RBTofDoom<Key>::Iterator find_first(const A& cmpto) const; // find the smallest which is >= key
    //template<class A> typename RBTofDoom<Key>::Iterator find_last(const A& cmpto) const; // find the largest which is <= key



    Key& deref(unsigned int ite){return tree[( ite + ((1 << alloc_mag) - size) >= firstlone) ? ite + ((1 << alloc_mag) - size) + 1  : ite+1].first;}
    const Key& deref(unsigned int ite) const {return tree[ (ite + ((1 << alloc_mag) - size) >= firstlone) ? ite + ((1 << alloc_mag) - size) +1  :  ite+1].first;}

    Key& deref_index(unsigned int index){return tree[index].first;}
    const Key& deref_index(unsigned int index) const {return tree[index].first;}

    Key& operator()(const QueryIterator& ite){return tree[ (ite.path[0] ==0) ? 1 : ite.path[ite.path[0]] ].first;}
    Key operator()(const QueryIterator& ite) const {return tree[(ite.path[0] ==0) ? 1 : ite.path[ite.path[0]] ].first;}

    unsigned int getSize()const{return size; }

	Key getMin() const;
	Key getMax() const;


	//	void realloc(const Key* oldtree, Vector<Key> &batch, bool istoremove);


	// remove, can have an extra verificator!
	//	void remove(const Key &newf, bool (*isTheOne)(const Key &newf) = NULL);


	unsigned int find_index(const Key &newf)const ; // 0 for not found
    template<class O_KEY> unsigned int find_index(const O_KEY &newf)const ; // 0 for not found


	void findPN(const Key &newf, unsigned int &prev, unsigned int &next)const; // 0 for not found

	void insert(const Key &newf);
	void remove(const Key &todelete);
	void remove(typename RBTofDoom<Key,Comparator>::Iterator &todelete); // WARNING, iterator is illegal afterwards

	void remove_update(typename RBTofDoom<Key,Comparator>::Iterator &todelete);

    template<class A> void remove(const A &todelete); // A,Key needs to be comparable

    void batchInit(Vector<Key> &batch, bool needsort = true); // clears current tree!!!
	void removeRange(const Key &min, const Key &max); // inclusive!

	void removeRange(const Key &min, typename RBTofDoom<Key,Comparator>::Iterator &max); // iterator will point to previous when done (ASSUMES iterator is valid!)
	void removeRange(typename RBTofDoom<Key,Comparator>::Iterator &min, const Key &max); // iterator will point to next when done     (ASSUMES iterator is valid!)
	template<class OKEY> void removeRange(const OKEY &min, typename RBTofDoom<Key,Comparator>::Iterator &max);
	template<class OKEY> void removeRange(typename RBTofDoom<Key,Comparator>::Iterator &min, const OKEY &max);

	void remove_mm(typename RBTofDoom<Key,Comparator>::Iterator &todelete); // todelete become prev elem
	void remove_pp(typename RBTofDoom<Key,Comparator>::Iterator &todelete); // todelete become prev elem

    template<class OKEY> void overwriteRange(const OKEY &min, const OKEY &max, const RBTofDoom<Key,Comparator>& );
    void overwriteRange(const Key &min, const Key &max, const RBTofDoom<Key,Comparator>& );


    template<class OKEY> void removeRange(const OKEY &min, const OKEY &max); // inclusive!
    template<class OKEY> void removeRange_faster(const OKEY &min, const OKEY &max, bool doit =false); // inclusive!
    template<class OKEY, class FKEY> void removeRange(const OKEY &min, const OKEY &max, const FKEY& filter_func); // inclusive!

	void moveRangeIn(const Key &min, const Key &max, Vector<Key>& ); // inclusive!



	unsigned int inorderFirst(unsigned int *) const;
	unsigned int inorderNext(unsigned int *) const;

	bool isRed(unsigned int) const;
	bool isLeaf(unsigned int) const;

	// class needs to implement SETCMP_enum operator()(const Key &min, const Key &max)const;
	template<class C> void intersection(vector<Key> &_out, const C &query) const;
	void intersectionInterval(vector<Key> &_out, const Key &min, const Key &max) const;

	template<class C> void intersection(Vector<Key> &_out, const C &query  ) const;
	void intersectionInterval(Vector<Key> &_out, const Key &min, const Key &max) const;



	//	void orderIncrease(); // reallocate tree! balance it!
	//	void orderDecrease(); // reallocate tree! balance it!
	//	const Key& getKey(unsigned int n) const; // shifted by 1 ! node 0 does not exist!

    RBTofDoom<Key,Comparator>& operator=(const RBTofDoom<Key,Comparator>&);

    RBTofDoom<Key,Comparator>& toZero(){toMemfree(); return *this;}
    RBTofDoom<Key,Comparator>& toMemmove(RBTofDoom<Key>& source){firstlone = source.firstlone; min = source.min; max = source.max; alloc_mag = source.alloc_mag; path = source.path; tree = source.tree; size = source.size; source.size = 0; return *this;}
    RBTofDoom<Key,Comparator>& toMemfree(){if (size> 0){delete[](tree+1);delete[](path);size=0;}return *this;}

	void show(FILE* f = stdout, int level =0) const;

	ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f);

    void test()const;
};



template<class Key>
class RBTofDoom<Key, void>{
	public: class Iterator;
	private:
	void swap_with_next();
	void memswap_childrens();
	void memswap_torightchild();
	void memswap_toleftchild();

	void pop_out();

	bool reach(const Key &newf); // binarysearch check equals and stops
	void reach_par(const Key &newf); // binarysearch

	bool reach_continue(const Key &newf); // found a node, but we want to loop to find all equal nodes
	template<class A> bool reach(const A &newf); // binarysearch check equals and stops

	void removeSubTree();
	void premoveSubTree();

//	void compress_toppair();

void batchInit_routine(Vector<Key> &batch, bool needsort = true); // clears current tree!!!

template<class OKEY> void range_move_routine(const OKEY& min,const OKEY& max, Vector<Key> &moveout);

void resolve_doubleblack( typename RBTofDoom<Key,void>::Iterator &ite);
void remove_static_routine(typename RBTofDoom<Key,void>::Iterator &todel,  pair<Vector<unsigned int>, unsigned int> &list_and_lone);

bool areLonersHealty()const;
void makeEmptyTree(unsigned int size); // initialize a most balanced tree, that it to be filled!

public:

    void remove_subtree_routine( typename RBTofDoom<Key,void>::Iterator &ite, pair<Vector<unsigned int>, unsigned int> &list_and_lone); //  deletes parent of subtree as well, updates pointer to brother subtree
    void delete_nodes_routine(pair<Vector<unsigned int>, unsigned int> &list_and_lone); // listed nodes are assumed to be all the unreachable ones, *MAKES ANY ITERATOR INVALID*



	//	Key* top; // chunk-o-data, 2 structures
	//	Interval<Key> range;


	unsigned int size;
	unsigned char alloc_mag;
	uint32_t* path;

//	unsigned int toins_double;
//	unsigned int nb_holes;



	uint32_t min_offset, max_offset;
	// allocated pointer is shifted, dont access position 0!
	pair<Key, uint32_t> *tree; // the int points to the left or right child, depending if it is red or black
	// the tail of the array contains pointers to parents!

	unsigned int firstlone; // min range of lone nodes!, all lone nodes are right childs (>= than parent)

	// call this to create a perfect balanced tree, with ((2^n)-1) nodes and 2^(n+1) space, all nodes are black, but the parents of the leaves are red
	RBTofDoom(): size(0){}
	~RBTofDoom(){if (size > 0){delete[](tree+1);delete[](path);}	}
	RBTofDoom(const RBTofDoom<Key,void> &other);
    RBTofDoom(RBTofDoom<Key,void>&& other): size(other.size), tree(other.tree), path(other.path), alloc_mag(other.alloc_mag),min_offset(other.min_offset),max_offset(other.max_offset){other.size =0;}
    RBTofDoom<Key,void>& operator=(RBTofDoom<Key,void>&& source){firstlone = source.firstlone; min_offset = source.min_offset; max_offset = source.max_offset; alloc_mag = source.alloc_mag; path = source.path; tree = source.tree; size = source.size; source.size = 0; return *this;}
    RBTofDoom<Key,void>& operator=(const RBTofDoom<Key>&);

    RBTofDoom<Key,void>& toMemmove(RBTofDoom<Key>& source){firstlone = source.firstlone; min_offset = source.min_offset; max_offset = source.max_offset; alloc_mag = source.alloc_mag; path = source.path; tree = source.tree; size = source.size; source.size = 0; return *this;}




    class Iterator{
	//	void toMaxWithinBrotherSubtree_routine(int height);
	//	void toMinWithinBrotherSubtree_routine(int height);
    public:
        const RBTofDoom<Key,void>& target;
        unsigned int path[63];
        int depth;
        Iterator(const RBTofDoom<Key,void>& );

        Iterator(const RBTofDoom<Key,void>& target, const Key& key); // init iterator and findGE
        template <class O> Iterator(const RBTofDoom<Key,void>& target, const O& key); // init iterator and findGE
        Iterator(const RBTofDoom<Key,void>& target, const Key& key, SETCMP_enum cmp); // init iterator and find based on "cmp" {SETCMP_GE,SETCMP_LE,SETCMP_LT,SETCMP_GT}
        template <class O> Iterator(const RBTofDoom<Key,void>& target, const O& key, SETCMP_enum cmp); // init iterator and  find based on "cmp"


		Iterator& operator=(const Iterator& other);
        const typename RBTofDoom<Key,void>::Iterator& operator++(); // prefix
        const typename RBTofDoom<Key,void>::Iterator& operator--(); // (TODO) prefix
        const Key* operator->() const;
        bool isValid() const;
        Iterator& toInvalid();
        const Key& operator*() const; // MUST BE VALID!
        const Key& Prev() const;
        const Key& Next() const;
        unsigned int get_serial_index() const;
        unsigned int get_index() const;

        unsigned int Prev_index() const;
        unsigned int Next_index() const;
        bool operator==(const Iterator& other)const;
        bool operator!=(const Iterator& other)const{return !((*this) == other);}
        Key* ordering_preserving_change() const; // allows modification, beware!
        Key& ordering_preserving_reference() const; // allows modification, beware!

        bool hsPrev() const;
        bool hsNext() const;

        bool hasLeftChild() const{return (target->tree[get_index()].second - 2) < target->firstlone-1;}
        bool hasRightChild() const{return target->tree[get_index()].second > 1;}

		bool findFirst();
		bool findLast();

		bool findGE(const Key&);
		bool findLE(const Key&);
		bool findGT(const Key&);
		bool findLT(const Key&);
		bool findEQ(const Key&);

		template<class O> bool findGE(const O&);
		template<class O> bool findLE(const O&);
		template<class O> bool findGT(const O&);
		template<class O> bool findLT(const O&);
		template<class O> bool findEQ(const O&);


		template<class C> bool findFirst(const C &query, Tuple<Key, 2u> &mima);
		template<class C> bool findNext(const C &query, Tuple<Key, 2u> &mima);
    };

    class QueryIterator{
    public:
        const RBTofDoom<Key,void>& target;
        unsigned int* path;
        uint32_t depth;
        Key tmin;
		Key tmax;
        QueryIterator(const RBTofDoom<Key,void>& _target);
        ~QueryIterator();
        const Key& operator()()const{return target.tree[path[depth]].first;}
        template <class QUERY> bool first(const QUERY& q);
        template <class QUERY> bool next(const QUERY& q);
    };

	//void minmax_update();

	Key& orderpreserve_deref(const typename RBTofDoom<Key,void>::Iterator& ite); // warning, any modification must not perturb ordering!
	Key& orderpreserve_deref_index(uint32_t index); // warning, any modification must not perturb ordering!

    template<class F> bool first(QueryIterator&, const F&)const;
    template<class F> bool next(QueryIterator&, const F&)const;

    typename RBTofDoom<Key,void>::Iterator find(const Key& tval) const;

    //typename RBTofDoom<Key,void>::Iterator find_first(const Key& tval) const; // find the smallest which is >= key
    //typename RBTofDoom<Key>::Iterator find_last(const Key& tval) const; // find the largest which is <= key

    //template<class A> typename RBTofDoom<Key>::Iterator find_first(const A& cmpto) const; // find the smallest which is >= key
    //template<class A> typename RBTofDoom<Key>::Iterator find_last(const A& cmpto) const; // find the largest which is <= key

    uint32_t findEQ(const Key& what) const;
    uint32_t findGE(const Key& what) const;
    uint32_t findGT(const Key& what) const;
    uint32_t findLE(const Key& what) const;
    uint32_t findLT(const Key& what) const;
    template<class C> uint32_t findEQ(const C& what) const;
    template<class C> uint32_t findGE(const C& what) const;
    template<class C> uint32_t findGT(const C& what) const;
    template<class C> uint32_t findLE(const C& what) const;
    template<class C> uint32_t findLT(const C& what) const;


    Key& deref(unsigned int ite){return tree[( ite + ((1 << alloc_mag) - size) >= firstlone) ? ite + ((1 << alloc_mag) - size) + 1  : ite+1].first;}
    const Key& deref(unsigned int ite) const {return tree[ (ite + ((1 << alloc_mag) - size) >= firstlone) ? ite + ((1 << alloc_mag) - size) +1  :  ite+1].first;}



    Key& deref_index(unsigned int index){return tree[index].first;}
    const Key& deref_index(unsigned int index) const {return tree[index].first;}

    Key& operator()(const QueryIterator& ite){return tree[ (ite.path[0] ==0) ? 1 : ite.path[ite.path[0]] ].first;}
    Key operator()(const QueryIterator& ite) const {return tree[(ite.path[0] ==0) ? 1 : ite.path[ite.path[0]] ].first;}

    unsigned int getSize()const{return size; }
    RBTofDoom<Key,void>& setSize(uint32_t size){this->toMemfree(); this->makeEmptyTree(size); return *this;}


	Key getMin() const{return tree[min_offset].first;}
	Key getMax() const{return tree[max_offset].first;}

	//	void realloc(const Key* oldtree, Vector<Key> &batch, bool istoremove);


	// remove, can have an extra verificator!
	//	void remove(const Key &newf, bool (*isTheOne)(const Key &newf) = NULL);


	unsigned int find_index(const Key &newf)const ; // 0 for not found
    template<class O_KEY> unsigned int find_index(const O_KEY &newf)const ; // 0 for not found


	void findPN(const Key &newf, unsigned int &prev, unsigned int &next)const; // 0 for not found

	void insert(const Key &newf);

	void remove(const Key &todelete);
    template<class A> void remove(const A &todelete); // A,Key needs to be comparable
	void remove(typename RBTofDoom<Key,void>::Iterator &todelete); // WARNING, iterator is illegal afterwards
	void remove_update(typename RBTofDoom<Key,void>::Iterator &todelete);
	void remove_mm(typename RBTofDoom<Key,void>::Iterator &todelete); // todelete become prev elem
	void remove_pp(typename RBTofDoom<Key,void>::Iterator &todelete); // todelete become prev elem

    void batchInit(Vector<Key> &batch, bool needsort = true); // clears current tree!!!
	void removeRange(const Key &min, const Key &max); // inclusive!

	void removeRange(const Key &min, typename RBTofDoom<Key,void>::Iterator &max); // iterator will point to previous when done (ASSUMES iterator is valid!)
	void removeRange(typename RBTofDoom<Key,void>::Iterator &min, const Key &max); // iterator will point to next when done     (ASSUMES iterator is valid!)
	template<class OKEY> void removeRange(const OKEY &min, typename RBTofDoom<Key,void>::Iterator &max);
	template<class OKEY> void removeRange(typename RBTofDoom<Key,void>::Iterator &min, const OKEY &max);


    template<class OKEY> void overwriteRange(const OKEY &min, const OKEY &max, const RBTofDoom<Key,void>& );
    void overwriteRange(const Key &min, const Key &max, const RBTofDoom<Key,void>& );


    template<class OKEY> void removeRange(const OKEY &min, const OKEY &max); // inclusive!
    template<class OKEY> void removeRange_faster(const OKEY &min, const OKEY &max, bool doit =false); // inclusive!
    template<class OKEY, class FKEY> void removeRange(const OKEY &min, const OKEY &max, const FKEY& filter_func); // inclusive!

	void moveRangeIn(const Key &min, const Key &max, Vector<Key>& ); // inclusive!



	unsigned int inorderFirst(unsigned int *) const;
	unsigned int inorderNext(unsigned int *) const;
    unsigned int inorderLast(unsigned int *) const;


	bool isRed(unsigned int) const;
	bool isLeaf(unsigned int) const;

	// class needs to implement SETCMP_enum operator()(const Key &min, const Key &max)const;
	template<class C> void intersection(vector<Key> &_out, const C &query) const;
	void intersectionInterval(vector<Key> &_out, const Key &min, const Key &max) const;

	template<class C> void intersection(Vector<Key> &_out, const C &query  ) const;
	void intersectionInterval(Vector<Key> &_out, const Key &min, const Key &max) const;


    template<class Q, class R> void queryParallel(ThreadBase &tb, Q& query, R& receiver, FUNCTOREQUIRES_DECL2(Q, bool, const Key& , const Key& ), ACCEPTOR_DECL(R, Key)) const; //, FUNCTOREQUIRES_DECL(G, int, Key)


	//	void orderIncrease(); // reallocate tree! balance it!
	//	void orderDecrease(); // reallocate tree! balance it!
	//	const Key& getKey(unsigned int n) const; // shifted by 1 ! node 0 does not exist!
    RBTofDoom<Key,void>& toZero(){toMemfree(); return *this;}
    RBTofDoom<Key,void>& toMemfree(){if (size> 0){delete[](tree+1);delete[](path);size=0;}return *this;}
	void show(FILE* f = stdout, int level =0) const;
	ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f);
    void test()const;
};

template<class Key, class Data>
class RBTofApocalypse : public RBTofDoom<Key, void>{
    public:
    Data* data;

};



/*
template<class Key>
class RBTofApocalypse{
public:
	static bool const IsPOD = false;
	uint32_t size;
	unsigned char alloc_mag;
	void makeEmptyTree(unsigned int size); // initialize a most balanced tree, that it to be filled!
	uint32_t min_offset, max_offset;
	pair<Key, uint32_t> *tree; // the int points to the left or right child, depending if it is red or black
	unsigned int firstlone; // min range of lone nodes!, all lone nodes are right childs (>= than parent)
	RBTofApocalypse(): size(0){}
	~RBTofApocalypse(){if (size > 0){delete[](tree+1);}	}

};*/

// FUNCTOR has two type OUTTYPE SCOPETYPE defined and the following functions:
//
// OUTTYPE FUNCTOR(SCOPETYPE, KEY) // extract cached
// OUTTYPE FUNCTOR(SCOPETYPE, KEY, OUTTYPE) // extract cached and merge with 1
// OUTTYPE FUNCTOR(SCOPETYPE, KEY, OUTTYPE, OUTTYPE) // extract cached and merge with 2
// bool FUNCTOR(KEY,KEY) // is Greather or equal, for ordering

// any inner node will have the result of the functor cached, where the function has been called on subtrees
// can hold minimums, sums, variances
template<class KEY, class FUNCTOR>
class RBDynaTree{
    public:
    Tuple< pair<KEY, unsigned int> > tree;
    typename FUNCTOR::OUTTYPE *cached;
	uint32_t size;
	uint32_t firstlone;
    uint32_t min_offset;
	uint32_t max_offset;

    bool isRed(unsigned int) const;
	bool isLeaf(unsigned int) const;
	ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f);
};


//m = pr
//v = pr/(1-p)


template<class C>
class SparseTrianglix{
    C& acxEntry(uint32_t); // assumes non-zero entries are added this way
    bool hasEntry(uint32_t)const;
    uint32_t getFirstID(uint32_t) const;
    bool isAFirstEntry(uint32_t)const;
    void addEntry(uint32_t);
    void removeEntry(uint32_t);
    void partition_popswap_routine(uint32_t what);
    void partition_acyclic_routine(uint16_t source, uint16_t sink, uint16_t partitionID, bool does_add_edge);
    bool isCyclicFID(uint16_t firstID) const;

public:
    //RBTofDoom<KeyElem<uint32_t, uint32_t > > attrib; // labels nodes and edges, allows to list edges quickly
    myHashmap<uint32_t, C> data; // data
    Tuple<Vector<uint16_t > > neigh;

    // entry exists if edge is acyclic *and* adjacent to a cyclic component,
    // value is next edgeID also adjacent to neighboring component (linked list, ending with 0xFFFFFFFF)

    // if next_id(this) == this, all neighbors belongs to distinct cyclic groups,
    // partitions, where acyclic and cyclic components are also partitioned.

    Tuple<uint32_t> next_id; // aka circular linked lists
    Tuple<uint32_t> part_id; // partID at first entry, ID of first entry for the rest
    myHashmap<uint32_t, uint32_t> attrib;
    Vector<Tuple<uint32_t, 4u> > part_data; // first/size/nbedges/firstEdgeID(0xFFFFFFFF for Lone component), if size = nbedges-1, acyclic component!

    static uint32_t getPositionInTriangularMatrix(uint32_t x, uint32_t y);

    class PartitionIterator{
        const SparseTrianglix& target;
        public:
        PartitionIterator(const SparseTrianglix& _target): target(_target){}
        Tuple<uint32_t> part;
        bool isCyclic();
        uint32_t getSize(){return part.getSize();}
        bool first();
        bool next();
    };

//    ~SparseTrianglix();
//    SparseTrianglix& operator=(const SparseTrianglix&);
    SparseTrianglix<C>& toMemmove(SparseTrianglix<C>&);
    SparseTrianglix<C>& toMemfree();
    SparseTrianglix<C>& toZero();
    SparseTrianglix<C>& toOne();
    uint32_t getSize()const {return next_id.getSize();}
    void setSize(int _size);

    //void wrPartitionId(vector<uint32_t> &fout);

    C& operator[](uint32_t); // to access diagonal of matrix:
    C operator[](uint32_t)const; // to access diagonal of matrix:
    C& operator()(const Tuple<uint32_t, 2u> &); // assumes non-zero entries are added this way
    C operator()(const Tuple<uint32_t, 2u> &)const;
    bool hasEntry(const Tuple<uint32_t, 2u> &coor)const;

    void updateDiagoInverse(C* diagoinv, uint32_t edge, C target)const;
    bool wouldEdgeBeInCyclicGroup(uint32_t pos)const;
    Tuple<uint32_t, 2u> deref_key(uint32_t ite);

    operator Trianglix<C>()const;

    void wrPermutation(Tuple<uint32_t> &fout);


    template<class D> ERRCODE evalRegul(const Trianglix<C> &target, const Trianglix<D> &network,const  Tuple< Trianglix<C> > &cross_train,const  Tuple< Trianglix<C> > &cross_test, Tuple<double, 11u> &details,  SparseTrianglix<C> *cross_solution = NULL, Trianglix<char> *chkmatch = NULL);

    class SearchRegulTask;
    ERRCODE searchRegul2019(const Trianglix<C> &target, const  Tuple< Trianglix<C> > &cross_train, const  Tuple< Trianglix<C> > &cross_test, const Tuple<uint32_t> &cross_test_size, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 13u> > > &details, uint32_t &io_nbedges, Trianglix<char> *chkmatch = NULL, uint32_t maxcyclic = 0u, uint32_t lowersearchbound =0, bool force_selected_number_edge =false);
    ERRCODE searchRegul(const Trianglix<C> target,const  Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > > &details, double lambda = 0.0f, Trianglix<char> *chkmatch = NULL, uint32_t maxcyclic = 0u);
    void searchRegul2017(const Trianglix<C> target,const  Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 9u> > > &details, double lambda = 0.0f, Trianglix<char> *chkmatch = NULL);
    ERRCODE searchRegulHalf(const Trianglix<C> target,const  Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > > &details, double lambda = 0.0f, Trianglix<char> *chkmatch = NULL, uint32_t maxcyclic = 0u);

    Tuple<uint32_t> getPartition(int partID, Tuple<uint32_t>* opt_getEdges = NULL, bool opt_projected =true)const;
    Tuple<uint32_t> getEdgePartition(uint32_t edge, Tuple<uint32_t>* opt_getEdges = NULL, bool opt_projected =true)const;

    Tuple<uint32_t> getConnectedSet(int nodeID)const;
    Tuple<uint32_t> getConnectedPartIDs(int nodeID)const;
    Tuple<uint32_t> getWouldBeCyclicPartition(uint32_t fictiveEdge, Tuple<uint32_t>* opt_getEdges = NULL, bool opt_projected =true, Tuple<uint32_t>* opt_outerEdges = NULL) const;
    Tuple<uint32_t> getAcyclicPath(uint32_t source, uint32_t target)const;

    #ifdef Rcpp_hpp
    void toRegularizedInverseOf(const Trianglix<C> &target);
    void addEdgeToRegInv(int i, int j, const Trianglix<C> &target);
    void remEdgeToRegInv(int i, int j, const Trianglix<C> &target);

    template<class RCL> void wrMatrix(arma::Mat<RCL> &fout) const;
    template<class RCL> void wrDetMatrix(arma::Mat<RCL> &fout) const;
    template<class RCL> void wrSubMatrix(arma::Mat<RCL> &fout, const Tuple<uint32_t> &partit) const;



    //void rdMatrix(const Rcpp::NumericMatrix &where);
    void wrMatrix(Rcpp::NumericMatrix &where)const;

    #endif

    Trianglix<double> mkInverse() const;
    Tuple<C> solveInPartition(Tuple<uint32_t> partit, uint32_t column)const; // solve Px = e_column
    void solveConstainedML(uint32_t partID, const Trianglix<C> &target);
    void solveAllAcyclicEdges(const Trianglix<C> &target, bool diagisone = false); // assumes is optimal for cyclic components

//    void solveRegularization();

    double solveWithFictiveEdge(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error, Tuple<double> &wouldbe_P, bool return_det_only =false, bool debug = false) const;
    double solveWithMissingEdge(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error, bool return_det_only =false, bool debug = false) const;
    //double solveWithFictiveEdge2017(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, Tuple<double> &wouldbe_P, const Trianglix<C> target, double &error_out, bool return_det_only=false, bool debug=false) const;
    double solveWithFictiveEdge2017(uint32_t fictiveEdge,                                   Tuple<double> &wouldbe_P, const Trianglix<C> &target, double &error_out, bool return_det_only=false, bool debug=false) const;

    void wrTrace(C& fout) const;
    template<class O> void wrTrace(O& fout) const;
    void wrDeterminant(C& fout) const;
    template<class O> void wrDeterminant(O& fout) const;
    void wrLogDeterminant(C& fout) const;
    template<class O> void wrLogDeterminant(O& fout) const;

    void show(FILE* f= stdout, int level=0)const;
    void showShort(FILE* f= stdout, int level=0)const;
    void showWarn(FILE* f= stdout, int level=0)const;
};
template<class TAR, class OBJ>
class AsyncInserter{
    public:
    OBJ buffer[256];
    TAR& target;
    atomic<int> semawrite;
    int semaread;
    AsyncInserter(TAR& _target): target(_target),semawrite(0), semaread(0){}
    ERRCODE operator<<(OBJ &what);
    ERRCODE operator<<=(OBJ &what);
    ERRCODE insert(OBJ &what);
};



template<class C>
class SparseTrianglixHanddle{ // controls a SparseTrianglix, and maintains properties (bi-connected partition, diago of inverse, determinant)
    public:
    SparseTrianglix<C> trianglix;
    C* diago_inverse;
    C* part_determinant;
    SparseTrianglixHanddle():diago_inverse(NULL),part_determinant(NULL){}

    const SparseTrianglix<C>& get()const {return trianglix;}

    void setSize(uint32_t _size);
    const SparseTrianglix<C>& toZero();
    const SparseTrianglix<C>& toOne();
    SparseTrianglixHanddle<C>& toMemmove(SparseTrianglixHanddle<C>& other);

    void setMaint_InvDiago(bool _does_maintain =true);
    void setMaint_Determinant(bool _does_maintain =true);

};

template<class C, unsigned int nbdim, unsigned int nblead>
class HyperPosition : public Tuple<C,nbdim>{
public:
	typedef std::true_type IsPOD;
    typedef HyperPosition<C,nbdim,nblead> SAFETYPE;

	void getOrderAndLeadHyper(unsigned int &order, unsigned int &lead) const;

    static int pairBox_MinMax(HyperPosition<C,nbdim,nblead> *fout, const Tuple<C, nbdim> &coor_min,const Tuple<C, nbdim> &coor_max);
    static int pairBox_CenSiz(HyperPosition<C,nbdim,nblead> *fout, const Tuple<C, nbdim> &center,const Tuple<C, nbdim> &size, bool size_minus1 = false);

	HyperPosition(){}
	HyperPosition(const HyperPosition<C,nbdim,nblead>& other){for(unsigned int i=0;i<nbdim;i++) this->data[i] = other.data[i];}
	HyperPosition(C* in): Tuple<C,nbdim>(in){}

    HyperPosition(const Tuple<C, nbdim> &coor_cent, C &width);

	HyperPosition(const Tuple<C, nbdim> &coor_min,const Tuple<C, nbdim> &coor_max);
	HyperPosition(const Tuple<C, nbdim> &coor_in,unsigned short magnitude, char last_dim);

    HyperPosition<C,nbdim,nblead>& operator=(const HyperPosition<C,nbdim,nblead>&other){for(unsigned int i=0;i<nbdim;i++) this->data[i] = other.data[i]; return(*this);}
    HyperPosition<C,nbdim,nblead>& operator=(HyperPosition<C,nbdim,nblead>&&other){for(unsigned int i=0;i<nbdim;i++) this->data[i] = other.data[i]; return(*this);}
    template<class D> HyperPosition<C,nbdim,nblead>& operator=(const HyperPosition<D,nbdim,nblead>&other){for(unsigned int i=0;i<nbdim;i++) this->data[i] = other.data[i];return(*this);}

    unsigned int getBestXorb(C &xorb) const;

	Tuple<C, nbdim> getMin() const;
	Tuple<C, nbdim> getMax() const;
    void wrMinMax(Tuple<C, nbdim>& f_min,Tuple<C, nbdim>& f_max) const;

	Tuple<C, nbdim> getCenter() const;
	C getMin(unsigned int dir) const;
	C getMax(unsigned int dir) const;

	HyperPosition<C,nbdim,nblead> getMinBox() const;
	HyperPosition<C,nbdim,nblead> getMaxBox() const;
    HyperPosition<C,nbdim,nblead>& toLargestBoxPossible();

    void setFromMinMax(const Tuple<C, nbdim> &coor_min,const Tuple<C, nbdim> &coor_max);
    void setFromCenterWidth(const Tuple<C, nbdim> &coor_center, C width);

    bool isZero() const {unsigned int i; for(i=0;i<nbdim;i++) if (!ExOp::isZero(this->data[i])) break; return i == nbdim;}
    HyperPosition<C,nbdim,nblead>& toZero() {unsigned int i; for(i=0;i<nbdim;i++) ExOp::toZero(this->data[i]); return (*this);}
    HyperPosition<C,nbdim,nblead>& toOne() {unsigned int i; for(i=0;i<nbdim;i++) ExOp::toOne(this->data[i]); return (*this);}

    char comparesimple(const HyperPosition<C,nbdim,nblead>& other) const;
	SETCMP_enum compare(const HyperPosition<C,nbdim,nblead>& other) const;
    template<class D> SETCMP_enum compare(const HyperPosition<D,nbdim,nblead>& other) const; // assumes sizeof(D) > sizeof(C)
	bool operator>(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	bool operator<(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
	bool operator>=(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_GE);}
	bool operator<=(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_LE);}
	bool operator==(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
	bool operator!=(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}

	bool doesIntersect(const HyperPosition<C,nbdim,nblead>& other) const{return((compare(other) & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);}


    // equality if contained only
    SETCMP_enum compare(const Tuple<C,nbdim>& other) const;
	bool operator>(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	bool operator<(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
	bool operator>=(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_GE);}
	bool operator<=(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_LE);}
	bool operator==(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
	bool operator!=(const Tuple<C,nbdim>& other) const{return((compare(other) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}
	bool isContained(const Tuple<C, nbdim>& other) const{return((compare(other) & SETCMP_DISJOINT) == SETCMP_EQUAL);}

	template<class D> bool operator>(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_GT); else return((other.compare(*this) & SETCMP_CMP_T_MASK) == SETCMP_LT); }
	template<class D> bool operator<(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_CMP_T_MASK) == SETCMP_LT); else return((other.compare(*this) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	template<class D> bool operator>=(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_GE); else return((other.compare(*this) & SETCMP_CMP_E_MASK) == SETCMP_LE);}
	template<class D> bool operator<=(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_CMP_E_MASK) == SETCMP_LE); else return((other.compare(*this) & SETCMP_CMP_E_MASK) == SETCMP_GE);}
	template<class D> bool operator==(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_DISJOINT) == SETCMP_EQUAL); else return((other.compare(*this) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
	template<class D> bool operator!=(const HyperPosition<D,nbdim,nblead>& other) const{if (sizeof(D) > sizeof(C)) return((compare(other) & SETCMP_DISJOINT) == SETCMP_DISJOINT); else return((other.compare(*this) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}

	void show(FILE * f = stdout, int level = 0) const;
//    __attribute__((noinline)) __attribute__((used)) DebugBuffer<3u>& dbg() const;
    unsigned short getMag() const;

	HyperPosition<C,nbdim,nblead> mkCommonContainer(const HyperPosition<C,nbdim,nblead>&) const;


	HyperPosition<C,nbdim,nblead> getParent() const;
	HyperPosition<C,nbdim,nblead> getLeftChild() const;
	HyperPosition<C,nbdim,nblead> getRightChild() const;
	HyperPosition<C,nbdim,nblead> getBrother() const;

    double rayIntersect(const Tuple<C,nbdim> &origin, const Tuple<C,nbdim> &target)const;  // returns the closest to origin
};



template<class C,unsigned int nbdim, unsigned int nblead>
class HyperCursor{
public:
    HyperPosition<C,nbdim,nblead> hyperpos;
typedef HyperCursor<C,nbdim,nblead> SAFETYPE;
	unsigned short mag; // 256 * magnitude + dimenstion
	// smallest mag is 0, mag of 65535 is for NaN for MAX magnitude
	// SpacePartition owned Variable:
	unsigned short par_mag; // 256 * magnitude + 4 * dimenstion + 1 isnotmin + 2 isnotmax
	HyperCursor(){} // default
	HyperCursor(const HyperCursor<C,nbdim,nblead>& in_hyper);
    HyperCursor<C,nbdim,nblead>& operator=(const HyperCursor<C,nbdim,nblead> &other);
    HyperCursor<C,nbdim,nblead>& operator=(HyperCursor<C,nbdim,nblead> &&other);

	HyperCursor(const HyperPosition<C,nbdim,nblead>& in_hyper);

	HyperCursor(const Tuple<C,nbdim>& in_coor, unsigned short magnitude);

static void incr_mag(unsigned short &a){a += ((a & 255) == nbdim -1) ? 257 + nblead - nbdim  : 1;}
static void decr_mag(unsigned short &a){a -= ((a & 255) == nblead) ? 257 + nblead - nbdim  : 1;}

    void updateMagnitude();
	operator HyperPosition<C,nbdim,nblead>() const;

	int getSplitDirection(){return mag & 255;}
	SETCMP_enum compare(const HyperCursor<C,nbdim,nblead>& other) const;

	bool operator>(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	bool operator<(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
	bool operator>=(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_CMP_T_MASK) != SETCMP_LT);}
	bool operator<=(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_CMP_T_MASK) != SETCMP_GT);}
	bool operator==(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
	bool operator!=(const HyperCursor<C,nbdim,nblead>& other) const{return((hyperpos.compare(other.hyperpos) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}


	bool operator>(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_CMP_T_MASK) == SETCMP_GT);}
	bool operator<(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_CMP_T_MASK) == SETCMP_LT);}
	bool operator>=(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_CMP_E_MASK) == SETCMP_GE);}
	bool operator<=(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_CMP_E_MASK) == SETCMP_LE);}
	bool operator==(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_DISJOINT) == SETCMP_EQUAL);}
	bool operator!=(const HyperPosition<C,nbdim,nblead>& other) const{return((hyperpos.compare(other) & SETCMP_DISJOINT) == SETCMP_DISJOINT);}

	bool isRealParentDirect() const;
	HyperPosition<C,nbdim,nblead> getRealParent() const;

	HyperPosition<C,nbdim,nblead> getParent() const;
    HyperCursor<C,nbdim,nblead> mkParent() const;
	HyperPosition<C,nbdim,nblead> getLeftChild() const;
	HyperPosition<C,nbdim,nblead> getRightChild() const;
	HyperPosition<C,nbdim,nblead> getBrother() const;
	HyperCursor<C,nbdim,nblead> mkBrother() const;
	HyperCursor<C,nbdim,nblead> mkLeftChild() const;
	HyperCursor<C,nbdim,nblead> mkRightChild() const;

    HyperPosition<C,nbdim,nblead> getCousin() const;

    HyperCursor<C,nbdim,nblead> NextTinyInBox(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max)const; // smallest but greater box within bounds;
    HyperCursor<C,nbdim,nblead> NextTinyInBox2(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max)const; // smallest but greater box within bounds;
    bool isIntersectingRect(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max)const; // smallest but greater box within bounds;

    HyperCursor<C,nbdim,nblead> getTiny_Next() const{return getNext().getMinBox();} // gets the first box, which is higher in the ordering and not within current box;
    HyperCursor<C,nbdim,nblead> getTiny_Prev() const{return getPrev().getMaxBox();} // gets the first box, which is higher in the ordering and not within current box;

	SETCMP_enum compareMagnitude(const HyperCursor<C,nbdim,nblead>& other) const;
	SETCMP_enum compareMagnitude(const HyperPosition<C,nbdim,nblead>& other) const;

    bool isMagnitudeLE(const HyperPosition<C,nbdim,nblead>& other)const{return this->compareMagnitude(other) != SETCMP_GT;}

	HyperCursor<C,nbdim,nblead>& toMin();
	HyperCursor<C,nbdim,nblead>& toMax();

	void wrReducedHyperCursor(HyperCursor<C,nbdim-1,nblead>& fout)const;
	void wrReducedHyperCursor(HyperCursor<C,nbdim-1,nblead>& fout, uint32_t dir)const;

    unsigned short largest_Container_mag_rect(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max) const;

    unsigned short commonContainer_mag(const HyperCursor<C,nbdim,nblead> &other) const;


	HyperCursor<C,nbdim,nblead>& toLargeMinBoxContainer(); // makes the largest container sharing the same MinBox

	HyperCursor<C,nbdim,nblead>& toMinBoxExclusiveContainer(const Tuple<C,nbdim> &other); // make the largest container *within* with matching MinBox that does not intersect with "other"

//	HyperCursor<C,nbdim,nblead>& toLargeMinBoxExclusiveContainer(const Tuple<C,nbdim> &other); // will increase size if needed


	bool findNext(); // update to next, return if it exists

	HyperPosition<C,nbdim,nblead> commonContainer(const HyperCursor<C,nbdim,nblead> &other) const;

    HyperCursor<C,nbdim,nblead>& operator=(const HyperPosition<C,nbdim,nblead> &other);

    template<class D> HyperCursor<C,nbdim,nblead>& operator=(const HyperCursor<D,nbdim,nblead> &other);
	const HyperCursor<C,nbdim,nblead>& operator++(); // prefix: gets the first box, which is higher in the ordering and of same dimensions
	const HyperCursor<C,nbdim,nblead>& operator--(); // prefix: gets the first box, which is lower in the ordering and of same dimensions
    HyperCursor<C,nbdim,nblead> getNext() const; // gets the first box, which is higher in the ordering and of same dimensions
    HyperCursor<C,nbdim,nblead> getPrev() const; // gets the first box, which is lower in the ordering and of same dimensions

	bool doesIntersect(const HyperPosition<C,nbdim,nblead> &other) const;
	bool isContainedBy(const HyperPosition<C,nbdim,nblead> &other) const;
	bool isContainedBy_safe(const HyperPosition<C,nbdim,nblead> &other) const;

	bool isOlderBrother() const; // return *this > this->kmBrother();



	C getMin(unsigned char dim) const;
	C getMax(unsigned char dim) const;

	HyperCursor<C,nbdim,nblead> getMinBox() const;
	HyperCursor<C,nbdim,nblead> getMaxBox() const;


	Tuple<C, nbdim> getMin() const;
	Tuple<C, nbdim> getMax() const;
    Tuple<C, nbdim> getCenter() const;
    void wrMinMax(Tuple<C, nbdim>& f_min,Tuple<C, nbdim>& f_max) const;
    void wrMinMax(C& f_min,C& f_max, unsigned int dir) const;

	void setOrder(unsigned char n_mag, unsigned char n_dir);
	const HyperCursor<C,nbdim,nblead>& setMagnitude(unsigned short n_mag);

    bool strictly_contains(const HyperCursor<C,nbdim,nblead> &other)const; // not working
    bool contains(const HyperPosition<C,nbdim,nblead> &other)const;

	HyperCursor<C,nbdim,nblead>& toLower(int dir);
	HyperCursor<C,nbdim,nblead>& toHigher(int dir);
	HyperCursor<C,nbdim,nblead> mkLower(int dir)const;
	HyperCursor<C,nbdim,nblead> mkHigher(int dir)const;
	HyperCursor<C,nbdim,nblead>& toParent();
	HyperCursor<C,nbdim,nblead>& toLeftChild();
	HyperCursor<C,nbdim,nblead>& toRightChild();
    HyperCursor<C,nbdim,nblead>& toBrother();
    HyperCursor<C,nbdim,nblead>& toLargestBoxPossible();
	HyperCursor<C,nbdim,nblead>& toCousin();

	HyperCursor<C,nbdim,nblead>& toMagnitudeLeftChild();
	HyperCursor<C,nbdim,nblead>& toMagnitudeParent();


    string type_tostring()const {return string("HyperCursor<") + ExOp::type_tostring(hyperpos[0]) + string(",...>");}
	void show(FILE *out = stdout, int level = 0) const;
	void show(char *out, int level = 0) const;
//    __attribute__((noinline)) __attribute__((used)) DebugBuffer<3u> dbg() const;

    ERRCODE save(FILE *f) const;
    ERRCODE load(FILE *f);
};


    template<class C, unsigned int nbdim>
    class HyperRectKeyIterator{

        HyperCursor<C,nbdim,0u > cur;
        public:
        Tuple<C, nbdim> min;
        Tuple<C, nbdim> max;
        HyperRectKeyIterator(){}
        HyperRectKeyIterator(const Tuple<C, nbdim> &_min,const Tuple<C, nbdim> &_max): min(_min), max(_max){}
        HyperCursor<C,nbdim,0u > operator()()const{return cur;}
        void setMinandMax(const Tuple<C, nbdim>  &_min,const Tuple<C, nbdim> &_max ){min = _min; max = _max;}
        void setCenterWidth(const Tuple<C, nbdim>  &_center,const C width ){min = _center - ((width+1) >> 1); max = _center + (width >> 1);}




        bool first();
        bool next();

        bool find_first(const HyperCursor<C,nbdim,0u > & item); // find first box in 3d rect that contains "item"
        void show(FILE*f = stdout, int level=0)const{fprintf(f,"HyperRectIterator:\n"); min.show(f); max.show(f);}
    };

template<class INTCLASS, unsigned int BITS>
class HyperPosition3D{
public:
    static const unsigned int SIZE = 1 + ((3 * BITS -1) / sizeof(INTCLASS));
    INTCLASS pos[(1 + ((3 * BITS -1) / sizeof(INTCLASS)))];
    HyperPosition3D<INTCLASS,BITS>& operator=(const HyperPosition<unsigned short , 3,0> &);
};

template<class C, unsigned int nbdim,unsigned int nblead>
class HyperMagnitudeFilter{
public:
    unsigned short min_mag;
    HyperMagnitudeFilter(){}
    HyperMagnitudeFilter(unsigned short _min_mag):min_mag(_min_mag){}
    template<class INNER> bool operator()(const KeyElem<HyperCursor<C, nbdim, nblead >, INNER> &what) const {return(what.k.mag >= min_mag);}
    template<class INNER> bool operator()(const KeyElem<HyperPosition<C, nbdim, nblead >, INNER> &what) const {return(what.k.getMag() >= min_mag);}
};

template<class C, unsigned int nbdim,unsigned int nblead, unsigned int hint = 0>
class HyperOrthoPlaneFilter{
public:
    C value;
    unsigned int dir;
    unsigned short min_mag;
    HyperOrthoPlaneFilter(){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir),min_mag(0){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir, unsigned short _min_mag):value(_inv), dir(_dir),min_mag(_min_mag){}
    HyperOrthoPlaneFilter(const HyperCursor<C, nbdim, nblead > &_inv):value(_inv.getMax(dir)){}

    bool operator()(const Tuple<C, nbdim > &what) const{return(what[dir] == value);}

    bool operator()(const HyperCursor<C, nbdim, nblead > &what) const{
        return((what.getMin(dir) <= value)&&(what.getMax(dir) >= value)&&(what.mag >= min_mag));
        }

    bool operator()(const HyperPosition<C, nbdim, nblead > &what) const{
        if ((what.getMin(dir) > value)||(what.getMax(dir) < value)) return false;
        // check mag!
        return true;
        }

     template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}


    const HyperOrthoPlaneFilter& show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoPlaneFilter: (dir:%i, minmag:%i, value:",(int)dir,(int)min_mag); ExOp::show(value,f,1); fprintf(f,")\n");return *this;}
};

template<class C, unsigned int nbdim,unsigned int nblead> // is minimum
class HyperOrthoPlaneFilter<C, nbdim, nblead,1>{
public:
    C value;
    unsigned int dir;
    unsigned short min_mag;
    HyperOrthoPlaneFilter(){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir),min_mag(0){}
    HyperOrthoPlaneFilter(const HyperCursor<C, nbdim, nblead > &_inv):value(_inv.getMin(dir)){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir, unsigned short _min_mag):value(_inv), dir(_dir),min_mag(_min_mag){}
    bool operator()(const Tuple<C, nbdim > &what) const{return(what[dir] == value);}
    bool operator()(const HyperCursor<C, nbdim, nblead > &what) const {return((what.getMin(dir) == value)&&(what.mag >= min_mag));}
    bool operator()(const HyperPosition<C, nbdim, nblead > &what) const {return(what.getMin(dir) == value);}
    template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}
    void show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoPlaneFilter for minimum: (dir:%i, minmag:%i, value:",(int)dir,(int)min_mag); ExOp::show(value,f,1); fprintf(f,")\n");}

//    template<class INNER> bool operator()(const KeyElem<HyperPosition<C, nbdim, nblead >, INNER> &what) const {return((what.k.getMin(dir) == value)&&(what.k.getMag() >= min_mag));}
};

template<class C, unsigned int nbdim,unsigned int nblead> // is maximum
class HyperOrthoPlaneFilter<C, nbdim, nblead, 2>{
public:
    C value;
    unsigned int dir;
    unsigned short min_mag;
    HyperOrthoPlaneFilter(){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir),min_mag(0){}
    HyperOrthoPlaneFilter(const HyperCursor<C, nbdim, nblead > &_inv):value(_inv.getMin(dir)){}
    HyperOrthoPlaneFilter(C _inv, unsigned int _dir, unsigned short _min_mag):value(_inv), dir(_dir),min_mag(_min_mag){}
    bool operator()(const Tuple<C, nbdim > &what) const{return(what[dir] == value);}
    bool operator()(const HyperCursor<C, nbdim, nblead > &what) const {return((what.getMax(dir) == value)&&(what.mag >= min_mag));}
    bool operator()(const HyperPosition<C, nbdim, nblead > &what) const {return(what.getMax(dir) == value);}
    template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}
//    template<class INNER> bool operator()(const KeyElem<HyperPosition<C, nbdim, nblead >, INNER> &what) const {return((what.k.getMax(dir) == value)&&(what.k.getMag() >= min_mag));}
    void show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoPlaneFilter for maximum: (dir:%i, minmag:%i, value:",(int)dir,(int)min_mag); ExOp::show(value,f,1); fprintf(f,")\n");}
};

template<class C,unsigned int nblead, unsigned int hint = 0>
class HyperOrthoHalfPlaneFilter{
public:
    C value;
    unsigned int dir;
    HyperOrthoHalfPlaneFilter(){}
    HyperOrthoHalfPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir){}
    HyperOrthoHalfPlaneFilter(const HyperCursor<C, 3u, nblead > &_inv):value(_inv.getMax(dir)){}

    bool operator()(const Tuple<C, 3u > &what) const{return(what[dir] == value);}

    bool operator()(const HyperCursor<C, 3u, nblead > &what) const{
        return((what.getMin(dir) <= value)&&(what.getMax(dir) >= value));
        }

    bool operator()(const HyperPosition<C, 3u, nblead > &what) const{
        if ((what.getMin(dir) > value)||(what.getMax(dir) < value)) return false;
        // check mag!
        return true;
        }

     template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}


    void show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoPlaneFilter: (dir:%i, value:",(int)dir); ExOp::show(value,f,1); fprintf(f,")\n");}
};

template<class C,unsigned int nblead> // is minimum
class HyperOrthoHalfPlaneFilter<C, nblead,1>{
public:
    C value;
    unsigned int dir;
    HyperOrthoHalfPlaneFilter(){}
    HyperOrthoHalfPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir){}
    HyperOrthoHalfPlaneFilter(const HyperCursor<C, 3u, nblead > &_inv):value(_inv.getMin(dir)){}
    bool operator()(const Tuple<C, 3u > &what) const{return(what[dir] == value);}
    bool operator()(const HyperCursor<C, 3u, nblead > &what) const {return((what.getMin(dir) == value));}
    bool operator()(const HyperPosition<C, 3u, nblead > &what) const {return(what.getMin(dir) == value);}
    template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}
    void show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoHalfPlaneFilter for minimum: (dir:%i, value:",(int)dir); ExOp::show(value,f,1); fprintf(f,")\n");}

//    template<class INNER> bool operator()(const KeyElem<HyperPosition<C, 3u, nblead >, INNER> &what) const {return((what.k.getMin(dir) == value)&&(what.k.getMag() >= min_mag));}
};

template<class C,unsigned int nblead> // is maximum
class HyperOrthoHalfPlaneFilter<C, nblead, 2>{
public:
    C value;
    unsigned int dir;
    HyperOrthoHalfPlaneFilter(){}
    HyperOrthoHalfPlaneFilter(C _inv, unsigned int _dir):value(_inv), dir(_dir){}
    HyperOrthoHalfPlaneFilter(const HyperCursor<C, 3u, nblead > &_inv):value(_inv.getMin(dir)){}
    bool operator()(const Tuple<C, 3u > &what) const{return(what[dir] == value);}
    bool operator()(const HyperCursor<C, 3u, nblead > &what) const {return((what.getMax(dir) == value));}
    bool operator()(const HyperPosition<C, 3u, nblead > &what) const {return(what.getMax(dir) == value);}
    template<class I, class D> inline bool operator()(const KeyElem<I, D> &what) const { return( (*this)(what.k));}
//    template<class INNER> bool operator()(const KeyElem<HyperPosition<C, 3u, nblead >, INNER> &what) const {return((what.k.getMax(dir) == value)&&(what.k.getMag() >= min_mag));}
    void show(FILE*f = stdout, int level =0)const{fprintf(f, "HyperOrthoHalfPlaneFilter for maximum: (dir:%i, value:",(int)dir); ExOp::show(value,f,1); fprintf(f,")\n");}

};

template<class STORETYPE,unsigned int  DIMS, unsigned int  LEAD,class NODES>
class HyperPositionQuery_Cube{
	public:
	Tuple<STORETYPE,DIMS> min;
	Tuple<STORETYPE,DIMS> max;
	HyperPositionQuery_Cube(){}
    HyperPositionQuery_Cube(const Tuple<STORETYPE,DIMS> i_min, const Tuple<STORETYPE,DIMS> i_max);
	HyperPositionQuery_Cube(const HyperPosition<STORETYPE,DIMS,LEAD> dabox);

    void initFromGaussElem(const GaussElem<Tuple<double, DIMS> >&);

    void setMinandMax(const Tuple<STORETYPE, DIMS>  &_min,const Tuple<STORETYPE, DIMS> &_max ){min = _min; max = _max;}
    void setCenterWidth(const Tuple<STORETYPE, DIMS>  &_center,const STORETYPE width){min = _center - ((width+1) >> 1); max = _center + (width >> 1);} // minus 1? To check
    void setCRect(const Tuple<STORETYPE, DIMS>  &_center,const Tuple<STORETYPE, DIMS>  &_size, bool size_minus1 =false);

    SETCMP_enum operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &query) const;
	SETCMP_enum operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &min, const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>  &max) const;

	void show(FILE*f, int level)const;
    };
template<class STORETYPE,unsigned int  DIMS, unsigned int  LEAD, class NODES>
class Query_HyperCube{
	public:
	Tuple<STORETYPE,DIMS> min;
	Tuple<STORETYPE,DIMS> max;
	Query_HyperCube(){}
    Query_HyperCube(const Tuple<STORETYPE,DIMS> i_min, const Tuple<STORETYPE,DIMS> i_max);
	Query_HyperCube(const HyperPosition<STORETYPE,DIMS,LEAD> dabox);
    void setMinandMax(const Tuple<STORETYPE, DIMS>  &_min,const Tuple<STORETYPE, DIMS> &_max ){min = _min; max = _max;}
    void setCenterWidth(const Tuple<STORETYPE, DIMS>  &_center,const STORETYPE width){min = _center - ((width+1) >> 1); max = _center + (width >> 1);} // minus 1? To check
    void setCRect(const Tuple<STORETYPE, DIMS>  &_center,const Tuple<STORETYPE, DIMS>  &_size, bool size_minus1 =false);
    bool operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &query)const {return (*this)(query.k);}
	bool operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &min, const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>  &max)const {return (*this)(min.k,max.k);}
    bool operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &query)const;
	bool operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &min, const HyperCursor<STORETYPE,DIMS,LEAD>  &max)const;
	void show(FILE*f, int level)const;
};



template<class STORETYPE,int DIMS, int LEAD,class NODES>
class HyperPositionQuery_Cube_Touch{
	public:
	Tuple<STORETYPE,DIMS> min;
	Tuple<STORETYPE,DIMS> max;
	unsigned char t_dim;
	STORETYPE value;
	HyperPositionQuery_Cube_Touch(const Tuple<STORETYPE,DIMS> i_min, const Tuple<STORETYPE,DIMS> i_max);
	HyperPositionQuery_Cube_Touch(const HyperPosition<STORETYPE,DIMS,LEAD> dabox);
	SETCMP_enum operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>  &min, const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>  &max) const;
	};

template<class STORETYPE,int DIMS,int LEAD,class NODES>
class HyperPosition_Contains{
	public:
	Tuple<STORETYPE,DIMS> value;
	HyperPosition_Contains(const Tuple<STORETYPE,DIMS> _value):value(_value){}
    SETCMP_enum operator()(const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES> &query) const;
	SETCMP_enum operator()(const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES>  &min, const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES>  &max) const;
	};


template<class STORETYPE,class NODES>
class HyperPositionQuery_LinkedCubes{
	public:
	Tuple<STORETYPE,3u> min;
	Tuple<STORETYPE,3u> max;

	Tuple<STORETYPE,3u> c_mima[4];

	Tuple<float, 3u> tan_gle[4];
	Tuple<float, 3u> tan_gle_cmp[4];

	HyperPositionQuery_LinkedCubes(){}

	void setQuery(const Tuple<STORETYPE,3u> *mima);

    SETCMP_enum operator()(const KeyElem<HyperCursor<STORETYPE,3u,0u>, NODES> &query) const;
	SETCMP_enum operator()(const KeyElem<HyperCursor<STORETYPE,3u,0u>, NODES> &min, const KeyElem<HyperCursor<STORETYPE,3u,0u>, NODES>  &max) const;
	void show(FILE*f, int level)const;
    };

// Tree of boxes, which maintain links to surrounding boxes
template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES>
class HierarchicalTree{
	// subroutine, assumes ite is first elem after point
    inline unsigned short par_mag_uneq(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, bool _debug = false) const;
    Vector<HyperCursor<STORETYPE,DIMS,LEAD> > makeBoxCover(const Tuple<STORETYPE, DIMS> &min, const Tuple<STORETYPE, DIMS> &max) const;
public:
	RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > hierarchy;

	typedef typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator Iterator;

    void updateParMag(const HyperCursor<STORETYPE,DIMS,LEAD> &where, unsigned short mag, unsigned short oldmag);
    void rayIntersection(Vector<NODES> &fout, const Tuple<STORETYPE, DIMS> &origin, const Tuple<STORETYPE, DIMS> &target, double maxrange =1.0f) const;
    void rayIntersection(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, DIMS> &origin, const Tuple<STORETYPE, DIMS> &target, double maxrange =1.0f) const;

    // find the magnitude of the smallest box containing subbox
    inline unsigned short par_mag_spacepartition(const HyperCursor<STORETYPE,DIMS, LEAD> &pos) const;
    inline unsigned short par_mag(const HyperCursor<STORETYPE,DIMS, LEAD> &pos) const;



    unsigned int getSize()const{return hierarchy.getSize();}

    void query(Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite) const;
    void query(const HyperCursor<STORETYPE,DIMS,LEAD> &key, Vector<NODES> &fout );

    void query(Vector<NODES> &fout, const Tuple<STORETYPE, 3>&pos, unsigned int width =0) const;
    void query(Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3>&pos,  unsigned int width=0) const;
    void queryRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3>&min, const Tuple<STORETYPE, 3>&max) const;
    void queryRect(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3>&min, const Tuple<STORETYPE, 3>&max) const;
    void queryCRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3>&center, const Tuple<STORETYPE, 3>&size, bool size_minus1 = false) const;
    void queryCRect(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3>&center, const Tuple<STORETYPE, 3>&size, bool size_minus1 = false) const;

    void queryLRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3> *mimas) const;
    void queryLRect(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3> *mimas) const;


    void insertElem(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data);
    void insertElem(const HyperPosition<STORETYPE,DIMS,LEAD> &key, const NODES &data){this->insertElem(HyperCursor<STORETYPE,DIMS,LEAD>(key), data);}

//    void insertElem(const pair<Tuple<STORETYPE,DIMS>,Tuple<STORETYPE,DIMS> > &rect, const NODES &data);




	HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>& toMemfree(){ hierarchy.toMemfree();return *this;}
	HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>& toMemmove(HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>& other){hierarchy.toMemmove(other.hierarchy);return *this;}

	void insertCover(const Tuple<STORETYPE,DIMS> & mincoor,const Tuple<STORETYPE,DIMS> & maxcoor, const NODES &data);
	void removeCover(const Tuple<STORETYPE,DIMS> & mincoor,const Tuple<STORETYPE,DIMS> & maxcoor, const NODES &data);


    void insertPartition(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data);
    void insertPartition(const HyperPosition<STORETYPE,DIMS,LEAD> &key, const NODES &data){this->insertPartition(HyperCursor<STORETYPE,DIMS,LEAD>(key), data);}

 //   void expandAll();

	unsigned short findUniqueContainerMagnitude(const HyperCursor<STORETYPE,DIMS,LEAD> &key) const;

    void removeElem(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data);

	void remove(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite);
	void removePP(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite);

    NODES* operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &where);
    const NODES* operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;

  //  void insertElemAndMerge(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data);
  //  void insertElemAndMerge(const HyperPosition<STORETYPE,DIMS,LEAD> &key, const NODES &data);

    bool checkHierarchySimple() const;

    bool checkHierarchyAt(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite) const;




   // void show(FILE*, int level)const;

    ERRCODE save(FILE *f) const{return hierarchy.save(f);}
    ERRCODE load(FILE *f){return hierarchy.load(f);}
    void show(FILE* f = stdout, int level = 0) const{hierarchy.show(f,level);}
};

/*

template<class TYPE, class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES, class FUNC>
class HierarchicalMap{
	public:
	myHashmap<typename FUNC::TYPE  , HierarchicalTree<STORETYPE, DIMS, LEAD, NODES> > trees;



};
*/
/*
template<class STORETYPE, unsigned int DIMS, class NODES>
class RectBoxes{
public:
    HierarchicalTree<STORETYPE,DIMS,0u,NODES> boxes;

//   void

};*/



// use hyperposition, automatically merges similar hyperposition, assuming hierachical overwrite
// needs to maintain the level of each box ancestor in a char! (storing its height)

template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES, unsigned int SUPER_COMPRESSION>
class SpacePartition{
    void initMagnitudes_routine(); // init mag, from hyperposition alone
    void clearParMagnitudes_routine(); // init par_mag to 0xFFFF (no parent)
    // { TO DO
    void initParMagnitudes_routine(); // init par_mag, from hyperposition alone
    // } TO DO
    unsigned short get_parent_order(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> & );
    void demote(const HyperCursor<STORETYPE,DIMS,LEAD>&);
	public:
    typedef SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION> SAFETYPE;
    HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES > partition;

	// Should be sorted and compressed at all time
	class Segment : public Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > {
		public:

		void shiftSegment(STORETYPE value, unsigned int dir);
		void shiftSegment(const Tuple<STORETYPE, DIMS> &value);
	};

    class KeyIterator : public AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION> > {
		public:
        KeyIterator(): AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION> >(){}
        KeyIterator(const SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& target) : AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION> >(target){}
		bool first();
        bool next();
        template<class FUNCTOR> bool next(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)

   //     bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&) = NULL);
   //     bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&) = NULL);

        bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        template<class FUNCTOR> bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);
        template<class FUNCTOR> bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);

        bool prev(){return false;}
        bool last(){return false;}
	};
    unsigned int getSize()const{return partition.getSize();}
    KeyIterator getKeyIterator()const{return KeyIterator(*this);}


    // { TO DO
    SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& operator=(const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& no_compression);
    // } TO DO
	NODES getAt(const HyperPosition<STORETYPE,DIMS,LEAD> posis, bool &isPure) const;
//	unsigned int getContainerOf(const HyperPosition<STORETYPE,DIMS,LEAD> posis) const;
	KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES> getAtPoint(const Tuple<STORETYPE,DIMS> &posis) const; // guarrantied to bepure

	void intersection(const HyperPosition<STORETYPE,DIMS,LEAD> &center, vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const; // always returns at least 1 element!
	void intersection(const HyperPosition<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const; // always returns at least 1 element!

    SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& toZero(){partition.toMemfree(); return *this;}

    SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& toMemfree(){partition.toMemfree(); return *this;}
    SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& toMemmove(SpacePartition<STORETYPE, DIMS, LEAD, NODES, SUPER_COMPRESSION>& other){partition.toMemmove(other.partition); return *this;}

    void remove(const KeyElem<const HyperCursor<STORETYPE,DIMS, LEAD>, NODES> what);
	void show(FILE* f = stdout, int level =0)const { partition.show(f,level);}

    void partitionIntersection(const HyperCursor<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const;

    // get the separation of *pure* block, stores oriantation of contact in par_mag
    void getContacts( Vector< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &out, HyperCursor<STORETYPE,DIMS,LEAD> &where)const;

    ERRCODE save(FILE *f) const { return partition.save(f);}
    ERRCODE load(FILE *f) { return partition.load(f); }

    void save_mk2(FILE *f) const;
    void load_mk2(FILE *f, unsigned int s =0);

    void save_portion(const HyperCursor<STORETYPE,DIMS,LEAD>& box, FILE *f) const;
    void load_portion(const HyperCursor<STORETYPE,DIMS,LEAD>& box, FILE *f);

	// partition functions:

    void fillvoid(); // use connectivity to fill the complete volumes
    void fillvoid_directional(int direction, NODES before); // use connectivity to fill the complete volumes, can only propagate in the input direction
	void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES); // insert and erase
    void remove_partition(const KeyElem<const HyperCursor<STORETYPE,DIMS, LEAD>, NODES> what);

    bool detect_artifacts()const; // used to detect useless inner boxes
    void remove_artifacts(); // used to detect useless inner boxes

	// void decompress(const HyperCursor<STORETYPE,DIMS, LEAD> &pos);
	unsigned int batchInsert_safe(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const Vector<KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES > >& data, unsigned int optional_first_index = 0u);
	unsigned int batchInsert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const Vector<KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES > >& data, unsigned int optional_first_index = 0u);

	Segment mkSegment(const Tuple<STORETYPE, DIMS>& min,  const Tuple<STORETYPE, DIMS>& max) const;

};

// NO compression! (super simple!)
template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES>
class SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u> {
    unsigned short get_parent_order(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> & );
    void demote(const HyperCursor<STORETYPE,DIMS,LEAD>&);
	public:

    RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > > partition;
    mutable KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > voidnode; // if no node is found, return this node

    class KeyIterator : public AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u> > {
		public:
        typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator trite;

     //   KeyIterator(): AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u> >(){}
        KeyIterator(const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& target) : AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>  , const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u> >(target), trite(target.partition){}
		bool first();
        bool next();
        bool next_withvoid();


   //     bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&) = NULL);
   //     bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&) = NULL);

        bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool next_withvoid(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);


        NODES* deref(); // get the actual node, NULL if no data entry at position

        template<class FUNCTOR> bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);
        template<class FUNCTOR> bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);
        template<class FUNCTOR> bool next(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)

        bool prev(){return false;}
        bool last(){return false;}
	};


    unsigned int getSize()const{return partition.getSize();}
    KeyIterator getKeyIterator()const{return KeyIterator(*this);}

	NODES getAt(const HyperPosition<STORETYPE,DIMS,LEAD> posis, bool &isPure) const;
    bool isEmptyAt(const HyperPosition<STORETYPE,DIMS,LEAD> &posis) const;
    bool isEmptyAt(const HyperCursor<STORETYPE,DIMS,LEAD> &posis) const;


// unsigned int getContainerOf(const HyperPosition<STORETYPE,DIMS,LEAD> posis) const;
	KeyElem< HyperPosition<STORETYPE,DIMS,LEAD> , NODES> operator()(const Tuple<STORETYPE,DIMS> &posis) const;
	KeyElem< HyperPosition<STORETYPE,DIMS,LEAD> , NODES>& operator()(const Tuple<STORETYPE,DIMS> &posis); // guarrantied to bepure

	unsigned short getContainerMagOf(const Tuple<STORETYPE,DIMS> &posis) const;

	void intersection(const HyperPosition<STORETYPE,DIMS,LEAD> &center, vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const; // always returns at least 1 element!
	void intersection(const HyperPosition<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const; // always returns at least 1 element!

    SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& toZero(){partition.toMemfree(); return *this;}
    SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& toMemfree(){partition.toMemfree(); return *this;}
    SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& toMemmove(SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0>& other){partition.toMemmove(other.partition); return *this;}

    void remove(const KeyElem<const HyperCursor<STORETYPE,DIMS, LEAD>, NODES> what);
	void show(FILE*f = stdout, int level =0)const { partition.show(f,level);}

    void clearBoxArea(const HyperCursor<STORETYPE,DIMS, LEAD> &where);

    void partitionIntersection(const HyperCursor<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const;

    // get the separation of *pure* block, stores oriantation of contact in par_mag
    void getContacts( Vector< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &out, HyperCursor<STORETYPE,DIMS,LEAD> &where)const;

    ERRCODE save(FILE *f) const { return partition.save(f);} // can be improved!!!
    ERRCODE load(FILE *f) { return partition.load(f); }  // can be improved!!!

	// partition functions:

	void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES); // insert and erase
    void remove_partition(const KeyElem<const HyperCursor<STORETYPE,DIMS, LEAD>, NODES> what);

    bool detect_artifacts()const; // used to detect useless inner boxes
    void remove_artifacts(); // used to detect useless inner boxes
};

template<class STORETYPE>
class HeirarchicalSurface3{
    public:
        STORETYPE inter_coor;
        unsigned short dir;
        unsigned short expand(HyperCursor<STORETYPE,3u,0u>&)const;
        HyperCursor<STORETYPE,3u,0u> compress(HyperCursor<STORETYPE,3u,0u>&)const;
    };


// class NODES must have 2 functions:
// bool isMergeMutable(KeyElem< HyperCursor<,,>, NODES > &f_inout) const
// void split(Tuple<,2u>& fout, HyperCursor<STORETYPE,DIMS,LEAD>&)const;

// Rules:
// A: Neighboring nodes does not intersect (strictly larger/smaller)
// B: Neighboring nodes cannot be equal (the second is redundant if so)
// C: non-trivial nodes head intersects plane (make it trivial if so)
// D: brother nodes must disagree (they should be joined)
// E: trivial nodes cannot aggree with preceding nodes (the trivial node should be removed if so)


/*
template<class NODES, class STORETYPE = typename NODES::IntType, unsigned int DIMS = NODES::nb_dims, unsigned int LEAD = NODES::nb_lead, class MERGER = typename NODES::MergeScopeType>
class HierarchicalMergable{
	//void clearHyperPosition_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where);
	bool findContainerOf_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const HyperCursor<STORETYPE,DIMS,LEAD> where);


	void simplifyCheck(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite);

	void insertAfter_routine();

	bool writePrevNext_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what,typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, NODES &prevprop, NODES &nextprop);
	bool writePrevNextMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what,typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, NODES &prevprop, NODES &nextprop, bool debug = false);

	void changeNodeSize_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos);
	void changeNodeSizeMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos, bool debug = false);
	void HyperPositionInflate_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where);
	void HyperPositionInflateMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where, bool debug = false);
	void writeHyperPosition_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what);
	void writeHyperPositionMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, NODES &what, bool debug = false);
public:
	RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > > space;
	class KeyIterator : public AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES >  , const HierarchicalMergable<NODES, STORETYPE, DIMS, LEAD, MERGER> > {
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;
		Tuple<STORETYPE,DIMS> limit;
		NODES lastnode;
		public:
        KeyIterator(const HierarchicalMergable<NODES, STORETYPE, DIMS, LEAD, MERGER>& target) : AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES >  , const HierarchicalMergable<NODES, STORETYPE, DIMS, LEAD, MERGER> >(target), ite(target.space){}
		bool first();
        bool next();
        bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool prev(){return false;}
        bool last(){return false;}
	};

	HierarchicalMergable& toZero();
	HierarchicalMergable& toMemfree(){space.toMemfree();return *this;}
	HierarchicalMergable& toMemmove(HierarchicalMergable& other){space.toMemmove(other.space);return *this;}

	unsigned int getSize()const{return space.getSize();}

    KeyIterator getKeyIterator()const{return KeyIterator(*this);}
    ERRCODE save(FILE *f) const { return space.save(f);}
    ERRCODE load(FILE *f) { return space.load(f); }


	void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node); // insert and erase

	void insertMK2(HyperCursor<STORETYPE,DIMS, LEAD> pos, NODES n_node, bool debug = false); // insert and erase
    // nodes are simple or complex, merging is allowed if equal, or possible is one of the two node is simple


	void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node, const MERGER& merger); // insert and erase


	bool getAt(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	bool getAt_simplify(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	bool getAt_update(NODES& fout,HyperCursor<STORETYPE,DIMS,LEAD> &where) const;


	bool checkIntegrity()const;
	void insert_with_bypass(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node); // simple insert, may corrupt structure if overlap
	void show(FILE*f = stdout, int level =0)const{
	    return space.show(f,level);
	    }
};*/

template<class NODES, class STORETYPE = typename NODES::IntType, unsigned int DIMS = NODES::nb_dims, unsigned int LEAD = NODES::nb_lead>
class HierarchicalMergable{
	//void clearHyperPosition_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where);
	bool findContainerOf_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const HyperCursor<STORETYPE,DIMS,LEAD> where);


	void simplifyCheck(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite);

	void insertAfter_routine();

	bool writePrevNext_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what,typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, NODES &prevprop, NODES &nextprop);
	bool writePrevNextMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what,typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, NODES &prevprop, NODES &nextprop, bool debug = false);

	void changeNodeSize_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos);
	void changeNodeSizeMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos, bool debug = false);
	void HyperPositionInflate_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where);
	void HyperPositionInflateMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where, bool debug = false);
	void writeHyperPosition_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what);
	void writeHyperPositionMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, NODES &what, bool debug = false);


	template<class MERGER> void changeNodeSize_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos, const MERGER& merger);
	template<class MERGER> void HyperPositionInflate_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where, const MERGER& merger);
	template<class MERGER> void writeHyperPosition_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what, const MERGER& merger);


	template<class MERGER> bool writePrevNextMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what,typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, NODES &prevprop, NODES &nextprop, const MERGER& merger, bool debug = false);
	template<class MERGER> void changeNodeSizeMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, const Tuple<STORETYPE,DIMS> &minpos, const MERGER& merger, bool debug = false);
	template<class MERGER> void HyperPositionInflateMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &where, const MERGER& merger, bool debug = false);
	template<class MERGER> void writeHyperPositionMK2_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, NODES &what, const MERGER& merger,bool debug = false);


	//void writeHyperPositionTrivial_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &where, const NODES& what);



	// returns "does the container exist"? (it is a (ZERO) otherwise"
public:
	RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > > space;

    HierarchicalMergable(){}
    HierarchicalMergable(const HierarchicalMergable& other):space(other.space){}
    HierarchicalMergable(HierarchicalMergable&& other):space(std::move(other.space)) {}
    HierarchicalMergable& operator=(const HierarchicalMergable& other){space = other.space; return *this;}
    HierarchicalMergable& operator=(HierarchicalMergable&& other){space = std::move(other.space); return *this;}
    HierarchicalMergable& toMemmove(HierarchicalMergable& other){space.toMemmove(other.space);return *this;}

    class KeyIterator : public AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES >  , const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> > {
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;
		Tuple<STORETYPE,DIMS> limit;
		//HyperCursor<STORETYPE,DIMS,LEAD> cur;
		NODES lastnode;
		bool first_inrect_routine(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
		public:
//		KeyIterator(){}
        KeyIterator(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& target) : AbstractKeyIterator< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES >  , const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> >(target), ite(target.space){}
		bool first();
        bool next();
        template<class FUNCTOR> bool first(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)
        template<class FUNCTOR> bool next(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)

        // TODO
        template<class FUNCTOR> bool firstBin(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)
        template<class FUNCTOR> bool nextBin(const FUNCTOR &f); // bool (*filter)(const KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES>&)
        // TODO

        bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);
        bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box);

        template<class FUNCTOR> bool first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);
        template<class FUNCTOR> bool next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const FUNCTOR &f);

        bool first(const Tuple< Tuple<STORETYPE,DIMS>, 2u> &mima); // TODO
        bool next(const Tuple< Tuple<STORETYPE,DIMS>, 2u> &mima); // TODO

        bool prev(){return false;}
        bool last(){return false;}
	};

	class Iterator{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;

        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;

		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;
		HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		bool next();
		public:
        Iterator(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const HyperCursor<STORETYPE,DIMS,LEAD>& _in_box) : ite(_target.space), target(_target), in_box(_in_box){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};

    class IteratorMima{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;


        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;

		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;
		Tuple< Tuple<STORETYPE,DIMS>, 2u> mima;
		HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		bool next();
		public:
        IteratorMima(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const Tuple< Tuple<STORETYPE,DIMS>, 2u> &_mima) : ite(_target.space), target(_target), mima(_mima){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};
    template<class M> class IteratorMimerger{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;

        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;
        const M& mergescope;
        Tuple< Tuple<STORETYPE,DIMS>, 2u> mima;
        HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;

		bool next();
		public:
        IteratorMimerger(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const Tuple< Tuple<STORETYPE,DIMS>, 2u> &_mima, const M& _mergescope) : ite(_target.space), target(_target), mima(_mima), mergescope(_mergescope){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};
	template<class M> class IteratorMerger{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;

        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;
        const M& mergescope;

		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;
		HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		bool next();
		public:
        IteratorMerger(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const HyperCursor<STORETYPE,DIMS,LEAD>& _in_box, const M& _mergescope) : ite(_target.space), target(_target), in_box(_in_box), mergescope(_mergescope){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};
    template<class F> class IteratorFunc{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;

        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;
        const F& func;

		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;
		HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		bool next();
		public:
        IteratorFunc(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const HyperCursor<STORETYPE,DIMS,LEAD>& _in_box, const F& _func) : ite(_target.space), target(_target), in_box(_in_box),func(_func){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};
	template<class F,class M> class IteratorFuncMerger{
		typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite;
		unsigned short mag;

        HyperCursor<STORETYPE,DIMS,LEAD> curkey;
        const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD> &target;
        const M& mergescope;
        const F& func;

		Tuple<STORETYPE,DIMS> limit;
		NODES simplenode;
		NODES lastnode;
		HyperCursor<STORETYPE,DIMS,LEAD> in_box;
		bool next();
		public:
        IteratorFuncMerger(const HierarchicalMergable<NODES,STORETYPE,DIMS,LEAD>& _target, const HyperCursor<STORETYPE,DIMS,LEAD>& _in_box, const F& _func, const M& _mergescope) : ite(_target.space), target(_target), in_box(_in_box),func(_func), mergescope(_mergescope){}
        operator bool ();
        bool operator++(int);
        const HyperCursor<STORETYPE,DIMS,LEAD>& operator()()const{return curkey;}
        const NODES* operator->()const{return &simplenode;}
        const NODES& operator*()const{return simplenode;}
	};
	Iterator getIterator(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box) const {return Iterator(*this, in_box);}
	IteratorMima getIterator(const Tuple< Tuple<STORETYPE,DIMS>, 2u> &mima) const {return IteratorMima(*this, mima);}

	template<class M> IteratorMimerger<M> getIterator(const Tuple< Tuple<STORETYPE,DIMS>, 2u> &mima, const M& mergescope) const {return IteratorMimerger<M>(*this, mima, mergescope);}

	template<class M> IteratorMerger<M> getIteratorMerger(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const M& mergescope) const {return IteratorMerger<M>(*this, in_box, mergescope);}
	template<class F> IteratorFunc<F> getIterator(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const F& func) const {return IteratorFunc<F>(*this, in_box, func);}
	template<class F, class M> IteratorFuncMerger<F,M> getIteratorMerger(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box, const F& func, const M& mergescope) const {return IteratorFuncMerger<F,M>(*this, in_box, func, mergescope);}

	HierarchicalMergable& toZero();
	HierarchicalMergable& toMemfree(){space.toMemfree();return *this;}


	unsigned int getSize()const{return space.getSize();}

    KeyIterator getKeyIterator()const{return KeyIterator(*this);}
    ERRCODE save(FILE *f) const {return space.save(f);}
    ERRCODE load(FILE *f) {return space.load(f);}


	void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node); // insert and erase

	void insertMK2(HyperCursor<STORETYPE,DIMS, LEAD> pos, NODES n_node, bool debug = false); // insert and erase

    // nodes are simple or complex, merging is allowed if equal, or possible is one of the two node is simple
	template<class MERGER> void insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node, const MERGER& merger); // insert and erase

    // nodes are simple or complex, merging is allowed if equal, or possible is one of the two node is simple

	// NODES has 2 functions that predicts that mergeInto/Simplify are guarrantied to return false in all circonstances:
    // bool canMerge() const;
    // bool canSimplify() const;
    // merger class has 2 functions:
    // bool mergeInto(const NODES& this, NODES& brother_io, const HyperCursor<STORETYPE,DIMS, LEAD> &this_loc) const;
    // bool simplifyAt(const NODES& this, NODES &fout, const HyperCursor<STORETYPE,DIMS, LEAD> &location) const; (return false if this == fout)
	template<class MERGER> void insertMK2(HyperCursor<STORETYPE,DIMS, LEAD> pos, NODES n_node, const MERGER& merger, bool need_simplification= false, bool debug = false); // insert and erase



	void wrGetAt(HyperCursor<STORETYPE,DIMS,LEAD>& _out, NODES& _out_data, const Tuple<STORETYPE,DIMS> &where) const; // at point, always pure
	bool getAt(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	bool getAt_simplify(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	bool getAt_simplifyMK2(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	template<class M> bool getAt_simplifyMK2(NODES& fout,const HyperCursor<STORETYPE,DIMS,LEAD> &where, const M& merger) const;
	bool getAt_update(NODES& fout,HyperCursor<STORETYPE,DIMS,LEAD> &where) const;
	template<class M> bool getAt_update(NODES& fout,HyperCursor<STORETYPE,DIMS,LEAD> &where, const M& merger) const;


	bool checkIntegrity()const;
	void insert_with_bypass(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES& n_node); // simple insert, may corrupt structure if overlap
	void show(FILE*f = stdout, int level =0)const{
	    return space.show(f,level);
	    }
};

template<class STORETYPE, class ALIASTYPE, unsigned int DIMS>
class ContinuousRangeQueries{
public:
    class EntryData{
    public:
        Tuple<STORETYPE, DIMS> pos;
        Tuple<float, DIMS> velo;
        Tuple<float, DIMS> accel;
        uint32_t time;

        float radius;
        float queryradius;
    };

    HeapTree<KeyElem<uint32_t, ALIASTYPE> > tocheck;
    myHashmap<ALIASTYPE, EntryData> entries;
    HierarchicalTree<STORETYPE, DIMS, 0, ALIASTYPE> queries;
    HierarchicalTree<STORETYPE, DIMS, 0, ALIASTYPE> objects;
};

template<unsigned int DIMS>
class GaussCloudNode{
public:
    double weight;
    Tuple<double, DIMS> center;
    Tuple<double, DIMS> box_radii;
    void show(FILE* f= stdout, int level=0) const;
};

template<class NODE, unsigned int DIMS, class INTCLASS>
class GaussCloud{
public:
    HierarchicalTree<INTCLASS, DIMS, 0, NODE> hie;
    myHashmap<NODE, GaussCloudNode<DIMS> > data;
    TMatrix<double> proj_matrix;
    double Zcover;

    GaussCloud(){}
    GaussCloud(const Vector< KeyElem<NODE, Tuple<double> > > &points);

    //void initFrom(const Vector< KeyElem<NODE, Tuple<double> > > &points);
    void setProjectionMatrix(const TMatrix<double>& dam){ proj_matrix = dam;}

    void insert(const NODE &node, const Tuple<double> &coor, const Trianglix<double> &var, double weight = 1.0);
    uint32_t findClosest(uint32_t query, const Trianglix<double, DIMS> &proj) const;
    uint32_t wrWithinRange(uint32_t* wr_buf, double threshold, uint32_t query, const Trianglix<double, DIMS> &proj) const;

    HyperPosition<INTCLASS, 0u,0u> makeBox(const GaussElem<Tuple<double, 0u> > datapoint, double LL_bound) const;

    Vector<NODE> queryEllipsoid(ThreadBase &tb, const Tuple<double> &center, const Trianglix<double> &radius) const;

    template<class R>  void inflateZcover(ThreadBase &tb, double newZcover, R& receiver, ACCEPTOR_DECL(R, pair<NODE MAKE_COMMA NODE> )) const;

	Forest<uint32_t,3u>  cluster_likelihood(const TMatrix< double > &data, double scale = 0.125f);

    void show(FILE* f=stdout, int level=0)const;
};

// much simpler... only an ordering
template<class NODE, class INTCLASS>
class GaussCloud<NODE, 1u,INTCLASS>{
public:
    TMatrix<double> proj_matrix;
    myHashmap<NODE, GaussCloudNode<1u> > data;
    RBTofDoom<KeyElem<double, NODE> > hierach;
    ERRCODE setProjectionMatrix(const TMatrix<double>& dam);
    void insert(const NODE &node, const Tuple<double> &coor, const Trianglix<double> &var, double weight = 1.0);
    void remove(const NODE &node);
};

/*
// silly! no ordering
template<class NODE, class INTCLASS>
class GaussCloud<NODE, 0u,INTCLASS>{
public:
    myHashmap<NODE> data;
    ERRCODE setProjectionMatrix(const TMatrix<double>& dam){return 0;}
    void insert(const NODE &node, const Tuple<double> &coor, const Trianglix<double> &var, double weight = 1.0){data.addEntry(node);}
    void remove(const NODE &node){data.removeEntry(node);}
};*/


/*
template<class D> class DoClusterLikelihoodRatioScope : public Event<uint32_t>{
public:
    class ParallelScope{
        uint32_t colrange[2];
    };
    GaussCloud<uint32_t, 3u> hie;
    Tuple<ParallelScope> scp;
    HeapTree< KeyElem<double, Tuple<int,2> > > pairs;
    Tuple<HeapTree< KeyElem<double, unsigned int > > > loc_merge;
    GaussElem< Tuple<D, 0u > >* stats;
    HeapTree< KeyElem<double, unsigned int > > merge;
    unsigned int* links;
    unsigned int* tom;
    uint32_t totsize;
    uint32_t midsize;
    DoClusterLikelihoodRatioScope(uint32_t nbthread = 4u);
    ~DoClusterLikelihoodRatioScope();
    DoClusterLikelihoodRatioScope& toMemfree();
    DoClusterLikelihoodRatioScope& toSize(uint32_t nsize);
    uint32_t operator()(uint32_t thread_no);
    template<class C,int NBREL,class TB> void run(Forest<C,NBREL>& target, const Vector< GaussElem< Tuple<double, 0u > > > &data, TB &tb);
};*/


template<class C>
class ClientPointer{ // uses man in middle to track the object
inline void clear(){if (t == NULL) return;}
	inline void set(C * const target){t = target;}
public:
	C* t;
	ClientPointer(): t(NULL){}
	ClientPointer(C * const target){}
	~ClientPointer() {clear();}

	const ClientPointer& operator=(C * const target){clear(); set(target);}
};


//\Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2))) \sum w_i (x_i - \Hat\mu)^2
//\Hat\sigma^2 = ((\sum w_i)/((\sum w_i)^2 - \sum (w_i^2))) \sum w_i x_i^2 - 2*x_i\Hat\mu + \Hat\mu^2
//\Hat\sigma^2 = (1/((\sum w_i)^2 - \sum (w_i^2))) (((\sum w_i) * \sum w_i x_i^2) - (\sum w_i x_i)^2)
//

//\Hat\mu = (\sum w_i x_i) / (\sum w_i)




// event: any change of state

// abstract class

enum Event_TYPE{
    EVTYPE_SYNC,
    EVTYPE_ASYNC,
    EVTYPE_LOCATION,
};
/*
template<class ARG, class KEY>
class Event{
public:
    typedef KEY INDEX_TYPE;
	virtual KEY getIndex() const=0;
    virtual unsigned int operator()(ARG)=0; // output 0 to be deleted 0xFFFFFFFF to be ignored, and anything else to be re-inserted (with a increment)
    virtual ~Event(){}
    };*/

/*
template<class SCOPE, class FUNCTION,class OUTPUT, class >
class FunctionCall
*/

template<class MSG>
class Listener{
public:
    virtual void listen(const MSG&)=0;
    virtual unsigned int getAlias()const{return 0;}
};



/** \brief Buffer for messages sent from X threads to a single thread that reads messages. The order of sent messages dictates the order in which they are read.
 *
 *
 */

template<class C>
class FIFOQueue{
    public:
    typedef void INDEX_TYPE;
    Event<C>* buffer[256];
    C scope;
    std::atomic<int> semaphore;
    int async_read;

    FIFOQueue(C _scope): scope(_scope),semaphore(0), async_read(0){}
    ERRCODE insert(Event<C>* ev);
    void runAll();
    void runAll(C _scope);
    bool isEmpty()const{return semaphore == async_read;}
};

template< >
class FIFOQueue<void>{
    public:
    typedef void INDEX_TYPE;
    Event<void>* buffer[256];
    std::atomic<int> semaphore;
    int async_read;

    FIFOQueue(): semaphore(0), async_read(0){}
    ERRCODE insert(Event<void>* ev){int pos = semaphore.fetch_add(1);
    if (((async_read - pos) & 0xFF) == 1){
        semaphore--;
        return 1; // maximum reached!
    }
    buffer[pos & 0xFF] = ev;
    return 0;}
    void runAll();
    bool isEmpty()const{return semaphore == async_read;}
};

/** \brief Buffer for messages sent from X threads to a single thread that reads messages. The order of sent messages dictates the order in which they are read.
 *
 *
 */
template<class C, class ARG, class R = C>
class OrderedTask{
public:
	C time;
	std::function<R(ARG)> task;
	OrderedTask()=default;
	OrderedTask(C _time, std::function<R(ARG)> _task): time(_time), task(_task){}
	bool operator>(const OrderedTask& other)const{return time > other.time;}
	bool operator<(const OrderedTask& other)const{return time < other.time;}
	bool operator>=(const OrderedTask& other)const{return time >= other.time;}
	bool operator<=(const OrderedTask& other)const{return time <= other.time;}
	OrderedTask& toMemmove(OrderedTask& other){time = other.time;task = other.task; other.task = std::function<R(ARG)>(); return *this;}
	OrderedTask& toMemswap(OrderedTask& other){std::function<R(ARG)> tmp = task; task = other.task; other.task= tmp;C tmp2 = time; time=other.time;other.time=tmp2; return *this;}
	//bool operator>(const OrderedTask& other)const{return time > other.time;}
	//bool operator>(const OrderedTask& other)const{return time > other.time;}
	void show(FILE *f =stdout, int level=0)const{fprintf(f,"Scheduled task at time:"); ExOp::show(time,f,level);}
};


// integer with overflow unaware comparisons
class CircularInteger{
public:
	int32_t value;
	CircularInteger()=default;
	CircularInteger(int32_t v): value(v){}
	operator int32_t(){return value;}
	CircularInteger& operator=(int32_t a){value = a;return *this;}
	CircularInteger& operator=(uint32_t a){value = *(int32_t*)&a;return *this;}
	CircularInteger& operator+=(int32_t a){value += a;return *this;}
	CircularInteger& operator-=(int32_t a){value -= a;return *this;}
	CircularInteger& operator+=(uint32_t a){value += a;return *this;}
	CircularInteger& operator-=(uint32_t a){value -= a;return *this;}
	inline bool operator>(const CircularInteger& other)const{return ((other.value - value) & 0x80000000) != 0;}
	inline bool operator<(const CircularInteger& other)const{return ((value - other.value) & 0x80000000) != 0;}
	inline bool operator>=(const CircularInteger& other)const{return ((value - other.value) & 0x80000000) == 0;}
	inline bool operator<=(const CircularInteger& other)const{return ((other.value - value) & 0x80000000) == 0;}
	inline bool operator>(const int32_t other)const{return ((other - value) & 0x80000000) != 0;}
	inline bool operator<(const int32_t other)const{return ((value - other) & 0x80000000) != 0;}
	inline bool operator>=(const int32_t other)const{return ((value - other) & 0x80000000) == 0;}
	inline bool operator<=(const int32_t other)const{return ((other - value) & 0x80000000) == 0;}
	inline bool operator!=(const CircularInteger& other)const{return value != other.value;}
	inline bool operator==(const CircularInteger& other)const{return value == other.value;}
	inline bool operator!=(const int32_t other)const{return value != other;}
	inline bool operator==(const int32_t other)const{return value == other;}
	void show(FILE*f = stdout, int level=0)const{fprintf(f, "%i", value); if (level == 0) fprintf(f, "\n"); }
};


/** \brief Buffer for messages sent from X threads to a single thread that reads messages. The order of sent messages dictates the order in which they are read.
 *
 *
 */
template<class C, class ARG, int ASYNC_MAG>
class PriorityQueue{
    public:
    typedef C INDEX_TYPE;
    typedef ARG INNER_TYPE;
    AsyncBuffer< OrderedTask<C, ARG>, ASYNC_MAG> async_buf;
    HeapTree<OrderedTask<C, ARG> > heap;
    PriorityQueue(){}
    ERRCODE insertSync(C priority, std::function<C(ARG)> ev);
    ERRCODE insertAsync(C priority, std::function<C(ARG)> ev);
    bool run(ARG& scope);

    void run_until(ARG& scope, C priority);
};

template<class C, int ASYNC_MAG>
class PriorityQueue<C, void, ASYNC_MAG>{
    public:
    typedef C INDEX_TYPE;
    typedef void INNER_TYPE;
    AsyncBuffer< OrderedTask<C, void> , ASYNC_MAG> async_buf;
    HeapTree< OrderedTask<C, void> > heap;
    ERRCODE insertSync(C priority, std::function<C()> ev);
    ERRCODE insertAsync(C priority, std::function<C()> ev);
    bool run();
    void run_until(C priority);
};

template<int ASYNC_MAG>
class PriorityQueue<CircularInteger, void, ASYNC_MAG>{
    public:
    typedef CircularInteger INDEX_TYPE;
    typedef void INNER_TYPE;
    AsyncBuffer< OrderedTask<CircularInteger, int32_t, int32_t> , ASYNC_MAG> async_buf;
    HeapTree< OrderedTask<CircularInteger, int32_t, int32_t> > heap;
    ERRCODE insertSync(int32_t time, std::function<int32_t(const int32_t)> ev);
    ERRCODE insertAsync(int32_t time, std::function<int32_t(const int32_t)> ev);
    bool run_until(int32_t);
    void show(FILE*f = stdout, int level=0)const;
};

template<class C, int ASYNC_MAG>
class PriorityQueue<CircularInteger, C, ASYNC_MAG>{
    public:
    typedef CircularInteger INDEX_TYPE;
    typedef void INNER_TYPE;
    AsyncBuffer< OrderedTask<CircularInteger, C, int32_t> , ASYNC_MAG> async_buf;
    HeapTree< OrderedTask<CircularInteger, C, int32_t> > heap;
    ERRCODE insertSync(int32_t time, std::function<int32_t(C)> ev);
    ERRCODE insertAsync(int32_t time, std::function<int32_t(C)> ev);
    bool run_until(C arg, int32_t);
    void show(FILE*f = stdout, int level=0)const;
};



/** \brief Buffer for messages sent from X threads to a single thread that reads messages. The order of sent messages dictates the order in which they are read.
 *
 *
 */

template<class KEY>
class StochasticQueue{
    void runAtIndex_routine(int index);

    uint16_t getIndex_from_exptime_routine(uint32_t expect_time) const;
    public:
    CategoryHashmap<KEY, Event<KEY>*, uint16_t> events;

    int32_t time;
    // Lambas are of the form 2^(i/8) where i is a integer


    double mainLambda;
    uint32_t lambdasum[1];


    StochasticQueue(int32_t start_time);

    const Event<KEY>* operator[](const KEY &key) const{return events[key];}
    Event<KEY>* operator[](const KEY &key){return events[key];}

    void updateTime(const KEY &key, uint32_t new_expect_time); /* time: 256 - 16777216 */
    void insert(const KEY key, Event<KEY>* ev, uint32_t expect_time); /* time: 256 - 16777216 */

    void remove(const KEY &key);
    void runTo(int32_t n_time);

	void show(FILE* f = stdout, int level =0) const;
};


template<class MSG>
class ListenerPointer{
public:
    virtual void Listen(MSG&) const=0;
};

template<class A, class MSG, class MSG_FILTER>
class MakeitListen : public ListenerPointer<MSG>{
public:
    A* target;
    MSG_FILTER msg_filter;
    MakeitListen(){}
    MakeitListen(A* i_target): target(i_target){}
    MakeitListen(A* i_target, const MSG_FILTER* _msg_filter): target(i_target) {if (_msg_filter) msg_filter = *_msg_filter;}

    MakeitListen(const  MakeitListen<A,MSG,MSG_FILTER>& other){memcpy(this, &other, sizeof(MakeitListen<A,MSG,MSG_FILTER>));}

    MakeitListen<A,MSG,MSG_FILTER>& toZero(){this->target = NULL; return(*this);}
    MakeitListen& operator=(const  MakeitListen& other){memcpy(this, &other, sizeof(MakeitListen<A,MSG,MSG_FILTER>)); return *this;}
    virtual void Listen(MSG&) const;
};

template<class A, class MSG>
class MakeitListen<A,MSG,void> : public ListenerPointer<MSG>{
public:
    A* target;
    MakeitListen(){}
    MakeitListen(A* i_target): target(i_target){}
    MakeitListen(const  MakeitListen<A,MSG,void>& other){memcpy(this, &other, sizeof(MakeitListen<A,MSG,void>));}

    MakeitListen<A,MSG,void>& toZero(){this->target = NULL; return(*this);}
    MakeitListen& operator=(const  MakeitListen& other){memcpy(this, &other, sizeof(MakeitListen<A,MSG,void>)); return *this;}
    virtual void Listen(MSG&) const;
};

template<class A>
class MakeitDie : public ListenerPointer<int>{
public:
    A* target;
    MakeitDie(A* i_target):target(i_target){}
    virtual void Listen(int&) const{delete(this->target);}
};

// prototype message
class MessageExample{ // example of a message class
    public:
    typedef unsigned int KEY;
    typedef unsigned int MSG_FILTER; // filter data, stored on listenerbase, can be void
    bool validFilter(const MSG_FILTER &){return(true);}
    // or "typedef void MSG_FILTER;" so that "bool validFilter()" is not needed
    unsigned int operator()(){return 0;} // dispatches the message, hence the Message may be a event too.
};

class ListenerExample{
    public:
    void listen(MessageExample &){printf("base listener used!, WTF!!\n");}
    };

template<class MSG>
class ListenerBase_Hashtable : myHashmap< typename MSG::MSG_KEY , Vector< MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> >  >{
    public:
    Vector< MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> > universal_listeners;
    template<class A> void addListener(A* , typename MSG::MSG_KEY what, typename MSG::MSG_FILTER * filter = NULL);
    template<class A> void removeListener(A* , typename MSG::MSG_KEY what, typename MSG::MSG_FILTER * filter = NULL);

    template<class A> void addUniversal_Listener(A*, typename MSG::MSG_FILTER * filter = NULL);
    template<class A> void removeUniversal_Listener(A* , typename MSG::MSG_FILTER * filter = NULL);



    void dispatch(MSG& msg);
};
/*
template<class MSG>
class ListenerBase : public Vector< AliasPtr< Listener<MSG> > >{
public:
    void dispatch(const MSG& msg);
};*/

// should remove this class
template<class MSG> // a unique listener per message ONLY, identified with MSG::KEY
class ListenerBase_Alias : myHashmap< typename MSG::KEY , MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER>  >{
    public:
    template<class A> void addListener(A*); // and static ListenerBase_Alias* getAliasBase();

    void dispatch(MSG& msg);
    void* getTarget(const typename MSG::KEY&) const;

    void show(FILE* f = stdout, int level =0) const;
};

template<class MSG, class Resol, unsigned int nbdims>
class Listener_Hyperspace{
    HyperPosition<Resol, nbdims, 0u> listen_volume;
    public:
    virtual void listen(const MSG&)=0;
};

template<class MSG, class RES, unsigned int DIM> // a unique listener per message ONLY, identified with MSG::KEY
class ListenerBase_Hyperspace{
    public:

    HierarchicalTree<RES, DIM, 0u,  Listener_Hyperspace<MSG, RES, DIM>* > listeners;

    void addListener(Listener_Hyperspace<MSG,RES, DIM>*, HyperPosition<RES, DIM, 0u> &where); // class A needs a member: (MSG::KEY alias) and static ListenerBase_Alias* getAliasBase();
    void dispatch(Tuple<RES,DIM> const &point, MSG& msg);
    void show(FILE* f = stdout, int level =0) const;
};

class ExecutePointer{
public:
    virtual void execute() const=0;
};


template<class A>
class MakeitRun : public ExecutePointer{
public:
    A* target;
    unsigned int time;
    MakeitRun(){}
    MakeitRun(A* i_target): target(i_target){}
    MakeitRun(const  MakeitRun<A>& other){memcpy(this, &other, sizeof(MakeitRun<A>));}

    MakeitRun<A>& toZero(){this->target = NULL; return(*this);}
    MakeitRun& operator=(const  MakeitRun<A>& other){memcpy(this, &other, sizeof(MakeitRun<A>)); return(*this);}
    virtual void execute() const{target->execute(time);}

    bool operator>(const MakeitRun &o)const { return (time - o.time-1) < 0x80000000;}
	bool operator>=(const MakeitRun &o)const{ return (time - o.time) <= 0x80000000;}
	bool operator<(const MakeitRun &o)const{ return (o.time - time-1) < 0x80000000;}
	bool operator<=(const MakeitRun &o)const{ return (o.time - time) <= 0x80000000;}
	bool operator==(const MakeitRun &o)const{ return time != o.time; }
	bool operator!=(const MakeitRun &o)const{ return time == o.time; }
};

// objects are destroyed or unloaded(and have a load method)
/*
class AliasBank{
    public:
    myHashmap< unsigned int , pair< void*, unsigned int> > pointers;
    myHashmap< unsigned int , Vector< ListenerPointer<int> > > auto_destroy;
};*/


class AliasHost{
   // friend template<C> class AliasPtr<C>;
    public:
    AliasHost(bool autoAssign = true);
 //   AliasHost(const AliasHost&);
 //   AliasHost& operator=(const AliasHost&);

    void setAHAlias(unsigned int alias =0); // set to random if no argument

    virtual ~AliasHost();
};


// unit which allows to uncover host using an alias, and to lock the host for exclusive use (in multi-threading)
/*
class Alias{
    public:
    mutable unsigned int alias;
    Alias(): alias(0){}
    ~Alias();
    template<class TARG> inline bool operator==(const AliasPtr<TARG>& other)const {return alias == other.alias;}
};
*/


template<bool RUNNING = LFH_HAS_RUNNINGTIME_STATISTICS>
class RunningStatistics{
public:
  inline void begin(const char* )const{}
  inline void end(const char* )const{}
  inline void report( bool exit_with_error)const{}
};

template< >
class RunningStatistics<true>{
public:
  Vector<KeyElem<unsigned int, unsigned int> > scope;
  myHashmap<string, unsigned int> procmap;
  Vector<WeightElem<double, 2> > runs;
  Vector<string> procnames;
  void begin(const char* label);
  void end(const char* label);
  void report(bool exit_with_error = false);
};

extern RunningStatistics< > lfhst_stats;

class Object_with_Attributes{
public:
	virtual int stauto_nbattributes() const;
	virtual double stauto_attributes(unsigned int which) const;
};

class StateAutomata_State{
public:
	unsigned int query;
    unsigned int output; // on success, may be a transition
    unsigned int next; // on failure, transition
    unsigned int flags;
    unsigned int indexes[15];
    double weights[16];
};

template<class HOST>
class StateAutomata{
    public:
    unsigned int runAuto(unsigned int state, unsigned int reward, Tuple<unsigned int, 0u> &scope);
};


/*

class StateAutomata_Scope{
public:
	void fillStateScope(const unsigned int state, vector<Object_with_Attributes*> &objs);
};

template<class Scope, class State>
class StateAutomata_state{
public:
	int query;
	unsigned int fieldselector[8];
};

template<class Scope, class State>
class StateAutomata{
public:
	void run(Scope &what);
	void runwithTrace(Scope &what, FILE *f);
};

    class emptyclass{
    public:
        const emptyclass& toRand(){return(*this);}
        const emptyclass& toZero(){return(*this);}
        const emptyclass& toOne(){return(*this);}
    };


//

class StateAutomata_State{
public:
	unsigned int queries[8];
	unsigned int nextstate[4]; // threshold, true, false,
	char weights[28]; // last 6 char are for mult selections
	char nextarg_weights[22 * 8];
	unsigned int selector;
//	unsigned int default_arguments[16]; // default input
//	unsigned int operation[16]; // first 8 -bit determines operation
//	unsigned int extra_arguments[32];
	void randomize();
	static int app_weight(char weight, int value);
};



// Actor->event_target_listener , (event, event_target)
// starts at state 0, fail increments the state by 1, sucess add a value to current state, adding 0 or overflowing 8 times finishes the state automata.


template<int statemag>
class StateAutomata{
public:
	StateAutomata_State states[1 << statemag];
	short values[16 << 4]; // input values, for the neural network
	StateAutomata CrossOver(const StateAutomata&) const;
	const  StateAutomata&  mutate();
	void randomize();

	// Class C must implement both
	// void operator()(unsigned int* o, unsigned int* i, int a); query for a values, with label i to o, if this applies i[0] may change the meaning of the other i[x]
	template<class C> void process(int state, unsigned char* out_input_state, const C & , int* init_arg) const;
	// same as process, but mutates illegal values!
	//void mutationFix(); // fixes
	// template<class C> void wise_process(int state, int* out_input_state, const & C, int* init_arg = NULL);
};

template<int statemag, int population_mag =0>
class StateAutomata_Population{
	public:
	StateAutomata<statemag>  ais[1<< population_mag];

	unsigned int fitcount[1<< population_mag];
	double fit[1<< population_mag];

	StateAutomata<statemag>& operator[](int i);

	void registerFitness(double fitness,int which, int incr, int nbchromo);

	void randomize();

	void save(char* path);
	ERRCODE load(char* path);

	};

template <class Node, int base_resol, int incr_resol, int nb_resol, class intclass, int nbdim>
class MultiResolHashTable : public Vector<Node> {
	public:
		Vector<Tuple<unsigned int, nb_resol+1> > partitions;
		map<Tuple<intclass,nbdim> , unsigned int> hash_ressource[nb_resol+1];

		void insert(Tuple<intclass,nbdim> key, Node what);
	};
*/

template <class IntClass, int nbdim, class node, int flag>
class Spatial{
	public:
		map<Tuple<IntClass,nbdim> , unsigned int> links;

	};
// convolution windows! mean for const construction

	template<int length,  const double &relativeFrequency> class LowPassFilterConvulutionWindow : public Tuple<double, length>{LowPassFilterConvulutionWindow();}; // class meant for const initialisation

/*
class InterruptHost{
    public:
    map<unsigned int, >
};

template<class EVENTHOST>
class Listener{
    public:

    Listener(typename EVENTHOST::KEY &query);
    virtual void listen(typename const EVENTHOST::DATA& data);
};*/







class variint{
	public:
	int piv;
	int var;
	int eval(int time = clock());
	void setValue(int value, int derivative =0, int time = clock());
};
/*
class Time{
public:
	int time;
	double offset;
	double speed;
	Time() : time(0), offset(0.0f), speed(0.0f){}
	void operator()();
};*/

template <class Key>
class defaultcomparator{
	int operator()( Key const a ,  Key const b) const{
			if (a > b) return(1);
			else return( a==b ? 0 : -1);
		}
	};

template <class Key,class Comparator = defaultcomparator<Key> >
class RBTree{
public:
	RBTreeNode<Key,Comparator>* root;
	RBTreeNode<Key,Comparator>* first;
	RBTreeNode<Key,Comparator>* last;

	Comparator comp;

	RBTree();
	virtual ~RBTree();
	RBTreeNode<Key,Comparator>* Pop();
	RBTreeNode<Key,Comparator>* Top();
	bool empty();

	RBTreeNode<Key,Comparator>* find(Key where);

	void Insert(RBTreeNode<Key,Comparator>* what);
	void Delete(Key where);
	void Delete(RBTreeNode<Key,Comparator>* what);
	void Remove(RBTreeNode<Key,Comparator>* what);

	void MakeDBlack(RBTreeNode<Key,Comparator>* what); // routine for removal of nodes
	void display(FILE* where);
	void display(FILE* where,RBTreeNode<Key,Comparator> *mark);
	void flush();
};

template <class Key,class Comparator>
class RBTreeNode{
public:
	RBTreeNode<Key,Comparator>* l;
	RBTreeNode<Key,Comparator>* r;
	RBTreeNode<Key,Comparator>* p;
	RBTreeNode<Key,Comparator>* b;
	RBTreeNode<Key,Comparator>* n;
	bool red;
	RBTreeNode();
	virtual Key getIndex()=0;
	void recursiveDisplay(FILE* where,char* buffer);
	void recursiveDisplay(FILE* where,char* buffer,RBTreeNode<Key,Comparator> *mark);
	void clear();
};






//#include "Primitive_spatial.h"

/* Priority queue

class Event : public RBTreeNode<int,IntComparator>{
public:
	Event();
	virtual ~Event();
int getIndex();
virtual void operator()()=0; // stuff to run
virtual int getTime()=0; // time of insertion and confirmation
virtual bool autoDelete();
};







class PriorityQueue{
public:
//	PriorityQueue *metaqueue; // there is another queue, which allows this queue to work.
	RBTree<int,IntComparator> list;
	Time time;

	vector<Event*> toinsert;
	bool sema;

	PriorityQueue();
	void Insert(Event* ev);
	void dotilldone();
	void do_one();
	void start();
	void pause();
};*/

#define POLYTHING_DESIRED_ORDER 4
#define POLYTHING_DESIRED_ORDER_FACTORIAL 24
#define POLYTHING_DESIRED_ITERATION_NB 1000


class abstractMeshFace{
public:
	unsigned short points[3];
	unsigned short adjacent[3];
	int tag;
};
class abstractMesh{
public:
	vector<abstractMeshFace> mesh;
	unsigned short nbpts;
	abstractMesh();
	void insertTriangle(unsigned short a, unsigned short b, unsigned short c);
	void insertNew2dGrid(unsigned short x, unsigned short y);

	void insertNew2dGridHole(unsigned short x, unsigned short y, unsigned short hole);

	void fillajacents();
//	void pathFinding(PathFinderQueryFunction* func);
};

template<class TARG> class AliasPtr{
public:
    unsigned int alias;
    AliasPtr(): alias(0){}
    AliasPtr(const AliasPtr<TARG>& );
    AliasPtr(const TARG * );
    AliasPtr<TARG>& operator=(const AliasPtr<TARG>&);
    AliasPtr<TARG>& operator=(const TARG* );
    ~AliasPtr(){clear();}

    void clear();

    inline bool operator==(const AliasPtr<TARG>& other)const {return alias == other.alias;}
//    inline bool operator==(const Alias& other)const {return alias == other.alias;}

    const TARG* operator->() const; // read-only ptr, use with care
    const TARG* operator*() const; // read-only ptr, use with care
    void freePointer() const;

    AliasPtr<TARG>& memmove(AliasPtr<TARG>& from );

    bool isValid()const;
    // template<class HOST> void setDestroyOnClear(HOST* what); // if object pointed by alias is destroyed, destroy host as well!
};
template<class TARG, bool READONLY = false> class LockPtr{
    void initAlias(const unsigned int &alias); // called by contruction
    public:
    TARG* target; // locks target (cannot be moved or destroyed)
    LockPtr(): target(NULL){}
    LockPtr(const AliasPtr<TARG>&);
    LockPtr(const unsigned int&);

    ~LockPtr();
    template<class NTYPE> operator NTYPE* ()const;
    bool isValid()const;
    TARG* operator->() const;
    TARG& operator*() const;
};
template<class TARG> class LockPtr<TARG, true>{
    void initAlias(const unsigned int &alias); // called by contruction
    public:
    const TARG* target; // locks target (cannot be moved or destroyed)
    LockPtr(): target(NULL){}
    LockPtr(const AliasPtr<TARG>&);
    LockPtr(const unsigned int&);
    ~LockPtr();
    template<class NTYPE> operator const NTYPE* ()const;
    bool isValid()const;
    const TARG* operator->() const;
    const TARG& operator*() const;
};
template< > class LockPtr<void, true>{
    void initAlias(const unsigned int &alias); // called by contruction
    public:
    const void* target; // locks target (cannot be moved or destroyed)
    LockPtr(const AliasPtr<void>&);
    LockPtr(const unsigned int&);
    ~LockPtr();
    template<class NTYPE> operator const NTYPE* ()const;
    bool isValid()const{return target != NULL;}
    const void* operator->() const{return target;}
};
template< > class LockPtr<void, false>{
    void initAlias(const unsigned int &alias); // called by contruction
    public:
    void* target; // locks target (cannot be moved or destroyed)
    LockPtr(const AliasPtr<void>&);
    LockPtr(const unsigned int&);
    ~LockPtr();
    template<class NTYPE> operator NTYPE* ()const;
    bool isValid()const{return target != NULL;}
    void* operator->() const{return target;}
};
template<class KEY, class KEYHASH = defaultHashFnc<KEY> > class HashofAliases{
public:
    LFHPrimitive::myHashmap<KEY, unsigned int, KEYHASH> alias_hash;
    template<class C> LockPtr<C, false> operator[](const KEY &);
    template<class C> LockPtr<C, true> operator[](const KEY &);
};

// a ressource of type A must have:
// typedef ?ENUM? ressource_enum; // defineding the enum used
// A(?ENUM?) a constructor which uses that enum alone
// a function getManager() returning the address of a ResManager<A>
// a function randomKey() returning ressource_enum available for custom uses


template<class RES>
class RessourcePtr{
public:

    typedef RessourcePtr<RES> SAFETYPE;
	RES* target;
	RessourcePtr(): target(NULL) {}
	RessourcePtr(typename RES::ressource_enum const ID) : target(RES::getManager().get(ID)){}
    RessourcePtr(typename RES::ressource_enum const ID, unsigned int readonly_handdlechannel) : target(RES::getManager().get(ID, readonly_handdlechannel)){}
	RessourcePtr(RessourcePtr<RES> const &other) : target(other.target){if (target != NULL) RES::getManager().increment(target);}
	~RessourcePtr()     {if (target != NULL) RES::getManager().free(target); }
    RessourcePtr<RES>& operator=(typename RES::ressource_enum const ID) {ResManager<RES, RES::ressource_appends >& mana =RES::getManager();if (target != NULL) mana.free(target);target =  mana.get(ID);return *this;}
	RessourcePtr<RES>& operator=(RES* n_target)     {ResManager<RES, RES::ressource_appends>& mana =RES::getManager();if (target != NULL) mana.free(target);target = n_target; if (n_target != NULL) mana.put(n_target);return *this;}
	RessourcePtr<RES>& operator=(const RessourcePtr<RES>& n_target)     {if (target == n_target.target) return *this; ResManager<RES, RES::ressource_appends>& mana =RES::getManager();if (target != NULL) mana.free(target);target = n_target.target; if (target != NULL) RES::getManager().increment(target); return *this;}
	RessourcePtr<RES>& setTo(typename RES::ressource_enum const ID, unsigned int readonly_handdlechannel){ResManager<RES, RES::ressource_appends >& mana =RES::getManager();if (target != NULL) mana.free(target); target = RES::getManager().get(ID, readonly_handdlechannel);return *this;}

	inline operator bool()const     {return target != NULL; }
	inline RES* operator->()const     {return target; }
	inline RES& operator*()const     {return *target; }
	typename RES::ressource_enum getID()const      {return(RES::getManager().getId(target));}
    RessourcePtr<RES>& toZero(){if (target != NULL) RES::getManager().free(target); target = NULL; return *this;}
    RessourcePtr<RES>& toMemmove(RessourcePtr<RES>& other){target = other.target; other.target = NULL; return *this;}
    void show(FILE*f = stdout, int level=0) const{if (target == NULL) fprintf(f,"NULL%c", level == 0 ? '\n' : '\t'); else fprintf(f,"%X%c",(int)(intptr_t)target , level == 0 ? '\n' : '\t'); }
};

/*
template<class RES>
class RessourcePtr{
public:
    typedef RessourcePtr<RES> SAFETYPE;
	RES* target;
	RessourcePtr(): target(NULL) {}
	RessourcePtr(typename RES::ressource_enum const ID) : target(RES::manager.get(ID)){}
    RessourcePtr(typename RES::ressource_enum const ID, unsigned int readonly_handdlechannel) : target(RES::manager.get(ID, readonly_handdlechannel)){}
	RessourcePtr(RessourcePtr<RES> const &other) : target(other.target){if (target != NULL) RES::manager.increment(target);}
	~RessourcePtr()     {if (target != NULL) RES::manager.free(target); }
    RessourcePtr<RES>& operator=(typename RES::ressource_enum const ID) {ResManager<RES, RES::ressource_appends >& mana =RES::manager;if (target != NULL) mana.free(target);target =  mana.get(ID);return *this;}
	RessourcePtr<RES>& operator=(RES* n_target)     {ResManager<RES, RES::ressource_appends>& mana =RES::manager;if (target != NULL) mana.free(target);target = n_target; if (n_target != NULL) mana.put(n_target);return *this;}
	RessourcePtr<RES>& operator=(const RessourcePtr<RES>& n_target)     {ResManager<RES, RES::ressource_appends>& mana =RES::manager;if (target != NULL) mana.free(target);target = n_target.target; if (target != NULL) RES::manager.increment(target); return *this;}
	RessourcePtr<RES>& setTo(typename RES::ressource_enum const ID, unsigned int readonly_handdlechannel){ResManager<RES, RES::ressource_appends >& mana =RES::manager;if (target != NULL) mana.free(target); target = RES::manager.get(ID, readonly_handdlechannel);return *this;}

	operator bool()const     {return target != NULL; }
	RES* operator->()const     {return target; }
	typename RES::ressource_enum getID()const      {return(RES::manager.getId(target));}
    RessourcePtr<RES>& toMemmove(RessourcePtr<RES>& other){target = other.target; other.target = NULL; return *this;}
    void show(FILE*f = stdout, int level=0) const{if (target == NULL) fprintf(f,"NULL%c", level == 0 ? '\n' : '\t'); else fprintf(f,"%X%c",(intptr_t)target , level == 0 ? '\n' : '\t'); }
};*/

// Constructs an object with additionnal memory allocated after
// constructor
template<class A>
class AppendixPtr{
   // friend class ResManager<A,1u>;
    public:
    typedef AppendixPtr<A> SAFETYPE;
    char* data;
    AppendixPtr();
    ~AppendixPtr();
    AppendixPtr<A>& setExtraSize(unsigned int size);
    AppendixPtr<A>& toZero();
    AppendixPtr<A>& toMemfree();
    AppendixPtr<A>& toMemmove(AppendixPtr&);
    AppendixPtr<A>& operator=(const AppendixPtr&);
    A* operator->()const {return (A*) data;}
    void show(FILE* f = stdout, int level =0)const{((A*) data)->show(f,level);}
};

// a ressource of type A must have:
// typedef ?ENUM? ressource_enum; // defineding the enum used
// const int ressource_appends = 0 or 1; // allows the variable allocation
// a static function ResManager<?ENUM?, ressource_appends> getManager();
// a static FILE* getFileHanddle(?ENUM? id); // get a FILE* to data of interest (*wont close file*)
// a function void finisgload(FILE* f, unsigned int varsize)const;
// a static function ?ENUM? randomKey() returning ressource_enum available for custom uses



template<class RES, unsigned int ressource_appends>
class ResManager{
	myHashmap<EnumBox<typename RES::ressource_enum>, RES* > hash_ressource;
	myHashmap< RES* , pair<EnumBox<typename RES::ressource_enum>, unsigned int> > back_hash;
    typename RES::ressource_enum buffer[256];
    unsigned short buffer_indexes[2];
	public:
    ResManager(){buffer_indexes[0] = 0; buffer_indexes[1] =0;}
    RES* unsafe_get(typename RES::ressource_enum const ID) const{return hash_ressource[ID];} // crashes if not loaded, ressource can be deleted if not owned by at least 1 real ressource pointer
	RES* get(typename RES::ressource_enum const ID); // get ressource, and load it with default function if needed
	RES* get(typename RES::ressource_enum const ID, unsigned int init_flag); // get ressource, and load it with default function if needed

	void increment(RES* res){back_hash[res].second++;}
	typename RES::ressource_enum getId(RES* res); // get ID for a loaded ressource
	void free(RES* res ); // frees a loaded ressource
	typename RES::ressource_enum put(RES* res); // stores a custom ressource, needs an new custom ID!
	void put_at(RES* res, typename RES::ressource_enum key); // stores a custom ressource, needs an new custom ID!
    void show(FILE* f = stdout, int level =0)const{hash_ressource.show(f,level);}
};



template<class RES>
class ResManager<RES,1u>{
	myHashmap<EnumBox<typename RES::ressource_enum> , AppendixPtr<RES> > hash_ressource;
	myHashmap< void* , pair<EnumBox<typename RES::ressource_enum>, unsigned int> > back_hash;
    typename RES::ressource_enum buffer[256];
    unsigned short buffer_indexes[2];

	public:
    ResManager(){buffer_indexes[0] = 0; buffer_indexes[1] =0;}
    RES* unsafe_get(typename RES::ressource_enum const ID) const{return (RES*)(hash_ressource[ID].data);} // crashes if not loaded, ressource can be deleted if not owned by at least 1 real ressource pointer

	RES* get(typename RES::ressource_enum const ID); // get ressource, and load it with default function if needed
	RES* get(typename RES::ressource_enum const ID, unsigned int init_flag); // get ressource, and load it with default function if needed
	void increment(RES* res){back_hash[(char*)res].second++;}
	typename RES::ressource_enum getId(RES* res); // get ID for a loaded ressource
	void free(RES* res ); // frees a loaded ressource
	//typename RES::ressource_enum put(RES* res); // stores a custom ressource, needs an new custom ID!
	void put_at(RES* res, typename RES::ressource_enum key); // stores a custom ressource, needs an new custom ID!
    void show(FILE* f = stdout, int level =0)const{hash_ressource.show(f,level);}
	};

template<class RES>
class ResManager<RES,2u>{
	myHashmap<EnumBox<typename RES::ressource_enum>, RES* > hash_ressource;
	myHashmap< RES* , pair<EnumBox<typename RES::ressource_enum>, unsigned int> > back_hash;
    typename RES::ressource_enum buffer[256];
    unsigned short buffer_indexes[2];

	public:
    ResManager();
    RES* unsafe_get(typename RES::ressource_enum const ID) const{return hash_ressource[ID];} // crashes if not loaded, ressource can be deleted if not owned by at least 1 real ressource pointer
	RES* get(typename RES::ressource_enum const ID); // get ressource, and load it with default function if needed
	RES* get(typename RES::ressource_enum const ID, unsigned int init_flag); // get ressource, and load it with default function if needed
    void* ressourceArgs;

	void increment(RES* res){back_hash[res].second++;}
	typename RES::ressource_enum getId(RES* res); // get ID for a loaded ressource
	void free(RES* res ); // frees a loaded ressource
	typename RES::ressource_enum put(RES* res); // stores a custom ressource, needs an new custom ID!
	void put_at(RES* res, typename RES::ressource_enum key); // stores a custom ressource, needs an new custom ID!

    void show(FILE* f = stdout, int level =0)const{hash_ressource.show(f,level);}
};


template<class RES>
class ResManagerFileHanddle{
	myHashmap<EnumBox<typename RES::ressource_enum> , AppendixPtr<RES> > hash_ressource;
	myHashmap< void* , pair<EnumBox<typename RES::ressource_enum>, unsigned int> > back_hash;
    SerialStore<typename RES::ressource_enum> fileheap;
    public:
    ResManagerFileHanddle(){}
    RES* unsafe_get(typename RES::ressource_enum const ID) const{return (RES*)(hash_ressource[ID].data);} // crashes if not loaded, ressource can be deleted if not owned by at least 1 real ressource pointer

	RES* get(typename RES::ressource_enum const ID); // get ressource, and load it with default function if needed
	RES* get(typename RES::ressource_enum const ID, unsigned int init_flag); // get ressource, and load it with default function if needed
	void increment(RES* res){back_hash[(char*)res].second++;}
	typename RES::ressource_enum getId(RES* res); // get ID for a loaded ressource
	void free(RES* res ); // frees a loaded ressource
	typename RES::ressource_enum put(RES* res); // stores a custom ressource, needs an new custom ID!
    void show(FILE* f = stdout, int level =0)const{hash_ressource.show(f,level);}
};


// RES needs a typename KEY_CLASS for the KEY class which maintains the access to the manager instance
template<class RES>
class EntryPtr{
public:
	RES* target;
	EntryPtr(): target(NULL) {}
	EntryPtr(typename RES::KEY_CLASS const ID) {
		EntryManager< typename RES::KEY_CLASS >& mana = RES::getManager();
		target = mana.get(ID);
	}
	EntryPtr(EntryPtr< RES> const &other) : target(other.target){RES::getManager().increment(target);}
	~EntryPtr()     {if (target != NULL) RES::getManager().free(target); }
	void operator=(typename RES::KEY_CLASS const ID)     {EntryManager< typename RES::KEY_CLASS >& mana = RES::getManager(); if (target != NULL) mana.free(target);target =  mana.get(ID);}
	void operator=(RES* n_target)     {EntryManager< typename RES::KEY_CLASS >& mana = RES::getManager();if (target != NULL) mana.free(target);target = n_target; if (n_target != NULL) mana.put(n_target);}
	operator bool()const     {return target != NULL; }
	RES* operator->()const     {return target; }
	typename RES::ressource_enum getID()const {EntryManager< typename RES::KEY_CLASS >& mana = RES::getManager(); return(mana.getId(target));}
};


// RES needs a typename KEY defining the access to the manager instance
template<class KEY>
class EntryManager{
	map< KEY , void* > hash_ressource;
	map< void* , pair< KEY , unsigned int> > back_hash;
public:
	void* get(KEY const ID); // get ressource, or uses KEY to load it
	void increment(void* res){back_hash.find(res)->second.second++;}
	KEY getId(void* res); // get ID for a loaded ressource
	void free(void* res ); // frees a loaded ressource
	KEY put(void* res); // stores a custom ressource, needs an new custom ID!
};
template <class T>
class Ressource{
public:
	Ressource* next;
	Ressource* prev;
	T target;

	Ressource();
	Ressource(T);
};
class BaseManager {
public:
	int maxitem;
	int nbitem;
	bool useTag;
	Data* owner;
    virtual ~BaseManager(){}
};
template <class E, class T>
class DataManager : public BaseManager{
public:
	map<E, Ressource<T> > hash_ressource;

	DataManager();
	DataManager(int capacity);
	~DataManager();
	E getFreeKey(); // for int enums only, start with 0xFFFFFFFF downwards
	T getRessource(E key); // get it if available, load and free ressources it if necessessary
	T quickgetRessource(E key); // get it if available only, return NULL otherwise
	virtual Ressource<T> loadRessource(E key);
	void insertRessource(E key, T data);
	void freeRessource(E key);
	bool hasRessource(E key); // check if it is there, nothing else
	void flushall();
};
class Data{
public:
	vector<BaseManager*> managelist;
	bool creatinglist;
/*
	map<texture, int> textures;
	priority_queue<TextRessource, vector<TextRessource>, RessourceOrder> texturesOrder;
	map<texture, AnimRessource*> anims;
	priority_queue<AnimRessource*, vector<AnimRessource*>, RessourceOrder> animsOrder;
*/
	Data();
	void addManager(BaseManager*);
//	void draw(_Display_enum_Mesh what);
//	void setTexture(_Display_enum_Texture what,int slot);
//	unsigned int loadTexture(char* path,int size);
//	short* anim(Anim what);


};
/*
class HeapStore{
	void change_nb_records(int variation);
	void box_insert(int logsize); // makes space at the end of the array of logsiz
	void box_delete(int logsize); // removes space at the end of the array of logsize

	void box_enlarge(int startlogsize, int endlogsize);
	void box_shrink(int startlogsize, int endlogsize);
	FILE *f;
	// data stored in file
	int sizes[16];

	// data maintained for speed
	int nbrecords;
	int quickoffsets[16];
	int shiftoffsets[16];
public:
	HeapStore(char* path);
	~HeapStore();
	void clear();
	void erase(unsigned int saveId);
	void save(unsigned int saveId, int charsize, void* data);
	void savenew(unsigned int &saveId,int charsize, void* data);
	void load(unsigned int saveId, void* &data);

//	LFHDisplay::SerialStruct* (*initialiser)(int myclass);
//	void load(vector<int> &which, vector<LFHDisplay::SerialStruct*> &objs);
//	void update(vector<LFHDisplay::SerialStruct*> &objs);
};*/
class RemoteLink{
public:
	bool master;
	void broadcast(void* message, int length);
	void digest(void* message, int length);
};
class RKenvironnement{
		int step;
		int offset;
	friend class RKable;
	public:
		double timestep;
		RKenvironnement(): step(0), offset(0), timestep(1.0f){}
		void nextstep(){step =(step +1) & 3; offset = (step) & (4 - step);}
	};
class RKable{
	double value[3]; // current value, cmputed derivative
	friend class RKenvironnement;
	public:

		void setValue(double val){ // RK state MUST be 0
			value[0] = val;
		}
		double eval(RKenvironnement* RK = NULL){
			if (RK == NULL) return(value[0]);
			else return(value[RK->offset]);
		}
		void registerDerivative(double deriv, RKenvironnement* RK){
			double tmp;
			switch(RK->step){
				case 0:
					value[1] = value[0] + RK->timestep * deriv / 2;
				break;
				case 1:
					value[2] = value[0] + RK->timestep * deriv / 2;
				break;
				case 2:
					tmp = value[0];
					value[0] = value[1] - value[0];
					value[1] = tmp + RK->timestep * deriv;
				break;
				case 3:
					value[0] = (value[0] + value[2] * 2 + value[1] + RK->timestep * deriv / 2)/3;
				break;
				}
		}
	};

class OwnerRegisterCmp{
	public:
	int operator()(void*  const a , void* const b) const;
	};
class OwnerRegisterEntry: public RBTreeNode<void*, OwnerRegisterCmp>{
	public:
	void* key;
	void* endkey;
	void* target;
	OwnerRegisterEntry(){}
	OwnerRegisterEntry(void *owner, void* start, void* end);
	void* getIndex();
	};
	/*
class LinkEntry{
	public:
	int count; // negative if owned!
	void* link;

	LinkEntry();
	LinkEntry(void* _link, int _count);
	};*/


// hold a back pointer for to the onwer of an memory area ( array ),
// caches thing that are reused
/*
class LinkRegister: public RBTree<void*, OwnerRegisterCmp>{
	public:

		OwnerRegisterEntry cachedOwners[16];
		int nbcached;


	//	map<unsigned int, LinkEntry > ID_to_ADD;
	//	map<unsigned int, void* > ID_to_ADD2;
	//  map<void*, unsigned int > ADD_to_ID;

		LinkRegister(){nbcached =0;}


        //unsigned int getID(void* address);
		//void* getLink(unsigned int &virtID);
		//void incrLink(unsigned int &virtID);
		//void decrLink(unsigned int &virtID);
		//unsigned int incrgetID(void* address);

		void* removeOwner(void* start);
		void addOwner(void *owner, void* start, void* end); // assumes target has to owner
		void setOwner(void *owner, void* start, void* end); // assumes nothing

		void* operator()(void* address) const; // return link
		void* moveOwner(void* oldstart, void* newstart, void* newend);
	};*/

template <class Host, class Target, unsigned int Hoff, unsigned int Toff>
class DLink{
	Target* tar;
	public:
	DLink(): tar(NULL){}
	DLink(const DLink<Host, Target, Hoff, Toff> &){}
	Target* operator->()const{return(tar);}
	DLink<Target,Host,Toff,Hoff>& getOther()const{return( *((DLink<Target,Host,Toff,Hoff>*)(((char*)tar) + Target::consint[Toff])) );}
	DLink& memmove(DLink& source){
		tar = source.tar;
		source.tar = NULL;
		DLink<Target,Host,Toff,Hoff>& backlink = getOther();
		backlink.tar = (Host*) (((char*)tar) - Host::consint[Hoff]);
		return(*this);
		}


};

template <class C>
class Ptr{
	public:
	unsigned int Id;

	Ptr();
	~Ptr();
	C& operator->();
	C& operator*();

	operator C*();
	C* operator=(C const * const nval);
	C* operator=(Ptr<C> const & nval);
};

	// assumes that the class C extends WTarget!





// template<class C,int size> template <> Tuple<C,size>& Tuple<C,size>::operator+=(Tuple<C,size> const & other);

	/*
class movie{
	public:
		int sizex,sizey,nbframes;
		char* data;
		movie();
		void saveAVI(char* path);
	};*/
//class ThreeDState{
//public:
//	double pos[3];
//	SQuaternion r;
//	ThreeDState interpolate(float fraction, ThreeDState other);
//	ThreeDState interpolate(float timeoffset, float statetimedistance, ThreeDLinearState other);
//};

/*
class ThreeDLinearState{
public:
	double pos[3];
	SQuaternion r;
	double speed[3];
	double angspeed[3];
	ThreeDState interpolate(float timeoffset, float statetimedistance, ThreeDLinearState other);
};*/

// templates specializations:

	/*
template <class C, int nbdim>
class HyperPositionCube : public IntervalSet<HyperPosition<C,nbdim> >{
	public:
	HyperPosition<C,nbdim> min;
	HyperPosition<C,nbdim> max;
	HyperPositionCube(const HyperPosition<C,nbdim> & _min, const HyperPosition<C,nbdim> & _max): min(_min), max(_max){}
	const SetComparison& compare(const HyperPosition<C,nbdim> &other) const;
	const SetComparison& compareInterval(const HyperPosition<C,nbdim> & _min, const HyperPosition<C,nbdim> & _max)const;
	};*/

template <class C, int nbdim>
class HyperPositionProjectCube: public IntervalSet<HyperPosition<C,nbdim> > {
	public:
	HyperPosition<C,nbdim> min;
	C diam;
	const SetComparison& compare(const HyperPosition<C,nbdim> &other) const;
	const SetComparison& compareInterval(const HyperPosition<C,nbdim> & min, const HyperPosition<C,nbdim> & max)const;
	};


template <class Var = int, bool hasone = false>
class mantissa{
	public:
	Var data;

	operator double() {return((double)data);}

	mantissa<Var,hasone>& operator=(Var const & val){data = val; return(*this);}
	mantissa<Var,hasone>& operator=(double const & val){data = val; return(*this);}

	mantissa<Var>& operator+=(mantissa<Var> const & other) {data += other.data; return(*this);}
	mantissa<Var>& operator-=(mantissa<Var> const & other) {data -= other.data; return(*this);}

	template <class OVar> mantissa<Var>& operator+=(mantissa<OVar> const & other) {data = ((double)(*this) + ((double)other.data)); return(*this);}
	template <class OVar> mantissa<Var>& operator-=(mantissa<OVar> const & other) {data = ((double)(*this) - ((double)other.data)); other.data; return(*this);}
	template <class OVar> mantissa<Var>& operator*=(mantissa<OVar> const & other) {data = ((double)(*this) * ((double)other.data)); other.data; return(*this);}
	template <class OVar> mantissa<Var>& operator/=(mantissa<OVar> const & other) {data = ((double)(*this) / ((double)other.data)); other.data; return(*this);}

	LFHDECL_TRIVIAL_OPERATOR(+,mantissa<Var>)
	LFHDECL_TRIVIAL_OPERATOR(-,mantissa<Var>)
	LFHDECL_TRIVIAL_OPERATOR(*,mantissa<Var>)
	LFHDECL_TRIVIAL_OPERATOR(/,mantissa<Var>)
};



template < > mantissa<char,true>::operator double();
template < > mantissa<unsigned char,true>::operator double();
template < > mantissa<short,true>::operator double();
template < > mantissa<unsigned short,true>::operator double();
template < > mantissa<int,true>::operator double();
template < > mantissa<unsigned int,true>::operator double();

template < > mantissa<char,false>::operator double();
template < > mantissa<unsigned char,false>::operator double();
template < > mantissa<short,false>::operator double();
template < > mantissa<unsigned short,false>::operator double();
template < > mantissa<int,false>::operator double();
template < > mantissa<unsigned int,false>::operator double();

template < > mantissa<char,true>& mantissa<char,true>::operator=(double const & v);
template < > mantissa<unsigned char,true>& mantissa<unsigned char,true>::operator=(double const & v);
template < > mantissa<short,true>& mantissa<short,true>::operator=(double const & v);
template < > mantissa<unsigned short,true>& mantissa<unsigned short,true>::operator=(double const & v);
template < > mantissa<int,true>& mantissa<int,true>::operator=(double const & v);
template < > mantissa<unsigned int,true>& mantissa<unsigned int,true>::operator=(double const & v);

template < > mantissa<char,false>& mantissa<char,false>::operator=(double const & v);
template < > mantissa<unsigned char,false>& mantissa<unsigned char,false>::operator=(double const & v);
template < > mantissa<short,false>& mantissa<short,false>::operator=(double const & v);
template < > mantissa<unsigned short,false>& mantissa<unsigned short,false>::operator=(double const & v);
template < > mantissa<int,false>& mantissa<int,false>::operator=(double const & v);
template < > mantissa<unsigned int,false>& mantissa<unsigned int,false>::operator=(double const & v);

// y = F(x) * F(i) + F(j) * F( u + g[y] )




// group of orthogonal transformations in 3D (48 bases)

template<int nbdim>
class Parobject{
	public:
	virtual const mycomplex& eval(Tuple<mycomplex, nbdim> where) const =0;
	virtual const Tuple<mycomplex, nbdim>& evald(Tuple<mycomplex, nbdim> where) const =0;
	virtual const TMatrix<mycomplex, nbdim>& evaldd(Tuple<mycomplex, nbdim> where) const =0;
};

template<int nbdim>
class Sphere : public Parobject<nbdim>{
	public:
	const mycomplex& eval(Tuple<mycomplex, nbdim> where) const;
	const Tuple<mycomplex, nbdim>& evald(Tuple<mycomplex, nbdim> where) const;
	const TMatrix<mycomplex, nbdim>& evaldd(Tuple<mycomplex, nbdim> where) const;
};
/*
class applyfunction{
public:
	functionname which;
	applyfunction(functionname _op) : which(_op){}
	void operator()(double &what) const{
		switch(which){
			case FUNCNAME_NULL:what = 0.0f;break;
			case FUNCNAME_SQUARRE:what = what*what;break;
			case FUNCNAME_INVERSE:what = (what == 0.0f) ? 0.0f: 1 / what ;break;
			case FUNCNAME_SIN: what = cos(what); break;
			case FUNCNAME_COS: what = sin(what); break;
			case FUNCNAME_LNGAMMA: what = lngamma(what);break;
			case FUNCNAME_POLYGAMMA_0: what = polygamma0(what);break;
			case FUNCNAME_SINH: what = sinh(what);break;
			case FUNCNAME_COSH: what = cosh(what);break;
			case FUNCNAME_BESSEL_J0: what = BesselJ0(what);break;
			case FUNCNAME_BESSEL_J1: what = BesselJ1(what);break;
			case FUNCNAME_BESSEL_I0: what = BesselI0(what);break;
			case FUNCNAME_BESSEL_I1: what = BesselI1(what);break;
		}
	}
};*/

/*
template<unsigned int TSIZE, class C>
class RankNormalizationScope : public Tuple<HeapTree<C>, TSIZE>{
public:

	 Vector<Tuple<unsigned int, TSIZE> > regist(const Vector<Tuple<C, TSIZE> >&);

};*/




// Y = exp(x)

//exp(x) * exp(x) = exp(2x);


/*
template<class C>
class pointMax : public Oper2< C , C >{
	 public:
	 virtual void operator()(C const & input, C &) const;
};

template<class C, class ARG>
class pointArgMax : public Oper2< C ,pair<C,ARG> >{
	 public:
	 ARG argument;
	 pointArgMax(ARG aaa): argument(aaa){}
	 virtual void operator()(C const & input, pair<C,ARG> &) const;
};

template<class C, int size>
class Tuple_FiniteDifference : public Oper2< Tuple<C, size>, Tuple<C, size-1> >{
public:
	virtual void operator()(Tuple<C, size> const & input,Tuple<C, size-1> &) const;

};

template<class C, int size>
class Tuple_FiniteDifferenceNORM : public Oper2< Tuple<C, size>, Tuple<C, size-1> >{
public:
	virtual  void operator()(Tuple<C, size> const & input, Tuple<C, size-1> &) const;

};*/
/*
template<class C, int size, int pos>
class Tuple_InterpolatedMax : public Oper< Tuple<C, size>, double>{
public:
	virtual void operator()(Tuple<C, size> const & input, double &) const;

};


template<int size, int pos>
class Tuple_InterpolatedMaxS : public Oper< Tuple<double, size>, double>{
public:
	virtual void operator()(Tuple<double, size> const & input,double&) const;
};*/
/*
template<class C,int size>
class FTTuple : public Tuple<C,size>{
	public:


	template<class D> operator Tuple<D, size>() const; // Inverse Fourrier Transform!

};*/

/*
$r^2 = 4c^2d^2e^2+a^2c^4+2abc^2d^2+4ac^3de+4bcd^3e+b^2d^4+
y(-4abc^2d-4ac^3e-8c^2de^2-12bcd^2e-4b^2d^3)
x(-4a^2c^3-12ac^2de-8cd^2e^2-4abcd^2-4bd^3e)$\\
xy^2(-4abc-8ce^2-12bde)
x^2y(-4abd-12ace-8de^2)

y^2(+2abc^2+4c^2e^2+12bcde+6b^2d^2)
x^2(6a^2c^2+2abd^2+12acde+4d^2e^2)

x^3(-4a^2c-4ade)+y^3(-4b^2d-4bce)+xy^3(4be)+x^3y(4ae)
x^4a^2+y^4b^2+x^2y^2(+2ab+4e^2)
xy(+12ac^2e+12bd^2e+8abcd+16cde^2)
*/

// min sum (c + x.*b + xTAx )^2

// sum (2*(c + x.*b + xTAx ) * d/dx (c + x.*b + xTAx));

// sum (c + x.*b + xTAx) = 0;
// c = -1/n * sum (x.*b + xTAx);

// sum (c + x.*b + xTAx ) * x = 0;
// sum (x.*b + xTAx ) * x = -1/n * (sum (x.*b + xTAx)) * (sum x) ;


// min sum (x^TAx )^2

// 0 = sum (x^TAx ) * x x^T


// ../SharedLibraries/Advanced.h.gch: ../SharedLibraries/Advanced.hpp ../SharedLibraries/Advanced.h ../SharedLibraries/Hierachical.hpp primitive.o
//	g++ ${COMPILER_FLAGS}  -c ../SharedLibraries/Advanced.h -o ../SharedLibraries/Advanced.h.gch

class GaussianRegression{
    public:
    Tuple<Tuple<double, 0u>, 0u> means;
    Tuple<Tuple<double, 0u>, 0u> stds;


    void EMinit(unsigned int size, unsigned int nbstates);
    void EMregist(const Tuple<double, 0u>&);
    void EMfinit();

};

// Polymial in SIZE variables, where the maximum sum of degree is ORDER-1
template <class C, unsigned int SIZE, unsigned int ORDER=0u>
class PolyTuple{
public:
    Tuple<C, TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans> data;

    PolyTuple(){}

    PolyTuple(Tuple<C,SIZE> const & input, double weight = 1.0f);


    C operator[](const Tuple<C, SIZE>& in)const; // TODO

    C& cell(const Tuple<unsigned int, SIZE>& in); // TODO
    C cell(const Tuple<unsigned int, SIZE>& in) const; // TODO


    PolyTuple<C,SIZE,ORDER>& operator+=(PolyTuple<C,SIZE,ORDER> const & other){for(unsigned int i = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i!=0xFFFFFFFF;i--) data[i] += other.data[i];  return(*this);}
    PolyTuple<C,SIZE,ORDER>& operator-=(PolyTuple<C,SIZE,ORDER> const & other){for(unsigned int i = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i!=0xFFFFFFFF;i--) data[i] -= other.data[i];  return(*this);}


    PolyTuple<C,SIZE,ORDER>& shift(Tuple<C, SIZE> const & other){for(unsigned int k=0;k<SIZE;k++) shiftdir(other[k],k); return (*this);}
    PolyTuple<C,SIZE,ORDER>& shiftdir(C const & other, unsigned int d);
    PolyTuple<C,SIZE,ORDER>& statshift(Tuple<C, SIZE> const & other){for(unsigned int k=0;k<SIZE;k++) statshiftdir(other[k],k); return (*this);}
    PolyTuple<C,SIZE,ORDER>& statshiftdir(C const & other, unsigned int d);


    PolyTuple<C,SIZE,ORDER>& toZero(){for(unsigned int i = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i!=0xFFFFFFFF;i--) ExOp::toZero(data[i]); return(*this);}
    PolyTuple<C,SIZE,ORDER>& toOne(){for(unsigned int i = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i!=0;i--) ExOp::toZero(data[i]); ExOp::toOne(data[0]);return(*this);}

    PolyTuple<C,SIZE,ORDER*2-1> operator*(PolyTuple<C,SIZE,ORDER> const & other)const; //TODO

	PolyTuple<C,SIZE,ORDER> operator*(DBL_Weight const & other)const;
	PolyTuple<C,SIZE,ORDER>& operator*=(DBL_Weight const & other);


    PolyTuple<C,SIZE,ORDER-1> getDerivative(unsigned int direction)const; // DONE!

	template<unsigned int OSIZE> PolyTuple<C,OSIZE,ORDER> operator<<(TMatrix<C,OSIZE,SIZE> const & xform) const; // TODO


	class KeyIterator : public AbstractKeyIterator< Tuple<unsigned int, SIZE> , PolyTuple<C, SIZE, ORDER> >{
	public:
       unsigned int s;
       KeyIterator(const PolyTuple<C, SIZE, ORDER>& tar) :AbstractKeyIterator< Tuple<unsigned int, SIZE> , PolyTuple<C, SIZE, ORDER> >(tar){}
       bool first(); bool prev(); bool next(); bool last();
       bool first_oforder(unsigned int order);
    };

    void show(FILE* f = stdout, int level=0) const;
};
template <class C, unsigned int SIZE =0u>
class PolyThing{
	public:
	C data[SIZE];
	PolyThing<C,SIZE>& operator+=(const C& scalar){data[0] += scalar; return(*this);}
	PolyThing<C,SIZE>& operator+=(const PolyThing<C,SIZE>& other) {for(unsigned int i=0;i<SIZE;i++) data[i] += other.data[i];return(*this);}
	C operator[](double value) const;

	typename MT_IFTYPE<SIZE-2, PolyThing<C,SIZE-1> ,C>::TYPE mkDerivative() const;
};
template <class C>
class PolyThing<C, 0u>{
	public:
		int order;
		C* coef;
		PolyThing():coef(NULL){}
		PolyThing(PolyThing<C> const & clonefrom): order(clonefrom.order){
			coef = new C[order];  LFH_NICE_ALLOCERROR(coef ,"")
			int i;
			for(i=0;i<order;i++) coef[i] = clonefrom.coef[i];
			};
		~PolyThing(){if (coef) delete[](coef);}
		PolyThing<C>& operator=(PolyThing<C> const & clonefrom){
			if ((order != clonefrom.order)||(coef == NULL)){
				if (coef) delete[](coef);
				order = clonefrom.order;
				coef = new C[order]; LFH_NICE_ALLOCERROR(coef ,"")
				}
				int i;
				for(i=0;i<order;i++) coef[i] = clonefrom.coef[i];
			}
		void clear(){if (coef) memset(coef,'\0',sizeof(C)*(order));}
		void setOrder(int norder){
			if (coef) delete[](coef);
			coef = new C[norder+1]; LFH_NICE_ALLOCERROR(coef ,"")
			order = norder+1;
		}

		//	void PolyThing::ODEextends(LFHPrimitive::ODEfunctions ode, double factor, double coef);

		// Trivial scalar operations:

		PolyThing<C>& operator+=(int const & other){coef[0] += other;return(*this);}
		PolyThing<C>& operator+=(double const & other){coef[0] += other;return(*this);}

		template<class D> PolyThing<C>& operator+=(D const & other){coef[0] += other;return(*this);}

		PolyThing<C> secondDerivative();

		PolyThing<C>& operator-=(int const & other){coef[0] -= other;return(*this);}
		PolyThing<C>& operator-=(double const & other){coef[0] -= other;return(*this);}
		PolyThing<C>& operator*=(int const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] *= other;} return(*this);}
		PolyThing<C>& operator*=(double const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] *= other;} return(*this);}
		PolyThing<C>& operator/=(int const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] /= other;} return(*this);}
		PolyThing<C>& operator/=(double const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] /= other;} return(*this);}


		template<class D> PolyThing<C> operator+(D const & other) const {return( (PolyThing<C>(*this)) += other ); }
		template<class D> PolyThing<C> operator-(D const & other) const {return( (PolyThing<C>(*this)) -= other ); }
		template<class D> PolyThing<C> operator*(D const & other) const {return( (PolyThing<C>(*this)) *= other ); }
		template<class D> PolyThing<C> operator/(D const & other) const {return( (PolyThing<C>(*this)) /= other ); }

		// Simple poly operations:

		PolyThing<C>& operator+=(PolyThing<C> const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] += other.coef[i];} return(*this);}
		PolyThing<C>& operator-=(PolyThing<C> const & other) {int i; for(i= order-1;i>=0;i--) {coef[i] -= other.coef[i];} return(*this);}
		PolyThing<C> operator+(PolyThing<C> const & other) const {return( (PolyThing<C>(*this)) += other ); }
		PolyThing<C> operator-(PolyThing<C> const & other) const {return( (PolyThing<C>(*this)) -= other ); }


		double findMax(double minv =0.0f, double maxv =0.0f);
		Vector<mycomplex> getZeros();

		// complex poly operations:
		PolyThing<C> operator*(PolyThing<C> const & other) const{

		}
		C eval(double where);


		double newtonRoot(double guess, double min, double max) const;



		static PolyThing<C> interpolate(C* points, int nbpts);
		static PolyThing<C> interpolate_1free(C* points, int nbpts);
		static PolyThing<C> interpolate_2free(C* points, int nbpts);
		static PolyThing<C> interpolate_1free_middlefit(C* points, int nbpts); // error is totelated but for middle pts
		static PolyThing<C> interpolate_2free_middlefit(C* points, int nbpts); // error is totelated but for middle pts

	//	static PolyThing<C> interpolate_2free_middlefit(C* points, int nbpts); // error is totelated but for middle pts

	};
template < > double PolyThing<double>::findMax(double minv, double maxv);
template < > Vector<mycomplex> PolyThing<double>::getZeros();
template<class O, class I, class M, int _out_dim, int _in_dim, int _inner_dim>
class Chained_Function : ContinuousFunction<O,I,_out_dim,_in_dim>{
public:
	ContinuousFunction<M, I, _inner_dim, _in_dim>* first;
	ContinuousFunction<O, M, _out_dim, _inner_dim>* last;

	void operator()(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &) const;
//	void derivative(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	void derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;

};
template<class O, class I, int _out_dim, int _in_dim, int _parameter_dim>
class Function_Sum : public Oper2< Tuple<O, _out_dim> , Tuple<I, _in_dim > >{
	public:
	Vector< Tuple<I, _parameter_dim > > params;
	ContinuousFunction<O, I, _out_dim, _in_dim + _parameter_dim>* function;
	void operator()(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &) const;

	void derivative(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct) const;
	void derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &, Tuple<I, _in_dim > &) const;
};
template<int dims> // LLIkely <- [RVAR, RADIUS, C_0, C_1 , C_2...]
class SphereDistanceLikelyhood : public ContinuousFunction< double, double , 1, dims+2 >{
public:


	void operator()(Tuple<double, 1> &, Tuple<double, dims+2 > &) const;
	void derivative(Tuple<double, 1> &, Tuple<double, dims+2 > &, int in_direct) const;

	void doublederivative(Tuple<double, 1> &, Tuple<double, dims+2 > &, int in_direct, int in_direct2) const;
};
template<int dims> // LLIkely <- [RVAR, RADIUS, C_0, C_1 , C_2...]
class SphereDistanceLikelyhood_Derivative : public ContinuousFunction< double , double , dims+1 , dims+2 >{
public:

	void operator()(Tuple<double, dims+1> &, Tuple<double, dims+2 > &) const;
	void derivative(Tuple<double, dims+1> &, Tuple<double, dims+2 > &, int in_direct) const;
	void doublederivative(Tuple<double, dims+1> &, Tuple<double, dims+2 > &, int in_direct, int in_direct2) const;
};




//	v= (x/2 - (81.169*(x-255.5)*(asin(sqrt((x-255.5)*(x-255.5)+(y-255.5)*(y-255.5))/255.5)/sqrt((x-255.5)*(x-255.5)+(y-255.5)*(y-255.5)))))  * 4.75 - 478.125

//	v= sqrt(255.5*255.5 - (x-255.5)*(x-255.5)-(y-255.5)*(y-255.5))

// underbound on range [-1,1] for arbitrary fnctions
// the value matches on point -1,1 and the derivative

template<class C>
class QuadraticBound_old{
    C st, en;
    double err;
public:
    QuadraticBound_old(const C& start, const C& endpt, double mid_max_err): st(start), en(endpt), err(mid_max_err){}
    QuadraticBound_old<C>& operator+=(const QuadraticBound_old<C>& other){st += other.st;en += other.en; err += other.err;return(*this);}
    QuadraticBound_old<C>& operator-=(const QuadraticBound_old<C>& other){st -= other.st;en -= other.en; err += other.err;return(*this);}
    template<class O> QuadraticBound_old<C>& operator+=(const O& other){st += other;en += other;return(*this);}
    template<class O> QuadraticBound_old<C>& operator-=(const O& other){st += other;en += other;return(*this);}
    QuadraticBound_old<C>& operator*=(const double& other){st *= other;en *= other;err *= fabs(other);return(*this);}
    QuadraticBound_old<C>& operator/=(const double& other){st /= other;en /= other;err /= fabs(other);return(*this);}
    QuadraticBound_old<C>& operator*=(const QuadraticBound_old<C>& other){
        double n1 = ExOp::norm(st) * other.err + ExOp::norm(other.st) * err;
        double n2 = ExOp::norm(en) * other.err + ExOp::norm(other.en) * err;
        double tmperr = (n1 > n2) ? n1 : n2;
        double errp = err * other.err;
        double s1 = ExOp::norm(st - en) / errp;

        if ((ExOp::isValid(s1))&&(fabs(s1)< 1.0f)){
            // there are 2 alterative optima, at s1 and -s1
            //n1 =
        }
        s1 = ExOp::norm(other.st - other.en)/ errp;
        if ((ExOp::isValid(s1))&&(fabs(s1)< 1.0f)){
            // there are 2 alterative optima, at s1 and -s1

        }
        err = n1;

    st*= other.st;en *= other.en; return(*this);}

	double minRoot() const; // returns minimum time where ||F[x]|| = 0 is possible, if impossible, 1 is returned (since nothing is known >1)
};

template< >	double QuadraticBound_old<double>::minRoot() const;

template<class C, unsigned int order>
class PolyBound{
    C data[order];
    double err;
public:
    PolyBound(double mid_max_err =0.0f): err(mid_max_err){}
    PolyBound<C,order>& operator+=(const PolyBound<C,order>& other){for(unsigned int i=0;i<order;i++) data[i] += other.st; err += other.err;return(*this);}
    PolyBound<C,order>& operator-=(const PolyBound<C,order>& other){for(unsigned int i=0;i<order;i++) data[i] -= other.st; err += other.err;return(*this);}
    template<class O> PolyBound<C,order>& operator+=(const O& other){data[0] += other;return(*this);}
    template<class O> PolyBound<C,order>& operator-=(const O& other){data[0] -= other;return(*this);}
    PolyBound<C,order>& operator*=(const double& other){for(unsigned int i=0;i<order;i++) data[i] *= other;err *= fabs(other);return(*this);}
    PolyBound<C,order>& operator/=(const double& other){for(unsigned int i=0;i<order;i++) data[i] /= other;err /= fabs(other);return(*this);}
    PolyBound<C,order>& operator*=(const PolyBound<C,order>& other){PolyBound<C,order> tmp = (*this) * other; return ((*this) = tmp); }
    PolyBound<C,order> operator*(const PolyBound<C,order>& other) const{ PolyBound<C,order> fout;
        unsigned int i,j;
        fout.err = err * other.err;
        for(i=0;i<order;i++) {
            fout.data[i] = data[i] * other.data[0];
            fout.err += err * ExOp::norm(other.data[i]);
            fout.err += other.err * ExOp::norm(data[i]);
        }
        for(j=1;j<order;i++) for(i=0;i+j<order;i++) fout.data[j+i] += data[i] * other.data[j];
        // ez part done,


        return(fout);
    }


	double minRoot() const; // returns minimum time where ||F[x]|| = 0 is possible, if impossible, 1 is returned (since nothing is known >1)
};

// pointer to char, in UTF_8, does NOT own data
class UTF_8{
    public:
    const char* data;
    unsigned int strlen() const;
    unsigned int operator[](int pos) const;
};

// underbound on range [0,1] for arbitrary fnctions
// the value matches on point 0 and its first 3 derivatives
template<class C>
class QuarticBound{
public:
    C data[4];
    double err; // P(x) = a + bx + cx2 +dx3 +- err*x*x*x*x
    QuarticBound<C>& operator+=(const QuarticBound<C>& other){data[0] += other.data[0];data[1] += other.data[1];data[2] += other.data[2];data[3] += other.data[3]; err += other.err;return(*this);}
    QuarticBound<C>& operator-=(const QuarticBound<C>& other){data[0] -= other.data[0];data[1] -= other.data[1];data[2] -= other.data[2];data[3] -= other.data[3]; err += other.err;return(*this);}
    QuarticBound<C>& operator+=(const double& other){data[0] += other;return(*this);}
    QuarticBound<C>& operator-=(const double& other){data[0] -= other;return(*this);}
    QuarticBound<C>& operator*=(const double& other){data[0] *= other;data[1] *= other;data[2] *= other;data[3] *= other; err *= fabs(other);return(*this);}
    QuarticBound<C>& operator/=(const double& other){data[0] /= other;data[1] /= other;data[2] /= other;data[3] /= other; err /= fabs(other);return(*this);}

	double minRoot() const;
};

template< >	double QuarticBound<double>::minRoot() const;



	class FileIO{
		public:
		static void loadBMP3(const char * path,DataGrid<unsigned char,3> &im);
		static void saveBMP3(const char* const path,const DataGrid<unsigned char,3> &im);
		};


extern double timestep;
extern myHashmap< unsigned int , pair< void*, unsigned int> > AliasBank;
extern myHashmap< const void* , unsigned int > AliasOf;

	class LinkEntry;
	template <class Type> class Ptr;
	class LinkRegister;

	template <class Host, class Target, unsigned int Hoff, unsigned int T_off> class DLink;

	class OwnerRegisterEntry;
	class OwnerRegisterCmp;


class ratio{
public:
	int num;
	unsigned int den;

	ratio(){}
	ratio(int n, unsigned int d): num(n),den(d){}
	ratio(double val);
	operator double ()const{return(((double)num)/den);};


};

template<class C> class AbsTuple{
public:
	virtual int getSize() const =0; // return size
	virtual C& operator[](int const pos)=0;
	virtual C operator[](int const pos) const=0;
	virtual operator C* ()=0;

};
class GaussianProcess{
public:
	template<class C> static void covMatrix(DataGrid<double, 2> &f_out, Vector<C>& data, double (*metric)(const C&, const C&), double noiselevel =0.0f);
	template<class C> static void covMatrix(DataGrid<double, 2> &f_out, Vector<C>& data, Vector<C>& query, double (*metric)(const C&, const C&));
	template<class C> static void covMatrix(DataGrid<double, 2> &f_out, Vector<WeightElem<C,1> >& data, double (*metric)(const C&, const C&), double noiselevel =0.0f);
	template<class C> static void covMatrix(DataGrid<double, 2> &f_out, Vector<WeightElem<C,1> >& data, Vector<C>& query, double (*metric)(const C&, const C&));

	template<class K, class C, int NBDIM> static void GPweightsregression(DataGrid<C, NBDIM> &f_out, const DataGrid<C, NBDIM> &f_data, const Vector<K> &f_Keys , double (*dist)(const K &a, const K &b), double (*weight)(const K &a));
	template<class K, class C> static Vector<C> GPweightsregression(const Vector< KeyElem<K,C> > &f_data, const Vector<K> &f_Keys);
	};


class GaussianProcessMK2{
public:
	Tuple<double, 4u> data;
	Trianglix<double> scale;
	double& getMean(){return data[0];} double getMean()const{return data[0];}
	double& getNoise(){return data[1];} double getNoise()const{return data[1];}
	double& getSignal(){return data[2];} double getSignal()const{return data[2];}
	double& getScale(){return data[3];} double getScale()const{return data[3];}
	template<int ISCONST> class Iterator{
    public:
        typename MetaType<GaussianProcessMK2,ISCONST>::IS_CONST_REF target;
        typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
        uint32_t cur;
        Iterator(typename MetaType<GaussianProcessMK2,ISCONST>::IS_CONST_REF trg): target(trg){}
        operator bool (){cur = 0; return true;}
        bool operator++(int){return (++cur < target.scale.totsize()+3);}
        uint32_t operator()()const{return cur;}
        typename MetaType<double*,ISCONST>::IS_CONST operator->(){return  (cur < 3) ? & target.data[cur] : &target.scale.data[cur-3];}
        typename MetaType<double,ISCONST>::IS_CONST_RETURN_RVAL operator*(){return (cur < 3) ? target.data[cur] : target.scale.data[cur-3];}
        Iterator& mkIterator(){return *this;}
	};


	class LearnParaTask{
		/* Iterates over Gaussian Processes, which share identical time/state associated with a fit number of datapoint
		*/
		double LL;
	public:
		LearnParaTask():LL(0.0){}
		void operator<<=(LearnParaTask& other){LL += other.LL;}
		double operator()(){return LL;}
		template<class P, class S, class O> void operator()(P& gp_producer, S& state_iterator, O& observation_dico);
	};

	class CmpStateDerivative{
		/* Iterates over Gaussian Processes, which share identical time/state associated with a fit number of datapoint
		*/
		double LL;
	public:
		CmpStateDerivative():LL(0.0){}
		void operator<<=(CmpStateDerivative& other){LL += other.LL;}
		double operator()(){return LL;}
		template<class P, class S, class G> void operator()(P& observation_producer, S& state_iterator, G& gp_dico);
	};

	Iterator<0> mkParamIterator(){return Iterator<0>(*this);}
	Iterator<1> mkParamIterator()const{return Iterator<1>(*this);}

	mutable Trianglix<double, 0u> cached_kernel;
	mutable Tuple<double> cached_projection;

	GaussianProcessMK2& init(double mean = 0.0, double signal = 1.0, double noise = 0.125, double window= 1.0, uint32_t nbdims=2) {scale.setSize(nbdims).toOne() *= window; this->getMean() = mean; this->getNoise() = noise; this->getSignal() = signal; return *this;}
	GaussianProcessMK2& toZero(){data.toZero(); return *this;}

	template<class I, ENABLEIF_ARRAY(I, double) =true> GaussianProcessMK2& setKernel(I timestamp, I noisescale, uint32_t length=0);
	template<class I, class J, ENABLEIF_ARRAY(I, double) =true, ENABLEIF_ARRAY(J, double) =true> GaussianProcessMK2& setKernel2(I timestamp, J noisescale, uint32_t length=0);

	Tuple<double> setKernel3(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp) const;
	//template<class I, ENABLEIF_WRITE_ITERATOR_ANYKEY(I) =true> GaussianProcessMK2& setKernelWithScaledNoise(I timestamp, I noisescale, uint32_t length=0);

	template<class I, ENABLEIF_ARRAY(I, double) =true> double wrLLDerivative(Tuple<double, 4u>& fout, I observ, I timestamp, I noisescale, bool verbose = false) const;

	double gradientAscent(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, unsigned int nbsteps=1000, double *startLL = NULL);


	double gradientAscentMK2(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, unsigned int nbsteps=1000);

	double getLL(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp);

	double getMeanAndDeriv_p1(double extrastate, double weight, Tuple<double>& statederiv, RemoteMemory<const double,2u> delta, double* buffer)const;

	template<class I, typename RequireType<const I>::beIterator = true> void syntaxtest(I testtest){}
	template<class I, typename RequireType<I>::beIterator = true> void syntaxtest2(I testtest){}

	template<class I, class J, class K, ENABLEIF_ARRAY(I, double) =true, ENABLEIF_ARRAY(J, double) =true, ENABLEIF_ARRAY(K, double) =true>
		double wrLLDerivative2(Tuple<double>& fout, I observ, J timestamp, K noisescale, bool verbose = false) const;

	template<class I, class J, ENABLEIF_ARRAY(I, (Tuple<double, 2u>)) =true, ENABLEIF_ARRAY(J, double) =true>
		double wrLLDerivative3(Tuple<double>& fout, I observ, J timestamp, bool verbose = false) const; // obs are mean and vars

	template<class I, class J, class K, ENABLEIF_ARRAY(I, double) =true, ENABLEIF_ARRAY(J, double) =true, ENABLEIF_ARRAY(K, double) =true>
		double wrLLTimeDerivative(Tuple<double>& fout, I observ, J timestamp, K noisescale, bool verbose = false) const;

	double addLLTimeDerivative(Tuple<double>& fout, RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, bool verbose) const;

	//double addLLTimeDerivative2(Tuple<double>& fout, double value, RemoteMemory<const double,1u> &state, RemoteMemory<const double,2u> timestamp, bool verbose) const;

	void checkDomain(CurvatureSearchScope &learnscope);
	void checkDomain(AdaptiveLearningScope &learnscope);

	DataGrid<double, 3u> makeImage(RemoteMemory<const double,2u> observ, RemoteMemory<const double,2u> timestamp, int sizex = 256, int sizey =256, int borderx = 0, int bordery = 0, Tuple<double, 0u> *extraout = NULL) const;
	const GaussianProcessMK2& show(FILE* f = stdout, int level=0)const{fprintf(f,"GP with %e +- %e and signal %e over %i dims\n", getMean(), getNoise(), getSignal(), scale.getSize()); return *this;}
};

class GaussianProcessArray{
public:
	GaussianProcessArray& initUsingPCA(RemoteMemory<const double,2u> observ, uint32_t nbPCs);
	GaussianProcessArray& initUsingPCA(RemoteMemory<const double,3u> observ, uint32_t nbPCs){if (observ.dims[0] == 2) this->initUsingPCA(observ.selectSlice(0,0),nbPCs); return *this;}

	Tuple<GaussianProcessMK2> gps;
	TMatrix<double> hstate;

	template<int ISCONST> class Iterator{
    public:
        typename MetaType<GaussianProcessArray,ISCONST>::IS_CONST_REF target;
        typedef typename MetaType<uint32_t ,ISCONST>::IS_CONST KEYITERATOR_TYPE;
        uint32_t cur;
        uint32_t gpparamsize;
        Iterator(typename MetaType<GaussianProcessArray,ISCONST>::IS_CONST_REF trg): target(trg){gpparamsize = 3 + ((target.hstate.sizes[0] + (target.hstate.sizes[0] +1))>>1);}
        operator bool (){cur = 0; return target.gps.getSize() > 0;}
        bool operator++(int){return (++cur < target.gps.getSize() * gpparamsize + target.hstate.sizes[0] * target.hstate.sizes[1]);}
        uint32_t operator()()const{return cur;}
        typename MetaType<double*,ISCONST>::IS_CONST operator->(){
        	if (cur >= target.gps.getSize() * gpparamsize) return &(target.hstate.data[cur - target.gps.getSize() * gpparamsize]);
        	int res = cur % gpparamsize;
        	return (res < 3) ? & (target.gps[cur / gpparamsize].data[res]) : & (target.gps[cur / gpparamsize].data[res-3]);
		}
        typename MetaType<double,ISCONST>::IS_CONST_RETURN_RVAL operator*(){
        	if (cur >= target.gps.getSize() * gpparamsize) return target.hstate.data[cur - target.gps.getSize() * gpparamsize];
        	int res = cur % gpparamsize;
        	return (res < 3) ? target.gps[cur / gpparamsize].data[res] : target.gps[cur / gpparamsize].data[res-3];}
        Iterator& mkIterator(){return *this;}
	};

	Iterator<0> mkIterator(){return Iterator<0>(*this);}

	uint32_t getNbvalues()const{return gps.getSize();}

	GaussianProcessMK2& operator[](int i) {return gps[i];}
	const GaussianProcessMK2& operator[](int i) const {return gps[i];}
	double cacheObservations(RemoteMemory<const double,3u> observ);
	double learnThis(RemoteMemory<const double,3u> observ, uint32_t nbsteps = 10);


	TMatrix<double> learnPCAreduction(RemoteMemory<const double,4u> observ, uint32_t nb_latent = 2, uint32_t nbsteps = 10, uint32_t nbpc_per_ct=0);  // (mean,std) x (gp_obs) x (gene) x (celltype)
	double learnGPParams(RemoteMemory<const double,3u> observ); // requir3e set hidden state

	Tuple< TMatrix<double>, 3u > localPoissonRegression(RemoteMemory<const uint32_t,2u> observ) const;
	GaussianProcessArray& setSize(uint32_t nbvalues, uint32_t nbdims=2){gps.setSize(nbvalues); for(unsigned int i=0;i<nbvalues;i++) gps[i].init(); return *this;}

	TMatrix<double> makeGParammatrix()const;
};

	extern LinkRegister LinkMem;

template<class H, class R, int Hoffoff, int Roffoff>
class dxorptr{
public:
	LFH_address resprev;
	LFH_address hosnext;

	dxorptr();
	dxorptr(const LFH_address rp); // used for list init
	dxorptr(const LFH_address rp, const LFH_address hn);
	const H* clear();
	operator const H*() const;
	const H* operator->() const;
	dxorptr<H,R,Hoffoff,Roffoff>& operator=(const H* target);

};

// an array which assumes that another slave exists in the next memory slot!
template<class T>
class SlaveArray{
    public:
    T* data;
    T& operator[](unsigned int i){return data[i];}
    unsigned int getSize() const{ return ((&data)[1] - data) / sizeof(T);  }
};

template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
class dxorptrarray{
	public:
		dxorptr<H,R,Hoffoff,Roffoff> aptr[nbptr];

		const H* operator[](unsigned int) const;

		dxorptrarray();
		operator const dxorptr<H,R,Hoffoff,Roffoff>& () const;

//		operator const H*() const;
//		const H* operator->() const;
//		dxorptr<H,R,Hoffoff,Roffoff>& operator=(const H* target);

	};
// Semaphore semaphore;
// syntax if (auto access = semaphore.tryRead/Write()){
//    success
// }else{
//    failed
// }
//
class Semaphore{
    atomic<int> sema;
    friend class ReadAccess;
    friend class WriteAccess;
    friend class WaitForReadAccess;
    friend class WaitForWriteAccess;
    public:
    class ReadAccess{
        Semaphore& target;
        int state;
        public:
        ReadAccess(Semaphore& _target):target(_target){}
        ~ReadAccess(){target.sema.fetch_add(-state);}
        operator bool ();
    };
    class WriteAccess{
        Semaphore& target;
        int state;
        public:
        WriteAccess(Semaphore& _target):target(_target){}
        ~WriteAccess(){target.sema.fetch_add(-state);}
        operator bool ();
    };
    class WaitForReadAccess{
        Semaphore& target;
        int delay;
        public:
        WaitForReadAccess(Semaphore& _target, int waitdelay):target(_target),delay(waitdelay){}
        ~WaitForReadAccess(){target.sema.fetch_add(-1);}
        operator bool ();
    };
    class WaitForWriteAccess{
        Semaphore& target;
        int delay;
        public:
        WaitForWriteAccess(Semaphore& _target, int waitdelay):target(_target),delay(waitdelay){}
        ~WaitForWriteAccess(){target.sema.fetch_add(-256);}
        operator bool ();
    };

    Semaphore():sema(0){}
    ReadAccess tryRead(){return ReadAccess(*this);}
    WriteAccess tryWrite(){return WriteAccess(*this);}
    WaitForReadAccess waitRead(int waitdelay = 1000){return WaitForReadAccess(*this,waitdelay);}
    WaitForWriteAccess waitWrite(int waitdelay = 1000){return WaitForWriteAccess(*this,waitdelay);}
};



class MutexAccess{
    std::mutex mtx;
    std::atomic<int> sema;
    std::condition_variable cv;
    friend class ReadAccess;
    friend class WriteAccess;
    public:
    class ReadAccess{
        MutexAccess& target;
        int state;
        public:
        ReadAccess(MutexAccess& _target):target(_target){}
        ~ReadAccess(){}
        operator bool ();
    };
    class WriteAccess{
        std::unique_lock<std::mutex> lock;
        int state;
        public:
        WriteAccess(MutexAccess& _target);
        ~WriteAccess(){}

    };
    MutexAccess(){}
    ReadAccess tryRead(){return ReadAccess(*this);}
//    const WriteAccess& tryWrite(){return WriteAccess(*this);}
};

template<class D, class T>
class ParamQueue{
    HeapTree<T> heap;
    myHashmap<T,D> hmap;
    T curtime;
    MutexAccess mutx;
    KeyElem<D,T> asyncbuffer[256];
public:
    void insert_sync(const T& t, const D& d){heap.addEntry(t); hmap.addEntry(t,d);}
    void remove_sync(const T& t, const D& d){hmap.removeEntry(t,d);}
    void insert_async(const T& t, const D& d);
    void remove_async(const T& t, const D& d);
    bool get(const T& t, D& d);
    uint32_t getSize()const{return hmap.getSize();}
};

template<class H, class R, int Hoffoff, int Roffoff>
class dxorlist : public dxorptr<H,R,Hoffoff,Roffoff>{
public:
	dxorlist();
	dxorptr<H,R,Hoffoff,Roffoff>& operator=(const H* target);
	unsigned int getSize();
	operator const Vector<R*> () const;
};

template<class H, class R, int Hoffoff, int Roffoff>
class xptr{
public:
	H* target;
	xptr<H,R,Hoffoff,Roffoff>* next;
	LFH_address resprev;
	xptr();
	xptr(H* _target, xptr<H,R,Hoffoff,Roffoff>* _next, LFH_address _resprev): target(_target),next(_next),resprev(_resprev) {}
	~xptr(){if (target != NULL) exit(123123);}
	void clear();
	H* operator->() const;
	xptr<H,R,Hoffoff,Roffoff>& operator=(H const * const);
};

// for ez-syntax, must implement all ptr functions!
template<class Host, class Resident, int Hoff, int Roff>
class xptrref{
	public:
	Resident* res;
	xptr<Host,Resident,Hoff,Roff>& trg;
	xptrref(xptr<Host,Resident,Hoff,Roff>& _trg, Resident* _res):trg(_trg), res(_res) {}

		void clear();
		xptrref<Host,Resident,Hoff,Roff>& operator=(Host * value);
		Host * operator->();

		operator Host* ();
	};



class GraphAxisData{
public:
	char* axis_name;
	double range[2];
	unsigned int nbbins;// nb+labels on scales, or bins
	GraphAxisData(): axis_name(NULL) {}
	~GraphAxisData(){delete[](axis_name);}
	GraphAxisData(const GraphAxisData &other) {axis_name = cloneString(other.axis_name);}
	GraphAxisData& operator=(const GraphAxisData &other) {delete[](axis_name);axis_name = cloneString(other.axis_name);return(*this);}
};

class Graph{ // 2D representations!
public:
	char* title;
	Vector<GraphAxisData> axises;
	unsigned int extra_rect[4];
	unsigned int graph_size[2];
	Graph(): title(NULL) {}
	~Graph(){delete[](title);}
	Graph(const Graph &other) {title = cloneString(other.title);}
	Graph& operator=(const Graph &other) {delete[](title);title = cloneString(other.title);return(*this);}

	template<class O_T> const DataGrid<O_T,3>& render_Histogram(DataGrid< O_T, 3 > &f_out, const DataGrid<double, 1 > &image) const; // 1D histogram
	template<class O_T> const DataGrid<O_T,3>& render_Histogram(DataGrid< O_T, 3 > &f_out, const DataGrid<double, 2 > &image) const; // 1D histogram, multiple channels
	void setInnerSize(unsigned int width,unsigned int height){graph_size[0] =width; graph_size[1] =height; }
};

template <class C>
class IntervalSet{
public:
    virtual ~IntervalSet(){}
	virtual const SetComparison& compare(const C & other)const =0;
	virtual const SetComparison& compareInterval(const C & min, const C & max)const =0;
	template<class D> const SetComparison& compare(const KeyElem<C, D> & other)const  {return(compare(other.k));}
	template<class D> const SetComparison& compareInterval(const KeyElem<C, D> & min, const KeyElem<C, D> & max)const {return(compareInterval(min.k,max.k));}
};

template <class C>
class Interval : public IntervalSet<C>{
public:
	C min;
	C max;
	Interval(C _min, C _max);
	SetComparison compare(const C & other) const;
	SetComparison compareInterval(const C & min, const C & max)const;
};

template < > double PolyThing<double>::findMax(double minv, double maxv);
template < > Vector<mycomplex> PolyThing<double>::getZeros();

template<class C>
class Spline{
public:
	Vector< KeyElem<double, C> >	spldata;
	C* SplineExtra;
	int last;

	Spline();
	void setExtra();
	C operator()(const double& value);
	void addPoint(const double&  where, const C&  data);
};



template <class C,class I, unsigned int D, unsigned int L = 0u>
class HyperPaint{
public:
    typedef I IntType;
    typedef void MergeScopeType;
    static const int nb_dims = D;
    static const int nb_lead = L;

    C data;
    HyperPaint<C,I,D,L>& toZero();
    C operator()()const{return data;}

	bool canMerge() const{return false;}
    bool mergeInto(HyperPaint<C,I,D,L>& brother_io, const HyperCursor<I,D,L> &this_loc) const;
    bool canSimplify() const{ return false; }
    bool simplifyAt(HyperPaint<C,I,D,L>& fout, const HyperCursor<I,D,L> &location) const;

    bool operator==(const HyperPaint<C,I,D>& other) const{return data == other.data;}
    bool operator!=(const HyperPaint<C,I,D>& other) const{return data != other.data;}
    bool operator>(const HyperPaint<C,I,D>& other) const{return data > other.data;}
    bool operator<(const HyperPaint<C,I,D>& other) const{return data < other.data;}
    bool operator>=(const HyperPaint<C,I,D>& other) const{return data >= other.data;}
    bool operator<=(const HyperPaint<C,I,D>& other) const{return data <= other.data;}

    void show(FILE* f = stdout, int level =0)const;
    ERRCODE save(FILE* f) const;
	ERRCODE load(FILE* f);
};

// TrivialityMask 4 bits
// every channel can be independently trivial, trivial agrees with anything

// 1001
// 0001
// or
// 1001
// (a.d.intersection(a.k) ^ b.d.intersection(b.k)) & (a.d.nonTrivialMask() && b.d.nonTrivialMask() ) == 0 => merge possible, where_ever

// class holding a set of operations to evecute on a stack of "ambiguous" variables
// operations are defined as uin32t_t codes, some that are built-in, and some that needs to be resolved externally via a provided function

enum LFHSCRIPT_enum : uint16_t{
	LFHSCRIPT_IF_THENELSE = 0xFFFF,

	LFHSCRIPT_CHECK_EQUAL = 0xFFFE,
	LFHSCRIPT_CHECK_UNEQUAL = 0xFFFD,
	LFHSCRIPT_CHECK_LESSTHAN = 0xFFFC,
	LFHSCRIPT_CHECK_GRTRTHAN = 0xFFFB,
	LFHSCRIPT_CHECK_LESSEQTHAN = 0xFFFA,
	LFHSCRIPT_CHECK_GRTREQTHAN = 0xFFF9,

	LFHSCRIPT_64BITS_TO_F32BITS  = 0xFFF4,
	LFHSCRIPT_64BITS_TO_F64BITS  = 0xFFF5,
	LFHSCRIPT_F64BITS_TO_F32BITS = 0xFFF6,
	LFHSCRIPT_F32BITS_TO_F64BITS = 0xFFF7,
	LFHSCRIPT_LOAD_CONST_ZERO =    0xFFF0,
	LFHSCRIPT_LOAD_CONST_16BITS =  0xFFF1,
	LFHSCRIPT_LOAD_CONST_32BITS =  0xFFF2,
	LFHSCRIPT_LOAD_CONST_64BITS =  0xFFF3,
	LFHSCRIPT_ADD_8BITS = 0xFFE0,
	LFHSCRIPT_ADD_16BITS = 0xFFE1,
	LFHSCRIPT_ADD_32BITS = 0xFFE2,
	LFHSCRIPT_ADD_64BITS = 0xFFE3,
	LFHSCRIPT_LOGICAL_OR = 0xFFE4,
	LFHSCRIPT_OR = 0xFFE5,
	LFHSCRIPT_ADD_F32BITS = 0xFFE6,
	LFHSCRIPT_ADD_F64BITS = 0xFFE7,

	LFHSCRIPT_SUB_8BITS = 0xFFE8,
	LFHSCRIPT_SUB_16BITS = 0xFFE9,
	LFHSCRIPT_SUB_32BITS = 0xFFEA,
	LFHSCRIPT_SUB_64BITS = 0xFFEB,
	LFHSCRIPT_LOGICAL_XOR = 0xFFEC,
	LFHSCRIPT_XOR = 0xFFED,
	LFHSCRIPT_SUB_F32BITS = 0xFFEE,
	LFHSCRIPT_SUB_F64BITS = 0xFFEF,

	LFHSCRIPT_MULT_8BITS = 0xFFD0,
	LFHSCRIPT_MULT_16BITS = 0xFFD1,
	LFHSCRIPT_MULT_32BITS = 0xFFD2,
	LFHSCRIPT_MULT_64BITS = 0xFFD3,
	LFHSCRIPT_LOGICAL_AND = 0xFFD4,
	LFHSCRIPT_AND = 0xFFD5,
	LFHSCRIPT_MULT_F32BITS = 0xFFD6,
	LFHSCRIPT_MULT_F64BITS = 0xFFD7,

	LFHSCRIPT_DIVI_8BITS = 0xFFD8,
	LFHSCRIPT_DIVI_16BITS = 0xFFD9,
	LFHSCRIPT_DIVI_32BITS = 0xFFDA,
	LFHSCRIPT_DIVI_64BITS = 0xFFDB,
	LFHSCRIPT_LOGICAL_NOT = 0xFFDC,
	LFHSCRIPT_NOT = 0xFFDD,
	LFHSCRIPT_DIVI_F32BITS = 0xFFDE,
	LFHSCRIPT_DIVI_F64BITS = 0xFFDF,

	LFHSCRIPT_GOTO = 0xFFCF,
	LFHSCRIPT_GOTO_IFFALSE_TERNARY_FINAL = 0xFFCE, //
	LFHSCRIPT_GOTO_IFTRUE = 0xFFCD,
	LFHSCRIPT_GOTO_IFFALSE = 0xFFCC,
	LFHSCRIPT_GOTO_IFTRUE_TERNARY = 0xFFCB,
	LFHSCRIPT_GOTO_IFFALSE_TERNARY = 0xFFCA,
	LFHSCRIPT_GOTO_ = 0xFFC8, // range for gotos


	LFHSCRIPT_LAST_OPERATION = 0xFF00
};

class Script;


class ScriptScope{
public:
	AliasedHashmap<string, KeyElem<uint32_t, Ambiguous> > scope;
};

// Function type: return type arg type
class ScriptCompiler{
	int mkCommand_matchtypes_routine(StructuredHash<uint16_t, uint32_t>::TreeIterator<0> &ite, uint32_t dacodecode, uint32_t ssize);
	void errorNoArgumentMatch_routine(string fnc, const char* precise);
	void errorSet(const char * fnc);
	void errorSet(const char * fnc, char dachar);
	void errorSet(const char * fnc, const char *dastr);
public:
	bool printError(FILE* f = stdout);
	// type, 0x8000 bit is L-value, and
	AliasedHashmap<string, Tuple<uint32_t> >  fnccodes;
	Vector< std::function< char* (Ambiguous*&) > > externfnc;
	AliasedHashmap<string, KeyElem<uint32_t, Ambiguous> > scope; // L-value R-value
	char* errbuffer;
	AliasedHashmap<string, uint32_t> types; // data is storage type
	uint32_t nbdeclared;
	ScriptCompiler();
	Script mkCommand(const char* command);
	void initFnccodes();
	void clearError(){if (errbuffer) {delete[](errbuffer); errbuffer = NULL;} }
	uint16_t getTypeCode(const char* type) const{return (uint16_t) types.getAlias(string(type));}
	uint16_t getTypeCode(const string type) const{return (uint16_t) types.getAlias(type);}
	uint16_t addEnumType(string type);
	void addEnumType(uint16_t enumCode, const char* type, uint32_t value);

	void declareFunction(string name, Tuple<const char*> arg_types);
	void declareFunction(string name, Tuple<uint16_t> arg_types);
	//void addExternFunction(string name, Tuple<uint16_t> arg_types , uint32_t fnccode );
	void defineFunction(string name, std::function< char* (Ambiguous*&) > dafnc );
	//void addExternFunction(string name, Tuple<const char*> arg_types , uint32_t fnccode );

	void loadDeclarations(const Vector< KeyElem< KeyElem<string, uint32_t>, Tuple<uint32_t> > > & dacst);
	void saveDeclarations(Vector< KeyElem< KeyElem<string, uint32_t>, Tuple<uint32_t> > > & dacst) const;

	Ambiguous executeCode(const uint16_t* code, uint32_t length, ScriptScope* fncscopr = NULL);

	const ScriptCompiler& show(FILE *f =stdout, int level =0)const;
};




class Script{
public:
	Tuple<uint16_t> code;
	ScriptCompiler& funchost;
	Script(ScriptCompiler& _funchost): funchost(_funchost){}
	Ambiguous execute() const;
	//Script& toCommand(const char* command, const myHashmap<string, uint32_t> &fncodes);

	const Script& show(FILE* f = stdout, int level =1) const;
};

// TODO list hurray:

// integrate only carlos

// only carlos
// all epi

// bulk deconvolution

// deconvoluve lesions?


class ScriptNode{
    public:
    virtual void operator()(void* arglist)=0;
    virtual void nbargs(unsigned int &min, unsigned int &max) const=0;
    virtual void args_type(unsigned int offset, char* fout_buffer) const=0;
    virtual char nodetype() const=0;
    virtual ~ScriptNode(){}
};
class ScriptNode_Type : public ScriptNode{
    public:
    void operator()(void* arglist){};
    virtual void nbargs(unsigned int &min, unsigned int &max) const{};
    virtual void args_type(unsigned int offset, char* fout_buffer) const{};
    virtual char nodetype() const{return(1);};
};
class ScriptScopeOld{
    public:
    char *script;
    char *cursor;
    char *stack;
    char compbit;
    myHashmap<string, unsigned int> dico;

    void initdictionary();
    void finitdictionary();
    void compileScript(FILE* f);

    unsigned int operator()();
};
template<unsigned int DIMS>
class BBNetwork{
    class SampleScope{
        public:
        unsigned int routine_ite;
        HeapTree<KeyElem<double, Tuple< unsigned int, 2> > > proxi;
        Tuple<double> dist_buf;
        Tuple<unsigned int> dirty_buf;
        myHashmap<unsigned int, unsigned int> indexmap;
    };
public:
    Tuple< GaussElem< Tuple<double, DIMS> > > nodes;

    union{
    SampleScope* sample_scope;
    };

    // routine for 'k-mean' clustering for likelihood ratio metric
    void sampleInit(unsigned int nb_nodes);
    void sampleRegist(const GaussElem< Tuple<double, DIMS> > &inst);
    void sampleFinit();


	Tuple< KeyElem<double , unsigned int > > classmarginal(const GaussElem< Tuple<double, DIMS> > &inst, const double epsilon =0.0f) const;

    void EMinit();
    void EMregist(const GaussElem< Tuple<double, DIMS> > &inst);
    void EMfinit();
};
 // return Id of compiled function calls
union Instructions{
    void* p;
    unsigned int v;
};
class ScriptHead{
    public:
    myHashmap<string, ScriptNode*> parse_dictionary;
    ScriptNode_Type type;
    ScriptHead();
    Vector<Instructions> parseScript(const char*);

    void runScript( Vector<Instructions> &);
};


// Distribution over sparse vectors

template<unsigned int INPUTSIZE>
class NeuralWindow{
public:
	DataGrid<double, INPUTSIZE> weights;
	Tuple<unsigned int,0u> outputsize;

	void setOutputAndWindowSize(const Tuple<unsigned int,0u> &_outputsize,const Tuple<unsigned int, INPUTSIZE - 1> & _windowsize);
	template<class C, unsigned int IS, unsigned int OS> void wrForward(DataGrid<C, OS> &f_out, const DataGrid<C, IS> &f_in) const;

	static DataGrid<double, 3u> imageFeatureExtraction(const DataGrid<unsigned char, 3u>& image);
	static DataGrid<unsigned char, 3u> genImageFromFeatures(const DataGrid<double, 3u>& image);

	static DataGrid<double, 3u> imageHartWaveletTransform(const DataGrid<double, 3u>& image);
};



template<int MAG, int MODE, int HEADSIZE>
class MessageBuffer{
public:
	uint8_t buffer[(1 << MAG)];
	MessageBuffer();

	class NewMessageScope{
		NewMessageScope(MessageBuffer &_target):target(_target) {}
		friend class MessageBuffer;
		MessageBuffer &target;
		uint32_t nextpos;
	public:
		uint8_t* msg;
		~NewMessageScope();
		NewMessageScope(NewMessageScope&&o):msg(o.msg),target(o.target),nextpos(o.nextpos){o.msg=NULL;}
		NewMessageScope(const NewMessageScope&)=delete;
		NewMessageScope& operator=(const NewMessageScope&)=delete;
		NewMessageScope& operator=(NewMessageScope&&o){msg=o.msg;o.msg=NULL; nextpos = o.nextpos;return*this;}
		operator bool(){return msg != NULL;}
		uint8_t* operator+(uint32_t offset){return msg + offset + HEADSIZE;}
		uint8_t* operator()(){return msg;}
	};
	uint32_t getCurrentPosition()const {return *(uint32_t*)buffer;}
	NewMessageScope tryAlloc(uint16_t length, uint8_t nbreads);
	void printBufferStatus(FILE *f)const;
	uint8_t* readAt(uint32_t& position, uint16_t& length);
	void flushReader(uint32_t position);
};



/*
// hashmap which is accessed by several threads, so to possibly contain mutexes
template<class Key, class DATA, class HashFnc>
class MutexHashmap{
    typedef typename ExCo<Key>::TYPE KeyType;
    unsigned int hashpos(unsigned int seed) const;
    void swap_and_pop(unsigned int ite);
    void rehash(unsigned char _new_mag);
    DATA master;
    public:
    typedef unsigned int ITERATOR_TYPE;
    Vector< pair< KeyElem<Key, DATA> , unsigned int > > heap;
    unsigned int* hash;
    unsigned char hash_mag;

    MutexHashmap(): hash_mag(0), master(0){}
    bool grab(const Key &); // assumes entry exists, use "find" otherwise
    bool free(const Key &); // assumes entry exists, use "find" otherwise
};
*/
// maintains linked-list of categories, which are gradually ordered at the structure is queried
/*
template<class Key, class Data, class HashFnc>
class SortedHashmap{
    myHashmap<Category, unsigned int> categories;
    Vector< Vector<Data> > data;
public:

    Vector< pair< KeyElem<Key, Data> , pair<unsigned int, unsigned int> > > heap;
    unsigned int* hash;
    unsigned char hash_mag;

    CategoryHashmap(): hash_mag(0){}
    void insert(const Key&, const Data&, const Category&);
    void remove(const Key&);
    const Vector<Data>& getCategory() const; // requires *no* alteration to hashmap
};*/


	// or mantissa [0,1] range

	/*template<class C = int>
	class angle{

	public:
	C ang;

	angle();
	angle(C const & value);
	operator double() const;
	operator complex() const;

	double real(){return(cos((double)(*this)));}
	double imag(){return(sin((double)(*this)));}
	complex getcomplex() const;
};*/

#include "VLongint.hpp"
#include "RBTree.hpp"
#include "Hierachical.hpp"

#include "Advanced.hpp"

} // end of namespace

#endif // _defined_LFHAdvanced

