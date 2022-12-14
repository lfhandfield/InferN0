/*
 * HashMap.hpp
 *
 * Copyright (C) 2013 Louis-Francois Handfield
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
 */


#undef LFHTEMP
#define LFHTEMP template <class C>


LFHTEMP Vector<C>::Vector(const C* _data, int _size){
    unsigned int i;
    setSize_init(_size);
    for(i=0;i<_size;i++) new(darray + i) C(_data[i]);
}
LFHTEMP Vector<C>::Vector(std::initializer_list<C> ivl){setSize_init(ivl.size()); C* tar = darray; const C* sour = std::begin(ivl); for(uint32_t i=0;i<asize;i++) new(tar++) C(*(sour++));}
LFHTEMP template<class O> Vector<C>::Vector(std::initializer_list<O> ivl){setSize_init(ivl.size()); C* tar = darray; const O* sour = std::begin(ivl); for(uint32_t i=0;i<asize;i++) new(tar++) C(*(sour++));}

LFHTEMP Vector<C>& Vector<C>::toMemfree(){if (asize != 0) {delete[](darray); asize =0;} return(*this);}
LFHTEMP void Vector<C>::show(FILE* f_out, int l) const{
    unsigned int s = getSize();
    switch(l){
        case 1:
        case 0:
            for(unsigned int i=0;i<s;i++) {
                if (i != 0) fprintf(f_out,"\t");
                ExOp::show(darray[i],f_out,2);
            }
            if (l==0) fprintf(f_out,"\n");
            break;
        case 2:
            fprintf(f_out,"[");
            for(unsigned int i=0;i<s;i++) {
                if (i != 0) fprintf(f_out,";");
                ExOp::show(darray[i],f_out,3);
            }
            fprintf(f_out,"]");
            break;
        default:
            fprintf(f_out,"(");
            for(unsigned int i=0;i<s;i++) {
                if (i != 0) fprintf(f_out,",");
                ExOp::show(darray[i],f_out,4);
            }
            fprintf(f_out,")");
            break;
    }
}
/*LFHTEMP void Vector<C>::save(FILE* f) const {
	unsigned int s = getSize();
	int j;
	if ((j=fwrite(&s,sizeof(unsigned int), 1, f)) != 1) {printf("warning, got %i != 1 for fwrite",j);}
	if (ExCo<C>::IsPOD) if (j=fwrite(darray,sizeof(C),getSize(),f)) != getSize()) {printf("warning, got %i != %i for fwrite",j, getSize());}
	else{
		unsigned int s = getSize();
		for(unsigned int i=0;i<s;i++) ExOp::save(darray[i],f);
}	}*/
LFHTEMP template<class O> std::vector<O>& Vector<C>::wrStdVector(std::vector<O>& fout)const{
    fout.clear();
    fout.reserve(getSize());
    for(uint32_t i =0;i<getSize();i++) fout.push_back(darray[i]);
return fout;}
LFHTEMP ERRCODE Vector<C>::save(FILE* f) const {
	uint32_t s = getSize();
	ERRCODE fout = (fwrite(&s,sizeof(uint32_t), 1, f) != 1) ? 1 : 0;
	if (ExCo<C>::IsPOD::value) return fout | ((fwrite(darray,sizeof(C),getSize(),f) != 1) ? 1 : 0);
    //    if (j != getSize()) fprintf(stderr, "warning, got %i != %i for fwrite",j, getSize());
	else{
		unsigned int s = getSize();
		for(unsigned int i=0;i<s;i++) fout |= ExOp::save(darray[i],f);
		return fout;
}	}
LFHTEMP ERRCODE Vector<C>::load(FILE* f) {
	uint32_t s;
	ERRCODE fout = (fread(&s,sizeof(uint32_t), 1, f) != 1) ? 1 : 0;
	this->setSize(s);
	if (ExCo<C>::IsPOD::value){
        if (fread(darray,sizeof(C),s,f) != s) fout |= 1;
	}else for(unsigned int i=0;i<s;i++) fout |= ExOp::load(darray[i],f);
return fout;}


/*LFHTEMP void Vector<C>::setLinkMemory(void* new_mem){
	LinkAssert< ExCo<C>::NeedsAddLink == false > ass;

	if (asize == 0) darray = new_mem;
	else{
		LinkMem.setOwner(new_mem, darray, darray + sizeof(C) * (*this).allocatedSize() )(darray);
	}

}*/


LFHTEMP unsigned int Vector<C>::allocatedSize() const{

	uint32_t tmp = (asize-1) & 0x7FFFFFFF;
	if (tmp == 0) return(0);
	else tmp--;
	tmp |= tmp >> 8;
	tmp |= tmp >> 4;
	tmp |= tmp >> 2;
	tmp |= tmp >> 1;
	tmp++;
	return  (asize & 0x80000000) ? (tmp << 1) : tmp;
}


LFHTEMP C&  Vector<C>::push_back(){
	uint32_t rsize = (asize & 0x7FFFFFFF);
	uint32_t i;
	if (((asize + 1)^(asize))> rsize){
		if (asize & 0x80000000) asize &= 0x7FFFFFFf;
		else{
				// needs up alloc
				/*printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);
				C* hahaha = new C[256];
				printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);
				delete[](hahaha);

				printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);*/
				C* swp = new C[((rsize+1) << 1)];
				LFH_NICE_ALLOCERROR(swp, "Could not allocate for Vector<C>::push_back\n(size = %i)\n", (int)(((rsize+1) << 1)* sizeof(C)));

			   // printf("triedand got %i!\n", swp ); fflush(stdout);

				if (rsize> 0){

					for(i=0;i<rsize;i++){
						// printf("cpytime! %i\n", i); fflush(stdout);
						ExOp::toMemmove(swp[i], darray[i]);
						//swp[i]= darray[i];
					}
					delete[](darray);
				}
				darray = swp;

		}
	}
	return(darray[((asize++) & 0x7FFFFFFF)]);
}

LFHTEMP C&  Vector<C>::push_back(C && newitem){
	uint32_t rsize = (asize & 0x7FFFFFFF);
	uint32_t i;
	if (((asize + 1)^(asize))> rsize){
		if (asize & 0x80000000) asize &= 0x7FFFFFFf;
		else{
				// needs up alloc
				/*printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);
				C* hahaha = new C[256];
				printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);
				delete[](hahaha);

				printf("up-alloec in vector! %i planned\n", ((rsize+1) << 1) ); fflush(stdout);*/
				C* swp = new C[((rsize+1) << 1)];
				LFH_NICE_ALLOCERROR(swp, "Could not allocate for Vector<C>::push_back\n(size = %i)\n", (int)(((rsize+1) << 1)* sizeof(C)))

			   // printf("triedand got %i!\n", swp ); fflush(stdout);

				if (rsize> 0){

					for(i=0;i<rsize;i++){
						// printf("cpytime! %i\n", i); fflush(stdout);
						ExOp::toMemmove(swp[i], darray[i]);
						//swp[i]= darray[i];
					}
					delete[](darray);
				}
				darray = swp;
		}
	}
	return(darray[((asize++) & 0x7FFFFFFF)] = newitem);
}

LFHTEMP void Vector<C>::pop_back(){
	uint32_t rsize = (asize & 0x7FFFFFFF);
	if (((asize)^(asize-1))> rsize-1){
		if (rsize-1 == 0){
			delete[](darray);
			darray = NULL;
			asize =0;
			return;
		}else{
			if ((asize & 0x80000000) == 0) asize |= 0x80000000;
			else{
				// needs down alloc
				{//:
					C* swp = new C[(rsize) << 1];
					LFH_NICE_ALLOCERROR(swp,"Could not allocate for Vector<C>::push_back\n(size = %i)\n", (int)(((rsize) << 1)* sizeof(C)))
					uint32_t i;
					for(i=0;i<rsize;i++){
						ExOp::toMemmove(swp[i], darray[i]);
					}
					delete[](darray);
					//if (ExCo<C>::NeedsAddLink) LinkMem.moveOwner(darray, swp, swp+ ((rsize) << 1));
					darray = swp;
				}//:
			}
		}
	}
	asize--;
	ExOp::toMemfree(darray[asize & 0x7FFFFFFF]);
}

LFHTEMP void Vector<C>::insert_at_position(unsigned int pos, const C& entry){ // may be slow, shift most entries at any beyong "pos" by 1
uint32_t rsize = (asize & 0x7FFFFFFF);
uint32_t i;
if (((asize + 1)^(asize))> rsize){
	if (asize & 0x80000000) asize &= 0x7FFFFFFf;
	else{ // could be more clever, oh well
			// needs up alloc
			C* swp = new C[((rsize+1) << 1)];
			LFH_NICE_ALLOCERROR(swp,"Could not allocate for Vector<C>::push_back\n(size = %i)\n", (int)(((rsize+1) << 1)* sizeof(C)))
            if (rsize> 0){
				for(i=0;i<rsize;i++){
					ExOp::toMemmove(swp[i], darray[i]);
					//swp[i]= darray[i];
				}
				delete[](darray);
			}
			darray = swp;
	}
}
rsize = (asize & 0x7FFFFFFF);
for(rsize--; rsize != pos; rsize--) ExOp::toMemmove(darray[rsize], darray[rsize-1]);
darray[pos] = entry;
}
LFHTEMP void Vector<C>::pop_swap(unsigned int w){
	if (w+1 != (asize & 0x7FFFFFFF) ) ExOp::toMemmove(darray[w], darray[(asize & 0x7FFFFFFF)-1]);
	uint32_t rsize = (asize & 0x7FFFFFFF);
	if (((asize)^(asize-1))> rsize-1){
		if (rsize-1 == 0){
			delete[](darray);
			// if (ExCo<C>::NeedsAddLink) darray =(C*) LinkMem.removeOwner(darray); else
			darray = NULL;
			asize =0;
			return;
		}else{
			if ((asize & 0x80000000) == 0) asize |= 0x80000000;
			else{
				// needs down alloc
				{//:
					C* swp = new C[(rsize) << 1];

					LFH_NICE_ALLOCERROR(swp, "Could not allocate for Vector<C>::pop_swap\n")
					uint32_t i;
					for(i=0;i<rsize;i++){
						ExOp::toMemmove(swp[i], darray[i]);
					}
					delete[](darray);
					//if (ExCo<C>::NeedsAddLink) LinkMem.moveOwner(darray, swp, swp+ ((rsize) << 1));
					darray = swp;
				}//:
			}
		}
	}
	asize--;
	ExOp::toMemfree(darray[asize & 0x7FFFFFFF]);
}
LFHTEMP void Vector<C>::setSize_init(uint32_t nsize) {
    asize = nsize;
	if (nsize == 0) return;
	uint32_t i =1;
	while(i<= nsize) i = (i<< 1);
	darray = new C[i];
	LFH_NICE_ALLOCERROR(darray,"Could not allocate for Vector<C>::push_back\n(size = %i)\n",(int)(i*sizeof(C)))
}
LFHTEMP Vector<C>& Vector<C>::setSize(uint32_t nsize) {
	if (nsize == (asize &0x7FFFFFFF)) return *this;
	if (asize != 0) delete[](darray);
	asize = nsize;
	if (nsize == 0) return *this;
	uint32_t i =1;
	while(i<= nsize) i = (i<< 1);
	darray = new C[i];
	LFH_NICE_ALLOCERROR(darray,"Could not allocate for Vector<C>::push_back\n(size = %i)\n",(int)(i*sizeof(C)))
return *this;}
LFHTEMP Vector<C>& Vector<C>::DownSize(unsigned int nsize) {
	if (nsize == (unsigned int)(asize &0x7FFFFFFF)) return *this;

	if (nsize == 0) return *this;
	uint32_t i =1;
	while(i< nsize) i = (i<< 1);
	if (i < (asize & 0x7FFFFFFF)){
		C* tswap = new C[i];
		LFH_NICE_ALLOCERROR(tswap,"Could not allocate for Vector<C>::push_back\n(size = %i)\n", (int)i)
		asize = (int)nsize;
		for(i=0;i< nsize;i++) ExOp::toMemmove(tswap[i], darray[i]);
		delete[](darray);
		darray = tswap;
	}else{ // no resize
		if ((asize & 0x7FFFFFFF) < nsize) {fprintf(stderr,"Illegal Vector DownSize!"); LFH_exit(1);}
		asize = (nsize & 0x7FFFFFFF) | (asize & 0x80000000);
	}

	//printf("allocating (%x) (owner=%x, size=%i)\n",(int)darray, (int)this,i);fflush(stdout);
	//	printf("allocating (%x) (owner=%x)\n",(int)darray, (int)this);fflush(stdout);

return *this;}
LFHTEMP Vector<C>& Vector<C>::upSize(unsigned int nsize) {
    uint32_t tmp = asize & 0x7FFFFFFF;
	ExOp::toLeftFlood(tmp);
	if (nsize <= tmp) {asize = (asize & 0x80000000) | nsize; return *this;}
	if (asize & 0x80000000){
        if (nsize <= ((tmp<<1)|1)) {asize = nsize; return *this;}
	}
    // change of allocation needed
    //printf("Up aLLoc needed!\n"); fflush(stdout);
	tmp =1;
	while(tmp<= nsize) tmp = (tmp<< 1);
	C* tmpar = new C[tmp];
	LFH_NICE_ALLOCERROR(tmpar,"Could not allocate for Vector<C>::push_back\n(size = %i)\n",(int)(tmp*sizeof(C)))
    asize &= 0x7FFFFFFF;
    for(tmp=0;tmp< asize;tmp++) ExOp::toMemmove(tmpar[tmp], darray[tmp]);
    if (asize != 0) delete[](darray);
    asize = nsize;
    darray = tmpar;
return *this;}
LFHTEMP Vector<C>& Vector<C>::toAppend(const Vector<C>& other){
    uint32_t ts = asize & 0x7FFFFFFF;
    this->upSize(ts + other.getSize());
    for(uint32_t i=0;i<other.getSize();i++) darray[i+ts] = other[i];
return(*this);}
LFHTEMP Vector<C>& Vector<C>::toMemappend(Vector<C>& other){
    uint32_t ts = asize & 0x7FFFFFFF;
    this->upSize(ts + other.getSize());
    for(uint32_t i=0;i<other.getSize();i++) ExOp::toMemmove(darray[i+ts],other[i]);
return(*this);}
LFHTEMP Vector<C>& Vector<C>::toMemappendInverted(Vector<C>& other){
    uint32_t ts = asize & 0x7FFFFFFF;
    this->upSize(ts + other.getSize());
    for(uint32_t i=0;i<other.getSize();i++) ExOp::toMemmove(darray[other.getSize() - i - 1 +ts],other[i]);
return(*this);}


LFHTEMP uint32_t Vector<C>::size() const {return (uint32_t)(asize & 0x7FFFFFFF);}
LFHTEMP uint32_t Vector<C>::getSize() const {return (uint32_t)(asize & 0x7FFFFFFF);}
LFHTEMP	Vector<C>& Vector<C>::operator=(const Vector<C> & v){
	setSize(v.getSize());
	for(uint32_t i=0;i< (unsigned int)(asize &0x7FFFFFFF) ;i++) darray[i] = v.darray[i];
return(*this);}
LFHTEMP Vector<C>& Vector<C>::toMemmove(Vector<C> & v){
    if (asize != 0) delete[](darray);
	asize = v.asize;
	v.asize = 0;
	darray = v.darray;
return(*this);}
LFHTEMP Vector<C>& Vector<C>::toMemswap(Vector<C> & v){
    union{
        uint32_t asizeswp;
        C* darrayswp;
    };
	asizeswp = asize; asize = v.asize; v.asize = asizeswp;
	darrayswp = darray; darray = v.darray; v.darray = darrayswp;
return(*this);}
LFHTEMP template<class OC> Vector<C>& Vector<C>::operator=(const Vector<OC> & v){
	setSize(v.getSize());
	for(uint32_t i=0;i<(asize & 0x7FFFFFFF);i++) ExOp::toClone(darray[i],v.darray[i]);
return(*this);}
LFHTEMP template<class OC> Vector<C>& Vector<C>::toConvert(const Vector<OC> & v){
	//	printf("vector assignement\n");
	setSize(v.getSize());
	for(uint32_t i=0;i<(asize & 0x7FFFFFFF);i++)  ExOp::toConvert(darray[i],v.darray[i]);
return(*this);}
LFHTEMP
template<class D> Vector<C>& Vector<C>::operator+=(const Vector<D> & v){
	for(uint32_t i=0;i<(asize & 0x7FFFFFFF);i++) darray[i] += v[i];
return(*this);}
LFHTEMP template<class D> Vector<C>& Vector<C>::operator+=(const D & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] += v;
return(*this);}
LFHTEMP template<class D> Vector<C>& Vector<C>::operator-=(const Vector<D> & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] -= v[i];
return(*this);}
LFHTEMP template<class D> Vector<C>& Vector<C>::operator-=(const D & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] -= v;
return(*this);}
LFHTEMP
template<class D> Vector<C>& Vector<C>::operator*=(const Vector<D> & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] *= v[i];
	return(*this);
}
LFHTEMP
template<class D> Vector<C>& Vector<C>::operator*=(const D & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] *= v;
	return(*this);
}
LFHTEMP
template<class D> Vector<C>& Vector<C>::operator/=(const Vector<D> & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] /= v[i];
	return(*this);
}
LFHTEMP
template<class D> Vector<C>& Vector<C>::operator/=(const D & v){
	for(unsigned int i=0;i<(unsigned int)(asize & 0x7FFFFFFF);i++) darray[i] /= v;
	return(*this);
}
/** \brief get order of element for array, so that Data[order[x]] is ordered, O(n) if already sorted
 *
 * \return void
 *
 */
LFHTEMP template<class O> Vector<C>& Vector<C>::toOrderOf(const O* data){
    C swap;
    int i;
    int j;
    for(i=0;i< (asize & 0x7FFFFFFF);i++) darray[i] = (C)i;
    Vector<unsigned int> maxi;
    if ( (asize& 0x7FFFFFFF) == 0) return *this;
    maxi.push_back(asize & 0x7FFFFFFF);
    int mini = 0;
    while(maxi.getSize() > 0u){
        //		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
        switch(maxi.last() - mini){
            case 0:
            case 1:
                mini = maxi.last()+1; maxi.pop_back();
                break;
            case 2:
                if (ExOp::isGT(data[darray[mini]], data[darray[mini+1]])){
                    ExOp::toMemmove(swap, darray[mini+1]);
                    ExOp::toMemmove(darray[mini+1], darray[mini]);
                    ExOp::toMemmove(darray[mini], swap);
                }
                mini += 3;  maxi.pop_back();
                break;
            case 3:
                if (ExOp::isLT(data[darray[mini]], data[darray[mini+2]])){
                    if (ExOp::isLT(data[darray[mini+1]], data[darray[mini]])){
                        ExOp::toMemmove(swap,darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], darray[mini]);
                        ExOp::toMemmove(darray[mini], swap);
                    }else if (ExOp::isLT(data[darray[mini+2]] , data[darray[mini+1]])){
                        ExOp::toMemmove(swap,darray[mini+2]);
                        ExOp::toMemmove(darray[mini+2],darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], swap);
                    }
                }else{
                    ExOp::toMemmove(swap, darray[mini+2]);
                    if (ExOp::isLT(data[darray[mini+1]], data[darray[mini+2]])){
                        ExOp::toMemmove(darray[mini+2], darray[mini]);
                        ExOp::toMemmove(darray[mini], darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], swap);
                    }else{
                        if (ExOp::isLT(data[darray[mini]] , data[darray[mini+1]])){
                            ExOp::toMemmove(darray[mini+2], darray[mini+1]);
                            ExOp::toMemmove(darray[mini+1], darray[mini]);
                        }else{
                            ExOp::toMemmove(darray[mini+2], darray[mini]);

                        }
                        ExOp::toMemmove(darray[mini], swap);
                    }
                }
                mini += 4;  maxi.pop_back();
                break;
            default:
                j = maxi.last()-1;
                i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
                switch( ((ExOp::isLT(data[darray[mini]], data[darray[i]])) ? 1 : 0) | (ExOp::isLT(data[darray[j-1]], data[darray[j]]) ? 2 : 0) ){
                    case 0:
                        if (ExOp::isLT(data[darray[mini]] , data[darray[j-1]])){
                            // daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
                            ExOp::toMemmove(swap, darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[i]);
                            ExOp::toMemmove(darray[i], darray[j]);
                            ExOp::toMemmove(darray[j], darray[j-1]);
                        }else{
                            // daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
                            ExOp::toMemmove(swap, darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[j]);
                            ExOp::toMemmove(darray[j], swap);
                            ExOp::toMemmove(swap, darray[j-1]);
                        }
                        break;
                    case 1:

                        if (ExOp::isLT(data[darray[i]] ,data[darray[j-1]])){
                            // daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
                            ExOp::toMemmove(swap, darray[i]);
                            ExOp::toMemmove(darray[i], darray[j]);
                            ExOp::toMemmove(darray[j], darray[j-1]);
                        }else{
                            // daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
                            ExOp::toMemmove(swap, darray[i]);
                            ExOp::toMemmove(darray[i], darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[j]);
                            ExOp::toMemmove(darray[j], swap);
                            ExOp::toMemmove(swap, darray[j-1]);
                        }
                        break;
                    case 2:
                        if (ExOp::isLT(data[darray[mini]] , data[darray[j]])){
                            // daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
                            ExOp::toMemmove(swap, darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[i]);
                            ExOp::toMemmove(darray[i], darray[j-1]);
                        }else{
                            // daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
                            ExOp::toMemmove(swap, darray[j]);
                            ExOp::toMemmove(darray[j], darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[j-1]);
                        }
                        break;
                    case 3:
                        if (ExOp::isLT(data[darray[i]] , data[darray[j]])){
                            // daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
                            ExOp::toMemmove(swap, darray[i]);
                            ExOp::toMemmove(darray[i], darray[j-1]);
                        }else{
                            // daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
                            swap = darray[j];
                            ExOp::toMemmove(darray[j], darray[i]);
                            ExOp::toMemmove(darray[i], darray[mini]);
                            ExOp::toMemmove(darray[mini], darray[j-1]);
                        }
                        break;
                }
                //		printf("pivot! %i!\n", swap);
                i = mini+1;
                j--;
                //		printf("pospos! %i,%i!\n", i,j);
                while(true){
                    while(ExOp::isLE(data[darray[i]] , data[swap])) {
                        if ((++i) == j) break;
                        if (data[darray[i]] >= data[swap]) break;
                        if ((++i) == j) break;
                    }
                    //			printf("%i,%i mint\n", i,j);
                    if (i >= j) break;
                    ExOp::toMemmove(darray[j], darray[i]);
                    while(ExOp::isGE(data[darray[j]], data[swap])) {
                        if (i == (--j)) break;
                        if (data[darray[j]] <= data[swap]) break;
                        if (i == (--j)) break;
                    }
                    //			printf("%i,%i maxt\n", i,j);
                    if (i >= j) break;
                    ExOp::toMemmove(darray[i], darray[j]);
                }
                ExOp::toMemmove(darray[i], swap);
                maxi.push_back(i);
        }
    }
return *this;}
/** \brief sort array, O(n) if already sorted
 *
 * \return This Vector
 *
 */
LFHTEMP Vector<C>& Vector<C>::toRandPermute(){
	if (this->getSize() <= 1) return *this;
	C swappos;
	ExOp::toMemmove(swappos, darray[0]);
	uint32_t dasize = asize & 0x7FFFFFFF;
	uint32_t r,i;
	for(i=0; i < dasize - 1;i++){
		r = rand() % (dasize - i);
		if (r != 0) {
			ExOp::toMemmove(darray[i], darray[i+r]);
			ExOp::toMemmove(darray[i+r], darray[i+1]);
		}else {
			ExOp::toMemmove(darray[i], swappos);
			ExOp::toMemmove(swappos, darray[i+1]);
		}
	}
	ExOp::toMemswap(darray[i], swappos);
return *this;}

LFHTEMP Vector<C>& Vector<C>::sort(){return this->sortWithinRange(0,asize & 0x7FFFFFFF);}
/** \brief sort array, O(n) if already sorted, and computes a sum with number of tied entries within range
 *
 * \return Integer with value \sum t_i (t_i + 1) / 2
 *
 */
LFHTEMP uint32_t Vector<C>::sortWithinRangeAndSumPerm(int first, int limit){
    if ( first+1 >= limit) return 0;
    this->sortWithinRange(first, limit);
    int i,j;
    uint32_t fout =0;
    for(i = first,j= first+1; j < limit;j++){
        if (darray[j] != darray[i]) {
            fout += ((j-i) * (j-i -1))>>1;
            i = j;
        }
    }
    fout += ((j-i) * (j-i -1))>>1;
return fout;}
/** \brief sort array, O(n) if already sorted
 *
 * \return void
 *
 */
LFHTEMP Vector<C>& Vector<C>::sortWithinRange(int first, int limit){
    C swap;
    int i;
    int j;
    Vector<unsigned int> maxi;
    if ( first+1 >= limit) return *this;
    maxi.push_back(limit);
    int mini = first;
    while(maxi.getSize() > 0u){
        //		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
        switch(maxi.last() - mini){
		case 0:
		case 1:
			mini = maxi.last()+1; maxi.pop_back();
        break; case 2:
			if (ExOp::isGT(darray[mini], darray[mini+1])){
				ExOp::toMemmove(swap, darray[mini+1]);
				ExOp::toMemmove(darray[mini+1], darray[mini]);
				ExOp::toMemmove(darray[mini], swap);
			}
			mini += 3;  maxi.pop_back();
        break; case 3:
			if (ExOp::isLT(darray[mini], darray[mini+2])){
				if (darray[mini+1] < darray[mini]){
					ExOp::toMemmove(swap,darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], darray[mini]);
					ExOp::toMemmove(darray[mini], swap);
				}else if (ExOp::isLT(darray[mini+2] , darray[mini+1])){
					ExOp::toMemmove(swap,darray[mini+2]);
					ExOp::toMemmove(darray[mini+2],darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}
			}else{
				ExOp::toMemmove(swap, darray[mini+2]);
				if (ExOp::isLT(darray[mini+1] , darray[mini+2])){
					ExOp::toMemmove(darray[mini+2], darray[mini]);
					ExOp::toMemmove(darray[mini], darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}else{
					if (ExOp::isLT(darray[mini] , darray[mini+1])){
						ExOp::toMemmove(darray[mini+2], darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], darray[mini]);
					}else{
						ExOp::toMemmove(darray[mini+2], darray[mini]);

					}
					ExOp::toMemmove(darray[mini], swap);
				}
			}
			mini += 4;  maxi.pop_back();
        break; default:
			j = maxi.last()-1;
			i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
			switch( ((ExOp::isLT(darray[mini], darray[i])) ? 1 : 0) | (ExOp::isLT(darray[j-1], darray[j]) ? 2 : 0) ){
            case 0:
                if (ExOp::isLT(darray[mini] , darray[j-1])){
                    // daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
                    ExOp::toMemmove(swap, darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[i]);
                    ExOp::toMemmove(darray[i], darray[j]);
                    ExOp::toMemmove(darray[j], darray[j-1]);
                }else{
                    // daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
                    ExOp::toMemmove(swap, darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[j]);
                    ExOp::toMemmove(darray[j], swap);
                    ExOp::toMemmove(swap, darray[j-1]);
                }
            break; case 1:

                if (ExOp::isLT(darray[i] ,darray[j-1])){
                    // daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
                    ExOp::toMemmove(swap, darray[i]);
                    ExOp::toMemmove(darray[i], darray[j]);
                    ExOp::toMemmove(darray[j], darray[j-1]);
                }else{
                    // daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
                    ExOp::toMemmove(swap, darray[i]);
                    ExOp::toMemmove(darray[i], darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[j]);
                    ExOp::toMemmove(darray[j], swap);
                    ExOp::toMemmove(swap, darray[j-1]);
                }
            break; case 2:
                if (ExOp::isLT(darray[mini] , darray[j])){
                    // daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
                    ExOp::toMemmove(swap, darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[i]);
                    ExOp::toMemmove(darray[i], darray[j-1]);
                }else{
                    // daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
                    ExOp::toMemmove(swap, darray[j]);
                    ExOp::toMemmove(darray[j], darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[j-1]);
                }
            break; case 3:
                if (ExOp::isLT(darray[i] , darray[j])){
                    // daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
                    ExOp::toMemmove(swap, darray[i]);
                    ExOp::toMemmove(darray[i], darray[j-1]);
                }else{
                    // daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
                    swap = darray[j];
                    ExOp::toMemmove(darray[j], darray[i]);
                    ExOp::toMemmove(darray[i], darray[mini]);
                    ExOp::toMemmove(darray[mini], darray[j-1]);
                }
            break;}
			//		printf("pivot! %i!\n", swap);
			i = mini+1;
			j--;
			//		printf("pospos! %i,%i!\n", i,j);
			while(true){
				while(ExOp::isLE(darray[i] , swap)) {
					if ((++i) == j) break;
					if (darray[i] >= swap) break;
					if ((++i) == j) break;
				}
				//			printf("%i,%i mint\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[j], darray[i]);
				while(ExOp::isGE(darray[j], swap)) {
					if (i == (--j)) break;
					if (darray[j] <= swap) break;
					if (i == (--j)) break;
				}
				//			printf("%i,%i maxt\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[i], darray[j]);
			}
			ExOp::toMemmove(darray[i], swap);
			maxi.push_back(i);
        }
    }
return *this;}

/** \brief sort array, O(n) if already sorted, and report the number of inconsistencies in the ordering
 *
 * \return void
 *
 */
LFHTEMP Vector<C>& Vector<C>::sortAndCountPermutation(uint32_t &nbperm){
    C swap;
    int i,j,k;
    Vector<unsigned int> maxi;
    nbperm =0;
    maxi.push_back(asize & 0x7FFFFFFF);
    if ( maxi.last() < 2) return *this;
    int mini = 0;
    Vector<C> swapbuf;
    while(maxi.getSize() > 0u){
        //		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
        //printf("maxi.last() - mini is %i - %i = %i\n", maxi.last(), mini,maxi.last() - mini);
        //this->show();
        switch(maxi.last() - mini){
            case 0:
            case 1:
                mini = maxi.last()+1; maxi.pop_back();
                break;
            case 2:
                if (ExOp::isGT(darray[mini], darray[mini+1])){
                    ExOp::toMemmove(swap, darray[mini+1]);
                    ExOp::toMemmove(darray[mini+1], darray[mini]);
                    ExOp::toMemmove(darray[mini], swap);
                    nbperm++;
                }
                mini += 3;  maxi.pop_back();
                break;
            case 3:
                if (ExOp::isGT(darray[mini],darray[mini+2])){
                    ExOp::toMemmove(swap, darray[mini+2]); nbperm += 2;
                    if (ExOp::isGT(darray[mini+1],darray[mini+2])){
                        if (ExOp::isGT(darray[mini],darray[mini+1])){
                            ExOp::toMemmove(darray[mini+2], darray[mini]);
                            nbperm++;
                        }else{
                            ExOp::toMemmove(darray[mini+2], darray[mini+1]);
                            ExOp::toMemmove(darray[mini+1], darray[mini]);
                        }
                        ExOp::toMemmove(darray[mini], swap);
                    }else{
                        ExOp::toMemmove(darray[mini+2], darray[mini]);
                        ExOp::toMemmove(darray[mini], darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], swap);
                    }
                }else{
                    if (ExOp::isGT(darray[mini], darray[mini+1])){
                        ExOp::toMemmove(swap,darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], darray[mini]);
                        ExOp::toMemmove(darray[mini], swap);
                        nbperm++;
                    }else if (ExOp::isGT(darray[mini+1], darray[mini+2])){
                        ExOp::toMemmove(swap,darray[mini+2]);
                        ExOp::toMemmove(darray[mini+2],darray[mini+1]);
                        ExOp::toMemmove(darray[mini+1], swap);
                        nbperm++;
                    }
                }
                mini += 4;  maxi.pop_back();
            break; default:
                j = maxi.last()-1;
                // select pivot among mini, mini+1, mini+2
                // min is placed, swap is obtained, max is at top for swap buffer
                swapbuf.setSize(j-mini);
                if (ExOp::isGT(darray[mini],darray[mini+2])){
                    nbperm += 2;
                    if (ExOp::isGT(darray[mini+1],darray[mini+2])){ // 312
                        if (ExOp::isGT(darray[mini],darray[mini+1])){ // 231
                            ExOp::toMemmove(swap, darray[mini+1]);
                            ExOp::toMemmove(swapbuf[0], darray[mini]);
                            nbperm++;
                        }else{ // 321
                            ExOp::toMemmove(swap, darray[mini]);
                            ExOp::toMemmove(swapbuf[0], darray[mini+1]);
                        }
                        ExOp::toMemmove(darray[mini], darray[mini+2]);
                    }else{
                        ExOp::toMemmove(swap, darray[mini+2]);
                        ExOp::toMemmove(swapbuf[0], darray[mini]);
                        ExOp::toMemmove(darray[mini], darray[mini+1]);
                    }
                }else{
                    if (ExOp::isGT(darray[mini],darray[mini+1])){ // 213
                        ExOp::toMemmove(swap,darray[mini]);
                        ExOp::toMemmove(swapbuf[0], darray[mini+2]);
                        ExOp::toMemmove(darray[mini], darray[mini+1]);
                        nbperm++;
                    }else if (ExOp::isGT(darray[mini+1],darray[mini+2])){ // 132
                        ExOp::toMemmove(swap,darray[mini+2]);
                        ExOp::toMemmove(swapbuf[0],darray[mini+1]);
                        nbperm++;
                    }else{ // 123
                        ExOp::toMemmove(swap,darray[mini+1]);
                        ExOp::toMemmove(swapbuf[0],darray[mini+2]);
                    }
                }
                j=mini+1; k=1;
                for(i = mini +3; i <  maxi.last(); i++){
                    if (ExOp::isLT(darray[i] , swap)) {ExOp::toMemmove(darray[j++], darray[i]); nbperm+= k;}
                    else ExOp::toMemmove(swapbuf[k++], darray[i]);
                }
                //printf("median is at %i\n", j);
                nbperm += j - mini - 1;
                maxi.push_back(j);
                ExOp::toMemmove(darray[j++], swap);
                for(i=0;i<k;i++){
                    ExOp::toMemmove(darray[j++], swapbuf[i]);
                }

        }
    }
return *this;}
/** \brief sort the first N elements in array, O(n) if already sorted
 *
 * \param uint32_t (number of element to sort)
 * \return void
 *
 */
LFHTEMP Vector<C>& Vector<C>::sortFirstN(uint32_t _tmpsize){
C swap;
int i;
int j;
Vector<unsigned int> maxi;
if ( (asize& 0x7FFFFFFF) < 2) return *this;
maxi.push_back(_tmpsize);
int mini = 0;
while(maxi.getSize() > 0u){
	//		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
	switch(maxi.last() - mini){
		case 0:
		case 1:
			mini = maxi.last()+1; maxi.pop_back();
			break;
		case 2:
			if (ExOp::isGT(darray[mini], darray[mini+1])){
				ExOp::toMemmove(swap, darray[mini+1]);
				ExOp::toMemmove(darray[mini+1], darray[mini]);
				ExOp::toMemmove(darray[mini], swap);
			}
			mini += 3;  maxi.pop_back();
			break;
		case 3:
			if (ExOp::isLT(darray[mini], darray[mini+2])){
				if (darray[mini+1] < darray[mini]){
					ExOp::toMemmove(swap,darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], darray[mini]);
					ExOp::toMemmove(darray[mini], swap);
				}else if (ExOp::isLT(darray[mini+2] , darray[mini+1])){
					ExOp::toMemmove(swap,darray[mini+2]);
					ExOp::toMemmove(darray[mini+2],darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}
			}else{
				ExOp::toMemmove(swap, darray[mini+2]);
				if (ExOp::isLT(darray[mini+1] , darray[mini+2])){
					ExOp::toMemmove(darray[mini+2], darray[mini]);
					ExOp::toMemmove(darray[mini], darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}else{
					if (ExOp::isLT(darray[mini] , darray[mini+1])){
						ExOp::toMemmove(darray[mini+2], darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], darray[mini]);
					}else{
						ExOp::toMemmove(darray[mini+2], darray[mini]);

					}
					ExOp::toMemmove(darray[mini], swap);
				}
			}
			mini += 4;  maxi.pop_back();
			break;
		default:
			j = maxi.last()-1;
			i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
			switch( ((ExOp::isLT(darray[mini], darray[i])) ? 1 : 0) | (ExOp::isLT(darray[j-1], darray[j]) ? 2 : 0) ){
				case 0:
					if (ExOp::isLT(darray[mini] , darray[j-1])){
						// daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[i]);
						ExOp::toMemmove(darray[i], darray[j]);
						ExOp::toMemmove(darray[j], darray[j-1]);
					}else{
						// daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j]);
						ExOp::toMemmove(darray[j], swap);
						ExOp::toMemmove(swap, darray[j-1]);
					}
					break;
				case 1:

					if (ExOp::isLT(darray[i] ,darray[j-1])){
						// daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[j]);
						ExOp::toMemmove(darray[j], darray[j-1]);
					}else{
						// daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j]);
						ExOp::toMemmove(darray[j], swap);
						ExOp::toMemmove(swap, darray[j-1]);
					}
					break;
				case 2:
					if (ExOp::isLT(darray[mini] , darray[j])){
						// daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[i]);
						ExOp::toMemmove(darray[i], darray[j-1]);
					}else{
						// daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
						ExOp::toMemmove(swap, darray[j]);
						ExOp::toMemmove(darray[j], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j-1]);
					}
					break;
				case 3:
					if (ExOp::isLT(darray[i] , darray[j])){
						// daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[j-1]);
					}else{
						// daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
						swap = darray[j];
						ExOp::toMemmove(darray[j], darray[i]);
						ExOp::toMemmove(darray[i], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j-1]);
					}
					break;
			}
			//		printf("pivot! %i!\n", swap);
			i = mini+1;
			j--;
			//		printf("pospos! %i,%i!\n", i,j);
			while(true){
				while(ExOp::isLE(darray[i] , swap)) {
					if ((++i) == j) break;
					if (darray[i] >= swap) break;
					if ((++i) == j) break;
				}
				//			printf("%i,%i mint\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[j], darray[i]);
				while(ExOp::isGE(darray[j], swap)) {
					if (i == (--j)) break;
					if (darray[j] <= swap) break;
					if (i == (--j)) break;
				}
				//			printf("%i,%i maxt\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[i], darray[j]);
			}
			ExOp::toMemmove(darray[i], swap);
			maxi.push_back(i);


	}
}
return *this;}
LFHTEMP template<class COMP> Vector<C>& Vector<C>::sort_comp(const COMP& comp){
C swap;
int i;
int j;
Vector<unsigned int> maxi;
if ( (asize& 0x7FFFFFFF) < 2) return *this;
maxi.push_back(asize & 0x7FFFFFFF);
int mini = 0;
while(maxi.getSize() > 0u){
	//		printf("doing (%i, %i)!\n", mini, maxi.top()-1);
	switch(maxi.last() - mini){
		case 0:
		case 1:
			mini = maxi.last()+1; maxi.pop_back();
			break;
		case 2:
			if (comp(darray[mini], darray[mini+1])>0){
				ExOp::toMemmove(swap, darray[mini+1]);
				ExOp::toMemmove(darray[mini+1], darray[mini]);
				ExOp::toMemmove(darray[mini], swap);
			}
			mini += 3;  maxi.pop_back();
			break;
		case 3:
			if (comp(darray[mini], darray[mini+2])<0){
				if (darray[mini+1] < darray[mini]){
					ExOp::toMemmove(swap,darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], darray[mini]);
					ExOp::toMemmove(darray[mini], swap);
				}else if (ExOp::isLT(darray[mini+2] , darray[mini+1])){
					ExOp::toMemmove(swap,darray[mini+2]);
					ExOp::toMemmove(darray[mini+2],darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}
			}else{
				ExOp::toMemmove(swap, darray[mini+2]);
				if (comp(darray[mini+1] , darray[mini+2])<0){
					ExOp::toMemmove(darray[mini+2], darray[mini]);
					ExOp::toMemmove(darray[mini], darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], swap);
				}else{
					if (comp(darray[mini] , darray[mini+1])<0){
						ExOp::toMemmove(darray[mini+2], darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], darray[mini]);
					}else{
						ExOp::toMemmove(darray[mini+2], darray[mini]);

					}
					ExOp::toMemmove(darray[mini], swap);
				}
			}
			mini += 4;  maxi.pop_back();
			break;
		default:
			j = maxi.last()-1;
			i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
			switch( ((comp(darray[mini], darray[i]) <0) ? 1 : 0) | (ExOp::isLT(darray[j-1], darray[j]) ? 2 : 0) ){
				case 0:
					if (comp(darray[mini] , darray[j-1])<0){
						// daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[i]);
						ExOp::toMemmove(darray[i], darray[j]);
						ExOp::toMemmove(darray[j], darray[j-1]);
					}else{
						// daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j]);
						ExOp::toMemmove(darray[j], swap);
						ExOp::toMemmove(swap, darray[j-1]);
					}
					break;
				case 1:

					if (comp(darray[i] ,darray[j-1])<0){
						// daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[j]);
						ExOp::toMemmove(darray[j], darray[j-1]);
					}else{
						// daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j]);
						ExOp::toMemmove(darray[j], swap);
						ExOp::toMemmove(swap, darray[j-1]);
					}
					break;
				case 2:
					if (comp(darray[mini] , darray[j])<0){
						// daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
						ExOp::toMemmove(swap, darray[mini]);
						ExOp::toMemmove(darray[mini], darray[i]);
						ExOp::toMemmove(darray[i], darray[j-1]);
					}else{
						// daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
						ExOp::toMemmove(swap, darray[j]);
						ExOp::toMemmove(darray[j], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j-1]);
					}
					break;
				case 3:
					if (comp(darray[i] , darray[j])<0){
						// daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
						ExOp::toMemmove(swap, darray[i]);
						ExOp::toMemmove(darray[i], darray[j-1]);
					}else{
						// daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
						swap = darray[j];
						ExOp::toMemmove(darray[j], darray[i]);
						ExOp::toMemmove(darray[i], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[j-1]);
					}
					break;
			}
			//		printf("pivot! %i!\n", swap);
			i = mini+1;
			j--;
			//		printf("pospos! %i,%i!\n", i,j);
			while(true){
				while(comp(darray[i] , swap)<=0) {
					if ((++i) == j) break;
					if (darray[i] >= swap) break;
					if ((++i) == j) break;
				}
				//			printf("%i,%i mint\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[j], darray[i]);
				while(comp(darray[j], swap)>=0) {
					if (i == (--j)) break;
					if (darray[j] <= swap) break;
					if (i == (--j)) break;
				}
				//			printf("%i,%i maxt\n", i,j);
				if (i >= j) break;
				ExOp::toMemmove(darray[i], darray[j]);
			}
			ExOp::toMemmove(darray[i], swap);
			maxi.push_back(i);


	}
}
return *this;}
LFHTEMP Vector<C>& Vector<C>::sort_unique(){
	//printf("got in sortunique! %i\n", asize & 0x7FFFFFFF);
	if ((asize & 0x7FFFFFFF) < 2) return *this;
	this->sort();
	//printf("got in sorted! %i\n", asize & 0x7FFFFFFF);
	C *hare, *turtle;
	hare = darray;
	for(hare = darray+1; (uint32_t)(hare - darray) < (asize & 0x7FFFFFFF) ; hare++){
		//printf("cmp %i and %i\n", (int)(hare - darray),(int)(hare - darray -1));
		if (ExOp::isEQ(hare[0],hare[-1])) break;
	}
	if ((uint32_t)(hare - darray) < (asize & 0x7FFFFFFF)){
		//printf("found turtle!\n");
		turtle = hare; turtle--;
		for(hare++; (uint32_t)(hare - darray) < (asize & 0x7FFFFFFF) ; hare++){
			if (ExOp::isNQ(*turtle,*hare)) {
				turtle++;ExOp::toMemmove(*turtle, *hare);
			}
		}
		turtle++;
		//printf("downsize %i\n",turtle - darray);
		this->DownSize(turtle - darray);
	} // else all were unique
return *this;}
LFHTEMP uint32_t Vector<C>::searchFor(const C& what, bool _is_sorted) const{
    if ((asize & 0x7FFFFFFF) == 0) return 0xFFFFFFFF;
    uint32_t a =0;
    uint32_t b = asize & 0x7FFFFFFF;
    if (_is_sorted){
        b--;
        uint32_t c;
        while(a != b){
            c = (a+b)>>1;
            if (ExOp::isGT(darray[c], what)) a = c+1;
            else b = c;
        }
        return (ExOp::isEQ(darray[a], what)) ? a : 0xFFFFFFFF;
    }else{
        for(;a<b;a++) if (ExOp::isEQ(darray[a], what)) return a;
        return 0xFFFFFFFF;
    }
}
LFHTEMP Vector<C>& Vector<C>::reverse(){
	if (asize == 0) return *this;
	C swap;
	C* i;
	C* j;
	for(i=darray,j=darray + ((asize & 0x7FFFFFFF)-1); i<j; i++,j--){
		ExOp::toMemmove(swap,*i);
		ExOp::toMemmove(*i,*j);
		ExOp::toMemmove(*j,swap);
	}
return *this;}
LFHTEMP Vector<C>& Vector<C>::reverseLast(uint32_t nbpermuted){
	if (asize == 0) return *this;
	C swap;
	C* i;
	C* j;
	for(i=darray + ((asize & 0x7FFFFFFF)-nbpermuted) ,j=darray + ((asize & 0x7FFFFFFF)-1); i<j; i++,j--){
		ExOp::toMemmove(swap,*i);
		ExOp::toMemmove(*i,*j);
		ExOp::toMemmove(*j,swap);
	}
return *this;}
LFHTEMP bool Vector<C>::issorted()const{
	int i;
	for(i=(asize & 0x7FFFFFFF)-2;i>=0;i--) if (ExOp::isGT(darray[i], darray[i+1])) break;
	return(i < 0);
}
LFHTEMP bool Vector<C>::issorted_decr()const{
	int i;
	for(i=(asize & 0x7FFFFFFF)-2;i>=0;i--) if (ExOp::isLT(darray[i],darray[i+1])) break;
	return(i < 0);
}
LFHTEMP template<class COMP> bool Vector<C>::issorted_comp(const COMP& comp) const{
	int i;
	for(i=(asize & 0x7FFFFFFF)-2;i>=0;i--) if (comp(darray[i],darray[i+1])>0) break;
	return(i < 0);
}

LFHTEMP Vector<C>& Vector<C>::sort_decr(){return(sortWithinRange_decreasing(0, asize & 0x7FFFFFFF));}
LFHTEMP Vector<C>& Vector<C>::sortWithinRange_decreasing(int first, int limit){
	C swap;
	int i;
	int j;
	Vector<unsigned int> maxi;
	if (limit - first < 2) return *this;

	maxi.push_back(limit);
	int mini = first;
	while(maxi.getSize()>0u){
		//		printf("doing (%i, %i)!\n", mini, maxi.last()-1);
		switch(maxi.last() - mini){
			case 0:
			case 1:
				mini = maxi.last()+1; maxi.pop_back();
				break;
			case 2:
				if (ExOp::isLT(darray[mini], darray[mini+1])){
					ExOp::toMemmove(swap, darray[mini+1]);
					ExOp::toMemmove(darray[mini+1], darray[mini]);
					ExOp::toMemmove(darray[mini], swap);
				}
				mini += 3;  maxi.pop_back();
				break;
			case 3:
				if (ExOp::isGT(darray[mini], darray[mini+2])){
					if (darray[mini+1] < darray[mini]){
						ExOp::toMemmove(swap,darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], darray[mini]);
						ExOp::toMemmove(darray[mini], swap);
					}else if (ExOp::isGT(darray[mini+2] , darray[mini+1])){
						ExOp::toMemmove(swap,darray[mini+2]);
						ExOp::toMemmove(darray[mini+2],darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], swap);
					}
				}else{
					ExOp::toMemmove(swap, darray[mini+2]);
					if (ExOp::isGT(darray[mini+1] , darray[mini+2])){
						ExOp::toMemmove(darray[mini+2], darray[mini]);
						ExOp::toMemmove(darray[mini], darray[mini+1]);
						ExOp::toMemmove(darray[mini+1], swap);
					}else{
						if (ExOp::isGT(darray[mini] , darray[mini+1])){
							ExOp::toMemmove(darray[mini+2], darray[mini+1]);
							ExOp::toMemmove(darray[mini+1], darray[mini]);
						}else{
							ExOp::toMemmove(darray[mini+2], darray[mini]);

						}
						ExOp::toMemmove(darray[mini], swap);
					}
				}
				mini += 4;  maxi.pop_back();
				break;
			default:
				j = maxi.last()-1;
				i = (mini + j) >> 1; // printf("%i\t%i\n", i, j);
				switch( ((ExOp::isGT(darray[mini], darray[i])) ? 1 : 0) | (ExOp::isLT(darray[j-1], darray[j]) ? 2 : 0) ){
					case 0:
						if (ExOp::isGT(darray[mini] , darray[j-1])){
							// daray[i] < darray[mini] < darray[j-1]  ?  darray[j]
							ExOp::toMemmove(swap, darray[mini]);
							ExOp::toMemmove(darray[mini], darray[i]);
							ExOp::toMemmove(darray[i], darray[j]);
							ExOp::toMemmove(darray[j], darray[j-1]);
						}else{
							// daray[j] < darray[j-1] < darray[mini]  ?  darray[i]
							ExOp::toMemmove(swap, darray[mini]);
							ExOp::toMemmove(darray[mini], darray[j]);
							ExOp::toMemmove(darray[j], swap);
							ExOp::toMemmove(swap, darray[j-1]);
						}
						break;
					case 1:

						if (ExOp::isGT(darray[i] ,darray[j-1])){
							// daray[mini] < darray[i] < darray[j-1]  ?  darray[j]
							ExOp::toMemmove(swap, darray[i]);
							ExOp::toMemmove(darray[i], darray[j]);
							ExOp::toMemmove(darray[j], darray[j-1]);
						}else{
							// daray[j] < darray[j-1] < darray[i]  ?  darray[mini]
							ExOp::toMemmove(swap, darray[i]);
							ExOp::toMemmove(darray[i], darray[mini]);
							ExOp::toMemmove(darray[mini], darray[j]);
							ExOp::toMemmove(darray[j], swap);
							ExOp::toMemmove(swap, darray[j-1]);
						}
						break;
					case 2:
						if (ExOp::isGT(darray[mini] , darray[j])){
							// daray[i] < darray[mini] < darray[j]  ?  darray[j-1]
							ExOp::toMemmove(swap, darray[mini]);
							ExOp::toMemmove(darray[mini], darray[i]);
							ExOp::toMemmove(darray[i], darray[j-1]);
						}else{
							// daray[j-1] < darray[j] < darray[mini]  ?  darray[i]
							ExOp::toMemmove(swap, darray[j]);
							ExOp::toMemmove(darray[j], darray[mini]);
							ExOp::toMemmove(darray[mini], darray[j-1]);
						}
						break;
					case 3:
						if (ExOp::isGT(darray[i] , darray[j])){
							// daray[mini] < darray[i] < darray[j]  ?  darray[j-1]
							ExOp::toMemmove(swap, darray[i]);
							ExOp::toMemmove(darray[i], darray[j-1]);
						}else{
							// daray[j-1] < darray[j] < darray[i]  ?  darray[mini]
							swap = darray[j];
							ExOp::toMemmove(darray[j], darray[i]);
							ExOp::toMemmove(darray[i], darray[mini]);
							ExOp::toMemmove(darray[mini], darray[j-1]);
						}
						break;
				}
				//		printf("pivot! %i!\n", swap);
				i = mini+1;
				j--;
				//		printf("pospos! %i,%i!\n", i,j);
				while(true){
					while(ExOp::isGE(darray[i] , swap)) {
						if ((++i) == j) break;
						if (darray[i] >= swap) break;
						if ((++i) == j) break;
					}
					//			printf("%i,%i mint\n", i,j);
					if (i >= j) break;
					ExOp::toMemmove(darray[j], darray[i]);
					while(ExOp::isLE(darray[j], swap)) {
						if (i == (--j)) break;
						if (darray[j] <= swap) break;
						if (i == (--j)) break;
					}
					//			printf("%i,%i maxt\n", i,j);
					if (i >= j) break;
					ExOp::toMemmove(darray[i], darray[j]);
				}
				ExOp::toMemmove(darray[i], swap);
				maxi.push_back(i);
    }	}


return *this;}
LFHTEMP uint32_t Vector<C>::findLE(const C& what)const{
    uint32_t a =0;
	uint32_t b = asize & 0x7FFFFFFF;
	if (b == 0) return 0xFFFFFFFF;
    b--;
    uint32_t c;
    while(a != b){
        c = (a+b)>>1;
        if (ExOp::isGT(darray[c], what)) a = c+1;
        else b = c;
    }
    return (ExOp::isGT(darray[a], what)) ? a : 0xFFFFFFFF;
}
LFHTEMP template<class D> uint32_t Vector<C>::findLE(const D& what)const{
    uint32_t a =0;
	uint32_t b = asize & 0x7FFFFFFF;
	if (b == 0) return 0xFFFFFFFF;
    b--;
    uint32_t c;
    while(a != b){
        c = (a+b)>>1;
        if (ExOp::isGT(darray[c], what)) a = c+1;
        else b = c;
    }
    return (ExOp::isGT(darray[a], what)) ? a : 0xFFFFFFFF;
}
LFHTEMP void Vector<C>::random_permute(){
	unsigned int i;
	unsigned int j = asize & 0x7FFFFFFF;
	C swap;
	if (j < 2) return;
	do{
		i = rand() % (j--);
		if (i == j) continue;
		ExOp::toMemmove(swap, darray[j]);
		ExOp::toMemmove(darray[j], darray[i]);
		ExOp::toMemmove(darray[i], swap);
	}while (j > 1);
}
LFHTEMP Vector<mycomplex> Vector<C>::bluesteinWindow(unsigned int size){
	unsigned int t = (size >> 2);
	int i;
	for(i=0;t != 0;i++) t >>=1;

	printf("%i is order\n",i); fflush(stdout);


	return(Vector<mycomplex>());
	switch(i){
		case 1: return(Vector<mycomplex>(Tuple<mycomplex,16>::bluesteinWindow(size)));
		case 2: return(Vector<mycomplex>(Tuple<mycomplex,32>::bluesteinWindow(size)));
		case 3: return(Vector<mycomplex>(Tuple<mycomplex,64>::bluesteinWindow(size)));
		case 4: return(Vector<mycomplex>(Tuple<mycomplex,128>::bluesteinWindow(size)));
		case 5: return(Vector<mycomplex>(Tuple<mycomplex,256>::bluesteinWindow(size)));
		case 6: return(Vector<mycomplex>(Tuple<mycomplex,512>::bluesteinWindow(size)));
		case 7: return(Vector<mycomplex>(Tuple<mycomplex,1024>::bluesteinWindow(size)));
		case 8: return(Vector<mycomplex>(Tuple<mycomplex,2048>::bluesteinWindow(size)));
		case 9: return(Vector<mycomplex>(Tuple<mycomplex,4096>::bluesteinWindow(size)));
		case 10: return(Vector<mycomplex>(Tuple<mycomplex,8192>::bluesteinWindow(size)));
		case 11: return(Vector<mycomplex>(Tuple<mycomplex,16384>::bluesteinWindow(size)));
		case 12: return(Vector<mycomplex>(Tuple<mycomplex,32768>::bluesteinWindow(size)));
		case 13: return(Vector<mycomplex>(Tuple<mycomplex,65536>::bluesteinWindow(size)));
		default:
			return(Vector<mycomplex>(Tuple<mycomplex,8>::bluesteinWindow(size)));

	}

}
/*
 LFHTEMP void Vector<C>::bluesteinWindow(complex *& fout, unsigned int size){
 unsigned int t = (size >> 2);
 int i;
 for(i=0;t != 0;i++) t >>=1;
 printf("%i is order\n",i); fflush(stdout);
 switch(i){
 case 1: Tuple<complex,16>::bluesteinWindow(fout,size);
 case 2: Tuple<complex,32>::bluesteinWindow(fout,size);
 case 3:Tuple<complex,64>::bluesteinWindow(fout,size);
 case 4: Tuple<complex,128>::bluesteinWindow(fout,size);
 case 5: Tuple<complex,256>::bluesteinWindow(fout,size);
 case 6: Tuple<complex,512>::bluesteinWindow(fout,size);
 case 7: Tuple<complex,1024>::bluesteinWindow(fout,size);
 case 8: Tuple<complex,2048>::bluesteinWindow(fout,size);
 case 9:Tuple<complex,4096>::bluesteinWindow(fout,size);
 case 10: Tuple<complex,8192>::bluesteinWindow(fout,size);
 case 11: Tuple<complex,16384>::bluesteinWindow(fout,size);
 case 12: Tuple<complex,32768>::bluesteinWindow(fout,size);
 case 13: Tuple<complex,65536>::bluesteinWindow(fout,size);
 default:
 Tuple<complex,8>::bluesteinWindow(fout,size);
 }
 }*/

LFHTEMP Vector<C> Vector<C>::fourierTransform(const Vector<mycomplex>& bluewindow) const{
	unsigned int t = (bluewindow.getSize())>>4;
	int i;
	for(i=0;t != 0;i++) t >>=1;

	switch(i){
		case 1: return(fourierTransform_routine<16>(bluewindow));
		case 2: return(fourierTransform_routine<32>(bluewindow));
		case 3: return(fourierTransform_routine<64>(bluewindow));
		case 4: return(fourierTransform_routine<128>(bluewindow));
		case 5: return(fourierTransform_routine<256>(bluewindow));
		case 6: return(fourierTransform_routine<512>(bluewindow));
		case 7: return(fourierTransform_routine<1024>(bluewindow));
		case 8: return(fourierTransform_routine<2048>(bluewindow));
		case 9: return(fourierTransform_routine<4096>(bluewindow));
		case 10: return(fourierTransform_routine<8192>(bluewindow));
		case 11: return(fourierTransform_routine<16384>(bluewindow));
		case 12: return(fourierTransform_routine<32768>(bluewindow));
		case 13: return(fourierTransform_routine<65536>(bluewindow));
		default:
			return(fourierTransform_routine<8>(bluewindow));

	}

}
LFHTEMP Vector<C> Vector<C>::invfourierTransform(const Vector<mycomplex>& bluewindow) const{
	unsigned int t = (bluewindow.getSize())>>4;
	int i;
	for(i=0;t != 0;i++) t >>=1;

	switch(i){
		case 1: return(invfourierTransform_routine<16>(bluewindow));
		case 2: return(invfourierTransform_routine<32>(bluewindow));
		case 3: return(invfourierTransform_routine<64>(bluewindow));
		case 4: return(invfourierTransform_routine<128>(bluewindow));
		case 5: return(invfourierTransform_routine<256>(bluewindow));
		case 6: return(invfourierTransform_routine<512>(bluewindow));
		case 7: return(invfourierTransform_routine<1024>(bluewindow));
		case 8: return(invfourierTransform_routine<2048>(bluewindow));
		case 9: return(invfourierTransform_routine<4096>(bluewindow));
		case 10: return(invfourierTransform_routine<8192>(bluewindow));
		case 11: return(invfourierTransform_routine<16384>(bluewindow));
		case 12: return(invfourierTransform_routine<32768>(bluewindow));
		case 13: return(invfourierTransform_routine<65536>(bluewindow));
		default:
			return(invfourierTransform_routine<8>(bluewindow));

	}
}



LFHTEMP template<int supersize> Vector<C> Vector<C>::fourierTransform_routine(const Vector<mycomplex>& bluewindow) const{

	Tuple<C, supersize> tmp;
	int i;
	int j=0;
	int s;
	int vsize = getSize();
	for(i=0;i<vsize;i++){
		double ang = -i * i * M_PI / vsize;
		mycomplex factor = mycomplex(cos(ang),sin(ang));
		tmp[j] = darray[i];
		tmp[j] *= factor;
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}
	for(;i<supersize;i++){
		tmp[j] = ExCo<C>::zero();
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}

	//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	tmp.fourierTransform_routine();
	//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

	Tuple<C, supersize> tmp2;
	//		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();

	//		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

	for(i=0;i<supersize;i++){
		tmp2[j]	= tmp[i];
		tmp2[j]	*= bluewindow[i];
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}

	//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	tmp2.invfourierTransform_routine();
	//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

	Vector<C> _out;
	_out.setSize(vsize);
	for(i=0;i<vsize;i++) {
		double ang = -i * i * M_PI / vsize;
		mycomplex factor = mycomplex(cos(ang),sin(ang));
		_out[i] = tmp2[i];
		_out[i] *= factor;

	}
	return(_out);
}


LFHTEMP template<int supersize> Vector<C> Vector<C>::invfourierTransform_routine(const Vector<mycomplex>& bluewindow) const{

	// implement Bluestein algorithm!

	// B windows : 2 {1.0+1.0i, 1.0-1.0i }
	// B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


	Tuple<C, supersize> tmp;
	int i;
	int j=0;
	int s;
	int vsize = getSize();

	for(i=0;i<vsize;i++){
		double ang = i * i * M_PI / vsize;
		mycomplex factor = mycomplex(cos(ang),sin(ang));
		tmp[j] = darray[i];
		tmp[j] *= factor;
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}
	for(;i<supersize;i++){
		tmp[j] = ExCo<C>::zero();
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}

	//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
	tmp.fourierTransform_routine();
	//				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

	Tuple<C, supersize> tmp2;
	//		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();

	//		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

	for(i=0;i<supersize;i++){
		tmp2[j]	= tmp[i];
		mycomplex factor =bluewindow[i];
		factor[1] = -factor[1];
		tmp2[j]	*= factor;
		for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
		j |= s;
	}

	//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
	tmp2.invfourierTransform_routine();
	//		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

	Vector<C> _out;
	_out.setSize(vsize);
	for(i=0;i<vsize;i++) {
		double ang = i * i * M_PI / vsize;
		mycomplex factor = mycomplex(cos(ang)/vsize,sin(ang)/vsize);
		_out[i] = tmp2[i];
		_out[i] *= factor;

	}
	return(_out);

}

/*
 LFHTEMP void Vector<C>::fourierTransform_semiroutine(const Vector<mycomplex>&  bluewindow) {

 unsigned int supersize = getSize();
 int i;
 int j=0;
 int s;
 int vsize = getSize();
 for(i=0;i<vsize;i++){
 double ang = -i * i * M_PI / vsize;
 complex factor = complex(cos(ang),sin(ang));
 tmp[j] = darray[i];
 tmp[j] *= factor;
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }
 for(;i<supersize;i++){
 tmp[j] = ExCo<C>::zero();
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }


 int step;
 int cur;
 int icur;
 complex multip[size];
 multip[0] = complex(1.0f,0.0f);

 for(step=1;step<size;step<<=1){
 for(icur=1;icur<(step<<1);icur++){
 double tmpdouble =  ((double)icur * M_PI) / step;
 multip[icur] = complex(cos(tmpdouble),sin(tmpdouble));
 }
 for(cur =0;cur<size;cur+= (step<< 1)){
 for(icur=0;icur<step;icur++){
 //					printf("%f\t%f\t",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
 C tmp = data[cur | step | icur];
 tmp *= multip[step | icur];
 //			printf("%f\t%f\t",tmp[0],tmp[1]);
 tmp += data[cur | icur];
 //		printf("%f\t%f\t",tmp[0],tmp[1]);
 data[cur | step | icur] *= multip[icur];
 data[cur | icur] += data[cur | step | icur];
 data[cur | step | icur] = tmp;
 //						printf("%f\t%f\n",((complex)data[cur | step | icur])[0],((complex)data[cur | step | icur])[1]);
 }
 }
 }



 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
 //	tmp.fourierTransform_routine();
 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

 Tuple<C, supersize> tmp2;
 //		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();

 //		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

 for(i=0;i<supersize;i++){
 tmp2[j]	= tmp[i];
 tmp2[j]	*= bluewindow[i];
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }

 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
 //	tmp2.invfourierTransform_routine();
 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

 Vector<C> _out;
 _out.setSize(vsize);
 for(i=0;i<vsize;i++) {
 double ang = -i * i * M_PI / vsize;
 complex factor = complex(cos(ang),sin(ang));
 _out[i] = tmp2[i];
 _out[i] *= factor;

 }



 }


 LFHTEMP void Vector<C>::invfourierTransform_routine(Vector<C> & fout, const complex * bluewindow, unsigned char mag) const{

 // implement Bluestein algorithm!

 // B windows : 2 {1.0+1.0i, 1.0-1.0i }
 // B windows : 4 {0.0+0.0i, 2.0-0.0i,(1 + i) * sqrt(2)/2, (1 + i) * sqrt(2)/2}


 Tuple<C, supersize> tmp;
 int i;
 int j=0;
 int s;
 int vsize = getSize();

 for(i=0;i<vsize;i++){
 double ang = i * i * M_PI / vsize;
 complex factor = complex(cos(ang),sin(ang));
 tmp[j] = darray[i];
 tmp[j] *= factor;
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }
 for(;i<supersize;i++){
 tmp[j] = ExCo<C>::zero();
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }

 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);
 tmp.fourierTransform_routine();
 //				for(i=0;i<supersize;i++) printf("1:%i: %f\t%f\n",i,tmp[i][0],tmp[i][1]);

 Tuple<C, supersize> tmp2;
 //		Tuple<C, supersize,Cflag> tmpb = bluesteinWindow<supersize>();

 //		for(i=0;i<supersize;i++) printf("b:%i: %f\t%f\n",i,tmpb[i][0],tmpb[i][1]);

 for(i=0;i<supersize;i++){
 tmp2[j]	= tmp[i];
 complex factor =bluewindow[i];
 factor[1] = -factor[1];
 tmp2[j]	*= factor;
 for(s = supersize >> 1;(j & s);s >>=1) j ^= s;
 j |= s;
 }

 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);
 tmp2.invfourierTransform_routine();
 //		for(i=0;i<supersize;i++) printf("2:%i: %f\t%f\n",i,tmp2[i][0],tmp2[i][1]);

 Vector<C> _out;
 _out.setSize(vsize);
 for(i=0;i<vsize;i++) {
 double ang = i * i * M_PI / vsize;
 complex factor = complex(cos(ang)/vsize,sin(ang)/vsize);
 _out[i] = tmp2[i];
 _out[i] *= factor;

 }
 return(_out);

 }*/


// uses interval queries D must be comparable to C !
LFHTEMP template<class D>
void Vector<C>::getIntersection(const IntervalSet<D> &query, Vector<C> &out) const{
	int node = asize & 0x7FFFFFF;
	stack<unsigned int> indexes;
	int imin = 0;
	int imax = node -1;
	int imid;
	SetComparison cmp, cmp2;
	cmp = query.compareInterval((*this)[imin],(*this)[imax]);
	if (cmp.areMonotonicDisjoint()) return;
	while(true){
		while(true){
			//printf("%i,%i\n",imin,imax); fflush(stdout);
			if (imax - imin  < 2){
				for(imid = imin; imid <= imax;imid++){
					cmp = query.compare((*this)[imid]);
					if (cmp.doesContain()) out.push_back((*this)[imid]);
				}
				break;
			}else{
				imid = (imin + imax) >> 1;
				cmp = query.compareInterval((*this)[imin], (*this)[imid]);
				cmp2 = query.compareInterval((*this)[imid+1], (*this)[imax]);
				if (cmp.areDisjoint()){
					if (cmp2.areDisjoint()) break;
					else{
						imin = imid+1;
					}
				}else{
					if (!(cmp2.areDisjoint())){
						indexes.push(imid+1);
						indexes.push(imax);
					}
					imax = imid;
				}
			}
		}

		if (indexes.empty()) break;
		imax = indexes.top();indexes.pop();
		imin = indexes.top();indexes.pop();
	}
}
LFHTEMP template<class I> void Vector<C>::operator() (Oper1<I> const & op){ // not a match
	int i = asize -1;
	for(;i>=0;i--) darray[i](op);
}
LFHTEMP void Vector<C>::operator() (Oper1< C> const & op){ // match
	int i = asize -1;
	for(;i>=0;i--) op(darray[i], darray[i]);
}
LFHTEMP template<class A_1, class A_2, class C_I> void Vector<C>::operator() (Oper2<A_1,A_2> const & op, const Vector<C_I> & a_2 ){ // not a match
	//setSize(_in.getSize());
	int i = asize -1;
	for(;i>=0;i--) (darray[i])(op,a_2.darray[i]);
}
LFHTEMP template<class C_I> void Vector<C>::operator() (Oper2<C,C_I> const & op, const Vector<C_I> & a_2){ // match
	//setSize(_in.getSize());
	int i = asize -1;
	for(;i>=0;i--) op(darray[i], a_2.darray[i]);
}
LFHTEMP template<class A_1, class A_2, class C_I> void Vector<C>::operator() (Oper2<A_1,A_2> const & op, Vector<C_I> & a_2 ){ // not a match
	//setSize(_in.getSize());
	int i = asize -1;
	for(;i>=0;i--) (darray[i])(op,a_2.darray[i]);
}
LFHTEMP template<class C_I> void Vector<C>::operator() (Oper2<C,C_I> const & op, Vector<C_I> & a_2){ // match
	//setSize(_in.getSize());
	int i = asize -1;
	for(;i>=0;i--) op(darray[i], a_2.darray[i]);
}
LFHTEMP template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> const & a_2, Vector<C_3> const & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
}
LFHTEMP template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> const & a_2, Vector<C_3> const & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
}
LFHTEMP template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> & a_2, Vector<C_3> const & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
}
LFHTEMP template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> & a_2, Vector<C_3> const & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
}
LFHTEMP template<class A_1, class A_2, class A_3, class C_2, class C_3> void Vector<C>::operator() (Oper3<A_1,A_2,A_3> const & op, Vector<C_2> & a_2, Vector<C_3> & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) (darray[i])(op,a_2.darray[i],a_3.darray[i]);
}
LFHTEMP template<class C_2, class C_3> void Vector<C>::operator() (Oper3<C,C_2,C_3> const & op, Vector<C_2> & a_2, Vector<C_3> & a_3 ){ // not a match
	int i = asize -1;
	for(;i>=0;i--) op(darray[i],a_2.darray[i],a_3.darray[i]);
}


/*
LFHTEMP
void Vector<C>::GP_Covariance(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), double noise_variance) const{

	unsigned int coor[2];
	coor[1] =coor[0] = getSize();
	f_out.setSizes(coor);

	double tmp;

	typedef class DataGrid<double,2>::KeyIterator itetype;
	itetype ite = f_out.getKeyIterator();

	if (ite.first()) do{

		if (ite()[0] == ite()[1]) {
			f_out(ite()) = 1.0f + noise_variance;
		}else{
			tmp = metric((*this)[ite()[0]], (*this)[ite()[1]]);
			f_out(ite()) = exp(-0.5f * tmp * tmp);
		}

	} while(ite.next());
}*/
/*
LFHTEMP  void Vector<C>::GP_Covariance_cross(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&), Vector<C>& query) const{

	unsigned int coor[2];
	coor[0] = getSize();
	coor[1] = query.getSize();

	f_out.setSizes(coor);

	double tmp;

	typedef class DataGrid<double,2>::KeyIterator itetype;
	itetype ite = f_out.getKeyIterator();

	if (ite.first()) do{
		tmp = metric((*this)[ite()[0]], query[ite()[1]]);
		f_out(ite()) = exp(-0.5f * tmp * tmp);
	} while(ite.next());



}

LFHTEMP
void Vector<C>::GP_fromCorrel(DataGrid<double, 2> &f_out, double (*metric)(const C&, const C&)) const{

	unsigned int coor[2];
	coor[1] =coor[0] = getSize();
	f_out.setSizes(coor);

	double tmp;

	typedef class DataGrid<double,2>::KeyIterator itetype;
	itetype ite = f_out.getKeyIterator();

	if (ite.first()) do{

		if (ite()[0] == ite()[1]) {
			f_out(ite()) = metric((*this)[ite()[0]], (*this)[ite()[1]]);
		}else{
			f_out(ite()) = metric((*this)[ite()[0]], (*this)[ite()[1]]);
		}

	} while(ite.next());



}

LFHTEMP void Vector<C>::GP_fromCorrel_cross(DataGrid<double, 2> &f_out, Vector<C>& query, double (*correl)(const C&, const C&)) const{

	unsigned int coor[2];
	coor[0] = getSize();
	coor[1] = query.getSize();

	f_out.setSizes(coor);

	double tmp;

	typedef class DataGrid<double,2>::KeyIterator itetype;
	itetype ite = f_out.getKeyIterator();

	if (ite.first()) do{
		tmp = metric((*this)[ite()[0]], query[ite()[1]]);
		f_out(ite()) = exp(-0.5f * tmp * tmp);
	} while(ite.next());



}

LFHTEMP void Vector<C>::GP_fromCorrel_complete(DataGrid<double, 2> &f_out, Vector<C>& query, double (*corr_metric)(const C&, const C&)) const{

	unsigned int coor[2];
	coor[0] = getSize();
	coor[1] = getSize() + query.getSize();

	f_out.setSizes(coor);
	int s = getSize();

	typedef class DataGrid<double,2>::KeyIterator itetype;
	itetype ite = f_out.getKeyIterator();

	if (ite.first()) do{
		if (ite()[1] < s){
			f_out(ite()) = corr_metric((*this)[ite()[0]], (*this)[ite()[1]]);
		}else{
			f_out(ite()) = corr_metric((*this)[ite()[0]], query[ite()[1] - s]);
		}
	} while(ite.next());
}*/

LFHTEMP Matrix<C> Vector<C>::mkouterprod() const{unsigned int s = asize & 0x7FFFFFFF; Matrix<C> fout(s, s);
	unsigned i,j;
	C* cur = fout.data;
	for(j=0;j<s;j++) for(i=0;i<s;i++) *(cur++) = darray[i] * darray[j];
	return(fout);
}
/*
VectorArray::VectorArray(): asize(0){}

LFHTEMP C& VectorArray::pushIn(int index){
	int32_t rsize = (asize & 0x7FFFFFFF);
	C* swp;
	int i,j;
	if (((asize + 1)^(asize))> rsize){
        if (asize & 0x80000000) asize &= 0x7FFFFFFf;
		else{
			swp = new C[((rsize+1) << 1)];
			if (rsize> 0){
                for(i=0;i<sizes.getSize();i++){
                    for(j=0;j<sizes[];j++){
                        ExOp::memmove(swp[j + (indexes[i] << 1)],darray[j + indexes[i]]);
                    }
                    indexes[i] <<= 1;
                }
                delete[](darray);
			}else{
                for(i=0;i<sizes.getSize();i++){
                    sizes[i] = 0;
                    indexes[i] = 0;
                }
			}
			darray = swp;
		}
	}

	i = (index == sizes.getSize() -1) ? rsize : indexes[index-1];
    asize++;
    sizes[index]++;



	if (indexes[index] + sizes[index] >= i) {
        i = (index == 0) ? 0 : indexes[index-1] + sizes[index-1];
        if (i<indexes[index]) {
            indexes[index]--;
            j= indexes[index];
        }else{

        }
	}else j = indexes[index] + sizes[index];
	return j;
}

void VectorArray::show(FILE* f, int level) const{
    fprintf(f, "VectorArray with size %i\n", asize & 0x7FFFFFFF);

}*/

#ifdef Rcpp_hpp
LFHTEMP template<int D> void Vector<C>::wrRCPP(Rcpp::Vector<D> &fout) const{
    fout = Rcpp::NumericVector(this->getSize());
    for(uint32_t i=0;i<this->getSize();i++) fout[i] = darray[i];
}

LFHTEMP template<int D> void Vector<C>::rdRCPP(const Rcpp::Vector<D> &fin){
    this->setSize(fin.length());
    for(uint32_t i=0;i<this->getSize();i++) darray[i] = Rcpp::as<C>(fin[i]);
}

#endif

LFHTEMP template<int ISC> Vector<C>::MetaIterator<ISC>::MetaIterator(const Vector<C>::MetaIterator<ISC>& trg): target(trg.target), i(trg.i), maxval(trg.maxval){
	subite = (trg.subite == NULL) ? NULL : new decltype(C().mkIterator())(*(trg.subite));
}
LFHTEMP template<int ISC> Vector<C>::MetaIterator<ISC>::operator bool (){
	if (target.asize == 0) return false;
	i=0;
	if (subite) delete(subite);
	do{
		subite = new decltype(C().mkIterator())(target.darray[i].mkIterator());
		if (*subite) return true;
		i++;delete(subite);
	}while(i< target.getSize());
	subite = NULL;
return false;}
LFHTEMP template<int ISC> bool Vector<C>::MetaIterator<ISC>::operator++(int){
	if ((*subite)++) return true;
	i++;delete(subite);
	while(i< target.getSize()){
		subite = new decltype(C().mkIterator())(target.darray[i].mkIterator());
		if (*subite) return true;
		i++;delete(subite);
	}
	subite = NULL;
return false;}
LFHTEMP template<int ISC> typename Vector<C>::template MetaIterator<ISC>& Vector<C>::MetaIterator<ISC>::operator=(const Vector<C>::MetaIterator<ISC>& ite){
	i = ite.i;
	maxval = ite.maxval;
	if (subite) delete(subite);
	subite = (ite.subite == NULL) ? NULL : new decltype(C().mkIterator())(*(ite.subite));
return *this;}

/*
LFHTEMP template<int ISC> Vector<C>::MetaIteratorIsaPtr<ISC>::operator bool (){
	if (target.asize == 0) return false;
	i=0;
	do{
		subite = new decltype(C()->mkIterator())(target.darray[i]->mkIterator());
		if (*subite) return true;
		i++;delete(subite);
	}while(i< target.getSize());
return false;}
LFHTEMP template<int ISC> bool Vector<C>::MetaIteratorIsaPtr<ISC>::operator++(int){
	if ((*subite)++) return true;
	i++;delete(subite);
	while(i< target.getSize()){
		subite = new decltype(C()->mkIterator())(target.darray[i]->mkIterator());
		if (*subite) return true;
		i++;delete(subite);
	}
return false;}

*/
#undef LFHTEMP
#define LFHTEMP template<class Key, class Data, class HashFunc>

LFHTEMP myHashmap<Key,Data,HashFunc>::myHashmap(const myHashmap<Key, Data, HashFunc>& other): heap(other.heap),hash_mag(other.hash_mag){
    if (hash_mag > 0){
        hash = new unsigned int[1 << other.hash_mag];
        LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc> constructor\n(size = %i)\n", (int)((1 << other.hash_mag)* sizeof(unsigned int)))
        memcpy(hash, other.hash, sizeof(unsigned int) << hash_mag);
    }
}
LFHTEMP myHashmap<Key,Data,HashFunc>::myHashmap(myHashmap<Key, Data, HashFunc>&& other): heap(std::move(other.heap)), hash(other.hash),hash_mag(other.hash_mag){
    other.hash_mag =0;
}
LFHTEMP myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::operator=(const myHashmap<Key, Data, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap = other.heap;
    hash_mag = other.hash_mag;
    if (hash_mag != 0){
        hash = new unsigned int[1 << hash_mag];
        LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc>::operator=\n(size = %i)\n",(int)(sizeof(unsigned int) << hash_mag))
        memcpy(hash, other.hash , sizeof(unsigned int) <<  hash_mag);
    }
return *this;}
LFHTEMP myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::operator=(myHashmap<Key, Data, HashFunc>&& other){
    heap = std::move(other.heap);
    hash = other.hash;
    hash_mag = other.hash_mag;
    other.hash_mag =0;
return *this;}
LFHTEMP myHashmap<Key,Data,HashFunc>::~myHashmap(){if (hash_mag) delete[](hash);}
LFHTEMP Data& myHashmap<Key,Data,HashFunc>::operator[](const Key & key){
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF){

        //printf("insert\n");
        //showLinks();
       	uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1u)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            //printf("up alloc in progress!\n"); fflush(stdout);
            heap.push_back();
            heap[j].k.k = key;
           // printf("rehash!\n"); fflush(stdout);
            rehash(hash_mag + 1);
        }else{
            heap.push_back();
            heap[j].k.k = key;
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].d = hash[i];
            hash[i] = j;
        }
        ExOp::toZero(heap[j].k.d);
    }
return heap[j].k.d;}
LFHTEMP KeyElem<Key, Data>& myHashmap<Key,Data,HashFunc>::getKeyElem(const Key &key){
 unsigned int j = find(key);
    if (j == 0xFFFFFFFF){

        //printf("insert\n");
        //showLinks();
       	uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1u)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            //printf("up alloc in progress!\n"); fflush(stdout);
            heap.push_back();
            heap[j].k.k = key;
           // printf("rehash!\n"); fflush(stdout);
            rehash(hash_mag + 1);
        }else{
            heap.push_back();
            heap[j].k.k = key;
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].d = hash[i];
            hash[i] = j;
        }
        ExOp::toZero(heap[j].k.d);
    }
return heap[j].k;}
LFHTEMP Data myHashmap<Key,Data,HashFunc>::operator[](const Key & key) const{
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF) {Data fout; ExOp::toZero(fout); return fout;}
return heap[j].k.d;}
LFHTEMP Data& myHashmap<Key,Data,HashFunc>::addEntry(const Key & key){
    uint32_t rsize = (heap.asize & 0x7FFFFFFF);
    uint32_t j = heap.size();
    if ((((heap.asize + 1u)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
		//printf("would up alloc...\n");
        heap.push_back();
        heap[j].k.k = key;
		//printf("would rehash...\n");
        rehash(hash_mag + 1);
    }else{
        heap.push_back();
        heap[j].k.k = key;
        uint32_t  i = hashpos( HashFunc::makeSeed(key) );
        heap[j].d = hash[i];
        hash[i] = j;
    }
return heap[j].k.d;}
LFHTEMP template<class Data_COMP> bool myHashmap<Key,Data,HashFunc>::remEntry(const Key & k, const Data_COMP& datacomp){
    uint32_t* ite = &hash[hashpos( HashFunc::makeSeed(k) )];
    while(*ite != 0xFFFFFFFF){
        if ((heap[*ite].k.k == k)&&(heap[*ite].k.d == datacomp)){
            uint32_t tmp = *ite;
            *ite = heap[*ite].d;
            swap_and_pop(tmp);
            return true;
        }
        ite = &(heap[*ite].d);
    }
return false;}
LFHTEMP myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::toMemmove(myHashmap<Key, Data, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap.toMemmove(other.heap);
    hash = other.hash; hash_mag = other.hash_mag; other.hash_mag = 0;
return *this;}
LFHTEMP myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::toMemswap(myHashmap<Key, Data, HashFunc>& other){
    heap.toMemswap(other.heap);
    union{
        unsigned int* swhash;
        unsigned char swhash_mag;
    };
    swhash = other.hash; other.hash = hash; hash = swhash;
    swhash_mag = other.hash_mag; other.hash_mag = hash_mag; hash_mag = swhash_mag;
return *this;}

LFHTEMP template<class D> myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::operator=(const myHashmap<Key, D, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap = other.heap;
    hash_mag = other.hash_mag;
    if (hash_mag != 0){
        hash = new unsigned int[1 << hash_mag];
        LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc>::operator=\n(size = %i)\n",(int)(sizeof(unsigned int) << hash_mag))
        memcpy(hash, other.hash , sizeof(unsigned int) <<  hash_mag);
    }
return *this;}
LFHTEMP template<class D> myHashmap<Key, Data, HashFunc>& myHashmap<Key,Data,HashFunc>::toConvert(const myHashmap<Key, D, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap.toConvert(other.heap);
    hash_mag = other.hash_mag;
    if (hash_mag != 0){
        hash = new unsigned int[1 << hash_mag];
        LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc> constructor\n(size = %i)\n", (int)((1 << hash_mag)* sizeof(unsigned int)))
        memcpy(hash, other.hash , sizeof(unsigned int) <<  hash_mag);
    }
return *this;}
LFHTEMP void myHashmap<Key,Data,HashFunc>::remake_links_routine(){
    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);
    unsigned int i;
    for(i=0;i< heap.getSize();i++){
        unsigned int j = hashpos( HashFunc::makeSeed(heap[i].k.k) );
        heap[i].d = hash[j];
        hash[j] = i;
    }
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::sort(){
    if (hash_mag == 0) return;
    Compareclass comp;
    heap.sort_comp(comp);
    remake_links_routine();
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::sort_decr(){
    if (hash_mag == 0) return;
    Compareclass_decr comp;
    heap.sort_comp(comp);
    remake_links_routine();
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::keysort(){
    if (hash_mag == 0) return;
    Compareclass_key comp;
    heap.sort_comp(comp);
    remake_links_routine();
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::keysort_decr(){
    if (hash_mag == 0) return;
    Compareclass_keydecr comp;
    heap.sort_comp(comp);
    remake_links_routine();
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::keysort_byteswap(){
    if (hash_mag == 0) return;
    Compareclass_key comp;
	int i;
	for(i=0;i<heap.getSize();i++) ExOp::bytereverse(heap[i].k.k);
    heap.sort_comp(comp);
	for(i=0;i<heap.getSize();i++) ExOp::bytereverse(heap[i].k.k);
    remake_links_routine();
}

/*
LFHTEMP  int  myHashmap<Key,Data,HashFunc>::Compareclass_keybyteswap::operator()(const KeyElem< KeyElem<Key, Data> , unsigned int >& a, const KeyElem< KeyElem<Key, Data> , unsigned int >&b)const{
	if (ExOp::isEQ(a.k.k, b.k.k)) return 0;
	Key as = a.k.k;
	Key bs = b.k.k;
	ExOp::bytereverse(as);
	ExOp::bytereverse(bs);
	return (ExOp::isGT(as, bs)) ? 1 : -1;
}*/
LFHTEMP template<class COMP> void myHashmap<Key,Data,HashFunc>::sortKeyElem(const COMP &comp){
    if (hash_mag == 0) return;
	heap.sort_comp(comp);
    remake_links_routine();
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::insert_at_position(unsigned int pos, const Key &key, const Data& data){
    unsigned int j = find(key);
    //printf("insert at pos\n");
    //showLinks();

    if (j != 0xFFFFFFFF){
       	fprintf(stderr,"entry already there! (myhashexception)\n"); LFH_exit(1);
    }else{
        uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1u)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            heap.insert_at_position(pos, pair< KeyElem<Key, Data> , unsigned int >( KeyElem<Key, Data>(key,data) , 0) );
            rehash(hash_mag + 1);
        }else{
            heap.insert_at_position(pos, pair< KeyElem<Key, Data> , unsigned int >( KeyElem<Key, Data>(key,data) , 0) );
            for(;pos<heap.getSize();pos++){
                unsigned int i =  hashpos( HashFunc::makeSeed(heap[pos].k.k) );
                heap[pos].d = hash[i];
                hash[i] = pos;
            }
        }
    }
}
LFHTEMP unsigned int myHashmap<Key,Data,HashFunc>::find(const Key & k )const{
    if (hash_mag == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( HashFunc::makeSeed(k) )];
    //printf("pos %i is ...%i, mag is %i gives %i\n", heap.getSize(), hashpos( HashFunc::makeSeed(k) ), hash_mag, j); fflush(stdout);
    /*if ((heap.getSize() == 64)&&(7 == hashpos(HashFunc::makeSeed(k)))){
        printf("breaktime\n");
    }*/
    for(;j != 0xFFFFFFFF;j = heap[j].d){
       // printf("%i and %i\n", j, 1 << hash_mag);
        if (k == heap[j].k.k) break;
    }
return(j);}
LFHTEMP void myHashmap<Key,Data,HashFunc>::showLinks(FILE* f) const{
    unsigned int k;
    fprintf(f, "Hash (Address %X), hash Range:\n", (int) this);
    if (hash_mag == 0) {fprintf(f, "Empty!\n"); return;}
    for(k =0 ;(k >> hash_mag)  == 0;k++) fprintf(f, "%X: %i\n", k, hash[k] );
    fprintf(f, "Vector Range:\n");
    for(k =0 ; k < heap.size();k++) fprintf(f,"%i: %i (hashpos: %X)\n", k, heap[k].d, hashpos( HashFunc::makeSeed(heap[k].k.k)) );
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::swap_and_pop(unsigned int ite){
    // swap and pop
   // printf("pop swap %i, ,%i\n", ite, heap.size());fflush(stdout);
    uint32_t rsize = (heap.asize & 0x7FFFFFFF);
    if (ite == heap.size()-1){
		if ((((heap.asize)^(heap.asize-1))> rsize-1)&&(( (heap.asize & 0x80000000) != 0)||(rsize == 1) )  ) {
			heap.pop_back();
			rehash((rsize == 1) ? 0 : hash_mag-1);
		}else{
			heap.pop_back();
		}
    }else{
        // printf("spop not last!\n"); fflush(stdout);

		if ((((heap.asize)^(heap.asize-1))> rsize-1)&&( (heap.asize & 0x80000000) != 0)  ) { // going to down_alloc
			heap.pop_swap(ite);
			rehash( hash_mag-1);
		} else{
			unsigned int i = heap.size()-1;
			unsigned int *j = hash + hashpos( HashFunc::makeSeed(heap[i].k.k) );
			while(*j != i) j = &(heap[*j].d);
			*j = ite;
			heap.pop_swap(ite);
		}
    }
}


LFHTEMP void myHashmap<Key,Data,HashFunc>::erase(const Key &k){
    //printf("erase\n");
    //showLinks();
    unsigned int i = hashpos( HashFunc::makeSeed(k) );
    unsigned int *j = hash + i;
    i = hash[i];
    for(;i != 0xFFFFFFFF;i = heap[i].d){
        if (k == heap[i].k.k) break;
        j = &(heap[i].d);
    }
    if (i == 0xFFFFFFFF) return;
    *j = heap[i].d;
    swap_and_pop(i);
}
LFHTEMP void myHashmap<Key,Data,HashFunc>::erase_from_iterator(unsigned int ite){
    //printf("erase from ite\n");
    //showLinks();

    if (ite == 0xFFFFFFFF) return;
    if (ite >= heap.getSize()) {fprintf(stderr, "tried to delete illegal iterator: %i\n", ite); LFH_exit(1);}
    unsigned int i = hashpos( HashFunc::makeSeed(heap[ite].k.k) );
    unsigned int *j = hash + i;

   // printf("%i and %i\n",*j, ite ); fflush(stdout);
    while(*j != ite) {
        j = &(heap[*j].d);
        if (*j == 0xFFFFFFFF) {fprintf(stderr, "tried to delete illegal iterator: %i\n", ite); LFH_exit(1);}
    }
    *j = heap[*j].d;


    // swap and pop
    //printf("spop in!\n"); fflush(stdout);

    swap_and_pop(ite);
  //  printf("spop out!\n"); fflush(stdout);
//showLinks();

}

LFHTEMP void myHashmap<Key,Data,HashFunc>::changeKey(unsigned int ite, const Key &nk){
    unsigned int i = hashpos( HashFunc::makeSeed(heap[ite].k.k) );
    unsigned int *j = hash + i;

   // printf("%i and %i\n",*j, ite ); fflush(stdout);
    while(*j != ite) {

        j = &(heap[*j].d);

        if (*j == 0xFFFFFFFF) {fprintf(stderr, "tried to delete illegal iterator: %i\n", ite); LFH_exit(1);}
    }
    *j = heap[*j].d;

    heap[ite].k.k = nk;
    i = hashpos( HashFunc::makeSeed(heap[ite].k.k) );

    heap[ite].d = hash[i];
    hash[i] = ite;
}

LFHTEMP unsigned int myHashmap<Key,Data,HashFunc>::hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}

LFHTEMP void myHashmap<Key,Data,HashFunc>::rehash(unsigned char _new_mag){
    //if (hash_mag == 0) printf("Hash (Address %X) init hashing!\n", (int)this);
    //else if (_new_mag == 0) printf("Hash (Address %X) flush hashing!\n", (int)this);
    //else printf("Hash (Address %X) init hashing!\n", (int)this);
    //showLinks();
    //fflush(stdout);
    //printf("rehash %i -> %i\n", hash_mag , _new_mag); fflush(stdout);

    if (hash_mag) delete[](hash);
    hash_mag = _new_mag;
 //   printf("rehashing to size %i\n", 1<<_new_mag); fflush(stdout);
    if (_new_mag == 0) {hash = nullptr; return;}

    hash = new unsigned int[1 << hash_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc>::rehash\n(size = %i)\n",(int)(sizeof(unsigned int) << hash_mag))


    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);

    unsigned int i;
    for(i=0;i< heap.size();i++){
        unsigned int j = hashpos( HashFunc::makeSeed(heap[i].k.k) );
        heap[i].d = hash[j];
        hash[j] = i;
    }
}


LFHTEMP void myHashmap<Key,Data,HashFunc>::swapEntries(uint32_t iteA, uint32_t iteB){
    unsigned int* ia = hash + hashpos( HashFunc::makeSeed(heap[iteA].k.k) );
    unsigned int* ib = hash + hashpos( HashFunc::makeSeed(heap[iteB].k.k) );
    if (ia == ib) heap[iteA].k.toMemswap(heap[iteB].k); // lucky, no need to change hashmap
    else{
        while(*ia != iteA) ia = &(heap[*ia].d);
        while(*ib != iteB) ib = &(heap[*ib].d);
        heap[iteA].toMemswap(heap[iteB]);
        *ia = iteB;
        *ib = iteA;
    }
}

LFHTEMP ERRCODE myHashmap<Key,Data,HashFunc>::save(FILE*f)const{ ERRCODE fout;
    uint32_t s = heap.getSize();
    fout = (fwrite(&s, sizeof(uint32_t),1,f) == 1) ? 0 : 1;
    for(unsigned int i=0;i<s;i++) fout |= ExOp::save(heap[i].k,f);
return fout;}

LFHTEMP ERRCODE myHashmap<Key,Data,HashFunc>::load(const uint8_t * &chunk, const uint8_t* const endpos){
	uint32_t s,i;
	if (endpos - chunk < sizeof(uint32_t)) return 1;
	s = *(uint32_t*) chunk; chunk += sizeof(uint32_t);
	if (s == 0) {this->toMemfree(); return 0;}
	heap.setSize(s);
	ERRCODE fout =0;
	for(i=0;i<s;i++) {
		fout |= ExOp::load(heap[i].k,chunk, endpos);
	}
	if (fout != 0) {heap.toMemfree(); return 1;}
    for(i=1;i<32;i++) if ((s >> i) == 0) break;
    this->rehash(i);
return 0;}

LFHTEMP ERRCODE myHashmap<Key,Data,HashFunc>::load(FILE*f){ ERRCODE fout;
    uint32_t s;
    fout = (fread(&s, sizeof(uint32_t),1,f) != 1) ? 1 : 0;
    if ((fout == 1)||(s == 0)) this->toMemfree();
    else{
    	heap.setSize(s);
		unsigned int i;
		for(i=0;i<s;i++) fout |= ExOp::load(heap[i].k,f);
		for(i=1;i<32;i++) if ((s >> i) == 0) break;
		this->rehash(i);
    }
return fout;}

LFHTEMP void myHashmap<Key,Data,HashFunc>::show(FILE*f, int level)const{
    if (hash_mag == 0) { fprintf(f, "Un-initialized HashMap, nbitems=%i\n", heap.size());}
    else{
        fprintf(f, "HashMap hashsize=%i, nbitems=%i\n", 1 << hash_mag, heap.size());
        unsigned int i;
        unsigned int ttt =0;
        for(i=0;(i>> hash_mag) == 0;i++){
               printf("%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].d) {printf("\t");ExOp::show(heap[j].k,f,2);ttt++;}
               printf("\n");
            }
        if (ttt != heap.size()) {printf("ill state! %i %i\n", heap.size(), ttt); LFH_exit(1);}
    }
}


#undef LFHTEMP
#define LFHTEMP template<class K, class D>

LFHTEMP AliasedHashmap<K,D>::AliasedHashmap(const AliasedHashmap<K,D>& o): asize(o.asize), hash_mask(o.hash_mask), emptypos(o.emptypos){
	if (asize != 0){
		main_hash = new uint32_t[hash_mask +1];
		memcpy(main_hash, o.main_hash,sizeof(uint32_t) * (hash_mask +1));
		heap = new pair< KeyElem< K, D>, Tuple<uint32_t, 2u> >[hash_mask +1];
		for(uint32_t i = 0 ; i <= hash_mask; i++ ) heap[i] = o.heap[i];
	}
}
LFHTEMP AliasedHashmap<K,D>::~AliasedHashmap(){
    if (asize) {delete[](main_hash); delete[](heap);}
}
LFHTEMP AliasedHashmap<K,D>& AliasedHashmap<K,D>::operator=(const AliasedHashmap<K,D>& o){
	if (asize) {delete[](main_hash); delete[](heap);}
	asize = o.asize;
	hash_mask = o.hash_mask;
	emptypos = o.emptypos;
	if (asize != 0){
		main_hash = new uint32_t[hash_mask +1];
		memcpy(main_hash, o.main_hash,sizeof(uint32_t) * (hash_mask +1));
		heap = new pair< KeyElem< K, D>, Tuple<uint32_t, 2u> >[hash_mask +1];
		for(uint32_t i = 0 ; i <= hash_mask; i++ ) heap[i] = o.heap[i];
	}
return *this;}
LFHTEMP AliasedHashmap<K,D>& AliasedHashmap<K,D>::operator=(AliasedHashmap<K,D>&& o){
	if (asize) {delete[](main_hash); delete[](heap);}
	asize = o.asize; hash_mask = o.hash_mask; heap = o.heap; main_hash = o.main_hash; emptypos = o.emptypos;
	o.asize=0;
return *this;}

LFHTEMP uint32_t AliasedHashmap<K,D>::prepare_slot_routine(uint32_t &alias, uint8_t protected_bits){
	if (heap[emptypos].second[1] == emptypos) this->rehash_routine((hash_mask << 1) | 1);
	uint32_t j = heap[emptypos].second[1];
	bool freetouse;
	if (protected_bits < 32) {
		alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		freetouse = (((alias ^ j) & hash_mask) >> (31 - protected_bits)) < 2;
		if (!freetouse){ // need to make sure it is not already used
			while(this->hasAlias(alias)) alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		}
	}else if ((freetouse = ((heap[alias & hash_mask].second[0]+1) >= 0x80000001))){
		j = alias & hash_mask;
		if ((alias & hash_mask) == emptypos) emptypos = heap[emptypos].second[1];
	}
	if (freetouse) { // check if alias can be selected to fit
		heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
		heap[heap[j].second[1] ].second[0] = heap[j].second[0];
		alias = (alias & (0xFFFFFFFF ^ hash_mask)) | j;
		heap[j].second[1] = alias | 0x80000000;
	}else{ // alias is not allowed to be set to the evident empty spot, need to check it the natural position is occupied or not.

		if ((heap[alias & hash_mask].second[0]+1) >= 0x80000001) { // 1-forced position is free :-)
			j = alias & hash_mask;
			//printf("it is... Free :-)! s[%i]=%X\n", j,  heap[j].second[0]+1 );

			heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
			heap[heap[j].second[1] ].second[0] = heap[j].second[0];
			if (emptypos == j) emptypos = heap[j].second[1];
			heap[j].second[1] = alias | 0x80000000;
		}else if ((heap[alias & hash_mask].second[1] & 0x80000000) != 0){ // 2-forced position is used by same hashID entry :-/
			//printf("it is... not me first! s[%i]=%X\n", alias & hash_mask,  heap[alias & hash_mask].second[1]);
			heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
			heap[heap[j].second[1] ].second[0] = heap[j].second[0];

			//printf("%X ^  %X -> %X\n", alias & hash_mask, j , (((alias ^ j) & hash_mask) >> (32 - protected_bits)) );
			//printf("I am forced! to %i but %i [A=%X] is there\n", j, heap[alias & hash_mask].first.k , heap[alias & hash_mask].second[1]);
			heap[j].second[1] = (alias & (0x7FFFFFFF ^ hash_mask)) | (heap[alias & hash_mask].second[1]  & hash_mask);
			heap[alias & hash_mask].second[1] = (heap[alias & hash_mask].second[1]& (0xFFFFFFFF ^ hash_mask)) | j;
		}else{ // 3- forced position is used by an another hashID   :-(
			// first move that foreign away
			//printf("it is... me First! get out! D: s[%i]=%X\n", alias & hash_mask,  heap[alias & hash_mask].second[1] );
			j = alias & hash_mask;
			uint32_t i = heap[j].second[1] & hash_mask;
			while((heap[i].second[1] & hash_mask) != j) i = heap[i].second[1] & hash_mask;
			heap[i].second[1] = (heap[i].second[1] & (0xFFFFFFFF ^ hash_mask)) | heap[emptypos].second[1];
			i = heap[emptypos].second[1];
			heap[emptypos].second[1] = heap[i].second[1];
			heap[heap[i].second[1] ].second[0] = emptypos | 0x80000000;

			ExOp::toMemmove(heap[i].first, heap[j].first);
			heap[i].second = heap[j].second;

			uint32_t k = ExOp::mkHashValue(heap[j].first.k) & hash_mask;
			if (main_hash[k] == j) main_hash[k] = i;
			else {
				k = main_hash[k];
				while(heap[k].second[0] != j) k = heap[k].second[0];
				heap[k].second[0] = i;
			}
			// ok, use it then it is all yours
			heap[j].second[1] = alias | 0x80000000;
		}
	}
return j;}

LFHTEMP D& AliasedHashmap<K,D>::createEntry(const K &key, uint32_t &alias, bool insert_at_end, uint8_t protected_bits){
	uint32_t i,j;
	if (asize == 0){
		if (protected_bits < 32) alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		hash_mask = 1;
		heap = new pair< KeyElem< K, D>, Tuple<uint32_t, 2u> >[2];
		main_hash = new uint32_t[2];
		j = ExOp::mkHashValue(key) & 1;
		main_hash[j] = alias & 1;
		main_hash[j^1] = 0xFFFFFFFF;

		j = alias & 1;
		heap[j].first.k = key;
		heap[j].second[0] = 0xFFFFFFFF;
		heap[j].second[1] = alias | 0x80000000;
		emptypos = j ^ 1;
		heap[emptypos].second[1] = emptypos;
		heap[emptypos].second[0] = 0x80000000 | emptypos;
		asize = 1;
	}else{
		// find empty position
		j = prepare_slot_routine(alias, protected_bits);
		heap[j].first.k = key;
		//(hash_mask >> (32 - protected_bits)) == 0
		i = ExOp::mkHashValue(ExOp::getIndex(heap[j].first)) & hash_mask;
		if ((!insert_at_end)||(main_hash[i] == 0xFFFFFFFF)){
			heap[j].second[0] = main_hash[i];
			main_hash[i] = j;
		}else{
			heap[j].second[0] = 0xFFFFFFFF;
			i = main_hash[i];
			while(heap[i].second[0] != 0xFFFFFFFF) i = heap[i].second[0];
			heap[i].second[0] = j;
		}
		asize++;
	}
return heap[j].first.d;}

LFHTEMP D& AliasedHashmap<K,D>::createEntryAfter(const U32_Alias &prev, uint32_t &alias, bool insert_at_end, uint8_t protected_bits){
	if (asize == 0) {printf("Error, trying to insert after unexisting node, is empty as well\n"); exit(1);}
	if (heap[emptypos].second[1] == emptypos) this->rehash_routine((hash_mask << 1) | 1);
	uint32_t ite = this->find(prev);
	if (ite == 0xFFFFFFFF) {printf("Error, trying to insert after unexisting node\n"); exit(1);}
	uint32_t j = prepare_slot_routine(alias, protected_bits);
	heap[j].first.k = heap[ite].first.k;
	heap[j].second[0] = heap[ite].second[0];
	heap[ite].second[0] = j;
	asize++;
return heap[j].first.d;}
LFHTEMP D& AliasedHashmap<K,D>::operator[](const K & k ){
    uint32_t ite = this->find(k);
    if (ite != 0xFFFFFFFF) return heap[ite].first.d;
	printf("Warning, forcefully accessing missing item with key:"); ExOp::show(k);
return heap[0].first.d;}
LFHTEMP D& AliasedHashmap<K,D>::operator[](const U32_Alias & alias ){
	uint32_t ite = this->find(alias);
    if (ite != 0xFFFFFFFF) return heap[ite].first.d;
	printf("Warning, forcefully accessing missing item with alias %X!\n", alias);
return heap[0].first.d;}
LFHTEMP uint32_t AliasedHashmap<K,D>::find(const K & k )const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = main_hash[ExOp::mkHashValue(k) & hash_mask];
    for(;j != 0xFFFFFFFF;j = heap[j].second[0]){
        if (k == heap[j].first.k) break;
    }
return(j);}
LFHTEMP unsigned int AliasedHashmap<K,D>::find(const U32_Alias &alias)const{
	if (asize == 0) return 0xFFFFFFFF;
	uint32_t j = alias.value & hash_mask;
	if ((heap[j].second[1] & 0x80000000) == 0) return 0xFFFFFFFF;
	if (((alias.value ^ heap[j].second[1]) & 0x7FFFFFFF) <= hash_mask)	return j;
	uint32_t c = heap[j].second[1] & hash_mask;
	while(((alias.value ^ heap[c].second[1]) & 0x7FFFFFFF) > hash_mask){
		c = heap[c].second[1] & hash_mask;
		if (j == c) return 0xFFFFFFFF;
	}
return c;}
LFHTEMP bool AliasedHashmap<K,D>::hasAlias(uint32_t alias)const{
	if (asize == 0) return false;
	uint32_t j = alias & hash_mask;
	if ((heap[j].second[1] & 0x80000000) == 0) return false;
	if (((alias ^ heap[j].second[1]) & 0x7FFFFFFF) <= hash_mask) return true;
	uint32_t c = heap[j].second[1] & hash_mask;
	while(((alias ^ heap[c].second[1]) & 0x7FFFFFFF) > hash_mask){
		c = heap[c].second[1] & hash_mask;
		if (j == c) return false;
	}
return true;}

LFHTEMP void AliasedHashmap<K,D>::erase_from_iterator(uint32_t ite, uint32_t *ite_to_maintain){
    if (ite == 0xFFFFFFFF) return;
    // step 1 unlink from main hash
    uint32_t j = ExOp::mkHashValue(heap[ite].first.k) & hash_mask;
    if (main_hash[j] == ite) main_hash[j] = heap[ite].second[0];
    else{
		j = main_hash[j];
		while(heap[j].second[0] != ite) j = heap[j].second[0];
		heap[j].second[0] = heap[ite].second[0];
    }

    // step 2: unlink the alias hash
    j = heap[ite].second[1] & hash_mask;
    if (j != ite){
		if ((heap[ite].second[1] & 0x80000000) != 0){ //Move next in line D:
			// update hash
			uint32_t i = ExOp::mkHashValue(heap[j].first.k) & hash_mask;
			if (main_hash[i] == j) main_hash[i] = ite;
			else{
				i = main_hash[i];
				while(heap[i].second[0] != j) i = heap[i].second[0];
				heap[i].second[0] = ite;
			}
			ExOp::toMemmove(heap[ite].first,heap[j].first);
			heap[ite].second[0] = heap[j].second[0];
			heap[ite].second[1] = heap[j].second[1] | 0x80000000;

			ite = j; // flag for removal
		}else{ // wont be missed :D
			while( (heap[j].second[1] & hash_mask) != ite) j = heap[j].second[1] & hash_mask;
			heap[j].second[1] = (heap[j].second[1] & (0xFFFFFFFF ^ hash_mask)) | (heap[ite].second[1] & hash_mask);
		}
    }// else: the last one... so not sad!

    // step 3: join ite to the empty bunch



    heap[ite].second[0] = emptypos | 0x80000000;
    heap[heap[emptypos].second[1]].second[0] = ite | 0x80000000;
    heap[ite].second[1] = heap[emptypos].second[1];
    heap[emptypos].second[1] = ite;

/*
    unsigned int *j = main_hash + hashpos( ExOp::mkHashValue(heap[ite].first.k.first) );
There is no repository for such regions, I would not venture to say I would know how to even define that. I think the best we can do is to list examples that work and abide by certain rules (stable and periodic solutions), since sharing partial solution on the forum will not yield working conduit unless it is a near miss, but listing catalysts needed by search algorithms would.

I am aware that we do have a chicken-and-egg problem here, but the idea is that one should only report conduit section that can be shown to be linkable to pre-existing conduits, otherwise is it is an oddity one should but in the "useless discovery" section or the newly created discussion threads, which still are looked at by all and might be fixed to a known conduit, if one presume a extended search would be successful, at least I cant see how to fix that O-to-B without periodic circuits. I did note your new P-to-G in the other thread, I believe there is no stable circuit for that, but only a periodic one (somewhere, sadly these results are all over the place) so I think one needs to showcase that it can be used. If the solution is only periodic, I think one could still put it in the stable conduit section, still it would be nice to show one can use it without a 2-glider collision (but even that would good if all else fails). Otherwise, I would be bragging about my marvelous syringe v2.0 discovery:

[code]x = 44, y = 36, rule = LifeHistory
8$22.2A$22.A$14.2A3.2A.A$14.A4.A.A$15.3A.A$17.A.2A2.5B$17.13B$18.2B.
10B$14.2A6.2C.7B$14.A3.B2.2B2CB2C4B$11.2A.A3.B3.C3B2C4B$11.A2.2A2.B8.
5B$12.2A4.A$14.5A$14.A4.2A$15.2A2.A$17.A.A$17.A.A.A2.A$18.2A.4A$20.A$
20.A.A$21.2A![/code]

That's unclear how to best to it, but I did try to organize the periodic solutions somehow, which was never done before, and the hope would be to make a collection too, but these should work (connect to other conduits, in and out) and low repeat time <300 as well.

    while(*j != ite) j = &(heap[*j].second.first);
    *j = heap[*j].second.first;


    j = main_hash + (hashpos( ExOp::mkHashValue(heap[ite].first.k.second) ) | (1 << hash_mag));

    while(*j != ite) j = &(heap[*j].second.second);
    *j = heap[*j].second.second;

    // swap and pop
    swap_and_pop(ite);*/


	asize--;
    if ((asize << 2) <= hash_mask) {
		this->rehash_routine( (asize == 0) ? 0 : (hash_mask >> 1));
	}
}
LFHTEMP uint32_t AliasedHashmap<K,D>::getAliasAt_routine(uint32_t pos) const{
	if ((heap[pos].second[1] & 0x80000000) != 0) return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (pos & hash_mask);
	uint32_t i,j;
	i = heap[pos].second[1] & hash_mask;
	while((heap[i].second[1] & 0x80000000) == 0){
		j = heap[i].second[1] & hash_mask;
		if ((heap[j].second[1] & 0x80000000) != 0) return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (j & hash_mask);
		i = heap[j].second[1] & hash_mask;
	}
return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (i & hash_mask);}
LFHTEMP AliasedHashmap<K,D>& AliasedHashmap<K,D>::toMemfree(){if (asize) {asize = 0; delete[](main_hash); delete[](heap); } return *this;}
LFHTEMP void AliasedHashmap<K,D>::rehash_routine(uint32_t _new_mask){
    if (_new_mask == 0) {hash_mask = 0; delete[](heap); delete[](main_hash); return; }
    //printf("allcatint %X was %X\n", _new_mask + 1, hash_mask); fflush(stdout);
    pair< KeyElem< K, D>, Tuple<uint32_t, 2u>  > *newheap = new pair< KeyElem< K, D>, Tuple<uint32_t, 2u> >[_new_mask+1];
    if (newheap == NULL) {printf("memory full! 32bit curse eh?\n"); exit(0);}
	// mark emptypos as such

	uint32_t i,j,k,l,m,n;

	for( i=0; i<= _new_mask ; i++) newheap[i].second[1] = 0;
	// insert non-empty
	Vector<KeyElem<uint32_t, uint32_t> > oldinds;

	Tuple<uint32_t> daorder; daorder.setSize(asize);
	for(j=0, i=0; i<= hash_mask ; i++){
		k = main_hash[i];
		while( k != 0xFFFFFFFF){
			daorder[j++] = k;
			k = heap[k].second[0];
		}
	}
	for( m=daorder.getSize(); m != 0 ;){
		i = daorder[--m];
		if ((((heap[i].second[1] ^ i) & hash_mask) != 0)&&( ((heap[heap[i].second[1] & hash_mask].second[1] ^ heap[i].second[1]) & hash_mask)  != 0)) { // replace with actual alias

			k = this->getAliasAt_routine(i)  | 0x80000000;
			j = heap[i].second[1] & hash_mask;
			heap[i].second[1] = k;
			k = (k & hash_mask) | 0x80000000;
			do{
				l = heap[j].second[1] & hash_mask;
				heap[j].second[1] = (heap[j].second[1] & (0x7FFFFFFF ^ hash_mask)) | k;
				if (l == i) break;
				j = heap[l].second[1] & hash_mask;
				heap[l].second[1] = (heap[l].second[1] & (0x7FFFFFFF ^ hash_mask)) | k;
			}while(j !=i);
			k = heap[i].second[1];
		}else k = heap[i].second[1];

		if (newheap[k & _new_mask].second[1] == 0) {
			newheap[k & _new_mask].second[1] = k | 0x80000000;
			k = k & _new_mask;
			ExOp::toMemmove(newheap[k].first, heap[i].first);
		}else oldinds.push_back( KeyElem<uint32_t, uint32_t>(k & _new_mask,i));
	}

	delete[](main_hash);
	main_hash = new uint32_t[_new_mask +1];
	if (main_hash == NULL) {printf("memory full! 32bit curse eh?\n"); exit(0);}
	memset(main_hash, '\xFF',sizeof(uint32_t) * (_new_mask +1));

	k = 0; n =0;
	//for(i=0; i <= _new_mask;i++) printf("%i:%X\n", i, newheap[i].second[1]);
	for( m=daorder.getSize(); m != 0 ;){
		i = heap[daorder[--m]].second[1] & _new_mask;
		//printf("cmp cmp %X vs %X\n", newheap[i].second[1], heap[daorder[m]].second[1] | 0x80000000);
		if ((newheap[i].second[1] ^ (heap[daorder[m]].second[1] | 0x80000000)) <= _new_mask){
			//printf("inserting that normal at %i\n", i);
			l = ExOp::mkHashValue(newheap[i].first.k) & _new_mask;
			newheap[i].second[0] = main_hash[l];
			main_hash[l] = i;
		}else{
			while(newheap[n].second[1] != 0) n++;
			//printf("inserting that extra at %i!\n", i);
			ExOp::toMemmove(newheap[n].first, heap[oldinds[k].d].first);
			l = ExOp::mkHashValue(newheap[n].first.k) & _new_mask;
			//printf("%i\t%i! dahask %i %X\n", oldinds[k].k , oldinds[k].d, l, heap[oldinds[k].d].second[1]);
			newheap[n].second[0] = main_hash[l];
			main_hash[l] = n;
			newheap[n].second[1] =  (heap[oldinds[k].d].second[1] & (0x7FFFFFFF ^ _new_mask)) |  (newheap[i].second[1] & _new_mask);
			newheap[i].second[1] =  (newheap[i].second[1] & (0xFFFFFFFF ^ _new_mask)) | n;
			k++; n++;
		}
	}

	while(newheap[n].second[1] != 0) n++;
	i = n++;
	emptypos = i;
	j = i;
	for(;n <= _new_mask; n++){
		if (newheap[n].second[1] != 0) continue;
		newheap[j].second[1] = n;
		newheap[n].second[0] = j | 0x80000000;
		j = n;
	}
	newheap[j].second[1] = emptypos;
	newheap[emptypos].second[0] = j | 0x80000000;

	delete[](heap); heap = newheap;
	hash_mask = _new_mask;
}
LFHTEMP ERRCODE AliasedHashmap<K,D>::save(FILE*f)const{
    ERRCODE fout = (fwrite(&asize,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
    if (asize == 0) return fout;
    Tuple<uint32_t> daorder; daorder.setSize(asize);
	uint32_t i,j,k;
	for(j=0, i=0; i<= hash_mask ; i++){
		k = main_hash[i];
		while( k != 0xFFFFFFFF){
			daorder[j++] = k;
			k = heap[k].second[0];
		}
	}
	for( i=daorder.getSize(); i-- != 0 ;){
		fout |= ExOp::save(heap[daorder[i]].first,f);
		j = getAliasAt_routine(daorder[i]);
		fout |= (fwrite(&j,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
	}
return(fout);}

LFHTEMP ERRCODE AliasedHashmap<K,D>::load(FILE*f){
    this->toMemfree();
    ERRCODE fout = (fread(&asize,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
    if ((fout)||(asize == 0)) return fout;
    hash_mask = asize + 1;
    ExCo<uint32_t>::toLeftFlood(hash_mask);

    Tuple<uint32_t> daorder; daorder.setSize(asize);
	main_hash = new uint32_t[hash_mask +1];
	memset(main_hash, '\xFF',sizeof(uint32_t) * (hash_mask +1));
	heap = new pair< KeyElem< K, D>, Tuple<uint32_t, 2u> >[hash_mask+1];
	uint32_t i,m,n,l;
	for( i=0; i<= hash_mask ; i++) heap[i].second[1] = 0;
	Vector<pair< KeyElem< K, D> , uint32_t> > oldinds;
	pair< KeyElem< K, D> , uint32_t> tmpelem;
	for(uint32_t i=0;i < asize ;i++){
		fout |= ExOp::load(tmpelem.first,f);
		fout |= (fread(&tmpelem.second,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
		if (heap[tmpelem.second & hash_mask].second[1] == 0) {
			heap[tmpelem.second & hash_mask].second[1] = tmpelem.second | 0x80000000;
			ExOp::toMemmove(heap[tmpelem.second & hash_mask].first, tmpelem.first);
			daorder[i] = (tmpelem.second & hash_mask) | 0x80000000;
		}else {daorder[i] = oldinds.getSize(); oldinds.push_back(tmpelem);}
	}
	n =0;
	for( m=0; m < asize ; m++){
		if ((daorder[m] & 0x80000000) != 0){
			i = daorder[m] & hash_mask;
			l = ExOp::mkHashValue(heap[i].first.k) & hash_mask;
			heap[i].second[0] = main_hash[l];
			main_hash[l] = i;
		}else{
			i = oldinds[daorder[m]].second & hash_mask;
			while(heap[n].second[1] != 0) n++;
			ExOp::toMemmove(heap[n].first, oldinds[daorder[m]].first);
			l = ExOp::mkHashValue(heap[n].first.k) & hash_mask;
			heap[n].second[0] = main_hash[l];
			main_hash[l] = n;
			heap[n].second[1] =  ( oldinds[daorder[m]].second  & (0x7FFFFFFF ^ hash_mask)) |  (heap[i].second[1] & hash_mask);
			heap[i].second[1] =  (heap[i].second[1] & (0xFFFFFFFF ^ hash_mask)) | n;
			n++;
		}
	}
	while(heap[n].second[1] != 0) n++;
	i = n++;
	emptypos = i;
	l = i;
	for(;n <= hash_mask; n++){
		if (heap[n].second[1] != 0) continue;
		heap[l].second[1] = n;
		heap[n].second[0] = l | 0x80000000;
		l = n;
	}
	heap[l].second[1] = emptypos;
	heap[emptypos].second[0] = l | 0x80000000;

	if (fout != 0) this->toMemfree();
return(fout);}
LFHTEMP void AliasedHashmap<K,D>::show(FILE*f, int level)const{
    if (asize == 0) { fprintf(f, "Un-initialized AliasedHashmap\n");}
    else{
        fprintf(f, "AliasedHashmap size=%i hashsize=%i \n",asize, hash_mask + 1);
        unsigned int i, ttt;
        for(i=0, ttt =0; i <= hash_mask;i++){
		   fprintf(f,"%i:",i);
		   for(unsigned int j= main_hash[i]; j != 0xFFFFFFFF;) {
				fprintf(f,"\t");ExOp::show(heap[j].first,f,2);
				fprintf(f,"\tAlias=%X,Pos=%i",getAliasAt_routine(j), j & hash_mask);
				ttt++;
				j = heap[j].second[0];
		   }
		   fprintf(f,"\n");
		}
		fprintf(f, "set of empty pos: %i", emptypos);
		for(i = heap[emptypos].second[1]; i != emptypos ; i = heap[i].second[1]) fprintf(f, ", %i", i);
		fprintf(f, "\n");
    }
}
LFHTEMP void AliasedHashmap<K,D>::test(FILE*f){
	Vector< KeyElem<K, Tuple< U32_Alias >   > > list; // list of entries
	KeyElem<K, Tuple< U32_Alias > > tmpelem;
	fprintf(f,"Dual hash map test, randomly insert and delete element, with objective size alternating between 8 and  2000 every 4096 iterations for million operations\n");
	unsigned int i;
	unsigned int r,j;
	for(i=0;i<10000;i++){
		ExOp::toRand(r); r >>= (i & 2047) ? 28: 20;
		//this->show();
		//if (asize != 0) for(int k=0; k <= hash_mask;k++) printf("%i:%u\t%e\t%X\t%X\n", k, heap[k].first.k, heap[k].first.d, heap[k].second[0], heap[k].second[1]);
		if (r >= list.size()){ // insert new entry
			j=rand() & 31;
			if (j & 16) j = 1;
			else if (j & 8) j = 2;
			else if (j == 0) j = 3;
			tmpelem.d.setSize(j);
			do{	ExOp::toRand(tmpelem.k);
			}while(find(tmpelem.k) != 0xFFFFFFFF);
			for(j =0; j < tmpelem.d.getSize();j++) {
				ExOp::toRand(tmpelem.d[tmpelem.d.getSize()-j-1].value);
				tmpelem.d[tmpelem.d.getSize()-j-1].value &= 0x7FFFFFFF;
				if (rand() & 4096) ExOp::toRand(this->createEntry(tmpelem.k, tmpelem.d[tmpelem.d.getSize()-j-1].value, 3));
				else {
					while(this->hasAlias(tmpelem.d[tmpelem.d.getSize()-j-1])){
						ExOp::toRand(tmpelem.d[tmpelem.d.getSize()-j-1].value);
						tmpelem.d[tmpelem.d.getSize()-j-1].value &= 0x7FFFFFFF;
					}
					ExOp::toRand(this->createEntryWithAlias(tmpelem.k, tmpelem.d[tmpelem.d.getSize()-j-1].value));
				}
				printf("%X\n", tmpelem.d[tmpelem.d.getSize()-j-1].value);
			}
			list.push_back(tmpelem);
			fprintf(f,"Insert :"); ExOp::show(tmpelem.k,f,1);fprintf(f,"\t%e\t%i elems\n", (*this)[tmpelem.k],tmpelem.d.getSize());
		}else{ // remove random entry
			ExOp::toRand(r);
			j = r % list.size();
			fprintf(f,"Delete key: %u with %i\n", list[j].k, list[j].d.getSize());
			for(int k=0; k < list[j].d.getSize();k++) {
				r = find(list[j].k);
				if ((r == 0xFFFFFFFF)||(((getAliasAt_routine(r) ^ list[j].d[k].value) & 0x7FFFFFFF)!=0)) {printf("finding %u, %X:\n", list[j].k , list[j].d[k].value); fprintf(f,"Could not find a used key-data by key, %X al=%X  exiting\n", r, r == 0xFFFFFFFF ? 0 : getAliasAt_routine(r));
					this->show();
					for(i=0; i <= hash_mask;i++) printf("%i:%u\t%e\t%X\t%X\n", i, heap[i].first.k, heap[i].first.d, heap[i].second[0], heap[i].second[1]);
					LFH_exit(1);
				}
				r = find(list[j].d[k]);
				if ((r == 0xFFFFFFFF)||(deref_key(r) != list[j].k)) {printf("finding %u, %X:\n", list[j].k , list[j].d[k].value); fprintf(f,"Could not find a used key-data by alias (needs to be first in line), %X  %u exiting\n", r, r ==0xFFFFFFFF ? 0 : deref_key(r));
					this->show();
					for(i=0; i <= hash_mask;i++) printf("%i:%u\t%e\t%X\t%X\n", i, heap[i].first.k, heap[i].first.d, heap[i].second[0], heap[i].second[1]);
					LFH_exit(1);}
				erase_from_iterator(r);
			}
			list.pop_swap(j);
		}
		if ((rand() & 255) == 0){
			FILE *f = fopen("dathing", "wb+");
			this->save(f);
			fclose(f);
			f = fopen("dathing", "rb+");
			this->load(f);
			fclose(f);
		}
	}
}

#undef LFHTEMP
#define LFHTEMP template<class K>



LFHTEMP AliasedHashmap<K,void>::AliasedHashmap(const AliasedHashmap<K,void>& o): asize(o.asize), hash_mask(o.hash_mask), emptypos(o.emptypos){
	if (asize != 0){
		main_hash = new uint32_t[hash_mask +1];
		memcpy(main_hash, o.main_hash,sizeof(uint32_t) * (hash_mask +1));
		heap = new pair< K, Tuple<uint32_t, 2u> >[hash_mask +1];
		for(uint32_t i = 0 ; i <= hash_mask; i++ ) heap[i] = o.heap[i];
	}
}
LFHTEMP AliasedHashmap<K,void>::~AliasedHashmap(){
    if (asize) {delete[](main_hash); delete[](heap);}
}
LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::operator=(const AliasedHashmap<K,void>& o){
	if (asize) {delete[](main_hash); delete[](heap);}
	asize = o.asize;
	hash_mask = o.hash_mask;
	emptypos = o.emptypos;
	if (asize != 0){
		main_hash = new uint32_t[hash_mask +1];
		memcpy(main_hash, o.main_hash,sizeof(uint32_t) * (hash_mask +1));
		heap = new pair< K, Tuple<uint32_t, 2u> >[hash_mask +1];
		for(uint32_t i = 0 ; i <= hash_mask; i++ ) heap[i] = o.heap[i];
	}
return *this;}
LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::operator=(AliasedHashmap<K,void>&& o){
	if (asize) {delete[](main_hash); delete[](heap);}
	asize = o.asize; hash_mask = o.hash_mask; heap = o.heap; main_hash = o.main_hash; emptypos = o.emptypos;
	o.asize=0;
return *this;}

LFHTEMP uint32_t AliasedHashmap<K,void>::prepare_slot_routine(uint32_t &alias, uint8_t protected_bits){
	if (heap[emptypos].second[1] == emptypos) this->rehash_routine((hash_mask << 1) | 1);
	uint32_t j = heap[emptypos].second[1];
	bool freetouse;
	if (protected_bits < 32) {
		alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		freetouse = (((alias ^ j) & hash_mask) >> (31 - protected_bits)) < 2;
		if (!freetouse){ // need to make sure it is not already used
			while(this->hasAlias(alias)) alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		}
	}else if ((freetouse = ((heap[alias & hash_mask].second[0]+1) >= 0x80000001))){
		j = alias & hash_mask;
		if ((alias & hash_mask) == emptypos) emptypos = heap[emptypos].second[1];
	}
	if (freetouse) { // check if alias can be selected to fit
		heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
		heap[heap[j].second[1] ].second[0] = heap[j].second[0];
		alias = (alias & (0xFFFFFFFF ^ hash_mask)) | j;
		heap[j].second[1] = alias | 0x80000000;
	}else{ // alias is not allowed to be set to the evident empty spot, need to check it the natural position is occupied or not.

		if ((heap[alias & hash_mask].second[0]+1) >= 0x80000001) { // 1-forced position is free :-)
			j = alias & hash_mask;
			//printf("it is... Free :-)! s[%i]=%X\n", j,  heap[j].second[0]+1 );

			heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
			heap[heap[j].second[1] ].second[0] = heap[j].second[0];
			if (emptypos == j) emptypos = heap[j].second[1];
			heap[j].second[1] = alias | 0x80000000;
		}else if ((heap[alias & hash_mask].second[1] & 0x80000000) != 0){ // 2-forced position is used by same hashID entry :-/
			//printf("it is... not me first! s[%i]=%X\n", alias & hash_mask,  heap[alias & hash_mask].second[1]);
			heap[heap[j].second[0] & hash_mask].second[1] = heap[j].second[1];
			heap[heap[j].second[1] ].second[0] = heap[j].second[0];

			//printf("%X ^  %X -> %X\n", alias & hash_mask, j , (((alias ^ j) & hash_mask) >> (32 - protected_bits)) );
			//printf("I am forced! to %i but %i [A=%X] is there\n", j, heap[alias & hash_mask].first.k , heap[alias & hash_mask].second[1]);
			heap[j].second[1] = (alias & (0x7FFFFFFF ^ hash_mask)) | (heap[alias & hash_mask].second[1]  & hash_mask);
			heap[alias & hash_mask].second[1] = (heap[alias & hash_mask].second[1]& (0xFFFFFFFF ^ hash_mask)) | j;
		}else{ // 3- forced position is used by an another hashID   :-(
			// first move that foreign away
			//printf("it is... me First! get out! D: s[%i]=%X\n", alias & hash_mask,  heap[alias & hash_mask].second[1] );
			j = alias & hash_mask;
			uint32_t i = heap[j].second[1] & hash_mask;
			while((heap[i].second[1] & hash_mask) != j) i = heap[i].second[1] & hash_mask;
			heap[i].second[1] = (heap[i].second[1] & (0xFFFFFFFF ^ hash_mask)) | heap[emptypos].second[1];
			i = heap[emptypos].second[1];
			heap[emptypos].second[1] = heap[i].second[1];
			heap[heap[i].second[1] ].second[0] = emptypos | 0x80000000;

			ExOp::toMemmove(heap[i].first, heap[j].first);
			heap[i].second = heap[j].second;

			uint32_t k = ExOp::mkHashValue(ExOp::getIndex(heap[j].first)) & hash_mask;
			if (main_hash[k] == j) main_hash[k] = i;
			else {
				k = main_hash[k];
				while(heap[k].second[0] != j) k = heap[k].second[0];
				heap[k].second[0] = i;
			}
			// ok, use it then it is all yours
			heap[j].second[1] = alias | 0x80000000;
		}
	}
return j;}

LFHTEMP K& AliasedHashmap<K,void>::createEntry(const K &key, uint32_t &alias, bool insert_at_end, uint8_t protected_bits){
	uint32_t i,j;
	if (asize == 0){
		if (protected_bits < 32) alias ^= (ExCo<uint32_t>::mkRand() >> protected_bits);
		hash_mask = 1;
		heap = new pair< K, Tuple<uint32_t, 2u> >[2];
		main_hash = new uint32_t[2];
		j = ExOp::mkHashValue(key) & 1;
		main_hash[j] = alias & 1;
		main_hash[j^1] = 0xFFFFFFFF;

		j = alias & 1;
		heap[j].first = key;
		heap[j].second[0] = 0xFFFFFFFF;
		heap[j].second[1] = alias | 0x80000000;
		emptypos = j ^ 1;
		heap[emptypos].second[1] = emptypos;
		heap[emptypos].second[0] = 0x80000000 | emptypos;
		asize = 1;
	}else{
		// find empty position
		j = prepare_slot_routine(alias, protected_bits);
		heap[j].first = key;
		//(hash_mask >> (32 - protected_bits)) == 0
		i = ExOp::mkHashValue(ExOp::getIndex(heap[j].first)) & hash_mask;
		if ((!insert_at_end)||(main_hash[i] == 0xFFFFFFFF)){
			heap[j].second[0] = main_hash[i];
			main_hash[i] = j;
		}else{
			heap[j].second[0] = 0xFFFFFFFF;
			i = main_hash[i];
			while(heap[i].second[0] != 0xFFFFFFFF) i = heap[i].second[0];
			heap[i].second[0] = j;
		}
		asize++;
	}
return heap[j].first;}
LFHTEMP K& AliasedHashmap<K,void>::createEntryAfter(const U32_Alias &prev, uint32_t &alias, bool insert_at_end, uint8_t protected_bits){
	if (asize == 0) {printf("Error, trying to insert after unexisting node, is empty as well\n"); exit(1);}
	if (heap[emptypos].second[1] == emptypos) this->rehash_routine((hash_mask << 1) | 1);
	uint32_t ite = this->find(prev);
	if (ite == 0xFFFFFFFF) {printf("Error, trying to insert after unexisting node\n"); exit(1);}
	uint32_t j = prepare_slot_routine(alias, protected_bits);
	heap[j].first = heap[ite].first;
	heap[j].second[0] = heap[ite].second[0];
	heap[ite].second[0] = j;
	asize++;
return heap[j].first;}



LFHTEMP K& AliasedHashmap<K,void>::operator[](const typename ExCo<K>::INDEX_TYPE & k ){
    uint32_t ite = this->find(k);
    if (ite != 0xFFFFFFFF) return heap[ite].first.d;
return heap[0].first.d;}
LFHTEMP K& AliasedHashmap<K,void>::operator[](const U32_Alias & alias ){
	uint32_t ite = this->find(alias);
    if (ite != 0xFFFFFFFF) return heap[ite].first.d;
return heap[0].first.d;}


LFHTEMP uint32_t AliasedHashmap<K,void>::find(const typename ExCo<K>::INDEX_TYPE & k )const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = main_hash[ExOp::mkHashValue(k) & hash_mask];
    for(;j != 0xFFFFFFFF;j = heap[j].second[0]){
        if (k == ExOp::getIndex(heap[j].first)) break;
    }
return(j);}

LFHTEMP unsigned int AliasedHashmap<K,void>::find(const U32_Alias &alias)const{
	if (asize == 0) return 0xFFFFFFFF;
	uint32_t j = alias.value & hash_mask;
	if ((heap[j].second[1] & 0x80000000) == 0) return 0xFFFFFFFF;
	if (((alias.value ^ heap[j].second[1]) & 0x7FFFFFFF) <= hash_mask)	return j;
	uint32_t c = heap[j].second[1] & hash_mask;
	while(((alias.value ^ heap[c].second[1]) & 0x7FFFFFFF) > hash_mask){
		c = heap[c].second[1] & hash_mask;
		if (j == c) return 0xFFFFFFFF;
	}
return c;}
LFHTEMP bool AliasedHashmap<K,void>::hasAlias(uint32_t alias)const{
	if (asize == 0) return false;
	uint32_t j = alias & hash_mask;
	if ((heap[j].second[1] & 0x80000000) == 0) return false;
	if (((alias ^ heap[j].second[1]) & 0x7FFFFFFF) <= hash_mask) return true;
	uint32_t c = heap[j].second[1] & hash_mask;
	while(((alias ^ heap[c].second[1]) & 0x7FFFFFFF) > hash_mask){
		c = heap[c].second[1] & hash_mask;
		if (j == c) return false;
	}
return true;}

LFHTEMP void AliasedHashmap<K,void>::erase_from_iterator(uint32_t ite){
    if (ite == 0xFFFFFFFF) return;
    // step 1 unlink from main hash
    uint32_t j = ExOp::mkHashValue(ExOp::getIndex(heap[ite].first)) & hash_mask;
    if (main_hash[j] == ite) main_hash[j] = heap[ite].second[0];
    else{
		j = main_hash[j];
		while(heap[j].second[0] != ite) j = heap[j].second[0];
		heap[j].second[0] = heap[ite].second[0];
    }

    // step 2: unlink the alias hash
    j = heap[ite].second[1] & hash_mask;
    if (j != ite){
		if ((heap[ite].second[1] & 0x80000000) != 0){ //Move next in line D:
			// update hash
			uint32_t i = ExOp::mkHashValue(ExOp::getIndex(heap[j].first)) & hash_mask;
			if (main_hash[i] == j) main_hash[i] = ite;
			else{
				i = main_hash[i];
				while(heap[i].second[0] != j) i = heap[i].second[0];
				heap[i].second[0] = ite;
			}
			ExOp::toMemmove(heap[ite].first,heap[j].first);
			heap[ite].second[0] = heap[j].second[0];
			heap[ite].second[1] = heap[j].second[1] | 0x80000000;

			ite = j; // flag for removal
		}else{ // wont be missed :D
			while( (heap[j].second[1] & hash_mask) != ite) j = heap[j].second[1] & hash_mask;
			heap[j].second[1] = (heap[j].second[1] & (0xFFFFFFFF ^ hash_mask)) | (heap[ite].second[1] & hash_mask);
		}
    }// else: the last one... so not sad!

    // step 3: join ite to the empty bunch



    heap[ite].second[0] = emptypos | 0x80000000;
    heap[heap[emptypos].second[1]].second[0] = ite | 0x80000000;
    heap[ite].second[1] = heap[emptypos].second[1];
    heap[emptypos].second[1] = ite;

/*
    unsigned int *j = main_hash + hashpos( ExOp::mkHashValue(heap[ite].first.k.first) );
There is no repository for such regions, I would not venture to say I would know how to even define that. I think the best we can do is to list examples that work and abide by certain rules (stable and periodic solutions), since sharing partial solution on the forum will not yield working conduit unless it is a near miss, but listing catalysts needed by search algorithms would.

I am aware that we do have a chicken-and-egg problem here, but the idea is that one should only report conduit section that can be shown to be linkable to pre-existing conduits, otherwise is it is an oddity one should but in the "useless discovery" section or the newly created discussion threads, which still are looked at by all and might be fixed to a known conduit, if one presume a extended search would be successful, at least I cant see how to fix that O-to-B without periodic circuits. I did note your new P-to-G in the other thread, I believe there is no stable circuit for that, but only a periodic one (somewhere, sadly these results are all over the place) so I think one needs to showcase that it can be used. If the solution is only periodic, I think one could still put it in the stable conduit section, still it would be nice to show one can use it without a 2-glider collision (but even that would good if all else fails). Otherwise, I would be bragging about my marvelous syringe v2.0 discovery:

[code]x = 44, y = 36, rule = LifeHistory
8$22.2A$22.A$14.2A3.2A.A$14.A4.A.A$15.3A.A$17.A.2A2.5B$17.13B$18.2B.
10B$14.2A6.2C.7B$14.A3.B2.2B2CB2C4B$11.2A.A3.B3.C3B2C4B$11.A2.2A2.B8.
5B$12.2A4.A$14.5A$14.A4.2A$15.2A2.A$17.A.A$17.A.A.A2.A$18.2A.4A$20.A$
20.A.A$21.2A![/code]

That's unclear how to best to it, but I did try to organize the periodic solutions somehow, which was never done before, and the hope would be to make a collection too, but these should work (connect to other conduits, in and out) and low repeat time <300 as well.

    while(*j != ite) j = &(heap[*j].second.first);
    *j = heap[*j].second.first;


    j = main_hash + (hashpos( ExOp::mkHashValue(heap[ite].first.k.second) ) | (1 << hash_mag));

    while(*j != ite) j = &(heap[*j].second.second);
    *j = heap[*j].second.second;

    // swap and pop
    swap_and_pop(ite);*/


	asize--;
    if ((asize << 2) <= hash_mask) {
		this->rehash_routine( (asize == 0) ? 0 : (hash_mask >> 1));
	}
}
LFHTEMP uint32_t AliasedHashmap<K,void>::getAliasAt_routine(uint32_t pos) const{
	if ((heap[pos].second[1] & 0x80000000) != 0) return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (pos & hash_mask);
	uint32_t i,j;
	i = heap[pos].second[1] & hash_mask;
	while((heap[i].second[1] & 0x80000000) == 0){
		j = heap[i].second[1] & hash_mask;
		if ((heap[j].second[1] & 0x80000000) != 0) return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (j & hash_mask);
		i = heap[j].second[1] & hash_mask;
	}
return (heap[pos].second[1] & (0x7FFFFFFF ^ hash_mask)) | (i & hash_mask);}
LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::toMemfree(){if (asize) {asize = 0; delete[](main_hash); delete[](heap); } return *this;}
LFHTEMP void AliasedHashmap<K,void>::rehash_routine(uint32_t _new_mask){
    if (_new_mask == 0) {hash_mask = 0; delete[](heap); delete[](main_hash); return; }

    pair< K, Tuple<uint32_t, 2u>  > *newheap = new pair< K, Tuple<uint32_t, 2u> >[_new_mask+1];
	// mark emptypos as such
	uint32_t i,j,k,l,m,n;

	for( i=0; i<= _new_mask ; i++) newheap[i].second[1] = 0;
	// insert non-empty
	Vector<KeyElem<uint32_t, uint32_t> > oldinds;

	Tuple<uint32_t> daorder; daorder.setSize(asize);
	for(j=0, i=0; i<= hash_mask ; i++){
		k = main_hash[i];
		while( k != 0xFFFFFFFF){
			daorder[j++] = k;
			k = heap[k].second[0];
		}
	}
	for( m=daorder.getSize(); m != 0 ;){
		i = daorder[--m];
		if ((((heap[i].second[1] ^ i) & hash_mask) != 0)&&( ((heap[heap[i].second[1] & hash_mask].second[1] ^ heap[i].second[1]) & hash_mask)  != 0)) { // replace with actual alias

			k = this->getAliasAt_routine(i)  | 0x80000000;
			j = heap[i].second[1] & hash_mask;
			heap[i].second[1] = k;
			k = (k & hash_mask) | 0x80000000;
			do{
				l = heap[j].second[1] & hash_mask;
				heap[j].second[1] = (heap[j].second[1] & (0x7FFFFFFF ^ hash_mask)) | k;
				if (l == i) break;
				j = heap[l].second[1] & hash_mask;
				heap[l].second[1] = (heap[l].second[1] & (0x7FFFFFFF ^ hash_mask)) | k;
			}while(j !=i);
			k = heap[i].second[1];
		}else k = heap[i].second[1];

		if (newheap[k & _new_mask].second[1] == 0) {
			newheap[k & _new_mask].second[1] = k | 0x80000000;
			k = k & _new_mask;
			ExOp::toMemmove(newheap[k].first, heap[i].first);
		}else oldinds.push_back( KeyElem<uint32_t, uint32_t>(k & _new_mask,i));
	}

	delete[](main_hash);
	main_hash = new uint32_t[_new_mask +1];
	printf("%p and %p\n",newheap, main_hash);
	memset(main_hash, '\xFF',sizeof(uint32_t) * (_new_mask +1));

	k = 0; n =0;
	//for(i=0; i <= _new_mask;i++) printf("%i:%X\n", i, newheap[i].second[1]);
	for( m=daorder.getSize(); m != 0 ;){
		i = heap[daorder[--m]].second[1] & _new_mask;
		//printf("cmp cmp %X vs %X\n", newheap[i].second[1], heap[daorder[m]].second[1] | 0x80000000);
		if ((newheap[i].second[1] ^ (heap[daorder[m]].second[1] | 0x80000000)) <= _new_mask){
			//printf("inserting that normal at %i\n", i);
			l = ExOp::mkHashValue(ExOp::getIndex(newheap[i].first)) & _new_mask;
			newheap[i].second[0] = main_hash[l];
			main_hash[l] = i;
		}else{
			while(newheap[n].second[1] != 0) n++;
			//printf("inserting that extra at %i!\n", i);
			l = oldinds[k].d; // l = newheap[i].second[1] & _new_mask;
			ExOp::toMemmove(newheap[n].first, heap[l].first);
			l = ExOp::mkHashValue(ExOp::getIndex(newheap[n].first)) & _new_mask;
			//printf("%i\t%i! dahask %i %X\n", oldinds[k].k , oldinds[k].d, l, heap[oldinds[k].d].second[1]);
			newheap[n].second[0] = main_hash[l];
			main_hash[l] = n;
			newheap[n].second[1] =  (heap[oldinds[k].d].second[1] & (0x7FFFFFFF ^ _new_mask)) |  (newheap[i].second[1] & _new_mask);
			newheap[i].second[1] =  (newheap[i].second[1] & (0xFFFFFFFF ^ _new_mask)) | n;
			k++; n++;
		}
	}

	while(newheap[n].second[1] != 0) n++;
	i = n++;
	emptypos = i;
	j = i;
	for(;n <= _new_mask; n++){
		if (newheap[n].second[1] != 0) continue;
		newheap[j].second[1] = n;
		newheap[n].second[0] = j | 0x80000000;
		j = n;
	}
	newheap[j].second[1] = emptypos;
	newheap[emptypos].second[0] = j | 0x80000000;

	delete[](heap); heap = newheap;
	hash_mask = _new_mask;
}
LFHTEMP ERRCODE AliasedHashmap<K,void>::save(FILE*f)const{
    ERRCODE fout = (fwrite(&asize,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
    if (asize == 0) return fout;
    Tuple<uint32_t> daorder; daorder.setSize(asize);
	uint32_t i,j,k;
	for(j=0, i=0; i<= hash_mask ; i++){
		k = main_hash[i];
		while( k != 0xFFFFFFFF){
			daorder[j++] = k;
			k = heap[k].second[0];
		}
	}
	for( i=daorder.getSize(); i-- != 0 ;){
		fout |= ExOp::save(heap[daorder[i]].first,f);
		j = getAliasAt_routine(daorder[i]);
		fout |= (fwrite(&j,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
	}
return(fout);}

LFHTEMP ERRCODE AliasedHashmap<K,void>::load(FILE*f){
    this->toMemfree();
    ERRCODE fout = (fread(&asize,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
    if ((fout)||(asize == 0)) return fout;
    hash_mask = asize + 1;
    ExCo<uint32_t>::toLeftFlood(hash_mask);

    Tuple<uint32_t> daorder; daorder.setSize(asize);
	main_hash = new uint32_t[hash_mask +1];
	memset(main_hash, '\xFF',sizeof(uint32_t) * (hash_mask +1));
	heap = new pair< K, Tuple<uint32_t, 2u> >[hash_mask+1];
	uint32_t i,m,n,l;
	for( i=0; i<= hash_mask ; i++) heap[i].second[1] = 0;
	Vector<pair< K , uint32_t> > oldinds;
	pair< K, uint32_t> tmpelem;
	for(uint32_t i=0;i < asize ;i++){
		fout |= ExOp::load(tmpelem.first,f);
		fout |= (fread(&tmpelem.second,sizeof(uint32_t),1,f) == 1) ? 0 : 1;
		if (heap[tmpelem.second & hash_mask].second[1] == 0) {
			heap[tmpelem.second & hash_mask].second[1] = tmpelem.second | 0x80000000;
			ExOp::toMemmove(heap[tmpelem.second & hash_mask].first, tmpelem.first);
			daorder[i] = (tmpelem.second & hash_mask) | 0x80000000;
		}else {daorder[i] = oldinds.getSize(); oldinds.push_back(tmpelem);}
	}
	n =0;
	for( m=0; m < asize ; m++){
		if ((daorder[m] & 0x80000000) != 0){
			i = daorder[m] & hash_mask;
			l = ExOp::mkHashValue(ExOp::getIndex(heap[i].first)) & hash_mask;
			heap[i].second[0] = main_hash[l];
			main_hash[l] = i;
		}else{
			i = oldinds[daorder[m]].second & hash_mask;
			while(heap[n].second[1] != 0) n++;
			ExOp::toMemmove(heap[n].first, oldinds[daorder[m]].first);
			l = ExOp::mkHashValue(ExOp::getIndex(heap[n].first)) & hash_mask;
			heap[n].second[0] = main_hash[l];
			main_hash[l] = n;
			heap[n].second[1] =  ( oldinds[daorder[m]].second  & (0x7FFFFFFF ^ hash_mask)) |  (heap[i].second[1] & hash_mask);
			heap[i].second[1] =  (heap[i].second[1] & (0xFFFFFFFF ^ hash_mask)) | n;
			n++;
		}
	}
	while(heap[n].second[1] != 0) n++;
	i = n++;
	emptypos = i;
	l = i;
	for(;n <= hash_mask; n++){
		if (heap[n].second[1] != 0) continue;
		heap[l].second[1] = n;
		heap[n].second[0] = l | 0x80000000;
		l = n;
	}
	heap[l].second[1] = emptypos;
	heap[emptypos].second[0] = l | 0x80000000;

	if (fout != 0) this->toMemfree();
return(fout);}
LFHTEMP void AliasedHashmap<K,void>::show(FILE*f, int level)const{
    if (asize == 0) { fprintf(f, "Un-initialized AliasedHashmap\n");}
    else{
        fprintf(f, "AliasedHashmap size=%i hashsize=%i \n",asize, hash_mask + 1);
        unsigned int i;
        unsigned int ttt =0;
        for(i=0; i <= hash_mask;i++){
		   fprintf(f,"%i:",i);
		   for(unsigned int j= main_hash[i]; j != 0xFFFFFFFF;j = heap[j].second[0]) {fprintf(f,"\t");ExOp::show(heap[j].first,f,2);fprintf(f,"\tAlias=%X",getAliasAt_routine(j)); ttt++;}
		   fprintf(f,"\n");
		}
		fprintf(f, "set of empty pos: %i", emptypos);
		for(i = heap[emptypos].second[1]; i != emptypos ; i = heap[i].second[1]) fprintf(f, ", %i", i);
		fprintf(f, "\n");
    }
}

LFHTEMP void AliasedHashmap<K,void>::test(FILE*f){
	Vector< KeyElem<K, Tuple< U32_Alias >   > > list; // list of entries
	KeyElem<K, Tuple< U32_Alias > > tmpelem;
	fprintf(f,"Dual hash map test, randomly insert and delete element, with objective size alternating between 8 and  2000 every 4096 iterations for million operations\n");
	unsigned int i;
	unsigned int r,j;
	for(i=0;i<10000;i++){
		ExOp::toRand(r); r >>= (i & 2047) ? 28: 20;
		//this->show();
		//if (asize != 0) for(int k=0; k <= hash_mask;k++) printf("%i:%u\t%e\t%X\t%X\n", k, heap[k].first.k, heap[k].first.d, heap[k].second[0], heap[k].second[1]);
		if (r >= list.size()){ // insert new entry
			j=rand() & 31;
			if (j & 16) j = 1;
			else if (j & 8) j = 2;
			else if (j == 0) j = 3;
			tmpelem.d.setSize(j);
			do{	ExOp::toRand(tmpelem.k);
			}while(find(tmpelem.k) != 0xFFFFFFFF);
			for(j =0; j < tmpelem.d.getSize();j++) {
				ExOp::toRand(tmpelem.d[tmpelem.d.getSize()-j-1].value);
				tmpelem.d[tmpelem.d.getSize()-j-1].value &= 0x7FFFFFFF;
				if (rand() & 4096) this->createEntry(tmpelem.k, tmpelem.d[tmpelem.d.getSize()-j-1].value, 3);
				else {
					while(this->hasAlias(tmpelem.d[tmpelem.d.getSize()-j-1])){
						ExOp::toRand(tmpelem.d[tmpelem.d.getSize()-j-1].value);
						tmpelem.d[tmpelem.d.getSize()-j-1].value &= 0x7FFFFFFF;
					}
					this->createEntryWithAlias(tmpelem.k, tmpelem.d[tmpelem.d.getSize()-j-1].value);
				}
				printf("%X\n", tmpelem.d[tmpelem.d.getSize()-j-1].value);
			}
			list.push_back(tmpelem);
			fprintf(f,"Insert :"); ExOp::show(tmpelem.k,f,1);fprintf(f,"\t%i elems\n",tmpelem.d.getSize());
		}else{ // remove random entry
			ExOp::toRand(r); j = r % list.size();
			fprintf(f,"Delete key: %u with %i\n", list[j].k, list[j].d.getSize());
			for(int k=0; k < list[j].d.getSize();k++) {
				r = find(list[j].k);
				if ((r == 0xFFFFFFFF)||(((getAliasAt_routine(r) ^ list[j].d[k].value) & 0x7FFFFFFF)!=0)) {printf("finding %u, %X:\n", list[j].k , list[j].d[k].value); fprintf(f,"Could not find a used key-data by key, %X al=%X  exiting\n", r, r == 0xFFFFFFFF ? 0 : getAliasAt_routine(r));
					this->show();
					for(i=0; i <= hash_mask;i++) printf("%i:%u\t%X\t%X\n", i, heap[i].first, heap[i].second[0], heap[i].second[1]);
					LFH_exit(1);
				}
				r = find(list[j].d[k]);
				if ((r == 0xFFFFFFFF)||(deref(r) != list[j].k)) {printf("finding %u, %X:\n", list[j].k , list[j].d[k].value); fprintf(f,"Could not find a used key-data by alias (needs to be first in line), %X  %u exiting\n", r, r ==0xFFFFFFFF ? 0 : deref(r));
					this->show();
					for(i=0; i <= hash_mask;i++) printf("%i:%u\t%X\t%X\n", i, heap[i].first, heap[i].second[0], heap[i].second[1]);
					LFH_exit(1);}
				erase_from_iterator(r);
			}
			list.pop_swap(j);
		}
		if ((rand() & 255) == 0){
			FILE *f = fopen("dathing", "wb+");
			this->save(f);
			fclose(f);
			f = fopen("dathing", "rb+");
			this->load(f);
			fclose(f);
		}
	}
}
LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::moveNextEntry(uint32_t cur, const K& new_key, int pos){
	uint32_t ncur = cur;
	uint32_t tmpnext = heap[ncur].second[0];
	while(true){
		if (tmpnext == 0xFFFFFFFF) return *this; // no next... nothing to move!
		if (heap[tmpnext].first == heap[cur].first)  break;
		ncur = tmpnext;
		tmpnext = heap[ncur].second[0];
	}

	heap[ncur].second[0] = heap[tmpnext].second[0];
	heap[tmpnext].first = new_key;
	// OK found that thing, move time!
	uint32_t i = ExOp::getIndex(new_key) & hash_mask;
	ncur = main_hash[i ];
	if (ncur == 0xFFFFFFFF) {
		// very easy, nothing here insert no matter pos
		 main_hash[i & hash_mask] = tmpnext;
		 heap[tmpnext].second[0] = ncur;
		 return *this;
	}
	while(pos < 0){
		while(ncur != 0xFFFFFFFF){
			if (heap[ncur].first == new_key) pos++;
			ncur = heap[ncur].second[0];
		}
		pos++;
		ncur = main_hash[i];
	}
	if (pos > 0){ // insert in front
		while(true){
			if (heap[ncur].first == new_key) if (pos-- == 1) break;
			ncur = heap[ncur].second[0];
		}
		heap[tmpnext].second[0] = heap[ncur].second[0];
		heap[ncur].second[0] = tmpnext;
	}else{ // insert in back
		 main_hash[i] = tmpnext;
		 heap[tmpnext].second[0] = ncur;
	}
return *this;}
LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::changeKeyForNextAfter(uint32_t ite, const K &new_key, uint32_t nbitem){
	if (nbitem == 0) return *this;
	uint32_t i = ExOp::getIndex(new_key);
	//this->show();
	//if (asize != 0) for(int kss=0; kss <= hash_mask;kss++) printf("%i:%u\t%X\t%X\n", kss, heap[kss].first, heap[kss].second[0], heap[kss].second[1]);
	if (((ExOp::getIndex(heap[ite].first) ^ i) & hash_mask) == 0){ // same hash position... super lucky?
		//printf("so lucky...!\n");
		for(i = heap[ite].second[0]; nbitem > 0;i = heap[i].second[0]){
			if (heap[ite].first == heap[i].first) {heap[i].first = ite; nbitem--;}
		}
	}else{
		uint32_t j = ite;
		//printf("%X == %X at %i right?\n",getAliasAt_routine(j), new_key );
		while(heap[ite].first != heap[heap[j].second[0]].first) j = heap[j].second[0];
		uint32_t curnext = main_hash[i & hash_mask];
		main_hash[i & hash_mask] = heap[j].second[0];
		//printf("main hash for %i, got to set H[%i]=%i first (used to be %X) \n", new_key, i & hash_mask, heap[j].second[0], curnext);
		i = j;
		//printf("set key %i -> %i\n", j, heap[j].second[0]);
		j = heap[j].second[0];
		heap[i].second[0] = heap[j].second[0];
		//printf("set n[%i]=%i\n", i, heap[j].second[0]);
		heap[j].first = new_key;
		//printf("set key %i at %i\n", new_key, j);
		for(nbitem--;nbitem > 0;){
			if (heap[ite].first == heap[ heap[i].second[0] ].first) {
				heap[j].second[0] = heap[i].second[0];
				j = heap[j].second[0];
				heap[i].second[0] = heap[j].second[0];
				heap[j].first = new_key;
				nbitem--;
			}else i = heap[i].second[0];
		}   // N->2->3->4->6->danewhash
		heap[j].second[0] = curnext;
	}
	//this->show();
	//if (asize != 0) for(int kss=0; kss <= hash_mask;kss++) printf("%i:%u\t%X\t%X\n", kss, heap[kss].first, heap[kss].second[0], heap[kss].second[1]);
return *this;}

//LFHTEMP AliasedHashmap<K,void>& AliasedHashmap<K,void>::changeKeyForFirst(const K &key, const K &new_key, uint32_t nbitem){return *this;}


//  int >
/*LFHTEMP template<int ISCONST> typename std::enable_if<ISCONST == 0, typename AliasedHashmap<K,void>::SlotIterator<ISCONST >& >::type
AliasedHashmap<K,void>::SlotIterator<ISCONST>::changeKeyForNext(typename ExCo<K>::INDEX_TYPE newkey, uint32_t nbnexts){

return *this;}*/

#undef LFHTEMP
#define LFHTEMP template<class K, class D >


LFHTEMP StructuredHash<K,D>::StructuredHash(){
	links.createEntryWithAlias(0,0);
}

LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::insertElem(const K& key, const D& ndata){
	//printf("Insert unrooted:\n");
	uint32_t dalias;
	links.createEntry(0,dalias,true);
	data.createEntryWithAlias(key,dalias) = ndata;
return *this;}

LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::insertElemUnder(const K& key, const D& ndata, const K& parent){
	//printf("Insert:\n");
	uint32_t dalias=0;
	uint32_t paralias = data.find(parent);
	if (paralias == 0xFFFFFFFF) {printf("warning: Could not find parent, inserting as unrooted"); paralias =0;}
	else paralias = data.deref_alias(paralias);
	links.createEntry(paralias,dalias,true,1);
	data.createEntryWithAlias(key,dalias) = ndata;
return *this;}

LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::insertElemAfter(const K& key, const D& ndata, const K& prev){
	//printf("Insert:\n");
	uint32_t dalias=0;
	uint32_t prealias = data.find(prev);
	if (prealias == 0xFFFFFFFF) {printf("warning: Could not find parent, inserting as unrooted"); prealias =0;}
	else prealias = data.deref_alias(prealias);
	links.createEntryAfter(U32_Alias(prealias),dalias,true,1);
	data.createEntryWithAlias(key,dalias) = ndata;
return *this;}

LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::insertParentNodeAt(const K& key, const D& ndata, U32_Alias first, uint32_t nbsub){
	uint32_t lcur = links.find(first);
	uint32_t pcur = links.main_hash[ links.heap[lcur].first & links.hash_mask];
	uint32_t nalias =0;
	if (pcur == lcur){
		links.createEntry(links.heap[lcur].first,nalias,false,1);
		// very first...
	}else{
		while(links.heap[pcur].second[0] != lcur) pcur = links.heap[pcur].second[0];
		links.createEntryAfter( links.getAliasAt_routine(pcur), nalias,false,1);
	}
	data.createEntryWithAlias(key,nalias) = ndata;
	if (nbsub == 0) return *this;
	// steal next, put under (some might mismatch!)
	links.changeKeyForNextAfter(links.find(U32_Alias(nalias)),nalias,nbsub);
return *this;}

LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::setChild(const K& child_key, const K& parent_key){
	// cannot assume anything! hurray!


return *this;}
LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::setNext(const K& left_key, const K& right_key){

return *this;}
LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::addChild(const K& child_key, const K& parent_key){

return *this;}
LFHTEMP StructuredHash<K,D>& StructuredHash<K,D>::addNext(const K& left_key, const K& right_key){

return *this;}

LFHTEMP ERRCODE StructuredHash<K,D>::save(FILE*f)const{
	ERRCODE fout = voidparents.save(f);
	fout |= links.save(f);
	fout |= data.save(f);
return fout;}
LFHTEMP ERRCODE StructuredHash<K,D>::load(FILE*f){
	ERRCODE fout = voidparents.load(f);
	fout |= links.load(f);
	fout |= data.load(f);
return fout;}

LFHTEMP void StructuredHash<K,D>::show(FILE*f, int level) const{
	Vector<uint32_t> todo;


	char funbuf[256];
	memset(funbuf,' ', sizeof(char) * 256);
	funbuf[255] = '\0';

	if (auto ite = links.mkIterator(0)) do{
		if (ite() == 0) continue;
		todo.push_back(0);
		todo.push_back(ite());
	}while(ite++);
	todo.reverse();
	fprintf(f, "Tree structure:\n");
	uint32_t nbspace, subite;
	U32_Alias daalias;
	uint32_t c;
	while(todo.getSize() != 0){
		nbspace = todo.last(); todo.pop_back();
		funbuf[nbspace] = '\0';
		fprintf(f, "%s", funbuf);
		daalias = todo.last(); todo.pop_back();
		if ((subite = data.find(daalias)) != 0xFFFFFFFF) {
			ExOp::show( data.deref_key(subite), f,1);
			fprintf(f, ":");
			ExOp::show( data.deref(subite), f,1);
			fprintf(f, " {Alias=%X}\n", daalias.value);
			c = todo.getSize();
			if (auto ite = links.mkIterator(daalias)) do{
				todo.push_back(nbspace +2);
				todo.push_back(ite());
			}while(ite++);
			todo.reverseLast(todo.getSize() - c);
		}else {printf("missing entry!\n"); exit(1);}
		funbuf[nbspace] = ' ';
	}
}

LFHTEMP template<int IC> typename std::enable_if<IC == 0, uint32_t>::type StructuredHash<K,D>::TreeIterator<IC>::moveNextForward(uint32_t nbpos){
	uint32_t ncur;
	uint32_t tmpnext = cur;
	uint32_t tmpnext2;
	//for(int kss=0; kss <= target.links.hash_mask;kss++) printf("%i:%u\t%X\t%X\n", kss, target.links.heap[kss].first, target.links.heap[kss].second[0], target.links.heap[kss].second[1]);
	do{
		ncur =  target.links.heap[tmpnext].second[0];
		//printf("vis %i->%i\t %X vs %X\n", tmpnext, ncur, target.links.heap[ncur].first, target.links.heap[cur].first);
		if (target.links.heap[ncur].first == target.links.heap[cur].first) break;
		tmpnext = ncur;
	}while(true);
	//printf("found da pre ! %i\n", tmpnext);
	while(true){
		tmpnext2 =  target.links.heap[ncur].second[0];
		//printf("vis %i->%i\t %X vs %X\n", ncur, tmpnext2, target.links.heap[tmpnext2].first, target.links.heap[cur].first);
		if (target.links.heap[tmpnext2].first == target.links.heap[cur].first) {
			if (nbpos-- == 1) break;
		}
		ncur = tmpnext2;
	}
	tmpnext2 = ncur;
	//printf("found da prenext ! %i\n", tmpnext2);
	ncur = target.links.heap[tmpnext2].second[0];
	tmpnext2 = target.links.heap[tmpnext2].second[0] = target.links.heap[tmpnext].second[0];
	target.links.heap[tmpnext].second[0] = ncur;
	//for(int kss=0; kss <= target.links.hash_mask;kss++) printf("%i:%u\t%X\t%X\n", kss, target.links.heap[kss].first, target.links.heap[kss].second[0], target.links.heap[kss].second[1]);
	//printf("%i\t%i\t%i\t%i\t%i\n", ncur, tmpnext, tmpnext2, target.links.heap[ncur].second[0], target.links.heap[tmpnext2].second[0]);
	tmpnext = target.links.heap[ncur].second[0];
	target.links.heap[ncur].second[0] = target.links.heap[tmpnext2].second[0];
	target.links.heap[tmpnext2].second[0] = tmpnext;
	//for(int kss=0; kss <= target.links.hash_mask;kss++) printf("%i:%u\t%X\t%X\n", kss, target.links.heap[kss].first, target.links.heap[kss].second[0], target.links.heap[kss].second[1]);
return 0;}

LFHTEMP template<int IC> typename std::enable_if<IC == 0, uint32_t>::type StructuredHash<K,D>::TreeIterator<IC>::removeNext(){
	uint32_t ncur = cur;
	uint32_t tmpnext = target.links.heap[ncur].second[0];
	while(true){
		if (tmpnext == 0xFFFFFFFF) return 1;
		if (target.links.heap[tmpnext].first == target.links.heap[cur].first)  break;
		ncur = tmpnext;
		tmpnext = target.links.heap[ncur].second[0];
	}
	// OK found that thing, delete time!

	// erase_from_iterator
	U32_Alias daalias2(target.links.deref_alias(tmpnext));
	target.data.erase_from_iterator(target.data.find(daalias2));

	// update d ptr
	U32_Alias daalias(target.links.deref_alias(cur));
	dcur = target.data.find(daalias);

	target.links.erase_from_iterator(tmpnext);
	cur = target.links.find(daalias);


return 0;}
LFHTEMP template<int IC> typename std::enable_if<IC == 0, uint32_t>::type StructuredHash<K,D>::TreeIterator<IC>::moveNextUnder(int pos){
	U32_Alias daalias(target.links.deref_alias(cur));
	target.links.moveNextEntry(cur, daalias.value, pos);
	cur = target.links.find(daalias);
return 0;}
#undef LFHTEMP
#define LFHTEMP template<class Key, class HashFunc >

LFHTEMP myHashmap<Key,void,HashFunc>::myHashmap(const myHashmap<Key,void, HashFunc>& other): heap(other.heap),hash_mag(other.hash_mag),hash(new unsigned int[1 << other.hash_mag]){
    LFH_NICE_ALLOCERROR(hash,"Could not allocate for myHashmap<Key,Data,HashFunc>::rehash\n(size = %i)\n",(int)(sizeof(unsigned int) << other.hash_mag))
    memcpy(hash, other.hash, sizeof(unsigned int) << hash_mag);
}
LFHTEMP myHashmap<Key,void,HashFunc>::myHashmap(myHashmap<Key,void, HashFunc>&& other): heap(std::move(other.heap)),hash_mag(other.hash_mag), hash(other.hash){
    other.hash_mag =0;
}
LFHTEMP myHashmap<Key, void, HashFunc>& myHashmap<Key,void,HashFunc>::operator=(const myHashmap<Key, void, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap = other.heap;
    hash_mag = other.hash_mag;
    if (hash_mag != 0){
        hash = new unsigned int[1 << hash_mag];
        LFH_NICE_ALLOCERROR(hash, "Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))
        memcpy(hash, other.hash , sizeof(unsigned int) <<  hash_mag);
    }
return *this;}
LFHTEMP myHashmap<Key, void, HashFunc>& myHashmap<Key, void,HashFunc>::operator=(myHashmap<Key,void, HashFunc>&& other){
    heap = std::move(other.heap);
    hash = other.hash;
    hash_mag = other.hash_mag;
    other.hash_mag =0;
return *this;}

LFHTEMP myHashmap<Key,void,HashFunc>::~myHashmap(){
    if (hash_mag) delete[](hash);
}
LFHTEMP unsigned int myHashmap<Key,void,HashFunc>::hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}
LFHTEMP Key& myHashmap<Key,void,HashFunc>::operator[](const typename ExCo<Key>::INDEX_TYPE & key){
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF){
        //printf("insert\n");
        //showLinks();
       	uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            //printf("up alloc in progress!\n"); fflush(stdout);
            heap.push_back();
            ExOp::getIndex(heap[j].first) = key;
           // printf("rehash!\n"); fflush(stdout);
            rehash(hash_mag + 1);
        }else{
            heap.push_back();
            ExOp::getIndex(heap[j].first) = key;
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].second = hash[i];
            hash[i] = j;
        }
    }
    return heap[j].first;
}
LFHTEMP Key myHashmap<Key,void,HashFunc>::operator[](const typename ExCo<Key>::INDEX_TYPE & key) const{
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF) {Key fout; ExOp::toZero(fout); return fout;}
    return heap[j].first;
}
LFHTEMP void myHashmap<Key,void,HashFunc>::addEntry(const Key & key){


    unsigned int j = find(ExOp::getIndex(key));
    if (j == 0xFFFFFFFF){

        //printf("insert\n");
        //showLinks();
       	uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            //printf("up alloc in progress!\n"); fflush(stdout);
            heap.push_back();
            heap[j].first = key;
           // printf("rehash!\n"); fflush(stdout);
            rehash(hash_mag + 1);
        }else{
            heap.push_back();
            heap[j].first = key;
            unsigned int i = hashpos( HashFunc::makeSeed(ExOp::getIndex(key)) );
            heap[j].second = hash[i];
            hash[i] = j;
        }
    }
}
LFHTEMP typename ExCo<Key>::INDEX_TYPE myHashmap<Key,void,HashFunc>::deref_key(unsigned int ite) const{return ExOp::getIndex(heap[ite].first);}
LFHTEMP Key& myHashmap<Key,void,HashFunc>::createEntry(const typename ExCo<Key>::INDEX_TYPE &key){
	unsigned int j = find(key);
    if (j == 0xFFFFFFFF){

        //printf("insert\n");
        //showLinks();
       	uint32_t rsize = (heap.asize & 0x7FFFFFFF);
        j = heap.size();
        if ((((heap.asize + 1)^(heap.asize))> rsize)&&((heap.asize & 0x80000000) == 0)) { // going to up_alloc
            //printf("up alloc in progress!\n"); fflush(stdout);
            heap.push_back();
            heap[j].first.getIndex() = key;

           // printf("rehash!\n"); fflush(stdout);
            rehash(hash_mag + 1);
        }else{
            heap.push_back();
            heap[j].first.getIndex() = key;
            unsigned int i = hashpos( HashFunc::makeSeed(key) );
            heap[j].second = hash[i];
            hash[i] = j;
        }
    }
	return heap[j].first;

}

LFHTEMP template<class F> myHashmap<Key, void, HashFunc>& myHashmap<Key,void,HashFunc>::toConvert(F& f, ITERABLE_DEF(F) ){
    this->toMemfree();
    if (auto ite = f.getIterator()) do{
        this->addEntry(*ite);
    }while(ite++);
return(*this);}

LFHTEMP void myHashmap<Key,void,HashFunc>::showLinks(FILE* f) const{
    unsigned int k;
    fprintf(f, "Hash Range:\n");
    for(k =0 ;(k >> hash_mag)  == 0;k++) fprintf(f, "%X: %i\n", k, hash[k] );
    fprintf(f, "Vector Range:\n");
    for(k =0 ; k < heap.size();k++) fprintf(f,"%i: %i (hashpos: %X)\n", k, heap[k].second, hashpos( HashFunc::makeSeed(heap[k].first.k)) );
}
LFHTEMP void myHashmap<Key,void,HashFunc>::swap_and_pop(unsigned int ite){
    // swap and pop
   // printf("pop swap %i, ,%i\n", ite, heap.size());fflush(stdout);
    uint32_t rsize = (heap.asize & 0x7FFFFFFF);
    if (ite == heap.size()-1){
        //    printf("spop last!\n"); fflush(stdout);

    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&(( (heap.asize & 0x80000000) != 0)||(rsize == 1) )  ) {
        heap.pop_back();rehash((rsize == 1) ? 0 : hash_mag-1);
    }else{
        heap.pop_back();
    }

    }else{

    if ((((heap.asize)^(heap.asize-1))> rsize-1)&&( (heap.asize & 0x80000000) != 0)  ) { // going to down_alloc
           heap.pop_swap(ite);
        rehash( hash_mag-1);
    } else{
        //showLinks();
        unsigned int i = heap.size()-1;
        unsigned int *j = hash + hashpos( HashFunc::makeSeed(ExOp::getIndex(heap[i].first)) );

        //printf("fix! %i %i %i\n", *j, i , ite); fflush(stdout);
        while(*j != i) j = &(heap[*j].second);
        *j = ite;
       heap.pop_swap(ite);
    }
    }
  //  ExOp::show(heap);
}
LFHTEMP unsigned int myHashmap<Key,void,HashFunc>::find(const typename ExCo<Key>::INDEX_TYPE & k)const{
    if (hash_mag == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( HashFunc::makeSeed(k) )];
    for(;j != 0xFFFFFFFF;j = heap[j].second){
        if (k == ExOp::getIndex(heap[j].first) ) break;
    }
    return(j);
}
LFHTEMP myHashmap<Key, void, HashFunc>& myHashmap<Key,void,HashFunc>::toMemmove(myHashmap<Key, void, HashFunc>& other){
    if (hash_mag != 0) delete[](hash);
    heap.toMemmove(other.heap);
    hash = other.hash; hash_mag = other.hash_mag; other.hash_mag = 0;
return *this;}
LFHTEMP myHashmap<Key, void, HashFunc>& myHashmap<Key,void,HashFunc>::toMemswap(myHashmap<Key, void, HashFunc>& other){
    heap.toMemswap(other.heap);
    union{
        unsigned int* swhash;
        unsigned char swhash_mag;
    };
    swhash = other.hash; other.hash = hash; hash = swhash;
    swhash_mag = other.hash_mag; other.hash_mag = hash_mag; hash_mag = swhash_mag;
return *this;}

LFHTEMP void myHashmap<Key,void,HashFunc>::erase(const Key &k){
    unsigned int i = hashpos( HashFunc::makeSeed(k) );
    unsigned int *j = hash + i;
    i = hash[i];
    for(;i != 0xFFFFFFFF;i = heap[i].second){

        if (k == heap[i].first.k) break;
        j = &(heap[i].second);
        }

    if (i == 0xFFFFFFFF) return;


    *j = heap[i].second;
    swap_and_pop(i);
}
LFHTEMP void myHashmap<Key,void,HashFunc>::removeEntry(const typename ExCo<Key>::INDEX_TYPE &k){
return this->erase_from_iterator(this->find(k));}
LFHTEMP void myHashmap<Key,void,HashFunc>::erase_from_iterator(unsigned int ite){
    if (ite == 0xFFFFFFFF) return;
    unsigned int i = hashpos( HashFunc::makeSeed(ExOp::getIndex(heap[ite].first)) );
    unsigned int *j = hash + i;

    while(*j != ite) j = &(heap[*j].second);
    *j = heap[ite].second;

    // swap and pop
    swap_and_pop(ite);
}
LFHTEMP void myHashmap<Key,void,HashFunc>::rehash(unsigned char _new_mag){
    if (hash_mag) delete[](hash);
    hash_mag = _new_mag;
 //   printf("rehashing to size %i\n", 1<<_new_mag); fflush(stdout);
    if (_new_mag == 0) {hash = nullptr; return;}
    hash = new unsigned int[1 << hash_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))
    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);

    unsigned int i;
   for(i=0;i< heap.size();i++){
        unsigned int j = hashpos( HashFunc::makeSeed(ExOp::getIndex(heap[i].first)) );
        heap[i].second = hash[j];
        hash[j] = i;
    }
}
LFHTEMP ERRCODE myHashmap<Key,void,HashFunc>::save(FILE*f)const{ ERRCODE fout;
    unsigned int s = heap.getSize();
    fout = (fwrite(&s,sizeof(unsigned int),1,f) == 1) ? 0 : 1;
    for(unsigned int i=0;i<s;i++) fout |= ExOp::save(heap[i].first,f);
    return fout;
}

LFHTEMP ERRCODE myHashmap<Key,void,HashFunc>::load(FILE*f){
    unsigned int s;
    if (fread(&s, sizeof(unsigned int),1 ,f) != 1) return 1;
    heap.setSize(s);
    unsigned int i;
    for(i=0;i<s;i++) ExOp::load(heap[i].first,f);
    for(i=1;i<32;i++) if ((s >> i) == 0) break;
    this->rehash(i);
	return 0;
}

LFHTEMP void myHashmap<Key,void,HashFunc>::show(FILE*f, int level)const{
    if (hash_mag == 0) { fprintf(f, "Un-initialized HashMap\n");}
    else{
        fprintf(f, "HashMap hashsize=%i, nbitems=%i\n", 1 << hash_mag, heap.size());
        unsigned int i;
        unsigned int ttt =0;
        for(i=0;(i>> hash_mag) == 0;i++){
               printf("%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].second) {printf("\t");ExOp::show(heap[j].first,f,2);ttt++;}
               printf("\n");
            }
        if (ttt != heap.size()) {printf("ill state! %i %i\n", heap.size(), ttt); LFH_exit(1);}
    }
}

#undef LFHTEMP
#define LFHTEMP template<class K,class H>

LFHTEMP unsigned int myHashmap<K,Anything,H>::hashpos(unsigned int seed) const{
return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}
LFHTEMP void myHashmap<K,Anything,H>::swap_and_pop(unsigned int ite){  // moves entry in vector, into an unlinked location
}
LFHTEMP uint32_t myHashmap<K,Anything,H>::find(const K &) const{
return 0xFFFFFFFF;}
LFHTEMP myHashmap<K,Anything,H>& myHashmap<K,Anything,H>::toMemmove(myHashmap<K,Anything,H>& other){
    if (hash_mag != 0) delete[](hash);
    heap.toMemmove(other.heap);
    data.toMemmove(other.data);
    hash = other.hash; hash_mag = other.hash_mag; other.hash_mag = 0;
return *this;}
LFHTEMP myHashmap<K,Anything,H>& myHashmap<K,Anything,H>::toMemswap(myHashmap<K,Anything,H>& other){
    heap.toMemswap(other.heap);
    data.toMemmove(other.data);
    union{
        unsigned int* swhash;
        unsigned char swhash_mag;
    };
    swhash = other.hash; other.hash = hash; hash = swhash;
    swhash_mag = other.hash_mag; other.hash_mag = hash_mag; hash_mag = swhash_mag;
return *this;}
LFHTEMP myHashmap<K,Anything,H>& myHashmap<K,Anything,H>::operator=(const myHashmap<K,Anything,H>& other){
}
LFHTEMP AnythingConst myHashmap<K,Anything,H>::operator[](const K &key) const{
    uint32_t ite = this->find(key);
return data[ite];}
LFHTEMP Anything myHashmap<K,Anything,H>::operator[](const K &key){
    uint32_t ite = this->find(key);
return data[ite];}
LFHTEMP Anything myHashmap<K,Anything,H>::addEntry(const K &key){
    data.push_back();
return data[data.getSize()-1];}
LFHTEMP void myHashmap<K,Anything,H>::erase(const K &){
}
LFHTEMP void myHashmap<K,Anything,H>::rmEntry(const K &){
}
LFHTEMP void myHashmap<K,Anything,H>::erase_from_iterator(unsigned int ite){
}
LFHTEMP void myHashmap<K,Anything,H>::rehash(unsigned char _new_mag){
    if (hash_mag) delete[](hash);
    hash_mag = _new_mag;
    if (_new_mag == 0) return;
    hash = new unsigned int[1 << hash_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))
    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);
    unsigned int i;
    for(i=0;i< heap.size();i++){
        unsigned int j = hashpos( H::makeSeed(ExOp::getIndex(heap[i].first)) );
        heap[i].second = hash[j];
        hash[j] = i;
    }
}
LFHTEMP ERRCODE myHashmap<K,Anything,H>::save(FILE*f)const{ ERRCODE fout;
    uint32_t s = heap.getSize();
    fout = (fwrite(&s,sizeof(unsigned int),1,f) == 1) ? 0 : 1;
    for(unsigned int i=0;i<s;i++) fout |= ExOp::save(heap[i].first,f);
    return fout | data.save(f);
}
LFHTEMP ERRCODE myHashmap<K,Anything,H>::load(FILE*f){
    uint32_t s;
    if (fread(&s, sizeof(uint32_t),1 ,f) != 1) return 1;
    heap.setSize(s);
    unsigned int i;
    for(i=0;i<s;i++) ExOp::load(heap[i].first,f);
    for(i=1;i<32;i++) if ((s >> i) == 0) break;
    this->rehash(i);
	return data.load(f);
}


#undef LFHTEMP
#define LFHTEMP template<class K, class D, class C, class F, class CF>

LFHTEMP CategoryHashmap<K,D,C,F,CF>::~CategoryHashmap(){
    if (asize) delete[](hash);
}

LFHTEMP unsigned int CategoryHashmap<K,D,C,F,CF>::hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}


LFHTEMP void CategoryHashmap<K,D,C,F,CF>::swap_and_pop(unsigned int ite){
   //TODO
}


LFHTEMP uint32_t CategoryHashmap<K,D,C,F,CF>::getNewSlotForCategory_routine(const C& cat) {
    uint32_t fout;
    uint32_t j = categories.find(cat);
    unsigned int b;
    if (j == 0xFFFFFFFF){
        if (categories.getSize() == 0) fout=nb_categoryless;
        else{
            b=categories.getSize()-1;
            fout = categories.deref(b).first + categories.deref(b).second;
        }
        categories[cat] = pair<unsigned int, unsigned int >(fout,1);
    }else{
        b = (j == 0) ? nb_categoryless : categories.deref(j-1).first + categories.deref(j-1).second;
        if (b == categories.deref(j).first)  fout = enlargeback_routine(j);
        else fout = (--categories.deref(j).first);
        categories.deref(j).second++;
    }
    return fout;
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::move_up_routine(unsigned int from, unsigned int unto){

	ExOp::toMemmove(heap[unto], heap[from]);
	unsigned int tmp = hashpos( F::makeSeed(ExOp::getIndex(heap[unto].first.k)) );
	if (hash[tmp] == from) hash[tmp] = unto;
	else{
		tmp = hash[tmp];
		while(heap[tmp].second != from) tmp = heap[tmp].second;
		heap[tmp].second = unto;
	}
}
LFHTEMP unsigned int CategoryHashmap<K,D,C,F,CF>::find(const K &key) const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( F::makeSeed(key) )];
    for(;j != 0xFFFFFFFF;j = heap[j].second){
        if (key == heap[j].first.k) break;
    }
    return(j);
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::remove_links_routine(uint32_t target, uint32_t catite){
    unsigned int *j = &(hash[hashpos( F::makeSeed(heap[target].first.k) )]);
    while((*j) != target) j = &heap[*j].second;
    (*j) = heap[target].second;
    if (catite == 0xFFFFFFFF) {nb_categoryless--; this->memmove_routine(target, nb_categoryless);}
    else{categories.deref(catite).second--; this->memmove_routine(target, categories.deref(catite).first + categories.deref(catite).second);}
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::memmove_routine(uint32_t target, uint32_t source){
    if (target == source) return;
    printf("called memmove %i->%i\n", source,target);
    unsigned int *j = &(hash[hashpos( F::makeSeed(heap[source].first.k) )]);
    while((*j) != source) j = &heap[*j].second;
    *j = target;
    heap[target].second = heap[source].second;
    heap[target].first.toMemmove(heap[source].first);
}



LFHTEMP unsigned int CategoryHashmap<K,D,C,F,CF>::enlargeback_routine(unsigned int catite){
    uint32_t curc =catite;
    while(curc+1 < categories.getSize()){
        if (categories.deref(curc).first + categories.deref(curc).second < categories.deref(curc+1).first) break;
        curc++;
    }
    if (categories.deref(curc).first + categories.deref(curc).second == (1u << hash_mag)){
        for(curc = catite-1;curc>0; curc--){
            if (categories.deref(curc-1).first + categories.deref(curc-1).second < categories.deref(curc+1).first) break;
        }
        while(curc < catite){
            categories.deref(curc).first--;
            memmove_routine(categories.deref(curc).first, categories.deref(curc).first + categories.deref(curc).second);
            curc++;
        }
        categories.deref(curc).first--;
        return categories.deref(curc).first;
    }
    while (curc > catite){
        memmove_routine(categories.deref(curc).first + categories.deref(curc).second, categories.deref(curc).first);
        categories.deref(curc).first++;
        curc--;
    }
    return categories.deref(catite).first + categories.deref(catite).second;
    /*     unsigned int fout = categories.deref(catite).first + categories.deref(catite).second;

	unsigned int to_move = asize & 0x7FFFFFF;
	while(fout != to_move){
		for(unsigned int i=0;;i++){
			if (i>=categories.getSize()){printf("error in CategoryHash\n");LFH_exit(1);}
			if (categories.deref(i).first + categories.deref(i).second == to_move) {
				move_up_routine(categories.deref(i).first,to_move);
				to_move = categories.deref(i).first++;
				break;
			}
		}
	}
    return fout;*/
    /*
    unsigned int fout = categories.deref(catite).first + categories.deref(catite).second;
	unsigned int to_move = asize & 0x7FFFFFF;
	uint32_t prec;

	//while


	uint32_t curcat = categories.getSize()-1;
	while(curcat != catite){
        prec = categories.deref(catite).first;
        memmove_routine(to_move, prec);
        to_move = prec;
        categories.deref(catite).first++;
        catite--;
	}
    return fout;*/

}
LFHTEMP D& CategoryHashmap<K,D,C,F,CF>::operator[](const K & key){
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF){
        // adding categoryless;
        if (((asize - 1)^(asize)) == ((asize << 1) -1)){
            if (asize & 0x80000000) asize &= 0x7FFFFFFf;
            else rehash( (asize == 0) ? 0 : hash_mag + 1);
        }
        if ((categories.getSize() != 0)&&(categories.deref(0).first == nb_categoryless)){
            uint32_t curc =0;
            while(curc+1 < categories.getSize()){
                if (categories.deref(curc).first + categories.deref(curc).second < categories.deref(curc+1).first) break;
                curc++;
            }
            while (curc > 0){
                memmove_routine(categories.deref(curc).first + categories.deref(curc).second, categories.deref(curc).first);
                categories.deref(curc).first++;
                curc--;
            }
            memmove_routine(categories.deref(curc).first + categories.deref(curc).second, categories.deref(curc).first);
            categories.deref(curc).first++;
        }

        j = nb_categoryless++;
        asize++;
        heap[j].first.k = key;
        unsigned int i = hashpos( F::makeSeed(key) );
        heap[j].second = hash[i];
        //printf("%i has %i as next!\n", j, heap[j].second);
        hash[i] = j;
        ExOp::toZero(heap[j].first.d);
        //printf("just inserted an category less %X\n", (int)this);
        //this->show();
        //fflush(stdout);
    }
    return heap[j].first.d;
}
LFHTEMP D CategoryHashmap<K,D,C,F,CF>::operator[](const K & key) const{
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF) {D fout; ExOp::toZero(fout); return fout;}
    return heap[j].first.d;
}
LFHTEMP void CategoryHashmap<K,D,C,F,CF>::getCategoryIndexes(const C& c, unsigned int& first, unsigned int& size)const{
    //ScopeVerbose sv("category index nv");
    unsigned int j = categories.find(c);
	if (j == 0xFFFFFFFF){
		size = 0;
	}else{
		first = categories.deref(j).first;
		size = categories.deref(j).second;
	}
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::insert(const K& key, const D& data, const C& category){
   // ScopeVerbose sv("insert category nv");
    if (((asize - 1)^(asize)) == ((asize << 1) -1)){
        if (asize & 0x80000000) asize &= 0x7FFFFFFf;
        else rehash( (asize == 0) ? 0 : hash_mag + 1);
    }
    unsigned int p = getNewSlotForCategory_routine(category);
    asize++;
    heap[p].first.k = key;
    unsigned int i = hashpos( F::makeSeed(key) );
    heap[p].second = hash[i];
    hash[i] = p;
    heap[p].first.d = data;

  //  printf("end insert: %X\n", (int) this);
  //  this->show();
}


LFHTEMP uint32_t CategoryHashmap<K,D,C,F,CF>::findCategory(uint32_t ite)const{
    if (ite < nb_categoryless) return 0xFFFFFFFF;
    int a = 0;
    int b = categories.getSize() - 1;
    int c;
    while(a < b){
        c = (a + b + 1) >> 1;
        if (ite >= categories.deref(c).first) a = c;
        else b = c - 1;
    }
    return a;
}
LFHTEMP C CategoryHashmap<K,D,C,F,CF>::getCategory(const K &key)const{
    uint32_t catite = this->findCategory(this->find(key));
    if (catite == 0xFFFFFFFF) {printf("HAS NO CATEGORY... illegal! \n"); LFH_exit(1);}
    return categories.deref_key(catite);
}


LFHTEMP C CategoryHashmap<K,D,C,F,CF>::deref_category(uint32_t ite)const{
    uint32_t catite = this->findCategory(ite);
    if (catite == 0xFFFFFFFF) {printf("HAS NO CATEGORY... illegal! \n"); LFH_exit(1);}
    return categories.deref_key(catite);
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::erase(const K& key){
    uint32_t ite = this->find(key);
    if (ite == 0xFFFFFFFF) {printf("does not exist in categoryhash!\n"); return;}
    else return this->erase_from_iterator(ite);
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::erase_from_iterator(uint32_t ite){
    remove_links_routine(ite, this->findCategory(ite));
    uint32_t rsize = (asize & 0x7FFFFFFF);
    if (((asize)^(asize-1))> rsize-1){
		if (rsize-1 == 0){
			asize =0;
			delete[](heap);
			return;
		}else{
			if ((asize & 0x80000000) == 0) asize |= 0x80000000;
			else{
				// needs down alloc, size is of from 2^i - 1
				{//:
					pair<KeyElem<K,D>, unsigned int >* swp = new pair<KeyElem<K,D>, unsigned int >[(rsize) << 1];
                    LFH_NICE_ALLOCERROR(swp,"Could not allocate %i bytes (critical)\n",(int)(sizeof(pair<KeyElem<K,D>, unsigned int >) * rsize*2))

					int i,j;


					for(i=0; i< nb_categoryless;i++) ExOp::toMemmove(swp[i].first, heap[i].first);
					uint32_t k=nb_categoryless<<1;
                    for(j=0;j<categories.getSize();j++){
                        for(i=0; i< categories.deref(j).second;i++){
                            ExOp::toMemmove(swp[k+i].first, heap[categories.deref(j).first + i].first);
                        }
                        categories.deref(j).first = k;
                        k+= (categories.deref(j).second << 1);
                    }
					delete[](heap);
					heap = swp;
					hash_mag--;
					delete[](hash);
                    hash = new unsigned int[1 << hash_mag];
                    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))

                    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);

                    for(i=0; i< nb_categoryless;i++){
                        unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[i].first) ) );
                        heap[i].second = hash[s];
                        hash[s] = i;
                    }

                    for(j=0;j<categories.getSize();j++){
                        for(i=0; i< categories.deref(j).second;i++){
                            unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[categories.deref(j).first + i].first)) );
                            heap[categories.deref(j).first + i].second = hash[s];
                            hash[s] = categories.deref(j).first + i;
                        }
                    }
				}//:
			}
		}
	}
    asize--;
    // alloc down check? bah...
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::setCategory(const K& key, const C& category){
    uint32_t ite = this->find(key);
    uint32_t catite = this->findCategory(ite);
    if (category == categories.deref_key(catite)) return; // no change!
    D swapbuffer;
    ExOp::toMemmove(swapbuffer, heap[ite].first.d);
    remove_links_routine(ite, catite);
    ite = getNewSlotForCategory_routine(category);
    heap[ite].first.k = key;
    unsigned int i = hashpos( F::makeSeed(key) );
    heap[ite].second = hash[i];
    hash[i] = ite;
    ExOp::toMemmove(heap[ite].first.d, swapbuffer);
}

LFHTEMP void CategoryHashmap<K,D,C,F,CF>::rehash(unsigned char _new_mag){
 //   ScopeVerbose sv("rehash nv");
    pair< KeyElem<K,D> , unsigned int >* old_heap = heap;
    heap = new pair< KeyElem<K,D> , unsigned int >[1 << _new_mag];
    LFH_NICE_ALLOCERROR(heap,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))

    unsigned int i,j,k;
    if (asize){ // was non-empty,moves entries
        for(i=0; i< nb_categoryless;i++) ExOp::toMemmove(heap[i], old_heap[i]);

        if (_new_mag == hash_mag+1){
            for(j=0;j<categories.getSize();j++){
                k = categories.deref(j).first << 1;
                for(i=0; i< categories.deref(j).second;i++){
                    ExOp::toMemmove(heap[k+i], old_heap[categories.deref(j).first + i]);
                }
                categories.deref(j).first = k;
            }
        }else{
            for(j=0;j<categories.getSize();j++){
                k = categories.deref(j).first + (1 << (_new_mag -1));
                for(i=0; i< categories.deref(j).second;i++){
                    ExOp::toMemmove(heap[k+i], old_heap[categories.deref(j).first + i]);
                }
                categories.deref(j).first = k;
            }
        }
        delete[](old_heap);
        delete[](hash);
    } else nb_categoryless =0;
    hash = new unsigned int[1 << _new_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << _new_mag))

    memset(hash, '\xFF',sizeof(unsigned int) << _new_mag);

    hash_mag = _new_mag;


    for(i=0; i< nb_categoryless;i++){
        unsigned int s = hashpos( F::makeSeed(heap[i].first.k) );
        heap[i].second = hash[s];
        hash[s] = i;
    }

    for(j=0;j<categories.getSize();j++){
        for(i=0; i< categories.deref(j).second;i++){
            unsigned int s = hashpos( F::makeSeed(heap[categories.deref(j).first + i].first.k) );
            heap[categories.deref(j).first + i].second = hash[s];
            hash[s] = categories.deref(j).first + i;
        }
    }
}

LFHTEMP CategoryHashmap<K,D,C,F,CF>& CategoryHashmap<K,D,C,F,CF>::toMemmove(CategoryHashmap<K,D,C,F,CF>& other){
    if (asize) delete[](hash);
    categories.toMemmove(other.categories);
    heap = other.heap;
    hash = other.hash;
    hash_mag = other.hash_mag;
    asize = other.asize;
    other.asize =0;
    return *this;
}
LFHTEMP CategoryHashmap<K,D,C,F,CF>& CategoryHashmap<K,D,C,F,CF>::operator=(const CategoryHashmap<K,D,C,F,CF>& other){
    if (asize) delete[](hash);
    categories = other.categories;
    hash_mag = other.hash_mag;
    hash = new unsigned int[1 << hash_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))
    memcpy(hash,other.hash,sizeof(unsigned int)<< hash_mag);
    return *this;
}
LFHTEMP void CategoryHashmap<K,D,C,F,CF>::show(FILE*f, int level)const{
    if (asize == 0) { fprintf(f, "Un-initialized CategoryHashMap\n");}
    else{
        fprintf(f, "CategoryHashmap, nbitems=%i and nbcategories=%i\n", asize & 0x7FFFFFFF, categories.getSize());
        unsigned int i;
        unsigned int ttt;
        for(i=0;(i>> hash_mag) == 0;i++){
               fprintf(f,"%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].second) {fprintf(f,"\t");ExOp::show(heap[j].first,f,2);}
               fprintf(f,"\n");
            }
        fprintf(f, "its categories are:\n");
        for(i=0;i<categories.heap.getSize();i++){
            ExOp::show(categories.heap[i].k.k,f,1);
            for(ttt= 0; ttt < categories.deref(i).second;ttt++){
                fprintf(f,"\t");ExOp::show(heap[ttt + categories.deref(i).first].first,f,2);
            }
            fprintf(f,"\n");
        }
    }
}

#undef LFHTEMP
#define LFHTEMP template<class K, class C, class F,class H>

LFHTEMP CategoryHashmap<K,void,C,F,H>::~CategoryHashmap(){
    if (asize) delete[](hash);
}

LFHTEMP unsigned int CategoryHashmap<K,void,C,F,H>::hashpos(unsigned int seed) const{return (seed & ((0xFFFFFFFF << hash_mag)^ 0xFFFFFFFF));}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::swap_and_pop(unsigned int ite){
   //TODO
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::move_up_routine(unsigned int from, unsigned int unto){
	ExOp::toMemmove(heap[unto], heap[from]);
	unsigned int tmp = hashpos( F::makeSeed(ExOp::getIndex(heap[unto].first)) );
	if (hash[tmp] == from) hash[tmp] = unto;
	else{
		tmp = hash[tmp];
		while(heap[tmp].second != from) tmp = heap[tmp].second;
		heap[tmp].second = unto;
	}
}


LFHTEMP unsigned int CategoryHashmap<K,void,C,F,H>::find(const typename ExCo<K>::INDEX_TYPE  &key) const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = hash[hashpos( F::makeSeed(key) )];
    for(;j != 0xFFFFFFFF;j = heap[j].second){
        if (key == ExOp::getIndex(heap[j].first)) break;
    }
    return(j);
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::remove_links_routine(uint32_t target, uint32_t catite){
    unsigned int *j = &(hash[hashpos( F::makeSeed(heap[target].first.getIndex()) )]);
    while((*j) != target) j = &heap[*j].second;
    (*j) = heap[target].second;
    if (catite == 0xFFFFFFFF) {nb_categoryless--; this->memmove_routine(target, nb_categoryless);}
    else{categories.deref(catite).second--; this->memmove_routine(target, categories.deref(catite).first + categories.deref(catite).second);}
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::memmove_routine(uint32_t target, uint32_t source){
    if (target == source) return;
    unsigned int *j = &(hash[hashpos( F::makeSeed(heap[source].first.getIndex()) )]);
    while((*j) != source) j = &heap[*j].second;
    *j = target;
    heap[target].second = heap[source].second;
    ExOp::toMemmove(heap[target].first, heap[source].first);
}

LFHTEMP unsigned int CategoryHashmap<K,void,C,F,H>::enlargeback_routine(unsigned int catite){
    uint32_t curc =catite;
    while(curc+1 < categories.getSize()){
        if (categories.deref(curc).first + categories.deref(curc).second < categories.deref(curc+1).first) break;
        curc++;
    }
    if (categories.deref(curc).first + categories.deref(curc).second == (1u << hash_mag)){
        for(curc = catite-1;curc>0; curc--){
            if (categories.deref(curc-1).first + categories.deref(curc-1).second < categories.deref(curc+1).first) break;
        }
        while(curc < catite){
            categories.deref(curc).first--;
            memmove_routine(categories.deref(curc).first, categories.deref(curc).first + categories.deref(curc).second);
            curc++;
        }
        categories.deref(curc).first--;
        return categories.deref(curc).first;
    }
    while (curc > catite){
        memmove_routine(categories.deref(curc).first + categories.deref(curc).second, categories.deref(curc).first);
        categories.deref(curc).first++;
        curc--;
    }
    return categories.deref(catite).first + categories.deref(catite).second;
}
LFHTEMP K& CategoryHashmap<K,void,C,F,H>::operator[](const typename ExCo<K>::INDEX_TYPE &key){
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF){
		printf("Could not find an entry in hashtable for "); ExOp::show(key,stderr);
        LFH_exit(1);
    }
    return heap[j].first;
}
LFHTEMP K CategoryHashmap<K,void,C,F,H>::operator[](const typename ExCo<K>::INDEX_TYPE &key) const{
    unsigned int j = find(key);
    if (j == 0xFFFFFFFF) {K fout; ExOp::toZero(fout); return fout;}
    return heap[j].first;
}
LFHTEMP void CategoryHashmap<K,void,C,F,H>::getCategoryIndexes(const C& c, unsigned int& first, unsigned int& size)const{
    //ScopeVerbose sv("category index void");
    unsigned int j = categories.find(c);
	if (j == 0xFFFFFFFF) {
		size =0;
	}else{
		first = categories.deref(j).first;
		size = categories.deref(j).second;
	}
}



LFHTEMP unsigned int CategoryHashmap<K,void,C,F,H>::insert(const K& key, const C& category){
    //ScopeVerbose sv("insert category void");
    if (((asize - 1)^(asize)) == ((asize << 1) -1)){
        if (asize & 0x80000000) asize &= 0x7FFFFFFf;
        else rehash( (asize == 0) ? 0 : hash_mag + 1);
    }
    unsigned int j = categories.find(category);
    unsigned int p,b;
    if (j == 0xFFFFFFFF){ // category does not exist...
        if (categories.getSize() == 0) p=nb_categoryless;
        else{
            b=categories.getSize()-1;
            p = categories.deref(b).first + categories.deref(b).second;
        }
        categories[category] = pair<unsigned int, unsigned int >(p,1);
    }else{
        b = (j == 0) ? nb_categoryless : categories.deref(j-1).first + categories.deref(j-1).second;
        if (b == categories.deref(j).first)  p = enlargeback_routine(j);
        else p = (--categories.deref(j).first);
        categories.deref(j).second++;
    }

    asize++;
    heap[p].first = key;
    unsigned int i = hashpos( F::makeSeed(ExOp::getIndex(key)) );
    heap[p].second = hash[i];
    hash[i] = p;
	return p;
}



LFHTEMP uint32_t CategoryHashmap<K,void,C,F,H>::findCategory(uint32_t ite)const{
    if (ite < nb_categoryless) return 0xFFFFFFFF;
    int a = 0;
    int b = categories.getSize() - 1;
    int c;
    while(a < b){
        c = (a + b + 1) >> 1;
        if (ite >= categories.deref(c).first) a = c;
        else b = c - 1;
    }
    return a;
}
LFHTEMP C CategoryHashmap<K,void,C,F,H>::getCategory(const typename ExCo<K>::INDEX_TYPE &key)const{
    uint32_t catite = this->findCategory(this->find(key));
    if (catite == 0xFFFFFFFF) {printf("HAS NO CATEGORY... illegal! \n"); LFH_exit(1);}
    return categories.deref_key(catite);
}
LFHTEMP C CategoryHashmap<K,void,C,F,H>::deref_category(uint32_t ite)const{
    uint32_t catite = this->findCategory(ite);
    if (catite == 0xFFFFFFFF) {printf("HAS NO CATEGORY... illegal! \n"); LFH_exit(1);}
    return categories.deref_key(catite);
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::erase(const typename ExCo<K>::INDEX_TYPE &key){
    uint32_t ite = this->find(key);
    if (ite == 0xFFFFFFFF) {printf("does not exist in categoryhash!\n"); return;}
    else return this->erase_from_iterator(ite);
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::erase_from_iterator(uint32_t ite){
    remove_links_routine(ite, this->findCategory(ite));
    uint32_t rsize = (asize & 0x7FFFFFFF);
    if (((asize)^(asize-1))> rsize-1){
		if (rsize-1 == 0){
			asize =0;
			delete[](heap);
			return;
		}else{
			if ((asize & 0x80000000) == 0) asize |= 0x80000000;
			else{
				// needs down alloc, size is of from 2^i - 1
				{//:
					pair<K, unsigned int >* swp = new pair< K , unsigned int >[(rsize) << 1];
					uint32_t i,j;
                    LFH_NICE_ALLOCERROR(swp,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) * rsize*2))


					for(i=0; i< nb_categoryless;i++) ExOp::toMemmove(swp[i].first, heap[i].first);
					uint32_t k=nb_categoryless<<1;
                    for(j=0;j<categories.getSize();j++){
                        for(i=0; i< categories.deref(j).second;i++){
                            ExOp::toMemmove(swp[k+i].first, heap[categories.deref(j).first + i].first);
                        }
                        categories.deref(j).first = k;
                        k+= (categories.deref(j).second << 1);
                    }
					delete[](heap);
					heap = swp;
					hash_mag--;
					delete[](hash);
                    hash = new unsigned int[1 << hash_mag];
                    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))

                    memset(hash, '\xFF',sizeof(unsigned int) << hash_mag);

                     for(i=0; i< nb_categoryless;i++){
                            unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[i].first) ) );
                            heap[i].second = hash[s];
                            hash[s] = i;
                        }

                    for(j=0;j<categories.getSize();j++){
                        for(i=0; i< categories.deref(j).second;i++){
                            unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[categories.deref(j).first + i].first)) );
                            heap[categories.deref(j).first + i].second = hash[s];
                            hash[s] = categories.deref(j).first + i;
                        }
                    }
				}//:
			}
		}
	}
    asize--;
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::rehash(unsigned char _new_mag){
    pair< K , unsigned int >* old_heap = heap;
    heap = new pair< K , unsigned int >[1 << _new_mag];
    LFH_NICE_ALLOCERROR(heap, "Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << _new_mag))

    unsigned int i,j,k;
    if (asize){ // was non-empty,moves entries

        for(i=0; i< nb_categoryless;i++) ExOp::toMemmove(heap[i], old_heap[i]);

        if (_new_mag == hash_mag+1){
            for(j=0;j<categories.getSize();j++){
                k = categories.deref(j).first << 1;
                for(i=0; i< categories.deref(j).second;i++){
                    ExOp::toMemmove(heap[k+i], old_heap[categories.deref(j).first + i]);
                }
                categories.deref(j).first = k;
            }
        }else{
            for(j=0;j<categories.getSize();j++){
                k = categories.deref(j).first + (1 << (_new_mag -1));
                for(i=0; i< categories.deref(j).second;i++){
                    ExOp::toMemmove(heap[k+i], old_heap[categories.deref(j).first + i]);
                }
                categories.deref(j).first = k;
            }
        }
        delete[](old_heap);
        delete[](hash);
    } else nb_categoryless =0;
    hash = new unsigned int[1 << _new_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << _new_mag))

    memset(hash, '\xFF',sizeof(unsigned int) << _new_mag);

    hash_mag = _new_mag;


    for(i=0; i< nb_categoryless;i++){
        unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[i].first) ) );
        heap[i].second = hash[s];
        hash[s] = i;
    }

    for(j=0;j<categories.getSize();j++){
        for(i=0; i< categories.deref(j).second;i++){
            unsigned int s = hashpos( F::makeSeed(ExOp::getIndex(heap[categories.deref(j).first + i].first)) );
            heap[categories.deref(j).first + i].second = hash[s];
            hash[s] = categories.deref(j).first + i;
        }
    }
}

LFHTEMP CategoryHashmap<K,void,C,F,H>& CategoryHashmap<K,void,C,F,H>::toMemmove(CategoryHashmap<K,void,C,F,H>& other){
    if (asize != 0) delete[](hash);
    categories.toMemmove(other.categories);
    heap = other.heap;
    hash = other.hash;
    hash_mag = other.hash_mag;
    asize = other.asize;
    other.asize =0;
    return *this;
}
LFHTEMP CategoryHashmap<K,void,C,F,H>& CategoryHashmap<K,void,C,F,H>::operator=(const CategoryHashmap<K,void,C,F,H>& other){
    if (asize != 0) delete[](hash);
    categories = other.categories;
    hash_mag = other.hash_mag;
    hash = new unsigned int[1 << hash_mag];
    LFH_NICE_ALLOCERROR(hash,"Could not allocate %i bytes (critical)\n",(int)(sizeof(unsigned int) << hash_mag))
    memcpy(hash,other.hash,sizeof(unsigned int)<< hash_mag);
    return *this;
}

LFHTEMP void CategoryHashmap<K,void,C,F,H>::show(FILE*f, int level)const{
    if (asize == 0) { fprintf(f, "Un-initialized CategoryHashMap\n");}
    else{
        fprintf(f, "CategoryHashmap, nbitems=%i and nbcategories=%i\n", asize & 0x7FFFFFFF, categories.getSize());
        unsigned int i;
        unsigned int ttt;
        for(i=0;(i>> hash_mag) == 0;i++){
               fprintf(f,"%i:",i);
               for(unsigned int j= hash[i]; j != 0xFFFFFFFF;j = heap[j].second) {fprintf(f,"\t");ExOp::show(heap[j].first,f,2);}
               fprintf(f,"\n");
            }
        fprintf(f, "its categories are:\n");
        for(i=0;i<categories.getSize();i++){
            ExOp::show(categories.deref_key(i),f,1);
            for(ttt= 0; ttt < categories.deref(i).second;ttt++){
                fprintf(f,"\t");ExOp::show(heap[ttt + categories.deref(i).first].first,f,2);
            }
            fprintf(f,"\n");
        }


    }
}

LFHTEMP bool CategoryHashmap<K,void,C,F,H>::QueryIterator::findFirst(const C &cat){
    unsigned int j = target.categories.find(cat);
	if (j == 0xFFFFFFFF) return false;
    ite = target.categories.deref(j).first;
	nb = target.categories.deref(j).second;
	return (nb > 0);
}

LFHTEMP bool CategoryHashmap<K,void,C,F,H>::QueryIterator::findNext(){ite++;return (--nb > 0);}
LFHTEMP const K* CategoryHashmap<K,void,C,F,H>::QueryIterator::operator->()const{return &(target.heap[ite].first);}


#undef LFHTEMP
#define LFHTEMP template<class IC, class D>

LFHTEMP void MaskMap<IC,D>::insert(const Tuple<IC>& key, const D& value){
    hashmap[key] = value;
}
LFHTEMP void MaskMap<IC,D>::getIntersection(Vector<D> &res, const Tuple<IC>& value, const Tuple<IC>& mask) const{
    uint32_t i,j;
    for(i=0;i<hashmap.getSize();i++){
        for(j=0;j<value.getSize();j++) if ((hashmap.deref_key(i)[j] ^ value[j]) & mask) break;
        if (j == value.getSize()) res.push_back(hashmap.deref(i));
    }
}

#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int A>

LFHTEMP void HeapTree<C,A>::insert(const C &what){
    asyncbuffer.insert(what);
}
LFHTEMP bool HeapTree<C,A>::pop(C& fout){
    typename AsyncBuffer<C,A>::ReadAccess racc(asyncbuffer);
    C target;
    while(racc.needsUpdate()) heaptree.insert(racc.update());
    if (heaptree.isEmpty()) return false;
    fout = heaptree.pop();
return true;}
LFHTEMP bool HeapTree<C,A>::top(C& fout){
    typename AsyncBuffer<C,A>::ReadAccess racc(asyncbuffer);
    C target;
    while(racc.needsUpdate()) heaptree.insert(racc.update());
    if (heaptree.isEmpty()) return false;
    fout = heaptree.top();
return true;}

LFHTEMP void HeapTree<C,A>::show(FILE*f, int level)const{
	fprintf(f, (hasunsorted) ? "Async HeapTree with 1 unsorted" : "Async HeapTree");
	asyncbuffer.show(f);
	if (data.getSize() == 0) fprintf(f, "(empty)\n");
	else{
		for( int i = (hasunsorted) ?  0 : 1; i< data.getSize();i++){
			fprintf(f,"%i:", i); ExOp::show(data[i],1);
			fprintf(f,"\n");
		}
	}
}

#undef LFHTEMP
#define LFHTEMP template <class C>
LFHTEMP void HeapTree<C,0>::insert(const C &what){
	if (hasunsorted){

	int cur = data.getSize();
	data.push_back();
	//push_back(KeyElem<Key, Node>());

	while(true){
	if ((cur == 1)||( data[cur >>1] < what )) break;
	ExOp::toMemmove(data[cur], data[cur  >>1]);
	cur  >>= 1;
	}
	data[cur] = what;
	}else{
		hasunsorted = true;
		data[0] = what; // quick insert!
	}
}
LFHTEMP HeapTree<C,0>& HeapTree<C>::operator<<=(C &what){ // memory insert!
    if (hasunsorted){
        int cur = data.getSize();
        data.push_back();
        //push_back(KeyElem<Key, Node>());
        while(true){
            if ((cur == 1)||( data[cur >>1] < what )) break;
            ExOp::toMemmove(data[cur], data[cur  >>1]);
            cur  >>= 1;
        }
        ExOp::toMemmove(data[cur], what);
	}else{
		hasunsorted = true;
		ExOp::toMemmove(data[0],what); // quick insert!
	}
    return *this;
}
LFHTEMP C HeapTree<C,0>::pop(){
	C swap;
	unsigned int s = data.getSize();
	if (!hasunsorted){
        if (s == 1) {fprintf(stderr,"poping form an empty heap! (required to check first)\n");LFH_exit(1);}
		ExOp::toMemmove(data[0],data[1]);
		ExOp::toMemmove(swap,data.last());
		data.pop_back();s--;
		if (s == 1) return(data[0]);
	}else{
		hasunsorted = false;
		if ((s ==1)||(data[0] < data[1])) {
			return(data[0]);
		}else{
			ExOp::toMemmove(swap,data[0]);
			ExOp::toMemmove(data[0], data[1]);
		}

	}

	unsigned int cur = 2;
	while(true){
		if (cur+1 < s) {
			if (data[cur+1] < data[cur]) cur++;
		}else{
			if (cur+1 !=  s) break;
		}
		if (data[cur] > swap) break;
		ExOp::toMemmove(data[cur >> 1], data[cur]);
		cur <<=1;
	}
	ExOp::toMemmove(data[cur >> 1],swap);

	return(data[0]);
	}
LFHTEMP void HeapTree<C,0>::pop_and_insert(const C &what){
	C swap;
	unsigned int s = data.getSize();
	if (!hasunsorted){
        if (s == 1) {fprintf(stderr,"poping form an empty heap! (required to check first)\n");LFH_exit(1);}
		ExOp::toMemmove(data[0],data[1]);
		ExOp::toMemmove(swap,data.last());
		data.pop_back();s--;
		data[0] = what;
		hasunsorted = true;
		if (s == 1) return;;
	}else{
		if ((s ==1)||(data[0] < data[1])) {
			data[0] = what;
			return;
		}else{
			ExOp::toMemmove(swap,data[0]);
			ExOp::toMemmove(data[0], data[1]);
		}
		data[0] = what;
	}

	unsigned int cur = 2;
	while(true){
		if (cur+1 < s) {
			if (data[cur+1] < data[cur]) cur++;
		}else{
			if (cur+1 !=  s) break;
		}
		if (data[cur] > swap) break;
		ExOp::toMemmove(data[cur >> 1], data[cur]);
		cur <<=1;
	}
	ExOp::toMemmove(data[cur >> 1],swap);
}

LFHTEMP void HeapTree<C,0>::update(unsigned int cur){// Elem was changed! update its position!
	if (cur == 0) return; // unsorted position, ignore!
	C swap;
	int s = this->size();
	if ((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *(*this)[cur])){
		// got smaller
		swap = (*this)[cur];
		do{
			(*this)[cur] = (*this)[cur >> 1];
			cur = cur >> 1;
		}while((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *swap));
		(*this)[cur] = swap;
		return;
	}
	if ((cur << 1) < s){
		if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
			if (*(*this)[cur << 1] < *(*this)[cur]){
				swap = (*this)[cur];
				(*this)[cur] = (*this)[cur << 1];
				cur = cur << 1;
			}else return;
		}else{
			if (*(*this)[1+(cur << 1)] < *(*this)[cur]){
				swap = (*this)[cur];
				(*this)[cur] = (*this)[1 + (cur << 1)];
				cur = 1 | (cur << 1);
			}else return;
		}
	} else return;
	while((cur << 1) < s){
		if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
			if (*(*this)[cur << 1] < *swap){
				(*this)[cur] = (*this)[cur << 1];
				cur = cur << 1;
			}else break;
		}else{
			if (*(*this)[1+(cur << 1)] < *swap){
				(*this)[cur] = (*this)[1 + (cur << 1)];
				cur = 1 | (cur << 1);
			}else break;
		}
	}
	(*this)[cur] = swap;
}

LFHTEMP void HeapTree<C,0>::show(FILE*f, int level)const{
	fprintf(f, (hasunsorted) ? "HeapTree with 1 unsorted" : "HeapTree");
	if (data.getSize() == 0) fprintf(f, "(empty)\n");
	else{
		for( int i = (hasunsorted) ?  0 : 1; i< data.getSize();i++){
			fprintf(f,"%i:", i); ExOp::show(data[i],f,1);
			fprintf(f,"\n");
		}
	}
}
	// pointer specialization!!!
#undef LFHTEMP
#define LFHTEMP template <class C>
	// pointer specialization!!!

	LFHTEMP HeapTree<C*,0>::HeapTree() : hasunsorted(false){
		this->setSize(1);
	}

	LFHTEMP bool HeapTree<C*,0>::isEmpty(){
		return((!hasunsorted) && (this->size() == 1));
	}

	//template <class Key, class Node>
	//Node& HeapTree<Key,Node>::operator()(const Key & where){ // yet a method to insert

	//	}
	LFHTEMP void HeapTree<C*,0>::insert(C* what){
		unsigned int cur;
		if (hasunsorted){

			cur = this->size();
			this->push_back();
			//push_back(KeyElem<Key, Node>());

			while(true){
				if ((cur == 1)||( *(*this)[cur >>1] < *what )) break;
				(*this)[cur] = (*this)[cur  >>1];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur  >>= 1;
			}
		}else{
			hasunsorted = true;
			cur = 0;// quick insert!
		}
		(*this)[cur] = what;

		if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(*what) = cur;

	}
	LFHTEMP C& HeapTree<C*,0>::pop(){
		C* swap;

		if (!hasunsorted){
			(*this)[0] = (*this)[1];
			swap = (this->last());
			this->pop_back();
		}else{
			hasunsorted = false;
			if (*(*this)[0] < *(*this)[1]) {
				return(*(*this)[0]);
			}else{
				swap = (*this)[0];
				(*this)[0] = (*this)[1];
			}

		}

		int s = this->size();
		int cur = 1;
		while((cur << 1) < s){
			if (((cur << 1) == s-1)||(*(*this)[cur<<1] < *(*this)[1 | (cur<< 1)])){
				if (*swap < *(*this)[cur<<1]) break;
				(*this)[cur] = (*this)[cur <<1];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur = cur << 1;
			}else{
				if (*swap < *(*this)[1 | (cur<<1)]) break;
				(*this)[cur] = (*this)[1 | (cur <<1)];
				if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
				cur = 1 | (cur << 1);
			}
		}

		(*this)[cur] = swap;
		if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;

		return(*((*this)[0]));
	}
	LFHTEMP void HeapTree<C*,0>::update(unsigned int cur){// Elem was changed! update its position!
		if (cur == 0) return; // unsorted position, ignore!
		C* swap;
		int s = this->size();
		if ((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *(*this)[cur])){
			// got smaller
			swap = (*this)[cur];
			do{
			(*this)[cur] = (*this)[cur >> 1];
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
			cur = cur >> 1;
			}while((cur & 0xFFFFFFFE)&&(*(*this)[cur >> 1] > *swap));
			(*this)[cur] = swap;
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
			return;
		}
			if ((cur << 1) < s){
				if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
					if (*(*this)[cur << 1] < *(*this)[cur]){
						swap = (*this)[cur];
						(*this)[cur] = (*this)[cur << 1];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = cur << 1;
					}else return;
				}else{
					if (*(*this)[1+(cur << 1)] < *(*this)[cur]){
						swap = (*this)[cur];
						(*this)[cur] = (*this)[1 + (cur << 1)];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = 1 | (cur << 1);
					}else return;
				}
			} else return;
			while((cur << 1) < s){
				if (((cur << 1) == s-1)||(*(*this)[cur << 1] < *(*this)[1 | (cur << 1)] )){
					if (*(*this)[cur << 1] < *swap){
						(*this)[cur] = (*this)[cur << 1];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = cur << 1;
					}else break;
				}else{
					if (*(*this)[1+(cur << 1)] < *swap){
						(*this)[cur] = (*this)[1 + (cur << 1)];
						if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;
						cur = 1 | (cur << 1);
					}else break;
				}
			}
			(*this)[cur] = swap;
			if (HeapTreeBackPointerOffset<C>::ans) HeapTreeBackPointerOffset<C>::getBackOffset(* (*this)[cur] ) = cur;

	}


#undef LFHTEMP
#define LFHTEMP template <class K, class D, class HF>

LFHTEMP void AsyncHashmap<K,D,HF>::rehash(uint32_t _hsize){
    if (asize != 0) delete[](dahash);
    hsize = _hsize;
    if (_hsize == 0) return;
    dahash = new uint32_t[hsize << 1];
    memset(dahash + hsize , '\xFF', hsize << 2);
    uint32_t i;
    for(i=0;i< asize;i++){
        uint32_t j = (HF::makeSeed(darray[i].k) & (hsize-1)) | hsize;
        dahash[i] = dahash[j];
        dahash[j] = i;
    }
    /*for(i=0;i< (hsize << 1);i++){
        printf("H[%i]: %i\n",i, dahash[i]);

    }*/
}
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::find_routine(const K &key) const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
    for(;j != 0xFFFFFFFF;j = dahash[j]){
        if (key == darray[j].k) break;
    }
return(j);}
LFHTEMP void AsyncHashmap<K,D,HF>::delayed_remove_routine(const K &key, std::unique_lock<std::mutex> &lck){
    // just kidding, not deleting ^^'
    //printf("delayed wait!\n"); fflush(stdout);
    uint32_t ind;
    foreach.accessFinalCondvar().wait(lck,[this, key, &ind]{
            if (asize == 0) {ind = 0xFFFFFFFF; return true;}
            if ((chunkIO & (3 << (((hsize >> 5) == 0) ?  ((asize-1) << 1) : (((asize-1) / (hsize >> 5)) & 0x1E)))) != 0) return false;
            uint32_t i = (HF::makeSeed(key) & (hsize-1)) | hsize;
            for(ind = dahash[i];;ind= dahash[ind]){
                if (ind == 0xFFFFFFFF) return true;
                if (key == darray[ind].k) break;
            }
            return ((chunkIO & (3 << (((hsize >> 5) == 0) ?  (ind << 1) : ((ind / (hsize >> 5)) & 0x1E)))) == 0);
        });
    if (ind == 0xFFFFFFFF) return;
    if (asize == 1){ // deleting the last one!
            //printf("delay del last!\n"); fflush(stdout);
        delete[](dahash);
        delete[](darray);
        asize =0;
    }else{
        //printf("delay del!"); fflush(stdout);
        if (key != darray[ind].k){
            printf("crazy, did not match!\n");
        }
        uint32_t i = (HF::makeSeed(key) & (hsize-1)) | hsize;
        uint32_t j = dahash[i];
        if (j == ind) dahash[i] = dahash[j];
        else{
            while (dahash[j] != ind) j = dahash[j];
            dahash[j] = dahash[ind];
        }
        moveEntry_routine(asize-1, ind);
        asize--;

        if (asize <= (hsize >> 2)){
            //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO==0)||(asize > (hsize >> 2));});
            //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)chunkIO);
            if (asize <= (hsize >> 2)){
                //printf("delay dadel downalloc!\n"); fflush(stdout);
                KeyElem<K, D>* swp = darray;
                darray = new KeyElem<K, D>[hsize >> 1];
                for(i=0;i< asize; i++) ExOp::toMemmove(darray[i], swp[i]);
                rehash(hsize >> 1);
                delete[](swp);
            }
        }
    }
return;}
LFHTEMP void AsyncHashmap<K,D,HF>::moveEntry_routine(uint32_t from, uint32_t to){
    if (from == to) return;
    uint32_t i;
    uint32_t j = dahash[i = ((HF::makeSeed(darray[from].k) & (hsize-1)) | hsize)];
    if (j == from) dahash[i] = to;
    else{
        while (dahash[j] != from) j = dahash[j];
        dahash[j] = to;
    }
    dahash[to] = dahash[from];
    darray[to] = std::move(darray[from]);
}
LFHTEMP AsyncHashmap<K,D,HF>::WriteAccess::~WriteAccess(){
	uint32_t tmp;
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); tmp = target.chunkIO.fetch_sub(exitcode); if ((tmp & exitcode) != exitcode) {printf("freed %X but %X is atomic previous value...\n", exitcode, tmp); exit(1);} foreach.accessFinalCondvar().notify_all();}
}
LFHTEMP void AsyncHashmap<K,D,HF>::WriteAccess::unlock(){
	uint32_t tmp;
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); tmp = target.chunkIO.fetch_sub(exitcode); if ((tmp & exitcode) != exitcode) {printf("freed %X but %X is atomic previous value...\n", exitcode, tmp); exit(1);} foreach.accessFinalCondvar().notify_all(); exitcode =0;}
}
LFHTEMP AsyncHashmap<K,D,HF>::ReadAccess::~ReadAccess(){
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); target.chunkIO.fetch_sub(exitcode); foreach.accessFinalCondvar().notify_all();}
}
LFHTEMP void AsyncHashmap<K,D,HF>::ReadAccess::unlock(){
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); target.chunkIO.fetch_sub(exitcode); foreach.accessFinalCondvar().notify_all(); exitcode =0;}
}
/*LFHTEMP bool AsyncHashmap<K,D,HF>::ReadAccess::find_only_routine(const K &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0; return false;}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == target.darray[j].k) break;
    }
    exitcode =0;
return(j != 0xFFFFFFFF);}*/
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::ReadAccess::find_read_syncroutine(const K &key){
    if (target.asize == 0) {exitcode = 0; return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == target.darray[j].k) break;
    }

    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        uint32_t tmp = (exitcode | (exitcode << 1));
        if ((target.chunkIO & tmp) == tmp) exitcode = 0;
        else {
            target.chunkIO.fetch_add(exitcode);
            elem = target.darray + j;
        }
    }
return(j);}
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::WriteAccess::try_write_routine(const K &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == target.darray[j].k) break;
    }
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 3 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        if ((target.chunkIO & exitcode) == 0) target.chunkIO.fetch_or(exitcode);
        else exitcode = 0;
    }
return(j);}
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::ReadAccess::block_routine(const K &key){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == target.darray[j].k) break;
    }
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        uint32_t tmp = (exitcode | (exitcode << 1));

        //if (((target.chunkIO | target.angryIO) & tmp) == tmp) {
            //foreach.fprintf_l(stdout,"Deadlock prone want %X got chunkIO %X\n", tmp, (uint32_t)target.chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this,tmp]{return ((target.chunkIO | target.angryIO) & tmp) != tmp;});
            //foreach.fprintf_l(stdout,"Woke up! got everything chunkIO %X\n", (uint32_t)target.chunkIO);
        //}

        target.chunkIO.fetch_add(exitcode);
    }
return(j);}
/* Delete entry from hashmap, needs to make sure no other reader exist
 */
LFHTEMP void AsyncHashmap<K,D,HF>::ReadAccess::remove(){
     if (exitcode != 0){
        // unlink from hash
        {std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
            uint32_t i;
            if (target.chunkIO & (exitcode << 1)){
                // cant write...
                target.chunkIO.fetch_sub(exitcode); exitcode =0;
                foreach.accessFinalCondvar().notify_all();
                target.delayed_remove_routine(elem->d, lck);
                elem = NULL;
                return;
            }

            // first try to swap/delete entry

            if (target.asize == 1){
                delete[](target.dahash);
                delete[](target.darray);
                target.asize =0;
                //printf("dadel del last! %i %p %p L=%c\n", target.asize, target.dahash, target.darray, lck ? 'L' : 'U'); fflush(stdout);
            }else{
                // needs to lock entry and end of hash!
                uint32_t endcode = 1 << (((target.hsize >> 5) == 0) ?  ((target.asize-1) << 1) : (((target.asize-1) / (target.hsize >> 5)) & 0x1E));
                if ((endcode == exitcode)||((target.chunkIO & (endcode | (endcode << 1)) ) == 0)) {
                    uint32_t ind = (uint32_t) (elem - target.darray);
                    //printf("dadel read! %i\n", ind); fflush(stdout);
                    if (ind != target.find_routine(elem->k)) printf("Disagree for position!\n", ind, target.find_routine(elem->k));
                    uint32_t j = target.dahash[i = ((HF::makeSeed(elem->k) & (target.hsize-1)) | target.hsize)];
                    if (j == ind) target.dahash[i] = target.dahash[j];
                    else{
                        while (target.dahash[j] != ind) j = target.dahash[j];
                        target.dahash[j] = target.dahash[ind];
                    }
                    target.moveEntry_routine(target.asize-1, ind);
                    target.asize--;
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                }else{ // cannot access end... need free it and to wait for it  and end of vector, implying remove is postponed!
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                    target.delayed_remove_routine(elem->k, lck);
                    elem = NULL;
                    return;
                }
                // cant exit if downalloc is needed
                if (target.asize <= (target.hsize >> 2)){
                    //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) target.chunkIO);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO==0)||(target.asize > (target.hsize >> 2));});
                    //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
                    if (target.asize <= (target.hsize >> 2)){
                        //printf("dadel downalloc!\n"); fflush(stdout);
                        KeyElem<K, D>* swp = target.darray;
                        target.darray = new KeyElem<K, D>[target.hsize >> 1];
                        for(i=0;i< target.asize; i++) target.darray[i] = std::move(swp[i]);
                        target.rehash(target.hsize >> 1);
                        delete[](swp);
                    }
                }
            }
            elem = NULL;
        }
    }
}
/* Delete entry from hashmap
 */
LFHTEMP void AsyncHashmap<K,D,HF>::WriteAccess::remove(){
    if (exitcode != 0){
        // unlink from hash
        {std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
            uint32_t i;

            // first try to swap/delete entry

            if (target.asize == 1){
                delete[](target.dahash);
                delete[](target.darray);
                target.asize =0;
            }else{
                // needs to lock entry and end of hash!
                uint32_t endcode = 3 << (((target.hsize >> 5) == 0) ?  ((target.asize-1) << 1) : (((target.asize-1) / (target.hsize >> 5)) & 0x1E));
                if ((endcode == exitcode)||((target.chunkIO & endcode) == 0)) {
                    uint32_t ind = (uint32_t) (elem - target.darray);
                    //printf("dadel write! %i\n", ind);
                    if (ind != target.find_routine(elem->k)) printf("Disagree for position!\n", ind, target.find_routine(elem->k));
                    uint32_t j = target.dahash[i = ((HF::makeSeed(elem->k) & (target.hsize-1)) | target.hsize)];
                    if (j == ind) target.dahash[i] = target.dahash[j];
                    else{
                        while (target.dahash[j] != ind) j = target.dahash[j];
                        target.dahash[j] = target.dahash[ind];
                    }
                    target.moveEntry_routine(target.asize-1, ind);
                    target.asize--;
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                }else{ // cannot access end... need free it and to wait for it  and end of vector, implying remove is postponed!
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                    target.delayed_remove_routine(elem->k, lck);
                    elem = NULL;
                    return;
                }
                // cant exit if downalloc is needed
                if (target.asize <= (target.hsize >> 2)){
                    //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) target.chunkIO);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO==0)||(target.asize > (target.hsize >> 2));});
                    //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
                    if (target.asize <= (target.hsize >> 2)){
                        KeyElem<K, D>* swp = target.darray;
                        target.darray = new KeyElem<K, D>[target.hsize >> 1];
                        for(i=0;i< target.asize; i++) ExOp::toMemmove(target.darray[i], swp[i]);
                        target.rehash(target.hsize >> 1);
                        delete[](swp);
                    }
                }
            }
            elem = NULL;
        }
    }
}
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::WriteAccess::block_routine(const K &key){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]) if (key == target.darray[j].k) break;
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        exitcode |= exitcode << 1;

        if ((target.chunkIO & exitcode) != 0) {
            // trying to stay calm ^_^
            //if ((target.angryIO & exitcode) != 0) foreach.accessFinalCondvar().wait(lck,[this]{return (target.angryIO & exitcode) == 0;});
            //target.angryIO |=  exitcode; // got angry!
            //if ((target.chunkIO & exitcode) != 0) {
               // foreach.fprintf_l(stdout,"writing block %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
                foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO & exitcode) == 0;});
                //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t) target.chunkIO);
            //}
            //target.angryIO ^=  exitcode; // is happy now!
        }
        target.chunkIO.fetch_or(exitcode);

    }
return(j);}

LFHTEMP AsyncHashmap<K,D,HF>::ConstIterator::~ConstIterator(){
    if (exitcode != 0){
        target.chunkIO.fetch_sub(exitcode);
        foreach.accessFinalCondvar().notify_all();
    }
}
LFHTEMP AsyncHashmap<K,D,HF>::ConstIterator::operator  bool (){
    if (target.asize == 0) return false;
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    // locks first and last (no delete and insert allowed, read/write can occur)
    exitcode = 1 | (0x55555555 << ((((target.hsize >> 5) == 0) ? target.asize - 1 : target.asize / (target.hsize >> 5)) & 0x1E));
   // if (exitcode & (target.angryIO |(target.chunkIO & (target.chunkIO >>1))) != 0){
        //foreach.fprintf_l(stdout,"iterator start %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
        foreach.accessFinalCondvar().wait(lck,[this]{return (exitcode & (target.angryIO | (target.chunkIO & (target.chunkIO >>1))) == 0);});
        //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
   //     target.chunkIO.fetch_add(exitcode);
   // }else
    target.chunkIO.fetch_add(exitcode);
    i = 0;
return true;}
LFHTEMP bool AsyncHashmap<K,D,HF>::ConstIterator::operator++(int){
    if (++i == target.asize){
        target.chunkIO.fetch_sub(exitcode); exitcode = 0;
        foreach.accessFinalCondvar().notify_all();
        return false;
    }
    uint32_t tmp = 1 << (((target.hsize >> 5) == 0) ?  (i << 1) : ((i / (target.hsize >> 5)) & 0x1E));
    if ((tmp & exitcode) == 0){
        target.chunkIO.fetch_sub(tmp >> 2);
        uint32_t  tmp2 = tmp | (tmp << 1);
        std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());

        if (((target.angryIO | target.chunkIO) & tmp2) == tmp2){
            foreach.accessFinalCondvar().notify_all();
            //foreach.fprintf_l(stdout,"Iterator move wants %X got chunkIO %X\n", tmp2, (uint32_t)target.chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this,tmp2]{return (((target.angryIO | target.chunkIO) & tmp2) != tmp2);});
            //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
            target.chunkIO.fetch_add(tmp);
        }else{
            target.chunkIO.fetch_add(tmp);
            foreach.accessFinalCondvar().notify_all();
        }
        exitcode ^= tmp | (tmp >> 2);
    }
return true;}
LFHTEMP void AsyncHashmap<K,D,HF>::remove(const K &key){
    // needs to lock entry and end of array!
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
return delayed_remove_routine(key, lck);}

LFHTEMP D& AsyncHashmap<K,D,HF>::createEntry_sync(const K key){
    uint32_t j;
	if (asize == 0) {
		hsize = 1;
		dahash = new uint32_t[2];
		memset(dahash + hsize , '\xFF', 4);
		this->rehash(1);
		darray = new KeyElem<K, D>[1];
		chunkIO = 0;
	}else if (asize == hsize){ // need to upalloc!!!
		this->rehash(hsize << 1);
		KeyElem<K, D>* swp = darray;
		darray = new KeyElem<K, D>[hsize];
		for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
		delete[](swp);
	}
    j = asize++;
    darray[j].k = key;
    uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[j] = dahash[tj];
    dahash[tj] = j;
return darray[j].d;}

LFHTEMP uint32_t AsyncHashmap<K,D,HF>::find_sync(const K key)const{
	if (asize == 0) return(0xFFFFFFFF);
    uint32_t j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
    for(;j != 0xFFFFFFFF;j = dahash[j]){
        if (key == ExOp::getIndex(darray[j])) break;
    }
return j;}


LFHTEMP uint32_t AsyncHashmap<K,D,HF>::try_create_routine(const K &key, uint32_t &exitcode){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    uint32_t j;
    if (asize == 0) {
        hsize = 1;
        darray = new KeyElem<K, D>[1];
        dahash = new uint32_t[2];
        chunkIO = 3;
        exitcode = 3;
        dahash[0] = 0xFFFFFFFF;
        dahash[1] = 0;
        return asize++;
    }else if (asize == hsize){
        // need to upalloc!!!
        if (chunkIO != 0) {exitcode = 0; return 0xFFFFFFFF;} // no can do...
        this->rehash(hsize << 1);
        KeyElem<K, D>* swp = darray;
        darray = new KeyElem<K, D>[hsize];
        for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
        delete[](swp);
    }else j = asize;
    exitcode = 3 << (((hsize >> 5) == 0) ?  (j << 1) : ((j / (hsize >> 5)) & 0x1E));
    if ((chunkIO & exitcode) != 0){exitcode = 0; return j;}
    chunkIO |= exitcode;
    darray[j].k = key;
    j = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[asize] = dahash[j];
    dahash[j] = asize;
    asize++;
return asize;}
LFHTEMP uint32_t AsyncHashmap<K,D,HF>::block_create_routine(const K &key, uint32_t &exitcode){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    uint32_t j;
    do{
        if (asize == 0) {
            //printf("created 1! %c\n", lck ? 'L' : 'U'); fflush(stdout);
            //fflush(stdout);
            this->rehash(1);
            darray = new KeyElem<K, D>[1];
            chunkIO = 0;
            exitcode = 3;
            break;
        }else if (asize == hsize){
            // need to upalloc!!!
            //printf("upalloc time! %c\n", lck ? 'L' : 'U'); fflush(stdout);
            do{
                if (chunkIO != 0) { // todo, make angry
                    //foreach.fprintf_l(stdout,"Deadlockprone want everything got chunkIO %X\n", (uint32_t)chunkIO);
                    //printf("sleep time! %c\n", lck ? 'L' : 'U'); fflush(stdout);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO==0)||(asize !=  hsize);});
                    //foreach.fprintf_l(stdout,"Woke up! got everything chunkIO %X\n", (uint32_t)chunkIO);
                    if (asize < hsize) {
                        //printf("something happened! %i\n", asize);
                        if (asize == 0){
                            this->rehash(1);
                            darray = new KeyElem<K, D>[1];
                            chunkIO = 0;
                        }
                        break;
                    } // some delete or up alloc happened
                }
                //printf("all good! %i %i %c\n", asize, hsize, lck ? 'L' : 'U'); fflush(stdout);
                this->rehash(hsize << 1);
                KeyElem<K, D>* swp = darray;
                darray = new KeyElem<K, D>[hsize];
                for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
                delete[](swp);
                j = asize++;
                exitcode = 3 << (((hsize >> 5) == 0) ?  (j << 1) : ((j / (hsize >> 5)) & 0x1E));
                chunkIO |= exitcode;

                darray[j].k = key;
                uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
                dahash[j] = dahash[tj];
                dahash[tj] = j;
                return j;
            }while(false);

        }
        exitcode = 3 << (((hsize >> 5) == 0) ?  (asize << 1) : ((asize / (hsize >> 5)) & 0x1E));
        if ((chunkIO & exitcode) == 0) break;
                    // trying to stay calm ^_^
            //if ((angryIO & exitcode) != 0) foreach.accessFinalCondvar().wait(lck,[this,exitcode]{return (angryIO & exitcode) == 0;});
            //angryIO |=  exitcode; // got angry!
        //    if ((chunkIO & exitcode) != 0) {
               // foreach.fprintf_l(stdout,"writing block %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
                foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO & (3 << (((hsize >> 5) == 0) ?  (asize << 1) : ((asize / (hsize >> 5)) & 0x1E)))) == 0;});

                //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t) target.chunkIO);
           //}
            //angryIO ^=  exitcode; // is happy now!
        //}
    }while(true);
    j = asize++;
    chunkIO |= exitcode;
    darray[j].k = key;
    uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[j] = dahash[tj];
    dahash[tj] = j;
    /*printf("created exity! %i %i %i %i %c\n", asize, hsize, tj, j, lck ? 'L' : 'U'); fflush(stdout);
    for(int i=0;i< (hsize << 1);i++){
        printf("H[%i]: %i\n",i, dahash[i]);
    }*/
    fflush(stdout);
return j;}
/*LFHTEMP uint32_t AsyncHashmap<K,D,HF>::find_create_routine(const K &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    unsigned int j;
    if (asize == 0) {
        j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
        for(;j != 0xFFFFFFFF;j = dahash[j]){
            if (key == ExOp::getIndex(darray[j]) ) return(j);
        }
    }

    uint32_t rsize = (asize & 0x7FFFFFFF);
    j = rsize;
    if ((((asize + 1u)^(asize))> rsize)&&((asize & 0x80000000) == 0)) { // going to up_alloc
        // TODO
       // heap.push_back();
       // heap[j].k.k = key;
       // rehash(hash_mag + 1);
    }else{
        uint32_t j = (HF::makeSeed(key) & (hsize-1)) | hsize;
        dahash[asize] = dahash[j];
        dahash[j] = asize;
    }
return asize++;}*/

LFHTEMP typename AsyncHashmap<K,D,HF>::WriteAccess AsyncHashmap<K,D,HF>::operator[](const K &key){ typename AsyncHashmap<K,D,HF>::WriteAccess fout(*this);

    // wait... that hash can be thrashed *but* the offset should be allowed to change... need to lock that first
    uint32_t ind = fout.block_routine(key);
    if (ind == 0xFFFFFFFF) ind = this->block_create_routine(key, fout.exitcode);
    //printf("daind%i%i%i%idniad\n", ind, ind, ind, ind);
    fout.elem = darray + ind;
return fout;}
LFHTEMP typename AsyncHashmap<K,D,HF>::ReadAccess AsyncHashmap<K,D,HF>::read(const K &key)const{ typename AsyncHashmap<K,D,HF>::ReadAccess fout(*this);
    uint32_t ind = fout.block_routine(key);
    if (ind != 0xFFFFFFFF) fout.elem = darray + ind;
return fout;}

LFHTEMP typename AsyncHashmap<K,D,HF>::ReadAccess AsyncHashmap<K,D,HF>::await(const K &key, atomic<int32_t>* can_block_flag) const{typename AsyncHashmap<K,D,HF>::ReadAccess fout(*this);
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    fout.find_read_syncroutine(key);
    if (fout.exitcode != 0) return fout;
    // does not exist, wait for it!
    if (0 == (int32_t) *can_block_flag) foreach.accessFinalCondvar().wait(lck,[&fout, &key, can_block_flag]{fout.find_read_syncroutine(key);return (fout.exitcode != 0)||(0 != (int32_t) *can_block_flag);});
return fout;}

LFHTEMP void AsyncHashmap<K,D,HF>::insert(const K key, D && newelem){ typename AsyncHashmap<K,D,HF>::WriteAccess fout(*this);
    uint32_t ind = this->block_create_routine(key, fout.exitcode);
    darray[ind].d = newelem;
}

LFHTEMP bool AsyncHashmap<K,D,HF>::doesContain(const K &key) const{
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
return this->find_routine(key)!= 0xFFFFFFFF;}
LFHTEMP typename AsyncHashmap<K,D,HF>::ReadAccess AsyncHashmap<K,D,HF>::tryRead(const K &key) const{ typename AsyncHashmap<K,D,HF>::ReadAccess fout(*this);
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    uint32_t ind = fout.find_read_syncroutine(key);
return fout;}
LFHTEMP typename AsyncHashmap<K,D,HF>::WriteAccess AsyncHashmap<K,D,HF>::tryWrite(const K &key){ typename AsyncHashmap<K,D,HF>::WriteAccess fout(*this);
    uint32_t ind = fout.try_write_routine(key);
    fout.elem = (fout.exitcode == 0) ? nullptr : darray + ind;
return fout;}
LFHTEMP typename AsyncHashmap<K,D,HF>::WriteAccess AsyncHashmap<K,D,HF>::tryCreate(const K &key){ typename AsyncHashmap<K,D,HF>::WriteAccess fout(*this);
    uint32_t ind = this->try_create_routine(key, fout.exitcode);
    fout.elem = (fout.exitcode == 0) ? nullptr : darray + ind;
return fout;}

/*
LFHTEMP typename AsyncHashmap<K,D,HF>::WriteAccess  AsyncHashmap<K,D,HF>::addEntry(K &key){
    return AsyncHashmap<K,D,HF>::WriteAccess(*this);
}
LFHTEMP void AsyncHashmap<K,D,HF>::removeEntry(K &key){
}*/
LFHTEMP void AsyncHashmap<K,D,HF>::show(FILE *f, int level)const{
    printf("Async Hashmap with %i items\n", this->getSize());
    if (asize == 0) return;
    /*if (auto ite = this->getIterator()) do{
        ExOp::show(*ite, f, level);
    }while(ite++);*/
    uint32_t j;
    for(int i=0;i < hsize;i++){
        uint32_t j = dahash[i | hsize];
        printf("\n%i:", i);
        while(j != 0xFFFFFFFF) { ExOp::show(darray[j], f, 1); j = dahash[j];}
    }
}

LFHTEMP void AsyncHashmap<K,D,HF>::testThisDataStructure(uint32_t nbthreads, uint32_t nbitems, uint32_t noisythread){
    class crazyTask{
        public:
        LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> > &target;
        uint32_t nbitems;
        uint32_t mask;
        uint32_t donetask;
        uint32_t nbnastythreads;
        uint32_t iteratortrap;
        uint32_t noisythread;
        crazyTask(LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> >& _target, uint32_t _nbitems): target(_target), nbitems(_nbitems){}
        uint32_t operator()(uint32_t threadID){
            myHashmap<uint32_t> dabuffer;
            uint32_t r,fail;
            bool doit;
            for(fail =nbitems+1; fail != nbitems;){
                if ((iteratortrap == threadID) && ((rand() & 255) == 0)){
                    r =0;
                    uint32_t sum =0;
                    if (auto ite = target.getIterator()) do{
                        r++;
                        sum += ite->d;
                        if (ite() > nbitems * 7) foreach.fprintf(stdout, "strange %i!!!\n", ite());
                        if (((rand() & 1023) == 0)&&(!target.doesContain(ite()))) foreach.fprintf(stdout, "failed to reach %i!!!\n", ite());
                    }while(ite++);
                    if (r < 16){
                        if (auto ite = target.getIterator()) do{
                            foreach.fprintf(stdout, "T[%i]: F[%i] = %i\n",iteratortrap , ite(), ite->d);
                        }while(ite++);
                    }else foreach.fprintf(stdout, "T[%i] = %i (%i items) (curfail %i)\n",iteratortrap , sum, r, fail);
                    iteratortrap = (iteratortrap + 1) % nbnastythreads;
                }
                r = (fail > nbitems) ? (rand() % nbitems) : fail;
                r = r * 7;
                if ((noisythread == threadID)&&((rand() & 255) == 0)) foreach.fprintf(stdout,"%i r is %i %i\n", fail, r, (r - threadID) & mask);
                if (dabuffer.find(r) == 0xFFFFFFFF){
                    if (((r ^ threadID ) & mask) == 0){
                        if (threadID <= mask) { // in charge of creation
                            if ((rand() & 15) == 0){
                                target[r]->d = threadID+1;
                                dabuffer.addEntry(r); foreach.updateProgress(threadID);
                            }else{
                                if (auto danew = target.tryCreate(r)) {
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                    danew->d = threadID+1;
                                    //std::this_thread::sleep_for(std::chrono::seconds((rand() & 3) + 1)); // hold it for fun ^^
                                }else if (noisythread == threadID) foreach.fprintf(stdout," could create %i!\n", r);
                            }
                        }else{ // in charge of deletion
                            if (auto daacc = target.tryWrite(r)) {
                                 if (daacc->d >= ((threadID * (threadID + 1)) >> 1)) {
                                    daacc.remove();
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                 }else{
                                    daacc.unlock();
                                    //if (fail <= nbitems) std::this_thread::sleep_for(std::chrono::seconds(1)); // give a chance to other threads if scanning
                                 }
                            }
                        }
                    }else{
                        if ((threadID > mask)&&((threadID & mask) > (r & mask))) {dabuffer.addEntry(r); foreach.updateProgress(threadID); continue;} // do not mind this one
                        if ((rand() & 15) == 0){
                            if (target.doesContain(r)){
                                doit = (target.read(r)->d >= ((threadID * (threadID + 1)) >> 1));
                                if (doit) {
                                    target[r]->d += threadID+1;
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                }else{
                                    //if (fail <= nbitems) std::this_thread::sleep_for(std::chrono::seconds(1)); // give a chance to other threads if scanning
                                    if (noisythread == threadID) foreach.fprintf(stdout," could not update, f{%i] = %i is val!\n", r, target[r]->d );
                                }
                            }
                        }else{
                            if (auto danew = target.tryRead(r)) {
                                doit = (danew->d >= ((threadID * (threadID + 1)) >> 1));
                                if ((!doit)&&(noisythread == threadID)) foreach.fprintf(stdout," could not update try, f{%i] = %i is val!\n", r, danew->d );
                                danew.unlock();
                                //if ((!doit)&&(fail <= nbitems)) std::this_thread::sleep_for(std::chrono::seconds(1));
                            }else doit = false;

                            if (doit){
                                if (auto daacc = target.tryWrite(r)){
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                    daacc->d += threadID +1;

                                }
                            }// give a chance to other threads if scanning
                        }
                    }

                   // if ((fail < nbitems)&&(threadID == 0))  printf("%i\n", fail);
                }else {
                    //if ((rand() & 1023) == 0) std::this_thread::sleep_for(std::chrono::seconds(1));
                    fail++;
                    if (fail >= nbitems * 10) {
                        fail = 0;
                        foreach.fprintf(stdout," is now scanning to finish!!!\n", threadID);
                        /*if (auto ite = target.getIterator()) do{
                            if ((rand() & 255) == 0) printf("F[%i] = %i\n", ite(), ite->d);
                        }while(ite++);*/
                    }
                }
            }
            foreach.finishProgress(threadID);
            foreach.fprintf(stdout,"\n%i %X exit!\n", threadID, target.queryChunkIO());
           /* if (auto ite = target.getIterator()) do{
                foreach.printf("F[%i] = %i\n", ite(), ite->d);
            }while(ite++);
            foreach.printf("\n%i %X\n", threadID, target.queryChunkIO());*/
        return 0;}

        void iterMadness(){
            /*for(int i =0;i<10;i++){
                if (auto ite = target.getIterator()) do{
                    if ((rand() & 255) == 0) printf("F[%i] = %i\n", ite(), ite->d);
                }while(ite++);
                printf("sleepy\n");
                std::this_thread::sleep_for(std::chrono::seconds(2));
            }*/
            printf("im done\n");
        }
        std::function<void()> fncIterMadness(){return std::bind(&crazyTask::iterMadness,*this);}
    };


    LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> > test;
    crazyTask lotask(test, nbitems);
    lotask.nbnastythreads = nbthreads;
    lotask.iteratortrap =0;
    lotask.noisythread = noisythread;
    lotask.mask = (lotask.nbnastythreads >> 1)-1; ExOp::toLeftFlood(lotask.mask);
    foreach.startThreadArray(nbthreads);
    printf("Doing this with %i threads (mask %X)\n",foreach.getThreadArraySize(), lotask.mask);

    foreach.startProgress("Lets do this!", lotask.nbnastythreads * lotask.nbitems);
    for(int i = lotask.nbnastythreads-1 ; i>0;i--) {
        foreach.submit(lotask, i);
        //std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    foreach.submit_ThenWait(lotask, 0);
   // foreach.submit(lotask.fncIterMadness());
}



#undef LFHTEMP
#define LFHTEMP template <class C, class HF>
LFHTEMP void AsyncHashmap<C,void,HF>::rehash(uint32_t _hsize){
    if (asize != 0) delete[](dahash);
    hsize = _hsize;
    if (_hsize == 0) return;
    dahash = new uint32_t[hsize << 1];
    memset(dahash + hsize , '\xFF', hsize << 2);
    uint32_t i;
    for(i=0;i< asize;i++){
        uint32_t j = (HF::makeSeed(ExOp::getIndex(darray[i])) & (hsize-1)) | hsize;
        dahash[i] = dahash[j];
        dahash[j] = i;
    }
    /*for(i=0;i< (hsize << 1);i++){
        printf("H[%i]: %i\n",i, dahash[i]);

    }*/
}
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::find_routine(const typename ExCo<C>::INDEX_TYPE &key) const{
    if (asize == 0) return(0xFFFFFFFF);
    unsigned int j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
    for(;j != 0xFFFFFFFF;j = dahash[j]){
        if (key == ExOp::getIndex(darray[j])) break;
    }
return(j);}
LFHTEMP void AsyncHashmap<C,void,HF>::delayed_remove_routine(const typename ExCo<C>::INDEX_TYPE &key, std::unique_lock<std::mutex> &lck){
    // just kidding, not deleting ^^'
    //printf("delayed wait!\n"); fflush(stdout);
    uint32_t ind;
    foreach.accessFinalCondvar().wait(lck,[this, key, &ind]{
            if (asize == 0) {ind = 0xFFFFFFFF; return true;}
            if ((chunkIO & (3 << (((hsize >> 5) == 0) ?  ((asize-1) << 1) : (((asize-1) / (hsize >> 5)) & 0x1E)))) != 0) return false;
            uint32_t i = (HF::makeSeed(key) & (hsize-1)) | hsize;
            for(ind = dahash[i];;ind= dahash[ind]){
                if (ind == 0xFFFFFFFF) return true;
                if (key == ExOp::getIndex(darray[ind])) break;
            }
            return ((chunkIO & (3 << (((hsize >> 5) == 0) ?  (ind << 1) : ((ind / (hsize >> 5)) & 0x1E)))) == 0);
        });
    if (ind == 0xFFFFFFFF) return;
    if (asize == 1){ // deleting the last one!
            //printf("delay del last!\n"); fflush(stdout);
        delete[](dahash);
        delete[](darray);
        asize =0;
    }else{
        //printf("delay del!"); fflush(stdout);
        if (key != ExOp::getIndex(darray[ind])){
            printf("crazy, did not match!\n");
        }
        uint32_t i = (HF::makeSeed(key) & (hsize-1)) | hsize;
        uint32_t j = dahash[i];
        if (j == ind) dahash[i] = dahash[j];
        else{
            while (dahash[j] != ind) j = dahash[j];
            dahash[j] = dahash[ind];
        }
        moveEntry_routine(asize-1, ind);
        asize--;

        if (asize <= (hsize >> 2)){
            //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO==0)||(asize > (hsize >> 2));});
            //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)chunkIO);
            if (asize <= (hsize >> 2)){
                //printf("delay dadel downalloc!\n"); fflush(stdout);
                C* swp = darray;
                darray = new C[hsize >> 1];
                for(i=0;i< asize; i++) ExOp::toMemmove(darray[i], swp[i]);
                rehash(hsize >> 1);
                delete[](swp);
            }
        }
    }
return;}
LFHTEMP void AsyncHashmap<C,void,HF>::moveEntry_routine(uint32_t from, uint32_t to){
    if (from == to) return;
    uint32_t i;
    uint32_t j = dahash[i = ((HF::makeSeed(ExOp::getIndex(darray[from])) & (hsize-1)) | hsize)];
    if (j == from) dahash[i] = to;
    else{
        while (dahash[j] != from) j = dahash[j];
        dahash[j] = to;
    }
    dahash[to] = dahash[from];
    darray[to] = std::move(darray[from]);
}
LFHTEMP AsyncHashmap<C,void,HF>::WriteAccess::~WriteAccess(){
	uint32_t tmp;
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); tmp = target.chunkIO.fetch_sub(exitcode); if ((tmp & exitcode) != exitcode) {printf("freed %X but %X is atomic previous value...\n", exitcode, tmp); exit(1);} foreach.accessFinalCondvar().notify_all();}
}
LFHTEMP void AsyncHashmap<C,void,HF>::WriteAccess::unlock(){
	uint32_t tmp;
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); tmp = target.chunkIO.fetch_sub(exitcode); if ((tmp & exitcode) != exitcode) {printf("freed %X but %X is atomic previous value...\n", exitcode, tmp); exit(1);} foreach.accessFinalCondvar().notify_all(); exitcode =0;}
}
LFHTEMP AsyncHashmap<C,void,HF>::ReadAccess::~ReadAccess(){
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());  target.chunkIO.fetch_sub(exitcode); foreach.accessFinalCondvar().notify_all();}
}
LFHTEMP void AsyncHashmap<C,void,HF>::ReadAccess::unlock(){
    if (exitcode != 0){std::lock_guard<std::mutex> lck(foreach.accessFinalMutex()); target.chunkIO.fetch_sub(exitcode); foreach.accessFinalCondvar().notify_all(); exitcode =0;}
}
/*LFHTEMP bool AsyncHashmap<C,void,HF>::ReadAccess::find_only_routine(const typename ExCo<C>::INDEX_TYPE &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0; return false;}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == ExOp::getIndex(target.darray[j])) break;
    }
    exitcode =0;
return(j != 0xFFFFFFFF);}*/
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::ReadAccess::find_read_syncroutine(const typename ExCo<C>::INDEX_TYPE &key){
    if (target.asize == 0) {exitcode = 0; return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == ExOp::getIndex(target.darray[j])) break;
    }

    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        uint32_t tmp = (exitcode | (exitcode << 1));
        if ((target.chunkIO & tmp) == tmp) exitcode = 0;
        else {
            target.chunkIO.fetch_add(exitcode);
            elem = target.darray + j;
        }
    }
return(j);}
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::WriteAccess::try_write_routine(const typename ExCo<C>::INDEX_TYPE &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == ExOp::getIndex(target.darray[j])) break;
    }
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 3 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        if ((target.chunkIO & exitcode) == 0) target.chunkIO.fetch_or(exitcode);
        else exitcode = 0;
    }
return(j);}
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::ReadAccess::block_routine(const typename ExCo<C>::INDEX_TYPE &key){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]){
        if (key == ExOp::getIndex(target.darray[j])) break;
    }
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        uint32_t tmp = (exitcode | (exitcode << 1));

        //if (((target.chunkIO | target.angryIO) & tmp) == tmp) {
            //foreach.fprintf_l(stdout,"Deadlock prone want %X got chunkIO %X\n", tmp, (uint32_t)target.chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this,tmp]{return ((target.chunkIO | target.angryIO) & tmp) != tmp;});
            //foreach.fprintf_l(stdout,"Woke up! got everything chunkIO %X\n", (uint32_t)target.chunkIO);
        //}

        target.chunkIO.fetch_add(exitcode);
    }
return(j);}
/* Delete entry from hashmap, needs to make sure no other reader exist
 */
LFHTEMP void AsyncHashmap<C,void,HF>::ReadAccess::remove(){
     if (exitcode != 0){
        // unlink from hash
        {std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
            uint32_t i;
            if (target.chunkIO & (exitcode << 1)){
                // cant write...
                target.chunkIO.fetch_sub(exitcode); exitcode =0;
                foreach.accessFinalCondvar().notify_all();
                target.delayed_remove_routine(ExOp::getIndex(*elem), lck);
                elem = NULL;
                return;
            }

            // first try to swap/delete entry

            if (target.asize == 1){
                delete[](target.dahash);
                delete[](target.darray);
                target.asize =0;
                //printf("dadel del last! %i %p %p L=%c\n", target.asize, target.dahash, target.darray, lck ? 'L' : 'U'); fflush(stdout);
            }else{
                // needs to lock entry and end of hash!
                uint32_t endcode = 1 << (((target.hsize >> 5) == 0) ?  ((target.asize-1) << 1) : (((target.asize-1) / (target.hsize >> 5)) & 0x1E));
                if ((endcode == exitcode)||((target.chunkIO & (endcode | (endcode << 1)) ) == 0)) {
                    uint32_t ind = (uint32_t) (elem - target.darray);
                    //printf("dadel read! %i\n", ind); fflush(stdout);
                    if (ind != target.find_routine(ExOp::getIndex(*elem))) printf("Disagree for position!\n", ind, target.find_routine(ExOp::getIndex(*elem)));
                    uint32_t j = target.dahash[i = ((HF::makeSeed(ExOp::getIndex(*elem)) & (target.hsize-1)) | target.hsize)];
                    if (j == ind) target.dahash[i] = target.dahash[j];
                    else{
                        while (target.dahash[j] != ind) j = target.dahash[j];
                        target.dahash[j] = target.dahash[ind];
                    }
                    target.moveEntry_routine(target.asize-1, ind);
                    target.asize--;
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                }else{ // cannot access end... need free it and to wait for it  and end of vector, implying remove is postponed!
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                    target.delayed_remove_routine(ExOp::getIndex(*elem), lck);
                    elem = NULL;
                    return;
                }
                // cant exit if downalloc is needed
                if (target.asize <= (target.hsize >> 2)){
                    //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) target.chunkIO);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO==0)||(target.asize > (target.hsize >> 2));});
                    //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
                    if (target.asize <= (target.hsize >> 2)){
                        //printf("dadel downalloc!\n"); fflush(stdout);
                        C* swp = target.darray;
                        target.darray = new C[target.hsize >> 1];
                        for(i=0;i< target.asize; i++) target.darray[i] = std::move(swp[i]);
                        target.rehash(target.hsize >> 1);
                        delete[](swp);
                    }
                }
            }
            elem = NULL;
        }
    }
}
/* Delete entry from hashmap
 */
LFHTEMP void AsyncHashmap<C,void,HF>::WriteAccess::remove(){
    if (exitcode != 0){
        // unlink from hash
        {std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
            uint32_t i;

            // first try to swap/delete entry

            if (target.asize == 1){
                delete[](target.dahash);
                delete[](target.darray);
                target.asize =0;
            }else{
                // needs to lock entry and end of hash!
                uint32_t endcode = 3 << (((target.hsize >> 5) == 0) ?  ((target.asize-1) << 1) : (((target.asize-1) / (target.hsize >> 5)) & 0x1E));
                if ((endcode == exitcode)||((target.chunkIO & endcode) == 0)) {
                    uint32_t ind = (uint32_t) (elem - target.darray);
                    //printf("dadel write! %i\n", ind);
                    if (ind != target.find_routine(ExOp::getIndex(*elem))) printf("Disagree for position!\n", ind, target.find_routine(ExOp::getIndex(*elem)));
                    uint32_t j = target.dahash[i = ((HF::makeSeed(ExOp::getIndex(*elem)) & (target.hsize-1)) | target.hsize)];
                    if (j == ind) target.dahash[i] = target.dahash[j];
                    else{
                        while (target.dahash[j] != ind) j = target.dahash[j];
                        target.dahash[j] = target.dahash[ind];
                    }
                    target.moveEntry_routine(target.asize-1, ind);
                    target.asize--;
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                }else{ // cannot access end... need free it and to wait for it  and end of vector, implying remove is postponed!
                    target.chunkIO.fetch_sub(exitcode); exitcode =0;
                    foreach.accessFinalCondvar().notify_all();
                    target.delayed_remove_routine(ExOp::getIndex(*elem), lck);
                    elem = NULL;
                    return;
                }
                // cant exit if downalloc is needed
                if (target.asize <= (target.hsize >> 2)){
                    //foreach.fprintf_l(stdout,"Wants end (downalloc) chunkIO %X\n", (uint32_t) target.chunkIO);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO==0)||(target.asize > (target.hsize >> 2));});
                    //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
                    if (target.asize <= (target.hsize >> 2)){
                        C* swp = target.darray;
                        target.darray = new C[target.hsize >> 1];
                        for(i=0;i< target.asize; i++) ExOp::toMemmove(target.darray[i], swp[i]);
                        target.rehash(target.hsize >> 1);
                        delete[](swp);
                    }
                }
            }
            elem = NULL;
        }
    }
}
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::WriteAccess::block_routine(const typename ExCo<C>::INDEX_TYPE &key){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    if (target.asize == 0) {exitcode = 0;return(0xFFFFFFFF);}
    uint32_t j = target.dahash[(HF::makeSeed(key) & (target.hsize-1)) | target.hsize];
    for(;j != 0xFFFFFFFF;j = target.dahash[j]) if (key == ExOp::getIndex(target.darray[j])) break;
    if (j == 0xFFFFFFFF) exitcode =0;
    else{
        exitcode = 1 << (((target.hsize >> 5) == 0) ?  (j << 1) : ((j / (target.hsize >> 5)) & 0x1E));
        exitcode |= exitcode << 1;

        if ((target.chunkIO & exitcode) != 0) {
            // trying to stay calm ^_^
            //if ((target.angryIO & exitcode) != 0) foreach.accessFinalCondvar().wait(lck,[this]{return (target.angryIO & exitcode) == 0;});
            //target.angryIO |=  exitcode; // got angry!
            //if ((target.chunkIO & exitcode) != 0) {
               // foreach.fprintf_l(stdout,"writing block %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
                foreach.accessFinalCondvar().wait(lck,[this]{return (target.chunkIO & exitcode) == 0;});
            //}
            //target.angryIO ^=  exitcode; // is happy now!
        }
        if ((target.chunkIO.fetch_or(exitcode) & exitcode) != 0){
			printf("Impossible!!! required access without holding mutex!\n");
			exit(1);
        }
    }
return(j);}

LFHTEMP AsyncHashmap<C,void,HF>::ConstIterator::~ConstIterator(){
    if (exitcode != 0){
        target.chunkIO.fetch_sub(exitcode);
        foreach.accessFinalCondvar().notify_all();
    }
}
LFHTEMP AsyncHashmap<C,void,HF>::ConstIterator::operator  bool (){
    if (target.asize == 0) return false;
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    // locks first and last (no delete and insert allowed, read/write can occur)
    exitcode = 1 | (0x55555555 << ((((target.hsize >> 5) == 0) ? target.asize - 1 : target.asize / (target.hsize >> 5)) & 0x1E));
   // if (exitcode & (target.angryIO |(target.chunkIO & (target.chunkIO >>1))) != 0){
        //foreach.fprintf_l(stdout,"iterator start %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
        foreach.accessFinalCondvar().wait(lck,[this]{return (exitcode & (target.angryIO | (target.chunkIO & (target.chunkIO >>1))) == 0);});
        //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
   //     target.chunkIO.fetch_add(exitcode);
   // }else
    target.chunkIO.fetch_add(exitcode);
    i = 0;
return true;}
LFHTEMP bool AsyncHashmap<C,void,HF>::ConstIterator::operator++(int){
    if (++i == target.asize){
        target.chunkIO.fetch_sub(exitcode); exitcode = 0;
        foreach.accessFinalCondvar().notify_all();
        return false;
    }
    uint32_t tmp = 1 << (((target.hsize >> 5) == 0) ?  (i << 1) : ((i / (target.hsize >> 5)) & 0x1E));
    if ((tmp & exitcode) == 0){
        target.chunkIO.fetch_sub(tmp >> 2);
        uint32_t  tmp2 = tmp | (tmp << 1);
        std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());

        if (((target.angryIO | target.chunkIO) & tmp2) == tmp2){
            foreach.accessFinalCondvar().notify_all();
            //foreach.fprintf_l(stdout,"Iterator move wants %X got chunkIO %X\n", tmp2, (uint32_t)target.chunkIO);
            foreach.accessFinalCondvar().wait(lck,[this,tmp2]{return (((target.angryIO | target.chunkIO) & tmp2) != tmp2);});
            //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t)target.chunkIO);
            target.chunkIO.fetch_add(tmp);
        }else{
            target.chunkIO.fetch_add(tmp);
            foreach.accessFinalCondvar().notify_all();
        }
        exitcode ^= tmp | (tmp >> 2);
    }
return true;}
LFHTEMP void AsyncHashmap<C,void,HF>::remove(const typename ExCo<C>::INDEX_TYPE &key){
    // needs to lock entry and end of array!
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
return delayed_remove_routine(key, lck);}
LFHTEMP uint32_t AsyncHashmap<C,void,HF>::try_create_routine(const typename ExCo<C>::INDEX_TYPE &key, uint32_t &exitcode){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    uint32_t j;
    if (asize == 0) {
        hsize = 1;
        darray = new C[1];
        dahash = new uint32_t[2];
        chunkIO = 3;
        angryIO = 0;
        exitcode = 3;
        dahash[0] = 0xFFFFFFFF;
        dahash[1] = 0;
        return asize++;
    }else if (asize == hsize){
        // need to upalloc!!!
        if (chunkIO != 0) {exitcode = 0; return 0xFFFFFFFF;} // no can do...
        this->rehash(hsize << 1);
        C* swp = darray;
        darray = new C[hsize];
        for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
        delete[](swp);
    }else j = asize;
    exitcode = 3 << (((hsize >> 5) == 0) ?  (j << 1) : ((j / (hsize >> 5)) & 0x1E));
    if ((chunkIO & exitcode) != 0){exitcode = 0; return j;}
    chunkIO |= exitcode;
    ExOp::getIndex(darray[j]) = key;
    j = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[asize] = dahash[j];
    dahash[j] = asize;
    asize++;
return asize;}

LFHTEMP C& AsyncHashmap<C,void,HF>::createEntry_sync(const typename ExCo<C>::INDEX_TYPE key){
    uint32_t j;
	if (asize == 0) {
		hsize = 1;
		dahash = new uint32_t[2];
		memset(dahash + hsize , '\xFF', 4);
		this->rehash(1);
		darray = new C[1];
		chunkIO = 0;
	}else if (asize == hsize){ // need to upalloc!!!
		this->rehash(hsize << 1);
		C* swp = darray;
		darray = new C[hsize];
		for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
		delete[](swp);
	}
    j = asize++;
    ExOp::getIndex(darray[j]) = key;
    uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[j] = dahash[tj];
    dahash[tj] = j;
return darray[j];}

LFHTEMP uint32_t AsyncHashmap<C,void,HF>::find_sync(const typename ExCo<C>::INDEX_TYPE key)const{
	if (asize == 0) return(0xFFFFFFFF);
    uint32_t j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
    for(;j != 0xFFFFFFFF;j = dahash[j]){
        if (key == ExOp::getIndex(darray[j])) break;
    }
return j;}


LFHTEMP uint32_t AsyncHashmap<C,void,HF>::block_create_routine(const typename ExCo<C>::INDEX_TYPE &key, uint32_t &exitcode){
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    printf("block_create_routine...\n");
    uint32_t j;
    do{
        if (asize == 0) {
            //printf("created 1! %c\n", lck ? 'L' : 'U'); fflush(stdout);
            //fflush(stdout);
            this->rehash(1);
            darray = new C[1];
            chunkIO = 0;
            exitcode = 3;
            break;
        }else if (asize == hsize){
            // need to upalloc!!!
            //printf("upalloc time! %c\n", lck ? 'L' : 'U'); fflush(stdout);
            do{
                if (chunkIO != 0) { // todo, make angry
                    //foreach.fprintf_l(stdout,"Deadlockprone want everything got chunkIO %X\n", (uint32_t)chunkIO);
                    //printf("sleep time! %c\n", lck ? 'L' : 'U'); fflush(stdout);
                    foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO==0)||(asize !=  hsize);});
                    //foreach.fprintf_l(stdout,"Woke up! got everything chunkIO %X\n", (uint32_t)chunkIO);
                    if (asize < hsize) {
                        //printf("something happened! %i\n", asize);
                        if (asize == 0){
                            this->rehash(1);
                            darray = new C[1];
                            chunkIO = 0;
                        }
                        break;
                    } // some delete or up alloc happened
                }
                //printf("all good! %i %i %c\n", asize, hsize, lck ? 'L' : 'U'); fflush(stdout);
                this->rehash(hsize << 1);
                C* swp = darray;
                darray = new C[hsize];
                for(j=0;j< asize;j++) darray[j] = std::move(swp[j]);
                delete[](swp);
                j = asize++;
                exitcode = 3 << (((hsize >> 5) == 0) ?  (j << 1) : ((j / (hsize >> 5)) & 0x1E));
                chunkIO |= exitcode;

                ExOp::getIndex(darray[j]) = key;
                uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
                dahash[j] = dahash[tj];
                dahash[tj] = j;
                return j;
            }while(false);

        }
        exitcode = 3 << (((hsize >> 5) == 0) ?  (asize << 1) : ((asize / (hsize >> 5)) & 0x1E));
        if ((chunkIO & exitcode) == 0) break;
                    // trying to stay calm ^_^
            //if ((angryIO & exitcode) != 0) foreach.accessFinalCondvar().wait(lck,[this,exitcode]{return (angryIO & exitcode) == 0;});
            //angryIO |=  exitcode; // got angry!
        //    if ((chunkIO & exitcode) != 0) {
               // foreach.fprintf_l(stdout,"writing block %X got chunkIO %X\n", exitcode, (uint32_t)target.chunkIO);
                foreach.accessFinalCondvar().wait(lck,[this]{return (chunkIO & (3 << (((hsize >> 5) == 0) ?  (asize << 1) : ((asize / (hsize >> 5)) & 0x1E)))) == 0;});

                //foreach.fprintf_l(stdout,"Woke up! chunkIO %X\n", (uint32_t) target.chunkIO);
           //}
            //angryIO ^=  exitcode; // is happy now!
        //}
    }while(true);
    j = asize++;
    chunkIO |= exitcode;
    ExOp::getIndex(darray[j]) = key;
    uint32_t tj = (HF::makeSeed(key) & (hsize-1)) | hsize;
    dahash[j] = dahash[tj];
    dahash[tj] = j;
    /*printf("created exity! %i %i %i %i %c\n", asize, hsize, tj, j, lck ? 'L' : 'U'); fflush(stdout);
    for(int i=0;i< (hsize << 1);i++){
        printf("H[%i]: %i\n",i, dahash[i]);
    }*/
    fflush(stdout);
return j;}
/*LFHTEMP uint32_t AsyncHashmap<C,void,HF>::find_create_routine(const typename ExCo<C>::INDEX_TYPE &key){
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    unsigned int j;
    if (asize == 0) {
        j = dahash[(HF::makeSeed(key) & (hsize-1)) | hsize];
        for(;j != 0xFFFFFFFF;j = dahash[j]){
            if (key == ExOp::getIndex(darray[j]) ) return(j);
        }
    }

    uint32_t rsize = (asize & 0x7FFFFFFF);
    j = rsize;
    if ((((asize + 1u)^(asize))> rsize)&&((asize & 0x80000000) == 0)) { // going to up_alloc
        // TODO
       // heap.push_back();
       // heap[j].k.k = key;
       // rehash(hash_mag + 1);
    }else{
        uint32_t j = (HF::makeSeed(key) & (hsize-1)) | hsize;
        dahash[asize] = dahash[j];
        dahash[j] = asize;
    }
return asize++;}*/

LFHTEMP typename AsyncHashmap<C,void,HF>::WriteAccess AsyncHashmap<C,void,HF>::operator[](const typename ExCo<C>::INDEX_TYPE &key){ typename AsyncHashmap<C,void,HF>::WriteAccess fout(*this);

    // wait... that hash can be thrashed *but* the offset should be allowed to change... need to lock that first
    uint32_t ind = fout.block_routine(key);
    if (ind == 0xFFFFFFFF) {
    	printf("I thought this is impossible now!\n");
		exit(0);
		ind = this->block_create_routine(key, fout.exitcode);
	}
    //printf("daind%i%i%i%idniad\n", ind, ind, ind, ind);
    fout.elem = darray + ind;
return fout;}
LFHTEMP typename AsyncHashmap<C,void,HF>::ReadAccess AsyncHashmap<C,void,HF>::read(const typename ExCo<C>::INDEX_TYPE &key)const{ typename AsyncHashmap<C,void,HF>::ReadAccess fout(*this);
    uint32_t ind = fout.block_routine(key);
    if (ind != 0xFFFFFFFF) fout.elem = darray + ind;
return fout;}

LFHTEMP typename AsyncHashmap<C,void,HF>::ReadAccess AsyncHashmap<C,void,HF>::await(const typename ExCo<C>::INDEX_TYPE &key, atomic<int32_t>* can_block_flag) const{typename AsyncHashmap<C,void,HF>::ReadAccess fout(*this);
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
    fout.find_read_syncroutine(key);
    if (fout.exitcode != 0) return fout;
    // does not exist, wait for it!
    if (0 == (int32_t) *can_block_flag) foreach.accessFinalCondvar().wait(lck,[&fout, &key, can_block_flag]{fout.find_read_syncroutine(key);return (fout.exitcode != 0)||(0 != (int32_t) *can_block_flag);});
return fout;}

LFHTEMP void AsyncHashmap<C,void,HF>::insert(C && newelem){ typename AsyncHashmap<C,void,HF>::WriteAccess fout(*this);
    uint32_t ind = this->block_create_routine(ExOp::getIndex(newelem), fout.exitcode);
    darray[ind] = newelem;
}

LFHTEMP bool AsyncHashmap<C,void,HF>::doesContain(const typename ExCo<C>::INDEX_TYPE &key) const{
    std::unique_lock<std::mutex> lck(foreach.accessFinalMutex());
return this->find_routine(key)!= 0xFFFFFFFF;}
LFHTEMP typename AsyncHashmap<C,void,HF>::ReadAccess AsyncHashmap<C,void,HF>::tryRead(const typename ExCo<C>::INDEX_TYPE &key) const{ typename AsyncHashmap<C,void,HF>::ReadAccess fout(*this);
    std::lock_guard<std::mutex> lck(foreach.accessFinalMutex());
    uint32_t ind = fout.find_read_syncroutine(key);
return fout;}
LFHTEMP typename AsyncHashmap<C,void,HF>::WriteAccess AsyncHashmap<C,void,HF>::tryWrite(const typename ExCo<C>::INDEX_TYPE &key){ typename AsyncHashmap<C,void,HF>::WriteAccess fout(*this);
    uint32_t ind = fout.try_write_routine(key);
    fout.elem = (fout.exitcode == 0) ? nullptr : darray + ind;
return fout;}
LFHTEMP typename AsyncHashmap<C,void,HF>::WriteAccess AsyncHashmap<C,void,HF>::tryCreate(const typename ExCo<C>::INDEX_TYPE &key){ typename AsyncHashmap<C,void,HF>::WriteAccess fout(*this);
    uint32_t ind = this->try_create_routine(key, fout.exitcode);
    fout.elem = (fout.exitcode == 0) ? nullptr : darray + ind;
return fout;}

/*
LFHTEMP typename AsyncHashmap<C,void,HF>::WriteAccess  AsyncHashmap<C,void,HF>::addEntry(typename ExCo<C>::INDEX_TYPE &key){
    return AsyncHashmap<C,void,HF>::WriteAccess(*this);
}
LFHTEMP void AsyncHashmap<C,void,HF>::removeEntry(typename ExCo<C>::INDEX_TYPE &key){
}*/
LFHTEMP void AsyncHashmap<C,void,HF>::show(FILE *f, int level)const{
    printf("Async Hashmap with %i items\n", this->getSize());
    if (asize == 0) return;
    /*if (auto ite = this->getIterator()) do{
        ExOp::show(*ite, f, level);
    }while(ite++);*/
    uint32_t j;
    for(int i=0;i < hsize;i++){
        uint32_t j = dahash[i | hsize];
        printf("\n%i:", i);
        while(j != 0xFFFFFFFF) { ExOp::show(darray[j], f, 1); j = dahash[j];}
    }
}

LFHTEMP void AsyncHashmap<C,void,HF>::testThisDataStructure(uint32_t nbthreads, uint32_t nbitems, uint32_t noisythread){
    class crazyTask{
        public:
        LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> > &target;
        uint32_t nbitems;
        uint32_t mask;
        uint32_t donetask;
        uint32_t nbnastythreads;
        uint32_t iteratortrap;
        uint32_t noisythread;
        crazyTask(LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> >& _target, uint32_t _nbitems): target(_target), nbitems(_nbitems){}
        uint32_t operator()(uint32_t threadID){
            myHashmap<uint32_t> dabuffer;
            uint32_t r,fail;
            bool doit;
            for(fail =nbitems+1; fail != nbitems;){
                if ((iteratortrap == threadID) && ((rand() & 255) == 0)){
                    r =0;
                    uint32_t sum =0;
                    if (auto ite = target.getIterator()) do{
                        r++;
                        sum += ite->d;
                        if (ite() > nbitems * 7) foreach.fprintf(stdout, "strange %i!!!\n", ite());
                        if (((rand() & 1023) == 0)&&(!target.doesContain(ite()))) foreach.fprintf(stdout, "failed to reach %i!!!\n", ite());
                    }while(ite++);
                    if (r < 16){
                        if (auto ite = target.getIterator()) do{
                            foreach.fprintf(stdout, "T[%i]: F[%i] = %i\n",iteratortrap , ite(), ite->d);
                        }while(ite++);
                    }else foreach.fprintf(stdout, "T[%i] = %i (%i items) (curfail %i)\n",iteratortrap , sum, r, fail);
                    iteratortrap = (iteratortrap + 1) % nbnastythreads;
                }
                r = (fail > nbitems) ? (rand() % nbitems) : fail;
                r = r * 7;
                if ((noisythread == threadID)&&((rand() & 255) == 0)) foreach.fprintf(stdout,"%i r is %i %i\n", fail, r, (r - threadID) & mask);
                if (dabuffer.find(r) == 0xFFFFFFFF){
                    if (((r ^ threadID ) & mask) == 0){
                        if (threadID <= mask) { // in charge of creation
                            if ((rand() & 15) == 0){
                                target[r]->d = threadID+1;
                                dabuffer.addEntry(r); foreach.updateProgress(threadID);
                            }else{
                                if (auto danew = target.tryCreate(r)) {
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                    danew->d = threadID+1;
                                    //std::this_thread::sleep_for(std::chrono::seconds((rand() & 3) + 1)); // hold it for fun ^^
                                }else if (noisythread == threadID) foreach.fprintf(stdout," could create %i!\n", r);
                            }
                        }else{ // in charge of deletion
                            if (auto daacc = target.tryWrite(r)) {
                                 if (daacc->d >= ((threadID * (threadID + 1)) >> 1)) {
                                    daacc.remove();
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                 }else{
                                    daacc.unlock();
                                    //if (fail <= nbitems) std::this_thread::sleep_for(std::chrono::seconds(1)); // give a chance to other threads if scanning
                                 }
                            }
                        }
                    }else{
                        if ((threadID > mask)&&((threadID & mask) > (r & mask))) {dabuffer.addEntry(r); foreach.updateProgress(threadID); continue;} // do not mind this one
                        if ((rand() & 15) == 0){
                            if (target.doesContain(r)){
                                doit = (target.read(r)->d >= ((threadID * (threadID + 1)) >> 1));
                                if (doit) {
                                    target[r]->d += threadID+1;
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                }else{
                                    //if (fail <= nbitems) std::this_thread::sleep_for(std::chrono::seconds(1)); // give a chance to other threads if scanning
                                    if (noisythread == threadID) foreach.fprintf(stdout," could not update, f{%i] = %i is val!\n", r, target[r]->d );
                                }
                            }
                        }else{
                            if (auto danew = target.tryRead(r)) {
                                doit = (danew->d >= ((threadID * (threadID + 1)) >> 1));
                                if ((!doit)&&(noisythread == threadID)) foreach.fprintf(stdout," could not update try, f{%i] = %i is val!\n", r, danew->d );
                                danew.unlock();
                                //if ((!doit)&&(fail <= nbitems)) std::this_thread::sleep_for(std::chrono::seconds(1));
                            }else doit = false;

                            if (doit){
                                if (auto daacc = target.tryWrite(r)){
                                    dabuffer.addEntry(r); foreach.updateProgress(threadID);
                                    daacc->d += threadID +1;

                                }
                            }// give a chance to other threads if scanning
                        }
                    }

                   // if ((fail < nbitems)&&(threadID == 0))  printf("%i\n", fail);
                }else {
                    //if ((rand() & 1023) == 0) std::this_thread::sleep_for(std::chrono::seconds(1));
                    fail++;
                    if (fail >= nbitems * 10) {
                        fail = 0;
                        foreach.fprintf(stdout," is now scanning to finish!!!\n", threadID);
                        /*if (auto ite = target.getIterator()) do{
                            if ((rand() & 255) == 0) printf("F[%i] = %i\n", ite(), ite->d);
                        }while(ite++);*/
                    }
                }
            }
            foreach.finishProgress(threadID);
            foreach.fprintf(stdout,"\n%i %X exit!\n", threadID, target.queryChunkIO());
           /* if (auto ite = target.getIterator()) do{
                foreach.printf("F[%i] = %i\n", ite(), ite->d);
            }while(ite++);
            foreach.printf("\n%i %X\n", threadID, target.queryChunkIO());*/
        return 0;}

        void iterMadness(){
            /*for(int i =0;i<10;i++){
                if (auto ite = target.getIterator()) do{
                    if ((rand() & 255) == 0) printf("F[%i] = %i\n", ite(), ite->d);
                }while(ite++);
                printf("sleepy\n");
                std::this_thread::sleep_for(std::chrono::seconds(2));
            }*/
            printf("im done\n");
        }
        std::function<void()> fncIterMadness(){return std::bind(&crazyTask::iterMadness,*this);}
    };


    LFHPrimitive::AsyncHashmap<KeyElem<uint32_t, uint32_t> > test;
    crazyTask lotask(test, nbitems);
    lotask.nbnastythreads = nbthreads;
    lotask.iteratortrap =0;
    lotask.noisythread = noisythread;
    lotask.mask = (lotask.nbnastythreads >> 1)-1; ExOp::toLeftFlood(lotask.mask);
    foreach.startThreadArray(nbthreads);
    printf("Doing this with %i threads (mask %X)\n",foreach.getThreadArraySize(), lotask.mask);

    foreach.startProgress("Lets do this!", lotask.nbnastythreads * lotask.nbitems);
    for(int i = lotask.nbnastythreads-1 ; i>0;i--) {
        foreach.submit(lotask, i);
        //std::this_thread::sleep_for(std::chrono::seconds(1));
    }
    foreach.submit_ThenWait(lotask, 0);
   // foreach.submit(lotask.fncIterMadness());
}

