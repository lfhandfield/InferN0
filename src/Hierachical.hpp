/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */


#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int nbdim,unsigned int nblead>

LFHTEMP	HyperPosition<C,nbdim,nblead>::HyperPosition(const Tuple<C, nbdim> &coor_min,const Tuple<C, nbdim> &coor_max){setFromMinMax(coor_min,coor_max);}
LFHTEMP	HyperPosition<C,nbdim,nblead>::HyperPosition(const Tuple<C, nbdim> &coor_cent, C &width){
    Tuple<C, nbdim> coor_min = coor_cent;
    Tuple<C, nbdim> coor_max = coor_cent;
    for(unsigned int i=nblead;i<nbdim;i++) coor_min[i] -= (width>>1)+ (width &1);
    for(unsigned int i=nblead;i<nbdim;i++) coor_max[i] += (width>>1);
    setFromMinMax(coor_min,coor_max);
}
LFHTEMP	int HyperPosition<C,nbdim,nblead>::pairBox_CenSiz(HyperPosition<C,nbdim,nblead> *fout, const Tuple<C, nbdim> &center,const Tuple<C, nbdim> &size, bool size_minus1){
	Tuple<C, nbdim> coor_min, coor_max;
	coor_min = center;
	for(unsigned int i =0 ; i< nbdim;i++){
		coor_min[i] -= (size[i]>> 1);
		coor_max[i] = coor_min[i] + size[i];
		if (!size_minus1) coor_max[i]--;
	}
return pairBox_MinMax(fout,coor_min, coor_max);}

LFHTEMP	int HyperPosition<C,nbdim,nblead>::pairBox_MinMax(HyperPosition<C,nbdim,nblead> *fout, const Tuple<C, nbdim> &coor_min, const Tuple<C, nbdim> &coor_max){
    C xorb,xorc;
    //C dabit;
    unsigned int best =0;
    unsigned int cur=0;
    best = nbdim -1;
    xorb = coor_min[best] ^ coor_max[best];
    for(cur= best-1; cur+1 > nblead;cur--){
        xorc = coor_min[cur] ^ coor_max[cur];
        if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
    }

    // found maximal box, rafine left and right
    // needs to find the top bit in xorb
    xorb |=  xorb >>1;
    xorb |=  xorb >>2;
    xorb |=  xorb >>4;
    if (sizeof(C)>1) xorb |=  xorb >>8;
    if (sizeof(C)>2) xorb |=  xorb >>16;
	C sepval = coor_max[best] & (0xFFFFFFFF ^ (xorb >> 1)) ;
    cur = nbdim -1;
    C xorbs = (cur == best) ?  coor_max[best] ^ sepval  : coor_min[cur] ^ coor_max[cur];
    unsigned int sbest = cur;
    for(cur--; cur+1 > nblead;cur--){
        xorc = (cur == best) ?  coor_max[best] ^ sepval  : coor_min[cur] ^ coor_max[cur];
        if ((xorbs <= xorc)&&((xorbs & xorc) <= (xorbs ^ xorc))) {sbest = cur; xorbs = xorc;}
    }
    xorbs |=  xorbs >>1;
    xorbs |=  xorbs >>2;
    xorbs |=  xorbs >>4;
    if (sizeof(C)>1) xorbs |=  xorbs >>8;
    if (sizeof(C)>2) xorbs |=  xorbs >>16;


    if (sizeof(C)==4) xorbs ^= 0xFFFFFFFF;
    else if (sizeof(C)==2) xorbs ^= 0xFFFF;

	for(cur = nblead; cur < sbest;cur++) fout[1][cur] = coor_max[cur]  & xorbs;

	if (xorbs) xorbs |= (xorbs >>1);
	else xorbs = (1 << ((sizeof(C) << 3)-1));
    for(; cur < nbdim;cur++) fout[1][cur] = coor_max[cur]  & xorbs;

	sepval--;

    cur = nbdim -1;
    xorbs = (cur == best) ?  coor_min[best] ^ sepval : coor_min[cur] ^ coor_max[cur];
    sbest = cur;
    for(cur--; cur+1 > nblead;cur--){
        xorc = (cur == best) ?  coor_min[best] ^ sepval : coor_min[cur] ^ coor_max[cur];
        if ((xorbs <= xorc)&&((xorbs & xorc) <= (xorbs ^ xorc))) {sbest = cur; xorbs = xorc;}
    }



    xorbs |=  xorbs >>1;
    xorbs |=  xorbs >>2;
    xorbs |=  xorbs >>4;
    if (sizeof(C)>1) xorbs |=  xorbs >>8;
    if (sizeof(C)>2) xorbs |=  xorbs >>16;

    if (sizeof(C)==4) xorbs ^= 0xFFFFFFFF;
    else if (sizeof(C)==2) xorbs ^= 0xFFFF;

    for(cur = nblead; cur < sbest;cur++) fout[0][cur] = ((cur == best)? sepval : coor_max[cur]) & xorbs;

	if (xorbs) xorbs |= (xorbs >>1);
	else xorbs = (1 << ((sizeof(C) << 3)-1));
      //  }
       // (*this)[cur] = coor_max[cur] & xorb;
    for(; cur < nbdim;cur++) fout[0][cur] = ((cur == best)? sepval : coor_max[cur]) & xorbs;





/*
//  CODE TO DETECT ERRORS
    Tuple<C, nbdim> q_min[2];
    Tuple<C, nbdim> q_max[2];


    q_min[0] = HyperCursor<C,nbdim,nblead>(fout[0]).getMin();
    q_min[1] = HyperCursor<C,nbdim,nblead>(fout[1]).getMin();
    q_max[0] = HyperCursor<C,nbdim,nblead>(fout[0]).getMax();
    q_max[1] = HyperCursor<C,nbdim,nblead>(fout[1]).getMax();

    for(cur =0;cur< 3;cur++){
		if (cur == best){
			if (q_max[0][cur]+1 != q_min[1][cur]) break;
			if (q_min[0][cur]> coor_min[cur]) break;
			if (q_max[1][cur]< coor_min[cur]) break;
		}else{
			if (q_min[0][cur]> coor_min[cur]) break;
			if (q_min[1][cur]> coor_min[cur]) break;
			if (q_max[0][cur]< coor_max[cur]) break;
			if (q_max[1][cur]< coor_max[cur]) break;
		}
    }

    printf("found error ? %c\n", (cur < 3) ? 'Y' : 'N');
 //   if (cur < 3){
		printf("[%X-%X,%X-%X,%X-%X] mainsep: d:%i sepval: %X\n", coor_min[0],coor_max[0],coor_min[1],coor_max[1],coor_min[2],coor_max[2] , best, sepval);
		printf("bounds:\n");
		printf("LowBox : [%X-%X,%X-%X,%X-%X]\n", q_min[0][0],q_max[0][0],q_min[0][1],q_max[0][1],q_min[0][2],q_max[0][2]);
		printf("HighBox: [%X-%X,%X-%X,%X-%X]\n", q_min[1][0],q_max[1][0],q_min[1][1],q_max[1][1],q_min[1][2],q_max[1][2]);


       printf("best %i, xorb %X\n", best, xorb);

        printf("base error min %X\t%X\t%X\n",q_min[0][0] ^ coor_min[0],q_min[0][1] ^ coor_min[1],q_min[0][2] ^ coor_min[2]);

        printf("base error max %X\t%X\t%X\n",q_max[0][0] ^ coor_max[0],q_max[0][1] ^ coor_max[1],q_max[0][2] ^ coor_max[2]);
        printf("bmin %X\t%X\t%X\n",q_min[0][0],q_min[0][1],q_min[0][2]);
        printf("bmax %X\t%X\t%X\n",q_max[0][0],q_max[0][1],q_max[0][2]);

        printf("min %X\t%X\t%X\n",coor_min[0],coor_min[1],coor_min[2]);
        printf("max %X\t%X\t%X\n",coor_max[0],coor_max[1],coor_max[2]);
        printf("fout %X\t%X\t%X\n",fout[0][0],fout[0][1],fout[0][2]);

        printf("input xor %X\t%X\t%X\n",coor_min[0] ^ coor_max[0],coor_min[1] ^ coor_max[1],coor_min[2] ^ coor_max[2]);
        printf("input andxor %X\t%X\t%X\n",coor_max[0] & (coor_min[0] ^ coor_max[0]),coor_max[1] & (coor_min[1] ^ coor_max[1]), coor_max[2] & (coor_min[2] ^ coor_max[2]));
   //     LFH_exit(1);
  //  }
*/

	if (fout[0].getBrother() != fout[1]) return 2;
	fout[0] = fout[1].getParent();
    return 1;
}

LFHTEMP	void HyperPosition<C,nbdim,nblead>::setFromMinMax(const Tuple<C, nbdim> &coor_min,const Tuple<C, nbdim> &coor_max){
		C xorb;
		C xorc;
		unsigned int best =0;
		unsigned int cur=0;
		for(cur=0;cur<nblead;cur++) {(*(Tuple<C,nbdim>*)this)[cur] = coor_min[cur]; if ((coor_min[cur]) !=  coor_max[cur]) {static_warning_handdle << LFH_WARNING_UNEXPECTED_INPUT;}}
		if (nbdim != nblead) {
/*
		xorb = (coor_min[cur]) ^ (coor_max[cur]);
		for(cur = 1; cur < nbdim;cur++){
			xorc = (coor_min[cur]) ^ (coor_max[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}
		xorb >>= 1;
		for(cur = best+1; cur < nbdim;cur++){
			xorc = (coor_min[cur]) ^ (coor_max[cur]);
			if ((xorb > xorc)||((xorb & xorc) > (xorb ^ xorc))) break;
			}
		best = cur-1;*/

        best = nbdim -1;
		xorb = coor_min[best] ^ coor_max[best];
		for(cur= best-1; cur+1 > nblead;cur--){
			xorc = coor_min[cur] ^ coor_max[cur];
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

		if (xorb){
			xorb |=  xorb >>1;
			xorb |=  xorb >>2;
			xorb |=  xorb >>4;
			if (sizeof(C)>1) xorb |=  xorb >>8;
			if (sizeof(C)>2) xorb |=  xorb >>16;
	//		if (sizeof(C)>4) xorb |=  xorb >>32;
			if (sizeof(C)==4) xorb ^= 0xFFFFFFFF;
            else if (sizeof(C)==2) xorb ^= 0xFFFF;

         //   if (cur < best){
			for(cur = nblead; cur < best;cur++) (*this)[cur] =  coor_max[cur] & xorb;
			xorb |= xorb >>1;
          //  }
           // (*this)[cur] = coor_max[cur] & xorb;
			for(; cur < nbdim;cur++) (*this)[cur] =  coor_max[cur] & xorb;
		//	printf("%x\n",xorb);

		}else{ // smallest possible!

			}

		}
}


LFHTEMP	void HyperPosition<C,nbdim,nblead>::setFromCenterWidth(const Tuple<C, nbdim> &coor_center, C width){
		C xorb;
		C xorc;
		unsigned int best =0;
		unsigned int cur=0;
		if (nbdim != nblead) {

        best = nbdim -1;
		xorb = (coor_center[best] -((width+1) >> 1)) ^ (coor_center[best] +(width >> 1));
		for(cur= best-1; cur+1 > nblead;cur--){
			xorc = (coor_center[cur] -((width+1) >> 1)) ^ (coor_center[cur] +(width >> 1));
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

		if (xorb){
			xorb |=  xorb >>1;
			xorb |=  xorb >>2;
			xorb |=  xorb >>4;
			if (sizeof(C)>1) xorb |=  xorb >>8;
			if (sizeof(C)>2) xorb |=  xorb >>16;
	//		if (sizeof(C)>4) xorb |=  xorb >>32;
			if (sizeof(C)==4) xorb ^= 0xFFFFFFFF;
            else if (sizeof(C)==2) xorb ^= 0xFFFF;

         //   if (cur < best){
			for(cur = nblead; cur < best;cur++) (*this)[cur] =  (coor_center[cur] +(width >> 1)) & xorb;
			xorb |= xorb >>1;
          //  }
           // (*this)[cur] = coor_max[cur] & xorb;
			for(; cur < nbdim;cur++) (*this)[cur] =  (coor_center[cur] +(width >> 1)) & xorb;
		//	printf("%x\n",xorb);

		}else{ // smallest possible!

			}

		}
}
	// init Brick!
LFHTEMP	HyperPosition<C,nbdim,nblead>::HyperPosition(const Tuple<C, nbdim> &coor_in,unsigned short magnitude, char last_dim){
	C mask;
	if (sizeof(C) == 4) mask = 0xFFFFFFFF;

	if (magnitude & 0x0010) mask <<= 16;
	if (magnitude & 0x0008) mask <<= 8;
	if (magnitude & 0x0004) mask <<= 4;
	if (magnitude & 0x0002) mask <<= 2;
	if (magnitude & 0x0001) mask <<= 1;

	int i;
	for(i=0;i< nbdim-last_dim;i++){
		(*(Tuple<C, nbdim>*)this)[i] = coor_in[i] & mask;
		}
		mask |= mask >>1;
	for(;i< nbdim;i++){
		(*(Tuple<C, nbdim>*)this)[i] = coor_in[i] & mask;
		}


	}
LFHTEMP char HyperPosition<C,nbdim,nblead>::comparesimple(const HyperPosition<C,nbdim,nblead>& other) const{
		C xorb;
		C xorc;
		unsigned int best;
		unsigned int cur=0;
		for(cur=0;cur<nblead;cur++) if (((*this)[cur]) !=  other[cur]) return((*this)[cur] < other[cur] ? 2 : 1);
		if (nbdim == nblead) return(0);

	//	xorb = ((*this)[cur]) ^ (other[cur]);
	//	for(cur = 1; cur < nbdim;cur++){
	//		xorc = ((*this)[cur]) ^ (other[cur]);
	//		if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
	//	}
		best = nbdim -1;
		xorb = ((*this)[best]) ^ (other[best]);

		for(cur = nbdim-2; cur+1!= nblead;cur--){
			xorc = ((*this)[cur]) ^ (other[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

		if (xorb) return((*this)[best] < other[best] ? 2 : 1);
		else return(0);

	}
LFHTEMP SETCMP_enum HyperPosition<C,nbdim,nblead>::compare(const HyperPosition<C,nbdim,nblead>& other) const{
		C xorb;
		C xorc;
		unsigned int best;
		unsigned int cur=0;
		for(cur=0;cur<nblead;cur++) if (((*this)[cur]) !=  other[cur]) return((*this)[cur] < other[cur] ? SETCMP_LT : SETCMP_GT);
		if (nbdim == nblead) return(SETCMP_EQUAL);


	//	xorb = ((*this)[cur]) ^ (other[cur]);
	//	for(cur = 1; cur < nbdim;cur++){
	//		xorc = ((*this)[cur]) ^ (other[cur]);
	//		if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
	//	}
		best = nbdim -1;
		xorb = ((*this)[best]) ^ (other[best]);

		for(cur = nbdim-2; cur+1 != nblead;cur--){
			xorc = ((*this)[cur]) ^ (other[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

		if (xorb) {
			return((*this)[best] < other[best] ? SETCMP_LT : SETCMP_GT);
		}else return(SETCMP_EQUAL);
	}
LFHTEMP template<class D> SETCMP_enum HyperPosition<C,nbdim,nblead>::compare(const HyperPosition<D,nbdim,nblead>& other) const{ // assumes sizeof(D) > sizeof(C)
	C xorb;
	C xorc;
	C mmmm; ExOp::toMax(mmmm);
	unsigned int best;
	unsigned int cur=0;
	for(cur=0;cur<nblead;cur++) if ( (mmmm & (other[cur] ^ ((*this)[cur]))) !=  0 ) return((*this)[cur] < other[cur] ? SETCMP_LT : SETCMP_GT);
	if (nbdim == nblead) return(SETCMP_EQUAL);


//	xorb = ((*this)[cur]) ^ (other[cur]);
//	for(cur = 1; cur < nbdim;cur++){
//		xorc = ((*this)[cur]) ^ (other[cur]);
//		if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
//	}
	best = nbdim -1;
	xorb = mmmm & ((other[best]) ^ ((*this)[best]));

	for(cur = nbdim-2; cur+1 != nblead;cur--){
		xorc = mmmm & ((other[cur]) ^ ((*this)[cur]));
		if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
	}

	if (xorb) {
		return(  (*this)[best] < (mmmm & other[best]) ? SETCMP_LT : SETCMP_GT);
	}else return(SETCMP_EQUAL);
}
LFHTEMP SETCMP_enum HyperPosition<C,nbdim,nblead>::compare(const Tuple<C, nbdim>& query) const{
    C xorb;
    C xorc;
    C mmmm; ExOp::toMax(mmmm);
    unsigned int best;
    unsigned int cur=0;
    for(cur=0;cur<nblead;cur++) if ( (mmmm & (query[cur] ^ ((*this)[cur]))) !=  0 ) return((*this)[cur] < query[cur] ? SETCMP_LT : SETCMP_GT);
    if (nbdim == nblead) return(SETCMP_EQUAL);

    best = nbdim -1;
    xorb = mmmm & ((query[best]) ^ ((*this)[best]));
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    if (best == lead) order++;
    bool inside = (xorb >> order) ==0;

    for(cur = nbdim-2; cur+1 != nblead;cur--){
        if (cur == lead) order++;
        xorc = mmmm & ((query[cur]) ^ ((*this)[cur]));
        inside &= (xorc >> order) ==0;
        if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
    }

    if (inside) return SETCMP_SUP_BOUNDED;
    else return(  (*this)[best] < (mmmm & query[best]) ? SETCMP_LT : SETCMP_GT);

}

LFHTEMP HyperPosition<C,nbdim,nblead>& HyperPosition<C,nbdim,nblead>::toLargestBoxPossible(){
    unsigned int i;
    for(i =0;i< nbdim-1;i++) ExOp::toZero((*this)[i]);
    (*this)[i] = 1 << ((sizeof(C)<< 3) -1);
return(*this);}
/*
LFHTEMP SETCMP_enum HyperPosition<C,nbdim,nblead>::compare(const Tuple<C, nbdim>& query) const{
    unsigned int order,lead,i;
    getOrderAndLeadHyper(order,lead);
    printf("%X and %X  %c%c\n", (query[0] ^ this->data[0])>>order,(query[1] ^ this->data[1])>>order,this->data[0] < query[0] ? 'L' : 'G',this->data[1] < query[1] ? 'L' : 'G');
    printf("q %X,%X\n",query[0],query[1]);
    printf("t %X,%X\n",this->data[0],this->data[1]);
    SETCMP_enum o = SETCMP_EQUAL;
    Tuple<unsigned short, 2> mima[2];
    mima[0] = this->getMin();
    mima[1] = this->getMax();
    printf("i %X,%X\n",mima[0][0],mima[0][1]);
    printf("a %X,%X\n",mima[1][0],mima[1][1]);
 //   for(i=0;i < nblead;i++) if (query[i] != this->data[i]) o = ((this->data[i] < query[i]) ? SETCMP_LT : SETCMP_GT;
    i=nbdim-1;
    for(;i >= lead;i--) if ((query[i] ^ this->data[i]) >> order) {o = (this->data[i] < query[i]) ? SETCMP_LT : SETCMP_GT; break;}
    if (o == SETCMP_EQUAL){
    order++;
    for(;i != nblead-1;i--) if ((query[i] ^ this->data[i]) >> order) {o = (this->data[i] < query[i]) ? SETCMP_LT : SETCMP_GT; break;}
    }
    if (o == SETCMP_EQUAL) {printf("E\n");return SETCMP_SUP_BOUNDED;}
    printf("%C\n", o == SETCMP_LT ? 'L' : 'G');
    return o;
}*/

/*
LFHTEMP void HyperPosition<C,nbdim,nblead>::getOrderAndLeadHyper(unsigned int &order, unsigned int &lead) const{
    unsigned int m = sizeof(nbdim) * 4;
    order=0;
    while(m >0){
	    if ( ((*this)[nblead] << (order+m)) != 0) order+= m;
	    m >>= 1;
	    }
	    order++;
	lead =nblead;
	if (order == 32) return;
	unsigned int i;
	for(i=nblead+1;i<nbdim;i++){
		if ( ((*this)[i] << order) != 0) {
			lead =i;order++;
			while((((*this)[i] << order) != 0)&&(order < 32)) order++;
			if (order == 32) return;
			}

		}
	}
*/
LFHTEMP void HyperPosition<C,nbdim,nblead>::wrMinMax(Tuple<C, nbdim>& f_min,Tuple<C, nbdim>& f_max) const{
    unsigned int i;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    for(i=0;i < nbdim;i++) f_min[i] = this->data[i];
    f_min[lead] -= (1 << order);
    for(i=0;i < nblead;i++) f_max[i] = this->data[i];
    for(;i < lead;i++) f_max[i] = (this->data[i] + (2 << order))-1;
    for(;i < nbdim;i++) f_max[i] = (this->data[i] + (1 << order))-1;
}
LFHTEMP Tuple<C, nbdim> HyperPosition<C,nbdim,nblead>::getMin() const{
    Tuple<C, nbdim> f_out;
    unsigned int i;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    for(i=0;i < nbdim;i++) f_out[i] = this->data[i];
    f_out[lead] -= (1 << order);
return(f_out);}
LFHTEMP Tuple<C, nbdim> HyperPosition<C,nbdim,nblead>::getMax() const{
    Tuple<C, nbdim> f_out;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    unsigned int i;
    for(i=0;i < nblead;i++) f_out[i] = this->data[i];
    for(;i < lead;i++) f_out[i] = (this->data[i] + (2 << order))-1;
    for(;i < nbdim;i++) f_out[i] = (this->data[i] + (1 << order))-1;
return(f_out);}
LFHTEMP Tuple<C, nbdim> HyperPosition<C,nbdim,nblead>::getCenter() const{
    Tuple<C, nbdim> f_out;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    unsigned int i;
    for(i=0;i < nblead;i++) f_out[i] = this->data[i];
    for(;i < lead;i++) f_out[i] = (this->data[i] + (1 << order));
    f_out[i] = this->data[i];
    if (order > 0)	for(i++ ;i < nbdim;i++) f_out[i] = (this->data[i] + (1 << (order-1)));
return(f_out);}
LFHTEMP C HyperPosition<C,nbdim,nblead>::getMin(unsigned int dir) const{
    C f_out;
    unsigned int i;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    f_out = this->data[dir];
    if (lead == dir) f_out -= (1 << order);
return(f_out);}
LFHTEMP C HyperPosition<C,nbdim,nblead>::getMax(unsigned int dir) const{
    C f_out;
    unsigned int order,lead;
    getOrderAndLeadHyper(order,lead);
    unsigned int i;
    f_out = this->data[dir];
    f_out += (((dir < lead) ? 2 : 1) << order);
    f_out--;
return(f_out);}
LFHTEMP unsigned int HyperPosition<C,nbdim,nblead>::getBestXorb(C &xorb) const{
    unsigned int m = sizeof(nbdim) * 4;
    unsigned int i;
    unsigned int o=0;
    while(m){
        for(i=nblead;i<nbdim;i++) if (this->data[i] & (0xFFFFFFFF >> (m+o))) break;
        if (i != nbdim) o+= m;
        m >>=1;
    }
    xorb = 0xFFFFFFFF >> o;
    if (i == nbdim) for(i=nblead;i<nbdim;i++) if (this->data[i] & (0xFFFFFFFF >> o)) break;
    if ((this->data[i] & (xorb >> 1))||(i == nbdim)) {printf("failure %i!\n",i); LFH_exit(1);}
return i;}

/*
LFHTEMP void HyperPosition<C,nbdim,nblead>::getOrderAndLeadHyper(unsigned int &order, unsigned int &lead) const{
    unsigned int m = sizeof(nbdim) * 4;
    order=32;
    while(m >0){
	    if ( ((*this)[nblead] << (32-order+m)) != 0) order-= m;
	    m >>= 1;
    }
    order--;
	lead =nblead;
	if (order == 0) return;
	unsigned int i;
	for(i=nblead+1;i<nbdim;i++){
		if ( ((*this)[i] << (32-order)) != 0) {
			lead =i;order--;
			while((((*this)[i] << (32-order)) != 0)&&(order > 0)) order--;
			if (order == 0) return;
        }
    }
}*/

LFHTEMP void HyperPosition<C,nbdim,nblead>::getOrderAndLeadHyper(unsigned int &order, unsigned int &lead) const{
    unsigned int i = nblead;
	if (nblead == nbdim) return;
	C msk = (*this)[i];
	for(i++;i<nbdim;i++) msk |= (*this)[i];
	order=0;
//	if ((TEMPLATE_SIZECHECK<C,4>::isGT)&&((msk & 0xFFFFFFFF) == 0) ) {msk >>= 32; order += 32;}
	if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((msk & 0x0000FFFF) == 0) ) {msk >>= 16; order += 16;}
	if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((msk & 0x00FF) == 0) ) {msk >>= 8; order += 8;}
	if ((msk & 0x0F) == 0) {msk >>= 4; order += 4;}
	if ((msk & 0x03) == 0) {msk >>= 2; order += 2;}
	if ((msk & 0x01) == 0) order++;
	for(i=nblead; ;i++) if ((*this)[i] & (1 << order)) break;
	lead = i;
}

LFHTEMP unsigned short HyperPosition<C,nbdim,nblead>::getMag() const{
	unsigned int order, lead;
	this->getOrderAndLeadHyper(order,lead);
	return (order<< 8) | lead;
}

LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getParent() const{
		unsigned int order,lead;
		getOrderAndLeadHyper(order,lead);
		HyperPosition<C,nbdim,nblead> f_out = (*this);
		if (lead == nbdim-1){
			f_out[lead] -= 1<< order;
			f_out[nblead] |= 1<< (1+order);
			}else{
			f_out[lead] -= 1<< order;
			f_out[lead+1] |= 1<< order;

		}
		return(f_out);
	}




LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getLeftChild() const{
		unsigned int order,lead;
		getOrderAndLeadHyper(order,lead);

		HyperPosition<C,nbdim,nblead> f_out = (*this);
		f_out[lead] -= 1<< order;
		if (lead == nblead){
			if (order > 0) f_out[nbdim-1] |= 1<< (order - 1);
			}else{
			f_out[lead-1] |= 1<< order;

		}
	return(f_out);
	}
LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getRightChild() const{
		unsigned int order,lead;
		this->getOrderAndLeadHyper(order,lead);
		HyperPosition<C,nbdim,nblead> f_out = (*this);
		if (lead == nblead){
			if (order > 0) f_out[nbdim-1] |= 1<< (order - 1);
			}else{
			f_out[lead-1] |= 1<< order;
		}
		return(f_out);
	}

LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getBrother() const{
		unsigned int order,lead;
		this->getOrderAndLeadHyper(order,lead);
		HyperPosition<C,nbdim,nblead> f_out = (*this);
		if (lead == nbdim-1){
			f_out[nblead] ^= 1<< (1+order);
			}else{
			f_out[lead+1] ^= 1<< order;
		}
		return(f_out);
	}


LFHTEMP void HyperPosition<C,nbdim,nblead>::show(FILE * f, int level) const{
    //HyperCursor<C,nbdim,nblead>(*this).show(f,level);
    unsigned int i=0;
    if (level == 0) for(i=0;i<nbdim;i++) fprintf(f,"%X%c", (*this)[i], i +1 == nbdim ? '\n' : '\t');
    else{for(i=0;i<nbdim-1;i++) fprintf(f,"%X;", (*this)[i]);
        fprintf(f,"%X", (*this)[i]);
    }
    //((Tuple<C,nbdim>*)this)->show(f, level);
}
/*LFHTEMP DebugBuffer<3u>& HyperPosition<C,nbdim,nblead>::dbg()const{ DebugBuffer<3u>& buffer = alloced_dbg_buffer; char* cur;
    //unsigned int i=0;
    //sprintf(buffer[0], "HyperPosition\n");
    //cur = buffer[1] + sprintf(buffer[1], "val=[";
    //for(i=0;i<nbdim;i++) cur += sprintf(cur,"%X%c", (*this)[i], i +1 == nbdim ? ']' : '\t');
    //cur = buffer[2] + sprintf(buffer[1], "rect=[";
    //for(i=0;i<nbdim;i++) cur += sprintf(cur,"%X-%X%c", (*this)[i], i +1 == nbdim ? ']' : ';');
return(buffer);}*/


    LFHTEMP double HyperPosition<C,nbdim,nblead>::rayIntersect(const Tuple<C,nbdim> &origin, const Tuple<C,nbdim> &target)const{  // returns the closest to origin
    unsigned int order;
    unsigned int lead;
    getOrderAndLeadHyper(order, lead);


    Tuple<C,nbdim> middle;
    Tuple<double,nbdim> tro, trd ;
    unsigned int i;
    for(i=0; i< nbdim;i++) {middle[i] = (*this)[i]; tro[i] = ((double)middle[i]) - origin[i]; trd[i] = ((double)target[i]) - origin[i];}
    printf("order %i  %i\n", order, lead);
    for(i= 0; i< lead;i++){tro[i] *= pow(0.5f, 33 -order);trd[i] *= pow(0.5f, 33 - order);}
    for(; i< nbdim;i++){tro[i] *= pow(0.5f, 32 -order);trd[i] *= pow(0.5f, 32 -order);}

    double polyno[3];

    polyno[0] = tro[0] * tro[0] - 4.0f;
    polyno[1] = trd[0] * tro[0]; // needs a *2
    polyno[2] = trd[0] * trd[0];

    for(i = 1; i< nbdim;i++){
    polyno[0] += tro[i] * tro[i];
    polyno[1] += trd[i] * tro[i];
    polyno[2] += trd[i] * trd[i];
    }
    // closest point at : num / den;

    polyno[0] = polyno[1]*polyno[1] - polyno[0] * polyno[2];

    if (polyno[0] < 0) return(-1);
    else{
        polyno[0] = sqrt(polyno[0]);
        if (fabs(polyno[1]) < fabs(polyno[0])) return(0.0f);
        else return ((polyno[1] + ((polyno[1] < 0.0f) ? -polyno[0] : polyno[0])) / polyno[2]);
    }
}

	LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getMinBox() const{
		unsigned int order,lead;    getOrderAndLeadHyper(order, lead);
		HyperPosition<C,nbdim,nblead>  f_out = *this;
		f_out[lead] ^= (1 << order);
		f_out[nblead] |= 1;
		return(f_out);
		}

	LFHTEMP HyperPosition<C,nbdim,nblead> HyperPosition<C,nbdim,nblead>::getMaxBox() const{
		HyperPosition<C,nbdim,nblead> f_out;
		unsigned int i;
		unsigned int order,lead;
		getOrderAndLeadHyper(order, lead);
		for(i=0;i < nblead;i++) f_out.data[i] = this->data[i];

		for(;i < lead;i++) f_out.data[i] = this->data[i] +(2 << order)-1;
		for(;i < nbdim;i++) f_out.data[i] =  this->data[i] + (1 << order)-1;
		//f_out.data[nblead] |= 1;
		return(f_out);
		}


#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int nbdim,unsigned int nblead>

LFHTEMP	HyperCursor<C,nbdim,nblead>::HyperCursor(const HyperPosition<C,nbdim,nblead>& in_hyper): hyperpos(in_hyper){
		unsigned int order,lead;
		in_hyper.getOrderAndLeadHyper(order,lead);
		//printf("convert: %i,%i\n", order, lead);
		mag = (order<< 8) | lead;
	}

LFHTEMP	HyperCursor<C,nbdim,nblead>::HyperCursor(const Tuple<C,nbdim>& in_coor, unsigned short magnitude){
        for(unsigned int i;i< nbdim;i++) hyperpos[i] = in_coor[i];
        setMagnitude(magnitude);
    }


LFHTEMP	HyperCursor<C,nbdim,nblead>::HyperCursor(const HyperCursor<C,nbdim,nblead>& in_hyper) :  hyperpos(in_hyper.hyperpos), mag(in_hyper.mag), par_mag(in_hyper.par_mag){}


LFHTEMP	SETCMP_enum HyperCursor<C,nbdim,nblead>::compareMagnitude(const HyperCursor<C,nbdim,nblead>& other) const{
	return (mag > other.mag) ? SETCMP_GT : (mag == other.mag) ? SETCMP_EQUAL: SETCMP_LT;
}

LFHTEMP	SETCMP_enum HyperCursor<C,nbdim,nblead>::compareMagnitude(const HyperPosition<C,nbdim,nblead>& other) const{
	unsigned int order,lead;
	other.getOrderAndLeadHyper(order,lead);
	unsigned short omag = (order<< 8) | lead; // safe for now...
	return (mag > omag) ? SETCMP_GT : (mag == omag) ? SETCMP_EQUAL: SETCMP_LT;
}
LFHTEMP	void HyperCursor<C,nbdim,nblead>::updateMagnitude(){
    unsigned int i = nblead;
	if (nblead == nbdim) return;
	C msk = hyperpos[i];
	for(i++;i<nbdim;i++) msk |= hyperpos[i];
	mag=0;
//	if ((TEMPLATE_SIZECHECK<C,4>::isGT)&&((msk & 0xFFFFFFFF) == 0) ) {msk >>= 32; order += 32;}
	if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((msk & 0x0000FFFF) == 0) ) {msk >>= 16; mag += 16;}
	if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((msk & 0x00FF) == 0) ) {msk >>= 8; mag += 8;}
	if ((msk & 0x0F) == 0) {msk >>= 4; mag += 4;}
	if ((msk & 0x03) == 0) {msk >>= 2; mag += 2;}
	if ((msk & 0x01) == 0) mag++;
	;
	for(i=nblead; ;i++) if (hyperpos[i] & (1 << mag)) break;
	mag = (mag << 8) |  i;
}
LFHTEMP	HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator=(const HyperCursor<C,nbdim,nblead> &other){
    hyperpos = other.hyperpos;
    mag = other.mag;
    par_mag= other.par_mag;
return *this;}
LFHTEMP	HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator=(HyperCursor<C,nbdim,nblead> &&other){
    hyperpos = other.hyperpos;
    mag = other.mag;
    par_mag= other.par_mag;
return *this;}

LFHTEMP	template<class D> HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator=(const HyperCursor<D,nbdim,nblead> &other){
    hyperpos = other.hyperpos;
    mag = other.mag;
    par_mag= other.par_mag;
    return (*this);
    }

LFHTEMP	HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator=(const HyperPosition<C,nbdim,nblead> &other){
    hyperpos = other; unsigned int order,lead;
    other.getOrderAndLeadHyper(order,lead);
	mag = (order<< 8) | lead;
    return (*this);
}

LFHTEMP HyperCursor<C,nbdim,nblead>::operator HyperPosition<C,nbdim,nblead>() const{return(hyperpos);}

	LFHTEMP SETCMP_enum HyperCursor<C,nbdim,nblead>::compare(const HyperCursor<C,nbdim,nblead>& other) const{
		return( hyperpos.compare(other.hyperpos));
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toMin(){
    unsigned int i;
    mag = 0;
    for(i=0;i<nbdim;i++) ExOp::toMin(hyperpos[i]);
    if (nblead != nbdim) ++hyperpos[nbdim-1];
    return(*this);
}
LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toMax(){
    unsigned int i;
    for(i=0;i<nbdim;i++)ExOp::toMax(hyperpos[i]);
    mag = 0;
    return(*this);}
LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getParent() const{
		int order = mag >> 8;
		int lead = mag & 255;
		HyperPosition<C,nbdim,nblead> f_out = hyperpos;
		f_out[lead] &= 0xFFFFFFFF ^ (1<< order);
		if (lead == nbdim-1){
			f_out[lead] -= 1<< order;
			f_out[nblead] |= 2<< order;
			}else{
			f_out[lead] -= 1<< order;
			f_out[lead+1] |= 1<< order;

		}
		return(f_out);
	}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkParent() const{ HyperCursor<C,nbdim,nblead> f_out;
    int order = mag >> 8;
    int lead = mag & 255;
    f_out.hyperpos = hyperpos;
    f_out.hyperpos[lead] &= 0xFFFFFFFF ^ (1<< order);
    if (lead == nbdim-1){
        f_out.hyperpos[lead] -= 1<< order;
        f_out.hyperpos[nblead] |= 2<< order;
        f_out.mag = ((mag +256) & 0xFF00) ;
        }else{
        f_out.hyperpos[lead] -= 1<< order;
        f_out.hyperpos[lead+1] |= 1<< order;
        f_out.mag = mag +1;
    }
    return(f_out);
}

LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkLeftChild() const{ HyperCursor<C,nbdim,nblead> f_out = hyperpos;
    int order = mag >> 8;
    int lead =  (mag & 255);
    f_out.hyperpos[lead] -= 1<< order;
    if (lead == nblead){
        f_out.hyperpos[nbdim-1] |= 1<< (order-1);
        f_out.mag = mag - (0x0101 + (nblead - nbdim));
    }else{
        f_out.hyperpos[lead-1] |= 1<< order;
        f_out.mag = mag -1;
    }
return(f_out);}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkRightChild() const{ HyperCursor<C,nbdim,nblead> f_out = hyperpos;
    int order = mag >> 8;
    int lead =  (mag & 255);
    if (lead == nblead){
        f_out.hyperpos[nbdim-1] |= 1<< (order-1);
        f_out.mag = mag - (0x0101 + (nblead - nbdim));
        }else{
        f_out.hyperpos[lead-1] |= 1<< order;
        f_out.mag = mag -1;
    }

return(f_out);}
LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getLeftChild() const{
    int order = mag >> 8;
    int lead =  (mag & 255);

    HyperPosition<C,nbdim,nblead> f_out = hyperpos;
    f_out[lead] -= 1<< order;
    if (lead == nblead){
        f_out[nbdim-1] |= 1<< (order-1);
    }else{
        f_out[lead-1] |= 1<< order;
    }
    return(f_out);
}
LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getRightChild() const{
    int order = mag >> 8;
    int lead =  (mag & 255);
    HyperPosition<C,nbdim,nblead> f_out = hyperpos;
    if (lead == nblead){
        f_out[nbdim-1] |= 1<< (order-1);
        }else{
        f_out[lead-1] |= 1<< order;
    }
    return(f_out);
}
LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getBrother() const{
    int order = mag >> 8;
    int lead =  (mag & 255);
    HyperPosition<C,nbdim,nblead> f_out = hyperpos;
    if (lead == nbdim-1){
        f_out[nblead] ^= 2<< order;
        }else{
        f_out[lead+1] ^=(1<< order);
    }
    return(f_out);
}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkBrother() const{
    int order = mag >> 8;
    int lead =  (mag & 255);
    HyperCursor<C,nbdim,nblead> f_out;
    f_out.hyperpos = hyperpos;
    f_out.mag = mag;
    if (lead == nbdim-1){
        f_out.hyperpos[nblead] ^= 2<< order;
        }else{
        f_out.hyperpos[lead+1] ^=(1<< order);
    }
    return(f_out);
}
LFHTEMP bool HyperCursor<C,nbdim,nblead>::isOlderBrother()const{
	int order = mag >> 8;
	int lead =  (mag & 255);
	if (lead == nbdim-1) return hyperpos[nblead] & (2<< order);
	else return hyperpos[lead+1] & (1<< order);
}

LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getCousin()const{
        HyperPosition<C,nbdim,nblead> f_out = hyperpos;
		if ((mag & 255) ==  nbdim - 1){
            if (nblead +1 == nbdim) f_out[ nblead ] ^= 4 << (mag >> 8);
			else f_out[ nblead + 1 ] ^= 2 << (mag >> 8);
		}else if ((mag & 255) ==  nbdim - 2){
            f_out[ nblead ] ^= 2 << (mag >> 8);
        }else f_out[ ((mag+2) & 255) ] ^= 1 << (mag >> 8);
        return(f_out);
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toParent(){
		hyperpos[ (mag & 255) ] ^= 1 << (mag >> 8);
		if ((mag & 255) ==  nbdim - 1){
			mag += 257 - nbdim + nblead ;
			}else mag++;
		hyperpos[ (mag & 255) ] |= 1 << (mag >> 8);
        return(*this);
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toBrother(){
		if ((mag & 255) ==  nbdim - 1){
			hyperpos[ nblead ] ^= 2 << (mag >> 8);
		}else hyperpos[ ((mag+1) & 255) ] ^= 1 << (mag >> 8);
        return(*this);
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toCousin(){
		if ((mag & 255) ==  nbdim - 1){
            if (nblead +1 == nbdim) hyperpos[ nblead ] ^= 4 << (mag >> 8);
			else hyperpos[ nblead + 1 ] ^= 2 << (mag >> 8);
		}else if ((mag & 255) ==  nbdim - 2){
            hyperpos[ nblead ] ^= 2 << (mag >> 8);
        }else hyperpos[ ((mag+2) & 255) ] ^= 1 << (mag >> 8);
        return(*this);
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toLeftChild(){
		hyperpos[ (mag & 255) ] ^= 1 << (mag >> 8);
		if (nblead == (int) (mag & 255)){
			mag -= 257 - nbdim + nblead ;
			}else mag--;
		hyperpos[ (mag & 255) ] |= 1 << (mag >> 8);
        return(*this);
	}
LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toRightChild(){
		if (nblead == (int) (mag & 255)){
			mag -= 257 - nbdim + nblead ;
			}else mag--;
		hyperpos[ (mag & 255) ] |= 1 << (mag >> 8);
        return(*this);
	}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toMagnitudeLeftChild(){
	if (mag < 256) return (*this);
	mag -= 256;
	hyperpos[ (mag & 255) ] ^= 3 << (mag >> 8);
	return(*this);
}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toMagnitudeParent(){
	int dir = nblead;
	for(;dir < (mag & 255); dir++) hyperpos[dir] = hyperpos[dir] ^ (hyperpos[dir] & (1 << (mag >> 8)));
	hyperpos[ (mag & 255) ] ^= 1 << (mag >> 8);
	mag += 256;
	hyperpos[ (mag & 255) ] |= 1 << (mag >> 8);
	for(dir++ ;dir < nbdim; dir++) hyperpos[dir] = hyperpos[dir] ^ (hyperpos[dir] & (1 << (mag >> 8)));
return(*this);}
LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toLower(int dir){
    int order = mag >> 8;
    if (dir > (mag &255)) hyperpos[dir] -= 1<< order;
    else  hyperpos[dir] -= 2<< order;
return(*this);}
LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toHigher(int dir){
    int order = mag >> 8;
    if (dir > (mag &255)) hyperpos[dir] += 1<< order;
    else  hyperpos[dir] += 2<< order;
return(*this);}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkLower(int dir)const{HyperCursor<C,nbdim,nblead> fout;
    int order;
    for(order=0;order< nbdim;order++) fout.hyperpos[order] = hyperpos[order];
    fout.mag = mag;
    order = mag >> 8;
    if (dir > (mag &255)) fout.hyperpos[dir] -= (1<< order);
    else  fout.hyperpos[dir] -= (2<< order);
return fout;}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::mkHigher(int dir)const{HyperCursor<C,nbdim,nblead> fout;
    int order;
    for(order=0;order< nbdim;order++) fout.hyperpos[order] = hyperpos[order];
    fout.mag = mag;
    order = mag >> 8;
    if (dir > (mag &255)) fout.hyperpos[dir] += (1<< order);
    else fout.hyperpos[dir] += (2<< order);
return fout;}
LFHTEMP void HyperCursor<C,nbdim,nblead>::wrReducedHyperCursor(HyperCursor<C,nbdim-1,nblead>& fout)const{
	unsigned int i;
	for(i =0;i<nbdim-1;i++) fout.hyperpos[i] = this->hyperpos[i];
	i--;
	fout.mag = this->mag;
	if ((this->mag & 0xFF) == nbdim-1){
		fout.mag--;
		fout.hyperpos[i] ^= (1 << (mag >> 8));
	}
}

LFHTEMP void HyperCursor<C,nbdim,nblead>::wrReducedHyperCursor(HyperCursor<C,nbdim-1,nblead>& fout, uint32_t dir)const{
	unsigned int i;
	fout.mag = this->mag;
    for(i =0;i<dir;i++) fout.hyperpos[i] = this->hyperpos[i];
    for(i++;i<nbdim;i++) fout.hyperpos[i-1] = this->hyperpos[i];
	if ((this->mag & 0xFF) == dir){
        if (dir == 0){
            fout.mag = (this->mag & 0xFF00) - 0x0102 + nbdim;
            fout.hyperpos[nbdim-2] ^= (1 << (mag >> 8));
        }else{
            fout.mag--;
            fout.hyperpos[dir-1] ^= (1 << (mag >> 8));
        }
	}else if ((this->mag & 0xFF) > dir) fout.mag--;
}


LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toLargestBoxPossible(){
    unsigned int i;
    for(i =0;i< nbdim-1;i++) ExOp::toZero(hyperpos[i]);
    hyperpos[i] = 1 << ((sizeof(C)<< 3) -1);
    mag = (nbdim - 1) | (((sizeof(C)<< 3) -1) << 8);
return(*this);}

LFHTEMP C HyperCursor<C,nbdim,nblead>::getMin(unsigned char dir) const{
return(hyperpos[dir] ^ ( dir == (mag & 255) ? 1 << (mag >> 8) : 0 ) );}
LFHTEMP C HyperCursor<C,nbdim,nblead>::getMax(unsigned char dir) const{
return( (hyperpos[dir] + ((dir < nblead) ? 0u : (((dir < (mag & 255) ) ? 2 : 1) << (mag >> 8) ) -1u )) );}


LFHTEMP Tuple<C, nbdim> HyperCursor<C,nbdim,nblead>::getMin() const{
    Tuple<C, nbdim> f_out;
    unsigned int i;
    for(i=0;i < nbdim;i++) f_out[i] = hyperpos[i];
    f_out[(mag & 255)] ^= (1 << (mag >> 8));
return(f_out);}
LFHTEMP Tuple<C, nbdim> HyperCursor<C,nbdim,nblead>::getMax() const{
    Tuple<C, nbdim> f_out;
    unsigned int order = (mag >> 8);
    unsigned int lead =  (mag & 255);
    unsigned int i;
    for(i=0;i < nblead;i++) f_out[i] = hyperpos[i];
    for(;i < lead;i++) f_out[i] = (hyperpos[i] + (2 << order))-1;
    for(;i < nbdim;i++) f_out[i] = (hyperpos[i] + (1 << order))-1;
return(f_out);}
	LFHTEMP Tuple<C, nbdim> HyperCursor<C,nbdim,nblead>::getCenter() const{
		Tuple<C, nbdim> f_out;
		unsigned int order = (mag >> 8);
		unsigned int lead =  (mag & 255);

		unsigned int i;
		for(i=0;i < nblead;i++) f_out[i] = hyperpos[i];
		for(;i < lead;i++) f_out[i] = (hyperpos[i] + (1 << order));
        f_out[i] = hyperpos[i];
        if (order > 0)	for(i++ ;i < nbdim;i++) f_out[i] = (hyperpos[i] + (1 << (order-1)));
		return(f_out);
		}

LFHTEMP void HyperCursor<C,nbdim,nblead>::wrMinMax(Tuple<C, nbdim>& f_min,Tuple<C, nbdim>& f_max) const{
    unsigned int order = (mag >> 8);
    unsigned int lead =  (mag & 255);
    unsigned int i;
    for(i=0;i < nblead;i++) f_min[i] = f_max[i] = hyperpos[i];
    for(;i < lead;i++) f_max[i] = (hyperpos[i] + (2 << order))-1;
    for(;i < nbdim;i++) f_max[i] = (hyperpos[i] + (1 << order))-1;
    for(i=nblead;i < nbdim;i++) f_min[i] = hyperpos[i];
    f_min[lead] -= (1 << order);
}
LFHTEMP void HyperCursor<C,nbdim,nblead>::wrMinMax(C& f_min,C& f_max, unsigned int dir) const{
    unsigned int order = (mag >> 8);
    unsigned int lead =  (mag & 255);
    f_min = hyperpos[dir];
    if (dir< nblead) f_max = hyperpos[dir];
    else if (dir < lead) f_max = (hyperpos[dir] + (2 << order))-1;
    else {
        f_max = (hyperpos[dir] + (1 << order))-1;
        if (dir == lead) f_min -= (1 << order);
    }
}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getMinBox() const{
	HyperCursor<C,nbdim,nblead>  f_out;
	f_out.hyperpos = hyperpos;
	f_out.hyperpos[(mag & 255)] ^= (1 << (mag >> 8));
	f_out.hyperpos[nblead] |= 1;
	f_out.mag =0;
	return(f_out);
}
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getMaxBox() const{
	HyperCursor<C,nbdim,nblead>  f_out;
	unsigned int i;
	unsigned int order = mag >> 8;
	unsigned int dir = mag & 255 ;
	for(i=0;i < nblead;i++) f_out.hyperpos[i] = hyperpos[i];
	for(;i < dir;i++) f_out.hyperpos[i] = hyperpos[i] +(2 << order)-1;
	for(;i < nbdim;i++) f_out.hyperpos[i] =  hyperpos[i] + (1 << order)-1;
	f_out.mag =0;
	return(f_out);
}

LFHTEMP bool HyperCursor<C,nbdim,nblead>::doesIntersect(const HyperPosition<C,nbdim,nblead>& other)const{
	Tuple<C,nbdim> mimas[4];
	mimas[0] = this->getMin();
	mimas[1] = this->getMax();
	mimas[2] = other.getMin();
	mimas[3] = other.getMax();
	for(unsigned int i=0;i<nbdim;i++){
		if (mimas[0][i] > mimas[3][i]) return false;
		if (mimas[1][i] < mimas[2][i]) return false;
	}
	return true;
}

LFHTEMP bool HyperCursor<C,nbdim,nblead>::isContainedBy_safe(const HyperPosition<C,nbdim,nblead>& other)const{
	Tuple<C,nbdim> mimas[4];
	mimas[0] = this->getMin();
	mimas[1] = this->getMax();
	mimas[2] = other.getMin();
	mimas[3] = other.getMax();
//	for(unsigned int i=0;i<4;i++){
//		printf("%X\t%X\t%X\n",mimas[i][0],mimas[i][1],mimas[i][2]);
//	}
	for(unsigned int i=0;i<nbdim;i++){
		if (mimas[0][i] < mimas[2][i]) return false;
		if (mimas[1][i] > mimas[3][i]) return false;
	}
	return true;
}

LFHTEMP bool HyperCursor<C,nbdim,nblead>::isContainedBy(const HyperPosition<C,nbdim,nblead>& other)const{
//	bool safe_ans = this->isContainedBy_safe(other);
	C msk;
	uint32_t i = 0;
	msk = other[0];
	for(i=1;i<nbdim;i++) msk |= other[i];
	int omag=0;
//	if ((TEMPLATE_SIZECHECK<C,4>::isGT)&&((msk & 0xFFFFFFFF) == 0) ) {msk >>= 32; omag += 32;}
	if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((msk & 0x0000FFFF) == 0) ) {msk >>= 16; omag += 16;}
	if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((msk & 0x00FF) == 0) ) {msk >>= 8; omag += 8;}
	if ((msk & 0x0F) == 0) {msk >>= 4; omag += 4;}
	if ((msk & 0x03) == 0) {msk >>= 2; omag += 2;}
	if ((msk & 0x01) == 0) omag++;
	if (omag < (mag >> 8)){
		//if (safe_ans != false) {printf("Wrong!"); LFH_exit(1);}
		return false;
	}

	for(i=0;;i++){
		if ((other[i]^ hyperpos[i]) >> (omag+1)) {
			//if (safe_ans != false) {printf("Wrong3!"); LFH_exit(1);}
			return false;
		}
		if (other[i] & (1 << omag)) break;
	}
	if ((omag == (mag >> 8))&&(i < (mag & 255))){
		//if (safe_ans != false) {	printf("Wrong2!"); LFH_exit(1);	}
		return false;
	}
	for(i++;i<nbdim;i++){
		if ((other[i]^ hyperpos[i]) >> omag) {
			//if (safe_ans != false) {printf("Wrong4!"); LFH_exit(1);}
			return false;
		}
	}
	//if (safe_ans != true) {	printf("Wrong5!"); LFH_exit(1);	}
	return true;
}


LFHTEMP unsigned short HyperCursor<C,nbdim,nblead>::commonContainer_mag(const HyperCursor<C,nbdim,nblead> &other) const{
		C xorb;
		C xorc;

		unsigned int cur=0;

		unsigned int best = nbdim -1;
		xorb = ((*this).hyperpos[best]) ^ (other.hyperpos[best]);

		for(cur = nbdim-2; cur != ((unsigned int)nblead) - 1u;cur--){
			xorc = ((*this).hyperpos[cur]) ^ (other.hyperpos[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}
        if (xorb == 0) return mag; // they are equal
    //    printf("dadif: %X\n", xorb);
        unsigned short damag =0;
		if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((xorb & 0xFFFF0000) != 0)) {xorb = xorb >>  16; damag += 16;}
		if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((xorb & 0xFF00) != 0)) {xorb = xorb >>  8; damag += 8;}
		if ((xorb & 0xF0) != 0) {xorb = xorb >>  4; damag += 4;}
		if ((xorb & 0x0C) != 0) {xorb = xorb >>  2; damag += 2;}
		if ((xorb & 0x02) != 0) {xorb = xorb >>  1; damag++;}
        if ((xorb & 0x01) == 0) damag--;
        damag = (damag << 8) + best;

        if (damag < mag) return(mag);
        else if (damag < other.mag) return(other.mag);
        else return(damag);
		}



LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toMinBoxExclusiveContainer(const Tuple<C,nbdim> &other){
	hyperpos[(mag & 255)] -= (1 << (mag >> 8));
	C xorb, xorc;
	unsigned int best = nbdim -1;
	unsigned int cur;
	xorb = ((*this).hyperpos[best] ^ other[best]);
	for(cur = nbdim-2; cur != ((unsigned int)nblead) - 1u;cur--){
		xorc = ((*this).hyperpos[cur]) ^ (other[cur]);
		if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
	}
	//printf("b%i cmp %X   %X + %X\n", best, other[best], hyperpos[best],  (((mag & 255) >= best ? 2 : 1) << (mag >> 8)) );
	if (other[best] >= hyperpos[best] + (((mag & 255) >= best ? 2 : 1) << (mag >> 8))){ // is way far off this one is ok
		hyperpos[(mag & 255)] += (1 << (mag >> 8));
		//printf("excape!\n");
		return *this;
	}

	if (xorb == 0) {hyperpos[nblead] ^= 1; mag = 0xFFFF; return *this;} // they are equal, no intersection, return illegal Hyperpos
	unsigned short damag =0;
	xorc = xorb;
	if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((xorb & 0xFFFF0000) != 0)) {xorb = xorb >>  16; damag += 16;}
	if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((xorb & 0xFF00) != 0)) {xorb = xorb >>  8; damag += 8;}
	if ((xorb & 0xF0) != 0) {xorb = xorb >>  4; damag += 4;}
	if ((xorb & 0x0C) != 0) {xorb = xorb >>  2; damag += 2;}
	if ((xorb & 0x02) != 0) damag++;

//	if ((1 << damag) == xorc){ // the only difference is the half in hyperpos, return illegal
//		mag = 0xFFFF; return *this;
//	}

	if (best == nblead) mag = ((damag -1) << 8) + nbdim -1;
	else mag = (damag << 8) + best -1;
//	printf("got max mag %X\n", mag);
	while(hyperpos[(mag & 255)] & (1 << (mag >> 8))){
		if ((mag & 255) != nblead) mag--;
		else if (mag == nblead) break;
		else mag -= 257 + nblead - nbdim;
	}
	// check if parent is "other", if so we are done!

//	printf("got final mag %X\n", mag);
	hyperpos[(mag & 255)] += (1 << (mag >> 8));
	return *this;
}

LFHTEMP HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::toLargeMinBoxContainer(){
	hyperpos[(mag & 255)] ^= (1 << (mag >> 8));
	C xorb;
	unsigned int best = nbdim-1;
	xorb = (*this).hyperpos[best];
	for(best--; best != ((unsigned int)nblead) - 1u;best--) xorb |=  ((*this).hyperpos[best]);
	unsigned short damag =0;

	if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((xorb & 0xFFFF) == 0)) {xorb = xorb >>  16; damag += 16;}
	if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((xorb & 0xFF) == 0)) {xorb = xorb >>  8; damag += 8;}
	if ((xorb & 0xF) == 0) {xorb = xorb >>  4; damag += 4;}
	if ((xorb & 0x3) == 0) {xorb = xorb >>  2; damag += 2;}
	if ((xorb & 0x1) == 0) damag++;

	if (xorb == 0){
		// xorb was 0
		mag = ((sizeof(C) << 11 ) -256) | (nbdim - 1);
	}else{
		for(best = nblead; best <= nbdim;best++){
			if ((*this).hyperpos[best] & (1 << damag)) break;
		}
		if (best == nblead) mag = ((damag -1) << 8) + nbdim -1;
		else mag = (damag << 8) + best -1;
	}
	hyperpos[(mag & 255)] ^= (1 << (mag >> 8));
	return *this;
}


LFHTEMP bool HyperCursor<C,nbdim,nblead>::findNext(){
	unsigned int lead = (mag & 255);
	unsigned int order = mag >> 8;
	hyperpos[lead] ^= (1<< order);
	do{
		if (lead == nbdim-1){
			lead = nblead;
			order++;
			if (order == (sizeof(C) <<3)) {mag = 0xFFFF; return false;}
		}else lead++;
		hyperpos[lead] ^= (1<< order);
	}while((hyperpos[lead] & (1<< order)) == 0);
	if (lead == nblead){
		lead = nbdim-1;
		order--;
	}else lead--;
	hyperpos[lead] ^= (1<< order);
	mag = (order <<8) | lead;
    return true;
}
/*
// Finds the largest box that contains cursor, that is bounded by a given rectangle
    LFHTEMP unsigned short HyperCursor<C,nbdim,nblead>::largest_Container_mag_rect(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max) const{
 		C xorb;
		C xorc;


      //  printf("%X\t%X\t%X\n", min[0],min[1],min[2]);
     //   printf("%X\t%X\t%X\n", max[0],max[1],max[2]);

        C damax; ExOp::toMax(damax);
        unsigned int best;
        for(best=0 ;best < nbdim;best++) if ((min[best] > 0)||(max[best] < damax)) break;
        if (best == nbdim) return (sizeof(C) << 11) + nbdim - 257;// silly rectangle: whole area!



        xorb = (min[best] > 0) ? ((*this).hyperpos[best]) ^ (min[best]-1) : damax;
    //    printf("%X b\n",xorb);
        xorc = (max[best] < damax) ? ((*this).hyperpos[best]) ^ (max[best]+1): damax;
   //     printf("%X b\n",xorc);
        if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) xorb = xorc;
        unsigned int cur = best;
        for(cur=1; cur < nbdim;cur++) {
            xorc = (min[cur] > 0) ?  ((*this).hyperpos[cur]) ^ (min[cur]-1) : damax;
     //               printf("%X c\n",xorc);
            if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
            xorc = (max[cur] < damax) ? ((*this).hyperpos[cur]) ^ (max[cur]+1): damax;
             printf("%X c %X\n",xorc, xorb);
            if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
        }


        printf("xorb %X", xorb);

        unsigned short damag =0;
		if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((xorb & 0xFFFF0000) != 0)) {xorb = xorb >>  16; damag += 16;}
		if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((xorb & 0xFF00) != 0)) {xorb = xorb >>  8; damag += 8;}
		if ((xorb & 0xF0) != 0) {xorb = xorb >> 4; damag += 4;}
		if ((xorb & 0x0C) != 0) {xorb = xorb >> 2; damag += 2;}
		if ((xorb & 0x02) != 0) {xorb = xorb >> 1; damag++;}
        if ((xorb & 0x01) == 0) damag--;
        printf("  mag cmp %i\n", damag);

        printf("%i nn\n", best);
        damag = (damag << 8) + best;
        HyperCursor<C,nbdim,nblead>::decr_mag(damag);
        return damag;
        }*/
    // Finds the largest box that contains cursor, that is bounded by a given rectangle
LFHTEMP unsigned short HyperCursor<C,nbdim,nblead>::largest_Container_mag_rect(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max) const{
    C xorb;
    C xorc;
    /*
    unsigned int best = nbdim -1;
    xorb = ((*this).hyperpos[best]) ^ (min[best]-1);
    xorc = ((*this).hyperpos[best]) ^ (max[best]+1);

    if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) xorb = xorc;
    unsigned int cur = best;
    while(cur-- > nblead){
        xorc = ((*this).hyperpos[cur]) ^ (min[cur]-1);
        if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
        xorc = ((*this).hyperpos[cur]) ^ (max[cur]+1);
        if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
    }
    */

  //  printf("%X\t%X\t%X\n", min[0],min[1],min[2]);
 //   printf("%X\t%X\t%X\n", max[0],max[1],max[2]);

    C damax; ExOp::toMax(damax);
    unsigned int best;
    for(best=0 ;best < nbdim;best++) if ((min[best] > 0)||(max[best] < damax)) break;
    if (best == nbdim) return (sizeof(C) << 11) + nbdim - 257;// silly rectangle: whole area!



    xorb = (min[best] > 0) ? ((*this).hyperpos[best]) ^ (min[best]-1) : damax;
//    printf("%X b\n",xorb);
    xorc = (max[best] < damax) ? ((*this).hyperpos[best]) ^ (max[best]+1): damax;
//     printf("%X b\n",xorc);
    if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) xorb = xorc;
    unsigned int cur;
    for(cur=1+best; cur < nbdim;cur++) {
        xorc = (min[cur] > 0) ?  ((*this).hyperpos[cur]) ^ (min[cur]-1) : damax;
 //               printf("%X c\n",xorc);
        if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
        xorc = (max[cur] < damax) ? ((*this).hyperpos[cur]) ^ (max[cur]+1): damax;
  //       printf("%X c %X\n",xorc, xorb);
        if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
    }


//    printf("xorb %X", xorb);

    unsigned short damag =0;
    if ((TEMPLATE_SIZECHECK<C,2>::isGT)&&((xorb & 0xFFFF0000) != 0)) {xorb = xorb >>  16; damag += 16;}
    if ((TEMPLATE_SIZECHECK<C,1>::isGT)&&((xorb & 0xFF00) != 0)) {xorb = xorb >>  8; damag += 8;}
    if ((xorb & 0xF0) != 0) {xorb = xorb >> 4; damag += 4;}
    if ((xorb & 0x0C) != 0) {xorb = xorb >> 2; damag += 2;}
    if ((xorb & 0x02) != 0) {xorb = xorb >> 1; damag++;}
    if ((xorb & 0x01) == 0) damag--;
//    printf("  mag cmp %i\n", damag);
   /*   best = 0;
    xorb = ((*this).hyperpos[best]) ^ (min[best]-1);
    xorc = ((*this).hyperpos[best]) ^ (max[best]+1);
    if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) xorb = xorc;
    cur = best;
    for(cur=1; cur < nbdim;cur++) {
        xorc = ((*this).hyperpos[cur]) ^ (min[cur]-1);
        if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
        xorc = ((*this).hyperpos[cur]) ^ (max[cur]+1);
        if ((xorb > xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {xorb = xorc; best = cur;}
    }

    damag =0;
    if ((sizeof(C) >= 4)&&((xorb & 0xFFFF0000) != 0)) {xorb = xorb >>  16; damag += 16;}
    if ((sizeof(C) >= 2)&&((xorb & 0xFF00) != 0)) {xorb = xorb >>  8; damag += 8;}
    if ((xorb & 0xF0) != 0) {xorb = xorb >>  4; damag += 4;}
    if ((xorb & 0x0C) != 0) {xorb = xorb >>  2; damag += 2;}
    if ((xorb & 0x02) != 0) {xorb = xorb >>  1; damag++;}
    if ((xorb & 0x01) == 0) damag--;
    printf("mag cur %i, %X\n", damag, xorb);*/
    damag = (damag << 8) + best;
   /* if (damag != 0){
        for(cur=0; cur < nbdim;cur++) {
            xorc = (min[cur] > 0) ?  ((*this).hyperpos[cur]) ^ (min[cur]-1) : damax;
            printf("min%i %X c!!\n",cur,xorc);
            xorc = (max[cur] < damax) ? ((*this).hyperpos[cur]) ^ (max[cur]+1): damax;
            printf("max%i %X c!!\n",xorc);
        }
    }
    printf("%i nn\n", best);*/

    if (damag != 0) HyperCursor<C,nbdim,nblead>::decr_mag(damag);
return damag;}

LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::NextTinyInBox(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max) const{
    HyperCursor<C,nbdim,nblead> fout = *this;
    Tuple<C, nbdim> toxor; ExOp::toZero(toxor);
    unsigned int lead = (mag & 255);
    unsigned int order = mag >> 8;
    fout.hyperpos[lead] ^= (1<< order); // remove size indicator
    if (lead == nbdim-1){
        lead = nblead;
        order++;
    }else lead++;

    while(order < 32){
        if (fout.hyperpos[lead] & (1<< order)) {
            if (((fout.hyperpos[lead] ^ min[lead]) >> order) != 0) fout.hyperpos[lead] ^= (1<< order);
            else toxor[lead] ^= (1<< order);
        }else{
            if (((max[lead] ^ fout.hyperpos[lead]) >> order) != 0){
                fout.hyperpos[lead] ^= (1<< order) ^ toxor[lead];
                break;
            }
        }
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
    }
    fout.mag = (order == 32) ? 0xFFFF : 0;
    for(lead=nblead;lead<nbdim;lead++) if (min[lead] > fout.hyperpos[lead]) fout.hyperpos[lead] = min[lead];
    fout.hyperpos[0] |= 1;
return fout;}
// does not assume that current box is inside
LFHTEMP HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::NextTinyInBox2(const Tuple<C,nbdim>& min, const Tuple<C,nbdim>& max) const{
    HyperCursor<C,nbdim,nblead> fout = *this; // maintains min-coor
    uint16_t flag =0;
    unsigned int lead = (mag & 255);
    unsigned int order = mag >> 8;
    fout.hyperpos[lead] ^= (1<< order); // maintains min-coor
    for(uint32_t i=0;i< nbdim;i++){
        if (fout.hyperpos[i] > max[i]) flag |= (1 << i);
        else if (fout.hyperpos[i] + (( (i <= lead ? 2u : 1u) << order) -1u) < min[i]) flag |= (256 << i);
    }
    //printf("iniflag: %X\n", flag);
    while(true){
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
        if (order == sizeof(C)*8) {
            fout.mag = 0xFFFF;
            return fout;
        }
        if (fout.hyperpos[lead] & (1<< order)) {
            fout.hyperpos[lead] ^= (1<< order);
            if ((flag & (1 << lead))&&(fout.hyperpos[lead] <= max[lead])) flag ^= (1 << lead);
        }else{
            if (flag & (256 << lead)){
                //printf("nm %X vs %X\n", fout.hyperpos[lead] + ((2 << order)-1) , min[lead] );
                if (fout.hyperpos[lead] + ((2 << order)-1) >= min[lead]) flag ^= (256 << lead);
            }
            if ((flag == 0)&&((fout.hyperpos[lead] | (1<< order))  <= max[lead])) break;

        }
    }
    //printf("final: l%i o%i\n",lead, order);
    fout.hyperpos[lead] ^= (1<< order);



    fout.mag = 0;
    for(lead=nblead;lead<nbdim;lead++) if (min[lead] > fout.hyperpos[lead]) fout.hyperpos[lead] = min[lead];
    fout.hyperpos[0] |= 1;
    return fout;
}
LFHTEMP bool HyperCursor<C,nbdim,nblead>::isIntersectingRect(const Tuple<C,nbdim>& r_min, const Tuple<C,nbdim>& r_max)const{
    unsigned int order = (mag >> 8);
    unsigned int lead =  (mag & 255);
    unsigned int i;
    for(i=0; i < lead;i++){
        if (r_min[i] > (this->hyperpos[i] + (2 << order))) return false;
        if (r_max[i] < this->hyperpos[i]) return false;
    }
    if (r_min[i] > (this->hyperpos[i] + (1 << order))) return false;
    if (r_max[i] < this->hyperpos[i] - (1 << order)) return false;
    for(i++; i < nbdim;i++){
        if (r_min[i] > (this->hyperpos[i] + (1 << order))) return false;
        if (r_max[i] < this->hyperpos[i]) return false;
    }
    return true;
}


	LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::commonContainer(const HyperCursor<C,nbdim,nblead> &other) const{HyperCursor<C,nbdim,nblead>  f_out = *this;
        unsigned int da_mag = this->commonContainer_mag(other);
        if (da_mag < mag) f_out.setMagnitude(da_mag);
	return(f_out);}

	LFHTEMP bool HyperCursor<C,nbdim,nblead>::isRealParentDirect() const{return(false);}

LFHTEMP HyperPosition<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getRealParent() const{
    HyperPosition<C,nbdim,nblead> fout;
    int i;
    for(i=0;i < nblead;i++) fout[i] = hyperpos[i];
    if (par_mag == 65535){
        for(;i < nbdim;i++) ExOp::toZero(fout[i]);
    }else{
    unsigned int mask = 0xFFFFFFFF << (1 +(par_mag >> 8));
    int dir = (par_mag >> 2) & 63;
    for(;i < dir;i++) fout[i] = hyperpos[i] & mask;
    mask |= (1 << (par_mag >> 8));
    for(;i < nbdim;i++) fout[i] =  hyperpos[i] & mask;
    fout[dir] |= (1 << (par_mag >> 8));
    }
return(fout);}
LFHTEMP void HyperCursor<C,nbdim,nblead>::show(FILE * f, int level) const{
    Tuple<C,nbdim> tmp[2];
	unsigned int order, lead;
	this->hyperpos.getOrderAndLeadHyper(order,lead);
	if (mag != ((order << 8) | lead)){
        printf("%i\t%i\t%i\n", hyperpos[0], hyperpos[1], hyperpos[2]);
        printf("coor: [");
        tmp[0] = this->getMin(); tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) printf("%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) printf("%.4X-%.4X%c\n", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
		printf("error detected! Magnitude should be %X, not %X\n", ((order << 8) | lead) , mag);
		LFH_exit(1);
	}

    switch(level){
    case 0:
        fprintf(f,"coor: [");
        tmp[0] = this->getMin(); tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) fprintf(f,"%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) fprintf(f,"%.4X-%.4X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        fprintf(f,"cur %.4X, par %.4X\n", mag& 0xFFFF, par_mag & 0xFFFF );
    break;
    default:
        fprintf(f,"coor: [");
        tmp[0] = this->getMin();tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) fprintf(f,"%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) fprintf(f,"%.4X-%.4X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        fprintf(f," (T%.4XP%.4X) ", mag& 0xFFFF, par_mag & 0xFFFF);
    }
}
/*LFHTEMP DebugBuffer<3u>  HyperCursor<C,nbdim,nblead>::dbg() const{DebugBuffer<3u> buffer; char* cur;
    try {
    cur = buffer[0] + sprintf(buffer[0], "HyperCursor is mag");// %X and pmag %X", mag, par_mag);
    cur[0] = ' ';
    unsigned int i=0;
    cur = buffer[1] + sprintf(buffer[1], "val=[");
    //for(i=0;i<nbdim;i++) cur += sprintf(cur,"%X%c", hyperpos[i], i +1 == nbdim ? ']' : '\t');
    cur[0] = ' ';
    cur = buffer[2] + sprintf(buffer[2], "rect=[");
    unsigned int order = (mag >> 8);
    unsigned int lead =  (mag & 255);
    for(i=0;i < nblead;i++) cur += sprintf(cur,"%X%c", hyperpos[i], i +1 == nbdim ? ']' : ';');
    for(;i < lead;i++) cur += sprintf(cur,"%X-%X%c", hyperpos[i], hyperpos[i]+ (2 << order)-1, i +1 == nbdim ? ']' : ';');
    cur += sprintf(cur,"%X-%X%c", hyperpos[i] - (1 << order), hyperpos[i] + (1 << order)-1, i +1 == nbdim ? ']' : ';');
    for(i++;i < nbdim;i++) cur += sprintf(cur,"%X-%X%c", hyperpos[i], hyperpos[i] + (1 << order)-1, i +1 == nbdim ? ']' : ';');
    }catch(const std::exception&){
    }
return(buffer);} */
LFHTEMP void HyperCursor<C,nbdim,nblead>::show(char * b, int level) const{
    Tuple<C,nbdim> tmp[2];
	unsigned int order, lead;
	this->hyperpos.getOrderAndLeadHyper(order,lead);
	if (mag != ((order << 8) | lead)){
        sprintf("%i\t%i\t%i\n", hyperpos[0], hyperpos[1], hyperpos[2]);
        printf("coor: [");
        tmp[0] = this->getMin(); tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) printf("%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) printf("%.4X-%.4X%c\n", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
		printf("error detected! Magnitude should be %X, not %X\n", ((order << 8) | lead) , mag);
		LFH_exit(1);
	}

    switch(level){
    case 0:
        sprintf(b,"coor: [");
        tmp[0] = this->getMin(); tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) sprintf(b,"%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) sprintf(b,"%.4X-%.4X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        sprintf(b,"cur %.4X, par %.4X\n", mag& 0xFFFF, par_mag & 0xFFFF );
    break;
    default:
        sprintf(b,"coor: [");
        tmp[0] = this->getMin();tmp[1] = this->getMax();
        if (TEMPLATE_SIZECHECK<C,4>::isEQ) for(unsigned int i=0;i<nbdim;i++) sprintf(b,"%.8X-%.8X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        else for(unsigned int i=0;i<nbdim;i++) sprintf(b,"%.4X-%.4X%c", tmp[0][i], tmp[1][i], i+1 == nbdim ? ']' : ';');
        sprintf(b," (T%.4XP%.4X) ", mag& 0xFFFF, par_mag & 0xFFFF);
    }
}

LFHTEMP ERRCODE HyperCursor<C,nbdim,nblead>::save(FILE *f) const {return (fwrite(this, sizeof(HyperCursor<C,nbdim,nblead>) ,1, f) != 1) ? 1 : 0;}
LFHTEMP ERRCODE HyperCursor<C,nbdim,nblead>::load(FILE *f){return (fread(this, sizeof(HyperCursor<C,nbdim,nblead>) ,1, f) != 1) ? 1 : 0;}
LFHTEMP  const HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator++(){ // prefix
    unsigned int lead =  (mag & 255);
    unsigned int order = mag >> 8;
    if (lead == nbdim-1){
        lead = nblead;
        order++;
    }else lead++;
    while(order < 32){
        hyperpos[lead] ^= (1<< order);
        if (hyperpos[lead] & (1<< order)) break;
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
    }
return(*this);}
LFHTEMP  const HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::operator--(){ // prefix
    unsigned int lead =  (mag & 255);
    unsigned int order = mag >> 8;
    if (lead == nbdim-1){
        lead = nblead;
        order++;
    }else lead++;
    while(order < 32){
        hyperpos[lead] ^= (1<< order);
        if (!(hyperpos[lead] & (1<< order))) break;
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
    }
return(*this);}
LFHTEMP  HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getNext() const{ HyperCursor<C,nbdim,nblead> fout = *this;
    unsigned int lead =  (mag & 255);
    unsigned int order = mag >> 8;
    if (lead == nbdim-1){
        lead = nblead;
        order++;
    }else lead++;
    while(order < 32){
        fout.hyperpos[lead] ^= (1<< order);
        if (fout.hyperpos[lead] & (1<< order)) break;
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
    }
return(fout);}
LFHTEMP  HyperCursor<C,nbdim,nblead> HyperCursor<C,nbdim,nblead>::getPrev() const { HyperCursor<C,nbdim,nblead> fout = *this;
    unsigned int lead =  (mag & 255);
    unsigned int order = mag >> 8;
    if (lead == nbdim-1){
        lead = nblead;
        order++;
    }else lead++;
    while(order < 32){
        fout.hyperpos[lead] ^= (1<< order);
        if (!(fout.hyperpos[lead] & (1<< order))) break;
        if (lead == nbdim-1){
            lead = nblead;
            order++;
        }else lead++;
    }
return(fout);}

	/*

LFHTEMP const SetComparison& HyperPositionCube<C,nbdim>::compare(const HyperPosition<C,nbdim,0> & other) const {
	int i;
	if (other < min){
		return(SetComparison(3+8));
	}else if (other > max){
		return(SetComparison(3+4));
	}else{
		if (min == max) return(SetComparison(0));
		for(i=0;i<nbdim;i++){
			if ((other[i]< min[i])||(other[i]> max[i])) break;
			}
		if (i == nbdim) return(SetComparison(1));
		else return(SetComparison(3));
	}



	}

LFHTEMP const SetComparison& HyperPositionCube<C,nbdim>::compareInterval(const HyperPosition<C,nbdim,0> & q_min,const HyperPosition<C,nbdim,0> & q_max) const {
	int i;
	if (q_max < min){
		return(SetComparison(3+8));
	}else if (q_min > max){
		return(SetComparison(3+4));
	}else{
		if (min == max) return(SetComparison(0));
//		for(i=0;i<nbdim;i++){
//			if ((other[i]< min[i])||(other[i]> max[i])) break;
//			}
//		if (i == nbdim) return(SetComparison(1));
//		else

		C xorb;
		C xorc;
		int best =0;
		int cur =1;
		xorb = (q_min[0]) ^ (q_max[0]);
		for(cur = 1; cur < nbdim;cur++){
			xorc = (q_min[cur]) ^ (q_max[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}
		for(cur = 1 << (nbdimof(C)+2);cur >0;cur >>= 1) xorb |= xorb >> cur;
		xorc = ~xorb;
		for(cur = nbdim-1; cur >= 0;cur--){
			if (cur-1 == best) {xorb >>= 1;xorc = ~xorb;}

			if (((q_min[cur] | xorb) < min[cur])||((q_min[cur] & xorc) > max[cur])) break;
		}
		if (cur >= 0) return(SetComparison(3)); // no intersection!
		else return(SetComparison(3)); // loosy!
	}

	}

*/
    LFHTEMP const HyperCursor<C,nbdim,nblead>& HyperCursor<C,nbdim,nblead>::setMagnitude(unsigned short nmag){

        if (((unsigned short)sizeof(C) * 8) < (nmag >> 8)) {printf("HyperCursor error! Tried to make magnitude %x! critical!\n",nmag); LFH_exit(1);}
        C d_max = ExCo<C>::mkMaximum();
        d_max <<=  (nmag >> 8);
        unsigned int i;
  //      for(i= nblead; i< (nmag & 255);i++) hyperpos[i] &= d_max;
  //      d_max <<=  1;
  //      for(; i< nbdim;i++) hyperpos[i] &= d_max;
  //       printf("%i\t%i\t%i\n",  hyperpos[0],  hyperpos[1],  hyperpos[2]);

        for(i= nbdim-1; i > (nmag & 255);i--) hyperpos[i] &= d_max;
        d_max <<=  1;
        for(; i != (((unsigned int)nblead)-1);i--) hyperpos[i] &= d_max;
        hyperpos[nmag & 255] |= 1 << (nmag >> 8);
     //   printf("%i\t%i\t%i\n",  hyperpos[0],  hyperpos[1],  hyperpos[2]);
        mag = nmag;

      //  printf("has mag %i\n", mag);

   //     unsigned int j;
   //    hyperpos.getOrderAndLeadHyper(i, j);

       // printf("%i,%i,   %i\n",i,j, (i<<8) + j);


        return(*this);
   }

LFHTEMP bool HyperCursor<C,nbdim,nblead>::strictly_contains(const HyperCursor<C,nbdim,nblead> &other) const{
    if (other.mag > mag) return(false);
}

LFHTEMP bool HyperCursor<C,nbdim,nblead>::contains(const HyperPosition<C,nbdim,nblead> &other) const{
    int i;
    for(i=0;i<= (mag & 255);i++) if (((other.data[i] ^ this->hyperpos[i]) >>(1+(mag >> 8))) != 0) return false;
    for(;i< nbdim;i++) if (((other.data[i] ^ this->hyperpos[i]) >> (mag >> 8)) != 0) return false;
return true;}


#undef LFHTEMP
#define LFHTEMP template<class I, unsigned int D, unsigned int L, class N>

// finds the find element intersecting a ray, the magnitude of direction

// insert, mag = where.mag,  oldmag = par_mag
// delete, mag =  par_mag,  oldmag = mag = where.mag

LFHTEMP Vector<HyperCursor<I,D,L> > HierarchicalTree<I,D,L,N>::makeBoxCover(const Tuple<I,D> &min, const Tuple<I,D> &max) const{
    Vector<HyperCursor<I,D,L> > fout;
    // trying to minimize K * volume
    Tuple<int32_t, D> magn;
    for(uint32_t i=0;i<D;i++) magn[i] = ExOp::mkMagnitude(min[i] ^ max[i]);
    //Tuple<uint32_t> ord = magn.mkOrdering();
return(fout);}

#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES>

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::updateParMag(const HyperCursor<STORETYPE,DIMS,LEAD> &where, unsigned short mag, unsigned short oldmag){
	//printf("is updating par mag %X -> %X at:\n", oldmag, mag);
    //where.show();fflush(stdout);
	HyperCursor<STORETYPE,DIMS,LEAD> cur = where.getMinBox();
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy,cur);
    HyperCursor<STORETYPE,DIMS,LEAD> max = where.getMaxBox();
    while((ite.isValid())&&(ite->k <= max)){
        while(ite->k.par_mag < oldmag){ // too small, uses the par_mag
            if (ite->k.mag == ite->k.par_mag) {printf("illegal par_mag! equal to self! %X\n", ite->k.par_mag);fflush(stdout);LFH_exit(1);}
            cur = ite->k;
            cur.setMagnitude(ite->k.par_mag);
            if (cur < ite->k){
                ite->k.show();cur.show();
                printf("moved backward... impossible! %X\n", ite->k.mag); fflush(stdout);
				ite.findGE(where);
				this->par_mag_uneq(where,ite,true);
				LFH_exit(1);
            }
            ite.findGE(cur);

            if ((!ite.isValid())||(ite->k.mag != cur.mag)){
                if (ite->k < cur) {

					printf("super impossible in par mag!\n");
					printf("found:\n");
					ite->k.show();
					printf("as greater than:\n");
					cur.show();
					printf("but it's the other way around!\n");
					LFH_exit(1);
				}
                printf("got:   ");ite->k.show();
                printf("wanted:");cur.show();
                printf("detected ill, expect mag %X got %X\n", cur.mag, ite->k.mag); fflush(stdout); LFH_exit(1);
                }
            if (ite->k > max){
                printf("parent outside! %X\n", ite->k.mag); fflush(stdout); LFH_exit(1);
                }
        }
        if (ite->k == where){ // ignore top

            //printf("skip equal\n"); // should not occur...
            do {++ite;}while((ite.isValid())&&(ite->k == where)); //
        }else{
            if (ite->k.par_mag > oldmag) {
                printf("IM mag %X  par_mag %X !!!\n", ite->k.mag , ite->k.par_mag ); ite->k.show();
                printf("ILLEGAL par mag\n"); LFH_exit(1);
            }
			//printf("Affected:"); ite->k.show(); fflush(stdout);
            hierarchy.orderpreserve_deref(ite).k.par_mag = mag;
            cur = ite->k.getTiny_Next();
            if (cur <= ite->k) break;
            ite.findGE(cur);
        }
    }
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::insertElem(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data){
    HyperCursor<STORETYPE,DIMS,LEAD> ney(key);

//	printf("finding first...\n"); fflush(stdout);
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy,key);
//	printf("done...\n"); fflush(stdout);

    if ((ite.isValid())&&(ite->k == key)) {
        ney.par_mag = ite->k.par_mag;
    }else{
        ney.par_mag = par_mag_uneq(key,ite);
        this->updateParMag(key,key.mag,ney.par_mag);
    }
	//printf("Inserted:");
	//ney.show();
//		printf("relay insertion...\n"); fflush(stdout);
    hierarchy.insert(KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >(ney,data));
//		printf("should be done...\n"); fflush(stdout);
}
/*
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::insertElem(const pair<Tuple<STORETYPE,DIMS>,Tuple<STORETYPE,DIMS> > &rect, const NODES &data){

    HyperPosition<STORETYPE, DIMS, LEAD> pairout[2];
    HyperPosition<STORETYPE, DIMS, LEAD>::pairBox_MinMax(pairout, rect.first, rect.second);
    {
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy,pairout[0]);
    HyperCursor<STORETYPE,DIMS,LEAD> ney(pairout[0]);
    if ((ite.isValid())&&(ite->k == pairout[0])) {
        ney.par_mag = ite->k.par_mag;
    }else{
        ney.par_mag = par_mag(pairout[0]);
        this->updateParMag(pairout[0],ney.mag,ney.par_mag);
    }
    hierarchy.insert(KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >(ney,data));
    }

    {
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2(hierarchy, pairout[1]);
    HyperCursor<STORETYPE,DIMS,LEAD> ney2(pairout[1]);
    if ((ite2.isValid())&&(ite2->k == pairout[1])) {
        ney2.par_mag = ite2->k.par_mag;
    }else{
        ney2.par_mag = par_mag(pairout[1]);
        this->updateParMag(pairout[1],ney2.mag,ney2.par_mag);
    }
    hierarchy.insert(KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >(ney2,data));
    }
}
*/
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::insertPartition(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data){
    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> ney = KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>(key, data);
  //  printf("at %X ", this); ExOp::show(this->max);
   data.trivialContract(ney.k);  // try to merge

 //   printf("ip start!\n"); fflush(stdout);
  unsigned short tmag = data.spiritBrother(ney.k);
  // printf("sb done!\n"); fflush(stdout);
   typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = this->find(ney);

   while(ite.isValid()){
        //printf("youppi! mag : dir %i, %X -> %X\n", 3-(data & 3), ney.k.mag,  tmag); fflush(stdout);
     //   this->removeElem(ney.k,data);
   typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2 = this->find_first(ney.k);

        if (ite2->d == data) {
            ++ite2;
            if ((ite2.isValid())&&(ite2->k == ney.k)) ite2.ordering_preserving_change()->k.par_mag = ite->k.par_mag;
            else this->updateParMag(ite->k,ite->k.par_mag,ite->k.mag);
        }
        this->remove(ite);
        ney.k.setMagnitude(tmag);
        tmag = data.spiritBrother(ney.k);

//
       ite = this->find(ney); ite2 = this->find_first(ney.k);
   }
 //   printf("p1 done!\n"); fflush(stdout);
   data.spiritBrother(ney.k);
    ite = this->find_first(ney.k);
//printf("sB done!\n"); fflush(stdout);
    if ((ite.isValid())&&(ite->k == ney.k)) {
        ney.k.par_mag = ite->k.par_mag;
    }else{
   //    printf("gpm todo!\n"); fflush(stdout);
        ney.k.par_mag = this->par_mag(ney.k);
   //     printf("tou todo!\n"); fflush(stdout);
        this->updateParMag(ney.k,ney.k.mag,ney.k.par_mag);
    }
//   printf("p2 done!\n"); fflush(stdout);
    if (ney.k.hyperpos[0] & 0xFFFF0000) {printf("ins ");ExOp::show(ney.k);}
    hierarchy.insert(KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >(ney.k,data));
 //   printf("finish!\n"); fflush(stdout);
}

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::insertCover(const Tuple<STORETYPE,DIMS> &min,const Tuple<STORETYPE,DIMS> &max, const NODES &data){
    Vector<HyperCursor<STORETYPE,DIMS,LEAD> > cover = this->makeBoxCover(min,max);
    for(uint32_t i=0;i< cover.getSize();i++) this->insertElem(cover[0],data);
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::removeCover(const Tuple<STORETYPE,DIMS> &min,const Tuple<STORETYPE,DIMS> &max, const NODES &data){
    Vector<HyperCursor<STORETYPE,DIMS,LEAD> > cover = this->makeBoxCover(min,max);
    for(uint32_t i=0;i< cover.getSize();i++) this->removeElem(cover[0],data);
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::removeElem(const HyperCursor<STORETYPE,DIMS,LEAD> &key, const NODES &data){
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy);
	ite.findGE(key);

    if (ite.isValid()){
		if (ite->k != key) {printf("Trying to delete non-existant!\n"); return;}
		if (ite->d != data){
			for(++ite;ite.isValid();++ite){
				if (ite->k != key) {printf("Trying to delete non-existant!\n"); return;}
				if (ite->d == data) break;
            }
			hierarchy.remove(ite); // par_mag can be ignored!
		}else this->remove(ite);
    } else printf("Could not find expected\n");
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::removePP(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite){
	if (&(ite.target) != &hierarchy) {printf("Trying to remove item within HierarchicalTree with erroneous Iterator"); LFH_exit(1);}
	HyperCursor<STORETYPE,DIMS,LEAD> thatkey = ite->k;
	this->remove(ite);
	ite.findGE(thatkey);
}

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::remove(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite){
	HyperCursor<STORETYPE,DIMS,LEAD> key = ite->k;
	//printf("Deleted:");
	//key.show();
	--ite;
	if (ite.isValid()){
		if (ite->k == key) {++ite; return hierarchy.remove(ite);}
		++ite;
	}else ite.findFirst();

	++ite;
	if (ite.isValid()){
		if (ite->k == key) ite.ordering_preserving_change()->k.par_mag = key.par_mag;
		else this->updateParMag(key,key.par_mag,key.mag);
		--ite;
	}else{
		this->updateParMag(key,key.par_mag,key.mag);
		ite.findLast();
	}
	hierarchy.remove(ite);
}


LFHTEMP NODES* HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &where){
typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy);
	ite.findGE(where);
    if ((ite.isValid())&&(ite->k == where)) return &(ite.ordering_preserving_change()->d);
    return NULL;
}
LFHTEMP const NODES* HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &where) const{
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = this->find_first(where);
    if ((ite.isValid())&&(ite->k == where)) return &(ite->d);
    return NULL;
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::query(const HyperCursor<STORETYPE,DIMS,LEAD> &key, Vector<NODES> &fout ){
    HyperCursor<STORETYPE,DIMS,LEAD> cur = key.getMinBox();
    HyperCursor<STORETYPE,DIMS,LEAD> max = key.getMaxBox();
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = this->find_first(cur);
    while((ite.isValid())&&(ite->k <= max)){
        fout.push_back(ite->d);
        ++ite;
    }

    unsigned short nmag = par_mag(key);
    while(nmag != 0xFFFF){
        max.setMagnitude(nmag);
        ite = this->find_first(max);
        if ((!ite.isValid())||(ite->k != max)) {printf("illegal parent pointer!\n"); max.show(); max = key; max.setMagnitude(nmag); max.show(); LFH_exit(1);}
        fout.push_back(ite->d);
        nmag = ite->k.par_mag;
        for(++ite;(ite.isValid())&&(ite->k == max);++ite) fout.push_back(ite->d);
        }
    fout.sort_unique();
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::query(Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite) const{
   // printf("offsets %i %i %i\n", hierarchy.min_offset,hierarchy.max_offset, 1 << hierarchy.alloc_mag );

    hierarchy.intersection(fout, box_ite); // find inner items
    HyperRectKeyIterator<STORETYPE, DIMS> box_iterator(box_ite.min, box_ite.max);
    HyperCursor<unsigned int, 3,0> cur,old;
/*    unsigned short nmag;
    // find outer items
    if (box_iterator.first()){
        old = box_iterator();
        nmag = this->par_mag(old);
        while(nmag != 0xFFFF){
            old.setMagnitude(nmag);
            break;
            }

    }*/

    fout.sort_unique();
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::query(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3> &pos, unsigned int width) const{

    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;

    box_ite.setCenterWidth(pos,width);

    this->query(box_ite, fout);

    fout.sort_unique();

}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::query(Vector<NODES> &fout, const Tuple<STORETYPE, 3> &pos, unsigned int width) const{

    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;

    box_ite.setCenterWidth(pos,width);
    Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > inner;
    this->query(inner, box_ite);

    fout.setSize(inner.getSize());
    for(unsigned int i=0;i<inner.getSize();i++) fout[i] = inner[i].d;
    fout.sort_unique();

}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryRect(Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3>&min, const Tuple<STORETYPE, 3>&max) const{
    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;

    box_ite.setMinandMax(min,max);

    this->query(fout, box_ite);
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3>&min, const Tuple<STORETYPE, 3>&max) const{

    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;

    box_ite.setMinandMax(min,max);
    Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > inner;
    this->query(inner, box_ite);

    for(unsigned int i=0;i<inner.getSize();i++) fout.push_back(inner[i].d);
    fout.sort_unique();

}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryCRect(Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3>&center, const Tuple<STORETYPE, 3>&size, bool size_minus1) const{
    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;
    box_ite.setCRect(center,size,size_minus1);
    this->query(fout, box_ite);
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryCRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3>&center, const Tuple<STORETYPE, 3>&size, bool size_minus1) const{

    HyperPositionQuery_Cube<STORETYPE,DIMS, LEAD, NODES> box_ite;

    box_ite.setCRect(center,size,size_minus1);
    Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > inner;
//	box_ite.show();
    this->query(inner, box_ite);

    for(unsigned int i=0;i<inner.getSize();i++) fout.push_back(inner[i].d);
    fout.sort_unique();

}

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryLRect(Vector<NODES> &fout, const Tuple<STORETYPE, 3> *mimas) const{
    HyperPositionQuery_LinkedCubes<STORETYPE, NODES> box_ite;
    box_ite.setQuery(mimas);
    this->query(fout, box_ite);
}

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::queryLRect(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &fout, const Tuple<STORETYPE, 3> *mimas) const{
    HyperPositionQuery_LinkedCubes<STORETYPE,NODES> box_ite;

    box_ite.setQuery(mimas);
    Vector<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > inner;

    this->query(inner, box_ite);

    for(unsigned int i=0;i<inner.getSize();i++) fout.push_back(inner[i].d);
    fout.sort_unique();

}


LFHTEMP unsigned short HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::par_mag_spacepartition(const HyperCursor<STORETYPE,DIMS, LEAD> &pos) const{
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy,pos);
    unsigned short tmag = pos.mag;
    unsigned short tmagp = pos.mag;

    unsigned int daind;
    HyperCursor<STORETYPE,DIMS,LEAD> t_pos;
    if (ite.isValid()){
     //   printf("path AA\n");
        tmag = pos.commonContainer_mag((*ite).k);
        if (pos == (*ite).k) { // got a perfect match!
           // printf("out!\n");fflush(stdout);
            //printf("Query Parent got %X! route A\n", (*ite).k.par_mag); fflush(stdout);
            return((*ite).k.par_mag);
        } else {
            daind = ite.Prev_index(); //printf("path AB %i %c\n", daind, (ite.hsPrev())? 'Y' :'N');
            if (daind != 0) {
                tmagp = pos.commonContainer_mag(hierarchy.deref_index(daind).k);
                t_pos = hierarchy.deref_index(daind).k; // printf("path AC\n");
                if (hierarchy.deref_index(daind).k.mag >= tmagp) tmagp = hierarchy.deref_index(daind).k.mag;
                else{
                while(hierarchy.deref_index(daind).k.par_mag < tmagp) {t_pos.setMagnitude(hierarchy.deref_index(daind).k.par_mag);daind = hierarchy.find_index(t_pos);}//printf("mag to par: %i->%i \n", this->deref_index(daind).k.mag, this->deref_index(daind).k.par_mag);}
                tmagp = hierarchy.deref_index(daind).k.par_mag;
                }

            } else tmagp = 0xFFFF;
        }    // first is outside! need to check both sides!
        t_pos = (*ite).k;
        if ((*ite).k.mag >= tmag) tmag = (*ite).k.mag;
        else{
        while((*ite).k.par_mag < tmag) {t_pos.setMagnitude((*ite).k.par_mag); ite.findGE(t_pos);}// printf("mag TO par: %i %i->%i \n", t_pos.k.mag, (*ite).k.mag, (*ite).k.par_mag);}
        tmag = (*ite).k.par_mag;
        tmag = (*ite).k.par_mag;
        }
        //printf("out! %x\t%x\n", tmagp, tmag); fflush(stdout);
        //printf("Query Parent got %X! route B \n", (tmagp < tmag) ? tmagp : tmag); fflush(stdout);
        return (tmagp < tmag) ? tmagp : tmag;            // if ((*ite).k)> pos) other = pos.prev
    }else if (hierarchy.getSize() == 0) return(0xFFFF);
    else{ // not valid, use the max
        //if (this->size == 0) {printf("Query Parent got %X! route C\n", 0xFFFF); fflush(stdout); return(0xFFFF);}
        daind = hierarchy.find_index(hierarchy.getMax().k);

        if (daind == 0){
            printf("supposed to be impossible (max of Tree corrupted!)\n"); fflush(stdout);
            printf("min:"); ExOp::show(hierarchy.getMin());
            printf("max:"); ExOp::show(hierarchy.getMax());

            //((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->show();
            LFH_exit(1);
        }
        tmag = pos.commonContainer_mag(hierarchy.getMax().k);

        if (hierarchy.deref_index(daind).k.mag >= tmag) return hierarchy.deref_index(daind).k.mag;
      //  printf("in rare! %i, %i\n", daind, tmag);fflush(stdout);
        t_pos = hierarchy.getMax().k;
       // this->show();unsigned int ccc =0;
        while(hierarchy.deref_index(daind).k.par_mag < tmag) {t_pos.setMagnitude(hierarchy.deref_index(daind).k.par_mag);daind = hierarchy.find_index(t_pos);}//  printf("daind: %i\n",daind);fflush(stdout); printf("%i->%i \n", (this->deref_index(daind).k.par_mag) >> 8, (this->deref_index(daind).k.par_mag) &255); fflush(stdout);if (ccc++ > 10) LFH_exit(1);}
     //   printf("out!\n");fflush(stdout);
        //printf("Query Parent got %X! route E\n", this->deref_index(daind).k.par_mag); fflush(stdout);
        return hierarchy.deref_index(daind).k.par_mag;
    }
}

LFHTEMP unsigned short HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::par_mag(const HyperCursor<STORETYPE,DIMS, LEAD> &pos) const{
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy,pos);
	return ((ite.isValid())&&(pos == (*ite).k)) ? ((*ite).k.par_mag) : this->par_mag_uneq(pos,ite);
}

LFHTEMP unsigned short HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::par_mag_uneq(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite, bool _debug) const{
 //   printf("g\n");fflush(stdout);
  //  printf("h\n");fflush(stdout);
    HyperCursor<STORETYPE,DIMS,LEAD> p_pos;
    HyperCursor<STORETYPE,DIMS,LEAD> n_pos;

	if (ite.isValid()) {n_pos = ite->k; --ite;
	}else if (hierarchy.getSize() == 0) return(0xFFFF);
	else{
		n_pos.mag = 0xFFFF;
		ite.findLast();
	}
	if (_debug) {printf("next found:"); n_pos.show();}



	if (!ite.isValid()) p_pos.mag = 0xFFFF;
	else{
		p_pos = ite->k;
		--ite;
		while((ite.isValid())&&(ite->k == p_pos)){
			p_pos.par_mag = ite->k.par_mag;
			--ite;
		}
		if (_debug) {printf("prev found:"); p_pos.show();}

		while(pos.commonContainer_mag(p_pos) != p_pos.mag){
			if (_debug) printf("set mag to %X\n", p_pos.par_mag);
			if (p_pos.par_mag == 0xFFFF) {p_pos.mag = 0xFFFF; break;}
			p_pos.setMagnitude(p_pos.par_mag);
			ite.findGE(p_pos);
			if (!ite.isValid()){ printf("did not find parent!\n");  fflush(stdout); LFH_exit(1);}
			p_pos = ite->k;
		}
	}
	if (_debug) printf("next now\n");
	if (n_pos.mag == 0xFFFF) return p_pos.mag;
	while(pos.commonContainer_mag(n_pos) != n_pos.mag){
		if (_debug) printf("set mag to %X\n", n_pos.par_mag);
		if (n_pos.par_mag == 0xFFFF)  {n_pos.mag = 0xFFFF; break;}
		n_pos.setMagnitude(n_pos.par_mag);
		ite.findGE(n_pos);
		if (!ite.isValid()){ printf("did not find parent!\n"); fflush(stdout); LFH_exit(1);}
		n_pos = ite->k;
	}
	if (_debug) {printf("looped to %X %X\n", p_pos.mag, n_pos.mag);}

	return (p_pos.mag < n_pos.mag) ? p_pos.mag : n_pos.mag;
}

/*

    if (ite.isValid()){
//printf("i\n");fflush(stdout);
     //      printf("f\n");fflush(stdout);
        if (pos == (*ite).k) return((*ite).k.par_mag); // got a perfect match!
            //printf("Query Parent got %X! route A\n", (*ite).k.par_mag); fflush(stdout);
		tmag = pos.commonContainer_mag((*ite).k);

   //     printf("iza here!\n");fflush(stdout);
		daind = ite.Prev_index(); //printf("path AB %i %c\n", daind, (ite.hsPrev())? 'Y' :'N');
		if (daind != 0) {
			tmagp = pos.commonContainer_mag(hierarchy.deref_index(daind).k);
			t_pos = hierarchy.deref_index(daind).k; // printf("path AC\n");
			if (hierarchy.deref_index(daind).k.mag >= tmagp) tmagp = hierarchy.deref_index(daind).k.mag;
			else{
			ite.findGE(hierarchy.deref_index(daind).k);
			while( ite->k.par_mag < tmagp) {t_pos.setMagnitude(ite->k.par_mag);tmagp = ite->k.par_mag; ite.findGE(t_pos); }//printf("mag to par: %i->%i \n", this->deref_index(daind).k.mag, this->deref_index(daind).k.par_mag);}
			tmagp = ite->k.par_mag;
			ite.findGE(pos);
			}
		} else tmagp = 0xFFFF;
            // first is outside! need to check both sides!
    //    printf("j\n");fflush(stdout);
        t_pos = (*ite).k;
        if ((*ite).k.mag >= tmag) tmag = (*ite).k.mag;
        else{
        while((*ite).k.par_mag < tmag) {t_pos.setMagnitude((*ite).k.par_mag); tmag = (*ite).k.par_mag; ite.findGE(t_pos);  }// printf("mag TO par: %i %i->%i \n", t_pos.k.mag, (*ite).k.mag, (*ite).k.par_mag);}
        tmag = (*ite).k.par_mag;
        }
    //    printf("out! %x\t%x\n", tmagp, tmag); fflush(stdout);
      //  printf("Query Parent got %X! route B \n", (tmagp < tmag) ? tmagp : tmag); fflush(stdout);
        return (tmagp < tmag) ? tmagp : tmag;            // if ((*ite).k)> pos) other = pos.prev

     }else if (hierarchy.getSize() == 0) return(0xFFFF);
     else{ // not valid, use the max
     //   printf("ee\n");
     //   printf("p-in rare! %i, %i\n", daind, tmag);fflush(stdout);
        ite.findGE(hierarchy.getMax().k); // guarrantied to work, finds the maximum
        if (!ite.isValid()){
            printf("supposed to be impossible (max of Tree corrupted!)\n"); fflush(stdout);
            printf("min:"); ExOp::show(hierarchy.getMin().k);
            printf("max:"); ExOp::show(hierarchy.getMax().k);

            //((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->show();
            LFH_exit(1);
        }
        tmag = pos.commonContainer_mag(hierarchy.getMax().k);
     //   printf("rr\n");
        if (ite->k.mag >= tmag) return ite->k.mag;
     //   printf("in rare! %i\n", tmag);fflush(stdout);
        t_pos = hierarchy.getMax().k;
       // this->show();unsigned int ccc =0;
        while(ite->k.par_mag < tmag) {t_pos.setMagnitude(ite->k.par_mag);ite.findGE(t_pos);}//  printf("daind: %i\n",daind);fflush(stdout); printf("%i->%i \n", (this->deref_index(daind).k.par_mag) >> 8, (this->deref_index(daind).k.par_mag) &255); fflush(stdout);if (ccc++ > 10) LFH_exit(1);}
     //   printf("Query Parent got %X! route E\n", ite->k.par_mag); fflush(stdout);
        return ite->k.par_mag;
    }*/

LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::rayIntersection(Vector<NODES> &f_out,  const Tuple<STORETYPE, DIMS> &origin, const Tuple<STORETYPE, DIMS> &target, double maxrange) const{
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		HyperCursor<STORETYPE,DIMS,LEAD> tmin = this->getMin().k;
		HyperCursor<STORETYPE,DIMS,LEAD> tmax = this->getMax().k;
        double intermed;

		if (this->size ==0) return;

		while( true ){

            intermed =  tmin.commonContainer(tmax).rayIntersect(origin,target);
			while(  (intermed >= 0)&&(intermed <= maxrange) ){
				if (this->tree[path[depth]].second - 2 > this->firstlone-2) {
                    intermed =  this->tree[path[depth]].first.k.hyperpos.rayIntersect(origin,target);
					if  ((intermed >= 0)&&(intermed <= maxrange)) f_out.push_back(this->tree[path[depth]].first.d);

					if (this->tree[path[depth]].second > this->firstlone){
                        intermed =  this->tree[this->tree[path[depth]].second].first.k.hyperpos.rayIntersect(origin,target);
						if ((intermed >= 0)&&(intermed <= maxrange))   f_out.push_back(this->tree[this->tree[path[depth]].second].first.d);
					}

					break;
				}
				tmax = this->tree[path[depth]].first.k;
				path[depth+1] = this->tree[path[depth]].second & 0xFFFFFFFE;
				depth++;

				intermed =  tmin.commonContainer(tmax).rayIntersect(origin,target);
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
            intermed =  this->tree[this->tree[path[depth]].second].first.k.hyperpos.rayIntersect(origin,target);
			if ((intermed >= 0)&&(intermed <= maxrange))  f_out.push_back(this->tree[path[depth]].first.d);

			tmin = this->tree[path[depth]].first.k;
			if (depth == 0) tmax = this->getMax().k;
			else tmax = this->tree[path[depth-1]].first.k;

			// erase current parent in history, no longer needed
			path[depth] = this->tree[path[depth]].second | 1;




		}
}
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::rayIntersection(Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >  > &f_out,  const Tuple<STORETYPE, DIMS> &origin, const Tuple<STORETYPE, DIMS> &target, double maxrange) const{

		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		HyperCursor<STORETYPE,DIMS,LEAD> tmin = this->min.k;
		HyperCursor<STORETYPE,DIMS,LEAD> tmax = this->max.k;
        double intermed;

		if (this->size ==0) return;

		while( true ){

            intermed =  tmin.commonContainer(tmax).rayIntersect(origin,target);
			while(  (intermed >= 0)&&(intermed <= maxrange) ){
				if (this->tree[path[depth]].second - 2 > this->firstlone-2) {
                    intermed =  this->tree[path[depth]].first.k.hyperpos.rayIntersect(origin,target);
					if  ((intermed >= 0)&&(intermed <= maxrange)) f_out.push_back(this->tree[path[depth]].first);

					if (this->tree[path[depth]].second > this->firstlone){
                        intermed =  this->tree[this->tree[path[depth]].second].first.k.hyperpos.rayIntersect(origin,target);
						if ((intermed >= 0)&&(intermed <= maxrange))   f_out.push_back(this->tree[this->tree[path[depth]].second].first);
					}

					break;
				}
				tmax = this->tree[path[depth]].first.k;
				path[depth+1] = this->tree[path[depth]].second & 0xFFFFFFFE;
				depth++;

				intermed =  tmin.commonContainer(tmax).rayIntersect(origin,target);
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
            intermed =  this->tree[this->tree[path[depth]].second].first.k.hyperpos.rayIntersect(origin,target);
			if ((intermed >= 0)&&(intermed <= maxrange))  f_out.push_back(this->tree[path[depth]].first);

			tmin = this->tree[path[depth]].first.k;
			if (depth == 0) tmax = this->max.k;
			else tmax = this->tree[path[depth-1]].first.k;

			// erase current parent in history, no longer needed
			path[depth] = this->tree[path[depth]].second | 1;




		}

}
LFHTEMP bool HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::checkHierarchySimple()const{
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(hierarchy);
	ite.findFirst();

    HyperCursor<STORETYPE,DIMS,LEAD> tmp;
    unsigned int qq;
    for(;ite.isValid();++ite){
        tmp = (*ite).k;
        if ((*ite).k.mag >= (*ite).k.par_mag) return false;
        if (tmp.par_mag != 0xFFFF){
			tmp.setMagnitude(tmp.par_mag);
			qq = hierarchy.find_index(tmp);
			if (qq == 0) return false;
		}
		continue;
		tmp = (*ite).k;
		tmp.toParent();
		while(tmp.mag != tmp.par_mag){
			qq = hierarchy.find_index(tmp);
			if (qq != 0){
				printf("Error, par_mag %X, but %X exists!", tmp.par_mag, tmp.mag);
				printf("Node");ite->k.show();
				printf("Intercept");tmp.show();
				return false;
			}
			tmp.toParent();
			if (tmp.mag > 0x0C00) break;
		}
    }
    return true;
    }// class HyperPositionQuery_Cube






LFHTEMP bool HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::checkHierarchyAt(typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator &ite) const{
	HyperCursor<STORETYPE,DIMS,LEAD> tmp = ite->k;
	--ite;
	if (ite.isValid()){
		if (ite->k == tmp) {++ite;return true;}
		++ite;
	}else ite.findFirst();
    unsigned int qq;
	if ((*ite).k.mag >= (*ite).k.par_mag) return false;
	if (tmp.par_mag != 0xFFFF){
		tmp.setMagnitude(tmp.par_mag);
		qq = hierarchy.find_index(tmp);
		if (qq == 0) return false;
	}
	tmp = (*ite).k;
	tmp.toParent();
	while(tmp.mag != tmp.par_mag){
		qq = hierarchy.find_index(tmp);
		if (qq != 0){
			printf("Error, par_mag %X, but %X exists!", tmp.par_mag, tmp.mag);
			printf("Node");ite->k.show();
			printf("Intercept");tmp.show();
			return false;
		}
		tmp.toParent();
		if (tmp.mag > 0x0C00) break;
	}
	return true;
}





/*
LFHTEMP void HierarchicalTree<STORETYPE, DIMS, LEAD, NODES>::expandAll(){
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = this->first();

    for(;ite.isValid();++ite){
        ite->d.trivialExpand(ite.ordering_preserving_change()->k);
    }
}*/

#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int  DIMS, unsigned int  LEAD, class NODES>
	LFHTEMP HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::HyperPositionQuery_Cube(const Tuple<STORETYPE,DIMS> i_min, const Tuple<STORETYPE,DIMS> i_max): min(i_min),max(i_max){

	}
	LFHTEMP HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::HyperPositionQuery_Cube(const HyperPosition<STORETYPE,DIMS,LEAD> dabox){
		HyperCursor<STORETYPE,DIMS,LEAD> dibox(dabox);
		min = dibox.getMin();
		max = dibox.getMax();
	}

	LFHTEMP void HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::setCRect(const Tuple<STORETYPE, DIMS>  &_center,const Tuple<STORETYPE, DIMS>  &_size, bool size_minus1){
		min = _center;
		for(unsigned int i =0; i<DIMS;i++){
			min[i] -= (_size[i] >> 1);
			max[i] = min[i] + _size[i];
			if (!size_minus1) max[i]--;
		}
	}

	LFHTEMP SETCMP_enum HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &q_min, const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &q_max) const{
		// check if the intervals are: disjoint, intersect, or contained!
   //     q_min.show();
   //     q_max.show();
		/*
		printf("Range Query:\n");
		q_min.show();
		q_max.show();*/


       // return SETCMP_EQUAL;

		STORETYPE xorb;
		STORETYPE xorc;
		STORETYPE xorq;
		int best;
		int cur=0;
		for(cur=0;cur<LEAD;cur++) if ((q_min.k.hyperpos[cur]) !=  q_max.k.hyperpos[cur]) break;
		if (cur < LEAD)return SETCMP_DISJOINT; // some day... fix this// disjoint... for now
		//	if (DIMS == LEAD) (...) min == max

		int qbest;
		if (q_min.k.mag > q_max.k.mag) {xorq = 1 << (q_min.k.mag>>8); qbest = (q_min.k.mag & 0xFF);}
        else {xorq = 1 << (q_max.k.mag>>8); qbest = (q_min.k.mag & 0xFF);}

		cur = DIMS-1;
		best = DIMS-1;
		xorb = (q_min.k.hyperpos[cur]) ^ (q_max.k.hyperpos[cur]);
		if (cur == qbest){
			if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) xorb = xorq;
		}
		for(cur--; cur >= (int)LEAD;cur--){
			xorc = (q_min.k.hyperpos[cur]) ^ (q_max.k.hyperpos[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			if (cur == qbest){
				if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) {best = cur; xorb = xorq;}
			}
		}

		ExOp::toLeftFlood(xorb);
        ExOp::toMax(xorc);


		for(cur= LEAD;cur <= best;cur++){
       //     printf("%X and %X\n" , max[cur], (q_max.k.hyperpos[cur] & (xorb ^xorc)) );
            if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
        xorb = xorb >> 1;
		for(;cur < DIMS;cur++){
        //    printf("%X and %X\n" , max[cur], (q_max.k.hyperpos[cur] & (xorb ^xorc)) );
		 	if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
	//	printf("TRUE!\n");
		return SETCMP_EQUAL;

	}
	LFHTEMP SETCMP_enum HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::operator()(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> &query) const{
		// check if the intervals are: disjoint, intersect, or contained!
        Tuple<unsigned int, 3> minmax = query.k.getMin();

        for(unsigned int i=0; i<DIMS;i++) if( minmax[i] > max[i]) return SETCMP_DISJOINT;
        minmax = query.k.getMax();
        for(unsigned int i=0; i<DIMS;i++) if( minmax[i] < min[i]) return SETCMP_DISJOINT;

        return SETCMP_EQUAL;
	/*	STORETYPE xorb;
		STORETYPE xorc;
		int best;
		int cur=0;
		for(cur=0;cur<LEAD;cur++) if ((q_min.k.hyperpos[cur]) !=  q_max.k.hyperpos[cur]) break;
		if (cur < LEAD){
			// some day... fix this
			return SETCMP_DISJOINT; // disjoint... for now
		}
		//	if (DIMS == LEAD) (...) min == max

		best = DIMS -1;
		xorb = (q_min.k.hyperpos[best]) ^ (q_max.k.hyperpos[best]);

		for(cur = DIMS-2; cur >= LEAD;cur--){
			xorc = (q_min.k.hyperpos[cur]) ^ (q_max.k.hyperpos[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

		// most distant direction found!
		// xorb contain most distant bit;

		if (sizeof(STORETYPE) >= 4) xorb |= (xorb >> 16);
		if (sizeof(STORETYPE) >= 2) xorb |= (xorb >> 8);
		xorb |= (xorb >> 4);
		xorb |= (xorb >> 2);
		xorb |= (xorb >> 1);
        if (sizeof(STORETYPE) == 4) xorc = 0xFFFFFFFF;
        else if (sizeof(STORETYPE) == 2) xorc = 0xFFFF;

		for(cur= LEAD;cur < best;cur++){
            if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
        xorb = xorb >> 1;
		for(;cur < DIMS;cur++){
			if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
		return SETCMP_EQUAL;
		*/
	}



    LFHTEMP void HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::initFromGaussElem(const GaussElem<Tuple<double, DIMS> >& input){
        unsigned int i;
        Tuple<double, DIMS> m = input.getMean();
        Trianglix<double, DIMS> cov = input.getVar_biaised();
        double r;
        for(i=0;i<DIMS;i++){
            r = sqrt(cov[i]);
            max[i] = (STORETYPE)(1367130551.1528631795645895562207f * atan(m[i] + r));
            min[i] = (STORETYPE)(1367130551.1528631795645895562207f * atan(m[i] - r));
        }
    }
    LFHTEMP void HyperPositionQuery_Cube<STORETYPE,DIMS,LEAD,NODES>::show(FILE*f, int level)const{
        fprintf(f, "HyperPositionQuery, Cube Area\n min");
        min.show(f,level+1);
        fprintf(f, " max ");
        max.show(f,level+1);
        fprintf(f, "\n");
    }
#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int  DIMS, unsigned int  LEAD,class NODES>
LFHTEMP bool Query_HyperCube<STORETYPE,DIMS,LEAD,NODES>::operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &q_min, const HyperCursor<STORETYPE,DIMS,LEAD> &q_max)const{
    STORETYPE xorb,xorc,xorq;
    int best;
    int cur=0;
    for(cur=0;cur<LEAD;cur++) if ((q_min.hyperpos[cur]) !=  q_max.hyperpos[cur]) break;
    if (cur < LEAD) return false; // some day... fix this// disjoint... for now
    //	if (DIMS == LEAD) (...) min == max

    int qbest;
    if (q_min.mag > q_max.mag) {xorq = 1 << (q_min.mag>>8); qbest = (q_min.mag & 0xFF);}
    else {xorq = 1 << (q_max.mag>>8); qbest = (q_min.mag & 0xFF);}

    cur = DIMS-1;
    best = DIMS-1;
    xorb = (q_min.hyperpos[cur]) ^ (q_max.hyperpos[cur]);
    if (cur == qbest){
        if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) xorb = xorq;
    }
    for(cur--; cur >= LEAD;cur--){
        xorc = (q_min.hyperpos[cur]) ^ (q_max.hyperpos[cur]);
        if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
        if (cur == qbest){
            if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) {best = cur; xorb = xorq;}
        }
    }
    if (sizeof(STORETYPE) >= 4) xorb |= (xorb >> 16);
    if (sizeof(STORETYPE) >= 2) xorb |= (xorb >> 8);
    xorb |= (xorb >> 4);
    xorb |= (xorb >> 2);
    xorb |= (xorb >> 1);
    ExOp::toMax(xorc);

    for(cur= LEAD;cur <= best;cur++){
        if (max[cur] < (q_max.hyperpos[cur] & (xorb ^xorc))) return false;
        if (min[cur] > (q_max.hyperpos[cur] | xorb )) return false;
    }
    xorb = xorb >> 1;
    for(;cur < DIMS;cur++){
        if (max[cur] < (q_max.hyperpos[cur] & (xorb ^xorc))) return false;
        if (min[cur] > (q_max.hyperpos[cur] | xorb )) return false;
    }
return true;}
LFHTEMP bool Query_HyperCube<STORETYPE,DIMS,LEAD,NODES>::operator()(const HyperCursor<STORETYPE,DIMS,LEAD> &query)const{
    Tuple<unsigned int, 3> minmax = query.getMin();
    for(unsigned int i=0; i<DIMS;i++) if( minmax[i] > max[i]) return false;
    minmax = query.getMax();
    for(unsigned int i=0; i<DIMS;i++) if( minmax[i] < min[i]) return false;
return true;}
#undef LFHTEMP
#define LFHTEMP template<class STORETYPE,int DIMS,int LEAD,class NODES>
LFHTEMP SETCMP_enum HyperPosition_Contains<STORETYPE,DIMS,LEAD,NODES>::operator()(const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES> &query) const{
    return query.k.isContained(this->value) ? SETCMP_EQUAL : SETCMP_DISJOINT;
}
LFHTEMP SETCMP_enum HyperPosition_Contains<STORETYPE,DIMS,LEAD,NODES>::operator()(const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES>  &min, const KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES>  &max) const{
    printf("2v query:"); ExOp::show(min); ExOp::show(max);
    return SETCMP_EQUAL;
}

#undef LFHTEMP
#define LFHTEMP template<class S,unsigned int B>
LFHTEMP HyperPosition3D<S,B>& HyperPosition3D<S,B>::operator=(const HyperPosition<unsigned short, 3,0u> &other){
    S buf[HyperPosition3D<S,B>::SIZE];



    return *this;
}

#undef LFHTEMP
#define LFHTEMP template<class S, class N>

LFHTEMP void HyperPositionQuery_LinkedCubes<S,N>::setQuery(const Tuple<S,3u> *mima){
    c_mima[0] = mima[0];
    c_mima[1] = mima[1];
    c_mima[2] = mima[2];
    c_mima[3] = mima[3];
    for(int i=0;i<3;i++){
        min[i] = c_mima[(c_mima[0][i] < c_mima[2][i]) ? 0 : 2][i];
        max[i] = c_mima[(c_mima[1][i] > c_mima[3][i]) ? 1 : 3][i];
    }

    for(int i=0;i<3;i++){
        if (mima[0][(i+1)%3] == mima[2][(i+1)%3]){
            tan_gle[0][i] = 0.0f;
            tan_gle_cmp[0][i] = 0.0f;
            tan_gle[1][i] = 0.0f;
            tan_gle_cmp[1][i] = 0.0f;
        }else{
            if ((mima[0][i] == mima[2][i])||( (mima[0][(i+1)%3] < mima[2][(i+1)%3])^ (mima[0][i] < mima[2][i]))){
                tan_gle[0][i] = 0.0f;
                tan_gle_cmp[0][i] = 0.0f;
            }else{


            }
            if ((mima[1][i] == mima[3][i])||( (mima[0][(i+1)%3] < mima[2][(i+1)%3])^ (mima[1][i] < mima[3][i]))){
                tan_gle[1][i] = 0.0f;
                tan_gle_cmp[1][i] = 0.0f;
            }else{

            }

        }

        if (mima[1][(i+1)%3] == mima[3][(i+1)%3]){
            tan_gle[2][i] = 0.0f;
            tan_gle_cmp[2][i] = 0.0f;
            tan_gle[3][i] = 0.0f;
            tan_gle_cmp[3][i] = 0.0f;
        }else{
            if ((mima[0][i] == mima[2][i])||( (mima[1][(i+1)%3] < mima[3][(i+1)%3])^ (mima[0][i] < mima[2][i]))){
                tan_gle[2][i] = 0.0f;
                tan_gle_cmp[2][i] = 0.0f;
            }else{

            }
            if ((mima[1][i] == mima[3][i])||( (mima[1][(i+1)%3] < mima[3][(i+1)%3])^ (mima[1][i] < mima[3][i]))){
                tan_gle[3][i] = 0.0f;
                tan_gle_cmp[3][i] = 0.0f;
            }else{

            }

        }

    }
}

LFHTEMP SETCMP_enum HyperPositionQuery_LinkedCubes<S,N>::operator()(const KeyElem<HyperCursor<S,3u,0u>, N> &query) const{
    Tuple<unsigned int, 3> mima = query.k.getMin();
    for(unsigned int i=0; i<3;i++) if( mima[i] > max[i]) return SETCMP_DISJOINT;
    mima = query.k.getMax();
    for(unsigned int i=0; i<3;i++) if( mima[i] < min[i]) return SETCMP_DISJOINT;
    return SETCMP_EQUAL;
}
LFHTEMP SETCMP_enum HyperPositionQuery_LinkedCubes<S,N>::operator()(const KeyElem<HyperCursor<S,3u,0u>, N> &q_min, const KeyElem<HyperCursor<S,3u,0u>, N>  &q_max) const{

		S xorb, xorc, xorq;
		int best;
		int cur=0;


		int qbest;
		if (q_min.k.mag > q_max.k.mag) {xorq = 1 << (q_min.k.mag>>8); qbest = (q_min.k.mag & 0xFF);}
        else {xorq = 1 << (q_max.k.mag>>8); qbest = (q_min.k.mag & 0xFF);}

		cur = 2;
		best = 2;
		xorb = (q_min.k.hyperpos[cur]) ^ (q_max.k.hyperpos[cur]);
		if (cur == qbest){
			if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) xorb = xorq;
		}
		for(cur--; cur >= 0;cur--){
			xorc = (q_min.k.hyperpos[cur]) ^ (q_max.k.hyperpos[cur]);
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
			if (cur == qbest){
				if ((xorb <= xorq)&&((xorb & xorq) <= (xorb ^ xorq))) {best = cur; xorb = xorq;}
			}
		}
		/*
        if (best >= DIMS){ // one box contains the other...
            best -= DIMS;
			printf("b%i: M:%X querybased\n", best, xorb);
        } else*/
	//	printf("b%i: M:%X\n", best, xorb);
		// most distant direction found!
		// xorb contain most distant bit;

		if (sizeof(S) >= 4) xorb |= (xorb >> 16);
		if (sizeof(S) >= 2) xorb |= (xorb >> 8);
		xorb |= (xorb >> 4);
		xorb |= (xorb >> 2);
		xorb |= (xorb >> 1);
        ExOp::toMax(xorc);


		for(cur= 0;cur <= best;cur++){
       //     printf("%X and %X\n" , max[cur], (q_max.k.hyperpos[cur] & (xorb ^xorc)) );
            if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
        xorb = xorb >> 1;
		for(;cur < 3;cur++){
        //    printf("%X and %X\n" , max[cur], (q_max.k.hyperpos[cur] & (xorb ^xorc)) );
		 	if (max[cur] < (q_max.k.hyperpos[cur] & (xorb ^xorc))) return SETCMP_DISJOINT;
            if (min[cur] > (q_max.k.hyperpos[cur] | xorb )) return SETCMP_DISJOINT;
		}
	//	printf("TRUE!\n");
		return SETCMP_EQUAL;

}


LFHTEMP void HyperPositionQuery_LinkedCubes<S,N>::show(FILE*f, int level)const{

}


#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES, unsigned int COMP>
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::demote(const HyperCursor<STORETYPE,DIMS,LEAD>& what){


	unsigned int ql = partition.hierarchy.find_index(what);
    KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmp;
    tmp.k = what;
    if (ql){
        tmp.d = partition.hierarchy.tree[ql].first.d;
        tmp.k.par_mag = partition.hierarchy.tree[ql].first.k.par_mag;
        partition.hierarchy.remove(tmp);
        while(true){
            if (tmp.k.mag == 0) { ql = 0; break;}
        tmp.k.toLeftChild();
        ql = partition.hierarchy.find_index(tmp.k);
        tmp.k.toBrother();
        unsigned int qr = partition.hierarchy.find_index(tmp.k);
      //  printf("%c\t%c\n", (ql) ? 'Y' : 'N', (qr) ? 'Y' : 'N' );
        if (ql){
           partition.hierarchy.tree[ql].first.k.par_mag = tmp.k.par_mag;
           if (qr){
                partition.hierarchy.tree[qr].first.k.par_mag = tmp.k.par_mag;
                break;
           }
        }else{
           if (qr){
                partition.hierarchy.tree[qr].first.k.par_mag = tmp.k.par_mag;
                tmp.k.toBrother();
           }else tmp.k.toParent(); break;
        }
        }
        if (ql == 0) partition.hierarchy.insert(tmp);
    }
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){
	HyperCursor<STORETYPE,DIMS,LEAD> hcur(pos);
	Vector< HyperCursor<STORETYPE,DIMS,LEAD> > in_nodes;
    // Step 1: clear smaller boxes

	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > bi;
	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ba;
	 bi.k =  pos.getMinBox();
	 ba.k =  pos.getMaxBox();
     partition.hierarchy.removeRange(bi.k,ba.k);


	// Step 2: check brother, if same type, upgrade it and exit

	ba.k = hcur.getBrother();
	unsigned int qq = partition.hierarchy.find_index(ba.k);
	if (qq){ // use its parent!
		 bi.k.par_mag = partition.hierarchy.tree[qq].first.k.par_mag;

	if (partition.hierarchy.tree[qq].first.d == i_n) {

		hcur.toParent();
        ba.d = i_n;
		partition.hierarchy.remove(ba);

        typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite(partition.hierarchy, ba.k.getMinBox());
        for(; (daite.isValid())&&(ba.k.getMaxBox() >= (*daite).k); ++daite ){
               if ((*daite).k.par_mag == ba.k.mag) in_nodes.push_back((*daite).k);
            }

	 ba.k = hcur.getBrother();
     daite.findGE(ba.k);
	 while((daite.isValid())&&(ba.k == (*daite).k)){

		if ((*daite).d != i_n) break;
        // deleting brother, need its hierarchical childs
 		partition.hierarchy.remove(ba);

        for(daite.findGE(ba.k.getMinBox()); (daite.isValid())&&(ba.k.getMaxBox() >= (*daite).k); ++daite ){
               if ((*daite).k.par_mag == ba.k.mag) in_nodes.push_back((*daite).k);
            }


		hcur.toParent();

            ba.k = hcur.getBrother();
	 		daite.findGE(ba.k);

	 		}

	 ba.k = hcur;
	 ba.k.par_mag = bi.k.par_mag;
	 ba.d = i_n;

	partition.hierarchy.insert(ba);

	// Update the parent pointer of ex-brother subboxes
        for(qq=0;qq<in_nodes.size();qq++){
            //printf("updateParent!\n"); fflush(stdout);
            //in_nodes[qq].show();
            daite.findGE(in_nodes[qq]);
            if ((!daite.isValid())||((*daite).k != in_nodes[qq] )) {printf("Hierarchical tree is corrupted!\n"); LFH_exit(1);}
            partition.hierarchy.orderpreserve_deref(daite).k.par_mag = ba.k.mag;
        }
    ba.k.toParent();
    if (ba.k.mag == ba.k.par_mag) demote(ba.k);
	return;
	}
	}

    // Step 3: find parent
    ba.k = pos;
    bi.k = pos;
    bi.d = i_n;
    // bi.k.par_mag = 0xFFFF;
    bi.k.par_mag = partition.par_mag_spacepartition(bi.k);

    // check if parent is the same, dont inster a thing!

    // step 4 check if the parent is the same!
    if (bi.k.par_mag != 0xFFFF){
    ba.k = bi.k;
    ba.k.setMagnitude(bi.k.par_mag);
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator da_ite(partition.hierarchy,ba.k);
    if (ba.k.mag != bi.k.par_mag) {printf("erroneous paren magnitude!\n");}
    if ((*da_ite).d == bi.d) {/*printf("Meaningless insert!\n");*/return;}


    partition.hierarchy.insert(bi);
    bi.k.toParent();
    if (bi.k.mag == bi.k.par_mag) demote(bi.k);

    }else{
        ExOp::toZero(ba.d);
        if (ba.d != bi.d) partition.hierarchy.insert(bi);
        else {/*printf("Meaningless null insert!\n");*/return;}
    }
        if (! partition.checkHierarchySimple()) {printf("ERROR IN HIERARCHICAL!\n");fflush(stdout);}
	}
/*
//SETCMP_enum gagahaga(const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &mi,const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &ma );

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::insert_partition(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){

	HyperCursor<STORETYPE,DIMS,LEAD> hcur(pos);

	vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > in_nodes;

	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > bi;
	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ba;

	 bi.k =  pos.getMinBox();
	 ba.k =  pos.getMaxBox();
	// hcur.show();

	//  printf("%c\t%c\n",hcur.getBrother() > b[0].k ? 'Y' : 'N', hcur.getBrother() < b[1].k ? 'Y' : 'N');

	 // Step 1: remove sub or equal boxes!

 //    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite =
 //   ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(bi.k);
 //   if(ite.isValid()){
//		 if ((*ite).par_mag< bi.k.par_mag) bi.k.par_mag = (*ite).par_mag;
//		 ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ite);
//    }

    ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->removeRange(bi,ba);

//	 this->intersectionInterval(in_nodes,bi,ba);
//	 bi.k.par_mag =65535;
//	 for(int i=0 ;i< in_nodes.size();i++) {
//		 if ((*ite).par_mag< bi.k.par_mag) bi.k.par_mag = (*ite).par_mag;
//		 ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ite);
//	 }

	 // Step 2: check brother, if same type, upgrade it and exit

	 ba.k = hcur.getBrother();
	 unsigned int qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_index(ba.k);
	if (qq){ // use its parent!
		 bi.k.par_mag = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.k.par_mag;

	if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.d == i_n) {
		hcur.toParent();
		((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);



	 ba.k = hcur.getBrother();
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(ba);
	 while((daite.isValid())&&(ba.k == (*daite).k)){


		 if ((*daite).d != i_n) break;
		((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);
		hcur.toParent();

	 	ba.k = hcur.getBrother();

	 		 daite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(ba);

	 		}

	 ba.k = hcur;
	 ba.k.par_mag = bi.k.par_mag;

	 ba.d = i_n;
	((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->insert(ba);



	// TODO, update the parent pointer of ex-brother subboxes

	return;
	}
	}



	 // Step 3: find parent



	// Step 4: if parent is 1 layer up, downgrade/erase it, ( and update parent

//	if (bi.k.par_mag != 65535) {
//	ba.k = bi.k.getRealParent();
//	qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find(ba);

//	if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.d != i_n) return; // ministep if parent exist and is of same type, exit

//	if (bi.k.isRealParentDirect()){
		// todo
	//	}
//	}

	// Step 5: insert box at last!

	 	bi.d = i_n;
		bi.k = pos;
        bi.k.par_mag = 0xFFFF;
//	 if (in_nodes.size() == 0){
		 // we did not find anything, need to recover hierac parent!

		 // TODO!

//		 }

	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->insert(bi);


	}
//SETCMP_enum gagahaga(const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &mi,const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &ma );
*/
LFHTEMP typename SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::Segment SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::mkSegment(const Tuple<STORETYPE, DIMS>& min,  const Tuple<STORETYPE, DIMS>& max) const{ SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::Segment fout;
	HyperRectKeyIterator<STORETYPE, DIMS> rect_ite;	rect_ite.setMinandMax(min, max);
	typename LFHPrimitive::SpacePartition<STORETYPE, DIMS, LEAD, NODES>::KeyIterator box_ite = this->getKeyIterator();
	//printf("Making segment!\n");
	//ExOp::show(min);
	//ExOp::showf(max);

	if (rect_ite.first()) do{
		if (box_ite.first(rect_ite())) do {
			fout.push_back(box_ite());
		}while(box_ite.next(rect_ite()));
	}while(rect_ite.next());


	return fout;
}
LFHTEMP unsigned int SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::batchInsert_safe(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const Vector<KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES > >& data, unsigned int optional_first_index){
	unsigned int ite = optional_first_index;
	if (ite == 0u){
		HyperCursor<STORETYPE,DIMS, LEAD> cursor = pos.getMinBox();
		for(;ite<data.getSize();ite++){
			if (data[ite].k >= cursor) break;
		}
	}
	if (ite >= data.getSize()) return ite;
	if (pos == data[ite].k){
		this->insert(data[ite].k, data[ite].d);
		ite++;
	}else{
		KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES> da_void; ExOp::toZero(da_void.d);
		da_void.k = pos;
		HyperCursor<STORETYPE,DIMS, LEAD> da_max = pos.getMaxBox();
		this->insert(pos, da_void.d);
		for(;ite<data.getSize();ite++){
			if (data[ite].k > da_max) break;
			this->insert(data[ite].k, data[ite].d);
		}
	}
	// needs to fix par_mag!
	return ite;
}
LFHTEMP unsigned int SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::batchInsert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const Vector<KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES > >& data, unsigned int optional_first_index){
	unsigned int ite = optional_first_index;
	if (ite == 0u){
		HyperCursor<STORETYPE,DIMS, LEAD> cursor = pos.getMinBox();
		for(;ite<data.getSize();ite++){
			if (data[ite].k >= cursor) break;
		}
	}
	if (ite >= data.getSize()) return ite;
	if (pos == data[ite].k){
		this->insert(data[ite].k, data[ite].d);
		ite++;
	}else{
		RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > new_part;
		HyperCursor<STORETYPE,DIMS, LEAD> da_max = pos.getMaxBox();
		HyperCursor<STORETYPE,DIMS, LEAD> cursor = pos.getMinBox();
		KeyElem<HyperCursor<STORETYPE,DIMS, LEAD>, NODES> da_void; ExOp::toZero(da_void.d);
		typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator da_ite = new_part.find_first(cursor);
		for(;ite<data.getSize();ite++){
			if (data[ite].k > da_max) break;
			da_void.k = data[ite].k.commonContainer(cursor);
			if (da_void.k != data[ite].k){ // some zero found!
				da_void.k.toLeftChild();
				new_part.insert(da_void);
			}else{

			}
		}
		partition.hierarchy.overwriteRange(pos.getMinBox(), da_max, pos.getMaxBox(),new_part);
	}
	// needs to fix par_mag!
	return ite;
}
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::detect_artifacts() const{

    typename RBTofDoom<KeyElem<HyperCursor<STORETYPE, DIMS, LEAD>,  NODES > >::Iterator ite = partition.hierarchy.first();

    HyperCursor<STORETYPE, DIMS, LEAD> cur;
    for(;ite.isValid();++ite){
        cur = ite->k.getMinBox();
            typename RBTofDoom<KeyElem<HyperCursor<STORETYPE, DIMS, LEAD>,  NODES > >::Iterator ite2(partition.hierarchy, cur);
            if (ite2 != ite){ printf("testing!\n");
                do{
                    if (cur < ite2->k.getMinBox()) {++ite2;
                        if (ite2 == ite) ++ite2;
                    }else{
                    cur = ite2->k.getTiny_Next();
                    if (cur > ite->k.getMaxBox()) break; // useless here!
                    ite2.findGE(cur);
                    if (ite2 == ite) ++ite2;
                    }
                } while(ite2.isValid()); // exit?alright some space remained
                if (cur > ite->k.getMaxBox()) return true;
            } // else first half is empty, all good!
        }
    return false;
    }
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::remove_artifacts(){

    typename RBTofDoom<KeyElem<HyperCursor<STORETYPE, DIMS, LEAD>,  NODES > >::Iterator ite(partition.hierarchy);
	ite.findFirst();

    HyperCursor<STORETYPE, DIMS, LEAD> cur;
    for(;ite.isValid();++ite){
        cur = ite->k.getMinBox();
            typename RBTofDoom<KeyElem<HyperCursor<STORETYPE, DIMS, LEAD>,  NODES > >::Iterator ite2(partition.hierarchy, cur);
            if (ite2 != ite){ printf("testing!\n");
                do{
                    if (cur < ite2->k.getMinBox()) {++ite2;
                        if (ite2 == ite) ++ite2;
                    }else{
                    cur = ite2->k.getTiny_Next();
                    if (cur > ite->k.getMaxBox()) break; // useless here!
                    ite2.findGE(cur);
                    if (ite2 == ite) ++ite2;
                    }
                } while(ite2.isValid()); // exit?alright some space remained
                if (cur > ite->k.getMaxBox()) {
                    partition.hierarchy.updateParMag(ite->k,ite->k.par_mag, ite->k.mag );
                    partition.hierarchy.remove(ite);
                }
            } // else first half is empty, all good!

        }
    }
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::partitionIntersection(const HyperCursor<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const{

    // for now return what is in interval;

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > query;



}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::fillvoid(){ // proceedure for initialization

    const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> minimum;
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->first();

    HyperCursor<STORETYPE,DIMS,LEAD> other;


    Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > redo;
    unsigned short mmag,tmag;

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmpin;
    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmpin2;

    while(ite.isValid()){

        typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2 = ite;

        tmpin = *ite;
        ++ite;
        --ite2;
        if ((ite.isValid())&&((*ite).d == tmpin.d)){
            //  multiple blocks!
            redo.push_back(tmpin);
            while((ite.isValid())&&((*ite).d == tmpin.d)) {
                tmpin2 = *ite;
                ++ite;
            }
            if ((!ite.isValid())&&(!ite2.isValid())){
                // all blocks are of the same type!
                LFH_exit(1);
            }

            if (ite.isValid()){
                mmag = tmpin.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite2).k);
                }
            }

            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
            if  (mmag > tmpin.k.mag) {
                tmpin.k.setMagnitude(mmag);
            }
            redo.push_back(tmpin);

            if (ite.isValid()){
                mmag = tmpin2.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin2.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin2.k.commonContainer_mag((*ite2).k);
                }
            }

            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);

            if  (mmag > tmpin2.k.mag) {
                tmpin2.k.setMagnitude(mmag);
            }

            if (tmpin.k != tmpin2.k) {
                redo.push_back(tmpin2);
                // need to fill the space between  tmpin2 and  tmpin


                while(true){
                ++(tmpin.k);
                if (ite.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite).k);

                    if (ite2.isValid()){
                        tmag = tmpin.k.commonContainer_mag((*ite2).k);
                    if (tmag < mmag) mmag = tmag;
                }
                }else{
                    if (ite2.isValid()){
                        mmag = tmpin.k.commonContainer_mag((*ite2).k);
                    }else break;
                }
                HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
                if  (mmag > tmpin.k.mag) {
                    tmpin.k.setMagnitude(mmag);
                    redo.push_back(tmpin);
                } else break;


                }

                while(true){
                --(tmpin2.k);

                if (ite.isValid()){
                    mmag = tmpin2.k.commonContainer_mag((*ite).k);

                    if (ite2.isValid()){
                        tmag = tmpin2.k.commonContainer_mag((*ite2).k);
                    if (tmag < mmag) mmag = tmag;
                }
                }else{
                    if (ite2.isValid()){
                        mmag = tmpin2.k.commonContainer_mag((*ite2).k);
                    }else break;
                }

                HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
                if  (mmag > tmpin2.k.mag) {
                    tmpin2.k.setMagnitude(mmag);
                    if (tmpin.k == tmpin2.k) break;
                    else redo.push_back(tmpin2);
                } else break;


                }
            }

        }else{
            // single block!

            if (ite.isValid()){
                mmag = tmpin.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite2).k);
                }else{redo.push_back(tmpin); continue;}
            }
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
            if  (mmag > tmpin.k.mag) {
                tmpin.k.setMagnitude(mmag);
            }
            redo.push_back(tmpin);
        }
    }

    ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->batchInit(redo);

}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::fillvoid_directional(int direction, NODES before){ // use connectivity to fill the complete volumes, can only propagate in the input direction
    // fill void using the elems on the top

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > cur = this->max;
    Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > f_out;
    bool isPure;
    KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tcur;

    KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES > cmpress_cur;


    unsigned int i,j;
    j=0;

    NODES majority[32 * DIMS];

    cmpress_cur.k.toMax();
    cmpress_cur.k.hyperpos.getOrderAndLeadHyper(i,j); printf("%i\t%i\n",i,j);
/*
    cmpress_cur.k.show();
unsigned short tmp_mag = cmpress_cur.k.commonContainer_mag(cur.k); cmpress_cur.k.setMagnitude(tmp_mag);
    this->insert_partition(cmpress_cur.k, tcur.d);
    cmpress_cur.k.show();*/

    for(i=0;i<100000000;i++){






        if ((i % 100000) == 0) printf("tic\n");
        tcur.k = cur.k;
        tcur.k.toLower(2);

        if (cur.k.hyperpos[2] > tcur.k.hyperpos[2]){
     //   this->intersection(cur.k.hyperpos,f_out);
     //   if (f_out.size() != 0) break;
        this->insert_partition(tcur.k, cur.d);
        }
        /*
        tcur.k = cur.k;
        tcur.k.toHigher(2);
        tcur.d = this->getAt(tcur.k,isPure);
        if ((isPure)&&(ExOp::isZero(tcur.d))){
            this->insert_partition(tcur.k, tcur.d);
        }*/

        typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = this->find_first(cur.k);

        if (ite.isValid()) {


        if (((*ite).k != cur.k)&&((*ite).k.commonContainer_mag(cur.k) == (*ite).k.mag)){
        }else {
            --ite;
            if (!ite.isValid()) break;
        }
        /*
        if ((*ite).k == cur.k){
            --ite;
            if (!ite.isValid()) break;
        }*/
        cur = *ite;
        } else cur = this->max;


    }

/*

    const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> minimum;
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)&_input)->last();

    HyperCursor<STORETYPE,DIMS,LEAD> other;


    Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > redo;
    unsigned short mmag,tmag;
*/
}
LFHTEMP NODES SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::getAt(const HyperPosition<STORETYPE,DIMS,LEAD> posis, bool &isPure) const{NODES f_out;
	SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator kite(*this);
    isPure = true;
    if (kite.first(posis)){
        if (kite().k.hyperpos == posis){
        f_out = kite().d;
        }else isPure = false;
    }else ExOp::toZero(f_out);

    return(f_out);
	}


    // get the separation of *pure* block, stores oriantation of contact in par_mag
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::save_portion(const HyperCursor<STORETYPE,DIMS,LEAD>& box, FILE *f) const{

}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::load_portion(const HyperCursor<STORETYPE,DIMS,LEAD>& box, FILE *f){

}

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::getContacts( Vector< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &out, HyperCursor<STORETYPE,DIMS,LEAD> &where)const{

}
/*
LFHTEMP unsigned int SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::getContainerOf(const HyperPosition<STORETYPE,DIMS,LEAD> posis) const{
	KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,NODES> query;  query.k = posis;
	unsigned int prev,next;
	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,NODES> >*)this)->findPN(query,prev,next);
	return 0;
	}*/
LFHTEMP KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES> SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::getAtPoint(const Tuple<STORETYPE,DIMS> &posis) const{ KeyElem< HyperCursor<STORETYPE,DIMS,LEAD> , NODES> fout;
    fout.k.hyperpos[0] = posis[0] | 1;
    for(unsigned int i=1;i<DIMS;i++) fout.k.hyperpos[i] = posis[i];

    fout.k.mag = 0;
    bool isPure;

    unsigned short tmag = this->partition.par_mag_spacepartition(fout.k);
    if (tmag != 0xFFFF){
    fout.k.setMagnitude(tmag);
    fout.d = this->getAt(fout.k, isPure);
    } else ExOp::toZero(fout.d);
    return fout;
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::intersection(const HyperPosition<STORETYPE,DIMS, LEAD> &pos, vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &f_out) const{

//	int i;
	Tuple<STORETYPE, DIMS> min[2];

	min[0][0] = 0; min[0][1] = 0; min[0][2] = 0;
	min[1][0] = 0; min[1][1] = 15; min[1][2] = 0;

	HyperPositionQuery_Cube<STORETYPE, DIMS, LEAD,  NODES > cq( min[0], min[1]);

	(( RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> >const *const )this)->intersection(f_out, cq );



	}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::intersection(const HyperPosition<STORETYPE,DIMS, LEAD> &pos, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const{


	Tuple<STORETYPE, DIMS> min[2];

	min[0][0] = 0; min[0][1] = 0; min[0][2] = 0;
	min[1][0] = 0; min[1][1] = 15; min[1][2] = 0;

	HyperPositionQuery_Cube<STORETYPE, DIMS, LEAD, NODES > cq( min[0], min[1]);

	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,  NODES > >const *const )this)->intersection(f_out, cq );



	}
	// SpacePartition private routines
LFHTEMP unsigned short SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::get_parent_order(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>& from){
    // get closest Item
     KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmp_key(from.k.getMinBox(), from.d);
     typename RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key.k);
    unsigned short f_out;

	tmp_key.k = from.k.getMaxBox();
    if (ite.isValid()){
        if ((*ite) < tmp_key){
            // found an item within box, iterate up!
        //    if (*ite).
            while( (*ite).k.par_mag != ExCo<unsigned short>::mkMaximum() ){
                if ((*ite).k.par_mag > from.k.mag) break;
                tmp_key.k.setMagnitude((*ite).k.par_mag);
                ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key.k);
                if ((!ite.isValid())||((*ite).k.mag != tmp_key.k.mag)){
                    printf("Should never happen!\n");
                }
            }
            return( (*ite).k.par_mag );
        }

    }else{
        if (((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->size == 0) return(ExCo<unsigned short>::mkMaximum());
        ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->last();
    }

    // check knowing that no elements exists inside

    // the one pointed to, and the the precepding are to find the best quandidates if existing

    tmp_key.k = (*ite).k;
    while(!tmp_key.k.strictly_contains(from.k)){


    }

    typename RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2 = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key.k);

    return(0);
}
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::first(){
        const HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >& tree = this->target->partition;
        if (tree.hierarchy.size == 0) {
            (*this).curkey.k.toLargestBoxPossible();
            ExOp::toZero((*this).curkey.d);
            return true;
        }

        HyperCursor<STORETYPE,DIMS,LEAD> smallest; smallest.toMin();

        typename HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >::Iterator iten(tree.hierarchy); iten.findFirst();
        if (!iten.isValid()) return false;
        unsigned short mag = smallest.commonContainer_mag((*iten).k);


        if (mag == (*iten).k.mag) {
            (*this).curkey = (*iten);
                ++iten;
                if (iten.isValid()){
                    mag = (*this).curkey.k.commonContainer_mag((*iten).k);
                    if (mag == (*this).curkey.k.mag) (*this).curkey.k.toLeftChild();
                }
        }else{ // iterator in preceeding void
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mag);
            smallest.setMagnitude(mag);
            (*this).curkey.k = smallest;
            ExOp::toZero((*this).curkey.d);
        }
        return true;
    }
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
        const HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >& tree = this->target->partition;
        if (tree.hierarchy.getSize() == 0) {
            (*this).curkey.k = in_box;
            ExOp::toZero((*this).curkey.d);
            return true;
        }
        //printf("first\n"); fflush(stdout);
        HyperCursor<STORETYPE,DIMS,LEAD> smallest = in_box.getMinBox();
        typename HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >::Iterator iten(tree.hierarchy,smallest);

        unsigned short mag;

        if (iten.isValid()){

        mag = smallest.commonContainer_mag((*iten).k);
        if (mag == (*iten).k.mag) {
            (*this).curkey = (*iten);
            if (mag >= in_box.mag) (*this).curkey.k = in_box;
            else{
                ++iten;
                if (iten.isValid()){
                    mag = (*this).curkey.k.commonContainer_mag((*iten).k);
                    if (mag == (*this).curkey.k.mag) (*this).curkey.k.toLeftChild();
                }
            }
        }else{
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mag);
            if (mag > in_box.mag) mag = in_box.mag;
            smallest.setMagnitude(mag);
            unsigned short magn= tree.par_mag_spacepartition(smallest);
            (*this).curkey.k = smallest;
            if (0xFFFF !=  magn){ //printf("path B\n");
                 smallest.setMagnitude(magn);
                 iten.findGE(smallest);
                 (*this).curkey.d = (*iten).d;
            }else{ //printf("path C\n");
            ExOp::toZero((*this).curkey.d);
            }}
        } else{
            // mag = smallest.commonContainer_mag(tree.max.k);
            mag = tree.par_mag_spacepartition(smallest);
            if (0xFFFF != mag){
                (*this).curkey.k = in_box;
                smallest.setMagnitude(mag);
                iten.findGE(smallest);
                (*this).curkey.d = (*iten).d;
            }else{
                (*this).curkey.k = in_box;
                ExOp::toZero((*this).curkey.d);
            }
        }
        return true;
    }
// might be not working. depending on filter
LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box,const FUNCTOR &filter){
        const HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >& tree = this->target->partition;
        if (tree.getSize() == 0) {
            (*this).curkey.k = in_box;
            ExOp::toZero((*this).curkey.d);
            return filter((*this).curkey);
        }

    //    printf("first\n"); fflush(stdout);
        HyperCursor<STORETYPE,DIMS,LEAD> smallest = in_box.getMinBox();
        typename HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >::Iterator iten(tree.hierarchy,smallest);

     //   if (iten.isValid()){
     //       printf("%c:", filter(*iten) ? 'Y' : 'N' ); ExOp::show(*iten);
     //   }

       /* while((iten.isValid())&&(!filter(*iten))) {
            ++iten;
            if ((!iten.isValid())||(filter(*iten))) break;
            ++iten;
            if ((!iten.isValid())||(filter(*iten))) break;
            ++iten;
            if ((!iten.isValid())||(filter(*iten))) break;
            ++iten;
            if ((!iten.isValid())||(filter(*iten))) break;
            if (in_box.getMaxBox() >= (*iten).k) {iten.toInvalid(); break; }
            ++iten;
        }*/

        while((iten.isValid())&&(!filter(*iten))) ++iten;

        unsigned short mag;
        (*this).curkey.k = in_box;
      //  if (!filter((*this).curkey)) {printf("self filter!\n"); in_box.show();}


        if (iten.isValid()){
      //  printf("id?\n");
      //  (*iten).show();
      //  in_box.show();

        mag = smallest.commonContainer_mag((*iten).k);
      //  printf("%X\t%X mag  (%i)\n", mag, (*iten).k.mag, (*iten).d);
        if (mag == (*iten).k.mag) { // found box contains smallest
            (*this).curkey = (*iten);
            if (mag >= in_box.mag) (*this).curkey.k = in_box;
            else{ // found box is smaller than query, maybe is not full
                ++iten;
                while((iten.isValid())&&(!filter(*iten))) ++iten;
                if (iten.isValid()){
                    mag = (*this).curkey.k.commonContainer_mag((*iten).k);
                    if (mag == (*this).curkey.k.mag) (*this).curkey.k.toLeftChild();
                }
            }
          //  printf("path A\n");
        }else{ // paritalbox, needs to find parent
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mag);
            if (mag > in_box.mag) mag = in_box.mag;
            smallest.setMagnitude(mag);
            unsigned short magn= tree.par_mag_spacepartition(smallest);
            (*this).curkey.k = smallest;
            if (0xFFFF !=  magn){// printf("path B\n");
                 smallest.setMagnitude(magn);
                 iten.findGE(smallest);
                 (*this).curkey.d = (*iten).d;
            }else{ //printf("path C\n");
            ExOp::toZero((*this).curkey.d);
            }}
        } else{ // no valid found area, needs to find parent of area
            // mag = smallest.commonContainer_mag(tree.max.k);
            mag = tree.par_mag_spacepartition(smallest);
            while(mag < in_box.mag){
                smallest.setMagnitude(mag);
                mag = tree.par_mag_spacepartition(smallest);
            }
            if (0xFFFF != mag){
                (*this).curkey.k = in_box;
                smallest.setMagnitude(mag);
                iten.findGE(smallest);
                (*this).curkey.d = (*iten).d;
            } else {
                (*this).curkey.k = in_box;
                ExOp::toZero((*this).curkey.d);// printf("gothere!\n");
                return filter((*this).curkey);
            }
        }

        if (!filter((*this).curkey)) {
       //     printf("first relay!\n");fflush(stdout);
            if (!this->next(filter)) { // empty then... try zero
                (*this).curkey.k = in_box;
                ExOp::toZero((*this).curkey.d);
                return filter((*this).curkey);
            }return (in_box.getMaxBox() >= (*this).curkey.k);
        }


     //   printf("first out!\n");fflush(stdout);
        return (in_box.getMaxBox() >= (*this).curkey.k);
    }
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::next(){
    const HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >& tree = this->target->partition;
    if (tree.hierarchy.getSize() == 0) return false;
    //printf("next!\n");fflush(stdout);
    HyperCursor<STORETYPE,DIMS,LEAD> smallest;
    smallest = (*this).curkey.k.getTiny_Next();

    if (smallest < (*this).curkey.k) return(false);

    typename HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >::Iterator iten(tree.hierarchy,smallest);

    unsigned short magn = (*this).curkey.k.commonContainer_mag(smallest);
    HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(magn);

    if (iten.isValid()){
     // printf("valid!\n");fflush(stdout);
        unsigned short magp = smallest.commonContainer_mag((*iten).k);
        if (magp <= magn) {
            if (magp == (*iten).k.mag) {
                (*this).curkey = (*iten);
                ++iten;
                if (iten.isValid()){
                    magp = (*this).curkey.k.commonContainer_mag((*iten).k);
                    if (magp == (*this).curkey.k.mag) (*this).curkey.k.toLeftChild();
                }
                return(true);
            } else {magn = magp;HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(magn);}
        }

    } //printf("validout!\n");fflush(stdout);
    smallest.setMagnitude(magn);
    (*this).curkey.k = smallest;
    magn = tree.par_mag_spacepartition(smallest);//printf("parvalidout!\n");fflush(stdout);

    if (magn == 0xFFFF) ExOp::toZero((*this).curkey.d);
    else {
    smallest.setMagnitude(magn);
    iten.findGE(smallest);
        if (iten.isValid()) { (*this).curkey.d = (*iten).d; if ((*iten).k.mag != magn) printf("WTF!!!!\n");}
        else {ExOp::toZero((*this).curkey.d);printf("WTF!\n");}
    }
    return true;
}
LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::next(const FUNCTOR& filter ){
    const HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >& tree = this->target->partition;
    if (tree.hierarchy.size == 0) return false;

    HyperCursor<STORETYPE,DIMS,LEAD> smallest;
   // printf("next!\n");fflush(stdout);

    do{
    smallest = (*this).curkey.k.getTiny_Next();
    if (smallest < (*this).curkey.k) return(false);
    typename HierarchicalTree< STORETYPE, DIMS,  LEAD, NODES >::Iterator iten(tree.hierarchy,smallest);

    while((iten.isValid())&&(!filter(*iten))) ++iten;

    unsigned short magn = (*this).curkey.k.commonContainer_mag(smallest);
    HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(magn);

    if (iten.isValid()){
    //  printf("valid!\n");fflush(stdout);
        unsigned short magp = smallest.commonContainer_mag((*iten).k);
        if (magp <= magn) {
            if (magp == (*iten).k.mag) {
                (*this).curkey = (*iten);
                ++iten;
                while((iten.isValid())&&(!filter(*iten))) ++iten;
                if (iten.isValid()){
                    magp = (*this).curkey.k.commonContainer_mag((*iten).k);
                    if (magp == (*this).curkey.k.mag) (*this).curkey.k.toLeftChild();
                }
                if (!filter((*this).curkey)) continue;

             //   printf("next out!\n");fflush(stdout);
                return(true);
            } else {magn = magp;HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(magn);}
        }

    } //printf("validout!\n");fflush(stdout);
    smallest.setMagnitude(magn);
    (*this).curkey.k = smallest;
    magn = tree.par_mag_spacepartition(smallest);//printf("parvalidout!\n");fflush(stdout);

    if (magn == 0xFFFF) ExOp::toZero((*this).curkey.d);
    else {
    smallest.setMagnitude(magn);
    iten.findGE(smallest);
        if (iten.isValid()) { (*this).curkey.d = (*iten).d; if ((*iten).k.mag != magn) printf("WTF!!!!\n");}
        else {ExOp::toZero((*this).curkey.d);printf("WTF!\n");}
    }
        if (filter((*this).curkey)) break;
    } while (true);
          //  printf("next out!\n");fflush(stdout);
    return true;
}
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
    if (!next()) return false;
    return (in_box.getMaxBox() >= (*this).curkey.k);
}
LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::KeyIterator::next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box,const FUNCTOR& filter){
    if (!next(filter)) return false;

    return (in_box.getMaxBox() >= (*this).curkey.k);
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::save_mk2(FILE *f) const{
		unsigned int buffer[1024];
		fwrite(&partition.size ,sizeof(unsigned int),1 ,f);
		if (partition.size == 0) return;

		for(unsigned int cur = partition.inorderFirst(buffer); cur != 0xFFFFFFFF;cur = partition.inorderNext(buffer)){
			ExOp::save(partition.tree[cur].first.k.hyperpos, f);
		}

		for(unsigned int cur = partition.inorderFirst(buffer); cur != 0xFFFFFFFF;cur = partition.inorderNext(buffer)){
			ExOp::save(partition.tree[cur].first.d, f);
		}
    }
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::load_mk2(FILE *f, unsigned int s){
		unsigned int c,i_size;
        this->toMemfree();
        if (1 != fread(&i_size ,sizeof(unsigned int),1 ,f)) {fprintf(stderr, "Corrupt File! cant load Space Partition!\n"); LFH_exit(1);}

        if (i_size == 0) return;
		partition.makeEmptyTree(i_size);

        // fixed known size!
        unsigned int buffer[1024];
        unsigned int cur;

        for(cur = partition.inorderFirst(buffer),c = i_size ;(c--) != 0; cur = partition.inorderNext(buffer)){
            ExOp::load(partition.tree[cur].first.k.hyperpos, f);
        }

        for(cur = partition.inorderFirst(buffer),c = i_size ;(c--) != 0; cur = partition.inorderNext(buffer)){
            ExOp::load(partition.tree[cur].first.d, f);
        }

        this->initMagnitudes_routine();
        this->initParMagnitudes_routine();

        max = partition.tree[cur].first;
        cur = partition.inorderFirst(buffer);
        min = partition.tree[cur].first;
    }
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::initMagnitudes_routine(){
    if (partition.size == 0) return;
    unsigned int i,s;
    unsigned int order,lead;

    s = partition.size - (1 << partition.alloc_mag) + partition.firstlone;

    for(i=1;i<= s;i++){
        partition.tree[i].first.k.hyperposr.getOrderAndLeadHyper(order,lead);
        partition.tree[i].first.k.mag = (order<< 8) | lead;
    }

    // lone leaves
    for(i = partition.firstlone +1; (i >> partition.alloc_mag) == 0; i++){
        partition.tree[i].first.k.hyperposr.getOrderAndLeadHyper(order,lead);
        partition.tree[i].first.k.mag = (order<< 8) | lead;
    }
    if ((i & 1) == 0) {
        partition.tree[i].first.k.hyperposr.getOrderAndLeadHyper(order,lead);
        partition.tree[i].first.k.mag = (order<< 8) | lead;
    }
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::clearParMagnitudes_routine(){
    if (partition.size == 0) return;
    unsigned int i,s;
    unsigned int order,lead;

    s = partition.size - (1 << partition.alloc_mag) + partition.firstlone;

    for(i=1;i<= s;i++){
        ExOp::toMax(partition.tree[i].first.k.par_mag);
    }

    // lone leaves
    for(i = partition.firstlone +1; (i >> partition.alloc_mag) == 0; i++){
        ExOp::toMax(partition.tree[i].first.k.par_mag);
    }
    if ((i & 1) == 0) {
        ExOp::toMax(partition.tree[i].first.k.par_mag);
    }
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::initParMagnitudes_routine(){
    if (partition.size == 0) return;
    this->clearParMagnitudes_routine();
    unsigned int buffer[256];
    unsigned int cur;
    typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator i_ite(partition);

    unsigned int c = partition.size;
    HyperCursor<STORETYPE,DIMS, LEAD> mima[2];
    for(cur = partition.inorderFirst(buffer) ;(c--) != 0; cur = partition.inorderNext(buffer)){
        mima[0] = partition.tree[cur].first.k.getMinBox();
        mima[1] = partition.tree[cur].first.k.getMaxBox();
        i_ite = partition.find_first(mima[0]);
    }
}
LFHTEMP SpacePartition<STORETYPE,DIMS,LEAD,NODES,COMP>& SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::operator=(const SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>& no_compression){
    this->initMagnitudes_routine();
    this->clearParMagnitudes_routine();
    return *this;
}

// Class Segment
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::Segment::shiftSegment(const Tuple<STORETYPE, DIMS> &value){
	for(unsigned int i=0;i<DIMS;i++) this->shiftSegment(value[i],i);
}
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::Segment::shiftSegment(STORETYPE value, unsigned int dir){
	if (this->getSize() == 0) return;
	unsigned int i;
	SpacePartition<STORETYPE, DIMS, LEAD, NODES,COMP>::Segment bufout;
	HeapTree< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > purbuf;
	HyperRectKeyIterator<STORETYPE, DIMS> rect_ite;
	Tuple<STORETYPE, DIMS> mima[2];

	unsigned int mag = value;
	mag |= mag << 1;
	mag |= mag << 2;
	mag |= mag << 4;
	if (sizeof(STORETYPE) > 1)  mag |= mag << 8;
	if (sizeof(STORETYPE) > 2)  mag |= mag << 16;
	mag = mag ^ (((STORETYPE)0)-1);


	HyperCursor<STORETYPE,DIMS,LEAD> cursor = (*this)[0].k.getMinBox();
	cursor.hyperpos[dir] += value;

	for(i=0;i<this->getSize();i++){
		if (((*this)[i].k.hyperpos[dir] & mag) == 0){ // no split, ez mode!
			(*this)[i].k.hyperpos[dir] += value;
			purbuf.insert((*this)[i]);
		}else{
			mima[0] = (*this)[i].k.getMin();
			mima[1] = (*this)[i].k.getMax();
			mima[0][dir]+= value;
			mima[1][dir]+= value;
			rect_ite.setMinandMax(mima[0], mima[1]);
			if (rect_ite.first()) do{
				purbuf.insert( KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES >( rect_ite(), (*this)[i].d));
			}while(rect_ite.next());
		}
	}

	while(!purbuf.isEmpty()){
		bufout.push_back(purbuf.pop());
	}

	this->toMemmove(bufout);
}

#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES>

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES, 0u>::insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){
	HyperCursor<STORETYPE,DIMS,LEAD> hcur(pos);

    // Step 1: clear smaller boxes

	KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > bi;
	KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ba;

	bi.k =  pos.getMinBox();
	ba.k =  pos.getMaxBox();
    partition.removeRange(bi.k.hyperpos,ba.k.hyperpos);

	// Step 2: check brother, if same type, upgrade it and exit

	ba.k = hcur.getBrother();
	unsigned int qq = partition.find_index(ba.k.hyperpos);
	if (qq){ // use its parent!

        if (partition.tree[qq].first.d == i_n) { // same type, merging
            hcur.toParent();
            ba.d = i_n;

            partition.remove(ba);

            ba.k = hcur.getBrother();
            typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite(partition,ba.k);
            while((daite.isValid())&&(ba.k.hyperpos == (*daite).k)){

                if ((*daite).d != i_n) break;
                // deleting brother, need its hierarchical childs
                partition.remove( ba );

                hcur.toParent();

                ba.k = hcur.getBrother();
                daite.findGE(ba.k);

            }
            partition.insert( KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES >(hcur.hyperpos, i_n) );
            return;
        }
	}else{
        // need to check parent to cut down if different
        ba.k = hcur;
        typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite2(partition, ba.k);
        if (daite2.isValid()){
//            if daite2().

        }
        qq = partition.find_index(ba.k.hyperpos);
    }


    partition.insert( KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES >(pos.hyperpos, i_n) );
}

LFHTEMP KeyElem< HyperPosition<STORETYPE,DIMS,LEAD> , NODES> SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::operator()(const Tuple<STORETYPE,DIMS> &posis) const{
    typename RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = partition.find_first(posis);
    if (!ite.isValid()){
        return voidnode;
    }else if (ite->k == posis){
        return (*ite);
    }else{
        return voidnode;
    }
}
LFHTEMP KeyElem< HyperPosition<STORETYPE,DIMS,LEAD> , NODES>& SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::operator()(const Tuple<STORETYPE,DIMS> &posis){ // guarrantied to bepure
    typename RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(partition,posis);
    STORETYPE xora, xorb, xorc;
    unsigned int cur, best, aest;
    if (!ite.isValid()){ // *no* larger exists, use max as bound
        cur=0;
		best = DIMS -1;
		xorb = (partition.getMax().k[best]) ^ posis[best];
		for(cur = DIMS-2; cur+1u != LEAD;cur--){
			xorc = (partition.getMax().k[cur]) ^ posis[cur];
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}
        //printf("query %X\t%X\n", posis[0],posis[1]);
        //pri = ite->k.getMin(); printf("max mi %X\t%X\n", pri[0],pri[1]);
        //pri = ite->k.getMax(); printf("max ma %X\t%X\n", pri[0],pri[1]);
    }else if (ite->k == posis){
     //   printf("B\n");
        return ite.ordering_preserving_reference();
    }else{
		cur=0;
		best = DIMS -1;
		xorb = (ite->k[best]) ^ posis[best];
		for(cur = DIMS-2; cur+1u != LEAD;cur--){
			xorc = (ite->k[cur]) ^ posis[cur];
			if ((xorb <= xorc)&&((xorb & xorc) <= (xorb ^ xorc))) {best = cur; xorb = xorc;}
		}

       // Tuple<unsigned short, 2> pri;
        //printf("query %X\t%X\n", posis[0],posis[1]);
       // pri = ite->k.getMin(); printf("max mi %X\t%X\n", pri[0],pri[1]);
        //pri = ite->k.getMax(); printf("max ma %X\t%X\n", pri[0],pri[1]);

        --ite;
        if (ite.isValid()){
         //   pri = ite->k.getMin(); printf("min mi %X\t%X\n", pri[0],pri[1]);
         //   pri = ite->k.getMax(); printf("min ma %X\t%X\n", pri[0],pri[1]);
            cur=0;
            aest = DIMS -1;
            xora = (ite->k[best]) ^ posis[best];
            for(cur = DIMS-2; cur+1u != LEAD;cur--){
                xorc = (ite->k[cur]) ^ posis[cur];
                if ((xora <= xorc)&&((xora & xorc) <= (xora ^ xorc))) {aest = cur; xora = xorc;}
            }
            // min order wins!
            if ((xorb & xora) <= (xorb ^ xora)) {
                if (xora < xorb){best = aest; xorb = xora;}
            } else if (aest < best) best = aest;
        }
    }

    if (sizeof(STORETYPE) >= 4) xorb |= xorb >>  16;
    if (sizeof(STORETYPE) >= 2) xorb |= xorb >>  8;
    xorb |= xorb >> 4;
    xorb |= xorb >> 2;
    xorb |= xorb >> 1;

    //printf("%i and %X mask, %X \n", best, xorb, ExCo<STORETYPE>::mkMaximum());
    if (best == LEAD) {xorb >>= 1; best = DIMS -1;}
    else best--;
    xorb ^= ExCo<STORETYPE>::mkMaximum();

    for(cur=0;cur<best;cur++) voidnode.k[cur] =  posis[cur] & xorb;
    voidnode.k[cur] = (posis[cur] & xorb);
    xorb |= xorb >>1;
    voidnode.k[cur] |= xorb ^ (xorb <<1);
    for(cur++;cur<DIMS;cur++) voidnode.k[cur] =  posis[cur] & xorb;
//        pri = voidnode.k.getMin(); printf("res mi %X\t%X\n", pri[0],pri[1]);
//        pri = voidnode.k.getMax(); printf("res ma %X\t%X\n", pri[0],pri[1]);
    ExOp::toZero(voidnode.d);
    return voidnode;
}


LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::clearBoxArea(const HyperCursor<STORETYPE,DIMS, LEAD> &where){



}


LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::first(){
    trite.findFirst();
    if (trite.isValid()){
        (*this).curkey.k = trite->k;
        (*this).curkey.d = trite->d;
        return true;
        }else{
        return false;
    }
}

/*
LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
        const RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >& tree = this->target->partition;
        trite = tree.find_first(in_box.getMinBox());
        if ((trite.isValid())&&(trite->k <= in_box.getMaxBox())) {
            (*this).curkey.k = trite->k;
            (*this).curkey.d = trite->d;
            return true;
            }
        return false;
    }*/


LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
    const RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >& tree = this->target->partition;
    trite.findGE(in_box.getMinBox());
    if (trite.isValid()){
        if (trite->k <= in_box.getMaxBox()){
            (*this).curkey.k = trite->k;
            (*this).curkey.d = trite->d;
            return true;
        }
        if (trite->k.getMinBox() > in_box.getMaxBox()){
            --trite;
            if ((!trite.isValid())||(trite->k.getMaxBox() < in_box.getMinBox())) return false;
        }
        (*this).curkey.d = trite->d;

        //printf("mega!\n");
        //in_box.show();
        //trite->k.show();
    }else if ((tree.size != 0)&&(tree.getMax().k.getMaxBox() > in_box.getMinBox())){
        (*this).curkey.d = tree.getMax().d;
        //printf("megamax!\n");
        //in_box.show();
        //tree.max.k.show();
    } else return false;
    (*this).curkey.k = in_box;
    return true;
}

LFHTEMP NODES* SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::deref(){

}

LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::first(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box,const FUNCTOR &filter){
        const RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >& tree = this->target->partition;
        trite = tree.find_first(in_box.getMinBox());
        while((trite.isValid())&&(!filter(*trite))) ++trite;
        if ((!trite.isValid())||(trite->k > in_box.getMaxBox()))  return false;
        (*this).curkey.k = trite->k;
        (*this).curkey.d = trite->d;
        return true;
    }

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next(){
    ++trite;
    if (!trite.isValid()) return false;
    (*this).curkey.k = trite->k;
    (*this).curkey.d = trite->d;
    return true;
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next_withvoid(){
    ++trite;
    if (!trite.isValid()) return false;
    (*this).curkey.k = trite->k;
    (*this).curkey.d = trite->d;
    return true;
}

LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next(const FUNCTOR& filter ){
    ++trite; if (!trite.isValid()) return false;
    while(!filter(*trite)) {
        ++trite;
        if (!trite.isValid()) return false;
    }
    (*this).curkey.k = trite->k;
    (*this).curkey.d = trite->d;
    return true;
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
    if (!next()) return false;
    return (in_box.getMaxBox() >= (*this).curkey.k);
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next_withvoid(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box){
    if (!next_withvoid()) return false;
    return (in_box.getMaxBox() >= (*this).curkey.k);
}


LFHTEMP template<class FUNCTOR> bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box,const FUNCTOR& filter){
    if (!next(filter)) return false;
    return (in_box.getMaxBox() >= (*this).curkey.k);
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::isEmptyAt(const HyperPosition<STORETYPE,DIMS,LEAD> &in_box) const{
    typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = partition.find_first(in_box);
    if (ite.isValid()){
        if ((in_box.compare(ite->k) & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT) return false;
        --ite;
        if (!ite.isValid()) return true;
        return ((in_box.compare(ite->k) & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);
    }else return ((in_box.compare(partition.getMax().k) & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::isEmptyAt(const HyperCursor<STORETYPE,DIMS,LEAD> &in_box) const{
    typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = partition.find_first(in_box);
    if (ite.isValid()){
        if ((in_box.compare(ite->k) & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT) return false;
        --ite;
        if (!ite.isValid()) return true;
        return ((in_box.compare(ite->k) & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);
    }else return ((in_box.compare(partition.getMax().k) & SETCMP_DISJOINT_BIT) == SETCMP_DISJOINT_BIT);
}

LFHTEMP NODES SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::getAt(const HyperPosition<STORETYPE,DIMS,LEAD> posis, bool &isPure) const{ NODES f_out;

    typename RBTofDoom< KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite(partition, posis);
    printf("got it %i\n", partition.getSize());
    if (ite.isValid()){

        if ((posis.compare(ite->k) & SETCMP_DISJOINT_BIT) != SETCMP_DISJOINT_BIT){
            isPure = true;
             printf("A\n");
        }else{
            --ite;
            printf("B\n");
            if (!(ite.isValid())){
                isPure = true;printf("Ex2\n");
                ExOp::toZero<NODES>(f_out);
                return(f_out);
            }
        }
    }else{
        // TODO
       isPure = true;printf("Ex\n");
       ExOp::toZero<NODES>(f_out);
       return(f_out);
    }
    f_out = ite->d;printf("Done\n");
    return(f_out);
	}




#undef LFHTEMP
#define LFHTEMP template <class C,unsigned int nbdim>


LFHTEMP bool HyperRectKeyIterator<C, nbdim>::first(){
    unsigned int i;
    for(i=0;i<nbdim;i++) cur.hyperpos[i]= min[i];
    cur.hyperpos[0] |= 1; cur.mag =0; // smallest box
    unsigned short mag;
    mag = cur.largest_Container_mag_rect(min, max);
    if (mag > 0x8FFF) return(false);
    cur.setMagnitude(mag);
    return true;
}

LFHTEMP bool HyperRectKeyIterator<C, nbdim>::next(){
    // parent of current must intersect with box
    // so if current is locally a left node, either the whole right node is void or some is to be covered
    HyperCursor<C,nbdim,0u> smallest = cur.NextTinyInBox(min, max);

    if (smallest < cur) return false; // overflowed!

    if (smallest.mag == 0xFFFF) return(false);
    unsigned short mag = smallest.commonContainer_mag(cur); // max magnitude, due to preceding
    HyperCursor<C,nbdim,0u>::decr_mag(mag);
    unsigned short tmag = smallest.largest_Container_mag_rect(min, max);
    cur = smallest;
    cur.setMagnitude((tmag < mag) ? tmag : mag);
    return true;
}

LFHTEMP bool HyperRectKeyIterator<C, nbdim>::find_first(const HyperCursor<C,nbdim,0u > & item){ // find first box in 3d rect that contains "item"
    unsigned int i;
    for(i=0;i<nbdim;i++) cur.hyperpos[i]= min[i];
    cur.hyperpos[0] |= 1; cur.mag =0; // smallest box
    unsigned short mag;
    mag = cur.largest_Container_mag_rect(min, max);
    if (mag == 0xFFFF) return(false);
    cur.setMagnitude(mag);
    return true;
}


#undef LFHTEMP
#define LFHTEMP template<class STORETYPE>
LFHTEMP unsigned short HeirarchicalSurface3<STORETYPE>::expand(HyperCursor<STORETYPE,3u,0u>& box)const{
    return(box.mag);
}
LFHTEMP HyperCursor<STORETYPE,3u,0u> HeirarchicalSurface3<STORETYPE>::compress(HyperCursor<STORETYPE,3u,0u>& box)const{
    return(box);
}

#undef LFHTEMP
#define LFHTEMP template<class N, class S, unsigned int D, unsigned int L>

LFHTEMP bool HierarchicalMergable<N,S,D,L>::findContainerOf_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const HyperCursor<S,D,L> where){
	//printf("search\n");
	ite = space.find_first(where);
	Tuple<S,D> pos[2];
	bool no_intersect;
	if (!ite.isValid()){
		//printf("hit pos end...\n");
		ite = space.find_last(where);
		if (!ite.isValid()) return true;
		else if (where.isMagnitudeLE(ite->k)) return false;
		no_intersect = (ite->k < where.getMinBox());
	}else{
		if (where.isContainedBy(ite->k)) return false;
		no_intersect = (ite->k > where.getMaxBox());
		--ite;
	//	printf("trying previous\n");
		if (!ite.isValid()) return no_intersect;
		if (where.isContainedBy(ite->k)) return false;
		no_intersect = no_intersect & (ite->k < where.getMinBox());
	}
	//printf("Neither contains!\n");
	ite.toInvalid();
	return no_intersect;
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::wrGetAt(HyperCursor<S,D,L>& _out, N& _out_data, const Tuple<S,D> &where) const{
    for(uint32_t i=0;i<D;i++) _out.hyperpos[i] = where[i];
    _out.mag =0;
    _out.hyperpos[0] |=1;
    getAt_update(_out_data,_out);
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::getAt(N& fout, const HyperCursor<S,D,L> &where) const{
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(space,where);
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(fout); return true;}
		fout = ite->d;
		return (where.getMinBox() > ite->k);
	}
	if (where.isContainedBy(ite->k)) {fout = ite->d; return true;}
	Tuple<S,D> limit = ite->k.getMin();
	HyperCursor<S,D,L> cursor;
	--ite;
	if (ite.isValid()) {
		fout = ite->d;
		if (where.isContainedBy(ite->k)) return true;
		if (ite->k >= where.getMinBox()) return false;
	} else ExOp::toZero(fout);
	cursor = where;
	cursor.toMinBoxExclusiveContainer(limit);
	if (where.isContainedBy(cursor.hyperpos)) return true;
	return false;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::getAt_simplify(N& fout, const HyperCursor<S,D,L> &where) const{
	if (!this->getAt(fout,where)) return false;
	switch(fout.planeIntersection(where)){
		case 1: fout.toZero(); break;
		case 2: fout.toOne(); break;
		default: break;
	}
return true;}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::getAt_simplifyMK2(N& fout, const HyperCursor<S,D,L> &where) const{
	if (!this->getAt(fout,where)) return false;
	if (fout.canSimplify()) fout.simplifyAt(fout,where);
return true;}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::getAt_simplifyMK2(N& fout, const HyperCursor<S,D,L> &where, const M& merger) const{
	if (!this->getAt(fout,where)) return false;
	if (fout.canSimplify()) merger.simplifyAt(fout, fout,where);
return true;}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::getAt_update(N& fout, HyperCursor<S,D,L> &where) const{
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(space,where); // init iterator and findGE
	Tuple<S,D> limit;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(fout); where.toLargestBoxPossible();return true;}
		fout = ite->d;
		if (where.getMinBox() > ite->k){ // 2019june swap
			limit = ite->k.getMin();
			where.toLargestBoxPossible();
			where.toMinBoxExclusiveContainer(limit);
			return true;
		}else return false;
	}
	if (where.isContainedBy(ite->k)) {fout = ite->d; where = HyperCursor<S,D,L> (ite->k); return true;}
	limit = ite->k.getMin();
	HyperCursor<S,D,L> cursor = where;
	--ite;
	if (ite.isValid()) {
		fout = ite->d;
		where = HyperCursor<S,D,L> (ite->k);
		if (cursor.isContainedBy(ite->k)) return true;
		if (ite->k >= cursor.getMinBox()) return false;
	} else {ExOp::toZero(fout);where.toLargestBoxPossible();}
	where.toMinBoxExclusiveContainer(limit);
	do{
		if (cursor.isContainedBy(where.hyperpos)) {
            N tmp;
            if (fout.simplifyAt(tmp, where)) fout = tmp;
            return true;
        }
		if (!where.findNext()) return false;
		where.toMinBoxExclusiveContainer(limit);
	}while(where.mag != 0xFFFF);
	return false;
}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::getAt_update(N& fout, HyperCursor<S,D,L> &where, const M& merger) const{
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(space,where); // init iterator and findGE
	Tuple<S,D> limit;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(fout); where.toLargestBoxPossible();return true;}
		fout = ite->d;
		if (where.getMinBox() > ite->k){ // 2019june swap
			limit = ite->k.getMin();
			where.toLargestBoxPossible();
			where.toMinBoxExclusiveContainer(limit);
			return true;
		}else return false;
	}
	if (where.isContainedBy(ite->k)) {fout = ite->d; where = HyperCursor<S,D,L> (ite->k); return true;}
	limit = ite->k.getMin();
	HyperCursor<S,D,L> cursor = where;
	--ite;
	if (ite.isValid()) {
		fout = ite->d;
		where = HyperCursor<S,D,L> (ite->k);
		if (cursor.isContainedBy(ite->k)) return true;
		if (ite->k >= cursor.getMinBox()) return false;
	} else {ExOp::toZero(fout);where.toLargestBoxPossible();}
	where.toMinBoxExclusiveContainer(limit);
	do{
		if (cursor.isContainedBy(where.hyperpos)) {
		    N tmp;
            if (merger.simplifyAt(fout, tmp, where)) fout = tmp;
            return true;
        }
		if (!where.findNext()) return false;
		where.toMinBoxExclusiveContainer(limit);
	}while(where.mag != 0xFFFF);
	return false;
}




/*
LFHTEMP bool HierarchicalMergable<D,S,L,N>::getAt(N& fout, HyperCursor<S,D,L> &where) const{
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(&(this->space));
	ite = space.find_last(where);
	if (!ite.isValid()){
		ite = space.first();
		ExOp::toZero(fout);
		if (!ite.isValid()) return true;
		return (where.getMaxBox() < ite->k);
	}
	if (where.isContainedBy(ite->k)) {fout = ite->d; return true;}
	++ite;
	if (!ite.isValid()) return false;
	if (where.isContainedBy(ite->k)) {fout = ite->d; return true;}
	return false;
}*/

// if trivial, must desagree with preceding box, if non-trivial, an intersection with plane of interest must exist in the *first* box

LFHTEMP void HierarchicalMergable<N,S,D,L>::simplifyCheck(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite){
	typename HierarchicalMergable<N,S,D,L>::KeyIterator kite = this->getKeyIterator();
	if ((!kite.first(HyperCursor<S,D,L>(ite->k)))||(kite().d.isTrivial())) {printf("simplify should not be called!\n");LFH_exit(1);}
	N tmp = kite().d;
	uint32_t border = 0;
	do{
		border |= kite().d.planeIntersection(kite().k);
		if (border == 3u) return;
	}while((kite.next())&&(kite().d == tmp));
	// should be simplified!
//	printf("Performing a simplification!\n"); fflush(stdout);

	if (border == 2)  ExOp::toZero(tmp);
	else ExOp::toOne(tmp);

	--ite;
	if (!ite.isValid()){ ite = space.first();
		if (border == 2) space.remove_update(ite);
		else {ite.ordering_preserving_reference().d = tmp; ++ite;}
	}else{
		if (ite->d == tmp) {++ite; space.remove_update(ite);}
		else {++ite;ite.ordering_preserving_reference().d = tmp;++ite;}
	}

	if (ite.isValid()) if (ite->d == tmp) space.remove(ite);

}
LFHTEMP void HierarchicalMergable<N,S,D,L>::insertAfter_routine(){


}
LFHTEMP void HierarchicalMergable<N,S,D,L>::changeNodeSize_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const Tuple<S,D> &minpos){
	HyperCursor<S,D,L> result = ite->k;
	result.toMinBoxExclusiveContainer(minpos);
	if (result.mag == 0xFFFF) {
		//printf("reducing to nothingness!\n");
		space.remove(ite); ite.findLE(result); return;
	}
	N prev = ite->d;
	if (prev.isTrivial()) {ite.ordering_preserving_reference().k = result.hyperpos; return;}
	unsigned int tmppl = prev.planeIntersection(result);
	//printf("super reduction time: %i at:\n", tmppl);
	//ExOp::show(result);

	if (tmppl == 3u) {ite.ordering_preserving_reference().k = result.hyperpos; return;}
	// complex becomes simple
	N newins;
	newins.toTrivial(tmppl);
	--ite;
	if (!ite.isValid()){
		//printf("hit an invalid...!\n");
		N tmptmp; tmptmp.toZero();
		ite.findFirst();
		if (newins == tmptmp) {space.remove(ite); // zero in front, delete and beware!
			ite.findFirst();
			while((ite.isValid())&&(ite->d == tmptmp)) space.remove(ite);
			ite.toInvalid();
		}else{
			ite.ordering_preserving_reference().k = result.hyperpos;
			ite.ordering_preserving_reference().d = newins;
		}
		return;
	}else{
		if ((ite->d == newins)||(newins.trivialIsCoherrentWith(ite->d.planeIntersection(result)))){
			++ite;
			space.remove(ite);
			ite.findLE(result);
		}else{
			++ite;
			ite.ordering_preserving_reference().k = result;
			ite.ordering_preserving_reference().d = newins;
		}
	}
	if (!result.findNext()) {
            printf("SHOULD NEVER BE A!\n"); LFH_exit(1);
    }
	result.toMinBoxExclusiveContainer(minpos);
	while (result.mag != 0xFFFF){
		unsigned int tmppl2 = prev.planeIntersection(result);
		if (tmppl2 == 3u) {
			space.insert(KeyElem<HyperPosition<S,D,L>, N >(result, prev));
			ite.findGE(result);
			return;
		}

		if (tmppl != tmppl2){
			tmppl = tmppl2;
			newins.toTrivial(tmppl);
			space.insert(KeyElem<HyperPosition<S,D,L>, N >(result,newins));
			ite.findGE(result);
		}

		if (!result.findNext()) {
                printf("SHOULD NEVER BE B!!!\n"); LFH_exit(1);
        }
		result.toMinBoxExclusiveContainer(minpos);
	}
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::changeNodeSizeMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const Tuple<S,D> &minpos, bool debug){
	HyperCursor<S,D,L> cursor = ite->k;
	cursor.toMinBoxExclusiveContainer(minpos);
	if (debug) {printf("Deflating:\n"); minpos.show();}
	N cur = ite->d;
	N prev;

	if (cursor.mag == 0xFFFF) {
		//printf("reducing to nothingness!\n");
		space.remove(ite); ite.findLE(cursor);
		return;
	}
    if ((cursor.hyperpos[0] == 0xA400)&&(cursor.hyperpos[1] == 0x8200)&&(cursor.hyperpos[2] == 0x8300)){
        minpos.show();
        debug = true; cur.show(); printf("Old Support\n"); HyperCursor<S,D,L>(ite->k).show(); printf("To be:\n"); cursor.show();
    }

	N simplyfied[2];

	if (!cur.simplifyAt(simplyfied[0],cursor)) {ite.ordering_preserving_reference().k = cursor.hyperpos; return;}

	if (debug){
        printf("Deflating:\n");
        cur.show();
        simplyfied[0].show();
	}

	// got simplified, can be eaten by previous, maybe
	--ite;
	if (ite.isValid()){
        prev = ite->d;
        if (debug){ printf("previous is:\n"); HyperCursor<S,D,L>(ite->k).show(); ite->d.show();}
        ++ite;
    }else {prev.toZero();ite.findFirst();}

    if (debug) ite->d.show();

	while(true){
        prev.simplifyAt(simplyfied[1], cursor);
        if (simplyfied[0] != simplyfied[1]) break;
        if (!cursor.findNext()) {
                printf("SHOULD NEVER BE C!\n"); LFH_exit(1);
        }
        cursor.toMinBoxExclusiveContainer(minpos);

        if (cursor.mag == 0xFFFF){
            // fully consumed...
            space.remove(ite);
            ite.findLE(cursor);
            return;
        }
        if (!cur.simplifyAt(simplyfied[0], cursor)){
                ite.ordering_preserving_reference().k = cursor.hyperpos;return; // next do not simplify, wont be eaten
        }
	}

    if (debug) {printf("Replacement:\n"); cursor.show(); simplyfied[0].show();}
	ite.ordering_preserving_reference().k = cursor.hyperpos;
	ite.ordering_preserving_reference().d = simplyfied[0];
    prev = simplyfied[0];
    HyperCursor<S,D,L> precur = cursor;
	// fix till mnipos
	while(cursor.findNext()){
	    cursor.toMinBoxExclusiveContainer(minpos);
        if (cursor.mag == 0xFFFF) break;
        if (!cur.simplifyAt(simplyfied[1], cursor)){
            space.insert(KeyElem<HyperPosition<S,D,L>, N >(cursor, cur));
            ite.findGE(cursor);
            return;
        }
        prev.simplifyAt(simplyfied[0], cursor);
        if (simplyfied[1] != simplyfied[0]){
            space.insert(KeyElem<HyperPosition<S,D,L>, N >(cursor, simplyfied[1]));
            prev = simplyfied[1];
            precur = cursor;
        }
	}
	ite.findGE(precur);
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::changeNodeSizeMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const Tuple<S,D> &minpos, const M& merger, bool debug){
	HyperCursor<S,D,L> cursor = ite->k;
	cursor.toMinBoxExclusiveContainer(minpos);
	if (debug) printf("Deflating:\n");
	N cur = ite->d;
	N prev;

	if (cursor.mag == 0xFFFF) {
		//printf("reducing to nothingness!\n");
		space.remove(ite); ite.findLE(cursor);
		return;
	}


	N simplyfied[2];

	if (!merger.simplifyAt(cur,simplyfied[0],cursor)) {ite.ordering_preserving_reference().k = cursor.hyperpos; return;}

	if (debug){
        printf("Deflating:\n");
        cur.show();
        simplyfied[0].show();
	}

	// got simplified, can be eaten by previous, maybe
	--ite;
	if (ite.isValid()){prev = ite->d;++ite;}
	else {prev.toZero();ite.findFirst();}

	while(true){
        merger.simplifyAt(prev,simplyfied[1], cursor);
        if (simplyfied[0] != simplyfied[1]) break;
        if (!cursor.findNext()) {
                printf("Surprise in deflating:\n");
                cur.show();
                simplyfied[0].show();
                prev.show();
                printf("SHOULD NEVER BE D!\n"); LFH_exit(1); // does makee= sense?
        }
        cursor.toMinBoxExclusiveContainer(minpos);

        if (cursor.mag == 0xFFFF){
            // fully consumed...
            space.remove(ite);
            ite.findLE(cursor);
            return;
        }
        if (!merger.simplifyAt(cur,simplyfied[0], cursor)){ite.ordering_preserving_reference().k = cursor.hyperpos;return;}
	}


	ite.ordering_preserving_reference().k = cursor.hyperpos;
	ite.ordering_preserving_reference().d = simplyfied[0];
    prev = simplyfied[0];
    HyperCursor<S,D,L> precur = cursor;
	// fix till mnipos
	while(cursor.findNext()){
	    cursor.toMinBoxExclusiveContainer(minpos);
        if (cursor.mag == 0xFFFF) break;
        if (!merger.simplifyAt(cur,simplyfied[1], cursor)){
            space.insert(KeyElem<HyperPosition<S,D,L>, N >(cursor, cur));
            ite.findGE(cursor);
            return;
        }
        merger.simplifyAt(prev,simplyfied[0], cursor);
        if (simplyfied[1] != simplyfied[0]){
            space.insert(KeyElem<HyperPosition<S,D,L>, N >(cursor, simplyfied[1]));
            prev = simplyfied[1];
            precur = cursor;
        }
	}
	ite.findGE(precur);
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::writePrevNext_routine(const HyperCursor<S,D,L> &where, const N& what,typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, N &prevprop, N &nextprop){
	HyperCursor<S,D,L> cursor;	cursor = where.getMaxBox();
	ite.findLE(cursor);

	if (!ite.isValid())	{ // is empty and first box
		ite.findFirst();
		if (!ite.isValid())	{
			//printf("hh h haha!\n"); fflush(stdout);
			nextprop.toZero();
			prevprop.toZero();
		}else{
			//printf("hh h haha2!\n");
			if (where.isContainedBy(ite->k)){
				nextprop = ite->d;
				prevprop = ite->d;
				//printf("da routine!\n"); fflush(stdout);
				changeNodeSize_routine(ite,where.getMin());// reduces box!
			//	this->show();
				//printf("da routine done!\n"); fflush(stdout);

			}else{
				nextprop.toZero();
				prevprop.toZero();
			}
			++ite;
		}
	}else{
		//printf("hh h hehe!\n");
		//printf("has ttfound!\n");
		// need to check if neighbor overlap with current contains query...
		//ite->k.show();
		++ite;
		if ((ite.isValid())&&(where.isContainedBy(ite->k))){
			nextprop = ite->d;
			//printf("hh h hoho! %X\n", ite->d());
			changeNodeSize_routine(ite,where.getMin());// reduces box!

			//space.show();
		}else{
			//printf("hh h hihi!\n", ite.isValid()?'Y':'N');
			if (!ite.isValid()) ite.findLast();
			else --ite;
			nextprop = ite->d;
			if (ite->k == where.hyperpos){
				if (ite->d == what) {/*printf("detected pointless!\n");*/return true;} // pointless insert!
				space.remove_mm(ite);
			}else if (where.isContainedBy(ite->k)){
				//printf("got here! large\n");
				changeNodeSize_routine(ite,where.getMin());// reduces box!
			}else{
			//	printf("range rem here!\n");
			//	where.show();
			//	ite->k.show();
				cursor = where.getMinBox();
				space.removeRange(cursor.hyperpos, ite); // points to last that is not removed!
			}
		}
		if (ite.isValid()) {prevprop = ite->d; ++ite;}
		else {prevprop.toZero(); ite.findFirst();}
	}
	return false;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::writePrevNextMK2_routine(const HyperCursor<S,D,L> &where, const N& what,typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, N &prevprop, N &nextprop, bool debug){
	HyperCursor<S,D,L> cursor;	cursor = where.getMaxBox();
	ite.findLE(cursor);

	if (!ite.isValid())	{ // is empty and first box
		ite.findFirst();
		if (!ite.isValid())	{ // empty world!
			nextprop.toZero();
			prevprop.toZero();
		}else{
			//printf("hh h haha2!\n");
			if (where.isContainedBy(ite->k)){
				nextprop = ite->d;
				//printf("da routine!\n"); fflush(stdout);
				changeNodeSizeMK2_routine(ite,where.getMin(),debug);// reduces box!
				if (ite.isValid()) {prevprop = ite->d;++ite;}
				else {prevprop.toZero(); ite.findFirst();}
			//	this->show();
				//printf("da routine done!\n"); fflush(stdout);
			}else{
				nextprop.toZero();
				prevprop.toZero();
				++ite;
			}
		}
	}else{
		//printf("hh h hehe!\n");
		//printf("has ttfound!\n");
		// need to check if neighbor overlap with current contains query...
		//ite->k.show();
		++ite;
		if ((ite.isValid())&&(where.isContainedBy(ite->k))){
			nextprop = ite->d;
			//printf("hh h hoho! %X\n", ite->d());
			changeNodeSizeMK2_routine(ite,where.getMin(),debug);// reduces box!
		}else{
			//printf("hh h hihi!\n", ite.isValid()?'Y':'N');
			if (!ite.isValid()) ite.findLast();
			else --ite;
			nextprop = ite->d;
			if (ite->k == where.hyperpos){
				if (ite->d == what) {/*printf("detected pointless!\n");*/return true;} // pointless insert!
				space.remove_mm(ite);
			}else if (where.isContainedBy(ite->k)){
				//printf("got here! large\n");
				changeNodeSizeMK2_routine(ite,where.getMin(),debug);// reduces box!
			}else{
			//	printf("range rem here!\n");
			//	where.show();
			//	ite->k.show();
				cursor = where.getMinBox();
				space.removeRange(cursor.hyperpos, ite); // points to last that is not removed!
			}
		}
		if (ite.isValid()) {prevprop = ite->d; ++ite;}
		else {prevprop.toZero(); ite.findFirst();}
	}
	return false;
}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::writePrevNextMK2_routine(const HyperCursor<S,D,L> &where, const N& what,typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, N &prevprop, N &nextprop,  const M& merger, bool debug){
	HyperCursor<S,D,L> cursor;	cursor = where.getMaxBox();
	ite.findLE(cursor);

	if (!ite.isValid())	{ // is empty and first box
		ite.findFirst();
		if (!ite.isValid())	{ // empty world!
			nextprop.toZero();
			prevprop.toZero();
		}else{
			//printf("hh h haha2!\n");
			if (where.isContainedBy(ite->k)){
				nextprop = ite->d;
				//printf("da routine!\n"); fflush(stdout);
				changeNodeSizeMK2_routine(ite,where.getMin(),merger,debug);// reduces box!
				if (ite.isValid()) {prevprop = ite->d;++ite;}
				else {prevprop.toZero(); ite.findFirst();}
			//	this->show();
				//printf("da routine done!\n"); fflush(stdout);
			}else{
				nextprop.toZero();
				prevprop.toZero();
				++ite;
			}
		}
	}else{
		//printf("hh h hehe!\n");
		//printf("has ttfound!\n");
		// need to check if neighbor overlap with current contains query...
		//ite->k.show();
		++ite;
		if ((ite.isValid())&&(where.isContainedBy(ite->k))){
			nextprop = ite->d;
			//printf("hh h hoho! %X\n", ite->d());
			changeNodeSizeMK2_routine(ite,where.getMin(),merger,debug);// reduces box!
		}else{
			//printf("hh h hihi!\n", ite.isValid()?'Y':'N');
			if (!ite.isValid()) ite.findLast();
			else --ite;
			nextprop = ite->d;
			if (ite->k == where.hyperpos){
				if (ite->d == what) {/*printf("detected pointless!\n");*/return true;} // pointless insert!
				space.remove_mm(ite);
			}else if (where.isContainedBy(ite->k)){
				//printf("got here! large\n");
				changeNodeSizeMK2_routine(ite,where.getMin(),merger,debug);// reduces box!
			}else{
			//	printf("range rem here!\n");
			//	where.show();
			//	ite->k.show();
				cursor = where.getMinBox();
				space.removeRange(cursor.hyperpos, ite); // points to last that is not removed!
			}
		}
		if (ite.isValid()) {prevprop = ite->d; ++ite;}
		else {prevprop.toZero(); ite.findFirst();}
	}
	return false;
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::HyperPositionInflate_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite){
	HyperCursor<S,D,L> cursor = ite->k;
	bool _debug_ = false;
	N toexp = ite->d;
	unsigned short mag;
	if (_debug_) {printf("Inflating %X: ",toexp());	cursor.show();}
	++ite;
	cursor.toLargeMinBoxContainer();

	bool hasNext = ite.isValid();
	Tuple<S,D> minOfNext;
	uint32_t tmppl;

	if (hasNext) {
		if (_debug_) {printf("NextLimit %X: ", ite->d()); HyperCursor<S,D,L>(ite->k).show();}
		minOfNext = HyperCursor<S,D,L>(ite->k).getMin();
		cursor.toMinBoxExclusiveContainer(minOfNext);
		--ite;
	}else ite.findLast();
	--ite;
	do{
		if (cursor.mag == 0xFFFF) printf("should be impossible... got intersecting nodes! fdadf\n");
		HyperCursor<S,D,L> brother = cursor.getBrother();

		if (_debug_) printf("inner! continue %c\n", (brother <= cursor)  ? 'Y' : 'N');
		if (brother > cursor) break;

		tmppl = toexp.planeIntersection(brother);


		if (ite.isValid()){
		//	if (toexp == ite->d){ //
		//		cursor = ite->k;
		//		++ite;
		//		space.remove_mm(ite);
		//	}else{

			HyperCursor<S,D,L> prev(ite->k);
			if (_debug_) {printf("PrevFound: ");  prev.show();}
			//ite->k.show();
			mag = brother.commonContainer_mag(prev);
			if (_debug_) {
				printf("Cur: ");cursor.hyperpos.show();
				printf("mag for cur... %X\n", cursor.mag);
				printf("Cur Bro: ");  brother.show();
				printf("PE %X %X %X\n", mag, brother.mag, prev.mag);
			}
			if (mag > brother.mag){
				if (_debug_) printf("HH\n");
				if (tmppl == 3u) break;
				unsigned int tmppl2 = ite->d.planeIntersection(brother);
				if (tmppl2 == tmppl){
					if (_debug_) printf("aggree... %c\n", (toexp.isTrivial())? 'Y' : 'N');
					if (toexp.isTrivial()){

						if (ite->d.planeIntersection(cursor) == tmppl){
						//	printf("LOOOOPY!\n");
						//	ite->d.show();
						//	cursor.show();
							cursor.toParent();
							return writeHyperPosition_routine(cursor, ite->d);
						}
					}
				} else if ((tmppl2 != 3u)&&(toexp.isTrivial())){
					N tmptmp; tmptmp.toTrivial(tmppl2);
					tmptmp.TrivialMerge(cursor, toexp);
					cursor.toParent();
					if (_debug_) printf("did an trivial merge path3! %X\n", toexp());
					return writeHyperPosition_routine(cursor, toexp);
				} else {if (_debug_) printf("disaggree...\n"); break;}
				cursor.toParent();
				if (_debug_) printf("PFF\n");
			}else if (brother.hyperpos == ite->k){
				if (_debug_) printf("PFCdd\n");
				if (ite->d.isTrivial()){
					if (ite->d.trivialIsCoherrentWith(tmppl)){
					}else if (toexp.isTrivial()) {
						ite->d.TrivialMerge(cursor, toexp);
						++ite; ite.ordering_preserving_reference().k = cursor.hyperpos;
						if (_debug_) printf("did an trivial merge! %X\n", toexp());
						cursor.toParent();
						return writeHyperPosition_routine(cursor, toexp);
					}else break;
					cursor.toParent();
					space.remove_mm(ite);
					if ((ite.isValid())&&(ite->d == toexp)) {
						// previous matches! use it
						cursor = ite->k;
						++ite;
						if (_debug_) printf("auto remove previous, matches!\n");
						space.remove_mm(ite);
						--ite;
					}
				}else{
					// previous is brother, hence ite->k.planeIntersection(brother) == 3
					 // DONE!
					break;
				}
			}else break;
		//	}
		}else if (tmppl == 3u) break;
		else{
			// previous is non-existant (trivial 0)
			// *TODO*
			//HyperCursor<S,D,L> prev;
			printf("PFWTF\n");
			break;
		}
		//cursor.show();
		cursor.toLargeMinBoxContainer();
		if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
		if (cursor.mag == 0xFFFF){
			printf("got cornered!\n");
			ite->k.show();
			ite->d.show();
			if (hasNext) {printf("nextis:"); ExOp::show(minOfNext);}
			else printf("no next\n");
			LFH_ALIVE; LFH_exit(1);
		}
	}while(true);
	if (ite.isValid()) ++ite;
	else ite.findFirst();
	//printf("to %X\n", cursor.mag);
	ite.ordering_preserving_reference().k = cursor.hyperpos;
	++ite;
	if (!ite.isValid()) return;
	if (toexp.isTrivial()){


		HyperCursor<S,D,L> ncursor = HyperCursor<S,D,L>(ite->k);
		//if (ite->d.isTrivial()){
		//	if (ite->d.trivialIsCoherrentWith(toexp.planeIntersection(ncursor))) return writeHyperPosition_routine(ncursor, toexp);
		//}

		if (_debug_) {printf("recurmerge!\n"); ncursor.show();}
		HyperCursor<S,D,L> brother = ncursor.getBrother();
		if ((brother > ncursor)||(brother.getMinBox() < cursor.getMinBox())) return;
		tmppl = ite->d.planeIntersection(brother);
		if (tmppl == 3u) return;
		if (ite->d.isTrivial()) ite->d.TrivialMerge(brother, toexp); // merging opposites!
		else if (toexp.trivialIsCoherrentWith(tmppl)) toexp = ite->d;
		else return;
		if (_debug_) printf("recurmerge what? %X\n", toexp());
		ncursor.toParent();
		return writeHyperPosition_routine(ncursor, toexp);
	}else{
		Tuple<S,D> damin;
		N tmptmp;
		bool change = false;
		if (_debug_) printf("end check\n");
		if (ite.isValid()) while(ite->d.isTrivial()){
			if (ite->d.trivialIsCoherrentWith(toexp.planeIntersection(cursor = HyperCursor<S,D,L>(ite->k)))){
				tmptmp = ite->d;
				++ite;
				if (ite.isValid()) damin = HyperCursor<S,D,L>(ite->k).getMin();
				cursor.findNext();
				if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
				change = false;
				while(cursor.mag != 0xFFFF) {
					if (!tmptmp.trivialIsCoherrentWith(toexp.planeIntersection(cursor))){
						if (ite.isValid()) --ite;
						else ite.findLast();
						ite.ordering_preserving_reference().k = cursor.hyperpos;
						if (_debug_) {printf("could not destroy %X at\n", ite->d()); cursor.show();}
						change = true;
						break;
					}
					cursor.findNext();
					if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
				}
				if (change) break;
				if (ite.isValid()) --ite;
				else ite.findLast();
				if (_debug_) {printf("Destroying %X\n", ite->d());	HyperCursor<S,D,L>(ite->k).show();}
				space.remove_pp(ite);
				change = true;
				if (!ite.isValid()) break;
				if (ite->d == toexp) {
					space.remove_pp(ite);
					if (!ite.isValid()) break;
				}
			} else break;
		}
		if (change){
			if ((hasNext = ite.isValid())) {damin = HyperCursor<S,D,L>(ite->k).getMin();	--ite;}
			else ite.findLast();
			cursor = HyperCursor<S,D,L>(ite->k);
			cursor.toLargeMinBoxContainer();
			if (hasNext) cursor.toMinBoxExclusiveContainer(damin);
			ite.ordering_preserving_reference().k = cursor.hyperpos;
			return HyperPositionInflate_routine(ite);
		}
		if (_debug_) printf("nothing!\n");
	}
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::HyperPositionInflateMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, bool debug){
	HyperCursor<S,D,L> cursor = ite->k;
	N toexp = ite->d;
	if (!ite->d.canSimplify()) return;
	N simplified[2];
	if (debug) {printf("Inflating: ");toexp.show();	cursor.show();}
    bool hasNext;
	Tuple<S,D> minOfNext;

	/*
	++ite;
	cursor.toLargeMinBoxContainer();



	if (hasNext) {
		if (debug) {printf("NextLimit %X: ", ite->d()); HyperCursor<S,D,L>(ite->k).show();}
		minOfNext = HyperCursor<S,D,L>(ite->k).getMin();
		cursor.toMinBoxExclusiveContainer(minOfNext);
		--ite;
	}else ite.findLast();
	--ite;

	do{
		if (cursor.mag == 0xFFFF) printf("should be impossible... got intersecting nodes! fdadf\n");
		HyperCursor<S,D,L> brother = cursor.getBrother();

		if (debug) printf("inner! continue %c\n", (brother <= cursor)  ? 'Y' : 'N');
		if (brother > cursor) break;

		toexp.simplifyAt(simplified[0],brother);


		if (ite.isValid()){
		//	if (toexp == ite->d){ //
		//		cursor = ite->k;
		//		++ite;
		//		space.remove_mm(ite);
		//	}else{

			HyperCursor<S,D,L> prev(ite->k);
			if (debug) {printf("PrevFound: ");  prev.show();}
			//ite->k.show();
			mag = brother.commonContainer_mag(prev);
			if (debug) {
				printf("Cur: ");cursor.hyperpos.show();
				printf("mag for cur... %X\n", cursor.mag);
				printf("Cur Bro: ");  brother.show();
				printf("PE %X %X %X\n", mag, brother.mag, prev.mag);
			}
			if (mag > brother.mag){
				if (debug) printf("HH\n");
				if (!simplified[0].canMerge()) break;
				ite->d.simplifyAt(simplified[1],brother);
				if (simplified[0] == simplified[1]){
					if (debug) printf("aggree... %c\n", (toexp.canMerge())? 'Y' : 'N');
					if (toexp.canMerge()){
                        ite->d.simplifyAt(simplified[1], cursor);
						if (simplified[1] == simplified[0]){
							cursor.toParent();
							return writeHyperPositionMK2_routine(cursor, ite->d,debug);
						}
					}
				} else if ((simplified[1].canMerge())&&(toexp.canMerge())){
					//NODES tmptmp; tmptmp.toTrivial(tmppl2);
					//tmptmp.TrivialMerge(cursor, toexp);
					cursor.toParent();
					if (debug) printf("did an trivial merge path3! %X\n", toexp());
					return writeHyperPositionMK2_routine(cursor, toexp,debug);
				} else {if (debug) printf("disaggree...\n"); break;}
				cursor.toParent();
				if (debug) printf("PFF\n");
			}else if (brother.hyperpos == ite->k){
				if (debug) printf("PFCdd\n");
				if (ite->d.canMerge()){
					if (ite->d == simplified[0]){
					}else if (ite->d.mergeInto(toexp,cursor)) {
						++ite; ite.ordering_preserving_reference().k = cursor.hyperpos;
						if (debug) printf("did an trivial merge! %X\n", toexp());
						cursor.toParent();
						return writeHyperPositionMK2_routine(cursor, toexp,debug);
					}else break;
					cursor.toParent();
					space.remove_mm(ite);
					if ((ite.isValid())&&(ite->d == toexp)) {
						// previous matches! use it
						cursor = ite->k;
						++ite;
						if (debug) printf("auto remove previous, matches!\n");
						space.remove_mm(ite);
						--ite;
					}
				}else{
					// previous is brother, hence ite->k.planeIntersection(brother) == 3
					 // DONE!
					break;
				}
			}else break;
		//	}
		}else if (toexp == simplified[0]) break;
		else{
			// previous is non-existant (trivial 0)
			// *TODO*
			HyperCursor<S,D,L> prev();
			printf("PFWTF\n");
			break;
		}
		//cursor.show();
		cursor.toLargeMinBoxContainer();
		if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
		if (cursor.mag == 0xFFFF){
			printf("got cornered!\n");
			ite->k.show();
			ite->d.show();
			if (hasNext) {printf("nextis:"); ExOp::show(minOfNext);}
			else printf("no next\n");
			LFH_exit(1);
		}
	}while(true);
	if (ite.isValid()) ++ite;
	else ite.findFirst();*/
	//printf("to %X\n", cursor.mag);
	//ite.ordering_preserving_reference().k = cursor.hyperpos;
	++ite;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	while(ite.isValid()){
        if (debug){
            printf("Trying to eat: "); ite->d.show();
            printf("Sitting at:"); HyperCursor<S,D,L>(ite->k).show();
        }
        toexp.simplifyAt(simplified[0],ite->k);
        if (simplified[0] != ite->d) break;
        cursor = ite->k;
        ++ite;
        if ((hasNext = ite.isValid())) {minOfNext = HyperCursor<S,D,L>(ite->k).getMin(); --ite;}
        else ite.findLast();
        do{
            cursor.findNext();
            if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
            if (cursor.mag == 0xFFFF){
                // next is vanishing!
                space.remove_pp(ite);
                if ((ite.isValid())&&(ite->d == toexp)) space.remove_pp(ite);
                hasNext = false;
                break;
            }
            toexp.simplifyAt(simplified[0],cursor);
            ite->d.simplifyAt(simplified[1],cursor);
            if (simplified[0] != simplified[1]){
                ite.ordering_preserving_reference().k = cursor.hyperpos;
                if (ite->d != simplified[1]){
                    // needs to insert this, and fix after, ya, real fun
                    toexp = ite->d;
                    ite.ordering_preserving_reference().d = simplified[1];
                    toins.d = simplified[1];
                    HyperCursor<S,D,L> cursorB = ite->k;
                    if (debug) {printf("downgraded next, fixin'\n");cursor.show();}
                    cursor.findNext();
                    if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
                    while(cursor.mag != 0xFFFF){
                        if (!toexp.simplifyAt(simplified[0],cursor)) {
                            toins.d = toexp; toins.k = cursor.hyperpos; space.insert(toins);
                            break;
                        }
                        toins.d.simplifyAt(simplified[1],cursor);
                        if (simplified[1] != simplified[0]){
                            toins.d = simplified[0];
                            toins.k = cursor.hyperpos;
                            space.insert(toins);
                        }
                        cursor.findNext();
                        if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
                    }
                    ite.findLE(cursorB);
                }
                hasNext = true;
                break;
            }
        }while(true);
        if (hasNext) break;
        /*
        if (toexp.canMerge()){


            HyperCursor<S,D,L> ncursor = HyperCursor<S,D,L>(ite->k);
            //if (ite->d.isTrivial()){
            //	if (ite->d.trivialIsCoherrentWith(toexp.planeIntersection(ncursor))) return writeHyperPosition_routine(ncursor, toexp);
            //}

            if (debug) {printf("recurmerge!\n"); ncursor.show();}
            HyperCursor<S,D,L> brother = ncursor.getBrother();
            if ((brother > ncursor)||(brother.getMinBox() < cursor.getMinBox())) return;
            if (!ite->d.simplifyAt(simplified[0], brother)) return;
            if (!ite->d.mergeInto(toexp, brother)) { // merging opposites!
                if (toexp == simplified[0]) toexp = ite->d;
                else return;
            }
            ncursor.toParent();
            return writeHyperPositionMK2_routine(ncursor, toexp,debug);
        }else{
            Tuple<STORETYPE,DIMS> damin;
            NODES tmptmp;
            bool change = false;
            if (debug) printf("end check\n");
            if (ite.isValid()) while(ite->d.canMerge()){
                toexp.simplifyAt(simplified[1], cursor = HyperCursor<S,D,L>(ite->k) );
                if (ite->d == simplified[1]){
                    tmptmp = ite->d;
                    ++ite;
                    if (ite.isValid()) damin = HyperCursor<S,D,L>(ite->k).getMin();
                    cursor.findNext();
                    if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
                    change = false;
                    while(cursor.mag != 0xFFFF) {
                        toexp.simplifyAt(simplified[1], cursor);
                        if (tmptmp != simplified[1]){
                            if (ite.isValid()) --ite;
                            else ite.findLast();
                            ite.ordering_preserving_reference().k = cursor.hyperpos;
                            if (debug) {printf("could not destroy %X at\n", ite->d()); cursor.show();}
                            change = true;
                            break;
                        }
                        cursor.findNext();
                        if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
                    }
                    if (change) break;
                    if (ite.isValid()) --ite;
                    else ite.findLast();
                    if (debug) {printf("Destroying %X\n", ite->d());	HyperCursor<S,D,L>(ite->k).show();}
                    space.remove_pp(ite);
                    change = true;
                    if (!ite.isValid()) break;
                    if (ite->d == toexp) {
                        space.remove_pp(ite);
                        if (!ite.isValid()) break;
                    }
                } else break;
            }
            if (change){
                if ((hasNext = ite.isValid()) == true) {damin = HyperCursor<S,D,L>(ite->k).getMin();	--ite;}
                else ite.findLast();
                cursor = HyperCursor<S,D,L>(ite->k);
                cursor.toLargeMinBoxContainer();
                if (hasNext) cursor.toMinBoxExclusiveContainer(damin);
                ite.ordering_preserving_reference().k = cursor.hyperpos;
                return HyperPositionInflateMK2_routine(ite,debug);
            }
            if (debug) printf("nothing!\n");
        }*/
	}

	//update magnitude!]
    if ((hasNext = ite.isValid())){
        minOfNext = HyperCursor<S,D,L>(ite->k).getMin();
        if (debug) HyperCursor<S,D,L>(ite->k).show();
        --ite;
    }else ite.findLast();

    cursor = ite->k;
    if (debug) {ExOp::show(minOfNext);printf("inflated from mag %X to", cursor.mag);}
    cursor.toLargeMinBoxContainer();

	if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
	while(cursor.isOlderBrother()){
	    HyperCursor<S,D,L> cursorC = cursor.mkBrother();
        if (!this->getAt_update(toins.d, cursorC)) break;
        ite->d.simplifyAt(simplified[0],cursorC);
        if (simplified[0] != toins.d) break;
        cursor.toParent();cursor.toLargeMinBoxContainer();
        if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
	}
	if (debug) cursor.show();
	//

	--ite;
	if (!ite.isValid()){
	    ite.findFirst();
        ite.ordering_preserving_reference().k = cursor.hyperpos;
        return;

	}

    if (cursor.getMinBox() > ite->k){
        ++ite;
        ite.ordering_preserving_reference().k = cursor.hyperpos;
        return;
	}
	++ite;
    toins.d = ite->d;
    toins.k = cursor.hyperpos;
    space.removeRange(cursor.getMinBox(), cursor.getMaxBox());
    space.insert(toins);

}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::HyperPositionInflateMK2_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const M& merger, bool debug){
	HyperCursor<S,D,L> cursor = ite->k;
	N toexp = ite->d;
	if (!ite->d.canSimplify()) return;
	N simplified[2];
	if (debug) {printf("Inflating: ");toexp.show();	cursor.show();}
    bool hasNext;
	Tuple<S,D> minOfNext;
	++ite;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	while(ite.isValid()){
        if (debug){
            printf("Trying to eat: "); merger.showThat(ite->d);
            printf("Sitting at:"); HyperCursor<S,D,L>(ite->k).show();
        }
        merger.simplifyAt(toexp, simplified[0],ite->k, debug);
        if (debug){  printf("simplified to: ");merger.showThat(simplified[0]);}
        if (simplified[0] != ite->d) break;
        cursor = ite->k;
        ++ite;
        if ((hasNext = ite.isValid())) {minOfNext = HyperCursor<S,D,L>(ite->k).getMin(); --ite;}
        else ite.findLast();
        do{
            cursor.findNext();
            if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
            if (cursor.mag == 0xFFFF){
                // next is vanishing!
                space.remove_pp(ite);
                if ((ite.isValid())&&(ite->d == toexp)) space.remove_pp(ite);
                hasNext = false;
                break;
            }
            merger.simplifyAt(toexp,simplified[0],cursor, debug);
            merger.simplifyAt(ite->d,simplified[1],cursor, debug);
            if (simplified[0] != simplified[1]){
                ite.ordering_preserving_reference().k = cursor.hyperpos;
                if (ite->d != simplified[1]){
                    // needs to insert this, and fix after, ya, real fun
                    toexp = ite->d;
                    ite.ordering_preserving_reference().d = simplified[1];
                    toins.d = simplified[1];
                    HyperCursor<S,D,L> cursorB = ite->k;
                    if (debug) {printf("downgraded next, fixin'\n");cursor.show();}
                    cursor.findNext();
                    if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
                    while(cursor.mag != 0xFFFF){
                        if (!merger.simplifyAt(toexp,simplified[0],cursor, debug)) {
                            toins.d = toexp; toins.k = cursor.hyperpos; space.insert(toins);
                            break;
                        }
                        merger.simplifyAt(toins.d,simplified[1],cursor, debug);
                        if (simplified[1] != simplified[0]){
                            toins.d = simplified[0]; toins.k = cursor.hyperpos; space.insert(toins);
                        }
                        cursor.findNext();
                        if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
                    }
                    ite.findLE(cursorB);
                }
                hasNext = true;
                break;
            }
        }while(true);
        if (hasNext) break;
	}

    if ((hasNext = ite.isValid())){
        minOfNext = HyperCursor<S,D,L>(ite->k).getMin();
        if (debug) HyperCursor<S,D,L>(ite->k).show();
        --ite;
    }else ite.findLast();

    cursor = ite->k;
    if (debug) {ExOp::show(minOfNext);printf("inflated from mag %X to", cursor.mag);}
    cursor.toLargeMinBoxContainer();

	if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
	if (debug){
        printf("oldermerge!\n");
        cursor.show();
        toins.d.show();
	}

	while(cursor.isOlderBrother()){
	    HyperCursor<S,D,L> cursorC = cursor.mkBrother();
        if (!this->getAt_update(toins.d, cursorC,merger)) break;
        merger.simplifyAt(ite->d,simplified[0],cursorC,debug);
        if (simplified[0] != toins.d) break;
        cursor.toParent();cursor.toLargeMinBoxContainer();
        if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
	}
	if (debug) cursor.show();

	--ite;
	if (!ite.isValid()){
	    ite.findFirst();
	    if (debug){ite.ordering_preserving_reference().d.show(); printf("changed to:"); cursor.show(); }
        ite.ordering_preserving_reference().k = cursor.hyperpos;
        return;

	}

    if (cursor.getMinBox() > ite->k){
        ++ite;
        if (debug) {ite.ordering_preserving_reference().d.show(); printf("chanGed to:"); cursor.show(); }
        ite.ordering_preserving_reference().k = cursor.hyperpos;
        return;
	}
	++ite;
    toins.d = ite->d;
    toins.k = cursor.hyperpos;
    if (debug) {ite.ordering_preserving_reference().d.show(); printf("would be inserted at:"); cursor.show(); }

    space.removeRange(cursor.getMinBox(), cursor.getMaxBox());
    space.insert(toins);
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::writeHyperPosition_routine(const HyperCursor<S,D,L> &where, const N& what){
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(this->space);
	HyperCursor<S,D,L> cursor;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	N prevprop;
	N nextprop;
	cursor = where.getMaxBox();
	//printf("Query:");
	//cursor.show();

	bool _debug_ = false;
	if (_debug_) {printf("Actual insertion: %X  mag %X\n", what(), where.mag);where.show();}

	if (this->writePrevNext_routine(where,what,ite,prevprop,nextprop)) return;

	if (_debug_) printf("Prev %X, Next %X, ToIns %X (cursize %i)\n", prevprop(), nextprop(), what(), space.getSize());

	if (nextprop != what){
		cursor = where;
		if (cursor.findNext()){
			if (_debug_) printf("gota\n");
			toins.d = nextprop;
			Tuple<S,D> damin;
			if (ite.isValid()){
				damin = HyperCursor<S,D,L>(ite->k).getMin();
			//	damin.show();
			//	ite->k.show();cursor.show();
				//if (cursor.mag == cursor.commonContainer_mag(tmp))
				cursor.toMinBoxExclusiveContainer(damin);
			}
			if (cursor.mag == 0xFFFF){ // no next, merge is possible!
				if (_debug_) printf("NO NEXT!!!\n");
				if (ite.isValid()){
					if (ite->d == what) { // merging to next!
						if (prevprop == what){
							//printf("remnext\n");
							space.remove_mm(ite);
                            if (!ite.isValid()){
                                if (_debug_) printf("end that end?\n");
                                return;
                            }
						}else{
							ite.ordering_preserving_reference().k = where.hyperpos;
						}
						if (_debug_) printf("INF call A\n");
						//this->checkIntegrity();
					}else if (prevprop != what){
						toins.d = what;
						toins.k = where.hyperpos;
						space.insert(toins);
						ite.findLE(where);
						if (_debug_) printf("INF call hAh\n");
					}else{
						--ite;
						//printf("INF call B\n");
						if (!ite.isValid()) return;
					}
					HyperPositionInflate_routine(ite);
				}else if (prevprop != what){
					if (_debug_) printf("super espacedd!\n");
					toins.d = what;
					toins.k = where.hyperpos;
					space.insert(toins);
				} else{ if (_debug_) printf("miiiiised!\n");}
				return;
			}else{
				if (_debug_) printf("gotD %X\n", nextprop());
				if (nextprop.isTrivial()) {
					toins.k = cursor.hyperpos;
					space.insert(toins); //printf("nxt in: "); toins.k.show();//if (_debug_)
                    //    printf("got %X inserted\n", toins.d());
				}else{
					bool startin = true;

					do{
						if (_debug_) printf("Duper fun! %c\n", startin ? 'Y' : 'N');
						//space.show();
						unsigned int tmppl = nextprop.planeIntersection(cursor);
						if (tmppl == 3u) {
							//printf("%X gave int3\n",  nextprop());
							//printf("insert %X at:\n", toins.d); cursor.show();
							//printf("insert %X yes\n", toins.d);
							toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); startin = false; break;
						}
						if (tmppl != what.planeIntersection(cursor)){
							//printf("duper disagree! %i and %i\n", tmppl,  what.planeIntersection(cursor));
							toins.d.toTrivial(tmppl);
							toins.k = cursor.hyperpos;
							startin = false;
							//printf("Duper fun! %c\n", startin ? 'Y' : 'N');
							if (!startin) space.insert(toins);
							cursor.findNext();
							if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							do{
								if (cursor.mag == 0xFFFF) {
										ite.findGE(toins); // was just inserted, should be there
										++ite;
										if ((ite.isValid())&&(tmppl == ite->d.planeIntersection())) space.remove(ite);
									break;
								}
								unsigned int tmppl2 = nextprop.planeIntersection(cursor);
								if (tmppl2 == 3u) {toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); break;}
								if (tmppl2 != tmppl){
									tmppl = tmppl2;
									toins.d.toTrivial(tmppl);
									toins.k = cursor.hyperpos;
									//printf("Duper fun! N=%c\n", startin ? 'Y' : 'N');
									space.insert(toins);
								}
								cursor.findNext();
								if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							}while(true);
							break;
						}
						cursor.findNext();
						if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
					} while(cursor.mag != 0xFFFF);
					if (_debug_) {printf("after duper\n");space.show();}
					if (startin == true){ // ended up deleting nodes, not inserting
						if (_debug_) {printf("got hhheehheh %X\n", ite->d());	ExOp::show(ite->k);}

						if (ite.isValid()){
							if (ite->d == what) { // merging to next!
								if (prevprop == what){
									if (_debug_) printf("remnext... reallyss?\n");
									space.remove_mm(ite);
									if (!ite.isValid()){
                                        if (_debug_) printf("end that end?\n");
                                        return;
									}
								}else if ((what.isTrivial())&&(what.trivialIsCoherrentWith(prevprop.planeIntersection(where.hyperpos)))){
									if (_debug_) printf("simply! not eqq!\n");
									return writeHyperPosition_routine(ite->k, prevprop);
								}else{
									ite.ordering_preserving_reference().k = where.hyperpos;
								}
								if (_debug_) printf("INF call A_hehe! does inflate!\n");
								//this->checkIntegrity();
							}else if (prevprop != what){
								if (_debug_) printf("not eqq!\n");
								toins.d = what;
								toins.k = where.hyperpos; if (_debug_) printf("INF call Bsds\n");
								space.insert(toins);
								ite.findLE(where);
							}else{
								--ite;
								if (_debug_) printf("INF call B\n");
								if (!ite.isValid()) return;
							}
							HyperPositionInflate_routine(ite);
						}else if (prevprop != what){
							toins.d = what;
							toins.k = where.hyperpos;if (_debug_) printf("INF call Bh\n");
							space.insert(toins);
						}else{ if (_debug_) printf("got nowhre!\n");}
						return;
					}
				}
			}

		}
	}else{ // inflated node might be needed, but not done yet
		if (_debug_) printf("gotb\n");
		if (prevprop != what){
			toins.d = what;
			toins.k = where.hyperpos;
			space.insert(toins);
		}
		ite.findLE(where);
		//printf("INF call C\n");

		if (ite.isValid()) HyperPositionInflate_routine(ite);
		return;
	}
	if (_debug_) printf("gotc\n");
	if (prevprop != what){
		//printf("add last %X != %X\n", what.anbox, prevprop.anbox );
//		if ((what.isTrivial())&&(what.trivialIsCoherrentWith(prevprop.planeIntersection(where)))){

//		}else{
		toins.d = what;
		toins.k = where.hyperpos;
		space.insert(toins);
		if (_debug_) {printf("again insert %X\n", toins.d());		space.show();}
//		}
		ite.findLE(toins.k);
	}else{
		if (_debug_) printf("would upgrade...\n");
		// prev node might need to be enlarged!
		ite.findLE(where);
		if (!ite.isValid()) return;
	}
	return HyperPositionInflate_routine(ite);
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::writeHyperPositionMK2_routine(const HyperCursor<S,D,L> &where, N& what, bool debug){
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(this->space);
	HyperCursor<S,D,L> cursor;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	N prevprop;
	N nextprop;
	N simplified[2];
	cursor = where.getMaxBox();
    bool _debug_ = false;
	if (_debug_) {printf("Actual insertion: %X  mag %X\n", what(), where.mag);where.show();}
	if (this->writePrevNextMK2_routine(where,what,ite,prevprop,nextprop,debug)) return;
	if (debug){
        printf("This: "); what.show();
        printf("Prev: "); prevprop.show();
        printf("Next: "); nextprop.show();
	}
	if (_debug_) printf("Prev %X, Next %X, ToIns %X (cursize %i)\n", prevprop(), nextprop(), what(), space.getSize());

	// if agrees with previous, assume it *will* like previous
	prevprop.simplifyAt(simplified[0], where);
	if (simplified[0] == what) what = prevprop;

	if (nextprop != what){
		cursor = where;
		if (cursor.findNext()){
			if (debug) printf("got a\n");
			toins.d = nextprop;
			Tuple<S,D> damin;
			if (ite.isValid()){
				damin = HyperCursor<S,D,L>(ite->k).getMin();
			//	damin.show();
			//	ite->k.show();cursor.show();
				//if (cursor.mag == cursor.commonContainer_mag(tmp))
				cursor.toMinBoxExclusiveContainer(damin);
			}
			if (cursor.mag == 0xFFFF){ // no space for next, merge with actual next is possible!
				if (debug) printf("NO NEXT!!!\n");
				if (ite.isValid()){
					if (ite->d == what) { // merging to next if equal, inflating should do the rest!
						space.remove(ite);
						/*if (prevprop == what){
							//printf("remnext\n");
							space.remove_mm(ite);
						}else{
							ite.ordering_preserving_reference().k = where.hyperpos;
						}
						if (debug) printf("INF call A\n");*/
						//this->checkIntegrity();
					}/*else if (prevprop != what){
						toins.d = what;
						toins.k = where.hyperpos;
						space.insert(toins);
						ite.findLE(where);
						if (debug) printf("INF call hAh\n");
					}else{
						--ite;
						//printf("INF call B\n");
						if (!ite.isValid()) return;
					}*/
					//HyperPositionInflateMK2_routine(ite,debug);
				}/*else if (prevprop != what){
					if (debug) printf("super espacedd!\n");
					toins.d = what;
					toins.k = where.hyperpos;
					space.insert(toins);
				} else{ if (debug) printf("miiiiised!\n");}
				return;*/
			}else{
				if (debug) printf("gotD %X\n", nextprop());
				toins.d = what;
				if (nextprop.canSimplify()) {
					bool startin = true;
					do{
						if (debug) printf("Duper fun! %c\n", startin ? 'Y' : 'N');
						//space.show();
						if (!nextprop.simplifyAt(simplified[0],cursor)) {
							//printf("%X gave int3\n",  nextprop());
							//printf("insert %X at:\n", toins.d); cursor.show();
							//printf("insert %X yes\n", toins.d);
							toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); startin = false;
							break;
						}
						toins.d.simplifyAt(simplified[1],cursor);
						if (simplified[1] != simplified[0]){
                            if (debug){
                                printf("duper disagree!\n");
                                nextprop.show();
                                simplified[0].show();
                                toins.d.show();
                                simplified[1].show();
                                cursor.show();
							}
							toins.d = simplified[0];
							toins.k = cursor.hyperpos;
							startin = false;
							//printf("Duper fun! %c\n", startin ? 'Y' : 'N');
							space.insert(toins);
							cursor.findNext();
							if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							if (debug){
                                printf("fixing next: first with:");  toins.d.show();
                                damin.show();
							}
							do{
								if (cursor.mag == 0xFFFF) {
                                    ite.findGE(toins); // was just inserted, should be there
                                    ++ite;
                                    if (debug){printf("danext is..."); ite->d.show();}
                                    if ((ite.isValid())&&(ite->d == toins.d)) space.remove(ite);
									break;
								}

								if (!nextprop.simplifyAt(simplified[0],cursor)) {
                                    toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins);
                                    break;
                                }

                                toins.d.simplifyAt(simplified[1],cursor);

								if (simplified[1] != simplified[0]){
									toins.d = simplified[0];
									toins.k = cursor.hyperpos;
									if (debug){printf("inserted an:"); toins.d.show();}
									space.insert(toins);
								}else if (debug) {simplified[1].show();cursor.show();}
								cursor.findNext();
								if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							}while(true);
							break;
						}
						cursor.findNext();
						if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
					} while(cursor.mag != 0xFFFF);
					if (debug) {printf("after duper\n");}
					if (startin == true){ // ended up deleting nodes, not inserting
						if (ite.isValid()){
                            if (debug) {printf("got hhheehheh %X\n", ite->d());	ExOp::show(ite->k);}
							if (ite->d == what) { // merging to next!
								if (prevprop == what){
									if (debug) printf("remnext... realddly?\n");
									space.remove_mm(ite);
									if (!ite.isValid()){
                                        if (debug) printf("end that end?\n");
                                        return;
									}
								}else{
								    prevprop.simplifyAt(simplified[1], where);
                                    if ((what.canMerge())&&(what == simplified[1])){
                                    if (debug) printf("simply! not eqq!\n");
                                        return writeHyperPositionMK2_routine(ite->k, prevprop,debug);
                                    }else{
                                        ite.ordering_preserving_reference().k = where.hyperpos;
                                    }
								}
								if (debug) printf("INF call A_hehe! does inflate!\n");
								//this->checkIntegrity();
							}else if (prevprop != what){
								if (debug) printf("not eqq!\n");
								toins.d = what;
								toins.k = where.hyperpos; if (_debug_) printf("INF call Bsds\n");
								space.insert(toins);
								ite.findLE(where);
							}else{
								--ite;
								if (debug) printf("INF call B\n");
								if (!ite.isValid()) return;
							}
							HyperPositionInflateMK2_routine(ite,debug);
						}else if (prevprop != what){
							toins.d = what;
							toins.k = where.hyperpos;if (debug) printf("INF call Bh\n");
							space.insert(toins);
						}else{ if (debug) printf("got nowhre!\n");}

						return;
                    }
                }else{
                    toins.d = nextprop;
					toins.k = cursor.hyperpos;
					space.insert(toins);
					if (debug){
                        printf("Next is trivial:"); toins.d.show();
                        printf("Inserted at:"); cursor.show();
					}
		}   }	}
	}/*else{ // inflated node might be needed, but not done yet
		if (debug) printf("gotb\n");
		if (prevprop != what){
			toins.d = what;
			toins.k = where.hyperpos;
			space.insert(toins);
		}
		ite.findLE(where);
		//printf("INF call C\n");

		if (ite.isValid()) HyperPositionInflateMK2_routine(ite,debug);
		return;
	}*/
	if (debug) {printf("gotc\n");}

	if (prevprop != what){
		toins.d = what;
		toins.k = where.hyperpos;
		space.insert(toins);
		if (debug) {printf("again insert %X\n", toins.d());}
		ite.findLE(toins.k);
	}else{
		if (debug) printf("would upgrade...\n");
		ite.findLE(where);
		if (!ite.isValid()) return;
	}
	return HyperPositionInflateMK2_routine(ite,debug);
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::writeHyperPositionMK2_routine(const HyperCursor<S,D,L> &where, N& what, const M& merger, bool debug){
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(this->space);
	HyperCursor<S,D,L> cursor;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	N prevprop;
	N nextprop;
	N simplified[2];
	cursor = where.getMaxBox();
    bool _debug_ = false;
	if (_debug_) {printf("Actual insertion: %X  mag %X\n", what(), where.mag);where.show();}
	if (this->writePrevNextMK2_routine(where,what,ite,prevprop,nextprop,merger,debug)) return;
	if (debug){
        printf("This: "); merger.showThat(what);
        printf("Prev: "); merger.showThat(prevprop);
        printf("Next: "); merger.showThat(nextprop);
	}
	if (_debug_) printf("Prev %X, Next %X, ToIns %X (cursize %i)\n", prevprop(), nextprop(), what(), space.getSize());

	// if agrees with previous, assume it *will* like previous
	merger.simplifyAt(prevprop,simplified[0], where, debug);
	if (simplified[0] == what) what = prevprop;

	if (nextprop != what){
		cursor = where;
		if (cursor.findNext()){
			if (debug) {
                printf("got a\n");
                printf("Inserted at:"); where.show();
                cursor = where; cursor.findNext();
                printf("nextwas...?"); cursor.show();
            }

			toins.d = nextprop;
			Tuple<S,D> damin;
			if (ite.isValid()){
				damin = HyperCursor<S,D,L>(ite->k).getMin();
			//	damin.show();
			//	ite->k.show();cursor.show();
				//if (cursor.mag == cursor.commonContainer_mag(tmp))
				cursor.toMinBoxExclusiveContainer(damin);
			} // else no next...
			if (cursor.mag == 0xFFFF){ // no space for next, merge with actual next is possible!
				if (debug) printf("NO NEXT!!!\n");
				if (ite.isValid()){
					if (ite->d == what) { // merging to next if equal, inflating should do the rest!
						space.remove(ite);
						/*if (prevprop == what){
							//printf("remnext\n");
							space.remove_mm(ite);
						}else{
							ite.ordering_preserving_reference().k = where.hyperpos;
						}
						if (debug) printf("INF call A\n");*/
						//this->checkIntegrity();
					}/*else if (prevprop != what){
						toins.d = what;
						toins.k = where.hyperpos;
						space.insert(toins);
						ite.findLE(where);
						if (debug) printf("INF call hAh\n");
					}else{
						--ite;
						//printf("INF call B\n");
						if (!ite.isValid()) return;
					}*/
					//HyperPositionInflateMK2_routine(ite,debug);
				}/*else if (prevprop != what){
					if (debug) printf("super espacedd!\n");
					toins.d = what;
					toins.k = where.hyperpos;
					space.insert(toins);
				} else{ if (debug) printf("miiiiised!\n");}
				return;*/
			}else{
				if (debug) printf("gotD %X\n", nextprop());
				toins.d = what;
				if (nextprop.canSimplify()) {
					bool startin = true;
					do{
						if (debug) printf("Duper fun! %c\n", startin ? 'Y' : 'N');
						//space.show();
						if (!merger.simplifyAt(nextprop,simplified[0],cursor)) {
							//printf("%X gave int3\n",  nextprop());
							//printf("insert %X at:\n", toins.d); cursor.show();
							//printf("insert %X yes\n", toins.d);
							toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); startin = false;
							break;
						}
						merger.simplifyAt(toins.d,simplified[1],cursor);
						if (simplified[1] != simplified[0]){
                            if (debug){
                                printf("duper disagree!\n");
                                nextprop.show();
                                simplified[0].show();
                                toins.d.show();
                                simplified[1].show();
                                cursor.show();
							}
							toins.d = simplified[0];
							toins.k = cursor.hyperpos;
							startin = false;
							//printf("Duper fun! %c\n", startin ? 'Y' : 'N');
							space.insert(toins);
							cursor.findNext();
							if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							if (debug){
                                printf("fixing next: first with:");  toins.d.show();
                                damin.show();
							}
							do{
								if (cursor.mag == 0xFFFF) {
                                    ite.findGE(toins); // was just inserted, should be there
                                    ++ite;
                                    if (debug){  printf("danext is..."); if (ite.isValid()) ite->d.show(); else printf(" invalid\n");}
                                    if ((ite.isValid())&&(ite->d == toins.d)) space.remove(ite);
									break;
								}

								if (!merger.simplifyAt(nextprop,simplified[0],cursor)) {
                                    toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins);
                                    break;
                                }

                                merger.simplifyAt(toins.d,simplified[1],cursor);

								if (simplified[1] != simplified[0]){
									toins.d = simplified[0];
									toins.k = cursor.hyperpos;
									if (debug){printf("inserted an:"); toins.d.show();}
									space.insert(toins);
								}else if (debug) {simplified[1].show();cursor.show();}
								cursor.findNext();
								if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							}while(true);
							break;
						}
						cursor.findNext();
						if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
					} while(cursor.mag != 0xFFFF);
					if (debug) {printf("after duper\n");}
					if (startin == true){ // ended up deleting nodes, not inserting
						if (ite.isValid()){
                            if (debug) {printf("got hhheehheh %X\n", ite->d());	ExOp::show(ite->k);}
							if (ite->d == what) { // merging to next!
								if (prevprop == what){
									if (debug) printf("remnext... reallyee?\n");
									space.remove_mm(ite);
                                    if (!ite.isValid()){
                                        if (debug) printf("end that end?\n");
                                        return;
									}
								}else{
								    merger.simplifyAt(prevprop,simplified[1], where);
                                    if ((what.canMerge())&&(what == simplified[1])){
                                    if (debug) printf("simply! not eqq!\n");
                                        return writeHyperPositionMK2_routine(ite->k, prevprop,merger,debug);
                                    }else{
                                        ite.ordering_preserving_reference().k = where.hyperpos;
                                    }
								}
								if (debug) printf("INF call A_hehe! does inflate!\n");
								//this->checkIntegrity();
							}else if (prevprop != what){
								if (debug) printf("not eqq!\n");
								toins.d = what;
								toins.k = where.hyperpos; if (_debug_) printf("INF call Bsds\n");
								space.insert(toins);
								ite.findLE(where);
							}else{
								--ite;
								if (debug) printf("INF call B\n");
								if (!ite.isValid()) return;
							}
							HyperPositionInflateMK2_routine(ite,merger,debug);
						}else if (prevprop != what){
							toins.d = what;
							toins.k = where.hyperpos;if (debug) printf("INF call Bh\n");
							space.insert(toins);
						}else{ if (debug) printf("got nowhre!\n");}

						return;
                    }
                }else{
                    toins.d = nextprop;
					toins.k = cursor.hyperpos;
					space.insert(toins);
					if (debug){
                        printf("Next is trivial:"); toins.d.show();
                        printf("Inserted at:"); cursor.show();
					}
		}   }	}
	}/*else{ // inflated node might be needed, but not done yet
		if (debug) printf("gotb\n");
		if (prevprop != what){
			toins.d = what;
			toins.k = where.hyperpos;
			space.insert(toins);
		}
		ite.findLE(where);
		//printf("INF call C\n");

		if (ite.isValid()) HyperPositionInflateMK2_routine(ite,debug);
		return;
	}*/
	if (debug) {printf("gotc\n");}

	if (prevprop != what){
		toins.d = what;
		toins.k = where.hyperpos;
		space.insert(toins);
		if (debug) {printf("again insert\n"); toins.k.show();toins.d.show();}
		ite.findLE(toins.k);
	}else{
		if (debug) printf("would upgrade...\n");
		ite.findLE(where);
		if (!ite.isValid()) return;
	}
	return HyperPositionInflateMK2_routine(ite,merger,debug);
}
LFHTEMP HierarchicalMergable<N,S,D,L>& HierarchicalMergable<N,S,D,L>::toZero(){
	space.toMemfree();
	return *this;
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::insert(const HyperCursor<S,D,L> &pos, const N& n_node){
	N query;
	HyperCursor<S,D,L> tpos;
	if (!n_node.isTrivial()) {query = n_node; tpos = pos;}
	else{ // trying to insert trivial block

		//this->clearHyperPosition_routine(pos);
		//printf("hello trivial being!\n");
		bool isPure = this->getAt(query,pos);
		if (isPure){
			if (n_node.trivialIsCoherrentWith(query.planeIntersection(pos))) {
			//	printf("super coherent!\n");
				return;
			}
		}

		tpos = pos.getBrother();
		if (tpos == pos){ // is maximum box? (bro is itself)
			space.toMemfree();
			query.toZero();
			if (query != n_node) space.insert(KeyElem<HyperPosition<S,D,L>, N >(pos.hyperpos,n_node));
			return;
		}
		if (!this->getAt(query, tpos)) return writeHyperPosition_routine(pos, n_node);
		//printf("super found pure at bro!\n");
		if (!query.isTrivial()){
			//if (!(n_node.trivialIsCoherrentWith(query.planeIntersection(pos)))) query = n_node; // return writeHyperPosition_routine(pos, n_node);
			//printf("brother is not trivial\n");
			if (!(n_node.trivialIsCoherrentWith(query.planeIntersection(pos)))){
				//printf("no! not coherrent!\n");
				if (n_node.otherTrivial().trivialIsCoherrentWith(query.planeIntersection(tpos))){
					//printf("brother holds opposite!\n");
					query = n_node.otherTrivial();
					n_node.TrivialMerge(tpos, query);
					tpos.toParent();
					/*n_node.show();
					query.show();
					tpos.show();*/
				} else {query = n_node; tpos = pos;}// printf("brother disaggrees!\n");}
			} else tpos = pos;
		}else{
			if (query != n_node) n_node.TrivialMerge(tpos, query);
			tpos.toParent();
		}
	}
	writeHyperPosition_routine(tpos, query);
}

// first step, ensures that the insertion is a complex node whenever possible
LFHTEMP void HierarchicalMergable<N,S,D,L>::insertMK2(HyperCursor<S,D,L> pos, N n_node, bool debug){
	N query, simpler, simpler2;
	HyperCursor<S,D,L> tpos = pos;

	if ((n_node.canSimplify())&&(n_node.simplifyAt(simpler2, tpos))){
        printf("Silly! tried to insert an un-simplified!\n");
        n_node.show();
        simpler2.show();
        tpos.show();
        n_node = simpler2;

	}

	// check if insertion is pointless
	if (this->getAt_update(query, tpos)){
        if (pos.mag == tpos.mag){
            if (query == n_node) return;
        }else{
            query.simplifyAt(simpler, pos);
            if (simpler == n_node) return;
        }
	}
    if (debug) printf("startInsert\n");
	// Combines with brother(s) if possible
	while(true){
		tpos = pos.getBrother();
		if (tpos == pos){ // is maximum box? (bro is itself)
			space.toMemfree();
			query.toZero();
			if (query != n_node) space.insert(KeyElem<HyperPosition<S,D,L>, N >(pos.hyperpos,n_node));
			if (debug) printf("end Insert\n");
			return;
		}
		if (!this->getAt_simplifyMK2(query, tpos)) break;
		n_node.simplifyAt(simpler2, tpos);
		if (debug) {printf("this simple is"); simpler2.show();}
		if (simpler2 != query){
            if (debug) {printf("at"); tpos.show(); printf("simple bro is:");query.show();}
            query.simplifyAt(simpler, pos);
            if (simpler != n_node){
                if (!n_node.canMerge()) break;
                if (!n_node.mergeInto(query, pos)) break;
                if (debug) {printf("mergbro is:");query.show();}
            }
            n_node = query;
		}
		pos.toParent();
	}

	if (debug){
        printf("Inserting ");  n_node.show();
        printf("At: "); pos.show();
	}
	writeHyperPositionMK2_routine(pos, n_node,debug);
	if (debug) printf("end Insert\n");
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::changeNodeSize_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const Tuple<S,D> &minpos, const M& merger){
	HyperCursor<S,D,L> result = ite->k;
	result.toMinBoxExclusiveContainer(minpos);
	if (result.mag == 0xFFFF) {
		//printf("reducing to nothingness!\n");
		space.remove(ite); ite.findLE(result); return;
	} // has matching minbox
	N prev = ite->d;
	if (prev.isTrivial()) {ite.ordering_preserving_reference().k = result.hyperpos; return;}
	unsigned int tmppl = merger.planeIntersection(prev,result);
	//printf("super reduction time: %i at:\n", tmppl);
	//ExOp::show(result);

	if (tmppl == 3u) {ite.ordering_preserving_reference().k = result.hyperpos; return;}
	// complex becomes simple
	N newins;
	newins.toTrivial(tmppl);
	--ite;
	if (!ite.isValid()){
		//printf("hit an invalid...!\n");
		N tmptmp; tmptmp.toZero();
		ite.findFirst();
		if (newins == tmptmp) {space.remove(ite); // zero in front, delete and beware!
			ite.findFirst();
			while((ite.isValid())&&(ite->d == tmptmp)) space.remove(ite);
			ite.toInvalid();
		}else{
			ite.ordering_preserving_reference().k = result.hyperpos;
			ite.ordering_preserving_reference().d = newins;
		}
		return;
	}else{
		if ((ite->d == newins)||(newins.trivialIsCoherrentWith(merger.planeIntersection(ite->d,result)))){
			++ite;
			space.remove(ite);
			ite.findLE(result);
		}else{
			++ite;
			ite.ordering_preserving_reference().k = result;
			ite.ordering_preserving_reference().d = newins;
		}
	}
	if (!result.findNext()) {
            printf("SHOULD NEVER BE E!\n"); LFH_exit(1);
    }
	result.toMinBoxExclusiveContainer(minpos);
	while (result.mag != 0xFFFF){
		unsigned int tmppl2 = merger.planeIntersection(prev,result);
		if (tmppl2 == 3u) {
			space.insert(KeyElem<HyperPosition<S,D,L>, N >(result, prev));
			ite.findGE(result);
			return;
		}

		if (tmppl != tmppl2){
			tmppl = tmppl2;
			newins.toTrivial(tmppl);
			space.insert(KeyElem<HyperPosition<S,D,L>, N >(result,newins));
			ite.findGE(result);
		}

		if (!result.findNext()) {
                printf("SHOULD NEVER BE !!! F\n"); LFH_exit(1);
        }
		result.toMinBoxExclusiveContainer(minpos);
	}
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::HyperPositionInflate_routine(typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator &ite, const M& merger){
	HyperCursor<S,D,L> cursor = ite->k;
	bool _debug_ = false;
	N toexp = ite->d;
	unsigned short mag;
	if (_debug_) {printf("Inflating %X: ",toexp());	cursor.show();}
	++ite;
	cursor.toLargeMinBoxContainer();

	bool hasNext = ite.isValid();
	Tuple<S,D> minOfNext;
	unsigned int tmppl;

	if (hasNext) {
		if (_debug_) {printf("NextLimit %X: ", ite->d()); HyperCursor<S,D,L>(ite->k).show();}
		minOfNext = HyperCursor<S,D,L>(ite->k).getMin();
		cursor.toMinBoxExclusiveContainer(minOfNext);
		--ite;
	}else ite.findLast();
	--ite;
	do{
		if (cursor.mag == 0xFFFF) printf("should be impossible... got intersecting nodes! fdadf\n");
		HyperCursor<S,D,L> brother = cursor.getBrother();

		if (_debug_) printf("inner! continue %c\n", (brother <= cursor)  ? 'Y' : 'N');
		if (brother > cursor) break;

		tmppl = merger.planeIntersection(toexp,brother);


		if (ite.isValid()){
		//	if (toexp == ite->d){ //
		//		cursor = ite->k;
		//		++ite;
		//		space.remove_mm(ite);
		//	}else{

			HyperCursor<S,D,L> prev(ite->k);
			if (_debug_) {printf("PrevFound: ");  prev.show();}
			//ite->k.show();
			mag = brother.commonContainer_mag(prev);
			if (_debug_) {
				printf("Cur: ");cursor.hyperpos.show();
				printf("mag for cur... %X\n", cursor.mag);
				printf("Cur Bro: ");  brother.show();
				printf("PE %X %X %X\n", mag, brother.mag, prev.mag);
			}
			if (mag > brother.mag){
				if (_debug_) printf("HH\n");
				if (tmppl == 3u) break;
				unsigned int tmppl2 = merger.planeIntersection(ite->d,brother);
				if (tmppl2 == tmppl){
					if (_debug_) printf("aggree... %c\n", (toexp.isTrivial())? 'Y' : 'N');
					if (toexp.isTrivial()){
						if (merger.planeIntersection(ite->d,cursor) == tmppl){
							cursor.toParent();
							return writeHyperPosition_routine(cursor, ite->d);
						}
					}
				} else if ((tmppl2 != 3)&&(toexp.isTrivial())){
					N tmptmp; tmptmp.toTrivial(tmppl2);
					tmptmp.TrivialMerge(cursor, toexp);
					cursor.toParent();
					if (_debug_) printf("did an trivial merge path3! %X\n", toexp());
					return writeHyperPosition_routine(cursor, toexp);
				} else {if (_debug_) printf("disaggree...\n"); break;}
				cursor.toParent();
				if (_debug_) printf("PFF\n");
			}else if (brother.hyperpos == ite->k){
				if (_debug_) printf("PFCdd\n");
				if (ite->d.isTrivial()){
					if (ite->d.trivialIsCoherrentWith(tmppl)){
					}else if (toexp.isTrivial()) {
						ite->d.TrivialMerge(cursor, toexp);
						++ite; ite.ordering_preserving_reference().k = cursor.hyperpos;
						if (_debug_) printf("did an trivial merge! %X\n", toexp());
						cursor.toParent();
						return writeHyperPosition_routine(cursor, toexp);
					}else break;
					cursor.toParent();
					space.remove_mm(ite);
					if ((ite.isValid())&&(ite->d == toexp)) {
						// previous matches! use it
						cursor = ite->k;
						++ite;
						if (_debug_) printf("auto remove previous, matches!\n");
						space.remove_mm(ite);
						--ite;
					}
				}else{
					// previous is brother, hence ite->k.planeIntersection(brother) == 3
					 // DONE!
					break;
				}
			}else break;
		//	}
		}else if (tmppl == 3u) break;
		else{
			// previous is non-existant (trivial 0)
			// *TODO*
			//HyperCursor<S,D,L> prev();
			printf("PFWTF\n");
			break;
		}
		//cursor.show();
		cursor.toLargeMinBoxContainer();
		if (hasNext) cursor.toMinBoxExclusiveContainer(minOfNext);
		if (cursor.mag == 0xFFFF){
			printf("got cornered!\n");
			ite->k.show();
			ite->d.show();
			if (hasNext) {printf("nextis:"); ExOp::show(minOfNext);}
			else printf("no next\n");
			LFH_ALIVE;LFH_exit(1);
		}
	}while(true);
	if (ite.isValid()) ++ite;
	else ite.findFirst();
	//printf("to %X\n", cursor.mag);
	ite.ordering_preserving_reference().k = cursor.hyperpos;
	++ite;
	if (!ite.isValid()) return;
	if (toexp.isTrivial()){


		HyperCursor<S,D,L> ncursor = HyperCursor<S,D,L>(ite->k);
		//if (ite->d.isTrivial()){
		//	if (ite->d.trivialIsCoherrentWith(toexp.planeIntersection(ncursor))) return writeHyperPosition_routine(ncursor, toexp);
		//}

		if (_debug_) {printf("recurmerge!\n"); ncursor.show();}
		HyperCursor<S,D,L> brother = ncursor.getBrother();
		if ((brother > ncursor)||(brother.getMinBox() < cursor.getMinBox())) return;
		tmppl = merger.planeIntersection(ite->d,brother);
		if (tmppl == 3u) return;
		if (ite->d.isTrivial()) ite->d.TrivialMerge(brother, toexp); // merging opposites!
		else if (toexp.trivialIsCoherrentWith(tmppl)) toexp = ite->d;
		else return;
		if (_debug_) printf("recurmerge what? %X\n", toexp());
		ncursor.toParent();
		return writeHyperPosition_routine(ncursor, toexp);
	}else{
		Tuple<S,D> damin;
		N tmptmp;
		bool change = false;
		if (_debug_) printf("end check\n");
		if (ite.isValid()) while(ite->d.isTrivial()){
			if (ite->d.trivialIsCoherrentWith(merger.planeIntersection(toexp,cursor = HyperCursor<S,D,L>(ite->k)))){
				tmptmp = ite->d;
				++ite;
				if (ite.isValid()) damin = HyperCursor<S,D,L>(ite->k).getMin();
				cursor.findNext();
				if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
				change = false;
				while(cursor.mag != 0xFFFF) {
					if (!tmptmp.trivialIsCoherrentWith(merger.planeIntersection(toexp,cursor))){
						if (ite.isValid()) --ite;
						else ite.findLast();
						ite.ordering_preserving_reference().k = cursor.hyperpos;
						if (_debug_) {printf("could not destroy %X at\n", ite->d()); cursor.show();}
						change = true;
						break;
					}
					cursor.findNext();
					if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
				}
				if (change) break;
				if (ite.isValid()) --ite;
				else ite.findLast();
				if (_debug_) {printf("Destroying %X\n", ite->d());	HyperCursor<S,D,L>(ite->k).show();}
				space.remove_pp(ite);
				change = true;
				if (!ite.isValid()) break;
				if (ite->d == toexp) {
					space.remove_pp(ite);
					if (!ite.isValid()) break;
				}
			} else break;
		}
		if (change){
			if ((hasNext = ite.isValid())) {damin = HyperCursor<S,D,L>(ite->k).getMin();	--ite;}
			else ite.findLast();
			cursor = HyperCursor<S,D,L>(ite->k);
			cursor.toLargeMinBoxContainer();
			if (hasNext) cursor.toMinBoxExclusiveContainer(damin);
			ite.ordering_preserving_reference().k = cursor.hyperpos;
			return HyperPositionInflate_routine(ite, merger);
		}
		if (_debug_) printf("nothing!\n");
	}
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::writeHyperPosition_routine(const HyperCursor<S,D,L> &where, const N& what, const M& merger){
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(this->space);
	HyperCursor<S,D,L> cursor;
	KeyElem<HyperPosition<S,D,L>, N > toins;
	N prevprop;
	N nextprop;
	cursor = where.getMaxBox();

	bool _debug_ = false;
	if (_debug_) {printf("Actual insertion: %X  mag %X\n", what(), where.mag);where.show();}


	ite.findLE(cursor);

	if (!ite.isValid())	{ // is empty and first box
		ite.findFirst();
		if (!ite.isValid())	{
			//printf("hh h haha!\n"); fflush(stdout);
			nextprop.toZero();
			prevprop.toZero();
		}else{
			//printf("hh h haha2!\n");
			if (where.isContainedBy(ite->k)){
				nextprop = ite->d;
				prevprop = ite->d;
				//printf("da routine!\n"); fflush(stdout);
				changeNodeSize_routine(ite,where.getMin(),merger);// reduces box!
				//this->show();
				//printf("da routine done!\n"); fflush(stdout);
			}else{
				nextprop.toZero();
				prevprop.toZero();
			}
		}
	}else{
		//printf("hh h hehe!\n");
		//printf("has ttfound!\n");
		// need to check if neighbor overlap with current contains query...
		//ite->k.show();
		++ite;
		if ((ite.isValid())&&(where.isContainedBy(ite->k))){
			nextprop = ite->d;
			//printf("hh h hoho! %X\n", ite->d());
			changeNodeSize_routine(ite,where.getMin(),merger);// reduces box!

			//space.show();
		}else{
			//printf("hh h hihi!\n", ite.isValid()?'Y':'N');
			if (!ite.isValid()) ite.findLast();
			else --ite;
			nextprop = ite->d;
			if (ite->k == where.hyperpos){
				if (ite->d == what) {/*printf("detected pointless!\n");*/return;} // pointless insert!
				space.remove_mm(ite);
			}else if (where.isContainedBy(ite->k)){
				//printf("got here! large\n");
				changeNodeSize_routine(ite,where.getMin(),merger);// reduces box!
			}else{
			//	printf("range rem here!\n");
			//	where.show();
			//	ite->k.show();
				cursor = where.getMinBox();
				space.removeRange(cursor.hyperpos, ite); // points to last that is not removed!
			}
		}
		if (ite.isValid()) {prevprop = ite->d; ++ite;}
		else {prevprop.toZero(); ite.findFirst();}
	}

	cursor = where.getMaxBox();



	if (_debug_) printf("Prev %X, Next %X, ToIns %X (cursize %i)\n", prevprop(), nextprop(), what(), space.getSize());

	if (nextprop != what){
		cursor = where;
		if (cursor.findNext()){
			if (_debug_) printf("gota\n");
			toins.d = nextprop;
			Tuple<S,D> damin;
			if (ite.isValid()){
				damin = HyperCursor<S,D,L>(ite->k).getMin();
				//if (cursor.mag == cursor.commonContainer_mag(tmp))
				cursor.toMinBoxExclusiveContainer(damin);
			}
			if (cursor.mag == 0xFFFF){ // no next, merge is possible!
				if (_debug_) printf("NO NEXT!!!\n");
				if (ite.isValid()){
					if (ite->d == what) { // merging to next!
						if (prevprop == what){
							//printf("remnext\n");
							space.remove_mm(ite);
                            if (!ite.isValid()){
                                if (_debug_) printf("end that end?\n");
                                return;
                            }
						}else{
							ite.ordering_preserving_reference().k = where.hyperpos;
						}
						if (_debug_) printf("INF call A\n");
						//this->checkIntegrity();
					}else if (prevprop != what){
						toins.d = what;
						toins.k = where.hyperpos;
						space.insert(toins);
						ite.findLE(where);
						if (_debug_) printf("INF call hAh\n");
					}else{
						--ite;
						//printf("INF call B\n");
						if (!ite.isValid()) return;
					}
					HyperPositionInflate_routine(ite,merger);
				}else if (prevprop != what){
					if (_debug_) printf("super espacedd!\n");
					toins.d = what;
					toins.k = where.hyperpos;
					space.insert(toins);
				} else{ if (_debug_) printf("miiiiised!\n");}
				return;
			}else{
				if (_debug_) printf("gotD %X\n", nextprop());
				if (nextprop.isTrivial()) {
					toins.k = cursor.hyperpos;
					space.insert(toins);// printf("nxt in: "); toins.k.show();
					if (_debug_) printf("got %X inserted\n", toins.d());
				}else{
					bool startin = true;

					do{
						if (_debug_) printf("Duper fun! %c\n", startin ? 'Y' : 'N');
						//space.show();
						unsigned int tmppl = merger.planeIntersection(nextprop,cursor);
						if (tmppl == 3u) {
							//printf("%X gave int3\n",  nextprop());
							//printf("insert at:"); cursor.show();
							toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); startin = false; break;
						}
						if (tmppl != merger.planeIntersection(what,cursor)){
							//printf("duper disagree! %i and %i\n", tmppl,  what.planeIntersection(cursor));
							toins.d.toTrivial(tmppl);
							toins.k = cursor.hyperpos;
							startin = false;
							//printf("Duper fun! %c\n", startin ? 'Y' : 'N');
							if (!startin) space.insert(toins);
							cursor.findNext();
							if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							do{
								if (cursor.mag == 0xFFFF) {
										ite.findGE(toins); // was just inserted, should be there
										++ite;
										if ((ite.isValid())&&(ite->d.planeIntersection() == tmppl)) space.remove(ite);
									break;
								}
								unsigned int tmppl2 = merger.planeIntersection(nextprop,cursor);
								if (tmppl2 == 3u) {toins.d = nextprop; toins.k = cursor.hyperpos; space.insert(toins); break;}
								if (tmppl2 != tmppl){
									tmppl = tmppl2;
									toins.d.toTrivial(tmppl);
									toins.k = cursor.hyperpos;
									//printf("Duper fun! N=%c\n", startin ? 'Y' : 'N');
									space.insert(toins);
								}
								cursor.findNext();
								if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
							}while(true);
							break;
						}
						cursor.findNext();
						if (ite.isValid()) cursor.toMinBoxExclusiveContainer(damin);
					} while(cursor.mag != 0xFFFF);
					if (_debug_) {printf("after duper\n");space.show();}
					if (startin == true){ // ended up deleting nodes, not inserting
						if (_debug_) {printf("got hhheehheh %X\n", ite->d());	ExOp::show(ite->k);}

						if (ite.isValid()){
							if (ite->d == what) { // merging to next!
								if (prevprop == what){
									if (_debug_) printf("remnext... reallytt?\n");
									space.remove_mm(ite);
                                    if (!ite.isValid()){
                                        if (_debug_) printf("end that end?\n");
                                        return;
									}
								}else if ((what.isTrivial())&&(what.trivialIsCoherrentWith(merger.planeIntersection(prevprop,where.hyperpos)))){
									if (_debug_) printf("simply! not eqq!\n");
									return writeHyperPosition_routine(ite->k, prevprop,merger);
								}else{
									ite.ordering_preserving_reference().k = where.hyperpos;
								}
								if (_debug_) printf("INF call A_hehe! does inflate!\n");
								//this->checkIntegrity();
							}else if (prevprop != what){
								if (_debug_) printf("not eqq!\n");
								toins.d = what;
								toins.k = where.hyperpos; if (_debug_) printf("INF call Bsds\n");
								space.insert(toins);
								ite.findLE(where);
							}else{
								--ite;
								if (_debug_) printf("INF call B\n");
								if (!ite.isValid()) return;
							}
							HyperPositionInflate_routine(ite,merger);
						}else if (prevprop != what){
							toins.d = what;
							toins.k = where.hyperpos;if (_debug_) printf("INF call Bh\n");
							space.insert(toins);
						}else{ if (_debug_) printf("got nowhre!\n");}
						return;
					}
				}
			}

		}
	}else{ // inflated node might be needed, but not done yet
		if (_debug_) printf("gotb\n");
		if (prevprop != what){
			toins.d = what;
			toins.k = where.hyperpos;
			space.insert(toins);
		}
		ite.findLE(where);
		//printf("INF call C\n");

		if (ite.isValid()) HyperPositionInflate_routine(ite,merger);
		return;
	}
	if (_debug_) printf("gotc\n");
	if (prevprop != what){
		//printf("add last %X != %X\n", what.anbox, prevprop.anbox );
//		if ((what.isTrivial())&&(what.trivialIsCoherrentWith(prevprop.planeIntersection(where)))){

//		}else{
		toins.d = what;
		toins.k = where.hyperpos;
		space.insert(toins);
		if (_debug_) {printf("again insert %X\n", toins.d());		space.show();}
//		}
		ite.findLE(toins.k);
	}else{
		if (_debug_) printf("would upgrade...\n");
		// prev node might need to be enlarged!
		ite.findLE(where);
		if (!ite.isValid()) return;
	}
	return HyperPositionInflate_routine(ite,merger);
}
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::insert(const HyperCursor<S,D,L> &pos, const N& n_node, const M& merger){
	N query;
	HyperCursor<S,D,L> tpos;
	if (!merger.shouldInsert(pos)) return;
	if (!n_node.isTrivial()) {query = n_node; tpos = pos;}
	else{
		//printf("inserting trivial%i\n", n_node());
		bool isPure = this->getAt(query,pos);
		if (isPure){

			if (n_node.trivialIsCoherrentWith(merger.planeIntersection(query,pos))) {
                //printf("super coherent!\n");
				return;
			}
		}

		tpos = pos.getBrother();
		if (tpos == pos){ // is maximum box? (bro is itself)
			space.toMemfree();
			query.toZero();
			if (query != n_node) space.insert(KeyElem<HyperPosition<S,D,L>, N >(pos.hyperpos,n_node));
			return;
		}
		if (!this->getAt(query, tpos)) return writeHyperPosition_routine(pos, n_node, merger);
		//printf("super found pure at bro!\n");
		if (!query.isTrivial()){
			//printf("brother is found %i\n", query());
			tpos.toParent();
			//if (!(n_node.trivialIsCoherrentWith(query.planeIntersection(pos)))) query = n_node; // return writeHyperPosition_routine(pos, n_node);
			//printf("brother is not trivial\n");
			/*if (!(n_node.trivialIsCoherrentWith(merger.planeIntersection(query,pos)))){
				//printf("no! not coherrent!\n");
				if (n_node.otherTrivial().trivialIsCoherrentWith(merger.planeIntersection(query,tpos))){
					//printf("brother holds opposite!\n");
					query = n_node.otherTrivial();
					n_node.TrivialMerge(tpos, query);
					tpos.toParent();
				} else {query = n_node; tpos = pos;}// printf("brother disaggrees!\n");}
			} else tpos = pos;*/
		}else{
			printf("brother trivial is found %i\n", query());
			if (query != n_node) n_node.TrivialMerge(tpos, query);
			tpos.toParent();
		}
	}
	writeHyperPosition_routine(tpos, query, merger);
}
// first step, ensures that the insertion is a complex node whenever possible
LFHTEMP template<class M> void HierarchicalMergable<N,S,D,L>::insertMK2(HyperCursor<S,D,L> pos, N n_node, const M& merger, bool need_simplification, bool debug){
	N query, simpler, simpler2;
	HyperCursor<S,D,L> tpos = pos;

	if ((need_simplification)&&(merger.simplifyAt(n_node, simpler, pos))){n_node = simpler;}

	// check if insertion is pointless
	if (this->getAt_update(query, tpos, merger)){
        //if (pos.mag == tpos.mag){
        //    if (query == n_node) return; june2019
        //}else{
            merger.simplifyAt(query,simpler, pos, debug);
            if (simpler == n_node) return;

        //}
	}
	// Combines with brother(s) if possible
	if (debug) printf("startInsert\n");
	while(n_node.canMerge()){
		tpos = pos.getBrother();
		if (tpos == pos){ // is maximum box? (bro is itself)
			space.toMemfree();
			query.toZero();
			if (query != n_node) space.insert(KeyElem<HyperPosition<S,D,L>, N >(pos.hyperpos,n_node));
			if (debug) printf("end  Insert\n");
			return;
		}
		if (!this->getAt_simplifyMK2(query, tpos,merger)) break;
		merger.simplifyAt(n_node,simpler2, tpos, debug); if (debug) {tpos.show(); merger.showThat(query); merger.showThat(simpler2);}
		if (query != simpler2) {
            merger.simplifyAt(query,simpler,pos, debug);
            if (simpler != n_node){
                if (!merger.mergeInto(n_node,query, pos, debug)) break;
            }
            n_node = query;
		}
		pos.toParent();
	}
	if (debug){
        printf("Inserting ");  merger.showThat(n_node);
        printf("At: "); pos.show();
	}
	writeHyperPositionMK2_routine(pos, n_node,merger,debug);
	if (debug) printf("end  Insert\n");
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::checkIntegrity()const{
	// checks that boxes spans entire space exactly:
	typename RBTofDoom< KeyElem<HyperPosition<S,D,L>, N > >::Iterator ite(this->space);
	HyperCursor<S,D,L> mima;
	ite.findFirst();
	N prev;
	if (ite.isValid()){
		mima = HyperCursor<S,D,L>(ite->k).getMaxBox();
		prev = ite->d;
		++ite;
		while(ite.isValid()){
			HyperCursor<S,D,L> tmp = HyperCursor<S,D,L>(ite->k);
			if (tmp.getMinBox() <= mima){
				printf("Found error: Intersection!\n");
				--ite;
				HyperCursor<S,D,L>(ite->k).show();
				tmp.show();
				ExOp::show(ite->d);
				++ite;
				HyperCursor<S,D,L>(ite->k).getMaxBox().show();
				tmp.getMinBox().show();
				ExOp::show(ite->d);
				return false;
			}else if (prev == ite->d){
				printf("Found error: Double\n");
				--ite;
				HyperCursor<S,D,L>(ite->k).show();
				ExOp::show(ite->d);
				++ite;
				HyperCursor<S,D,L>(ite->k).show();
				ExOp::show(ite->d);
				return false;
			}
			prev = ite->d;
			mima = tmp.getMaxBox();
			++ite;
		}
	}
	return true;
}
LFHTEMP void HierarchicalMergable<N,S,D,L>::insert_with_bypass(const HyperCursor<S,D,L> &pos, const N& n_node){
	KeyElem<HyperPosition<S,D,L>, N > to_ins(pos,n_node);
	space.insert(to_ins);
	}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::first(){
	ite.findFirst();
	this->curkey.k.toLargestBoxPossible();
	ExOp::toZero(this->curkey.d);
	lastnode = this->curkey.d;
	if (!ite.isValid()) {mag = 0xFFFF; return true;}
	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	//tmp.show();
	//HyperCursor<S,D,L>(this->curkey.k).show();
	this->curkey.k.toMinBoxExclusiveContainer(limit);


	//printf("%i is mag\n", this->curkey.k.mag);
	if (this->curkey.k.mag == 0xFFFF){
		this->curkey.k = ite->k;
		this->curkey.d = ite->d;
		lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.k.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
		//	printf("CHANGE! new mag %X\n", mag);
		//	ite->k.show();
		}else{mag = 0xFFFF;} // printf("HIT END\n");
	}
	return true;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::next(){
	HyperCursor<S,D,L> cursor;
	if (!this->curkey.k.findNext()) return false;
	if ((this->curkey.k.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.k.toMinBoxExclusiveContainer(limit); mag = 0;}
	if (this->curkey.k.mag == 0xFFFF){ // no exclusive! use next
		this->curkey.k = ite->k;
		this->curkey.d = ite->d;
		lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.k.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
			//printf("CHANGE! new mag %X\n", mag);
			//ite->k.show();
		}else{mag = 0xFFFF;} // printf("HIT END\n");

	}else{
		switch(lastnode.planeIntersection(this->curkey.k)){
			case 1: this->curkey.d.toZero(); break;
			case 2: this->curkey.d.toOne(); break;
			default: this->curkey.d = lastnode;
		}
	}
	return true;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::first(const HyperCursor<S,D,L> &in_box){
	ite.findGE(in_box.getMinBox());
	this->curkey.k = in_box.hyperpos;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) ExOp::toZero(this->curkey.d);
		else{
			lastnode = ite->d;
			switch(lastnode.planeIntersection(this->curkey.k)){
				case 1: this->curkey.d.toZero(); break;
				case 2: this->curkey.d.toOne(); break;
				default: this->curkey.d = lastnode;
			}
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
		lastnode = ite->d;
		switch(lastnode.planeIntersection(in_box)){
			case 1: this->curkey.d.toZero(); break;
			case 2: this->curkey.d.toOne(); break;
			default: this->curkey.d = lastnode;
		}
		mag = 0xFFFF;
		return true;
	}
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(this->curkey.d);
		lastnode = this->curkey.d;
		ite.findFirst();
	}else{
		this->curkey.d = ite->d;
		lastnode = ite->d;
		if (in_box.isContainedBy(ite->k)) {
			switch(lastnode.planeIntersection(in_box)){
				case 1: this->curkey.d.toZero(); break;
				case 2: this->curkey.d.toOne(); break;
				default: this->curkey.d = lastnode;
			}
			mag = 0xFFFF;
			return true;
		}
		++ite;
	}

	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	mag =0;
	this->curkey.k.toMinBoxExclusiveContainer(limit);
	if (this->curkey.k.mag == 0xFFFF){
		this->curkey.k = ite->k;
		this->curkey.d = ite->d;
		lastnode = ite->d;
		++ite; // guarrantied to be valid
        if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag =this->curkey.k.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else{
		switch(lastnode.planeIntersection(this->curkey.k)){
			case 1: this->curkey.d.toZero(); break;
			case 2: this->curkey.d.toOne(); break;
			default: this->curkey.d = lastnode;
		}
	}
	return true;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::next(const HyperCursor<S,D,L> &in_box){
	if (!this->next()) return false;
	return (in_box.getMaxBox() >= (*this).curkey.k);
}
/*
LFHTEMP bool HierarchicalMergable<N,S,D,L>::Iterator::first(){
	ite.findFirst();
	this->curkey.k.toLargestBoxPossible();
	ExOp::toZero(this->curkey.d);
	lastnode = this->curkey.d;
	if (!ite.isValid()) {mag = 0xFFFF; return true;}
	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	//tmp.show();
	//HyperCursor<S,D,L>(this->curkey.k).show();
	this->curkey.k.toMinBoxExclusiveContainer(limit);


	//printf("%i is mag\n", this->curkey.k.mag);
	if (this->curkey.k.mag == 0xFFFF){
		this->curkey.k = ite->k;
		this->curkey.d = ite->d;
		lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.k.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
		}else{mag = 0xFFFF;} // printf("HIT END\n");
	}
	return true;
}*/
LFHTEMP bool HierarchicalMergable<N,S,D,L>::Iterator::next(){
	HyperCursor<S,D,L> cursor;
	if (!this->curkey.findNext()) return false;
	if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}
	if (this->curkey.mag == 0xFFFF){ // no exclusive! use next
		this->curkey = ite->k;
		simplenode = lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else lastnode.simplifyAt(simplenode, this->curkey);
	return true;
}
LFHTEMP HierarchicalMergable<N,S,D,L>::Iterator::operator bool (){
	ite.findGE(in_box.getMinBox());
	this->curkey = in_box.hyperpos;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(simplenode); ExOp::toZero(lastnode);}
		else{
			lastnode = ite->d;
			lastnode.simplifyAt(simplenode, this->curkey);
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
		lastnode = ite->d;
		lastnode.simplifyAt(simplenode, in_box);
		mag = 0xFFFF;
		return true;
	}
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(simplenode);
		lastnode = simplenode;
		ite.findFirst();
	}else{
		lastnode = simplenode = ite->d;
		if (in_box.isContainedBy(ite->k)) {
            lastnode.simplifyAt(simplenode, in_box);
			mag = 0xFFFF;
			return true;
		}
		++ite;
	}

	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	mag =0;
	this->curkey.toMinBoxExclusiveContainer(limit);
	if (this->curkey.mag == 0xFFFF){
		this->curkey = ite->k;
		lastnode = simplenode = ite->d;
		++ite; // guarrantied to be valid
		if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag =this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else lastnode.simplifyAt(simplenode, this->curkey);
	return true;
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::Iterator::operator++(int){
	if (!this->next()) return false;
	return (in_box.getMaxBox() >= (*this).curkey);
}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::IteratorMerger<M>::next(){
	HyperCursor<S,D,L> cursor;
	if (!this->curkey.findNext()) return false;
	if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}
	if (this->curkey.mag == 0xFFFF){ // no exclusive! use next
		this->curkey = ite->k;
		simplenode = lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else mergescope.simplifyAt(lastnode,simplenode, this->curkey);
return true;}
LFHTEMP template<class M> HierarchicalMergable<N,S,D,L>::IteratorMerger<M>::operator bool (){
	ite.findGE(in_box.getMinBox());
	this->curkey = in_box.hyperpos;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(simplenode); ExOp::toZero(lastnode);}
		else{
			lastnode = ite->d;
			mergescope.simplifyAt(lastnode,simplenode, this->curkey);
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
		lastnode = ite->d;
		mergescope.simplifyAt(lastnode,simplenode, in_box);
		mag = 0xFFFF;
		return true;
	}
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(simplenode);
		lastnode = simplenode;
		ite.findFirst();
	}else{
		lastnode = simplenode = ite->d;
		if (in_box.isContainedBy(ite->k)) {
            mergescope.simplifyAt(lastnode,simplenode, in_box);
			mag = 0xFFFF;
			return true;
		}
		++ite;
	}

	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	mag =0;
	this->curkey.toMinBoxExclusiveContainer(limit);
	if (this->curkey.mag == 0xFFFF){
		this->curkey = ite->k;
		lastnode = simplenode = ite->d;
		++ite;
        if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag =this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else mergescope.simplifyAt(lastnode,simplenode, this->curkey);
return true;}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::IteratorMerger<M>::operator++(int){
	if (!this->next()) return false;
return (in_box.getMaxBox() >= (*this).curkey);}

LFHTEMP template<class F> bool HierarchicalMergable<N,S,D,L>::IteratorFunc<F>::next(){
	while(true){
        HyperCursor<S,D,L> cursor;
        if (!this->curkey.findNext()) return false;
        if (in_box.getMaxBox() < (*this).curkey) return false;
        if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}
        if (this->curkey.mag == 0xFFFF){ // no exclusive! use next
            this->curkey = ite->k;
            simplenode = lastnode = ite->d;
            ++ite;
            if (ite.isValid()) {
                HyperCursor<S,D,L> tmp(ite->k);
                mag = this->curkey.commonContainer_mag(tmp);
                HyperCursor<S,D,L>::decr_mag(mag);
                limit = tmp.getMin();
            }else{mag = 0xFFFF;}
        }else lastnode.simplifyAt(simplenode, this->curkey);
        if (func( this->curkey)) return true;
	}
}
LFHTEMP template<class F> HierarchicalMergable<N,S,D,L>::IteratorFunc<F>::operator bool (){
	ite.findGE(in_box.getMinBox());
	this->curkey = in_box.hyperpos;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(simplenode); ExOp::toZero(lastnode);}
		else{
			lastnode = ite->d;
			lastnode.simplifyAt(simplenode, this->curkey);
		}
		mag = 0xFFFF;
		if (!func( this->curkey)) return (*this)++;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
		lastnode = ite->d;
		lastnode.simplifyAt(simplenode, in_box);
		mag = 0xFFFF;
		if (!func( this->curkey)) return (*this)++;
		return true;
	}
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(simplenode);
		lastnode = simplenode;
		ite.findFirst();
	}else{
		lastnode = simplenode = ite->d;
		if (in_box.isContainedBy(ite->k)) {
            lastnode.simplifyAt(simplenode, in_box);
			mag = 0xFFFF;
			if (!func( this->curkey)) return (*this)++;
			return true;
		}
		++ite;
	}

	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	mag =0;
	this->curkey.toMinBoxExclusiveContainer(limit);
	if (this->curkey.mag == 0xFFFF){
		this->curkey = ite->k;
		lastnode = simplenode = ite->d;
		++ite;
        if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag =this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else lastnode.simplifyAt(simplenode, this->curkey);

	if (!func( this->curkey)) return (*this)++;
return true;}
LFHTEMP template<class F> bool HierarchicalMergable<N,S,D,L>::IteratorFunc<F>::operator++(int){
	if (!this->next()) return false;
return (in_box.getMaxBox() >= (*this).curkey);}

LFHTEMP template<class F,class M> bool HierarchicalMergable<N,S,D,L>::IteratorFuncMerger<F,M>::next(){
	HyperCursor<S,D,L> cursor;
	if (!this->curkey.findNext()) return false;
	if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}
	if (this->curkey.mag == 0xFFFF){ // no exclusive! use next
		this->curkey = ite->k;
		simplenode = lastnode = ite->d;
		++ite;
		if (ite.isValid()) {
			HyperCursor<S,D,L> tmp(ite->k);
			mag = this->curkey.commonContainer_mag(tmp);
			HyperCursor<S,D,L>::decr_mag(mag);
			limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else mergescope.simplifyAt(lastnode,simplenode, this->curkey);
return true;}
LFHTEMP template<class F,class M> HierarchicalMergable<N,S,D,L>::IteratorFuncMerger<F,M>::operator bool (){
	ite.findGE(in_box.getMinBox());
	this->curkey = in_box.hyperpos;
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) {ExOp::toZero(simplenode); ExOp::toZero(lastnode);}
		else{
			lastnode = ite->d;
			mergescope.simplifyAt(lastnode,simplenode, this->curkey);
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
		lastnode = ite->d;
		mergescope.simplifyAt(lastnode,simplenode, in_box);
		mag = 0xFFFF;
		return true;
	}
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(simplenode);
		lastnode = simplenode;
		ite.findFirst();
	}else{
		lastnode = simplenode = ite->d;
		if (in_box.isContainedBy(ite->k)) {
            mergescope.simplifyAt(lastnode,simplenode, in_box);
			mag = 0xFFFF;
			return true;
		}
		++ite;
	}

	HyperCursor<S,D,L> tmp(ite->k);
	limit = tmp.getMin();
	mag =0;
	this->curkey.toMinBoxExclusiveContainer(limit);
	if (this->curkey.mag == 0xFFFF){
		this->curkey = ite->k;
		lastnode = simplenode = ite->d;
		++ite;
        if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag =this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
		}else{mag = 0xFFFF;}
	}else mergescope.simplifyAt(lastnode,simplenode, this->curkey);
return true;}
LFHTEMP template<class F, class M> bool HierarchicalMergable<N,S,D,L>::IteratorFuncMerger<F,M>::operator++(int){
	if (!this->next()) return false;
return (in_box.getMaxBox() >= (*this).curkey);}

LFHTEMP bool HierarchicalMergable<N,S,D,L>::IteratorMima::next(){
    ite.findGE(in_box);
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) ExOp::toZero(this->simplenode);
		else{
			lastnode = ite->d;
            this->simplenode = lastnode;
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
   //     printf("iscontained!\n");
        //cur = ite->k;
		lastnode = ite->d;
        this->simplenode = lastnode;
        this->curkey = ite->k;
        ++ite;
		if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag = this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
        }else{mag = 0xFFFF;}
		return true;
	}
	HyperCursor<S,D,L> tmp2(ite->k);
	limit = tmp2.getMin();
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(this->simplenode);
		lastnode = this->simplenode;
		ite.findFirst();
		for(uint32_t i=0;i< D;i++) this->curkey.hyperpos[i];
		this->curkey.hyperpos[0] |=1;
		this->curkey.mag=0;
		this->curkey.toMinBoxExclusiveContainer(limit);
	}else{
 //       printf("backtrack to\n");
		this->simplenode = ite->d;
		lastnode = ite->d;
		this->curkey = ite->k;
		ite->k.show();
		++ite;
	}
    mag = this->curkey.commonContainer_mag(tmp2);
    HyperCursor<S,D,L>::decr_mag(mag);
	return true;
}
LFHTEMP HierarchicalMergable<N,S,D,L>::IteratorMima::operator bool (){
    for(uint32_t i=0;i<D;i++) in_box.hyperpos[i] = mima[0][i];
    in_box.hyperpos[0] |= 1; in_box.mag =0;
    uint16_t mag_box = in_box.largest_Container_mag_rect(mima[0], mima[1]); // box *must* share minimum, but...
    in_box.setMagnitude(mag_box);
    return (this->next());
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::IteratorMima::operator++(int){
    while (this->curkey.findNext()){
        if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}

        if (this->curkey.mag != 0xFFFF){ // no exclusive! use next
            if (!this->curkey.isIntersectingRect(mima[0],mima[1])) continue;
            lastnode.simplifyAt(simplenode, this->curkey);
            return true;
        }
        // check intersection!
        this->curkey = ite->k;
        if (this->curkey.isIntersectingRect(mima[0],mima[1])){
            this->simplenode = ite->d;
            lastnode = ite->d;
            ++ite;
            if (ite.isValid()) {
                HyperCursor<S,D,L> tmp(ite->k);
                mag = this->curkey.commonContainer_mag(tmp);
                HyperCursor<S,D,L>::decr_mag(mag);
                limit = tmp.getMin();
            }else{mag = 0xFFFF;}
            return true;
        }
        break;
    }

    if (this->curkey.mag == 0xFFFF) return false;
    HyperCursor<S,D,L> cur = this->curkey;
    in_box = cur.NextTinyInBox2(mima[0], mima[1]);

    //TODO!
    if (in_box.mag == 0xFFFF) return(false);
	return (this->next());
}

LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::IteratorMimerger<M>::next(){
    ite.findGE(in_box);
	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) ExOp::toZero(this->simplenode);
		else{
			lastnode = ite->d;
            this->simplenode = lastnode;
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
   //     printf("iscontained!\n");
        //cur = ite->k;
		lastnode = ite->d;
        this->simplenode = lastnode;
        this->curkey = ite->k;
        ++ite;
		if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag = this->curkey.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
        }else{mag = 0xFFFF;}
		return true;
	}
	HyperCursor<S,D,L> tmp2(ite->k);
	limit = tmp2.getMin();
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(this->simplenode);
		lastnode = this->simplenode;
		ite.findFirst();
		for(uint32_t i=0;i< D;i++) this->curkey.hyperpos[i];
		this->curkey.hyperpos[0] |=1;
		this->curkey.mag=0;
		this->curkey.toMinBoxExclusiveContainer(limit);
	}else{
 //       printf("backtrack to\n");
		this->simplenode = ite->d;
		lastnode = ite->d;
		this->curkey = ite->k;
		ite->k.show();
		++ite;
	}
    mag = this->curkey.commonContainer_mag(tmp2);
    HyperCursor<S,D,L>::decr_mag(mag);
	return true;
}

LFHTEMP template<class M> HierarchicalMergable<N,S,D,L>::IteratorMimerger<M>::operator bool (){
    for(uint32_t i=0;i<D;i++) in_box.hyperpos[i] = mima[0][i];
    in_box.hyperpos[0] |= 1; in_box.mag =0;
    uint16_t mag_box = in_box.largest_Container_mag_rect(mima[0], mima[1]); // box *must* share minimum, but...
    in_box.setMagnitude(mag_box);
    return (this->next());
}
LFHTEMP template<class M> bool HierarchicalMergable<N,S,D,L>::IteratorMimerger<M>::operator++(int){
    while (this->curkey.findNext()){
        if ((this->curkey.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.toMinBoxExclusiveContainer(limit); mag = 0;}

        if (this->curkey.mag != 0xFFFF){ // no exclusive! use next
            if (!this->curkey.isIntersectingRect(mima[0],mima[1])) continue;
            mergescope.simplifyAt(lastnode,simplenode, this->curkey);
            return true;
        }
        // check intersection!
        this->curkey = ite->k;
        if (this->curkey.isIntersectingRect(mima[0],mima[1])){
            this->simplenode = ite->d;
            lastnode = ite->d;
            ++ite;
            if (ite.isValid()) {
                HyperCursor<S,D,L> tmp(ite->k);
                mag = this->curkey.commonContainer_mag(tmp);
                HyperCursor<S,D,L>::decr_mag(mag);
                limit = tmp.getMin();
            }else{mag = 0xFFFF;}
            return true;
        }
        break;
    }

    if (this->curkey.mag == 0xFFFF) return false;
    HyperCursor<S,D,L> cur = this->curkey;
    in_box = cur.NextTinyInBox2(mima[0], mima[1]);

    //TODO!
    if (in_box.mag == 0xFFFF) return(false);
	return (this->next());
}

LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::first(const FUNCTOR &f){
	if (!this->first()) return false;
	if (f((*this).curkey)) return true;
	else return this->next(f);
}
LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::next(const FUNCTOR &f){
	do {
		if (!this->next()) return false;
	}while(!f((*this).curkey));
	return true;
}
LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::first(const HyperCursor<S,D,L> &in_box, const FUNCTOR &f){
	if (!this->first(in_box)) return false;
	if (f((*this).curkey)) return true;
	else return this->next(in_box,f);
}
LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::next(const HyperCursor<S,D,L> &in_box, const FUNCTOR &f){
	if (!this->next(f)) return false;
	return (in_box.getMaxBox() >= (*this).curkey.k);
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::first(const Tuple< Tuple<S,D>, 2u> &mima){
    HyperCursor<S,D,L> smallest;
    for(uint32_t i=0;i<D;i++) smallest.hyperpos[i] = mima[0][i];
    smallest.hyperpos[0] |= 1; smallest.mag =0;

    uint16_t mag_box = smallest.largest_Container_mag_rect(mima[0], mima[1]); // box *must* share minimum, but...
    smallest.setMagnitude(mag_box);
    return (this->first_inrect_routine(smallest));
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::next(const Tuple< Tuple<S,D>, 2u> &mima){
    while (this->curkey.k.findNext()){
        if ((this->curkey.k.mag >= mag)&&(mag != 0xFFFF)) {this->curkey.k.toMinBoxExclusiveContainer(limit); mag = 0;}

        if (this->curkey.k.mag != 0xFFFF){ // no exclusive! use next
            if (!this->curkey.k.isIntersectingRect(mima[0],mima[1])) continue;
            switch(lastnode.planeIntersection(this->curkey.k)){
                case 1: this->curkey.d.toZero(); break;
                case 2: this->curkey.d.toOne(); break;
                default: this->curkey.d = lastnode;
            }
            return true;
        }
        // check intersection!
        this->curkey.k = ite->k;
        if (this->curkey.k.isIntersectingRect(mima[0],mima[1])){
            this->curkey.d = ite->d;
            lastnode = ite->d;
            ++ite;
            if (ite.isValid()) {
                HyperCursor<S,D,L> tmp(ite->k);
                mag = this->curkey.k.commonContainer_mag(tmp);
                HyperCursor<S,D,L>::decr_mag(mag);
                limit = tmp.getMin();
            }else{mag = 0xFFFF;}
            return true;
        }
        break;
    }

    if (this->curkey.k.mag == 0xFFFF) return false;
    HyperCursor<S,D,L> cur = this->curkey.k;
    //printf("  get next: "); cur.show();
    HyperCursor<S,D,L> smallest = cur.NextTinyInBox2(mima[0], mima[1]);

    //TODO!
    if (smallest.mag == 0xFFFF) {
        //printf("got the end!\n");
        return(false);
    }
    //printf("got: "); smallest.show();
 /*   printf("[%X-%X,%X-%X,%X-%X]\n", mima[0][0], mima[1][0], mima[0][1], mima[1][1], mima[0][2], mima[1][2]);
    printf("rein!\n");
    smallest.show();*/
/*
    unsigned short mag = smallest.commonContainer_mag(cur); // max magnitude, due to preceding
    //printf("next: mag is %i\n", mag);
    HyperCursor<S,D,L>::decr_mag(mag);
    unsigned short tmag = smallest.largest_Container_mag_rect(mima[0], mima[1]);
    cur = smallest;
    cur.setMagnitude((tmag < mag) ? tmag : mag);*/
	return (this->first_inrect_routine(smallest));
}
LFHTEMP bool HierarchicalMergable<N,S,D,L>::KeyIterator::first_inrect_routine(const HyperCursor<S,D,L> &in_box){

    ite.findGE(in_box);

	if (!ite.isValid()){
		ite.findLast();
		if (!ite.isValid()) ExOp::toZero(this->curkey.d);
		else{
			lastnode = ite->d;
            this->curkey.d = lastnode;
		}
		mag = 0xFFFF;
		return true;
	}

	if (in_box.isContainedBy(ite->k)){
   //     printf("iscontained!\n");
        //cur = ite->k;
		lastnode = ite->d;
        this->curkey.d = lastnode;
        this->curkey.k = ite->k;
        ++ite;
		if (ite.isValid()) {
            HyperCursor<S,D,L> tmp(ite->k);
            mag = this->curkey.k.commonContainer_mag(tmp);
            HyperCursor<S,D,L>::decr_mag(mag);
            limit = tmp.getMin();
        }else{mag = 0xFFFF;}
		return true;
	}
	HyperCursor<S,D,L> tmp2(ite->k);
	limit = tmp2.getMin();
	--ite;
	if (!ite.isValid()) {
		ExOp::toZero(this->curkey.d);
		lastnode = this->curkey.d;
		ite.findFirst();
		for(uint32_t i=0;i< D;i++) this->curkey.k.hyperpos[i];
		this->curkey.k.hyperpos[0] |=1;
		this->curkey.k.mag=0;
		this->curkey.k.toMinBoxExclusiveContainer(limit);
	}else{
 //       printf("backtrack to\n");
		this->curkey.d = ite->d;
		lastnode = ite->d;
		this->curkey.k = ite->k;
		ite->k.show();
		++ite;
	}
    mag = this->curkey.k.commonContainer_mag(tmp2);
    HyperCursor<S,D,L>::decr_mag(mag);
	return true;
}
LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::firstBin(const FUNCTOR &f){
    return false;
    /*if (!ite->findFirst(f,mima)) false;

	if (!this->first()) return false;
	if (f((*this).curkey)) return true;
	else return this->next(f);*/
}
LFHTEMP template<class FUNCTOR> bool HierarchicalMergable<N,S,D,L>::KeyIterator::nextBin(const FUNCTOR &f){
	return false;
	/*do {
		if (!this->next()) return false;
	}while(!f((*this).curkey));
	return true;*/
}

#undef LFHTEMP
#define LFHTEMP template<class C,class I, unsigned int D, unsigned int L>


LFHTEMP HyperPaint<C,I,D,L>& HyperPaint<C,I,D,L>::toZero(){ExOp::toZero(data);return *this;}

LFHTEMP bool HyperPaint<C,I,D,L>::mergeInto(HyperPaint<C,I,D,L>& brother_io, const HyperCursor<I,D,L> &this_loc) const{
return false;}

LFHTEMP bool HyperPaint<C,I,D,L>::simplifyAt(HyperPaint<C,I,D,L>& fout, const HyperCursor<I,D,L> &location) const{
    fout = *this;
return false;}

LFHTEMP ERRCODE HyperPaint<C,I,D,L>::save(FILE* f) const{return ExOp::save(data,f);}
LFHTEMP ERRCODE HyperPaint<C,I,D,L>::load(FILE* f){return ExOp::load(data, f);}
LFHTEMP void HyperPaint<C,I,D,L>::show(FILE* f, int level)const{ExOp::show(data,f,level);}

#undef LFHTEMP
#define LFHTEMP template<unsigned int D>
LFHTEMP void GaussCloudNode<D>::show(FILE* f, int level) const{
    fprintf(f,"w=%f\t", weight);
    ExOp::show(center, f,level);
    fprintf(f, "\t+-\t");
    ExOp::show(box_radii, f,level);
}


#undef LFHTEMP
#define LFHTEMP template<class N, unsigned int D, class I>
LFHTEMP GaussCloud<N,D,I>::GaussCloud(const Vector< KeyElem<N, Tuple<double> > >  &points){
    if (points.getSize() == 0) return;
    data.setSize(points.getSize());
    GaussElem<Tuple<double,D> > gaussian; gaussian.toZero();
    unsigned int i;
    gaussian = new GaussElem<Tuple<double,D> >(points[0]);
    for(i=1;i<points.getSize();i++) gaussian +=  new GaussElem<Tuple<double,D> >(points.data[i]);
}

LFHTEMP uint32_t GaussCloud<N,D,I>::findClosest(uint32_t query, const Trianglix<double, D> &proj) const{
    // needs to find all such that threshold > (mu_i-query)^t proj (mu_i-query) / (1.0 + query_w / w_i)
    // also finds the minimum for all i, if no none are below the threshold


return 0;}

LFHTEMP uint32_t GaussCloud<N,D,I>::wrWithinRange(uint32_t* wr_buf, double threshold, uint32_t query, const Trianglix<double, D> &proj) const{
    // needs to find all such that threshold > (mu_i-query)^t proj (mu_i-query) / (1.0 + query_w / w_i)
    // also finds the minimum for all i, if no none are below the threshold


return 0;}
/*
LFHTEMP void GaussCloud<N,D,I>::initFrom(const Vector< KeyElem<N, Tuple<double> > > &points){
    Vector<KeyElem<HyperCursor<uint32_t, D>, uint32_t> > treeinit;



}*/

LFHTEMP void GaussCloud<N,D,I>::insert(const N& node, const Tuple<double>& coor, const Trianglix<double> &var, double weight){
    GaussCloudNode<D>& newelem = data.addEntry(node);

    newelem.center = proj_matrix * coor;
    Trianglix<double, 0u> cov = var.mkInnerMult(proj_matrix);
    printf("well\n");
    var.show();
    cov.show();
    double tmp[3];
    if (D < 4){ // x cov^{-1} x = r^2
        if (D == 1) newelem.box_radii[0] = cov.data[0];
        else if (D == 2){
            tmp[0] = cov.data[1] * cov.data[1];
            tmp[1] = cov.data[2] * cov.data[0] - tmp[0];
            newelem.box_radii[0] = (tmp[1])/ (cov.data[2] - tmp[0] / cov.data[0]);
            newelem.box_radii[1] = (tmp[1])/ (cov.data[0] - tmp[0] / cov.data[2]);
        }else{
            // TODO !

        }

    } // else todo!
    newelem.weight = weight;

    Tuple<uint32_t,D> projmima[2];
    HyperCursor<I, D> projbox;

    for(uint32_t i =0;i<D;i++){
        tmp[0] = erf(newelem.center[i] - Zcover * newelem.box_radii[i]);
        tmp[1] = erf(newelem.center[i] + Zcover * newelem.box_radii[i]);
        projmima[0][i] = (0x80000000 + (int32_t)(tmp[0] * pow(2.0, 32.0)));
        projmima[1][i] = (0x80000000 + (int32_t)(tmp[0] * pow(2.0, 32.0)));
    }
    hie.insertCover(projmima[0],projmima[1],node);
}

LFHTEMP Vector<N> GaussCloud<N,D,I>::queryEllipsoid(ThreadBase &tb, const Tuple<double> &center, const Trianglix<double> &radius) const{
    class Receiver{
        std::mutex mut;
    public:
    Query_HyperCube<I,D,0u,uint32_t> daquery;
    const Tuple<double>& center;
    const Trianglix<double> &radius;
    Vector<N> & target;
    const GaussCloud<N,D,I> & _this;
    Receiver(const GaussCloud<N,D,I> &__this, Vector<N>& _target, const Tuple<double> &_center, const Trianglix<double> &_radius): _this(__this),target(_target), center(_center),radius(_radius){}
        void operator()(KeyElem<HyperCursor<I,D>, uint32_t> node){
            /*if (!daquery(node)) return; // check rect,
            if ((radius.mkInnerMult(_this.data[node.d].center - center)) >= 1.0) return; // check higher-dim ellipse
            std::lock_guard<std::mutex> guard(mut);
            target.push_back(_this.data[node.d].k);*/
        }
    };
    Vector<N> target;
    Receiver curreceive(*this, target,center,radius);
    //hie.hierarchy.queryParallel(tb, curreceive.daquery, curreceive);
return target;}

LFHTEMP HyperPosition<I, 0u,0u> GaussCloud<N,D,I>::makeBox(const GaussElem<Tuple<double, 0u> > datapoint, double LL_bound) const{ HyperPosition<I, 0u,0u> fout;
    fout.setSize(3);


return fout;}


LFHTEMP Forest<uint32_t,3u> GaussCloud<N,D,I>::cluster_likelihood(const TMatrix< double > &data, double prior_scale){ Forest<uint32_t,3u> fout;
class Task{
    public:
    HeapTree<KeyElem<double,uint32_t> > mastertree;
    Tuple< HeapTree<KeyElem<double,uint32_t> > > nodetree;
    Tuple<uint32_t> thranges;
    uint32_t task;
    double shared_minimum;
    TMatrix<double> projection_denum;
    TMatrix<double> projection;
    TMatrix<double> brojection;
    Tuple< void* > curoutput;
    Tuple< Trianglix<double> > output_covar;
    const GaussCloud<N,D,I>& cloud;
    Tuple< GaussElem<Tuple<double> > > out_gausselem;
    const TMatrix<double>* tmatrix_input;

    Task(const GaussCloud<N,D,I>& _cloud):cloud(_cloud){ // reference initializations for shared scope
    }
    uint32_t operator()(uint32_t threadID){
        switch(task){
        case 0:{ // generate a random matrix
            Tuple<double,0u, TUPLE_FLAG_REMOTE_MEMORY> rdata; rdata.tup_size = projection.sizes[0];
            output_covar[threadID].setSize(rdata.tup_size).toZero();
            for(uint32_t i = thranges[threadID]; i <thranges[threadID+1];i++){ // , tb.updateProgress(threadID)
                rdata = projection.data + i * rdata.tup_size;
                for(uint32_t j = 0; j < rdata.tup_size ;j++) rdata[j] =sampleGaussian();
                for(uint32_t j = 0; j < rdata.tup_size ;j++)
                output_covar[threadID] += Trianglix<double>(rdata);
            }
        }break;
        case 1:{
            Tuple<double,0u, TUPLE_FLAG_REMOTE_MEMORY> rdata; rdata.tup_size = brojection.sizes[0];
            for(uint32_t i = thranges[threadID]; i <thranges[threadID+1];i++){ // , tb.updateProgress(threadID)
                rdata = brojection.data + i * rdata.tup_size;
                rdata.toZero();
                //rdata = projection_denum * rdata;

                //if (auto ite = data.data[i].getIterator()) do{
                  //  for(uint32_t j = 0 ; j< brojection.getSize();j++)

                //}while(ite++);
                // output_covar[threadID] += Trianglix<double>(rdata); to check
            }
        }break;
        case 10:{ // compute prior covariance
            out_gausselem[threadID].setSize(tmatrix_input->sizes[0]).toZero();
            const double *target;
            if (auto rite = foreach.getRangeIterator(threadID, tmatrix_input->sizes[1],true)){
                target = tmatrix_input->data() + rite() * tmatrix_input->sizes[0];
                do{
                    out_gausselem[threadID].addSample(target);
                    target += tmatrix_input->sizes[0];
                }while(rite++);
            }
        }break;
        //tb.finishProgress(threadID);
        }
    return 0;}
	} lotask(*this);
    foreach.startThreadArray();
    uint32_t thite;
    // Task 1, find the prior covariance
    lotask.out_gausselem.setSize(foreach.getThreadArraySize());
    lotask.tmatrix_input = &data;
    lotask.task = 10;
    printf("mat size %i x %i\n", data.sizes[0], data.sizes[1]);
    foreach.startProgress("Finding prior scale", data.sizes[1]);
    for(thite = foreach.getThreadArraySize()-1; thite> 0;thite--) foreach.submit(lotask, thite);
    foreach.submit_ThenWait(lotask, thite);
    foreach.printf("done!\n");
    for(thite = foreach.getThreadArraySize()-1; thite> 0;thite--) lotask.out_gausselem[0] += lotask.out_gausselem[thite];
    lotask.out_gausselem[0].getKernelDensity_Silverman(prior_scale).show();

    TMatrix<double> projection; projection.setSizes(data.sizes[1], data.sizes[1]).toOne();
    this->setProjectionMatrix(projection);

    foreach.stopThreadArray();
return fout;}


LFHTEMP void GaussCloud<N,D,I>::show(FILE* f, int level)const{
    fprintf(f,"Gausscloud if Zcover %f:\n", Zcover);
    data.show(f);
    fprintf(f,"Tree structure:\n");
    hie.show(f);
}


#undef LFHTEMP
#define LFHTEMP template<class N, class I>

LFHTEMP ERRCODE GaussCloud<N,1u,I>::setProjectionMatrix(const TMatrix<double>& dam){
    if (dam.sizes[1] != 1) {printf("Expected 1 dim only\n"); return 1;} proj_matrix = dam;
return 0;}
LFHTEMP void GaussCloud<N,1u,I>::insert(const N &node, const Tuple<double> &coor, const Trianglix<double> &var, double weight){

}
LFHTEMP void GaussCloud<N,1u,I>::remove(const N &node){

}



/*
#undef LFHTEMP
#define LFHTEMP template<class D>

LFHTEMP DoClusterLikelihoodRatioScope<D>::DoClusterLikelihoodRatioScope(uint32_t _nbthread):loc_merge(NULL){
        scp.setSize(_nbthread);
    }
LFHTEMP DoClusterLikelihoodRatioScope<D>::~DoClusterLikelihoodRatioScope(){toMemfree();}
LFHTEMP DoClusterLikelihoodRatioScope<D>& DoClusterLikelihoodRatioScope<D>::toMemfree(){
    if (loc_merge.getSize() != 0){
        loc_merge.toMemfree();
        delete[](stats);delete[](tom);delete[](links);
	}
return *this;}
LFHTEMP DoClusterLikelihoodRatioScope<D>& DoClusterLikelihoodRatioScope<D>::toSize(uint32_t nsize){
    if (loc_merge.getSize() == nsize) return *this;
    if (loc_merge.getSize() != 0){
        delete[](stats);delete[](tom);delete[](links);
	}
	loc_merge.setSize(nsize);
    stats = new GaussElem< Tuple<D, 0u > >[nsize];
    links = new unsigned int[nsize*2 -1];
    tom = new unsigned int[nsize];
    totsize = nsize;
    midsize = nsize >> 1;
return *this;}

LFHTEMP uint32_t DoClusterLikelihoodRatioScope<D>::operator()(uint32_t thread_no){
    uint32_t col, row;*/
   /* if (isInit){
        for(col=scp[thread_no].colrange[0];col < scp[thread_no].colrange[0];col++){
            // row > col loop
            if (col >= midsize) break;
            for(row=col+1;row<totsize;row++){
                if (row == midsize) row = totsize - col - 1;

            }
        }
        for(;col < scp[thread_no].colrange[0];col++){
            // row < col loop
            for(row=0;row<col;row++){
                if (row == ) row = midsize;
            }
        }
    }else{


    }*/
/*return 0;}

LFHTEMP template<class C,int NBREL,class TB> void DoClusterLikelihoodRatioScope<D>::run(Forest<C,NBREL>& target, const Vector< GaussElem< Tuple<double, 0u > > > &data, TB &tb){
    uint32_t i;
    this->toSize(data.getSize());

    // populate cloud in projected dimentions
    // find Schur complement associated with each nodes

    // for all node, list all neighboring node that are within a threshold, or is the closest among all.

    for(uint32_t i = 0;i< scp.getSize();i++){
        scp[i].colrange[0] = (i * data.getSize()) / scp.getSize();
        scp[i].colrange[1] = ((i+1) * data.getSize()) / scp.getSize();
    }

    for(i=scp.getSize()-1;i>0;i--) tb.startEvent(this, i);
    tb.startEvent_ThenWait(this,i);


}*/

/*
#undef LFHTEMP
#define LFHTEMP template<class STORETYPE, unsigned int DIMS, unsigned int LEAD, class NODES>

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::insert(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){
    printf("WTH!\n");
     // find parent



    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > to_in(pos,i_n);

  //  to_in.k.par_mag = get_parent_order( to_in);

    ((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->insert(to_in);


	}


LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::insert_partition(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){

	HyperCursor<STORETYPE,DIMS,LEAD> hcur(pos);

	vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > in_nodes;

	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > bi;
	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ba;

	 bi.k =  pos.getMinBox();
	 ba.k =  pos.getMaxBox();
	// hcur.show();

	//  printf("%c\t%c\n",hcur.getBrother() > b[0].k ? 'Y' : 'N', hcur.getBrother() < b[1].k ? 'Y' : 'N');

	 // Step 1: remove sub or equal boxes!


    partition.removeRange(bi,ba);

//	 this->intersectionInterval(in_nodes,bi,ba);
//	 bi.k.par_mag =65535;
//	 for(int i=0 ;i< in_nodes.size();i++) {
//		 if ((*ite).par_mag< bi.k.par_mag) bi.k.par_mag = (*ite).par_mag;
//		 ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ite);
//	 }

	 // Step 2: check brother, if same type, upgrade it and exit

	ba.k = hcur.getBrother();
	unsigned int qq = partition.find_index(ba.k);
	if (qq){ // use its parent!
		 bi.k.par_mag = partition.tree[qq].first.k.par_mag;

	if (partition.tree[qq].first.d == i_n) {
		hcur.toParent();
		partition.remove(ba);

		//if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.k.isRealParentDirect()){
		//	((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);
		//	ba.k = hcur;
		//	qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find(ba);
		//	bi.k.par_mag = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.k.par_mag;
		//	((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);
		//	} else ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);



	ba.k = hcur.getBrother();
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite = partition.find_first(ba);
	 while((daite.isValid())&&(ba.k == (*daite).k)){


		 if ((*daite).d != i_n) break;
		partition.remove(ba);
		hcur.toParent();

	 	ba.k = hcur.getBrother();

	 		 daite = partition.find_first(ba);

	 		}

	 ba.k = hcur;
	 ba.k.par_mag = bi.k.par_mag;

	 ba.d = i_n;
	partition.insert(ba);



	// TODO, update the parent pointer of ex-brother subboxes

	return;
	}
	}


    // Step 3: find parent
    ba.k = hcur;
     bi.k = pos;
    bi.d = i_n;
     bi.k.par_mag = 0xFFFF;



	// Step 4: if parent is 1 layer up, downgrade/erase it, ( and update parent

//	if (bi.k.par_mag != 65535) {
//	ba.k = bi.k.getRealParent();
//	qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find(ba);

//	if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.d != i_n) return; // ministep if parent exist and is of same type, exit

//	if (bi.k.isRealParentDirect()){
		// todo
	//	}
//	}



//	 if (in_nodes.size() == 0){
		 // we did not find anything, need to recover hierac parent!

		 // TODO!

//		 }

	partition.insert(bi);


	 b[0].k = hcur.getBrother();
	 unsigned int qq = ((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->find(b[0]);

	if ((qq)&&(((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->tree[qq].first.d == i_n)) {

	((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->remove(b[0]);
	hcur.toParent();
	 b[0].k = hcur.getBrother(); qq = ((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->find(b[0]);

	 while(qq != 0){

		 if (((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->tree[qq].first.d != i_n) break;

		((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->remove(b[0]);
		hcur.toParent();
	 	b[0].k = hcur.getBrother();

	 		 qq = ((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, pair<NODES, unsigned char> > >*)this)->find(b[0]);

	 		}

	 b[0].k = (HyperPosition<STORETYPE,DIMS,LEAD>)hcur;
	 b[0].d = i_n;
	((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES> >*)this)->insert(b[0]);

	}else{

	// needs to check if ancestor is the same

	((RBTofDoom<KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>, NODES> >*)this)->insert(KeyElem<HyperPosition<STORETYPE,DIMS,LEAD>,NODES>(pos, i_n));

	}

	}
//SETCMP_enum gagahaga(const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &mi,const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &ma );

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::insert_partition(const HyperCursor<STORETYPE,DIMS, LEAD> &pos, const NODES i_n){

	HyperCursor<STORETYPE,DIMS,LEAD> hcur(pos);

	vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > in_nodes;

	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > bi;
	 KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ba;

	 bi.k =  pos.getMinBox();
	 ba.k =  pos.getMaxBox();
	// hcur.show();

	//  printf("%c\t%c\n",hcur.getBrother() > b[0].k ? 'Y' : 'N', hcur.getBrother() < b[1].k ? 'Y' : 'N');

	 // Step 1: remove sub or equal boxes!

 //    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite =
 //   ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(bi.k);
 //   if(ite.isValid()){
//		 if ((*ite).par_mag< bi.k.par_mag) bi.k.par_mag = (*ite).par_mag;
//		 ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ite);
//    }

    ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->removeRange(bi,ba);

//	 this->intersectionInterval(in_nodes,bi,ba);
//	 bi.k.par_mag =65535;
//	 for(int i=0 ;i< in_nodes.size();i++) {
//		 if ((*ite).par_mag< bi.k.par_mag) bi.k.par_mag = (*ite).par_mag;
//		 ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ite);
//	 }

	 // Step 2: check brother, if same type, upgrade it and exit

	 ba.k = hcur.getBrother();
	 unsigned int qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_index(ba.k);
	if (qq){ // use its parent!
		 bi.k.par_mag = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.k.par_mag;

	if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.d == i_n) {
		hcur.toParent();
		((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);



	 ba.k = hcur.getBrother();
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator daite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(ba);
	 while((daite.isValid())&&(ba.k == (*daite).k)){


		 if ((*daite).d != i_n) break;
		((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->remove(ba);
		hcur.toParent();

	 	ba.k = hcur.getBrother();

	 		 daite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(ba);

	 		}

	 ba.k = hcur;
	 ba.k.par_mag = bi.k.par_mag;

	 ba.d = i_n;
	((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->insert(ba);



	// TODO, update the parent pointer of ex-brother subboxes

	return;
	}
	}



	 // Step 3: find parent



	// Step 4: if parent is 1 layer up, downgrade/erase it, ( and update parent

//	if (bi.k.par_mag != 65535) {
//	ba.k = bi.k.getRealParent();
//	qq = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find(ba);

//	if (((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->tree[qq].first.d != i_n) return; // ministep if parent exist and is of same type, exit

//	if (bi.k.isRealParentDirect()){
		// todo
	//	}
//	}

	// Step 5: insert box at last!

	 	bi.d = i_n;
		bi.k = pos;
        bi.k.par_mag = 0xFFFF;
//	 if (in_nodes.size() == 0){
		 // we did not find anything, need to recover hierac parent!

		 // TODO!

//		 }

	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->insert(bi);


	}
//SETCMP_enum gagahaga(const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &mi,const KeyElem<HyperPosition<unsigned int,3,0>, unsigned int> &ma );


LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::partitionIntersection(const HyperCursor<STORETYPE,DIMS,LEAD> &center, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const{

    // for now return what is in interval;

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > query;



}

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::fillvoid(){ // proceedure for initialization

    const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> minimum;
    typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->first();

    HyperCursor<STORETYPE,DIMS,LEAD> other;


    Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > redo;
    unsigned short mmag,tmag;

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmpin;
    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmpin2;

    while(ite.isValid()){

        typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2 = ite;

        tmpin = *ite;
        ++ite;
        --ite2;
        if ((ite.isValid())&&((*ite).d == tmpin.d)){
            //  multiple blocks!
            redo.push_back(tmpin);
            while((ite.isValid())&&((*ite).d == tmpin.d)) {
                tmpin2 = *ite;
                ++ite;
            }
            if ((!ite.isValid())&&(!ite2.isValid())){
                // all blocks are of the same type!
                LFH_exit(1);
            }

            if (ite.isValid()){
                mmag = tmpin.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite2).k);
                }
            }

            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
            if  (mmag > tmpin.k.mag) {
                tmpin.k.setMagnitude(mmag);
            }
            redo.push_back(tmpin);

            if (ite.isValid()){
                mmag = tmpin2.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin2.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin2.k.commonContainer_mag((*ite2).k);
                }
            }

            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);

            if  (mmag > tmpin2.k.mag) {
                tmpin2.k.setMagnitude(mmag);
            }

            if (tmpin.k != tmpin2.k) {
                redo.push_back(tmpin2);
                // need to fill the space between  tmpin2 and  tmpin


                while(true){
                ++(tmpin.k);
                if (ite.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite).k);

                    if (ite2.isValid()){
                        tmag = tmpin.k.commonContainer_mag((*ite2).k);
                    if (tmag < mmag) mmag = tmag;
                }
                }else{
                    if (ite2.isValid()){
                        mmag = tmpin.k.commonContainer_mag((*ite2).k);
                    }else break;
                }
                HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
                if  (mmag > tmpin.k.mag) {
                    tmpin.k.setMagnitude(mmag);
                    redo.push_back(tmpin);
                } else break;


                }

                while(true){
                --(tmpin2.k);

                if (ite.isValid()){
                    mmag = tmpin2.k.commonContainer_mag((*ite).k);

                    if (ite2.isValid()){
                        tmag = tmpin2.k.commonContainer_mag((*ite2).k);
                    if (tmag < mmag) mmag = tmag;
                }
                }else{
                    if (ite2.isValid()){
                        mmag = tmpin2.k.commonContainer_mag((*ite2).k);
                    }else break;
                }

                HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
                if  (mmag > tmpin2.k.mag) {
                    tmpin2.k.setMagnitude(mmag);
                    if (tmpin.k == tmpin2.k) break;
                    else redo.push_back(tmpin2);
                } else break;


                }
            }

        }else{
            // single block!

            if (ite.isValid()){
                mmag = tmpin.k.commonContainer_mag((*ite).k);

                if (ite2.isValid()){
                    tmag = tmpin.k.commonContainer_mag((*ite2).k);
                if (tmag < mmag) mmag = tmag;
            }
            }else{
                if (ite2.isValid()){
                    mmag = tmpin.k.commonContainer_mag((*ite2).k);
                }else{redo.push_back(tmpin); continue;}
            }
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mmag);
            if  (mmag > tmpin.k.mag) {
                tmpin.k.setMagnitude(mmag);
            }
            redo.push_back(tmpin);
        }
    }

    ((RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->batchInit(redo);

}


LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::fillvoid_directional(int direction, NODES before){ // use connectivity to fill the complete volumes, can only propagate in the input direction
    // fill void using the elems on the top

    KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > cur = this->max;
    Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > f_out;
    bool isPure;
    KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tcur;

    KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES > cmpress_cur;


    unsigned int i,j;
    j=0;

    NODES majority[32 * DIMS];

    cmpress_cur.k.toMax();
    cmpress_cur.k.hyperpos.getOrderAndLeadHyper(i,j); printf("%i\t%i\n",i,j);


    for(i=0;i<100000000;i++){






        if ((i % 100000) == 0) printf("tic\n");
        tcur.k = cur.k;
        tcur.k.toLower(2);

        if (cur.k.hyperpos[2] > tcur.k.hyperpos[2]){
     //   this->intersection(cur.k.hyperpos,f_out);
     //   if (f_out.size() != 0) break;
        this->insert_partition(tcur.k, cur.d);
        }


        typename RBTofDoom< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > ,1u>::Iterator ite = this->find_first(cur);

        if (ite.isValid()) {


        if (((*ite).k != cur.k)&&((*ite).k.commonContainer_mag(cur.k) == (*ite).k.mag)){
        }else {
            --ite;
            if (!ite.isValid()) break;
        }

        cur = *ite;
        } else cur = this->max;


    }


}



    // get the separation of *pure* block, stores oriantation of contact in par_mag
LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::getContacts( Vector< KeyElem< HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &out, HyperCursor<STORETYPE,DIMS,LEAD> &where)const{




}

LFHTEMP unsigned int SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::getContainerOf(const HyperPosition<STORETYPE,DIMS,LEAD> posis) const{
	KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,NODES> query;  query.k = posis;

	unsigned int prev,next;
	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,NODES> >*)this)->findPN(query,prev,next);




	return 0;
	}





LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::intersection(const HyperPosition<STORETYPE,DIMS, LEAD> &pos, vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> > &f_out) const{

//	int i;
	Tuple<STORETYPE, DIMS> min[2];

	min[0][0] = 0; min[0][1] = 0; min[0][2] = 0;
	min[1][0] = 0; min[1][1] = 15; min[1][2] = 0;

	HyperPositionQuery_Cube<STORETYPE, DIMS, LEAD,  NODES > cq( min[0], min[1]);

	(( RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES> >const *const )this)->intersection(f_out, cq );



	}

LFHTEMP void SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::intersection(const HyperPosition<STORETYPE,DIMS, LEAD> &pos, Vector< KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > > &f_out) const{


	Tuple<STORETYPE, DIMS> min[2];

	min[0][0] = 0; min[0][1] = 0; min[0][2] = 0;
	min[1][0] = 0; min[1][1] = 15; min[1][2] = 0;

	HyperPositionQuery_Cube<STORETYPE, DIMS, LEAD, NODES > cq( min[0], min[1]);

	((RBTofDoom<KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>,  NODES > >const *const )this)->intersection(f_out, cq );



	}

	// SpacePartition private routines

LFHTEMP unsigned short SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::get_parent_order(const KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES>& from){
    // get closest Item
     KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > tmp_key(from.k.getMinBox(), from.d);
     typename RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key);
    unsigned short f_out;

	tmp_key.k = from.k.getMaxBox();
    if (ite.isValid()){
        if ((*ite) < tmp_key){
            // found an item within box, iterate up!
        //    if (*ite).
            while( (*ite).k.par_mag != ExCo<unsigned short>::mkMaximum() ){
                if ((*ite).k.par_mag > from.k.mag) break;
                tmp_key.k.setMagnitude((*ite).k.par_mag);
                ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key);
                if ((!ite.isValid())||((*ite).k.mag != tmp_key.k.mag)){
                    printf("Should never happen!\n");
                }
            }
            return( (*ite).k.par_mag );
        }

    }else{
        if (((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->size == 0) return(ExCo<unsigned short>::mkMaximum());
        ite = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->last();
    }

    // check knowing that no elements exists inside

    // the one pointed to, and the the precepding are to find the best quandidates if existing



    tmp_key.k = (*ite).k;
    while(!tmp_key.k.strictly_contains(from.k)){


    }

    typename RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2 = ((RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >*)this)->find_first(tmp_key);

    return(0);
}

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::first(){
        if (partition.size() == 0)return false;

        HyperCursor<STORETYPE,DIMS,LEAD> smallest; smallest.toMin();
        unsigned short mag = smallest.commonContainer_mag(partition.min.k);
        if (mag == partition.min.k.mag) {
            (*this).curkey = partition.min;
        }else{
            HyperCursor<STORETYPE,DIMS,LEAD>::decr_mag(mag);
            smallest.setMagnitude(mag);
            (*this).curkey = smallest;
            ExOp::toZero((*this).curkey.d);
        }
        return true;
    }

LFHTEMP bool SpacePartition<STORETYPE, DIMS, LEAD, NODES,0u>::KeyIterator::next(){
    typename RBTofDoom<  KeyElem<HyperCursor<STORETYPE,DIMS,LEAD>, NODES > >::Iterator ite2;

    if (ite2.first()){
    }else{

    }
return false;
    }

*/
 // end of namespace
