/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */


#undef LFHTEMP
#define LFHTEMP template <int M, int O, int H>

LFHTEMP MessageBuffer<M,O,H>::MessageBuffer(){
	((uint32_t*)buffer)[0] = sizeof(uint32_t) << 1;
	((uint32_t*)buffer)[1] = sizeof(uint32_t) << 1;
}

LFHTEMP void MessageBuffer<M,O,H>::printBufferStatus(FILE *f)const{
	if (((uint32_t*)buffer)[0] >= ((uint32_t*)buffer)[1]){
		fprintf(f, "(%i / %i bytes of buffer available)\n", (1 << M) + ((uint32_t*)buffer)[1] - ((uint32_t*)buffer)[0], (1 << M));
	}else{
		fprintf(f, "(%i / %i bytes of buffer available)\n", ((uint32_t*)buffer)[1] - ((uint32_t*)buffer)[0], (1 << M));
	}
}
/** \brief Tries to make space for new message
 *
 * \param length of message required
 * \param nbreads exact number of reads expected (message is deleted after)
 *
 * \return address of new message or NULL pointer
 */
LFHTEMP typename MessageBuffer<M,O,H>::NewMessageScope MessageBuffer<M,O,H>::tryAlloc(uint16_t length, uint8_t nbreads){ typename MessageBuffer<M,O,H>::NewMessageScope fout(*this);
	//printf("trying %p... state is: %i %i\n", (void*)this, nbreads, length);
	if ((nbreads == 0)||(length >= 2047)) fout.msg = NULL;
	else {
		if (((uint32_t*)buffer)[0] >= ((uint32_t*)buffer)[1]){
	//		printf("pre %i >= %i\n", ((uint32_t*)buffer)[0], ((uint32_t*)buffer)[1]);
			if ((((*(uint16_t*)buffer) + length + 4) >> M) != 0){
				if (((uint32_t*)buffer)[1] <= (sizeof(uint32_t) << 1) + 2 + length) {fout.msg = NULL;return fout;}
				// loop around time!!
				*(uint16_t*)(buffer + ((uint32_t*)buffer)[0])  = 0;
				fout.nextpos = (sizeof(uint32_t) << 1) ;
			}else fout.nextpos = ((uint32_t*)buffer)[0] ;
		}else{
	//		printf("post  %i < %i\n", ((uint32_t*)buffer)[0], ((uint32_t*)buffer)[1] );
			if (((uint32_t*)buffer)[0] + (length + 2) >= ((uint32_t*)buffer)[1]) {fout.msg = NULL;return fout;}
		    fout.nextpos = ((uint32_t*)buffer)[0];
		}
		*(uint16_t*)(buffer + fout.nextpos) = length | ((nbreads - 1) << 11);
		fout.msg = buffer +( fout.nextpos+ 2);
		fout.nextpos += length + 2;
	}
return fout;}


/** \brief Return the location of the latest message
 *
 * \param position current read position (will be updated)
 * \param length current message length
 * (undefined behavior warning) assumes that reader knows that a message is available!
 * \return address of message
 */
LFHTEMP uint8_t* MessageBuffer<M,O,H>::readAt(uint32_t& position, uint16_t& length){
	if (position == ((uint32_t*)buffer)[0]) return NULL;
	length = *((uint16_t*)(buffer + position));
	if (length == 0) { // wrap aroung!
		length = (*(uint16_t*)(buffer + (sizeof(uint32_t) << 1)));
		position = (sizeof(uint32_t) << 1);
	}
	if (length & 0xF800) {
			*((uint16_t*)(buffer + position)) -= 0x0800;
			length &= 0x07FF;
	}else ((uint32_t*)buffer)[1] = position + (length + 2);
	position += 2;
	uint8_t* fout = buffer + position;
	position += length;
return fout;}

/** \brief Notifies that message after will not be read by one source
 *
 * \param position current reader position (to flush
 *
 * \return address of message
 */
LFHTEMP void MessageBuffer<M,O,H>::flushReader(uint32_t position){
	while(position != ((uint32_t*)buffer)[0]){
		uint16_t length = *((uint16_t*)(buffer + position));
		if (length == 0){
			position = (sizeof(uint32_t) << 1);
			length = *((uint16_t*)(buffer + position));
		}
		if (length & 0xF800) {
			*((uint16_t*)(buffer + position)) -= 0x0800;
			position += ((uint32_t)2) + (length & 0x07FFF);
		}else position = ((uint32_t*)buffer)[1] = position + (length + 2);
	}
}


LFHTEMP MessageBuffer<M,O,H>::NewMessageScope::~NewMessageScope(){
	if (msg != NULL) ((uint32_t*)target.buffer)[0] = nextpos;
}

 //static double systemSolveSph(double* a, double* b, PolyVector &c){
/*	double argd[5];
	Polything arg[16];
	argd[0] = b[0] - a[0]; //
	argd[1] = b[1] - a[1]; //
	argd[2] = b[2] - a[2]; //
	argd[3] = b[3] - a[3]; //
	arg[4] = c.pos[0] - a[0];
	arg[5] = c.pos[1] - a[1];
	arg[6] = c.pos[2] - a[2];
	arg[7] = c.pos[3] + a[3];

	argd[4] = argd[0]*argd[0] + argd[1]*argd[1] + argd[2]*argd[2];

	arg[9] = arg[4]*argd[0] +arg[5]*argd[1] +arg[6]*argd[2]; // c-a * b-a

	arg[10] =  arg[5]*argd[2]- arg[6]*argd[1];
	arg[11] =  arg[6]*argd[0]- arg[4]*argd[2];
	arg[12] =  arg[4]*argd[1]- arg[5]*argd[0];
	arg[13] =  arg[10]* arg[10] + arg[11] * arg[11] + arg[12] * arg[12];
	arg[14] =  arg[3]* arg[3] - argd[4];
	arg[15] = arg[9] * arg[3] + arg[7]*argd[4];


	Polything f = (arg[15]*arg[15] + arg[14]*arg[13]) * -1.0f;*/

	//return(f.solve());
//	return(0.0f);
//}


/* Class Ressource
 *
 * stored header of ressources which mainly order use ressources by their
 * last time used
 *
 */

/*
void showDoubleMap(char* path, double* data, int resolution, int sizex, int sizey, int nbchannels){
	image im;
	double* maxs = new double[nbchannels];
	memset(maxs,'\0',sizeof(double)*nbchannels);
	int x,y,i,j;
	for(y=0;y<sizey;y++) for (x=0;x<sizey;x++){
		for(j=0;j<nbchannels;j++){
		if (data[j + nbchannels * (x + y * sizex)] > maxs[j]) maxs[j] = data[j + nbchannels * (x + y * sizex)];
	}
}

	for(j=0;j<nbchannels;j++) maxs[j] = 255.0f / maxs[j];
	im.data = new unsigned char[resolution*resolution*sizex*sizey*3];

	im.sizex = sizex * resolution;
	im.sizey = sizey * resolution;


	int w;
	int v;
	for(y=0;y<sizey;y++) for (x=0;x<sizey;x++){

		for(j=0;j<resolution;j++) for (i=0;i<resolution;i++){
			w= (i * nbchannels) / resolution;
			v =(int) data[w + nbchannels * (x + y * sizex)];
			im.data[  3*(i+resolution*(x+sizex*(j+resolution*y)))] = (unsigned char)((v * maxs[w]) * ((w & 1) ? 0: 1));
			im.data[1+3*(i+resolution*(x+sizex*(j+resolution*y)))] = (unsigned char)((v * maxs[w]) * ((w & 4) ? 0: 1));
			im.data[2+3*(i+resolution*(x+sizex*(j+resolution*y)))] = (unsigned char)((v * maxs[w]) * ((w & 2) ? 0: 1));
	}}
	im.save(path);
	delete[](maxs);
}*/
/*
template <class C>
Ptr<C>::Ptr():Id(0){}

template <class C>
Ptr<C>::~Ptr(){
	if (Id != 0) LR.decrLink(Id);
	}

template <class C>
C& Ptr<C>::operator->(){
	return(*((C*)LinkMem.getLink(Id)));
	}
template <class C>
C& Ptr<C>::operator*(){
	return(*((C*)LinkMem.getLink(Id)));
	}
template <class C>
Ptr<C>::operator C*(){
	return(((C*)LinkMem.getLink(Id)));
	}
template <class C>
C* Ptr<C>::operator=(C const * const nval){
	if (Id != 0){
		// decrement link
		LinkMem.decrLink(Id);
		}
	Id = LinkMem.incrgetID((void*)nval);
	}
template <class C>
C* Ptr<C>::operator=(Ptr<C> const & nval){
	if (Id != 0){
		// decrement link
		LinkMem.decrLink(Id);
		}
	Id = nval.Id;
	LinkMem.incrLink(Id);
	}*/







// class PolyThing


#undef LFHTEMP
#define LFHTEMP template <class C>

LFHTEMP ERRCODE FIFOQueue<C>::insert(Event<C>* ev){
    int pos = semaphore.fetch_add(1);
    if (((async_read - pos) & 0xFF) == 1){
        semaphore--;
        return 1; // maximum reached!
    }
    buffer[pos & 0xFF] = ev;
    return 0;
}

LFHTEMP void FIFOQueue<C>::runAll(){
    while(semaphore != async_read){
        (*buffer[async_read])(scope);
        delete(buffer[async_read]);
        async_read++;
        if (async_read & 0x100){
            async_read &= 0xFF;
            semaphore -= 0x100;
        }
    }
}

LFHTEMP void FIFOQueue<C>::runAll(C scope){
    while(semaphore != async_read){
        (*buffer[async_read])(scope);
        delete(buffer[async_read]);
        async_read++;
        if (async_read & 0x100){
            async_read &= 0xFF;
            semaphore -= 0x100;
        }
    }
}

#undef LFHTEMP
#define LFHTEMP template <class T, class O>

LFHTEMP ERRCODE AsyncInserter<T,O>::operator<<(O &what){
    int pos = semawrite.fetch_add(0x201); // reserving writing slot
    semawrite.fetch_and(0xFEFF);
    if (pos != semaread){
        if (((semaread - pos) & 0xFF) == 1){
            semawrite.fetch_add(-0x201);
            return 1; // maximum reached!
        }
        ExOp::toMemmove(buffer[pos & 0xFF],what);
        do{
            pos = semawrite.fetch_add(-0x200);
            if ((pos & 0xFE00) != 0x200) return 0;
            pos = semawrite.fetch_add(0x200);
        }while((pos & 0xFE00) != 0); // well, someone else just arrived...
    }else{ // writing to first unread element... insert instead!
        target << what;
    }
    semaread = (semaread +1) & 0xFF;
    while(true){
        pos = semawrite.fetch_add(-0x200);
        if (((pos & 0xFE00) != 0x200)||(pos == semaread + 0x200)) return 0; // DONE , or someone is busy writing... will take your place!)
        target << buffer[semaread & 0xFF];
        semawrite.fetch_add(0x200);
        semaread = (semaread +1) & 0xFF;
    }
    return 0;
}

LFHTEMP ERRCODE AsyncInserter<T,O>::operator<<=(O &what){
    int pos = semawrite.fetch_add(0x201); // reserving writing slot
    semawrite.fetch_and(0xFEFF);
    if (pos != semaread){
        if (((semaread - pos) & 0xFF) == 1){
            semawrite.fetch_add(-0x201);
            return 1; // maximum reached!
        }
        ExOp::toMemmove(buffer[pos & 0xFF],what);
        do{
            pos = semawrite.fetch_add(-0x200);
            if ((pos & 0xFE00) != 0x200) return 0;
            pos = semawrite.fetch_add(0x200);
        }while((pos & 0xFE00) != 0); // well, someone else just arrived...
    }else{ // writing to first unread element... insert instead!
        target <<= what;
    }
    semaread = (semaread +1) & 0xFF;
    while(true){
        pos = semawrite.fetch_add(-0x200);
        if (((pos & 0xFE00) != 0x200)||(pos == semaread + 0x200)) return 0; // DONE , or someone is busy writing... will take your place!)
        target <<= buffer[semaread & 0xFF];
        semawrite.fetch_add(0x200);
        semaread = (semaread +1) & 0xFF;
    }
    return 0;
}

LFHTEMP ERRCODE AsyncInserter<T,O>::insert(O &what){
    int pos = semawrite.fetch_add(0x201); // reserving writing slot
    semawrite.fetch_and(0xFEFF);
    if (pos != semaread){
        if (((semaread - pos) & 0xFF) == 1){
            semawrite.fetch_add(-0x201);
            return 1; // maximum reached!
        }
        ExOp::toMemmove(buffer[pos & 0xFF],what);
        do{
            pos = semawrite.fetch_add(-0x200);
            if ((pos & 0xFE00) != 0x200) return 0;
            pos = semawrite.fetch_add(0x200);
        }while((pos & 0xFE00) != 0); // well, someone else just arrived...
    }else{ // writing to first unread element... insert instead!
        target.insert(what);
    }
    semaread = (semaread +1) & 0xFF;
    while(true){
        pos = semawrite.fetch_add(-0x200);
        if (((pos & 0xFE00) != 0x200)||(pos == semaread + 0x200)) return 0; // DONE , or someone is busy writing... will take your place!)
        target.insert(buffer[semaread & 0xFF]);
        semawrite.fetch_add(0x200);
        semaread = (semaread +1) & 0xFF;
    }
    return 0;
}




#undef LFHTEMP
#define LFHTEMP template <class C, class A, int M>

LFHTEMP ERRCODE PriorityQueue<C,A,M>::insertSync(C priority, std::function<C(A)> ev){return async_buf.wrSync(OrderedTask<C, A>(priority, ev));}
LFHTEMP ERRCODE PriorityQueue<C,A,M>::insertAsync(C priority, std::function<C(A)> ev){return async_buf.wrAsync(OrderedTask<C, A>(priority, ev));}
LFHTEMP bool PriorityQueue<C,A,M>::run(A& scope){
    OrderedTask<C, A> tmp;
    while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
    if (heap.isEmpty()) return false;
    tmp = heap.pop();
    C fout = tmp.task(scope);
return true;}

LFHTEMP void PriorityQueue<C,A,M>::run_until(A& scope, C priority){
	OrderedTask<C, A> tmp;
	while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
	if (heap.isEmpty()) return;
	while( heap.top() < priority){
		tmp = heap.pop();
		C fout = tmp.task(scope);
	}
}

#undef LFHTEMP
#define LFHTEMP template <class C, int M>

LFHTEMP ERRCODE PriorityQueue<C,void,M>::insertSync(C priority, std::function<C()>  ev){return async_buf.wrSync(OrderedTask<C, void>(priority, ev));}
LFHTEMP ERRCODE PriorityQueue<C,void,M>::insertAsync(C priority, std::function<C()>  ev){return async_buf.wrAsync(OrderedTask<C, void>(priority, ev));}
LFHTEMP bool PriorityQueue<C,void,M>::run(){
    OrderedTask<C, void> tmp;
    while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
    if (heap.isEmpty()) return false;
    tmp = heap.pop();
    C fout = tmp.task();
return true;}

LFHTEMP void PriorityQueue<C,void,M>::run_until(C priority){
	OrderedTask<C, void> tmp;
	while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
	if (heap.isEmpty()) return;
	while( heap.top() < priority){
		tmp = heap.pop();
		C fout = tmp.task();
	}
}
#undef LFHTEMP
#define LFHTEMP template <int M>

LFHTEMP ERRCODE PriorityQueue<CircularInteger,void,M>::insertSync(int32_t time, std::function<int32_t(int32_t)>  ev){return async_buf.wrSync(OrderedTask<CircularInteger, int32_t,int32_t>(time, ev));}
LFHTEMP ERRCODE PriorityQueue<CircularInteger,void,M>::insertAsync(int32_t time, std::function<int32_t(int32_t)>  ev){return async_buf.wrAsync(OrderedTask<CircularInteger, int32_t,int32_t>(time, ev));}
LFHTEMP bool PriorityQueue<CircularInteger,void,M>::run_until(int32_t time){
    OrderedTask<CircularInteger, int32_t, int32_t> tmp;
    while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
    if (heap.isEmpty()) return false;
    tmp = heap.top();
	while (tmp.time < time) {
		CircularInteger fout = tmp.task(tmp.time);
		if (fout != 0) {tmp.time += fout; heap.pop_and_insert(tmp);}
		else heap.pop();
		while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
		if (heap.isEmpty()) return false;
		tmp = heap.top();
	}
return true;}
LFHTEMP void PriorityQueue<CircularInteger,void,M>::show(FILE*f, int level)const{
	fprintf(f, "Async PriorityQueue:\n");
    async_buf.show(f);
    heap.show(f);
}

#undef LFHTEMP
#define LFHTEMP template <class C, int M>

LFHTEMP ERRCODE PriorityQueue<CircularInteger,C,M>::insertSync(int32_t time, std::function<int32_t(C)>  ev){return async_buf.wrSync(OrderedTask<CircularInteger, C,int32_t>(CircularInteger(time), ev));}
LFHTEMP ERRCODE PriorityQueue<CircularInteger,C,M>::insertAsync(int32_t time, std::function<int32_t(C)>  ev){return async_buf.wrAsync(OrderedTask<CircularInteger, C,int32_t>(CircularInteger(time), ev));}
LFHTEMP bool PriorityQueue<CircularInteger,C,M>::run_until(C arg, int32_t time){
    OrderedTask<CircularInteger, C, int32_t> tmp;
    while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
    if (heap.isEmpty()) return false;
    tmp = heap.top();
	while (tmp.time < time) {
		CircularInteger fout = tmp.task(arg);
		if (fout != 0) {tmp.time += fout; heap.pop_and_insert(tmp);}
		else heap.pop();
		while(async_buf.rdSync(tmp)==0) heap.insert(tmp);
		if (heap.isEmpty()) return false;
		tmp = heap.top();
	}
return true;}
LFHTEMP void PriorityQueue<CircularInteger,C,M>::show(FILE*f, int level)const{
	fprintf(f, "Async PriorityQueue:\n");
    async_buf.show(f);
    heap.show(f);
}

#undef LFHTEMP
#define LFHTEMP template <class KEY>
LFHTEMP StochasticQueue<KEY>::StochasticQueue(int32_t start_time): time(start_time){
    for(int i = 0 ;i < 1; i++){
        lambdasum[i] = 0;
    }
    mainLambda = 0.0f;
}

LFHTEMP uint16_t StochasticQueue<KEY>::getIndex_from_exptime_routine(uint32_t expect_time) const{
    int exponent;
    int mantissa = (int)(log(frexp((double)expect_time, &exponent)) * 11.541560327111707258879397448015);

    exponent -= 7;
    if (exponent < 0) exponent = 0;
    mantissa = 0; // for now

    return (15-exponent) | (mantissa << 4);
}

LFHTEMP void StochasticQueue<KEY>::runAtIndex_routine(int index){
    uint32_t i; ExOp::toRand(i);
    uint32_t within = i % lambdasum[index];
    uint16_t catindex = 15 | (index << 4);
    for(i=15;i!=0;i--,catindex--){
        if (within < (events.getSizeWithinCategory(catindex) << i)) break;
        else within -= (events.getSizeWithinCategory(catindex) << i);
    }
    uint32_t ite = events.getFirstWithinCategory(catindex);
    ite += (within >> i);
    if ((within >> i) > events.getSizeWithinCategory(catindex)) printf("%i > %i !!%i!!\n",(within >> i), events.getSizeWithinCategory(catindex), catindex);
    else (*events.deref(ite))(events.deref_key(ite));
}

LFHTEMP void StochasticQueue<KEY>::runTo(int32_t n_time){
    uint32_t r;
    while(time - n_time < 0){
        if (lambdasum[0] == 0) {time = n_time; return;}
        ExOp::toRand(r);
        double mtime = (2097152.0f / lambdasum[0]) * (22.180709777918249901351427886662f - log((double)r));
        runAtIndex_routine(0);
        time += (int32_t)mtime;
    }
}
LFHTEMP void StochasticQueue<KEY>::insert(const KEY key, Event<KEY>* ev, uint32_t expect_time){
    uint16_t index = getIndex_from_exptime_routine(expect_time);
    events.insert(key, ev, index);
    lambdasum[(index >> 4)] += (1 << (index & 15));
}
LFHTEMP void StochasticQueue<KEY>::updateTime(const KEY &key, uint32_t new_expect_time){
    uint16_t index = events.getCategory(key);
    lambdasum[(index >> 4)] -= (1 << (index & 15));
    index = getIndex_from_exptime_routine(new_expect_time);
    lambdasum[(index >> 4)] += (1 << (index & 15));
    events.setCategory(key, index);
}
LFHTEMP void StochasticQueue<KEY>::remove(const KEY &key){
    uint16_t index = events.getCategory(key);
    lambdasum[(index >> 4)] -= (1 << (index & 15));
    events.erase(key);
}

LFHTEMP void StochasticQueue<KEY>::show(FILE* f, int level) const{
    fprintf(f,"Stochastic queue:\n");
  //  events.show(f,level);
    fprintf(f,"Lambda sum:\n");
    ExOp::show(lambdasum);
    uint16_t catindex = 0;
    uint32_t ite;
    int j;
    for(int i = 0;i<16;i++, catindex++){
        j = events.getSizeWithinCategory(catindex);
        ite = events.getFirstWithinCategory(catindex);
        fprintf(f,"%i[%i at %i]: ",i, j, ite);
        while(j != 0){
            ExOp::show(events.deref_key(ite),f,1);
            j--; ite++;
            fprintf(f,"\t");
        }
        fprintf(f,"\n");
    }
}

#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int SIZE, unsigned int ORDER>

LFHTEMP PolyTuple<C,SIZE,ORDER>::PolyTuple(Tuple<C,SIZE> const & input, double weight){
      //  for(unsigned int i = 0;i<ORDER;i++) ExOp::intpow(input[0],i);


    data[1] = input[0];
    unsigned int i,s,k;
    s = ORDER;
    ExOp::toOne(data[0]);

    if (SIZE == 1){

    for(i=1;i<s;i++) data[i] = ExOp::mkPowInt(input[0],i);
    }else{
        Tuple<unsigned int, SIZE> ite;
        ExOp::toZero(ite);
        Tuple<C, SIZE> tmp; ExOp::toOne(tmp);k=1;
    while(true){
        if (ite[0] != 0) {
        ite[1]++; ite[0]--;
        if (SIZE == 2) data[k++] = ExOp::mkPowInt(input[0],ite[0]) * ExOp::mkPowInt(input[1],ite[1]);
        else data[k++] = tmp[2] * ExOp::mkPowInt(input[0],ite[0]) * ExOp::mkPowInt(input[1],ite[1]);
    }else{
        unsigned int i;
        for(i=1;i<SIZE-1;i++) if (ite[i] != 0) break;
        if (i == SIZE-1) {
            ite[0] = ite[SIZE-1] +1;
            if (ite[0] == ORDER) break;
            ExOp::toOne(tmp);
            ite[SIZE-1] =0;s++;
            data[k++] = ExOp::mkPowInt(input[0],ite[0]);
        }else{
            ite[0] = ite[i] -1;
            ite[i] =0;
            ite[i+1]++;
            if (i == SIZE -2) {tmp[SIZE -1] = ExOp::mkPowInt(input[SIZE -1],ite[SIZE -1]); i--;}
            for(i++;i>0;i--) tmp[i] = tmp[i+1] * ExOp::mkPowInt(input[i],ite[i]);
            data[k++] = tmp[1] * ExOp::mkPowInt(input[0],ite[0]);
        }
    }

    }
    }

	data *= weight;
}


LFHTEMP C PolyTuple<C,SIZE,ORDER>::operator[](const Tuple<C, SIZE>& in) const{
    Tuple<C,SIZE> buffer;
    unsigned int i;
    if (ORDER == 1) return(data[0]);
    if (SIZE == 1){
        i = SIZE-1; buffer[0] = data[i];
        for(i--; i != 0xFFFFFFFF;i--) buffer[0] = buffer[0] * in[0] + data[i];
        return buffer[0];
    }else{
        i = 0;
    }
}

LFHTEMP C& PolyTuple<C,SIZE,ORDER>::cell(const Tuple<unsigned int, SIZE>& in){
    if (SIZE == 1) return data[in[0]];
}

LFHTEMP C PolyTuple<C,SIZE,ORDER>::cell(const Tuple<unsigned int, SIZE>& in) const{
    if (SIZE == 1) return data[in[0]];
}

LFHTEMP PolyTuple<C,SIZE,ORDER-1> PolyTuple<C,SIZE,ORDER>::getDerivative(unsigned int direction)const{
    PolyTuple<C,SIZE,ORDER-1> fout;
    if (ORDER == 1) fout.data[0] = data[direction+1];
    else{
       bool bite;
       typename PolyTuple<C,SIZE,ORDER>::KeyIterator rite(*this);
        unsigned int k,l;
       k=0;l=0;
       if ((rite.first())) do{
           while(rite()[direction] == 0) {if (!rite.next()) return fout; l++;}
           fout.data[k++] = (*this).data[l++] * rite()[direction];
       } while (rite.next());
    }
    return fout;
    }



 LFHTEMP PolyTuple<C,SIZE,ORDER>& PolyTuple<C,SIZE,ORDER>::shiftdir(C const & other, unsigned int d){
    Tuple<unsigned int, ORDER> offsets;
    unsigned int l,k;
    unsigned int coef;
    // init offset to subgroups
    offsets[0] =0;
    if (ORDER > 1){
    offsets[1] =1;
    for(l=2;l<ORDER;l++) offsets[l] = (offsets[l-1] * (SIZE +l-1)) / (l-1);
    }
    typename PolyTuple<C,SIZE,ORDER>::KeyIterator ite(*this);
    C tmp;

    if (ite.first_oforder(ORDER-1)) do{

        for(k = ORDER - ite()[d] -1; k< ORDER-1;k++){
            coef = 1;
            tmp = data[offsets[ORDER-1]];
            for(l= ORDER-2;l>k;l--){
                coef *= 1+l - k;
                coef /= ORDER - l -1;

                tmp = (tmp * other) + (data[offsets[l]] * coef);
            }
            data[offsets[k]] += tmp*other;
        }
        for(k = ite()[d]; k< ORDER;k--) offsets[ORDER-k-1]++;
    }while (ite.next());

    return (*this);
    }

 LFHTEMP PolyTuple<C,SIZE,ORDER>& PolyTuple<C,SIZE,ORDER>::statshiftdir(C const & other, unsigned int d){
    Tuple<unsigned int, ORDER> offsets;
    unsigned int l,k;
    unsigned int coef;
    // init offset to subgroups
    offsets[0] =0;
    if (ORDER > 1){
    offsets[1] =1;
    for(l=2;l<ORDER;l++) offsets[l] = (offsets[l-1] * (SIZE +l-1)) / (l-1);
    }
    typename PolyTuple<C,SIZE,ORDER>::KeyIterator ite(*this);
    C tmp;
    if (ite.first_oforder(ORDER-1)) do{


        for(k = ORDER-1 ; k > ORDER - ite()[d]-1;k--){
            coef= 1;
            tmp = data[offsets[ORDER-ite()[d]-1]];
            for(l= ORDER-ite()[d];l<k;l++){
              //  printf("%f pow %i is %f\n", -other, l - (ORDER-k-2),  ExOp::mkPowInt(-other,l - (ORDER-k-2)) * coef);
              //  printf("at %i += %i, coef %i\n%f", offsets[ORDER-k-2], offsets[l], coef, data[offsets[ORDER-k-2]]); fflush(stdout);
              //  data[offsets[ORDER-k-2]] -= data[offsets[l]] * (ExOp::mkPowInt(-other,l - (ORDER-k-2)) * coef);

             //   printf(" -= %f\n", data[offsets[l]] * (ExOp::mkPowInt(-other,l - (ORDER-k-2)) * coef)); fflush(stdout);
                coef *= 1+k - l;
                coef /= l-ORDER+ite()[d]+1;
                tmp = (tmp * other) + (data[offsets[l]] * coef);
            }
            data[offsets[k]] += tmp * other;
        }
        for(k = ite()[d]; k< ORDER;k--) offsets[ORDER-k-1]++;
    }while (ite.next());

    return (*this);
    }

LFHTEMP PolyTuple<C,SIZE,ORDER> PolyTuple<C,SIZE,ORDER>::operator*(DBL_Weight const & other)const{ PolyTuple<C,SIZE,ORDER> fout;
	unsigned int max = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans;
	for( unsigned int i=0;i<max;i++) fout.data[i] = ExOp::mkMult(data[i], other);
	return fout;
	}
LFHTEMP PolyTuple<C,SIZE,ORDER>& PolyTuple<C,SIZE,ORDER>::operator*=(DBL_Weight const & other){
	unsigned int max = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans;
	for( unsigned int i=0;i<max;i++) ExOp::toMult(data[i], other);
	return *this;
	}

LFHTEMP PolyTuple<C,SIZE,ORDER*2-1> PolyTuple<C,SIZE,ORDER>::operator*(PolyTuple<C,SIZE,ORDER> const & other)const{
        PolyTuple<C,SIZE,ORDER*2-1> fout; fout.toZero();
        for(unsigned int i = TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i!=0xFFFFFFFF;i--) data[i] -= other.data[i];  return(*this);
        return fout;
    }


LFHTEMP template<unsigned int OSIZE> PolyTuple<C,OSIZE,ORDER>  PolyTuple<C,SIZE,ORDER>::operator<<(TMatrix<C,OSIZE,SIZE> const & xform) const{PolyTuple<C,OSIZE,ORDER> fout; ExOp::toZero(fout);
	typename PolyTuple<C,OSIZE,ORDER>::KeyIterator oite(fout);
	typename PolyTuple<C,SIZE,ORDER>::KeyIterator ite(*this);
	bool iiout;
	oite.first();
	fout.data[0] = this->data[0];
	unsigned int l,k;
	l=0;
	while(oite.next()){l++;
		ite.first(); k=1;
		while(ite.next()){k++;
			 while(oite.s != ite.s) {k++;if (!ite.next()) {oite.s =0;break;}}
			 if (oite.s == 0) break;
			 fout.data[l] += this->data[k];
		}
	}

	return fout;
}


LFHTEMP void PolyTuple<C,SIZE,ORDER>::show(FILE* f, int level) const{
        unsigned int i;
        switch(level){
        case 0:
        for(i = 0;i!=TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans;i++) {ExOp::show(data[i],f,level+1); fprintf(f,"\n");}
        break;
        case 1:
        for(i = 0;i!=TEMPLATE_TRIANGLE_NUMBER<ORDER,SIZE>::ans-1;i++) {ExOp::show(data[i],f,level+1); fprintf(f,"\n");}
        ExOp::show(data[i],f,level+1);
        break;
        }
        }


LFHTEMP bool PolyTuple<C,SIZE,ORDER>::KeyIterator::first(){ExOp::toZero((*this).curkey);s=0;return(true);}
LFHTEMP bool PolyTuple<C,SIZE,ORDER>::KeyIterator::next(){
    if (SIZE == 1){
        (*this).curkey[0]++;
        return (*this).curkey[0] != ORDER-1;
    }else if ((*this).curkey[0] != 0) {
        (*this).curkey[1]++; (*this).curkey[0]--;
    }else{
        unsigned int i;
        for(i=1;i<SIZE-1;i++) if ((*this).curkey[i] != 0) break;
        if (i == SIZE-1) {
            (*this).curkey[0] = (*this).curkey[SIZE-1] +1;
            if ((*this).curkey[0] == ORDER) return false;
            (*this).curkey[SIZE-1] =0;s++;
        }else{
            (*this).curkey[0] = (*this).curkey[i] -1;
            (*this).curkey[i] =0;
            (*this).curkey[i+1]++;
        }
    }
    return true;
    }

LFHTEMP bool PolyTuple<C,SIZE,ORDER>::KeyIterator::last(){
    unsigned int i=0;
    for(i=0;i<SIZE-1;i++) (*this).curkey[i]=0;
    (*this).curkey[i] = ORDER-1; s = ORDER-1;
    return(true);
    }
LFHTEMP bool PolyTuple<C,SIZE,ORDER>::KeyIterator::prev(){
    if (SIZE == 1){
        (*this).curkey[0]--;
        return (*this).curkey[0] != 0xFFFFFFFF;
    }else if ((*this).curkey[1] != 0) {
        (*this).curkey[1]--; (*this).curkey[0]++;
    }else{
        unsigned int i;
        for(i=2;i<SIZE;i++) if ((*this).curkey[i] != 0) break;
        if (i == SIZE) {
            (*this).curkey[SIZE-1] = (*this).curkey[0];
            if (((*this).curkey[SIZE-1]--) == 0) return false;
            (*this).curkey[0] = 0;s--;
        }else{
            (*this).curkey[i]--;
            (*this).curkey[i-1] = 1 + (*this).curkey[0];
            (*this).curkey[0] = 0;
        }
    }
    return true;
    }

LFHTEMP bool PolyTuple<C,SIZE,ORDER>::KeyIterator::first_oforder(unsigned int order){
    if (order >= ORDER) return false;
    ExOp::toZero((*this).curkey);s=order;
    (*this).curkey[0] = order;
    return(true);
}

// class PolyThing
#undef LFHTEMP
#define LFHTEMP template <class C, unsigned int SIZE>

LFHTEMP typename MT_IFTYPE<SIZE-2, PolyThing<C,SIZE-1> ,C>::TYPE PolyThing<C,SIZE>::mkDerivative() const{typename MT_IFTYPE<SIZE-2, PolyThing<C,SIZE-1> ,C>::TYPE fout;
	if (SIZE == 2) fout = data[1];
	else for(unsigned int i=0; i<SIZE-1;i++) fout.data[i] = data[i+1];
	return fout;
}

// class PolyThing
#undef LFHTEMP
#define LFHTEMP template <class C>

LFHTEMP C PolyThing<C,0u>::eval(double where){
	C _out = coef[order-1];
	int i = order -2;
	for(;i>=0;i--) _out = coef[i] + where * _out;
	//C &__out = _out;
return(_out);}

LFHTEMP PolyThing<C,0u> PolyThing<C,0u>::secondDerivative(){
	PolyThing<C,0u> _out;
	_out.setOrder(order-3);
	int i;
	for(i=0;i<order-2;i++) _out.coef[i] = coef[i+2] * (i*i+3*i+2);
	return(_out);
	}

LFHTEMP double PolyThing<C,0u>::newtonRoot(double guess, double min, double max) const{
	C moms[2];
	int i;
	bool err =false;
	double tmp;
	double val[2];
	val[0] = ExCo<double>::max();
	val[1] = ExCo<double>::min();
	int j;
	for(j=0;j<100;j++){
		i =order;
		moms[0] = coef[i];
		moms[1] = coef[i] * (double)i;
		for(i--;i>2;i--){
			moms[0] = (coef[i]) + (moms[0] * guess);
			moms[1] = (coef[i] * (double)i) + (moms[1] * guess);
		}
		val[ (moms[0] < ExCo<C>::zero()) ? 0 : 1] = guess;
		if (moms[1] != ExCo<C>::zero()){
		tmp = guess - ExCo<C>::doubleratio(moms[0],moms[1]);
			if (tmp == guess) break;
		}else err = true;

		if ((tmp < min)||(tmp >max)||(err)){
			err =0;
			if (val[0] < val[1]){
				guess = (val[0] + val[1])/2;
			}else{
				switch((rand()) & 15){
					case 0: tmp = -guess;
					case 1: tmp = ExCo<C>::invert(guess);
					default: tmp *= exp( ((double)(rand() - (RAND_MAX >> 1))) /(RAND_MAX >> 4));
				}
				if (tmp < min) guess = (tmp + min) *0.5f;
				else if (tmp > max) guess = (tmp + max) *0.5f;
				else guess = tmp;
			}
		}
	}
	return(guess);
}


LFHTEMP	double PolyThing<C,0u>::findMax(double imin, double imax){
		PolyThing<C,0u> nd = secondDerivative();
		Vector<mycomplex> zeros = nd.getZeros();

		double best,bestv;
		double tmp,tmpv;
		int i;
		best =imin;
		bestv = eval(imin);

		if (imin != imax){
			tmpv = eval(imax);
			if (tmpv > bestv){
				tmpv = bestv;
				best = imax;
			}
		}
		for(i=zeros.getSize()-1;i>=0;i--){
			if (zeros[i][1] == 0.0f) {
				tmpv = eval(zeros[i][0]);
				if (tmpv > bestv){
					tmpv = bestv;
					best = zeros[i][0];
				}
			}
		}
		return(best);
	}

LFHTEMP Vector<mycomplex> PolyThing<C,0u>::getZeros(){
	double tmp, tmp2;
	Vector<mycomplex> _out;
	switch(order){
		case 2:
			tmp = coef[1] / (coef[2] *2);
			tmp2 = tmp*tmp +  coef[0] / coef[2];
			if (tmp2 < 0) {
				tmp2 = sqrt(-tmp2);
				_out.push_back(mycomplex(tmp, tmp2));
				_out.push_back(mycomplex(tmp, -tmp2));
			}else {
				tmp2 = sqrt(tmp2);
				_out.push_back(mycomplex(tmp + tmp2,0.0f));
				_out.push_back(mycomplex(tmp - tmp2,0.0f));
			}
		break;

	}
	return(_out);
	}

LFHTEMP PolyThing<C,0u> PolyThing<C,0u>::interpolate(C* pt, int nbpts){

		PolyThing<C,0u> _out;
		_out.setOrder(nbpts-1);
		switch(nbpts){
			case 1: _out.coef[0] = pt[0]; break;
			case 2: _out.coef[0] = (pt[0] + pt[1]) *0.5f; _out.coef[1] = (pt[1] - pt[0]) *0.5f;  break;
			case 3: _out.coef[0] = pt[1]; _out.coef[1] = (pt[2] - pt[0]) *0.25f; _out.coef[2] = (pt[2] + pt[0]) *0.125f - pt[1]*0.25f;  break;
			case 4:{// DONE
				C tmp[2];
				tmp[0] = (pt[3] - pt[0])/ (3 * 16.0f);
				tmp[1] = (pt[2] - pt[1])/16.0f;
			 _out.coef[1] = -tmp[0] + 9.0f *tmp[1];
			 _out.coef[3] = tmp[0] - tmp[1];
			 tmp[0] = (pt[3] + pt[0])/ 16.0f;
			 tmp[1] = (pt[2] + pt[1]) /(16.0f);
			 _out.coef[0] = -tmp[0] + 9.0f *tmp[1];
			 _out.coef[2] = tmp[0] - tmp[1];
		 }
			break;
			case 5:{
				C tmp[2];
				tmp[0] = pt[4] - pt[0];
				tmp[1] = pt[3] - pt[1];
			 _out.coef[1] = -tmp[0] + tmp[2];
			 _out.coef[3] = tmp[0] - tmp[2];
			 tmp[0] = pt[4] + pt[0];
			 tmp[1] = pt[3] + pt[1];
			 _out.coef[0] = pt[2];
			 _out.coef[2] = -tmp[0] + tmp[2] - pt[2];
			 _out.coef[4] = tmp[0] - tmp[2] + pt[2];
		 }
		 break;
			case 6:{// DONE
				C tmp[3];
			  tmp[0] = (pt[5] - pt[0]) / (5.0f * 768.0f);
			  tmp[1] = (pt[4] - pt[1]) / (-6.0f * 128.0f);
			  tmp[2] = (pt[3] - pt[2]) / 384.0f;
			 _out.coef[1] = 9*tmp[0] +25*tmp[1] + 225*tmp[2];
			 _out.coef[3] = -10*tmp[0] -26*tmp[1] -34*tmp[2];
			 _out.coef[5] = tmp[0] +tmp[1] + tmp[2];
			 tmp[0] = (pt[5] + pt[0]) / (768.0f);
			 tmp[1] = (pt[4] + pt[1]) / (-2.0f * 128.0f);
			 tmp[2] = (pt[3] + pt[2])/ (384.0f);
			 _out.coef[0] = 9*tmp[0] +25*tmp[1] +225*tmp[2] ;
			 _out.coef[2] = -10*tmp[0] -26*tmp[1] -34*tmp[2];
			 _out.coef[4] = tmp[0] + tmp[1] +tmp[2];
		 		}



			 break;
			case 7:{
				C tmp[3];
			 tmp[0] = (pt[6] - pt[0]);
			 tmp[1] = (pt[5] - pt[1]);
			 tmp[2] = (pt[4] - pt[2]);
			 _out.coef[1] = (tmp[0] -9 * tmp[1] + 45.0f * tmp[2])/60.0f;
			 _out.coef[3] = (-15.0f *tmp[0] +120.0f * tmp[1] - 195.0f * tmp[2])/720.0f;
			 _out.coef[5] = (3.0f *tmp[0] -12.0f * tmp[1] + 15.0f * tmp[2])/720.0f;
			 tmp[0] = (pt[6] + pt[0]);
			 tmp[1] = (pt[5] + pt[1]);
			 tmp[2] = (pt[4] + pt[2]);
			 _out.coef[0] = pt[3];
			 _out.coef[2] = (4.0f *tmp[0] -108.0f *  tmp[1] + 540.0f * tmp[2] - 980.0f * pt[3] )/720.0f;
			 _out.coef[4] = (-5.0f *tmp[0] +60.0f *  tmp[1] - 195.0f * tmp[2] + 280.0f * pt[3])/720.0f;
			 _out.coef[6] = (1.0f *tmp[0] -6.0f *  tmp[1] + 15.0f * tmp[2] - 20.0f * pt[3])/720.0f;
		 }
				break;
			case 8:{ // Done!
				C tmp[4];
			  tmp[0] = (pt[7] - pt[0]) / (-7*92160);
			  tmp[1] = (pt[6] - pt[1]) / (5*18432);
			  tmp[2] = (pt[5] - pt[2]) / (-3*10240.0f);
			  tmp[3] = (pt[4] - pt[3]) / (9216.0f*2);
			 _out.coef[1] = 225*tmp[0]+441*tmp[1]+1225.0f*tmp[2]+ 11025.0f * tmp[3];
			 _out.coef[3] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2]+ -1891.0f * tmp[3];
			 _out.coef[5] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ 83.0f*tmp[3];
			 _out.coef[7] = -tmp[0]-tmp[1]-tmp[2]-tmp[3];
			 tmp[0] = (pt[7] + pt[0])/ (-92160);
			 tmp[1] = (pt[6] + pt[1])/ (18432);
			 tmp[2] = (pt[5] + pt[2])/ (-10240.0f);
			 tmp[3] = (pt[4] + pt[3])/ (9216.0f *2);
			 _out.coef[0] = 225*tmp[0]+441*tmp[1]+1225.0f*tmp[2]+11025.0f * tmp[3];
			 _out.coef[2] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2]+-1891.0f * tmp[3];
			 _out.coef[4] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ 83.0f* tmp[3];
			 _out.coef[6] = -tmp[0]-tmp[1]-tmp[2]-tmp[3];
		 		}






			 break;
			 case 10:{ // done!
				C tmp[5];
			  tmp[0] = (pt[9] - pt[0]) / (9*20643840);
			  tmp[1] = (pt[8] - pt[1]) / (-7*2949120);
			  tmp[2] = (pt[7] - pt[2]) / (5*1032192);
			  tmp[3] = (pt[6] - pt[3]) / (-3*737280.0f);
			  tmp[4] = (pt[5] - pt[4]) / 1474560.0f;
			 _out.coef[1] = 11025*tmp[0]+18225*tmp[1]+35721*tmp[2]+99225.0f * tmp[3]+ 893025.0f * tmp[4];
			 _out.coef[3] = -12916*tmp[0]-21204*tmp[1]-40860*tmp[2]-106444.0f * tmp[3]-164196* tmp[4];
			 _out.coef[5] = 1974*tmp[0]+3094*tmp[1]+5278*tmp[2]+7374.0f*tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[7] = -84*tmp[0]-116*tmp[1]-140*tmp[2]-156.0f*tmp[3]-164* tmp[4];
			 _out.coef[9] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4];
			 tmp[0] = (pt[9] + pt[0]) / (20643840);
			 tmp[1] = (pt[8] + pt[1])/ (-2949120);
			 tmp[2] = (pt[7] + pt[2]) / (1032192);
			 tmp[3] = (pt[6] + pt[3])  / (-737280.0f);
			 tmp[4] = (pt[5] + pt[4]) / 1474560.0f;
			 _out.coef[0] = 11025*tmp[0]+18225*tmp[1]+35721*tmp[2]+99225.0f * tmp[3] + 893025.0f * tmp[4];
			 _out.coef[2] = -12916*tmp[0]-21204*tmp[1]-40860*tmp[2]-106444.0f * tmp[3] -164196* tmp[4];
			 _out.coef[4] = 1974*tmp[0]+3094*tmp[1]+5278*tmp[2]+7374.0f* tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[6] = -84*tmp[0]-116*tmp[1]-140*tmp[2]-156.0f*tmp[3]-164* tmp[4];
			 _out.coef[8] = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4];
		 		}


			 break;
		}
		return(_out);
	}

LFHTEMP	PolyThing<C,0u> PolyThing<C,0u>::interpolate_1free(C* pt, int nbpts){

		PolyThing<C,0u> _out;
		_out.setOrder(nbpts-2);
		switch(nbpts){
			case 1: exit(1); break;
			case 2: _out.coef[0] = (pt[0] + pt[1]) *0.5f;  break;
			case 3: _out.coef[0] = (pt[2] + pt[1] + pt[0]) /3.0f; _out.coef[1] = (pt[2] - pt[0]) *0.25f; break;
			case 4:{// DONE
				C tmp[2];
				tmp[0] = pt[3] - pt[0];
				tmp[1] = pt[2] - pt[1];
			 _out.coef[1] =  (3.0f*tmp[0] + tmp[1])/20.0f; // 82/10 coef[3]
			 tmp[0] = (pt[3] + pt[0])/ 16.0f;
			 tmp[1] = (pt[2] + pt[1]) /(16.0f);
			 _out.coef[0] = -tmp[0] + 9.0f *tmp[1];
			 _out.coef[2] = tmp[0] - tmp[1];
		 }
			break;
			case 5:{
				C tmp[2];
				tmp[0] = pt[4] - pt[0];
				tmp[1] = pt[3] - pt[1];
			 _out.coef[1] = -tmp[0] + tmp[2];
			 _out.coef[3] = tmp[0] - tmp[2];
			 tmp[0] = pt[4] + pt[0];
			 tmp[1] = pt[3] + pt[1];
			 _out.coef[0] = pt[2];
			 _out.coef[2] = -tmp[0] + tmp[2] - pt[2];
			 _out.coef[4] = tmp[0] - tmp[2] + pt[2];
		 }
		 break;
			case 6:{// DONE
				C tmp[3];
			  tmp[0] = (pt[5] - pt[0]) / (10.0f);
			  tmp[1] = (pt[4] - pt[1]) / -2.0f;
			  tmp[2] = (pt[3] - pt[2]);

			 _out.coef[1] = (-1375*tmp[0]	-1249*tmp[1] + 326*tmp[2])/3024.0f; // -11567 / 63 coef[5]
			 _out.coef[3] = (175*tmp[0] +49*tmp[1] -14*tmp[2])/3024.0f; // +290/9 coef[5]
//			 _out.coef[5] = ;
			 tmp[0] = (pt[5] + pt[0]) / (768.0f);
			 tmp[1] = (pt[4] + pt[1]) / (-2.0f * 128.0f);
			 tmp[2] = (pt[3] + pt[2])/ (384.0f);
			 _out.coef[0] = 9*tmp[0] +25*tmp[1] +225*tmp[2] ;
			 _out.coef[2] = -10*tmp[0] -26*tmp[1] -34*tmp[2];
			 _out.coef[4] = tmp[0] + tmp[1] +tmp[2];
		 		}


			 break;
			case 7:{
				C tmp[3];
			 tmp[0] = (pt[6] - pt[0]);
			 tmp[1] = (pt[5] - pt[1]);
			 tmp[2] = (pt[4] - pt[2]);
			 _out.coef[1] = (tmp[0] -9 * tmp[1] + 45.0f * tmp[2])/60.0f;
			 _out.coef[3] = (-15.0f *tmp[0] +120.0f * tmp[1] - 195.0f * tmp[2])/720.0f;
			 _out.coef[5] = (3.0f *tmp[0] -12.0f * tmp[1] + 15.0f * tmp[2])/720.0f;
			 tmp[0] = (pt[6] + pt[0]);
			 tmp[1] = (pt[5] + pt[1]);
			 tmp[2] = (pt[4] + pt[2]);
			 _out.coef[0] = pt[3];
			 _out.coef[2] = (4.0f *tmp[0] -108.0f *  tmp[1] + 540.0f * tmp[2] - 980.0f * pt[3] )/720.0f;
			 _out.coef[4] = (-5.0f *tmp[0] +60.0f *  tmp[1] - 195.0f * tmp[2] + 280.0f * pt[3])/720.0f;
			 _out.coef[6] = (1.0f *tmp[0] -6.0f *  tmp[1] + 15.0f * tmp[2] - 20.0f * pt[3])/720.0f;
		 }
				break;
			case 8:{
				C tmp[4];
			  tmp[0] = (pt[7] - pt[0]);
			  tmp[1] = (pt[6] - pt[1]);
			  tmp[2] = (pt[5] - pt[2]) / (3*10240.0f);
			  tmp[3] = (pt[4] - pt[3]) / (9216.0f*2);
			 _out.coef[1] = 1225.0f*tmp[2]+ 11025.0f * tmp[3];
			 _out.coef[3] = -1299.0f*tmp[2]+ -1891.0f * tmp[3];
			 _out.coef[5] = 75.0f*tmp[2]+ 83.0f*tmp[3];
			 _out.coef[7] = -tmp[2]-tmp[3];
			 tmp[0] = (pt[7] + pt[0]);
			 tmp[1] = (pt[6] + pt[1]);
			 tmp[2] = (pt[5] + pt[2])/ 10240.0f;
			 tmp[3] = (pt[4] + pt[3])/ (9216.0f *2);
			 _out.coef[0] = 1225.0f*tmp[2]+11025.0f * tmp[3];
			 _out.coef[2] = -1299.0f*tmp[2]+-1891.0f * tmp[3];
			 _out.coef[4] = 75.0f*tmp[2]+ 83.0f* tmp[3];
			 _out.coef[6] = -tmp[2]-tmp[3];
		 		}





			 break;
			 case 10:{
				C tmp[5];
			  tmp[0] = (pt[9] - pt[0]);
			  tmp[1] = (pt[8] - pt[1]);
			  tmp[2] = (pt[7] - pt[2]);
			  tmp[3] = (pt[6] - pt[3]) / (3*737280.0f);
			  tmp[4] = (pt[5] - pt[4]) / 1474560.0f;
			 _out.coef[1] = 99225.0f * tmp[3]+ 893025.0f * tmp[4];
			 _out.coef[3] = -106444.0f * tmp[3]-164196* tmp[4];
			 _out.coef[5] = 7374.0f*tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[7] = -156.0f*tmp[3]-164* tmp[4];
			 _out.coef[9] = tmp[3]+tmp[4];
			 tmp[0] = (pt[9] + pt[0]);
			 tmp[1] = (pt[8] + pt[1]);
			 tmp[2] = (pt[7] + pt[2]);
			 tmp[3] = (pt[6] + pt[3])  / 737280.0f;
			 tmp[4] = (pt[5] + pt[4]) / 1474560.0f ;
			 _out.coef[0] = 99225.0f * tmp[3] + 893025.0f * tmp[4];
			 _out.coef[2] = -106444.0f * tmp[3] -164196* tmp[4];
			 _out.coef[4] = 7374.0f* tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[6] = -156.0f*tmp[3]-164* tmp[4];
			 _out.coef[8] = tmp[3]+tmp[4];
		 		}

			 break;
		}
		return(_out);
	}

LFHTEMP	PolyThing<C,0u> PolyThing<C,0u>::interpolate_2free(C* pt, int nbpts){

		PolyThing<C,0u> _out;
		_out.setOrder(nbpts-3);
		switch(nbpts){
			case 1:
			case 2: exit(1); break;
			case 3: _out.coef[0] = (pt[2] + pt[1] + pt[0]) /3.0f; break;
			case 4:{// DONE
				C tmp[2];
				tmp[0] = pt[3] - pt[0];
				tmp[1] = pt[2] - pt[1];
			 _out.coef[1] = ( 3.0f*tmp[0] + tmp[1])/20.0f; // 82/10 coef[3]
			 _out.coef[0] = (pt[3] + pt[0] + pt[2] + pt[1])/(4.0f); // 5 coef[2]
		 }
			break;
			case 5:{
				C tmp[2];
				tmp[0] = pt[4] - pt[0];
				tmp[1] = pt[3] - pt[1];
			 _out.coef[1] = -tmp[0] + tmp[2];
			 _out.coef[3] = tmp[0] - tmp[2];
			 tmp[0] = pt[4] + pt[0];
			 tmp[1] = pt[3] + pt[1];
			 _out.coef[0] = pt[2];
			 _out.coef[2] = -tmp[0] + tmp[2] - pt[2];
			 _out.coef[4] = tmp[0] - tmp[2] + pt[2];
		 }
		 break;
			case 6:{ // DONE
				C tmp[3];
			  tmp[0] = (pt[5] - pt[0]) / (10.0f);
			  tmp[1] = (pt[4] - pt[1]) / -2.0f;
			  tmp[2] = (pt[3] - pt[2]);

			 _out.coef[1] = (-1375*tmp[0]	-1249*tmp[1] + 326*tmp[2])/3024.0f; // -11567 / 63 coef[5]
			 _out.coef[3] = (175*tmp[0] +49*tmp[1] -14*tmp[2])/3024.0f; // +290/9 coef[5]
//			 _out.coef[5] = ;
			 tmp[0] = (pt[5] + pt[0]) / (768.0f);
			 tmp[1] = (pt[4] + pt[1]) / (-2.0f * 128.0f);
			 tmp[2] = (pt[3] + pt[2])/ (384.0f);
			 _out.coef[0] = -72*tmp[0] -56*tmp[1] +144*tmp[2];// 190/7
			 _out.coef[2] = (120*tmp[0] +8*tmp[1] -48*tmp[2])/7;  // -81
			 //_out.coef[4] = tmp[0] + tmp[1] +tmp[2];
		 		}


			 break;
			case 7:{
				C tmp[3];
			 tmp[0] = (pt[6] - pt[0]);
			 tmp[1] = (pt[5] - pt[1]);
			 tmp[2] = (pt[4] - pt[2]);
			 _out.coef[1] = (tmp[0] -9 * tmp[1] + 45.0f * tmp[2])/60.0f;
			 _out.coef[3] = (-15.0f *tmp[0] +120.0f * tmp[1] - 195.0f * tmp[2])/720.0f;
			 _out.coef[5] = (3.0f *tmp[0] -12.0f * tmp[1] + 15.0f * tmp[2])/720.0f;
			 tmp[0] = (pt[6] + pt[0]);
			 tmp[1] = (pt[5] + pt[1]);
			 tmp[2] = (pt[4] + pt[2]);
			 _out.coef[0] = pt[3];
			 _out.coef[2] = (4.0f *tmp[0] -108.0f *  tmp[1] + 540.0f * tmp[2] - 980.0f * pt[3] )/720.0f;
			 _out.coef[4] = (-5.0f *tmp[0] +60.0f *  tmp[1] - 195.0f * tmp[2] + 280.0f * pt[3])/720.0f;
			 _out.coef[6] = (1.0f *tmp[0] -6.0f *  tmp[1] + 15.0f * tmp[2] - 20.0f * pt[3])/720.0f;
		 }
				break;
			case 8:{
				C tmp[4];
			  tmp[0] = (pt[7] - pt[0]);
			  tmp[1] = (pt[6] - pt[1]);
			  tmp[2] = (pt[5] - pt[2]) / (3*10240.0f);
			  tmp[3] = (pt[4] - pt[3]) / (9216.0f*2);
			 _out.coef[1] = 1225.0f*tmp[2]+ 11025.0f * tmp[3];
			 _out.coef[3] = -1299.0f*tmp[2]+ -1891.0f * tmp[3];
			 _out.coef[5] = 75.0f*tmp[2]+ 83.0f*tmp[3];
			 _out.coef[7] = -tmp[2]-tmp[3];
			 tmp[0] = (pt[7] + pt[0]);
			 tmp[1] = (pt[6] + pt[1]);
			 tmp[2] = (pt[5] + pt[2])/ 10240.0f;
			 tmp[3] = (pt[4] + pt[3])/ (9216.0f *2);
			 _out.coef[0] = 1225.0f*tmp[2]+11025.0f * tmp[3];
			 _out.coef[2] = -1299.0f*tmp[2]+-1891.0f * tmp[3];
			 _out.coef[4] = 75.0f*tmp[2]+ 83.0f* tmp[3];
			 _out.coef[6] = -tmp[2]-tmp[3];
		 		}





			 break;
			 case 10:{
				C tmp[5];
			  tmp[0] = (pt[9] - pt[0]);
			  tmp[1] = (pt[8] - pt[1]);
			  tmp[2] = (pt[7] - pt[2]);
			  tmp[3] = (pt[6] - pt[3]) / (3*737280.0f);
			  tmp[4] = (pt[5] - pt[4]) / 1474560.0f;
			 _out.coef[1] = 99225.0f * tmp[3]+ 893025.0f * tmp[4];
			 _out.coef[3] = -106444.0f * tmp[3]-164196* tmp[4];
			 _out.coef[5] = 7374.0f*tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[7] = -156.0f*tmp[3]-164* tmp[4];
			 _out.coef[9] = tmp[3]+tmp[4];
			 tmp[0] = (pt[9] + pt[0]);
			 tmp[1] = (pt[8] + pt[1]);
			 tmp[2] = (pt[7] + pt[2]);
			 tmp[3] = (pt[6] + pt[3])  / 737280.0f;
			 tmp[4] = (pt[5] + pt[4]) / 1474560.0f ;
			 _out.coef[0] = 99225.0f * tmp[3] + 893025.0f * tmp[4];
			 _out.coef[2] = -106444.0f * tmp[3] -164196* tmp[4];
			 _out.coef[4] = 7374.0f* tmp[3]+ 8614.0f * tmp[4];
			 _out.coef[6] = -156.0f*tmp[3]-164* tmp[4];
			 _out.coef[8] = tmp[3]+tmp[4];
		 		}

			 break;
		}
		return(_out);
	}

LFHTEMP	PolyThing<C,0u> PolyThing<C,0u>::interpolate_1free_middlefit(C* pt, int nbpts){

		PolyThing<C,0u> _out;
		_out.setOrder(nbpts-2);
		switch(nbpts){
			case 1:  exit(1);
			case 2: _out.coef[0] = (pt[0] + pt[1]) * 0.5f; break; // is ILL defined!
			case 3: _out.coef[0] = pt[1]; _out.coef[1] = (pt[2] - pt[0]) *0.25f; break;
			case 4:{// DONE
				C tmp[2];
//				tmp[0] = (pt[3] - pt[0])/ (-3.0f * 16.0f);
				//tmp[1] = (pt[2] - pt[1])/ 2.0f;
			 _out.coef[1] = (pt[2] - pt[1])/ 2.0f;
			 tmp[0] = (pt[3] + pt[0]) /(-16.0f);
			 tmp[1] = (pt[2] + pt[1]) /(16.0f);
			 _out.coef[0] = tmp[0] + 9.0f *tmp[1];
			 _out.coef[2] = -tmp[0] - tmp[1];
		 }
			break;
			case 5:{

		 }
		 break;
			case 6:{// DONE
				C tmp[3];
			  tmp[0] = (pt[5] - pt[0]) / (5.0f * 768.0f);
			  tmp[1] = (pt[4] - pt[1]) / (-6.0f * 128.0f);
			  tmp[2] = (pt[3] - pt[2]) / 384.0f;
			 _out.coef[1] = (-200*tmp[0] +8*tmp[1] + 2608*tmp[2])/13.0f; //
			 _out.coef[3] = (200*tmp[0] -8*tmp[1] -112*tmp[2]) /13.0f; // -330/13

			 tmp[0] = (pt[5] + pt[0]) / (768.0f);
			 tmp[1] = (pt[4] + pt[1]) / (-2.0f * 128.0f);
			 tmp[2] = (pt[3] + pt[2])/ (384.0f);
			 _out.coef[0] = 9*tmp[0] +25*tmp[1] +225*tmp[2] ;
			 _out.coef[2] = -10*tmp[0] -26*tmp[1] -34*tmp[2];
			 _out.coef[4] = tmp[0] + tmp[1] +tmp[2];
		 		}



			 break;
			case 7:{


		 }
				break;
			case 8:{
				C tmp[4];

		 		}






			 break;
			 case 10:{

		 		}


			 break;
		}
		return(_out);
	}

LFHTEMP	PolyThing<C,0u> PolyThing<C,0u>::interpolate_2free_middlefit(C* pt, int nbpts){

		PolyThing<C,0u> _out;
		_out.setOrder(nbpts-3);
		switch(nbpts){
			case 1:
			case 2: exit(1);
			case 3: _out.coef[0] = pt[1];  break; // LAME!
			case 4:{// (lame, ignores the 2 other pts!)DONE
			 _out.coef[1] = (pt[2] - pt[1])/ 2.0f;
			 _out.coef[0] = (pt[2] + pt[1]) /2.0f;
		 }
			break;
			case 5:{

		 }
		 break;
			case 6:{// DONE
				C tmp[3];
			  tmp[0] = (pt[5] - pt[0]) / (5.0f * 768.0f);
			  tmp[1] = (pt[4] - pt[1]) / (-6.0f * 128.0f);
			  tmp[2] = (pt[3] - pt[2]) / 384.0f;
			 _out.coef[1] = (-200*tmp[0] +8*tmp[1] + 2608*tmp[2])/13.0f; //
			 _out.coef[3] = (200*tmp[0] -8*tmp[1] -112*tmp[2]) /13.0f; // -330/13

			 tmp[0] = (pt[5] + pt[0]) / (768.0f);
			 tmp[1] = (pt[4] + pt[1]) / (-2.0f * 128.0f);
			 tmp[2] = (pt[3] + pt[2])/ (384.0f);
			 _out.coef[0] = 9*tmp[0] +25*tmp[1] +225*tmp[2] - (117.0f/5)*(tmp[0] + tmp[1] +tmp[2]) ; //117/5
			 _out.coef[2] = -10*tmp[0] -26*tmp[1] -34*tmp[2]+ (122.0f/5)*(tmp[0] + tmp[1] +tmp[2]); // -122/5
		 		}



			 break;
			case 7:{

		 }
				break;
			case 8:{

		 		}






			 break;
			 case 10:{

		 		}


			 break;
		}
		return(_out);
	}




	// DOUBLE XOR LIST!

	template<class H, class R, int Hoffoff, int Roffoff>
	dxorptr<H,R,Hoffoff,Roffoff>::dxorptr() {
		hosnext = ((LFH_address)this);
		resprev = hosnext ^ (((LFH_address)this) - R::member_offsets[Roffoff]);
	}

	template<class H, class R, int Hoffoff, int Roffoff>
	dxorptr<H,R,Hoffoff,Roffoff>::dxorptr(const LFH_address dual) : resprev(dual), hosnext(dual){}
	template<class H, class R, int Hoffoff, int Roffoff>
	dxorptr<H,R,Hoffoff,Roffoff>::dxorptr(const LFH_address rp, const LFH_address hn) : resprev(rp), hosnext(hn){}

	template<class H, class R, int Hoffoff, int Roffoff>
	const H* dxorptr<H,R,Hoffoff,Roffoff>::operator->() const{
		   return((const H*)(*this));
		} // - H::member_offsets[Hoffoff]

	template<class H, class R, int Hoffoff, int Roffoff>
	 dxorptr<H,R,Hoffoff,Roffoff>::operator const H*() const{
		dxorptr<H,R,Hoffoff,Roffoff>* prev = (dxorptr<H,R,Hoffoff,Roffoff>*)(resprev ^ (((LFH_address)this) - R::member_offsets[Roffoff]));
		if (prev == this) return(NULL);
		return((const H*)(prev->hosnext ^ ((LFH_address)this)));
		}

	template<class H, class R, int Hoffoff, int Roffoff>
	dxorptr<H,R,Hoffoff,Roffoff>& dxorptr<H,R,Hoffoff,Roffoff>::operator=(const H* target){
		clear();
		if (target == NULL) return(*this);
		dxorptr<H,R,Hoffoff,Roffoff>* link = (dxorptr<H,R,Hoffoff,Roffoff>*)(((LFH_address)target) + H::member_offsets[Hoffoff]);
		dxorptr<H,R,Hoffoff,Roffoff>* next = (dxorptr<H,R,Hoffoff,Roffoff>*)(link->hosnext ^ ((LFH_address)target));

		hosnext = ((LFH_address)next) ^ ((LFH_address)target);
		resprev = ((LFH_address)link) ^ ((LFH_address)(((LFH_address)this) - R::member_offsets[Roffoff]) );

		link->hosnext ^= ((LFH_address)this)^((LFH_address)next);
		next->resprev ^= ((LFH_address)this)^((LFH_address)link);
		return(*this);
		}


	template<class H, class R, int Hoffoff, int Roffoff>
	const H* dxorptr<H,R,Hoffoff,Roffoff>::clear(){

		dxorptr<H,R,Hoffoff,Roffoff>* prev = (dxorptr<H,R,Hoffoff,Roffoff>*)(resprev ^ (((LFH_address)this) - R::member_offsets[Roffoff]));
		if (prev == this) return(NULL);
		H* tar = (H*)(prev->hosnext ^ ((LFH_address)this));
		dxorptr<H,R,Hoffoff,Roffoff>* next = (dxorptr<H,R,Hoffoff,Roffoff>*)(hosnext ^ ((LFH_address)tar));
		next->resprev ^= (((LFH_address)this) ^ ((LFH_address)prev) );
		prev->hosnext ^= (((LFH_address)this) ^ ((LFH_address)next) );

		return(NULL);
	}


	template<class H, class R, int Hoffoff, int Roffoff>
	dxorlist<H,R,Hoffoff,Roffoff>::dxorlist() : dxorptr<H,R,Hoffoff,Roffoff>(((LFH_address)this) ^ (((LFH_address)this) - H::member_offsets[Hoffoff])) {}

	template<class H, class R, int Hoffoff, int Roffoff>
	unsigned int dxorlist<H,R,Hoffoff,Roffoff>::getSize(){
			unsigned int i= 0;
			dxorptr<H,R,Hoffoff,Roffoff>* cur = (dxorptr<H,R,Hoffoff,Roffoff>*)(((LFH_address)((dxorptr<H,R,Hoffoff,Roffoff>*)this)->hosnext) ^ (((LFH_address)this) - H::member_offsets[Hoffoff]));
			while(cur != this) {
				cur = (dxorptr<H,R,Hoffoff,Roffoff>*)(((LFH_address)cur->hosnext) ^ (((LFH_address)this) - H::member_offsets[Hoffoff]));
				i++;
			}
			return(i);
		}

	template<class H, class R, int Hoffoff, int Roffoff>
	dxorlist<H,R,Hoffoff,Roffoff>::operator const Vector<R*> () const{
		Vector<R*> _out;
			dxorptr<H,R,Hoffoff,Roffoff>* prev = (dxorptr<H,R,Hoffoff,Roffoff>*)this;
			dxorptr<H,R,Hoffoff,Roffoff>* cur = (dxorptr<H,R,Hoffoff,Roffoff>*)(((LFH_address)prev->hosnext) ^ (((LFH_address)this) - H::member_offsets[Hoffoff]));
			while(cur != this) {
				_out.push_back((R*)(((LFH_address)cur->resprev) ^ ((LFH_address)prev)));
				prev = cur;
				cur = (dxorptr<H,R,Hoffoff,Roffoff>*)(((LFH_address)cur->hosnext) ^ (((LFH_address)this) - H::member_offsets[Hoffoff]));
			}
		return(_out);
		}

	template<class H, class R, int Hoffoff, int Roffoff>
	dxorptr<H,R,Hoffoff,Roffoff>& dxorlist<H,R,Hoffoff,Roffoff>::operator=(const H* target){
		exit(1234); // not allowed!
	}

	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::dxorptrarray(){
		int i;
		LFH_address hosnext = (LFH_address)aptr;
		for(i=0;i<nbptr;i++) {
			aptr[0] = dxorptr<H,R,Hoffoff,Roffoff>(hosnext, hosnext ^ (((LFH_address)this) - R::member_offsets[Roffoff]));
			hosnext += sizeof(dxorptr<H,R,Hoffoff,Roffoff>);
		}
		}

	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	const H* dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator[](unsigned int which) const{
		return(NULL);
		}

	/*template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	dxorptr<H,R,Hoffoff,Roffoff>& dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator[](unsigned int which) const{
		return(aptr[which]);
	}
	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	dxorptr<H,R,Hoffoff,Roffoff>& dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator[](unsigned int which){
		return(aptr[which]);
	}*/






	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator const dxorptr<H,R,Hoffoff,Roffoff>& () const{
		return((*this)[0]);
	}

	/*

	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	const H* dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator->() const{
		   return( ((*this)[0]) );
		} // - H::member_offsets[Hoffoff]

	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	 dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator const H*() const{
		   return((H*)((*this)[0]));
		} // - H::member_offsets[Hoffoff]

	template<class H, class R, int Hoffoff, int Roffoff, int nbptr>
	dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>& dxorptrarray<H,R,Hoffoff,Roffoff,nbptr>::operator=(const H* target){
		(*this)[0] = (void*) target;


		}*/

		template<class H, class R, int Hoff, int Roff>
		xptr<H,R,Hoff, Roff>::xptr(): target(NULL){}

		template<class H, class R, int Hoff, int Roff>
		H* xptr<H,R,Hoff, Roff>::operator->() const{return(target);}

		template<class H, class R, int Hoff, int Roff>
		void xptr<H,R,Hoff, Roff>::clear(){}
		template<class H, class R, int Hoff, int Roff>
		xptr<H,R,Hoff, Roff>& xptr<H,R,Hoff, Roff>::operator=(H const * const newt){
			clear();
			target = newt;
            return(*this);
			}

		template<class H, class R, int Hoff, int Roff>
		xptrref<H,R,Hoff, Roff>& xptrref<H,R,Hoff,Roff>::operator=(H * value){
			clear();
			if (value){
			trg.target = value;
			xptr<H,R,Hoff,Roff>* head = (xptr<H,R,Hoff,Roff>*)(((LFH_address)value) + H::member_offsets[Hoff]);
			xptr<H,R,Hoff,Roff>* prev = (xptr<H,R,Hoff,Roff>*)((head->resprev) ^ ((LFH_address)value));
			prev->next = &trg;
			trg.next = head;
			trg.resprev = ((LFH_address)prev) ^ ((LFH_address)res);
			head->resprev = ((LFH_address)&trg) ^ ((LFH_address)value);
			}
			return(*this);
			}
		template<class H, class R, int Hoff, int Roff>
		H * xptrref<H,R,Hoff,Roff>::operator->(){return(trg.target);}

		template<class H, class R, int Hoff, int Roff>
		xptrref<H,R,Hoff,Roff>::operator H* (){return(trg.target);}




template <class T>
Ressource<T>::Ressource(): next(NULL), prev(NULL), target(NULL){}

template <class T>
Ressource<T>::Ressource(T tar): next(NULL), prev(NULL), target(tar){}

template <class E,class T>
DataManager<E,T>::DataManager(): BaseManager() {}

template <class E,class T>
DataManager<E,T>::DataManager(int capacity): BaseManager() {}

template <class E,class T>
DataManager<E,T>::~DataManager() {hash_ressource.clear();}

template <class E,class T>
T DataManager<E,T>::getRessource(E key){
	if (hash_ressource.find(key) == hash_ressource.end()) return((hash_ressource[key] = loadRessource(key)).target);
	return(hash_ressource[key].target);
}

template <class E,class T>
Ressource<T> DataManager<E,T>::loadRessource(E key){
	return(T());
}

template <class E,class T>
void DataManager<E,T>::freeRessource(E key){

	if (hash_ressource.find(key) != hash_ressource.end()) hash_ressource.erase(hash_ressource.find(key));

}

template <class E,class T>
E DataManager<E,T>::getFreeKey(){
	int guess = 0xFFFFFFFF;

	while( hash_ressource.find((E)guess) != hash_ressource.end()) guess--;
	return((E)guess);
}


template <class E,class T>
T DataManager<E,T>::quickgetRessource(E key){

	if (hash_ressource.find(key) == hash_ressource.end()) return((T)NULL);
	return(hash_ressource[key].target);
}

template <class E,class T>
bool DataManager<E,T>::hasRessource(E key){
	return(hash_ressource.find(key) != hash_ressource.end());
}


template <class E,class T>
void DataManager<E,T>::flushall(){
	while(hash_ressource.end() != hash_ressource.begin()){
		hash_ressource.erase(hash_ressource.begin());
	}
}

template <class E,class T>
void DataManager<E,T>::insertRessource(E key, T data){
	hash_ressource[key] = Ressource<T>(data);
}


#undef LFHTEMP
#define LFHTEMP template<class RES,unsigned int ressource_appends>

LFHTEMP	RES* LFHPrimitive::ResManager<RES,ressource_appends>::get(typename RES::ressource_enum const ID){ // get ressource, and load it with default function if needed
	unsigned int ite = hash_ressource.find(ID);
	RES* f_out;
	if (ite == 0xFFFFFFFF){
		f_out = new RES(ID); LFH_NICE_ALLOCERROR(f_out ,"")
//		printf("ressource %i:%X created! %i %i\n", (int)ID, (intptr_t)f_out, back_hash.getSize(), hash_ressource.getSize() );
		hash_ressource[ID] = f_out;
		back_hash[f_out] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
//		printf("ressource %i:%X created! %i %i\n", (int)ID, (intptr_t)f_out, back_hash.getSize(), hash_ressource.getSize() );
		}else{
		f_out =  hash_ressource.deref(ite);
	   // printf("ASSUMES!\n");
	   // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
		(back_hash[f_out].second)++;
	}
	return f_out;
}
LFHTEMP	RES* LFHPrimitive::ResManager<RES,ressource_appends>::get(typename RES::ressource_enum const ID, unsigned int init_flag){ // get ressource, and load it with default function if needed
	unsigned int ite = hash_ressource.find(ID);
	RES* f_out;
	if (ite == 0xFFFFFFFF){
		f_out = new RES(ID, init_flag); LFH_NICE_ALLOCERROR(f_out ,"")
		//printf("ressource %i:%X created!\n", (int)ID, (intptr_t)f_out);
		hash_ressource[ID] = f_out;
		back_hash[f_out] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
		}else{
		f_out =  hash_ressource.deref(ite);
	   // printf("ASSUMES!\n");
	   // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
		(back_hash[f_out].second)++;
	}
	return f_out;
}
LFHTEMP	typename RES::ressource_enum LFHPrimitive::ResManager<RES,ressource_appends>::getId(RES* res){ // get ID for a loaded ressource
	int ite = back_hash.find(res);
	return (typename RES::ressource_enum) ((ite == 0xFFFFFFFF) ? 0 : back_hash.deref(ite).first());
}
LFHTEMP	void LFHPrimitive::ResManager<RES,ressource_appends>::free(RES* res ){ // frees a loaded ressource
		int ite = back_hash.find(res);
		if (--(back_hash.deref(ite).second) == 0){
			unsigned int ite2 = hash_ressource.find(back_hash.deref(ite).first);
			delete(res);
			//printf("resource %i:%X is deleted!\n", (int)back_hash.deref(ite).first ,(intptr_t)res);
			hash_ressource.erase_from_iterator(ite2);

          //  printf("PRE-DEL!\n");
          //  for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
			back_hash.erase_from_iterator(ite);
			}
		}

LFHTEMP	typename RES::ressource_enum LFHPrimitive::ResManager<RES,ressource_appends>::put(RES* res){ // stores a custom ressource, needs an new custom ID!
	pair<EnumBox<typename RES::ressource_enum>, unsigned int> tkey;
	tkey.second = 1;
	unsigned int ite;
	do{
	tkey.first = RES::randomKey();
	ite = hash_ressource.find(tkey.first);
	} while(ite != 0xFFFFFFFF);

	hash_ressource[tkey.first] = res;
	back_hash[res] = tkey;
return(tkey.first.data);}

LFHTEMP	void LFHPrimitive::ResManager<RES,ressource_appends>::put_at(RES* res, typename RES::ressource_enum key){ // stores a custom ressource, needs an new custom ID!
	pair<EnumBox<typename RES::ressource_enum>, unsigned int> tkey;
	tkey.second = 1;
	unsigned int ite;
	tkey.first = key;
	if (hash_ressource.find(tkey.first) != 0xFFFFFFFF) {printf("tried to load to an existing slot!\n"); return;}
	hash_ressource[tkey.first] = res;
	back_hash[res] = tkey;
}

#undef LFHTEMP
#define LFHTEMP template<class RES>

LFHTEMP	RES* LFHPrimitive::ResManager<RES,1u>::get(typename RES::ressource_enum const ID){ // get ressource, and load it with default function if needed
    uint32_t ite = hash_ressource.find(ID);
    if (ite == 0xFFFFFFFF){
        FILE* f = RES::getFileHanddle(ID);
        if (1 != fread(&ite,sizeof(uint32_t), 1, f)) {
            if (ID >= RES::getNbInstances()) printf("Tring to load entry outside range (%i >= %i)\n", ID, RES::getNbInstances());
            else printf("Saved entry for %i is corrupted!\n", ID);
            exit(1);
        }
        //printf("read %i for var size %i/%i\n",ite, ID, RES::getNbInstances() ); fflush(stdout);
        AppendixPtr<RES>& t_out = hash_ressource[ID].setExtraSize(ite);
        back_hash[t_out.data] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
        t_out->finishload(f, ite);
        return (RES*) hash_ressource[ID].data;
        }else{
        char* f_out = hash_ressource.deref(ite).data;
       // printf("ASSUMES app!\n");
       // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
        (back_hash[f_out].second)++;
        return (RES*)f_out;
    }
}
LFHTEMP	RES* LFHPrimitive::ResManager<RES,1u>::get(typename RES::ressource_enum const ID, unsigned int init_flag){ // get ressource, and load it with default function if needed
    uint32_t ite = hash_ressource.find(ID);
    if (ite == 0xFFFFFFFF){
        FILE* f = RES::getFileHanddle(ID, init_flag);
        if (1 != fread(&ite,sizeof(uint32_t), 1, f)) {
            if (ID >= RES::getNbInstances()) printf("Tring to load entry outside range (%i >= %i)\n", ID, RES::getNbInstances());
            else printf("Saved entry for %i is corrupted!\n", ID);
            exit(1);
        }
       // printf("read %i for var size %i/%i\n",ite, ID, RES::getNbInstances() ); fflush(stdout);
        AppendixPtr<RES>& t_out = hash_ressource[ID].setExtraSize(ite);
        back_hash[t_out.data] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
        t_out->finishload(f, ite,init_flag);
        return (RES*) hash_ressource[ID].data;
        }else{
        char* f_out = hash_ressource.deref(ite).data;
       // printf("ASSUMES app!\n");
       // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
        (back_hash[f_out].second)++;
        return (RES*)f_out;
    }
}
LFHTEMP	typename RES::ressource_enum LFHPrimitive::ResManager<RES,1u>::getId(RES* res){ // get ID for a loaded ressource
    unsigned int ite = back_hash.find( (char*)res);
    return (ite == 0xFFFFFFFF) ?  (typename RES::ressource_enum) 0 : back_hash.deref(ite).first();
}
LFHTEMP	void LFHPrimitive::ResManager<RES,1u>::free(RES* res ){ // frees a loaded ressource
		unsigned int ite = back_hash.find( (char*)res);
		if (--(back_hash.deref(ite).second) == 0){
			unsigned int ite2 = hash_ressource.find(back_hash.deref(ite).first);
			if (ite2 == 0xFFFFFFFF) {fprintf(stderr, "tried to free an unexisting ressource!\n"); exit(1);}
       //     printf("super deldel!\n"); fflush(stdout);
       //     res->show();

       //     hash_ressource.show();
         //   back_hash.show();
         //   fflush(stdout);

			hash_ressource.erase_from_iterator(ite2);

           // printf("PRE-DEL app!\n");
           // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
			back_hash.erase_from_iterator(ite);

            // delete(res);

			}
		}
LFHTEMP	void LFHPrimitive::ResManager<RES,1u>::put_at(RES* res, typename RES::ressource_enum key){ // stores a custom ressource, needs an new custom ID!
	pair<EnumBox<typename RES::ressource_enum>, unsigned int> tkey;
	tkey.second = 1;
	unsigned int ite;
	tkey.first = key;
	if (hash_ressource.find(tkey.first) != 0xFFFFFFFF) {printf("tried to load to an existing slot!\n"); return;}
	hash_ressource[tkey.first] = res;
	back_hash[res] = tkey;
}

LFHTEMP	 LFHPrimitive::ResManager<RES,2u>::ResManager(){buffer_indexes[0] = 0; buffer_indexes[1] =0;}// ressourceStream = tmpfile();}

LFHTEMP	RES* LFHPrimitive::ResManager<RES,2u>::get(typename RES::ressource_enum const ID){ // get ressource, and load it with default function if needed
	unsigned int ite = hash_ressource.find(ID);
	RES* f_out;
	if (ite == 0xFFFFFFFF){
		f_out = new RES(ID, ressourceArgs); LFH_NICE_ALLOCERROR(f_out ,"")
//		printf("ressource %i:%X created! %i %i\n", (int)ID, (intptr_t)f_out, back_hash.getSize(), hash_ressource.getSize() );
		hash_ressource[ID] = f_out;
		back_hash[f_out] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
//		printf("ressource %i:%X created! %i %i\n", (int)ID, (intptr_t)f_out, back_hash.getSize(), hash_ressource.getSize() );
		}else{
		f_out =  hash_ressource.deref(ite);
	   // printf("ASSUMES!\n");
	   // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
		(back_hash[f_out].second)++;
	}
return f_out;}
LFHTEMP	RES* LFHPrimitive::ResManager<RES,2u>::get(typename RES::ressource_enum const ID, unsigned int init_flag){ // get ressource, and load it with default function if needed
	unsigned int ite = hash_ressource.find(ID);
	RES* f_out;
	if (ite == 0xFFFFFFFF){
		f_out = new RES(ID, init_flag); LFH_NICE_ALLOCERROR(f_out ,"")
		//printf("ressource %i:%X created!\n", (int)ID, (intptr_t)f_out);
		hash_ressource[ID] = f_out;
		back_hash[f_out] =  pair<EnumBox<typename RES::ressource_enum>, unsigned int>(ID,1);
		}else{
		f_out =  hash_ressource.deref(ite);
	   // printf("ASSUMES!\n");
	   // for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
		(back_hash[f_out].second)++;
	}
return f_out;}
LFHTEMP	typename RES::ressource_enum LFHPrimitive::ResManager<RES,2u>::getId(RES* res){ // get ID for a loaded ressource
	int ite = back_hash.find(res);
return (typename RES::ressource_enum) ((ite == 0xFFFFFFFF) ? 0 : back_hash.deref(ite).first());}
LFHTEMP	void LFHPrimitive::ResManager<RES,2u>::free(RES* res ){ // frees a loaded ressource
    int ite = back_hash.find(res);
    if (--(back_hash.deref(ite).second) == 0){
        unsigned int ite2 = hash_ressource.find(back_hash.deref(ite).first);
        delete(res);
        //printf("resource %i:%X is deleted!\n", (int)back_hash.deref(ite).first ,(intptr_t)res);
        hash_ressource.erase_from_iterator(ite2);

      //  printf("PRE-DEL!\n");
      //  for(unsigned int i=0;i< back_hash.heap.getSize();i++) printf("%i: %X\n", i, (int) back_hash.heap[i].first.k );
        back_hash.erase_from_iterator(ite);
    }
}
LFHTEMP	typename RES::ressource_enum LFHPrimitive::ResManager<RES,2u>::put(RES* res){ // stores a custom ressource, needs an new custom ID!
    pair<EnumBox<typename RES::ressource_enum>, unsigned int> tkey;
    tkey.second = 1;
    unsigned int ite;
    do{
    tkey.first = RES::randomKey();
    ite = hash_ressource.find(tkey.first);
    } while(ite != 0xFFFFFFFF);

    hash_ressource[tkey.first] = res;
    back_hash[res] = tkey;
return(tkey.first.data);}

LFHTEMP	void LFHPrimitive::ResManager<RES,2u>::put_at(RES* res, typename RES::ressource_enum key){ // stores a custom ressource, needs an new custom ID!
	pair<EnumBox<typename RES::ressource_enum>, unsigned int> tkey;
	tkey.second = 1;
	unsigned int ite;
	tkey.first = key;
	if (hash_ressource.find(tkey.first) != 0xFFFFFFFF) {printf("tried to load to an existing slot!\n"); return;}
	hash_ressource[tkey.first] = res;
	back_hash[res] = tkey;
}

LFHTEMP	AppendixPtr<RES>::AppendixPtr() : data(NULL){}
LFHTEMP	AppendixPtr<RES>::~AppendixPtr(){if (data != NULL) {ExOp::toMemfree(*(RES*) data); delete[](data);}}
LFHTEMP	AppendixPtr<RES>& AppendixPtr<RES>::setExtraSize(unsigned int extra) {data = new char[sizeof(RES) + extra]; /*printf("Allocated [%X-%X]\n", (int)data, ((int)data) + sizeof(RES) + extra);*/ return *this;}
LFHTEMP	AppendixPtr<RES>& AppendixPtr<RES>::toZero(){return this->toMemfree();}
LFHTEMP	AppendixPtr<RES>& AppendixPtr<RES>::toMemfree(){ if (data != NULL){ExOp::toMemfree(*(RES*) data); /*printf("free %X\n", (int)data);*/ delete[](data); data = NULL;} return *this;}
LFHTEMP	AppendixPtr<RES>& AppendixPtr<RES>::operator=(const AppendixPtr<RES>&){  fprintf(stderr, "any assignment for AppendixPtr<> is illegal, use memmove\n");*(int*)(8 >> 4) = 89;exit(1);return *this;}
LFHTEMP	AppendixPtr<RES>& AppendixPtr<RES>::toMemmove(AppendixPtr& other){data = other.data;other.data = NULL; return *this;}

#undef LFHTEMP
#define LFHTEMP template<class MSG>


LFHTEMP template<class A> void ListenerBase_Hashtable<MSG>::addListener(A* target, typename MSG::MSG_KEY whatkey, typename MSG::MSG_FILTER * whatfilter){
        unsigned int j = this->find(whatkey);
        if (j == 0xFFFFFFFF){
            (*this)[whatkey] = Vector<MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> >();
            j =  this->find(whatkey);
        }
        Vector< MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> > &funfun = this->deref(j);
        funfun.push_back();
        MakeitListen<A,MSG, typename MSG::MSG_FILTER> testest(target,whatfilter);
        memcpy(&(funfun[funfun.getSize()-1]),&testest, sizeof(MakeitListen<ListenerExample, unsigned int, typename MSG::MSG_FILTER>)); // MOUHAHAHAHAH !!!
    }

LFHTEMP template<class A> void ListenerBase_Hashtable<MSG>::removeListener(A* target, typename MSG::MSG_KEY whatkey, typename MSG::MSG_FILTER * whatfilter){
        unsigned int j = this->find(whatkey);
        if (j == 0xFFFFFFFF) LFH_exit("could not remove listernerlink, already gone!\n");
        Vector< MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> > &funfun = this->deref(j);
        unsigned int i;
        for(i=0;i< funfun.getSize();i++) {if (target == (A*)funfun[i].target) break;}
        if (i < funfun.getSize()) funfun.pop_swap(i);
        else LFH_exit("could not remove listernerlink, already gone!\n");
    }

LFHTEMP template<class A> void ListenerBase_Hashtable<MSG>::addUniversal_Listener(A* target, typename MSG::MSG_FILTER * whatfilter){
        universal_listeners.push_back();
        MakeitListen<A,MSG, typename MSG::MSG_FILTER> testest(target,whatfilter);
        memcpy(&(universal_listeners[universal_listeners.getSize()-1]),&testest, sizeof(MakeitListen<ListenerExample, unsigned int, typename MSG::MSG_FILTER>)); // MOUHAHAHAHAH !!!
    }

LFHTEMP template<class A> void ListenerBase_Hashtable<MSG>::removeUniversal_Listener(A* target, typename MSG::MSG_FILTER * whatfilter){
        unsigned int i;
        for(i=0;i<universal_listeners.getSize();i++) if (universal_listeners[i].target == target) break;
        if (i == universal_listeners.getSize()) LFH_exit("could not remove listernerlink, already gone!\n");
        else{
            if (i != universal_listeners.getSize()-1) memcpy(&(universal_listeners[i]), &(universal_listeners[universal_listeners.getSize()-1]),sizeof(MakeitListen<ListenerExample, unsigned int, typename MSG::MSG_FILTER>)); // MOUHAHAHAHAH !!!

        }
    }


LFHTEMP void ListenerBase_Hashtable<MSG>::dispatch(MSG& msg){
        unsigned int j = this->find(msg.msg_key);
        if (j != 0xFFFFFFFF){
            Vector<MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> > &funfun = this->deref(j);
            for(unsigned int i=0;i< funfun.getSize();i++) ((ListenerPointer<MSG>*)(&(funfun[i])))->Listen(msg);
            }
        for(j=0;j<universal_listeners.getSize();j++) ((ListenerPointer<MSG>*)(&(universal_listeners[j])))->Listen(msg);
    }

LFHTEMP template<class A> void ListenerBase_Alias<MSG>::addListener(A* target){
        // find valid aliaz
        unsigned int j = this->find(target->alias);
        while((j != 0xFFFFFFFF)||(ExOp::isZero(target->alias))){
            ExOp::toRand(target->alias);
            j = this->find(target->alias);
        }
        (*this)[target->alias] = MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER>();
        j = this->find(target->alias);
        MakeitListen<ListenerExample, unsigned int, typename MSG::MSG_FILTER> &funfun = this->deref(j);
        MakeitListen<A,MSG, typename MSG::MSG_FILTER> testest(target);
        memcpy(&funfun,&testest, sizeof(MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER>)); // MOUHAHAHAHAH !!!
    }

LFHTEMP void ListenerBase_Alias<MSG>::dispatch(MSG& msg){
       unsigned int j = this->find(msg.key);
        if (j != 0xFFFFFFFF){
            MakeitListen<ListenerExample, MessageExample, typename MSG::MSG_FILTER> &funfun = this->deref(j);
            ((ListenerPointer<MSG>*)(&funfun))->Listen(msg);
        }
    }


LFHTEMP void ListenerBase_Alias<MSG>::show(FILE*f, int level) const{
    printf("ListenerBase_Alias, %i entries!\n", this->size());
    }

template<class A, class MSG> void MakeitListen<A,MSG,void>::Listen(MSG& msg) const{ ((A*)this->target)->listen(msg);}
template<class A, class MSG, class MSG_FILTER> void MakeitListen<A,MSG,MSG_FILTER>::Listen(MSG& msg) const{if (msg.validFilter(msg_filter)) ((A*)this->target)->listen(msg);}

#undef LFHTEMP
#define LFHTEMP template<class MSG, class RES, unsigned int DIM>

LFHTEMP void ListenerBase_Hyperspace<MSG,RES,DIM>::addListener(Listener_Hyperspace<MSG,RES, DIM>* listener, HyperPosition<RES, DIM, 0u> &where){
    // check if it is already inside!
   // typename HierarchicalTree<RES, DIM, 0u,  Listener_Hyperspace<MSG, RES, DIM>* >::Iterator ite = listeners.find_first();

    // then insert!
    listeners.insert( KeyElem< HyperCursor<RES,DIM, 0u>, Listener_Hyperspace<MSG,RES, DIM>* >(where,listener));
}

LFHTEMP void ListenerBase_Hyperspace<MSG,RES,DIM>::dispatch(Tuple<RES,DIM> const &point, MSG& msg){
    typename HierarchicalTree<RES, DIM, 0u,  Listener_Hyperspace<MSG, RES, DIM>* >::Iterator ite = listeners.first();
    for(;ite.isValid(); ++ite) (*ite).d->listen(msg);

    }

#undef LFHTEMP
#define LFHTEMP template<class D, class T>

LFHTEMP void ParamQueue<D,T>::insert_async(const T& t, const D& d){
   /* if (auto accesss = mutx.tryWrite()){
        heap.addEntry(t);
        hmap.addEntry(t,d);
    }else{

    }*/
}

LFHTEMP void ParamQueue<D,T>::remove_async(const T& t, const D& d){
    /*if (auto accesss = mutx.tryWrite()){
        hmap.removeEntry(t,d);
    }else{

    }
    hmap.addEntry(t,d);*/
}


LFHTEMP bool ParamQueue<D,T>::get(const T& t, D& d){
    while((!heap.isEmpty())&&((heap.top() - t) <= 0)){
        uint32_t ite = hmap.find(heap.pop());
        if (ite != 0xFFFFFFFF){
            d = hmap.deref(ite);
            hmap.erase_from_iterator(ite);
            return true;
        }
    }
return false;}


template<class prec, unsigned int nbchan>
void loadWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im){
	FILE* in = fopen(path,"rb+");
	if (in == NULL){
		LFH_STOPPOPUP("Close da file!");
		in = fopen(path,"rb+");
	}
	if (in == NULL) LFH_ENDPOPUP("Grr you did not!");

	char buffer[65536*4];
	fread(buffer,sizeof(char),44,in);

	if (nbchan != *((short*)(buffer + 22))) LFH_ENDPOPUP("Number of channels of opened file is unexpected!");

	unsigned int i = (*((short*)(buffer + 34))) >> 3;
	i *= nbchan;
	i = *((int*)(buffer + 40)) / i; // computed length

//	printf("%i\t%i\n",*((int*)(buffer + 40)) , (*((int*)(buffer + 16))));
//	printf("%i\n",i);
	im.setSizes(&i);

	fread(&(im.data[0]),sizeof(prec) * nbchan, im.dims[0], in);

	fclose(in);



	}
template<class prec, unsigned int nbchan>
void saveWAV(const char* path,DataGrid<Tuple<prec, nbchan>,1> &im){
	FILE* out = fopen(path,"wb+");
	if (out == NULL){
		LFH_STOPPOPUP("Close da file!");
		out = fopen(path,"wb+");
	}
	if (out == NULL) LFH_ENDPOPUP("Grr you did not!");
	char buffer[65536*4];
	memcpy(buffer,"RIFF", sizeof(char)*4);
	memcpy(buffer+8,"WAVEfmt ", sizeof(char)*8);
	memcpy(buffer+36,"data", sizeof(char)*4);
	*((int*)(buffer + 4)) = 36 + sizeof(prec)*nbchan* im.dims[0];
	*((int*)(buffer + 16)) = 16;
	*((short*)(buffer + 20)) = 1;
	*((short*)(buffer + 22)) = nbchan; // stereo
	*((int*)(buffer + 24)) =44100; // samplerate
	*((int*)(buffer + 28)) = 44100 * sizeof(prec) * nbchan; // bitrate
	*((short*)(buffer + 32)) = sizeof(prec) * nbchan; // BlockAlign
	*((short*)(buffer + 34)) = sizeof(prec) * 8; // bitpersample
	*((int*)(buffer + 40)) = sizeof(prec) * nbchan * im.dims[0];
	fwrite(buffer, sizeof(char),44, out);
	fwrite(&(im.data[0]),sizeof(prec) * nbchan, im.dims[0], out);
	fclose(out);
	}

template<class prec, int nbchan>
void lowpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out, DataGrid<Tuple<prec, nbchan>,1> &sour, double freq){
	unsigned int i = sour.dims[0]; // computed length
	out.setSizes(&i);
double factor = 2.0f * freq * M_PI / 44100.0f;
	int hsize = 5 * 44100.0f / freq;
//	if (hsize > 150) hsize = 150;
	int wsize = (hsize<< 1) | 1;
	double* window = new double[wsize]; LFH_NICE_ALLOCERROR(window ,"")
	unsigned int x,y;
	int k;
			double tmp[nbchan];

	double hartor = 2*M_PI / (wsize +1);tmp[0]=0;
	for(k=0;k<wsize;k++){
			window[k] = factor * sinc( (k - hsize)* factor) * (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor))/ (M_PI);
			tmp[0] += window[k];
	//		printf("%i\t%f\t%f\n",k,window[k],tmp[0]);
		}
		tmp[0] = 1.0f / tmp[0];
	for(k=0;k<wsize;k++) window[k] *= tmp[0];

		for(y=0;y<i;y++){
			memset(tmp,'\0',sizeof(double) * nbchan);
			for(x=(y - hsize >=0) ? y - hsize : 0; (x <= hsize + y) && (x < i); x++){
					for(k=0;k<nbchan;k++) tmp[k] += window[x+hsize - y] * sour(&x)[k];
			}
		for(k=0;k<nbchan;k++) out(&y)[k] = (prec)tmp[k];
	}
}

template<class prec, int nbchan>
void highpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out, DataGrid<Tuple<prec, nbchan>,1> &sour, double freq){
	unsigned int i = sour.dims[0]; // computed length
	out.setSizes(&i);
double factor = 2.0f * freq * M_PI / 44100.0f;
	int hsize = 5 * 44100.0f / freq;
//	if (hsize > 150) hsize = 150;
	int wsize = (hsize<< 1) | 1;
	double* window = new double[wsize]; LFH_NICE_ALLOCERROR(window ,"")
	unsigned int x,y;
	int k;
			double tmp[nbchan];

	double hartor = 2*M_PI / (wsize +1);tmp[0]=0;
	for(k=0;k<wsize;k++){
			window[k] = -factor * sinc( (k - hsize)* factor) * (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor))/ (M_PI);
			tmp[0] += window[k];
//			printf("%i\t%f\t%f\n",k,window[k],tmp[0]);
		}
		tmp[0] = -1.0f / tmp[0];
	for(k=0;k<wsize;k++) window[k] *= tmp[0];
		window[hsize] += 1.0f;
		for(y=0;y<i;y++){
			memset(tmp,'\0',sizeof(double) * nbchan);
			for(x=(y - hsize >=0) ? y - hsize : 0; (x <= hsize + y) && (x < i); x++){
					for(k=0;k<nbchan;k++) tmp[k] += window[x+hsize - y] * sour(&x)[k];
			}
		for(k=0;k<nbchan;k++) out(&y)[k] = (prec)tmp[k];
	}
}




template<class prec, int nbchan>
void lowandhighpassfilter(DataGrid<Tuple<prec, nbchan>,1> &out, DataGrid<Tuple<prec, nbchan>,1> &sour, double under, double above){
	unsigned int i = sour.dims[0]; // computed length
	out.setSizes(&i);
double factoru = 2.0f * under * M_PI / 44100.0f;
double factora = 2.0f * above * M_PI / 44100.0f;
	int hsize = (int)(5 * 44100.0f / ((under < above) ? under : above));
//	if (hsize > 150) hsize = 150;
	int wsize = (hsize<< 1) | 1;
	double* window = new double[wsize]; LFH_NICE_ALLOCERROR(window ,"")
	unsigned int x,y;
	int k;
	double tmp[nbchan];

	double hartor = 2*M_PI / (wsize +1);tmp[0]=0;
	double tmpd[2];
	tmpd[0] =0;
	tmpd[1] =0;
	for(k=0;k<wsize;k++){
		tmp[0] = (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
		tmpd[0] += tmp[0] * sinc( (k - hsize)* factoru);
		tmpd[1] += tmp[0] * sinc( (k - hsize)* factora);
		}

		if (under > above){
	tmpd[0] = -1.0f / tmpd[0];
	tmpd[1] = 1.0f / tmpd[1];
	}else{
	tmpd[0] = 1.0f / tmpd[0];
	tmpd[1] = -1.0f / tmpd[1];
		}
	for(k=0;k<wsize;k++){
			window[k] = (tmpd[1] * sinc( (k - hsize)* factora) + tmpd[0] * sinc( (k - hsize)* factoru))*(0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
		}

		if (under > above) window[hsize] += 1.0f;

		for(y=0;y<i;y++){
			memset(tmp,'\0',sizeof(double) * nbchan);
			for(x=(y - hsize >=0) ? y - hsize : 0; (x <= hsize + y) && (x < i); x++){
					for(k=0;k<nbchan;k++) tmp[k] += window[x+hsize - y] * sour(&x)[k];
			}
		for(k=0;k<nbchan;k++) out(&y)[k] = (prec)tmp[k];
	}
}



template<class prec, int nbchan, Tuple_flag flag> static void fourierfilter(DataGrid<Tuple<prec, nbchan, flag>,1> &out,DataGrid<Tuple<prec, nbchan, flag>,1> &sour, double freq){
	unsigned int i = sour.dims[0]; // computed length


//	Tuple<Tuple<prec, nbchan, flag>, i,0>* track;

//	Tuple<mycomplex, super,0> bluewindow = (*track).bluesteinWindow<super>();
//	Tuple<mycomplex, super,0> FT = (*track).fourierTransform(bluewindow);
//	printf("%i\n");
	Tuple<mycomplex,16> haha = Tuple<mycomplex,16>::bluesteinWindow(i);
//	Vector<mycomplex,flag> blue = Vector<mycomplex,flag>(haha);

//	FFout2 = FFout.invfourierTransform(bluewindow);


	}


template<class prec, int nbchan>
void rescale(DataGrid<Tuple<prec, nbchan>,1> &out,DataGrid<Tuple<prec, nbchan>,1> &sour, double freq){
	unsigned int i = sour.dims[0]; // computed length

	double tmp = i * freq;
	out.setSizes(&i);


	}


template<class prec, int nbchan>
void extractFreqfromWAV(DataGrid<Tuple<prec, nbchan>,1> &out, DataGrid<Tuple<prec, nbchan>,1> &sour, double freq){
	unsigned int i = sour.dims[0]; // computed length
	out.setSizes(&i);

	double quart_period = 11025.0f / freq;
	int nbp = 4 + (int)(i / quart_period);
	Tuple<double, nbchan>* buffer = new Tuple<double, nbchan>[nbp];

	memset(buffer,'\0', sizeof(Tuple<double, nbchan>) * nbp);
	double sincos[4];
	unsigned int j;
	int off[2];
	off[0] =2;
	off[1] =1;
	int k,l;
	j=0;


	double factor = 2.0f * freq * M_PI / 44100.0f;
	double envfact = 0.5f / sin(factor * 0.5f);

/*
	sincos[0] = (sin((j + 0.5f) * factor)) * envfact;
	sincos[1] = (cos((j + 0.5f) * factor) - 1.0f) * envfact;

	for(k=0;k<nbchan;k++) {
		sincos[2] = sour(&j)[k];
		buffer[0][k] = sincos[0]*sincos[2];//sour(&j)[k];
		buffer[1][k] = sincos[1]*sincos[2];//sour(&j)[k];
	}
	*/


//	printf("%f mmm!\n", quart_period);
//	exit(1);

	double interpobuf[nbchan][8];
	double tmp[4];
	double polynomial[8];

	double hsize = 24.0f; // nb quarter window size

	memset(interpobuf,'\0',sizeof(double)*nbchan*8);

	for(j = 0;j<3;j++) {
		for(k=0;k<nbchan;k++) interpobuf[k][j] = (double) (sour(&j)[k]);
		}
	for(j = 3;j<i+3;j+=4){
		if (j < i) for(k=0;k<nbchan;k++) interpobuf[k][j & 7] = (double) (sour(&j)[k]);
		else for(k=0;k<nbchan;k++) interpobuf[k][j & 7] = (double) 0.0f;
		j -= 3;

			for(k=0;k<nbchan;k++){
				tmp[0] = (interpobuf[k][j & 7] + interpobuf[k][(j+1) & 7]) / (-92160.0f);
				tmp[1] = (interpobuf[k][(j+7) & 7] + interpobuf[k][(j+2) & 7]) / (18432.0f);
				tmp[2] = (interpobuf[k][(j+6) & 7] + interpobuf[k][(j+3) & 7]) / (-10240.0f);
				tmp[3] = (interpobuf[k][(j+5) & 7] + interpobuf[k][(j+4) & 7]) / (18432.0f) ;
				polynomial[0] = 225*tmp[0]+441*tmp[1]+1225.0f*tmp[2]+11025.0f * tmp[3];
				polynomial[2] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2]+-1891.0f * tmp[3];
				polynomial[4] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ 83.0f* tmp[3];
				polynomial[6] = -tmp[0]-tmp[1]-tmp[2]-tmp[3];
				tmp[0] = (interpobuf[k][j & 7] - interpobuf[k][(j+1) & 7]) / (-7*92160);
				tmp[1] = (interpobuf[k][(j+7) & 7] - interpobuf[k][(j+2) & 7]) / (5*18432);
				tmp[2] = (interpobuf[k][(j+6) & 7] - interpobuf[k][(j+3) & 7]) / (-3*10240.0f);
				tmp[3] = (interpobuf[k][(j+5) & 7] - interpobuf[k][(j+4) & 7]) / (18432.0f);
				polynomial[1] = 225*tmp[0]+441*tmp[1]+1225.0f*tmp[2]+ 11025.0f * tmp[3];
				polynomial[3] = -259*tmp[0]-499*tmp[1]-1299.0f*tmp[2]+ -1891.0f * tmp[3];
				polynomial[5] = 35*tmp[0]+59*tmp[1]+75.0f*tmp[2]+ 83.0f*tmp[3];
				polynomial[7] = -tmp[0]-tmp[1]-tmp[2]-tmp[3];


				/* works!
				printf("1 :\t%f-> %f\n", interpobuf[k][(j+5) & 7], polynomial[0] + polynomial[1]+ polynomial[2]+ polynomial[3]+ polynomial[4]+ polynomial[5]+ polynomial[6]+ polynomial[7]);
				printf("-1:\t%f-> %f\n", interpobuf[k][(j+4) & 7], polynomial[0] - polynomial[1]+ polynomial[2]- polynomial[3]+ polynomial[4]- polynomial[5]+ polynomial[6]- polynomial[7]);
				printf("3 :\t%f-> %f\n", interpobuf[k][(j+6) & 7], polynomial[0] + 3*(polynomial[1]+ 3*(polynomial[2]+ 3*(polynomial[3]+ 3*(polynomial[4]+ 3*(polynomial[5]+ 3*(polynomial[6]+ 3*(polynomial[7]))))))));
				printf("-3:\t%f-> %f\n", interpobuf[k][(j+3) & 7], polynomial[0] - 3*(polynomial[1]- 3*(polynomial[2]- 3*(polynomial[3]- 3*(polynomial[4]- 3*(polynomial[5]- 3*(polynomial[6]- 3*(polynomial[7]))))))));
				*/


				l = (int)((j / quart_period) - hsize)+1;
				if (l < 0) l=0;
				for( ;(l < nbp)&&(l< (int)((j / quart_period) + hsize));l++){
						tmp[0] = 0.125f * M_PI * ((j / quart_period) - l);
//						if (fabs(tmp[0]) > 3.8f) printf("%f\n",tmp[0]);
						buffer[l][k] += polynomial[0] * (0 * sinc( 2* tmp[0]) + 3 * sinc(3 * tmp[0])) *((0.42f/ (M_PI)) - (0.5f/ (M_PI)) * cos( tmp[0] /3) + (0.08f/ (M_PI)) * cos(2*tmp[0] /3));

						// (tmpd[1] * sinc( (k - hsize)* factora) + tmpd[0] * sinc( (k - hsize)* factoru))*(0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
					}

				}




		}



	/*

		sincos[0] = 0;
		sincos[1] = 0;
		sincos[2] = 4.0f * freq * M_PI / (3*44100.0f);
		sincos[3] = 8.0f * freq * M_PI / (3*44100.0f);
		for(k=0;k<wsize;k++){
			sincos[0] += sinc( (k - hsize)* sincos[2]) * (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
			sincos[1] += sinc( (k - hsize)* sincos[3]) * (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
		}
		sincos[0] = 1.0f / sincos[0];
		sincos[1] = -1.0f / sincos[1];
		for(k=0;k<wsize;k++){
			window[k] = (sincos[0] * sinc( (k - hsize)* sincos[2])+ sincos[1] * sinc( (k - hsize)* sincos[3])) * (0.42f - 0.5f * cos( (k+1) * hartor) + 0.08f * cos(2 * (k+1) * hartor));
		}
		window[hsize] += 1;
		LFH_ALIVE; exit(1); // exit in main

*/

	// envfact = 1.0f / (envfact * factor);
//	k=0;	for(j=10;j<nbp;j++) printf("%f\n",buffer[j][k] / (2.0f * envfact) );

		Tuple<Tuple<double, 2>, nbchan>* buffer_phase = new Tuple<Tuple<double, 2>, nbchan>[nbp-5];
		Tuple<Tuple<double, 2>, nbchan>* buffer_ampli = new Tuple<Tuple<double, 2>, nbchan>[nbp-5];
		double vals[6];

		for(j=0;j<nbp-5;j++) for(k=0;k<nbchan;k++){
			vals[0] = buffer[j+2][k] - buffer[j][k];
			vals[1] = buffer[j+3][k] - buffer[j+1][k];
			vals[2] = buffer[j+4][k] - buffer[j+2][k];
			vals[3] = buffer[j+5][k] - buffer[j+3][k];


			vals[4] = vals[1] * vals[1] - vals[0] * vals[2];
			vals[5] = vals[2] * vals[2] - vals[1] * vals[3];

			if (vals[4] * vals[5] <= 0.0f){

				if (vals[4] * vals[5] == 0.0f){
				if ((vals[1] * vals[2]) < 0){
					buffer_phase[j][k][0] = M_PI;
					buffer_phase[j][k][1] = (vals[1] < 0) ? M_PI/2 : -M_PI/2;
				}else{
					buffer_phase[j][k][0] = 0.0f;
					buffer_phase[j][k][1] = 0.0f;
				}
				buffer_ampli[j][k][0] = sqrt(fabs(vals[2] * vals[1]));
				buffer_ampli[j][k][1] = sqrt(fabs(vals[2] / vals[1]));
				}else{
					buffer_phase[j][k][0] = 0.0f;
					buffer_phase[j][k][1] = 0.0f;
					buffer_ampli[j][k][0] = 0.0f;
					buffer_ampli[j][k][1] = 1.0f;
				}

			}else{
				sincos[0] = vals[5] / vals[4];
				buffer_ampli[j][k][1] = sqrt(sincos[0]);
				sincos[1] = vals[1] * vals[2] - vals[3] * vals[0];


				do{

				if (sincos[1] == 0) {buffer_phase[j][k][0] = 0.0f; buffer_phase[j][k][1] = 0.0f; buffer_ampli[j][k][0] = 0.0f; buffer_ampli[j][k][1] = 1.0f;  break;}
				// both roots are shared! => cos(w)^2 = cos(3w)^2
					sincos[2] = vals[1] * buffer_ampli[j][k][1] + vals[2] ;
					sincos[3] = sincos[2] * sincos[2];
					sincos[2] = vals[0] * sincos[0] + vals[3] / buffer_ampli[j][k][1] ;
					sincos[3] -= sincos[2] * sincos[2];
					sincos[3] /= sincos[1] * buffer_ampli[j][k][1]*4;


					sincos[1] = 2 * sincos[3] -1;

					if (fabs(sincos[1]) > 1) {buffer_phase[j][k][0] = 0.0f; buffer_phase[j][k][1] = 0.0f; buffer_ampli[j][k][0] = 0.0f; buffer_ampli[j][k][1] = 1.0f;  break;}


					buffer_phase[j][k][0] = acos(sincos[1])/2;

					sincos[1] = vals[0] * vals[2]* sincos[0] - vals[1] * vals[3];
				if (sincos[1] == 0) {buffer_phase[j][k][0] = 0.0f; buffer_phase[j][k][1] = 0.0f; buffer_ampli[j][k][0] = 0.0f; buffer_ampli[j][k][1] = 1.0f;  break;}
					// both roots are shared! => cos(t+w/2)^2 = cos(t-w/2)^2
					sincos[2] = vals[0] * sincos[0] + vals[2];
					sincos[3] = sincos[2] * sincos[2];
					sincos[2] =  vals[1] * buffer_ampli[j][k][1] + vals[3] / buffer_ampli[j][k][1] ;
					sincos[3] -= sincos[2] * sincos[2];
					sincos[3] /= sincos[1] * 4;


					buffer_ampli[j][k][0] = sqrt(sqrt(vals[4] * vals[5]) / (1 - sincos[3]));

					sincos[1] =  sqrt(sincos[3]);

					sincos[2] = vals[0] * sincos[0] + vals[1] * buffer_ampli[j][k][1];
					sincos[3] = sincos[2] * sincos[2];
					sincos[2] =  vals[2]  + vals[3] / buffer_ampli[j][k][1] ;
					sincos[3] -= sincos[2] * sincos[2];
					sincos[3] /= 4 * (vals[0] * vals[1]* sincos[0]*buffer_ampli[j][k][1] - vals[2] * vals[3] / buffer_ampli[j][k][1]);
					if (sincos[3] < 0.5f) sincos[1] = -sincos[1];

					if (fabs(sincos[1]) > 1) {buffer_phase[j][k][0] = 0.0f; buffer_phase[j][k][1] = 0.0f; buffer_ampli[j][k][0] = 0.0f; buffer_ampli[j][k][1] = 1.0f;  break;}
					buffer_phase[j][k][1] = acos(sincos[1]);


					sincos[0] = M_PI - buffer_phase[j][k][0];
					sincos[2] = buffer_ampli[j][k][0] * cos(buffer_phase[j][k][0] + buffer_phase[j][k][1] /2) ;
					if (vals[2] * sincos[2] < 0.0f) {buffer_phase[j][k][0] += M_PI; sincos[2] += vals[2] / sqrt(buffer_ampli[j][k][1]);}
					else sincos[2] = (vals[2]/ sqrt(buffer_ampli[j][k][1])) - sincos[2];
					sincos[3] = sincos[2] * sincos[2];

					sincos[2] = vals[1] * sqrt(buffer_ampli[j][k][1]) - buffer_ampli[j][k][0] * cos(buffer_phase[j][k][0] - buffer_phase[j][k][1] /2) ;
					sincos[3] += sincos[2] * sincos[2];

					sincos[2] = buffer_ampli[j][k][0] * cos(sincos[0] + buffer_phase[j][k][1] /2)  ;
					if (vals[2] * sincos[2] < 0.0f) {sincos[0] += M_PI; sincos[2] += vals[2]/ sqrt(buffer_ampli[j][k][1]);}
					else sincos[2] = (vals[2]/ sqrt(buffer_ampli[j][k][1])) - sincos[2];
					sincos[1] = sincos[2] * sincos[2];
					sincos[2] = vals[1] * sqrt(buffer_ampli[j][k][1]) - buffer_ampli[j][k][0] * cos(sincos[0] - buffer_phase[j][k][1] /2);
					sincos[1] += sincos[2] * sincos[2];
					if (sincos[1] < sincos[3]) buffer_phase[j][k][0] = sincos[0];




				}while(false);









			}
		//	printf("%f\t%f\t%f\t%f\n", buffer_ampli[j][k][0], buffer_ampli[j][k][1], buffer_phase[j][k][0], buffer_phase[j][k][1]);
		//	printf("no1: %f\t%f\n", vals[0], buffer_ampli[j][k][0] * pow(buffer_ampli[j][k][1], -1.5f) * cos(buffer_phase[j][k][0] - 1.5f * buffer_phase[j][k][1]));
		//	printf("no2: %f\t%f\n", vals[1], buffer_ampli[j][k][0] * pow(buffer_ampli[j][k][1], -0.5f) * cos(buffer_phase[j][k][0] - 0.5f * buffer_phase[j][k][1]));
		//	printf("no3: %f\t%f\n", vals[2], buffer_ampli[j][k][0] * pow(buffer_ampli[j][k][1], 0.5f) * cos(buffer_phase[j][k][0] + 0.5f * buffer_phase[j][k][1]));
		//	printf("no4: %f\t%f\n", vals[3], buffer_ampli[j][k][0] * pow(buffer_ampli[j][k][1], 1.5f) * cos(buffer_phase[j][k][0] + 1.5f * buffer_phase[j][k][1]));



		}



											off[0] =0;
					sincos[1] = 0.99	/ (M_PI * envfact);
			for(j=0;j<i;j++){
					if ((off[0] +2.0f) * quart_period < j) off[0]++;
					sincos[0] = (j / quart_period) - off[0] -1.5f;
					for(k=0;k<nbchan;k++) out(&j)[k] = (prec) (sincos[1] * buffer_ampli[off[0]][k][0] * pow(buffer_ampli[off[0]][k][1], sincos[0]) * sin(buffer_phase[off[0]][k][0] + sincos[0] * buffer_phase[off[0]][k][1]));
				}

	}
		/*

	for(j=1;j<i;j++){
			if (off[0] * quart_period < j){
				if (off[0] & 2){



				sincos[0] = (-1.0f - cos((j - 0.5f) * factor)) * envfact;
				sincos[1] = (cos((j + 0.5f) * factor) - -1.0f) * envfact;
				}else{


				sincos[0] = (1.0f -cos((j - 0.5f) * factor)) * envfact;
				sincos[1] = (cos((j + 0.5f) * factor) - 1.0f) * envfact;
				}
				for(k=0;k<nbchan;k++){
					 sincos[2] = sour(&j)[k];
					 buffer[off[0]][k] += sincos[2] * sincos[0];
					 buffer[off[0]+2][k] += sincos[2] * sincos[1];
				 }
				if (off[0] & 2) for(k=0;k<nbchan;k++) buffer[off[0]][k] *= -1.0f;

				off[0] += 2;

			}else{
				sincos[0] = -sin(j * factor);
				for(k=0;k<nbchan;k++) buffer[off[0]][k] +=  sour(&j)[k] * sincos[0]; // * cos(factor) / factor

			}

			if (off[1] * quart_period < j){

				if (off[1] & 2){

//				if ((off[1] & 65532) == 0)		printf("%i done!  %f \n", j, sin(j * factor) );

				sincos[0] = (-1.0f - sin((j - 0.5f) * factor)) * envfact;
				sincos[1] = (sin((j + 0.5f) * factor) - -1.0f) * envfact;
				}else{
				sincos[0] = (1.0f -sin((j - 0.5f) * factor)) * envfact;
				sincos[1] = (sin((j + 0.5f) * factor) - 1.0f) * envfact;
				}
				for(k=0;k<nbchan;k++){
					 sincos[2] = sour(&j)[k];

					 buffer[off[1]][k] += sincos[2] * sincos[0];
					 buffer[off[1]+2][k] += sincos[2] * sincos[1];
				 }
				if (off[1] & 2) for(k=0;k<nbchan;k++) buffer[off[1]][k] *= -1.0f;
				off[1] += 2;

			}else{
				sincos[0] = cos(j * factor);
				for(k=0;k<nbchan;k++) buffer[off[1]][k] +=  sour(&j)[k]* sincos[0]; // * cos(factor) / factor

			}
		}

	sincos[0] = sin((j - 0.5f) * factor);
	if (off[1] & 2) sincos[0] -=  1.0f;
	else sincos[0] +=  1.0f;
	sincos[1] = cos((j - 0.5f) * factor);
	if (off[0] & 2) sincos[1] -= 1.0f ;
	else sincos[1] +=  1.0f;
	for(k=0;k<nbchan;k++) buffer[1][k] *= 2.0f;
	if (sincos[1] != 0.0f) {
		sincos[1] = fabs(2.0f / sincos[1]);
		for(k=0;k<nbchan;k++) buffer[off[0]][k] *= sincos[1];
		}
	if (sincos[0] != 0.0f) {
		sincos[0] = fabs(2.0f / sincos[0]);
		for(k=0;k<nbchan;k++) buffer[off[1]][k] *= sincos[0];
	}


	*/

			// Special Constructors!!!

// class Spline
#undef LFHTEMP
#define LFHTEMP template <class C>
LFHTEMP	Spline<C>::Spline(): SplineExtra(NULL){

	}
LFHTEMP	void Spline<C>::setExtra(){
	spldata.sort();
	SplineExtra = new C[spldata.getSize()];

	C diffbuf[2];
	int i;
	diffbuf[0] = (spldata[1].d - spldata[0].d) * (3.0f / (spldata[1].k - spldata[0].k));
	for(i=1;i<spldata.getSize();i++){
		diffbuf[(i & 1)] = (spldata[i+1].d - spldata[i].d) * (6.0f / (spldata[i+1].k - spldata[i].k));

		}


	}
LFHTEMP	C Spline<C>::operator()(const double& value){
	if (SplineExtra == NULL) setExtra();
	if ((spldata[last].k < value)&&(value < spldata[last+1])){
		// same interval!

	}



	}
LFHTEMP	void Spline<C>::addPoint(const double& where, const C& data){
	if (SplineExtra != NULL) {delete[](SplineExtra); SplineExtra = NULL;}
	spldata.push_back(KeyElem<double,C>(where,data));
	}




#undef LFHTEMP
#define LFHTEMP template<class O, class I, int _out_dim, int _in_dim>




LFHTEMP void ContinuousFunction<O,I,_out_dim,_in_dim>::derivative(Tuple<O, _out_dim> &, Tuple<I, _in_dim > &, int in_direct)const{

	}



LFHTEMP void ContinuousFunction<O,I,_out_dim,_in_dim>::derivative_from_difference(Tuple<O, _out_dim> &_out, Tuple<I, _in_dim > &_in , int in_direct)const{
			Tuple<I, _in_dim > query1 = _in;

			query1[in_direct] += 0.01f;

			Tuple<O, _out_dim > answer1;
			Tuple<O, _out_dim > answer2;
			(*this)(answer1, query1);
			query1 = _in;
			query1[in_direct] -= 0.01f;
			(*this)(answer2, query1);
			_out = answer1 - answer2;
			_out *= 50.0f;
		}

	LFHTEMP void ContinuousFunction<O,I,_out_dim,_in_dim>::derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &_out, Tuple<I, _in_dim > &_in) const{
		derivativeMatrix_default(_out,_in);
	}

	LFHTEMP void ContinuousFunction<O,I,_out_dim,_in_dim>::derivativeMatrix_default(TMatrix<O, _in_dim, _out_dim> &_out, Tuple<I, _in_dim > &_in) const{
		// calls the derivative function iteratively
		int i,j;
		Tuple<O, _out_dim> deriv;
		for(i=0;i<_in_dim;i++){
			this->derivative(deriv, _in,i);
			for(j=0;j<_out_dim;j++) _out.data[i + j * _in_dim] = deriv[j];
		}
	}

LFHTEMP void ContinuousFunction<O,I,_out_dim,_in_dim>::newtonStep(Tuple<I, _in_dim > &){



	}

#undef LFHTEMP
#define LFHTEMP template<class O, class I, class M, int _out_dim, int _in_dim, int _inner_dim>
/*
LFHTEMP void Chained_Function<O,I,M, _out_dim, _in_dim, _inner_dim>::operator()(Tuple<O, _out_dim> &_out, Tuple<I, _in_dim > &_in)const{
		Tuple<M, _inner_dim> _inner;
		(*first)(_inner, _in);
		(*last)(_out, _inner);
	}

LFHTEMP void Chained_Function<O,I,M, _out_dim, _in_dim, _inner_dim>::derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &_out, Tuple<I, _in_dim > &_in)const{


	TMatrix<M, _in_dim, _inner_dim> _f_out;
	first->derivativeMatrix(_f_out,_in);
	Tuple<M, _inner_dim> _inner;
	(*first)(_inner, _in);

	TMatrix<O, _inner_dim, _out_dim> _s_out;
	last->derivativeMatrix(_s_out,_inner);

	_out = _s_out * _f_out ; // chain-rule joy!
}


#undef LFHTEMP
#define LFHTEMP template<class O, class I, int _out_dim, int _in_dim, int _parameter_dim>

LFHTEMP void Function_Sum<O,I,_out_dim, _in_dim, _parameter_dim>::operator()(Tuple<O, _out_dim> &_out, Tuple<I, _in_dim > &_in) const{

		Tuple<I , _in_dim + _parameter_dim> complete;
		int i;
		int j = params.getSize();
		Tuple<O, _out_dim> tmp;
		if (j == 0){
			for(i=0;i< _out_dim;i++) ExCo<O>::toZero(_out[i]);
		}else{
			j--;
			for(i=0;i< _in_dim;i++) complete[i] = _in[i];

			for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
			(*function)(_out, complete);
			for(j--;j>=0;j--){
				for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
				(*function)(tmp, complete);
				_out += tmp;
			}
		}
	}

LFHTEMP void Function_Sum<O,I,_out_dim, _in_dim, _parameter_dim>::derivative(Tuple<O, _out_dim> &_out, Tuple<I, _in_dim > &_in, int in_direct) const{

	Tuple<I , _in_dim + _parameter_dim> complete;
	int i;
	int j = params.getSize();
	Tuple<O, _out_dim> tmp;
	if (j == 0){
		for(i=0;i< _out_dim;i++) ExCo<O>::toZero(_out[i]);
	}else{
		j--;
		for(i=0;i< _in_dim;i++) complete[i] = _in[i];

		for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
		function->derivative(_out, complete,in_direct);
		for(j--;j>=0;j--){
			for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
			function->derivative(tmp, complete, in_direct);
			_out += tmp;
		}
	}

	}

	template<class A,class B>
	void Convert<A,B>::operator()(A & a , B & b) const{a = (A)b;}

	template<class A,class B>
	void Convert<A,B>::operator()(A & a , const B & b) const{a = (A)b;}

	template<class A>
	void getNorm<A>::operator()(double & _out , A & a) const{_out = ExCo<A>::norm(a);}


	template<class A>
	void getPositiveNorm<A>::operator()(double & _out , A & a) const{_out = (a > ExOp::mkzero<A>()) ? ExOp::norm(a) : 0.0f;}



LFHTEMP void Function_Sum<O,I,_out_dim, _in_dim, _parameter_dim>::derivativeMatrix(TMatrix<O, _in_dim, _out_dim> &_out, Tuple<I, _in_dim > &_in) const{

	Tuple<I , _in_dim + _parameter_dim> complete;
	int i;
	int j = params.getSize();
	TMatrix<O, _in_dim, _out_dim> tmp;
	if (j == 0){
		for(i=0;i< _out_dim * _in_dim ;i++) ExCo<O>::toZero(_out.data[i]);
	}else{
		j--;
		for(i=0;i< _in_dim;i++) complete[i] = _in[i];

		for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
		function->derivativeMatrix(_out, complete);
		for(j--;j>=0;j--){
			for(i=0;i<  _parameter_dim;i++) complete[i + _in_dim] = params[j][i];
			function->derivativeMatrix(tmp, complete);
			_out += tmp;
		}
	}


	}



#undef LFHTEMP
#define LFHTEMP template<int dims>

	LFHTEMP	void SphereDistanceLikelyhood<dims>::operator()(Tuple<double, 1> &_out, Tuple<double, dims+2 > &_in) const{
		double sum, tmp;
		int i;
		sum = _in[2]*_in[2];
		for(i=3;i<dims+2;i++) sum += _in[i]*_in[i];
		sum = (_in[0] - sqrt(sum))/ _in[1] ;
		_out[0] = - log(_in[1]) -0.5f * (log(2 * M_PI) + sum * sum);
	}


	LFHTEMP	void SphereDistanceLikelyhood<dims>::derivative(Tuple<double, 1> &_out, Tuple<double, dims+2> &_in, int in_direct) const{
		double sum, tmp;
		int i;
		sum = _in[2]*_in[2];
		for(i=3;i<dims+2;i++) sum += _in[i]*_in[i];


		switch(in_direct){
			case 0:

				break;
			case 1:
				_out[0] = -(1.0f / _in[1]) + (  (_in[0]* (_in[0] - 2.0f * sqrt(sum)) + sum ) / (_in[1]*_in[1]*_in[1]));

				break;
			default:
				_out[0] =  _in[in_direct] * (((_in[0] / sqrt(sum)) - 1.0f) / (_in[1] *_in[1])) ;

		}


	}

	LFHTEMP	void SphereDistanceLikelyhood_Derivative<dims>::operator()(Tuple<double, dims+1> &_out, Tuple<double, dims+2 > &_in) const{
		double sum, tmp;
		int i;
		sum = _in[2]*_in[2];
		for(i=3;i<dims+2;i++) sum += _in[i]*_in[i];
		_out[0] = (((_in[0] / sqrt(sum)) - 1.0f) / (_in[1] *_in[1]));
		for(i = 0;i<dims;i++) _out[i+1] = _in[i+2] * _out[0];
		_out[0] = (-1.0f + (_in[0]* (_in[0] - 2.0f * sqrt(sum)) + sum) / (_in[1]*_in[1])) / _in[1];
	}


	LFHTEMP	void SphereDistanceLikelyhood_Derivative<dims>::derivative(Tuple<double, dims+1> &_out, Tuple<double, dims+2> &_in, int in_direct) const{
		double sum, tmp;
		int i;
		sum = _in[2]*_in[2];
		for(i=3;i<dims+2;i++) sum += _in[i]*_in[i];




	}*/


template<class I, ENABLEIF_ARRAY(I,double)> GaussianProcessMK2& GaussianProcessMK2::setKernel(I timestamp, I noisescale, uint32_t length){
	if (length == 0) LFH_exit("give a length plz");
	cached_kernel.setSize(length);
	int i,j,k;
	double tmp;
	for(k=0,i=0;i<cached_kernel.getSize();i++){
		for(j=0;j<i;j++){
			tmp = timestamp[j] - timestamp[i];
			cached_kernel.data[k++] = this->getSignal() * exp( -tmp * tmp / getScale());
		}
		cached_kernel.data[k++] = this->getSignal() + this->getNoise() * noisescale[i];
	}
return *this;}
template<class I, class J, ENABLEIF_ARRAY(I, double), ENABLEIF_ARRAY(J, double)> GaussianProcessMK2& GaussianProcessMK2::setKernel2(I timestamp, J noisescale, uint32_t length){
	if (length == 0) LFH_exit("give a length plz");
	cached_kernel.setSize(length);
	int i,j,k,l;
	Tuple<double> tmp; tmp.setSize(scale.getSize());
	for(k=0,i=0;i<length;i++){
		for(j=0;j<i;j++){
			for(l=0;l<scale.getSize();l++) tmp[l] = timestamp[j * scale.getSize() + l] - timestamp[i * scale.getSize() + l];
			cached_kernel.data[k++] = this->getSignal() * exp( -scale.Xformed_inner_product(tmp));
		}
		cached_kernel.data[k++] = this->getSignal() + this->getNoise() * noisescale[i];
	}
return *this;}
template<class I, ENABLEIF_ARRAY(I,double)> double GaussianProcessMK2::wrLLDerivative(Tuple<double, 4u>& fout, I observ, I timestamp, I noisescale, bool verbose)const{
	fout.toZero();
	int i,j,k;
	Tuple<double> dif; dif.setSize(cached_kernel.getSize());
	for(i=0;i<dif.getSize();i++) dif[i] = observ[i] - this->getMean();

	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	if (verbose) {
		printf("LL %e\n",LL);
		dif.show();
	}
	LL += cached_kernel.log_determinant();
	LL *= -0.5;
	Tuple<double> old = dif;
	dif = cached_kernel.leftDivision(dif);
	double altLL = 0.0;

	// derivative in Signal
	Trianglix<double> der, der2, der3; der.setSize(cached_kernel.getSize()); der2.setSize(cached_kernel.getSize()); der3.setSize(cached_kernel.getSize());
	for(k=0,i=0;i<cached_kernel.getSize();i++){
		for(j=0;j<i;j++,k++){
			double tmp = timestamp[j] - timestamp[i];
			der.data[k] = cached_kernel.data[k] / getSignal();
			der2.data[k] = cached_kernel.data[k] * tmp * tmp / (getScale() * getScale());
			der3.data[k] = 0.0;
		}
		altLL += dif[i] * old[i];
		fout[0] += dif[i];
		fout[1] += noisescale[i] * dif[i] * dif[i]; // + 1.0 / cached_kernel.data[k]);
		der3.data[k] = noisescale[i];
		der.data[k] = 1.0;
		der2.data[k++] = 0.0;
	}
	if (verbose) {
		printf("%e %i\n", altLL,k);
		dif.show();
		der.show();
		der2.show();
		der3.show();
	}
	fout[3] = (der2.Xformed_inner_product(dif) - der2.trace_of_division(cached_kernel)) * 0.5;
	fout[2] = (der.Xformed_inner_product(dif) - der.trace_of_division(cached_kernel)) * 0.5;
	fout[1] -= der3.trace_of_division(cached_kernel);
	fout[1] *= 0.5;
return LL;}



template<class I, class J, class K, ENABLEIF_ARRAY(I, double), ENABLEIF_ARRAY(J, double), ENABLEIF_ARRAY(K, double)>
double GaussianProcessMK2::wrLLDerivative2(Tuple<double>& fout, I observ, J timestamp, K noisescale, bool verbose)const{
	fout.setSize(3 + scale.totsize()).toZero();
	int i,j,k,l,m,n;
	Tuple<double> dif; dif.setSize(cached_kernel.getSize());
	for(i=0;i<dif.getSize();i++) dif[i] = observ[i] - this->getMean();

	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	if (verbose) {
		cached_kernel.show();
		dif.show();
		printf("LL %e\n",LL);
		dif.show();
	}
	LL += cached_kernel.log_determinant();
	LL *= -0.5;
	Tuple<double> old = dif;
	dif = cached_kernel.leftDivision(dif);
	double altLL = 0.0;

	Trianglix< Tuple<double> > ders; ders.setSize(cached_kernel.getSize());

	Tuple<double> tmp; tmp.setSize(scale.getSize());
	// derivative in Signal
	//Trianglix<double> der, der2; der.setSize(cached_kernel.getSize()); der2.setSize(cached_kernel.getSize());
	for(k=0,i=0;i<cached_kernel.getSize();i++){
		for(j=0;j<i;j++,k++){
			ders.data[k].setSize(fout.getSize() - 1);
			for(l=0;l<scale.getSize();l++) tmp[l] = timestamp[j * scale.getSize()  + l] - timestamp[i * scale.getSize() + l];
			ders.data[k][0] = cached_kernel.data[k] / getSignal();
			ders.data[k][1] = 0.0;
			//der.data[k] = cached_kernel.data[k] / getSignal();
			//der2.data[k] = 0.0;
			for(l=2,m=0;m<scale.getSize();m++){
				for(n=0;n<m;n++) ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[n] * 2.0;
				ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[m];
			}
		}
		altLL += dif[i] * old[i];
		fout[0] += dif[i];
		//fout[1] -= noisescale[i] * (dif[i] * dif[i] + 1.0 / cached_kernel.data[k]);
		ders.data[k].setSize(fout.getSize() - 1).toZero();
		ders.data[k][0] = 1.0;
		ders.data[k][1] = noisescale[i];
		//der2.data[k] = noisescale[i];
		//der.data[k] = 1.0;
		k++;
	}


	Tuple<double> trofdif = ders.Xformed_inner_product(dif);
	if (verbose) {
		printf("%e %i\n", altLL,k);
		dif.show();
		for(k=0;k<10;k++){
		printf("v %i:",k);
		ders.data[k].show();}
		printf("those 2\n");
		//der.show();
		//der2.show();
		printf("new xform %e %e\n", trofdif[0],trofdif[1]);
		//printf("old xform %e %e\n", der.Xformed_inner_product(dif),der2.Xformed_inner_product(dif));
		trofdif = ders.trace_of_division(cached_kernel);
		printf("new ofdiv %e %e\n", trofdif[0],trofdif[1]);
		//printf("old ofdiv %e %e\n", der.trace_of_division(cached_kernel),der2.trace_of_division(cached_kernel));
		trofdif = ders.Xformed_inner_product(dif);
	}
	trofdif -= ders.trace_of_division(cached_kernel);

	for(l=0;l< trofdif.getSize();l++) fout[l+1] = trofdif[l] * 0.5;
return LL;}

template<class I, class J, ENABLEIF_ARRAY(I, (Tuple<double, 2u>)), ENABLEIF_ARRAY(J, double)>
double GaussianProcessMK2::wrLLDerivative3(Tuple<double>& fout, I observ, J timestamp, bool verbose)const{
	fout.setSize(3 + scale.totsize()).toZero();
	int i,j,k,l,m,n;
	Tuple<double> dif; dif.setSize(cached_kernel.getSize());
	for(i=0;i<dif.getSize();i++) dif[i] = observ[i][0] - this->getMean();

	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	if (verbose) {
		cached_kernel.show();
		dif.show();
		printf("LL %e\n",LL);
		dif.show();
	}
	LL += cached_kernel.log_determinant();
	LL *= -0.5;
	Tuple<double> old = dif;
	dif = cached_kernel.leftDivision(dif);
	double altLL = 0.0;

	Trianglix< Tuple<double> > ders; ders.setSize(cached_kernel.getSize());

	Tuple<double> tmp; tmp.setSize(scale.getSize());
	// derivative in Signal
	//Trianglix<double> der, der2; der.setSize(cached_kernel.getSize()); der2.setSize(cached_kernel.getSize());
	for(k=0,i=0;i<cached_kernel.getSize();i++){
		for(j=0;j<i;j++,k++){
			ders.data[k].setSize(fout.getSize() - 1);
			for(l=0;l<scale.getSize();l++) tmp[l] = timestamp[j * scale.getSize()  + l] - timestamp[i * scale.getSize() + l];
			ders.data[k][0] = cached_kernel.data[k] / getSignal();
			ders.data[k][1] = 0.0;
			//der.data[k] = cached_kernel.data[k] / getSignal();
			//der2.data[k] = 0.0;
			for(l=2,m=0;m<scale.getSize();m++){
				for(n=0;n<m;n++) ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[n] * 2.0;
				ders.data[k][l++] = cached_kernel.data[k] * tmp[m] * tmp[m];
			}
		}
		altLL += dif[i] * old[i];
		fout[0] += dif[i];
		//fout[1] -= noisescale[i] * (dif[i] * dif[i] + 1.0 / cached_kernel.data[k]);
		ders.data[k].setSize(fout.getSize() - 1).toZero();
		ders.data[k][0] = 1.0;
		ders.data[k][1] = observ[i][1]; // noise
		//der2.data[k] = noisescale[i];
		//der.data[k] = 1.0;
		k++;
	}


	Tuple<double> trofdif = ders.Xformed_inner_product(dif);
	if (verbose) {
		printf("%e %i\n", altLL,k);
		dif.show();
		for(k=0;k<10;k++){
		printf("v %i:",k);
		ders.data[k].show();}
		printf("those 2\n");
		//der.show();
		//der2.show();
		printf("new xform %e %e\n", trofdif[0],trofdif[1]);
		//printf("old xform %e %e\n", der.Xformed_inner_product(dif),der2.Xformed_inner_product(dif));
		trofdif = ders.trace_of_division(cached_kernel);
		printf("new ofdiv %e %e\n", trofdif[0],trofdif[1]);
		//printf("old ofdiv %e %e\n", der.trace_of_division(cached_kernel),der2.trace_of_division(cached_kernel));
		trofdif = ders.Xformed_inner_product(dif);
	}
	trofdif -= ders.trace_of_division(cached_kernel);

	for(l=0;l< trofdif.getSize();l++) fout[l+1] = trofdif[l] * 0.5;
return LL;}
template<class I, class J, class K, ENABLEIF_ARRAY(I, double), ENABLEIF_ARRAY(J, double), ENABLEIF_ARRAY(K, double)>
double GaussianProcessMK2::wrLLTimeDerivative(Tuple<double>& fout, I observ, J timestamp, K noisescale, bool verbose)const{
	int i,j,k,l,m,n;
	Tuple<double> dif; dif.setSize(cached_kernel.getSize());
	for(i=0;i<dif.getSize();i++) dif[i] = observ[i] - this->getMean();

	double LL = cached_kernel.Xformed_inner_product_of_inverse(dif);
	if (verbose) {
		cached_kernel.show();
		dif.show();
	}
	LL += cached_kernel.log_determinant();
	LL *= -0.5;
	Tuple<double> old = dif;
	dif = cached_kernel.leftDivision(dif);
	double altLL = 0.0;

	Tuple< Tuple<double> > ders; ders.setSize(cached_kernel.getSize() * scale.getSize());

	Tuple<double> tmp; tmp.setSize(scale.getSize());
	for(l=0;l<ders.getSize();l++) ders[l].setSize(cached_kernel.getSize()).toZero();
	for(k=0,i=0;i<cached_kernel.getSize();i++){
		for(j=0;j<i;j++,k++){ // i and j are sampleIDs
			for(l=0;l<scale.getSize();l++) tmp[l] = timestamp[j * scale.getSize()  + l] - timestamp[i * scale.getSize() + l];
			tmp = scale * tmp; // projection
			for(m=0;m<scale.getSize();m++){ // n and m are dimensions
				ders[i * scale.getSize() + m][j] -= cached_kernel.data[k] * tmp[m];
				ders[j * scale.getSize() + m][i] += cached_kernel.data[k] * tmp[m];
			}
		}
		k++;
	}

	if (verbose) {
		for(l=0;l<ders.getSize();l++) {
			printf("%i,%i: ",l / cached_kernel.getSize(), l % cached_kernel.getSize()); ders[l].show();
		}
	}
	fout.setSize(ders.getSize());
	for(l=0;l<ders.getSize();l++) {
		fout[l] =  dif[0] * ders[l][0];
		for(i=1;i<cached_kernel.getSize();i++) fout[l] += dif[i] * ders[l][i];
	}
	if (verbose) {
		printf("and\n");
		fout.show();
		printf("and\n");
		dif.show();
		printf("and\n");
	}
	for(l=0;l<ders.getSize();l++) {
		fout[l] *= dif[(l / scale.getSize())];
		ders[l] = cached_kernel.leftDivision(ders[l]);
		fout[l] -= ders[l][(l / scale.getSize())];
		fout[l] *= 2.0;
	}
return LL;}

template<class P, class S, class O> void GaussianProcessMK2::LearnParaTask::operator()(P& gp_producer, S& state_iterator, O& observation_dico){
	CurvatureSearchScope dacurv;
	int j;
	Tuple<double> curderiv;
    LL = 0.0;
	if (auto jobite = gp_producer()) do{
		dacurv.init();
		for(j = 0 ;j < 1000;j++){
			jobite->wrLLDerivative3(curderiv, observation_dico[jobite()], state_iterator);
		}
		if (j == 1000) dacurv.wrFinalGuess( jobite->mkParamIterator());
		LL +=  dacurv.getLastValue();
	}while(jobite++);
}

template<class P, class S, class G> void GaussianProcessMK2::CmpStateDerivative::operator()(P& observation_producer, S& state_iterator, G& gp_dico){
	if (auto jobite = observation_producer()) do{

	}while(jobite++);
}



	/*
	// Gaussian Process!

	template<class C, LFHVectorflag d_flag>

	void GaussianProcess::covMatrix(DataGrid<double, 2,> &f_out, Vector<C,d_flag>& data, double (*metric)(const C&, const C&), double noiselevel){
		unsigned int coor[2];
		coor[1] =coor[0] = data.getSize();
		f_out.setSizes(coor);
		double tmp;
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();

		if (ite.first()) do{

			if (ite()[0] == ite()[1]) {
				f_out(ite()) = 1.0f + noiselevel;
			}else{
				tmp = metric(data[ite()[0]], data[ite()[1]]);
				f_out(ite()) = exp(-0.5f * tmp * tmp);
			}

		} while(ite.next());
	}

	template<class C, LFHVectorflag d_flag, LFHVectorflag q_flag>
	void GaussianProcess::covMatrix(DataGrid<double, 2,> &f_out, Vector<C,d_flag>& data, Vector<C,q_flag>& query, double (*metric)(const C&, const C&)){
		unsigned int coor[2];
		coor[0] = data.getSize();
		coor[1] = query.getSize();
		f_out.setSizes(coor);
		double tmp;
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		if (ite.first()) do{
			tmp = metric(data[ite()[0]], query[ite()[1]]);
			f_out(ite()) = exp(-0.5f * tmp * tmp);
		} while(ite.next());
	}

	template<class C, LFHVectorflag d_flag>
	void GaussianProcess::covMatrix(DataGrid<double, 2> &f_out, Vector<WeightElem<C,1>,d_flag>& data, double (*metric)(const C&, const C&), double noiselevel){
		unsigned int coor[2];
		coor[1] =coor[0] = data.getSize();
		f_out.setSizes(coor);
		double tmp;
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();

		if (ite.first()) do{

			if (ite()[0] == ite()[1]) {
				f_out(ite()) = 1.0f + noiselevel;
			}else{
				tmp = metric(data[ite()[0]].e[0], data[ite()[1]].e[0]);
				f_out(ite()) = exp(-0.5f * tmp * tmp) * data[ite()[0]].w[0] * data[ite()[1]].w[0];
			}

		} while(ite.next());
	}

	template<class C, LFHVectorflag d_flag, LFHVectorflag q_flag>
	void GaussianProcess::covMatrix(DataGrid<double, 2> &f_out, Vector<WeightElem<C,1>,d_flag>& data, Vector<C,q_flag>& query, double (*metric)(const C&, const C&)){
		unsigned int coor[2];
		coor[0] = data.getSize();
		coor[1] = query.getSize();
		f_out.setSizes(coor);
		double tmp;
		typedef class DataGrid<double,2>::KeyIterator itetype;
		itetype ite = f_out.getKeyIterator();
		if (ite.first()) do{
			tmp = metric(data[ite()[0]].e[0], query[ite()[1]]);
			f_out(ite()) = exp(-0.5f * tmp * tmp) * data[ite()[0]].w[0];
		} while(ite.next());
	}

	template <class C, Tuple_flag flag>  double twoddistcorr(const Tuple<C,3,flag> *pts, const Tuple<unsigned int ,2,flag> &indexes){
		if (indexes[0] == indexes[1]) return(2.0f);
		Tuple<C,3> diff = pts[indexes[0]] - pts[indexes[1]];
		return(pts[indexes[0]][2]*pts[indexes[1]][2]* exp(-1.25f * ((double) diff[0]* diff[0] + diff[1]*diff[1])));
		}

	template <class C, Tuple_flag flag>  double twoddistcorr_noweight(const Tuple<C,3,flag> *pts, const Tuple<unsigned int ,2,flag> &indexes){
		if (indexes[0] == indexes[1]) return(2.0f);
		Tuple<C,3> diff = pts[indexes[0]] - pts[indexes[1]];
		return(exp(-1.25f * ((double) diff[0]* diff[0] + diff[1]*diff[1])));
		}

	 template<class K, class C, int NBDIM, LFHVectorflag d_flag, LFHVectorflag d_flag_2>
	 void GaussianProcess::GPweightsregression(DataGrid<C, NBDIM> &f_out, const DataGrid<C, NBDIM, d_flag> &p_data, const Vector<K, d_flag_2> &p_keys , double (*p_dist)(const K &a, const K &b), double (*p_weight)(const K &a)){
*/





	// n2 linear regression
	template<class K, class C>
	Vector<C> GaussianProcess::GPweightsregression(const Vector< KeyElem<K,C> > &p_data, const Vector<K> &p_keys){
		Vector<C> f_out;

		DataGrid< double, 2 > front_stuff;

		DataGrid< C, 2 > back_stuff;

		// fill front

		unsigned int dims[2];
		unsigned int coor[2];
		dims[0] = p_data.getSize();
		dims[1] = p_keys.getSize();
		front_stuff.setSizes(dims);
		double norm;
		for(coor[1]=0;coor[1]<dims[1];coor[1]++){
			for(coor[0]=0;coor[0]<dims[0];coor[0]++){
				norm = p_data[coor[0]].k.getNorm(p_keys[coor[1]]);
				front_stuff(coor) = norm;
			}
		}


		dims[0] = 1;
		dims[1] = p_data.getSize();
		back_stuff.setSizes(dims);
		coor[0]=0;
		for(coor[1]=0;coor[1]<dims[1];coor[1]++){
			back_stuff(coor) = p_data[coor[1]].d;
		}

		ExOp::show(front_stuff);
		ExOp::show(back_stuff);

		// compute K * COV-1 * Y


		return(f_out);

		/*	LinkAssert< (NBDIM>1) > ass;
		unsigned int dims[2];

		DataGrid<C, NBDIM> residual;
		FuncGrid<double, K,   , 2> corr;



		dims[0] = p_data.dims[0];
		dims[1] = p_data.dims[1];
		residual.setSizes(dims);
		// sqrt grouping n^2
		// local regressions, sqrt * sqrt^3
		// iterative converge n^2

		// STEP 1: Cluster instances!

		// O(N^2)


		Vector<int> partition;
		Vector<int> p_maxes;
		int s = p_keys.getSize();
		partition.setSize(s);
		int i,j,k,l;
		for(i=0;i<dims[1];i++) corr.data.push_back(p_keys[i]);
		{//:
		Forest< unsigned int ,3> hierach;
		unsigned int cur;
		hierach.partition_squarreroot_nbcluster(p_keys, p_dist, (int) pow(s, 0.6f) );


		stack<unsigned int> nodes;
		p_maxes.setSize(hierach.nbroots);j=0;
		for(i=0;i<hierach.nbroots;i++){
			nodes.push(i);
	//		printf("nroot\n");
			while(!nodes.empty()){
				cur = nodes.top();nodes.pop();
	//			printf("visitnroot %i\n",cur);
				if (hierach.hasLeft(cur)){
					nodes.push(hierach.getLeft(cur));
					nodes.push(hierach.getRight(cur));
					}else{
					partition[j++] = hierach[cur].second;
	//				printf("hello %i\n", cur);
				}
				}
				p_maxes[i] =j;
			}
		}//:




		// STEP 2: Com'e inverses of partial linear systems!
		// O(N^2)

		j=0;



		int maxsize = 1+ (int) sqrt((double)s);
		int nbsubl = p_maxes.getSize();
		DataGrid<double, 2>* HH = new DataGrid<double, 2>[nbsubl+1];

		DataGrid<double, 2> tmpsys;

		unsigned int coor[2];
		unsigned int cooralt[2];

		// fill modified housederler matrices
		int lastsize =0;

		corr.fun = &twoddistcorr;
		k=0;
		double tmp, sum, average;

		for(i=0;i<nbsubl;i++){
			k = (i ==0) ? 0 : p_maxes[i -1];
			dims[1] = dims[0] = p_maxes[i] -k;
			tmpsys.setSizes(dims);
			for(coor[1]=0;coor[1]<dims[0];coor[1]++) {
				cooralt[1] = partition[coor[1]+k];
			for(coor[0]=0;coor[0]<dims[0];coor[0]++) {
				cooralt[0] = partition[coor[0]+k];
				tmpsys(coor) = corr(cooralt);
			}}

			HH[i].makeInverse(tmpsys);

	//		printf("inverses!\n");
	//		HH[i].show();
			}
		dims[1] = dims[0] = nbsubl;

		tmpsys.setSizes(dims);
		ExOp::toZero(tmpsys);

		cooralt[0] =0;
		cooralt[1] =0;
		for(coor[1]=0;coor[1]<corr.data.getSize();coor[1]++){
			while(coor[1]>=p_maxes[cooralt[1]]) cooralt[1]++;
			cooralt[0] =cooralt[1];
			coor[0] = coor[1];
			tmpsys(cooralt) += corr(coor);
			for(coor[0]++;coor[0]<corr.data.getSize();coor[0]++){
				while(coor[0]>=p_maxes[cooralt[0]]) cooralt[0]++;
				if (cooralt[0] == cooralt[1]) tmpsys(cooralt) += corr(coor) * 2.0f;
				else tmpsys(cooralt) +=  corr(coor);
				}
			}

		for(coor[1]=0; coor[1]<nbsubl;coor[1]++){
			cooralt[0] = coor[1];
			coor[0]=coor[1]+1;
			cooralt[1] = coor[0];
			for(;coor[0]<nbsubl;coor[0]++,cooralt[1]++){
				tmpsys(cooralt) = tmpsys(coor);
			}
		}

	//	printf("inverses!\n");
	//		tmpsys.show();
		HH[nbsubl].makeInverse(tmpsys);
	//	HH[nbsubl].show();
		DataGrid<double, 2>* RD = new DataGrid<double, 2>[nbsubl+1];
		dims[0] =p_data.dims[0];
		for(i=0;i<=nbsubl;i++){
				dims[1] = HH[i].dims[0];
				RD[i].setSizes(dims);
			}

		int loop;

		for(i=0;i<p_maxes.getSize();i++) printf("%i%c", p_maxes[i], i == p_maxes.getSize()-1 ? '\n':'\t');
		for(i=0;i<=nbsubl;i++) printf("%i%c", HH[i].dims[0], i == nbsubl ? '\n':'\t');
		// STEP 3: Iterative refinement using partial inverses!

		double* cachedweights = new double[maxsize];
		for(loop =0; loop< 25;loop++){
			// get residual
			residual.LinearSystem_residual(f_out,p_data, corr); // O(N * D)

			printf("res %i\n",loop);
			residual.show(stdout);


	//		if (loop ==1) exit(1);

			// Xform to partial residual
			for(i=0;i<nbsubl;i++){
				coor[1] = i;
				coor[0] =0,cooralt[0]=0;
				j = k = (i ==0) ? 0 : p_maxes[i -1];
				cooralt[1] =  partition[k];
				cachedweights[0] = sum = p_weight(p_keys[partition[k]]);
				for(coor[0] = 0;coor[0] <p_data.dims[0];coor[0]++) {
						cooralt[0] = coor[0];  RD[nbsubl](coor) = residual(cooralt) * cachedweights[0];}
				for(k++;k<p_maxes[i];k++){
					sum += cachedweights[k-j] = p_weight(p_keys[partition[k]]);
					cooralt[1] =  partition[k];
					for(coor[0] = 0;coor[0] <p_data.dims[0];coor[0]++) {
						cooralt[0] = coor[0];  RD[nbsubl](coor) += residual(cooralt) * cachedweights[k-j];}
					}
				sum = 1.0f / sum;
				for(coor[0] = 0;coor[0] <p_data.dims[0];coor[0]++) {
					cooralt[0] = coor[0];
					average = RD[nbsubl](coor) * sum;
					j = k = (i ==0) ? 0 : p_maxes[i -1];
					RD[nbsubl](coor) = (double)(p_maxes[i]-k) * average;
					for(;k<p_maxes[i];k++){
						cooralt[1] =  partition[k];
						tmp = residual(cooralt) - average;
						cooralt[1] =  k-j;
						RD[i](cooralt) = tmp;
						}
					}
				}


			// partial sysmtem solve!
			for(i=0;i<=nbsubl;i++) {
		//		printf("component %i (size = )\n",i); fflush(stdout);
		//		RD[i].show(stdout);
				RD[i].Matrixleftmultiply(HH[i]);fflush(stdout);
		//		printf("becomes\n");
		//		RD[i].show(stdout);
				}


			// xform and update to complete

			for(i=0;i<nbsubl;i++){
				j = k = (i ==0) ? 0 : p_maxes[i -1];
				for(;k<p_maxes[i];k++){
					cooralt[1] =  k;
					cachedweights[k-j] = p_weight(p_keys[partition[k]]);
					cooralt[1] =  partition[k];
					for(coor[0] = cooralt[0] = 0;coor[0] <p_data.dims[0];coor[0] = ++(cooralt[0])) {
						coor[1] = k-j;
						f_out(cooralt) += RD[i](coor);
						coor[1] = i;
						f_out(cooralt) += RD[nbsubl](coor);
					}	}
				}
			}

		if (loop == 25){
			// fail to converge!
			static_warning_handdle << LFH_WARNING_MAXITE;

			}

		delete[](RD);
		delete[](HH);

		delete[](cachedweights);*/

		}



template<class O_T> const DataGrid< O_T, 3 >& Graph::render_Histogram(DataGrid< O_T, 3 > &f_out, const DataGrid<double, 1 > &image) const{
	Tuple<unsigned int,3> coor;
	coor[0] = 3;
	coor[1] = extra_rect[0] + extra_rect[2] + graph_size[0];
	coor[2] = extra_rect[1] + extra_rect[3] + graph_size[1];
	f_out.setSizes(coor);


		return(f_out);
	}
template<class O_T> const DataGrid< O_T, 3 >& Graph::render_Histogram(DataGrid< O_T, 3 > &f_out, const DataGrid<double, 2 > &image) const{
	Tuple<unsigned int,3> coor;
	coor[0] = 3;
	coor[1] = extra_rect[0] + extra_rect[2] + graph_size[0];
	coor[2] = extra_rect[1] + extra_rect[3] + graph_size[1];

	Tuple<unsigned int, 2> i_coor;


	DataGrid<double,2> trimage = image.resize_dir(1 , graph_size[0]);
	trimage *= ((double)image.dims[1]) / graph_size[0];

//	for(i_coor[0] = extra_rect[0]; i_coor[0]  < extra_rect[0] + graph_size[1];i_coor[0]){
//	}

	f_out.setSizes(coor);
	return(f_out);
	}

template <class C>
Interval<C>::Interval(C _min, C _max) : min(_min), max(_max){}
template <class C>
SetComparison Interval<C>::compare(const C & other) const{
	if (other < min){
		return(SetComparison(3+8));
	}else if (other > max){
		return(SetComparison(3+4));
	}else{
		if (min == max) return(SetComparison(0));
		else return(SetComparison(1));
	}
}

template <class C>
SetComparison Interval<C>::compareInterval(const C & q_min, const C & q_max) const{
	if (q_min  < min){
		if (q_max < max){
			if (q_max < min) return(SetComparison(3+8));
			else return(SetComparison(3));
		}else{
			if ((q_min == min)&&(q_max == max)) return(SetComparison(0));
			else return(SetComparison(3));
		}
	}else{
		if (q_max > max){
			if (q_min > max) return(SetComparison(3+4));
			else return(SetComparison(3));
		}else{
			return(SetComparison(3));
		}
	}
}

template<class I> I& OrthoRef3::applyTo(I& target)const{ // assume I is a int class, modifies last 6bit only
    switch(((r >> 3) & 7) | (target & 56)){
        case 0: case 1: case 2: case 3: case 4: case 5:
            target ^= r; // simplest of all, inherit directed and flip mismatching directions
        break;
        case 9:
        case 8: target ^= ((r&57)|((r&2) << 1)| ((r&4) >> 1));break;
        case 16: target ^= (((r&6)>>1)|((r&1) << 2));break;
        case 27:
        case 24: target ^= ((r&60)|((r&1) << 1)| ((r&2) >> 1));break;
        case 32: target ^= (((r&3)<<1)|((r&4) >> 2)); break;
        case 45:
        case 40:target ^= ((r&58)|((r&1) << 2)| ((r&4) >> 2));break;

        case 18:target ^= 48 ^ (((r&6)>>1)|((r&1) << 2));break;

        case 36:target ^= 48 ^ (((r&3)<<1)|((r&4) >> 2));break;

        case 10:target ^= 16^((r&1)|((r&2) << 1)| ((r&4) >> 1));break;
        case 11:target ^= 24^((r&1)|((r&2) << 1)| ((r&4) >> 1));break;
        case 12:target ^= 32^((r&1)|((r&2) << 1)| ((r&4) >> 1));break;
        case 13:target ^= 40^((r&1)|((r&2) << 1)| ((r&4) >> 1));break;

        case 26:target ^= 48^((r&4)|((r&1) << 1)| ((r&2) >> 1));break;
        case 25:target ^= 56^((r&4)|((r&1) << 1)| ((r&2) >> 1));break;
        case 28:target ^= 16^((r&4)|((r&1) << 1)| ((r&2) >> 1));break;
        case 29:target ^= 8^((r&4)|((r&1) << 1)| ((r&2) >> 1));break;

        case 42:target ^= 32^((r&2)|((r&1) << 2)| ((r&4) >> 2));break;
        case 43:target ^= 8^((r&2)|((r&1) << 2)| ((r&4) >> 2));break;
        case 44:target ^= 48^((r&2)|((r&1) << 2)| ((r&4) >> 2));break;
        case 41:target ^= 56^((r&2)|((r&1) << 2)| ((r&4) >> 2));break;

        case 17:target ^= 56^(((r&6)>>1)|((r&1) << 2));break;
        case 19:target ^= 24^(((r&6)>>1)|((r&1) << 2));break;
        case 20:target ^= 16^(((r&6)>>1)|((r&1) << 2));break;
        case 21:target ^= 8^(((r&6)>>1)|((r&1) << 2));break;

        case 33:target ^= 56^(((r&3)<<1)|((r&4) >> 2));break;
        case 34:target ^= 32^(((r&3)<<1)|((r&4) >> 2));break;
        case 35:target ^= 8^(((r&3)<<1)|((r&4) >> 2));break;
        case 37:target ^= 40^(((r&3)<<1)|((r&4) >> 2));break;
    }
return target;}
template <class C> OrthoRef3::operator TMatrix<C, 3,3>() const{ TMatrix<C,3,3> fout;
    ExOp::toZero(fout);
	unsigned short dir = (r >> 4) & 3;
	if (dir == 3) fprintf(stderr,"Not a direction!\n");
	else{
		ExOp::toOne(fout.data[(dir % 3)*3]); if (r & 1) ExOp::toNegative(fout.data[(dir % 3) *3]);
		ExOp::toOne(fout.data[(((dir+ ((r & 8) ? 2 : 1)) % 3)*3 ) +1]); if (r & 2) ExOp::toNegative(fout.data[(((dir+ ((r & 8) ? 2 : 1)) % 3)*3 ) +1]);
		ExOp::toOne(fout.data[(((dir+ ((r & 8) ? 1 : 2)) % 3)*3 ) +2]); if (r & 4) ExOp::toNegative(fout.data[(((dir+ ((r & 8) ? 1 : 2)) % 3)*3 ) +2]);
	}

    return fout;
}

template <class C> OrthoRef3::operator TMatrix<C, 4,4>() const{ TMatrix<C,4,4> fout;
    ExOp::toZero(fout);
	unsigned short dir = (r >> 4) & 3;
	if (dir == 3) fprintf(stderr,"Not a direction!\n");
	else{
		ExOp::toOne(fout.data[((dir % 3)<< 2)]); if (r & 1) ExOp::toNegative(fout.data[(dir % 3) << 2]);
		ExOp::toOne(fout.data[(((dir+ ((r & 8) ? 2 : 1)) % 3) << 2 ) +1]); if (r & 2) ExOp::toNegative(fout.data[(((dir+ ((r & 8) ? 2 : 1)) % 3) << 2 ) +1]);
		ExOp::toOne(fout.data[(((dir+ ((r & 8) ? 1 : 2)) % 3) << 2 ) +2]); if (r & 4) ExOp::toNegative(fout.data[(((dir+ ((r & 8) ? 1 : 2)) % 3) << 2 ) +2]);
	}
	fout.data[15] = 1.0f;

    return fout;
}



template<class C> Tuple<C,3u> OrthoRef3::multAbs(const Tuple<C,3u> &_in) const{ Tuple<C,3u> fout; // no reflexion, rotation only (preserve sign
	switch(r >> 3){
		case 1: fout[0] = _in[0]; fout[1] = _in[2]; fout[2] = _in[1];break;
		case 2: fout[0] = _in[1]; fout[1] = _in[2]; fout[2] = _in[0];break;
		case 3: fout[0] = _in[1]; fout[1] = _in[0]; fout[2] = _in[2];break;
		case 4: fout[0] = _in[2]; fout[1] = _in[0]; fout[2] = _in[1];break;
		case 5: fout[0] = _in[2]; fout[1] = _in[1]; fout[2] = _in[0];break;
		default: fout[0] = _in[0]; fout[1] = _in[1]; fout[2] = _in[2];
	}
	return fout;
}

#undef LFHTEMP
#define LFHTEMP template<unsigned int DIMS>


LFHTEMP void BBNetwork<DIMS>::sampleInit(unsigned int nb_nodes){
    nodes.setSize(nb_nodes);
    sample_scope = new BBNetwork::SampleScope;
    sample_scope->routine_ite =0;
    sample_scope->dist_buf.setSize(nb_nodes);
    sample_scope->dirty_buf.setSize(nb_nodes); ExOp::toZero(sample_scope->dirty_buf);
}

LFHTEMP void BBNetwork<DIMS>::sampleRegist(const GaussElem< Tuple<double, DIMS> > &inst){
    Tuple< unsigned int, 2> dap;
    unsigned int mind,oth;
    if (sample_scope->routine_ite < nodes.getSize()){
        dap[1] = sample_scope->routine_ite++;
        nodes[dap[1]] = inst;
        for(dap[0] =0; dap[0] < dap[1];dap[0]++)  sample_scope->proxi.insert(KeyElem<double, Tuple< unsigned int, 2> >(nodes[dap[1]].likelihoodratio_dist(nodes[dap[0]]), dap));
        sample_scope->indexmap[dap[1]] = dap[1];
        sample_scope->dirty_buf[dap[1]] = dap[1];
    }else{

        for(dap[0] =0; dap[0] < nodes.getSize();dap[0]++) {
            sample_scope->dist_buf[dap[0]] = nodes[dap[0]].likelihoodratio_dist(inst);
            if (sample_scope->dist_buf[dap[0]] < sample_scope->proxi.top().k) break;
        }
        if (dap[0] == nodes.getSize()){ // all merge are worse

            mind = sample_scope->indexmap[ sample_scope->proxi.top().d[0] ] ;
            nodes[ mind] +=   nodes[ sample_scope->indexmap[ sample_scope->proxi.top().d[1] ] ];
            nodes[ sample_scope->indexmap[ sample_scope->proxi.top().d[1] ] ] = inst;
            dap[1] = sample_scope->routine_ite++;
            for(oth=0; oth < nodes.getSize();oth++) {
                if ((oth == sample_scope->indexmap[ sample_scope->proxi.top().d[0] ]) ||(oth == sample_scope->indexmap[ sample_scope->proxi.top().d[1]] )) continue;
                dap[0] = sample_scope->dirty_buf[oth];
                sample_scope->proxi.insert(KeyElem<double, Tuple< unsigned int, 2> >(sample_scope->dist_buf[oth], dap));
            }

            sample_scope->indexmap[dap[1]] = sample_scope->indexmap[ sample_scope->proxi.top().d[1] ];
            sample_scope->indexmap.erase(sample_scope->proxi.top().d[1]);
            dap[1] = sample_scope->routine_ite++;
            sample_scope->indexmap[dap[1]] = sample_scope->indexmap[ sample_scope->proxi.top().d[0] ];
            sample_scope->indexmap.erase(sample_scope->proxi.top().d[0]);

            for(oth=0; oth < nodes.getSize();oth++) {
                if (oth == mind ) continue;
                dap[0] = sample_scope->dirty_buf[oth];
                sample_scope->proxi.insert(KeyElem<double, Tuple< unsigned int, 2> >(nodes[ mind].likelihoodratio_dist(nodes[oth]), dap));
            }

        }else{
            mind=dap[0];
            for(dap[0]++; dap[0] < nodes.getSize();dap[0]++) {
                sample_scope->dist_buf[dap[0]] = nodes[dap[0]].likelihoodratio_dist(inst);
                if (sample_scope->dist_buf[dap[0]] < sample_scope->dist_buf[mind]) mind = dap[0];
            }

            nodes[mind] += inst;
            dap[1] = sample_scope->routine_ite++;
            sample_scope->indexmap.erase(sample_scope->dirty_buf[mind]);
            sample_scope->indexmap[dap[1]] = mind;

            sample_scope->dirty_buf[mind] = dap[1];
            for(oth =0; oth < nodes.getSize();oth++) {
                dap[0] = sample_scope->dirty_buf[oth];
                if (mind == oth) continue;
                sample_scope->proxi.insert(KeyElem<double, Tuple< unsigned int, 2> >(nodes[mind].likelihoodratio_dist(nodes[oth]), dap));
            }
        }

        while(true){
            mind = sample_scope->indexmap.find(sample_scope->proxi.top().d[0]);
            if (mind != 0xFFFFFFFF){
                mind = sample_scope->indexmap.find(sample_scope->proxi.top().d[1]);
                if (mind != 0xFFFFFFFF) break;
            }
            sample_scope->proxi.pop();
        }

    }
}


LFHTEMP void BBNetwork<DIMS>::sampleFinit(){
    delete(sample_scope);
}

LFHTEMP Tuple< KeyElem<double , unsigned int > > BBNetwork<DIMS>::classmarginal(const GaussElem< Tuple<double, DIMS> > &inst, const double epsilon) const{ Tuple< KeyElem<double , unsigned int > > fout;
	double total =0.0f;
	unsigned int i;
	for(i=0;i<nodes.getSize();i++){

	}
	return fout;
}

LFHTEMP void BBNetwork<DIMS>::EMinit(){

}

LFHTEMP void BBNetwork<DIMS>::EMregist(const GaussElem< Tuple<double, DIMS> > &inst){

}

LFHTEMP void BBNetwork<DIMS>::EMfinit(){

}


#undef LFHTEMP
#define LFHTEMP template<unsigned int I>
LFHTEMP void NeuralWindow<I>::setOutputAndWindowSize(const Tuple<unsigned int, 0u > &_outputsize,const Tuple<unsigned int, I - 1> &windowsize){
	outputsize = _outputsize;
	Tuple<unsigned int, I> ndim;
	unsigned int i;
	for(i=0;i< I-1;i++) ndim[i] = windowsize[i];
	ndim[I-1] = _outputsize[_outputsize.getSize()-1];
	weights.setDims(ndim);
	weights.toRandom();
}

LFHTEMP template<class C, unsigned int IS, unsigned int OS> void NeuralWindow<I>::wrForward(DataGrid<C, OS> &f_out, const DataGrid<C, IS> &f_in) const{
	f_out.setSizes(outputsize);
}


LFHTEMP DataGrid<double, 3u> NeuralWindow<I>::imageFeatureExtraction(const DataGrid<unsigned char, 3u>& image){ DataGrid<double, 3u> fout;
	Tuple<unsigned int, 3u> coor = image.dims;
	coor[0] = 3;
	fout.setSizes(coor);
	Tuple<double, 3> rgb;
	Tuple<double, 3> hcy;
	double tmp;
	for(coor[2] = 0;coor[2]< image.dims[2];coor[2]++){
		for(coor[1] = 0;coor[1]< image.dims[1];coor[1]++){
			coor[0]=0; rgb[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			coor[0]++; rgb[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			coor[0]++; rgb[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			hcy = HCYfromRGB(rgb);
			hcy[2] = pow(hcy[2], 1.0f / 2.2f) * 2.0f - 1.0f;
			coor[0]=0; fout(coor) = hcy[2];
			hcy[1] *= sqrt(1.0f - hcy[2] * hcy[2]);
			coor[0]++; fout(coor) = cos(hcy[0]) * hcy[1];
			coor[0]++; fout(coor) = sin(hcy[0]) * hcy[1];
		}
	}
	return fout;
}

LFHTEMP DataGrid<unsigned char, 3u> NeuralWindow<I>::genImageFromFeatures(const DataGrid<double, 3u>& image){ DataGrid<unsigned char, 3u> fout;
	Tuple<unsigned int, 3u> coor = image.dims;
	coor[0] = 3;
	fout.setSizes(coor);
	Tuple<double, 3> rgb;
	Tuple<double, 3> hcy;
	double tmp;
	for(coor[2] = 0;coor[2]< image.dims[2];coor[2]++){
		for(coor[1] = 0;coor[1]< image.dims[1];coor[1]++){
			coor[0]=0; hcy[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			coor[0]++; hcy[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			coor[0]++; hcy[coor[0]] = pow(image(coor) /255.0f, 2.2f);
			hcy = HCYfromRGB(rgb);
			hcy[2] = pow(hcy[2], 1.0f / 2.2f) * 2.0f - 1.0f;
			coor[0]=0; fout(coor) = hcy[2];
			hcy[1] *= sqrt(1.0f - hcy[2] * hcy[2]);
			coor[0]++; fout(coor) = cos(hcy[0]) * hcy[1];
			coor[0]++; fout(coor) = sin(hcy[0]) * hcy[1];
		}
	}
	return fout;
}

LFHTEMP DataGrid<double, 3u> NeuralWindow<I>::imageHartWaveletTransform(const DataGrid<double, 3u>& image){
	Tuple<unsigned int, 3u> coor = image.dims;
	coor[1] >>=1;
	coor[2] >>=1;
	DataGrid<double, 3u> fout;
	fout.setSizes(coor);
	double *d = image.data;
	double *c = fout.data;
	double *mark;
	for(coor[2]=0; coor[2] < fout.dims[2]; coor[2]++){
		mark = c;
		for(coor[1]=0;coor[1] < fout.dims[1];coor[1]++){
			for(coor[0]=0; coor[0] < image.dims[0]; coor[0]++) *(c++) = *(d++);
			c-= image.dims[0];
			for(coor[0]=0; coor[0] < image.dims[0]; coor[0]++) *(c++) += *(d++);
		}
		if (image.dims[1] & 1) d+= image.dims[0];
		c = mark;
		for(coor[1]=0;coor[1] < fout.dims[1];coor[1]++){
			for(coor[0]=0; coor[0] < image.dims[0]; coor[0]++) *(c++) += *(d++);
			c-= image.dims[0];
			for(coor[0]=0; coor[0] < image.dims[0]; coor[0]++) {*c += *(d++); *(c++) /= 4.0f;}
		}
	}
	return fout;
}

#undef LFHTEMP
#define LFHTEMP template <class C>
/*LFHTEMP SparseTrianglix<C>::SparseTrianglix(){}
LFHTEMP SparseTrianglix<C>::~SparseTrianglix(){}
LFHTEMP SparseTrianglix& SparseTrianglix<C>::operator=(const SparseTrianglix&){

}*/
LFHTEMP SparseTrianglix<C>& SparseTrianglix<C>::toMemfree(){
    next_id.toMemfree(); part_id.toMemfree(); neigh.toMemfree(); attrib.toMemfree(); data.toMemfree();
    return *this;
}

LFHTEMP SparseTrianglix<C>& SparseTrianglix<C>::toMemmove(SparseTrianglix& other){
    attrib.toMemmove(other.attrib);
    data.toMemmove(other.data);
    next_id.toMemmove(other.next_id);
    part_id.toMemmove(other.part_id);
    part_data.toMemmove(other.part_data);
    neigh.toMemmove(other.neigh);
    return *this;
}


LFHTEMP void SparseTrianglix<C>::setSize(int _size){
    next_id.setSize(_size);
    part_id.setSize(_size);
    neigh.setSize(_size);
    unsigned int i;
    for(i=0;i<_size;i++) {next_id[i] = i; part_id[i] = 0xFFFFFFFF;}
    }
LFHTEMP bool SparseTrianglix<C>::hasEntry(uint32_t pos)const{return(data.find(pos) != 0xFFFFFFFF);}
LFHTEMP C& SparseTrianglix<C>::acxEntry(uint32_t pos){
    uint32_t ite = data.find(pos);
    if (ite != 0xFFFFFFFF) return data.deref(ite);
    this->addEntry(pos);
    return data[pos];
}
LFHTEMP uint32_t SparseTrianglix<C>::getFirstID(uint32_t which) const{
    if (which >= next_id.getSize()) printf("Warning, trying to find first of %X >= %i\n", which, next_id.getSize());
    if (next_id[which] == which) return which;
    which = (next_id[which] > which) ? next_id[which] : which;
    return part_id[which];
}
LFHTEMP bool SparseTrianglix<C>::isAFirstEntry(uint32_t which)const{
    if (which >= next_id.getSize()) printf("Warning, trying to find first of %X >= %i\n", which, next_id.getSize());
    if (next_id[which] == which) return true;
    if (next_id[which] < which) return false;
    return (which == part_id[next_id[which]]);
}

LFHTEMP void SparseTrianglix<C>::partition_popswap_routine(uint32_t what){
  //  printf("is popswap %i\n",  what);
    if (what < part_data.getSize()-1){
      //   ExOp::show( part_data.last());
      //   printf("set %i part ;n=%i", part_data.last()[0], next_id[part_data.last()[0]]);
         part_id[part_data.last()[0]] = what;
         part_data.pop_swap(what);
      //  ExOp::show( part_data[what]);
      //  printf("done\n");
    }else part_data.pop_back();
}

LFHTEMP bool SparseTrianglix<C>::isCyclicFID(uint16_t firstID)const{return (part_data[part_id[firstID]][1] <= part_data[part_id[firstID]][2]);}

/* finds paths between source and sink, and remove residual nodes not found in the path into disconnected partitions
 * source and sink must belong to the same component with partitionID, which must be acyclic.
 */
LFHTEMP void SparseTrianglix<C>::partition_acyclic_routine(uint16_t source, uint16_t sink, uint16_t partitionID, bool does_add_edge){
    myHashmap<uint16_t, uint16_t> visit;
    myHashmap<uint16_t, uint16_t> subvisit;
    visit[source] = source;
    uint32_t i,j,ite, oldedge,sv_ite,k;
    uint16_t fs,prev;
    uint32_t tmpkey;

    //printf("acyclic partition %i %i  (part %i)\n", source, sink, partitionID);
    if (source == sink){
        if (part_data[partitionID][1] == 1) return; // nothing to do!
      //  printf("simple sink!\n");
        next_id[source] = source;
        part_id[source] = partitionID;
        oldedge = part_data[partitionID][3];
        part_data[partitionID][3] = 0xFFFFFFFF;

        for(j=0;j< neigh[source].getSize();j++){
            if (attrib.find((source << 16) | neigh[source][j]) != 0xFFFFFFFF) continue;
            part_data.push_back();
            subvisit[source] = source;
            subvisit[neigh[source][j]] = neigh[source][j];
            //printf("new partition starting with %i\n",neigh[source][j] );
            part_data.last()[0] = neigh[source][j];
            part_data.last()[3] = (neigh[source][j]<< 16) | source;
            attrib[part_data.last()[3]] = 0xFFFFFFFF;
            attrib[(source << 16) | neigh[source][j]] = part_data[partitionID][3];
            part_data[partitionID][3] = (source << 16) | neigh[source][j];

            for(sv_ite=1; sv_ite < subvisit.getSize();sv_ite++){

                if (sv_ite > 1){
                    tmpkey = subvisit.deref_key(sv_ite-1); // dangerous
                    next_id[subvisit.deref_key(sv_ite)] = tmpkey;
                   // printf("subset n[%i] = %i\n", subvisit.deref_key(sv_ite), subvisit.deref_key(sv_ite-1));

                    if (subvisit.deref_key(sv_ite) < part_data.last()[0]) part_data.last()[0] = subvisit.deref_key(sv_ite);
                }

                for(k=0;k< neigh[subvisit.deref_key(sv_ite)].getSize();k++){
                    if (attrib.find((subvisit.deref_key(sv_ite) << 16) | neigh[subvisit.deref_key(sv_ite)][k]) != 0xFFFFFFFF) continue;
                    ite = subvisit.find(neigh[subvisit.deref_key(sv_ite)][k]);
                    if (ite == 0xFFFFFFFF){
                        tmpkey = subvisit.deref_key(sv_ite); // dangerous
                        subvisit[neigh[subvisit.deref_key(sv_ite)][k]] = tmpkey;
                    }
                }
            }
            tmpkey = subvisit.deref_key(subvisit.getSize()-1); // dangerous
            next_id[subvisit.deref_key(1)] = tmpkey;
       //     printf("subset n[%i] = %i\n", subvisit.deref_key(1), subvisit.deref_key(subvisit.getSize()-1));
            for(sv_ite=1; sv_ite< subvisit.getSize();sv_ite++) part_id[subvisit.deref_key(sv_ite)] = (subvisit.deref_key(sv_ite) == part_data.last()[0]) ? part_data.getSize()-1 : part_data.last()[0];
            part_data.last()[1] = subvisit.getSize() - 1;
            part_data.last()[2] = part_data.last()[1] - 1;
            //printf("da new partition %i has %i nodes\n", part_data.last()[0] , part_data.last()[1]);
            subvisit.toMemfree();
        }
        part_data[partitionID][0] = source;
        part_data[partitionID][1] = 1;
        part_data[partitionID][2] = 0;
    }else{
        for(i=0;i<visit.getSize();i++){
            //printf("visiting %i\n", visit.deref_key(i));
            for(j=0;j< neigh[visit.deref_key(i)].getSize();j++){
                if (attrib.find((visit.deref_key(i) << 16) | neigh[visit.deref_key(i)][j]) != 0xFFFFFFFF) continue;
                //printf("considering %i\n", neigh[visit.deref_key(i)][j]);
                ite = visit.find(neigh[visit.deref_key(i)][j]);
                if (ite == 0xFFFFFFFF){
                     // outside acyclic partition, do not bother
                    tmpkey = visit.deref_key(i);
                    visit[neigh[visit.deref_key(i)][j]] = tmpkey;
                    if (neigh[visit.deref_key(i)][j] == sink){
                        // found the path! done!
                        // partition time!
                        fs = sink;
                        i = sink;
                        oldedge = part_data[partitionID][3];
                        part_data[partitionID][3] = 0xFFFFFFFF;
                        prev = i;
                        do{
                            //printf("Spliting from %i\n",i);
                            for(j=0;j< neigh[i].getSize() ;j++){

                                if ((neigh[i][j] != visit[i])&&(prev != neigh[i][j])){
                                    if (attrib.find((i << 16) | neigh[i][j]) != 0xFFFFFFFF) continue; // that's another partition nothing to do for now
                                    // assign to new partition
                                    part_data.push_back();
                                    subvisit[i] = i;
                                    subvisit[neigh[i][j]] = neigh[i][j];
                                    //printf("New group starting with %i\n", neigh[i][j]);
                                    part_data.last()[0] = neigh[i][j];
                                    part_data.last()[3] = (neigh[i][j]<< 16) | i;
                                    //printf("edges are: %X and %X\n", (neigh[i][j]<< 16) | i, (i << 16) | neigh[i][j]);
                                    attrib[part_data.last()[3]] = 0xFFFFFFFF;

                                    attrib[(i << 16) | neigh[i][j]] = part_data[partitionID][3];
                                    part_data[partitionID][3] = (i << 16) | neigh[i][j];

                                    for(sv_ite=1; sv_ite < subvisit.getSize();sv_ite++){

                                        if (sv_ite > 1){
                                            next_id[subvisit.deref_key(sv_ite)] = subvisit.deref_key(sv_ite-1);
                     //                       printf("subset n[%i] = %i\n", subvisit.deref_key(sv_ite), subvisit.deref_key(sv_ite-1));

                                            if (subvisit.deref_key(sv_ite) < part_data.last()[0]) part_data.last()[0] = subvisit.deref_key(sv_ite);
                                        }
                                        /*if (subvisit.deref_key(sv_ite) < part_data.last()[0]){

                                            part_id[part_data.last()[0]] =

                                        }else{

                                        }*/

                                        for(k=0;k< neigh[subvisit.deref_key(sv_ite)].getSize();k++){
                                            if (attrib.find((subvisit.deref_key(sv_ite) << 16) | neigh[subvisit.deref_key(sv_ite)][k]) != 0xFFFFFFFF) continue;
                                            ite = subvisit.find(neigh[subvisit.deref_key(sv_ite)][k]);
                                            if (ite == 0xFFFFFFFF){
                                                tmpkey = subvisit.deref_key(sv_ite);
                                                subvisit[neigh[subvisit.deref_key(sv_ite)][k]] = tmpkey;
                                            }
                                        }
                                    }
                                    next_id[subvisit.deref_key(1)] = subvisit.deref_key(subvisit.getSize()-1);
                                 //   printf("subset n[%i] = %i\n", subvisit.deref_key(1), subvisit.deref_key(subvisit.getSize()-1));
                                    for(sv_ite=1; sv_ite< subvisit.getSize();sv_ite++) part_id[subvisit.deref_key(sv_ite)] = (subvisit.deref_key(sv_ite) == part_data.last()[0]) ? part_data.getSize()-1 : part_data.last()[0];
                                    part_data.last()[1] = subvisit.getSize() - 1;
                                    part_data.last()[2] = part_data.last()[1] - 1;
                                  //  printf("da new partition %i has %i nodes\n", part_data.last()[0] , part_data.last()[1]);
                                    subvisit.toMemfree();
                                }
                            }
                            //attrib.erase(i < visit[i] ? (i << 16) | visit[i]  : i | (visit[i]<< 16));
                            next_id[i] = (visit[i] == i) ? sink : visit[i];
                      //      printf("set n[%i] = %i\n", i, next_id[i]);
                            prev = i;
                            if (i == visit[i]) break;
                            i = visit[i];
                            if (i < fs) fs = i;
                        }while(true);
                        part_data[partitionID][0] = fs;
                        //printf("updated cycle to %i\n", fs);
                        part_data[partitionID][1] = 0;
                        i = sink;
                        do{
                            part_data[partitionID][1]++;
                            part_id[i] = (i == fs) ? partitionID : fs;
                            if (i == visit[i]) break;
                            next_id[i] = visit[i];
                            i = visit[i];
                        }while(true);
                        next_id[i] = sink;
                        part_data[partitionID][2] = (does_add_edge) ? part_data[partitionID][1] : part_data[partitionID][1] - 1;
                        visit.toMemfree();
                        break;
                    }
                }
            }
        }
    }
    while(oldedge != 0xFFFFFFFF){
        sv_ite = oldedge;
        oldedge = attrib[oldedge];
        //printf("metafix %X ->", sv_ite);
        k = part_id[getFirstID(sv_ite >> 16)];
        attrib[sv_ite] = part_data[k][3];
        //printf("%i\n", k);
        part_data[k][3] = sv_ite;
    }
}

LFHTEMP SparseTrianglix<C>& SparseTrianglix<C>::toZero(){
    neigh.toMemfree(); attrib.toMemfree(); data.toMemfree();
    unsigned int i;
    for(i=0;i<next_id.getSize();i++) {next_id[i] = i; part_id[i] = 0xFFFFFFFF;}
    for(i=0;i<neigh.getSize();i++) neigh[i].toMemfree();
    return *this;
}
LFHTEMP SparseTrianglix<C>& SparseTrianglix<C>::toOne(){
     attrib.toMemfree(); data.toMemfree();
    unsigned int i;
    for(i=0;i<next_id.getSize();i++) {next_id[i] = i; part_id[i] = 0xFFFFFFFF; ExOp::toOne(data[i | (i << 16)]);}
    for(i=0;i<neigh.getSize();i++) neigh[i].toMemfree();
    return *this;
}

LFHTEMP void SparseTrianglix<C>::addEntry(uint32_t pos){
    int ff = getFirstID(pos & 0xFFFF);
    int fs = getFirstID(pos >> 16);
    uint32_t i,j,prev,k;
    uint32_t ite;
    uint32_t oldedge;
    bool doesjoin = true;
   // printf("%i/%i is it %i is head\n", part_id[getFirstID(100)], part_data.getSize(), (part_id[getFirstID(100)] == 0xFFFFFFFF) ? 100 : part_data[part_id[getFirstID(100)]][0]);
    uint32_t sv_ite;
   // printf("%i/%i\n",part_id[getFirstID(62)], part_id[getFirstID(918)] );

    if (ff != fs){ // connecting distinct components
        if (part_id[ff] == 0xFFFFFFFF){ // is lone
            if (part_id[fs] == 0xFFFFFFFF){ // simple pairing!
                //printf("Easy pairing...\n");
                if (ff < fs){
                    part_id[ff] = part_data.getSize();
                    part_data.push_back();
                    part_data.last()[0] = ff;
                    part_id[fs] = ff;
                }else{
                    part_id[fs] = part_data.getSize();
                    part_data.push_back();
                    part_data.last()[0] = fs;
                    part_id[ff] = fs;
                }
                part_data.last()[1] = 2;
                part_data.last()[2] = 1;
                part_data.last()[3] = 0xFFFFFFFF;
            }else{
                if (isCyclicFID(fs)){
                    //printf("Easy mode1\n");
                    part_id[ff] = part_data.getSize();
                    part_data.push_back();
                    part_data.last()[0] = ff;
                    part_data.last()[1] = 1;
                    part_data.last()[2] = 0;
                    part_data.last()[3] = (pos >> 16) | (pos << 16);
                    attrib[part_data.last()[3]] = 0xFFFFFFFF;
                    attrib[pos] = part_data[part_id[fs]][3];
                    part_data[part_id[fs]][3] = pos;
                    doesjoin = false;
                }else{ // joining lone to acyclic component, all good
                    prev = part_id[fs];
                    if (ff < fs){
                        part_id[ff] = prev;
                        for(i=0,j = fs;i< part_data[prev][1];i++){
                            part_id[j] = ff;
                            j = next_id[j];
                        }
                        part_data[prev][0] = ff;
                    }else part_id[ff] = fs;
                    part_data[prev][1]++;
                    part_data[prev][2]++;
                    //printf("Silly mode1\n");
                }
            }
        }else{
            if (part_id[fs] == 0xFFFFFFFF){
                if (isCyclicFID(ff)){
                   // printf("Easy mode2\n");
                    part_id[fs] = part_data.getSize();
                    part_data.push_back();
                    part_data.last()[0] = fs;
                    part_data.last()[1] = 1;
                    part_data.last()[2] = 0;
                    part_data.last()[3] = pos;
                    attrib[part_data.last()[3]] = 0xFFFFFFFF;

                    attrib[(pos >> 16) | (pos << 16)] = part_data[part_id[ff]][3];
                    part_data[part_id[ff]][3] = (pos >> 16) | (pos << 16);
                    doesjoin = false;
                }else{ // joining lone to acyclic component, all good
                    prev = part_id[ff];
                    if (fs < ff){
                        part_id[fs] = prev;
                        for(i=0,j = ff;i< part_data[prev][1];i++){
                            part_id[j] = fs;
                            j = next_id[j];
                        }
                        part_data[prev][0] = fs;
                    }else part_id[fs] = ff;
                    part_data[prev][1]++;
                    part_data[prev][2]++;
                   // printf("Silly mode2\n");
                }
            }else{
                //printf("trying to join distinct components... %i to %i\n", fs, ff);
                // find if meta_path exists or not...
                Vector<uint32_t> scope;
                myHashmap<uint16_t, uint32_t> visit;

                scope.push_back(part_data[part_id[fs]][3]);
                visit[fs] = 0xFFFFFFFF;
                while(scope.getSize() !=0){
                    if (scope.last() == 0xFFFFFFFF) {scope.pop_back();continue;}

                    i = getFirstID(scope.last() & 0xFFFF);
                    ite = visit.find(i);

                    if (ite == 0xFFFFFFFF){ // not going backward
                       // printf("proc Edge %X  (%i-%i)\n", scope.last(), getFirstID(scope.last() >> 16), getFirstID(scope.last() & 0xFFFF));
                        visit[i] = scope.last();
                        if (i == ff) break;
                        scope.last() = attrib[scope.last()];
                       // printf("%i: %X at head\n", part_id[i], part_data[part_id[i]][3]);
                        if (part_data[part_id[i]][3] != 0xFFFFFFFF) scope.push_back(part_data[part_id[i]][3]);
                    }else scope.last() = attrib[scope.last()];
                }
                if (scope.getSize() !=0){ // found loop!
                  // printf("found loop!\n");
                   /* printf("there is a loop! merging all found in a single cyclic group... the nightmare\n");
                    i = ff;
                    do{
                        printf("%i--%X->", i, visit[i]);
                        i = getFirstID(visit[i] >> 16);
                    }while(visit[i] != 0xFFFFFFFF);
                    printf("%i\n", i);*/
                    doesjoin = false;

                    Vector<uint32_t> spart_id;
                    myHashmap<uint32_t, uint32_t> edges_to_kill;
                    //for(i=0;i<visit.getSize();i++) printf("visit[%X] = %X\n", visit.deref_key(i), visit.deref(i));
                    i = ff;
                    prev = (pos & 0xFFFF);
                    do{
                        spart_id.push_back(part_id[i]);
                        //printf("mmm %i %i %X %i\n", prev, (visit[i] == 0xFFFFFFFF) ? (pos >> 16) : (visit[i] & 0xFFFF), visit[i],i );
                        if (!this->isCyclicFID(i)){
                            partition_acyclic_routine(prev, (visit[i] == 0xFFFFFFFF) ? (pos >> 16) : (visit[i] & 0xFFFF), part_id[i], false);
                        }
                        if (visit[i] == 0xFFFFFFFF) break;
                        prev = (visit[i] >> 16);
                        edges_to_kill[visit[i]] =1;
                        edges_to_kill[(visit[i] >> 16) | (visit[i] << 16)] = 1;
                        i = getFirstID(prev);
                    }while(true);
                    spart_id.sort();
                    //printf("partID to merge:");
                    //spart_id.show();

                    k = part_data[spart_id[0]][0];
                    for(j=1; j< spart_id.getSize();j++) if (part_data[spart_id[j]][0] < k) k = part_data[spart_id[j]][0];

                    prev = part_data[spart_id[0]][0];
                    for(j=1; j< part_data[spart_id[0]][1];j++) {
                            part_id[prev] = k;
                            prev = next_id[prev];
                    }
                    part_id[prev] = k;
                    //printf("set N[%i] = %i\n",prev, part_data[spart_id.last()][0]);
                    next_id[prev] = part_data[spart_id.last()][0];
                    part_data[spart_id[0]][2] += spart_id.getSize();

                    oldedge = part_data[spart_id[0]][3];
                    part_data[spart_id[0]][3] = 0xFFFFFFFF;
                    int nbkill=0;
                    while(oldedge != 0xFFFFFFFF) {
                        prev = attrib[oldedge];
                        if (edges_to_kill.find(oldedge) == 0xFFFFFFFF){
                            attrib[oldedge] = part_data[spart_id[0]][3];
                            part_data[spart_id[0]][3] = oldedge;
                        }else nbkill++;
                        oldedge = prev;
                    }
                    //printf("Killed %i/%i metaedges after processing %i\n", nbkill, edges_to_kill.getSize(),spart_id[0]);


                    for(i=spart_id.getSize()-1;i>0;i--){
                        prev = part_data[spart_id[i]][0];
                        j=0;
                        while(true){
                            part_id[prev] = k;
                            j++;
                            if (j == part_data[spart_id[i]][1]) break;
                            prev = next_id[prev];
                        }
                        if (next_id[prev] != part_data[spart_id[i]][0]) {
                            printf("Warning! expected first of partition%i(s=%i) %i, found %i\n", spart_id[i], part_data[spart_id[i]][1], part_data[spart_id[i]][0], next_id[prev] );
                            /*prev = part_data[spart_id[i]][0];
                            j=0;
                            while(true){
                                j++;
                                if (j == part_data[spart_id[i]][1]) break;
                                printf("%i->%i\n", prev, next_id[prev]);
                                prev = next_id[prev];
                            }*/
                        }


                        next_id[prev] = part_data[spart_id[i-1]][0];
                        part_data[spart_id[0]][1] += part_data[spart_id[i]][1];
                        part_data[spart_id[0]][2] += part_data[spart_id[i]][2];

                        while (part_data[spart_id[i]][3] != 0xFFFFFFFF){
                            prev = attrib[part_data[spart_id[i]][3]];
                            if (edges_to_kill.find(part_data[spart_id[i]][3]) == 0xFFFFFFFF){
                                attrib[part_data[spart_id[i]][3]] = part_data[spart_id[0]][3];
                                part_data[spart_id[0]][3] = part_data[spart_id[i]][3];
                            }else nbkill++;
                            part_data[spart_id[i]][3] = prev;
                        }
                        //printf("Killed %i/%i metaedges after processing %i\n", nbkill, edges_to_kill.getSize(),spart_id[i]);

                        this->partition_popswap_routine(spart_id[i]);
                    }
                    if (nbkill != edges_to_kill.getSize()){
                        for(i=0;i<edges_to_kill.getSize();i++) printf("E%X\n", edges_to_kill.deref_key(i));
                    }
                    part_data[spart_id[0]][0] = k;
                    part_id[part_data[spart_id[0]][0]] = spart_id[0];
                    doesjoin = false;
                }else{ // found no loop!
                    //printf("there is no loop!\n");
                    if ((!this->isCyclicFID(ff))&&(!this->isCyclicFID(fs))){
                        //printf("Merge %i and %i\n", fs,ff);
                        if (fs < ff){i =fs; fs = ff; ff =i;}
                        if ((oldedge = part_data[part_id[ff]][3]) == 0xFFFFFFFF) part_data[part_id[ff]][3] = part_data[part_id[fs]][3];
                        else{
                            while(attrib[oldedge] != 0xFFFFFFFF) oldedge = attrib[oldedge];
                            attrib[oldedge] = part_data[part_id[fs]][3];
                        }
                        prev = part_id[fs];
                        for(i =0,j = part_data[prev][0];i<part_data[prev][1];i++){
                            part_id[j] = ff;
                            j = next_id[j];
                        }
                        part_data[part_id[ff]][1] += part_data[prev][1];
                        part_data[part_id[ff]][2] = part_data[part_id[ff]][1]-1;

                        this->partition_popswap_routine(prev);
                    }else{ // at least 1 cyclic, add silly link...
                        //printf("adding a link...\n");
                        prev = (pos >> 16) | (pos << 16);
                        attrib[prev] = part_data[part_id[ff]][3];
                        part_data[part_id[ff]][3] = prev;
                        attrib[pos] = part_data[part_id[fs]][3];
                        part_data[part_id[fs]][3] = pos;
                        doesjoin = false;
                    }
                }
                // two distinct components are joined... FUN!
                /*part_data[part_id[ff]][1] += part_data[part_id[fs]][1];
                part_data[part_id[ff]][2] += part_data[part_id[fs]][2] +1;
                part_data[part_id[fs]] = part_data.last();
                part_id[part_data[part_id[fs]][0]] = part_id[fs];
                part_data.pop_back();

                part_id[fs] = ff;
                for(i = next_id[fs]; fs != i; i = next_id[i]) part_id[i] = ff;
                */
            }
        }
        if (doesjoin){
            i = next_id[ff];
            next_id[ff] = next_id[fs];
            next_id[fs] = i;
        }
    }else if ((pos & 0xFFFF) != (pos >> 16)) {
        // connection within component
        if (this->isCyclicFID(ff)) part_data[part_id[ff]][2]++; // connection within cyclic component,... nothing to do
        else{
           // printf("Connecting within... %i ->", part_data.getSize());
            partition_acyclic_routine(pos & 0xFFFF, pos >> 16, part_id[ff], true);
           // printf("%i\n", part_data.getSize());
            // connecting within acylic component, need to find the one cycle made, then partition the partition...
        }
    }
    if ((pos & 0xFFFF) != (pos >> 16)) {
        neigh[pos >> 16].push_back( (uint16_t)(pos & 0xFFFF));
        neigh[pos & 0xFFFF].push_back((uint16_t)(pos >> 16));
    }
}
LFHTEMP void SparseTrianglix<C>::removeEntry(uint32_t pos){
    // TODO, cannot be done fast...
}

LFHTEMP bool SparseTrianglix<C>::wouldEdgeBeInCyclicGroup(uint32_t pos)const{
    int ff = getFirstID(pos & 0xFFFF);
    int fs = getFirstID(pos >> 16);
    uint32_t i,j,prev,k;
    uint32_t ite;
    uint32_t oldedge;
    bool doesjoin = true;
    uint32_t sv_ite;
    if (part_id[ff] == 0xFFFFFFFF) return false;
    if (ff == fs) return true;
    if (part_id[fs] == 0xFFFFFFFF) return false;
    Vector<uint32_t> scope;
    myHashmap<uint16_t, uint32_t> visit;
    scope.push_back(part_data[part_id[fs]][3]);
    visit[fs] = 0xFFFFFFFF;
    while(scope.getSize() !=0){
        if (scope.last() == 0xFFFFFFFF) {scope.pop_back();continue;}
        i = getFirstID(scope.last() & 0xFFFF);
        ite = visit.find(i);
        if (ite == 0xFFFFFFFF){ // not going backward
           // printf("proc Edge %X  (%i-%i)\n", scope.last(), getFirstID(scope.last() >> 16), getFirstID(scope.last() & 0xFFFF));
            visit[i] = scope.last();
            if (i == ff) break;
            scope.last() = attrib[scope.last()];
           // printf("%i: %X at head\n", part_id[i], part_data[part_id[i]][3]);
            if (part_data[part_id[i]][3] != 0xFFFFFFFF) scope.push_back(part_data[part_id[i]][3]);
        }else scope.last() = attrib[scope.last()];
    }
    return (scope.getSize() !=0);
}

LFHTEMP C& SparseTrianglix<C>::operator[](uint32_t pos){
    pos |= pos << 16;
    uint32_t ite = data.find(pos | (pos << 16));
    if (ite != 0xFFFFFFFF) return data.deref(ite);
    this->addEntry(pos);
    return data[pos];
}

LFHTEMP C SparseTrianglix<C>::operator[](uint32_t pos)const{
    pos |= pos << 16;
    uint32_t ite = data.find(pos | (pos << 16));
    if (ite != 0xFFFFFFFF) return data.deref(ite);
    C fout; ExOp::toZero(fout);
    return fout;
}

LFHTEMP C& SparseTrianglix<C>::operator()(const Tuple<uint32_t, 2u> &coor){
    uint32_t pos = (coor[0] < coor[1]) ? (coor[1] |  (coor[0] << 16)) : (coor[0] |  (coor[1] << 16));
    uint32_t ite = data.find(pos);
    if (ite != 0xFFFFFFFF) return data.deref(ite);
    this->addEntry(pos);
    return data[pos];
}
LFHTEMP C SparseTrianglix<C>::operator()(const Tuple<uint32_t, 2u> &coor)const{
    uint32_t pos = (coor[0] < coor[1]) ? (coor[1] |  (coor[0] << 16)) : (coor[0] |  (coor[1] << 16));
    uint32_t ite = data.find(pos);
    if (ite != 0xFFFFFFFF) return data.deref(ite);
    C input; ExOp::toZero(input);
    return input;
}
LFHTEMP bool SparseTrianglix<C>::hasEntry(const Tuple<uint32_t, 2u> &coor)const{
    uint32_t pos = (coor[0] < coor[1]) ? (coor[1] |  (coor[0] << 16)) : (coor[0] |  (coor[1] << 16));

    return(data.find(pos) != 0xFFFFFFFF);
}
LFHTEMP SparseTrianglix<C>::operator Trianglix<C>()const{
    Trianglix<C> fout; fout.setSize(this->getSize()); fout.toZero();
    uint32_t i,p;
    for(i=0;i<data.getSize();i++){
        p = (data.deref_key(i) & 0xFFFF);
        p = ((p * (p+1)) >> 1) + (data.deref_key(i) >> 16);
        fout.data[p] = data.deref(i);
    }
    return fout;
}
LFHTEMP Tuple<uint32_t, 2u> SparseTrianglix<C>::deref_key(uint32_t ite){ Tuple<uint32_t, 2u> fout;
    uint32_t pos = data.deref_key(ite);
    fout[0] = pos & 0xFFFF;
    fout[1] = pos >> 16;
    return fout;
}
LFHTEMP void SparseTrianglix<C>::wrPermutation(Tuple<uint32_t> &fout){
    fout.setSize(this->getSize());
    unsigned int cursor, i;
    PartitionIterator pp(*this);
    cursor = 0;
    if (pp.first()) do{
        for(i=0;i<pp.getSize();i++) fout[cursor++] = pp.part[i];
    }while(pp.next());
}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getPartition(int partID, Tuple<uint32_t>* opt_getEdges, bool opt_projected)const{Tuple<uint32_t> fout;
    fout.setSize(part_data[partID][1]);
    fout[0] = part_data[partID][0];
    uint32_t i,j,k;
    for(i = 1u; i < fout.getSize();i++) fout[i] = next_id[fout[i-1u]];
    if (opt_getEdges != NULL){
        opt_getEdges->setSize(part_data[partID][2]);
        if (opt_projected){
            myHashmap<uint16_t, uint32_t> visit;
            for(j=0;j<fout.getSize();j++) visit[fout[j]] = j;
            for(j=0,k=0;j<fout.getSize();j++){
                for(i=0;i<neigh[fout[j]].getSize();i++){
                    if (getFirstID(neigh[fout[j]][i]) != part_data[partID][0]) continue;
                    if (visit[neigh[fout[j]][i]] < visit[fout[j]]) (*opt_getEdges)[k++] = visit[fout[j]] | (visit[neigh[fout[j]][i]] << 16);
                }
            }
        }else{
            for(j=0,k=0;j<fout.getSize();j++){
                for(i=0;i<neigh[fout[j]].getSize();i++){
                    if (getFirstID(neigh[fout[j]][i]) != part_data[partID][0]) continue;
                    if (neigh[fout[j]][i] < fout[j]) (*opt_getEdges)[k++] = fout[j] | (neigh[fout[j]][i] << 16);
                }
            }
        }
    }
    return fout;
}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getEdgePartition(uint32_t edge, Tuple<uint32_t>* opt_getEdges, bool opt_projected)const{
    uint32_t ite = data.find(edge);
    if (ite == 0xFFFFFFFF) return this->getWouldBeCyclicPartition(edge, opt_getEdges,opt_projected); // edge does not exist, return partition if it was added
    int ff = getFirstID(edge & 0xFFFF);
    int fs = getFirstID(edge >> 16);
    if ((ff == fs)&&(isCyclicFID(ff))) return getPartition(part_id[ff], opt_getEdges, opt_projected);
    // is acyclic edge! custom partition
    Tuple<uint32_t> fout; fout.setSize(2);
    if (opt_projected){
        fout[0] = 0; fout[1] = 1;
        if (opt_getEdges) opt_getEdges->setSize(1)[0] = 1;
    }else{
        fout[0] = edge & 0xFFFF; fout[1] = edge >> 16;
        if (opt_getEdges) opt_getEdges->setSize(1)[0] = edge;
    }
return fout;}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getConnectedSet(int member)const{Tuple<uint32_t> fout;
    int partID = part_id[getFirstID(member)];
    if (partID == 0xFFFFFFFF){fout.setSize(1); fout[0] = member; return fout;}
    Vector<uint32_t> pids;
    Vector<uint32_t> edges;

    uint32_t tsize = part_data[partID][1];
    pids.push_back(partID);
    edges.push_back(0);
    uint32_t curedge;
    for(int i = 0; i < pids.getSize(); i++){
        curedge = part_data[pids[i]][3];
        while (curedge != 0xFFFFFFFF){
            if (curedge != edges[i]) {
                edges.push_back((curedge >> 16) | (curedge << 16) );
                partID = part_id[getFirstID(curedge & 0xFFFF)];
                pids.push_back(partID);
                tsize += part_data[partID][1];
            }
            curedge = attrib[curedge];
        }
    }
    fout.setSize(tsize);
    tsize =0;
    for(int i = 0; i < pids.getSize(); i++){
        fout[tsize++] = part_data[pids[i]][0];
        for(curedge = next_id[part_data[pids[i]][0]]; curedge != part_data[pids[i]][0]; curedge =next_id[curedge]) fout[tsize++] = curedge;
    }
    return fout;
}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getConnectedPartIDs(int nodeID)const{Tuple<uint32_t> fout;
    uint32_t partID = part_id[getFirstID(nodeID)];
    if (partID == 0xFFFFFFFF) {fout.setSize(1); fout[0] = 0xFFFFFFFF; return fout;}
    Vector<uint32_t> pids;
    Vector<uint32_t> edges;
    pids.push_back(partID);
    edges.push_back(0);
    uint32_t curedge;
    for(int i = 0; i < pids.getSize(); i++){
        curedge = part_data[pids[i]][3];
        while (curedge != 0xFFFFFFFF){
            if (curedge != edges[i]) {
                edges.push_back((curedge >> 16) | (curedge << 16) );
                partID = part_id[getFirstID(curedge & 0xFFFF)];
                pids.push_back(partID);
            }
            curedge = attrib[curedge];
        }
    }
    return fout.toMemmove(pids);
}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getAcyclicPath(uint32_t source, uint32_t sink)const{ Tuple<uint32_t> fout;
    if (source == sink) {fout.setSize(1); fout[0] =source; return fout;}
    myHashmap<uint16_t, uint16_t> visit;
    visit[source] = source;
    uint32_t i,j, ite, c;
    for(i=0;i<visit.getSize();i++){
        c = visit.deref_key(i);
        //printf("visiting %i\n", c);
        for(j=0;j< neigh[c].getSize();j++){
            if (attrib.find((c << 16) | neigh[c][j]) != 0xFFFFFFFF) continue;
            //printf("considering %i\n", neigh[c][j]);
            ite = visit.find(neigh[c][j]);
            if (ite == 0xFFFFFFFF){
                //visit.show();
                //printf("set set%i! %i an %i\n", j, neigh[c][j],c);
                visit[neigh[c][j]] = (uint32_t)c;
                //visit.show();
                if (neigh[c][j] == sink){
                    //printf("found target rebuild time!\n");
                    //visit.show();
                    //fflush(stdout);
                    i = 0;
                    j = sink;
                    while(visit[j] != source) {
                        //printf("F[%i] = %i visit\n", i, visit[j]); fflush(stdout);
                        //if (i > 1000) exit(0);
                        i++; j = visit[j];
                    }
                    fout.setSize(i+2);
                    j = sink;
                    for(i++;i != 0;j = visit[j]) fout[i--] = j;
                    fout[0] = source;
                    return fout;
                }
            }
        }
    }
    LFH_exit("Did not find a acyclic path, as should be expected");
    fflush(stdout);
return fout;}

LFHTEMP Tuple<uint32_t> SparseTrianglix<C>::getWouldBeCyclicPartition(uint32_t fictiveEdge, Tuple<uint32_t>* opt_getEdges, bool opt_projected, Tuple<uint32_t>* opt_outerEdges)const{Tuple<uint32_t> fout;
    Vector<uint32_t> scope;
    myHashmap<uint16_t, uint32_t> visit;
    foreach.printf("fictive %i\n", fictiveEdge); fflush(stdout);
    int ff = getFirstID(fictiveEdge & 0xFFFF);
    int fs = getFirstID(fictiveEdge >> 16);
    //printf("%i, %i ss%c%c\n", ff, fs, isCyclicFID(ff)? 'C' : 'A', isCyclicFID(fs) ? 'C' : 'A'); fflush(stdout);
    uint32_t i,j,k;
    if (!opt_projected) {printf("Warning, order is not preserved...\n");}
    if (ff == fs){
        if (isCyclicFID(ff)){
            fout = this->getPartition(part_id[ff]);
            if (opt_getEdges){
                opt_getEdges->setSize(part_data[part_id[ff]][2] +1);
                k =0;
                if (opt_projected){
                    for(j=0;j<fout.getSize();j++) visit[fout[j]] = j;
                    for(j=0;j<fout.getSize();j++){
                        for(i=0;i<neigh[fout[j]].getSize();i++){
                            if (getFirstID(neigh[fout[j]][i]) != ff) continue;
                            if (visit[neigh[fout[j]][i]] < visit[fout[j]]) (*opt_getEdges)[k++] = visit[fout[j]] | (visit[neigh[fout[j]][i]] << 16);
                        }
                    }
                    (*opt_getEdges)[k] = (visit[fictiveEdge >> 16] < visit[fictiveEdge & 0xFFFF]) ? visit[fictiveEdge & 0xFFFF] | (visit[fictiveEdge >> 16] << 16) : visit[fictiveEdge >> 16] | (visit[fictiveEdge & 0xFFFF] << 16);
                }else{
                    for(j=0;j<fout.getSize();j++){
                        for(i=0;i<neigh[fout[j]].getSize();i++){
                            if (getFirstID(neigh[fout[j]][i]) != ff) continue;
                            if (neigh[fout[j]][i] < fout[j]) (*opt_getEdges)[k++] = fout[j] | (neigh[fout[j]][i] << 16);
                        }
                    }
                    (*opt_getEdges)[k] = fictiveEdge;
                }
            }
            if (opt_outerEdges){ // Ezy-pezy
                j = 0;
                for(i = part_data[part_id[ff]][3];i != 0xFFFFFFFF; i = attrib[i]) j++;
                opt_outerEdges->setSize(j);
                j = 0;
                for(i = part_data[part_id[ff]][3];i != 0xFFFFFFFF; i = attrib[i]) (*opt_outerEdges)[j++] = i;
            }
        }else{// within acyclic component!, it must be a single loop!
            fout = this->getAcyclicPath(fictiveEdge & 0xFFFF,fictiveEdge >> 16);
            if (opt_getEdges){
                opt_getEdges->setSize(fout.getSize());
                if (opt_projected){
                    for(i=0;i<fout.getSize()-1;i++) (*opt_getEdges)[i] = (i+1) | (i << 16);
                    (*opt_getEdges)[i] = i;
                }else{
                    for(i=0;i<fout.getSize()-1;i++) (*opt_getEdges)[i] = fout[i+1] > fout[i] ? fout[i+1] | (fout[i] << 16) : fout[i] | (fout[i+1] << 16);
                    (*opt_getEdges)[i] = (fout[i] > fout[0]) ? fout[i] | (fout[0] << 16) : fout[0] | (fout[i] << 16) ;
                }
            }
            if (opt_outerEdges){
                scope.toMemfree();
                for(i=0,j=0;j<neigh[fout[i]].getSize();j++){
                    if (neigh[fout[i]][j] == fout[i+1]) continue;
                    else scope.push_back((fout[i] << 16) | neigh[fout[i]][j]);
                }
                for(i++;i<fout.getSize()-1;i++){
                    for(j=0;j<neigh[fout[i]].getSize()-1;j++){
                        if ((neigh[fout[i]][j] == fout[i+1])||(neigh[fout[i]][j] == fout[i-1])) continue;
                        else scope.push_back((fout[i] << 16) | neigh[fout[i]][j]);
                    }
                }
                for(j=0;j<neigh[fout[i]].getSize();j++){
                    if (neigh[fout[i]][j] == fout[i-1]) continue;
                    else scope.push_back((fout[i] << 16) | neigh[fout[i]][j]);
                }
                (*opt_outerEdges).toMemmove(scope);
            }
        }
    }else{
        Vector<Tuple<uint32_t> > spart_id;
        scope.push_back(part_data[part_id[fs]][3]);
        visit[fs] = 0xFFFFFFFF;
        while(scope.getSize() !=0){
            if (scope.last() == 0xFFFFFFFF) {scope.pop_back();continue;}
            i = getFirstID(scope.last() & 0xFFFF);
            k = visit.find(i);
            if (k == 0xFFFFFFFF){ // not going backward
               // printf("proc Edge %X  (%i-%i)\n", scope.last(), getFirstID(scope.last() >> 16), getFirstID(scope.last() & 0xFFFF));
                visit[i] = scope.last();
                if (i == ff) break;
                scope.last() = attrib[scope.last()];
               // printf("%i: %X at head\n", part_id[i], part_data[part_id[i]][3]);
                if (part_data[part_id[i]][3] != 0xFFFFFFFF) scope.push_back(part_data[part_id[i]][3]);
            }else scope.last() = attrib[scope.last()];
        }
        if (scope.getSize() ==0){ printf("illegal call to getWouldBeCyclicPartition, no cycle would be created %i %i! (check 'wouldEdgeBeInCyclicGroup' for true)\n", fictiveEdge & 0xFFFF, (fictiveEdge >> 16)); exit(1);}

        //for(i=0;i<visit.getSize();i++) printf("visit[%X] = %X\n", visit.deref_key(i), visit.deref(i));
        i = ff;
        j = (fictiveEdge & 0xFFFF);
        do{
            //printf("mmm %i %i %X %i\n", prev, (visit[i] == 0xFFFFFFFF) ? (pos >> 16) : (visit[i] & 0xFFFF), visit[i],i );
            if (this->isCyclicFID(i)) spart_id.push_back(this->getPartition(part_id[i]));
            else spart_id.push_back(this->getAcyclicPath(j,(visit[i] == 0xFFFFFFFF) ? (fictiveEdge >> 16) : visit[i] & 0xFFFF));
            if (visit[i] == 0xFFFFFFFF) break;
            j = (visit[i] >> 16);
            i = getFirstID(j);
        }while(true);
        for(j=0,i=0;i < spart_id.getSize();i++) j+= spart_id[i].getSize();
        fout.setSize(j);
        for(i=0,k=0;i < spart_id.getSize();i++){
            for(j=0;j<spart_id[i].getSize();j++) fout[k++] =  spart_id[i][j];
        }

        if (opt_getEdges){
            int l,m;
            if (opt_projected){
                myHashmap<uint16_t, uint16_t> mapmap;
                for(k=0,i=0,m=spart_id.getSize();i < spart_id.getSize();i++){
                    j = getFirstID(spart_id[i][0]);
                    for(l=0;l<spart_id[i].getSize();l++) mapmap[spart_id[i][l]] = k++;
                    m += (this->isCyclicFID(j)) ? part_data[part_id[j]][2] : spart_id[i].getSize()-1;
                }
                opt_getEdges->setSize(m);
                m = 0;
                for(k=0;;k++){
                    i = getFirstID(spart_id[k][0]);
                    if (!this->isCyclicFID(i)) {
                        l = mapmap[spart_id[k][0]];
                        for(j=0;j<spart_id[k].getSize()-1;j++,l++) (*opt_getEdges)[m++] = (l + 1) | (l << 16);
                    }else{
                        for(j=0;j<spart_id[k].getSize();j++){
                            for(l=0;l<neigh[spart_id[k][j]].getSize();l++){
                                if (getFirstID(neigh[spart_id[k][j]][l]) != i) continue;
                                if (mapmap[neigh[spart_id[k][j]][l]] < mapmap[spart_id[k][j]]) (*opt_getEdges)[m++] = mapmap[spart_id[k][j]] | (mapmap[neigh[spart_id[k][j]][l]] << 16);
                            }
                        }
                    }
                    if (visit[i] == 0xFFFFFFFF) break;
                    (*opt_getEdges)[m++] = mapmap[visit[i] & 0xFFFF]  | (mapmap[visit[i] >> 16] << 16);
                }
                (*opt_getEdges)[m++] = mapmap[fictiveEdge & 0xFFFF] | (mapmap[fictiveEdge >> 16] << 16);
            }else{
                for(i=0,l=spart_id.getSize();i < spart_id.getSize();i++){
                    j = getFirstID(spart_id[i][0]);
                    l += (this->isCyclicFID(j)) ? part_data[part_id[j]][2] : spart_id[i].getSize()-1;
                }
                opt_getEdges->setSize(l);
                for(k=0,m = 0;;k++){
                    i = getFirstID(spart_id[k][0]);
                    if (!this->isCyclicFID(i)) {
                        for(j=0;j<spart_id[k].getSize()-1;j++) (*opt_getEdges)[m++] = (spart_id[k][j] < spart_id[k][j+1]) ? spart_id[k][j+1] | (spart_id[k][j] << 16) : spart_id[k][j] | (spart_id[k][j+1] << 16);
                    }else{
                        for(j=0;j<spart_id[k].getSize();j++){
                            for(l=0;l<neigh[spart_id[k][j]].getSize();l++){
                                if (getFirstID(neigh[spart_id[k][j]][l]) != i) continue;
                                if (neigh[spart_id[k][j]][l] < spart_id[k][j]) (*opt_getEdges)[m++] = spart_id[k][j] | (neigh[spart_id[k][j]][l] << 16);
                            }
                        }
                    }
                    if (visit[i] == 0xFFFFFFFF) break;
                    (*opt_getEdges)[m++] = (visit[i] & 0xFFFF) > (visit[i] >> 16) ? visit[i] : (visit[i] >> 16) | (visit[i]<<16);
                }
                (*opt_getEdges)[m++] = fictiveEdge;
            }
        }
        if (opt_outerEdges) {
            scope.toMemfree();
            for(i=0;i < spart_id.getSize();i++){
                j = getFirstID(spart_id[i][0]);
                if (this->isCyclicFID(j)) {
                    for(k = part_data[part_id[j]][3]; k != 0xFFFFFFFF; k = attrib[k]){
                        if ((i>0)&&(getFirstID(spart_id[i-1][0]) == getFirstID(k & 0xFFFF))) continue;
                        if ((i < spart_id.getSize()-1)&&(getFirstID(spart_id[i+1][0]) == getFirstID(k & 0xFFFF))) continue;
                        scope.push_back(k);
                    }
                }else{
                    if (spart_id[i].getSize() == 1){
                        for(k=0,j=0;j<neigh[spart_id[i][k]].getSize();j++){
                            if ((i!=spart_id.getSize()-1)&&(getFirstID(spart_id[i+1][0]) == getFirstID(neigh[spart_id[i][k]][j]))) continue;
                            else if ((i!=0)&&(getFirstID(spart_id[i-1][0]) == getFirstID(neigh[spart_id[i][k]][j]))) continue;
                            scope.push_back((spart_id[i][k] << 16) | neigh[spart_id[i][k]][j]);
                        }
                    }else{
                        for(k=0,j=0;j<neigh[spart_id[i][k]].getSize();j++){
                            if (neigh[spart_id[i][k]][j] == spart_id[i][k+1]) continue;
                            else if ((i!=0)&&(getFirstID(spart_id[i-1][0]) == getFirstID(neigh[spart_id[i][k]][j]))) continue;
                            scope.push_back((spart_id[i][k] << 16) | neigh[spart_id[i][k]][j]);
                        }
                        for(k++;k<spart_id[i].getSize()-1;k++){
                            for(j=0;j<neigh[spart_id[i][k]].getSize();j++){
                                if ((neigh[spart_id[i][k]][j] == spart_id[i][k+1])||(neigh[spart_id[i][k]][j] == spart_id[i][k-1])) continue;
                                else scope.push_back((spart_id[i][k] << 16) | neigh[spart_id[i][k]][j]);
                            }
                        }
                        for(j=0;j<neigh[spart_id[i][k]].getSize();j++){
                            if (neigh[spart_id[i][k]][j] == spart_id[i][k-1]) continue;
                            else if ((i!=spart_id.getSize()-1)&&(getFirstID(spart_id[i+1][0]) == getFirstID(neigh[spart_id[i][k]][j]))) continue;
                            scope.push_back((spart_id[i][k] << 16) | neigh[spart_id[i][k]][j]);
                        }
                    }
                }
            }
            (*opt_outerEdges).toMemmove(scope);
        }
    }
return fout;}

// assumes edge inclusion either did not change partition for all cyclic groups
LFHTEMP void SparseTrianglix<C>::updateDiagoInverse(C* diagoinv, uint32_t edge, C target)const{
    int firstID = getFirstID(edge & 0xFFFF);
    Tuple<C> inv_row;
    Tuple<uint32_t> partit;
    C tsqr = ExOp::mkPowInt(target,-2);
    if (getFirstID(edge >> 16) != firstID){ // this is a meta edge ez!
        partit = getPartition(part_id[firstID]);
        inv_row = this->solveInPartition(edge & 0xFFFF);
        C factor = tsqr * (1.0 - diagoinv[edge & 0xFFFF]) - (1.0 - diagoinv[(edge >> 16)]);
        for(int i=0;i< partit.getSize();i++) {
            if (partit[i] == (edge & 0xFFFF)) continue;
            diagoinv[partit[i]] -= (inv_row[i] * inv_row[i]) / factor;
        }
        partit = getPartition(part_id[getFirstID(edge >> 16)]);
        inv_row = this->solveInPartition(edge >> 16);
        factor = tsqr * (1.0 - diagoinv[edge >> 16]) - (1.0 - diagoinv[(edge & 0xFFFF)]);
        for(int i=0;i< partit.getSize();i++) {
            if (partit[i] == (edge & 0xFFFF)) continue;
            diagoinv[partit[i]] -= (inv_row[i] * inv_row[i]) / factor;
        }
    }else if (isCyclicFID(firstID)) {


    }else{


    }
}
LFHTEMP template<class D> ERRCODE SparseTrianglix<C>::evalRegul(const Trianglix<C> &target, const Trianglix<D> &network,const  Tuple< Trianglix<C> > &cross_train,const  Tuple< Trianglix<C> > &cross_test, Tuple<double, 11u> &details, SparseTrianglix<C> *cross_solution, Trianglix<char> *chkmatch){
    this->setSize(target.getSize());
    int i,j,k,l;
    double norm_constant =0.0;
    Tuple<uint32_t, 4u> match; match.toZero();
    if (network.getSize() != target.getSize()) {printf("network and target size mismatch!\n"); return 1;}
    if (cross_test.getSize() != cross_train.getSize()) {printf("number of training set and test set must match!\n"); return 1;}
    if (chkmatch){
        for(k=0, i=0;i<target.getSize();i++,k++){
            for(j=0;j<i;j++) {
                if (chkmatch->data[k++] != 0) match[0]++;
            }
        }
        match[1] = (target.getSize() * (target.getSize()-1)) >> 1;
        match[1] -= match[0];
    }
    Tuple<C, 0u> eigen = target.getEigenValues();
    for(i=0;i<eigen.getSize();i++){
        if ((!ExOp::isValid(eigen[i]))||(eigen[i] <=0)) break;
    }
    if (i != eigen.getSize()){
        printf("Warning, got non-positive eigen values in input covariance! \n");
    }



    Tuple<uint32_t> part;
    SparseTrianglix<C> *trainsparse = (cross_train.getSize() == 0) ? NULL : new SparseTrianglix<C>[cross_train.getSize()];

    // Init to diagonal
    for(i=0;i< target.getSize();i++) this->data[i | (i<<16)] = 1.0 / target[i];
    for(j=0;j< cross_train.getSize();j++) {
        trainsparse[j].setSize(target.getSize());
        for(i=0;i< target.getSize();i++) trainsparse[j].data[i | (i<<16)] = 1.0 / cross_train[j][i];
    }

    for(k=0, i=0;i<target.getSize();i++,k++){
        for(j=0;j<i;j++) {
            if (network.data[k++] != (D)0) {
                this->addEntry( (j << 16) | i);
                for(l=0;l< cross_train.getSize();l++) trainsparse[l].addEntry( (j << 16) | i);
    }   }   }

    for(i=0;i< part_data.getSize();i++){
        if (part_data[i][1] <= part_data[i][2]) {
            this->solveConstainedML(i,target);
            for(j=0;j< cross_train.getSize();j++) trainsparse[j].solveConstainedML(i, cross_train[j]);
        }
    }
    this->solveAllAcyclicEdges(target);
    for(j=0;j< cross_train.getSize();j++) trainsparse[j].solveAllAcyclicEdges(cross_train[j]);

    details.toZero();
    //details[1] = this->log_determinant();
    for(i=0;i<data.getSize();i++){
        j = data.deref_key(i);
        j = ((j & 0xFFFF) < (j >> 16)) ? (j & 0xFFFF) + (((j >> 16) *((j >> 16)+1))>>1) : (j >> 16) + (((j & 0xFFFF) *((j & 0xFFFF)+1))>>1);
        details[2] += data.deref(i) * target.data[j];
    }

    if (cross_solution){
        for(j=0;j< cross_train.getSize();j++) cross_solution[j].toMemmove(trainsparse[j]);
    }
return 0;}

LFHTEMP ERRCODE SparseTrianglix<C>::searchRegul(const Trianglix<C> target, const Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > > &details, double lambda, Trianglix<char> *chkmatch, uint32_t maxcyclic){
    this->setSize(target.getSize());
    uint32_t totsize = target.getSize();
    int partID;
    Tuple<uint32_t> partit[4];
    uint32_t init_t = clock();

    KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > details_input;

    HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > > heap_cyclic_merges;
    KeyElem<double, KeyElem<Tuple<uint32_t, 2u>, Tuple<double> > > cycinput;
    Tuple<double> theotherone;



    HeapTree< KeyElem<double, uint32_t> > heap_acyclic_merges;
    HeapTree< KeyElem<double, uint32_t> > heuristic_merges;
    Tuple<uint32_t> dirty_buf; dirty_buf.setSize(totsize); dirty_buf.toZero();
    uint32_t dirtyID = 0;
    Tuple< KeyElem<double, uint32_t> > current_rep; current_rep.setSize(totsize);
    Tuple<uint32_t, 2u> coor; // (coor[0] > coor[1])
    Tuple<uint32_t> bestedge_buf; bestedge_buf.setSize(totsize);

    Tuple< myHashmap<uint32_t, C> > cross_P; cross_P.setSize(cross_pairs.getSize()); // lambda finder scope


    uint32_t bestedge_size;
    KeyElem<double, uint32_t > input;


    // step 1: normalize target;
    Trianglix<C> ntrgt; ntrgt.setSize(totsize);




    int i,j,k,l;

    double norm_constant =0.0;
    Tuple<uint32_t, 4u> match; match.toZero();
    if (chkmatch){
        for(k=0, i=0;i<target.getSize();i++,k++){
            for(j=0;j<i;j++) {
                if (chkmatch->data[k++] != 0) match[0]++;
            }
        }
        match[1] = (target.getSize() * (target.getSize()-1)) >> 1;
        match[1] -= match[0];
    }

    Tuple< Trianglix<C> > locross_pairs = cross_pairs;
    double tmp;
    for(k=0, i=0;i<target.getSize();i++){
        for(j=0;j<i;j++) {
            ntrgt.data[k] = target.data[k];
            tmp = ExOp::mkPowInvInt(target[j] * target[i], -2);
            for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] *= tmp;
            ntrgt.data[k++] *= tmp;
        }
        norm_constant += log(target[i]);
        for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] /= target[i];
        ExOp::toOne(ntrgt.data[k++]);
    }
    Tuple<C> invup;


/*
    Tuple<C, 0u> eigen = ntrgt.getEigenValues();




    // makes sure matrix is not singular...
    for(k=0;k<100;k++){
        eigen = ntrgt.getEigenValues();
        j = 0;
        for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
        if  (eigen[j] >0.0) break;
        if (eigen[j] > -0.0000001) tmp= 0.999999;
        else tmp = 1.0 / (1.0 - eigen[j]);
        if (k == 0) printf("warning: singular covariance, adding prior\n");
        for(i=0,l=0;i<ntrgt.getSize();i++,l++){
            for(j=0;j<i;j++) ntrgt.data[l++] *= tmp;
        }
    }
    if (k != 0){
        ntrgt.show();
        printf("final eigenvalues:");
        eigen.show();
        if (k ==100) {printf("failed to make the matrix non-singular... really?\n"); return 1;}
    }

    for(k=0;k<100;k++){
        for(partID=0;partID< locross_pairs.getSize();partID++){
            eigen = locross_pairs[partID].getEigenValues();
            j = 0;
            for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
            if ((partID == 0)||(eigen[j] < tmp)) tmp = eigen[j];
        }
        if (tmp > 0.0) break;
        if (k == 0) printf("warning: singular covariance for cross validation, adding prior\n");
        if (tmp > -0.0000001) tmp = 1.0 - 0.000001 * k;
        else tmp = 1.0 / (1.0 - tmp * (k+1));

        for(i=0,l=0;i<ntrgt.getSize();i++,l++){
            for(j=0;j<i;j++,l++) {
                for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[l] *= tmp;
            }
        }
    }
    if (k != 0){
        printf("final eigenvalues:");
        for(partID=0;partID< locross_pairs.getSize();partID++){
            eigen = locross_pairs[partID].getEigenValues();
            eigen.show();
        }
        if (k ==100) {printf("failed to make cthe matrix non-singular for cross validation... really?\n"); return 1;}
    }*/



    Trianglix<C> trgtinv = ntrgt.mkInverse();
    this->toOne(); // initialize P to one

    Tuple<C> cross_deter; cross_deter.setSize(locross_pairs.getSize()); ExOp::toZero(cross_deter);
    Tuple<C> cross_trace; cross_trace.setSize(locross_pairs.getSize()); ExOp::toZero(cross_trace);

    Tuple<C> cross_deter_ac; cross_deter_ac.setSize(locross_pairs.getSize()); ExOp::toZero(cross_deter_ac);
    Tuple<C> cross_trace_ac; cross_trace_ac.setSize(locross_pairs.getSize()); ExOp::toZero(cross_trace_ac);

    myHashmap< uint32_t, Tuple<C> > cross_acyclic_comp;
    Tuple<C> cross_acyclic;

    cross_deter += norm_constant; // includes normalization for comparison
    C deter;
    C trace;
    C acyclic;
    ExOp::toZero(deter); ExOp::toZero(acyclic); ExOp::toOne(trace); trace *= target.getSize();

    WeightElem<Tuple<double ,3u> ,2u> crossstat;
    Tuple<double ,3u> crossstat_input;


    Tuple<double ,2u> tmpdoubles;
    ExOp::toZero(details_input.k);
    ExOp::toZero(details_input.d);
    //details_input.d[2] -= target.getSize() * locross_pairs.getSize();
    //details_input.d[1] += norm_constant * locross_pairs.getSize();
    for(j= 0; j < target.getSize(); j++){
        k = j | (j << 16);
        for(i= 0; i < locross_pairs.getSize(); i+=2){
            crossstat_input[0] = log(cross_P[i][k] = ExOp::mkInverse(locross_pairs[i][j]) );
            cross_deter[i] += crossstat_input[0];
            crossstat_input[1] = log(cross_P[i|1][k] = ExOp::mkInverse(locross_pairs[i|1][j]) );
            cross_deter[i|1] += crossstat_input[1];
      //      printf("logdet: %e %e\n", crossstat_input[0], crossstat_input[1]);
            details_input.d[1] += crossstat_input[0] + crossstat_input[1];

            cross_trace[i]   += (crossstat_input[0] = cross_P[i][k] * locross_pairs[i|1][j] );
            cross_trace[i|1] += (crossstat_input[1] = cross_P[i|1][k] * locross_pairs[i][j] );
            details_input.d[2] += crossstat_input[0] + crossstat_input[1];
      //      printf("trace: %e %e\n", crossstat_input[0], crossstat_input[1]);

        }
    }

    crossstat.toZero();

    for(i= 0; i < locross_pairs.getSize(); i++){
        crossstat_input[0] = cross_deter[i];
        crossstat_input[1] = cross_trace[i];
        crossstat_input[2] = cross_deter[i] - cross_trace[i];
        crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
    }

    details_input.d[4] = sqrt(crossstat.getVar()[2]);
    details.push_back(details_input);
    details_input.d[3] = details_input.d[2] - details_input.d[1];

 //   crossstat_input = crossstat.getMean();
 //   printf("0\t%e-%e (%e-%e  ", deter, trace, crossstat_input[0], crossstat_input[1]);
 //   crossstat_input = crossstat.getVar();
 //   printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

    for(k = 1; k < totsize; k++){
        for(input.d =k; input.d < (k << 16); input.d+= 0x10000){
            input.k = ntrgt.data[(input.d >> 16) + ((k * (k+1)) >> 1)];
            input.k = log(1.0f - input.k * input.k);
            heap_acyclic_merges.insert(input);
        }
    }
    C tmp2[6];

    while((!heap_acyclic_merges.isEmpty())||(!heap_cyclic_merges.isEmpty())){


        if ((!heap_cyclic_merges.isEmpty())&&((heap_acyclic_merges.isEmpty())||(heap_acyclic_merges.top().k > heap_cyclic_merges.top().k))){
            //printf("Current\n");
            //((Trianglix<C>)*this).show();
            //printf("Difference\n");
            //(this->mkInverse() - ntrgt).show();

            details_input.k[0] = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
            details_input.k[1] = heap_cyclic_merges.top().d.k[0] >> 16;
            details_input.d[0] -= heap_cyclic_merges.top().k;


      //      printf("Cyclic %i,%i is %e\n", heap_cyclic_merges.top().d.k[0] & 0xFFFF, heap_cyclic_merges.top().d.k[0] >> 16, -heap_cyclic_merges.top().k);
            partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0], &partit[1], true, &partit[2]);

            details_input.d[9] = partit[0].getSize();
            for(i=0;i<partit[0].getSize();i++){
                if (this->isCyclicFID(this->getFirstID(partit[0][i]))) continue;
                for(j=0;j<partit[0].getSize();j++){
                    l = (partit[0][j] < partit[0][i]) ? partit[0][i] | (partit[0][j] << 16) : partit[0][j] | (partit[0][i] << 16);
                    if ((k = cross_acyclic_comp.find(l)) != 0xFFFFFFFF){
                        for(l=0;l<locross_pairs.getSize();l++){
                            cross_trace_ac[l] -= cross_acyclic_comp.deref(k)[l];
                            cross_deter_ac[l] -= cross_acyclic_comp.deref(k)[l+locross_pairs.getSize()];
                            cross_deter[l] += cross_acyclic_comp.deref(k)[l+locross_pairs.getSize()];
                            details_input.d[2] -= cross_acyclic_comp.deref(k)[l];
                        }
                        cross_acyclic_comp.deref(k).toMemfree();
                        cross_acyclic_comp.erase_from_iterator(k);
                        k = (partit[0][j] < partit[0][i]) ? partit[0][j] + ((partit[0][i] * (partit[0][i]+1)) >> 1) : partit[0][i] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        acyclic += log(1.0 - ntrgt.data[k] * ntrgt.data[k]);
                        deter -= log(1.0 - ntrgt.data[k] * ntrgt.data[k]);

                    }
                }
            }

            this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],cycinput.d.d, ntrgt, details_input.d[5], false,false);//,true);,data

            //printf("got %i acyclic to fix too!\n", partit[3].getSize());
            // replace by edgeID and triindex

            partit[3].setSize(partit[1].getSize() + partit[0].getSize());
            for(i=0;i<partit[1].getSize();i++) {
                partit[3][i] = (partit[0][partit[1][i] & 0xFFFF] < partit[0][partit[1][i] >> 16]) ? partit[0][partit[1][i] >> 16] | (partit[0][partit[1][i] & 0xFFFF] << 16) : partit[0][partit[1][i] & 0xFFFF] | (partit[0][partit[1][i] >> 16] << 16);
            }
            for(i=0;i<partit[0].getSize();i++) partit[3][i + partit[1].getSize()] = partit[0][i] | (partit[0][i] << 16);
            partit[1].setSize(partit[3].getSize());
            for(i=0;i<partit[3].getSize();i++) partit[1][i] = (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

//            for(i=0;i< partit[3].getSize();i++) partit[1][i+partit[3].getSize()] =  (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

            // adding edge in crosses!

            details_input.d[6] = 0.0;
            for(l =0; l< locross_pairs.getSize();l+=2){
                this->data.toMemswap(cross_P[l]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],cycinput.d.d, locross_pairs[l], crossstat_input[1], true,false);//,cross_P[l|1],true
                this->data.toMemswap(cross_P[l]);
                cross_deter[l] += crossstat_input[0];
                details_input.d[1] += crossstat_input[0];
                details_input.d[6] += crossstat_input[1];

                this->data.toMemswap(cross_P[l|1]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],theotherone, locross_pairs[l|1], crossstat_input[1],true,false);//,cross_P[l|1],true
                this->data.toMemswap(cross_P[l|1]);

                cross_deter[l|1] += crossstat_input[0];
                details_input.d[1] += crossstat_input[0];
                details_input.d[6] += crossstat_input[1];

                crossstat_input.toZero();
                //tmpdoubles.toZero();
                if (partit[3].getSize() != theotherone.getSize()) printf("adding edge problem! theotherone is has wrong length! %i != %i\n", partit[3].getSize(),theotherone.getSize() );
                if (partit[3].getSize() != cycinput.d.d.getSize()) printf("adding edge problem! cycinput.d.d is has wrong length! %i != %i\n", partit[3].getSize(),cycinput.d.d.getSize() );
                for(i=0;i< partit[3].getSize() - partit[0].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];


                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }
                crossstat_input *= 2.0f; // sym matrix right?
                //tmpdoubles *= 2.0f;

                for(;i<partit[3].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }

                cross_trace[l] += crossstat_input[0];
                cross_trace[l|1] += crossstat_input[1];
                //cross_deter[l] += tmpdoubles[0];
                //cross_deter[l|1] += tmpdoubles[1];

                details_input.d[2] += crossstat_input[0] + crossstat_input[1];
                //details_input.d[1] += tmpdoubles[1] + tmpdoubles[0];
            }

            crossstat_input[0] = 0.0f;
            for(i=0;i<partit[3].getSize() - partit[0].getSize();i++){
                if ((j = data.find(partit[3][i])) != 0xFFFFFFFF){
                    crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data.deref(j)) * ntrgt.data[partit[1][i]];
                    data.deref(j) = heap_cyclic_merges.top().d.d[i];
                }else{
          //          printf("Adding that entry! P[%i,%i] = %e\n",partit[3][i] & 0xFFFF, partit[3][i] >> 16, heap_cyclic_merges.top().d.d[i] );
                    this->addEntry(partit[3][i]);
                    crossstat_input[0] += heap_cyclic_merges.top().d.d[i] * ntrgt.data[partit[1][i]];
                    data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
                }
            }
            crossstat_input[0] *= 2.0f;
            for(;i<partit[3].getSize();i++){
                crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data[partit[3][i]]) * ntrgt.data[partit[1][i]];
                data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
            }

            trace += crossstat_input[0];
            deter -= heap_cyclic_merges.top().k - crossstat_input[0];

            // get a fresh taste within cycle group! at last
       //     printf("Heuristic quandidates Cycle\n"); fflush(stdout);
            // finding new candidates within cycle or leading to new cycles

            Vector<uint32_t> edges;
            i=0;
            k = part_id[getFirstID(partit[0][0])];
            j = part_data[k][0];
            do{bestedge_buf[i++] = j; j = next_id[j];
            }while(j != part_data[k][0]);
            partit[1].toMemfree();
            edges.push_back(part_data[k][3]);
            for(;edges[0] != 0xFFFFFFFF;edges[0] = attrib[ edges[0]]){
                partit[1].push_back(i);
                k = part_id[getFirstID(edges[0] & 0xFFFF)];
                edges.push_back(part_data[k][3]);
                j = part_data[k][0];
                do{bestedge_buf[i++] = j; j = next_id[j];
                }while(j != part_data[k][0]);
                k = part_id[getFirstID(partit[0][0])];
                while(true){
                    if (edges.last() == 0xFFFFFFFF) {
                        edges.pop_back();
                        if ((k = edges.getSize()) == 1) break;
                        k = part_id[getFirstID(edges[k-2] >>16 )];
                        edges.last() = attrib[edges.last()];
                    }else{
                        j = part_id[getFirstID(edges.last() & 0xFFFF)];
                        if (j != k){
                            k = part_data[j][0];
                            do{bestedge_buf[i++] = k; k = next_id[k];
                            }while(k != part_data[j][0]);
                            k = part_id[getFirstID(edges.last() >> 16)];
                            edges.push_back(part_data[j][3]);
                        }else edges.last() = attrib[edges.last()];
                    }
                }
            }
            bestedge_size = i;
            edges.toMemfree();
            if (partit[1].getSize() == 0) partit[1].push_back(bestedge_size);
            for(i=0;i < partit[1][0];i++){
                for(j=i+1;j<bestedge_size;j++){
                    l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                    if (this->hasEntry(l)) continue;
                    heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                }
            }
            for(k=1;k<partit[1].getSize();k++){
                for(;i < partit[1][k];i++){
                    for(j=partit[1][k];j<bestedge_size;j++){
                        l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                        if (this->hasEntry(l)) continue;
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                    }
                }
            }

            if (bestedge_size > heuristic_merges.getSize()) bestedge_size = heuristic_merges.getSize();
            for(i=0;i<bestedge_size;) bestedge_buf[i++] = heuristic_merges.pop().d;
            heuristic_merges.toMemfree();
            dirtyID++;
            for(i=0;i<partit[0].getSize();i++) dirty_buf[partit[0][i]] = dirtyID;

            if (chkmatch){
                i = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
                j = heap_cyclic_merges.top().d.k[0] >> 16;
                i = (i < j) ? i + ((j * (j+1)) >> 1) : j + ((i * (i+1)) >> 1);
                match[ (chkmatch->data[i] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }
            // makes sure next cyclic merge is not dirty, and not of a too large size!
            for(heap_cyclic_merges.pop();!heap_cyclic_merges.isEmpty();heap_cyclic_merges.pop()){
                partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0]);
                if ((maxcyclic != 0)&& (partit[0].getSize() > maxcyclic)) continue;
                for(i=0;i<partit[0].getSize();i++){
                    if (dirty_buf[partit[0][i]] != heap_cyclic_merges.top().d.k[1]) break;
                }
                if (i == partit[0].getSize()) break;
            }
            //printf("Current\n");((Trianglix<C>)*this).show();printf("Difference\n");(this->mkInverse() - ntrgt).show();
        }else{
            k = heap_acyclic_merges.top().d;
            //if (wouldEdgeBeInCyclicGroup(k)) LFH_exit("ImpPPPPosible");
            details_input.d[9] = 0;


            details_input.k[0] = k & 0xFFFF;
            details_input.k[1] = k >> 16;
            details_input.d[0] -= heap_acyclic_merges.top().k;

            coor[0] = k & 0xFFFF;
            coor[1] = (k >> 16);
       //     printf("Acyclic %i,%i is %e\n", coor[0], coor[1], -heap_acyclic_merges.top().k);

            // find best candidates for *new* cyclic edges!!
            partit[0] = this->getConnectedSet(coor[0]);
            partit[1] = this->getConnectedSet(coor[1]);
      //      printf("Heuristic quandidates Acycle\n"); fflush(stdout);
            if ((partit[0].getSize() > 1)&&(partit[1].getSize() > 1)){
                for(i=0;i<partit[0].getSize();i++){
                    for(j=0;j<partit[1].getSize();j++){
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(target(partit[0][i], partit[1][j])), partit[0][i] | (partit[1][j] << 16)));
                    }
                }
                bestedge_size = partit[0].getSize() + partit[1].getSize() - 1;
                for(i=0;i<bestedge_size;){
                    k = (heuristic_merges.top().d >> 16);
                    j = (heuristic_merges.top().d & 0xFFFF);
                    if (k > j){
                        if ((coor[0] != k)||(coor[1] != j)) bestedge_buf[i++] = (heuristic_merges.top().d >> 16) | (heuristic_merges.top().d << 16);
                    }else{
                        if ((coor[0] != j)||(coor[1] != k)) bestedge_buf[i++] = heuristic_merges.top().d;
                    }
                    heuristic_merges.pop();
                }
                heuristic_merges.toMemfree();
            }else if (partit[0].getSize() == 1){
                bestedge_size = partit[1].getSize() -1;
                for(i=0,j=0;i<partit[1].getSize();i++){
                    if (partit[1][i] == coor[1]) j++;
                    else bestedge_buf[i-j] = (partit[1][i] < partit[0][0]) ? partit[0][0] | (partit[1][i] << 16) : partit[1][i] | (partit[0][0] << 16);
                }
            }else{
                bestedge_size = partit[0].getSize() -1;
                for(i=0,j=0;i<partit[0].getSize();i++){
                    if (partit[0][i] == coor[0]) j++;
                    else bestedge_buf[i-j] = (partit[0][i] < partit[1][0]) ? partit[1][0] | (partit[0][i] << 16) : partit[0][i] | (partit[1][0] << 16);
                }
            }
            k = coor[1] + ((coor[0] * (coor[0] + 1)>>1));
            // does not need to update !!invP!! lols

            //tmp2[3] = ntrgt.data[k];
            //tmp2[2] = tmp2[3] * tmp2[3];
            //tmp2[1] = -tmp2[2] +1.0f;
            //tmp2[0] = tmp2[2] / tmp2[1];

            /*
            partit = this->getConnectedSet(coor[0]);
            invup = this->solveInPartition(partit, coor[0]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];

            partit = this->getConnectedSet(coor[1]);
            invup = this->solveInPartition(partit,coor[1]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];*/

            //(*this)[coor[0]] += tmp2[0];
            //(*/this)[coor[1]] += tmp2[0];
            //(*this)(coor) = -ntrgt.data[k] / tmp2[1]; // ADDING THAT EDGE!
            (*this)(coor) = 0.0f; // ADDING THAT EDGE!
            //trace += heap_acyclic_merges.top().k - tmp2[0];
            //deter -= heap_acyclic_merges.top().k;
            acyclic -= heap_acyclic_merges.top().k;

            cross_acyclic.setSize(locross_pairs.getSize()*2);
            for(i =0; i< locross_pairs.getSize();i+= 2){
                tmp2[0] = locross_pairs[i][coor[0]] * locross_pairs[i][coor[1]];
                tmp2[1] = locross_pairs[i].data[k] * locross_pairs[i].data[k];
                cross_acyclic[locross_pairs.getSize() + i] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                details_input.d[1] += cross_acyclic[locross_pairs.getSize() + i];
                cross_deter_ac[i] += cross_acyclic[locross_pairs.getSize() + i];
                tmp2[5] = tmp2[0] - tmp2[1];
                tmp2[4] = tmp2[1] / tmp2[5];
                tmp2[0] = locross_pairs[i|1][coor[0]] * locross_pairs[i|1][coor[1]];
                tmp2[1] = locross_pairs[i|1].data[k] * locross_pairs[i|1].data[k];
                cross_acyclic[locross_pairs.getSize() + (i|1)] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                details_input.d[1] += cross_acyclic[locross_pairs.getSize() + (i|1)];
                cross_deter_ac[i|1] += cross_acyclic[locross_pairs.getSize() + (i|1)] ;

//                printf("logdet accyclic: %e %e\n", cross_acyclic[locross_pairs.getSize() + i], cross_acyclic[locross_pairs.getSize() + (i|1)]);

                tmp2[3] = tmp2[0] - tmp2[1];
                tmp2[2] = tmp2[1] / tmp2[3];

                tmp2[0] = tmp2[4] / locross_pairs[i][coor[0]];
                cross_acyclic[i] = tmp2[0] * locross_pairs[i|1][coor[0]];
                //cross_P[i][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[4] / locross_pairs[i][coor[1]];
                cross_acyclic[i] += tmp2[0] * locross_pairs[i|1][coor[1]];
                //cross_P[i][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -locross_pairs[i].data[k] / tmp2[5];
                cross_acyclic[i] += tmp2[0] * locross_pairs[i|1].data[k] *2.0;
                //cross_P[i][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i] += cross_acyclic[i];



                tmp2[0] = tmp2[2] / locross_pairs[i|1][coor[0]];
                cross_acyclic[i|1] = tmp2[0] * locross_pairs[i][coor[0]];
                //cross_P[i|1][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[2] / locross_pairs[i|1][coor[1]];
                cross_acyclic[i|1] += tmp2[0] * locross_pairs[i][coor[1]];
                //cross_P[i|1][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -locross_pairs[i|1].data[k] / tmp2[3];
                cross_acyclic[i|1] += tmp2[0] * locross_pairs[i].data[k] *2.0;
                //cross_P[i|1][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i|1] += cross_acyclic[i|1];
                details_input.d[2] += cross_acyclic[i] + cross_acyclic[i|1];

 //               printf("accyclic trace: %e %e\n", cross_acyclic[i], cross_acyclic[i|1]);
            }
            cross_acyclic_comp[heap_acyclic_merges.top().d].toMemmove(cross_acyclic);

            if (chkmatch){
                match[ (chkmatch->data[coor[1] + ((coor[0] * (coor[0] + 1)>>1))] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }else {details_input.d[7] = details_input.d[8] = 0.0;}

            for(heap_acyclic_merges.pop();!heap_acyclic_merges.isEmpty();heap_acyclic_merges.pop()){
                if (!wouldEdgeBeInCyclicGroup(heap_acyclic_merges.top().d)) break;
            }


        }
        details_input.d[3] += details_input.d[1] -  details_input.d[2];

        crossstat.toZero();
        for(i= 0; i < locross_pairs.getSize(); i++){
            crossstat_input[0] = cross_deter[i] + cross_deter_ac[i];
            crossstat_input[1] = cross_trace[i] + cross_trace_ac[i];
            crossstat_input[2] = crossstat_input[0] - crossstat_input[1];
            crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
        }

        details_input.d[4] = sqrt(crossstat.getVar()[2]);
        details_input.d[10] = clock() - init_t;
        details.push_back(details_input);
        details_input.d[3] = details_input.d[2] - details_input.d[1];
    //    printf("Start fill cyclic quandidates\n"); fflush(stdout);
        cycinput.d.k[1] = dirtyID;
        for(i=0;i<bestedge_size;i++){
            cycinput.d.k[0] = bestedge_buf[i];
            cycinput.k = -this->solveWithFictiveEdge2017(bestedge_buf[i],cycinput.d.d, ntrgt,crossstat_input[0] ,false,false);// ,data,false
            heap_cyclic_merges <<= cycinput;
        }
    //    printf("End fill cyclic quandidates %i \n", heap_cyclic_merges.getSize()); fflush(stdout);
/*  checking zone!



        // check trace is behaving!

        crossstat_input[0] =0.0;
        for(i=0; i< data.getSize();i++){
            j = data.deref_key(i); k =(j & 0xFFFF);
            k = (j >> 16) + ((k * (k+1)) >> 1);
            crossstat_input[0] += data.deref(i) * ntrgt.data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        this->wrDeterminant(crossstat_input[1]);
        printf("trace equal: %e ?= %e\n", trace, crossstat_input[0]);
        printf("deter equal: %e ?= %e\n", deter, log(crossstat_input[1]));

        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[0].getSize();i++){
                j = cross_P[0].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[0].deref(i) * locross_pairs[1].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[0], crossstat_input[0]);
        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[1].getSize();i++){
                j = cross_P[1].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[1].deref(i) * locross_pairs[0].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[1], crossstat_input[0]);

        printf("Acyclic traces: %e & %e\n", cross_trace_ac[0], cross_trace_ac[1]);
        printf("Acyclic determinant: %e & %e\n", cross_deter_ac[0], cross_deter_ac[1]);


        this->data.toMemswap(cross_P[0]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[0]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[0]);
        this->data.toMemswap(cross_P[1]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[1]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[1]);


        crossstat_input = crossstat.getMean();
        printf("%i\t%e-%e (%e-%e  ", data.getSize() - totsize, deter, trace, crossstat_input[0], crossstat_input[1]);
        crossstat_input = crossstat.getVar();
        printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

        printf("%i: %e\n", data.getSize() - totsize, deter + acyclic - trace);
        crossstat_input.toZero();
        for(i=0; i< cross_P.getSize();i++){
            crossstat_input[0] += cross_deter_ac[i] - cross_trace_ac[i];
            crossstat_input[1] += cross_deter[i] - cross_trace[i];
        }
        printf("vs: %e\n", (crossstat_input[0]  + crossstat_input[1]) / cross_P.getSize());
        */
        printf("%i edges searched!\n",(data.getSize() - totsize));
        if ((data.getSize() - totsize) > lambda) break;
    }

    // formally add acyclic edges
    for(i=0;i<part_data.getSize();i++){
        if (part_data[i][1] > part_data[i][2]){
            partit[0] = this->getPartition(i);
            for(j=0;j<partit[0].getSize();j++){
                for(k=0;k<neigh[partit[0][j]].getSize();k++){
                    if ((getFirstID(neigh[partit[0][j]][k]) != partit[0][0])||(partit[0][j] < neigh[partit[0][j]][k])){
                        l = neigh[partit[0][j]][k] > partit[0][j] ? partit[0][j] + ((neigh[partit[0][j]][k] * (neigh[partit[0][j]][k]+1)) >> 1) :  neigh[partit[0][j]][k] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        tmp2[3] = ntrgt.data[l];
                        tmp2[2] = tmp2[3] * tmp2[3];
                        tmp2[1] = -tmp2[2] +1.0f;
                        tmp2[0] = tmp2[2] / tmp2[1];
                        (*this)[partit[0][j]] += tmp2[0];
                        (*this)[neigh[partit[0][j]][k]] += tmp2[0];
                        coor[0] = neigh[partit[0][j]][k] > partit[0][j] ? neigh[partit[0][j]][k] | (partit[0][j] << 16) : partit[0][j] | (neigh[partit[0][j]][k] << 16);
                        (*this).data[coor[0]] = -ntrgt.data[l] / tmp2[1];
                    }
                }
            }
        }
    }


    /*printf("Error\n");
    ntrgt -= this->mkInverse();
    for(i=0,k=0;i<target.getSize();i++,k++){
        for(j=0;j<i;j++,k++){
            if (data.find(i |(j << 16)) == 0xFFFFFFFF) ntrgt.data[k] = 0.0f;
        }
    }*/
   // ntrgt.show();

    // scale precision matrix to match prior to normalization
    cycinput.d.d.setSize(target.getSize());
    for(i=0;i<target.getSize();i++) cycinput.d.d[i] = ExOp::mkPowInvInt(target.data[(i * (i+3))>> 1] , -2); // 1/sqrt
    for(i=0;i<data.getSize();i++) data.deref(i) *= cycinput.d.d[(data.deref_key(i)& 0xFFFF)] * cycinput.d.d[(data.deref_key(i) >> 16)];

   // printf("Final matrix\n");
   // ((Trianglix<C>)*this).show();


    for(i=0;i<details.getSize();i++){
        details[i].d[1] /= locross_pairs.getSize();
        details[i].d[2] /= locross_pairs.getSize();
        details[i].d[3] /= locross_pairs.getSize();
    }
/*
        Tuple<D, TSIZE > tmp,mean;
        unsigned int totsize = data.getSize();

		this->setSize(totsize*2 -1);

		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;


        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;

        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];


		unsigned int i;
		KeyElem<double, unsigned int > toins;

        printf("Initializing Clustering: "); fflush(stdout);
        for(i=0;i<totsize;i++) stats[i] = data[i];


		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){

                toins.k = stats[links[i]].likelihoodratio_dist(stats[links[toins.d]]);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}

            if (false == loc_merge[links[i]].isEmpty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%'); if ((i & 7) ==0)fflush(stdout);
                toins.d= i;
                merge.insert(toins);
			}
		}

        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isEmpty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isEmpty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();

                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');if ((nbroots & 7) ==0)fflush(stdout);
                (*this)[nbroots].second = toins.k ; //report( stats[links[i]], stats[links[toins.d]]);

                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].toMemfree();
                loc_merge[links[i]].toMemfree();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;


                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].likelihoodratio_dist(stats[links[toins.d]]);
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }

                if (!loc_merge[links[nbroots]].isEmpty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }

        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;

        printf("(DONE)\n"); fflush(stdout);

		//	delete[](counts[0].second);
		//	delete[](counts);
		*/
return 0;}

LFHTEMP void SparseTrianglix<C>::searchRegul2017(const Trianglix<C> target,const  Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 9u> > > &details, double lambda, Trianglix<char> *chkmatch){
    this->setSize(target.getSize());
    uint32_t totsize = target.getSize();
    int partID;
    Tuple<uint32_t> partit[4];

    KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 9u> > details_input;

    HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > > heap_cyclic_merges;
    KeyElem<double, KeyElem<Tuple<uint32_t, 2u>, Tuple<double> > > cycinput;
    Tuple<double> theotherone;



    HeapTree< KeyElem<double, uint32_t> > heap_acyclic_merges;
    HeapTree< KeyElem<double, uint32_t> > heuristic_merges;
    Tuple<uint32_t> dirty_buf; dirty_buf.setSize(totsize); dirty_buf.toZero();
    uint32_t dirtyID = 0;
    Tuple< KeyElem<double, uint32_t> > current_rep; current_rep.setSize(totsize);
    Tuple<uint32_t, 2u> coor; // (coor[0] > coor[1])
    Tuple<uint32_t> bestedge_buf; bestedge_buf.setSize(totsize);

    Tuple< myHashmap<uint32_t, C> > cross_P; cross_P.setSize(cross_pairs.getSize()); // lambda finder scope


    uint32_t bestedge_size;
    KeyElem<double, uint32_t > input;


    // step 1: normalize target;
    Trianglix<C> ntrgt; ntrgt.setSize(totsize);
    int i,j,k,l;

    double norm_constant =0.0;
    Tuple<uint32_t, 4u> match; match.toZero();
    if (chkmatch){
        for(k=0, i=0;i<target.getSize();i++,k++){
            for(j=0;j<i;j++) {
                if (chkmatch->data[k++] != 0) match[0]++;
            }
        }
        match[1] = (target.getSize() * (target.getSize()-1)) >> 1;
        match[1] -= match[0];
    }


    for(k=0, i=0;i<target.getSize();i++){
        for(j=0;j<i;j++) {
            ntrgt.data[k] = target.data[k];
            ntrgt.data[k++] *= ExOp::mkPowInvInt(target[j] * target[i], -2);
        }
        norm_constant += log(target[i]);
        ExOp::toOne(ntrgt.data[k++]);
    }
    Tuple<C> invup;
    double tmp;
    // initialize P to one


    Trianglix<C> trgtinv = ntrgt.mkInverse();


    this->toOne();

    Tuple<C> cross_deter; cross_deter.setSize(cross_pairs.getSize()); ExOp::toZero(cross_deter);
    Tuple<C> cross_trace; cross_trace.setSize(cross_pairs.getSize()); ExOp::toZero(cross_trace);

    Tuple<C> cross_deter_ac; cross_deter_ac.setSize(cross_pairs.getSize()); ExOp::toZero(cross_deter_ac);
    Tuple<C> cross_trace_ac; cross_trace_ac.setSize(cross_pairs.getSize()); ExOp::toZero(cross_trace_ac);

    myHashmap< uint32_t, Tuple<C> > cross_acyclic_comp;
    Tuple<C> cross_acyclic;

    cross_deter += norm_constant; // includes normalization for comparison
    C deter;
    C trace;
    C acyclic;
    ExOp::toZero(deter); ExOp::toZero(acyclic); ExOp::toOne(trace); trace *= target.getSize();

    WeightElem<Tuple<double ,3u> ,2u> crossstat;
    Tuple<double ,3u> crossstat_input;


    Tuple<double ,2u> tmpdoubles;
    ExOp::toZero(details_input.k);
    ExOp::toZero(details_input.d);
    details_input.d[2] -= target.getSize() * cross_pairs.getSize();
    details_input.d[1] += norm_constant * cross_pairs.getSize();
    for(j= 0; j < target.getSize(); j++){
        k = j | (j << 16);
        for(i= 0; i < cross_pairs.getSize(); i+=2){
            crossstat_input[0] = log(cross_P[i][k] = ExOp::mkInverse(cross_pairs[i][j]) );
            cross_deter[i] += crossstat_input[0];
            crossstat_input[1] = log(cross_P[i|1][k] = ExOp::mkInverse(cross_pairs[i|1][j]) );
            cross_deter[i|1] += crossstat_input[1];
            details_input.d[1] += crossstat_input[0] + crossstat_input[1];

            cross_trace[i]   += (crossstat_input[0] = cross_P[i][k] * cross_pairs[i|1][j]);
            cross_trace[i|1] += (crossstat_input[1] = cross_P[i|1][k] * cross_pairs[i][j]);
            details_input.d[2] += crossstat_input[0] + crossstat_input[1];
        }
    }

    crossstat.toZero();

    for(i= 0; i < cross_pairs.getSize(); i++){
        crossstat_input[0] = cross_deter[i];
        crossstat_input[1] = cross_trace[i];
        crossstat_input[2] = cross_deter[i] - cross_trace[i];
        crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
    }
    details_input.d[4] = sqrt(crossstat.getVar()[2]);
    details.push_back(details_input);
    details_input.d[3] = details_input.d[2] - details_input.d[1];

 //   crossstat_input = crossstat.getMean();
 //   printf("0\t%e-%e (%e-%e  ", deter, trace, crossstat_input[0], crossstat_input[1]);
 //   crossstat_input = crossstat.getVar();
 //   printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

    for(k = 1; k < totsize; k++){
        for(input.d =k; input.d < (k << 16); input.d+= 0x10000){
            input.k = ntrgt.data[(input.d >> 16) + ((k * (k+1)) >> 1)];
            input.k = input.k * input.k;
            input.k = log(1.0f - input.k); // - 4.0f * input.k / (1.0f - input.k );
            heap_acyclic_merges.insert(input);
        }
    }
    C tmp2[6];

    while((!heap_acyclic_merges.isEmpty())||(!heap_cyclic_merges.isEmpty())){


        if ((!heap_cyclic_merges.isEmpty())&&((heap_acyclic_merges.isEmpty())||(heap_acyclic_merges.top().k > heap_cyclic_merges.top().k))){
            //printf("Current\n");
            //((Trianglix<C>)*this).show();
            //printf("Difference\n");
            //(this->mkInverse() - ntrgt).show();

            details_input.k[0] = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
            details_input.k[1] = heap_cyclic_merges.top().d.k[0] >> 16;
            details_input.d[0] -= heap_cyclic_merges.top().k;


      //      printf("Cyclic %i,%i is %e\n", heap_cyclic_merges.top().d.k[0] & 0xFFFF, heap_cyclic_merges.top().d.k[0] >> 16, -heap_cyclic_merges.top().k);
            partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0], &partit[1], true, &partit[2]);

            for(i=0;i<partit[0].getSize();i++){
                if (this->isCyclicFID(this->getFirstID(partit[0][i]))) continue;
                for(j=0;j<partit[0].getSize();j++){
                    l = (partit[0][j] < partit[0][i]) ? partit[0][i] | (partit[0][j] << 16) : partit[0][j] | (partit[0][i] << 16);
                    if ((k = cross_acyclic_comp.find(l)) != 0xFFFFFFFF){
                        for(l=0;l<cross_pairs.getSize();l++){
                            cross_trace_ac[l] -= cross_acyclic_comp.deref(k)[l];
                            cross_deter_ac[l] -= cross_acyclic_comp.deref(k)[l+cross_pairs.getSize()];
                            cross_deter[l] += cross_acyclic_comp.deref(k)[l+cross_pairs.getSize()];
                            details_input.d[2] -= cross_acyclic_comp.deref(k)[l];
                        }
                        cross_acyclic_comp.deref(k).toMemfree();
                        cross_acyclic_comp.erase_from_iterator(k);
                        k = (partit[0][j] < partit[0][i]) ? partit[0][j] + ((partit[0][i] * (partit[0][i]+1)) >> 1) : partit[0][i] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        acyclic += log(1.0 - ntrgt.data[k] * ntrgt.data[k]);
                        deter -= log(1.0 - ntrgt.data[k] * ntrgt.data[k]);

                    }
                }
            }

            this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],cycinput.d.d, ntrgt, details_input.d[5], false);

            //printf("got %i acyclic to fix too!\n", partit[3].getSize());
            // replace by edgeID and triindex

            partit[3].setSize(partit[1].getSize() + partit[0].getSize());
            for(i=0;i<partit[1].getSize();i++) {
                partit[3][i] = (partit[0][partit[1][i] & 0xFFFF] < partit[0][partit[1][i] >> 16]) ? partit[0][partit[1][i] >> 16] | (partit[0][partit[1][i] & 0xFFFF] << 16) : partit[0][partit[1][i] & 0xFFFF] | (partit[0][partit[1][i] >> 16] << 16);
            }
            for(i=0;i<partit[0].getSize();i++) partit[3][i + partit[1].getSize()] = partit[0][i] | (partit[0][i] << 16);
            partit[1].setSize(partit[3].getSize());
            for(i=0;i<partit[3].getSize();i++) partit[1][i] = (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

//            for(i=0;i< partit[3].getSize();i++) partit[1][i+partit[3].getSize()] =  (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

            // adding edge in crosses!

            details_input.d[6] = 0.0;
            for(l =0; l< cross_pairs.getSize();l+=2){
                this->data.toMemswap(cross_P[l]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],cycinput.d.d, cross_pairs[l], crossstat_input[1], true);
                this->data.toMemswap(cross_P[l]);
                cross_deter[l] += crossstat_input[0];
                details_input.d[1] += crossstat_input[0];
                details_input.d[6] += crossstat_input[1];

                this->data.toMemswap(cross_P[l|1]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],theotherone, cross_pairs[l|1], crossstat_input[1],true);
                this->data.toMemswap(cross_P[l|1]);

                cross_deter[l|1] += crossstat_input[0];
                details_input.d[1] += crossstat_input[0];
                details_input.d[6] += crossstat_input[1];

                crossstat_input.toZero();
                //tmpdoubles.toZero();
                if (partit[3].getSize() != theotherone.getSize()) printf("adding edge problem! theotherone is has wrong length! %i != %i\n", partit[3].getSize(),theotherone.getSize() );
                if (partit[3].getSize() != cycinput.d.d.getSize()) printf("adding edge problem! cycinput.d.d is has wrong length! %i != %i\n", partit[3].getSize(),cycinput.d.d.getSize() );

                for(i=0;i< partit[3].getSize() - partit[0].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * cross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * cross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * cross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * cross_pairs[l|1].data[partit[1][i]];


                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }
                crossstat_input *= 2.0f; // sym matrix right?
                //tmpdoubles *= 2.0f;

                for(;i<partit[3].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * cross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * cross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * cross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * cross_pairs[l|1].data[partit[1][i]];
                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }

                cross_trace[l] += crossstat_input[0];
                cross_trace[l|1] += crossstat_input[1];
                //cross_deter[l] += tmpdoubles[0];
                //cross_deter[l|1] += tmpdoubles[1];

                details_input.d[2] += crossstat_input[0] + crossstat_input[1];
                //details_input.d[1] += tmpdoubles[1] + tmpdoubles[0];
            }

            crossstat_input[0] = 0.0f;
            for(i=0;i<partit[3].getSize() - partit[0].getSize();i++){
                if ((j = data.find(partit[3][i])) != 0xFFFFFFFF){
                    crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data.deref(j)) * ntrgt.data[partit[1][i]];
                    data.deref(j) = heap_cyclic_merges.top().d.d[i];
                }else{
          //          printf("Adding that entry! P[%i,%i] = %e\n",partit[3][i] & 0xFFFF, partit[3][i] >> 16, heap_cyclic_merges.top().d.d[i] );
                    this->addEntry(partit[3][i]);
                    crossstat_input[0] += heap_cyclic_merges.top().d.d[i] * ntrgt.data[partit[1][i]];
                    data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
                }
            }
            crossstat_input[0] *= 2.0f;
            for(;i<partit[3].getSize();i++){
                crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data[partit[3][i]]) * ntrgt.data[partit[1][i]];
                data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
            }

            trace += crossstat_input[0];
            deter -= heap_cyclic_merges.top().k - crossstat_input[0];

            // get a fresh taste within cycle group! at last
       //     printf("Heuristic quandidates Cycle\n"); fflush(stdout);


            Vector<uint32_t> edges;
            i=0;
            k = part_id[getFirstID(partit[0][0])];
            j = part_data[k][0];
            do{bestedge_buf[i++] = j; j = next_id[j];
            }while(j != part_data[k][0]);
            partit[1].toMemfree();
            edges.push_back(part_data[k][3]);
            for(;edges[0] != 0xFFFFFFFF;edges[0] = attrib[ edges[0]]){
                partit[1].push_back(i);
                k = part_id[getFirstID(edges[0] & 0xFFFF)];
                edges.push_back(part_data[k][3]);
                j = part_data[k][0];
                do{bestedge_buf[i++] = j; j = next_id[j];
                }while(j != part_data[k][0]);
                k = part_id[getFirstID(partit[0][0])];
                while(true){
                    if (edges.last() == 0xFFFFFFFF) {
                        edges.pop_back();
                        if ((k = edges.getSize()) == 1) break;
                        k = part_id[getFirstID(edges[k-2] >>16 )];
                        edges.last() = attrib[edges.last()];
                    }else{
                        j = part_id[getFirstID(edges.last() & 0xFFFF)];
                        if (j != k){
                            k = part_data[j][0];
                            do{bestedge_buf[i++] = k; k = next_id[k];
                            }while(k != part_data[j][0]);
                            k = part_id[getFirstID(edges.last() >> 16)];
                            edges.push_back(part_data[j][3]);
                        }else edges.last() = attrib[edges.last()];
                    }
                }
            }
            bestedge_size = i;
            edges.toMemfree();
            if (partit[1].getSize() == 0) partit[1].push_back(bestedge_size);
            for(i=0;i < partit[1][0];i++){
                for(j=i+1;j<bestedge_size;j++){
                    l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                    if (this->hasEntry(l)) continue;
                    heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                }
            }
            for(k=1;k<partit[1].getSize();k++){
                for(;i < partit[1][k];i++){
                    for(j=partit[1][k];j<bestedge_size;j++){
                        l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                        if (this->hasEntry(l)) continue;
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                    }
                }
            }

            if (bestedge_size > heuristic_merges.getSize()) bestedge_size = heuristic_merges.getSize();
            for(i=0;i<bestedge_size;) bestedge_buf[i++] = heuristic_merges.pop().d;
            heuristic_merges.toMemfree();
            dirtyID++;
            for(i=0;i<partit[0].getSize();i++) dirty_buf[partit[0][i]] = dirtyID;

            if (chkmatch){
                i = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
                j = heap_cyclic_merges.top().d.k[0] >> 16;
                i = (i < j) ? i + ((j * (j+1)) >> 1) : j + ((i * (i+1)) >> 1);
                match[ (chkmatch->data[i] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }


            for(heap_cyclic_merges.pop();!heap_cyclic_merges.isEmpty();heap_cyclic_merges.pop()){
                partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0]);
                for(i=0;i<partit[0].getSize();i++){
                    if (dirty_buf[partit[0][i]] != heap_cyclic_merges.top().d.k[1]) break;
                }
                if (i == partit[0].getSize()) break;
            }
            //printf("Current\n");
            //((Trianglix<C>)*this).show();
            //printf("Difference\n");
            //(this->mkInverse() - ntrgt).show();



        }else{
            k = heap_acyclic_merges.top().d;
            //if (wouldEdgeBeInCyclicGroup(k)) LFH_exit("ImpPPPPosible");


            details_input.k[0] = k & 0xFFFF;
            details_input.k[1] = k >> 16;
            details_input.d[0] -= heap_acyclic_merges.top().k;

            coor[0] = k & 0xFFFF;
            coor[1] = (k >> 16);
       //     printf("Acyclic %i,%i is %e\n", coor[0], coor[1], -heap_acyclic_merges.top().k);

            // find best candidates for *new* cyclic edges!!
            partit[0] = this->getConnectedSet(coor[0]);
            partit[1] = this->getConnectedSet(coor[1]);
      //      printf("Heuristic quandidates Acycle\n"); fflush(stdout);
            if ((partit[0].getSize() > 1)&&(partit[1].getSize() > 1)){
                for(i=0;i<partit[0].getSize();i++){
                    for(j=0;j<partit[1].getSize();j++){
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(target(partit[0][i], partit[1][j])), partit[0][i] | (partit[1][j] << 16)));
                    }
                }
                bestedge_size = partit[0].getSize() + partit[1].getSize() - 1;
                for(i=0;i<bestedge_size;){
                    k = (heuristic_merges.top().d >> 16);
                    j = (heuristic_merges.top().d & 0xFFFF);
                    if (k > j){
                        if ((coor[0] != k)||(coor[1] != j)) bestedge_buf[i++] = (heuristic_merges.top().d >> 16) | (heuristic_merges.top().d << 16);
                    }else{
                        if ((coor[0] != j)||(coor[1] != k)) bestedge_buf[i++] = heuristic_merges.top().d;
                    }
                    heuristic_merges.pop();
                }
                heuristic_merges.toMemfree();
            }else if (partit[0].getSize() == 1){
                bestedge_size = partit[1].getSize() -1;
                for(i=0,j=0;i<partit[1].getSize();i++){
                    if (partit[1][i] == coor[1]) j++;
                    else bestedge_buf[i-j] = (partit[1][i] < partit[0][0]) ? partit[0][0] | (partit[1][i] << 16) : partit[1][i] | (partit[0][0] << 16);
                }
            }else{
                bestedge_size = partit[0].getSize() -1;
                for(i=0,j=0;i<partit[0].getSize();i++){
                    if (partit[0][i] == coor[0]) j++;
                    else bestedge_buf[i-j] = (partit[0][i] < partit[1][0]) ? partit[1][0] | (partit[0][i] << 16) : partit[0][i] | (partit[1][0] << 16);
                }
            }
            k = coor[1] + ((coor[0] * (coor[0] + 1)>>1));
            // does not need to update !!invP!! lols

            //tmp2[3] = ntrgt.data[k];
            //tmp2[2] = tmp2[3] * tmp2[3];
            //tmp2[1] = -tmp2[2] +1.0f;
            //tmp2[0] = tmp2[2] / tmp2[1];

            /*
            partit = this->getConnectedSet(coor[0]);
            invup = this->solveInPartition(partit, coor[0]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];

            partit = this->getConnectedSet(coor[1]);
            invup = this->solveInPartition(partit,coor[1]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];*/

            //(*this)[coor[0]] += tmp2[0];
            //(*/this)[coor[1]] += tmp2[0];
            //(*this)(coor) = -ntrgt.data[k] / tmp2[1]; // ADDING THAT EDGE!
            (*this)(coor) = 0.0f; // ADDING THAT EDGE!
            //trace += heap_acyclic_merges.top().k - tmp2[0];
            //deter -= heap_acyclic_merges.top().k;
            acyclic -= heap_acyclic_merges.top().k;

            cross_acyclic.setSize(cross_pairs.getSize()*2);
            for(i =0; i< cross_pairs.getSize();i+= 2){
                tmp2[0] = cross_pairs[i][coor[0]] * cross_pairs[i][coor[1]];
                tmp2[1] = cross_pairs[i].data[k] * cross_pairs[i].data[k];
                cross_acyclic[cross_pairs.getSize() + i] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                details_input.d[1] += cross_acyclic[cross_pairs.getSize() + i];
                cross_deter_ac[i] += cross_acyclic[cross_pairs.getSize() + i];
                tmp2[5] = tmp2[0] - tmp2[1];
                tmp2[4] = tmp2[1] / tmp2[5];
                tmp2[0] = cross_pairs[i|1][coor[0]] * cross_pairs[i|1][coor[1]];
                tmp2[1] = cross_pairs[i|1].data[k] * cross_pairs[i|1].data[k];
                cross_acyclic[cross_pairs.getSize() + (i|1)] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                details_input.d[1] += cross_acyclic[cross_pairs.getSize() + (i|1)];
                cross_deter_ac[i|1] += cross_acyclic[cross_pairs.getSize() + (i|1)] ;
                tmp2[3] = tmp2[0] - tmp2[1];
                tmp2[2] = tmp2[1] / tmp2[3];

                tmp2[0] = tmp2[4] / cross_pairs[i][coor[0]];
                cross_acyclic[i] = tmp2[0] * cross_pairs[i|1][coor[0]];
                //cross_P[i][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[4] / cross_pairs[i][coor[1]];
                cross_acyclic[i] += tmp2[0] * cross_pairs[i|1][coor[1]];
                //cross_P[i][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -cross_pairs[i].data[k] / tmp2[5];
                cross_acyclic[i] += tmp2[0] * cross_pairs[i|1].data[k] *2.0;
                //cross_P[i][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i] += cross_acyclic[i];



                tmp2[0] = tmp2[2] / cross_pairs[i|1][coor[0]];
                cross_acyclic[i|1] = tmp2[0] * cross_pairs[i][coor[0]];
                //cross_P[i|1][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[2] / cross_pairs[i|1][coor[1]];
                cross_acyclic[i|1] += tmp2[0] * cross_pairs[i][coor[1]];
                //cross_P[i|1][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -cross_pairs[i|1].data[k] / tmp2[3];
                cross_acyclic[i|1] += tmp2[0] * cross_pairs[i].data[k] *2.0;
                //cross_P[i|1][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i|1] += cross_acyclic[i|1];
                details_input.d[2] += cross_acyclic[i] + cross_acyclic[i|1];
            }
            cross_acyclic_comp[heap_acyclic_merges.top().d].toMemmove(cross_acyclic);

            if (chkmatch){
                match[ (chkmatch->data[coor[1] + ((coor[0] * (coor[0] + 1)>>1))] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }else {details_input.d[7] = details_input.d[8] = 0.0;}

            for(heap_acyclic_merges.pop();!heap_acyclic_merges.isEmpty();heap_acyclic_merges.pop()){
                if (!wouldEdgeBeInCyclicGroup(heap_acyclic_merges.top().d)) break;
            }


        }
        details_input.d[3] += details_input.d[1] -  details_input.d[2];

        crossstat.toZero();
        for(i= 0; i < cross_pairs.getSize(); i++){
            crossstat_input[0] = cross_deter[i] + cross_deter_ac[i];
            crossstat_input[1] = cross_trace[i] + cross_trace_ac[i];
            crossstat_input[2] = crossstat_input[0] - crossstat_input[1];
            crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
        }

        details_input.d[4] = sqrt(crossstat.getVar()[2]);
        details.push_back(details_input);
        details_input.d[3] = details_input.d[2] - details_input.d[1];
    //    printf("Start fill cyclic quandidates\n"); fflush(stdout);
        cycinput.d.k[1] = dirtyID;
        for(i=0;i<bestedge_size;i++){
            cycinput.d.k[0] = bestedge_buf[i];
            cycinput.k = -this->solveWithFictiveEdge2017(bestedge_buf[i],cycinput.d.d, ntrgt,crossstat_input[0] ,false);
            heap_cyclic_merges <<= cycinput;
        }
    //    printf("End fill cyclic quandidates %i \n", heap_cyclic_merges.getSize()); fflush(stdout);
/*  checking zone!



        // check trace is behaving!

        crossstat_input[0] =0.0;
        for(i=0; i< data.getSize();i++){
            j = data.deref_key(i); k =(j & 0xFFFF);
            k = (j >> 16) + ((k * (k+1)) >> 1);
            crossstat_input[0] += data.deref(i) * ntrgt.data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        this->wrDeterminant(crossstat_input[1]);
        printf("trace equal: %e ?= %e\n", trace, crossstat_input[0]);
        printf("deter equal: %e ?= %e\n", deter, log(crossstat_input[1]));

        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[0].getSize();i++){
                j = cross_P[0].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[0].deref(i) * cross_pairs[1].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[0], crossstat_input[0]);
        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[1].getSize();i++){
                j = cross_P[1].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[1].deref(i) * cross_pairs[0].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[1], crossstat_input[0]);

        printf("Acyclic traces: %e & %e\n", cross_trace_ac[0], cross_trace_ac[1]);
        printf("Acyclic determinant: %e & %e\n", cross_deter_ac[0], cross_deter_ac[1]);


        this->data.toMemswap(cross_P[0]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[0]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[0]);
        this->data.toMemswap(cross_P[1]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[1]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[1]);


        crossstat_input = crossstat.getMean();
        printf("%i\t%e-%e (%e-%e  ", data.getSize() - totsize, deter, trace, crossstat_input[0], crossstat_input[1]);
        crossstat_input = crossstat.getVar();
        printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

        printf("%i: %e\n", data.getSize() - totsize, deter + acyclic - trace);
        crossstat_input.toZero();
        for(i=0; i< cross_P.getSize();i++){
            crossstat_input[0] += cross_deter_ac[i] - cross_trace_ac[i];
            crossstat_input[1] += cross_deter[i] - cross_trace[i];
        }
        printf("vs: %e\n", (crossstat_input[0]  + crossstat_input[1]) / cross_P.getSize());
        */
        printf("%i edges considered!\n",data.getSize() - totsize);

        if ((data.getSize() - totsize) > lambda) break;
    }

    // formally add acyclic edges
    for(i=0;i<part_data.getSize();i++){
        if (part_data[i][1] > part_data[i][2]){
            partit[0] = this->getPartition(i);
            for(j=0;j<partit[0].getSize();j++){
                for(k=0;k<neigh[partit[0][j]].getSize();k++){
                    if ((getFirstID(neigh[partit[0][j]][k]) != partit[0][0])||(partit[0][j] < neigh[partit[0][j]][k])){
                        l = neigh[partit[0][j]][k] > partit[0][j] ? partit[0][j] + ((neigh[partit[0][j]][k] * (neigh[partit[0][j]][k]+1)) >> 1) :  neigh[partit[0][j]][k] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        tmp2[3] = ntrgt.data[l];
                        tmp2[2] = tmp2[3] * tmp2[3];
                        tmp2[1] = -tmp2[2] +1.0f;
                        tmp2[0] = tmp2[2] / tmp2[1];
                        (*this)[partit[0][j]] += tmp2[0];
                        (*this)[neigh[partit[0][j]][k]] += tmp2[0];
                        coor[0] = neigh[partit[0][j]][k] > partit[0][j] ? neigh[partit[0][j]][k] | (partit[0][j] << 16) : partit[0][j] | (neigh[partit[0][j]][k] << 16);
                        (*this).data[coor[0]] = -ntrgt.data[l] / tmp2[1];
                    }
                }
            }
        }
    }

    /*printf("Error\n");
    ntrgt -= this->mkInverse();
    for(i=0,k=0;i<target.getSize();i++,k++){
        for(j=0;j<i;j++,k++){
            if (data.find(i |(j << 16)) == 0xFFFFFFFF) ntrgt.data[k] = 0.0f;
        }
    }*/
   // ntrgt.show();

    // scale precision matrix to match prior to normalization
    cycinput.d.d.setSize(target.getSize());
    for(i=0;i<target.getSize();i++) cycinput.d.d[i] = ExOp::mkPowInvInt(target.data[(i * (i+3))>> 1] , -2); // 1/sqrt
    for(i=0;i<data.getSize();i++) data.deref(i) *= cycinput.d.d[(data.deref_key(i)& 0xFFFF)] * cycinput.d.d[(data.deref_key(i) >> 16)];

   // printf("Final matrix\n");
   // ((Trianglix<C>)*this).show();


    for(i=0;i<details.getSize();i++){
        details[i].d[1] /= cross_pairs.getSize();
        details[i].d[2] /= cross_pairs.getSize();
        details[i].d[3] /= cross_pairs.getSize();
    }
/*
        Tuple<D, TSIZE > tmp,mean;
        unsigned int totsize = data.getSize();

		this->setSize(totsize*2 -1);

		HeapTree< KeyElem<C, Tuple<int,2> > > pairs;


        GaussElem< Tuple<D, TSIZE > >* stats = new GaussElem< Tuple<D, TSIZE > >[totsize];
        HeapTree< KeyElem<double, unsigned int > > merge;

        unsigned int* links = new unsigned int[totsize*2 -1];
        unsigned int* tom = new unsigned int[totsize];


		unsigned int i;
		KeyElem<double, unsigned int > toins;

        printf("Initializing Clustering: "); fflush(stdout);
        for(i=0;i<totsize;i++) stats[i] = data[i];


		for(i=totsize-1;i<totsize*2-1;i++){
			links[i] = i+ 1- totsize;
            tom[links[i]] = i;
            for(toins.d=totsize-1;toins.d<i;toins.d++){

                toins.k = stats[links[i]].likelihoodratio_dist(stats[links[toins.d]]);
				if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                loc_merge[links[i]].insert(toins);
			}

            if (false == loc_merge[links[i]].isEmpty()){
                toins = loc_merge[links[i]].top();
                printf("partial best : (%i,%i) with %e (%.1f%c done!)\n",i, toins.d, toins.k,  ((50.0f*(i-totsize+1))*(i-totsize+2)) / (double)(totsize * (totsize+1)),'%'); if ((i & 7) ==0)fflush(stdout);
                toins.d= i;
                merge.insert(toins);
			}
		}

        printf("(DONE)\nClustering: "); fflush(stdout);
        nbroots =totsize-2;
        while(!merge.isEmpty()){
            toins = merge.pop();
            i = toins.d;
            if (links[i] == 0xFFFFFFFF) continue;
			//printf("popinmerg %i, %i \n",i,links[i]);fflush(stdout);
            toins = loc_merge[links[i]].pop();
			//          printf("popinmerg\n");fflush(stdout);
            if (links[toins.d] == 0xFFFFFFFF){
                while(!loc_merge[links[i]].isEmpty()){
                    toins = loc_merge[links[i]].top();
					//       printf("popsinmerg\n");fflush(stdout);
                    if (links[toins.d] == 0xFFFFFFFF) loc_merge[links[i]].pop();

                    else {toins.d= i;merge.insert(toins); break;}
					//printf("popsinmerg\n");fflush(stdout);
                }
            }else{
				//      printf("damerge!\n");
				//	printf("%f\t%f\n", stats[links[i]].w,stats[links[toins.d]].w);
                if (stats[links[i]].w > stats[links[toins.d]].w){
                	makeLeftof(i, nbroots);
                    makeRightof(toins.d, nbroots);
                }else{
                	makeLeftof(toins.d, nbroots);
                    makeRightof(i, nbroots);
                }
                printf("%i <- %i,%i with %e (%.1f%c done!)\n",nbroots,i, toins.d, toins.k,  100.0f - ((50.0f*(nbroots+1))*(nbroots)) / (double)(totsize * (totsize+1)),'%');if ((nbroots & 7) ==0)fflush(stdout);
                (*this)[nbroots].second = toins.k ; //report( stats[links[i]], stats[links[toins.d]]);

                stats[links[i]] += stats[links[toins.d]];
                links[nbroots] = links[i];
                loc_merge[links[toins.d]].toMemfree();
                loc_merge[links[i]].toMemfree();
                tom[links[toins.d]] = 0xFFFFFFFF;
                tom[links[i]] = nbroots;
                links[i] = 0xFFFFFFFF;
                links[toins.d] = 0xFFFFFFFF;


                for(i=0;i<totsize;i++){
                    if (tom[i] == 0xFFFFFFFF) continue;
                    if (tom[i] == nbroots) continue;
                    toins.d = tom[i]; // printf("%i is avail\n", tom[i]);
                    toins.k = stats[links[nbroots]].likelihoodratio_dist(stats[links[toins.d]]);
					if (!(ExOp::isValid(toins.k))) ExOp::toMax(toins.k);
                    loc_merge[links[nbroots]].insert(toins);
                }

                if (!loc_merge[links[nbroots]].isEmpty()){
                    toins = loc_merge[links[nbroots]].top();
                    toins.d= nbroots;
                    merge.insert(toins);
                }
                nbroots--;
            }
        }

        delete[](loc_merge);
        delete[](stats);
        delete[](tom);
        delete[](links);
        nbroots =1;

        printf("(DONE)\n"); fflush(stdout);

		//	delete[](counts[0].second);
		//	delete[](counts);
		*/



	}





LFHTEMP class SparseTrianglix<C>::SearchRegulTask{
    public:
    uint32_t nbthreads;
    const Trianglix<C> &target;
    uint32_t totsize;

    Tuple<Tuple<double, 4u> > para_crossstats;
    const  Tuple< Trianglix<C> > &cross_train;
    const  Tuple< Trianglix<C> > &cross_test;
    const  Tuple< uint32_t > &cross_test_size;
    uint32_t total_nb_observations;

    Tuple< myHashmap<uint32_t, C> > cross_P;
    Tuple<C> cross_trace, cross_deter, cross_error, cross_trainace;
    SparseTrianglix<C> current;

    uint32_t dirtyID;

    Tuple<uint32_t> partit[4];
    uint32_t init_t;
    const Trianglix<C> *n;
    Trianglix<C> ntrgt;

    Tuple<uint32_t> bestedge_buf;
    Tuple<double> delta_buf;
    uint32_t bestedge_size;
    //HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > > heap_cyclic_merges;
    //AsyncInserter<HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > > > async;
    HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > >, 3 > heap_cyclic_merges;
    uint32_t bestedge;
    uint32_t curtask;
    double test_size_factor;

    SearchRegulTask(uint32_t _nbthreads, const Trianglix<C> &_target, const  Tuple< Trianglix<C> > &_cross_train,const  Tuple< Trianglix<C> > &_cross_test, const  Tuple< uint32_t > &_cross_test_size ): nbthreads(_nbthreads), target(_target),cross_train(_cross_train), cross_test(_cross_test), cross_test_size(_cross_test_size),totsize(_target.getSize()){ //, async(heap_cyclic_merges){
        init_t = clock();
        bestedge_buf.setSize(totsize);
        para_crossstats.setSize(_nbthreads);
        cross_deter.setSize(_cross_train.getSize()).toZero();
        cross_trace.setSize(_cross_train.getSize()).toZero();
        cross_trainace.setSize(_cross_train.getSize()).toZero();
        cross_error.setSize(_cross_train.getSize()).toZero();
        total_nb_observations =0;
        for(int i =0;i < _cross_test_size.getSize();i++) total_nb_observations += _cross_test_size[i];
        test_size_factor = 1.0 / total_nb_observations;
    }
    uint32_t operator()(uint32_t threadID){
        uint32_t maxranage,i;
        switch(curtask){ case 0:{
            KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > cycinput;
            double err;
            cycinput.d.k[1] = dirtyID;
            if (bestedge_size < nbthreads){
                i = threadID;
                maxranage = i + 1;
            }else{
                i=(threadID * bestedge_size) / nbthreads;
                maxranage = ((threadID + 1) * bestedge_size) / nbthreads;
            }
            foreach.printf("Thread %i running %i-%i\n", threadID, i, maxranage-1);
            fflush(stdout);
            for(;i<maxranage;i++){
                if (bestedge_size > 1) foreach.updateProgress(threadID);
                foreach.printf("%i doing E[%i] = %X\n", threadID, i, bestedge_buf[i]);
                cycinput.d.k[0] = bestedge_buf[i];
                cycinput.k = -current.solveWithFictiveEdge(bestedge_buf[i], current.data,  ntrgt,err ,cycinput.d.d,false,threadID==0);
                foreach.printf("%i insert output!\n", threadID);
                heap_cyclic_merges.insert(cycinput);
            }
        break;}case 1:{
            Tuple<double ,2u> crossstat_input;
            Tuple<double> wouldbe_P;
            uint32_t i;
            maxranage = ((threadID + 1) * cross_train.getSize()) / nbthreads;
            for(uint32_t l = (threadID* cross_train.getSize()) / nbthreads ; l < maxranage; l++){
                crossstat_input[0] = current.solveWithFictiveEdge(bestedge,cross_P[l], cross_train[l], crossstat_input[1],wouldbe_P, true,false);
                cross_deter[l] += crossstat_input[0];
                cross_error[l] = crossstat_input[1];
                para_crossstats[threadID][1] += crossstat_input[0] * cross_test_size[l];
                para_crossstats[threadID][0] += crossstat_input[1];
                foreach.printf("%i: got %e\t%e\n", crossstat_input[0] * cross_test_size[l], crossstat_input[1]);

                crossstat_input.toZero();
                for(i=0;i< partit[3].getSize() - partit[0].getSize();i++){
//                    if (l == 0) printf("read from %X is %e\n", partit[1][i], cross_test[l].data[partit[1][i]]);
                    crossstat_input[0] += (wouldbe_P[i] - cross_P[l][partit[3][i]]) * cross_test[l].data[partit[1][i]] ;
                    crossstat_input[1] += (wouldbe_P[i] - cross_P[l][partit[3][i]]) * cross_train[l].data[partit[1][i]];
                    cross_P[l][partit[3][i]] = wouldbe_P[i];
                }
                crossstat_input *= 2.0f;
                for(;i<partit[3].getSize();i++){
                    crossstat_input[0] += (wouldbe_P[i] - cross_P[l][partit[3][i]]) * cross_test[l].data[partit[1][i]];
                    crossstat_input[1] += (wouldbe_P[i] - cross_P[l][partit[3][i]]) * cross_train[l].data[partit[1][i]];
                    cross_P[l][partit[3][i]] = wouldbe_P[i];
                }
                cross_trace[l] += crossstat_input[0];
                cross_trainace[l] += crossstat_input[1];
                foreach.printf("%i: got %e\t%e\n", crossstat_input[0], crossstat_input[1]);
                para_crossstats[threadID][2] += crossstat_input[0];
                para_crossstats[threadID][4] += crossstat_input[1];
            }
        break;} case 2:
            if (bestedge_size < nbthreads){
                i = threadID;
                maxranage = i + 1;
            }else{
                i=(threadID * bestedge_size) / nbthreads;
                maxranage = ((threadID + 1) * bestedge_size) / nbthreads;
            }
            double err;
            for(;i<maxranage;i++){ foreach.updateProgress(threadID);
                delta_buf[i] = -current.solveWithMissingEdge(bestedge_buf[i], current.data, ntrgt, err);
            }
        }
        foreach.printf("Thread DONE %i!\n", threadID);
    return 0;}
};
LFHTEMP ERRCODE SparseTrianglix<C>::searchRegul2019(const Trianglix<C> &target, const  Tuple< Trianglix<C> > &cross_train,const  Tuple< Trianglix<C> > &cross_test, const Tuple<uint32_t> &cross_test_size, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 13u> > > &details, uint32_t &io_nbedges, Trianglix<char> *chkmatch, uint32_t maxcyclic, uint32_t lowersearchbound, bool force_selected_number_edge){
    if (foreach.getThreadArraySize() == 0) {printf("forgot to start threads!\n"); return 1;}
    else printf("running with %i threads", foreach.getThreadArraySize());
    SparseTrianglix<C>::SearchRegulTask lotask(foreach.getThreadArraySize(), target, cross_train, cross_test,cross_test_size);
    int partID;
    this->setSize(target.getSize());
    lotask.current.setSize(target.getSize());
    if (lowersearchbound == 0) lowersearchbound = target.getSize();
    uint32_t init_t = clock();
    KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 13u> > details_input;
    KeyElem<double, KeyElem<Tuple<uint32_t, 2u>, Tuple<double> > > cycinput;
    Tuple<double> theotherone;
    uint32_t best =0u;
    HeapTree< KeyElem<double, uint32_t> > heap_acyclic_merges;
    HeapTree< KeyElem<double, uint32_t> > heuristic_merges;
    Tuple<uint32_t> dirty_buf; dirty_buf.setSize(target.getSize()); dirty_buf.toZero();
    lotask.dirtyID = 0u;
    Tuple< KeyElem<double, uint32_t> > current_rep; current_rep.setSize(target.getSize());
    Tuple<uint32_t, 2u> coor; // (coor[0] > coor[1])

    Vector<Tuple<double> > Xstats;

    if (cross_train.getSize() != cross_test.getSize()) {printf("Number of test and training sets mismatches!\n"); return 1;}
    lotask.cross_P.setSize(cross_train.getSize()); // lambda finder scope
    int i,j,k,l;
    KeyElem<double, uint32_t > input;

    // step 1: normalize target;

    double tmp;
    lotask.ntrgt.setSize(target.getSize());
    double norm_constant =0.0;

    for(k=0, i=0;i<target.getSize();i++){
        for(j=0;j<i;j++) {
            lotask.ntrgt.data[k] = target.data[k];
            tmp = ExOp::mkPowInvInt(target[j] * target[i], -2);
        //    for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] *= tmp;
            lotask.ntrgt.data[k++] *= tmp;
        }
        norm_constant += log(target[i]);
        //for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] /= target[i];
        ExOp::toOne(lotask.ntrgt.data[k++]);
    }

    Tuple<uint32_t, 4u> match; match.toZero();
    if (chkmatch){
        for(k=0, i=0;i<target.getSize();i++,k++){
            for(j=0;j<i;j++) {
                if (chkmatch->data[k++] != 0) match[0]++;
            }
        }
        match[1] = (target.getSize() * (target.getSize()-1)) >> 1;
        match[1] -= match[0];
    }
    Trianglix<C> trgtinv = target.mkInverse();

    // check that matrices are not singular!
    Tuple<C, 0u> eigen = lotask.ntrgt.getEigenValues();
    for(i=0;i<eigen.getSize();i++){
        if ((!ExOp::isValid(eigen[i]))||(eigen[i] <=0)) break;
    }

    if (i != eigen.getSize()) printf("Warning, got non-positive eigen values in input covariance! \n");
    for(j=0;j<cross_train.getSize();j++) {
        eigen = cross_train[j].getEigenValues();
        for(i=0;i<eigen.getSize();i++){
            if ((!ExOp::isValid(eigen[i]))||(eigen[i] <=0)) break;
        }
        if (i != eigen.getSize()){
            printf("Warning, got non-positive values in training set covariance %i\n", j);
        }
    }

    lotask.current.toOne();
    Tuple<C> cross_deter_ac; cross_deter_ac.setSize(cross_train.getSize()); ExOp::toZero(cross_deter_ac);
    Tuple<C> cross_trace_ac; cross_trace_ac.setSize(cross_train.getSize()); ExOp::toZero(cross_trace_ac);

    myHashmap< uint32_t, Tuple<C> > cross_acyclic_comp;
    Tuple<C> cross_acyclic;

    C deter, trace, acyclic;
    ExOp::toZero(deter); ExOp::toZero(acyclic); ExOp::toOne(trace); trace *= target.getSize();

    WeightElem<Tuple<double ,3u> ,2u> crossstat;
    Tuple<double ,3u> crossstat_input;
    Tuple<double ,2u> tmpdoubles;
    ExOp::toZero(details_input.k);
    ExOp::toZero(details_input.d); /*LL,  , mean logDeterminant, mean lo */
    for(j= 0; j < target.getSize(); j++){
        details_input.d[0] -= log(target[j]);
        k = j | (j << 16);
        for(i= 0; i < cross_train.getSize(); i++){
            crossstat_input[0] = log(lotask.cross_P[i][k] = ExOp::mkInverse(cross_train[i][j]) );
            lotask.cross_deter[i] += crossstat_input[0];
            details_input.d[3] += crossstat_input[0] * cross_test_size[i];
            crossstat_input[0] = lotask.cross_P[i][k] * cross_test[i][j];
            lotask.cross_trace[i] += crossstat_input[0];
            //details_input.d[4] += crossstat_input[0];
        }
    }

    crossstat.toZero();
    Xstats.push_back();
    Xstats.last().setSize(cross_train.getSize()*3);
    details_input.d[3] = 0;details_input.d[4] = 0;
    for(i= 0; i < cross_train.getSize(); i++){
        crossstat_input[0] = lotask.cross_deter[i];
        details_input.d[3] += lotask.cross_deter[i] * cross_test_size[i];
        crossstat_input[1] = lotask.cross_trace[i];
        details_input.d[4] += crossstat_input[1];
        crossstat_input[2] = lotask.cross_deter[i] - lotask.cross_trace[i];
        crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
        Xstats.last()[i] = lotask.cross_deter[i];
        Xstats.last()[i + cross_train.getSize()] = lotask.cross_trace[i];
        Xstats.last()[i + cross_train.getSize()*2] = 0.0;
    }
    details_input.d[1] = (crossstat.getMean()[2] - log(M_2PI))*0.5; // mean of LLS
    details_input.d[6] = sqrt(crossstat.getVar()[2])*0.5; // STD of LLs
    details_input.d[3] /= lotask.total_nb_observations; details.push_back(details_input);
    details_input.d[5] = details_input.d[2];
 //   crossstat_input = crossstat.getMean();
 //   printf("0\t%e-%e (%e-%e  ", deter, trace, crossstat_input[0], crossstat_input[1]);
 //   crossstat_input = crossstat.getVar();
 //   printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));
    for(k = 1; k < target.getSize(); k++){
        for(input.d =k; input.d < (k << 16); input.d+= 0x10000){
            input.k = lotask.ntrgt.data[(input.d >> 16) + ((k * (k+1)) >> 1)];
            input.k = log(1.0f - input.k * input.k);
            heap_acyclic_merges.insert(input);
        }
    }
    C tmp2[6];

    KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > cyclictop;
    bool has_cyclic;
    while((has_cyclic = lotask.heap_cyclic_merges.top(cyclictop))||(!heap_acyclic_merges.isEmpty())){
        #ifdef Rcpp_hpp
            //if (Rcpp::R_TopLevelExec == NULL) printf("function is not defined...\n");
            Rcpp::checkUserInterrupt();
        #endif // Rcpp_hpp

        if ((has_cyclic)&&((heap_acyclic_merges.isEmpty())||(heap_acyclic_merges.top().k > cyclictop.k))){

            details_input.k[0] = cyclictop.d.k[0] & 0xFFFF;
            details_input.k[1] = cyclictop.d.k[0] >> 16;
            details_input.d[0] -= cyclictop.k;

            lotask.partit[0] = lotask.current.getWouldBeCyclicPartition(cyclictop.d.k[0]);

      //      printf("Cyclic %i,%i is %e\n", heap_cyclic_merges.top().d.k[0] & 0xFFFF, heap_cyclic_merges.top().d.k[0] >> 16, -heap_cyclic_merges.top().k);
            lotask.partit[0] = lotask.current.getWouldBeCyclicPartition(cyclictop.d.k[0], &lotask.partit[1], true, &lotask.partit[2]);

            details_input.d[11] = lotask.partit[0].getSize();
            for(i=0;i<lotask.partit[0].getSize();i++){
                if (lotask.current.isCyclicFID(lotask.current.getFirstID(lotask.partit[0][i]))) continue;
                for(j=0;j<lotask.partit[0].getSize();j++){
                    l = (lotask.partit[0][j] < lotask.partit[0][i]) ? lotask.partit[0][i] | (lotask.partit[0][j] << 16) : lotask.partit[0][j] | (lotask.partit[0][i] << 16);
                    if ((k = cross_acyclic_comp.find(l)) != 0xFFFFFFFF){
                        for(l=0;l<cross_train.getSize();l++){
                            cross_trace_ac[l] -= cross_acyclic_comp.deref(k)[l];
                            cross_deter_ac[l] -= cross_acyclic_comp.deref(k)[l+cross_train.getSize()];
                            lotask.cross_deter[l] += cross_acyclic_comp.deref(k)[l+cross_train.getSize()];
                            //details_input.d[4] -= cross_acyclic_comp.deref(k)[l];
                        }
                        cross_acyclic_comp.deref(k).toMemfree();
                        cross_acyclic_comp.erase_from_iterator(k);
                        k = (lotask.partit[0][j] < lotask.partit[0][i]) ? lotask.partit[0][j] + ((lotask.partit[0][i] * (lotask.partit[0][i]+1)) >> 1) : lotask.partit[0][i] + ((lotask.partit[0][j] * (lotask.partit[0][j]+1)) >> 1);
                        acyclic += log(1.0 - lotask.ntrgt.data[k] * lotask.ntrgt.data[k]);
                        deter -= log(1.0 - lotask.ntrgt.data[k] * lotask.ntrgt.data[k]);
                    }
                }
            }

            //printf("%e == %e right?\n", -cyclictop.k,
            lotask.current.solveWithFictiveEdge(cyclictop.d.k[0],lotask.current.data,lotask.ntrgt, details_input.d[7], cycinput.d.d, false,false);

            lotask.partit[3].setSize(lotask.partit[1].getSize() + lotask.partit[0].getSize());
            for(i=0;i<lotask.partit[1].getSize();i++) {
                lotask.partit[3][i] = (lotask.partit[0][lotask.partit[1][i] & 0xFFFF] < lotask.partit[0][lotask.partit[1][i] >> 16]) ? lotask.partit[0][lotask.partit[1][i] >> 16] | (lotask.partit[0][lotask.partit[1][i] & 0xFFFF] << 16) : lotask.partit[0][lotask.partit[1][i] & 0xFFFF] | (lotask.partit[0][lotask.partit[1][i] >> 16] << 16);
            }
            for(i=0;i<lotask.partit[0].getSize();i++) lotask.partit[3][i + lotask.partit[1].getSize()] = lotask.partit[0][i] | (lotask.partit[0][i] << 16);
            lotask.partit[1].setSize(lotask.partit[3].getSize());
            for(i=0;i<lotask.partit[3].getSize();i++) lotask.partit[1][i] = (lotask.partit[3][i] >> 16) + (((lotask.partit[3][i] & 0xFFFF) * ((lotask.partit[3][i] & 0xFFFF)+1)) >> 1);

//            for(i=0;i< lotask.partit[3].getSize();i++) lotask.partit[1][i+lotask.partit[3].getSize()] =  (lotask.partit[3][i] >> 16) + (((lotask.partit[3][i] & 0xFFFF) * ((lotask.partit[3][i] & 0xFFFF)+1)) >> 1);

            // adding edge in crosses!
            details_input.d[8] = 0.0;

            lotask.curtask = 1; lotask.para_crossstats.toZero();
            lotask.bestedge = cyclictop.d.k[0];
            for(i=foreach.getThreadArraySize()-1;i>0;i--) foreach.submit(lotask, i);
            foreach.submit_ThenWait(lotask, i);
            for(i=0;i<foreach.getThreadArraySize();i++) {
                //details_input.d[4] += lotask.para_crossstats[i][2]; // test trace
                //details_input.d[3] += lotask.para_crossstats[i][1]; // delta-log-determinant
                details_input.d[8] += lotask.para_crossstats[i][0]; // error
            }
            crossstat_input[0] = 0.0f;
            for(i=0;i<lotask.partit[3].getSize() - lotask.partit[0].getSize();i++){
                if ((j = lotask.current.data.find(lotask.partit[3][i])) != 0xFFFFFFFF){
                    crossstat_input[0] += (cycinput.d.d[i] - lotask.current.data.deref(j)) * lotask.ntrgt.data[lotask.partit[1][i]];
                    lotask.current.data.deref(j) = cycinput.d.d[i];
                }else{
                    lotask.current.addEntry(lotask.partit[3][i]);
                    crossstat_input[0] += cycinput.d.d[i] * lotask.ntrgt.data[lotask.partit[1][i]];
                    lotask.current.data[lotask.partit[3][i]] = cycinput.d.d[i];
                }
            }

            crossstat_input[0] *= 2.0f;
            for(;i<lotask.partit[3].getSize();i++){
                crossstat_input[0] += (cycinput.d.d[i] - lotask.current.data[lotask.partit[3][i]]) * lotask.ntrgt.data[lotask.partit[1][i]];
                lotask.current.data[lotask.partit[3][i]] = cycinput.d.d[i];
            }

            trace += crossstat_input[0];
            deter -= cyclictop.k - crossstat_input[0];

            // get a fresh taste within cycle group! at last
       //     printf("Heuristic quandidates Cycle\n"); fflush(stdout);
            // finding new candidates within cycle or leading to new cycles

            Vector<uint32_t> edges;
            i=0;
            k = lotask.current.part_id[lotask.current.getFirstID(lotask.partit[0][0])];
            j = lotask.current.part_data[k][0];
            do{lotask.bestedge_buf[i++] = j; j = lotask.current.next_id[j];
            }while(j != lotask.current.part_data[k][0]);
            lotask.partit[1].toMemfree();
            edges.push_back(lotask.current.part_data[k][3]);
            for(;edges[0] != 0xFFFFFFFF;edges[0] = lotask.current.attrib[ edges[0]]){
                lotask.partit[1].push_back(i);
                k = lotask.current.part_id[lotask.current.getFirstID(edges[0] & 0xFFFF)];
                edges.push_back(lotask.current.part_data[k][3]);
                j = lotask.current.part_data[k][0];
                do{lotask.bestedge_buf[i++] = j; j = lotask.current.next_id[j];
                }while(j != lotask.current.part_data[k][0]);
                k = lotask.current.part_id[lotask.current.getFirstID(lotask.partit[0][0])];
                while(true){
                    if (edges.last() == 0xFFFFFFFF) {
                        edges.pop_back();
                        if ((k = edges.getSize()) == 1) break;
                        k = lotask.current.part_id[lotask.current.getFirstID(edges[k-2] >>16 )];
                        edges.last() = lotask.current.attrib[edges.last()];
                    }else{
                        j = lotask.current.part_id[lotask.current.getFirstID(edges.last() & 0xFFFF)];
                        if (j != k){
                            k = lotask.current.part_data[j][0];
                            do{lotask.bestedge_buf[i++] = k; k = lotask.current.next_id[k];
                            }while(k != lotask.current.part_data[j][0]);
                            k = lotask.current.part_id[lotask.current.getFirstID(edges.last() >> 16)];
                            edges.push_back(lotask.current.part_data[j][3]);
                        }else edges.last() = lotask.current.attrib[edges.last()];
                    }
                }
            }
            lotask.bestedge_size = i;
            edges.toMemfree();
            if (lotask.partit[1].getSize() == 0) lotask.partit[1].push_back(lotask.bestedge_size);
            for(i=0;i < lotask.partit[1][0];i++){
                for(j=i+1;j<lotask.bestedge_size;j++){
                    l = (lotask.bestedge_buf[i] < lotask.bestedge_buf[j]) ? lotask.bestedge_buf[j] | (lotask.bestedge_buf[i] << 16) : lotask.bestedge_buf[i] | (lotask.bestedge_buf[j] << 16);
                    if (lotask.current.hasEntry(l)) continue;
                    heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(lotask.bestedge_buf[i], lotask.bestedge_buf[j])), l));
                }
            }
            for(k=1;k<lotask.partit[1].getSize();k++){
                for(;i < lotask.partit[1][k];i++){
                    for(j=lotask.partit[1][k];j<lotask.bestedge_size;j++){
                        l = (lotask.bestedge_buf[i] < lotask.bestedge_buf[j]) ? lotask.bestedge_buf[j] | (lotask.bestedge_buf[i] << 16) : lotask.bestedge_buf[i] | (lotask.bestedge_buf[j] << 16);
                        if (lotask.current.hasEntry(l)) continue;
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(lotask.bestedge_buf[i], lotask.bestedge_buf[j])), l));
                    }
                }
            }

            if (lotask.bestedge_size > heuristic_merges.getSize()) lotask.bestedge_size = heuristic_merges.getSize();
            for(i=0;i<lotask.bestedge_size;) lotask.bestedge_buf[i++] = heuristic_merges.pop().d;
            heuristic_merges.toMemfree();
            lotask.dirtyID++;
            for(i=0;i<lotask.partit[0].getSize();i++) dirty_buf[lotask.partit[0][i]] = lotask.dirtyID;

            if (chkmatch){
                i = cyclictop.d.k[0] & 0xFFFF;
                j = cyclictop.d.k[0] >> 16;
                i = (i < j) ? i + ((j * (j+1)) >> 1) : j + ((i * (i+1)) >> 1);
                match[ (chkmatch->data[i] != 0) ? 2 : 3]++;
                details_input.d[9] = ((double)match[2]) / match[0];
                details_input.d[10] = ((double)match[3]) / match[1];
            }
            // makes sure next cyclic merge is not dirty, and not of a too large size!

            lotask.heap_cyclic_merges.pop(cyclictop);
            has_cyclic = lotask.heap_cyclic_merges.top(cyclictop);

            //printf("Current\n");((Trianglix<C>)*this).show();printf("Difference\n");(lotask.current.mkInverse() - lotask.ntrgt).show();
        }else{
            printf("in acyclic\n"); fflush(stdout);
            k = heap_acyclic_merges.top().d;
            //if (wouldEdgeBeInCyclicGroup(k)) LFH_exit("ImpPPPPosible");
            details_input.d[11] = 0;


            details_input.k[0] = k & 0xFFFF;
            details_input.k[1] = k >> 16;
            details_input.d[0] -= heap_acyclic_merges.top().k;

            coor[0] = k & 0xFFFF;
            coor[1] = (k >> 16);
       //     printf("Acyclic %i,%i is %e\n", coor[0], coor[1], -heap_acyclic_merges.top().k);

            // find best candidates for *new* cyclic edges!!
            lotask.partit[0] = lotask.current.getConnectedSet(coor[0]);
            lotask.partit[1] = lotask.current.getConnectedSet(coor[1]);
      //      printf("Heuristic quandidates Acycle\n"); fflush(stdout);
            if ((lotask.partit[0].getSize() > 1)&&(lotask.partit[1].getSize() > 1)){
                for(i=0;i<lotask.partit[0].getSize();i++){
                    for(j=0;j<lotask.partit[1].getSize();j++){
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(target(lotask.partit[0][i], lotask.partit[1][j])), lotask.partit[0][i] | (lotask.partit[1][j] << 16)));
                    }
                }
                lotask.bestedge_size = lotask.partit[0].getSize() + lotask.partit[1].getSize() - 1;
                for(i=0;i<lotask.bestedge_size;){
                    k = (heuristic_merges.top().d >> 16);
                    j = (heuristic_merges.top().d & 0xFFFF);
                    if (k > j){
                        if ((coor[0] != k)||(coor[1] != j)) lotask.bestedge_buf[i++] = (heuristic_merges.top().d >> 16) | (heuristic_merges.top().d << 16);
                    }else{
                        if ((coor[0] != j)||(coor[1] != k)) lotask.bestedge_buf[i++] = heuristic_merges.top().d;
                    }
                    heuristic_merges.pop();
                }
                heuristic_merges.toMemfree();
            }else if (lotask.partit[0].getSize() == 1){
                lotask.bestedge_size = lotask.partit[1].getSize() -1;
                for(i=0,j=0;i<lotask.partit[1].getSize();i++){
                    if (lotask.partit[1][i] == coor[1]) j++;
                    else lotask.bestedge_buf[i-j] = (lotask.partit[1][i] < lotask.partit[0][0]) ? lotask.partit[0][0] | (lotask.partit[1][i] << 16) : lotask.partit[1][i] | (lotask.partit[0][0] << 16);
                }
            }else{
                lotask.bestedge_size = lotask.partit[0].getSize() -1;
                for(i=0,j=0;i<lotask.partit[0].getSize();i++){
                    if (lotask.partit[0][i] == coor[0]) j++;
                    else lotask.bestedge_buf[i-j] = (lotask.partit[0][i] < lotask.partit[1][0]) ? lotask.partit[1][0] | (lotask.partit[0][i] << 16) : lotask.partit[0][i] | (lotask.partit[1][0] << 16);
                }
            }
            k = coor[1] + ((coor[0] * (coor[0] + 1)>>1));

            lotask.current(coor) = 0.0f;
            acyclic -= heap_acyclic_merges.top().k;

            cross_acyclic.setSize(cross_train.getSize()*2);
            for(i =0; i< cross_train.getSize();i++){
                tmp2[0] = cross_train[i][coor[0]] * cross_train[i][coor[1]];
                tmp2[1] = cross_train[i].data[k] * cross_train[i].data[k];
                cross_acyclic[cross_train.getSize() + i] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                cross_deter_ac[i] += cross_acyclic[cross_train.getSize() + i];
                tmp2[5] = tmp2[0] - tmp2[1];
                tmp2[4] = tmp2[1] / tmp2[5];
                tmp2[0] = tmp2[4] / cross_train[i][coor[0]];
                cross_acyclic[i] = tmp2[0] * cross_test[i][coor[0]];
                tmp2[0] = tmp2[4] / cross_train[i][coor[1]];
                cross_acyclic[i] += tmp2[0] * cross_test[i][coor[1]];
                tmp2[0] = -cross_train[i].data[k] / tmp2[5];
                cross_acyclic[i] += tmp2[0] * cross_test[i].data[k] *2.0;
                cross_trace_ac[i] += cross_acyclic[i];
            }
            cross_acyclic_comp[heap_acyclic_merges.top().d].toMemmove(cross_acyclic);

            if (chkmatch){
                match[ (chkmatch->data[coor[1] + ((coor[0] * (coor[0] + 1)>>1))] != 0) ? 2 : 3]++;
                details_input.d[9] = ((double)match[2]) / match[0];
                details_input.d[10] = ((double)match[3]) / match[1];
            }else {details_input.d[9] = details_input.d[10] = 0.0;}

            for(heap_acyclic_merges.pop();!heap_acyclic_merges.isEmpty();heap_acyclic_merges.pop()){
                if (!lotask.current.wouldEdgeBeInCyclicGroup(heap_acyclic_merges.top().d)) break;
            }
            printf("out acyclic\n"); fflush(stdout);
        }

        // make sure top is not dirty!
        if (has_cyclic) do{
            lotask.partit[0] = lotask.current.getWouldBeCyclicPartition(cyclictop.d.k[0]);
            if ((maxcyclic != 0)&& (lotask.partit[0].getSize() > maxcyclic)) continue;
            for(i=0;i<lotask.partit[0].getSize();i++){
                if (dirty_buf[lotask.partit[0][i]] > cyclictop.d.k[1]) break;
            }
            if (i == lotask.partit[0].getSize()) break;
        }while(true==(has_cyclic = lotask.heap_cyclic_merges.pop(cyclictop)));




        //details_input.d[5] += (details_input.d[3] -  details_input.d[4]) * 0.5;

        crossstat.toZero();
        Xstats.push_back();
        Xstats.last().setSize(cross_train.getSize()*3);
        details_input.d[3] = 0; details_input.d[4] = 0;
        for(i= 0; i < cross_train.getSize(); i++){
            crossstat_input[0] = lotask.cross_deter[i] + cross_deter_ac[i]; // LL
            details_input.d[3] += crossstat_input[0] * cross_test_size[i];
            crossstat_input[1] = lotask.cross_trace[i] + cross_trace_ac[i];
            details_input.d[4] += crossstat_input[1];
            crossstat_input[2] = crossstat_input[0] - crossstat_input[1];
            crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
            crossstat_input.show();
            Xstats.last()[i] = lotask.cross_deter[i] + cross_deter_ac[i];
            Xstats.last()[i + cross_train.getSize()] = lotask.cross_trace[i] + cross_trace_ac[i];
            Xstats.last()[i + cross_train.getSize()*2] = lotask.cross_error[i];
        }

        /*details_input.d[5] = 0;
        for(i= 0; i < cross_train.getSize(); i++) details_input.d[5] += lotask.cross_trace[i] + cross_trace_ac[i];
        printf("Trace sums... %e and %e\n", details_input.d[3], details_input.d[5]);
        details_input.d[5] = 0;
        for(i= 0; i < cross_train.getSize(); i++) details_input.d[5] += lotask.cross_deter[i] + cross_deter_ac[i];
        printf("logdeter sums... %e and %e\n", details_input.d[2], details_input.d[5]);*/

        // details_input.d[0] log det
        // details_input.d[1] LL training
        // details_input.d[2] LL test set

        // details_input.d[3] weithed average log det
        // details_input.d[4] sum trace

        // details_input.d[5] deltaLL for test
        // details_input.d[6] std crxLL
        tmp = details_input.d[3] - details_input.d[4];
        details_input.d[1] = (details_input.d[0] - (log(M_2PI) + 1.0)*  target.getSize()) * 0.5 * lotask.total_nb_observations;
        details_input.d[2] = (tmp  - log(M_2PI) * target.getSize() * lotask.total_nb_observations) * 0.5;

        details_input.d[6] = sqrt(crossstat.getVar()[2]) * 0.5;
        details_input.d[12] = clock() - init_t;
        details_input.d[3] /= lotask.total_nb_observations;
        details_input.d[5] = details_input.d[2] - details_input.d[5];
        if (force_selected_number_edge){
            if (details.getSize() <= io_nbedges){
                printf("%i left to check (TrainLL %e TestLL %e)\n", io_nbedges - details.getSize() , details_input.d[2] , details_input.d[1]);
            }else{
                best = details.getSize();
                *this = lotask.current;
                details.push_back(details_input);
                break;
            }
        }else{
            if ((ExOp::isValid(details_input.d[5]))&&(details_input.d[5] > details[best].d[3] - details[best].d[4])){
                best = details.getSize();
                *this = lotask.current; // stores best
                printf("%i edges currently predicted! (TrainLL %e TestLL %e)\n", best, details_input.d[1] , details_input.d[2]);
            }else{
                uint32_t tmp = io_nbedges + best;
                if (tmp < lowersearchbound) tmp = lowersearchbound;
                printf("%i edges still (%i left to check) (TrainLL %e TestLL %e)\n", best,tmp - details.getSize(), details_input.d[1] , details_input.d[2]);
                if (( ((int32_t)details.getSize()) - best >= io_nbedges)&&(details.getSize() >= lowersearchbound)) {
                    details.push_back(details_input);
                    break;
                }
            }
        }
        details.push_back(details_input); details_input.d[5] = details_input.d[2];
    //    printf("Start fill cyclic quandidates\n"); fflush(stdout);
        lotask.curtask = 0;
        printf("run %i new cyclic edge candidates! (for %i threads)\n", lotask.bestedge_size, foreach.getThreadArraySize());

        /*if (lotask.bestedge_size >0){
            Tuple<uint32_t> hehe,haha;
            lotask.current.show();
            printf("Ex Edge %X\n", lotask.bestedge_buf[0]);fflush(stdout);
            lotask.current.getWouldBeCyclicPartition(lotask.bestedge_buf[0], &hehe,false,&haha);
            hehe.show();
            haha.show();
            fflush(stdout);
        }*/
        if (lotask.bestedge_size > foreach.getThreadArraySize()) {
            foreach.startProgress("Search for next edge", lotask.bestedge_size);
            for(i=foreach.getThreadArraySize()-1;i>0;i--) foreach.submit(lotask, i);
            foreach.submit_ThenWait(lotask, i);
        }else if (lotask.bestedge_size > 1){
            foreach.startProgress("Search for next edge", lotask.bestedge_size);
            for(i=lotask.bestedge_size-1;i>0;i--) foreach.submit(lotask, i);
            foreach.submit_ThenWait(lotask, i);
        }else if (lotask.bestedge_size == 1) lotask(0);
        printf("Init!\n");
        SparseTrianglix<double> chkchk = *this;
        chkchk.solveAllAcyclicEdges(lotask.ntrgt,true);
        double dares; chkchk.wrLogDeterminant(dares);
        printf("%i: %e\t%e are logdets\n", details.getSize(), details_input.d[0], dares);

        if (details.getSize() >= 100) break;
    }
    io_nbedges = best;

    lotask.bestedge_size = io_nbedges;
    lotask.delta_buf.setSize(lotask.bestedge_size);
    lotask.curtask = 2;
    if (lotask.bestedge_size < foreach.getThreadArraySize()) {
        if (lotask.bestedge_size != 0){
            foreach.startProgress("Get edges delta likelihood", lotask.bestedge_size);
            for(i=lotask.bestedge_size-1;i>0;i--) foreach.submit(lotask, i);
            foreach.submit_ThenWait(lotask, i);
        }
    }else{
        foreach.startProgress("Get edges delta likelihood", lotask.bestedge_size);
        for(i=foreach.getThreadArraySize()-1;i>0;i--) foreach.submit(lotask, i);
        foreach.submit_ThenWait(lotask, i);
    }


    // fill delta likelihoods


    this->solveAllAcyclicEdges(lotask.ntrgt,true);

    TMatrix<double> tmpshow;
    tmpshow.setSizes(Xstats.getSize(), cross_train.getSize());
    if (auto ite = tmpshow.getIterator()) do{
        *ite = Xstats[ite()[0]][ite()[1]];
    }while(ite++);
    printf("cross deter\n");
    tmpshow.show();
   if (auto ite = tmpshow.getIterator()) do{
        *ite = Xstats[ite()[0]][ite()[1]+ cross_train.getSize()];
    }while(ite++);
    printf("cross trace\n");
    tmpshow.show();
   if (auto ite = tmpshow.getIterator()) do{
        *ite = Xstats[ite()[0]][ite()[1]+ cross_train.getSize()*2];
    }while(ite++);
    printf("cross error\n");
    tmpshow.show();
    /*   if (auto ite = tmpshow.getIterator()) do{
        *ite = Xstats[ite()[0]][ite()[1]+ cross_train.getSize()*3];
    }while(ite++);
    printf("cross train trace...\n");
    tmpshow.show();*/

    // scale output

    cycinput.d.d.setSize(target.getSize());
    for(i=0;i<target.getSize();i++) cycinput.d.d[i] = ExOp::mkPowInvInt(target.data[(i * (i+3))>> 1] , -2); // 1/sqrt
    for(i=0;i<this->data.getSize();i++) this->data.deref(i) *= cycinput.d.d[(this->data.deref_key(i)& 0xFFFF)] * cycinput.d.d[(this->data.deref_key(i) >> 16)];

return 0;}
/** \brief Search Regulatory Network, where the training and test sets are lists of binary partitions for the same observations (50%/50%)
 *
 * \param
 * \param
 * \return
 *
 */

LFHTEMP ERRCODE SparseTrianglix<C>::searchRegulHalf(const Trianglix<C> target, const Tuple< Trianglix<C> > &cross_pairs, Vector< KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > > &details, double lambda, Trianglix<char> *chkmatch, uint32_t maxcyclic){
    this->setSize(target.getSize());
    uint32_t totsize = target.getSize();
    int partID;
    Tuple<uint32_t> partit[4];
    uint32_t init_t = clock();
    printf("called half original seach!\n"); fflush(stdout);

    KeyElem<Tuple<uint32_t, 2u> , Tuple<double, 11u> > details_input;

    HeapTree< KeyElem<double, KeyElem< Tuple<uint32_t, 2u>, Tuple<double> > > > heap_cyclic_merges;
    KeyElem<double, KeyElem<Tuple<uint32_t, 2u>, Tuple<double> > > cycinput;
    Tuple<double> theotherone;



    HeapTree< KeyElem<double, uint32_t> > heap_acyclic_merges;
    HeapTree< KeyElem<double, uint32_t> > heuristic_merges;
    Tuple<uint32_t> dirty_buf; dirty_buf.setSize(totsize); dirty_buf.toZero();
    uint32_t dirtyID = 0;
    Tuple< KeyElem<double, uint32_t> > current_rep; current_rep.setSize(totsize);
    Tuple<uint32_t, 2u> coor; // (coor[0] > coor[1])
    Tuple<uint32_t> bestedge_buf; bestedge_buf.setSize(totsize);

    Tuple< myHashmap<uint32_t, C> > cross_P; cross_P.setSize(cross_pairs.getSize()); // lambda finder scope


    uint32_t bestedge_size;
    KeyElem<double, uint32_t > input;


    // step 1: normalize target;
    Trianglix<C> ntrgt; ntrgt.setSize(totsize);



    int i,j,k,l;

    double norm_constant =0.0;
    Tuple<uint32_t, 4u> match; match.toZero();
    if (chkmatch){
        for(k=0, i=0;i<target.getSize();i++,k++){
            for(j=0;j<i;j++) {
                if (chkmatch->data[k++] != 0) match[0]++;
            }
        }
        match[1] = (target.getSize() * (target.getSize()-1)) >> 1;
        match[1] -= match[0];
    }

    Tuple< Trianglix<C> > locross_pairs = cross_pairs;
    double tmp;
    for(k=0, i=0;i<target.getSize();i++){
        for(j=0;j<i;j++) {
            ntrgt.data[k] = target.data[k];
            tmp = ExOp::mkPowInvInt(target[j] * target[i], -2);
            for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] *= tmp;
            ntrgt.data[k++] *= tmp;
        }
        norm_constant += log(target[i]);
        for(partID=0;partID< locross_pairs.getSize();partID++) locross_pairs[partID].data[k] /= target[i];
        ExOp::toOne(ntrgt.data[k++]);
    }
    Tuple<C> invup;

    // initialize P to one

    Tuple<C, 0u> eigen = ntrgt.getEigenValues();

    // makes sure matrix is not singular...
   for(k=0;k<100;k++){
       eigen = ntrgt.getEigenValues();
       j = 0;
       for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
       if  (eigen[j] >0.0) break;
       if (eigen[j] > -0.0000001) tmp= 0.999999;
       else tmp = 1.0 / (1.0 - eigen[j]);
       if (k == 0) {printf("warning: non-positive definite covariance, adding prior, start eigenval:\n");eigen.show();}
       for(i=0,l=0;i<ntrgt.getSize();i++,l++){
           for(j=0;j<i;j++) ntrgt.data[l++] *= tmp;
       }
   }
   if (k != 0){
       ntrgt.show();
       printf("final eigenvalues:");
       eigen.show();
       if (k ==100) {printf("failed to make the matrix positive definite ... really?\n"); return 1;}
   }

   for(k=0;k<100;k++){
       for(partID=0;partID< locross_pairs.getSize();partID+=2){
           eigen = locross_pairs[partID].getEigenValues();
           j = 0;
           for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
           if ((partID == 0)||(eigen[j] < tmp)) tmp = eigen[j];
       }
       if (tmp > 0.0) break;
       if (k == 0) {printf("warning: non-positive definite  covariance for cross validation, adding prior\n");eigen.show();}
       if (tmp > -0.0000001) tmp = 1.0 - 0.000001 * k;
       else tmp = 1.0 / (1.0 - tmp * (k+1));

       for(i=0,l=0;i<ntrgt.getSize();i++,l++){
           for(j=0;j<i;j++,l++) {
               for(partID=0;partID< locross_pairs.getSize();partID+=2) locross_pairs[partID].data[l] *= tmp;
           }
       }
   }
   if (k != 0){
       printf("final eigenvalues:");
       for(partID=0;partID< locross_pairs.getSize();partID++){
           eigen = locross_pairs[partID].getEigenValues();
           eigen.show();
       }
       if (k ==100) {printf("failed to make matrix positive definite... really?\n"); return 1;}
   }

    Trianglix<C> trgtinv = ntrgt.mkInverse();
    SparseTrianglix<C> cur;
    double best_LL;
    this->toOne();
    cur = (*this);

    Tuple<C> cross_deter; cross_deter.setSize(locross_pairs.getSize()); ExOp::toZero(cross_deter);
    Tuple<C> cross_trace; cross_trace.setSize(locross_pairs.getSize()); ExOp::toZero(cross_trace);

    Tuple<C> cross_deter_ac; cross_deter_ac.setSize(locross_pairs.getSize()); ExOp::toZero(cross_deter_ac);
    Tuple<C> cross_trace_ac; cross_trace_ac.setSize(locross_pairs.getSize()); ExOp::toZero(cross_trace_ac);

    myHashmap< uint32_t, Tuple<C> > cross_acyclic_comp;
    Tuple<C> cross_acyclic;

    cross_deter += norm_constant; // includes normalization for comparison
    C deter;
    C trace;
    C acyclic;
    ExOp::toZero(deter); ExOp::toZero(acyclic); ExOp::toOne(trace); trace *= target.getSize();

    WeightElem<Tuple<double ,3u> ,2u> crossstat;
    Tuple<double ,3u> crossstat_input;


    Tuple<double ,2u> tmpdoubles;
    ExOp::toZero(details_input.k);
    ExOp::toZero(details_input.d);
    details_input.d[2] -= target.getSize() * locross_pairs.getSize();
    details_input.d[1] += norm_constant * locross_pairs.getSize();
    for(j= 0; j < target.getSize(); j++){
        k = j | (j << 16);
        for(i= 0; i < locross_pairs.getSize(); i+=2){
            crossstat_input[0] = log(cross_P[i][k] = ExOp::mkInverse(locross_pairs[i][j]) );
            cross_deter[i] += crossstat_input[0];
            crossstat_input[1] = log(cross_P[i|1][k] = ExOp::mkInverse(locross_pairs[i|1][j]) );
            cross_deter[i|1] += crossstat_input[1];
            details_input.d[1] += crossstat_input[0];

            cross_trace[i]   += (crossstat_input[0] = cross_P[i][k] * locross_pairs[i|1][j]);
            cross_trace[i|1] += (crossstat_input[1] = cross_P[i|1][k] * locross_pairs[i][j]);
            details_input.d[2] += crossstat_input[0];
        }
    }

    crossstat.toZero();

    for(i= 0; i < locross_pairs.getSize(); i++){
        crossstat_input[0] = cross_deter[i];
        crossstat_input[1] = cross_trace[i];
        crossstat_input[2] = cross_deter[i] - cross_trace[i];
        crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
    }
    details_input.d[4] = sqrt(crossstat.getVar()[2]);
    best_LL = details_input.d[1] - details_input.d[2];
    details.push_back(details_input);
    details_input.d[3] = details_input.d[2] - details_input.d[1];

 //   crossstat_input = crossstat.getMean();
 //   printf("0\t%e-%e (%e-%e  ", deter, trace, crossstat_input[0], crossstat_input[1]);
 //   crossstat_input = crossstat.getVar();
 //   printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

    for(k = 1; k < totsize; k++){
        for(input.d =k; input.d < (k << 16); input.d+= 0x10000){
            input.k = ntrgt.data[(input.d >> 16) + ((k * (k+1)) >> 1)];
            input.k = log(1.0f - input.k * input.k);
            heap_acyclic_merges.insert(input);
        }
    }
    C tmp2[6];

    printf("start!\n");fflush(stdout);

    while((!heap_acyclic_merges.isEmpty())||(!heap_cyclic_merges.isEmpty())){
        #ifdef Rcpp_hpp
            Rcpp::checkUserInterrupt();
        #endif // Rcpp_hpp

        if ((!heap_cyclic_merges.isEmpty())&&((heap_acyclic_merges.isEmpty())||(heap_acyclic_merges.top().k > heap_cyclic_merges.top().k))){
            //printf("Current\n");
            //((Trianglix<C>)*this).show();
            //printf("Difference\n");
            //(this->mkInverse() - ntrgt).show();

            if (heap_cyclic_merges.top().k <= -100.0f){
                printf("tried to add edge that did not get computed properly...\n");
                break;
            }

            details_input.k[0] = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
            details_input.k[1] = heap_cyclic_merges.top().d.k[0] >> 16;
            details_input.d[0] -= heap_cyclic_merges.top().k;


      //      printf("Cyclic %i,%i is %e\n", heap_cyclic_merges.top().d.k[0] & 0xFFFF, heap_cyclic_merges.top().d.k[0] >> 16, -heap_cyclic_merges.top().k);
            partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0], &partit[1], true, &partit[2]);

            details_input.d[9] = partit[0].getSize();
            for(i=0;i<partit[0].getSize();i++){
                if (this->isCyclicFID(this->getFirstID(partit[0][i]))) continue;
                for(j=0;j<partit[0].getSize();j++){
                    l = (partit[0][j] < partit[0][i]) ? partit[0][i] | (partit[0][j] << 16) : partit[0][j] | (partit[0][i] << 16);
                    if ((k = cross_acyclic_comp.find(l)) != 0xFFFFFFFF){
                        for(l=0;l<locross_pairs.getSize();l++){
                            cross_trace_ac[l] -= cross_acyclic_comp.deref(k)[l];
                            cross_deter_ac[l] -= cross_acyclic_comp.deref(k)[l+locross_pairs.getSize()];
                            cross_deter[l] += cross_acyclic_comp.deref(k)[l+locross_pairs.getSize()];
                            if ((l & 1) == 0) details_input.d[2] -= cross_acyclic_comp.deref(k)[l];
                        }
                        cross_acyclic_comp.deref(k).toMemfree();
                        cross_acyclic_comp.erase_from_iterator(k);
                        k = (partit[0][j] < partit[0][i]) ? partit[0][j] + ((partit[0][i] * (partit[0][i]+1)) >> 1) : partit[0][i] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        acyclic += log(1.0 - ntrgt.data[k] * ntrgt.data[k]);
                        deter -= log(1.0 - ntrgt.data[k] * ntrgt.data[k]);
                    }
                }
            }

            this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],cycinput.d.d, ntrgt, details_input.d[5], false,false);//,true);

            //printf("got %i acyclic to fix too!\n", partit[3].getSize());
            // replace by edgeID and triindex

            partit[3].setSize(partit[1].getSize() + partit[0].getSize());
            for(i=0;i<partit[1].getSize();i++) {
                partit[3][i] = (partit[0][partit[1][i] & 0xFFFF] < partit[0][partit[1][i] >> 16]) ? partit[0][partit[1][i] >> 16] | (partit[0][partit[1][i] & 0xFFFF] << 16) : partit[0][partit[1][i] & 0xFFFF] | (partit[0][partit[1][i] >> 16] << 16);
            }
            for(i=0;i<partit[0].getSize();i++) partit[3][i + partit[1].getSize()] = partit[0][i] | (partit[0][i] << 16);
            partit[1].setSize(partit[3].getSize());
            for(i=0;i<partit[3].getSize();i++) partit[1][i] = (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

//            for(i=0;i< partit[3].getSize();i++) partit[1][i+partit[3].getSize()] =  (partit[3][i] >> 16) + (((partit[3][i] & 0xFFFF) * ((partit[3][i] & 0xFFFF)+1)) >> 1);

            // adding edge in crosses!

            details_input.d[6] = 0.0;
            for(l =0; l< locross_pairs.getSize();l+=2){
                this->data.toMemswap(cross_P[l]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0] ,cycinput.d.d, locross_pairs[l], crossstat_input[1], true,true);//,true);,cross_P[l]
                this->data.toMemswap(cross_P[l]);
                cross_deter[l] += crossstat_input[0];
                details_input.d[1] += crossstat_input[0];
                details_input.d[6] += crossstat_input[1];

                this->data.toMemswap(cross_P[l|1]);
                crossstat_input[0] = this->solveWithFictiveEdge2017(heap_cyclic_merges.top().d.k[0],theotherone, locross_pairs[l|1], crossstat_input[1],true,true);//,true);,cross_P[l|1]
                this->data.toMemswap(cross_P[l|1]);

                cross_deter[l|1] += crossstat_input[0];
                //details_input.d[1] += crossstat_input[0];
                //details_input.d[6] += crossstat_input[1];

                crossstat_input.toZero();
                //tmpdoubles.toZero();
                printf("edge is %X, %i and %i, %i\n", heap_cyclic_merges.top().d.k[0] , partit[3].getSize(), cycinput.d.d.getSize(), theotherone.getSize());

                for(i=0;i< partit[3].getSize() - partit[0].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];


                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }
                crossstat_input *= 2.0f; // sym matrix right?
                //tmpdoubles *= 2.0f;

                for(;i<partit[3].getSize();i++){
                    crossstat_input[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    //tmpdoubles[0] += (cycinput.d.d[i] - cross_P[l][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    crossstat_input[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l].data[partit[1][i]];
                    //tmpdoubles[1] += (theotherone[i] - cross_P[l|1][partit[3][i]]) * locross_pairs[l|1].data[partit[1][i]];
                    cross_P[l][partit[3][i]] = cycinput.d.d[i];
                    cross_P[l|1][partit[3][i]] = theotherone[i];
                }

                cross_trace[l] += crossstat_input[0];
                cross_trace[l|1] += crossstat_input[1];
                //cross_deter[l] += tmpdoubles[0];
                //cross_deter[l|1] += tmpdoubles[1];

                details_input.d[2] += crossstat_input[0];// + crossstat_input[1];
                //details_input.d[1] += tmpdoubles[1] + tmpdoubles[0];
            }

            crossstat_input[0] = 0.0f;
            for(i=0;i<partit[3].getSize() - partit[0].getSize();i++){
                if ((j = data.find(partit[3][i])) != 0xFFFFFFFF){
                    crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data.deref(j)) * ntrgt.data[partit[1][i]];
                    data.deref(j) = heap_cyclic_merges.top().d.d[i];
                }else{
          //          printf("Adding that entry! P[%i,%i] = %e\n",partit[3][i] & 0xFFFF, partit[3][i] >> 16, heap_cyclic_merges.top().d.d[i] );
                    this->addEntry(partit[3][i]);
                    crossstat_input[0] += heap_cyclic_merges.top().d.d[i] * ntrgt.data[partit[1][i]];
                    data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
                }
            }
            crossstat_input[0] *= 2.0f;
            for(;i<partit[3].getSize();i++){
                crossstat_input[0] += (heap_cyclic_merges.top().d.d[i] - data[partit[3][i]]) * ntrgt.data[partit[1][i]];
                data[partit[3][i]] = heap_cyclic_merges.top().d.d[i];
            }

            trace += crossstat_input[0];
            deter -= heap_cyclic_merges.top().k - crossstat_input[0];

            // get a fresh taste within cycle group! at last
       //     printf("Heuristic quandidates Cycle\n"); fflush(stdout);
            // finding new candidates within cycle or leading to new cycles

            Vector<uint32_t> edges;
            i=0;
            k = part_id[getFirstID(partit[0][0])];
            j = part_data[k][0];
            do{bestedge_buf[i++] = j; j = next_id[j];
            }while(j != part_data[k][0]);
            partit[1].toMemfree();
            edges.push_back(part_data[k][3]);
            for(;edges[0] != 0xFFFFFFFF;edges[0] = attrib[ edges[0]]){
                partit[1].push_back(i);
                k = part_id[getFirstID(edges[0] & 0xFFFF)];
                edges.push_back(part_data[k][3]);
                j = part_data[k][0];
                do{bestedge_buf[i++] = j; j = next_id[j];
                }while(j != part_data[k][0]);
                k = part_id[getFirstID(partit[0][0])];
                while(true){
                    if (edges.last() == 0xFFFFFFFF) {
                        edges.pop_back();
                        if ((k = edges.getSize()) == 1) break;
                        k = part_id[getFirstID(edges[k-2] >>16 )];
                        edges.last() = attrib[edges.last()];
                    }else{
                        j = part_id[getFirstID(edges.last() & 0xFFFF)];
                        if (j != k){
                            k = part_data[j][0];
                            do{bestedge_buf[i++] = k; k = next_id[k];
                            }while(k != part_data[j][0]);
                            k = part_id[getFirstID(edges.last() >> 16)];
                            edges.push_back(part_data[j][3]);
                        }else edges.last() = attrib[edges.last()];
                    }
                }
            }
            bestedge_size = i;
            edges.toMemfree();
            if (partit[1].getSize() == 0) partit[1].push_back(bestedge_size);
            for(i=0;i < partit[1][0];i++){
                for(j=i+1;j<bestedge_size;j++){
                    l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                    if (this->hasEntry(l)) continue;
                    heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                }
            }
            for(k=1;k<partit[1].getSize();k++){
                for(;i < partit[1][k];i++){
                    for(j=partit[1][k];j<bestedge_size;j++){
                        l = (bestedge_buf[i] < bestedge_buf[j]) ? bestedge_buf[j] | (bestedge_buf[i] << 16) : bestedge_buf[i] | (bestedge_buf[j] << 16);
                        if (this->hasEntry(l)) continue;
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(trgtinv(bestedge_buf[i], bestedge_buf[j])), l));
                    }
                }
            }

            if (bestedge_size > heuristic_merges.getSize()) bestedge_size = heuristic_merges.getSize();
            for(i=0;i<bestedge_size;) bestedge_buf[i++] = heuristic_merges.pop().d;
            heuristic_merges.toMemfree();
            dirtyID++;
            for(i=0;i<partit[0].getSize();i++) dirty_buf[partit[0][i]] = dirtyID;

            if (chkmatch){
                i = heap_cyclic_merges.top().d.k[0] & 0xFFFF;
                j = heap_cyclic_merges.top().d.k[0] >> 16;
                i = (i < j) ? i + ((j * (j+1)) >> 1) : j + ((i * (i+1)) >> 1);
                match[ (chkmatch->data[i] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }
            // makes sure next cyclic merge is not dirty, and not of a too large size!
            heap_cyclic_merges.pop();

            //printf("Current\n");((Trianglix<C>)*this).show();printf("Difference\n");(this->mkInverse() - ntrgt).show();
        }else{
            k = heap_acyclic_merges.top().d;
            //if (wouldEdgeBeInCyclicGroup(k)) LFH_exit("ImpPPPPosible");
            details_input.d[9] = 0;

            details_input.k[0] = k & 0xFFFF;
            details_input.k[1] = k >> 16;
            details_input.d[0] -= heap_acyclic_merges.top().k;

            coor[0] = k & 0xFFFF;
            coor[1] = (k >> 16);
       //     printf("Acyclic %i,%i is %e\n", coor[0], coor[1], -heap_acyclic_merges.top().k);

            // find best candidates for *new* cyclic edges!!
            partit[0] = this->getConnectedSet(coor[0]);
            partit[1] = this->getConnectedSet(coor[1]);
      //      printf("Heuristic quandidates Acycle\n"); fflush(stdout);
            if ((partit[0].getSize() > 1)&&(partit[1].getSize() > 1)){
                for(i=0;i<partit[0].getSize();i++){
                    for(j=0;j<partit[1].getSize();j++){
                        heuristic_merges.insert(KeyElem<double, uint32_t>(-fabs(target(partit[0][i], partit[1][j])), partit[0][i] | (partit[1][j] << 16)));
                    }
                }
                bestedge_size = partit[0].getSize() + partit[1].getSize() - 1;
                for(i=0;i<bestedge_size;){
                    k = (heuristic_merges.top().d >> 16);
                    j = (heuristic_merges.top().d & 0xFFFF);
                    if (k > j){
                        if ((coor[0] != k)||(coor[1] != j)) bestedge_buf[i++] = (heuristic_merges.top().d >> 16) | (heuristic_merges.top().d << 16);
                    }else{
                        if ((coor[0] != j)||(coor[1] != k)) bestedge_buf[i++] = heuristic_merges.top().d;
                    }
                    heuristic_merges.pop();
                }
                heuristic_merges.toMemfree();
            }else if (partit[0].getSize() == 1){
                bestedge_size = partit[1].getSize() -1;
                for(i=0,j=0;i<partit[1].getSize();i++){
                    if (partit[1][i] == coor[1]) j++;
                    else bestedge_buf[i-j] = (partit[1][i] < partit[0][0]) ? partit[0][0] | (partit[1][i] << 16) : partit[1][i] | (partit[0][0] << 16);
                }
            }else{
                bestedge_size = partit[0].getSize() -1;
                for(i=0,j=0;i<partit[0].getSize();i++){
                    if (partit[0][i] == coor[0]) j++;
                    else bestedge_buf[i-j] = (partit[0][i] < partit[1][0]) ? partit[1][0] | (partit[0][i] << 16) : partit[0][i] | (partit[1][0] << 16);
                }
            }
            k = coor[1] + ((coor[0] * (coor[0] + 1)>>1));
            // does not need to update !!invP!! lols

            //tmp2[3] = ntrgt.data[k];
            //tmp2[2] = tmp2[3] * tmp2[3];
            //tmp2[1] = -tmp2[2] +1.0f;
            //tmp2[0] = tmp2[2] / tmp2[1];

            /*
            partit = this->getConnectedSet(coor[0]);
            invup = this->solveInPartition(partit, coor[0]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];

            partit = this->getConnectedSet(coor[1]);
            invup = this->solveInPartition(partit,coor[1]);
            for(i=0;i<partit.getSize();i++) invP[partit[i]] -= invup[i] * invup[i] * tmp2[0];*/

            //(*this)[coor[0]] += tmp2[0];
            //(*/this)[coor[1]] += tmp2[0];
            //(*this)(coor) = -ntrgt.data[k] / tmp2[1]; // ADDING THAT EDGE!
            (*this)(coor) = 0.0f; // ADDING THAT EDGE!
            //trace += heap_acyclic_merges.top().k - tmp2[0];
            //deter -= heap_acyclic_merges.top().k;
            acyclic -= heap_acyclic_merges.top().k;

            cross_acyclic.setSize(locross_pairs.getSize()*2);
            for(i =0; i< locross_pairs.getSize();i+= 2){
                tmp2[0] = locross_pairs[i][coor[0]] * locross_pairs[i][coor[1]];
                tmp2[1] = locross_pairs[i].data[k] * locross_pairs[i].data[k];
                cross_acyclic[locross_pairs.getSize() + i] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                details_input.d[1] += cross_acyclic[locross_pairs.getSize() + i];
                cross_deter_ac[i] += cross_acyclic[locross_pairs.getSize() + i];
                tmp2[5] = tmp2[0] - tmp2[1];
                tmp2[4] = tmp2[1] / tmp2[5];
                tmp2[0] = locross_pairs[i|1][coor[0]] * locross_pairs[i|1][coor[1]];
                tmp2[1] = locross_pairs[i|1].data[k] * locross_pairs[i|1].data[k];
                cross_acyclic[locross_pairs.getSize() + (i|1)] = log(tmp2[0]) - log(tmp2[0] - tmp2[1]);
                //details_input.d[1] += cross_acyclic[locross_pairs.getSize() + (i|1)];
                cross_deter_ac[i|1] += cross_acyclic[locross_pairs.getSize() + (i|1)] ;
                tmp2[3] = tmp2[0] - tmp2[1];
                tmp2[2] = tmp2[1] / tmp2[3];

                tmp2[0] = tmp2[4] / locross_pairs[i][coor[0]];
                cross_acyclic[i] = tmp2[0] * locross_pairs[i|1][coor[0]];
                //cross_P[i][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[4] / locross_pairs[i][coor[1]];
                cross_acyclic[i] += tmp2[0] * locross_pairs[i|1][coor[1]];
                //cross_P[i][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -locross_pairs[i].data[k] / tmp2[5];
                cross_acyclic[i] += tmp2[0] * locross_pairs[i|1].data[k] *2.0;
                //cross_P[i][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i] += cross_acyclic[i];



                tmp2[0] = tmp2[2] / locross_pairs[i|1][coor[0]];
                cross_acyclic[i|1] = tmp2[0] * locross_pairs[i][coor[0]];
                //cross_P[i|1][coor[0] | (coor[0] << 16)] += tmp2[0];
                tmp2[0] = tmp2[2] / locross_pairs[i|1][coor[1]];
                cross_acyclic[i|1] += tmp2[0] * locross_pairs[i][coor[1]];
                //cross_P[i|1][coor[1] | (coor[1] << 16)] += tmp2[0];
                tmp2[0] = -locross_pairs[i|1].data[k] / tmp2[3];
                cross_acyclic[i|1] += tmp2[0] * locross_pairs[i].data[k] *2.0;
                //cross_P[i|1][coor[0] | (coor[1] << 16)] = tmp2[0];
                cross_trace_ac[i|1] += cross_acyclic[i|1];
                details_input.d[2] += cross_acyclic[i];// + cross_acyclic[i|1];
            }
            cross_acyclic_comp[heap_acyclic_merges.top().d].toMemmove(cross_acyclic);

            if (chkmatch){
                match[ (chkmatch->data[coor[1] + ((coor[0] * (coor[0] + 1)>>1))] != 0) ? 2 : 3]++;
                details_input.d[7] = ((double)match[2]) / match[0];
                details_input.d[8] = ((double)match[3]) / match[1];
            }else {details_input.d[7] = details_input.d[8] = 0.0;}

            for(heap_acyclic_merges.pop();!heap_acyclic_merges.isEmpty();heap_acyclic_merges.pop()){
                if (!wouldEdgeBeInCyclicGroup(heap_acyclic_merges.top().d)) break;
            }
        }

        for(;!heap_cyclic_merges.isEmpty();heap_cyclic_merges.pop()){
            partit[0] = getWouldBeCyclicPartition(heap_cyclic_merges.top().d.k[0]);
            if ((maxcyclic != 0)&&(partit[0].getSize() > maxcyclic)) continue;
            for(i=0;i<partit[0].getSize();i++){
                if (dirty_buf[partit[0][i]] > heap_cyclic_merges.top().d.k[1]) break;
            }
            if (i == partit[0].getSize()) break;
        }

        details_input.d[3] += details_input.d[1] -  details_input.d[2];

        if (best_LL < details_input.d[1] -  details_input.d[2]){
            best_LL = details_input.d[1] -  details_input.d[2];
            cur = (*this);
        }

        crossstat.toZero();
        for(i= 0; i < locross_pairs.getSize(); i++){
            crossstat_input[0] = cross_deter[i] + cross_deter_ac[i];
            crossstat_input[1] = cross_trace[i] + cross_trace_ac[i];
            crossstat_input[2] = crossstat_input[0] - crossstat_input[1];
            crossstat += WeightElem<Tuple<double ,3u> ,2u >(crossstat_input);
        }

        details_input.d[4] = sqrt(crossstat.getVar()[2]);
        details_input.d[10] = clock() - init_t;
        details.push_back(details_input);
        details_input.d[3] = details_input.d[2] - details_input.d[1];
    //    printf("Start fill cyclic quandidates\n"); fflush(stdout);
        cycinput.d.k[1] = dirtyID;
        for(i=0;i<bestedge_size;i++){
            cycinput.d.k[0] = bestedge_buf[i];
            cycinput.k = -this->solveWithFictiveEdge2017(bestedge_buf[i],cycinput.d.d, ntrgt,crossstat_input[0] ,false,false); // ,false
            heap_cyclic_merges <<= cycinput;
        }
    //    printf("End fill cyclic quandidates %i \n", heap_cyclic_merges.getSize()); fflush(stdout);
/*  checking zone!



        // check trace is behaving!

        crossstat_input[0] =0.0;
        for(i=0; i< data.getSize();i++){
            j = data.deref_key(i); k =(j & 0xFFFF);
            k = (j >> 16) + ((k * (k+1)) >> 1);
            crossstat_input[0] += data.deref(i) * ntrgt.data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        this->wrDeterminant(crossstat_input[1]);
        printf("trace equal: %e ?= %e\n", trace, crossstat_input[0]);
        printf("deter equal: %e ?= %e\n", deter, log(crossstat_input[1]));

        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[0].getSize();i++){
                j = cross_P[0].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[0].deref(i) * locross_pairs[1].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[0], crossstat_input[0]);
        crossstat_input[0] =0.0;
        for(i=0; i< cross_P[1].getSize();i++){
                j = cross_P[1].deref_key(i);  k =(j & 0xFFFF);
                k = (j >> 16) + ((k * (k+1)) >> 1);
                crossstat_input[0] += cross_P[1].deref(i) * locross_pairs[0].data[k] * (((j >> 16) == (j & 0xFFFF)) ? 1.0 : 2.0) ;
        }
        printf("trace should be equal: %e ?= %e\n", cross_trace[1], crossstat_input[0]);

        printf("Acyclic traces: %e & %e\n", cross_trace_ac[0], cross_trace_ac[1]);
        printf("Acyclic determinant: %e & %e\n", cross_deter_ac[0], cross_deter_ac[1]);


        this->data.toMemswap(cross_P[0]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[0]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[0]);
        this->data.toMemswap(cross_P[1]);
        this->wrDeterminant(crossstat_input[0]);
        this->data.toMemswap(cross_P[1]);
        printf("Cylic Determinant %e vs %e\n", log(crossstat_input[0]), cross_deter[1]);


        crossstat_input = crossstat.getMean();
        printf("%i\t%e-%e (%e-%e  ", data.getSize() - totsize, deter, trace, crossstat_input[0], crossstat_input[1]);
        crossstat_input = crossstat.getVar();
        printf("+- %e-%e)\n", sqrt(crossstat_input[0]), sqrt(crossstat_input[1]));

        printf("%i: %e\n", data.getSize() - totsize, deter + acyclic - trace);
        crossstat_input.toZero();
        for(i=0; i< cross_P.getSize();i++){
            crossstat_input[0] += cross_deter_ac[i] - cross_trace_ac[i];
            crossstat_input[1] += cross_deter[i] - cross_trace[i];
        }
        printf("vs: %e\n", (crossstat_input[0]  + crossstat_input[1]) / cross_P.getSize());
        */

        printf("%i edge considered so far\n", data.getSize() - totsize);
        if ((data.getSize() - totsize) > lambda) break;
    }

    (*this) = cur;

    // formally add acyclic edges
    for(i=0;i<part_data.getSize();i++){
        if (part_data[i][1] > part_data[i][2]){
            partit[0] = this->getPartition(i);
            for(j=0;j<partit[0].getSize();j++){
                for(k=0;k<neigh[partit[0][j]].getSize();k++){
                    if ((getFirstID(neigh[partit[0][j]][k]) != partit[0][0])||(partit[0][j] < neigh[partit[0][j]][k])){
                        l = neigh[partit[0][j]][k] > partit[0][j] ? partit[0][j] + ((neigh[partit[0][j]][k] * (neigh[partit[0][j]][k]+1)) >> 1) :  neigh[partit[0][j]][k] + ((partit[0][j] * (partit[0][j]+1)) >> 1);
                        tmp2[3] = ntrgt.data[l];
                        tmp2[2] = tmp2[3] * tmp2[3];
                        tmp2[1] = -tmp2[2] +1.0f;
                        tmp2[0] = tmp2[2] / tmp2[1];
                        (*this)[partit[0][j]] += tmp2[0];
                        (*this)[neigh[partit[0][j]][k]] += tmp2[0];
                        coor[0] = neigh[partit[0][j]][k] > partit[0][j] ? neigh[partit[0][j]][k] | (partit[0][j] << 16) : partit[0][j] | (neigh[partit[0][j]][k] << 16);
                        (*this).data[coor[0]] = -ntrgt.data[l] / tmp2[1];
                    }
                }
            }
        }
    }


    /*printf("Error\n");
    ntrgt -= this->mkInverse();
    for(i=0,k=0;i<target.getSize();i++,k++){
        for(j=0;j<i;j++,k++){
            if (data.find(i |(j << 16)) == 0xFFFFFFFF) ntrgt.data[k] = 0.0f;
        }
    }*/
   // ntrgt.show();

    // scale precision matrix to match prior to normalization
    cycinput.d.d.setSize(target.getSize());
    for(i=0;i<target.getSize();i++) cycinput.d.d[i] = ExOp::mkPowInvInt(target.data[(i * (i+3))>> 1] , -2); // 1/sqrt
    for(i=0;i<data.getSize();i++) data.deref(i) *= cycinput.d.d[(data.deref_key(i)& 0xFFFF)] * cycinput.d.d[(data.deref_key(i) >> 16)];

   // printf("Final matrix\n");
   // ((Trianglix<C>)*this).show();

    for(i=0;i<details.getSize();i++){
        details[i].d[1] /= locross_pairs.getSize();
        details[i].d[2] /= locross_pairs.getSize();
        details[i].d[3] /= locross_pairs.getSize();
    }
return 0;}

LFHTEMP void SparseTrianglix<C>::wrTrace(C& fout) const{
    if (next_id.getSize() == 0) {ExOp::toZero(fout); return;}
    fout = data[0];
    for(int i=1;i<next_id.getSize();i++) fout += (*this)[i];
}
LFHTEMP template<class O> void SparseTrianglix<C>::wrTrace(O& fout) const{
    if (next_id.getSize() == 0) {ExOp::toZero(fout); return;}
    O tmp;
    ExOp::wrTrace(data[0], tmp);
    fout = tmp;
    for(int i=1;i<next_id.getSize();i++) {ExOp::wrTrace(data[i], tmp); fout += tmp;}
}

LFHTEMP void SparseTrianglix<C>::solveAllAcyclicEdges(const Trianglix<C> &target, bool diagisone){
    uint32_t i,j,k,l,m;
    Tuple<uint32_t> partit;
    double tmp[4];
    if (diagisone){
        for(i=0;i<part_data.getSize();i++){
            if (part_data[i][1] > part_data[i][2]){
                partit = this->getPartition(i);
                for(j=0;j<partit.getSize();j++){
                    for(k=0;k<neigh[partit[j]].getSize();k++){
                        if ((getFirstID(neigh[partit[j]][k]) != partit[0])||(partit[j] < neigh[partit[j]][k])){
                            l = neigh[partit[j]][k] > partit[j] ? partit[j] + ((neigh[partit[j]][k] * (neigh[partit[j]][k]+1)) >> 1) :  neigh[partit[j]][k] + ((partit[j] * (partit[j]+1)) >> 1);
                            tmp[3] = target.data[l];
                            tmp[2] = tmp[3] * tmp[3];
                            tmp[1] = -tmp[2] +1.0f;
                            tmp[0] = tmp[2] / tmp[1];
                            (*this)[partit[j]] += tmp[0];
                            (*this)[neigh[partit[j]][k]] += tmp[0];
                            m = neigh[partit[j]][k] > partit[j] ? neigh[partit[j]][k] | (partit[j] << 16) : partit[j] | (neigh[partit[j]][k] << 16);
                            (*this).data[m] = -target.data[l] / tmp[1];
        }   }   }   }   }
    }else{
        for(i=0;i<part_data.getSize();i++){
            if (part_data[i][1] > part_data[i][2]){
                partit = this->getPartition(i);
                for(j=0;j<partit.getSize();j++){
                    for(k=0;k<neigh[partit[j]].getSize();k++){
                        if ((getFirstID(neigh[partit[j]][k]) != partit[0])||(partit[j] < neigh[partit[j]][k])){
                            l = neigh[partit[j]][k] > partit[j] ? partit[j] + ((neigh[partit[j]][k] * (neigh[partit[j]][k]+1)) >> 1) :  neigh[partit[j]][k] + ((partit[j] * (partit[j]+1)) >> 1);
                            tmp[0] = target[neigh[partit[j]][k]] * target[partit[j]];
                            tmp[1] = target.data[l]* target.data[l];
                            tmp[3] = tmp[0] - tmp[1];
                            tmp[2] = tmp[1] / tmp[3];
                            (*this)[partit[j]] += tmp[2] / target[partit[j]];
                            (*this)[neigh[partit[j]][k]] += tmp[2] / target[neigh[partit[j]][k]];
                            m = (neigh[partit[j]][k] > partit[j]) ? neigh[partit[j]][k] | (partit[j] << 16) : partit[j] | (neigh[partit[j]][k] << 16);
                            (*this).data[m] = -target.data[l] / tmp[3];
    }   }   }   }   }   }
}

#ifdef Rcpp_hpp

LFHTEMP template<class RCL> void SparseTrianglix<C>::wrMatrix(arma::Mat<RCL> &fout) const{
    fout.zeros(next_id.getSize(),next_id.getSize());
    uint32_t i,j,k;
    uint32_t ite;
    for(i=0;i<data.getSize();i++){
        j = data.deref_key(i);
        k = (j >> 16); j &= 0xFFFF;
        fout.at(j,k) = data.deref(i);
        if (k != j) fout.at(k,j) = data.deref(i);
    }
}

LFHTEMP void SparseTrianglix<C>::wrMatrix(Rcpp::NumericMatrix &fout) const{
    fout = Rcpp::NumericMatrix(next_id.getSize(),next_id.getSize());
    uint32_t i,j,k;
    uint32_t ite;
    for(i=0;i<data.getSize();i++){
        j = data.deref_key(i);
        k = (j >> 16); j &= 0xFFFF;
        fout(j,k) = data.deref(i);
        if (k != j) fout(k,j) = data.deref(i);
    }
}

LFHTEMP template<class RCL> void SparseTrianglix<C>::wrDetMatrix(arma::Mat<RCL> &fout) const{
    fout.zeros(next_id.getSize(),next_id.getSize());
    uint32_t i,j,k;
    uint32_t ite;
    for(i=0;i<data.getSize();i++){
        j = data.deref_key(i);
        k = (j >> 16); j &= 0xFFFF;
        ExOp::wrDeterminant(data.deref(i), fout.at(j,k));
        if (k != j) fout.at(k,j) = fout.at(j,k);
    }
}


LFHTEMP template<class RCL> void SparseTrianglix<C>::wrSubMatrix(arma::Mat<RCL> &fout, const Tuple<uint32_t> &partit) const{
    fout.zeros(partit.getSize(),partit.getSize());
    myHashmap<uint16_t, uint16_t> backward;
    uint32_t i,j,k;
    for(i=0;i<partit.getSize();i++) {
        backward[partit[i]] = i;
        fout.at(i,i) = (*this)[partit[i]];
    }
    uint32_t ite;
    for(i=0;i<partit.getSize();i++){
        for(j= 0;j<neigh[partit[i]].getSize();j++) {
            if (neigh[partit[i]][j] > partit[i]) continue;
            if ((ite = backward.find(neigh[partit[i]][j])) != 0xFFFFFFFF){
                k = backward.deref(ite);
                ite = data.find(partit[i]  | (neigh[partit[i]][j] << 16));
                fout.at(i,k) = data.deref(ite);
                fout.at(k,i) = data.deref(ite);
            }
        }
    }
}

LFHTEMP Trianglix<double> SparseTrianglix<C>::mkInverse()const{Trianglix<double> fout;
    arma::Mat<double> armat[2];
    this->wrMatrix(armat[0]);
    if (!arma::inv(armat[1], armat[0])) {fout.setSize(next_id.getSize()); ExOp::toUndefined(fout);}
    else fout.rdMatrix(armat[1]);
    return fout;
}

LFHTEMP void SparseTrianglix<C>::toRegularizedInverseOf(const Trianglix<C> &target){
    // needs to minimize J = log |P| - tr(SP),
    // at optimum, dJ/dP = 0, for all non-zero entries
    // that is: 0 = S - P^-1, only for all non-zero entries of P
    // hence, if P is highly sparse, individual components can be solved individually
    uint32_t i,x,loop;
    uint32_t y[2];
    CurvatureSearchScope cssc;
    double* guess;
    double* g_target;
    double* deriv;
    double value;
    uint32_t* trindexes;
    uint32_t* axes;
    int totsize;
    int diagosize;
    Trianglix<double> sub, dainv;

    PartitionIterator pp(*this);
    if (pp.first()) do{
        //if (pp.getSize() < 4){
            if (pp.getSize() < 3){
                if (pp.getSize() == 1){
                    this->acxEntry(pp.part[0] | (pp.part[0] << 16)) = ExOp::mkInverse(target(pp.part[0],pp.part[0]));
                }else{
                    C denum = ExOp::mkInverse(target(pp.part[0], pp.part[0]) * target(pp.part[1], pp.part[1]) - ExOp::mkTrjuProd(target(pp.part[0], pp.part[1])));
                    this->acxEntry(pp.part[0] | (pp.part[0] << 16)) = target(pp.part[1], pp.part[1]) * denum;
                    this->acxEntry(pp.part[1] | (pp.part[0] << 16)) = -target(pp.part[0], pp.part[1]) * denum;
                    this->acxEntry(pp.part[1] | (pp.part[1] << 16)) = target(pp.part[0], pp.part[0]) * denum;
                }
            }else{
                totsize = part_data[pp.part[0]][2] + part_data[pp.part[0]][1];
                diagosize = part_data[pp.part[0]][2];

                guess = new double[totsize];
                deriv = new double[totsize];
                g_target = new double[totsize];
                trindexes = new uint32_t[totsize];
                axes = new uint32_t[diagosize * 2];


                myHashmap<uint32_t,uint32_t> backmap;
                uint32_t ite;

                for(i= 0;i<part_data[pp.part[0]][1];i++){
                    trindexes[i+diagosize] = (i  * (i + 3))>> 1;
                  //  printf("trindex[%i] = %i\n",i+diagosize, trindexes[i+diagosize] );
                    backmap[pp.part[i]] = i;
                    g_target[i+diagosize] = target(pp.part[i],pp.part[i]);
                }

                for(i=0,x=0;x< pp.getSize();x++){
                    y[0] = backmap[pp.part[x]];
                    for(ite = 0;ite < neigh[pp.part[x]].getSize();ite++){
                        if (pp.part[x] < neigh[pp.part[x]][ite]) continue;
                        y[1] = backmap[neigh[pp.part[x]][ite]];
                        g_target[i] =  target(pp.part[x],neigh[pp.part[x]][ite]);
                        axes[(i<<1)] = y[0] + diagosize;
                        axes[(i<<1)|1] = y[1] + diagosize;
                        trindexes[i] = (y[0] > y[1]) ? y[1] + ((y[0]*(y[0]+1))>> 1) : y[0] + ((y[1]*(y[1]+1))>> 1);
                       // printf("trindex[%i] = %i [a%i,%i]\n",i, trindexes[i], axes[(i<<1)] , axes[(i<<1)|1] );
                        i++;
                    }
                }

                // initial guess;
                sub.setSize(part_data[pp.part[0]][1]); sub.toZero();
                for(i=0;i<totsize;i++) sub.data[trindexes[i]] = g_target[i];

                //dainv = sub.mkInverse();
                arma::Mat<double> armat[2];
                sub.wrMatrix(armat[0]);
                armat[1] = armat[0].i();
                dainv.rdMatrix(armat[1]);

                for(i=diagosize;i<totsize;i++){
                    guess[i] = sqrt(sub.data[trindexes[i]]);
                }
                for(;i<diagosize;i++){
                    guess[i] = sub.data[trindexes[i]];
                }
                cssc.init(0.01f, totsize);
                for(loop=0;loop<250;loop++){
                    // value = log |P| - tr(SP),

                    value = 0.0f;
                    for(i=0;i<diagosize;i++) {sub.data[trindexes[i]] = guess[axes[(i<<1)]]*guess[axes[(i<<1)|1]]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
                    for(;i<totsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * guess[i]* guess[i];} //


                    //sub.log_determinant(&i);
                    //if (i != 0) value = nan("");
                    //printf("Eigen\n");
                    //ExOp::show(sub.getEigenValues());
                    //dainv = sub.mkInverse();
                    arma::Mat<double> armat[2];
                    sub.wrMatrix(armat[0]);

                    armat[1] = armat[0].i();
                    dainv.rdMatrix(armat[1]);
                    value += arma::log_det(armat[0]).real();
                    //(dainv * sub).show();
                    //sub.show();
                    //dainv.show();

                    for(i=diagosize;i<totsize;i++) deriv[i] = guess[i] * (dainv.data[trindexes[i]] - g_target[i]);
                    for(i=0;i<diagosize;i++) {
                        deriv[i] = tanh(guess[i]);
                        deriv[i] = guess[axes[(i<<1)]]*guess[axes[(i<<1)|1]]* (dainv.data[trindexes[i]] - g_target[i]) * (1.0f - deriv[i] * deriv[i]);
                        deriv[axes[(i<<1)]] += (sub.data[trindexes[i]] / guess[axes[(i<<1)]]) * (dainv.data[trindexes[i]] - g_target[i]);
                        deriv[axes[(i<<1)|1]] += (sub.data[trindexes[i]] / guess[axes[(i<<1)|1]]) * (dainv.data[trindexes[i]] - g_target[i]);
                    }

                    value *= 0.5f;
                    //if (cssc.updateAscent(value,guess,deriv)) break;
                    if (loop < 10){
                        if (cssc.updateAscent(value,guess,deriv)) break;
                    }else{
                        if (loop == 10) cssc.init(0.01f, totsize);
                        if (cssc.checkDerivative(value,guess,deriv) < 0.0000000001f) break;
                    }
                    for(i=0;i<diagosize;i++) printf("%e\t", guess[i]);
                    for(;i<totsize;i++) printf("%e\t", guess[i]*guess[i]);

                    printf("\nF[%i] = %f\n", loop, value);
                }


                for(i=0,x=0;x< pp.getSize();x++){
                    for(ite = 0;ite < neigh[pp.part[x]].getSize();ite++){
                        if (pp.part[x] < neigh[pp.part[x]][ite]) continue;
                        this->acxEntry((pp.part[x] << 16) | neigh[pp.part[x]][ite] ) = guess[i++];
                    }
                }
                for(;i<totsize;i++){
                    this->acxEntry(pp.part[i-diagosize] | (pp.part[i-diagosize] << 16)) = guess[i] * guess[i];
                }

                for(i=0;i<totsize;i++){
                    sub.data[trindexes[i]] = g_target[i];
                    dainv.data[trindexes[i]] = i < diagosize ? guess[i]  : guess[i] * guess[i];
                }
                sub.show();
                dainv.show();
                printf("is ID!\n");

                (dainv * sub).show();

                delete[](guess); delete[](deriv); delete[](g_target); delete[](trindexes); delete[](axes);





                // either fully connected or 1 missing
                /*if (part_data[pp.part[0]][2] == 3){
                    // full matrix;
                }else{

                }
                if (this->hasEntry(pp.part[0] | (pp.part[1] << 16))){
                    if (this->hasEntry(pp.part[0] | (pp.part[2] << 16))){
                        if (this->hasEntry(pp.part[1] | (pp.part[2] << 16))){

                        }else{

                        }
                    }else{

                    }
                }else{

                }*/


            }
       // }else{
       //     exit(1);
       // }
    }while(pp.next());
}

LFHTEMP void SparseTrianglix<C>::addEdgeToRegInv(int i, int j, const Trianglix<C> &target){
    if (getFirstID[i] == getFirstID[j]){
        // edge is within a component

    }else{
        // edge is across two connected component

    }
    uint32_t pos = (i < j) ? (j |  (i << 16)) : (i |  (j << 16));
    this->addEntry(pos);
}

LFHTEMP void SparseTrianglix<C>::remEdgeToRegInv(int i, int j, const Trianglix<C> &target){
    if (getFirstID[i] == getFirstID[j]){
        // edge is within a component

    }else{
        // edge is across two connected component

    }

    uint32_t pos = (i < j) ? (j |  (i << 16)) : (i |  (j << 16));
    this->removeEntry(pos);
}

LFHTEMP Tuple<C> SparseTrianglix<C>::solveInPartition(Tuple<uint32_t> partit, uint32_t column)const{ Tuple<C> fout; // solve Px = e_column
    fout.setSize(partit.getSize());
    int i=0;
    while(partit[i] !=column) i++;
    arma::Mat<C> mat;
    wrSubMatrix(mat, partit);
    arma::Col<C> vec = arma::zeros<arma::Col<C> >(partit.getSize());
    vec.at(i) = 1.0f;
    arma::Col<C> sol = arma::solve(mat,vec);
    for(i=0;i<partit.getSize();i++) fout[i] = sol.at(i);
    return fout;
}

LFHTEMP void SparseTrianglix<C>::wrDeterminant(C& fout) const{ // assumes C is double!
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
    arma::Mat<double> damat;
    this->wrMatrix(damat);
    fout = arma::det(damat);
}
LFHTEMP template<class O> void SparseTrianglix<C>::wrDeterminant(O& fout) const{
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
    arma::Mat<double> damat;
    this->wrDetMatrix(damat);
    fout = arma::det(damat);
}

LFHTEMP void SparseTrianglix<C>::wrLogDeterminant(C& fout) const{ // assumes C is double!
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
    arma::Mat<double> damat;
    this->wrMatrix(damat);
    double sign;
    arma::log_det(fout, sign, damat);
    if (sign < 0) ExOp::toUndefined(fout);
}
LFHTEMP template<class O> void SparseTrianglix<C>::wrLogDeterminant(O& fout) const{
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
    arma::Mat<double> damat;
    this->wrDetMatrix(damat);
    double sign;
    arma::log_det(fout, sign, damat);
    if (sign < 0) ExOp::toUndefined(fout);
}

LFHTEMP void SparseTrianglix<C>::solveConstainedML(uint32_t partID, const Trianglix<C> &target){
    Tuple<uint32_t> wouldedge;
    Tuple<uint32_t> partit = getPartition(partID, &wouldedge);
    uint32_t tsize = wouldedge.getSize() + partit.getSize();
    bool debug = false;
    bool insist = true;
    double error_out;
    Tuple<double> guess; guess.setSize(tsize);
    Tuple<double> deriv; deriv.setSize(tsize);
    double* ondiago = &(guess[0]) + wouldedge.getSize();
    double* onderiv = &(deriv[0]) + wouldedge.getSize();

    Tuple<double> g_target; g_target.setSize(tsize);

    double value;
    int ofdiagosize = wouldedge.getSize();
    int loop,i;
    uint32_t x;
    uint32_t y[2];
    double tmp[4];


    myHashmap<uint32_t,uint32_t> backmap;
    uint32_t ite;

    Tuple<uint32_t> trindexes; trindexes.setSize(tsize);
    Trianglix<double> sub, dainv; sub.setSize(partit.getSize()); sub.toZero();





    for(i= 0;i<partit.getSize();i++){
        trindexes[i+ofdiagosize] = (i  * (i + 3))>> 1;
        ondiago[i] = pow(target[partit[i]], -0.5);
        onderiv[i] = -target[partit[i]] * target[partit[i]]; // = -((P^-1)_{ii})^2 == -target[partit[i]]^2 at start!...
        g_target[i+ofdiagosize] = target[partit[i]];
    }

    // (P^-1 - S) * tanh" + (pii^2+2pij^2+pjj^2)*(tanh')^2
    for(i= 0;i<wouldedge.getSize();i++){
        trindexes[i] = (wouldedge[i] >> 16) > (wouldedge[i] & 0xFFFF) ? (wouldedge[i] & 0xFFFF) + (((wouldedge[i] >>16) * ((wouldedge[i] >>16)+1))>>1) : (wouldedge[i] >> 16) + (((wouldedge[i] &0xFFFF) * ((wouldedge[i] &0xFFFF)+1))>>1);
        guess[i] = 0.0f;
        deriv[i] = 1.0f;
        g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);
    }

    if ((!ExOp::isValid(guess))||(!ExOp::isValid(deriv))||(!ExOp::isValid(g_target))){
        printf("got NAs... nose!\n");
        ExOp::show(guess);
        ExOp::show(deriv);
        ExOp::show(g_target);
        this->show();
        exit(1);
    }

    // get double derivative to start!

    CurvatureSearchScope cssc;
    cssc.init(tsize, deriv());
    arma::Mat<double> armat[2];

    Tuple<double,0u> eigenv;

    for(loop=0;loop<25000;loop++){
        // value = log |P| - tr(SP),  // unconstrained maximum is P = S^-1  => = -log |S| - n

        //if ((loop & 63) == 63) cssc.updateRandom(guess());
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) {sub.data[trindexes[i]] = ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
        for(;i<tsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * guess[i]* guess[i];} //



        sub.wrMatrix(armat[0]);
        if (!arma::inv(armat[1], armat[0])){
            ExOp::toUndefined(value);
            if (loop == 0) {sub.show();printf("started iterative search with Nan, fix this!\n");exit(1);}
            cssc.updateGotNAN(guess());
        }else{

            dainv.rdMatrix(armat[1]);





            //sub.show();
            //dainv = sub.mkInverse();
            //(dainv * sub).show();


            for(i=ofdiagosize;i<tsize;i++) deriv[i] = guess[i] * (dainv.data[trindexes[i]] - g_target[i]);
            for(i=0;i<ofdiagosize;i++) {
                deriv[i] = ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i]>> 16]* (dainv.data[trindexes[i]] - g_target[i]) * d_tanh(guess[i]);
                onderiv[wouldedge[i] &0xFFFF] += (sub.data[trindexes[i]] / ondiago[wouldedge[i] &0xFFFF]) * (dainv.data[trindexes[i]] - g_target[i]);
                onderiv[wouldedge[i]>> 16] += (sub.data[trindexes[i]] / ondiago[wouldedge[i]>> 16]) * (dainv.data[trindexes[i]] - g_target[i]);
            }

            if (loop == 0){
                if (debug) {
                    printf("start subP value =%e ldet =%e\n", value, tmp[0]);
                    sub.show();
                    printf("start deriv\n");
                    deriv.show();
                }
                if ((!ExOp::isValid(value))||(!ExOp::isValid(tmp[0]))){
                    // beyond repair...
                    printf("got the impossible here!");
                    printf("start subP value =%e ldet =%e\n", value, tmp[0]);
                    sub.show();
                    printf("start deriv\n");
                    deriv.show();

                    return;

                }

            }

            tmp[0] = 0.0f; // arma::log_det(armat[0]).real();
            eigenv = sub.getEigenValues();
            for(i=0;i<eigenv.getSize();i++) if (eigenv[i] <= 0.0) break; else tmp[0] += log(eigenv[i]);

            value += tmp[0];
            value *= 0.5f; // there is a 2.0 factor for deriv and dderiv anyways

        //if (loop < 10){
            if ((i < eigenv.getSize())||(!ExOp::isValid(value))) cssc.updateGotNAN(guess());
            else if (cssc.updateAscent(value,guess(),deriv())) {
                if (value > 0.00001f){
                    printf("after %i: got stuck with error=%e\n", loop , value);
                    if (!insist) break;
                }else break;
            }
         }
        //}else{
        //    if (loop == 10) cssc.init(tsize, 0.01f );
        //    if (cssc.checkDerivative(value,guess,deriv) < 0.0000000001f) break;
        // }
        /*for(i=0;i<ofdiagosize;i++) printf("%e\t", guess[i]);
        for(;i<tsize;i++) printf("%e\t", guess[i]*guess[i]);

        printf("\nF[%i] = %f\n", loop, value);*/
    }
    cssc.wrFinalGuess(guess());

    for(i=0;i<wouldedge.getSize();i++) {
        uint32_t j = (partit[wouldedge[i] &0xFFFF] > partit[wouldedge[i] >> 16]) ? (partit[wouldedge[i] >> 16] << 16) | partit[wouldedge[i] &0xFFFF] : (partit[wouldedge[i] &0xFFFF] << 16) | partit[wouldedge[i]>> 16 ];
        this->data[j] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i] >> 16];
    }

    for(i=0;i<partit.getSize();i++) this->data[partit[i] | (partit[i] << 16)] = ondiago[i] * ondiago[i];

    for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16];
    for(;i<tsize;i++) sub.data[trindexes[i]] = guess[i]*guess[i];

    sub.wrMatrix(armat[0]);
    if (!arma::inv(armat[1], armat[0])) {ExOp::toUndefined(value); ExOp::toUndefined(error_out);}
    else{
        tmp[0] = arma::log_det(armat[0]).real();
        dainv.rdMatrix(armat[1]);
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);
        for(;i<tsize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);

        error_out = value;
        if (debug) {
            printf("Current sub-P:\n");
            sub.show();
            for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = g_target[i];
            for(;i<tsize;i++) sub.data[trindexes[i]] = g_target[i];
            printf("Current target:\n");
            sub.show();
            printf("Current sub-error:\n");
            (((Trianglix<C>)dainv)- sub).show();
        }
    }
return;}
// solve modified system, where "datamap" contains the current values of P that minimizes det(P) - tr(target * P)
LFHTEMP double SparseTrianglix<C>::solveWithMissingEdge(uint32_t missingEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error_out, bool return_det_only, bool debug) const{
    Tuple<uint32_t> wouldedge;
    Tuple<uint32_t> partit = getEdgePartition(missingEdge, &wouldedge);
    int loop,i;
    double tmp[4];
    if (wouldedge.getSize() == 1){ // removing an acyclic edge... use close form!
        loop = (partit[wouldedge[i] &0xFFFF] > partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] + ((partit[wouldedge[i] &0xFFFF]*(partit[wouldedge[i] &0xFFFF]+1))>>1) : partit[wouldedge[i] &0xFFFF] +  ((partit[wouldedge[i] >> 16]*(partit[wouldedge[i] >> 16]+1)) >> 1);
        return log_x_p1( target.data[loop] * target.data[loop] / (target.data[missingEdge & 0xFFFF] * target.data[missingEdge >> 16]));
    }

    uint32_t x;
    for( x=0;x < wouldedge.getSize(); x++) if ((partit[wouldedge[x] & 0xFFFF] | (partit[wouldedge[x] >>16] << 16) ) == missingEdge) break;
    wouldedge.removeAt(x);


    uint32_t tsize = wouldedge.getSize() + partit.getSize();
    foreach.printf("%i,%i are sizes\n", partit.getSize(), wouldedge.getSize()); fflush(stdout);
    Tuple<double> guess; guess.setSize(tsize);
    Tuple<double> deriv; deriv.setSize(tsize);
    double* ondiago = &(guess[0]) + wouldedge.getSize();
    double* onderiv = &(deriv[0]) + wouldedge.getSize();
    Tuple<double> g_target; g_target.setSize(tsize);

    arma::Col<C> eigval;

    Tuple<double> datmptmp;

    double value, start_value, start_det;
    int ofdiagosize = wouldedge.getSize();


    uint32_t y[2];

    double sillyvalgrind;
    double daerror;

    myHashmap<uint32_t,uint32_t> backmap;
    uint32_t ite;


    Tuple<uint32_t> trindexes; trindexes.setSize(tsize);
    Trianglix<double> sub, dainv; sub.setSize(partit.getSize()); sub.toZero();
    Trianglix<double> startsub;

    for(i= 0;i<partit.getSize();i++){
        trindexes[i+ofdiagosize] = (i  * (i + 3))>> 1;
        g_target[i+ofdiagosize] = datamap[partit[i] | (partit[i] << 16)];
    }

    for(i= 0;i<wouldedge.getSize();i++){ // add edge if acyclic *and* not the fictive edge
        trindexes[i] = (wouldedge[i] >> 16) > (wouldedge[i] & 0xFFFF) ? (wouldedge[i] & 0xFFFF) + (((wouldedge[i] >>16) * ((wouldedge[i] >>16)+1))>>1) : (wouldedge[i] >> 16) + (((wouldedge[i] &0xFFFF) * ((wouldedge[i] &0xFFFF)+1))>>1);
        loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
        g_target[i] = datamap[loop];
    }

    foreach.printf("filling!\n");
    datmptmp = g_target;
    for(i= 0;i<partit.getSize();i++){
        ondiago[i] = sqrt(g_target[i+ofdiagosize]);
        onderiv[i] = -target[partit[i]] * target[partit[i]]; // = -((P^-1)_{ii})^2 == -target[partit[i]]^2 at start!...
        g_target[i+ofdiagosize] = target[partit[i]];
    }
    for(i= 0;i<wouldedge.getSize()-1;i++){
        guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] )  );
        if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
        deriv[i] = (g_target[i]*g_target[i] + 1.0) * ExOp::mkSquare(d_tanh(guess[i])); // assumes dJ(P)/dP = 0
        if (!ExOp::isValid(deriv[i])){
            g_target[i] = fabs(g_target[i]);
            deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
            if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
        }
        g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);
    }
    // the new edge:
    guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] ) );
    if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
    deriv[i] = (- 1.0 - g_target[i]*g_target[i] ) * ExOp::mkSquare(d_tanh(guess[i]));
    // deriv[i] += d2_tanh(guess[i]) * P^{-1}{ij} - S^{ij} // is zero if acyclic edge... but if we are here...
    if (!ExOp::isValid(deriv[i])){
        g_target[i] = fabs(g_target[i]);
        deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
        if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
    }
    g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);

    if ((!ExOp::isValid(guess))||(!ExOp::isValid(deriv))||(!ExOp::isValid(g_target))){
        if (!ExOp::isValid(guess)) {foreach.printf("got NAs in guess... nose!\n"); ExOp::show(guess);}
        if (!ExOp::isValid(deriv)) {foreach.printf("got NAs in deriv... nose!\n"); ExOp::show(deriv);}
        if (!ExOp::isValid(g_target)) {foreach.printf("got NAs in target... nose!\n"); ExOp::show(g_target);}
        if (debug) foreach.printf("out solve!\n");
        return -100.0f;
    }

    CurvatureSearchScope cssc;
    cssc.init(tsize, deriv());
    arma::Mat<double> armat[2];
    double sign;
    // get double derivative to start!
    int maxloop = 2500;
    foreach.printf("Start looping!\n");


return 0;}

/*LFHTEMP double SparseTrianglix<C>::solveWithFictiveEdge(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error_out, Tuple<double> &wouldbe_P, bool return_det_only, bool debug) const{
    Tuple<uint32_t> wouldedge;
    foreach.printf("in solve!\n"); fflush(stdout);
    Tuple<uint32_t> partit = getWouldBeCyclicPartition(fictiveEdge, &wouldedge);
    uint32_t tsize = wouldedge.getSize() + partit.getSize();
    foreach.printf("%i,%i are sizes\n", partit.getSize(), wouldedge.getSize()); fflush(stdout);
    Tuple<double> guess; guess.setSize(tsize);
    Tuple<double> deriv; deriv.setSize(tsize);
    double* ondiago = &(guess[0]) + wouldedge.getSize();
    double* onderiv = &(deriv[0]) + wouldedge.getSize();

    Tuple<double> g_target; g_target.setSize(tsize);
    arma::Col<C> eigval;

    Tuple<double> datmptmp;

    double value, start_value, start_det;
    int ofdiagosize = wouldedge.getSize();
    int loop,i;
    uint32_t x;
    uint32_t y[2];
    double tmp[4];
    double sillyvalgrind;
    double daerror;

    myHashmap<uint32_t,uint32_t> backmap;
    uint32_t ite;


    Tuple<uint32_t> trindexes; trindexes.setSize(tsize);
    Trianglix<double> sub, dainv; sub.setSize(partit.getSize()); sub.toZero();
    Tuple<unsigned int> wasacyclic; wasacyclic.setSize(partit.getSize());
    Trianglix<double> startsub;

    //if (debug) printf("start diago prior to acyclic:\n");
    for(i= 0;i<partit.getSize();i++){
        trindexes[i+ofdiagosize] = (i  * (i + 3))>> 1;
        wasacyclic[i] = isCyclicFID(getFirstID(partit[i])) ? 1 : 0;
        g_target[i+ofdiagosize] = datamap[partit[i] | (partit[i] << 16)];
        //if (debug) printf("%e\t", g_target[i+ofdiagosize]);
    }
    //if (debug) printf("\nstart offdiago elems:\n");
    for(i= 0;i<wouldedge.getSize();i++){ // add edge if acyclic *and* not the fictive edge
        trindexes[i] = (wouldedge[i] >> 16) > (wouldedge[i] & 0xFFFF) ? (wouldedge[i] & 0xFFFF) + (((wouldedge[i] >>16) * ((wouldedge[i] >>16)+1))>>1) : (wouldedge[i] >> 16) + (((wouldedge[i] &0xFFFF) * ((wouldedge[i] &0xFFFF)+1))>>1);

        if (i == wouldedge.getSize() - 1) g_target[i] = 0.0f;
        else if (wasacyclic[wouldedge[i] &0xFFFF] * wasacyclic[wouldedge[i] >> 16] == 0){
            // add acyclic edge
            loop = (partit[wouldedge[i] &0xFFFF] > partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] + ((partit[wouldedge[i] &0xFFFF]*(partit[wouldedge[i] &0xFFFF]+1))>>1) : partit[wouldedge[i] &0xFFFF] +  ((partit[wouldedge[i] >> 16]*(partit[wouldedge[i] >> 16]+1)) >> 1);
            tmp[0] = target[partit[wouldedge[i] &0xFFFF]] * target[partit[wouldedge[i] >> 16]] ;
            tmp[1] = target.data[loop] * target.data[loop];
            tmp[3] = tmp[0] - tmp[1];
            tmp[2] = tmp[1] / tmp[3];
            g_target[ofdiagosize + (wouldedge[i] >> 16)] += tmp[2] / target[partit[wouldedge[i] >> 16]];
            g_target[ofdiagosize + (wouldedge[i] & 0xFFFF)] += tmp[2] / target[partit[wouldedge[i] & 0xFFFF]];
            g_target[i] = -target.data[loop] / tmp[3];
            //if (debug) printf("AE(%i,%i): %e\n", wouldedge[i] &0xFFFF, wouldedge[i] >> 16,g_target[i]);
        }else{
            loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
            g_target[i] = datamap[loop];
            //if (debug) printf("CE(%i,%i): %e\n", wouldedge[i] &0xFFFF, wouldedge[i] >> 16,datamap[loop]);
        }
    }
    foreach.printf("filling!\n");
    datmptmp = g_target;
    for(i= 0;i<partit.getSize();i++){
       // backmap[partit[i]] = i;
        ondiago[i] = sqrt(g_target[i+ofdiagosize]);
        onderiv[i] = -target[partit[i]] * target[partit[i]]; // = -((P^-1)_{ii})^2 == -target[partit[i]]^2 at start!...
        g_target[i+ofdiagosize] = target[partit[i]];
    }

    // (P^-1 - S) * tanh" + (pii^2+2pij^2+pjj^2)*(tanh')^2
    for(i= 0;i<wouldedge.getSize()-1;i++){
        // updating entries that are already inserted edges
        guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] )  );
        if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
        deriv[i] = (g_target[i]*g_target[i] + 1.0) * ExOp::mkSquare(d_tanh(guess[i])); // assumes dJ(P)/dP = 0
        if (!ExOp::isValid(deriv[i])){
            g_target[i] = fabs(g_target[i]);
            deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
            if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
        }
        g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);
    }
    // the new edge:
    guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] ) );
    if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
    deriv[i] = (- 1.0 - g_target[i]*g_target[i] ) * ExOp::mkSquare(d_tanh(guess[i]));
    // deriv[i] += d2_tanh(guess[i]) * P^{-1}{ij} - S^{ij} // is zero if acyclic edge... but if we are here...
    if (!ExOp::isValid(deriv[i])){
        g_target[i] = fabs(g_target[i]);
        deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
        if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
    }
    g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);

    if ((!ExOp::isValid(guess))||(!ExOp::isValid(deriv))||(!ExOp::isValid(g_target))){
        if (!ExOp::isValid(guess)) {foreach.printf("got NAs in guess... nose!\n"); ExOp::show(guess);}
        if (!ExOp::isValid(deriv)) {foreach.printf("got NAs in deriv... nose!\n"); ExOp::show(deriv);}
        if (!ExOp::isValid(g_target)) {foreach.printf("got NAs in target... nose!\n"); ExOp::show(g_target);}
        //printf("current:\n");
        //datamap.show();
        if (debug) foreach.printf("out solve!\n");
        return -100.0f;
    }

    // get double derivative to start!

    CurvatureSearchScope cssc;
    cssc.init(tsize, deriv());
    arma::Mat<double> armat[2];
    double sign;
    // get double derivative to start!
    int maxloop = 2500;
    foreach.printf("Start looping!\n");
    for(loop=0;loop<maxloop;loop++){
        // value = log |P| - tr(SP),  // unconstrained maximum is P = S^-1  => = -log |S| - n
        Rcpp::checkUserInterrupt(); // hello!!
        //if ((loop & 63) == 63) cssc.updateRandom(guess());
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) {sub.data[trindexes[i]] = ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
        for(;i<tsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * guess[i]* guess[i];} //



        sub.wrMatrix(armat[0]);
        foreach.printf("armacall...\n");
        if (!arma::inv(armat[1], armat[0])){
            foreach.printf("armacall... done!\n");
            ExOp::toUndefined(value);
            if (loop == 0) break;
        }else{
            foreach.printf("armacall... done!\n");

            dainv.rdMatrix(armat[1]);
            foreach.printf("armacall2...\n");
            arma::eig_sym( eigval, armat[0] );
            foreach.printf("armacall2... done!\n");
            tmp[0] = 0.0;
            for(i=0;i < dainv.getSize(); i++){
                if (eigval[i] <= 0.0) break;
                tmp[0] += log(eigval[i]);
            }
            if ((debug)&&(loop == 0)){
                printf("Determinant at start:\n");
                for(i=0;i < dainv.getSize();i++) printf("%e%c",eigval[i], (i+1 ==  dainv.getSize()) ? '\n' : '\t');
                printf("For %i edges in:\n", wouldedge.getSize());
                partit.show();
                //sub.show();
            }
            if ((i < dainv.getSize())||(!ExOp::isValid(tmp[0]))){
                if (loop == 0) {
                    foreach.printf("problematic determinant at start:\n");
                    for(i=0;i < dainv.getSize();i++) foreach.printf("%e%c",eigval[i], (i+1 ==  dainv.getSize()) ? '\n' : '\t');
                    foreach.printf("For %i edges in:\n", wouldedge.getSize());
                    partit.show();
                    //sub.show();
                    printf("dathat start values:\n");
                    datmptmp.show();
                    break;
                }else ExOp::toUndefined(tmp[0]); // went too far!
            }


            for(i=ofdiagosize;i<tsize;i++) deriv[i] = guess[i] * (dainv.data[trindexes[i]] - g_target[i]);
            for(i=0;i<ofdiagosize;i++) {
                deriv[i] = ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i]>> 16]* (dainv.data[trindexes[i]] - g_target[i]) * d_tanh(guess[i]);
                onderiv[wouldedge[i] &0xFFFF] += (sub.data[trindexes[i]] / ondiago[wouldedge[i] &0xFFFF]) * (dainv.data[trindexes[i]] - g_target[i]);
                onderiv[wouldedge[i]>> 16] += (sub.data[trindexes[i]] / ondiago[wouldedge[i]>> 16]) * (dainv.data[trindexes[i]] - g_target[i]);
            }

            if (loop == 0){
                if (debug) foreach.printf("start subP value =%e ldet =%e\n", value, tmp[0]);
                if (!ExOp::isValid(value)) break;
                start_value = (return_det_only) ? tmp[0] : value + tmp[0];
            }

            value += tmp[0];
            value *= 0.5f; // there is a 2.0 factor for deriv and dderiv anyways


            if (cssc.updateAscent(value,guess(),deriv()))  break;
            // if (debug) printf("F[%i] = %e\n",loop, value);
         }
        //}else{
        //    if (loop == 10) cssc.init(tsize, 0.01f );
        //    if (cssc.checkDerivative(value,guess,deriv) < 0.0000000001f) break;
        // }
    }
    foreach.printf("End looping!\n");
    if (maxloop == loop) foreach.printf("warning, failed to converge\n");
    wouldbe_P.setSize(tsize);
    if (loop == 0){
        foreach.printf("Started a search with invalid guess\n");
        foreach.printf("\nstart subP value =%e ldet =%e\n", value, tmp[0]);
        guess.show();
        g_target.show();
        foreach.printf("out solve!\n"); fflush(stdout);
        return -100.0;
    }
    if (debug) {
        foreach.printf("final eigenvals:\n");
        for(i=0;i < dainv.getSize(); i++) foreach.printf("%e\t",eigval[i]);
        foreach.printf("\nFor %i edges in:\n", wouldedge.getSize());
        partit.show();
        //sub.show();
    }
    cssc.wrFinalGuess(guess());




    for(i=0;i<wouldedge.getSize();i++) wouldbe_P[i] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i] >> 16];
    for(i=0;i<partit.getSize();i++) wouldbe_P[i+ofdiagosize] = ondiago[i] * ondiago[i];
    if (debug) {foreach.printf("finalize search!\n"); fflush(stdout); wouldbe_P.show();}

    for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16];
    for(;i<tsize;i++) sub.data[trindexes[i]] = guess[i]*guess[i];

    sub.wrMatrix(armat[0]);
            foreach.printf("armacall3...\n");
    if (!arma::inv(armat[1], armat[0])) {foreach.printf("armacall...done\n");ExOp::toUndefined(value); ExOp::toUndefined(error_out);}
    else{
        foreach.printf("armacall3...done\n");
        tmp[0] = arma::log_det(armat[0]).real();
        foreach.printf("armacall3...doneagain\n");
        dainv.rdMatrix(armat[1]);
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);
        for(;i<tsize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);

        error_out = value;
    }

    if (return_det_only) start_value = tmp[0] - start_value;
    else start_value = 2.0f * cssc.getLastValue() - start_value;

    if (debug) {
        foreach.printf("Proposed update (from->to) %e (s=%i,%i nbstep = %i, Error: %e)\n",  start_value, partit.getSize(), wouldedge.getSize(), loop, value);
    }

    foreach.printf("out solve!\n"); fflush(stdout);
return ((ExOp::isValid(value))&&(value < 0.1f)) ? start_value : -100.0;}*/

LFHTEMP double SparseTrianglix<C>::solveWithFictiveEdge2017(uint32_t fictiveEdge, Tuple<double> &wouldbe_P, const Trianglix<C> &target, double &error_out, bool return_det_only, bool debug) const{
    Tuple<uint32_t> wouldedge;
    Tuple<uint32_t> partit = getWouldBeCyclicPartition(fictiveEdge, &wouldedge);
    uint32_t tsize = wouldedge.getSize() + partit.getSize();

    Tuple<double> guess; guess.setSize(tsize);
    Tuple<double> deriv; deriv.setSize(tsize);
    double* ondiago = &(guess[0]) + wouldedge.getSize();
    double* onderiv = &(deriv[0]) + wouldedge.getSize();

    Tuple<double> g_target; g_target.setSize(tsize);

    double value, start_value, start_det;
    int ofdiagosize = wouldedge.getSize();
    int loop,i;
    uint32_t x;
    uint32_t y[2];
    double tmp[4];
    double sillyvalgrind;
    double daerror;

    myHashmap<uint32_t,uint32_t> backmap;
    uint32_t ite;

    Tuple<uint32_t> trindexes; trindexes.setSize(tsize);
    Trianglix<double> sub, dainv; sub.setSize(partit.getSize()); sub.toZero();
    Trianglix<double> startsub;

    Tuple<unsigned int> wasacyclic; wasacyclic.setSize(partit.getSize());
    for(i= 0;i<partit.getSize();i++){
        trindexes[i+ofdiagosize] = (i  * (i + 3))>> 1;
        wasacyclic[i] = isCyclicFID(getFirstID(partit[i])) ? 1 : 0;
        g_target[i+ofdiagosize] = (*this)[partit[i]];
    }
    for(i= 0;i<wouldedge.getSize();i++){ // add edge if acyclic *and* not the fictive edge
        trindexes[i] = (wouldedge[i] >> 16) > (wouldedge[i] & 0xFFFF) ? (wouldedge[i] & 0xFFFF) + (((wouldedge[i] >>16) * ((wouldedge[i] >>16)+1))>>1) : (wouldedge[i] >> 16) + (((wouldedge[i] &0xFFFF) * ((wouldedge[i] &0xFFFF)+1))>>1);

        if (i == wouldedge.getSize()-1) g_target[i] = 0.0f;
        else if (wasacyclic[wouldedge[i] &0xFFFF] * wasacyclic[wouldedge[i] >> 16] == 0){
            // add acyclic edge
            loop = (partit[wouldedge[i] &0xFFFF] > partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] + ((partit[wouldedge[i] &0xFFFF]*(partit[wouldedge[i] &0xFFFF]+1))>>1) : partit[wouldedge[i] &0xFFFF] +  ((partit[wouldedge[i] >> 16]*(partit[wouldedge[i] >> 16]+1)) >> 1);
            tmp[0] = target[partit[wouldedge[i] &0xFFFF]] * target[partit[wouldedge[i] >> 16]] ;
            tmp[1] = target.data[loop] * target.data[loop];
            tmp[3] = tmp[0] - tmp[1];
            tmp[2] = tmp[1] / tmp[3];
            g_target[ofdiagosize + (wouldedge[i] >> 16)] += tmp[2] / target[partit[wouldedge[i] >> 16]];
            g_target[ofdiagosize + (wouldedge[i] & 0xFFFF)] += tmp[2] / target[partit[wouldedge[i] & 0xFFFF]];
            g_target[i] = -target.data[loop] / tmp[3];
        }else{
            loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
            g_target[i] = data[loop];
        }
    }
    for(i= 0;i<partit.getSize();i++){
       // backmap[partit[i]] = i;
        ondiago[i] = sqrt(g_target[i+ofdiagosize]);
        onderiv[i] = -target[partit[i]] * target[partit[i]]; // = -((P^-1)_{ii})^2 == -target[partit[i]]^2 at start!...
        g_target[i+ofdiagosize] = target[partit[i]];
    }

    // (P^-1 - S) * tanh" + (pii^2+2pij^2+pjj^2)*(tanh')^2
    for(i= 0;i<wouldedge.getSize()-1;i++){
        // updating entries that are already inserted edges
        guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] )  );
        if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
        deriv[i] = (g_target[i]*g_target[i] + 1.0) * ExOp::mkSquare(d_tanh(guess[i])); // assumes dJ(P)/dP = 0
        g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);
    }
    // the new edge:
    guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] ) );
    if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
    deriv[i] = (- 1.0 - g_target[i]*g_target[i] ) * ExOp::mkSquare(d_tanh(guess[i]));
    // deriv[i] += d2_tanh(guess[i]) * P^{-1}{ij} - S^{ij} // is zero if acyclic edge... but if we are here...
    g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);

    // get double derivative to start!

    CurvatureSearchScope cssc;
    cssc.init(tsize, deriv());
    arma::Mat<double> armat[2];
    wouldbe_P.setSize(tsize); // hurray, wont crash... but that's unitilizalized all right
    double sign;
    arma::Col<C> eigval;

    int maxloop = 2500;
    for(loop=0;loop<maxloop;loop++){
        // value = log |P| - tr(SP),

        //if ((loop & 63) == 63) cssc.updateRandom(guess());
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) {sub.data[trindexes[i]] = ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
        for(;i<tsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * sub.data[trindexes[i]];} //

        if (loop == 0) startsub = sub;

        sub.wrMatrix(armat[0]);
        if (!arma::inv(armat[1], armat[0])){
            ExOp::toUndefined(value);
            if (loop == 0) {printf("warning: Started iterative search singular matrix!\n");break;}
        }else{

            dainv.rdMatrix(armat[1]);
            //arma::log_det(,sign,armat[0]);
            arma::eig_sym( eigval, armat[0] );
            tmp[0] = 0.0;
            for(i=0;i < dainv.getSize(); i++){
                if (eigval[i] <= 0.0) break;
                tmp[0] += log(eigval[i]);
            }
            if (i < dainv.getSize()){
                if (loop == 0) break;
                else ExOp::toUndefined(tmp[0]); // went too far!
            }
            //sub.show();
            //dainv = sub.mkInverse();
            //(dainv * sub).show();


            for(i=ofdiagosize;i<tsize;i++) deriv[i] = guess[i] * (dainv.data[trindexes[i]] - g_target[i]);
            for(i=0;i<ofdiagosize;i++) {
                deriv[i] = ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i]>> 16]* (dainv.data[trindexes[i]] - g_target[i]) * d_tanh(guess[i]);
                onderiv[wouldedge[i] &0xFFFF] += (sub.data[trindexes[i]] / ondiago[wouldedge[i] &0xFFFF]) * (dainv.data[trindexes[i]] - g_target[i]);
                onderiv[wouldedge[i]>> 16] += (sub.data[trindexes[i]] / ondiago[wouldedge[i]>> 16]) * (dainv.data[trindexes[i]] - g_target[i]);
            }

            if (loop == 0){
                // if (debug) { printf("start subP value =%e ldet =%e\n", value, tmp[0]); sub.show(); printf("start deriv\n"); deriv.show(); }
                if ((!ExOp::isValid(value))||(!ExOp::isValid(tmp[0]))) break;
                start_det = tmp[0];
                start_value = value;
            }
            value += tmp[0];
            value *= 0.5f; // there is a 2.0 factor for deriv and dderiv anyways
        }

        //if (loop < 10){
            // if (debug) {ExOp::show(guess); printf("%i:%e\n",loop,value);}

            if (cssc.updateAscent(value,guess(),deriv())) break;
        //}else{
        //    if (loop == 10) cssc.init(tsize, 0.01f );
        //    if (cssc.checkDerivative(value,guess,deriv) < 0.0000000001f) break;
        // }
        //for(i=0;i<ofdiagosize;i++) printf("%e\t", guess[i]);
        //for(;i<tsize;i++) printf("%e\t", guess[i]*guess[i]);
        //printf("F[%i] = %e   (log det %e, modif %e)\n", loop, value, tmp[0], cssc.getRelErrorBound());
        //printf("F[%i] = %e\n", loop, value);
    }
    if (loop == 0){
        printf("Started a search with invalid input\n");
        sub.show();
        printf("start deriv\n");
        deriv.show();
        printf("Covar:\n"); dainv.show();
        return -100.0f;
    }
    //printf("F[%i] = %e   (log det %e, modif %e)\n", loop, value, tmp[0], cssc.getRelErrorBound());

    cssc.wrFinalGuess(guess());
    if (loop == maxloop) {
        printf("failed to converge!\n");
        return -100.0;
    }
    for(i=0;i<wouldedge.getSize();i++) wouldbe_P[i] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i] >> 16];
    for(i=0;i<partit.getSize();i++) wouldbe_P[i+ofdiagosize] = ondiago[i] * ondiago[i];

    for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16];
    for(;i<tsize;i++) sub.data[trindexes[i]] = guess[i]*guess[i];

    sub.wrMatrix(armat[0]);
    if (!arma::inv(armat[1], armat[0])) {ExOp::toUndefined(error_out); ExOp::toUndefined(error_out);}
    else{
        tmp[0] = arma::log_det(armat[0]).real();
        dainv.rdMatrix(armat[1]);
        error_out = 0.0f;
        for(i=0;i<ofdiagosize;i++) error_out += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);
        for(;i<tsize;i++) error_out += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);

        if ((error_out)&&((!ExOp::isValid(error_out))||(error_out >= 0.1f))) {
            printf("Got problematic error %e... %e -> %e\n", error_out, start_value, value);
            printf("Got log determinant %e -> %e\n", start_det, tmp[0]);
            printf("Current start-sub-P:\n"); startsub.show();
            printf("Current sub-P:\n"); sub.show();
            for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = g_target[i];
            for(;i<tsize;i++) sub.data[trindexes[i]] = g_target[i];
            printf("Current target:\n"); sub.show();
            printf("Current Inverse:\n"); dainv.show();
            printf("Current Sub-error:\n"); (((Trianglix<C>)dainv)- sub).show();


        }
    }

    if (return_det_only) start_value = tmp[0] - start_value;
    else start_value = 2.0f * cssc.getLastValue() - start_value - start_det;

    if (debug) {
        printf("Proposed update (from->to) %e (s=%i,%i nbstep = %i, Error: %e)\n",  start_value, partit.getSize(), wouldedge.getSize(), loop, error_out);
    }

    /*
    for(i= 0;i<partit.getSize();i++) ondiago[i] = (*this)[partit[i]];
    for(i= 0;i<wouldedge.getSize();i++){
        loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
        guess[i] = data[loop];
    }
    guess.show();
    wouldbe_P.show();*/

    return ((ExOp::isValid(value))&&(value < 0.1f)) ? start_value : -100.0;
}

#else

LFHTEMP Trianglix<double> SparseTrianglix<C>::mkInverse()const{Trianglix<double> fout;
    fout.setSize(next_id.getSize());
    return fout;
}

LFHTEMP Tuple<C> SparseTrianglix<C>::solveInPartition(Tuple<uint32_t> partit, uint32_t column)const{ Tuple<C> fout; // solve Px = e_column
    fout.setSize(partit.getSize());
    fout.toUndefined();
    return fout;
}
LFHTEMP void SparseTrianglix<C>::wrDeterminant(C& fout) const{
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
}
LFHTEMP template<class O> void SparseTrianglix<C>::wrDeterminant(O& fout) const{
    if (next_id.getSize() == 0) {ExOp::toOne(fout); return;}
}
LFHTEMP void SparseTrianglix<C>::wrLogDeterminant(C& fout) const{ // assumes C is double!
    if (next_id.getSize() == 0) {ExOp::toZero(fout); return;}
}
LFHTEMP template<class O> void SparseTrianglix<C>::wrLogDeterminant(O& fout) const{
    if (next_id.getSize() == 0) {ExOp::toZero(fout); return;}
}
LFHTEMP double SparseTrianglix<C>::solveWithMissingEdge(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error, bool return_det_only, bool debug) const{

return 0;}

#endif

LFHTEMP double SparseTrianglix<C>::solveWithFictiveEdge(uint32_t fictiveEdge, const myHashmap<uint32_t, C> &datamap, const Trianglix<C> &target, double &error, Tuple<double> &wouldbe_P,  bool return_det_only, bool debug) const{
    Tuple<uint32_t> wouldedge;
    foreach.printf("in solve!\n"); fflush(stdout);
    Tuple<uint32_t> partit = getWouldBeCyclicPartition(fictiveEdge, &wouldedge);
    uint32_t tsize = wouldedge.getSize() + partit.getSize();
    foreach.printf("%i,%i are sizes\n", partit.getSize(), wouldedge.getSize()); fflush(stdout);
    Tuple<double> guess; guess.setSize(tsize);
    Tuple<double> deriv; deriv.setSize(tsize);
    double* ondiago = &(guess[0]) + wouldedge.getSize();
    double* onderiv = &(deriv[0]) + wouldedge.getSize();

    Tuple<double> g_target; g_target.setSize(tsize);
    Tuple<double> datmptmp;

    double value, start_value, start_det;
    int ofdiagosize = wouldedge.getSize();
    int loop,i;
    uint32_t x;
    uint32_t y[2];
    double tmp[4];
    double sillyvalgrind;
    double daerror;
    double error_out;
    myHashmap<uint32_t,uint32_t> backmap;
    uint32_t ite;


    Tuple<uint32_t> trindexes; trindexes.setSize(tsize);
    Trianglix<double> sub, dainv; sub.setSize(partit.getSize()); sub.toZero();
    Tuple<unsigned int> wasacyclic; wasacyclic.setSize(partit.getSize());
    Trianglix<double> startsub;

    //if (debug) printf("start diago prior to acyclic:\n");
    for(i= 0;i<partit.getSize();i++){
        trindexes[i+ofdiagosize] = (i  * (i + 3))>> 1;
        wasacyclic[i] = isCyclicFID(getFirstID(partit[i])) ? 1 : 0;
        g_target[i+ofdiagosize] = datamap[partit[i] | (partit[i] << 16)];
        //if (debug) printf("%e\t", g_target[i+ofdiagosize]);
    }
    //if (debug) printf("\nstart offdiago elems:\n");
    for(i= 0;i<wouldedge.getSize();i++){ // add edge if acyclic *and* not the fictive edge
        trindexes[i] = (wouldedge[i] >> 16) > (wouldedge[i] & 0xFFFF) ? (wouldedge[i] & 0xFFFF) + (((wouldedge[i] >>16) * ((wouldedge[i] >>16)+1))>>1) : (wouldedge[i] >> 16) + (((wouldedge[i] &0xFFFF) * ((wouldedge[i] &0xFFFF)+1))>>1);

        if (i == wouldedge.getSize() - 1) g_target[i] = 0.0f;
        else if (wasacyclic[wouldedge[i] &0xFFFF] * wasacyclic[wouldedge[i] >> 16] == 0){
            // add acyclic edge
            loop = (partit[wouldedge[i] &0xFFFF] > partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] + ((partit[wouldedge[i] &0xFFFF]*(partit[wouldedge[i] &0xFFFF]+1))>>1) : partit[wouldedge[i] &0xFFFF] +  ((partit[wouldedge[i] >> 16]*(partit[wouldedge[i] >> 16]+1)) >> 1);
            tmp[0] = target[partit[wouldedge[i] &0xFFFF]] * target[partit[wouldedge[i] >> 16]] ;
            tmp[1] = target.data[loop] * target.data[loop];
            tmp[3] = tmp[0] - tmp[1];
            tmp[2] = tmp[1] / tmp[3];
            g_target[ofdiagosize + (wouldedge[i] >> 16)] += tmp[2] / target[partit[wouldedge[i] >> 16]];
            g_target[ofdiagosize + (wouldedge[i] & 0xFFFF)] += tmp[2] / target[partit[wouldedge[i] & 0xFFFF]];
            g_target[i] = -target.data[loop] / tmp[3];
            //if (debug) printf("AE(%i,%i): %e\n", wouldedge[i] &0xFFFF, wouldedge[i] >> 16,g_target[i]);
        }else{
            loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
            g_target[i] = datamap[loop];
            //if (debug) printf("CE(%i,%i): %e\n", wouldedge[i] &0xFFFF, wouldedge[i] >> 16,datamap[loop]);
        }
    }
    /*if (debug) {
        printf("start diago: after acyclic\n");
        for(i= 0;i<partit.getSize();i++) printf("%e\t", g_target[i+ofdiagosize]);
    }*/
    foreach.printf("filling!\n");
    datmptmp = g_target;
    for(i= 0;i<partit.getSize();i++){
       // backmap[partit[i]] = i;
        ondiago[i] = sqrt(g_target[i+ofdiagosize]);
        onderiv[i] = -target[partit[i]] * target[partit[i]]; // = -((P^-1)_{ii})^2 == -target[partit[i]]^2 at start!...
        g_target[i+ofdiagosize] = target[partit[i]];
    }

    // (P^-1 - S) * tanh" + (pii^2+2pij^2+pjj^2)*(tanh')^2
    for(i= 0;i<wouldedge.getSize()-1;i++){
        // updating entries that are already inserted edges
        guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] )  );
        if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
        deriv[i] = (g_target[i]*g_target[i] + 1.0) * ExOp::mkSquare(d_tanh(guess[i])); // assumes dJ(P)/dP = 0
        if (!ExOp::isValid(deriv[i])){
            g_target[i] = fabs(g_target[i]);
            deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
            if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
        }
        g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);
    }
    // the new edge:
    guess[i] = atanh( g_target[i] / (ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16] ) );
    if ((!ExOp::isValid(guess[i]))||(fabs(guess[i]) > 25.0f)) guess[i] = (guess[i] > 0.0f) ? 25.0f : -25.0f; // we dont like plateaus
    deriv[i] = (- 1.0 - g_target[i]*g_target[i] ) * ExOp::mkSquare(d_tanh(guess[i]));
    // deriv[i] += d2_tanh(guess[i]) * P^{-1}{ij} - S^{ij} // is zero if acyclic edge... but if we are here...
    if (!ExOp::isValid(deriv[i])){
        g_target[i] = fabs(g_target[i]);
        deriv[i] = 4.0 * exp( -4.0 * g_target[i] + log(g_target[i])* 2.0);
        if (!ExOp::isValid(deriv[i])) deriv[i] = 0.0f;
    }
    g_target[i] = target(partit[wouldedge[i] &0xFFFF], partit[wouldedge[i] >> 16]);

    if ((!ExOp::isValid(guess))||(!ExOp::isValid(deriv))||(!ExOp::isValid(g_target))){
        if (!ExOp::isValid(guess)) {foreach.printf("got NAs in guess... nose!\n"); ExOp::show(guess);}
        if (!ExOp::isValid(deriv)) {foreach.printf("got NAs in deriv... nose!\n"); ExOp::show(deriv);}
        if (!ExOp::isValid(g_target)) {foreach.printf("got NAs in target... nose!\n"); ExOp::show(g_target);}
        //printf("current:\n");
        //datamap.show();
        if (debug) foreach.printf("out solve!\n");
        return -100.0f;
    }
    /*
    value = 0.0f;
    for(i=0;i<ofdiagosize;i++) {sub.data[trindexes[i]] = ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
    for(;i<tsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * guess[i]* guess[i];}
    sub.wrMatrix(armat[0]);
    arma::eig_sym( eigval, armat[0] );
    tmp[0] = 0.0;
    for(i=0;i < dainv.getSize(); i++){
        if (eigval[i] <= 0.0) break;
        tmp[0] += log(eigval[i]);
    }
    if ((i < dainv.getSize())||(!ExOp::isValid(tmp[0]))){


    }else if (!ExOp::isValid(value)){
        printf("Invalid Start...\n");
        return -100.0f
    }*/

    // get double derivative to start!

    CurvatureSearchScope cssc;
    cssc.init(tsize, deriv());
    Magscaled<double> deter;
    short posdef;
    double sign;
    // get double derivative to start!
    int maxloop = 2500;
    foreach.printf("Start looping!\n");
    for(loop=0;loop<maxloop;loop++){
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) {sub.data[trindexes[i]] = ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16]* tanh(guess[i]); value -= 2.0f * g_target[i] * sub.data[trindexes[i]];}
        for(;i<tsize;i++) {sub.data[trindexes[i]] = guess[i]*guess[i]; value -= g_target[i] * guess[i]* guess[i];} //


        dainv = sub.mkInverse(&deter, &posdef);

        foreach.printf("armacall...\n");
        if (!posdef){
            if ((debug)&&(loop == 0)){
                printf("Not Positive definite at start! Determinant at start: %e * e^%e\n", deter.value , deter.exponent);
            }
            foreach.printf("armacall... done!\n");
            ExOp::toUndefined(value);
            if (loop == 0) break;
        }else{
            foreach.printf("armacall... done!\n");
            tmp[0] = log(deter.value) + deter.exponent;
            for(i=ofdiagosize;i<tsize;i++) deriv[i] = guess[i] * (dainv.data[trindexes[i]] - g_target[i]);
            for(i=0;i<ofdiagosize;i++) {
                deriv[i] = ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i]>> 16]* (dainv.data[trindexes[i]] - g_target[i]) * d_tanh(guess[i]);
                onderiv[wouldedge[i] &0xFFFF] += (sub.data[trindexes[i]] / ondiago[wouldedge[i] &0xFFFF]) * (dainv.data[trindexes[i]] - g_target[i]);
                onderiv[wouldedge[i]>> 16] += (sub.data[trindexes[i]] / ondiago[wouldedge[i]>> 16]) * (dainv.data[trindexes[i]] - g_target[i]);
            }

            if (loop == 0){
                if (debug) foreach.printf("start subP value =%e ldet =%e\n", value, tmp[0]);
                if (!ExOp::isValid(value)) break;
                start_value = (return_det_only) ? tmp[0] : value + tmp[0];
            }

            value += tmp[0];
            value *= 0.5f; // there is a 2.0 factor for deriv and dderiv anyways


            if (cssc.updateAscent(value,guess(),deriv()))  break;
            // if (debug) printf("F[%i] = %e\n",loop, value);
            /*
            printf("curvaturescope start\n");
            if ((i < eigenv.getSize())||(!ExOp::isValid(value))) cssc.updateGotNAN(guess());
            else if (cssc.updateAscent(value,guess(),deriv())) {
                if (value > 0.00001f){
                    printf("after %i: got stuck with error=%e\n", loop , value);
                    if (!insist) break;
                }else break;
            }
            printf("curvaturescope end\n");*/
         }
        //}else{
        //    if (loop == 10) cssc.init(tsize, 0.01f );
        //    if (cssc.checkDerivative(value,guess,deriv) < 0.0000000001f) break;
        // }
        /*for(i=0;i<ofdiagosize;i++) printf("%e\t", guess[i]);
        for(;i<tsize;i++) printf("%e\t", guess[i]*guess[i]);

        printf("\nF[%i] = %f\n", loop, value);*/
    }
    foreach.printf("End looping!\n");
    if (maxloop == loop) foreach.printf("warning, failed to converge\n");
    wouldbe_P.setSize(tsize);
    if (loop == 0){
        foreach.printf("Started a search with invalid guess\n");
        foreach.printf("\nstart subP value =%e ldet =%e\n", value, tmp[0]);
        guess.show();
        g_target.show();
        foreach.printf("out solve!\n"); fflush(stdout);
        return -100.0;
    }
    if (debug) {
        foreach.printf("\nFor %i edges in:\n", wouldedge.getSize());
        partit.show();
        //sub.show();
    }
    cssc.wrFinalGuess(guess());




    for(i=0;i<wouldedge.getSize();i++) wouldbe_P[i] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF] * ondiago[wouldedge[i] >> 16];
    for(i=0;i<partit.getSize();i++) wouldbe_P[i+ofdiagosize] = ondiago[i] * ondiago[i];
    if (debug) {foreach.printf("finalize search!\n"); fflush(stdout); wouldbe_P.show();}

    for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = tanh(guess[i]) * ondiago[wouldedge[i] &0xFFFF]*ondiago[wouldedge[i]>> 16];
    for(;i<tsize;i++) sub.data[trindexes[i]] = guess[i]*guess[i];

    dainv = sub.mkInverse(&deter, &posdef);
    if (posdef != 1) {foreach.printf("armacall...done\n");ExOp::toUndefined(value); ExOp::toUndefined(error_out);}
    else{
        tmp[0] = log(deter.value) + deter.exponent;
        value = 0.0f;
        for(i=0;i<ofdiagosize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);
        for(;i<tsize;i++) value += ExOp::mkSquare(dainv.data[trindexes[i]] - g_target[i]);

        error_out = value;
        /*if (debug) {
            printf("Current sub-P:\n");
            sub.show();
            for(i=0;i<ofdiagosize;i++) sub.data[trindexes[i]] = g_target[i];
            for(;i<tsize;i++) sub.data[trindexes[i]] = g_target[i];
            printf("Current target:\n");
            sub.show();
            printf("Current sub-error:\n");
            (((Trianglix<C>)dainv)- sub).show();
        }*/
    }

    if (return_det_only) start_value = tmp[0] - start_value;
    else start_value = 2.0f * cssc.getLastValue() - start_value;

    if (debug) {
        foreach.printf("Proposed update (from->to) %e (s=%i,%i nbstep = %i, Error: %e)\n",  start_value, partit.getSize(), wouldedge.getSize(), loop, value);
    }

    /*
    for(i= 0;i<partit.getSize();i++) ondiago[i] = (*this)[partit[i]];
    for(i= 0;i<wouldedge.getSize();i++){
        loop = (partit[wouldedge[i] &0xFFFF] < partit[wouldedge[i] >> 16]) ? partit[wouldedge[i] >> 16] |  (partit[wouldedge[i] &0xFFFF]<<16) : partit[wouldedge[i] &0xFFFF] |  (partit[wouldedge[i] >> 16]<<16);
        guess[i] = data[loop];
    }
    guess.show();
    wouldbe_P.show();*/
    foreach.printf("out solve!\n"); fflush(stdout);
return ((ExOp::isValid(value))&&(value < 0.1f)) ? start_value : -100.0;}




LFHTEMP void SparseTrianglix<C>::show(FILE* f, int level)const{
    fprintf(f,"Sparse matrix with %i factored dimensions:", this->getSize());
    PartitionIterator pp(*this);
    int i,j;
    int nblone=0;
    if (pp.first()) do{
        if (pp.isCyclic()){
            if (part_data[part_id[pp.part[0]]][2] == ((part_data[part_id[pp.part[0]]][1] * (part_data[part_id[pp.part[0]]][1] - 1)) >>1)) fprintf(f,"\n Clique: {");
            else fprintf(f,"\n Cyclic(E=%i): {", part_data[part_id[pp.part[0]]][2] );
        }else{
            if (pp.getSize() != 1) fprintf(f,"\n Acyclic(E=%i): {", part_data[part_id[pp.part[0]]][2]);
            else if (part_id[pp.part[0]] == 0xFFFFFFFF) {nblone++; continue;}
            else fprintf(f,"\n Satellite {");
        }
        for(i=0;i<pp.getSize();i++) fprintf(f,"%i%c", pp.part[i], ((i+1) == pp.getSize()) ? '}' : ',');
        if (next_id[pp.part.last()] != pp.part[0]) fprintf(f,"\nWARNING! corrupted next or size for %i (s = %i?)\n", pp.part[0], pp.getSize());
        if (part_data[part_id[pp.part[0]]][3] != 0xFFFFFFFF){
            j = part_data[part_id[pp.part[0]]][3];
            if (getFirstID(j >> 16) != pp.part[0]) fprintf(f,"\nWarning! Illegal meta-edge %X for %i! (should be in %i %i)\n", j, pp.part[0], getFirstID(j >> 16), part_id[j >> 16] );
            fprintf(f," - {%i", getFirstID(j & 0xFFFF));
            j = attrib[j];
            while(j != 0xFFFFFFFF){
                if (getFirstID(j >> 16) != pp.part[0]) fprintf(f,"\nWarning! Illegal meta-edge %X for %i! (should be in %i %i)\n", j, pp.part[0], getFirstID(j >> 16), part_id[j >> 16] );
                fprintf(f,",%i", getFirstID(j & 0xFFFF));
                j = attrib[j];
            }
            fprintf(f,"}");
        }
    }while(pp.next());
    fprintf(f,"\n number of lone nodes: %i\n", nblone);
}

LFHTEMP void SparseTrianglix<C>::showShort(FILE* f, int level)const{
    fprintf(f,"Sparse matrix with %i factored dimensions:", this->getSize());
    PartitionIterator pp(*this);
    int i,j;
    int nblone=0;
    int nbsat=0;
    int nbcyclic=0;
    int nbacyclic=0;
    int nbclique=0;

    if (pp.first()) do{
        if (pp.isCyclic()){
            if (part_data[part_id[pp.part[0]]][2] == ((part_data[part_id[pp.part[0]]][1] * (part_data[part_id[pp.part[0]]][1] - 1)) >>1)) nbclique++;
            else nbcyclic++;
        }else{
            if (pp.getSize() != 1) nbacyclic++;
            else if (part_id[pp.part[0]] == 0xFFFFFFFF) nblone++;
            else nbsat++;
        }
    }while(pp.next());
    fprintf(f,"lone %i sat %i acyclic %i cyclic %i clique %i\n", nblone, nbsat, nbacyclic, nbcyclic, nbclique);
}

LFHTEMP void SparseTrianglix<C>::showWarn(FILE* f, int level)const{
    //fprintf(f,"Sparse matrix with %i factored dimensions:", this->getSize());
    PartitionIterator pp(*this);
    int i,j;
    int nblone=0;
    int nbcycle=0;
    int nbacycle=0;
    bool warn= false;
    bool has = false;
    if (pp.first()) do{
        if (pp.isCyclic()) nbcycle++;
        else{
            if (pp.getSize() != 1) nbacycle++;
            else if (part_id[pp.part[0]] == 0xFFFFFFFF) {nblone++; continue;}
            else nbacycle++;
        }
        /*if (pp.isCyclic()){
            if (part_data[part_id[pp.part[0]]][2] == ((part_data[part_id[pp.part[0]]][1] * (part_data[part_id[pp.part[0]]][1] - 1)) >>1)) fprintf(f,"\n Clique: {");
            else fprintf(f,"\n Cyclic(E=%i): {", part_data[part_id[pp.part[0]]][2] );
        }else{
            if (pp.getSize() != 1) fprintf(f,"\n Acyclic(E=%i): {", part_data[part_id[pp.part[0]]][2]);
            else if (part_id[pp.part[0]] == 0xFFFFFFFF) {nblone++; continue;}
            else fprintf(f,"\n Satellite {", pp.part[0]);
        }*/

        has |= (pp.part[0] == 2);
        for(i=1;i<pp.getSize();i++) {
            has |= (pp.part[i] == 2);
            if (isAFirstEntry(pp.part[i])){
                warn=true;
                printf("\nWARNING! corrupted partID, %i > %i thinks it is smallest\n", pp.part[i], pp.part[0]);
            }
        }

        if (next_id[pp.part.last()] != pp.part[0]) {warn=true;printf("\nWARNING! corrupted next or size for %i (s = %i?, s[%i] = %i)\n", pp.part[0], pp.getSize(), pp.getSize(), next_id[pp.part.last()]); }
        if (part_data[part_id[pp.part[0]]][3] != 0xFFFFFFFF){
            j = part_data[part_id[pp.part[0]]][3];
            if (getFirstID(j >> 16) != pp.part[0]) {warn=true;fprintf(f,"\nWarning! Illegal meta-edge %X for %i! (should be in %i %i)\n", j, pp.part[0], getFirstID(j >> 16), part_id[j >> 16] );}
            //fprintf(f," - {%i", getFirstID(j & 0xFFFF));
            j = attrib[j];
            while(j != 0xFFFFFFFF){
                if (getFirstID(j >> 16) != pp.part[0]) {warn=true; fprintf(f,"\nWarning! Illegal meta-edge %X for %i! (should be in %i %i)\n", j, pp.part[0], getFirstID(j >> 16), part_id[j >> 16] );}
                //fprintf(f,",%i", getFirstID(j & 0xFFFF));
                j = attrib[j];
            }
            //fprintf(f,"}");
        }
    }while(pp.next());
    //fprintf(f,"\n number of lone nodes: %i\n", nblone);

    for(i=0;i<attrib.getSize();i++){
        if ((attrib.deref(i) != 0xFFFFFFFF)&&(attrib.find(attrib.deref(i)) == 0xFFFFFFFF)) {warn=true; fprintf(f,"\nWarning! Illegal Meta-edge pointer entry %X->%X\n", attrib.deref_key(i), attrib.deref(i));}
        if (((attrib.deref_key(i) & 0xFFFF) >= neigh.getSize())||((attrib.deref_key(i) >> 16 ) >= neigh.getSize()))  {warn=true; fprintf(f,"\nWarning! Illegal Meta-edge key entry %X->%X\n", attrib.deref_key(i), attrib.deref(i));}

        if (part_id[getFirstID(attrib.deref_key(i) & 0xFFFF)] == 0xFFFFFFFF) {warn=true; printf("edge %X cannot be, %i->%i has no partition data...\n", attrib.deref_key(i), attrib.deref_key(i) & 0xFFFF, getFirstID(attrib.deref_key(i) & 0xFFFF));}
        if (part_id[getFirstID(attrib.deref_key(i) >> 16)] == 0xFFFFFFFF) {warn=true; printf("edge %X cannot be, %i->%i has no partition data...\n", attrib.deref_key(i), attrib.deref_key(i) >> 16, getFirstID(attrib.deref_key(i) >> 16));}
    }

    if (warn) printf("l%i c%i a%i\n", nblone, nbcycle, nbacycle);
    printf("*%i*\n", has ? 0xFFFF : 0xFFFFFF);
}

LFHTEMP bool SparseTrianglix<C>::PartitionIterator::first(){
    if (target.getSize() == 0) return false;
    part.setSize((target.part_id[0] == 0xFFFFFFFF) ? 1 : target.part_data[target.part_id[0]][1]);
    part[0] = 0;
    for(int i = 1; i < part.getSize();i++) part[i] = target.next_id[part[i-1]];
    return true;
}

LFHTEMP bool SparseTrianglix<C>::PartitionIterator::next(){
    uint32_t nextcur = part[0];
    for(nextcur++;nextcur< target.getSize();nextcur++) if (target.isAFirstEntry(nextcur)) break;
    //printf("found next %i\n", nextcur);
    if (nextcur >= target.getSize()) return false;
    part.setSize((target.part_id[nextcur] == 0xFFFFFFFF) ? 1 : target.part_data[target.part_id[nextcur]][1]);
    part[0] = nextcur;
    for(int i = 1; i < part.getSize();i++) part[i] = target.next_id[part[i-1]];
    return true;
}
LFHTEMP bool SparseTrianglix<C>::PartitionIterator::isCyclic(){
    if (part.getSize() < 3) return false;
    return (target.part_data[target.part_id[part[0]]][1] <= target.part_data[target.part_id[part[0]]][2] );
}

LFHTEMP void SparseTrianglixHanddle<C>::setSize(uint32_t _size){
    if (trianglix.getSize() == _size) return;;
    if (diago_inverse) {
        delete[](diago_inverse);
        diago_inverse = new C[_size];
    }
    if (part_determinant) {
        delete[](part_determinant);
        part_determinant = new C[_size];
    }
}
LFHTEMP const SparseTrianglix<C>& SparseTrianglixHanddle<C>::toZero(){
    // TODO...
}
LFHTEMP const SparseTrianglix<C>& SparseTrianglixHanddle<C>::toOne(){
    int i;
    if (part_determinant) for(i=0;i< trianglix.getSize();i++){
        ExOp::toOne(part_determinant[i]);
    }
}

LFHTEMP SparseTrianglixHanddle<C>& SparseTrianglixHanddle<C>::toMemmove(SparseTrianglixHanddle<C>& other){
    trianglix.toMemmove(other.trianglix);
    diago_inverse = other.diago_inverse; other.diago_inverse = NULL;
    part_determinant = other.part_determinant; other.part_determinant = NULL;
    return *this;
}


LFHTEMP void SparseTrianglixHanddle<C>::setMaint_InvDiago(bool _does_maintain){
    if (!_does_maintain) {delete[](diago_inverse); diago_inverse = NULL; return;}
    if (diago_inverse) return;
    if (trianglix.getSize() == 0){
        diago_inverse = new C[1];
        return;
    }
    // TODO!
}
LFHTEMP void SparseTrianglixHanddle<C>::setMaint_Determinant(bool _does_maintain){
    if (!_does_maintain) {delete[](part_determinant); part_determinant = NULL; return;}
    if (part_determinant) return;
    if (trianglix.getSize() == 0){
        part_determinant = new C[1];
        return;
    }
    // TODO!
}

/*
#undef LFHTEMP
#define LFHTEMP template <class TARG,bool READONLY>
LFHTEMP void LockPtr<TARG,READONLY>::initAlias(const unsigned int &alias){
   unsigned int tridmask = ((unsigned int)Controlstate::ThreadID_mask[SDL_ThreadID()]) << 28;
   unsigned int ite =  AliasBank.find(alias);// printf("found %i\n", ite);
    if (ite == 0xFFFFFFFF) target =NULL;
    else{
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        if ((AliasBank.deref(ite).second & 0xF0000000) == 0xF0000000) { // read-only lock!
            target = NULL;
        }else if ((AliasBank.deref(ite).second & 0xFF000000) != 0) {
            if  ((AliasBank.deref(ite).second & 0xF0000000) == tridmask) {AliasBank.deref(ite).second += 0x01000000; target = (TARG*) AliasBank.deref(ite).first;}
            else target = NULL;
        }else{
            AliasBank.deref(ite).second += 0x01000000;
            if ((AliasBank.deref(ite).second & 0xFF000000) == 0x01000000) {target = (TARG*) AliasBank.deref(ite).first; AliasBank.deref(ite).second |= tridmask;}
            else{ AliasBank.deref(ite).second -= 0x01000000;target = NULL;}
        }
        }
}

 LFHTEMP LockPtr<TARG,READONLY>::LockPtr(const AliasPtr<TARG> &tar){initAlias(tar.alias);}
 LFHTEMP LockPtr<TARG,READONLY>::LockPtr(const unsigned int &alias){initAlias(alias);}


 LFHTEMP LockPtr<TARG,READONLY>::~LockPtr(){if (target != NULL) {
        unsigned int ite = AliasOf.find(target);
        ite = AliasBank.find(AliasOf.deref(ite));
        AliasBank.deref(ite).second -= 0x01000000;
        if ((AliasBank.deref(ite).second & 0x0F000000) == 0) AliasBank.deref(ite).second &= 0x0FFFFFFF;
    }
    }
 LFHTEMP TARG* LockPtr<TARG,READONLY>::operator->()const{return target;}
 LFHTEMP TARG& LockPtr<TARG,READONLY>::operator*()const{return *target;}

 LFHTEMP bool LockPtr<TARG,READONLY>::isValid()const{return target != NULL;}

 LFHTEMP template<class NTYPE> LockPtr<TARG,READONLY>::operator NTYPE* ()const{return (NTYPE*) target;}
#undef LFHTEMP
#define LFHTEMP template <class TARG>
LFHTEMP void LockPtr<TARG,true>::initAlias(const unsigned int &alias){
   unsigned int tridmask = ((unsigned int)Controlstate::ThreadID_mask[SDL_ThreadID()]) << 28;
   unsigned int ite =  AliasBank.find(alias);
   printf("ThreadID %i\n",SDL_ThreadID());
    if (ite == 0xFFFFFFFF) target =NULL;
    else{
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        AliasBank.deref(ite).second += 0x01000000;
       // printf("mass found %i\n", (AliasBank.deref(ite).second & 0xFF000000) >> 24);
        if ((AliasBank.deref(ite).second & 0xFF000000) == 0x01000000) target = (TARG*) AliasBank.deref(ite).first;
        else{AliasBank.deref(ite).second -= 0x01000000;
        target = NULL;}
    }
    }

 LFHTEMP LockPtr<TARG, true>::LockPtr(const AliasPtr<TARG> &tar){initAlias(tar.alias);}
 LFHTEMP LockPtr<TARG, true>::LockPtr(const unsigned int &alias){initAlias(alias);}
 LFHTEMP LockPtr<TARG, true>::~LockPtr(){if (target != NULL) {
        unsigned int ite = AliasOf.find(target);
        ite = AliasBank.find(AliasOf.deref(ite));
        AliasBank.deref(ite).second -= 0x01000000;}
    }
 LFHTEMP const TARG* LockPtr<TARG, true>::operator->()const{return target;}
 LFHTEMP const TARG& LockPtr<TARG, true>::operator*()const{return *target;}
 LFHTEMP bool LockPtr<TARG, true>::isValid()const{return target != NULL;}

 LFHTEMP template<class NTYPE> LockPtr<TARG, true>::operator const NTYPE* ()const{return (const NTYPE*) target;}


template<class NTYPE> LockPtr<void, true>::operator const NTYPE* ()const{return (const NTYPE*) target;}
template<class NTYPE> LockPtr<void, false>::operator NTYPE* ()const{return (NTYPE*) target;}
*/
/*
 LFHTEMP template<class HOST> void AliasPtr<TARG>::setDestroyOnClear(HOST* what){
    if (alias == 0) {fprintf(stderr,"alias based AutoDelete on void target, is illegal!\n");exit(1);}
    unsigned int ite =  static_alias_bank.auto_destroy.find(alias);
    if (ite == 0xFFFFFFFF) {
        static_alias_bank.auto_destroy[alias] = Vector< ListenerPointer<int> >();
        ite =  static_alias_bank.auto_destroy.find(alias);
        }
//    static_alias_bank.auto_destroy.deref(ite).push_back(MakeitDie<HOST>(what));
}*/

#undef LFHTEMP
#define LFHTEMP template<class K, class D, class C, class H>


/*
LFHTEMP void CategoryHashmap<K,D,C,H>::insert(const K&, const D&, const C&){

    }

LFHTEMP void CategoryHashmap<K,D,C,H>::remove(const K&){

    }*/

/*template <class C> angle<C>::angle(C const & value) : ang(value){}
template <class C>
angle<C>::operator double() const{
	return( (double)ang);
}
template <class C>
angle<C>::operator complex() const{
	double a = (double)(*this);
	return(complex(cos(a),sin(a)));
}
template <class C>
complex angle<C>::getcomplex() const{
	double a = (double)(*this);
	return(complex(cos(a),sin(a)));
}*/
