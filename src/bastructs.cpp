/*
 * bastructs.cpp
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

#include "bastructs.h"

namespace LFHPrimitive{

LFHPrimitive::ThreadBase foreach;

char* cloneString(const char* const what){
    unsigned int i = strlen(what) +1;
    char* _out = new char[i];
    memcpy(_out,what,sizeof(char)*i);
return(_out);}

ProgressBarPrint::ProgressBarPrint(uint32_t _length):length(_length),state(0){}
void ProgressBarPrint::start(const char* text){
    unsigned int i;
    for(i=0;i<length;i++) putc(' ', stdout);
    fprintf(stdout, ": %s", text);
    for(i=length + strlen(text)+2;i != 0;i--) putc('\b', stdout);fflush(stdout);
    lasttime = clock();
    state =0;
}
void ProgressBarPrint::update(double fraction){
    if (((unsigned int)(fraction * length * 16)) > state) {
    if (state & 15) putc('\b', stdout);
    state &= 0xFFFFFFF0;
        while (((unsigned int)(fraction * length * 16)) > (state | 15)) {state += 16; putc('#', stdout);}
        state |= ((unsigned int)(fraction * length * 16)) & 15;
        switch(state & 15){
            case 1: putc('.', stdout);break;
            case 2: putc(',', stdout);break;
            case 3: putc(':', stdout);break;
            case 4: putc(';', stdout);break;
            case 5: putc('i', stdout);break;
            case 6: putc('j', stdout);break;
            case 7: putc('l', stdout);break;
            case 8: putc('!', stdout);break;
            case 9: putc('|', stdout);break;
            case 10: putc('L', stdout);break;
            case 11: putc('C', stdout);break;
            case 12: putc('G', stdout);break;
            case 13: putc('O', stdout);break;
            case 14: putc('Q', stdout);break;
            case 15: putc('@', stdout);break;
        }
        if ((clock() - lasttime) & 0xFFFFFC00 ) {fflush(stdout); lasttime =clock();}
    }
}
void ProgressBarPrint::finish(){
    if (state & 15) putc('\b', stdout);
    while (length > (state >> 4)) {state+= 16;putc('#', stdout);}
    fprintf(stdout,"\n");fflush(stdout);
}

ERRCODE Workflow::operator()(){
	tb.startThreadArray();
	uint32_t orderite =0;
	uint32_t dirty;
    tb.stopThreadArray();
return 0;}
void Workflow::show(FILE *f, int level) const{
	for(int i =0 ; i<task_order.getSize();i++){
		if (task_order[i][1] == 0) fprintf(f, "Wait for %s\n", tb.tasknodes[task_order[i][0]].taskname.c_str());
		else fprintf(f, "Submit task %i for %s\n", i, tb.tasknodes[task_order[i][0]].taskname.c_str());
	}
}

Workflow ThreadBase::makeWorkflow(uint32_t target){ Workflow wfout(*this);

return(wfout);}

ThreadBase::ThreadBase():nbactive(0), progb(20),async_progress_maintain(0xFFFFFFFF),nbthreads(0),running(true){
}

ThreadBase::~ThreadBase(){this->toMemfree(); stopThreadArray(); joinThreads();}
ThreadBase& ThreadBase::toMemfree(){while(nbactive) {thrds[--nbactive]->join();delete(thrds[nbactive]);}
    if (nbthreads) {delete[](thrds); nbthreads = 0;}
return *this;}

ERRCODE ThreadBase::startEvent(Event<void>* ev){
    if (nbactive == nbthreads) return 1;
    //thrds[nbactive++] = new std::thread(Event<void>::operator(), ev);
    thrds[nbactive] = new std::thread(ThreadBase::callThatEventVoid, ev); LFH_NICE_ALLOCERROR(thrds[nbactive],"")
    nbactive++;
    return 0;
}
void ThreadBase::waitForAllThreads(){
    if (async_progress_maintain != 0xFFFFFFFF){
        while(nbactive>0){async_progress_maintain = nbactive--; thrds[nbactive]->join(); delete(thrds[nbactive]);}
        progb.finish();
        async_progress_maintain = 0xFFFFFFFF;
    }else while(nbactive>0){async_progress_maintain = nbactive--; thrds[nbactive]->join(); delete(thrds[nbactive]);}
}
uint32_t ThreadBase::callThatEventVoid(Event<void>* ev) {return (*ev)();}
void ThreadBase::submit(std::function<void()> what){
    {
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(what);
        #endif // STDBINDFILTER
    }
    condvar.notify_one();
}
void ThreadBase::submit_ThenWait(std::function<void()> what){
    /*{
        std::lock_guard<std::mutex> lock(mut);
        #ifdef STDBINDFILTER
        todolist.push_front(what);
        #endif // STDBINDFILTERwaitForAllThreads
    }

    std::unique_lock<std::mutex> innerlock(mut);
    while(true){
        condvar.wait(innerlock,[this]{return true;});
        if (todolist.empty()) break;
        std::function<void()> fout = std::move(todolist.back());
        todolist.pop_back();
        nbactivethr++;
        innerlock.unlock();
    }
    innerlock.unlock();*/
    what();
    {
        std::unique_lock<std::mutex> lock(mut);
        main_condvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
    }
    if (async_progress_maintain != 0xFFFFFFFF){
        progb.finish();
        async_progress_maintain = 0xFFFFFFFF;
    }
    flushMsgs();
}


int ThreadBase::printLog(const char *__format, ...) {
    std::unique_lock<std::mutex> lck(accessFinalMutex());
    va_list __va; va_start(__va, __format);
    // std::fprintf(log, "thr%i:", std::this_thread::get_id());
    int __res = std::vfprintf(log, __format, __va);
    va_end(__va);
    fflush(log);
return __res;}
int ThreadBase::printLogF(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    std::fprintf(log, "thr%X:", (uint32_t)std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(log, __format, __va);
    va_end(__va);
    fflush(log);
return __res;}
int ThreadBase::printf(const char *__format, ...) {
    std::unique_lock<std::mutex> lck(accessFinalMutex());
    va_list __va; va_start(__va, __format);
    std::fprintf(stdout,"thr%X:", (uint32_t)std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout, __format, __va);
    va_end(__va);
return __res;}
int ThreadBase::printf_l(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    std::fprintf(stdout,"thr%X:", (uint32_t)std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout,__format, __va);
    va_end(__va);
return __res;}
int ThreadBase::fprintf(FILE* f, const char *__format, ...) {
    if (f == NULL) f = log;
    std::unique_lock<std::mutex> lck(accessFinalMutex());
    va_list __va; va_start(__va, __format);
    std::fprintf(f, "thr%X:", (uint32_t)std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(stdout, __format, __va);
    va_end(__va);
return __res;}
int ThreadBase::fprintf_l(FILE* f,const char *__format, ...) {
    if (f == NULL) f = log;
    va_list __va; va_start(__va, __format);
    std::fprintf(f, "thr%X:", (uint32_t)std::hash<std::thread::id>()(std::this_thread::get_id()) & 65535);
    int __res = std::vfprintf(f, __format, __va);
    va_end(__va);
    fflush(f);
return __res;}

int ThreadBase::print(const char *__format, ...) {
    va_list __va; va_start(__va, __format);
    char buffer[65536];
    int __res = std::vsprintf(buffer, __format, __va);
    va_end(__va);
    msgs.insert(string(buffer));
return __res;}


void ThreadBase::terminate(const char* __format, ...){
    running = false;
    va_list __va; va_start(__va, __format);
    char buffer[65536];
    int __res = std::vsprintf(buffer, __format, __va);
    va_end(__va);
    msgs.insert(string(buffer));
exit(1);}


std::function<void()> ThreadBase::getFunc(){
    std::unique_lock<std::mutex> lock(mut);
    condvar.wait(lock,[this]{return !todolist.empty(); });
    std::function<void()> fout = std::move(todolist.back());
    todolist.pop_back();
return(fout);}

void ThreadBase::operator()(){
    while(true){
        std::unique_lock<std::mutex> lock(mut); // ,std::defer_lock
        condvar.wait(lock,[this]{return !todolist.empty()||(!running); });
        if (!running) break;
        std::function<void()> fout = std::move(todolist.back());
        nbactivethr++;
        todolist.pop_back();
        //printf("^^^ %i\n", nbactivethr);
        //printf("Thread Awoke! lets do this %i\n", nbactivethr);
        lock.unlock();
        fout();
//        try{ fout();} catch(...) {std::throw_with_nested( std::runtime_error("Inside Task %p") );}
        lock.lock();
        nbactivethr--;
		//printf("vvv %i\n", nbactivethr);

		if (nbactivethr == 0){
			main_condvar.notify_one();
		//	printf("notification done! %i %i\n", nbactivethr, todolist.size());
        }
    }
    condvar.notify_one();
}
void ThreadBase::joinAll(){
    std::unique_lock<std::mutex> lock(mut);
    main_condvar.wait(lock,[this]{return ((nbactivethr==0)&&(todolist.empty())) ||(!running); });
    /*while(true){
        mut.lock();
        if (!running) {mut.unlock(); break;}
        if (todolist.empty()){
            mut.unlock();
            {
                std::unique_lock<std::mutex> lock(endmut);
                endcondvar.wait(lock,[this]{return (nbactivethr==0)||(!running)||(!todolist.empty()); });
                lock.unlock();
            }
            if ((!running)||((nbactivethr==0)&&(todolist.empty()))) break;
        }else{
            std::function<void()> fout = std::move(todolist.back());
            todolist.pop_back();
            mut.unlock();
            fout();
        }
    }*/
}
uint32_t ThreadBase::startThread(uint32_t thr_input, std::function< int(uint32_t) > fnc){
    uint32_t r;
    do{
        ExOp::toRand(r);
    }while((r == 0)||(dedicated.find(r) != 0xFFFFFFFF));
    dedicated[r] = new std::thread(std::bind(fnc, thr_input));
return r;}
void ThreadBase::joinThread(uint32_t thrID){
    if (dedicated.find(thrID) == 0xFFFFFFFF){
        this->printf("warning! trying to join an un-existing thread\n");
    }else{
        dedicated[thrID]->join(); delete(dedicated[thrID]); dedicated.erase(thrID);
    }
}
void ThreadBase::startThreadArray(uint32_t _nbthreads){
	if (nbthreads !=0) {printf("ThreadArray is already started...\n"); exit(1);}
	printf("Init threadArray for length %i\n", _nbthreads);
    running = true;
    nbactivethr =0;
    nbthreads = _nbthreads;
    mainthreadid = std::this_thread::get_id();
   // #ifdef STDBINDFILTER
    for(int i=0;i<_nbthreads;i++) futures.emplace_back(std::bind(&ThreadBase::operator(), std::ref(*this)));
  //  #endif // STDBINDFILTER
}
void ThreadBase::stopThreadArray(){
	printf("Stopping threadArray for length %i\n", nbthreads);
    running = false;
    condvar.notify_all();
    fflush(stdout);
    for(int i=0;i<futures.size();i++) futures[i].join();
    futures.clear();
    flushMsgs();
    nbthreads =0;
	printf("Done\n");
}
void ThreadBase::joinThreads(){
    if (auto ite = dedicated.mkIterator()) do{
        (*ite)->join();
        delete(*ite);
    }while(ite++);
    dedicated.toMemfree();
    flushMsgs();
}
void ThreadBase::flushMsgs(FILE *f){
    int cnt =0;
    string msgstar,old;
    if (msgs.pop(old)){
        while (msgs.pop(msgstar)){
        if (old == msgstar) cnt++;
        else {
        if (cnt) {std::fprintf(f,"(x%i) %s\n", cnt+1, old.c_str()); cnt=0;}
        else std::fprintf(f,"%s\n", old.c_str());
        old = std::move(msgstar);
        }
        }
        if (cnt) std::fprintf(f,"(x%i) %s\n", cnt+1, old.c_str());
        else std::fprintf(f,"%s\n", old.c_str());
    }
}
void ThreadBase::startEvent_ThenWait(Event<void>* ev){(*ev)();this->waitForAllThreads();}
void ThreadBase::startProgress(const char* text, uint32_t totalsteps){
    async_progress_maintain = 0u;
    async_progress = 0;
    async_progress_max = totalsteps;
    progb.start(text);
    Vector<double> hehe;
}
void ThreadBase::updateProgress(uint32_t threadID){
    if (threadID != async_progress_maintain) async_progress++;
    else progb.update(((double)async_progress++) / async_progress_max);
}
void ThreadBase::finishProgress(uint32_t threadID){
    if (threadID == async_progress_maintain) {async_progress_maintain = 0xFFFFFFFF; progb.finish();}
}

void ThreadBase::initEqualRanges(Tuple<uint32_t> &fout, uint32_t nbelems, uint32_t nb_threads) const{
    if (nb_threads == 0) nb_threads = nbthreads;
    fout.setSize(nb_threads+1);
    for(uint32_t i =0 ; i <= nb_threads;i++) fout[i] = (i *  nbelems) / nb_threads;
}

bool ThreadBase::ThreadArrayScope::getRangeVal(uint32_t& val) {
    bool fout = (val = range_iterator.fetch_add(1)) < range_limit;
    std::unique_lock<std::mutex> lck(tb.mut, std::defer_lock_t());
    if (fout) {
        tb.async_progress++;
        if (lck.try_lock()) tb.progb.update(((double)tb.async_progress) / tb.async_progress_max);
    }else if (val == range_limit){
        lck.lock();
        tb.progb.finish();
    }
return fout;}

HiddenClasses::DictionaryBase::DictionaryBase(): flags(0){}
HiddenClasses::DictionaryBase::~DictionaryBase(){this->toMemfree();}
HiddenClasses::DictionaryBase::DictionaryBase(const HiddenClasses::DictionaryBase& o): flags(o.flags), word_links(o.word_links), hashes(o.hashes){
    entries.setSize(o.entries.getSize());
    for(uint32_t i=0;i< entries.getSize();i++) {uint32_t j = strlen(o.entries[i])+1; entries[i] = new char[j]; memcpy(entries[i], o.entries[i],j);}
}
HiddenClasses::DictionaryBase& HiddenClasses::DictionaryBase::operator=(const HiddenClasses::DictionaryBase&o){flags = o.flags; word_links = word_links; hashes = o.hashes;
    entries.setSize(o.entries.getSize());
    for(uint32_t i=0;i< entries.getSize();i++) {uint32_t j = strlen(o.entries[i])+1; entries[i] = new char[j]; memcpy(entries[i], o.entries[i],j);}
return *this;}
HiddenClasses::DictionaryBase& HiddenClasses::DictionaryBase::toMemfree(){
    for(uint32_t i=0;i< entries.getSize();i++) delete[](entries[i]);
    entries.toMemfree();
    word_links.toMemfree();
    hashes.toMemfree();
return *this;}
uint32_t HiddenClasses::DictionaryBase::findWord(const char* word, int word_lenght) const{
	//printf("find lenght: %i\n", word_lenght);

	int proc = 0;
	uint32_t key;
	uint32_t selhash = 0;
	uint32_t ite;
    if (hashes.getSize() == 0) return 0xFFFFFFFF;
	while(word_lenght - proc >= 4){
		key = (((int)word[proc]) << 24) | (((int)word[proc | 1]) << 16) | (((int)word[proc | 2]) << 8) | ((int)word[proc | 3]);
		ite = hashes[selhash].find(key);
		if (ite == 0xFFFFFFFF) return 0xFFFFFFFF;
		selhash = hashes[selhash].deref(ite);
		//printf("%i %i tmp\n", selhash, proc);
		//hashes[selhash].show();
		proc+= 4;
	}
	if (word_lenght - proc >= 2){
		key = (((int)word[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)word[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)word[proc]) << 24);
	ite = hashes[selhash].find(key);
	//if (ite == 0xFFFFFFFF) printf("cant find %s\n", word);
	//else printf("found %s at %X\n", word,hashes[selhash].deref(ite));
return (ite == 0xFFFFFFFF) ? 0xFFFFFFFF : hashes[selhash].deref(ite);}
uint32_t HiddenClasses::DictionaryBase::insertWord(const char* word, int word_lenght){
	//printf("inserting word %s\n", word); fflush(stdout);
	uint32_t fout = word_links.getSize();
	word_links.push_back();
//	word_links[fout].k = new char[word_lenght+1];
//	memcpy(word_links[fout].k,word, word_lenght);
//	word_links[fout].k[word_lenght] = '\0';
	int proc = 0;
	unsigned int ite;
	uint32_t selhash = 0;
	if (hashes.getSize() ==0 ) hashes.setSize(1);
	uint32_t key;
	while(word_lenght - proc >= 4){
		key = (((int)word[proc]) << 24) | (((int)word[proc | 1]) << 16) | (((int)word[proc | 2]) << 8) | ((int)word[proc | 3]);
		if ((ite = hashes[selhash].find(key)) != 0xFFFFFFFF) {
            selhash =  hashes[selhash].deref(ite);
          //  printf("proc%i got to %i\n", proc, selhash);
		}else{
		   // printf("creating new hash for %X\n", key);
			hashes[selhash][key] = hashes.getSize();
			hashes[selhash].keysort();
			selhash = hashes.getSize();
			hashes.push_back();
		}
		proc+= 4;
	}

	if (word_lenght - proc >= 2){
		key = (((int)word[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)word[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)word[proc]) << 24);
	ite = hashes[selhash].find(key);
	if (ite != 0xFFFFFFFF) {fprintf(stderr, "DICO error: trying to insert word that already exist!\n"); LFH_exit(1);}
//	printf("final entry within %X, (index %i)", key, selhash);
	hashes[selhash][key] = fout;
	hashes[selhash].keysort();
return fout;}
Tuple<uint32_t, 4u> HiddenClasses::DictionaryBase::findPrefix(const char* prefix, int word_lenght) const{
	Tuple<uint32_t, 4u> fout;
	int proc = 0;
	uint32_t key;
	uint32_t selhash = 0;
	uint32_t ite;
	if (hashes.getSize() ==0 ) {fout[0]= 0xFFFFFFFF; return fout;}

	while(word_lenght - proc >= 4){
		key = (((int)prefix[proc]) << 24) | (((int)prefix[proc | 1]) << 16) | (((int)prefix[proc | 2]) << 8) | ((int)prefix[proc | 3]);
		ite = hashes[selhash].find(key);
		if (ite == 0xFFFFFFFF) {fout[0]= 0xFFFFFFFF; return fout;}
		selhash = hashes[selhash].deref(ite);
		proc+= 4;
	}
	if (word_lenght - proc >= 2){
		key = (((int)prefix[proc | 1]) << 16) | ((word_lenght - proc == 3) ?  (((int)prefix[proc | 2]) << 8) : 0);
	} else key = 0;
	if (word_lenght - proc > 0) key |= (((int)prefix[proc]) << 24);

	fout[0] = 0;
	fout[1] = hashes[selhash].getSize()-1;
	//printf("binary search!\n");

	while(fout[0] != fout[1]){
		if (hashes[selhash].deref_key((fout[0]+fout[1]) >> 1) >= key) fout[1] = (fout[0]+fout[1]) >> 1;
		else fout[0] = ((fout[0]+fout[1]) >> 1) +1;
	}


	fout[2] = key;
	fout[3] = (word_lenght == proc) ? 0 : 0xFFFFFFFF << ((4 - word_lenght + proc) * 8);

	//printf("%X and %X xor : %X & %X\n", key, hashes[selhash].deref_key(fout[1]), (hashes[selhash].deref_key(fout[1]) ^ key), fout[3]);

	ite = hashes[selhash].deref_key(fout[1]) ^ key;



	fout[0] = (ite & fout[3]) ? 0xFFFFFFFF : selhash;

	//if (ite == 0xFFFFFFFF) printf("cant find %s\n", word);
	//else printf("found %s at %X\n", word,hashes[selhash].deref(ite));
return fout;}
ERRCODE HiddenClasses::DictionaryBase::save(FILE*f) const{
    uint32_t i = entries.getSize();
    if (fwrite(&flags,sizeof(uint32_t),1,f) != 1) return 1;
    if (fwrite(&i,sizeof(uint32_t),1,f) != 1) return 1;
    for(i=0;i<entries.getSize();i++){
        if (fwrite(entries[i],sizeof(char),strlen(entries[i])+1,f) != strlen(entries[i])+1) return 1;
    }
    word_links.save(f);
return hashes.save(f);}
ERRCODE HiddenClasses::DictionaryBase::load(FILE*f){
    char buffer[256];
    char* cur;
    uint32_t i;
    if (fread(&flags,sizeof(uint32_t),1,f) != 1) return 1;
    if (fread(&i,sizeof(uint32_t),1,f) != 1) return 1;
    entries.setSize(i);
    for(i=0;i<entries.getSize();i++){
        cur = buffer;
        do{
            if (fread(cur,sizeof(char),1,f) != 1) return 1;
        }while (*(cur++) != '\0');
        entries[i] = new char[(int)(cur - buffer)];
        strcpy(entries[i], buffer);
    }
    word_links.load(f);
    hashes.load(f);
return 0;}
void HiddenClasses::DictionaryBase::addEntry_routine(const char* str, uint32_t strlenght){ // , bool allow_duplicate
	uint32_t i;
	char* newstr = new char[strlenght+1];
	memcpy(newstr, str, strlenght+1);
	entries.push_back(newstr);
	int c;
	if (flags & 1){ // fix Capitalize
		for(c=0;newstr[c] != '\0';c++){
			if (((newstr[c] & 0xC0) == 0x40)&&(newstr[c] <= 'Z')) newstr[c] |= 0x20;
		}
	}
	const char* c_start = newstr;
	const char* c_end = newstr;
	uint32_t word_index;
	while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			word_index = this->findWord(c_start, (int)(c_end - c_start));
			if (word_index == 0xFFFFFFFF) word_index = insertWord(c_start, (int)(c_end - c_start));
			for(i=0;i<word_links[word_index].getSize();i++) if (word_links[word_index][i] == entries.getSize()-1) break;
			if (i==word_links[word_index].getSize()) word_links[word_index].push_back(entries.getSize()-1);
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}
	if (flags & 1) memcpy(newstr, str, strlenght+1); // fix Capitalize
}

void HiddenClasses::DictionaryBase::findIndexes(Vector<uint32_t> &fout, const char* query, uint32_t max_requested) const{
	uint32_t i = strlen(query);
	char* mquery;
	const char* c_start;
	const char* c_end;
	int c;
	if (flags & 1){
		mquery = new char[i+1];
		memcpy(mquery, query, i+1);
		for(c=0;mquery[c] != '\0';c++){
			if (((mquery[c] & 0xC0) == 0x40)&&(mquery[c] <= 'Z')) mquery[c] |= 0x20;
		}
		c_start = mquery;
		c_end = mquery;
	}else {
		c_start = query;
		c_end = query;
	}

	Vector< Tuple<uint32_t, 4u> > found_indexes;


	while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			found_indexes.push_back(this->findPrefix(c_start, (int)(c_end - c_start)));
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}
	uint32_t j;
	int k,l;
	Vector< KeyElem<uint32_t, uint32_t> > scope;

	if (found_indexes.getSize() == 1){
		//printf("found indexes: %X %X %X %X\n", found_indexes[0][0], found_indexes[0][1], found_indexes[0][2],found_indexes[0][3]);
		if (found_indexes[0][0] != 0xFFFFFFFF)
		for(i=found_indexes[0][1];i<hashes[found_indexes[0][0]].getSize();i++){
			if ((hashes[found_indexes[0][0]].deref_key(i) ^ found_indexes[0][2]) & found_indexes[0][3]) break; // mismatch
			k = hashes[found_indexes[0][0]].deref(i);
			if (hashes[found_indexes[0][0]].deref_key(i) & 0xFF){ // sub is hash!
				scope.push_back( KeyElem<uint32_t, uint32_t>(k,0) );
				while( (l = scope.size()) > 0){
					l--;
					if (hashes[scope[l].k].deref_key(scope[l].d) & 0xFF){
						scope.push_back( KeyElem<uint32_t, uint32_t>(hashes[scope[l].k].deref(scope[l].d),0) );
					}else{
						k = hashes[scope[l].k].deref(scope[l].d);
						for(j=0 ;j <word_links[k].getSize();j++) {
							fout.push_back(word_links[k][j]);
							if (fout.getSize() == max_requested) {
								fout.sort_unique();
								if (fout.getSize() == max_requested){
									if (flags & 1) delete[](mquery);
									return;
								}
							}
						}
						while(true){
							scope[l].d++;
							if ( scope[l].d < hashes[scope[l].k].getSize()) break;
							else{
								scope.pop_back();
								if (l ==0) break;
								l--;
							}
						}
					}
				}
			}else{

				for(j=0 ;j <word_links[k].getSize();j++) {

					fout.push_back(word_links[k][j]);
					if (fout.getSize() == max_requested) {
						fout.sort_unique();
						if (fout.getSize() == max_requested){
							if (flags & 1) delete[](mquery);
							return;
						}
					}
				}
			}
		}
	}else{
		/*
		int *subind = new int[word_indexes.getSize()];
		uint32_t curmin;
		uint32_t curmin_index;

		for(i=0;i<word_indexes.getSize();i++) subind [i] =0;*/
	}

	if (fout.getSize() != 0) fout.sort_unique();
	if (flags & 1) delete[](mquery);
}
void HiddenClasses::DictionaryBase::findKeys(Vector<const char*> &fout, const char* query, uint32_t max_requested)const{
	Vector< uint32_t> foutindex;
	this->findIndexes(foutindex,query,max_requested);
	for(uint32_t i=0;i<foutindex.getSize();i++) fout.push_back(entries[foutindex[i]]);
}

uint32_t HiddenClasses::DictionaryBase::findEntry(const char* query) const{
    // if query is a word, ez-pez;
    const char* c_start = query;
	const char* c_end = query;

    Vector<uint32_t> wordind;
    while(true){
		while(true){
			if ((*c_end == '\0')||(*c_end == ' ')||(*c_end == '\t')) break;
			c_end++;
		}
		if (c_end != c_start){
			wordind.push_back(this->findWord(c_start,(int)(c_end - c_start)));
			if (wordind.last() == 0xFFFFFFFF) return 0xFFFFFFFF;
		}
		if (*c_end == '\0') break;
		c_end++;
		c_start = c_end;
	}
    uint32_t i,j;
    for(i=0,j=1;j<wordind.getSize();j++){
        if (word_links[wordind[j]].getSize() < word_links[wordind[i]].getSize()) i=j;
    }
    j = wordind[i];
    for(i=0;i< word_links[j].getSize();i++) if (strcmp(query, entries[word_links[j][i]]) == 0) break;
return (i< word_links[j].getSize()) ? word_links[j][i] : 0xFFFFFFFF;}

Vector<uint32_t> HiddenClasses::DictionaryBase::findEntries(const vector<string> &queries)const{Vector<uint32_t> fout;
    fout.setSize(queries.size());
    uint32_t i,j;
    for(i=0,j=0;i< fout.getSize();i++)  if ((fout[i] = this->findEntry(queries[i+j].c_str())) == 0xFFFFFFFF) {printf("Warning, did not find %s\n", queries[i+j].c_str()); i--; j++; fout.pop_back();}
return fout;}

//uint32_t Dictionary::find(const char* what) const{
//    return findWord(what, strlen(what));
    /*uint32_t key;
    uint32_t selhash = 0;
    if (strlen(what) <3){
        if (strlen(what) == 2) ((((uint32_t)what[0]) << 8) | what[1]) << 16;
        else key = (((uint32_t)what[0]) << 24);
    }else key = (((((((uint32_t)what[0]) << 8) | what[1]) << 8) | what[2])) | what[3];
    ite = hashes[selhash].find(key);
    while(ite != 0xFFFFFFFF){
        selhash =  hashes[selhash].deref(ite);

    }
    return 0xFFFFFFFF;*/
//}

//void Dictionary::registerSubstrings(int maxlen){
//	unsigned int i;
/*,j,cur;
	maxlen = (maxlen + 3) >>2;

	unsigned int ite;
	unsigned int curlen;

	char cap[6];
	cap[2] ='\0';
	uint32_t key;
	uint32_t nbsub;
	hashes.setSize(1);

	for(i=0;i< entries.size();i++){
		key =0;
		curlen = strlen(entries[i]);
		if (curlen > 0){
			cap[1] = entries[i][curlen-1];
			cap[3] = entries[i][0];
			if (curlen > 1){
				cap[0] = entries[i][curlen-2];
				cap[4] = entries[i][1];
				cap[5] = (curlen > 2) ? entries[i][2] : '\0';
			}else{
				cap[4] = '\0';
				cap[5] = '\0';
			}
		}else{
			cap[3] = '\0';
			cap[4] = '\0';
			cap[5] = '\0';
		}

		for(cur=0; entries[i][cur] != '\0';cur++){
			key =  *(uint32_t*)(curlen - cur > 2 ? entries[i] + cur : cap + 3 + cur - curlen);
			ite =hashes[0].find(key);
			if (ite == 0xFFFFFFFF){
				nbsub = hashes.getSize();
				hashes[0][key] = nbsub | 0x80000000;
				key = *(uint32_t*)(curlen - cur > 6 ? entries[i] + cur+4 : curlen > cur+4 ? cap + 7 + cur - curlen: cap+2);
				hashes.push_back();
				hashes[nbsub][key] = i;
			}else{
				nbsub = hashes[0].deref(ite);
				key = *(uint32_t*)(curlen - cur > 6 ? entries[i] + cur+4 : curlen > cur+4 ? cap + 7 + cur - curlen: cap+2);
				ite = hashes[nbsub].find(key);
				if (ite ==  0xFFFFFFFF) hashes[nbsub][key] = i;
				// else TODO!
			}
		}
	}*/

	//Compare_key ck;
//	for(i=0;i<hashes.size();i++) hashes[i].keysort();
//}







} // end of namespace

