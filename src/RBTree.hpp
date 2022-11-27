/*

 Copyright (C) 2013 Louis-Francois Handfield
 All rights reserved.

 */




#undef LFHTEMP
#define LFHTEMP template <class Key, class Comparator>

LFHTEMP RBTree<Key,Comparator>::RBTree(): root(NULL),first(NULL),last(NULL){}
LFHTEMP RBTree<Key,Comparator>::~RBTree(){flush();}
LFHTEMP void RBTree<Key,Comparator>::flush(){
RBTreeNode<Key,Comparator>* cur = first;
RBTreeNode<Key,Comparator>* tmp;
while(cur != NULL){
	tmp = cur->n;
	delete(cur);
	cur = tmp;
}
first = NULL;
root = NULL;
last = NULL;
}
LFHTEMP bool RBTree<Key,Comparator>::empty(){ return(root == NULL);}
LFHTEMP RBTreeNode<Key,Comparator>* RBTree<Key,Comparator>::find(Key where){
RBTreeNode<Key,Comparator>* tmp =root;
int tmpc;
while((tmp != NULL)&&(tmpc = comp(where,tmp->getIndex()))){
if (tmpc < 0) tmp = tmp->l;
else if (tmpc != 0) tmp = tmp->r;
else return(tmp);
}
return(tmp);
}
LFHTEMP void RBTree<Key,Comparator>::Insert(RBTreeNode<Key,Comparator>* ins){
	RBTreeNode<Key,Comparator>* cur;
	Key insindex = ins->getIndex();
	ins->red = true;
	ins->l = NULL;
	ins->r = NULL;
	if (root == NULL){
		first = last = root =ins;
		ins->red = false;
	}else{
		cur = root;
		while(true){
			if (comp(cur->getIndex(), insindex)>=0){// to the left
//			if (false){
				if (cur->l == NULL){
					cur->l = ins;
					if (cur->b != NULL) cur->b->n = ins;
					ins->b = cur->b;
					cur->b = ins;
					ins->n = cur;
					break;
				}else cur = cur->l;
			}else{
				if (cur->r == NULL){
					cur->r = ins;
					if (cur->n != NULL) cur->n->b = ins;
					ins->n = cur->n;
					cur->n = ins;
					ins->b = cur;
					break;
				}else cur = cur->r;
			}
		}
		ins->p = cur;
		// need to maitain order  and tree balance

		while((cur != NULL)&&(cur->red)){
			if (cur->p == NULL) {cur->red = false; break;}
			if (ins == cur->l){
				if (cur->p->l == cur){ // rigth rotation
					ins->red = false;
					cur = cur->p;
					ins = ins->p;
					ins->p = cur->p;
					cur->p = ins;
					cur->l = ins->r;
					ins->r = cur;
					if (cur->l != NULL) cur->l->p = cur;
					if (ins->p == NULL) root = ins;
					else if (ins->p->l == cur) ins->p->l = ins;
					else ins->p->r = ins;
					cur = ins->p;
				}else{ // double rotation
					cur->red =false;
					cur->l = ins->r;
					if (ins->r != NULL) ins->r->p = cur;
					ins->r = cur;
					cur->p->r = ins->l;
					if (ins->l != NULL) ins->l->p = cur->p;
					ins->l = cur->p;
					ins->p = ins->l->p;
					cur->p = ins;
					cur = ins->p;
					ins->l->p = ins;
					if (cur == NULL) root = ins;
					else if (cur->l == ins->l) cur->l = ins;
					else cur->r = ins;
				}
			}else{
				if (cur->p->l == cur){ // double rotation
					cur->red =false;
					cur->r = ins->l;
					if (ins->l != NULL) ins->l->p = cur;
					ins->l = cur;
					cur->p->l = ins->r;
					if (ins->r != NULL) ins->r->p = cur->p;
					ins->r = cur->p;
					ins->p = ins->r->p;
					cur->p = ins;
					cur = ins->p;
					ins->r->p = ins;
					if (cur == NULL) root = ins;
					else if (cur->l == ins->r) cur->l = ins;
					else cur->r = ins;
				}else{ // left rotation
					ins->red = false;
					cur = cur->p;
					ins = ins->p;
					ins->p = cur->p;
					cur->p = ins;
					cur->r = ins->l;
					ins->l = cur;
					if (cur->r != NULL) cur->r->p = cur;
					if (ins->p == NULL) root = ins;
					else if (ins->p->l == cur) ins->p->l = ins;
					else ins->p->r = ins;
					cur = ins->p;
				}
			}
		}
		root->red = false;
		if (first->b != NULL) first = first->b;
		if (last->n != NULL) last = last->n;
		while(root->p != NULL) root= root->p;
	}

}
LFHTEMP RBTreeNode<Key,Comparator>* RBTree<Key,Comparator>::Top(){
	return(first);
}
LFHTEMP RBTreeNode<Key,Comparator>* RBTree<Key,Comparator>::Pop(){
	RBTreeNode<Key,Comparator> *tmp = first;
	this->Remove(tmp);
	return(tmp);
}
LFHTEMP void RBTree<Key,Comparator>::Delete(Key where){

	RBTreeNode<Key,Comparator>* what = find(where);
	if (what != NULL) {
	this->Remove(what);
	delete(what);
	}
}
LFHTEMP void RBTree<Key,Comparator>::Remove(RBTreeNode<Key,Comparator>* what){
//	display(Dbg);
//	fflush(Dbg);
	RBTreeNode<Key,Comparator>* tofree;

	if ((what->l == NULL)&&(what->r == NULL)) tofree = what;
	else if ((what->l == NULL)||(what->r == NULL)){
		if (what->l == NULL) tofree = what->r;
		else tofree = what->l;
	}else if (rand() & 64){
		if (what->n != NULL) tofree = what->n;
		else tofree = what->b;
	}else{
		if (what->b != NULL) tofree = what->b;
		else tofree = what->n;
	}

	if ((tofree->l != NULL)||(tofree->r != NULL)){
		// not a close leaf, must do a swap first, but it's already red so not need for coloring
		if (tofree->r != NULL){
			tofree->r->red = false;
			tofree->r->p = tofree->p;
			if (tofree->p == NULL) root = tofree->r;
			else if (tofree->p->l == tofree) tofree->p->l = tofree->r;
			else tofree->p->r = tofree->r;
		}else{
			tofree->l->red = false;
			tofree->l->p = tofree->p;
			if (tofree->p == NULL) root = tofree->l;
			else if (tofree->p->l == tofree) tofree->p->l = tofree->l;
			else tofree->p->r = tofree->l;
		}
	}else {
		MakeDBlack(tofree);
	if (tofree->p == NULL) root = NULL;
	else if (tofree->p->l == tofree) tofree->p->l = NULL;
	else tofree->p->r = NULL;
		}
	if (tofree != what) {
		// swap
		tofree->red = what->red;
		tofree->p =what->p;
		if (what->p == NULL) root = tofree;
		else if (what->p->l == what) what->p->l = tofree;
		else what->p->r = tofree;


		tofree->l = what->l;
		if (what->l != NULL) what->l->p = tofree;
		tofree->r = what->r;
		if (what->r != NULL) what->r->p = tofree;


	}
	if (what->b != NULL) what->b->n = what->n;
	else first = what->n;
	if (what->n != NULL) what->n->b = what->b;
	else last = what->b;

}
LFHTEMP void RBTree<Key,Comparator>::MakeDBlack(RBTreeNode<Key,Comparator>* tofree){
	RBTreeNode<Key,Comparator>* tmp;

	while(true){
		tmp = tofree->p;
		if ((tofree->red == true)||(tmp == NULL)) {tofree->red = false; break;}
		if (tmp->l == tofree){
			if (tmp->r->red == true) {
				// we want it black, rotating!
				tmp->red = true;
				tmp = tmp->r;
				tmp->red = false;
				tmp->p = tofree->p->p;
				if (tofree->p->p == NULL) root = tmp;
				else if (tofree->p->p->l == tofree->p) tofree->p->p->l = tmp;
				else tofree->p->p->r = tmp;
				tofree->p->r = tmp->l;
				if (tmp->l != NULL) tmp->l->p = tofree->p;
				tofree->p->p = tmp;
				tmp->l = tofree->p;
			}
			tofree = tofree->p;
			tmp = tofree->r;
			if ( ((tmp->r == NULL)||(tmp->r->red == false)) && ((tmp->l == NULL)||(tmp->l->red == false)) ){
				// no rotation! not done thought
			//	if (tmp != NULL)
					tmp->red = true;
			}else if ((tmp->r != NULL)&&(tmp->r->red)){
				// a rotation, and we are done
				tmp->red = tofree->red;
				tmp->r->red= false;
				tofree->red = false;
				tmp->p = tofree->p;
				if (tofree->p == NULL) root = tmp;
				else if (tofree->p->l == tofree) tofree->p->l = tmp;
				else tofree->p->r = tmp;
				tofree->r = tmp->l;
				if (tmp->l != NULL) tmp->l->p = tofree;
				tmp->l = tofree;
				tofree->p = tmp;
				break;
			}else{
				// double rotation, and we are done >_<
				tmp->l->red = tofree->red;
				tofree->red= false;

				tmp->l->p = tofree->p;
				if (tofree->p == NULL) root = tmp->l;
				else if (tofree->p->l == tofree) tofree->p->l = tmp->l;
				else tofree->p->r = tmp->l;

				tofree->r = tmp->l->l;
				if (tofree->r != NULL) tofree->r->p = tofree;
				tmp->l->l = tofree;
				tofree->p = tmp->l;
				tmp->p = tofree->p;
				tmp->l = tofree->p->r;
				if (tmp->l != NULL) tmp->l->p = tmp;
				tofree->p->r = tmp;
				break;
			}
		}else{
			if (tmp->l->red == true) {
				// we want it black, rotating!
				tmp->red = true;
				tmp = tmp->l;
				tmp->red = false;
				tmp->p = tofree->p->p;
				if (tofree->p->p == NULL) root = tmp;
				else if (tofree->p->p->l == tofree->p) tofree->p->p->l = tmp;
				else tofree->p->p->r = tmp;
				tofree->p->l = tmp->r;
				if (tmp->r != NULL) tmp->r->p = tofree->p;
				tofree->p->p = tmp;
				tmp->r = tofree->p;
			}
			tofree = tofree->p;
			tmp = tofree->l;
			if ( ((tmp->r == NULL)||(tmp->r->red == false)) && ((tmp->l == NULL)||(tmp->l->red == false)) ){
				// no rotation! not done thought
			//	if (tmp != NULL)
					tmp->red = true;
			}else if ((tmp->l != NULL)&&(tmp->l->red)){
				// a rotation, and we are done :-/
				tmp->red = tofree->red;
				tmp->l->red= false;
				tofree->red = false;
				tmp->p = tofree->p;
				if (tofree->p == NULL) root = tmp;
				else if (tofree->p->l == tofree) tofree->p->l = tmp;
				else tofree->p->r = tmp;
				tofree->l = tmp->r;
				if (tmp->r != NULL) tmp->r->p = tofree;
				tofree->p = tmp;
				tmp->r = tofree;
				break;
			}else{
				// double rotation, and we are done >_<
				tmp->r->red = tofree->red;
				tofree->red= false;
				tmp->r->p = tofree->p;
				if (tofree->p == NULL) root = tmp->r;
				else if (tofree->p->l == tofree) tofree->p->l = tmp->r;
				else tofree->p->r = tmp->r;
				tofree->l = tmp->r->r;
				if (tofree->l != NULL) tofree->l->p = tofree;
				tmp->r->r = tofree;
				tofree->p = tmp->r;
				tmp->p = tofree->p;
				tmp->r = tofree->p->l;
				if (tmp->r != NULL) tmp->r->p = tmp;
				tofree->p->l = tmp;
				break;
			}
		}
	}
}
LFHTEMP void RBTree<Key,Comparator>::display(FILE* where){

	char buffer[256];
	memset(buffer,' ',256*sizeof(char));
	buffer[255] = '\0';
	if (root == NULL){
		fprintf(where,"empty Tree\n");
	}else{
		root->recursiveDisplay(where,&buffer[255]);
	}

}
LFHTEMP void RBTree<Key,Comparator>::display(FILE* where, RBTreeNode<Key,Comparator> *mark){

	char buffer[256];
	memset(buffer,' ',256*sizeof(char));
	buffer[255] = '\0';
	if (root == NULL){
		fprintf(where,"empty Tree\n");
	}else{
		root->recursiveDisplay(where,&buffer[255],mark);
	}

}
LFHTEMP RBTreeNode<Key,Comparator>::RBTreeNode(): l(NULL),r(NULL),p(NULL),b(NULL),n(NULL),red(true){
}
LFHTEMP void RBTreeNode<Key,Comparator>::clear(){
	red =true;
	p =NULL;
	b =NULL;
	n =NULL;
	l =NULL;
	r =NULL;
}
LFHTEMP void RBTreeNode<Key,Comparator>::recursiveDisplay(FILE* where,char* buffer){
	if (l != NULL) l->recursiveDisplay(where,buffer-1);
	if ( ((l != NULL)&&(l->p != this)) || ((r != NULL)&&(r->p != this)) || ((b != NULL)&&(b->n != this)) || ((n != NULL)&&(n->b != this)) ) fprintf(where,"%s%c bugged\n",buffer,((red) ? 'R' : 'B'));
	else fprintf(where,"%s%c\n",buffer,((red) ? 'R' : 'B'), ((int)this) );
	if (r != NULL) r->recursiveDisplay(where,buffer-1);
}
LFHTEMP void RBTreeNode<Key,Comparator>::recursiveDisplay(FILE* where,char* buffer, RBTreeNode<Key,Comparator> *mark){
	if (l != NULL) l->recursiveDisplay(where,buffer-1);
	if (this == mark){
	if (((l != NULL)&&(l->p != this))||((r != NULL)&&(r->p != this)) || ((b != NULL)&&(b->n != this)) || ((n != NULL)&&(n->b != this))) fprintf(where,"%s%c target bugged\n",buffer,((red) ? 'R' : 'B'));
	else fprintf(where,"%s%c target\n",buffer,((red) ? 'R' : 'B'), ((int)this) );
	}else{
	if (((l != NULL)&&(l->p != this))||((r != NULL)&&(r->p != this)) || ((b != NULL)&&(b->n != this)) || ((n != NULL)&&(n->b != this))) fprintf(where,"%s%c bugged\n",buffer,((red) ? 'R' : 'B'));
	else fprintf(where,"%s%c\n",buffer,((red) ? 'R' : 'B'), ((int)this) );
	}
	if (r != NULL) r->recursiveDisplay(where,buffer-1);
}

#undef LFHTEMP
#define LFHTEMP template<class Key, class C>

LFHTEMP RBTofDoom<Key,C>::RBTofDoom(const RBTofDoom<Key,C>& other){
    min = other.min;
    max = other.max;
    makeEmptyTree(other.getSize());
    unsigned int buffer[256];
    unsigned int obuffer[256];
    tree[inorderFirst(buffer)].first = other.tree[other.inorderFirst(obuffer)].first;
    unsigned int i;
    for(i=1;i < size;i++) tree[inorderNext(buffer)].first = other.tree[other.inorderNext(obuffer)].first;
}
LFHTEMP RBTofDoom<Key,C>& RBTofDoom<Key,C>::operator=(const RBTofDoom<Key,C>& other){
    this->toMemfree();
    min = other.min;
    max = other.max;
    makeEmptyTree(other.getSize());
    unsigned int buffer[256];
    unsigned int obuffer[256];
    tree[inorderFirst(buffer)].first = other.tree[other.inorderFirst(obuffer)].first;
    unsigned int i;
    for(i=1;i < size;i++) tree[inorderNext(buffer)].first = other.tree[other.inorderNext(obuffer)].first;
return *this;}

LFHTEMP void RBTofDoom<Key,C>::makeEmptyTree(unsigned int nbelem){ // initialize a most balanced tree, that is to be filled!
		size = nbelem;
		alloc_mag = ExOp::upperbound_pow_of_2(nbelem+2);
		if (alloc_mag < 3){

			switch(nbelem){
				case 0: tree = NULL; path =NULL;
				break;
				case 1:
					tree = new pair<Key, unsigned int>[4]; tree--;
					path = new unsigned int[10];
					tree[1].second = 0;
                    firstlone =4;
				break;
				case 2:
					tree = new pair<Key, unsigned int>[4]; tree--;
					path = new unsigned int[10];
					tree[1].second = 4;
					tree[4].second = 1;
                    firstlone =3;
				break;
				}

			}else{

		//       B
		//   B       B
		// R   R   R   R
		//B B B B B B B B
		// ????

		tree = new pair<Key, unsigned int>[1 << alloc_mag];
        if (tree == NULL) {printf("RBTofDoom: could not allocate!\n"); LFH_exit(1);}
		path = new unsigned int[alloc_mag*2+6];
        if (path == NULL) {printf("RBTofDoom: could not allocate!\n"); LFH_exit(1);}
		firstlone = 1 << alloc_mag;
		tree--;
		unsigned int c=1;
		unsigned int m=2;
		unsigned int r,y;

	//	unsigned int cum;
		unsigned int inc;

		for(unsigned char h=0;h<alloc_mag-2;h++){
			r = ((c != 1)&&( ((int)h) % 3 == (((int)alloc_mag) % 3)  )) ? 1 : 0;
			for(;c<m;c++) tree[c].second = (c<<1) | r  ;
			m <<=1;
			}



		unsigned int nbadd = 1 + nbelem - (1 << (alloc_mag-1)); //printf("ndadd = %i\n", nbadd);

		inc =0;
		if (nbadd <= (((unsigned int)1) << (alloc_mag-2))){
			// only loners!
		//	toins_double = 4 << (alloc_mag-3);
            //    printf("hello!\n");
				for(;c<m;c++) {
					inc += nbadd;
			//		printf("%i : %i %i tfund\n",c, inc >> (alloc_mag-2), inc);
					if (inc >> (alloc_mag-2)){
					tree[c].second = firstlone;
					tree[firstlone].second = c;
					firstlone--;
                    inc -= 1 << (alloc_mag-2);
					} else tree[c].second = 0;
					}

			}else{
			nbadd -= (1 << (alloc_mag-2));
			y = (4 << (alloc_mag-3));
		//	toins_double = y+ 2 * nbadd;
				for(;c<m;c++) {
					inc += nbadd;
				//	printf("%i : %i %i fun\n",c, inc >> (alloc_mag-2), inc);
					if (inc >> (alloc_mag-2)){
					tree[c].second = y;
					tree[y].second = 1;
					tree[y+1].second = 1;
					y+= 2;
                    inc -= 1 << (alloc_mag-2);
					} else{
					tree[c].second = firstlone;
					tree[firstlone].second = c;
					firstlone--;
					}
					}

		}
    }
}



LFHTEMP	Key RBTofDoom<Key,C>::getMin() const{return min;}

LFHTEMP	Key RBTofDoom<Key,C>::getMax() const{return max;}

LFHTEMP	void RBTofDoom<Key,C>::batchInit(Vector<Key> &batch, bool needsort){
    this->toMemfree();
    this->batchInit_routine(batch, needsort);
}

LFHTEMP	void RBTofDoom<Key,C>::batchInit_routine(Vector<Key> &batch, bool needsort){
    if (needsort) batch.sort();

    this->makeEmptyTree(batch.size());
    unsigned int buffer[1024];

    min = tree[inorderFirst(buffer)].first = batch[0];
    unsigned int i;
    for(i=1;i < size-1;i++) tree[inorderNext(buffer)].first = batch[i];
    max = tree[inorderNext(buffer)].first = batch[size-1];
}
LFHTEMP	bool RBTofDoom<Key,C>::reach(const Key &q){ // stops on equality
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (q != tree[cur].first) {
					if (tree[tree[cur].second].first == q){
						path[++(path[0])] =tree[cur].second;
						return(true);
					} else return false;
					}
					return(true);

				}
				if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
				else if (q != tree[cur].first) cur = tree[cur].second | 1;
				else return true;
				path[++(path[0])] =cur;
			}
			return (q == tree[cur].first);
		}
		return false;
	}
// loop throught the remaining equal nodes, (ignore the highest in tree, it was already found)
LFHTEMP	bool RBTofDoom<Key,C>::reach_continue(const Key &q){ // stops on equality
			// go right!
		/*	for(cur = path[0] ; cur >0 ;cur--) if (tree[path[cur]].first == q) break;
			if ((cur == 0)&&(tree[0].first != q)){
				//reach minimum equal in tree

			}else{

			}*/
			unsigned int cur = path[0] == 0 ? 1 : path[path[0]];
			if (tree[cur].second < 2) return false;
			cur = tree[cur].second | 1;
			path[++(path[0])] =cur;
				while(tree[cur].second > 1){
					if (tree[cur].second > firstlone) {

						if (q != tree[cur].first) {
							if (tree[tree[cur].second].first == q){
								path[++(path[0])] =tree[cur].second;
								return(true);
							} else return false;
						}
						return(true);

					}
					if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
						else if (q != tree[cur].first) cur = tree[cur].second | 1;
							else return true;
					path[++(path[0])] =cur;
				}
				return (q == tree[cur].first);

			return false;
		}
LFHTEMP	void RBTofDoom<Key,C>::reach_par(const Key &q){ // ignores equalities!
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) break;
				if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
				else cur = tree[cur].second | 1;
				path[++(path[0])] =cur;
			}
		} else path[0] = 0;
	}
LFHTEMP	template<class A> bool RBTofDoom<Key,C>::reach(const A &q){ // stops on equality
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (q != tree[cur].first) {
					if (ExOp::isEQ(tree[tree[cur].second].first,q)){
						path[++(path[0])] =tree[cur].second;
						return(true);
					} else return false;
					}
					return(true);

				}
				if (ExOp::isLT(q, tree[cur].first)) cur = tree[cur].second & 0xFFFFFFFE;
				else if (ExOp::isNQ(q, tree[cur].first)) cur = tree[cur].second | 1;
				else return true;
				path[++(path[0])] =cur;
			}
			return (ExOp::isEQ(q,tree[cur].first));
		}
		return false;
	}
LFHTEMP	unsigned int RBTofDoom<Key,C>::find_index(const Key &what) const{ // 0 for not found
    if (size > 0){
        unsigned int cur=1;
        while(tree[cur].second > 1){
            if (tree[cur].second > firstlone) {

                if (what != tree[cur].first) {
                    if (tree[tree[cur].second].first == what){
                        return tree[cur].second;
                    } else return 0;
                }
                return cur;
            }
            if (what < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
            else if (what != tree[cur].first) cur = tree[cur].second | 1;
            else return cur;
        }

        return (what == tree[cur].first) ? cur : 0;
    }
    return 0;
    }
LFHTEMP	template<class O_KEY> unsigned int RBTofDoom<Key,C>::find_index(const O_KEY &what) const{ // 0 for not found
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (tree[cur].first == what) return cur;
					if (tree[tree[cur].second].first == what){
						return tree[cur].second;
					} else return 0;


				}
				if (tree[cur].first > what) cur = tree[cur].second & 0xFFFFFFFE;
				else if (tree[cur].first != what) cur = tree[cur].second | 1;
				else return cur;
			}

			return (tree[cur].first == what) ? cur : 0;
		}
		return 0;
		}
	// regarless if "what" exists, find the previous and next non-equal element to it
LFHTEMP	void RBTofDoom<Key,C>::findPN(const Key &what, unsigned int &prev, unsigned int &next)const{
			prev = next = 0;
			if (size > 0){
			unsigned int cur=1;

			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

						if (what < tree[cur].first) next = cur;
						else{
							if (what != tree[cur].first) prev = cur;
							if (what < tree[tree[cur].second].first) next =tree[cur].second;
							else if (tree[tree[cur].second].first != what) prev = tree[cur].second;
						}
					return;
				}

				if (what < tree[cur].first) {next = cur;cur = tree[cur].second & 0xFFFFFFFE;  }
				else if (what != tree[cur].first) {prev = cur;cur = tree[cur].second | 1;  }
				else { // found item... need to look into both subtrees
					// todo
					next = tree[cur].second | 1;
				//	while(tree[next].second > 1){

				//		}

					prev = tree[cur].second & 0xFFFFFFFE;
					// todo
				//	while(tree[prev].second > 1){

				//		}
					return;
					}
			}

			if (what < tree[cur].first) {next = cur;cur = tree[cur].second & 0xFFFFFFFE;  }
			else if (what != tree[cur].first) {prev = cur;cur = tree[cur].second | 1;  }


		}


		}
LFHTEMP	void RBTofDoom<Key,C>::insert(const Key &newf){
		unsigned int i;
		bool cached;
		Key swapk;
		if (size > 0){
			if (newf < min) min = newf ;
			else if (newf > max)  max = newf;

			if (size+2 == (1u << alloc_mag) ){
				Vector<Key> batch;
				batch.setSize(size+1);

				//	show(stdout);
				//		printf("inserting %i, %i\n", newf, tree[inorderFirst(path)].first);
				i =0;
				batch[i] = tree[inorderFirst(path)].first;
				while(true){
					if (newf < batch[i]){
						batch[i+1] = batch[i];
						batch[i] = newf;
						i++;
						break;
					}
					i++;
					if (i == size){
						batch[size] = newf;
						break;
					}
					batch[i] = tree[inorderNext(path)].first;
				}
				for(i++;i<= size;i++) { batch[i] = tree[inorderNext(path)].first;}

				delete[](tree+1);
				delete[](path);
//printf("reloaded updated AAA!\n");
				batchInit_routine(batch,false);
			}else{
				reach_par(newf);
				unsigned int cur = path[0] == 0 ? 1 : path[path[0]];
				unsigned int cur2,cur3;
				size++;
				if (tree[cur].second >= firstlone){

					cur2 = tree[cur].second;
					cur3 = size - (1 << alloc_mag) + firstlone;

					tree[tree[cur2].second].second = cur3; // parent should be black

					LFH_ASSERT(cur3 +1 < firstlone, "Corrupted RBTree!\n");
                    LFH_ASSERT((cur3 & 1) == 0, "Corrupted RBTree!\n");
					// get next free slot!

					tree[cur3].second = 1;
					tree[cur3+1].first = tree[cur2].first;
					tree[cur3+1].second = 1;




					firstlone++;
					if (cur2 != firstlone) {
						tree[cur2] = tree[firstlone];
						tree[tree[cur2].second].second = cur2;
					}

					if (newf < tree[cur].first){

						tree[cur3].first = newf;

					}else{
						tree[cur3].first =  tree[cur].first;
						if (newf < tree[cur3+1].first){
							tree[cur].first = newf;

						}else{
							tree[cur].first = tree[cur3+1].first;
							tree[cur3+1].first = newf;
						}
					}

				}else{ // add loner
					//	printf("helloB!");
					cached = (tree[cur].second == 1);
					tree[firstlone].second = cur;
					tree[cur].second = firstlone;

					if (newf < tree[cur].first){
						tree[firstlone].first = tree[cur].first;
						tree[cur].first = newf;
					}else tree[firstlone].first = newf;
					firstlone--;
					if (cached){

						tree[cur ^ 1].second = 0; // paint the other children in black;
						path[0]--;
						if (path[0] >0){
							tree[path[path[0]]].second |=1;
							path[0]--;
							while((path[0]>0)&&(tree[path[path[0]]].second & 1)){
								if ((tree[path[path[0]]^ 1].second & 1)&&(tree[path[path[0]]^ 1].second <=firstlone)){
									tree[path[path[0]]^ 1].second &= 0xFFFFFFFE;
									tree[path[path[0]]].second &= 0xFFFFFFFE;

									path[0]--;
									if (path[0] !=0) {tree[path[path[0]]].second |= 1; path[0]--;}
								}else{
									if ((path[path[0]+1] & 1)^(path[path[0]] & 1)){ // double rotation
										// 		printf("double rot!\n");show(stdout);
										if ((tree[path[path[0]+1]].second < 2)||(tree[path[path[0]+1]].second > firstlone)){
											printf("OMFG!!!!!\n");

										}else{

											path[0]--;
									//			printf("here wego !\n");
											if (path[0] == 0) {
												swapk = tree[1].first;
												tree[1].first = tree[path[2]].first;
												tree[1].second = tree[path[2]].second & 0xFFFFFFFE;
											}else {
												swapk = tree[path[path[0]]].first;
												tree[path[path[0]]]= tree[path[path[0]+2]];
											}

											//		printf("%c\n", (path[path[0]+2] & 1) ? 'Y' : 'N');
											if (path[path[0]+2] & 1) path[path[0]+3] = tree[ path[path[0]+2] ].second | 1;
											else path[path[0]+3] = tree[ path[path[0]+2] ].second  & 0xFFFFFFFE;

											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+2]].first = tree[path[path[0]+3]  ^ 1].first;
											tree[path[path[0]+2]].second = tree[path[path[0]+3]  ^ 1].second;
											if (tree[path[path[0]+2]].second > firstlone) tree[tree[path[path[0]+2]].second].second = path[path[0]+2];

											tree[path[path[0]+3]  ^ 1].first = tree[path[path[0]+1]].first;
											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+1]].first = tree[ path[path[0]+3] ].first;
											tree[path[path[0]+1]].second = tree[ path[path[0]+3] ].second;
											if (tree[path[path[0]+1]].second > firstlone) tree[tree[path[path[0]+1]].second].second = path[path[0]+1];

											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+3] ].first = swapk;

											tree[ path[path[0]+3] ].second = path[path[0]+1] & 0xFFFFFFFE;
											tree[ path[path[0]+3] ^ 1].second = path[path[0]+2] & 0xFFFFFFFE;

											//			printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);




										}


										if (path[0] != 0) { path[0]--;}

										// 		printf("is that fine? %i, %i, %i\n", path[0], path[1], path[2] );show(stdout);




									}else{ // rotation
										//printf("here wego rot !\n");
										path[0]--;
										if (path[0] == 0) {
											swapk = tree[1].first;
											tree[1].first = tree[path[1]].first;
											tree[1].second = path[path[0]+2] &0xFFFFFFFE;
										}else {
											swapk = tree[path[path[0]]].first;
											tree[path[path[0]]].first = tree[path[path[0]+1]].first;
											tree[path[path[0]]].second = path[path[0]+2] | 1;
										}


										tree[path[path[0]+1]].first = tree[path[path[0]+2] ^ 1].first;

										if (tree[path[path[0]+2] ^ 1].second > firstlone){
											tree[path[path[0]+1]].second = tree[path[path[0]+2] ^ 1].second;
											tree[tree[path[path[0]+2] ^ 1].second].second = path[path[0]+1];
										}else tree[path[path[0]+1]].second = tree[path[path[0]+2] ^ 1].second & 0xFFFFFFFE;
										//	if (tree[path[path[0]+2] ^ 1].second > firstlone) tree[tree[path[path[0]+1]].second].second = path[path[0]+1];

										tree[path[path[0]+2]^ 1].first = swapk;
										tree[path[path[0]+2]^ 1].second = path[path[0]+1] & 0xFFFFFFFE;

										//}
										tree[path[path[0]+2]].second &= 0xFFFFFFFE;
										if (path[0] != 0) { path[0]--;}
										// 	show(stdout);
										//


									}}}}}





				}

			}
		}else{
			alloc_mag = 2;
			tree = new pair<Key, unsigned int>[4]; 			path = new unsigned int[10];

			tree--;
			tree[1].first = newf;
			tree[1].second = 0;
			size =1;
			firstlone=1 << alloc_mag;
	//		toins_double = 2;
			min = newf;
			max = newf;
		}

	}


LFHTEMP	void RBTofDoom<Key,C>::remove(const Key &todelete){ // remove 1 element which is equal, if found
    if ((size > 0)&&(reach(todelete))) remove_maintain_routine(todelete);
}

LFHTEMP	template<class A> void RBTofDoom<Key,C>::remove(const A &todelete){ // remove 1 element which is equal, if found
    if ((size > 0)&&(reach(todelete))) {
        this->remove_maintain_routine(tree[path[0] == 0 ? 1 : path[path[0]]].first);
    }
}

LFHTEMP void RBTofDoom<Key,C>::remove_maintain_routine(const Key &todelete){
    bool prob;
    unsigned int cur,cur2,cur3;
    unsigned int ipath;
    Key swap;
    if (size == 1) this->toMemfree();
    else{
    cur = path[0] == 0 ? 1 : path[path[0]];
    if (todelete == min){
        typename RBTofDoom<Key,C>::Iterator ite(*this); ite.findFirst();
        ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
        if (ipath == cur){
            if (tree[cur].second > firstlone){
                min = tree[tree[cur].second].first;
            }else min = tree[ path[0] == 1 ? 1 : path[path[0]-1] ].first;
        }
    }else if (todelete == max) {
        typename RBTofDoom<Key,C>::Iterator ite(*this); ite.findLast();
        ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
        if (ipath == cur) max = tree[ path[0] == 1 ? 1 : path[path[0]-1] ].first;
    }
    if (size == (((unsigned int) 1) << (alloc_mag-2)) ){ // Downsize, complete rebalance
        Vector<Key> batch;
        batch.setSize(size-1);
        cur =0;
        batch[cur] = tree[inorderFirst(path)].first;
        while(true){
            if (todelete == batch[cur]) break;
            cur++;
            if (cur == size-1) break;
            batch[cur] = tree[inorderNext(path)].first;
        }

        for(;cur< size-1;cur++) { batch[cur] = tree[inorderNext(path)].first;}
        delete[](tree+1);
        delete[](path);
        batchInit_routine(batch,false);
    }else{
    size--;
    cur = path[0] == 0 ? 1 : path[path[0]];
    if (!isLeaf(cur)){ // swap node... preserves structure
        // swap node with a leaf!
        cur2 = cur;
        if (!isLeaf(cur = tree[cur].second)) {
            cur = cur & 0xFFFFFFFE; // left child!
            while(tree[cur].second > 1){
                path[++path[0]] = cur;
                if (tree[cur].second > firstlone) {
                    cur = tree[cur].second;
                    path[++path[0]] = cur;
                    break;
                }
                cur = (tree[cur].second) | 1; // right child!
            }
        }
        path[++path[0]] = cur;
        ExOp::toMemmove(tree[cur2].first,tree[cur].first);
    }
    if (cur > firstlone){
        firstlone++;
        tree[tree[cur].second].second = 0;
        if (firstlone != cur) {
        tree[cur] = tree[firstlone];
        tree[tree[firstlone].second].second = cur;
        }
    }else{
        cur2 = path[0] == 1 ? 1 : path[path[0]-1];
        if ((tree[cur].second == 0)&&(tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // annoying red brother... ROTATE!
            cur3 = (cur & 1) ? tree[cur^1].second: tree[cur^1].second & 0xFFFFFFFE;

            swap = tree[cur2].first;
            ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
            ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
            tree[cur3].first = swap;

            tree[cur ^ 1].second = tree[cur3].second;
            if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

            tree[cur2].second = cur3 & 0xFFFFFFFE;
            tree[cur3].second = cur | 1;

            path[path[0]+1] = path[path[0]];
            path[path[0]] = cur3;
            path[0]++;

            cur = path[path[0]];
            cur2 = path[path[0]-1];
        }



        if (tree[cur^1].second > firstlone){ // brother has lone child, use him!

            firstlone++;
            cur3= tree[cur^1].second;

            if (cur & 1){
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
            }else{
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
            }
            if (cur3 != firstlone){
            ExOp::toMemmove(tree[cur3].first,tree[firstlone].first);
            tree[cur3].second = tree[firstlone].second;
            tree[tree[firstlone].second].second = cur3;
            }
            tree[cur ^ 1].second =0;

        }else{
            // deleting a inside leaf...
    //		printf("DELIN!!! %i \n", tree[cur^1].second);
            cur2 = path[0] == 1 ? 1 : path[path[0]-1];




            if (tree[cur^1].second > 1){ // brother has siblings, do a rotation
    //			printf("NOT oK!!!! %i \n", tree[cur^1].second);
                prob = false; // done after this
                ipath = tree[cur^1].second & 0xFFFFFFFE; // to fill later

                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                if (cur & 1){
                    ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second | 1 ].first);
                    ExOp::toMemmove(tree[firstlone].first,tree[cur^1].first);
                    ExOp::toMemmove(tree[cur^1].first,tree[tree[cur^1].second].first);
                }else{
                    ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second ].first);
                    ExOp::toMemmove(tree[firstlone].first,tree[tree[cur^1].second | 1].first);
                }
                tree[cur^1].second = firstlone;
                tree[firstlone].second = cur^1;
            }else{
            prob = (tree[cur2].second & 1) == 0; // is black?
            // send brother to loners
            tree[firstlone].second = cur2;
            tree[cur2].second = firstlone;

            if (cur & 1){
                // needs to swap parent
                ExOp::toMemmove(tree[firstlone].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
            }else ExOp::toMemmove(tree[firstlone].first,tree[cur ^ 1].first);

                ipath = cur & 0xFFFFFFFE; // to fill later
                path[0]--;
            }
            firstlone--;
            if ((prob)&&(tree[cur].second == 0)){ // both child and parent were black, needs rotations! (cur2 is double-black)
                if (path[0] != 0)
                while(true){
                    // rotate to guarranty black brother
                    if (tree[path[path[0]] ^ 1].second & 1){
                        cur = path[0] == 1 ? 1 : path[path[0]-1];
                        swap = tree[cur].first;
                        ExOp::toMemmove(tree[cur].first,tree[path[path[0]] ^ 1].first);
                        cur2 = (path[path[0]] & 1) ? tree[path[path[0]] ^ 1].second | 1 : tree[path[path[0]] ^ 1].second & 0xFFFFFFFE;
                        ExOp::toMemmove(tree[path[path[0]] ^ 1].first, tree[cur2].first);
                        tree[cur2].first = swap;

                        tree[path[path[0]] ^ 1].second = tree[cur2].second;
                        tree[cur].second = cur2 &0xFFFFFFFE;
                        tree[cur2].second =  path[path[0]] | 1;

                        path[path[0]+1] = path[path[0]];
                        path[path[0]] = cur2;
                        path[0]++;
                    }

                    cur = (path[path[0]] & 1) ? tree[path[path[0]] ^ 1].second | 1 : tree[path[path[0]] ^ 1].second & 0xFFFFFFFE ;
                    if ((tree[cur^1].second & 1)&&(tree[cur^1].second  <=firstlone)) {
                        // one rotation, then done!

                        cur2 = path[0] == 1 ? 1 : path[path[0]-1];
                        swap = tree[cur2].first;
                        ExOp::toMemmove(tree[cur2].first,tree[path[path[0]] ^ 1].first);
                        ExOp::toMemmove(tree[path[path[0]] ^ 1].first, tree[cur].first);
                        tree[cur].first = swap;

                        tree[path[path[0]] ^ 1].second = tree[cur].second;
                        if (tree[cur].second > firstlone) tree[tree[cur].second].second = path[path[0]] ^ 1;
                        tree[cur2].second = (tree[cur2].second & 1) | (cur & 0xFFFFFFFE);
                        tree[cur].second =  path[path[0]] & 0xFFFFFFFE;

                        tree[cur^1].second &= 0xFFFFFFFE;

                        break;
                    }else if ((tree[cur].second & 1)&&(tree[cur].second  <=firstlone)){
                        // double rotation, then done!
                        cur2 = path[0] == 1 ? 1 : path[path[0]-1];
                        swap = tree[cur2].first;
                        ExOp::toMemmove(tree[cur2].first,tree[cur].first);
                        cur3 = (path[path[0]] & 1) ? tree[cur].second | 1 : tree[cur].second & 0xFFFFFFFE ;
                        ExOp::toMemmove(tree[cur].first, tree[cur3 ^ 1].first );
                        ExOp::toMemmove(tree[cur3 ^ 1].first, tree[path[path[0]] ^ 1].first );
                        ExOp::toMemmove(tree[path[path[0]] ^ 1].first, tree[cur3].first);
                        tree[cur3].first = swap;

                        tree[path[path[0]] ^ 1].second = tree[cur3].second;
                        if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = path[path[0]] ^ 1;
                        tree[cur].second = tree[cur3^1].second;
                        if (tree[cur3^1].second > firstlone) tree[tree[cur3^1].second].second = cur;

                        tree[cur3].second = path[path[0]] & 0xFFFFFFFE;
                        tree[cur3^1].second = cur & 0xFFFFFFFE;
                        tree[cur2].second = (tree[cur2].second & 1) | (cur3 & 0xFFFFFFFE);



                        break;
                        }
                    // no rotation... not finished D:
                    //	cur3 = tree[path[path[0]]].first;
                    //	tree[path[path[0]]].first = -1;
                    //	show();
                    //	tree[path[path[0]]].first = cur3;

                        tree[path[path[0]] ^ 1].second |= 1; // paint brother red
                        path[0]--;

                        if (path[0] == 0) break;
                        if (tree[path[path[0]]].second & 1){
                            tree[path[path[0]]].second &= 0xFFFFFFFE;
                            break;
                            }
                    }
            }

            cur3 = size - (1 << alloc_mag) + firstlone + 1; // last pair position

            if (cur3 != ipath){

                typename RBTofDoom<Key,C>::Iterator ite(*this, tree[cur3].first );
                for(; (ite.path[0] == 0) || (ite.path[ite.path[0]] != cur3) ; ++ite){
                }


                ExOp::toMemmove(tree[ipath].first, tree[cur3].first );
                tree[ipath].second = tree[cur3].second;
                ExOp::toMemmove(tree[ipath | 1].first, tree[cur3 | 1].first );
                tree[ipath| 1].second = tree[cur3| 1].second;

                if (tree[ipath].second > firstlone) tree[tree[ipath].second].second = ipath;
                if (tree[ipath|1].second > firstlone) tree[tree[ipath|1].second].second = ipath | 1;

                cur3 = (ite.path[0] == 1) ? 1 : ite.path[ite.path[0]-1];
                tree[cur3].second = (tree[cur3].second & 1) | ipath;
            }
        }
    }
    }
    }
}


LFHTEMP	void RBTofDoom<Key,C>::remove_static_routine(typename RBTofDoom<Key,C>::Iterator &todel,  pair<Vector<unsigned int>, unsigned int> &list_and_lone){ // delete pointed node and update iterator which will target some element
    if (!todel.isValid()) return;
    //	display(Dbg);
    //	fflush(Dbg);

    LFH_ASSERT(((size - (1 << alloc_mag) + firstlone) & 1) == 1, "Not legal deleting");

  //  printf("single del on "); ExOp::show(*todel);
  //  printf("firstlone is  %i\n", firstlone);

    bool prob;
    unsigned int cur,cur2,cur3;
    unsigned int ipath;
    Key swap;
    // Maintain min-max
    if (size == 1) {
        list_and_lone.first.push_back(0); list_and_lone.second++;
    }else{
        cur = todel.path[0] == 0 ? 1 : todel.path[todel.path[0]];
        // maintain iterator! this is crazy!
        if (!isLeaf(cur)){ // swap node with previous node, preserves structure
            // swap node with a leaf!
            cur2 = cur;
            if (!isLeaf(cur = tree[cur].second)) {
                cur = cur & 0xFFFFFFFE; // left child!
                while(tree[cur].second > 1){
                    todel.path[++(todel.path[0])] = cur;
                    if (tree[cur].second > firstlone) {
                        cur = tree[cur].second;
                        todel.path[++(todel.path[0])] = cur;
                        break;
                    }
                    cur = (tree[cur].second) | 1; // right child!

                }
            }
            todel.path[++(todel.path[0])] = cur;
            ExOp::toMemmove(tree[cur2].first,tree[cur].first);
        }

        if (cur > firstlone){
             list_and_lone.first.push_back(cur); list_and_lone.second++;
            tree[tree[cur].second].second = 0;
        }else{
            cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
            if ((tree[cur].second == 0)&&(tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // annoying red brother... ROTATE!
                cur3 = (cur & 1) ? tree[cur^1].second: tree[cur^1].second & 0xFFFFFFFE;
                swap = tree[cur2].first;
                ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
                ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
                tree[cur3].first = swap;

                tree[cur ^ 1].second = tree[cur3].second;
                if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

                tree[cur2].second = cur3 & 0xFFFFFFFE;
                tree[cur3].second = cur | 1;

                todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
                todel.path[todel.path[0]] = cur3;
                todel.path[0]++;

                cur = todel.path[todel.path[0]];
                cur2 = todel.path[todel.path[0]-1];
            }



            if (tree[cur^1].second > firstlone){ // brother has lone child, use him!
                cur3= tree[cur^1].second;
                if (cur & 1){
                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
                }else{
                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                    ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
                }
                list_and_lone.first.push_back(cur3); list_and_lone.second++;
                tree[cur ^ 1].second =0;

            }else{
                // deleting a inside leaf...
                if (list_and_lone.second != 0) {
                    for(ipath=0; firstlone>=list_and_lone.first[ipath];ipath++);
                    cur3 = list_and_lone.first[ipath];
                    list_and_lone.first.pop_swap(ipath); list_and_lone.second--;
                }else{cur3 = firstlone--;size++; } // printf("WaRnInG! this is bad news (RBTofDOOm) %i is size\n", size); fflush(stdout);

                cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];

                if (tree[cur^1].second > 1){ // brother has siblings, do a rotation
        //			printf("NOT oK!!!! %i \n", tree[cur^1].second);
                    prob = false; // done after this
                    ipath = tree[cur^1].second & 0xFFFFFFFE; // to fill later

                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);

                    if (cur & 1){
                        ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second | 1 ].first);

                        ExOp::toMemmove(tree[cur3].first,tree[cur^1].first);
                        ExOp::toMemmove(tree[cur^1].first,tree[tree[cur^1].second].first);

                    }else{

                        ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second ].first);
                        ExOp::toMemmove(tree[cur3].first,tree[tree[cur^1].second | 1].first);

                    }
                    tree[cur^1].second = cur3;
                    tree[cur3].second = cur^1;

                }else{
                    prob = (tree[cur2].second & 1) == 0; // is black?

                    // send brother to loners


                    tree[cur3].second = cur2;
                    tree[cur2].second = cur3;

                    if (cur & 1){
                        // needs to swap parent
                        ExOp::toMemmove(tree[cur3].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                    }else ExOp::toMemmove(tree[cur3].first,tree[cur ^ 1].first);
                    ipath = cur & 0xFFFFFFFE; // to fill later
                    todel.path[0]--;
                }

                if ((prob)&&(tree[cur].second == 0)) this->resolve_doubleblack(todel); // both child and parent were black, needs rotations! (cur2 is double-black)

                list_and_lone.first.push_back(ipath);
}   }   }	}
LFHTEMP	void RBTofDoom<Key,C>::remove_update(typename RBTofDoom<Key,C>::Iterator &todelete){
	Key tmp = (*todelete);
	this->remove(todelete);
	todelete = this->find_first(tmp);
}
LFHTEMP	void RBTofDoom<Key,C>::remove_mm(typename RBTofDoom<Key,C>::Iterator &todelete){
	Key thatkey = (*todelete);
	this->remove(todelete);
	todelete.findLE(thatkey);
}
LFHTEMP	void RBTofDoom<Key,C>::remove_pp(typename RBTofDoom<Key,C>::Iterator &todelete){
	Key thatkey = (*todelete);
	this->remove(todelete);
	todelete.findGE(thatkey);
}
LFHTEMP	void RBTofDoom<Key,C>::remove(typename RBTofDoom<Key,C>::Iterator &todel){ // delete pointed node and update iterator which will target some element
    if (!todel.isValid()) return;

    //	display(Dbg);
    //	fflush(Dbg);
    bool prob;
    unsigned int cur,cur2,cur3;
    unsigned int ipath;
    Key swapbuf;
    // Maintain min-max
    if (size == 1) {
        delete[](tree+1);
        tree = NULL;
        size=0;
        todel.path[0] = ExCo<unsigned int>::mkMaximum(); // make unvalid
    }else{
        cur = todel.path[0] == 0 ? 1 : todel.path[todel.path[0]];
        // MAINTAIN MIN/MAX
        if (tree[cur].first == min){
            typename RBTofDoom<Key,C>::Iterator ite(*this); ite.findFirst();
            ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
            if (ipath == cur){
                if (tree[cur].second > firstlone){
                    min = tree[tree[cur].second].first;
                }else min = tree[ todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1] ].first;
            }
        }else if (tree[cur].first == max) {
            typename RBTofDoom<Key,C>::Iterator ite(*this); ite.findLast();
            ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
            if (ipath == cur) max = tree[ todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1] ].first;
        }

        if (size == (((unsigned int) 1) << (alloc_mag-2)) ){ // Downsize, complete rebalance
            Vector<Key> batch;
            batch.setSize(size-1);
            ipath =0;
            cur2 = inorderFirst(path);
            while(true){
                if (cur == cur2) {ExOp::toMemmove(swapbuf,tree[cur2].first); break;}
                ExOp::toMemmove(batch[ipath++], tree[cur2].first);
                cur2 =inorderNext(path);
            }
            while(ipath< size-1) batch[ipath++] = tree[inorderNext(path)].first;
            delete[](tree+1);
            delete[](path);
            batchInit_routine(batch,false);

            todel.findGE(swapbuf);

        }else{ // RBT updates
            size--;
            // maintain iterator! this is crazy!
            if (!isLeaf(cur)){ // swapbuf node with previous node, preserves structure
                // swapbuf node with a leaf!
                cur2 = cur;
                if (!isLeaf(cur = tree[cur].second)) {
                    cur = cur & 0xFFFFFFFE; // left child!
                    while(tree[cur].second > 1){
                        todel.path[++(todel.path[0])] = cur;
                        if (tree[cur].second > firstlone) {
                            cur = tree[cur].second;
                            todel.path[++(todel.path[0])] = cur;
                            break;
                        }
                        cur = (tree[cur].second) | 1; // right child!

                    }
                }
                todel.path[++(todel.path[0])] = cur;
                ExOp::toMemmove(tree[cur2].first,tree[cur].first);
            }

            if (cur > firstlone){
                firstlone++;
                tree[tree[cur].second].second = 0;
                if (firstlone != cur) {
                    tree[cur] = tree[firstlone];
                    tree[tree[firstlone].second].second = cur;
                }
            }else{

                cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
                if ((tree[cur].second == 0)&&(tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // annoying red brother... ROTATE!
                    cur3 = (cur & 1) ? tree[cur^1].second: tree[cur^1].second & 0xFFFFFFFE;

                    ExOp::toMemmove(swapbuf,tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
                    ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
                    ExOp::toMemmove(tree[cur3].first,swapbuf);

                    tree[cur ^ 1].second = tree[cur3].second;
                    if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

                    tree[cur2].second = cur3 & 0xFFFFFFFE;
                    tree[cur3].second = cur | 1;

                    todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
                    todel.path[todel.path[0]] = cur3;
                    todel.path[0]++;

                    cur = todel.path[todel.path[0]];
                    cur2 = todel.path[todel.path[0]-1];
                }



                if (tree[cur^1].second > firstlone){ // brother has lone child, use him!

                    firstlone++;
                    cur3= tree[cur^1].second;

                    if (cur & 1){
                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
                    }else{
                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                        ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
                    }
                    if (cur3 != firstlone){
                    ExOp::toMemmove(tree[cur3].first,tree[firstlone].first);
                    tree[cur3].second = tree[firstlone].second;
                    tree[tree[firstlone].second].second = cur3;
                    }
                    tree[cur ^ 1].second =0;

                }else{
                    // deleting a inside leaf...
                    cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];

                    if (tree[cur^1].second > 1){ // brother has siblings, do a rotation
            //			printf("NOT oK!!!! %i \n", tree[cur^1].second);
                        prob = false; // done after this
                        ipath = tree[cur^1].second & 0xFFFFFFFE; // to fill later

                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        if (cur & 1){
                            ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second | 1 ].first);

                            ExOp::toMemmove(tree[firstlone].first,tree[cur^1].first);
                            ExOp::toMemmove(tree[cur^1].first,tree[tree[cur^1].second].first);

                        }else{

                            ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second ].first);
                            ExOp::toMemmove(tree[firstlone].first,tree[tree[cur^1].second | 1].first);

                        }
                        tree[cur^1].second = firstlone;
                        tree[firstlone].second = cur^1;


        //					firstlone--;show();	firstlone++;
                    }else{
                        prob = (tree[cur2].second & 1) == 0; // is black?

                        // send brother to loners


                        tree[firstlone].second = cur2;
                        tree[cur2].second = firstlone;

                        if (cur & 1){
                            // needs to swap parent
                            ExOp::toMemmove(tree[firstlone].first,tree[cur2].first);
                            ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                        }else ExOp::toMemmove(tree[firstlone].first,tree[cur ^ 1].first);
                        ipath = cur & 0xFFFFFFFE; // to fill later
                        todel.path[0]--;
                    }


                    firstlone--;


                    if ((prob)&&(tree[cur].second == 0)) this->resolve_doubleblack(todel); // both child and parent were black, needs rotations! (cur2 is double-black)

                    cur3 = size - (1 << alloc_mag) + firstlone + 1; // last pair position

                    // move last pair inside
                    if (cur3 != ipath){
                        ExOp::toMemmove(tree[ipath].first, tree[cur3].first );
                        tree[ipath].second = tree[cur3].second;
                        ExOp::toMemmove(tree[ipath | 1].first, tree[cur3 | 1].first );
                        tree[ipath| 1].second = tree[cur3| 1].second;

                        if (tree[ipath].second > firstlone) tree[tree[ipath].second].second = ipath;
                        if (tree[ipath|1].second > firstlone) tree[tree[ipath|1].second].second = ipath | 1;


                        typename RBTofDoom<Key,C>::Iterator ite(*this, tree[cur3].first );
                        for(; (ite.path[0] == 0) || (ite.path[ite.path[0]] != cur3) ; ++ite){

                        }
                        cur3 = (ite.path[0] == 1) ? 1 : ite.path[ite.path[0]-1];
                        tree[cur3].second = (tree[cur3].second & 1) | ipath;
}   }   }   }   }	}

LFHTEMP	void RBTofDoom<Key,C>::resolve_doubleblack(typename RBTofDoom<Key,C>::Iterator &todel){
Key swapbuf;
unsigned int cur, cur2, cur3;

if (todel.path[0] != 0)
while(true){
    // rotate to guarranty black brother
    if (tree[todel.path[todel.path[0]] ^ 1].second & 1){
        cur = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        ExOp::toMemmove(swapbuf, tree[cur].first);
        ExOp::toMemmove(tree[cur].first,tree[todel.path[todel.path[0]] ^ 1].first);
        cur2 = (todel.path[todel.path[0]] & 1) ? tree[todel.path[todel.path[0]] ^ 1].second | 1 : tree[todel.path[todel.path[0]] ^ 1].second & 0xFFFFFFFE;
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur2].first);
        ExOp::toMemmove(tree[cur2].first, swapbuf);

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur2].second;
        tree[cur].second = cur2 &0xFFFFFFFE;
        tree[cur2].second =  todel.path[todel.path[0]] | 1;

        todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
        todel.path[todel.path[0]] = cur2;
        todel.path[0]++;
    }

    cur = (todel.path[todel.path[0]] & 1) ? tree[todel.path[todel.path[0]] ^ 1].second | 1 : tree[todel.path[todel.path[0]] ^ 1].second & 0xFFFFFFFE ;
    if ((tree[cur^1].second & 1)&&(tree[cur^1].second  <=firstlone)) { // brother is red
        // one rotation, then done!

        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        ExOp::toMemmove(swapbuf, tree[cur2].first);
        ExOp::toMemmove(tree[cur2].first,tree[todel.path[todel.path[0]] ^ 1].first);
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur].first);
        ExOp::toMemmove(tree[cur].first, swapbuf);

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur].second;
        if (tree[cur].second > firstlone) tree[tree[cur].second].second = todel.path[todel.path[0]] ^ 1;
        tree[cur2].second = (tree[cur2].second & 1) | (cur & 0xFFFFFFFE);
        tree[cur].second =  todel.path[todel.path[0]] & 0xFFFFFFFE;

        tree[cur^1].second &= 0xFFFFFFFE;

        break;
    }else if ((tree[cur].second & 1)&&(tree[cur].second  <=firstlone)){
        // double rotation, then done!
        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        ExOp::toMemmove(swapbuf, tree[cur2].first);
        ExOp::toMemmove(tree[cur2].first,tree[cur].first);
        cur3 = (todel.path[todel.path[0]] & 1) ? tree[cur].second | 1 : tree[cur].second & 0xFFFFFFFE ;
        ExOp::toMemmove(tree[cur].first, tree[cur3 ^ 1].first );
        ExOp::toMemmove(tree[cur3 ^ 1].first, tree[todel.path[todel.path[0]] ^ 1].first );
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur3].first);
        ExOp::toMemmove(tree[cur3].first, swapbuf);

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur3].second;
        if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = todel.path[todel.path[0]] ^ 1;
        tree[cur].second = tree[cur3^1].second;
        if (tree[cur3^1].second > firstlone) tree[tree[cur3^1].second].second = cur;

        tree[cur3].second = todel.path[todel.path[0]] & 0xFFFFFFFE;
        tree[cur3^1].second = cur & 0xFFFFFFFE;
        tree[cur2].second = (tree[cur2].second & 1) | (cur3 & 0xFFFFFFFE);



        break;
        }

    // no rotation... not finished D:
    //	cur3 = tree[path[path[0]]].first;
    //	tree[path[path[0]]].first = -1;
    //	show();
    //	tree[path[path[0]]].first = cur3;

        tree[todel.path[todel.path[0]] ^ 1].second |= 1; // paint brother red
        todel.path[0]--;

        if (todel.path[0] == 0) break;
        if (tree[todel.path[todel.path[0]]].second & 1){
            tree[todel.path[todel.path[0]]].second &= 0xFFFFFFFE;
            break;
            }
    }

}

LFHTEMP	bool RBTofDoom<Key,C>::areLonersHealty()const{
    for(unsigned int i=firstlone; (i++ >> alloc_mag)==0;){
        if (tree[i].second > firstlone) return false;
        if (tree[tree[i].second].second != i) return false;
    }
    return true;
}


LFHTEMP	void RBTofDoom<Key,C>::delete_nodes_routine(pair<Vector<unsigned int>, unsigned int> &list_and_lone){ // listed nodes are assumed to be all the unreachable ones
    unsigned int cur;
    unsigned int cur2,cur3;
    unsigned int i = size - list_and_lone.first.getSize() * 2 + list_and_lone.second;
  //      printf("swapdel start! %c \n", i < (((unsigned int) 1) << (alloc_mag-2)) ? 'Y' : 'N'); fflush(stdout);
    if (i < (((unsigned int) 1) << (alloc_mag-2)) ){ // Downsize, complete rebalance
        if (i == 0) {this->toMemfree(); return;}
        Vector<Key> batch;
        batch.setSize(i);

        cur2 = inorderFirst(path);
        batch[0] = tree[cur2].first;
        for(i=1;i< batch.getSize();i++) { batch[i] = tree[inorderNext(path)].first;}
        delete[](tree+1);
        delete[](path);
        batchInit_routine(batch,false);
       // printf("new size %i, %i\n", size, firstlone);
    }else{
        list_and_lone.first.sort();
       // printf("%i is size (%i >= %i right?)\n", size, firstlone, size - (1 << alloc_mag) + firstlone + 1);
       // printf("got sorted! (first lone was %i)\n", firstlone); fflush(stdout);
       // this->show();
        // swapping loners
        cur = list_and_lone.first.getSize() - list_and_lone.second;
        cur2 = list_and_lone.first.getSize();
        while(cur< cur2){
            firstlone++;
            if (list_and_lone.first[cur] == firstlone) { cur++; continue;} // no swapping needed
            // swap with cur2
            cur2--;
            if (firstlone != tree[tree[firstlone].second].second){
                printf(" %i is nbrem (%i cur cur %i)\n", list_and_lone.second,cur  ,  cur2);
                printf("%i != %i  (p=%i) \n",list_and_lone.first[cur2], tree[tree[list_and_lone.first[cur2]].second].second, tree[list_and_lone.first[cur2]].second);
                ExOp::show(list_and_lone.first);

                this->show();
            }
            LFH_ASSERT(firstlone == tree[tree[firstlone].second].second, "corrupted lone leaf\n");

            ExOp::toMemmove(tree[list_and_lone.first[cur2]].first, tree[firstlone].first);
            tree[list_and_lone.first[cur2]].second = tree[firstlone].second;
            tree[tree[firstlone].second].second = list_and_lone.first[cur2];


            }
        //printf("loners done!\n"); fflush(stdout);
        LFH_ASSERT(this->areLonersHealty(), "Unhealthy\n");


        // swapping inner nodes
        cur = 0;
        cur2 = list_and_lone.first.getSize() - list_and_lone.second -1;
        list_and_lone.second = size - (1 << alloc_mag) + firstlone + 1 - list_and_lone.second;
        size = i;
//        printf("dacur %i %i (firstlo = %i)\n", cur, cur2, firstlone); fflush(stdout);
        while(cur < (cur2+1)){
//             printf("dacur %i %i %i\n", cur, cur2,  list_and_lone.second); fflush(stdout);
             list_and_lone.second -= 2;
             if (((list_and_lone.first[cur2] ^ list_and_lone.second) & 0xFFFFFFFE) == 0 ) { cur2--; continue;} // no swapping needed

//            printf("time to find %i!\n", tree[list_and_lone.second].first  ); fflush(stdout);
            typename RBTofDoom<Key,C>::Iterator ite = this->find_first( tree[list_and_lone.second].first );
            if (!ite.isValid()) {printf("Could not reach %i (at %i)!\n", tree[list_and_lone.second].first, list_and_lone.second ); this->show(); LFH_exit(1);}
            for(; (ite.path[0] == 0) || ((ite.path[ite.path[0]] ^ list_and_lone.second) & 0xFFFFFFFE) != 0 ; ){
                ++ite;
                if (!ite.isValid()) {printf("Could not reach %i (at %i)!\n", tree[list_and_lone.second].first, list_and_lone.second ); this->show(); LFH_exit(1);}
            }

            cur3 = (ite.path[0] == 1) ? 1 : ite.path[ite.path[0]-1];
//            printf("replace %i->%i\n", list_and_lone.second & 0xFFFFFFFE , list_and_lone.first[cur] & 0xFFFFFFFE); fflush(stdout);
            tree[cur3].second = (tree[cur3].second & 1) | (list_and_lone.first[cur] & 0xFFFFFFFE);
            cur3 = list_and_lone.first[cur] & 0xFFFFFFFE;

            ExOp::toMemmove(tree[cur3].first, tree[list_and_lone.second].first );
            ExOp::toMemmove(tree[cur3 | 1].first, tree[list_and_lone.second | 1].first );
            tree[cur3].second = tree[list_and_lone.second].second;
            tree[cur3| 1].second = tree[list_and_lone.second| 1].second;


            if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur3;
            if (tree[cur3|1].second > firstlone) tree[tree[cur3|1].second].second = cur3 | 1;

            cur++;
        }
        //printf("pairs done!\n"); fflush(stdout);
        LFH_ASSERT(this->areLonersHealty(), "Unhealthy\n");

    }
    //printf("swapdel end!\n"); fflush(stdout);
}

LFHTEMP	void RBTofDoom<Key,C>::remove_subtree_routine( typename RBTofDoom<Key,C>::Iterator &todel, pair<Vector<unsigned int>, unsigned int> &list_and_lone){
    if (todel.path[0] == 0) {todel.path[0] = 0xFFFFFFFF; this->toMemfree(); return;}
    //printf("remove_subtree_routine called on"); ExOp::show(*todel);
    bool prob;
    unsigned int cur = todel.path[todel.path[0]];
    Key swapbuf;
    unsigned int cur2 = (todel.path[0] == 1) ? 1 : todel.path[todel.path[0]-1];
    unsigned int cur3;
    if (cur > firstlone){// printf("remrout: lone\n"); fflush(stdout);
        --todel;
        if (todel.path[0] == 0) {todel.path[0] = 0xFFFFFFFF; this->toMemfree(); return;}
        list_and_lone.first.push_back(cur); list_and_lone.second++;
        cur = todel.path[todel.path[0]];
        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        tree[cur].second = 0;

        if ((tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // Red brother... double rotate!
            //   B
            // B   R
            //    B B
            //   ? ?
            cur3 = (cur & 1) ? tree[cur^1].second : tree[cur^1].second & 0xFFFFFFFE;

            ExOp::toMemmove(swapbuf, tree[cur2].first);
            ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
            ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
            ExOp::toMemmove(tree[cur3].first, swapbuf);

            tree[cur ^ 1].second = tree[cur3].second;
            if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

            tree[cur2].second = cur3 & 0xFFFFFFFE;
            tree[cur3].second = cur | 1;

            todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
            todel.path[todel.path[0]] = cur3;
            todel.path[0]++;

            cur = todel.path[todel.path[0]];
            cur2 = todel.path[todel.path[0]-1];
        }

        if (tree[cur^1].second > firstlone){ // brother has lone child, use him!
            //   X
            // B   B
            //    R
            cur3= tree[cur^1].second;

            if (cur & 1){
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
            }else{
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
            }
            list_and_lone.first.push_back(cur3); list_and_lone.second++;
            tree[cur ^ 1].second =0;
        }else{
            cur3 = tree[cur^1].second;
            if ( cur3 == 0) {
                prob = (tree[cur2].second & 1) == 0;
                //   X
                // B   B
           //     printf("in here %i %i (fl%i)\n", list_and_lone.first.last(), cur, firstlone);
           //     printf("in here %i %i\n", tree[list_and_lone.first.last()].first, tree[cur].first);
           //     printf("in here %i %i\n", tree[cur2].first, tree[cur^1].first);


                tree[cur2].second = list_and_lone.first.last();
                tree[list_and_lone.first.last()].second = cur2;

                if (cur & 1){
                //    printf("haha\n");
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first, tree[cur^1].first);
                }else{
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur^1].first);
                }


                list_and_lone.first.last() = cur;
                todel.path[0]--;
                if (prob) this->resolve_doubleblack(todel);
           //     this->
;
            } else {
                //   X
                // B   B
                //    R R
                tree[cur^1].second = list_and_lone.first.last();
                tree[list_and_lone.first.last()].second = cur^1;
                ExOp::toMemmove(tree[cur].first, tree[cur2].first);
                if (cur & 1){
                    ExOp::toMemmove(tree[cur2].first, tree[cur3 | 1].first);
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur^1].first);
                    ExOp::toMemmove(tree[cur^1].first, tree[cur3 &0xFFFFFFFE].first);
                }else{
                    ExOp::toMemmove(tree[cur2].first, tree[cur3 &0xFFFFFFFE].first);
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur3 | 1].first);
                }
                list_and_lone.first.last() = cur3;
            }
            list_and_lone.second--;

        //    printf("dblack is at %i\n",  todel.path[0] == 0 ? 1 :  todel.path[todel.path[0]]);

        }
    }else{
        cur2 = (todel.path[0] == 1) ? 1 : todel.path[todel.path[0]-1];
        //if (!delete_parent_as_well) ExOp::toMemmove(swap: tree[cur2]);
        prob = ((tree[cur2].second & 1) == 0)&&((tree[cur^1].second > firstlone)||((tree[cur^1].second & 1) == 0));  // are parent and brother both black?

        ExOp::toMemmove(tree[cur2].first, tree[cur^1].first);
        if (tree[cur^1].second > firstlone){
            tree[cur2].second = tree[cur^1].second;
            tree[tree[cur^1].second].second = cur2;
        }else tree[cur2].second = tree[cur^1].second & 0xFFFFFFFE;
        unsigned int i;
        list_and_lone.first.push_back(cur);
        i=list_and_lone.first.getSize();
        if (tree[cur].second > 1) {
            list_and_lone.first.push_back(tree[cur].second);
            for(;i< list_and_lone.first.getSize();i++){
                if (list_and_lone.first[i] > firstlone) {list_and_lone.second++; continue;}
                if (tree[list_and_lone.first[i]].second > 1) list_and_lone.first.push_back(tree[list_and_lone.first[i]].second);
                if (tree[list_and_lone.first[i]^1].second > 1) list_and_lone.first.push_back(tree[list_and_lone.first[i]^1].second);
            }
        }
        todel.path[0]--;
        //printf("does res %c\n", prob ? 'Y' : 'N');
        if (prob) {
            this->resolve_doubleblack(todel);

        }
    }
}



    LFHTEMP	template<class SORTER> void RBTofDoom<Key,C>::intersection(vector<Key> &f_out, const SORTER &Query) const{

		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = min;
		Key tmax = max;

		if (size ==0) return;
		while( true ){


			while(!(SETCMP_DISJOINT & Query(tmin,tmax))){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if (!(SETCMP_DISJOINT & Query(tree[tree[path[depth]].second].first,tree[tree[path[depth]].second].first))) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
			if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = max;
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;

		}



	}
	LFHTEMP	void RBTofDoom<Key,C>::intersectionInterval(vector<Key> &f_out, const Key &q_min, const Key &q_max) const{
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = min;
		Key tmax = max;

		if (size ==0) return;
		while( true ){


			while((tmin <= q_max)&&(tmax >= q_min)){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if ((q_min <= tree[tree[path[depth]].second].first)&&(q_max >=tree[tree[path[depth]].second].first)) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;

			if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = max;
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;
		}


		}
	LFHTEMP	template<class SORTER> void RBTofDoom<Key,C>::intersection(Vector<Key> &f_out, const SORTER& Query) const{
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = min;
		Key tmax = max;

		if (size ==0) return;
		while( true ){


			while(!(SETCMP_DISJOINT & Query(tmin,tmax))){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if (!(SETCMP_DISJOINT & Query(tree[tree[path[depth]].second].first))) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
			if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = max;
			else tmax = tree[path[depth-1]].first;
			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;
		}
	}
	LFHTEMP	void RBTofDoom<Key,C>::intersectionInterval(Vector<Key> &f_out, const Key &q_min, const Key &q_max) const{
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = min;
		Key tmax = max;


		if (size ==0) return;
		while( true ){


			while((tmin <= q_max)&&(tmax >= q_min)){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if ((q_min <= tree[tree[path[depth]].second].first)&&(q_max >=tree[tree[path[depth]].second].first)) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
			if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = max;
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;

		}
    }


LFHTEMP	void RBTofDoom<Key,C>::removeRange(const Key &q_min, const Key &q_max){
//	vector<Key> fin;
//	intersection(fin, q_min,q_max);
//	for(int i=0 ;i< fin.size();i++) remove(fin[i]);
    while (true){
    typename RBTofDoom<Key,C>::Iterator ite = this->find_first(q_min);
        if (!ite.isValid()) break;
        if ((*ite) > q_max) break;
        this->remove(ite);
    }
}

LFHTEMP	void RBTofDoom<Key,C>::removeRange(const Key &q_min, typename RBTofDoom<Key,C>::Iterator &i_max){
	Key maxkey = (*i_max);
    while((*i_max) >= q_min){
        this->remove(i_max);
		i_max.findLE(maxkey);
		if (!i_max.isValid()) break;
    }
}

LFHTEMP	void RBTofDoom<Key,C>::removeRange(typename RBTofDoom<Key,C>::Iterator &i_min, const Key &q_max){
	Key minkey = (*i_min);
    while ((*i_min) <= q_max){
        this->remove(i_min);
        i_min.finGE(minkey);
        if (!i_min.isValid()) break;
    }
}


LFHTEMP template<class OKEY> void RBTofDoom<Key,C>::removeRange(const OKEY &q_min, typename RBTofDoom<Key,C>::Iterator &i_max){
	Key maxkey = (*i_max);
    while((*i_max) >= q_min){
        this->remove(i_max);
		i_max.findLE(maxkey);
		if (!i_max.isValid()) break;
    }
}

LFHTEMP	template<class OKEY> void RBTofDoom<Key,C>::removeRange(typename RBTofDoom<Key,C>::Iterator &i_min, const OKEY &q_max){
	Key minkey = (*i_min);
    while ((*i_min) <= q_max){
        this->remove(i_min);
        i_min.findGE(minkey);
        if (!i_min.isValid()) break;
    }
}

LFHTEMP	template<class OKEY> void RBTofDoom<Key,C>::range_move_routine(const OKEY& q_min,const OKEY& q_max, Vector<Key> &moveout){

    typename RBTofDoom<Key,C>::Iterator ite(*this, q_min);
    while (true){
        if (!ite.isValid()) break;
        if ((*ite) > q_max) break;
       moveout.mempush_back(*ite);
    }

}

// can be improved...
LFHTEMP	void RBTofDoom<Key,C>::overwriteRange(const Key &q_min, const Key &q_max, const RBTofDoom<Key,C>& replacement){
    if (replacement.getSize() == 0) return removeRange(q_min,q_max);
    if ((replacement.min < q_min)||(replacement.max > q_max)){
        printf("illegal overwriteRange, inserted elements must be within range\n");
        ExOp::show(replacement.min);
        ExOp::show(replacement.max);
        ExOp::show(q_min);
        ExOp::show(q_max);
        LFH_exit(1);
    }

    typename RBTofDoom<Key,C>::Iterator ite = this->find_first(q_min);
    typename RBTofDoom<Key,C>::Iterator ite_in = replacement.first();
    unsigned int flag =0;
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) flag |= 1;
        if (!ite_in.isValid()) flag |= 2;
        if (flag != 0) break;
        ++ite;    ++ite_in;
    }
    Key  cur;
    switch(flag){
        case 1: // has more to insert!
                ite = this->find_first(q_min);
                ite_in = replacement.first();
                while (true){
                    if ((!ite.isValid())||((*ite) > q_max)) break;
                    if (!ite_in.isValid()) break;
                    ite.ordering_preserving_reference() = (*ite_in);
                    ++ite;    ++ite_in;
                }

                do{
                this->insert(*ite_in);
                ++ite_in;
                } while(ite_in.isValid());
            return;
        case 2: // has more to delete!
            cur = (*ite);
            while (true){
                this->remove(ite);
                ite = this->find_first(cur);
                if (!ite.isValid()) break;
                if ((*ite) > q_max) break;
            }
            break;
        default:// exact replacement
            break;
    }
    ite = this->find_first(q_min);
    ite_in = replacement.first();
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) break;
        if (!ite_in.isValid()) break;
        ite.ordering_preserving_reference() = (*ite_in);
        ++ite;    ++ite_in;
    }


}


// can be improved...
LFHTEMP	template<class OKEY> void RBTofDoom<Key,C>::overwriteRange(const OKEY &q_min, const OKEY &q_max, const RBTofDoom<Key,C>& replacement){
    if (replacement.getSize() == 0) return removeRange(q_min,q_max);
    if ((replacement.min < q_min)||(replacement.max > q_max)){
        printf("illegal overwriteRange, inserted elements must be within range\n");
        ExOp::show(replacement.min);
        ExOp::show(replacement.max);
        ExOp::show(q_min);
        ExOp::show(q_max);
        LFH_exit(1);
    }

    typename RBTofDoom<Key,C>::Iterator ite = this->find_first(q_min);
    typename RBTofDoom<Key,C>::Iterator ite_in = replacement.first();
    unsigned int flag =0;
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) flag |= 1;
        if (!ite_in.isValid()) flag |= 2;
        if (flag != 0) break;
        ++ite;    ++ite_in;
    }
    Key  cur;
    switch(flag){
        case 1: // has more to insert!
                ite = this->find_first(q_min);
                ite_in = replacement.first();
                while (true){
                    if ((!ite.isValid())||((*ite) > q_max)) break;
                    if (!ite_in.isValid()) break;
                    ite.ordering_preserving_reference() = (*ite_in);
                    ++ite;    ++ite_in;
                }

                do{
                this->insert(*ite_in);
                ++ite_in;
                } while(ite_in.isValid());
            return;
        case 2: // has more to delete!
            cur = (*ite);
            while (true){
                this->remove(ite);
                ite = this->find_first(cur);
                if (!ite.isValid()) break;
                if ((*ite) > q_max) break;
            }
            break;
        default:// exact replacement
            break;
    }
    ite = this->find_first(q_min);
    ite_in = replacement.first();
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) break;
        if (!ite_in.isValid()) break;
        ite.ordering_preserving_reference() = (*ite_in);
        ++ite;    ++ite_in;
    }


}


LFHTEMP	template<class OKEY> void RBTofDoom<Key,C>::removeRange_faster(const OKEY &q_min, const OKEY &q_max, bool doit){
    // remove red nodes in range (should not rotate at all)
    pair<Vector<unsigned int>, unsigned int> list_and_lone; list_and_lone.second =0;

    typename RBTofDoom<Key,C>::Iterator ite(this);
    unsigned int cur;
    int flag =0;

    if ((max < q_min)||(min > q_max)) return;

    if (max <= q_max) flag |= 1;
    if (min >= q_min) flag |= 2;

    if (flag == 3) {this->toMemfree(); return;}
    unsigned int cflag;
    while(true){
        if (doit) this->show();
        ite.path[0] = 0;
        switch(flag){
            case 0:
            for(cflag=0; cflag == 0;){
                while(ExOp::isLT(tree[ite.get_index()].first,q_min)){
                    if (!ite.hasRightChild()) {delete_nodes_routine(list_and_lone);return;}
                    if ((cur =tree[ite.get_index()].second) > firstlone){
                        if ((ExOp::isGE(tree[cur].first, q_min))&&(ExOp::isLE(tree[cur].first, q_max))) {
                            ite.path[++ite.path[0]] = cur;
                            remove_static_routine(ite,list_and_lone);
                        }
                        delete_nodes_routine(list_and_lone);return;
                    }
                    ite.path[++ite.path[0]] = cur | 1;
                }
                if (ExOp::isLE(tree[ite.get_index()].first,q_max)) {cflag = ite.path[0]; break;}
                while(ExOp::isGT(tree[ite.get_index()].first,q_max)){
                    if (!ite.hasLeftChild()) {delete_nodes_routine(list_and_lone);return;}
                    cur = (tree[ite.get_index()].second & 0xFFFFFFFE); // should exist, since q_min < q_max < max
                    ite.path[++ite.path[0]] = cur;
                }
                if (ExOp::isGE(tree[ite.get_index()].first,q_min)) {cflag = ite.path[0]; break;}
            }
            break;
            case 1:
           //     printf("hahaha\n");
                if (ExOp::isLT(tree[1].first, q_min)) {
                    cur =1;
                    do{
                        if ((cur = tree[cur].second) > firstlone){ // deleting that last lone
                            ite.path[++ite.path[0]] = cur; remove_static_routine(ite,list_and_lone);  delete_nodes_routine(list_and_lone);
                            ite = this->last();
                            max = *ite;
                            return;
                        }
                        cur |=1;
                        if (cur == 1){
                            delete_nodes_routine(list_and_lone);
                            ite = this->last();
                            max = *ite;
                            return;
                        }
                        ite.path[++ite.path[0]] = cur;
                    }while(ExOp::isLT(tree[cur].first, q_min));
                }
                cflag = ite.path[0];
            break;
            case 2:
            //   printf("hehehe\n");
                if (ExOp::isGT(tree[1].first, q_max)) {
                    cur =1;
                    do{
                        cur = tree[cur].second;
                        if ((cur > firstlone)||((cur & 0xFFFFFFFE) == 0)) {
                             delete_nodes_routine(list_and_lone);
                             ite = this->first();
                             min = *ite;
                             return;
                        }
                        cur &= 0xFFFFFFFE;
                      //  printf("%i ttt\n", cur);
                        ite.path[++ite.path[0]] = cur;
                    }while(ExOp::isGT(tree[cur].first, q_max));
                }
                cflag = ite.path[0];
            break;
        }
      //  printf("%i\n", cflag);
      //  printf("highest in is"); ExOp::show(*ite);
        // trying Right
        if (ite.hasRightChild()){
            if ((cur =tree[ite.get_index()].second) > firstlone){
                if (ExOp::isLE(tree[cur].first, q_max)) {ite.path[++ite.path[0]] = cur; remove_subtree_routine(ite, list_and_lone);}
                else remove_static_routine(ite,list_and_lone);
                continue;
            }else{
                ite.path[++ite.path[0]] = cur | 1;
                while(ExOp::isGT(tree[ite.get_index()].first,q_max)){
                    if (!ite.hasLeftChild()) {ite.path[0] = cflag;break;}
                    cur = (tree[ite.get_index()].second & 0xFFFFFFFE); // should exist, since q_min < q_max < max
                    ite.path[++ite.path[0]] = cur;
                }
                if (ite.path[0] != cflag){
                    if (!ite.hasLeftChild()) remove_static_routine(ite,list_and_lone);
                    else {
                        cur = tree[ite.get_index()].second & 0xFFFFFFFE;
                        ite.path[++ite.path[0]] = cur;
                        remove_subtree_routine(ite, list_and_lone);
                    }
                    continue;
                }
               // printf("continue! on "); ExOp::show(*ite); fflush(stdout);
            }
        }

        // trying left
        if (ite.hasLeftChild()){
            cur = (tree[ite.get_index()].second & 0xFFFFFFFE);
            ite.path[++ite.path[0]] = cur;
            while(ExOp::isLT(tree[ite.get_index()].first,q_min)){
                if (!ite.hasRightChild()) {ite.path[0] = cflag;break;}
                if ((cur =tree[ite.get_index()].second) > firstlone){
                    if (ExOp::isGE(tree[cur].first,q_min)) {
                        ite.path[++ite.path[0]] = cur;
                        remove_static_routine(ite,list_and_lone);
                    }
                    ite.path[0] = cflag;break;
                }
                cur = (tree[ite.get_index()].second | 1); // should exist, since q_min < q_max < max
                ite.path[++ite.path[0]] = cur;
            }
            if (ite.path[0] != cflag){
                if (!ite.hasRightChild()) remove_static_routine(ite,list_and_lone);
                else{
                    cur = tree[ite.get_index()].second;
                    ite.path[++ite.path[0]] = cur | (cur > firstlone ? 0 : 1);
                    remove_subtree_routine(ite, list_and_lone);
                }
                continue;
            }
          //   printf("continue!\n"); fflush(stdout);
        }

        //printf("final!\n");    printf("firstlone is  %i\n", firstlone);fflush(stdout);
        remove_static_routine(ite,list_and_lone);
        break;
    }
   // this->show();

    delete_nodes_routine(list_and_lone);

    //printf("recoverytime\n"); fflush(stdout);

    switch(flag){
        case 1:
            ite = this->last();
            max = *ite;
        break;
        case 2:
            ite = this->first();
            min = *ite;
        break;

    }
    return;
}



LFHTEMP	template<class OKEY> void RBTofDoom<Key,C>::removeRange(const OKEY &q_min, const OKEY &q_max){
    while (true){
    typename RBTofDoom<Key,C>::Iterator ite(*this, q_min);
        if (!ite.isValid()) break;
        if (ExOp::isGT((*ite),q_max)) break;

        this->remove(ite);
    }

}

LFHTEMP	template<class OKEY, class FF> void RBTofDoom<Key,C>::removeRange(const OKEY &q_min, const OKEY &q_max, const FF& filter_func){
    OKEY tmpmin = q_min;
    typename RBTofDoom<Key,C>::Iterator ite(*this,q_min);
    while (true){
        if (!ite.isValid()) break;
        if (ExOp::isGT((*ite),q_max)) break;
        if (FF(*ite)) {this->remove(ite); ite = this->find_first(tmpmin);}
        else {tmpmin = q_min; ++ite;}
    }
}

LFHTEMP	unsigned int RBTofDoom<Key,C>::inorderFirst(unsigned int *path) const{
    unsigned int cur = 1;
    path[0] =0;
    while(tree[cur].second > 1){
        if (tree[cur].second > firstlone) break;
        else cur = tree[cur].second & 0xFFFFFFFE;
        path[++(path[0])] =cur;
    }
    return(cur);
}

	LFHTEMP	unsigned int RBTofDoom<Key,C>::inorderNext(unsigned int *path) const{

		unsigned int cur;
		bool branch;
		if (path[0] == 0) {
			if (tree[1].second == 0) {return(0xFFFFFFFF);}
			branch = false;
		} else{
			branch = ((tree[path[path[0]]].second < 2)||(path[path[0]] > firstlone))  ;
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no right child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			if ((path[0]>0)&&((tree[path[path[0]]].second > firstlone)||((tree[path[path[0]]].second | 1) == path[path[0]+1]))) {path[0]--;
				while((path[0]>0)&&((tree[path[path[0]]].second | 1) == path[path[0]+1])) path[0]--;
			}

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
				return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}else return(path[path[0]]);

		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (tree[cur].second > firstlone) { cur = tree[cur].second; path[++(path[0])] =cur;}
			else {cur = tree[cur].second | 1;
				path[0]++;
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);
				path[path[0]] =cur;
				while(tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (tree[cur].second > firstlone) break;
					else cur = tree[cur].second & 0xFFFFFFFE;
					path[++(path[0])] =cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
			}
		}

		return(cur);
	}
	LFHTEMP	bool RBTofDoom<Key,C>::isRed(unsigned int w) const{

		return((w > firstlone)||((tree[w].second & 1)&&(tree[w].second <= firstlone)));

	}
    LFHTEMP	bool RBTofDoom<Key,C>::isLeaf(unsigned int w) const{

			return((w > firstlone)||(tree[w].second< 2));
		}
	LFHTEMP	void RBTofDoom<Key,C>::show(FILE* f, int level) const{
		unsigned int buffer[1024];
		char funbuf[256];
		memset(funbuf,' ', sizeof(char) * 256);
		funbuf[255] = '\0';
		int loop =0;
		if (size == 0) {fprintf(f,"Empty RBTree!\n"); return;}
		fprintf(f,"RBTree with %i nodes on interval [", size);
		ExOp::show(min,f,1);
		fprintf(f," ; ");
		ExOp::show(max,f,1);
		fprintf(f,"] allocsize = %i\n", 1 << alloc_mag);
		int fun;

		for(unsigned int cur = inorderFirst(buffer); cur != 0xFFFFFFFF;cur = inorderNext(buffer)){
			if (cur > firstlone) {
				fun = 0;
				for (loop = buffer[0]-1;loop > 0;loop--) if (!isRed(buffer[loop])) fun++;
				fprintf(f,"%i%sR{%i}%c",fun , funbuf + 255 - buffer[0]+1,cur, (tree[tree[cur].second].second != cur) ? 'E' : ' ');

			}else fprintf(f,"%s%c{%i}  ", funbuf + 255 - buffer[0], isRed(cur) ? 'R' : 'B',cur);
			ExOp::show(tree[cur].first,f,1);
			fprintf(f,"\n");
		}
	}
LFHTEMP	ERRCODE RBTofDoom<Key,C>::save(FILE* f) const{
    //saves an ordered array!
    unsigned int buffer[1024];
    fwrite(&size ,sizeof(unsigned int),1 ,f);
    if (size == 0) return 0;

    for(unsigned int cur = inorderFirst(buffer); cur != 0xFFFFFFFF;cur = inorderNext(buffer)){
        ExOp::save(tree[cur].first, f);
    }
    return 0;
}
LFHTEMP	ERRCODE RBTofDoom<Key,C>::load(FILE* f){
	unsigned int c,i_size;
	this->toMemfree();
	if (1 != fread(&i_size ,sizeof(unsigned int),1 ,f)) return 1;

	if (i_size == 0) return 0;
	this->makeEmptyTree(i_size);

	// fixed known size!
	unsigned int buffer[1024];
	unsigned int cur = inorderFirst(buffer);
	ExOp::load(tree[cur].first, f);
	min = tree[cur].first;c = i_size -1;
	if (i_size > 1) while((c--) != 0){
		cur = inorderNext(buffer);
		ExOp::load(tree[cur].first, f);
	}
	max = tree[cur].first;
	return 0;
}
LFHTEMP	typename RBTofDoom<Key,C>::Iterator RBTofDoom<Key,C>::find(const Key& tval) const{ // find the smallest which is >= key
    typename RBTofDoom<Key,C>::Iterator fout = find_first(tval);
    if ((fout.isValid())&&(*fout != tval)) fout.path[0] = 0xFFFFFFFF;
    return(fout);
    }



LFHTEMP void RBTofDoom<Key,C>::test() const{
    printf("testing... ");
    typename RBTofDoom<Key,C>::Iterator itemin = this->first();
    unsigned int nbitem;
    for(nbitem=0;itemin.isValid();++itemin){
        typename RBTofDoom<Key,C>::Iterator itemax = this->find(*itemin);
        nbitem++;

   //     for(++itemax;itemax.isValid();++itemax){
            if ((*itemin) > (*itemax)) {printf("failed!\n"); return;}
   //     }

    }
    if (nbitem != size) printf("failed!\n");
    else printf("passed!\n");
}



LFHTEMP template<class F> bool RBTofDoom<Key,C>::first(RBTofDoom<Key,C>::QueryIterator& ite, const F &f)const{
    SETCMP_enum sc;
    if (size < 2){
        if (size == 0) return false;
        ite.path[0] =0;
        return f(tree[1].first);
    }else if (f(min,max) & SETCMP_DISJOINT) return false;
    unsigned int cur=1;
    ite.path[0] =0;
    ite.iaindex[0] =0; ite.iaindex[1] =32;
    ite.iaval[0] = min; ite.iaval[32] = max;
    while(tree[cur].second > 1){
        if (tree[cur].second > firstlone) {
            // lone children
            if ((f(tree[cur].first)& SETCMP_DISJOINT) == 0) return true;
            if ((f(tree[tree[cur].second].first)& SETCMP_DISJOINT) == 0){
                path[++(path[0])] =tree[cur].second;
                return true;
            }
            break;
            // needs to back up!

        } else if (((f(ite.iaval[ite.iaindex[0]] ,tree[cur].first)) & SETCMP_DISJOINT) == 0){
            ite.iaindex[1]++; ite.iaval[ite.iaindex[1]] = tree[cur].first;
            cur = tree[cur].second & 0xFFFFFFFE;
            path[++(path[0])] =cur;
        }else if (((f(tree[cur].first, ite.iaval[ite.iaindex[1]] )) & SETCMP_DISJOINT) == 0){
            ite.iaindex[0]++; ite.iaval[ite.iaindex[0]] = tree[cur].first;
            cur = tree[cur].second | 1;
            path[++(path[0])] =cur;
        }else break;// else needs to back up!, but *should* not happen
    }
    return false;
}

LFHTEMP template<class F> bool RBTofDoom<Key,C>::next(RBTofDoom<Key,C>::QueryIterator&, const F&)const{
    return true;
}


LFHTEMP	Key& RBTofDoom<Key,C>::orderpreserve_deref(const typename RBTofDoom<Key,C>::Iterator& ite){return tree[ite.get_index()].first;}



LFHTEMP	RBTofDoom<Key,C>::Iterator::Iterator(const RBTofDoom<Key,C>& i_target) : target(i_target){}
LFHTEMP	RBTofDoom<Key,C>::Iterator::Iterator(const RBTofDoom<Key,C>& _target, const Key& key): target(_target){this->findGE(key);}
LFHTEMP	template <class O> RBTofDoom<Key,C>::Iterator::Iterator(const RBTofDoom<Key,C>& _target, const O& key): target(_target){this->findGE(key);}
LFHTEMP	RBTofDoom<Key,C>::Iterator::Iterator(const RBTofDoom<Key,C>& _target, const Key& key, SETCMP_enum q): target(_target){
	switch(q){
		case SETCMP_GE: this->findGE(key); break;
		case SETCMP_GT: this->findGT(key); break;
		case SETCMP_LE: this->findLE(key); break;
		case SETCMP_LT: this->findLT(key); break;
		default: this->find(key); break;
	}
}

LFHTEMP	template <class O> RBTofDoom<Key,C>::Iterator::Iterator(const RBTofDoom<Key,C>& _target, const O& key, SETCMP_enum q): target(_target){
	switch(q){
		case SETCMP_GE: this->findGE(key); break;
		case SETCMP_GT: this->findGT(key); break;
		case SETCMP_LE: this->findLE(key); break;
		case SETCMP_LT: this->findLT(key); break;
		default: this->find(key); break;
	}
}

LFHTEMP	const typename RBTofDoom<Key,C>::Iterator& RBTofDoom<Key,C>::Iterator::operator++(){
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return(*this);}
			branch = false;
		} else{
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no right child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			if ((path[0]>0)&&((target.tree[path[path[0]]].second > target.firstlone)||((target.tree[path[path[0]]].second | 1) == path[path[0]+1]))) {path[0]--;
				while((path[0]>0)&&((target.tree[path[path[0]]].second | 1) == path[path[0]+1])) path[0]--;
			}

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
            if (((target.tree[1].second & 0xFFFFFFFE)!= path[1])||(target.tree[1].second > target.firstlone)) path[0] = ExCo<unsigned int>::mkMaximum();
		//	return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}

		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (target.tree[cur].second > target.firstlone) { cur = target.tree[cur].second; path[++(path[0])] =cur;}
			else {cur = target.tree[cur].second | 1;
				path[0]++;

				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);

				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);
				path[path[0]] =cur;




				while(target.tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (target.tree[cur].second > target.firstlone) break;
					else cur = target.tree[cur].second & 0xFFFFFFFE;
					path[++(path[0])] =cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
			}
		}
		return(*this);
}
LFHTEMP	const typename RBTofDoom<Key,C>::Iterator& RBTofDoom<Key,C>::Iterator::operator--(){
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) {path[0] = ExCo<unsigned int>::mkMaximum(); return(*this);}
			branch = false;
		} else if (path[path[0]] > target.firstlone) {
		    path[0]--; return(*this);
		}else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no left child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			while((path[0]>0)&&((target.tree[path[path[0]]].second & 0xFFFFFFFE) == path[path[0]+1])) path[0]--;


			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
            if ((target.tree[1].second | 1) != path[1]) path[0] = ExCo<unsigned int>::mkMaximum();
		//	return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}

		}else{
			// has left
			cur = (path[0] ==0) ? 1 : path[path[0]];
			cur = target.tree[cur].second & 0xFFFFFFFE; path[++(path[0])] =cur;


				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);





				while(target.tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (target.tree[cur].second > target.firstlone) {
					    cur = target.tree[cur].second; path[++(path[0])] =cur;
					    break;
					    }
					cur = target.tree[cur].second | 1;
					path[++(path[0])] = cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);

		}
return(*this);
}


LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findFirst(){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur = 1;
	path[0] =0;
	//  printf("%i.l = %i\n",cur, tree[cur].second); fflush(stdout);
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) break;
		else cur = target.tree[cur].second & 0xFFFFFFFE;
		path[++(path[0])] =cur;
	}
	return true;
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findLast(){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur = 1;
	path[0] =0;
	//  printf("%i.l = %i\n",cur, tree[cur].second); fflush(stdout);
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			path[++(path[0])] =target.tree[cur].second;
			break;
		}else cur = target.tree[cur].second | 1;
		path[++(path[0])] =cur;
	}
	return true;
}


LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findGE(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLT(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLT(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLT(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findLE(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGT(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGT(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGT(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findGT(const Key& key){
if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLE(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLE(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLE(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findLT(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGE(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGE(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGE(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::findEQ(const Key& key){*this = target.find_first(key); return this->isValid();}

LFHTEMP	template<class O> bool RBTofDoom<Key,C>::Iterator::findEQ(const O& key){*this = target.find_first(key); return this->isValid();}

LFHTEMP	template<class O> bool RBTofDoom<Key,C>::Iterator::findGE(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLT(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLT(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLT(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	template<class O> bool RBTofDoom<Key,C>::Iterator::findLE(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGT(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGT(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGT(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}


LFHTEMP	template<class O> bool RBTofDoom<Key,C>::Iterator::findGT(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLE(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLE(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLE(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	template<class O> bool RBTofDoom<Key,C>::Iterator::findLT(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	unsigned int height =0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGE(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGE(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGE(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}


/*
LFHTEMP	void RBTofDoom<Key,C>::Iterator::toMaxWithinBrotherSubtree_routine(unsigned int height){
	if height == 0
}
LFHTEMP	void RBTofDoom<Key,C>::Iterator::toMinWithinBrotherSubtree_routine(unsigned int height){

}
*/

//prev LT
//next GT


LFHTEMP	typename RBTofDoom<Key,C>::Iterator& RBTofDoom<Key,C>::Iterator::operator=(const typename RBTofDoom<Key,C>::Iterator& other){
	target = other.target;
	memcpy(path,other.path, sizeof(int) *63);
return *this;}

LFHTEMP	bool RBTofDoom<Key,C>::Iterator::hsNext() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(false); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) return(false);
			branch = false;
		} else{
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
		}
		if (branch) {
			// no right child!
			cur = path[0]-1;
			if ((cur>0)&&((target.tree[path[cur]].second > target.firstlone)||((target.tree[path[cur]].second | 1) == path[cur+1]))) {cur--;
				while((cur>0)&&((target.tree[path[cur]].second | 1) == path[cur+1])) cur--;
			}

			if (cur==0) {
            if (((target.tree[1].second & 0xFFFFFFFE)!= path[1])||(target.tree[1].second > target.firstlone)) return false;
			}

		}
		return(true);
}
LFHTEMP	bool RBTofDoom<Key,C>::Iterator::hsPrev() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(false); // invalid!
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) return(false);
			branch = false;
		} else if (path[path[0]] > target.firstlone) {return true;
		} else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
		}
		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			if (cur==0) {
                if ((target.tree[1].second | 1) != path[1]) return(false);
			}
		}
    return(true);
}

LFHTEMP	const Key& RBTofDoom<Key,C>::Iterator::Next() const{
    unsigned int j = Next_index();
    return target.tree[(j == 0) ? 1: j].first;
    }

LFHTEMP	unsigned int RBTofDoom<Key,C>::Iterator::Next_index() const{
        unsigned int cur;
		bool branch;
	//	if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) return 0;
			branch = false;
		} else {
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
		}

		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			if (cur==0) {
                if ((target.tree[1].second | 1) != path[1]) return(false);
			}
		}

		if (branch) {
			cur = path[0]-1;
			if ((cur>0)&&((target.tree[path[cur]].second > target.firstlone)||((target.tree[path[cur]].second | 1) == path[cur+1]))) {cur--;
				while((cur>0)&&((target.tree[path[cur]].second | 1) == path[cur+1])) cur--;
			}
			return (cur==0) ? 0 : path[cur];
		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (target.tree[cur].second > target.firstlone) { cur = target.tree[cur].second;}
			else {cur = target.tree[cur].second | 1;
				while(target.tree[cur].second > 1){
					if (target.tree[cur].second > target.firstlone) break;
					else cur = target.tree[cur].second & 0xFFFFFFFE;
				}
			}
		}
	return cur;
}

LFHTEMP	const Key& RBTofDoom<Key,C>::Iterator::Prev() const{
    unsigned int j = Prev_index();
    return target.tree[(j == 0) ? 1: j].first;
    }

LFHTEMP	unsigned int RBTofDoom<Key,C>::Iterator::Prev_index() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return 0;
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) return 0;
			branch = false;
		} else if (path[path[0]] > target.firstlone) {
            return( path[0] == 1 ? 1 : path[path[0]-1] );
		}else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
		}
    //    printf("%c brbr\n", branch ? 'Y' : 'N');
		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			return (cur==0) ? (((target.tree[1].second | 1) == path[1]) ? 1 : 0 )  : path[cur] ;
		}else{
			// has left
			cur = (path[0] ==0) ? 1 : path[path[0]];
			cur = target.tree[cur].second & 0xFFFFFFFE;
				while(target.tree[cur].second > 1){
					if (target.tree[cur].second > target.firstlone) {
					    cur = target.tree[cur].second;
					    break;
				    }
				cur = target.tree[cur].second | 1;
			}
		}
	return(cur);
}


LFHTEMP unsigned int RBTofDoom<Key,C>::Iterator::get_index() const{return (path[0]) ? path[path[0]] : 1 ; }
LFHTEMP bool RBTofDoom<Key,C>::Iterator::isValid() const{return ((path[0] &0xFFFF0000) == 0);}
LFHTEMP typename RBTofDoom<Key,C>::Iterator& RBTofDoom<Key,C>::Iterator::toInvalid(){path[0] = 0xFFFFFFFF; return *this;}

LFHTEMP const Key& RBTofDoom<Key,C>::Iterator::operator*() const{ return  target.tree[(path[0]) ? path[path[0]] : 1 ].first;  }
LFHTEMP const Key* RBTofDoom<Key,C>::Iterator::operator->() const{return  &(target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP Key* RBTofDoom<Key,C>::Iterator::ordering_preserving_change() const{return  &(target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP Key& RBTofDoom<Key,C>::Iterator::ordering_preserving_reference() const{return  (target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP bool RBTofDoom<Key,C>::Iterator::operator==(const RBTofDoom<Key,C>::Iterator& other)const{
    if ((path[0] &0xFFFF0000) != 0) return false;
    unsigned int i;
    for(i=0;i<path[0];i++) if (path[i] != other.path[i]) return false;
    return true;
}



#undef LFHTEMP
#define LFHTEMP template<class Key>

LFHTEMP RBTofDoom<Key,void>::RBTofDoom(const RBTofDoom<Key,void>& other): size(other.getSize()){
    makeEmptyTree(other.getSize());
    unsigned int buffer[256];
    unsigned int obuffer[256];
    min_offset = inorderFirst(buffer);
    tree[min_offset].first = other.tree[other.inorderFirst(obuffer)].first;
    unsigned int i;
    for(i=1;i < size-1;i++) tree[inorderNext(buffer)].first = other.tree[other.inorderNext(obuffer)].first;
    if (size == 1) max_offset = min_offset;
    else {max_offset = inorderNext(buffer);  tree[max_offset].first = other.tree[other.inorderNext(obuffer)].first;}
}
LFHTEMP RBTofDoom<Key,void>& RBTofDoom<Key,void>::operator=(const RBTofDoom<Key,void>& other){
    this->toMemfree();
    makeEmptyTree(other.getSize());
    if (size == 0) return 0;
    unsigned int buffer[256];
    unsigned int obuffer[256];
    min_offset = inorderFirst(buffer);
    tree[min_offset].first = other.tree[other.inorderFirst(obuffer)].first;
    unsigned int i;
    for(i=1;i < size-1;i++) tree[inorderNext(buffer)].first = other.tree[other.inorderNext(obuffer)].first;
    if (size == 1) max_offset = min_offset;
    else {max_offset = inorderNext(buffer);  tree[max_offset].first = other.tree[other.inorderNext(obuffer)].first;}
return *this;}

LFHTEMP void RBTofDoom<Key,void>::makeEmptyTree(unsigned int nbelem){ // initialize a most balanced tree, that is to be filled!
		size = nbelem;
		alloc_mag = ExOp::upperbound_pow_of_2(nbelem+2);
		if (alloc_mag < 3){

			switch(nbelem){
				case 0: tree = NULL; path =NULL;
				break;
				case 1:
					tree = new pair<Key, unsigned int>[4];  tree--;
					path = new unsigned int[10];
					tree[1].second = 0;
                    firstlone =4;
				break;
				case 2:
					tree = new pair<Key, unsigned int>[4]; tree--;
					path = new unsigned int[10];
					tree[1].second = 4;
					tree[4].second = 1;
                    firstlone =3;
				break;
				}

			}else{

		//       B
		//   B       B
		// R   R   R   R
		//B B B B B B B B
		// ????

		tree = new pair<Key, unsigned int>[1 << alloc_mag];
        if (tree == NULL) {printf("RBTofDoom: could not allocate!\n"); LFH_exit(1);}
		path = new unsigned int[alloc_mag*2+6];
        if (path == NULL) {printf("RBTofDoom: could not allocate!\n"); LFH_exit(1);}
		firstlone = 1 << alloc_mag;
		tree--;
		unsigned int c=1;
		unsigned int m=2;
		unsigned int r,y;

	//	unsigned int cum;
		unsigned int inc;

		for(unsigned char h=0;h<alloc_mag-2;h++){
			r = ((c != 1)&&( ((int)h) % 3 == (((int)alloc_mag) % 3)  )) ? 1 : 0;
			for(;c<m;c++) tree[c].second = (c<<1) | r  ;
			m <<=1;
			}



		unsigned int nbadd = 1 + nbelem - (1 << (alloc_mag-1)); //printf("ndadd = %i\n", nbadd);

		inc =0;
		if (nbadd <= (((unsigned int)1) << (alloc_mag-2))){
			// only loners!
		//	toins_double = 4 << (alloc_mag-3);
            //    printf("hello!\n");
				for(;c<m;c++) {
					inc += nbadd;
			//		printf("%i : %i %i tfund\n",c, inc >> (alloc_mag-2), inc);
					if (inc >> (alloc_mag-2)){
					tree[c].second = firstlone;
					tree[firstlone].second = c;
					firstlone--;
                    inc -= 1 << (alloc_mag-2);
					} else tree[c].second = 0;
					}

			}else{
			nbadd -= (1 << (alloc_mag-2));
			y = (4 << (alloc_mag-3));
		//	toins_double = y+ 2 * nbadd;
				for(;c<m;c++) {
					inc += nbadd;
				//	printf("%i : %i %i fun\n",c, inc >> (alloc_mag-2), inc);
					if (inc >> (alloc_mag-2)){
					tree[c].second = y;
					tree[y].second = 1;
					tree[y+1].second = 1;
					y+= 2;
                    inc -= 1 << (alloc_mag-2);
					} else{
					tree[c].second = firstlone;
					tree[firstlone].second = c;
					firstlone--;
					}
					}

		}
    }
}


LFHTEMP	void RBTofDoom<Key,void>::batchInit(Vector<Key> &batch, bool needsort){
    this->toMemfree();
    this->batchInit_routine(batch, needsort);
}

LFHTEMP	void RBTofDoom<Key,void>::batchInit_routine(Vector<Key> &batch, bool needsort){
    if (needsort) batch.sort();
    this->makeEmptyTree(batch.size());
    unsigned int buffer[1024];

    min_offset = inorderFirst(buffer);
    tree[min_offset].first = batch[0];

    unsigned int i;
    for(i=1;i < size-1;i++) tree[inorderNext(buffer)].first = batch[i];
    if (batch.size() > 1){
        max_offset = inorderNext(buffer);
        tree[max_offset].first = batch[size-1];
    }else max_offset = min_offset;
}

LFHTEMP	bool RBTofDoom<Key,void>::reach(const Key &q){ // stops on equality
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (q != tree[cur].first) {
					if (tree[tree[cur].second].first == q){
						path[++(path[0])] =tree[cur].second;
						return(true);
					} else return false;
					}
					return(true);

				}
				if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
				else if (q != tree[cur].first) cur = tree[cur].second | 1;
				else return true;
				path[++(path[0])] =cur;
			}
			return (q == tree[cur].first);
		}
		return false;
	}
// loop throught the remaining equal nodes, (ignore the highest in tree, it was already found)
LFHTEMP	bool RBTofDoom<Key,void>::reach_continue(const Key &q){ // stops on equality
			// go right!
		/*	for(cur = path[0] ; cur >0 ;cur--) if (tree[path[cur]].first == q) break;
			if ((cur == 0)&&(tree[0].first != q)){
				//reach minimum equal in tree

			}else{

			}*/
			unsigned int cur = path[0] == 0 ? 1 : path[path[0]];
			if (tree[cur].second < 2) return false;
			cur = tree[cur].second | 1;
			path[++(path[0])] =cur;
				while(tree[cur].second > 1){
					if (tree[cur].second > firstlone) {

						if (q != tree[cur].first) {
							if (tree[tree[cur].second].first == q){
								path[++(path[0])] =tree[cur].second;
								return(true);
							} else return false;
						}
						return(true);

					}
					if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
						else if (q != tree[cur].first) cur = tree[cur].second | 1;
							else return true;
					path[++(path[0])] =cur;
				}
				return (q == tree[cur].first);

			return false;
		}
LFHTEMP	void RBTofDoom<Key,void>::reach_par(const Key &q){ // ignores equalities!
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) break;
				if (q < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
				else cur = tree[cur].second | 1;
				path[++(path[0])] =cur;
			}
		} else path[0] = 0;
	}
LFHTEMP	template<class A> bool RBTofDoom<Key,void>::reach(const A &q){ // stops on equality
		path[0] =0;
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (q != tree[cur].first) {
					if (ExOp::isEQ(tree[tree[cur].second].first,q)){
						path[++(path[0])] =tree[cur].second;
						return(true);
					} else return false;
					}
					return(true);

				}
				if (ExOp::isLT(q, tree[cur].first)) cur = tree[cur].second & 0xFFFFFFFE;
				else if (ExOp::isNQ(q, tree[cur].first)) cur = tree[cur].second | 1;
				else return true;
				path[++(path[0])] =cur;
			}
			return (ExOp::isEQ(q,tree[cur].first));
		}
		return false;
	}
LFHTEMP	unsigned int RBTofDoom<Key,void>::find_index(const Key &what) const{ // 0 for not found
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (what != tree[cur].first) {
						if (tree[tree[cur].second].first == what){
							return tree[cur].second;
						} else return 0;
					}
					return cur;
				}
				if (what < tree[cur].first) cur = tree[cur].second & 0xFFFFFFFE;
				else if (what != tree[cur].first) cur = tree[cur].second | 1;
				else return cur;
			}

			return (what == tree[cur].first) ? cur : 0;
		}
		return 0;
		}
LFHTEMP	template<class O_KEY> unsigned int RBTofDoom<Key,void>::find_index(const O_KEY &what) const{ // 0 for not found
		if (size > 0){
			unsigned int cur=1;
			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

					if (tree[cur].first == what) return cur;
					if (tree[tree[cur].second].first == what){
						return tree[cur].second;
					} else return 0;


				}
				if (tree[cur].first > what) cur = tree[cur].second & 0xFFFFFFFE;
				else if (tree[cur].first != what) cur = tree[cur].second | 1;
				else return cur;
			}

			return (tree[cur].first == what) ? cur : 0;
		}
		return 0;
		}
	// regarless if "what" exists, find the previous and next non-equal element to it
LFHTEMP	void RBTofDoom<Key,void>::findPN(const Key &what, unsigned int &prev, unsigned int &next)const{
			prev = next = 0;
			if (size > 0){
			unsigned int cur=1;

			while(tree[cur].second > 1){
				if (tree[cur].second > firstlone) {

						if (what < tree[cur].first) next = cur;
						else{
							if (what != tree[cur].first) prev = cur;
							if (what < tree[tree[cur].second].first) next =tree[cur].second;
							else if (tree[tree[cur].second].first != what) prev = tree[cur].second;
						}
					return;
				}

				if (what < tree[cur].first) {next = cur;cur = tree[cur].second & 0xFFFFFFFE;  }
				else if (what != tree[cur].first) {prev = cur;cur = tree[cur].second | 1;  }
				else { // found item... need to look into both subtrees
					// todo
					next = tree[cur].second | 1;
				//	while(tree[next].second > 1){

				//		}

					prev = tree[cur].second & 0xFFFFFFFE;
					// todo
				//	while(tree[prev].second > 1){

				//		}
					return;
					}
			}

			if (what < tree[cur].first) {next = cur;cur = tree[cur].second & 0xFFFFFFFE;  }
			else if (what != tree[cur].first) {prev = cur;cur = tree[cur].second | 1;  }


		}


		}
LFHTEMP	void RBTofDoom<Key,void>::insert(const Key &newf){
		unsigned int i;
		bool cached;
		Key swapk;
		if (size > 0){
			//if (newf < min) min = newf ;
			//else if (newf > max)  max = newf;

			if (size+2 == (1u << alloc_mag) ){

				Vector<Key> batch;
				batch.setSize(size+1);

				//	show(stdout);
				//		printf("inserting %i, %i\n", newf, tree[inorderFirst(path)].first);
				i =0;
				batch[i] = tree[inorderFirst(path)].first;
				while(true){
					if (newf < batch[i]){
						batch[i+1] = batch[i];
						batch[i] = newf;
						i++;
						break;
					}
					i++;
					if (i == size){
						batch[size] = newf;
						break;
					}
					batch[i] = tree[inorderNext(path)].first;
				}
				for(i++;i<= size;i++) { batch[i] = tree[inorderNext(path)].first;}

				delete[](tree+1);
				delete[](path);
//printf("reloaded updated AAA!\n");
				batchInit_routine(batch,false);
			}else{
				reach_par(newf);
				unsigned int cur = path[0] == 0 ? 1 : path[path[0]];
				unsigned int cur2,cur3;
				size++;
				if (tree[cur].second >= firstlone){

					cur2 = tree[cur].second;
					cur3 = size - (1 << alloc_mag) + firstlone;

					tree[tree[cur2].second].second = cur3; // parent should be black

					LFH_ASSERT(cur3 +1 < firstlone, "Corrupted RBTree!\n");
                    LFH_ASSERT((cur3 & 1) == 0, "Corrupted RBTree!RB\n");
					// get next free slot!

					tree[cur3].second = 1; // red leaf
					tree[cur3+1].first = tree[cur2].first;
					tree[cur3+1].second = 1; // red leaf
					if (max_offset == cur2) max_offset = cur3+1;
					if (min_offset == cur) min_offset = cur3;

					firstlone++;
					if (cur2 != firstlone) {
						tree[cur2] = tree[firstlone];
						tree[tree[cur2].second].second = cur2;
						if (max_offset == firstlone) max_offset = cur2;
					}

					if (newf < tree[cur].first){
						tree[cur3].first = newf;

					}else{
						tree[cur3].first =  tree[cur].first;
						if (newf < tree[cur3+1].first){
							tree[cur].first = newf;

						}else{
							tree[cur].first = tree[cur3+1].first;
							tree[cur3+1].first = newf;
						}
					}

				}else{ // add loner

					//	printf("helloB!");
					cached = (tree[cur].second == 1); // is node red?
					tree[firstlone].second = cur;
					tree[cur].second = firstlone;
                    if (cur == max_offset) max_offset = firstlone;

					if (newf < tree[cur].first){
						tree[firstlone].first = tree[cur].first;
						tree[cur].first = newf;
					}else tree[firstlone].first = newf;
					firstlone--;
					if (cached){

						tree[cur ^ 1].second = 0; // paint the other children in black;
						path[0]--;
						if (path[0] >0){
							tree[path[path[0]]].second |=1;
							path[0]--;
							while((path[0]>0)&&(tree[path[path[0]]].second & 1)){
								if ((tree[path[path[0]]^ 1].second & 1)&&(tree[path[path[0]]^ 1].second <=firstlone)){
									tree[path[path[0]]^ 1].second &= 0xFFFFFFFE;
									tree[path[path[0]]].second &= 0xFFFFFFFE;

									path[0]--;
									if (path[0] !=0) {tree[path[path[0]]].second |= 1; path[0]--;}
								}else{
									if ((path[path[0]+1] & 1)^(path[path[0]] & 1)){ // double rotation
										// 		printf("double rot!\n");show(stdout);
										if ((tree[path[path[0]+1]].second < 2)||(tree[path[path[0]+1]].second > firstlone)){
											printf("OMFG!!!!!\n");

										}else{

											path[0]--;
									//			printf("here wego !\n");
											if (path[0] == 0) {
												swapk = tree[1].first;
												tree[1].first = tree[path[2]].first;
												tree[1].second = tree[path[2]].second & 0xFFFFFFFE;
											}else {
												swapk = tree[path[path[0]]].first;
												tree[path[path[0]]]= tree[path[path[0]+2]];
											}

											//		printf("%c\n", (path[path[0]+2] & 1) ? 'Y' : 'N');
											if (path[path[0]+2] & 1) path[path[0]+3] = tree[ path[path[0]+2] ].second | 1;
											else path[path[0]+3] = tree[ path[path[0]+2] ].second  & 0xFFFFFFFE;

											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+2]].first = tree[path[path[0]+3]  ^ 1].first;
											tree[path[path[0]+2]].second = tree[path[path[0]+3]  ^ 1].second;
											if (tree[path[path[0]+2]].second > firstlone) tree[tree[path[path[0]+2]].second].second = path[path[0]+2];

											tree[path[path[0]+3]  ^ 1].first = tree[path[path[0]+1]].first;
											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+1]].first = tree[ path[path[0]+3] ].first;
											tree[path[path[0]+1]].second = tree[ path[path[0]+3] ].second;
											if (tree[path[path[0]+1]].second > firstlone) tree[tree[path[path[0]+1]].second].second = path[path[0]+1];

											//		printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);
											tree[path[path[0]+3] ].first = swapk;

											tree[ path[path[0]+3] ].second = path[path[0]+1] & 0xFFFFFFFE;
											tree[ path[path[0]+3] ^ 1].second = path[path[0]+2] & 0xFFFFFFFE;

											//			printf("%i tmptmp\n",tree[ path[path[0]+3] ].first);




										}


										if (path[0] != 0) { path[0]--;}

										// 		printf("is that fine? %i, %i, %i\n", path[0], path[1], path[2] );show(stdout);




									}else{ // rotation
										//printf("here wego rot !\n");
										path[0]--;
										if (path[0] == 0) {
											swapk = tree[1].first;
											tree[1].first = tree[path[1]].first;
											tree[1].second = path[path[0]+2] &0xFFFFFFFE;
										}else {
											swapk = tree[path[path[0]]].first;
											tree[path[path[0]]].first = tree[path[path[0]+1]].first;
											tree[path[path[0]]].second = path[path[0]+2] | 1;
										}


										tree[path[path[0]+1]].first = tree[path[path[0]+2] ^ 1].first;

										if (tree[path[path[0]+2] ^ 1].second > firstlone){
											tree[path[path[0]+1]].second = tree[path[path[0]+2] ^ 1].second;
											tree[tree[path[path[0]+2] ^ 1].second].second = path[path[0]+1];
										}else tree[path[path[0]+1]].second = tree[path[path[0]+2] ^ 1].second & 0xFFFFFFFE;
										//	if (tree[path[path[0]+2] ^ 1].second > firstlone) tree[tree[path[path[0]+1]].second].second = path[path[0]+1];

										tree[path[path[0]+2]^ 1].first = swapk;
										tree[path[path[0]+2]^ 1].second = path[path[0]+1] & 0xFFFFFFFE;

										//}
										tree[path[path[0]+2]].second &= 0xFFFFFFFE;
										if (path[0] != 0) { path[0]--;}
										// 	show(stdout);
										//
									}}}}}
				}
			}
		}else{
			alloc_mag = 2;
			tree = new pair<Key, unsigned int>[4];
			path = new unsigned int[10];

			tree--;
			tree[1].first = newf;
			tree[1].second = 0;
			size =1;
			firstlone=1 << alloc_mag;
	//		toins_double = 2;
			min_offset = 1;
			max_offset = 1;
		}

    if (min_offset != inorderFirst(path)) {
        printf("min got corrupted in insert!\n");
        printf("real min%i:", inorderFirst(path)); ExOp::show(tree[inorderFirst(path)].first);
        printf("wrong min%i:", min_offset); ExOp::show(getMin());
        this->show();
        LFH_exit(1);
    }
    if (max_offset != inorderLast(path)) {
        printf("max got corrupted in insert!\n");
        printf("real max%i:", inorderLast(path)); ExOp::show(tree[inorderLast(path)].first);
        printf("wrong max%i:", max_offset); ExOp::show(getMax());
        this->show();


        LFH_exit(1);
    }

}


LFHTEMP	void RBTofDoom<Key,void>::remove(const Key &todelete){ // remove 1 element which is equal, if found
    typename RBTofDoom<Key,void>::Iterator ite(*this);
    ite.findEQ(todelete);
    if (ite.isValid()) remove(ite);
}

LFHTEMP	template<class A> void RBTofDoom<Key,void>::remove(const A &todelete){ // remove 1 element which is equal, if found
    typename RBTofDoom<Key,void>::Iterator ite(*this);
    ite.findEQ(todelete);
    if (ite.isValid()) remove(ite);
}

LFHTEMP	void RBTofDoom<Key,void>::remove_static_routine(typename RBTofDoom<Key,void>::Iterator &todel,  pair<Vector<unsigned int>, unsigned int> &list_and_lone){ // delete pointed node and update iterator which will target some element
    if (!todel.isValid()) return;
    //	display(Dbg);
    //	fflush(Dbg);

    LFH_ASSERT(((size - (1 << alloc_mag) + firstlone) & 1) == 1, "Not legal deleting");

  //  printf("single del on "); ExOp::show(*todel);
  //  printf("firstlone is  %i\n", firstlone);

    bool prob;
    unsigned int cur,cur2,cur3;
    unsigned int ipath;
    Key swap;
    // Maintain min-max
    if (size == 1) {
        list_and_lone.first.push_back(0); list_and_lone.second++;
    }else{
        cur = todel.path[0] == 0 ? 1 : todel.path[todel.path[0]];
        // maintain iterator! this is crazy!
        if (!isLeaf(cur)){ // swap node with previous node, preserves structure
            // swap node with a leaf!
            cur2 = cur;
            if (!isLeaf(cur = tree[cur].second)) {
                cur = cur & 0xFFFFFFFE; // left child!
                while(tree[cur].second > 1){
                    todel.path[++(todel.path[0])] = cur;
                    if (tree[cur].second > firstlone) {
                        cur = tree[cur].second;
                        todel.path[++(todel.path[0])] = cur;
                        break;
                    }
                    cur = (tree[cur].second) | 1; // right child!

                }
            }
            todel.path[++(todel.path[0])] = cur;
            ExOp::toMemmove(tree[cur2].first,tree[cur].first);
        }

        if (cur > firstlone){
             list_and_lone.first.push_back(cur); list_and_lone.second++;
            tree[tree[cur].second].second = 0;
        }else{
            cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
            if ((tree[cur].second == 0)&&(tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // annoying red brother... ROTATE!
                cur3 = (cur & 1) ? tree[cur^1].second: tree[cur^1].second & 0xFFFFFFFE;
                swap = tree[cur2].first;
                ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
                ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
                tree[cur3].first = swap;

                tree[cur ^ 1].second = tree[cur3].second;
                if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

                tree[cur2].second = cur3 & 0xFFFFFFFE;
                tree[cur3].second = cur | 1;

                todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
                todel.path[todel.path[0]] = cur3;
                todel.path[0]++;

                cur = todel.path[todel.path[0]];
                cur2 = todel.path[todel.path[0]-1];
            }



            if (tree[cur^1].second > firstlone){ // brother has lone child, use him!
                cur3= tree[cur^1].second;
                if (cur & 1){
                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
                }else{
                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                    ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
                }
                list_and_lone.first.push_back(cur3); list_and_lone.second++;
                tree[cur ^ 1].second =0;

            }else{
                // deleting a inside leaf...
                if (list_and_lone.second != 0) {
                    for(ipath=0; firstlone>=list_and_lone.first[ipath];ipath++);
                    cur3 = list_and_lone.first[ipath];
                    list_and_lone.first.pop_swap(ipath); list_and_lone.second--;
                }else{cur3 = firstlone--;size++; } // printf("WaRnInG! this is bad news (RBTofDOOm) %i is size\n", size); fflush(stdout);

                cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];

                if (tree[cur^1].second > 1){ // brother has siblings, do a rotation
        //			printf("NOT oK!!!! %i \n", tree[cur^1].second);
                    prob = false; // done after this
                    ipath = tree[cur^1].second & 0xFFFFFFFE; // to fill later

                    ExOp::toMemmove(tree[cur].first,tree[cur2].first);

                    if (cur & 1){
                        ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second | 1 ].first);

                        ExOp::toMemmove(tree[cur3].first,tree[cur^1].first);
                        ExOp::toMemmove(tree[cur^1].first,tree[tree[cur^1].second].first);

                    }else{

                        ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second ].first);
                        ExOp::toMemmove(tree[cur3].first,tree[tree[cur^1].second | 1].first);

                    }
                    tree[cur^1].second = cur3;
                    tree[cur3].second = cur^1;

                }else{
                    prob = (tree[cur2].second & 1) == 0; // is black?

                    // send brother to loners


                    tree[cur3].second = cur2;
                    tree[cur2].second = cur3;

                    if (cur & 1){
                        // needs to swap parent
                        ExOp::toMemmove(tree[cur3].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                    }else ExOp::toMemmove(tree[cur3].first,tree[cur ^ 1].first);
                    ipath = cur & 0xFFFFFFFE; // to fill later
                    todel.path[0]--;
                }

                if ((prob)&&(tree[cur].second == 0)) this->resolve_doubleblack(todel); // both child and parent were black, needs rotations! (cur2 is double-black)

                list_and_lone.first.push_back(ipath);
}   }   }	}




LFHTEMP	void RBTofDoom<Key,void>::remove_update(typename RBTofDoom<Key,void>::Iterator &todelete){
	Key tmp = (*todelete);
	this->remove(todelete);
	todelete = this->find_first(tmp);
}


LFHTEMP	void RBTofDoom<Key,void>::remove_mm(typename RBTofDoom<Key,void>::Iterator &todelete){
	Key thatkey = (*todelete);
	this->remove(todelete);
	todelete.findLE(thatkey);
}

LFHTEMP	void RBTofDoom<Key,void>::remove_pp(typename RBTofDoom<Key,void>::Iterator &todelete){
	Key thatkey = (*todelete);
	this->remove(todelete);
	todelete.findGE(thatkey);
}


LFHTEMP	void RBTofDoom<Key,void>::remove(typename RBTofDoom<Key,void>::Iterator &todel){ // delete pointed node and update iterator which will target some element
    if (!todel.isValid()) return;
    //	display(Dbg);
    //	fflush(Dbg);
    bool prob;
    unsigned int cur,cur2,cur3;
    unsigned int ipath;
    Key swap;
    // Maintain min-max
    if (size == 1) {
        delete[](tree+1);
        tree = NULL;
        size=0;
        todel.path[0] = ExCo<unsigned int>::mkMaximum(); // make unvalid
    }else{
        cur = todel.path[0] == 0 ? 1 : todel.path[todel.path[0]];
        // MAINTAIN MIN/MAX
        /*
        if (tree[cur].first == min){
            typename RBTofDoom<Key,void>::Iterator ite(*this); ite.findFirst();
            ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
            if (ipath == cur){
                if (tree[cur].second > firstlone){
                    min = tree[tree[cur].second].first;
                }else min = tree[ todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1] ].first;
            }
        }else if (tree[cur].first == max) {
            typename RBTofDoom<Key,void>::Iterator ite(*this); ite.findLast();
            ipath = ite.path[0] == 0 ? 1 : ite.path[ite.path[0]];
            if (ipath == cur) max = tree[ todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1] ].first;
        }*/

        if (size == (((unsigned int) 1) << (alloc_mag-2)) ){ // Downsize, complete rebalance
            Vector<Key> batch;
            batch.setSize(size-1);
            ipath =0;
            cur2 = inorderFirst(path);
            while(true){
                if (cur == cur2) {ExOp::toMemmove(swap, tree[cur2].first); break;}
                ExOp::toMemmove(batch[ipath++],tree[cur2].first);
                cur2 =inorderNext(path);
            }
            while(ipath< size-1) batch[ipath++] = tree[inorderNext(path)].first;
            delete[](tree+1);
            delete[](path);
            batchInit_routine(batch,false);

            todel.findGE(swap);

        }else{ // RBT updates
            size--;
            // maintain iterator! this is crazy!
            if (!isLeaf(cur)){ // swap node with previous node, preserves structure
                // swap node with a leaf!
                cur2 = cur;
                if (!isLeaf(cur = tree[cur].second)) {
                    cur = cur & 0xFFFFFFFE; // left child!
                    while(tree[cur].second > 1){
                        todel.path[++(todel.path[0])] = cur;
                        if (tree[cur].second > firstlone) {
                            cur = tree[cur].second;
                            todel.path[++(todel.path[0])] = cur;
                            break;
                        }
                        cur = (tree[cur].second) | 1; // right child!

                    }
                }
                todel.path[++(todel.path[0])] = cur;
                ExOp::toMemmove(tree[cur2].first,tree[cur].first);
            }

            if (cur > firstlone){
                if (max_offset == cur) max_offset = tree[cur].second;
                firstlone++;
                tree[tree[cur].second].second = 0;
                if (firstlone != cur) {
                    tree[cur] = tree[firstlone];
                    tree[tree[firstlone].second].second = cur;
                    if (max_offset == firstlone) max_offset = cur;
                }
            }else{

                cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
                if ((tree[cur].second == 0)&&(tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // annoying red brother... ROTATE!


                    cur3 = (cur & 1) ? tree[cur^1].second: tree[cur^1].second & 0xFFFFFFFE;

                    swap = tree[cur2].first;
                    ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
                    ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
                    tree[cur3].first = swap;

                    tree[cur ^ 1].second = tree[cur3].second;
                    if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

                    tree[cur2].second = cur3 & 0xFFFFFFFE;
                    tree[cur3].second = cur | 1;

                    todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
                    todel.path[todel.path[0]] = cur3;
                    todel.path[0]++;

                    cur = todel.path[todel.path[0]];
                    cur2 = todel.path[todel.path[0]-1];
                }



                if (tree[cur^1].second > firstlone){ // brother has lone child, use him!
                    if (tree[cur^1].second == max_offset) max_offset = cur^1;
                    firstlone++;
                    cur3= tree[cur^1].second;

                    if (cur & 1){
                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
                    }else{
                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                        ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
                    }
                    if (cur3 != firstlone){
                    ExOp::toMemmove(tree[cur3].first,tree[firstlone].first);
                    tree[cur3].second = tree[firstlone].second;
                    tree[tree[firstlone].second].second = cur3;
                    if (max_offset == firstlone) max_offset = cur3;
                    }
                    tree[cur ^ 1].second =0;

                }else{
                    // deleting a inside leaf...
                    cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];

                    if (tree[cur^1].second > 1){ // brother has siblings, do a rotation
                        prob = false; // done after this
                        ipath = tree[cur^1].second & 0xFFFFFFFE; // to fill later
                        ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                        if (cur & 1){
                            ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second | 1 ].first);
                            ExOp::toMemmove(tree[firstlone].first,tree[cur^1].first);
                            ExOp::toMemmove(tree[cur^1].first,tree[tree[cur^1].second].first);
                            if (ipath == min_offset) min_offset = cur ^1;
                        }else{
                            ExOp::toMemmove(tree[cur2].first,tree[ tree[cur^1].second ].first);
                            ExOp::toMemmove(tree[firstlone].first,tree[tree[cur^1].second | 1].first);
                            if ((ipath |1) == max_offset) max_offset = firstlone;
                        }
                        tree[cur^1].second = firstlone;
                        tree[firstlone].second = cur^1;
        //					firstlone--;show();	firstlone++;
                    }else{
                        prob = (tree[cur2].second & 1) == 0; // is black?

                        // send brother to loners

                        tree[firstlone].second = cur2;
                        tree[cur2].second = firstlone;

                        if (cur & 1){
                            // needs to swap parent
                            ExOp::toMemmove(tree[firstlone].first,tree[cur2].first);
                            ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                        }else ExOp::toMemmove(tree[firstlone].first,tree[cur ^ 1].first);
                        if ((cur&0xFFFFFFFE) == min_offset) min_offset = cur2;
                        if ((cur|1) == max_offset) max_offset = firstlone;
                        ipath = cur & 0xFFFFFFFE; // to fill later
                        todel.path[0]--;
                    }


                    firstlone--;


                    if ((prob)&&(tree[cur].second == 0)) this->resolve_doubleblack(todel); // both child and parent were black, needs rotations! (cur2 is double-black)

                    cur3 = size - (1 << alloc_mag) + firstlone + 1; // last pair position

                    // move last pair inside
                    if (cur3 != ipath){
                        if (cur3 == min_offset) min_offset = ipath & 0xFFFFFFFE;
                        if ((cur3|1) == max_offset) max_offset = ipath | 1;
                        ExOp::toMemmove(tree[ipath].first, tree[cur3].first );
                        tree[ipath].second = tree[cur3].second;
                        ExOp::toMemmove(tree[ipath | 1].first, tree[cur3 | 1].first );
                        tree[ipath| 1].second = tree[cur3| 1].second;

                        if (tree[ipath].second > firstlone) tree[tree[ipath].second].second = ipath;
                        if (tree[ipath|1].second > firstlone) tree[tree[ipath|1].second].second = ipath | 1;


                        typename RBTofDoom<Key,void>::Iterator ite(*this, tree[cur3].first );
                        for(; (ite.path[0] == 0) || (ite.path[ite.path[0]] != cur3) ; ++ite){

                        }
                        cur3 = (ite.path[0] == 1) ? 1 : ite.path[ite.path[0]-1];
                        tree[cur3].second = (tree[cur3].second & 1) | ipath;
}   }   }   }   }

    if (size == 0) return;
    if (min_offset != inorderFirst(path)) {
        printf("min got corrupted in delete!\n");
        printf("real min%i:", inorderFirst(path)); ExOp::show(tree[inorderFirst(path)].first);
        printf("wrong min%i:", min_offset); ExOp::show(getMin());
        this->show();
        LFH_exit(1);
    }
    if (max_offset != inorderLast(path)) {
        printf("max got corrupted in delete!\n");
        printf("real max%i:", inorderLast(path)); ExOp::show(tree[inorderLast(path)].first);
        printf("wrong max%i:", max_offset); ExOp::show(getMax());
        this->show();
        LFH_exit(1);
    }
}


LFHTEMP	void RBTofDoom<Key,void>::resolve_doubleblack(typename RBTofDoom<Key,void>::Iterator &todel){
Key swap;
unsigned int cur, cur2, cur3;

if (todel.path[0] != 0)
while(true){
    // rotate to guarranty black brother
    if (tree[todel.path[todel.path[0]] ^ 1].second & 1){
        cur = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        swap = tree[cur].first;
        ExOp::toMemmove(tree[cur].first,tree[todel.path[todel.path[0]] ^ 1].first);
        cur2 = (todel.path[todel.path[0]] & 1) ? tree[todel.path[todel.path[0]] ^ 1].second | 1 : tree[todel.path[todel.path[0]] ^ 1].second & 0xFFFFFFFE;
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur2].first);
        tree[cur2].first = swap;

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur2].second;
        tree[cur].second = cur2 &0xFFFFFFFE;
        tree[cur2].second =  todel.path[todel.path[0]] | 1;

        todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
        todel.path[todel.path[0]] = cur2;
        todel.path[0]++;
    }

    cur = (todel.path[todel.path[0]] & 1) ? tree[todel.path[todel.path[0]] ^ 1].second | 1 : tree[todel.path[todel.path[0]] ^ 1].second & 0xFFFFFFFE ;
    if ((tree[cur^1].second & 1)&&(tree[cur^1].second  <=firstlone)) { // brother is red
        // one rotation, then done!

        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        swap = tree[cur2].first;
        ExOp::toMemmove(tree[cur2].first,tree[todel.path[todel.path[0]] ^ 1].first);
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur].first);
        tree[cur].first = swap;

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur].second;
        if (tree[cur].second > firstlone) tree[tree[cur].second].second = todel.path[todel.path[0]] ^ 1;
        tree[cur2].second = (tree[cur2].second & 1) | (cur & 0xFFFFFFFE);
        tree[cur].second =  todel.path[todel.path[0]] & 0xFFFFFFFE;

        tree[cur^1].second &= 0xFFFFFFFE;

        break;
    }else if ((tree[cur].second & 1)&&(tree[cur].second  <=firstlone)){
        // double rotation, then done!
        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        swap = tree[cur2].first;
        ExOp::toMemmove(tree[cur2].first,tree[cur].first);
        cur3 = (todel.path[todel.path[0]] & 1) ? tree[cur].second | 1 : tree[cur].second & 0xFFFFFFFE ;
        ExOp::toMemmove(tree[cur].first, tree[cur3 ^ 1].first );
        ExOp::toMemmove(tree[cur3 ^ 1].first, tree[todel.path[todel.path[0]] ^ 1].first );
        ExOp::toMemmove(tree[todel.path[todel.path[0]] ^ 1].first, tree[cur3].first);
        tree[cur3].first = swap;

        tree[todel.path[todel.path[0]] ^ 1].second = tree[cur3].second;
        if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = todel.path[todel.path[0]] ^ 1;
        tree[cur].second = tree[cur3^1].second;
        if (tree[cur3^1].second > firstlone) tree[tree[cur3^1].second].second = cur;

        tree[cur3].second = todel.path[todel.path[0]] & 0xFFFFFFFE;
        tree[cur3^1].second = cur & 0xFFFFFFFE;
        tree[cur2].second = (tree[cur2].second & 1) | (cur3 & 0xFFFFFFFE);



        break;
        }

    // no rotation... not finished D:
    //	cur3 = tree[path[path[0]]].first;
    //	tree[path[path[0]]].first = -1;
    //	show();
    //	tree[path[path[0]]].first = cur3;

        tree[todel.path[todel.path[0]] ^ 1].second |= 1; // paint brother red
        todel.path[0]--;

        if (todel.path[0] == 0) break;
        if (tree[todel.path[todel.path[0]]].second & 1){
            tree[todel.path[todel.path[0]]].second &= 0xFFFFFFFE;
            break;
            }
    }

}

LFHTEMP	bool RBTofDoom<Key,void>::areLonersHealty()const{
    for(unsigned int i=firstlone; (i++ >> alloc_mag)==0;){
        if (tree[i].second > firstlone) return false;
        if (tree[tree[i].second].second != i) return false;
    }
    return true;
}


LFHTEMP	void RBTofDoom<Key,void>::delete_nodes_routine(pair<Vector<unsigned int>, unsigned int> &list_and_lone){ // listed nodes are assumed to be all the unreachable ones
    unsigned int cur;
    unsigned int cur2,cur3;
    unsigned int i = size - list_and_lone.first.getSize() * 2 + list_and_lone.second;
  //      printf("swapdel start! %c \n", i < (((unsigned int) 1) << (alloc_mag-2)) ? 'Y' : 'N'); fflush(stdout);
    if (i < (((unsigned int) 1) << (alloc_mag-2)) ){ // Downsize, complete rebalance
        if (i == 0) {this->toMemfree(); return;}
        Vector<Key> batch;
        batch.setSize(i);

        cur2 = inorderFirst(path);
        batch[0] = tree[cur2].first;
        for(i=1;i< batch.getSize();i++) { batch[i] = tree[inorderNext(path)].first;}
        delete[](tree+1);
        delete[](path);
        batchInit_routine(batch,false);
       // printf("new size %i, %i\n", size, firstlone);
    }else{
        list_and_lone.first.sort();
       // printf("%i is size (%i >= %i right?)\n", size, firstlone, size - (1 << alloc_mag) + firstlone + 1);
       // printf("got sorted! (first lone was %i)\n", firstlone); fflush(stdout);
       // this->show();
        // swapping loners
        cur = list_and_lone.first.getSize() - list_and_lone.second;
        cur2 = list_and_lone.first.getSize();
        while(cur< cur2){
            firstlone++;
            if (list_and_lone.first[cur] == firstlone) { cur++; continue;} // no swapping needed
            // swap with cur2
            cur2--;
            if (firstlone != tree[tree[firstlone].second].second){
                printf(" %i is nbrem (%i cur cur %i)\n", list_and_lone.second,cur  ,  cur2);
                printf("%i != %i  (p=%i) \n",list_and_lone.first[cur2], tree[tree[list_and_lone.first[cur2]].second].second, tree[list_and_lone.first[cur2]].second);
                ExOp::show(list_and_lone.first);

                this->show();
            }
            LFH_ASSERT(firstlone == tree[tree[firstlone].second].second, "corrupted lone leaf\n");

            ExOp::toMemmove(tree[list_and_lone.first[cur2]].first, tree[firstlone].first);
            tree[list_and_lone.first[cur2]].second = tree[firstlone].second;
            tree[tree[firstlone].second].second = list_and_lone.first[cur2];


            }
        //printf("loners done!\n"); fflush(stdout);
        LFH_ASSERT(this->areLonersHealty(), "Unhealthy\n");


        // swapping inner nodes
        cur = 0;
        cur2 = list_and_lone.first.getSize() - list_and_lone.second -1;
        list_and_lone.second = size - (1 << alloc_mag) + firstlone + 1 - list_and_lone.second;
        size = i;
//        printf("dacur %i %i (firstlo = %i)\n", cur, cur2, firstlone); fflush(stdout);
        while(cur < (cur2+1)){
//             printf("dacur %i %i %i\n", cur, cur2,  list_and_lone.second); fflush(stdout);
             list_and_lone.second -= 2;
             if (((list_and_lone.first[cur2] ^ list_and_lone.second) & 0xFFFFFFFE) == 0 ) { cur2--; continue;} // no swapping needed

//            printf("time to find %i!\n", tree[list_and_lone.second].first  ); fflush(stdout);
            typename RBTofDoom<Key,void>::Iterator ite = this->find_first( tree[list_and_lone.second].first );
            if (!ite.isValid()) {printf("Could not reach %i (at %i)!\n", tree[list_and_lone.second].first, list_and_lone.second ); this->show(); LFH_exit(1);}
            for(; (ite.path[0] == 0) || ((ite.path[ite.path[0]] ^ list_and_lone.second) & 0xFFFFFFFE) != 0 ; ){
                ++ite;
                if (!ite.isValid()) {printf("Could not reach %i (at %i)!\n", tree[list_and_lone.second].first, list_and_lone.second ); this->show(); LFH_exit(1);}
            }

            cur3 = (ite.path[0] == 1) ? 1 : ite.path[ite.path[0]-1];
//            printf("replace %i->%i\n", list_and_lone.second & 0xFFFFFFFE , list_and_lone.first[cur] & 0xFFFFFFFE); fflush(stdout);
            tree[cur3].second = (tree[cur3].second & 1) | (list_and_lone.first[cur] & 0xFFFFFFFE);
            cur3 = list_and_lone.first[cur] & 0xFFFFFFFE;

            ExOp::toMemmove(tree[cur3].first, tree[list_and_lone.second].first );
            ExOp::toMemmove(tree[cur3 | 1].first, tree[list_and_lone.second | 1].first );
            tree[cur3].second = tree[list_and_lone.second].second;
            tree[cur3| 1].second = tree[list_and_lone.second| 1].second;


            if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur3;
            if (tree[cur3|1].second > firstlone) tree[tree[cur3|1].second].second = cur3 | 1;

            cur++;
        }
        //printf("pairs done!\n"); fflush(stdout);
        LFH_ASSERT(this->areLonersHealty(), "Unhealthy\n");

    }
    //printf("swapdel end!\n"); fflush(stdout);
    min_offset = inorderFirst(path);
    max_offset = inorderLast(path);
}

LFHTEMP	void RBTofDoom<Key,void>::remove_subtree_routine( typename RBTofDoom<Key,void>::Iterator &todel, pair<Vector<unsigned int>, unsigned int> &list_and_lone){
    if (todel.path[0] == 0) {todel.path[0] = 0xFFFFFFFF; this->toMemfree(); return;}
    //printf("remove_subtree_routine called on"); ExOp::show(*todel);
    bool prob;
    unsigned int cur = todel.path[todel.path[0]];
    Key swap;
    unsigned int cur2 = (todel.path[0] == 1) ? 1 : todel.path[todel.path[0]-1];
    unsigned int cur3;
    if (cur > firstlone){// printf("remrout: lone\n"); fflush(stdout);
        --todel;
        if (todel.path[0] == 0) {todel.path[0] = 0xFFFFFFFF; this->toMemfree(); return;}
        list_and_lone.first.push_back(cur); list_and_lone.second++;
        cur = todel.path[todel.path[0]];
        cur2 = todel.path[0] == 1 ? 1 : todel.path[todel.path[0]-1];
        tree[cur].second = 0;

        if ((tree[cur^1].second <= firstlone)&&(tree[cur^1].second & 1)){ // Red brother... double rotate!
            //   B
            // B   R
            //    B B
            //   ? ?
            cur3 = (cur & 1) ? tree[cur^1].second : tree[cur^1].second & 0xFFFFFFFE;

            swap = tree[cur2].first;
            ExOp::toMemmove(tree[cur2].first,tree[cur^1].first);
            ExOp::toMemmove(tree[cur^1].first, tree[cur3].first);
            tree[cur3].first = swap;

            tree[cur ^ 1].second = tree[cur3].second;
            if (tree[cur3].second > firstlone) tree[tree[cur3].second].second = cur^1;

            tree[cur2].second = cur3 & 0xFFFFFFFE;
            tree[cur3].second = cur | 1;

            todel.path[todel.path[0]+1] = todel.path[todel.path[0]];
            todel.path[todel.path[0]] = cur3;
            todel.path[0]++;

            cur = todel.path[todel.path[0]];
            cur2 = todel.path[todel.path[0]-1];
        }

        if (tree[cur^1].second > firstlone){ // brother has lone child, use him!
            //   X
            // B   B
            //    R
            cur3= tree[cur^1].second;

            if (cur & 1){
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur3].first);
            }else{
                ExOp::toMemmove(tree[cur].first,tree[cur2].first);
                ExOp::toMemmove(tree[cur2].first,tree[cur ^ 1].first);
                ExOp::toMemmove(tree[cur^1].first,tree[cur3].first);
            }
            list_and_lone.first.push_back(cur3); list_and_lone.second++;
            tree[cur ^ 1].second =0;
        }else{
            cur3 = tree[cur^1].second;
            if ( cur3 == 0) {
                prob = (tree[cur2].second & 1) == 0;
                //   X
                // B   B
           //     printf("in here %i %i (fl%i)\n", list_and_lone.first.last(), cur, firstlone);
           //     printf("in here %i %i\n", tree[list_and_lone.first.last()].first, tree[cur].first);
           //     printf("in here %i %i\n", tree[cur2].first, tree[cur^1].first);


                tree[cur2].second = list_and_lone.first.last();
                tree[list_and_lone.first.last()].second = cur2;

                if (cur & 1){
                //    printf("haha\n");
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur2].first);
                    ExOp::toMemmove(tree[cur2].first, tree[cur^1].first);
                }else{
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur^1].first);
                }


                list_and_lone.first.last() = cur;
                todel.path[0]--;
                if (prob) this->resolve_doubleblack(todel);
           //     this->
;
            } else {
                //   X
                // B   B
                //    R R
                tree[cur^1].second = list_and_lone.first.last();
                tree[list_and_lone.first.last()].second = cur^1;
                ExOp::toMemmove(tree[cur].first, tree[cur2].first);
                if (cur & 1){
                    ExOp::toMemmove(tree[cur2].first, tree[cur3 | 1].first);
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur^1].first);
                    ExOp::toMemmove(tree[cur^1].first, tree[cur3 &0xFFFFFFFE].first);
                }else{
                    ExOp::toMemmove(tree[cur2].first, tree[cur3 &0xFFFFFFFE].first);
                    ExOp::toMemmove(tree[list_and_lone.first.last()].first, tree[cur3 | 1].first);
                }
                list_and_lone.first.last() = cur3;
            }
            list_and_lone.second--;

        //    printf("dblack is at %i\n",  todel.path[0] == 0 ? 1 :  todel.path[todel.path[0]]);

        }
    }else{
        cur2 = (todel.path[0] == 1) ? 1 : todel.path[todel.path[0]-1];
        //if (!delete_parent_as_well) ExOp::toMemmove(swap: tree[cur2]);
        prob = ((tree[cur2].second & 1) == 0)&&((tree[cur^1].second > firstlone)||((tree[cur^1].second & 1) == 0));  // are parent and brother both black?

        ExOp::toMemmove(tree[cur2].first, tree[cur^1].first);
        if (tree[cur^1].second > firstlone){
            tree[cur2].second = tree[cur^1].second;
            tree[tree[cur^1].second].second = cur2;
        }else tree[cur2].second = tree[cur^1].second & 0xFFFFFFFE;
        unsigned int i;
        list_and_lone.first.push_back(cur);
        i=list_and_lone.first.getSize();
        if (tree[cur].second > 1) {
            list_and_lone.first.push_back(tree[cur].second);
            for(;i< list_and_lone.first.getSize();i++){
                if (list_and_lone.first[i] > firstlone) {list_and_lone.second++; continue;}
                if (tree[list_and_lone.first[i]].second > 1) list_and_lone.first.push_back(tree[list_and_lone.first[i]].second);
                if (tree[list_and_lone.first[i]^1].second > 1) list_and_lone.first.push_back(tree[list_and_lone.first[i]^1].second);
            }
        }
        todel.path[0]--;
        //printf("does res %c\n", prob ? 'Y' : 'N');
        if (prob) {
            this->resolve_doubleblack(todel);

        }
    }
    min_offset = inorderFirst(todel.path);
    max_offset = inorderLast(todel.path);
}



    LFHTEMP	template<class SORTER> void RBTofDoom<Key,void>::intersection(vector<Key> &f_out, const SORTER &Query) const{
		if (size ==0) return;

		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = this->getMin();
		Key tmax = this->getMax();

		while( true ){


			while(!(SETCMP_DISJOINT & Query(tmin,tmax))){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if (!(SETCMP_DISJOINT & Query(tree[tree[path[depth]].second].first,tree[tree[path[depth]].second].first))) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
			if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = this->getMax();
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;
		}
	}
	LFHTEMP	void RBTofDoom<Key,void>::intersectionInterval(vector<Key> &f_out, const Key &q_min, const Key &q_max) const{
		if (size ==0) return;
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = this->getMin();
		Key tmax = this->getMax();

		while( true ){


			while((tmin <= q_max)&&(tmax >= q_min)){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if ((q_min <= tree[tree[path[depth]].second].first)&&(q_max >=tree[tree[path[depth]].second].first)) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;

			if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = this->getMax();
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;
		}


		}
	LFHTEMP	template<class SORTER> void RBTofDoom<Key,void>::intersection(Vector<Key> &f_out, const SORTER& Query) const{
		if (size == 0) return;
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = this->getMin();
		Key tmax = this->getMax();

		while( true ){
			while(!(SETCMP_DISJOINT & Query(tmin,tmax))){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if (!(SETCMP_DISJOINT & Query(tree[tree[path[depth]].second].first))) f_out.push_back(tree[tree[path[depth]].second].first);
					}
					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}
			depth--;
			if (depth == 0xFFFFFFFF) break;
			if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = this->getMax();
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;
		}
	}
	LFHTEMP	void RBTofDoom<Key,void>::intersectionInterval(Vector<Key> &f_out, const Key &q_min, const Key &q_max) const{
		if (size == 0) return;
		unsigned int path[65536];
		unsigned int depth =0;
		path[0] = 1;
		//	unsigned int depth_min=0;
		Key tmin = this->getMin();
		Key tmax = this->getMax();


		while( true ){


			while((tmin <= q_max)&&(tmax >= q_min)){
				if (tree[path[depth]].second - 2 > firstlone-2) {
					if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);
					if (tree[path[depth]].second > firstlone){
						if ((q_min <= tree[tree[path[depth]].second].first)&&(q_max >=tree[tree[path[depth]].second].first)) f_out.push_back(tree[tree[path[depth]].second].first);
					}

					break;
				}
				tmax = tree[path[depth]].first;
				path[depth+1] = tree[path[depth]].second & 0xFFFFFFFE;
				depth++;
			}



			depth--;
			if (depth == 0xFFFFFFFF) break;
			if ((q_min <= tree[path[depth]].first)&&(q_max >=tree[path[depth]].first)) f_out.push_back(tree[path[depth]].first);

			tmin = tree[path[depth]].first;
			if (depth == 0) tmax = this->getMax();
			else tmax = tree[path[depth-1]].first;

			// erase current parent in history, no longer needed
			path[depth] = tree[path[depth]].second | 1;

		}
    }


LFHTEMP	void RBTofDoom<Key,void>::removeRange(const Key &q_min, const Key &q_max){
//	vector<Key> fin;
//	intersection(fin, q_min,q_max);
//	for(int i=0 ;i< fin.size();i++) remove(fin[i]);
    while (true){
    typename RBTofDoom<Key,void>::Iterator ite = this->find_first(q_min);
        if (!ite.isValid()) break;
        if ((*ite) > q_max) break;
        this->remove(ite);
    }
}

LFHTEMP	void RBTofDoom<Key,void>::removeRange(const Key &q_min, typename RBTofDoom<Key,void>::Iterator &i_max){
	Key maxkey = (*i_max);
    while((*i_max) >= q_min){
        this->remove(i_max);
		i_max.findLE(maxkey);
		if (!i_max.isValid()) break;
    }
}

LFHTEMP	void RBTofDoom<Key,void>::removeRange(typename RBTofDoom<Key,void>::Iterator &i_min, const Key &q_max){
	Key minkey = (*i_min);
    while ((*i_min) <= q_max){
        this->remove(i_min);
        i_min.finGE(minkey);
        if (!i_min.isValid()) break;
    }
}


LFHTEMP template<class OKEY> void RBTofDoom<Key,void>::removeRange(const OKEY &q_min, typename RBTofDoom<Key,void>::Iterator &i_max){
	Key maxkey = (*i_max);
    while((*i_max) >= q_min){
        this->remove(i_max);
		i_max.findLE(maxkey);
		if (!i_max.isValid()) break;
    }
}

LFHTEMP	template<class OKEY> void RBTofDoom<Key,void>::removeRange(typename RBTofDoom<Key,void>::Iterator &i_min, const OKEY &q_max){
	Key minkey = (*i_min);
    while ((*i_min) <= q_max){
        this->remove(i_min);
        i_min.findGE(minkey);
        if (!i_min.isValid()) break;
    }
}

LFHTEMP	template<class OKEY> void RBTofDoom<Key,void>::range_move_routine(const OKEY& q_min,const OKEY& q_max, Vector<Key> &moveout){

    typename RBTofDoom<Key,void>::Iterator ite(*this, q_min);
    while (true){
        if (!ite.isValid()) break;
        if ((*ite) > q_max) break;
       moveout.mempush_back(*ite);
    }

}

// can be improved...
LFHTEMP	void RBTofDoom<Key,void>::overwriteRange(const Key &q_min, const Key &q_max, const RBTofDoom<Key,void>& replacement){
    if (replacement.getSize() == 0) return removeRange(q_min,q_max);
    if ((replacement.getMin() < q_min)||(replacement.getMax() > q_max)){
        printf("illegal overwriteRange, inserted elements must be within range\n");
        ExOp::show(replacement.min);
        ExOp::show(replacement.max);
        ExOp::show(q_min);
        ExOp::show(q_max);
        LFH_exit(1);
    }

    typename RBTofDoom<Key,void>::Iterator ite = this->find_first(q_min);
    typename RBTofDoom<Key,void>::Iterator ite_in = replacement.first();
    unsigned int flag =0;
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) flag |= 1;
        if (!ite_in.isValid()) flag |= 2;
        if (flag != 0) break;
        ++ite;    ++ite_in;
    }
    Key  cur;
    switch(flag){
        case 1: // has more to insert!
                ite = this->find_first(q_min);
                ite_in = replacement.first();
                while (true){
                    if ((!ite.isValid())||((*ite) > q_max)) break;
                    if (!ite_in.isValid()) break;
                    ite.ordering_preserving_reference() = (*ite_in);
                    ++ite;    ++ite_in;
                }

                do{
                this->insert(*ite_in);
                ++ite_in;
                } while(ite_in.isValid());
            return;
        case 2: // has more to delete!
            cur = (*ite);
            while (true){
                this->remove(ite);
                ite = this->find_first(cur);
                if (!ite.isValid()) break;
                if ((*ite) > q_max) break;
            }
            break;
        default:// exact replacement
            break;
    }
    ite = this->find_first(q_min);
    ite_in = replacement.first();
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) break;
        if (!ite_in.isValid()) break;
        ite.ordering_preserving_reference() = (*ite_in);
        ++ite;    ++ite_in;
    }


}


// can be improved...
LFHTEMP	template<class OKEY> void RBTofDoom<Key,void>::overwriteRange(const OKEY &q_min, const OKEY &q_max, const RBTofDoom<Key,void>& replacement){
    if (replacement.getSize() == 0) return removeRange(q_min,q_max);
    if ((replacement.getMin() < q_min)||(replacement.getMax() > q_max)){
        printf("illegal overwriteRange, inserted elements must be within range\n");
        ExOp::show(replacement.min);
        ExOp::show(replacement.max);
        ExOp::show(q_min);
        ExOp::show(q_max);
        LFH_exit(1);
    }

    typename RBTofDoom<Key,void>::Iterator ite = this->find_first(q_min);
    typename RBTofDoom<Key,void>::Iterator ite_in = replacement.first();
    unsigned int flag =0;
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) flag |= 1;
        if (!ite_in.isValid()) flag |= 2;
        if (flag != 0) break;
        ++ite;    ++ite_in;
    }
    Key  cur;
    switch(flag){
        case 1: // has more to insert!
                ite = this->find_first(q_min);
                ite_in = replacement.first();
                while (true){
                    if ((!ite.isValid())||((*ite) > q_max)) break;
                    if (!ite_in.isValid()) break;
                    ite.ordering_preserving_reference() = (*ite_in);
                    ++ite;    ++ite_in;
                }

                do{
                this->insert(*ite_in);
                ++ite_in;
                } while(ite_in.isValid());
            return;
        case 2: // has more to delete!
            cur = (*ite);
            while (true){
                this->remove(ite);
                ite = this->find_first(cur);
                if (!ite.isValid()) break;
                if ((*ite) > q_max) break;
            }
            break;
        default:// exact replacement
            break;
    }
    ite = this->find_first(q_min);
    ite_in = replacement.first();
    while (true){
        if ((!ite.isValid())||((*ite) > q_max)) break;
        if (!ite_in.isValid()) break;
        ite.ordering_preserving_reference() = (*ite_in);
        ++ite;    ++ite_in;
    }


}


LFHTEMP	template<class OKEY> void RBTofDoom<Key,void>::removeRange_faster(const OKEY &q_min, const OKEY &q_max, bool doit){
    // remove red nodes in range (should not rotate at all)
    pair<Vector<unsigned int>, unsigned int> list_and_lone; list_and_lone.second =0;

    typename RBTofDoom<Key,void>::Iterator ite(this);
    unsigned int cur;
    int flag =0;

    if ((this->getMax() < q_min)||(this->getMin() > q_max)) return;

    if (max <= q_max) flag |= 1;
    if (min >= q_min) flag |= 2;

    if (flag == 3) {this->toMemfree(); return;}
    unsigned int cflag;
    while(true){
        if (doit) this->show();
        ite.path[0] = 0;
        switch(flag){
            case 0:
            for(cflag=0; cflag == 0;){
                while(ExOp::isLT(tree[ite.get_index()].first,q_min)){
                    if (!ite.hasRightChild()) {delete_nodes_routine(list_and_lone);return;}
                    if ((cur =tree[ite.get_index()].second) > firstlone){
                        if ((ExOp::isGE(tree[cur].first, q_min))&&(ExOp::isLE(tree[cur].first, q_max))) {
                            ite.path[++ite.path[0]] = cur;
                            remove_static_routine(ite,list_and_lone);
                        }
                        delete_nodes_routine(list_and_lone);return;
                    }
                    ite.path[++ite.path[0]] = cur | 1;
                }
                if (ExOp::isLE(tree[ite.get_index()].first,q_max)) {cflag = ite.path[0]; break;}
                while(ExOp::isGT(tree[ite.get_index()].first,q_max)){
                    if (!ite.hasLeftChild()) {delete_nodes_routine(list_and_lone);return;}
                    cur = (tree[ite.get_index()].second & 0xFFFFFFFE); // should exist, since q_min < q_max < max
                    ite.path[++ite.path[0]] = cur;
                }
                if (ExOp::isGE(tree[ite.get_index()].first,q_min)) {cflag = ite.path[0]; break;}
            }
            break;
            case 1:
           //     printf("hahaha\n");
                if (ExOp::isLT(tree[1].first, q_min)) {
                    cur =1;
                    do{
                        if ((cur = tree[cur].second) > firstlone){ // deleting that last lone
                            ite.path[++ite.path[0]] = cur; remove_static_routine(ite,list_and_lone);  delete_nodes_routine(list_and_lone);
                            ite = this->last();
                            max = *ite;
                            return;
                        }
                        cur |=1;
                        if (cur == 1){
                            delete_nodes_routine(list_and_lone);
                            ite = this->last();
                            max = *ite;
                            return;
                        }
                        ite.path[++ite.path[0]] = cur;
                    }while(ExOp::isLT(tree[cur].first, q_min));
                }
                cflag = ite.path[0];
            break;
            case 2:
            //   printf("hehehe\n");
                if (ExOp::isGT(tree[1].first, q_max)) {
                    cur =1;
                    do{
                        cur = tree[cur].second;
                        if ((cur > firstlone)||((cur & 0xFFFFFFFE) == 0)) {
                             delete_nodes_routine(list_and_lone);
                             ite = this->first();
                             min = *ite;
                             return;
                        }
                        cur &= 0xFFFFFFFE;
                      //  printf("%i ttt\n", cur);
                        ite.path[++ite.path[0]] = cur;
                    }while(ExOp::isGT(tree[cur].first, q_max));
                }
                cflag = ite.path[0];
            break;
        }
      //  printf("%i\n", cflag);
      //  printf("highest in is"); ExOp::show(*ite);
        // trying Right
        if (ite.hasRightChild()){
            if ((cur =tree[ite.get_index()].second) > firstlone){
                if (ExOp::isLE(tree[cur].first, q_max)) {ite.path[++ite.path[0]] = cur; remove_subtree_routine(ite, list_and_lone);}
                else remove_static_routine(ite,list_and_lone);
                continue;
            }else{
                ite.path[++ite.path[0]] = cur | 1;
                while(ExOp::isGT(tree[ite.get_index()].first,q_max)){
                    if (!ite.hasLeftChild()) {ite.path[0] = cflag;break;}
                    cur = (tree[ite.get_index()].second & 0xFFFFFFFE); // should exist, since q_min < q_max < max
                    ite.path[++ite.path[0]] = cur;
                }
                if (ite.path[0] != cflag){
                    if (!ite.hasLeftChild()) remove_static_routine(ite,list_and_lone);
                    else {
                        cur = tree[ite.get_index()].second & 0xFFFFFFFE;
                        ite.path[++ite.path[0]] = cur;
                        remove_subtree_routine(ite, list_and_lone);
                    }
                    continue;
                }
               // printf("continue! on "); ExOp::show(*ite); fflush(stdout);
            }
        }

        // trying left
        if (ite.hasLeftChild()){
            cur = (tree[ite.get_index()].second & 0xFFFFFFFE);
            ite.path[++ite.path[0]] = cur;
            while(ExOp::isLT(tree[ite.get_index()].first,q_min)){
                if (!ite.hasRightChild()) {ite.path[0] = cflag;break;}
                if ((cur =tree[ite.get_index()].second) > firstlone){
                    if (ExOp::isGE(tree[cur].first,q_min)) {
                        ite.path[++ite.path[0]] = cur;
                        remove_static_routine(ite,list_and_lone);
                    }
                    ite.path[0] = cflag;break;
                }
                cur = (tree[ite.get_index()].second | 1); // should exist, since q_min < q_max < max
                ite.path[++ite.path[0]] = cur;
            }
            if (ite.path[0] != cflag){
                if (!ite.hasRightChild()) remove_static_routine(ite,list_and_lone);
                else{
                    cur = tree[ite.get_index()].second;
                    ite.path[++ite.path[0]] = cur | (cur > firstlone ? 0 : 1);
                    remove_subtree_routine(ite, list_and_lone);
                }
                continue;
            }
          //   printf("continue!\n"); fflush(stdout);
        }

        //printf("final!\n");    printf("firstlone is  %i\n", firstlone);fflush(stdout);
        remove_static_routine(ite,list_and_lone);
        break;
    }
   // this->show();

    delete_nodes_routine(list_and_lone);

    //printf("recoverytime\n"); fflush(stdout);

    switch(flag){
        case 1:
            ite = this->last();
            max = *ite;
        break;
        case 2:
            ite = this->first();
            min = *ite;
        break;

    }
    return;
}



LFHTEMP	template<class OKEY> void RBTofDoom<Key,void>::removeRange(const OKEY &q_min, const OKEY &q_max){
    while (true){
    typename RBTofDoom<Key,void>::Iterator ite(*this, q_min);
        if (!ite.isValid()) break;
        if (ExOp::isGT((*ite),q_max)) break;

        this->remove(ite);
    }

}

LFHTEMP	template<class OKEY, class FF> void RBTofDoom<Key,void>::removeRange(const OKEY &q_min, const OKEY &q_max, const FF& filter_func){
    OKEY tmpmin = q_min;
    typename RBTofDoom<Key,void>::Iterator ite(*this,q_min);
    while (true){
        if (!ite.isValid()) break;
        if (ExOp::isGT((*ite),q_max)) break;
        if (FF(*ite)) {this->remove(ite); ite = this->find_first(tmpmin);}
        else {tmpmin = q_min; ++ite;}
    }
}

LFHTEMP	unsigned int RBTofDoom<Key,void>::inorderFirst(unsigned int *path) const{
    unsigned int cur = 1;
    path[0] =0;
    while(tree[cur].second > 1){
        if (tree[cur].second > firstlone) break;
        else cur = tree[cur].second & 0xFFFFFFFE;
        path[++(path[0])] =cur;
    }
    return(cur);
}

LFHTEMP	unsigned int RBTofDoom<Key,void>::inorderLast(unsigned int *path) const{
    unsigned int cur = 1;
    path[0] =0;
    while(tree[cur].second > 1){
        if (tree[cur].second > firstlone) {cur = tree[cur].second; path[++(path[0])] =cur; return cur;}
        else cur = tree[cur].second | 1;
        path[++(path[0])] =cur;
    }
    return(cur);
}

	LFHTEMP	unsigned int RBTofDoom<Key,void>::inorderNext(unsigned int *path) const{

		unsigned int cur;
		bool branch;
		if (path[0] == 0) {
			if (tree[1].second == 0) {return(0xFFFFFFFF);}
			branch = false;
		} else{
			branch = ((tree[path[path[0]]].second < 2)||(path[path[0]] > firstlone))  ;
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no right child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			if ((path[0]>0)&&((tree[path[path[0]]].second > firstlone)||((tree[path[path[0]]].second | 1) == path[path[0]+1]))) {path[0]--;
				while((path[0]>0)&&((tree[path[path[0]]].second | 1) == path[path[0]+1])) path[0]--;
			}

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
				return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}else return(path[path[0]]);

		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (tree[cur].second > firstlone) { cur = tree[cur].second; path[++(path[0])] =cur;}
			else {cur = tree[cur].second | 1;
				path[0]++;
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);
				path[path[0]] =cur;
				while(tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (tree[cur].second > firstlone) break;
					else cur = tree[cur].second & 0xFFFFFFFE;
					path[++(path[0])] =cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
			}
		}

		return(cur);
	}
	LFHTEMP	bool RBTofDoom<Key,void>::isRed(unsigned int w) const{

		return((w > firstlone)||((tree[w].second & 1)&&(tree[w].second <= firstlone)));

	}
    LFHTEMP	bool RBTofDoom<Key,void>::isLeaf(unsigned int w) const{

			return((w > firstlone)||(tree[w].second< 2));
		}

LFHTEMP	uint32_t RBTofDoom<Key,void>::findGE(const Key& what) const{
    Iterator ite(this);
    if (ite.findGE(what)) return ite.get_index();
	return 0;
}
LFHTEMP	template<class C> uint32_t RBTofDoom<Key,void>::findGE(const C& what) const{
    Iterator ite(this);
    if (ite.findGE(what)) return ite.get_index();
	return 0;
}
LFHTEMP	uint32_t RBTofDoom<Key,void>::findEQ(const Key& what) const{
    Iterator ite(*this);
    if (ite.findEQ(what)) return ite.get_index();
	return 0;
}
LFHTEMP	template<class C> uint32_t RBTofDoom<Key,void>::findEQ(const C& what) const{
    Iterator ite(*this);
    if (ite.findEQ(what)) return ite.get_index();
	return 0;
}
LFHTEMP	uint32_t RBTofDoom<Key,void>::findGT(const Key& what) const{
    Iterator ite(*this);
    if (ite.findGT(what)) return ite.get_index();
	return 0;
}
LFHTEMP	template<class C> uint32_t RBTofDoom<Key,void>::findGT(const C& what) const{
    Iterator ite(*this);
    if (ite.findGT(what)) return ite.get_index();
	return 0;
}
LFHTEMP	uint32_t RBTofDoom<Key,void>::findLE(const Key& what) const{
    Iterator ite(*this);
    if (ite.findLE(what)) return ite.get_index();
	return 0;
}
LFHTEMP	template<class C> uint32_t RBTofDoom<Key,void>::findLE(const C& what) const{
    Iterator ite(*this);
    if (ite.findLE(what)) return ite.get_index();
	return 0;
}
LFHTEMP	uint32_t RBTofDoom<Key,void>::findLT(const Key& what) const{
    Iterator ite(*this);
    if (ite.findLT(what)) return ite.get_index();
	return 0;
}
LFHTEMP	template<class C> uint32_t RBTofDoom<Key,void>::findLT(const C& what) const{
    Iterator ite(*this);
    if (ite.findLT(what)) return ite.get_index();
	return 0;
}


LFHTEMP template<class Q, class R> void RBTofDoom<Key,void>::queryParallel(ThreadBase &tb, Q& query, R& receiver,
FUNCTOREQUIRES_DEF2(Q, bool, const Key& , const Key&), ACCEPTOR_DEF(R, Key)) const{ //
class Task{
public:
    const RBTofDoom<Key,void>& _this;
    ThreadBase &tb;
    Q& query;
    R& receiver;
    Task(const RBTofDoom<Key,void>& __this, ThreadBase &_tb,Q& _query,R& _receiver): _this(__this), tb(_tb),query(_query), receiver(_receiver){}
    void operator()(uint32_t node, uint32_t minoff, uint32_t maxoff){
       // printf("sta: %i,%i,%i\n", _this.tree[node].first, _this.tree[minoff].first, _this.tree[maxoff].first);
        while(true){
            receiver(_this.tree[node].first);
            if (_this.tree[node].second == 0) break;
            else if (_this.tree[node].second > _this.firstlone) {
                receiver(_this.tree[_this.tree[node].second].first);
                break;
            }else if (query(_this.tree[minoff].first, _this.tree[node].first)){
                if (query(_this.tree[node].first, _this.tree[maxoff].first)){
                    //printf("sub: %i,%i,%i\n", _this.tree[_this.tree[node].second | 1].first, _this.tree[node].first, _this.tree[maxoff].first);
                    tb.submit(*this, _this.tree[node].second | 1, node, maxoff);

                }
                maxoff = node; node = _this.tree[node].second & 0xFFFFFFFE;
            }else if (query(_this.tree[node].first, _this.tree[maxoff].first)){
                minoff = node; node = _this.tree[node].second | 1;
            }else break;
        }
    }
};
    Task locask(*this, tb,query,receiver);
    if (this->getSize() <2){
        if (this->getSize() ==0) return;
        receiver(tree[1].first); // so root is not a leaf
    }else{
        locask(1,min_offset, max_offset);
        tb.joinAll();
    }
}

	LFHTEMP	void RBTofDoom<Key,void>::show(FILE* f, int level) const{
		unsigned int buffer[1024];
		char funbuf[256];
		memset(funbuf,' ', sizeof(char) * 256);
		funbuf[255] = '\0';
		int loop =0;
		if (size == 0) {fprintf(f,"Empty RBTree!\n"); return;}
		fprintf(f,"RBTree with %i nodes on interval [", size);
		ExOp::show(this->getMin(),f,1);
		fprintf(f," ; ");
		ExOp::show(this->getMax(),f,1);
		fprintf(f,"] allocsize = %i\n", 1 << alloc_mag);
		int fun;

		for(unsigned int cur = inorderFirst(buffer); cur != 0xFFFFFFFF;cur = inorderNext(buffer)){
			if (cur > firstlone) {
				fun = 0;
				for (loop = buffer[0]-1;loop > 0;loop--) if (!isRed(buffer[loop])) fun++;
				fprintf(f,"%i%sR{%i}%c",fun , funbuf + 255 - buffer[0]+1,cur, (tree[tree[cur].second].second != cur) ? 'E' : ' ');

			}else fprintf(f,"%s%c{%i}  ", funbuf + 255 - buffer[0], isRed(cur) ? 'R' : 'B',cur);
			ExOp::show(tree[cur].first,f,1);
			fprintf(f,"\n");
		}
	}
LFHTEMP	ERRCODE RBTofDoom<Key,void>::save(FILE* f) const{
    //saves an ordered array!
    unsigned int buffer[1024];
    ERRCODE fout = (fwrite(&size ,sizeof(unsigned int),1 ,f) == 1) ? 0 : 1;
    if (size == 0) return fout;
    for(unsigned int cur = inorderFirst(buffer); cur != 0xFFFFFFFF;cur = inorderNext(buffer)){
        fout |= ExOp::save(tree[cur].first, f);
    }
    return fout;
}
LFHTEMP	ERRCODE RBTofDoom<Key,void>::load(FILE* f){
	unsigned int c,i_size;
	this->toMemfree();
	if (1 != fread(&i_size ,sizeof(unsigned int),1 ,f)) return 1;

	if (i_size == 0) return 0;
	this->makeEmptyTree(i_size);

	// fixed known size!
	unsigned int buffer[1024];
	unsigned int cur = inorderFirst(buffer);
	min_offset = cur;
	ExOp::load(tree[cur].first, f);
	c = i_size -1;
	if (i_size > 1) while((c--) != 0){
		cur = inorderNext(buffer);
		ExOp::load(tree[cur].first, f);
	}
	max_offset = cur;
	return 0;
}
LFHTEMP	typename RBTofDoom<Key,void>::Iterator RBTofDoom<Key,void>::find(const Key& tval) const{ // find the smallest which is >= key
    typename RBTofDoom<Key,void>::Iterator fout = find_first(tval);
    if ((fout.isValid())&&(*fout != tval)) fout.path[0] = 0xFFFFFFFF;
    return(fout);
    }

LFHTEMP void RBTofDoom<Key,void>::test() const{
    printf("testing... ");
    typename RBTofDoom<Key,void>::Iterator itemin = this->first();
    unsigned int nbitem;
    for(nbitem=0;itemin.isValid();++itemin){
        typename RBTofDoom<Key,void>::Iterator itemax = this->find(*itemin);
        nbitem++;

   //     for(++itemax;itemax.isValid();++itemax){
            if ((*itemin) > (*itemax)) {printf("failed!\n"); return;}
   //     }

    }
    if (nbitem != size) printf("failed!\n");
    else printf("passed!\n");
}



LFHTEMP template<class F> bool RBTofDoom<Key,void>::first(RBTofDoom<Key,void>::QueryIterator& ite, const F &f)const{
    SETCMP_enum sc;
    if (size < 2){
        if (size == 0) return false;
        ite.path[0] =0;
        return f(tree[1].first);
    }else if (f(min,max) & SETCMP_DISJOINT) return false;
    unsigned int cur=1;
    ite.path[0] =0;
    ite.iaindex[0] =0; ite.iaindex[1] =32;
    ite.iaval[0] = this->getMin(); ite.iaval[32] = this->getMax();
    while(tree[cur].second > 1){
        if (tree[cur].second > firstlone) {
            // lone children
            if ((f(tree[cur].first)& SETCMP_DISJOINT) == 0) return true;
            if ((f(tree[tree[cur].second].first)& SETCMP_DISJOINT) == 0){
                path[++(path[0])] =tree[cur].second;
                return true;
            }
            break;
            // needs to back up!

        } else if (((f(ite.iaval[ite.iaindex[0]] ,tree[cur].first)) & SETCMP_DISJOINT) == 0){
            ite.iaindex[1]++; ite.iaval[ite.iaindex[1]] = tree[cur].first;
            cur = tree[cur].second & 0xFFFFFFFE;
            path[++(path[0])] =cur;
        }else if (((f(tree[cur].first, ite.iaval[ite.iaindex[1]] )) & SETCMP_DISJOINT) == 0){
            ite.iaindex[0]++; ite.iaval[ite.iaindex[0]] = tree[cur].first;
            cur = tree[cur].second | 1;
            path[++(path[0])] =cur;
        }else break;// else needs to back up!, but *should* not happen
    }
    return false;
}

LFHTEMP template<class F> bool RBTofDoom<Key,void>::next(RBTofDoom<Key,void>::QueryIterator&, const F&)const{
    return true;
}


LFHTEMP	Key& RBTofDoom<Key,void>::orderpreserve_deref(const typename RBTofDoom<Key,void>::Iterator& ite){return tree[ite.get_index()].first;}
LFHTEMP	Key& RBTofDoom<Key,void>::orderpreserve_deref_index(uint32_t index){return tree[index].first;}


LFHTEMP	RBTofDoom<Key,void>::Iterator::Iterator(const RBTofDoom<Key,void>& i_target) : target(i_target){}
LFHTEMP	RBTofDoom<Key,void>::Iterator::Iterator(const RBTofDoom<Key,void>& _target, const Key& key): target(_target){this->findGE(key);}
LFHTEMP	template <class O> RBTofDoom<Key,void>::Iterator::Iterator(const RBTofDoom<Key,void>& _target, const O& key): target(_target){this->findGE(key);}
LFHTEMP	RBTofDoom<Key,void>::Iterator::Iterator(const RBTofDoom<Key,void>& _target, const Key& key, SETCMP_enum q): target(_target){
	switch(q){
		case SETCMP_GE: this->findGE(key); break;
		case SETCMP_GT: this->findGT(key); break;
		case SETCMP_LE: this->findLE(key); break;
		case SETCMP_LT: this->findLT(key); break;
		default: this->find(key); break;
	}
}

LFHTEMP	template <class O> RBTofDoom<Key,void>::Iterator::Iterator(const RBTofDoom<Key,void>& _target, const O& key, SETCMP_enum q): target(_target){
	switch(q){
		case SETCMP_GE: this->findGE(key); break;
		case SETCMP_GT: this->findGT(key); break;
		case SETCMP_LE: this->findLE(key); break;
		case SETCMP_LT: this->findLT(key); break;
		default: this->find(key); break;
	}
}

LFHTEMP	const typename RBTofDoom<Key,void>::Iterator& RBTofDoom<Key,void>::Iterator::operator++(){
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return(*this);}
			branch = false;
		} else{
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no right child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			if ((path[0]>0)&&((target.tree[path[path[0]]].second > target.firstlone)||((target.tree[path[path[0]]].second | 1) == path[path[0]+1]))) {path[0]--;
				while((path[0]>0)&&((target.tree[path[path[0]]].second | 1) == path[path[0]+1])) path[0]--;
			}

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
            if (((target.tree[1].second & 0xFFFFFFFE)!= path[1])||(target.tree[1].second > target.firstlone)) path[0] = ExCo<unsigned int>::mkMaximum();
		//	return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}

		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (target.tree[cur].second > target.firstlone) { cur = target.tree[cur].second; path[++(path[0])] =cur;}
			else {cur = target.tree[cur].second | 1;
				path[0]++;

				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);

				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);
				path[path[0]] =cur;




				while(target.tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (target.tree[cur].second > target.firstlone) break;
					else cur = target.tree[cur].second & 0xFFFFFFFE;
					path[++(path[0])] =cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
			}
		}
		return(*this);
}
LFHTEMP	const typename RBTofDoom<Key,void>::Iterator& RBTofDoom<Key,void>::Iterator::operator--(){
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) {path[0] = ExCo<unsigned int>::mkMaximum(); return(*this);}
			branch = false;
		} else if (path[path[0]] > target.firstlone) {
		    path[0]--; return(*this);
		}else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
			//	 printf("%i.(l) = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
		}
		//	printf("has %c\n", branch? 'N': 'Y'); fflush(stdout);
		if (branch) {
			// no left child!

			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);
			path[0]--;
			while((path[0]>0)&&((target.tree[path[path[0]]].second & 0xFFFFFFFE) == path[path[0]+1])) path[0]--;


			//	printf("%i.c = %i\n", path[path[0]], tree[path[path[0]]].second); fflush(stdout);

			if (path[0]==0) {
            if ((target.tree[1].second | 1) != path[1]) path[0] = ExCo<unsigned int>::mkMaximum();
		//	return( (tree[1].second == path[1])&&(tree[1].second <= firstlone) ?  1 : 0xFFFFFFFF);
			}

		}else{
			// has left
			cur = (path[0] ==0) ? 1 : path[path[0]];
			cur = target.tree[cur].second & 0xFFFFFFFE; path[++(path[0])] =cur;


				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
				//	printf("bed %i = %i; %i= %i\n",path[0],path[1],path[2],path[3]);





				while(target.tree[cur].second > 1){
					//		printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);
					if (target.tree[cur].second > target.firstlone) {
					    cur = target.tree[cur].second; path[++(path[0])] =cur;
					    break;
					    }
					cur = target.tree[cur].second | 1;
					path[++(path[0])] = cur;

				}
				//	printf("%i.c = %i\n", cur, tree[cur].second); fflush(stdout);

		}
return(*this);
}


LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findFirst(){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur = 1;
	path[0] =0;
	//  printf("%i.l = %i\n",cur, tree[cur].second); fflush(stdout);
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) break;
		else cur = target.tree[cur].second & 0xFFFFFFFE;
		path[++(path[0])] =cur;
	}
	return true;
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findLast(){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur = 1;
	path[0] =0;
	//  printf("%i.l = %i\n",cur, tree[cur].second); fflush(stdout);
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			path[++(path[0])] =target.tree[cur].second;
			break;
		}else cur = target.tree[cur].second | 1;
		path[++(path[0])] =cur;
	}
	return true;
}


LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findGE(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLT(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLT(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLT(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findLE(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGT(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGT(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGT(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findGT(const Key& key){
if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLE(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLE(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLE(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findLT(const Key& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGE(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGE(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGE(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::findEQ(const Key& key){this->findLE(key); if (!this->isValid()) return false; if (key != *(*this)) this->toInvalid(); return this->isValid();}
LFHTEMP	template<class O> bool RBTofDoom<Key,void>::Iterator::findEQ(const O& key){this->findLE(key); if (!this->isValid()) return false; if (ExOp::isNQ(key,*(*this))) this->toInvalid(); return this->isValid();}
LFHTEMP	template<class O> bool RBTofDoom<Key,void>::Iterator::findGE(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLT(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLT(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLT(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	template<class O> bool RBTofDoom<Key,void>::Iterator::findLE(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGT(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGT(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGT(target.tree[cur].first, key)) --(*this);
	return this->isValid();
}


LFHTEMP	template<class O> bool RBTofDoom<Key,void>::Iterator::findGT(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}

	unsigned int cur=1;
	path[0] = 0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isLE(target.tree[cur].first, key)) {cur = target.tree[cur].second; path[++(path[0])] = cur ; break;}
			else return true;
			}
			if (ExOp::isLE(target.tree[cur].first, key)) cur = target.tree[cur].second | 1;
			else cur = target.tree[cur].second & 0xFFFFFFFE;
			path[++(path[0])] = cur;
		}

	if (ExOp::isLE(target.tree[cur].first, key)) ++(*this);
	return this->isValid();
}
LFHTEMP	template<class O> bool RBTofDoom<Key,void>::Iterator::findLT(const O& key){
	if (target.size == 0) {path[0] = ExCo<unsigned int>::mkMaximum(); return false;}
	unsigned int cur=1;
	path[0] = 0;
	unsigned int height =0;
	while(target.tree[cur].second > 1){
		if (target.tree[cur].second > target.firstlone) {
			if (ExOp::isGE(target.tree[target.tree[cur].second].first, key)) break;
			else{cur = target.tree[cur].second; path[++(path[0])] = cur; return true;}
			}
			if (ExOp::isGE(target.tree[cur].first, key)) cur = target.tree[cur].second & 0xFFFFFFFE;
			else cur = target.tree[cur].second | 1;
			path[++(path[0])] = cur;
		}
	if (ExOp::isGE(target.tree[cur].first, key)) --(*this);
return this->isValid();}
LFHTEMP	template<class C> bool RBTofDoom<Key,void>::Iterator::findFirst(const C &query, Tuple<Key, 2u> &mima){
    if (target.size ==0) return false;
    mima[0] = min;
    mima[1] = max;
    path[0] =0;
return this->findNext(query,mima);}
LFHTEMP	template<class C> bool RBTofDoom<Key,void>::Iterator::findNext(const C &query, Tuple<Key, 2u> &mima){
    unsigned int cur;
    while( true ){
        while(!(SETCMP_DISJOINT & Query(mima[0],mima[1]))){
            if (target.tree[path[depth]].second - 2 > target.firstlone-2) {
                /*if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);
                if (tree[path[depth]].second > firstlone){
                    if (!(SETCMP_DISJOINT & Query(tree[tree[path[depth]].second].first,tree[tree[path[depth]].second].first))) f_out.push_back(tree[tree[path[depth]].second].first);
                }*/
                break;
            }
            cur = (path[0] == 0) ? 1 : path[path[0]];
            mima[1] = target.tree[cur].first;
            path[0]++;
            path[path[0]] = target.tree[cur].second & 0xFFFFFFFE;
        }


        if (path[0] < 2) {
            if (path[0] == 1) return false;
            mima[0] = target.tree[cur].first;
            mima[1] = max;
            cur = 0;
        }else{
            path[0]--;
            mima[1] = target.tree[(path[0] == 1) ? 1 : path[path[0]-1]].first;
            cur = path[path[0]];
        }
      //  if (!(SETCMP_DISJOINT & Query(tree[path[depth]].first,tree[path[depth]].first))) f_out.push_back(tree[path[depth]].first);



        // erase current parent in history, no longer needed
        path[depth] = target.tree[cur].second | 1;
    }
}


/*
LFHTEMP	void RBTofDoom<Key,void>::Iterator::toMaxWithinBrotherSubtree_routine(unsigned int height){
	if height == 0
}
LFHTEMP	void RBTofDoom<Key,void>::Iterator::toMinWithinBrotherSubtree_routine(unsigned int height){

}
*/

//prev LT
//next GT


LFHTEMP	typename RBTofDoom<Key,void>::Iterator& RBTofDoom<Key,void>::Iterator::operator=(const typename RBTofDoom<Key,void>::Iterator& other){
	target = other.target;
	memcpy(path,other.path, sizeof(int) *63);
return *this;}

LFHTEMP	bool RBTofDoom<Key,void>::Iterator::hsNext() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(false); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) return(false);
			branch = false;
		} else{
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
		}
		if (branch) {
			// no right child!
			cur = path[0]-1;
			if ((cur>0)&&((target.tree[path[cur]].second > target.firstlone)||((target.tree[path[cur]].second | 1) == path[cur+1]))) {cur--;
				while((cur>0)&&((target.tree[path[cur]].second | 1) == path[cur+1])) cur--;
			}

			if (cur==0) {
            if (((target.tree[1].second & 0xFFFFFFFE)!= path[1])||(target.tree[1].second > target.firstlone)) return false;
			}

		}
		return(true);
}
LFHTEMP	bool RBTofDoom<Key,void>::Iterator::hsPrev() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return(false); // invalid!
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) return(false);
			branch = false;
		} else if (path[path[0]] > target.firstlone) {return true;
		} else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
		}
		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			if (cur==0) {
                if ((target.tree[1].second | 1) != path[1]) return(false);
			}
		}
    return(true);
}

LFHTEMP	const Key& RBTofDoom<Key,void>::Iterator::Next() const{
    unsigned int j = Next_index();
    return target.tree[(j == 0) ? 1: j].first;
    }

LFHTEMP	unsigned int RBTofDoom<Key,void>::Iterator::Next_index() const{
        unsigned int cur;
		bool branch;
	//	if (path[0] == ExCo<unsigned int>::mkMaximum()) return(*this); // invalid!
		if (path[0] == 0) {
			if (target.tree[1].second == 0) return 0;
			branch = false;
		} else {
			branch = ((target.tree[path[path[0]]].second < 2)||(path[path[0]] > target.firstlone))  ;
		}

		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			if (cur==0) {
                if ((target.tree[1].second | 1) != path[1]) return(false);
			}
		}

		if (branch) {
			cur = path[0]-1;
			if ((cur>0)&&((target.tree[path[cur]].second > target.firstlone)||((target.tree[path[cur]].second | 1) == path[cur+1]))) {cur--;
				while((cur>0)&&((target.tree[path[cur]].second | 1) == path[cur+1])) cur--;
			}
			return (cur==0) ? 0 : path[cur];
		}else{
			// has right
			cur = (path[0] ==0) ? 1 : path[path[0]];
			if (target.tree[cur].second > target.firstlone) { cur = target.tree[cur].second;}
			else {cur = target.tree[cur].second | 1;
				while(target.tree[cur].second > 1){
					if (target.tree[cur].second > target.firstlone) break;
					else cur = target.tree[cur].second & 0xFFFFFFFE;
				}
			}
		}
	return cur;
}

LFHTEMP	const Key& RBTofDoom<Key,void>::Iterator::Prev() const{
    unsigned int j = Prev_index();
    return target.tree[(j == 0) ? 1: j].first;
    }

LFHTEMP	unsigned int RBTofDoom<Key,void>::Iterator::Prev_index() const{
		unsigned int cur;
		bool branch;
		if (path[0] == ExCo<unsigned int>::mkMaximum()) return 0;
		if (path[0] == 0) {
			if ((target.tree[1].second == 0)|| (target.tree[1].second > target.firstlone)) return 0;
			branch = false;
		} else if (path[path[0]] > target.firstlone) {
            return( path[0] == 1 ? 1 : path[path[0]-1] );
		}else{
			branch = (target.tree[path[path[0]]].second < 2)|| (target.tree[path[path[0]]].second > target.firstlone);
		}
    //    printf("%c brbr\n", branch ? 'Y' : 'N');
		if (branch) {
			// no left child!
			cur = path[0]-1;
			while((cur>0)&&((target.tree[path[cur]].second & 0xFFFFFFFE) == path[cur+1])) cur--;
			return (cur==0) ? (((target.tree[1].second | 1) == path[1]) ? 1 : 0 )  : path[cur] ;
		}else{
			// has left
			cur = (path[0] ==0) ? 1 : path[path[0]];
			cur = target.tree[cur].second & 0xFFFFFFFE;
				while(target.tree[cur].second > 1){
					if (target.tree[cur].second > target.firstlone) {
					    cur = target.tree[cur].second;
					    break;
				    }
				cur = target.tree[cur].second | 1;
			}
		}
	return(cur);
}


LFHTEMP unsigned int RBTofDoom<Key,void>::Iterator::get_index() const{return (path[0]) ? path[path[0]] : 1 ; }
LFHTEMP unsigned int RBTofDoom<Key,void>::Iterator::get_serial_index() const{
    uint32_t tmp = (path[0]) ? path[path[0]] : 1;
    return tmp; // return something else, to avoid confusion...
}

LFHTEMP bool RBTofDoom<Key,void>::Iterator::isValid() const{return ((path[0] &0xFFFF0000) == 0);}
LFHTEMP typename RBTofDoom<Key,void>::Iterator& RBTofDoom<Key,void>::Iterator::toInvalid(){path[0] = 0xFFFFFFFF; return *this;}

LFHTEMP const Key& RBTofDoom<Key,void>::Iterator::operator*() const{ return  target.tree[(path[0]) ? path[path[0]] : 1 ].first;  }
LFHTEMP const Key* RBTofDoom<Key,void>::Iterator::operator->() const{return  &(target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP Key* RBTofDoom<Key,void>::Iterator::ordering_preserving_change() const{return  &(target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP Key& RBTofDoom<Key,void>::Iterator::ordering_preserving_reference() const{return  (target.tree[(path[0]) ? path[path[0]] : 1 ].first);}
LFHTEMP bool RBTofDoom<Key,void>::Iterator::operator==(const RBTofDoom<Key,void>::Iterator& other)const{
    if ((path[0] &0xFFFF0000) != 0) return false;
    unsigned int i;
    for(i=0;i<path[0];i++) if (path[i] != other.path[i]) return false;
return true;}


LFHTEMP RBTofDoom<Key,void>::QueryIterator::QueryIterator(const RBTofDoom<Key,void>& _target): target(_target), path((_target.size == 0) ? NULL : new uint32_t[3 << _target.alloc_mag]){}
LFHTEMP RBTofDoom<Key,void>::QueryIterator::~QueryIterator(){delete[](path);}
LFHTEMP template <class QRY> bool RBTofDoom<Key,void>::QueryIterator::first(const QRY& query){
    if (path == NULL) return false;
    depth =0;
    path[0] = 1;
    tmin = target.getMin();
    tmax = target.getMax();

    if (target.size == 0) return false;

    while( true ){
        while(!(SETCMP_DISJOINT & query(tmin,tmax))){
            if (target.tree[path[depth]].second - 2 > target.firstlone-2) {

                if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;
                if (target.tree[path[depth]].second > target.firstlone){
                    if (!(SETCMP_DISJOINT & query(target.tree[target.tree[path[depth]].second].first,target.tree[target.tree[path[depth]].second].first))) {
                        path[depth+1] = target.tree[path[depth]].second;
                        depth++;
                        return true;
                    }

                }

                break;
            }
            tmax = target.tree[path[depth]].first;
            path[depth+1] = target.tree[path[depth]].second & 0xFFFFFFFE;
            depth++;
        }
        depth--;
        if (depth == 0xFFFFFFFF) break;
        if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;

        tmin = target.tree[path[depth]].first;
        if (depth == 0) tmax = target.getMax();
        else tmax = target.tree[path[depth-1]].first;

        // erase current parent in history, no longer needed
        path[depth] = target.tree[path[depth]].second | 1;
    }
    return false;
}
LFHTEMP template <class QRY> bool RBTofDoom<Key,void>::QueryIterator::next(const QRY& query){
    if (path[depth] > target.firstlone) {
        depth-= 2;
        if (depth == 0xFFFFFFFF) return false;
        if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;
    }else if (target.tree[path[depth]].second - 2 > target.firstlone-2) {
        if (target.tree[path[depth]].second > target.firstlone){
            if (!(SETCMP_DISJOINT & query(target.tree[target.tree[path[depth]].second].first,target.tree[target.tree[path[depth]].second].first))) {
                path[depth+1] = target.tree[path[depth]].second;
                depth++;
                return true;
            }
        }
        depth--;
        if (depth == 0xFFFFFFFF) return false;
        if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;
    }

    while( true ){
        // erase current parent in history, no longer needed
        tmin = target.tree[path[depth]].first;
        if (depth == 0) tmax = target.getMax();
        else tmax = target.tree[path[depth-1]].first;
        path[depth] = target.tree[path[depth]].second | 1;
        while(!(SETCMP_DISJOINT & query(tmin,tmax))){
            if (target.tree[path[depth]].second - 2 > target.firstlone-2) {
                if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;
                if (target.tree[path[depth]].second > target.firstlone){
                    if (!(SETCMP_DISJOINT & query(target.tree[target.tree[path[depth]].second].first,target.tree[target.tree[path[depth]].second].first))) {
                        path[depth+1] = target.tree[path[depth]].second;
                        depth++;
                        return true;
                    }
                }
                break;
            }
            tmax = target.tree[path[depth]].first;
            path[depth+1] = target.tree[path[depth]].second & 0xFFFFFFFE;
            depth++;
        }
        depth--;
        if (depth == 0xFFFFFFFF) break;
        if (!(SETCMP_DISJOINT & query(target.tree[path[depth]].first,target.tree[path[depth]].first))) return true;
    }
    return false;
}




#undef LFHTEMP
#define LFHTEMP template <class C, class F>

LFHTEMP bool RBDynaTree<C,F>::isRed(uint32_t w)const{return (w > firstlone)||((tree[w].second&1)&&(tree[w].second <= firstlone));}
LFHTEMP bool RBDynaTree<C,F>::isLeaf(uint32_t w)const{return((w > firstlone)||(tree[w].second <2));}
LFHTEMP ERRCODE RBDynaTree<C,F>::save(FILE* f) const{
    ERRCODE fout = tree.save(f);
    if (tree.getSize()>0){
        fwrite(&size,sizeof(uint32_t),1,f);
        fwrite(&firstlone,sizeof(uint32_t),1,f);
        fwrite(&min_offset,sizeof(uint32_t),1,f);
        fwrite(&max_offset,sizeof(uint32_t),1,f);
        if (ExCo<typename F::OUTTYPE>::IsPOD) fwrite(cached,sizeof(typename F::OUTTYPE),tree.getSize(),f);
        else for(uint32_t i=0;i<tree.getSize();i++) ExOp::save(cached[i],f);
    }
return fout;}
LFHTEMP ERRCODE RBDynaTree<C,F>::load(FILE* f){
    if (tree.getSize()>0) delete[](cached);
    ERRCODE fout = tree.load(f);
    if (tree.getSize()>0){
        fread(&size,sizeof(uint32_t),1,f);
        fread(&firstlone,sizeof(uint32_t),1,f);
        fread(&min_offset,sizeof(uint32_t),1,f);
        fread(&max_offset,sizeof(uint32_t),1,f);
        cached = new typename F::OUTTYPE[tree.getSize()];
        if (ExCo<typename F::OUTTYPE>::IsPOD) fwrite(cached,sizeof(typename F::OUTTYPE),tree.getSize(),f);
        else for(uint32_t i=0;i<tree.getSize();i++) ExOp::load(cached[i],f);
    }
return fout;}

 // end of namespace
