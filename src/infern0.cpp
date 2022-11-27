// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// dropoutrate * (log2(average) + 1)
// we only include RcppArmadillo.h which pulls Rcpp.h in for us

#define ARMA_DONT_PRINT_ERRORS

#include "RcppArmadillo.h"
#include "infern0.hpp"

using namespace LFHPrimitive;
using namespace arma;

template< uint32_t size> double batta_metric(const GaussElem<Tuple<double,size> > &l,const GaussElem<Tuple<double,size> > &r){return l.likelihoodratio_dist(r);}

Infern0scope::Infern0scope(){

    uint32_t ite = lfhstatic_ressources.find((uint32_t)LFHSTATIC_RESSOURCE_ANNOTATIONS);
    if (ite == 0xFFFFFFFF){
        lfhstatic_ressources[(uint32_t)LFHSTATIC_RESSOURCE_ANNOTATIONS] = (void*)&annot;
    }

}

Infern0scope::~Infern0scope(){}


uint32_t Infern0scope::getNBgenes()const{return annot.dicos[(int)INFRN0_AXES_GENE_NAMES].getSize();}
uint32_t Infern0scope::getNBcells()const{return annot.dicos[(int)INFRN0_AXES_CELL_NAMES].getSize();}


/*
void Infern0scope::fitNegativeBinomial(){
    int nb_features = rawdata.getNBrows();
    double shape_guess = 0.5f;
    cell_scale.toSize(rawdata.getNBcols()).toZero();
    gene_scale.toSize(rawdata.getNBrows()).toZero();
    r_features.setSize(rawdata.getNBrows());
    drop_features.toSize(rawdata.getNBrows()).toZero();
    drop_features_empirical.setSize(rawdata.getNBrows());

    Tuple<unsigned int, 2u> coor;
    for(coor[1] = 0 ; coor[1]< nb_features; coor[1]++) {r_features[coor[1]] = 0.5;}

    double tmp,sum, sum2;
    double stat[2];
    // first get initial guess by simple mean/var fits.
    WeightElem<double, 2u> hehe;
    ProgressBarPrint ppri(20);


    shape_guess = shape_guess / (1.0f - shape_guess);

    sum2 = 0.0f;
    sum =0;
    uint32_t rite;
    // initialize
    for(coor[0] = 0 ; coor[0]< cell_scale.getSize(); coor[0]++){
      for(rite = 0; rite < rawdata.data[coor[0]].getSize();rite++){
        drop_features[rawdata.data[coor[0]].deref_key(rite)] += 1.0f;
        sum2 += rawdata.data[coor[0]].deref(rite);
        gene_scale[rawdata.data[coor[0]].deref_key(rite)] += rawdata.data[coor[0]].deref(rite);
        cell_scale[coor[0]] += rawdata.data[coor[0]].deref(rite);
      }
      if (cell_scale[coor[0]] > sum) sum = cell_scale[coor[0]];
    }
    for(coor[1] = 0 ; coor[1]< nb_features; coor[1]++){
      drop_features_empirical[coor[1]] = drop_features[coor[1]] = 1.0f - (drop_features[coor[1]] / cell_scale.getSize());
    }
    for(coor[0] = 0 ; coor[0]< cell_scale.getSize(); coor[0]++) {
        cell_scale[coor[0]] /= 1.125f * sum; // 0 <= scale_cell <= 8/9
        cell_scale[coor[0]] = pow(cell_scale[coor[0]], 0.125f);
     }
    tmp = ((double)cell_scale.getSize()) / sum2;
    for(coor[1] = 0 ; coor[1]< nb_features; coor[1]++) {
        gene_scale[coor[1]] *= tmp; //= log(gene_scale[coor[1]]*tmp);
        r_features[coor[1]] = gene_scale[coor[1]] * 2.0f;
        if (r_features[coor[1]] < 0.1f) r_features[coor[1]] += 0.1f;
        gene_scale[coor[1]] = 0.5f;
    }

    CurvatureSearchScope css;

    double guess[2];
    double trguess[2];
    double deriv[2];
    double value, startval;
    double denum;
    int l;
    uint32_t total;
    uint32_t ite;
    char buffer[256];
    int maxouter =5;
    double total_ll;
    uint32_t total_cnt;
    uint32_t* indexes = new uint32_t[rawdata.getNBrows() > rawdata.getNBcols() ? rawdata.getNBrows()  : rawdata.getNBcols()];
    for(int outer =0; outer < maxouter;outer++){

        sprintf(buffer, "Normalize features (%i/%i)", outer + 1, maxouter);
        ppr.start(buffer);total_ll =0.0f; total_cnt =0;
        for(coor[0] = 0; coor[0]< nb_features; coor[0]++){ppr.update(((double)coor[0]) / nb_features);
            for(coor[1] = 0, total =0; coor[1]< cell_scale.getSize(); coor[1]++) {
                rite = rawdata.data[coor[1]].find(coor[0]);
                if ((rite == 0xFFFFFFFF)||(cell_scale[coor[1]] == 0)) continue;
                indexes[total++] = coor[1];
            }
            if (total == 0) continue;
            guess[0] = 2.0f * atanh(1.0f - 2.0f * gene_scale[coor[0]]);
            guess[1] = sqrt(r_features[coor[0]]- 0.05f);
            if (!ExOp::isValid(guess[0])) guess[0] =0.0f;
            if (!ExOp::isValid(guess[1])) guess[1] =0.0f;
            css.init(2, 0.001f);
            for(l=0;l<10 * (1 << outer);l++){
                trguess[0] = 1.0f / (1.0f + exp(guess[0]));
                trguess[1] = guess[1] * guess[1] + 0.0625f;
                value = -(lngamma(trguess[1]) * total);
                deriv[0] = 0;
                deriv[1] = -(d_lngamma_dx(trguess[1]) * total);


                for(ite =0 ; ite < total; ite++){
                    coor[1] = indexes[ite];
                    denum = 1.0f / (1.0f - trguess[0]* (1.0f - cell_scale[coor[1]]));
                    sum = -trguess[1] * log((1.0f - trguess[0]) * denum);
//                    if (!ExOp::isValid(log_exp_m1(sum))) printf("Nan for %e\n", sum);
//                    if (!ExOp::isValid(d_log_exp_m1_dx(sum))) printf("Nan dfor %e\n", sum);
                    value += lngamma(trguess[1]+rawdata(coor)) -log_exp_m1(sum)+ rawdata(coor) * log((trguess[0] * cell_scale[coor[1]]) * denum);
                    sum = d_log_exp_m1_dx(sum);
                    //if (!ExOp::isValid(sum)) printf("dx %e gave %e\n", -trguess[1] * log((1.0f - cell_scale[coor[1]]) * denum), sum);

//                    deriv[0] += -sum * trguess[1] * denum * cell_scale[coor[1]] + rawdata(coor) * denum * (1.0f - cell_scale[coor[1]]) / trguess[0];
                    deriv[0] -= (sum *  trguess[1] * cell_scale[coor[1]] / (1.0f - trguess[0]) - rawdata(coor) / trguess[0]) * denum;

                    deriv[1] += d_lngamma_dx(trguess[1]+rawdata(coor)) + sum * log((1.0f - trguess[0]) * denum);
                }
                deriv[1] *= 2.0f * guess[1];
                deriv[0] *= trguess[0] * (trguess[0] - 1.0f);

                if ((l==0)&&( (!ExOp::isValid(value)) || (!ExOp::isValid(deriv[0])) || (!ExOp::isValid(deriv[1])) )){
                  printf("error inbound!\n");
                  ExOp::show(value);
                  ExOp::show(deriv);
                  ExOp::show(guess);
                  ExOp::show(trguess);
                  exit(1);
                }


                if (l ==0) startval = value;
                if ((denum = css.updateAscent(value,guess,deriv)) < 0.0000000000001f) break;
                //if (l < 5) css.updateAscent(value,guess,deriv);else{if (l ==5) css.init(0.001f, 2);css.checkDerivative(value,guess,deriv);}
                //printf("%i: (%f,%f) %e   %e   %e %e\n", l, trguess[0], trguess[1], value, denum, guess[0], guess[1]);
            }
            css.wrFinalGuess(guess);
            total_ll += css.getLastValue();
            total_cnt += total;
            //printf("converged in %i %e\n", l, trguess[0]);
         //   printf("init guess %e vs final %e   %e\n", log(gene_scale[coor[0]] * 2.0f), guess[1], log(gene_scale[coor[0]] * 2.0f) - guess[1]);
            if (l == 1000){printf("%i: (%f,%f) %e     %f %e\n", l, trguess[0], trguess[1], value, denum, guess[0]);}
            sum = pow(0.5f, (double)outer);
            sum =0.0f;
            gene_scale[coor[0]] = 0.125f * sum  + (1.0f - 0.25f * sum) / (1.0f + exp(guess[0]));
            r_features[coor[0]] = guess[1] * guess[1] + 0.0625f;
            drop_features[coor[0]] = l;
            if (outer == 0) drop_features_empirical[coor[0]] = (css.getLastValue() - startval) / total;
        }ppr.finish();
        printf("Total LL is %e (%i)\n", total_ll/ total_cnt, total_cnt);
        if (outer+1 == maxouter) break;

        sprintf(buffer, "Normalize cells (%i/%i)", outer + 1, maxouter);
        ppr.start(buffer);total_ll =0.0f; total_cnt =0;
        for(coor[1] = 0; coor[1]< cell_scale.getSize(); coor[1]++) {ppr.update(((double)coor[1]) / cell_scale.getSize());
            guess[0] = 2.0f * atanh(1.0f - 2.0f * cell_scale[coor[1]]);
            if (!ExOp::isValid(guess[0])) guess[0] =0.0f;
            css.init(1,0.001f);
            for(rite = 0,total =0; rite< rawdata.data[coor[1]].getSize();  rite++) {
                if (gene_scale[rawdata.data[coor[1]].deref_key(rite)] < 0.0000001f) continue;
                indexes[total++] = rawdata.data[coor[1]].deref_key(rite);
            }
            if (total == 0) continue;
            for(l=0;l<10* (1 << outer);l++){
                trguess[0] = 1.0f / (1.0f + exp(guess[0]));// printf("guess %e\n", trguess[0]);
                value = 0;
                deriv[0] = 0;
                for(ite =0 ; ite < total; ite++){
                    coor[0] = indexes[ite];
                    denum = 1.0f / (1.0f - gene_scale[coor[0]] * (1.0f - trguess[0]));
                    sum = -r_features[coor[0]] * log((1.0f - gene_scale[coor[0]]) * denum);
                    value += -log_exp_m1(sum) +  rawdata(coor) * log (trguess[0] * gene_scale[coor[0]] * denum);
                    deriv[0] -= (d_log_exp_m1_dx(sum) * r_features[coor[0]] * gene_scale[coor[0]] + data(coor)  * (gene_scale[coor[0]] - 1.0f) / trguess[0]) * denum;
                }
                deriv[0] *= trguess[0] * (trguess[0] - 1.0f);

                if ((denum = css.updateAscent(value,guess,deriv)) < 0.0000000000001f) break;
            }
            css.wrFinalGuess(guess);

            total_ll += css.getLastValue();
            total_cnt += total;
            for(ite =0 ; ite < total; ite++){
                coor[0] = indexes[ite];
                total_ll += lngamma(r_features[coor[0]]+data(coor)) - lngamma(r_features[coor[0]]);
            }
            sum = pow(0.5f, (double)outer);
            cell_scale[coor[1]] = 0.125f * sum  + (1.0f - 0.25f * sum) / (1.0f + exp(guess[0]));
        }ppr.finish();
        printf("Total LL is %e (%i)\n", total_ll/ total_cnt, total_cnt);

        // transfer scale_cell to gene_scale
        sum =0;
        for(coor[1] = 0; coor[1]< cell_scale.getSize(); coor[1]++) {
            if (cell_scale[coor[1]] > sum) sum = cell_scale[coor[1]];
        }
        if (sum >= 0.99f) continue;
        //exit(0);
    }
    delete[](indexes);
}*/

void Infern0scope::cmpDeviation(){
  /*double scale,shape,prob_zero;
  SparseMatrix<uint32_t>::KeyIterator kite(rawdata);
  data.setDims(rawdata.getDims());
  if (kite.first()) do{
    scale = gene_scale[kite.getRow()] * cell_scale[kite.getCol()] / (1.0 - gene_scale[kite.getRow()] * (1.0f - cell_scale[kite.getCol()]));
    shape = r_features[kite.getRow()];
    prob_zero = pow(1.0- scale, shape);
    data(kite()) = LogPvalue_to_stdnorm((LogPvalue_NBdistrib_LH(*kite,shape,scale) - prob_zero) / (1.0f -prob_zero));
  }while(kite.next());*/
}

ERRCODE Infern0scope::wrDeviationMatrix(Rcpp::NumericMatrix &fout, const myHashmap<uint32_t> &rows, const myHashmap<uint32_t> &cols, int mode) const{
	fout = Rcpp::NumericMatrix(rows.getSize(), cols.getSize());
	Tuple< double > curproj_R, curproj_M;
	Tuple<double, 3u> coeff;
	Tuple<double, 2u> RM;
	uint32_t daite;
	bool use_cached = (mode == 0);
	if (auto itec = cols()) do{
		curproj_R = par_R * cell_state.getColumn(*itec);
		curproj_M = par_M * cell_state.getColumn(*itec);
		coeff = cell_scale[*itec];
		if (auto iter = rows()) do{
			if ((daite = data.data[*itec].find(*iter)) != 0xFFFFFFFF) {
				if (use_cached) {fout(iter.getOffset(), itec.getOffset()) = data.data[*itec].deref(daite); continue;}
				else daite = rawdata(*iter,*itec);
			}else daite =0;
			RM[0] = curproj_R.mkInnerProd(	gene_state.getColumn(*iter)) + coeff[0] + gene_scale[*iter][0];
			RM[1] = curproj_M.mkInnerProd(	gene_state.getColumn(*iter)) + coeff[1] + gene_scale[*iter][1];
			if ((!ExOp::isValid(RM[0]))||(!ExOp::isValid(RM[1]))) fout(iter.getOffset(), itec.getOffset()) = nan("");
			else switch(mode) {
				case 0: case 5: fout(iter.getOffset(), itec.getOffset()) = LogPvalue_to_stdnorm(LogPvalue_NBdistrib_exppara_LH_exact(daite,RM[0],-RM[1])); break;
				case 1: fout(iter.getOffset(), itec.getOffset()) = exp(RM[0]) * (1.0 + exp(-RM[1])) / (1.0 + exp(RM[1])); break; // mean
				case 2: fout(iter.getOffset(), itec.getOffset()) = exp(RM[0]) * (1.0 + exp(-RM[1])) * (1.0 + exp(-RM[1])) / (1.0 + exp(RM[1])); break; // var
				case 3: fout(iter.getOffset(), itec.getOffset()) = exp(RM[0]); break; // r param
				case 4: fout(iter.getOffset(), itec.getOffset()) = 1.0 / (1.0 + exp(RM[1])); break;
			}
		}while(iter++);
	}while(itec++);
return 0;}

Tuple<uint32_t,3> Infern0scope::testfunc(int value){Tuple<uint32_t,3> fout;
	fout[0]= 1;
	fout[1]= value;
	fout[2]= 42;
return fout;}

void Infern0scope::dothattest(){
//    using std::placeholders::_1;
//    std::function<void()> that = std::bind(&Infern0scope::cmpDeviation, this);
//    std::function<Tuple<uint32_t,3>(int)> thath = std::bind(&Infern0scope::testfunc, this,_1);
}


Tuple<double> Infern0scope::getBiodropoutProb(uint32_t col){Tuple<double> fout;



  fout.setSize(gene_scale.getSize());
  Tuple<double> curprojR = par_R * cell_state.getColumn(col);
  Tuple<double > curprojM = par_M * cell_state.getColumn(col);
  double rvb,mvb,rv,mv;
  if (ExOp::isValid(cell_scale[col][0])){
    rvb =cell_scale[col][0];
    mvb =cell_scale[col][1];
  }else{
    rvb =0.0f;
    mvb =0.0f;
  }

  for(uint32_t row =0;row<gene_scale.getSize();row++){
    rv = rvb + curprojR.mkInnerProd(gene_state.getColumn(row));
    mv = mvb + curprojM.mkInnerProd(gene_state.getColumn(row));
    if (ExOp::isValid(gene_scale[row][0])){
      rv += gene_scale[row][0];
      mv += gene_scale[row][1];
    }
    fout[row] = pow(1+exp(-mv) , -exp(rv));
  }

return fout;}
Tuple<double> Infern0scope::getBiodropoutProb(const Vector<uint32_t> &rows, uint32_t col){Tuple<double> fout;
  fout.setSize(rows.getSize());
  Tuple<double > curprojR = par_R * cell_state.getColumn(col);
  Tuple<double > curprojM = par_M * cell_state.getColumn(col);
  double rvb,mvb,rv,mv;
  if (ExOp::isValid(cell_scale[col][0])){
    rvb =cell_scale[col][0];
    mvb =cell_scale[col][1];
  }else{
    rvb =0.0f;
    mvb =0.0f;
  }
  for(uint32_t row =0;row<rows.getSize();row++){
    rv = rvb + curprojR.mkInnerProd(gene_state.getColumn(rows[row]));
    mv = mvb + curprojM.mkInnerProd(gene_state.getColumn(rows[row]));
    if (ExOp::isValid(gene_scale[rows[row]][0])){
      rv += gene_scale[rows[row]][0];
      mv += gene_scale[rows[row]][1];
    }
    fout[row] = pow(1+exp(-mv) , -exp(rv));
  }
return fout;}
double Infern0scope::getBiodropoutProb(uint32_t row,uint32_t col){double fout;
  Tuple<double > curprojR = par_R * cell_state.getColumn(col);
  Tuple<double > curprojM = par_M * cell_state.getColumn(col);
  double rvb,mvb,rv,mv;
  if (ExOp::isValid(cell_scale[col][0])){
    rvb =cell_scale[col][0];
    mvb =cell_scale[col][1];
  }else{
    rvb =0.0f;
    mvb =0.0f;
  }

    rv = rvb + curprojR.mkInnerProd(gene_state.getColumn(row));
    mv = mvb + curprojM.mkInnerProd(gene_state.getColumn(row));
    if (ExOp::isValid(gene_scale[row][0])){
      rv += gene_scale[row][0];
      mv += gene_scale[row][1];
    }
    fout = pow(1+exp(-mv) , -exp(rv));

return fout;}

ERRCODE Infern0scope::getGeneCellIndices(myHashmap<uint32_t> &rows, myHashmap<uint32_t>&cols, const std::vector<std::string> &genelist, const std::vector<std::string> &celllist, int& missinggenes, int& missingcells, bool printwarning) const{
	missinggenes=0;
	missingcells=0;
	int k;
	if ((genelist.size() == 0)||(genelist[0].length()==0)) for(int i=0; i < this->getNBgenes();i++) rows.addEntry(i);
	else for(int i=0;i<genelist.size();i++) {
		k = annot.dicos[INFRN0_AXES_GENE_NAMES].findEntry(genelist[i].c_str());
		if (k != 0xFFFFFFFF) rows.addEntry(k);
		else missinggenes++;
	}
	if ((celllist.size() == 0)||(celllist[0].length()==0)) for(int i=0; i < this->getNBcells();i++) cols.addEntry(i);
	else for(int i=0;i<celllist.size();i++) {
		k = annot.dicos[INFRN0_AXES_CELL_NAMES].findEntry(celllist[i].c_str());
		if (k != 0xFFFFFFFF) cols.addEntry(k);
		else missingcells++;
	}
	if (((missinggenes == genelist.size())&&(genelist.size() != 0))||((missingcells == celllist.size())&&(celllist.size() != 0))) {
		if (printwarning) printf("error: %i cells and %i genes queried were not found.\n", missingcells,missinggenes);
	return 1;}
	if ((printwarning)&&((missingcells != 0)||(missinggenes != 0))) printf("warning: %i cells and %i genes queried were not found.\n", missingcells,missinggenes);
return 0;}


void Infern0scope::wrSubMatrix(SparseMatrix<double> &target, const std::vector<std::string> genes, const std::vector<std::string> cells, Rcpp::CharacterVector *vecpair, int mincell){
  Vector<uint32_t> inds = annot.dicos[INFRN0_AXES_GENE_NAMES].findEntries(genes);
  Tuple<uint32_t> coverage; coverage.setSize(inds.getSize()); coverage.toZero();
  Vector<uint32_t> cinds;
  uint32_t i,j,k;

  if (cells.size() == 0){
    for(i=0;i<data.data.getSize();i++){
      for(j=0;j< inds.getSize();j++){
        if ((k = data.data[i].find(inds[j])) != 0xFFFFFFFF) coverage[j]++;
      }
    }
    if (vecpair) annot.getAxeNamesAl(vecpair[1], (int)INFRN0_AXES_CELL_NAMES);
  }else{
    cinds = annot.dicos[INFRN0_AXES_CELL_NAMES].findEntries(cells);
    for(i=0;i< cinds.getSize();i++){
      for(j=0;j< inds.getSize();j++){
        if ((k = data.data[cinds[i]].find(inds[j])) != 0xFFFFFFFF) coverage[j]++;
      }
    }
    if (vecpair) {
      for(i=0;i< cinds.getSize();i++) vecpair[1].push_back(annot.dicos[INFRN0_AXES_CELL_NAMES].entries[cinds[i]]);
    }
  }

  myHashmap<uint32_t, void> curind;
  for(i=0,j=0;j<inds.getSize();j++){
    if (coverage[j] < mincell){
      printf("warning, %s has too little cell coverage (filtered)\n", annot.dicos[INFRN0_AXES_GENE_NAMES].entries[inds[j]]);
    }else if (curind.find(inds[j]) != 0xFFFFFFFF){
      printf("warning, %s is found more than once\n", annot.dicos[INFRN0_AXES_GENE_NAMES].entries[inds[j]]);
    }else {inds[i++] = inds[j]; curind.addEntry(inds[j]);}
  }
  while( i < inds.getSize()) inds.pop_back();
  if (cells.size() == 0){
    target.setNBcols(data.getNBcols());
    for(i=0;i<data.data.getSize();i++){
      for(j=0;j< inds.getSize();j++){
        if ((k = data.data[i].find(inds[j])) != 0xFFFFFFFF) target.data[i][j] = data.data[i].deref(k);
      }
    }
  }else{
    target.setNBcols(cinds.getSize());
    for(i=0;i< target.getNBcols();i++){
      for(j=0;j< inds.getSize();j++){
        if ((k = data.data[cinds[i]].find(inds[j])) != 0xFFFFFFFF) target.data[i][j] = data.data[cinds[i]].deref(k);
      }
    }
    cinds.toMemfree();
  }
  if (vecpair){
    for(j=0;j< inds.getSize();j++) vecpair[0].push_back(annot.dicos[INFRN0_AXES_GENE_NAMES].entries[inds[j]]);
  }
}

void Infern0scope::findDiscriminatoryGenes(){
return;}

void Infern0scope::clusterHierarchical(bool do_cluster_cells, bool do_overwrite){
  Vector< GaussElem< Tuple<double, 0u > > > to_clust;

  uint32_t i,j,k,ind;
  Tuple<double, 0u> prior[3];
  prior[0].setSize(par_R.sizes[0] + par_R.sizes[1]).toZero();
  prior[1].setSize(par_R.sizes[0] + par_R.sizes[1]).toZero();
  prior[2].setSize(par_R.sizes[0] + par_R.sizes[1]);

  int nbval;
  Vector<uint32_t> perm;

  if (do_cluster_cells){
    nbval = this->populateProfile(to_clust, prior, true);
    double dist;


    for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++){
      prior[0][j] = (prior[1][j] - (prior[0][j] *prior[0][j] / nbval))/ nbval;
      prior[0][j] *= 0.0001;
    }
    printf("total prior:");
    prior[0].show();

    // adding prior, some oil for LL ratio
    for(i=0;i<cell_scale.getSize();i++){
      for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++) to_clust[i].cov.data[((j*(j+3))>>1)] += prior[0][j];
      to_clust[i].determinant = to_clust[i].getCovariance_biased().log_determinant();
    }

    for(i=1;i<cell_scale.getSize();i++){
       dist = (to_clust[i] + to_clust[i-1]).getCovariance_biased().log_determinant();
       //printf("dain %e %e and %e  -> %e?\n",to_clust[i-1].determinant,to_clust[i].determinant,dist,2.0 * dist - to_clust[i].determinant - to_clust[i-1].determinant);
       if ((2.0 * dist < to_clust[i].determinant + to_clust[i-1].determinant)||(!ExOp::isValid(2.0 * dist - to_clust[i].determinant + to_clust[i-1].determinant))){
          printf("GotNegative...\n");
          //to_clust[i-1].show();
          //to_clust[i].show();
          //(to_clust[i] + to_clust[i-1]).show();
       }
    }


    cell_hierarch.cluster_likelihood_ratio(to_clust);
    if (do_overwrite){cell_hierarch.wrAsPermutation(perm); cell_ordering.toMemmove(perm);}

  }else{
    nbval = this->populateProfile(to_clust, prior, false);
    double dist;


    for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++){
      prior[0][j] = (prior[1][j] - (prior[0][j] *prior[0][j] / nbval))/ nbval;
      prior[0][j] *= 0.0001;
    }
    printf("total prior:");
    prior[0].show();

    // adding prior, some oil for LL ratio
    for(i=0;i<gene_scale.getSize();i++){
      for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++) to_clust[i].cov.data[((j*(j+3))>>1)] += prior[0][j];
      to_clust[i].determinant = to_clust[i].getCovariance_biased().log_determinant();
    }

    for(i=1;i<gene_scale.getSize();i++){
       dist = (to_clust[i] + to_clust[i-1]).getCovariance_biased().log_determinant();
      // printf("dain %e %e and %e  -> %e?\n",to_clust[i-1].determinant,to_clust[i].determinant,dist,2.0 * dist - to_clust[i].determinant - to_clust[i-1].determinant);
       if ((2.0 * dist < to_clust[i].determinant + to_clust[i-1].determinant)||(!ExOp::isValid(2.0 * dist - to_clust[i].determinant + to_clust[i-1].determinant))){
          printf("GotNegative...\n");
       }
    }
    gene_hierarch.cluster_likelihood_ratio(to_clust);

    if (do_overwrite){gene_hierarch.wrAsPermutation(perm); gene_ordering.toMemmove(perm);}
  }

  //Rcpp::XPtr<Forest<double,2> > da_clust(new Forest<double,2>);
  //da_clust.cluster_likelihood_ratio(to_clust, batta_metric);
}

uint32_t Infern0scope::populateProfile(Vector< GaussElem< Tuple<double, 0u > > > &to_clust, Tuple<double, 0u> (&prior)[3], bool is_cell_profile){
  GaussElem< Tuple<double, 0u > > to_clust_input;
  GaussElem< Tuple<double, 0u > > to_clust_input_small;
  uint32_t i,j,ind,k;
  KeyElem<double, SparseTuple<double> > tmp;tmp.k = 1.0;
  ProgressBarPrint ppr(20);
  uint32_t  nbval=0u;
  if (is_cell_profile){
    to_clust_input_small.setSize(par_R.sizes[0]);
    ppr.start("Generating Profiles");
    for(i=0;i<data.data.getSize();i++){ ppr.update(((double)i)/cell_scale.getSize());

      to_clust_input_small.toZero();

      for(j=0;j< data.data[i].getSize();j++){
        tmp.d.toMemfree();
        if (auto ite = gene_state.data[data.data[i].deref_key(j)]()) do{
          tmp.d[ite()] = (*ite) * data.data[i].deref(j);
        }while(ite++);
        to_clust_input_small += tmp;
      }

      to_clust_input_small.setWeight(1.0);
      to_clust_input.setSize(par_R.sizes[1]+par_R.sizes[0]);
      to_clust_input.w = 1.0;
      to_clust_input.w2 = 0.0;
  /*
      varw.toZero();
      for(j=0;j< par_R.sizes[0];j++) {
        if (to_clust_input_small.cov[j] > 0.0) varw += WeightElem<double, 1>(log(to_clust_input_small.cov[j]));
      }*/

      for(j=0;j< par_R.sizes[0];j++) to_clust_input.mean[j] = to_clust_input_small.mean[j];
      for(j=0;j< par_R.sizes[1];j++) to_clust_input.mean[j+par_R.sizes[0]] = ((ind = cell_state.data[i].find(j)) == 0xFFFFFFFF) ? 0.0 : cell_state.data[i].deref(ind);
      k = (par_R.sizes[0] * (par_R.sizes[0] +1)) >> 1;
      for(j=0;j<k;j++) to_clust_input.cov.data[j] = to_clust_input_small.cov.data[j];
      for(ind=par_R.sizes[0];ind<par_R.sizes[1]+par_R.sizes[0];ind++){
        for(j=0;j<ind;j++) to_clust_input.cov.data[k++] = to_clust_input.mean[j] * to_clust_input.mean[ind];
        to_clust_input.cov.data[k++] = to_clust_input.mean[ind] * to_clust_input.mean[ind];
      }

      for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++) prior[2][j] = to_clust_input.cov.data[((j*(j+3))>>1)];

      if ((ExOp::isValid(to_clust_input.mean[0]))&&(ExOp::isValid(prior[2]))){
        prior[0] += to_clust_input.mean;
        prior[1] += prior[2];
        nbval++;
      }

      to_clust.push_back().toMemmove(to_clust_input);
    }ppr.finish();
  }else{
    to_clust_input_small.setSize(par_R.sizes[1]);

    ppr.start("Generating Profiles");
    SparseMatrix<double> trdata = data.mkTranspose();
    nbval=0;

    for(i=0;i<trdata.data.getSize();i++){ ppr.update(((double)i)/gene_scale.getSize());

      to_clust_input_small.toZero();

      for(j=0;j< trdata.data[i].getSize();j++){
        tmp.d.toMemfree();
        if (auto ite = cell_state.data[trdata.data[i].deref_key(j)]()) do{
          tmp.d[ite()] = (*ite) * trdata.data[i].deref(j);
        }while(ite++);
        to_clust_input_small += tmp;
      }

      to_clust_input_small.setWeight(1.0);
      to_clust_input.setSize(par_R.sizes[1]+par_R.sizes[0]);
      to_clust_input.w = 1.0;
      to_clust_input.w2 = 0.0;

      for(j=0;j< par_R.sizes[1];j++) to_clust_input.mean[j] = to_clust_input_small.mean[j];
      for(j=0;j< par_R.sizes[0];j++) to_clust_input.mean[j+par_R.sizes[1]] = ((ind = gene_state.data[i].find(j)) == 0xFFFFFFFF) ? 0.0 : gene_state.data[i].deref(ind);
      k = (par_R.sizes[1] * (par_R.sizes[1] +1)) >> 1;
      for(j=0;j<k;j++) to_clust_input.cov.data[j] = to_clust_input_small.cov.data[j];
      for(ind=par_R.sizes[1];ind<par_R.sizes[1]+par_R.sizes[0];ind++){
        for(j=0;j<ind;j++) to_clust_input.cov.data[k++] = to_clust_input.mean[j] * to_clust_input.mean[ind];
        to_clust_input.cov.data[k++] = to_clust_input.mean[ind] * to_clust_input.mean[ind];
      }

      for(j=0;j< par_R.sizes[0]+par_R.sizes[1];j++) prior[2][j] = to_clust_input.cov.data[((j*(j+3))>>1)];
      if ((ExOp::isValid(to_clust_input.mean[0]))&&(ExOp::isValid(prior[2]))){
        prior[0] += to_clust_input.mean;
        prior[1] += prior[2];
        nbval++;
      }
      to_clust.push_back().toMemmove(to_clust_input);
    }ppr.finish();
  }
return nbval;}

void Infern0scope::saveHierarchical(const char* path){
  char buffer[1024]; int p_l = strlen(path);
  memcpy(buffer,path,p_l);
  uint32_t i,j,k,ind;
  FILE* fcdt;
  Vector<char*> header;
  Vector< GaussElem< Tuple<double, 0u > > > to_clust;
  KeyElem<double, SparseTuple<double> > tmp;tmp.k = 1.0;

  Tuple<double, 0u> prior[3];
  prior[0].setSize(par_R.sizes[0] + par_R.sizes[1]).toZero();
  prior[1].setSize(par_R.sizes[0] + par_R.sizes[1]).toZero();
  prior[2].setSize(par_R.sizes[0] + par_R.sizes[1]);
  int nbval;
  if (cell_hierarch.getSize() > 0){
    strcpy(buffer+p_l,"_cell.gtr");
    cell_hierarch.saveGTRfile(buffer, "GENE");
    strcpy(buffer+p_l,"_cell.cdt");
    fcdt = fopen(buffer,"w+");

    fprintf(fcdt, "GID\tID\tNAME\tNB_P\tNB_R\tGWEIGHT");
    for(i=0;i< par_R.sizes[0];i++) fprintf(fcdt, "\tGDev%i",i);
    for(i=0;i< par_R.sizes[1];i++) fprintf(fcdt, "\tCSta%i",i);
    for(i=0;i< par_R.sizes[0];i++) fprintf(fcdt, "\tGDevV%i",i);
    for(i=0;i< par_R.sizes[1];i++) fprintf(fcdt, "\tCStaV%i",i);
    fprintf(fcdt, "\nEWEIGHT\t\t\t\t\t\t\t\t"); for(int i=0;i<(par_R.sizes[0]+ par_R.sizes[1] * 2);i++) fprintf(fcdt, "\t1.0");
    fprintf(fcdt, "\n");

    header.setSize(cell_scale.getSize());
    for(i=0;i<cell_scale.getSize();i++){
        sprintf(buffer, "GENE%iX\t%s\t%s\t%f\t%f\t1.0", i+1, annot.dicos[(uint32_t)INFRN0_AXES_CELL_NAMES].entries[i], annot.dicos[(uint32_t)INFRN0_AXES_CELL_NAMES].entries[i], cell_scale[i][0], cell_scale[i][1]);
        int j =strlen(buffer)+1;
        header[i] = new char[j]; memcpy(header[i],buffer,j);
    }

    // rere compute... (sigh)
    this->populateProfile(to_clust, prior, true);
    // rere compute end ... (sigh)

    cell_hierarch.saveCDTfile_W(fcdt, header, to_clust);
    for(i=0;i<header.getSize();i++) delete[](header[i]);
    fclose(fcdt);
    to_clust.toMemfree();
  }

  if (gene_hierarch.getSize() > 0){
    memcpy(buffer,path,p_l);
    strcpy(buffer+p_l,"_gene.gtr");
    gene_hierarch.saveGTRfile(buffer, "GENE");
    strcpy(buffer+p_l,"_gene.cdt");
    fcdt = fopen(buffer,"w+");

    fprintf(fcdt, "GID\tID\tNAME\tNB_P\tNB_R\tGWEIGHT");
    for(i=0;i< par_R.sizes[0];i++) fprintf(fcdt, "\tCDev%i",i);
    for(i=0;i< par_R.sizes[1];i++) fprintf(fcdt, "\tGSta%i",i);
    for(i=0;i< par_R.sizes[0];i++) fprintf(fcdt, "\tCDevV%i",i);
    for(i=0;i< par_R.sizes[1];i++) fprintf(fcdt, "\tGStaV%i",i);
    fprintf(fcdt, "\nEWEIGHT\t\t\t\t\t\t\t\t"); for(int i=0;i<(par_R.sizes[0]+ par_R.sizes[1] * 2);i++) fprintf(fcdt, "\t1.0");
    fprintf(fcdt, "\n");

    header.setSize(gene_scale.getSize());
    for(i=0;i<gene_scale.getSize();i++){
        sprintf(buffer, "GENE%iX\t%s\t%s\t%f\t%f\t1.0", i+1, annot.dicos[(uint32_t)INFRN0_AXES_GENE_NAMES].entries[i], annot.dicos[(uint32_t)INFRN0_AXES_GENE_NAMES].entries[i], gene_scale[i][0], gene_scale[i][1]);
        int j =strlen(buffer)+1;
        header[i] = new char[j]; memcpy(header[i],buffer,j);
    }
    // rere compute... (sigh)
    this->populateProfile(to_clust, prior, false);
    // rere compute end ... (sigh)

    gene_hierarch.saveCDTfile_W(fcdt, header, to_clust);
    for(i=0;i<header.getSize();i++) delete[](header[i]);
    fclose(fcdt);
    to_clust.toMemfree();
  }
}
ERRCODE Infern0scope::save(FILE*f)const{
    ERRCODE fout = 0;
    fout |= rawdata.save(f);
    fout |= gene_scale.save(f);
    fout |= cell_scale.save(f);
    fout |= par_R.save(f);
    fout |= par_M.save(f);
    fout |= par_D.save(f);
    fout |= par_S.save(f);

    fout |= gene_state.save(f);
    fout |= cell_state.save(f);

    fout |= data.save(f);
    fout |= pval_data.save(f);
    fout |= annot.save(f);

    fout |= gene_hierarch.save(f);
    fout |= cell_hierarch.save(f);
    fout |= gene_ordering.save(f);
    fout |= cell_ordering.save(f);

    fout |= gene_clusterID.save(f);
    fout |= cell_clusterID.save(f);

    fout |= metadata.save(f);

    return fout;
}
ERRCODE Infern0scope::load(FILE*f){
    ERRCODE fout = 0;
    fout |= rawdata.load(f);
    fout |= gene_scale.load(f);
    fout |= cell_scale.load(f);
    fout |= par_R.load(f);
    fout |= par_M.load(f);
    fout |= par_D.load(f);
    fout |= par_S.load(f);

    fout |= gene_state.load(f);
    fout |= cell_state.load(f);

    fout |= data.load(f);
    fout |= pval_data.load(f);
    fout |= annot.load(f);

    fout |= gene_hierarch.load(f);

    fout |= cell_hierarch.load(f);

    fout |= gene_ordering.load(f);

    fout |= cell_ordering.load(f);


    fout |= gene_clusterID.load(f);

    fout |= cell_clusterID.load(f);
    if (fout |= metadata.load(f)) return fout;


return fout;}

// return derivation and weight for a dropout
Tuple<double, 2u> Infern0scope::DropoutDerivation::operator()(const Tuple<uint32_t,2u>&) const{ Tuple<double,2u> fout;


return fout;}

void Infern0scope::runTaskPartcorrel(TMatrix<double,0u,0u> &out_correl, Vector<uint32_t> &out_genelist, const Vector<uint32_t> &gind, const Vector<uint32_t> &cind, uint32_t nb_threads, uint32_t indimax, uint32_t mincell,const double* minthr) const{
	TaskPartcorrel task(*this,out_correl,true);
	uint32_t i;
	task.ranges.setSize(nb_threads+1);
  task.part.setSizes(gene_scale.getSize(), gind.size());
	for(i=0;i<= nb_threads;i++){
		task.ranges[i] = (i * gene_scale.getSize()) / nb_threads;
	}
	for(i=0;i<gind.getSize();i++) task.mapping[gind[i]] =i;
	if (gind.getSize() != task.mapping.getSize()){
		printf("Detected duplicated entries in input, abort."); return;
	}
	task.topass.setSize(gene_scale.getSize());
	task.asad.setSize(gene_scale.getSize());

	task.trdata.setNBcols(gene_scale.getSize());
	task.mincell = mincell;
	task.minthr = minthr;

	if (cind.getSize() == 0){
    for(uint32_t c=0;c < cell_scale.getSize();c++){
       for(uint32_t ite=0;ite< data.data[c].getSize();ite++) task.trdata.data[data.data[c].deref_key(ite)][c] = data.data[c].deref(ite);
    }
	}else{
		for(uint32_t c=0;c < cind.getSize();c++){
			for(uint32_t ite=0;ite< data.data[cind[c]].getSize();ite++) task.trdata.data[data.data[cind[c]].deref_key(ite)][c] = data.data[cind[c]].deref(ite);
		}
	}
  	ThreadBase tb;

	tb.toSize(nb_threads-1);
	printf("\n\n\n\n\n\n\n\n\n\n\nstart runrun!\n");
	for(i=nb_threads-1;i>0;i--) tb.startEvent(&task, i);
        tb.startEvent_ThenWait(&task, i);
	printf("done runrun\n"); fflush(stdout);

	// very sad
	for(i=0;i<gene_scale.getSize();i++){
		if (task.asad[i].d != 0xFFFFFFFF) task.heap_correl <<= task.asad[i];
	}
	// very sad
	KeyElem<double, uint32_t> heapin;
	Tuple<uint32_t> coucount; coucount.setSize(gind.getSize()).toZero();
	for(i=0;i<out_genelist.getSize();i++){
		if (task.heap_correl.isEmpty()) break;
		heapin = task.heap_correl.pop();
		if (task.topass[heapin.d] >= gind.getSize()){
			printf("warning, invalid top associated %i to %i\n",task.topass[heapin.d],heapin.d);
			i--;
			continue;
		}
		if ((indimax != 0u)&&(coucount[task.topass[heapin.d]]++ >= indimax)) {i--;  continue;}
		printf("%e with %i\n", heapin.k, task.topass[heapin.d]);
		out_genelist[i] = heapin.d;
	}
	while (i < out_genelist.getSize()) out_genelist.pop_back();

return;}
uint32_t Infern0scope::TaskPartcorrel::operator()(uint32_t threadID){
	uint32_t inds[2];
	uint32_t ite;
	GaussElem< Tuple<double, 2u>, 0u > correlbuf;
	Tuple<double, 2u> correlinput;
	Tuple<uint32_t, 2u> coroutcoor;
	KeyElem<double, uint32_t> toins;
	uint32_t best;
	double bestv;
	ProgressBarPrint progbar(20);
	if (threadID == 0) progbar.start("Computing Correlation");

	for(toins.d = ranges[threadID];toins.d < ranges[threadID+1];toins.d++){
		if (threadID == 0) progbar.update(((double)toins.d)/ranges[1]);
		if (mapping.find(toins.d) != 0xFFFFFFFF) continue;
		toins.k = 0.0;
		coroutcoor[0] = toins.d;
		bestv = -1.0;
		for(uint32_t x=0;x<mapping.getSize();x++){
			coroutcoor[1] = mapping.deref_key(x);
			if (trdata.data[coroutcoor[1]].getSize() < trdata.data[coroutcoor[1]].getSize()) {inds[0] = coroutcoor[1] ; inds[1] = toins.d;}
			else {inds[1] = coroutcoor[1]; inds[0] = toins.d;}
			correlbuf.toZero();
			for(uint32_t c=0;c<trdata.data[inds[0]].getSize();c++){
				if ((ite = trdata.data[inds[1]].find(trdata.data[inds[0]].deref_key(c))) != 0xFFFFFFFF){
					correlinput[0] = trdata.data[inds[0]].deref(c);
					correlinput[1] = trdata.data[inds[1]].deref(ite);
					correlbuf += GaussElem< Tuple<double, 2u>, 0u >(correlinput);
				}
			}
			Trianglix<double, 2u> cov = correlbuf.getCovariance_biased();
			double tmp = cov.data[1] / sqrt(cov.data[0] * cov.data[2]);
			coroutcoor[1] =x;
			if ((((int)correlbuf.w) >= mincell)&&(ExOp::isValid(tmp))) {
				part(coroutcoor) = tmp;
				tmp *= tmp;
				if (((int)correlbuf.w) >= 3){
					if (tmp > minthr[((int)correlbuf.w)-3]){ // needs to be higher than threshold to count...
						toins.k -= tmp - minthr[((int)correlbuf.w)-3];
						if (bestv < tmp) {best = x;bestv=tmp;}
					}
				}
			} else  part(coroutcoor) = 0.0;
		}
		if (bestv != -1.0) {
			topass[toins.d] = best;
			//async <<= toins;
			asad[toins.d] = toins;
		}else asad[toins.d].d = 0xFFFFFFFF;

	}
  if (threadID == 0) progbar.finish();

return 0;}


void Infern0scope::getHypergeometricDropout(vector<double> &logpval, vector<uint32_t> &maxcluster){
  /*uint32_t i,j;
  if (!annot.hasDico((int)INFRN0_AXES_CELL_CLUSTER_NAMES)) return;
  uint32_t nbclass = annot.dicos[(int)INFRN0_AXES_CELL_CLUSTER_NAMES].entries.getSize();
  if (!annot.hasDico((int)INFRN0_AXES_GENE_NAMES)) return;
  uint32_t nbgene = annot.dicos[(int)INFRN0_AXES_GENE_NAMES].entries.getSize();

  Tuple<uint32_t> clustersizes; clustersizes.toSize(nbclass).toZero();
  Tuple<uint32_t> cellcov; cellcov.toSize(nbclass+1).toZero();

  for(uint32_t j=0;j<cell_scale.getSize();j++) clustersizes[cell_clusterID[i]]++;

  for(uint32_t i=0;i<nbgene;i++){
    for(uint32_t j=0;j<cell_scale.getSize();j++) if (data.data.find(i) != 0xFFFFFFFF) {cellcov[cell_clusterID[i]]++; cellcov[nbclass]++;}


    cellcov.toZero();
  }*/
}


// predecated

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// predecated
// [[Rcpp::export]]
Rcpp::List InferN0_identifyL0Network(const Rcpp::NumericMatrix & x, const double & y) {
	//printf("%i and %i\n", x.n_cols, x.n_rows);

	arma::mat xmat = Rcpp::as<arma::mat>(x);

	SparseTrianglix<double> hehe;
	Trianglix<double> target; target.setSize(xmat.n_cols); target.toZero();

	Trianglix< Tuple<double, 4u> > target_w; target_w.setSize(target.getSize()); target_w.toZero();
	Tuple< Trianglix< Tuple<double, 4u> >, 9u> cross_target_w;
	Tuple<double, 4u> input_w; input_w[0] = 1.0f;

	for(int i = 0 ; i < 8;i++){
	cross_target_w[i].setSize(target.getSize());
	cross_target_w[i].toZero();
	}

	Tuple<double> input; input.setSize(target.getSize());


	unsigned int ite,i,j,k;
	double tmpdouble;
	for(ite = 0 ; ite < xmat.n_rows;ite++){
	for(j =0,k=0; j < target.getSize();j++,k++){
	input[j] = xmat.at(ite,j);
	input_w[1] = input[j];
	if (input[j] == 0.0f){k+=j; continue;}
	for(i=0;i<j;i++,k++){
	if (input[i] != 0.0f) {
	input_w[2] = input[i];
	input_w[3] = input_w[2] * input_w[1];
	target_w.data[k]+= input_w;
	cross_target_w[ite & 7].data[k]+= input_w;
	}
	}
	input_w[3] = input_w[1] * input_w[1];
	target_w.data[k]+= input_w;
	cross_target_w[ite & 7].data[k]+= input_w;
	}
	}


	Tuple< Trianglix<double>, 0u> cross; cross.setSize(30);
	arma::colvec mu_out = arma::colvec(target.getSize());
	arma::colvec std_out = arma::colvec(target.getSize());


	for(j =0,k=0; j < input.getSize();j++,k++){
	for(i=0;i<j;i++,k++){
	target.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][2] / target_w.data[k][0]) / target_w.data[k][0];
	}
	target.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][1] / target_w.data[k][0]) / target_w.data[k][0];
	mu_out.at(j) = target_w.data[k][1] / target_w.data[k][0];
	std_out.at(j) = sqrt(target.data[k]);
	}

	int indexes[9];
	for(ite = 0; ite< 8;ite++) indexes[ite] = ite;

	for(ite = 0; ite< cross.getSize();ite+=2){ // up to 70
	//  printf("%i%i%i%i,%i%i%i%i\n", indexes[0], indexes[1], indexes[2], indexes[3], indexes[4], indexes[5], indexes[6], indexes[7]);
	cross[ite].setSize(target.getSize());
	cross[ite|1].setSize(target.getSize());

	cross_target_w[8] = cross_target_w[indexes[0]];
	cross_target_w[8]+= cross_target_w[indexes[1]];
	cross_target_w[8]+= cross_target_w[indexes[2]];
	cross_target_w[8]+= cross_target_w[indexes[3]];
	//  cross_target_w[8].show();
	for(j =0,k=0; j < input.getSize();j++,k++){
	for(i=0;i<j;i++,k++){
	cross[ite].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][2]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	cross[ite].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][1]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	//   cross[ite].show();
	cross_target_w[8] = cross_target_w[indexes[4]];
	cross_target_w[8]+= cross_target_w[indexes[5]];
	cross_target_w[8]+= cross_target_w[indexes[6]];
	cross_target_w[8]+= cross_target_w[indexes[7]];
	//   cross_target_w[8].show();
	for(j =0,k=0; j < input.getSize();j++,k++){
	for(i=0;i<j;i++,k++){
	cross[ite|1].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][2]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	cross[ite|1].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][1]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	//    cross[ite|1].show();
	if ((ite & 2) == 0){
	indexes[8] = indexes[7];
	indexes[7] = indexes[3];
	indexes[3] = indexes[8];
	}else switch(ite >> 2){
	case 0: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 1: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;
	case 2: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 3: indexes[8] = indexes[4];  indexes[4] = indexes[5]; indexes[5] = indexes[6]; indexes[6] = indexes[7]; indexes[7] = indexes[8]; indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8];break;
	case 4: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 5: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;
	case 6: indexes[8] = indexes[4];  indexes[4] = indexes[5]; indexes[5] = indexes[6]; indexes[6] = indexes[7]; indexes[7] = indexes[8]; indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8];break;
	case 7: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 8: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;

	//  case 3:
	//  case 1: indexes[1] = 5; indexes[5] = 1; break;
	//  case 2: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	}

	}

	//   hehe.setSize(target.getSize());
	//   hehe.toOne();

	Vector< KeyElem<Tuple<uint32_t, 2u>, Tuple<double, 11u> > > details;

	hehe.searchRegul(target, cross, details, y);

	arma::mat crdmat = arma::mat(details.getSize(),2);
	for(int i =0;i<details.getSize();i++){
	crdmat.at(i,0) = details[i].k[0];
	crdmat.at(i,1) = details[i].k[1];
	}

	Rcpp::NumericMatrix llincrmat(details.getSize(),11);

	for(int i =0;i<details.getSize();i++){
	for(int j=0;j<11;j++) llincrmat(i,j) = details[i].d[j];
	}



	Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("LL","CrX_deter","CrX_trace","CrX_LLgradient","CrX_LLstdev","Error","CrX_error","TPrate","FPrate", "InsIntoSize", "TimeMilli");



	Rcpp::NumericMatrix fout;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());
	hehe.wrMatrix(fout);
	Rcpp::NumericMatrix dacov;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());
	target.wrMatrix(dacov);
	Rcpp::colnames(dacov) = Rcpp::colnames(x);
	Rcpp::rownames(dacov) = Rcpp::colnames(x);
	Rcpp::colnames(fout) = Rcpp::colnames(x);
	Rcpp::rownames(fout) = Rcpp::colnames(x);

	return Rcpp::List::create(Rcpp::Named("output") = fout,
	Rcpp::Named("covar") = dacov,
	Rcpp::Named("EdgeList") = crdmat,
	Rcpp::Named("LL_incrs") = llincrmat,
	Rcpp::Named("mean") = mu_out,
	Rcpp::Named("stddev") = std_out);
}

// report on precision and false positive rate instead of error in LL_incr
//
// [[Rcpp::export]]
Rcpp::List InferN0_identifyL0NetworkGold(const Rcpp::NumericMatrix & x, const arma::mat & y, const double & z) {
	//printf("%i and %i\n", x.n_cols, x.n_rows);

	arma::mat xmat = Rcpp::as<arma::mat>(x);
	SparseTrianglix<double> hehe;
	Trianglix<double> target; target.setSize(xmat.n_cols); target.toZero();

	Trianglix< Tuple<double, 4u> > target_w; target_w.setSize(target.getSize()); target_w.toZero();
	Tuple< Trianglix< Tuple<double, 4u> >, 9u> cross_target_w;
	Tuple<double, 4u> input_w; input_w[0] = 1.0f;

	for(int i = 0 ; i < 8;i++){
	cross_target_w[i].setSize(target.getSize());
	cross_target_w[i].toZero();
	}

	Tuple<double> input; input.setSize(target.getSize());


	unsigned int ite,i,j,k;
	double tmpdouble;
	for(ite = 0 ; ite < xmat.n_rows;ite++){
	  for(j =0,k=0; j < target.getSize();j++,k++){
	    input[j] = xmat.at(ite,j);
	    input_w[1] = input[j];
	    if (input[j] == 0.0f){k+=j; continue;}
      for(i=0;i<j;i++,k++){
        if (input[i] != 0.0f) {
	        input_w[2] = input[i];
	        input_w[3] = input_w[2] * input_w[1];
	        target_w.data[k]+= input_w;
	        cross_target_w[ite & 7].data[k]+= input_w;
	      }
	    }
	    input_w[3] = input_w[1] * input_w[1];
	    target_w.data[k]+= input_w;
	    cross_target_w[ite & 7].data[k]+= input_w;
	  }
	}



	Tuple< Trianglix<double>, 0u> cross; cross.setSize(30);
	arma::colvec mu_out = arma::colvec(target.getSize());
	arma::colvec std_out = arma::colvec(target.getSize());


	for(j =0,k=0; j < input.getSize();j++,k++){
	  for(i=0;i<j;i++,k++){
	    target.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][2] / target_w.data[k][0]) / target_w.data[k][0];
	  }
	  target.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][1] / target_w.data[k][0]) / target_w.data[k][0];
	  mu_out.at(j) = target_w.data[k][1] / target_w.data[k][0];
	  std_out.at(j) = sqrt(target.data[k]);
	}

	int indexes[9];
	for(ite = 0; ite< 8;ite++) indexes[ite] = ite;

	for(ite = 0; ite< cross.getSize();ite+=2){ // up to 70
	//  printf("%i%i%i%i,%i%i%i%i\n", indexes[0], indexes[1], indexes[2], indexes[3], indexes[4], indexes[5], indexes[6], indexes[7]);
	cross[ite].setSize(target.getSize());
	cross[ite|1].setSize(target.getSize());

	cross_target_w[8] = cross_target_w[indexes[0]];
	cross_target_w[8]+= cross_target_w[indexes[1]];
	cross_target_w[8]+= cross_target_w[indexes[2]];
	cross_target_w[8]+= cross_target_w[indexes[3]];
	//  cross_target_w[8].show();
	for(j =0,k=0; j < input.getSize();j++,k++){
	for(i=0;i<j;i++,k++){
	cross[ite].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][2]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	cross[ite].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][1]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	//   cross[ite].show();
	cross_target_w[8] = cross_target_w[indexes[4]];
	cross_target_w[8]+= cross_target_w[indexes[5]];
	cross_target_w[8]+= cross_target_w[indexes[6]];
	cross_target_w[8]+= cross_target_w[indexes[7]];
	//   cross_target_w[8].show();
	for(j =0,k=0; j < input.getSize();j++,k++){
	for(i=0;i<j;i++,k++){
	cross[ite|1].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][2]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	cross[ite|1].data[k] = (cross_target_w[8].data[k][3] - cross_target_w[8].data[k][1]*cross_target_w[8].data[k][1]/cross_target_w[8].data[k][0]) / cross_target_w[8].data[k][0];
	}
	//    cross[ite|1].show();
	if ((ite & 2) == 0){
	indexes[8] = indexes[7];
	indexes[7] = indexes[3];
	indexes[3] = indexes[8];
	}else switch(ite >> 2){
	case 0: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 1: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;
	case 2: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 3: indexes[8] = indexes[4];  indexes[4] = indexes[5]; indexes[5] = indexes[6]; indexes[6] = indexes[7]; indexes[7] = indexes[8]; indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8];break;
	case 4: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 5: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;
	case 6: indexes[8] = indexes[4];  indexes[4] = indexes[5]; indexes[5] = indexes[6]; indexes[6] = indexes[7]; indexes[7] = indexes[8]; indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8];break;
	case 7: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	case 8: indexes[8] = indexes[1];  indexes[1] = indexes[5]; indexes[5] = indexes[8]; break;

	//  case 3:
	//  case 1: indexes[1] = 5; indexes[5] = 1; break;
	//  case 2: indexes[8] = indexes[2];  indexes[2] = indexes[6]; indexes[6] = indexes[8]; break;
	}

	}

	//   hehe.setSize(target.getSize());
	//   hehe.toOne();

	Vector< KeyElem<Tuple<uint32_t, 2u>, Tuple<double, 11u> > > details;

	Trianglix<char> gold; gold.setSize(target.getSize()); gold.toZero();
	for(j =0,k=0; j < target.getSize();j++,k++){
	  for(i=0;i<j;i++,k++){
	    if (y.at(i,j) != 0.0f) gold.data[k] = 65;
	  }
	}

  printf("Start Search!\n");
	hehe.searchRegul(target, cross, details, z, &gold);
  printf("End Search!\n");

	arma::mat crdmat = arma::mat(details.getSize()-1,2);
	Rcpp::CharacterMatrix namededges = Rcpp::CharacterMatrix(details.getSize()-1,2);
	for(int i =0;i<details.getSize()-1;i++){
	  crdmat.at(i,0) = details[i+1].k[0];
	  crdmat.at(i,1) = details[i+1].k[1];
	  namededges(i,0) = Rcpp::CharacterVector::create("LL")[0];//Rcpp::colnames(x)[details[i+1].k[0]];
	  namededges(i,1) = Rcpp::CharacterVector::create("LLfd")[0];//Rcpp::colnames(x)[details[i+1].k[1]];
	}


	Rcpp::NumericMatrix llincrmat(details.getSize(),13);

	int best =0;
	double bestv = details[0].d[1] - details[0].d[2];
	double tmp;
	for(int i =0;i<details.getSize();i++){
	for(int j=0;j<13;j++) llincrmat(i,j) = details[i].d[j];
	if ((ExOp::isValid(details[i].d[5]))&&(tmp > bestv)){best = i;bestv = details[i].d[5];}
	}
	Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("LL","LL2","CrX_deter","CrX_trace","CrX_LLgradient","CrX_LLmean","CrX_LLstdev","Error","CrX_error","TPrate","FPrate", "InsIntoSize", "TimeMilli");



	Rcpp::NumericMatrix fout;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());
	hehe.wrMatrix(fout);
	Rcpp::NumericMatrix dacov;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());
	target.wrMatrix(dacov);
	Rcpp::colnames(dacov) = Rcpp::colnames(x);
	Rcpp::rownames(dacov) = Rcpp::colnames(x);
	Rcpp::colnames(fout) = Rcpp::colnames(x);
	Rcpp::rownames(fout) = Rcpp::colnames(x);

	return Rcpp::List::create(Rcpp::Named("PrecisionMatrix") = fout,
	Rcpp::Named("covar") = dacov,
	Rcpp::Named("EdgeListCoor") = crdmat,
	Rcpp::Named("EdgeListName") = namededges,
	Rcpp::Named("LL_incrs") = llincrmat,
	Rcpp::Named("mean") = mu_out,
	Rcpp::Named("stddev") = std_out,
	Rcpp::Named("nbedges") = best );
}

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List InferN0_printstringlist(SEXP x) {
	std::vector<std::string> hoho = Rcpp::as<std::vector<std::string> >(x);
	int i;
	for(i=0;i<hoho.size();i++) printf("%i:%s\n",i, hoho[i].c_str());
return(Rcpp::List::create(Rcpp::Named("output") = Rcpp::wrap((double)0.0f)));}


// library(InferN0) ; cov <- readRDS("/lustre/scratch117/cellgen/team218/vk8/scRNAseqData/Pluripotency/Inf_raw_covar.rds"); out <- InferN0IdentifyNetworkOlderVersion(c(),method=cov, check.nb=20)

// report on precision and false positive rate instead of error in LL_incr
//
// [[Rcpp::export]]
Rcpp::List InferN0_identifyL0NetworkGoldResurrected(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean){
  SparseMatrix<double> inputSP; inputSP.rdRcppdgCMatrix(input);
  Rcpp::List davec = input.slot("Dimnames");
  Rcpp::CharacterVector genenames = davec[0];
  Rcpp::CharacterVector cellnames = davec[1];
  printf("Run on %i x %i\n", (int)genenames.length(), (int) cellnames.length())  ; fflush(stdout);
  Rcpp::NumericMatrix covar = Rcpp::as< Rcpp::NumericMatrix >(list[0]);
  Rcpp::NumericVector mean = Rcpp::as< Rcpp::NumericVector >(list[1]);
  Rcpp::Vector<INTSXP> cell_partition = Rcpp::as< Rcpp::Vector<INTSXP> >(list[2]);
  uint32_t maxcyclic = Rcpp::as< uint32_t >(list[3]);
  uint32_t nbthreads = Rcpp::as< uint32_t >(list[4]);
  uint32_t nbcheck = Rcpp::as< uint32_t >(list[5]);
  uint32_t mincheck = Rcpp::as< uint32_t >(list[6]);
  Rcpp::NumericMatrix rcpptarget = Rcpp::as< Rcpp::NumericMatrix >(list[7]);
  uint32_t flags = Rcpp::as< uint32_t >(list[8]);
  Trianglix<char> target; target.rdMatrix(rcpptarget);
  uint32_t nbedges = nbcheck;
  Tuple<Trianglix<double> > training; training.setSize(traincov.length());
  uint32_t i,j,k;
  Trianglix<double> totalcovar; totalcovar.rdMatrix(covar);
  for(i =0; i< training.getSize();i++) training[i].rdMatrix(Rcpp::as< Rcpp::NumericMatrix >(traincov[i]));
  Trianglix<double> wholetestset; wholetestset.setSize(totalcovar.getSize()).toZero();
  Tuple<Trianglix<double> > testset; testset.setSize(traincov.length());
  Tuple< Rcpp::NumericVector  > dameans; dameans.setSize(traincov.length());
  for(i =0; i< training.getSize();i++) {
    testset[i].setSize(totalcovar.getSize()).toZero();
    dameans[i] = Rcpp::as<Rcpp::NumericVector>(trainmean[i]);
  }
  double tmp, tmp2;
  Tuple<uint32_t> partcount;partcount.setSize(training.getSize()).toZero();
  for(i =0; i< inputSP.getNBcols();i++){
    partcount[cell_partition[i]]++;
    if (auto ite = inputSP.data[i]()) do{
      k = (ite() * (ite() +1)) >> 1;
	    tmp = (*ite) - dameans[cell_partition[i]][ite()];
      tmp2 = (*ite) - mean[ite()];
      if (auto ite2 = inputSP.data[i]()) do{
        if (ite2()<=ite()){
          testset[cell_partition[i]].data[k+ite2()] += tmp * ((*ite2) - dameans[cell_partition[i]][ite2()]);
          wholetestset.data[k+ite2()]  += tmp2 * ((*ite2) - mean[ite2()]);
        }
      }while(ite2++);
    }while(ite++);
  }
  for(i =0; i< training.getSize();i++) testset[i] /= partcount[i];
  Vector< KeyElem<Tuple<uint32_t, 2u>, Tuple<double, 11u> > > details;
  SparseTrianglix<double> precision;
  printf("Run\n")  ; fflush(stdout);


  Vector< KeyElem<Tuple<uint32_t, 2u>, Tuple<double, 9u> > > details9;
  Tuple< Trianglix<double> > cross; cross.setSize(testset.getSize()*2);
  Tuple<double> eigen;
  for(i=0;i<testset.getSize();i++){
    cross[(i<<1) | 0].toMemmove(training[i]);
    cross[(i<<1) | 1].toMemmove(testset[i]);
    eigen = cross[(i<<1) | 0].getEigenValues();
    for(j=0;j<eigen.getSize();j++) if (eigen[j] <= 0.0) break;
    if (j != eigen.getSize()) {printf("Training set%i is not positive definite! Eigenvalues:\n",i);eigen.show();}
  }
  eigen = totalcovar.getEigenValues();
  for(j=0;j<eigen.getSize();j++) if (eigen[j] <= 0.0) break;
  if (j != eigen.getSize()) {printf("Input Sample Covariance is not positive definite! Eigenvalues:\n");eigen.show();}

  printf("Start Search! for %i edges (use all %c)\n", nbcheck, (flags & 4) ? 'T' : 'F'); fflush(stdout);

//	if (flags & 4) precision.searchRegul(totalcovar, cross, details, nbcheck, &target, 10); // precision.searchRegul2017(totalcovar, cross, details9, nbcheck, &target);
// else
  precision.searchRegulHalf(totalcovar, cross, details, nbcheck, target.getSize() == totalcovar.getSize() ? &target : NULL, 10);

  printf("End Search!\n");fflush(stdout);

int best =0;
  double bestv;
  Rcpp::NumericMatrix llincrmat;

	arma::mat crdmat;
	Rcpp::CharacterMatrix namededges;

  /*if (flags & 4) {
	  crdmat = arma::mat(details9.getSize()-1,2);
	  namededges = Rcpp::CharacterMatrix(details9.getSize()-1,2);
	  for(int i =0;i<details9.getSize()-1;i++){
	    crdmat.at(i,0) = details9[i+1].k[0];
	    crdmat.at(i,1) = details9[i+1].k[1];
	    namededges(i,0) = Rcpp::CharacterVector::create("LL")[0];//Rcpp::colnames(x)[details[i+1].k[0]];
	    namededges(i,1) = Rcpp::CharacterVector::create("LLfd")[0];//Rcpp::colnames(x)[details[i+1].k[1]];
	  }

	  llincrmat = Rcpp::NumericMatrix(details9.getSize(),9);
	  bestv = details9[0].d[1] - details9[0].d[2];
	  for(int i =0;i<details9.getSize();i++){
	    for(int j=0;j<9;j++) llincrmat(i,j) = details9[i].d[j];
      tmp = details9[i].d[1] - details9[i].d[2];
	    if ((ExOp::isValid(tmp))&&(tmp > bestv)){best = i;bestv = details9[i].d[5];}
	  }
	  Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("LL","CrX_deter","CrX_trace","CrX_LLgradient","CrX_LLstdev","Error","CrX_error","TPrate","FPrate");
  }else{*/
	  crdmat = arma::mat(details.getSize()-1,2);
	  namededges = Rcpp::CharacterMatrix(details.getSize()-1,2);
	  for(int i =0;i<details.getSize()-1;i++){
	    crdmat.at(i,0) = details[i+1].k[0];
	    crdmat.at(i,1) = details[i+1].k[1];
	    namededges(i,0) = Rcpp::CharacterVector::create("LL")[0];//Rcpp::colnames(x)[details[i+1].k[0]];
	    namededges(i,1) = Rcpp::CharacterVector::create("LLfd")[0];//Rcpp::colnames(x)[details[i+1].k[1]];
	  }

	  llincrmat = Rcpp::NumericMatrix(details.getSize(),11);
	  bestv = details[0].d[1] - details[0].d[2];
	  for(int i =0;i<details.getSize();i++){
	    for(int j=0;j<11;j++) llincrmat(i,j) = details[i].d[j];
      tmp = details[i].d[1] - details[i].d[2];
	    if ((ExOp::isValid(tmp))&&(tmp > bestv)){best = i;bestv = tmp;}
	  }
    if (target.getSize() == totalcovar.getSize()) Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("LL","CrX_deter","CrX_trace","CrX_LLgradient","CrX_LLstdev","Error","CrX_error","TPrate","FPrate", "InsIntoSize", "TimeMilli");
	  else Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("LL","CrX_deter","CrX_trace","CrX_LLgradient","CrX_LLstdev","Error","CrX_error","Empty","Empty2", "InsIntoSize", "TimeMilli");
  //}
	Rcpp::NumericMatrix fout;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());
	precision.wrMatrix(fout);
	Rcpp::colnames(fout) = genenames;
	Rcpp::rownames(fout) = genenames;
return Rcpp::List::create(Rcpp::Named("PrecisionMatrix") = fout,Rcpp::Named("EdgeListCoor") = crdmat,Rcpp::Named("EdgeListName") = namededges,Rcpp::Named("LL_incrs") = llincrmat,Rcpp::Named("nbedges") = best );}

Rcpp::NumericMatrix wilcox(const SparseMatrix<uint32_t> &inputSP,  uint32_t flag, Rcpp::List list, Rcpp::CharacterVector genenames){
	vector<uint32_t> alist = Rcpp::as<vector<uint32_t> > (list[0]);
	vector<uint32_t> blist = Rcpp::as<vector<uint32_t> > (list[1]);
	Tuple<uint32_t> groups[2];
	groups[0].setSize(alist.size());
	for(uint32_t i=0;i<alist.size(); i++) groups[0][i] = alist[i]-1;
	groups[1].setSize(blist.size());
	for(uint32_t i=0;i<blist.size();i++) groups[1][i] = blist[i]-1;
	Tuple<double> logitauroc;
	Tuple<double> seenlogitauroc;
	Tuple<double> output = inputSP.UTestZScore(groups[0], groups[1], &logitauroc, &seenlogitauroc, (flag & 8) == 8,false,(flag & 16) == 16);
	Rcpp::NumericMatrix fout(output.getSize(),(flag & 2) ? 6 : 3);
	Rcpp::CharacterVector cvec;
	cvec.push_back("Zscore");
	cvec.push_back("Logit_AUROC");
	cvec.push_back("Logit_AUROC_in_obs");
	if (flag & 2){
		WeightElem<double, 2u> average;
		uint32_t j,k;
		for(uint32_t i=0;i< groups[0].getSize();i++){
			average.toZero();
			for(uint32_t j=0;j<inputSP.getNBcols();j++){
				if ((k = inputSP.data[j].find(groups[0][j])) != 0xFFFFFFFF) average += WeightElem<double, 2u>( inputSP.data[j].deref(k));
			}
			fout(i,3) = average.getMean();
		}
		for(uint32_t i=0;i< groups[1].getSize();i++){
			average.toZero();
			for(uint32_t j=0;j<inputSP.getNBcols();j++){
				if ((k = inputSP.data[j].find(groups[1][j])) != 0xFFFFFFFF) average += WeightElem<double, 2u>( inputSP.data[j].deref(k));
			}
			fout(i,4) = average.getMean();
			fout(i,5) = (log(fout(i,3)) - log(fout(i,4)))/log(2);
		}
		cvec.push_back("mean_in_obs_pos");
		cvec.push_back("mean_in_obs_neg");
		cvec.push_back("Fold_change_in_obs");
	}
	for(uint32_t i=0;i< output.getSize();i++){
		fout(i,0) = output[i];
		fout(i,1) = logitauroc[i];
		fout(i,2) = seenlogitauroc[i];
	}
	Rcpp::colnames(fout) = cvec;
	Rcpp::rownames(fout) =  genenames;

return fout;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_wilcox_geneset)]]
Rcpp::List infernal_wilcox_geneset(SEXP scope, Rcpp::List list){
	Rcpp::XPtr< Infern0scope > scp(scope);
	vector<uint32_t> alist = Rcpp::as<vector<uint32_t> > (list[0]);
	vector<uint32_t> blist = Rcpp::as<vector<uint32_t> > (list[1]);
	vector<uint32_t> olist = Rcpp::as<vector<uint32_t> > (list[2]);
	vector<uint32_t> plist = Rcpp::as<vector<uint32_t> > (list[3]);
	uint32_t nbquartiles = Rcpp::as<uint32_t > (list[4]);
return Rcpp::List::create();}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_wilcox_scope)]]
Rcpp::NumericMatrix infernal_wilcox_scope(SEXP scope, Rcpp::List list){
	Rcpp::XPtr< Infern0scope > scp(scope);
	uint32_t flag = Rcpp::as<uint32_t > (list[2]);
	Rcpp::CharacterVector genenames = Rcpp::CharacterVector();
	scp->annot.getAxeNamesAl(genenames, (int)INFRN0_AXES_CELL_NAMES);
return wilcox(scp->rawdata,flag,list,genenames);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_wilcox_matrix)]]
Rcpp::NumericMatrix infernal_wilcox_matrix(Rcpp::S4 input, Rcpp::List list){
	uint32_t flag = Rcpp::as<uint32_t > (list[2]);
	SparseMatrix<uint32_t> inputSP; inputSP.rdRcppdgCMatrix(input, (flag & 1) == 1);
	Rcpp::List davec = input.slot("Dimnames");
        Rcpp::CharacterVector genenames = davec[0];
	if (inputSP.getNBcols() < genenames.length()) inputSP.data.upSize(genenames.length());
return wilcox(inputSP,flag,list,genenames);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_find_DE)]]
Rcpp::List infernal_find_DE(SEXP scope, Rcpp::List list){
	Rcpp::XPtr< Infern0scope > scp(scope);
	Rcpp::CharacterVector genenames = Rcpp::CharacterVector();
	scp->annot.getAxeNamesAl(genenames, (int)INFRN0_AXES_CELL_NAMES);
	Tuple<uint32_t> ordering, partit;
	Tuple<Tuple<uint32_t> > groups; groups.setSize(2);
	vector<uint32_t> alist = Rcpp::as<vector<uint32_t> > (list[0]);
	vector<uint32_t> blist = Rcpp::as<vector<uint32_t> > (list[1]);
	vector<uint32_t> olist = Rcpp::as<vector<uint32_t> > (list[2]);
	vector<uint32_t> plist = Rcpp::as<vector<uint32_t> > (list[3]);
	uint32_t nbquartiles = Rcpp::as<uint32_t > (list[4]);
	uint32_t nbthreads = Rcpp::as<uint32_t > (list[5]);
	uint32_t flag = Rcpp::as<uint32_t > (list[6]);
	uint32_t i,j;
	groups[0].setSize(alist.size());
	for(uint32_t i=0;i<alist.size(); i++) groups[0][i] = alist[i]-1;
	groups[1].setSize(blist.size());
	for(uint32_t i=0;i<blist.size();i++) groups[1][i] = blist[i]-1;
	TMatrix<double> output2, output3, output4, output5;

	if ((groups[0].getSize() == 0)||(groups[1].getSize() == 0)) return(Rcpp::List::create());
//	printf("start!; %i vs %i (in %i)\n", groups[0].getSize(),groups[1].getSize(), scp->rawdata.getNBcols());

	SparseMatrix<uint32_t> subdata;
	// if (glist.size() != 0) subdata = scp->rawdata.subseddtRows(myglist);
	const SparseMatrix<uint32_t> &inputref = scp->rawdata;

	ordering.setSize(olist.size());
	for(i=0;i<olist.size();i++) ordering[i] = olist[i];

	TMatrix<Zscore> output = inputref.UTestZScoreSplit(groups, ordering, nbthreads, nbquartiles, partit,  &output2, &output3, &output4, flag & 1, flag & 2, flag & 64);
	Rcpp::NumericMatrix zscores;
	Rcpp::NumericMatrix logitauroc;
	Rcpp::NumericMatrix mean;
	Rcpp::NumericMatrix drop_enrich;
	Rcpp::NumericMatrix zweight;

	Rcpp::CharacterVector cvec = Rcpp::CharacterVector();
    	scp->annot.getAxeNamesAl(cvec, (int)INFRN0_AXES_CELL_NAMES);

	output2.wrMatrix(logitauroc);
	output3.wrMatrix(mean);
	output4.wrMatrix(drop_enrich);
	output5.setSizes(output2.nbRows(), output2.nbCols());

	if (auto ite = output.getIterator()) do{
		output5(ite()) = ite->weight;
	}while(ite++);
        output5.wrMatrix(zweight);

	if (auto ite = output.getIterator()) do{
		output5(ite()) = (double)*ite;
	}while(ite++);
        output5.wrMatrix(zscores);

	Rcpp::rownames(zscores) = cvec;
	Rcpp::rownames(logitauroc) = cvec;
	Rcpp::rownames(mean) = cvec;
	Rcpp::rownames(drop_enrich) = cvec;
	Rcpp::rownames(zweight) = cvec;

	Rcpp::CharacterVector cvec2 = Rcpp::CharacterVector(nbquartiles *2);
	char buffer[256];
	for(i=0;i<nbquartiles;i++){sprintf(buffer,"Pos_q%i", i+1);cvec2[i << 1] = buffer;	sprintf(buffer,"Neg_q%i", i+1);cvec2[(i << 1)|1] = buffer;}
	Rcpp::colnames(zscores) = cvec2;
	Rcpp::colnames(logitauroc) = cvec2;
	Rcpp::colnames(mean) = cvec2;
	Rcpp::colnames(drop_enrich) = cvec2;
	Rcpp::colnames(zweight) = cvec2;

	double daval;


	if (flag & 128){
		// overwite data with DE value, assumes the number of partitions is 1
		if (scp->data.getNBcols() != scp->rawdata.getNBcols()) scp->data.data.setSize(scp->rawdata.getNBcols());
		/*for(j=0;j<2;j++){
		for(i = 0; i< groups[j].getSize();i++){
			if (auto ite = scp->rawdata.data[groups[j][i]]()) do{
				scp->data.data[groups[j][i]][ite()] = 1.0 - 2.0 / (1.0 - exp(output2(ite(),1)));
			}while(ite++);
		}i
		}*/

		myHashmap<uint32_t> dalistlist;
		for(j=0;j<2;j++) for(i = 0; i< groups[j].getSize();i++) dalistlist.addEntry(groups[j][i]);
		for(i = 0; i< scp->rawdata.getNBcols();i++){
			if (fabs(output5(i,0)) > 1.95) {
				daval = 1.0 - 2.0 / (1.0 + exp(output2(i,0)));
				if (auto ite = scp->rawdata.data[i]()) do{
					if (dalistlist.find(ite()) != 0xFFFFFFFF) scp->data.data[i][ite()] = daval;
				}while(ite++);
			}else{
				//if (auto ite = scp->rawdata.data[i]()) do{
				//	if (dalistlist.find(ite()) != 0xFFFFFFFF) scp->data.data[i].erase(ite());
				//}while(ite++);
			}
		}
	}
return Rcpp::List::create(Rcpp::Named("Zscore") = zscores,Rcpp::Named("LogitAuroc") = logitauroc,Rcpp::Named("Mean") = mean, Rcpp::Named("CoverageEnrichment") = drop_enrich, Rcpp::Named("Weight") = zweight);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_find_markers)]]
Rcpp::List infernal_find_markers(SEXP scope, Rcpp::List rlist){
	Rcpp::XPtr< Infern0scope > scp(scope);
	uint32_t flag = Rcpp::as<uint32_t > (rlist[3]);
	Rcpp::CharacterVector genenames = Rcpp::CharacterVector();
	scp->annot.getAxeNamesAl(genenames, (int)INFRN0_AXES_CELL_NAMES);
	Tuple<uint32_t> ordering,partit;
	Tuple<Tuple<uint32_t> > groups;
	Vector<Tuple<uint32_t> > pregroups;
	vector<uint32_t> partition = Rcpp::as<vector<uint32_t> > (rlist[0]);
	vector<uint32_t> olist = Rcpp::as<vector<uint32_t> > (rlist[1]);
	uint32_t nbquartiles = Rcpp::as<uint32_t > (rlist[2]);
	vector<bool> celluse = Rcpp::as<vector<bool> > (rlist[4]);

	for(uint32_t i=0;i<partition.size();i++){
		if (!celluse[i]) continue;
		if (partition[i] >= pregroups.getSize()) {
			pregroups.upSize(partition[i]+1);
		}
		pregroups[partition[i]].push_back(i);
	}
	groups.toMemmove(pregroups);
	ordering.setSize(partition.size());
	if (nbquartiles != 1){
		if (olist.size() != ordering.getSize()){
			ordering = scp->gene_ordering;
			printf("used gene ordering! %i %i\n", scp->gene_ordering.getSize(), (int) partition.size());
		}else{
			for(uint32_t i=0;i<olist.size();i++) ordering[i] = olist[i]-1;
		}
	}

	for(uint32_t i=0;i<olist.size();i++) ordering[i] = olist[i]-1;
	TMatrix<double> output2, output3, output4;

	TMatrix<Zscore> output = scp->rawdata.UTestZScoreSplit(groups, ordering, 1, nbquartiles, partit, &output2, &output3, &output4 , flag & 1, (flag & 2) == 0, (flag & 64));

	Rcpp::NumericMatrix tfidf;
        scp->rawdata.computeTFIDF(groups).wrMatrix(tfidf);


	printf("output  has sizes %i, %i\n", output.sizes[0], output.sizes[1]);
	printf("output2 has sizes %i, %i\n", output2.sizes[0], output2.sizes[1]);
	printf("output3 has sizes %i, %i\n", output3.sizes[0], output3.sizes[1]);
	printf("output4 has sizes %i, %i\n", output4.sizes[0], output4.sizes[1]);

	TMatrix<double> possize, negsize;
	possize.setSizes(output.sizes[0],groups.getSize()).toZero();
	negsize.setSizes(output.sizes[0],groups.getSize()).toZero();
	Tuple<uint32_t> totals; totals.setSize(groups.getSize()).toZero();
	Tuple<uint32_t,2> coor;

	for(uint32_t i=0;i<partition.size();i++){
		if (celluse[i]) totals[partition[i]]++;
	}
	for(uint32_t i=0;i<scp->rawdata.data.getSize();i++){
		coor[0] = i;
		if (auto ite = scp->rawdata.data[i]()) do{
			coor[1] = ite();
			if (!celluse[ite()]) continue;
			coor[1] = partition[ite()];
			possize(coor) += 1.0;
		}while(ite++);
	}
	uint32_t datotal=0;
	uint32_t stotals;
	for(int i =0; i < totals.getSize(); i++) datotal += totals[i];

	ExOp::show(totals);
	for(coor[0]=0; coor[0] < output.sizes[0];coor[0]++){
		coor[1] = 0;
		stotals = 0;
		for(coor[1] =0; coor[1] < totals.getSize(); coor[1]++) stotals += possize(coor);
		for(coor[1] =0; coor[1] < totals.getSize(); coor[1]++) {
			negsize(coor) = (stotals - possize(coor)) / (datotal - totals[coor[1]]);
	//		possize(coor) /= totals[coor[1]];
		}
	}
	ExOp::show(totals);

	Rcpp::NumericMatrix zscores;
	Rcpp::NumericMatrix logitauroc,mean,covenrich,zweight,poscov,negcov;

	Rcpp::CharacterVector cvec = Rcpp::CharacterVector();
    	scp->annot.getAxeNamesAl(cvec, (int)INFRN0_AXES_CELL_NAMES);

	output2.wrMatrix(logitauroc);
	output3.wrMatrix(mean);
	output4.wrMatrix(covenrich);
	possize.wrMatrix(poscov);
	negsize.wrMatrix(negcov);


	if (auto ite = output.getIterator()) do{
		output2(ite()) = (double)*ite;
	}while(ite++);
        output2.wrMatrix(zscores);
	if (auto ite = output.getIterator()) do{
		output2(ite()) = ite->weight;
	}while(ite++);
        output2.wrMatrix(zweight);

	Rcpp::rownames(zscores) = cvec;
	Rcpp::rownames(logitauroc) = cvec;
	Rcpp::rownames(mean) = cvec;
	Rcpp::rownames(covenrich) = cvec;
	Rcpp::rownames(zweight) = cvec;
	Rcpp::rownames(poscov) = cvec;
	Rcpp::rownames(negcov) = cvec;
	Rcpp::rownames(tfidf) = cvec;


	Rcpp::CharacterVector cvec2 = Rcpp::CharacterVector(nbquartiles *groups.getSize());
	char buffer[256];
	for(unsigned int i=0;i<nbquartiles;i++){
		for(int j=0;j< groups.getSize();j++){
		sprintf(buffer,"Cl%i_q%i", j, i+1); cvec2[i * groups.getSize() +j] = buffer;
		}
	}
	Rcpp::colnames(zscores) = cvec2;
	Rcpp::colnames(logitauroc) = cvec2;
	Rcpp::colnames(mean) = cvec2;
	Rcpp::colnames(covenrich) = cvec2;
	Rcpp::colnames(zweight) = cvec2;

	Rcpp::CharacterVector cvec3 = Rcpp::CharacterVector(groups.getSize());
	for(int j=0;j< groups.getSize();j++) {sprintf(buffer,"Cl%i", j); cvec3[j] = buffer;}
	Rcpp::colnames(tfidf) = cvec3;

	strcpy(buffer, "Zscore"); //, (flag & 1) ? "LogitPval": "Zscore");

return Rcpp::List::create(Rcpp::Named(buffer) = zscores,Rcpp::Named("LogitAuroc") = logitauroc,	Rcpp::Named("Mean") = mean, Rcpp::Named("CoverageEnrichment") = covenrich, Rcpp::Named("PosCoverage") = poscov, Rcpp::Named("NegCoverage") = negcov, Rcpp::Named("Weight") = zweight, Rcpp::Named("TfIdf") = tfidf);}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_sharedDE)]]
Rcpp::List infernal_sharedDE(Rcpp::S4 input){
  SparseMatrix<double> dainp;
  dainp.rdRcppdgCMatrix(input);

  TMatrix<double> output;

  Rcpp::NumericMatrix zscores;
  output.wrMatrix(zscores);

  return Rcpp::List::create(Rcpp::Named("Zscores") = zscores);
}

// d <- 14; tmp <- svd(matrix(runif(d*d),d,d)); tmp <- tmp$u %*% diag((runif(d) - 0.5) * exp(10 * (runif(d)-0.5))) %*% t(tmp$u);
// d <- 14; tmp <- matrix(runif(d*d),d,d)
// out <- InferN0:::.infernal_cmpInverses(tmp); sum(abs(out$Arma %*% out$Input -diag(d))) / (d*d); sum(abs(out$LFH %*% out$Input - diag(d))) / (d*d)
//  (out$Diff / out$Arma)/(ncol(out$Diff)-1)

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_cmpInverses)]]
Rcpp::List infernal_cmpInverses(Rcpp::NumericMatrix input){
  Trianglix<double> target;
  target.rdMatrix(input);
  uint32_t milli, milli2, milli3;
  milli = clock();
  Trianglix<double> res1 = target.mkInverse();
  milli2 = clock();
  Magscaled<double> opt;
  short ispos;
  Trianglix<double> res2 = target.mkInverse(&opt, &ispos);
  milli3 = clock();

  double det1 =  target.determinant();
  double det2 = target.log_determinant();

  printf("is definite matrix... code %i!!!\n", ispos);
  printf("Arma vs LFH  %e and %e secs\n",((double)milli2- milli) / CLOCKS_PER_SEC, ((double)milli3-milli2) / CLOCKS_PER_SEC);
  printf("%e and %e !!!\n", det1, det2);
  Rcpp::NumericMatrix fout, fout2, fout3, fout4, fout5;
  target.wrMatrix(fout4);
  res1.wrMatrix(fout);
  res2.wrMatrix(fout2);
  res2 -= res1;
  res2.wrMatrix(fout3);
  printf("%e, log det %e\n", opt.value * exp(opt.exponent), log(fabs(opt.value))  + opt.exponent);
return Rcpp::List::create(Rcpp::Named("Input") = fout4,Rcpp::Named("Arma") = fout,Rcpp::Named("LFH") = fout2, Rcpp::Named("Diff") = fout3, Rcpp::Named("det") = log(abs(opt.value)) + opt.exponent);}


// Find network in a small gene subset
//
//' param list of genes to organize into a network (<30 genes)
//' param list of cells that are used, (empty list means all cells)
//
//' export
/*
Rcpp::List InferN0_findNetwork(SEXP genes, SEXP cells) {
    PartialGaussElem<double> target[9];
    vector<std::string> genames = Rcpp::as<std::vector<std::string> >(genes);

    printf("nbdico %i\n data sizes: %i,%i\n", stscope.annot.dicos.getSize(), stscope.data.row_sizes.getSize(), stscope.data.data.getSize());


    uint32_t ite = stscope.annot.findDico(stscope.data,1);
    if (ite == 0xFFFFFFFF){
        printf("Error, data is not loaded\n");
        return Rcpp::List::create();
    }
    Dictionnary& genedico = stscope.annot.dicos[ite];

    Vector<uint32_t> goffsets;
    for(uint32_t i ;i < genames.size();i++){
        ite = genedico.find(genames[i].c_str());
  //      if (ite != 0xFFFFFFFF) goffsets.push_back(ite);
  //      else printf("Warning, did not find gene %s\n", genames[i].c_str());
    }
    printf("Found %i/%i\n", goffsets.getSize(), genames.size());
    uint32_t i,j;
    for(i=0;i<9;i++) target[i].setSize(goffsets.getSize());

    for(i=0;i<stscope.data.data.getSize();i++){
        target[i&7].toAddFiltCol(stscope.data.data[i], &(goffsets[0]));
    }
    target[8] = target[0];
    for(i++;i<8;i++) target[8] += target[i];



return(Rcpp::List::create(Rcpp::Named("output") = Rcpp::wrap((double)0.0f)));}*/

//' Save InferN0 Scope
//'
//' param path file path for the scope to save
//'
//' export
/*
void saveScope(SEXP path){
    try{
        std::string spath = Rcpp::as<std::string>(path);
        FILE *f = fopen(spath.c_str(), "w+");
        if (f == NULL) {printf("Could not open %s for wrinting.\n", spath.c_str()); return;}
        if (stscope.save(f)) {printf("Error found while writing scope.\n", spath.c_str());}
        fclose(f);
    }catch(std::exception &ex){
       printf("Expected \"path\" argument would be a string\n");
    }
}*/


//' Load InferN0 Scope
//'
//' param path file path to a saved scope
//'
//' export
/*
void loadScope(SEXP path){
    try{
        std::string spath = Rcpp::as<std::string>(path);
        FILE *f = fopen(spath.c_str(), "r+");
        if (f == NULL) {printf("Could not open %s for reading.\n", spath.c_str()); return;}
        if (stscope.load(f)) {printf("Error found while reading scope.\n", spath.c_str());}
        fclose(f);
    }catch(std::exception &ex){
       printf("Expected \"path\" argument would be a string");
    }
}*/

//' Check S4 Class Object
//'
//' @param S4 Class object
//' @param string name of queried slot
//'
//' @export
// [[Rcpp::export]]
void testSlot(Rcpp::S4 obj, SEXP slot){
    std::string spath = Rcpp::as<std::string>(slot);
    printf("Object has slot %s? %c\n", spath.c_str(), obj.hasSlot(spath.c_str()) ? 'Y' : 'N');
}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_serialize)]]
Rcpp::RawVector infernal_serialize(Rcpp::S4 obj){
	Rcpp::XPtr< Infern0scope > scp(obj);
	FILE* f = tmpfile();
	scp->save(f);
	uint64_t length = ftell(f);
	Rcpp::RawVector fout(length);
	fseek(f,SEEK_SET, 0u);
	if (fread(&fout[0], 1 , length, f) != length) printf("did not read serialized data fully!\n");
return fout;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_deserialize)]]
void infernal_deserialize(Rcpp::S4 obj, Rcpp::RawVector src){
	Rcpp::XPtr< Infern0scope > scp(obj);
	FILE* f = tmpfile();
	//uint64_t cur =0;
	//while( src.size() -  cur > 65535){
	//	fwrite(&src[cur], 1 , 65536, f);
//		cur += 65536;
//	}
//	fwrite(&src[cur], 1 , src.size() -  cur 65536, f);
	fwrite(&src[0], 1 , src.size(), f);
	fseek(f,SEEK_SET, 0u);
	// TODO convert to vector stream
	scp->load(f);
return;}





/*

loopy <- 2
param <- c(1000,5,20,0.1,0.5,0)

graph <- GetGoldStandardNetwork()
for(i in 1:loopy){
  data <-GenerateSyntheticData2(GetGoldStandardNetwork(),param[1],mu.min=param[2],mu.max=param[3],p.min=param[4],p.max=param[5], s=param[6])
  srcout <- InferN0_identifyL0NetworkGold(data,GetGoldStandardNetwork(),50)
  if (i == 1){
    summy <- srcout$LL_incrs[,c("TPrate","FPrate")]
    best <- srcout$LL_incrs[(srcout$nbedges+1),c("TPrate","FPrate")]
  }else{
    summy <- summy + srcout$LL_incrs[,c("TPrate","FPrate")]
    best <- best + srcout$LL_incrs[(srcout$nbedges+1),c("TPrate","FPrate")]
  }
}
summy <- summy / loopy
best <- best/ loopy
saveRDS(list(best=best,summy=summy,nb=loopy,graph=graph,param=param,paramnames=c("nbcell","mu.min","mu.max","p.min","p.max","noise.mu")),"Infernoresult.rds")
quit("no")
*/



