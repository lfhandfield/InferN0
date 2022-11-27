// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// Internal functions, no documentation provided, use with care...

#define ARMA_DONT_PRINT_ERRORS

#include "RcppArmadillo.h"
#include "infern0.hpp"

using namespace LFHPrimitive;
using namespace arma;


// list of internal functions that are called by the R wrapper, see R wrapper code instead for documentation

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_createScope)]]
SEXP infernal_createScope(){Rcpp::XPtr< Infern0scope > p(new Infern0scope(),true); return wrap(p);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_readScope)]]
SEXP infernal_readScope(SEXP path){Rcpp::XPtr< Infern0scope > p(new Infern0scope(),true);
  std::string pathvar = Rcpp::as<std::string>(path);
  printf("reading from '%s'\n",pathvar.c_str());
  FILE* f = fopen(pathvar.c_str(), "rb");
  if (f == NULL) {
	printf("could not open %s\n", pathvar.c_str());
        return Rcpp::wrap(NULL);
  }else{
    p->load(f);
    fclose(f);
  }
return wrap(p);}
//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_readScope2)]]
SEXP infernal_readScope2(SEXP path){Rcpp::XPtr< Infern0scope > p(new Infern0scope(),true);
  std::string pathvar = Rcpp::as<std::string>(path);
  printf("reading from '%s'\n",pathvar.c_str());
  FILE* f = fopen(pathvar.c_str(), "rb");
  if (f == NULL) {
	printf("could not open %s\n", pathvar.c_str());
        return Rcpp::wrap(NULL);
  }
    uint32_t length;
    if (fread(&length, sizeof(int32_t), 1, f) != 1) return Rcpp::List::create();
    Rcpp::RawVector serialized(length);
    if (fread((char*)&serialized[0], 1 , length, f) != length) return Rcpp::List::create();
    p->load(f);
    fclose(f);
return Rcpp::List::create(Rcpp::Named("ptr")=wrap(p), Rcpp::Named("serial") = wrap(serialized));}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_writeScope)]]
void infernal_writeScope(SEXP infscope, SEXP path){
  Rcpp::XPtr< Infern0scope > scptr(infscope);
  std::string pathvar = Rcpp::as<std::string>(path);
  FILE* f = fopen(pathvar.c_str(), "wb+"); if (f == NULL) {printf("could not save to '%s'\n",pathvar.c_str()); return;}
  printf("saving to '%s'\n",pathvar.c_str());
  scptr->save(f);
  fclose(f);
return;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_writeScope2)]]
void infernal_writeScope2(SEXP infscope, Rcpp::RawVector serialized, SEXP path){
  Rcpp::XPtr< Infern0scope > scptr(infscope);
  std::string pathvar = Rcpp::as<std::string>(path);
  FILE* f = fopen(pathvar.c_str(), "wb+"); if (f == NULL) {printf("could not save to '%s'\n",pathvar.c_str()); return;}
  uint32_t length = serialized.size();
  if (fwrite(&length, sizeof(int32_t), 1, f) != 1) {  printf("could not save to '%s'\n",pathvar.c_str()); return;}
  printf("saving to '%s'\n",pathvar.c_str());
  fwrite((char*)&serialized[0], 1 , length, f);
  scptr->save(f);
  fclose(f);
return;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_readMatrixFolder)]]
Rcpp::List infernal_readMatrixFolder(Rcpp::List list){
  char buffer[512];
  char sep;
  std::string path = Rcpp::as< std::string >(list[0]);
  uint32_t umithresthold = Rcpp::as< uint32_t >(list[1]);
  uint32_t genethresthold = Rcpp::as< uint32_t >(list[2]);
  std::vector<std::string> ignore = Rcpp::as<std::vector<std::string> >(list[3]);
  double soupfraction = Rcpp::as< double >(list[4]);
  double souplog10slope = Rcpp::as< double >(list[5]);
  std::vector<uint32_t> soup_samples = Rcpp::as<std::vector<uint32_t> >(list[6]);
  uint32_t do_progress = Rcpp::as< uint32_t >(list[7]);


  uint32_t l = path.length();

  uint32_t buf[3];
  double dbuf;
  uint32_t nbgenes;
  uint32_t nbentries;

  myHashmap<uint32_t,uint32_t> allfiltered, allunfiltered;
  myHashmap<uint32_t,double> meatfraction;

  SparseMatrix<uint32_t> histogram; histogram.data.setSize(5);

  SparseMatrix<uint32_t> data;
  SparseMatrix<double> ddata;
  SparseMatrix<uint32_t> softfiltered_data;
  SparseMatrix<uint32_t> filt_data[3];
  if (umithresthold > 1) filt_data[0].data.setSize(umithresthold);
  if (genethresthold > 1) filt_data[1].data.setSize(genethresthold-1);

  Rcpp::CharacterVector genenames,cellnames,filtcellnames, geneacc_names;
  uint32_t fscout;

  myHashmap<string,uint32_t> toignore;
  myHashmap<uint32_t,void> toignore_index;
  for(fscout=0;fscout<ignore.size();fscout++) toignore[ignore[fscout].c_str()]  = 7;

  memcpy(buffer, path.c_str(), l);
  strcpy(buffer+l, "matrix.mtx");
  FILE* g = fopen(buffer,"r");
  if (g == NULL) {fprintf(stderr, "could not open '%s'!\n", buffer); return Rcpp::List::create();}

  strcpy(buffer+l, "genes.tsv");
  FILE* f = fopen(buffer,"r");
  bool isfeat;
  if (isfeat = (f == NULL)) {
	fprintf(stderr, "could not open '%s'!\n", buffer);
	fprintf(stderr, "trying features.tsv instead!\n");
	strcpy(buffer+l, "features.tsv");
	f = fopen(buffer,"r");
        if (f == NULL) {fprintf(stderr, "could not open '%s'!\n", buffer); return Rcpp::List::create();}
 }

  while(1 == fscanf(g," %1[%]", buf)) {
	if (fscanf(g,"%*[^\n]\n") != 0 ) return R_NilValue; // comment line detected
  }
  if (3 != fscanf(g,"%i %i %i\n", &nbgenes,buf+1,&nbentries)) {fprintf(stderr, "Failed to read Header of file '%s' ", buffer); return Rcpp::List::create();}
  printf("%i genes, %i cells, %i entries\n", nbgenes , buf[1],nbentries); fflush(stdout);

  uint32_t nbannot=1;
  if (fscanf(f,"%*[^\t\n]%1[\t\n]", &sep) != 1) return R_NilValue;
  if (sep == '\t'){
    do{
    nbannot++;
    if (fscanf(f,"%*[^\t\n]%1[\t\n]", &sep) != 1) return R_NilValue;
    }while(sep != '\n');
  }
  fclose(f);
  f = fopen(buffer,"r");
  uint32_t UMIc =0;

  genenames = Rcpp::CharacterVector(nbgenes);

  if (nbannot > 2){
    geneacc_names = Rcpp::CharacterVector(nbgenes);
    printf("Reading Gene Accesion and Gene Names..."); fflush(stdout);
    while(2 == (fscout = fscanf(f,"%[^ \t]%*1[ \t]%[^ \t]%*1[ \t]%*[^\n]\n", buffer,buffer + 256))){
      if (UMIc < genenames.length()){genenames[UMIc] = buffer+256;geneacc_names[UMIc] = buffer;
      }else{genenames.push_back(buffer+256);geneacc_names.push_back(buffer);}
      if ((fscout = toignore.find(buffer+256)) != 0xFFFFFFFF) toignore_index.addEntry(UMIc);
      UMIc++;
    }
  }else if (nbannot == 2){ // assumes columns are geneaccession and genenames respectively
    geneacc_names = Rcpp::CharacterVector(nbgenes);
    printf("Reading Gene Accesion and Gene Names..."); fflush(stdout);
    while(2 == (fscout = fscanf(f,"%[^ \t]%*1[ \t]%[^\n]\n", buffer,buffer + 256))){
      if (UMIc < genenames.length()){genenames[UMIc] = buffer+256;geneacc_names[UMIc] = buffer;
      }else{genenames.push_back(buffer+256);geneacc_names.push_back(buffer);}
      if ((fscout = toignore.find(buffer+256)) != 0xFFFFFFFF) toignore_index.addEntry(UMIc);
      UMIc++;
    }
  }else{ // use everything
    printf("Reading Gene Names..."); fflush(stdout);
    while(1 == (fscout = fscanf(f,"%[^\n]\n", buffer))){
      if (UMIc < genenames.length()) genenames[UMIc] = buffer;
      else genenames.push_back(buffer);
      if ((fscout = toignore.find(buffer)) != 0xFFFFFFFF) toignore_index.addEntry(UMIc);
      UMIc++;
    }
  }
  fclose(f);


  printf("(Done), %i genes match filter list\n",toignore_index.getSize()); // library(InferN0); data <- InferN0ReadRangerMatrixCPP("/warehouse/team218_wh01/lh20/references/Zhong2018/", do.soup.regression=F)

  uint32_t old =1; // starts at 1 eh?
  UMIc =0;
  uint32_t UMItc =0;
  uint32_t UMItc_max =0;
  uint32_t samplethr[2];


  Vector<uint32_t> selcol;
  myHashmap<uint32_t, uint32_t> curcol;
  Rcpp::S4 fout;
	ProgressBarPrint progbar(20);
	if (do_progress & 1) progbar.start("Reading Sparse Matrix");
  uint32_t cur=0;

  if (do_progress & 4){ // floating point and not filetering!
     while(true){
      if ((do_progress & 1)&&(((cur++) & 63) == 0)) progbar.update(((double)cur) / nbentries);
      if (3 != fscanf(g,"%i %i %lf\n", buf,buf+1,&dbuf)) {
        if (do_progress & 1) progbar.finish();
        ddata.wrRcppdgCMatrix(fout, nbgenes);
        printf("Reading Cell Names..."); fflush(stdout);

        memcpy(buffer, path.c_str(), l);
        strcpy(buffer+l, "barcodes.tsv");
        f = fopen(buffer,"r+");
        if (f == NULL) {fprintf(stderr, "could not open '%s'!\n", buffer); return Rcpp::List::create();}
        for(buf[0] =0,buf[1] =0; 1 == (fscout = fscanf(g,"%[^\n]\n", buffer)); buf[0]++) cellnames.push_back(buffer);
        fclose(f);
        printf("(Done)\n");

        if (do_progress & 8) { // tries to recover raw counts: assumes unique denominator was used for all counts
          dbuf =0;
          for(old=0;old < ddata.data.getSize();old++){
            for(l=0;l <  ddata.data[old].getSize();l++) dbuf += ddata.data[old].deref(l);
          }
          printf("total is %f /(= 1000000?)\n", dbuf);
          dbuf /= 1000000.0f;
          for(old=0;old < ddata.data.getSize();old++){
            for(l=0;l <  ddata.data[old].getSize();l++) ddata.data[old].deref(l) *= dbuf;
          }
        }

        fout.slot("Dimnames") = Rcpp::List::create(genenames,cellnames);
        return Rcpp::List::create(Rcpp::Named("data") = fout);
      }
      while(buf[1] > ddata.data.getSize()) ddata.data.push_back();
      ddata.data[buf[1]-1][buf[0]-1] = dbuf; // starts at 1 eh?
    }
  }else{

    if (soup_samples.size() != 0) filt_data[2].data.setSize(soup_samples[2]);
    else{
        soup_samples.resize(2);
        soup_samples[0] = 10;
        soup_samples[1] = 0;
    }


    while(true){
      if ((do_progress & 1)&&(((cur++) & 63) == 0)) progbar.update(((double)cur) / nbentries);
      fscout = fscanf(g,"%i %i %i\n", buf,buf+1,buf+2);
      if ((fscout != 3)||(old != buf[1])){
	if ((UMIc >= umithresthold)&&(curcol.getSize() >= genethresthold)){
          if (auto ite = curcol()) do{
            allunfiltered[ite()] += *ite;
          }while(ite++);
          data.data.mempush_back(curcol);
          selcol.push_back(old-1);
          histogram.data[2][UMIc]++;
          histogram.data[3][UMItc]++;
        }else{
          if (auto ite = curcol()) do{
            allfiltered[ite()] += *ite;
          }while(ite++);
          if ((UMIc >= soup_samples[0])&&(UMIc <= soup_samples[1])){
		uint32_t sss = rand() % soup_samples[2];
		if (auto ite = curcol()) do{filt_data[2].data[sss][ite()] += *ite;}while(ite++);
	  }
          if (UMIc < umithresthold){
            if (auto ite = curcol()) do{
              filt_data[0].data[UMIc][ite()] += *ite;
            }while(ite++);
            histogram.data[0][UMIc]++;
            histogram.data[1][UMItc]++;
          }
          if (curcol.getSize() < genethresthold){
            if (auto ite = curcol()) do{
              filt_data[1].data[curcol.getSize()-1][ite()] += *ite;
            }while(ite++);
            histogram.data[4][curcol.getSize()-1]++;
          }
          curcol.toMemfree();
        }
        if (UMItc_max < UMItc) UMItc_max = UMItc;
        if (fscout != 3) break;
        old = buf[1]; UMIc =0; UMItc =0;
      }
      curcol[buf[0]-1] = buf[2]; // starts at 1 eh?
      if (toignore_index.find(buf[0]) == 0xFFFFFFFF) UMIc += buf[2];
      UMItc += buf[2];
    }
  }
	if (do_progress & 1) progbar.finish();
  fclose(g);
  uint32_t i,j,k;
  Tuple<double, 3u> ccc;



  myHashmap<uint32_t, void> tofilter2;
  Vector<double> soupvec, soupfvec;

  if (do_progress & 2){ // soup regression threshold
    // base threhold defines underbound for soup frequency of transcripts
    // now cells needs to pass for a nUMI nGENE for same thresholds, but now weighted down by fraction of SOUP!
    k = 0;
    while(true){
      if (auto ite = allunfiltered.mkIterator()) do{
        if ((i = allfiltered.find(ite())) == 0xFFFFFFFF) meatfraction[ite()] = 1.0;
        else meatfraction[ite()] = ((double)(*ite)) / ((*ite) + allfiltered.deref(i));
      }while(ite++);



      for(i=0;i<data.data.getSize();i++){
        if (tofilter2.find(i) != 0xFFFFFFFF) continue;
        ccc.toZero();
        if (auto ite = data.data[i]()) do{
	  if ((toignore_index.find(ite())) == 0xFFFFFFFF) { // must be a non-prior bg transcript to count!
            ccc[0] += meatfraction[ite()];
            ccc[1] += meatfraction[ite()] * (*ite);
	  }
        }while(ite++);
        if ((ccc[0] < genethresthold)||(ccc[1] < umithresthold)) tofilter2.addEntry(i);
      }
      if (k == tofilter2.getSize()) break;
      printf("adaptive got %i/%i more filtered\n", tofilter2.getSize(), data.data.getSize());

      // update soup fraction!
      for(i=k;i<tofilter2.getSize();i++){
        if (auto ite = data.data[tofilter2.deref(i)]()) do{
          allunfiltered[ite()] -= *ite;
          allfiltered[ite()] += *ite;
        }while(ite++);
      }
      k = tofilter2.getSize();
    }

    printf("adaptive got %i/%i more filtered\n", tofilter2.getSize(), data.data.getSize());

    for(i=0;i<data.data.getSize();i++){
      ccc.toZero();
      if (auto ite = data.data[i]()) do{
	if ((toignore_index.find(ite())) == 0xFFFFFFFF) { // must be a non-prior bg transcript to count!
          ccc[0] += meatfraction[ite()] * (*ite);
	  ccc[2] += meatfraction[ite()];
	}
        ccc[1] += (double)(*ite);
      }while(ite++);
      double tmp = 1.0 - (ccc[0] / ccc[1]);
      double tmps = soupfraction + souplog10slope * log(ccc[2]) / log(10);
      if (tofilter2.find(i) != 0xFFFFFFFF) soupfvec.push_back(tmp);
      else if (tmp < tmps) soupvec.push_back(tmp);
      else {tofilter2.addEntry(i); soupfvec.push_back(tmp);}
    }
    printf("adaptive got %i/%i filtered based on soup threshold\n", tofilter2.getSize(), data.data.getSize());
  }

  printf("soupvec %i\n", soupvec.getSize());
  printf("soupvecf %i\n", soupfvec.getSize());


  Rcpp::S4 fout2, fout3, foutf, foutsam;


  printf("Reading Cell Names..."); fflush(stdout);

  memcpy(buffer, path.c_str(), l);
  strcpy(buffer+l, "barcodes.tsv");
  f = fopen(buffer,"r");
  if (f == NULL) {fprintf(stderr, "could not open '%s'!\n", buffer); return Rcpp::List::create();}

  cellnames = Rcpp::CharacterVector(selcol.getSize() - tofilter2.getSize());
  filtcellnames = Rcpp::CharacterVector(tofilter2.getSize());
  for(buf[0] =0,buf[1] =0,k=0; 1 == (fscout = fscanf(f,"%[^\n]\n", buffer)); buf[0]++){
    if (selcol[buf[1]] == buf[0]) {
      if (tofilter2.find(buf[1]++) == 0xFFFFFFFF) cellnames[k++] = buffer;
      else filtcellnames[buf[1]-k-1] = buffer;
      if (buf[1] >= selcol.getSize()) break;
    }
  }
  fclose(f);

  printf("(Done) cells = %i\n", (int)cellnames.length()); fflush(stdout);
  if (cellnames.length() < 2) return Rcpp::List::create(); // ALL cells got filtered... or 1... R type conversion is annoying when only 1 cell is detected... sorry
  if (do_progress  & 16){
	  data.wrRcppdgCMatrix(fout, nbgenes);
	  fout.slot("Dimnames") = Rcpp::List::create(genenames,cellnames);
	return Rcpp::List::create(Rcpp::Named("data") = fout);
  }

  if (tofilter2.getSize() != 0){
    //remove data at 2nd step
    softfiltered_data.setNBcols(tofilter2.getSize());
    for(i=0; tofilter2.find(i) == 0xFFFFFFFF; i++) ;
    k=0;
    softfiltered_data.data[k++].toMemmove(data.data[i]);
    j=i;
    for(i++;i < data.data.getSize();i++){
      if (tofilter2.find(i) == 0xFFFFFFFF) data.data[j++].toMemmove(data.data[i]);
      else softfiltered_data.data[k++].toMemmove(data.data[i]);
    }
    data.data.DownSize(j);
    softfiltered_data.wrRcppdgCMatrix(foutf, nbgenes);
    foutf.slot("Dimnames") = Rcpp::List::create(genenames,filtcellnames);
  }else if (do_progress & 2){
     softfiltered_data.setNBcols(0);
     softfiltered_data.wrRcppdgCMatrix(foutf, nbgenes);
     foutf.slot("Dimnames") = Rcpp::List::create(genenames,filtcellnames);
  }
  printf("charveclen: %i %i %i\n",(int)genenames.length(), (int)cellnames.length(),(int)filtcellnames.length() );
  printf("souplen %i %i\n", soupvec.getSize() ,  soupfvec.getSize() ) ; fflush(stdout);

  data.wrRcppdgCMatrix(fout, nbgenes);
  fout.slot("Dimnames") = Rcpp::List::create(genenames,cellnames);

  if (filt_data[2].data.getSize() > 0){
 	 Rcpp::CharacterVector samnames = Rcpp::CharacterVector(filt_data[2].data.getSize());
	 for(fscout=0;fscout<filt_data[2].data.getSize();fscout++){sprintf(buffer,"%i",fscout); samnames[fscout] = buffer;}
	 filt_data[2].wrRcppdgCMatrix(foutsam, nbgenes);
	 foutsam.slot("Dimnames") = Rcpp::List::create(genenames,samnames);
  }


  Vector<double> soupgene;
  Rcpp::List fouttt;
  if (umithresthold > 0){
    filt_data[0].wrRcppdgCMatrix(fout2, nbgenes);
    Rcpp::CharacterVector umicnames = Rcpp::CharacterVector(umithresthold);
    for(fscout=0;fscout<umithresthold ;fscout++){sprintf(buffer,"%i",fscout); umicnames[fscout] = buffer;}
    fout2.slot("Dimnames") = Rcpp::List::create(genenames,umicnames);
  }
  if (do_progress & 2){
    histogram.wrRcppdgCMatrix(fout3, UMItc_max);
    Rcpp::CharacterVector historow,histocol;
    histocol.push_back("Filteredcells_UMI_total_count");
    histocol.push_back("Filteredcells_unfilteredgene_UMI_total_count");
    histocol.push_back("UMI_total_count");
    histocol.push_back("unfilteredgene_UMI_total_count");
    histocol.push_back("Filteredcells_GENE_coverage");
     historow = Rcpp::CharacterVector(UMItc_max);
    for(fscout=0;fscout<UMItc_max ;fscout++){sprintf(buffer,"%i",fscout+1); historow[fscout] = buffer;}
    fout3.slot("Dimnames") = Rcpp::List::create(historow,histocol);
    soupgene.setSize(genenames.length());
    for(i=0;i<genenames.length();i++){
      if ((j = meatfraction.find(i)) == 0xFFFFFFFF) soupgene[i] = 1.0;
      else soupgene[i] = 1.0 - meatfraction.deref(j);
    }
    Rcpp::NumericVector fout_soupgene; soupgene.wrRCPP(fout_soupgene);  fout_soupgene.names() = genenames;
    Rcpp::NumericVector fout_soup; soupvec.wrRCPP(fout_soup); fout_soup.names() = cellnames;
    Rcpp::NumericVector foutf_soup; soupfvec.wrRCPP(foutf_soup); foutf_soup.names() = filtcellnames;

    fouttt = Rcpp::List::create(Rcpp::Named("data") = fout,
    Rcpp::Named("filtered.data") = foutf,
    Rcpp::Named("gene_accession") = Rcpp::wrap(geneacc_names),
    Rcpp::Named("umi_histogram") = fout3,
    Rcpp::Named("soupfraction_cell") = fout_soup,
    Rcpp::Named("soupfraction_filteredcell")= foutf_soup,
    Rcpp::Named("soupfraction_gene") = fout_soupgene);
  }else if (umithresthold > 1){
    fouttt = Rcpp::List::create(Rcpp::Named("data") = fout,
    Rcpp::Named("filtered.data") = foutf,
    Rcpp::Named("gene_accession") = Rcpp::wrap(geneacc_names),
    Rcpp::Named("umi_histogram") = fout3);
  }else{
    fouttt = Rcpp::List::create(Rcpp::Named("data") = fout,
    Rcpp::Named("filtered.data") = foutf,
    Rcpp::Named("gene_accession") = Rcpp::wrap(geneacc_names));
  }
  if (filt_data[2].data.getSize() > 0) fouttt["soup.samples"] = foutsam;
  if (umithresthold > 0) fouttt["hardfiltered.UMIs"] = fout2;

return fouttt;}



//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_convertDenseTextMatrixToSparseFormat)]]
void infernal_convertDenseTextMatrixToSparseFormat(Rcpp::List list){
  char buffer[512];
  std::vector< std::string >  paths = Rcpp::as< std::vector< std::string > >(list[0]);
  std::string outprefix = Rcpp::as< std::string >(list[1]);
  std::vector< std::string >  nameprefix = Rcpp::as< std::vector< std::string > >(list[2]);
  int flags = Rcpp::as< int >(list[3]);
  SparseMatrix<double> sp;
  SparseMatrix<double> sptmp;
  Vector<string> rownames,colnames, curcolnames;
  string tnamebuf;
  uint32_t fite;
  for(fite=0;fite<paths.size();fite++){
    sptmp.readDenseTable(paths[fite].c_str(),rownames,curcolnames, tnamebuf, "do.progress", flags & 1, ((flags & 2) ? ',' : '\t'));
    if (nameprefix[fite].length() != 0){
      for(uint32_t i = 0; i < curcolnames.size();i++) curcolnames[i] = nameprefix[fite] + curcolnames[i];
    }
    if (fite == 0){
      sp.toMemmove(sptmp);
      colnames.toMemmove(curcolnames);
    }else{
      sp.toColMemappend(sptmp);
      colnames.toMemappend(curcolnames);
    }
  }
  sp.writeMTXTable(outprefix.c_str(),rownames,colnames, tnamebuf,nullptr, flags & 1,true);
  return;
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_exportDataAsMtx)]]
void infernal_exportDataAsMtx(SEXP infscope, Rcpp::List list){
  Rcpp::XPtr< Infern0scope > scptr(infscope);
  char buffer[512];
  std::string pathprefix = Rcpp::as< std::string >(list[0]);
  std::string matrixname = Rcpp::as< std::string >(list[1]);
  int flags = Rcpp::as< int >(list[2]);

  Vector<string> rownames = Vector<string>(scptr->annot.dicos[INFRN0_AXES_GENE_NAMES].entries);
  Vector<string> colnames = Vector<string>(scptr->annot.dicos[INFRN0_AXES_CELL_NAMES].entries);
  printf("Exporting data to %s with %i genes and %i cells\n",pathprefix.c_str(), rownames.getSize(),colnames.getSize());
  if (flags & 1) scptr->data.writeMTXTable(pathprefix.c_str(),rownames,colnames, matrixname,nullptr, true,true);
  else scptr->rawdata.writeMTXTable(pathprefix.c_str(),rownames,colnames, matrixname,nullptr, false,true);
}


void reverseTCM(SparseMatrix<uint32_t> &targ, const SparseMatrix<double> &sp, bool logtr){
  uint32_t i,j,k,l;
  double sum, smallest;
  double ssmallest[2];
  targ.setNBcols(sp.getNBcols());
  for(j=0;j<sp.getNBcols();j++){
    targ.data[j].toMemfree();
    if (sp.data[j].getSize() == 0) continue;
    smallest = sp.data[j].deref(0);
    for(i=1;i<sp.data[j].getSize();i++) if (sp.data[j].deref(i) != smallest) break;
    if (i == sp.data[j].getSize()){ // all equal... assumes they are 1s!
      for(i=0;i<sp.data[j].getSize();i++) targ.data[j][sp.data[j].deref_key(i)] = 1.0;
      continue;
    }
    if (sp.data[j].deref(i) < smallest){
      ssmallest[0] = smallest; smallest = sp.data[j].deref(i);
    }else ssmallest[0] = sp.data[j].deref(i);

    for(i++;i<sp.data[j].getSize();i++) if ((sp.data[j].deref(i) != smallest)&&(sp.data[j].deref(i) != ssmallest[0])) break;

    if (i == sp.data[j].getSize()) ssmallest[1] = ssmallest[0]; // only 2 value... so with them!
    else{

         if (sp.data[j].deref(i) > ssmallest[0]) ssmallest[1] = sp.data[j].deref(i);
	 else{
	        ssmallest[1] = ssmallest[0];
		if (sp.data[j].deref(i) >= smallest) ssmallest[0] = sp.data[j].deref(i);
	        else {ssmallest[0] = smallest; smallest = sp.data[j].deref(i);}
	 }
	for(i++;i<sp.data[j].getSize();i++) {
         if (sp.data[j].deref(i) >= ssmallest[1]) continue;
         if (sp.data[j].deref(i) > ssmallest[0]) ssmallest[1] = sp.data[j].deref(i);
	 else{
	        ssmallest[1] = ssmallest[0];
		if (sp.data[j].deref(i) >= smallest) ssmallest[0] = sp.data[j].deref(i);
	        else {ssmallest[0] = smallest; smallest = sp.data[j].deref(i);}
	 }
      }
    }

    if (logtr) {smallest = exp(smallest); ssmallest[0] = exp(ssmallest[0]);ssmallest[1] = exp(ssmallest[1]);}

    ssmallest[0] /= smallest;
    ssmallest[1] /= smallest;
     for(l=1;l<1000; l++){
	if (fabs(ssmallest[0] * l - floor((ssmallest[0] * l) + 0.5)) + fabs(ssmallest[1] * l - floor((ssmallest[1] * l) + 0.5)) < 0.01) break;
     }
     if (l == 1000) {
	if ((j % 100) == 0){
 	    for(l=1;l<20; l++){
		printf("%i:%e and %e\n", l, ssmallest[0] * l - floor((ssmallest[0] * l) + 0.5),  ssmallest[1] * l - floor((ssmallest[1] * l) + 0.5));
     	    }
	    l = 1000;
	}
	printf("warning!, could not figure denominator for cell no. '%i' (2 smallest ratios %e and %e) using 1000\n",j, ssmallest[0], ssmallest[1]);
	}

    if (logtr){
      ssmallest[0] = log((double)l) - log(smallest);
      for(i=0;i<sp.data[j].getSize();i++) targ.data[j][sp.data[j].deref_key(i)] = floor(exp(sp.data[j].deref(i) + ssmallest[0])  + 0.5);
    }else{
      ssmallest[0] = ((double)l) / smallest;
      for(i=0;i<sp.data[j].getSize();i++) targ.data[j][sp.data[j].deref_key(i)] = floor(sp.data[j].deref(i) * ssmallest[0] + 0.5);
    }
  }
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_reverseTCM)]]
void infernal_reverseTCM(SEXP rscp, Rcpp::List list){
  char buffer[512];
  uint32_t flag = Rcpp::as< uint32_t >(list[0]);
  Rcpp::XPtr< Infern0scope > scp(rscp);
  reverseTCM(scp->rawdata, scp->data, (flag & 1) == 1);
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_reverseTCM_inMTX)]]
void infernal_reverseTCM_inMTX(Rcpp::List list){
  char buffer[512];
  uint32_t flag = Rcpp::as< uint32_t >(list[2]);
  std::string inprefix = Rcpp::as< std::string >(list[0]);
  std::string outprefix = Rcpp::as< std::string >(list[1]);

  SparseMatrix<double> sp;
  SparseMatrix<uint32_t> spi;
  Vector<string> rownames,colnames;
  string tnamebuf;
  sp.readMTXTable(inprefix.c_str(),rownames,colnames, tnamebuf, true, true);
  reverseTCM(spi, sp, (flag & 1) == 1 );
  spi.writeMTXTable(outprefix.c_str(),rownames,colnames, tnamebuf, nullptr, false, true);
}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_readSparseMatrix)]]
void infernal_readSparseMatrix(SEXP scope, Rcpp::S4 x, bool isRaw){ Rcpp::XPtr< Infern0scope > scp(scope);
  if (isRaw) scp->rawdata.rdRcppdgCMatrix(x);
  else scp->data.rdRcppdgCMatrix(x);
  scp->annot.rdAxesNamesAl(x,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES);

}
//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_readDenseMatrix)]]
void infernal_readDenseMatrix(SEXP scope, Rcpp::NumericMatrix x, bool isRaw){ Rcpp::XPtr< Infern0scope > scp(scope);
	if (isRaw) scp->rawdata.rdRcpp(x);
  else scp->data.rdRcpp(x);
  scp->annot.rdAxesNamesAl(x,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES);
}
//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_createOutput)]]
Rcpp::List infernal_createOutput(SEXP data){
  Rcpp::XPtr< Infern0scope > scp(data);
  scp->data.show();
  funtest = funtest +1;
  return(Rcpp::List::create(Rcpp::Named("output") = Rcpp::wrap((double)666.0f + funtest)));
}
//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_fitNegativeBinomial)]]
Rcpp::List infernal_fitNegativeBinomial(SEXP data){ Rcpp::XPtr< Infern0scope > scp(data);
  //scp->fitNegativeBinomial();
  return(Rcpp::List::create(Rcpp::Named("output") = Rcpp::wrap((double)666.0f + scp->sillyint)));
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_getParameterList)]]
Rcpp::CharacterVector infernal_getParameterList(SEXP rscp){
  Rcpp::CharacterVector fout;
  Rcpp::XPtr< Infern0scope > scp(rscp);
  if (scp->rawdata.getNBcols() != 0)  {fout.push_back("cell.coverage"); fout.push_back("cell.nbcounts"); fout.push_back("gene.coverage"); fout.push_back("gene.nbcounts"); fout.push_back("raw.data");  }
  if (scp->data.getNBcols() != 0)  fout.push_back("data");
  if (scp->cell_ordering.getSize() != 0) fout.push_back("cell.order");
  if (scp->gene_ordering.getSize() != 0) fout.push_back("gene.order");
  if (scp->cell_clusterID.getSize() != 0) fout.push_back("cell.cluster");
  if (scp->gene_clusterID.getSize() != 0) fout.push_back("gene.cluster");
  if (scp->annot.hasDico((int)INFRN0_AXES_CELL_NAMES)) fout.push_back("cell.names");
  if (scp->annot.hasDico((int)INFRN0_AXES_GENE_NAMES)) fout.push_back("gene.names");
  if (scp->cell_scale.getSize() != 0) fout.push_back("cell.scale");
  if (scp->gene_scale.getSize() != 0) fout.push_back("gene.scale");
  if (scp->cell_state.getNBcols() != 0) {fout.push_back("cell.state"); fout.push_back("cell.state.transition");}
  if (scp->gene_state.getNBcols() != 0) {fout.push_back("gene.state"); fout.push_back("gene.state.transition");}
return(fout);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_setParameter)]]
uint32_t infernal_setParameter(SEXP rscp, uint32_t which, SEXP value){ Rcpp::XPtr< Infern0scope > scp(rscp);
  char buffer[256];
  std::vector<std::string> names;
  uint32_t i;
  switch((INFRN0_enum)which){
    case INFRN0_CELL_CLUSTER:{
      int total=0;
      Vector<uint32_t> groupId_vec(Rcpp::as<std::vector<uint32_t> >(value));
      for(i=0;i<groupId_vec.getSize();i++) if (groupId_vec[i] >= total) total = groupId_vec[i]+1;
      scp->cell_clusterID.toMemmove(groupId_vec);


      for(i=0;i< total;i++) {sprintf(buffer, "cluster%i", i); names.push_back(std::string(buffer));}
      Rcpp::CharacterVector genenames = Rcpp::wrap(names);
      scp->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_CELL_CLUSTER_NAMES);
      names.clear();

    }return 0;
    case INFRN0_GENE_CLUSTER:{
      int total=0;
      Vector<uint32_t> groupId_vec(Rcpp::as<std::vector<uint32_t> >(value));
      for(i=0;i<groupId_vec.getSize();i++) if (groupId_vec[i] >= total) total = groupId_vec[i]+1;
      scp->gene_clusterID.toMemmove(groupId_vec);

      for(i=0;i< total;i++) {sprintf(buffer, "cluster%i", i); names.push_back(std::string(buffer));}
      Rcpp::CharacterVector genenames = Rcpp::wrap(names);
      scp->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_GENE_CLUSTER_NAMES);
      names.clear();
    }return 0;
    case INFRN0_CELL_ORDER:{
      int total=0;
      Vector<uint32_t> groupId_vec(Rcpp::as<std::vector<uint32_t> >(value));
      for(i=0;i<groupId_vec.getSize();i++) if (groupId_vec[i] >= total) total = groupId_vec[i]+1;
      scp->cell_ordering.toMemmove(groupId_vec);
    }return 0;
    case INFRN0_GENE_ORDER:{
      int total=0;
      Vector<uint32_t> groupId_vec(Rcpp::as<std::vector<uint32_t> >(value));
      for(i=0;i<groupId_vec.getSize();i++) if (groupId_vec[i] >= total) total = groupId_vec[i]+1;
      scp->gene_ordering.toMemmove(groupId_vec);
    }return 0;
    case INFRN0_CELL_NAMES:{
      Rcpp::CharacterVector groupId_vec = Rcpp::as<Rcpp::CharacterVector >(value);
      scp->annot.setAxeNamesAl(groupId_vec,(int)INFRN0_AXES_CELL_NAMES);
    }return 0;
    case INFRN0_GENE_NAMES:{
      Rcpp::CharacterVector groupId_vec = Rcpp::as<Rcpp::CharacterVector >(value);
      scp->annot.setAxeNamesAl(groupId_vec,(int)INFRN0_AXES_GENE_NAMES);
    }return 0;
  }
return 0;}
//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_getParameter)]]
SEXP infernal_getParameter(SEXP rscp, uint32_t which){Rcpp::XPtr< Infern0scope > scp(rscp);
  std::vector<uint32_t> groupId;
  Rcpp::S4 fout;
  SEXP foute;

  uint32_t i,j;
  switch((INFRN0_enum)which){
  case INFRN0_CELL_CLUSTER:{
    ExOp::toConvert(groupId,scp->cell_clusterID);
    return Rcpp::wrap(groupId);
  }break;
  case INFRN0_GENE_CLUSTER:{
    ExOp::toConvert(groupId,scp->gene_clusterID);
    return Rcpp::wrap(groupId);
  }break;
  case INFRN0_CELL_SCALE:{
    Rcpp::NumericMatrix target(scp->cell_scale.getSize(),3);
    for(i=0;i<scp->cell_scale.getSize();i++){
      target(i,0) = scp->cell_scale[i][0];
      target(i,1) = scp->cell_scale[i][1];
      target(i,2) = scp->cell_scale[i][2];
    }
    colnames(target) = Rcpp::CharacterVector::create("NB_Rscale", "NB_Pscale", "Dropout_scale");
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_CELL_NAMES);
    rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_GENE_SCALE:{
    Rcpp::NumericMatrix target(scp->gene_scale.getSize(),3);
    for(i=0;i<scp->gene_scale.getSize();i++){
      target(i,0) = scp->gene_scale[i][0];
      target(i,1) = scp->gene_scale[i][1];
      target(i,2) = scp->gene_scale[i][2];
    }
    colnames(target) = Rcpp::CharacterVector::create("NB_Rscale", "NB_Pscale", "Dropout_scale");
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_GENE_NAMES);
    rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_CELL_STATE: scp->cell_state.wrRcppdgCMatrix(fout, scp->annot.getDicoSize(INFRN0_AXES_CELL_STATE_NAMES) ); scp->annot.wrAxesNamesAl(fout,  (int)INFRN0_AXES_CELL_STATE_NAMES,  (int)INFRN0_AXES_CELL_NAMES); return Rcpp::wrap(fout);
  case INFRN0_GENE_STATE: scp->gene_state.wrRcppdgCMatrix(fout, scp->annot.getDicoSize(INFRN0_AXES_GENE_STATE_NAMES) ); scp->annot.wrAxesNamesAl(fout,  (int)INFRN0_AXES_GENE_STATE_NAMES,  (int)INFRN0_AXES_GENE_NAMES); return Rcpp::wrap(fout);
  case INFRN0_CELL_NAMES:{
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_CELL_NAMES);
    return Rcpp::wrap(vec);
  }break;
  case INFRN0_GENE_NAMES:{
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_GENE_NAMES);
    return Rcpp::wrap(vec);
  }break;
  case INFRN0_GENE_COVERAGE:{
    Rcpp::NumericMatrix target(scp->rawdata.getNBcols(),1);
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_CELL_NAMES);
    for(i=0;i<scp->rawdata.getNBcols();i++) target(i,0) = scp->rawdata.data[i].getSize();
    Rcpp::CharacterVector cvec =  Rcpp::CharacterVector(); cvec.push_back("Gene Coverage");
    Rcpp::colnames(target) = cvec;
    Rcpp::rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_CELL_COVERAGE:{
    Vector<uint32_t> row_size = scp->rawdata.compute_rowsizes();
    Rcpp::NumericMatrix target(row_size.getSize(),1);
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_GENE_NAMES);
    for(i=0;i<row_size.getSize();i++) target(i,0) = row_size[i];
    Rcpp::CharacterVector cvec =  Rcpp::CharacterVector(); cvec.push_back("Cell Coverage");
    Rcpp::colnames(target) = cvec;
    Rcpp::rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_CELL_SIZE:{
    Rcpp::NumericMatrix target(scp->rawdata.getNBcols(),1);
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_CELL_NAMES);
    for(i=0;i<scp->rawdata.getNBcols();i++) {
	if (auto ite = scp->rawdata.data[i]()) {
		target(i,0) = *ite;
		while(ite++) target(i,0) += *ite;
	}
    }
    Rcpp::CharacterVector cvec =  Rcpp::CharacterVector(); cvec.push_back("Cell nbcounts");
    Rcpp::colnames(target) = cvec;
    Rcpp::rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_GENE_SIZE:{
    int nbgene = scp->getNBgenes();
    Rcpp::NumericMatrix target(nbgene,1);
    Rcpp::CharacterVector vec = Rcpp::CharacterVector();
    scp->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_GENE_NAMES);
    for(i=0;i<nbgene;i++) target(i,0) = 0;
    for(i=0;i<scp->rawdata.getNBcols();i++) {
	if (auto ite = scp->rawdata.data[i]()) {
		do{ target(ite(),0) += *ite;}while(ite++);
	}
    }
    Rcpp::CharacterVector cvec =  Rcpp::CharacterVector(); cvec.push_back("Gene nbcounts");
    Rcpp::colnames(target) = cvec;
    Rcpp::rownames(target) = vec;
    return Rcpp::wrap(target);
  }break;
  case INFRN0_CELL_ORDER:{
    ExOp::toConvert(groupId,scp->cell_ordering);
    return Rcpp::wrap(groupId);
  }break;
  case INFRN0_GENE_ORDER:{
    ExOp::toConvert(groupId,scp->gene_ordering);
    return Rcpp::wrap(groupId);
  }break;
  case INFRN0_RAW_DATA: scp->rawdata.wrRcppdgCMatrix(fout, scp->annot.getDicoSize(INFRN0_AXES_GENE_NAMES) ); scp->annot.wrAxesNamesAl(fout,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES); return Rcpp::wrap(fout);
  case INFRN0_DATA: scp->data.wrRcppdgCMatrix(fout, scp->annot.getDicoSize(INFRN0_AXES_GENE_NAMES) ); scp->annot.wrAxesNamesAl(fout,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES); return Rcpp::wrap(fout);
  case INFRN0_CELL_STATE_MEAN_DEVIATION:
    /*q()Rcpp::NumericMatrix target(scp->gene_scale.getSize(), scp->par_R.sizes[1] );
    Tuple<double> denum; denum.setSize(scp->par_R.sizes[1]);
    for(i=0;i<scp->gene_scale.getSize();i++){
      for(j=0;j<scp->par_R.sizes[1];j++) target(i,j) = 0;
      denum.toZero();
      if (auto ite = scp->data.data[i]()) do{
          if (auto ite2 <- scp->cell_state.data[ite()]()) do{
            target(i,ite2()) += (*ite2) * (*ite);
            denum += (*ite2);
          }while(ite2++);
      }while(ite++);
      for(j=0;j<scp->par_R.sizes[0];j++) target(i,j) /= denum[j];
    }
    scp->annot.wrAxesNamesAl(target,  (int)INFRN0_AXES_CELL_NAMES,  (int)INFRN0_AXES_GENE_STATE_NAMES);
    return Rcpp::wrap(target);*/
  break;
  case INFRN0_GENE_STATE_MEAN_DEVIATION:{
    Rcpp::NumericMatrix target(scp->cell_scale.getSize(), scp->par_R.sizes[0] );
    Tuple<double> denum; denum.setSize(scp->par_R.sizes[0]);
    for(i=0;i<scp->cell_scale.getSize();i++){
      for(j=0;j<scp->par_R.sizes[0];j++) target(i,j) = 0;
      denum.toZero();
      if (auto ite = scp->data.data[i]()) do{
          if (auto ite2 = scp->gene_state.data[ite()]()) do{
            target(i,ite2()) += (*ite2) * (*ite);
            denum += (*ite2);
          }while(ite2++); }while(ite++);
      for(j=0;j<scp->par_R.sizes[0];j++) {
        target(i,j) /= denum[j];
        if (!ExOp::isValid(target(i,j))) target(i,j) =0.0;
      }
    }
    scp->annot.wrAxesNamesAl(target,  (int)INFRN0_AXES_CELL_NAMES,  (int)INFRN0_AXES_GENE_STATE_NAMES);
    return Rcpp::wrap(target);
  }
  case INFRN0_GENE_STATE_MATRICES:{
    Rcpp::NumericMatrix matR(scp->par_R.sizes[0],scp->par_R.sizes[1]);
    Rcpp::NumericMatrix matM(scp->par_M.sizes[0],scp->par_M.sizes[1]);
    Rcpp::NumericMatrix matD(scp->par_D.sizes[0],scp->par_D.sizes[1]);
    scp->par_R.wrMatrix(matR);
    scp->par_M.wrMatrix(matM);
    scp->par_D.wrMatrix(matD);
    scp->annot.wrAxesNamesAl(matR,  (int)INFRN0_AXES_GENE_STATE_NAMES,  (int)INFRN0_AXES_CELL_STATE_NAMES);
    scp->annot.wrAxesNamesAl(matM,  (int)INFRN0_AXES_GENE_STATE_NAMES,  (int)INFRN0_AXES_CELL_STATE_NAMES);
    scp->annot.wrAxesNamesAl(matD,  (int)INFRN0_AXES_GENE_STATE_NAMES,  (int)INFRN0_AXES_CELL_STATE_NAMES);
    return Rcpp::List::create(Rcpp::Named("R")= matR,Rcpp::Named("M")= matM,Rcpp::Named("D")= matD);
  }break;
	case INFRN0_PVAL_DATA:scp->pval_data.wrRcppdgCMatrix(fout); scp->annot.wrAxesNamesAl(fout,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES); return Rcpp::wrap(fout);
	case INFRN0_BIODROP_PROB:{
    Rcpp::NumericMatrix hugebio(scp->gene_scale.getSize(),scp->cell_scale.getSize());
    for(uint32_t col=0;col< scp->cell_scale.getSize();col++){
      Tuple<double> fout = scp->getBiodropoutProb(col);
      for(uint32_t row=0;row< scp->cell_scale.getSize();row++) hugebio(row,col) = fout[row];
    }
    scp->annot.wrAxesNamesAl(hugebio,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES);
    return Rcpp::wrap(hugebio);
  }break;
  case INFRN0_CELL_TRANSITION:{
    Rcpp::NumericMatrix tr;
    Rcpp::NumericVector vc = Rcpp::NumericVector(scp->par_R.sizes[1]);
    LFHPrimitive::MPReachDistribution tolearn; tolearn.setDefaultRatesAndPrior(scp->par_R.sizes[1], 2.0);
    auto lrn = tolearn.mkEMscope();
    lrn.init();
    for(i=0;i<scp->cell_state.getNBcols();i++) lrn.EMregist(scp->cell_state.data[i]);
    tolearn.learn(lrn, 1.0);
    if (auto ite = tolearn.log_transition.getIterator()) do {*ite = exp(*ite);} while(ite++);
    if (auto ite = tolearn.log_start_prior.getIterator()) do {vc[ite()] = exp(*ite);} while(ite++);
    tolearn.log_transition.wrMatrix(tr);
    for(int i=0; i<scp->par_R.sizes[1];i++) for(int j= i+1;j < scp->par_R.sizes[1];j++) tr(i,j) =0;


    scp->annot.wrAxesNamesAl(tr,INFRN0_AXES_CELL_STATE_NAMES,INFRN0_AXES_CELL_STATE_NAMES);

     return Rcpp::List::create(Rcpp::Named("Start")= vc,Rcpp::Named("Transition")= tr);
  }break;
  case INFRN0_GENE_TRANSITION:{
    Rcpp::NumericMatrix tr;
    Rcpp::NumericVector vc = Rcpp::NumericVector(scp->par_R.sizes[0]);
    LFHPrimitive::MPReachDistribution tolearn; tolearn.setDefaultRatesAndPrior(scp->par_R.sizes[0], 2.0);
    auto lrn = tolearn.mkEMscope();
    lrn.init();
    for(i=0;i<scp->gene_state.getNBcols();i++) lrn.EMregist(scp->gene_state.data[i]);
    tolearn.learn(lrn, 1.0);
    if (auto ite = tolearn.log_transition.getIterator()) do {*ite = exp(*ite);} while(ite++);
    if (auto ite = tolearn.log_start_prior.getIterator()) do {vc[ite()] = exp(*ite);} while(ite++);
    tolearn.log_transition.wrMatrix(tr);
    scp->annot.wrAxesNamesAl(tr,INFRN0_AXES_GENE_STATE_NAMES,INFRN0_AXES_GENE_STATE_NAMES);
    for(int i=0; i<scp->par_R.sizes[0];i++) for(int j= i+1;j < scp->par_R.sizes[0];j++) tr(i,j) =0;

     return Rcpp::List::create(Rcpp::Named("Start")= vc,Rcpp::Named("Transition")= tr);
  } break;


default: printf("warning, item %X does not exist!", which);
  }
  return Rcpp::wrap(groupId);
}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_getDataRow)]]
SEXP infernal_getDataRow(SEXP rscp, Rcpp::List list){
	Rcpp::XPtr< Infern0scope > scp(rscp);
	std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
	std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
	uint32_t which = Rcpp::as< uint32_t >(list[2]);

	if (genelist[0].length() == 0) genelist.clear();
	if (celllist[0].length() == 0) celllist.clear();

	Rcpp::NumericMatrix fout;
	uint32_t ite;
	if (which == 0){
		if (genelist.size() == 1){
			SEXP vout;
			ite = scp->annot.dicos[INFRN0_AXES_GENE_NAMES].findEntry(genelist[0].c_str());
			if (ite != 0xFFFFFFFF)	{
				scp->rawdata.wrRcppRow(ite, vout);
			}else return R_NilValue;
			return vout;
		}else if (celllist.size() == 1){
			SEXP vout;
			ite = scp->annot.dicos[INFRN0_AXES_CELL_NAMES].findEntry(celllist[0].c_str());
			if (ite != 0xFFFFFFFF) {
				scp->rawdata.wrRcppCol(ite, vout, scp->getNBgenes());
			}else return R_NilValue;
			return vout;
		}else return R_NilValue;
	}else{
		myHashmap<uint32_t> rows;
		myHashmap<uint32_t> cols;
		int cnf =0;
		int gnf =0;
		if (genelist.size() == 0) for(int i=0; i < scp->getNBgenes();i++) rows.addEntry(i);
		else for(int i=0;i<genelist.size();i++) {
			ite = scp->annot.dicos[INFRN0_AXES_GENE_NAMES].findEntry(genelist[i].c_str());
			if (ite != 0xFFFFFFFF) rows.addEntry(ite);
			else gnf++;
		}
		if (celllist.size() == 0) for(int i=0; i < scp->getNBcells();i++) cols.addEntry(i);
		else for(int i=0;i<celllist.size();i++) {
			ite = scp->annot.dicos[INFRN0_AXES_CELL_NAMES].findEntry(celllist[i].c_str());
			if (ite != 0xFFFFFFFF) cols.addEntry(ite);
			else cnf++;
		}
		if (((gnf == genelist.size())&&(genelist.size() != 0))||((cnf == celllist.size())&&(celllist.size() != 0))) {printf("error: %i cells and %i genes queried were not found.\n", gnf,cnf); return R_NilValue;}
		if ((cnf != 0)||(gnf != 0)) printf("warning: %i cells and %i genes queried were not found.\n", gnf,cnf);
		scp->wrDeviationMatrix(fout, rows, cols, which - 1);
		return fout;
	}
}





//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_runTestFunction)]]
void infernal_runTestFunction(SEXP data){
  Rcpp::XPtr< Infern0scope > scp(data);
  scp->dothattest();
return;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_HierarchicalClustering)]]
void infernal_HierarchicalClustering(SEXP infscope, Rcpp::List list){
  Rcpp::XPtr< Infern0scope > scptr(infscope);
  bool do_cluster_cells = Rcpp::as<bool>(list[0]);
  bool do_overwrite = Rcpp::as<bool>(list[1]);
  scptr->clusterHierarchical(do_cluster_cells, do_overwrite);
return;}





//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_computePartialCorrelation)]]
Rcpp::List infernal_computePartialCorrelation(SEXP infscope, Rcpp::List list){
  Rcpp::XPtr< Infern0scope > scptr(infscope);
  std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
  std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
  uint32_t nbthreads = Rcpp::as< uint32_t >(list[2]);
  uint32_t nbout = Rcpp::as< uint32_t >(list[3]);
  uint32_t nbindimax = Rcpp::as< uint32_t >(list[4]);
  uint32_t mincell = Rcpp::as< uint32_t >(list[5]);
  std::vector<double> minpthr = Rcpp::as< std::vector<double> >(list[6]);


  Vector<uint32_t> cinds = scptr->annot.dicos[INFRN0_AXES_CELL_NAMES].findEntries(celllist);
  Vector<uint32_t> ginds = scptr->annot.dicos[INFRN0_AXES_GENE_NAMES].findEntries(genelist);

  TMatrix<double, 0u,0u> covout;
  Vector<uint32_t> gout; gout.setSize(nbout);

  scptr->runTaskPartcorrel(covout,gout,ginds,cinds, nbthreads, nbindimax, mincell, &(minpthr[0]));
  Rcpp::NumericMatrix cov;
  covout.wrMatrix(cov);
  colnames(cov) =  Rcpp::as< Rcpp::CharacterVector >(list[0]);
  Rcpp::CharacterVector vec = Rcpp::CharacterVector();
  scptr->annot.getAxeNamesAl(vec, (int)INFRN0_AXES_GENE_NAMES);
  rownames(cov) = vec;
  Rcpp::CharacterVector gvec = Rcpp::CharacterVector();
  for(uint32_t i=0;i< gout.getSize();i++) gvec.push_back(vec.at(gout[i]));
return Rcpp::List::create(Rcpp::Named("correlation") = cov,Rcpp::Named("gene.list") = gvec);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_IdentifyNetwork)]]
Rcpp::List infernal_IdentifyNetwork(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean){
  SparseMatrix<double> inputSP; inputSP.rdRcppdgCMatrix(input);
  Rcpp::List davec = input.slot("Dimnames");
  Rcpp::CharacterVector genenames = davec[0];
  Rcpp::CharacterVector cellnames = davec[1];
  printf("Run on %i x %i\n", (int)genenames.length(), (int)cellnames.length())  ; fflush(stdout);
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
  Tuple< Tuple < double> > delta_mean;

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
  // produces Suffitient statistics for test sets centered at each training set mean (and at global mean, actual sample variance)
  for(i =0; i< inputSP.getNBcols();i++){
    partcount[cell_partition[i]]++;
    if (auto ite = inputSP.data[i]()) do{
      k = (ite() * (ite() +1)) >> 1;
      tmp = (*ite) - dameans[cell_partition[i]][ite()];
      tmp2 = (*ite) - mean[ite()];
      if (auto ite2 = inputSP.data[i]()) do{
        if (ite2()<=ite()){
          testset[cell_partition[i]].data[k+ite2()] += tmp * ((*ite2) - dameans[cell_partition[i]][ite2()]);
          wholetestset.data[k+ite2()] += tmp2 * ((*ite2) - mean[ite2()]);
        }
      }while(ite2++);
    }while(ite++);
  }
  Vector< KeyElem<Tuple<uint32_t, 2u>, Tuple<double, 13u> > > details;
  SparseTrianglix<double> precision;
  printf("Run\n")  ; fflush(stdout);
  Tuple<double> eigen;
  for(i=0;i<testset.getSize();i++){
    eigen = training[i].getEigenValues();
    for(j=0;j<eigen.getSize();j++) if (eigen[j] <= 0.0) break;
    if (j != eigen.getSize()) {
	printf("Training set%i is not positive definite! Eigenvalues:\n",i);eigen.show();
        return R_NilValue;
    }
  }

  //srtsk.run(precision,tb,training, testset, details, nbedges, target.getSize() == mean.length() ? &target : NULL, maxcyclic, mincheck, (flags & 1));
  foreach.startThreadArray(nbthreads);
  precision.searchRegul2019(totalcovar,training, testset, partcount, details, nbedges, target.getSize() == mean.length() ? &target : NULL, maxcyclic, mincheck);
  foreach.stopThreadArray();
  arma::mat crdmat = arma::mat(details.getSize(),2);
  for(int i =0;i<details.getSize();i++){crdmat.at(i,0) = details[i].k[0];crdmat.at(i,1) = details[i].k[1];}
  Rcpp::NumericMatrix prec, suff;
  precision.wrMatrix(prec);
  wholetestset.wrMatrix(suff);
  Rcpp::colnames(prec) = genenames;
  Rcpp::rownames(prec) = genenames;
  Rcpp::colnames(suff) = genenames;
  Rcpp::rownames(suff) = genenames;

  Rcpp::NumericMatrix llincrmat(details.getSize(),13);
  for(int i =0;i<details.getSize();i++){
    for(int j=0;j<13;j++) llincrmat(i,j) = details[i].d[j];
  }
  if (target.getSize() == mean.length()){
	  Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("Log_Determinent","TrainSet_LL","TestSet_LL","CrxSum_Log_det","CrXSum_Trace","CrXSum_gradient","CrX_LLstdev","Error","CrX_error","TPrate","FPrate", "InsIntoSize", "TimeMilli");
  }else{
	  Rcpp::colnames(llincrmat) = Rcpp::CharacterVector::create("Log_Determinent","TrainSet_LL","TestSet_LL","CrxSum_Log_det","CrXSum_Trace","CrXSum_gradient","CrX_LLstdev","Error","CrX_error","NoData1","NoData2", "InsIntoSize", "TimeMilli");
  }
return Rcpp::List::create(Rcpp::Named("PrecisionMatrix") = prec, Rcpp::Named("SuffitientMatrix") = suff, Rcpp::Named("EdgeListCoor") = crdmat,	Rcpp::Named("LL_incrs") = llincrmat,Rcpp::Named("nbedges") = nbedges);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_IdentifyConstrainedCovar)]]
Rcpp::List infernal_IdentifyConstrainedCovar(Rcpp::S4 input, Rcpp::List list, Rcpp::List traincov, Rcpp::List trainmean){
  SparseMatrix<double> inputSP; inputSP.rdRcppdgCMatrix(input);

  Rcpp::List davec = input.slot("Dimnames");
  Rcpp::CharacterVector genenames = davec[0];
  Rcpp::CharacterVector cellnames = davec[1];

  Rcpp::NumericMatrix covar = Rcpp::as< Rcpp::NumericMatrix >(list[0]);
  Rcpp::NumericVector mean = Rcpp::as< Rcpp::NumericVector >(list[1]);
  Rcpp::Vector<INTSXP> cell_partition = Rcpp::as< Rcpp::Vector<INTSXP> >(list[2]);
  Rcpp::NumericMatrix constraint = Rcpp::as< Rcpp::NumericMatrix >(list[3]);


  Tuple<Trianglix<double> > training; training.setSize(traincov.length());

  uint32_t i,k;
  Trianglix<double> totalcovar; totalcovar.rdMatrix(covar);
  Trianglix<char> constr; constr.rdMatrix(constraint);
  for(i =0; i< training.getSize();i++) training[i].rdMatrix(Rcpp::as< Rcpp::NumericMatrix >(traincov[i]));

  Trianglix<double> wholetestset; wholetestset.setSize(totalcovar.getSize()).toZero();

  Tuple<Trianglix<double> > testset; testset.setSize(traincov.length());
  Tuple< Rcpp::NumericVector  > dameans; dameans.setSize(traincov.length());

  for(i =0; i< training.getSize();i++) {
    testset[i].setSize(totalcovar.getSize()).toZero();
    dameans[i] = Rcpp::as<Rcpp::NumericVector>(trainmean[i]);
  }
  double tmp, tmp2;
  for(i =0; i< inputSP.getNBcols();i++){
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

	Tuple<double, 11u> details;
  SparseTrianglix<double> precision;
  SparseTrianglix<double>* train_precision = (traincov.length() == 0) ? NULL : new SparseTrianglix<double>[traincov.length()];
  precision.evalRegul(totalcovar,constr,training,testset, details,train_precision, NULL);

  Rcpp::NumericMatrix prec, suff;
  precision.wrMatrix(prec);
  wholetestset.wrMatrix(suff);
  Rcpp::colnames(prec) = genenames;
  Rcpp::rownames(prec) = genenames;
  Rcpp::colnames(suff) = genenames;
  Rcpp::rownames(suff) = genenames;
  arma::mat outprec;
  Rcpp::List crosslist, testsuffitient;
  crosslist = Rcpp::List(traincov.length());
  testsuffitient = Rcpp::List(traincov.length());
  for(i=0;i< traincov.length();i++){
    train_precision[i].wrMatrix(outprec);
    crosslist[i] = outprec;
   // Rcpp::colnames(crosslist[i]) = genenames;
   // Rcpp::rownames(crosslist[i]) = genenames;
    testset[i].wrMatrix(outprec);
    testsuffitient[i] = outprec;
   // Rcpp::colnames(testsuffitient[i]) = genenames;
    //Rcpp::rownames(testsuffitient[i]) = genenames;
  }

 	return Rcpp::List::create(Rcpp::Named("PrecisionMatrix") = prec, Rcpp::Named("SuffitientMatrix") = suff, Rcpp::Named("TrainingPrecisionMatrices") = crosslist, Rcpp::Named("TestSetSuffitientMatrices") = testsuffitient);
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_saveHierarchical)]]
void infernal_saveHierarchical(SEXP data, SEXP path){ Rcpp::XPtr< Infern0scope > scp(data);
  std::string pathvar = Rcpp::as<std::string>(path);
  scp->saveHierarchical(pathvar.c_str());
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_show)]]
void infernal_show(SEXP data){ Rcpp::XPtr< Infern0scope > scp(data);
  printf("Summary: Raw data: %i genes %i cells\n", scp->rawdata.computeNBrows(), scp->rawdata.getNBcols());
}


Rcpp::Vector<INTSXP> makeRandomPartition(uint32_t nbpairs, int nbcols, uint32_t seed){Rcpp::Vector<INTSXP> fout(nbcols);
  if (seed != 0xFFFFFFFF) srand(seed);
  uint32_t i,j,k;
  Tuple<uint32_t> randomizer;
  randomizer.setSize(nbpairs);
  for(i=0;i< nbpairs;i++) randomizer[i] = i;
  for(i=0;i< nbcols;i++){
    j = nbpairs - 1 - (i % nbpairs);
    if (j == 0) fout[i] = randomizer[0];
    else{
      k = rand() % (j + 1);
      fout[i] = randomizer[k];
      if (k != j) ExOp::toMemswap(randomizer[k], randomizer[j]);
    }
  }
  // not yet clever...
return(fout);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_computeCovar)]]
Rcpp::List infernal_computeCovar(SEXP scptrexp, Rcpp::List list){
  std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
  std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
  uint32_t nbpairs = Rcpp::as< uint32_t >(list[2]);
  uint32_t method = Rcpp::as< uint32_t >(list[3]); // 0 or 1
  uint32_t seed = Rcpp::as< uint32_t >(list[4]); // 0 or 1

  uint32_t inputClass=0;
  if (Rf_isMatrix(scptrexp)) inputClass=1;
  else if (Rf_isS4(scptrexp)) inputClass=2;

  Rcpp::CharacterVector axes[2];
  SparseMatrix<double> data;
  uint32_t nb_gene;
  switch(inputClass){
    case 0:{
    Rcpp::XPtr< Infern0scope > scptr(scptrexp);
    scptr->wrSubMatrix(data, genelist, celllist, axes, 3);
    nb_gene = axes[0].size();
    }break;
    case 2:{
    Rcpp::S4 dadata(scptrexp);
    data.rdRcppdgCMatrix(dadata);
    nb_gene = data.computeNBrows();
    Rcpp::List davec = dadata.slot("Dimnames");
    axes[0] = davec[0];
    axes[1] = davec[1];
    }break;
  }


  uint32_t col,i,j,k,l;


	Trianglix< Tuple<double, 4u> > target_w; target_w.setSize(nb_gene).toZero();
	Trianglix< Tuple<double, 4u> > tmp_target_w;
	Tuple< Trianglix< Tuple<double, 4u> > > cross_target_w;
	Trianglix< double > tmp;

	//Tuple<double> input; input.setSize(x.ncol());
	Tuple<double, 4u> input_w;
	arma::colvec mu_out = arma::colvec(nb_gene);
	arma::mat target_test(nb_gene,nb_gene);
	Rcpp::NumericMatrix target(nb_gene,nb_gene);
	Rcpp::List crosslist_mean, crosslist;
	Rcpp::Vector<INTSXP> rcpartition;

  Tuple<double, 0u> eigen;
  double tmpei;

  Rcpp::S4 input_out; data.wrRcppdgCMatrix(input_out, nb_gene);
  input_out.slot("Dimnames") = Rcpp::List::create(axes[0],axes[1]);

  if (nbpairs != 0) {
    rcpartition = makeRandomPartition(nbpairs, data.getNBcols(),seed);
    cross_target_w.setSize(nbpairs);
    for(i=0;i<nbpairs;i++) cross_target_w[i].setSize(nb_gene);
    cross_target_w.toZero();
    crosslist = Rcpp::List(nbpairs);
    crosslist_mean = Rcpp::List(nbpairs);
  }
  if (method == 0){
    {
     Tuple<uint32_t> partsize;
     input_w[0] = 1.0f;
     for(col = 0 ; col < data.getNBcols();col++){
        if (auto ite = data.data[col]()) do{
          k = (ite() * (ite() +1)) >> 1;
    	    input_w[1] = *ite;
          if (auto ite2 = data.data[col]()) do{
            if (ite2()<=ite()){
              input_w[2] = *ite2;
              input_w[3] = (*ite2) * (*ite);
              target_w.data[k + ite2()]+= input_w;
              if (nbpairs > 0) cross_target_w[rcpartition[col]].data[k+ite2()] += input_w;
            }
          }while(ite2++);
        }while(ite++);
      }
      if (nbpairs != 0) {
        partsize.setSize(nbpairs).toZero();
        for(i =0; i < data.getNBcols();i++) partsize[rcpartition[col]]++;
      }

      for(l=0;l<nbpairs;l++){
        cross_target_w[l] = target_w - cross_target_w[l];
        partsize[l] = data.getNBcols() - partsize[l];
        for(j =0,k=0; j < nb_gene;j++,k++){
          for(i=0;i<j;i++,k++) target_test.at(i,j) = target_test.at(j,i) = ( cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][2] / partsize[l]) /  partsize[l];
	        target_test.at(j,j) = (cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][1] / cross_target_w[l].data[k][0]) / partsize[l];
	        mu_out.at(j) = cross_target_w[l].data[k][1] / partsize[l];
        }
        crosslist[l] = target_test;
        crosslist_mean[l] = mu_out;
      }
	    for(j =0,k=0; j < nb_gene;j++,k++){
	      for(i=0;i<j;i++,k++) target(i,j) = target(j,i) = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][2] / data.getNBcols()) / data.getNBcols();
	      target(j,j) = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][1] / data.getNBcols()) / data.getNBcols();
	      mu_out.at(j) = target_w.data[k][1] / data.getNBcols();
	    }
    }

  }else{  // method == 1 or 2
     input_w[0] = 1.0f;
     tmp.setSize(nb_gene);
     for(col = 0 ; col < data.getNBcols();col++){
        if (auto ite = data.data[col]()) do{
          k = (ite() * (ite() +1)) >> 1;
    	    input_w[1] = *ite;
          if (auto ite2 = data.data[col]()) do{
            if (ite2()<=ite()){
              input_w[2] = *ite2;
              input_w[3] = (*ite2) * (*ite);
              target_w.data[k + ite2()]+= input_w;
              if (nbpairs > 0) cross_target_w[rcpartition[col]].data[k+ite2()] += input_w;
            }
          }while(ite2++);
        }while(ite++);
      }
      for(l=0;l<nbpairs;l++){
        cross_target_w[l] = target_w - cross_target_w[l];
        for(j =0,k=0; j < nb_gene;j++,k++){
          for(i=0;i<j;i++,k++) {
            tmp.data[k] = ( cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][2] / cross_target_w[l].data[k][0]) /  cross_target_w[l].data[k][0];
            if (!ExOp::isValid(tmp.data[k])) tmp.data[k] =0;
          }
	        tmp.data[k] = (cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][1] / cross_target_w[l].data[k][0]) / cross_target_w[l].data[k][0];
	        mu_out.at(j) = cross_target_w[l].data[k][1] / cross_target_w[l].data[k][0];
        }
        if (method == 2){
            for(k=0;k<100;k++){
                eigen = tmp.getEigenValues();
                j = 0;
                for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
                if  (eigen[j] >0.0) break;
                tmpei = 1.0 / (1.0 - eigen[j] * k);
                if (!ExOp::isValid(tmpei)) {
                   for(j=0;j<eigen.getSize();j++) if (eigen[j] != 0.0) break;
                   if (j == eigen.getSize()) tmp.toOne();
                   else{
                      for(i=j+1;i<eigen.getSize();i++) if ((eigen[i] <= eigen[j])&&(eigen[i] != 0.0)) j =i;
                      for(i=0;i<nb_gene;i++) tmp[i] += eigen[j] * 0.000001;
                   }
                }else {
                  tmpei = 1.0 / (1.0 - eigen[j] * k);
                  uint32_t ii;
                  for(i=0,ii=0;i<nb_gene;i++,ii++){
                      for(j=0;j<i;j++) tmp.data[ii++] *= tmpei;
                  }
                }
            }
            if (k ==100) printf("failed to make the matrix positive definite...\n");
        }
        tmp.wrMatrix(target_test);

        crosslist[l] = target_test;
        crosslist_mean[l] = mu_out;
      }
	    for(j =0,k=0; j < nb_gene;j++,k++){
	      for(i=0;i<j;i++,k++) {
          tmp.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][2] / target_w.data[k][0]) / target_w.data[k][0];
          if (!ExOp::isValid(tmp.data[k])) tmp.data[k] =0;
        }
	      tmp.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][1] / target_w.data[k][0]) / target_w.data[k][0];
	      mu_out.at(j) = target_w.data[k][1] / target_w.data[k][0];
	    }
      if (method == 2){
          for(k=0;k<100;k++){
              eigen = tmp.getEigenValues();
              j = 0;
              for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
              if  (eigen[j] >0.0) break;
              tmpei = 1.0 / (1.0 - eigen[j] * k);
              if (!ExOp::isValid(tmpei)){
                 for(j=0;j<eigen.getSize();j++) if (eigen[j] != 0.0) break;
                 if (j == eigen.getSize()) tmp.toOne();
                 else{
                    for(i=j+1;i<eigen.getSize();i++) if ((eigen[i] <= eigen[j])&&(eigen[i] != 0.0)) j =i;
                    for(i=0;i<nb_gene;i++) tmp[i] += eigen[j] * 0.000001;
                 }
              }else {

                uint32_t ii;
                for(i=0,ii=0;i<nb_gene;i++,ii++){
                    for(j=0;j<i;j++) tmp.data[ii++] *= tmpei;
                }
              }
          }
          if (k ==100) printf("failed to make the matrix positive definite...\n");
      }
      tmp.wrMatrix(target);
  }


  if (nbpairs != 0) {
    return Rcpp::List::create(Rcpp::Named("Covar")= target,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out, Rcpp::Named("Training.Covar") =  crosslist, Rcpp::Named("Training.Mean") =  crosslist_mean, Rcpp::Named("Testset.Partition") = rcpartition, Rcpp::Named("RAND_SEED") = seed);
  }else{
    return Rcpp::List::create(Rcpp::Named("Covar")= target,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out);
  }
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_getFullDeviations)]]
Rcpp::NumericMatrix infernal_getFullDeviations(SEXP InferN0scp, Rcpp::List list){
	Rcpp::XPtr< Infern0scope > scp(InferN0scp);
	std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
	std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
	myHashmap<uint32_t> rows, cols;
	int missing[2];
	if (scp->getGeneCellIndices(rows,cols,genelist,celllist,missing[0], missing[1],true)){
		return R_NilValue; // got nothing to output
	}
	printf("Wmmiissould produce a %i x %i matrix\n",missing[0], missing[1]); fflush(stdout);
	printf("Would produce a %i x %i matrix\n",rows.getSize(),cols.getSize()); fflush(stdout);
	Rcpp::NumericMatrix fout(rows.getSize(),cols.getSize());

	Tuple< double > curproj_R, curproj_M;
	Tuple<double, 3u> coeff;
	uint32_t k;
	double RM[2];
	if (auto itec = cols()) do{
		curproj_R = scp->par_R * scp->cell_state.getColumn(*itec);
		curproj_M = scp->par_M * scp->cell_state.getColumn(*itec);
		coeff = scp->cell_scale[*itec];

		if (auto iter = rows()) do{
			if ((k = scp->data.data[*itec].find(*iter)) != 0xFFFFFFFF) {
				fout(iter.getOffset(),itec.getOffset()) = scp->data.data[*itec].deref(k);
			}else{
				RM[0] = curproj_R.mkInnerProd(	scp->gene_state.getColumn(*iter)) + coeff[0] + scp->gene_scale[*iter][0];
				RM[1] = curproj_M.mkInnerProd(	scp->gene_state.getColumn(*iter)) + coeff[1] + scp->gene_scale[*iter][1];
				fout(iter.getOffset(),itec.getOffset()) = LogPvalue_to_stdnorm(LogPvalue_NBdistrib_exppara_LH_exact(0,RM[0],-RM[1]));
			}
		}while(iter++);
	}while(itec++);
	Rcpp::CharacterVector genenames = Rcpp::CharacterVector(rows.getSize());
	Rcpp::CharacterVector cellnames = Rcpp::CharacterVector(cols.getSize());
	for(k =0; k < rows.getSize();k++) genenames[k] = scp->annot.dicos[INFRN0_AXES_GENE_NAMES].entries[rows.deref(k)];
	for(k =0; k < cols.getSize();k++) cellnames[k] = scp->annot.dicos[INFRN0_AXES_CELL_NAMES].entries[cols.deref(k)];
	Rcpp::colnames(fout) = 	cellnames;
	Rcpp::rownames(fout) = 	genenames;
return(fout);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_cmpModeledVariance)]]
Rcpp::List infernal_cmpModeledVariance(SEXP infscope, Rcpp::List list){ Rcpp::XPtr< Infern0scope > scp(infscope);
	std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
	std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
	uint32_t nbthreads = Rcpp::as< uint32_t >(list[2]);
	uint32_t nbpartition = Rcpp::as< uint32_t >(list[3]);
	uint32_t seed = Rcpp::as< uint32_t >(list[4]);
	uint32_t method = Rcpp::as< uint32_t >(list[5]);
	uint32_t i,j,k,l;

	if ((genelist.size() != 0)&&(genelist[0].length() == 0)) genelist.clear();
	if ((celllist.size() != 0)&&(celllist[0].length() == 0)) celllist.clear();
	myHashmap<uint32_t> rows;
	myHashmap<uint32_t> cols;
	int cnf,gnf;
	if (scp->getGeneCellIndices(rows, cols, genelist, celllist, cnf, gnf,true)){
		return R_NilValue; // got nothing to output
	}
	int pcnf=0;
	for(int i=0;i<cols.getSize();) {
		if (ExOp::isValid(scp->cell_scale[cols.deref(i)])) i++;
		else {pcnf++; cols.erase_from_iterator(i);}
	}

	int pgnf=0;
	for(int i=0;i<rows.getSize();) {
		if (ExOp::isValid(scp->gene_scale[rows.deref(i)])) i++;
		else {pgnf++; rows.erase_from_iterator(i);}
	}
	if ((pcnf != 0)||(pgnf != 0)) printf("warning: %i cells and %i genes are filtered as they have too little expression (only 1-counts is whole row/column).\n", pcnf,pgnf);
	if (rows.getSize() <= 2) {printf("error: only %i genes, needs 3 or more.\n", rows.getSize()); return R_NilValue;}
	if (cols.getSize() <= rows.getSize()) {printf("error: got more genes than cells (%i cells and %i genes).\n", cols.getSize() ,rows.getSize()); return R_NilValue;}
	else if (cols.getSize() <= rows.getSize()*2) printf("warning: got little ammount of cells for the queried number of genes (%i cells and %i genes).\n", cols.getSize() ,rows.getSize());

  	Rcpp::Vector<INTSXP> rcpartition;
	Rcpp::List crosslist_mean, crosslist;
	Tuple< Trianglix< Tuple<double, 4u> > > cross_target_w;
	if (nbpartition > 0){
		rcpartition= makeRandomPartition(nbpartition, cols.getSize(),seed);
		crosslist = Rcpp::List(nbpartition);
		crosslist_mean = Rcpp::List(nbpartition);
    		cross_target_w.setSize(nbpartition);
		for(i=0;i<nbpartition;i++) cross_target_w[i].setSize(rows.getSize());
		cross_target_w.toZero();
	}

	Rcpp::NumericMatrix input_out = Rcpp::NumericMatrix(rows.getSize(), cols.getSize());

	Trianglix< Tuple<double, 4u> > target_w; target_w.setSize(rows.getSize()).toZero();
	Trianglix< Tuple<double, 4u> > tmp_target_w;
	Trianglix< double > tmp;
	Tuple<double, 4u> input_w;
	arma::colvec mu_out = arma::colvec(rows.getSize());
	Rcpp::NumericMatrix target_test; // (rows.getSize(),rows.getSize());
	Rcpp::NumericMatrix target; // (rows.getSize(),rows.getSize());
	Tuple<double, 0u> eigen;
	double tmpei;


	Tuple< double > curproj_R, curproj_M;
	Tuple<double, 3u> coeff;
	Tuple<double, 2u> RM;
	uint32_t daite;

	Tuple<double> curvec; curvec.setSize(rows.getSize());
	if (auto itec = cols()) do{
		curproj_R = scp->par_R * scp->cell_state.getColumn(*itec);
		curproj_M = scp->par_M * scp->cell_state.getColumn(*itec);
		coeff = scp->cell_scale[*itec];
		if (auto iter = rows()) do{
			if ((k = scp->data.data[*itec].find(*iter)) != 0xFFFFFFFF) {
				k = scp->rawdata(*iter,*itec);
			}else k =0;
			RM[0] = curproj_R.mkInnerProd(	scp->gene_state.getColumn(*iter)) + coeff[0] + scp->gene_scale[*iter][0];
			RM[1] = curproj_M.mkInnerProd(	scp->gene_state.getColumn(*iter)) + coeff[1] + scp->gene_scale[*iter][1];
			if (!ExOp::isValid(RM)) {printf("unexpected error"); return R_NilValue;}
			else curvec[iter.getOffset()] = LogPvalue_to_stdnorm(LogPvalue_NBdistrib_exppara_LH_exact(k,RM[0],-RM[1]));
		}while(iter++);
		input_w[0] = 1.0;
		for(i=0,k=0; i < curvec.getSize();i++){
			input_out(i,itec.getOffset()) = curvec[i];
			input_w[1] = curvec[i];
			for(j=0;j <= i;j++){
				input_w[2] = curvec[j];
				input_w[3] = curvec[i] * curvec[j];
				if (nbpartition > 0) cross_target_w[rcpartition[itec.getOffset()]].data[k] += input_w;
				target_w.data[k++] += input_w;
			}
		}
	}while(itec++);

	tmp.setSize(rows.getSize());

	// Rcpp::CharacterVector cellnames;
	Rcpp::CharacterVector genenames = Rcpp::CharacterVector(rows.getSize());
	Rcpp::CharacterVector cellnames = Rcpp::CharacterVector(cols.getSize());
	for(k =0; k < rows.getSize();k++){
		genenames[k] = scp->annot.dicos[INFRN0_AXES_GENE_NAMES].entries[rows.deref(k)];
	}
	for(k =0; k < cols.getSize();k++){
		cellnames[k] = scp->annot.dicos[INFRN0_AXES_CELL_NAMES].entries[cols.deref(k)];
	}

      for(l=0;l<nbpartition;l++){
        cross_target_w[l] = target_w - cross_target_w[l];
        for(j =0,k=0; j < rows.getSize();j++,k++){
          for(i=0;i<j;i++,k++) {
            tmp.data[k] = ( cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][2] / cross_target_w[l].data[k][0]) /  cross_target_w[l].data[k][0];
            if (!ExOp::isValid(tmp.data[k])) tmp.data[k] =0;
          }
	        tmp.data[k] = (cross_target_w[l].data[k][3] - cross_target_w[l].data[k][1]*cross_target_w[l].data[k][1] / cross_target_w[l].data[k][0]) / cross_target_w[l].data[k][0];
	        mu_out.at(j) = cross_target_w[l].data[k][1] / cross_target_w[l].data[k][0];
        }
        if (method != 0){
            for(k=0;k<100;k++){
                eigen = tmp.getEigenValues();
                j = 0;
                for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
                if  (eigen[j] >0.0) break;
                tmpei = 1.0 / (1.0 - eigen[j] * k);
                if (!ExOp::isValid(tmpei)) {
                   for(j=0;j<eigen.getSize();j++) if (eigen[j] != 0.0) break;
                   if (j == eigen.getSize()) tmp.toOne();
                   else{
                      for(i=j+1;i<eigen.getSize();i++) if ((eigen[i] <= eigen[j])&&(eigen[i] != 0.0)) j =i;
                      for(i=0;i<rows.getSize();i++) tmp[i] += eigen[j] * 0.000001;
                   }
                }else {
                  tmpei = 1.0 / (1.0 - eigen[j] * k);
                  uint32_t ii;
                  for(i=0,ii=0;i<rows.getSize();i++,ii++){
                      for(j=0;j<i;j++) tmp.data[ii++] *= tmpei;
                  }
                }
            }
            if (k ==100) printf("failed to make the matrix positive definite...\n");
        }
        tmp.wrMatrix(target_test);
        colnames(target_test) = genenames;
	rownames(target_test) = genenames;
	crosslist[l] = Rcpp::NumericMatrix(target_test);
        crosslist_mean[l] = mu_out;
      }

	    for(j =0,k=0; j < rows.getSize();j++,k++){
	      for(i=0;i<j;i++,k++) {
          tmp.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][2] / target_w.data[k][0]) / target_w.data[k][0];
          if (!ExOp::isValid(tmp.data[k])) tmp.data[k] =0;
        }
	      tmp.data[k] = (target_w.data[k][3] - target_w.data[k][1]*target_w.data[k][1] / target_w.data[k][0]) / target_w.data[k][0];
	      mu_out.at(j) = target_w.data[k][1] / target_w.data[k][0];
	    }
      if (method != 0){
          for(k=0;k<100;k++){
              eigen = tmp.getEigenValues();
              j = 0;
              for(i=1;i<eigen.getSize();i++) if (eigen[i] <= eigen[j]) j =i;
              if  (eigen[j] >0.0) break;
              tmpei = 1.0 / (1.0 - eigen[j] * k);
              if (!ExOp::isValid(tmpei)){
                 for(j=0;j<eigen.getSize();j++) if (eigen[j] != 0.0) break;
                 if (j == eigen.getSize()) tmp.toOne();
                 else{
                    for(i=j+1;i<eigen.getSize();i++) if ((eigen[i] <= eigen[j])&&(eigen[i] != 0.0)) j =i;
                    for(i=0;i<rows.getSize();i++) tmp[i] += eigen[j] * 0.000001;
                 }
              }else {

                uint32_t ii;
                for(i=0,ii=0;i<rows.getSize();i++,ii++){
                    for(j=0;j<i;j++) tmp.data[ii++] *= tmpei;
                }
              }
          }
          if (k ==100) printf("failed to make the matrix positive definite...\n");
      }
      tmp.wrMatrix(target);
      colnames(target) = genenames;
      rownames(target) = genenames;
      colnames(input_out) = cellnames;
	rownames(input_out) = genenames;

//	 target.slot("Dimnames") = Rcpp::List::create(genenames, genenames);
//      input_out.slot("Dimnames") = Rcpp::List::create(genenames, cellnames);

  if (nbpartition != 0) {
    return Rcpp::List::create(Rcpp::Named("Covar")= target,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out, Rcpp::Named("Training.Covar") =  crosslist, Rcpp::Named("Training.Mean") =  crosslist_mean, Rcpp::Named("Testset.Partition") = rcpartition, Rcpp::Named("RAND_SEED") = seed);
  }else{
    return Rcpp::List::create(Rcpp::Named("Covar")= target,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out);
  }


}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_cmpVariance)]]
Rcpp::List infernal_cmpVariance(SEXP infscope, Rcpp::List list){
  std::vector<std::string> genelist = Rcpp::as< std::vector<std::string> >(list[0]);
  std::vector<std::string> celllist = Rcpp::as< std::vector<std::string> >(list[1]);
  uint32_t nbthreads = Rcpp::as< uint32_t >(list[2]);
  uint32_t nbpartition = Rcpp::as< uint32_t >(list[3]);
  uint32_t nbiteration = Rcpp::as< uint32_t >(list[4]);
  double outlier_frequency = Rcpp::as< double >(list[5]);
  uint32_t does_include_uncertainty = Rcpp::as< uint32_t >(list[6]);
  uint32_t seed = Rcpp::as< uint32_t >(list[7]);


  uint32_t inputClass=0;
  if (Rf_isMatrix(infscope)) inputClass=1;
  else if (Rf_isS4(infscope)) inputClass=2;

  SparseMatrix<double> data;

  ThreadBase tb;
	tb.toSize(nbthreads-1);
  Rcpp::CharacterVector axes[2];
  uint32_t nb_gene;

  switch(inputClass){
    case 0:{
      Rcpp::XPtr< Infern0scope > scptr(infscope);
	scptr->wrSubMatrix(data, genelist, celllist, axes, 3);
      nb_gene = axes[0].size();
	printf("%i %i\n", (int)axes[0].length(), (int)axes[1].length());
    }break;
    case 1:{
    Rcpp::NumericMatrix dadata(infscope);
    data.rdRcpp(dadata);
    nb_gene = data.computeNBrows();
    axes[0] = Rcpp::rownames(dadata);
    axes[1] = Rcpp::colnames(dadata);
    }break;
    case 2:{
    Rcpp::S4 dadata(infscope);
    data.rdRcppdgCMatrix(dadata);
    nb_gene = data.computeNBrows();
    Rcpp::List davec = dadata.slot("Dimnames");
    axes[0] = davec[0];
    axes[1] = davec[1];
    }break;
  }
  Rcpp::Vector<INTSXP> rcpartition;
  Rcpp::List crosslist_mean, crosslist;
  Vector<double> daweight;

  if (nbpartition > 0){
    daweight.setSize(data.getNBcols());
    rcpartition= makeRandomPartition(nbpartition, data.getNBcols(),seed);
    crosslist = Rcpp::List(nbpartition);
    crosslist_mean = Rcpp::List(nbpartition);
  }

  GaussElem<Tuple<double> > toinf;
  SparseMatrix<double> partdata;
  uint32_t i,j,k;
  arma::colvec mu_out = arma::colvec(nb_gene);
  Rcpp::NumericMatrix cov;
  arma::mat training;


  for(i=0;i< nbpartition;i++){
    for(j=0;j<data.getNBcols();j++) daweight[j] = (rcpartition[j] == i) ? 0.0 : 1.0;
	  toinf.runTaskPartialObservationMaximize(tb, data, nb_gene,nbiteration,outlier_frequency,NULL,does_include_uncertainty != 0, &daweight);
    toinf.cov.wrMatrix(training);
    crosslist[i] = training;
    for(j=0;j<nb_gene;j++) mu_out.at(j) = toinf.mean[j];
    crosslist_mean[i] = mu_out;
  }
	toinf.runTaskPartialObservationMaximize(tb, data, nb_gene,nbiteration,outlier_frequency,NULL,does_include_uncertainty != 0,NULL);

  toinf.cov.wrMatrix(cov);
  for(i=0;i<nb_gene;i++) mu_out.at(i) = toinf.mean[i];

  colnames(cov) = Rcpp::as< Rcpp::CharacterVector >(axes[0]);
  rownames(cov) = Rcpp::as< Rcpp::CharacterVector >(axes[0]);


  Rcpp::S4 input_out; data.wrRcppdgCMatrix(input_out, nb_gene);
  input_out.slot("Dimnames") = Rcpp::List::create(axes[0],axes[1]);
  if (nbpartition > 0) {
    return Rcpp::List::create(Rcpp::Named("Covar") = cov,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out, Rcpp::Named("Training.Covar") =  crosslist, Rcpp::Named("Training.Mean") =  crosslist_mean, Rcpp::Named("Testset.Partition") = rcpartition);
  }else{
    return Rcpp::List::create(Rcpp::Named("Covar") = cov,Rcpp::Named("Mean")= mu_out,Rcpp::Named("Input")= input_out);
  }
}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.internal_Infern0_hierarchicalClustering_CustomGroups)]]
Rcpp::List internal_Infern0_hierarchicalClustering(SEXP data, SEXP cellgroup, int nbclusters){
  Rcpp::XPtr< LFHPrimitive::SparseMatrix<double> > matptr(data);
  std::vector<uint32_t> groupId = Rcpp::as<std::vector<uint32_t> >(cellgroup);

  Vector< Tuple< WeightElem<double, 2> , 0u > > to_clust;
  uint32_t i,j;
  Vector<uint32_t> rowsize = matptr->compute_rowsizes();
  to_clust.setSize(rowsize.getSize());
  for(i=0;i<rowsize.getSize();i++) to_clust[i].setSize(nbclusters).toZero();
  for(j=0;j<matptr->getNBcols();j++){
      for(i=0;i<matptr->data.getSize();i++){
         to_clust[matptr->data[j].deref_key(i)][groupId[j]] = WeightElem<double, 2>(matptr->data[j].deref(i));
      }
  }
  for(i=0;i<rowsize.getSize();i++){
    for(j=0;j<nbclusters;j++) to_clust[i][j].replaceNan(0.0f,0.01f,true);
  }
  Rcpp::XPtr<Forest<double,2> > da_clust(new Forest<double,2>);
  //da_clust.cluster_likelihood_ratio(to_clust, batta_metric);

return Rcpp::List::create(Rcpp::Named("VectorPTR")= wrap(da_clust));}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_mergeSparseMatrices)]]
Rcpp::S4 infernal_mergeSparseMatrices(Rcpp::S4 matrixA, Rcpp::S4 matrixB){

  Rcpp::List axeA = matrixA.slot("Dimnames");

  Rcpp::List axeB = matrixB.slot("Dimnames");

  myHashmap<string, uint32_t> rowmap;
  myHashmap<string, uint32_t> colmap;
  myHashmap<uint32_t, uint32_t> rowmap2;
  myHashmap<uint32_t, uint32_t> colmap2;


  Rcpp::CharacterVector curcv = axeA[0];
  for(uint32_t i=0;i<curcv.size();i++) rowmap[string(curcv[i])] = i;
  curcv = axeA[1];
  for(uint32_t i=0;i<curcv.size();i++) colmap[string(curcv[i])] = i;
  uint32_t j;
  curcv = axeB[0];

  for(uint32_t i=0;i<curcv.size();i++) {
    if ((j = rowmap.find(string(curcv[i]))) == 0xFFFFFFFF) {
      j = rowmap.getSize(); rowmap2[i] = j;  rowmap[string(curcv[i])]= j;
    }else rowmap2[i] = rowmap.deref(j);
  }

  curcv = axeB[1];
  for(uint32_t i=0;i<curcv.size();i++) {
    if ((j = colmap.find(string(curcv[i]))) == 0xFFFFFFFF) {
    j = colmap.getSize(); colmap2[i] = j; colmap[string(curcv[i])]= j;
    }else colmap2[i] = colmap.deref(j);
  }

  SparseMatrix<double> merged;
  merged.rdRcppdgCMatrix(matrixA);
  while(merged.data.getSize() < colmap.getSize()) merged.data.push_back();

  Rcpp::CharacterVector rownames,colnames;
  for(uint32_t i=0;i<rowmap.getSize();i++) rownames.push_back(rowmap.deref_key(i).c_str());
  for(uint32_t i=0;i<colmap.getSize();i++) colnames.push_back(colmap.deref_key(i).c_str());

  const int RTYPE= Rcpp::traits::r_sexptype_traits<double>::rtype;
  Rcpp::IntegerVector Bi = matrixB.slot("i");
  Rcpp::IntegerVector Bp = matrixB.slot("p");
  Rcpp::Vector<RTYPE> Bx = matrixB.slot("x");

  uint32_t pite =1;
  for(j=0;j<Bi.size();j++){
    while((pite < Bp.size())&&(Bp[pite] == j)) pite++;
    merged.data[colmap2[pite-1]][rowmap2[Bi[j]]] = Bx[j];
  }


  Rcpp::S4 fout;
  merged.wrRcppdgCMatrix(fout, rowmap.getSize());
  fout.slot("Dimnames") = Rcpp::List::create(rownames,colnames);
return fout;}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_genSynthetic)]]
Rcpp::List infernal_genSynthetic(SEXP infscope, Rcpp::List data){
    Rcpp::XPtr< Infern0scope > scptr(infscope);

    int32_t nbhiddenrows = Rcpp::as<int32_t>(data[0]); // ignored for now...
    int32_t nbhiddencols = Rcpp::as<int32_t>(data[1]); // ignored for now...

    Tuple<uint32_t, 2u> coor;
    coor[0] = Rcpp::as<int32_t>(data[2]);
    coor[1] = Rcpp::as<int32_t>(data[3]);
    double mixed_weight = Rcpp::as<double>(data[4]);

    TMatrix< double > synt_R; synt_R.setSizes(coor).toRand();
    TMatrix< double > synt_M; synt_M.setSizes(coor).toRand();
    TMatrix< double > synt_D; synt_D.setSizes(coor).toRand();
    synt_M -= 0.75;
    synt_M *= 4.0;
    if (auto ite = synt_R.getIterator()) do{
      (*ite) = log(0.0001 + 4.0 * (*ite) );
    }while(ite++);


    coor[1] = coor[0] = 2u;
    Tuple< Tuple<double ,2u> > par_colscale;
    Vector< Tuple<double ,2u> > randgroup; randgroup.setSize(40);

    TMatrix<double> colhid;
    coor[0] = synt_R.sizes[1];
    coor[1] = coor[0] * randgroup.getSize(); colhid.setSizes(coor).toRand();

    for(coor[0] = 0; coor[0] < synt_R.sizes[1];coor[0]++){
      randgroup.toRand();
      randgroup.sort();
      for(coor[1]=0;coor[1] < randgroup.getSize();coor[1]++) par_colscale.push_back(randgroup[coor[1]]);
    }
    for(coor[0] = 0; coor[0] < par_colscale.getSize();coor[0]++){
      par_colscale[coor[0]][0] =3.0 * par_colscale[coor[0]][0] - 1.5;
      par_colscale[coor[0]][1] =1.0 * par_colscale[coor[0]][1] - 0.5;
    }

    double sum;
    for(coor[1]=0;coor[1]< colhid.sizes[1];coor[1]++) {
        coor[0] = (synt_R.sizes[1] * coor[1]) / colhid.sizes[1];
        colhid(coor) += (1.0 / mixed_weight) * synt_R.sizes[1];
        sum=0.0;
        for(coor[0]=0; coor[0] <colhid.sizes[0];coor[0]++) sum+= colhid(coor);
        for(coor[0]=0; coor[0] <colhid.sizes[0];coor[0]++) colhid(coor) /= sum;
    }
    int tmp;
    TMatrix<double> rowhid;
    Tuple< Tuple<double ,2u> > par_rowscale;
    coor[0] = synt_R.sizes[0];
    coor[1] = coor[0] * randgroup.getSize();
    rowhid.setSizes(coor).toRand();
    //par_rowscale.toSize(coor[1]).toRand();

    for(coor[0] = 0; coor[0] < synt_R.sizes[0];coor[0]++){
      randgroup.toRand();
      randgroup.sort();
      for(coor[1]=0;coor[1] < randgroup.getSize();coor[1]++) par_rowscale.push_back(randgroup[coor[1]]);
    }
    for(coor[0] = 0; coor[0] < par_rowscale.getSize();coor[0]++){
      par_rowscale[coor[0]][0] =3.0 * par_rowscale[coor[0]][0] - 1.5;
      par_rowscale[coor[0]][1] =1.0 * par_rowscale[coor[0]][1] - 0.5;
    }

    for(coor[1]=0;coor[1]< rowhid.sizes[1];coor[1]++) {
        coor[0] = (synt_R.sizes[0] * coor[1]) / rowhid.sizes[1];
        rowhid(coor) += (1.0 / mixed_weight) * synt_R.sizes[0];
        sum=0.0;
        for(coor[0]=0; coor[0] <rowhid.sizes[0];coor[0]++) sum+= rowhid(coor);
        for(coor[0]=0; coor[0] <rowhid.sizes[0];coor[0]++) rowhid(coor) /= sum;
    }



    scptr->rawdata.setNBcols(par_colscale.getSize());
    Tuple<double, 0u> proj_R,proj_M,proj_D;
    double val_r,val_m,val_d;
    double bufm[2];
    WeightElem<double, 2u> dropdrop; dropdrop.toZero();
    for(coor[1]=0;coor[1]< par_colscale.getSize();coor[1]++) {
        proj_R = synt_R * colhid.getColumn(coor[1]);
        proj_M = synt_M * colhid.getColumn(coor[1]);
        proj_D = synt_D * colhid.getColumn(coor[1]);
        for(coor[0]=0;coor[0]< par_rowscale.getSize();coor[0]++) {
            val_r = proj_R.mkInnerProd(rowhid.getColumn(coor[0]));
            val_r += par_rowscale[coor[0]][1];
            val_r += par_colscale[coor[1]][1];
            val_m = proj_M.mkInnerProd(rowhid.getColumn(coor[0]));
            val_m += par_rowscale[coor[0]][0];
            val_m += par_colscale[coor[1]][0];
            val_d = proj_D.mkInnerProd(rowhid.getColumn(coor[0]));
            dropdrop += WeightElem<double, 2u>(val_d);
            //if (rand() < (val_d * RAND_MAX)) scptr->rawdata(coor) = tmp = (int)(exp(10.0*val_m));
            if (rand() < (val_d * RAND_MAX)) {
              deflatedNB_getMeanAndVar(bufm, val_r,val_m);
              scptr->rawdata(coor) = (bufm[0] < 1000) ? deflatedNB_sample(val_r,val_m,2000) : (uint32_t)(bufm[0] + sampleGaussian() * sqrt(bufm[1]));
              if (scptr->rawdata(coor) > 100000) scptr->rawdata(coor) = (uint32_t) bufm[0];
            }
        }
    }
  Rcpp::NumericMatrix nummatrix;
  scptr->rawdata.wrRcpp(nummatrix);

  char buffer[256];
  std::vector<std::string> names;
  for(coor[0]=0;coor[0]< par_rowscale.getSize();coor[0]++) {sprintf(buffer, "gene%i", coor[0]); names.push_back(std::string(buffer));}
  Rcpp::CharacterVector genenames = Rcpp::wrap(names);
  scptr->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_GENE_NAMES);
  names.clear();
  for(coor[0]=0;coor[0]< par_colscale.getSize();coor[0]++) {sprintf(buffer, "cell%i", coor[0]); names.push_back(std::string(buffer));}
  genenames = Rcpp::wrap(names);
  scptr->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_CELL_NAMES);
  names.clear();

  scptr->annot.wrAxesNamesAl(nummatrix,  (int)INFRN0_AXES_GENE_NAMES,  (int)INFRN0_AXES_CELL_NAMES);


return Rcpp::List::create(Rcpp::Named("data")= nummatrix,Rcpp::Named("axename")= genenames);}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_initModelWithClustering)]]
void infernal_initModelWithClustering(SEXP infscope, Rcpp::List data){ Rcpp::XPtr< Infern0scope > scp(infscope);
	int32_t nbhiddenrows = Rcpp::as<int32_t>(data[0]);
	int32_t nbhiddencols = Rcpp::as<int32_t>(data[1]);
	vector<uint32_t> cluster = Rcpp::as<vector<uint32_t> >(data[2]);

	scp->par_D.setSizes(nbhiddenrows, nbhiddencols);
	scp->par_R.setSizes(nbhiddenrows, nbhiddencols);
	scp->par_M.setSizes(nbhiddenrows, nbhiddencols);
	scp->cell_state.data.setSize(scp->getNBcells());
	for(int i=0;i<scp->rawdata.getNBcols();i++) scp->cell_state.data[i][cluster[i]] = 1.0;

	MPReachDistribution sampler;
	sampler.setDefaultRatesAndPrior(nbhiddenrows, 2.0);
	scp->gene_state.data = std::move(sampler.mkSamples(scp->getNBgenes()));

}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_modelhidden)]]
Rcpp::List infernal_modelhidden(SEXP infscope, Rcpp::List data){
  Rcpp::XPtr< Infern0scope > scptr(infscope);

  int32_t nbhiddenrows = Rcpp::as<int32_t>(data[0]);
  int32_t nbhiddencols = Rcpp::as<int32_t>(data[1]);
  int32_t nb_step = Rcpp::as<int32_t>(data[2]);
  int32_t nb_thread = Rcpp::as<int32_t>(data[3]);
  bool outputtiff = Rcpp::as<bool>(data[4]);

  BiClassifier<uint32_t> dafun;

  ThreadBase tb; tb.toSize(nb_thread-1);
  printf("Randomize\n");fflush(stdout);

  Vector<uint32_t> rowsize = scptr->rawdata.compute_rowsizes();

  dafun.setRandom(scptr->annot.dicos[(int)INFRN0_AXES_GENE_NAMES].getSize(), scptr->annot.dicos[(int)INFRN0_AXES_CELL_NAMES].getSize(), nbhiddenrows,nbhiddencols);

  if (scptr->cell_scale.getSize() == scptr->rawdata.getNBcols()){
    printf("reusing scales\n");
    dafun.par_rowscale = scptr->gene_scale;
    dafun.par_colscale = scptr->cell_scale;
    if (scptr->cell_state.getNBcols() == scptr->rawdata.getNBcols()){
       printf("reusing states\n");
       dafun.rowhid = scptr->gene_state;
       dafun.colhid = scptr->cell_state;
    }
  }else{
    dafun.initScale(scptr->rawdata);
  }

  if (scptr->par_D.data.getSize() != 0){
    if ((scptr->par_D.sizes[0] <= dafun.par_D.sizes[0])&&(scptr->par_D.sizes[1] <= dafun.par_D.sizes[1])){
      printf("reusing matrices\n");
      if (auto ite = scptr->par_D.getIterator()) do{
         dafun.par_D(ite()) = *ite;
         dafun.par_RM(ite())[0] = scptr->par_R(ite());
         dafun.par_RM(ite())[1] = scptr->par_M(ite());
      }while(ite++);
    }else if ((scptr->par_D.sizes[0] >= dafun.par_D.sizes[0])&&(scptr->par_D.sizes[1] >= dafun.par_D.sizes[1])){
      printf("reusing parts of matrices\n");
      if (auto ite = scptr->par_D.getIterator()) do{
         *ite = dafun.par_D(ite());
         dafun.par_RM(ite())[0] = scptr->par_R(ite());
         dafun.par_RM(ite())[1] = scptr->par_M(ite());
      }while(ite++);
    }
  }

  Vector<uint32_t> excl;
  printf("Running EM\n");fflush(stdout);
  dafun.run2D_EM_v2(scptr->rawdata, excl, tb, nb_step, nb_thread);
  printf("Finished Successfully\n");fflush(stdout);

  if (outputtiff){
    printf("Saving\n");fflush(stdout);
    dafun.exportTiff("/nfs/users/nfs_l/lh20/");
    printf("Saving2\n");fflush(stdout);
    dafun.exportZscore("/nfs/users/nfs_l/lh20/", scptr->rawdata);
  }

  scptr->data = dafun.normalize(scptr->rawdata);

  printf("Finished Normalization Successfully\n");fflush(stdout);

   bool hasresized;
   if (hasresized = (scptr->gene_scale.getSize() != dafun.par_rowscale.getSize())) scptr->gene_scale.setSize(dafun.par_rowscale.getSize());
   	for(int i=0;i< dafun.par_rowscale.getSize();i++){
		scptr->gene_scale[i][0] = dafun.par_rowscale[i][0];
		scptr->gene_scale[i][1] = dafun.par_rowscale[i][1];
		if (hasresized) scptr->gene_scale[i][2] = 0.0;
	}
   if (hasresized = (scptr->cell_scale.getSize() != dafun.par_colscale.getSize())) scptr->cell_scale.setSize(dafun.par_colscale.getSize());
   for(int i=0;i< dafun.par_colscale.getSize();i++){
	scptr->cell_scale[i][0] = dafun.par_colscale[i][0];
	scptr->cell_scale[i][1] = dafun.par_colscale[i][1];
	if (hasresized) scptr->cell_scale[i][2] = 0.0;
	}
  scptr->gene_state.toMemmove(dafun.rowhid);
  scptr->cell_state.toMemmove(dafun.colhid);
  scptr->par_D.toMemmove(dafun.par_D);
  scptr->par_R.setSizes(dafun.par_RM.sizes[0],dafun.par_RM.sizes[1]);
  scptr->par_M.setSizes(dafun.par_RM.sizes[0],dafun.par_RM.sizes[1]);
  if (auto ite = dafun.par_RM.getIterator()) do{
    scptr->par_R(ite()) = (*ite)[0]; scptr->par_M(ite()) = (*ite)[1];
  }while(ite++);

  printf("Populated Scope Successfully\n");fflush(stdout);

  Rcpp::S4 normalized,hiddenrows, hiddencols;
  scptr->data.wrRcppdgCMatrix(normalized);
  dafun.rowhid.wrRcppdgCMatrix(hiddenrows);
  dafun.colhid.wrRcppdgCMatrix(hiddencols);

  char buffer[256];
  uint32_t i;
  std::vector<std::string> names;
  for(i=0;i< dafun.par_RM.sizes[0];i++) {sprintf(buffer, "gene state %i", i); names.push_back(std::string(buffer));}
  Rcpp::CharacterVector genenames = Rcpp::wrap(names);
  scptr->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_GENE_STATE_NAMES);
  names.clear();
  for(i=0;i< dafun.par_RM.sizes[1];i++) {sprintf(buffer, "cell state %i", i); names.push_back(std::string(buffer));}
  genenames = Rcpp::wrap(names);
  scptr->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_CELL_STATE_NAMES);
  names.clear();

  return Rcpp::List::create(Rcpp::Named("data")= Rcpp::wrap(normalized),Rcpp::Named("hiddenrows")= Rcpp::wrap(hiddenrows), Rcpp::Named("hiddencols")= Rcpp::wrap(hiddencols));
  //return Rcpp::List::create();
}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_modelhidden_v2)]]
Rcpp::List infernal_modelhidden_v2(SEXP infscope, Rcpp::List data){
  Rcpp::XPtr< Infern0scope > scp(infscope);

  int32_t nbhiddenrows = Rcpp::as<int32_t>(data[0]);
  int32_t nbhiddencols = Rcpp::as<int32_t>(data[1]);
  int32_t nb_step = Rcpp::as<int32_t>(data[2]);
  int32_t nb_thread = Rcpp::as<int32_t>(data[3]);
  int32_t flag = Rcpp::as<bool>(data[4]);
  vector<uint32_t> batch = Rcpp::as< vector<uint32_t> >(data[5]);

  BiClassifier<uint32_t> dafun;

  Vector<uint32_t> rowsize = scp->rawdata.compute_rowsizes();
  if (dafun.initHidden(nbhiddenrows,nbhiddencols,scp->getNBgenes(),scp->getNBcells(), 2.0f,2.0f) != 0 ){
	return Rcpp::List::create();
  }

  bool isIncohent=false;
  if (scp->cell_scale.getSize() == scp->rawdata.getNBcols()) dafun.par_colscale = scp->cell_scale;
  else dafun.initScale(scp->rawdata,false,true);
  if (scp->cell_scale.getSize() == scp->rawdata.getNBcols()) dafun.par_rowscale = scp->gene_scale;
  else dafun.initScale(scp->rawdata,true,false);

  if (scp->par_D.data.getSize() != 0){
    if ((scp->par_D.sizes[0] <= dafun.par_D.sizes[0])&&(scp->par_D.sizes[1] <= dafun.par_D.sizes[1])){
      if (scp->cell_state.getNBcols() == scp->rawdata.getNBcols()) {dafun.colhid = scp->cell_state; dafun.rowhid = scp->gene_state;}
      printf("reusing states\n");
      if (auto ite = scp->par_D.getIterator()) do{
         dafun.par_D(ite()) = *ite;
         dafun.par_RM(ite())[0] = scp->par_R(ite());
         dafun.par_RM(ite())[1] = scp->par_M(ite());
      }while(ite++);
    } else if ((scp->par_D.sizes[0] >= dafun.par_D.sizes[0])&&(scp->par_D.sizes[1] >= dafun.par_D.sizes[1])){
      printf("reusing parts of matrices\n");
      if (auto ite = scp->par_D.getIterator()) do{
         *ite = dafun.par_D(ite());
         dafun.par_RM(ite())[0] = scp->par_R(ite());
         dafun.par_RM(ite())[1] = scp->par_M(ite());
      }while(ite++);
    }
  }

  Vector<uint32_t> excl;
  foreach.startThreadArray(nb_thread);

  if (batch.size() == scp->cell_scale.getSize()){
  	Vector<uint32_t> batch_ids(batch);
	dafun.run2D_EM_v5(scp->rawdata, excl, batch_ids, nb_step, (flag & 1) != 0);
  }else{
  	dafun.run2D_EM_v4(scp->rawdata, excl, nb_step, (flag & 1) != 0);
  }
  printf("Finished Successfully\n");fflush(stdout);
/*  if (outputtiff){
    printf("Saving\n");fflush(stdout);
    dafun.exportTiff("/nfs/users/nfs_l/lh20/");
    printf("Saving2\n");fflush(stdout);
    dafun.exportZscore("/nfs/users/nfs_l/lh20/", scp->rawdata);
  }*/

  scp->data = dafun.normalizeMK2(scp->rawdata); // tmptmp

  printf("Finished Normalization Successfully\n");fflush(stdout);

  foreach.stopThreadArray();

  printf("stopped the threads");fflush(stdout);



   bool hasresized;
   if (hasresized = (scp->gene_scale.getSize() != dafun.par_rowscale.getSize())) scp->gene_scale.setSize(dafun.par_rowscale.getSize());
   	for(int i=0;i< dafun.par_rowscale.getSize();i++){
		scp->gene_scale[i][0] = dafun.par_rowscale[i][0];
		scp->gene_scale[i][1] = dafun.par_rowscale[i][1];
		if (hasresized) scp->gene_scale[i][2] = 0.0;
	}
   if (hasresized = (scp->cell_scale.getSize() != dafun.par_colscale.getSize())) scp->cell_scale.setSize(dafun.par_colscale.getSize());
   for(int i=0;i< dafun.par_colscale.getSize();i++){
	scp->cell_scale[i][0] = dafun.par_colscale[i][0];
	scp->cell_scale[i][1] = dafun.par_colscale[i][1];
	if (hasresized) scp->cell_scale[i][2] = 0.0;
	}
  scp->gene_state.toMemmove(dafun.rowhid);
  scp->cell_state.toMemmove(dafun.colhid);
  scp->par_D.toMemmove(dafun.par_D);
  scp->par_R.setSizes(dafun.par_RM.sizes[0],dafun.par_RM.sizes[1]);
  scp->par_M.setSizes(dafun.par_RM.sizes[0],dafun.par_RM.sizes[1]);
  if (auto ite = dafun.par_RM.getIterator()) do{
    scp->par_R(ite()) = (*ite)[0]; scp->par_M(ite()) = (*ite)[1];
  }while(ite++);

  printf("Populated Scope Successfully\n");fflush(stdout);

  Rcpp::S4 normalized,hiddenrows, hiddencols;
  scp->data.wrRcppdgCMatrix(normalized);
  dafun.rowhid.wrRcppdgCMatrix(hiddenrows);
  dafun.colhid.wrRcppdgCMatrix(hiddencols);

  char buffer[256];
  uint32_t i;
  std::vector<std::string> names;
  for(i=0;i< dafun.par_RM.sizes[0];i++) {sprintf(buffer, "gene state %i", i); names.push_back(std::string(buffer));}
  Rcpp::CharacterVector genenames = Rcpp::wrap(names);
  scp->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_GENE_STATE_NAMES);
  names.clear();
  for(i=0;i< dafun.par_RM.sizes[1];i++) {sprintf(buffer, "cell state %i", i); names.push_back(std::string(buffer));}
  genenames = Rcpp::wrap(names);
  scp->annot.setAxeNamesAl(genenames,  (uint32_t)INFRN0_AXES_CELL_STATE_NAMES);
  names.clear();

  return Rcpp::List::create(Rcpp::Named("data")= Rcpp::wrap(normalized),Rcpp::Named("hiddenrows")= Rcpp::wrap(hiddenrows), Rcpp::Named("hiddencols")= Rcpp::wrap(hiddencols));
}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernal_exportpcs)]]
Rcpp::List infernal_exportpcs(SEXP infscope, Rcpp::List rlist){
	Rcpp::XPtr< Infern0scope > scptr(infscope);
	vector<int32_t> reduceddims = Rcpp::as<vector<int32_t> >(rlist[0]);
	int32_t nb_threads = Rcpp::as<int32_t>(rlist[1]);
	int32_t nb_step = Rcpp::as<int32_t>(rlist[2]);
	int32_t flags = Rcpp::as<int32_t>(rlist[3]);

	Rcpp::S4 normalized;
	// var = w * sum x x^t
  	// x_ij ~ nb()
	/*
	GaussElem<Tuple<double> > data;	data.setSize(reduceddims.size());
	data.toZero();
	Tuple<double> data_input; data_input.setSize(reduceddims.size());

	uint32_t i,j,ind;
	for(i=0;i< scptr->data.getNBcols();i++){
		for(j=0;j<reduceddims.size();j++){
			if ((ind = scptr->data.data[i].find(reduceddims[j])) != 0xFFFFFFFF) data_input[j] =  scptr->data.data[i].deref(ind);
			else{

			}
		}
		data += GaussElem<Tuple<double> >(data_input);
	}*/

 	 BiClassifier<uint32_t> dascp;


	 if (scptr->cell_scale.getSize() == scptr->rawdata.getNBcols()){
	    printf("reusing scales\n");
	    dascp.par_rowscale = scptr->gene_scale;
	    dascp.par_colscale = scptr->cell_scale;
	    if (scptr->cell_state.getNBcols() == scptr->rawdata.getNBcols()){
      		 printf("reusing states\n");
		 dascp.rowhid = scptr->gene_state;
	      	 dascp.colhid = scptr->cell_state;
	    }
		 dascp.par_RM.setSizes(scptr->par_R.sizes[0],scptr->par_R.sizes[1]);
		 if (auto ite = dascp.par_RM.getIterator()) do{
		    (*ite)[0] = scptr->par_R(ite());
		    (*ite)[1] = scptr->par_M(ite());
		 }while(ite++);

	  }else{
	    dascp.initScale(scptr->rawdata);
	  }


  	ThreadBase tb; tb.toSize(nb_threads-1);

	Tuple< double > fout_rowdrop, fout_coldrop;
	dascp.run2D_EM_dropout(tb, fout_rowdrop, fout_coldrop, scptr->rawdata, nb_step);

	if (flags & 1){ // renormalize with dropout!
		scptr->data = dascp.normalize_dropout(scptr->rawdata);
	}

	bool hasresized;
	if (hasresized = (scptr->gene_scale.getSize() != fout_rowdrop.getSize())) scptr->gene_scale.setSize(fout_rowdrop.getSize());
   	for(int i=0;i< fout_rowdrop.getSize();i++){
		scptr->gene_scale[i][2] = fout_rowdrop[i];
		if (hasresized) scptr->gene_scale[i][0] =  scptr->gene_scale[i][1] = 0.0;
	}
	if (hasresized = (scptr->cell_scale.getSize() != fout_coldrop.getSize())) scptr->cell_scale.setSize(fout_coldrop.getSize());
   	for(int i=0;i< fout_coldrop.getSize();i++){
		scptr->cell_scale[i][2] = fout_coldrop[i];
		if (hasresized) scptr->cell_scale[i][0] = scptr->cell_scale[i][1] = 0.0;
	}



	return(Rcpp::List::create());
}





// library(InferN0); scptr <- Infern0GenerateSyntheticData(0,0,8,5); fakeout <- Infern0ModelData(scptr,8,5,25,1); Infern0_HierarchicalClustering(scptr) ; Infern0_HierarchicalClustering(scptr, cluster.cells=FALSE); Infern0saveHierarchical(scptr, "/lustre/scratch117/cellgen/team218/lh20/test");

// library(InferN0); scptr <- Infern0GenerateSyntheticData(0,0,8,5); fakeout <- Infern0ModelData(scptr,8,5,25,1); testout <- Infern0_IdentifyNetwork(scptr, c("gene10","gene12","gene13","gene18","gene20","gene90","gene100","gene120","gene130","gene180","gene200","gene220","gene230","gene280"))

// library(InferN0); scptr <- Infern0GenerateSyntheticData(0,0,8,5);
// fakeout <- Infern0ModelData(scptr,8,5,25,1) ; saveInferN0file(dascp,"/lustre/scratch117/cellgen/team218/lh20/infn0.scp");
// heptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/infn0.scp") ;
// Infern0ModelData(heptr,8,5,25,1) ;

// sro <- readRDS("neuropilot.sro.rds")
// initInferN0(sro@raw.data)

// library(InferN0); sro<- readRDS("/lustre/scratch117/cellgen/team218/lh20/neuropilot2.sro.rds"); scptr <- initInferN0(sro@raw.data); fakeout <- Infern0ModelData(scptr,8,5,100,1); testout <- Infern0_IdentifyNetwork(scptr, c("GRIA1","CALB2","ROBO1","ROBO2","GAD2","DLX5","CNTNAP2","SIX3","MEIS2","MEF2C","NPAS3","GRIK2","CCND2","SOX4"))
//c("Gria1","Calb2","Robo1","Robo2","Gad2","Dlx5","Cntnap2","Six3","Meis2","Mef2c","Npas3","Grik2","Ccnd2","Sox4")

// raw<- Infern0get(scptr,"raw.data"); data <- Infern0get(scptr,"data"); hid<- Infern0get(scptr,"cell.state");
// library(Seurat); sro <- CreateSeuratObject(raw.data = rbind(hid,raw) ,min.cell = 0, min.genes =0); sro@scale.data <- rbind(hid,data)

// library(InferN0); scptr <- initInferN0(readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/BrainNeurons_FACS.rds"))
// library(InferN0); scptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/inferN0output/zeisel.ifn.scp"); testout <- Infern0_IdentifyNetwork(scptr,c("Gria1","Calb2","Robo1","Robo2","Gad2","Dlx5","Cntnap2","Six3","Meis2","Mef2c","Npas3","Grik2","Ccnd2","Sox4"))

// library(InferN0); scptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/inferN0output/tabulamuris_brainneuronsF.ifn.scp"); testout <- Infern0_IdentifyNetwork(scptr,c("Gria1","Calb2","Robo1","Robo2","Gad2","Dlx5","Cntnap2","Six3","Meis2","Mef2c","Npas3","Grik2","Ccnd2","Sox4"))


// detach("package:InferN0", unload=T); library(InferN0);
// mmm <- readRDS("/lustre/scratch117/cellgen/team218/lh20/lung3.sre.rds")
// library(InferN0); lungptr <- initInferN0(mmm@raw.data)
// saveInferN0file(lungptr,"/lustre/scratch117/cellgen/team218/lh20/lung3.scp")
// mmm <- Infern0get(lungptr,"raw.data")


// library(InferN0); lungptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/lung3.scp")
// fakeout <- Infern0ModelData(lungptr,12,8,25,12);

// library(InferN0); scptr <- Infern0GenerateSyntheticData(0,0,8,5); fakeout <- Infern0ModelData(scptr,8,5,25,1);
// Infern0performHierarchical(scptr); Infern0saveHierarchical(scptr, "/lustre/scratch117/cellgen/team218/lh20/test")


// Infern0performHierarchical(lungptr); Infern0saveHierarchical(lungptr, "/lustre/scratch117/cellgen/team218/lh20/lung");
//

// DoHeatmap(object = sro3,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
// 2635190




// sh ~/scRNASeq/DO_Rscript.sh ~/scRNASeq/runInferno.Rscript 20,10,100 /lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/BrainNeurons_FACS.rds /lustre/scratch117/cellgen/team218/lh20/inferN0output/tabulamuris_brainneuronsF.ifn.scp

// sh ~/scRNASeq/DO_Rscript.sh ~/scRNASeq/runInferno.Rscript 20,10,100 /lustre/scratch117/cellgen/team218/lh20/lung3.sre.rds /lustre/scratch117/cellgen/team218/lh20/inferN0output/lung3.ifn.scp

// sh ~/scRNASeq/DO_Rscript.sh ~/scRNASeq/runInferno.Rscript 20,10,100 /lustre/scratch117/cellgen/team218/lh20/neuropilot2.sro.rds /lustre/scratch117/cellgen/team218/lh20/inferN0output/neuropilot2.ifn.scp



// c("BCL11B","CUX1","CUX2", "ETV1","FEZF2", "NR4A2","GRIA1","CALB2","ROBO1","ROBO2","GAD2","DLX5","CNTNAP2","SIX3","MEIS2","MEF2C","NPAS3","GRIK2","CCND2","SOX4")



// library(InferN0);  scptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/inferN0output/prefrontalcortex.ifn.scp"); mmm <- Infern0_HierarchicalClustering(scptr,cluster.cells=F); saveInferN0file(scptr"/lustre/scratch117/cellgen/team218/lh20/inferN0output/prefrontalcortex.ifn.scp");


// library(InferN0);  scptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/inferN0outpuprefrontalcortex.ifn.scp");

// library(InferN0);  scptr <- loadInferN0file("/lustre/scratch117/cellgen/team218/lh20/inferN0output/neuropilot2.ifn.scp");genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4","MEF2C","NPAS3","GRIA1","NLGN1","DLX6-AS1","VIM","APOE","APP","PSEN1","STMN2","ROBO1","TOP2A","HIST1H4C","H3F3B","MALAT1","RPLP0","RPL41","RPL12","MT-ND2");genecolors= c(rep("#FF4400",12),rep("#0088FF",9), rep("#FF00FF", 11), rep("#00FF00", 8)); out <-InferN0identifyNetworkExtended(scptr, gene.list=genelist, nb.threads=8, individual.coexp.max=2, color = genecolors,extendedcolor="#AAAAAA", check.min = 40, cell.cluster=10)

/*
  cn <- read.csv("/warehouse/team218_wh01/lh20/SLX-15201/human_MTG_2018-06-14_samples-columns.csv")
  rn <- read.csv("/warehouse/team218_wh01/lh20/SLX-15201/human_MTG_2018-06-14_genes-rows.csv")
  ex <-readRDS("/lustre/scratch117/cellgen/team218/lh20/MDSexspr.rds")
  library("SingleCellExperiment")
  sce <- SingleCellExperiment(assays=list(counts=ex[, colnames(ex) != "X"],logcounts=ex[, colnames(ex) != "X"]))
  sce@assays[["logcounts"]]@x <- log2(sce@assays[["logcounts"]]@x + 1);
  colData(sce)$feature_symbol <- colnames(ex[, colnames(ex) != "X"])
  rowData(sce)$feature_symbol <- rownames(rn)
  library(scmap)
  sce <- selectFeatures(sce, 500,F)
  saveRDS(sce,"MDS.sce.rds")
  colData(sce)$class <- cn$class[colnames(ex) != 'X'] ; colData(sce)$cluster <- cn$cluster[colnames(ex) != 'X']




*/

// library(InferN0); scptr <- InferN0LoadFile("/lustre/scratch117/cellgen/team218/lh20/inferN0output/frozenhigh.ifn.scp");subcell <- readRDS("tmp.rds");







//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernalDensity)]]
Rcpp::List infernalDensity(Rcpp::NumericMatrix data, Rcpp::List inlist){
	uint32_t flags = Rcpp::as<uint32_t> (inlist[0]);
	uint32_t echan =0;
	uint32_t wchan =0;
	uint32_t dchan =0;

	bool has_weight, has_data, normalize_bandwidth;
	if (has_weight=(flags & 1)) wchan=echan++;
	if (has_data=(flags & 2)) dchan=echan++;
	normalize_bandwidth = (flags & 4);

	if (data.ncol() > 2+ echan) return Rcpp::List::create(Rcpp::Named("Error") =  "Excepted Input should have 1 or 2 collumns");
	uint32_t i;

	if (data.ncol() == 2+echan){
		wchan += 2;
		dchan += 2;
		vector<int32_t> dims = Rcpp::as<vector<int32_t> > (inlist[1]);
		vector<double> rect = Rcpp::as<vector<double> > (inlist[2]);
		vector<double> covar = Rcpp::as<vector<double> > (inlist[3]);
		DataGrid<double , 2> buffer;
		DataGrid<double , 2> duffer;
		Tuple<unsigned int,2> cor;
		cor[0] = dims[0]; cor[1] = dims[0];
		buffer.setSizes(cor); buffer.toZero();
		duffer.setSizes(cor); duffer.toZero();

		NormalDistribution<2> dadamat;
		Trianglix<double, 2u> cov;


		KeyElem<Tuple<double,2> , double > valin;
		double pofacts[4];
		pofacts[0] = (-rect[0] * (dims[0]-1))/(rect[2] - rect[0]);
		pofacts[1] = (rect[3] * (dims[0]-1))/(rect[3] - rect[1]);
		pofacts[2] = ((double)(dims[0]-1))/(rect[2] - rect[0]);
		pofacts[3] = -((double)(dims[0]-1))/(rect[3] - rect[1]);
		double totweight=0;
		for(i=0;i< data.nrow();i++){
			valin.k[1] = pofacts[0] + data.at(i,0) *  pofacts[2];
			if ((valin.k[1] < 0)||(valin.k[1] >= dims[0]-1)) continue;
			valin.k[0] = pofacts[1] + data.at(i,1) *  pofacts[3];
			if ((valin.k[0] < 0)||(valin.k[0] >= dims[0]-1)) continue;
			valin.d = (has_weight) ? data.at(i,wchan) : 1.0;
			buffer += valin;
			totweight += valin.d;
			if (has_data){
				valin.d *= data.at(i,dchan);
				duffer += valin;
			}
		}
		buffer /= totweight;
		duffer /= totweight;
		if (normalize_bandwidth) totweight= pow(totweight, -1.0/3.0) * (((double)dims[0]) * dims[0]);
		else totweight= (((double)dims[0]) * dims[0]);
		cov.data[0] = covar[2] * totweight / ((rect[3] - rect[1]) * (rect[3] - rect[1]));
		cov.data[1] = -covar[1] * totweight / ((rect[2] - rect[0]) * (rect[3] - rect[1]));
		cov.data[2] = covar[0] * totweight / ((rect[2] - rect[0]) * (rect[2] - rect[0]));
		// PUTBACK dadamat.setCovariance(cov);
		dadamat.mean.toZero();


		/* = correl[k];
		dadamat.mean[0] =0.0f;
		dadamat.mean[1] =0.0f;

		dascale[0] = (mmm.ims / (mmm.rect[1] - mmm.rect[0])) * pow(correl[k].weight, -1.0f / 6.0f);
		dascale[1] = -(mmm.ims / (mmm.rect[3] - mmm.rect[2])) * pow(correl[k].weight, -1.0f / 6.0f); // (negative because of Y-reversal)

		dadamat *= dascale;*/
//PUTBACK		buffer.blur_zeroborder(dadamat.bindLogLikelihood());

		Rcpp::NumericMatrix outim(dims[0],dims[0]);
		vector<double> outdens; outdens.resize(data.nrow());
		if (!has_data){
			for(cor[1]=0;cor[1] < dims[0]; cor[1]++){
				for(cor[0]=0;cor[0] < dims[0]; cor[0]++) outim.at(cor[0],cor[1]) = buffer(cor) <0 ? 0.0 : buffer(cor);
			}

			for(i=0;i< data.nrow();i++){
				valin.k[1] = pofacts[0] + data.at(i,0) *  pofacts[2];
				if ((valin.k[1] < 0)||(valin.k[1] >= dims[0]-1)) {outdens[i] =0; continue;}
				valin.k[0] = pofacts[1] + data.at(i,1) *  pofacts[3];
				if ((valin.k[0] < 0)||(valin.k[0] >= dims[0]-1)) {outdens[i] =0; continue;}
				outdens[i] = buffer.linearInterpolation(valin.k);
			}
			return(Rcpp::List::create(Rcpp::Named("Query") = Rcpp::wrap(outdens), Rcpp::Named("Density") = outim));
		}else{
			Rcpp::NumericMatrix outim2(dims[0],dims[0]);
			// PUTBACK duffer.blur_zeroborder(dadamat.bindLogLikelihood());
			for(cor[1]=0;cor[1] < dims[0]; cor[1]++){
				for(cor[0]=0;cor[0] < dims[0]; cor[0]++) {
					if (buffer(cor) < pow(0.5,32.0f)){
						outim.at(cor[0],cor[1]) = 0.0;
						outim2.at(cor[0],cor[1]) = 0.0;
					}else{
						outim.at(cor[0],cor[1]) = duffer(cor) / buffer(cor);
						outim2.at(cor[0],cor[1]) = buffer(cor);
					}
				}
			}
			for(i=0;i< data.nrow();i++){
				valin.k[1] = pofacts[0] + data.at(i,0) *  pofacts[2];
				if ((valin.k[1] < 0)||(valin.k[1] >= dims[0]-1)) {outdens[i] =0; continue;}
				valin.k[0] = pofacts[1] + data.at(i,1) *  pofacts[3];
				if ((valin.k[0] < 0)||(valin.k[0] >= dims[0]-1)) {outdens[i] =0; continue;}
				outdens[i] = duffer.linearInterpolation(valin.k) / buffer.linearInterpolation(valin.k);

			}
			return(Rcpp::List::create(Rcpp::Named("Query") = Rcpp::wrap(outdens), Rcpp::Named("Image") = outim, Rcpp::Named("Density") = outim2));
		}
	}
return R_NilValue;}



/*
Rcpp::List principalComponents(const Rcpp::S4 &obj, Rcpp::List inlist){
	uint32_t nbcomp = Rcpp::as<uint32_t> (inlist[0]);
	uint32_t flags = Rcpp::as<uint32_t> (inlist[1]);
	SparseMatrix<double> damat;
	damat.rdRcppdgCMatrix(obj);
	TMatrix<double> out = damat.getPrincipalComponents(nbcomp);
	Rcpp::NumericMatrix output;
	out.wrMatrix(output);
return Rcpp::wrap(output);}
*/


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernalTest)]]
Rcpp::List infernalTest(){

    Forest<uint32_t, 3u> daforest;
    SparseMatrix<double> dadata; dadata.setNBcols(1024);

    uint32_t i,j;
    double da[32];
    for(j=0;j<32;j++) da[j] = sampleGaussian();
    for(i=0; i<256;i++){
        for(j=0;j< 512;j++){
            if (rand() & 8) dadata.data[i][j] = da[j % 31] + da[j % 19] * sampleGaussian();
        }
    }
    for(j=0;j<32;j++) da[j] = sampleGaussian();
    for(; i<512;i++){
        for(j=0;j< 512;j++){
            if (rand() & 8) dadata.data[i][j] = da[j % 31] + da[j % 19] * sampleGaussian();
        }
    }
    for(j=0;j<32;j++) da[j] = sampleGaussian();
    for(; i<768;i++){
        for(j=0;j< 512;j++){
            if (rand() & 8) dadata.data[i][j] = da[j % 31] + da[j % 19] * sampleGaussian();
        }
    }
    for(j=0;j<32;j++) da[j] = sampleGaussian();
    for(; i<1024;i++){
        for(j=0;j< 512;j++){
            if (rand() & 8) dadata.data[i][j] = da[j % 31] + da[j % 19] * sampleGaussian();
        }
    }



    daforest.cluster_likelihood(dadata);

return(Rcpp::List::create());}

//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernalKendallCols)]]
Rcpp::List infernalKendallCols(Rcpp::S4 x){
	SparseMatrix<double> damat;
	damat.rdRcppdgCMatrix(x);
	TMatrix<double> zscores(damat.UTestPaired(true));
	Rcpp::NumericMatrix fout;
	zscores.wrMatrix(fout);
return(Rcpp::List::create(Rcpp::Named("Zscore"), fout));}


//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernalCumulant)]]
Rcpp::List infernalCumulant(Rcpp::List inlist){
	vector<double> data = Rcpp::as<vector<double> > (inlist[0]);
	uint32_t op = Rcpp::as<uint32_t> (inlist[1]);
	uint32_t flags = Rcpp::as<uint32_t> (inlist[2]);
	uint32_t order = Rcpp::as<uint32_t> (inlist[3]);

	Cumulant dacum;
	Rcpp::NumericVector v(16);
	Rcpp::NumericVector stdm(4);

	if (op == 0) dacum.setNormal(data[0], data[1]);
	else if (op == 1) dacum.setNegativeBinomial(data[0], data[1],order);
	else if (op == 2) dacum.setGamma(data[0], data[1], order);
	else if (op == 3) dacum.setMannWhitney((int)data[0], (int)data[1]);
	for(int i =0; i < dacum.data.getSize();i++) v[i] = dacum.data[i];
	Tuple<double, 4u> stdmom = dacum.getStdMoments();
	for(int i = 0;i < 4;i++) stdm[i] = stdmom[i];

	double logQ;
	if (flags !=0){
		if (flags < 100) {
			Tuple<double> tbuf;
			tbuf.setSize(flags);
			for(int i=0;i<flags;i++) tbuf[i] = dacum.data[i];
			dacum.data.setSize(flags);
			for(int i=0;i<flags;i++)  dacum.data[i] = tbuf[i];
			logQ = dacum.evalLogQuantile(data[3], true);
		}else if (flags < 1200) {
			Tuple<double> tbuf;
			tbuf.setSize(flags - 10);
			for(int i=0;i<tbuf.getSize();i++) tbuf[i] = dacum.data[i];
			dacum.data.setSize(flags - 10);
			for(int i=0;i<tbuf.getSize();i++)  dacum.data[i] = tbuf[i];
			logQ = (-dacum).evalLogQuantileSaddle(-data[3]);
		}else{
			logQ = (-dacum).evalLogQuantileOld(-data[3]);
		}
	}else logQ = dacum.evalLogQuantile(data[3], true);

//	printf("Sanity check quantile %e == %e ? (val = %e)\n", data[2], dacum.evalQuantile(dacum.fromQuantile(data[2])), dacum.fromQuantile(data[2]));
//	printf("Sanity check logquantile %e == %e ? (val = %e)\n", log(data[2]), dacum.evalLogQuantile(dacum.fromLogQuantile(log(data[2]))), dacum.fromLogQuantile(log(data[2])));
//	printf("Sanity check value %e == %e ? (quant = %e)\n", data[3], dacum.fromQuantile(dacum.evalQuantile(data[3])), dacum.evalQuantile(data[3]));
//	printf("Sanity check value %e == %e ? (logquant = %e)\n", data[3], dacum.fromLogQuantile(dacum.evalLogQuantile(data[3])), dacum.evalLogQuantile(data[3]));

	double target = exp(dacum.evalLogLikelihood(data[3]));
//	printf("%e is the target\n", target);
//	for(int i=0;i<10;i++){
//		double diff = (dacum.evalQuantile(data[3]+ exp(-4.0 - i)) -  dacum.evalQuantile(data[3] - exp(-4.0-i))) / (2.0 * exp(-4.0-i));
//		printf("%e is tangent\n", diff);
//	}
	return Rcpp::List::create(Rcpp::Named("K") = Rcpp::wrap(v), Rcpp::Named("sample") = Rcpp::wrap(dacum.evalQuantile(data[2])), Rcpp::Named("stdmom") = stdm, Rcpp::Named("density") = Rcpp::wrap(dacum.evalDensity(data[3])), Rcpp::Named("LogLike") = Rcpp::wrap(dacum.evalLogLikelihood(data[3])), Rcpp::Named("Quantile") = Rcpp::wrap(exp(logQ)), Rcpp::Named("LogQuantile") = Rcpp::wrap(logQ));
}



//' Internal Function, see R code for documentation instead
//' @export
// [[Rcpp::export(.infernalDeconvolve2D)]]
Rcpp::List infernalDeconvolve2D(Rcpp::S4 scdata_in, Rcpp::S4 bulk_in, Rcpp::List inlist){ class Task{ public:
	SparseMatrix<uint32_t> scdata;
 	SparseMatrix<uint32_t> tr_scdata;

	SparseMatrix<double> sc_ct_state;
	SparseMatrix<double> sc_ct_state_deriv;
	SparseMatrix<double> sc_ct_state_zerobuf_M;
	SparseMatrix<double> sc_ct_state_zerobuf_V;

	SparseMatrix<uint32_t> bkdata;
	std::vector<uint32_t> ct_index;
	std::vector<uint32_t> dn_index;
	Rcpp::CharacterVector tmp_annotation;
	Rcpp::CharacterVector ct_names;
	Rcpp::CharacterVector dn_names;
	Rcpp::CharacterVector gene_names;
	Rcpp::CharacterVector vgene_names;
	vector<uint32_t> genefilter;

	vector<double> LLoutput_buffer;
	Rcpp::CharacterVector bulknames;
	double averagedepth;
	uint32_t curtask;
	Tuple<double> bk_scale;
	Vector<TMatrix<double> > llout;
	int nbcelltypes, nbgenotypes, nbgenes;
	Tuple< LFHPrimitive::WeightElem<double, 2> > dasuperlight_prior;
	Tuple<uint32_t> sc_depth, sc_genecount;

	Vector< Tuple<double, 2u> > lonely_MV_param;
	Vector< TMatrix<double> > M_param;
	Vector< TMatrix<double> > V_param;
	Vector< TMatrix<double> > M_param_deriv;
	Vector< TMatrix<double> > V_param_deriv;
	Vector< TMatrix<double> > outer;


	myHashmap<string, uint32_t> genemap;


	//x_i ~ N( mu(s_i), \sigma(s_i))  -> c + B x_i ~ N(c + Bmu , B Sigma B^t )


    GaussianProcessArray dagp;
    DataGrid<double, 4> dagpinput; // PCA-reduced

	MemoryLimitedOptimizationScope adl;
	uint32_t gradient_search_steps;
	double zerodderbuf[6];

	Tuple< uint32_t> thrbuf_warn;
	Tuple< Tuple<double, 2u> > thrbuf_LL;
	Tuple< Tuple<double, 2u> > thrbuf_sum_mv;
	Tuple< Tuple<TMatrix<double>, 2u> > thrbuf_totsum_mv; // matrix sum
	Tuple< TMatrix<double> > thrbuf_totsum_df; // outerproduct sum

	bool doload, useold,useinter, useold2 ,usegp;

	int mainloopcount;
	Tuple<uint32_t> dasumcheck;
	myHashmap< uint32_t, Tuple< GaussianProcessMK2> > gp_mv; // gene x celltype
	myHashmap< uint32_t, Tuple< Tuple< double> > > gp_mv_dx; // gene x celltype
	Vector< TMatrix<double> > M_serror; // noise

	TMatrix<double> gpstate;
	Tuple< TMatrix<double> > thrbuf_gpstate_dx;

	ThreadBase::ThreadArrayScope thrarray;

	TMatrix<double> dn_deconv_results, ct_deconv_results;
	vector<double> intercept;
	Task():curtask(0), nbcelltypes(0), averagedepth(0.0),thrarray(foreach,16){}

Rcpp::List run(int nbsteps, const char* savepath){
	if (ct_index.size() != scdata.getNBcols()){
		printf("unexpected annotation length %i (!= %i)\n", (int)ct_index.size() ,scdata.getNBcols());
		return(Rcpp::List());
	}

	SparseMatrix<uint32_t> annot_count;
	Tuple<uint32_t, 2> coor, sccoor;


	for(int i =0;i<ct_index.size();i++){
		if ( (coor[1] = dn_index[i]) >= annot_count.data.getSize()) annot_count.data.upSize(dn_index[i]+1);
		if ( (coor[0] = ct_index[i]) >= nbcelltypes) nbcelltypes = ct_index[i]+1;
		annot_count(coor) += 1;
	}
	nbgenotypes = annot_count.data.getSize();

	int nbsinglecell = scdata.getNBcols();


	// Step 1, retrieve NB parameter for all genes, donors and celltypes
	scdata.purgeValue(0);
	tr_scdata = scdata.mkTranspose();

	sc_depth = scdata.sumRows();
	sc_genecount = scdata.sumCols();
	printf("some depths... %i %i %i len = %i\n", sc_depth[0], sc_depth[1], sc_depth[2], sc_depth.getSize());
	printf("some moredepths... %i %i %i len = %i\n", sc_genecount[0], sc_genecount[1], sc_genecount[2], sc_genecount.getSize());
	nbgenes = sc_genecount.getSize();
	double	prior_weight[3];
	prior_weight[1] =0.0;
	for(int j=0;j<sc_genecount.getSize();j++) if (genefilter[j] == 1) prior_weight[1] += 1.0;

	printf("Deconvoling %i bulks with %i/%i genes using %i single cells from  %i donors and with %i celltypes\n", bkdata.getNBcols(),(int)prior_weight[1],nbgenes,scdata.getNBcols(),nbgenotypes, nbcelltypes);

	if (doload){
		FILE *f = fopen(savepath, "rb+");
		if (f == NULL) {printf("could not open %s\n", savepath);return Rcpp::List::create(); }
		M_param.load(f);
		V_param.load(f);
		sc_ct_state.load(f);
		lonely_MV_param.load(f);
		fclose(f);
	}else{
		M_param.setSize(nbgenes); V_param.setSize(nbgenes);lonely_MV_param.setSize(nbgenes);sc_ct_state.setNBcols(scdata.getNBcols());
		for(int i = 0; i< nbgenes;i++){
			M_param[i].setSizes(nbcelltypes+1,nbgenotypes);
			V_param[i].setSizes(nbcelltypes+1,nbgenotypes);
		}
	}
	outer.setSize(nbgenes);
	M_param_deriv.setSize(nbgenes); V_param_deriv.setSize(nbgenes);
	for(int i = 0; i< nbgenes;i++) outer[i].setSizes(nbcelltypes+1,nbgenotypes);


	dasuperlight_prior.setSize(nbgenes).toZero();
	prior_weight[1] = 0.000001;

	for(int j =0;j<scdata.getNBcols();j++) averagedepth += log((double)sc_depth[j]);
	averagedepth = exp(averagedepth / scdata.getNBcols()); // use geometric average has reference

	sc_ct_state_deriv.setNBcols(scdata.getNBcols()); sc_ct_state_zerobuf_M.setNBcols(scdata.getNBcols()); sc_ct_state_zerobuf_V.setNBcols(scdata.getNBcols());
	for(int j =0;j<scdata.getNBcols();j++){
		if (auto ite =scdata.data[j]()) do { dasuperlight_prior[ite()] += LFHPrimitive::WeightElem<double,2>( ((double)*ite) * (averagedepth / sc_depth[j]),  prior_weight[1]); }while(ite++);
		if (doload == false){
			sc_ct_state.data[j][ct_index[j]] = (averagedepth / sc_depth[j]) * 0.9;
			sc_ct_state.data[j][nbcelltypes] = (averagedepth / sc_depth[j]) * 0.1;
		}
		sc_ct_state_deriv.data[j][ct_index[j]] = ((averagedepth / sc_depth[j]) + 1.0) / 2;
		sc_ct_state_deriv.data[j][nbcelltypes] = (sc_ct_state.data[j][nbcelltypes] + 1.0)  /2;
	}
	for(int i=0; i < dasuperlight_prior.getSize(); i++) dasuperlight_prior[i].w[0] = prior_weight[1] * scdata.getNBcols();

	Tuple<uint32_t> bk_depth = bkdata.sumRows();

	bk_scale.setSize(bk_depth.getSize());
	for(int i=0;i<bk_depth.getSize();i++) bk_scale[i] = ((double) bk_depth[i]) / averagedepth;

	llout.setSize(bkdata.getNBcols());
	if (doload == false){
		thrarray.initRange(nbcelltypes * nbgenotypes, "Init Parameters Heuristic");
		curtask = 0; thrarray.submit(*this);
	}
	thrbuf_totsum_mv.setSize(thrarray.nbthreads); thrbuf_totsum_df.setSize(thrarray.nbthreads);
	dasumcheck.setSize(thrarray.nbthreads);

	thrbuf_warn.setSize(thrarray.nbthreads);

	thrbuf_LL.setSize(thrarray.nbthreads);
	thrbuf_sum_mv.setSize(thrarray.nbthreads);
	for(int i=0;i < thrbuf_totsum_mv.getSize();i++) {thrbuf_totsum_mv[i][0].setSizes(nbcelltypes+1,nbgenotypes);thrbuf_totsum_mv[i][1].setSizes(nbcelltypes+1,nbgenotypes);thrbuf_totsum_df[i].setSizes(nbcelltypes+1,nbgenotypes);}

	curtask = 1; thrarray.initRange(nbgenes, "Init Parameters For Ambiant RNA and prior"); thrarray.submit(*this);
	thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum(); thrbuf_LL[0] = thrbuf_LL.mkSum();
	printf("prior LL %e -> %e\n", thrbuf_LL[0][0], thrbuf_LL[0][1]);

	int valtocheck = 0;
	double tocheckLL;

	int failcount=0;

/*
	printf("comp that silly %i \t %i\n", nbgenes, tr_scdata.getNBcols());
	adl.initVariable("Mmat", M_param_deriv.mkPartitionedMetaIterator(),1).initVariable("Vmat", V_param_deriv.mkPartitionedMetaIterator(),1);
	adl.initVariable("cellstate", sc_ct_state_deriv.mkIterator(),0).initLearn();
	for(int ss =0; ; ss++){
		Rcpp::checkUserInterrupt();

//		if ((adl.step & 1) == 1){
//			curtask = 7; thrarray.initRange(nbgenes, "update");thrarray.submit(*this); thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum();
//			curtask = 2; thrarray.initRange(scdata.getNBcols(), "EM step: cell-wise"); thrarray.submit(*this);
//			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_sum_mv[0] = thrbuf_sum_mv.mkSum(); thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum();
//
//			thrbuf_sum_mv[0][0] = thrbuf_totsum_mv[0][0].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][0];
//			thrbuf_sum_mv[0][1] = thrbuf_totsum_mv[0][1].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][1];
//			NBmeanvar_logprob(zerodderbuf, 0u, thrbuf_sum_mv[0][0], thrbuf_sum_mv[0][1]);
//
//			printf("Ok, model param for zeros is m=%e v=%e: LL=%e\n", thrbuf_sum_mv[0][0] ,thrbuf_sum_mv[0][1], zerodderbuf[0]);
//			thrbuf_LL[0][0] += zerodderbuf[0];
//		}else{
//			curtask = 6; thrarray.initRange(scdata.getNBcols(), "update");thrarray.submit(*this); thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum();
//			curtask = 3; thrarray.initRange(nbgenes, "EM step: gene-wise"); thrarray.submit(*this);
//			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_sum_mv[0] = thrbuf_sum_mv.mkSum(); thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum();
//
//
//			thrbuf_sum_mv[0][0] = thrbuf_totsum_mv[0][0].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][0];
//			thrbuf_sum_mv[0][1] = thrbuf_totsum_mv[0][1].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][1];
//			NBmeanvar_logprob(zerodderbuf, 0u, thrbuf_sum_mv[0][0], thrbuf_sum_mv[0][1]);
//
//			printf("Ok, model param for zeros is m=%e v=%e: LL=%e\n", thrbuf_sum_mv[0][0] ,thrbuf_sum_mv[0][1], zerodderbuf[0]);
//			thrbuf_LL[0][0] += zerodderbuf[0];
//		}
//		printf("daother LL %e\n", thrbuf_LL[0][0]); M_param_deriv.toZero(); V_param_deriv.toZero(); if (auto ite = sc_ct_state_deriv.mkIterator()) do{*ite = 0.0;} while(ite++);

		if ((adl.step & 1) == 0){
			curtask = 7; thrarray.initRange(nbgenes, "update");thrarray.submit(*this); thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum();
			curtask = 2; thrarray.initRange(scdata.getNBcols(), "EM step: cell-wise"); thrarray.submit(*this);
			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_sum_mv[0] = thrbuf_sum_mv.mkSum(); thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum();

			thrbuf_sum_mv[0][0] = thrbuf_totsum_mv[0][0].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][0];
			thrbuf_sum_mv[0][1] = thrbuf_totsum_mv[0][1].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][1];
			NBmeanvar_logprob(zerodderbuf, 0u, thrbuf_sum_mv[0][0], thrbuf_sum_mv[0][1]);

			printf("Ok, model param for zeros is m=%e v=%e: LL=%e\n", thrbuf_sum_mv[0][0] ,thrbuf_sum_mv[0][1], zerodderbuf[0]);
			curtask = 4; thrarray.initRange(scdata.getNBcols(), "EM step: cell-wise zeros");thrarray.submit(*this);
			thrbuf_LL[0][0] += zerodderbuf[0];
		}else{
			curtask = 6; thrarray.initRange(scdata.getNBcols(), "update");thrarray.submit(*this); thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum();
			curtask = 3; thrarray.initRange(nbgenes, "EM step: gene-wise"); thrarray.submit(*this);
			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_sum_mv[0] = thrbuf_sum_mv.mkSum(); thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum();


			thrbuf_sum_mv[0][0] = thrbuf_totsum_mv[0][0].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][0];
			thrbuf_sum_mv[0][1] = thrbuf_totsum_mv[0][1].vectorize().mkInnerProd(thrbuf_totsum_df[0].vectorize()) - thrbuf_sum_mv[0][1];
			NBmeanvar_logprob(zerodderbuf, 0u, thrbuf_sum_mv[0][0], thrbuf_sum_mv[0][1]);


			printf("Ok, model param for zeros is m=%e v=%e: LL=%e\n", thrbuf_sum_mv[0][0] ,thrbuf_sum_mv[0][1], zerodderbuf[0]);
			curtask = 5; thrarray.initRange(nbgenes, "EM step: gene-wise zeros");thrarray.submit(*this);
			thrbuf_LL[0][0] += zerodderbuf[0];
		}
		printf("sum count is %i\n", dasumcheck.mkSum());

		if (adl.doesAccept(thrbuf_LL[0][0])){
			failcount=0;

			printf("LL[%i]: %e\n",adl.step, thrbuf_LL[0][0]);
			if (ss >= nbsteps) break; // All good
			if ((adl.step & 1) == 1) {
				printf("cellstate max rel change: %e\n", adl.updatePositiveVariable("cellstate",sc_ct_state.mkIterator(), sc_ct_state_deriv.mkIterator()));
			}else{
				printf("Vmat max rel change: %e\n", adl.updatePositiveVariable("Vmat",V_param.mkPartitionedMetaIterator(), V_param_deriv.mkPartitionedMetaIterator()));
				printf("Mmat max rel change: %e\n", adl.updatePositiveVariable("Mmat",M_param.mkPartitionedMetaIterator(), M_param_deriv.mkPartitionedMetaIterator()));
			}

		}else{
			if (++failcount > 5){
				printf("too many fails\n");
				if ((adl.step & 1) == 1) adl.markVariableToCheck("Mmat", 1);
				else adl.markVariableToCheck("cellstate", 0);
			}

			if (ss >= nbsteps+10) break; // tryagain...
			printf("LL[%i]: %e (Rejected!)\n",adl.step, thrbuf_LL[0][0]);
			if ((adl.step & 1) == 1) {
				adl.backtrackVariable("cellstate",sc_ct_state.mkIterator());
				M_param_deriv.toZero(); V_param_deriv.toZero();
			}else{
				adl.backtrackVariable("Mmat",M_param.mkPartitionedMetaIterator());
				adl.backtrackVariable("Vmat",V_param.mkPartitionedMetaIterator());
				if (auto ite = sc_ct_state_deriv.mkIterator()) do{*ite = 0.0;} while(ite++);
			}
		}
	}*/


	//curtask =0; thrbase.submit_Array(*this, 8);
	//curtask =1; thrbase.submit_Array(*this, 1);

	//printf("deconvulve for real now!\n");
	//if (!useold){
		for(mainloopcount=0;mainloopcount<nbsteps;mainloopcount++){
			Rcpp::checkUserInterrupt();
			curtask = 12; thrarray.initRange(scdata.getNBcols(), "EM step: cell-wise"); thrarray.submit(*this);
			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum(); thrbuf_warn[0] = thrbuf_warn.mkSum();
			if (!ExOp::isValid(thrbuf_totsum_df[0])){
				printf("grg\n")	;
				thrbuf_totsum_df[0].show();
			}
			printf("%i.A: LL %e -> %e (%i warnings)\n", mainloopcount, thrbuf_LL[0][0], thrbuf_LL[0][1], thrbuf_warn[0]);

			Rcpp::checkUserInterrupt();
			curtask = 14; thrarray.initRange(nbgenes, "EM step: gene-wise"); thrarray.submit(*this);
			thrbuf_LL[0] = thrbuf_LL.mkSum(); thrbuf_totsum_mv[0] = thrbuf_totsum_mv.mkSum(); thrbuf_warn[0] = thrbuf_warn.mkSum();
			if (!ExOp::isValid(thrbuf_totsum_mv[0])){
				printf("grg M \n");	thrbuf_totsum_mv[0][0].show();
				printf("grg V\n");thrbuf_totsum_mv[0][1].show();
			}
			printf("%i.B: LL %e -> %e (%i warnings)\n", mainloopcount, thrbuf_LL[0][0], thrbuf_LL[0][1], thrbuf_warn[0]);

		}
	//}
	//printf("did exit...\n");
	//exit(1);

	if (doload == false){
        if (strlen(savepath) != 0){
            FILE *f = fopen(savepath, "wb+");
            if (f == NULL) {printf("could not open %s\n", savepath);return Rcpp::List::create(); }
            M_param.save(f);
            V_param.save(f);
            sc_ct_state.save(f);
            lonely_MV_param.save(f);
            fclose(f);
        }
	}
	dn_deconv_results.setSizes(nbgenotypes,bulknames.size());
	ct_deconv_results.setSizes(nbcelltypes,bulknames.size());
	intercept.resize(bulknames.size());
	mainloopcount=0;
	for(int i=0;i<M_param.getSize();i++){
		bool isok = false;
		if ((ExOp::isValid(M_param[i]) )&&(ExOp::isValid(V_param[i]) )) {
			// check if gene is 100% bg
			if (auto ite = M_param[i].getIterator()){
				do{
					if (ite()[0] == nbcelltypes) continue;
					if (*ite != 0.0) {isok = true; break;}
				}while(ite++);
			}
			if (auto ite = V_param[i].getIterator()){
				do{
					if (ite()[0] == nbcelltypes) continue;
					if (*ite != 0.0) {isok = true; break;}
				}while(ite++);
			}
		}
		printf("%i: %i (vs %i)\n", i, tr_scdata.data[i].getSize() , scdata.getNBcols() / 2);
		if (tr_scdata.data[i].getSize() * 2 < scdata.getNBcols()) isok = false;
		if (!isok) {if (genefilter[i] != 0) mainloopcount++; genefilter[i] = 0;}
		if (genefilter[i] != 0) {
            genemap[string(gene_names[i])] = i;
            vgene_names.push_back(gene_names[i]);
        }
	}
	printf("filtered %i genes with undefined param, remains %i\n",mainloopcount,genemap.getSize());
	LLoutput_buffer.resize(bkdata.getNBcols());
	//curtask =2; thrbase.submit_Array(*this, 8);

	Rcpp::List lfout = Rcpp::List::create();
	Rcpp::NumericMatrix fout;//= arma::eye<arma::mat>(hehe.getSize(), hehe.getSize());

	if (usegp == false){
	//	myHashmap<uint32_t, uint32_t> show_these_counts;
	//	if (auto ite = bkdata.getIterator()) do{
	//		show_these_counts[*ite]++;
	//	}while (ite++);
	//	show_these_counts.show();
		if(useold2){
			thrarray.initRange(bkdata.getNBcols(), "Deconvonvolve Bulk");
			curtask=8;thrarray.submit(*this);
		}else{
			thrarray.initRange(bkdata.getNBcols(), "Deconvonvolve Bulk");
			curtask=9;thrarray.submit(*this);
		}
	}else{

		curtask = 6; thrarray.initRange(scdata.getNBcols(), "Recompute NBcell per combination"); thrarray.submit(*this);
		thrbuf_totsum_df[0] = thrbuf_totsum_df.mkSum();
        thrbuf_totsum_df[0].show();

        Tuple<unsigned int, 4u> coor4; coor4[0] = 2; coor4[1] = nbgenotypes; coor4[2] = genemap.getSize(); coor4[3] = nbcelltypes;
        dagpinput.setSizes(coor4);

        printf("%i genes x %i celltype\n", genemap.getSize() , nbcelltypes); fflush(stdout);
        curtask = 23; thrarray.initRange(genemap.getSize(), "Init GP input"); thrarray.submit(*this);
        printf("done\n"); fflush(stdout);

        TMatrix<double> daclone((*dagpinput).vectorizeOntoDimension(2).selectSlice(0,0));
        daclone.wrMatrix(fout);
        Rcpp::rownames(fout) = dn_names;
        lfout["GPinput_mean"] = fout;

        daclone = TMatrix<double>((*dagpinput).vectorizeOntoDimension(2).selectSlice(1,0));
        daclone.wrMatrix(fout);
        Rcpp::rownames(fout) = dn_names;
        lfout["GPinput_var"] = fout;

        //dagp.setSize(coor3[2], 2);
        //dagp.initUsingPCA(dagpinput,2);

        //dagp.hstate.wrMatrix(fout);
        //Rcpp::colnames(fout) = dn_names;
        //Rcpp::rownames(fout) = Rcpp::CharacterVector({"PC_1","PC_2"});
        //lfout["Donor_PCA_State"] = fout;

        printf("reduction start\n"); fflush(stdout);
        dagp.learnPCAreduction(dagpinput,2, 10, 10);

        printf("Learned hidden state for reference:\n"); dagp.hstate.show();

        dagp.hstate.wrMatrix(fout);
        Rcpp::colnames(fout) = dn_names;
        Rcpp::rownames(fout) = Rcpp::CharacterVector({"GpLat_1","GpLat_2"});
        lfout["GP_State_Reference"] = fout;


        dagp.learnGPParams((*dagpinput).vectorizeOntoDimension(2));


        daclone = dagp.makeGParammatrix();
        for(int i=0;i<ct_names.size();i++){
            TMatrix<double> somematrix((*daclone).splitDimension(0,genemap.getSize()).selectSlice(i,1));
            somematrix.wrMatrix(fout);
            Rcpp::colnames(fout) = Rcpp::CharacterVector({"Mean","Noise","Signal","K00","K01","K11"});
            Rcpp::rownames(fout) = vgene_names;
            lfout[string("GP_Parameters_") + ct_names[i]] = fout;
        }

     /*
        uint32_t dagenegetto = genemap.find(string("PAEP"));
        if (dagenegetto == 0xFFFFFFFF) printf("did not find gene\n");
        else printf("%i is gene offset\n", dagenegetto);
        if (dagenegetto != 0xFFFFFFFF){
            TiffFile tf("GPput.tif", true);
            DataGrid<double, 3u> daimage;
            Tuple<double> extra;
            for(int iiii=0; iiii < nbcelltypes;iiii++){
                daimage = dagp[iiii + nbcelltypes * dagenegetto].makeImage(dagpinput.getPixelsOnRow(iiii + nbcelltypes * dagenegetto),dagp.hstate,5u<<6,5u<<6,2u<<6,2u<<6, &extra);
                printf("%i: extra is \n", iiii); extra.show(); printf("make image:\n");
                tf.put(daimage.applyFunctionOnPixels<unsigned char>([r = extra()+4](const double* o, unsigned char* pix){
                                    double phase = sin(((o[0] - r[0]) / (r[1] - r[0])) * M_PI - M_PI * 0.5);
                                    double scstd = (r[3] - o[1]) / (r[3] - r[2]);
                                    double hue = phase * (2.125 - scstd*scstd * 1.625);
                                    Tuple<double, 3u> pix2 = RGBfromHCY_gamma((hue + 2.0) / 6.0, abs(hue) < 1.0 ? abs(hue): 1.0, scstd*scstd + (1.0 - scstd*scstd) * 0.25 * phase * phase);
                                    pix[0] = (unsigned char)(pow(pix2[0],0.454f) * 255.0); pix[1] = (unsigned char)(pow(pix2[1],0.454f) * 255.0); pix[2] = (unsigned char)(pow(pix2[2],0.454f) * 255.0);}, 3));
            }
        }*/

        TMatrix<uint32_t> bulksubset;
        bulksubset.setSizes(genemap.getSize(),bkdata.getNBcols()).toZero();
        printf("Write Matrix Subset %i x %i\n", genemap.getSize() , bkdata.getNBcols());
        for(coor[0]=0; coor[0] < genemap.getSize();coor[0]++){
            int currow = genemap.deref(coor[0]);
            int tmpoff;
            for(coor[1] = 0 ; coor[1] < bkdata.getNBcols();coor[1]++){
                if ((tmpoff = bkdata.data[coor[1]].find(currow)) != 0xFFFFFFFF) bulksubset(coor) = bkdata.data[coor[1]].deref(tmpoff);
            }
        }

        printf("Run regressions\n");
        Tuple< TMatrix<double>, 3u > dares = dagp.localPoissonRegression(bulksubset);

        printf("Write that one then!\n");

        dares[0].wrMatrix(fout);
        Rcpp::colnames(fout) = bulknames;
        Rcpp::rownames(fout) = Rcpp::CharacterVector({"GpLat_1","GpLat_2","GpLat_prt","LL"});
        lfout["Donor_deconv"] = fout;
        dares[1].wrMatrix(fout);
        Rcpp::colnames(fout) = bulknames;
        Rcpp::rownames(fout) = ct_names;
        lfout["CellType_deconv"] = fout;
        dares[2].wrMatrix(fout);
        Rcpp::colnames(fout) = bulknames;
        Rcpp::rownames(fout) = vgene_names;
        lfout["Inferred_mean"] = fout;


        printf("DONE!\n");

        /*
        M_serror.setSize(M_param.getSize());
		for(int i=0;i<M_param.getSize();i++){
			if (genefilter[i] != 0 ) {
				gp_mv[i].setSize(nbcelltypes);
				gp_mv_dx[i].setSize(nbcelltypes);
			}
		}
		thrbuf_gpstate_dx.setSize(thrarray.nbthreads);

		gpstate.setSizes(2,nbgenotypes).toRand();


		thrarray.initRange(M_param.getSize(), "Init GP scope"); thrbuf_LL[0] = thrbuf_LL.mkSum();
		curtask =20; thrarray.submit(*this);
        printf("so?\n");
		for(mainloopcount=0;mainloopcount<4;mainloopcount++){
			Rcpp::checkUserInterrupt();
            printf("so??\n");
			thrarray.initRange(M_param.getSize(), "Update GP params");
			curtask =21; thrarray.submit(*this); thrbuf_LL[0] = thrbuf_LL.mkSum();
			printf("%i.C: LL %e -> %e (%i warnings)\n", mainloopcount, thrbuf_LL[0][0], thrbuf_LL[0][1], thrbuf_warn[0]);

            printf("so???\n");
			Rcpp::checkUserInterrupt();
			thrarray.initRange(M_param.getSize(), "Update GP state");
			curtask =22; thrarray.submit(*this);
			thrbuf_gpstate_dx[0] = thrbuf_gpstate_dx.mkSum();thrbuf_LL[0] = thrbuf_LL.mkSum();
			printf("%i.D: LL %e -> %e (%i warnings)\n", mainloopcount, thrbuf_LL[0][0], thrbuf_LL[0][1], thrbuf_warn[0]);
		}

		gpstate.wrMatrix(fout);
		Rcpp::colnames(fout) = dn_names;
		lfout["Donor_GP_state"] = fout;
		*/

		//samplingnoise.wrMatrix(fout);
		//Rcpp::colnames(fout) = dn_names;
		//Rcpp::rownames(fout) = ct_names;
		//lfout["SamplingVariance"] = fout;


	}

	printf("alldone!!!\n"); fflush(stdout);
	foreach.stopThreadArray();
	TMatrix<double> summy[2];
	curtask =0;
	while (genefilter[curtask] ==0) curtask++;
	summy[0] = M_param[curtask];
	summy[1] = V_param[curtask];

	for(curtask++;curtask < nbgenes;curtask++){
		if (genefilter[curtask] == 0) continue;
		summy[0] += M_param[curtask];
		summy[1] += V_param[curtask];
	}
	printf("show mean\n"); summy[0].show();
	printf("show var\n"); summy[1].show();
	ct_names.push_back("Background");
	summy[0].wrMatrix(fout);
	Rcpp::colnames(fout) = dn_names; Rcpp::rownames(fout) = ct_names;
	lfout["Mean_Param"] = fout;
	summy[1].wrMatrix(fout);
	Rcpp::colnames(fout) = dn_names; Rcpp::rownames(fout) = ct_names;
	lfout["Var_Param"] = fout;
	char danum[256];

	if (usegp == false){

        dn_deconv_results.wrMatrix(fout);
        Rcpp::colnames(fout) = bulknames;
        Rcpp::rownames(fout) = dn_names;
        lfout["Donor_deconv"] = fout;
        ct_deconv_results.wrMatrix(fout);
        Rcpp::colnames(fout) = bulknames;
        Rcpp::rownames(fout) = ct_names;
        lfout["CellType_deconv"] = fout;


        TMatrix<double> massive; massive.setSizes(nbgenes,(nbcelltypes+1) *nbgenotypes);
        for(coor[0]=0;coor[0]<nbgenes;coor[0]++){
            coor[1]=0;
            if (auto ite = M_param[coor[0]].getIterator()) do{
                massive(coor) = *ite;
                coor[1]++;
            }while(ite++);
        }
        massive.wrMatrix(fout);
        Rcpp::rownames(fout) = gene_names;
        lfout["Learned_Mean"] = fout;
    }
/*
    if (usegp != false){
    	massive.toZero();
        for(coor[0]=0;coor[0]<nbgenes;coor[0]++){
            coor[1]=0;
            if (auto ite = M_serror[coor[0]].getIterator()) do{
                massive(coor) = *ite;
                coor[1]++;
            }while(ite++);
        }
        massive.wrMatrix(fout);
        Rcpp::rownames(fout) = gene_names;
        lfout["Sampling_Error"] = fout;
    }*/


    /*
	lfout["Bulk_LogLikelihood"] = Rcpp::wrap(LLoutput_buffer);
	printf("warp\n");
	lfout["Bulk_Intercept"] = Rcpp::wrap(intercept);
	vector<double> sc_bgfrac; sc_bgfrac.resize(sc_ct_state.getNBcols());
	double updown[2];
	for(int i = 0 ; i <sc_ct_state.getNBcols();i++){
		if (auto ite = sc_ct_state.data[i]()) do{
			if (ite() < nbcelltypes) updown[0] = *ite;
			else updown[1] = *ite;
		}while(ite++);
		sc_bgfrac[i] = updown[1] / (updown[1] + updown[0]);
	}

	lfout["sc_bg_fraction"] = Rcpp::wrap(sc_bgfrac);
	*/

	/*
	printf("danamelen %i %i %i\n",ct_names.size(), dn_names.size(), bulknames.size());
	for(int i = 0 ; i< llout.getSize();i++){
		llout[i].wrMatrix(fout);
		Rcpp::colnames(fout) = dn_names;
		Rcpp::rownames(fout) = ct_names;
		printf("%s\n", Rcpp::as<string>(bulknames[i]).c_str());
		lfout[Rcpp::as<string>(bulknames[i]).c_str()] = Rcpp::wrap(fout);
	}*/
return lfout;}

   	operator std::function<uint32_t(uint32_t)> (){using std::placeholders::_1; return std::bind(&Task::operator(),this,_1);}

	uint32_t operator()(uint32_t threadID){
		uint32_t jobID;
		double dderbuf[3];
		CurvatureSearchScope dacurv;
		int innerloop;
		double LL;
		uint32_t i;
		switch(curtask){
			case 0:{
			Vector< uint32_t > dacelllist[3];
			double prior_weight[3];  prior_weight[0] = 1.0;
			while(thrarray.getRangeVal(jobID)){
			//if (auto jobite = thrbase.getRangeIterator(threadID, nbcelltypes * nbgenotypes,true)) do{
				Tuple<uint32_t, 2u> coor; coor[0] = jobID % (nbcelltypes); coor[1] = jobID / (nbcelltypes);
				for(i =0 ; i < scdata.getNBcols(); i++ ){
					if (ct_index[i] == coor[0]) {
						if (dn_index[i] == coor[1]) dacelllist[0].push_back(i);
						else dacelllist[1].push_back(i);
					}else if (dn_index[i] == coor[1]) dacelllist[2].push_back(i);
				}
				prior_weight[1] = 0.05;
				prior_weight[2] = 0.05;
				Tuple< LFHPrimitive::WeightElem<double, 2> > daelem_buff = dasuperlight_prior;
				for(int j=0;j<3;j++){
					for(int i =0; i < dacelllist[j].getSize();i++) {
						if (auto ite =scdata.data[dacelllist[j][i]]()) do {
							daelem_buff[ite()] += LFHPrimitive::WeightElem<double,2>( ((double)*ite) * (averagedepth / sc_depth[ dacelllist[j][i] ]),  prior_weight[j]);
							if ((*ite) == 0) printf("found a zero!\n");
							if ((*ite) >= 0x0FFFFFFF) printf("found an extreme val %i or %X!\n", (*ite) , (*ite) );
						}while(ite++);
					}
				}
				prior_weight[1] =  prior_weight[2] *  dacelllist[2].getSize() +  prior_weight[1] * dacelllist[1].getSize()  + dacelllist[0].getSize() + dasuperlight_prior[0].w[0];
				prior_weight[2] =0;

				for(int j=0;j<daelem_buff.getSize();j++) {
					//if (sc_genecount[j] < mincount_pergene) {
					//	M_param[j](coor) = 0.0;
					//	V_param[j](coor) = 0.0;
					//	continue;
					//}
					daelem_buff[j].w[0] = prior_weight[1];
					M_param[j](coor) = daelem_buff[j].getMean();
					if (!ExOp::isValid(M_param[j](coor))){ // gene has ZERO counts!
						M_param[j](coor) = 0.0;
						V_param[j](coor) = 0.0;
						continue;
					}
					V_param[j](coor) = daelem_buff[j].getVar_biaised();
					if (!ExOp::isValid(V_param[j](coor))){
						V_param[j](coor) = M_param[j](coor) * 2.0;
						continue;
					}
					if (V_param[j](coor) < M_param[j](coor)) M_param[j](coor) = V_param[j](coor);
				}
				dacelllist[0].toMemfree();
				dacelllist[1].toMemfree();
				dacelllist[2].toMemfree();
			}//while(jobite++);
		}
		break; case 1: { // finit inital guesses for M V matrices
			double buf;
			thrbuf_totsum_mv[threadID].toZero();
			Tuple<double, 2u> guess, deriv, param;
			myHashmap<uint32_t, uint32_t> dacounts;
			double LL, LLinit;
			thrbuf_LL[threadID].toZero();
			int maxstep = 1000;
			while(thrarray.getRangeVal(jobID)){
				buf =0;
				if (doload ==false){
					if (auto ite = M_param[jobID].getIterator()) do{
						if (ite()[0] != nbcelltypes) buf += *ite;
						else {*ite = buf / nbcelltypes; buf = 0.0;}
					}while(ite++);
					if (auto ite = V_param[jobID].getIterator()) do{
						if (ite()[0] != nbcelltypes) buf += *ite;
						else {*ite = buf / nbcelltypes; buf = 0.0;}
					}while(ite++);
				}
				guess[0] = 4.0 * (((double)tr_scdata.data[jobID].getSize()) / scdata.getNBcols());
				guess[1] = 4.0 * (((double)tr_scdata.data[jobID].getSize()) / scdata.getNBcols());


				if (tr_scdata.data[jobID].getSize() == 0) lonely_MV_param[jobID].toZero();
				else{
				dacounts.toMemfree();
				dacounts[0] = scdata.getNBcols() - tr_scdata.data[jobID].getSize();
				if (auto ite = tr_scdata.data[jobID]()) do {
					dacounts[*ite]++;
				}while(ite++);
				dacurv.init();deriv.toZero();
				for(innerloop=0;innerloop<maxstep;innerloop++){

					if (auto ite = dacounts()) {
						NBmeanvar_logprob(dderbuf, ite(), guess[0], guess[1]);
						LL = dderbuf[0] * *ite; deriv[0] = dderbuf[1] * *ite; deriv[1] = dderbuf[2] * *ite;
						while(ite++){
							NBmeanvar_logprob(dderbuf, ite(), guess[0], guess[1]);
							LL += dderbuf[0] * *ite; deriv[0] += dderbuf[1] * *ite; deriv[1] += dderbuf[2] * *ite;
						}
					}
					if (innerloop == 0) thrbuf_LL[threadID][0] += LLinit = LL;
					if (dacurv.updateAscentPositiveDomain(LL,guess.mkIterator(),deriv.mkIterator())) break;
				}
				if (innerloop == maxstep) dacurv.wrFinalGuess(guess.mkIterator()); // did not converge, get the best one
				thrbuf_LL[threadID][1] += dacurv.getLastValue();
				//printf("LL%i %e -> %e in %i  %e %e\n", jobID, LLinit, dacurv.getLastValue(), innerloop,  guess[0],  guess[1]);
				lonely_MV_param[jobID] = guess;
					if (innerloop == 0) dacounts.show();
				}
				M_param_deriv[jobID] = M_param[jobID] * 0.9;
				V_param_deriv[jobID] = V_param[jobID] * 0.95;
				thrbuf_totsum_mv[threadID][0] += M_param[jobID];
				thrbuf_totsum_mv[threadID][1] += V_param[jobID];
			}
		}
		break; case 2:{ // cell-wise update f_i
			thrbuf_LL[ threadID].toZero(); thrbuf_sum_mv[threadID].toZero(); thrbuf_totsum_df[threadID].toZero();
			Tuple<double, 2u> involres;
			dasumcheck[threadID]=0;
			while(thrarray.getRangeVal(jobID)){
				if (auto site = sc_ct_state_deriv.data[jobID]()) do{
					sc_ct_state_zerobuf_M.data[jobID][site()] = thrbuf_totsum_mv[0][0](site(),dn_index[jobID]);
					sc_ct_state_zerobuf_V.data[jobID][site()] = thrbuf_totsum_mv[0][1](site(),dn_index[jobID]);
				}while(site++);
				if (auto ite = scdata.data[jobID]()) do { // loop over genes
					involres[0] = M_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]);
					involres[1] = V_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]);
					thrbuf_sum_mv[threadID] += involres;
					dasumcheck[threadID] += *ite;
					if (!ExOp::isValid(involres)){
						printf("Problems... at %i,%i  (dims %i,%i)\n", jobID, ite(), tr_scdata.data.getSize(), scdata.data.getSize());
						involres.show();
						sc_ct_state.data[jobID].show();
						M_param[ite()].selectColumn(dn_index[jobID]).show();
						V_param[ite()].selectColumn(dn_index[jobID]).show();
					}else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
				//	if ((!ExOp::isValid(dderbuf[0]))||(!ExOp::isValid(dderbuf[1]))||(!ExOp::isValid(dderbuf[2]))){
				//		printf("got that... %i %e %e -> %e %e %e\n", *ite, involres[0], involres[1], dderbuf[0],dderbuf[1],dderbuf[2]);
				//	}

					thrbuf_LL[threadID][0] += dderbuf[0];
					if (auto site = sc_ct_state_deriv.data[jobID]()) do{
						*site += M_param[ite()](site() ,dn_index[jobID]) * dderbuf[1] + V_param[ite()](site() ,dn_index[jobID]) * dderbuf[2];
						sc_ct_state_zerobuf_M.data[jobID][site()] -= M_param[ite()](site() ,dn_index[jobID]);
						sc_ct_state_zerobuf_V.data[jobID][site()] -= V_param[ite()](site() ,dn_index[jobID]);
					}while(site++);
				}while(ite++);
				if (auto site = sc_ct_state.data[jobID]()) do{
					thrbuf_totsum_df[threadID](site(),dn_index[jobID]) += *site;
				}while(site++);
			}
		}
		break; case 3:{ // gene-wise update ambiant/mix parameters M and V
			Tuple<uint32_t, 2> sccoor;
			thrbuf_LL[ threadID].toZero(); thrbuf_sum_mv[threadID].toZero(); thrbuf_totsum_mv[threadID].toZero();
			Tuple<double, 2u> involres;
			dasumcheck[threadID]=0;

			while(thrarray.getRangeVal(jobID)){
				outer[jobID] = thrbuf_totsum_df[0];
				if ((M_param_deriv[jobID].mkSum() != 0.0)) {
					printf("dirty\n");
					ExOp::show(M_param_deriv[jobID].mkSum());
					//ExOp::show(V_param_deriv[jobID].mkSum());||(V_param_deriv[jobID].mkSum() != 0.0)
					continue;
				}
				if (auto ite = tr_scdata.data[jobID]()) do { // loop over cells
					involres[0] = M_param[jobID].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]);
					involres[1] = V_param[jobID].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]);
					thrbuf_sum_mv[threadID] += involres;
					dasumcheck[threadID] += *ite;

					if (!ExOp::isValid(involres)){
						printf("Problems... at %i,%i  (dims %i,%i)\n", jobID, ite(), tr_scdata.data.getSize(), scdata.data.getSize());
						involres.show();
						sc_ct_state.data[ite()].show();
						M_param[jobID].selectColumn(dn_index[ite()]).show();
						V_param[jobID].selectColumn(dn_index[ite()]).show();
					}else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
				//	if ((!ExOp::isValid(dderbuf[0]))||(!ExOp::isValid(dderbuf[1]))||(!ExOp::isValid(dderbuf[2]))){
				//		printf("got that... %i %e %e -> %e %e %e\n", *ite, involres[0], involres[1], dderbuf[0],dderbuf[1],dderbuf[2]);
				//	}
					thrbuf_LL[threadID][0] +=  dderbuf[0];
					if (auto site = sc_ct_state.data[ite()]()) do{
						M_param_deriv[jobID](site(),dn_index[ite()]) += *site * dderbuf[1];
						V_param_deriv[jobID](site(),dn_index[ite()]) += *site * dderbuf[2];
						outer[jobID](site(),dn_index[ite()]) -= *site;
					}while(site++);
				}while(ite++);
				thrbuf_totsum_mv[threadID][0] += M_param[jobID];
				thrbuf_totsum_mv[threadID][1] += V_param[jobID];
			/*	if (jobID ==6){
					TMatrix<double> fun = outer[jobID]; fun.toZero();
					TMatrix<double> fun2 = outer[jobID]; fun2.toZero();

					for(int i =0;i<scdata.getNBcols();i++){
						if (tr_scdata.data[jobID].find(i) == 0xFFFFFFFF) { // gene has not that cell!
							if (auto site = sc_ct_state.data[i]()) do{
								fun(site(),dn_index[i]) += *site;
								fun2(site(),dn_index[i]) += *site;
							}while(site++);
						}else{
							if (auto site = sc_ct_state.data[i]()) do{
								fun2(site(),dn_index[i]) += *site;
							}while(site++);
						}
					}
					printf("check %i\n",jobID);
					outer[jobID].show();
					printf("check %i\n",jobID);
					fun.show();
					if (jobID == 6){
						printf("toto\n");
						thrbuf_totsum_df[0].show();
						fun2.show();
					}

				}*/

			}
		}
		break; case 4:{ // gene-wise update ambiant/mix parameters M and V
			while(thrarray.getRangeVal(jobID)){
				if (auto site = sc_ct_state_deriv.data[jobID]()) do{
					*site += sc_ct_state_zerobuf_M.data[jobID][site()] * zerodderbuf[1] + sc_ct_state_zerobuf_V.data[jobID][site()] * zerodderbuf[2];
				}while(site++);
			}
		}
		break; case 14:{ // gene-wise update update M_j V_j
			Tuple<uint32_t, 2> sccoor;
			thrbuf_LL[ threadID].toZero(); thrbuf_totsum_mv[threadID].toZero(); thrbuf_warn[threadID] = 0;
			Tuple<double, 2u> involres;
			Tuple<double, 2u> zeroparam;
			Vector< TMatrix<double> > tar_guess; tar_guess.setSize(2);
			Vector< TMatrix<double> > tar_deriv; tar_deriv.setSize(2);
			TMatrix<double> daouter;
			int kkk,kkkk;
			bool issafe;
			int maxstep = 5 + mainloopcount * 5 + ((mainloopcount == 0) ? 15 : 0); //first time need to guess curvature so more steps are needed
			while(thrarray.getRangeVal(jobID)){
				outer[jobID] = thrbuf_totsum_df[0];
				if (tr_scdata.data[jobID].getSize() == 0) continue; // gene has no expression... ignore
//				if (jobID < 1000) {
				tar_guess[0].toMemmove(M_param[jobID]);
				tar_guess[1].toMemmove(V_param[jobID]);
				tar_deriv[0].toMemmove(M_param_deriv[jobID]);
				tar_deriv[1].toMemmove(V_param_deriv[jobID]);
				if ((mainloopcount == 0)||(!ExOp::isValid(tar_deriv))) dacurv.init();
				else dacurv.initInvCurv(tar_deriv.mkPartitionedMetaIterator());
				daouter = thrbuf_totsum_df[0];

				if (auto ite = tr_scdata.data[jobID]()) do { // loop over cells
					daouter.selectColumn(dn_index[ite()]) -= sc_ct_state.data[ite()]();
				}while(ite++);

				for(innerloop=0;innerloop<maxstep;innerloop++){
					zeroparam[0] = tar_guess[0].vectorize().mkInnerProd(daouter.vectorize());
					zeroparam[1] = tar_guess[1].vectorize().mkInnerProd(daouter.vectorize());
					//if (jobID == 0) printf("zero%i %i %e %e\n", jobID, innerloop, zeroparam[0], zeroparam[1]);
					if ((!ExOp::isValid(zeroparam[0]))||(!ExOp::isValid(zeroparam[1]))) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else if ((zeroparam[0] < 0.0)||(zeroparam[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else NBmeanvar_logprob(dderbuf, 0, zeroparam[0], zeroparam[1]);
					// printf("dazero: %e,%e = %e\n", zeroparam[0], zeroparam[1], dderbuf[0]);
					LL = dderbuf[0];
					tar_deriv[0] = daouter * dderbuf[1];
					tar_deriv[1] = daouter * dderbuf[2];
					if ((!ExOp::isValid(dderbuf[1]))||(!ExOp::isValid(dderbuf[2]))) {
						printf("NBWIERD: %i %e %e %e %e %e\n", 0, zeroparam[0], zeroparam[1], dderbuf[0],dderbuf[1],dderbuf[2]);
					}
					if (auto ite = tr_scdata.data[jobID]()) do { // loop over cells
						involres[0] = tar_guess[0].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]());
						involres[1] = tar_guess[1].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]());
						if ((!ExOp::isValid(involres[0]))||(!ExOp::isValid(involres[1]))) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else if ((involres[0] < 0.0)||(involres[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
						LL += dderbuf[0];
						if ((!ExOp::isValid(dderbuf[1]))||(!ExOp::isValid(dderbuf[2]))) {
							printf("NBWIERD: %i %e %e %e %e %e\n", *ite, involres[0], involres[1], dderbuf[0],dderbuf[1],dderbuf[2]);
						}

						if (auto site = sc_ct_state.data[ite()]()) do{
							tar_deriv[0](site(),dn_index[ite()]) += *site * dderbuf[1];
							tar_deriv[1](site(),dn_index[ite()]) += *site * dderbuf[2];
						}while(site++);
					}while(ite++);
					//kkk=0; kkkk=0;
					//if (auto ite = tar_guess.mkPartitionedMetaIterator()) do{kkk++;	}while(ite++);
					//if (auto ite = tar_deriv.mkPartitionedMetaIterator()) do{kkkk++;}while(ite++);

					//if (jobID == 0) printf("dasize %i %i at %i LL = %e %i %i\n" ,dacurv.getSize() ,dacurv.curvature.getSize(), innerloop, LL, kkk, kkkk);
					if (innerloop == 0) thrbuf_LL[threadID][0] += LL;
					if (dacurv.updateAscentPositiveDomain(LL,tar_guess.mkPartitionedMetaIterator(),tar_deriv.mkPartitionedMetaIterator())) break;
					//if (jobID == 0) printf("hello\n");
				}
				if (((!ExOp::isValid(LL))&&(innerloop <maxstep))||(LL == 0.0)){
					printf("got %i %e critical %i m:\n", innerloop,LL, tr_scdata.data[jobID].getSize());
					tar_guess[0].show();
					printf("got critical v:\n");
					tar_guess[1].show();
					printf("got critical o:\n");
					daouter.show();
					zeroparam[0] = tar_guess[0].vectorize().mkInnerProd(daouter.vectorize());
					zeroparam[1] = tar_guess[1].vectorize().mkInnerProd(daouter.vectorize());
					//if (jobID == 0) printf("zero%i %i %e %e\n", jobID, innerloop, zeroparam[0], zeroparam[1]);
					if ((!ExOp::isValid(zeroparam[0]))||(!ExOp::isValid(zeroparam[1]))) {ExOp::show(zeroparam);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else if ((zeroparam[0] < 0.0)||(zeroparam[1] <0.0)) {ExOp::show(zeroparam);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else NBmeanvar_logprob(dderbuf, 0, zeroparam[0], zeroparam[1]);
					printf("dazero: %e,%e = %e %e %e\n", zeroparam[0], zeroparam[1], dderbuf[0], dderbuf[1], dderbuf[2]);
					LL = dderbuf[0];
					if (auto ite = tr_scdata.data[jobID]()) do { // loop over cells
						involres[0] = tar_guess[0].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]());
						involres[1] = tar_guess[1].selectColumn(dn_index[ite()]).mkInnerProd(sc_ct_state.data[ite()]());
						if ((!ExOp::isValid(involres[0]))||(!ExOp::isValid(involres[1]))) {ExOp::show(involres);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else if ((involres[0] < 0.0)||(involres[1] <0.0)) {ExOp::show(involres);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
						if (!ExOp::isValid(dderbuf)) {
							sc_ct_state.data[ite()].show();
							printf("%i:%i %e %e -> %e %e %e\n", ite(), *ite, involres[0], involres[1],dderbuf[0], dderbuf[1],dderbuf[2]);
						}
						LL += dderbuf[0];
					}while(ite++);
					printf("got LL recomputed %e:\n", LL);
				}

				//if (jobID == 0) printf("hello\n");
				if (innerloop == maxstep) dacurv.wrFinalGuess(tar_guess.mkPartitionedMetaIterator()); // did not converge, get the best one
				dacurv.wrInvCurv(tar_deriv.mkPartitionedMetaIterator());

				issafe = true;
				if (auto ite = tar_guess.mkPartitionedMetaIterator()) do{
					if (!ExOp::isValid(*ite)) issafe = false;
					else if (*ite < 0.0) issafe = false;
				}while(ite++);
				if (!issafe){
					printf("terrible %i M nan %i\n", jobID, innerloop);
					tar_guess[0].show();
					printf("terrible V nan %i\n", jobID);
					tar_guess[1].show();
					if (auto ite = tar_guess.mkPartitionedMetaIterator()) do{
						if (!ExOp::isValid(*ite)) *ite =0.0;
						else if (*ite < 0.0) *ite =fabs(*ite);
					}while(ite++);
				}

				M_param[jobID].toMemmove(tar_guess[0]);
				V_param[jobID].toMemmove(tar_guess[1]);
				M_param_deriv[jobID].toMemmove(tar_deriv[0]);
				V_param_deriv[jobID].toMemmove(tar_deriv[1]);
				//printf("test nan again %i %i %i %e\n", jobID, maxstep, innerloop, dacurv.getLastValue());
				//M_param[jobID].show();
								//printf("done%i\n",jobID);
				thrbuf_LL[threadID][1] += dacurv.getLastValue();
			//	}
				if (ExOp::isValid(M_param[jobID])) thrbuf_totsum_mv[threadID][0] += M_param[jobID];
				if (ExOp::isValid(V_param[jobID])) thrbuf_totsum_mv[threadID][1] += V_param[jobID];
			}


		}
		break; case 12:{ // cell-wise update f_i
			thrbuf_LL[ threadID].toZero(); thrbuf_totsum_df[threadID].toZero(); thrbuf_warn[threadID] = 0;
			Tuple<double, 2u> involres;
			TMatrix<double> daouter[2];
			Tuple<double, 2u> zeroparam;
			int maxstep = 5 + mainloopcount * 5 + ((mainloopcount == 0) ? 15 : 0);
			//printf("this is new\n");
			bool issafe = true;

			while(thrarray.getRangeVal(jobID)){
				issafe = true;
				if (scdata.data[jobID].getSize() == 0) continue; // cell has no umi... strange
				if (auto ite = sc_ct_state.data[jobID]()) do{
					if (!ExOp::isValid(*ite)) issafe = false;
					else if (*ite < 0.0) issafe = false;
				}while(ite++);

				//if (jobID < 1000) {
				if ((mainloopcount == 0)||(!ExOp::isValid(sc_ct_state_deriv.data[jobID]))) dacurv.init();
				else dacurv.initInvCurv(sc_ct_state_deriv.data[jobID]());
				daouter[0] = thrbuf_totsum_mv[0][0];
				daouter[1] = thrbuf_totsum_mv[0][1];
				if (auto site = sc_ct_state.data[jobID]()) do{ // loop over genes
					if (tr_scdata.data[site()].getSize() == 0) continue; // ignore genes with no expression
					if (site() >= M_param.getSize()) {printf("oh no!\n"); continue;}
					daouter[0] -= M_param[site()];
					daouter[1] -= V_param[site()];
				}while(site++);
				for (innerloop=0;innerloop<maxstep;innerloop++){
					zeroparam[0] = daouter[0].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					zeroparam[1] = daouter[1].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					//if (jobID == 0) printf("zero%i %i %e %e\n", jobID, innerloop, zeroparam[0], zeroparam[1]);
					if ((!ExOp::isValid(zeroparam[0]))||(!ExOp::isValid(zeroparam[1]))) {ExOp::show(zeroparam);ExOp::toZero(dderbuf); thrbuf_warn[threadID]++;}
					else if ((zeroparam[0] < 0.0)||(zeroparam[1] <0.0)) {ExOp::show(zeroparam);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else NBmeanvar_logprob(dderbuf, 0, zeroparam[0], zeroparam[1]);
					LL = dderbuf[0];
					if (auto site = sc_ct_state_deriv.data[jobID]()) do{
						*site = daouter[0](site(),dn_index[jobID]) * dderbuf[1] + daouter[1](site(),dn_index[jobID]) * dderbuf[2];
					}while(site++);

					if (auto ite = scdata.data[jobID]()) do { // loop over genes
						involres[0] = M_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						involres[1] = V_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						if ((!ExOp::isValid(involres[0]))||(!ExOp::isValid(involres[1]))) {ExOp::show(involres);ExOp::toZero(dderbuf); thrbuf_warn[threadID]++;}
						else if ((involres[0] < 0.0)||(involres[1] <0.0)) {ExOp::show(involres);ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
						LL += dderbuf[0];

						if (auto site = sc_ct_state_deriv.data[jobID]()) do{
							*site += M_param[ite()](site() ,dn_index[jobID]) * dderbuf[1] + V_param[ite()](site() ,dn_index[jobID]) * dderbuf[2];
						}while(site++);

					}while(ite++);
					//if (jobID == 0) printf("LL %e size %i %i\n", LL, dacurv.getSize(), dacurv.curvature.getSize());
					if (innerloop == 0) thrbuf_LL[threadID][0] += LL;
					if (dacurv.updateAscentPositiveDomain(LL,sc_ct_state.data[jobID](),sc_ct_state_deriv.data[jobID]())) break;
					//if (jobID == 0) printf("hello\n");
				}
				if (((!ExOp::isValid(LL))&&(innerloop <maxstep))||(LL == 0.0)){
					printf("got cr %eitical m %i:\n", LL, scdata.data[jobID].getSize());
					daouter[0].selectColumn(dn_index[jobID]).show();
					printf("got critical %i v:\n",dn_index[jobID]);
					daouter[1].selectColumn(dn_index[jobID]).show();
					printf("got critical state:\n");
					sc_ct_state.data[jobID].show();
					zeroparam[0] = daouter[0].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					zeroparam[1] = daouter[1].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					//if (jobID == 0) printf("zero%i %i %e %e\n", jobID, innerloop, zeroparam[0], zeroparam[1]);
					if ((!ExOp::isValid(zeroparam[0]))||(!ExOp::isValid(zeroparam[1]))) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else if ((zeroparam[0] < 0.0)||(zeroparam[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else NBmeanvar_logprob(dderbuf, 0, zeroparam[0], zeroparam[1]);
					LL = dderbuf[0];
					printf("dazero: %e,%e = %e %e %e\n", zeroparam[0], zeroparam[1], dderbuf[0], dderbuf[1], dderbuf[2]);
					if (auto ite = scdata.data[jobID]()) do { // loop over cells
						involres[0] = M_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						involres[1] = V_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						if ((!ExOp::isValid(involres[0]))||(!ExOp::isValid(involres[1]))) {ExOp::toZero(dderbuf); thrbuf_warn[threadID]++;}
						else if ((involres[0] < 0.0)||(involres[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
						if (!ExOp::isValid(dderbuf)) {
							printf("%i:%i %e %e -> %e %e %e\n", ite(), *ite, involres[0], involres[1],dderbuf[0], dderbuf[1],dderbuf[2]);
							M_param[ite()].selectColumn(dn_index[jobID]).show();
							V_param[ite()].selectColumn(dn_index[jobID]).show();
						}
						LL += dderbuf[0];
					}while(ite++);
					printf("%i recompute %i:%e\n", M_param.getSize(),jobID, LL);
				}
				//if (jobID == 0) printf("hello\n");
				if (innerloop == maxstep) dacurv.wrFinalGuess(sc_ct_state.data[jobID]()); // did not converge, get the best one
				dacurv.wrInvCurv(sc_ct_state_deriv.data[jobID]());
				thrbuf_LL[threadID][1] += dacurv.getLastValue();
				issafe = true;
				if (auto ite = sc_ct_state.data[jobID]()) do{
					if (!ExOp::isValid(*ite)) issafe = false;
					else if (*ite < 0.0) issafe = false;
				}while(ite++);
				if (!issafe) {
					printf("terrible %i state nan %i:", innerloop, jobID);
					printf("got cr %eitical m:\n", LL);
					daouter[0].selectColumn(dn_index[jobID]).show();
					printf("got critical %i v:\n",dn_index[jobID]);
					daouter[1].selectColumn(dn_index[jobID]).show();
					printf("got critical state:\n");
					sc_ct_state.data[jobID].show();
					zeroparam[0] = daouter[0].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					zeroparam[1] = daouter[1].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
					//if (jobID == 0) printf("zero%i %i %e %e\n", jobID, innerloop, zeroparam[0], zeroparam[1]);
					if ((!ExOp::isValid(zeroparam[0]))||(!ExOp::isValid(zeroparam[1]))) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else if ((zeroparam[0] < 0.0)||(zeroparam[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
					else NBmeanvar_logprob(dderbuf, 0, zeroparam[0], zeroparam[1]);
					LL = dderbuf[0];
					if (auto site = sc_ct_state_deriv.data[jobID]()) do{
						*site = daouter[0](site(),dn_index[jobID]) * dderbuf[1] + daouter[1](site(),dn_index[jobID]) * dderbuf[2];
					}while(site++);

					printf("dazero: %e,%e = %e\n", zeroparam[0], zeroparam[1], dderbuf[0]);
					printf("zerozero: %i %e %e %e %e %e\n", 0, zeroparam[0], zeroparam[1], dderbuf[0],dderbuf[1],dderbuf[2]);
					printf("tmpderiv:"); sc_ct_state_deriv.data[jobID].show();
					if (auto ite = scdata.data[jobID]()) do { // loop over cells
						involres[0] = M_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						involres[1] = V_param[ite()].selectColumn(dn_index[jobID]).mkInnerProd(sc_ct_state.data[jobID]());
						if ((!ExOp::isValid(involres[0]))||(!ExOp::isValid(involres[1]))) {ExOp::toZero(dderbuf); thrbuf_warn[threadID]++;}
						else if ((involres[0] < 0.0)||(involres[1] <0.0)) {ExOp::toZero(dderbuf);thrbuf_warn[threadID]++;}
						else NBmeanvar_logprob(dderbuf, *ite, involres[0], involres[1]);
						if (!ExOp::isValid(dderbuf)) {
							printf("%i:%i %e %e -> %e %e %e\n", ite(), *ite, involres[0], involres[1],dderbuf[0], dderbuf[1],dderbuf[2]);
							M_param[ite()].selectColumn(dn_index[jobID]).show();
							V_param[ite()].selectColumn(dn_index[jobID]).show();
						}
						printf("yesyesy: %i %e %e %e %e %e\n", *ite, involres[0], involres[1], dderbuf[0],dderbuf[1],dderbuf[2]);
						LL += dderbuf[0];
						printf("tmpderiv:"); sc_ct_state_deriv.data[jobID].show();
					}while(ite++);
					printf("%i recompute %i:%e\n", M_param.getSize(),jobID, LL);


					if (auto ite = sc_ct_state.data[jobID]()) do{
						printf("%i:%e\t", ite(),*ite);
						if (!ExOp::isValid(*ite)) *ite = 0.0;
						else if (*ite < 0.0) *ite = 0.0;
					}while(ite++);
					printf(" and fix\n");


				}
				//}
				issafe = true;
				if (auto ite = sc_ct_state.data[jobID]()) do{
					if (!ExOp::isValid(*ite)) issafe = false;
				}while(ite++);

				if (issafe) thrbuf_totsum_df[threadID].selectColumn(dn_index[jobID]) += sc_ct_state.data[jobID]();
			}
		}
		break;	case 5:{
			while(thrarray.getRangeVal(jobID)){
				//if (jobID == 6){printf("deriv update with %e\n",zerodderbuf[1]);M_param_deriv[jobID].show();}
				M_param_deriv[jobID] += outer[jobID] * zerodderbuf[1];
				V_param_deriv[jobID] += outer[jobID] * zerodderbuf[2];
				//if (jobID == 6){printf("gives\n");M_param_deriv[jobID].show();}
			}
		}
		break; case 6:{ // compute sum of states/donors
			thrbuf_totsum_df[threadID].toZero();
			while(thrarray.getRangeVal(jobID)){
				if (auto site = sc_ct_state.data[jobID]()) do{
					thrbuf_totsum_df[threadID](site(),dn_index[jobID]) += *site;
				}while(site++);
			}
		}
		break; case 7:{ // compute sum of MV params
			thrbuf_totsum_mv[threadID].toZero();
			while(thrarray.getRangeVal(jobID)){
				thrbuf_totsum_mv[threadID][0] += M_param[jobID];
				thrbuf_totsum_mv[threadID][1] += V_param[jobID];
			}
		}
		break; case 8:{
			Tuple<double> curg; curg.setSize(nbcelltypes + nbgenotypes);
			Tuple<double> derv; derv.setSize(nbcelltypes + nbgenotypes);

			Tuple<uint32_t, 2> sccoor;
			Tuple<double> ctselect, dnselect;
			Tuple<double> dLL_ctselect, dLL_dnselect;
			Tuple<double> derbuf[2];
			Tuple<double> projection[2][2];


			int i,j;
			ctselect.setSize(nbcelltypes + 1); ctselect[nbcelltypes] =0.0; // do not use ambiant RNA
			dnselect.setSize(nbgenotypes);
			derbuf[0].setSize(nbcelltypes + 1);
			derbuf[1].setSize(nbgenotypes);
			MemoryLimitedOptimizationScope adllocal;
			uint32_t dac;

			//if (auto jobite = thrbase.getRangeIterator(threadID, bkdata.getNBcols(),true)) do{
			while(thrarray.getRangeVal(jobID)){
				sccoor[1] = jobID;//ite();

				curg[0] = sqrt(sqrt(bk_scale[sccoor[1]]) / nbcelltypes);
				for(i=1;i<nbcelltypes;i++) curg[i] = curg[0];
				curg[nbcelltypes] = sqrt(sqrt(bk_scale[sccoor[1]]) / nbgenotypes);
				for(i++;i<curg.getSize();i++) curg[i] = curg[nbcelltypes];
				double LL, initLL, initLL2;
				dacurv.init(nbcelltypes + nbgenotypes,0.0000001);
				int minstep = 0;
				for(j =0; j<1000;j++){
					for(i=0;i<nbcelltypes;i++) ctselect[i] = curg[i] * curg[i];
					for(;i<curg.getSize();i++) dnselect[i-nbcelltypes] = curg[i] * curg[i];
					LL = 0; derbuf[0].toZero(); derbuf[1].toZero();
					for(sccoor[0] = 0; sccoor[0]< nbgenes ; sccoor[0]++) {
						if (genefilter[sccoor[0]] == 0) continue;
						projection[0][0] =  M_param[sccoor[0]] * dnselect;
						projection[0][1] =  V_param[sccoor[0]] * dnselect;
						projection[1][0] =  M_param[sccoor[0]].mkBackMult(ctselect);
						projection[1][1] =  V_param[sccoor[0]].mkBackMult(ctselect);
						dac = bkdata.data[sccoor[1]].find(sccoor[0]);
						dac = (dac == 0xFFFFFFFF) ? 0 : bkdata.data[sccoor[1]].deref(dac);
						NBmeanvar_logprob(dderbuf, dac, projection[0][0].mkInnerProd(ctselect), projection[0][1].mkInnerProd(ctselect));
						LL += dderbuf[0];
						derbuf[0] += projection[0][0] *  dderbuf[1] + projection[0][1] *  dderbuf[2];
						derbuf[1] += projection[1][0] *  dderbuf[1] + projection[1][1] *  dderbuf[2];
					}
					for(i=0;i<nbcelltypes;i++) derv[i] = derbuf[0][i] * curg[i] * 2.0;
					for(i++;i<curg.getSize();i++) derv[i] = derbuf[1][i-nbcelltypes] * curg[i]* 2.0;
					//printf("ll[%i]: %e\n", j, LL);
					//printf("deriv:"); ExOp::show(derv);
					if (j < 2){ if (j== 0) initLL = LL; else initLL2 = LL;}
					else{
						if (!ExOp::isValid(initLL2)) printf("LL_%i[%i]->%e\n",jobID,j,LL);
						if ((ExOp::isValid(LL))&&(LL >= dacurv.getLastValue())) minstep++;
					}
					if ((dacurv.updateAscent(LL,curg(),derv()))&&(minstep > 25))   break;

				}
				printf("col%i: LL %e,%e -> %e in %i steps\n", jobID, initLL,initLL2, LL,j);
				curg.show();
				Tuple<uint32_t,2> coor; coor[1] = jobID;
				for(coor[0]=0;coor[0]<nbcelltypes;coor[0]++) ct_deconv_results(coor) = curg[coor[0]] * curg[coor[0]];
				for(coor[0]=0;coor[0]<nbgenotypes;coor[0]++) dn_deconv_results(coor) = curg[coor[0]+nbcelltypes] * curg[coor[0] +nbcelltypes];

			/*	if ((!ExOp::isValid(LL))||(LL < initLL)){
					printf("%i failed... stragery no2.!\n", jobID);
				ctselect[0] = sqrt(bk_scale[sccoor[1]]) / nbcelltypes;
				for(i=1;i<nbcelltypes;i++) ctselect[i] = ctselect[0];
				dnselect[0] = sqrt(bk_scale[sccoor[1]]) / nbgenotypes;
				for(i=1;i<nbgenotypes;i++) dnselect[i] = dnselect[0];


				dLL_ctselect = ctselect * 0.9;
				dLL_dnselect = dnselect * 0.9;

				adllocal.initVariable("Cvec", dLL_ctselect.mkIterator() ,0).initVariable("Dvec",dLL_dnselect.mkIterator() ,0).initLearn();

				double LL, initLL, initLL2, datmp, datmp2;
				int minstep = 0;
				bool daexitflag=false;
				for(j =0; j<2500;j++){
					LL = 0.0;
					for(sccoor[0] = 0; sccoor[0]< nbgenes ; sccoor[0]++) {
						if (genefilter[sccoor[0]] == 0) continue;
						projection[0][0] =  M_param[sccoor[0]] * dnselect;
						projection[0][1] =  V_param[sccoor[0]] * dnselect;
						projection[1][0] =  M_param[sccoor[0]].mkBackMult(ctselect);
						projection[1][1] =  V_param[sccoor[0]].mkBackMult(ctselect);
						dac = bkdata.data[sccoor[1]].find(sccoor[0]);
						dac = (dac == 0xFFFFFFFF) ? 0 : bkdata.data[sccoor[1]].deref(dac);

						NBmeanvar_logprob(dderbuf, dac, projection[0][0].mkInnerProd(ctselect), projection[0][1].mkInnerProd(ctselect));
						LL += dderbuf[0];
						dLL_ctselect += projection[0][0] *  dderbuf[1] + projection[0][1] *  dderbuf[2];
						dLL_dnselect += projection[1][0] *  dderbuf[1] + projection[1][1] *  dderbuf[2];
					}

					if (j < 2){ if (j== 0) initLL = LL; else initLL2 = LL;}
					else{
						if (!ExOp::isValid(initLL2)) printf("LL_%i[%i]->%e\n",jobID,j,LL);
					}
					dLL_ctselect[nbcelltypes] = 0.0; // ignore that parameter...
					if (jobID == 0) {
						printf("job0 start%i\n", j);
						dLL_ctselect.show();dLL_dnselect.show();
						ctselect.show();dnselect.show();
						printf("job0 end%i\n", j);
					}
					if (adllocal.doesAccept(LL)){
						if ((daexitflag)&&(j>50)) break;
						if (jobID == 0) dnselect.show();
						daexitflag = (datmp = adllocal.updatePositiveVariable("Cvec",ctselect.mkIterator(), dLL_ctselect.mkIterator())) < 0.0001;
						daexitflag &=(datmp2 = adllocal.updatePositiveVariable("Dvec",dnselect.mkIterator(), dLL_dnselect.mkIterator())) < 0.0001;
						if (jobID == 0) dnselect.show();
						//printf("%i_%i: %e %e\n", jobID, j, datmp ,datmp2);
					}else{
						//printf("LL[%i]: %e (Rejected!)\n",adl.step, thrbuf_LL[0]);
						adllocal.backtrackVariable("Cvec",ctselect.mkIterator());
						adllocal.backtrackVariable("Dvec",dnselect.mkIterator());
						dLL_ctselect.toZero();
						dLL_dnselect.toZero();
					}
				}

				printf("col%i: LL %e,%e -> %e in %i steps  (%e is zero? :D)\n", jobID, initLL,initLL2, LL,j, ctselect[nbcelltypes]);
				Tuple<uint32_t,2> coor; coor[1] = jobID;
				for(coor[0]=0;coor[0]<nbcelltypes;coor[0]++) ct_deconv_results(coor) = ctselect[coor[0]];
				for(coor[0]=0;coor[0]<nbgenotypes;coor[0]++) dn_deconv_results(coor) = dnselect[coor[0]];
				}*/
			}

		}
		break; case 9:{

			Tuple<double> furg, dfurg; furg.setSize((nbcelltypes+1) * nbgenotypes);
			dfurg.setSize(furg.getSize());
			Tuple<double> curg; curg.setSize(nbcelltypes + nbgenotypes + (useinter ? 1 : 0));
			Tuple<double> derv; derv.setSize(curg.getSize());

			Tuple<uint32_t, 2> sccoor;
			Tuple<double> ctselect, dnselect;
			Tuple<double> dLL_ctselect, dLL_dnselect;
			Tuple<double> derbuf[2];
			Tuple<double> projection[2][2];


			int i,j,k;
			ctselect.setSize(nbcelltypes + 1); ctselect[nbcelltypes] =0.0; // do not use ambiant RNA
			dnselect.setSize(nbgenotypes);
			derbuf[0].setSize(nbcelltypes + 1);
			derbuf[1].setSize(nbgenotypes);
			double dadebg[2];
			uint32_t dac, tdac;
			useinter = false;
			SparseMatrix<double> dadbgmast; dadbgmast.setNBcols(15);
			//if (auto jobite = thrbase.getRangeIterator(threadID, bkdata.getNBcols(),true)) do{
			while(thrarray.getRangeVal(jobID)){
				sccoor[1] = jobID;//ite();
				if (jobID != 0) continue;
				curg[0] = sqrt(bk_scale[sccoor[1]]) / nbcelltypes;
				double LL, initLL, initLL2;

				LL = bk_scale[sccoor[1]] / (nbcelltypes * nbgenotypes);
				for(j=0,k=0;j<nbgenotypes;j++){
					for(i=0;i<nbcelltypes;i++) furg[k++] = LL;
					furg[k++] = 0;
				}

				//if (jobID == 0) dacurv.init( curg.getSize(), 0.0001);
				//else i
				dacurv.init(pow(0.5, 15.0));
				int minstep = 0;

				// unconstrained search!

				for(j=0;j<1000;j++){
					LL = 0.0; dfurg.toZero();
					for(sccoor[0] = 0; sccoor[0]< nbgenes ; sccoor[0]++) {
						if (genefilter[sccoor[0]] == 0) continue;

						dac = bkdata.data[sccoor[1]].find(sccoor[0]);
						dac = ((dac == 0xFFFFFFFF) ? 0 : bkdata.data[sccoor[1]].deref(dac));
						dadebg[0] = furg.mkInnerProd(M_param[sccoor[0]].vectorize().mkIterator());
						dadebg[1] = furg.mkInnerProd(V_param[sccoor[0]].vectorize().mkIterator());
						NBmeanvar_logprob(dderbuf, dac, dadebg[0], dadebg[1]);
						for(k=0;k<dfurg.getSize();k++) dfurg[k] += M_param[sccoor[0]].data[k] * dderbuf[1] + V_param[sccoor[0]].data[k] * dderbuf[2];
						LL += dderbuf[0];

					}
					if (j== 0) initLL = LL;
					printf("LL %e\n", LL); furg.show();
					for(k=0;k<nbgenotypes;k++) dfurg[((k+1) *nbcelltypes)-1] = 0.0;
					if (dacurv.checkDerivative(LL, furg.mkModuloRangeIterator(nbcelltypes+1,0,nbcelltypes), dfurg.mkModuloRangeIterator(nbcelltypes+1,0,nbcelltypes))) break;
					if (j == 100) break;
					//if (dacurv.updateAscentPositiveDomain(LL,furg.mkModuloRangeIterator(nbcelltypes+1,0,nbcelltypes),dfurg.mkModuloRangeIterator(nbcelltypes+1,0,nbcelltypes),&dadbgmast,j)) break;
				}
				dadbgmast.show();
				printf("unconstr %i in %i steps: %e -> %e\n", dacurv.getSize(),j, initLL, LL);
				furg.show();
				// use unconstrained solution as starting guess...
				int nbstep_unconstr = j;
				LL = 0;curg.toZero();
				for(j=0,k=0;j<nbgenotypes;j++){
					for(i=0;i<nbcelltypes;i++){
						curg[j+nbcelltypes] += furg[k];
						curg[i] += furg[k];
						LL += furg[k++];
					}
					k++;
				}
				curg /= sqrt(LL);

				//if (!ExOp::isValid(curg[0])) curg[0] = 1.0;
				//for(i=1;i<nbcelltypes;i++) curg[i] = curg[0];
				//curg[nbcelltypes] = sqrt(bk_scale[sccoor[1]]) / nbgenotypes;
				//if (!ExOp::isValid(curg[nbcelltypes])) curg[nbcelltypes] = 1.0;
				//for(i=1;i<nbgenotypes;i++) curg[i+nbcelltypes] = curg[nbcelltypes];
				if (useinter) curg.last() = 0.01 * curg[0] * curg[nbcelltypes];

				dacurv.init();
				for(j =0; j<1000;j++){
					for(i=0;i<nbcelltypes;i++) ctselect[i] = curg[i];
					for(i=0;i<nbgenotypes;i++) dnselect[i] = curg[i+ nbcelltypes];
					LL = 0.0; derbuf[0].toZero(); derbuf[1].toZero();
					if (useinter) derv.last() = 0.0;
					tdac =0;
					for(sccoor[0] = 0; sccoor[0]< nbgenes ; sccoor[0]++) {
						if (genefilter[sccoor[0]] == 0) continue;

						projection[0][0] =  M_param[sccoor[0]] * dnselect;
						projection[0][1] =  V_param[sccoor[0]] * dnselect;
						projection[1][0] =  M_param[sccoor[0]].mkBackMult(ctselect);
						projection[1][1] =  V_param[sccoor[0]].mkBackMult(ctselect);
						dac = bkdata.data[sccoor[1]].find(sccoor[0]);
						//printf("finding %i in %i give %i\n", sccoor[0] , sccoor[1], dac);
						dac = ((dac == 0xFFFFFFFF) ? 0 : bkdata.data[sccoor[1]].deref(dac));
						//printf("and hence %i\n",dac);

						tdac += dac;
						//if (useinter) NBmeanvar_logprob(dderbuf, dac, projection[0][0].mkInnerProd(ctselect) + curg.last() * lonely_MV_param[sccoor[0]][0] , projection[0][1].mkInnerProd(ctselect)+ curg.last() * lonely_MV_param[sccoor[0]][1]);
						//else i
						//printf("%i %i param %e %e\n", j, dac, projection[0][0].mkInnerProd(ctselect), projection[0][1].mkInnerProd(ctselect));
						//printf("%i %i param %e %e\n", j, dac, projection[1][0].mkInnerProd(dnselect), projection[1][1].mkInnerProd(dnselect));
						/*curg.show();
						printf("ct:");ctselect.show();
						printf("dn:");dnselect.show();

						printf("00:");projection[0][0].show();
						printf("01:");projection[0][1].show();
						printf("10:");projection[1][0].show();
						printf("11:");projection[1][1].show();*/
						NBmeanvar_logprob(dderbuf, dac, projection[0][0].mkInnerProd(ctselect), projection[0][1].mkInnerProd(ctselect));

						LL += dderbuf[0];
						//ExOp::show(dderbuf);

						if (auto ite = derbuf[0].getIterator()) do{
							*ite += projection[0][0][ite()] * dderbuf[1] + projection[0][1][ite()] * dderbuf[2];
						}while(ite++);
						if (auto ite = derbuf[1].getIterator()) do{
							*ite += projection[1][0][ite()] * dderbuf[1] + projection[1][1][ite()] * dderbuf[2];
						}while(ite++);
						//derbuf[0] = projection[0][0] * dderbuf[1];// + projection[0][1] * dderbuf[2];
						//derbuf[1] = projection[1][0] * dderbuf[1];// + projection[1][1] * dderbuf[2];
						/*ExOp::show(projection[0][0] * dderbuf[1]);
						ExOp::show(projection[1][0] * dderbuf[1]);
						ExOp::show(derbuf[0]);
						ExOp::show(derbuf[1]);*/
						if (useinter) derv.last() += lonely_MV_param[sccoor[0]][0] * dderbuf[1] + lonely_MV_param[sccoor[0]][1] * dderbuf[2];
					}
					for(i=0;i<nbcelltypes;i++) derv[i] = derbuf[0][i];
					for(i=0;i< nbgenotypes;i++) derv[i+nbcelltypes] = derbuf[1][i];
					//printf("ll[%i]: %e\n", j, LL);
					//printf("deriv %i %i:", nbcelltypes, nbgenotypes);
					//ExOp::show(derv);
					//printf("LL %e\n", LL);
					if (j== 0) initLL = LL;
					if (j== 100) {printf("startder%i: ", tdac);  derv.show();}
					//if (jobID == 0) dacurv.checkDerivative(LL,curg(),derv());
					if (dacurv.updateAscentPositiveDomain(LL,curg.getIterator(),derv.getIterator())) break;
					//printf("dadelta = %e\n", dacurv.curvature[0]);
				}
				if (j == 1000) dacurv.wrFinalGuess(curg.getIterator());
				if (j == 0){
					printf("mmm, NAN on first try...\n");
					curg.show();
					/*printf("and then\n");
					for(sccoor[0] = 0; sccoor[0]< nbgenes ; sccoor[0]++) {
						if (genefilter[sccoor[0]] == 0) continue;
						projection[0][0] =  M_param[sccoor[0]] * dnselect;
						projection[0][1] =  V_param[sccoor[0]] * dnselect;
						projection[1][0] =  M_param[sccoor[0]].mkBackMult(ctselect);
						projection[1][1] =  V_param[sccoor[0]].mkBackMult(ctselect);
						dac = bkdata.data[sccoor[1]].find(sccoor[0]);
						dac = (dac == 0xFFFFFFFF) ? 0 : bkdata.data[sccoor[1]].deref(dac);
						NBmeanvar_logprob(dderbuf, dac, projection[0][0].mkInnerProd(ctselect) + curg.last() * lonely_MV_param[sccoor[0]][0] , projection[0][1].mkInnerProd(ctselect)+ curg.last() * lonely_MV_param[sccoor[0]][1]);
						if (!ExOp::isValid(dderbuf[0])){
							M_param[sccoor[0]].show();
							V_param[sccoor[0]].show();
							printf("%i: %i with %e\t%e LL=%e\n", sccoor[0] , dac, projection[0][0].mkInnerProd(ctselect), projection[0][1].mkInnerProd(ctselect), dderbuf[0]);
						}
					}*/

				}
				LLoutput_buffer[jobID] = dacurv.getLastValue();
				printf("col%i: LL %e -> %e in %i steps %e\n", jobID, initLL, dacurv.getLastValue(),j, curg.last());
				curg.show();
				Tuple<uint32_t,2> coor; coor[1] = jobID;
				double dasum = 0;
				for(coor[0]=0;coor[0]<nbcelltypes;coor[0]++) dasum += curg[coor[0]];
				for(coor[0]=0;coor[0]<nbcelltypes;coor[0]++) ct_deconv_results(coor) = curg[coor[0]] / dasum;
				dasum = 0;
				for(coor[0]=0;coor[0]<nbgenotypes;coor[0]++) dasum += curg[coor[0]+nbcelltypes];
				for(coor[0]=0;coor[0]<nbgenotypes;coor[0]++) dn_deconv_results(coor) = curg[coor[0]+nbcelltypes] / dasum;
				intercept[jobID] = curg.last();
			}
		}
		break; case 20:{ // prepare GP inference ( needs thrbuf_totsum_df[0] defined
			int i;
			thrbuf_gpstate_dx[threadID].setSizes(2,nbgenotypes).toRand();
			thrbuf_LL[ threadID].toZero();
			while(thrarray.getRangeVal(jobID)){
                if (genefilter[jobID] == 0) continue;
                M_serror[jobID].setSizes(nbcelltypes+1,nbgenotypes);
                if (auto ite = thrbuf_totsum_df[0].getIterator()) do{
                    M_serror[jobID](ite()) = (M_param[jobID](ite()) + V_param[jobID](ite())) / *ite;
                    if (!ExOp::isValid(M_serror[jobID](ite()))) M_serror[jobID](ite()) = 1000000.0 * (M_param[jobID](ite()) + V_param[jobID](ite()));
                }while( ite++);
				for(int curct =0; curct < nbcelltypes;curct++) gp_mv[jobID][curct].init();
            }
		}
		break; case 21:{ // fit MV-GP, for every genes
			thrbuf_LL[ threadID].toZero();
			Tuple<double> curderiv;
			int j;
			printf("lets go! %i\n", threadID);
			while(thrarray.getRangeVal(jobID)){
				if (genefilter[jobID] == 0) continue;
				if (gp_mv.find(jobID) == 0xFFFFFFFF) printf("did not find!!!\n");
				gp_mv[jobID][0].show();
				if (M_param[jobID].selectRow(0).getShape()[0] != nbgenotypes){
                    printf("not the right size! %i vs %i\n", nbgenotypes , M_param[jobID].selectRow(0).getShape()[0]);
                    exit(1);
				}
				for(int curct=0; curct < nbcelltypes;curct++){
                    dacurv.init();
					//if ((mainloopcount == 0)||(!ExOp::isValid(gp_mv_dx[jobID][curct].getIterator())))
					//else dacurv.initInvCurv(gp_mv_dx[jobID][curct].getIterator());
					for(j =0; j<1000;j++){
						LL = gp_mv[jobID][curct].wrLLDerivative2(curderiv, M_param[jobID].selectRow(curct), gpstate.vectorize(),M_serror[jobID].selectRow(curct), threadID == 0);
						if (j == 0) thrbuf_LL[threadID][0] += LL;
						if (dacurv.updateAscent(LL,gp_mv[jobID][curct].mkParamIterator(),curderiv.getIterator())) break;
						gp_mv[jobID][curct].checkDomain(dacurv);
						if (threadID == 0) gp_mv[jobID][curct].show();
					}
					printf("done in %i steps\n", j);
					if (j == 1000) dacurv.wrFinalGuess(gp_mv[jobID][curct].mkParamIterator());
					//dacurv.wrInvCurv(gp_mv_dx[jobID][curct].getIterator());
					thrbuf_LL[threadID][1] += dacurv.getLastValue();
				}
			}
		}
		break; case 22:{ // fit donor-states, cummul for every genes
			thrbuf_gpstate_dx[threadID].toZero();
			thrbuf_LL[ threadID].toZero();
			thrbuf_gpstate_dx[threadID].toZero(); thrbuf_LL[ threadID].toZero();
			Tuple<double> daderiv;
			while(thrarray.getRangeVal(jobID)){
                if (genefilter[jobID] == 0) continue;
				for(int curct=0; curct < nbcelltypes;curct++){
					thrbuf_LL[ threadID][0] += gp_mv[jobID][curct].wrLLTimeDerivative(daderiv, M_param[jobID].selectRow(curct), gpstate.vectorize(), M_serror[jobID].selectRow(curct));
					thrbuf_gpstate_dx[threadID].vectorize() += daderiv.mkIterator();
				}
			}
		}
		break; case 23:{ // Define mean and std error in the input for the GP model, use prior count of a fictive 0.01 expression cell is all cases.
            Tuple<unsigned int, 4u> coor4;
            Tuple<unsigned int, 2> coor2;
            double kappa, one_over_theta;
            while(thrarray.getRangeVal(jobID)){
                uint32_t dag = genemap.deref(jobID);
                coor4[2] = jobID;
                for(coor2[0] =0; coor2[0] < nbcelltypes;coor2[0]++){
                    coor4[3] = coor2[0];
                    for(coor4[1]=0; coor4[1] < nbgenotypes;coor4[1]++){
                        coor2[1] = coor4[1];
                        kappa = M_param[dag](coor2) * thrbuf_totsum_df[0](coor2) + 1.0; // 1 / 0.01 value for prior count
                        one_over_theta = thrbuf_totsum_df[0](coor2) + 1.0;
                        coor4[0] = 0; dagpinput(coor4) = d_lngamma_dx(kappa) - log(one_over_theta);
                        if (!ExOp::isValid(dagpinput(coor4))) dagpinput(coor4) =0.0;
                        coor4[0] = 1; dagpinput(coor4) = d2_lngamma_dx2(kappa);
                        if (!ExOp::isValid(dagpinput(coor4))) dagpinput(coor4) =100000.0;
                    }
                }
            }

		}break; default: printf("unknown task %i for %i\n", curtask, threadID);
		}
		//printf("scope happy! %i\n", threadID); fflush(stdout);
	return 0;}
};	Task ltk;

	ltk.scdata.rdRcppdgCMatrix(scdata_in);
	ltk.bkdata.rdRcppdgCMatrix(bulk_in);

	ltk.ct_index = Rcpp::as<std::vector<uint32_t> >(inlist[0]);
	ltk.dn_index = Rcpp::as<std::vector<uint32_t> >(inlist[1]);
	ltk.ct_names = Rcpp::as<Rcpp::CharacterVector >(inlist[2]);
	ltk.dn_names = Rcpp::as<Rcpp::CharacterVector >(inlist[3]);
	ltk.gene_names = Rcpp::as<Rcpp::CharacterVector >(inlist[4]);
	int nbsteps =  Rcpp::as<int >(inlist[6]);
	ltk.genefilter = Rcpp::as<vector<uint32_t> >(inlist[5]);
	string dapath = Rcpp::as<string >(inlist[7]);
	int flag = Rcpp::as<int >(inlist[8]);
	ltk.doload = (flag & 1) != 0;
	ltk.useold = (flag & 2) != 0;
	ltk.useinter = (flag & 4) != 0;
	ltk.useold2 = (flag & 8) != 0;
	ltk.usegp = (flag & 16) != 0;
	Rcpp::List davec = bulk_in.slot("Dimnames");
	ltk.bulknames = davec[1]; // Rcpp::as<Rcpp::CharacterVector>(Rcpp::colnames(bulk_in));
return ltk.run(nbsteps, dapath.c_str());}

