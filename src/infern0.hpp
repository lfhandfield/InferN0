#include "require11.h"
#include "bastructs.h"
#include "primitive.h"
#include "Advanced.h"

using namespace LFHPrimitive;

enum INFRN0_enum{
	INFRN0_CELL_CLUSTER=0,
	INFRN0_GENE_CLUSTER=1,
	INFRN0_CELL_SCALE=2,
	INFRN0_GENE_SCALE=3,
	INFRN0_CELL_STATE=4,
	INFRN0_GENE_STATE=5,
	INFRN0_RAW_DATA=6,
	INFRN0_DATA=7,
	INFRN0_CELL_ORDER=8,
	INFRN0_GENE_ORDER=9,
	INFRN0_CELL_NAMES=10,
	INFRN0_GENE_NAMES=11,
	INFRN0_CELL_STATE_MEAN_DEVIATION=12,
	INFRN0_GENE_STATE_MEAN_DEVIATION=13,
	INFRN0_GENE_STATE_MATRICES=14,
	INFRN0_CELL_COVERAGE=15, // aka 1 - dropout rate
	INFRN0_GENE_COVERAGE=16,
	INFRN0_PVAL_DATA=17,	
	INFRN0_BIODROP_PROB=18,
	INFRN0_CELL_COLOR=19,
	INFRN0_GENE_COLOR=20,
	INFRN0_CELL_SIZE=21,
	INFRN0_GENE_SIZE=22,
	INFRN0_CELL_TRANSITION=23,
	INFRN0_GENE_TRANSITION=24

};

enum INFRN0_AXES_enum{
	INFRN0_AXES_GENE_NAMES=0,
	INFRN0_AXES_CELL_NAMES=1,
	INFRN0_AXES_GENE_STATE_NAMES=2,
	INFRN0_AXES_CELL_STATE_NAMES=3,
	INFRN0_AXES_GENE_CLUSTER_NAMES=4,
	INFRN0_AXES_CELL_CLUSTER_NAMES=5
};

class Infern0scope{
    public:
    Annotation annot;
    int sillyint;
    SparseMatrix<uint32_t> rawdata; // Raw Data
    SparseMatrix<uint32_t> soupdata;

		// Normalization/reduction scope:

    Tuple<Tuple<double,3u> > gene_scale;
    Tuple<Tuple<double,3u> > cell_scale;

    TMatrix<double> par_R;
    TMatrix<double> par_M;
    TMatrix<double> par_D;
    TMatrix<double> par_S;
    SparseMatrix<double> gene_state;
    SparseMatrix<double> cell_state;

    myHashmap<string, Vector<Anything > > metadata; 



    SparseMatrix<double> pval_data; // Transformed Data
    SparseMatrix<double> data; // Transformed Data
    Trianglix<double> cell_distance;
    Forest<double,2> gene_hierarch;
    Forest<double,2> cell_hierarch;
    Tuple<uint32_t> cell_ordering;
    Tuple<uint32_t> gene_ordering;
    Tuple<uint32_t> cell_color;
    Tuple<uint32_t> gene_color;
    Tuple<uint32_t> cell_clusterID;
    Tuple<uint32_t> gene_clusterID;

    Tuple<double> r_features;


    Tuple<double> P_mean_features;
    Tuple<double> P_var_features;

    DataGrid<double, 2u> scale_class_feature;

    Tuple<double> nbread_features;
    Tuple<double> drop_features;

    

    Tuple<double> drop_features_empirical;


		Infern0scope();
		~Infern0scope();
		//void fitNegativeBinomial();

		uint32_t getNBgenes()const;
		uint32_t getNBcells()const;

		Tuple<double> getBiodropoutProb(const Vector<uint32_t> &rows, uint32_t col);
		Tuple<double> getBiodropoutProb(uint32_t col);
		double getBiodropoutProb(uint32_t row,uint32_t col);
		uint32_t populateProfile(Vector< GaussElem< Tuple<double, 0u > > > &to_clust, Tuple<double, 0u> (&prior)[3], bool is_cell_profile);

		void findDiscriminatoryGenes();
		Tuple<uint32_t, 3u> testfunc(int value);
		void dothattest();
		
		ERRCODE getGeneCellIndices(myHashmap<uint32_t> &rows, myHashmap<uint32_t>&cols, const std::vector<std::string> &genelist, const std::vector<std::string> &celllist, int& missinggene, int& missingcells, bool printwarning) const;

		void clusterHierarchical(bool do_cluster_cells, bool do_overwrite_order);
		void saveHierarchical(const char* path);

		void cmpDeviation();

		void wrSubMatrix(SparseMatrix<double> &target, const std::vector<std::string> genes, const std::vector<std::string> cells, Rcpp::CharacterVector *vecpair, int mincell=3);
//		Tuple<uint32_t> findCleverRandomPartition(uint32_t nb_partition, const SparseMatrix<double> &data, uint32_t seed=0xFFFFFFFF);

		void getHypergeometricDropout(vector<double> &pval, vector<uint32_t> &maxcluster);


		Rcpp::List identifyNetwork(const vector<string> genes, const vector<string> cells, uint32_t nb_thread, uint32_t maxcyclicsize,uint32_t nbcheck,uint32_t mincheck);

    ERRCODE save(FILE*f) const;
    ERRCODE load(FILE*f);
    void chkInit();
	//	void cmpVariance(Trianglix<double> &var, Tuple<double> &mean, const Vector<uint32_t> &gene_list, const Vector<uint32_t> &cell_list, uint32_t nb_threads, unsigned int itemax) const;


	class TaskPartcorrel :public Event<uint32_t>{
        public:
        const Infern0scope& target;
	SparseMatrix<double> trdata;
	Tuple< uint32_t > ranges;
        TMatrix<double,0u,0u> &part;
        myHashmap<uint32_t, uint32_t> mapping;
	const double* minthr;
	bool method;
	uint32_t mincell;
	Tuple<uint32_t> topass;
	HeapTree< KeyElem<double, uint32_t > > heap_correl;
	Tuple< KeyElem<double, uint32_t> > asad;
    	AsyncInserter<HeapTree< KeyElem<double, uint32_t> > > async;

	TaskPartcorrel(const Infern0scope& _target, TMatrix<double,0u,0u> &_cort, bool _met = true ):target(_target),async(heap_correl),part(_cort), method(_met){}
        uint32_t operator()(uint32_t threadID);
	};
	void runTaskPartcorrel(TMatrix<double,0u,0u> &out_correl, Vector<uint32_t> &out_genelist, const Vector<uint32_t> &gene_list, const Vector<uint32_t> &cell_list, uint32_t nbthread=4u, uint32_t indvm=0u, uint32_t mincell =3u,const double* minthr=NULL)const;


	ERRCODE wrDeviationMatrix(Rcpp::NumericMatrix &fout, const myHashmap<uint32_t> &rows, const myHashmap<uint32_t> &cols, int mode) const; // mode: 0 Normalized value, 1 Expected value, 2 Variance, 3 R parameter, 4 M parameter , 5 recomputed Normalized Value

	class DropoutDerivation{
	public:
		const Infern0scope& target;
		DropoutDerivation(const Infern0scope& _target) : target(_target){}
		Tuple<double, 2u> operator()(const Tuple<uint32_t,2u>&) const;
	};
};



static int funtest;

class Parameters{
    public:
    std::vector<double> row_factor;
    std::vector<double> col_factor;
};


