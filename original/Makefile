LFLAG = -lgsl -lgslcblas 
SRC = BFGSWrapper.C ClusterManager.C Expert.C Framework.C Kmeans.C  MotifManager.C MotifRegressor.C Randomizer.C TestData.C  Distance.C HyperGeomPval.C GeneExpManager.C common/Variable.C common/VariableManager.C common/Evidence.C common/EvidenceManager.C common/Error.C  common/Matrix.C
LIBPATH = ../../thirdparty/gsl_install/lib
INCLPATH1 = ../../thirdparty/gsl_install/include
INCLPATH2 = common

CC=g++
CFLAGS = -g

learnMoE: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o learnMoE 

