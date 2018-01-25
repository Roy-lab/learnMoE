#include <fstream>
#include <iostream>
#include <string.h>

#include <unistd.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "MotifManager.H"
#include "BFGSWrapperData.H"
#include "BFGSWrapper.H"
#include "TestData.H"
#include "MotifRegressor.H"
#include "Framework.H"

Framework::Framework()
{
	epsThreshold=-1;
	untransformedData[0]='\0';
}

Framework::~Framework()
{
}

//We will use getopt here
//The options are 
//-m modelname
//-o outputdir
//-e epsilon to control the number of standard deviations above random
//-k maxfactorsize
//-s number of samples for approximate information estimation
//-x k for which we approximate information
//-n cnt of the top candidate MBs to save

Error::ErrorCode
Framework::init(int argc, char** argv)
{
	evManager.setVariableManager(&varManager);
	int optret='-';
	opterr=1;
	int oldoptind=optind;
	cvfolds=1;
	gmmClusterFName[0]='\0';
	while(optret=getopt(argc,argv,"m:o:k:e:x:p:t:l:b:i:f:d:v:g:u:")!=-1)
	{
		if(optret=='?')
		{
			cout <<"Option error " << optopt << endl;
			return Error::UNKNOWN;
		}
		char c;
		char* my_optarg=NULL;
		c=*(argv[oldoptind]+1);
		if(optind-oldoptind ==2)
		{
			my_optarg=argv[oldoptind+1];	
		}
		else
		{
			my_optarg=argv[oldoptind]+2;
		}
		switch(c)
		{
			case 'm':
			{
				char fName[256];
				sprintf(fName,"%s.model",my_optarg);
				Error::ErrorCode eCode=varManager.readVariables(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				sprintf(fName,"%s.data",my_optarg);
				eCode=evManager.loadEvidenceFromFile_Continuous(fName);
				if(eCode!=Error::SUCCESS)
				{
					cout << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				break;
			}
			case 'o':
			{
				strcpy(outputDir,my_optarg);
				break;
			}
			case 'l':
			{
				expertCnt=atoi(my_optarg);
				break;
			}
			case  'i':
			{
				motifManager.readMotifs(my_optarg);
				break;
			}
			case 'p':
			{
				strcpy(projectFName,my_optarg);
				break;
			}
			case 'f':
			{
				motifManager.readMotifTFMap(my_optarg);
				break;
			}
			case 'e':
			{
				if(strcmp(my_optarg,"kmeans")==0)
				{
					initType=MotifRegressor::KMEANS;
				}
				else if(strcmp(my_optarg,"gmm")==0)
				{
					initType=MotifRegressor::GMM;
				}
				else
				{
					initType=MotifRegressor::RAND;
				}
				break;
			}
			case 'd':
			{
				testData.readData(my_optarg);
				break;
			}
			case 'v':
			{
				cvfolds=atoi(my_optarg);
				break;
			}
			case 'g':
			{
				strcpy(gmmClusterFName,my_optarg);
				break;
			}
			case 'u':
			{
				strcpy(untransformedData,my_optarg);
				break;
			}
			default:
			{
				cout <<"Unhandled option " << c  << endl;
				return Error::UNKNOWN;
			}
		}
		oldoptind=optind;
	}
	if(initType==MotifRegressor::GMM)
	{
		if(strlen(gmmClusterFName)==0)
		{
			cout <<"No GMM Cluster file " << endl;
			exit(0);
		}
	}
	return Error::SUCCESS;
}


int 
Framework::start()
{
	MotifRegressor mlearner;
	mlearner.setInitType(initType);
	if(initType==MotifRegressor::GMM)
	{
		mlearner.setInitClusterFile(gmmClusterFName);
	}
	mlearner.setExpertCnt(expertCnt);
	mlearner.setVariableManager(&varManager);
	mlearner.setEvidenceManager(&evManager);
	mlearner.setMotifManager(&motifManager);
	mlearner.setOutputDir(outputDir);
	mlearner.setUntransformedData(untransformedData);
	//mlearner.learnMoE();
	mlearner.learnMoE_CrossValidation(cvfolds);
	//mlearner.showMoE(outputDir);
	//mlearner.showMoEParameters();
	//mlearner.showClusterAssignment(outputDir);
	mlearner.dispTFsPerCluster();
	//mlearner.showGenatomyModule();
	if(testData.getSize()>0)
	{
		mlearner.predictTestData(testData.getDataSet());
	}
        return 0;
}


int
main(int argc, char* argv[])
{
	if(argc<2)
	{
		cout <<"factorGraphInf " <<  endl
			<<"-m modelname " << endl
			<< "-o outputdir " << endl
			 << "-p penalty"<< endl
			<< "-k maxfactorsize " << endl
			 << "-x maxfactorsize_approx" << endl
			 << "-t convergence_threshold" << endl
			 << "-l expertcnt"<< endl
			 << "-b priornet"<< endl
			 << "-i motifinstance" << endl
			 << "-f tfmotifmap"<< endl 
			<< "-e [kmeans|random|gmm]" << endl
			<< "-g gmmclusterfile" << endl
			<< "-d testdata" << endl
			<< "-v crossvalidation_folds" << endl
			<< "-u untransformed_data"<< endl;

		return 0;
	}
	Framework fw;
	if(fw.init(argc,argv)!=Error::SUCCESS)
	{
		return 0;
	}
	fw.start();
	return 0;
	
}

