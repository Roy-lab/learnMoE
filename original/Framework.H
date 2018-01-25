/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

//All the managers need forward declarations
class Error;
class VariableManager;
class EvidenceManager;
class MotifManager;
class MotifRegressor;
class TestData;
class Framework
{
	public:
		Framework();
		~Framework();
		Error::ErrorCode init(int, char**);
		int start();
		int predictTestData(const char*);
	private:
		VariableManager varManager;
		EvidenceManager evManager;
		MotifManager motifManager;
		TestData testData;
		double epsThreshold;
		char actionStr[256];
		char outputDir[1024];
		double penalty;
		int expertCnt;
		char projectFName[1024];
		MotifRegressor::InitType initType;
		int cvfolds;
		char gmmClusterFName[1024];
		char untransformedData[1024];
};
#endif
