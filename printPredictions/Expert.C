#include <math.h>
#include <gsl/gsl_randist.h>
#include "Matrix.H"
#include "Expert.H"


Expert::Expert()
{
	covariance=NULL;
	mean=NULL;
	invCovariance=NULL;
}

Expert::~Expert()
{
	if(covariance!=NULL)
	{
		delete covariance;
	}
	if(mean!=NULL)
	{
		delete mean;
	}
	if(invCovariance!=NULL)
	{
		delete invCovariance;
	}
}


int 
Expert::setMean(Matrix* m)
{
	if(mean!=NULL)
	{
		delete mean;
	}
	mean=m;
	return 0;
}

int 
Expert::setCovariance(Matrix* c)
{
	if(covariance!=NULL)
	{
		delete covariance;
	}
	covariance=c;
	//c->showMatrix();
	if(invCovariance!=NULL)
	{
		delete invCovariance;
	}
	//invCovariance=covariance->invMatrix();
	invCovariance=new Matrix(c->getRowCnt(),c->getColCnt());
	invCovariance->setAllValues(0);
	double det=1;
	double n=((double)c->getRowCnt());
	for(int i=0;i<n;i++)
	{
		double invval=c->getValue(i,i);
		invCovariance->setValue(1.0/invval,i,i);
		det=det+log(invval);
		//det=det*(c->getValue(i,i));
	}
	double testdet=c->detMatrix();
	//normFactor=pow(2*PI,n)*det;
	normFactor=(n*log(2*PI)) + det;
	//normFactor=sqrt(normFactor);
	normFactor=normFactor/2;
	return 0;
}


int
Expert::updateCovariance()
{
	if(invCovariance!=NULL)
	{
		delete invCovariance;
	}
	invCovariance=covariance->invMatrix();
	//double det=covariance->detMatrix();
	double n=((double)covariance->getRowCnt());
	//normFactor=pow(2*PI,n)*det;
	//normFactor=sqrt(normFactor);
	/*invCovariance=new Matrix(covariance->getRowCnt(),covariance->getColCnt());
	invCovariance->setAllValues(0);*/
	double det=0;
	for(int i=0;i<n;i++)
	{
		double invval=covariance->getValue(i,i);
		//invCovariance->setValue(1.0/invval,i,i);
		det=det+log(invval);
		//det=det*covariance->getValue(i,i);
	}
	
	//normFactor=pow(2*PI,n)*det;
	normFactor=(n*log(2*PI)) + det;
	//normFactor=sqrt(normFactor);
	normFactor=normFactor/2;
	return 0;
}


int
Expert::setFeatureWeight(int fId, double weight)
{
	if(featureWeights.find(fId)!=featureWeights.end())
	{
	//	cout <<"Old feature wt " << featureWeights[fId] << " " << weight << endl;
	}
	featureWeights[fId]=weight;
	return 0;	
}

double 
Expert::getOutputPDF(Matrix* y)
{
	double pdf=0;
	Matrix diffMat(y->getRowCnt(),y->getColCnt());
	int colCnt=y->getColCnt();
	for(int i=0;i<colCnt;i++)
	{
		double diff=y->getValue(0,i)-mean->getValue(0,i);
		diffMat.setValue(diff,0,i);
	}
	Matrix* p1=diffMat.multiplyMatrix(invCovariance);
	double sum=0;
	for(int i=0;i<colCnt;i++)
	{
		double s=p1->getValue(0,i)*diffMat.getValue(0,i);
		sum=sum+s;
	}
	//pdf=exp(-0.5*sum)/normFactor;
	double lpdf=(-0.5*sum)-normFactor;
	pdf=exp(pdf);
	if(pdf<1e-300)
	{
		cout << "Correcting " << pdf << " to " << 1e-80 << endl;
		pdf=1e-300;
	}
	delete p1;
	return lpdf;
//	return pdf;
}

double 
Expert::getOutputPDF_Fast(Matrix* y)
{
	double pdf=0;
	Matrix diffMat(y->getRowCnt(),y->getColCnt());
	int colCnt=y->getColCnt();
	double sum=0;
	for(int i=0;i<colCnt;i++)
	{
		double diff=y->getValue(0,i)-mean->getValue(0,i);
		double corre=invCovariance->getValue(i,i);
		double s=diff*corre*diff;
		sum=sum+s;
	}
	//pdf=exp(-0.5*sum)/normFactor;
	double lpdf=(-0.5*sum)-normFactor;
	pdf=exp(pdf);
	if(pdf<1e-300)
	{
		cout << "Correcting " << pdf << " to " << 1e-80 << endl;
		pdf=1e-300;
	}
	return lpdf;
//	return pdf;
}


double 
Expert::getOutputPDF_Nocov(Matrix* y)
{
	double pdf=0;
	Matrix diffMat(y->getRowCnt(),y->getColCnt());
	int colCnt=y->getColCnt();
	double lpdf=0;
	for(int i=0;i<colCnt;i++)
	{
		double mean_i=mean->getValue(0,i);
		double var_i=covariance->getValue(i,i);
		double val_i=y->getValue(0,i);
		double gpdf=gsl_ran_gaussian_pdf(val_i-mean_i,sqrt(var_i));
		if(gpdf<1e-30)
		{
			gpdf=1e-30;
		}
		lpdf=lpdf+log(gpdf);
	}
	pdf=exp(lpdf);
	if(pdf<1e-300)
	{
		cout << "Correcting " << pdf  << " logval="  <<lpdf << " to " << 1e-80 << endl;
		pdf=1e-300;
	}
	return pdf;
}




double 
Expert::getMixtureWeight(INTDBLMAP* x)
{
	double unnormWeight=0;
	for(INTDBLMAP_ITER wIter=featureWeights.begin();wIter!=featureWeights.end();wIter++)
	{
		if(x->find(wIter->first)!=x->end())
		{
			double val=wIter->second*(*x)[wIter->first];
			unnormWeight=unnormWeight+val;
		}
	}
	return unnormWeight;
}

Matrix*
Expert::getMean()
{
	return mean;
}

Matrix*
Expert::getCovariance()
{
	return covariance;
}


INTDBLMAP& 
Expert::getFeatureWeights()
{
	return featureWeights;
}


int 
Expert::sortFeatures()
{
	for(INTDBLMAP_ITER fIter=featureWeights.begin();fIter!=featureWeights.end();fIter++)
	{
		sortedFeatureIDs.push_back(fIter->first);
	}
	for(int i=0;i<sortedFeatureIDs.size();i++)
	{
		for(int j=i+1;j<sortedFeatureIDs.size();j++)
		{
				double aVal=fabs(featureWeights[sortedFeatureIDs[i]]);
				double bVal=fabs(featureWeights[sortedFeatureIDs[j]]);
				if(aVal<bVal)
				{
					int aId=sortedFeatureIDs[i];
					sortedFeatureIDs[i]=sortedFeatureIDs[j];
					sortedFeatureIDs[j]=aId;
				}
		}
	}
	return 0;
}


int 
Expert::sortFeatures_Enrichment(INTDBLMAP& fEnr)
{
	sortedFeatureIDs.clear();
	for(INTDBLMAP_ITER fIter=fEnr.begin();fIter!=fEnr.end();fIter++)
	{
		sortedFeatureIDs.push_back(fIter->first);
	}
	for(int i=0;i<sortedFeatureIDs.size();i++)
	{
		for(int j=i+1;j<sortedFeatureIDs.size();j++)
		{
			double aVal=fEnr[sortedFeatureIDs[i]];
			double bVal=fEnr[sortedFeatureIDs[j]];
			if(aVal>bVal)
			{
				int aId=sortedFeatureIDs[i];
				sortedFeatureIDs[i]=sortedFeatureIDs[j];
				sortedFeatureIDs[j]=aId;
			}
		}
	}
	return 0;
}


vector<int>& 
Expert::getSortedFeatures()
{
	return sortedFeatureIDs;
}

int
Expert::assignGeneToExpert(int g)
{
	geneSet[g]=0;
	return 0;
}

map<int,int>&
Expert::getGeneSet()
{
	return geneSet;
}

int
Expert::resetAssignedGenes()
{
	geneSet.clear();
}

double
Expert::getEntropy()
{
	double determinant=covariance->detMatrix();
	//covariance->showMatrix();
	//cout <<"Determinant: " << determinant << endl;
	double n=((double)covariance->getRowCnt())/2.0;
	double commFact=1+log(2*3.1472);
	double jointEntropy=0.5*((n*commFact) + log(determinant));
	//cout <<"Joint entropy: " << jointEntropy << endl;
	return jointEntropy;
}


int 
Expert::setPrior(double p)
{
	prior=p;
	return 0;
}

double 
Expert::getPrior()
{
	return prior;
}
