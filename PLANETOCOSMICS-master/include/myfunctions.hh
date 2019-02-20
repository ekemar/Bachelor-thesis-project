#ifndef MYFUNC_h
#define MYFUNC_h 1
#include <fstream>
#include <complex>
#include <vector>

#ifndef USE_ANALYSIS_ROOT
#define MYFUNCH1D IHistogram1D
namespace AIDA
{class IHistogram1D;
}
using namespace AIDA;
#else
#define MYFUNCH1D TH1F
class TH1F;
#endif 
namespace myfunc{
//integration of y(x) from x1 to x2 by considering linear variation of y with x
double IntegrationOfY(const std::vector< double > x, 
                      const std::vector< double > y,
		      double x1, double x2); //integration of y vs x from x1 to x2
//integration of y(x) from x1 to x2 by considering exponential variation of y with x
double IntegrationOfY_exp(const std::vector< double > x, 
                      const std::vector< double > y,
		      double x1, double x2);
//integration of y(x) from x1 to x2 by considering power law variation of y with x
double IntegrationOfY_pow(const std::vector< double > x, 
                      const std::vector< double > y,
		      double x1, double x2);		      
// mean value of y(x) from x1 to x2		      		     
double MeanValueOfY(const  std::vector< double> x,
                    const  std::vector< double> y,
		    double x1,double x2);

//linear interpolation of y(x) at x1 
double LinearInterpolation(const  std::vector< double> x,
                           const  std::vector< double> y,
			   const double x1);
double LogInterpolation(const  std::vector< double> x,
                           const  std::vector< double> y,
			   const double x1);			   
//Bilinear extra/interpolation of z(x,y) at x1,y1  
double BilinearInterpolation(const  std::vector< double> x,
                            const  std::vector< double> y,
                            const std::vector< std::vector< double> > z,
			    const double x1, const double y1);
//locate x1 in a vector x
//return i such that x is in the segment [x[i],x[i+1]] or [x[i+1],x[i]]
// x is monotically increasing or decreasing vector
// if x is out of range it returns 0 or x.size()-2 and out_of_range is true



unsigned int locate(const  std::vector< double> x, double x1, bool& out_of_range);
 			    

// dsitribute linearly a weight over  an [x1,x2] intio the histogram aHistogram
void LinearDistributionInHistogram(MYFUNCH1D* anHisto, 
                                   double x1,
				   double x2,
				   double weight);

}
#endif
