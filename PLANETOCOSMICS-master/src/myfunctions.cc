#include"myfunctions.hh"
#include "globals.hh"
#ifdef USE_ANALYSIS_ROOT
#define HISTO1D TH1F
#include "TH1.h"
#include "TAxis.h"
#else
#define HISTO1D IHistogram1D
#include "AIDA/IHistogram1D.h"
#include "AIDA/IAxis.h"
#endif

double myfunc::IntegrationOfY(const std::vector< double > x, 
                              const std::vector< double > y,
		              double x1, double x2)
{//y is always positive and x is monotically decreasing or increasing
 // we consider the absolute value of the integration
   
   unsigned int i1,i2,nx;
   bool out_of_range1,out_of_range2 ;
   i1 = myfunc::locate(x,x1,out_of_range1);
   i2 = myfunc::locate(x,x2,out_of_range2);
   nx=x.size();
  
   double y1,y2;
   y1=y[i1] + (x1-x[i1])*(y[i1+1]-y[i1])/(x[i1+1]-x[i1]);
   y2=y[i2] + (x2-x[i2])*(y[i2+1]-y[i2])/(x[i2+1]-x[i2]);
   
   int ii1,ii2;
   ii1=i1;
   ii2=i2;
   std::vector <double> x_vec,y_vec;
   
   //direction of vectors from x1 or x2 or from x2 to x1  
   
   if  ((x[0]-x1)/(x[0]-x[1]) < (x[0]-x2)/(x[0]-x[1])){
     	x_vec.push_back(x1);
      	y_vec.push_back(y1);
      	if ( !out_of_range1 || i1 ==0) {
        	ii1=i1+1;
	 	if  (out_of_range1) ii1=i1; // below 0
	 	ii2=i2;
	 	if  (out_of_range2 && i2 > 0) ii2=i2+1;
	 	else if  ((x[i2]-x2)/(x[i2]-x[i2+1]) <= 0) ii2=i2-1;
	 	for (int j=ii1;j<ii2+1;j++){
	     		x_vec.push_back(x[j]);
	      		y_vec.push_back(y[j]);
	     	}
	}
       	x_vec.push_back(x2);
       	y_vec.push_back(y2);
   }		       
   else {
  	x_vec.push_back(x2);
       	y_vec.push_back(y2);
       	if (!out_of_range2 || i2 ==0){ 
        	ii2=i2+1;
	 	if  (out_of_range2) ii2=i2;
	 	ii1=i1;
	 	if  (out_of_range1 && i1 > 0) ii1=i1+1;
	 	if  ((x[i1]-x1)/(x[i1]-x[i1+1]) <= 0) ii1=i1-1;
	 	for (int j=ii2;j<ii1+1;j++){
	     		x_vec.push_back(x[j]);
	      		y_vec.push_back(y[j]);
	     	}
	}
	x_vec.push_back(x1);
	y_vec.push_back(y1);
   }
    
    //now comes the integration 
       
   double integration_y=0.;
   for (unsigned int j=0;j<x_vec.size()-1;j++)
         	integration_y+=(x_vec[j+1]-x_vec[j])
	                      		* (y_vec[j+1] + y_vec[j])/2.;
       
   return std::abs(integration_y);	 
}
////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////		      
double myfunc::IntegrationOfY_exp(const std::vector< double > x, 
                                  const std::vector< double > y,
		                  double x1, double x2)
{  
   //return myfunc::IntegrationOfY(x,y,x1,x2);

 //y is always positive and x is monotically decreasing or increasing
 // we consider the absolute value of the integration
   
   unsigned int i1,i2,nx;
   bool out_of_range1,out_of_range2 ;
   i1 = myfunc::locate(x,x1,out_of_range1);
   i2 = myfunc::locate(x,x2,out_of_range2);
   
   
   nx=x.size();
  
   double y1,y2;
   y1=y[i1] *std::pow(y[i1+1]/y[i1] ,(x1-x[i1])/(x[i1+1]-x[i1]));
   y2=y[i2] *std::pow(y[i2+1]/y[i2] ,(x2-x[i2])/(x[i2+1]-x[i2]));
   
   int ii1,ii2;
   ii1=i1;
   ii2=i2;
   std::vector <double> x_vec,y_vec;
   
   //direction of vectors from x1 or x2 or from x2 to x1  
   
   if  ((x[0]-x1)/(x[0]-x[1]) < (x[0]-x2)/(x[0]-x[1]))
     {x_vec.push_back(x1);
      y_vec.push_back(y1);
      if ( !out_of_range1 || i1 ==0) 
        {ii1=i1+1;
	 if  (out_of_range1) ii1=i1; // below 0
	 ii2=i2;
	 if  (out_of_range2 && i2 > 0) ii2=i2+1;
	 else if  ((x[i2]-x2)/(x[i2]-x[i2+1]) <= 0) ii2=i2-1;
	 for (int j=ii1;j<ii2+1;j++)
	     {x_vec.push_back(x[j]);
	      y_vec.push_back(y[j]);
	     }
	}
       x_vec.push_back(x2);
       y_vec.push_back(y2);
      }		       
   else{x_vec.push_back(x2);
       y_vec.push_back(y2);
       if (!out_of_range2 || i2 ==0) 
        {ii2=i2+1;
	 if  (out_of_range2) ii2=i2;
	 ii1=i1;
	 if  (out_of_range1 && i1 > 0) ii1=i1+1;
	 if  ((x[i1]-x1)/(x[i1]-x[i1+1]) <= 0) ii1=i1-1;
	 for (int j=ii2;j<ii1+1;j++)
	     {x_vec.push_back(x[j]);
	      y_vec.push_back(y[j]);
	     }
	 }
	 x_vec.push_back(x1);
	 y_vec.push_back(y1);
       }
   
       
   //now comes the integration 
       double integration_y=0.;
       for (unsigned int j=0;j<x_vec.size()-1;j++)
        {if (std::abs (y_vec[j+1] -y_vec[j] ) >0.0000000001*y_vec[j] &&
	                                  y_vec[j] >= 0.00000000000001)
          {integration_y+=
	        std::abs(x_vec[j+1]-x_vec[j])* (y_vec[j+1] - y_vec[j])
	                                  /   std::log(y_vec[j+1]/y_vec[j]);
	  }
	else //we should use the lineear expression
	  {integration_y+=(x_vec[j+1]-x_vec[j])
	                      * (y_vec[j+1] + y_vec[j])/2.;
	  }
	}    				  
       
       return std::abs(integration_y);	
}
////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////		      
double myfunc::IntegrationOfY_pow(const std::vector< double > x, 
                                  const std::vector< double > y,
		                  double x1, double x2)
{//y is always positive and x is monotically decreasing or increasing
 // we consider the absolute value of the integration
   
   unsigned int i1,i2,nx;
   bool out_of_range1,out_of_range2 ;
   i1 = myfunc::locate(x,x1,out_of_range1);
   i2 = myfunc::locate(x,x2,out_of_range2);
   nx=x.size();
  
   double y1,y2,log1,log2;
   double dlogy_dlogx= (std::log(y[i1+1])-std::log(y[i1])) /
                                (std::log(x[i1+1])-std::log(x[i1]));
   log1 =std::log(y[i1]) + (std::log(x1)-std::log(x[i1]))*dlogy_dlogx;				       
   y1= std::exp(log1);
   
   dlogy_dlogx= (std::log(y[i2+1])-std::log(y[i2])) /
                                (std::log(x[i2+1])-std::log(x[i2]));
   log2 =std::log(y[i2]) + (std::log(x2)-std::log(x[i2]))*dlogy_dlogx;				       
   y2= std::exp(log2);
  
   
   
  
   
   int ii1,ii2;
   ii1=i1;
   ii2=i2;
   std::vector <double> x_vec,y_vec;
   
   //direction of vectors from x1 or x2 or from x2 to x1  
   
   if  ((x[0]-x1)/(x[0]-x[1]) < (x[0]-x2)/(x[0]-x[1]))
     {x_vec.push_back(x1);
      y_vec.push_back(y1);
      if ( !out_of_range1 || i1 ==0) 
        {ii1=i1+1;
	 if  (out_of_range1) ii1=i1; // below 0
	 ii2=i2;
	 if  (out_of_range2 && i2 > 0) ii2=i2+1;
	 else if  ((x[i2]-x2)/(x[i2]-x[i2+1]) <= 0) ii2=i2-1;
	 for (int j=ii1;j<ii2+1;j++)
	     {x_vec.push_back(x[j]);
	      y_vec.push_back(y[j]);
	     }
	}
       x_vec.push_back(x2);
       y_vec.push_back(y2);
      }		       
   else
      {x_vec.push_back(x2);
       y_vec.push_back(y2);
       if (!out_of_range2 || i2 ==0) 
        {ii2=i2+1;
	 if  (out_of_range2) ii2=i2;
	 ii1=i1;
	 if  (out_of_range1 && i1 > 0) ii1=i1+1;
	 if  ((x[i1]-x1)/(x[i1]-x[i1+1]) <= 0) ii1=i1-1;
	 for (int j=ii2;j<ii1+1;j++)
	     {x_vec.push_back(x[j]);
	      y_vec.push_back(y[j]);
	     }
	 }
	 x_vec.push_back(x1);
	 y_vec.push_back(y1);
       }
   
       
   //now comes the integration 
       double integration_y=0.;
       for (unsigned int j=0;j<x_vec.size()-1;j++)
        {double b= (std::log(y_vec[j+1])-std::log(y_vec[j]))
	                / (std::log(x_vec[j+1])-std::log(x_vec[j]));
	
	//std::cout<<b<<std::endl;
	//std::cout<<x_vec[j+1]<<std::endl;
	//std::cout<<y_vec[j+1]<<std::endl;
	 if ( std::abs( b + 1) >0.0000001)
	    integration_y += (x_vec[j+1]*y_vec[j+1] - x_vec[j]*y_vec[j])/(b+1);
	 else    
	    integration_y += x_vec[j+1]*y_vec[j+1] * 
	                         (std::log(x_vec[j+1]) -std::log(x_vec[j]));
	 
	}    				  
       
       return std::abs(integration_y);	
}



////////////////////////////////////////////////////////////////////////
double myfunc::MeanValueOfY(const  std::vector< double> x,
                            const  std::vector< double> y,
		            double x1,double x2)
{return IntegrationOfY(x,y,x1,x2)/std::abs(x1-x2);
}
//////////////////////////////////////////////////////////////////////////  
double myfunc::LinearInterpolation(const  std::vector< double> x,
                                     const  std::vector< double> y,
                                     const double x1)
{
 bool out_of_range;
 unsigned int ix = myfunc::locate(x,x1,out_of_range);

 
 double y0 = y[ix];
 double y1 = y[ix+1];
 double t = (x1 - x[ix])/(x[ix+1]-x[ix]);
 
 double res = y0 + (y1-y0) *t;    
 
 return res;
}
//////////////////////////////////////////////////////////////////////////  
double myfunc::LogInterpolation(const  std::vector< double> x,
                                     const  std::vector< double> y,
                                     const double x1)
{
 bool out_of_range;
 unsigned int ix = myfunc::locate(x,x1,out_of_range);

 
 double logy0 = std::log(y[ix]);
 double logy1 = std::log(y[ix+1]);
 double logx0 = std::log(x[ix]);
 double logx1 = std::log(x[ix+1]);
 double logy = logy0 +(std::log(x1)-logx0)*(logy1 -logy0)/(logx1-logx0);
 
 double res = std::exp(logy);    
 
 return res;
}


			  
//////////////////////////////////////////////////////////////////////////  
double myfunc::BilinearInterpolation(const  std::vector< double> x,
                                     const  std::vector< double> y,
                                     const std::vector< std::vector< double> > z,
			             const double x1, const double y1)
{bool out_of_range;
 unsigned int ix = myfunc::locate(x,x1,out_of_range);
 unsigned int iy = myfunc::locate(y,y1,out_of_range);
 
 double z0,z1,z2,z3;
 z0 = z[ix][iy];
 z1 = z[ix+1][iy];
 z2 = z[ix][iy+1];
 z3 = z[ix+1][iy+1];
 
 double t = (x1 - x[ix])/(x[ix+1]-x[ix]);
 double u = (y1 - y[iy])/(y[iy+1]-y[iy]);
 
 double res = (1. - t) * (1. -u) * z0 + 
                  t    * (1. -u) * z1 +
	      (1. - t) *    u    * z2 +
	          t    *    u    * z3;	   
 
 return res;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


unsigned int myfunc::locate( const  std::vector< double> x, double x1, bool& out_of_range)
{// copy and arranged from Numerical recipies in C++ pp 120.
 int ju, jm, jl;
 bool increase;
 int n= int (x.size());
 
 if (n<2)
   {//std::cout<<"the size of x is lower than 2"<<std::endl;
    return 0;   
   }
 jl=-1;
 jm=0;
 ju= n;
 increase = (x[n-1]>=x[0]);
 out_of_range = false;
 
 if (x1 == x[0]) jl=0;
 else if (x1 == x[n-1]) jl = n-2;
 else
   {while (ju-jl >1)
     {jm =(ju+jl) >>1;//shifting of one bit to the right is like devided by 2
      if ( (x1 >= x[jm]) == increase)  jl =jm;
      else                             ju =jm;
     }
   } 
 
 // check if x1 is out of range
 if (jl == -1 || jl == n-1  ) 
    {jl=std::abs(jl) -1;
     out_of_range = true;
    }  

      

 return (unsigned int) jl; 
}
///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//works for an histogram with increasing linear binning scheme 
#ifdef USE_ANALYSIS_ROOT
#define fill Fill
#endif
void myfunc::LinearDistributionInHistogram(HISTO1D* anHisto, 
                                   double x1,
				   double x2,
				   double weight)
{ if (x1==x2){
  	anHisto->fill(x1,weight);
   	return;
  }

  double xmin,xmax,low,up;
  int nbins;
  int n1=-1;
  int n2=-1;
  double weight_per_bin=0.;
  double weight_n1 =0.;
  double weight_n2 =0.;
#ifndef USE_ANALYSIS_ROOT 
  low=anHisto->axis().lowerEdge();
  up=anHisto->axis().upperEdge();
  nbins=anHisto->axis().bins();
#else
  low=anHisto->GetXaxis()->GetXmin();
  up=anHisto->GetXaxis()->GetXmax();
  nbins=anHisto->GetXaxis()->GetNbins();
 
#endif 
  double dx_bin= (up-low)/nbins;
 
  xmin=std::max(low,std::min(up,std::min(x1,x2)))+0.0000000000000001;
  xmax=std::min(up,std::max(low,std::max(x1,x2)))-0.0000000000000001;
 
  if (xmax>xmin){   
#ifndef USE_ANALYSIS_ROOT    
    	n1=anHisto->coordToIndex(xmin);
    	//std::cout<<"n1"<<n1<<std::endl;
    	n2=anHisto->coordToIndex(xmax);
#else
    	n1=anHisto->FindBin(xmin)-1;
    	n2=anHisto->FindBin(xmax)-1;	
#endif    
    	//std::cout<<"n2"<<n2<<std::endl;
    	weight_per_bin=weight*dx_bin/(x1-x2);
    	if (weight_per_bin <0) weight_per_bin = -weight_per_bin;
    	//std::cout<<"weight_per_bin "<<weight_per_bin<<std::endl;
    	double xx1=low+(n1+1)*dx_bin;
    	double xx2=low+(n2)*dx_bin;
    	//std::cout<<"xx1"<<xx1<<std::endl;
    	//std::cout<<"xx2"<<xx2<<std::endl;
    
    	weight_n1=weight*(xx1-xmin)/(x1-x2);
    	if (weight_n1 <0) weight_n1=-weight_n1; 
    	weight_n2=weight*(xmax-xx2)/(x1-x2);
    	if (weight_n2 <0) weight_n2=-weight_n2;
  }
  if (n1>-1){
     	if (n1==n2) {
		anHisto->fill(xmin,weight);
		if (weight <0) G4cout<<weight<<std::endl; 
	}	
      	else {
		anHisto->fill(xmin,weight_n1);
	 	//std::cout<<"xmin "<<xmin<<"weight n1 "<<weight_n1<<std::endl;
	 	//std::cout<<"xmax "<<xmax<<"weight n2 "<<weight_n2<<std::endl;
	 	anHisto->fill(xmax,weight_n2);
	 	for (int j=n1+1;j<n2;j++){
	     		double x=low+ (j+0.5)*dx_bin;
	      		//std::cout<<"x"<<x<<"weight "<<weight_per_bin<<std::endl;
	       		anHisto->fill(x,weight_per_bin);
	     	}				     				     
	}
  }
  else return;
}   
