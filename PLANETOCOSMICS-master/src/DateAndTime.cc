#include "DateAndTime.hh"
#include "G4ios.hh"
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime()
{year=2000;
 month=1;
 day=1;
 hour=0;
 min=0;
 sec=0;
 msec=0;
 leap_seconds.clear();
 julday_for_leap_seconds.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d)
{year=y;
 month=mo;
 day=d;
 hour=0;
 min=0;
 sec=0;
 msec=0;
 leap_seconds.clear();
 julday_for_leap_seconds.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d, int h, int mi, int s)
{year=y;
 month=mo;
 day=d;
 hour=h;
 min=mi;
 sec=s;
 msec=0;
 leap_seconds.clear();
 julday_for_leap_seconds.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(int y, int mo, int d, int h, int mi, int s, int ms)
{year=y;
 month=mo;
 day=d;
 hour=h;
 min=mi;
 sec=s;
 msec=ms;
 leap_seconds.clear();
 julday_for_leap_seconds.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::DateAndTime(double julian_date)
{ConvertToJulianDate(julian_date);
 leap_seconds.clear();
 julday_for_leap_seconds.clear();
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime::~DateAndTime()
{;
}	

//Public methods
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInDays(const DateAndTime aDate) const
{return JulianDate()-aDate.JulianDate();
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInHours(const DateAndTime aDate) const
{return DifferenceInDays(aDate) * 24.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInMinutes(const DateAndTime aDate) const
{return DifferenceInDays(aDate) * 24. * 60.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInSeconds(const DateAndTime aDate) const
{return DifferenceInDays(aDate) * 24. * 3600.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::DifferenceInMilliseconds(const DateAndTime aDate)const
{return DifferenceInDays(aDate) * 24. *  3600000.;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::SecondFromDayStart() const
{return DifferenceInSeconds(DateAndTime(year,month,day));
}
////////////////////////////////////////////////////////////////////////////////
//
long DateAndTime::JulianDay() const
{// The julian day 0 is 1 January 4713 BC    
 // Algorithm for computing julian day from Year Month and Day of Gregorian
 // calendar, taken from 
 //"explanatory supplement to the Astronomical Almanach", 1992, p.604 
 
  long ly=long(year);
  long lm=long(month);
  long ld=long(day);
  
  // Adjust BC years
  if ( ly < 0 ) ly++;
  
  long JD= ld - 32075L +
           1461L * (ly + 4800L + (lm - 14L) / 12L) / 4L +
           367L * (lm - 2L - 12L * ((lm - 14L) / 12L)) / 12L -
           3L * ((ly + 4900L + (lm - 14L) / 12L) / 100L) / 4L;

  return JD;
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::JulianDate() const
{// The julian date 0.0 starts at 1 January 4713 BC at Greenwhich mean noon   
 
 
 double dms=double(msec);
 double dsec=double(sec);
 double dmin=double(min);
 double dhour=double(hour);
 double fraction_day = 
          (((dms / 1000. + dsec) / 60. + dmin) / 60. + dhour) / 24. - 0.5;
  	  
 double JulDate = double (JulianDay()) + fraction_day;
 return JulDate;
}
////////////////////////////////////////////////////////////////////////////////
//
int DateAndTime::DayOfYear() const
{// return the day of the year with 1st january =day 0
   int day_of_year = int (DifferenceInDays(DateAndTime(year,1,1))); 
   return day_of_year;  
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::ConvertToJulianDate(double aJulianDate)
{// set the date according to a given julian date

 //Julian day
 // 0.5 is because Julian Date are calciulated form noon on 1st January BC 4713
  
  long JD=long(aJulianDate + 0.5);
  double fraction_of_day =  (aJulianDate + 0.5) - double(JD);
  
  //calculating year month day with algorithm from
  //"explanatory supplement to the Astronomical Almanach", 1992, p.604 
   
   long L = JD + 68569L;
   long N = (4 * L) / 146097L;
   L = L - ((146097 * N + 3) / 4L);
   long I = (4000L * (L + 1))/1461001L;
   L = L - (1461L * I) / 4L + 31L;
   long J = (80L * L) / 2447L;
   day = int ( L - (2447L * J) / 80L );
   L  = J / 11L;
   month = int (J + 2 - (12 * L));
   year = int (100L * (N-49) + I + L);
   
  //calculating hour, min, sec , msec
  
    double dhour  = fraction_of_day * 24.;
    hour = int (dhour);
    double dmin = (dhour - double (hour)) * 60.;
    min = int (dmin);
    double dsec = (dmin - double(min)) * 60.;
    sec = int(dsec);
    msec = int ((dsec - double (sec)) * 1000.); 
    
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::ConvertToSpenvisModifiedJulianDay(double aModJulianDay)
{ ConvertToJulianDate(aModJulianDay + 2433282.5);    
}
////////////////////////////////////////////////////////////////////////////////
//
DateAndTime DateAndTime::TDB_from_UTC()
{double  UTC_juldate=JulianDate();
 double leap_sec =LeapSecond();
 double TAI_juldate = UTC_juldate+ leap_sec/86400.;
 double TDT_juldate = TAI_juldate + 32.184/86400.;
 // first approximation for TDB from Explanatory supplement to the astronomical 
 // almanach (1992) p.42
 double JD = double(int(TDT_juldate*100.))/100.;
 double g_deg= 357.53 + 0.9856003*(JD-2451545.0);
 double g = g_deg*std::acos(0.)/90.;
 double TDB_juldate = TDT_juldate +(0.001658*std::sin(g) + 
                                    0.000014*std::sin(2.*g))/86400.;
 DateAndTime aDate;
 aDate.ConvertToJulianDate(TDB_juldate);				    
 return aDate;				    
}
////////////////////////////////////////////////////////////////////////////////
//
/*DateAndTime DateAndTime::UTC_from_TDB()
{
}*/
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::operator=(const DateAndTime aDate)
{year =aDate.year;
 month=aDate.month;
 day=aDate.day;
 hour=aDate.hour;
 min=aDate.min;
 sec=aDate.sec;
 msec=aDate.msec;
}
////////////////////////////////////////////////////////////////////////////////
//
void DateAndTime::BuildLeapSecondVectors()
{ leap_seconds.clear();
  julday_for_leap_seconds.clear();
  
  
  //leap secons taken form Fraenz, M. and D., Harper, Heliospheric Coordinate
  // Systems, corrected version of Plan. Space Sci., 50, 217,2002
  julday_for_leap_seconds.push_back(DateAndTime(1972,1,1).JulianDate());
  leap_seconds.push_back(10.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1972,7,1).JulianDate());
  leap_seconds.push_back(11.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1973,1,1).JulianDate());
  leap_seconds.push_back(12.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1974,1,1).JulianDate());
  leap_seconds.push_back(13.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1975,1,1).JulianDate());
  leap_seconds.push_back(14.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1976,1,1).JulianDate());
  leap_seconds.push_back(15.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1977,1,1).JulianDate());
  leap_seconds.push_back(16.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1978,1,1).JulianDate());
  leap_seconds.push_back(17.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1979,1,1).JulianDate());
  leap_seconds.push_back(18.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1980,1,1).JulianDate());
  leap_seconds.push_back(19.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1981,7,1).JulianDate());
  leap_seconds.push_back(20.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1982,7,1).JulianDate());
  leap_seconds.push_back(21.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1983,7,1).JulianDate());
  leap_seconds.push_back(22.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1985,7,1).JulianDate());
  leap_seconds.push_back(23.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1988,1,1).JulianDate());
  leap_seconds.push_back(24.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1990,1,1).JulianDate());
  leap_seconds.push_back(25.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1991,1,1).JulianDate());
  leap_seconds.push_back(26.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1992,7,1).JulianDate());
  leap_seconds.push_back(27.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1993,7,1).JulianDate());
  leap_seconds.push_back(28.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1994,7,1).JulianDate());
  leap_seconds.push_back(29.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1996,1,1).JulianDate());
  leap_seconds.push_back(30.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1997,7,1).JulianDate());
  leap_seconds.push_back(31.);
  
  julday_for_leap_seconds.push_back(DateAndTime(1999,1,1).JulianDate());
  leap_seconds.push_back(32.);
}
////////////////////////////////////////////////////////////////////////////////
//
double DateAndTime::LeapSecond()
{ if (leap_seconds.size()<1) BuildLeapSecondVectors();
  double jul_date = JulianDate();
  unsigned int index =0;
  double leap_second=0.;
  while (jul_date >= julday_for_leap_seconds[index] && 
                                      index < leap_seconds.size()){
  	leap_second=leap_seconds[index];
	index++;
  } 
  return leap_second; 
  	 
}

