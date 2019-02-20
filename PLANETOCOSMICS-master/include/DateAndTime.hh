#ifndef DATEANDTIME_HH
#define DATEANDTIME_HH
#include <vector>
class DateAndTime 
{
public:
	//constructor destructor
	DateAndTime();
	DateAndTime(int y, int mo, int day);
	DateAndTime(int y, int mo, int day, int h, int mi, int s);
	DateAndTime(int y, int mo, int day, int h, int mi, int s, int msec);
	DateAndTime(double julian_date);
	~DateAndTime();
	
	//public methods
	
	double DifferenceInDays(const DateAndTime) const;
	double DifferenceInHours(const DateAndTime) const;
	double DifferenceInMinutes(const DateAndTime) const;
	double DifferenceInSeconds(const DateAndTime) const;
	double DifferenceInMilliseconds(const DateAndTime)const;
	long   JulianDay() const;
	inline double   SpenvisModifiedJulianDay()const{return (JulianDate() - 2433282.5);}
        double JulianDate() const;
	int   DayOfYear() const;
	void ConvertToJulianDate(double aJulianDate); 
	void ConvertToSpenvisModifiedJulianDay(double aModJulianDay);
	//conversion from UTC time to barycentric dynamical time TDB
	double SecondFromDayStart() const;
	DateAndTime TDB_from_UTC();
	//conversion from barycentric dynamical time TDB to UTC time
	//DateAndTime UTC_from_TDB();
	
	//operators
	void operator=(const DateAndTime);
	
	
private:
       //private methods
       void BuildLeapSecondVectors();
       double LeapSecond();

//attributes
public: //public attributes
        int year,month,day,hour,min,sec,msec;
	
private: 
	//leap second vectors
	std::vector<double> leap_seconds;
	std::vector<double> julday_for_leap_seconds;
	
		
	
	
};
#endif
