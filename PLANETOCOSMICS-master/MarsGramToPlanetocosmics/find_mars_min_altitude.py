#!/usr/bin/python
import os
import sys
import string
import math
import time


def read_mola_data():
	mola_file=open("/home/laurent/planets/mars/marsgram_model/UnxFiles/txtFiles/molatoph.txt",'r')
	data=mola_file.readlines()
	longitude=[]
	latitude=[]
	altitude=[]
	min_alt=99999999999.
	lat_min=0.
	long_min=0.
	lat_max=0.
	long_max=0.
	max_alt=-99999999999.
	for line in data[1:]:
		words=line.split()
		
		longi=float(words[0])
		lat=float(words[1])
		longitude+=[longi]
		latitude+=[lat]
		#read the altitude in km
		alt=float(words[4])/1000.
		if min_alt > alt:
			min_alt=alt
			lat_min=lat
			long_min=longi
		if max_alt<alt:
			max_alt=alt
			lat_max=lat
			long_max=longi
			
		altitude+=[alt]
	return [latitude,longitude,altitude,min_alt,lat_min,long_min,max_alt,lat_max,long_max]	
	
	

res=read_mola_data()
latitude=res[0]
longitude=res[1]
altitude=res[2]
min_alt=res[3]
lat_min=res[4]
long_min=res[5]
max_alt=res[6]
lat_max=res[7]
long_max=res[8]
print long_min,lat_min,min_alt
print long_max,lat_max,max_alt



