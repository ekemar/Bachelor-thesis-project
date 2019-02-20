#!/usr/bin/python
import os
import sys
import string
import math
import time

def replace_strings_in_lines(old_strings,new_strings,lines):
	new_lines=lines
	i=0
	for new_line in new_lines :
		j=0
		for str1 in old_strings:
            		str2=new_strings[j]
			new_line=string.replace(new_line,str1,str2) 
			j+=1
		new_lines[i]=new_line
		i+=1
	return new_lines
	
def read_lines_in_file(file_name):
	myfile_in = open(file_name,'r')
       	lines = myfile_in.readlines()
       	myfile_in.close()
	return lines
def write_lines_in_file(lines,file_name):
	myfile_out = open(file_name,'w')
	for line in lines:
		myfile_out.write(line)
	myfile_out.close()		

def replace_strings_in_file(old_strings,new_strings, file_name):
  	lines = read_lines_in_file(file_name)
	return replace_strings_in_lines(old_strings,new_strings,lines)


def write_marsgram_name_list(year,month,day,hour,minute,second,latitude,longitude,template_nml,output_nml):
	new_strings = [str(year),str(month),str(day),str(hour),str(minute),str(second),str(latitude),str(longitude) ]
	old_strings = ["PYEAR","PMONTH", "PDAY", "PHOUR", "PMIN", "PSEC", "PLAT", "PLONG" ]
	lines = replace_strings_in_file(old_strings,new_strings, template_nml)
	write_lines_in_file(lines,output_nml)
	
	

def get_parameter_value_to_name_list(name_list,tag_name):
	file_nml = open (name_list,'r')
	lines =file_nml.readlines()
	res='unfound'
	end_tag_reached = 0
	for line in lines:
		index = string.find(line,"$END")
		if not (index  == -1):
			end_tag_reached = 1
		if (end_tag_reached == 0):	 
			index = string.find(line,tag_name)
			if not (index  == -1): 
				words = string.split(line)
				res = words[2]
	
	file_nml.close()
	return res
def compute_table_for_specific_conditions(year,month,day,hour,minute,sec,lat,longi,output_file_name):
	write_marsgram_name_list(str(int(year)),str(int(month)),str(int(day)),str(int(hour)),str(int(minute)),str(int(sec)),str(lat),str(longi),"marsgram_template.nml","run.nml")	
	os.system("echo run.nml | ./marsgram.out")
	output_file=open (output_file_name,'w')
	name_list="run.nml"
	mgram_output_file=get_parameter_value_to_name_list(name_list,"OUTFL")
	year = get_parameter_value_to_name_list(name_list,"MYEAR")
	month = get_parameter_value_to_name_list(name_list,"MONTH")
	day = get_parameter_value_to_name_list(name_list,"MDAY")
	hour = get_parameter_value_to_name_list(name_list,"IHR")
	minute = get_parameter_value_to_name_list(name_list,"IMIN")
	sec = get_parameter_value_to_name_list(name_list,"SEC")
	latitude= get_parameter_value_to_name_list(name_list,"FLAT")
	longitude = get_parameter_value_to_name_list(name_list,"FLON")



	file_day_data = open ("DayData.txt",'r')
	lines = file_day_data.readlines()
	words=string.split(lines[1])
	altitude = str(float(words[0]) - float(words[1]))

	
	output_file.write("\\comments\n")
	output_file.write("\t Dayly averaged Mars atmospheric profiles computed with MarsGram2001 \n")  
	output_file.write("\tdate :\t"+year+"/"+month+"/"+day+" "+hour+":"+minute+":"+sec+"\n")
	output_file.write("\tground altitude :\t"+altitude+" km\n")
	output_file.write("\tlatitude :\t"+latitude+"\n")
	output_file.write("\tlongitude :\t"+longitude+"\n")

	output_file.write("\\definition\n")
	output_file.write("\t \\type_of_composition{mass_composition}\n")
	output_file.write("\t \\mass_density_unit{kg/m3}\n")
	output_file.write("\t \\pressure_unit{Pa}\n")
	output_file.write("\t \\altitude_unit{km}\n")
	output_file.write("\t \\ground_altitude{"+altitude+"}\n")
	output_file.write("\t \\top_altitude{90}\n")


	output_file.write("\\data\n")
	output_file.write("\taltitude\ttemperature\tpressure\tdensity\tCO2\tN2\tAr\n")

	lines_newfile=["\n"]

	for line in lines[1:1000]:
		word=string.split(line)
		height= float(word[1])
		altitude = float(word[0])
		if (height >= 0 and altitude<=100) :
			new_line="\t"+word[0]+"\t"+word[2]+"\t"+word[3]+"\t"
			new_line+=word[4]+"\t"+str(0.957)+"\t"+str(0.027)+"\t"+str(0.016)+"\n"
			lines_newfile=[new_line]+lines_newfile

	for line in  lines_newfile:  
		output_file.write(line)

def compute_table_for_a_grid(year,month,day,hour,minute,sec,lat0,dlat,nlat,long0,dlong,nlong,pre_file_name):
	for i in range(nlat):
		for j in range(nlong):	
			latitude=lat0+i*(dlat)
			longitude=long0+j*(dlong)
			post_str='_%.2f_%.2f.ascii' %(latitude,longitude)
			file_name=pre_file_name+post_str
			print file_name
			compute_table_for_specific_conditions(year,month,day,hour,minute,sec,latitude,longitude,file_name)
	
	
	
PYEAR=float(sys.argv[1])
PMONTH=	float(sys.argv[2])
PDAY=float(sys.argv[3])			 
PHOUR=float(sys.argv[4])
PMIN=float(sys.argv[5])
PSEC=float(sys.argv[6])
output_file_name=sys.argv[7]
compute_table_for_a_grid(PYEAR,PMONTH,PDAY,PHOUR,PMIN,PSEC,-85.,10.,18,0.,90.,4, output_file_name)
