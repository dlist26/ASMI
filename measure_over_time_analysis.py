import csv
import sys
import time
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.optimize import curve_fit
import os

print("This is the analysis program for the measure_over_time script. This will not work for data collected with the "
      "measure and custom_measure programs")
bad_name = True
files = os.listdir("/home/robot/ASMI_KABlab") #change to correct directory for your device!!!
while bad_name:
    filename = input(
        "Please enter the name of the the file you would like to analyze the data from. The name is case"
        " sensitive, so please enter the name in exactly. To see a list of all files, press Ctrl + C and"
        " type ls in the command window, then rerun this program. Note, please do not type in the .csv, just enter "
        "the name of the file. ")
    filename = filename + "_results_timed.csv"
    if filename in files:
        print("File was found")
        bad_name = False
    else:
        print("File was not found, please try again")
print(filename)
with open(filename, 'r') as file:
    reader = csv.reader(file)
    data = list(reader)
cleaned_data = []
for i in range(0, len(data)):
    if data[i] != []:
        cleaned_data.append(data[i])

cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
invalid = True
while invalid: #obtain wells to be tested
   well = input("Enter a well: ")
   if well[0] not in cols or well.lstrip("ABCDEFGH") not in rows: #error check if invalid well is entered
       print("Error, please try again")
   else:
       invalid = False

well_data = []
#print(cleaned_data)
for i in range(0, len(cleaned_data)):
    if cleaned_data[i][0] == well:
        values = cleaned_data[i]
        well_data.append(values)
#print(well_data)
#print("\n")
if len(well_data) == 0:
    print("Well was not tested")
    sys.exit()

times = []
e_mods = []
for j in range(0, len(well_data)):
    date_time = well_data[j][3].split(" ")
    time = date_time[1]
    print(f"At time {time}, the elastic modulus of well {well_data[j][0]} was {well_data[j][1]} N/m^2 with an uncertainty of {well_data[j][2]} N/m^2")
    times.append(time)
    e_mods.append(well_data[j][1])
    if (j % 52) - 50 == 1:
        bad_answer = True
        while bad_answer:
            cont = input("The maximum number of lines that can fit on the screen has been reached, but there is still more data. Please input Y when you are ready to continue: ")
            if cont.strip() == 'y' or cont.strip() == 'Y':
                bad_answer = False
                print("Okay, displaying next datapoints")
                time.sleep(2)
            else:
                print("Invalid response, please try again")
pyplot.plot(times, e_mods)
pyplot.xlabel("Time")
pyplot.ylabel("Elastic Modulus (N/m^2)")
pyplot.title(f"Elastic modulus of well {well} over time")
pyplot.show()



