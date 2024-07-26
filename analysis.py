"""
This is the analysis script to display results from measurements made using the measure and custom_measure programs
"""

import csv
import sys
import time
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.optimize import curve_fit
import os

files = os.listdir("/home/robot/ASMI_KABlab") #change to correct directory for your device!!!


def load_csv(): #load data from csv file
   bad_name = True
   while bad_name: #select file to obtain data from
        filename = input(
            "Please enter the name of the the file you would like to analyze the data from. The name is case"
            " sensitive, so please enter the name in exactly. To see a list of all files, press Ctrl + C and"
            " type ls in the command window, then rerun this program. Note, please do not type in the .csv, just enter "
            "the name of the file. ")
        filename = filename + ".csv"
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
   return cleaned_data


def collect_run_data(data, well): #collect data for specific run from csv file
    well_data = []
    no_contact = []
    run_array = []
    forces = []
    for i in range(0, len(data)):
        if data[i][0] == well: #collect data from specified well
            values = [data[i][1], data[i][2]]
            well_data.append(values)
    #print(well_data)
    #print("\n")
    if len(well_data) == 0:
        print("Well was not tested")
        sys.exit()
    for l in range(1, len(well_data)):
        if float(well_data[l][1]) <= -1*float(well_data[0][0]) + 2*float(well_data[0][1]): #determine which measurements correspond to contact
            no_contact.append(l)
        run_array.append([well_data[l][0], well_data[l][1]])
    #print(run_array)
    #print("\n")
    #print(no_contact)
    #print("\n")
    if len(run_array) - int(no_contact[len(no_contact) - 1]) <= 10: #check if no or not enough data was collected for well
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        sys.exit()
    if len(no_contact) > 0: #find index of first continuous contact measurement
        start_val = int(no_contact[len(no_contact)-1]+1)
    else:
        start_val = 0
    #print(start_val)
    #print(len(well_data))
    #print(run_array[start_val][0])
    for k in range(0, len(run_array)): #format data for analysis
        run_array[k][0] = round(-1*(float(run_array[k][0]) - float(well_data[start_val][0])), 2) #set indentation depths relative to initial contact height
        run_array[k][1] = float(run_array[k][1]) + float(well_data[0][0]) #zero forces
        forces.append(run_array[k][1])
    if forces == [] or max(forces)-min(forces) < 0.04: #check that force measurements were large enough to make proper measurement
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        sys.exit()
    get_ratio = True
    while get_ratio: #obtain Poisson's ratios for specified sample
        p_ratio = input("What is the approximate Poisson's Ratio of the sample? Value should be between 0.3-0.5. ")
        cleaned_input = p_ratio.replace(".", "1")
        if cleaned_input.isnumeric() and float(p_ratio) >= 0.3 and float(p_ratio) <= 0.5:
            p_ratio = float(p_ratio)
            get_ratio = False
        else:
            print("Improper Poisson's Ratio, please try again.")
    return run_array, p_ratio

def approximate_height(run_array): #find height of sample to determine correction equation used
    depths = []
    for i in range(0, len(run_array)):
        depths.append(run_array[i][0])
    for j in range(0, len(depths)):
        depths[j] = abs(depths[j])
    zero = min(depths)
    num = depths.index(zero)
    z_pos = (num * 0.02) + 3
    approx_height = 15 - z_pos
    #print(approx_height)
    return approx_height

def split(run_array): #splits data from well into separate depth and force arrays
    depths = []
    forces = []
    for i in range(0, len(run_array)):
        depths.append(run_array[i][0])
        forces.append(run_array[i][1])
    return depths, forces

def find_d_and_f_in_range(run_array): #select data within desired depth range to determine elastic modulus
    forces = []
    depths = []
    for i in range(0, len(run_array)):
        if run_array[i][0] >= 0.24 and run_array[i][0] <= 0.5: #.04, .3
            forces.append(run_array[i][1])
            depths.append(run_array[i][0])
            #print(i)
    return depths, forces

def correct_force(depths, forces, p_ratio, approx_height): #add correction factor based on simulation data since samples are not ideal shapes
    new_array = []
    for i in range(0, len(depths)):
        if p_ratio < 0.325:
            if approx_height >= 9.5:
                b = 0.13
                c = 1.24
            elif approx_height >= 8.5 and approx_height < 9.5:
                b = 0.131
                c = 1.24
            elif approx_height >= 7.5 and approx_height < 8.5:
                b = 0.133
                c = 1.25
            elif approx_height >= 6.5 and approx_height < 7.5:
                b = 0.132
                c = 1.24
            elif approx_height >= 5.5 and approx_height < 6.5:
                b = 0.132
                c = 1.24
            elif approx_height >= 4.5 and approx_height < 5.5:
                b = 0.139
                c = 1.27
            elif approx_height >= 3.5 and approx_height < 4.5:
                b = 0.149
                c = 1.3
            else:
                b = 0.162
                c = 1.38
        elif p_ratio >= 0.325 and p_ratio < 0.375:
            if approx_height >= 9.5:
                b = 0.132
                c = 1.25
            elif approx_height >= 8.5 and approx_height < 9.5:
                b = 0.132
                c = 1.25
            elif approx_height >= 7.5 and approx_height < 8.5:
                b = 0.134
                c = 1.25
            elif approx_height >= 6.5 and approx_height < 7.5:
                b = 0.136
                c = 1.26
            elif approx_height >= 5.5 and approx_height < 6.5:
                b = 0.126
                c = 1.25
            elif approx_height >= 4.5 and approx_height < 5.5:
                b = 0.133
                c = 1.27
            elif approx_height >= 3.5 and approx_height < 4.5:
                b = 0.144
                c = 1.32
            else:
                b = 0.169
                c = 1.42
        elif p_ratio >= 0.375 and p_ratio < 0.425:
            if approx_height >= 9.5:
                b = 0.181
                c = 1.33
            elif approx_height >= 8.5 and approx_height < 9.5:
                b = 0.182
                c = 1.34
            elif approx_height >= 7.5 and approx_height < 8.5:
                b = 0.183
                c = 1.34
            elif approx_height >= 6.5 and approx_height < 7.5:
                b = 0.183
                c = 1.34
            elif approx_height >= 5.5 and approx_height < 6.5:
                b = 0.194
                c = 1.38
            elif approx_height >= 4.5 and approx_height < 5.5:
                b = 0.198
                c = 1.4
            elif approx_height >= 3.5 and approx_height < 4.5:
                b = 0.203
                c = 1.44
            else:
                b = 0.176
                c = 1.46
        elif p_ratio >= 0.425 and p_ratio < 0.475:
            if approx_height >= 9.5:
                b = 0.156
                c = 1.35
            elif approx_height >= 8.5 and approx_height < 9.5:
                b = 0.152
                c = 1.34
            elif approx_height >= 7.5 and approx_height < 8.5:
                b = 0.156
                c = 1.35
            elif approx_height >= 6.5 and approx_height < 7.5:
                b = 0.161
                c = 1.37
            elif approx_height >= 5.5 and approx_height < 6.5:
                b = 0.153
                c = 1.37
            elif approx_height >= 4.5 and approx_height < 5.5:
                b = 0.166
                c = 1.42
            elif approx_height >= 3.5 and approx_height < 4.5:
                b = 0.179
                c = 1.47
            else:
                b = 0.205
                c = 1.59
        else:
            if approx_height >= 9.5:
                b = 0.203
                c = 1.58
            elif approx_height >= 8.5 and approx_height < 9.5:
                b = 0.207
                c = 1.6
            elif approx_height >= 7.5 and approx_height < 8.5:
                b = 0.212
                c = 1.62
            elif approx_height >= 6.5 and approx_height < 7.5:
                b = 0.217
                c = 1.65
            elif approx_height >= 5.5 and approx_height < 6.5:
                b = 0.21
                c = 1.64
            elif approx_height >= 4.5 and approx_height < 5.5:
                b = 0.22
                c = 1.68
            elif approx_height >= 3.5 and approx_height < 4.5:
                b = 0.17
                c = 1.58
            else:
                b = 0.182
                c = 1.64
        val = (forces[i])/(c*pow(depths[i], b))
        new_array.append(val)
    return new_array


def adjust_depth(run_array, d0): #using curve fit, adjust depth so that at zero force, depth is 0
    for i in range(0, len(run_array)):
        run_array[i][0] = run_array[i][0]-d0
    return run_array


def find_E(A, p_ratio): #determine elastic modulus from curve fit
    r_sphere = 0.0025
    sphere_p_ratio = 0.28
    sphere_E = 1.8 * pow(10, 11)
    polymer_p_ratio = p_ratio
    actual_A = A * pow(1000, 1.5)
    E_star = (actual_A * 0.75)/pow(r_sphere, 0.5)
    E_inv = 1/(E_star * (1 - pow(polymer_p_ratio, 2))) - (1 - pow(sphere_p_ratio, 2))/(sphere_E * (1 - pow(polymer_p_ratio, 2)))
    E_polymer = 1/E_inv
    return E_polymer

def adjust_E(E): #an empirical correction factor for softer samples which causes issues with getting proper data at small indentation depths
    if E < 660000:
        factor = 457 * pow(E, -0.457)
        E = E/factor
    return E



data = load_csv()
#print(data)
cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
invalid = True
while invalid: #obtain well to be analyzed
   well = input("Enter a well: ")
   if well[0] not in cols or well.lstrip("ABCDEFGH") not in rows: #error check if invalid well is entered
       print("Error, please try again")
   else:
       invalid = False
run_array, p_ratio = collect_run_data(data, well)
well_data = run_array
#print(run_array)
depths, forces = split(run_array)
height = approximate_height(run_array)
#pyplot.scatter(depths, forces)
#pyplot.show()
#print(depths)
#print(forces)
well_depths = depths
well_forces = forces
depth_in_range, force_in_range = find_d_and_f_in_range(run_array)
#print(depth_in_range)
#print(force_in_range)
adjusted_forces = correct_force(depth_in_range, force_in_range, p_ratio, height)
#print(adjusted_forces)
depth_in_range = np.asarray(depth_in_range)
adjusted_forces = np.asarray(adjusted_forces)


def Hertz_func(depth, A, d0):
  F = A * pow(depth - d0, 1.5)
  return F



try: #fit data to Hertzian contact mechanics
    parameters, covariance = curve_fit(Hertz_func, depth_in_range, adjusted_forces, p0=[2, 0.03])
except:
    print("Data could not be analyzed")
    sys.exit()
else:

    #print(depth_in_range)
    #print(force_in_range)


    fit_A = float(parameters[0])
    fit_d0 = float(parameters[1])
    #print(fit_A)
    #print(fit_d0)
    #pyplot.scatter(depth_in_range, adjusted_forces)
    #y_var = []
    #for i in range(0, len(depth_in_range)):
        #y_var.append(fit_A * pow(depth_in_range[i], 1.5))
    #pyplot.plot(depth_in_range, y_var)
    #pyplot.xlabel("Depth (mm)")
    #pyplot.ylabel("Force (N)")
    #pyplot.title("Force vs. Indentation Depth")
    #pyplot.show()

    count = 0 #adjust if approximate initial depth was incorrect
    continue_to_adjust = True
    if abs(fit_d0) < 0.01:
        continue_to_adjust = False
    min_d0 = 100
    while continue_to_adjust:
        count = count + 1
        old_d0 = fit_d0
        run_array = adjust_depth(run_array, fit_d0)
        height = approximate_height(run_array)
        depth_in_range, force_in_range = find_d_and_f_in_range(run_array)
        #print(depth_in_range)
        adjusted_forces = correct_force(depth_in_range, force_in_range, p_ratio, height)
        #print(adjusted_forces)
        depth_in_range = np.asarray(depth_in_range)
        adjusted_forces = np.asarray(adjusted_forces)
        try:
            parameters, covariance = curve_fit(Hertz_func, depth_in_range, adjusted_forces, p0=[2, 0.03])
        except:
            print("Data could not be analyzed")
            sys.exit()
        else:
            fit_A = float(parameters[0])
            fit_d0 = float(parameters[1])
            #print(fit_A)
            #print(fit_d0)
            if abs(fit_d0) < min_d0:
                min_d0 = abs(fit_d0)
            #print(f"min {min_d0}")
            #pyplot.scatter(depth_in_range, adjusted_forces)
            #y_var = []
            #for i in range(0, len(depth_in_range)):
                #y_var.append(fit_A * pow(depth_in_range[i], 1.5))
            #pyplot.plot(depth_in_range, y_var)
            #pyplot.xlabel("Depth (mm)")
            #pyplot.ylabel("Force (N)")
            #pyplot.title("Force vs. Indentation Depth")
            #pyplot.show()
            if abs(round(old_d0, 5)) == abs(round(fit_d0, 5)): #if fit continues to converge to improper value
                fit_d0 = -0.75 * fit_d0
            elif abs(fit_d0) < 0.01:
                continue_to_adjust = False
                break
            elif count > 100 and count < 200:
                if abs(round(fit_d0, 2)) == round(min_d0, 2):
                    break
            elif count >= 200 and count < 300:
                if abs(round(fit_d0, 1)) == round(min_d0, 1):
                    break
            elif count == 300:
                print("Fit quality may be poor")
                break

    E = find_E(fit_A, p_ratio) #determine elastic modulus from measurements
    #print(E)
    E = adjust_E(E)
    E = round(E)
    ##print(E)
    if round(max(depth_in_range), 2) < 0.4:
        print("Sample was not indented far enough")
        print(
            f"The range the measurement was made with was {round(min(depth_in_range), 2)} mm to {round(max(depth_in_range), 2)} mm")
    err = np.sqrt(np.diag(covariance))
    #print(covariance[0][0])
    std_dev = round(find_E(err[0], p_ratio))
    ##print(std_dev)
    print(f"Well {well}: E = {E} N/m^2, Uncertainty = {std_dev} N/m^2")
    pyplot.scatter(depth_in_range, adjusted_forces)
    y_var = []
    for i in range(0, len(depth_in_range)):
       y_var.append(fit_A * pow(depth_in_range[i], 1.5))
    pyplot.plot(depth_in_range, y_var) #plot data and curve fit
    pyplot.xlabel("Depth (mm)")
    pyplot.ylabel("Force (N)")
    pyplot.title(f"Force vs. Indentation Depth of Well {well}")
    pyplot.show()
