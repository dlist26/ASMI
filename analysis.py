"""
This is a script to analyze one specific well from the last time either the measure or custom_measure scripts
were used.

"""



import csv
import sys
import time
import numpy as np
import matplotlib.pyplot as pyplot
from scipy.optimize import curve_fit




def load_csv(): #load data from csv file containing data from only the most recent round of testing
   with open('measurements.csv', 'r') as file:
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
    for i in range(0, len(data)): #first, gather the data from only the specified well
        if data[i][0] == well:
            values = [data[i][1], data[i][2]]
            well_data.append(values)
    #print(well_data)
    #print("\n")
    if len(well_data) == 0:
        print("Well was not tested")
        sys.exit()
    for l in range(1, len(well_data)): #determine at what z-height first contact was made
        if float(well_data[l][1]) <= -1*float(well_data[0][0]) + 2*float(well_data[0][1]):
            no_contact.append(l)
        run_array.append([well_data[l][0], well_data[l][1]])
    #print(run_array)
    #print("\n")
    #print(no_contact)
    #print("\n")
    if len(run_array) - int(no_contact[len(no_contact) - 1]) <= 10: #do not analyze if not enough measurements were made
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        sys.exit()
    if len(no_contact) > 0:
        start_val = int(no_contact[len(no_contact)-1]+1)
    else:
        start_val = 0
    #print(start_val)
    #print(len(well_data))
    #print(run_array[start_val][0])
    for k in range(0, len(run_array)): #adjust forces by average "zeroed" force and adjust z heights to represent indentation depth
        run_array[k][0] = round(-1*(float(run_array[k][0]) - float(well_data[start_val][0])), 2)
        run_array[k][1] = float(run_array[k][1]) + float(well_data[0][0])
        forces.append(run_array[k][1])
    if forces == [] or max(forces)-min(forces) < 0.04: #do not analyze if sample is too soft
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        sys.exit()
    #print(run_array)
    #print("\n")
    return run_array

def split(run_array): #get force and depth data in separate arrays to adjust individually
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

def correct_force(depths, forces): #add correction factor based on simulation data since samples are not ideal shapes
    new_array = []
    for i in range(0, len(depths)):
        val = (forces[i])/(1.5364*pow(depths[i], 0.1113))
        new_array.append(val)
    return new_array


def adjust_depth(run_array, d0): #using curve fit, adjust depth so that at zero force, depth is 0
    for i in range(0, len(run_array)):
        run_array[i][0] = run_array[i][0]-d0
    return run_array


def find_E(A): #determine elastic modulus from curve fit
    r_sphere = 0.0025
    sphere_p_ratio = 0.28
    sphere_E = 1.8 * pow(10, 11)
    polymer_p_ratio = 0.5
    actual_A = A * pow(1000, 1.5)
    E_star = (actual_A * 0.75)/pow(r_sphere, 0.5)
    E_inv = 1/(E_star * (1 - pow(polymer_p_ratio, 2))) - (1 - pow(sphere_p_ratio, 2))/(sphere_E * (1 - pow(polymer_p_ratio, 2)))
    E_polymer = 1/E_inv
    return E_polymer

#def adjust_E(E): #a correction factor based on calibrated data that can be uncommented and adjusted to your own needs
    #factor = 457*pow(E, -0.457)
    #E = E/factor
    #return E



data = load_csv()
#print(data)
cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
invalid = True
while invalid: #obtain wells to be tested
   well = input("Enter a well: ")
   if well[0] not in cols or well.lstrip("ABCDEFGH") not in rows: #error check if invalid well is entered
       print("Error, please try again")
   else:
       invalid = False
run_array = collect_run_data(data, well) #get data for specified well
well_data = run_array
#print(run_array)
depths, forces = split(run_array)
#pyplot.scatter(depths, forces)
#pyplot.show()
#print(depths)
#print(forces)
well_depths = depths
well_forces = forces
depth_in_range, force_in_range = find_d_and_f_in_range(run_array) #find desired depths and forces at those depths
#print(depth_in_range)
#print(force_in_range)
adjusted_forces = correct_force(depth_in_range, force_in_range) #correct force based on non ideal shape
#print(adjusted_forces)
depth_in_range = np.asarray(depth_in_range)
adjusted_forces = np.asarray(adjusted_forces)


def Hertz_func(depth, A, d0):
  F = A * pow(depth - d0, 1.5)
  return F



try: #fit data to a Hertzian contact mechanics function
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

    count = 0
    continue_to_adjust = True #if improper starting depth was determined using only force data, adjust based on curve fit
    if abs(fit_d0) < 0.01:
        continue_to_adjust = False
    min_d0 = 100
    while continue_to_adjust:
        count = count + 1
        old_d0 = fit_d0
        run_array = adjust_depth(run_array, fit_d0) #adjust depths by curve fit offset
        depth_in_range, force_in_range = find_d_and_f_in_range(run_array) #find new forces and depth in the desired range
        #print(depth_in_range)
        adjusted_forces = correct_force(depth_in_range, force_in_range) #adjust original force values in range by updated depths
        #print(adjusted_forces)
        depth_in_range = np.asarray(depth_in_range)
        adjusted_forces = np.asarray(adjusted_forces)
        try: #refit to curve
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
            if abs(round(old_d0, 5)) == abs(round(fit_d0, 5)): #sometimes curve fit converges to improper value, nudges it away from convergence
                fit_d0 = -0.75 * fit_d0
            elif abs(fit_d0) < 0.01:
                continue_to_adjust = False
                break
            elif count > 100 and count < 200: #if optimal value is greater than 0.01, find the lowest obtained optimal value
                if abs(round(fit_d0, 2)) == round(min_d0, 2):
                    break
            elif count >= 200 and count < 300:
                if abs(round(fit_d0, 1)) == round(min_d0, 1):
                    break
            elif count == 300:
                fit_d0 = min_d0
                print("Possible error in data analysis")
                print(f"Optimal offset depth was {fit_d0}, normally it is < 0.01")
                print("if data is plotted, curve fit will likely not match data, but it is based on data that was used from a previous adjustment")
                break

    E = find_E(fit_A) #find elastic modulus
    #print(E)
    #if E < 2000000: #uncomment if you need to adjust calculated value based on calobration data
       #E = adjust_E(E)
    E = round(E)
    ##print(E)
    if round(max(depth_in_range), 2) < 0.4:
        print("Sample was not indented far enough")
        print(
            f"The range the measurement was made with was {round(min(depth_in_range), 2)} mm to {round(max(depth_in_range), 2)} mm")
    err = np.sqrt(np.diag(covariance))
    #print(covariance[0][0])
    std_dev = round(find_E(err[0]))
    ##print(std_dev)
    print(f"Well {well}: E = {E} N/m^2, Uncertainty = {std_dev} N/m^2")
    pyplot.scatter(depth_in_range, adjusted_forces)
    y_var = []
    for i in range(0, len(depth_in_range)):
       y_var.append(fit_A * pow(depth_in_range[i], 1.5))
    pyplot.plot(depth_in_range, y_var)
    pyplot.xlabel("Depth (mm)")
    pyplot.ylabel("Force (N)")
    pyplot.title(f"Force vs. Indentation Depth")# of Well {well}")
    pyplot.show()
