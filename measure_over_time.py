import serial
import time
from threading import Event
from godirect import GoDirect
import statistics
import csv
import logging
import numpy as np
import math
import warnings
import os
from datetime import datetime

import matplotlib.pyplot as pyplot

warnings.filterwarnings("ignore", category=RuntimeWarning)
from scipy.optimize import curve_fit
logging.basicConfig()

BAUD_RATE = 115200
GRBL_port_path = "/dev/ttyUSB0" #"/dev/ttyUSB0"
x_init = 17.4  # set position of well A1
y_init = 51.1  # set position of well A1
offset = 8.95  # set distance between wells

godirect = GoDirect(use_ble=True, use_usb=True)
device = godirect.get_device(threshold=-100)
lowest = -11
height_offset = 4
files = os.listdir("/home/robot/ASMI_KABlab") #change to correct directory for your device!!!

def remove_comment(string):
    if (string.find(';') == -1):
        return string
    else:
        return string[:string.index(';')]


def remove_eol_chars(string):
    # removed \n or traling spaces
    return string.strip()


def send_wake_up(ser):
    # Wake up
    # Hit enter a few times to wake the Printrbot
    ser.write(str.encode("\r\n\r\n"))
    time.sleep(2)  # Wait for Printrbot to initialize
    ser.flushInput()  # Flush startup text in serial input


def wait_for_movement_completion(ser, cleaned_line):
    Event().wait(1)

    if cleaned_line != '$X' or '$$':

        idle_counter = 0

        while True:

            # Event().wait(0.01)
            ser.reset_input_buffer()
            command = str.encode('?' + '\n')
            ser.write(command)
            grbl_out = ser.readline()
            grbl_response = grbl_out.strip().decode('utf-8')
            #print(grbl_response)

            if grbl_response != 'ok':

                if grbl_response.find('Idle') > 0:
                    idle_counter += 1

            if idle_counter > 0:
                break
    return


def move_gcode(GRBL_port_path, gcode, home, x, y, z):
    # with contect opens file/connection and closes it if function(with) scope is left
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home:
            ##print(f'G92 X{x} Y{y} Z{z}\n')
            command = str.encode(f'G92 X{x} Y{y} Z{z}\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            x = 0
            y = 0
            z = 0
        if cleaned_line:  # checks if string is empty
            ##print("Sending gcode:" + str(cleaned_line))
                # converts string to byte encoded string and append newline
            command = str.encode(gcode + '\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]

        ##print('End of gcode')
        with open("position.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

def get_start_stats(well, filename, trial):
    measurements = []
    device.start()
    sensors = device.get_enabled_sensors()
    ##print("reading force")
    for i in range(0, 10):
        if device.read():
            for sensor in sensors:
                measurements.append(sensor.values[0])
                #print(measurements)
                ##print(sensor.values)
                sensor.clear()
        time.sleep(0.05)
    device.stop()
    average = statistics.mean(measurements)
    standard_dev = statistics.stdev(measurements)
    ##print(f"avg = {average}")
    ##print(f"stdev = {standard_dev}")
    with open(filename, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        row = [well, str(average), str(standard_dev), trial]
        csvwriter.writerow(row)
    return average, standard_dev

def get_measurement():
    device.start()
    val = 0
    sensors = device.get_enabled_sensors()
    if device.read():
        for sensor in sensors:
            val = sensor.values
            val = val[0]
            ##print(val)
            sensor.clear()
    device.stop()
    #val = val[0]
    val = val
    return val


def stream_gcode(GRBL_port_path, gcode, x, y, well, filename, trial):
    # with contect opens file/connection and closes it if function(with) scope is left
    stiff = False
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        avg, stdev = get_start_stats(well, filename, trial)
        down = True
        up = False
        z = -1*height_offset+1
        z_max = lowest - 2
        contact = False
        measurements = []
        for line in gcode:
            ##print(z)
            value = get_measurement()
            if value < avg - 2 * stdev:
                ##print("contact")
                ##print(z_max)
                measurements.append(value * -1)
                if z == z_max or value <= -20:
                    if value <= -20:
                        print("Sample too stiff to analyze")
                        stiff = True
                    return measurements, z, stiff
                if len(measurements) == 1:
                    ##print("creating z_max")
                    z_max = round(z - 1, 2)
                contact = True
                z = round(z-0.02, 2)
            elif contact == True and value >= avg - 2 * stdev:
                ##print("False alarm")
                z_max = lowest - 2
                contact = False
                measurements = []
                z = round(z-0.02, 2)
                #Start here!
            else:
                z = round(z-0.02, 2)
            with open(filename, 'a') as csvfile:
                csvwriter = csv.writer(csvfile)
                row = [str(well), str(z), str(value * -1), trial]
                csvwriter.writerow(row)
            cleaned_line = remove_eol_chars(remove_comment(line))
            if cleaned_line:  # checks if string is empty
                ##print("Sending gcode:" + str(cleaned_line))
                # converts string to byte encoded string and append newline
                command = str.encode(line + '\n')
                ser.write(command)  # Send g-code

                wait_for_movement_completion(ser, cleaned_line)

                grbl_out = ser.readline()  # Wait for response with carriage return
                ##print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]
            with open("position.csv", 'w') as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(position)

        return measurements, z, stiff

def load_csv(filename): #load data from csv file
   with open(filename, 'r') as file:
       reader = csv.reader(file)
       data = list(reader)
   cleaned_data = []
   for i in range(0, len(data)):
       if data[i] != []:
           cleaned_data.append(data[i])
   return cleaned_data


def collect_run_data(data, well, stiff, trial): #collect data for specific run from csv file
    well_data = []
    no_contact = []
    run_array = []
    forces = []
    if stiff:
        return run_array
    for i in range(0, len(data)):
        if data[i][0] == well and data[i][3] == trial:
            values = [data[i][1], data[i][2]]
            well_data.append(values)
    #print(well_data)
    #print("\n")
    for l in range(1, len(well_data)):
        if float(well_data[l][1]) <= -1*float(well_data[0][0]) + 2*float(well_data[0][1]):
            no_contact.append(l)
        run_array.append([well_data[l][0], well_data[l][1]])
    #print(run_array)
    #print("\n")
    #print(no_contact)
    #print("\n")
    if len(run_array) - int(no_contact[len(no_contact) - 1]) <= 10:
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        return run_array
    if len(no_contact) > 0:
        start_val = int(no_contact[len(no_contact)-1]+1)
    else:
        start_val = 0
    #print(start_val)
    #print(len(well_data))
    #print(run_array[start_val][0])
    for k in range(0, len(run_array)):
        run_array[k][0] = round(-1*(float(run_array[k][0]) - float(well_data[start_val][0])), 2)
        run_array[k][1] = float(run_array[k][1]) + float(well_data[0][0])
        forces.append(run_array[k][1])
    if forces == [] or max(forces)-min(forces) < 0.04:
        print("Either well was not tested or no data was collected, either because sample was too short or too soft")
        run_array = []
        return run_array
    #print(run_array)
    #print("\n")
    return run_array

def split(run_array):
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
    return depths, forces

def approximate_height(run_array):
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

def adjust_E(E):
    if E < 660000:
        factor = 457*pow(E, -0.457)
        E = E/factor
    return E


def go_home(GRBL_port_path):
    n = 1
    with open('position.csv', mode='r') as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            if n == 1:
                position = lines
                #print(position)
            n = n + 1
    #print(position)
    x = position[0]
    y = position[1]
    z = position[2]
    if x !=0 or y != 0 or z !=0:
        #print(x)
        #print(y)
        #print(z)
        home = True
        gcode = f"G01 Z0 F500"
        #print(f"G01 Z0 F500")
        home_z(GRBL_port_path, gcode, home, x, y, z)
        gcode = f"G01 X0 Y0 F500"
        #print(f"G01 X0 Y0 F500")
        home_xy(GRBL_port_path, gcode, home, x, y, z)

def home_z(GRBL_port_path, gcode, home, x, y, z):
    # with contect opens file/connection and closes it if function(with) scope is left
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home:
            #print(f'G92 Z{z}\n')
            command = str.encode(f'G92 Z{z}\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            #print(" : ", grbl_out.strip().decode('utf-8'))
            z = 0
        if cleaned_line:  # checks if string is empty
            #print("Sending gcode:" + str(cleaned_line))
                # converts string to byte encoded string and append newline
            command = str.encode(gcode + '\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            #print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]

        #print('End of gcode')
        with open("position.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

def home_xy(GRBL_port_path, gcode, home, x, y, z):
    # with contect opens file/connection and closes it if function(with) scope is left
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home:
            #print(f'G92 X{x} Y{y}\n')
            command = str.encode(f'G92 X{x} Y{y}\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            #print(" : ", grbl_out.strip().decode('utf-8'))
            x = 0
            y = 0
            z = 0
        if cleaned_line:  # checks if string is empty
            #print("Sending gcode:" + str(cleaned_line))
            # converts string to byte encoded string and append newline
            command = str.encode(gcode + '\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            #print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]

        #print('End of gcode')
        with open("position.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

if __name__ == "__main__":
    # GRBL_port_path = '/dev/tty.usbserial-A906L14X'
    #gcode_path = 'gcode.gcode'
    #serial.Serial(GRBL_port_path, BAUD_RATE).close()

    #print("USB Port: ", GRBL_port_path)
    #print("Gcode file: ", gcode_path)
    device.open(auto_start=False)

    f = open("measurements.csv", "w")
    f.truncate()
    f.close()

    go_home(GRBL_port_path)

    entry_wells = []
    cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
    rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    x = {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": "", "H": ""}
    y = {"1": "", "2": "", "3": "", "4": "", "5": "", "6": "", "7": "", "8": "", "9": "", "10": "", "11": "", "12": ""}
    z_up = "-2.50"
    height_offset = 4  # set starting distance between indenter and wells
    # y_disp = 0.1 #well plate is not precisely aligned, should get fixed in future iterations
    h_speed = "500"  # speed sensor moves between wells
    v_speed = "100"  # speed sensor moves while testing sample
    results = []

    for i in range(0, 8):  # load x values into x dictionairy
        x[cols[i]] = str(x_init + offset * i)

    for j in range(0, 12):  # load y values into y dictionairy
        y[rows[j]] = str(y_init + offset * j)


    bad_input = True
    while bad_input:
        num_tests = input("Please enter how many tests you would like to be run: ")
        time_between = input("Please enter how many minutes you would like between set of tests to the nearest whole number: ")
        if num_tests.isnumeric() and time_between.isnumeric():
            print(f"Valid inputs, a total of {num_tests} tests will be run with {time_between} minutes in between each test.")
            bad_input = False
        else:
            print("Invalid input, please try again")

    num_tests = int(num_tests)
    time_between = int(time_between)
    time_between = round(time_between * 60)

    bad_name = True
    while bad_name:
        filename = input(
            "Please enter the name of the the file you would like to save the data to. The format should have"
            " no spaces like the following: File_name ")
        count = 0
        for i in range(0, len(filename)):
            if filename[i] == " ":
                count = count + 1
        if count > 0:
            print("Improper filename, file contained a space, please try again")
            retry = True
        else:
            retry = False
        if not retry:
            results_filename = filename + "_results_timed.csv"
            filename = filename + ".csv"
            if filename in files:
                bad_answer = True
                while bad_answer:
                    overwrite = input("A file with that name already exists, would you like to overwrite it? Y or N: ")
                    if overwrite.strip() == 'y' or overwrite.strip() == 'Y':
                        bad_answer = False
                        bad_name = False
                        print("Okay, overwriting file")
                    elif overwrite.strip() == "n" or overwrite.strip() == "N":
                        bad_answer = False
                        print("Okay, restarting file naming process")
                    else:
                        print("Invalid response, please try again")
            else:
                bad_name = False
    # print(filename)
    # print(files)

    print("Wells can be entered in combinations of different modes. Please note wells will be tested in the order they are entered.")
    print("Also, wells that are entered more than once will be tested only once.")
    modes = ["1", "2", "3", "4"]
    entered = False
    while not entered:
        mode = input("Please enter 1 to enter individual wells, 2 to enter a row/rows, 3 to enter a column/columns,\n"
                     " or 4 to enter the entire plate: ")
        while mode not in modes:
            print("Invalid mode, please try again")
            mode = input("Please enter 1 to enter individual wells, 2 to enter a row/rows, 3 to enter a column/columns,\n"
                         " or 4 to enter the entire plate: ")
        if mode == "1":
            more = True
            while more:
                well = input("Enter a well: ")
                if well == '':  # exit condition if all wells are entered
                    more = False
                elif well[0] not in cols or well.lstrip("ABCDEFGH") not in rows:  # error check if invalid well is entered
                    print("Error, please try again")
                else:
                    entry_wells.append(well)  # if valid well is entered, add to wells to be tested
                    print("Valid well, press enter if all wells are entered")
        if mode == "2":
            more = True
            while more:
                row = input("Enter a row (\"1\", \"2\", ... : ")
                if row == '':  # exit condition if all wells are entered
                    more = False
                elif row not in rows:
                    print("Error, please try again")
                else:
                    for i in range(0, len(cols)):
                        well = cols[i] + row
                        entry_wells.append(well)
                    print("Valid row, press enter if all rows are entered")
        if mode == "3":
            more = True
            while more:
                col = input("Enter a column (\"A\", \"B\", ... : ")
                if col == '':  # exit condition if all wells are entered
                    more = False
                elif col not in cols:
                    print("Error, please try again")
                else:
                    for i in range(0, len(rows)):
                        well = col + rows[i]
                        entry_wells.append(well)
                    print("Valid column, press enter if all columns are entered")
        if mode == "4":
            for j in range(0, len(cols)):
                for i in range(0, len(rows)):
                    well = cols[j] + rows[i]
                    entry_wells.append(well)
                print("Valid entry")
        print(entry_wells)
        invalid_answer = True
        while invalid_answer:
            correct = input("Please confirm that these are the correct wells, \"Y\" or \"N\" ")
            #print(correct)
            if correct.strip() == 'y' or correct.strip() == 'Y':
                invalid_answer = False
            elif correct.strip() == "n" or correct.strip() == "N":
                invalid_answer = False
                entry_wells = []
                print("Okay, restarting well entry process")
            else:
                print("Invalid response, please try again")
        different_mode = False
        while not different_mode:
            correct = input("Are all wells entered?, \"Y\" or \"N\" ")
            if correct.strip() == "y" or correct.strip() == "Y":
                different_mode = True
                entered = True
            elif correct.strip() == "n" or correct.strip() == "N":
                different_mode = True
            else:
                print("Invalid response, please try again")

    wells = []
    for well in entry_wells:
        if well not in wells:
            wells.append(well)

    p_ratios = []
    invalid_answer = True
    while invalid_answer:
        same_ratio = input("Do all of the samples have approximately the same Poisson's Ratio?, \"Y\" or \"N\" ")
        # print(correct)
        if same_ratio.strip() == 'y' or same_ratio.strip() == 'Y':
            get_ratio = True
            while get_ratio:
                p_ratio = input(
                    "What is the approximate Poisson's Ratio of the samples? Value should be between 0.3-0.5. ")
                cleaned_input = p_ratio.replace(".", "1")
                if cleaned_input.isnumeric() and float(p_ratio) >= 0.3 and float(p_ratio) <= 0.5:
                    p_ratio = float(p_ratio)
                    for i in range(0, len(wells)):
                        p_ratios.append(p_ratio)
                    get_ratio = False
                else:
                    print("Improper Poisson's Ratio, please try again.")
            invalid_answer = False
        elif same_ratio.strip() == "n" or same_ratio.strip() == "N":
            print("Okay, for each well, please input the approximate Poisson's Ratio")
            for well in wells:
                print(f"What is the approximate Poisson's Ratio of well {well}? ")
                get_ratio = True
                while get_ratio:
                    p_ratio = input(
                        "What is the approximate Poisson's Ratio of the sample? Value should be between 0.3-0.5. ")
                    cleaned_input = p_ratio.replace(".", "1")
                    if cleaned_input.isnumeric() and float(p_ratio) >= 0.3 and float(p_ratio) <= 0.5:
                        p_ratio = float(p_ratio)
                        p_ratios.append(p_ratio)
                        get_ratio = False
                    else:
                        print("Improper Poisson's Ratio, please try again.")
            invalid_answer = False
        else:
            print("Invalid response, please try again")

    print(f"Okay, testing wells: {wells}")

    max_time = (len(wells) * 9.2 * 60) * num_tests + ((num_tests - 1) * time_between)
    hrs = math.floor(max_time / 3600)
    mins = round(max_time / 60) - hrs * 60
    print(f"Estimated time to run tests is {hrs} hours and {mins} minutes, though times may vary")

    for trial in range(0, num_tests):
        #Test machine homing
        curr_x = 0
        curr_y = 0
        home = False
        for n in range(0, len(wells)):
            max_time = (len(wells)-n) * 9.2 * 60
            hrs = math.floor(max_time / 3600)
            mins = round(max_time / 60) - hrs * 60
            print(f"Testing well {wells[n]}")
            #print(f"Estimated time remaining is {hrs} hours and {mins} minutes, though times may vary")
            well = wells[n]
            col = wells[n][0]
            X = x[col]
            X = float(X)
            #print(f"X: {X}")
            new_X = X-curr_x
            #print(f"new_X: {new_X}")
            new_X = str(new_X)
            ind = cols.index(col) + 1
            row = wells[n].lstrip("ABCDEFGH")
            Y = str(float(y[row]))
            Y = float(Y)
            #print(f"Y: {Y}")
            new_Y = Y-curr_y
            #print(f"new_Y: {new_Y}")
            new_Y = str(new_Y)
            #print(f"G01 X{new_X} Y{new_Y} F{h_speed}")
            z_int = round(-1*height_offset+1, 2)
            Z = z_int
            if n == 0:
                gcode = f"G01 X{new_X} Y{new_Y} Z{z_int} F{h_speed}"
            else:
                gcode = f"G01 X{new_X} Y{new_Y} F{h_speed}"
            move_gcode(GRBL_port_path, gcode, home, X, Y, Z)
            gcode = []
            z = -0.02
            while z >= lowest+1:
                gcode.append(f"G01 Z{z} F{v_speed}")
                z = round(z-0.02, 2)
            measurements, z, stiff = stream_gcode(GRBL_port_path, gcode, X, Y, well, filename, trial)
            gcode = f"G01 Z{-z+Z} F{v_speed}"
            move_gcode(GRBL_port_path, gcode, home, X, Y, Z)
            curr_x = X
            #print(f"curr_X: {curr_x}")
            curr_y = Y
            #print(f"curr_Y: {curr_y}")
            #Save measurements to csv, probably do in function so depth values can be added
            #Analysis
            data = load_csv(filename)
            # print(data)
            cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
            rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
            well = wells[n]
            run_array = collect_run_data(data, well, stiff, trial)
            if run_array != []:
                well_data = run_array
                # print(run_array)
                height = approximate_height(run_array)
                depths, forces = split(run_array)
                # pyplot.scatter(depths, forces)
                # pyplot.show()
                # print(depths)
                # print(forces)
                well_depths = depths
                well_forces = forces
                depth_in_range, force_in_range = find_d_and_f_in_range(run_array)
                # print(depth_in_range)
                # print(force_in_range)
                p_ratio = p_ratios[n]
                adjusted_forces = correct_force(depth_in_range, force_in_range, p_ratio, height)
                # print(adjusted_forces)
                depth_in_range = np.asarray(depth_in_range)
                adjusted_forces = np.asarray(adjusted_forces)


                def Hertz_func(depth, A, d0):
                    F = A * pow(depth - d0, 1.5)
                    return F

                try:
                    parameters, covariance = curve_fit(Hertz_func, depth_in_range, adjusted_forces, p0=[2, 0.03])
                except:
                    print("Data could not be analyzed")
                    error = True
                else:

                    # print(depth_in_range)
                    # print(force_in_range)

                    fit_A = float(parameters[0])
                    fit_d0 = float(parameters[1])
                    # print(fit_A)
                    # print(fit_d0)
                    # pyplot.scatter(depth_in_range, adjusted_forces)
                    # y_var = []
                    # for i in range(0, len(depth_in_range)):
                    # y_var.append(fit_A * pow(depth_in_range[i], 1.5))
                    # pyplot.plot(depth_in_range, y_var)
                    # pyplot.xlabel("Depth (mm)")
                    # pyplot.ylabel("Force (N)")
                    # pyplot.title("Force vs. Indentation Depth")
                    # pyplot.show()

                    count = 0
                    continue_to_adjust = True
                    if abs(fit_d0) < 0.01:
                        continue_to_adjust = False
                    min_d0 = 100
                    error = False
                    while continue_to_adjust:
                        count = count + 1
                        old_d0 = fit_d0
                        run_array = adjust_depth(run_array, fit_d0)
                        depth_in_range, force_in_range = find_d_and_f_in_range(run_array)
                        # print(depth_in_range)
                        height = approximate_height(run_array)
                        adjusted_forces = correct_force(depth_in_range, force_in_range, p_ratio, height)
                        # print(adjusted_forces)
                        depth_in_range = np.asarray(depth_in_range)
                        adjusted_forces = np.asarray(adjusted_forces)
                        try:
                            parameters, covariance = curve_fit(Hertz_func, depth_in_range, adjusted_forces, p0=[2, 0.03])
                        except:
                            print("Data could not be analyzed")
                            error = True
                        else:
                            fit_A = float(parameters[0])
                            fit_d0 = float(parameters[1])
                            # print(fit_A)
                            # print(fit_d0)
                            if abs(fit_d0) < min_d0:
                                min_d0 = abs(fit_d0)
                            # print(f"min {min_d0}")
                            # pyplot.scatter(depth_in_range, adjusted_forces)
                            # y_var = []
                            # for i in range(0, len(depth_in_range)):
                            # y_var.append(fit_A * pow(depth_in_range[i], 1.5))
                            # pyplot.plot(depth_in_range, y_var)
                            # pyplot.xlabel("Depth (mm)")
                            # pyplot.ylabel("Force (N)")
                            # pyplot.title("Force vs. Indentation Depth")
                            # pyplot.show()
                            if abs(round(old_d0, 5)) == abs(round(fit_d0, 5)):
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
                                print("Error in data analysis")
                                error = True
                                break

                if not error:
                    E = find_E(fit_A, p_ratio)
                    #print(E)
                    E = adjust_E(E)
                    ##print(E)
                    E = round(E)
                    if round(max(depth_in_range), 2) < 0.4:
                        print("Sample was not indented far enough")
                        print(f"The range the measurement was made with was {round(min(depth_in_range), 2)} mm to {round(max(depth_in_range), 2)} mm")
                    err = np.sqrt(np.diag(covariance))
                    # print(covariance[0][0])
                    std_dev = round(find_E(err[0], p_ratio))
                    ##(std_dev)
                    now = datetime.now()
                    row = [wells[n], E, std_dev, now]
                    results.append(row)
                    print(f"Well {wells[n]}: E = {E} N/m^2, Uncertainty = {std_dev} N/m^2 at {now}")
                    with open(results_filename, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        csvwriter.writerow(row)
                else:
                    now = datetime.now()
                    row = [wells[n], "no data", "no data", now]
                    results.append(row)
                    with open(results_filename, 'a') as csvfile:
                        csvwriter = csv.writer(csvfile)
                        csvwriter.writerow(row)
            else:
                now = datetime.now()
                row = [wells[n], "no data", "no data", now]
                results.append(row)
                with open(results_filename, 'a') as csvfile:
                    csvwriter = csv.writer(csvfile)
                    csvwriter.writerow(row)


        home = True
        gcode = f"G01 X0 Y0 Z0 F500"
        ##print(f"G01 X0 Y0 Z0 F500")
        move_gcode(GRBL_port_path, gcode, home, round(curr_x, 2), round(curr_y, 2), round(Z, 2))
        #print("Here are the results:")
        #for l in range(0, len(results)):
            #print(f"Well {results[l][0]}: E = {results[l][1]} N/m^2, Uncertainty = {results[l][2]} N/m^2")
        print("Testing in progress! Please do not touch anything!")
        time.sleep(time_between)

    print("Testing is complete! Please us the analysis_over_time program to view the results for each well.")




