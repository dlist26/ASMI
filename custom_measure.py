"""
This is a script which can be used to make measurements if samples are in non-standard well plate with 96 wells. It allows
the user to determine the location of three wells in order to map out every well before testing begins.


Parts of this code are adapted from the godirect-examples github repository found here:

https://github.com/VernierST/godirect-examples/tree/main
Copyright (c) 2018, Vernier Software and Technology, @VernierST

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This acknowledgement does not constitute an endorsement or promotion


Part of this code are also adapted from the python_to_GRBL github repository found here:

https://github.com/Sam-Freitas/python_to_GRBL/tree/main
Copyright (c) 2016 Florent Gallaire <fgallaire@gmail.com>
Copyright (c) 2014 Gabriele Cirulli

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

"""


import serial
import time
from threading import Event
from godirect import GoDirect
import statistics
import csv
import logging
import numpy as np
import math

import matplotlib.pyplot as pyplot
from scipy.optimize import curve_fit
logging.basicConfig()

BAUD_RATE = 115200
GRBL_port_path = "/dev/ttyUSB0" #Change this to the desired serial port!

godirect = GoDirect(use_ble=True, use_usb=True)
device = godirect.get_device(threshold=-100)
lowest = -11
height_offset = 4

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
    # Hit enter a few times to wake the cnc
    ser.write(str.encode("\r\n\r\n"))
    time.sleep(2)  # Wait for cnc to initialize
    ser.flushInput()  # Flush startup text in serial input


def wait_for_movement_completion(ser, cleaned_line): #wait for cnc to reach destination before sending new movement
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
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home: #if device is being moved home, run this to reset coordinate system
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
        with open("position.csv", 'w') as csvfile: #update position of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

def get_start_stats(well): #take measurements with no contact above each well so force measurements can be zeroed
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
    with open("measurements.csv", 'a') as csvfile: #first entry for each well is average non-contact force
        csvwriter = csv.writer(csvfile)
        row = [well, str(average), str(standard_dev)]
        csvwriter.writerow(row)
    return average, standard_dev

def get_measurement(): #takes 1 measuremnt
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


def stream_gcode(GRBL_port_path, gcode, x, y, well): #moves device based multiple lines of g-code to indent well
    stiff = False
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        avg, stdev = get_start_stats(well) #get non contact measurements for the current well
        down = True
        up = False
        z = -1*height_offset+1 #sets starting z-height of indenter
        z_max = lowest - 2 #sets the maximum z-height the indenter will go
        contact = False
        measurements = []
        for line in gcode: #for each line of gcode:
            ##print(z)
            value = get_measurement() #take a measurement
            if value < avg - 2 * stdev: #determine if contact was made based on non-contact measurements
                ##print("contact")
                ##print(z_max)
                measurements.append(value * -1) #multipy by -1 since pushing sensor is defined as a negative force
                if z == z_max or value <= -20: #exit condition if sample was indented far enough or force is too high that it might damage the indenter or sensor
                    if value <= -20:
                        print("Sample too stiff to analyze")
                        stiff = True
                    return measurements, z, stiff
                if len(measurements) == 1: #sets lowest z-height indenter needs to go to if contact was made
                    ##print("creating z_max")
                    z_max = round(z - 1, 2)
                contact = True
                z = round(z-0.02, 2)
            elif contact == True and value >= avg - 2 * stdev: #if contact was determined bu not sustained due to random fluctuations
                ##print("False alarm")
                z_max = lowest - 2 #reset lowest z height
                contact = False
                measurements = []
                z = round(z-0.02, 2)
            else:
                z = round(z-0.02, 2)
            with open("measurements.csv", 'a') as csvfile: #add data to csv file
                csvwriter = csv.writer(csvfile)
                row = [str(well), str(z), str(value * -1)]
                csvwriter.writerow(row)
            cleaned_line = remove_eol_chars(remove_comment(line))
            if cleaned_line:  # checks if string is empty, moves cnc
                ##print("Sending gcode:" + str(cleaned_line))
                # converts string to byte encoded string and append newline
                command = str.encode(line + '\n')
                ser.write(command)  # Send g-code

                wait_for_movement_completion(ser, cleaned_line)

                grbl_out = ser.readline()  # Wait for response with carriage return
                ##print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]
            with open("position.csv", 'w') as csvfile: #update the position of the cnc
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(position)

        return measurements, z

def load_csv(): #load data from csv file
   with open('measurements.csv', 'r') as file:
       reader = csv.reader(file)
       data = list(reader)
   cleaned_data = []
   for i in range(0, len(data)):
       if data[i] != []:
           cleaned_data.append(data[i])
   return cleaned_data


def collect_run_data(data, well, stiff): #collect data for specific run from csv file
    well_data = []
    no_contact = []
    run_array = []
    forces = []
    if stiff:
        return run_array
    for i in range(0, len(data)): #first, gather the data from only the specified well
        if data[i][0] == well:
            values = [data[i][1], data[i][2]]
            well_data.append(values)
    #print(well_data)
    #print("\n")
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
        return run_array
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
        return run_array
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

def adjust_E(E): #an empirical correction factor for softer samples which causes issues with getting proper data at small indentation depths
    if E < 660000:
        factor = 457 * pow(E, -0.457)
        E = E / factor
    return E


def go_home(GRBL_port_path): #ensures device starts at its home position, will not work is device is moved externally
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
    #print(x)
    #print(y)
    #print(z)
    home = True
    gcode = f"G01 Z0 F500"
    #print(f"G01 Z0 F500")
    home_z(GRBL_port_path, gcode, home, x, y, z) #move z up first to prevent accidental contact
    gcode = f"G01 X0 Y0 F500"
    #print(f"G01 X0 Y0 F500")
    home_xy(GRBL_port_path, gcode, home, x, y, z) #then move x and y to home position

def home_z(GRBL_port_path, gcode, home, x, y, z):
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
        with open("position.csv", 'w') as csvfile: #update position of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

def home_xy(GRBL_port_path, gcode, home, x, y, z):
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
        with open("position.csv", 'w') as csvfile: #update position of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

if __name__ == "__main__":
    # GRBL_port_path = '/dev/tty.usbserial-A906L14X'
    print("This function is for nonstandard well plates or if the sensor continues to run into the walls of the plate.")
    print("Move in the x direction by using the a (move sensor left) and d (move sensor right) keys")
    print("Move in the y direction by using the w (move tray forward) and s (move tray back) keys")
    print("Move in the z direction by using the r (move sensor down) and f (move sensor up) keys")
    print("Specify the distance of steps before each set of movements between .01 mm and 100 mm")
    print("Change step size by inputting exit")
    print("For reference, a standard well plate is approximately 85 mm across in the x direction")
    print("To start, please jog the machine so that the indenter sits exactly centered above the bottom left well A1")
    #gcode_path = 'gcode.gcode'
    #serial.Serial(GRBL_port_path, BAUD_RATE).close()

    #print("USB Port: ", GRBL_port_path)
    #print("Gcode file: ", gcode_path)
    device.open(auto_start=False)

    x_custom = 0
    y_custom = 0
    z_custom = 0
    entry_wells = []
    cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
    rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    x = {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": "", "H": ""}
    y = {"1": "", "2": "", "3": "", "4": "", "5": "", "6": "", "7": "", "8": "", "9": "", "10": "", "11": "", "12": ""}
    z_up = "-2.50"
    height_offset = 4  # set starting distance between indenter and wells
    h_speed = "500"  # speed sensor moves between wells
    v_speed = "100"  # speed sensor moves while testing sample
    x_init = -100  # set position of well A1
    y_init = -100  # set position of well A1
    x_offset = -100  # set distance between wells
    y_offset = -100

    set_coordinates = True
    while set_coordinates: #initiates coordinate setting process
        if x_init == -100:
            print("Please center the indenter above well A1 (bottom left corner)")
        elif x_offset == -100:
            print("Please center the indenter above well H1 (bottom right corner)")
        elif y_offset == -100:
            print("Please center the indenter above well H12 (top right corner)")
        step_size = input("Please enter the size of the step in mm for the next series of movements: ")
        if step_size.replace(".", "").isnumeric() and float(step_size) <= 100 and float(step_size) > 0.01:
            step_size = float(step_size) #sets size of step
            print(f"Valid step size, moving in increments of {step_size} mm")
            jog = True
            while jog: #moves device in desired direction by step size
                direction = input("Please move the device (a, d, w, s, r, f, or exit): ")
                if direction == "a":
                    if x_custom - step_size >= 0:
                        x_custom = x_custom - step_size
                        gcode = f"G01 X{-1*step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "d":
                    if x_custom + step_size <= 200:
                        x_custom = x_custom + step_size
                        gcode = f"G01 X{step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "w":
                    if y_custom + step_size <= 150:
                        y_custom = y_custom + step_size
                        gcode = f"G01 Y{step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "s":
                    if y_custom - step_size >= 0:
                        y_custom = y_custom - step_size
                        gcode = f"G01 Y{-1*step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "r":
                    if z_custom + step_size <= 0:
                        z_custom = z_custom + step_size
                        gcode = f"G01 Z{step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "f":
                    if z_custom - step_size >= -15:
                        z_custom = z_custom - step_size
                        gcode = f"G01 Z{-1*step_size} F{h_speed}"
                        move_gcode(GRBL_port_path, gcode, False, x_custom, y_custom, z_custom)
                    else:
                        print("Movement would exceed bound of machine, please try again")
                elif direction == "exit": #if finished making movements that are the specified size, new size will be asked if not in desired location
                    if x_init == -100: #first sets position of first well
                        invalid = True
                        while invalid:
                            positioned = input("Is the sensor centered above well A1? (Y or N): ")
                            if positioned == 'Y' or positioned == 'y':
                                x_init = x_custom
                                y_init = y_custom
                                print("Okay, setting initial x and y coordinates")
                                #print(x_init)
                                #print(y_init)
                                invalid = False
                            elif positioned == 'N' or positioned == 'n':
                                invalid = False
                            else:
                                print("Invalid input, please try again")
                    elif x_offset == -100: #next determines offset in x direction
                        invalid = True
                        while invalid:
                            positioned = input("Is the sensor centered above well H1? (Y or N): ")
                            if positioned == 'Y' or positioned == 'y':
                                x_offset = round((x_custom - x_init)/7, 2)
                                print("Okay, setting offsets in x direction")
                                #print(x_offset)
                                invalid = False
                            elif positioned == 'N' or positioned == 'n':
                                invalid = False
                            else:
                                print("Invalid input, please try again")
                    elif y_offset == -100: #lastly determines offset in y direction
                        invalid = True
                        while invalid:
                            positioned = input("Is the sensor centered above well H12? (Y or N): ")
                            if positioned == 'Y' or positioned == 'y':
                                y_offset = round((y_custom - y_init)/11, 2)
                                print("Okay, setting offsets in y direction")
                                #print(y_offset)
                                invalid = False
                                set_coordinates = False
                                break
                            elif positioned == 'N' or positioned == 'n':
                                invalid = False
                            else:
                                print("Invalid input, please try again")
                    jog = False
                else:
                    print("Invalid input, please try again")
            #print(set_coordinates)
            if y_offset != -100:
                break
        else:
            print("Invalid step size, please try again")

    print("Thank you for setting coordinates, moving device home to begin taking measurements")

    #print(x_init)
    #print(y_init)
    #print(x_offset)
    #print(y_offset)

    f = open("measurements.csv", "w")
    f.truncate()
    f.close()

    go_home(GRBL_port_path)

    results = []

    for i in range(0, 8):  # load x values into x dictionary
        x[cols[i]] = str(x_init + x_offset * i)

    for j in range(0, 12):  # load y values into y dictionary
        y[rows[j]] = str(y_init + y_offset * j)

    print("Wells can be entered in combinations of different modes. Please note wells will be tested in the order they are entered.")
    print("Also, wells that are entered more than once will be tested only once.")
    modes = ["1", "2", "3", "4"]
    entered = False
    while not entered: #allows wells to be entered in a combination of different modes
        mode = input("Please enter 1 to enter individual wells, 2 to enter a row/rows, 3 to enter a column/columns,\n"
                     " or 4 to enter the entire plate: ")
        while mode not in modes:
            print("Invalid mode, please try again")
            mode = input("Please enter 1 to enter individual wells, 2 to enter a row/rows, 3 to enter a column/columns,\n"
                         " or 4 to enter the entire plate: ")
        if mode == "1": #allows wells to be entered one at a time
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
        if mode == "2": #allows entire rows of wells to be entered at once
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
        if mode == "3": #allows entire columns of wells to be entered at once
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
                    print("Valid column, press enter if all column are entered")
            if mode == "4": #allows every well in the well plate to be entered at once
                for j in range(0, len(cols)):
                    for i in range(0, len(rows)):
                        well = cols[j] + rows[i]
                        entry_wells.append(well)
                    print("Valid entry")
        print(entry_wells)
        invalid_answer = True
        while invalid_answer: #allows user to check that the correct wells will be tested
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
        while not different_mode: #checks if more wells need to be entered
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

    #Test machine homing
    curr_x = 0
    curr_y = 0
    home = False
    for n in range(0, len(wells)):
        max_time = (len(wells)-n) * 9.2 * 60
        hrs = math.floor(max_time / 3600)
        mins = round(max_time / 60) - hrs * 60
        print(f"Testing well {wells[n]}")
        print(f"Estimated time remaining is {hrs} hours and {mins} minutes, though times may vary")
        well = wells[n]
        col = wells[n][0]
        X = x[col] #determines coordinates machine needs to move to in x and y direction to reach certain well
        X = float(X)
        #print(f"X: {X}")
        new_X = X-curr_x #determines distance cnc need to move to reach that well
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
        move_gcode(GRBL_port_path, gcode, home, X, Y, Z) #moves cnc to x and y coordinates of current well being tested
        gcode = []
        z = -0.02
        while z >= lowest+1: #creates lines of gcode to indent sample
            gcode.append(f"G01 Z{z} F{v_speed}")
            z = round(z-0.02, 2)
        measurements, z, stiff = stream_gcode(GRBL_port_path, gcode, X, Y, well)
        gcode = f"G01 Z{-z+Z} F{v_speed}"
        move_gcode(GRBL_port_path, gcode, home, X, Y, Z) #moves indenter back up to original z height
        curr_x = X
        #print(f"curr_X: {curr_x}")
        curr_y = Y
        #print(f"curr_Y: {curr_y}")
        #Save measurements to csv, probably do in function so depth values can be added
        #Analysis
        data = load_csv() #analyzes data from last well that was tested
        # print(data)
        cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
        rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
        well = wells[n]
        run_array = collect_run_data(data, well, stiff) #get data for specified well
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

            try: #fit data to a Hertzian contact mechanics function
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
                continue_to_adjust = True #if improper starting depth was determined using only force data, adjust based on curve fit
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
                    try: #refit to curve
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
                            break

            if not error:
                E = find_E(fit_A) #find elastic modulus
                #print(E)
                E = adjust_E(E) #adjustemnt based on empirical data for softer samples which may cause difficultly in making initial measurements
                ##print(E)
                E = round(E)
                if round(max(depth_in_range), 2) < 0.4:
                    print("Sample was not indented far enough")
                    print(f"The range the measurement was made with was {round(min(depth_in_range), 2)} mm to {round(max(depth_in_range), 2)} mm")
                err = np.sqrt(np.diag(covariance))
                # print(covariance[0][0])
                std_dev = round(find_E(err[0]))
                #print(std_dev)
                row = [wells[n], E, std_dev]
                results.append(row)
                print(f"Well {wells[n]}: E = {E} N/m^2, Uncertainty = {std_dev} N/m^2")
            else:
                row = [wells[n], "no data", "no data"]
                results.append(row)
        else:
            row = [wells[n], "no data", "no data"]
            results.append(row)


    home = True #move machine home
    gcode = f"G01 X0 Y0 Z0 F500"
    #print(f"G01 X0 Y0 Z0 F500")
    move_gcode(GRBL_port_path, gcode, home, round(curr_x, 2), round(curr_y, 2), round(Z, 2))
    print("Here are the results:")
    for l in range(0, len(results)):
        print(f"Well {results[l][0]}: E = {results[l][1]} N/m^2, Uncertainty = {results[l][2]} N/m^2")




