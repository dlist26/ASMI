"""
This is a script which moves the cnc to the location of inputted wells. This is a good check to make sure that the cnc
is set up and connected properly


Part of this code are adapted from the python_to_GRBL github repository found here:

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
import csv

BAUD_RATE = 115200
GRBL_port_path = "/dev/ttyUSB0" #Change this to the desired serial port!
x_init = 17.4  # set x position of well A1, change to fit to your device!
y_init = 51.1  # set y position of well A1, change to fit to your device!
offset = 8.95  # set distance between wells, change to fit to your device!


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

            if grbl_response != 'ok':

                if grbl_response.find('Idle') > 0:
                    idle_counter += 1

            if idle_counter > 10:
                break
    return


def stream_gcode(GRBL_port_path, gcode, home, x, y):
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home: #if device is being moved home, run this to reset coordinate system
            ##print(f'G92 X{x} Y{y} Z{z}\n')
            command = str.encode(f'G92 X{x} Y{y} \n')
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


if __name__ == "__main__":
    # GRBL_port_path = '/dev/tty.usbserial-A906L14X'
    #gcode_path = 'gcode.gcode'
    #serial.Serial(GRBL_port_path, BAUD_RATE).close()

    #print("USB Port: ", GRBL_port_path)
    #print("Gcode file: ", gcode_path)
    wells = []
    cols = ["A", "B", "C", "D", "E", "F", "G", "H"]
    rows = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    x = {"A": "", "B": "", "C": "", "D": "", "E": "", "F": "", "G": "", "H": ""}
    y = {"1": "", "2": "", "3": "", "4": "", "5": "", "6": "", "7": "", "8": "", "9": "", "10": "", "11": "", "12": ""}
    z_up = "-2.50"
    height_offset = 1.5  # set starting distance between indenter and wells
    # y_disp = 0.1 #well plate is not precisely aligned, should get fixed in future iterations
    h_speed = "500"  # speed sensor moves between wells
    v_speed = "10"  # speed sensor moves while testing sample

    for i in range(0, 8):  # load x values into x dictionary
        x[cols[i]] = str(x_init + offset * i)

    for j in range(0, 12):  # load y values into y dictionary
        y[rows[j]] = str(y_init + offset * j)

    f = open("wells.csv", "w")
    f.truncate()
    f.close()
    more = True
    while more:
        well = input("Enter a well: ")
        if well == '':  # exit condition if all wells are entered
            more = False
        elif well[0] not in cols or well.lstrip("ABCDEFGH") not in rows:  # error check if invalid well is entered
            print("Error, please try again")
        else:
            wells.append(well)  # if valid well is entered, add to wells to be tested
            print("Valid well, press enter if all wells are entered")

    print(f"Okay, moving to wells {wells}")
    curr_x = 0
    curr_y = 0
    curr_z = 0
    home = False
    for n in range(0, len(wells)):
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
        gcode = f"G01 X{new_X} Y{new_Y} F{h_speed}"
        stream_gcode(GRBL_port_path, gcode, home, 0, 0) #moves cnc to x and y coordinates of current well being tested
        curr_x = X
        #print(f"curr_X: {curr_x}")
        curr_y = Y
        #print(f"curr_Y: {curr_y}")

    home = True #move machine home
    gcode = f"G01 X0 Y0 F500"
    #print(f"G01 X0 Y0 F500")
    stream_gcode(GRBL_port_path, gcode, home, curr_x, curr_y)

    #with open("gcode.gcode", "w") as f:
     #   f.truncate()
      #  print(f"G91 X{str(-curr_x)} Y{str(-curr_y)}", file=f)
       # f.close()
    #stream_gcode(GRBL_port_path, gcode_path)


    #print('EOF')