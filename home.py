"""
This is a script which moves the cnc back to its home position. This will only work properly if the device has not
been moved externally


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


def home_z(GRBL_port_path, gcode, home, x, y, z):
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home: #reset z coordinate system of device
            ##print(f'G92 Z{z}\n')
            command = str.encode(f'G92 Z{z}\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            z = 0
        if cleaned_line:  # checks if string is empty, moves cnc to home z height
            ##print("Sending gcode:" + str(cleaned_line))
                # converts string to byte encoded string and append newline
            command = str.encode(gcode + '\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]

        ##print('End of gcode')
        with open("position.csv", 'w') as csvfile: #update postion of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)

def home_xy(GRBL_port_path, gcode, home, x, y, z):
    with serial.Serial(GRBL_port_path, BAUD_RATE) as ser:
        send_wake_up(ser)
        cleaned_line = remove_eol_chars(remove_comment(gcode))
        if home: #reset x and y coordinate system of device
            ##print(f'G92 X{x} Y{y}\n')
            command = str.encode(f'G92 X{x} Y{y}\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            x = 0
            y = 0
            z = 0
        if cleaned_line:  # checks if string is empty, moves cnc to home x and y position
            ##print("Sending gcode:" + str(cleaned_line))
            # converts string to byte encoded string and append newline
            command = str.encode(gcode + '\n')
            ser.write(command)  # Send g-code

            wait_for_movement_completion(ser, cleaned_line)

            grbl_out = ser.readline()  # Wait for response with carriage return
            ##print(" : ", grbl_out.strip().decode('utf-8'))
            position = [x, y, z]

        ##print('End of gcode')
        with open("position.csv", 'w') as csvfile: #updates position of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)


if __name__ == "__main__":
    n = 1
    with open('position.csv', mode='r') as file: #finds position of cnc
        csvFile = csv.reader(file)
        for lines in csvFile:
            if n == 1:
                position = lines
                #print(position)
            n = n + 1
    ##print(position)
    x = position[0]
    y = position[1]
    z = position[2]
    ##print(x)
    ##print(y)
    ##print(z)
    home = True
    gcode = f"G01 Z0 F500"
    ##print(f"G01 Z0 F500")
    home_z(GRBL_port_path, gcode, home, x, y, z) #moves cnc to home z height first
    gcode = f"G01 X0 Y0 F500"
    ##print(f"G01 X0 Y0 F500")
    home_xy(GRBL_port_path, gcode, home, x, y, z) #moves cnc to home x and y position

    #with open("gcode.gcode", "w") as f:
     #   f.truncate()
      #  print(f"G91 X{str(-curr_x)} Y{str(-curr_y)}", file=f)
       # f.close()
    #stream_gcode(GRBL_port_path, gcode_path)


    ##print('EOF')
    print("Successfully moved home")