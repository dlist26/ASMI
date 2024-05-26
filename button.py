"""
This is a script which can be used to check that the force sensor and cnc device are connected properly so that data can
be read and movements can be sent.


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


Parts of this code are also adapted from the python_to_GRBL github repository found here:

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
import csv

import logging
logging.basicConfig()

BAUD_RATE = 115200
GRBL_port_path = "/dev/ttyUSB0" #Change this to the desired serial port!

godirect = GoDirect(use_ble=True, use_usb=True)
device = godirect.get_device(threshold=-100)

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
            command = str.encode(f'G92 X{x} Y{y}\n')
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
            z = 0
            position = [x, y, z]

        ##print('End of gcode')
        with open("position.csv", 'w') as csvfile: #update position of cnc
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)


if __name__ == "__main__":
    # GRBL_port_path = '/dev/tty.usbserial-A906L14X'
    home = False
    right = True
    measure = True
    device.open(auto_start=False)
    while measure:
        device.start()
        sensors = device.get_enabled_sensors()
        if device.read(): #read measurements from force sensor
            for sensor in sensors:
                print(str(sensor.values))
                if abs(sensor.value) > 0.1: #if sensor is contacted, move cnc
                    if right:
                        gcode = f"G01 X50 F500"
                        right = False
                        x = 50;
                    else:
                        gcode = f"G01 X-50 F500"
                        right = True
                        x = 0
                    stream_gcode(GRBL_port_path, gcode, home, x, 0)
                    sensor.clear()
                else:
                    sensor.clear()
        time.sleep(0.1)
        device.stop()
    device.close()

