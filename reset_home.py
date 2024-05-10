"""
This is a script to reset the coordinates of the cnc tracked in a csv to (0,0,0). This should only be used after the
device has externally been moved back to its home position
"""

import sys
import csv

invalid = True
reset = input("Are you sure you would like to reset the home position? This should only be done if the device was moved"
              " externally. (Y or N): ")
while invalid: #ensures user really wants to reset the coordinate
    if reset == 'Y' or reset == 'y':
        position = [0, 0, 0]
        with open("position.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(position)
        print("Okay, resetting position of CNC to 0, 0, 0")
        invalid = False
    elif reset == 'N' or reset == 'n':
        invalid = False
        print("Okay, position will not be reset")
        sys.exit()
    else:
        print("Invalid input, please try again")
        reset = input(
            "Are you sure you would like to reset the home position? This should only be done if the device was moved"
            " externally. (Y or N): ")