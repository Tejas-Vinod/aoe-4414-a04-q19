# ecef_to_eci.py
#
# Usage: python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km
#  Converts ECEF to ECI vector components
#
# Parameters:
#  year: year
#  month: month
#  day: day
#  hour: hour
#  minute: minute
#  second: second
#  ecef_x_km: ECEF x-component in km
#  ecef_y_km: ECEF y-component in km
#  ecef_z_km: ECEF z-component in km
# Output:
#  Prints ECI components in km
#
# Written by Tejas Vinod
# Other contributors: 
#
# This work is licensed under CC BY-SA 4.0
#
# Test 
# python3 ecef_to_eci.py 2054 4 29 11 29 3.3 6778.137000 -0.028633 3838.027968 
# 5870.038832 3389.068500 3838.027968

# import Python modules
import math # math module
import sys  # argv

# "constants"
w = 7.292115 * 10**(-5) # rad/s

# helper functions

## matrix multiplication
def matrix_multiply(matrix1, matrix2):
    # Get the dimensions of the matrices
    rows_matrix1 = len(matrix1)
    cols_matrix1 = len(matrix1[0])
    rows_matrix2 = len(matrix2)
    cols_matrix2 = len(matrix2[0])
    
    # Check if matrices are able to be multiplied
    if cols_matrix1 != rows_matrix2:
        raise ValueError("Matrix multiplication is not possible. Columns of matrix1 must equal rows of matrix2.") 
    
    # Perform the multiplication
    result = [[0 for _ in range(cols_matrix2)] for _ in range(rows_matrix1)]
    for i in range(rows_matrix1):
        for j in range(cols_matrix2):
            for k in range(cols_matrix1):
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    
    return result

# initialize script arguments
year = float('nan') 
month = float('nan') 
day = float('nan') 
hour = float('nan') 
minute = float('nan') 
second = float('nan') 
ecef_x_km = float('nan') # ECEF x-component in km
ecef_y_km = float('nan') # ECEF y-component in km
ecef_z_km = float('nan') # ECEF z-component in km

# parse script arguments
if len(sys.argv)==10:
  year = float(sys.argv[1])
  month = float(sys.argv[2])
  day = float(sys.argv[3])
  hour = float(sys.argv[4]) 
  minute = float(sys.argv[5])
  second = float(sys.argv[6])
  ecef_x_km = float(sys.argv[7])
  ecef_y_km = float(sys.argv[8])
  ecef_z_km = float(sys.argv[9])
else:
  print(\
   'Usage: '\
   'python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km'\
  )
  exit()

# write script below this line
# Find GMST angle
JD = day - 32075 + int((1461 * (year + 4800 + int((month - 14) / 12))) / 4) + int((367 * (month - 2 - 12 * int((month - 14) / 12))) / 12) - int((3 * int(int(year + 4900 + (month - 14) / 12) / 100)) / 4)
JD_midnight = JD - 0.5
D_fractional = (second + 60 * (minute + 60 * hour)) / 86400
JD_fractional = JD_midnight + D_fractional
T_UT1 = (JD_fractional-2451545.0)/36525.0
GMST_sec = (67310.54841 + (876600*60*60 + 8640184.812866) * T_UT1 + 0.093104 * T_UT1**2 - 6.2*10**(-6) * T_UT1**3+86400)%(86400)
GMST_rad = (GMST_sec/240)*math.pi/180
GMST_rad = GMST_sec*w

# Rotation Matrix
r_ecef = [[ecef_x_km], [ecef_y_km], [ecef_z_km]]

Rzinv = [ 
      [ math.cos(-GMST_rad), math.sin(-GMST_rad), 0],
      [-math.sin(-GMST_rad), math.cos(-GMST_rad), 0],
      [                  0,                    0, 1] ] # rzinv(theta)^-1

r_eci = matrix_multiply(Rzinv, r_ecef)
                 
print(r_eci[0][0])
print(r_eci[1][0])
print(r_eci[2][0])   