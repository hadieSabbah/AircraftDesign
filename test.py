import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from sympy import symbols, Eq, solve
mass = 1.7
a = 0.71 # table 6.3
c = 0.48 # table 6.3
g = 9.81
ib_kg_conversion = 2.2
kg_ib_conversion = 1/ib_kg_conversion
weight_pound = mass*kg_ib_conversion
weight_total = mass*g
fuslage_length = a*weight_pound**c #Table 6.3
thrust_weight = 0.75
thrust = thrust_weight*weight_total
print("The Thrust is %.3f, and the total weight is %.3f\n" %(thrust,weight_total))
print("The fuslage length is %.3f meters\n The total weight of the aircarft is %.3f\n"%(fuslage_length,weight_total))
## Solving for s_wing 
s_wing = symbols("s_wing")
WCL = 5 # Based on mission design parameters(Long distance and not fast)
expr= (weight_total/((s_wing)**(3/2))) - WCL
sol = solve(expr)
s_wing = sol[0]
print('The wing area is %.3f meters sqared'%(s_wing))
wing_loading = weight_total/s_wing
## Power System preliminary selection
watt_ratio = 125 #units in watt per pound. This is based on the YouTube Video: "How to design your motor + esc_ battery Combination"
motor_watt = watt_ratio*mass*ib_kg_conversion
watt_ratio = 125 #units in watt per pound. This is based on the YouTube Video: "How to design your motor + esc_ battery Combination"
ib_kg_conversion = 2.2
kg_ib_conversion = 1/ib_kg_conversion
motor_watt = watt_ratio*mass*ib_kg_conversion
## Power System preliminary selection
ib_kg_conversion = 2.2
kg_ib_conversion = 1/ib_kg_conversion
power = watt_ratio*mass*ib_kg_conversion
print("The power needed for you power motor is %.3f\n" %(power))
#Dependent on the number of cells. Here the cell is at minimum of 3.7v per cell. The units here are all in volts. The variation can be graphed for a better optimization rate
battery_2 = 7.4
battery_3 = 11.1
battery_4 = 14.8 
battery_6 = 22.2 
battery_list = [battery_2,battery_3,battery_4,battery_6]
current_max = []
for battery in battery_list:
    current_max.append(power/battery)
current_desirable = np.divide(current_max,0.8) #Generally, you do not want the current to be going at max during flight. Heating might cause problems and will wear out the esc
print(current_desirable)
plt.plot(battery_list,current_desirable)
plt.grid()
plt.title("Battery Voltage with Current")
plt.xlabel("Voltage")
plt.ylabel("Current")
plt.show()
#### Create a solidworks model that automatically inputs values found into this program. 
#### Create a rough sketch of the model and determine the airfoil to determine K 
#### Airfoil optimization programs are found in the book by Alexander on airfoil theory. 
#### Look for manfacturing process and how that would apply with designs. 
### The airfoil decided on is E423 Airfoil. The decision was made based upon previously made aircrafts. 
AR_wing = 9
b_wing = np.sqrt((AR_wing*s_wing)) #Span(both wings)
taper = 0.45 #Taper Ratio of wing.Number based on Raymer's book for optimal taper for subsonic flight
c_r = (2*s_wing)/(b_wing*(1+taper)) #Root chord of wing. To be determined by CFD simulations
c_t = taper*c_r #Tip chord wing. Equation from Raymer's book. Layout and configuration chapter 7
c_bar = 2/3*(c_r)*((1+taper+taper**2)/(1+taper)) #M.A.C of the wing
y_bar = (b_wing/2)*(1/3)*((1+2*taper)/(1+taper)) #M.A.C location from the tip chord. 
rho = 1.255 #Expected altitude is 60 meters.Thus, assuming sea level. Units are also meteric
n_p = 0.8 #the average power efficiency is 0.8. However, this must be researched even further for the particular motor used
velocity_cruise = (power/thrust)*n_p #n_p is efficiency of engine, and power refers to the motor.
nu = 1.81e-5 #based on  Standard sea level measurements
Re = (rho*velocity_cruise*c_bar)/(nu) #Reynolds Number
print("The Reynold's Number is %.3f\n"%(Re))