### WCL and Thurst-to-Weight ratios 
### Determing the wing loading and thrust to weigth ratio is an important part of any aircraft design. 
### A high WCL displays that an aircarft will be quick with stubber and smaller wing, while a lower WCL...
### displays  a bigger with(higher Aspect ratio) that is generally slower potentially more manueverable 
### To estimate these parameters, the estimated weight of the aircraft must be computed.
### Additoinally, the wing loading must be estimated from historical data and design parameter
### If the aircraft that is being designed is for gliding, a lower wing loading can be assumed, and if...
### a quicker dog-fight aircarft is being designed, then a higher wing loading must be assumed. 
### The initial number estimation can be found from pervious designs. D. Rayemer: Aircraft Design:Conceptual Approach...
### has a table(5.5 page 84) that shows different wing loadings from different aircraft types. 
### The thrust to weight ratio is also another parameter that can be decided based on the intended misison design
### However, this will affect the type of motor required, which might increase the weight, hence increasing the wing loading.


### This process is highly iterative process, due to the fact that weight is subject to change
### Through this process, the area of the wing be calcualted. 
### The thrust will be calculated, Cl will be calculated, and drag will be calcualted.
### There will be data representation to illustrate the changes in AoA with different parameters(thrust,Cl,etc.)
### XLR5 will be utilized to compute the airfoil parameter design. 

## Weight estimation
import numpy as np
import math
import matplotlib.pyplot as plt
from PIL import Image
from sympy import symbols, Eq, solve
import scipy 
#mass_servo = []
#mass_battery = []
#mass_sec = []
#mass_motor = []
#mass_prop = []
#mass_receiver = []
#mass_control_horn_rods = []
#mass_fuselage = []
#mass_total =weight_servo+weight_battery+weight_sec+weight_motor+weight_prop+weight_receiver+weight_control_horn_rods+weight_fuselage
mass = 1.7
a = 0.71 # table 6.3
c = 0.48 # table 6.3
g = 9.81
ib_kg_conversion = 2.2
kg_ib_conversion = 1/ib_kg_conversion
weight_pound = mass*kg_ib_conversion
weight_total = mass*g
fuselage_length = a*weight_pound**c #Table 6.3 from Conceptual Approach
thrust_weight = 0.75
thrust = thrust_weight*weight_total
print("The Thrust is %.3f, and the total weight is %.3f\n" %(thrust,weight_total))
## Solving for s_wing 
s_wing = symbols("s_wing")
WCL = 5 # Based on mission design parameters(Long distance and not fast)
expr= (weight_total/((s_wing)**(3/2))) - WCL
sol = solve(expr)
s_wing = sol[0]
wing_loading = weight_total/s_wing
watt_ratio = 125 #units in watt per pound. This is based on the YouTube Video: "How to design your motor + esc_ battery Combination"
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
print("You need a %.3f watt motor, a %.3f volt battery, and an ampere esc of %.3f\n" %(power,battery_4,current_desirable[3]))

## Determining Prop Sizing 
#Input data from a variety of motors for the particular wattage and determine the most efficient motor for your design 

    

#### For stability and control of an aircraft, a thourough simulation of the aircraft must be made. 
#### The aircraft must be stable in pitch, roll, and yaw axis. 


#### For an aircrat to be statically stable in the pitch axis, the aircraft must have a positive Static Margin.
#### Kn referes to h_n - h, where h_n is the N.P ratio relative to M.A.C and h is the C.G. ratio relative to M.A.C
#### Another way to determine longitudinal stability is by observing the CmvsCl graph. 
#### For static stability, the slope must be negative at a positive Cl for upright flight 
#### Same goes with CmvsAlpha, where alpha is the angle of attack at the zero lift line
#### When these conditions are satisified, the airplane is considered to have Positive pitch Stiffness.
#### At trim, Cm = 0
#### To improve performance of Rc Model the aerodynamics of the RC plane must be taken into consideration.
#### Mass should be minimized as much as possible. Attempts to minimze mass will be included in this code
#### For lateral stability, cl_beta must be negative, and for stability in the yawing axis, Cn_beta must be positive.
#### More about dynamic stability in the next section 


#Assumptions: Fixed wing, with no aeroelasticity. Conventinal Aircraft configuration.
#Determining Vertical tail sizing based on wing sizing. 

## WING GEOMETRY CALCULATIONS 
AR_wing = 9
b_wing = math.sqrt(AR_wing*s_wing) #Span(both wings)
taper = 0.45 #Taper Ratio of wing.Number based on Raymer's book for optimal taper for subsonic flight
c_r = (2*s_wing)/(b_wing*(1+taper)) #Root chord of wing. To be determined by CFD simulations
c_t = taper*c_r #Tip chord wing. Equation from Raymer's book. Layout and configuration chapter 7
c_bar = 2/3*(c_r)*((1+taper+taper**2)/(1+taper)) #M.A.C of the wing
y_bar = (b_wing/2)*(1/3)*((1+2*taper)/(1+taper)) #M.A.C location from the tip chord. 

## Layout Sizing(Area calculations of HT and VT)
VH_bar = 0.6 #Based on table 6.4 page 112 from Raymer's book
l_bar = 0.6*fuselage_length #Distance from wing A.C. to tail A.C... To be determine by the weight required for operation
s_tail = (VH_bar*s_wing*c_bar)/l_bar
Vv_bar = 0.04
s_vtail = (Vv_bar*s_wing*b_wing)/l_bar
AR_wetted = .8 ## Estimated from Fig 3.6 
Sref_Swet = AR_wetted/AR_wing
L_D_design = 14 #Estimated from figure 3.6 as well for subsonic aircrafts
## TAIL CALCULATIONS(Horizontal tail and Vertical Tail)
AR_tail = 2/3*(AR_wing)
b_tail = math.sqrt(AR_tail*s_tail)#Span(both HT tails)
taper_tail = 0.45 #Taper ratio of tail
c_r_tail = (2*s_tail)/(b_tail*(1+taper_tail)) #Root chord of tail
c_t_tail = taper_tail*c_r_tail #Tip chord of tail
c_bar_tail = 2/3*(c_r_tail)*((1+taper_tail+taper_tail**2)/(1+taper_tail)) #M.A.C of tail
y_bar_tail = (b_tail/2)*(1/3)*((1+2*taper_tail)/(1+taper_tail)) #M.A.C location from the tip chord.
AR_vtail = 1.5 #Table 4.3
b_vtail = math.sqrt(AR_vtail*s_vtail)
taper_vtail = 0.45 #Table 4.3
c_r_vtail = (2*s_vtail)/(b_vtail*(1+taper_vtail)) #Root chord of tail
c_t_vtail = taper_vtail*c_r_vtail #Tip chord of tail
c_bar_vtail = 2/3*(c_r_vtail)*((1+taper_vtail+taper_vtail**2)/(1+taper_vtail)) #M.A.C of the wing
y_bar_vtail = (b_vtail/2)*(1/3)*((1+2*taper_vtail)/(1+taper_vtail)) #M.A.C location from the tip chord.
## Calculation of aerodynamics forces 

#### Create a solidworks model that automatically inputs values found into this program. 
#### Create a rough sketch of the model and determine the airfoil to determine K 
#### Airfoil optimization programs are found in the book by Alexander on airfoil theory. 
#### Look for manfacturing process and how that would apply with designs. 
### The airfoil decided on is E423 Airfoil. The decision was made based upon previously made aircrafts. 
rho = 1.255 #Expected altitude is 60 meters.Thus, assuming sea level. Units are also meteric
n_p = 0.8 #the average power efficiency is 0.8. However, this must be researched even further for the particular motor used
velocity_cruise = (power/thrust)*n_p #n_p is efficiency of engine, and power refers to the motor.
print("My Velocity Cruise is %.3f meteres per second"%(velocity_cruise))
nu = 1.81e-5 #based on  Standard sea level measurements
Re = (rho*velocity_cruise*c_bar)/(nu) #Reynolds Number
speed_sound = 343 #meters/second
#Calculating a_w 
M = velocity_cruise/speed_sound #mach number 
beta = np.sqrt(1-M**2)
mid_chord_angle = [] #in radians
#cl_alpha were obtained through XFLR5 computation 
alpha_1 = 4
alpha_2 = 5
cl_1 =  1.52161 
cl_2 =  1.52161 
cl_alpha_2D =((cl_2-cl_1)/(alpha_2-alpha_1))*(180/np.pi) #in radians 
k = (beta*cl_alpha_2D)/(2*np.pi)
####a_w = (2*np.pi*AR_wing)/(2+ math.sqrt((((AR_wing**2)*(beta**2)/(k**2))*(1+ ((math.tan(mid_chord_angle))/beta)+4)))) #Cl_alpha_wing. Formula found in Bernard's book for stability and control
#Calculating Downwash per angle from Appendix B.5 Downwash from Brenard's book(Dynamics of Flight)
k_a = 1/AR_wing - 1/(1+AR_wing**1.7)
k_lambda = (10-(3*taper))/7
h_h = 0.1 #in meters. This is the hieght difference between the mac and the tail mac... Optimization shall be conducted here for a low downwash angle per degree
k_h = (1-abs(h_h/b_wing))/((((2*l_bar)/b_wing))**1/3)
quarter_chord_ang = ... #in radians
###epsilion_alpha = 4.4*(k_a*k_lambda*k_h*(math.cos(quarter_chord_ang))**(1/2))**1.19
#Calculating cl_alphas for different components of the aircraft
fitness_ratio = 3 #This is the ratio of fuselage length per max diameter. This is based on Chapter 6 of Raymer's Textbook on Fuselage
fuselage_diameter = fitness_ratio*fuselage_length
body_diameter = fuselage_diameter/b_wing
print(body_diameter)