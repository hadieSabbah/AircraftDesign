import math
import numpy as np
AR_wing = 5
power = 5
thrust = 15
c_bar = 1.2
rho = 1.255 #Expected altitude is 60 meters.Thus, assuming sea level. Units are also meteric
n_p = 0.8 #the average power efficiency is 0.8. However, this must be researched even further for the particular motor used
velocity_cruise = (power/thrust)*n_p #n_p is efficiency of engine, and power refers to the motor.
print("My Velocity Cruise is %.3f meteres per second"%(velocity_cruise))
nu = 1.81e-5 #based on  Standard sea level measurements
Re = (rho*velocity_cruise*c_bar)/(nu) #Reynolds Number
speed_sound = 343 #meters/second
#Calculating a_w 
M = velocity_cruise/speed_sound #mach number 
print("The mach is %.3f\n"%(M))
beta = np.sqrt(1-M**2)
print("Beta is %.3f\n"%(beta))
mid_chord_angle_list = []
for mid_chord_angle in range(0,20,5):
    mid_chord_angle_list.append(mid_chord_angle)
print(mid_chord_angle_list)
mid_chord_angle_list = [number*(math.pi/180) for number in mid_chord_angle_list ] # in radians
#cl_alpha were obtained through XFLR5 computation 
alpha_1 = 4
alpha_2 = 8
cl_1 =  1.52161 
cl_2 =  1.87636 
cl_alpha_2D =((cl_2-cl_1)/(alpha_2-alpha_1))*(180/np.pi) #in radians 
k = (beta*cl_alpha_2D)/(2*np.pi)
a_w = []
mid_chord_angle = 20
a_w= ((2*np.pi*AR_wing)/(2 + math.sqrt(((AR_wing**2 * beta**2)/(k**2)) *(1+ (math.tan(mid_chord_angle)**2)/(beta**2))+4))) #Cl_alpha_wing. Formula found in Bernard's book for stability and control
print(mid_chord_angle)
print(a_w)
print(k)
print(cl_alpha_2D)
print(math.sqrt(((AR_wing**2 * beta**2)/(k**2)) *(1+ (math.tan(mid_chord_angle)**2)/(beta**2))+4))
