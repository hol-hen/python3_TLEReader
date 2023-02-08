#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 17:05:48 2021

@author: hollywilson
"""

'''
Required commands for the program:

    1) When prompted to "Enter TLE text file of object: ", type the file name of the text file that is the TLE (eg. ISS.txt)
    2) Enter the simple Y/N responses for the other questions

    To run program open in shell and change directory to folder with TLEreader.py and TLE file. Type python3 TLEreader.py
    
'''

#Load packages
import datetime
import calendar
from math import pi, cos, sin, sqrt, acos, floor
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#Initial data
G = 6.67408e-11   #Universal gravitational constant
mu = 3.986004418e14 #Earth gravitational constant with m^3
mu2 = 3.986004418e5 #Earth grav constant with km^3
re = 6378.1370 #Earth radius km
m_earth = 5.972e24  #Mass of Earth


#Define TLE class
class TLE:
    
    def __init__(self, i, RAAN, e, w, n, M, juliandate, year, sat_id): #constructor
        self.i = i * pi/180 #inclination converted to radians from degrees
        self.RAAN = RAAN * pi/180 #Right Ascension of the ascending node converted to radians from degrees
        self.e = e #eccentricity
        self.w = w * pi/180 #Argument of perigree converted to radians from degrees
        self.n =  n * 2*pi/86400 #mean motion converted to rad/s from rev/day
        self.M = M * pi/180 #Mean anomaly converted to radians from degrees
        self.juliandate = juliandate #julian date number
        self._year = year #last 2 digitis of year launched
        self.sat_id = sat_id #NORAD id
    
    def compute_year_launched(self): #year in TLE is just the last 2 digits
        if self._year > 60:
            self.year_launched = self._year + 1900 #add 1900 for years in 1900s (no launches before 1961)
        elif self._year <= 60:
            self.year_launched = self._year + 2000 #add 2000 for years in 2000s

    def get_year_launched(self): 
        self.compute_year_launched()
        return self.year_launched
    
    def compute_period(self):
        self.period = 2*pi/self.n #calculate orbital period in seconds
    
    def get_period(self):
        self.compute_period()
        return self.period
    
    def compute_a(self):
        self.compute_period()
        self.a = (self.period**2 * mu / (4 * pi**2) )**(1/3)/1000 #calculate semi major axis in km
    
    def get_a(self):
        self.compute_a()
        return self.a
    
    def compute_altitude(self):
        self.compute_a()
        self.altitude = self.a - re #altitude is semi-major axis- radius of earth
        
    def get_altitude(self):
        self.compute_altitude()
        return self.altitude
        
    def compute_orbit_type(self):
        self.compute_altitude()

        #label orbit type depending on altitude
        if self.altitude < 2000:
            self.orbit_type = 'LEO (Low Earth Orbit)' 
        elif self.altitude >= 2000 and self.altitude < 35700:
            self.orbit_type = 'MEO (Medium Earth Orbit)'
        elif self.altitude >= 35700 and self.altitude <= 35820:
            self.orbit_type = 'GEO (Geostationary Orbit)'
        elif self.altitude < 35820:
            self.orbit_type = 'HEO (High Earth Orbit)'
    
    def get_orbit_type(self):
        self.compute_orbit_type()
        return self.orbit_type
    
    def compute_orbits_per_day(self):
        self.compute_period()
        self.orbits_per_day = 24*60*60/self.period #calculate orbits per day
        
    def get_orbits_per_day(self):
        self.compute_orbits_per_day()
        return self.orbits_per_day
    
    def compute_speed(self):
        self.compute_a()
        self.speed = sqrt(G*m_earth/(self.a*1000))/1000 #calculate orbital speed
        
    def get_speed(self):
        self.compute_speed()
        return self.speed
    
    def compute_theta(self):
        self.E = kepler_E(self.e, self.M)
        self.theta = acos( (cos(self.E) - self.e) / (1 - self.e*cos(self.E)) ) #true anomaly in radians
        
    def get_theta(self):
        self.compute_theta()
        return self.theta
    
    def compute_r(self):
        self.compute_a()
        self.compute_theta()
        
        self.r_0 = self.a*(1-self.e**2)/(1+self.e*cos(self.theta)) #calculate position at time of TLE

        #calculate x,y,z components
        self.x = (cos(self.RAAN) * cos(self.w) - sin(self.RAAN) * sin(self.w) * cos(self.i)) * self.r_0 * cos(self.theta) + (-cos(self.RAAN) * sin(self.w) - sin(self.RAAN) * cos(self.w) * cos(self.i)) * self.r_0 * sin(self.theta)
        self.y = (sin(self.RAAN) * cos(self.w) + cos(self.RAAN) * sin(self.w) * cos(self.i)) *self. r_0 * cos(self.theta) + (-sin(self.RAAN) * sin(self.w) + cos(self.RAAN) * cos(self.w) * cos(self.i)) * self.r_0 * sin(self.theta)
        self.z = (sin(self.w) * sin(self.i)) * self.r_0 * cos(self.theta) + (cos(self.w) * sin(self.i)) * self.r_0 * sin(self.theta)
        
        self.r = np.array([self.x *1000,self.y*1000,self.z*1000]) #store as position vector
    
    def get_r(self):
        self.compute_r()
        return self.r

    def compute_v(self):
        self.compute_a()
        self.compute_theta()

        #calculate x,y,z components of velocity
        self.vx = (cos(self.RAAN) * cos(self.w) - sin(self.RAAN) * sin(self.w) * cos(self.i)) * (-(mu2/(self.a*(1-self.e **2)))**(1/2) * sin(self.theta)) + (-cos(self.RAAN) * sin(self.w) - sin(self.RAAN) * cos(self.w) * cos(self.i)) * (mu2/(self.a*(1-self.e **2)))**(1/2) * (self.e + cos(self.theta))
        self.vy = (sin(self.RAAN) * cos(self.w) + cos(self.RAAN) * sin(self.w) * cos(self.i)) * (-(mu2/(self.a*(1-self.e **2)))**(1/2) * sin(self.theta)) + (-sin(self.RAAN) * sin(self.w) + cos(self.RAAN) * cos(self.w) * cos(self.i)) * (mu2/(self.a*(1-self.e **2)))**(1/2) * (self.e + cos(self.theta))
        self.vz = (sin(self.w) * sin(self.i)) * (-(mu2/(self.a*(1-self.e **2)))**(1/2) * sin(self.theta)) + (cos(self.w) * sin(self.i)) * (mu2/(self.a*(1-self.e **2)))**(1/2) * (self.e + cos(self.theta))
        self.v = np.array([self.vx*1000,self.vy*1000,self.vz*1000]) #store as velocity vector

    
    def get_v(self):
        self.compute_v()
        return self.v


#Function to convert from julian time to gregorian time and get a string
def datestr(date_julian): 
    date = datetime.datetime.strptime(str(floor(date_julian)), '%y%j').date() #removes decimal using floor to only get date and then converts date to gregorian (YYYY-MM-DD)
    date = str(date).split('-') #splits date into year, month, day
    month = calendar.month_name[int(date[1])] #gets month as the character name
                                      
    UTCtime = float(str(date_julian)[5:]) #get the time of day as decimal point
    UTChour = '{0}'.format(str(floor(UTCtime*24)).zfill(2)) #get hour
    UTCminutes = '{0}'.format(str(floor(UTCtime*1440 % 24)).zfill(2)) #get mintue
    UTCsecond = '{0}'.format(str(round(UTCtime*86400 %60)).zfill(2)) #get second
       
    date_string = '{}:{}:{} UTC {} {} {}'.format(UTChour, UTCminutes, UTCsecond, date[2], month, date[0]) #put info into a date string with UTC as time zone
    
    return date_string

#function to read and extract information from the TLE
def readTLE(textfile):
    #open TLE file and read lines
    file = open(textfile, 'r')
    lines = file.readlines() #read lines and store as lines list
    file.close()
    
    objectname = " ".join(lines[0][:-1].split())  #remove any blank spaces and rejoin to get satellite object name
    data = lines[1][:-1].split() + lines[2][:-1].split() #store TLE as list (exclude name)
    
    #extract info from the data list
    ideg = float(data[11])   #inclination in degrees
    RAANdeg = float(data[12]) #RAAN in degrees
    e = float('0.{}'.format(data[13])) #eccentricity
    wdeg = float(data[14]) #AOP in degrees
    ndeg =  float(data[16]) #Mean motion in degrees per rev
    Mdeg = float(data[15]) # Mean anomaly in degrees
    julian_date = float(data[3]) #date of TLE observation in julian time
    date_string = datestr(julian_date) #call datestr function to get year as string
    year = float(data[2][:2])#last 2 digits of year launched
    sat_id = float(data[10]) #NORAD ID
    
    return objectname, ideg, RAANdeg, e, wdeg, ndeg, Mdeg, julian_date, date_string, year, sat_id


#Function to calculate eccentric anomaly using Newton's method
def kepler_E(e,M): 
    step_size = 1e-8 #integration step
    if M < pi:#check which quadrant Mean anomaly in (because it is an angle)
        E = M + e/2
    else:
        E = M - e/2
    
    ratio = 1

    #find convergence point = eccentric anomaly
    while abs(ratio) > step_size:
            ratio = (E - e*sin(E) - M)/(1 - e*cos(E))
            E = E - ratio
    return E


#function to create information text for object
def text(name, date_string, TLE): 
    text = [
    #General info  
    '------------------------------ {} ------------------------------'.format(name),
    'NORAD ID:                  {:0.0f}'.format(TLE.sat_id),
    'Year Launched:             {}'.format(int(TLE.get_year_launched())),
    'Orbit Type:                {}'.format(TLE.get_orbit_type()),
    
    '\n',
    
    #TLE observation info
    'TLE observation time:      {}'.format(date_string),
    'Altitude at time:          {:0.2f} km'. format(TLE.get_altitude()),
    'Speed at time:             {:0.2f} km/s'. format(TLE.get_speed()),
    
    '\n',
    
    #Orbit info
    'Orbital Period:            {:0.2f} mins'. format(TLE.get_period()/60),
    'Number of orbits per day:  {:0.2f}'. format(TLE.get_orbits_per_day()),
    'Orbit inclination:         {:0.2f} degrees'. format(TLE.i*180/pi),
    'Orbit eccentricity:        {:0.5f}'. format(TLE.e),
    'Orbit RAAN:                {:0.2f} degrees'.format(TLE.RAAN*180/pi),
    'Orbit AOP:                 {:0.2f} degrees'.format(TLE.w*180/pi)]
    return text


#function to print information on console
def text_print(name, date_string, TLE):
    a = text(name, date_string, TLE)  #call text() function to get required text
    
    print('\n') 
    for i in a:
        print(i) #print info on console


#function to save the text file        
def save_txt_file(name, date_string, TLE):
    a = text(name, date_string, TLE) #call text() function to get required text
    
    f = open('{}_orbit.txt'.format(name), 'w') #open file
    
    for i in a:
        f.write(i) #write each on a new line
        f.write('\n')
        
    f.close() #close file


#Function to plot the orbit of the satellite and save as a .png file
def orbit_plot(name, TLE):

    # define initial values in orbit
    r_0 = TLE.get_r()
    v_0 = TLE.get_v()
    t_0 = 0
    t_n = TLE.get_period()
    dt = t_n / 100000 #integration step 
    r = r_0.copy()
    v = v_0.copy()
    t = t_0
    orbit = [r.copy()]
    vel = [v.copy()]
    times = [t]
    
    #Euler Loop
    while t < t_n:
        #calculate change in r and v vectors for each integration step
        r3 = sum((-r)**2)**1.5
        r = r + dt * v
        v = v - dt * G * m_earth * (r) / r3
        t = t + dt
        #append orbit lists for plotting
        orbit.append(r.copy())
        vel.append(v.copy())
        times.append(t)
    
    #Convert into coordinates to be plotted
    x,y,z = np.asarray(orbit).T
    
    #Create a figure with 3D axes
    fig = plt.figure(1)
    ax = plt.axes(projection='3d') #create 3d plot

    #define surface of Earth
    theta1 = np.linspace(0, 2*pi, 100)
    theta2 = np.linspace(0, pi, 100)
    x_s = re * np.outer(np.cos(theta1), np.sin(theta2))
    y_s = re * np.outer(np.sin(theta1), np.sin(theta2))
    z_s = re * np.outer(np.ones(np.size(theta1)), np.cos(theta2))
       
    #plot surface of Earth
    ax.plot_wireframe(x_s, y_s, z_s, color='#23395d', alpha = 0.2) #changed transparency of wireframe Earth because matplotlib does not integrate well with the plotted lines and would not be able to see orbit at all otherwise.
    
    # Plot body 2 orbit
    ax.plot3D(x/1000,y/1000,z/1000, 'r'); #plotting in km 
    ax.scatter(r[0]/1000,r[1]/1000,r[2]/1000, c='r', marker='o'); #mark satellite point at TLE
    ax.text2D(0.05, 0.93, "Orbit of {}".format(name) , fontsize = 14, transform=ax.transAxes) #title

    #create axis limits, to ensure scale is right for all orbits
    ax.set_xlim3d(-TLE.get_a(), TLE.get_a()) #define in terms of apogee (furthest point form earth) to make all sats fit
    ax.set_ylim3d(-TLE.get_a(), TLE.get_a())
    ax.set_zlim3d(-TLE.get_a(), TLE.get_a())

    #add x, y, z to axis labels for orientation
    ax.set_xlabel("x", fontsize = 12)
    ax.set_ylabel("y", fontsize = 12)
    ax.set_zlabel("z", fontsize = 12)

    #remove axis numbers (takes away from the visualisation)
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.zaxis.set_ticklabels([])

    #save graph as .png file
    plt.savefig('{}_orbit.png'.format(name))


#Function to ask to save text file and then do so if yes
def ask_save_graph(name, TLE):
    graph = 0
    while graph not in ("Y", "y", "N", 'n'):
        graph = input('Do you want to get a visualisation of the orbit? (Y/N): ') #ask for user response

    if graph == 'Y' or graph == 'y': #if yes call function orbit_plot() to run plot sequence
        orbit_plot(name, TLE)
        print('Your graph has been saved as {}_orbit.png'.format(name)) #tell user what it has been saved as


#Function to ask to save text file and then do so if yes
def ask_save_txt(name, date_string, TLE):
    text = 0
    while text not in ("Y", "y", "N", 'n'):
        text = input('Do you want to save the data of the orbit as text file? (Y/N): ') #ask for user response

    if text == 'Y' or text == 'y':
            save_txt_file(name, date_string, TLE) #if yes call function save_txt_file() to run save sequence
            print('Your file has been saved as {}_orbit.txt'.format(name)) #tell user what it has been saved as


#Function to restart or exit program
def exit_restart():
    run_exit = 0
    
    while run_exit not in('', 's', 'S'):
        run_exit = input('- - - - - Press enter to exit or S to start again- - - - -') #ask for user response
               
    if run_exit == '':
        os._exit(os.EX_OK) #ends program
    elif run_exit == 'S' or 's':
        main()   #restarts program


#Main script
def main():
    
    #Title and instructions
    print('\n')
    print('- - - - - - - - - - TLE Reader - - - - - - - - - -')
    print("This program analyses the TLE of an object")
    print('\n')

    #Get file name of TLE as input    
    objectfile = input('Enter TLE text file of object: ')
    
    #read TLE file by calling function READtle()
    name1, ideg1, RAANdeg1, e1, wdeg1, ndeg1, Mdeg1, juliandate1, date_string1, year1, sat_id1 = readTLE(objectfile) 
    
    #instantiate class TLE by using output of readTLE()
    TLE1 = TLE(i = ideg1, RAAN = RAANdeg1, e = e1,
                   w = wdeg1, n = ndeg1, M = Mdeg1,
                   juliandate = juliandate1, year = year1, 
                   sat_id = sat_id1)
    
    
    text_print(name1, date_string1, TLE1) #print information by calling text_print() function

    print('\n')

    ask_save_graph(name1, TLE1) #ask and save graph by calling ask_save_graph() function

    print('\n')

    ask_save_txt(name1, date_string1, TLE1) #ask and save text file by ask_save_txt() function

    print('\n')

    exit_restart() #ask to exit or restart by calling exit_restart() function


if __name__ == "__main__":
   main() #run main function at start of program
