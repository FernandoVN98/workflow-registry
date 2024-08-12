#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Unified CyberShake Workflow (UnifiedCSWFlow)
#
# Author:  Juan Esteban Rodríguez, Josep de la Puente
# Contact: juan.rodriguez@bsc.es, josep.delapuente@bsc.es
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

################################################################################
# Module imports
import os
import sqlite3
import numpy as np
from tqdm import tqdm
import pandas as pd
from turfpy import measurement
from geojson import Point, Feature
from turfpy.measurement import destination

from modules.cybershake_workflows.utils.DAL import SQLiteHandler

################################################################################
# Methods and classes

# Rupture generator
class Ruptures_FM():
    # Initialization method
    def __init__(self, opath, database):
        # Set the output path
        self.opath = opath
        
        # Initialze focal mechanism
        self.fm = {"dip": 0.0, "strike": 0.0, "rake": 0.0}
        
        # Create DB  sconnection        
        self.dal = SQLiteHandler(database, True)
        
    # Set the focal mechanism
    def setFocalMechanism(self, dip, strike, rake):
        self.fm = {"dip": dip, "strike": strike, "rake": rake}
        
    # Generate ruptures
    def generateRuptures(self, id, site, fms):

        # Obtain all ruptures in the DB
        rows = self.dal.getRupturesDict()
        count = len(list(rows))
                
        # For each Rupture found
        pbar = tqdm(rows, desc='[' + site +'] Generating rupture files', total=count, position=id, leave=False)
        df = pd.read_csv(fms)
        cont = 0
        for row in pbar:
            print('cont',cont)
            # Creating output directory
            path = "%s/%s/%s" % (self.opath, row['Source_ID'], row['Rupture_ID']);
            self.fm['dip'] = df['dip'].iloc[cont]
            self.fm['rake'] = df['rake'].iloc[cont]
            self.fm['strike'] = df['strike'].iloc[cont]

            if os.path.isdir(path):
               continue
            os.makedirs(path, exist_ok=True)

            # Building the filename
            file = "%s_%s.txt"% (row['Source_ID'], row['Rupture_ID'])

            #pbar.set_description("Generating rupture files --> %s/%s" % (path, file))
            #print("Generating file: " + path + "/" + file)

            # Generate the rupture
            with open(path + "/" + file, 'w') as rup:

                slat = float(row['Start_Lat'])
                elat = float(row['End_Lat'])
                slon = float(row['Start_Lon'])
                elon = float(row['End_Lon'])
                sdepth = float(row['Start_Depth'])
                edepth = float(row['End_Depth'])
                gridspa = float(row['Grid_Spacing'])


		# Computing the lower corners to compute the width of the fault depending on the width
                origin = Feature(geometry=Point([slon,slat]))
                r = (edepth-sdepth)/np.sin(self.fm['dip']*np.pi/180)
                y_end = (edepth-sdepth)
                distance = r*np.cos(self.fm['dip']*np.pi/180)
                if (self.fm['strike'] < 0) or (180 < self.fm['strike'] < 360):
                    bearing = self.fm['strike'] + 90
                elif (self.fm['strike'] > 0) and (0 < self.fm['strike'] < 180):
                    bearing = self.fm['strike'] - 90    
                elif self.fm['strike'] == 0:
                    bearing = self.fm['strike'] + 90
                elif self.fm['strike'] == 180:
                    bearing = self.fm['strike'] - 90                


                options = {'units': 'km'}
                desti_tmp = destination(origin,distance,bearing,options)
                Lon_C1_origin, Lat_C1_origin = desti_tmp.geometry.coordinates
	
                # Computing the Number of cols and rows according to the dip-dependent faul-width

                start = Feature(geometry=Point([slon,slat]))
                end = Feature(geometry=Point((Lon_C1_origin, Lat_C1_origin)))
                distance_start_c1 = measurement.distance(start,end,"km")
                
                dist_Width = r    # fault width 
                nrows = int(dist_Width/gridspa)+1  # along dip


                # Determine the length of the Fault 
                start = Feature(geometry=Point((slon, slat)))
                end = Feature(geometry=Point((elon,elat)))
                dist_Length = measurement.distance(start,end,"km")
                ncols = int(dist_Length/gridspa) + 1  # along strike

                # Dumping header information
                rup.write("Probability = " + str(row['Prob']) + "\n")
                rup.write("Magnitude = " + str(row['Mag']) + "\n")
                rup.write("GridSpacing = " + str(row['Grid_Spacing']) + "\n")
               # rup.write("NumRows = " + str(row['Num_Rows']) + "\n")
               # rup.write("NumCols = " + str(row['Num_Columns']) + "\n")
                rup.write("NumRows = " + str(nrows) + "\n")
                rup.write("NumCols = " + str(ncols) + "\n")
                rup.write("#   Lat         Lon         Depth      Rake    Dip     Strike\n")

                #ncols = int(row['Num_Columns'])-1
                #nrows = int(row['Num_Rows'])-1

                #sphereDist = haversine((x2, y1), (x1, y2)) * 1000
                #print(str(round((sphereDist/200) + 0.5)) + " " + str(row['Num_Columns']))

                print(slat,elat,slon,elon,sdepth,edepth,ncols, nrows, self.fm['strike'], self.fm['dip'], self.fm['rake'], gridspa,cont)

                ## Generate the 3 axis points only for faults pure strike-slip
                #x = self._generateAxis(slat, elat, 0, ncols, 1e-10)
                #y = self._generateAxis(slon, elon, 0, ncols, 1e-10)
                #z = self._generateAxis(sdepth, edepth, 0, nrows, 1e-2)

		# Generate the 3 axis considering the dip and strike of the fault
                vecTestDepth = self._generateAxisFMs(slat, elat, slon, elon, sdepth, edepth, ncols, nrows, self.fm['strike'], self.fm['dip'], self.fm['rake'], gridspa)
                x = vecTestDepth[:,0] # Latitude
                y = vecTestDepth[:,1] # Longitude
                z = vecTestDepth[:,2] # Depth


                # Write to file the ruptures
                s =  "    "
                
                for k in np.arange(len(vecTestDepth)):
                    rup.write(str(x[k]) + s + str(y[k]) + s + str(z[k]) + s + str(self.fm['rake'])\
                     + s + str(self.fm['dip']) + s + str(self.fm['strike']) + "\n")

                # This for loop was before in the original code
                #for j in z:
                #    print('j',j)
                #    for i in np.arange(0, ncols):
                #        rup.write(str(x[i]) + s + str(y[i]) + s + str(j) + s + str(self.fm['rake'])\
                #         + s + str(self.fm['dip']) + s + str(self.fm['strike']) + "\n")
            cont += 1
            
    # Generate axis
    def _generateAxis(self, ini, end, offset, num, prec):
        # Calculate the step
        step = abs(end-ini) / num

        prec = 1/prec
        step = int(step *prec)/prec

        # Generate the axis
        if ini == end:
            axis = np.repeat(ini, num)
        elif ini < end:
            axis = np.arange(ini+offset, end, step)
        else:
            axis = np.arange(end+offset, ini, step)

        return axis        
            
    
    # Generate axis including the dip
    def _generateAxisFMs(self, slat, elat, slon, elon, sdepth, edepth, ncols, nrows, stk, dip, rake, gridspa):
        VecTestDepth = []
        for j in np.arange(nrows):  # loop along dip (depth)
            for i in np.arange(ncols):   # loop along strike
                origin = Feature(geometry=Point([slon,slat]))
                distance = j*gridspa*np.cos(dip*np.pi/180)
                if (stk < 0) or (180 < stk < 360):
                    bearing = stk + 90
                elif (stk > 0) and (0 < stk < 180):
                    bearing = stk - 90    
                elif stk == 0:
                    bearing = stk + 90
                elif stk == 180:
                    bearing = stk - 90

                    # bearing = 0
                options = {'units': 'km'}
                desti_tmp = destination(origin,distance,bearing,options)
                Lon_tmp_origin, Lat_tmp_origin = desti_tmp.geometry.coordinates

                origin = Feature(geometry=Point([Lon_tmp_origin,Lat_tmp_origin]))
                distance = i*gridspa
                bearing = stk
                options = {'units': 'km'}
                desti_tmp = destination(origin,distance,bearing,options)
                Lon_tmp_C1, Lat_tmp_C1 = desti_tmp.geometry.coordinates
                VecTestDepth = np.append(VecTestDepth,[Lat_tmp_C1, Lon_tmp_C1, sdepth+(j*gridspa*np.sin(dip*np.pi/180))])                

        VecTestDepth = VecTestDepth.reshape(int(len(VecTestDepth)/3),3)
 
        return VecTestDepth
      
        
    
        
        
    
          
    
