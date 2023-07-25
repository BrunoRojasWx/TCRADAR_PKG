# Class that contains all the tools to process TDR data from TC-RADAR for Tropycal merge
# TC-RADAR is available here: https://www.aoml.noaa.gov/ftp/pub/hrd/data/radar/level3/
# Also to use FLIGHT+ dataset?
# Bruno S. Rojas    - PhD Student,          Penn State,              bsr5234@psu.edu
# Michael Fischer   - Assistant Scientist,  UM/CIMAS NOAA/AOML/HRD,  michael.fischer@noaa.gov

import os
srcpath=os.getcwd()
os.chdir(srcpath)
import netCDF4 as NC
import numpy as np
import datetime as dt
import glob


"""
Data coverage figure - Flight+, SHIPS shear (maybe get from TC-RADAR, but need to interpolate)
Types of variables:
    Reflectivity
    Tan Wind
    Rad wind
    Vort
    Vertical vel
    Convergence/Div
    Momentum (?)
Plan views
Plan views with shear dir/mag and quadrant lines
Azimuthal averages (including data coverage contours 50% and 25%)

Shear quadrant averages
Swath mean sections
CFADs
    Over whole area
    Within a radius
    Annulus
Shear Quadrant CFADs
Reflectivity CFADs
VV CFADs
CFAD quadrant anomalies
Vertical mass transport plots (need same capability as CFADS (quadrants, annuli))
"""

class TCR_toolkit:

    # def __init__(self, TCRADARpath):
    #     os.chdir(TCRADARpath)
    # def __init__(self,mergedfilename,swathsfilename):
        # self.latest_version='v3j'
        # swathsfilename="tc_radar_"+self.latest_version+"_combined_xy_rel_swath_ships.nc"
        # mergedfilename="tc_radar_"+self.latest_version+"_combined_xy_rel_merge_ships.nc"
        
        # if mergedfilename==None and swathsfilename==None:
        #     print('*** Needs at least one swath or merged dataset ***')

        # try:    
        #     self.swaths = NC.Dataset(swathsfilename)
        # except:
        #     self.swaths = None
        #     print('No swaths dataset imported.')
        # try:    
        #     self.merged = NC.Dataset(mergedfilename)
        # except:
        #     self.merged = None
        #     print('No merged dataset imported.')
        # print(f"Using version {self.latest_version}")

    def __init__(self,datafilename):

        try:    
            self.data = NC.Dataset(datafilename)
        except:
            self.data = None
            print('No dataset imported.')
        
        self.levels = self.data.variables['level']
        self.allmissionID = list(self.data.variables['mission_ID'][:])
        self.stormnames = self.data.variables['storm_name']
        self.uniquestormnames=list(dict.fromkeys(self.stormnames)) #remove duplicates
        self.missionsperstorm={i:list(self.stormnames[:]).count(i) for i in list(self.stormnames[:])}
        self.number_of_swaths = self.data.variables['number_of_swaths'][:]
        # self.levels = self.data.variables['level']
        # self.levels = self.data.variables['level']
        # self.levels = self.data.variables['level']
        # self.levels = self.data.variables['level']
        # self.levels = self.data.variables['level']
        # self.levels = self.data.variables['level']

    def get_levels(self):
        return self.levels[:]

    def get_vbl(self, variablekey):
        return self.data.variables[variablekey][:]

        # swath_filenames=glob.glob("tc_radar_*_combined_xy_rel_swath_ships.nc")
        # while len(swath_filenames)>1:
        #     TCRversionspresent = []
        #     for TCRversionfile in swath_filenames:
        #         TCRversion=TCRversionfile.split('_')[2]
        #         TCRversionspresent.append(TCRversion)
        #     print("More than one version of the swath file is present.")
        #     print(f"The {len(TCRversionspresent)} versions available are: {TCRversionspresent}")
        #     selectedversion=input("Enter selected version: ")
        #     swath_filenames=["tc_radar_"+selectedversion+"_combined_xy_rel_swath_ships.nc"]
                
        # if len(swath_filenames)==1:
        #     TCRversion=swath_filenames[0].split('_')[2]
        #     swathsfilename = swath_filenames[0]
        #     self.swaths = NC.Dataset(swathsfilename)
        #     print(f"Using version {TCRversion}")



        # print(fname)
        # swathsfilename = "tc_radar_v3j_combined_xy_rel_swath_ships.nc"
        # mergedfilename = "tc_radar_v3j_combined_xy_rel_merge_ships.nc"

        # self.swaths = NC.Dataset(swathsfilename)
        # self.merged = NC.Dataset(mergedfilename)

        # self.filename=filename
        # self.data = NC.Dataset(filename)
        # self.allvblkeys = self.data.variables.keys()


class missiondata:

    def __init__(self,radardataset, stormID):
        self.radardata = radardataset
        self.stormindex = radardataset.allmissionID.index(stormID)
    
    def get_stormvbl(self, variablekey):
        return self.radardata.get_vbl(variablekey)[self.stormindex]
    
    def get_datetime(self):
        import datetime as dt
        self.year = self.radardata.data.variables['merge_year'][self.stormindex]
        self.month = self.radardata.data.variables['merge_month'][self.stormindex]
        self.day = self.radardata.data.variables['merge_day'][self.stormindex]
        self.hour = self.radardata.data.variables['merge_hour'][self.stormindex]
        self.minute = self.radardata.data.variables['merge_min'][self.stormindex]

        self.datetime= dt.datetime(year=self.year, month=self.month, day=self.day,
                                    hour=self.hour, minute=self.minute)
        return self.datetime
    
    def get_field(self, fieldtype, recentered=False, total_recentered=False):
        '''
        Returns 3-dimensional grid of storm radar data (latitude/meridional 
        displacement, longitude/zonal displacement, vertical level)

        Field options:

        zonal_wind, meridional_wind, vertical_velocity, reflectivity,
        wind_speed, radial_wind, tangential_wind, earth_relative_zonal_wind,
        earth_relative_meridional_wind, relative_vorticity, divergence
        '''
        if recentered==False and total_recentered==False:
            field = self.get_stormvbl('merged_{}'.format(fieldtype))
        elif recentered==True and total_recentered==False:
            field = self.get_stormvbl('recentered_{}'.format(fieldtype))
        elif total_recentered==True:
            field = self.get_stormvbl('total_recentered_{}'.format(fieldtype))            
            
        return field