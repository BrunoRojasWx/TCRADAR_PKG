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
    def __init__(self):
        self.default_version='v3j'
        swathsfilename="tc_radar_"+self.default_version+"_combined_xy_rel_swath_ships.nc"
        mergedfilename="tc_radar_"+self.default_version+"_combined_xy_rel_merge_ships.nc"
        self.swaths = NC.Dataset(swathsfilename)
        self.merged = NC.Dataset(mergedfilename)
        print(f"Using version {self.default_version}")



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


