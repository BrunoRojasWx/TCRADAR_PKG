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
from matplotlib.dates import date2num


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
    '''
    Creates an instance for a particular mission

    Takes the radar dataset instance(swath/merged), and missionID
    Attributes and functions will provide data for the given mission
    '''
    def __init__(self,radardataset, missionID):
        self.radardata = radardataset
        self.stormindex = radardataset.allmissionID.index(missionID)
    
    def get_stormvbl(self, variablekey):
        '''
        Returns any variable present in the database for the storm instance
        '''
        return self.radardata.get_vbl(variablekey)[self.stormindex]
    
    def get_datetime(self):
        '''
        Returns mission date and time in datetime format
        '''
        import datetime as dt
        self.year = self.radardata.data.variables['merge_year'][self.stormindex]
        self.month = self.radardata.data.variables['merge_month'][self.stormindex]
        self.day = self.radardata.data.variables['merge_day'][self.stormindex]
        self.hour = self.radardata.data.variables['merge_hour'][self.stormindex]
        self.minute = self.radardata.data.variables['merge_min'][self.stormindex]

        self.datetime= dt.datetime(year=self.year, month=self.month, day=self.day,
                                    hour=self.hour, minute=self.minute)
        return self.datetime

    def get_SHIPS_times(self):
        '''
        Returns a list of datetimes for the times of the SHIPS data for the mission.

        The datetimes are +/- 48 hours from the nearest 6-hly time to the mission time

        Index 8 corresponds to the nearest 6-hly time to the mission time.
        '''
        missiontime = self.get_datetime()   #pull the mission time

        #calculate the nearest ships time (nearest synoptic 6-hly time)
        base_time = dt.datetime.combine(missiontime.date(), dt.datetime.min.time())
        sixh_interval = dt.timedelta(hours=6)

        # Calculate the time difference from the base time
        time_difference = missiontime - base_time
        
        # Round to the nearest 6-hour interval
        rounded_hours = round(time_difference.total_seconds() / sixh_interval.total_seconds())
        
        # Calculate the nearest 6-hour time
        nearest_time = base_time + (rounded_hours * sixh_interval)

        times_forward = [nearest_time + dt.timedelta(hours=i*6) for i in range(1, 9)]
        times_backward = [nearest_time - dt.timedelta(hours=i*6) for i in range(1, 9)]
        times_backward.reverse()  # Reverse the backward list to maintain chronological order
        
        return times_backward + [nearest_time] + times_forward

    def get_shear(self):
        '''
        Returns: shear direction (degrees meteorological, vector heading 90 to the east, westerly), shear magnitude (kts)
        Interpolated to the mission time from SHIPS.

        Uses the SHRD and SHTD SHIPS parameters:
        SHRD: 850-200 hPa shear magnitude (kt * 10) vs time (200-800 km)
        SHTD: Heading (deg) of above shear vector. Westerly shear has a value of 90 deg

        SHIPS variables documented here:
        https://rammb.cira.colostate.edu/research/tropical_cyclones/ships/data/ships_predictor_file_2022.pdf
        '''
        #Pull the SHIPS shear variables
        shrDir=self.get_stormvbl('shtd_ships')[:]  #shear heading degrees
        shrMag=self.get_stormvbl('shrd_ships')[:]  #shear magnitude knots
        print(shrMag)

        #convert all the datetimes to number format for interpolation
        numbertime_missiondatetime=date2num(self.get_datetime())    #get the mission time
        ships_times_numberformat=date2num(self.get_SHIPS_times())   #get the SHIPS times

        #interpolate to the mission time
        shear_direction=np.interp(numbertime_missiondatetime, ships_times_numberformat, shrDir)
        shear_magnitude=np.interp(numbertime_missiondatetime, ships_times_numberformat, shrMag)

        #calculate the U and V components of the shear
        cart_angle=-(shear_direction-90)%360 #cartesian angle
        ush=shear_magnitude*np.cos(np.deg2rad(cart_angle))
        vsh=shear_magnitude*np.sin(np.deg2rad(cart_angle))

        return shear_direction, shear_magnitude, ush, vsh

    def get_field(self, fieldtype, recentered=False, total_recentered=False):
        '''
        Returns 3-dimensional grid of storm radar data (latitude/meridional 
        displacement, longitude/zonal displacement, vertical level) (200, 200, 37)

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
    
    def get_radiusgrid(self, center_altitude=2000):
        longitude = self.get_stormvbl('merged_longitudes')
        latitude = self.get_stormvbl('merged_latitudes')
        tc_ctr_longitude = self.get_stormvbl('tc_center_longitudes')
        tc_ctr_latitude = self.get_stormvbl('tc_center_latitudes')

        center_altitude_index = list(self.radardata.get_levels()[:]).index(int(center_altitude/1000))   #pull the corresponding z-index given the altitude
        if np.isnan(tc_ctr_latitude[center_altitude_index])==True:  #check if the center has a real value
            raise Exception('Azimuthal average ERROR:\nCenter data missing at requested altitude: {} m\nSelect a new altitude'.format(center_altitude))        
        
        from tdr_tc_centering_with_example import distance
        #calculate the radius of each gridbox
        self.radiusgrid = distance(tc_ctr_latitude[center_altitude_index],tc_ctr_longitude[center_altitude_index],latitude,longitude)
        return self.radiusgrid
    
    def azimuthal_average(self, field_variable, center_altitude=2000, reflectivity_flag=False, maximum_radius=200, radial_bin_size=1):
        '''
        Returns a 2-Dimensional array of the azimuthal storm average for 
        the requested field. (height, radius) (37, 200)
        
        Field options:

        zonal_wind, meridional_wind, vertical_velocity, reflectivity,
        wind_speed, radial_wind, tangential_wind, earth_relative_zonal_wind,
        earth_relative_meridional_wind, relative_vorticity, divergence

        Set reflectivity_flag = True (False by default) when averaging
        reflectivity values. This is important to average the reflectivity in linear
        units of Z and then reconvert back to dbz. Setting to True will handle
        this conversion to Z and back to dbz.

        maximum_radius sets the maximum radius in kilometers over which averaging will 
        be performed and affects the shape of the returned array. 200 km by default.

        radial_bin_size sets the radial bin size for the azimuthal mean in km
        default is radial_bin_size = 1 km

        The bins can be retrieved for plotting by calling the self.range_bins attribute:
        e.g.: mission_class_name.range_bins
        note: this will be written over if this function is called again with a different bin size

        Data coverage can be retrieved as an attribute:
        e.g.: mission_class_name.datacoverage
        Expressed as an azimuthal percentage at each radius.
        '''
        datavariable = field_variable

        radiusgrid = self.get_radiusgrid()
        radiusgrid_rounded = np.round(radiusgrid / radial_bin_size) * radial_bin_size #rounds radii to integers for binning
        raveled_radii=radiusgrid_rounded.ravel() # ravel(flatten) the array of radii
        self.range_bins = np.arange(0,maximum_radius+radial_bin_size-1,radial_bin_size)
        gridboxes_byradius = np.histogram(raveled_radii, bins=self.range_bins)    #gives the amount of grid boxes at each radius
        gridlen = len(self.range_bins)
        aziavg=[]
        self.datacoverage=[]
        for lvindex, lv in enumerate(datavariable[0,0,:]):  #iterate over each vertical level
            raveled_data=datavariable[:,:,lvindex].ravel() #flattens the array in to one long 1D array
            keep = ~np.isnan(raveled_data) #using this as an index allows to only look at existing data values
            if reflectivity_flag==True:
                raveled_data[keep] = 10**(raveled_data[keep]/10)  #convert from dbz to Z
                # conversion needs to be done here,
                # otherwise if done earlier, all masked values will be assigned a value 
                # that will end up at 10dbz following the averaging and reconversion to dbz
            level_aziavg = np.histogram(raveled_radii[keep], weights=raveled_data[keep], bins=self.range_bins) #creates a sum of all the existing data values at each radius
            radial_bintotals = np.histogram(raveled_radii[keep], bins=self.range_bins) #gives the amount of grid boxes with data at each radius
            azimean = level_aziavg[0] / radial_bintotals[0] #takes the summed data values and divides them by the amount of grid boxes that have data at each radius
            coverage_byradius = radial_bintotals[0] / gridboxes_byradius[0]
            aziavg.append(azimean[0:gridlen])  #add the level to the azimuthal mean stack, and trim down to the max radius
            self.datacoverage.append(coverage_byradius[0:gridlen+1])
        if reflectivity_flag==True:
            aziavg = 10 * np.log10(aziavg)  # convert Z back to dbz
        self.range_bins=self.range_bins[0:-1]   #since these are bin edges, we trim the last element to have the same shape as the data to use as an x-axis
        return aziavg

    def quadrant_average(self, field_variable, center_altitude=2000, shD='none', halftype='none', reflectivity_flag=False):
        '''
        Returns 4 two-dimensional arrays of shear-relative quadrant means 
        of the given field. Quadrants are returned in the order: UR,DR,DL,UL.
        (Quadrant, height, radius) (4, 37, 200)

        Shear direction is chosed based on the SHIPS shrd variable. This is 
        an annulus mean shear between 200-800km radius. The synoptic 6-hly SHIPS
        data that is closest to the mission time is chosen. (e.g. a mission at
        2000 UTC will use the 1800 UTC SHIPS time for the shear direction)

        shD keyword argument takes a shear direction as a heading in meteorological azimuth.
        If not specified, will default to the nearest SHIPS shear direction.

        halftype keyword argument allows to average halves of the storm, either:
        'halfRL' for Right and Left of shear
        'halfDU' for Upshear and Downshear halves
        This will return a shape of (half, height, radius) (2, 37, 200)

        Field options:

        zonal_wind, meridional_wind, vertical_velocity, reflectivity,
        wind_speed, radial_wind, tangential_wind, earth_relative_zonal_wind,
        earth_relative_meridional_wind, relative_vorticity, divergence
        '''

        if shD=='none': #if no shear is explicitly passed, use the nearest SHIPS shear direction
            # Pull the shear data at the closest SHIPS time to the mission
            shD=self.get_stormvbl('shrd_ships')[8]  #index 8 is +/- 0 hours from nearest time to mission
        ups=(shD+180) % 360 #Upshear direction
        dws=shD #Downshear Direction
        rts=(shD+90) % 360 #right-of-shear 
        lfs=(shD-90) % 360 #left-of-shear

        if halftype=='halfRL':#use for R/L of shear comparison
            RoShr=[dws,ups,"Right of Shear"]#DownShear Right Quadrant
            LoShr=[ups,dws,"Left of Shear"]#DownShear Right Quadrant
            self.SQ=[RoShr,LoShr]   #right and left of shear halves

        elif halftype=='halfDU':   #use for upshear/downshear comparison
            DwShr=[ups,dws,"Downshear"]
            UpShr=[dws,ups,"Upshear"]
            self.SQ=[DwShr,UpShr]   #downshear and upshear halves

        else:   #for normal shear quadrant average
            UR=[rts,ups,"UpShear Right"]#UpShear Right Quadrant
            DR=[dws,rts,"DownShear Right"]#DownShear Right Quadrant
            DL=[lfs,dws,"DownShear Left"]#DownShear Left Quadrant
            UL=[ups,lfs,"UpShear Left"]#UpShear Left Quadrant
            self.SQ=[UR,DR,DL,UL]#All shear quadrants

        from tdr_tc_centering_with_example import get_bearing
        #pull the relevant lat/lon grid and storm center for azimuth calculation
        longitude = self.get_stormvbl('merged_longitudes')
        latitude = self.get_stormvbl('merged_latitudes')
        tc_ctr_longitude = self.get_stormvbl('tc_center_longitudes')
        tc_ctr_latitude = self.get_stormvbl('tc_center_latitudes')

        center_altitude_index = list(self.radardata.get_levels()[:]).index(int(center_altitude/1000))   #pull the corresponding z-index given the altitude
        if np.isnan(tc_ctr_latitude[center_altitude_index])==True:  #check if the center has a real value
            raise Exception('Azimuthal average ERROR:\nCenter data missing at requested altitude: {} m\nSelect a new altitude'.format(center_altitude))        
        
        # calculate the azimuth of each grid box
        azimuthgrid = np.rad2deg(get_bearing(tc_ctr_latitude[center_altitude_index],tc_ctr_longitude[center_altitude_index],latitude,longitude))
        azimuthgrid = ((-azimuthgrid+90)%360).astype(np.int32) #converts to meteo azimuth, rounds radii to integers for binning

        q_means = []
        import copy
        for shrquad in self.SQ:  # loop over each quadrant
            data = copy.copy(field_variable)    #duplicate the input field so it does not get overwritten
            Alim=shrquad[0]#edges of the desired quadrant
            Blim=shrquad[1]
            #check and ignore data outside the desired quadrant
            if Alim < Blim: #if the quadrant does not include 0 degrees azimuth (North)
                for i in range(len(latitude)):
                    for j in range(len(latitude)):
                        if (azimuthgrid[i][j] < Alim or azimuthgrid[i][j] > Blim):
                            data[i][j] = np.nan
            else:           #if the quadrant does include 0 degrees azimuth (North)
                for i in range(len(latitude)):
                    for j in range(len(latitude)):
                        if (azimuthgrid[i][j] > Alim and azimuthgrid[i][j] > Blim):
                            data[i][j] = np.nan
            q_means.append(self.azimuthal_average(data,reflectivity_flag=reflectivity_flag))    #add each quadrant to a list
        return q_means  #list of four two-dimensional arrays (quadrant, height, radius)

