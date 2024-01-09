# Examples on how to use the TC-RADAR toolkit
# Created: Aug 31st 2023 (7/31/23)
# Last Updated: Jan 9th 2024 (1/09/24)

import os
from TCR_Toolkit import TCR_toolkit


# change directory to the location of the data in not in the same directory
os.chdir('/rita/s1/nrb171/TCRADAR') 

# currently the code is built to handle merged data, not swath data
mergefname="tc_radar_v3k_1997_2019_xy_rel_merge_ships.nc" # merged file name

# create a TCRADAR instance of the merged data using the TCR_toolkit class
TCRADAR = TCR_toolkit(mergefname)

# here we can display all variable keys avaliable in the dataset
print(TCRADAR.data.variables.keys())

# a particular variable and relevat metadata (such as units and dimensions)
# can be accessed by calling the variable directly:
variablekey = 'total_recentered_reflectivity'
print(TCRADAR.data.variables[variablekey])

# We can also pull the values directly from the dataset scale class more concisely:
print(TCRADAR.get_vbl(variablekey))

# using the dataset instance/class we can access some metadata more easily:
# What storms are available in the dataset?
print(TCRADAR.uniquestormnames[:])

# How many missions are there in each storm?
print(TCRADAR.missionsperstorm)

# How many missions for a particular storm?
print(TCRADAR.missionsperstorm['IVAN'])


from TCR_Toolkit import missiondata

# Here we imported the missiondata class. This allows us to look at
# a particular mission's data more easily. We can create an instance
# of the mission class given a mission ID:

# We can create a mission object using the dataset instance created earlier
# and a mission ID of the mission we are interested in:
ivan14H = missiondata(TCRADAR, '040914H1')

# Now that we have our mission instance we can pull information using the
# mission instance functions, such as the date of the mission in datetime 
# format. This saves time parsing and constructing the datetime variable.
mission_datetime = ivan14H.get_datetime()
print(mission_datetime)

# We can pull the reflectivity data with the recentered data option.
# dbz = ivan14H.get_field('reflectivity', total_recentered=True)
dbz = ivan14H.get_field('reflectivity')

'''
# Scripts for plotting the data will eventually be in an optional 
# companion class built around this package and data format.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
altitude_idx = 3
altitude = TCRADAR.get_levels()[altitude_idx]

cbrange_dbz=np.arange(5,60,1)
longitude = ivan14H.get_stormvbl('merged_longitudes')
latitude = ivan14H.get_stormvbl('merged_latitudes')
PT=plt.contourf(longitude, latitude, dbz[:,:,altitude_idx], 30,levels=cbrange_dbz,
                 cmap=cm.gist_ncar, origin='lower', extend='both')
cl=plt.colorbar()
cl.ax.set_title('dBZ')

plt.xlabel('Longitude')
plt.ylabel('Latitude')

plt.title("Reflectivity at {0:.1f} km \n Mission time: {1} UTC".format(altitude, mission_datetime))

os.chdir('/rita/s0/bsr5234/TCRADAR_Package')
plt.savefig('14Hdbz%i.png' %(altitude*1000))
'''

# We can pull other variables using a the get_field() function:
div = ivan14H.get_field('divergence')

# If we want to calulate the azimuthal average of the divergence field, 
# we can feed the divergence field into the azimuthal average function.
div_aziavg = ivan14H.azimuthal_average(div)

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
PT=plt.contourf(np.arange(200), TCRADAR.get_levels()[:], div_aziavg, levels=30,
                 cmap=cm.seismic, origin='lower', extend='both')
cl=plt.colorbar()
cl.ax.set_title(r's^{-1}')

plt.xlabel('Horizontal distance from storm center (km)')
plt.title("Azimuthal Mean of Divergence\n Time: %s" %(mission_datetime))

os.chdir('/rita/s0/bsr5234/TCRADAR_Package')
plt.savefig('14H_divaziavg.png')
'''

# We can also calculate the shear-relative quadrant means. 
# This is possible without importing shear data because 
# TCRADAR includes SHIPS shear data natively.
divergence_quadrantmeans = ivan14H.quadrant_average(div)

'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
fig, axs = plt.subplots(1, 4, sharex=True, sharey=True, constrained_layout=True)
cbrange_div=np.arange(-0.001,0.0010001,.000025) #TURN ON WHEN USING CONTOUR_FILL
print(np.shape(divergence_quadrantmeans))

PTG=plt.gcf()
PTG.set_size_inches(16,4)
for qinx, ax in enumerate(axs):
    ax.contourf(np.arange(200), TCRADAR.get_levels()[:], divergence_quadrantmeans[qinx], 30, levels=cbrange_div, origin='lower', cmap=cm.seismic, extend='both') #for shaded
    ax.set_title(ivan14H.SQ[qinx][2])

os.chdir('/rita/s0/bsr5234/TCRADAR_Package')
plt.savefig('14H_div_Qavg.png')
'''

#We can calculate vertical cross sections at a specified azimuth
#creating a height-radius cross section from the center of the storm

Xsec, cross_section_distances = ivan14H.calculate_cross_section(dbz, 95)
radarheights = TCRADAR.get_levels()

'''
import matplotlib.pyplot as plt
plt.contourf(cross_section_distances, radarheights, Xsec)
os.chdir('/rita/s0/bsr5234/TCRADAR_Package') 
plt.savefig('final_crossectiontest.png')
'''