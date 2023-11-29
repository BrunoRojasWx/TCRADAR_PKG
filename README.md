# TC-RADAR Toolkit

The TC-RADAR Toolkit is designed to provide functions for meteorological analysis of the NOAA P3 airborne doppler radar dataset, Tropical Cyclone Radar Archive of Doppler Analyses with Recentering [(TC-RADAR)](https://www.aoml.noaa.gov/ftp/pub/hrd/data/radar/level3/TC_RADAR_README.pdf). The TC-RADAR Toolkit is specifically developed for TC-RADAR and the dataset can be accessed [here](https://www.aoml.noaa.gov/ftp/pub/hrd/data/radar/level3/) at NOAA Hurricane Research Division's website.

The original motivation for this project was to organize the set of functions used for analysis and vizulaization of the tail doppler radar analyses into one project. The consolidation of all the tools used makes it easy to use and reference.

## Installation

The TC-RADAR Toolkit can be installed for use by simply cloning this repository to your system.

```bash
git clone https://github.com/BrunoRojasWx/TCRADAR_PKG
```

Alternatively, a zipped folder of the repository can be downloaded from the green "Code" dropdown menu button here on the github page.

## Usage

TC-RADAR Toolkit makes it easy to read and access variables in TC-RADAR and access metadata about the mission of interest. This toolkit is designed to be used as a set of functions to perform calculations on the radar data, such as azimuthal means and shear-relative analyses.

A more compresenive set of examples is available in [usage_examples_tutorial.py](https://github.com/BrunoRojasWx/TCRADAR_PKG/blob/master/usage_examples_tutorial.py).

## Roadmap and Features

Currently the TC-RADAR Toolkit is in active development and not all features are implemented. 

Features present in the TC-RADAR Toolkit include:
- [x] Accessing merged analyses (mission scale)
- [x] Pulling merged analyses variables
- [x] Accessing datetime of the mission
- [x] Accessing SHIPS shear information at the time of the mission
- [x] Calculate azimuthal means of any variable (with variable radial bin size)
- [x] Correct reflectivity averaging by conversion to Z and back to dBZ
- [x] Shear-relative quadrant and half azimuthal means

We hope to implement the list of features below:

-Swath analyses function class
    -swath mean slices
    -contoured frequency by altitude diagrams (CFADs)
    -vetrical mass transport
    -CFADs by quadrants
    -CFAD anomalies
    -convective stratiform separation algorithm

-Cross Sections
-Storm context information using SHIPS data (shear history, intensity, relative humidity)
-Vortex tilt
-Wavenumber analysis
-Data vizualization class to plot figures of the data natively

## Contributing

Contributions from other meteorologists, hurricane scientists, and researchers are welcome!

## License

[MIT](https://choosealicense.com/licenses/mit/)