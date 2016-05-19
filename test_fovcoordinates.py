#! /usr/bin/env python
""" short script to test FoV coordinates from fits export """

import os
import sys
import numpy as np
from astropy.table import Table


# read in file and events hdu

if not len(sys.argv)==2:
    sys.exit("Usage: python test_coordinates.py event-fits-file.fits")

eventfile = sys.argv[1]

if not os.path.isfile(eventfile):
  raise RuntimeError('\n file %s not found ...' %(infile))
else:
  print 'reading file %s...' %eventfile
  print '   '

table = Table.read(eventfile, hdu='EVENTS')


# definitons and settings
atol=1e-6
atol_nomsys=1e-1

DEG_TO_RAD = np.pi/180.
RAD_TO_DEG = 180./np.pi


###########################################################################
#  test FOV_RADEC_LON/LAT
###########################################################################

# get coordinates of tangential RADec system
skyx_radec = table['SKYX_RADEC'].data
skyy_radec = table['SKYY_RADEC'].data

# get lon/lat fov coordinates in radec system
lon_radec = table['FOV_RADEC_LON'].data
lat_radec = table['FOV_RADEC_LAT'].data

# calculate vector from radec lon/lat coordinates
print 'Initialize vector from RADec lon/lat FoV coordinates'
x = np.cos(lat_radec*DEG_TO_RAD)*np.cos(lon_radec*DEG_TO_RAD)
y = np.cos(lat_radec*DEG_TO_RAD)*np.sin(lon_radec*DEG_TO_RAD)
z = np.sin(lat_radec*DEG_TO_RAD)

# test norm
print '>> test norm of initialized vector'
isnormed = np.allclose((x**2+y**2+z**2),1.)

if isnormed:
  print '>> succeeded...'
else:
  print '*** WARNING: could not verify norm' 

# calculate projection on tangential plane, use axis convention from tangential radec sys
print 'Calculate projection on tangential plane'
proj_x =  np.divide(-z,x) * RAD_TO_DEG
proj_y =  np.divide(y,x) * RAD_TO_DEG

# compare to coordinates of tangential RADec system
print 'Compare to coordinates of tangential radec system'
if np.allclose(skyx_radec, proj_x, atol=atol)and np.allclose(skyy_radec, proj_y, atol=atol):
  print'>> successfully verified lon/lat RADec FoV coordinates (agreement to coordinates of tangential RADec sys at level atol=%s)'%atol
else:
  print'*** WARNING: no agreement between lon/lat RADec FoV coordinates and coordinates of tangential RADEC sys at level atol=%s)'%atol
print '   '


###########################################################################
#  test FOV_RADEC_THETA/PHI
###########################################################################

# get theta/phi fov coordinates in radec system
theta_radec = table['FOV_RADEC_THETA'].data
phi_radec = table['FOV_RADEC_PHI'].data

# calculate vector from radec theta/phi coordinates
print 'Initialize vector from RADec offset FoV coordinates (THETA,PHI)'
x = np.sin(theta_radec*DEG_TO_RAD)*np.cos(phi_radec*DEG_TO_RAD)
y = np.sin(theta_radec*DEG_TO_RAD)*np.sin(phi_radec*DEG_TO_RAD)
z = np.cos(theta_radec*DEG_TO_RAD)

# test norm
print '>> test norm of initialized vector'
isnormed = np.allclose((x**2+y**2+z**2),1.)

if isnormed:
  print '>> succeeded...'
else:
  print '*** WARNING: could not verify norm' 

# calculate projection on tangential plane, use axis convention from tangential radec sys
print 'Calculate projection on tangential plane'
proj_x = np.divide(-y,z) * RAD_TO_DEG
proj_y = np.divide(-x,z) * RAD_TO_DEG

# compare to coordinates of tangential RADec system
print 'Compare to coordinates of tangential radec system'
if np.allclose(skyx_radec, proj_x, atol=atol)and np.allclose(skyy_radec, proj_y, atol=atol):
  print'>> successfully verified offset RADec FoV coordinates (agreement to coordinates of tangential RADec sys at level atol=%s)'%atol
else:
  print'*** WARNING: no agreement between offset RADec FoV coordinates and coordinates of tangential RADEC sys at level atol=%s)'%atol
print '   '
print '---------------------------------------------------------------------------------------------------------------------------'
print '   '


###########################################################################
#  test FOV_ALTAZ_LON/LAT
###########################################################################

# get coordinates of tangential RADec system
detx = table['DETX'].data
dety = table['DETY'].data

# get lon/lat fov coordinates in altaz system
lon_altaz = table['FOV_ALTAZ_LON'].data
lat_altaz = table['FOV_ALTAZ_LAT'].data

# calculate vector from radec lon/lat coordinates
print 'Initialize vector from AltAz lon/lat FoV coordinates'
x = np.cos(lat_altaz*DEG_TO_RAD)*np.cos(lon_altaz*DEG_TO_RAD)
y = np.cos(lat_altaz*DEG_TO_RAD)*np.sin(lon_altaz*DEG_TO_RAD)
z = np.sin(lat_altaz*DEG_TO_RAD)

# test norm
print '>> test norm of initialized vector'
isnormed = np.allclose((x**2+y**2+z**2),1.)

if isnormed:
  print '>> succeeded...'
else:
  print '*** WARNING: could not verify norm' 

# calculate projection on tangential plane, use axis convention from tangential radec sys
print 'Calculate projection on tangential plane, switch sign due to conventions of nominal system'
proj_x_lonlat = -(np.divide(-z,x) * RAD_TO_DEG)
proj_y_lonlat = -(np.divide(y,x) * RAD_TO_DEG)

# compare to coordinates of tangential RADec system
print 'Compare to coordinates of nominal system'
if np.allclose(detx, proj_x_lonlat, atol=atol_nomsys) and np.allclose(dety, proj_y_lonlat, atol=atol_nomsys):
  print'>> verified lon/lat AltAz FoV coordinates (agreement to coordinates of nominal sys at level atol=%s)'%atol_nomsys
else:
  print'*** WARNING: no agreement between lon/lat RADec FoV coordinates and coordinates of nominal sys at level atol=%s)'%atol_nomsys
print '   '


###########################################################################
#  test FOV_ALTAZ_THETA/PHI
###########################################################################

# get theta/phi fov coordinates in radec system
theta_altaz = table['FOV_ALTAZ_THETA'].data
phi_altaz = table['FOV_ALTAZ_PHI'].data

# calculate vector from radec theta/phi coordinates
print 'Initialize vector from AltAz offset FoV coordinates (THETA,PHI)'
x = np.sin(theta_altaz*DEG_TO_RAD)*np.cos(phi_altaz*DEG_TO_RAD)
y = np.sin(theta_altaz*DEG_TO_RAD)*np.sin(phi_altaz*DEG_TO_RAD)
z = np.cos(theta_altaz*DEG_TO_RAD)

# test norm
print '>> test norm of initialized vector'
isnormed = np.allclose((x**2+y**2+z**2),1.)

if isnormed:
  print '>> succeeded...'
else:
  print '*** WARNING: could not verify norm' 

# calculate projection on tangential plane, use axis convention from tangential radec sys
print 'Calculate projection on tangential plane, switch sign du to conventions of nominal system'
proj_x = -(np.divide(-y,z) * RAD_TO_DEG)
proj_y = -(np.divide(-x,z) * RAD_TO_DEG)

# compare to coordinates of tangential RADec system
print 'Compare to coordinates of nominal system'
if np.allclose(detx, proj_x, atol=atol_nomsys) and np.allclose(dety, proj_y, atol=atol_nomsys):
  print'>> verified offset AltAz FoV coordinates (agreement to coordinates of nominal sys at level atol=%s)'%atol_nomsys
else:
  print'*** WARNING: no agreement between offset AltAz FoV coordinates and coordinates of nominal sys at level atol=%s)'%atol_nomsys
print '   '



###########################################################################
#  test agreement between FoV AltAz coordinates directly
###########################################################################

print 'Test agreement of projected coordinates calculated from FOV_ALTAZ_LON/LAT and FOV_ALTAZ_THETA/PHI' 
print '(agreement to nominal system is only supposed at low level...)'
if np.allclose(proj_x, proj_x_lonlat, atol=atol) and np.allclose(proj_y, proj_y_lonlat, atol=atol):
  print '>> succesfully verified agreement between AltAz FoV coordinates at level %s'%atol
else:
  print '***WARNING: could not verify agreement at level atol=%s'%atol
print '   '

