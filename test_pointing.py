"""Checks of the pointing table from HESS.
"""
from gammapy.data import PointingInfo


def check_pointing_info():
    pointing_info = PointingInfo.read('hess_event_list.fits')
    print(pointing_info)

    # Check altaz coordinates from HESS software against Astropy
    altaz_hess = pointing_info.altaz
    altaz_astropy = pointing_info.radec.transform_to(pointing_info.altaz_frame)

    sky_diff = altaz_astropy.separation(altaz_hess).to('arcsec')
    az_diff = (altaz_astropy.az - altaz_hess.az).to('arcsec')
    alt_diff = (altaz_astropy.alt - altaz_hess.alt).to('arcsec')
    hourangle_diff = az_diff.to('hourangle')

    print('sky_diff:  ', sky_diff.min(), sky_diff.max())
    print(' az_diff:  ', az_diff.min(), az_diff.max())
    print('alt_diff:  ', alt_diff.min(), alt_diff.max())

    print('time diff: ', hourangle_diff.min(), hourangle_diff.max())

    # import IPython; IPython.embed()

"""
Current output of check_pointing_info.

The cause of the difference is unknown.

Discussion of this issue is here: https://github.com/gammasky/hess-host-analyses/issues/47

Pointing info:

Location:     (<Longitude 16.5002222222222 deg>, <Latitude -23.2717777777778 deg>, <Quantity 1834.999999999783 m>)
MJDREFI, MJDREFF, TIMESYS = (51910, 0.000742870370370241, 'TT')
Time ref:     2001-01-01T00:01:04.184(TT)
Time ref:     51910.00074287037 MJD (TT)
Duration:     1586.0000000044238 sec = 0.4405555555567844 hours
Table length: 100

START:
Time:  2004-01-21T19:50:02.184(TT)
Time:  53025.826414166666 MJD (TT)
RADEC: 83.6333 24.5144 deg
ALTAZ: 11.2043 41.3792 deg


END:
Time:  2004-01-21T20:16:28.184(TT)
Time:  53025.844770648146 MJD (TT)
RADEC: 83.6333 24.5144 deg
ALTAZ: 3.18474 42.1431 deg


sky_diff:   697.907arcsec 697.909arcsec
 az_diff:   911.484arcsec 939.596arcsec
alt_diff:   -137.97arcsec -40.2693arcsec
time diff:  0h01m00.7656s 0h01m02.6397s
"""


if __name__ == '__main__':
    check_pointing_info()
