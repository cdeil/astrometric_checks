"""
Basic class and checks for pointing table.

TODO: move to Gammapy.

Usage: py.test test_pointing.py
"""
from astropy.utils import lazyproperty
from astropy.units import Quantity
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz
from gammapy.time.utils import time_ref_from_dict


def _earth_location_from_dict(meta):
    """TODO: re-use `gammapy.data.utils._earth_location_from_dict`

    At the moment the keys are different: GEOLAT vs. ALTITUDE.
    Support both or just the one from the spec?
    """
    lon = Angle(meta['GEOLON'], 'deg')
    lat = Angle(meta['GEOLAT'], 'deg')
    height = Quantity(meta['GEOALT'], 'meter')
    return EarthLocation(lon=lon, lat=lat, height=height)


class PointingInfo(object):
    """IACT array pointing info.

    TODO: link to open spec.
    TODO: share code with `EventList` class.

    This class has many cached properties.
    Should be used as read-only.

    Parameters
    ----------
    table : `~astropy.table.Table`
        Table (with meta header info) on pointing
    """

    def __init__(self, table):
        self.table = table

    @classmethod
    def read(cls, filename, hdu=None):
        """Read from file.
        """
        if hdu is None:
            hdu = 'POINTING'

        table = Table.read(filename, hdu=hdu)
        return cls(table=table)

    def __str__(self):
        """Basic info."""
        ss = 'Pointing info:\n\n'
        ss += 'Location:     {}\n'.format(self.location.geodetic)
        m = self.table.meta
        ss += 'MJDREFI, MJDREFF, TIMESYS = {}\n'.format((m['MJDREFI'], m['MJDREFF'], m['TIMESYS']))
        ss += 'Time ref:     {}\n'.format(self.time_ref.fits)
        ss += 'Time ref:     {} MJD (TT)\n'.format(self.time_ref.mjd)
        sec = self.duration.to('second').value
        hour = self.duration.to('hour').value
        ss += 'Duration:     {} sec = {} hours\n'.format(sec, hour)
        ss += 'Table length: {}\n'.format(len(self.table))

        ss += '\nSTART:\n' + self._str_for_index(0) + '\n'
        ss += '\nEND:\n' + self._str_for_index(-1) + '\n'

        return ss

    def _str_for_index(self, idx):
        """Information for one point in the pointing table"""
        ss = 'Time:  {}\n'.format(self.time[idx].fits)
        ss += 'Time:  {} MJD (TT)\n'.format(self.time[idx].mjd)
        ss += 'RADEC: {} deg\n'.format(self.radec[idx].to_string())
        ss += 'ALTAZ: {} deg\n'.format(self.altaz[idx].to_string())
        return ss

    @lazyproperty
    def time_ref(self):
        """Time reference (`~astropy.time.Time`)"""
        # For debugging ... change TIMESYS
        # self.table.meta['TIMESYS'] = 'utc'
        return time_ref_from_dict(self.table.meta)

    @lazyproperty
    def duration(self):
        """Duration (`~astropy.time.Time`)"""
        return self.time[-1] - self.time[0]

    @lazyproperty
    def time(self):
        """Time array (`~astropy.time.Time`)"""
        met = Quantity(self.table['TIME'].astype('float64'), 'second')
        time = self.time_ref + met
        return time.tt

    @lazyproperty
    def location(self):
        """Location (`~astropy.coordinates.EarthLocation`)"""
        return _earth_location_from_dict(self.table.meta)

    @lazyproperty
    def radec(self):
        """RA / DEC position (`~astropy.coordinates.SkyCoord`)"""
        lon = self.table['RA_PNT'].astype('float64')
        lat = self.table['DEC_PNT'].astype('float64')
        return SkyCoord(lon, lat, unit='deg', frame='icrs')

    @lazyproperty
    def altaz_frame(self):
        """AltAz frame object."""
        return AltAz(obstime=self.time, location=self.location)

    @lazyproperty
    def altaz(self):
        """ALT / AZ position (`~astropy.coordinates.SkyCoord`)"""
        lon = self.table['AZ_PNT'].astype('float64')
        lat = self.table['ALT_PNT'].astype('float64')
        return SkyCoord(lon, lat, unit='deg', frame=self.altaz_frame)


class TestPointingInfo:
    """Test `PointingInfo` class.
    """

    def setup(self):
        self.pointing_info = PointingInfo.read('hess_event_list.fits')

    def test_str(self):
        ss = str(self.pointing_info)
        print(ss)
        assert 'todo' in ss


def check_pointing_info():
    pointing_info = PointingInfo.read('hess_event_list.fits')
    print(pointing_info)

    # Check altaz coordinates from HESS software against Astropy
    altaz_hess = pointing_info.altaz
    altaz_astropy = pointing_info.radec.transform_to(pointing_info.altaz_frame)

    sky_diff = altaz_astropy.separation(altaz_hess).to('arcsec')
    az_diff = (altaz_astropy.az - altaz_hess.az).to('arcsec')
    alt_diff = (altaz_astropy.alt - altaz_astropy.alt).to('arcsec')
    hourangle_diff = az_diff.to('hourangle')

    print('sky_diff:  ', sky_diff.min(), sky_diff.max())
    print(' az_diff:  ', az_diff.min(), az_diff.max())
    print('alt_diff:  ', alt_diff.min(), alt_diff.max())

    print('time diff: ', hourangle_diff.min(), hourangle_diff.max())

    # import IPython; IPython.embed()

"""
Current output of check_pointing_info.

Looks like the times are off by about one minute.
Changing to TIMESYS=UTC in the header gives a different time offset (-20 sec)


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
alt_diff:   0arcsec 0arcsec
time diff:  0h01m00.7656s 0h01m02.6397s
"""


if __name__ == '__main__':
    check_pointing_info()
