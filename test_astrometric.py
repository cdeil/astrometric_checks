"""Tests for the Astropy Astrometric class.
"""
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle


def copy_test_event_list():
    """Copy HESS event list as test data file.

    This isn't public data, so we cut out the photons from the Crab nebula.
    """
    filename = '/Users/deil/work/hess-host-analyses/checks/coordinates_check/run_0018406_std_fullEnclosure_eventlist.fits'
    table = Table.read(filename)
    source_pos = SkyCoord(table.meta['RA_OBJ'], table.meta['DEC_OBJ'], unit='deg')
    event_pos = SkyCoord(table['RA'], table['DEC'], unit='deg')
    sep = source_pos.separation(event_pos)
    mask = (sep > Angle(0.3, 'deg'))
    table = table[mask]
    filename = 'hess_event_list.fits'
    print('Writing {}'.format(filename))
    table.write(filename)


if __name__ == '__main__':
    copy_test_event_list()
