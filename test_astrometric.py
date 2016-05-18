"""Tests for the Astropy Astrometric class against the HESS software.

http://docs.astropy.org/en/latest/coordinates/matchsep.html?highlight=astrometric%20frame#astrometric-frames
https://github.com/astropy/astropy/pull/4909
https://github.com/astropy/astropy/issues/4931
https://github.com/astropy/astropy/pull/4941
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


def test_fov_radec():
    """Test FOV RADEC coordinate transformations.
    """
    # Set up test data and astrometric frame (a.k.a. FOV frame)
    # centered on the telescope pointing position
    table = Table.read('hess_event_list.fits')
    center = SkyCoord(table.meta['RA_PNT'], table.meta['DEC_PNT'], unit='deg')
    aframe = center.astrometric_frame()

    # Transform: RADEC -> FOV_RADEC
    event_radec = SkyCoord(table['RA'], table['DEC'], unit='deg')
    event_fov = event_radec.transform_to(aframe)
    table['FOV_RADEC_LON_ASTROPY'] = event_fov.data.lon.wrap_at('180 deg').to('deg')
    table['FOV_RADEC_LAT_ASTROPY'] = event_fov.data.lat.to('deg')
    table['FOV_RADEC_LON_DIFF'] = Angle(table['FOV_RADEC_LON_ASTROPY'] - table['FOV_RADEC_LON'], 'deg').to('arcsec')
    table['FOV_RADEC_LAT_DIFF'] = Angle(table['FOV_RADEC_LAT_ASTROPY'] - table['FOV_RADEC_LAT'], 'deg').to('arcsec')

    # Transform: FOV_RADEC -> RADEC    #
    event_fov = SkyCoord(table['FOV_RADEC_LON_ASTROPY'], table['FOV_RADEC_LAT_ASTROPY'], unit='deg', frame=aframe)
    event_radec = event_fov.transform_to('icrs')
    table['RA_ASTROPY'] = event_radec.data.lon.to('deg')
    table['DEC_ASTROPY'] = event_radec.data.lat.to('deg')
    table['RA_DIFF'] = Angle(table['RA_ASTROPY'] - table['RA'], 'deg').to('arcsec')
    table['DEC_DIFF'] = Angle(table['DEC_ASTROPY'] - table['DEC'], 'deg').to('arcsec')

    # Check results
    table.info('stats')
    """
    FOV_RADEC_LON_ASTROPY    0.0443400478952    1.58878237603   -35.1578248604   21.7907000503
    FOV_RADEC_LAT_ASTROPY   -0.0177905539829    1.57634999964   -28.7981936822     17.18667566
       FOV_RADEC_LON_DIFF     -2.07926024446   0.792609733345   -8.00939051063   6.52532690282
       FOV_RADEC_LAT_DIFF      -16.225735491   0.608440988983   -24.6793550923  -6.86899057409
               RA_ASTROPY      83.6820322695    1.72309019655    52.9889995666   106.496595676
              DEC_ASTROPY      24.4867638715    1.58698917572   -8.10738650931   41.7000435121
                  RA_DIFF  -0.00228730168484 0.00720260391234 -0.0173224899129 0.0125271212539
                 DEC_DIFF -0.000683806319612 0.00343228179454 -0.0148859721222 0.0125804299074

    Conclusions:
    * currently results for RADEC -> FOV_RADEC are off by this much for unknown reasons:
        -8 to +6 arcsec in LON
        -24 to -6 arcsec in LAT
    * the Astropy transformation does roundtrip with this accuracy
      * good enough for us ... won't investigate further.
      * could be limited by float32 and switching to float64 would improve accuracy?
      * 0.01 arcsec in LON and LAT
    """
    # import IPython; IPython.embed()


if __name__ == '__main__':
    # copy_test_event_list()
    test_fov_radec()
