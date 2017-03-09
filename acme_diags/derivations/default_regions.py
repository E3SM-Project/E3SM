import cdutil

regions_specs = {
    'NHEX': {'domain': cdutil.region.domain(latitude=(30., 90, 'ccb'))},
    'SHEX': {'domain': cdutil.region.domain(latitude=(-90., -30, 'ccb'))},
    'TROPICS': {'domain': cdutil.region.domain(latitude=(-30., 30, 'ccb'))},
    "global": {},
    'TRMM_region': {'domain': cdutil.region.domain(latitude=(-38., 38, 'ccb'))},
    '90S50S': {'domain': cdutil.region.domain(latitude=(-90., -50, 'ccb'))},
    '50S20S': {'domain': cdutil.region.domain(latitude=(-50., -20, 'ccb'))},
    '20S20N': {'domain': cdutil.region.domain(latitude=(-20., 20, 'ccb'))},
    '20N50N': {'domain': cdutil.region.domain(latitude=(20., 50, 'ccb'))},
    '50N90N': {'domain': cdutil.region.domain(latitude=(50., 90, 'ccb'))},
    'ocean': {'value': 0.65, },
    'land': {'value': 0.65, },
    'ocean_TROPICS': {'value': 0.65, 'domain': cdutil.region.domain(latitude=(-30., 30, 'ccb'))},
    #'land_NHEX': {'value': 100, 'domain': cdutil.region.domain(latitude=(30., 90, 'ccb'))},
    #'land_SHEX': {'value': 100, 'domain': cdutil.region.domain(latitude=(-90., -30, 'ccb'))},
    #'land_TROPICS': {'value': 100, 'domain': cdutil.region.domain(latitude=(-30., 30, 'ccb'))},
    #"land": {'value': 100, },
    #'ocean_NHEX': {'value': 0, 'domain': cdutil.region.domain(latitude=(30., 90, 'ccb'))},
    #'ocean_SHEX': {'value': 0, 'domain': cdutil.region.domain(latitude=(-90., -30, 'ccb'))},
    # Below is for modes of variability
    "NAM": {'domain': cdutil.region.domain(latitude=(20., 90, 'ccb'), longitude=(-180, 180, 'ccb'))},
    "NAO": {'domain': cdutil.region.domain(latitude=(20., 80, 'ccb'), longitude=(-90, 40, 'ccb'))},
    "SAM": {'domain': cdutil.region.domain(latitude=(-20., -90, 'ccb'), longitude=(0, 360, 'ccb'))},
    "PNA": {'domain': cdutil.region.domain(latitude=(20., 85, 'ccb'), longitude=(120, 240, 'ccb'))},
    "PDO": {'domain': cdutil.region.domain(latitude=(20., 70, 'ccb'), longitude=(110, 260, 'ccb'))},
    # Below is for monsoon domains
    # All monsoon domains
    'AllM': {'domain': cdutil.region.domain(latitude=(-45., 45., 'ccb'), longitude=(0., 360., 'ccb'))},
    # North American Monsoon
    'NAMM': {'domain': cdutil.region.domain(latitude=(0., 45., 'ccb'), longitude=(210., 310., 'ccb'))},
    # South American Monsoon
    'SAMM': {'domain': cdutil.region.domain(latitude=(-45., 0., 'ccb'), longitude=(240., 330., 'ccb'))},
    # North African Monsoon
    'NAFM': {'domain': cdutil.region.domain(latitude=(0., 45., 'ccb'), longitude=(310., 60., 'ccb'))},
    # South African Monsoon
    'SAFM': {'domain': cdutil.region.domain(latitude=(-45., 0., 'ccb'), longitude=(0., 90., 'ccb'))},
    # Asian Summer Monsoon
    'ASM': {'domain': cdutil.region.domain(latitude=(0., 45., 'ccb'), longitude=(60., 180., 'ccb'))},
    # Australian Monsoon
    'AUSM': {'domain': cdutil.region.domain(latitude=(-45., 0., 'ccb'), longitude=(90., 160., 'ccb'))},
}

default_regions = ['global', 'NHEX', 'SHEX', 'TROPICS']

