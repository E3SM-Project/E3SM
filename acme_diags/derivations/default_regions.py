import cdutil

regions_specs = {
    "NHEX": {"domain": cdutil.region.domain(latitude=(30.0, 90, "ccb"))},
    "SHEX": {"domain": cdutil.region.domain(latitude=(-90.0, -30, "ccb"))},
    "TROPICS": {"domain": cdutil.region.domain(latitude=(-30.0, 30, "ccb"))},
    "global": {},
    "TRMM_region": {
        "domain": cdutil.region.domain(latitude=(-38.0, 38, "ccb"))
    },
    "90S50S": {"domain": cdutil.region.domain(latitude=(-90.0, -50, "ccb"))},
    "50S20S": {"domain": cdutil.region.domain(latitude=(-50.0, -20, "ccb"))},
    "20S20N": {"domain": cdutil.region.domain(latitude=(-20.0, 20, "ccb"))},
    "50S50N": {"domain": cdutil.region.domain(latitude=(-50.0, 50, "ccb"))},
    "5S5N": {"domain": cdutil.region.domain(latitude=(-5.0, 5, "ccb"))},
    "20N50N": {"domain": cdutil.region.domain(latitude=(20.0, 50, "ccb"))},
    "50N90N": {"domain": cdutil.region.domain(latitude=(50.0, 90, "ccb"))},
    "60S90N": {"domain": cdutil.region.domain(latitude=(-60.0, 90, "ccb"))},
    "ocean": {
        "value": 0.65,
    },
    "ocean_seaice": {
        "value": 0.65,
    },
    "land": {
        "value": 0.65,
    },
    "land_60S90N": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(-60.0, 90, "ccb")),
    },
    "ocean_TROPICS": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(-30.0, 30, "ccb")),
    },
    "land_NHEX": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(30.0, 90, "ccb")),
    },
    "land_SHEX": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(-90.0, -30, "ccb")),
    },
    "land_TROPICS": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(-30.0, 30, "ccb")),
    },
    "ocean_NHEX": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(30.0, 90, "ccb")),
    },
    "ocean_SHEX": {
        "value": 0.65,
        "domain": cdutil.region.domain(latitude=(-90.0, -30, "ccb")),
    },
    # follow AMWG polar range,more precise selector
    "polar_N": {"domain": cdutil.region.domain(latitude=(50.0, 90.0, "ccb"))},
    "polar_S": {
        "domain": cdutil.region.domain(latitude=(-90.0, -55.0, "ccb"))
    },
    # To match AMWG results, the bounds is not as precise in this case
    # 'polar_N_AMWG':{'domain': Selector(latitude=(50., 90.))},
    # 'polar_S_AMWG':{'domain': Selector(latitude=(-90., -55.))},
    # Below is for modes of variability
    "NAM": {
        "domain": cdutil.region.domain(
            latitude=(20.0, 90, "ccb"), longitude=(-180, 180, "ccb")
        )
    },
    "NAO": {
        "domain": cdutil.region.domain(
            latitude=(20.0, 80, "ccb"), longitude=(-90, 40, "ccb")
        )
    },
    "SAM": {
        "domain": cdutil.region.domain(
            latitude=(-20.0, -90, "ccb"), longitude=(0, 360, "ccb")
        )
    },
    "PNA": {
        "domain": cdutil.region.domain(
            latitude=(20.0, 85, "ccb"), longitude=(120, 240, "ccb")
        )
    },
    "PDO": {
        "domain": cdutil.region.domain(
            latitude=(20.0, 70, "ccb"), longitude=(110, 260, "ccb")
        )
    },
    # Below is for monsoon domains
    # All monsoon domains
    "AllM": {
        "domain": cdutil.region.domain(
            latitude=(-45.0, 45.0, "ccb"), longitude=(0.0, 360.0, "ccb")
        )
    },
    # North American Monsoon
    "NAMM": {
        "domain": cdutil.region.domain(
            latitude=(0.0, 45.0, "ccb"), longitude=(210.0, 310.0, "ccb")
        )
    },
    # South American Monsoon
    "SAMM": {
        "domain": cdutil.region.domain(
            latitude=(-45.0, 0.0, "ccb"), longitude=(240.0, 330.0, "ccb")
        )
    },
    # North African Monsoon
    "NAFM": {
        "domain": cdutil.region.domain(
            latitude=(0.0, 45.0, "ccb"), longitude=(310.0, 60.0, "ccb")
        )
    },
    # South African Monsoon
    "SAFM": {
        "domain": cdutil.region.domain(
            latitude=(-45.0, 0.0, "ccb"), longitude=(0.0, 90.0, "ccb")
        )
    },
    # Asian Summer Monsoon
    "ASM": {
        "domain": cdutil.region.domain(
            latitude=(0.0, 45.0, "ccb"), longitude=(60.0, 180.0, "ccb")
        )
    },
    # Australian Monsoon
    "AUSM": {
        "domain": cdutil.region.domain(
            latitude=(-45.0, 0.0, "ccb"), longitude=(90.0, 160.0, "ccb")
        )
    },
    # Below is for NINO domains.
    "NINO3": {
        "domain": cdutil.region.domain(
            latitude=(-5.0, 5.0, "ccb"), longitude=(210.0, 270.0, "ccb")
        )
    },
    "NINO34": {
        "domain": cdutil.region.domain(
            latitude=(-5.0, 5.0, "ccb"), longitude=(190.0, 240.0, "ccb")
        )
    },
    "NINO4": {
        "domain": cdutil.region.domain(
            latitude=(-5.0, 5.0, "ccb"), longitude=(160.0, 210.0, "ccb")
        )
    },
    # Below is for additional domains for diurnal cycle of precipitation
    "W_Pacific": {
        "domain": cdutil.region.domain(
            latitude=(-20.0, 20.0, "ccb"), longitude=(90.0, 180.0, "ccb")
        )
    },
    "CONUS": {
        "domain": cdutil.region.domain(
            latitude=(25.0, 50.0, "ccb"), longitude=(-125.0, -65.0, "ccb")
        )
    },
    "Amazon": {
        "domain": cdutil.region.domain(
            latitude=(-20.0, 5.0, "ccb"), longitude=(-80.0, -45.0, "ccb")
        )
    },
    # Below is for RRM(regionally refined model) domains.
    #'CONUS_RRM': {'domain': cdutil.region.domain(latitude=(20., 50., 'ccb'), longitude=(-125., -65., 'ccb'))},For RRM dataset, negative value won't work
    "CONUS_RRM": {
        "domain": cdutil.region.domain(
            latitude=(20.0, 50.0, "ccb"), longitude=(235.0, 295.0, "ccb")
        )
    },
    # Below is for debugging. A smaller latitude range reduces processing time.
    "DEBUG": {"domain": cdutil.region.domain(latitude=(-2.0, 2, "ccb"))},
}

points_specs = {
    # ARM sites coordinates, select nearest grid poit to ARM site coordinates
    # Each point is supplied with [latitude, longitude ,select method, description of the point]
    "sgp": [36.4, -97.5, "cob", "97.5W 36.4N Oklahoma ARM"],
    "nsa": [71.3, -156.6, "cob", "156.6W 71.3N Barrow ARM"],
    "twpc1": [-2.1, 147.4, "cob", "147.4E 2.1S Manus ARM"],
    "twpc2": [-0.5, 166.9, "cob", "166.9E 0.5S Nauru ARM"],
    "twpc3": [-12.4, 130.9, "cob", "130.9E 12.4S Darwin ARM"],
}
