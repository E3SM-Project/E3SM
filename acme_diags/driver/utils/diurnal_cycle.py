import MV2
import numpy
import numpy.ma as ma


def composite_diurnal_cycle(var, season, fft=True):
    """
    Compute the composite diurnal cycle for var for the given season.
    Return mean + amplitudes and times-of-maximum of the first Fourier harmonic component as three transient variables.
    """
    season_idx = {
        "01": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "02": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "03": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "04": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        "05": [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        "06": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        "07": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        "08": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        "09": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        "10": [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        "11": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        "12": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        "DJF": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        "MAM": [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
        "JJA": [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
        "SON": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
        "ANN": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    }

    site = False
    if var.getLatitude() is None and var.getLongitude() is None:
        site = True
        lat = var.lat
        lon = var.lon
    # Redefine time to be in the middle of the time interval
    var_time = var.getTime()
    if var_time is None:
        # Climo cannot be run on this variable.
        return var

    #    tbounds = var_time.getBounds()
    #    var_time[:] = 0.5*(tbounds[:,0]+tbounds[:,1]) #time bounds for h1-h4 are problematic
    var_time_absolute = var_time.asComponentTime()
    time_freq = int(
        24 / (var_time_absolute[1].hour - var_time_absolute[0].hour)
    )  # This only valid for time interval >= 1hour
    start_time = var_time_absolute[0].hour
    print("start_time", var_time_absolute[0], var_time_absolute[0].hour)
    print("var_time_freq={}".format(time_freq))

    # Convert to masked array
    v = var.asma()

    # Select specified seasons:
    if season == "ANNUALCYCLE":  # Not supported yet!
        cycle = [
            "01",
            "02",
            "03",
            "04",
            "05",
            "06",
            "07",
            "08",
            "09",
            "10",
            "11",
            "12",
        ]
    elif season == "SEASONALCYCLE":  # Not supported yet!
        cycle = ["DJF", "MAM", "JJA", "SON"]
    else:
        cycle = [season]

    ncycle = len(cycle)
    # var_diurnal has shape i.e. (ncycle, ntimesteps, [lat,lon]) for lat lon data
    var_diurnal = ma.zeros([ncycle] + [time_freq] + list(numpy.shape(v))[1:])
    for n in range(ncycle):
        # Get time index for each month/season.
        idx = numpy.array(
            [
                season_idx[cycle[n]][var_time_absolute[i].month - 1]
                for i in range(len(var_time_absolute))
            ],
            dtype=numpy.int,
        ).nonzero()
        var_diurnal[n,] = ma.average(
            numpy.reshape(
                v[idx],
                (int(v[idx].shape[0] / time_freq), time_freq)
                + v[idx].shape[1:],
            ),
            axis=0,
        )

    # Convert GMT to local time
    if site:
        nlat = 1
        nlon = 1
        # lat = [36.6]
        # lon = [262.5]
        lat = [
            lat,
        ]
        lon = [
            lon,
        ]
    else:
        nlat = var.shape[1]
        nlon = var.shape[2]
        lat = var.getLatitude()
        lon = var.getLongitude()
        var_diurnal = numpy.squeeze(var_diurnal)

    nt = time_freq
    lst = numpy.zeros((nt, nlat, nlon))
    for it, itime in enumerate(numpy.arange(0, 24, int(24 / nt))):
        for ilon in range(nlon):
            lst[it, :, ilon] = (
                itime + start_time + lon[ilon] / 360 * 24
            ) % 24  # convert GMT to LST

    # Compute mean, amplitude and max time of the first three Fourier components.
    if not fft:
        return var_diurnal, lst

    else:
        cycmean, maxvalue, tmax = fastAllGridFT(var_diurnal, lst)

        # Save phase, amplitude, and mean for the first homonic,
        amplitude = MV2.zeros((nlat, nlon))
        amplitude[:, :] = maxvalue[0]
        amplitude.id = var.id + "_diurnal_amplitude"
        amplitude.longname = "Amplitude of diurnal cycle of " + var.id
        amplitude.units = var.units
        amplitude.setAxis(0, lat)
        amplitude.setAxis(1, lon)

        maxtime = MV2.zeros((nlat, nlon))
        maxtime[:, :] = tmax[0]
        maxtime.id = var.id + "_diurnal_phase"
        maxtime.longname = "Phase of diurnal cycle of " + var.id
        maxtime.units = "hour"
        maxtime.setAxis(0, lat)
        maxtime.setAxis(1, lon)

        cmean = MV2.zeros((nlat, nlon))
        cmean[:, :] = cycmean
        cmean.id = var.id + "_diurnal_cycmean"
        cmean.longname = "Mean of diurnal cycle of " + var.id
        cmean.units = var.units
        cmean.setAxis(0, lat)
        cmean.setAxis(1, lon)

        return cmean, amplitude, maxtime


def fastAllGridFT(x, t):
    """
    This version of fastFT does all gridpoints at once.
    Use a Numerical Python function to compute a FAST Fourier transform -- which should give the same result as a simple
    SLOW Fourier integration via the trapezoidal rule.
    Return mean + amplitudes and times-of-maximum of the first three Fourier harmonic components of a time series x(t).
    Do NOT detrend the time series first, in order to retain the "sawtooth" frequency implied by the inumpy.t length of the
    time series (e.g. the 24-hour period from a composite-diurnal cycle).
    On inumpy.t: x[k,i,j] = values      at each gridpoint (i,j) for N times (k), e.g. N = 8 for a 3-hr composite-diurnal cycle
          t[k,i,j] = timepoints  at each gridpoint (i,j) for N times (k), e.g. Local Standard Times
    On output: c[i,j] = mean value at each gridpoint (i,j) in the time series ("zeroth" term in Fourier series)
           maxvalue[n,i,j] = amplitude       at each gridpoint (i,j) for each Fourier harmonic (n)
           tmax    [n,i,j] = time of maximum at each gridpoint (i,j) for each Fourier harmonic (n)
                Curt Covey, PCMDI/LLNL                                      December 2016
    """

    # Creating output arrays
    if len(x.shape) == 1:
        nx = 1
        ny = 1
    else:
        nx = x.shape[1]
        ny = x.shape[2]
    # time  of maximum for nth component (n=0 => diurnal, n=1 => semi...)
    tmax = numpy.zeros((3, nx, ny))
    # value of maximum for nth component (= 1/2 peak-to-peak amplitude)
    maxvalue = numpy.zeros((3, nx, ny))

    print(
        "Calling numpy FFT function and converting from complex-valued FFT to real-valued amplitude and phase"
    )
    X = numpy.fft.ifft(x, axis=0)
    print("FFT output shape={}".format(X.shape))

    # Converting from complex-valued FFT to real-valued amplitude and phase
    a = X.real
    b = X.imag
    S = numpy.sqrt(a ** 2 + b ** 2)
    c = S[0]  # Zeroth harmonic = mean-value "constant term" in Fourier series.
    for n in range(3):
        # Adding first + last terms, second + second-to-last, ...
        maxvalue[n] = S[n + 1] + S[-n - 1]
        tmax[n] = numpy.arctan2(b[n + 1], a[n + 1])
        tmax[n] = tmax[n] * 12.0 / (numpy.pi * (n + 1))  # Radians to hours
        tmax[n] = tmax[n] + t[0]  # GMT to LST
        tmax[n] = tmax[n] % (24 / (n + 1))
    return c, maxvalue, tmax
