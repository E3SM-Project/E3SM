import plot_set_5_vcs
import plot_set_5_cartopy

def plot(reference, test, reference_regrid, test_regrid, parameter):

    if parameter.backend.lower() == 'vcs':
        plot_set_5_vcs.plot(reference, test, reference_regrid, test_regrid, parameter)

    elif parameter.backend.lower() == 'cartopy':
        plot_set_5_cartopy.plot(reference, test, reference_regrid, test_regrid, parameter)

    else:
        print("Unknown backend '%s'" % parameter.backend)
        quit()

