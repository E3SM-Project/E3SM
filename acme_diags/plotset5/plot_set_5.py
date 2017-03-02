import plot_set_5_vcs
import plot_set_5_cartopy
#import plot_set_7_cartopy

def plot(reference, test, reference_regrid, test_regrid, parameter):

    if parameter.backend.lower() == 'vcs':
        plot_set_5_vcs.plot(reference, test, reference_regrid, test_regrid, parameter)

    elif parameter.backend.lower() == 'cartopy':
        plot_set_5_cartopy.plot(reference, test, reference_regrid, test_regrid, parameter)
        #plot_set_7_cartopy.plot(reference, test, reference_regrid, test_regrid, parameter, 'N')
        #plot_set_7_cartopy.plot(reference, test, reference_regrid, test_regrid, parameter, 'S')

    else:
        print("Unknown backend '%s'" % parameter.backend)
        quit()

