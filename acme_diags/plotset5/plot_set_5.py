import plot_set_5_vcs

def plot(reference, test, reference_regrid, test_regrid, parameter):

    if parameter.backend.lower() == 'vcs':
        plot_set_5_vcs.plot(reference, test, reference_regrid, test_regrid, parameter)

    else:
        print("Unknown backend '%s'" % parameter.backend)
        quit()

