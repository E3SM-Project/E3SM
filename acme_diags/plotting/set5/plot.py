import acme_diags.plotting.set5.set5vcs
import acme_diags.plotting.set5.set5cartopy

def plot(reference, test, reference_regrid, test_regrid, parameter):
    if parameter.backend.lower() == 'vcs':
        acme_diags.plotting.set5.set5vcs.plot(reference, test, reference_regrid, test_regrid, parameter)

    elif parameter.backend.lower() == 'cartopy':
        acme_diags.plotting.set5.set5cartopy.plot(reference, test, reference_regrid, test_regrid, parameter)

    else:
        print("Unknown backend '%s'" % parameter.backend)
        quit()
