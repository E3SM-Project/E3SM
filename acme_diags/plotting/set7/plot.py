import acme_diags.plotting.set7.set7vcs
import acme_diags.plotting.set7.set7cartopy

def plot(reference, test, reference_regrid, test_regrid, parameter):
    if parameter.backend.lower() == 'vcs':
        acme_diags.plotting.set7.set7vcs.plot(reference, test, reference_regrid, test_regrid, parameter)

    elif parameter.backend.lower() == 'cartopy':
        acme_diags.plotting.set7.set7cartopy.plot(reference, test, reference_regrid, test_regrid, parameter)

    else:
        print("Unknown backend '%s'" % parameter.backend)
        quit()
