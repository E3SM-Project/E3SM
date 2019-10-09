Building CIME on an UNSUPPORTED local machine
---------------------------------------------

These directions are for a Mac OS X 10.9 or 10.10 laptop using
homebrew or macports to install the required software. The procedure
is similar for a linux workstation or cluster, you will just use
different package management tools to install the third party
libraries.

Setup
=====

  - install xcode, including the command line tools. Failure to
    install the command line tools is the most likely cause if you
    get an error about the compilers not being able to create
    executables.

  - install third party libraries from homebrew or macports.

      - home brew

        Install science tap : <https://github.com/Homebrew/homebrew-science>

            brew install gcc --without-multilib cmake mpich hdf5 --enable-fortran netcdf --enable-fortran


      - macports

            sudo port install mpich +gcc48 hdf5-18 +mpich netcdf-fortran +gcc48 +mpich cmake

        Note: If you see an error while running create_newcase that
        indicates perl can't find XML::LibXML, you may need to install
        p5-xml-libxml as well.


  - Some of the shell scripts used by cesm hard code "gmake" instead
    of using the GMAKE variable from env_build.xml. To work around
    this, you should install gnu make, or simply create a link from
    make to gmake in you path.

        mkdir -p ${HOME}/local/bin
        ln -s `whereis make` ${HOME}/local/bin/gmake
        cat >> ${HOME}/.bashrc <<EOF
        export PATH=${PATH}:${HOME}/local/bin
        EOF

  - Create a directory for the local copy of the cesm input data:

        mkdir -p ~/projects/cesm-inputdata

    Plan on downloading ~30 GB of data before the first build.

  - build cprnc :

    homebrew :

        cd ${CIMEROOT}/tools/cprnc
        ${CIMEROOT}/tools/configure --macros-format CMake
        mkdir build
	cd build
        cmake \
            -DCMAKE_Fortran_COMPILER=/usr/local/bin/mpif90 \
            -DHDF5_DIR=/usr/local \
            -DNetcdf_INCLUDE_DIR=/usr/local/include ..
        make


    macports :

        cd ${CESMROOT}/tools/cprnc
        mkdir build
        cd build
        cmake \
            -DCMAKE_Fortran_COMPILER=/opt/local/bin/mpif90-mpich-gcc48 \
            -DHDF5_DIR=/opt/local \
            -DNetcdf_INCLUDE_DIR=/opt/local/include ..
        make

  - copy the template directory:

        cp ${CESMROOT}/machines/userdefined_laptop_template/* ${HOME}/.cime

    NOTE: it is highly reccommend that you place ~/.cime under version
    control with your favorite tool.

  - make sure your machine
    is listed in config\_machines.xml, and optionally config\_compilers.xml,
    config\_pes.xml.

    - config_machines.xml must have a section for this machine and all
    the paths are correct. Especially note that the example puts input
    data, scratch, and archive into ${HOME}/projects

    - config\_compilers.xml - set the correct path's for your system.

    - config\_pes.xml

    NOTE: that the selected pe layout is order dependent, and the
    standard config_pes.xml is read first, then over ridden by this
    file.



First Use
=========

  - Use should be the same on a laptop or workstation as a supercomputer.

  - create a case:

        cd ~/projects/my-cesm-sandbox/scripts
        ./create_newcase \
            --machine durango \
            --compiler gnu \
            --case ../junk-1x1_brazil \
            --compset ICLM45CN \
            --res 1x1_brazil

  - run `case_setup`, `case.build`, `case.run` as normal.

  - Test a global simulation, use:
        -res f45_f45

 - Running a test suite or individual test should work like on any supported platform:

        ./create_test \
            --machine durango --compiler gnu \
            --xml-mach durango --xml-compiler gnu \
            --testroot ../junk-root \
            --xml-category aux_pop_obgc_se \
            --nobatch
