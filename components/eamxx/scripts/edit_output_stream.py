import pathlib

from utils import expect, ensure_yaml

ensure_yaml()
import yaml

###############################################################################
def generate_empty_yaml(filename,overwrite):
###############################################################################
    """
    Generate a yaml file with basic fields set, but containing empty or invalid values
    >>> fname = "__test__.yaml"
    >>> generate_empty_yaml(fname,False)
    >>> generate_empty_yaml(fname,False)
    Traceback (most recent call last):
    SystemExit: ERROR: YAML file already exist. Re-run with -O/--overwrite to overwrite existing file
    >>> generate_empty_yaml(fname,True)
    >>> data = yaml.load(open(fname,'r'),Loader=yaml.SafeLoader)
    >>> len(data)
    4
    >>> data["filename_prefix"]
    'UNSET'
    >>> data["averaging_type"]
    'INVALID'
    >>> len(data["fields"])
    0
    >>> oc = data["output_control"]
    >>> len(oc)
    2
    >>> oc["frequency"]
    -1
    >>> oc["frequency_units"]
    'never'
    >>> # Clean up the file
    >>> pathlib.Path(fname).unlink()
    """
    file = pathlib.Path(filename).resolve()

    expect (overwrite or not file.exists(),
            "YAML file already exist. Re-run with -O/--overwrite to overwrite existing file")

    if file.exists():
        file.unlink()

    data = {}
    data["filename_prefix"] = "UNSET"
    data["averaging_type"] = "INVALID"
    data["fields"] = {}
    data["output_control"] = {}
    data["output_control"]["skip_t0_output"] = "false"
    data["output_control"]["frequency"] = -1
    data["output_control"]["frequency_units"] = "never"

    with open(file,'w') as fd:
        yaml.dump(data,fd,Dumper=yaml.SafeDumper,explicit_start=True,explicit_end=True,version=(1,2))

###############################################################################
def edit_output_stream_impl(filename,prefix=None,generate=False,overwrite=False,
                            avg_type=None,skip_t0_output=False,freq_units=None,freq=None,
                            grid=None,fields=[],reset=None,io_grid=None,
                            horiz_remap_file=None,vertical_remap_file=None):
###############################################################################
    """
    Apply the requested changes to the output stream yaml file
    >>> fname = '__test__.yaml'
    >>> # Create the file
    >>> edit_output_stream_impl(fname,generate=True,prefix='foo')
    >>> # Set some basic options, and then check
    >>> edit_output_stream_impl(fname,avg_type='max',freq_units='ndays',freq=10)
    >>> data = yaml.load(open(fname,'r'),Loader=yaml.SafeLoader)
    >>> data['filename_prefix']
    'foo'
    >>> data['averaging_type']
    'max'
    >>> data['output_control']['frequency']
    10
    >>> data['output_control']['frequency_units']
    'ndays'
    >>> # Set fields options, and then check
    >>> edit_output_stream_impl(fname,fields=['a','b'],grid='my_grid',io_grid='other_grid')
    >>> data = yaml.load(open(fname,'r'),Loader=yaml.SafeLoader)
    >>> f = data['fields']['my_grid']['field_names']
    >>> f.sort()
    >>> f
    ['a', 'b']
    >>> data['fields']['my_grid']['io_grid_name']
    'other_grid'
    >>> # No remap if online remap (io_grid_name) is set
    >>> edit_output_stream_impl(fname,horiz_remap_file='blah')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot use online remap and horiz/vert remap at the same time.
    >>> edit_output_stream_impl(fname,vertical_remap_file='blah')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot use online remap and horiz/vert remap at the same time.
    >>> # Remove io grid and fields
    >>> edit_output_stream_impl(fname,reset=['fields','io-grid'])
    Traceback (most recent call last):
    SystemExit: ERROR: fields reset requested, but no grid name provided. Re-run with --grid GRID_NAME
    >>> edit_output_stream_impl(fname,reset=['fields','io-grid'],grid='my_grid')
    >>> data = yaml.load(open(fname,'r'),Loader=yaml.SafeLoader)
    >>> 'my_grid' in data['fields'].keys()
    False
    >>> # Set remap options, and then check
    >>> edit_output_stream_impl(fname,horiz_remap_file='blah1',vertical_remap_file='blah2')
    >>> data = yaml.load(open(fname,'r'),Loader=yaml.SafeLoader)
    >>> data['horiz_remap_file']
    'blah1'
    >>> data['vertical_remap_file']
    'blah2'
    >>> # Clean up the file
    >>> pathlib.Path(fname).unlink()
    """

    if generate:
        generate_empty_yaml(filename,overwrite)

    file = pathlib.Path(filename).resolve()

    expect (file.exists(),
            "YAML file does not exist. Re-run with -g/--generate to create")

    data = yaml.load(open(file,"r"),Loader=yaml.SafeLoader)

    # Before adding new options, process all the reset requests
    if reset is not None:
        for s in reset:
            if s=="avg-type":
                data["averaging_type"] = "INVALID"
            elif s=="skip_t0_output":
                data["skip_t0_output"] = "false"
            elif s=="preifx":
                data["filename_prefix"] = "UNSET"
            elif s=="freq":
                data["output_control"]["frequency"] = -1
            elif s=="freq_units":
                data["output_control"]["frequency_units"] = "never"
            elif s=="horiz_remap_file":
                del data["horiz_remap_file"]
            elif s=="vert_remap_file":
                del data["vert_remap_file"]
            elif s=="fields":
                expect (grid is not None,
                        "fields reset requested, but no grid name provided. Re-run with --grid GRID_NAME")
                if grid in data["fields"].keys():
                    data["fields"][grid]["field_names"] = []

                    # Remove this grid if there are no other options set
                    if len(data["fields"][grid])==1:
                        del data["fields"][grid]
            elif s=="io-grid":
                expect (grid is not None,
                        "IO grid reset requested, but no grid name provided. Re-run with --grid GRID_NAME")
                if grid in data["fields"].keys():
                    del data["fields"][grid]["io_grid_name"]

                    # Remove this grid if there's not other options set other than
                    # fields names, and field names is an empty list
                    if len(data["fields"][grid])==1 and len(data["fields"][grid]["field_names"])==0:
                        del data["fields"][grid]

    if prefix is not None:
        data["filename_prefix"] = prefix

    if avg_type is not None:
        data["averaging_type"] = avg_type

    if skip_t0_output is not None:
        data["skip_t0_output"] = skip_t0_output

    if freq is not None:
        expect (freq.lstrip('-+').isnumeric() or freq=='hist_n',
                f"Invalid value '{freq}' for --freq. Valid options are\n"
                 " - an integer\n"
                 " - HIST_N\n")
        data["output_control"]["frequency"] = int(freq) if freq.lstrip('-+').isnumeric() else f"${{{freq.upper()}}}"

    if freq_units is not None:
        explicit  = ['nsteps','nsecs','nmins','nhours','ndays','nmonths','nyears']
        expect (freq_units in explicit or freq_units=='hist_option',
                f"Invalid value '{freq_units}' for --freq-units. Valid options are (case insensitive)\n"
                 " - explicit values: 'nsteps','nsecs','nmins','nhours','ndays','nmonths','nyears'\n"
                 " - CIME variables : 'HIST_OPTION'\n")

        data["output_control"]["frequency_units"] = freq_units if freq_units in explicit else f"${{{freq_units.upper()}}}"

    if horiz_remap_file is not None:
        data["horiz_remap_file"] = horiz_remap_file

    if vertical_remap_file is not None:
        data["vertical_remap_file"] = vertical_remap_file

    if len(fields)>0 or io_grid is not None:
        expect (grid is not None,
                "fields list specified, but no grid name provided. Re-run with --grid GRID_NAME")

        section = data["fields"].setdefault(grid,{})
        if "field_names" not in section.keys():
            section["field_names"] = []

        fnames = section["field_names"]
        fnames += fields
        fnames = list(set(fnames))
        section["field_names"] = fnames

        if io_grid is not None:
            section["io_grid_name"] = io_grid
            # If not already present, add an empty list of field names
            section.setdefault("field_names",[])

        data["fields"][grid] = section

    # We cannot do online remap (typically dyn->physgll) if horiz or vert remap is used
    has_online_remap = False
    for k,v in data["fields"].items():
        has_online_remap = has_online_remap or "io_grid_name" in v.keys();
    has_vert_remap   = "vertical_remap_file" in data.keys()
    has_horiz_remap  = "horiz_remap_file" in data.keys()
    expect (not has_online_remap or (not has_vert_remap and not has_horiz_remap),
            "Cannot use online remap and horiz/vert remap at the same time.")

    with open(file,'w') as fd:
        yaml.dump(dict(data),fd,Dumper=yaml.SafeDumper,explicit_start=True,explicit_end=True,version=(1,2))
