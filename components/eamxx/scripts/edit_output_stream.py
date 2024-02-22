import pathlib

from utils import expect, ensure_yaml

ensure_yaml()
import yaml

###############################################################################
def generate_empty_yaml(filename,overwrite):
###############################################################################
    """
    Generate a yaml file with basic fields set, but containing empty or invalid values
    """
    file = pathlib.Path(filename).resolve()

    expect (overwrite or not file.exists(),
            "YAML file already exist. Re-run with -O/--overwrite to overwrite existing file")

    if file.exists():
        file.unlink()

    data = {}
    data["filename_prefix"] = "UNSET"
    data["Averaging Type"] = "INVALID"
    data["Fields"] = {}
    data["output_control"] = {}
    data["output_control"]["Frequency"] = -1
    data["output_control"]["frequency_units"] = "never"

    with open(file,'w') as fd:
        yaml.dump(data,fd,Dumper=yaml.SafeDumper,explicit_start=True,explicit_end=True,version=(1,2))

###############################################################################
def edit_output_stream_impl(filename,generate,overwrite,avg_type,freq_units,
                            freq,grid,fields,reset,io_grid,horiz_remap_file,vertical_remap_file):
###############################################################################
    """
    Apply the requested changes to the output stream yaml file
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
                data["Averaging Type"] = "INVALID"
            elif s=="freq":
                data["output_control"]["Frequency"] = -1
            elif s=="freq_units":
                data["output_control"]["frequency_units"] = "never"
            elif s=="horiz_remap_file":
                del data["horiz_remap_file"]
            elif s=="vert_remap_file":
                del data["vert_remap_file"]
            elif s=="fields":
                expect (grid is not None,
                        "Fields reset requested, but no grid name provided. Re-run with --grid GRID_NAME")
                if grid in data["Fields"].keys():
                    data["Fields"][grid]["Field Names"] = []

                    # Remove this grid if there are no other options set
                    if len(data["Fields"][grid])==1:
                        del data["Fields"][grid]
            elif s=="io-grid":
                expect (grid is not None,
                        "IO grid reset requested, but no grid name provided. Re-run with --grid GRID_NAME")
                if grid in data["Fields"].keys():
                    del data["Fields"][grid]["IO Grid Name"]

                    # Remove this grid if there's not other options set other than
                    # fields names, and field names is an empty list
                    if len(data["Fields"][grid])==1 and len(data["Fields"][grid]["Field Names"])==0:
                        del data["Fields"][grid]

    if avg_type is not None:
        data["Averaging Type"] = avg_type

    if freq is not None:
        data["output_control"]["Frequency"] = freq

    if freq_units is not None:
        data["output_control"]["frequency_units"] = freq_units

    if horiz_remap_file is not None:
        data["horiz_remap_file"] = horiz_remap_file

    if vertical_remap_file is not None:
        data["vertical_remap_file"] = vertical_remap_file

    if len(fields)>0 or io_grid is not None:
        expect (grid is not None,
                "Fields list specified, but no grid name provided. Re-run with --grid GRID_NAME")

        section = data["Fields"].setdefault(grid,{})
        if "Field Names" not in section.keys():
            section["Field Names"] = []

        fnames = section["Field Names"]
        fnames += fields
        fnames = list(set(fnames))
        section["Field Names"] = fnames

        if io_grid is not None:
            section["IO Grid Name"] = io_grid
            # If not already present, add an empty list of field names
            section.setdefault("Field Names",[])
        
        data["Fields"][grid] = section

    # We cannot do online remap (typically dyn->physGLL) if horiz or vert remap is used
    has_online_remap = False
    for k,v in data["Fields"].items():
        has_online_remap = has_online_remap or "IO Grid Name" in v.keys();
    has_vert_remap   = "vertical_remap_file" in data.keys()
    has_horiz_remap  = "horiz_remap_file" in data.keys()
    expect (not has_online_remap or (not has_vert_remap and not has_horiz_remap),
            "Cannot use online remap when horiz/vert remap is used.")

    with open(file,'w') as fd:
        yaml.dump(dict(data),fd,Dumper=yaml.SafeDumper,explicit_start=True,explicit_end=True,version=(1,2))
    
    return True
