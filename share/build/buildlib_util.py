
###############################################################################
def extract_from_macros(macro_dump, varname):
###############################################################################
    look_for = f"{varname} :="
    for line in macro_dump.splitlines():
        if look_for in line:
            return line.split(":=")[-1].strip()

    return ""

