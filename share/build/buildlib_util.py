
###############################################################################
def extract_from_macros(macro_dump, varname):
###############################################################################
    look_for = f"{varname} :="
    for line in macro_dump.splitlines():
        if line.startswith(look_for):
            return line.split(":=")[-1].strip()

    return ""

