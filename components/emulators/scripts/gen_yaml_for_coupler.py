import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Optional

from spel.scripts.fortran_parser.spel_ast import (
    DoLoop,
    Expression,
    ExpressionStatement,
    FuncExpression,
    Identifier,
    IfConstruct,
    InfixExpression,
    Statement,
    StringLiteral,
    SubCallStatement,
)
from spel.scripts.fortran_parser.spel_parser import Parser
from spel.scripts.types import LineTuple, LogicalLineIterator

FieldtoComponentMap = dict[str, list[str]]


# @dataclass
# class VariableStack:
#     name: str
#     values: list[Any] = field(default_factory=list)
#
#     def pop(self):
#         assert len(self.values) > 0
#         return self.values.pop()

VariableStack = dict[str, list[Any]]


@dataclass
class MetadataEntry:
    attname: str
    longname: str
    stdname: str
    units: str


def check_expression_statement(
    expr: Expression,
    vtable: VariableStack,
):
    assert isinstance(
        expr, InfixExpression
    ), f"Unexpected expression statement:\n{expr}"

    if expr.operator == "=":
        assert isinstance(expr.left_expr, Identifier)
        key = expr.left_expr.value
        val = expr.right_expr
        vtable[key].append(str(val))
    else:
        sys.exit(f"Unexpected Infix operator in Expression Statement\n{expr}")

    return


def add_field(
    func_expr: FuncExpression,
    fmap: FieldtoComponentMap,
    vtable: VariableStack,
):
    """

    Ex calls:
    call seq_flds_add(dom_other,'frac')
    call seq_flds_add(a2x_states,"Sa_z")
    call seq_flds_add(x2l_states,"Sa_z")

    result after all 3:
        fmap = {'frac': ['dom_other'], 'Sa_z':['a2x_states','x2l_states']}
    """
    # assumptions
    assert len(func_expr.args) == 2, f"Encountered optional args:\n{func_expr}"

    component_arg, field_arg = func_expr.args

    assert isinstance(
        component_arg, Identifier
    ), f"Arguments took unexpected form\n{func_expr}"

    field_str = get_arg_value(field_arg, vtable)
    fmap[field_str].append(str(component_arg))
    return


def get_arg_value(_arg: Expression, vtable: VariableStack) -> str:
    if isinstance(_arg, Identifier):
        if str(_arg) not in vtable:
            raise ValueError(f"{_arg} not found in variable table")
        return vtable[str(_arg)][-1]
    elif isinstance(_arg, StringLiteral):
        return str(_arg)
    elif isinstance(_arg, InfixExpression):
        assert _arg.operator == "//"
        return str(_arg.left_expr) + "_{i}"
    elif isinstance(_arg, FuncExpression):
        fn = _arg.function.value
        assert fn == "trim"
        return get_arg_value(_arg.args[0], vtable)
    raise ValueError(f"Arg is unexpected type: {_arg}")


def get_metadata(
    expr: FuncExpression,
    vtable: VariableStack,
) -> MetadataEntry:
    """
    Ex calls:
        call metadata_set(attname, longname, stdname, units)

    """
    assert (
        len(expr.args) == 4
    ), f"Unexpected number of arguments to metadata_set:\n{expr}"
    arg = expr.args[0]
    assert isinstance(arg, Identifier), str(expr)
    attname = get_arg_value(arg, vtable) 

    arg = expr.args[1]
    longname = get_arg_value(arg, vtable)
    arg = expr.args[2]
    stdname = get_arg_value(arg, vtable)
    arg = expr.args[3]
    units = get_arg_value(arg, vtable)

    return MetadataEntry(
        attname=attname,
        longname=longname,
        stdname=stdname,
        units=units,
    )


def main():
    fn = "seq_flds_mod.F90"
    with open(fn, "r") as ifile:
        lines = ifile.readlines()

    lpairs = [LineTuple(ln=i, line=line) for i, line in enumerate(lines)]
    parser = Parser(lines=lpairs)
    program = parser.parse_program()

    variable_table: VariableStack = defaultdict(list)
    field_to_component_map: dict[str, list[str]] = defaultdict(list)
    metadata_map: dict[str,MetadataEntry] = {}

    def walk(stmts: list[Statement]):
        for stmt in stmts:
            match stmt:
                case ExpressionStatement():
                    expr = stmt.expression
                    check_expression_statement(expr, variable_table)
                case SubCallStatement():
                    expr = stmt.function
                    func_name = expr.function.value
                    if func_name == "seq_flds_add":
                        add_field(expr, field_to_component_map, variable_table)
                    elif func_name == "metadata_set":
                        entry = get_metadata(expr, variable_table)
                        metadata_map[entry.attname] = entry
                    else:
                        print(f"skipping subroutine call to: {expr}")
                case DoLoop():
                    block_stmts = stmt.body.statements
                    walk(block_stmts)
                case IfConstruct():
                    # TODO: track variables guarded by ifs
                    # blockstmts = stmt.consequence.statements
                    # walk(blockstmts)
                    if stmt.else_ifs:
                        print("Need to add else ifs")
                    if stmt.else_:
                        print("need to add else")
        return

    walk(program.statements)

    from pprint import pprint

    pprint(field_to_component_map)
    # pprint(metadata_map)

    to_yaml: list[YAML_Entry] = []
    for field, components in field_to_component_map.items():
        f_entry = decode_field(field,components)
        if f_entry:
            metadata = metadata_map[field]
            to_yaml.append(YAML_Entry(metadata=metadata,field=f_entry))


@dataclass
class FieldEntry:
    kind: str
    name: str
    merge_type: str
    sources: str
    destinations: list[str] = field(default_factory=list)


@dataclass
class YAML_Entry:
    metadata: MetadataEntry
    field: FieldEntry


ABBREV_TO_COMPONENT_MAP = {
    "a": "atm",
    "l": "lnd",
    "o": "ocean",
    "i": "ice",
    "x": "cpl",
    "g": "glc",
    "r": "rof",
    "w": "wav",
    "z": "iac",
    "f":"atm",
}


def get_state_destinations(source: str, structures: list[str]) -> list[str]:
    dests: list[str] = []
    if source == "cpl":
        for s in structures:
            dests.append(ABBREV_TO_COMPONENT_MAP[s[0]])
    else:
        for s in structures:
            if s[0] == "x":
                dests.append(ABBREV_TO_COMPONENT_MAP[s[2]])
    return dests


def decode_field(field_name: str, structures: list[str]) -> Optional[FieldEntry]:
    if field_name[0] == "s" and '_' in field_name:
        kind = 'state'
    elif field_name[0] == "f" and '_' in field_name: 
        kind = "flux"
    else: 
        return None

    if kind == "state":
        assert field_name[2] == "_"
        source = ABBREV_TO_COMPONENT_MAP.get(field_name[1] )
        if not source:
            sys.exit(f"Error incorrect state field {field_name}")
        dests = get_state_destinations(source, structures)
    elif kind == 'flux':
        assert field_name[0] == "f" and field_name[4] == "_", f"Malformed Flux field {field_name}"
        source = ABBREV_TO_COMPONENT_MAP[field_name[3]]
        dests = [ABBREV_TO_COMPONENT_MAP[field_name[2]]]

    if source == "cpl":
        merge_type = "by_name"
    else:
        merge_type = "direct"

    return FieldEntry(
        kind=kind,
        name=field_name.split('_')[-1],
        sources=source,
        destinations=dests,
        merge_type=merge_type,
    )


if __name__ == "__main__":
    main()
