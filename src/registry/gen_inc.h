// Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
// and the University Corporation for Atmospheric Research (UCAR).
//
// Unless noted otherwise source code is licensed under the BSD license.
// Additional copyright and license information can be found in the LICENSE file
// distributed with this code, or at http://mpas-dev.github.com/license.html
//
void gen_namelists(struct namelist *);
void gen_history_attributes(char * modelname, char * corename, char * version);
void gen_field_defs(struct group_list * groups, struct variable *, struct dimension *);
void gen_reads(struct group_list * groups, struct variable *, struct dimension *);
void gen_writes(struct group_list * groups, struct variable *, struct dimension *, struct namelist *);
