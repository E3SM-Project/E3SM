#! /usr/bin/env python

# Copyright (C) 2005, 2006, 2007, 2009, 2010
# Glimmer-CISM contributors - see AUTHORS file for list of contributors
#
# This file is part of Glimmer-CISM.
#
# Glimmer-CISM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# Glimmer-CISM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
#
# Glimmer-CISM is hosted on BerliOS.de:
# https://developer.berlios.de/projects/glimmer-cism/

# python script used to generate source code files given a variable definition file

import ConfigParser, sys, time, string,re, os.path

NOATTRIB = ['name','dimensions','dimlen','data','factor','load','hot','type','average','coordinates']
dimensions = {}
module = {}

AVERAGE_SUFFIX = 'tavg'

def dimid(name):
    return '%s_dimid'%name

def is_dimvar(var):
    """Return True if variable is associated with a dimension.

    this is assumed to be the case if no time dim is present
    """

    if len(string.split(var['dimensions'],',')) == 1 and 'dimlen' in var:
        return True
    else:
        return False
    
class Variables(dict):
    """Dictionary containing variable definitions."""

    def __init__(self,filename):
        """Initialise Variable class.

        filename: name of file containing variable definitions."""

        dict.__init__(self)

        # reading variable configuration file
        vars = ConfigParser.ConfigParser()
        vars.readfp(open(filename))

        self.__have_avg = False

        for v in vars.sections():
            if v == 'VARSET':
                for (name, value) in vars.items(v):
                    module[name]=value
                continue
            vardef = {}
            vardef['name'] = v
            for (name, value) in vars.items(v):
                vardef[name] = value
            if 'type' not in vardef:
                vardef['type'] = 'float'
            if 'average' in vardef:
                if vardef['average'].lower() in ['1','true','t']:
                    vardef['average'] = True
                    self.__have_avg = True
                else:
                    vardef['average'] = False
            else:
                vardef['average'] = False
            # handle dims
            for d in vardef['dimensions'].split(','):
                d=d.strip()
                if 'dimlen' in vardef:
                    dimensions[d] = vardef['dimlen']
                if d not in dimensions:
                    dimensions[d] = '-1'
            self.__setitem__(v,vardef)

            # handle average
            if vardef['average']:
                # copy data structure
                vardef_avg = vardef.copy()
                vardef_avg['average'] = False
                vardef_avg['load'] = '0'
                vardef_avg['data'] = '%s_%s'%(vardef_avg['data'],AVERAGE_SUFFIX)
                vardef_avg['name'] = '%s_%s'%(vardef_avg['name'],AVERAGE_SUFFIX)
                if 'long_name' in vardef_avg:
                    vardef_avg['long_name'] = '%s (time average)'%vardef_avg['long_name']
                if 'cell_methods' in vardef_avg:
                    vardef_avg['cell_methods'] = '%s time: mean over years'%vardef_avg['cell_methods']
                else:
                    vardef_avg['cell_methods'] = 'time: mean over years'
                vardef_avg['avg_factor'] = 'tavgf'
                # and add to dictionary
                self.__setitem__('%s_%s'%(v,AVERAGE_SUFFIX),vardef_avg)
                

    def keys(self):
        """Reorder standard keys alphabetically."""
        dk = []
        vk = []
        for v in dict.keys(self):
            if is_dimvar(self.__getitem__(v)):
                dk.append(v)
            else:
                vk.append(v)
        dk.sort()
        vk.sort()
        return dk+vk

    def get_avg(self):
        return self.__have_avg
    have_avg = property(get_avg)

class PrintVars:
    """Base class for printing variables."""
    canhandle = None
    comment = '!'

    def __init__(self,filename,outname=None):
        """Initialise.

        filename: name of file to be processed."""
        if os.path.basename(filename) != self.canhandle:
            raise NotImplementedError, 'Can only handle %s'%self.canhandle

        self.infile = open(filename,'r')
        if outname==None:
            self.stream = open(self.canhandle[:-3],'w')
        else:
            self.stream = open(outname,'w')

        self.handletoken = {'!GENVARS!' : self.print_var}

    def print_warning(self):
        """Write a warning message to stream"""

        self.stream.write("%s\n"%(80*self.comment))
        self.stream.write("%s WARNING: this file was automatically generated on\n%s %s\n%s from %s\n"%(self.comment,
                                                                                                       self.comment,time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()),
                                                                                                       self.comment, self.canhandle))
        self.stream.write("%s\n\n"%(80*self.comment))
        
    def print_var(self, var):
        """Template for writing single variable"""

        raise NotImplementedError, 'You should use one of the derived classes'

    def write(self,vars):
        """Merge file with definitions"""

        self.print_warning()
        for l in self.infile.readlines():
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

class PrintDoc(PrintVars):
    """Process varlist.tex"""
    canhandle = 'varlist.tex.in'
    comment = '%'

    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename,'%s_varlist.tex'%module['name'])
        
    def print_var(self, var):
        """Write single variable block to stream for ncdf_params."""

        # skip variables associated with dimension 
        load = ''
        if 'load' in var:
            if var['load'].lower() in ['1','true','t']:
                load = '$^\\ast$'

        self.stream.write("\\texttt{%s}%s & %s & %s\\\\\n"%(var['name'].replace('_','\_'),load,var['long_name'].replace('_','\_'),
                                                            var['units'].replace('_','\_')))
        if 'standard_name' in var:
            self.stream.write("&CF name: \\texttt{%s}&\\\\\n"%(var['standard_name'].replace('_','\_')))
        self.stream.write("\\hline\n")

class PrintNC_template(PrintVars):
    """Process ncdf_template.F90.in"""
    canhandle = 'ncdf_template.F90.in'
    
    def __init__(self,filename):
        """Initialise.

        filename: name of file to be processed."""

        PrintVars.__init__(self,filename,'%s_io.F90'%module['name'])
        self.numvars = 0
        self.handletoken['!GENVAR_VARDEF!'] = self.print_vardef
        self.handletoken['!GENVAR_WRITE!'] = self.print_var_write
        self.handletoken['!GENVAR_READ!'] = self.print_var_read
        self.handletoken['!GENVAR_ACCESSORS!'] = self.print_var_accessor
        self.handletoken['!GENVAR_CALCAVG!'] = self.print_var_avg_accu
        self.handletoken['!GENVAR_RESETAVG!'] = self.print_var_avg_reset

    def write(self,vars):
        """Merge ncdf.F90.in with definitions."""

        numvars = 0
        for v in vars:
            if vars[v]['dimensions'] != v:
                numvars = numvars + 1

        self.thisvar = 1

        self.print_warning()
        for l in self.infile.readlines():
            for k in module.keys():
                l = l.replace(k.upper(),module[k])
            for token in self.handletoken:
                if string.find(l,token) is not -1:
                    break
            if string.find(l,token) is not -1:
                for v in vars.keys():
                    self.handletoken[token](vars[v])
            elif '!GENVAR_DIMS!' in l:
                self.print_dimensions()
            elif '!GENVAR_CHECKDIM!' in l:
                self.print_checkdims()
            elif '!GENVAR_HAVE_AVG!' in l:
                self.print_have_avg(vars.have_avg)
            elif 'AVG_SUFF' in l:
                self.stream.write("%s"%l.replace('AVG_SUFF','\"_%s\"'%AVERAGE_SUFFIX))
            else:
                self.stream.write("%s"%l)
        self.infile.close()
        self.stream.close()

    def print_have_avg(self,have_avg):
        """define whether we have time averaged vars or not"""

        if have_avg:
            self.stream.write("#define HAVE_AVG 1\n")

        
    def print_vardef(self,var):
        """Write single variable block to stream for ncdf_file."""

        dims = string.split(var['dimensions'],',')
        dims.reverse()
        dimstring = dimid(dims[0].strip())
        for d in dims[1:]:
            dimstring = '%s, %s'%(dimstring,dimid(d.strip()))
        
        self.stream.write("    !     %s -- %s\n"%(var['name'],var['long_name'])) # writing comment
        spaces = 0
        idstring = 'varid'
        if not is_dimvar(var):
            spaces=3
            self.stream.write("    pos = index(NCO%%vars,' %s ')\n"%var['name'])
            self.stream.write("    status = parallel_inq_varid(NCO%%id,'%s',varid)\n"%var['name'])
            self.stream.write("    if (pos.ne.0) then\n")
            self.stream.write("      NCO%%vars(pos+1:pos+%d) = '%s'\n"%(len(var['name']),len(var['name'])*' '))
            self.stream.write("    end if\n")
            self.stream.write("    if (pos.ne.0 .and. status.eq.nf90_enotvar) then\n")
        else:
            spaces=3
            self.stream.write("    if (.not.outfile%append) then\n")
        self.stream.write("%s    call write_log('Creating variable %s')\n"%(spaces*' ',var['name']))
        self.stream.write("%s    status = parallel_def_var(NCO%%id,'%s',get_xtype(outfile,NF90_%s),(/%s/),%s)\n"%(spaces*' ',
                                                                                                              var['name'],
                                                                                                              var['type'].upper(), 
                                                                                                              dimstring,
                                                                                                              idstring
                                                                                                              ))
        self.stream.write("%s    call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces*' '))
        if 'factor' in var:
            if var['factor'] == 'noscale':
                self.stream.write("%s    status = parallel_put_att(NCO%%id, %s, 'scale_factor',1.0)\n"%(spaces*' ',idstring))
            else:
                self.stream.write("%s    status = parallel_put_att(NCO%%id, %s, 'scale_factor',(%s))\n"%(spaces*' ',idstring,var['factor']))
        for attrib in var:
            if attrib not in NOATTRIB:
                self.stream.write("%s    status = parallel_put_att(NCO%%id, %s, '%s', '%s')\n"%(spaces*' ',
                                                                                            idstring,
                                                                                            attrib,
                                                                                            var[attrib]))
        if not is_dimvar(var):
            self.stream.write("%s    if (glimmap_allocated(model%%projection)) then\n"%(spaces*' '))
            self.stream.write("%s       status = parallel_put_att(NCO%%id, %s, 'grid_mapping',glimmer_nc_mapvarname)\n"%(spaces*' ',idstring))
            attrib='coordinates'
            if attrib in var:
                self.stream.write("%s       status = parallel_put_att(NCO%%id, %s, '%s', '%s')\n"%(spaces*' ',idstring,attrib,var[attrib]))
            self.stream.write("%s    end if\n"%(spaces*' '))
            self.stream.write("%s  end if\n"%(spaces*' '))
        else:
            self.stream.write("%s  end if\n"%(spaces*' '))
        self.stream.write("\n")

    def print_dimensions(self):
        """Set up dimensions."""

        dims = dimensions.keys()
        dims.sort()
        # generate list of dimension ids
        for d in dims:
            self.stream.write("    integer :: %s_dimid\n"%d)
        # get dimension ids
        self.stream.write("\n    ! defining dimensions\n")
        for d in dims:
            if dimensions[d]!='-1': # create a new dimension
                self.stream.write("    if (.not.outfile%append) then\n")
                self.stream.write("       status = parallel_def_dim(NCO%%id,'%s',%s,%s)\n"%(d,dimensions[d],dimid(d)))
                self.stream.write("    else\n")
                self.stream.write("       status = parallel_inq_dimid(NCO%%id,'%s',%s)\n"%(d,dimid(d)))
                self.stream.write("    endif\n")
            else:
                self.stream.write("    status = parallel_inq_dimid(NCO%%id,'%s',%s)\n"%(d,dimid(d)))
            self.stream.write("    call nc_errorhandle(__FILE__,__LINE__,status)\n")

    def print_checkdims(self):
        """Produce code for checking dimension sizes"""

        dims = dimensions.keys()
        dims.sort()
        for d in dims:
            if dimensions[d]!='-1':
                self.stream.write("    status = parallel_inq_dimid(NCI%%id,'%s',dimid)\n"%(d))
                self.stream.write("    if (dimid.gt.0) then\n")
                self.stream.write("       status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)\n")
                self.stream.write("       if (dimsize.ne.%s) then\n"%dimensions[d])
                self.stream.write("          write(message,*) 'Error, reading file ',trim(NCI%%filename),' size %s does not match: ', &\n               %s\n"%(d,dimensions[d]))
                self.stream.write("          call write_log(message,GM_FATAL)\n")
                self.stream.write("       end if\n")
                self.stream.write("    end if\n")
                
        
    def print_var_write(self,var):
        """Write single variable block to stream for ncdf_file."""

        # skip variables associated with dimension 
        if not is_dimvar(var):
            dims = string.split(var['dimensions'],',')
            dims.reverse()
            for i in range(0,len(dims)):
                dims[i] = dims[i].strip()
            self.stream.write("    status = parallel_inq_varid(NCO%%id,'%s',varid)\n"%var['name'])
            self.stream.write("    if (status .eq. nf90_noerr) then\n")
            
            dimstring = ''
            spaces = ''
            for i in range(0,len(dims)):
                if i > 0:
                    dimstring = dimstring + ','
                if dims[i] == 'time':
                    dimstring = dimstring + 'outfile%timecounter'
                elif dims[i] == 'level':
                    dimstring = dimstring + 'up'

                #*sfp* added to deal w/ writing of vars associated w/ stag vert coord 
                elif dims[i] == 'staglevel':
                    dimstring = dimstring + 'up'
                #*MJH* added to deal w/ writing of vars associated w/ stag vert coord w/bnd
                elif dims[i] == 'stagwbndlevel':
                    dimstring = dimstring + 'up+1'  # goes to index up+1
                else:
                    dimstring = dimstring + '1'
                
            if  'level' in dims:
                # handle 3D fields
                spaces = ' '*3
                self.stream.write("       do up=1,NCO%nlevel\n")

            #*sfp* added to handle writing of vars associated w/ stag vert coord
            if  'staglevel' in dims:
                # handle 3D fields
                spaces = ' '*3
                self.stream.write("       do up=1,NCO%nstaglevel\n")
                        
            #*MJH* added to handle writing of vars associated w/ stag vert coord w/ bnd
            if  'stagwbndlevel' in dims:
                # handle 3D fields
                spaces = ' '*3
                self.stream.write("       do up=0,NCO%nstagwbndlevel\n")  # starts with index 0

            data = var['data']
            if 'avg_factor' in var:
                data = '(%s)*(%s)'%(var['avg_factor'],data)
            self.stream.write("%s       status = distributed_put_var(NCO%%id, varid, &\n%s            %s, (/%s/))\n"%(spaces,
                                                                                                               spaces,data, dimstring))
            self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))

            if  'level' in dims:
                self.stream.write("       end do\n")

            #*sfp* added to handle writing of vars associated w/ stag vert coord
            if  'staglevel' in dims:
                self.stream.write("       end do\n")

            #*MJH* added to handle writing of vars associated w/ stag vert coord w/ bnd
            if  'stagwbndlevel' in dims:
                self.stream.write("       end do\n")

            # remove self since it's not time dependent
            if 'time' not in dims:
                self.stream.write("       NCO%%do_var(%s) = .False.\n"%(var_type(var)))
                
            self.stream.write("    end if\n\n")

    def print_var_read(self,var):
        """Write single variable block to stream for reading netCDF data."""

        if 'load' in var:
            if var['load'].lower() in ['1','true','t']:
                dims = string.split(var['dimensions'],',')
                dims.reverse()
                for i in range(0,len(dims)):
                    dims[i] = dims[i].strip()
                self.stream.write("    status = parallel_inq_varid(NCI%%id,'%s',varid)\n"%var['name'])
                self.stream.write("    if (status .eq. nf90_noerr) then\n")
                self.stream.write("       call write_log('  Loading %s')\n"%var['name'])
                dimstring = ''
                spaces = ''
                for i in range(0,len(dims)):
                    if i > 0:
                        dimstring = dimstring + ','
                    if dims[i] == 'time':
                        dimstring = dimstring + 'infile%current_time'
                    elif dims[i] == 'level':
                        dimstring = dimstring + 'up'
                    #*sfp* added to deal w/ writing of vars associated w/ stag vert coord 
                    elif dims[i] == 'staglevel':
                        dimstring = dimstring + 'up'
                    #*MJH* added to deal w/ writing of vars associated w/ stag vert coord w/ bnd
                    elif dims[i] == 'stagwbndlevel':
                        dimstring = dimstring + 'up+1'   # goes to index up+1
                    else:
                        dimstring = dimstring + '1'

                if  'level' in dims:
                    # handle 3D fields
                    spaces = ' '*3
                    self.stream.write("       do up=1,NCI%nlevel\n")

                #*sfp* added to handle writing of vars associated w/ stag vert coord
                if  'staglevel' in dims:
                    # handle 3D fields
                    spaces = ' '*3
                    self.stream.write("       do up=1,NCI%nstaglevel\n")

                #*MJH* added to handle writing of vars associated w/ stag vert coord w/ bnd
                if  'stagwbndlevel' in dims:
                    # handle 3D fields
                    spaces = ' '*3
                    self.stream.write("       do up=0,NCI%nstagwbndlevel\n")  # starts at index 0

                self.stream.write("%s       status = distributed_get_var(NCI%%id, varid, &\n%s            %s, (/%s/))\n"%(spaces,
                                                                                                               spaces,var['data'], dimstring))
                self.stream.write("%s       call nc_errorhandle(__FILE__,__LINE__,status)\n"%(spaces))
                self.stream.write("%s       status = parallel_get_att(NCI%%id, varid,'scale_factor',scaling_factor)\n"%(spaces))
                self.stream.write("%s       if (status.ne.NF90_NOERR) then\n"%(spaces))
                if 'factor' in var:
                    self.stream.write("%s          scaling_factor = 1.0d0/(%s)\n"%(spaces,var['factor']))
                else:
                    self.stream.write("%s          scaling_factor = 1.0d0\n"%(spaces))
                if 'factor' in var:
                    self.stream.write("%s       else\n"%(spaces))
                    self.stream.write("%s          scaling_factor = scaling_factor/(%s)\n"%(spaces,var['factor']))
                self.stream.write("%s       end if\n"%(spaces))
                self.stream.write("%s       if (abs(scaling_factor-1.0d0).gt.1.d-17) then\n"%(spaces))
                self.stream.write("%s          call write_log(\"scaling %s\",GM_DIAGNOSTIC)\n"%(spaces,var['name']))
                self.stream.write("%s          %s = %s*scaling_factor\n"%(spaces,var['data'],var['data']))
                self.stream.write("%s       end if\n"%(spaces))

                if  'level' in dims:
                    self.stream.write("       end do\n")

                #*sfp* added to handle writing of vars associated w/ stag vert coord
                if  'staglevel' in dims:
                    self.stream.write("       end do\n")

                #*MJH* added to handle writing of vars associated w/ stag vert coord w/ bnd
                if  'stagwbndlevel' in dims:
                    self.stream.write("       end do\n")
                
                self.stream.write("    end if\n\n")

    def print_var_accessor(self,var):
        """Write accessor function to stream."""

        dims = string.split(var['dimensions'],',')
        dimlen = len(dims)-1
        if dimlen>0:
            dimstring = ", dimension(:"+",:"*(dimlen-1)+")"
        else:
            dimstring = ""
        if not is_dimvar(var) and dimlen<3 and AVERAGE_SUFFIX not in var['name']:
            # get
            self.stream.write("  subroutine %s_get_%s(data,outarray)\n"%(module['name'],var['name']))
            self.stream.write("    use glimmer_scales\n")
            self.stream.write("    use glimmer_paramets\n")
            self.stream.write("    use %s\n"%module['datamod'])
            self.stream.write("    implicit none\n")
            self.stream.write("    type(%s) :: data\n"%module['datatype'])
            if var['type'] == 'int':
                vtype = 'integer'
            else:
#WHL - Changing the default type to real(dp)
#                vtype = 'real'
                vtype = 'real(dp)'
            self.stream.write("    %s%s, intent(out) :: outarray\n\n"%(vtype,dimstring))
            if 'factor' in var:
                if var["factor"] == 'noscale':
                    self.stream.write("    outarray = %s\n"%(var['data']))
                else:
                    self.stream.write("    outarray = (%s)*(%s)\n"%(var['factor'], var['data']))
            else:
                self.stream.write("    outarray = %s\n"%(var['data']))
            self.stream.write("  end subroutine %s_get_%s\n\n"%(module['name'],var['name']))
            # set
            # only creating set routine if the variable is not derived
            if len(var['data'].split('data'))<3:
                self.stream.write("  subroutine %s_set_%s(data,inarray)\n"%(module['name'],var['name']))
                self.stream.write("    use glimmer_scales\n")
                self.stream.write("    use glimmer_paramets\n")
                self.stream.write("    use %s\n"%module['datamod'])
                self.stream.write("    implicit none\n")
                self.stream.write("    type(%s) :: data\n"%module['datatype'])
                if var['type'] == 'int':
                    vtype = 'integer'
                else:
#WHL - Changing the default type to real(dp)
#                    vtype = 'real'
                    vtype = 'real(dp)'
                self.stream.write("    %s%s, intent(in) :: inarray\n\n"%(vtype,dimstring))
                if 'factor' in var:
                    if var['factor'] == 'noscale':
                       self.stream.write("!  no rescaling here\n")
                    else:
                      self.stream.write("    %s = inarray/(%s)\n"%(var['data'], var['factor']))
                else:
                    self.stream.write("    %s = inarray\n"%(var['data']))
                self.stream.write("  end subroutine %s_set_%s\n\n"%(module['name'],var['name']))

    def print_var_avg_accu(self,var):
        """Take average of a single variable"""

        if var['average']:
            avgname = '%s_%s'%(var['name'],AVERAGE_SUFFIX)
            avgdata = '%s_%s'%(var['data'],AVERAGE_SUFFIX)
            self.stream.write("    ! accumulate %s\n"%var['name'])
            self.stream.write("    status = parallel_inq_varid(NCO%%id,'%s',varid)\n"%avgname)
            self.stream.write("    if (status .eq. nf90_noerr) then\n")
            self.stream.write("       %s = %s + factor * %s\n"%(avgdata,avgdata,var['data']))
            self.stream.write("    end if\n\n")

    def print_var_avg_reset(self,var):
        """Reset average variables"""

        if var['average']:
            avgname = '%s_%s'%(var['name'],AVERAGE_SUFFIX)
            avgdata = '%s_%s'%(var['data'],AVERAGE_SUFFIX)
            self.stream.write("    ! reset %s\n"%var['name'])
            self.stream.write("    status = parallel_inq_varid(NCO%%id,'%s',varid)\n"%avgname)
            self.stream.write("    if (status .eq. nf90_noerr) then\n")
            self.stream.write("       %s = 0.\n"%avgdata)
            self.stream.write("    end if\n\n")

def usage():
    """Short help message."""

    print 'Usage generate_ncvars.py vardef [outfile0.in [,... [,outfileN.in]]]'
    print 'generate source code files given a variable definition file'
    print ''
    print 'vardef: file containing variable definition'
    print 'outfile.in: output template to be processed'
    print 'print variables if no templates are given'

HandleFile={}
HandleFile['ncdf_template.F90.in'] = PrintNC_template
HandleFile['varlist.tex.in'] = PrintDoc

if __name__ == '__main__':

    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    vars = Variables(sys.argv[1])

    if len(sys.argv) == 2:
        for v in vars.keys():
            print v
            for o in vars[v]:
                print '%s: %s'%(o, vars[v][o])
            print ''
            print module
    else:
        for f in sys.argv[2:]:
            base_f = os.path.basename(f)
            if base_f in HandleFile:
                handle = HandleFile[base_f](f)
                handle.write(vars)
            else:
                handle = PrintNCDF_IO(f)
                handle.write(vars)
                     
