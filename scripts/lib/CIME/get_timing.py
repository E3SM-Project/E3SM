#!/usr/bin/env python

"""
Library for implementing getTiming tool which gets timing
information from a run.
"""

from CIME.XML.standard_module_setup import *

import datetime, shutil, re

logger = logging.getLogger(__name__)

class _GetTimingInfo:
    def __init__(self, name):
        self.name = name
        self.tmin = 0
        self.tmax = 0
        self.adays = 0

class _TimingParser:
    def __init__(self, case, lid="999999-999999"):
        self.case = case
        self.caseroot = case.get_value("CASEROOT")
        self.lid = lid
        self.finlines = None
        self.fout = None
        self.adays=0
        self.models = {}

    def write(self, text):
        self.fout.write(text)

    def prttime(self, label, offset=None, div=None, coff=-999):
        if offset is None:
            offset=self.models['CPL'].offset
        if div is None:
            div = self.adays
        datalen = 20
        cstr = "<---->"
        clen = len(cstr)

        minval, maxval, found = self.gettime(label)
        if div >= 1.0:
            mind = minval/div
            maxd = maxval/div
        else:
            mind = minval
            maxd = maxval

        pstrlen = 25
        if mind >= 0 and maxd >= 0 and found:
            if coff >= 0:
                zoff = pstrlen + coff + int((datalen-clen)/2)
                csp = offset - coff - int((datalen-clen)/2)
                self.write(" {label:<{width1}}{cstr:<{width2}} "
                           "{minv:8.3f}:{maxv:8.3f} \n".format(label=label,
                                                               width1=zoff,
                                                               cstr=cstr,
                                                               width2=csp,
                                                               minv=mind,
                                                               maxv=maxd))
            else:
                zoff = pstrlen + offset
                self.write(" {label:<{width1}} {minv:8.3f}:{maxv:8.3f} \n"
                           .format(label=label, width1=zoff,
                                   minv=mind, maxv=maxd))

    def gettime2(self, heading_padded):
        nprocs = 0
        ncount = 0

        heading = '"' + heading_padded.strip() + '"'
        for line in self.finlines:
            m = re.match(r'\s*%s\s*(\d+)\s*\d+\s*(\S+)'%heading, line)
            if m:
                nprocs = int(float(m.groups()[0]))
                ncount = int(float(m.groups()[1]))
                return (nprocs, ncount)
            else:
                m = re.match(r'\s*%s\s+(\d+)\s'%heading, line)
                if m:
                    nprocs = 1
                    ncount = int(float(m.groups()[0]))
                    return (nprocs, ncount)
        return (0, 0)

    def gettime(self, heading_padded):
        found = False
        heading = '"' + heading_padded.strip() + '"'
        minval = 0
        maxval = 0

        for line in self.finlines:
            m = re.match(r'\s*%s\s*\d+\s*\d+\s*\S+\s*\S+\s*(\d*\.\d+)'
                         r'\s*\(.*\)\s*(\d*\.\d+)\s*\(.*\)'%heading, line)
            if m:
                maxval = float(m.groups()[0])
                minval = float(m.groups()[1])
                found = True
                return (minval, maxval, found)
        return (0, 0, False)

    def getTiming(self):
        components=self.case.get_values("COMP_CLASSES")
        for s in components:
            self.models[s] = _GetTimingInfo(s)
        atm = self.models['ATM']
        lnd = self.models['LND']
        rof = self.models['ROF']
        ice = self.models['ICE']
        ocn = self.models['OCN']
        glc = self.models['GLC']
        cpl = self.models['CPL']
        cime_model = self.case.get_value("MODEL")
        caseid = self.case.get_value("CASE")
        mach = self.case.get_value("MACH")
        user = self.case.get_value("USER")
        continue_run = self.case.get_value("CONTINUE_RUN")
        rundir = self.case.get_value("RUNDIR")
        run_type = self.case.get_value("RUN_TYPE")
        ncpl_base_period = self.case.get_value("NCPL_BASE_PERIOD")
        atm_ncpl = self.case.get_value("ATM_NCPL")
        ocn_ncpl = self.case.get_value("OCN_NCPL")
        compset = self.case.get_value("COMPSET")
        grid = self.case.get_value("GRID")
        run_type = self.case.get_value("RUN_TYPE")
        stop_option = self.case.get_value("STOP_OPTION")
        stop_n = self.case.get_value("STOP_N")
        cost_pes = self.case.get_value("COST_PES")
        totalpes = self.case.get_value("TOTALPES")
        pes_per_node = self.case.get_value("PES_PER_NODE")
        smt_factor = max(1,int(self.case.get_value("MAX_TASKS_PER_NODE") / pes_per_node))

        if cost_pes > 0:
            pecost = cost_pes
        else:
            pecost = totalpes

        for m in self.models.values():
            for key in ["NTASKS", "ROOTPE", "PSTRID", "NTHRDS", "NINST"]:
                if key == "NINST" and m.name == "CPL":
                    m.ninst = 1
                else:
                    setattr(m, key.lower(),
                            int(self.case.get_value("%s_%s" % (key, m.name))))

            m.comp = self.case.get_value("COMP_%s" % (m.name))
            m.pemax = m.rootpe + m.ntasks * m.pstrid - 1

        now = datetime.datetime.ctime(datetime.datetime.now())
        inittype = "FALSE"
        if (run_type == "startup" or run_type == "hybrid") and \
                not continue_run:
            inittype = "TRUE"

        binfilename = os.path.join(rundir, "timing", "model_timing_stats")
        finfilename = os.path.join(self.caseroot, "timing",
                                   "%s_timing_stats.%s" % (cime_model, self.lid))
        foutfilename = os.path.join(self.caseroot, "timing",
                                    "%s_timing.%s.%s" % (cime_model, caseid, self.lid))

        timingDir = os.path.join(self.caseroot, "timing")
        if not os.path.isdir(timingDir):
            os.makedirs(timingDir)

        try:
            shutil.copyfile(binfilename, finfilename)
        except Exception, e:
            if not os.path.isfile(binfilename):
                logger.critical("File %s not found" % binfilename)
            else:
                logger.critical("Unable to cp %s to %s" %
                                (binfilename, finfilename))
            raise e

        os.chdir(self.caseroot)
        try:
            fin = open(finfilename, "r")
            self.finlines = fin.readlines()
            fin.close()
        except Exception, e:
            logger.critical("Unable to open file %s" % finfilename)
            raise e

        tlen = 1.0
        if ncpl_base_period == "decade":
            tlen = 3650.0
        elif ncpl_base_period == "year":
            tlen = 365.0
        elif ncpl_base_period == "day":
            tlen = 1.0
        elif ncpl_base_period == "hour":
            tlen = 1.0/24.0
        else:
            logger.warning("Unknown NCPL_BASE_PERIOD=%s" % ncpl_base_period)

        nprocs, ncount = self.gettime2('CPL:CLOCK_ADVANCE ')
        nsteps = ncount / nprocs
        adays = nsteps*tlen/atm_ncpl
        odays = adays
        if inittype == "TRUE":
            odays = adays - (tlen/ocn_ncpl)
        peminmax = max([m.rootpe for m in self.models.values()])+1
        if ncpl_base_period in ["decade","year","day"] and int(adays) > 0:
            adays = int(adays)
            if tlen % ocn_ncpl == 0:
                odays = int(odays)
        self.adays = adays
        maxoffset = 40
        extraoff = 20
        for m in self.models.values():
            m.offset = int((maxoffset*m.rootpe)/peminmax) + extraoff
        cpl.offset = 0
        try:
            self.fout = open(foutfilename, "w")
        except Exception, e:
            logger.critical("Could not open file for writing: %s"
                            % foutfilename)
            raise e

        self.write("---------------- TIMING PROFILE ---------------------\n")

        self.write("  Case        : %s\n" % caseid)
        self.write("  LID         : %s\n" % self.lid)
        self.write("  Machine     : %s\n" % mach)
        self.write("  Caseroot    : %s\n" % self.caseroot)
        self.write("  Timeroot    : %s/Tools\n" % self.caseroot)
        self.write("  User        : %s\n" % user)
        self.write("  Curr Date   : %s\n" % now)

        self.write("  grid        : %s\n" % grid)
        self.write("  compset     : %s\n" % compset)
        self.write("  run_type    : %s, continue_run = %s "
                   "(inittype = %s)\n" %
                   (run_type, str(continue_run).upper(), inittype))
        self.write("  stop_option : %s, stop_n = %s\n" %
                   (stop_option, stop_n))
        self.write("  run_length  : %s days (%s for ocean)\n\n" %
                   (adays, odays))

        self.write("  component       comp_pes    root_pe   tasks  "
                   "x threads"
                   " instances (stride) \n")
        self.write("  ---------        ------     -------   ------   "
                   "------  ---------  ------  \n")
        maxthrds = 0
        for k in self.case.get_values("COMP_CLASSES"):
            m = self.models[k]
            self.write("  %s = %-8s   %-6u      %-6u   %-6u x %-6u  "
                       "%-6u (%-6u) \n"
                       % (m.name.lower(), m.comp, (m.ntasks*m.nthrds *smt_factor), m.rootpe,
                          m.ntasks, m.nthrds, m.ninst, m.pstrid))
            if m.nthrds > maxthrds:
                maxthrds = m.nthrds
        nmax  = self.gettime(' CPL:INIT ')[1]
        tmax  = self.gettime(' CPL:RUN_LOOP ')[1]
        wtmax = self.gettime(' CPL:TPROF_WRITE ')[1]
        fmax  = self.gettime(' CPL:FINAL ')[1]
        for k in components:
            if k != "CPL":
                m = self.models[k]
            m.tmin, m.tmax, _ = self.gettime(' CPL:%s_RUN ' %  m.name)
        otmin, otmax, _ = self.gettime(' CPL:OCNT_RUN ')

        # pick OCNT_RUN for tight coupling
        if otmax > ocn.tmax:
            ocn.tmin = otmin
            ocn.tmax = otmax

        cpl.tmin, cpl.tmax, _ = self.gettime(' CPL:RUN ')
        xmax = self.gettime(' CPL:COMM ')[1]
        ocnwaittime = self.gettime(' CPL:C2O_INITWAIT')[0]

        if odays != 0:
            ocnrunitime = ocn.tmax * (adays/odays - 1.0)
        else:
            ocnrunitime = 0.0


        correction = max(0, ocnrunitime - ocnwaittime)

        tmax = tmax + wtmax + correction
        ocn.tmax += ocnrunitime

        for m in self.models.values():
            m.tmaxr = 0
            if m.tmax > 0:
                m.tmaxr = adays*86400.0/(m.tmax*365.0)
        xmaxr = 0
        if xmax > 0:
            xmaxr = adays*86400.0/(xmax*365.0)
        tmaxr = 0
        if tmax > 0:
            tmaxr = adays*86400.0/(tmax*365.0)

        self.write("\n")
        self.write("  total pes active           : %s \n" % (totalpes*maxthrds*smt_factor ))
        self.write("  pes per node               : %s \n" % pes_per_node)
        self.write("  pe count for cost estimate : %s \n" % pecost)
        self.write("\n")

        self.write("  Overall Metrics: \n")
        self.write("    Model Cost:         %10.2f   pe-hrs/simulated_year \n"
                   % ((tmax*365.0*pecost)/(3600.0*adays)))
        self.write("    Model Throughput:   %10.2f   simulated_years/day \n"
                   % ((86400.0*adays)/(tmax*365.0)))
        self.write("\n")

        self.write("    Init Time   :  %10.3f seconds \n" % nmax)
        self.write("    Run Time    :  %10.3f seconds   %10.3f seconds/day \n"
                   % (tmax, tmax/adays))
        self.write("    Final Time  :  %10.3f seconds \n" % fmax)

        self.write("\n")

        self.write("    Actual Ocn Init Wait Time     :  %10.3f seconds \n"
                   % ocnwaittime)
        self.write("    Estimated Ocn Init Run Time   :  %10.3f seconds \n"
                   % ocnrunitime)
        self.write("    Estimated Run Time Correction :  %10.3f seconds \n"
                   % correction)
        self.write("      (This correction has been applied to the ocean and"
                   " total run times) \n")

        self.write("\n")
        self.write("Runs Time in total seconds, seconds/model-day, and"
                   " model-years/wall-day \n")
        self.write("CPL Run Time represents time in CPL pes alone, "
                   "not including time associated with data exchange "
                   "with other components \n")
        self.write("\n")

        self.write("    TOT Run Time:  %10.3f seconds   %10.3f seconds/mday   "
                   "%10.2f myears/wday \n" % (tmax, tmax/adays, tmaxr))
        for k in self.case.get_values("COMP_CLASSES"):
            m = self.models[k]
            self.write("    %s Run Time:  %10.3f seconds   "
                       "%10.3f seconds/mday   %10.2f myears/wday \n"
                       % (k, m.tmax, m.tmax/adays, m.tmaxr))
        self.write("    CPL COMM Time: %10.3f seconds   %10.3f seconds/mday   "
                   "%10.2f myears/wday \n" % (xmax, xmax/adays, xmaxr))

        self.write("\n\n---------------- DRIVER TIMING FLOWCHART "
                   "--------------------- \n\n")

        pstrlen = 25
        hoffset = 1
        self.write("   NOTE: min:max driver timers (seconds/day):   \n")

        for k in ["CPL", "OCN", "LND", "ROF", "ICE", "ATM", "GLC", "WAV"]:
            m = self.models[k]
            xspace = (pstrlen+hoffset+m.offset) * ' '
            self.write(" %s %s (pes %u to %u) \n"
                       % (xspace, k, m.rootpe, m.pemax))
        self.write("\n")

        self.prttime(' CPL:CLOCK_ADVANCE ')
        self.prttime(' CPL:OCNPRE1_BARRIER ')
        self.prttime(' CPL:OCNPRE1 ')
        self.prttime(' CPL:ATMOCN1_BARRIER ')
        self.prttime(' CPL:ATMOCN1 ')
        self.prttime(' CPL:OCNPREP_BARRIER ')
        self.prttime(' CPL:OCNPREP ')
        self.prttime(' CPL:C2O_BARRIER ', offset=ocn.offset, div=odays,
                     coff=cpl.offset)
        self.prttime(' CPL:C2O ', offset=ocn.offset, div=odays, coff=cpl.offset)
        self.prttime(' CPL:LNDPREP_BARRIER ')
        self.prttime(' CPL:LNDPREP ')
        self.prttime(' CPL:C2L_BARRIER ', offset=lnd.offset, coff=cpl.offset)
        self.prttime(' CPL:C2L ', offset=lnd.offset, coff=cpl.offset)
        self.prttime(' CPL:ICEPREP_BARRIER ')
        self.prttime(' CPL:ICEPREP ')
        self.prttime(' CPL:C2I_BARRIER ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:C2I ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:WAVPREP_BARRIER ')
        self.prttime(' CPL:WAVPREP ')
        self.prttime(' CPL:C2W_BARRIER ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:C2W ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:ROFPREP_BARRIER ')
        self.prttime(' CPL:ROFPREP ')
        self.prttime(' CPL:C2R_BARRIER ', offset=rof.offset, coff=cpl.offset)
        self.prttime(' CPL:C2R ', offset=rof.offset, coff=cpl.offset)
        self.prttime(' CPL:ICE_RUN_BARRIER ', offset=ice.offset)
        self.prttime(' CPL:ICE_RUN ', offset=ice.offset)
        self.prttime(' CPL:LND_RUN_BARRIER ', offset=lnd.offset)
        self.prttime(' CPL:LND_RUN ', offset=lnd.offset)
        self.prttime(' CPL:ROF_RUN_BARRIER ', offset=rof.offset)
        self.prttime(' CPL:ROF_RUN ', offset=rof.offset)
        self.prttime(' CPL:WAV_RUN_BARRIER ', offset=rof.offset)
        self.prttime(' CPL:WAV_RUN ', offset=rof.offset)
        self.prttime(' CPL:OCNT_RUN_BARRIER ', offset=ocn.offset, div=odays)
        self.prttime(' CPL:OCNT_RUN ', offset=ocn.offset, div=odays)
        self.prttime(' CPL:O2CT_BARRIER ', offset=ocn.offset, div=odays,
                     coff=cpl.offset)
        self.prttime(' CPL:O2CT ', offset=ocn.offset, div=odays,
                     coff=cpl.offset)
        self.prttime(' CPL:OCNPOSTT_BARRIER ')
        self.prttime(' CPL:OCNPOSTT ')
        self.prttime(' CPL:ATMOCNP_BARRIER ')
        self.prttime(' CPL:ATMOCNP ')
        self.prttime(' CPL:L2C_BARRIER ', offset=lnd.offset, coff=cpl.offset)
        self.prttime(' CPL:L2C ', offset=lnd.offset, div=cpl.offset)
        self.prttime(' CPL:LNDPOST_BARRIER ')
        self.prttime(' CPL:LNDPOST ')
        self.prttime(' CPL:GLCPREP_BARRIER ')
        self.prttime(' CPL:GLCPREP ')
        self.prttime(' CPL:C2G_BARRIER ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:C2G ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:R2C_BARRIER ', offset=rof.offset, coff=cpl.offset)
        self.prttime(' CPL:R2C ', offset=rof.offset, coff=cpl.offset)
        self.prttime(' CPL:ROFPOST_BARRIER ')
        self.prttime(' CPL:ROFPOST ')
        self.prttime(' CPL:BUDGET1_BARRIER ')
        self.prttime(' CPL:BUDGET1 ')
        self.prttime(' CPL:I2C_BARRIER ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:I2C ', offset=ice.offset, coff=cpl.offset)
        self.prttime(' CPL:ICEPOST_BARRIER ')
        self.prttime(' CPL:ICEPOST ')
        self.prttime(' CPL:FRACSET_BARRIER ')
        self.prttime(' CPL:FRACSET ')
        self.prttime(' CPL:ATMOCN2_BARRIER ')
        self.prttime(' CPL:ATMOCN2 ')
        self.prttime(' CPL:OCNPRE2_BARRIER ')
        self.prttime(' CPL:OCNPRE2 ')
        self.prttime(' CPL:C2O2_BARRIER ', offset=ocn.offset, div=odays,
                     coff=cpl.offset)
        self.prttime(' CPL:C2O2 ', offset=ocn.offset, div=odays,
                     coff=cpl.offset)
        self.prttime(' CPL:ATMOCNQ_BARRIER')
        self.prttime(' CPL:ATMOCNQ ')
        self.prttime(' CPL:ATMPREP_BARRIER ')
        self.prttime(' CPL:ATMPREP ')
        self.prttime(' CPL:C2A_BARRIER ', offset=atm.offset, coff=cpl.offset)
        self.prttime(' CPL:C2A ', offset=atm.offset, coff=cpl.offset)
        self.prttime(' CPL:OCN_RUN_BARRIER ', offset=ocn.offset, div=odays)
        self.prttime(' CPL:OCN_RUN ', offset=ocn.offset, div=odays)
        self.prttime(' CPL:ATM_RUN_BARRIER ', offset=atm.offset)
        self.prttime(' CPL:ATM_RUN ', offset=atm.offset)
        self.prttime(' CPL:GLC_RUN_BARRIER ', offset=glc.offset)
        self.prttime(' CPL:GLC_RUN ', offset=glc.offset)
        self.prttime(' CPL:W2C_BARRIER ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:W2C ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:WAVPOST_BARRIER ')
        self.prttime(' CPL:WAVPOST ', cpl.offset)
        self.prttime(' CPL:G2C_BARRIER ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:G2C ', offset=glc.offset, coff=cpl.offset)
        self.prttime(' CPL:GLCPOST_BARRIER ')
        self.prttime(' CPL:GLCPOST ')
        self.prttime(' CPL:A2C_BARRIER ', offset=atm.offset, coff=cpl.offset)
        self.prttime(' CPL:A2C ', offset=atm.offset, coff=cpl.offset)
        self.prttime(' CPL:ATMPOST_BARRIER ')
        self.prttime(' CPL:ATMPOST ')
        self.prttime(' CPL:BUDGET2_BARRIER ')
        self.prttime(' CPL:BUDGET2 ')
        self.prttime(' CPL:BUDGET3_BARRIER ')
        self.prttime(' CPL:BUDGET3 ')
        self.prttime(' CPL:BUDGETF_BARRIER ')
        self.prttime(' CPL:BUDGETF ')
        self.prttime(' CPL:O2C_BARRIER ', offset=ocn.offset,
                     div=odays, coff=cpl.offset)
        self.prttime(' CPL:O2C ', offset=ocn.offset, div=odays, coff=cpl.offset)
        self.prttime(' CPL:OCNPOST_BARRIER ')
        self.prttime(' CPL:OCNPOST ')
        self.prttime(' CPL:RESTART_BARRIER ')
        self.prttime(' CPL:RESTART')
        self.prttime(' CPL:HISTORY_BARRIER ')
        self.prttime(' CPL:HISTORY ')
        self.prttime(' CPL:TSTAMP_WRITE ')
        self.prttime(' CPL:TPROF_WRITE ')
        self.prttime(' CPL:RUN_LOOP_BSTOP ')

        self.write("\n\n")
        self.write("More info on coupler timing:\n")

        self.write("\n")
        self.prttime(' CPL:OCNPRE1 ')
        self.prttime(' CPL:ocnpre1_atm2ocn ')

        self.write("\n")
        self.prttime(' CPL:OCNPREP ')
        self.prttime(' CPL:OCNPRE2 ')
        self.prttime(' CPL:ocnprep_avg ')
        self.prttime(' CPL:ocnprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:LNDPREP ')
        self.prttime(' CPL:lndprep_atm2lnd ')
        self.prttime(' CPL:lndprep_mrgx2l ')
        self.prttime(' CPL:lndprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:ICEPREP ')
        self.prttime(' CPL:iceprep_ocn2ice ')
        self.prttime(' CPL:iceprep_atm2ice ')
        self.prttime(' CPL:iceprep_mrgx2i ')
        self.prttime(' CPL:iceprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:WAVPREP ')
        self.prttime(' CPL:wavprep_atm2wav ')
        self.prttime(' CPL:wavprep_ocn2wav ')
        self.prttime(' CPL:wavprep_ice2wav ')
        self.prttime(' CPL:wavprep_mrgx2w ')
        self.prttime(' CPL:wavprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:ROFPREP ')
        self.prttime(' CPL:rofprep_l2xavg ')
        self.prttime(' CPL:rofprep_lnd2rof ')
        self.prttime(' CPL:rofprep_mrgx2r ')
        self.prttime(' CPL:rofprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:GLCPREP ')
        self.prttime(' CPL:glcprep_avg ')
        self.prttime(' CPL:glcprep_lnd2glc ')
        self.prttime(' CPL:glcprep_mrgx2g ')
        self.prttime(' CPL:glcprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:ATMPREP ')
        self.prttime(' CPL:atmprep_xao2atm ')
        self.prttime(' CPL:atmprep_ocn2atm ')
        self.prttime(' CPL:atmprep_alb2atm ')
        self.prttime(' CPL:atmprep_ice2atm ')
        self.prttime(' CPL:atmprep_lnd2atm ')
        self.prttime(' CPL:atmprep_mrgx2a ')
        self.prttime(' CPL:atmprep_diagav ')

        self.write("\n")
        self.prttime(' CPL:ATMOCNP ')
        self.prttime(' CPL:ATMOCN1 ')
        self.prttime(' CPL:ATMOCN2 ')
        self.prttime(' CPL:atmocnp_ice2ocn ')
        self.prttime(' CPL:atmocnp_wav2ocn ')
        self.prttime(' CPL:atmocnp_fluxo ')
        self.prttime(' CPL:atmocnp_fluxe ')
        self.prttime(' CPL:atmocnp_mrgx2o ')
        self.prttime(' CPL:atmocnp_accum ')
        self.prttime(' CPL:atmocnp_ocnalb ')

        self.write("\n")
        self.prttime(' CPL:ATMOCNQ ')
        self.prttime(' CPL:atmocnq_ocn2atm ')
        self.prttime(' CPL:atmocnq_fluxa ')
        self.prttime(' CPL:atmocnq_atm2ocnf ')

        self.write("\n")
        self.prttime(' CPL:OCNPOSTT ')
        self.prttime(' CPL:OCNPOST ')
        self.prttime(' CPL:ocnpost_diagav ')

        self.write("\n")
        self.prttime(' CPL:LNDPOST ')
        self.prttime(' CPL:lndpost_diagav ')
        self.prttime(' CPL:lndpost_acc2lr ')
        self.prttime(' CPL:lndpost_acc2lg ')

        self.write("\n")
        self.prttime(' CPL:ROFOST ')
        self.prttime(' CPL:rofpost_diagav ')
        self.prttime(' CPL:rofpost_histaux ')
        self.prttime(' CPL:rofpost_rof2lnd ')
        self.prttime(' CPL:rofpost_rof2ice ')
        self.prttime(' CPL:rofpost_rof2ocn ')

        self.write("\n")
        self.prttime(' CPL:ICEPOST ')
        self.prttime(' CPL:icepost_diagav ')

        self.write("\n")
        self.prttime(' CPL:WAVPOST ')
        self.prttime(' CPL:wavpost_diagav ')

        self.write("\n")
        self.prttime(' CPL:GLCPOST ')
        self.prttime(' CPL:glcpost_diagav ')
        self.prttime(' CPL:glcpost_glc2lnd ')
        self.prttime(' CPL:glcpost_glc2ice ')
        self.prttime(' CPL:glcpost_glc2ocn ')

        self.write("\n")
        self.prttime(' CPL:ATMPOST ')
        self.prttime(' CPL:atmpost_diagav ')

        self.write("\n")
        self.prttime(' CPL:BUDGET ')
        self.prttime(' CPL:BUDGET1 ')
        self.prttime(' CPL:BUDGET2 ')
        self.prttime(' CPL:BUDGET3 ')
        self.prttime(' CPL:BUDGETF ')
        self.write("\n\n")

        self.fout.close()

def get_timing(case, lid):
    parser = _TimingParser(case, lid)
    parser.getTiming()
