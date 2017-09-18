#!/usr/bin/env python
"""
Optimization algorithm for solving MILP from timing data.
"""

import sys
import copy
import logging
import operator
import importlib
try:
    import pulp
except ImportError, e:
    sys.stderr.write("pulp library not installed or located. "
                     "Try pip install [--user] pulp\n")
    raise e

logger = logging.getLogger(__name__)

def solver_factory(data):
    """
    load data either from a json file or dictionary
    """
    if data.has_key('description'):
        description = data['description']

    if data.has_key('totaltasks'):
        maxtasks = data['totaltasks']
    else:
        logger.critical("ERROR: totaltasks not found in data")
        raise KeyError("totaltasks not found in data")

    layout = data['layout']
    if layout in globals():
        solverclass = globals()[layout]
    else:
        sp = layout.rsplit('.', 1)
        if len(sp) > 1:
            try:
                layout_module = importlib.import_module(sp[0])
                layout = sp[1]
            except ImportError, e:
                logger.critical("ERROR: cannot import %s\n" % sp[0])
                raise e
        try:
            solverclass = getattr(layout_module, layout)
        except KeyError, e:
            logger.critical("ERROR: layout class %s not found\n" % layout)
            raise e

    solver = solverclass()

    for c in solver.get_required_components():
        assert(data.has_key(c), "ERROR: component %s not found in data" % c)

    solver.set_data(data)
    return solver

class ModelData:
    """
    Convert dictionary data entry into usable object
    """
    def __init__(self, name, model_dict):
        self.name = name
        self.blocksize = model_dict['blocksize']
        self.nthrds = model_dict['nthrds'][0]
        ntasks = copy.deepcopy(model_dict['ntasks'])
        cost = copy.deepcopy(model_dict['cost'])
        assert len(ntasks) == len(cost), "ntasks data not same length as cost for %s" % name
        # sort smallest ntasks to largest
        tup = zip(*sorted(zip(cost, ntasks),
                          key=operator.itemgetter(1)))
        self.cost = list(tup[0])
        self.ntasks = list(tup[1])
        for j in self.ntasks:
            if j % self.blocksize:
                logger.warning("WARNING: %s pe %d not divisible by "
                               "blocksize %d. Results may be invalid\n"
                               % (name, j, self.blocksize))

class OptimizeModel(object):
    def __init__(self):
        self.models = {}
        self.status = pulp.constants.LpStatusUndefined

    def set_data(self, data_dict):
        """
        Add data to the model.
        data_dict is dictionary of components with their data
        example: {'totaltasks':64
                  'ICE': {'ntasks': [2,4,8],
                          'costs':  [10.0,6.0,4.0],
                          'nthrds': [1,1,1],
                          'blocksize': 8}
                  'LND': {...}
                 }

        data is extrapolated as needed for n=1 and n=totaltasks
        """
        # get deep copy, because we need to divide ntasks by blocksize
        self.maxtasks = data_dict['totaltasks']

        for key in data_dict:
            if isinstance(data_dict[key], dict) and 'ntasks' in data_dict[key]:
                self.models[key] = ModelData(key, data_dict[key])


        # extrapolate for n=1 and n=maxtasks
        for m in self.models.values():
            m.extrapolated = [False] * len(m.cost)

            # add in data for ntasks=1 if not provided
            if m.ntasks[0] > 1:
                m.cost.insert(0, m.ntasks[0] * m.cost[0])
                m.ntasks.insert(0, 1)
                m.extrapolated.insert(0, True)

            # add in data for maxtasks if not available
            # assume same scaling factor as previous interval
            if m.ntasks[-1] < self.maxtasks:
                if m.cost[-2] <= 0.0:
                    factor = 1.0
                elif len(m.ntasks) > 1:
                    factor = (1.0 - m.cost[-1]/m.cost[-2]) / \
                             (1.0 - 1. * m.ntasks[-2] / m.ntasks[-1])
                else:
                    # not much information to go on ...
                    factor = 1.0
                m.cost.append(m.cost[-1] * (1.0 - factor +
                                        factor * m.ntasks[-1] / self.maxtasks))
                m.ntasks.append(self.maxtasks)
                m.extrapolated.append(True)

        self.check_requirements()
        self.status = pulp.constants.LpStatusUndefined

    def add_model_constraints(self):
        """
        Build constraints based on the cost vs ntask models
        This should be the same for any layout so is provided in base class
        Assumes cost variables are 'Txxx' and ntask variables are 'Nxxx'
        """
        for k in self.get_required_components():
            m = self.models[k]
            tk = 'T' + k.lower() # cost(time) key
            nk = 'N' + k.lower() # nprocs key
            for i in range(0, len(m.cost) - 1):
                slope = (m.cost[i+1] - m.cost[i]) / (1. * m.ntasks[i+1] - m.ntasks[i])
                self.constraints.append([self.X[tk] - slope * self.X[nk] >= \
                                         m.cost[i] - slope * m.ntasks[i],
                                         "T%s - %f*N%s >= %f" % \
                                         (k.lower(), slope, k.lower(),
                                          m.cost[i] - slope * m.ntasks[i])])
                if slope > 0:
                    logger.warning("WARNING: Nonconvex cost function for "
                                   "model %s. Please")
                    logger.warning("Review costs to ensure data is correct "
                                   "(--graph_models or --print_models)")
                    break
                if slope == 0:
                    break

        self.status = pulp.constants.LpStatusNotSolved

    def get_required_components(self):
        """
        Should be overridden by derived class
        """
        return []

    def check_requirements(self):
        """
        Check to make sure that each element of the subclass's list of
        required components has some data provided.
        """
        for r in self.get_required_components():
            if r not in self.models:
                logger.critical("Data for component %s not available" % r)

    def write_timings(self, fd=sys.stdout, level=logging.DEBUG):
        """
        Print out the data used for the ntasks/cost models.
        Can be used to check that the data provided to the
        model is reasonable. Also see graph_costs()
        """
        for k in self.models:
            m = self.models[k]
            message = "***%s***" % k
            if fd is not None:
                fd.write("\n" + message + "\n")
            logger.log(level, message)

            for i in range(len(m.cost)):
                extra = ""
                if m.extrapolated[i]:
                    extra = " (extrapolated)"
                message = "%4d: %f%s" % \
                           (m.ntasks[i], m.cost[i], extra)
                if fd is not None:
                    fd.write(message + "\n")
                logger.log(level, message)

    def graph_costs(self):
        """
        Use matplotlib to graph the ntasks/cost data.
        This provides a quick visual to check that the
        data used for the optimization is reasonable.

        If matplotlib is not available, nothing will happen
        """
        try:
            import matplotlib.pyplot as pyplot
        except ImportError, e:
            logger.info("matplotlib not found, skipping graphs")
            return

        nplots = len(self.models)
        nrows = (nplots + 1) / 2
        ncols = 2
        fig, ax = pyplot.subplots(nrows, ncols)
        row = 0; col = 0
        for k in self.models:
            m = self.models[k]
            p = ax[row, col]
            p.loglog(m.ntasks, m.cost, 'k-')
            for i in range(len(m.ntasks)):
                if not m.extrapolated[i]:
                    p.plot(m.ntasks[i], m.cost[i], 'bx')
                else:
                    p.plot(m.ntasks[i], m.cost[i], 'rx')
            p.set_title(m.name)
            p.set_xlabel('ntasks')
            p.set_ylabel('cost (s/mday)')
            p.set_xlim([1, self.maxtasks])
            row += 1
            if row == nrows:
                row = 0
                col += 1

        fig.suptitle("log-log plot of Cost/mday vs ntasks for designated "
                     "components.\nPerfectly scalable components would have a "
                     "straight line. Blue 'X's designate points\nfrom data, "
                     "red 'X's designate extrapolated data. Areas above the "
                     "line plots represent\nthe feasible region. Global "
                     "optimality of solution depends on the convexity of "
                     "these line plots.\nClose graph to continue on to solve.")
        fig.tight_layout()
        fig.subplots_adjust(top=0.75)
        logger.info("close graph window to continue")
        pyplot.show()

    def optimize(self):
        """
        Run the optimization.
        Must set self.status using LpStatus object
           1  "Optimal"       LpStatusOptimal
           0  "Not Solved"    LpStatusNotSolved
          -1  "Infeasible"    LpStatusInfeasible
          -2  "Unbounded"     LpStatusUnbounded
          -3  "Undefined"     LpStatusUndefined
        Returns LpStatus

        If solved, then solution will be stored in self.X dictionary, indexed
        by variable name. Suggested convention:
        'Tice', 'Tlnd', ... for cost per component
        'Nice', 'Nlnd', ... for ntasks per component
        'NBice', 'NBlnd', ... for number of blocks per component

        The default implementation of get_solution() returns a dictionary
        of these variable keys and their values.
        """
        raise NotImplementedError

    def get_solution(self):
        """
        Return a dictionary of the solution variables, can be overridden
        """
        retval = {}
        for k in self.X:
            retval[k] = self.X[k].varValue
        return retval

    def get_all_variables(self):
        return self.prob.X

    def write_pe_file(self, pefilename):
        raise NotImplementedError

    def write_xml_changes(self, outfile):
        """
        Write out a list of xmlchange commands to implement
        the optimal layout
        """
        raise NotImplementedError

    def write_pe_template(self, pefilename, ntasks, nthrds, roots):
        from distutils.spawn import find_executable
        from xml.etree import ElementTree as ET
        from CIME.utils import run_cmd
        logger.info("Writing pe node info to %s" % pefilename)
        root = ET.Element('config_pes')
        grid = ET.SubElement(root, 'grid')
        grid.set('name', 'any')
        mach = ET.SubElement(grid, 'mach')
        mach.set('name', 'any')
        pes = ET.SubElement(mach, 'pes')
        pes.set('compset', 'any')
        pes.set('pesize', '')
        ntasks_node = ET.SubElement(pes, 'ntasks')
        for k in ntasks:
            node = ET.SubElement(ntasks_node, 'ntasks_' + k)
            node.text = str(ntasks[k])
        nthrds_node = ET.SubElement(pes, 'nthrds')
        for k in nthrds:
            node = ET.SubElement(nthrds_node, 'nthrds_' + k)
            node.text = str(nthrds[k])
        rootpe_node = ET.SubElement(pes, 'rootpe')
        for k in roots:
            node = ET.SubElement(rootpe_node, 'rootpe_' + k)
            node.text = str(roots[k])
        xmllint = find_executable("xmllint")
        import pdb; pdb.set_trace()
        if xmllint is not None:
            run_cmd("%s --format --output %s -" % (xmllint, pefilename),
                    input_str=ET.tostring(root))

class IceLndAtmOcn(OptimizeModel):
    """
    Optimized the problem based on the Layout
              ____________________
             | ICE  |  LND  |     |
             |______|_______|     |
             |              | OCN |
             |    ATM       |     |
             |______________|_____|

      Min T
      s.t.  T[ice]      <= T1
            T[lnd]      <= T1
            T1 + T[atm] <= T
            T[ocn]      <= T

            NB[c]        >= 1 for c in [ice,lnd,ocn,atm]
            NB[ice] + NB[lnd] <= NB[atm]
            atm_blocksize*NB[atm] + ocn_blocksize*NB[ocn] <= TotalTasks
            (NB[*] is number of processor blocks)

            T[c]        >= C[c]_{i} - NB[c]_{i} *
                       (C[c]_{i+1} - C[c]_{i}) / (NB[c]_{i+1} - NB[c]_{i})
                       + NB[c] * (C[c]_{i+1} - C[c]_{i})
                                               / (NB[c]_{i+1} - NB[c]_{i}),
                        i=1..ord(NB), c in [ice,lnd,ocn,atm]

    These assumptions are checked when solver is initialized
      . Assuming cost is monotonic decreasing vs ntasks
      . Assuming perfect scalability for ntasks < tasks[0]
      . Assuming same scalability factor for ntasks > ntasks[last] as for
                              last two data points

    Returns LpStatus
    If solved, then solution will be stored in self.X dictionary, indexed
    by variable name. Suggested convention:
        'Tice', 'Tlnd', ... for cost per component
        'Nice', 'Nlnd', ... for ntasks per component
        'NBice', 'NBlnd', ... for number of blocks per component
    """
    def get_required_components(self):
        return ['LND', 'ICE', 'ATM', 'OCN']

    def optimize(self):
        """
        Run the optimization.
        Must set self.status using LpStatus object
           1  "Optimal"       LpStatusOptimal
           0  "Not Solved"    LpStatusNotSolved
          -1  "Infeasible"    LpStatusInfeasible
          -2  "Unbounded"     LpStatusUnbounded
          -3  "Undefined"     LpStatusUndefined
        """
        assert(self.status == pulp.constants.LpStatusUndefined,
               "set_data() must be called before optimize()!")
        self.atm = self.models['ATM']
        self.lnd = self.models['LND']
        self.ice = self.models['ICE']
        self.ocn = self.models['OCN']
        self.real_variables = ['TotalTime', 'T1', 'Tice', 'Tlnd', 'Tatm',
                               'Tocn']
        self.integer_variables = ['NBice', 'NBlnd', 'NBatm', 'NBocn',
                                  'Nice', 'Nlnd', 'Natm', 'Nocn']
        self.X = {}
        X = self.X
        self.prob = pulp.LpProblem("Minimize ACME time cost", pulp.LpMinimize)
        for rv in self.real_variables:
            X[rv] = pulp.LpVariable(rv, lowBound=0)

        for iv in self.integer_variables:
            X[iv] = pulp.LpVariable(iv, lowBound=1, cat=pulp.LpInteger)


        # cost function
        self.prob += X['TotalTime']

        #constraints
        self.constraints = []
        # Layout-dependent constraints. Choosing another layout to model
        # will require editing these constraints
        self.constraints.append([X['Tice'] - X['T1'] <= 0, "Tice - T1 == 0"])
        self.constraints.append([X['Tlnd'] - X['T1'] <= 0, "Tlnd - T1 == 0"])
        self.constraints.append([X['T1'] + X['Tatm'] - X['TotalTime'] <= 0,
                                 "T1 + Tatm - TotalTime <= 0"])
        self.constraints.append([X['Tocn'] - X['TotalTime'] <= 0,
                                 "Tocn - TotalTime == 0"])
        self.constraints.append([X['Nice'] + X['Nlnd'] - X['Natm'] == 0,
                                 "Nice + Nlnd - Natm == 0"])
        self.constraints.append([X['Natm'] + X['Nocn'] == self.maxtasks,
                                 "Natm + Nocn <= %d" % (self.maxtasks)])
        self.constraints.append([self.atm.blocksize * X['NBatm'] - X['Natm'] == 0,
                                 "Natm = %d * NBatm" % self.atm.blocksize])
        self.constraints.append([self.ice.blocksize * X['NBice'] - X['Nice'] == 0,
                                 "Nice = %d * NBice" % self.ice.blocksize])
        self.constraints.append([self.lnd.blocksize * X['NBlnd'] - X['Nlnd'] == 0,
                                 "Nlnd = %d * NBlnd" % self.lnd.blocksize])
        self.constraints.append([self.ocn.blocksize * X['NBocn'] - X['Nocn'] == 0,
                                 "Nocn = %d * NBocn" % self.ocn.blocksize])

        # These are the constraints based on the timing data.
        # They should be the same no matter what the layout of the components.
        self.add_model_constraints()

        for c, s in self.constraints:
            self.prob += c, s

        # Write the program to file and solve (using coin-cbc)
        self.prob.writeLP("IceLndAtmOcn_model.lp")
        self.prob.solve()
        self.status = self.prob.status
        return pulp.LpStatus[self.prob.status]

    def get_solution(self):
        """
        Return a dictionary of the solution variables.
        """
        return {'NBLOCKS_ICE':self.X['NBice'].varValue,
                'NBLOCKS_LND':self.X['NBlnd'].varValue,
                'NBLOCKS_ATM':self.X['NBatm'].varValue,
                'NBLOCKS_OCN':self.X['NBocn'].varValue,
                'NTASKS_ICE':self.X['Nice'].varValue,
                'NTASKS_LND':self.X['Nlnd'].varValue,
                'NTASKS_ATM':self.X['Natm'].varValue,
                'NTASKS_OCN':self.X['Nocn'].varValue,
                'NTASKS_TOTAL':self.maxtasks,
                'COST_ICE':self.X['Tice'].varValue,
                'COST_LND':self.X['Tlnd'].varValue,
                'COST_ATM':self.X['Tatm'].varValue,
                'COST_OCN':self.X['Tocn'].varValue,
                'COST_TOTAL':self.X['TotalTime'].varValue}

    def write_pe_file(self, pefilename):
        """
        Write out a pe_file that can be used to implement the
        optimized layout
        """
        natm = int(self.X['Natm'].varValue)
        nlnd = int(self.X['Nlnd'].varValue)
        nice = int(self.X['Nice'].varValue)
        nocn = int(self.X['Nocn'].varValue)
        ntasks = {'atm':natm, 'lnd':nlnd, 'rof':1, 'ice':nice,
                  'ocn':nocn, 'glc':1, 'wav':1, 'cpl':1}
        roots = {'atm':0, 'lnd':nice, 'rof':0, 'ice':0,
                 'ocn':natm, 'glc':0, 'wav':0, 'cpl':0}
        nthrds = {}
        for c in ['atm', 'lnd', 'rof', 'ice', 'ocn', 'glc', 'wav', 'cpl']:
            if c.upper() in self.models:
                nthrds[c] = self.models[c.upper()].nthrds
            else:
                nthrds[c] = 1
        self.write_pe_template(pefilename, ntasks, nthrds, roots)
