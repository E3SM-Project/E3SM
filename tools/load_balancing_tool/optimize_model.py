#!/usr/bin/env python
"""
Optimization algorithm for solving MILP from timing data.
"""

import sys
import copy
import logging
import operator
import importlib
from CIME.utils import expect
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
    expect(data.has_key('totaltasks'),"totaltasks not found in data")

    layout = data['layout']
    sp = layout.rsplit('.', 1)
    try:
        if len(sp) > 1:
            layout_module = importlib.import_module(sp[0])
            layout = sp[1]
        else:
            import layouts
            layout_module = layouts
    except ImportError:
        expect(False,"cannot import %s\n")

    try:
        solverclass = getattr(layout_module, layout)
    except KeyError:
        expect(False, "layout class %s not found in %s\n",
               layout, layout_module)

    solver = solverclass()

    for c in solver.get_required_components():
        assert data.has_key(c), "ERROR: component %s not found in data" % c

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
            if j > 1 and j % self.blocksize:
                logger.warning("WARNING: %s pe %d not divisible by "
                               "blocksize %d. Results may be invalid\n",
                               name, j, self.blocksize)

class OptimizeModel(object):
    STATE_UNDEFINED = 0
    STATE_UNSOLVED = 1
    STATE_SOLVED_OK = 2
    STATE_SOLVED_BAD = 3
    states = ['Undefined', 'Unsolved', 'Solved', 'No Solution']

    def __init__(self):
        self.models = {}
        self.state = self.STATE_UNDEFINED
        self.X = {}
        self.constraints = []
        self.maxtasks = 0

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
        sets state to STATE_UNSOLVED
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

            if len(m.ntasks) > 1 and m.ntasks[-1] < self.maxtasks:
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
        self.state = self.STATE_UNSOLVED

    def add_model_constraints(self):
        """
        Build constraints based on the cost vs ntask models
        This should be the same for any layout so is provided in base class
        Assumes cost variables are 'Txxx' and ntask variables are 'Nxxx'
        """
        assert self.state != self.STATE_UNDEFINED,\
               "set_data() must be called before add_model_constraints()"
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
                    logger.warning("WARNING: Nonconvex cost function for model "
                                   "%s. Review costs to ensure data is correct "
                                   "(--graph_models or --print_models)", k)

                    break
                if slope == 0:
                    break

    def get_required_components(self):
        """
        Should be overridden by derived class. Return a list of required
        components (capitalized) used in the layout.

        Example: return ['ATM', 'LND', 'ICE']
        """
        return []

    def check_requirements(self):
        """
        Check to make sure that each element of the subclass's list of
        required components has some data provided.
        """
        for r in self.get_required_components():
            if r not in self.models:
                logger.critical("Data for component %s not available", r)

    def write_timings(self, fd=sys.stdout, level=logging.DEBUG):
        """
        Print out the data used for the ntasks/cost models.
        Can be used to check that the data provided to the
        model is reasonable. Also see graph_costs()
        """
        assert self.state != self.STATE_UNDEFINED,\
               "set_data() must be called before write_timings()"
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
        assert self.state != self.STATE_UNDEFINED,\
               "set_data() must be called before graph_costs()"
        try:
            import matplotlib.pyplot as pyplot
        except ImportError:
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
        Must set self.state using LpStatus object:
              LpStatusOptimal    -> STATE_SOLVED_OK
              LpStatusNotSolved  -> STATE_UNSOLVED
              LpStatusInfeasible -> STATE_SOLVED_BAD
              LpStatusUnbounded  -> STATE_SOLVED_BAD
              LpStatusUndefined  -> STATE_UNDEFINED
              -- use self.set_state(lpstatus) --
        Returns state

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
        Return a dictionary of the solution variables, can be overridden.
        Default implementation returns values in self.X
        """
        assert self.state == self.STATE_SOLVED_OK,\
               "solver failed, no solution available"
        retval = {}
        if hasattr(self,'X') and isinstance(self.X, dict):
            for k in self.X:
                retval[k] = self.X[k].varValue
        return retval

    def set_state(self, lpstatus):
        if lpstatus == pulp.constants.LpStatusOptimal:
            self.state = self.STATE_SOLVED_OK
        elif lpstatus == pulp.constants.LpStatusNotSolved:
            self.state = self.STATE_UNSOLVED
        elif lpstatus == pulp.constants.LpStatusUndefined:
            self.state = self.STATE_UNDEFINED
        else:
            self.state = self.STATE_SOLVED_BAD

    def get_state(self):
        return self.state

    def get_state_string(self, state):
        return self.states[state]

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
        logger.info("Writing pe node info to %s", pefilename)
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
        if xmllint is not None:
            run_cmd("%s --format --output %s -" % (xmllint, pefilename),
                    input_str=ET.tostring(root))
