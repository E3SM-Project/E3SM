import sys, logging
import pulp
import optimize_model
logger = logging.getLogger(__name__)

class AtmLnd(optimize_model.OptimizeModel):
    def get_required_components(self):
        return ['ATM', 'LND', 'ROF', 'ICE', 'CPL', 'OCN']

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
        """
        assert self.state != self.STATE_UNDEFINED,\
            "set_data() must be called before optimize()!"
        self.atm = self.models['ATM']
        self.lnd = self.models['LND']
        self.ice = self.models['ICE']
        self.ocn = self.models['OCN']
        self.rof = self.models['ROF']
        self.cpl = self.models['CPL']

        self.real_variables = ['TotalTime', 'Tice', 'Tlnd', 'Tatm',
                               'Tocn', 'Trof', 'Tcpl']
        self.integer_variables = ['NBice', 'NBlnd', 'NBatm', 'NBocn',
                                  'NBrof', 'NBcpl', 'Nrof', 'Ncpl',
                                  'Nice', 'Nlnd', 'Natm', 'Nocn', 'N1']
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
        self.constraints.append([X['Tatm'] + X['Trof'] + X['Tcpl'] - X['TotalTime'] <= 0, "Tatm + Trof + Tcpl - TotalTime <= 0"])
        self.constraints.append([X['Tlnd'] + X['Tice'] + X['Tocn'] - X['TotalTime'] <= 0, "Tlnd + Tice + Tocn - TotalTime <= 0"])

        self.constraints.append([X['Natm'] - X['N1'] == 0,
                            "Natm - N1 <= 0"])
        self.constraints.append([X['Nrof'] - X['N1'] == 0,
                            "Nrof - N1 <= 0"])
        self.constraints.append([X['Ncpl'] - X['N1'] == 0,
                            "Ncpl - N1 <= 0"])

        self.constraints.append([X['Nlnd'] + X['N1'] == self.maxtasks,
                                 "Nlnd + N1 <= MAXN"])
        self.constraints.append([X['Nice'] + X['N1'] == self.maxtasks,
                                 "Nice + N1 <= MAXN"])
        self.constraints.append([X['Nocn'] + X['N1'] == self.maxtasks,
                                 "Nocn + N1 <= MAXN"])

        self.constraints.append([self.atm.blocksize * X['NBatm'] - X['Natm'] == 0,
                            "Natm = %d * NBatm" % self.atm.blocksize])
        self.constraints.append([self.ice.blocksize * X['NBice'] - X['Nice'] == 0,
                            "Nice = %d * NBice" % self.ice.blocksize])
        self.constraints.append([self.lnd.blocksize * X['NBlnd'] - X['Nlnd'] == 0,
                            "Nlnd = %d * NBlnd" % self.lnd.blocksize])
        self.constraints.append([self.ocn.blocksize * X['NBocn'] - X['Nocn'] == 0,
                            "Nocn = %d * NBocn" % self.ocn.blocksize])
        self.constraints.append([self.rof.blocksize * X['NBrof'] - X['Nrof'] == 0,
                            "Nrof = %d * NBrof" % self.rof.blocksize])
        self.constraints.append([self.cpl.blocksize * X['NBcpl'] - X['Ncpl'] == 0,
                            "Ncpl = %d * NBcpl" % self.cpl.blocksize])

        # These are the constraints based on the timing data.
        # They should be the same no matter what the layout of the components.
        self.add_model_constraints()

        for c, s in self.constraints:
            self.prob += c, s

        # Write the program to file and solve (using glpk)
        self.prob.writeLP("IceLndAtmOcn_model.lp")
        self.prob.solve()
        self.set_state(self.prob.status)
        return self.state
    

    def write_pe_file(self, pefilename):
        """
        Write out a pe_file that can be used to implement the 
        optimized layout
        """
        natm = self.X['Natm'].varValue
        nlnd = self.X['Nlnd'].varValue
        nice = self.X['Nice'].varValue
        nocn = self.X['Nocn'].varValue
        ncpl = self.X['Ncpl'].varValue
        nrof = self.X['Nrof'].varValue
        npart = max(natm, nrof, ncpl)
        ntasks = {'atm':natm, 'lnd':nldn, 'rof':nrof, 'ice':nice,
                  'ocn':nocn, 'glc':1, 'wav':1, 'cpl':ncpl}
        roots = {'atm':0, 'lnd':npart, 'rof':0, 'ice':npart,
                'ocn':npart, 'glc':0, 'wav':0, 'cpl':0}
        nthrds = {}
        for c in ['atm', 'lnd', 'rof', 'ice', 'ocn', 'glc', 'wav', 'cpl']:
            nthrds[c] = self.models[c.upper()].nthrds
        
        self.write_pe_template(pefilename, ntasks, nthrds, roots)
