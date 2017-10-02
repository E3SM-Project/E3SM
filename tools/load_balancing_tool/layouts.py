import optimize_model
import pulp

class IceLndAtmOcn(optimize_model.OptimizeModel):
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

    Returns state (STATE_SOLVED_OK, STATE_SOLVED_BAD, STATE_UNSOLVED)
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
        set solution in self.X
        set state STATE_SOLVED_OK if solved,
        otherwise STATE_SOLVED_BAD
        """
        assert self.state != self.STATE_UNDEFINED,\
               "set_data() must be called before optimize()!"
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
        self.set_state(self.prob.status)
        return self.state

    def get_solution(self):
        """
        Return a dictionary of the solution variables.
        """
        assert self.state == self.STATE_SOLVED_OK,\
               "solver failed, no solution available"
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
        assert self.state == self.STATE_SOLVED_OK,\
               "solver failed, no solution available"
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

class IceLndWavAtmOcn(optimize_model.OptimizeModel):
    """
    Optimized the problem based on the Layout
              __________________________
             | ICE  |  LND  | WAV |     |
             |______|_______|_____|     |
             |                    | OCN |
             |    ATM             |     |
             |____________________|_____|

      Min T
      s.t.  T[ice]      <= T1
            T[lnd]      <= T1
            T[wav]      <= T1
            T1 + T[atm] <= T
            T[ocn]      <= T

            NB[c]        >= 1 for c in [ice,lnd,wav,ocn,atm]
            NB[ice] + NB[lnd] + NB[wav] <= NB[atm]
            atm_blocksize*NB[atm] + ocn_blocksize*NB[ocn] <= TotalTasks
            (NB[*] is number of processor blocks)

            T[c]        >= C[c]_{i} - NB[c]_{i} *
                       (C[c]_{i+1} - C[c]_{i}) / (NB[c]_{i+1} - NB[c]_{i})
                       + NB[c] * (C[c]_{i+1} - C[c]_{i})
                                               / (NB[c]_{i+1} - NB[c]_{i}),
                        i=1..ord(NB), c in [ice,lnd,wav,ocn,atm]

    These assumptions are checked when solver is initialized
      . Assuming cost is monotonic decreasing vs ntasks
      . Assuming perfect scalability for ntasks < tasks[0]
      . Assuming same scalability factor for ntasks > ntasks[last] as for
                              last two data points
      . Assuming components are capable of running on ntasks

    Returns state (STATE_SOLVED_OK, STATE_SOLVED_BAD, STATE_UNSOLVED)
    If solved, then solution will be stored in self.X dictionary, indexed
    by variable name. Suggested convention:
        'Tice', 'Tlnd', ... for cost per component
        'Nice', 'Nlnd', ... for ntasks per component
        'NBice', 'NBlnd', ... for number of blocks per component
    """

    def __init__(self):
        self.models = {}


    def get_required_components(self):
        return ['LND', 'ICE', 'WAV', 'ATM', 'OCN']

    def optimize(self):
        """
        Run the optimization.
        set solution in self.X
        set state STATE_SOLVED_OK if solved,
        otherwise STATE_SOLVED_BAD
        """
        assert self.state != self.STATE_UNDEFINED,\
               "set_data() must be called before optimize()!"
        self.atm = self.models['ATM']
        self.lnd = self.models['LND']
        self.ice = self.models['ICE']
        self.ocn = self.models['OCN']
        self.wav = self.models['WAV']
        self.real_variables = ['TotalTime', 'T1', 'Tice', 'Tlnd', 'Tatm',
                               'Tocn', 'Twav']
        self.integer_variables = ['NBice', 'NBlnd', 'NBatm', 'NBocn', 'NBwav',
                                  'Nice', 'Nlnd', 'Natm', 'Nocn', 'Nwav']
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
        self.constraints.append([X['Twav'] - X['T1'] <= 0, "Twav - T1 == 0"])
        self.constraints.append([X['T1'] + X['Tatm'] - X['TotalTime'] <= 0,
                                 "T1 + Tatm - TotalTime <= 0"])
        self.constraints.append([X['Tocn'] - X['TotalTime'] <= 0,
                                 "Tocn - TotalTime == 0"])
        self.constraints.append([X['Nice'] + X['Nlnd'] + X['Nwav'] - X['Natm'] == 0,
                                 "Nice + Nlnd + Nwav - Natm == 0"])
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
        self.constraints.append([self.wav.blocksize * X['NBwav'] - X['Nwav'] == 0,
                                 "Nwav = %d * NBwav" % self.wav.blocksize])

        # These are the constraints based on the timing data.
        # They should be the same no matter what the layout of the components.
        self.add_model_constraints()

        for c, s in self.constraints:
            self.prob += c, s

        # Write the program to file and solve (using coin-cbc)
        self.prob.writeLP("IceLndWavAtmOcn_model.lp")
        self.prob.solve()
        self.set_state(self.prob.status)
        return self.state

    def get_solution(self):
        """
        Return a dictionary of the solution variables.
        """
        assert self.state == self.STATE_SOLVED_OK,\
               "solver failed, no solution available"
        return {'NBLOCKS_ICE':self.X['NBice'].varValue,
                'NBLOCKS_LND':self.X['NBlnd'].varValue,
                'NBLOCKS_WAV':self.X['NBwav'].varValue,
                'NBLOCKS_ATM':self.X['NBatm'].varValue,
                'NBLOCKS_OCN':self.X['NBocn'].varValue,
                'NTASKS_ICE':self.X['Nice'].varValue,
                'NTASKS_LND':self.X['Nlnd'].varValue,
                'NTASKS_WAV':self.X['Nwav'].varValue,
                'NTASKS_ATM':self.X['Natm'].varValue,
                'NTASKS_OCN':self.X['Nocn'].varValue,
                'NTASKS_TOTAL':self.maxtasks,
                'COST_ICE':self.X['Tice'].varValue,
                'COST_LND':self.X['Tlnd'].varValue,
                'COST_WAV':self.X['Twav'].varValue,
                'COST_ATM':self.X['Tatm'].varValue,
                'COST_OCN':self.X['Tocn'].varValue,
                'COST_TOTAL':self.X['TotalTime'].varValue}

    def write_pe_file(self, pefilename):
        """
        Write out a pe_file that can be used to implement the
        optimized layout
        """
        assert self.state == self.STATE_SOLVED_OK,\
               "solver failed, no solution available"
        natm = int(self.X['Natm'].varValue)
        nlnd = int(self.X['Nlnd'].varValue)
        nice = int(self.X['Nice'].varValue)
        nocn = int(self.X['Nocn'].varValue)
        nwav = int(self.X['Nwav'].varValue)

        ntasks = {'atm':natm, 'lnd':nlnd, 'rof':1, 'ice':nice,
                  'ocn':nocn, 'glc':1, 'wav':nwav, 'cpl':1}
        roots = {'atm':0, 'lnd':0, 'rof':0, 'ice':nlnd,
                 'ocn':natm, 'glc':0, 'wav':nlnd+nice, 'cpl':0}
        nthrds = {}
        for c in ['atm', 'lnd', 'rof', 'ice', 'ocn', 'glc', 'wav', 'cpl']:
            if c.upper() in self.models:
                nthrds[c] = self.models[c.upper()].nthrds
            else:
                nthrds[c] = 1
        self.write_pe_template(pefilename, ntasks, nthrds, roots)
