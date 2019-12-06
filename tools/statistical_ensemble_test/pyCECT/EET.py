#!/usr/bin/python

import re
import argparse
import itertools

class exhaustive_test(object):

    def __init__(self):
        super(exhaustive_test, self).__init__()

    def file_to_sets(self, compfile):
        set_dict = {}
        with open(compfile, 'r') as f:
            for line in f:
                line.strip
                key, failset = line.replace(' ', '').split(';', 1)

                try:
                    failset = map(int, failset.split(','))
                    failset = set(failset)

                except:
                    failset = set()

                set_dict[key] = failset

        return set_dict

    def test_combinations(self, dictionary, runsPerTest=3, nRunFails=2):
        sims = dictionary.keys()

        passed = failed = 0
        for compset in itertools.combinations(sims, runsPerTest):
            # This block is slightly slower than manually 
            # specifying the pairs, but it generalizes
            # easily.
            failsets = [dictionary[s] for s in compset]
            # The following three lines are adapted from 
            # user doug's answer in
            # http://stackoverflow.com/questions/27369373/pairwise-set-intersection-in-python
            pairs = itertools.combinations(failsets, 2)
            isect = lambda a, b: a.intersection(b)
            isect_list = [isect(t[0], t[1]) for t in pairs]
            isect_tot = set()
            [isect_tot.update(x) for x in isect_list]

            if len(isect_tot) > nRunFails:
                # print statements for debugging
                # print "this set failed"
                # print compset
                failed += 1
            else:
                # print "this set passed"
                # print compset
                passed +=1

        return passed, failed

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to calculate all combinations of ensemble tests")
    parser.add_argument("-f", dest="compfile",
    help="compfile location", metavar="PATH")

    args = parser.parse_args()

    eet = exhaustive_test()
    compare_dict = eet.file_to_sets(args.compfile)
    print "failure percent is %s" % eet.test_combinations(compare_dict)
