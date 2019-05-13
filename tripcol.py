# Python program for triple collocation
#
# Version 1.0  01-05-2019
# Jur Vogelzang (KNMI)
# with help from Jos de Kloe and Jeroen Verspeek (KNMI)
#

import argparse
from triple_collocation_module import do_tc

# routine for printing the available command line arguments
#
def usage():
    print("tc:")
    print("tc:  usage: python tripcol.py <-i IN> [OPTIONS]")
    print("tc:  with < > mandatory arguments and [ ] free options")
    print("tc:")
    print("tc:  mandatory arguments:")
    print("tc:  -i  <IN>          : read collocations from ASCII file IN")
    print("tc:  --input <IN>      : same as -i")
    print("tc:")
    print("tc:  free arguments:")
    print("tc:  -f <F>            : set sigma test factor to F (default: 4.0)")
    print("tc:  --f_sigma <F>     : same as -f")
    print("tc:")
    print("tc:  -m <M>            : set maximum number of iterations to M (default: 20)")
    print("tc:  --maxiter <M>     : same as -m")
    print("tc:")
    print("tc:  -p <EPS           : set precision to EPS (default: 0.00001)")
    print("tc:  --precision <EPS> : same as -p")
    print("tc:")
    print("tc:  -r <R2>           : set representativeness error variance to R2 (default: 0.0)")
    print("tc:  --reprerr <R2>    : same as -r")
    print("tc:")


# MAIN PROGRAM
# define settings
#
infile  = ""
fs      = float(4.0)
maxiter = int(20)
re      = float(0.0)
prec    = float(0.00001)
verb    = int(1)

print("tc:")
print("tc:  program tripcol.py - Python program for triple collocation")

# define command line arguments
#
parser = argparse.ArgumentParser()
parser.add_argument("-i" , "--input"     , help="Read collocations from file INPUT (mandatory)")
parser.add_argument("-f" , "--f_sigma"   , help="Set sigma test factor to F_SIGMA (default: 4.0)")
parser.add_argument("-m" , "--maxiter"   , help="Set maximum number of iterations to MAXITER (default: 20)")
parser.add_argument("-r" , "--reprerr"   , help="Set representativeness error to REPRERR (default: 0.0)")
parser.add_argument("-p" , "--precision" , help="Set precision to PRECISION (default: 0.00001)")
parser.add_argument("-v" , "--verbosity" , help="The higher, the more output (default: 0)")

# read command line arguments
#
args = parser.parse_args()
if not args.input:
    print("tc:")
    print("tc:  ERROR: no file with collocations given")
    usage()
    print("tc:")
    print("tc:  program tripcol.py aborted")
    quit()
infile = args.input

if args.f_sigma:
    fs = float(args.f_sigma)
if args.maxiter:
    maxiter = int(args.maxiter)
if args.reprerr:
    re = float(args.reprerr)
if args.precision:
    prec = float(args.precision)
if args.verbosity:
    verb = int(args.verbosity)

# do triple collocation
#
result = do_tc(infile, f_sigma=fs, max_nr_of_iterations=maxiter, repr_err=re, precision=prec, verbosity=verb)
