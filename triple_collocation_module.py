# Python program for triple collocation
#
# Version 2.0  27-07-2024
# Jur Vogelzang (KNMI, Netherlands)
# with help from Jos de Kloe and Jeroen Verspeek (KNMI, Netherlands)
# Special thanks to Weicheng Ni (NUDT, China) for correct handling of
# representativeness errors and Federico Cossu (CSIC, Spain) for
# improving the sigma test (implemented here as a distance squared test)
#

from __future__ import print_function
import math

class TripColProcess:
    """ do triple collocation """

    def __init__(self):
        self.b         = [0.0, 0.0, 0.0]                                       # calibration biases
        self.a         = [1.0, 1.0, 1.0]                                       # calibration scalings
        self.M1        = [0.0, 0.0, 0.0]                                       # first-order moments (averages)
        self.M2        = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]   # second-order moments
        self.C         = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]   # covariances
        self.db        = [0.0, 0.0, 0.0]                                       # increment in calibration bias (additive)
        self.da        = [1.0, 1.0, 1.0]                                       # increment in calibration scaling (multiplicative)
        self.d2max     = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]   # maximum allowed distance squared between each pair
        self.t2        = 0.0                                                   # common variance of all three systems
        self.errvar    = [0.0, 0.0, 0.0]                                       # error variances
        self.accepted  = 0                                                     # number of accepted collocations
        self.rejected  = 0                                                     # number of rejected collocations (failed sigma test)
        self.converged = False                                                 # has triple collocation converged?


    def update_distances(self, input_file, f_sigma):
    #
    #   reads the collocations and calculates the distances squared multiplied by
    #   f_sigma^2 for the sigma test
    #
        self.d2max = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        f = open(input_file , "r")
        n = 0

        # read collocations from file
        #
        for line in f:
            x = line.split( )

            # calibrate the collocation data with the current calibration coefficients
            # note: the error model is x = a*t + b + delta, with delta the error
            # here we need an estimate for t (the calibrated value),
            # so we must apply the inverse formula t = x_cal = (x - b)/a
            #
            xcal = []
            for i in range(3):
                xcal.append((float(x[i]) - self.b[i])/self.a[i])

            # update distances
            #            
            n = n + 1
            R = 1.0/float(n)
            for i in range(3):
                for j in range(i+1,3):
                    dist2 = (xcal[i] - xcal[j])**2
                    self.d2max[i][j] = self.d2max[i][j] + R*(dist2 - self.d2max[i][j])
        f.close()

        # multiply with f_sigma^2 and symmetrize
        #
        for i in range(3):
            for j in range(i+1,3):
                self.d2max[i][j] = f_sigma*f_sigma * self.d2max[i][j]
                self.d2max[j][i] = self.d2max[i][j] 

    
    def update_moments(self, input_file):
    #
    #   reads the collocations, applies the sigma-test, and updates the moments
    #
        self.M1 = [0.0, 0.0, 0.0]
        self.M2 = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        f = open(input_file , "r")
        n = 0
        m = 0

        # read collocations from file
        #
        for line in f:
            x = line.split( )

            # calibrate the collocation data with the current calibration coefficients
            # note: the error model is x = a*t + b + delta, with delta the error
            # here we need an estimate for t (the calibrated value),
            # so we must apply the inverse formula t = x_cal = (x - b)/a
            #
            xcal = []
            for i in range(3):
                xcal.append((float(x[i]) - self.b[i])/self.a[i])

            # sigma test (implemented here as a distance-squared test)
            #
            accept = True
            for i in range(3):
                for j in range(i+1,3):
                    dist2 = (xcal[i] - xcal[j])**2
                    if dist2 > self.d2max[i][j]:
                        accept = False
                        
            # update moments if the collocation passed the sigma test
            #
            if accept:
                n = n + 1
                R = 1.0/float(n)

                for i in range(3):
                    self.M1[i] = self.M1[i] + R*(xcal[i] - self.M1[i])
                    for j in range(3):
                        self.M2[i][j] = self.M2[i][j] + R*(xcal[i]*xcal[j] - self.M2[i][j])
            else:
                m = m + 1
        f.close()

        # set and check number of accepted collocations
        #
        self.accepted  = n
        self.rejected = m

        if n < 2:
            print("tc:")
            print("tc:  ERROR in routine update_moments: number of accepted collocations equals",n)
            print("tc:  insufficient number of collocations passed the sigma test")
            print("tc:  triple collocation aborted")
            print("tc:")
            quit()

        
    def print_moments(self, verbosity):
        if verbosity > 4:
            print("tc:")
            print("tc:  - first moments (averages)")
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.M1[0], self.M1[1], self.M1[2]))
            print("tc:")
            print("tc:  - second moments")
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.M2[0][0], self.M2[0][1], self.M2[0][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.M2[1][0], self.M2[1][1], self.M2[1][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.M2[2][0], self.M2[2][1], self.M2[2][2]))
        if verbosity > 5: 
            print("tc:")
            print("tc:  - maximum distances squared for sigma test")
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.d2max[0][0], self.d2max[0][1], self.d2max[0][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.d2max[1][0], self.d2max[1][1], self.d2max[1][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.d2max[2][0], self.d2max[2][1], self.d2max[2][2]))


    def solve(self, repr_err, precision):
        #
        # calculate covariances
        #
        for i in range(3):
            for j in range(3):
                self.C[i][j] = self.M2[i][j] - self.M1[i]*self.M1[j]

        self.C[0][0] = self.C[0][0] - repr_err
        self.C[0][1] = self.C[0][1] - repr_err
        self.C[1][0] = self.C[1][0] - repr_err
        self.C[1][1] = self.C[1][1] - repr_err

        # solve covariance equations
        # the first moments and the covariances are from calibrated data
        #
        self.t2 = self.C[1][0] * self.C[2][0] / self.C[2][1]                     # common variance

        self.da[1] = self.C[2][1] / self.C[2][0]                                 # calibration scaling increment system 1
        self.da[2] = self.C[2][1] / self.C[1][0]                                 # calibration scaling increment system 2

        self.db[1] = self.M1[1] - self.da[1]*self.M1[0]                          # calibration bias increment system 1
        self.db[2] = self.M1[2] - self.da[2]*self.M1[0]                          # calibration bias increment system 2

        self.errvar[0] = self.C[0][0] - self.C[1][0]*self.C[2][0]/self.C[2][1]   # error variance system 0
        self.errvar[1] = self.C[1][1] - self.C[1][0]*self.C[2][1]/self.C[2][0]   # error variance system 1
        self.errvar[2] = self.C[2][2] - self.C[2][0]*self.C[2][1]/self.C[1][0]   # error variance system 2

        # update calibration scalings and biases; check if calculation has converged
        #
        self.converged = True

        for i in range(1,3):
            self.a[i] = self.a[i]*self.da[i]
            if abs(self.da[i] - 1.0) > precision:
                self.converged = False

            self.b[i] = self.b[i] + self.db[i]
            if abs(self.db[i]) > precision:
                self.converged = False


    def print_results(self, verbosity):
        if verbosity > 2:
            print("tc:")
            print("tc:  - covariances")
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.C[0][0], self.C[0][1], self.C[0][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.C[1][0], self.C[1][1], self.C[1][2]))
            print("tc:    " , "{:12.6f}{:12.6f}{:12.6f}".format(self.C[2][0], self.C[2][1], self.C[2][2]))
            print("tc:")
            print("tc:  - increments of calibration parameters")
            print("tc:                     system 0    system 1    system 2")
            print("tc:    -------------------------------------------------")
            print("tc:    - scalings :" , "{:12.6f}{:12.6f}{:12.6f}".format(self.da[0], self.da[1], self.da[2]))
            print("tc:    - biases   :" , "{:12.6f}{:12.6f}{:12.6f}".format(self.db[0], self.db[1], self.db[2]))
        if verbosity > 3:
            print("tc:")
            print("tc:  - intermediate triple collocation results:")
            print("tc:                                system 0    system 1    system 2")
            print("tc:    ------------------------------------------------------------")
            print("tc:    - calibration scalings:" , "{:12.6f}{:12.6f}{:12.6f}".format(self.a[0], self.a[1], self.a[2]))
            print("tc:    - calibration biases  :" , "{:12.6f}{:12.6f}{:12.6f}".format(self.b[0], self.b[1], self.b[2]))
            print("tc:    - error variances     :" , "{:12.6f}{:12.6f}{:12.6f}".format(self.errvar[0], self.errvar[1], self.errvar[2]))
            print("tc:    - common variance     :" , "{:12.6f}".format(self.t2))


def print_settings(input_file, f_sigma, max_nr_of_iterations, repr_err, precision, verbosity):
    #
    # print current settings and parameters
    #
    print("tc:")
    print("tc:  settings for triple collocation")
    print("tc:  - input collocation file            :" , input_file)
    print("tc:  - sigma test factor                 :" , "{:12.6f}".format(f_sigma))
    print("tc:  - maximum number of iterations      :" , "{:12d}".format(max_nr_of_iterations))
    print("tc:  - precision                         :" , "{:12.6f}".format(precision))
    print("tc:  - representativeness error variance :" , "{:12.6f}".format(repr_err))
    print("tc:  - verbosity level                   :" , "{:12d}".format(verbosity))


def do_tc(input_file, f_sigma=4.0, max_nr_of_iterations=20, repr_err=0.0, precision=0.00001, verbosity=1):
    #
    # do triple collocation
    #
    if verbosity > 0:
        print_settings(input_file, f_sigma, max_nr_of_iterations, repr_err, precision, verbosity)

    # start triple collocation calculation; start iteration
    #
    tc = TripColProcess()

    for i in range(max_nr_of_iterations):
        iteration = i + 1
        if verbosity > 1:
            print("tc:")
            print("tc:  iteration",iteration)
            print("tc:  - accepted collocations       " , "{:12d}".format(tc.accepted))
            print("tc:  - rejected collocations       " , "{:12d}".format(tc.rejected))
            print("tc:  - total number of aollocations" , "{:12d}".format(tc.accepted + tc.rejected))

        # update distances for sigma test; update moments
        #
        tc.update_distances(input_file, f_sigma)

        tc.update_moments(input_file)
        tc.print_moments(verbosity)

        # solve covariance equations
        #
        tc.solve(repr_err, precision)
        tc.print_results(verbosity)

        # print output if needed when iteration has converged and stop iteration loop
        #
        if tc.converged:
            stddev = [0.0 , 0.0 , 0.0]
            for j in range(3):
                if tc.errvar[j] > 0.0:
                    stddev[j] = math.sqrt(tc.errvar[j])

            if verbosity > 0:
                print("tc:")
                print("tc:  triple collocation converged at iteration",iteration)
                print("tc:  final results, calibration in the form of t = (x - b)/a")
                print("tc:                                      system 0    system 1    system 2")
                print("tc:  --------------------------------------------------------------------")
                print("tc:  - calibration scalings a      :" , "{:12.6f}{:12.6f}{:12.6f}".format(tc.a[0], tc.a[1], tc.a[2]))
                print("tc:  - calibration biases b        :" , "{:12.6f}{:12.6f}{:12.6f}".format(tc.b[0], tc.b[1], tc.b[2]))
                print("tc:  - error variances             :" , "{:12.6f}{:12.6f}{:12.6f}".format(tc.errvar[0], tc.errvar[1], tc.errvar[2]))
                print("tc:  - error standard deviations   :" , "{:12.6f}{:12.6f}{:12.6f}".format(stddev[0], stddev[1], stddev[2]))
                print("tc:")
                print("tc:  - common variance             :" , "{:12.6f}".format(tc.t2))
                print("tc:  - accepted collocations       :" , "{:12d}".format(tc.accepted))
                print("tc:  - rejected collocations       :" , "{:12d}".format(tc.rejected))
                print("tc:  - total number of collocations:" , "{:12d}".format(tc.accepted + tc.rejected))
                print("tc:")
                print("tc:  triple collocation completed succesfully")
                print("tc:")
            break
    
    # iteration loop has ended; print warning if it has not converged
    #
    if not tc.converged:
        print("tc:")
        print("tc:  WARNING: triple collocation did not converge")
        print("tc:")

    # define result and return it
    #
    result = []
    result.append(tc.a)
    result.append(tc.b)
    result.append(tc.errvar)
    result.append(tc.t2)
    result.append(tc.accepted)
    result.append(tc.rejected)

    return result

