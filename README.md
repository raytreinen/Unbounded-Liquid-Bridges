# Unbounded-Liquid-Bridges

Codes associated with my paper "Discussion of a uniqueness result in "Equilibrium Configurations for a Floating Drop""

GPL-3.0-or-later

As of June 28th, these codes support the numerical portions of the paper.  Updates may happen in the future to make repository more general-use.

Each of these programs needs Chebfun. Please install Chebfun from chebfun.org.

# unbd_plotter

This script file plots a generating curve for an unbounded liquid bridge.  One can set the radius sigma of the vertial point (sigma,T(sigma)) along the curve.  One may also set the upper angle phi to a value in the range (0,pi/2).  The current range of sigma values allowed is in [0.085, 2].  The upper limit can be changed, but changing the lower limit is likely to lead to poor performance.

This file is slow to run as it computes the curve (sigma,T(sigma)) in order to build the unbounded liquid bridge.  If this code is used frequently, one can run the first line once and save the data in memory or a file, then the remaining lines can be run using that data.

# symmetric_capillary_unbounded

This function file is called by unbd_plotter and generates the data needed, including computing (sigma,T(sigma)) over an interval.

# symmetric_capillary_unbounded_paper

This function file is an alteration of the previous file that also generates the figures in the paper.
