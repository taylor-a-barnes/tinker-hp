c
c     Sorbonne University
c     Washington University in Saint Louis
c     University of Texas at Austin
c
c     ####################################################################
c     ##                                                                ##
c     ##  module shunt  --  polynomial switching function coefficients  ##
c     ##                                                                ##
c     ####################################################################
c
c
c     off    distance at which the potential energy goes to zero
c     off2   square of distance at which the potential goes to zero
c     cut    distance at which switching of the potential begins
c     cut2   square of distance at which the switching begins
c     c0     zeroth order coefficient of multiplicative switch
c     c1     first order coefficient of multiplicative switch
c     c2     second order coefficient of multiplicative switch
c     c3     third order coefficient of multiplicative switch
c     c4     fourth order coefficient of multiplicative switch
c     c5     fifth order coefficient of multiplicative switch
c     f0     zeroth order coefficient of additive switch function
c     f1     first order coefficient of additive switch function
c     f2     second order coefficient of additive switch function
c     f3     third order coefficient of additive switch function
c     f4     fourth order coefficient of additive switch function
c     f5     fifth order coefficient of additive switch function
c     f6     sixth order coefficient of additive switch function
c     f7     seventh order coefficient of additive switch function
c
c
#include "tinker_macro.h"
      module shunt
      implicit none
      real(t_p) off,off2,cut,cut2
      real(t_p) c0,c1,c2,c3,c4,c5
      real(t_p) f0,f1,f2,f3,f4,f5,f6,f7
      save
      end
