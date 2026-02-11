c
      double precision e1,el,ein,en
c
c  E1  = current  CM energy
c  EL  = current lab energy
c  EIN = array of CM energies
c  EN  = working array for reading of energies (lab or CM)
c
      common /ener  /  e1,el,ein(nemax),en(nemax)
c-----------------------------------------------------------------------
