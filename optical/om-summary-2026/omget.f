      PROGRAM omget

      USE om_retrieve 
    
      integer nnn
      real*8  REnergies (500)      

      OPEN(unit=5,file='ominput.inp')

      READ(5,*) Number_Energies
      
      IF (Number_Energies.gt.0) 
     >     READ(5,*) (REnergies(n),n=1,Number_Energies)
      Def_Index = 1      
      READ(5,*) Ztarget,Atarget,RIPL_Index,Calc_Type
      IF (iabs(Calc_Type).eq.3) READ(5,*) Def_Index

      do nnn=1,Number_Energies
        Energies(nnn) = REnergies(nnn)
	enddo

      CALL retrieve
       
      CLOSE(5)
                   
      END PROGRAM omget
