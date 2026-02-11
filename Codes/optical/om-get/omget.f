      PROGRAM omget
     
      USE om_retrieve 
     
      integer n
      
      OPEN(unit=5,file='ominput.inp')

      READ(5,*) Number_Energies
      
      IF (Number_Energies.gt.0) 
     >     READ(5,*) (Energies(n),n=1,Number_Energies)
      Def_Index = 1      
      READ(5,*) Ztarget,Atarget,RIPL_Index,Calc_Type
      IF (iabs(Calc_Type).eq.3) READ(5,*) Def_Index

      CALL retrieve
       
      CLOSE(5)
                   
      END PROGRAM omget
