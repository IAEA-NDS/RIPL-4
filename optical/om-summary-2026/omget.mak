# Microsoft Developer Studio Generated NMAKE File, Format Version 4.20
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=omget - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to omget - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "omget - Win32 Release" && "$(CFG)" != "omget - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "omget.mak" CFG="omget - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "omget - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "omget - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "omget - Win32 Debug"
F90=fl32.exe
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "omget - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\omget.exe"

CLEAN : 
	-@erase "$(INTDIR)\om_retrieve.obj"
	-@erase "$(INTDIR)\omget.obj"
	-@erase "$(OUTDIR)\omget.exe"

# ADD BASE F90 /Ox /c /nologo
# ADD F90 /Ox /c /nologo
F90_PROJ=/Ox /c /nologo 
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/omget.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/omget.pdb" /machine:I386 /out:"$(OUTDIR)/omget.exe" 
LINK32_OBJS= \
	"$(INTDIR)\om_retrieve.obj" \
	"$(INTDIR)\omget.obj"

"$(OUTDIR)\omget.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "omget - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\omget.exe"

CLEAN : 
	-@erase "$(INTDIR)\om_retrieve.obj"
	-@erase "$(INTDIR)\omget.obj"
	-@erase "$(OUTDIR)\omget.exe"
	-@erase "$(OUTDIR)\omget.ilk"
	-@erase "$(OUTDIR)\omget.pdb"

# ADD BASE F90 /Zi /c /nologo
# ADD F90 /Zi /c /nologo
F90_PROJ=/Zi /c /nologo /Fd"omget.pdb" 
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/omget.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/omget.pdb" /debug /machine:I386 /out:"$(OUTDIR)/omget.exe" 
LINK32_OBJS= \
	"$(INTDIR)\om_retrieve.obj" \
	"$(INTDIR)\omget.obj"

"$(OUTDIR)\omget.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for.obj:
   $(F90) $(F90_PROJ) $<  

.f.obj:
   $(F90) $(F90_PROJ) $<  

.f90.obj:
   $(F90) $(F90_PROJ) $<  

CPP_PROJ=/ML 

.c.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx.obj:
   $(CPP) $(CPP_PROJ) $<  

.c.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx.sbr:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Target

# Name "omget - Win32 Release"
# Name "omget - Win32 Debug"

!IF  "$(CFG)" == "omget - Win32 Release"

!ELSEIF  "$(CFG)" == "omget - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\omget.f
DEP_F90_OMGET=\
	".\om_retrieve.mod"\
	

"$(INTDIR)\omget.obj" : $(SOURCE) $(DEP_F90_OMGET) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\om_retrieve.f

!IF  "$(CFG)" == "omget - Win32 Release"

F90_MODOUT=\
	"om_retrieve"


"$(INTDIR)\om_retrieve.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ELSEIF  "$(CFG)" == "omget - Win32 Debug"

F90_MODOUT=\
	"om_retrieve"


"$(INTDIR)\om_retrieve.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


!ENDIF 

# End Source File
# End Target
# End Project
################################################################################
