C------------------------------------------------------------
C
C Rothamsted Carbon and Nitrogen Turnover Model
C Example WRAPPER
C
C Adapted from SUNDIAL (MAGEC)
C by Jo Smith & Kevin Coleman
C 02/03/05
C
C Modified for highly organic soils
C by Jo Smith, Bente Foereid & Matt Aitkenhead
C Started 01/04/05
C
C Modified for tropical soils
C by Jo Smith & Pia Gottschalk
C Started 01/04/05
C
C Modified for spatial applications
C by Jo Smith 
C Started 01/08/06
C
C------------------------------------------------------------
	IMPLICIT NONE

C
C Variables passed to/from calling subroutine
C Mode variables
C
	INTEGER MODE				! IN:Mode of model run (site or spatial)
	INTEGER GISMODE				! IN:Spatial simulation of cells
	INTEGER ISWAIT				! OUT:Code to wait for key press (1) or not (0)
	INTEGER SITEMODE			! IN:Site specific simulation
	INTEGER LIMMODE				! IN:Limited data site specific simulation
	INTEGER TESTMODE			! IN:Test mode for VSD
C
C Wait for key stroke on noerror?
C
      ISWAIT=1
C
C Title on screen
C
	PRINT*,'*********************************************************'
	PRINT*,'                                                         '
	PRINT*,'                                                         '
	PRINT*,'  E____________C____O_______S______S___________E         '
	PRINT*,'  Estimator of C in Organic Soils: Sequestrn & Emissions '
	PRINT*,'                                                         '
	PRINT*,'               _______________________                   '
	PRINT*,'                                                         '
	PRINT*,'                ECOSSE VERSION 6.1.a1                    '
	PRINT*,'                     ___________                         '
	PRINT*,'                                                         '
	PRINT*,'                Modular SUNDIAL-MAGEC                    '
	PRINT*,'                      + 5cm layers for all soil states   '
	PRINT*,'*********************************************************'
C
C Get mode of model run
C
      CALL GETMODE(MODE,SITEMODE,GISMODE,LIMMODE)
C
C-----------------------------------------------------------------------
C Site specific simulation
C-----------------------------------------------------------------------
	IF(MODE.EQ.SITEMODE)THEN
	  CALL ECOSSE_SITE_RUN(ISWAIT) 
C-----------------------------------------------------------------------
C Spatial simulation of cells
C-----------------------------------------------------------------------
	ELSEIF(MODE.EQ.GISMODE)THEN
	  CALL ECOSSE_GIS_RUN(ISWAIT)
C-----------------------------------------------------------------------
C Site specific simulation using limited data
C-----------------------------------------------------------------------
	ELSEIF(MODE.EQ.LIMMODE)THEN
	  CALL ECOSSE_LIM_RUN(ISWAIT)
C-----------------------------------------------------------------------
C Test run of VSD
C-----------------------------------------------------------------------
c	ELSEIF(MODE.EQ.4)THEN
c	  CALL TEST_ECOSSE_FOR_VSD(ISWAIT)
	ENDIF
401   END
C
C--------------------------------------------------------
C
	SUBROUTINE GETMODE(MODE,SITEMODE,GISMODE,LIMMODE)
C
C Variables passed to / from this routine
C
	INTEGER MODE				! IN:Mode of model run (site or spatial)
	INTEGER GISMODE				! IN:Spatial simulation of cells
	INTEGER SITEMODE			! IN:Site specific simulation
	INTEGER LIMMODE				! IN:Limited data site specific simulation
C
C Set value of SITEMODE and GISMODE
C
	SITEMODE=1
	GISMODE=2
	LIMMODE=3
	TESTMODE=4
C
C Ask user for mode of simulation
C
101   CONTINUE
      WRITE(*,10)
10	FORMAT('Choose mode of model run'/
     &       ' 1 = Site specific'/
     &       ' 2 = Spatial simulation of cells'/
     &       ' 3 = Limited data site simulation'/
     &       ' 4 = Test run of VSD')
      MODE=2
	READ(*,*)MODE
C
C Check entered number is OK
C
	IF(MODE.NE.SITEMODE.AND.MODE.NE.GISMODE.AND.MODE.NE.LIMMODE.AND.
     &   MODE.NE.TESTMODE)THEN
	  PRINT*,'Incorrect mode entered'
	  GOTO 101
      ENDIF
C
C Echo mode back to user
C 
      IF(MODE.EQ.SITEMODE)THEN
	  PRINT*,'*******************************************************'
	  PRINT*,' SITE SPECIFIC SIMULATION'
	  PRINT*,'*******************************************************'
      ELSEIF(MODE.EQ.GISMODE)THEN
	  PRINT*,'*******************************************************'
	  PRINT*,' SPATIAL SIMULATION OF CELLS'
	  PRINT*,'*******************************************************'
      ELSEIF(MODE.EQ.LIMMODE)THEN
	  PRINT*,'*******************************************************'
	  PRINT*,' SITE SPECIFIC SIMULATION USING LIMITED DATA'
	  PRINT*,'*******************************************************'
      ELSEIF(MODE.EQ.TESTMODE)THEN
	  PRINT*,'*******************************************************'
	  PRINT*,' LIMITED DATA RUN - TEST FOR VSD IMPLEMENTATION'
	  PRINT*,'*******************************************************'
	ENDIF
      END
C
C--------------------------------------------------------
C
      SUBROUTINE TELL_MODEL(DMODEL,ICMODEL,INMODEL,CNMODEL,DOCMODEL,
     &                      SPARMODEL,EQMODEL,IMFUNC,ITFUNC,CH4MODEL,
     &                      ISPINUP,ISWAIT)

C
C Subroutine to tell the user about the model
C      
      IMPLICIT NONE
C
C Variables local to this subroutine
C
	INTEGER CH4_OFF				! CH4 model off
	INTEGER CH4_RICHARDS	    ! Richards CH4 model on		
	INTEGER CH4_AITKENHEAD      ! Aitkenhead CH4 model on
	INTEGER CNFOEREID			! C:N ratio obtained by method of Foereid
	INTEGER CNMAGEC				! C:N ratio obtained by method of MAGEC
	INTEGER DBRADBURY			! Bradbury model for denitrification
	INTEGER DNEMIS				! Nemis model for denitrification
	INTEGER DOC_OFF				! DOC model off
	INTEGER DOC_ON				! DOC model on
	INTEGER EQNPP				! Model initialised using plant inputs 
								!    from measured NPP
	INTEGER EQTOC				! Model initialised using measured TOC
	INTEGER EQNPPTOC			! Model initialised using both plant inputs
								!    from measured NPP and measured TOC
	INTEGER EQHILLIER			! Model initialised using TOC and
								!    the Hillier solver to get C pools and plant inputs
	INTEGER ICFIXED				! Initialisation of C is fixed
	INTEGER ICROTHCEQ			! Initialisation of C by RothC 
								!  equilibrium run
	DATA ICFIXED,ICROTHCEQ /1,2/
	INTEGER ISPINUP_OFF			! N limitation spin-up is not used
	INTEGER ISPINUP_ON			! N limitation spin-up is used
	DATA ISPINUP_OFF,ISPINUP_ON /0,1/
      INTEGER INSTABLECN			! Initisalisation of N by assuming 
								!  steady state (after Bradbury)
	INTEGER ISWAIT				! IN:Code to wait for key press (1) or not (0)
	INTEGER INPASSCN			! Initialisation of N by passing the
								!  C:N ratio of DPM and RPM
	DATA INSTABLECN,INPASSCN /1,2/
	INTEGER IMFUNC_ROTHC		! ROTHC moisture rate modifier
	INTEGER IMFUNC_HADLEY		! HADLEY moisture rate modifier
	INTEGER ITFUNC_ROTHC		! ROTHC temperature rate modifier
	INTEGER ITFUNC_HADLEY		! HADLEY temperature rate modifier
	INTEGER SPARFILE			! Soil parameters read in from file
	INTEGER SPARCALC			! Soil parameters calculated from TOC etc

	DATA CH4_OFF,CH4_RICHARDS,CH4_AITKENHEAD /0,1,2/
	DATA CNMAGEC,CNFOEREID /1,2/
	DATA DBRADBURY,DNEMIS /1,2/
	DATA DOC_ON,DOC_OFF /1,2/
	DATA EQNPP,EQTOC,EQNPPTOC,EQHILLIER /1,2,3,4/
	DATA ICFIXED,ICROTHCEQ /1,2/
      DATA IMFUNC_ROTHC,IMFUNC_HADLEY /0,1/
	DATA INSTABLECN,INPASSCN /1,2/
      DATA ITFUNC_ROTHC,ITFUNC_HADLEY /0,1/
	DATA SPARFILE,SPARCALC /1,2/

C
C Variables passed to/from calling subroutine
C
	INTEGER CH4MODEL			! IN:Methane model (off, Richards or Aitkenhead model on) 
	INTEGER CNMODEL				! IN:Type of Model for C:N ratio of DPM and 
	                            !     RPM (MAGEC or Foereid)
	INTEGER DMODEL				! IN:Type denitrification model chosen 
	                            ! (Bradbury or NEMIS)
	INTEGER DOCMODEL			! IN:DOC model (on or off)
	INTEGER EQMODEL				! IN:Type of equilibrium run 
	                            ! (NPP, TOC or both)
	INTEGER IMFUNC				! IN:Choice of moisture rate modifier 
								! (ROTHC or HADLEY)
	INTEGER ICMODEL				! IN:Type of C initialisation model 
      INTEGER ISPINUP				! IN:Is N limitation spin-up used?
	INTEGER INMODEL				! IN:Type of N initialisation model
	INTEGER ITFUNC				! IN:Choice of temperature rate modifier 
								! (ROTHC or HADLEY)
	INTEGER SPARMODEL			! IN:Soil parameter model (from file or calc)
C
C Tell the user about the model
C
	IF(DMODEL.EQ.DNEMIS)THEN
      PRINT*,'                      + NEMIS denitrification			 '
	PRINT*,'                      + partitioning of gaseous emissions'
	ELSE
      PRINT*,'                      + Bradbury denitrification	     '
	ENDIF
	IF(ICMODEL.EQ.ICROTHCEQ)THEN
	PRINT*,'                      + C initialisation by ROTHC eq.run'
	ELSEIF(ICMODEL.EQ.ICFIXED)THEN
	PRINT*,'                      + C initialisation fixed'
	ENDIF 
	IF(ISPINUP.EQ.ISPINUP_ON)THEN
	PRINT*,'                      + spin up for N limitation'
	ENDIF 
	IF(INMODEL.EQ.INSTABLECN)THEN
	PRINT*,'                      + N init.assuming steady state'
	ELSEIF(ICMODEL.EQ.INPASSCN)THEN
	PRINT*,'                      + N init.by passed C:N ratio'
	ENDIF 
	IF(CNMODEL.EQ.CNFOEREID)THEN
	PRINT*,'                      + Foereid C:N model                ' 
	ELSE
	PRINT*,'                      + MAGEC C:N model                  '
	ENDIF 
      IF(DOCMODEL.EQ.DOC_ON)THEN
	PRINT*,'                      + Aitkenhead DOC model             '
	ENDIF
      IF(SPARMODEL.EQ.SPARFILE)THEN
	PRINT*,'                      + soil parameters from file        '
	ELSE
	PRINT*,'                      + soil parameters from TOC & %clay '
	ENDIF
      IF(EQMODEL.EQ.EQNPP)THEN
	PRINT*,'                      + initialised to measured NPP'
      ELSEIF(EQMODEL.EQ.EQTOC)THEN
	PRINT*,'                      + initialised to measured TOC'
	ELSEIF(EQMODEL.EQ.EQNPPTOC)THEN
	PRINT*,'                      + initialised to measured NPP & TOC'
	ELSEIF(EQMODEL.EQ.EQHILLIER)THEN
	PRINT*,'                      + initialised by Hillier solver'
	ENDIF
      IF(IMFUNC.EQ.IMFUNC_ROTHC)THEN
	PRINT*,'                      + ROTHC moisture modifiers         '
	ELSE
	PRINT*,'                      + Hadley Centre moisture modifiers '
	ENDIF
      IF(ITFUNC.EQ.ITFUNC_ROTHC)THEN
	PRINT*,'                      + ROTHC temperature modifiers      '
	ELSE
	PRINT*,'                      + Hadley Centre temp.modifiers     '
	ENDIF
      IF(CH4MODEL.EQ.CH4_RICHARDS)THEN
	PRINT*,'                      + Richards methane emissions model '
	ELSEIF(CH4MODEL.EQ.CH4_AITKENHEAD)THEN
	PRINT*,'                      + Aitkenhead methane emissions model'
	ENDIF
	PRINT*,'                                                         '
	PRINT*,'          Last modified on 25/11/09 by Jo Smith          '
	PRINT*,'*********************************************************'
	PRINT*,''
      IF(ISWAIT.EQ.1)THEN
	PRINT*,'                             ...Press any key to continue'
	READ(*,*)
	ENDIF
C
C Leave TELL_MODEL
C
      END
