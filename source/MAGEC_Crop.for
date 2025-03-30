C------------------------------------------------------------
C
C Rothamsted Carbon and Nitrogen Turnover Model
C Soil Crop Routines
C
C Adapted from SUNDIAL (MAGEC)
C by Jo Smith & Kevin Coleman
C 02/03/05
C
C-------------------------------------------------------------
C
C EXTERNAL SUBROUTINES
C 1. RUN1_MAGEC_CROP
C 2. RUN2_MAGEC_CROP
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN1_MAGEC_CROP(I_TS,IEND,IYEAR,IK,SECONDS,
     &                  CACTOT,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,
     &                  ISOWN,IANTHES,NSOW,ICOVER,CULTIVATE,NCULT)
C
C Subroutine to Calculate the total C returns for next crop and set up crop at sowing
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
C
C Variables passed to/from calling subroutine
C
	INTEGER MAXCROP				! Max.no.of crop types allows
      PARAMETER (MAXCROP=36)
      INTEGER CULTIVATE			! IN(CALL): Code to cultivate this week (1=yes, 0=no)
	INTEGER NCULT				! IN(CALL): Number of cultivations
	INTEGER I_TS,ISTART,IEND,IYEAR,IK,ICROP,N_STEPS,ISOWN,IANTHES,NSOW
	REAL SECONDS
	REAL CINP,ANINP,CACTOT,CREQN,TRNIN,RNIN,RNIN15,CATOT15
      REAL SXORGN,SXORGN15,C_TS,SXORGC,TCINP,XORGC,XORGN,ORGC,ORGN
	REAL CRATE(0:MAXCROP),ANRATE(0:MAXCROP),MCROP,CONVER_F
	INTEGER ICOVER
C
C Set times for MAGEC crop
C
      CALL SETTIME_MAGEC_CROP(SECONDS,CONVER_F,N_STEPS)
C
C Count number of weeks since sowing
C
      CALL COUNTW(IYEAR,IK,ISTART,N_STEPS,ISOWN,IANTHES,I_TS,NSOW)
C
C Add any crop debris
C 
      CALL RESIDU(I_TS,ISTART,IEND,IYEAR,CINP,ANINP,IK,
     &                  CACTOT,CREQN,TRNIN,RNIN,RNIN15,CATOT15,
     &                  ICROP,SXORGN,SXORGN15,C_TS,SXORGC,TCINP,
     &                  XORGC,XORGN,ORGC,ORGN,CRATE,ANRATE,ICROP,
     &                  CONVER_F,CULTIVATE,NCULT)
C
C Work out if there is crop cover or not
C
      IF(IYEAR.GT.0.AND.IK.GT.NSOW+ANINT(4/CONVER_F).AND.IK.LE.IEND)THEN
	  ICOVER=1
      ELSE
        ICOVER=0
      ENDIF
C
C Leave RUN1_MAGEC_CROP
C
	END
C
C-------------------------------------------------------------
C
      SUBROUTINE RUN2_MAGEC_CROP(CRIT,IYEAR,ICROP,IK,NSOW,MEND,C_TS,
     &                          RNIN,DTR,RDD,DVP,VP,TMAX,
     &                          TMMX,TMIN,TMMN,DAVTMP,NAVTMP,
     &                          DTRJM2,TAVS,BBRADS,SVPS,SLOPES,
     &                          RLWNS,NRADS,PENMRS,WDF,WN,PENMD,
     &                          PE,AE,SOILW,EVAPW,DDAYS,IDAG,
     &                          WMAX,N_TS,RAIN,
     &                          LAT,SOILN,AMMN,CRITMIN,IROCKJ,INSTRAW,
     &                          EVAPO,CACT,SEEDN_S,SEED,MCROP,NSOIL)
C
C Subroutine to run crop Calculate the total C returns for next crop
C
      IMPLICIT NONE
C
C Variables local to this subroutine
C
      REAL DELT, BOLTZM, LHVAP, VHCA, PSYCH
	REAL CO2A, COEFN,COEFR,COEFT,COEFV	
	PARAMETER (DELT = 1.)
      PARAMETER (BOLTZM=5.668E-8,LHVAP=2.4E6,VHCA=1200.,PSYCH=0.067)
      PARAMETER (CO2A=350.,COEFN=1.,COEFR=1.,COEFT=0.,COEFV=1.)
      INTEGER MAXSOIL,MAXCROP,MAXLAYER	
	PARAMETER (MAXSOIL=50,MAXCROP=36,MAXLAYER=60)
C
C Variables passed to/from calling subroutine
C
      INTEGER IYEAR, ICROP, IK, NSOW, MEND, IDAG, INSTRAW, MCROP 
	INTEGER NSOIL, N_TS, IROCKJ

C REAL passed from calling subroutine
      REAL CRIT(MAXSOIL,MAXLAYER)
	REAL SOILW(MAXLAYER), SOILN(MAXLAYER), AMMN(MAXLAYER)
	REAL SEED(0:MAXCROP)
      REAL C_TS, RNIN, DTR, RDD, DVP, VP
      REAL TMAX, TMMX, TMIN, TMMN, DAVTMP, NAVTMP, DTRJM2, TAVS
      REAL BBRADS, SVPS, SLOPES, RLWNS, NRADS, PENMRS
      REAL WDF, WN, PENMD, PE, AE, EVAPW, DDAYS
      REAL WMAX(12),  RAIN, LAT
      REAL CRITMIN,  EVAPO, CACT, SEEDN_S

C
      critmin = crit(nsoil,1)
      IF(IYEAR.EQ.0.OR.ICROP.EQ.0.OR.(IYEAR.GT.0.AND.IK.LT.NSOW)
     &	    .or.(IYEAR.GT.0.AND.IK.GT.MEND))THEN
        c_ts = 0.0
	  rnin = 0.0

        DTR    = (RDD/1.E+6) * COEFR
        DVP    = VP          * COEFV
        TMAX   = TMMX + COEFT
        TMIN   = TMMN + COEFT
        DAVTMP = 0.29*TMIN + 0.71*TMAX
        NAVTMP = 0.71*TMIN + 0.29*TMAX
        DTRJM2 = DTR * 1.E6
        TAVS   = (DAVTMP+NAVTMP)/2.
        BBRADS = BOLTZM * (TAVS+273.)**4 * 86400.
        SVPS   = 0.611 * EXP(17.4 * TAVS / (TAVS + 239.))
        SLOPES = 4158.6 * SVPS / (TAVS + 239.)**2
        RLWNS  = BBRADS * MAX(0., 0.55*(1.-DVP/SVPS))
        NRADS  = DTRJM2 * (1.-0.15) - RLWNS
        PENMRS = NRADS * SLOPES/(SLOPES+PSYCH)
      
        WDF    = 2.63*(1.0 + 0.54 * WN)
        PENMD  = LHVAP*WDF*MAX(0.,(SVPS-DVP))*PSYCH/(SLOPES+PSYCH)

        PE     = MAX(0., (PENMRS + PENMD)/LHVAP)
        AE     = MIN (PE, soilw(1))
            
        EVAPW  = AE
        ddays = 0.0
        IDAG = 0 
	else
	  ddays = ddays + (TMMN + TMMX)/2
ckwc: 16/3/00            if(ik.eq.nsow)write(76,*) n_ts, ' sowing date '
	  if(ddays.lt.45.0) then
c         write(76,*)n_ts,' sowing - germ '
	  else 
c	    write(76,*)n_ts, ' germ - harvest ' 
		call cropmodel(DELT, BOLTZM, CO2A, COEFN, COEFR, COEFT,  
     &              COEFV,LHVAP, PSYCH, VHCA,wmax,
     &              ik, nsow, mend, n_ts, rain, rdd, tmmn, tmmx, 
     &              vp, wn, lat, soiln,ammn, critmin, irockj, icrop, 
     &              instraw, soilw,c_ts, rnin, evapo, cact, SEEDN_S,
     &              IDAG, iyear) 
          SEED(MCROP) = SEEDN_S 
	    EVAPW = evapo
	  endif
	endif
C
C Leave RUN_MAGEC_CROP
C
	END
C
C-------------------------------------------------------------
C
C INTERNAL ROUTINES
C
C-------------------------------------------------------------
C
C Plant model
C written by Xinyou Yin 27 Sept 1999
C Linked to Sundial October 1999 
C
      subroutine cropmodel(DELT, BOLTZM, CO2A, COEFN, COEFR, COEFT, 
     &  COEFV, LHVAP, PSYCH, VHCA, maxwhc,
     &  IK, NSOW, MEND, DOY, RAIN, RDD, TMMN, TMMX, 
     &  VP, WN, LAT, soilno3, soilnh4, critmin, iRtRest, cropID, 
     &  instraw, soilh2o, c_return, n_return, evapo,crop_off, SEEDN_S, 
     &  IDAG, iyear)
     
      IMPLICIT NONE 

C INTEGER passed from calling subroutine
      INTEGER IK, NSOW, MEND, DOY, iRtRest, cropID, IDAG, iyear, instraw

C REAL passed from calling subroutine
      REAL soilno3(12), soilnh4(3), soilh2o(12), maxwhc(12)
	REAL DELT, BOLTZM, CO2A, COEFN, COEFR, COEFT
      REAL COEFV, LHVAP, PSYCH, VHCA
      REAL RAIN, RDD, TMMN, TMMX
      REAL VP, WN, LAT,  critmin  
      REAL c_return, n_return, evapo, crop_off, SEEDN_S

C INTEGER local to subroutine
      INTEGER nlayers, I

C REAL local to subroutine
 
C*     State variables, initial values and rates
      REAL DS, ZERO, DVR
      REAL WLV, WLVI, RWLV
      REAL WLVD, LWLV
      REAL WST, RWST
      REAL WSO, RWSO
      REAL WRT, WRTI, RWRT
      REAL WRTD, LWRT
      REAL WSTR, LSTR
      REAL WLVDS, LVDS
      REAL NTOT, NTOTI, RNTOT
      REAL NSH, NSHI, RNSH
      REAL NRT, NRTI, RNRT
      REAL NST, RNST
      REAL NLV, NLVI, RNLV
      REAL NSO, RNSO
      REAL SLNB, SLNBI, RSLNB
      REAL LAI, LAII, RLAI
      REAL TSOIL, TSOILI, RTSOIL
      REAL TSUM, HU

C*     Crop parameters 
      REAL DSR, EG, FSTR, FVPD, HNCMIN, HTMAX, K, KN, KW, LWMAX, NPL 
      REAL PCC, RDMAX, RNCMIN, SEEDNC, SEEDW, SLNMIN, TBD, TOD, LNC0
      REAL TSUMR, TSUMV, WSUP, GAMMA, C3C4
 

C*     Other calculated variables
      REAL ACO2I, ADIF, AE, APCAN, APCANN, ARCICA, ARSW, ASSA, AT
      REAL ATLEAF, BBRADS, CO2I, DAVTMP, DAYL, DEC, DERI, DLAI
      REAL DTR, DTRJM2, DVP, DVRR, DVRV, FCO2I, FDIF, FGC, FNRADC, FNSH
      REAL FPCAN, FPT, FRCICA, FRSW, FWLV, FWSH, FWSO, FWST, GBCROP
      REAL GBLEAF, GC, HNC, HT, INCN, LAIN, LITC, LITN, LNC, LNLV, LNRT
      REAL LWIDTH, NAVTMP, NDEM, NNC, NONC, NRADC, NRADS, NSUP, NUPT
      REAL ONC, PCAN, PCANA, PCANNA, PE, PENMD, PENMRS, PI, PNC, PT
      REAL RAD, RBH, RBW, RCICA, RD, RGR, RGRRT, RGRSH, RLWNS, RNC
      REAL RRTSH, RSW, RT, RTSA, SDIF, SHSA, SHSAN, SLA, SLN, SLNBC
      REAL SLNNT, SLNT, SLOPE, SLOPEL, SLOPES, SVP, SVPA, SVPL, SVPS
      REAL TAVS, TAVSD, TAVST, TLEAF, TMAX, TMIN, VPD, VPDA, VPDL, WDF
      REAL WND, WSH, WTOT, Y, BETA, PT1

      REAL MaxNO3(12), MaxNH4(3)
	REAL NdemNO3, NdemNH4
	REAL maxdem

	REAL TotAvNO3, TotAvNH4

c FRWN: fraction that N supply is reduced when water is low
c it only makes a big difference when the water content is low. 
	REAL frwn(12), h2o50, maxwhc50, FRWN50

	REAL c_root, c_st_lv, c_annual, n_root, n_st_lv, n_annual

      character*15 cropname(30)
	REAL cropnum(30), ARfvpd(30), ARseedw(30),
     &     ARseednc(30), AReg(30), ARpcc(30), ARtbd(30), ARtod(30),       
     &     ARlnc0(30), ARtsumv(30), ARtsumr(30), ARhtmax(30),  
     &     ARlwmax(30),ARrdmax(30),ARnpl(30), ARfstr(30), ARdsr(30),
     &     ARhncmin(30), ARrncmin(30),ARslnmin(30), ARk(30), ARkn(30), 
     &     ARkw(30),ARgamma(30), ARc3c4(30)
 

c UptNO3 & UptNH4 are temporary variables  
	REAL UptNO3, UptNH4
  
c      save

      ZERO = 1e-10


c      if(ik.eq.nsow)then  
       if (IDAG.eq.0) then
ckwc: 16/3/00	   write(76,*) DOY, ' Germination '
	   IDAG = 1
		    
	   if(instraw.eq.0)then 
	      beta = 0.3
	   else 
	      beta = 1.0
	   endif

	   open(41,file='crop_par.dat',status='unknown')
         do 10 i=1,30
	     read(41,*)cropnum(i), cropname(i)
		   read(41,*)ARfvpd(i), ARseedw(i), ARseednc(i), 
     &               AReg(i), ARpcc(i), ARtbd(i), ARtod(i),ARlnc0(i)
           read(41,*)ARtsumv(i), ARtsumr(i), ARhtmax(i), ARlwmax(i), 
     &               ARrdmax(i), ARnpl(i), ARfstr(i), ARdsr(i)
           read(41,*)ARhncmin(i), ARrncmin(i), ARslnmin(i), 
     &               ARk(i), ARkn(i), ARkw(i), ARgamma(i), ARc3c4(i)
 10      continue
         close(41)
	   FVPD    =  ARfvpd(cropID)
	   SEEDW   =  ARseedw(cropID)
	   SEEDNC  =  ARseednc(cropID)
         EG      =  AReg(cropID)
	   PCC     =  ARpcc(cropID)
	   TBD     =  ARtbd(cropID)
	   TOD     =  ARtod(cropID)
	   LNC0    =  ARlnc0(cropID)
	   TSUMV   =  ARtsumv(cropID)
	   TSUMR   =  ARtsumr(cropID)
	   HTMAX   =  ARhtmax(cropID)
	   LWMAX   =  ARlwmax(cropID)
         RDMAX   =  ARrdmax(cropID)
	   NPL     =  ARnpl(cropID)
	   FSTR    =  ARfstr(cropID)
	   DSR     =  ARdsr(cropID)
         HNCMIN  =  ARhncmin(cropID)
	   RNCMIN  =  ARrncmin(cropID)
	   SLNMIN  =  ARslnmin(cropID)
         K       =  ARk(cropID)
	   KN      =  ARkn(cropID)
	   KW      =  ARkw(cropID)
	   GAMMA   =  ARgamma(cropID)
	   C3C4    =  ARc3c4(cropID)

         WLVI   = NPL * SEEDW * EG * 0.5
         WRTI   = NPL * SEEDW * EG * 0.5
 
         NTOTI  = NPL * SEEDW * EG * 0.5 * (0.00415 + 1.514*LNC0)

	   SEEDN_S=(SEEDNC -0.5*EG*(0.00415+1.514*LNC0))*NPL*SEEDW
         NSHI   =(NTOTI*WLVI-0.00405*WLVI*WRTI)/(WLVI+0.514*WRTI)
         NRTI   = NTOTI - NSHI
         NLVI   = NSHI
 
         LAII   = WLVI*(HNCMIN/SLNMIN)
         SLNBI  = NLVI/LAII
 
         TSOILI = 10.
 
C*        Initialize state variables
         DS    = ZERO
         WLV   = WLVI
         WLVD  = ZERO
         WST   = ZERO
         WSO   = ZERO
         WRT   = WRTI
         WRTD  = ZERO
         WSTR  = ZERO
         WLVDS = ZERO
         NTOT  = NTOTI
         NSH   = NSHI
         NRT   = NRTI
         NST   = ZERO
         NLV   = NLVI
         NSO   = ZERO
         SLNB  = SLNBI
         LAI   = LAII
         TSOIL = TSOILI
         TSUM  = ZERO

	   maxdem = ZERO
 
C*        Dynamic calculations

c         TIME = STTIME
c         FNSH = NSHI/NTOTI
         SLN = SLNB  
	 
      endif
c      IF (FNSH.LT.0.995) THEN
       IF ((DS.LT.2.0).and.(SLN.GT.SLNMIN)) THEN
c      IF (SLN.GT.SLNMIN) THEN
         if(DS.LE.1.0)then
ckwc: 16/3/00	     write(76,*) DOY, ' Germ - grain fill '
	   else
ckwc: 16/3/00	  	 write(76,*) DOY, ' grain fill - maturity '
         endif

         TMAX   = TMMX + COEFT
         TMIN   = TMMN + COEFT
         DAVTMP = 0.29*TMIN + 0.71*TMAX
         NAVTMP = 0.71*TMIN + 0.29*TMAX
         DTR    = (RDD/1.E+6) * COEFR
         DVP    = VP          * COEFV
         WND    = 1.33* MAX (0.2, WN)
 
         PI     = 3.1416
 
         DVRV   = 1./TSUMV
         DVRR   = 1./TSUMR
 
         WSH    = WLV + WST + WSO
         Y      = 12./44.* 0.6/PCC

         FWLV   = 1. - DS
         IF (FWLV.LT.0.) THEN
           FWLV = 0.
         ENDIF
         
         FWSO   = 1./(1.-DSR)-1./(1.-DSR)*DS
         IF (FWSO.LT.0.) THEN
           FWSO = 0.
         ENDIF
         IF (FWSO.GT.1.) THEN
           FWSO = 1.
         ENDIF
 
         LNC    = NLV / WLV
         ONC    = NSO / MAX(1E-8, WSO)
         RNC    = NRT / WRT
 
         SLN    = NLV/LAI
         SLNT   = NLV*KN             /(1.-EXP(-KN*LAI))
         SLNBC  = NLV*KN*EXP(-KN*LAI)/(1.-EXP(-KN*LAI))
 
         LAIN   = LOG(1.+KN*NLV/SLNMIN)/KN
 
         SLA    = LAI/WLV
         IF (SLA.LT.0.02) THEN
            SLA = 0.02
         ENDIF
         IF (SLA.GT.0.05) THEN
            SLA = 0.05
         ENDIF
 
         WDF    = 2.63*(1.0 + 0.54 * WN)
 
         HT     = HTMAX/(1.+298.*EXP(-0.009*TSUM))
         LWIDTH = LWMAX/(1.+     EXP(-0.005*TSUM))

C stop the roots growing if you have rooting restrictions 
   	   if(iRtRest.eq.1)then
	     RD     = MIN(50., RDMAX/LOG(0.05)*LOG(MAX(0.05, 1.-DS)))
	   elseif (iRtRest.eq.2)then
           RD     = MIN(100., RDMAX/LOG(0.05)*LOG(MAX(0.05, 1.-DS)))
	   else
	     RD     = RDMAX/LOG(0.05)*LOG(MAX(0.05, 1.-DS))
	   endif

	   nlayers = int(rd/5) + 1
 
         LWLV   = (LAI - MIN(LAI,LAIN)) /SLA
 
         FWST   = 1. - FWLV - FWSO
         RAD    = PI/180.
         WTOT   = WSH + WRT
         HNC    = NSH / WSH
 
         RRTSH  = WRT / WSH
         RSLNB  = SLNBC - SLNB
         DLAI   = (WLVD-WLVDS)*SLA
 
         RT     = 0.74*(LOG((2.-0.7*HT)/(0.1*HT)))**2/(0.4**2*WND)
 
         GBLEAF = 0.01*SQRT(WND/LWIDTH)
         SVP    = 0.611 * EXP(17.4 * DAVTMP / (DAVTMP + 239.))
 
         DTRJM2 = DTR * 1.E6
         TAVS   = (TSOIL + (DAVTMP+NAVTMP)/2.) / 2.
         LNLV   = LWLV*SLNMIN*SLA
         DEC    = -ASIN(SIN(23.45*RAD)*COS(2.*PI*(DOY+10.)/365.))
         PNC    = NTOT/ WTOT
         LWRT   = (max(0., -0.02+ 0.035*DS)) * WRT
         GBCROP = (1.-EXP(-0.5*KW*LAI))/(0.5*KW)*GBLEAF
         VPD    = MAX (0., SVP  - DVP)
         SLOPE  = 4158.6 * SVP /(DAVTMP + 239.)**2
 
         BBRADS = BOLTZM * (TAVS+273.)**4 * 86400.
         SVPS   = 0.611 * EXP(17.4 * TAVS / (TAVS + 239.))
         INCN   = 0.001*PNC
         LNRT   = LWRT*RNCMIN
         RBH    = 1./GBCROP
 
         FRCICA = 0.9 - FVPD * VPD
         DAYL   = 0.5*(1.+2.*ASIN(TAN(RAD*LAT)*TAN(DEC))/PI)
         SLOPES = 4158.6 * SVPS / (TAVS + 239.)**2
         RLWNS  = BBRADS * MAX(0.,0.55*(1.-DVP/SVPS))
         SLNNT  = (NLV+NLV/NTOT*INCN*WTOT)*KN /(1.-EXP(-KN*LAI))
         RBW    = 0.93 * RBH
         FCO2I  = FRCICA * CO2A
         NRADS  = DTRJM2 * (1.-0.15) - RLWNS
         PENMD  = LHVAP * WDF * MAX(0.,(SVPS-DVP))*PSYCH/(SLOPES+PSYCH)
 
         CALL PCANF(C3C4,FCO2I,SLNMIN,SLNT,K,LAI,DAVTMP,DAYL,DTR,FPCAN)
         PENMRS = NRADS * SLOPES/(SLOPES+PSYCH)
 
         FGC    = 0.9*FPCAN*(273.+DAVTMP)/0.536341/86400./DAYL/(CO2A-
     $    FCO2I)
 
         PE     = MAX(0., EXP(-0.5*(LAI+DLAI))*(PENMRS + PENMD)/LHVAP)
         FRSW   = MAX(1E-10, 1./FGC - RBW*1.3 - RT)/1.6
 
         CALL 
     $    PTRAN(FRSW,RT,RBW,RBH,BOLTZM,LHVAP,VHCA,PSYCH,LAI,DTRJM2,DAVTM
     $    P,DAYL,DVP,SLOPE,SVP,VPD,FPT,FNRADC)
         CALL DIFLA (LHVAP,VHCA,FNRADC,FPT,DAYL,RBH,RT,FDIF)
 
         TLEAF  = DAVTMP + FDIF
         SVPL   = 0.611 * EXP(17.4 * TLEAF / (TLEAF + 239.))
         SLOPEL = (SVPL - SVP) /(TLEAF - DAVTMP)
 
         VPDL   = MAX  (0., SVPL - DVP)
         RCICA  = 0.9 - FVPD * VPDL
         CO2I   = RCICA * CO2A
 
         CALL PCANF(C3C4,CO2I,SLNMIN,SLNT,K,LAI,TLEAF,DAYL,DTR,PCAN)
 
         GC     = 0.9*PCAN *(273.+TLEAF)/0.536341/86400./DAYL/(CO2A-
     $    CO2I)
 
         RSW    = MAX(1E-10, 1./GC - RBW*1.3 - RT)/1.6
 
         CALL PTRAN 
     $    (RSW,RT,RBW,RBH,BOLTZM,LHVAP,VHCA,PSYCH,LAI,DTRJM2,TLEAF,DAYL,
     $    DVP,SLOPEL,SVPL,VPD,PT,NRADC)
         CALL DIFLA (LHVAP,VHCA,NRADC,PT,DAYL,RBH,RT,SDIF)

c ADD WATER IN EACH LAYER HERE 
         WSUP = 0.0
	   h2o50 = 0.0
         maxwhc50 = 0.0

	   soilh2o(1) = soilh2o(1) + 0.01
	   WSUP = WSUP + 0.01

	   do 25 i = 1,10
	     maxwhc50 = maxwhc50 + maxwhc(i)
	     h2o50 = h2o50 + soilh2o(i)
  25     continue

         if (nlayers.le.10)then
		 do 30 i = 1, nlayers
             WSUP = WSUP + soilh2o(i)  
  30       continue
	   elseif(nlayers.le.20)then
           do 40 i = 1,11
	       WSUP = WSUP + soilh2o(i)
  40       continue
	   else
	     do 50 i = 1,12
	       WSUP = WSUP + soilh2o(i)
  50       continue
	   endif 

C If Evaporation in all layers used the two statements below 
c	    AT     = MIN (PT, PT/(PT+PE)*WSUP)
c         AE     = MIN (PE, PE/(PT1+PE)*WSUP)

C If Evaporation only in the first layer used the three statements below
C AE proportion from layer 1 only
C AT proportion from layer 1 and all other rooting layers 

c         if(wsup.ne.0.0)then
	     PT1    = soilh2o(1)/WSUP*PT
c	   else
c	     PT1 = 0.0
c	   endif
         
c	   if(pt1+pe.ne.0.0)then
	     AT     = MIN (PT, (PT1/(PT1+PE)*soilh2o(1))+wsup-soilh2o(1))
           AE     = MIN (PE, PE/(PT1+PE)*soilh2o(1))
c	   else
c	     AT = 0.0
c	     AE = 0.0
c	   endif

	   evapo = AE + AT

	   write(64,6401)doy, evapo, AE,AT
 6401    format(1x, i6, 3f8.2)

         CALL DIFLA (LHVAP,VHCA,NRADC, AT,DAYL,RBH,RT,ADIF)
c         if (at.ne.0.0)then
		   ARSW   = (PT-AT)*(SLOPEL*(RBH+RT)+PSYCH*(RBW+RT))/
     $              AT/PSYCH+PT/AT*RSW
c	   else
c	     ARSW = RSW
c	   endif
C KWC NSUP is obtained from the soil : see below
c         NSUP   = COEFN*AT/PT*(0.8-0.39*DS)
         CALL SUBDD  (TMAX, TMIN, ADIF, DAYL, TBD ,TOD, HU )
         ATLEAF = DAVTMP + ADIF
 
         NONC =SEEDNC
c         NONC = MAX(0.55*SEEDNC, MIN(SEEDNC,0.75*LNC+0.25*SEEDNC))

         SVPA   = 0.611 * EXP(17.4 * ATLEAF / (ATLEAF + 239.))

         TAVSD  = (1.-EXP(-0.5*(LAI+DLAI)))*ATLEAF + EXP(-0.5*(LAI+DLAI)
     $    ) *DAVTMP
         LVDS   = (WLVD-WLVDS)/10.*HU/(TOD-TBD)
         CALL PHENOL (DS, nsow, ik, DVRV, DVRR, HU, DVR)
 
         VPDA   = MAX  (0., SVPA - DVP)
         TAVST  = (TAVSD + NAVTMP)/2.

         IF (DS.GT.DSR) THEN
           LSTR   = (WST+WSTR)*DVR*FSTR/(2.-DSR)
         ELSE
           LSTR   = 0.
         ENDIF
 
         LITC   = ((LWRT + LVDS)* PCC) * 10.
         LITN   =  (LNRT + LVDS * SLNMIN*SLA) * 10.

         ARCICA = 0.9 - FVPD * VPDA
         RTSOIL = (TAVST - TSOIL)/5.
         ACO2I  = ARCICA * CO2A
         CALL PCANF(C3C4,ACO2I,SLNMIN,SLNT, K,LAI,ATLEAF,DAYL,DTR,APCAN)
 
         CALL PCANF(C3C4,ACO2I,SLNMIN,SLNNT,K,LAI,ATLEAF,DAYL,DTR,
     &             APCANN)

         PCANA  = (1.6*RSW+1.3*RBW+RT)/(1.6*ARSW+1.3*RBW+RT)*APCAN
         PCANNA = (1.6*RSW+1.3*RBW+RT)/(1.6*ARSW+1.3*RBW+RT)*APCANN
 
         ASSA   = PCANA + 0.947*LSTR*PCC*44./12.
 
         RGR    = PCANA*Y/WTOT

         SHSA   = 12./44. * 0.6*PCANA  / WSH
         SHSAN  = 12./44. * 0.6*PCANNA / WSH
         DERI   = (SHSAN - SHSA)  / INCN 
         RTSA   = 1./PCC * SHSA**2/ DERI

c         NDEM   = MIN(RTSA * WRT, 0.5)
         NDEM   = MIN(MAX(RTSA * WRT, maxdem), 1.0)
c         NDEM   = RTSA * WRT

C To calculate Max Available in each layer
C NO3 0-5cm, ... 45-50cm (kg N / ha)
         do 60 i= 1, 10
	     FRWN(i) = 1. / (1.+9.*exp(-15.3*soilh2o(i)/maxwhc(i)))
c           FRWN(i) = 1.0
           MaxNO3(i)= FRWN(i)*(soilno3(i)-critmin/10.)
  60     continue
C NO3 50-100cm, 100-150cm (kg N / ha)
	   FRWN(11) = 1. / (1.+9.*exp(-15.3*soilh2o(11)/maxwhc(11)))
c         FRWN(11) = 1.0
         MaxNO3(11) = FRWN(11)*(soilno3(11) - critmin)
	   FRWN(12) = 1. / (1.+9.*exp(-15.3*soilh2o(12)/maxwhc(12)))
c         FRWN(12) = 1.0
	   MaxNO3(12) = FRWN(12)*(soilno3(12) - critmin)
C NH4 0-50cm, 50-100cm, 100-150cm (kg n / ha)
	   FRWN50 = 1. / (1.+9.*exp(-15.3*h2o50/maxwhc50))
c         FRWN50 = 1.0
	   MaxNH4(1) = FRWN50*(soilnh4(1) - critmin)
	   MaxNH4(2) = FRWN(11)*(soilnh4(2) - critmin)
	   MaxNH4(3) = FRWN(12)*(soilnh4(3) - critmin)

	   TotAvNO3 = 0.
	   TotAvNH4 = 0.

C Sum available NO3 (MaxNO3) & NH4 (MaxNH4) down to the 
C rooting depth (RD)
C NB : This is in kg N / ha
         if (nlayers.le.10)then
		 do 70 i = 1, nlayers
             TotAvNO3 = TotAvNO3 + MaxNO3(i)
  70       continue
           TotAvNH4 = MaxNH4(1)
	   elseif(nlayers.le.20)then
           do 80 i = 1,10
	       TotAvNO3 = TotAvNO3 + MaxNO3(i)
  80       continue
           TotAvNO3 = TotAvNO3 +  MaxNO3(11) 
      	 TotAvNH4 = MaxNH4(1) + MaxNH4(2)

	   else
	     do 90 i = 1,11
	       TotAvNO3 = TotAvNO3 + MaxNO3(i)
  90       continue
           TotAvNO3 = TotAvNO3 + MaxNO3(12)
	     TotAvNH4 = MaxNH4(1) + MaxNH4(2)+ MaxNH4(3)
	   endif   
         
C Sum Available NO3 & available NH4 to give NSUP (g N / m2)
	   if(rd.eq.0.)then
           NSUP   = COEFN*AT/PT*(0.8-0.39*DS)
	   else
           NSUP = (TotAvNO3 + TotAvNH4)
c          convert NSUP from kg N / ha to g N /m2
           NSUP = NSUP / 10.
	   endif
C
C check this I think it should be max(1.0, min(NDEM,NSUP))
C
         NUPT   = MAX(0.0, MIN(NDEM, NSUP))

         if (ndem.ge.nsup)then
	     maxdem = ndem
	   endif

         ndem     = ndem * 10.
	   nsup     = nsup * 10. 
	   nupt     = nupt * 10.
	   maxdem   = maxdem * 10.
	   crop_off = nupt
C stop the N 
         if (ndem.ge.nsup)then
c           take NDEM from the soil
c NB : This is in kg N / ha
           if (nlayers.le.10)then
		   do 100 i = 1, nlayers
	         soilno3(i) = soilno3(i) - MaxNO3(i) 
 100         continue
             soilnh4(1) = soilnh4(1) - MaxNH4(1)

		 elseif(nlayers.le.20)then
             do 110 i = 1,10
	         soilno3(i) = soilno3(i) - MaxNO3(i)
 110         continue
             soilno3(11) = soilno3(11) - MaxNO3(11) 
             soilnh4(1) = soilnh4(1) - MaxNH4(1)
	       soilnh4(2) = soilnh4(2) - MaxNH4(2)
		 else
	       do 120 i = 1,11
	         soilno3(i) = soilno3(i) - MaxNO3(i)
  120         continue
             soilno3(12) = soilno3(12) - MaxNO3(12)

		   soilnh4(1) = soilnh4(1) - MaxNH4(1)
		   soilnh4(2) = soilnh4(2) - MaxNH4(2)
		   soilnh4(3) = soilnh4(3) - MaxNH4(3)

		 endif   
         else
c           take demand from the soil leaving a bit in the soil
       		 NdemNO3 = ndem * TotAvNO3 / (TotAvNO3 + TotAvNH4)
	     NdemNH4 = ndem * TotAvNH4 / (TotAvNO3 + TotAvNH4)

           if (nlayers.le.10)then
		   do 130 i = 1, nlayers
             if (TotAvNO3.GT.0.)then
	         UptNO3 = MaxNO3(i)*(NdemNO3/TotAvNO3)
	       else 
	         UptNO3 = 0.
	       endif
	       soilno3(i) = soilno3(i) - UptNO3
 130         continue

             if (TotAvNH4.GT.0.)then
               UptNH4 = MaxNH4(1)*(NdemNH4/TotAvNH4)
	       else 
	         UptNH4 = 0.
	       endif
             soilnh4(1) = soilnh4(1) - UptNH4

		 elseif(nlayers.le.20)then
             do 140 i = 1,10
	         if (TotAvNO3.GT.0.)then
	           UptNO3 = MaxNO3(i)*(NdemNO3/TotAvNO3)
	         else 
	           UptNO3 = 0.
	         endif
	         soilno3(i) = soilno3(i) - UptNO3
 140         continue

	       if (TotAvNO3.GT.0.)then
               UptNO3 = MaxNO3(11)*(NdemNO3/TotAvNO3)
	       else 
	         UptNO3 = 0.
	       endif
             soilno3(11) = soilno3(11) -  UptNO3

	       if (TotAvNH4.GT.0.)then
               UptNH4 = MaxNH4(1)*(NdemNH4/TotAvNH4)
	       else 
	         UptNH4 = 0.
	       endif
             soilnh4(1) = soilnh4(1) - UptNH4

	       if(TotAvNH4.GT.0.)then 
	         UptNH4 = MaxNH4(2)*(NdemNH4/TotAvNH4)
             else
	         UptNH4 = 0.
       endif

	       soilnh4(2) = soilnh4(2) - UptNH4
		 else
	       do 150 i = 1,11
	       if (TotAvNO3.GT.0.)then
	         UptNO3 = MaxNO3(i)*(NdemNO3/TotAvNO3)
	       else 
	         UptNO3 = 0.
	       endif
             soilno3(i) = soilno3(i) - UptNO3
 150         continue

             if(TotAvNO3.GT.0.)then
               UptNO3 = MaxNO3(12)*(NdemNO3/TotAvNO3)
	       else 
	         UptNO3 = 0.
	       endif
             soilno3(12) = soilno3(12) - UptNO3

	       if(TotAvNH4.GT.0.)then 
               UptNH4 = MaxNH4(1)*(NdemNH4/TotAvNH4)
	       else
	         UptNH4 = 0.
	       endif
		   soilnh4(1) = soilnh4(1) - UptNH4

	       if(TotAvNH4.GT.0.)then 
               UptNH4 = MaxNH4(2)*(NdemNH4/TotAvNH4)
	       else
	         UptNH4 = 0.
	       endif
		   soilnh4(2) = soilnh4(2) - UptNH4

	       if(TotAvNH4.GT.0.)then 
	         UptNH4 = MaxNH4(3)*(NdemNH4/TotAvNH4)
	       else
	         UptNH4 = 0.
	       endif
		   soilnh4(3) = soilnh4(3) - UptNH4
       	 endif   
         endif

         ndem   = ndem / 10.
	   nsup   = nsup / 10. 
	   nupt   = nupt / 10.     
		 maxdem = maxdem / 10.  

         NNC    = NUPT/(PCANA*Y)
 
         RNTOT  = NUPT - LNLV - LNRT
 
         FWSH   = 1./(1.+NNC*DERI/SHSA)
 
         FNSH   = 1./(1.+NNC*DERI/SHSA*WSH/WRT*NRT/NSH)
         RGRRT  = PCANA*Y*(1.-FWSH)/WRT
 
         RWLV   = FWLV * FWSH * Y * ASSA - LWLV
         RNSO   = FWSO * FWSH * Y * ASSA * NONC
         RNST   =(FWST * FWSH * Y * ASSA - LSTR)*HNCMIN
         RWST   = FWST * FWSH * Y * ASSA - LSTR
         RWSO   = FWSO * FWSH * Y * ASSA
         RWRT   = (1. -  FWSH)* Y * ASSA - LWRT
         RGRSH  = PCANA*Y*FWSH/WSH
         RNRT   = (1.-FNSH)*NUPT - LNRT - 
     &            (NRT-wrt*rncmin)/(NLV+(NRT-wrt*rncmin))*RNSO
         RNSH   =    FNSH *NUPT - LNLV + 
     &               (NRT-wrt*rncmin)/(NLV+(NRT-wrt*rncmin))*RNSO
         RNLV   = RNSH - RNST - RNSO
         CALL GLA(DS,HNCMIN/SLNMIN,RWLV,LAI,KN,NLV,RNLV,SLNB,RSLNB,RLAI)
 
         DS    = DS    + DELT*DVR
         WLV   = WLV   + DELT*RWLV
         WLVD  = WLVD  + DELT*LWLV
         WST   = WST   + DELT*RWST
         WSO   = WSO   + DELT*RWSO
         WRT   = WRT   + DELT*RWRT
         WRTD  = WRTD  + DELT*LWRT
         WSTR  = WSTR  + DELT*LSTR
         WLVDS = WLVDS + DELT*LVDS
         NTOT  = NTOT  + DELT*RNTOT
         NSH   = NSH   + DELT*RNSH
         NRT   = NRT   + DELT*RNRT
         NST   = NST   + DELT*RNST
         NLV   = NLV   + DELT*RNLV
         NSO   = NSO   + DELT*RNSO
         SLNB  = SLNB  + DELT*RSLNB
         LAI   = LAI   + DELT*RLAI
         TSOIL = TSOIL + DELT*RTSOIL
         TSUM  = TSUM  + DELT*HU 

      write(52,5201)DOY, FNSH, FWSH, PT, AE, PCANA,SLA,NDEM,NSUP, NUPT            
 5201 format(1x, i5, 2f7.3, 2F6.1, F7.1, 4F8.4)     

	write(53, 5301) DOY, WLV, WLVD, WST, WSO, WRT, WRTD,
     &   WSTR, WLVDS,WSH
 5301 format(1x, i5, 9f7.1)      
       
	write(54, 5401) DOY, DS, NTOT, NSO, NRT, NSH, NLV, NST, LAI, 
     &    LAIN, SLN, TSOIL, TSUM 
 5401 format(1x, i5, f6.2, 8f5.1, f6.2, 2f7.1)
      
	write(55,5501)DOY, PNC, HNC, RNC, LNC, ONC, NNC, LITC, LITN, 
     &    DERI,RTSA, SHSA, ADIF
 5501 format(1x, i5, 9f6.2, f6.3, f6.2,f6.2)

         c_root     = pcc*wrt
         c_st_lv  = pcc*(wlv+wst+wlvd-wlvds)*GAMMA*BETA
	   c_annual = (c_root+c_st_lv)*10.

	   n_root = wrt*rnc
	   n_st_lv = (wlv*lnc + wst*hncmin + (wlvd-wlvds)*slnmin*sla)*
     &          GAMMA*BETA
	   n_annual = (n_root + n_st_lv)*10.
      
	   c_return = litc
	   n_return = litn

	   if(ik.eq.mend)then
	      c_return = c_annual
	      n_return = n_annual
	      c_annual = 0.
	      n_annual = 0.  
        endif
      else
ckwc: 16/3/00	  write(76,*) DOY, ' Maturity - Harvest '
        DTR    = (RDD/1.E+6) * COEFR
        DVP    = VP          * COEFV
        TMAX   = TMMX + COEFT
        TMIN   = TMMN + COEFT
        DAVTMP = 0.29*TMIN + 0.71*TMAX
        NAVTMP = 0.71*TMIN + 0.29*TMAX
        
	  DTRJM2 = DTR * 1.E6
        TAVS   = (DAVTMP+NAVTMP)/2.
        BBRADS = BOLTZM * (TAVS+273.)**4 * 86400.
        SVPS   = 0.611 * EXP(17.4 * TAVS / (TAVS + 239.))
        SLOPES = 4158.6 * SVPS / (TAVS + 239.)**2
        RLWNS  = BBRADS * MAX(0., 0.55*(1.-DVP/SVPS))
        NRADS  = DTRJM2 * (1.-0.15) - RLWNS
        PENMRS = NRADS * SLOPES/(SLOPES+PSYCH)
        WDF    = 2.63*(1.0 + 0.54 * WN)
        PENMD  = LHVAP*WDF*MAX(0.,(SVPS-DVP))*PSYCH/(SLOPES+PSYCH)

        PE     = MAX(0., (PENMRS + PENMD)/LHVAP)
	  PE   = PE * (exp(-0.5*(LAI+DLAI)))
        AE     = MIN (PE, soilh2o(1))

        evapo = AE

	  c_return = 0.0
	  n_return = 0.0

	  if(ik.eq.mend)then
          c_return = c_annual
	    n_return = n_annual
	    c_annual = 0.
	    n_annual = 0.  
        endif

      END IF
 
      RETURN
      END
 
*----------------------------------------------------------------------*
*  SUBROUTINE SUBDD                                                    *
*  Purpose: This subroutine calculates the daily amount of heat units. *
*----------------------------------------------------------------------*
      SUBROUTINE SUBDD(TMAX,TMIN,DIF,DAYL,TBD,TOD,HU)

      IMPLICIT NONE

C REAL passed from calling subroutine
      REAL TMAX,TMIN,DIF,DAYL,TBD,TOD,HU

C INTEGER local to subroutine
      INTEGER I

C REAL local to subroutine
      REAL SUNRIS, SUNSET, TM, TT, TD

      SAVE
 
      SUNRIS = 12. - 12.*DAYL
      SUNSET = 12. + 12.*DAYL
 
      TM    = (TMAX + TMIN)/2.
      TT    = 0.
      DO 10 I = 1, 24
        IF (I.GE.SUNRIS .AND. I.LE.SUNSET) THEN
          TD = TM + DIF + 0.5*ABS(TMAX-TMIN)*COS(0.2618*FLOAT(I-14))
        ELSE
          TD = TM       + 0.5*ABS(TMAX-TMIN)*COS(0.2618*FLOAT(I-14))
        ENDIF

        IF (TD.LT.TBD)THEN
	    TD = TBD
	  ELSE        
	    TD = MIN(TD, TOD)
	  ENDIF

        TT = TT + (TD-TBD)/24.
 
  10  CONTINUE

      HU = TT
 
      RETURN
      END
 
*----------------------------------------------------------------------*
*  SUBROUTINE PHENOL                                                   *
*  Purpose: This subroutine calculates phenological development rate.  *
*----------------------------------------------------------------------*
      SUBROUTINE PHENOL (DS,nsow,ik,DVRV,
     &                   DVRR,HU,DVR)

      IMPLICIT NONE

C INTEGER passed from calling subroutine
	INTEGER nsow, ik

C REAL passed from calling subroutine
	REAL DS, DVR, DVRV, DVRR, HU

      SAVE
 
      IF (nsow.eq.ik) THEN
         DVR    = 0.
      END IF
 
      IF (DS.GE.0. .AND. DS.LT.1.0) THEN
          DVR   = DVRV*HU
      ELSE
          DVR   = DVRR*HU
      ENDIF
 
      RETURN
      END
 
*----------------------------------------------------------------------*
* SUBROUTINE PCANF                                                     *
* Purpose: To compute canopy photosynthesis based on Farquhar's        *
*          biochemical leaf-photosynthesis model                       *
*----------------------------------------------------------------------*
      SUBROUTINE PCANF(C3C4,CO2I,SLNMIN,SLN,K,LAI,TLEAF,DAYL,DTR,PCAN)

      IMPLICIT NONE


C REAL passed from calling subroutine
      REAL C3C4,CO2I,SLNMIN,SLN,K,LAI,TLEAF,DAYL,DTR,PCAN

C REAL local to subroutine
      REAL O2, JMUMOL, KC25, KMC25, KMO25, KOKC, EAVCMX, EAKMC, EAKMO
      REAL PAR, PINT, RUBISC, VCMAX, KMC, KMO, GAMMAX
      REAL CO2IP, EFF, PMAX

      O2     = 21.
      JMUMOL = 4.56
      KC25   = 138.
      KMC25  = 460.
      KMO25  = 33.
      KOKC   = 0.21
      EAVCMX = 68000.
      EAKMC  = 65800.
      EAKMO  = 1400.
 
      PAR    = 0.5*DTR
      PINT   = 1.-EXP(-K*LAI)
 
      RUBISC = MAX (0., 2.21*(SLN-SLNMIN))
      VCMAX  = RUBISC*KC25*EXP((1./298.-1./(TLEAF+273.))*EAVCMX/8.314)
      KMC    = KMC25*EXP((1./298.-1./(TLEAF+273.))*EAKMC/8.314)
      KMO    = KMO25*EXP((1./298.-1./(TLEAF+273.))*EAKMO/8.314)
      GAMMAX = 0.5 * KOKC * KMC * O2 / KMO
 
      IF (C3C4.GE.0.) THEN
        CO2IP = CO2I 
		EFF    = 44.*JMUMOL/2.1*(CO2IP-GAMMAX)/(4.5*CO2IP+10.5*GAMMAX)
      ELSE
        CO2IP = 1500.
        EFF    = 44.*JMUMOL/2.1/7.5
      ENDIF

      PMAX   = VCMAX*(CO2IP-GAMMAX)/(CO2IP+KMC*(1.+O2/KMO))
      
      PCAN   = EFF*PMAX/(EFF*K*PAR/DAYL+PMAX)*PAR*PINT
 
      RETURN
      END
 
*----------------------------------------------------------------------*
* SUBROUTINE PTRAN                                                     *
* Purpose: To compute canopy transpiration, using Penman-Monteith      *
*          equation                                                    *
*----------------------------------------------------------------------*
      SUBROUTINE PTRAN (RSW,RT,RBW,RBH,BOLTZM,LHVAP,VHCA,PSYCH,LAI,
     &                DTRJM2,TLEAF,DAYL,DVP,SLOPE,SVP,VPD,PT,NRADC)

      IMPLICIT NONE

C REAL passed from calling subroutine       
	REAL RSW,RT,RBW,RBH,BOLTZM,LHVAP,VHCA,PSYCH,LAI
      REAL DTRJM2,TLEAF,DAYL,DVP,SLOPE,SVP,VPD,PT,NRADC

C REAL local to subroutine
      REAL BBRAD, RLWN, PSR, PTR, PTD   
 
      BBRAD  = BOLTZM * (TLEAF +273.)**4 * 86400.*DAYL
      RLWN   = BBRAD  * MAX(0., 0.55*(1.-DVP/SVP))
      NRADC = (1.-EXP(-0.5*LAI))*(DTRJM2*(1.-0.25) - RLWN)
 
      PSR    = PSYCH*(RBW+RT+RSW)/(RBH+RT)
 
      PTR    = NRADC*SLOPE        /(SLOPE+PSR)/LHVAP
      PTD    = (VHCA*VPD/(RBH+RT))/(SLOPE+PSR)/LHVAP*86400.*DAYL
 
      PT     = PTR + PTD
 
      RETURN
      END
 
*----------------------------------------------------------------------*
*  SUBROUTINE DIFLA                                                    *
*  Purpose: This subroutine calculates leaf(canopy)-air temperature    *
*           differential (oC).                                         *
*----------------------------------------------------------------------*
      SUBROUTINE DIFLA (LHVAP,VHCA,NRADC,PT,DAYL,RBH,RT,DIF)

      IMPLICIT NONE

C REAL passed from calling subroutine    
      REAL LHVAP,VHCA,NRADC,PT,DAYL,RBH,RT,DIF

      SAVE
 
      DIF   = (NRADC-LHVAP*PT)/86400./DAYL*(RBH+RT)/VHCA
 
      RETURN
      END
 
*----------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine calculates the daily increase of leaf      *
*           area index (m2 leaf/m2 ground/day).                        *
*----------------------------------------------------------------------*
      SUBROUTINE GLA(DS,SLA,RWLV,LAI,KN,NLV,RNLV,SLNB,RSLNB, RLAI)

      IMPLICIT NONE

C REAL passed from calling subroutine
      REAL DS,SLA,RWLV,LAI,KN,NLV,RNLV,SLNB,RSLNB, RLAI

      SAVE
 
      RLAI = SLA*RWLV
 
      IF ((LAI.LT.1.) .AND. (DS.LT.1.)) THEN
        RLAI = (SLNB*RNLV-NLV*RSLNB)/SLNB/(SLNB+KN*NLV)
      ENDIF
 
      RETURN
      END

*----------------------------------------------------------------------*
*  SUBROUTINE SETTIME_MAGEC_CROP                                       *
*  Purpose: This subroutine Set time factors from seconds / timestep   *
*----------------------------------------------------------------------*
      SUBROUTINE SETTIME_MAGEC_CROP(SECONDS,CONVER_F,N_STEPS)
      IMPLICIT NONE
C
C Variables passed to/from calling subroutine
C
      INTEGER N_STEPS
      REAL CONVER_F,SECONDS
C
C Set time factors from SECONDS
C
      N_STEPS=(365.25*24*60*60)/SECONDS
	CONVER_F=SECONDS/(7*24*60*60)
	END
