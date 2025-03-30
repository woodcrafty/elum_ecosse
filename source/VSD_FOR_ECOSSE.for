
      subroutine  Run_VSD(start, cNO3in, cNH4in, cDOCin,  
     &                    waterin,waterout, thick,pH, cAl, SoilSeries,
     &                    VSD_states)

! Declare VSD parameters 
      logical(4)         start
      integer(4)         modexc

      real(4)            thick, rho, theta, Cpool, ps, pCO2, fde, CEC
      real(4)            expAl, KAlox, KAlBc, KHBc, coacid, pKpar(3)
      real(4)            Nimacc, Caw, Mgw, Kw, Naw, Nu, Cau, Mgu, Ku
      real(4)            Sdep, Ndep, Cadep, Mgdep, Kdep, Nadep, Cldep
      real(4)            CNmin, CNmax, CNrat, Nim, Nde
      real(4)            cH, cSO4, cNO3, cBc, cNa, cCl
      real(4)            cAl, cOrg, cHCO3, cANC, EBc, EAl, EH

      byte               errbyte(64)

      real(4)		 cNO3in, cNH4in, cDOCin, waterin, waterout, pH
      integer(4)         SoilSeries, i
      real 		VSD_states
      dimension 	VSD_states(10)! array to store VSD states 

! Declare VSD subroutine
      external VSD

! Assign values to VSD parameters
! Values taken from VSD run for gridcell 285154 acid grassland
       thick =1.0
       rho =1.248539
       theta =0.3
       pCO2 =20.0
       expAl =3.0
       KAlox =8.5
       CEC =124.2585
       modexc =2
       KAlBc =1.549773
       KHBc =4.370109
       fde =0.1552648
       ps = 0.698
       coacid =0.025
       pKpar(1) =0.96
       pKpar(2) =0.9
       pKpar(3) =0.039
       Caw = 0.035
       Mgw = 0.105
       Kw = 0.105
       Naw = 0.105
       Nu =0.0
       Cau =0.0
       Mgu =0.0
       Ku =0.0
       Sdep = 0.0706
       Ndep = 0.0339
       Cadep = 0.0100
       Mgdep = 0.0118
       Kdep = 0.0089
       Nadep = 0.4051
       Cldep = 0.4710
       Nimacc =71.0
       Cpool =13294.0
       CNmin =7.5
       CNmax =20.8
       CNrat =15.02122
       Nim =2.0
       Nde =0.0
       cH = 0.01
       cSO4 =0.01
       cNO3 =0.01
       cNH4 =0.005
       cBc =0.08
       cNa =0.10
       cCl =0.11
       cOrg =0.005
       cHCO3 =0.02
       cANC =0.055
       EBc =0.001
       EAl =0.1
       EH =0.05
      do i=1,64
      		errbyte(i) = 0
      enddo
 
      call VSD(start,thick,rho,theta,pCO2,expAl,KAlox,CEC, 
     &             modexc,KAlBc,KHBc,fde,ps,coacid,pKpar, 
     &             Caw,Mgw,Kw,Naw,Nu,Cau,Mgu,Ku,Sdep,Ndep,Cadep, 
     &             Mgdep,Kdep,Nadep,Cldep,Nimacc,Cpool,CNmin, 
     &             CNmax,CNrat,Nim,Nde,cH,cSO4,cNO3,cBc,cNa,cCl, 
     &             cAl,cOrg,cHCO3,cANC,EBc,EAl,EH,errbyte,
     &		   cNO3in, cNH4in,VSD_states)

!      write(6,'("Cpool after VSD",f10.3)')Cpool     
 
      pH = -1 * log(cH*0.001)	! cH is in eq/m3, pH uses eq/dm3
      
!      write(6,'(" ")')
!      write(6,'("end of Call_VSD")')
!      write(6,'("cNO3in ",f10.3)')cNO3in
!      write(6,'("cNH4in ",f10.3)')cNH4in
!      write(6,'("cDOCin ",f10.3)')cDOCin
!      write(6,'("waterin ",f10.3)')waterin
!      write(6,'("waterout ",f10.3)')waterout
!      write(6,'("pH ",f10.3)')pH
!      write(6,'("cAl ",f10.3)')cAl

      END ! sub Call_Run_VSD


