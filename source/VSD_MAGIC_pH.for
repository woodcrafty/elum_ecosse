
      subroutine  VSD(start,thick,rho,theta,pCO2,expAl,KAlox, 
     &             CEC,modexc,KAlBc,KHBc,fde,ps,coacid,pKpar,  
     &             Caw,Mgw,Kw,Naw,Nu,Cau,Mgu,Ku,Sdep,Ndep,Cadep,   
     &             Mgdep,Kdep,Nadep,Cldep,Nimacc,Cpool,CNmin,   
     &             CNmax,CNrat,Nim,Nde,cH,cSO4,cNO3,cBc,cNa,cCl,  
     &             cAl,cOrg,cHCO3,cANC,EBc,EAl,EH,errbyte,
     &		   cNO3in, cNH4in,VSD_states)
!
!     Very Simple Dynamic (VSD) [soil acidification] model:
!     Calculates for a single time step a number of soil chemical
!     parameters depending on element inputs and chemical soil status.
!
!     INPUT variables:
!

!     start ..... start flag; if true: initialise certain variables for
!                 use in simulation; set to false afterwards (see Notes)
!     thick ..... soil thickness (m)
!     rho ....... bulk density (g/cm3)
!     theta ..... volumetric water content of the soil (m3/m3)
!     pCO2 ...... partial CO2-pressure in soil solution (atm)
!     expAl ..... exponent in [Al] = KAlox*[H]^expAl
!     KAlox ..... Al equilibrium constant ((eq/m3)^(1-expAl))
!     CEC ....... cation exchange capacity of the soil (meq/kg)
!     modexc .... cation exchange model option: 1=Gaines-Thomas, 2=Gapon
!     KAlBc ..... selectivity constant for Al-Bc exchange
!     KHBc ...... selectivity constant for H-Bc exchange
!     fde ....... denitrification factor (0<=fde<=1)
!     ps ........ precipitation surplus (m/yr)
!     coacid .... total concentration of organic acids (m*DOC) (mol/m3)
!     pKpar() ... 1-3 parameters of (Oliver-type) mono-protic organics model:
!                 pK = par(21)+par(22)*pH-par(23)*pH^2
!     Caw ....... weathering rate of Ca (eq/m3/yr)
!     Mgw ....... weathering rate of Mg (eq/m3/yr)
!     Kw ........ weathering rate of K (eq/m3/yr)
!     Naw ....... weathering rate of Na (eq/m3/yr)
!     Nu ........ net uptake of N (eq/m2/yr)
!     Cau ....... net uptake of Ca (eq/m2/yr)
!     Mgu ....... net uptake of Mg (eq/m2/yr)
!     Ku ........ net uptake of K (eq/m2/yr)
!     Sdep ...... S  deposition (eq/m2/yr)
!     Ndep ...... N  deposition (eq/m2/yr)
!     Cadep ..... Ca deposition (eq/m2/yr)
!     Mgdep ..... Mg deposition (eq/m2/yr)
!     Kdep ...... K  deposition (eq/m2/yr)
!     Nadep ..... Na deposition (eq/m2/yr)
!     Cldep ..... Cl deposition (eq/m2/yr)
!     Nimacc .... constant rate of N immobilization (eq/m2/yr)
!     CNmin ..... C:N ratio (g/g) below which no time-dependent Nimm
!     CNmax ..... C:N ratio (g/g) above which all available N is immobilised
!
!     OUTPUT variables:
!
!    +Cpool ..... amount of C in topsoil (g/m2)
!    +CNrat ..... C:N ratio in topsoil (g/g)
!     Nim ....... N immobilized (eq/m2/yr)
!     Nde ....... N denitrified (eq/m2/yr)
!     cH ........ H+ concentration (eq/m3)
!    *cSO4 ...... SO4 concentration (eq/m3)
!    *cNO3 ...... NO3 concentration (eq/m3)

! Ed's additions!!!!!!!!!!
!    *cNH4 ...... NH4 concentration (eq/m3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!    *cBc ....... Ca+Mg+K ('divalent') base cation concentration (eq/m3)
!    *cNa ....... Na concentration (eq/m3)
!    *cCl ....... Cl concentration (eq/m3)
!     cAl ....... aluminium concentration (eq/m3)
!     cOrg ...... concentration of organic anions (eq/m3)
!     cHCO3 ..... bicarbonate concentration (eq/m3)
!     cANC ...... ANC concentration (as function of [H]) (eq/m3)
!   +*EBc ....... exchangeable fraction of Ca+Mg+K ('base saturation')
!     EAl ....... exchangeable fraction of Al
!     EH ........ exchangeable fraction of H (EBc+EAl+EH=1)
!     errbyte ... byte-vector holding error message from 'root01':
!                 if char(errbyte(1)) not blank, there IS an error(message)!
!
!     NOTES:
!    +These variables have to hold the values from the previous timestep;
!     i.e. at the first call they have to be the initial values.
!    *These variables have to hold the values from the previous timestep,
!     UNLESS start=.true.: then they are initialized with equilibrium values
!                          (EBc only if <0).
!
!DEC$ATTRIBUTES DLLEXPORT :: VSD
!
      implicit           none
!
      integer(4)         LE
      parameter          (LE=64) ! length of 'errbyte'/'errmsg'
      real(4)            eps, Hmin, third
      parameter          (eps=1.e-6,Hmin=eps,third=1./3.)
!
      byte               errbyte(LE)
      logical(4)         start
      integer(4)         modexc
      real(4)            thick, rho, theta, Cpool, ps, pCO2, fde, CEC
      real(4)            expAl, KAlox, KAlBc, KHBc, coacid, pKpar(3)
      real(4)            Nimacc, Caw, Mgw, Kw, Naw, Nu, Cau, Mgu, Ku
      real(4)            Sdep, Ndep, Cadep, Mgdep, Kdep, Nadep, Cldep
      real(4)            CNmin, CNmax, CNrat, Nim, Nde
      real(4)            cH, cSO4, cNO3, cBc, cNa, cCl
      real(4)            cAl, cOrg, cHCO3, cANC, EBc, EAl, EH
!
!U    USES chbal, fAlH, fANC, fgtn, fOliver, giant0, giant, root01
!
!      logical(4)         do_CN
      character(len=LE)  errmsg
      integer(4)         n
      real(4)            par(40), SO4in, NO3in, Bcin, Nain, Clin
      real(4)            dt, KHCO3, cBcmin, fcap, CECa, wat
      real(4)            cOth, cSO4o, cNO3o, cNao, cClo, cBco, EBco
!      real(4)            Npool, fCN, Nav, Niminf, Nimt
      real(4)            chbal, fAlH, fANC, fgtn, fOliver
      real(4)            giant0, giant, root01
!
      external           chbal, fAlH, fANC, fgtn, fOliver
      external           giant0, giant, root01
!
! Ed's additions!!!!!!!!!
      real 		VSD_states
      real(4)           cNO3in, cNH4in, cNH4, NH4in
      dimension 	VSD_states(10)! array to store VSD states 
!!!!!!!!!!!!!!!!!!!!!!!!!

!      save               do_CN
!
      data  dt /1./        ! 1 year time-step
      data  KHCO3 /0.02/   ! = 10^-1.7 (mol/m3)^2/atm (at 8oC)
      data  cBcmin /0.001/ ! eq/m3
!
      do n = 1,LE
        errmsg(n:n) = ' '
      end do

! Ed's additions!!!!!!!!!
	if (start.eqv..false.) then ! Re-assign values from previous timestep 
		Cpool 	= VSD_states(1)
		CNrat 	= VSD_states(2) 
		cSO4 	= VSD_states(3)
		cNO3 	= VSD_states(4)
		cBc 	= VSD_states(5)
		cNa	    = VSD_states(6)
		cCl 	= VSD_states(7)
		EBc 	= VSD_states(8)
		cNH4 	= VSD_states(9)
	end if 
!!!!!!!!!!!!!!!!!!!!!!!!!

      fcap = theta*thick
      wat = fcap+ps*dt
      CECa = CEC*thick*rho
!
      par(1) = KHCO3*pCO2
      par(2) = expAl
      par(3) = KAlox
      par(4) = wat/CECa
      par(6) = KHBc
      par(7) = KAlBc
      par(17) = real(modexc)
      par(19) = coacid
      do n = 1,3
        par(20+n) = pKpar(n)
      end do
!
! - N immobilization and C:N ratio:
!      Nav = max(0.,Ndep-Nu)
!      Npool = 0.
!      if (CNrat > 0.) Npool = (Cpool/CNrat)/14. ! eq/m2
!      Niminf = 0.
!      if (Nimacc > 0.) then ! constant N immobilisation:
!        Niminf = min(Nimacc,Nav)
!        Cpool = Cpool+14.*CNrat*Niminf ! update C-pool to keep CNrat the same!
!        Npool = Npool+Niminf
!      end if
!      Nimt = 0.
!      if (start .and. Cpool > 0.) do_CN = .true.
!      if (do_CN) then
!        if (CNrat > CNmin) then ! C:N-dependent N immobilisation:
!          fCN = min(max(0.,(CNrat-CNmin)/(CNmax-CNmin)),1.)
!          Nimt = fCN*max(0.,Nav-Niminf)
!          Npool = Npool+Nimt
!          if (Npool > 0.) CNrat = Cpool/(14.*Npool) ! g/g
!        end if ! (no update of C-pool!)
!      end if
!      Nim = Niminf+Nimt
!
! - denitrification:
!      Nde = fde*max(0.,Nav-Nim)
!
! - net input fluxes:
      SO4in = Sdep
      NO3in = cNO3in*ps
      NH4in = cNH4in*ps 

      Bcin  = max(ps*cBcmin,max(0.,Cadep+Caw*thick-Cau)
     &                     +max(0.,Mgdep+Mgw*thick-Mgu)
     &                     +max(0.,Kdep +Kw *thick-Ku))
      Nain  = Nadep+Naw*thick
      Clin  = Cldep
!


      if (start) then   ! initialise (with equilibrium):
        cSO4o = SO4in/ps
        cNO3o = NO3in/ps
        cNao  = Nain/ps
        cClo  = Clin/ps
        cOth  = cSO4o+cNO3o-cNao+cClo-cNH4
        par(11) = cOth
        if (EBc < 0.) then ! initialize (also) EBc:
          cBc = Bcin/ps
          par(12) = cBc
          errmsg = ' chbal'
          cH = root01(Hmin,10.,par,eps,chbal,errmsg)
          if (errmsg(1:1) /= ' ')       goto 99
          cAl = fAlH(cH,par)
          if (modexc == 1) then ! Gaines-Thomas
            par(13) = cH*sqrt(KHBc/cBc)
            par(14) = cAl*sqrt(KAlBc/cBc)/cBc
            errmsg = ' fgtn'
            EBc = root01(0.,1.,par,eps,fgtn,errmsg)
            if (errmsg(1:1) /= ' ')     goto 99
          else ! Gapon
            EBc = 1./(1.+(KHBc*cH+KAlBc*cAl**third)/sqrt(cBc))
          end if
        else   ! use given EBc
          par(13) = EBc
          errmsg = ' giant0'
          cH = root01(Hmin,10.,par,eps,giant0,errmsg)
          if (errmsg(1:1) /= ' ')       goto 99
        end if
        cANC = fANC(cH,par)
        cBco = max(0.,cOth+cANC)
        EBco = EBc
        start = .false.
      else   ! use variables from previous time step:
        cSO4o = cSO4
        cNO3o = cNO3
        cBco = cBc
        cNao = cNa
        cClo = cCl
        EBco = EBc
      end if
!
! - Update concentrations in soil solution:
!
      cSO4 = (fcap*cSO4o+SO4in*dt)/wat
      cNO3 = (fcap*cNO3o+NO3in*dt)/wat
      cNa  = (fcap*cNao +Nain *dt)/wat
      cCl  = (fcap*cClo +Clin *dt)/wat
      cOth = cSO4+cNO3+cCl-cNa-cNH4
      par(11) = cOth
!
      par(14) = EBco+(fcap*cBco+Bcin*dt)/CECa
      errmsg = ' giant'
      cH = root01(Hmin,10.,par,eps,giant,errmsg)
      
      if (errmsg(1:1) /= ' ')           goto 99
      cANC = fANC(cH,par)
      cBc = max(0.,cOth+cANC)
      EBc = min(max(0.,par(14)-cBc*par(4)),1.)
      if (modexc == 1) then ! Gaines-Thomas
        EH = min(cH*sqrt(EBc*KHBc/cBc),1.-EBc)
      else   ! Gapon
        EH = min(cH*EBc*KHBc/sqrt(cBc),1.-EBc)
      end if
      EAl = max(0.,1.-EBc-EH)


      cAl = fAlH(cH,par)
      cOrg = fOliver(cH,par)
 
      cHCO3 = par(1)/cH

99    continue
      do n = 1,LE-1
 !!!The following line commented out as causes crashes 
 !       errbyte(n) = ichar(errmsg(n:n))
      end do
		
! Ed's additions!!!!!!!!!
! store values at end of current timestep
	VSD_states(1) = Cpool
	VSD_states(2) = CNrat
	VSD_states(3) = cSO4 
	VSD_states(4) = cNO3 	
	VSD_states(5) = cNa
	VSD_states(6) = cBc 	
	VSD_states(7) = cCl 
	VSD_states(8) = EBc
	VSD_states(9) = cNH4

!print out states
c      write(6,'(" ")')
c      write(6,'("Cpool : ",f10.6)')VSD_states(1)
c      write(6,'("CNrat : ",f10.6)')VSD_states(2)
c      write(6,'("cSO4  : ",f10.6)')VSD_states(3)
c      write(6,'("cNO3  : ",f10.6)')VSD_states(4)
c      write(6,'("cNa   : ",f10.6)')VSD_states(5)
c      write(6,'("cBc   : ",f10.6)')VSD_states(6)
c      write(6,'("cCl   : ",f10.6)')VSD_states(7)
c      write(6,'("EBc   : ",f10.6)')VSD_states(8)
c      write(6,'("cNH4  : ",f10.6)')VSD_states(9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!The following line commented out as causes crashes 
 !      errbyte(LE) = 0
                                          return
      end subroutine  VSD
!
!
      real function  fgtn  (EBc,par)
!
!
!     Solving Gaines-Thomas equations for EBc:
!
!     par(13) = cH*sqrt(KHBc/cBc)
!     par(14) = cAl*sqrt(KAlBc/cBc)/cBc
!
      implicit           none
!
      real(4)            EBc, par(*)
!
      fgtn = 1.-EBc-sqrt(EBc)*(par(13)+EBc*par(14))
                                         return
      end function  fgtn
!
!
      real function  fAlH  (cH,par)
!
!
!     par(2)  = expAl
!     par(3)  = KAlox
!
      implicit           none
!
      real(4)            cH, par(*)
!
      fAlH = 0.
      if (cH > 0.) fAlH = par(3)*cH**par(2) ! = cAl
                                        return
      end function  fAlH
!
!
      real function  fOliver  (cH,par)
!
!
!     Returns concentration of mono-protic organic anions in mol/m3=eq/m3
!     calculated with an Oliver-type model.
!
!     cH .......... [H+] (eq/m3=mol/m3)
!     par(19) ..... total concentration of organic acids (m*DOC) (mol/m3)
!     par(21)+par(22)*pH-par(23)*pH^2 = pK ... Oliver (type) model
!       If par(22)=par(23)=0: simple mono-protic model with pK = par(21).
!
      implicit           none
!
      real(4)            cH, par(*)
      real(4)            pH, pKa, dis
!
      fOliver = 0.
      if (par(19) <= 0. .or. cH <= 0.)  return
      pH = 3.-alog10(cH)
      pKa = par(21)+pH*(par(22)-pH*par(23))
      dis = 10.**(3.-pKa)
      fOliver = par(19)*dis/(dis+cH)
                                        return
      end function  fOliver
!
!
      real function  fANC  (cH,par)
!
!
!     [ANC] (eq/m3) as a function of cH=[H+]:
!
!     par(1) = KHCO3*pCO2
!
      implicit           none
!
      real(4)            cH, par(*)
      real(4)            fAlH, fOliver, cAl, cOrg, cHCO3
!
      external           fAlH, fOliver
!
      cAl = fAlH(cH,par)
      cOrg = fOliver(cH,par)
      cHCO3 = par(1)/cH
      fANC = cHCO3+cOrg-cH-cAl
                                         return
      end function  fANC
!
!
      real function  chbal  (cH,par)
!
!
!     Charge balance equation expressed in cH=[H+]:
!
!     par(11) = cOth = [SO4--]+[NO3-]-[Na+]+[Cl-]-[NH4]
!     par(12) = cBc = [Ca+Mg+K]
!
      implicit           none
!
      real(4)            cH, par(*)
      real(4)            fANC, cANC
!
      external           fANC
!
      cANC = fANC(cH,par)
      chbal = par(12)-par(11)-cANC
                                         return
      end function  chbal
!
!
      real function  giant0  (cH,par)
!
!
!     par(6) = KHBc
!     par(7) = KAlBc
!     par(11) = cOth = [SO4]+[NO3]-[Na]+[Cl]-[NH4]
!     par(13) = EBc
!     par(17) = modexc (1=Gaines-Thomas, 2=Gapon)
!
      implicit           none
!
      real(4)            third
      parameter          (third=1./3.)
!
      real(4)            cH, par(*)
      real(4)            cAl, HBC, sHBC, bsat, sKHBc, sKAlBc, cAlth
      real(4)            fAlH, fANC
!
      external           fAlH, fANC
!
      cAl = fAlH(cH,par)
      HBC = cH*(par(11)+fANC(cH,par)) ! cH*[BC] from charge balance
      HBC = max(0.,HBC)
      sHBC = sqrt(HBC)
      if (nint(par(17)) == 1) then ! Gaines-Thomas
        bsat = cH*par(13)
        sKHBc = sqrt(bsat*par(6))
        sKAlBc = sqrt(bsat*par(7))
        giant0 = (1.-par(13))*HBC*sHBC-cH*sKHBc*HBC-sKAlBc*cAl*bsat
      else   ! Gapon
        cAlth = cAl**third
        giant0 = sHBC-par(13)*(sHBC+sqrt(cH)*(par(6)*cH+par(7)*cAlth))
      end if
                                        return
      end function  giant0
!
!
      real function  giant  (cH,par)
!
!
!     par(4) = (fcap+ps*dt)/CECa
!     par(6) = KHBc
!     par(7) = KAlBc
!     par(11) = cOth = [SO4]+[NO3]-[Na]+[Cl]-[NH4]
!     par(14) = EBco+(fcap*cBco+Bcin*dt)/CECa
!     par(17) = modexc (1=Gaines-Thomas, 2=Gapon)
!
      implicit           none
!
      real(4)            third
      parameter          (third=1./3.)
!
      real(4)            cH, par(*)
      real(4)            fAlH, fANC, cAl
      real(4)            HBC, sHBC, bsat, sKHBc, sKAlBc, cAlth
!
      external           fAlH, fANC
!
      cAl = fAlH(cH,par)
      HBC = cH*(par(11)+fANC(cH,par)) ! cH*[BC] from charge balance
      HBC = max(0.,HBC)
      sHBC = sqrt(HBC)
      bsat = max(0.,cH*par(14)-HBC*par(4)) ! = cH*EBc
      if (nint(par(17)) == 1) then ! Gaines-Thomas
        sKHBc = sqrt(bsat*par(6))
        sKAlBc = sqrt(bsat*par(7))
        giant = (cH-bsat)*sHBC*HBC-cH*(cH*sKHBc*HBC+sKAlBc*cAl*bsat)
      else   ! Gapon
        cAlth = cAl**third
        giant = (cH-bsat)*sHBC-bsat*sqrt(cH)*(par(6)*cH+par(7)*cAlth)
      end if
                                        return
      end function  giant
!
!
      real function  root01  (x1,x2,par,tol,func,errmsg)
!
!
!     Returns the zero of the function 'func' -- with parameter
!     vector par() --, known to lie between x1 and x2, with an accuracy 'tol'
!     using bisection, limited by 'itmax' iterations.
!     errmsg ... error message:
!                errmsg(1:1) = ' ' if no error.
!                errmsg(1:1) = 'E' & errmsg(10:) = ...
!     1. = 'ROOT01: root not bracketed.', i.e. func(a)*func(b)&gt;0
!          (no or &gt;1 zeros: the value of root01 is undetermined)
!     2. = 'ROOT01: exceeding maximum iterations.', i.e. itmax exceeded
!          (root01 returns the best approximation it could find)
!
      implicit           none
!
      integer(4)         itmax
      parameter          (itmax=100)
!
      character*(*)      errmsg
      real(4)            x1, x2, par(*), tol, func
!
      external           func
!
!U    USES NOTHING
!
      integer(4)         n
      real(4)            a, b, c, fa, fb, fc
!
      a = x1
      b = x2
      fa = func(a,par)
      fb = func(b,par)
      if (fa*fb > 0.) then
        errmsg(1:1) = 'E'
        errmsg(10:) = ' ROOT01: root not bracketed.'
                                        return
      end if
      do n = 1,itmax
        c = 0.5*(a+b)
        fc = func(c,par)
        if (abs(b-a) < tol) then
          root01 = c
                                        return
        else if (fa*fc > 0.) then
          a = c
          fa = fc
        else ! fb*fc > 0
          b = c
!          fb = fc
        end if
      end do
      errmsg(1:1) = 'E'
      errmsg(10:) = ' ROOT01: exceeding maximum iterations.'
      root01 = 0.5*(a+b)
                                        return
      end function  root01
!</PRE></BODY></HTML>

