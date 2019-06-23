      program o2profile
      !! try to calculate o2 profile in a simple system
      !! full implicit method
      !! from ver 8 reflecting ph change in next step of pyrite weathering and assuming always wet surface 
      !! from ver 9 deterministic pH calculation; pH is not iteratively sought 
      !! from ver 9.5 pH is iteratively calculated to be consistent with pyrite; for the case excluding aqueous Fe2 oxidation 
      !! from ver 9.6 calculation is conducted within the whole domain 
      implicit none
      
      !-----------------------------
      
      double precision :: ztot = 300.0d0 ! m
      double precision dz
      integer, parameter :: nz = 300 
      double precision z(nz)
      double precision ze(nz+1)
      double precision :: ph = 5.0d0
      double precision :: tc = 15.0d0 ! deg celsius
      double precision dt  ! yr 
      integer, parameter :: nt = 50000000
      double precision time
      integer, parameter :: nsp = 4
      integer, parameter :: nsp3= 2
      
      double precision :: rg = 8.3d-3   ! kJ mol^-1 K^-1
      double precision :: rg2 = 8.2d-2  ! L mol^-1 atm K^-1
      
      double precision :: po2i = 0.21d0 ! atm
C       double precision :: po2i = 0.6d-1 ! atm
      double precision :: pco2i = 10.0d0**(-2.5d0) ! atm
C       double precision :: pco2i = 10.0d0**(-1.0d0) ! atm
      double precision :: ci = 0d0 
      double precision :: c2i = 0d0
      double precision :: so4i = 0d0
      double precision :: nai = 0d0
      
      double precision :: redsldi = 0.56d0 ! wt%
C       double precision :: redsldi = 1.12d0 ! wt%
C       double precision :: redsldi = 2.8d0 ! wt%
C       double precision :: redsldi = 2.24d0 ! wt%
C       double precision :: redsldi = 3.36d0 ! wt%
      
      double precision :: silwti = 30d0 ! wt%
C       double precision :: silwti = 24d0 ! wt%
C       double precision :: silwti = 15d0 ! wt%
      
      double precision sat(nz), poro(nz), torg(nz), tora(nz), deff(nz)
      double precision :: dgas = 6.09d2 ! m^2 yr^-1
      double precision :: daq = 5.49d-2 ! m^2 yr^-1
      double precision :: poroi = 0.1d0
      double precision :: sati = 0.50d0
      double precision :: satup = 0.10d0
      
C       double precision :: zsat = 30d0
      double precision :: zsat = 5d0
      
      double precision :: dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
      double precision :: dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
      double precision :: dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
      double precision :: dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
      
      double precision, parameter :: w = 5.0d-5 ! m yr^-1, uplift rate
C       double precision, parameter :: w = 1.0d-4 ! m yr^-1, uplift rate
      
      double precision, parameter :: vcnst = 1.0d1 ! m yr^-1, advection
C       double precision, parameter :: qin = 5d-3 ! m yr^-1, advection (m3 water / m2 profile / yr)

      double precision :: qin = 1d-1 ! m yr^-1, advection (m3 water / m2 profile / yr)
      double precision v(nz), q
      
      double precision :: hr = 1d5 ! m^2 m^-3, reciprocal of hydraulic radius
C       double precision :: hr = 1d4 ! m^2 m^-3, reciprocal of hydraulic radius
      
      double precision po2(nz), redsld(nz), redaq(nz), ms(nz), c(nz)
      double precision po2x(nz), msx(nz), cx(nz)
      double precision msi
      double precision msili
      double precision msil(nz),msilx(nz)
      
      double precision c2(nz), c2x(nz),ctmp, po2tmp
      double precision so4(nz), so4x(nz)
      double precision na(nz), nax(nz), naeq(nz), silsat(nz) 
      double precision pro(nz), prox(nz), dumreal(nz), dprodna(nz)
      
      double precision koxa(nz), kdis(nz), koxs(nz),koxs2(nz)
      double precision koxai(nz), kdisi(nz), koxsi(nz),koxs2i(nz)
      double precision ksil(nz) 
      
      double precision kho, ucv
      double precision kco2,k1, keqsil, kw
      
      integer iz, ie, it, ie2
      
      integer col, row
      
      double precision imbr
      
      integer, parameter :: nmx = nsp*nz
      integer, parameter :: nmx2 = 1*nz
      integer, parameter :: nmx3 = nsp3*nz
      
      double precision amx(nmx,nmx),ymx(nmx)
      double precision amx2(nmx2,nmx2),ymx2(nmx2)
      double precision amx3(nmx3,nmx3),ymx3(nmx3)
      double precision emx3(nz)
!       double precision y2mx(3*(nz-1)),a2mx(3*(nz-1),3*(nz-1)) 
!       double precision, allocatable :: y2mx(:),a2mx(:,:) 
      
      integer ipiv(nmx)
      integer ipiv2(nz)
      integer ipiv3(nmx3)
!       integer, allocatable :: ipiv(:)
      integer info
      
      external DGESV
      
      double precision error, error2
      double precision :: tol = 1d-6
      
      integer :: spc = 50
      
C       integer, parameter :: nrec = 22
      integer, parameter :: nrec = 20
      integer reclis(nrec)
      double precision rectime(nrec)
      character(3) chr
      character(256) runname,workdir, chrz(3), chrq(3),base,fname
      integer irec, iter, idum, iter2
      
      double precision :: swad = 1.0d0! 1.0 when advection included 0.0d0 when not
      double precision dt2, swpe  ! physical erosion 
      
!       integer, parameter :: nt2 = int(nt*w/v)
!       double precision time2(nt2)
      
      integer zlis(3*(nz-1))
      integer nel, imx, imx2
      
      double precision :: po2th = 1.0d-20
      double precision minpo2
      
      double precision :: cth = 1.0d-20
      double precision :: c2th = 1.0d-20
      double precision :: so4th = 1.0d-20
      double precision :: proth = 1.0d-20
      double precision :: nath = 1.0d-20
      double precision :: msth = 1.0d-300
      double precision :: msilth = 1.0d-300
      
      double precision prepo2
      double precision :: stoxs = 15.0d0/4.0d0  ! 15/4 py => Fe-oxide + sulfate; 7/2 py => Fe++ + sulfate
      double precision :: stoxa = 1.0d0/4.0d0  ! stoichiomety of oxidation in aq
      double precision :: swoxa = 0.0d0   ! switch for oxidation in aq
      double precision :: swoxs2 = 1.0d0  ! switch for oxidation in solid phase
      
      double precision :: swbr = 0.0d0  ! switch for biological respiration
      double precision :: vmax = 0.71d0 ! mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
      double precision :: mo2 = 0.121d0 ! Michaelis, Davidson et al. (2012)
      
      double precision :: swex = 0.0d0 ! switch for explicit
      double precision :: frex = 0.0d0 ! fraction of explicit
      
      double precision :: swadvmass = 0d0 ! switch; 1 when calculating q from advection mass balance
      double precision :: waterfluc = 0d0 ! switch: 1 when fluctuating water flow
      
      double precision pyoxflx, feoxflx, respflx, diflx, advflx, o2tflx
      double precision feadvflx(2), feox2flx(2), fepy1flx(2)
     $  , fepy2flx(2), fetflx(2), fediflx(2)
      double precision pyadvflx, pyox1flx, pyox2flx, pytflx
      double precision siladvflx, sildisflx,  siltflx
      double precision naadvflx, nadisflx,  natflx, nadiflx
      double precision so4advflx, so4disflx,  so4tflx, so4diflx
      double precision o2flxsum,feflxsum(2),pyflxsum, so4flxsum
      double precision silflxsum,naflxsum
      double precision dummy, zdum(nz)
      
      double precision :: maxdt = 10d0
      
      integer izdum
C       logical :: pre_calc = .false.
      logical :: pre_calc = .true.
      
      logical :: read_data = .false.
C       logical :: read_data = .true.
      
      data rectime /1d1,3d1,1d2,3d2,1d3,3d3,1d4,3d4
     $   ,1d5,2d5,3d5,4d5,5d5,6d5,7d5,8d5,9d5,1d6,1.1d6,1.2d6/
C       data rectime /-1d6,0d6,1d6,2d6,3d6,4d6,5d6,6d6,7d6,8d6
C      $   ,9d6,10d6,11d6,12d6,13d6,14d6,15d6,16d6,17d6,18d6,19d6,20d6/
C       data rectime /21d6,22d6,23d6,24d6,25d6,26d6,27d6,28d6,29d6,30d6
C      $   ,31d6,32d6,33d6,34d6,35d6,36d6,37d6,38d6,39d6,40d6,41d6,42d6/
      logical :: dir_exist
      double precision o2out,zrxn,dms_max
      double precision :: beta = 3.7d19/0.21d0
      double precision :: o2in = 1d13
      double precision :: area = 1.5d14
      double precision pyrxn(nz)
      double precision :: pi = 4.0d0*datan(1d0)
	  
      logical :: O2_evolution = .false.
C       logical :: O2_evolution = .true.

      logical :: flgback = .false.
      double precision :: zab(3), zpy(3) 
      integer :: oxj
      !-------------------------
      
C       rectime =rectime/1d4
C       rectime =rectime/2d0

#if var1==1
      qin=10d0**(-3.0d0)
#elif var1==2
      qin=10d0**(-2.8d0)
#elif var1==3
      qin=10d0**(-2.6d0)
#elif var1==4
      qin=10d0**(-2.4d0)
#elif var1==5
      qin=10d0**(-2.2d0)
#elif var1==6
      qin=10d0**(-2.0d0)
#elif var1==7
      qin=10d0**(-1.8d0)
#elif var1==8
      qin=10d0**(-1.6d0)
#elif var1==9
      qin=10d0**(-1.4d0)
#elif var1==10
      qin=10d0**(-1.2d0)
#elif var1==11
      qin=10d0**(-1.0d0)
#elif var1==12
      qin=10d0**(-0.8d0)
#elif var1==13
      qin=10d0**(-0.6d0)
#elif var1==14
      qin=10d0**(-0.4d0)
#elif var1==15
      qin=10d0**(-0.2d0)
#elif var1==16
      qin=10d0**(0.0d0)
#endif

#if var2==1
      zsat=1d0
#elif var2==2
      zsat=5d0
#elif var2==3
      zsat=10d0
#elif var2==4
      zsat=15
#elif var2==5
      zsat=20d0
#elif var2==6
      zsat=25d0
#elif var2==7
      zsat=30d0
#elif var2==8
      zsat=35d0
#elif var2==9
      zsat=40d0
#elif var2==10
      zsat=45d0
#elif var2==11
      zsat=50d0
#endif

      write(chrq(1),'(i0)') int(qin/(10d0**(floor(log10(qin)))))
      write(chrq(2),'(i0)') floor(log10(qin))
      chrq(3) = trim(adjustl(chrq(1)))//'E'//trim(adjustl(chrq(2)))
      write(chrz(3),'(i0)') int(zsat)
      
      vmax = vmax * 1d0  !!  vmax is increased by a factor of 100 (cf., soil respiration in Liu and Zhou 2006)
      
       mo2 = mo2*po2i/0.21d0     !! mo2 is assumed to proportional to po2i
      
C       do irec = 1,nrec 
C         rectime(irec) = 19d6 + irec*1d6
C       enddo
      
      write(workdir,*) 'C:/cygwin64/home/YK/PyWeath/'     
      write(base,*) '_w1e-5_msx1_msil%100_S1e5_z300'     
C       write(workdir,*) write(workdir,*) '/home/latruffe/PyWeath/'         
      write(runname,*) 'Fe2+SO4+sil+ph_wet_iter'
     $   //'---q'//trim(adjustl(chrq(3)))//'_z'
     $   //trim(adjustl(chrz(3)))//trim(adjustl(base))
C      $   //'_co21e-1-o26e-2'
#ifdef pyweath
      write(runname,*) 'Fe2+SO4_wet_iter'
     $   //'---q'//trim(adjustl(chrq(3)))//'_z'
     $   //trim(adjustl(chrz(3)))//trim(adjustl(base))
C      $   //'_co21e-1-o26e-2'
#endif
#ifdef silweath
      write(runname,*) 'sil+ph_wet_iter'
     $   //'---q'//trim(adjustl(chrq(3)))//'_z'
     $   //trim(adjustl(chrz(3)))//trim(adjustl(base))
C      $   //'_co21e-1-o26e-2'
#endif
      
      call system ('mkdir -p '//trim(adjustl(runname)))
      
C       if (.not.read_data) call system ('rm '//trim(adjustl(workdir))//
C      $   trim(adjustl(runname))//
C      $    '/*')
      
      open(65, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(o2flx).txt', 
     $                              status='unknown', action = 'write')
      open(67, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(fe2flx).txt', 
     $                              status='unknown', action = 'write')
      open(68, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(fe3flx).txt', 
     $                              status='unknown', action = 'write')
      open(69, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(pyflx).txt', 
     $                              status='unknown', action = 'write')
      open(55, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(naflx).txt', 
     $                              status='unknown', action = 'write')
      open(56, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(silflx).txt', 
     $                              status='unknown', action = 'write')
      open(57, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(so4flx).txt', 
     $                              status='unknown', action = 'write')
      
      open(71, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(o2flx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(73, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(fe2flx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(74, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(fe3flx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(75, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(pyflx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(76, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(naflx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(77, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(silflx-alltime).txt', 
     $                              status='unknown', action = 'write')
      open(78, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(so4flx-alltime).txt', 
     $                              status='unknown', action = 'write')
     
     
      open(95, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2sense-res.txt', 
     $                              status='unknown', action = 'write')
      open(97, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2sense-res(alltime).txt', 
     $                              status='unknown', action = 'write')
      
      
      do iz=1,nz+1
        ze(iz)=(iz-1)*ztot/nz
      end do
      
      z(:) = 0.5d0*(ze(2:nz+1)+ze(1:nz))
      
C       print *,ze;stop
C       do iz=1,nz
C         z(iz)=(iz-1)*ztot/(nz-1)
C       end do
      
      dz=z(2)-z(1)
      
      stoxs = stoxs - swoxa*stoxa
            
C       print *, stoxs
            
      do irec=1,nrec
        reclis(irec)=irec*nt/nrec
      end do
            
      po2 = po2i
C       po2 = po2th
      sat = sati
C       sat = min(1.0d0,0.90d0*(z-ztot/2d0)/ztot/2d0 + 1.d0)
C       zsat = ztot/6.0d0
C       zsat = 5d0
      sat = min(1.0d0,(1d0-satup)*z/zsat + satup)
      v = vcnst
      v = qin/poroi/sat
      poro = poroi
      torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
      tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
      deff = torg*dgas + tora*daq
      c(1) = 0.0d0
      c(2:) = 0.0d0
      c2(1) = 0.0d0
      c2(2:) = 0.0d0
      so4 = 0d0
      na = 0d0
      pro = 10.0d0**(-ph)
      
      dt = maxdt
      dt = 1d-2
      
      
      if (swex == 1.0d0) then 
      
      dt = dz/vcnst
      dt2 = dz/w
      
!       do ie = 1, nt2
!         time2(ie) = dt2*ie
!       end do
      
      end if
      
      kho=10.0d0**(-2.89d0)*exp(13.2d0*(1.0d0/(273.0d0+tc)
     $  -1.0d0/(273.0d0+25.0d0))/rg)
     
      kco2 = 10.0d0**(-1.34d0)  ! 15C Kanzaki Murakami 2015
      k1 = 10.0d0**(-6.42d0)  ! 15C Kanzaki Murakami 2015
      kw = 10.0d0**(-14.35d0)  ! 15C Kanzaki Murakami 2015
      
      pro = sqrt(kco2*k1*pco2i+kw)
      
C       print *,pro(1);stop
      
      
      ksil = 1.31d-9*1d4  ! mol/m2/yr  ! from Li et al., 2014
      
      keqsil = 3.412182823d0 - 0.5d0* 8.310989613d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
      keqsil = 10.0d0**(keqsil)
      
C       print *,keqsil;stop
      
      ucv = 1.0d0/(rg2*(273.0d0+tc))
      
      redsld = redsldi*1d-2/120.0d0*23.94d0*2.7d0*(1.0d0-poro)  
      ! 120 g mol^-1 and 23.94 cm^3 mol^-1 for pyrite and 2.7 g cm^-3 for solid density
      ! redsld is fraction and unitless
      
      msi = redsldi*1d-2/120.0d0*2.7d0*(1.0d0-poroi)*1d6 
      
      if (O2_evolution) msi =  4d0*o2in/15d0/w/area     
      ! mol m^-3       
      ms = msi
      
      msili = silwti*1d-2/262.2d0*2.7d0*(1.0d0-poroi)*1d6   
      ! 262.2 g mol^-1 and 100.07 cm^3 mol^-1 (Robie et al., 1967) for albite and 2.7 g cm^-3 for solid density
      
      msil = msili
      
C       print *,msili;stop
      
      koxsi = 10.0d0**(-8.19d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (po2**0.5)
     $   *(kho)**(0.50d0)/((10.0d0**(-ph))**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)
     
      koxai = swoxa*8.0d13*60.0d0*24.0d0*365.0d0   !  excluding the term (c*po2)
     $   *(10.0d0**(-14.0d0+ph))**2.0d0               ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)
     
      koxs2i = swoxs2*10.0d0**(-6.07d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (fe3**0.93/fe2**0.40)
      !! mol m^-2 yr^-1, Williamson and Rimstidt (1994)
      
      ! when adv o2 flux == adv pyrite flx (q*kho*po2i*1d3 = stoxs*msi*w)
      if (swadvmass == 1d0) then 
      q = 15d0/4d0*msi*w/kho/po2i/1d3
      v = q/poroi/sat
      endif 
      ! *** the above must be commented out when using arbitrary q value
      
      
      open(22, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-bsd.txt', 
     $                              status='unknown', action = 'write')
      
      do iz = 1,nz
        write(22,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz)
      enddo
      
      close(22)
      
      
      !  --------- read -----
      if (read_data) then 
      open(255, file=trim(adjustl(workdir))//
     $  'po2-sense_z1800_qmbs_dt2e-2_wt30_ver6'//
     $                  '-go2_v2_cnt7'//'/'//
     $   'o2profile-res-018.txt', 
     $    status ='old',action='read')
      do iz = 1,nz
        read(255,*) zdum(iz),po2x(iz), cx(iz), msx(iz),c2x(iz), time
      enddo
      close(255)
      
      if (zdum(nz) > ztot) then  ! interpolating to finer scale
         izdum = 1
         do iz = 1,nz-1
           do while (z(izdum)<=zdum(iz+1))  
              po2(izdum) = (po2x(iz)-po2x(iz+1))/(zdum(iz)-zdum(iz+1))
     $                       *(z(izdum)-zdum(iz+1))+po2x(iz+1)
              c(izdum) = (cx(iz)-cx(iz+1))/(zdum(iz)-zdum(iz+1))
     $                       *(z(izdum)-zdum(iz+1))+cx(iz+1)
              ms(izdum) = (msx(iz)-msx(iz+1))/(zdum(iz)-zdum(iz+1))
     $                       *(z(izdum)-zdum(iz+1))+msx(iz+1)
              c2(izdum) = (c2x(iz)-c2x(iz+1))/(zdum(iz)-zdum(iz+1))
     $                       *(z(izdum)-zdum(iz+1))+c2x(iz+1)
              izdum = izdum+1
              if (izdum>=nz) then
                po2(izdum) = po2x(nz)
                c(izdum) = cx(nz)
                ms(izdum) = msx(nz)
                c2(izdum) = c2x(nz)
                exit
              endif
           enddo
           if (izdum>=nz) exit
         enddo
      endif
      
      if (zdum(nz) == ztot) then 
        po2 = po2x
        c = cx
        ms = msx
        c2 = c2x
      endif

#ifdef display      
      Print *,'==== printing read data (1) ===='
      print *, (po2(iz),iz=1,Nz,50)
      print *, (c(iz),iz=1,Nz,50)
      print *, (c2(iz),iz=1,Nz,50)
      print *, (ms(iz),iz=1,Nz,50)
      print *, ''
#endif      
C       c = cth
C       c2 = c2th
      
      open(255, file=trim(adjustl(workdir))//
     $  'po2-sense_z1800_qmbs_dt2e-2_wt30_ver6'//
     $                  '-go2_v2_cnt7'//'/'//
     $   'o2sense-res(alltime).txt', 
     $    status ='old',action='read')
C       do iz = 1,21
      dummy = 0d0
      do while (.not.dummy == time)
        read(255,*) dummy,o2in, o2out, po2i,msi,zrxn
        if (dummy==time) exit
      enddo
      close(255)
#ifdef display      
      print *,'**** printing read data (2) *****'
      print *, 'dummy,o2in, o2out, po2i,msi,zrxn'
      print *, dummy,o2in, o2out, po2i,msi,zrxn
      print *, '~~~ End reading ~~~'
#endif      
      endif
      
      
      ! --------- loop -----
      it = 0
C       if (.not.read_data)time = 0
      
C       do irec = 1,nrec-1
C         if ((rectime(irec)-time)*(rectime(irec+1)-time)<0d0) then 
C             idum=irec
C             exit
C         endif
C         idum = 0
C       enddo
      
      idum =0
      if (O2_evolution) then 
        if (.not.read_data) time = -1.1d6  
      endif
      
      irec = idum
          
          
      nel = 0
      zlis = 0
      
      
!! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@
      
 100  do while (it<nt)
#ifdef display 
      print *, it, time
#endif
      if (time>rectime(nrec)) exit
      
      if (waterfluc == 1d0) then
      
      if (swadvmass == 1d0) then 
      q = 15d0/4d0*msi*w/kho/po2i/1d3
      else 
      q = qin
      endif
      
      if (mod(it,2)==0) then 
      q = q + q*0.5d0
      else if (mod(it,2)==1) then 
      q = q - q*0.50d0
      endif
      
      v = q/poroi/sat
      
      endif 
      
      swpe = 0.0d0 ! physical erosion
      
!       if (any((time+dt-time2)*(time-time2)<0.0d0)) then
!         swpe = 1.0d0
!         print *, "===================="
!       end if
       
       ! -------- modifying dt --------
       
C        if ((iter <= 10).and.(dt<1d1)) then
       if (dt<maxdt) then
         dt = dt*1.01d0
         if (dt>maxdt) dt = maxdt
       endif
       if (iter > 300) then
         dt = dt/1.05d0
       end if

      ! ======== modifying maxdt ===============
C       if (time >1d4) maxdt = 10d0
      ! ========================================
#ifdef pHiter      
!  ############## pH interation #######################
       
      error2 = 1d4
      iter2 = 0
      do while (error2 > tol ) 
#endif
      
      
      koxsi = 10.0d0**(-8.19d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (po2**0.5)
     $   *(kho)**(0.50d0)/(pro**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)
     
      koxai = max(swoxa*8.0d13*60.0d0*24.0d0*365.0d0   !  excluding the term (c*po2)
     $   *(kw/pro)**2.0d0               ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)
     $  , swoxa*1d-7*60.0d0*24.0d0*365.0d0)
      
      koxs = koxsi
      koxa = koxai
      
      koxs2 = koxs2i
      
      ! forcing and oxygen calculation 
      if (O2_evolution) then 
      pyrxn = 0d0
      do iz=1,nz
      pyrxn(iz)=koxs(iz)*poro(iz)*hr*23.94d0*1d-6*ms(iz)*
     $  po2(iz)**0.50d0
     $  *merge(0.0d0,1.0d0,po2(iz)<po2th)
      enddo
      if (time < 0d0) then 
         o2in = 1d13 
         o2out = 1d13
         msi = 4d0*o2in/15d0/w/area
         po2i = 0.21d0
         mo2 = mo2*po2i/0.21d0     !! mo2 is assumed to proportional to po2i
         if (it == 0) then 
            zrxn = 0d0
         else 
            dms_max = 0d0
            do iz = 1,nz
               if (abs(pyrxn(iz))>dms_max) then 
                 zrxn = z(Iz)
                 dms_max = max(dms_max,abs(pyrxn(iz)))
               endif
            enddo 
         endif
      else        
      ! water table
C       zsat = 30d0 + 20d0*sin(2d0*pi*time/1d4)  ! fluctuation
      zsat = 50d0 ! step
C       zsat = 30d0 + 20d0*(time/20d6)  ! linear increase 
      sat = min(1.0d0,(1d0-satup)*z/zsat + satup)
      if (swadvmass == 1d0) then 
      q = 15d0/4d0*msi*w/kho/po2i/1d3
      v = q/poroi/sat
      endif 
      poro = poroi
      torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
      tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
      deff = torg*dgas + tora*daq
      
      
C         o2in = 2.0d13                 !  doubling
C         o2in = 1.5d13                 !  1.5 times
C         o2in = 0.5d13                 !  halving
        o2in = 1.0d13                 !  halving
C         o2in = 2d13/2d7*time+1d13   ! case if o2in is continuously changed
        if (.not.it==0) o2out = area*pyoxflx
C         msi = 4d0*o2in/15d0/w/area  ! case if msi is modified with o2in
        po2i = po2i + (o2in - o2out)*dt/beta 
        dms_max = 0d0
        mo2 = mo2*po2i/0.21d0     !! mo2 is assumed to proportional to po2i
        do iz = 1,nz-1
           if (abs(pyrxn(iz))>dms_max) then 
             zrxn = z(iz)
             dms_max = max(dms_max,abs(pyrxn(iz)))
           endif
        enddo        
      endif
      
      endif
      
      
C       po2(1) = po2i
C       ms(nz) = msi
C       msil(nz) = msili
      
      po2x = po2
      cx = c
      msx = ms
      
      c2x = c2
      
      so4x=so4
      
      nax = na
      msilx = msil
      
      prox = pro
      
      error = 1d4
      iter=0
      
!       if (swex /= 1.0d0) then
      
      if (pre_calc) then 
      
      do iz = 1, nz
      
      if (po2x(iz)>=po2th) cycle
      
      if (iz/=nz) po2tmp = po2(iz+1)
      if (iz==nz) po2tmp = po2(iz)
      
      if (iz/=1) then 
      po2x(iz) = max(0.0d0,-(dt/
     $ (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))*
     $ ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(-po2(iz))/dt
     $  -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2tmp
     $  +po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0)
     $   -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2(iz)-po2(iz-1))/(dz**2.0d0)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad
C      $  +poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*po2(iz)/dz
     $                        ))
C      $  +(poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))
C      $           *v(iz)*kho*1d3*(po2(iz))/dz
C      $  +poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*po2(iz)/dz
      else ! iz == 1
      po2x(iz) = max(0.0d0,-(dt/
     $ (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))*
     $ ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(-po2(iz))/dt
     $  -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(-2.0d0*po2(iz) + po2(iz+1)+po2i)
     $                  /(dz**2.0d0)
C      $   -1d3*(ucv*(poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(po2(iz)-po2i)/(dz**2.0d0)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*swad
C      $  +poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*po2(iz)/dz
     $                        ))
C      $  +(poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))
C      $           *v(iz)*kho*1d3*(po2(iz))/dz
C      $  +poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*po2(iz)/dz
      endif
      
      end do 
      
      do iz = 1, nz
      
      if (msx(iz)>=msth) cycle
      
      if (iz/=nz) then 
      msx(iz) = max(0d0,
     $    ms(iz) +
     $         dt*(
     $   w*(ms(iz+1)-ms(iz))/dz
     $             )
     $             )
      else  
      msx(iz) = max(0d0,
     $    ms(iz) +
     $         dt*(
     $   w*(msi-ms(iz))/dz
     $             )
     $             )
      endif
      
      enddo
      
      do iz = 1, nz
      
      if (msilx(iz)>=msilth) cycle
      
      if (iz/=nz) then 
      msilx(iz) = max(0d0,
     $    msil(iz) +
     $         dt*(
     $   w*(msil(iz+1)-msil(iz))/dz
     $             )
     $             )
      else 
      msilx(iz) = max(0d0,
     $    msil(iz) +
     $         dt*(
     $   w*(msili-msil(iz))/dz
     $             )
     $             )
      endif 
      
      enddo
      
      
      if (swoxa == 1d0) then 
      
      do iz = 1, nz
      
      if (c2x(iz)>=c2th) cycle
      
      if (iz/=nz) ctmp = c2(iz+1)
      if (iz==nz) ctmp = c2(iz)
      if (iz/=1) then 
      c2x(iz) = max(0.0d0,
     $    c2(iz) +
     $         dt*(
     $   -v(iz)*(c2(iz)-c2(iz-1))/dz
     $   +dfe3*tora(iz)*(ctmp
     $   +c2(iz-1)-2d0*c2(iz))/(dz**2d0)
     $   +dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2(iz)-c2(iz-1))/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +koxa(iz)*po2x(iz)*c(iz)
     $             )
     $             )
      else
      c2x(iz) = max(0.0d0,
     $    c2(iz) +
     $         dt*(
     $   -v(iz)*(c2(iz)-c2i)/dz
     $   +dfe3*tora(iz)*(ctmp
     $   +c2i-2d0*c2(iz))/(dz**2d0)
C      $   +dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c2(iz)-c2i)/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +koxa(iz)*po2x(iz)*c(iz)
     $             )
     $             )
      end if
!       print *, iz,c2x(iz)
      
      if (cx(iz)>=cth) cycle
      
      if (iz/=nz) ctmp = c(iz+1)
      if (iz==nz) ctmp = c(iz)
      if (iz/=1) then 
      cx(iz) = max(0.0d0,
     $    c(iz) +
     $         dt*(
     $   -v(iz)*(c(iz)-c(iz-1))/dz
     $   +dfe2*tora(iz)*(ctmp
     $   +c(iz-1)-2d0*c(iz))/(dz**2d0)
     $   +dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c(iz)-c(iz-1))/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +15d0*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*c(iz)**(-0.40d0)
     $    *1d-3
     $  +koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0
     $   *1d-3
     $             )
     $             )
      else 
      cx(iz) = max(0.0d0,
     $    c(iz) +
     $         dt*(
     $   -v(iz)*(c(iz)-ci)/dz
     $   +dfe2*tora(iz)*(ctmp
     $   +ci-2d0*c(iz))/(dz**2d0)
C      $   +dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c(iz)-ci)/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +15d0*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*c(iz)**(-0.40d0)
     $    *1d-3
     $  +koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0
     $   *1d-3
     $             )
     $             )
      endif
      
      if (so4x(iz)>=so4th) cycle
      
      if (iz/=nz) ctmp = so4(iz+1)
      if (iz==nz) ctmp = so4(iz)
      if (iz/=1) then 
      so4x(iz) = max(0.0d0,
     $    so4(iz) +
     $         dt*(
     $   -v(iz)*(so4(iz)-so4(iz-1))/dz
     $   +dso4*tora(iz)*(ctmp
     $   +so4(iz-1)-2d0*so4(iz))/(dz**2d0)
     $   +dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(so4(iz)-so4(iz-1))/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +2d0*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*c(iz)**(-0.40d0)
     $    *1d-3
     $  +2d0*koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0
     $   *1d-3
     $             )
     $             )
      else 
      so4x(iz) = max(0.0d0,
     $    so4(iz) +
     $         dt*(
     $   -v(iz)*(so4(iz)-so4i)/dz
     $   +dso4*tora(iz)*(ctmp
     $   +so4i-2d0*so4(iz))/(dz**2d0)
C      $   +dso4/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(so4(iz)-so4i)/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +2d0*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*c(iz)**(-0.40d0)
     $    *1d-3
     $  +2d0*koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0
     $   *1d-3
     $             )
     $             )
      endif 
      
      if (nax(iz)>=nath) cycle
      
      if (iz/=nz) ctmp = na(iz+1)
      if (iz==nz) ctmp = na(iz)
      if (iz/=1) then 
      nax(iz) = max(0.0d0,
     $    na(iz) +
     $         dt*(
     $   -v(iz)*(na(iz)-na(iz-1))/dz
     $   +dna*tora(iz)*(ctmp
     $   +na(iz-1)-2d0*na(iz))/(dz**2d0)
     $   +dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(na(iz)-na(iz-1))/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)*1d-3
     $             )
     $             )
      else 
      nax(iz) = max(0.0d0,
     $    na(iz) +
     $         dt*(
     $   -v(iz)*(na(iz)-nai)/dz
     $   +dna*tora(iz)*(ctmp
     $   +nai-2d0*na(iz))/(dz**2d0)
C      $   +dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(na(iz)-nai)/(dz**2d0)
C      $   -(v(iz)-v(iz-1))*c2(iz)/dz
     $  +ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)*1d-3
     $             )
     $             )
      endif 
      end do
      
      end if 
      
C       nax(2:nz) = naeq(2:nz)
      
      if (any(isnan(po2x)).or.any(isnan(cx)).or.any(isnan(c2x))) then 
        print*, 'error in precalc'
        stop
      endif
      
      end if
      
      if (O2_evolution) then 
      if (swoxa == 1d0) then 
        c = 0.1d0*cth
        c2 = 0.1d0*c2th
        cx = c
        c2x = c2
      endif
      endif
      
!       end if
#ifndef silweath      
      do while ((.not.isnan(error)).and.(error > tol))
      
      amx=0.0d0
      ymx=0.0d0
      
      if (it == 0 .and. iter == 0) then
        cx(:) = 1.0d2
        c2x(:) = 1.0d2
        so4x(:) = 1.0d2
        nax(:) = 1.0d2
      end if
      
C       print*,msilx
C       stop
      
      
      do iz = 1, nz  !================================
      
        row = nsp*(iz-1)+1
        
        if (iz/=nz) then
        
          amx(row,row) = (1.0d0/dt    !  z
     $      + w/dz *(1.0d0-swex)   
     $     + (1.0d0-frex)*
     $  koxs(iz)*poro(iz)*hr*23.94d0*1d-6
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $
     $    + (1.0d0-frex)*merge(0.0d0,
     $      + koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $   cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $                )
C      $      * msx(iz)
     $      * merge(1.0d0,msx(iz),msx(iz)<msth)
          
          amx(row,row+nsp) = (-w/dz)*(1.0d0-swex) 
C      $       *msx(iz+1)      ! z + dz
     $       *merge(1.0d0,msx(iz+1),msx(iz)<msth)
          
          ymx(row) = (
     $     (msx(iz)-ms(iz))/dt
     $      -w*(msx(iz+1)-msx(iz))/dz*(1.0d0-swex) 
     $      -w*(ms(iz+1)-ms(iz))/dz*swex/dt*dt2*swpe
     $     + (1.0d0-frex)*
     $      koxs(iz)*poro(iz)*hr
     $          *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $     + frex*
     $      koxs(iz)*poro(iz)*hr
     $          *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))
     $       
     $   + (1.0d0-frex)*merge(0.0d0,
     $   koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $      cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $   + frex*merge(0.0d0,
     $   koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $      c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))
     $                       )
     $  *merge(0.0d0,1d0,msx(iz)<msth)
     
        else if (iz==nz) then
        
          amx(row,row) = (1.0d0/dt 
     $      + w/dz*(1.0d0-swex)
     $    + (1.0d0-frex)*
     $     koxs(iz)*poro(iz)*hr*23.94d0*1d-6
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $         
     $     + (1.0d0-frex)*merge(0.0d0,
     $    koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $          cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $               )
C      $        *msx(iz)
     $        *merge(1.0d0,msx(iz),msx(iz)<msth)
          
          ymx(row) = (
     $     (msx(iz)-ms(iz))/dt
     $      -w*(msi-msx(iz))/dz*(1.0d0-swex)
     $      -w*(msi-ms(iz))/dz*swex*dt2/dt*swpe
     $      + (1.0d0-frex)*
     $      koxs(iz)*poro(iz)*hr
     $             *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $      + frex*
     $      koxs(iz)*poro(iz)*hr
     $             *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))
     $          
     $    + (1.0d0-frex)*merge(0.0d0,
     $      koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $          cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $    + frex*merge(0.0d0,
     $      koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $          c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))
     $                       )
     $  *merge(0.0d0,1d0,msx(iz)<msth)
        end if 
        
C         if (iz /= 1) then
        
          amx(row,row + 2 ) = (
     $    + (1.0d0-frex)*merge(0.0d0,
     $   koxs(iz)*poro(iz)*hr*23.94d0*1d-6*msx(iz)
     $     *0.50d0*(po2x(iz)**(-0.50d0)),
     $       po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))
     $               )
     $        *po2x(iz)
     $  *merge(0.0d0,1d0,msx(iz)<msth)
      
          amx(row,row  + 1 ) = (
     $    + (1.0d0-frex)*merge(0.0d0,
     $     koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *msx(iz)*c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-1.40d0),
     $        cx(iz)<cth
     $   .or.isnan(c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-1.40d0)))
     $                     )
     $        *cx(iz)
     $  *merge(0.0d0,1d0,msx(iz)<msth)
          
          amx(row,row  + 3 ) = (
     $    + (1.0d0-frex)*merge(0.0d0,
     $    koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $      *msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),
     $      (c2x(iz)<c2th).or.(cx(iz)<cth)
     $   .or.isnan(c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0)))
     $                     )
     $        *c2x(iz)
     $  *merge(0.0d0,1d0,msx(iz)<msth)
     
C         end if 
        
      end do  !================================
      
      do iz = 1, nz
        
        row = nsp*(iz-1)+2
        
        if (.not.((iz == 1).or.(iz==nz))) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(-2d0)
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $  + (1.0d0-frex)*koxa(iz)*po2x(iz)
     $  + (1.0d0-frex)*merge(0.0d0,
     $  -(15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),
     $     cx(iz)<cth
     $  .or.isnan(c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0)))
     $   *1d-3
     $                       )
C      $    *cx(iz)
     $    *merge(1.0d0,cx(iz),cx(iz)<cth)
     
          amx(row,row-nsp) = (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
     $      *cx(iz-1)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
          amx(row,row+nsp) = (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
     $      *cx(iz+1)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
          ymx(row) = (
     $  (cx(iz)-c(iz))/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(cx(iz)-cx(iz-1))/(dz**2d0))
     $   +swex*(-dfe2*tora(iz)*(c(iz+1)+c(iz-1)-2d0*c(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c(iz)-c(iz-1))/(dz**2d0))
     $  + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(c(iz)-c(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $  + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  + frex*koxa(iz)*c(iz)*po2(iz)
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
        else if (iz == 1) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(-2d0)
     $    /(dz**2d0)
C      $   -dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(1d0)/(dz**2d0)
     $     )
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $  + koxa(iz)*po2x(iz)*(1.0d0-frex)
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),
     $     cx(iz)<cth
     $  .or.isnan(c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0)))
     $    *1d-3
     $                       )
C      $   *cx(iz)
     $    *merge(1.0d0,cx(iz),cx(iz)<c2th)
     
          amx(row,row+nsp) = (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
     $      *cx(iz+1)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
          ymx(row) = (
     $  (cx(iz)-c(iz))/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))
     $    /(dz**2d0)
C      $   -dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(cx(iz)-ci)/(dz**2d0)
     $     )
     $   +swex*(-dfe2*tora(iz)*(c(iz+1)+ci-2d0*c(iz))
     $    /(dz**2d0)
C      $   -dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c(iz)-ci)/(dz**2d0)
     $    )
     $  + v(iz)*(cx(iz)-ci)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*cx(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c(iz)-ci)/dz*swex
C      $  + (v(iz)-v(iz-1))*c(iz)/dz*swex
     $  + koxa(iz)*cx(iz)*po2x(iz)*(1.0d0-frex)
     $  + koxa(iz)*c(iz)*po2(iz)*frex
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
        
        else if (iz == nz) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(-1d0)
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $  + (1.0d0-frex)*koxa(iz)*po2x(iz)
     $  + (1.0d0-frex)*merge(0.0d0,
     $  -(15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),
     $     cx(iz)<cth)*1d-3
     $                       )
C      $    *cx(iz)
     $    *merge(1.0d0,cx(iz),cx(iz)<cth)
     
          amx(row,row-nsp) = (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
     $      *cx(iz-1)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
          ymx(row) = (
     $  (cx(iz)-c(iz))/dt 
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(cx(iz)-cx(iz-1))/(dz**2d0))
     $   +swex*(-dfe2*tora(iz)*(c(iz-1)-1d0*c(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c(iz)-c(iz-1))/(dz**2d0))
     $  + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(c(iz)-c(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $  + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  + frex*koxa(iz)*c(iz)*po2(iz)
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $     ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
        end if 
        
!         if (row == 188) write (*,*) row, it, iz, ymx(row)
        
        amx(row,row+1) = (
     $   + (1.0d0-frex)*merge(0.0d0,koxa(iz)*cx(iz),po2x(iz)<po2th)
     $   - (1.0d0-frex)*merge(0.0d0,
     $  koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*0.50d0
     $    *po2x(iz)**(-0.50d0),
     $   po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))*1d-3 
     $                    )
     $    *po2x(iz)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
      
        amx(row,row+2) = (
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $   (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),
     $    (c2x(iz)<c2th).or.(cx(iz)<cth)
     $  .or.isnan(c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0)))*1d-3
     $                    )
     $    *c2x(iz)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
!         if (it == 64 .and. iter == 1 .and. row == 137) then
!           print *, amx(row,row+1),po2x(iz),po2x(iz)**(-0.50d0),
!      $     po2x(iz)**(0.50d0)
!         end if
     
C         if (iz /= nz) then
        
          amx(row,row  - 1) = (     
     $   - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6
     $   *merge(0d0,po2x(iz)**(0.50d0),
     $    po2x(iz) <po2th.or. isnan(po2x(iz)**(0.50d0)))*1d-3 
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $     cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $                            )
     $        *msx(iz)
     $    *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)
     
C         end if 
        
      end do  ! ==============================
      
      do iz = 1, nz
        
        row = nsp*(iz-1)+4  !! fe+++
        
        if (.not.((iz == 1).or.(iz==nz))) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(-2d0)
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz *(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),
     $   (c2x(iz)<c2th).or.(cx(iz)<cth)
     $  .or.isnan(c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0)))*1d-3
     $                       )
C      $    *c2x(iz)
     $    *merge(1.0d0,c2x(iz),c2x(iz)<c2th)
     
          amx(row,row - nsp) = (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - v(iz)/dz*(1.0d0-swex)
     $                       )
C      $      *c2x(iz-1)
     $    *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)
     
          amx(row,row+nsp) = (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(1d0)
     $    /(dz**2d0))
     $            )
C      $      *c2x(iz+1)
     $    *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)

     
          ymx(row) = (
     $  (c2x(iz)-c2(iz))/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2x(iz)-c2x(iz-1))/(dz**2d0))
     $   +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2(iz-1)-2d0*c2(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2(iz)-c2(iz-1))/(dz**2d0))
     $  + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2(iz-1))/dz*swex
C      $  + (v(iz)-v(iz-1))*c2(iz)/dz*swex
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)
     
        else if (iz == 1) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(-2d0)
     $    /(dz**2d0)
C      $   -dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(1d0)/(dz**2d0)
     $    )
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $  +(1.0d0-frex)*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),
     $     (c2x(iz)<c2th).or.(cx(iz)<cth)
     $ .or.isnan(c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0)))*1d-3
     $                       )
C      $   *c2x(iz)
     $    *merge(1.0d0,c2x(iz),c2x(iz)<c2th)
     
          amx(row,row+nsp) = (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(1d0)
     $    /(dz**2d0))
     $            )
C      $      *c2x(iz+1)
     $    *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)
     
          ymx(row) = (
     $  (c2x(iz)-c2(iz))/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))
     $    /(dz**2d0)
C      $   -dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c2x(iz)-c2i)/(dz**2d0)
     $   )
     $   +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2i-2d0*c2(iz))
     $    /(dz**2d0)
C      $   -dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c2(iz)-c2i)/(dz**2d0)
     $   )
     $  + v(iz)*(c2x(iz)-c2i)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2i)/dz*swex
C      $  + (v(iz)-v(iz-1))*c2(iz)/dz*swex
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $  +(1.0d0-frex)*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)  ! commented out 
        
        else if (iz ==nz) then
        
          amx(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(-1d0)
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz *(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),
     $   (c2x(iz)<c2th).or.(cx(iz)<cth)
     $  .or.isnan(c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0)))
     $    *1d-3
     $                       )
C      $    *c2x(iz)
     $    *merge(1.0d0,c2x(iz),c2x(iz)<c2th)
     
          amx(row,row - nsp) = (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - v(iz)/dz*(1.0d0-swex)
     $                       )
C      $      *c2x(iz-1)
     $    *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)

     
          ymx(row) = (
     $  (c2x(iz)-c2(iz))/dt 
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2x(iz)-c2x(iz-1))/(dz**2d0))
     $   +swex*(-dfe3*tora(iz)*(c2(iz-1)-1d0*c2(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2(iz)-c2(iz-1))/(dz**2d0))
     $  + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2(iz-1))/dz*swex
C      $  + (v(iz)-v(iz-1))*c2(iz)/dz*swex
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)
        
        end if 
        
!         if (row == 188) write (*,*) row, it, iz, ymx(row)
        
        amx(row,row-1) = (
     $   - (1.0d0-frex)*koxa(iz)*cx(iz)
     $                    )
C      $    *po2x(iz)
     $    *merge(0.0d0,po2x(iz),c2x(iz)<c2th)
     
        amx(row,row-2) = (
     $   - (1.0d0-frex)*koxa(iz)*po2x(iz)
     $  +(1.0d0-frex)*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $       c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),
     $    cx(iz)<cth.or.
     $  isnan(c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0)))*1d-3
     $                    )
C      $    *cx(iz)
     $    *merge(0.0d0,cx(iz),c2x(iz)<c2th)
     
!         if (it == 64 .and. iter == 1 .and. row == 137) then
!           print *, amx(row,row+1),po2x(iz),po2x(iz)**(-0.50d0),
!      $     po2x(iz)**(0.50d0)
!         end if
     
C         if (iz /= nz) then
        
          amx(row,row  - 3) = (
     $   +(1.0d0-frex)*merge(0.0d0,
     $    (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $      cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $                            )
C      $        *msx(iz)
     $    *merge(0.0d0,msx(iz),c2x(iz)<c2th)
     
C         end if 
      
      end do
      
      do iz = 1, nz
        
        row = 3 + nsp*(iz-1)
        
        if (iz == 1) then
        
          amx(row,row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt
     $ +(1.0d0-frex)*
     $  2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
C      $   -(1.0d0-frex)*1d3*(ucv*
C      $   (poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(1.0d0)/(dz**2.0d0)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz*(1.0d0-swex)
C      $   +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $  po2x(iz)<po2th)
     $  +(1.0d0-frex)*merge(0.d0,
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $   po2x(iz)<po2th.or.isnan(msx(iz)*0.50d0*po2x(iz)**(-0.50d0)))
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
     
!      $  +
!      $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
!      $  
!      $  +
!      $   stoxs*koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0)
!      $   
!      $  +
!      $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0
!      $   
     $                    )
!      $        *po2x(iz)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $        *merge(1.0d0,po2x(iz),po2x(iz)<po2th)
     
      if (isnan(amx(row,row))) then 
        print*,'nan in oxygen',iz,
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt
     $ +(1.0d0-frex)*
     $  2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
C      $   -(1.0d0-frex)*1d3*(ucv*
C      $   (poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(1.0d0)/(dz**2.0d0)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz*(1.0d0-swex)
C      $   +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $  po2x(iz)<po2th)
     $  +(1.0d0-frex)*merge(0.d0,
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $   po2x(iz)<po2th.or.isnan(msx(iz)*0.50d0*po2x(iz)**(-0.50d0)))
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
      endif
     
          amx(row,row+nsp) = (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $                    )
!      $       *po2x(iz+1)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $        *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)
     
       if (isnan(amx(row,row+nsp))) then 
         print *,'error in oxygen +',iz,
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
       endif
       
          ymx(row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0)
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz+1)+po2i-2.0d0*po2(iz))/(dz**2.0d0)
C      $   -(1.0d0-frex)*1d3*(ucv*
C      $  (poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(po2x(iz)-po2i)/(dz**2.0d0)
C      $   -(frex)*1d3*(ucv*(poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(po2(iz)-po2i)/(dz**2.0d0)
     $ +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz
     $     *(1.0d0-swex)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*(swex)
C      $  +(1.0d0-swex)*
C      $ (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2x(iz))/dz
C      $  +(swex)*
C      $ (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $    *(po2(iz))/dz
C      $  +(1.0d0-swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
C      $  +(swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
!      $  +merge(0.0d0,
!      $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz),
!      $  po2x(iz)<po2th)
!      $  +merge(0.0d0,
!      $  stoxs*koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*msx(iz)*po2x(iz)**(0.50d0),
!      $  po2x(iz) < po2th)
!      $  +merge(0.0d0,
!      $  swbr*vmax*po2x(iz)/(po2x(iz)+mo2),
!      $   po2x(iz) < po2th)
     
     $  +(1.0d0-frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)
     $  +(frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $  
     $  +(1.0d0-frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $     ,(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0))))
     $  +(frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $     ,(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0))))
     $  
     $  +(1.0d0-frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $  +(frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $   
     $                    )
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $         *merge(0.0d0,1.0d0,po2x(iz)<po2th)
        
        else if (iz == nz) then
        
          amx(row,row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt
     $  +(1.0d0-frex)*
     $ 1.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(1.0d0)/(dz**2.0d0)
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3*(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $  po2x(iz) <po2th)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $   po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
     
!      $  +
!      $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
!      $  
!      $  +
!      $  stoxs*koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0)
!      $   
!      $  +
!      $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0
!      $  
     $                    )
!      $        *po2x(iz)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $        *merge(1.0d0,po2x(iz),po2x(iz)<po2th)
      
      if (isnan(amx(row,row))) then 
        print*,'nan in oxygen', iz,
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt,
     $  +(1.0d0-frex)*
     $ 1.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0),
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(1.0d0)/(dz**2.0d0),
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz,
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3*(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $  po2x(iz) <po2th),
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $   po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))),
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
      endif
      
          amx(row,row-nsp) = (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(-1.0d0)/(dz**2.0d0)
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz
     $                    )
!      $     *po2x(iz-1)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $     *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)
     
      if (isnan(amx(row,row-nsp))) then 
        print*,'nan in oxygen - ',iz,
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0),
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(-1.0d0)/(dz**2.0d0),
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz
       endif
     
          ymx(row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0)
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz-1)-1.0d0*po2(iz))/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2x(iz)-po2x(iz-1))/(dz**2.0d0)
     $   -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2(iz)-po2(iz-1))/(dz**2.0d0)
     $  +(1.0d0-swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz
     $  +(swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $     *(po2x(iz))/dz
C      $  +(swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $     *(po2(iz))/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
C      $  +(swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
!      $  +merge(0.0d0,
!      $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz),
!      $  po2x(iz) < po2th)
!      $  +merge(0.0d0,
!      $  stoxs*koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*msx(iz)*po2x(iz)**(0.50d0),
!      $  po2x(iz) < po2th)
!      $  +merge(0.0d0,
!      $  swbr*vmax*po2x(iz)/(po2x(iz)+mo2),
!      $  po2x(iz) < po2th)
     
     $  +(1.0d0-frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
     $   *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz)))
     $  
     $  +(1.0d0-frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6
     $   *merge(0d0,msx(iz)*po2x(iz)**(0.50d0),po2x(iz) <po2th
     $        .or.isnan(msx(iz)*po2x(iz)**(0.50d0)))
     $  
     $  +(1.0d0-frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $
     $  +(frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $  
     $  +(frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6
     $   *merge(0d0,ms(iz)*po2(iz)**(0.50d0),po2x(iz) <po2th
     $       .or.isnan(ms(iz)*po2(iz)**(0.50d0)))
     $  
     $  +(frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $  
     $                    )
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)
     $         *merge(0.0d0,1.0d0,po2x(iz)<po2th)
        
        else
        
          amx(row,row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt
     $  +(1.0d0-frex)*
     $  2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(1.0d0)/(dz**2.0d0)
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $     *(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $   po2x(iz) <po2th)
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $    po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
     
!      $  +
!      $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
!      $   
!      $  +
!      $  stoxs*koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0)
!      $    
!      $  +
!      $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0
!      $   
     $                    )
!      $           *po2x(iz)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $           *merge(1.0d0,po2x(iz),po2x(iz)<po2th)
     
      if (isnan(amx(row,row))) then 
        print*,'nan in oxygen ',iz,
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt,
     $  +(1.0d0-frex)*
     $  2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $  +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0),
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(1.0d0)/(dz**2.0d0),
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz,
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $     *(1.0d0)/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(1d0)/dz
     $  +(1.0d0-frex)*merge(0.0d0,
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),
     $   po2x(iz) <po2th),
     $  +(1.0d0-frex)*merge(0.0d0,
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0),
     $    po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))),
     $  +(1.0d0-frex)*merge(0.0d0,
     $   swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0,
     $   (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0)))
      endif
      
      
          amx(row,row+nsp) = (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $                    )
!      $     *po2x(iz+1)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $     *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)
      
      if (isnan(amx(row,row+nsp))) then 
        print*,'nan in oxygen +',iz,
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
      endif
      
!             if (abs(amx(row,row+3))>1d50) then
!               print *, amx(row,row+3),po2x(iz+1),po2x(iz)
!             end if
     
          amx(row,row-nsp) = (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(-1.0d0)/(dz**2.0d0)
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz
     $                    )
!      $       *po2x(iz-1)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)
     $       *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)
     
      if (isnan(amx(row,row-nsp))) then 
        print*,'nan in oxygen -',iz,
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0),
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(-1.0d0)/(dz**2.0d0),
     $  +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz
       endif
     
          ymx(row) = (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2x(iz)-po2x(iz-1))/(dz**2.0d0)
     $  +(1.0d0-swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz*swad
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $    *(po2x(iz))/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
     $
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz+1)+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0)
     $   -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2(iz)-po2(iz-1))/(dz**2.0d0)
     $  +(swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad
C      $  +(swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2(iz))/dz
C      $  +(swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
!      $  +merge(0.0d0,
!      $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz),
!      $   po2x(iz) < po2th)
!      $  +merge(0.0d0,
!      $   stoxs*koxs(iz)*poro(iz)*hr*sat(iz)
!      $      *23.94d0*1d-6*msx(iz)*po2x(iz)**(0.50d0),
!      $   po2x(iz) < po2th)
!      $  +merge(0.0d0,
!      $   swbr*vmax*po2x(iz)/(po2x(iz)+mo2),
!      $    po2x(iz) < po2th)
     
     $  +(1.0d0-frex)*
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
     &     *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz)))
     $  +(1.0d0-frex)*
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th
     $    .or.isnan(po2x(iz)**(0.50d0)))
     $   
     $  +(1.0d0-frex)*
     $   swbr*vmax
     $   *merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $  +(frex)*
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $  +(frex)*
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**(0.50d0)))
     $   
     $  +(frex)*
     $   swbr*vmax
     $   *merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $                    )
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $         *merge(0.0d0,1.0d0,po2x(iz)<po2th)
        
        end if 
        
!         if (row == 144) write (*,*) row, it, iz, ymx(row),amx(row,row)
        
        amx(row,row-1) = (
     $     +(1.0d0-frex)*poro(iz)*sat(iz)*1d3*koxa(iz)*po2x(iz)*stoxa              
     $                    )
!      $   *cx(iz)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $   *merge(0.0d0,cx(iz),po2x(iz)<po2th)
        
        if (isnan(amx(row,row-1))) then 
          print*,'nan in oxygen for fe2',iz,
     $     +(1.0d0-frex)*poro(iz)*sat(iz)*1d3*koxa(iz)*po2x(iz)*stoxa 
        endif
        
        
C         if (iz /= nz) then
          
          amx(row,row  - 2) = (
!      $      +merge(0.0d0,
!      $   koxs(iz)*poro(iz)*hr
!      $      *23.94d0*1d-6*po2x(iz)**(0.50d0),
!      $     po2x(iz)< po2th)
     
     $   (1.0d0-frex)*koxs(iz)*poro(iz)*hr*stoxs
     $      *23.94d0*1d-6
     $   *merge(0d0,po2x(iz)**(0.50d0),
     $   po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0)))
     $                            )
!      $       *msx(iz)
!      $         /merge(po2x(iz),1.0d0,po2x(iz)<po2th)

     $       *merge(0.0d0,msx(iz),po2x(iz)<po2th)
        
         if (isnan(amx(row,row - 2))) then 
           print*,'nan in oxygen for pyrite',iz,
     $   (1.0d0-frex)*koxs(iz)*poro(iz)*hr*stoxs
     $      *23.94d0*1d-6
     $   *merge(0d0,po2x(iz)**(0.50d0),
     $   po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0)))
        endif
        
C         end if 
        
      end do 
      
      ymx=-1.0d0*ymx
      
!       allocate (y2mx(3*(nz-1)-nel))
!       allocate (a2mx(3*(nz-1)-nel,3*(nz-1)-nel))
!       
!       allocate (ipiv(3*(nz-1)-nel))
      
!       write (*,*) "after allocation",nel,3*(nz-1)-nel
      
!       do ie = 1,3*(nz-1)
!         if (zlis(ie) == 1) write (*,*) ie
!       end do

!       y2mx = 0.0d0
!       a2mx = 0.0d0
!       
!       if (nel /= 0) then
!       
!       imx = 1
!       
!       do ie = 1, 3*(nz-1)
!       
!         if (zlis(ie) == 0) then
!           
!           y2mx(imx) = ymx(ie)
!           
!           imx2 = 1
!           
!           do ie2 = 1, 3*(nz-1)
!             
!             if (zlis(ie2) == 0) then 
!             
!               a2mx(imx,imx2) = amx(ie,ie2)
!               
!               imx2 = imx2 + 1
!             
!             end if
!           
!           end do
!           
!           imx = imx + 1
!           
!         end if
!       
!       end do
      
!       print *, "imx, imx2",imx, imx2
      
!       else if (nel == 0) then
!         y2mx = ymx
!         a2mx = amx
!       end if
      
!       if (it == 3 .and. iter == 17) then
!       
!       open (29,file='test.txt',status='replace')
!       do iz = 1,3*(Nz-1)
!       write (29,*) (amx(iz,ie),ie=1,3*(Nz-1))
!       end do
!       close(29)
!       
!       open (343,file='test2.txt',status='replace')
!       do iz = 1,3*(Nz-1)
!       write (343,*) (ymx(iz))
!       end do
!       close(343)
!       
!       end if
      
      if (any(isnan(amx)).or.any(isnan(ymx))) then 
        print*,'error in mtx'
        print*,'any(isnan(amx)),any(isnan(ymx))'
        print*,any(isnan(amx)),any(isnan(ymx))
        
        if (any(isnan(ymx))) then 
          do ie = 1,nmx
            if (isnan(ymx(ie))) then 
              print*,'NAN is here...',ie,mod(ie,300)
            endif
          enddo
        endif
        
        
        if (any(isnan(amx))) then 
          do ie = 1,nmx
            do ie2 = 1,nmx
              if (isnan(amx(ie,ie2))) then 
                print*,'NAN is here...',ie,mod(ie,nz),ie2,mod(ie2,nz)
              endif
            enddo
          enddo
        endif
        
        stop
      endif
      
      call DGESV(nmx,int(1),amx,
     $           nmx,IPIV,ymx,nmx,INFO) 
      
!       call DGESV(3*(Nz-1)-nel,int(1),a2mx,
!      $           3*(Nz-1)-nel,IPIV,y2mx,3*(Nz-1)-nel,INFO) 
      
      
!       if (it == 3 .and. iter == 17) then
!       
!       open (53,file='test-res.txt',status='replace')
!       do iz = 1,3*(Nz-1)
!       write (53,*) (ymx(iz))
!       end do
!       close(53)
!       
!       end if
      
!       ymx = 0.0d0
!       imx = 1
!       do ie = 1, 3*(nz-1)
!         if (zlis(ie) == 0) then
!           ymx(ie) = y2mx(imx)
!           imx = imx + 1
!         end if 
!       end do 
      
!       deallocate  (ipiv)
!       deallocate  (a2mx)
!       deallocate  (y2mx)
      
      if (any(isnan(ymx))) then
        print*,'error in soultion'
      endif
      
      do iz = 1, nz
        row = 1 + nsp*(iz-1)
        
!         if (msx(iz) /= 0.0d0) then
!         y2mx(row) = ymx(row)
!         else if (msx(iz) == 0.0d0) then
!         y2mx(row) = 0.0d0
!         end if 
        
        if (isnan(ymx(row))) then 
          print *,'nan at', iz,z(iz),zrxn,'pyrite'
          if (z(iz)<zrxn) then 
            ymx(row)=0d0
C             msx(iz) = msx(iz)*0.5d0
          endif
        endif
        
        if ((.not.isnan(ymx(row))).and.ymx(row) >10d0) then 
          msx(iz) = msx(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
          msx(iz) = msx(iz)*0.50d0
        else   
        
        
        msx(iz) = msx(iz)*exp(ymx(row))
        
        endif
        
!         y2mx(row) = abs(ymx(row)/msx(iz))
!         msx(iz) = msx(iz) + ymx(row)
!         
!         if (msx(iz) < 0.0d0) msx(iz) = msx(iz) - ymx(row)
        
      end do
      
      do iz = 1, nz
        row = 2 + nsp*(iz-1)
        
!         if (cx(iz) /= 0.0d0) then
!         y2mx(row) = ymx(row)
!         else if (cx(iz) == 0.0d0) then
!         y2mx(row) = 0.0d0
!         end if 
        
!         if (po2x(iz) /= 0.0d0) then
!         y2mx(row+1) = ymx(row+1)
!         else if (po2x(iz) == 0.0d0) then
!         y2mx(row+1) = 0.0d0
!         end if 
        if (isnan(ymx(row))) then 
          print *,'nan at', iz,z(iz),zrxn,'Fe2'
          if (z(iz)<zrxn) then 
            ymx(row)=0d0
            cx(iz) = 0.1d0*cth
          endif
        endif
        if (isnan(ymx(row+1))) then 
          print *,'nan at', iz,z(iz),zrxn,'oxygen'
          if (z(iz)<zrxn) then 
            ymx(row+1)=0d0
            po2x(iz) = 0.1d0*po2th
          endif
        endif
        if (isnan(ymx(row+2))) then 
          print *,'nan at', iz,z(iz),zrxn,'Fe3'
          if (z(iz)<zrxn) then 
            ymx(row+2)=0d0
            c2x(iz) = 0.1d0*c2th
          endif
        endif
        
        if ((.not.isnan(ymx(row))).and.ymx(row) >10d0) then 
          cx(iz) = cx(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
          cx(iz) = cx(iz)*0.50d0
        else
        cx(iz) = cx(iz)*exp(ymx(row))
        endif
        
        if ((.not.isnan(ymx(row+1))).and.ymx(row+1) >10d0) then 
          po2x(iz) = po2x(iz)*1.5d0
        else if (ymx(row+1) < -10d0) then 
          po2x(iz) = po2x(iz)*0.50d0
        else
        po2x(iz) = po2x(iz)*exp(ymx(row+1))
        endif
        if ((.not.isnan(ymx(row+2))).and.ymx(row+2) >10d0) then 
          c2x(iz) = c2x(iz)*1.5d0
        else if (ymx(row+2) < -10d0) then 
          c2x(iz) = c2x(iz)*0.50d0
        else
        c2x(iz) = c2x(iz)*exp(ymx(row+2))
        endif

        if (swoxa==0d0) then  ! ignoring Fe2 and Fe3
          ymx(row)=0d0
          ymx(row+2) = 0d0
        endif
        
!           y2mx(row) = abs(ymx(row)/cx(iz))
!           y2mx(row+1) = abs(ymx(row+1)/po2x(iz))
!           
!           cx(iz) = cx(iz) + ymx(row)
!           po2x(iz) = po2x(iz) + ymx(row+1)
!           
!           if (po2x(iz) < 0.0d0) po2x(iz) = po2x(iz) - ymx(row+1)
!           if (cx(iz) < 0.0d0) cx(iz) = cx(iz) - ymx(row + 1)
        
!         if (po2x(iz) < 1.0d-50) then
!           po2x(iz) = 1.0d-50
!           ymx(row+1) = 0.0d0
!         end if 
      end do 
      
!       error = maxval(y2mx)

      error = maxval(exp(abs(ymx))) - 1.0d0
      
      if (isnan(error).or.info/=0 
     $  .or. any(isnan(cx)) .or. any(isnan(c2x))
     $  .or. any(isnan(po2x)) .or. any(isnan(msx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to'
     &  //' those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(cx)), any(isnan(c2x))'
     $    //', any(isnan(po2x)),any(isnan(msx)))'
        print*,isnan(error),info/=0,any(isnan(cx)), any(isnan(c2x))
     $     , any(isnan(po2x)),any(isnan(msx))
        stop
        cx = c
        po2x =po2
        msx = ms
        c2x = c2
C         dt = dt/1.05d0
        iter = iter + 1
        cycle
      endif

!       print *, (po2x(iz),iz=1,nz, spc)
!       print *, (cx(iz),iz=1,nz, spc)
!       print *, (msx(iz),iz=1,nz, spc)
!       print *, (c2x(iz),iz=1,nz, spc)
      
      do iz = 1, nz
        row = 1 + nsp*(iz-1)
        
        if (msx(iz) < 0.0d0) then
!           msx(iz) = msx(iz) - ymx(row)
          msx(iz) = msx(iz)/exp(ymx(row))*0.5d0
          error = 1.0d0
        end if
      end do
      
!       koxs = koxsi
!       koxa = koxai
      
!       zlis = 0
!       nel = 0
      
      do iz = 1, nz
        row = 2 + nsp*(iz-1)
          
          if (po2x(iz) < 0.0d0) then
!             po2x(iz) = po2x(iz) - ymx(row+1)
            po2x(iz) = po2x(iz)/exp(ymx(row+1))*0.5d0
            error = 1.0d0
          end if
          if (cx(iz) < 0.0d0) then
!             cx(iz) = cx(iz) - ymx(row)
            cx(iz) = cx(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
          end if 
       
!         if (po2x(iz) < po2th) then
!           koxa(iz) = 0.0d0
!           koxs(iz) = 0.0d0
!           po2x(iz) = 0.0d0
!           write (*,*) "zero po2", row + 1
!           zlis(row+1) = 1
!           nel = nel + 1
!           error = 1.0d0
!         end if 

      end do 
      
#ifdef display      
      print *, 'error',error,info,nel
#endif      
      iter = iter + 1
      
!       if (iter>=50) then
!         dt = dt * 0.99d0
!       end if 

       if (iter > 300) then
         dt = dt/1.01d0
         if (dt==0d0) stop
C          go to 100 
       end if
      
      end do
      
      !!  so4 calculation start 
      
C       error = 1d4
C       iter=0
C       do while ((.not.isnan(error)).and.(error > tol))
      
      amx2=0.0d0
      ymx2=0.0d0
      
      do iz = 1, nz
        
        row = iz
        
        if (.not.((iz == 1).or.(iz==nz))) then
        
          amx2(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dso4*tora(iz)*(-2d0)
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $                       )
C      $    *cx(iz)
C      $    *merge(1.0d0,so4x(iz),so4x(iz)<so4th)
     
          amx2(row,row-1) = (
     $   +(1d0-swex)*(-dso4*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
C      $      *so4x(iz-1)
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
          amx2(row,row+1) = (
     $   +(1d0-swex)*(-dso4*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
C      $      *so4x(iz+1)
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
          ymx2(row) = (
     $  (-so4(iz))/dt 
     $   +swex*(-dso4*tora(iz)*(so4(iz+1)+so4(iz-1)-2d0*so4(iz))
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(so4(iz)-so4(iz-1))/(dz**2d0))
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*so4x(iz)/dz
     $  + swex*v(iz)*(so4(iz)-so4(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
        else if (iz == 1) then
        
          amx2(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dso4*tora(iz)*(-2d0)
     $    /(dz**2d0)
C      $   -dso4/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(1d0)/(dz**2d0)
     $  )
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $                       )
C      $   *so4x(iz)
C      $    *merge(1.0d0,so4x(iz),so4x(iz)<so4th)
     
          amx2(row,row+1) = (
     $   +(1d0-swex)*(-dso4*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
C      $      *so4x(iz+1)
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
          ymx2(row) = (
     $  (-so4(iz))/dt 
     $   +swex*(-dso4*tora(iz)*(so4(iz+1)+so4i-2d0*so4(iz))
     $    /(dz**2d0)
C      $   -dso4/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(so4(iz)-so4i)/(dz**2d0)
     $   )
C      $  + (v(iz)-v(iz-1))*so4x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(so4(iz)-so4i)/dz*swex
C      $  + (v(iz)-v(iz-1))*so4(iz)/dz*swex
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
        
        else if (iz == nz) then
        
          amx2(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dso4*tora(iz)*(-1d0)
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $                       )
C      $    *cx(iz)
C      $    *merge(1.0d0,cx(iz),cx(iz)<cth)
     
          amx2(row,row-1) = (
     $   +(1d0-swex)*(-dso4*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
C      $      *so4x(iz-1)
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
          ymx2(row) = (
     $  (-so4(iz))/dt 
     $   +swex*(-dso4*tora(iz)*(so4(iz-1)-1d0*so4(iz))
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(so4(iz)-so4(iz-1))/(dz**2d0))
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*so4x(iz)/dz
     $  + swex*v(iz)*(so4(iz)-so4(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*so4(iz)/dz
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $     ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )
C      $    *merge(0.0d0,1.0d0,so4x(iz)<so4th)   ! commented out (is this necessary?)
     
        end if 
        
      end do  ! ==============================
      
      ymx2=-1.0d0*ymx2
      
      if (any(isnan(amx2)).or.any(isnan(ymx2))) then 
        print*,'error in mtx'
        print*,'any(isnan(amx2)),any(isnan(ymx2))'
        print*,any(isnan(amx2)),any(isnan(ymx2))
        
        if (any(isnan(ymx2))) then 
          do ie = 1,(nz)
            if (isnan(ymx2(ie))) then 
              print*,'NAN is here...',ie
            endif
          enddo
        endif
        
        
        if (any(isnan(amx2))) then 
          do ie = 1,(nz)
            do ie2 = 1,(nz)
              if (isnan(amx2(ie,ie2))) then 
                print*,'NAN is here...',ie,ie2
              endif
            enddo
          enddo
        endif
        
        stop
      endif
      
      call DGESV((Nz),int(1),amx2,
     $           (Nz),IPIV2,ymx2,(Nz),INFO) 
      
      if (any(isnan(ymx2))) then
        print*,'error in soultion'
      endif
      
      do iz = 1, nz
        row = iz
        
        if (isnan(ymx2(row))) then 
          print *,'nan at', iz,z(iz),zrxn,'so4'
          if (z(iz)<zrxn) then 
            ymx2(row)=0d0
            so4x(iz) = 0.1d0*so4th
          endif
        endif
        
        if ((.not.isnan(ymx2(row))).and.(ymx2(row)>=0d0)) then 
C           if (so4x(iz)/=0d0) emx2(row) = abs(ymx2(row)/so4x(iz)) 
          so4x(iz) = ymx2(row)
        else
          stop
          so4x(iz) = 0d0
        endif

        if (O2_evolution .and. swoxa==0d0) then  ! ignoring so4
          ymx2(row)=0d0
        endif
      end do 
      
!       print *, (po2x(iz),iz=1,nz, spc)
!       print *, (cx(iz),iz=1,nz, spc)
!       print *, (msx(iz),iz=1,nz, spc)
!       print *, (c2x(iz),iz=1,nz, spc)
C        print *, (so4x(iz),iz=1,nz, spc)
       
      !!!  =-=-=-=-=-=-=-=-=- END of SO4 calculation  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif       
#ifndef pyweath       
      !!! =-=-=-=-=-=-=-=-=- START of calculation for Na and albite  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      
      error = 1d4
      iter = 0
      
C       prox = 1d-1

#ifdef silweath
      so4x = so4th
      cx = cth 
      c2x = c2th 
      
      if (it == 0 .and. iter == 0) nax(1:) = 1.0d2
      
#endif      
      
      if (it ==0) then 
      prox(:) = 0.5d0* (
     $   -1d0*(2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) 
     $ + sqrt((2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0 
     $ + 4d0*kco2*k1*pco2i)
     $                          )
      else 
      
      prox(:) = 0.5d0* (
     $   -1d0*(nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) 
     $ + sqrt((nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0 
     $ + 4d0*kco2*k1*pco2i)
     $                          )
      
      endif
      
      do while ((.not.isnan(error)).and.(error > tol))
      
      amx3=0.0d0
      ymx3=0.0d0
      
C       if (iter/=0) then 
      
      prox(:) = 0.5d0* (
     $   -1d0*(nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) 
     $ + sqrt((nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0 
     $ + 4d0*kco2*k1*pco2i)
     $                          )
      
      dprodna(:) = 0.5d0* (
     $   -1d0 
     $ + 0.5d0/sqrt(
     $  (nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0 
     $ + 4d0*kco2*k1*pco2i)*2d0
     $  *(nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))
     $                          )
      
C       endif
      
      do iz = 1, nz  !================================
      
        row = nsp3*(iz-1)+1
        
        if (iz/=nz) then
        
          amx3(row,row) = (1.0d0/dt    !  z
     $      + w/dz *(1.0d0-swex)   
     $     + (1.0d0-frex)*
     $  ksil(iz)*poro(iz)*hr*100.07d0*1d-6
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $                )
C      $      * msx(iz)
     $      * merge(1.0d0,msilx(iz),msilx(iz)<msth)
C      $      *1d-10
          
          amx3(row,row+nsp3) = (-w/dz)*(1.0d0-swex) 
C      $       *msx(iz+1)      ! z + dz
     $       *merge(1.0d0,msilx(iz+1),msilx(iz)<msilth)
C      $      *1d-10
          
          ymx3(row) = (
     $     (msilx(iz)-msil(iz))/dt
     $      -w*(msilx(iz+1)-msilx(iz))/dz*(1.0d0-swex) 
     $      -w*(msil(iz+1)-msil(iz))/dz*swex/dt*dt2*swpe
     $     + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $     + frex*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $                       )
     $  *merge(0.0d0,1d0,msilx(iz)<msilth)
C      $      *1d-10
     
        else if (iz==nz) then
        
          amx3(row,row) = (1.0d0/dt 
     $      + w/dz*(1.0d0-swex)
     $    + (1.0d0-frex)*
     $     ksil(iz)*poro(iz)*hr*100.07d0*1d-6
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $               )
C      $        *msx(iz)
     $        *merge(1.0d0,msilx(iz),msilx(iz)<msilth)
C      $      *1d-10
          
          ymx3(row) = (
     $     (msilx(iz)-msil(iz))/dt
     $      -w*(msili-msilx(iz))/dz*(1.0d0-swex)
     $      -w*(msili-msil(iz))/dz*swex*dt2/dt*swpe
     $      + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $      + frex*
     $      ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $                       )
     $  *merge(0.0d0,1d0,msilx(iz)<msilth)
C      $      *1d-10
        end if 
        
C         if (iz /= 1) then
      
          amx3(row,row + 1 ) = (
     $    + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msilx(iz)
     $   *(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $    + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msilx(iz)
     $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $                     )
     $        *nax(iz)
     $  *merge(0.0d0,1d0,msilx(iz)<msilth)
C      $      *1d-10
      
C           amx3(row,row-nsp3 + 2 ) = (
C      $    + (1.0d0-frex)*
C      $      ksil(iz)*poro(iz)*hr
C      $          *100.07d0*1d-6*msilx(iz)
C      $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)
C      $   *merge(0d0,1d0
C      $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
C      $                     )
C      $        *prox(iz)
C      $  *merge(0.0d0,1d0,msilx(iz)<msilth)
C      $      *1d-6
     
C         end if 
        
      end do  !================================
      
      do iz = 1, nz
        
        row = nsp3*(iz-1)+2
        
        if (.not.((iz == 1).or.(iz==nz))) then
        
          amx3(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(-2d0)
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $  + (1.0d0-frex)*(
     $  -(1.0d0)*ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           )
     $   *1d-3
     $  + (1.0d0-frex)*(
     $  -(1.0d0)*ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           )
     $   *1d-3
     $                       )
C      $    *cx(iz)
     $    *merge(1.0d0,nax(iz),nax(iz)<nath)
     
          amx3(row,row-nsp3) = (
     $   +(1d0-swex)*(-dna*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
     $      *nax(iz-1)
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
          amx3(row,row+nsp3) = (
     $   +(1d0-swex)*(-dna*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
     $      *nax(iz+1)
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
          ymx3(row) = (
     $  (nax(iz)-na(iz))/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nax(iz-1)-2d0*nax(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(nax(iz)-nax(iz-1))/(dz**2d0))
     $   +swex*(-dna*tora(iz)*(na(iz+1)+na(iz-1)-2d0*na(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(na(iz)-na(iz-1))/(dz**2d0))
     $  + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(na(iz)-na(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
        else if (iz == 1) then
        
          amx3(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(-2d0)
     $    /(dz**2d0)
C      $   -dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(1d0)/(dz**2d0)
     $   )
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))/dz*(1.0d0-swex)
     $  -(1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  -(1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )
C      $   *cx(iz)
     $    *merge(1.0d0,nax(iz),nax(iz)<nath)
     
          amx3(row,row+nsp3) = (
     $   +(1d0-swex)*(-dna*tora(iz)*(1d0)
     $    /(dz**2d0))
     $                )
     $      *nax(iz+1)
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
          ymx3(row) = (
     $  (nax(iz)-na(iz))/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nai-2d0*nax(iz))
     $    /(dz**2d0)
C      $   -dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(nax(iz)-nai)/(dz**2d0)
     $    )
     $   +swex*(-dna*tora(iz)*(na(iz+1)+nai-2d0*na(iz))
     $    /(dz**2d0)
C      $   -dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(na(iz)-nai)/(dz**2d0)
     $   )
     $  + v(iz)*(nax(iz)-nai)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*cx(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(na(iz)-nai)/dz*swex
C      $  + (v(iz)-v(iz-1))*c(iz)/dz*swex
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
        
        else if (iz == nz) then
        
          amx3(row,row) = (
     $  1.0d0/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(-1d0)
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(1d0)/(dz**2d0))
     $   + v(iz)/dz*(1.0d0-swex)
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))/dz
     $  + (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  + (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )
C      $    *cx(iz)
     $    *merge(1.0d0,nax(iz),nax(iz)<nath)
     
          amx3(row,row-nsp3) = (
     $   +(1d0-swex)*(-dna*tora(iz)*(1d0)
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(-1d0)/(dz**2d0))
     $  - (1.0d0-swex)*v(iz)/dz
     $                       )
     $      *nax(iz-1)
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
          ymx3(row) = (
     $  (nax(iz)-na(iz))/dt 
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz-1)-1d0*nax(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(nax(iz)-nax(iz-1))/(dz**2d0))
     $   +swex*(-dna*tora(iz)*(na(iz-1)-1d0*na(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(na(iz)-na(iz-1))/(dz**2d0))
     $  + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(na(iz)-na(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
        end if 
        
C           amx3(row,row + 1) = (     
C      $   - (1.0d0-frex)*
C      $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
C      $   *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)
C      $   *merge(0d0,1d0
C      $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
C      $           *1d-3 
C      $                            )
C      $        *prox(iz)
C      $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
C         if (iz /= nz) then
        
          amx3(row,row  - 1) = (     
     $   - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*1d0
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3 
     $                            )
     $        *msilx(iz)
     $    *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
     
C         end if 
        
      end do  ! ==============================
      
      ymx3=-1.0d0*ymx3
      
      if (any(isnan(amx3)).or.any(isnan(ymx3))) then 
        print*,'error in mtx'
        print*,'any(isnan(amx3)),any(isnan(ymx3))'
        print*,any(isnan(amx3)),any(isnan(ymx3))
        
        if (any(isnan(ymx3))) then 
          do ie = 1,nsp3*(nz)
            if (isnan(ymx3(ie))) then 
              print*,'NAN is here...',ie
            endif
          enddo
        endif
        
        
        if (any(isnan(amx3))) then 
          do ie = 1,nsp3*(nz)
            do ie2 = 1,nsp3*(nz)
              if (isnan(amx3(ie,ie2))) then 
                print*,'NAN is here...',ie,ie2
              endif
            enddo
          enddo
        endif
        
        stop
      endif
      
      call DGESV(nsp3*(Nz),int(1),amx3,
     $           nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 
     
     
      if (any(isnan(ymx3))) then
        print*,'error in soultion'
      endif
      
      do iz = 1, nz
        row = 1 + nsp3*(iz-1)
        
        if (isnan(ymx3(row))) then 
          print *,'nan at', iz,z(iz),'albite'
          if (z(iz)<zrxn) then 
            ymx3(row)=0d0
C             msx(iz) = msx(iz)*0.5d0
          endif
        endif
        
        if ((.not.isnan(ymx3(row))).and.ymx3(row) >10d0) then 
          msilx(iz) = msilx(iz)*1.5d0
        else if (ymx3(row) < -10d0) then 
          msilx(iz) = msilx(iz)*0.50d0
        else   
        
        
        msilx(iz) = msilx(iz)*exp(ymx3(row))
        
        endif
        
      end do
      
      do iz = 1, nz
        row = 2 + nsp3*(iz-1)
        
        if (isnan(ymx3(row))) then 
          print *,'nan at', iz,z(iz),'Na'
          if (z(iz)<zrxn) then 
            ymx3(row)=0d0
            nax(iz) = 0.1d0*nath
          endif
        endif
        
        if ((.not.isnan(ymx3(row))).and.ymx3(row) >1d0) then 
          nax(iz) = nax(iz)*1.5d0
        else if (ymx3(row) < -1d0) then 
          nax(iz) = nax(iz)*0.50d0
        else
        nax(iz) = nax(iz)*exp(ymx3(row))
        endif
        
C         if (isnan(ymx3(row+1))) then 
C           print *,'nan at', iz,z(iz),'pH'
C           if (z(iz)<zrxn) then 
C             ymx3(row+1)=0d0
C             prox(iz) = 0.1d0*proth
C           endif
C         endif
        
C         if ((.not.isnan(ymx3(row+1))).and.ymx3(row+1) >10d0) then 
C           prox(iz) = prox(iz)*1.5d0
C         else if (ymx3(row+1) < -10d0) then 
C           prox(iz) = prox(iz)*0.50d0
C         else
C         prox(iz) = prox(iz)*exp(ymx3(row+1))
C         endif
        
C         if (O2_evolution) then 
C         if (swoxa==0d0) then  ! ignoring Na
C           ymx3(row)=0d0
C         endif
C         endif
      end do 

      error = maxval(exp(abs(ymx3))) - 1.0d0
      
      if (isnan(error).or.info/=0 
     $  .or. any(isnan(nax)) .or. any(isnan(msilx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to'
     &  //' those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(nax))'
     &   //',any(isnan(msilx)))'
        print*,isnan(error),info/=0,any(isnan(nax)),any(isnan(msilx))
        stop
        nax = na
        msilx = msil
        prox = pro
C         dt = dt/1.05d0
        iter = iter + 1
        cycle
      endif
      
      prox(:) = 0.5d0* (
     $   -1d0*(nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) 
     $ + sqrt((nax(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0 
     $ + 4d0*kco2*k1*pco2i)
     $                          )
      
      naeq = (keqsil *prox/4d0)**(1d0/3d0)
      silsat = nax/naeq

!       print *, (po2x(iz),iz=1,nz, spc)
!       print *, (cx(iz),iz=1,nz, spc)
!       print *, (msx(iz),iz=1,nz, spc)
!       print *, (c2x(iz),iz=1,nz, spc)
!       print *, (so4x(iz),iz=1,nz, spc)
C        print *, (nax(iz),iz=1,nz, spc)
C        print *, (msilx(iz),iz=1,nz, spc)
C        print *, (prox(iz),iz=1,nz, spc)
      
      do iz = 1, nz
        row = 1 + nsp3*(iz-1)
        
        if (msilx(iz) < 0.0d0) then
!           msx(iz) = msx(iz) - ymx3(row)
          msilx(iz) = msilx(iz)/exp(ymx3(row))*0.5d0
          error = 1.0d0
        end if
      end do
      
      do iz = 1, nz
        row = 2 + nsp3*(iz-1)
        
          if (nax(iz) < 0.0d0) then
!             cx(iz) = cx(iz) - ymx3(row)
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
          end if 
        
C           if (prox(iz) < 0.0d0) then
C !             cx(iz) = cx(iz) - ymx3(row)
C             prox(iz) = prox(iz)/exp(ymx3(row+1))*0.5d0
C             error = 1.0d0
C           end if 

      end do 
      
#ifdef display      
      print *, 'error',error,info,nel
#endif      
      iter = iter + 1
      
!       if (iter>=50) then
!         dt = dt * 0.99d0
!       end if 

       if (iter > 300) then
         dt = dt/1.01d0
         if (dt==0d0) stop
C          go to 100
         flgback = .true.
       end if
      
      
      enddo
      
      if (flgback) then 
        flgback = .false. 
        go to 100
      endif     
      
#ifdef pHiter
      error2 = maxval( abs(pro-prox)/pro)
#ifdef display
      print*,"=========pH iteration========"
      print*,"iter2,error2",iter2,error2
      print*,"============================="
#endif
      iter2 = iter2 + 1
      
      pro = prox
C       pro = sqrt(prox*pro)
      
       if (iter2 > 300) then
C          dt = dt/1.01d0
         dt = dt/10d0
         if (dt==0d0) stop
C          go to 100
C          flgback = .true.
       end if
      
      enddo     
! ######################## end of pH iteration #######################
#endif

#endif  
      !  endif of ifndef pyweath
#ifdef display
       print *,'o2:', (po2x(iz),iz=1,nz, spc)
       print *,'fe2:', (cx(iz),iz=1,nz, spc)
       print *,'py:', (msx(iz),iz=1,nz, spc)
       print *, 'fe3:', (c2x(iz),iz=1,nz, spc)
       print *, 'so4:', (so4x(iz),iz=1,nz, spc)
       print *, 'na:', (nax(iz),iz=1,nz, spc)
       print *, 'sil:', (msilx(iz),iz=1,nz, spc)
       print *, 'ph:', (prox(iz),iz=1,nz, spc)
#endif      
      pyoxflx =0d0 
      feoxflx = 0d0
      respflx = 0d0
      diflx = 0d0 
      advflx = 0d0 
      o2tflx = 0d0
      
      feadvflx = 0d0
      feox2flx = 0d0
      fepy1flx = 0d0
      fepy2flx = 0d0
      fetflx = 0d0
      fediflx = 0d0
      
      pyadvflx = 0d0
      pyox1flx = 0d0
      pyox2flx = 0d0
      pytflx = 0d0
      
      siladvflx = 0d0 
      sildisflx = 0d0  
      siltflx = 0d0
      
      naadvflx = 0d0
      nadisflx = 0d0  
      natflx = 0d0 
      nadiflx = 0d0
      
      so4advflx = 0d0
      so4disflx = 0d0  
      so4tflx = 0d0 
      so4diflx = 0d0
      
      do iz = 1, nz
        
        if (iz == 1) then
          
          o2tflx = o2tflx +(
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $   )*dz
          
          diflx = diflx  + (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0)
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz+1)+po2i-2.0d0*po2(iz))/(dz**2.0d0)
C      $   -(1.0d0-frex)*1d3*(ucv*
C      $  (poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(po2x(iz)-po2i)/(dz**2.0d0)
C      $   -(frex)*1d3*(ucv*(poro(iz+1)*(1.0d0-sat(iz+1))*torg(iz+1)*dgas
C      $      -poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas)
C      $    +(poro(iz+1)*sat(iz+1)*kho*tora(iz+1)*daq
C      $             -poro(iz)*sat(iz)*kho*tora(iz)*daq))
C      $       *(po2(iz)-po2i)/(dz**2.0d0)
     $              )*dz
          
          advflx = advflx +(
     $ +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz
     $     *(1.0d0-swex)
     $  +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*(swex)
C      $  +(1.0d0-swex)*
C      $ (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2x(iz))/dz
C      $  +(swex)*
C      $ (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2(iz))/dz
C      $  +(1.0d0-swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
C      $  +(swex)*
C      $ poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
     $            ) *dz
     
          feoxflx = feoxflx +(
     $  +(1.0d0-frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)
     $  +(frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $       )*dz
     
          pyoxflx = pyoxflx  + (
     $  +(1.0d0-frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $     ,(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0))))
     $  +(frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $     ,(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0))))
     $    )*dz
     
          respflx = respflx  + (
     $  +(1.0d0-frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $  +(frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $    )*dz
     
        
        else if (iz == nz) then
          
          o2tflx = o2tflx  + (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $    )*dz
      
          diflx = diflx  + (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0)
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz-1)-1.0d0*po2(iz))/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2x(iz)-po2x(iz-1))/(dz**2.0d0)
     $   -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2(iz)-po2(iz-1))/(dz**2.0d0)
     $     )*dz
     
          advflx = advflx  + (
     $  +(1.0d0-swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz
     $  +(swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2x(iz))/dz
C      $  +(swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2(iz))/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
C      $  +(swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
     $    )*dz
      
          feoxflx = feoxflx  + (
     $  +(1.0d0-frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
     $   *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz)))
     $  +(frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $          )*dz
     
          pyoxflx = pyoxflx  + (
     $  +(1.0d0-frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6
     $   *merge(0d0,msx(iz)*po2x(iz)**(0.50d0),po2x(iz) <po2th
     $        .or.isnan(msx(iz)*po2x(iz)**(0.50d0)))
     $  +(frex)*
     $  stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6
     $   *merge(0d0,ms(iz)*po2(iz)**(0.50d0),po2x(iz) <po2th
     $       .or.isnan(ms(iz)*po2(iz)**(0.50d0)))
     $      )*dz
     
          respflx = respflx  + (
     $  +(frex)*
     $  stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $  +(frex)*
     $  swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $     )*dz
          
          
        
        else
      
          o2tflx = o2tflx +  (
     $   (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)
     $       *(po2x(iz)-po2(iz))/dt
     $               )*dz
     
          diflx = diflx   + (
     $  -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0)
     $   -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2x(iz)-po2x(iz-1))/(dz**2.0d0)
     $  -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas
     $    +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)
     $       *(po2(iz+1)+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0)
     $   -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas
     $      -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)
     $    +(poro(iz)*sat(iz)*kho*tora(iz)*daq
     $             -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))
     $       *(po2(iz)-po2(iz-1))/(dz**2.0d0)
     $             )*dz
     
          advflx = advflx  + (
     $  +(1.0d0-swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz*swad
C      $  +(1.0d0-swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2x(iz))/dz
C      $  +(1.0d0-swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2x(iz))/dz
     $  +(swex)*
     $  poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad
C      $  +(swex)*
C      $  (poro(iz)*sat(iz)-poro(iz-1)*sat(iz-1))*v(iz)*kho*1d3
C      $   *(po2(iz))/dz
C      $  +(swex)*
C      $  poro(iz)*sat(iz)*(v(iz)-v(iz-1))*kho*1d3*(po2(iz))/dz
     $              )*dz
      
          feoxflx = feoxflx  + (
     $  +(1.0d0-frex)*
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)
     &     *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz)))
     $  +(frex)*
     $   stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)
     $            )*dz
     
          pyoxflx = pyoxflx  + (
     $  +(1.0d0-frex)*
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th
     $    .or.isnan(po2x(iz)**(0.50d0)))
     $  +(frex)*
     $   stoxs*koxs(iz)*poro(iz)*hr
     $      *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**(0.50d0)
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**(0.50d0)))
     $               )*dz
     
          respflx = respflx  +(
     $  +(1.0d0-frex)*
     $   swbr*vmax
     $   *merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $  +(frex)*
     $   swbr*vmax
     $   *merge(0d0,po2x(iz)/(po2x(iz)+mo2),
     $         (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))
     $           )*dz
        
        end if 
        
      end do 
      
      o2flxsum = o2tflx + advflx + diflx + pyoxflx + feoxflx + respflx 
      
      do iz = 1, nz  !================================
        
        if (iz/=nz) then
        
          pytflx = pytflx + (
     $     (msx(iz)-ms(iz))/dt
     $           )*dz
     
          pyadvflx = pyadvflx + (
     $      -w*(msx(iz+1)-msx(iz))/dz*(1.0d0-swex) 
     $      -w*(ms(iz+1)-ms(iz))/dz*swex/dt*dt2*swpe
     $              )*dz
          
          pyox1flx = pyox1flx + (
     $     + (1.0d0-frex)*
     $      koxs(iz)*poro(iz)*hr
     $          *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $     + frex*
     $      koxs(iz)*poro(iz)*hr
     $          *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))
     $               )*dz
     
          pyox2flx = pyox2flx + (
     $   + (1.0d0-frex)*merge(0.0d0,
     $   koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $      cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $   + frex*merge(0.0d0,
     $   koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $      c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))
     $              )*dz
     
        else if (iz==nz) then
        
          pytflx = pytflx + (
     $     (msx(iz)-ms(iz))/dt
     $              )*dz
          
         pyadvflx = pyadvflx + ( 
     $      -w*(msi-msx(iz))/dz*(1.0d0-swex)
     $      -w*(msi-ms(iz))/dz*swex*dt2/dt*swpe
     &                 )*dz
          
          pyox1flx = pyox1flx + (
     $      + (1.0d0-frex)*
     $      koxs(iz)*poro(iz)*hr
     $             *23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))
     $      + frex*
     $      koxs(iz)*poro(iz)*hr
     $             *23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))
     $              )*dz
     
         pyox2flx = pyox2flx + (
     $    + (1.0d0-frex)*merge(0.0d0,
     $      koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $          cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))
     $    + frex*merge(0.0d0,
     $      koxs2(iz)*poro(iz)*hr*23.94d0*1d-6
     $         *ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $          c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))
     $                   ) *dz
        end if 
        
      end do  !================================
      
      pyflxsum = pytflx + pyadvflx + pyox1flx + pyox2flx 
      
      do iz = 1, nz
        
        if (.not.((iz == 1).or.(iz==nz))) then
          
          fetflx(1) = fetflx(1) + (
     $  (cx(iz)-c(iz))/dt 
     $           )*dz*poro(iz)*sat(iz)*1d3
      
          fediflx(1) = fediflx(1) + (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(cx(iz)-cx(iz-1))/(dz**2d0))
     $   +swex*(-dfe2*tora(iz)*(c(iz+1)+c(iz-1)-2d0*c(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c(iz)-c(iz-1))/(dz**2d0))
     $               )*dz*poro(iz)*sat(iz)*1d3
          
          feadvflx(1) = feadvflx(1) + (
     $  + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(c(iz)-c(iz-1))/dz
     $          )*dz*poro(iz)*sat(iz)*1d3
          
          feox2flx(1) = feox2flx(1) + (
     $  + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  + frex*koxa(iz)*c(iz)*po2(iz)
     $              )*dz*poro(iz)*sat(iz)*1d3
     
          fepy1flx(1) = fepy1flx(1) + (
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $                 )*dz*poro(iz)*sat(iz)*1d3
     
          fepy2flx(1) = fepy2flx(1) + (
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                   )*dz*poro(iz)*sat(iz)*1d3
     
        else if (iz == 1) then
        
          fetflx(1) = fetflx(1) + (
     $  (cx(iz)-c(iz))/dt 
     $          )*dz*poro(iz)*sat(iz)*1d3
      
          fediflx(1) = fediflx(1) + (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))
     $    /(dz**2d0)
C      $   -dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(cx(iz)-ci)/(dz**2d0)
     $     )
     $   +swex*(-dfe2*tora(iz)*(c(iz+1)+ci-2d0*c(iz))
     $    /(dz**2d0)
C      $   -dfe2/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c(iz)-ci)/(dz**2d0)
     $    )
     $               )*dz*poro(iz)*sat(iz)*1d3
          
          feadvflx(1) = feadvflx(1) + (
     $  + v(iz)*(cx(iz)-ci)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*cx(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c(iz)-ci)/dz*swex
     $             )*dz*poro(iz)*sat(iz)*1d3
     
          feox2flx(1) = feox2flx(1) + (
     $  + koxa(iz)*cx(iz)*po2x(iz)*(1.0d0-frex)
     $  + koxa(iz)*c(iz)*po2(iz)*frex
     $              )*dz*poro(iz)*sat(iz)*1d3
     
          fepy1flx(1) = fepy1flx(1) + (
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $              )*dz*poro(iz)*sat(iz)*1d3
     
          fepy2flx(1) = fepy2flx(1) + (
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $              )*dz*poro(iz)*sat(iz)*1d3
        
        else if (iz ==nz ) then
          
          fetflx(1) = fetflx(1) + (
     $  (cx(iz)-c(iz))/dt 
     $           )*dz*poro(iz)*sat(iz)*1d3
      
          fediflx(1) = fediflx(1) + (
     $   +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(cx(iz)-cx(iz-1))/(dz**2d0))
     $   +swex*(-dfe2*tora(iz)*(c(iz-1)-1d0*c(iz))
     $    /(dz**2d0)
     $   -dfe2/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c(iz)-c(iz-1))/(dz**2d0))
     $               )*dz*poro(iz)*sat(iz)*1d3
          
          feadvflx(1) = feadvflx(1) + (
     $  + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(c(iz)-c(iz-1))/dz
     $          )*dz*poro(iz)*sat(iz)*1d3
          
          feox2flx(1) = feox2flx(1) + (
     $  + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  + frex*koxa(iz)*c(iz)*po2(iz)
     $              )*dz*poro(iz)*sat(iz)*1d3
     
          fepy1flx(1) = fepy1flx(1) + (
     $  - (1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $     ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $                 )*dz*poro(iz)*sat(iz)*1d3
     
          fepy2flx(1) = fepy2flx(1) + (
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (15.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                   )*dz*poro(iz)*sat(iz)*1d3
     
        
        end if 
        
      end do  ! ==============================
      
      feflxsum(1) = fetflx(1) + fediflx(1) +  feadvflx(1) + feox2flx(1)
     $   + fepy1flx(1) + fepy2flx(1)
      
      do iz = 1, nz
        
        if (.not.((iz == 1).or.(iz==nz))) then
        
          fetflx(2) = fetflx(2) + (
     $  (c2x(iz)-c2(iz))/dt 
     $  )*dz*poro(iz)*sat(iz)*1d3
      
          fediflx(2) = fediflx(2) + (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2x(iz)-c2x(iz-1))/(dz**2d0))
     $   +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2(iz-1)-2d0*c2(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2(iz)-c2(iz-1))/(dz**2d0))
     $             )*dz*poro(iz)*sat(Iz)*1d3
     
          feadvflx(2) = feadvflx(2) + (
     $  + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2(iz-1))/dz*swex
     $               )*dz*poro(iz)*sat(Iz)*1d3
     
          feox2flx(2) = feox2flx(2) + (
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $                )*dz*poro(iz)*sat(iz)*1d3
     
          fepy2flx(2) = fepy2flx(2) + (
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)
     
        else if (iz == 1) then
          
          fetflx(2) = fetflx(2) + (
     $  (c2x(iz)-c2(iz))/dt 
     &           )*dz*poro(iz)*sat(Iz)*1d3
      
          fediflx(2) = fediflx(2) + (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))
     $    /(dz**2d0)
C      $   -dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c2x(iz)-c2i)/(dz**2d0)
     $   )
     $   +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2i-2d0*c2(iz))
     $    /(dz**2d0)
C      $   -dfe3/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(c2(iz)-c2i)/(dz**2d0)
     $   )
     $             )*dz*poro(iz)*sat(Iz)*1d3
     
          feadvflx(2) = feadvflx(2) + (
     $  + v(iz)*(c2x(iz)-c2i)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2i)/dz*swex
     $               )*dz*poro(iz)*sat(iz)*1d3
      
          feox2flx(2) = feox2flx(2) + (
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $                 )*dz*poro(iz)*sat(Iz)*1d3
     
          fepy2flx(2) = fepy2flx(2) + (
     $  +(1.0d0-frex)*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $   (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                )*dz*poro(iz)*sat(iz)*1d3
C      $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)
        
        else if (iz ==nz) then
        
          fetflx(2) = fetflx(2) + (
     $  (c2x(iz)-c2(iz))/dt 
     $  )*dz*poro(iz)*sat(iz)*1d3
      
          fediflx(2) = fediflx(2) + (
     $   +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2x(iz)-c2x(iz-1))/(dz**2d0))
     $   +swex*(-dfe3*tora(iz)*(c2(iz-1)-1d0*c2(iz))
     $    /(dz**2d0)
     $   -dfe3/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(c2(iz)-c2(iz-1))/(dz**2d0))
     $             )*dz*poro(iz)*sat(Iz)*1d3
     
          feadvflx(2) = feadvflx(2) + (
     $  + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*c2x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(c2(iz)-c2(iz-1))/dz*swex
     $               )*dz*poro(iz)*sat(Iz)*1d3
     
          feox2flx(2) = feox2flx(2) + (
     $  - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz)
     $  - frex*koxa(iz)*c(iz)*po2(iz)
     $                )*dz*poro(iz)*sat(iz)*1d3
     
          fepy2flx(2) = fepy2flx(2) + (
     $  +(1.0d0-frex)*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  +frex*merge(0.0d0,
     $  (14.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     $    *merge(0.0d0,1.0d0,c2x(iz)<c2th)
     
        
        end if 
      
      end do
      
      feflxsum(2) = fetflx(2) + fediflx(2) +  feadvflx(2) + feox2flx(2)
     $   + fepy2flx(2)
      
      do iz = 1, nz  !================================
        
        if (iz/=nz) then
          
          siltflx = siltflx  + (
     $     (msilx(iz)-msil(iz))/dt)*dz
          
          siladvflx = siladvflx + (
     $      -w*(msilx(iz+1)-msilx(iz))/dz*(1.0d0-swex) 
     $      -w*(msil(iz+1)-msil(iz))/dz*swex/dt*dt2*swpe
     $                       )*dz
          
          sildisflx = sildisflx + (
     $     + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $     + frex*
     $      ksil(iz)*poro(iz)*hr
     $          *100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $                       )*dz
     
        else if (iz==nz) then
          
          siltflx = siltflx + (
     $     (msilx(iz)-msil(iz))/dt
     $                       )*dz
          
          siladvflx = siladvflx + (
     $      -w*(msili-msilx(iz))/dz*(1.0d0-swex)
     $      -w*(msili-msil(iz))/dz*swex*dt2/dt*swpe
     $                       )*dz
         
          sildisflx = sildisflx + (
     $      + (1.0d0-frex)*
     $      ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $      + frex*
     $      ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $                       )*dz
     
        end if 
        
      end do  !================================
      
      silflxsum = siltflx + siladvflx + sildisflx 
      
      do iz = 1, nz
      
        if (.not.((iz == 1).or.(iz==nz))) then
     
          natflx = natflx + (
     $  (nax(iz)-na(iz))/dt
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          naadvflx = naadvflx + (
     $  + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(na(iz)-na(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadiflx = nadiflx + (
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nax(iz-1)-2d0*nax(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(nax(iz)-nax(iz-1))/(dz**2d0))
     $   +swex*(-dna*tora(iz)*(na(iz+1)+na(iz-1)-2d0*na(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(na(iz)-na(iz-1))/(dz**2d0))
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadisflx = nadisflx + (
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
        else if (iz == 1) then
     
          natflx = natflx + (
     $  (nax(iz)-na(iz))/dt 
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          naadvflx = naadvflx + (
     $  + v(iz)*(nax(iz)-nai)/dz*(1.0d0-swex)
C      $  + (v(iz)-v(iz-1))*cx(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(na(iz)-nai)/dz*swex
C      $  + (v(iz)-v(iz-1))*c(iz)/dz*swex
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadiflx = nadiflx + (
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nai-2d0*nax(iz))
     $    /(dz**2d0)
C      $   -dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(nax(iz)-nai)/(dz**2d0)
     $   )
     $   +swex*(-dna*tora(iz)*(na(iz+1)+nai-2d0*na(iz))
     $    /(dz**2d0)
C      $   -dna/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(na(iz)-nai)/(dz**2d0)
     $    )
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadisflx = nadisflx + (
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
        
        else if (iz == nz) then
     
          natflx = natflx + (
     $  (nax(iz)-na(iz))/dt 
     $                       )*dz*poro(iz)*sat(iz)*1d3   ! commented out (is this necessary?)
     
          naadvflx = naadvflx + (
     $  + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*cx(iz)/dz
     $  + swex*v(iz)*(na(iz)-na(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*c(iz)/dz
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadiflx = nadiflx + (
     $   +(1d0-swex)*(-dna*tora(iz)*(nax(iz-1)-1d0*nax(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(nax(iz)-nax(iz-1))/(dz**2d0))
     $   +swex*(-dna*tora(iz)*(na(iz-1)-1d0*na(iz))
     $    /(dz**2d0)
     $   -dna/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(na(iz)-na(iz-1))/(dz**2d0))
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          nadisflx = nadisflx + (
     $  - (1.0d0-frex)*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $           *1d-3
     $  - frex*
     $   ksil(iz)/sat(iz)*hr*100.07d0*1d-6*msil(iz)
     $   *(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)
     $           *1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
        end if ! =========================================
      enddo 
      
      naflxsum = natflx + nadiflx + naadvflx + nadisflx 
      
      do iz = 1, nz
        
        if (.not.((iz == 1).or.(iz==nz))) then
     
          so4tflx = so4tflx+ (
     $  (so4x(iz)-so4(iz))/dt
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
     
          so4diflx = so4diflx+ (
     $   +(-dso4*tora(iz)*(so4x(iz+1)+so4x(iz-1)-2d0*so4x(iz))
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(so4x(iz)-so4x(iz-1))/(dz**2d0))
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4advflx = so4advflx+ (
     $  + v(iz)*(so4x(iz)-so4x(iz-1))/dz
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4disflx = so4disflx+ (
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
        else if (iz == 1) then
     
          so4tflx = so4tflx+ (
     $  (so4x(iz)-so4(iz))/dt 
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4diflx = so4diflx+ (
     $   +(-dso4*tora(iz)*(so4x(iz+1)+so4i-2d0*so4x(iz))
     $    /(dz**2d0)
C      $   -dso4/poro(iz)/sat(iz)
C      $   *(poro(iz+1)*sat(iz+1)*tora(iz+1)-poro(iz)*sat(iz)*tora(iz))
C      $   *(so4x(iz)-so4i)/(dz**2d0)
     $   )
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4advflx = so4advflx+ (
C      $  + (v(iz)-v(iz-1))*so4x(iz)/dz*(1.0d0-swex)
     $  + v(iz)*(so4x(iz)-so4i)/dz
C      $  + (v(iz)-v(iz-1))*so4(iz)/dz*swex
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4disflx = so4disflx+ (
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $    ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
        
        else if (iz == nz) then
     
          so4tflx = so4tflx + (
     $  (so4x(iz)-so4(iz))/dt 
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4diflx = so4diflx + (
     $   +(-dso4*tora(iz)*(so4x(iz-1)-1d0*so4x(iz))
     $    /(dz**2d0)
     $   -dso4/poro(iz)/sat(iz)
     $   *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))
     $   *(so4x(iz)-so4x(iz-1))/(dz**2d0))
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4advflx = so4advflx + (
C      $  + (1.0d0-swex)*(v(iz)-v(iz-1))*so4x(iz)/dz
     $  + v(iz)*(so4x(iz)-so4x(iz-1))/dz
C      $  + swex*(v(iz)-v(iz-1))*so4(iz)/dz
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
          so4disflx = so4disflx + (
     $  - 2d0*(1.0d0-frex)*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)
     $   *merge(0d0,po2x(iz)**0.50d0
     $     ,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3
     $  - 2d0*frex*
     $   koxs(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)
     $   *merge(0d0,po2(iz)**0.50d0
     $   ,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3
     $   
     $  -(1.0d0-frex)*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*msx(iz)*
     $           c2x(iz)**0.93d0*cx(iz)**(-0.40d0),
     $    cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3
     $  -frex*merge(0.0d0,
     $   (2.0d0)*koxs2(iz)/sat(iz)*hr*23.94d0*1d-6*ms(iz)*
     $           c2(iz)**0.93d0*c(iz)**(-0.40d0),
     $    c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3
     $                       )*dz*poro(iz)*sat(iz)*1d3
     
        end if 
        
      end do  ! ==============================
      
      so4flxsum = so4tflx + so4diflx + so4advflx + so4disflx 
      
      
      po2 = po2x
      c = cx
      ms = msx
      c2 = c2x
      so4 = so4x
      na = nax
      msil = msilx
      pro = prox
      
!       print *, (po2(iz),iz=1,Nz,spc)
!       print *, (c(iz),iz=1,Nz,spc)
!       print *, (ms(iz),iz=1,Nz,spc)

C       if (any(it==reclis)) then
      if (time>=rectime(irec+1)) then
        write(chr,'(i3.3)') irec+1
        open (58, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(rate)-'//chr//'.txt', 
     $                              status='replace')
        open (22, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res-'//chr//'.txt', 
     $     status='replace')
        open (29, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(nom)-'//chr//'.txt', 
     $   status='replace')
        
        open(30, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-bsd-'//chr//'.txt', 
     $   status='replace')
        
        
        do iz = 1, Nz
          write (58,*) z(iz),
     $ koxs(iz)*poro(iz)*hr*23.94d0*1d-6*msx(iz)
     $  *po2x(iz)**0.50d0
     $  *merge(0.0d0,1.0d0,po2x(iz)<po2th),
     $ merge(0.0d0,
     $      + koxs2(iz)*poro(iz)*hr*23.94d0*1d-6*msx(iz)
     $         *c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth)
     $  *merge(0.0d0,1.0d0,c2x(iz)<c2th),
     $  poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)
     $  *merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth)
     $  , swbr*vmax*po2x(iz)/(po2x(iz)+mo2)
     $  ,ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $  ,time
          write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz)
     $     ,na(iz),msil(iz),pro(iz),silsat(iz), time
          write (29,*) z(iz),po2(iz)/maxval(po2(:)),c(iz)/maxval(c(:))
     $     ,ms(iz)/maxval(ms(:)),c2(iz)/maxval(c2(:))
     $     , so4(iz)/maxval(so4(:)), na(iz)/maxval(na(:))
     $     , msil(iz)/maxval(msil(:)), pro(iz)/maxval(pro(:))
     $     , silsat(iz),time
          write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz)
        end do
        irec=irec+1
        
        close(58)
        close(22)
        close(29)
        close(30)
        
        write(65,*) time, o2tflx, diflx, advflx, pyoxflx, 
     $     feoxflx, respflx,o2flxsum 
        write(67,*) time, fetflx(1),fediflx(1), feadvflx(1),feox2flx(1), 
     $     fepy1flx(1), fepy2flx(1),feflxsum(1)
        write(68,*) time, fetflx(2),fediflx(2), feadvflx(2),feox2flx(2), 
     $     fepy1flx(2), fepy2flx(2),feflxsum(2) 
        write(69,*) time, pytflx, pyadvflx, pyox1flx, pyox2flx,
     $   pyflxsum        
        write(56,*) time, siltflx, siladvflx, sildisflx,
     $      silflxsum        
        write(55,*) time, natflx, naadvflx, nadiflx, nadisflx,
     $     naflxsum        
        write(57,*) time, so4tflx, so4advflx, so4diflx, so4disflx,
     $     so4flxsum
        
        write(95,*) time, o2in, o2out, po2i,msi,zrxn, zsat
        
      end if
      
      
      write(71,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx
     $  ,o2flxsum
      write(73,*) time, fetflx(1),fediflx(1), feadvflx(1), feox2flx(1), 
     $     fepy1flx(1), fepy2flx(1),feflxsum(1)
      write(74,*) time, fetflx(2),fediflx(2), feadvflx(2), feox2flx(2), 
     $     fepy1flx(2), fepy2flx(2),feflxsum(2)
      write(75,*) time, pytflx, pyadvflx, pyox1flx, pyox2flx,
     $  pyflxsum      
      write(77,*) time, siltflx, siladvflx, sildisflx,
     $  silflxsum      
      write(76,*) time, natflx, naadvflx, nadiflx, nadisflx,
     $  naflxsum      
      write(78,*) time, so4tflx, so4advflx, so4diflx, so4disflx,
     $  so4flxsum      
      
      write(97,*) time, o2in, o2out, po2i,msi,zrxn, zsat
      
      it = it + 1
      time = time + dt
      
!       if (iter<=2) then
!         dt = dt * 1.01d0
!       else if (iter>=5) then
!         dt = dt * 0.99d0
!       end if 
      
      end do
      
      write(chr,'(i3.3)') irec+1
      open (58, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(rate)-'//chr//'.txt', 
     $                              status='replace')
      open (22, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res-'//chr//'.txt', 
     $                             status='replace')
      open (29, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-res(nom)-'//chr//'.txt', 
     $                              status='replace')
        
        open(30, file=trim(adjustl(workdir))//
     $  trim(adjustl(runname))//'/'//
     $ 'o2profile-bsd-'//chr//'.txt', 
     $   status='replace')
        
      do iz = 1, Nz
        write (58,*) z(iz),
     $ koxs(iz)*poro(iz)*hr*23.94d0*1d-6*msx(iz)*
     $  po2x(iz)**0.50d0
     $  *merge(0.0d0,1.0d0,po2x(iz)<po2th),
     $ merge(0.0d0,
     $      + koxs2(iz)*poro(iz)*hr*23.94d0*1d-6*msx(iz)
     $         *c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth)
     $  *merge(0.0d0,1.0d0,c2x(iz)<c2th),
     $  poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)
     $  *merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth)
     $  , swbr*vmax*po2x(iz)/(po2x(iz)+mo2)
     $  ,ksil(iz)*poro(iz)*hr
     $             *100.07d0*1d-6*msilx(iz)
     $   *(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil)
     $   *merge(0d0,1d0
     $  ,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)
     $  ,time
        write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz),so4(iz)
     $   ,na(iz),msil(iz),pro(iz),silsat(iz), time
        write (29,*) z(iz),po2(iz)/maxval(po2(:)),c(iz)/maxval(c(:))
     $   ,ms(iz)/maxval(ms(:)),c2(iz)/maxval(c2(:))
     $   ,so4(iz)/maxval(so4(:)),na(iz)/maxval(na(:))
     $   ,msil(iz)/maxval(msil(:)),pro(iz)/maxval(pro(:))
     $   ,silsat(iz), time
          write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz)
      end do
      
      
      write(65,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx
     $   ,o2flxsum
        
      close(58)
      close(22)
      close(29)
      close(30)
      
      close(65)
      close(67)
      close(68)
      close(69)
      close(55)
      close(56)
      close(57)
      
      close(71)
      close(73)
      close(74)
      close(75)
      close(76)
      close(77)
      close(78)
      
      close(95)
      close(97)
      
      zpy = 0d0
      zab = 0d0
      
      do iz=1,nz    
         if ( zpy(1)==0d0 .and. 
     $    ms(iz)>=0.1d0*msi+0.9d0*ms(1)) 
     $    zpy(1) = z(iz)    
         if ( zpy(2)==0d0 .and. 
     $    ms(iz)>=0.5d0*msi+0.5d0*ms(1)) 
     $    zpy(2) = z(iz)    
         if ( zpy(3)==0d0 .and. 
     $   ms(iz)>=0.9d0*msi+0.1d0*ms(1)) 
     $    zpy(3) = z(iz)    
         if ( zab(1)==0d0 .and. 
     $   msil(iz)>=0.1d0*msili+0.9d0*msil(1)) 
     $    zab(1) = z(iz)    
         if ( zab(2)==0d0 .and. 
     $    msil(iz)>=0.5d0*msili+0.5d0*msil(1)) 
     $    zab(2) = z(iz)    
         if ( zab(3)==0d0 .and. 
     $    msil(iz)>=0.9d0*msili+0.1d0*msil(1)) 
     $    zab(3) = z(iz)    
      enddo
      
      fname=trim(adjustl(workdir))//'sense'
     $ //trim(adjustl(base))//'.txt'
      open(20,file=trim(adjustl(fname)),
     $ action='write',status='unknown',access='append')
      write(20,*) qin,zsat,zpy(1:3),zab(1:3)
      close(20)
#ifdef pyweath 
      fname=trim(adjustl(workdir))//'sense-py'
     $ //trim(adjustl(base))//'.txt'
      open(400,file=trim(adjustl(fname)),
     $ action='write',status='unknown',access='append')
      write(400,*) qin,zsat,zpy(1:3)
      close(400)
#endif      
#ifdef silweath 
      fname=trim(adjustl(workdir))//'sense-ab'
     $ //trim(adjustl(base))//'.txt'
      open(401,file=trim(adjustl(fname)),
     $ action='write',status='unknown',access='append')
      write(401,*) qin,zsat,zab(1:3)
      close(401)
#endif
        
      end
      
      subroutine checkfile(fname,oxj)
      implicit none
      character(*),intent(in)::fname
      integer,intent(out)::oxj
 
      open(998,file=trim(fname),status='old',err=999)
      close(998)
      write(6,'(3A)')"file '",fname,"' exist"
      oxj=1
      return
 
  999 continue
      close(998)
      write(6,'(3A)')"file '",fname,"' don't exist"
      oxj=0

      return
      end subroutine checkfile