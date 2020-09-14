program o2profile
!! try to calculate o2 profile in a simple system
!! full implicit method
!! from ver 8 reflecting ph change in next step of pyrite weathering and assuming always wet surface 
!! from ver 9 deterministic pH calculation; pH is not iteratively sought 
!! from ver 9.5 pH is iteratively calculated to be consistent with pyrite; for the case excluding aqueous Fe2 oxidation 
!! from ver 9.6 calculation is conducted within the whole domain 
!! adding HCO3 profile recording
!! adding calcium profile to examine the effect of CaCO3 dissolution
!! allowing porosity & surface area change       
!! converted to f90 from o2profile+silweath+o2_v9_7.f
implicit none

!-----------------------------

real(kind=8) :: ztot = 20.0d0 ! m
real(kind=8) dz
integer, parameter :: nz = 100 
real(kind=8) z(nz)
real(kind=8) ze(nz+1)
real(kind=8) :: ph = 5.0d0
real(kind=8) :: tc = 15.0d0 ! deg celsius
real(kind=8) dt  ! yr 
integer, parameter :: nt = 50000000
real(kind=8) time
integer, parameter :: nsp = 4
integer, parameter :: nsp3= 2

real(kind=8) :: rg = 8.3d-3   ! kJ mol^-1 K^-1
real(kind=8) :: rg2 = 8.2d-2  ! L mol^-1 atm K^-1

real(kind=8) :: po2i = 0.21d0 ! atm **default
! real(kind=8) :: po2i = 0.6d-1 ! atm
real(kind=8) :: pco2i = 10.0d0**(-2.5d0) ! atm **default 
! real(kind=8) :: pco2i = 10.0d0**(-1.0d0) ! atm
real(kind=8) :: ci = 0d0 
real(kind=8) :: c2i = 0d0
real(kind=8) :: so4i = 0d0
real(kind=8) :: nai = 0d0
real(kind=8) :: mgi = 0d0
real(kind=8) :: sii = 0d0

real(kind=8) :: mvfo = 43.79d0 ! cm3/mol; molar volume of Fo; Robie et al. 1978
real(kind=8) :: mvab = 100.07d0 ! cm3/mol; molar volume of Ab; Robie et al. 1978

real(kind=8) :: mwtfo = 140.694d0 ! g/mol; molar volume of Fo; Robie et al. 1978
real(kind=8) :: mwtab = 262.225d0 ! g/mol; molar volume of Ab; Robie et al. 1978

! real(kind=8) :: redsldi = 0.56d0 ! wt%  **default 
! real(kind=8) :: redsldi = 1.12d0 ! wt%  x2
real(kind=8) :: redsldi = 2.8d0 ! wt%   x5
! real(kind=8) :: redsldi = 2.24d0 ! wt%  x4
! real(kind=8) :: redsldi = 3.36d0 ! wt%  x6

real(kind=8) :: silwti = 30d0 ! wt%  **default
! real(kind=8) :: silwti = 45d0 ! wt%  
! real(kind=8) :: silwti = 24d0 ! wt%
! real(kind=8) :: silwti = 1d-10 ! wt%

real(kind=8)::rainpowder = 40d2 !  g/m2/yr corresponding to 40 t/ha/yr (40x1e3x1e3/1e4)
! real(kind=8)::rainpowder = 0.5d2 !  g/m2/yr corresponding to 0.5 t/ha/yr (0.5x1e3x1e3/1e4)

real(kind=8)::rainfrc_fo = 0.12d0 ! rain wt fraction for Fo (Beering et al 2020)
real(kind=8)::rainfrc_ab = 0.172d0 ! rain wt fraction for Ab; assuming 0.43 for La and 0.4 of wt of La is Ab (Beering et al 2020)

real(kind=8)::zsupp = 0.3d0 !  e-folding decrease

real(kind=8) sat(nz), poro(nz), torg(nz), tora(nz), deff(nz)
real(kind=8) :: dgas = 6.09d2 ! m^2 yr^-1
real(kind=8) :: daq = 5.49d-2 ! m^2 yr^-1
real(kind=8) :: poroi = 0.1d0
real(kind=8) :: sati = 0.50d0
real(kind=8) :: satup = 0.10d0

! real(kind=8) :: zsat = 30d0  ! water table depth [m] ** default 
real(kind=8) :: zsat = 5d0  ! water table depth [m] 
! real(kind=8) :: zsat = 5d0

real(kind=8) :: dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
real(kind=8) :: dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dmg  = 0.017218079d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dsi  = 0.03689712d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 

real(kind=8), parameter :: w = 5.0d-5 ! m yr^-1, uplift rate ** default 
! real(kind=8), parameter :: w = 1.0d-5 ! m yr^-1, uplift rate

real(kind=8), parameter :: vcnst = 1.0d1 ! m yr^-1, advection
! real(kind=8), parameter :: qin = 5d-3 ! m yr^-1, advection (m3 water / m2 profile / yr)

real(kind=8) :: qin = 1d-1 ! m yr^-1, advection (m3 water / m2 profile / yr)  ** default
! real(kind=8) :: qin = 2d-1 ! m yr^-1 
real(kind=8) v(nz), q

! real(kind=8) :: hr = 1d5 ! m^2 m^-3, reciprocal of hydraulic radius  ** default 
! real(kind=8) :: hr = 1d4 ! m^2 m^-3, reciprocal of hydraulic radius
real(kind=8) :: hri = 1d5

real(kind=8) :: hr(nz)

real(kind=8) po2(nz), redsld(nz), redaq(nz), ms(nz), c(nz)
real(kind=8) po2x(nz), msx(nz), cx(nz)
real(kind=8) msi
real(kind=8) msili
real(kind=8) mfoi,mabi
real(kind=8),dimension(nz)::msil,msilx,mfo,mfox,mab,mabx

real(kind=8) c2(nz), c2x(nz),ctmp, po2tmp
real(kind=8) so4(nz), so4x(nz)
real(kind=8) na(nz), nax(nz), naeq(nz), silsat(nz) 
real(kind=8) pro(nz), prox(nz), dumreal(nz), dprodna(nz)
real(kind=8) hco3(nz), ca(nz), co2(nz), co3(nz), dic(nz)
real(kind=8) porox(nz), dporodta(nz),dporodtg(nz)
real(kind=8) mg(nz),mgx(nz), si(nz), six(nz)

real(kind=8) :: caeq = 1d-3  ! mol/L equilibrium Ca conc. 
real(kind=8) :: delca = 0.5d0  ! m reaction front width  
real(kind=8) :: zca = 50d0   ! m depth of reaction front for calcite          

real(kind=8) koxa(nz), kdis(nz), koxs(nz),koxs2(nz)
real(kind=8) koxai(nz), kdisi(nz), koxsi(nz),koxs2i(nz)
real(kind=8) ksil(nz), msilsupp(nz)
real(kind=8) kfo(nz), mfosupp(nz), omega_fo(nz) 
real(kind=8) kab(nz), mabsupp(nz), omega_ab(nz) 

real(kind=8) kho, ucv
real(kind=8) kco2,k1, keqsil, kw, k2, kcceq, keqfo, keqab, keqgb

integer iz, ie, it, ie2

integer col, row

real(kind=8) imbr

integer, parameter :: nmx = nsp*nz
integer, parameter :: nmx2 = 1*nz
integer, parameter :: nmx3 = nsp3*nz

real(kind=8) amx(nmx,nmx),ymx(nmx)
real(kind=8) amx2(nmx2,nmx2),ymx2(nmx2)
real(kind=8) amx3(nmx3,nmx3),ymx3(nmx3)
real(kind=8) emx3(nz)
! real(kind=8) y2mx(3*(nz-1)),a2mx(3*(nz-1),3*(nz-1)) 
! real(kind=8), allocatable :: y2mx(:),a2mx(:,:) 

integer ipiv(nmx)
integer ipiv2(nz)
integer ipiv3(nmx3)
! integer, allocatable :: ipiv(:)
integer info

external DGESV

real(kind=8) error, error2
real(kind=8) :: tol = 1d-6

integer :: spc = 50

! integer, parameter :: nrec = 22
integer, parameter :: nrec = 20
integer reclis(nrec)
real(kind=8) rectime(nrec)
character(3) chr
character(256) runname,workdir, chrz(3), chrq(3),base,fname
integer irec, iter, idum, iter2

real(kind=8) :: swad = 1.0d0! 1.0 when advection included 0.0d0 when not
real(kind=8) dt2, swpe  ! physical erosion 

! integer, parameter :: nt2 = int(nt*w/v)
! real(kind=8) time2(nt2)

integer zlis(3*(nz-1))
integer nel, imx, imx2

real(kind=8) :: po2th = 1.0d-20
real(kind=8) minpo2

real(kind=8) :: cth = 1.0d-20
real(kind=8) :: c2th = 1.0d-20
real(kind=8) :: so4th = 1.0d-20
real(kind=8) :: proth = 1.0d-20
real(kind=8) :: nath = 1.0d-20
real(kind=8) :: mgth = 1.0d-20
real(kind=8) :: sith = 1.0d-20
real(kind=8) :: msth = 1.0d-300
real(kind=8) :: msilth = 1.0d-300
real(kind=8) :: mfoth = 1.0d-300
real(kind=8) :: mabth = 1.0d-300

real(kind=8) prepo2
real(kind=8) :: stoxs = 15.0d0/4.0d0  ! 15/4 py => Fe-oxide + sulfate; 7/2 py => Fe++ + sulfate
real(kind=8) :: stoxa = 1.0d0/4.0d0  ! stoichiomety of oxidation in aq
real(kind=8) :: swoxa = 0.0d0   ! switch for oxidation in aq
real(kind=8) :: swoxs2 = 1.0d0  ! switch for oxidation in solid phase
real(kind=8) :: swoxall = 0d0   ! switch when only consider overall oxidation, i.e., 15/4 py => Fe-oxide + sulfate

real(kind=8) :: swbr = 0.0d0  ! switch for biological respiration
real(kind=8) :: vmax = 0.71d0 ! mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
real(kind=8) :: mo2 = 0.121d0 ! Michaelis, Davidson et al. (2012)

real(kind=8) :: swex = 0.0d0 ! switch for explicit
real(kind=8) :: frex = 0.0d0 ! fraction of explicit

real(kind=8) :: swadvmass = 0d0 ! switch; 1 when calculating q from advection mass balance
real(kind=8) :: waterfluc = 0d0 ! switch: 1 when fluctuating water flow

real(kind=8) pyoxflx, feoxflx, respflx, diflx, advflx, o2tflx
real(kind=8) feadvflx(2), feox2flx(2), fepy1flx(2) &
&  , fepy2flx(2), fetflx(2), fediflx(2)
real(kind=8) pyadvflx, pyox1flx, pyox2flx, pytflx
real(kind=8) siladvflx, sildisflx,  siltflx
real(kind=8) naadvflx, nadisflx,  natflx, nadiflx
real(kind=8) so4advflx, so4disflx,  so4tflx, so4diflx
real(kind=8) o2flxsum,feflxsum(2),pyflxsum, so4flxsum
real(kind=8) silflxsum,naflxsum
real(kind=8) flx_t(nz),flx_adv(nz),flx_dif(nz),flx_oxaq(nz)
real(kind=8) flx_oxpy1(nz),flx_oxpy2(nz),flx_rem(nz)
real(kind=8) dummy, zdum(nz)

integer,parameter :: nflx = 6
integer  iflx
real(kind=8),dimension(nflx,nz)::flx_fo,flx_mg,flx_si,flx_ab,flx_na

! real(kind=8) :: maxdt = 10d0
real(kind=8) :: maxdt = 0.2d0 ! for basalt exp?

integer izdum
! logical :: pre_calc = .false.
logical :: pre_calc = .true.

logical :: read_data = .false.
! logical :: read_data = .true.

data rectime /1d1,3d1,1d2,3d2,1d3,3d3,1d4,3d4 &
    & ,1d5,2d5,3d5,4d5,5d5,6d5,7d5,8d5,9d5,1d6,1.1d6,1.2d6/
! data rectime /-1d6,0d6,1d6,2d6,3d6,4d6,5d6,6d6,7d6,8d6
! &,9d6,10d6,11d6,12d6,13d6,14d6,15d6,16d6,17d6,18d6,19d6,20d6/
! data rectime /21d6,22d6,23d6,24d6,25d6,26d6,27d6,28d6,29d6,30d6
! & ,31d6,32d6,33d6,34d6,35d6,36d6,37d6,38d6,39d6,40d6,41d6,42d6/
logical :: dir_exist
real(kind=8) o2out,zrxn,dms_max
real(kind=8) :: beta = 3.7d19/0.21d0
real(kind=8) :: o2in = 1d13
real(kind=8) :: area = 1.5d14
real(kind=8) pyrxn(nz)
real(kind=8) :: pi = 4.0d0*datan(1d0)

logical :: O2_evolution = .false.
! logical :: O2_evolution = .true.

logical :: flgback = .false.
real(kind=8) :: zab(3), zpy(3) 
integer :: oxj

real(kind=8), parameter:: infinity = huge(0d0)
!-------------------------

! rectime =rectime/1d4 
rectime =rectime/1d3 ! better with basalt exp?
! rectime =rectime/2d0
! rectime=rectime*3d0 

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


write(workdir,*) '../pyweath_output/'     
write(base,*) '_basalt_test_cpl_high_rain'     
#ifdef test 
write(base,*) '_test'
#endif     
write(runname,*) 'Fe2+SO4+sil+ph_wet_iter'//'---q'//trim(adjustl(chrq(3)))//'_zsat' &
    & //trim(adjustl(chrz(3)))//trim(adjustl(base))
#ifdef pyweath
write(runname,*) 'Fe2+SO4_wet_iter'//'---q'//trim(adjustl(chrq(3)))//'_zsat'  &
    & //trim(adjustl(chrz(3)))//trim(adjustl(base))
#endif
#ifdef silweath
write(runname,*) 'sil+ph_wet_iter'//'---q'//trim(adjustl(chrq(3)))//'_zsat'  &
    & //trim(adjustl(chrz(3)))//trim(adjustl(base))
#endif

call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))

open(65, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(o2flx).txt', &
    & status='unknown', action = 'write')
open(67, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe2flx).txt', &
    & status='unknown', action = 'write')
open(68, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe3flx).txt', &
    & status='unknown', action = 'write')
open(69, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(pyflx).txt',  &
    & status='unknown', action = 'write')
open(55, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx).txt',  &
    & status='unknown', action = 'write')
open(56, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(silflx).txt',  &
    & status='unknown', action = 'write')
open(57, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(so4flx).txt',  &
    & status='unknown', action = 'write')
open(58, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(mgflx).txt',  &
    & status='unknown', action = 'write')
open(59, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(siflx).txt',  &
    & status='unknown', action = 'write')
open(60, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(foflx).txt',  &
    & status='unknown', action = 'write')
open(61, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx2).txt',  &
    & status='unknown', action = 'write')
open(62, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(abflx).txt',  &
    & status='unknown', action = 'write')

open(71, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(o2flx-alltime).txt',  &
    & status='unknown', action = 'write')
open(73, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe2flx-alltime).txt',  &
    & status='unknown', action = 'write')
open(74, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe3flx-alltime).txt',  &
    & status='unknown', action = 'write')
open(75, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(pyflx-alltime).txt',  &
    & status='unknown', action = 'write')
open(76, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx-alltime).txt',  &
    & status='unknown', action = 'write')
open(77, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(silflx-alltime).txt',  &
    & status='unknown', action = 'write')
open(78, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(so4flx-alltime).txt', &
    & status='unknown', action = 'write')


open(95, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2sense-res.txt', & 
    & status='unknown', action = 'write')
open(97, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2sense-res(alltime).txt', & 
    & status='unknown', action = 'write')


do iz=1,nz+1
    ze(iz)=(iz-1)*ztot/nz
end do

z(:) = 0.5d0*(ze(2:nz+1)+ze(1:nz))

dz=z(2)-z(1)

if (swoxall==1d0) swoxa = 0d0

stoxs = stoxs - swoxa*stoxa
#ifdef poroevol 
stoxs = stoxs - (1d0-swoxall)*stoxa  ! 15/4 if overall oxidation is considered (swoxall = 1)
                           ! 7/2  if pyrite oxidation and aqFe2+ oxidation is separately considered (swoxall = 0)
#endif 

do irec=1,nrec
reclis(irec)=irec*nt/nrec
end do

if (tc/=15d0) then 
    dgas = dgas*exp(-4.18d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    daq = daq*exp(-20.07d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    dfe2=dfe2*exp(-19.615251d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    dfe3=dfe3*exp(-14.33659d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    dso4=dso4*exp(-20.67364d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    dna=dna*exp(-20.58566d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
    dmg=dmg*exp(-18.51979d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
endif 

po2 = po2i
! po2 = po2th
sat = sati
! sat = min(1.0d0,0.90d0*(z-ztot/2d0)/ztot/2d0 + 1.d0)
! zsat = ztot/6.0d0
! zsat = 5d0
sat = min(1.0d0,(1d0-satup)*z/zsat + satup)
#ifdef satconvex 
sat = min(1.0d0, satup+(1d0-satup)*(z/zsat)**2d0)
#endif 
#ifdef satconcave 
sat = min(1.0d0, 1d0-(1d0-satup)*(1d0-z/zsat)**2d0)
do iz=1,nz
    if (z(iz)>=zsat) sat(iz)=1d0
enddo 
#endif 
! zca = zsat 
ca = 0d0
#ifdef carbonate
ca = min(caeq,max(0d0,caeq*(z-zca)/delca))
#endif 

hr = hri
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
mg = 0d0
si = 0d0
pro = 10.0d0**(-ph)

dt = maxdt
dt = 1d-2  
dt = 1d-6 ! for basalt exp?


if (swex == 1.0d0) then 

    dt = dz/vcnst
    dt2 = dz/w

end if

kho=10.0d0**(-2.89d0)*exp(13.2d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+25.0d0))/rg)

kco2 = 10.0d0**(-1.34d0)  ! 15C Kanzaki Murakami 2015
k1 = 10.0d0**(-6.42d0)  ! 15C Kanzaki Murakami 2015
k2 = 10d0**(-10.43d0)      ! 15C Kanzaki Murakami 2015
kw = 10.0d0**(-14.35d0)  ! 15C Kanzaki Murakami 2015
kcceq = 10d0**(-8.43d0)  ! 15C Kanzaki Murakami 2015

if (tc==5d0) then

    kco2 = 10.0d0**(-1.19d0)  ! 15C Kanzaki Murakami 2015
    k1 = 10.0d0**(-6.52d0)  ! 15C Kanzaki Murakami 2015
    k2 = 10d0**(-10.55d0)      ! 15C Kanzaki Murakami 2015
    kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0+0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
    kw =10d0**kw
    kcceq = 10d0**(-8.39d0)  ! 15C Kanzaki Murakami 2015

elseif (tc==25d0) then 

    kco2 = 10.0d0**(-1.47d0)  ! 15C Kanzaki Murakami 2015
    k1 = 10.0d0**(-6.35d0)  ! 15C Kanzaki Murakami 2015
    k2 = 10d0**(-10.33d0)      ! 15C Kanzaki Murakami 2015
    kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0 &
        & +0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
    kw =10d0**kw
    kcceq = 10d0**(-8.48d0)  ! 15C Kanzaki Murakami 2015

endif 

pro = sqrt(kco2*k1*pco2i+kw)


ksil = 1.31d-9*1d4  ! mol/m2/yr  ! from Li et al., 2014

keqsil = 3.412182823d0 - 0.5d0* 8.310989613d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
keqsil = 10.0d0**(keqsil)

kab = ksil
keqab = keqsil

kfo = -10.64d0  ! mol/m2/sec ! from Beering et al 2020 (only neutral weathering)
kfo = 10d0**(kfo)*60d0*60d0*24d0*365d0 ! mol/m2/yr 
! kfo = ksil

keqfo= 27.8626d0  ! Sugimori et al. (2012) 
keqfo = 10d0**keqfo
! -208.5932252 

if (tc==5d0) then 
    ksil =  5.13d-10*1d4  ! mol/m2/yr  ! from Li et al., 2014
    keqsil = 3.751169218d0 - 0.5d0* 9.242748918d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
    keqsil = 10.0d0**(keqsil)
elseif (tc==25d0) then
    ksil =  3.15d-9*1d4  ! mol/m2/yr  ! from Li et al., 2014
    keqsil = 3.080792422d0 - 0.5d0* 7.434511289d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
    keqsil = 10.0d0**(keqsil)
endif       


!       print *,keqsil;stop

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

mabi = msili
mab = mabi

mfoi = 1d-10

mfo = mfoi

koxsi = 10.0d0**(-8.19d0)*60.0d0*60.0d0*24.0d0*365.0d0  &!! excluding the term (po2**0.5)
    & *(kho)**(0.50d0)/((10.0d0**(-ph))**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)

koxai = swoxa*8.0d13*60.0d0*24.0d0*365.0d0   &!  excluding the term (c*po2)
    & *(10.0d0**(-14.0d0+ph))**2.0d0               ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)

koxs2i = swoxs2*10.0d0**(-6.07d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (fe3**0.93/fe2**0.40)
!! mol m^-2 yr^-1, Williamson and Rimstidt (1994)

! when adv o2 flux == adv pyrite flx (q*kho*po2i*1d3 = stoxs*msi*w)
if (swadvmass == 1d0) then 
    q = 15d0/4d0*msi*w/kho/po2i/1d3
    v = q/poroi/sat
endif 
! *** the above must be commented out when using arbitrary q value


open(22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd.txt', &
    & status='unknown', action = 'write')

do iz = 1,nz
    write(22,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz)
enddo

close(22)


!  --------- read -----
if (read_data) then 
    open(255, file=trim(adjustl(workdir))//'po2-sense_z1800_qmbs_dt2e-2_wt30_ver6'//'-go2_v2_cnt7'//'/'// &
        & 'o2profile-res-018.txt',status ='old',action='read')
    do iz = 1,nz
        read(255,*) zdum(iz),po2x(iz), cx(iz), msx(iz),c2x(iz), time
    enddo
    close(255)

    if (zdum(nz) > ztot) then  ! interpolating to finer scale
        izdum = 1
        do iz = 1,nz-1
            do while (z(izdum)<=zdum(iz+1))  
                po2(izdum) = (po2x(iz)-po2x(iz+1))/(zdum(iz)-zdum(iz+1))*(z(izdum)-zdum(iz+1))+po2x(iz+1)
                c(izdum) = (cx(iz)-cx(iz+1))/(zdum(iz)-zdum(iz+1))*(z(izdum)-zdum(iz+1))+cx(iz+1)
                ms(izdum) = (msx(iz)-msx(iz+1))/(zdum(iz)-zdum(iz+1))*(z(izdum)-zdum(iz+1))+msx(iz+1)
                c2(izdum) = (c2x(iz)-c2x(iz+1))/(zdum(iz)-zdum(iz+1))*(z(izdum)-zdum(iz+1))+c2x(iz+1)
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
    !       c = cth
    !       c2 = c2th

    open(255, file=trim(adjustl(workdir))//'po2-sense_z1800_qmbs_dt2e-2_wt30_ver6'//'-go2_v2_cnt7'//'/'// &
    & 'o2sense-res(alltime).txt',status ='old',action='read')
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
!       if (.not.read_data)time = 0

!       do irec = 1,nrec-1
!         if ((rectime(irec)-time)*(rectime(irec+1)-time)<0d0) then 
!             idum=irec
!             exit
!         endif
!         idum = 0
!       enddo

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
print *, 'it, time = ',it, time
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

!        if ((iter <= 10).and.(dt<1d1)) then
if (dt<maxdt) then
    dt = dt*1.01d0
    if (dt>maxdt) dt = maxdt
endif
if (iter > 300) then
    dt = dt/1.05d0
end if

! ======== modifying maxdt ===============
!       if (time >1d4) maxdt = 10d0
! ========================================
#ifdef pHiter      
!  ############## pH interation #######################

error2 = 1d4
iter2 = 0
do while (error2 > tol ) 
#endif

koxsi = 10.0d0**(-8.19d0)*60.0d0*60.0d0*24.0d0*365.0d0  &!! excluding the term (po2**0.5)
    & *(kho)**(0.50d0)/(pro**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)

koxai = max(swoxa*8.0d13*60.0d0*24.0d0*365.0d0   &!  excluding the term (c*po2)
    & *(kw/pro)**2.0d0               &! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)
    & , swoxa*1d-7*60.0d0*24.0d0*365.0d0)

if (tc/=15d0) then 
    koxsi = koxsi*exp(-57d0*(1.0d0/(273.0d0+tc)-1.0d0/(273.0d0+15.0d0))/rg)
endif 

koxs = koxsi
koxa = koxai

koxs2 = koxs2i

! forcing and oxygen calculation 
if (O2_evolution) then 
    pyrxn = 0d0
    do iz=1,nz
        pyrxn(iz)=koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0 &
            & *merge(0.0d0,1.0d0,po2(iz)<po2th)
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
        ! zsat = 30d0 + 20d0*sin(2d0*pi*time/1d4)  ! fluctuation
        zsat = 50d0 ! step
        ! zsat = 30d0 + 20d0*(time/20d6)  ! linear increase 
        sat = min(1.0d0,(1d0-satup)*z/zsat + satup)
        if (swadvmass == 1d0) then 
            q = 15d0/4d0*msi*w/kho/po2i/1d3
            v = q/poroi/sat
        endif 
        poro = poroi
        torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
        tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
        deff = torg*dgas + tora*daq

        ! o2in = 2.0d13                 !  doubling
        ! o2in = 1.5d13                 !  1.5 times
        ! o2in = 0.5d13                 !  halving
        o2in = 1.0d13                 !  halving
        ! o2in = 2d13/2d7*time+1d13   ! case if o2in is continuously changed
        if (.not.it==0) o2out = area*pyoxflx
        ! msi = 4d0*o2in/15d0/w/area  ! case if msi is modified with o2in
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


!       po2(1) = po2i
!       ms(nz) = msi
!       msil(nz) = msili

po2x = po2
cx = c
msx = ms

c2x = c2

so4x=so4

nax = na
msilx = msil

six = si
mgx = mg
mfox = mfo

mabx = mab

prox = pro  

porox = poro

error = 1d4
iter=0



mfosupp = rainpowder*rainfrc_fo/mwtfo*exp(-z/zsupp)/zsupp 
mabsupp = rainpowder*rainfrc_ab/mwtab*exp(-z/zsupp)/zsupp 

msilsupp = 0d0


if (pre_calc) then 

    do iz = 1, nz

        if (po2x(iz)>=po2th) cycle

        if (iz/=nz) po2tmp = po2(iz+1)
        if (iz==nz) po2tmp = po2(iz)

        if (iz/=1) then 
            po2x(iz) = max(0.0d0  &
                & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))* &
                & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(-po2(iz))/dt &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2tmp+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2(iz)-po2(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad &
                & ) &
                & )
        else ! iz == 1
            po2x(iz) = max(0.0d0 &
                & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))* &
                & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3) &
                & *(-po2(iz))/dt-(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(-2.0d0*po2(iz) + po2(iz+1)+po2i) &
                & /(dz**2.0d0)+poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*swad  &
                & ) &
                & )
        endif

    end do 

    do iz = 1, nz

        if (msx(iz)>=msth) cycle

        if (iz/=nz) then 
            msx(iz) = max(0d0, &
                & ms(iz) +dt*(w*(ms(iz+1)-ms(iz))/dz) &
                & )
        else  
            msx(iz) = max(0d0, &
                & ms(iz) +dt*(w*(msi-ms(iz))/dz) &
                & )
        endif

        if (msilx(iz)>=msilth) cycle

        if (iz/=nz) then 
            msilx(iz) = max(0d0, &
                & msil(iz) +dt*(w*(msil(iz+1)-msil(iz))/dz + msilsupp(iz)) &
                & )
        else 
            msilx(iz) = max(0d0, &
                & msil(iz) + dt*(w*(msili-msil(iz))/dz + msilsupp(iz)) &
                & )
        endif 

        if (mfox(iz)>=mfoth) cycle

        if (iz/=nz) then 
            mfox(iz) = max(0d0, &
                & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz + mfosupp(iz)) &
                & )
        else 
            mfox(iz) = max(0d0, &
                & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz+ mfosupp(iz)) &
                & )
        endif 

    enddo


    ! pause

    if (swoxa == 1d0) then 

        do iz = 1, nz

            if (c2x(iz)>=c2th) cycle

            if (iz/=nz) ctmp = c2(iz+1)
            if (iz==nz) ctmp = c2(iz)
            
            if (iz/=1) then 
                c2x(iz) = max(0.0d0, &
                    & c2(iz) +dt*(-v(iz)*(c2(iz)-c2(iz-1))/dz+dfe3*tora(iz)*(ctmp+c2(iz-1)-2d0*c2(iz))/(dz**2d0) &
                    & +dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(c2(iz)-c2(iz-1))/(dz**2d0)+koxa(iz)*po2x(iz)*c(iz)) &
                    & )
            else
                c2x(iz) = max(0.0d0, &
                    & c2(iz) + dt*(-v(iz)*(c2(iz)-c2i)/dz+dfe3*tora(iz)*(ctmp+c2i-2d0*c2(iz))/(dz**2d0) &
                    & +koxa(iz)*po2x(iz)*c(iz) &
                    & ) &
                    & )
            end if
        enddo

        do iz = 1, nz

            if (cx(iz)>=cth) cycle

            if (iz/=nz) ctmp = c(iz+1)
            if (iz==nz) ctmp = c(iz)
            
            if (iz/=1) then 
                cx(iz) = max(0.0d0, &
                    & c(iz) + dt*(-v(iz)*(c(iz)-c(iz-1))/dz+dfe2*tora(iz)*(ctmp+c(iz-1)-2d0*c(iz))/(dz**2d0) &
                    & +dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))  &
                    & *(c(iz)-c(iz-1))/(dz**2d0)            &
                    & +15d0*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3  &
                    & +koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3  &
                    & )  &
                    & )
            else 
                cx(iz) = max(0.0d0, &
                    & c(iz) + dt*(-v(iz)*(c(iz)-ci)/dz+dfe2*tora(iz)*(ctmp+ci-2d0*c(iz))/(dz**2d0) &
                    & +15d0*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
                    & +koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
                    & ) &
                    & )
            endif
        enddo
        
        do iz = 1, nz
            if (so4x(iz)>=so4th) cycle

            if (iz/=nz) ctmp = so4(iz+1)
            if (iz==nz) ctmp = so4(iz)
            if (iz/=1) then 
                so4x(iz) = max(0.0d0, &
                    & so4(iz) +dt*(-v(iz)*(so4(iz)-so4(iz-1))/dz+dso4*tora(iz)*(ctmp+so4(iz-1)-2d0*so4(iz))/(dz**2d0) &
                    & +dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(so4(iz)-so4(iz-1))/(dz**2d0) &
                    & +2d0*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
                    & +2d0*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
                    & ) &
                    & )
            else 
                so4x(iz) = max(0.0d0, &
                    & so4(iz) +dt*(-v(iz)*(so4(iz)-so4i)/dz+dso4*tora(iz)*(ctmp+so4i-2d0*so4(iz))/(dz**2d0) &
                    & +2d0*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
                    & +2d0*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
                    & ) &
                    & )
            endif 
        enddo
    endif 

    do iz = 1, nz

        if (nax(iz)>=nath) cycle

        if (iz/=nz) ctmp = na(iz+1)
        if (iz==nz) ctmp = na(iz)
        if (iz/=1) then 
            nax(iz) = max(0.0d0, &
                & na(iz) +dt*(-v(iz)*(na(iz)-na(iz-1))/dz+dna*tora(iz)*(ctmp+na(iz-1)-2d0*na(iz))/(dz**2d0) &
                & +dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                & *(na(iz)-na(iz-1))/(dz**2d0) &
                & +ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*1d-3 &
                & ) &
                & )
        else 
            nax(iz) = max(0.0d0, &
                & na(iz) + dt*(-v(iz)*(na(iz)-nai)/dz+dna*tora(iz)*(ctmp+nai-2d0*na(iz))/(dz**2d0) &
                & +ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*1d-3 &
                & ) &
                & )
        endif 
    enddo
            
    do iz = 1, nz

        if (mgx(iz)>=mgth) cycle

        if (iz/=nz) ctmp = mg(iz+1)
        if (iz==nz) ctmp = mg(iz)
        if (iz/=1) then 
            mgx(iz) = max(0.0d0, &
                & mg(iz) +dt*(-v(iz)*(mg(iz)-mg(iz-1))/dz+dmg*tora(iz)*(ctmp+mg(iz-1)-2d0*mg(iz))/(dz**2d0) &
                & +dmg/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                & *(mg(iz)-mg(iz-1))/(dz**2d0) &
                & +2d0*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*1d-3 &
                & ) &
                & )
        else 
            mgx(iz) = max(0.0d0, &
                & mg(iz) + dt*(-v(iz)*(mg(iz)-mgi)/dz+dmg*tora(iz)*(ctmp+mgi-2d0*mg(iz))/(dz**2d0) &
                & +2d0*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*1d-3 &
                & ) &
                & )
        endif 
    enddo
            
            
    do iz = 1, nz
        if (six(iz)>=sith) cycle

        if (iz/=nz) ctmp = si(iz+1)
        if (iz==nz) ctmp = si(iz)
        if (iz/=1) then 
            six(iz) = max(0.0d0, &
                & si(iz) +dt*(-v(iz)*(si(iz)-si(iz-1))/dz+dsi*tora(iz)*(ctmp+si(iz-1)-2d0*si(iz))/(dz**2d0) &
                & +dsi/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                & *(si(iz)-si(iz-1))/(dz**2d0) &
                & +kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*1d-3 &
                & ) &
                & )
        else 
            six(iz) = max(0.0d0, &
                & si(iz) + dt*(-v(iz)*(si(iz)-sii)/dz+dsi*tora(iz)*(ctmp+sii-2d0*si(iz))/(dz**2d0) &
                & +kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*1d-3 &
                & ) &
                & )
        endif 
            
    end do

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

#ifndef silweath      

if (it == 0 .and. iter == 0) then
    ! cx(:) = 1.0d2
    ! c2x(:) = 1.0d2
    so4x(:) = 1.0d2
    nax(:) = 1.0d2
    mgx(:) = 1.0d2
    six(:) = 1.0d2
end if

call pyweath_1D( &
    & nz,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
    & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daq,dgas,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
    & ,swoxa,swoxall,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,cx,c2x,po2x,msx,flx_t,flx_adv,flx_dif,flx_oxaq,flx_oxpy1,flx_oxpy2,flx_rem &! output
    ) 
!!!!!!!!!!!!!!!!!!!!!!!!  so4 calculation start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call pyweath_1D_SO4( &
    & nz,c,c2,po2,ms,hr,po2th,poro,z,dz,koxs2,koxs,dso4,sat,dporodta  &! input
    & ,cth,tora,v,tol,zrxn,dt,cx,c2x,po2x,msx,so4,swoxa,O2_evolution,so4i,so4th &! input
    & ,so4x &! output
    & )

!!!  =-=-=-=-=-=-=-=-=- END of SO4 calculation  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif       
#ifndef pyweath       
!!! =-=-=-=-=-=-=-=-=- START of calculation for Na and albite  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#ifdef silweath
so4x = so4th
cx = cth 
c2x = c2th 

if (it == 0 .and. iter == 0) then 
    nax(1:) = 1.0d2
    mgx(1:) = 1.0d2
    six(1:) = 1.0d2
endif 

#endif      

! call abweath_1D( &
    ! & nz,na,msil,hr,poro,z,dz,w,ksil,keqsil,msilth,dna,sat,dporodta,pro,msili,msilsupp  &! input
    ! & ,kco2,k1,k2,dt2,nath,tora,v,tol,nsp3,zrxn,it,cx,c2x,so4x,ca,pco2i,nai,mgx &! input
    ! & ,iter,error,dt,flgback &! inout
    ! & ,nax,prox,co2,hco3,co3,naeq,silsat,dic,msilx &! output
    ! & )

! call basaltweath_1D( &
    ! & nz,mfo,mg,si,hr,poro,z,dz,w,kfo,keqfo,mfoth,dmg,dsi,sat,dporodta,pro,mfoi,mfosupp  &! input
    ! & ,kco2,k1,k2,mgth,sith,tora,v,tol,zrxn,it,nax,cx,c2x,so4x,ca,pco2i,mgi,sii,mvfo,nflx &! input
    ! & ,iter,error,dt,flgback &! inout
    ! & ,mgx,six,prox,co2,hco3,co3,dic,mfox,omega_fo,flx_fo,flx_mg,flx_si &! output
    ! & )
    
call silicate_dis_1D( &
    & nz,mfo,mab,na,mg,si,hr,poro,z,dz,w,kfo,kab,keqfo,keqab,mfoth,mabth,dmg,dsi,dna,sat,dporodta,pro,mfoi,mabi,mfosupp,mabsupp  &! input
    & ,kco2,k1,k2,mgth,sith,nath,tora,v,tol,zrxn,it,cx,c2x,so4x,ca,pco2i,mgi,sii,nai,mvfo,mvab,nflx &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,prox,co2,hco3,co3,dic,mfox,mabx,omega_fo,omega_ab,flx_fo,flx_mg,flx_si,flx_ab,flx_na &! output
    & )


dporodtg = 0d0
dporodta = 0d0
#ifdef poroevol    
poro = poroi + (msili-msilx)*(100.07d0 -  0.5d0*99.52d0)*1d-6  &
    & +(msi-msx)*(23.94d0*(1d0-swoxall)+(23.94d0-20.82d0)*swoxall)*1d-6
v = qin/poro/sat
torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
deff = torg*dgas + tora*daq
dporodtg = ( &
    & (ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3) &
    & -(ucv*porox*(1.0d0-sat)*1d3+porox*sat*kho*1d3) &
    & )/dt
dporodta = (poro*sat-porox*sat)/dt/(poro*sat)
hr = hri
#ifdef surfevol1 
hr = hri*((1d0-poro)/(1d0-poroi))**(2d0/3d0)
#endif 
#ifdef surfevol2 
hr = hri*(poro/poroi)**(2d0/3d0)
#endif 
#endif 

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

if (iter2 > 300) then
    dt = dt/10d0
    if (dt==0d0) then 
        print *. 'dt==0d0; stop'
        stop
    endif 
end if

enddo     
! ######################## end of pH iteration #######################
#endif

#endif  
!  endif of ifndef pyweath
#ifdef display
print *
print *,'-=-=-=-=-=-= o2 & pyrite -=-=-=-=-=-=-='
print *,'o2:', (po2x(iz),iz=1,nz, nz/5)
print *,'fe2:', (cx(iz),iz=1,nz, nz/5)
print *,'py:', (msx(iz),iz=1,nz, nz/5)
print *, 'fe3:', (c2x(iz),iz=1,nz, nz/5)
print *, 'so4:', (so4x(iz),iz=1,nz, nz/5)
print *
print *,'-=-=-=-=-=-= Na & albite -=-=-=-=-=-=-='
print *, 'na:', (nax(iz),iz=1,nz, nz/5)
print *, 'sil:', (msilx(iz),iz=1,nz, nz/5)
print *
print *,'-=-=-=-=-=-= Mg, Si, Na, Fo, Ab -=-=-=-=-=-=-='
print *, 'mg:', (mgx(iz),iz=1,nz, nz/5)
print *, 'si:', (six(iz),iz=1,nz, nz/5)
print *, 'si:', (nax(iz),iz=1,nz, nz/5)
print *, 'fo:', (mfox(iz),iz=1,nz, nz/5)
print *, 'ab:', (mabx(iz),iz=1,nz, nz/5)
print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
print *
print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
print *
#endif 

! stop

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

        o2tflx = o2tflx +( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
            & +dporodtg(iz)*po2x(iz)   &
            & )*dz

        diflx = diflx  + ( &
            & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
            & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2(iz+1)+po2i-2.0d0*po2(iz))/(dz**2.0d0) &
            & )*dz

        advflx = advflx +( &
            & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz*(1.0d0-swex) &
            & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*(swex) &
            & ) *dz

        feoxflx = feoxflx +( &
            & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz) &
            & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
            & )*dz

        pyoxflx = pyoxflx  + ( &
            & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2(iz)**(0.50d0),(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0)))) &
            & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**(0.50d0),(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0)))) &
            & )*dz

        respflx = respflx  + ( &
            & +(1.0d0-frex)*swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2), &
            & (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & +(frex)*swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2), &
            & (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & )*dz

    else if (iz == nz) then

        o2tflx = o2tflx  + ( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
            & +dporodtg(iz)*po2x(iz)   &
            & )*dz

        diflx = diflx  + ( &
            & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
            & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2(iz-1)-1.0d0*po2(iz))/(dz**2.0d0) &
            & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
            & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
            & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
            & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
            & -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
            & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
            & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
            & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2(iz)-po2(iz-1))/(dz**2.0d0) &
            & )*dz

        advflx = advflx  + ( &
            & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
            & +(swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz &
            & )*dz

        feoxflx = feoxflx  + ( &
            & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
            & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
            & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
            & )*dz

        pyoxflx = pyoxflx  + ( &
            & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6 &
            & *merge(0d0,msx(iz)*po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(msx(iz)*po2x(iz)**(0.50d0))) &
            & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6 &
            & *merge(0d0,ms(iz)*po2(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(ms(iz)*po2(iz)**(0.50d0))) &
            & )*dz

        respflx = respflx  + ( &
            & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
            & +(frex)*swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2), &
            & (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & )*dz

    else

        o2tflx = o2tflx +  ( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
            & +dporodtg(iz)*po2x(iz)   &
            & )*dz

        diflx = diflx   + ( &
            & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
            & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
            & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)+(poro(iz)*sat(iz)*kho*tora(iz)*daq &
            & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
            & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2(iz+1)+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0) &
            & -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
            & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2(iz)-po2(iz-1))/(dz**2.0d0) &
            & )*dz

        advflx = advflx  + ( &
            & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz*swad &
            & +(swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad &
            & )*dz

        feoxflx = feoxflx  + ( &
            & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
            & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
            & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
            & )*dz

        pyoxflx = pyoxflx  + ( &
            & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
            & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2(iz)**(0.50d0))) &
            & )*dz

        respflx = respflx  +( &
            & +(1.0d0-frex)*swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2), &
            & (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & +(frex)*swbr*vmax*merge(0d0,po2x(iz)/(po2x(iz)+mo2), &
            & (po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & )*dz

    end if 

end do 

o2flxsum = o2tflx + advflx + diflx + pyoxflx + feoxflx + respflx 

do iz = 1, nz  !================================

    if (iz/=nz) then

        pytflx = pytflx + ( &
            & (msx(iz)-ms(iz))/dt &
            & )*dz

        pyadvflx = pyadvflx + ( &
            & -w*(msx(iz+1)-msx(iz))/dz*(1.0d0-swex)  &
            & -w*(ms(iz+1)-ms(iz))/dz*swex/dt*dt2*swpe &
            & )*dz

        pyox1flx = pyox1flx + ( &
            & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &
            & + frex*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0)) &
            & )*dz

        pyox2flx = pyox2flx + ( &
            & + (1.0d0-frex)*merge(0.0d0 &
            & ,koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
            & + frex*merge(0.0d0, &
            & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0))) &
            & )*dz

    else if (iz==nz) then

        pytflx = pytflx + ( &
            & (msx(iz)-ms(iz))/dt &
            & )*dz

        pyadvflx = pyadvflx + (  &
            & -w*(msi-msx(iz))/dz*(1.0d0-swex) &
            & -w*(msi-ms(iz))/dz*swex*dt2/dt*swpe &
            & )*dz

        pyox1flx = pyox1flx + ( &
            & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &
            & + frex*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0)) &
            & )*dz

        pyox2flx = pyox2flx + ( &
            & + (1.0d0-frex)*merge(0.0d0, &
            & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
            & + frex*merge(0.0d0, &
            & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0))) &
            & ) *dz
    end if 

end do  !================================

pyflxsum = pytflx + pyadvflx + pyox1flx + pyox2flx 

do iz = 1, nz

    if (.not.((iz == 1).or.(iz==nz))) then

        fetflx(1) = fetflx(1) + ( &
            & (cx(iz)-c(iz))/dt  &
            & +dporodta(iz)*cx(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        fediflx(1) = fediflx(1) + ( &
            & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))/(dz**2d0) &
            & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
            & +swex*(-dfe2*tora(iz)*(c(iz+1)+c(iz-1)-2d0*c(iz))/(dz**2d0) &
            & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        feadvflx(1) = feadvflx(1) + ( &
            & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
            & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        feox2flx(1) = feox2flx(1) + ( &
            & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
            & + frex*koxa(iz)*c(iz)*po2(iz) &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy1flx(1) = fepy1flx(1) + ( &
            & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy2flx(1) = fepy2flx(1) + ( &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)* &
            & (1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)* &
            & (1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz == 1) then

        fetflx(1) = fetflx(1) + ( &
            & (cx(iz)-c(iz))/dt  &
            & +dporodta(iz)*cx(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        fediflx(1) = fediflx(1) + ( &
            & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))/(dz**2d0)) &
            & +swex*(-dfe2*tora(iz)*(c(iz+1)+ci-2d0*c(iz))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        feadvflx(1) = feadvflx(1) + ( &
            & + v(iz)*(cx(iz)-ci)/dz*(1.0d0-swex) &
            & + v(iz)*(c(iz)-ci)/dz*swex &
            & )*dz*poro(iz)*sat(iz)*1d3

        feox2flx(1) = feox2flx(1) + ( &
            & + koxa(iz)*cx(iz)*po2x(iz)*(1.0d0-frex) &
            & + koxa(iz)*c(iz)*po2(iz)*frex &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy1flx(1) = fepy1flx(1) + ( &
            & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy2flx(1) = fepy2flx(1) + ( &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)* &
            & (1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)* &
            & (1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz ==nz ) then

        fetflx(1) = fetflx(1) + ( &
            & (cx(iz)-c(iz))/dt  &
            & +dporodta(iz)*cx(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        fediflx(1) = fediflx(1) + ( &
            & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))/(dz**2d0) &
            & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
            & +swex*(-dfe2*tora(iz)*(c(iz-1)-1d0*c(iz))/(dz**2d0) &
            & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        feadvflx(1) = feadvflx(1) + ( &
            & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
            & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        feox2flx(1) = feox2flx(1) + ( &
            & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
            & + frex*koxa(iz)*c(iz)*po2(iz) &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy1flx(1) = fepy1flx(1) + ( &
            & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy2flx(1) = fepy2flx(1) + ( &
            &  -(1.0d0-frex)*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)* &
            & (1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)* &
            & (1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    end if 

end do  ! ==============================

! feflxsum(1) = fetflx(1) + fediflx(1) +  feadvflx(1) + feox2flx(1)+ fepy1flx(1) + fepy2flx(1)

fetflx(1) = sum(flx_t(:)*dz*poro(:)*sat(:)*1d3)
fediflx(1) =sum(flx_dif(:)*dz*poro(:)*sat(:)*1d3) 
feadvflx(1) =sum(flx_adv(:)*dz*poro(:)*sat(:)*1d3) 
feox2flx(1) =sum(flx_oxaq(:)*dz*poro(:)*sat(:)*1d3) 
fepy1flx(1) =sum(flx_oxpy1(:)*dz*poro(:)*sat(:)*1d3) 
fepy2flx(1) =sum(flx_oxpy2(:)*dz*poro(:)*sat(:)*1d3) 
feflxsum(1) =sum(flx_rem(:)*dz*poro(:)*sat(:)*1d3) 

do iz = 1, nz

    if (.not.((iz == 1).or.(iz==nz))) then

        fetflx(2) = fetflx(2) + ( &
            & (c2x(iz)-c2(iz))/dt  &
            & +dporodta(iz)*c2x(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        fediflx(2) = fediflx(2) + ( &
            & +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))/(dz**2d0) &
            & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0)) &
            & +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2(iz-1)-2d0*c2(iz))/(dz**2d0) &
            & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2(iz)-c2(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(Iz)*1d3

        feadvflx(2) = feadvflx(2) + ( &
            & + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex) &
            & + v(iz)*(c2(iz)-c2(iz-1))/dz*swex &
            & )*dz*poro(iz)*sat(Iz)*1d3

        feox2flx(2) = feox2flx(2) + ( &
            & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
            & - frex*koxa(iz)*c(iz)*po2(iz) &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy2flx(2) = fepy2flx(2) + ( &
            & +(1.0d0-frex)*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & +frex*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3 &
            & *merge(0.0d0,1.0d0,c2x(iz)<c2th)

    else if (iz == 1) then

        fetflx(2) = fetflx(2) + ( &
            & (c2x(iz)-c2(iz))/dt  &
            & +dporodta(iz)*c2x(iz)   &
            & )*dz*poro(iz)*sat(Iz)*1d3

        fediflx(2) = fediflx(2) + ( &
            & +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))/(dz**2d0)) &
            & +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2i-2d0*c2(iz))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(Iz)*1d3

        feadvflx(2) = feadvflx(2) + ( &
            & + v(iz)*(c2x(iz)-c2i)/dz*(1.0d0-swex) &
            & + v(iz)*(c2(iz)-c2i)/dz*swex &
            & )*dz*poro(iz)*sat(iz)*1d3

        feox2flx(2) = feox2flx(2) + ( &
            & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
            & - frex*koxa(iz)*c(iz)*po2(iz) &
            & )*dz*poro(iz)*sat(Iz)*1d3

        fepy2flx(2) = fepy2flx(2) + ( &
            & +(1.0d0-frex)*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & +frex*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz ==nz) then

        fetflx(2) = fetflx(2) + ( &
            & (c2x(iz)-c2(iz))/dt  &
            & +dporodta(iz)*c2x(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        fediflx(2) = fediflx(2) + ( &
            & +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))/(dz**2d0) &
            & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0)) &
            & +swex*(-dfe3*tora(iz)*(c2(iz-1)-1d0*c2(iz))/(dz**2d0) &
            & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2(iz)-c2(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(Iz)*1d3

        feadvflx(2) = feadvflx(2) + ( &
            & + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex) &
            & + v(iz)*(c2(iz)-c2(iz-1))/dz*swex &
            & )*dz*poro(iz)*sat(Iz)*1d3

        feox2flx(2) = feox2flx(2) + ( &
            & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
            & - frex*koxa(iz)*c(iz)*po2(iz) &
            & )*dz*poro(iz)*sat(iz)*1d3

        fepy2flx(2) = fepy2flx(2) + ( &
            & +(1.0d0-frex)*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & +frex*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3 &
            & *merge(0.0d0,1.0d0,c2x(iz)<c2th)

    end if 

end do

feflxsum(2) = fetflx(2) + fediflx(2) +  feadvflx(2) + feox2flx(2)+ fepy2flx(2)

do iz = 1, nz  !================================

    if (iz/=nz) then

        siltflx = siltflx  + ((msilx(iz)-msil(iz))/dt)*dz

        siladvflx = siladvflx + ( &
            & -w*(msilx(iz+1)-msilx(iz))/dz*(1.0d0-swex)  &
            & -w*(msil(iz+1)-msil(iz))/dz*swex/dt*dt2*swpe &
            & )*dz

        sildisflx = sildisflx + ( &
            & + (1.0d0-frex)* &
            & ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
            & + frex*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0) &
            & )*dz

    else if (iz==nz) then

        siltflx = siltflx + ((msilx(iz)-msil(iz))/dt)*dz

        siladvflx = siladvflx + ( &
            & -w*(msili-msilx(iz))/dz*(1.0d0-swex) &
            & -w*(msili-msil(iz))/dz*swex*dt2/dt*swpe &
            & )*dz

        sildisflx = sildisflx + ( &
            & + (1.0d0-frex)* &
            & ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
            & + frex*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0) &
            & )*dz

    end if 

end do  !================================

silflxsum = siltflx + siladvflx + sildisflx 

do iz = 1, nz

    if (.not.((iz == 1).or.(iz==nz))) then

        natflx = natflx + ( &
            & (nax(iz)-na(iz))/dt &
            & +dporodta(iz)*nax(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        naadvflx = naadvflx + ( &
            & + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz &
            & + swex*v(iz)*(na(iz)-na(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadiflx = nadiflx + ( &
            & +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nax(iz-1)-2d0*nax(iz))/(dz**2d0) &
            & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(nax(iz)-nax(iz-1))/(dz**2d0)) &
            & +swex*(-dna*tora(iz)*(na(iz+1)+na(iz-1)-2d0*na(iz))/(dz**2d0) &
            & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(na(iz)-na(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadisflx = nadisflx + ( &
            & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
            & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz == 1) then

        natflx = natflx + ( &
            & (nax(iz)-na(iz))/dt  &
            & +dporodta(iz)*nax(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3

        naadvflx = naadvflx + ( &
            & + v(iz)*(nax(iz)-nai)/dz*(1.0d0-swex) &
            & + v(iz)*(na(iz)-nai)/dz*swex &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadiflx = nadiflx + ( &
            & +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nai-2d0*nax(iz))/(dz**2d0)) &
            & +swex*(-dna*tora(iz)*(na(iz+1)+nai-2d0*na(iz))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadisflx = nadisflx + ( &
            & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
            & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz == nz) then

        natflx = natflx + ( &
            & (nax(iz)-na(iz))/dt  &
            & +dporodta(iz)*nax(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3   ! commented out (is this necessary?)

        naadvflx = naadvflx + ( &
            & + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz &
            & + swex*v(iz)*(na(iz)-na(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadiflx = nadiflx + ( &
            & +(1d0-swex)*(-dna*tora(iz)*(nax(iz-1)-1d0*nax(iz))/(dz**2d0) &
            & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(nax(iz)-nax(iz-1))/(dz**2d0)) &
            & +swex*(-dna*tora(iz)*(na(iz-1)-1d0*na(iz))/(dz**2d0) &
            & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(na(iz)-na(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        nadisflx = nadisflx + ( &
            & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
            & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    end if ! =========================================
enddo 

naflxsum = natflx + nadiflx + naadvflx + nadisflx 

do iz = 1, nz

    if (.not.((iz == 1).or.(iz==nz))) then

        so4tflx = so4tflx+ ( &
            & (so4x(iz)-so4(iz))/dt &
            & +dporodta(iz)*so4x(iz)   &
            & )*dz*poro(iz)*sat(iz)*1d3


        so4diflx = so4diflx+ ( &
            & +(-dso4*tora(iz)*(so4x(iz+1)+so4x(iz-1)-2d0*so4x(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(so4x(iz)-so4x(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4advflx = so4advflx+ ( &
            & + v(iz)*(so4x(iz)-so4x(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4disflx = so4disflx+ ( &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz == 1) then

        so4tflx = so4tflx+ ( &
            & (so4x(iz)-so4(iz))/dt  &
            & +dporodta(iz)*so4x(iz)  &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4diflx = so4diflx+ ( &
            & +(-dso4*tora(iz)*(so4x(iz+1)+so4i-2d0*so4x(iz))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4advflx = so4advflx+ ( &
            & + v(iz)*(so4x(iz)-so4i)/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4disflx = so4disflx+ ( &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

    else if (iz == nz) then

        so4tflx = so4tflx + ( &
            & (so4x(iz)-so4(iz))/dt  &
            & +dporodta(iz)*so4x(iz)  &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4diflx = so4diflx + ( &
            & +(-dso4*tora(iz)*(so4x(iz-1)-1d0*so4x(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(so4x(iz)-so4x(iz-1))/(dz**2d0)) &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4advflx = so4advflx + ( &
            & + v(iz)*(so4x(iz)-so4x(iz-1))/dz &
            & )*dz*poro(iz)*sat(iz)*1d3

        so4disflx = so4disflx + ( &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )*dz*poro(iz)*sat(iz)*1d3

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
mg = mgx
si = six
mfo = mfox
mab = mabx
pro = prox

if (time>=rectime(irec+1)) then
    write(chr,'(i3.3)') irec+1
    open (28, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(rate)-'//chr//'.txt',  &
        & status='replace')
    open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-'//chr//'.txt',  &
        & status='replace')
    open (29, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(nom)-'//chr//'.txt',  &
        & status='replace')
    open(30, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-'//chr//'.txt',  &
        & status='replace')


    do iz = 1, Nz
        write (28,*) z(iz), &
            & koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*po2x(iz)**0.50d0*merge(0.0d0,1.0d0,po2x(iz)<po2th), &
            & merge(0.0d0, &
            & + koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth) &
            & *merge(0.0d0,1.0d0,c2x(iz)<c2th), &
            & poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz) &
            & *merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth) &
            & , swbr*vmax*po2x(iz)/(po2x(iz)+mo2) &
            & ,ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
            & ,time
        write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),mg(iz),si(iz),msil(iz),mfo(iz),pro(iz) &
            & ,silsat(iz), omega_fo(iz),time
        write (29,*) z(iz),po2(iz)/maxval(po2(:)),c(iz)/maxval(c(:)),ms(iz)/maxval(ms(:)),c2(iz)/maxval(c2(:)) &
            & , so4(iz)/maxval(so4(:)), na(iz)/maxval(na(:)), msil(iz)/maxval(msil(:)), pro(iz)/maxval(pro(:)) &
            & , silsat(iz),co2(iz),hco3(iz),co3(iz),dic(iz),time
        write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
    end do
    irec=irec+1

    close(28)
    close(22)
    close(29)
    close(30)

    write(65,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx,o2flxsum 
    write(67,*) time, fetflx(1),fediflx(1), feadvflx(1),feox2flx(1), fepy1flx(1), fepy2flx(1),feflxsum(1)
    write(68,*) time, fetflx(2),fediflx(2), feadvflx(2),feox2flx(2), fepy1flx(2), fepy2flx(2),feflxsum(2) 
    write(69,*) time, pytflx, pyadvflx, pyox1flx, pyox2flx,pyflxsum        
    write(56,*) time, siltflx, siladvflx, sildisflx,silflxsum        
    write(55,*) time, natflx, naadvflx, nadiflx, nadisflx,naflxsum        
    write(57,*) time, so4tflx, so4advflx, so4diflx, so4disflx,so4flxsum
    
    do iflx = 1,nflx
        flx_mg(iflx,:) = flx_mg(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_si(iflx,:) = flx_si(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_na(iflx,:) = flx_na(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_fo(iflx,:) = flx_fo(iflx,:)*dz
        flx_ab(iflx,:) = flx_ab(iflx,:)*dz
    enddo 
    
    write(58,*) time,(sum(flx_mg(iflx,:)),iflx=1,nflx)
    write(59,*) time,(sum(flx_si(iflx,:)),iflx=1,nflx)
    write(60,*) time,(sum(flx_fo(iflx,:)),iflx=1,nflx)
    write(61,*) time,(sum(flx_na(iflx,:)),iflx=1,nflx)
    write(62,*) time,(sum(flx_ab(iflx,:)),iflx=1,nflx)

    write(95,*) time, o2in, o2out, po2i,msi,zrxn, zsat

end if


write(71,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx,o2flxsum
write(73,*) time, fetflx(1),fediflx(1), feadvflx(1), feox2flx(1), fepy1flx(1), fepy2flx(1),feflxsum(1)
write(74,*) time, fetflx(2),fediflx(2), feadvflx(2), feox2flx(2), fepy1flx(2), fepy2flx(2),feflxsum(2)
write(75,*) time, pytflx, pyadvflx, pyox1flx, pyox2flx,pyflxsum      
write(77,*) time, siltflx, siladvflx, sildisflx,silflxsum      
write(76,*) time, natflx, naadvflx, nadiflx, nadisflx,naflxsum      
write(78,*) time, so4tflx, so4advflx, so4diflx, so4disflx,so4flxsum      

write(97,*) time, o2in, o2out, po2i,msi,zrxn, zsat

it = it + 1
time = time + dt

end do

write(chr,'(i3.3)') irec+1
open (58, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(rate)-'//chr//'.txt', &
    & status='replace')
open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-'//chr//'.txt', &
    & status='replace')
open (29, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(nom)-'//chr//'.txt', &
    & status='replace')

open(30, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-'//chr//'.txt',  &
    & status='replace')

do iz = 1, Nz
    write (58,*) z(iz), &
        & koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*po2x(iz)**0.50d0*merge(0.0d0,1.0d0,po2x(iz)<po2th), &
        & merge(0.0d0, &
        & + koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth) &
        & *merge(0.0d0,1.0d0,c2x(iz)<c2th), &
        & poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)*merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth) &
        & , swbr*vmax*po2x(iz)/(po2x(iz)+mo2) &
        & ,ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
        & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
        & ,time
    write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),mg(iz),si(iz),msil(iz),mfo(iz),pro(iz) &
        & ,silsat(iz), omega_fo(iz),time
    write (29,*) z(iz),po2(iz)/maxval(po2(:)),c(iz)/maxval(c(:)),ms(iz)/maxval(ms(:)),c2(iz)/maxval(c2(:)) &
        & ,so4(iz)/maxval(so4(:)),na(iz)/maxval(na(:)),msil(iz)/maxval(msil(:)),pro(iz)/maxval(pro(:)) &
        & ,silsat(iz), co2(iz),hco3(iz),co3(iz),dic(iz),time
    write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
end do


write(65,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx,o2flxsum

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
    if ( zpy(1)==0d0 .and. ms(iz)>=0.1d0*msi+0.9d0*ms(1)) zpy(1) = z(iz)    
    if ( zpy(2)==0d0 .and. ms(iz)>=0.5d0*msi+0.5d0*ms(1)) zpy(2) = z(iz)    
    if ( zpy(3)==0d0 .and. ms(iz)>=0.9d0*msi+0.1d0*ms(1)) zpy(3) = z(iz)    
    if ( zab(1)==0d0 .and. msil(iz)>=0.1d0*msili+0.9d0*msil(1)) zab(1) = z(iz)    
    if ( zab(2)==0d0 .and. msil(iz)>=0.5d0*msili+0.5d0*msil(1)) zab(2) = z(iz)    
    if ( zab(3)==0d0 .and. msil(iz)>=0.9d0*msili+0.1d0*msil(1)) zab(3) = z(iz)    
enddo

fname=trim(adjustl(workdir))//'sense'//trim(adjustl(base))//'.txt'
open(20,file=trim(adjustl(fname)),action='write',status='unknown',access='append')
write(20,*) qin,zsat,zpy(1:3),zab(1:3)
close(20)
#ifdef pyweath 
fname=trim(adjustl(workdir))//'sense-py'//trim(adjustl(base))//'.txt'
open(400,file=trim(adjustl(fname)),action='write',status='unknown',access='append')
write(400,*) qin,zsat,zpy(1:3)
close(400)
#endif      
#ifdef silweath 
fname=trim(adjustl(workdir))//'sense-ab'//trim(adjustl(base))//'.txt'
open(401,file=trim(adjustl(fname)),action='write',status='unknown',access='append')
write(401,*) qin,zsat,zab(1:3)
close(401)
#endif

endprogram o2profile

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine pyweath_1D( &
    & nz,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
    & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daq,dgas,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
    & ,swoxa,swoxall,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,cx,c2x,po2x,msx,flx_t,flx_adv,flx_dif,flx_oxaq,flx_oxpy1,flx_oxpy2,flx_rem &! output
    ) 
    
implicit none 

integer,intent(in)::nz,nsp
real(kind=8),intent(in)::ci,c2i,po2i,msi,po2th,dz,w,msth,dfe2,dfe3,kho,dt2,cth,c2th,stoxa,daq,dgas &
    & ,swbr,mo2,stoxs,tol,zrxn,swoxa,swoxall,ucv,vmax
real(kind=8),dimension(nz),intent(in)::c,c2,po2,ms,hr,poro,z,koxs2,koxs,sat,dporodta,dporodtg &
    & ,tora,torg,v,koxa
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(out)::cx,c2x,po2x,msx,flx_t,flx_adv,flx_dif,flx_oxaq,flx_oxpy1,flx_oxpy2 &
    & ,flx_rem
integer,intent(inout)::iter,it
real(kind=8),intent(inout)::error,dt

integer iz,row,nmx,ie,ie2

real(kind=8)::swex = 0.0d0 ! switch for explicit
real(kind=8)::frex = 0.0d0 ! fraction of explicit
real(kind=8)::swpe = 0.0d0 ! physical erosion
real(kind=8)::swad = 1.0d0! 1.0 when advection included 0.0d0 when not
real(kind=8),parameter::infinity = huge(0d0)

real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV


nmx = nsp*nz

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0

    if (it == 0 .and. iter == 0) then
        cx(:) = 1.0d2
        c2x(:) = 1.0d2
        ! so4x(:) = 1.0d2
        ! nax(:) = 1.0d2
    end if

    !       print*,msilx
    !       stop

    if (any(abs(cx)>infinity).or.any(abs(c2x)>infinity)  &
        & .or.any(abs(po2x)>infinity) .or.any(abs(msx)>infinity) ) then 
        print *, '*** ERROR before starting the newton iteration'
        print *, '*** INFINITY detected '       
        print *,any(abs(cx)>infinity),any(abs(c2x)>infinity) ,any(abs(po2x)>infinity) ,any(abs(msx)>infinity)  
        print *, '*** INFINITY returned to the arbitrary maximum value'
        do iz=1,nz
            if ((abs(cx(iz))>infinity))cx(iz)=1d2
            if ((abs(c2x(iz))>infinity))c2x(iz)=1d2 
            if ((abs(po2x(iz))>infinity)) po2x(iz)=po2i 
            if ((abs(msx(iz))>infinity)) msx(iz)=msi 
        enddo
    endif 
    if (any(isnan(cx)).or.any(isnan(c2x)).or.any(isnan(po2x)).or.any(isnan(msx)) ) then 
        print *, '*** ERROR before starting the newton iteration'
        print *, '*** NAN detected '       
        print *,any(isnan(cx)),any(isnan(c2x)) ,any(isnan(po2x)) ,any(isnan(msx)) 
        print *, '*** NAN returned to the arbitrary maximum value'
        do iz=1,nz
            if (isnan(cx(iz)))cx(iz)=1d2
            if (isnan(c2x(iz)))c2x(iz)=1d2 
            if (isnan(po2x(iz))) po2x(iz)=po2i 
            if (isnan(msx(iz))) msx(iz)=msi 
        enddo
    endif 

    do iz = 1, nz  !================================

        row = nsp*(iz-1)+1

        if (iz/=nz) then

            amx(row,row) = (1.0d0/dt    &
                & + w/dz *(1.0d0-swex)    &
                & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6 &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &
                & + (1.0d0-frex)*merge(0.0d0, &
                & + koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                &  cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & ) &
                & * merge(1.0d0,msx(iz),msx(iz)<msth)

            amx(row,row+nsp) = (-w/dz)*(1.0d0-swex)*merge(1.0d0,msx(iz+1),msx(iz)<msth)

            ymx(row) = ( &
                & (msx(iz)-ms(iz))/dt  &
                & -w*(msx(iz+1)-msx(iz))/dz*(1.0d0-swex)  &
                & -w*(ms(iz+1)-ms(iz))/dz*swex/dt*dt2*swpe  &
                & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))  &
                & + frex*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0)) &
                & + (1.0d0-frex)*merge(0.0d0, &
                & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & + frex*merge(0.0d0,koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0))) &
                & ) &
                & *merge(0.0d0,1d0,msx(iz)<msth)

        else if (iz==nz) then

            amx(row,row) = (1.0d0/dt  &
                & + w/dz*(1.0d0-swex)  &
                & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))  &
                & + (1.0d0-frex)*merge(0.0d0,  &
                & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))  &
                & )  &
                & *merge(1.0d0,msx(iz),msx(iz)<msth)

            ymx(row) = (   &
                & (msx(iz)-ms(iz))/dt  & 
                & -w*(msi-msx(iz))/dz*(1.0d0-swex) &
                & -w*(msi-ms(iz))/dz*swex*dt2/dt*swpe &
                & + (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &
                & + frex*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))   &          
                & + (1.0d0-frex)*merge(0.0d0, &
                & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & + frex*merge(0.0d0, &
                & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0))) &
                & ) &
                & *merge(0.0d0,1d0,msx(iz)<msth)
        end if 

        amx(row,row + 2 ) = ( & 
            & + (1.0d0-frex)*merge(0.0d0, &
            & koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*0.50d0*(po2x(iz)**(-0.50d0)), &
            & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
            & ) &
            & *po2x(iz)&
            & *merge(0.0d0,1d0,msx(iz)<msth) 

        amx(row,row  + 1 ) = (  &
            & + (1.0d0-frex)*merge(0.0d0,  &
            & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-1.40d0), &
            & cx(iz)<cth .or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
            & ) &
            & *cx(iz)&
            & *merge(0.0d0,1d0,msx(iz)<msth)

        amx(row,row  + 3 ) = ( &
            & + (1.0d0-frex)*merge(0.0d0,  &
            & koxs2(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),  &
            & (c2x(iz)<c2th).or.(cx(iz)<cth).or.c2x(iz)<c2th.or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0))) &
            & ) &
            & *c2x(iz)  &
            & *merge(0.0d0,1d0,msx(iz)<msth)


    end do  !================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Fe++   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    flx_t = 0d0
    flx_adv = 0d0
    flx_dif = 0d0
    flx_oxaq = 0d0
    flx_oxpy1 = 0d0
    flx_oxpy2 = 0d0
    flx_rem = 0d0

    do iz = 1, nz

        row = nsp*(iz-1)+2

        if (.not.((iz == 1).or.(iz==nz))) then

            amx(row,row) = (  &
                & 1.0d0/dt   &
                & +dporodta(iz)   &
                & +(1d0-swex)*(-dfe2*tora(iz)*(-2d0)/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0))  &
                & + v(iz)/dz*(1.0d0-swex)  &
                & + (1.0d0-frex)*koxa(iz)*po2x(iz) &
                & + (1.0d0-frex)*merge(0.0d0,  &
                & -(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*  &
                & c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))  &
                & *1d-3  &
                &  )  &
                & *merge(1.0d0,cx(iz),cx(iz)<cth)

            amx(row,row-nsp) = (  &
                & +(1d0-swex)*(-dfe2*tora(iz)*(1d0)/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - (1.0d0-swex)*v(iz)/dz &
                & ) &
                & *cx(iz-1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            amx(row,row+nsp) = ( &
                & +(1d0-swex)*(-dfe2*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *cx(iz+1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt &
                & +dporodta(iz) *cx(iz)+(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz+1)+c(iz-1)-2d0*c(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
                & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
                & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & + frex*koxa(iz)*c(iz)*po2(iz)  &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3  &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)  &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3  &
                & -(1.0d0-frex)*merge(0.0d0,  &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0  &
                & ,(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0),  &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & )  &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_t(iz) = ( &
                & (cx(iz)-c(iz))/dt  &
                & +dporodta(iz) *cx(iz)  &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_dif(iz) = (  &
                & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz+1)+c(iz-1)-2d0*c(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
                & )  &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_adv(iz) = (  &
                & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxaq(iz) = (  & 
                & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & + frex*koxa(iz)*c(iz)*po2(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy1(iz) = (  &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy2(iz) = ( &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        else if (iz == 1) then

            amx(row,row) = (1.0d0/dt +dporodta(iz)  &
                & +(1d0-swex)*(-dfe2*tora(iz)*(-2d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & + koxa(iz)*po2x(iz)*(1.0d0-frex) &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)&
                & *c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,cx(iz),cx(iz)<c2th)

            amx(row,row+nsp) = ( &
                & +(1d0-swex)*(-dfe2*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *cx(iz+1)  &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt  &
                & +dporodta(iz)*cx(iz)  &
                & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz+1)+ci-2d0*c(iz))/(dz**2d0)) &
                & + v(iz)*(cx(iz)-ci)/dz*(1.0d0-swex) &
                & + v(iz)*(c(iz)-ci)/dz*swex &
                & + koxa(iz)*cx(iz)*po2x(iz)*(1.0d0-frex) &
                & + koxa(iz)*c(iz)*po2(iz)*frex &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_t(iz) = ( &
                & (cx(iz)-c(iz))/dt +dporodta(iz)*cx(iz)  &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_dif(iz) = ( &
                & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz+1)+ci-2d0*c(iz))/(dz**2d0)) &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_adv(iz) = ( &
                & + v(iz)*(cx(iz)-ci)/dz*(1.0d0-swex) &
                & + v(iz)*(c(iz)-ci)/dz*swex &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxaq(iz) = ( &
                & + koxa(iz)*cx(iz)*po2x(iz)*(1.0d0-frex)  &
                & + koxa(iz)*c(iz)*po2(iz)*frex &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy1(iz) = ( &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy2(iz) = (  &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        else if (iz == nz) then

            amx(row,row) = (  &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dfe2*tora(iz)*(-1d0)/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                & *(1d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & + (1.0d0-frex)*koxa(iz)*po2x(iz) &
                & + (1.0d0-frex)*merge(0.0d0, &
                & -(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)&
                & *c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,cx(iz),cx(iz)<cth)

            amx(row,row-nsp) = ( &
                & +(1d0-swex)*(-dfe2*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - (1.0d0-swex)*v(iz)/dz &
                & ) &
                & *cx(iz-1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt &
                & +dporodta(iz) *cx(iz) &
                & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz-1)-1d0*c(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
                & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
                & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & + frex*koxa(iz)*c(iz)*po2(iz) &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0 &
                & ,(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_t(iz) = ( &
                & (cx(iz)-c(iz))/dt  &
                & +dporodta(iz) *cx(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_dif(iz) = ( &
                & +(1d0-swex)*(-dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0)) &
                & +swex*(-dfe2*tora(iz)*(c(iz-1)-1d0*c(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c(iz)-c(iz-1))/(dz**2d0)) &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_adv(iz) = ( &
                & + (1.0d0-swex)*v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + swex*v(iz)*(c(iz)-c(iz-1))/dz &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxaq(iz) = (  &
                & + (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & + frex*koxa(iz)*c(iz)*po2(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy1(iz) = ( &
                & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall) &
                & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_oxpy2(iz) = ( &
                & -(1.0d0-frex)*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & -frex*merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*(1d0-swoxall)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        end if 

        amx(row,row+1) = (  &
            & + (1.0d0-frex)*merge(0.0d0,koxa(iz)*cx(iz),po2x(iz)<po2th) &
            & - (1.0d0-frex)*merge(0.0d0, &
            & koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*0.50d0*(1d0-swoxall)*po2x(iz)**(-0.50d0), &
            & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))*1d-3  &
            & ) &
            & *po2x(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        amx(row,row+2) = ( &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(1d0-swoxall)* &
            & (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
            & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *c2x(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        amx(row,row  - 1) = (     &
            & - (1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or. isnan(po2x(iz)**(0.50d0)))*1d-3  &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *msx(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)


        if (isnan(ymx(row))) then 
            print*,'**** NAN found for Fe++ at iz = ',iz ,iter
        endif 

    end do  ! ==============================

    flx_rem = flx_t + flx_adv + flx_dif + flx_oxaq + flx_oxpy1  + flx_oxpy2       

    do iz = 1, nz

        row = nsp*(iz-1)+4  !! fe+++

        if (.not.((iz == 1).or.(iz==nz))) then

            amx(row,row) = ( &
                & 1.0d0/dt & 
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dfe3*tora(iz)*(-2d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                & + v(iz)/dz *(1.0d0-swex) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row - nsp) = ( &
                & +(1d0-swex)*(-dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - v(iz)/dz*(1.0d0-swex) &
                & ) &
                & *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)

            amx(row,row+nsp) = ( &
                & +(1d0-swex)*(-dfe3*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)


            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt &
                & +dporodta(iz) *c2x(iz) &
                & +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))&
                & *(c2x(iz)-c2x(iz-1))/(dz**2d0)) &
                & +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2(iz-1)-2d0*c2(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))&
                & *(c2(iz)-c2(iz-1))/(dz**2d0)) &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex) &
                & + v(iz)*(c2(iz)-c2(iz-1))/dz*swex &
                & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & - frex*koxa(iz)*c(iz)*po2(iz) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & +frex*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)

        else if (iz == 1) then

            amx(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dfe3*tora(iz)*(-2d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row+nsp) = ( &
                & +(1d0-swex)*(-dfe3*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)

            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt  &
                & +dporodta(iz) *c2x(iz)+(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))/(dz**2d0)) &
                & +swex*(-dfe3*tora(iz)*(c2(iz+1)+c2i-2d0*c2(iz))/(dz**2d0)) &
                & + v(iz)*(c2x(iz)-c2i)/dz*(1.0d0-swex) &
                & + v(iz)*(c2(iz)-c2i)/dz*swex &
                & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & - frex*koxa(iz)*c(iz)*po2(iz) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & +frex*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0),&
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)  ! commented out 

        else if (iz ==nz) then

            amx(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz) &
                & +(1d0-swex)*(-dfe3*tora(iz)*(-1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                & + v(iz)/dz *(1.0d0-swex) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row - nsp) = ( &
                & +(1d0-swex)*(-dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - v(iz)/dz*(1.0d0-swex) &
                & ) &
                & *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)


            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt   &
                & +dporodta(iz) *c2x(iz) &
                & +(1d0-swex)*(-dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))&
                & *(c2x(iz)-c2x(iz-1))/(dz**2d0)) &
                & +swex*(-dfe3*tora(iz)*(c2(iz-1)-1d0*c2(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))&
                & *(c2(iz)-c2(iz-1))/(dz**2d0)) &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz*(1.0d0-swex) &
                & + v(iz)*(c2(iz)-c2(iz-1))/dz*swex &
                & - (1.0d0-frex)*koxa(iz)*cx(iz)*po2x(iz) &
                & - frex*koxa(iz)*c(iz)*po2(iz) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & +frex*merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
                & c(iz)<cth.or.c2(iz)<c2th.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)

        end if 

        amx(row,row-1) = ( &
            & - (1.0d0-frex)*koxa(iz)*cx(iz) &
            & ) &
            & *merge(0.0d0,po2x(iz),c2x(iz)<c2th)

        amx(row,row-2) = ( &
            & - (1.0d0-frex)*koxa(iz)*po2x(iz) &
            & +(1.0d0-frex)*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *merge(0.0d0,cx(iz),c2x(iz)<c2th)

        amx(row,row  - 3) = ( &
            & +(1.0d0-frex)*merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *merge(0.0d0,msx(iz),c2x(iz)<c2th)

        if (isnan(ymx(row))) then 
            print *, '***** NAN for Fe+++ at iz = ',iz , iter
        endif 

    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iz = 1, nz

        row = 3 + nsp*(iz-1)

        if (iz == 1) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +(1.0d0-frex)*2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz*(1.0d0-swex) &
                & +(1.0d0-frex)*merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz)<po2th) &
                & +(1.0d0-frex)*merge(0.d0 &
                & ,stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz)<po2th.or.isnan(msx(iz)*0.50d0*po2x(iz)**(-0.50d0))) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            if (isnan(amx(row,row))) then 
                print*,'nan in oxygen',iz
            endif

            amx(row,row+nsp) = ( &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            if (isnan(amx(row,row+nsp))) then 
                print *,'error in oxygen +',iz
            endif

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2(iz+1)+po2i-2.0d0*po2(iz))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz*(1.0d0-swex) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz*(swex) &
                & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz) &
                & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
                & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
                & *merge(0d0,po2(iz)**(0.50d0),(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0)))) &
                & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
                & *merge(0d0,po2(iz)**(0.50d0),(po2x(iz) <po2th).or.(isnan(po2(iz)**(0.50d0)))) &
                & +(1.0d0-frex)*swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & +(frex)*swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +(1.0d0-frex)*1.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(1.0d0)/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +(1.0d0-frex)*merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz) <po2th) &
                & +(1.0d0-frex)*merge(0.0d0 &
                & ,stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            if (isnan(amx(row,row))) then 
                print*,'nan in oxygen', iz
            endif

            amx(row,row-nsp) = ( &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            if (isnan(amx(row,row-nsp))) then 
                print*,'nan in oxygen - ',iz
            endif

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2(iz-1)-1.0d0*po2(iz))/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2(iz)-po2(iz-1))/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +(swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz &
                & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
                & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
                & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6 &
                & *merge(0d0,msx(iz)*po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(msx(iz)*po2x(iz)**(0.50d0))) &
                &  +(1.0d0-frex)*swbr*vmax & 
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz)  &
                & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6 &
                & *merge(0d0,ms(iz)*po2(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(ms(iz)*po2(iz)**(0.50d0))) &
                & +(frex)*swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2))))  &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt & 
                & +dporodtg(iz)  &
                & +(1.0d0-frex)*2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(1.0d0)/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +(1.0d0-frex)*merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz) <po2th) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
                & +(1.0d0-frex)*merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            if (isnan(amx(row,row))) then 
                print*,'nan in oxygen ',iz
            endif


            amx(row,row+nsp) = ( &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            if (isnan(amx(row,row+nsp))) then 
                print*,'nan in oxygen +',iz
            endif

            amx(row,row-nsp) = ( &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)  &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            if (isnan(amx(row,row-nsp))) then 
                print*,'nan in oxygen -',iz
            endif

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
                & +dporodtg(iz) *po2x(iz) &
                & -(1.0d0-frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -(1.0d0-frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +(1.0d0-swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz*swad &
                & -(frex)*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2(iz+1)+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0) &
                & -(frex)*1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(po2(iz)-po2(iz-1))/(dz**2.0d0) &
                & +(swex)*poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz*swad &
                & +(1.0d0-frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
                & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
                & +(1.0d0-frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
                & +(1.0d0-frex)*swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & +(frex)*stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*c(iz)*po2(iz) &
                & +(frex)*stoxs*koxs(iz)*poro(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
                & *merge(0d0,po2(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2(iz)**(0.50d0))) &
                & +(frex)*swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)

        end if 

        amx(row,row-1) = ( &
            & +(1.0d0-frex)*poro(iz)*sat(iz)*1d3*koxa(iz)*po2x(iz)*stoxa               &
            & ) &
            & *merge(0.0d0,cx(iz),po2x(iz)<po2th)

        if (isnan(amx(row,row-1))) then 
            print*,'nan in oxygen for fe2',iz
        endif

        amx(row,row  - 2) = ( &
            & (1.0d0-frex)*koxs(iz)*poro(iz)*hr(iz)*stoxs*23.94d0*1d-6 &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
            & ) &
            & *merge(0.0d0,msx(iz),po2x(iz)<po2th)

        if (isnan(amx(row,row - 2))) then 
            print*,'nan in oxygen for pyrite',iz
        endif

        if (isnan(ymx(row))) then 
            print*, '****** NAN found for po2 at iz =',iz,iter
        endif 

    end do 

    ymx=-1.0d0*ymx

    if (any(isnan(amx)).or.any(isnan(ymx))) then 
        print*,'**** error in mtx py ***********'
        print*,'any(isnan(amx)),any(isnan(ymx))'
        print*,any(isnan(amx)),any(isnan(ymx))

        if (any(isnan(ymx))) then 
            open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_pre.txt',  &
                & status='unknown', action = 'write')
            do ie = 1,nmx
                write(25,*) ymx(ie)
                if (isnan(ymx(ie))) then 
                    print*,'NAN is here...',ie,mod(ie,300)
                endif
            enddo
            close(25)
        endif


        if (any(isnan(amx))) then 
            open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx.txt',  &
                & status='unknown', action = 'write')
            do ie = 1,nmx
                write(25,*) (amx(ie,ie2),ie2=1,nmx)
                do ie2 = 1,nmx
                    if (isnan(amx(ie,ie2))) then 
                        print*,'NAN is here...',ie,mod(ie,nz),ie2,mod(ie2,nz)
                    endif
                enddo
            enddo
            close(25)
        endif

        stop
    endif

    call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,INFO) 

    if (any(isnan(ymx))) then
        print*,'error in soultion'
        open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_aft.txt',  &
        & status='unknown', action = 'write')
        do ie = 1,nmx
            write(25,*) ymx(ie)
            if (isnan(ymx(ie))) then 
                print*,'NAN is here...',ie,mod(ie,300)
            endif
        enddo
        close(25)
        stop
    endif

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (isnan(ymx(row))) then 
            print *,'nan at', iz,z(iz),zrxn,'pyrite'
            if (z(iz)<zrxn) then 
                ymx(row)=0d0
            endif
        endif

        if (ymx(row) >10d0) then 
            msx(iz) = msx(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            msx(iz) = msx(iz)*0.50d0
        else   
            msx(iz) = msx(iz)*exp(ymx(row))
        endif

    end do

    do iz = 1, nz
        row = 2 + nsp*(iz-1)
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

        if (ymx(row) >10d0) then 
            cx(iz) = cx(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            cx(iz) = cx(iz)*0.50d0
        else
            if (.not.abs(cx(iz)*exp(ymx(row)))>infinity) then 
                cx(iz) = cx(iz)*exp(ymx(row))
            endif 
        endif

        if (ymx(row+1) >10d0) then 
            po2x(iz) = po2x(iz)*1.5d0
        else if (ymx(row+1) < -10d0) then 
            po2x(iz) = po2x(iz)*0.50d0
        else
            po2x(iz) = po2x(iz)*exp(ymx(row+1))
        endif
        if (ymx(row+2) >10d0) then 
            c2x(iz) = c2x(iz)*1.5d0
        else if (ymx(row+2) < -10d0) then 
            c2x(iz) = c2x(iz)*0.50d0
        else
            c2x(iz) = c2x(iz)*exp(ymx(row+2))
        endif

        if (swoxa==0d0) then  ! ignoring Fe2 and Fe3
            ymx(row+2) = 0d0
        endif
    end do 

    error = maxval(exp(abs(ymx))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(cx)) .or. any(isnan(c2x)) &
        & .or. any(isnan(po2x)) .or. any(isnan(msx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(cx)), any(isnan(c2x)), any(isnan(po2x)),any(isnan(msx)))'
        print*,isnan(error),info/=0,any(isnan(cx)), any(isnan(c2x)), any(isnan(po2x)),any(isnan(msx))
        stop
        cx = c
        po2x =po2
        msx = ms
        c2x = c2
        iter = iter + 1
        cycle
    endif

    if (any(abs(cx)>infinity).or.any(abs(c2x)>infinity) .or.any(abs(po2x)>infinity) &
        & .or.any(abs(msx)>infinity) ) then 
        error = 1d3 
        print*, '**** ERROR; values are ininity'
        print*,any(abs(cx)>infinity),any(abs(c2x)>infinity) ,any(abs(po2x)>infinity),any(abs(msx)>infinity)
        do iz=1,nz
            if (abs(cx(iz))>infinity) cx(iz)=c(iz)
            if (abs(c2x(iz))>infinity) c2x(iz)=c2(iz)
            if (abs(po2x(iz))>infinity) po2x(iz)=po2(iz)
            if (abs(msx(iz))>infinity) msx(iz)=ms(iz)
        enddo 
        iter = iter + 1
        cycle
    endif 

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (msx(iz) < 0.0d0) then
            msx(iz) = msx(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if
    end do

    do iz = 1, nz
        row = 2 + nsp*(iz-1)

        if (po2x(iz) < 0.0d0) then
            po2x(iz) = po2x(iz)/exp(ymx(row+1))*0.5d0
            error = 1.0d0
        end if
        if (cx(iz) < 0.0d0) then
            cx(iz) = cx(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if 

    end do 

#ifdef display      
    print *, 'py error',error,info
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
    end if

end do

endsubroutine pyweath_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine pyweath_1D_SO4( &
    & nz,c,c2,po2,ms,hr,po2th,poro,z,dz,koxs2,koxs,dso4,sat,dporodta  &! input
    & ,cth,tora,v,tol,zrxn,dt,cx,c2x,po2x,msx,so4,swoxa,O2_evolution,so4i,so4th &! input
    & ,so4x &! output
    & )
    
implicit none 

integer,intent(in)::nz
logical,intent(in)::O2_evolution
real(kind=8),intent(in)::po2th,dz,dso4,cth,tol,zrxn,dt,swoxa,so4i,so4th
real(kind=8),dimension(nz),intent(in)::c,c2,po2,ms,hr,poro,z,koxs2,koxs,sat,dporodta,tora,v,cx,c2x,po2x,msx,so4
real(kind=8),dimension(nz),intent(out)::so4x

integer iz,row,ie,ie2

real(kind=8)::swex = 0.0d0 ! switch for explicit
real(kind=8)::frex = 0.0d0 ! fraction of explicit

real(kind=8) amx2(nz,nz),ymx2(nz)
integer ipiv2(nz)
integer info

external DGESV


amx2=0.0d0
ymx2=0.0d0

do iz = 1, nz

    row = iz

    if (.not.((iz == 1).or.(iz==nz))) then

        amx2(row,row) = ( &
            & 1.0d0/dt &
            & +dporodta(iz)   &
            & +(1d0-swex)*(-dso4*tora(iz)*(-2d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
            & *(1d0)/(dz**2d0)) &
            & + v(iz)/dz*(1.0d0-swex) &
            & )

        amx2(row,row-1) = ( &
            & +(1d0-swex)*(-dso4*tora(iz)*(1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz) &
            & *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
            & *(-1d0)/(dz**2d0)) &
            & - (1.0d0-swex)*v(iz)/dz &
            & )

        amx2(row,row+1) = ( &
            & +(1d0-swex)*(-dso4*tora(iz)*(1d0)/(dz**2d0)) &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt & 
            & +swex*(-dso4*tora(iz)*(so4(iz+1)+so4(iz-1)-2d0*so4(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
            & *(so4(iz)-so4(iz-1))/(dz**2d0)) &
            & + swex*v(iz)*(so4(iz)-so4(iz-1))/dz &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0 &
            & ,(2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )

    else if (iz == 1) then

        amx2(row,row) = ( &
            & 1.0d0/dt  &
            & +dporodta(iz)  &
            & +(1d0-swex)*(-dso4*tora(iz)*(-2d0)/(dz**2d0)) &
            & + v(iz)/dz*(1.0d0-swex) &
            & )

        amx2(row,row+1) = ( &
            & +(1d0-swex)*(-dso4*tora(iz)*(1d0)/(dz**2d0)) &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt  &
            & +swex*(-dso4*tora(iz)*(so4(iz+1)+so4i-2d0*so4(iz))/(dz**2d0)) &
            & + v(iz)*(so4(iz)-so4i)/dz*swex &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0 &
            & ,(2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )

    else if (iz == nz) then

        amx2(row,row) = ( &
            & 1.0d0/dt  &
            & +dporodta(iz)  &
            & +(1d0-swex)*(-dso4*tora(iz)*(-1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
            & + v(iz)/dz*(1.0d0-swex) &
            & )

        amx2(row,row-1) = ( &
            & +(1d0-swex)*(-dso4*tora(iz)*(1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
            & - (1.0d0-swex)*v(iz)/dz &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt  &
            & +swex*(-dso4*tora(iz)*(so4(iz-1)-1d0*so4(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(so4(iz)-so4(iz-1))/(dz**2d0)) &
            & + swex*v(iz)*(so4(iz)-so4(iz-1))/dz &
            & - 2d0*(1.0d0-frex)*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - 2d0*frex*koxs(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz) &
            & *merge(0d0,po2(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2(iz)**0.50d0))*1d-3 &
            & -(1.0d0-frex)*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & -frex*merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*23.94d0*1d-6*ms(iz)*c2(iz)**0.93d0*c(iz)**(-0.40d0), &
            & c(iz)<cth.or.isnan(c2(iz)**0.93d0*c(iz)**(-0.40d0)))*1d-3 &
            & )

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

call DGESV((Nz),int(1),amx2,(Nz),IPIV2,ymx2,(Nz),INFO) 

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
        so4x(iz) = ymx2(row)
    else
        if ((.not.isnan(ymx2(row))).and.(abs(ymx2(row))<=tol)) then
            so4x(iz) = 0d0
        else 
            print *,'so4 nan? or negative; stop'
            print*,iz,ymx2(row)
            stop
        endif 
        so4x(iz) = 0d0
    endif

    if (O2_evolution .and. swoxa==0d0) then  ! ignoring so4
        ymx2(row)=0d0
    endif
end do 

endsubroutine pyweath_1D_SO4

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine abweath_1D( &
    & nz,na,msil,hr,poro,z,dz,w,ksil,keqsil,msilth,dna,sat,dporodta,pro,msili,msilsupp  &! input
    & ,kco2,k1,k2,dt2,nath,tora,v,tol,nsp3,zrxn,it,cx,c2x,so4x,ca,pco2i,nai,mgx &! input
    & ,iter,error,dt,flgback &! inout
    & ,nax,prox,co2,hco3,co3,naeq,silsat,dic,msilx &! output
    & )
    
implicit none 

integer,intent(in)::nz,nsp3
real(kind=8),intent(in)::dz,w,msilth,dt2,tol,zrxn,dna,nath,pco2i,kco2,k1,k2,msili,nai,keqsil
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,na,msil,ksil,cx,c2x,so4x,ca,pro,msilsupp,mgx
real(kind=8),dimension(nz),intent(out)::nax,prox,co2,hco3,co3,naeq,silsat,dic,msilx
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer iz,row,nmx,ie,ie2

real(kind=8)::swex = 0.0d0 ! switch for explicit
real(kind=8)::frex = 0.0d0 ! fraction of explicit
real(kind=8)::swpe = 0.0d0 ! physical erosion
real(kind=8)::swad = 1.0d0 ! 1.0 when advection included 0.0d0 when not
real(kind=8),dimension(nz)::dprodna

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV



error = 1d4
iter = 0

if (it ==0) then 
    prox(:) = 0.5d0* (  &
        & -1d0*(2d0*cx(:)+2d0*ca(:)+2d0*mgx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((2d0*cx(:)+2d0*ca(:)+2d0*mgx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )
else 

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

endif

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    dprodna(:) = 0.5d0* ( &
        &  -1d0  &
        & + 0.5d0/sqrt( &
        & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i)*2d0 &
        & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) &
        & )

    do iz = 1, nz  !================================

        row = nsp3*(iz-1)+1

        if (iz==nz) then

            amx3(row,row) = (1.0d0/dt  &
                & + w/dz*(1.0d0-swex) &
                & + (1.0d0-frex)*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & ) &
                & *merge(1.0d0,msilx(iz),msilx(iz)<msilth)

            ymx3(row) = ( &
                & (msilx(iz)-msil(iz))/dt &
                & -w*(msili-msilx(iz))/dz*(1.0d0-swex) &
                & -w*(msili-msil(iz))/dz*swex*dt2/dt*swpe &
                & + (1.0d0-frex)*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & + frex*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0) &
                & -msilsupp(iz)  &
                & ) &
                & *merge(0.0d0,1d0,msilx(iz)<msilth)

        else

            amx3(row,row) = (1.0d0/dt     &
                & + w/dz *(1.0d0-swex)    &
                & + (1.0d0-frex)*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & ) &
                & * merge(1.0d0,msilx(iz),msilx(iz)<msilth)

            amx3(row,row+nsp3) = (-w/dz)*(1.0d0-swex) *merge(1.0d0,msilx(iz+1),msilx(iz)<msilth)

            ymx3(row) = ( &
                & (msilx(iz)-msil(iz))/dt &
                & -w*(msilx(iz+1)-msilx(iz))/dz*(1.0d0-swex)  &
                & -w*(msil(iz+1)-msil(iz))/dz*swex/dt*dt2*swpe &
                & + (1.0d0-frex)*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & + frex*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0) &
                & -msilsupp(iz)  &
                & ) &
                & *merge(0.0d0,1d0,msilx(iz)<msilth)

        end if 

        amx3(row,row + 1 ) = ( &
            & + (1.0d0-frex)* &
            & ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
            & + (1.0d0-frex)*ksil(iz)*poro(iz)*hr(iz)*100.07d0*1d-6*msilx(iz) &
            & *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
            & ) &
            & *nax(iz) &
            & *merge(0.0d0,1d0,msilx(iz)<msilth)

    end do  !================================

    do iz = 1, nz

        row = nsp3*(iz-1)+2

        if (.not.((iz == 1).or.(iz==nz))) then

            amx3(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dna*tora(iz)*(-2d0)/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & + (1.0d0-frex)*( &
                & -(1.0d0)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & )*1d-3 &
                & + (1.0d0-frex)*( &
                & -(1.0d0)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz) &
                & *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0) &
                & )*1d-3 &
                & ) &
                & *merge(1.0d0,nax(iz),nax(iz)<nath)

            amx3(row,row-nsp3) = ( &
                & +(1d0-swex)*(-dna*tora(iz)*(1d0)/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - (1.0d0-swex)*v(iz)/dz &
                & ) &
                & *nax(iz-1) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

            amx3(row,row+nsp3) = ( &
                & +(1d0-swex)*(-dna*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *nax(iz+1) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

            ymx3(row) = ( &
                & (nax(iz)-na(iz))/dt  &
                & +dporodta(iz) *nax(iz) &
                & +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nax(iz-1)-2d0*nax(iz))/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(nax(iz)-nax(iz-1))/(dz**2d0)) &
                & +swex*(-dna*tora(iz)*(na(iz+1)+na(iz-1)-2d0*na(iz))/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(na(iz)-na(iz-1))/(dz**2d0)) &
                & + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz &
                & + swex*v(iz)*(na(iz)-na(iz-1))/dz &
                & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

        else if (iz == 1) then

            amx3(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dna*tora(iz)*(-2d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & -(1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & -(1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)&
                & *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & ) &
                & *merge(1.0d0,nax(iz),nax(iz)<nath)

            amx3(row,row+nsp3) = ( &
                & +(1d0-swex)*(-dna*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *nax(iz+1) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

            ymx3(row) = ( &
                & (nax(iz)-na(iz))/dt  &
                & +dporodta(iz) *nax(iz) &
                & +(1d0-swex)*(-dna*tora(iz)*(nax(iz+1)+nai-2d0*nax(iz))/(dz**2d0)) &
                & +swex*(-dna*tora(iz)*(na(iz+1)+nai-2d0*na(iz))/(dz**2d0)) &
                & + v(iz)*(nax(iz)-nai)/dz*(1.0d0-swex) &
                & + v(iz)*(na(iz)-nai)/dz*swex &
                & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

        else if (iz == nz) then

            amx3(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & +(1d0-swex)*(-dna*tora(iz)*(-1d0)/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                & + v(iz)/dz*(1.0d0-swex) &
                & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(-4d0*3d0*nax(iz)**2d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz) &
                & *(-4d0*nax(iz)**3d0*(-1d0)/(prox(iz)**2d0)/keqsil)*dprodna(iz) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & ) &
                & *merge(1.0d0,nax(iz),nax(iz)<nath)

            amx3(row,row-nsp3) = ( &
                & +(1d0-swex)*(-dna*tora(iz)*(1d0)/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - (1.0d0-swex)*v(iz)/dz &
                & ) &
                & *nax(iz-1) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

            ymx3(row) = ( &
                & (nax(iz)-na(iz))/dt  &
                & +dporodta(iz) *nax(iz) &
                & +(1d0-swex)*(-dna*tora(iz)*(nax(iz-1)-1d0*nax(iz))/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(nax(iz)-nax(iz-1))/(dz**2d0)) &
                & +swex*(-dna*tora(iz)*(na(iz-1)-1d0*na(iz))/(dz**2d0) &
                & -dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(na(iz)-na(iz-1))/(dz**2d0)) &
                & + (1.0d0-swex)*v(iz)*(nax(iz)-nax(iz-1))/dz &
                & + swex*v(iz)*(na(iz)-na(iz-1))/dz &
                & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msilx(iz)*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3 &
                & - frex*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*msil(iz)*(1d0-4d0*na(iz)**3d0/pro(iz)/keqsil) &
                & *merge(0d0,1d0,1d0-4d0*na(iz)**3d0/pro(iz)/keqsil < 0d0)*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)

        end if 

        amx3(row,row  - 1) = (     & 
            & - (1.0d0-frex)*ksil(iz)/sat(iz)*hr(iz)*100.07d0*1d-6*1d0*(1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil) &
            & *merge(0d0,1d0,1d0-4d0*nax(iz)**3d0/prox(iz)/keqsil < 0d0)*1d-3  &
            & ) &
            & *msilx(iz) &
            & *merge(0.0d0,1.0d0,nax(iz)<nath)   ! commented out (is this necessary?)
        
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

    call DGESV(nsp3*(Nz),int(1),amx3,nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 

    if (any(isnan(ymx3))) then
        print*,'error in soultion'
    endif

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'albite'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
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
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(nax)) .or. any(isnan(msilx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(nax)),any(isnan(msilx)))'
        print*,isnan(error),info/=0,any(isnan(nax)),any(isnan(msilx))
        stop
        nax = na
        msilx = msil
        prox = pro
        iter = iter + 1
        cycle
    endif

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    naeq = (keqsil *prox/4d0)**(1d0/3d0)
    silsat = nax/naeq

    co2 = kco2*pco2i
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (msilx(iz) < 0.0d0) then
            msilx(iz) = msilx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
    end do

    do iz = 1, nz
        row = 2 + nsp3*(iz-1)

        if (nax(iz) < 0.0d0) then
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'ab error',error,info
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
    end if


enddo

endsubroutine abweath_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine basaltweath_1D( &
    & nz,mfo,mg,si,hr,poro,z,dz,w,kfo,keqfo,mfoth,dmg,dsi,sat,dporodta,pro,mfoi,mfosupp  &! input
    & ,kco2,k1,k2,mgth,sith,tora,v,tol,zrxn,it,nax,cx,c2x,so4x,ca,pco2i,mgi,sii,mvfo,nflx &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,prox,co2,hco3,co3,dic,mfox,omega_fo,flx_fo,flx_mg,flx_si &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::dz,w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,pco2i,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo
real(kind=8),dimension(nz),intent(in)::nax,hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si
real(kind=8),dimension(nz),intent(out)::mgx,six,prox,co2,hco3,co3,dic,mfox,omega_fo
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 3
integer iz,row,nmx,ie,ie2,isp
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,domega_fo_disp,caq_tmp_prev,st_fo

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! print *, mgx
! print *, six
! print *, mfosupp
! stop

error = 1d4
iter = 0

if (it ==0) then 
    prox(:) = 0.5d0* (  &
        & -1d0*(2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )
else 

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

endif

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_fo = 0d0
    flx_mg = 0d0
    flx_si = 0d0

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    dprodna(:) = 0.5d0* ( &
        &  -1d0  &
        & + 0.5d0/sqrt( &
        & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i)*2d0 &
        & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) &
        & )

    dprodmg(:) = 0.5d0* ( &
        &  -2d0  &
        & + 0.5d0/sqrt( &
        & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i)*2d0 &
        & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))*2d0 &
        & )
    
    omega_fo(:) = mgx(:)**2d0*six(:)/(prox(:)**4d0)/keqfo
    domega_fo_dmg(:) = 2d0*mgx(:)*six(:)/(prox(:)**4d0)/keqfo+ 2d0*mgx(:)*six(:)*(-4d0)/(prox(:)**5d0)*dprodmg(:)/keqfo
    domega_fo_dsi(:) = mgx(:)**2d0/(prox(:)**4d0)/keqfo
    
    ! omega_fo(:) = mg(:)**2d0*si(:)/(pro(:)**4d0)/keqfo
    ! domega_fo_dmg(:) = 0d0
    ! domega_fo_dsi(:) = 0d0
    
    ! print *,omega_fo
    ! print *,domega_fo_dmg
    ! print *,domega_fo_dsi
    ! stop

    do iz = 1, nz  !================================

        row = nsp3*(iz-1)+1

        if (iz==nz) then

            amx3(row,row) = (1.0d0/dt  &
                & + w/dz  &
                & + kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                & ) &
                & *merge(1.0d0,mfox(iz),mfox(iz)<mfoth)

            ymx3(row) = ( &
                & (mfox(iz)-mfo(iz))/dt &
                & -w*(mfoi-mfox(iz))/dz &
                & + kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                & -mfosupp(iz)  &
                & ) &
                & *merge(0.0d0,1d0,mfox(iz)<mfoth)
            
            flx_fo(iadv,iz) = (&
                & -w*(mfoi-mfox(iz))/dz &
                & )

        else

            amx3(row,row) = (1.0d0/dt     &
                & + w/dz    &
                & + kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                & ) &
                & * merge(1.0d0,mfox(iz),mfox(iz)<mfoth)

            amx3(row,row+nsp3) = (-w/dz) *merge(1.0d0,mfox(iz+1),mfox(iz)<mfoth)

            ymx3(row) = ( &
                & (mfox(iz)-mfo(iz))/dt &
                & -w*(mfox(iz+1)-mfox(iz))/dz  &
                & + kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                & -mfosupp(iz)  &
                & ) &
                & *merge(0.0d0,1d0,mfox(iz)<mfoth)
            
            flx_fo(iadv,iz) = (&
                & -w*(mfox(iz+1)-mfox(iz))/dz  &
                & )

        end if 

        amx3(row,row + 1 ) = ( &
            & kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dmg(iz)) &
            & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
            & ) &
            & *mgx(iz) &
            & *merge(0.0d0,1d0,mfox(iz)<mfoth)

        amx3(row,row + 2 ) = ( &
            & kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz)) &
            & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
            & ) &
            & *six(iz) &
            & *merge(0.0d0,1d0,mfox(iz)<mfoth)
            
        flx_fo(itflx,iz) = (&
            & (mfox(iz)-mfo(iz))/dt &
            & )
            
        flx_fo(irxn_fo,iz) = (&
            & + kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
            & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
            & )
            
        flx_fo(irain,iz) = (&
            & -mfosupp(iz)  &
            & )
        flx_fo(ires,iz) = sum(flx_fo(:,iz))

    end do  !================================

    do iz = 1, nz
        
        do isp = 1, 2

            row = nsp3*(iz-1)+1 + isp
            
            if (isp==1) then ! mg 
                d_tmp = dmg
                caq_tmp = mgx(iz)
                caq_tmp_prev = mg(iz)
                caq_tmp_p = mgx(min(nz,iz+1))
                caq_tmp_n = mgx(max(1,iz-1))
                caqth_tmp = mgth
                caqi_tmp = mgi
                domega_fo_disp = domega_fo_dmg(iz)
                st_fo = 2d0
            elseif (isp==2) then  ! si
                d_tmp = dsi
                caq_tmp = six(iz)
                caq_tmp_prev = si(iz)
                caq_tmp_p = six(min(nz,iz+1))
                caq_tmp_n = six(max(1,iz-1))
                caqth_tmp = sith
                caqi_tmp = sii
                domega_fo_disp = domega_fo_dsi(iz)
                st_fo = 1d0
            endif 

            if (.not.((iz == 1).or.(iz==nz))) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-2d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                    & + v(iz)/dz  &
                    & -st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_disp) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                    & - v(iz)/dz &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                amx3(row,row+nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0)) &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                    & -st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == 1) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-2d0)/(dz**2d0)) &
                    & + v(iz)/dz &
                    & -st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_disp) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row+nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0)) &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == nz) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                    & + v(iz)/dz &
                    & -st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_disp) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                    & - v(iz)/dz &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 


            end if 

            amx3(row,row  - isp) = (     & 
                & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*1d0*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                & ) &
                & *mfox(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
            if (isp==1) then 
                flx_mg(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_mg(irxn_fo,iz) = (&
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) 
                flx_mg(ires,iz) = sum(flx_mg(:,iz))
            elseif (isp==2) then 
                flx_si(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_si(irxn_fo,iz) = (&
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & ) 
                flx_si(ires,iz) = sum(flx_si(:,iz))
            endif 
        
        enddo 
        
    end do  ! ==============================
    
    ymx3=-1.0d0*ymx3

    if (any(isnan(amx3)).or.any(isnan(ymx3))) then 
    ! if (.true.) then 
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
        
        open(unit=11,file='amx.txt',status = 'replace')
        open(unit=12,file='ymx.txt',status = 'replace')
        do ie = 1,nsp3*(nz)
            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
            write(12,*) ymx3(ie)
        enddo 
        close(11)
        close(12)       

        stop
    endif

    call DGESV(nsp3*(Nz),int(1),amx3,nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 

    if (any(isnan(ymx3))) then
        print*,'error in soultion'
    endif

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'albite'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >10d0) then 
            mfox(iz) = mfox(iz)*1.5d0
        else if (ymx3(row) < -10d0) then 
            mfox(iz) = mfox(iz)*0.50d0
        else   
            mfox(iz) = mfox(iz)*exp(ymx3(row))
        endif

    end do

    do iz = 1, nz
        row = 2 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Na'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                mgx(iz) = 0.1d0*mgth
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >1d0) then 
            mgx(iz) = mgx(iz)*1.5d0
        else if (ymx3(row) < -1d0) then 
            mgx(iz) = mgx(iz)*0.50d0
        else
            mgx(iz) = mgx(iz)*exp(ymx3(row))
        endif
        
        
        row = 3 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Na'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                six(iz) = 0.1d0*sith
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >1d0) then 
            six(iz) = six(iz)*1.5d0
        else if (ymx3(row) < -1d0) then 
            six(iz) = six(iz)*0.50d0
        else
            six(iz) = six(iz)*exp(ymx3(row))
        endif
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)) .or. any(isnan(mfox))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox))
        stop
        mgx = mg
        six = si
        mfox = mfo
        prox = pro
        iter = iter + 1
        cycle
    endif

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    co2 = kco2*pco2i
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (mfox(iz) < 0.0d0) then
            mfox(iz) = mfox(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
    end do

    do iz = 1, nz
        row = 2 + nsp3*(iz-1)

        if (mgx(iz) < 0.0d0) then
            mgx(iz) = mgx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 3 + nsp3*(iz-1)

        if (six(iz) < 0.0d0) then
            six(iz) = six(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'basalt error',error,info
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
    end if
    
enddo

endsubroutine basaltweath_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine silicate_dis_1D( &
    & nz,mfo,mab,na,mg,si,hr,poro,z,dz,w,kfo,kab,keqfo,keqab,mfoth,mabth,dmg,dsi,dna,sat,dporodta,pro,mfoi,mabi,mfosupp,mabsupp  &! input
    & ,kco2,k1,k2,mgth,sith,nath,tora,v,tol,zrxn,it,cx,c2x,so4x,ca,pco2i,mgi,sii,nai,mvfo,mvab,nflx &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,prox,co2,hco3,co3,dic,mfox,mabx,omega_fo,omega_ab,flx_fo,flx_mg,flx_si,flx_ab,flx_na &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::dz,w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,pco2i,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp
real(kind=8),dimension(nz),intent(out)::mgx,six,prox,co2,hco3,co3,dic,mfox,omega_fo,nax,mabx,omega_ab
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 5
integer iz,row,nmx,ie,ie2,isp
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! print *, mgx
! print *, six
! print *, mfosupp
! stop

error = 1d4
iter = 0

if (it ==0) then 
    prox(:) = 0.5d0* (  &
        & -1d0*(2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )
else 

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

endif

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_fo = 0d0
    flx_mg = 0d0
    flx_si = 0d0
    flx_ab = 0d0
    flx_na = 0d0

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    dprodna(:) = 0.5d0* ( &
        &  -1d0  &
        & + 0.5d0/sqrt( &
        & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i)*2d0 &
        & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) &
        & )

    dprodmg(:) = 0.5d0* ( &
        &  -2d0  &
        & + 0.5d0/sqrt( &
        & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i)*2d0 &
        & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))*2d0 &
        & )
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    omega_fo(:) = mgx(:)**2d0*six(:)/(prox(:)**4d0)/keqfo
    domega_fo_dmg(:) = 2d0*mgx(:)*six(:)/(prox(:)**4d0)/keqfo+ 2d0*mgx(:)*six(:)*(-4d0)/(prox(:)**5d0)*dprodmg(:)/keqfo
    domega_fo_dsi(:) = mgx(:)**2d0/(prox(:)**4d0)/keqfo
    
    ! omega_fo(:) = mg(:)**2d0*si(:)/(pro(:)**4d0)/keqfo
    ! domega_fo_dmg(:) = 0d0
    ! domega_fo_dsi(:) = 0d0
    
    ! print *,omega_fo
    ! print *,domega_fo_dmg
    ! print *,domega_fo_dsi
    ! stop
    
    ! Ab + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)
    omega_ab(:) = nax(:)*six(:)**2d0/prox(:)/keqab
    domega_ab_dna(:) = six(:)**2d0/prox(:)/keqab + nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodna(:)
    domega_ab_dsi(:) = nax(:)*(2d0)*six(:)/prox(:)/keqab

    do iz = 1, nz  !================================
        
        do isp = 1, 2
        
            row = nsp3*(iz-1)+isp
            
            if (isp==1) then 
                k_tmp = kfo(iz)
                mv_tmp = mvfo
                omega_tmp = omega_fo(iz)
                m_tmp = mfox(iz)
                mth_tmp = mfoth 
                mi_tmp = mfoi
                mp_tmp = mfox(min(nz,iz+1))
                msupp_tmp = mfosupp(iz)
                mprev_tmp = mfo(iz)
            elseif (isp==2)then 
                k_tmp = kab(iz)
                mv_tmp = mvab
                omega_tmp = omega_ab(iz)
                m_tmp = mabx(iz)
                mth_tmp = mabth 
                mi_tmp = mabi
                mp_tmp = mabx(min(nz,iz+1))
                msupp_tmp = mabsupp(iz)
                mprev_tmp = mab(iz)
            endif 

            if (iz==nz) then

                amx3(row,row) = (1.0d0/dt  &
                    & + w/dz  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & *merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp)/dt &
                    & -w*(mi_tmp-m_tmp)/dz &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & -msupp_tmp  &
                    & ) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                
                if (isp==1) then 
                    flx_fo(iadv,iz) = (&
                        & -w*(mi_tmp-m_tmp)/dz &
                        & )
                elseif (isp==2) then 
                    flx_ab(iadv,iz) = (&
                        & -w*(mi_tmp-m_tmp)/dz &
                        & )
                endif 

            else

                amx3(row,row) = (1.0d0/dt     &
                    & + w/dz    &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & * merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                amx3(row,row+nsp3) = (-w/dz) *merge(1.0d0,mp_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp)/dt &
                    & -w*(mp_tmp-m_tmp)/dz  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & -msupp_tmp  &
                    & ) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                
                if (isp==1) then 
                    flx_fo(iadv,iz) = (&
                        & -w*(mp_tmp-m_tmp)/dz  &
                        & )
                elseif (isp==2) then 
                    flx_ab(iadv,iz) = (&
                        & -w*(mp_tmp-m_tmp)/dz  &
                        & )
                endif 

            end if 
            
            if (isp==1) then 
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_fo(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_fo(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & )
                    
                flx_fo(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_fo(ires,iz) = sum(flx_fo(:,iz))
                
            elseif (isp==2) then 
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_ab(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_ab(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp < 0d0) &
                    & )
                    
                flx_ab(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_ab(ires,iz) = sum(flx_ab(:,iz))
            endif 
        enddo 
    end do  !================================

    do iz = 1, nz
        
        do isp = 1, 3

            row = nsp3*(iz-1)+2 + isp
            
            if (isp==1) then ! mg 
                d_tmp = dmg
                caq_tmp = mgx(iz)
                caq_tmp_prev = mg(iz)
                caq_tmp_p = mgx(min(nz,iz+1))
                caq_tmp_n = mgx(max(1,iz-1))
                caqth_tmp = mgth
                caqi_tmp = mgi
                st_fo = 2d0
                st_ab = 0d0
                rxn_tmp = st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 
                drxndisp_tmp = st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dmg(iz)) & 
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3
            elseif (isp==2) then  ! si
                d_tmp = dsi
                caq_tmp = six(iz)
                caq_tmp_prev = si(iz)
                caq_tmp_p = six(min(nz,iz+1))
                caq_tmp_n = six(max(1,iz-1))
                caqth_tmp = sith
                caqi_tmp = sii
                st_fo = 1d0
                st_ab = 2d0
                rxn_tmp = st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & + st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
                drxndisp_tmp = st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3 &
                    & + st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
            elseif (isp==3) then  ! na
                d_tmp = dna
                caq_tmp = nax(iz)
                caq_tmp_prev = na(iz)
                caq_tmp_p = nax(min(nz,iz+1))
                caq_tmp_n = nax(max(1,iz-1))
                caqth_tmp = nath
                caqi_tmp = nai
                st_fo = 0d0
                st_ab = 1d0
                rxn_tmp = st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
                drxndisp_tmp = st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
            endif 

            if (.not.((iz == 1).or.(iz==nz))) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-2d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                    & + v(iz)/dz  &
                    & -drxndisp_tmp &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                    & - v(iz)/dz &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                amx3(row,row+nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0)) &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                    & - rxn_tmp &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == 1) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-2d0)/(dz**2d0)) &
                    & + v(iz)/dz &
                    & - drxndisp_tmp &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row+nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0)) &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                    & - rxn_tmp &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == nz) then

                amx3(row,row) = ( &
                    & 1.0d0/dt  &
                    & +dporodta(iz)  &
                    & +(-d_tmp*tora(iz)*(-1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0)) &
                    & + v(iz)/dz &
                    & - drxndisp_tmp &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                    & - v(iz)/dz &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                    & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                    & - rxn_tmp &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 


            end if 
            
            amx3(row,row  - isp - 1) = (     & 
                & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*1d0*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                & ) &
                & *mfox(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp ) = (     & 
                & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*1d0*(1d0-omega_ab(iz)) &
                & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                & ) &
                & *mabx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            if (isp==1) then 
                flx_mg(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_mg(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_mg(ires,iz) = sum(flx_mg(:,iz))
            elseif (isp==2) then 
                flx_si(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_si(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_si(ires,iz) = sum(flx_si(:,iz))
            elseif (isp==3) then 
                flx_na(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_na(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_na(ires,iz) = sum(flx_na(:,iz))
            endif 
        
        enddo 
        
    end do  ! ==============================
    
    ymx3=-1.0d0*ymx3

    if (any(isnan(amx3)).or.any(isnan(ymx3))) then 
    ! if (.true.) then 
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
        
        open(unit=11,file='amx.txt',status = 'replace')
        open(unit=12,file='ymx.txt',status = 'replace')
        do ie = 1,nsp3*(nz)
            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
            write(12,*) ymx3(ie)
        enddo 
        close(11)
        close(12)       

        stop
    endif

    call DGESV(nsp3*(Nz),int(1),amx3,nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 

    if (any(isnan(ymx3))) then
        print*,'error in soultion'
    endif

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Fo'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >10d0) then 
            mfox(iz) = mfox(iz)*1.5d0
        else if (ymx3(row) < -10d0) then 
            mfox(iz) = mfox(iz)*0.50d0
        else   
            mfox(iz) = mfox(iz)*exp(ymx3(row))
        endif
        
        row = 2 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Ab'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >10d0) then 
            mabx(iz) = mabx(iz)*1.5d0
        else if (ymx3(row) < -10d0) then 
            mabx(iz) = mabx(iz)*0.50d0
        else   
            mabx(iz) = mabx(iz)*exp(ymx3(row))
        endif
        
        row = 3 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'mg'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                mgx(iz) = 0.1d0*mgth
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >1d0) then 
            mgx(iz) = mgx(iz)*1.5d0
        else if (ymx3(row) < -1d0) then 
            mgx(iz) = mgx(iz)*0.50d0
        else
            mgx(iz) = mgx(iz)*exp(ymx3(row))
        endif
        
        row = 4 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'si'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                six(iz) = 0.1d0*sith
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >1d0) then 
            six(iz) = six(iz)*1.5d0
        else if (ymx3(row) < -1d0) then 
            six(iz) = six(iz)*0.50d0
        else
            six(iz) = six(iz)*exp(ymx3(row))
        endif
        
        row = 5 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'na'
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
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx))
        stop
        mgx = mg
        six = si
        nax = na
        mfox = mfo
        mabx = mab
        prox = pro
        iter = iter + 1
        cycle
    endif

    prox(:) = 0.5d0* ( &
        & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        & + 4d0*kco2*k1*pco2i) &
        & )

    co2 = kco2*pco2i
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (mfox(iz) < 0.0d0) then
            mfox(iz) = mfox(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 2 + nsp3*(iz-1)

        if (mabx(iz) < 0.0d0) then
            mabx(iz) = mabx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 3 + nsp3*(iz-1)

        if (mgx(iz) < 0.0d0) then
            mgx(iz) = mgx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 4 + nsp3*(iz-1)

        if (six(iz) < 0.0d0) then
            six(iz) = six(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 5 + nsp3*(iz-1)

        if (nax(iz) < 0.0d0) then
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'silicate_dis error',error,info
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
    end if
    
enddo

endsubroutine silicate_dis_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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