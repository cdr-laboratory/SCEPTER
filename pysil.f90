program weathering
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

real(kind=8) :: ztot = 10.0d0 ! m
! real(kind=8) dz
integer, parameter :: nz = 100 
real(kind=8) z(nz),dz(nz)
real(kind=8) ze(nz+1)
real(kind=8) :: ph = 5.0d0
real(kind=8) :: tc = 15.0d0 ! deg celsius
real(kind=8) dt  ! yr 
integer, parameter :: nt = 50000000
real(kind=8) time
integer, parameter :: nsp = 4
integer, parameter :: nsp3= 2

real(kind=8),parameter :: rg = 8.3d-3   ! kJ mol^-1 K^-1
real(kind=8) :: rg2 = 8.2d-2  ! L mol^-1 atm K^-1

real(kind=8),parameter :: tempk_0 = 273d0
real(kind=8),parameter :: sec2yr = 60d0*60d0*24d0*365d0

real(kind=8) :: po2i = 0.21d0 ! atm **default
! real(kind=8) :: po2i = 0.6d-1 ! atm
real(kind=8) :: pco2i = 10.0d0**(-3.5d0) ! atm **default 
! real(kind=8) :: pco2i = 10.0d0**(-1.0d0) ! atm
real(kind=8) :: ci = 0d0 
real(kind=8) :: c2i = 0d0
real(kind=8) :: so4i = 0d0
real(kind=8) :: nai = 0d0
real(kind=8) :: mgi = 0d0
real(kind=8) :: sii = 0d0
real(kind=8) :: cai = 0d0

real(kind=8) :: mvkaol = 99.52d0 ! cm3/mol; molar volume of kaolinite; Robie et al. 1978
real(kind=8) :: mvfo = 43.79d0 ! cm3/mol; molar volume of Fo; Robie et al. 1978
real(kind=8) :: mvab = 100.07d0 ! cm3/mol; molar volume of Ab(NaAlSi3O8); Robie et al. 1978 
real(kind=8) :: mvan = 100.79d0 ! cm3/mol; molar volume of An (CaAl2Si2O8); Robie et al. 1978
real(kind=8) :: mvcc = 36.934d0 ! cm3/mol; molar volume of Cc (CaCO3); Robie et al. 1978
real(kind=8) :: mvpy = 23.94d0 ! cm3/mol; molar volume of Pyrite (FeS2); Robie et al. 1978

real(kind=8) :: mwtfo = 140.694d0 ! g/mol; molar volume of Fo; Robie et al. 1978
real(kind=8) :: mwtab = 262.225d0 ! g/mol; molar volume of Ab; Robie et al. 1978
real(kind=8) :: mwtan = 278.311d0 ! g/mol; molar volume of An; Robie et al. 1978
real(kind=8) :: mwtcc = 100.089d0 ! g/mol; molar volume of Cc; Robie et al. 1978

! real(kind=8) :: redsldi = 0.56d0 ! wt%  **default 
! real(kind=8) :: redsldi = 1.12d0 ! wt%  x2
real(kind=8) :: redsldi = 2.8d0 ! wt%   x5
! real(kind=8) :: redsldi = 2.24d0 ! wt%  x4
! real(kind=8) :: redsldi = 3.36d0 ! wt%  x6

real(kind=8) :: silwti = 30d0 ! wt%  **default
! real(kind=8) :: silwti = 45d0 ! wt%  
! real(kind=8) :: silwti = 24d0 ! wt%
! real(kind=8) :: silwti = 1d-10 ! wt%

! real(kind=8)::rainpowder = 40d2 !  g/m2/yr corresponding to 40 t/ha/yr (40x1e3x1e3/1e4)
! real(kind=8)::rainpowder = 0.5d2 !  g/m2/yr corresponding to 0.5 t/ha/yr (0.5x1e3x1e3/1e4)
real(kind=8)::rainpowder = 10d2 !  g/m2/yr 
! real(kind=8)::rainpowder = 10d2 !  g/m2/yr corresponding to 10 t/ha/yr (0.5x1e3x1e3/1e4)

real(kind=8)::rainfrc_fo = 0.12d0 ! rain wt fraction for Fo (Beering et al 2020)
real(kind=8)::rainfrc_ab = 0.172d0 ! rain wt fraction for Ab; assuming 0.43 for La and 0.4 of wt of La is Ab (Beering et al 2020)
! real(kind=8)::rainfrc_ab = 0d0 ! rain wt fraction for Ab; assuming 0.43 for La and 0.4 of wt of La is Ab (Beering et al 2020)
real(kind=8)::rainfrc_an = 0.258d0 ! rain wt fraction for An; assuming 0.43 for La and 0.6 of wt of La is An (Beering et al 2020)
! real(kind=8)::rainfrc_an = 0d0 ! rain wt fraction for An; assuming 0.43 for La and 0.6 of wt of La is An (Beering et al 2020)
real(kind=8)::rainfrc_cc = 0d0 ! rain wt fraction for Cc; None (Beering et al 2020)

real(kind=8)::zsupp = 0.3d0 !  e-folding decrease

real(kind=8) sat(nz), poro(nz), torg(nz), tora(nz), deff(nz)
real(kind=8) :: dgas = 6.09d2 ! m^2 yr^-1
real(kind=8) :: daq = 5.49d-2 ! m^2 yr^-1
real(kind=8) :: dgasc = 441.504d0 ! m^2 yr^-1 (Assuming 0.14 cm2/sec)
real(kind=8) :: daqc = 0.022459852 ! m^2 yr^-1 (for C32- from Li and Gregory 1974)
real(kind=8) :: poroi = 0.1d0
real(kind=8) :: sati = 0.50d0
real(kind=8) :: satup = 0.10d0

! real(kind=8) :: zsat = 30d0  ! water table depth [m] ** default 
real(kind=8) :: zsat = 5d0  ! water table depth [m] 
! real(kind=8) :: zsat = 15d0

real(kind=8) :: dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
real(kind=8) :: dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dmg  = 0.017218079d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dsi  = 0.03689712d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dca  = 0.019023312d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 

real(kind=8) :: w = 5.0d-5 ! m yr^-1, uplift rate ** default 
! real(kind=8), parameter :: w = 1.0d-4 ! m yr^-1, uplift rate

real(kind=8), parameter :: vcnst = 1.0d1 ! m yr^-1, advection
! real(kind=8), parameter :: qin = 5d-3 ! m yr^-1, advection (m3 water / m2 profile / yr)

! real(kind=8) :: qin = 1d-1 ! m yr^-1, advection (m3 water / m2 profile / yr)  ** default
real(kind=8) :: qin = 10d-1 ! m yr^-1 
! real(kind=8) :: qin = 0.1d-1 ! m yr^-1 
real(kind=8) v(nz), q

! real(kind=8) :: hr = 1d5 ! m^2 m^-3, reciprocal of hydraulic radius  ** default 
! real(kind=8) :: hr = 1d4 ! m^2 m^-3, reciprocal of hydraulic radius
real(kind=8) :: hri = 1d5

real(kind=8) :: hr(nz)

real(kind=8) msi,msili,mfoi,mabi,mani,mcci,ctmp,po2tmp
real(kind=8),dimension(nz)::msil,msilx,mfo,mfox,mab,mabx,man,manx,mcc,mccx &
    & ,po2,redsld,ms,c,po2x,msx,cx,resp,c2,c2x,so4,so4x,na,nax,naeq,silsat &
    & ,pro,prox,hco3,ca,co2,co3,dic,cax,porox,dporodta,dporodtg,dporodtgc,khco2  &
    & ,mg,mgx,si,six,pco2,pco2x,poroprev,khco2x,pco2x_prev,torgprev,toraprev

real(kind=8) :: caeq = 1d-3  ! mol/L equilibrium Ca conc. 
real(kind=8) :: delca = 0.5d0  ! m reaction front width  
real(kind=8) :: zca = 50d0   ! m depth of reaction front for calcite          

real(kind=8),dimension(nz)::koxa,koxs,koxs2,ksil,msilsupp,kfo,mfosupp,omega_fo &
    & ,kab,mabsupp,omega_ab,kan,mansupp,omega_an,kcc,mccsupp,omega_cc,preccc 

real(kind=8) kho,ucv,kco2,k1,keqsil,kw,k2,keqfo,keqab,keqgb,khco2i,keqan,keqcc

integer iz, ie, it, ie2, iter_co2

real(kind=8) error, error2, error_poro,error_co2
real(kind=8) :: tol = 1d-6

! integer, parameter :: nrec = 22
integer, parameter :: nrec = 20
integer reclis(nrec)
real(kind=8) rectime(nrec)
character(3) chr
character(256) runname,workdir, chrz(3), chrq(3),base,fname, runname_save,chrrain
integer irec, iter, idum, iter2

real(kind=8) dt2, swpe  ! physical erosion 

real(kind=8) :: po2th = 1.0d-20
real(kind=8) :: pco2th = 1.0d-20

real(kind=8) :: cth = 1.0d-20
real(kind=8) :: c2th = 1.0d-20
real(kind=8) :: so4th = 1.0d-20
real(kind=8) :: proth = 1.0d-20
real(kind=8) :: nath = 1.0d-20
real(kind=8) :: mgth = 1.0d-20
real(kind=8) :: sith = 1.0d-20
real(kind=8) :: cath = 1.0d-20
real(kind=8) :: msth = 1.0d-300
real(kind=8) :: msilth = 1.0d-300
real(kind=8) :: mfoth = 1.0d-300
real(kind=8) :: mabth = 1.0d-300
real(kind=8) :: manth = 1.0d-300
real(kind=8) :: mccth = 1.0d-300

real(kind=8) :: stoxs = 15.0d0/4.0d0  ! 15/4 py => Fe-oxide + sulfate; 7/2 py => Fe++ + sulfate
real(kind=8) :: stoxa = 1.0d0/4.0d0  ! stoichiomety of oxidation in aq
real(kind=8) :: swoxa = 0.0d0   ! switch for oxidation in aq
real(kind=8) :: swoxs2 = 1.0d0  ! switch for oxidation in solid phase
real(kind=8) :: swoxall = 0d0   ! switch when only consider overall oxidation, i.e., 15/4 py => Fe-oxide + sulfate

real(kind=8) :: swbr = 0.0d0  ! switch for biological respiration
real(kind=8) :: vmax = 0.71d0 ! mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
real(kind=8) :: mo2 = 0.121d0 ! Michaelis, Davidson et al. (2012)

real(kind=8) :: swadvmass = 0d0 ! switch; 1 when calculating q from advection mass balance
real(kind=8) :: waterfluc = 0d0 ! switch: 1 when fluctuating water flow

integer,parameter :: nflx = 6, nflx_py = 8
integer  iflx
real(kind=8),dimension(nflx,nz)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_o2,flx_an,flx_ca,flx_cc,flx_co2
real(kind=8),dimension(nflx_py,nz)::flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2,flx_py_so4
! real(kind=8) :: maxdt = 10d0
real(kind=8) :: maxdt = 0.2d0 ! for basalt exp?

! logical :: pre_calc = .false.
logical :: pre_calc = .true.

logical :: read_data = .false.
! logical :: read_data = .true.

logical :: initial_ss = .false.
! logical :: initial_ss = .true.

logical :: rain_wave = .false.
! logical :: rain_wave = .true.

logical :: co2_iteration = .false.
! logical :: co2_iteration = .true.

real(kind=8) :: wave_tau = 2d0 ! yr periodic time for wave 
real(kind=8) :: rain_norm = 0d0

data rectime /1d1,3d1,1d2,3d2,1d3,3d3,1d4,3d4 &
    & ,1d5,2d5,3d5,4d5,5d5,6d5,7d5,8d5,9d5,1d6,1.1d6,1.2d6/
! data rectime /-1d6,0d6,1d6,2d6,3d6,4d6,5d6,6d6,7d6,8d6
! &,9d6,10d6,11d6,12d6,13d6,14d6,15d6,16d6,17d6,18d6,19d6,20d6/
! data rectime /21d6,22d6,23d6,24d6,25d6,26d6,27d6,28d6,29d6,30d6
! & ,31d6,32d6,33d6,34d6,35d6,36d6,37d6,38d6,39d6,40d6,41d6,42d6/
real(kind=8) :: savetime = 1d3
real(kind=8) :: dsavetime = 1d3

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

real(kind=8), parameter:: infinity = huge(0d0)
real(kind=8) k_arrhenius
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

! write(chrq(1),'(i0)') int(qin/(10d0**(floor(log10(qin)))))
! write(chrq(2),'(i0)') floor(log10(qin))
! chrq(3) = trim(adjustl(chrq(1)))//'E'//trim(adjustl(chrq(2)))
write(chrq(3),'(E10.2)') qin
write(chrz(3),'(i0)') nint(zsat)
write(chrrain,'(E10.2)') rainpowder

vmax = vmax * 1d-1  !!  vmax is increased by a factor of 100 (cf., soil respiration in Liu and Zhou 2006)

mo2 = mo2*po2i/0.21d0     !! mo2 is assumed to proportional to po2i


write(workdir,*) '../pyweath_output/'     
write(base,*) 'test_cpl_rain-'//trim(adjustl(chrrain))    
 
#ifdef poroevol 
base = trim(adjustl(base))//'_pevol'
#endif 
#if defined(surfevol1)
base = trim(adjustl(base))//'_sevol1'
#elif defined(surfevol2)
base = trim(adjustl(base))//'_sevol2'
#endif 

#ifndef regulargrid
base = trim(adjustl(base))//'_irr'
#endif 

if (rain_wave)then 
    write(chrrain,'(E10.2)') wave_tau
    base = trim(adjustl(base))//'_rwave-'//trim(adjustl(chrrain))
endif 
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
write(runname,*) trim(adjustl(base))//'_q-'//trim(adjustl(chrq(3)))//'_zsat-'  &
    & //trim(adjustl(chrz(3)))
#endif

call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))

open(65, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(o2flx).txt', &
    & status='replace')
open(67, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe2flx).txt', &
    & status='replace')
open(68, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe3flx).txt', &
    & status='replace')
open(69, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(pyflx).txt',  &
    & status='replace')
open(55, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx).txt',  &
    & status='replace')
open(56, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(silflx).txt',  &
    & status='replace')
open(57, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(so4flx).txt',  &
    & status='replace')
open(58, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(mgflx).txt',  &
    & status='replace')
open(59, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(siflx).txt',  &
    & status='replace')
open(60, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(foflx).txt',  &
    & status='replace')
open(61, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx2).txt',  &
    & status='replace')
open(62, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(abflx).txt',  &
    & status='replace')
open(63, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(co2flx).txt', &
    & status='replace')
open(64, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(anflx).txt',  &
    & status='replace')
open(54, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(caflx).txt', &
    & status='replace')
open(53, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(ccflx).txt', &
    & status='replace')
open(52, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'rain.txt', &
    & status='replace')


close(65);close(67);close(68);close(69);close(55);close(56);close(57);close(58);close(59);close(60)
close(61);close(62);close(63);close(64);close(54);close(53);close(52)



!!!  MAKING GRID !!!!!!!!!!!!!!!!! 
beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
beta = 1.00005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z)


if (swoxall==1d0) swoxa = 0d0

stoxs = stoxs - swoxa*stoxa
#ifdef poroevol 
stoxs = stoxs - (1d0-swoxall)*stoxa  ! 15/4 if overall oxidation is considered (swoxall = 1)
                           ! 7/2  if pyrite oxidation and aqFe2+ oxidation is separately considered (swoxall = 0)
#endif 

do irec=1,nrec
    reclis(irec)=irec*nt/nrec
end do

po2 = po2i
pco2 = pco2i
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
ca = 0d0
pro = 10.0d0**(-ph)

resp = 0d0

dt = maxdt
if (.not.initial_ss) then 
    dt = 1d-2  
else
    dt = 1d-300 ! for basalt exp?
endif 

pro = 1d-5
! call calc_pH( &
    ! & nz,na+2d0*(mg+ca-so4),pco2,kw,kco2,k1,k2 &! input 
    ! & ,pro &! output
    ! & ) 
! print*,-log10(pro)
! stop

call coefs( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
    & ,pro,dgas,dgasc,daq,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
    & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
    & )
    
call calc_pH( &
    & nz,na+2d0*(mg+ca-so4),pco2,kw,kco2,k1,k2 &! input 
    & ,pro &! output
    & )


!       print *,keqsil;stop


redsld = redsldi*1d-2/120.0d0*mvpy*2.7d0*(1.0d0-poro)  
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

! mfoi = mfoth*0.1d0
! mani = manth*0.1d0

mfoi = 1d-10
mani = 1d-10
mcci = 1d-10

mfo = mfoi
man = mani
mcc = mcci

poroprev = poro

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
    runname_save = 'test_cpl_rain-0.50E+02_pevol_sevol1_q-0.10E+00_zsat-5' ! specifiy the file where restart data is stored 
    ! runname_save = runname  ! the working folder has the restart data 
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'o2profile-res-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt')
    open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt',  &
        & status ='old',action='read')

    do iz = 1, Nz
        read (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),mab(iz),man(iz),mfo(iz),mcc(iz) &
            & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_cc(iz),pco2(iz),pro(iz),time
    enddo 
    close(22) 
    pro = 10d0**(-pro) ! read data is -log10 (pro)
    time = 0d0
    
    ! do iz=1,nz
        ! if (po2(iz)<po2th) po2(iz)= po2th*0.1d0
        ! if (c(iz)<cth) c(iz)= cth*0.1d0
        ! if (ms(iz)<msth) ms(iz)= msth*0.1d0
        ! if (c2(iz)<c2th) c2(iz)= c2th*0.1d0
        ! if (so4(iz)<so4th) so4(iz)= so4th*0.1d0
        ! if (na(iz)<nath) na(iz)= nath*0.1d0
        ! if (ca(iz)<cath) ca(iz)= cath*0.1d0
        ! if (mg(iz)<mgth) mg(iz)= mgth*0.1d0
        ! if (si(iz)<sith) si(iz)= sith*0.1d0
        ! if (mab(iz)<mabth) mab(iz)= mabth*0.1d0
        ! if (man(iz)<manth) man(iz)= manth*0.1d0
        ! if (mfo(iz)<mfoth) mfo(iz)= mfoth*0.1d0
        ! if (mcc(iz)<mccth) mcc(iz)= mccth*0.1d0
        ! if (pco2(iz)<pco2th) pco2(iz)= pco2th*0.1d0
    ! enddo
    
    ! man = mani
    ! mfo = mfoi
    ! manx = man
    ! mfox = mfo
    mcc = mcci ! to precipirate calcite, crystal seeds are necessary 
    mccx = mcc
    
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'o2profile-bsd-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-save.txt')
    open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-save.txt',  &
        & status ='old',action='read')
    do iz=1,nz
        read(22,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
    enddo 
    close(22) 
    torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
    tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
        
#ifdef display      
    Print *,'==== printing read data  ===='
    print *
    print *,'-=-=-=-=-=-= o2 & pyrite -=-=-=-=-=-=-='
    print *,'o2:', (po2(iz),iz=1,nz, nz/5)
    print *,'fe2:', (c(iz),iz=1,nz, nz/5)
    print *,'py:', (ms(iz),iz=1,nz, nz/5)
    print *, 'fe3:', (c2(iz),iz=1,nz, nz/5)
    print *, 'so4:', (so4(iz),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An -=-=-=-=-=-=-='
    print *, 'mg:', (mg(iz),iz=1,nz, nz/5)
    print *, 'si:', (si(iz),iz=1,nz, nz/5)
    print *, 'na:', (na(iz),iz=1,nz, nz/5)
    print *, 'ca:', (ca(iz),iz=1,nz, nz/5)
    print *, 'fo:', (mfo(iz),iz=1,nz, nz/5)
    print *, 'ab:', (mab(iz),iz=1,nz, nz/5)
    print *, 'an:', (man(iz),iz=1,nz, nz/5)
    print *, 'cc:', (mcc(iz),iz=1,nz, nz/5)
    print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(pro(iz)),iz=1,nz, nz/5)
    print *
#endif      
endif

call coefs( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
    & ,pro,dgas,dgasc,daq,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
    & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
    & ) 
    
! --------- loop -----
it = 0
irec = 0

!! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@

do while (it<nt)
#ifdef display 
    print *, 'it, time = ',it, time
#endif
    ! if (time>rectime(nrec)) exit
    if (initial_ss.and.time>rectime(nrec)) exit
        
    if (.not.initial_ss .and. time > ztot/w*2d0) then 
        initial_ss = .true.
        time = 0
        dt = 1d-300
        ! pause
        
        ! man = mani
        ! mfo = mfoi
        mcc = mcci
        ! manx = man
        ! mfox = mfo
        mccx = mcc
    endif 

    if (.not.initial_ss) then 
        ! if (time > ztot/w*2d0 *0.99d0) then 
            ! maxdt = 1d0
        ! else  
            maxdt = 1d2
        ! endif 
    else 
        ! maxdt = 0.2d0
        maxdt = 0.02d0 ! when calcite is included smaller time step must be assumed 
        maxdt = 0.005d0 ! when calcite is included smaller time step must be assumed 
        maxdt = 0.002d0 ! when calcite is included smaller time step must be assumed 
    endif 

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
        ! dt = dt*1.01d0
        dt = dt*10d0
        if (dt>maxdt) dt = maxdt
    endif
    if (iter > 300) then
        ! dt = dt/1.05d0
        dt = dt/10d0
    end if

    ! incase temperature&ph change
    call coefs( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
        & ,pro,dgas,dgasc,daq,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
        & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
        & ) 

    ! ======== modifying maxdt ===============
    !       if (time >1d4) maxdt = 10d0
    ! ========================================
    ! #ifdef pHiter      
     ! ############## pH interation #######################

    ! error2 = 1d4
    ! iter2 = 0
    ! do while (error2 > tol ) 
    ! #endif


    error = 1d4
    iter=0

100 print *, iter,time

    po2x = po2
    pco2x = pco2
    cx = c
    msx = ms

    c2x = c2

    so4x=so4

    nax = na
    msilx = msil

    six = si
    mgx = mg
    cax = ca
    mfox = mfo

    mabx = mab
    manx = man
    mccx = mcc

    prox = pro  

    porox = poro


    if (initial_ss) then 
        mfosupp = rainpowder*rainfrc_fo/mwtfo*exp(-z/zsupp)/zsupp
        mabsupp = rainpowder*rainfrc_ab/mwtab*exp(-z/zsupp)/zsupp 
        mansupp = rainpowder*rainfrc_an/mwtan*exp(-z/zsupp)/zsupp
        mccsupp = 0d0
        ! kcc = k_arrhenius(10d0**(-5.81d0)*sec2yr,25d0+tempk_0,tc+tempk_0,23.5d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
        ! mccsupp = kcc*poro*hr*mvcc*1d-6*mccx*(cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc - 1d0) &
            ! & *merge(1d0,0d0,cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc - 1d0 > 0d0) 
            
        ! kcc = 0d0
        if (rain_wave) then 
            mfosupp = mfosupp*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau)) 
            mabsupp = mabsupp*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau)) 
            mansupp = mansupp*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
            if (time==0d0 .or. rain_norm /= merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))) then
                open(52, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'rain.txt', &
                    & status='old',action='write',access='append')
                write(52,*) time-dt,rain_norm
                write(52,*) time,merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                rain_norm = merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                close(52)
            endif 
        endif 
    else 
        mfosupp = 0d0
        mabsupp = 0d0
        mansupp = 0d0
        mccsupp = 0d0
    endif 

    msilsupp = 0d0


    if (pre_calc) then 
        
        call precalc_po2_v2( &
            & nz,po2th,dt,ucv,kho,dz,dgas,daq,po2i,poro,sat,po2,torg,tora,v &! input 
            & ,po2x &! output 
            & )
        
        call precalc_pco2_v2( &
            & nz,pco2th,dt,ucv,khco2,dz,dgasc,daqc,pco2i,poro,sat,pco2,torg,tora,v,resp &! input 
            & ,pco2x &! output 
            & )
        
        call precalc_slds( &
            & nz,msth,dt,w,dz,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth   &! input
            & ,ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp &! input
            & ,msx,msilx,mfox,mabx,manx,mccx &! output
            & )

        ! pause

        if (swoxa == 1d0) then 
        
            call precalc_pw_py( &
                & nz,dt,cth,c2th,v,c2,dz,dfe3,poro,sat,tora,koxa,po2x,c,dfe2,koxs,koxs2,hr,mvpy &! input
                & ,ms,po2,so4,so4th,ci,c2i,so4i,dso4,msx &! input 
                & ,c2x,cx,so4x  &! output
                & )
            
        endif 
        
        call precalc_pw_sil_v2( &
            & nz,nath,mgth,cath,sith,dt,v,na,ca,mg,si,dz,dna,dsi,dmg,dca,tora,poro,sat,nai,mgi,cai,sii &! input 
            & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx &! input 
            & ,nax,six,cax,mgx &! output
            & )

        if (any(isnan(po2x)).or.any(isnan(cx)).or.any(isnan(c2x)).or.any(isnan(pco2x))) then 
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
        cax(:) = 1.0d2
    end if

    call pyweath_1D( &
        & nz,nflx_py,mvpy,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
        & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daq,dgas,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
        & ,swoxa,swoxall,ucv,vmax  &! inpput
        & ,iter,error,dt &! inout
        & ,cx,c2x,po2x,msx,flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2 &! output
        ) 
    !!!!!!!!!!!!!!!!!!!!!!!!  so4 calculation start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call pyweath_1D_SO4( &
        & nz,nflx_py,mvpy,c,c2,po2,ms,hr,po2th,poro,z,dz,koxs2,koxs,dso4,sat,dporodta  &! input
        & ,cth,tora,v,tol,zrxn,dt,cx,c2x,po2x,msx,so4,swoxa,O2_evolution,so4i,so4th &! input
        & ,so4x,flx_py_so4 &! output
        & )

    !!!  =-=-=-=-=-=-=-=-=- END of SO4 calculation  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#endif       
#ifndef pyweath       
    !!! =-=-=-=-=-=-=-=-=- START of calculation for Na and albite  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#ifdef silweath
    so4x = so4th
    cx = cth 
    c2x = c2th 

    if ((.not.read_data) .and. it == 0 .and. iter == 0) then 
        nax(1:) = 1.0d2
        mgx(1:) = 1.0d2
        six(1:) = 1.0d2
        cax(1:) = 1.0d2
    endif 

#endif      

    error_co2 = 1d4
    iter_co2 = 0 
    do while (error_co2> 1e-3) 
        ! call silicate_dis_1D( &
            ! & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
            ! & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp  &! input
            ! & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
            ! & ,iter,error,dt,flgback &! inout
            ! & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
            ! & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc &! output
            ! & )
        call silicate_dis_1D_v2( &
            & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
            & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
            & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
            & ,iter,error,dt,flgback &! inout
            & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
            & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc &! output
            & )

        ! call oxygen_resp_1D( &
            ! & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
            ! & ,kho,tora,torg,daq,dgas,v,mo2,tol,runname,workdir,zrxn,ucv,vmax  &! inpput
            ! & ,iter,error,dt &! inout
            ! & ,po2x,flx_o2,resp &! output
            ! & ) 
            
        call oxygen_resp_1D_v2( &
            & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
            & ,kho,tora,torg,daq,dgas,v,mo2,tol,runname,workdir,zrxn,ucv,vmax,poroprev  &! inpput
            & ,dt,flgback &! inout
            & ,po2x,flx_o2,resp &! output
            & ) 

        pco2x_prev = pco2x
        preccc = flx_cc(4,:)
        ! preccc = 0d0
        ! call CO2_1D( &
            ! & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
            ! & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,flx_cc(4,:)  &! inpput
            ! & ,iter,error,dt &! inout
            ! & ,pco2x,flx_co2,khco2 &! output
            ! & )
        call CO2_1D_v2( &
            & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
            & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,preccc,poroprev,pro  &! inpput
            & ,dt,flgback &! inout
            & ,pco2x,flx_co2 &! output
            & ) 
        ! call CO2_1D_v3( &
            ! & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
            ! & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,flx_cc(4,:),poroprev,pro  &! inpput
            ! & ,dt &! inout
            ! & ,pco2x,flx_co2 &! output
            ! & ) 
        error_co2 = 0d0
        if (co2_iteration) then
            do iz = 1,nz
                if (pco2x(iz)>pco2th.and.pco2x_prev(iz)>pco2th) then 
                    error_co2 = max(error_co2,abs((pco2x(iz)-pco2x_prev(iz))/pco2x(iz)))
                endif 
            enddo
        endif 
        iter_co2 = iter_co2 + 1
        
        print *, iter_co2,error_co2
        if (iter_co2 > 100) then 
            dt = dt/10d0
            go to 100
        endif 
    enddo 
    
    ! call silicate_dis_co2_1D( &
        ! & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
        ! & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,resp,poroprev,daqc,dgasc  &! input
        ! & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
        ! & ,khco2i,pco2i,pco2th,torg,ucv &! input
        ! & ,iter,error,dt,flgback &! inout
        ! & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
        ! & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,pco2x,flx_co2,khco2x &! output
        ! & )

    if (flgback) then 
        flgback = .false. 
        go to 100
    endif    

    dporodtg = 0d0
    dporodtgc = 0d0
    dporodta = 0d0
#ifdef poroevol   
    poroprev = poro
    torgprev = torg
    toraprev = tora
    poro = poroi + (mabi-mabx)*(mvab -  0.5d0*99.52d0)*1d-6  &
        & +(mfoi-mfox)*(mvfo)*1d-6 &
        & +(mani-manx)*(mvan - mvkaol)*1d-6 &
        & +(mcci-mccx)*(mvcc )*1d-6 
    if (any(poro<0d0)) then 
        print*,'negative porosity: stop'
        print*,poro
        ! w = w*2d0
        ! go to 100
        stop
    endif 
    v = qin/poro/sat
    torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
    tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
    deff = torg*dgas + tora*daq
    dporodtg = ( &
        & (ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3) &
        & -(ucv*porox*(1.0d0-sat)*1d3+porox*sat*kho*1d3) &
        & )/dt
    dporodtgc = ( &
        & (ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2*1d3) &
        & -(ucv*porox*(1.0d0-sat)*1d3+porox*sat*khco2*1d3) &
        & )/dt
    dporodta = (poro*sat-porox*sat)/dt/(poro*sat)
    hr = hri
#ifdef surfevol1 
    hr = hri*((1d0-poro)/(1d0-poroi))**(2d0/3d0)
#endif 
#ifdef surfevol2 
    hr = hri*(poro/poroi)**(2d0/3d0)  ! SA increases with porosity 
#endif 
#ifdef display
    print *
    print *,'-=-=-=-=-=-= Porosity & SA -=-=-=-=-=-=-='
    print *, 'phi:', (poro(iz),iz=1,nz, nz/5)
    print *, 'SA:', (hr(iz),iz=1,nz, nz/5)
    print *
#endif 
#endif  

    ! #ifdef pHiter
    ! error2 = maxval( abs(pro-prox)/pro)
    ! #ifdef display
    ! print*,"=========pH iteration========"
    ! print*,"iter2,error2",iter2,error2
    ! print*,"============================="
    ! #endif
    ! iter2 = iter2 + 1

    ! pro = prox

    ! if (iter2 > 300) then
        ! dt = dt/10d0
        ! if (dt==0d0) then 
            ! print *. 'dt==0d0; stop'
            ! stop
        ! endif 
    ! end if

    ! enddo     
    ! ######################## end of pH iteration #######################
    ! #endif

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
    print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An -=-=-=-=-=-=-='
    print *, 'mg:', (mgx(iz),iz=1,nz, nz/5)
    print *, 'si:', (six(iz),iz=1,nz, nz/5)
    print *, 'na:', (nax(iz),iz=1,nz, nz/5)
    print *, 'ca:', (cax(iz),iz=1,nz, nz/5)
    print *, 'fo:', (mfox(iz),iz=1,nz, nz/5)
    print *, 'ab:', (mabx(iz),iz=1,nz, nz/5)
    print *, 'an:', (manx(iz),iz=1,nz, nz/5)
    print *, 'cc:', (mccx(iz),iz=1,nz, nz/5)
    print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= CO2 -=-=-=-=-=-=-='
    print *, 'co2:', (pco2x(iz),iz=1,nz, nz/5)
    print *
#endif 

    ! stop

    po2 = po2x
    pco2 = pco2x
    c = cx
    ms = msx
    c2 = c2x
    so4 = so4x
    na = nax
    ca = cax
    msil = msilx
    mg = mgx
    si = six
    mfo = mfox
    mab = mabx
    man = manx
    mcc = mccx
    pro = prox
        
    do iflx = 1,nflx_py
        flx_py(iflx,:) = flx_py(iflx,:)*dz
        flx_py_o2(iflx,:) = flx_py_o2(iflx,:)*dz
        flx_py_fe2(iflx,:) = flx_py_fe2(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_py_fe3(iflx,:) = flx_py_fe3(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_py_so4(iflx,:) = flx_py_so4(iflx,:)*dz*poro(:)*sat(:)*1d3
    enddo
        
    do iflx = 1,nflx
        flx_mg(iflx,:) = flx_mg(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_si(iflx,:) = flx_si(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_na(iflx,:) = flx_na(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_ca(iflx,:) = flx_ca(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_fo(iflx,:) = flx_fo(iflx,:)*dz
        flx_ab(iflx,:) = flx_ab(iflx,:)*dz
        flx_an(iflx,:) = flx_an(iflx,:)*dz
        flx_cc(iflx,:) = flx_cc(iflx,:)*dz
        flx_o2(iflx,:) = flx_o2(iflx,:)*dz
        flx_co2(iflx,:) = flx_co2(iflx,:)*dz
    enddo 

    if ((.not.initial_ss) .and. time > savetime) then 
        open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt',  &
            & status='replace')

        open(30, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-save.txt',  &
            & status='replace')

        do iz = 1, Nz
            write(22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),mab(iz),man(iz),mfo(iz),mcc(iz) &
                & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_cc(iz),pco2(iz),-log10(pro(iz)),time
            write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
        enddo 
        close(22)
        close(30)
        savetime = savetime + dsavetime
    endif 

    if (initial_ss .and. time>=rectime(irec+1)) then
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
                & koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*po2x(iz)**0.50d0*merge(0.0d0,1.0d0,po2x(iz)<po2th), &
                & merge(0.0d0, &
                & + koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th), &
                & poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth) &
                & , swbr*vmax*po2x(iz)/(po2x(iz)+mo2) &
                & ,kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
                & ,kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                & ,time
            write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),mab(iz),man(iz),mfo(iz),mcc(iz) &
            & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_cc(iz),pco2(iz),-log10(pro(iz)),time
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

        open(65, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(o2flx).txt' &
            & ,action='write',status='old',access='append')
        open(67, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe2flx).txt' &
            & ,action='write',status='old',access='append')
        open(68, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(fe3flx).txt' &
            & ,action='write',status='old',access='append')
        open(69, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(pyflx).txt'  &
            & ,action='write',status='old',access='append')
        open(55, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx).txt'  &
            & ,action='write',status='old',access='append')
        open(56, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(silflx).txt'  &
            & ,action='write',status='old',access='append')
        open(57, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(so4flx).txt'  &
            & ,action='write',status='old',access='append')
        open(58, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(mgflx).txt'  &
            & ,action='write',status='old',access='append')
        open(59, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(siflx).txt'  &
            & ,action='write',status='old',access='append')
        open(60, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(foflx).txt'  &
            & ,action='write',status='old',access='append')
        open(61, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(naflx2).txt'  &
            & ,action='write',status='old',access='append')
        open(62, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(abflx).txt'  &
            & ,action='write',status='old',access='append')
        open(63, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(co2flx).txt' &
            & ,action='write',status='old',access='append')
        open(64, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(anflx).txt'  &
            & ,action='write',status='old',access='append')
        open(54, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(caflx).txt' &
            & ,action='write',status='old',access='append')
        open(53, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(ccflx).txt' &
            & ,action='write',status='old',access='append')

        ! write(65,*) time, o2tflx, diflx, advflx, pyoxflx, feoxflx, respflx,o2flxsum 
        write(65,*) time,(sum(flx_o2(iflx,:)),iflx=1,nflx) 
        write(67,*) time,(sum(flx_py_fe2(iflx,:)),iflx=1,nflx_py) 
        write(68,*) time,(sum(flx_py_fe3(iflx,:)),iflx=1,nflx_py)  
        write(69,*) time,(sum(flx_py(iflx,:)),iflx=1,nflx_py)         
        write(56,*) time,(sum(flx_ab(iflx,:)),iflx=1,nflx)       
        write(55,*) time,(sum(flx_na(iflx,:)),iflx=1,nflx)      
        write(57,*) time,(sum(flx_py_so4(iflx,:)),iflx=1,nflx_py) 
        
        write(58,*) time,(sum(flx_mg(iflx,:)),iflx=1,nflx)
        write(59,*) time,(sum(flx_si(iflx,:)),iflx=1,nflx)
        write(60,*) time,(sum(flx_fo(iflx,:)),iflx=1,nflx)
        write(61,*) time,(sum(flx_na(iflx,:)),iflx=1,nflx)
        write(62,*) time,(sum(flx_ab(iflx,:)),iflx=1,nflx)
        write(63,*) time,(sum(flx_co2(iflx,:)),iflx=1,nflx)
        write(64,*) time,(sum(flx_an(iflx,:)),iflx=1,nflx)
        write(54,*) time,(sum(flx_ca(iflx,:)),iflx=1,nflx)
        write(53,*) time,(sum(flx_cc(iflx,:)),iflx=1,nflx)
        
        close(65);close(67);close(68);close(69);close(55);close(56);close(57);close(58);close(59);close(60)
        close(61);close(62);close(63);close(64);close(54);close(53)

    end if

    it = it + 1
    time = time + dt

end do

write(chr,'(i3.3)') irec+1
open (28, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(rate)-'//chr//'.txt', &
    & status='replace')
open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-'//chr//'.txt', &
    & status='replace')
open (29, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res(nom)-'//chr//'.txt', &
    & status='replace')

open(30, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-'//chr//'.txt',  &
    & status='replace')

do iz = 1, Nz
    write (28,*) z(iz), &
        & koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*po2x(iz)**0.50d0*merge(0.0d0,1.0d0,po2x(iz)<po2th), &
        & merge(0.0d0, &
        & + koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),cx(iz)<cth) &
        & *merge(0.0d0,1.0d0,c2x(iz)<c2th), &
        & poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz)*po2x(iz)*merge(0.0d0,1.0d0,po2x(iz)<po2th.or.cx(iz)<cth) &
        & , swbr*vmax*po2x(iz)/(po2x(iz)+mo2) &
        & ,kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
        & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
        & ,kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
        & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
        & ,time
    write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),mab(iz),man(iz),mfo(iz),mcc(iz) &
        & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_cc(iz),pco2(iz),-log10(pro(iz)),time
    write (29,*) z(iz),po2(iz)/maxval(po2(:)),c(iz)/maxval(c(:)),ms(iz)/maxval(ms(:)),c2(iz)/maxval(c2(:)) &
        & ,so4(iz)/maxval(so4(:)),na(iz)/maxval(na(:)),msil(iz)/maxval(msil(:)),pro(iz)/maxval(pro(:)) &
        & ,silsat(iz), co2(iz),hco3(iz),co3(iz),dic(iz),time
    write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
end do

close(28)
close(22)
close(29)
close(30)

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

#ifdef sense
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
#endif 

endprogram weathering

!**************************************************************************************************************************************
subroutine makegrid(beta,nz,ztot,dz,z)  !  making grid, after Hoffmann & Chiang, 2000
implicit none
integer(kind=4),intent(in) :: nz
real(kind=8),intent(in)::beta,ztot
real(kind=8),intent(out)::dz(nz),z(nz)
integer(kind=4) iz

do iz = 1, nz 
    z(iz) = iz*ztot/nz  ! regular grid 
    if (iz==1) then
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0))
    endif
    if (iz/=1) then 
        dz(iz) = ztot*log((beta+(z(iz)/ztot)**2d0)/(beta-(z(iz)/ztot)**2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(:iz-1))
    endif
enddo

#ifdef regulargrid
dz = ztot/nz  ! when implementing regular grid
#endif 

do iz=1,nz  ! depth is defined at the middle of individual layers 
    if (iz==1) z(iz)=dz(iz)*0.5d0  
    if (iz/=1) z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz)
enddo

endsubroutine makegrid
!**************************************************************************************************************************************

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine coefs( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
    & ,pro,dgas,dgasc,daq,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
    & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2th
real(kind=8),dimension(nz),intent(in)::pro,pco2
real(kind=8),intent(out)::dgas,dgasc,daq,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,keqsil &
    & ,keqab,keqfo,keqcc,keqan,ucv
real(kind=8),dimension(nz),intent(out)::ksil,kab,kfo,kan,kcc,koxs,koxa,koxs2,khco2

real(kind=8),dimension(nz):: koxsi,koxai,koxs2i
real(kind=8) k_arrhenius


ucv = 1.0d0/(rg2*(tempk_0+tc))

dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dmg  = 0.017218079d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dsi  = 0.03689712d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dca  = 0.019023312d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 

dgas = 6.09d2 ! m^2 yr^-1
daq = 5.49d-2 ! m^2 yr^-1
dgasc = 441.504d0 ! m^2 yr^-1 (Assuming 0.14 cm2/sec)
daqc = 0.022459852d0 ! m^2 yr^-1 (for C32- from Li and Gregory 1974)

dgas = dgas*exp(-4.18d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
daq = daq*exp(-20.07d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dfe2=dfe2*exp(-19.615251d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dfe3=dfe3*exp(-14.33659d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dso4=dso4*exp(-20.67364d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dna=dna*exp(-20.58566d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dmg=dmg*exp(-18.51979d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)

kho=10.0d0**(-2.89d0)*exp(13.2d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+25.0d0))/rg)

kco2 = 10.0d0**(-1.34d0)  ! 15C Kanzaki Murakami 2015
k1 = 10.0d0**(-6.42d0)  ! 15C Kanzaki Murakami 2015
k2 = 10d0**(-10.43d0)      ! 15C Kanzaki Murakami 2015
kw = 10.0d0**(-14.35d0)  ! 15C Kanzaki Murakami 2015

if (tc==5d0) then

    kco2 = 10.0d0**(-1.19d0)  ! 15C Kanzaki Murakami 2015
    k1 = 10.0d0**(-6.52d0)  ! 15C Kanzaki Murakami 2015
    k2 = 10d0**(-10.55d0)      ! 15C Kanzaki Murakami 2015
    kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0+0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
    kw =10d0**kw

elseif (tc==25d0) then 

    kco2 = 10.0d0**(-1.47d0)  ! 15C Kanzaki Murakami 2015
    k1 = 10.0d0**(-6.35d0)  ! 15C Kanzaki Murakami 2015
    k2 = 10d0**(-10.33d0)      ! 15C Kanzaki Murakami 2015
    kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0 &
        & +0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
    kw =10d0**kw

endif 


khco2i = kco2*(1d0+k1/sqrt(kco2*k1*pco2i)+k2/kco2/pco2i)
khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro)

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

! kan = 10d0**(-9.12d0)*60d0*60d0*24d0*365d0 &! mol/m2/yr  
    ! & *exp(-17.8d0/rg*(1d0/(tc+273d0)-1d0/(25d0+273d0)))
kan = k_arrhenius(10d0**(-9.12d0)*sec2yr,25d0+tempk_0,tc+tempk_0,17.8d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)

! anorthite (CaAl2Si2O8) + 2H+ + H2O --> kaolinite(Al2Si2O5(OH)4) + Ca2+ 
keqan = 28.8615308d0 - 8.310989613d0
keqan = 10d0**keqan

kcc = k_arrhenius(10d0**(-5.81d0)*sec2yr,25d0+tempk_0,tc+tempk_0,23.5d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
! kcc = kcc**merge(0.0d0,1.0d0,pco2<pco2th)
! kcc = 0d0

keqcc = 10d0**(-8.43d0) ! Kanzaki and Murakami 2015

if (tc==5d0) then 
    ksil =  5.13d-10*1d4  ! mol/m2/yr  ! from Li et al., 2014
    keqsil = 3.751169218d0 - 0.5d0* 9.242748918d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
    keqsil = 10.0d0**(keqsil)
elseif (tc==25d0) then
    ksil =  3.15d-9*1d4  ! mol/m2/yr  ! from Li et al., 2014
    keqsil = 3.080792422d0 - 0.5d0* 7.434511289d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
    keqsil = 10.0d0**(keqsil)
endif       

koxsi = 10.0d0**(-8.19d0)*60.0d0*60.0d0*24.0d0*365.0d0  &!! excluding the term (po2**0.5)
    & *(kho)**(0.50d0)/(pro**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)

koxai = max(swoxa*8.0d13*60.0d0*24.0d0*365.0d0   &!  excluding the term (c*po2)
    & *(kw/pro)**2.0d0               &! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)
    & , swoxa*1d-7*60.0d0*24.0d0*365.0d0)

if (tc/=15d0) then 
    koxsi = koxsi*exp(-57d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
endif 

koxs = koxsi
koxa = koxai

koxs2 = koxs2i

endsubroutine coefs

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_po2( &
    & nz,po2th,dt,ucv,kho,dz,dgas,daq,po2i,poro,sat,po2,torg,tora,v &! input 
    & ,po2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::po2th,dt,ucv,kho,dz,dgas,daq,po2i
real(kind=8),dimension(nz),intent(in)::poro,sat,po2,torg,tora,v
real(kind=8),dimension(nz),intent(out)::po2x

integer iz
real(kind=8) po2tmp

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
            & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz &
            & ) &
            & )
    else ! iz == 1
        po2x(iz) = max(0.0d0 &
            & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))* &
            & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3) &
            & *(-po2(iz))/dt-(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(-2.0d0*po2(iz) + po2(iz+1)+po2i) &
            & /(dz**2.0d0)+poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2i)/dz  &
            & ) &
            & )
    endif

end do 

endsubroutine precalc_po2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_po2_v2( &
    & nz,po2th,dt,ucv,kho,dz,dgas,daq,po2i,poro,sat,po2,torg,tora,v &! input 
    & ,po2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::po2th,dt,ucv,kho,dgas,daq,po2i
real(kind=8),dimension(nz),intent(in)::poro,sat,po2,torg,tora,v,dz
real(kind=8),dimension(nz),intent(out)::po2x

integer iz
real(kind=8) po2tmp,edifi,ediftmp
real(kind=8),dimension(nz)::alpha,edif

alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3
edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgas +poro*sat*kho*1d3*tora*daq
edifi = edif(1)
edifi = ucv*1d3*dgas 

do iz = 1, nz

    if (po2x(iz)>=po2th) cycle

    po2tmp = po2(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) po2tmp = po2i
    if (iz==1) ediftmp = edifi

    po2x(iz) = max(0.0d0 &
        & ,dt/alpha(iz)* &
        & ( &
        & alpha(iz)*po2(iz)/dt &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(po2(min(nz,iz+1))-po2(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(po2(iz)-po2tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & - poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2tmp)/dz(iz)  &
        & ) &
        & )

end do 

endsubroutine precalc_po2_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_pco2( &
    & nz,pco2th,dt,ucv,khco2,dz,dgasc,daqc,pco2i,poro,sat,pco2,torg,tora,v,resp &! input 
    & ,pco2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::pco2th,dt,ucv,dz,dgasc,daqc,pco2i
real(kind=8),dimension(nz),intent(in)::khco2,poro,sat,pco2,torg,tora,v,resp
real(kind=8),dimension(nz),intent(out)::pco2x

integer iz
real(kind=8) pco2tmp

do iz = 1, nz

    if (pco2x(iz)>=pco2th) cycle

    if (iz/=nz) pco2tmp = pco2(iz+1)
    if (iz==nz) pco2tmp = pco2(iz)

    if (iz/=1) then 
        pco2x(iz) = max(0.0d0  &
            & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3))* &
            & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3)*(-pco2(iz))/dt &
            & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
            & *(pco2tmp+pco2(iz-1)-2.0d0*pco2(iz))/(dz**2.0d0) &
            & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
            & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz)*tora(iz-1)*daqc)) &
            & *(pco2(iz)-pco2(iz-1))/(dz**2.0d0) &
            & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2(iz)-pco2(iz-1))/dz &
            & ) &
            & +resp(iz)*dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) &
            & )
    else ! iz == 1
        pco2x(iz) = max(0.0d0 &
            & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3))* &
            & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) &
            & *(-pco2(iz))/dt-(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
            & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)*(-2.0d0*pco2(iz) + pco2(iz+1)+pco2i) &
            & /(dz**2.0d0)+poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2(iz)-pco2i)/dz  &
            & ) &
            & +resp(iz)*dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) &
            & )
    endif

end do 

endsubroutine precalc_pco2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_pco2_v2( &
    & nz,pco2th,dt,ucv,khco2,dz,dgasc,daqc,pco2i,poro,sat,pco2,torg,tora,v,resp &! input 
    & ,pco2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::pco2th,dt,ucv,dgasc,daqc,pco2i
real(kind=8),dimension(nz),intent(in)::khco2,poro,sat,pco2,torg,tora,v,resp,dz
real(kind=8),dimension(nz),intent(out)::pco2x

integer iz
real(kind=8) pco2tmp,edifi,ediftmp
real(kind=8),dimension(nz)::alpha,edif

alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2*1d3
edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc +poro*sat*khco2*1d3*tora*daqc
edifi = edif(1)
edifi = ucv*1d3*dgasc 

do iz = 1, nz

    if (pco2x(iz)>=pco2th) cycle

    pco2tmp = pco2(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) pco2tmp = pco2i
    if (iz==1) ediftmp = edifi

    pco2x(iz) = max(0.0d0 &
        & ,dt/alpha(iz)* &
        & ( &
        & alpha(iz)*pco2(iz)/dt &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(pco2(min(nz,iz+1))-pco2(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(pco2(iz)-pco2tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & - poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2(iz)-pco2tmp)/dz(iz)  &
        & + resp(iz) & 
        & ) &
        & )

end do 

endsubroutine precalc_pco2_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds( &
    & nz,msth,dt,w,dz,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth   &! input
    & ,ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp &! input
    & ,msx,msilx,mfox,mabx,manx,mccx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::msth,dt,w,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth
real(kind=8),dimension(nz),intent(in)::ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp,dz
real(kind=8),dimension(nz),intent(out)::msx,msilx,mfox,mabx,manx,mccx

integer iz

do iz = 1, nz

    if (msx(iz)>=msth) cycle

    if (iz/=nz) then 
        msx(iz) = max(0d0, &
            & ms(iz) +dt*(w*(ms(iz+1)-ms(iz))/dz(iz)) &
            & )
    else  
        msx(iz) = max(0d0, &
            & ms(iz) +dt*(w*(msi-ms(iz))/dz(iz)) &
            & )
    endif

    if (msilx(iz)>=msilth) cycle

    if (iz/=nz) then 
        msilx(iz) = max(0d0, &
            & msil(iz) +dt*(w*(msil(iz+1)-msil(iz))/dz(iz) + msilsupp(iz)) &
            & )
    else 
        msilx(iz) = max(0d0, &
            & msil(iz) + dt*(w*(msili-msil(iz))/dz(iz) + msilsupp(iz)) &
            & )
    endif 

    if (mfox(iz)>=mfoth) cycle

    if (iz/=nz) then 
        mfox(iz) = max(0d0, &
            & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
            & )
    else 
        mfox(iz) = max(0d0, &
            & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
            & )
    endif 

    if (mabx(iz)>=mabth) cycle

    if (iz/=nz) then 
        mabx(iz) = max(0d0, &
            & mab(iz) +dt*(w*(mab(iz+1)-mab(iz))/dz(iz) + mabsupp(iz)) &
            & )
    else 
        mabx(iz) = max(0d0, &
            & mab(iz) + dt*(w*(mabi-mab(iz))/dz(iz)+ mabsupp(iz)) &
            & )
    endif 

    if (manx(iz)>=manth) cycle

    if (iz/=nz) then 
        manx(iz) = max(0d0, &
            & man(iz) +dt*(w*(man(iz+1)-man(iz))/dz(iz) + mansupp(iz)) &
            & )
    else 
        manx(iz) = max(0d0, &
            & man(iz) + dt*(w*(mani-man(iz))/dz(iz)+ mansupp(iz)) &
            & )
    endif 

    if (mccx(iz)>=mccth) cycle

    if (iz/=nz) then 
        mccx(iz) = max(0d0, &
            & mcc(iz) +dt*(w*(mcc(iz+1)-mcc(iz))/dz(iz) + mccsupp(iz) &
            ! & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0) &
            ! & *merge(1d0,0d0,cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0 > 0d0) &
            & ) &
            & )
    else 
        mccx(iz) = max(0d0, &
            & mcc(iz) + dt*(w*(mcci-mcc(iz))/dz(iz)+ mccsupp(iz)&
            ! & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0) &
            ! & *merge(1d0,0d0,cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0 > 0d0) &
            & ) &
            & )
    endif 

enddo

endsubroutine precalc_slds

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds_wopy( &
    & nz,dt,w,dz,mfoi,mabi,mani,mcci,mabth,manth,mfoth,mccth   &! input
    & ,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp &! input
    & ,mfox,mabx,manx,mccx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w,mfoi,mabi,mani,mcci,mabth,manth,mfoth,mccth
real(kind=8),dimension(nz),intent(in)::mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp,dz
real(kind=8),dimension(nz),intent(out)::mfox,mabx,manx,mccx

integer iz

do iz = 1, nz

    if (mfox(iz)>=mfoth) cycle

    if (iz/=nz) then 
        mfox(iz) = max(0d0, &
            & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
            & )
    else 
        mfox(iz) = max(0d0, &
            & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
            & )
    endif 

    if (mabx(iz)>=mabth) cycle

    if (iz/=nz) then 
        mabx(iz) = max(0d0, &
            & mab(iz) +dt*(w*(mab(iz+1)-mab(iz))/dz(iz) + mabsupp(iz)) &
            & )
    else 
        mabx(iz) = max(0d0, &
            & mab(iz) + dt*(w*(mabi-mab(iz))/dz(iz)+ mabsupp(iz)) &
            & )
    endif 

    if (manx(iz)>=manth) cycle

    if (iz/=nz) then 
        manx(iz) = max(0d0, &
            & man(iz) +dt*(w*(man(iz+1)-man(iz))/dz(iz) + mansupp(iz)) &
            & )
    else 
        manx(iz) = max(0d0, &
            & man(iz) + dt*(w*(mani-man(iz))/dz(iz)+ mansupp(iz)) &
            & )
    endif 

    if (mccx(iz)>=mccth) cycle

    if (iz/=nz) then 
        mccx(iz) = max(0d0, &
            & mcc(iz) +dt*(w*(mcc(iz+1)-mcc(iz))/dz(iz) + mccsupp(iz) &
            ! & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0) &
            ! & *merge(1d0,0d0,cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0 > 0d0) &
            & ) &
            & )
    else 
        mccx(iz) = max(0d0, &
            & mcc(iz) + dt*(w*(mcci-mcc(iz))/dz(iz)+ mccsupp(iz)&
            ! & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0) &
            ! & *merge(1d0,0d0,cax(iz)*k1*k2*kco2*pco2x(iz)/(prox(iz)**2d0)/kcceq - 1d0 > 0d0) &
            & ) &
            & )
    endif 

enddo

endsubroutine precalc_slds_wopy

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_pw_py( &
    & nz,dt,cth,c2th,v,c2,dz,dfe3,poro,sat,tora,koxa,po2x,c,dfe2,koxs,koxs2,hr,mvpy &! input
    & ,ms,po2,so4,so4th,ci,c2i,so4i,dso4,msx &! input 
    & ,c2x,cx,so4x  &! output
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::dt,cth,c2th,dz,dfe3,dfe2,mvpy,so4th,ci,c2i,so4i,dso4
real(kind=8),dimension(nz),intent(in)::v,c2,poro,sat,tora,koxa,po2x,c,koxs,koxs2,hr,ms,po2,so4,msx
real(kind=8),dimension(nz),intent(out)::c2x,cx,so4x

integer iz
real(kind=8) ctmp

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
            & +15d0*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3  &
            & +koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3  &
            & )  &
            & )
    else 
        cx(iz) = max(0.0d0, &
            & c(iz) + dt*(-v(iz)*(c(iz)-ci)/dz+dfe2*tora(iz)*(ctmp+ci-2d0*c(iz))/(dz**2d0) &
            & +15d0*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
            & +koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
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
            & +2d0*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
            & +2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
            & ) &
            & )
    else 
        so4x(iz) = max(0.0d0, &
            & so4(iz) +dt*(-v(iz)*(so4(iz)-so4i)/dz+dso4*tora(iz)*(ctmp+so4i-2d0*so4(iz))/(dz**2d0) &
            & +2d0*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*c(iz)**(-0.40d0)*1d-3 &
            & +2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*ms(iz)*po2(iz)**0.50d0*1d-3 &
            & ) &
            & )
    endif 
enddo

endsubroutine precalc_pw_py

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_pw_sil( &
    & nz,nath,mgth,cath,sith,dt,v,na,ca,mg,si,dz,dna,dsi,dmg,dca,tora,poro,sat,nai,mgi,cai,sii &! input 
    & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx &! input 
    & ,nax,six,cax,mgx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::nath,mgth,cath,sith,dt,dz,dna,dsi,dmg,dca,nai,mgi,cai,sii,mvab,mvan,mvfo,mvcc
real(kind=8),dimension(nz),intent(in)::v,na,ca,mg,si,tora,poro,sat,kab,kan,kcc,kfo,hr,mabx,manx,mfox,mccx
real(kind=8),dimension(nz),intent(out)::nax,six,cax,mgx

integer iz
real(kind=8) ctmp

do iz = 1, nz

    if (nax(iz)>=nath) cycle

    if (iz/=nz) ctmp = na(iz+1)
    if (iz==nz) ctmp = na(iz)
    if (iz/=1) then 
        nax(iz) = max(0.0d0, &
            & na(iz) +dt*(-v(iz)*(na(iz)-na(iz-1))/dz+dna*tora(iz)*(ctmp+na(iz-1)-2d0*na(iz))/(dz**2d0) &
            & +dna/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
            & *(na(iz)-na(iz-1))/(dz**2d0) &
            & +kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*1d-3 &
            & ) &
            & )
    else 
        nax(iz) = max(0.0d0, &
            & na(iz) + dt*(-v(iz)*(na(iz)-nai)/dz+dna*tora(iz)*(ctmp+nai-2d0*na(iz))/(dz**2d0) &
            & +kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*1d-3 &
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
            & +2d0*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*1d-3 &
            & ) &
            & )
    else 
        six(iz) = max(0.0d0, &
            & si(iz) + dt*(-v(iz)*(si(iz)-sii)/dz+dsi*tora(iz)*(ctmp+sii-2d0*si(iz))/(dz**2d0) &
            & +kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*1d-3 &
            & +2d0*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*1d-3 &
            & ) &
            & )
    endif 
        
end do
        
        
do iz = 1, nz
    if (cax(iz)>=cath) cycle

    if (iz/=nz) ctmp = ca(iz+1)
    if (iz==nz) ctmp = ca(iz)
    if (iz/=1) then 
        cax(iz) = max(0.0d0, &
            & ca(iz) +dt*(-v(iz)*(ca(iz)-ca(iz-1))/dz+dca*tora(iz)*(ctmp+ca(iz-1)-2d0*ca(iz))/(dz**2d0) &
            & +dca/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
            & *(ca(iz)-ca(iz-1))/(dz**2d0) &
            & +kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*1d-3 &
            & ) &
            & )
    else 
        cax(iz) = max(0.0d0, &
            & ca(iz) + dt*(-v(iz)*(ca(iz)-cai)/dz+dca*tora(iz)*(ctmp+cai-2d0*ca(iz))/(dz**2d0) &
            & +kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*1d-3 &
            & ) &
            & )
    endif 
        
end do

endsubroutine precalc_pw_sil

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_pw_sil_v2( &
    & nz,nath,mgth,cath,sith,dt,v,na,ca,mg,si,dz,dna,dsi,dmg,dca,tora,poro,sat,nai,mgi,cai,sii &! input 
    & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx &! input 
    & ,nax,six,cax,mgx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::nath,mgth,cath,sith,dt,dna,dsi,dmg,dca,nai,mgi,cai,sii,mvab,mvan,mvfo,mvcc
real(kind=8),dimension(nz),intent(in)::v,na,ca,mg,si,tora,poro,sat,kab,kan,kcc,kfo,hr,mabx,manx,mfox,mccx,dz
real(kind=8),dimension(nz),intent(out)::nax,six,cax,mgx

integer iz
real(kind=8) ctmp,edifi,ediftmp
real(kind=8),dimension(nz)::edif

edif = poro*sat*1d3*dna*tora
edifi = edif(1)

do iz = 1, nz

    if (nax(iz)>=nath) cycle

    ctmp = na(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) ctmp = nai
    if (iz==1) ediftmp = edifi
    
    nax(iz) = max(0.0d0, &
        & na(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
        & -poro(iz)*sat(iz)*1d3*v(iz)*(na(iz)-ctmp)/dz(iz) &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(na(min(nz,iz+1))-na(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(na(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & +kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz) &
        & ) &
        & )
        
enddo

edif = poro*sat*1d3*dmg*tora
edifi = edif(1)
        
do iz = 1, nz

    if (mgx(iz)>=mgth) cycle

    ctmp = mg(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) ctmp = mgi
    if (iz==1) ediftmp = edifi
    
    mgx(iz) = max(0.0d0, &
        & mg(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
        & -poro(iz)*sat(iz)*1d3*v(iz)*(mg(iz)-ctmp)/dz(iz) &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(mg(min(nz,iz+1))-mg(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(mg(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & +2d0*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz) &
        & ) &
        & )
enddo
        
edif = poro*sat*1d3*dsi*tora
edifi = edif(1)
        
do iz = 1, nz

    if (six(iz)>=sith) cycle

    ctmp = si(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) ctmp = sii
    if (iz==1) ediftmp = edifi
    
    six(iz) = max(0.0d0, &
        & si(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
        & -poro(iz)*sat(iz)*1d3*v(iz)*(si(iz)-ctmp)/dz(iz) &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(si(min(nz,iz+1))-si(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(si(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & +kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz) &
        & +2d0*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz) &
        & ) &
        & )
end do
        
edif = poro*sat*1d3*dca*tora
edifi = edif(1)
        
do iz = 1, nz


    if (cax(iz)>=cath) cycle

    ctmp = ca(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) ctmp = cai
    if (iz==1) ediftmp = edifi
    
    cax(iz) = max(0.0d0, &
        & ca(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
        & -poro(iz)*sat(iz)*1d3*v(iz)*(ca(iz)-ctmp)/dz(iz) &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(ca(min(nz,iz+1))-ca(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(ca(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & +kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz) &
        & ) &
        & )
        
end do

endsubroutine precalc_pw_sil_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine pyweath_1D( &
    & nz,nflx_py,mvpy,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
    & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daq,dgas,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
    & ,swoxa,swoxall,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,cx,c2x,po2x,msx,flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2 &! output
    ) 
    
implicit none 

integer,intent(in)::nz,nflx_py,nsp
real(kind=8),intent(in)::ci,c2i,po2i,msi,po2th,dz,w,msth,dfe2,dfe3,kho,dt2,cth,c2th,stoxa,daq,dgas &
    & ,swbr,mo2,stoxs,tol,zrxn,swoxa,swoxall,ucv,vmax,mvpy
real(kind=8),dimension(nz),intent(in)::c,c2,po2,ms,hr,poro,z,koxs2,koxs,sat,dporodta,dporodtg &
    & ,tora,torg,v,koxa
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(out)::cx,c2x,po2x,msx
real(kind=8),dimension(nflx_py,nz),intent(out)::flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2

integer::itflx,iadv,idif,irxn_pyo2,irxn_pyfe3,irxn_fe2o2,iresp,ires
data itflx,iadv,idif,irxn_pyo2,irxn_pyfe3,irxn_fe2o2,iresp,ires/1,2,3,4,5,6,7,8/

integer,intent(inout)::iter,it
real(kind=8),intent(inout)::error,dt

integer iz,row,nmx,ie,ie2

real(kind=8),parameter::infinity = huge(0d0)

real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV


nmx = nsp*nz

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0
    
    flx_py = 0d0
    flx_py_fe2 = 0d0
    flx_py_fe3 = 0d0
    flx_py_o2 = 0d0

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
                & + w/dz     &
                & + koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6 &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &
                & + merge(0.0d0, &
                & + koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                &  cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & ) &
                & * merge(1.0d0,msx(iz),msx(iz)<msth)

            amx(row,row+nsp) = (-w/dz)*merge(1.0d0,msx(iz+1),msx(iz)<msth)

            ymx(row) = ( &
                & (msx(iz)-ms(iz))/dt  &
                & -w*(msx(iz+1)-msx(iz))/dz &
                & + koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))  &
                & + merge(0.0d0, &
                & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & ) &
                & *merge(0.0d0,1d0,msx(iz)<msth)
            
            flx_py(iadv,iz) = (& 
                & -w*(msx(iz+1)-msx(iz))/dz &
                & )

        else if (iz==nz) then

            amx(row,row) = (1.0d0/dt  &
                & + w/dz  &
                & + koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))  &
                & + merge(0.0d0,  &
                & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))  &
                & )  &
                & *merge(1.0d0,msx(iz),msx(iz)<msth)

            ymx(row) = (   &
                & (msx(iz)-ms(iz))/dt  & 
                & -w*(msi-msx(iz))/dz &
                & + koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &         
                & + merge(0.0d0, &
                & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
                & ) &
                & *merge(0.0d0,1d0,msx(iz)<msth)
            
            flx_py(iadv,iz) = (& 
                & -w*(msi-msx(iz))/dz &
                & )
        end if 

        amx(row,row + 2 ) = ( & 
            & + merge(0.0d0, &
            & koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*0.50d0*(po2x(iz)**(-0.50d0)), &
            & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
            & ) &
            & *po2x(iz)&
            & *merge(0.0d0,1d0,msx(iz)<msth) 

        amx(row,row  + 1 ) = (  &
            & + merge(0.0d0,  &
            & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-1.40d0), &
            & cx(iz)<cth .or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
            & ) &
            & *cx(iz)&
            & *merge(0.0d0,1d0,msx(iz)<msth)

        amx(row,row  + 3 ) = ( &
            & + merge(0.0d0,  &
            & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0),  &
            & (c2x(iz)<c2th).or.(cx(iz)<cth).or.c2x(iz)<c2th.or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0))) &
            & ) &
            & *c2x(iz)  &
            & *merge(0.0d0,1d0,msx(iz)<msth)
            
        flx_py(itflx,iz) = (& 
            & (msx(iz)-ms(iz))/dt  & 
            & )
        flx_py(irxn_pyo2,iz) = (& 
            & + koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0)) &  
            & )
        flx_py(irxn_pyfe3,iz) = (& 
            & + merge(0.0d0, &
            & koxs2(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0))) &
            & )
        flx_py(ires,iz) = sum(flx_py(:,iz))


    end do  !================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Fe++   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do iz = 1, nz

        row = nsp*(iz-1)+2

        if (.not.((iz == 1).or.(iz==nz))) then

            amx(row,row) = (  &
                & 1.0d0/dt   &
                & +dporodta(iz)   &
                & +(-dfe2*tora(iz)*(-2d0)/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0))  &
                & + v(iz)/dz  &
                & + koxa(iz)*po2x(iz) &
                & + merge(0.0d0,  &
                & -(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)*  &
                & c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))  &
                & *1d-3  &
                &  )  &
                & *merge(1.0d0,cx(iz),cx(iz)<cth)

            amx(row,row-nsp) = (  &
                & +(-dfe2*tora(iz)*(1d0)/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0)) &
                & - v(iz)/dz &
                & ) &
                & *cx(iz-1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            amx(row,row+nsp) = ( &
                & +(-dfe2*tora(iz)*(1d0)/(dz**2d0)) &
                & ) &
                & *cx(iz+1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt &
                & +dporodta(iz) *cx(iz) &
                & -dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0) &
                & + v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + koxa(iz)*cx(iz)*po2x(iz) &
                & - koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)  &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3  &
                & - merge(0.0d0,  &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0),  &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & )  &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_py_fe2(idif,iz) = (  &
                & -dfe2*tora(iz)*(cx(iz+1)+cx(iz-1)-2d0*cx(iz))/(dz**2d0)  &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0) &
                & )  
            flx_py_fe2(iadv,iz) = (  &
                & + v(iz)*(cx(iz)-cx(iz-1))/dz &
                & ) 

        else if (iz == 1) then

            amx(row,row) = (1.0d0/dt +dporodta(iz)  &
                & - dfe2*tora(iz)*(-2d0)/(dz**2d0) &
                & + v(iz)/dz &
                & + koxa(iz)*po2x(iz) &
                & - merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)&
                & *c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,cx(iz),cx(iz)<c2th)

            amx(row,row+nsp) = ( &
                & -dfe2*tora(iz)*(1d0)/(dz**2d0) &
                & ) &
                & *cx(iz+1)  &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt  &
                & +dporodta(iz)*cx(iz)  &
                & -dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))/(dz**2d0) &
                & + v(iz)*(cx(iz)-ci)/dz &
                & + koxa(iz)*cx(iz)*po2x(iz) &
                & - koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & - merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_py_fe2(idif,iz) = (  &
                & -dfe2*tora(iz)*(cx(iz+1)+ci-2d0*cx(iz))/(dz**2d0) &
                & )  
            flx_py_fe2(iadv,iz) = (  &
                & + v(iz)*(cx(iz)-ci)/dz &
                & ) 

        else if (iz == nz) then

            amx(row,row) = (  &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & -dfe2*tora(iz)*(-1d0)/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0) &
                & + v(iz)/dz &
                & + koxa(iz)*po2x(iz) &
                & + merge(0.0d0, &
                & -(15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)&
                & *c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,cx(iz),cx(iz)<cth)

            amx(row,row-nsp) = ( &
                & -dfe2*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0) &
                & -v(iz)/dz &
                & ) &
                & *cx(iz-1) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            ymx(row) = (  &
                & (cx(iz)-c(iz))/dt &
                & +dporodta(iz) *cx(iz) &
                & -dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0) &
                & + v(iz)*(cx(iz)-cx(iz-1))/dz &
                & + koxa(iz)*cx(iz)*po2x(iz) &
                & - koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall) &
                & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
                & -merge(0.0d0, &
                & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

            flx_py_fe2(idif,iz) = (  &
                & -dfe2*tora(iz)*(cx(iz-1)-1d0*cx(iz))/(dz**2d0) &
                & -dfe2/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(cx(iz)-cx(iz-1))/(dz**2d0) &
                & )  
            flx_py_fe2(iadv,iz) = (  &
                & + v(iz)*(cx(iz)-cx(iz-1))/dz &
                & ) 

        end if 

        amx(row,row+1) = (  &
            & + merge(0.0d0,koxa(iz)*cx(iz),po2x(iz)<po2th) &
            & - merge(0.0d0, &
            & koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*0.50d0*(1d0-swoxall)*po2x(iz)**(-0.50d0), &
            & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0)))*1d-3  &
            & ) &
            & *po2x(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        amx(row,row+2) = ( &
            & - merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)* &
            & (0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
            & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *c2x(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        amx(row,row  - 1) = (     &
            & - koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or. isnan(po2x(iz)**(0.50d0)))*1d-3  &
            & - merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *msx(iz) &
            & *merge(0.0d0,1.0d0,cx(iz)<cth)   ! commented out (is this necessary?)

        flx_py_fe2(itflx,iz) = ( &
            & (cx(iz)-c(iz))/dt &
            & +dporodta(iz) *cx(iz) &
            & ) 
        flx_py_fe2(irxn_fe2o2,iz) = (  & 
            & + koxa(iz)*cx(iz)*po2x(iz) &
            & ) 
        flx_py_fe2(irxn_pyo2,iz) = (  &
            & - koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & ) 
        flx_py_fe2(irxn_pyfe3,iz) = ( &
            & -merge(0.0d0, &
            & (15.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(1d0-swoxall)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) 
        flx_py_fe2(ires,iz) = sum(flx_py_fe2(:,iz))

        if (isnan(ymx(row))) then 
            print*,'**** NAN found for Fe++ at iz = ',iz ,iter
        endif 

    end do  ! ==============================   

    do iz = 1, nz

        row = nsp*(iz-1)+4  !! fe+++

        if (.not.((iz == 1).or.(iz==nz))) then

            amx(row,row) = ( &
                & 1.0d0/dt & 
                & +dporodta(iz)  &
                & -dfe3*tora(iz)*(-2d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0) &
                & + v(iz)/dz  &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row - nsp) = ( &
                & -dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0) &
                & - v(iz)/dz &
                & ) &
                & *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)

            amx(row,row+nsp) = ( &
                & -dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & ) &
                & *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)


            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt &
                & +dporodta(iz) *c2x(iz) &
                & -dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0) &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz &
                & - koxa(iz)*cx(iz)*po2x(iz) &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)

            flx_py_fe3(idif,iz) = (  &
                & -dfe3*tora(iz)*(c2x(iz+1)+c2x(iz-1)-2d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0) &
                & )  
            flx_py_fe3(iadv,iz) = (  &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz &
                & ) 

        else if (iz == 1) then

            amx(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz)  &
                & -dfe3*tora(iz)*(-2d0)/(dz**2d0) &
                & + v(iz)/dz &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row+nsp) = ( &
                & -dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & ) &
                & *merge(0.0d0,c2x(iz+1),c2x(iz)<c2th)

            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt  &
                & +dporodta(iz) *c2x(iz) &
                & -dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))/(dz**2d0) &
                & + v(iz)*(c2x(iz)-c2i)/dz &
                & - koxa(iz)*cx(iz)*po2x(iz) &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)  ! commented out 

            flx_py_fe3(idif,iz) = (  &
                & -dfe3*tora(iz)*(c2x(iz+1)+c2i-2d0*c2x(iz))/(dz**2d0) &
                & )  
            flx_py_fe3(iadv,iz) = (  &
                & + v(iz)*(c2x(iz)-c2i)/dz &
                & ) 

        else if (iz ==nz) then

            amx(row,row) = ( &
                & 1.0d0/dt  &
                & +dporodta(iz) &
                & -dfe3*tora(iz)*(-1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0) &
                & + v(iz)/dz  &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*(0.93d0)*c2x(iz)**(0.93d0-1.0d0)*cx(iz)**(-0.40d0), &
                & (c2x(iz)<c2th).or.(cx(iz)<cth).or.isnan(c2x(iz)**(0.93d0)*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(1.0d0,c2x(iz),c2x(iz)<c2th)

            amx(row,row - nsp) = ( &
                & -dfe3*tora(iz)*(1d0)/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0) &
                & - v(iz)/dz &
                & ) &
                & *merge(0.0d0,c2x(iz-1),c2x(iz)<c2th)


            ymx(row) = ( &
                & (c2x(iz)-c2(iz))/dt   &
                & +dporodta(iz) *c2x(iz) &
                & -dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0) &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz &
                & - koxa(iz)*cx(iz)*po2x(iz) &
                & +merge(0.0d0, &
                & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
                & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
                & ) &
                & *merge(0.0d0,1.0d0,c2x(iz)<c2th)   ! commented out (is this necessary?)

            flx_py_fe3(idif,iz) = (  &
                & -dfe3*tora(iz)*(c2x(iz-1)-1d0*c2x(iz))/(dz**2d0) &
                & -dfe3/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(c2x(iz)-c2x(iz-1))/(dz**2d0) &
                & )  
            flx_py_fe3(iadv,iz) = (  &
                & + v(iz)*(c2x(iz)-c2x(iz-1))/dz &
                & ) 

        end if 

        amx(row,row-1) = ( &
            & - koxa(iz)*cx(iz) &
            & ) &
            & *merge(0.0d0,po2x(iz),c2x(iz)<c2th)

        amx(row,row-2) = ( &
            & - koxa(iz)*po2x(iz) &
            & + merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*(-0.40d0)*cx(iz)**(-0.40d0-1.0d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *merge(0.0d0,cx(iz),c2x(iz)<c2th)

        amx(row,row  - 3) = ( &
            & + merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) &
            & *merge(0.0d0,msx(iz),c2x(iz)<c2th)

        flx_py_fe3(itflx,iz) = ( &
            & (c2x(iz)-c2(iz))/dt   &
            & +dporodta(iz) *c2x(iz) &
            & ) 
        flx_py_fe3(irxn_fe2o2,iz) = (  & 
            & - koxa(iz)*cx(iz)*po2x(iz) &
            & ) 
        flx_py_fe3(irxn_pyfe3,iz) = ( &
            & +merge(0.0d0, &
            & (14.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.c2x(iz)<c2th.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & ) 
        flx_py_fe3(ires,iz) = sum(flx_py_fe3(:,iz))

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
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz)<po2th) &
                & +merge(0.d0 &
                & ,stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz)<po2th.or.isnan(msx(iz)*0.50d0*po2x(iz)**(-0.50d0))) &
                & +merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz &
                & + stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
                & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
                & +stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
                & +swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
                
            flx_py_o2(idif,iz) = (  &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & )  
            flx_py_o2(iadv,iz) = (  &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz &
                & ) 

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +1.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz) <po2th) &
                & +merge(0.0d0 &
                & ,stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
                & +merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & + poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
                & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
                & +stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
                &  +swbr*vmax & 
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
                
            flx_py_o2(idif,iz) = (  &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & )  
            flx_py_o2(iadv,iz) = (  &
                & + poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & ) 

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt & 
                & +dporodtg(iz)  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) *(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +merge(0.0d0,stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz),po2x(iz) <po2th) &
                & +merge(0.0d0, &
                & stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz)*0.50d0*po2x(iz)**(-0.50d0), &
                & po2x(iz) <po2th.or.isnan(po2x(iz)**(-0.50d0))) &
                & +merge(0.0d0, &
                & swbr*vmax*mo2/(po2x(iz)+mo2)**2.0d0, &
                & (po2x(iz) <po2th).or.(isnan(mo2/(po2x(iz)+mo2)**2.0d0))) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & - 1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)  &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & + poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
                & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
                & +stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
                & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
                & +swbr*vmax &
                & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
                
            flx_py_o2(idif,iz) = (  &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & )  
            flx_py_o2(iadv,iz) = (  &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & ) 

        end if 

        amx(row,row-1) = ( &
            & +poro(iz)*sat(iz)*1d3*koxa(iz)*po2x(iz)*stoxa               &
            & ) &
            & *merge(0.0d0,cx(iz),po2x(iz)<po2th)

        amx(row,row  - 2) = ( &
            & koxs(iz)*poro(iz)*hr(iz)*stoxs*mvpy*1d-6 &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
            & ) &
            & *merge(0.0d0,msx(iz),po2x(iz)<po2th)

        flx_py_o2(itflx,iz) = ( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
            & +dporodtg(iz) *po2x(iz) &
            & ) 
        flx_py_o2(irxn_fe2o2,iz) = (  & 
            & +stoxa*poro(iz)*sat(iz)*1d3*koxa(iz)*cx(iz) &
            & *merge(0d0,po2x(iz),po2x(iz)<po2th.or.isnan(po2x(iz))) &
            & ) 
        flx_py_o2(irxn_pyo2,iz) = ( &
            & +stoxs*koxs(iz)*poro(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**(0.50d0),po2x(iz) <po2th.or.isnan(po2x(iz)**(0.50d0))) &
            & ) 
        flx_py_o2(iresp,iz) = ( &
            & +swbr*vmax &
            & *merge(0d0,po2x(iz)/(po2x(iz)+mo2),(po2x(iz) <po2th).or.(isnan(po2x(iz)/(po2x(iz)+mo2)))) &
            & ) 
        flx_py_o2(ires,iz) = sum(flx_py_o2(:,iz))

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

subroutine oxygen_resp_1D( &
    & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
    & ,kho,tora,torg,daq,dgas,v,mo2,tol,runname,workdir,zrxn,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,po2x,flx_o2,resp &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::po2i,po2th,dz,kho,daq,dgas,mo2,tol,zrxn,ucv,vmax
real(kind=8),dimension(nz),intent(in)::po2,poro,z,sat,dporodtg,tora,torg,v
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(inout)::po2x
real(kind=8),dimension(nz),intent(out)::resp
real(kind=8),dimension(nflx,nz),intent(out)::flx_o2
integer,intent(inout)::iter
real(kind=8),intent(inout)::error,dt

integer iz,row,nmx,ie,ie2

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8),dimension(nz)::dresp_dpo2

integer::itflx,iadv,idif,iresp,irain,ires
data itflx,iadv,idif,iresp,irain,ires/1,2,3,4,5,6/

integer,parameter:: nsp = 1
real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV

nmx = nsp*nz

error = 1d4

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0
    
    flx_o2 = 0d0
    
    resp = vmax &
        & *merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
    
    dresp_dpo2 = merge(0.0d0, &
        & vmax*mo2/(po2x+mo2)**2.0d0, &
        & (po2x <po2th).or.(isnan(mo2/(po2x+mo2)**2.0d0)))
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iz = 1, nz

        row = 1 + nsp*(iz-1)

        if (iz == 1) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz &
                & +resp(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & )
            flx_o2(iadv,iz) = +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +resp(iz)  &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & ) 
            flx_o2(iadv,iz) = poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz 

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt & 
                & +dporodtg(iz)  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas)  &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +resp(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgas+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daq) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgas &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgas) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daq-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daq)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & ) 
            flx_o2(iadv,iz) = +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz

        end if 
        
        flx_o2(itflx,iz) = ( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
            & +dporodtg(iz) *po2x(iz) &
            & ) 
        flx_o2(iresp,iz) = resp(iz)
        flx_o2(ires,iz) = sum(flx_o2(:,iz))

    end do 

    ymx=-1.0d0*ymx
    ! pause

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
            po2x(iz) = po2x(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            po2x(iz) = po2x(iz)*0.50d0
        else   
            po2x(iz) = po2x(iz)*exp(ymx(row))
        endif

    end do

    error = maxval(exp(abs(ymx))) - 1.0d0
    ! pause
    if (isnan(error).or.info/=0 .or. any(isnan(po2x))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*,isnan(error), info/=0,any(isnan(po2x))
        print*,isnan(error),info/=0,any(isnan(po2x))
        stop
        po2x =po2
        iter = iter + 1
        cycle
    endif

    if (any(abs(po2x)>infinity)) then 
        error = 1d3 
        print*, '**** ERROR; values are ininity'
        do iz=1,nz
            if (abs(po2x(iz))>infinity) po2x(iz)=po2(iz)
        enddo 
        iter = iter + 1
        cycle
    endif 

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (po2x(iz) < 0.0d0) then
            po2x(iz) = po2x(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'o2 error',error,info
    ! pause
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

#ifdef display
print *
print *,'-=-=-=-=-=-= o2  -=-=-=-=-=-=-='
print *,'o2:', (po2x(iz),iz=1,nz,nz/5)
print *,'resp:', (resp(iz),iz=1,nz,nz/5)
print *,'dresp_dpo2:', (dresp_dpo2(iz),iz=1,nz,nz/5)
print *
! pause
#endif

endsubroutine oxygen_resp_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine oxygen_resp_1D_v2( &
    & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
    & ,kho,tora,torg,daq,dgas,v,mo2,tol,runname,workdir,zrxn,ucv,vmax,poroprev  &! inpput
    & ,dt,flgback &! inout
    & ,po2x,flx_o2,resp &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::po2i,po2th,kho,daq,dgas,mo2,tol,zrxn,ucv,vmax
real(kind=8),dimension(nz),intent(in)::po2,poro,z,sat,dporodtg,tora,torg,v,poroprev,dz
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(inout)::po2x
real(kind=8),dimension(nz),intent(out)::resp
real(kind=8),dimension(nflx,nz),intent(out)::flx_o2
real(kind=8),intent(inout)::dt
logical,intent(inout)::flgback

integer iz,row,nmx,ie,ie2,iter

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8),dimension(nz)::dresp_dpo2,kho2,kho2x,dkho2_do2,edif,dedif_do2,alpha,alphaprev,dalpha_do2
real(kind=8) edifi,kho2i,error

integer::itflx,iadv,idif,iresp,irain,ires
data itflx,iadv,idif,iresp,irain,ires/1,2,3,4,5,6/

integer,parameter:: nsp = 1
real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV

if (flgback) return 

nmx = nsp*nz

error = 1d4
iter = 0

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0
    
    flx_o2 = 0d0
    
    resp = vmax &
        & *merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
    
    dresp_dpo2 = merge(0.0d0, &
        & vmax*mo2/(po2x+mo2)**2.0d0, &
        & (po2x <po2th).or.(isnan(mo2/(po2x+mo2)**2.0d0)))
        
    kho2 = kho ! previous value; should not change through iterations 
    kho2x = kho
    kho2i = kho
    
    dkho2_do2 = 0d0
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgas+poro*sat*kho2x*1d3*tora*daq
    edifi = edif(1)
    edifi = ucv*1d3*dgas 
    
    dedif_do2 = poro*sat*dkho2_do2*1d3*tora*daq
    
    alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho2x*1d3
    alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*kho2*1d3
    
    dalpha_do2 = poro*sat*dkho2_do2*1d3
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iz = 1, nz

        row = nsp*(iz-1) + 1

        if (iz == 1) then

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_do2(iz)*po2x(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(dedif_do2(iz))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_do2(iz))*(po2x(iz)-po2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkho2_do2(iz)*po2x(iz))/dz(iz)*dt &
                & +dresp_dpo2(iz)*dt &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & + 0.5d0*(dedif_do2(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz) *dt&
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (alpha(iz)*po2x(iz)-alphaprev(iz)*po2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(po2x(iz)-po2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2i*po2i)/dz(iz)*dt &
                & +resp(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(po2x(iz)-po2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & )
            flx_o2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2i*po2i)/dz(iz) &
                & )

        else if (iz == nz) then

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_do2(iz)*po2x(iz)) &
                & -( - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_do2(iz))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkho2_do2(iz)*po2x(iz))/dz(iz)*dt &
                & +dresp_dpo2(iz)*dt &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_do2(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-kho2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkho2_do2(iz-1)*po2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (alpha(iz)*po2x(iz)-alphaprev(iz)*po2(iz)) &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) *dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2x(iz-1)*po2x(iz-1))/dz(iz)*dt &
                & +resp(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & ) 
            flx_o2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2x(iz-1)*po2x(iz-1))/dz(iz) &
                & )

        else

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_do2(iz)*po2x(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(dedif_do2(iz))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1))))  &
                & - 0.5d0*(dedif_do2(iz))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkho2_do2(iz)*po2x(iz))/dz(iz)*dt &
                & +dresp_dpo2(iz)*dt &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & + 0.5d0*(dedif_do2(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt  &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_do2(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-kho2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkho2_do2(iz-1)*po2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (alpha(iz)*po2x(iz)-alphaprev(iz)*po2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2x(iz-1)*po2x(iz-1))/dz(iz)*dt &
                & +resp(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(po2x(iz+1)-po2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(po2x(iz)-po2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & ) 
            flx_o2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(kho2x(iz)*po2x(iz)-kho2x(iz-1)*po2x(iz-1))/dz(iz) &
                & )

        end if 
        
        flx_o2(itflx,iz) = ( &
            & (alpha(iz)*po2x(iz)-alphaprev(iz)*po2(iz))/dt &
            & ) 
        flx_o2(iresp,iz) = resp(iz)
        flx_o2(ires,iz) = sum(flx_o2(:,iz))
        
        if (any(isnan(flx_o2(:,iz)))) then
            print *,flx_o2(:,iz)
        endif 

    end do 

    ymx=-1.0d0*ymx
    ! pause

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
            po2x(iz) = po2x(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            po2x(iz) = po2x(iz)*0.50d0
        else   
            po2x(iz) = po2x(iz)*exp(ymx(row))
        endif

    end do

    error = maxval(exp(abs(ymx))) - 1.0d0
    ! pause
    if (isnan(error).or.info/=0 .or. any(isnan(po2x))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*,isnan(error), info/=0,any(isnan(po2x))
        print*,isnan(error),info/=0,any(isnan(po2x))
        stop
        po2x =po2
        iter = iter + 1
        cycle
    endif

    if (any(abs(po2x)>infinity)) then 
        error = 1d3 
        print*, '**** ERROR; values are ininity'
        do iz=1,nz
            if (abs(po2x(iz))>infinity) po2x(iz)=po2(iz)
        enddo 
        iter = iter + 1
        cycle
    endif 

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (po2x(iz) < 0.0d0) then
            po2x(iz) = po2x(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'o2 error',error,info
    ! pause
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
        exit
    end if

end do

#ifdef display
print *
print *,'-=-=-=-=-=-= o2  -=-=-=-=-=-=-='
print *,'o2:', (po2x(iz),iz=1,nz,nz/5)
print *,'resp:', (resp(iz),iz=1,nz,nz/5)
print *,'dresp_dpo2:', (dresp_dpo2(iz),iz=1,nz,nz/5)
print *
! pause
#endif

endsubroutine oxygen_resp_1D_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine CO2_1D( &
    & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
    & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,preccc  &! inpput
    & ,iter,error,dt &! inout
    & ,pco2x,flx_co2,khco2 &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::pco2i,pco2th,dz,kco2,k1,k2,daqc,dgasc,tol,zrxn,ucv,khco2i
real(kind=8),dimension(nz),intent(in)::pco2,poro,z,sat,dporodtgc,tora,torg,v,resp,prox,preccc
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(inout)::pco2x
real(kind=8),dimension(nz),intent(out)::khco2
real(kind=8),dimension(nflx,nz),intent(out)::flx_co2
integer,intent(inout)::iter
real(kind=8),intent(inout)::error,dt

integer iz,row,nmx,ie,ie2

real(kind=8),parameter::infinity = huge(0d0)

integer::itflx,iadv,idif,iresp,irain,ires
data itflx,iadv,idif,iresp,irain,ires/1,2,3,4,5,6/

integer,parameter:: nsp = 1
real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV


nmx = nsp*nz

error = 1d4

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0
    
    flx_co2 = 0d0
    
    khco2 = kco2*(1d0+k1/prox + k1*k2/prox/prox)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iz = 1, nz

        row = 1 + nsp*(iz-1)

        if (iz == 1) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) &
                & +dporodtgc(iz)*dt  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2i)*1d3/dz*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3)*(pco2x(iz)-pco2(iz)) &
                & +dporodtgc(iz) *pco2x(iz)*dt &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz+1)+pco2i-2.0d0*pco2x(iz))/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2i)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2i)*1d3*pco2x(iz)/dz*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz+1)+pco2i-2.0d0*pco2x(iz))/(dz**2.0d0) &
                & )
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2i)/dz &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2i)*1d3*pco2x(iz)/dz &
                & )

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) &
                & +dporodtgc(iz)*dt  &
                & +(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(1.0d0)/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3/dz*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc &
                & -poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc))*(-1.0d0)/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(-1.0d0)/dz*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3)*(pco2x(iz)-pco2(iz)) &
                & +dporodtgc(iz) *pco2x(iz)*dt &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz-1)-1.0d0*pco2x(iz))/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(pco2x(iz)-pco2x(iz-1))/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2x(iz-1))/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3*pco2x(iz)/dz*dt &
                & -resp(iz)*dt  &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz-1)-1.0d0*pco2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(pco2x(iz)-pco2x(iz-1))/(dz**2.0d0) &
                & ) 
            flx_co2(iadv,iz) = ( &
                & poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2x(iz-1))/dz &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3*pco2x(iz)/dz &
                & )

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3) & 
                & +dporodtgc(iz)*dt  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(1.0d0)/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3/dz*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc &
                & +poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc)/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc)  &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc &
                & -poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc))*(-1.0d0)/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(-1.0d0)/dz*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3)*(pco2x(iz)-pco2(iz))  &
                & +dporodtgc(iz) *pco2x(iz)*dt &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz+1)+pco2x(iz-1)-2.0d0*pco2x(iz))/(dz**2.0d0)*dt &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(pco2x(iz)-pco2x(iz-1))/(dz**2.0d0)*dt &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2x(iz-1))/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3*pco2x(iz)/dz*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgasc+poro(iz)*sat(iz)*khco2(iz)*1d3*tora(iz)*daqc) &
                & *(pco2x(iz+1)+pco2x(iz-1)-2.0d0*pco2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgasc &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgasc) &
                & +(poro(iz)*sat(iz)*khco2(iz)*tora(iz)*daqc-poro(iz-1)*sat(iz-1)*khco2(iz-1)*tora(iz-1)*daqc)) &
                & *(pco2x(iz)-pco2x(iz-1))/(dz**2.0d0) &
                & ) 
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*khco2(iz)*1d3*(pco2x(iz)-pco2x(iz-1))/dz &
                & +poro(iz)*sat(iz)*v(iz)*(khco2(iz)-khco2(iz-1))*1d3*pco2x(iz)/dz &
                & )

        end if 
        
        flx_co2(itflx,iz) = ( &
            & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*khco2(iz)*1d3)*(pco2x(iz)-pco2(iz))/dt  &
            & +dporodtgc(iz) *pco2x(iz) &
            & ) 
        flx_co2(iresp,iz) = -resp(iz)
        flx_co2(irain,iz) = -preccc(iz)
        flx_co2(ires,iz) = sum(flx_co2(:,iz))
        
        if (any(isnan(flx_co2(:,iz)))) then
            print *,flx_co2(:,iz)
        endif 

    end do 

    ymx=-1.0d0*ymx
    ! pause

    if (any(isnan(amx)).or.any(isnan(ymx)).or.any(abs(amx)>infinity).or.any(abs(ymx)>infinity)) then 
        print*,'**** error in mtx co2 ***********'

        if (any(isnan(ymx))) then 
            open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_preCO2.txt',  &
                & status='replace', action = 'write')
            do ie = 1,nmx
                write(25,*) ymx(ie)
                if (isnan(ymx(ie))) then 
                    print*,'NAN is here...',ie,mod(ie,300)
                endif
            enddo
            close(25)
        endif


        if (any(isnan(amx))) then 
            open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx_preCO2.txt',  &
                & status='replace', action = 'write')
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
        print*,'error in soultion',info
        open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_aftCO2.txt',  &
        & status='replace', action = 'write')
        open(unit=26,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx_aftCO2.txt',  &
            & status='replace', action = 'write')
        do ie = 1,nmx
            write(25,*) ymx(ie)
            write(26,*) (amx(ie,ie2),ie2=1,nmx)
            if (isnan(ymx(ie))) then 
                print*,'NAN is here...',ie,mod(ie,300)
            endif
        enddo
        close(25)
        close(26)
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
            pco2x(iz) = pco2x(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            pco2x(iz) = pco2x(iz)*0.50d0
        else   
            pco2x(iz) = pco2x(iz)*exp(ymx(row))
        endif

    end do

    error = maxval(exp(abs(ymx))) - 1.0d0
    ! pause
    if (isnan(error).or.info/=0 .or. any(isnan(pco2x))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*,isnan(error), info/=0,any(isnan(pco2x))
        print*,isnan(error),info/=0,any(isnan(pco2x))
        stop
        pco2x =pco2
        iter = iter + 1
        cycle
    endif

    if (any(abs(pco2x)>infinity)) then 
        error = 1d3 
        print*, '**** ERROR; values are ininity'
        do iz=1,nz
            if (abs(pco2x(iz))>infinity) pco2x(iz)=pco2(iz)
        enddo 
        iter = iter + 1
        cycle
    endif 

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (pco2x(iz) < 0.0d0) then
            pco2x(iz) = pco2x(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'co2 error',error,info
    ! pause
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

#ifdef display
print *
print *,'-=-=-=-=-=-= co2  -=-=-=-=-=-=-='
print *,'co2:', (pco2x(iz),iz=1,nz,nz/5)
! pause
#endif

endsubroutine CO2_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine CO2_1D_v2( &
    & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
    & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,preccc,poroprev,pro  &! inpput
    & ,dt,flgback &! inout
    & ,pco2x,flx_co2 &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::pco2i,pco2th,kco2,k1,k2,daqc,dgasc,tol,zrxn,ucv,khco2i
real(kind=8),dimension(nz),intent(in)::pco2,poro,z,sat,dporodtgc,tora,torg,v,resp,prox,preccc,poroprev,pro,dz
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(inout)::pco2x
real(kind=8),dimension(nflx,nz),intent(out)::flx_co2
logical,intent(inout)::flgback
real(kind=8),intent(inout)::dt

integer iz,row,nmx,ie,ie2,iter

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8) edifi,error
real(kind=8),dimension(nz)::khco2,khco2x,dpreccc_dco2,dkhco2_dco2,edif,dedif_dco2,alpha,alphaprev,dalpha_dco2


integer::itflx,iadv,idif,iresp,irain,ires
data itflx,iadv,idif,iresp,irain,ires/1,2,3,4,5,6/

integer,parameter:: nsp = 1
real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV


if (flgback) return

nmx = nsp*nz

error = 1d4

do while ((.not.isnan(error)).and.(error > tol))

    amx=0.0d0
    ymx=0.0d0
    
    flx_co2 = 0d0
        
    khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
    khco2x = kco2*(1d0+k1/prox + k1*k2/prox/prox)
    
    dpreccc_dco2=0d0
    dkhco2_dco2 = 0d0
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc+poro*sat*khco2x*1d3*tora*daqc
    edifi = edif(1)
    edifi = ucv*1d3*dgasc 
    
    dedif_dco2 = poro*sat*dkhco2_dco2*1d3*tora*daqc
    
    alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2x*1d3
    alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khco2*1d3
    
    dalpha_dco2 = poro*sat*dkhco2_dco2*1d3
    
    if (any(isnan(alpha)).or.any(isnan(alphaprev))) then 
        print *, alpha 
        print *, alphaprev
        print *, poroprev
        print *, khco2
        print *, pro
        stop
    endif 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do iz = 1, nz

        row = nsp*(iz-1) + 1

        if (iz == 1) then

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)*pco2x(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(dedif_dco2(iz))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz(iz)*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row+nsp) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & + 0.5d0*(dedif_dco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt&
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(pco2x(iz)-pco2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2i*pco2i)/dz(iz)*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edifi)*(pco2x(iz)-pco2i)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  &
                & )
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2i*pco2i)/dz(iz) &
                & )

        else if (iz == nz) then

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)*pco2x(iz)) &
                & -( - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz(iz)*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row-nsp) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dco2(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & ) 
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz(iz) &
                & )

        else

            amx(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)*pco2x(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(dedif_dco2(iz))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1))))  &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz(iz)*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx(row,row+nsp) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & + 0.5d0*(dedif_dco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt  &
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            amx(row,row-nsp) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dco2(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            flx_co2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & ) 
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz(iz) &
                & )

        end if 
        
        flx_co2(itflx,iz) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz))/dt &
            & ) 
        flx_co2(iresp,iz) = -resp(iz)
        flx_co2(irain,iz) = -preccc(iz)
        flx_co2(ires,iz) = sum(flx_co2(:,iz))
        
        if (any(isnan(flx_co2(:,iz)))) then
            print *,flx_co2(:,iz)
        endif 

    end do 

    ymx=-1.0d0*ymx
    ! pause

    if (any(isnan(amx)).or.any(isnan(ymx)).or.any(abs(amx)>infinity).or.any(abs(ymx)>infinity)) then 
        print*,'**** error in mtx co2 ***********'

        open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_preCO2.txt',  &
            & status='replace', action = 'write')
        do ie = 1,nmx
            write(25,*) ymx(ie)
            if (isnan(ymx(ie))) then 
                print*,'NAN is here...',ie,mod(ie,300)
            endif
        enddo
        close(25)


        open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx_preCO2.txt',  &
            & status='replace', action = 'write')
        do ie = 1,nmx
            write(25,*) (amx(ie,ie2),ie2=1,nmx)
            do ie2 = 1,nmx
                if (isnan(amx(ie,ie2))) then 
                    print*,'NAN is here...',ie,mod(ie,nz),ie2,mod(ie2,nz)
                endif
            enddo
        enddo
        close(25)

        stop
    endif

    call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,INFO) 

    if (any(isnan(ymx))) then
        print*,'error in soultion',info
        open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_aftCO2.txt',  &
        & status='replace', action = 'write')
        open(unit=26,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx_aftCO2.txt',  &
            & status='replace', action = 'write')
        do ie = 1,nmx
            write(25,*) ymx(ie)
            write(26,*) (amx(ie,ie2),ie2=1,nmx)
            if (isnan(ymx(ie))) then 
                print*,'NAN is here...',ie,mod(ie,300)
            endif
        enddo
        close(25)
        close(26)
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
            pco2x(iz) = pco2x(iz)*1.5d0
        else if (ymx(row) < -10d0) then 
            pco2x(iz) = pco2x(iz)*0.50d0
        else   
            pco2x(iz) = pco2x(iz)*exp(ymx(row))
        endif

    end do

    error = maxval(exp(abs(ymx))) - 1.0d0
    ! pause
    if (isnan(error).or.info/=0 .or. any(isnan(pco2x))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*,isnan(error), info/=0,any(isnan(pco2x))
        print*,isnan(error),info/=0,any(isnan(pco2x))
        stop
        pco2x =pco2
        iter = iter + 1
        cycle
    endif

    if (any(abs(pco2x)>infinity)) then 
        error = 1d3 
        print*, '**** ERROR; values are ininity'
        do iz=1,nz
            if (abs(pco2x(iz))>infinity) pco2x(iz)=pco2(iz)
        enddo 
        iter = iter + 1
        cycle
    endif 

    do iz = 1, nz
        row = 1 + nsp*(iz-1)

        if (pco2x(iz) < 0.0d0) then
            pco2x(iz) = pco2x(iz)/exp(ymx(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'co2 error',error,info,iter
    ! pause
#endif      
    iter = iter + 1

    if (iter > 300) then
        dt = dt/1.01d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
        exit 
    end if

end do

#ifdef display
print *
print *,'-=-=-=-=-=-= co2  -=-=-=-=-=-=-='
print *,'co2:', (pco2x(iz),iz=1,nz,nz/5)
! pause
#endif

endsubroutine CO2_1D_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine CO2_1D_v3( &
    & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
    & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,preccc,poroprev,pro  &! inpput
    & ,dt &! inout
    & ,pco2x,flx_co2 &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::pco2i,pco2th,dz,kco2,k1,k2,daqc,dgasc,tol,zrxn,ucv,khco2i
real(kind=8),dimension(nz),intent(in)::pco2,poro,z,sat,dporodtgc,tora,torg,v,resp,prox,preccc,poroprev,pro
character(256),intent(in)::runname,workdir
real(kind=8),dimension(nz),intent(inout)::pco2x
real(kind=8),dimension(nflx,nz),intent(out)::flx_co2
real(kind=8),intent(inout)::dt

integer iz,row,nmx,ie,ie2

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8) edifi
real(kind=8),dimension(nz)::khco2,khco2x,dpreccc_dco2,dkhco2_dco2,edif,dedif_dco2,alpha,alphaprev,dalpha_dco2


integer::itflx,iadv,idif,iresp,irain,ires
data itflx,iadv,idif,iresp,irain,ires/1,2,3,4,5,6/

integer,parameter:: nsp = 1
real(kind=8) amx(nsp*nz,nsp*nz),ymx(nsp*nz)
integer ipiv(nsp*nz)
integer info

external DGESV


nmx = nsp*nz


amx=0.0d0
ymx=0.0d0

flx_co2 = 0d0
    
khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
khco2x = kco2*(1d0+k1/prox + k1*k2/prox/prox)

dpreccc_dco2=0d0
dkhco2_dco2 = 0d0

edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc+poro*sat*khco2x*1d3*tora*daqc
edifi = edif(1)
edifi = ucv*1d3*dgasc 

dedif_dco2 = poro*sat*dkhco2_dco2*1d3*tora*daqc

alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2x*1d3
alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khco2*1d3

dalpha_dco2 = poro*sat*dkhco2_dco2*1d3

if (any(isnan(alpha)).or.any(isnan(alphaprev))) then 
    print *, alpha 
    print *, alphaprev
    print *, poroprev
    print *, khco2
    print *, pro
    stop
endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do iz = 1, nz

    row = nsp*(iz-1) + 1

    if (iz == 1) then

        amx(row,row) = ( &
            & (alpha(iz)) &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/dz &
            & - 0.5d0*(edif(iz)+edifi)*(1d0)/dz )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
            & ) 

        amx(row,row+nsp) = ( &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/dz )/dz *dt&
            & ) 

        ymx(row) = ( &
            & (alpha(iz)*0d0-alphaprev(iz)*pco2(iz)) &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(0d0-0d00)/dz &
            & - 0.5d0*(edif(iz)+edifi)*(0d0-pco2i)/dz )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*0d0-khco2i*pco2i)/dz*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & ) 
        
        flx_co2(idif,iz) = ( &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
            & - 0.5d0*(edif(iz)+edifi)*(pco2x(iz)-pco2i)/dz )/dz  &
            & )
        flx_co2(iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2i*pco2i)/dz &
            & )

    else if (iz == nz) then

        amx(row,row) = ( &
            & (alpha(iz)) &
            & -( - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/dz  )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
            & ) 

        amx(row,row-nsp) = ( &
            & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/dz )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz*dt &
            & ) 

        ymx(row) = ( &
            & (alpha(iz)*0d0-alphaprev(iz)*pco2(iz)) &
            & -( 0d0 &
            & - 0.5d0*(edif(iz)+edif(iz-1))*(0d0-0d0)/dz )/dz *dt &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*0d0-khco2x(iz-1)*0d0)/dz*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & ) 
        
        flx_co2(idif,iz) = ( &
            & -( 0d0 &
            & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz  &
            & ) 
        flx_co2(iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz &
            & )

    else

        amx(row,row) = ( &
            & (alpha(iz)) &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/dz &
            & - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/dz  )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
            & ) 

        amx(row,row+nsp) = ( &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/dz )/dz*dt  &
            & ) 

        amx(row,row-nsp) = ( &
            & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/dz )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz*dt &
            & ) 

        ymx(row) = ( &
            & (alpha(iz)*0d0-alphaprev(iz)*pco2(iz)) &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(0d0-0d0)/dz &
            & - 0.5d0*(edif(iz)+edif(iz-1))*(0d0-0d0)/dz )/dz*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*0d0-khco2x(iz-1)*0d0)/dz*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & ) 
        
        flx_co2(idif,iz) = ( &
            & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
            & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz  &
            & ) 
        flx_co2(iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz &
            & )

    end if 
    
    flx_co2(itflx,iz) = ( &
        & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz))/dt &
        & ) 
    flx_co2(iresp,iz) = -resp(iz)
    flx_co2(irain,iz) = -preccc(iz)
    flx_co2(ires,iz) = sum(flx_co2(:,iz))
    
    if (any(isnan(flx_co2(:,iz)))) then
        print *,flx_co2(:,iz)
    endif 

end do 

ymx=-1.0d0*ymx
! pause

if (any(isnan(amx)).or.any(isnan(ymx)).or.any(abs(amx)>infinity).or.any(abs(ymx)>infinity)) then 
    print*,'**** error in mtx co2 ***********'

    open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/ymx_preCO2.txt',  &
        & status='replace', action = 'write')
    do ie = 1,nmx
        write(25,*) ymx(ie)
        if (isnan(ymx(ie))) then 
            print*,'NAN is here...',ie,mod(ie,300)
        endif
    enddo
    close(25)


    open(unit=25,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/amx_preCO2.txt',  &
        & status='replace', action = 'write')
    do ie = 1,nmx
        write(25,*) (amx(ie,ie2),ie2=1,nmx)
        do ie2 = 1,nmx
            if (isnan(amx(ie,ie2))) then 
                print*,'NAN is here...',ie,mod(ie,nz),ie2,mod(ie2,nz)
            endif
        enddo
    enddo
    close(25)

    stop
endif

call DGESV(nmx,int(1),amx,nmx,IPIV,ymx,nmx,INFO) 


do iz = 1, nz
    row = iz

    if (isnan(ymx(row))) then 
        print *,'nan at', iz,z(iz),zrxn,'pco2'
        if (z(iz)<zrxn) then 
            ymx(row)=0d0
            pco2x(iz) = 0.1d0*pco2th
        endif
    endif

    if ((.not.isnan(ymx(row))).and.(ymx(row)>=0d0)) then 
        pco2x(iz) = ymx(row)
    else
        if ((.not.isnan(ymx(row))).and.(abs(ymx(row))<=tol)) then
            pco2x(iz) = 0d0
        else 
            print *,'pco2 nan? or negative; stop'
            print*,iz,ymx(row)
            stop
        endif 
        pco2x(iz) = 0d0
    endif
end do 


#ifdef display
print *
print *,'-=-=-=-=-=-= co2  -=-=-=-=-=-=-='
print *,'co2:', (pco2x(iz),iz=1,nz,nz/5)
! pause
#endif

endsubroutine CO2_1D_v3

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine pyweath_1D_SO4( &
    & nz,nflx_py,mvpy,c,c2,po2,ms,hr,po2th,poro,z,dz,koxs2,koxs,dso4,sat,dporodta  &! input
    & ,cth,tora,v,tol,zrxn,dt,cx,c2x,po2x,msx,so4,swoxa,O2_evolution,so4i,so4th &! input
    & ,so4x,flx_py_so4 &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx_py
logical,intent(in)::O2_evolution
real(kind=8),intent(in)::po2th,dz,dso4,cth,tol,zrxn,dt,swoxa,so4i,so4th,mvpy
real(kind=8),dimension(nz),intent(in)::c,c2,po2,ms,hr,poro,z,koxs2,koxs,sat,dporodta,tora,v,cx,c2x,po2x,msx,so4
real(kind=8),dimension(nz),intent(out)::so4x
real(kind=8),dimension(nflx_py,nz),intent(out)::flx_py_so4

integer iz,row,ie,ie2

real(kind=8)::swex = 0.0d0 ! switch for explicit
real(kind=8)::frex = 0.0d0 ! fraction of explicit

integer::itflx,iadv,idif,irxn_pyo2,irxn_pyfe3,irxn_fe2o2,iresp,ires
data itflx,iadv,idif,irxn_pyo2,irxn_pyfe3,irxn_fe2o2,iresp,ires/1,2,3,4,5,6,7,8/

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
            & -dso4*tora(iz)*(-2d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0) &
            & + v(iz)/dz &
            & )

        amx2(row,row-1) = ( &
            & -dso4*tora(iz)*(1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz) &
            & *(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0) &
            & - v(iz)/dz &
            & )

        amx2(row,row+1) = ( &
            & -dso4*tora(iz)*(1d0)/(dz**2d0) &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt  &
            & - 2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & - merge(0.0d0 &
            & ,(2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & )

    else if (iz == 1) then

        amx2(row,row) = ( &
            & 1.0d0/dt  &
            & +dporodta(iz)  &
            & -dso4*tora(iz)*(-2d0)/(dz**2d0) &
            & + v(iz)/dz &
            & )

        amx2(row,row+1) = ( &
            & -dso4*tora(iz)*(1d0)/(dz**2d0) &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt  &
            & - 2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & -merge(0.0d0 &
            & ,(2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
            & )

    else if (iz == nz) then

        amx2(row,row) = ( &
            & 1.0d0/dt  &
            & +dporodta(iz)  &
            & -dso4*tora(iz)*(-1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0) &
            & + v(iz)/dz &
            & )

        amx2(row,row-1) = ( &
            & -dso4*tora(iz)*(1d0)/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0) &
            & -v(iz)/dz &
            & )

        ymx2(row) = ( &
            & (-so4(iz))/dt  &
            & - 2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
            & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
            & -merge(0.0d0, &
            & (2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
            & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
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


! flx calc
flx_py_so4 = 0d0

do iz = 1, nz

    if (.not.((iz == 1).or.(iz==nz))) then
        
        flx_py_so4(idif,iz) = ( & 
            & -dso4*tora(iz)*(so4x(iz-1)+so4x(iz+1)-2d0*so4x(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(so4x(iz)-so4x(iz-1))/(dz**2d0) &
            & )
        flx_py_so4(iadv,iz) = ( & 
            & + (so4x(iz)-so4x(iz-1))*v(iz)/dz &
            & )

    else if (iz == 1) then
        
        flx_py_so4(idif,iz) = ( & 
            & -dso4*tora(iz)*(so4x(iz+1)-2d0*so4x(iz))/(dz**2d0) &
            & )
        flx_py_so4(iadv,iz) = ( & 
            & + (so4x(iz))*v(iz)/dz &
            & )

    else if (iz == nz) then
        
        flx_py_so4(idif,iz) = ( & 
            & -dso4*tora(iz)*(so4x(iz-1)-so4x(iz))/(dz**2d0) &
            & -dso4/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(so4x(iz)-so4x(iz-1))/(dz**2d0) &
            & )
        flx_py_so4(iadv,iz) = ( & 
            & + (so4x(iz)-so4x(iz-1))*v(iz)/dz &
            & )

    end if 
        
    flx_py_so4(itflx,iz) = ( & 
        & (so4x(iz)-so4(iz))/dt  &
        & +dporodta(iz)*so4x(iz)   &
        & )
    flx_py_so4(irxn_pyo2,iz) = ( & 
        & - 2d0*koxs(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz) &
        & *merge(0d0,po2x(iz)**0.50d0,po2x(iz) <po2th.or.isnan(po2x(iz)**0.50d0))*1d-3 &
        & )
    flx_py_so4(irxn_pyfe3,iz) = ( & 
        & - merge(0.0d0 &
        & ,(2.0d0)*koxs2(iz)/sat(iz)*hr(iz)*mvpy*1d-6*msx(iz)*c2x(iz)**0.93d0*cx(iz)**(-0.40d0), &
        & cx(iz)<cth.or.isnan(c2x(iz)**0.93d0*cx(iz)**(-0.40d0)))*1d-3 &
        & )
    flx_py_so4(ires,iz) = sum(flx_py_so4(:,iz))

enddo 

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
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::dz,w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,pco2x,man,kan,mansupp,mcc,kcc,mccsupp
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx 
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 8
integer iz,row,nmx,ie,ie2,isp,iflx
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc
real(kind=8)::k1_fo = 10d0**(-6.85d0), E1_fo = 51.7d0, n1_fo = 0.5d0, k2_fo = 10d0**(-12.41d0),E2_fo = 38d0 &
    & ,k3_fo = 10d0**(-21.2d0),E3_fo = 94.1d0,n3_fo = -0.82d0  &
    & ,k1_ab = 10d0**(-10.16d0), E1_ab = 65d0, n1_ab = 0.457d0, k2_ab = 10d0**(-12.56d0),E2_ab = 69.8d8 &
    & ,k3_ab = 10d0**(-15.60d0),E3_ab = 71d0, n3_ab = -0.572d0 &
    & ,k1_an = 10d0**(-3.5d0),E1_an = 16.6d0,n1_an = 1.411d0,k2_an = 10d0**(-9.12d0), E2_an = 17.8d0 

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8)::dconc = 1d-6
real(kind=8)::threshold = 10d0
real(kind=8)::disonly = 0d0 ! for cc [1---yes, 0---no]
! real(kind=8)::disonly = 1d0 ! for cc 

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! print *, mgx
! print *, six
! print *, mfosupp
! stop

if (any(isnan(tora)))then 
    print*,tora
endif 

error = 1d4
iter = 0

if (it ==0) then 
    ! prox(:) = 0.5d0* (  &
        ! & -1d0*(2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,2d0*(cx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 
else 

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
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
    flx_an = 0d0
    flx_ca = 0d0
    flx_cc = 0d0

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*cax(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*cax(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )

    ! dprodna(:) = 0.5d0* ( &
        ! &  -1d0  &
        ! & + 0.5d0/sqrt( &
        ! & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:))*2d0 &
        ! & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) &
        ! & )

    ! dprodmg(:) = 0.5d0* ( &
        ! &  -2d0  &
        ! & + 0.5d0/sqrt( &
        ! & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:))*2d0 &
        ! & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))*2d0 &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 
    call calc_pH( &
        & nz,nax+dconc+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,dprodna &! output
        & ) 
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx+dconc-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,dprodmg &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodca = dprodmg
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    omega_fo(:) = mgx(:)**2d0*six(:)/(prox(:)**4d0)/keqfo
    domega_fo_dmg(:) = 2d0*mgx(:)*six(:)/(prox(:)**4d0)/keqfo + mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodmg(:)/keqfo
    domega_fo_dsi(:) = mgx(:)**2d0/(prox(:)**4d0)/keqfo
    domega_fo_dna(:) = mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodna(:)/keqfo
    domega_fo_dca(:) = mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodca(:)/keqfo
    
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
    domega_ab_dmg(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodmg(:)
    domega_ab_dca(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodca(:)
    
    ! An + 2H+ + H2O = kaolinite + Ca2+ 
    omega_an(:) = cax(:)/(prox(:)**2d0)/keqan
    domega_an_dca(:) = 1d0/(prox(:)**2d0)/keqan + cax(:)*(-2d0)/(prox(:)**3d0)*dprodca(:)/keqan
    domega_an_dna(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodna(:)/keqan
    domega_an_dmg(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodmg(:)/keqan
    domega_an_dsi(:) = 0d0
    
    ! Cc = Ca2+ + CO32- 
    omega_cc = cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc
    domega_cc_dca = k1*k2*kco2*pco2x/(prox**2d0)/keqcc + cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodca/keqcc
    domega_cc_dmg = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodmg/keqcc
    domega_cc_dna = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodna/keqcc
    domega_cc_dsi = 0d0

    do iz = 1, nz  !================================
        
        do isp = 1, 4
        
            row = nsp3*(iz-1)+isp
            
            if (isp==1) then  ! Fo
                k_tmp = kfo(iz)
                mv_tmp = mvfo
                omega_tmp = omega_fo(iz)
                omega_tmp_th = omega_tmp
                m_tmp = mfox(iz)
                mth_tmp = mfoth 
                mi_tmp = mfoi
                mp_tmp = mfox(min(nz,iz+1))
                msupp_tmp = mfosupp(iz)
                mprev_tmp = mfo(iz)
            elseif (isp==2)then ! Ab
                k_tmp = kab(iz)
                mv_tmp = mvab
                omega_tmp = omega_ab(iz)
                omega_tmp_th = omega_tmp
                m_tmp = mabx(iz)
                mth_tmp = mabth 
                mi_tmp = mabi
                mp_tmp = mabx(min(nz,iz+1))
                msupp_tmp = mabsupp(iz)
                mprev_tmp = mab(iz)
            elseif (isp==3)then ! An
                k_tmp = kan(iz)
                mv_tmp = mvan
                omega_tmp = omega_an(iz)
                omega_tmp_th = omega_tmp
                m_tmp = manx(iz)
                mth_tmp = manth 
                mi_tmp = mani
                mp_tmp = manx(min(nz,iz+1))
                msupp_tmp = mansupp(iz)
                mprev_tmp = man(iz)
            elseif (isp==4)then ! Cc
                k_tmp = kcc(iz)
                mv_tmp = mvcc
                omega_tmp = omega_cc(iz)
                omega_tmp_th = omega_tmp*disonly
                ! omega_tmp_th = 0d0
                m_tmp = mccx(iz)
                mth_tmp = mccth 
                mi_tmp = mcci
                mp_tmp = mccx(min(nz,iz+1))
                msupp_tmp = mccsupp(iz)
                mprev_tmp = mcc(iz)
            endif 

            if (iz==nz) then

                amx3(row,row) = (1.0d0/dt  &
                    & + w/dz  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp)/dt &
                    & -w*(mi_tmp-m_tmp)/dz &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
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
                elseif (isp==3) then 
                    flx_an(iadv,iz) = (&
                        & -w*(mi_tmp-m_tmp)/dz &
                        & )
                elseif (isp==4) then 
                    flx_cc(iadv,iz) = (&
                        & -w*(mi_tmp-m_tmp)/dz &
                        & )
                endif 

            else

                amx3(row,row) = (1.0d0/dt     &
                    & + w/dz    &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & * merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                amx3(row,row+nsp3) = (-w/dz) *merge(1.0d0,mp_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp)/dt &
                    & -w*(mp_tmp-m_tmp)/dz  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
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
                elseif (isp==3) then 
                    flx_an(iadv,iz) = (&
                        & -w*(mp_tmp-m_tmp)/dz  &
                        & )
                elseif (isp==4) then 
                    flx_cc(iadv,iz) = (&
                        & -w*(mp_tmp-m_tmp)/dz  &
                        & )
                endif 

            end if 
            
            if (isp==1) then 
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_fo(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_fo(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                    
                flx_fo(irain,iz) = (&
                    & - msupp_tmp  &
                    & )
                flx_fo(ires,iz) = sum(flx_fo(:,iz))
                if (isnan(flx_fo(ires,iz))) then 
                    print *,'fo',iz,(flx_fo(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==2) then 
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_ab(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_ab(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                    
                flx_ab(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_ab(ires,iz) = sum(flx_ab(:,iz))
                if (isnan(flx_ab(ires,iz))) then 
                    print *,'ab',iz,(flx_ab(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==3) then 
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_an(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_an(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                    
                flx_an(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_an(ires,iz) = sum(flx_an(:,iz))
                if (isnan(flx_an(ires,iz))) then 
                    print *,'an',iz,(flx_an(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==4) then 
                amx3(row,row + 1 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_cc(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                    
                flx_cc(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                    
                flx_cc(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_cc(ires,iz) = sum(flx_cc(:,iz))
                if (isnan(flx_cc(ires,iz))) then 
                    print *,'cc',iz,(flx_cc(iflx,iz),iflx=1,nflx)
                endif 
            endif 
        enddo 
    end do  !================================

    do iz = 1, nz
        
        do isp = 1, 4

            row = nsp3*(iz-1)+4 + isp
            
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
                st_an = 0d0
                st_cc = 0d0
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
                st_an = 0d0
                st_cc = 0d0
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
                st_an = 0d0
                st_cc = 0d0
                rxn_tmp = st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
                drxndisp_tmp = st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3 
            elseif (isp==4) then  ! ca
                d_tmp = dca
                caq_tmp = cax(iz)
                caq_tmp_prev = ca(iz)
                caq_tmp_p = cax(min(nz,iz+1))
                caq_tmp_n = cax(max(1,iz-1))
                caqth_tmp = cath
                caqi_tmp = cai
                st_fo = 0d0
                st_ab = 0d0
                st_an = 1d0
                st_cc = 1d0
                rxn_tmp = st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3 &
                    & + st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(1d0-omega_cc(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)*1d-3 
                drxndisp_tmp = ( &
                    & st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3 &
                    & + st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)*1d-3 &
                    & )
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
                
                if (any(isnan((/d_tmp,tora(iz),caq_tmp_p,caq_tmp_n,caq_tmp,dz &
                        & ,poro(iz),sat(iz),poro(iz-1),sat(iz-1),tora(iz-1)/)))) then 
                    print*,'nan iz isp',iz,isp
                    print *,(/d_tmp,tora(iz),caq_tmp_p,caq_tmp_n,caq_tmp,dz &
                        & ,poro(iz),sat(iz),poro(iz-1),sat(iz-1),tora(iz-1)/)
                endif 
                
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
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
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
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
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
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
                        & +(-d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp/poro(iz)/sat(iz)*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 


            end if 
            
            amx3(row,row  - isp - 3) = (     & 
                & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*1d0*(1d0-omega_fo(iz)) &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                & ) &
                & *mfox(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -2) = (     & 
                & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*1d0*(1d0-omega_ab(iz)) &
                & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                & ) &
                & *mabx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -1) = (     & 
                & - st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*1d0*(1d0-omega_an(iz)) &
                & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3  &
                & ) &
                & *manx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp) = (     & 
                & - st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*1d0*(1d0-omega_cc(iz)) &
                & *merge(0d0,1d0,1d0 < 0d0)*1d-3  &
                & ) &
                & *mccx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            if (isp==1) then 
            
                amx3(row,row  + 1) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 3) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_mg(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_mg(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_mg(ires,iz) = sum(flx_mg(:,iz))
                if (isnan(flx_mg(ires,iz))) then 
                    print *,'mg',iz,(flx_mg(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==2) then 
            
                amx3(row,row  - 1) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_fo*kfo(iz)/sat(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)*1d-3  &
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_si(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_si(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_si(ires,iz) = sum(flx_si(:,iz))
                if (isnan(flx_si(ires,iz))) then 
                    print *,'si',iz,(flx_si(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==3) then 
            
                amx3(row,row  - 2) = (     & 
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_ab*kab(iz)/sat(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)*1d-3  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_na(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_na(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_na(ires,iz) = sum(flx_na(:,iz))
                if (isnan(flx_na(ires,iz))) then 
                    print *,'na',iz,(flx_na(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==4) then 
            
                amx3(row,row  - 3) = (     & 
                    & - st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3  &
                    & - st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dmg(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)*1d-3  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3  &
                    & - st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)*1d-3  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_an*kan(iz)/sat(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)*1d-3  &
                    & - st_cc*kcc(iz)/sat(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)*1d-3  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_ca(itflx,iz) = (&
                    & (caq_tmp-caq_tmp_prev)/dt  &
                    & +dporodta(iz) *caq_tmp &
                    & ) 
                flx_ca(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_ca(ires,iz) = sum(flx_ca(:,iz))
                if (isnan(flx_ca(ires,iz))) then 
                    print *,'ca',iz,(flx_ca(iflx,iz),iflx=1,nflx)
                endif 
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

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mfox(iz) = mfox(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
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

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mabx(iz) = mabx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mabx(iz) = mabx(iz)*0.50d0
        else   
            mabx(iz) = mabx(iz)*exp(ymx3(row))
        endif
        
        row = 3 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'An'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            manx(iz) = manx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            manx(iz) = manx(iz)*0.50d0
        else   
            manx(iz) = manx(iz)*exp(ymx3(row))
        endif
        
        row = 4 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Cc'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mccx(iz) = mccx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mccx(iz) = mccx(iz)*0.50d0
        else   
            mccx(iz) = mccx(iz)*exp(ymx3(row))
        endif
        
        row = 5 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'mg'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                mgx(iz) = 0.1d0*mgth
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mgx(iz) = mgx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mgx(iz) = mgx(iz)*0.50d0
        else
            mgx(iz) = mgx(iz)*exp(ymx3(row))
        endif
        
        row = 6 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'si'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                six(iz) = 0.1d0*sith
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            six(iz) = six(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            six(iz) = six(iz)*0.50d0
        else
            six(iz) = six(iz)*exp(ymx3(row))
        endif
        
        row = 7 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'na'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                nax(iz) = 0.1d0*nath
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            nax(iz) = nax(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            nax(iz) = nax(iz)*0.50d0
        else
            nax(iz) = nax(iz)*exp(ymx3(row))
        endif
        
        row = 8 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'ca'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                cax(iz) = 0.1d0*cath
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            cax(iz) = cax(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            cax(iz) = cax(iz)*0.50d0
        else
            cax(iz) = cax(iz)*exp(ymx3(row))
        endif
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)).or. any(isnan(cax)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx)).or. any(isnan(manx)).or. any(isnan(mccx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx)) &
            & ,any(isnan(cax)),any(isnan(manx)),any(isnan(mccx))
        stop
        mgx = mg
        six = si
        nax = na
        cax = ca
        mfox = mfo
        mabx = mab
        manx = man
        mccx = mcc
        prox = pro
        iter = iter + 1
        cycle
    endif

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 

    co2 = kco2*pco2x
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

        if (manx(iz) < 0.0d0) then
            manx(iz) = manx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 4 + nsp3*(iz-1)

        if (mccx(iz) < 0.0d0) then
            mccx(iz) = mccx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 5 + nsp3*(iz-1)

        if (mgx(iz) < 0.0d0) then
            mgx(iz) = mgx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 6 + nsp3*(iz-1)

        if (six(iz) < 0.0d0) then
            six(iz) = six(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 7 + nsp3*(iz-1)

        if (nax(iz) < 0.0d0) then
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 8 + nsp3*(iz-1)

        if (cax(iz) < 0.0d0) then
            cax(iz) = cax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'silicate_dis error',error,info,iter,dt
#endif      
    iter = iter + 1 

    if (iter > 300) then
        ! dt = dt/1.01d0
        dt = dt/10d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
            stop
        endif 
        flgback = .true.
        ! exit 
    end if
    
enddo

endsubroutine silicate_dis_1D

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine silicate_dis_1D_v2( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,pco2x,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx 
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 8
integer iz,row,nmx,ie,ie2,isp,iflx
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p
real(kind=8)::k1_fo = 10d0**(-6.85d0), E1_fo = 51.7d0, n1_fo = 0.5d0, k2_fo = 10d0**(-12.41d0),E2_fo = 38d0 &
    & ,k3_fo = 10d0**(-21.2d0),E3_fo = 94.1d0,n3_fo = -0.82d0  &
    & ,k1_ab = 10d0**(-10.16d0), E1_ab = 65d0, n1_ab = 0.457d0, k2_ab = 10d0**(-12.56d0),E2_ab = 69.8d8 &
    & ,k3_ab = 10d0**(-15.60d0),E3_ab = 71d0, n3_ab = -0.572d0 &
    & ,k1_an = 10d0**(-3.5d0),E1_an = 16.6d0,n1_an = 1.411d0,k2_an = 10d0**(-9.12d0), E2_an = 17.8d0 

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-6
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 100d0
real(kind=8)::disonly = 0d0 ! for cc [1---yes, 0---no]
! real(kind=8)::disonly = 1d0 ! for cc 

integer,parameter :: iter_max = 100

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! print *, mgx
! print *, six
! print *, mfosupp
! stop

if (any(isnan(tora)))then 
    print*,tora
endif 

error = 1d4
iter = 0

if (it ==0) then 
    ! prox(:) = 0.5d0* (  &
        ! & -1d0*(2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((2d0*cx(:)+2d0*ca(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,2d0*(cx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 
else 

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
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
    flx_an = 0d0
    flx_ca = 0d0
    flx_cc = 0d0

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*cax(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*cax(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )

    ! dprodna(:) = 0.5d0* ( &
        ! &  -1d0  &
        ! & + 0.5d0/sqrt( &
        ! & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:))*2d0 &
        ! & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:)) &
        ! & )

    ! dprodmg(:) = 0.5d0* ( &
        ! &  -2d0  &
        ! & + 0.5d0/sqrt( &
        ! & (nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:))*2d0 &
        ! & *(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))*2d0 &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 
    call calc_pH( &
        & nz,nax+dconc+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,dprodna &! output
        & ) 
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx+dconc-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,dprodmg &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodca = dprodmg
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    omega_fo(:) = mgx(:)**2d0*six(:)/(prox(:)**4d0)/keqfo
    domega_fo_dmg(:) = 2d0*mgx(:)*six(:)/(prox(:)**4d0)/keqfo + mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodmg(:)/keqfo
    domega_fo_dsi(:) = mgx(:)**2d0/(prox(:)**4d0)/keqfo
    domega_fo_dna(:) = mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodna(:)/keqfo
    domega_fo_dca(:) = mgx(:)**2d0*six(:)*(-4d0)/(prox(:)**5d0)*dprodca(:)/keqfo
    
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
    domega_ab_dmg(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodmg(:)
    domega_ab_dca(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodca(:)
    
    ! An + 2H+ + H2O = kaolinite + Ca2+ 
    omega_an(:) = cax(:)/(prox(:)**2d0)/keqan
    domega_an_dca(:) = 1d0/(prox(:)**2d0)/keqan + cax(:)*(-2d0)/(prox(:)**3d0)*dprodca(:)/keqan
    domega_an_dna(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodna(:)/keqan
    domega_an_dmg(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodmg(:)/keqan
    domega_an_dsi(:) = 0d0
    
    ! Cc = Ca2+ + CO32- 
    omega_cc = cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc
    domega_cc_dca = k1*k2*kco2*pco2x/(prox**2d0)/keqcc + cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodca/keqcc
    domega_cc_dmg = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodmg/keqcc
    domega_cc_dna = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodna/keqcc
    domega_cc_dsi = 0d0

    do iz = 1, nz  !================================
        
        do isp = 1, 4
        
            row = nsp3*(iz-1)+isp
            
            if (isp==1) then  ! Fo
                k_tmp = kfo(iz)
                mv_tmp = mvfo
                omega_tmp = omega_fo(iz)
                omega_tmp_th = omega_tmp
                m_tmp = mfox(iz)
                mth_tmp = mfoth 
                mi_tmp = mfoi
                mp_tmp = mfox(min(nz,iz+1))
                msupp_tmp = mfosupp(iz)
                mprev_tmp = mfo(iz)
            elseif (isp==2)then ! Ab
                k_tmp = kab(iz)
                mv_tmp = mvab
                omega_tmp = omega_ab(iz)
                omega_tmp_th = omega_tmp
                m_tmp = mabx(iz)
                mth_tmp = mabth 
                mi_tmp = mabi
                mp_tmp = mabx(min(nz,iz+1))
                msupp_tmp = mabsupp(iz)
                mprev_tmp = mab(iz)
            elseif (isp==3)then ! An
                k_tmp = kan(iz)
                mv_tmp = mvan
                omega_tmp = omega_an(iz)
                omega_tmp_th = omega_tmp
                m_tmp = manx(iz)
                mth_tmp = manth 
                mi_tmp = mani
                mp_tmp = manx(min(nz,iz+1))
                msupp_tmp = mansupp(iz)
                mprev_tmp = man(iz)
            elseif (isp==4)then ! Cc
                k_tmp = kcc(iz)
                mv_tmp = mvcc
                omega_tmp = omega_cc(iz)
                omega_tmp_th = omega_tmp*disonly
                ! omega_tmp_th = 0d0
                m_tmp = mccx(iz)
                mth_tmp = mccth 
                mi_tmp = mcci
                mp_tmp = mccx(min(nz,iz+1))
                msupp_tmp = mccsupp(iz)
                mprev_tmp = mcc(iz)
            endif 
            
            if (iz==nz) mp_tmp = mi_tmp

            amx3(row,row) = (1.0d0     &
                & + w/dz(iz)*dt    &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp)*dt &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & ) &
                & * merge(1.0d0,m_tmp,m_tmp<mth_tmp)

            ymx3(row) = ( &
                & (m_tmp-mprev_tmp) &
                & -w*(mp_tmp-m_tmp)/dz(iz)*dt  &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp)*dt &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & -msupp_tmp*dt  &
                & ) &
                & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                
            if (iz/=nz) amx3(row,row+nsp3) = (-w/dz(iz))*dt *merge(1.0d0,mp_tmp,m_tmp<mth_tmp)
            
            if (isp==1) then 
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_fo(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_fo(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_fo(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                flx_fo(irain,iz) = (&
                    & - msupp_tmp  &
                    & )
                flx_fo(ires,iz) = sum(flx_fo(:,iz))
                if (isnan(flx_fo(ires,iz))) then 
                    print *,'fo',iz,(flx_fo(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==2) then 
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_ab(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_ab(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_ab(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                flx_ab(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_ab(ires,iz) = sum(flx_ab(:,iz))
                if (isnan(flx_ab(ires,iz))) then 
                    print *,'ab',iz,(flx_ab(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==3) then 
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_an(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_an(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_an(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                flx_an(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_an(ires,iz) = sum(flx_an(:,iz))
                if (isnan(flx_an(ires,iz))) then 
                    print *,'an',iz,(flx_an(iflx,iz),iflx=1,nflx)
                endif 
                
            elseif (isp==4) then 
                amx3(row,row + 1 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_cc(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_cc(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_cc(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                flx_cc(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_cc(ires,iz) = sum(flx_cc(:,iz))
                if (isnan(flx_cc(ires,iz))) then 
                    print *,'cc',iz,(flx_cc(iflx,iz),iflx=1,nflx)
                endif 
            endif 
        enddo 
    end do  !================================

    do iz = 1, nz
        
        do isp = 1, 4

            row = nsp3*(iz-1)+4 + isp
            
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
                st_an = 0d0
                st_cc = 0d0
                rxn_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) 
                drxndisp_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dmg(iz)) & 
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)
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
                st_an = 0d0
                st_cc = 0d0
                rxn_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) 
                drxndisp_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) 
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
                st_an = 0d0
                st_cc = 0d0
                rxn_tmp = st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) 
                drxndisp_tmp = st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) 
            elseif (isp==4) then  ! ca
                d_tmp = dca
                caq_tmp = cax(iz)
                caq_tmp_prev = ca(iz)
                caq_tmp_p = cax(min(nz,iz+1))
                caq_tmp_n = cax(max(1,iz-1))
                caqth_tmp = cath
                caqi_tmp = cai
                st_fo = 0d0
                st_ab = 0d0
                st_an = 1d0
                st_cc = 1d0
                rxn_tmp = st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(1d0-omega_cc(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0) 
                drxndisp_tmp = ( &
                    & st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0) &
                    & )
            endif 
            
            if (iz==1) caq_tmp_n = caqi_tmp
                
            edif_tmp = 1d3*poro(iz)*sat(iz)*tora(iz)*d_tmp
            edif_tmp_p = 1d3*poro(min(iz+1,nz))*sat(min(iz+1,nz))*tora(min(iz+1,nz))*d_tmp
            edif_tmp_n = 1d3*poro(max(iz-1,1))*sat(max(iz-1,1))*tora(max(iz-1,1))*d_tmp

            amx3(row,row) = ( &
                & (poro(iz)*sat(iz)*1d3*1d0)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*merge(0d0,-1d0,iz==nz)/( 0.5d0*(dz(iz)+dz(min(nz,iz+1))) ) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(1d0)/( 0.5d0*(dz(iz)+dz(max(1,iz-1))) ))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(1d0)/dz(iz)*dt &
                & -drxndisp_tmp*dt &
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*dt &
                & - rxn_tmp*dt &
                & ) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

            if (iz/=1) then 
                amx3(row,row-nsp3) = ( &
                    & -(-0.5d0*(edif_tmp +edif_tmp_n)*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(-1d0)/dz(iz)*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            endif 
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            endif 
            
            amx3(row,row  - isp - 3) = (     & 
                & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*1d0*(1d0-omega_fo(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                & ) &
                & *mfox(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -2) = (     & 
                & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*1d0*(1d0-omega_ab(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                & ) &
                & *mabx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -1) = (     & 
                & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*1d0*(1d0-omega_an(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                & ) &
                & *manx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp) = (     & 
                & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*1d0*(1d0-omega_cc(iz))*dt &
                & *merge(0d0,1d0,1d0 < 0d0)  &
                & ) &
                & *mccx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            if (isp==1) then 
            
                amx3(row,row  + 1) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 3) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_mg(itflx,iz) = (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_mg(iadv,iz) = (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_mg(idif,iz) = (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_mg(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_mg(ires,iz) = sum(flx_mg(:,iz))
                if (isnan(flx_mg(ires,iz))) then 
                    print *,'mg',iz,(flx_mg(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==2) then 
            
                amx3(row,row  - 1) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_si(itflx,iz) = (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_si(iadv,iz) = (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_si(idif,iz) = (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_si(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_si(ires,iz) = sum(flx_si(:,iz))
                if (isnan(flx_si(ires,iz))) then 
                    print *,'si',iz,(flx_si(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==3) then 
            
                amx3(row,row  - 2) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_na(itflx,iz) = (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_na(iadv,iz) = (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_na(idif,iz) = (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_na(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_na(ires,iz) = sum(flx_na(:,iz))
                if (isnan(flx_na(ires,iz))) then 
                    print *,'na',iz,(flx_na(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==4) then 
            
                amx3(row,row  - 3) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_ca(itflx,iz) = (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_ca(iadv,iz) = (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_ca(idif,iz) = (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_ca(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_ca(ires,iz) = sum(flx_ca(:,iz))
                if (isnan(flx_ca(ires,iz))) then 
                    print *,'ca',iz,(flx_ca(iflx,iz),iflx=1,nflx)
                endif 
            endif 
            
            ! amx3(row,:) = amx3(row,:)/(poro(iz)*sat(iz)*1d3)
            ! ymx3(row) = ymx3(row)/(poro(iz)*sat(iz)*1d3)
        
        enddo 
        
    end do  ! ==============================
    
    ymx3=-1.0d0*ymx3

    if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
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
        
        open(unit=11,file='amx.txt',status = 'replace')
        open(unit=12,file='ymx.txt',status = 'replace')
        do ie = 1,nsp3*(nz)
            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
            write(12,*) ymx3(ie)
        enddo 
        close(11)
        close(12)       
    endif

    do iz = 1, nz
        row = 1 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Fo'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mfox(iz) = mfox(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
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

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mabx(iz) = mabx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mabx(iz) = mabx(iz)*0.50d0
        else   
            mabx(iz) = mabx(iz)*exp(ymx3(row))
        endif
        
        row = 3 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'An'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            manx(iz) = manx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            manx(iz) = manx(iz)*0.50d0
        else   
            manx(iz) = manx(iz)*exp(ymx3(row))
        endif
        
        row = 4 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'Cc'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mccx(iz) = mccx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mccx(iz) = mccx(iz)*0.50d0
        else   
            mccx(iz) = mccx(iz)*exp(ymx3(row))
        endif
        
        row = 5 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'mg'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                mgx(iz) = 0.1d0*mgth
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mgx(iz) = mgx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mgx(iz) = mgx(iz)*0.50d0
        else
            mgx(iz) = mgx(iz)*exp(ymx3(row))
        endif
        
        row = 6 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'si'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                six(iz) = 0.1d0*sith
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            six(iz) = six(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            six(iz) = six(iz)*0.50d0
        else
            six(iz) = six(iz)*exp(ymx3(row))
        endif
        
        row = 7 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'na'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                nax(iz) = 0.1d0*nath
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            nax(iz) = nax(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            nax(iz) = nax(iz)*0.50d0
        else
            nax(iz) = nax(iz)*exp(ymx3(row))
        endif
        
        row = 8 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'ca'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                cax(iz) = 0.1d0*cath
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            cax(iz) = cax(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            cax(iz) = cax(iz)*0.50d0
        else
            cax(iz) = cax(iz)*exp(ymx3(row))
        endif
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)).or. any(isnan(cax)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx)).or. any(isnan(manx)).or. any(isnan(mccx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx)) &
            & ,any(isnan(cax)),any(isnan(manx)),any(isnan(mccx))
        stop
        mgx = mg
        six = si
        nax = na
        cax = ca
        mfox = mfo
        mabx = mab
        manx = man
        mccx = mcc
        prox = pro
        iter = iter + 1
        cycle
    endif

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        & ,prox &! output
        & ) 

    co2 = kco2*pco2x
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

        if (manx(iz) < 0.0d0) then
            manx(iz) = manx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 4 + nsp3*(iz-1)

        if (mccx(iz) < 0.0d0) then
            mccx(iz) = mccx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 5 + nsp3*(iz-1)

        if (mgx(iz) < 0.0d0) then
            mgx(iz) = mgx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 6 + nsp3*(iz-1)

        if (six(iz) < 0.0d0) then
            six(iz) = six(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 7 + nsp3*(iz-1)

        if (nax(iz) < 0.0d0) then
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 8 + nsp3*(iz-1)

        if (cax(iz) < 0.0d0) then
            cax(iz) = cax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if

    end do 

#ifdef display      
    print *, 'silicate_dis error',error,info,iter,dt
#endif      
    iter = iter + 1 

    if (iter > iter_Max) then
        ! dt = dt/1.01d0
        dt = dt/10d0
        if (dt==0d0) then 
            print *, 'dt==0d0; stop'
        
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
        flgback = .true.
        
        ! mgx = mg
        ! six = si
        ! nax = na
        ! cax = ca
        ! mfox = mfo
        ! mabx = mab
        ! manx = man
        ! mccx = mcc
        ! prox = pro
        
        ! call precalc_slds_wopy( &
            ! & nz,dt,w,dz,mfoi,mabi,mani,mcci,mabth,manth,mfoth,mccth   &! input
            ! & ,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp &! input
            ! & ,mfox,mabx,manx,mccx &! output
            ! & )
        
        ! call precalc_pw_sil_v2( &
            ! & nz,nath,mgth,cath,sith,dt,v,na,ca,mg,si,dz,dna,dsi,dmg,dca,tora,poro,sat,nai,mgi,cai,sii &! input 
            ! & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx &! input 
            ! & ,nax,six,cax,mgx &! output
            ! & )
        exit 
    end if
    
#ifdef dispiter
    print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An -=-=-=-=-=-=-='
    print *, 'mg:', (mgx(iz),iz=1,nz, nz/5)
    print *, 'si:', (six(iz),iz=1,nz, nz/5)
    print *, 'na:', (nax(iz),iz=1,nz, nz/5)
    print *, 'ca:', (cax(iz),iz=1,nz, nz/5)
    print *, 'fo:', (mfox(iz),iz=1,nz, nz/5)
    print *, 'ab:', (mabx(iz),iz=1,nz, nz/5)
    print *, 'an:', (manx(iz),iz=1,nz, nz/5)
    print *, 'cc:', (mccx(iz),iz=1,nz, nz/5)
    print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
#endif     
enddo

endsubroutine silicate_dis_1D_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH( &
    & nz,netcat,pco2x,kw,kco2,k1,k2 &! input 
    & ,prox &! output
    & ) 
! solving charge balance:
! [H+] + ZX[Xz+] - ZY[YZ-] - [HCO3] - 2[CO32-] - [OH-] = 0
! [H+] + ZX[Xz+] - ZY[YZ-] - k1kco2pCO2/[H+] - 2k2k1kco2pCO2/[H+]^2 - kw/[H+] = 0
! [H+]^3 + (ZX[Xz+] - ZY[YZ-])[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
! NetCat is defined as (ZX[Xz+] - ZY[YZ-])
! [H+]^3 + NetCat[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
implicit none
integer,intent(in)::nz
real(kind=8),intent(in)::kw,kco2,k1,k2
real(kind=8),dimension(nz),intent(in)::netcat,pco2x
real(kind=8),dimension(nz),intent(out)::prox

real(kind=8),dimension(nz)::df,f
real(kind=8) error,tol
integer iter,iz

error = 1d4
tol = 1d-6

prox = 1d0 
iter = 0
! print*,'calc_pH'
do while (error > tol)
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x 
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw)
    df = df*prox
    prox = prox*exp( -f/df )
    error = maxval(abs(exp( -f/df )-1d0))
    ! print*, iter,error
    ! print*,  (-log10(prox(iz)),iz=1,nz,nz/5)
    ! print*,  (-log10(f(iz)),iz=1,nz,nz/5)
    ! print*,  (-log10(df(iz)),iz=1,nz,nz/5)
    ! pause
    ! stop
    iter = iter + 1
enddo 

if (any(isnan(prox))) then     
    print *, (-log10(prox(iz)),iz=1,nz,nz/5)
    stop
endif 

endsubroutine calc_pH

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

function k_arrhenius(kref,tempkref,tempk,eapp,rg)
implicit none
real(kind=8) k_arrhenius,kref,tempkref,tempk,eapp,rg
k_arrhenius = kref*exp(-eapp/rg*(1d0/tempk-1d0/tempkref))
endfunction k_arrhenius
