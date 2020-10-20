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

integer nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nz
character(5),dimension(:),allocatable::chraq,chrsld,chrgas,chrrxn_ext 
character(500) sim_name
real(kind=8) ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80

call get_variables_num( &
    & nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! output
    & )

! print *,nsp_sld,nsp_aq,nsp_gas,nrxn_ext

allocate(chraq(nsp_aq),chrsld(nsp_sld),chrgas(nsp_gas),chrrxn_ext(nrxn_ext))

    
call get_variables( &
    & nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! input
    & ,chraq,chrgas,chrsld,chrrxn_ext &! output
    & ) 
    
print *,chraq
print *,chrsld 
print *,chrgas 
print *,chrrxn_ext 

! pause

sim_name = 'chkchk'

call get_bsdvalues( &
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name &! output
    & )
    
call weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,ttot  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name &! input
    & )

endprogram weathering

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
subroutine weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,ttot  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name &! input
    & )

implicit none

!-----------------------------

real(kind=8),intent(in) :: ztot != 3.0d0 ! m
real(kind=8),intent(in) :: ttot  ! yr
! real(kind=8) dz
integer,intent(in) :: nz != 30 
real(kind=8) z(nz),dz(nz)
real(kind=8) ze(nz+1)
real(kind=8) :: ph = 5.0d0
real(kind=8) :: tc = 15.0d0 ! deg celsius
real(kind=8) dt  ! yr 
integer, parameter :: nt = 50000000
real(kind=8) time
integer, parameter :: nsp = 4

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
real(kind=8) :: ali = 0d0

real(kind=8) :: mvka = 99.52d0 ! cm3/mol; molar volume of kaolinite; Robie et al. 1978
real(kind=8) :: mvfo = 43.79d0 ! cm3/mol; molar volume of Fo; Robie et al. 1978
real(kind=8) :: mvab = 100.07d0 ! cm3/mol; molar volume of Ab(NaAlSi3O8); Robie et al. 1978 
real(kind=8) :: mvan = 100.79d0 ! cm3/mol; molar volume of An (CaAl2Si2O8); Robie et al. 1978
real(kind=8) :: mvcc = 36.934d0 ! cm3/mol; molar volume of Cc (CaCO3); Robie et al. 1978
real(kind=8) :: mvpy = 23.94d0 ! cm3/mol; molar volume of Pyrite (FeS2); Robie et al. 1978
real(kind=8) :: mvgb = 31.956d0 ! cm3/mol; molar volume of Gibsite (Al(OH)3); Robie et al. 1978
real(kind=8) :: mvct = 108.5d0 ! cm3/mol; molar volume of Chrysotile (Mg3Si2O5(OH)4); Robie et al. 1978
real(kind=8) :: mvfa = 46.39d0 ! cm3/mol; molar volume of Fayalite (Fe2SiO4); Robie et al. 1978
real(kind=8) :: mvgt = 20.82d0 ! cm3/mol; molar volume of Goethite (FeO(OH)); Robie et al. 1978
real(kind=8) :: mvcabd = 129.77d0 ! cm3/mol; molar volume of Ca-beidellite (Ca(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
real(kind=8) :: mvdp = 66.09d0 ! cm3/mol; molar volume of Diopside (MgCaSi2O6);  Robie et al. 1978
real(kind=8) :: mvhb = 248.09d0/3.55d0 ! cm3/mol; molar volume of Hedenbergite (FeCaSi2O6); from a webpage
real(kind=8) :: mvkfs = 108.72d0 ! cm3/mol; molar volume of K-feldspar (KAlSi3O8); Robie et al. 1978

real(kind=8) :: mwtfo = 140.694d0 ! g/mol; formula weight of Fo; Robie et al. 1978
real(kind=8) :: mwtab = 262.225d0 ! g/mol; formula weight of Ab; Robie et al. 1978
real(kind=8) :: mwtan = 278.311d0 ! g/mol; formula weight of An; Robie et al. 1978
real(kind=8) :: mwtcc = 100.089d0 ! g/mol; formula weight of Cc; Robie et al. 1978
real(kind=8) :: mwtpy = 119.967d0 ! g/mol; formula weight of Py; Robie et al. 1978
real(kind=8) :: mwtka = 100.089d0 ! g/mol; formula weight of Ka; Robie et al. 1978
real(kind=8) :: mwtgb = 78.004d0 ! g/mol; formula weight of Gb; Robie et al. 1978
real(kind=8) :: mwtct = 277.113d0 ! g/mol; formula weight of Ct; Robie et al. 1978
real(kind=8) :: mwtfa = 203.778d0 ! g/mol; formula weight of Fa; Robie et al. 1978
real(kind=8) :: mwtgt = 88.854d0 ! g/mol; formula weight of Gt; Robie et al. 1978
real(kind=8) :: mwtcabd = 366.6252667d0 ! g/mol; formula weight of Cabd calculated from atmoic weight
real(kind=8) :: mwtdp = 216.553d0 ! g/mol;  Robie et al. 1978
real(kind=8) :: mwthb = 248.09d0 ! g/mol; from a webpage
real(kind=8) :: mwtkfs = 278.33d0 ! g/mol; formula weight of Kfs; Robie et al. 1978

! real(kind=8) :: redsldi = 0.56d0 ! wt%  **default 
! real(kind=8) :: redsldi = 1.12d0 ! wt%  x2
real(kind=8) :: redsldi = 2.8d0 ! wt%   x5
! real(kind=8) :: redsldi = 2.24d0 ! wt%  x4
! real(kind=8) :: redsldi = 3.36d0 ! wt%  x6

! real(kind=8) :: silwti = 30d0 ! wt%  **default
! real(kind=8) :: silwti = 45d0 ! wt%  
! real(kind=8) :: silwti = 24d0 ! wt%
real(kind=8) :: silwti = 1d-10 ! wt%

! real(kind=8)::rainpowder = 40d2 !  g/m2/yr corresponding to 40 t/ha/yr (40x1e3x1e3/1e4)
! real(kind=8)::rainpowder = 0.5d2 !  g/m2/yr corresponding to 0.5 t/ha/yr (0.5x1e3x1e3/1e4)
real(kind=8),intent(in)::rainpowder != 30d2 !  g/m2/yr 
! real(kind=8)::rainpowder = 10d2 !  g/m2/yr corresponding to 10 t/ha/yr (0.5x1e3x1e3/1e4)

real(kind=8)::rainfrc_fo = 0.12d0 ! rain wt fraction for Fo (Beering et al 2020)
real(kind=8)::rainfrc_ab = 0.172d0 ! rain wt fraction for Ab; assuming 0.43 for La and 0.4 of wt of La is Ab (Beering et al 2020)
! real(kind=8)::rainfrc_ab = 0d0 ! rain wt fraction for Ab; assuming 0.43 for La and 0.4 of wt of La is Ab (Beering et al 2020)
real(kind=8)::rainfrc_an = 0.258d0 ! rain wt fraction for An; assuming 0.43 for La and 0.6 of wt of La is An (Beering et al 2020)
! real(kind=8)::rainfrc_an = 0d0 ! rain wt fraction for An; assuming 0.43 for La and 0.6 of wt of La is An (Beering et al 2020)
real(kind=8)::rainfrc_cc = 0d0 ! rain wt fraction for Cc; None (Beering et al 2020)
real(kind=8)::rainfrc_ka = 0d0 ! rain wt fraction for Ka; None (Beering et al 2020)
real(kind=8)::rainfrc_gb = 0d0 ! rain wt fraction for Gb; None (Beering et al 2020)
real(kind=8)::rainfrc_ct = 0d0 ! rain wt fraction for ct; None (Beering et al 2020)
real(kind=8)::rainfrc_fa = 0.05d0 ! rain wt fraction for Fa (Beering et al 2020)
real(kind=8)::rainfrc_dp = 0.189d0 ! rain wt fraction for dp; 0.21 for augite and assuming 0.9 of augite is from diopside (Beering et al 2020)
real(kind=8)::rainfrc_hb = 0.021d0 ! rain wt fraction for hb; 0.21 for augite and assuming 0.1 of augite is from hedenbergite (Beering et al 2020)
real(kind=8)::rainfrc_kfs = 0.06d0 ! rain wt fraction for kfs (Beering et al 2020)

real(kind=8),intent(in)::zsupp != 0.3d0 !  e-folding decrease

real(kind=8) sat(nz), poro(nz), torg(nz), tora(nz), deff(nz)
real(kind=8) :: dgaso = 6.09d2 ! m^2 yr^-1
real(kind=8) :: daqo = 5.49d-2 ! m^2 yr^-1
real(kind=8) :: dgasc = 441.504d0 ! m^2 yr^-1 (Assuming 0.14 cm2/sec)
real(kind=8) :: daqc = 0.022459852 ! m^2 yr^-1 (for C32- from Li and Gregory 1974)

! real(kind=8) :: poroi = 0.1d0 !*** default
real(kind=8),intent(in) :: poroi != 0.5d0

real(kind=8) :: sati = 0.50d0
real(kind=8),intent(in) :: satup != 0.10d0

! real(kind=8) :: zsat = 30d0  ! water table depth [m] ** default 
real(kind=8),intent(in) :: zsat != 5d0  ! water table depth [m] 
! real(kind=8) :: zsat = 15d0

real(kind=8) :: dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
real(kind=8) :: dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dmg  = 0.017218079d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dsi  = 0.03689712d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dca  = 0.019023312d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
real(kind=8) :: dal  = 0.011656226d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 

real(kind=8),intent(in) :: w != 5.0d-5 ! m yr^-1, uplift rate ** default 
! real(kind=8), parameter :: w = 1.0d-4 ! m yr^-1, uplift rate

real(kind=8), parameter :: vcnst = 1.0d1 ! m yr^-1, advection
! real(kind=8), parameter :: qin = 5d-3 ! m yr^-1, advection (m3 water / m2 profile / yr)

! real(kind=8) :: qin = 1d-1 ! m yr^-1, advection (m3 water / m2 profile / yr)  ** default
real(kind=8),intent(in) :: qin != 10d-1 ! m yr^-1 
! real(kind=8) :: qin = 0.1d-1 ! m yr^-1 
real(kind=8) v(nz), q

! real(kind=8) :: hr = 1d5 ! m^2 m^-3, reciprocal of hydraulic radius  ** default 
! real(kind=8) :: hr = 1d4 ! m^2 m^-3, reciprocal of hydraulic radius
real(kind=8) :: hrii = 1d5

! real(kind=8) :: p80 = 10d-6 ! m (**default?)
real(kind=8),intent(in) :: p80 != 1d-6 ! m 

real(kind=8) msi,msili,mfoi,mabi,mani,mcci,ctmp,po2tmp,ssa_cmn,mvab_save,mvan_save,mvcc_save,mvfo_save &
    & ,mkai,mvka_save,mvgb_save,mcti
real(kind=8),dimension(nz)::msil,msilx,mfo,mfox,mab,mabx,man,manx,mcc,mccx,mka,mkax,mgb,mgbx &
    & ,po2,redsld,ms,c,po2x,msx,cx,resp,c2,c2x,so4,so4x,na,nax,naeq,silsat &
    & ,pro,prox,hco3,ca,co2,co3,dic,cax,porox,dporodta,dporodtg,dporodtgc,khco2  &
    & ,mg,mgx,si,six,pco2,pco2x,poroprev,khco2x,pco2x_prev,torgprev,toraprev,hr,rough,hri,al,alx

real(kind=8) :: caeq = 1d-3  ! mol/L equilibrium Ca conc. 
real(kind=8) :: delca = 0.5d0  ! m reaction front width  
real(kind=8) :: zca = 50d0   ! m depth of reaction front for calcite          

real(kind=8),dimension(nz)::koxa,koxs,koxs2,ksil,msilsupp,kfo,mfosupp,omega_fo,omega_ka,mkasupp,mgbsupp &
    & ,kab,mabsupp,omega_ab,kan,mansupp,omega_an,kcc,mccsupp,omega_cc,preccc,kcca,omega_cca,kka,kgb &
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,po2supp

real(kind=8) kho,ucv,kco2,k1,keqsil,kw,k2,keqfo,keqab,keqgb,khco2i,keqan,keqcc,k1si,k2si,keqcca &
    & ,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,k1al,k2al,k3al,k4al,keqka

integer iz, ie, it, ie2, iter_co2,ispa,ispg,isps,irxn,ispa2,ispg2,isps2

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
real(kind=8) :: alth = 1.0d-20
real(kind=8) :: msth = 1.0d-20
real(kind=8) :: msilth = 1.0d-20
real(kind=8) :: mfoth = 1.0d-20
real(kind=8) :: mabth = 1.0d-20
real(kind=8) :: manth = 1.0d-20
real(kind=8) :: mccth = 1.0d-20
real(kind=8) :: mkath = 1.0d-20
real(kind=8) :: mgbth = 1.0d-20
real(kind=8) :: mctth = 1.0d-20

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
real(kind=8),dimension(nflx,nz)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_o2,flx_an,flx_ca,flx_cc,flx_co2,flx_ka &
    & ,flx_al,flx_gb
real(kind=8),dimension(nflx_py,nz)::flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2,flx_py_so4
! real(kind=8) :: maxdt = 10d0
real(kind=8) :: maxdt = 0.2d0 ! for basalt exp?

logical :: pre_calc = .false.
! logical :: pre_calc = .true.

logical :: read_data = .false.
! logical :: read_data = .true.

! logical :: initial_ss = .false.
logical :: initial_ss = .true.

! logical :: incld_rough = .false.
logical :: incld_rough = .true.

! logical :: cplprec = .false.
logical :: cplprec = .true.

logical :: rain_wave = .false.
! logical :: rain_wave = .true.

logical :: co2_iteration = .false.
! logical :: co2_iteration = .true.

logical :: calcite_seed = .false.
! logical :: calcite_seed = .true.

logical :: al_inhibit = .false.
! logical :: al_inhibit = .true.

real(kind=8) :: authig = 0d0 ! 0 if not allowing authigenesis of CaCO3
! real(kind=8) :: authig = 1d0 ! 1 if allowing authigenesis of CaCO3

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
logical :: flgreducedt = .false.
logical :: flgreducedt_prev = .false.
real(kind=8) :: zab(3), zpy(3) 

real(kind=8), parameter:: infinity = huge(0d0)
real(kind=8) k_arrhenius

real(kind=8) time_start, time_fin, progress_rate, progress_rate_prev
integer count_dtunchanged 

integer,intent(in)::nsp_sld != 5
integer,parameter::nsp_sld_2 = 6
integer,parameter::nsp_sld_all = 14
integer ::nsp_sld_cnst != nsp_sld_all - nsp_sld
integer,intent(in)::nsp_aq != 5
integer,parameter::nsp_aq_ph = 9
integer,parameter::nsp_aq_all = 9
integer ::nsp_aq_cnst != nsp_aq_all - nsp_aq
integer,intent(in)::nsp_gas != 2
integer,parameter::nsp_gas_ph = 1
integer,parameter::nsp_gas_all = 2
integer ::nsp_gas_cnst != nsp_gas_all - nsp_gas
integer ::nsp3 != nsp_sld + nsp_aq + nsp_gas
integer,intent(in)::nrxn_ext != 1
integer,parameter::nrxn_ext_all = 2
character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_2)::chrsld_2
character(5),dimension(nsp_sld_all)::chrsld_all
character(5),dimension(nsp_sld_all - nsp_sld)::chrsld_cnst
character(5),dimension(nsp_aq),intent(in)::chraq
character(5),dimension(nsp_aq_ph)::chraq_ph
character(5),dimension(nsp_aq_all)::chraq_all
character(5),dimension(nsp_aq_all - nsp_aq)::chraq_cnst
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_ph)::chrgas_ph
character(5),dimension(nsp_gas_all)::chrgas_all
character(5),dimension(nsp_gas_all - nsp_gas)::chrgas_cnst
character(5),dimension(nrxn_ext),intent(in)::chrrxn_ext
character(5),dimension(nrxn_ext_all)::chrrxn_ext_all
real(kind=8),dimension(nsp_sld)::msldi,msldth,mv,rfrc_sld,mwt
real(kind=8),dimension(nsp_sld,nsp_aq)::staq
real(kind=8),dimension(nsp_sld,nsp_gas)::stgas
real(kind=8),dimension(nsp_sld,nz)::msldx,msld,ksld,omega,msldsupp
real(kind=8),dimension(nsp_sld,nflx,nz)::flx_sld
real(kind=8),dimension(nsp_aq)::maqi,maqth,daq
real(kind=8),dimension(nsp_aq,nz)::maqx,maq,rxnaq,maqsupp
real(kind=8),dimension(nsp_aq,nflx,nz)::flx_aq
real(kind=8),dimension(nsp_gas)::mgasi,mgasth,dgasa,dgasg,dmgas,khgasi,dgasi
real(kind=8),dimension(nsp_gas,nz)::mgasx,mgas,khgasx,khgas,dgas,agasx,agas,rxngas,mgassupp 
real(kind=8),dimension(nsp_gas,nflx,nz)::flx_gas 
real(kind=8),dimension(nrxn_ext,nz)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_gas)::stgas_ext,stgas_dext
real(kind=8),dimension(nrxn_ext,nsp_aq)::staq_ext,staq_dext
real(kind=8),dimension(nrxn_ext,nsp_sld)::stsld_ext,stsld_dext

real(kind=8),dimension(nsp_aq_all)::daq_all,maqi_all,maqth_all
real(kind=8),dimension(nsp_gas_all)::dgasa_all,dgasg_all,mgasi_all,mgasth_all
real(kind=8),dimension(nsp_gas_all,3)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2)::keqaq_c
real(kind=8),dimension(nsp_sld_all,nz)::ksld_all
real(kind=8),dimension(nsp_sld_all,nsp_aq_all)::staq_all
real(kind=8),dimension(nsp_sld_all,nsp_gas_all)::stgas_all
real(kind=8),dimension(nsp_sld_all)::keqsld_all,mv_all,msldi_all,msldth_all,rfrc_sld_all,mwt_all
real(kind=8),dimension(nrxn_ext_all,nz)::krxn1_ext_all
real(kind=8),dimension(nrxn_ext_all,nz)::krxn2_ext_all
real(kind=8),dimension(nrxn_ext_all,nsp_aq_all)::staq_ext_all,staq_dext_all
real(kind=8),dimension(nrxn_ext_all,nsp_gas_all)::stgas_ext_all,stgas_dext_all
real(kind=8),dimension(nrxn_ext_all,nsp_sld_all)::stsld_ext_all,stsld_dext_all

real(kind=8),dimension(nsp_aq_all - nsp_aq,nz)::maqc
real(kind=8),dimension(nsp_gas_all - nsp_gas,nz)::mgasc
real(kind=8),dimension(nsp_sld_all - nsp_sld,nz)::msldc

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer,dimension(nsp_aq)::iaqflx
integer,dimension(nsp_gas)::igasflx
integer,dimension(nsp_sld)::isldflx

integer,parameter::ibasaltrain = 15
integer ::isldprof != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 1
integer ::iaqprof != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 2
integer ::igasprof != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 3
integer ::isldsat != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 4
integer ::ibsd != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 5

real(kind=8) dt_prev

logical print_cb,ph_error
character(500) print_loc
character(500),intent(in):: sim_name

real(kind=8) def_dust,def_rain,def_pr
!-------------------------

nsp_sld_cnst = nsp_sld_all - nsp_sld
nsp_aq_cnst = nsp_aq_all - nsp_aq
nsp_gas_cnst = nsp_gas_all - nsp_gas
nsp3 = nsp_sld + nsp_aq + nsp_gas

isldprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 1
iaqprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 2
igasprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 3
isldsat = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 4
ibsd = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 5

! define all species and rxns definable in the model 
! note that rxns here exclude diss(/prec) of mineral 
! which are automatically included when associated mineral is chosen

chrsld_all = (/'fo   ','ab   ','an   ','cc   ','ka   ','gb   ','py   ','ct   ','fa   ','gt   ','cabd ' &
    & ,'dp   ','hb   ','kfs  '/)
chraq_all = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    '/)
chrgas_all = (/'pco2','po2 '/)
chrrxn_ext_all = (/'resp ','fe2o2'/)

! define the species and rxns explicitly simulated in the model in a fully coupled way
! should be chosen from definable species & rxn lists above 

! chrsld = (/'fo   ','ab   ','an   ','cc   ','ka   '/)
! chraq = (/'mg   ','si   ','na   ','ca   ','al   '/)
! chrgas = (/'pco2 ','po2  '/)
! chrrxn_ext = (/'resp '/)


! define solid species which can precipitate
! in default, all minerals only dissolve 
! should be chosen from the chrsld list
chrsld_2 = (/'cc   ','ka   ','gb   ','ct   ','gt   ','cabd '/) 

! below are species which are sensitive to pH 
chraq_ph = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    '/)
chrgas_ph = (/'pco2 '/)

if (nsp_aq_cnst .ne. 0) then 
    do ispa = 1, nsp_aq_cnst
        do ispa2=1,nsp_aq_all
            if (.not.any(chraq==chraq_all(ispa2)) .and. .not.any(chraq_cnst==chraq_all(ispa2))) then 
                chraq_cnst(ispa) = chraq_all(ispa2)
                exit 
            endif 
        enddo
    enddo 
    print *, chraq_cnst
    ! pause
endif 

if (nsp_gas_cnst .ne. 0) then 
    do ispg = 1, nsp_gas_cnst
        do ispg2=1,nsp_gas_all
            if (.not.any(chrgas==chrgas_all(ispg2)) .and. .not.any(chrgas_cnst==chrgas_all(ispg2))) then 
                chrgas_cnst(ispg) = chrgas_all(ispg2)
                exit 
            endif 
        enddo
    enddo 
    print *, chrgas_cnst
    ! pause 
endif 

if (nsp_sld_cnst .ne. 0) then 
    do isps = 1, nsp_sld_cnst
        do isps2=1,nsp_sld_all
            if (.not.any(chrsld==chrsld_all(isps2)) .and. .not.any(chrsld_cnst==chrsld_all(isps2))) then 
                chrsld_cnst(isps) = chrsld_all(isps2)
                exit 
            endif 
        enddo
    enddo 
    print *, chrsld_cnst
    ! pause 
endif 

! molar volume 

mv_all = (/mvfo,mvab,mvan,mvcc,mvka,mvgb,mvpy,mvct,mvfa,mvgt,mvcabd,mvdp,mvhb,mvkfs/)
mwt_all = (/mwtfo,mwtab,mwtan,mwtcc,mwtka,mwtgb,mwtpy,mwtct,mwtfa,mwtgt,mwtcabd,mwtdp,mwthb,mwtkfs/)

do isps = 1, nsp_sld 
    mv(isps) = mv_all(findloc(chrsld_all,chrsld(isps),dim=1))
    mwt(isps) = mwt_all(findloc(chrsld_all,chrsld(isps),dim=1))
enddo 

! initial values for all species 
! mgasi_all = (/pco2i,po2i/)

! msldi_all = 1d-10
! msldi_all(findloc(chrsld_all,'ab',dim=1)) =  silwti*1d-2/262.2d0*2.7d0*(1.0d0-poroi)*1d6  
! msldi_all(findloc(chrsld_all,'py',dim=1)) =  redsldi*1d-2/120.0d0*2.7d0*(1.0d0-poroi)*1d6  

! maqi_all = 0d0
    
def_rain = 0d0
def_pr = 1d-20
    
call get_rainwater( &
    & nsp_aq_all,chraq_all,def_rain &! input
    & ,maqi_all &! output
    & )
    
call get_parentrock( &
    & nsp_sld_all,chrsld_all,def_pr &! input
    & ,msldi_all &! output
    & )

msldi_all = msldi_all/mwt_all*2.7d0*(1d0-poroi)*1d6

call get_atm( &
    & nsp_gas_all,chrgas_all &! input
    & ,mgasi_all &! output
    & )

print*,maqi_all 
print*,mgasi_all 
print*,msldi_all

pause

! constant values are taken from the boundary values specified above 
do ispg = 1,nsp_gas_cnst
    mgasc(ispg,:) = mgasi_all(findloc(chrgas_all,chrgas_cnst(ispg),dim=1))
enddo 
do ispa = 1,nsp_aq_cnst
    maqc(ispa,:) = maqi_all(findloc(chraq_all,chraq_cnst(ispa),dim=1))
enddo 
do isps = 1,nsp_sld_cnst
    msldc(isps,:) = msldi_all(findloc(chrsld_all,chrsld_cnst(isps),dim=1))
enddo 

! threshould values 
mgasth_all = 1d-20
maqth_all = 1d-20
msldth_all = 1d-20


! passing initial and threshold values to explcit variables 
do isps = 1, nsp_sld    
    if (any(chrsld_all == chrsld(isps))) then 
        msldi(isps) = msldi_all(findloc(chrsld_all,chrsld(isps),dim=1))
        msldth(isps) = msldth_all(findloc(chrsld_all,chrsld(isps),dim=1))
    endif 
enddo 
do ispa = 1, nsp_aq    
    if (any(chraq_all == chraq(ispa))) then 
        maqi(ispa) = maqi_all(findloc(chraq_all,chraq(ispa),dim=1))
        maqth(ispa) = maqth_all(findloc(chraq_all,chraq(ispa),dim=1))
    endif 
enddo 
do ispg = 1, nsp_gas    
    if (any(chrgas_all == chrgas(ispg))) then 
        mgasi(ispg) = mgasi_all(findloc(chrgas_all,chrgas(ispg),dim=1))
        mgasth(ispg) = mgasth_all(findloc(chrgas_all,chrgas(ispg),dim=1))
    endif 
enddo 

print*,maqi 
print*,mgasi
print*,msldi

pause

! stoichiometry
! mineral dissolution(/precipitation)
staq_all = 0d0
stgas_all = 0d0
! Forsterite; Mg2SiO4
staq_all(findloc(chrsld_all,'fo',dim=1), findloc(chraq_all,'mg',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'fo',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
! Albite; NaAlSi3O8
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! K-feldspar; KAlSi3O8
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! Anothite; CaAl2Si2O8
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
! Calcite; CaCO3
staq_all(findloc(chrsld_all,'cc',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'cc',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
! Kaolinite; Al2Si2O5(OH)4
staq_all(findloc(chrsld_all,'ka',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'ka',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
! Gibbsite; Al(OH)3
staq_all(findloc(chrsld_all,'gb',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! Pyrite; FeS2
staq_all(findloc(chrsld_all,'py',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'py',dim=1), findloc(chraq_all,'so4',dim=1)) = 2d0
stgas_all(findloc(chrsld_all,'py',dim=1), findloc(chrgas_all,'po2',dim=1)) = -7d0/2d0
! Chrysotile; Mg3Si2O5(OH)4
staq_all(findloc(chrsld_all,'ct',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'ct',dim=1), findloc(chraq_all,'mg',dim=1)) = 3d0
! Fayalite; Fe2SiO4
staq_all(findloc(chrsld_all,'fa',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'fa',dim=1), findloc(chraq_all,'fe2',dim=1)) = 2d0
! Goethite; FeO(OH)
staq_all(findloc(chrsld_all,'gt',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
! Ca-beidellite; Ca(1/6)Al(7/3)Si(11/3)O10(OH)2
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0/6d0
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
! Diopside (MgCaSi2O6)
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
! Hedenbergite (FeCaSi2O6)
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0

staq = 0d0
stgas = 0d0

do isps = 1, nsp_sld
    if (any(chrsld_all == chrsld(isps))) then 
        do ispa = 1, nsp_aq 
            if (any(chraq_all == chraq(ispa))) then 
                staq(isps,ispa) = &
                    & staq_all(findloc(chrsld_all,chrsld(isps),dim=1), findloc(chraq_all,chraq(ispa),dim=1))
            endif 
        enddo 
        do ispg = 1, nsp_gas 
            if (any(chrgas_all == chrgas(ispg))) then 
                stgas(isps,ispg) = &
                    & stgas_all(findloc(chrsld_all,chrsld(isps),dim=1), findloc(chrgas_all,chrgas(ispg),dim=1))
            endif 
        enddo 
    endif 
enddo 

! external reactions
staq_ext_all = 0d0
stgas_ext_all = 0d0
stsld_ext_all = 0d0
! respiration 
stgas_ext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_ext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
! fe2 oxidation 
staq_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe2',dim=1)) = -1d0
staq_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0
stgas_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

! define 1 when a reaction is sensitive to a speces 
stgas_dext_all = 0d0
staq_dext_all = 0d0
stsld_dext_all = 0d0
! respiration 
stgas_dext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
! fe2 oxidation 
stgas_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
staq_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0

staq_ext = 0d0
stgas_ext = 0d0
stsld_ext = 0d0

do irxn = 1, nrxn_ext
    if (any(chrrxn_ext_all == chrrxn_ext(irxn))) then 
        do ispa = 1, nsp_aq 
            if (any(chraq_all == chraq(ispa))) then 
                staq_ext(irxn,ispa) = &
                    & staq_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chraq_all,chraq(ispa),dim=1))
                staq_dext(irxn,ispa) = &
                    & staq_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chraq_all,chraq(ispa),dim=1))
            endif 
        enddo 
        do ispg = 1, nsp_gas 
            if (any(chrgas_all == chrgas(ispg))) then 
                stgas_ext(irxn,ispg) = &
                    & stgas_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chrgas_all,chrgas(ispg),dim=1))
                stgas_dext(irxn,ispg) = &
                    & stgas_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chrgas_all,chrgas(ispg),dim=1))
            endif 
        enddo 
        do isps = 1, nsp_sld 
            if (any(chrsld_all == chrsld(isps))) then 
                stsld_ext(irxn,isps) = &
                    & stsld_ext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chrsld_all,chrsld(isps),dim=1))
                stsld_dext(irxn,isps) = &
                    & stsld_dext_all(findloc(chrrxn_ext_all,chrrxn_ext(irxn),dim=1) &
                    &   ,findloc(chrsld_all,chrsld(isps),dim=1))
            endif 
        enddo 
    endif 
enddo 

! rfrc_sld_all = 0d0 
! rfrc_sld_all(findloc(chrsld_all,'fo',dim=1)) = rainfrc_fo/mwtfo
! rfrc_sld_all(findloc(chrsld_all,'an',dim=1)) = rainfrc_an/mwtan
! rfrc_sld_all(findloc(chrsld_all,'ab',dim=1)) = rainfrc_ab/mwtab
! rfrc_sld_all(findloc(chrsld_all,'fa',dim=1)) = rainfrc_fa/mwtfa
! rfrc_sld_all(findloc(chrsld_all,'dp',dim=1)) = rainfrc_dp/mwtdp
! rfrc_sld_all(findloc(chrsld_all,'hb',dim=1)) = rainfrc_hb/mwthb
! rfrc_sld_all(findloc(chrsld_all,'kfs',dim=1)) = rainfrc_kfs/mwtkfs

def_dust = 0d0
    
call get_dust( &
    & nsp_sld_all,chrsld_all,def_dust &! input
    & ,rfrc_sld_all &! output
    & )

rfrc_sld_all = rfrc_sld_all/mwt_all


do isps = 1, nsp_sld 
    rfrc_sld(isps) = rfrc_sld_all(findloc(chrsld_all,chrsld(isps),dim=1))
enddo


do while (rectime(nrec)>ttot) 
    rectime = rectime/10d0
enddo 
do while (rectime(nrec)<ttot) 
    rectime = rectime*10d0
enddo 
! rectime =rectime/1d4 
! rectime =rectime/1d3 ! better with basalt exp? max 12000 yr 
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


write(workdir,*) '../pyweath_output/'     

! if (cplprec) then 
    ! write(base,*) 'test_cplp_test'
! else 
    ! write(base,*) 'test_cpl'
! endif 

base = trim(adjustl(sim_name))

if (al_inhibit) base = trim(adjustl(base))//'_alx'

base = trim(adjustl(base))//'_rain-'//trim(adjustl(chrrain))    
 
#ifdef poroevol 
base = trim(adjustl(base))//'_pevol'
#endif 
#if defined(surfevol1)
base = trim(adjustl(base))//'_sevol1'
#elif defined(surfevol2)
base = trim(adjustl(base))//'_sevol2'
#elif defined(surfssa)
base = trim(adjustl(base))//'_ssa'
#endif 

#ifndef regulargrid
base = trim(adjustl(base))//'_irr'
#endif 

if (rain_wave)then 
    write(chrrain,'(E10.2)') wave_tau
    base = trim(adjustl(base))//'_rwave-'//trim(adjustl(chrrain))
endif 

if (incld_rough)then 
    write(chrrain,'(E10.2)') p80
    base = trim(adjustl(base))//'_p80r-'//trim(adjustl(chrrain))
else
    write(chrrain,'(E10.2)') p80
    base = trim(adjustl(base))//'_p80-'//trim(adjustl(chrrain))
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

do isps = 1, nsp_sld 
    isldflx(isps) = ibasaltrain + isps
enddo 
    
do ispa = 1, nsp_aq 
    iaqflx(ispa) = ibasaltrain + nsp_sld  + ispa
enddo 

do ispg = 1, nsp_gas
    igasflx(ispg) = ibasaltrain + nsp_sld + nsp_aq + ispg
enddo 

call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))

do isps = 1,nsp_sld
    open(isldflx(isps), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', status='replace')
    write(isldflx(isps),*) ' time ',' itflx ',' iadv ',' idif ',' irxn ',' irain ',' ires '
    close(isldflx(isps))
enddo 

do ispa = 1,nsp_aq
    open(iaqflx(ispa), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', status='replace')
    write(iaqflx(ispa),*) ' time ',' itflx ',' iadv ',' idif ',' irxn ',' irain ',' ires '
    close(iaqflx(ispa))
enddo 

do ispg = 1,nsp_gas
    open(igasflx(ispg), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', status='replace')
    write(igasflx(ispg),*) ' time ',' itflx ',' iadv ',' idif ',' irxn ',' irain ',' ires '
    close(igasflx(ispg))
enddo 
open(ibasaltrain, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'rain.txt', &
    & status='replace')
close(ibasaltrain)



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

rough = 1d0
if (incld_rough) then 
    rough = 10d0**(3.3d0)*p80**0.33d0 ! from Navarre-Sitchler and Brantley (2007)
endif 

ssa_cmn = -4.4528d0*log10(p80*1d6) + 11.578d0 ! m2/g

hrii = 1d0/p80
hri = hrii

hr = hri*rough
v = vcnst
v = qin/poroi/sat
poro = poroi
torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
deff = torg*dgaso + tora*daqo

#ifdef surfssa
hri = ssa_cmn*1d6/poro
mvab_save = mvab
mvan_save = mvan
mvcc_save = mvcc
mvfo_save = mvfo
mvka_save = mvka
mvab = mwtab 
mvan = mwtan 
mvcc = mwtcc 
mvfo = mwtfo 
mvka = mwtka 
#endif 

dt = maxdt
if (.not.initial_ss) then 
    dt = 1d-20  
else
    dt = 1d-100 ! for basalt exp?
endif 

pro = 1d-5
    
call coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 

do ispa = 1, nsp_aq
    maq(ispa,:)=maqi(ispa)
enddo 
do ispg = 1, nsp_gas
    mgas(ispg,:)=mgasi(ispg)
enddo 
do isps = 1, nsp_sld
    msld(isps,:) = msldi(isps)
enddo 

omega = 0d0

print_cb = .false. 
print_loc = './ph.txt'

call calc_pH_v5( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
    & ,print_cb,print_loc,z &! input 
    & ,pro,ph_error &! output
    & ) 

poroprev = poro

! when adv o2 flux == adv pyrite flx (q*kho*po2i*1d3 = stoxs*msi*w)
if (swadvmass == 1d0) then 
    q = 15d0/4d0*msi*w/kho/po2i/1d3
    v = q/poroi/sat
endif 
! *** the above must be commented out when using arbitrary q value


!  --------- read -----
if (read_data) then 
    ! runname_save = 'test_cpl_rain-0.40E+04_pevol_sevol1_q-0.10E-01_zsat-5' ! specifiy the file where restart data is stored 
    runname_save = runname  ! the working folder has the restart data 
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'o2profile-res-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt')
    open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt',  &
        & status ='old',action='read')

    do iz = 1, Nz
        read (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),al(iz) &
            & ,mab(iz),man(iz),mfo(iz),mka(iz),mcc(iz) &
            & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_ka(iz),omega_cc(iz),pco2(iz),pro(iz),time
    enddo 
    close(22) 
    pro = 10d0**(-pro) ! read data is -log10 (pro)
    time = 0d0
    
    do iz=1,nz
        if (po2(iz)<po2th) po2(iz)= po2th*0.1d0
        if (c(iz)<cth) c(iz)= cth*0.1d0
        if (ms(iz)<msth) ms(iz)= msth*0.1d0
        if (c2(iz)<c2th) c2(iz)= c2th*0.1d0
        if (so4(iz)<so4th) so4(iz)= so4th*0.1d0
        if (na(iz)<nath) na(iz)= nath*0.1d0
        if (ca(iz)<cath) ca(iz)= cath*0.1d0
        if (mg(iz)<mgth) mg(iz)= mgth*0.1d0
        if (si(iz)<sith) si(iz)= sith*0.1d0
        if (al(iz)<sith) al(iz)= alth*0.1d0
        if (mab(iz)<mabth) mab(iz)= mabth*0.1d0
        if (man(iz)<manth) man(iz)= manth*0.1d0
        if (mfo(iz)<mfoth) mfo(iz)= mfoth*0.1d0
        if (mcc(iz)<mccth) mcc(iz)= mccth*0.1d0
        if (mka(iz)<mkath) mka(iz)= mkath*0.1d0
        if (pco2(iz)<pco2th) pco2(iz)= pco2th*0.1d0
    enddo
    
    ! man = mani
    ! mfo = mfoi
    ! manx = man
    ! mfox = mfo
    if (calcite_seed) then 
        mcc = mcci ! to precipirate calcite, crystal seeds are necessary 
        mccx = mcc
    endif 
    
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
    print *, 'al:', (al(iz),iz=1,nz, nz/5)
    print *, 'fo:', (mfo(iz),iz=1,nz, nz/5)
    print *, 'ab:', (mab(iz),iz=1,nz, nz/5)
    print *, 'an:', (man(iz),iz=1,nz, nz/5)
    print *, 'ka:', (mka(iz),iz=1,nz, nz/5)
    print *, 'cc:', (mcc(iz),iz=1,nz, nz/5)
    print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    print *, 'omega_ka:', (omega_ka(iz),iz=1,nz, nz/5)
    print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(pro(iz)),iz=1,nz, nz/5)
    print *
#endif      
endif
    
call coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
    
! --------- loop -----
print *, 'about to start time loop'
it = 0
irec = 0

count_dtunchanged = 0

!! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@

do while (it<nt)
    call cpu_time(time_start)
#ifdef display 
    print *, 'it, time = ',it, time
#endif
    dt_prev = dt
    ! if (time>rectime(nrec)) exit
    if (initial_ss.and.time>rectime(nrec)) exit
        
    if (.not.initial_ss .and. time > ztot/w*2d0) then 
        initial_ss = .true.
        time = 0
        it = 0
        dt = 1d-300
        ! pause
        
        ! man = mani
        ! mfo = mfoi
        ! manx = man
        ! mfox = mfo
        if (calcite_seed) then 
            mcc = mcci
            mccx = mcc
        endif 
    endif 

    if (.not.initial_ss .and. it==0) then 
    ! if (.not.initial_ss) then 
        maxdt = 1d2
        maxdt = 1d1
    else if (initial_ss .and. it==0) then
    ! else if (initial_ss) then
        maxdt = 0.2d0
        ! maxdt = 0.02d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.005d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.002d0 ! working with p80 = 10 um
        ! maxdt = 0.001d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.0005d0 ! working with p80 = 1 um
        
        ! if (time<1d-2) then  
            ! maxdt = 1d-6 
        ! elseif (time>=1d-2 .and. time<1d-1) then 
        
        ! if (time<1d-1) then  
            ! maxdt = 1d-5 
        ! elseif (time>=1d-1 .and. time<1d0) then  
            ! maxdt = 1d-4 
        ! elseif (time>=1d0 .and. time<1d1) then 
            ! maxdt = 1d-3 
        ! elseif (time>=1d1 .and. time<1d2) then 
            ! maxdt = 1d-2  
        ! elseif (time>=1d2 .and. time<1d3) then 
            ! maxdt = 1d-1 
        ! elseif (time>=1d3 .and. time<1d4) then 
            ! maxdt = 1d0 
        ! elseif (time>=1d4 .and. time<1d5) then 
            ! maxdt = 1d1 
        ! elseif (time>=1d5 ) then 
            ! maxdt = 1d2 
        ! endif 
        
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
    
    if (dt/=dt_prev) pre_calc = .true.

    ! incase temperature&ph change
    ! call coefs( &
        ! & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
        ! & ,pro,dgaso,dgasc,daqo,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
        ! & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
        ! & ,k1si,k2si,kcca,keqcca,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! output
        ! & ,k1al,k2al,k3al,k4al,keqka,kka,dal,keqgb,kgb &! output
        ! & ) 
        
    call coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
        & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
        & ) 
    
    do isps = 1, nsp_sld
        ksld(isps,:) = ksld_all(findloc(chrsld_all,chrsld(isps),dim=1),:)
        ! print *,chrsld(isps),ksld(isps,:)
    enddo
    
    do ispa = 1, nsp_aq 
        daq(ispa) = daq_all(findloc(chraq_all,chraq(ispa),dim=1))
        ! print *,chraq(ispa),daq(ispa)
    enddo 
    
    do ispg = 1, nsp_gas 
        dgasa(ispg) = dgasa_all(findloc(chrgas_all,chrgas(ispg),dim=1))
        dgasg(ispg) = dgasg_all(findloc(chrgas_all,chrgas(ispg),dim=1))
    enddo 
    kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)
    kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
    k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
    k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)
    khco2i = kco2*(1d0+k1/sqrt(kco2*k1*pco2i)+k2/kco2/pco2i)

    khgasi = (/khco2i,kho/)
    
    ! kinetic inhibition 
    if (al_inhibit) then 
        if (any(chraq == 'al')) then 
            do isps = 1, nsp_sld
                if (staq(isps,findloc(chraq,'al',dim=1)) .ne. 0d0) then 
                    ksld(isps,:) = ksld(isps,:) &
                        & *10d0**(-4.84d0)/(10d0**(-4.84d0)+maq(findloc(chraq,'al',dim=1),:)) 
                endif 
            enddo 
        endif 
    endif 

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

    mgasx = mgas
    msldx = msld
    maqx = maq
    
    prox = pro  

    porox = poro


    if (initial_ss) then 
        ! mfosupp = rainpowder*rainfrc_fo/mwtfo*exp(-z/zsupp)/zsupp
        ! mabsupp = rainpowder*rainfrc_ab/mwtab*exp(-z/zsupp)/zsupp 
        ! mansupp = rainpowder*rainfrc_an/mwtan*exp(-z/zsupp)/zsupp
        ! mccsupp = 0d0
        ! mkasupp = 0d0
        ! if (.not.cplprec)then
            ! mccsupp = kcc*poro*hr*mvcc*1d-6*mccx*(omega_cc - 1d0) &
                ! & *merge(1d0,0d0,omega_cc - 1d0 > 0d0) 
            ! mkasupp = kka*poro*hr*mvka*1d-6*mkax*(omega_ka - 1d0) &
                ! & *merge(1d0,0d0,omega_ka - 1d0 > 0d0) 
        ! endif 
            
        ! kcc = 0d0
        
        ! alsupp = -2d0*mkasupp
        ! mgsupp = 0d0
        ! nasupp = 0d0
        ! casupp = -mccsupp
        ! sisupp = -2d0*mkasupp
        
        ! pco2supp = -mccsupp
        ! po2supp = 0d0
        
        maqsupp = 0d0
        mgassupp = 0d0
        do isps = 1, nsp_sld
            msldsupp(isps,:) = rainpowder*rfrc_sld(isps)*exp(-z/zsupp)/zsupp
        enddo 
        
        
        if (rain_wave) then 
            do isps = 1, nsp_sld
                msldsupp(isps,:) = msldsupp(isps,:)*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
            enddo 
            if (time==0d0 .or. rain_norm /= merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))) then
                open(ibasaltrain, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'rain.txt', &
                    & status='old',action='write',access='append')
                write(ibasaltrain,*) time-dt,rain_norm
                write(ibasaltrain,*) time,merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                rain_norm = merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                close(ibasaltrain)
            endif 
        endif 
    else 
        mgassupp = 0d0
        msldsupp = 0d0
        maqsupp = 0d0
    endif 

    if (it==0) pre_calc = .true.
    pre_calc = .true.
    
    if (pre_calc) then 
    ! if (pre_calc .and. it ==0) then 
        pre_calc = .false.
        ! call precalc_po2_v2( &
            ! & nz,po2th,dt,ucv,kho,dz,dgaso,daqo,po2i,poro,sat,po2,torg,tora,v &! input 
            ! & ,po2x &! output 
            ! & )
        
        ! call precalc_pco2_v2( &
            ! & nz,pco2th,dt,ucv,khco2,dz,dgasc,daqc,pco2i,poro,sat,pco2,torg,tora,v,resp &! input 
            ! & ,pco2x &! output 
            ! & )

        ! mgas(findloc(chrgas,'pco2',dim=1),:)=pco2(:)
        ! mgas(findloc(chrgas,'po2',dim=1),:)=po2(:)
        
        call precalc_gases( &
            & nz,dt,ucv,dz,poro,sat,torg,tora,v,prox &! input 
            & ,nsp_gas,nsp_gas_all,chrgas,chrgas_all,keqgas_h,mgasi,mgasth,mgas &! input
            & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg &! input
            & ,mgasx &! output 
            & )
        
        ! call precalc_slds( &
            ! & nz,msth,dt,w,dz,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth   &! input
            ! & ,ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp,kcc,omega_cc,mvcc &! input
            ! & ,poro,hr,kcca,omega_cca,authig,sat,kka,mkai,mkath,omega_ka,mvka,mkasupp,mka &! input
            ! & ,msx,msilx,mfox,mabx,manx,mccx,mkax &! output
            ! & )
        
        ! msld(findloc(chrsld,'fo',dim=1),:)=mfo(:)
        ! msld(findloc(chrsld,'ab',dim=1),:)=mab(:)
        ! msld(findloc(chrsld,'an',dim=1),:)=man(:)
        ! msld(findloc(chrsld,'cc',dim=1),:)=mcc(:)
        ! msld(findloc(chrsld,'ka',dim=1),:)=mka(:)
        
        call precalc_slds_v2( &
            & nz,dt,w,dz,poro,hr,sat &! input
            & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
            & ,msldx &! output
            & )

        ! pause

        if (swoxa == 1d0) then 
        
            call precalc_pw_py( &
                & nz,dt,cth,c2th,v,c2,dz,dfe3,poro,sat,tora,koxa,po2x,c,dfe2,koxs,koxs2,hr,mvpy &! input
                & ,ms,po2,so4,so4th,ci,c2i,so4i,dso4,msx &! input 
                & ,c2x,cx,so4x  &! output
                & )
            
        endif 
        
        ! call precalc_pw_sil_v2( &
            ! & nz,nath,mgth,cath,sith,dt,v,na,ca,mg,si,dz,dna,dsi,dmg,dca,tora,poro,sat,nai,mgi,cai,sii &! input 
            ! & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx,alth,al,dal,ali,kka,mvka,mkax &! input 
            ! & ,nax,six,cax,mgx,alx &! output
            ! & )

        ! maq(findloc(chraq,'mg',dim=1),:)=mg(:)
        ! maq(findloc(chraq,'si',dim=1),:)=si(:)
        ! maq(findloc(chraq,'na',dim=1),:)=na(:)
        ! maq(findloc(chraq,'ca',dim=1),:)=ca(:)
        ! maq(findloc(chraq,'al',dim=1),:)=al(:)

        call precalc_aqs( &
            & nz,dt,v,dz,tora,poro,sat,hr &! input 
            & ,nsp_aq,nsp_sld,daq,maqth,maqi,maq,mv,msldx,ksld,staq &! input
            & ,nrxn_ext,staq_ext,rxnext &! input
            & ,maqx &! output
            & )

        ! if (any(isnan(po2x)).or.any(isnan(cx)).or.any(isnan(c2x)).or.any(isnan(pco2x))) then 
            ! print*, 'error in precalc'
            ! stop
        ! endif

        if (any(isnan(mgasx)).or.any(isnan(msldx)).or.any(isnan(maqx))) then 
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
        alx(:) = 1.0d2
    end if

    call pyweath_1D( &
        & nz,nflx_py,mvpy,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
        & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daqo,dgaso,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
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
        ! nax(1:) = 1.0d2
        ! mgx(1:) = 1.0d2
        ! six(1:) = 1.0d2
        ! cax(1:) = 1.0d2
        ! alx(1:) = 1.0d2
        maqx(:,1:) = 1d2
    endif 

#endif      

    error_co2 = 1d4
    iter_co2 = 0 
#ifdef two_way
    do while (error_co2> 1e-3) 
        call silicate_dis_1D_v2( &
            & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
            & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
            & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
            & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
            & ,iter,error,dt,flgback &! inout
            & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
            & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca &! output
            & )
            
        call oxygen_resp_1D_v2( &
            & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
            & ,kho,tora,torg,daqo,dgaso,v,mo2,tol,runname,workdir,zrxn,ucv,vmax,poroprev  &! inpput
            & ,dt,flgback &! inout
            & ,po2x,flx_o2,resp &! output
            & ) 

        pco2x_prev = pco2x
        preccc = flx_cc(4,:)
        ! preccc = 0d0
        call CO2_1D_v2_1( &
            & nz,nflx,pco2,pco2i,pco2th,poro,z,dz,sat,dporodtgc,v  &! input
            & ,kco2,k1,k2,tora,torg,daqc,dgasc,resp,tol,runname,workdir,zrxn,ucv,prox,khco2i,preccc,poroprev,pro  &! inpput
            & ,dt,flgback &! inout
            & ,pco2x,flx_co2 &! output
            & ) 
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
#endif 
    ! call silicate_dis_co2_1D_v2( &
        ! & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
        ! & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
        ! & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
        ! & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
        ! & ,pco2th,resp &! intput
        ! & ,iter,error,dt,flgback &! inout
        ! & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
        ! & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2 &! output
        ! & )
            
    ! call oxygen_resp_1D_v2( &
        ! & nz,nflx,po2,po2i,po2th,poro,z,dz,sat,dporodtg  &! input
        ! & ,kho,tora,torg,daqo,dgaso,v,mo2,tol,runname,workdir,zrxn,ucv,vmax,poroprev  &! inpput
        ! & ,dt,flgback &! inout
        ! & ,po2x,flx_o2,resp &! output
        ! & ) 
    ! call alsilicate_dis_co2_1D_v3( &
        ! & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
        ! & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
        ! & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
        ! & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
        ! & ,pco2th,resp,k1al,k2al,k3al,k4al,keqka,kka,al,mvka,ali,mkai,alth,mkath,dal,mkasupp,mka &! intput
        ! & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,cplprec &! input
        ! & ,iter,error,dt,flgback &! inout
        ! & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
        ! & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2,alx,flx_ka,flx_al,omega_ka,mkax &! output
        ! & )

    ! call alsilicate_aq_gas_1D( &
        ! & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
        ! & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
        ! & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
        ! & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
        ! & ,pco2th,k1al,k2al,k3al,k4al,keqka,kka,al,mvka,ali,mkai,alth,mkath,dal,mkasupp,mka &! intput
        ! & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,cplprec,po2,po2th,po2i,kho,po2supp,dgaso,daqo,vmax,mo2 &! input
        ! & ,iter,error,dt,flgback &! inout
        ! & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
        ! & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2,alx,flx_ka,flx_al,omega_ka,mkax,po2x,flx_o2,resp &! output
        ! & )

    ! call alsilicate_aq_gas_1D_v2( &
        ! & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
        ! & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
        ! & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
        ! & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
        ! & ,stgas_ext,stgas_dext,staq_ext,stsld_ext &
        ! & ,hr,poro,z,dz,w,keqfo,keqab,keqan,keqcc,sat,pro,poroprev &
        ! & ,kco2,k1,k2,tora,v,tol,it,so4x,nflx,kw,k1si,k2si & 
        ! & ,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,ucv,torg &
        ! & ,k1al,k2al,k3al,k4al,keqka,cplprec,kho,vmax,mo2  &
        ! & ,iter,error,dt,flgback &    
        ! & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,co2,hco3,co3,dic & 
        ! & )
    
    call alsilicate_aq_gas_1D_v3( &
        ! new input 
        & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
        & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
        & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
        & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
        & ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext &
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all &
        & ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &
        & ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &
        !  old inputs
        & ,hr,poro,z,dz,w,sat,pro,poroprev,tora,v,tol,it,nflx,kw & 
        & ,ucv,torg,cplprec  &
        ! old inout
        & ,iter,error,dt,flgback &    
        ! output 
        & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,co2,hco3,co3,dic & 
        & )

    ! if (iter > 75) then
        ! maxdt = maxdt/2d0
    ! end if
    ! if (iter<5) then 
        ! maxdt = maxdt*2d0
        ! if (maxdt >1d2) maxdt = 1d2
    ! endif 
    ! print*,iter,maxdt

    if (flgback) then 
        flgback = .false. 
        flgreducedt = .true.
        go to 100
    endif    
    
    dporodtg = 0d0
    dporodtgc = 0d0
    dporodta = 0d0
#ifdef poroevol   
    poroprev = poro
    torgprev = torg
    toraprev = tora
#ifdef surfssa
    mvab = mvab_save 
    mvan = mvan_save 
    mvcc = mvcc_save 
    mvfo = mvfo_save 
    mvka = mvka_save 
#endif 
    ! poro = poroi + (mabi-mabx)*(mvab)*1d-6  &
        ! & +(mfoi-mfox)*(mvfo)*1d-6 &
        ! & +(mani-manx)*(mvan)*1d-6 &
        ! & +(mcci-mccx)*(mvcc)*1d-6 &
        ! & +(mkai-mkax)*(mvka)*1d-6 
    poro = poroi
    do isps=1,nsp_sld
        poro = poro + (msldi(isps)-msldx(isps,:))*mv(isps)*1d-6
    enddo 
#ifdef surfssa
    mvab = mwtab 
    mvan = mwtan 
    mvcc = mwtcc 
    mvfo = mwtfo 
    mvka = mwtka 
#endif 
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
    deff = torg*dgaso + tora*daqo
    dporodtg = ( &
        & (ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3) &
        & -(ucv*porox*(1.0d0-sat)*1d3+porox*sat*kho*1d3) &
        & )/dt
    dporodtgc = ( &
        & (ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2*1d3) &
        & -(ucv*porox*(1.0d0-sat)*1d3+porox*sat*khco2*1d3) &
        & )/dt
    dporodta = (poro*sat-porox*sat)/dt/(poro*sat)
    hr = hri*rough
#ifdef surfevol1 
    hr = hri*rough*((1d0-poro)/(1d0-poroi))**(2d0/3d0)
#endif 
#ifdef surfevol2 
    hr = hri*rough*(poro/poroi)**(2d0/3d0)  ! SA increases with porosity 
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
    ! print *
    ! print *,'-=-=-=-=-=-= o2 & pyrite -=-=-=-=-=-=-='
    ! print *,'o2:', (po2x(iz),iz=1,nz, nz/5)
    ! print *,'fe2:', (cx(iz),iz=1,nz, nz/5)
    ! print *,'py:', (msx(iz),iz=1,nz, nz/5)
    ! print *, 'fe3:', (c2x(iz),iz=1,nz, nz/5)
    ! print *, 'so4:', (so4x(iz),iz=1,nz, nz/5)
    ! print *
    ! print *,'-=-=-=-=-=-= Na & albite -=-=-=-=-=-=-='
    ! print *, 'na:', (nax(iz),iz=1,nz, nz/5)
    ! print *, 'sil:', (msilx(iz),iz=1,nz, nz/5)
    ! print *
    ! print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An -=-=-=-=-=-=-='
    ! print *, 'mg:', (mgx(iz),iz=1,nz, nz/5)
    ! print *, 'si:', (six(iz),iz=1,nz, nz/5)
    ! print *, 'na:', (nax(iz),iz=1,nz, nz/5)
    ! print *, 'ca:', (cax(iz),iz=1,nz, nz/5)
    ! print *, 'al:', (alx(iz),iz=1,nz, nz/5)
    ! print *, 'fo:', (mfox(iz),iz=1,nz, nz/5)
    ! print *, 'ab:', (mabx(iz),iz=1,nz, nz/5)
    ! print *, 'an:', (manx(iz),iz=1,nz, nz/5)
    ! print *, 'ka:', (mkax(iz),iz=1,nz, nz/5)
    ! print *, 'cc:', (mccx(iz),iz=1,nz, nz/5)
    ! print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    ! print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    ! print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    ! print *, 'omega_ka:', (omega_ka(iz),iz=1,nz, nz/5)
    ! print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
    ! if (authig==1d0) print *, 'omega_cca:', (omega_cca(iz),iz=1,nz, nz/5)
    ! print *
    ! print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    ! print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    ! print *
    ! print *,'-=-=-=-=-=-= CO2 -=-=-=-=-=-=-='
    ! print *, 'co2:', (pco2x(iz),iz=1,nz, nz/5)
    ! print *
    
    
    print *
    print *,'-=-=-=-=-=-= Aq species -=-=-=-=-=-=-='
    do ispa = 1, nsp_aq
        print *, trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/5)
    enddo 
    print *
    print *,'-=-=-=-=-=-= Sld species -=-=-=-=-=-=-='
    do isps = 1, nsp_sld
        print *, trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *
    do isps = 1, nsp_sld
        print *, 'omega_'//trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *
    print *,'-=-=-=-=-=-= Gas species -=-=-=-=-=-=-='
    do ispg = 1, nsp_gas
        print *, trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/5)
    enddo 
    print *
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
#endif 

    ! stop
    
    mgas = mgasx
    maq = maqx
    msld = msldx
    
    pro = prox
        
    do iflx = 1,nflx_py
        flx_py(iflx,:) = flx_py(iflx,:)*dz
        flx_py_o2(iflx,:) = flx_py_o2(iflx,:)*dz
        flx_py_fe2(iflx,:) = flx_py_fe2(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_py_fe3(iflx,:) = flx_py_fe3(iflx,:)*dz*poro(:)*sat(:)*1d3
        flx_py_so4(iflx,:) = flx_py_so4(iflx,:)*dz*poro(:)*sat(:)*1d3
    enddo
        
    do iflx = 1,nflx
        flx_mg(iflx,:) = flx_mg(iflx,:)*dz
        flx_si(iflx,:) = flx_si(iflx,:)*dz
        flx_na(iflx,:) = flx_na(iflx,:)*dz
        flx_ca(iflx,:) = flx_ca(iflx,:)*dz
        flx_al(iflx,:) = flx_al(iflx,:)*dz
        flx_fo(iflx,:) = flx_fo(iflx,:)*dz
        flx_ab(iflx,:) = flx_ab(iflx,:)*dz
        flx_an(iflx,:) = flx_an(iflx,:)*dz
        flx_cc(iflx,:) = flx_cc(iflx,:)*dz
        flx_ka(iflx,:) = flx_ka(iflx,:)*dz
        flx_o2(iflx,:) = flx_o2(iflx,:)*dz
        flx_co2(iflx,:) = flx_co2(iflx,:)*dz
    enddo 

    if ((.not.initial_ss) .and. time > savetime) then 
        open (22, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-res-save.txt',  &
            & status='replace')

        open(30, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'o2profile-bsd-save.txt',  &
            & status='replace')

        do iz = 1, Nz
            write (22,*) z(iz),po2(iz),c(iz),ms(iz),c2(iz), so4(iz),na(iz),ca(iz),mg(iz),si(iz),al(iz) &
                & ,mab(iz),man(iz),mfo(iz),mka(iz),mcc(iz) &
                & ,omega_ab(iz), omega_fo(iz),omega_an(iz),omega_ka(iz),omega_cc(iz),pco2(iz),-log10(pro(iz)),time
            write(30,*) z(iz), poro(iz),sat(iz),v(iz),deff(iz),hr(iz)
        enddo 
        close(22)
        close(30)
        savetime = savetime + dsavetime
    endif 

    if (initial_ss .and. time>=rectime(irec+1)) then
        write(chr,'(i3.3)') irec+1
        
        
        print_cb = .true. 
        print_loc = trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'chrge_balance-'//chr//'.txt'

        call calc_pH_v5( &
            & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
            & ,print_cb,print_loc,z &! input 
            & ,prox,ph_error &! output
            & ) 
        
        open(isldprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_sld-'//chr//'.txt', status='replace')
        open(isldsat,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'sat_sld-'//chr//'.txt', status='replace')
        open(igasprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_gas-'//chr//'.txt', status='replace')
        open(iaqprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_aq-'//chr//'.txt', status='replace')
        open(ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
            & //'bsd-'//chr//'.txt', status='replace')
            
        write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(isldsat,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
        write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
        write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' time '

        do iz = 1, Nz
            write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
            write(isldsat,*) z(iz),(omega(isps,iz),isps = 1, nsp_sld),time
            write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
            write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
            write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
        end do
        irec=irec+1

        close(isldprof)
        close(isldsat)
        close(iaqprof)
        close(igasprof)
        close(ibsd)
        
        do isps=1,nsp_sld 
            open(isldflx(isps), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
                & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', action='write',status='old',access='append')
            write(isldflx(isps),*) time,(sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
            close(isldflx(isps))
        enddo 
        
        do ispa=1,nsp_aq 
            open(iaqflx(ispa), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
                & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', action='write',status='old',access='append')
            write(iaqflx(ispa),*) time,(sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
            close(iaqflx(ispa))
        enddo 
        
        do ispg=1,nsp_gas 
            open(igasflx(ispg), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
                & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', action='write',status='old',access='append')
            write(igasflx(ispg),*) time,(sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
            close(igasflx(ispg))
        enddo 

    end if

    it = it + 1
    time = time + dt
    count_dtunchanged = count_dtunchanged + 1
    
    progress_rate_prev = progress_rate
    
    call cpu_time(time_fin)
    
    ! progress_rate = dt/(time_fin-time_start)*sec2yr ! (model yr)/(computer yr)
    progress_rate = (time_fin-time_start) ! (computer sec)
    
    if (it/=1) then 
        if (flgreducedt) then 
            maxdt = maxdt/10d0
            flgreducedt = .false.
            count_dtunchanged = 0
        else
            ! maxdt = maxdt* (progress_rate/progress_rate_prev)**0.33d0
            maxdt = maxdt* (progress_rate/progress_rate_prev)**(-0.33d0)
            if (maxdt > 1d2) maxdt = 1d2
            if (dt < maxdt) count_dtunchanged = 0
            ! if (dt > maxdt) dt = maxdt
        endif 
        
        if (count_dtunchanged > 1000) then 
            maxdt = maxdt*10d0
            count_dtunchanged = 0
        endif 
    endif 
    
    print *,'progress_rate, maxdt, count_dtunchanged',progress_rate, maxdt, count_dtunchanged
    
    flgreducedt_prev = flgreducedt
    
end do

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

endsubroutine weathering_main

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_variables_num( &
    & nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! output
    & )
implicit none

integer,intent(out):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext
character(500) file_name

file_name = './slds.in'
call Console4(file_name,nsp_sld)
file_name = './solutes.in'
call Console4(file_name,nsp_aq)
file_name = './gases.in'
call Console4(file_name,nsp_gas)
file_name = './extrxns.in'
call Console4(file_name,nrxn_ext)

nsp_sld = nsp_sld - 1
nsp_aq = nsp_aq - 1
nsp_gas = nsp_gas - 1
nrxn_ext = nrxn_ext - 1

endsubroutine get_variables_num

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_variables( &
    & nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! input
    & ,chraq,chrgas,chrsld,chrrxn_ext &! output
    & )
implicit none

integer,intent(in):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext
character(5),dimension(nsp_sld),intent(out)::chrsld 
character(5),dimension(nsp_aq),intent(out)::chraq 
character(5),dimension(nsp_gas),intent(out)::chrgas 
character(5),dimension(nrxn_ext),intent(out)::chrrxn_ext 

character(500) file_name
integer ispa,ispg,isps,irxn

if (nsp_aq>=1) then 
    file_name = './solutes.in'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do ispa =1,nsp_aq
        read(50,*) chraq(ispa) 
    enddo 
    close(50)
endif 

if (nsp_sld>=1) then 
    file_name = './slds.in'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do isps =1,nsp_sld
        read(50,*) chrsld(isps) 
    enddo  
    close(50)
endif 

if (nsp_gas>=1) then 
    file_name = './gases.in'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do ispg =1,nsp_gas
        read(50,*) chrgas(ispg) 
    enddo 
    close(50)
endif 

if (nrxn_ext>=1) then 
    file_name = './extrxns.in'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do irxn =1,nrxn_ext
        read(50,*) chrrxn_ext(irxn) 
    enddo 
    close(50)
endif 

endsubroutine get_variables

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_bsdvalues( &
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name &! output
    & )
implicit none

integer,intent(out):: nz
real(kind=8),intent(out)::ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80
character(500),intent(out)::sim_name

character(500) file_name

file_name = './frame.in'
open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
read(50,*) ztot
read(50,*) nz
read(50,*) ttot
read(50,*) rainpowder
read(50,*) zsupp
read(50,*) poroi
read(50,*) satup
read(50,*) zsat
read(50,*) w
read(50,*) qin
read(50,*) p80
read(50,*) sim_name
close(50)

print*,'nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name'
print*,nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name

endsubroutine get_bsdvalues

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_rainwater( &
    & nsp_aq_all,chraq_all,def_rain &! input
    & ,rain_all &! output
    & )
implicit none

integer,intent(in):: nsp_aq_all
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
real(kind=8),dimension(nsp_aq_all),intent(out)::rain_all
real(kind=8),intent(in)::def_rain 
character(5) chr_tmp
real(kind=8) val_tmp

character(500) file_name
integer i,n_tmp

file_name = './rain.in'
call Console4(file_name,n_tmp)

n_tmp = n_tmp - 1

! in default 
rain_all = def_rain

if (n_tmp <= 0) return

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
do i =1,n_tmp
    read(50,*) chr_tmp,val_tmp
    if (any(chraq_all == chr_tmp)) then 
        rain_all(findloc(chraq_all,chr_tmp,dim=1)) = val_tmp
    endif 
enddo 
close(50)


endsubroutine get_rainwater

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_dust( &
    & nsp_sld_all,chrsld_all,def_dust &! input
    & ,dust_frct_all &! output
    & )
implicit none

integer,intent(in):: nsp_sld_all
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
real(kind=8),dimension(nsp_sld_all),intent(out)::dust_frct_all
real(kind=8),intent(in)::def_dust 
character(5) chr_tmp
real(kind=8) val_tmp

character(500) file_name
integer i,n_tmp

file_name = './dust.in'
call Console4(file_name,n_tmp)

n_tmp = n_tmp - 1

! in default 
dust_frct_all = def_dust

if (n_tmp <= 0) return

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
do i =1,n_tmp
    read(50,*) chr_tmp,val_tmp
    if (any(chrsld_all == chr_tmp)) then 
        dust_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
    endif 
enddo 
close(50)


endsubroutine get_dust

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_parentrock( &
    & nsp_sld_all,chrsld_all,def_pr &! input
    & ,parentrock_frct_all &! output
    & )
implicit none

integer,intent(in):: nsp_sld_all
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
real(kind=8),dimension(nsp_sld_all),intent(out)::parentrock_frct_all
real(kind=8),intent(in)::def_pr 
character(5) chr_tmp
real(kind=8) val_tmp

character(500) file_name
integer i,n_tmp

file_name = './parentrock.in'
call Console4(file_name,n_tmp)

n_tmp = n_tmp - 1

! in default 
parentrock_frct_all = def_pr

if (n_tmp <= 0) return

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
do i =1,n_tmp
    read(50,*) chr_tmp,val_tmp
    if (any(chrsld_all == chr_tmp)) then 
        parentrock_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
    endif 
enddo 
close(50)


endsubroutine get_parentrock

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_atm( &
    & nsp_gas_all,chrgas_all &! input
    & ,atm_all &! output
    & )
implicit none

integer,intent(in):: nsp_gas_all
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
real(kind=8),dimension(nsp_gas_all),intent(out)::atm_all
character(5) chr_tmp
real(kind=8) val_tmp

character(500) file_name
integer i,n_tmp

file_name = './atm.in'
call Console4(file_name,n_tmp)

n_tmp = n_tmp - 1

! in default 
atm_all(findloc(chrgas_all,'po2',dim=1)) = 0.21d0
atm_all(findloc(chrgas_all,'pco2',dim=1)) = 10d0**(-3.5d0)

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
do i =1,n_tmp
    read(50,*) chr_tmp,val_tmp
    if (any(chrgas_all == chr_tmp)) then 
        atm_all(findloc(chrgas_all,chr_tmp,dim=1)) = val_tmp
    endif 
enddo 
close(50)


endsubroutine get_atm

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine Console4(file_name,i)

implicit none

integer,intent(out) :: i
character(500),intent(in)::file_name

open(9, file =trim(adjustl(file_name)))

i = 0
do 
    read(9, *, end = 99)
    i = i + 1
end do 

! 99 print *, i
99 continue
close(9)

end subroutine Console4

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine coefs( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2,pco2th &! input
    & ,pro,dgaso,dgasc,daqo,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,ucv &! output
    & ,ksil,keqsil,kab,keqab,kfo,keqfo,kan,keqan,kcc,keqcc,koxs,koxa,koxs2,khco2 &! output
    & ,k1si,k2si,kcca,keqcca,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! output
    & ,k1al,k2al,k3al,k4al,keqka,kka,dal,keqgb,kgb &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0,pco2i,swoxa,pco2th
real(kind=8),dimension(nz),intent(in)::pro,pco2
real(kind=8),intent(out)::dgaso,dgasc,daqo,daqc,dfe2,dfe3,dso4,dna,dmg,dsi,dca,kho,kco2,k1,k2,kw,khco2i,keqsil &
    & ,keqab,keqfo,keqcc,keqan,ucv,k1si,k2si,keqcca,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka,dal,keqgb
real(kind=8),dimension(nz),intent(out)::ksil,kab,kfo,kan,kcc,koxs,koxa,koxs2,khco2,kcca,kka,kgb

real(kind=8),dimension(nz):: koxsi,koxai,koxs2i
real(kind=8) k_arrhenius
real(kind=8) :: cal2j = 4.184d0 


ucv = 1.0d0/(rg2*(tempk_0+tc))

dfe2 = 1.7016d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dfe3 = 1.5664d-2 ! m^2 yr^-1 ! at 15 C; Li and Gregory
dso4 = 2.54d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dna  = 3.19d-2   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dmg  = 0.017218079d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dsi  = 0.03689712d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dca  = 0.019023312d0   ! m^2 yr^-1 ! at 15 C; Li and Gregory 
dal = 0.011656226d0    ! m^2 yr^-1 ! at 15 C; Li and Gregory 

dgaso = 6.09d2 ! m^2 yr^-1
daqo = 5.49d-2 ! m^2 yr^-1
dgasc = 441.504d0 ! m^2 yr^-1 (Assuming 0.14 cm2/sec)
daqc = 0.022459852d0 ! m^2 yr^-1 (for C32- from Li and Gregory 1974)

dgaso = dgaso*exp(-4.18d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
daqo = daqo*exp(-20.07d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dfe2=dfe2*exp(-19.615251d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dfe3=dfe3*exp(-14.33659d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dso4=dso4*exp(-20.67364d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dna=dna*exp(-20.58566d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dmg=dmg*exp(-18.51979d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)
dal=dal*exp(-21.27788d0*(1.0d0/(tempk_0+tc)-1.0d0/(tempk_0+15.0d0))/rg)

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

kka = k_arrhenius(10d0**(-13.18d0)*sec2yr,25d0+tempk_0,tc+tempk_0,22.2d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
! kaolinite dissolution: Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
keqka = 8.310989613d0 ! gcw
keqka = k_arrhenius(10d0**(7.435d0),25d0+tempk_0,tc+tempk_0,-35.3d0*cal2j,rg) ! from PHREEQC.DAT  

! Al3+ + H2O = Al(OH)2+ + H+
k1al = k_arrhenius(10d0**(-5d0),25d0+tempk_0,tc+tempk_0,11.49d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 2H2O = Al(OH)2+ + 2H+
k2al = k_arrhenius(10d0**(-10.1d0),25d0+tempk_0,tc+tempk_0,26.90d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 3H2O = Al(OH)3 + 3H+
k3al = k_arrhenius(10d0**(-16.9d0),25d0+tempk_0,tc+tempk_0,39.89d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 4H2O = Al(OH)4- + 4H+
k4al = k_arrhenius(10d0**(-22.7d0),25d0+tempk_0,tc+tempk_0,42.30d0*cal2j,rg) ! from PHREEQC.DAT 

ksil = 1.31d-9*1d4  ! mol/m2/yr  ! from Li et al., 2014

keqsil = 3.412182823d0 - 0.5d0* 8.310989613d0   ! albite + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)  
keqsil = 10.0d0**(keqsil)

kab = ksil
kab = k_arrhenius(10d0**(-12.56d0)*sec2yr,25d0+tempk_0,tc+tempk_0,69.8d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
keqab = keqsil
! NaAlSi3O8 + 8 H2O = Na+ + Al(OH)4- + 3 H4SiO4
keqab = k_arrhenius(10d0**(-18.002d0),25d0+tempk_0,tc+tempk_0,25.896d0*cal2j,rg) ! from PHREEQC.DAT 
! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
keqab = 10d0**3.412182823d0

kfo = -10.64d0  ! mol/m2/sec ! from Beering et al 2020 (only neutral weathering)
kfo = 10d0**(kfo)*60d0*60d0*24d0*365d0 ! mol/m2/yr 
! kfo = ksil

! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O
keqfo= 27.8626d0  ! Sugimori et al. (2012) 
keqfo = 10d0**keqfo
! -208.5932252 

! kan = 10d0**(-9.12d0)*60d0*60d0*24d0*365d0 &! mol/m2/yr  
    ! & *exp(-17.8d0/rg*(1d0/(tc+273d0)-1d0/(25d0+273d0)))
kan = k_arrhenius(10d0**(-9.12d0)*sec2yr,25d0+tempk_0,tc+tempk_0,17.8d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)

! anorthite (CaAl2Si2O8) + 2H+ + H2O --> kaolinite(Al2Si2O5(OH)4) + Ca2+ 
keqan = 28.8615308d0 - 8.310989613d0
keqan = 10d0**keqan
! CaAl2Si2O8 + 8 H2O = Ca+2 + 2 Al(OH)4- + 2 H4SiO4
keqan = k_arrhenius(10d0**(-19.714d0),25d0+tempk_0,tc+tempk_0,11.580d0*cal2j,rg) ! from PHREEQC.DAT 
! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
keqan = 10d0**28.8615308d0

kcc = k_arrhenius(10d0**(-5.81d0)*sec2yr,25d0+tempk_0,tc+tempk_0,23.5d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
! kcc = kcc**merge(0.0d0,1.0d0,pco2<pco2th)
! kcc = 0d0

keqcc = 10d0**(-8.43d0) ! Kanzaki and Murakami 2015

kcca = kcc*1d-4 ! assuming 1e4 times slower kinetic but units are different (mol-1 m3) 
keqcca = keqcc*20d0 ! assuming more harder to precipitate

kgb = k_arrhenius(10d0**(-11.50d0)*sec2yr,25d0+tempk_0,tc+tempk_0,61.2d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
! Al(OH)3 + 3 H+ = Al+3 + 3 H2O
keqgb = k_arrhenius(10d0**(8.11d0),25d0+tempk_0,tc+tempk_0,-22.80d0*cal2j,rg) ! from PHREEQC.DAT 

k1si = k_arrhenius(10d0**(-9.83d0),25d0+tempk_0,tc+tempk_0,6.12d0*cal2j,rg) ! from PHREEQC.DAT 
k2si = k_arrhenius(10d0**(-23d0),25d0+tempk_0,tc+tempk_0,17.6d0*cal2j,rg) ! from PHREEQC.DAT 

! Mg2+ + H2O = Mg(OH)+ + H+
k1mg = k_arrhenius(10d0**(-11.44d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
! Mg2+ + CO32- = MgCO3 
k1mgco3 = k_arrhenius(10d0**(2.98d0),25d0+tempk_0,tc+tempk_0,2.713d0*cal2j,rg) ! from PHREEQC.DAT 
! Mg2+ + H+ + CO32- = MgHCO3
k1mghco3 = k_arrhenius(10d0**(11.399d0),25d0+tempk_0,tc+tempk_0,-2.771d0*cal2j,rg) ! from PHREEQC.DAT 

! Ca2+ + H2O = Ca(OH)+ + H+
k1ca = k_arrhenius(10d0**(-12.78d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
! (No delta_h is reported so used the same value for Mg)
! Ca2+ + CO32- = CaCO3 
k1caco3 = k_arrhenius(10d0**(3.224d0),25d0+tempk_0,tc+tempk_0,3.545d0*cal2j,rg) ! from PHREEQC.DAT 
! Ca2+ + H+ + CO32- = CaHCO3
k1cahco3 = k_arrhenius(10d0**(11.435d0),25d0+tempk_0,tc+tempk_0,-0.871d0*cal2j,rg) ! from PHREEQC.DAT 

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

subroutine coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0
real(kind=8),dimension(nz),intent(in)::pro
real(kind=8),dimension(nz)::oh
real(kind=8),intent(out)::ucv,kw

real(kind=8) k_arrhenius
real(kind=8) :: cal2j = 4.184d0 

integer,intent(in)::nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
real(kind=8),dimension(nsp_aq_all),intent(out)::daq_all
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all
real(kind=8),dimension(nsp_gas_all),intent(out)::dgasa_all,dgasg_all
real(kind=8),dimension(nsp_gas_all,3),intent(out)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(out)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_c
real(kind=8),dimension(nsp_sld_all,nz),intent(out)::ksld_all
real(kind=8),dimension(nsp_sld_all),intent(out)::keqsld_all
real(kind=8),dimension(nrxn_ext_all,nz),intent(out)::krxn1_ext_all
real(kind=8),dimension(nrxn_ext_all,nz),intent(out)::krxn2_ext_all

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

ucv = 1.0d0/(rg2*(tempk_0+tc))

! Aq species diffusion from Li and Gregory 1974 except for Si which is based on Rebreanu et al. 2008
daq_all(findloc(chraq_all,'fe2',dim=1))= k_arrhenius(1.7016d-2    , 15d0+tempk_0, tc+tempk_0, 19.615251d0, rg)
daq_all(findloc(chraq_all,'fe3',dim=1))= k_arrhenius(1.5664d-2    , 15d0+tempk_0, tc+tempk_0, 14.33659d0 , rg)
daq_all(findloc(chraq_all,'so4',dim=1))= k_arrhenius(2.54d-2      , 15d0+tempk_0, tc+tempk_0, 20.67364d0 , rg)
daq_all(findloc(chraq_all,'na',dim=1)) = k_arrhenius(3.19d-2      , 15d0+tempk_0, tc+tempk_0, 20.58566d0 , rg)
daq_all(findloc(chraq_all,'k',dim=1))  = k_arrhenius(4.8022699d-2 , 15d0+tempk_0, tc+tempk_0, 18.71816d0 , rg)
daq_all(findloc(chraq_all,'mg',dim=1)) = k_arrhenius(1.7218079d-2 , 15d0+tempk_0, tc+tempk_0, 18.51979d0 , rg)
daq_all(findloc(chraq_all,'si',dim=1)) = k_arrhenius(2.682396d-2  , 15d0+tempk_0, tc+tempk_0, 22.71378d0 , rg)
daq_all(findloc(chraq_all,'ca',dim=1)) = k_arrhenius(1.9023312d-2 , 15d0+tempk_0, tc+tempk_0, 20.219661d0, rg)
daq_all(findloc(chraq_all,'al',dim=1)) = k_arrhenius(1.1656226d-2 , 15d0+tempk_0, tc+tempk_0, 21.27788d0 , rg)

! values used in Kanzaki and Murakami 2016 for oxygen 
dgasa_all(findloc(chrgas_all,'po2',dim=1)) = k_arrhenius(5.49d-2 , 15d0+tempk_0, tc+tempk_0, 20.07d0 , rg)
dgasg_all(findloc(chrgas_all,'po2',dim=1)) = k_arrhenius(6.09d2  , 15d0+tempk_0, tc+tempk_0, 4.18d0  , rg)

! assuming a value of 0.14 cm2/sec and O2 gas activation energy for CO2 gas 
! and CO32- diffusion from Li and Greogy 1974 for aq CO2 
dgasa_all(findloc(chrgas_all,'pco2',dim=1)) = k_arrhenius(2.2459852d-2, 15d0+tempk_0, tc+tempk_0, 21.00564d0, rg)
dgasg_all(findloc(chrgas_all,'pco2',dim=1)) = k_arrhenius(441.504d0   , 15d0+tempk_0, tc+tempk_0, 4.18d0    , rg)

kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc**2d0+0.000000555d0*tc**3d0-0.0000000007581d0*tc**4d0  ! Murakami et al. 2011
kw = k_arrhenius(10d0**(-14.35d0), tempk_0+15.0d0, tempk_0+tc, 58.736742d0, rg) ! from Kanzaki and Murakami 2015

oh = kw/pro


keqgas_h = 0d0

! kho = k_arrhenius(10.0d0**(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg)
keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0) = &
    & k_arrhenius(10d0**(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg)

keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0) = &
    & k_arrhenius(10d0**(-1.34d0), tempk_0+15.0d0, tempk_0+tc, -21.33183d0, rg) ! from Kanzaki and Murakami 2015
keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1) = &
    & k_arrhenius(10d0**(-6.42d0), tempk_0+15.0d0, tempk_0+tc, 11.94453d0, rg) ! from Kanzaki and Murakami 2015
keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2) = &
    & k_arrhenius(10d0**(-10.43d0), tempk_0+15.0d0, tempk_0+tc, 17.00089d0, rg) ! from Kanzaki and Murakami 2015


keqaq_c = 0d0
keqaq_h = 0d0

! Al3+ + H2O = Al(OH)2+ + H+
keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h1) = &
    & k_arrhenius(10d0**(-5d0),25d0+tempk_0,tc+tempk_0,11.49d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 2H2O = Al(OH)2+ + 2H+
keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h2) = &
    & k_arrhenius(10d0**(-10.1d0),25d0+tempk_0,tc+tempk_0,26.90d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 3H2O = Al(OH)3 + 3H+
keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h3) = &
    & k_arrhenius(10d0**(-16.9d0),25d0+tempk_0,tc+tempk_0,39.89d0*cal2j,rg) ! from PHREEQC.DAT 
! Al3+ + 4H2O = Al(OH)4- + 4H+
keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h4) = &
    & k_arrhenius(10d0**(-22.7d0),25d0+tempk_0,tc+tempk_0,42.30d0*cal2j,rg) ! from PHREEQC.DAT 


keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h1) = &
    & k_arrhenius(10d0**(-9.83d0),25d0+tempk_0,tc+tempk_0,6.12d0*cal2j,rg) ! from PHREEQC.DAT 
keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h2) = &
    & k_arrhenius(10d0**(-23d0),25d0+tempk_0,tc+tempk_0,17.6d0*cal2j,rg) ! from PHREEQC.DAT 


! Mg2+ + H2O = Mg(OH)+ + H+
keqaq_h(findloc(chraq_all,'mg',dim=1),ieqaq_h1) = &
    & k_arrhenius(10d0**(-11.44d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
! Mg2+ + CO32- = MgCO3 
keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_co3) = &
    & k_arrhenius(10d0**(2.98d0),25d0+tempk_0,tc+tempk_0,2.713d0*cal2j,rg) ! from PHREEQC.DAT 
! Mg2+ + H+ + CO32- = MgHCO3
keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_hco3) = & 
    & k_arrhenius(10d0**(11.399d0),25d0+tempk_0,tc+tempk_0,-2.771d0*cal2j,rg) ! from PHREEQC.DAT 

! Ca2+ + H2O = Ca(OH)+ + H+
keqaq_h(findloc(chraq_all,'ca',dim=1),ieqaq_h1) =  &
    & k_arrhenius(10d0**(-12.78d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ! from PHREEQC.DAT 
! (No delta_h is reported so used the same value for Mg)
! Ca2+ + CO32- = CaCO3 
keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_co3) = &
    & k_arrhenius(10d0**(3.224d0),25d0+tempk_0,tc+tempk_0,3.545d0*cal2j,rg) ! from PHREEQC.DAT 
! Ca2+ + H+ + CO32- = CaHCO3
keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_hco3) = &
    & k_arrhenius(10d0**(11.435d0),25d0+tempk_0,tc+tempk_0,-0.871d0*cal2j,rg) ! from PHREEQC.DAT 
    
    
! Fe2+ + H2O = Fe(OH)+ + H+
keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1) = &
    & k_arrhenius(10d0**(-9.51d0),25d0+tempk_0,tc+tempk_0, 40.3d0,rg) ! from Kanzaki and Murakami 2016
! Fe2+ + CO32- = FeCO3 
keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3) = &
    & k_arrhenius(10d0**(5.69d0),25d0+tempk_0,tc+tempk_0, -45.6d0,rg) ! from Kanzaki and Murakami 2016
! Fe2+ + H+ + CO32- = FeHCO3
keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3) = &
    & k_arrhenius(10d0**(1.47d0),25d0+tempk_0,tc+tempk_0, -18d0,rg) &! from Kanzaki and Murakami 2016 
    & /keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2) 


! Fe3+ + H2O = Fe(OH)2+ + H+
keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h1) = &
    & k_arrhenius(10d0**(-2.19d0),25d0+tempk_0,tc+tempk_0,10.4d0*cal2j,rg) ! from PHREEQC.DAT 
! Fe3+ + 2H2O = Fe(OH)2+ + 2H+
keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h2) = &
    & k_arrhenius(10d0**(-5.67d0),25d0+tempk_0,tc+tempk_0,17.1d0*cal2j,rg) ! from PHREEQC.DAT 
! Fe3+ + 3H2O = Fe(OH)3 + 3H+
keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h3) = &
    & k_arrhenius(10d0**(-12.56d0),25d0+tempk_0,tc+tempk_0,24.8d0*cal2j,rg) ! from PHREEQC.DAT 
! Fe3+ + 4H2O = Fe(OH)4- + 4H+
keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h4) = &
    & k_arrhenius(10d0**(-21.6d0),25d0+tempk_0,tc+tempk_0,31.9d0*cal2j,rg) ! from PHREEQC.DAT 




!!! ----------- Solid phases ------------------------!!
ksld_all = 0d0 
keqsld_all = 0d0

ksld_all(findloc(chrsld_all,'ka',dim=1),:) = &
    & k_arrhenius(10d0**(-13.18d0)*sec2yr,25d0+tempk_0,tc+tempk_0,22.2d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**0.777d0*k_arrhenius(10d0**(-11.31d0)*sec2yr,25d0+tempk_0,tc+tempk_0,65.9d0,rg) &!(acid weathering from Palandri and Kharaka, 2004)
    & + pro**(-0.472d0)*k_arrhenius(10d0**(-17.05d0)*sec2yr,25d0+tempk_0,tc+tempk_0,17.9d0,rg) !(alkarine weathering from Palandri and Kharaka, 2004)
! kaolinite dissolution: Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
! keqka = 8.310989613d0 ! gcw
keqsld_all(findloc(chrsld_all,'ka',dim=1)) = &
    & k_arrhenius(10d0**(7.435d0),25d0+tempk_0,tc+tempk_0,-35.3d0*cal2j,rg) ! from PHREEQC.DAT  


ksld_all(findloc(chrsld_all,'ab',dim=1),:) = &
    & k_arrhenius(10d0**(-12.56d0)*sec2yr,25d0+tempk_0,tc+tempk_0,69.8d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**0.457d0*k_arrhenius(10d0**(-10.16d0)*sec2yr,25d0+tempk_0,tc+tempk_0,65d0,rg) &!(acid weathering from Palandri and Kharaka, 2004)
    & + pro**(-0.572d0)*k_arrhenius(10d0**(-15.6d0)*sec2yr,25d0+tempk_0,tc+tempk_0,71d0,rg) !(alkarine weathering from Palandri and Kharaka, 2004)

! NaAlSi3O8 + 8 H2O = Na+ + Al(OH)4- + 3 H4SiO4
keqsld_all(findloc(chrsld_all,'ab',dim=1)) = &
    & k_arrhenius(10d0**(-18.002d0),25d0+tempk_0,tc+tempk_0,25.896d0*cal2j,rg) ! from PHREEQC.DAT 
! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
keqsld_all(findloc(chrsld_all,'ab',dim=1)) = &
    & k_arrhenius(10d0**3.412182823d0,15d0+tempk_0,tc+tempk_0,-54.15042876d0,rg)   ! Kanzaki and Murakami 2018



ksld_all(findloc(chrsld_all,'kfs',dim=1),:) = &
    & k_arrhenius(10d0**(-12.41d0)*sec2yr,25d0+tempk_0,tc+tempk_0,9.08*cal2j,rg) &!(only neutral weathering from Brantley et al 2008)
    & + pro**0.5d0*k_arrhenius(10d0**(-10.06d0)*sec2yr,25d0+tempk_0,tc+tempk_0,12.4d0*cal2j,rg) &!(acid weathering from Brantley et al 2008)
    & + oh**0.823d0*k_arrhenius(10d0**(-9.68d0)*sec2yr,25d0+tempk_0,tc+tempk_0,22.5d0*cal2j,rg) !(alkarine weathering from Brantley et al 2008)
! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
keqsld_all(findloc(chrsld_all,'kfs',dim=1)) = &
    & k_arrhenius(10d0**0.227294204d0,15d0+tempk_0,tc+tempk_0,-26.30862098d0,rg)   ! Kanzaki and Murakami 2018


! ksld_all(findloc(chrsld_all,'fo',dim=1),:) = &
    ! & k_arrhenius(10d0**(-10.64d0)*sec2yr,25d0+tempk_0,tc+tempk_0,79d0,rg)  ! mol/m2/yr  from Beering et al 2020 (only neutral weathering)
ksld_all(findloc(chrsld_all,'fo',dim=1),:) = &
    & k_arrhenius(10d0**(-10.64d0)*sec2yr,25d0+tempk_0,tc+tempk_0,79d0,rg)  &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**0.47d0*k_arrhenius(10d0**(-6.85d0)*sec2yr,25d0+tempk_0,tc+tempk_0,67.2d0,rg)  !(acid weathering from Palandri and Kharaka, 2004)

! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O
! keqfo= 27.8626d0  ! Sugimori et al. (2012) 
! keqfo = 10d0**keqfo
keqsld_all(findloc(chrsld_all,'fo',dim=1)) = &
    & k_arrhenius(10d0**29.41364324d0,15d0+tempk_0,tc+tempk_0,-208.5932252d0,rg)   ! Kanzaki and Murakami 2018

! -208.5932252 


ksld_all(findloc(chrsld_all,'fa',dim=1),:) = &
    & k_arrhenius(10d0**(-12.80d0)*sec2yr,25d0+tempk_0,tc+tempk_0, 94.4d0, rg)  &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro*k_arrhenius(10d0**(-4.80d0)*sec2yr,25d0+tempk_0,tc+tempk_0, 94.4d0, rg)  !(acid weathering from Palandri and Kharaka, 2004)

! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O
keqsld_all(findloc(chrsld_all,'fa',dim=1)) = &
    & k_arrhenius(10d0**19.98781342d0,15d0+tempk_0,tc+tempk_0,-153.7676621d0,rg)   ! Kanzaki and Murakami 2018



! kan = 10d0**(-9.12d0)*60d0*60d0*24d0*365d0 &! mol/m2/yr  
    ! & *exp(-17.8d0/rg*(1d0/(tc+273d0)-1d0/(25d0+273d0)))
ksld_all(findloc(chrsld_all,'an',dim=1),:) = & 
    & k_arrhenius(10d0**(-9.12d0)*sec2yr,25d0+tempk_0,tc+tempk_0,17.8d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**1.411d0*k_arrhenius(10d0**(-3.5d0)*sec2yr,25d0+tempk_0,tc+tempk_0,16.6d0,rg) !(acid weathering from Palandri and Kharaka, 2004)

! anorthite (CaAl2Si2O8) + 2H+ + H2O --> kaolinite(Al2Si2O5(OH)4) + Ca2+ 
! keqan = 28.8615308d0 - 8.310989613d0
! keqan = 10d0**keqan
! CaAl2Si2O8 + 8 H2O = Ca+2 + 2 Al(OH)4- + 2 H4SiO4
keqsld_all(findloc(chrsld_all,'an',dim=1)) = & 
    & k_arrhenius(10d0**(-19.714d0),25d0+tempk_0,tc+tempk_0,11.580d0*cal2j,rg) ! from PHREEQC.DAT 
! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
keqsld_all(findloc(chrsld_all,'an',dim=1)) = & 
    & k_arrhenius(10d0**28.8615308d0,15d0+tempk_0,tc+tempk_0,-292.8769275d0,rg)   ! Kanzaki and Murakami 2018


ksld_all(findloc(chrsld_all,'cc',dim=1),:) = & 
    & k_arrhenius(10d0**(-5.81d0)*sec2yr,25d0+tempk_0,tc+tempk_0,23.5d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro*k_arrhenius(10d0**(-0.3d0)*sec2yr,25d0+tempk_0,tc+tempk_0,14.4d0,rg) !(acid weathering from Palandri and Kharaka, 2004)
! kcc = kcc**merge(0.0d0,1.0d0,pco2<pco2th)
! kcc = 0d0

keqsld_all(findloc(chrsld_all,'cc',dim=1)) = & 
    & k_arrhenius(10d0**(-8.43d0),15d0+tempk_0,tc+tempk_0,-8.028943471d0,rg) ! Kanzaki and Murakami 2015



ksld_all(findloc(chrsld_all,'gb',dim=1),:) = & 
    & k_arrhenius(10d0**(-11.50d0)*sec2yr,25d0+tempk_0,tc+tempk_0,61.2d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**0.992d0*k_arrhenius(10d0**(-7.65d0)*sec2yr,25d0+tempk_0,tc+tempk_0,47.5d0,rg) &!(acid weathering from Palandri and Kharaka, 2004)
    & + pro**(-0.784d0)*k_arrhenius(10d0**(-16.65d0)*sec2yr,25d0+tempk_0,tc+tempk_0,80.1d0,rg) !(alkarine weathering from Palandri and Kharaka, 2004)
! Al(OH)3 + 3 H+ = Al+3 + 3 H2O
keqsld_all(findloc(chrsld_all,'gb',dim=1)) = &
    & k_arrhenius(10d0**(8.11d0),25d0+tempk_0,tc+tempk_0,-22.80d0*cal2j,rg) ! from PHREEQC.DAT 


ksld_all(findloc(chrsld_all,'gt',dim=1),:) = & 
    & k_arrhenius(10d0**(-7.94d0)*sec2yr,25d0+tempk_0,tc+tempk_0,86.5d0,rg) !(only neutral weathering from Palandri and Kharaka, 2004)
! Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
keqsld_all(findloc(chrsld_all,'gt',dim=1)) = &
    & k_arrhenius(10d0**(0.5345d0),25d0+tempk_0,tc+tempk_0,-61.53703d0,rg) ! from Sugimori et al. 2012 


ksld_all(findloc(chrsld_all,'ct',dim=1),:) = & 
    & k_arrhenius(10d0**(-12d0)*sec2yr,25d0+tempk_0,tc+tempk_0,73.5d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**(-0.23d0)*k_arrhenius(10d0**(-13.58d0)*sec2yr,25d0+tempk_0,tc+tempk_0,73.5d0,rg) !(alkarine weathering from Palandri and Kharaka, 2004)
! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
keqsld_all(findloc(chrsld_all,'ct',dim=1)) = &
    & k_arrhenius(10d0**(32.2),25d0+tempk_0,tc+tempk_0,-46.800d0*cal2j,rg) ! from PHREEQC.DAT 


ksld_all(findloc(chrsld_all,'cabd',dim=1),:) = & 
    & k_arrhenius(10d0**(-12.78d0)*sec2yr,25d0+tempk_0,tc+tempk_0,35d0,rg) &!(only neutral weathering for smectite from Palandri and Kharaka, 2004)
    & + pro**0.34d0*k_arrhenius(10d0**(-10.98d0)*sec2yr,25d0+tempk_0,tc+tempk_0,23.6d0,rg) &!(only neutral weathering for smectite from Palandri and Kharaka, 2004)
    & + pro**(-0.4d0)*k_arrhenius(10d0**(-16.52d0)*sec2yr,25d0+tempk_0,tc+tempk_0,58.9d0,rg) !(only neutral weathering for smectite from Palandri and Kharaka, 2004)
! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
keqsld_all(findloc(chrsld_all,'cabd',dim=1)) = &
    & k_arrhenius(10d0**(7.269946518d0),15d0+tempk_0,tc+tempk_0,-157.0186168d0,rg) ! from Kanzaki & Murakami 2018


ksld_all(findloc(chrsld_all,'dp',dim=1),:) = & 
    & k_arrhenius(10d0**(-11.11d0)*sec2yr,25d0+tempk_0,tc+tempk_0,50.6d0,rg) &!(only neutral weathering from Palandri and Kharaka, 2004)
    & + pro**0.71d0*k_arrhenius(10d0**(-6.36d0)*sec2yr,25d0+tempk_0,tc+tempk_0,96.1d0,rg) !(acid weathering from Palandri and Kharaka, 2004)
! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
keqsld_all(findloc(chrsld_all,'dp',dim=1)) = &
    & k_arrhenius(10d0**(21.79853309d0),15d0+tempk_0,tc+tempk_0,-138.6020832d0,rg) ! from Kanzaki & Murakami 2018


ksld_all(findloc(chrsld_all,'hb',dim=1),:) = & 
    & k_arrhenius(10d0**(-11.97d0)*sec2yr,25d0+tempk_0,tc+tempk_0,78.0d0,rg) &!(only neutral weathering for augite from Palandri and Kharaka, 2004)
    & + pro**0.70d0*k_arrhenius(10d0**(-6.82d0)*sec2yr,25d0+tempk_0,tc+tempk_0,78.0d0,rg) !(acid weathering for augite from Palandri and Kharaka, 2004)
! Hedenbergite  + 4 H+  = 2 H2O  + 2 SiO2(aq)  + Fe++  + Ca++
keqsld_all(findloc(chrsld_all,'hb',dim=1)) = &
    & k_arrhenius(10d0**(20.20981116d0),15d0+tempk_0,tc+tempk_0,-128.5d0,rg) ! from Kanzaki & Murakami 2018




ksld_all(findloc(chrsld_all,'py',dim=1),:) = & 
    & k_arrhenius(10.0d0**(-8.19d0)*sec2yr,15d0+tempk_0,tc+tempk_0,57d0,rg)  !!! excluding po2 and ph dependence
    ! & *(kho)**(0.50d0)/(pro**0.11d0) ! mol m^-2 yr^-1, Williamson and Rimstidt (1994)


!--------- other reactions -------------! 
krxn1_ext_all = 0d0
krxn2_ext_all = 0d0

krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) = & 
    & max(8.0d13*60.0d0*24.0d0*365.0d0*(kw/pro)**2.0d0, 1d-7*60.0d0*24.0d0*365.0d0)    
    ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)excluding the term (c*po2)
     
krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.71d0 ! vmax mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = &
    & krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) *1d1 ! reducing a bit to be fitted with modern soil pco2

krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)

endsubroutine coefs_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_po2( &
    & nz,po2th,dt,ucv,kho,dz,dgaso,daqo,po2i,poro,sat,po2,torg,tora,v &! input 
    & ,po2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::po2th,dt,ucv,kho,dz,dgaso,daqo,po2i
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
            & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
            & *(po2tmp+po2(iz-1)-2.0d0*po2(iz))/(dz**2.0d0) &
            & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
            & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
            & *(po2(iz)-po2(iz-1))/(dz**2.0d0) &
            & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2(iz)-po2(iz-1))/dz &
            & ) &
            & )
    else ! iz == 1
        po2x(iz) = max(0.0d0 &
            & ,-(dt/(ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3))* &
            & ((ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3) &
            & *(-po2(iz))/dt-(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
            & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)*(-2.0d0*po2(iz) + po2(iz+1)+po2i) &
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
    & nz,po2th,dt,ucv,kho,dz,dgaso,daqo,po2i,poro,sat,po2,torg,tora,v &! input 
    & ,po2x &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::po2th,dt,ucv,kho,dgaso,daqo,po2i
real(kind=8),dimension(nz),intent(in)::poro,sat,po2,torg,tora,v,dz
real(kind=8),dimension(nz),intent(out)::po2x

integer iz
real(kind=8) po2tmp,edifi,ediftmp
real(kind=8),dimension(nz)::alpha,edif

alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3
edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgaso +poro*sat*kho*1d3*tora*daqo
edifi = edif(1)
edifi = ucv*1d3*dgaso 

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

subroutine precalc_gases( &
    & nz,dt,ucv,dz,poro,sat,torg,tora,v,prox &! input 
    & ,nsp_gas,nsp_gas_all,chrgas,chrgas_all,keqgas_h,mgasi,mgasth,mgas &! input
    & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg &! input
    & ,mgasx &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::dt,ucv
real(kind=8)::pco2th,pco2i,kco2,k1,k2,kho
real(kind=8),dimension(nz),intent(in)::poro,sat,torg,tora,v,dz,prox
real(kind=8),dimension(nz)::khco2,pco2,resp,pco2x

integer iz,ispg
real(kind=8) pco2tmp,edifi,ediftmp
real(kind=8),dimension(nz)::alpha,edif

integer,intent(in)::nsp_gas,nsp_gas_all,nrxn_ext
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
character(5),dimension(nrxn_ext),intent(in)::chrrxn_ext
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_gas),intent(in)::mgasi,mgasth,dgasa,dgasg
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas
real(kind=8),dimension(nrxn_ext,nz),intent(in)::rxnext
real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx



integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

if (nsp_gas == 0) return

kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)

khco2 = kco2*(1d0+k1/prox + k1*k2/prox/prox)

kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

do ispg = 1, nsp_gas

    select case(trim(adjustl(chrgas(ispg))))
        case('pco2')
            alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2*1d3
            if (any(chrrxn_ext=='resp')) then 
                resp = rxnext(findloc(chrrxn_ext,'resp',dim=1),:)
            else 
                resp = 0d0
            endif 
        case('po2')
            alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3
            resp = 0d0
    endselect 
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg) +poro*sat*khco2*1d3*tora*dgasa(ispg)
    edifi = edif(1)
    edifi = ucv*1d3*dgasg(ispg) 

    pco2x = mgasx(ispg,:)
    pco2 = mgas(ispg,:)
    
    pco2i = mgasi(ispg)
    pco2th = mgasth(ispg)
    
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
    
    mgasx(ispg,:) = pco2x
    
enddo

endsubroutine precalc_gases

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds( &
    & nz,msth,dt,w,dz,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth   &! input
    & ,ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp,kcc,omega_cc,mvcc &! input
    & ,poro,hr,kcca,omega_cca,authig,sat,kka,mkai,mkath,omega_ka,mvka,mkasupp,mka &! input
    & ,msx,msilx,mfox,mabx,manx,mccx,mkax &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::msth,dt,w,msili,msi,mfoi,mabi,mani,mcci,msilth,mabth,manth,mfoth,mccth,mvcc,authig &
    & ,mkai,mkath,mvka
real(kind=8),dimension(nz),intent(in)::ms,msil,msilsupp,mfo,mfosupp,mab,mabsupp,mansupp,man,mcc,mccsupp &
    & ,dz,kcc,omega_cc,poro,hr,kcca,omega_cca,sat,kka,omega_ka,mkasupp,mka
real(kind=8),dimension(nz),intent(out)::msx,msilx,mfox,mabx,manx,mccx,mkax

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
enddo 

do iz = 1, nz

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
enddo

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
enddo

do iz = 1, nz
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
enddo

do iz = 1, nz
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
enddo

do iz = 1, nz
    if (mccx(iz)>=mccth) cycle

    if (iz/=nz) then 
        mccx(iz) = max(0d0, &
            & mcc(iz) +dt*(w*(mcc(iz+1)-mcc(iz))/dz(iz) + mccsupp(iz) &
            & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(omega_cc(iz) - 1d0) &
            & *merge(1d0,0d0,omega_cc(iz) - 1d0 > 0d0) &
            ! & *merge(1d0,0d0,mccx(iz) > mccth) &
            & + kcca(iz)*poro(iz)*sat(iz)*(omega_cca(iz)-1d0)*authig &
            & *merge(1d0,0d0,omega_cca(iz)-1d0 > 0d0) &
            & ) &
            & )
    else 
        mccx(iz) = max(0d0, &
            & mcc(iz) + dt*(w*(mcci-mcc(iz))/dz(iz)+ mccsupp(iz)&
            & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(omega_cc(iz) - 1d0) &
            & *merge(1d0,0d0,omega_cc(iz) - 1d0 > 0d0) &
            ! & *merge(1d0,0d0,mccx(iz) > mccth) &
            & + kcca(iz)*poro(iz)*sat(iz)*(omega_cca(iz)-1d0)*authig &
            & *merge(1d0,0d0,omega_cca(iz)-1d0 > 0d0) &
            & ) &
            & )
    endif 

enddo

do iz = 1, nz
    if (mkax(iz)>=mkath) cycle

    if (iz/=nz) then 
        mkax(iz) = max(0d0, &
            & mka(iz) +dt*(w*(mka(iz+1)-mka(iz))/dz(iz) + mkasupp(iz) &
            & +kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(omega_ka(iz) - 1d0) &
            & *merge(1d0,0d0,omega_ka(iz) - 1d0 > 0d0) &
            ! & *merge(1d0,0d0,mkax(iz) > mkath) &
            & ) &
            & )
    else 
        mkax(iz) = max(0d0, &
            & mka(iz) + dt*(w*(mkai-mka(iz))/dz(iz)+ mkasupp(iz)&
            & +kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(omega_ka(iz) - 1d0) &
            & *merge(1d0,0d0,omega_ka(iz) - 1d0 > 0d0) &
            ! & *merge(1d0,0d0,mkax(iz) > mkath) &
            & ) &
            & )
    endif 

enddo
    
! print *,mccx

endsubroutine precalc_slds

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds_v2( &
    & nz,dt,w,dz,poro,hr,sat &! input
    & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
    & ,msldx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w
real(kind=8)::mfoi,mfoth
real(kind=8),dimension(nz),intent(in)::dz,poro,hr,sat
real(kind=8),dimension(nz)::mfo,mfosupp,mfox

integer iz,isps

integer,intent(in)::nsp_sld,nsp_sld_2
character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_2),intent(in)::chrsld_2
real(kind=8),dimension(nsp_sld),intent(in)::msldth,msldi,mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp,ksld,omega
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx

if (nsp_sld == 0) return

do isps = 1,nsp_sld

    mfox = msldx(isps,:)
    mfo = msld(isps,:)
    mfoth = msldth(isps)
    mfoi = msldi(isps)
    mfosupp = msldsupp(isps,:)
    
    if (any(chrsld_2 ==chrsld(isps))) then 
        do iz = 1, nz
            if (mfox(iz)>=mfoth) cycle

            if (iz/=nz) then 
                mfox(iz) = max(0d0, &
                    & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & )
            else 
                mfox(iz) = max(0d0, &
                    & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & )
            endif 
        enddo
    else
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
        enddo
    endif 
    
    msldx(isps,:) = mfox
    
enddo 

endsubroutine precalc_slds_v2

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
    & ,kab,kan,kcc,kfo,hr,mvab,mvan,mvfo,mvcc,mabx,manx,mfox,mccx,alth,al,dal,ali,kka,mvka,mkax &! input 
    & ,nax,six,cax,mgx,alx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::nath,mgth,cath,sith,dt,dna,dsi,dmg,dca,nai,mgi,cai,sii,mvab,mvan,mvfo,mvcc,alth,dal,ali &
    & ,mvka
real(kind=8),dimension(nz),intent(in)::v,na,ca,mg,si,tora,poro,sat,kab,kan,kcc,kfo,hr,mabx,manx,mfox,mccx,dz,al &
    & ,kka,mkax
real(kind=8),dimension(nz),intent(out)::nax,six,cax,mgx,alx

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
        & +3d0*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz) &
        & +2d0*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz) &
        & +2d0*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz) &
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
        & +kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz) &
        & ) &
        & )
        
end do
        
edif = poro*sat*1d3*dal*tora
edifi = edif(1)
        
do iz = 1, nz


    if (alx(iz)>=alth) cycle

    ctmp = al(max(1,iz-1))
    ediftmp = edif(max(1,iz-1))
    if (iz==1) ctmp = ali
    if (iz==1) ediftmp = edifi
    
    alx(iz) = max(0.0d0, &
        & al(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
        & -poro(iz)*sat(iz)*1d3*v(iz)*(al(iz)-ctmp)/dz(iz) &
        & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(al(min(nz,iz+1))-al(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
        & - 0.5d0*(edif(iz)+ediftmp)*(al(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
        & +2d0*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz) &
        & +kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz) &
        & +2d0*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz) &
        & ) &
        & )
        
end do

endsubroutine precalc_pw_sil_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_aqs( &
    & nz,dt,v,dz,tora,poro,sat,hr &! input 
    & ,nsp_aq,nsp_sld,daq,maqth,maqi,maq,mv,msldx,ksld,staq &! input
    & ,nrxn_ext,staq_ext,rxnext &! input
    & ,maqx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt
real(kind=8)::nath,dna,nai,rxn_tmp
real(kind=8),dimension(nz),intent(in)::v,tora,poro,sat,hr,dz
real(kind=8),dimension(nz)::na,nax

integer iz,ispa,isps,irxn
real(kind=8) ctmp,edifi,ediftmp
real(kind=8),dimension(nz)::edif

integer,intent(in)::nsp_aq,nsp_sld,nrxn_ext
real(kind=8),dimension(nsp_aq),intent(in)::daq,maqth,maqi
real(kind=8),dimension(nsp_aq,nz),intent(in)::maq
real(kind=8),dimension(nsp_sld),intent(in)::mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx,ksld
real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
real(kind=8),dimension(nsp_aq,nz),intent(inout)::maqx
real(kind=8),dimension(nrxn_ext,nz),intent(in)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_aq),intent(in)::staq_ext

if (nsp_aq == 0 ) return

do ispa = 1, nsp_aq

    dna = daq(ispa)
    nath = maqth(ispa)
    nai = maqi(ispa)
    
    na = maq(ispa,:)
    nax = maqx(ispa,:)

    edif = poro*sat*1d3*dna*tora
    edifi = edif(1)

    do iz = 1, nz

        if (nax(iz)>=nath) cycle

        ctmp = na(max(1,iz-1))
        ediftmp = edif(max(1,iz-1))
        if (iz==1) ctmp = nai
        if (iz==1) ediftmp = edifi
    
        rxn_tmp = 0d0
        
        do isps = 1, nsp_sld
            if (staq(isps,ispa)>0d0) then 
                rxn_tmp = rxn_tmp  &
                    & + staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)
            endif 
        enddo 
        
        do irxn = 1, nrxn_ext
            if (staq(irxn,ispa)>0d0) then 
                rxn_tmp = rxn_tmp  &
                    & + staq_ext(irxn,ispa)*rxnext(irxn,iz)
            endif 
        enddo 
        
        nax(iz) = max(0.0d0, &
            & na(iz) +dt/(poro(iz)*sat(iz)*1d3)*( &
            & -poro(iz)*sat(iz)*1d3*v(iz)*(na(iz)-ctmp)/dz(iz) &
            & +(0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(na(min(nz,iz+1))-na(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+ediftmp)*(na(iz)-ctmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
            & +rxn_tmp &
            & ) &
            & )
            
    enddo
    
    maqx(ispa,:) = nax 

enddo 

endsubroutine precalc_aqs

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine pyweath_1D( &
    & nz,nflx_py,mvpy,c,c2,ci,c2i,po2,po2i,ms,msi,hr,po2th,poro,z,dz,w,koxs2,koxs,msth,dfe2,dfe3,sat,dporodta,dporodtg  &! input
    & ,kho,koxa,dt2,cth,c2th,stoxa,tora,torg,daqo,dgaso,v,swbr,mo2,stoxs,tol,nsp,runname,workdir,zrxn,it &! input
    & ,swoxa,swoxall,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,cx,c2x,po2x,msx,flx_py,flx_py_fe2,flx_py_fe3,flx_py_o2 &! output
    ) 
    
implicit none 

integer,intent(in)::nz,nflx_py,nsp
real(kind=8),intent(in)::ci,c2i,po2i,msi,po2th,dz,w,msth,dfe2,dfe3,kho,dt2,cth,c2th,stoxa,daqo,dgaso &
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
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
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
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
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
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)*(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & )  
            flx_py_o2(iadv,iz) = (  &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz &
                & ) 

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +1.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(1.0d0)/(dz**2.0d0) &
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
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
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
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & )  
            flx_py_o2(iadv,iz) = (  &
                & + poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & ) 

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt & 
                & +dporodtg(iz)  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) *(1.0d0)/(dz**2.0d0) &
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
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & - (ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & - 1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso)  &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(-1.0d0)/(dz**2.0d0) &
                & + poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
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
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso-poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
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
    & ,kho,tora,torg,daqo,dgaso,v,mo2,tol,runname,workdir,zrxn,ucv,vmax  &! inpput
    & ,iter,error,dt &! inout
    & ,po2x,flx_o2,resp &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::po2i,po2th,dz,kho,daqo,dgaso,mo2,tol,zrxn,ucv,vmax
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
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz &
                & +resp(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2i-2.0d0*po2x(iz))/(dz**2.0d0) &
                & )
            flx_o2(iadv,iz) = +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2i)/dz

        else if (iz == nz) then

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt &
                & +dporodtg(iz)  &
                & +(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +resp(iz)  &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz-1)-1.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & ) 
            flx_o2(iadv,iz) = poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz 

        else

            amx(row,row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)/dt & 
                & +dporodtg(iz)  &
                & +2.0d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
                & *(1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3/dz &
                & +dresp_dpo2(iz) &
                & ) &
                & *merge(1.0d0,po2x(iz),po2x(iz)<po2th)

            amx(row,row+nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & ) &
                & *merge(0.0d0,po2x(iz+1),po2x(iz)<po2th)

            amx(row,row-nsp) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso &
                & +poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo)/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso)  &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo &
                & -poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo))*(-1.0d0)/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(-1.0d0)/dz &
                & ) &
                & *merge(0.0d0,po2x(iz-1),po2x(iz)<po2th)

            ymx(row) = ( &
                & (ucv*poro(iz)*(1.0d0-sat(iz))*1d3+poro(iz)*sat(iz)*kho*1d3)*(po2x(iz)-po2(iz))/dt  &
                & +dporodtg(iz) *po2x(iz) &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
                & *(po2x(iz)-po2x(iz-1))/(dz**2.0d0) &
                & +poro(iz)*sat(iz)*v(iz)*kho*1d3*(po2x(iz)-po2x(iz-1))/dz &
                & +resp(iz) &
                & ) &
                & *merge(0.0d0,1.0d0,po2x(iz)<po2th)
            
            flx_o2(idif,iz) = ( &
                & -(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(iz)*dgaso+poro(iz)*sat(iz)*kho*1d3*tora(iz)*daqo) &
                & *(po2x(iz+1)+po2x(iz-1)-2.0d0*po2x(iz))/(dz**2.0d0) &
                & -1d3*(ucv*(poro(iz)*(1.0d0-sat(iz))*torg(iz)*dgaso &
                & -poro(iz-1)*(1.0d0-sat(iz-1))*torg(iz-1)*dgaso) &
                & +(poro(iz)*sat(iz)*kho*tora(iz)*daqo-poro(iz-1)*sat(iz-1)*kho*tora(iz-1)*daqo)) &
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
    & ,kho,tora,torg,daqo,dgaso,v,mo2,tol,runname,workdir,zrxn,ucv,vmax,poroprev  &! inpput
    & ,dt,flgback &! inout
    & ,po2x,flx_o2,resp &! output
    & ) 
! only oxygen + soil respiration 
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::po2i,po2th,kho,daqo,dgaso,mo2,tol,zrxn,ucv,vmax
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
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgaso+poro*sat*kho2x*1d3*tora*daqo
    edifi = edif(1)
    edifi = ucv*1d3*dgaso 
    
    dedif_do2 = poro*sat*dkho2_do2*1d3*tora*daqo
    
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
iter = 0

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

subroutine CO2_1D_v2_1( &
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
real(kind=8) edifi,error,pco2n_tmp,khco2n_tmp,edifn_tmp
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
iter = 0

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
        
        pco2n_tmp = pco2x(max(1,iz-1))
        khco2n_tmp = khco2x(max(1,iz-1))
        edifn_tmp = edif(max(1,iz-1))
        if (iz == 1) then 
            pco2n_tmp = pco2i
            khco2n_tmp = khco2i
            edifn_tmp = edifi
        endif 

        amx(row,row) = ( &
            & (alpha(iz) + dalpha_dco2(iz)*pco2x(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & +0.5d0*(dedif_dco2(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
            & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dco2(iz)*dt &
            & ) &
            & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

        ymx(row) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & ) &
            & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
        
        
        if (iz/=nz) then 
        
            amx(row,row+nsp) = ( &
                    & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & + 0.5d0*(dedif_dco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt&
                    & ) &
                    & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)
        
        endif 
        
        if (iz/=1) then 

            amx(row,row-nsp) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dco2(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)
        endif 
        
        flx_co2(itflx,iz) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz))/dt &
            & )         
        flx_co2(idif,iz) = ( &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & )
        flx_co2(iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2n_tmp*pco2n_tmp)/dz(iz) &
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

endsubroutine CO2_1D_v2_1

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
    & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,k1si,k2si,keqcca,authig,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,pco2x,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz,kcca
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx 
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc,omega_cca
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
    & ,domega_cc_dsi,domega_cc_dmg,dprodsi,domega_fo_dpro,domega_ab_dpro,domega_an_dpro,domega_cc_dpro &
    & ,domega_cca_dmg,domega_cca_dca,domega_cca_dna,domega_cca_dsi,domega_cca_dpro
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,st_cca
real(kind=8)::k1_fo = 10d0**(-6.85d0), E1_fo = 51.7d0, n1_fo = 0.5d0, k2_fo = 10d0**(-12.41d0),E2_fo = 38d0 &
    & ,k3_fo = 10d0**(-21.2d0),E3_fo = 94.1d0,n3_fo = -0.82d0  &
    & ,k1_ab = 10d0**(-10.16d0), E1_ab = 65d0, n1_ab = 0.457d0, k2_ab = 10d0**(-12.56d0),E2_ab = 69.8d8 &
    & ,k3_ab = 10d0**(-15.60d0),E3_ab = 71d0, n3_ab = -0.572d0 &
    & ,k1_an = 10d0**(-3.5d0),E1_an = 16.6d0,n1_an = 1.411d0,k2_an = 10d0**(-9.12d0), E2_an = 17.8d0 

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 100d0
real(kind=8)::disonly = 0d0 ! for cc [1---yes, 0---no]
! real(kind=8)::disonly = 1d0 ! for cc 
! real(kind=8)::authig = 0d0 ! for cc whether to include authigenesis [1---yes, 0---no]
! real(kind=8)::authig = 1d0 ! 

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
    ! call calc_pH( &
        ! & nz,2d0*(cx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,prox &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,2d0*(cx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,si,k1si,k2si &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v3( &
        & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,prox &! output
        & ) 
else 

    ! prox(:) = 0.5d0* ( &
        ! & -1d0*(nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))  &
        ! & + sqrt((nax(:)+2d0*ca(:)+2d0*mgx(:)+2d0*cx(:)+3d0*c2x(:)-2d0*so4x(:))**2d0  &
        ! & + 4d0*kco2*k1*pco2x(:)) &
        ! & )
    ! call calc_pH( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,prox &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v3( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,prox &! output
        & ) 
endif

! print *, 'starting silciate calculation'

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
    ! call calc_pH( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,prox &! output
        ! & ) 
    ! call calc_pH( &
        ! & nz,nax+dconc+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,dprodna &! output
        ! & ) 
    ! call calc_pH( &
        ! & nz,nax+2d0*(cx+cax+mgx+dconc-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,dprodmg &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,prox &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+dconc+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,dprodna &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+mgx+dconc-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,dprodmg &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+dconc+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,dprodca &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six+dconc,k1si,k2si &! input 
        ! & ,dprodsi &! output
        ! & ) 
    call calc_pH_v3( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,prox &! output
        & ) 
    call calc_pH_v3( &
        & nz,nax,mgx+dconc,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,dprodmg &! output
        & ) 
    call calc_pH_v3( &
        & nz,nax+dconc,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,dprodna &! output
        & ) 
    call calc_pH_v3( &
        & nz,nax,mgx,cax+dconc,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,dprodca &! output
        & ) 
    call calc_pH_v3( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six+dconc,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,dprodsi &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodsi = (dprodsi-prox)/dconc
    dprodca = (dprodca-prox)/dconc
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
        & ,pco2x,cax,mgx,six,nax,prox,'fo' &! input 
        & ,omega_fo &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,'fo' &! input 
        & ,domega_fo_dmg &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &!input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,'fo' &! input 
        & ,domega_fo_dsi &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,'fo' &! input 
        & ,domega_fo_dna &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,'fo' &! input 
        & ,domega_fo_dca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox+dconc,'fo' &! input 
        & ,domega_fo_dpro &! output
        & )
    domega_fo_dmg = (domega_fo_dmg-omega_fo)/dconc
    domega_fo_dsi = (domega_fo_dsi-omega_fo)/dconc
    domega_fo_dna = (domega_fo_dna-omega_fo)/dconc
    domega_fo_dca = (domega_fo_dca-omega_fo)/dconc
    domega_fo_dpro = (domega_fo_dpro-omega_fo)/dconc
    
    domega_fo_dmg = domega_fo_dmg + domega_fo_dpro*dprodmg
    domega_fo_dca = domega_fo_dca + domega_fo_dpro*dprodca
    domega_fo_dna = domega_fo_dna + domega_fo_dpro*dprodna
    domega_fo_dsi = domega_fo_dsi + domega_fo_dpro*dprodsi
    
    
    ! omega_fo(:) = mgx(:)**2d0*six(:)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)/keqfo
    ! domega_fo_dmg(:) = 2d0*mgx(:)*six(:)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)/keqfo &
        ! & + mgx(:)**2d0*six(:)*(-1d0)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)**2d0 &
        ! & *(4d0*prox(:)**3d0+3d0*k1si*prox(:)**2d0+2d0*k2si*prox(:))*dprodmg(:)/keqfo
    ! domega_fo_dsi(:) = mgx(:)**2d0/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)/keqfo &
        ! & + mgx(:)**2d0*six(:)*(-1d0)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)**2d0 &
        ! & *(4d0*prox(:)**3d0+3d0*k1si*prox(:)**2d0+2d0*k2si*prox(:))*dprodsi(:)/keqfo
    ! domega_fo_dna(:) = mgx(:)**2d0*six(:)*(-1d0)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)**2d0 &
        ! & *(4d0*prox(:)**3d0+3d0*k1si*prox(:)**2d0+2d0*k2si*prox(:))*dprodna(:)/keqfo
    ! domega_fo_dca(:) = mgx(:)**2d0*six(:)*(-1d0)/(prox(:)**4d0+k1si*prox(:)**3d0+k2si*prox(:)**2d0)**2d0 &
        ! & *(4d0*prox(:)**3d0+3d0*k1si*prox(:)**2d0+2d0*k2si*prox(:))*dprodca(:)/keqfo
    
    ! omega_fo(:) = mg(:)**2d0*si(:)/(pro(:)**4d0)/keqfo
    ! domega_fo_dmg(:) = 0d0
    ! domega_fo_dsi(:) = 0d0
    
    ! print *,omega_fo
    ! print *,domega_fo_dmg
    ! print *,domega_fo_dsi
    ! stop
    
    ! Ab + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)
    
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,'ab' &! input 
        & ,omega_ab &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx+dconc,six,nax,prox,'ab' &! input 
        & ,domega_ab_dmg &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six+dconc,nax,prox,'ab' &! input 
        & ,domega_ab_dsi &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,'ab' &! input 
        & ,domega_ab_dna &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax+dconc,mgx,six,nax,prox,'ab' &! input 
        & ,domega_ab_dca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox+dconc,'ab' &! input 
        & ,domega_ab_dpro &! output
        & )
    domega_ab_dmg = (domega_ab_dmg-omega_ab)/dconc
    domega_ab_dsi = (domega_ab_dsi-omega_ab)/dconc
    domega_ab_dna = (domega_ab_dna-omega_ab)/dconc
    domega_ab_dca = (domega_ab_dca-omega_ab)/dconc
    domega_ab_dpro = (domega_ab_dpro-omega_ab)/dconc
    
    domega_ab_dmg = domega_ab_dmg + domega_ab_dpro*dprodmg
    domega_ab_dca = domega_ab_dca + domega_ab_dpro*dprodca
    domega_ab_dna = domega_ab_dna + domega_ab_dpro*dprodna
    domega_ab_dsi = domega_ab_dsi + domega_ab_dpro*dprodsi
    
    ! omega_ab(:) = nax(:)*six(:)**2d0/prox(:)/keqab
    ! domega_ab_dna(:) = six(:)**2d0/prox(:)/keqab + nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodna(:)
    ! domega_ab_dsi(:) = nax(:)*(2d0)*six(:)/prox(:)/keqab + nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodsi(:)
    ! domega_ab_dmg(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodmg(:)
    ! domega_ab_dca(:) = nax(:)*six(:)**2d0*(-1d0)/(prox(:)**2d0)/keqab*dprodca(:)
    
    ! An + 2H+ + H2O = kaolinite + Ca2+ 
    
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,'an' &! input 
        & ,omega_an &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx+dconc,six,nax,prox,'an' &! input 
        & ,domega_an_dmg &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,'an' &! input 
        & ,domega_an_dsi &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax+dconc,prox,'an' &! input 
        & ,domega_an_dna &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,'an' &! input 
        & ,domega_an_dca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,'an' &! input 
        & ,domega_an_dpro &! output
        & )
    domega_an_dmg = (domega_an_dmg-omega_an)/dconc
    domega_an_dsi = (domega_an_dsi-omega_an)/dconc
    domega_an_dna = (domega_an_dna-omega_an)/dconc
    domega_an_dca = (domega_an_dca-omega_an)/dconc
    domega_an_dpro = (domega_an_dpro-omega_an)/dconc
    
    domega_an_dmg = domega_an_dmg + domega_an_dpro*dprodmg
    domega_an_dca = domega_an_dca + domega_an_dpro*dprodca
    domega_an_dna = domega_an_dna + domega_an_dpro*dprodna
    domega_an_dsi = domega_an_dsi + domega_an_dpro*dprodsi
    
    ! omega_an(:) = cax(:)/(prox(:)**2d0)/keqan
    ! domega_an_dca(:) = 1d0/(prox(:)**2d0)/keqan + cax(:)*(-2d0)/(prox(:)**3d0)*dprodca(:)/keqan
    ! domega_an_dna(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodna(:)/keqan
    ! domega_an_dmg(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodmg(:)/keqan
    ! domega_an_dsi(:) = cax(:)*(-2d0)/(prox(:)**3d0)*dprodsi(:)/keqan
    
    ! Cc = Ca2+ + CO32- 
    
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,'cc' &! input 
        & ,omega_cc &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,'cc' &! input 
        & ,domega_cc_dmg &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,'cc' &! input 
        & ,domega_cc_dsi &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax+dconc,prox,'cc' &! input 
        & ,domega_cc_dna &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,'cc' &! input 
        & ,domega_cc_dca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,'cc' &! input 
        & ,domega_cc_dpro &! output
        & )
    domega_cc_dmg = (domega_cc_dmg-omega_cc)/dconc
    domega_cc_dsi = (domega_cc_dsi-omega_cc)/dconc
    domega_cc_dna = (domega_cc_dna-omega_cc)/dconc
    domega_cc_dca = (domega_cc_dca-omega_cc)/dconc
    domega_cc_dpro = (domega_cc_dpro-omega_cc)/dconc
    
    domega_cc_dmg = domega_cc_dmg + domega_cc_dpro*dprodmg
    domega_cc_dca = domega_cc_dca + domega_cc_dpro*dprodca
    domega_cc_dna = domega_cc_dna + domega_cc_dpro*dprodna
    domega_cc_dsi = domega_cc_dsi + domega_cc_dpro*dprodsi
    
    ! omega_cc = cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc
    ! domega_cc_dca = k1*k2*kco2*pco2x/(prox**2d0)/keqcc + cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodca/keqcc
    ! domega_cc_dmg = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodmg/keqcc
    ! domega_cc_dna = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodna/keqcc
    ! domega_cc_dsi = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodsi/keqcc
    
    ! authigenesis (prety much tentative parameterization)
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,'cc' &! input 
        & ,omega_cca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,'cc' &! input 
        & ,domega_cca_dmg &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,'cc' &! input 
        & ,domega_cca_dsi &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,'cc' &! input 
        & ,domega_cca_dna &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,'cc' &! input 
        & ,domega_cca_dca &! output
        & )
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,'cc' &! input 
        & ,domega_cca_dpro &! output
        & )
    domega_cca_dmg = (domega_cca_dmg-omega_cca)/dconc
    domega_cca_dsi = (domega_cca_dsi-omega_cca)/dconc
    domega_cca_dna = (domega_cca_dna-omega_cca)/dconc
    domega_cca_dca = (domega_cca_dca-omega_cca)/dconc
    domega_cca_dpro = (domega_cca_dpro-omega_cca)/dconc
    
    domega_cca_dmg = domega_cca_dmg + domega_cca_dpro*dprodmg
    domega_cca_dca = domega_cca_dca + domega_cca_dpro*dprodca
    domega_cca_dna = domega_cca_dna + domega_cca_dpro*dprodna
    domega_cca_dsi = domega_cca_dsi + domega_cca_dpro*dprodsi

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
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dmg(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dsi(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dna(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dca(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                
                if (authig == 1d0) then 
                    ymx3(row) = ymx3(row) + ( & 
                        & + kcca(iz)*poro(iz)*sat(iz)*(1d0-omega_cca(iz))*authig*dt &
                        & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                        & ) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                endif
                    
                flx_cc(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_cc(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_cc(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(1d0-omega_cca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
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
                st_cca = 0d0
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
                st_cca = 0d0
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
                st_cca = 0d0
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
                st_cca = 1d0
                rxn_tmp = st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(1d0-omega_cc(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0) &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(1d0-omega_cca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) 
                drxndisp_tmp = ( &
                    & st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0) &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
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
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dmg(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dsi(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dna(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
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
    ! call calc_pH( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2 &! input 
        ! & ,prox &! output
        ! & ) 
    ! call calc_pH_v2( &
        ! & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v3( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
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

subroutine calc_pH_v2( &
    & nz,netcat,pco2x,kw,kco2,k1,k2,six,k1si,k2si &! input 
    & ,prox &! output
    & ) 
! solving charge balance:
! [H+] + ZX[Xz+] - ZY[YZ-] - [HCO3-] - 2[CO32-] - [OH-] - [H3SiO4-] - 2[H2SiO42-] = 0
! [H+] + ZX[Xz+] - ZY[YZ-] - k1kco2pCO2/[H+] - 2k2k1kco2pCO2/[H+]^2 - kw/[H+] - [Si]/([H+]/k1si + 1 + k2si/k1si/[H+])
!       - 2[Si]/([H+]^2/k2si + [H+]k1si/k2si + 1) = 0
! [H+]^3 + (ZX[Xz+] - ZY[YZ-])[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
! NetCat is defined as (ZX[Xz+] - ZY[YZ-])
! [H+]^3 + NetCat[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
implicit none
integer,intent(in)::nz
real(kind=8),intent(in)::kw,kco2,k1,k2,k1si,k2si
real(kind=8),dimension(nz),intent(in)::netcat,pco2x,six
real(kind=8),dimension(nz),intent(inout)::prox

real(kind=8),dimension(nz)::df,f
real(kind=8) error,tol
integer iter,iz

error = 1d4
tol = 1d-6

prox = 1d0 
iter = 0

! print *,k1si,k2si

! print*,'calc_pH'
do while (error > tol)
    f = prox + netcat - k1*kco2*pco2x/prox - 2d0*k2*k1*kco2*pco2x/(prox**2d0) - kw/prox &
        & - six/(prox/k1si + 1d0 + k2si/k1si/prox) !- 2d0*six/(prox**2d0/k2si + prox*k1si/k2si + 1d0)
    df =1d0  - k1*kco2*pco2x*(-1d0)/(prox**2d0) - 2d0*k2*k1*kco2*pco2x*(-2d0)/(prox**3d0) - kw*(-1d0)/(prox**2d0) &
        & - six*(-1d0)/((prox/k1si + 1d0 + k2si/k1si/prox)**2d0)*(1d0/k1si + k2si/k1si*(-1d0)/(prox**2d0)) !&
        ! & - 2d0*six*(-1d0)/((prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0)*(2d0*prox/k2si + k1si/k2si ) 
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        & - six*prox**2d0/(prox/k1si + 1d0 + k2si/k1si/prox)  &
        & - 2d0*six*prox**2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0)
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw) &
        & - six*prox*2d0/(prox/k1si + 1d0 + k2si/k1si/prox) &
        & - six*prox**2d0*(-1d0)/(prox/k1si + 1d0 + k2si/k1si/prox)**2d0* (1d0/k1si + k2si/k1si*(-1d0)/prox**2d0) &
        & - 2d0*six*prox*2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & - 2d0*six*prox**2d0*(-1d0)/(prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0*(prox*2d0/k2si + k1si/k2si) 
    df = df*prox
    ! if (any(isnan(-f/df)) .or. any(isnan(exp(-f/df)))) then 
        ! print *,any(isnan(-f/df)),any(isnan(exp(-f/df)))
        ! print *,-f/df
    ! endif 
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

endsubroutine calc_pH_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH_v3( &
    & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
    & ,prox &! output
    & ) 
! solving charge balance:
! [H+] + ZX[Xz+] - ZY[YZ-] - [HCO3-] - 2[CO32-] - [OH-] - [H3SiO4-] - 2[H2SiO42-] = 0
! [H+] + ZX[Xz+] - ZY[YZ-] - k1kco2pCO2/[H+] - 2k2k1kco2pCO2/[H+]^2 - kw/[H+] - [Si]/([H+]/k1si + 1 + k2si/k1si/[H+])
!       - 2[Si]/([H+]^2/k2si + [H+]k1si/k2si + 1) = 0
! [H+]^3 + (ZX[Xz+] - ZY[YZ-])[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
! NetCat is defined as (ZX[Xz+] - ZY[YZ-])
! [H+]^3 + NetCat[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
implicit none
integer,intent(in)::nz
real(kind=8),intent(in)::kw,kco2,k1,k2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3
real(kind=8),dimension(nz),intent(in)::nax,mgx,cax,so4x,pco2x,six
real(kind=8),dimension(nz),intent(inout)::prox

real(kind=8),dimension(nz)::df,f,netcat
real(kind=8) error,tol
integer iter,iz

error = 1d4
tol = 1d-6

prox = 1d0 
iter = 0

! netcat = nax+2d0*cax-2d0*so4x
netcat = nax-2d0*so4x

! print *,k1si,k2si

! print*,'calc_pH'
do while (error > tol)
    f = prox + netcat - k1*kco2*pco2x/prox - 2d0*k2*k1*kco2*pco2x/(prox**2d0) - kw/prox &
        & - six/(prox/k1si + 1d0 + k2si/k1si/prox) !- 2d0*six/(prox**2d0/k2si + prox*k1si/k2si + 1d0)
        ! + 2d0*mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) 
        ! + mgx/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) 
    df =1d0  - k1*kco2*pco2x*(-1d0)/(prox**2d0) - 2d0*k2*k1*kco2*pco2x*(-2d0)/(prox**3d0) - kw*(-1d0)/(prox**2d0) &
        & - six*(-1d0)/((prox/k1si + 1d0 + k2si/k1si/prox)**2d0)*(1d0/k1si + k2si/k1si*(-1d0)/(prox**2d0)) !&
        ! & - 2d0*six*(-1d0)/((prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0)*(2d0*prox/k2si + k1si/k2si ) 
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        & - six*prox**2d0/(prox/k1si + 1d0 + k2si/k1si/prox)  &
        & - 2d0*six*prox**2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & + 2d0*mgx*prox**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + mgx*prox**2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & + 2d0*cax*prox**2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + cax*prox**2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) 
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw) &
        & - six*prox*2d0/(prox/k1si + 1d0 + k2si/k1si/prox) &
        & - six*prox**2d0*(-1d0)/(prox/k1si + 1d0 + k2si/k1si/prox)**2d0* (1d0/k1si + k2si/k1si*(-1d0)/prox**2d0) &
        & - 2d0*six*prox*2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & - 2d0*six*prox**2d0*(-1d0)/(prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0*(prox*2d0/k2si + k1si/k2si) &
        & + 2d0*mgx*prox*2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*mgx*prox**2d0*(-1d0) &
        &   /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + mgx*prox*2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & + mgx*prox**2d0*(-1d0) &
        &   /(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0)**2d0 & 
        &   *(1d0/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3*(-1d0)/prox**2d0)  &
        & + 2d0*cax*prox*2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*cax*prox**2d0*(-1d0) &
        &   /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + cax*prox*2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        & + cax*prox**2d0*(-1d0) &
        &   /(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0)**2d0 & 
        &   *(1d0/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3*(-1d0)/prox**2d0)  
    df = df*prox
    ! if (any(isnan(-f/df)) .or. any(isnan(exp(-f/df)))) then 
        ! print *,any(isnan(-f/df)),any(isnan(exp(-f/df)))
        ! print *,-f/df
    ! endif 
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

endsubroutine calc_pH_v3

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH_v4( &
    & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
    & ,alx,k1al,k2al,k3al,k4al &! input
    & ,prox &! output
    & ) 
! solving charge balance:
! [H+] + ZX[Xz+] - ZY[YZ-] - [HCO3-] - 2[CO32-] - [OH-] - [H3SiO4-] - 2[H2SiO42-] = 0
! [H+] + ZX[Xz+] - ZY[YZ-] - k1kco2pCO2/[H+] - 2k2k1kco2pCO2/[H+]^2 - kw/[H+] - [Si]/([H+]/k1si + 1 + k2si/k1si/[H+])
!       - 2[Si]/([H+]^2/k2si + [H+]k1si/k2si + 1) = 0
! [H+]^3 + (ZX[Xz+] - ZY[YZ-])[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
! NetCat is defined as (ZX[Xz+] - ZY[YZ-])
! [H+]^3 + NetCat[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
implicit none
integer,intent(in)::nz
real(kind=8),intent(in)::kw,kco2,k1,k2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,k1al,k2al,k3al,k4al
real(kind=8),dimension(nz),intent(in)::nax,mgx,cax,so4x,pco2x,six,alx
real(kind=8),dimension(nz),intent(inout)::prox

real(kind=8),dimension(nz)::df,f,netcat
real(kind=8) error,tol
integer iter,iz

error = 1d4
tol = 1d-6

prox = 1d0 
iter = 0

! netcat = nax+2d0*cax-2d0*so4x
netcat = nax-2d0*so4x

! print *,k1si,k2si

! print*,'calc_pH'
do while (error > tol)
    f = prox + netcat - k1*kco2*pco2x/prox - 2d0*k2*k1*kco2*pco2x/(prox**2d0) - kw/prox &
        & - six/(prox/k1si + 1d0 + k2si/k1si/prox) !- 2d0*six/(prox**2d0/k2si + prox*k1si/k2si + 1d0)
        ! + 2d0*mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) 
        ! + mgx/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) 
        ! + 3d0*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)
        ! + 2d0*alx*k1al/prox/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)
        ! + alx*k2al/prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)
        ! - alx*k4al/prox**4d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)
    df =1d0  - k1*kco2*pco2x*(-1d0)/(prox**2d0) - 2d0*k2*k1*kco2*pco2x*(-2d0)/(prox**3d0) - kw*(-1d0)/(prox**2d0) &
        & - six*(-1d0)/((prox/k1si + 1d0 + k2si/k1si/prox)**2d0)*(1d0/k1si + k2si/k1si*(-1d0)/(prox**2d0)) !&
        ! & - 2d0*six*(-1d0)/((prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0)*(2d0*prox/k2si + k1si/k2si ) 
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        & - six*prox**2d0/(prox/k1si + 1d0 + k2si/k1si/prox)  &
        & - 2d0*six*prox**2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & + 2d0*mgx*prox**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + mgx*prox**2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & + 2d0*cax*prox**2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + cax*prox**2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        & + 3d0*alx*prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 2d0*alx*k1al*prox/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + alx*k2al/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & - alx*k4al/prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw) &
        & - six*prox*2d0/(prox/k1si + 1d0 + k2si/k1si/prox) &
        & - six*prox**2d0*(-1d0)/(prox/k1si + 1d0 + k2si/k1si/prox)**2d0* (1d0/k1si + k2si/k1si*(-1d0)/prox**2d0) &
        & - 2d0*six*prox*2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & - 2d0*six*prox**2d0*(-1d0)/(prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0*(prox*2d0/k2si + k1si/k2si) &
        & + 2d0*mgx*prox*2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*mgx*prox**2d0*(-1d0) &
        &   /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + mgx*prox*2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & + mgx*prox**2d0*(-1d0) &
        &   /(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0)**2d0 & 
        &   *(1d0/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3*(-1d0)/prox**2d0)  &
        & + 2d0*cax*prox*2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*cax*prox**2d0*(-1d0) &
        &   /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + cax*prox*2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        & + cax*prox**2d0*(-1d0) &
        &   /(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0)**2d0 & 
        &   *(1d0/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3*(-1d0)/prox**2d0)   &
        & + 3d0*alx*prox*2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 3d0*alx*prox**2d0*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & + 2d0*alx*k1al/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 2d0*alx*k1al*prox*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & + alx*k2al*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & - alx*k4al*(-2d0)/prox**3d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & - alx*k4al/prox**2d0*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) 
    df = df*prox
    ! if (any(isnan(-f/df)) .or. any(isnan(exp(-f/df)))) then 
        ! print *,any(isnan(-f/df)),any(isnan(exp(-f/df)))
        ! print *,-f/df
    ! endif 
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

endsubroutine calc_pH_v4

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH_v5( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
    & ,print_cb,print_loc,z &! input 
    & ,prox,ph_error &! output
    & ) 
! solving charge balance:
! [H+] + ZX[Xz+] - ZY[YZ-] - [HCO3-] - 2[CO32-] - [OH-] - [H3SiO4-] - 2[H2SiO42-] = 0
! [H+] + ZX[Xz+] - ZY[YZ-] - k1kco2pCO2/[H+] - 2k2k1kco2pCO2/[H+]^2 - kw/[H+] - [Si]/([H+]/k1si + 1 + k2si/k1si/[H+])
!       - 2[Si]/([H+]^2/k2si + [H+]k1si/k2si + 1) = 0
! [H+]^3 + (ZX[Xz+] - ZY[YZ-])[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
! NetCat is defined as (ZX[Xz+] - ZY[YZ-])
! [H+]^3 + NetCat[H+]^2 - (k1kco2pCO2+kw)[H+] - 2k2k1kco2pCO2  = 0
implicit none
integer,intent(in)::nz
real(kind=8),intent(in)::kw
real(kind=8) kco2,k1,k2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,k1al,k2al,k3al,k4al &
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3
real(kind=8),dimension(nz)::nax,mgx,cax,so4x,pco2x,six,alx,fe2x,fe3x,kx
real(kind=8),dimension(nz),intent(in)::z
real(kind=8),dimension(nz),intent(inout)::prox
logical,intent(out)::ph_error

real(kind=8),dimension(nz)::df,f,netcat
real(kind=8) error,tol
integer iter,iz

integer,intent(in)::nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst
character(5),dimension(nsp_aq),intent(in)::chraq
character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

logical,intent(in)::print_cb
character(500),intent(in)::print_loc

error = 1d4
tol = 1d-6

prox = 1d0 
iter = 0


kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)
k1si = keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h1)
k2si = keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h2)
k1mg = keqaq_h(findloc(chraq_all,'mg',dim=1),ieqaq_h1)
k1mgco3 = keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_co3)
k1mghco3  = keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_hco3)
k1ca = keqaq_h(findloc(chraq_all,'ca',dim=1),ieqaq_h1)
k1caco3 = keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_co3)
k1cahco3 = keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_hco3)
k1al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h1)
k2al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h2)
k3al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h3)
k4al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h4)
k1fe2 = keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1)
k1fe2co3 = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3)
k1fe2hco3  = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3)
k1fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h1)
k2fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h2)
k3fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h3)
k4fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h4)

nax = 0d0
so4x = 0d0
kx = 0d0
if (any(chraq=='na')) then 
    nax = maqx(findloc(chraq,'na',dim=1),:)
elseif (any(chraq_cnst=='na')) then 
    nax = maqc(findloc(chraq_cnst,'na',dim=1),:)
endif 
if (any(chraq=='k')) then 
    kx = maqx(findloc(chraq,'k',dim=1),:)
elseif (any(chraq_cnst=='k')) then 
    kx = maqc(findloc(chraq_cnst,'k',dim=1),:)
endif 
if (any(chraq=='so4')) then 
    so4x = 2d0*maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = 2d0*maqc(findloc(chraq_cnst,'so4',dim=1),:)
endif 

netcat = nax + kx -2d0*so4x

six =0d0
cax =0d0
mgx =0d0
fe2x =0d0
alx =0d0
fe3x =0d0
pco2x =0d0

if (any(chraq=='si')) then 
    six = maqx(findloc(chraq,'si',dim=1),:)
elseif (any(chraq_cnst=='si')) then 
    six = maqc(findloc(chraq_cnst,'si',dim=1),:)
endif 
if (any(chraq=='ca')) then 
    cax = maqx(findloc(chraq,'ca',dim=1),:)
elseif (any(chraq_cnst=='ca')) then 
    cax = maqc(findloc(chraq_cnst,'ca',dim=1),:)
endif 
if (any(chraq=='mg')) then 
    mgx = maqx(findloc(chraq,'mg',dim=1),:)
elseif (any(chraq_cnst=='mg')) then 
    mgx = maqc(findloc(chraq_cnst,'mg',dim=1),:)
endif 
if (any(chraq=='fe2')) then 
    fe2x = maqx(findloc(chraq,'fe2',dim=1),:)
elseif (any(chraq_cnst=='fe2')) then 
    fe2x = maqc(findloc(chraq_cnst,'fe2',dim=1),:)
endif 
if (any(chraq=='al')) then 
    alx = maqx(findloc(chraq,'al',dim=1),:)
elseif (any(chraq_cnst=='al')) then 
    alx = maqc(findloc(chraq_cnst,'al',dim=1),:)
endif 
if (any(chraq=='fe3')) then 
    fe3x = maqx(findloc(chraq,'fe3',dim=1),:)
elseif (any(chraq_cnst=='fe3')) then 
    fe3x = maqc(findloc(chraq_cnst,'fe3',dim=1),:)
endif 
if (any(chrgas=='pco2')) then 
    pco2x = mgasx(findloc(chrgas,'pco2',dim=1),:)
elseif (any(chrgas_cnst=='pco2')) then 
    pco2x = mgasc(findloc(chrgas_cnst,'pco2',dim=1),:)
endif 

if (any(isnan(maqx)) .or. any(isnan(maqc))) then 
    print*,'nan in input aqueosu species'
    stop
endif 

200 continue

! print*,'calc_pH'
do while (error > tol)
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        ! si 
        ! & - six*prox**2d0/(prox/k1si + 1d0 + k2si/k1si/prox)  &
        ! & - 2d0*six*prox**2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & - six*k1si*prox/(1d0+k1si/prox+k2si/prox**2d0) &
        & - 2d0*six*k2si/(1d0+k1si/prox+k2si/prox**2d0) &
        ! mg
        & + 2d0*mgx*prox**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        ! & + mgx*prox**2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & + mgx*k1mg*prox/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + mgx*k1mghco3*k1*k2*kco2*pco2x*prox/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        ! ca
        & + 2d0*cax*prox**2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        ! & + cax*prox**2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        & + cax*k1ca*prox/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + cax*k1cahco3*k1*k2*kco2*pco2x*prox/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        ! al
        & + 3d0*alx*prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 2d0*alx*k1al*prox/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + alx*k2al/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & - alx*k4al/prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        ! fe2
        & + 2d0*fe2x*prox**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        ! & + fe2x*prox**2d0/(prox/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2co3/k1fe2hco3/prox+1d0) &
        & + fe2x*k1fe2*prox/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        & + fe2x*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        &       /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        ! fe3
        & + 3d0*fe3x*prox**2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & + 2d0*fe3x*k1fe3*prox/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & + fe3x*k2fe3/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & - fe3x*k4fe3/prox**2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) 
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw) &
        !
        ! si
        !
        ! & - six*prox*2d0/(prox/k1si + 1d0 + k2si/k1si/prox) &
        ! & - six*prox**2d0*(-1d0)/(prox/k1si + 1d0 + k2si/k1si/prox)**2d0* (1d0/k1si + k2si/k1si*(-1d0)/prox**2d0) &
        ! & - 2d0*six*prox*2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        ! & - 2d0*six*prox**2d0*(-1d0)/(prox**2d0/k2si + prox*k1si/k2si + 1d0)**2d0*(prox*2d0/k2si + k1si/k2si) &
        & - six*k1si*1d0/(1d0+k1si/prox+k2si/prox**2d0) &
        & - six*k1si*prox*(-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**2d0* (k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
        & - 2d0*six*k2si*(-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**2d0 *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
        !
        ! mg
        !
        & + 2d0*mgx*prox*2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*mgx*prox**2d0*(-1d0) &
        &   /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + mgx*k1mg*1d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + mgx*k1mg*prox*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + mgx*k1mghco3*k1*k2*kco2*pco2x*1d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
        & + mgx*k1mghco3*k1*k2*kco2*pco2x*prox*(-1d0) &
        &   /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        ! & + mgx*prox*2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        ! & + mgx*prox**2d0*(-1d0) &
        ! &   /(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0)**2d0 & 
        ! &   *(1d0/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3*(-1d0)/prox**2d0)  &
        ! 
        ! ca
        !
        & + 2d0*cax*prox*2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*cax*prox**2d0*(-1d0) &
        &   /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        ! & + cax*prox*2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        ! & + cax*prox**2d0*(-1d0) &
        ! &   /(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0)**2d0 & 
        ! &   *(1d0/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3*(-1d0)/prox**2d0)   &
        & + cax*k1ca*1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + cax*k1ca*prox*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + cax*k1cahco3*k1*k2*kco2*pco2x*1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
        & + cax*k1cahco3*k1*k2*kco2*pco2x*prox*(-1d0) &
        &   /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        !
        ! al
        !
        & + 3d0*alx*prox*2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 3d0*alx*prox**2d0*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & + 2d0*alx*k1al/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & + 2d0*alx*k1al*prox*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & + alx*k2al*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        & - alx*k4al*(-2d0)/prox**3d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
        & - alx*k4al/prox**2d0*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
        &   *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
        !
        ! fe2
        ! 
        & + 2d0*fe2x*prox*2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        & + 2d0*fe2x*prox**2d0*(-1d0) &
        &   /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        ! & + fe2x*prox*2d0/(prox/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2co3/k1fe2hco3/prox+1d0) &
        ! & + fe2x*prox**2d0*(-1d0) &
        ! &   /(prox/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2co3/k1fe2hco3/prox+1d0)**2d0 & 
        ! &   *(1d0/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2co3/k1fe2hco3*(-1d0)/prox**2d0)  &
        & + fe2x*k1fe2*1d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        & + fe2x*k1fe2*prox*(-1d0)/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & + fe2x*k1fe2hco3*k1*k2*kco2*pco2x*1d0 &
        &       /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
        & + fe2x*k1fe2hco3*k1*k2*kco2*pco2x*prox*(-1d0) &
        &   /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &   *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        !
        ! fe3 
        !
        & + 3d0*fe3x*prox*2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & + 3d0*fe3x*prox**2d0*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)**2d0 &
        &   *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0) &
        & + 2d0*fe3x*k1fe3/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & + 2d0*fe3x*k1fe3*prox*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)**2d0 &
        &   *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0) &
        & + fe3x*k2fe3*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)**2d0 &
        &   *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0) &
        & - fe3x*k4fe3*(-2d0)/prox**3d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
        & - fe3x*k4fe3/prox**2d0*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)**2d0 &
        &   *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0)
    df = df*prox
    if (any(isnan(-f/df)) .or. any(isnan(exp(-f/df)))) then 
        print *,any(isnan(-f/df)),any(isnan(exp(-f/df)))
        print *,-f/df
        print *,f
        print *,df
        print *,prox
        print * &
        & ,any(isnan(prox**3d0)), any(isnan(+ netcat*prox**2d0)), any(isnan(- (k1*kco2*pco2x+kw)*prox)) &
        & ,any(isnan( - 2d0*k2*k1*kco2*pco2x )) &
        & ,'si' &
        ! & - six*prox**2d0/(prox/k1si + 1d0 + k2si/k1si/prox)  &
        ! & - 2d0*six*prox**2d0/(prox**2d0/k2si + prox*k1si/k2si + 1d0) &
        & ,any(isnan(- six*k1si*prox/(1d0+k1si/prox+k2si/prox**2d0))) &
        & ,any(isnan(- 2d0*six*k2si/(1d0+k1si/prox+k2si/prox**2d0))) &
        & ,'mg' &
        & ,any(isnan(+ 2d0*mgx*prox**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox))) &
        ! & + mgx*prox**2d0/(prox/(k1mghco3*k1*k2*kco2*pco2x)+k1mg/(k1mghco3*k1*k2*kco2*pco2x)+k1mgco3/k1mghco3/prox+1d0) &
        & ,any(isnan(+ mgx*k1mg*prox/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox))) &
        & ,any(isnan(+ mgx*k1mghco3*k1*k2*kco2*pco2x*prox &
        &       /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox))) &
        & ,'ca' &
        & ,any(isnan(+ 2d0*cax*prox**2d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox))) &
        ! & + cax*prox**2d0/(prox/(k1cahco3*k1*k2*kco2*pco2x)+k1ca/(k1cahco3*k1*k2*kco2*pco2x)+k1caco3/k1cahco3/prox+1d0) &
        & ,any(isnan(+ cax*k1ca*prox/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox))) &
        & ,any(isnan(+ cax*k1cahco3*k1*k2*kco2*pco2x*prox &
        &       /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox))) &
        & ,'al' &
        & ,any(isnan(+ 3d0*alx*prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0))) &
        & ,any(isnan(+ 2d0*alx*k1al*prox/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0))) &
        & ,any(isnan(+ alx*k2al/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0))) &
        & ,any(isnan(- alx*k4al/prox**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0))) &
        & ,'fe2' &
        & ,any(isnan(+ 2d0*fe2x*prox**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox))) &
        ! & + fe2x*prox**2d0/(prox/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2/(k1fe2hco3*k1*k2*kco2*pco2x)+k1fe2co3/k1fe2hco3/prox+1d0) &
        & ,any(isnan(+ fe2x*k1fe2*prox/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox))) &
        & ,any(isnan(+ fe2x*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        &       /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox))) &
        & ,'fe3'  &
        & ,any(isnan(+ 3d0*fe3x*prox**2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0))) &
        & ,any(isnan(+ 2d0*fe3x*k1fe3*prox/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0))) &
        & ,any(isnan(+ fe3x*k2fe3/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0))) &
        & ,any(isnan(- fe3x*k4fe3/prox**2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)))  &
        & ,any(isnan(- fe3x*k4fe3/prox**2d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)))  
        ! stop
        ! prox = 10d0
        ! prox = -netcat
        ! goto 200
    endif 
    prox = prox*exp( -f/df )
    ! do iz = 1, nz
        ! if (-f(iz)/df(iz) > 100d0) then 
            ! prox(iz) = prox(iz)*1.5d0
        ! elseif (-f(iz)/df(iz) > -100d0) then 
            ! prox(iz) = prox(iz)*0.5d0
        ! else 
            ! prox(iz) = prox(iz)*exp( -f(iz)/df(iz) )
        ! endif 
    ! enddo
    error = maxval(abs(exp( -f/df )-1d0))
    if (isnan(error)) error = 1d4
    ! print*, iter,error
    ! print*,  (-log10(prox(iz)),iz=1,nz,nz/5)
    ! print*,  (-log10(f(iz)),iz=1,nz,nz/5)
    ! print*,  (-log10(df(iz)),iz=1,nz,nz/5)
    ! pause
    ! stop
    ! where (prox == 0)
        ! prox = 1d-14
    ! endwhere
    
    iter = iter + 1
enddo 
ph_error = .false.
if (any(isnan(prox))) then     
    print *, (-log10(prox(iz)),iz=1,nz,nz/5)
    print*,'ph is nan'
    ph_error = .true.
    ! stop
endif 

if (print_cb) then 
    open(88,file = trim(adjustl(print_loc)),status='replace')
    write(88,*) ' z ',' h+ ',' oh- ',' na+ ', ' k+ ',' so42- ', 'hco3- ', ' co32- ',' h4sio4 ',' h3sio4- ',' h2sio42- ' &
        & ,' mg2+ ', ' mg(oh)+ ', ' mgco3 ', 'mghco3+ ', ' ca2+ ', ' ca(oh)+ ', ' caco3 ', ' cahco3+ ' &
        & ,' al3+ ', ' al(oh)2+ ', ' al(oh)2+ ', ' al(oh)3 ', ' al(oh)4- ' &
        & , ' fe22+ ', ' fe2(oh)+ ', ' fe2co3 ', ' fe2hco3+ ' &
        & ,' fe3+ ', ' fe3(oh)2+ ', ' fe3(oh)2+ ', ' fe3(oh)3 ', ' fe3(oh)4- ', ' total_charge ' 
    do iz=1,nz
        write(88,*) z(iz) &
        &, prox(iz) &
        & ,kw/prox(iz) &
        & ,nax(iz) &
        & ,kx(iz) &
        & ,so4x(iz) &
        & ,k1*kco2*pco2x(iz)/prox(iz) &
        & ,k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0) &
        & ,six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0)*k1si/prox(iz) &
        & ,six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0)*k2si/prox(iz)**2d0 &
        & ,mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & ,mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mg/prox(iz)  &
        & ,mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & ,mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & ,cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & ,cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1ca/prox(iz) &
        & ,cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0) &
        & ,alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k1al/prox(iz) &
        & ,alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k2al/prox(iz)**2d0 &
        & ,alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k3al/prox(iz)**3d0 &
        & ,alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k4al/prox(iz)**4d0 &
        & ,fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & ,fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2/prox(iz) &
        & ,fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0) &
        & ,fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k1fe3/prox(iz) &
        & ,fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k2fe3/prox(iz)**2d0 &
        & ,fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k3fe3/prox(iz)**3d0 &
        & ,fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k4fe3/prox(iz)**4d0 &
        ! charge balance 
        & ,1d0*prox(iz) &
        & +(-1d0)*kw/prox(iz) &
        & +(1d0)*nax(iz) &
        & +(1d0)*kx(iz) &
        & +(-2d0)*so4x(iz) &
        & +(-1d0)*k1*kco2*pco2x(iz)/prox(iz) &
        & +(-2d0)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0) &
        & +(-1d0)*six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0)*k1si/prox(iz) &
        & +(-2d0)*six(iz)/(1d0 + k1si/prox(iz) + k2si/prox(iz)**2d0)*k2si/prox(iz)**2d0 &
        & +(2d0)*mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & +(1d0)*mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mg/prox(iz)  &
        & +(0d0)*mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & +(1d0)*mgx(iz)/(1d0+k1mg/prox(iz)+k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & +(2d0)*cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & +(1d0)*cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1ca/prox(iz) &
        & +(0d0)*cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*cax(iz)/(1d0+k1ca/prox(iz)+k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(3d0)*alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0) &
        & +(2d0)*alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k1al/prox(iz) &
        & +(1d0)*alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k2al/prox(iz)**2d0 &
        & +(0d0)*alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k3al/prox(iz)**3d0 &
        & +(-1d0)*alx(iz)/(1d0+k1al/prox(iz)+k2al/prox(iz)**2d0+k3al/prox(iz)**3d0+k4al/prox(iz)**4d0)*k4al/prox(iz)**4d0 &
        & +(2d0)*fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & +(1d0)*fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2/prox(iz) &
        & +(0d0)*fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*fe2x(iz)/(1d0+k1fe2/prox(iz)+k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(3d0)*fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0) &
        & +(2d0)*fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k1fe3/prox(iz) &
        & +(1d0)*fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k2fe3/prox(iz)**2d0 &
        & +(0d0)*fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k3fe3/prox(iz)**3d0 &
        & +(-1d0)*fe3x(iz)/(1d0+k1fe3/prox(iz)+k2fe3/prox(iz)**2d0+k3fe3/prox(iz)**3d0+k4fe3/prox(iz)**4d0)*k4fe3/prox(iz)**4d0 
    enddo 
    close(88)
endif 
            

endsubroutine calc_pH_v5

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_omega( &
    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
    & ,pco2x,cax,mgx,six,nax,prox,mineral &! input 
    & ,omega &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8),intent(in):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 
real(kind=8),dimension(nz),intent(in):: pco2x,cax,mgx,six,nax,prox
real(kind=8),dimension(nz),intent(out):: omega
character(2),intent(in):: mineral

select case(trim(adjustl(mineral)))
    case('fo')
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
        omega = mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo
        ! omega = mgx**2d0/(prox+k1mg+k1mgco3*k1*k2*kco2*pco2x/prox+k1mghco3*k1*k2*kco2*pco2x)**2d0 & 
            ! & *six/(prox**2d0+k1si*prox+k2si)/keqfo
    case('ab')
        omega = nax*six**2d0/prox/(1d0+k1si/prox+k2si/prox**2d0)**2d0/keqab
    case('an')
        omega = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & /(prox**2d0)/keqan
        ! omega = cax/(prox**2d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x+k1cahco3*k1*k2*kco2*pco2x*prox)/keqan
    case('cc')
        omega = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *k1*k2*kco2*pco2x/(prox**2d0)/keqcc
        ! omega = cax/(prox**2d0/(k1*k2*kco2*pco2x)+k1ca/prox/(k1*k2*kco2*pco2x)+k1caco3+k1cahco3*prox)/keqcc
endselect

endsubroutine calc_omega

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_omega_v2( &
    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
    & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,mineral &! input 
    & ,omega &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8),intent(in):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka
real(kind=8),dimension(nz),intent(in):: pco2x,cax,mgx,six,nax,prox,alx
real(kind=8),dimension(nz),intent(out):: omega
character(5),intent(in):: mineral

select case(trim(adjustl(mineral)))
    case('fo')
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
        omega = mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo
        ! omega = mgx**2d0/(prox+k1mg+k1mgco3*k1*k2*kco2*pco2x/prox+k1mghco3*k1*k2*kco2*pco2x)**2d0 & 
            ! & *six/(prox**2d0+k1si*prox+k2si)/keqfo
    case('ab')
    ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
        omega = nax*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0/prox**4d0/keqab
    case('an')
    ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
        omega = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**8d0/keqan
    case('cc')
        omega = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *k1*k2*kco2*pco2x/(prox**2d0)/keqcc
        ! omega = cax/(prox**2d0/(k1*k2*kco2*pco2x)+k1ca/prox/(k1*k2*kco2*pco2x)+k1caco3+k1cahco3*prox)/keqcc
    case('ka')
    ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
        omega = &
            & alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**6d0/keqka
endselect

endsubroutine calc_omega_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_omega_v3( &
    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
    & ,prox,mineral &! input 
    & ,omega &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka,keqgb,keqct,k1fe2,k1fe2co3,k1fe2hco3,keqfa,k1fe3,k2fe3,k3fe3,k4fe3,keqgt &
    & ,keqcabd,keqdp,keqhb,keqkfs
real(kind=8),dimension(nz),intent(in):: prox
real(kind=8),dimension(nz):: pco2x,cax,mgx,six,nax,alx,po2x,fe2x,fe3x,kx
real(kind=8),dimension(nz),intent(out):: omega
character(5),intent(in):: mineral

integer,intent(in)::nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst
character(5),dimension(nsp_aq),intent(in)::chraq
character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c
real(kind=8),dimension(nsp_sld_all),intent(in)::keqsld_all

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)
keqab = keqsld_all(findloc(chrsld_all,'ab',dim=1))
keqfo = keqsld_all(findloc(chrsld_all,'fo',dim=1))
keqan = keqsld_all(findloc(chrsld_all,'an',dim=1))
keqcc = keqsld_all(findloc(chrsld_all,'cc',dim=1))
k1si = keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h1)
k2si = keqaq_h(findloc(chraq_all,'si',dim=1),ieqaq_h2)
k1mg = keqaq_h(findloc(chraq_all,'mg',dim=1),ieqaq_h1)
k1mgco3 = keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_co3)
k1mghco3  = keqaq_c(findloc(chraq_all,'mg',dim=1),ieqaq_hco3)
k1ca = keqaq_h(findloc(chraq_all,'ca',dim=1),ieqaq_h1)
k1caco3 = keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_co3)
k1cahco3 = keqaq_c(findloc(chraq_all,'ca',dim=1),ieqaq_hco3)
k1al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h1)
k2al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h2)
k3al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h3)
k4al= keqaq_h(findloc(chraq_all,'al',dim=1),ieqaq_h4)
keqka = keqsld_all(findloc(chrsld_all,'ka',dim=1))
keqgb =  keqsld_all(findloc(chrsld_all,'gb',dim=1))
keqct =  keqsld_all(findloc(chrsld_all,'ct',dim=1))
k1fe2 = keqaq_h(findloc(chraq_all,'fe2',dim=1),ieqaq_h1)
k1fe2co3 = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_co3)
k1fe2hco3  = keqaq_c(findloc(chraq_all,'fe2',dim=1),ieqaq_hco3)
keqfa = keqsld_all(findloc(chrsld_all,'fa',dim=1))
k1fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h1)
k2fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h2)
k3fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h3)
k4fe3= keqaq_h(findloc(chraq_all,'fe3',dim=1),ieqaq_h4)
keqgt = keqsld_all(findloc(chrsld_all,'gt',dim=1))
keqcabd = keqsld_all(findloc(chrsld_all,'cabd',dim=1))
keqdp = keqsld_all(findloc(chrsld_all,'dp',dim=1))
keqhb = keqsld_all(findloc(chrsld_all,'hb',dim=1))
keqkfs = keqsld_all(findloc(chrsld_all,'kfs',dim=1))

nax = 0d0
kx = 0d0

if (any(chraq=='na')) then 
    nax = maqx(findloc(chraq,'na',dim=1),:)
elseif (any(chraq_cnst=='na')) then 
    nax = maqc(findloc(chraq_cnst,'na',dim=1),:)
endif 
if (any(chraq=='k')) then 
    kx = maqx(findloc(chraq,'k',dim=1),:)
elseif (any(chraq_cnst=='k')) then 
    kx = maqc(findloc(chraq_cnst,'k',dim=1),:)
endif 

six =0d0
cax =0d0
mgx =0d0
fe2x =0d0
fe3x =0d0
alx =0d0
pco2x =0d0

if (any(chraq=='si')) then 
    six = maqx(findloc(chraq,'si',dim=1),:)
elseif (any(chraq_cnst=='si')) then 
    six = maqc(findloc(chraq_cnst,'si',dim=1),:)
endif 
if (any(chraq=='ca')) then 
    cax = maqx(findloc(chraq,'ca',dim=1),:)
elseif (any(chraq_cnst=='ca')) then 
    cax = maqc(findloc(chraq_cnst,'ca',dim=1),:)
endif 
if (any(chraq=='mg')) then 
    mgx = maqx(findloc(chraq,'mg',dim=1),:)
elseif (any(chraq_cnst=='mg')) then 
    mgx = maqc(findloc(chraq_cnst,'mg',dim=1),:)
endif 
if (any(chraq=='fe2')) then 
    fe2x = maqx(findloc(chraq,'fe2',dim=1),:)
elseif (any(chraq_cnst=='fe2')) then 
    fe2x = maqc(findloc(chraq_cnst,'fe2',dim=1),:)
endif 
if (any(chraq=='al')) then 
    alx = maqx(findloc(chraq,'al',dim=1),:)
elseif (any(chraq_cnst=='al')) then 
    alx = maqc(findloc(chraq_cnst,'al',dim=1),:)
endif 
if (any(chraq=='fe3')) then 
    fe3x = maqx(findloc(chraq,'fe3',dim=1),:)
elseif (any(chraq_cnst=='fe3')) then 
    fe3x = maqc(findloc(chraq_cnst,'fe3',dim=1),:)
endif 
if (any(chrgas=='pco2')) then 
    pco2x = mgasx(findloc(chrgas,'pco2',dim=1),:)
elseif (any(chrgas_cnst=='pco2')) then 
    pco2x = mgasc(findloc(chrgas_cnst,'pco2',dim=1),:)
endif 
if (any(chrgas=='po2')) then 
    po2x = mgasx(findloc(chrgas,'po2',dim=1),:)
elseif (any(chrgas_cnst=='po2')) then 
    po2x = mgasc(findloc(chrgas_cnst,'po2',dim=1),:)
endif 

select case(trim(adjustl(mineral)))
    case('fo')
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
        omega = & 
            & mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo
        ! omega = mgx**2d0/(prox+k1mg+k1mgco3*k1*k2*kco2*pco2x/prox+k1mghco3*k1*k2*kco2*pco2x)**2d0 & 
            ! & *six/(prox**2d0+k1si*prox+k2si)/keqfo
    case('fa')
    ! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O 
        omega = & 
            & fe2x**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa
    case('ab')
    ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
        omega = & 
            & nax*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0/prox**4d0/keqab
    case('kfs')
    ! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
        omega = & 
            & kx &
            & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
            & /prox**4d0 &
            & /keqkfs 
    case('an')
    ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
        omega = & 
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**8d0/keqan
    case('cc')
        omega = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *k1*k2*kco2*pco2x/(prox**2d0)/keqcc
        ! omega = cax/(prox**2d0/(k1*k2*kco2*pco2x)+k1ca/prox/(k1*k2*kco2*pco2x)+k1caco3+k1cahco3*prox)/keqcc
    case('ka')
    ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
        omega = &
            & alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**6d0/keqka
    case('gb')
    ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O 
        omega = &
            & alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & /prox**3d0/keqgb
    case('gt')
    !  Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
        omega = &
            & fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
            & /prox**3d0/keqgt
    case('ct')
    ! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
        omega = &
            & mgx**3d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
            & /prox**6d0/keqct
    case('cabd')
    ! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
        omega = &
            & cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
            & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
            & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
            & /prox**(22d0/3d0)/keqcabd
    case('dp')
    ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
        omega = &
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
            & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
            & /prox**(4d0)/keqdp
    case('hb')
    ! Hedenbergite  + 4 H+  = 2 H2O  + 2 SiO2(aq)  + Fe++  + Ca++
        omega = &
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
            & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
            & /prox**(4d0)/keqhb
    case('py')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        omega = 1d0 - po2x**0.5d0
endselect

endsubroutine calc_omega_v3

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_rxn_ext( &
    & nz,vmax,mo2,po2th,po2x,rxn_name &! input 
    & ,rxn_ext &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8),intent(in):: vmax,mo2,po2th
real(kind=8),dimension(nz),intent(in):: po2x
real(kind=8),dimension(nz),intent(out):: rxn_ext
character(5),intent(in)::rxn_name

select case(trim(adjustl(rxn_name)))
    case('resp')
        rxn_ext = vmax*po2x/(po2x+mo2)
        ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
endselect

endsubroutine calc_rxn_ext

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_rxn_ext_v2( &
    & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
    & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
    & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
    & ,rxn_name &! input 
    & ,rxn_ext &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: po2th,fe2th
real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa
real(kind=8),dimension(nz),intent(in):: poro,sat
real(kind=8),dimension(nz),intent(out):: rxn_ext
character(5),intent(in)::rxn_name

integer,intent(in)::nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst

character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
character(5),dimension(nsp_aq),intent(in)::chraq
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst

real(kind=8),dimension(nsp_aq,nz),intent(in)::maqx
real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgasx
real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))


po2x = 0d0
if (any(chrgas=='po2')) then 
    po2x = mgasx(findloc(chrgas,'po2',dim=1),:)
elseif (any(chrgas_cnst=='po2')) then 
    po2x = mgasc(findloc(chrgas_cnst,'po2',dim=1),:)
endif 

fe2x = 0d0
if (any(chraq=='fe2')) then 
    fe2x = maqx(findloc(chraq,'fe2',dim=1),:)
elseif (any(chraq_cnst=='fe2')) then 
    fe2x = maqc(findloc(chraq_cnst,'fe2',dim=1),:)
endif 

select case(trim(adjustl(rxn_name)))
    case('resp')
        rxn_ext = vmax*po2x/(po2x+mo2)
        ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
    case('fe2o2')
        rxn_ext = poro*sat*1d3*koxa*fe2x*po2x &
            & *merge(0d0,1d0,po2x < po2th .or. fe2x < fe2th)
endselect

endsubroutine calc_rxn_ext_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
function k_arrhenius(kref,tempkref,tempk,eapp,rg)
implicit none
real(kind=8) k_arrhenius,kref,tempkref,tempk,eapp,rg
k_arrhenius = kref*exp(-eapp/rg*(1d0/tempk-1d0/tempkref))
endfunction k_arrhenius
!ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
