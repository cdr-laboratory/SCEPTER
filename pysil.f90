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
real(kind=8) ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain

call get_variables_num( &
    & nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! output
    & )

print *,nsp_sld,nsp_aq,nsp_gas,nrxn_ext

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
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain &! output
    & )
    
call weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,ttot,plant_rain  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name &! input
    & )

contains 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
subroutine weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,ttot,plant_rain  &! input
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
real(kind=8) :: mvom = 5.3d0 ! cm3/mol; molar volume of OM (CH2O); Lasaga and Ohmoto 2002
real(kind=8) :: mvomb = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8) :: mvg1 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8) :: mvg2 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8) :: mvg3 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8) :: mvamsi = 25.739d0 ! cm3/mol; molar volume of amorphous silica taken as cristobalite (SiO2); Robie et al. 1978

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
real(kind=8) :: mwtom = 30d0 ! g/mol; formula weight of CH2O
real(kind=8) :: mwtomb = 30d0 ! g/mol; formula weight of CH2O
real(kind=8) :: mwtg1 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8) :: mwtg2 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8) :: mwtg3 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8) :: mwtamsi = 60.085d0 ! g/mol; formula weight of amorphous silica

! real(kind=8) :: redsldi = 0.56d0 ! wt%  **default 
! real(kind=8) :: redsldi = 1.12d0 ! wt%  x2
real(kind=8) :: redsldi = 2.8d0 ! wt%   x5
! real(kind=8) :: redsldi = 2.24d0 ! wt%  x4
! real(kind=8) :: redsldi = 3.36d0 ! wt%  x6

! real(kind=8) :: silwti = 30d0 ! wt%  **default
! real(kind=8) :: silwti = 45d0 ! wt%  
! real(kind=8) :: silwti = 24d0 ! wt%
real(kind=8) :: silwti = 1d-10 ! wt%

real(kind=8) :: rho_grain = 2.7d0 ! g/cm3 as soil grain density 

! real(kind=8)::plant_rain = 1.4d-3 ! g C/g soil/yr; converted from 1.6d-4 mg C / g soil /hr from Georgiou et al. 2017 ! 
real(kind=8),intent(in)::plant_rain != 1d2 ! 1 t/ha/yr; approximate values from Vanveen et al. 1991 ! 
! real(kind=8)::plant_rain = 0.1d2 ! 

real(kind=8)::zsupp_plant = 0.3d0 !  e-folding decrease

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
    & ,pro,prox,ca,cax,porox,dporodta,dporodtg,dporodtgc,khco2  &
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

integer  iflx
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

logical :: timestep_fixed = .false.
! logical :: timestep_fixed = .true.

! logical :: no_biot = .false.
logical :: no_biot = .true.

logical :: biot_turbo2 = .false.
! logical :: biot_turbo2 = .true.

logical :: biot_labs = .false.
! logical :: biot_labs = .true.

! logical :: display = .false.
logical :: display = .true.

! logical :: method_precalc = .false.
logical :: method_precalc = .true.

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
integer count_dtunchanged, count_dtunchanged_Max  

integer,intent(in)::nsp_sld != 5
integer,parameter::nsp_sld_2 = 7
integer,parameter::nsp_sld_all = 20
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
integer,parameter::nrxn_ext_all = 7
integer :: nflx ! = 5 + nrxn_ext + nsp_sld 
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
real(kind=8),dimension(nsp_sld)::msldi,msldth,mv,rfrc_sld,mwt,rfrc_sld_plant
real(kind=8),dimension(nsp_sld,nsp_aq)::staq
real(kind=8),dimension(nsp_sld,nsp_gas)::stgas
real(kind=8),dimension(nsp_sld,nz)::msldx,msld,ksld,omega,msldsupp,nonprec
real(kind=8),dimension(nsp_sld,5 + nrxn_ext + nsp_sld,nz)::flx_sld
real(kind=8),dimension(nsp_aq)::maqi,maqth,daq
real(kind=8),dimension(nsp_aq,nz)::maqx,maq,rxnaq,maqsupp
real(kind=8),dimension(nsp_aq,5 + nrxn_ext + nsp_sld,nz)::flx_aq
real(kind=8),dimension(nsp_gas)::mgasi,mgasth,dgasa,dgasg,dmgas,khgasi,dgasi
real(kind=8),dimension(nsp_gas,nz)::mgasx,mgas,khgasx,khgas,dgas,agasx,agas,rxngas,mgassupp 
real(kind=8),dimension(nsp_gas,5 + nrxn_ext + nsp_sld,nz)::flx_gas 
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
real(kind=8),dimension(nsp_sld_all)::keqsld_all,mv_all,msldi_all,msldth_all,rfrc_sld_all,mwt_all,rfrc_sld_plant_all
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
integer ::isldprof2 != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 1
integer ::iaqprof != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 2
integer ::igasprof != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 3
integer ::isldsat != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 4
integer ::ibsd != ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 5

logical,dimension(nsp_sld)::turbo2,labs,nonlocal,nobio
real(kind=8),dimension(nz,nz,nsp_sld)::trans
real(kind=8) zml_ref,dbl_ref
integer izml

real(kind=8) dt_prev

logical print_cb,ph_error
character(500) print_loc
character(500),intent(in):: sim_name

real(kind=8) def_dust,def_rain,def_pr,def_OM_frc
character(5),dimension(5 + nrxn_ext + nsp_sld)::chrflx
character(10) chrfmt

integer::itflx,iadv,idif,irain,ires
data itflx,iadv,idif,irain/1,2,3,4/

integer,dimension(nsp_sld)::irxn_sld 
integer,dimension(nrxn_ext)::irxn_ext 
!-------------------------

nsp_sld_cnst = nsp_sld_all - nsp_sld
nsp_aq_cnst = nsp_aq_all - nsp_aq
nsp_gas_cnst = nsp_gas_all - nsp_gas
nsp3 = nsp_sld + nsp_aq + nsp_gas

isldprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 1
isldprof2 = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 2
iaqprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 3
igasprof = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 4
isldsat = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 5
ibsd = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 6

nflx = 5 + nrxn_ext + nsp_sld

do isps=1,nsp_sld
    irxn_sld(isps) = 4+isps
enddo 

do irxn=1,nrxn_ext
    irxn_ext(irxn) = 4+nsp_sld+irxn
enddo 

ires = nflx

chrflx(1:4) = (/'tflx ','adv  ','dif  ','rain '/)
if (nrxn_ext > 0) chrflx(irxn_sld(:)) = chrsld
if (nsp_sld > 0) chrflx(irxn_ext(:)) = chrrxn_ext
chrflx(nflx) = 'res  '

! print *,chrflx

! pause

! define all species and rxns definable in the model 
! note that rxns here exclude diss(/prec) of mineral 
! which are automatically included when associated mineral is chosen

chrsld_all = (/'fo   ','ab   ','an   ','cc   ','ka   ','gb   ','py   ','ct   ','fa   ','gt   ','cabd ' &
    & ,'dp   ','hb   ','kfs  ','om   ','omb  ','amsi ','g1   ','g2   ','g3   '/)
chraq_all = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    '/)
chrgas_all = (/'pco2','po2 '/)
chrrxn_ext_all = (/'resp ','fe2o2','omomb','ombto','g1dec','g2dec','g3dec'/)

! define the species and rxns explicitly simulated in the model in a fully coupled way
! should be chosen from definable species & rxn lists above 

! chrsld = (/'fo   ','ab   ','an   ','cc   ','ka   '/)
! chraq = (/'mg   ','si   ','na   ','ca   ','al   '/)
! chrgas = (/'pco2 ','po2  '/)
! chrrxn_ext = (/'resp '/)


! define solid species which can precipitate
! in default, all minerals only dissolve 
! should be chosen from the chrsld list
chrsld_2 = (/'cc   ','ka   ','gb   ','ct   ','gt   ','cabd ','amsi '/) 

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

mv_all = (/mvfo,mvab,mvan,mvcc,mvka,mvgb,mvpy,mvct,mvfa,mvgt,mvcabd,mvdp,mvhb,mvkfs,mvom,mvomb,mvamsi,mvg1,mvg2,mvg3/)
mwt_all = (/mwtfo,mwtab,mwtan,mwtcc,mwtka,mwtgb,mwtpy,mwtct,mwtfa,mwtgt,mwtcabd,mwtdp,mwthb,mwtkfs,mwtom,mwtomb,mwtamsi &
    & ,mwtg1,mwtg2,mwtg3/)

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

! print*,maqi_all 
! print*,mgasi_all 
print*,msldi_all

! pause

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
    print *, chrsld(isps)
    if (any(chrsld_all == chrsld(isps))) then 
        msldi(isps) = msldi_all(findloc(chrsld_all,chrsld(isps),dim=1))
        msldth(isps) = msldth_all(findloc(chrsld_all,chrsld(isps),dim=1))
        print *,msldi(isps),msldi_all(findloc(chrsld_all,chrsld(isps),dim=1))
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

! print*,maqi 
! print*,mgasi
print*,msldi

! pause

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
! Amorphous silica; SiO2
staq_all(findloc(chrsld_all,'amsi',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
! OMs; CH2O
! stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
! stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

! stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
! stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

! stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
! stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
! the above need to be modified to enable anoxic degradation 

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
! SOC assimilation by microbes 
stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'om',dim=1)) = -1d0
stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'omb',dim=1)) = 0.31d0
stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 0.69d0
stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'po2',dim=1)) = -0.69d0
! turnover of microbes 
stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'om',dim=1)) = 1d0
stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'omb',dim=1)) = -1d0
! OMs; CH2O -- G1
stsld_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1), findloc(chrsld_all,'g1',dim=1)) = -1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
! OMs; CH2O -- G2
stsld_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1), findloc(chrsld_all,'g2',dim=1)) = -1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
! OMs; CH2O -- G3
stsld_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1), findloc(chrsld_all,'g3',dim=1)) = -1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

! define 1 when a reaction is sensitive to a speces 
stgas_dext_all = 0d0
staq_dext_all = 0d0
stsld_dext_all = 0d0
! respiration 
stgas_dext_all(findloc(chrrxn_ext_all,'resp',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
! fe2 oxidation 
stgas_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
staq_dext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
! SOC assimilation by microbes 
stsld_dext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'om',dim=1)) = 1d0
stsld_dext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'omb',dim=1)) = 1d0
! turnover of microbes 
stsld_dext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'omb',dim=1)) = 1d0
! OMs; CH2O -- G1
stsld_dext_all(findloc(chrrxn_ext_all,'g1dec',dim=1), findloc(chrsld_all,'g1',dim=1)) = 1d0
stgas_dext_all(findloc(chrrxn_ext_all,'g1dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
! OMs; CH2O -- G2
stsld_dext_all(findloc(chrrxn_ext_all,'g2dec',dim=1), findloc(chrsld_all,'g2',dim=1)) = 1d0
stgas_dext_all(findloc(chrrxn_ext_all,'g2dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0
! OMs; CH2O -- G3
stsld_dext_all(findloc(chrrxn_ext_all,'g3dec',dim=1), findloc(chrsld_all,'g3',dim=1)) = 1d0
stgas_dext_all(findloc(chrrxn_ext_all,'g3dec',dim=1), findloc(chrgas_all,'po2',dim=1)) = 1d0

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
! rfrc_sld_all = rfrc_sld_all/sum(rfrc_sld_all)




! rfrc_sld_plant_all(findloc(chrsld_all,'om',dim=1)) = 1d0

! rfrc_sld_plant_all(findloc(chrsld_all,'g1',dim=1)) = 0.1d0
! rfrc_sld_plant_all(findloc(chrsld_all,'g2',dim=1)) = 0.8d0
! rfrc_sld_plant_all(findloc(chrsld_all,'g3',dim=1)) = 0.1d0

def_OM_frc = 0d0

call get_OM_rain( &
    & nsp_sld_all,chrsld_all,def_OM_frc &! input
    & ,rfrc_sld_plant_all &! output
    & )

rfrc_sld_plant_all = rfrc_sld_plant_all/mwt_all
rfrc_sld_plant_all = rfrc_sld_plant_all/sum(rfrc_sld_plant_all)



do isps = 1, nsp_sld 
    rfrc_sld(isps) = rfrc_sld_all(findloc(chrsld_all,chrsld(isps),dim=1))
    rfrc_sld_plant(isps) = rfrc_sld_plant_all(findloc(chrsld_all,chrsld(isps),dim=1))
enddo


call get_switches( &
    & no_biot,biot_turbo2,biot_labs,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc &! inout
    & )


do while (rectime(nrec)>ttot) 
    rectime = rectime/10d0
enddo 
do while (rectime(nrec)<ttot) 
    rectime = rectime*10d0
enddo 

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

write(runname,*) trim(adjustl(base))//'_q-'//trim(adjustl(chrq(3)))//'_zsat-'  &
    & //trim(adjustl(chrz(3)))

do isps = 1, nsp_sld 
    isldflx(isps) = ibasaltrain + isps
enddo 
    
do ispa = 1, nsp_aq 
    iaqflx(ispa) = ibasaltrain + nsp_sld  + ispa
enddo 

do ispg = 1, nsp_gas
    igasflx(ispg) = ibasaltrain + nsp_sld + nsp_aq + ispg
enddo 

! print*,workdir
! print*,runname
! pause

call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))

write(chrfmt,'(i0)') nflx+1

chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a))'

do isps = 1,nsp_sld
    open(isldflx(isps), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_sld-'//trim(adjustl(chrsld(isps)))//'.txt', status='replace')
    write(isldflx(isps),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
    close(isldflx(isps))
enddo 

do ispa = 1,nsp_aq
    open(iaqflx(ispa), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_aq-'//trim(adjustl(chraq(ispa)))//'.txt', status='replace')
    write(iaqflx(ispa),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
    close(iaqflx(ispa))
enddo 

do ispg = 1,nsp_gas
    open(igasflx(ispg), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_gas-'//trim(adjustl(chrgas(ispg)))//'.txt', status='replace')
    write(igasflx(ispg),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
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

pro = 1d-5
    
call coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 

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
        
    if (display) then
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
    endif      
endif
    
call coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
    
zml_ref = 1.5d0 ! mixed layer depth [m]
dbl_ref = 0d0  
labs = .false.
turbo2 = .false.
nobio = .false.
    
if (no_biot) nobio = .true.
if (biot_turbo2) turbo2 = .true.

call make_transmx(  &
    & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref  &! input
    & ,trans,nonlocal,izml  &! output 
    & )
    
! --------- loop -----
print *, 'about to start time loop'
it = 0
irec = 0

count_dtunchanged = 0

!! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@

do while (it<nt)
    call cpu_time(time_start)
    if (display) then 
        print *, 'it, time = ',it, time
    endif
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
    
    if (it == 0) then 
        if (.not.initial_ss) then 
            maxdt = 1d2
        else if (initial_ss) then
            maxdt = 0.2d0
        endif 
    endif 
    
    if (timestep_fixed) then 
        if (.not.initial_ss) then 
            maxdt = 1d2
            maxdt = 1d1
        ! else if (initial_ss .and. it==0) then
        else if (initial_ss) then
            maxdt = 0.2d0
            ! maxdt = 0.02d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.005d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.002d0 ! working with p80 = 10 um
            ! maxdt = 0.001d0 ! when calcite is included smaller time step must be assumed 
            ! maxdt = 0.0005d0 ! working with p80 = 1 um
            
            ! if (time<1d-2) then  
                ! maxdt = 1d-6 
            ! elseif (time>=1d-2 .and. time<1d-1) then 
            if (time<1d-3) then  
                maxdt = 1d-7 
            elseif (time>=1d-3 .and. time<1d-2) then  
                maxdt = 1d-6 
            elseif (time>=1d-2 .and. time<1d-1) then  
                maxdt = 1d-5 
            elseif (time>=1d-1 .and. time<1d0) then  
                maxdt = 1d-4 
            elseif (time>=1d0 .and. time<1d1) then 
                maxdt = 1d-3 
            elseif (time>=1d1 .and. time<1d2) then 
                maxdt = 1d-2  
            elseif (time>=1d2 .and. time<1d3) then 
                maxdt = 1d-1 
            elseif (time>=1d3 .and. time<1d4) then 
                maxdt = 1d0 
            elseif (time>=1d4 .and. time<1d5) then 
                maxdt = 1d1 
            elseif (time>=1d5 ) then 
                maxdt = 1d2 
            endif 
            
            ! maxdt = maxdt * 1d-1
        endif 
    endif 
    
    count_dtunchanged_Max = 1000
    ! if (dt<1d-5) then 
        ! count_dtunchanged_Max = 10
    ! elseif (dt>=1d-5 .and. dt<1d0) then
        ! count_dtunchanged_Max = 100
    ! elseif (dt>=1d0 ) then 
        ! count_dtunchanged_Max = 1000
    ! endif 

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
        
    call coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all &!input
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
    
    ! nobio = .true.
    call make_transmx(  &
        & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref  &! input
        & ,trans,nonlocal,izml  &! output 
        & )


    error = 1d4
    iter=0

100 print *, iter,time

    mgasx = mgas
    msldx = msld
    maqx = maq
    
    prox = pro  

    porox = poro


    if (initial_ss) then 
        
        maqsupp = 0d0
        mgassupp = 0d0
        do isps = 1, nsp_sld
            if (no_biot) then 
                msldsupp(isps,:) = rainpowder*rfrc_sld(isps)*exp(-z/zsupp)/zsupp
            else 
                msldsupp(isps,1) = rainpowder*rfrc_sld(isps)
            endif 
        enddo 
        
        
        if (rain_wave) then 
            do isps = 1, nsp_sld
                if (no_biot) then 
                    msldsupp(isps,:) = msldsupp(isps,:)*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                else 
                    msldsupp(isps,1) = msldsupp(isps,1)*merge(2d0,0d0,nint(time/wave_tau)==floor(time/wave_tau))
                endif 
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
        
        do isps = 1, nsp_sld
            if (no_biot) then 
                ! msldsupp(isps,:) = msldsupp(isps,:) &
                    ! & + plant_rain/12d0/((1d0-poroi)*rho_grain*1d6) &! converting g_C/g_soil/yr to mol_C/m3_soil/yr
                    ! & *1d0 &! assuming 1m depth to which plant C is supplied 
                    ! & *rfrc_sld_plant(isps) &
                    ! & *exp(-z/zsupp_plant)/zsupp_plant
                msldsupp(isps,:) = msldsupp(isps,:) &
                    & + plant_rain/12d0*rfrc_sld_plant(isps)*exp(-z/zsupp_plant)/zsupp_plant ! when plant_
            else 
                msldsupp(isps,1) = msldsupp(isps,1) &
                    & + plant_rain/12d0*rfrc_sld_plant(isps) ! when plant_rain is in g_C/m2/yr
            endif 
        enddo 
        
    else 
        mgassupp = 0d0
        msldsupp = 0d0
        maqsupp = 0d0
    endif 

    if (it==0) pre_calc = .true.
    if (method_precalc) pre_calc = .true.
    
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
            & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg,stgas_ext &! input
            & ,mgasx &! output 
            & )
            
        ! call precalc_gases_v2( &
            ! & nz,dt,ucv,dz,poro,sat,torg,tora,v,prox,hr &! input 
            ! & ,nsp_gas,nsp_gas_all,chrgas,chrgas_all,keqgas_h,mgasi,mgasth,mgas &! input
            ! & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg,stgas_ext &! input
            ! & ,nsp_sld,stgas,mv,ksld,msld,omega,nonprec &! input
            ! & ,mgasx &! output 
            ! & )
        
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
        
        ! call precalc_slds_v2( &
            ! & nz,dt,w,dz,poro,hr,sat &! input
            ! & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
            ! & ,nrxn_ext,rxnext,stsld_ext &!input
            ! & ,msldx &! output
            ! & )
            
        call precalc_slds_v2_1( &
            & nz,dt,w,dz,poro,hr,sat &! input
            & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
            & ,nrxn_ext,rxnext,stsld_ext &!input
            & ,labs,turbo2,trans &! input
            & ,msldx &! output
            & )
            
        ! call precalc_slds_v3( &
            ! & nz,dt,w,dz,poro,hr,sat &! input
            ! & ,nsp_sld,msldth,msldi,mv,msld,msldsupp,ksld,omega,nonprec &! input
            ! & ,msldx &! output
            ! & )
            
        ! call precalc_slds_v3_1( &
            ! & nz,dt,w,dz,poro,hr,sat &! input
            ! & ,nsp_sld,msldth,msldi,mv,msld,msldsupp,ksld,omega,nonprec &! input
            ! & ,labs,turbo2,trans &! input
            ! & ,msldx &! output
            ! & )

        ! pause
        
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
            
        ! call precalc_aqs_v2( &
            ! & nz,dt,v,dz,tora,poro,sat,hr &! input 
            ! & ,nsp_aq,nsp_sld,daq,maqth,maqi,maq,mv,msldx,ksld,staq,omega,nonprec &! input
            ! & ,nrxn_ext,staq_ext,rxnext &! input
            ! & ,maqx &! output
            ! & )

        if (any(isnan(mgasx)).or.any(isnan(msldx)).or.any(isnan(maqx))) then 
            print*, 'error in precalc'
            stop
        endif

    end if

    if ((.not.read_data) .and. it == 0 .and. iter == 0) then 
        do ispa = 1, nsp_aq
            if (chraq(ispa)/='so4') then
                maqx(ispa,1:) = 1d2
            endif 
        enddo
    endif 

    call alsilicate_aq_gas_1D_v3_1( &
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
        & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain &
        & ,turbo2,labs,trans,method_precalc,display &! input
        !  old inputs
        & ,hr,poro,z,dz,w,sat,pro,poroprev,tora,v,tol,it,nflx,kw & 
        & ,ucv,torg,cplprec  &
        ! old inout
        & ,iter,error,dt,flgback &    
        ! output 
        & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec & 
        & )

    ! if (iter > 75) then
        ! maxdt = maxdt/2d0
    ! end if
    ! if (iter<5) then 
        ! maxdt = maxdt*2d0
        ! if (maxdt >1d2) maxdt = 1d2
    ! endif 
    ! print*,iter,maxdt
    
    
    
    ! nobio = .false.
    ! nobio = .true.
    ! do isps =1,nsp_sld
        ! if (rfrc_sld(isps)>0d0 .or. rfrc_sld_plant(isps)>0d0) then 
            ! nobio(isps) = .false.
        ! else 
            ! nobio(isps) = .true.
        ! endif 
    ! enddo
    ! call make_transmx(  &
        ! & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref  &! input
        ! & ,trans,nonlocal,izml  &! output 
        ! & )
    ! print *,trans(:,:,1)
    ! print *
    ! print *,trans(:,:,2)
    
    ! call sld_biomix_sep( &
        ! & nz,nsp_sld,chrsld,turbo2,labs,trans,msldth &! input
        ! & ,poro,z,dz,tol,nflx,izml,nobio &! input 
        ! & ,dt,flgback,flx_sld &! inout    
        ! & ,msldx &! output 
        ! & )
    

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
    if (display) then 
        print *
        print *,'-=-=-=-=-=-= Porosity & SA -=-=-=-=-=-=-='
        print *, 'phi:', (poro(iz),iz=1,nz, nz/5)
        print *, 'SA:', (hr(iz),iz=1,nz, nz/5)
        print *
    endif 
#endif  


    if (display) then 
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
    endif 

    ! stop
    
    mgas = mgasx
    maq = maqx
    msld = msldx
    
    pro = prox

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
        open(isldprof2,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_sld(wt%)-'//chr//'.txt', status='replace')
        open(isldsat,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'sat_sld-'//chr//'.txt', status='replace')
        open(igasprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_gas-'//chr//'.txt', status='replace')
        open(iaqprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_aq-'//chr//'.txt', status='replace')
        open(ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
            & //'bsd-'//chr//'.txt', status='replace')
            
        write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(isldprof2,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(isldsat,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
        write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
        write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' time '

        do iz = 1, Nz
            write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
            write(isldprof2,*) z(iz),(msldx(isps,iz)*mwt(isps)*1d2/((1d0-poro(iz))*2.7d0*1d6),isps = 1, nsp_sld),time
            write(isldsat,*) z(iz),(omega(isps,iz),isps = 1, nsp_sld),time
            write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
            write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
            write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
        end do
        irec=irec+1

        close(isldprof)
        close(isldprof2)
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
    
    if (.not.timestep_fixed) then 
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
            
            if (count_dtunchanged > count_dtunchanged_Max) then 
                maxdt = maxdt*10d0
                count_dtunchanged = 0
            endif 
        endif 
    endif 
    
    print *,'progress_rate, maxdt, count_dtunchanged',progress_rate, maxdt, count_dtunchanged
    
    flgreducedt_prev = flgreducedt
    
end do

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
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain &! output
    & )
implicit none

integer,intent(out):: nz
real(kind=8),intent(out)::ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain
character(500),intent(out)::sim_name

character(500) file_name

file_name = './frame.in'
open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
read(50,*) ztot
read(50,*) nz
read(50,*) ttot
read(50,*) rainpowder
read(50,*) plant_rain
read(50,*) zsupp
read(50,*) poroi
read(50,*) satup
read(50,*) zsat
read(50,*) w
read(50,*) qin
read(50,*) p80
read(50,*) sim_name
close(50)

print*,'nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain'
print*,nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain

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

subroutine get_OM_rain( &
    & nsp_sld_all,chrsld_all,def_OM_frc &! input
    & ,OM_frct_all &! output
    & )
implicit none

integer,intent(in):: nsp_sld_all
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
real(kind=8),dimension(nsp_sld_all),intent(out)::OM_frct_all
real(kind=8),intent(in)::def_OM_frc 
character(5) chr_tmp
real(kind=8) val_tmp

character(500) file_name
integer i,n_tmp

file_name = './OM_rain.in'
call Console4(file_name,n_tmp)

n_tmp = n_tmp - 1

! in default 
OM_frct_all = def_OM_frc

if (n_tmp <= 0) return

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
do i =1,n_tmp
    read(50,*) chr_tmp,val_tmp
    if (any(chrsld_all == chr_tmp)) then 
        OM_frct_all(findloc(chrsld_all,chr_tmp,dim=1)) = val_tmp
    endif 
enddo 
close(50)


endsubroutine get_OM_rain

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

subroutine get_switches( &
    & no_biot,biot_turbo2,biot_labs,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc &! inout
    & )
implicit none

character(100) chr_tmp
logical,intent(inout):: no_biot,biot_turbo2,biot_labs,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc

character(500) file_name
integer i,n_tmp

file_name = './switches.in'

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')

read(50,*) no_biot,chr_tmp
read(50,*) biot_turbo2,chr_tmp
read(50,*) biot_labs,chr_tmp
read(50,*) display,chr_tmp
read(50,*) read_data,chr_tmp
read(50,*) incld_rough,chr_tmp
read(50,*) al_inhibit,chr_tmp
read(50,*) timestep_fixed,chr_tmp
read(50,*) method_precalc,chr_tmp

close(50)


endsubroutine get_switches

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

subroutine coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0
real(kind=8),dimension(nz),intent(in)::pro
real(kind=8),dimension(nz)::oh,po2
real(kind=8) kho,po2th
real(kind=8),intent(out)::ucv,kw

! real(kind=8) k_arrhenius
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

integer,intent(in)::nsp_gas,nsp_gas_cnst
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas
real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all

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
kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)

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


ksld_all(findloc(chrsld_all,'amsi',dim=1),:) = & 
    & k_arrhenius(10d0**(-12.23d0)*sec2yr,25d0+tempk_0,tc+tempk_0,74.5d0,rg) !(only neutral weathering for amsi from Palandri and Kharaka, 2004)
! SiO2 + 2 H2O = H4SiO4
keqsld_all(findloc(chrsld_all,'amsi',dim=1)) = &
    & k_arrhenius(10d0**(-2.71d0),25d0+tempk_0,tc+tempk_0,3.340d0*cal2j,rg) ! from PHREEQC.DAT 


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




po2 = 0d0
if (any(chrgas=='po2')) then 
    po2 = mgas(findloc(chrgas,'po2',dim=1),:)
elseif (any(chrgas_cnst=='po2')) then 
    po2 = mgasc(findloc(chrgas_cnst,'po2',dim=1),:)
endif 
! po2 = 1d0
po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

ksld_all(findloc(chrsld_all,'py',dim=1),:) = & 
    & k_arrhenius(10.0d0**(-8.19d0)*sec2yr,15d0+tempk_0,tc+tempk_0,57d0,rg)  &!!! excluding po2 and ph dependence
    & *(kho)**(0.50d0)/(pro**0.11d0) &! mol m^-2 yr^-1, Williamson and Rimstidt (1994)
    & *merge(0d0,1d0,po2<po2th)

!--------- other reactions -------------! 
krxn1_ext_all = 0d0
krxn2_ext_all = 0d0

krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) = & 
    & max(8.0d13*60.0d0*24.0d0*365.0d0*(kw/pro)**2.0d0, 1d-7*60.0d0*24.0d0*365.0d0) &   
    ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)excluding the term (c*po2)
    & *merge(0d0,1d0,po2<po2th)
     
krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.71d0 ! vmax mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = &
    & krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) !*1d1 ! reducing a bit to be fitted with modern soil pco2

krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
     
     
krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:) = 0.01d0*24d0*365d0 ! mg C mg-1 MBC yr-1
! converted from 0.01 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:) = 250d0 ! mg C g-1 soil  Georgiou et al. (2017)
     
     
krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:) = 0.00028d0*24d0*365d0 ! mg C mg-1 MBC yr-1
! converted from 0.00028 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:) = 2d0 ! beta value Georgiou et al. (2017)




!! OMs

krxn1_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:) = & 
    & 1d0/1d0 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1 year  
    & *merge(0d0,1d0,po2<po2th)
krxn2_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)




krxn1_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:) = & 
    & 1d0/30d0 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 30 year  
    & *merge(0d0,1d0,po2<po2th)
krxn2_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)



krxn1_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:) = & 
    & 1d0/1d3 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1k year  
    & *merge(0d0,1d0,po2<po2th)
krxn2_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:) = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)


endsubroutine coefs_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_gases( &
    & nz,dt,ucv,dz,poro,sat,torg,tora,v,prox &! input 
    & ,nsp_gas,nsp_gas_all,chrgas,chrgas_all,keqgas_h,mgasi,mgasth,mgas &! input
    & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg,stgas_ext &! input
    & ,mgasx &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::dt,ucv
real(kind=8)::pco2th,pco2i,kco2,k1,k2,kho
real(kind=8),dimension(nz),intent(in)::poro,sat,torg,tora,v,dz,prox
real(kind=8),dimension(nz)::khco2,pco2,resp,pco2x

integer iz,ispg,irxn
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
real(kind=8),dimension(nrxn_ext,nsp_gas),intent(in)::stgas_ext
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
            ! if (any(chrrxn_ext=='resp')) then 
                ! resp = rxnext(findloc(chrrxn_ext,'resp',dim=1),:)
            ! else 
                ! resp = 0d0
            ! endif 
        case('po2')
            alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3
            ! resp = 0d0
    endselect 
    
    resp = 0d0
    do irxn = 1, nrxn_ext
        if (stgas_ext(irxn,ispg)>0d0) then 
            resp = resp + stgas_ext(irxn,ispg)*rxnext(irxn,:)
        endif 
    enddo 
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg) +poro*sat*khco2*1d3*tora*dgasa(ispg)
    edifi = edif(1)
    edifi = ucv*1d3*dgasg(ispg) 

    pco2x = mgasx(ispg,:)
    pco2 = mgas(ispg,:)
    
    pco2i = mgasi(ispg)
    pco2th = mgasth(ispg)
    
    do iz = 1, nz

        ! if (pco2x(iz)>=pco2th) cycle

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

subroutine precalc_gases_v2( &
    & nz,dt,ucv,dz,poro,sat,torg,tora,v,prox,hr &! input 
    & ,nsp_gas,nsp_gas_all,chrgas,chrgas_all,keqgas_h,mgasi,mgasth,mgas &! input
    & ,nrxn_ext,chrrxn_ext,rxnext,dgasa,dgasg,stgas_ext &! input
    & ,nsp_sld,stgas,mv,ksld,msld,omega,nonprec &! input
    & ,mgasx &! output 
    & )
implicit none 
integer,intent(in)::nz
real(kind=8),intent(in)::dt,ucv
real(kind=8)::pco2th,pco2i,kco2,k1,k2,kho
real(kind=8),dimension(nz),intent(in)::poro,sat,torg,tora,v,dz,prox,hr
real(kind=8),dimension(nz)::khco2,pco2,resp,pco2x

integer iz,ispg,irxn,isps
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
real(kind=8),dimension(nrxn_ext,nsp_gas),intent(in)::stgas_ext
real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx

integer,intent(in)::nsp_sld
real(kind=8),dimension(nsp_sld),intent(in)::mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::ksld,msld,omega,nonprec
real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas

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
            ! if (any(chrrxn_ext=='resp')) then 
                ! resp = rxnext(findloc(chrrxn_ext,'resp',dim=1),:)
            ! else 
                ! resp = 0d0
            ! endif 
        case('po2')
            alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*kho*1d3
            ! resp = 0d0
    endselect 
    
    resp = 0d0
    do irxn = 1, nrxn_ext
        resp = resp + stgas_ext(irxn,ispg)*rxnext(irxn,:)
    enddo 
    
    do isps = 1,nsp_sld
        resp = resp + ( &
            & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msld(isps,:)*(1d0-omega(isps,:)) &
            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
            & )
    enddo 
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg) +poro*sat*khco2*1d3*tora*dgasa(ispg)
    edifi = edif(1)
    edifi = ucv*1d3*dgasg(ispg) 

    pco2x = mgasx(ispg,:)
    pco2 = mgas(ispg,:)
    
    pco2i = mgasi(ispg)
    pco2th = mgasth(ispg)
    
    do iz = 1, nz

        ! if (pco2x(iz)>=pco2th) cycle

        pco2tmp = pco2(max(1,iz-1))
        ediftmp = edif(max(1,iz-1))
        if (iz==1) pco2tmp = pco2i
        if (iz==1) ediftmp = edifi

        pco2x(iz) = max(pco2th*0.1d0 &
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

endsubroutine precalc_gases_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds_v2( &
    & nz,dt,w,dz,poro,hr,sat &! input
    & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
    & ,nrxn_ext,rxnext,stsld_ext &!input
    & ,msldx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w
real(kind=8)::mfoi,mfoth
real(kind=8),dimension(nz),intent(in)::dz,poro,hr,sat
real(kind=8),dimension(nz)::mfo,mfosupp,mfox,rxn_tmp

integer iz,isps,irxn

integer,intent(in)::nsp_sld,nsp_sld_2
character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_2),intent(in)::chrsld_2
real(kind=8),dimension(nsp_sld),intent(in)::msldth,msldi,mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp,ksld,omega
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx

integer,intent(in)::nrxn_ext
real(kind=8),dimension(nrxn_ext,nz),intent(in)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_sld),intent(in)::stsld_ext

if (nsp_sld == 0) return

do isps = 1,nsp_sld

    mfox = msldx(isps,:)
    mfo = msld(isps,:)
    mfoth = msldth(isps)
    mfoi = msldi(isps)
    mfosupp = msldsupp(isps,:)
    
    rxn_tmp = 0d0
    do irxn = 1, nrxn_ext
        if (stsld_ext(irxn,isps)>0d0) then 
            rxn_tmp = rxn_tmp + stsld_ext(irxn,isps)*rxnext(irxn,:)
        endif  
    enddo 
    
    if (any(chrsld_2 ==chrsld(isps))) then 
        do iz = 1, nz
            if (mfox(iz)>=mfoth) cycle

            if (iz/=nz) then 
                mfox(iz) = max(0d0, &
                    & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & +rxn_tmp(iz) &
                    & )
            else 
                mfox(iz) = max(0d0, &
                    & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & +rxn_tmp(iz) &
                    & )
            endif 
        enddo
    else
        do iz = 1, nz
            if (mfox(iz)>=mfoth) cycle

            if (iz/=nz) then 
                mfox(iz) = max(0d0, &
                    & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                    & +rxn_tmp(iz) &
                    & )
            else 
                mfox(iz) = max(0d0, &
                    & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                    & +rxn_tmp(iz) &
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

subroutine precalc_slds_v2_1( &
    & nz,dt,w,dz,poro,hr,sat &! input
    & ,nsp_sld,nsp_sld_2,chrsld,chrsld_2,msldth,msldi,mv,msld,msldsupp,ksld,omega &! input
    & ,nrxn_ext,rxnext,stsld_ext &!input
    & ,labs,turbo2,trans &! input
    & ,msldx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w
real(kind=8)::mfoi,mfoth
real(kind=8),dimension(nz),intent(in)::dz,poro,hr,sat
real(kind=8),dimension(nz)::mfo,mfosupp,mfox,rxn_tmp

integer iz,isps,irxn,iiz

integer,intent(in)::nsp_sld,nsp_sld_2
character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_2),intent(in)::chrsld_2
real(kind=8),dimension(nsp_sld),intent(in)::msldth,msldi,mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp,ksld,omega
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx

integer,intent(in)::nrxn_ext
real(kind=8),dimension(nrxn_ext,nz),intent(in)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_sld),intent(in)::stsld_ext

logical,dimension(nsp_sld),intent(in)::labs,turbo2
real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans

real(kind=8) trans_tmp(nz)

if (nsp_sld == 0) return

do isps = 1,nsp_sld

    mfox = msldx(isps,:)
    mfo = msld(isps,:)
    mfoth = msldth(isps)
    mfoi = msldi(isps)
    mfosupp = msldsupp(isps,:)
    
    rxn_tmp = 0d0
    do irxn = 1, nrxn_ext
        if (stsld_ext(irxn,isps)>0d0) then 
            rxn_tmp = rxn_tmp + stsld_ext(irxn,isps)*rxnext(irxn,:)
        endif  
    enddo 
    
    trans_tmp = 0d0
    do iz=1,nz
        do iiz=1,nz
            if (turbo2(isps) .or. labs(isps)) then 
                if (trans(iiz,iz,isps) >0d0) then 
                    trans_tmp(iz) = trans_tmp(iz)+ trans(iiz,iz,isps)/dz(iz)*dz(iiz)*mfo(iiz)
                endif 
            else
                if (trans(iiz,iz,isps) >0d0) then 
                    trans_tmp(iz) = trans_tmp(iz) + trans(iiz,iz,isps)/dz(iz)*mfo(iiz)
                endif 
            endif 
        enddo 
    enddo 
    
    if (any(chrsld_2 ==chrsld(isps))) then 
        do iz = 1, nz
            if (mfox(iz)>=mfoth) cycle

            if (iz/=nz) then 
                mfox(iz) = max(0d0, &
                    & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & +rxn_tmp(iz) &
                    & + trans_tmp(iz) &
                    & )
            else 
                mfox(iz) = max(0d0, &
                    & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                    & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(omega(isps,iz) - 1d0) &
                    & *merge(1d0,0d0,omega(isps,iz) - 1d0 > 0d0) &
                    & +rxn_tmp(iz) &
                    & + trans_tmp(iz) &
                    & )
            endif 
        enddo
    else
        do iz = 1, nz
            if (mfox(iz)>=mfoth) cycle

            if (iz/=nz) then 
                mfox(iz) = max(0d0, &
                    & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                    & +rxn_tmp(iz) &
                    & + trans_tmp(iz) &
                    & )
            else 
                mfox(iz) = max(0d0, &
                    & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                    & +rxn_tmp(iz) &
                    & + trans_tmp(iz) &
                    & )
            endif 
        enddo
    endif 
    
    msldx(isps,:) = mfox
    
enddo 

endsubroutine precalc_slds_v2_1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds_v3( &
    & nz,dt,w,dz,poro,hr,sat &! input
    & ,nsp_sld,msldth,msldi,mv,msld,msldsupp,ksld,omega,nonprec &! input
    & ,msldx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w
real(kind=8)::mfoi,mfoth
real(kind=8),dimension(nz),intent(in)::dz,poro,hr,sat
real(kind=8),dimension(nz)::mfo,mfosupp,mfox

integer iz,isps

integer,intent(in)::nsp_sld
real(kind=8),dimension(nsp_sld),intent(in)::msldth,msldi,mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp,ksld,omega,nonprec
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx

if (nsp_sld == 0) return

do isps = 1,nsp_sld

    mfox = msldx(isps,:)
    mfo = msld(isps,:)
    mfoth = msldth(isps)
    mfoi = msldi(isps)
    mfosupp = msldsupp(isps,:)
    
    do iz = 1, nz
        ! if (mfox(iz)>=mfoth) cycle

        if (iz/=nz) then 
            mfox(iz) = max(mfoth*0.1d0, &
                & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                & -ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(1d0-omega(isps,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & )
        else 
            mfox(iz) = max(mfoth*0.1d0, &
                & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                & -ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(1d0-omega(isps,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & )
        endif 
    enddo
    
    msldx(isps,:) = mfox
    
enddo 

endsubroutine precalc_slds_v3

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine precalc_slds_v3_1( &
    & nz,dt,w,dz,poro,hr,sat &! input
    & ,nsp_sld,msldth,msldi,mv,msld,msldsupp,ksld,omega,nonprec &! input
    & ,labs,turbo2,trans &! input
    & ,msldx &! output
    & )
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::dt,w
real(kind=8)::mfoi,mfoth
real(kind=8),dimension(nz),intent(in)::dz,poro,hr,sat
real(kind=8),dimension(nz)::mfo,mfosupp,mfox

integer iz,isps,iiz

integer,intent(in)::nsp_sld
real(kind=8),dimension(nsp_sld),intent(in)::msldth,msldi,mv
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp,ksld,omega,nonprec
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx
logical,dimension(nsp_sld)::labs,turbo2
real(kind=8),dimension(nz,nz,nsp_sld)::trans

real(kind=8) swnonloc

if (nsp_sld == 0) return

do isps = 1,nsp_sld

    mfox = msldx(isps,:)
    mfo = msld(isps,:)
    mfoth = msldth(isps)
    mfoi = msldi(isps)
    mfosupp = msldsupp(isps,:)
    
    swnonloc = 0d0
    if (turbo2(isps).or.labs(isps)) swnonloc = 1d0
    
    do iz = 1, nz
        ! if (mfox(iz)>=mfoth) cycle

        if (iz/=nz) then 
            mfox(iz) = max(mfoth*0.1d0, &
                & mfo(iz) +dt*(w*(mfo(iz+1)-mfo(iz))/dz(iz) + mfosupp(iz)) &
                & -ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(1d0-omega(isps,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & - sum(trans(:,iz,isps)/dz(iz)*mfo(:))*(1d0-swnonloc) &
                & - sum(trans(:,iz,isps)/dz(iz)*dz*mfo(:))*swnonloc &
                & )
        else 
            mfox(iz) = max(mfoth*0.1d0, &
                & mfo(iz) + dt*(w*(mfoi-mfo(iz))/dz(iz)+ mfosupp(iz)) &
                & -ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*mfox(iz)*(1d0-omega(isps,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & - sum(trans(:,iz,isps)/dz(iz)*mfo(:))*(1d0-swnonloc) &
                & - sum(trans(:,iz,isps)/dz(iz)*dz*mfo(:))*swnonloc &
                & )
        endif 
    enddo
    
    msldx(isps,:) = mfox
    
enddo 

endsubroutine precalc_slds_v3_1

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

subroutine precalc_aqs_v2( &
    & nz,dt,v,dz,tora,poro,sat,hr &! input 
    & ,nsp_aq,nsp_sld,daq,maqth,maqi,maq,mv,msldx,ksld,staq,omega,nonprec &! input
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
real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx,ksld,omega,nonprec
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

        ! if (nax(iz)>=nath) cycle

        ctmp = na(max(1,iz-1))
        ediftmp = edif(max(1,iz-1))
        if (iz==1) ctmp = nai
        if (iz==1) ediftmp = edifi
    
        rxn_tmp = 0d0
        
        do isps = 1, nsp_sld
            if (staq(isps,ispa)/=0d0) then 
                rxn_tmp = rxn_tmp  &
                    & + staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) 
            endif 
        enddo 
        
        do irxn = 1, nrxn_ext
            if (staq(irxn,ispa)/=0d0) then 
                rxn_tmp = rxn_tmp  &
                    & + staq_ext(irxn,ispa)*rxnext(irxn,iz)
            endif 
        enddo 
        
        nax(iz) = max(nath*0.1d0, &
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

endsubroutine precalc_aqs_v2

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
        ph_error = .true.
        return
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
    & ,keqcabd,keqdp,keqhb,keqkfs,keqamsi
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
keqamsi = keqsld_all(findloc(chrsld_all,'amsi',dim=1))

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
    case('amsi')
    ! SiO2 + 2 H2O = H4SiO4
        omega = &
            & six/(1d0+k1si/prox+k2si/prox**2d0) &
            & /keqamsi
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
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 - po2x**0.5d0*merge(0d0,1d0,po2x<1d-20)
        ! omega = 0d0
        
    case('om')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 
        ! omega = 0d0
    case('omb')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 
        ! omega = 0d0
    case default 
        omega = 1d0
endselect

if (any(isnan(omega))) then 
    print *,'nan in calc_omega_v3'
    stop
endif 

endsubroutine calc_omega_v3

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_omega_dev( &
    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
    & ,prox,mineral,sp_name &! input 
    & ,omega,domega_dmsp &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka,keqgb,keqct,k1fe2,k1fe2co3,k1fe2hco3,keqfa,k1fe3,k2fe3,k3fe3,k4fe3,keqgt &
    & ,keqcabd,keqdp,keqhb,keqkfs,keqamsi,keqg1,keqg2,keqg3
real(kind=8),dimension(nz),intent(in):: prox
real(kind=8),dimension(nz):: pco2x,cax,mgx,six,nax,alx,po2x,fe2x,fe3x,kx
real(kind=8),dimension(nz),intent(out):: domega_dmsp
real(kind=8),dimension(nz),intent(out):: omega
character(5),intent(in):: mineral,sp_name

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
keqamsi = keqsld_all(findloc(chrsld_all,'amsi',dim=1))
keqg1 = keqsld_all(findloc(chrsld_all,'g1',dim=1))
keqg2 = keqsld_all(findloc(chrsld_all,'g2',dim=1))
keqg3 = keqsld_all(findloc(chrsld_all,'g3',dim=1))

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
        omega = ( & 
            & mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo &
            & )
        ! omega = mgx**2d0/(prox+k1mg+k1mgco3*k1*k2*kco2*pco2x/prox+k1mghco3*k1*k2*kco2*pco2x)**2d0 & 
            ! & *six/(prox**2d0+k1si*prox+k2si)/keqfo
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & mgx**2d0*(-2d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                        & +k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo &
                    ! 
                    & +mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six*(-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**4d0/keqfo &
                    ! 
                    & +mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0) &
                    & *(-4d0)/prox**5d0/keqfo &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & 2d0*mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *1d0/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & mgx**2d0*(-2d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *(k1mgco3*k1*k2*kco2*1d0/prox**2d0+k1mghco3*k1*k2*kco2*1d0/prox) &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfo &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('fa')
    ! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O 
        omega = ( & 
            & fe2x**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
            & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa &
            & )
            
        ! copied and pasted from Fo case with mg changed with fe2 and fo changed with fa    
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & fe2x**2d0*(-2d0)/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                            & +k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa &
                    ! 
                    & +fe2x**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six*(-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**4d0/keqfa &
                    ! 
                    & +fe2x**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0) &
                    & *(-4d0)/prox**5d0/keqfa &
                    & )
            case('fe2')
                domega_dmsp = ( & 
                    & 2d0*fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & fe2x**2d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *1d0/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & fe2x**2d0*(-2d0)/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *(k1fe2co3*k1*k2*kco2*1d0/prox**2d0+k1fe2hco3*k1*k2*kco2*1d0/prox) &
                    & *six/(1d0+k1si/prox+k2si/prox**2d0)/prox**4d0/keqfa &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('ab')
    ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
        omega = ( & 
            & nax & 
            & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
            & /prox**4d0 &
            & /keqab &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & nax & 
                    & *alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    ! 
                    & +nax & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0*(-3d0)/(1d0+k1si/prox+k2si/prox**2d0)**4d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**4d0 &
                    & /keqab &
                    ! 
                    & +nax & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqab &
                    & )
            case('na')
                domega_dmsp = ( & 
                    & 1d0 & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & nax & 
                    & *1d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & nax & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *3d0*six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('kfs')
    ! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
        omega = ( & 
            & kx &
            & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
            & /prox**4d0 &
            & /keqkfs   &
            & )
            
        ! copied and pasted from ab case with na changed with k and ab changed with kfs  
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & kx & 
                    & *alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqkfs &
                    ! 
                    & +kx & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0*(-3d0)/(1d0+k1si/prox+k2si/prox**2d0)**4d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**4d0 &
                    & /keqkfs &
                    ! 
                    & +kx & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqkfs &
                    & )
            case('k')
                domega_dmsp = ( & 
                    & 1d0 & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqkfs &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & kx & 
                    & *1d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqkfs &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & kx & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *3d0*six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqkfs &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('an')
    ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
        omega = ( & 
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**8d0 &
            & /keqan &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    ! 
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *alx**2d0*(-2d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**3d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    ! 
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0*(-2d0)/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**8d0 &
                    & /keqan &
                    ! 
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & *(-8d0)/prox**9d0 &
                    & /keqan &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *2d0*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *1d0**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**8d0 &
                    & /keqan &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('cc')
        omega = ( &
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *k1*k2*kco2*pco2x/(prox**2d0) &
            & /keqcc &
            & )
        ! omega = cax/(prox**2d0/(k1*k2*kco2*pco2x)+k1ca/prox/(k1*k2*kco2*pco2x)+k1caco3+k1cahco3*prox)/keqcc
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ca++ dependence on pH is from 'an' case
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                                    & +k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    ! 
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *k1*k2*kco2*pco2x*(-2d0)/(prox**3d0) &
                    & /keqcc &
                    & )
            case('pco2')
                ! ca++ dependence on pCO2 is from 'an' case
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox) &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    ! 
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *k1*k2*kco2*1d0/(prox**2d0) &
                    & /keqcc &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('ka')
    ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
        omega = ( &
            & alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
            & /prox**6d0 &
            & /keqka &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependences of al and si are from an case
                domega_dmsp = ( & 
                    & alx**2d0*(-2d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**3d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    ! 
                    & +alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0*(-2d0)/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**6d0 &
                    & /keqka &
                    ! 
                    & +alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & *(-6d0)/prox**7d0 &
                    & /keqka &
                    ! 
                    & )
            case('al')
                domega_dmsp = ( & 
                    & 2d0*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & alx**2d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *six*2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('gb')
    ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O 
        omega = ( &
            & alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & /prox**3d0 &
            & /keqgb &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependence of al is from an case
                domega_dmsp = ( & 
                    & alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & /prox**3d0 &
                    & /keqgb &
                    ! 
                    & +alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *(-3d0)/prox**4d0 &
                    & /keqgb &
                    ! 
                    & )
            case('al')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & /prox**3d0 &
                    & /keqgb &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('amsi')
    ! SiO2 + 2 H2O = H4SiO4
        omega = ( &
            & six/(1d0+k1si/prox+k2si/prox**2d0) &
            & /keqamsi &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependence of si is from above
                domega_dmsp = ( & 
                    & six*(-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**2d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /keqamsi &
                    ! 
                    & )
            case('si')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1si/prox+k2si/prox**2d0) &
                    & /keqamsi &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('gt')
    !  Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
        omega = (&
            & fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
            & /prox**3d0 &
            & /keqgt &
            & )
            
        ! copied and pasted from gb case with replacing al and gb with fe3 and gt respectively
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0)**2d0 &
                    & *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0) &
                    & /prox**3d0 &
                    & /keqgt &
                    ! 
                    & +fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
                    & *(-3d0)/prox**4d0 &
                    & /keqgt &
                    ! 
                    & )
            case('fe3')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0) &
                    & /prox**3d0 &
                    & /keqgt &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('ct')
    ! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
        omega = ( &
            & mgx**3d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
            & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
            & /prox**6d0 &
            & /keqct &
            & )
            
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependence of mg and si are from above
                domega_dmsp = ( & 
                    & mgx**3d0*(-3d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**4d0 &
                    & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                                            & +k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    ! 
                    & +mgx**3d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *six**2d0*(-2d0)/(1d0+k1si/prox+k2si/prox**2d0)**3d0  &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**6d0 &
                    & /keqct &
                    ! 
                    & +mgx**3d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
                    & *(-6d0)/prox**7d0 &
                    & /keqct &
                    & )
            case('pco2')
                ! pco2 dependence of mg is from fo case
                domega_dmsp = ( & 
                    & mgx**3d0*(-3d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**4d0 &
                    & *(k1mgco3*k1*k2*kco2*1d0/prox**2d0+k1mghco3*k1*k2*kco2*1d0/prox) &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & 3d0*mgx**2d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *six**2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & mgx**3d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**3d0 &
                    & *six*2d0/(1d0+k1si/prox+k2si/prox**2d0)**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('cabd')
    ! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
        omega = ( &
            & cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
            & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
            & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
            & /prox**(22d0/3d0) &
            & /keqcabd &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &          
                    & cax**(1d0/6d0)*(-1d0/6d0) &
                            &/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0+1d0) &
                    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    !
                    & +cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *alx**(7d0/3d0)*(-7d0/3d0) &
                            & /(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0+1d0) &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    !
                    & +cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)*(-11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0+1d0) &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    !
                    & +cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & *(-22d0/3d0)/prox**(22d0/3d0-1d0) &
                    & /keqcabd &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & cax**(1d0/6d0)*(-1d0/6d0) &
                            & /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0+1d0) &
                    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & (1d0/6d0)*cax**(1d0/6d0-1d0) &
                            & /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *(7d0/3d0)*alx**(7d0/3d0-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *six**(11d0/3d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & cax**(1d0/6d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**(1d0/6d0) &
                    & *alx**(7d0/3d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**(7d0/3d0) &
                    & *(11d0/3d0)*six**(11d0/3d0-1d0)/(1d0+k1si/prox+k2si/prox**2d0)**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('dp')
    ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
        omega = ( &
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
            & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
            & /prox**(4d0) &
            & /keqdp &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)*(-2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(3d0) &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & *(-4d0)/prox**(5d0) &
                    & /keqdp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1mgco3*k1*k2*kco2*1d0/prox**2d0+k1mghco3*k1*k2*kco2*1d0/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *1d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox) &
                    & *six*(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('hb')
    ! Hedenbergite  + 4 H+  = 2 H2O  + 2 SiO2(aq)  + Fe++  + Ca++
        omega = &
            & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
            & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
            & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
            & /prox**(4d0)/keqhb
            
        ! copied and pasted from dp case with replacing mg and dp by fe2 and hb
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x*(-1d0)/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0 &
                                                        & +k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)*(-2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(3d0) &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & *(-4d0)/prox**(5d0) &
                    & /keqdp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    !
                    & +cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x*(-1d0)/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox)**2d0 &
                    & *(k1fe2co3*k1*k2*kco2*1d0/prox**2d0+k1fe2hco3*k1*k2*kco2*1d0/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('fe2')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *1d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six**(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox) &
                    & *fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox) &
                    & *six*(2d0)/(1d0+k1si/prox+k2si/prox**2d0)**(2d0) &
                    & /prox**(4d0) &
                    & /keqdp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('py')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 - po2x**0.5d0*merge(0d0,1d0,po2x<1d-20)
        ! omega = 0d0
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 0.5d0*po2x**(-0.5d0)*merge(0d0,1d0,po2x<1d-20) &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('om')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 
        domega_dmsp = 0d0
        ! omega = 0d0
        
    case('omb')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 
        domega_dmsp = 0d0
        ! omega = 0d0
        
        
    case('g1')
    ! omega is defined so that kg1*poro*hr*mvg1*1d-6*mg1x*(1d0-omega_g1) = kg1*poro*hr*mvg1*1d-6*mg1x*po2x/(po2x+mo2)
    ! i.e., 1.0 - omega_g1 = po2x/(po2x+mo2) 
        omega = 1d0 
        domega_dmsp = 0d0
        ! print *,mineral,sp_name
        ! omega = 1d0 - po2x/(po2x+mo2)
        
        
        
        ! omega = 1d0 - po2x/(po2x + keqg1)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! domega_dmsp = ( &
                    ! & - 1d0/(po2x + keqg1) &
                    ! & - po2x*(-1d0)/(po2x + keqg1)**2d0 &
                    ! & )
            ! case default 
                ! domega_dmsp = 0d0
        ! endselect 
        
        ! print*, omega
        
        
    case('g2')
    ! omega is defined so that kg2*poro*hr*mvg2*1d-6*mg2x*(1d0-omega_g2) = kg2*poro*hr*mvg2*1d-6*mg2x*po2x/(po2x+mo2)
    ! i.e., 1.0 - omega_g2 = po2x/(po2x+mo2) 
        ! omega = 1d0 - po2x/(po2x+mo2)
        ! print *,mineral,sp_name
        
        
        omega = 1d0 
        domega_dmsp = 0d0
        
        ! omega = 1d0 - po2x/(po2x + keqg2)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! domega_dmsp = ( &
                    ! & - 1d0/(po2x + keqg2) &
                    ! & - po2x*(-1d0)/(po2x + keqg2)**2d0 &
                    ! & )
            ! case default 
                ! domega_dmsp = 0d0
        ! endselect 
        
        ! print*, omega
        
        
    case('g3')
    ! omega is defined so that kg3*poro*hr*mvg3*1d-6*mg3x*(1d0-omega_g3) = kg3*poro*hr*mvg3*1d-6*mg3x*po2x/(po2x+mo2)
    ! i.e., 1.0 - omega_g3 = po2x/(po2x+mo2) 
        ! omega = 1d0 - po2x/(po2x+mo2)
        ! print *,mineral,sp_name
        omega = 1d0 
        domega_dmsp = 0d0
        
        ! omega = 1d0 - po2x/(po2x + keqg3)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! domega_dmsp = ( &
                    ! & - 1d0/(po2x + keqg3) &
                    ! & - po2x*(-1d0)/(po2x + keqg3)**2d0 &
                    ! & )
            ! case default 
                ! domega_dmsp = 0d0
        ! endselect 
        
        ! print*, omega
        
        
    case default 
        ! print *,'non-specified'
        omega = 1d0
        domega_dmsp = 0d0
        ! print *,omega
        
endselect

if (any(isnan(omega)) .or. any(isnan(domega_dmsp))) then 
    print *,'nan in calc_omega_dev'
    stop
endif 

endsubroutine calc_omega_dev

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_rxn_ext_v2( &
    & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
    & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
    & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
    & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
    & ,rxn_name &! input 
    & ,rxn_ext &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: po2th,fe2th,mwtom
real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa,vmax2,mom2,komb,beta,omx,ombx
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

integer,intent(in)::nsp_sld,nsp_sld_cnst

character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst

real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc

real(kind=8),intent(in)::rho_grain


vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))



vmax2 = krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)
mom2 = krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)

komb = krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)
beta = krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)


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

omx = 0d0
if (any(chrsld=='om')) then 
    omx = msldx(findloc(chrsld,'om',dim=1),:)
elseif (any(chraq_cnst=='om')) then 
    omx = msldc(findloc(chrsld_cnst,'om',dim=1),:)
endif 

ombx = 0d0
if (any(chrsld=='omb')) then 
    ombx = msldx(findloc(chrsld,'omb',dim=1),:)
elseif (any(chraq_cnst=='omb')) then 
    ombx = msldc(findloc(chrsld_cnst,'omb',dim=1),:)
endif 

select case(trim(adjustl(rxn_name)))
    case('resp')
        rxn_ext = vmax*po2x/(po2x+mo2)
        ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
    case('fe2o2')
        rxn_ext = poro*sat*1d3*koxa*fe2x*po2x &
            & *merge(0d0,1d0,po2x < po2th .or. fe2x < fe2th)
    case('omomb')
        rxn_ext = vmax2 & ! mg C / soil g /yr
            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
    case('ombto')
        rxn_ext = komb*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**beta &
            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
    case default 
        rxn_ext = 0d0
endselect

if (any(isnan(rxn_ext))) then 
    print *,'nan in calc_rxn_ext_v2'
    stop
endif 

endsubroutine calc_rxn_ext_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_rxn_ext_dev( &
    & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
    & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
    & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
    & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
    & ,rxn_name,sp_name &! input 
    & ,rxn_ext,drxnext_dmsp &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: po2th,fe2th,mwtom
real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa,vmax2,mom2,komb,beta,omx,ombx &
    & ,mo2g1,mo2g2,mo2g3,kg1,kg2,kg3,g1x,g2x,g3x
real(kind=8),dimension(nz),intent(in):: poro,sat
real(kind=8),dimension(nz),intent(out):: drxnext_dmsp
real(kind=8),dimension(nz),intent(out):: rxn_ext
character(5),intent(in)::rxn_name,sp_name

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

integer,intent(in)::nsp_sld,nsp_sld_cnst

character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst

real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc

real(kind=8),intent(in)::rho_grain


vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))



vmax2 = krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)
mom2 = krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)

komb = krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)
beta = krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)

kg1 = krxn1_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:)
mo2g1 = krxn2_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:)

kg2 = krxn1_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:)
mo2g2 = krxn2_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:)

kg3 = krxn1_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:)
mo2g3 = krxn2_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:)


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

omx = 0d0
if (any(chrsld=='om')) then 
    omx = msldx(findloc(chrsld,'om',dim=1),:)
elseif (any(chraq_cnst=='om')) then 
    omx = msldc(findloc(chrsld_cnst,'om',dim=1),:)
endif 

ombx = 0d0
if (any(chrsld=='omb')) then 
    ombx = msldx(findloc(chrsld,'omb',dim=1),:)
elseif (any(chraq_cnst=='omb')) then 
    ombx = msldc(findloc(chrsld_cnst,'omb',dim=1),:)
endif 

g1x = 0d0
if (any(chrsld=='g1')) then 
    g1x = msldx(findloc(chrsld,'g1',dim=1),:)
elseif (any(chraq_cnst=='g1')) then 
    g1x = msldc(findloc(chrsld_cnst,'g1',dim=1),:)
endif 

g2x = 0d0
if (any(chrsld=='g2')) then 
    g2x = msldx(findloc(chrsld,'g2',dim=1),:)
elseif (any(chraq_cnst=='g2')) then 
    g2x = msldc(findloc(chrsld_cnst,'g2',dim=1),:)
endif 

g3x = 0d0
if (any(chrsld=='g3')) then 
    g3x = msldx(findloc(chrsld,'g3',dim=1),:)
elseif (any(chraq_cnst=='g3')) then 
    g3x = msldc(findloc(chrsld_cnst,'g3',dim=1),:)
endif 

select case(trim(adjustl(rxn_name)))

    case('resp')
        rxn_ext = vmax*po2x/(po2x+mo2)
        ! rxn_ext = vmax*merge(0d0,po2x/(po2x+mo2),(po2x <po2th).or.(isnan(po2x/(po2x+mo2))))
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                drxnext_dmsp = (&
                    & vmax*1d0/(po2x+mo2) &
                    & +vmax*po2x*(-1d0)/(po2x+mo2)**2d0 &
                    & )
            case default
                drxnext_dmsp = 0d0
        endselect 
        
    case('fe2o2')
        rxn_ext = poro*sat*1d3*koxa*fe2x*po2x &
            & *merge(0d0,1d0,po2x < po2th .or. fe2x < fe2th)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                drxnext_dmsp = ( &
                    & poro*sat*1d3*koxa*fe2x*1d0 &
                    & ) &
                    & *merge(0d0,1d0,po2x < po2th .or. fe2x < fe2th)
            case('fe2')
                drxnext_dmsp = ( &
                    & poro*sat*1d3*koxa*1d0*po2x &
                    & ) &
                    & *merge(0d0,1d0,po2x < po2th .or. fe2x < fe2th)
            case default
                drxnext_dmsp = 0d0
        endselect
    
    case('omomb')
        rxn_ext = vmax2 & ! mg C / soil g /yr
            & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
            & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
            & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
        
        select case(trim(adjustl(sp_name)))
            case('om')
                drxnext_dmsp = ( &
                    & vmax2 & ! mg C / soil g /yr
                    & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                    & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                    & + vmax2 & ! mg C / soil g /yr
                    & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                    & *ombx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & *(-1d0)/(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3)**2d0 &
                    & *(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                    & ) 
            case('omb')
                drxnext_dmsp = ( &
                    & vmax2 & ! mg C / soil g /yr
                    & *omx*(1d0-poro)*rho_grain*1d6*12d0*1d3 &! mol/m3 converted to mg C/ soil g
                    & *1d0*(1d0-poro)*rho_grain*1d6*12d0*1d3 &
                    & /(mom2 + (omx*(1d0-poro)*rho_grain)*1d6*12d0*1d3) &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                    & ) 
            case default
                drxnext_dmsp = 0d0
        endselect
    
    case('ombto')
        rxn_ext = komb*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**beta &
            & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) ! converting mg_C/soil_g to mol_C/soil_m3
        
        select case(trim(adjustl(sp_name)))
            case('omb')
                drxnext_dmsp = ( &
                    & komb*beta*(ombx*(1d0-poro)*rho_grain*rho_grain*1d6)**(beta-1d0) &
                    & *(1d0-poro)*rho_grain*rho_grain*1d6 &
                    & *1d-3/12d0/((1d0-poro)*rho_grain*1d6) &! converting mg_C/soil_g to mol_C/soil_m3
                    & ) 
    
            case default
                drxnext_dmsp = 0d0
                
        endselect 
    
    case('g1dec')
        rxn_ext = kg1*g1x*po2x/(po2x+mo2g1)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                drxnext_dmsp = (&
                    & kg1*g1x*1d0/(po2x+mo2g1) &
                    & +kg1*g1x*po2x*(-1d0)/(po2x+mo2g1)**2d0 &
                    & )
            case('g1')
                drxnext_dmsp = (&
                    & kg1*1d0*po2x/(po2x+mo2g1) &
                    & )
            case default
                drxnext_dmsp = 0d0
        endselect 
    
    case('g2dec')
        rxn_ext = kg2*g2x*po2x/(po2x+mo2g2)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                drxnext_dmsp = (&
                    & kg2*g2x*1d0/(po2x+mo2g2) &
                    & +kg2*g2x*po2x*(-1d0)/(po2x+mo2g2)**2d0 &
                    & )
            case('g2')
                drxnext_dmsp = (&
                    & kg2*1d0*po2x/(po2x+mo2g2) &
                    & )
            case default
                drxnext_dmsp = 0d0
        endselect 
    
    case('g3dec')
        rxn_ext = kg3*g3x*po2x/(po2x+mo2g3)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                drxnext_dmsp = (&
                    & kg3*g3x*1d0/(po2x+mo2g3) &
                    & +kg3*g3x*po2x*(-1d0)/(po2x+mo2g3)**2d0 &
                    & )
            case('g3')
                drxnext_dmsp = (&
                    & kg3*1d0*po2x/(po2x+mo2g3) &
                    & )
            case default
                drxnext_dmsp = 0d0
        endselect 
                
                
                
    case default 
        rxn_ext = 0d0
        
endselect

if (any(isnan(rxn_ext)) .or. any(isnan(drxnext_dmsp))) then 
    print *,'nan in calc_rxn_ext_dev'
    stop
endif 

endsubroutine calc_rxn_ext_dev

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine make_transmx(  &
    & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref  &! input
    & ,trans,nonlocal,izml  &! output 
    & )
implicit none
integer,intent(in)::nsp_sld,nz
real(kind=8),intent(in)::dz(nz),poro(nz),z(nz),zml_ref,dbl_ref
real(kind=8)::sporo(nz)
logical,intent(in)::labs(nsp_sld),turbo2(nsp_sld),nobio(nsp_sld)
real(kind=8),intent(out)::trans(nz,nz,nsp_sld)
logical,intent(out)::nonlocal(nsp_sld)
integer,intent(out)::izml
integer iz,isp,iiz,izdbl
real(kind=8) :: translabs(nz,nz),dbio(nz),transdbio(nz,nz),transturbo2(nz,nz)
real(kind=8) :: zml(nsp_sld),probh,dbl

sporo = 1d0 - poro
trans = 0d0
!~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
if (any(labs)) then
    translabs = 0d0

    open(unit=88,file='../input/labs-mtx.txt',action='read',status='unknown')
    do iz=1,nz
        read(88,*) translabs(iz,:)  ! writing 
    enddo
    close(88)

endif

if (.true.) then  ! devided by the time duration when transition matrices are created in LABS and weakening by a factor
! if (.false.) then 
    ! translabs = translabs *365.25d0/10d0*1d0/3d0  
    ! translabs = translabs *365.25d0/10d0*1d0/15d0  
    translabs = translabs *365.25d0/10d0*1d0/10d0
endif
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
zml=zml_ref ! mixed layer depth assumed to be a reference value at first 

dbl = dbl_ref

nonlocal = .false. ! initial assumption 
do isp=1,nsp_sld
    if (turbo2(isp) .or. labs(isp)) nonlocal(isp)=.true. ! if mixing is made by turbo2 or labs, then nonlocal 
    
    dbio=0d0
    izdbl=0
    do iz = 1, nz
        if (z(iz) <= dbl) then 
            dbio(iz) = 0d0
            izdbl = iz
        elseif (dbl < z(iz) .and. z(iz) <=zml(isp)) then
            ! dbio(iz) =  0.15d-4   !  within mixed layer 150 cm2/kyr (Emerson, 1985) 
            dbio(iz) =  2d-4   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
            ! dbio(iz) =  2d-4*exp(z(iz)/0.1d0)   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
            ! dbio(iz) =  2d-7*exp(z(iz)/1d0)   !  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
            ! dbio(iz) =  2d-10   !  just a small value 
            izml = iz   ! determine grid of bottom of mixed layer 
        else
            dbio(iz) =  0d0 ! no biodiffusion in deeper depths 
        endif
    enddo

    transdbio = 0d0   ! transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above 
    ! do iz = max(1,izdbl), izml
        ! if (iz==max(1,izdbl)) then 
            ! transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! elseif (iz==izml) then 
            ! transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        ! else 
            ! transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                ! + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! endif
    ! enddo
    do iz = max(1,izdbl), izml
        if (iz==max(1,izdbl)) then 
            transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz))*(-1d0)/(0.5d0*(dz(iz)+dz(iz)))  &
                ! + 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        elseif (iz==izml) then 
            transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        else 
            transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                + 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        endif
    enddo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    ! transition matrix for random mixing 
    transturbo2 = 0d0
    ! ending up in upward mixing 
    probh = 0.0010d0
    transturbo2(max(1,izdbl):izml,max(1,izdbl):izml) = probh  ! arbitrary assumed probability 
    do iz=1,izml  ! when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
       transturbo2(iz,iz)=-probh*(izml-max(1,izdbl))  
    enddo
    ! trying real homogeneous 
    ! transturbo2 = 0d0
    ! probh = 0.0001d0
    ! do iz=1,izml 
        ! do iiz=1,izml
            ! if (iiz/=iz) then 
                ! transturbo2(iiz,iz) = probh*dz(iz)/dz(iiz)
                ! transturbo2(iiz,iiz) = transturbo2(iiz,iiz) - transturbo2(iiz,iz)
            ! endif 
        ! enddo
    ! enddo

    if (turbo2(isp)) translabs = transturbo2   ! translabs temporarily used to represents nonlocal mixing 

    trans(:,:,isp) = transdbio(:,:)  !  firstly assume local mixing implemented by dbio 

    if (nonlocal(isp)) trans(:,:,isp) = translabs(:,:)  ! if nonlocal, replaced by either turbo2 mixing or labs mixing 
    if (nobio(isp)) trans(:,:,isp) = 0d0  ! if assuming no bioturbation, transition matrix is set at zero  
enddo
! even when all are local Fickian mixing, mixing treatment must be the same as in case of nonlocal 
! if mixing intensity and depths are different between different species  
if (all(.not.nonlocal)) then  
    do isp=1,nsp_sld-1
        if (any(trans(:,:,isp+1)/=trans(:,:,isp))) nonlocal=.true.
    enddo
endif 

endsubroutine make_transmx

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine alsilicate_aq_gas_1D_v3_1( &
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
    & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain &
    & ,turbo2,labs,trans,method_precalc,display &! input
    !  old inputs
    & ,hr,poro,z,dz,w,sat,pro,poroprev,tora,v,tol,it,nflx,kw & 
    & ,ucv,torg,cplprec  &
    ! old inout
    & ,iter,error,dt,flgback &    
    ! output 
    & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec & 
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,tol,kw,ucv,rho_grain
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,tora,v,poroprev,dz,torg,pro
real(kind=8),dimension(nz),intent(out)::prox
integer,intent(inout)::iter,it
logical,intent(in)::cplprec,method_precalc,display
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,intent(in)::nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext
character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_2),intent(in)::chrsld_2
character(5),dimension(nsp_aq),intent(in)::chraq
character(5),dimension(nsp_aq_ph),intent(in)::chraq_ph
character(5),dimension(nsp_gas_ph),intent(in)::chrgas_ph
character(5),dimension(nsp_gas),intent(in)::chrgas
character(5),dimension(nrxn_ext),intent(in)::chrrxn_ext
real(kind=8),dimension(nsp_sld),intent(in)::msldi,msldth,mv
real(kind=8),dimension(nsp_aq),intent(in)::maqi,maqth,daq 
real(kind=8),dimension(nsp_gas),intent(in)::mgasi,mgasth,dgasa,dgasg,khgasi
real(kind=8),dimension(nsp_gas)::dgasi
real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,ksld,msldsupp 
real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx,omega,nonprec
real(kind=8),dimension(nsp_sld,nz)::domega_dpro,dmsld
real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq
real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas
real(kind=8),dimension(nsp_sld,nflx,nz),intent(out)::flx_sld
real(kind=8),dimension(nsp_aq,nz),intent(in)::maq,maqsupp
real(kind=8),dimension(nsp_aq,nz),intent(inout)::maqx 
real(kind=8),dimension(nsp_aq,nz)::dprodmaq,dmaq 
real(kind=8),dimension(nsp_aq,nflx,nz),intent(out)::flx_aq
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,mgassupp
real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx 
real(kind=8),dimension(nsp_gas,nz)::khgasx,khgas,dgas,agasx,agas,rxngas,dkhgas_dpro,dprodmgas,dmgas
real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
real(kind=8),dimension(nsp_gas,nflx,nz),intent(out)::flx_gas 
real(kind=8),dimension(nrxn_ext,nz),intent(inout)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_gas),intent(in)::stgas_ext,stgas_dext
real(kind=8),dimension(nrxn_ext,nsp_aq),intent(in)::staq_ext,staq_dext
real(kind=8),dimension(nrxn_ext,nsp_sld),intent(in)::stsld_ext,stsld_dext
real(kind=8),dimension(nrxn_ext,nsp_gas,nz)::drxnext_dmgas
real(kind=8),dimension(nrxn_ext,nsp_aq,nz)::drxnext_dmaq
real(kind=8),dimension(nrxn_ext,nsp_sld,nz)::drxnext_dmsld
logical,dimension(nsp_sld),intent(in)::labs,turbo2

integer,intent(in)::nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst,nsp_sld_cnst
character(5),dimension(nsp_aq_cnst),intent(in)::chraq_cnst
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
character(5),dimension(nsp_gas_cnst),intent(in)::chrgas_cnst
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all
real(kind=8),dimension(nsp_aq_cnst,nz),intent(in)::maqc
real(kind=8),dimension(nsp_gas_cnst,nz),intent(in)::mgasc
real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc
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

integer,intent(in)::nrxn_ext_all

character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all

real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

integer iz,row,ie,ie2,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn,isps2,iiz
integer::itflx,iadv,idif,irain,ires
data itflx,iadv,idif,irain/1,2,3,4/

integer,dimension(nsp_sld)::irxn_sld 
integer,dimension(nrxn_ext)::irxn_ext 

real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,omega_tmp_th,rxn_ext_tmp &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,khco2n_tmp,pco2n_tmp,edifn_tmp,caqsupp_tmp,kco2,k1,k2,kho,sw_red

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 3d0
real(kind=8),dimension(nz)::dummy,dummy2  

logical print_cb,ph_error
character(500) print_loc

integer,parameter :: iter_max = 50
! integer,parameter :: iter_max = 300

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

!-----------------------------------------------

sw_red = 1d0
if (.not.method_precalc) sw_red = -1d100

do isps=1,nsp_sld
    irxn_sld(isps) = 4+isps
enddo 

do irxn=1,nrxn_ext
    irxn_ext(irxn) = 4+nsp_sld+irxn
enddo 

ires = nflx

print_cb = .false. 
print_loc = './ph.txt'

kco2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h0)
k1 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h1)
k2 = keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)

kho = keqgas_h(findloc(chrgas_all,'po2',dim=1),ieqgas_h0)


! print *, staq
! print *, mgx
! print *, six
! print *, mfosupp
! stop

if (any(isnan(tora)))then 
    print*,tora
endif 

dummy = 0d0
dummy2 = 0d0

error = 1d4
iter = 0

if (it ==0) then 

    call calc_pH_v5( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqc,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
        & ,print_cb,print_loc,z &! input 
        & ,prox,ph_error &! output
        & ) 
    if (ph_error) then 
        dt = dt/10d0
        flgback = .true.
        return
    endif 
else 

    call calc_pH_v5( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
        & ,print_cb,print_loc,z &! input 
        & ,prox,ph_error &! output
        & ) 
    if (ph_error) then 
        dt = dt/10d0
        flgback = .true.
        return
    endif 
endif

! print *, 'starting silciate calculation'

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_sld = 0d0
    flx_aq = 0d0
    flx_gas = 0d0
    
    ! pH calculation and its derivative wrt aq and gas species
    
    call calc_pH_v5( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
        & ,print_cb,print_loc,z &! input 
        & ,prox,ph_error &! output
        & ) 
    if (ph_error) then 
        dt = dt/10d0
        flgback = .true.
        exit
    endif 
    do ispa=1,nsp_aq
        dmaq = 0d0
        dmgas = 0d0
        if (any (chraq_ph == chraq(ispa))) then 
            dmaq(ispa,:) = dconc
            call calc_pH_v5( &
                & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
                & ,maqx+dmaq,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
                & ,print_cb,print_loc,z &! input 
                & ,dummy,ph_error &! output
                & ) 
            dprodmaq(ispa,:) = (dummy - prox)/dconc
        endif 
    enddo 
    
    dprodmgas = 0d0
    do ispg=1,nsp_gas
        dmaq = 0d0
        dmgas = 0d0
        if (any (chrgas_ph == chrgas(ispg))) then 
            dmgas(ispg,:) = dconc
            call calc_pH_v5( &
                & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
                & ,maqx,maqc,mgasx+dmgas,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
                & ,print_cb,print_loc,z &! input 
                & ,dummy,ph_error &! output
                & ) 
            dprodmgas(ispg,:) = (dummy - prox)/dconc
        endif 
    enddo
    
    ! saturation state calc. and their derivatives wrt aq and gas species
    
    omega = 0d0
    domega_dpro = 0d0
    domega_dmaq = 0d0
    domega_dmgas = 0d0
    
    do isps =1, nsp_sld
        ! call calc_omega_v2( &
            ! & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            ! & ,mgasx(findloc(chrgas,'pco2',dim=1),:) &
            ! & ,maqx(findloc(chraq,'ca',dim=1),:) &
            ! & ,maqx(findloc(chraq,'mg',dim=1),:) &
            ! & ,maqx(findloc(chraq,'si',dim=1),:) &
            ! & ,maqx(findloc(chraq,'na',dim=1),:) &
            ! & ,prox &
            ! & ,maqx(findloc(chraq,'al',dim=1),:) &
            ! & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            ! & ,omega(isps,:) &! output
            ! & )
        ! call calc_omega_v3( &
            ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            ! & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
            ! & ,prox,chrsld(isps) &! input 
            ! & ,dummy &! output
            ! & )
        dummy = 0d0
        dummy2 = 0d0
        call calc_omega_dev( &
            & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
            & ,prox,chrsld(isps),'pro  ' &! input 
            & ,dummy,dummy2 &! output
            & )
        omega(isps,:) = dummy
        dmaq = 0d0
        dmgas = 0d0
        
        ! call calc_omega_v3( &
            ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            ! & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
            ! & ,prox+dconc,chrsld(isps) &! input 
            ! & ,domega_dpro(isps,:) &! output
            ! & )
        ! domega_dpro(isps,:) = (domega_dpro(isps,:)-omega(isps,:))/dconc
        
                
        ! call calc_omega_dev( &
            ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            ! & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
            ! & ,prox,chrsld(isps),'pro  ' &! input 
            ! & ,dummy &! output
            ! & )
        
        domega_dpro(isps,:) = dummy2
        
        do ispa = 1, nsp_aq
            dmaq = 0d0
            dmgas = 0d0
            ! if (any (chraq_ph == chraq(ispa))) then 
            if (any (chraq_ph == chraq(ispa)) .or. staq(isps,ispa)/=0d0 ) then 
                dmaq(ispa,:) = dconc
                ! call calc_omega_v3( &
                    ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    ! & ,maqx+dmaq,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
                    ! & ,prox,chrsld(isps) &! input 
                    ! & ,domega_dmaq(isps,ispa,:) &! output
                    ! & )
                ! domega_dmaq(isps,ispa,:) = (domega_dmaq(isps,ispa,:)-omega(isps,:))/dconc
                
                dummy = 0d0
                dummy2 = 0d0
                call calc_omega_dev( &
                    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
                    & ,prox,chrsld(isps),chraq(ispa) &! input 
                    & ,dummy,dummy2 &! output
                    & )
                
                domega_dmaq(isps,ispa,:) = dummy2 + domega_dpro(isps,:)*dprodmaq(ispa,:)
            endif 
        enddo
        do ispg = 1, nsp_gas
            dmaq = 0d0
            dmgas = 0d0
            ! print *,chrgas(ispg),stgas(isps,ispg)
            ! if (any (chrgas_ph == chrgas(ispg))) then 
            if (any (chrgas_ph == chrgas(ispg)) .or. stgas(isps,ispg)/=0d0) then 
                dmgas(ispg,:) = dconc
                ! print*,'derivation calc',chrsld(isps),chrgas(ispg)
                ! call calc_omega_v3( &
                    ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    ! & ,maqx,maqc,mgasx+dmgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
                    ! & ,prox,chrsld(isps) &! input 
                    ! & ,domega_dmgas(isps,ispg,:) &! output
                    ! & )
                ! domega_dmgas(isps,ispg,:) = (domega_dmgas(isps,ispg,:)-omega(isps,:))/dconc
                
                dummy = 0d0
                dummy2 = 0d0
                call calc_omega_dev( &
                    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
                    & ,prox,chrsld(isps),chrgas(ispg) &! input 
                    & ,dummy,dummy2 &! output
                    & )
                
                domega_dmgas(isps,ispg,:) = dummy2 + domega_dpro(isps,:)*dprodmgas(ispg,:)
            endif 
        enddo
    enddo 
    
    nonprec = 1d0 ! primary minerals only dissolve
    if (cplprec)then
        do isps = 1, nsp_sld
            if (any(chrsld_2 == chrsld(isps))) then  
                nonprec(isps,:) = 0d0 ! allowing precipitation for secondary phases
            endif 
        enddo
    endif 
    
    ! adding reactions that are not based on dis/prec of minerals
    rxnext = 0d0
    drxnext_dmaq = 0d0
    drxnext_dmgas = 0d0
    drxnext_dmsld = 0d0
    
    do irxn=1,nrxn_ext
        ! call calc_rxn_ext_v2( &
            ! & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
            ! & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
            ! & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
            ! & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
            ! & ,chrrxn_ext(irxn) &! input 
            ! & ,dummy &! output
            ! & )
        dummy = 0d0
        dummy2 = 0d0
        call calc_rxn_ext_dev( &
            & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
            & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
            & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
            & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
            & ,chrrxn_ext(irxn),'pro  ' &! input 
            & ,dummy,dummy2 &! output
            & )
        rxnext(irxn,:) = dummy
        dmgas = 0d0
        do ispg=1,nsp_gas
            if (stgas_dext(irxn,ispg)==0d0) cycle
            dmgas(ispg,:) = dconc
            
            ! call calc_rxn_ext_v2( &
                ! & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                ! & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                ! & ,poro,sat,maqx,maqc,mgasx+dmgas,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                ! & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
                ! & ,chrrxn_ext(irxn) &! input 
                ! & ,drxnext_dmgas(irxn,ispg,:) &! output
                ! & )
            ! drxnext_dmgas(irxn,ispg,:) = (drxnext_dmgas(irxn,ispg,:) -rxnext(irxn,:))/dconc
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
                & ,chrrxn_ext(irxn),chrgas(ispg) &! input 
                & ,dummy,dummy2 &! output
                & )
            drxnext_dmgas(irxn,ispg,:) = dummy2
        enddo 
        dmaq = 0d0
        do ispa=1,nsp_aq
            if (staq_dext(irxn,ispa)==0d0) cycle
            dmaq(ispa,:) = dconc
            
            ! call calc_rxn_ext_v2( &
                ! & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                ! & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                ! & ,poro,sat,maqx+dmaq,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                ! & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
                ! & ,chrrxn_ext(irxn) &! input 
                ! & ,drxnext_dmaq(irxn,ispa,:) &! output
                ! & )
            ! drxnext_dmaq(irxn,ispa,:) = (drxnext_dmaq(irxn,ispa,:) -rxnext(irxn,:))/dconc
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
                & ,chrrxn_ext(irxn),chraq(ispa) &! input 
                & ,dummy,dummy2 &! output
                & )
            drxnext_dmaq(irxn,ispa,:) = dummy2
        enddo 
        dmsld = 0d0
        do isps=1,nsp_sld
            if (stsld_dext(irxn,isps)==0d0) cycle
            dmsld(isps,:) = dconc
            
            ! call calc_rxn_ext_v2( &
                ! & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                ! & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                ! & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                ! & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx+dmsld,msldc,rho_grain &!input
                ! & ,chrrxn_ext(irxn) &! input 
                ! & ,drxnext_dmsld(irxn,isps,:) &! output
                ! & )
            ! drxnext_dmsld(irxn,isps,:) = (drxnext_dmsld(irxn,isps,:) -rxnext(irxn,:))/dconc
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
                & ,chrrxn_ext(irxn),chrsld(isps) &! input 
                & ,dummy,dummy2 &! output
                & )
            drxnext_dmsld(irxn,isps,:) = dummy2
        enddo 
    enddo 
            
    

    do iz = 1, nz  !================================
        
        do isps = 1, nsp_sld
        
            row = nsp3*(iz-1)+isps
            
            k_tmp = ksld(isps,iz)
            mv_tmp = mv(isps)
            omega_tmp = omega(isps,iz)
            omega_tmp_th = omega_tmp*nonprec(isps,iz)
            m_tmp = msldx(isps,iz)
            mth_tmp = msldth(isps) 
            mi_tmp = msldi(isps)
            mp_tmp = msldx(isps,min(nz,iz+1))
            msupp_tmp = msldsupp(isps,iz) 
            rxn_ext_tmp = sum(stsld_ext(:,isps)*rxnext(:,iz))
            mprev_tmp = msld(isps,iz)
            
            if (iz==nz) mp_tmp = mi_tmp

            amx3(row,row) = (1.0d0     &
                & + w/dz(iz)*dt    &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp)*dt &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps,iz))*dt &
                & ) &
                & * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red)

            ymx3(row) = ( &
                & (m_tmp-mprev_tmp) &
                & -w*(mp_tmp-m_tmp)/dz(iz)*dt  &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp)*dt &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & -msupp_tmp*dt  &
                & -rxn_ext_tmp*dt  &
                & ) &
                & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                
            if (iz/=nz) amx3(row,row+nsp3) = (-w/dz(iz))*dt *merge(1.0d0,mp_tmp,m_tmp<mth_tmp*sw_red)
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1)+nsp_sld + ispa
                
                amx3(row,col ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmaq(isps,ispa,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & - sum(stsld_ext(:,isps)*drxnext_dmaq(:,ispa,iz))*dt &
                    & ) &
                    & *maqx(ispa,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
            enddo 
            
            do ispg = 1, nsp_gas 
                col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                amx3(row,col) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmgas(isps,ispg,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & - sum(stsld_ext(:,isps)*drxnext_dmgas(:,ispg,iz))*dt &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
            enddo 
            
            do isps2 = 1,nsp_sld 
                if (isps2 == isps) cycle
                col = nsp3*(iz-1)+ isps2

                amx3(row,col) = ( &
                    & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps2,iz))*dt &
                    & ) &
                    & *msldx(isps2,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
            enddo 
            
            ! diffusion terms are filled with transition matrices 
            if (turbo2(isps).or.labs(isps)) then
                do iiz = 1, nz
                    col = nsp3*(iiz-1)+isps
                    ! col = 1 + (iiz-1)*nsp
                    if (trans(iiz,iz,isps)==0d0) cycle
                    amx3(row,col) = amx3(row,col) &
                        ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*msldx(iiz,isps)
                        & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)*dt
                    ymx3(row) = ymx3(row) &
                        ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*(1d0-poro(iiz))*msldx(iiz,isps)
                        & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)*dt
                        
                    flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                        & )
                enddo
            else
                do iiz = 1, nz
                    ! col = 1 + (iiz-1)*nsp
                    col = nsp3*(iiz-1)+isps
                    if (trans(iiz,iz,isps)==0d0) cycle
                    ! amx3(row,col) = amx3(row,col) -trans(iiz,iz,isps)/dz(iz)*msldx(iiz,isps)/(1d0-poro(iiz))
                    ! ymx3(row) = ymx3(row) - trans(iiz,iz,isps)/dz(iz)*msldx(iiz,isps)/(1d0-poro(iiz))
                        
                    ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        ! & - trans(iiz,iz,isps)/dz(iz)*msldx(iiz,isps)/(1d0-poro(iiz)) &
                        ! & )
                        
                    amx3(row,col) = amx3(row,col) -trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz)*dt &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    ymx3(row) = ymx3(row) - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz)*dt &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        
                    flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                        & )
                enddo
            endif
            
            flx_sld(isps,itflx,iz) = ( &
                & (m_tmp-mprev_tmp)/dt &
                & )
            flx_sld(isps,iadv,iz) = ( &
                & -w*(mp_tmp-m_tmp)/dz(iz)  &
                & )
            flx_sld(isps,irxn_sld(isps),iz) = ( &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & )
            flx_sld(isps,irain,iz) = (&
                & - msupp_tmp  &
                & )
            flx_sld(isps,irxn_ext(:),iz) = (&
                    & - stsld_ext(:,isps)*rxnext(:,iz)  &
                    & )
            flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz))
            if (isnan(flx_sld(isps,ires,iz))) then 
                print *,chrsld(isps),iz,(flx_sld(isps,iflx,iz),iflx=1,nflx)
            endif   
        enddo 
    end do  !================================

    do iz = 1, nz
        
        do ispa = 1, nsp_aq

            row = nsp3*(iz-1)+ nsp_sld + ispa
            
            d_tmp = daq(ispa)
            caq_tmp = maqx(ispa,iz)
            caq_tmp_prev = maq(ispa,iz)
            caq_tmp_p = maqx(ispa,min(nz,iz+1))
            caq_tmp_n = maqx(ispa,max(1,iz-1))
            caqth_tmp = maqth(ispa)
            caqi_tmp = maqi(ispa)
            caqsupp_tmp = maqsupp(ispa,iz) 
            rxn_ext_tmp = sum(staq_ext(:,ispa)*rxnext(:,iz))
            rxn_tmp = 0d0
            drxndisp_tmp = 0d0
            do isps = 1,nsp_sld
                ! if (staq(isps,ispa)==0d0) cycle
                rxn_tmp = rxn_tmp + ( &
                    & staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                    & ) 
                drxndisp_tmp = drxndisp_tmp + ( &
                    & staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,ispa,iz)) &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                    & ) 
            enddo
            
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
                & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa,iz))*dt &
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp*sw_red)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*dt &
                & - rxn_tmp*dt &
                & - caqsupp_tmp*dt &
                & - rxn_ext_tmp*dt &
                & ) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)

            if (iz/=1) then 
                amx3(row,row-nsp3) = ( &
                    & -(-0.5d0*(edif_tmp +edif_tmp_n)*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(-1d0)/dz(iz)*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            endif 
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            endif 
            
            do isps = 1, nsp_sld
                col = nsp3*(iz-1)+ isps
                
                amx3(row, col) = (     & 
                    & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                    & - sum(staq_ext(:,ispa)*drxnext_dmsld(:,isps,iz))*dt &
                    & ) &
                    & *msldx(isps,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            enddo 
            
            do ispa2 = 1, nsp_aq
                col = nsp3*(iz-1)+ nsp_sld + ispa2
                
                if (ispa2 == ispa) cycle
                
                amx3(row,col) = amx3(row,col) + (     & 
                    & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa2,iz))*dt &
                    & ) &
                    & *maqx(ispa2,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmaq(isps,ispa2,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *maqx(ispa2,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                enddo 
            enddo 
            
            do ispg = 1, nsp_gas
                col = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg
                
                amx3(row,col) = amx3(row,col) + (     & 
                    & - sum(staq_ext(:,ispa)*drxnext_dmgas(:,ispg,iz))*dt &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmgas(isps,ispg,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *mgasx(ispg,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                enddo 
            enddo 
                    
            flx_aq(ispa,itflx,iz) = (&
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                & ) 
            flx_aq(ispa,iadv,iz) = (&
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                & ) 
            flx_aq(ispa,idif,iz) = (&
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                & ) 
            flx_aq(ispa,irxn_sld(:),iz) = (& 
                & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
                & ) 
            flx_aq(ispa,irain,iz) = (&
                & - caqsupp_tmp &
                & ) 
            flx_aq(ispa,irxn_ext(:),iz) = (&
                & - staq_ext(:,ispa)*rxnext(:,iz) &
                & ) 
            flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz))
            if (isnan(flx_aq(ispa,ires,iz))) then 
                print *,chraq(ispa),iz,(flx_aq(ispa,iflx,iz),iflx=1,nflx)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/(poro(iz)*sat(iz)*1d3)
            ! ymx3(row) = ymx3(row)/(poro(iz)*sat(iz)*1d3)
        
        enddo 
        
    end do  ! ==============================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2 & pO2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    khgas = 0d0
    khgasx = 0d0
    dkhgas_dmaq = 0d0
    dkhgas_dmgas = 0d0
    
    dgas = 0d0
    ddgas_dmaq = 0d0
    ddgas_dmgas = 0d0
    
    agas = 0d0
    agasx = 0d0
    dagas_dmaq = 0d0
    dagas_dmgas = 0d0
    
    rxngas = 0d0
    drxngas_dmaq = 0d0
    drxngas_dmsld = 0d0
    drxngas_dmgas = 0d0
    
    do ispg = 1, nsp_gas
        
        select case (trim(adjustl(chrgas(ispg))))
            case('pco2')
                khgas(ispg,:) = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
                khgasx(ispg,:) = kco2*(1d0+k1/prox + k1*k2/prox/prox)
        
                dkhgas_dpro(ispg,:) = kco2*(k1*(-1d0)/prox**2d0 + k1*k2*(-2d0)/prox**3d0)
            case('po2')
                khgas(ispg,:) = kho ! previous value; should not change through iterations 
                khgasx(ispg,:) = kho
        
                dkhgas_dpro(ispg,:) = 0d0
        endselect 
        
        dgas(ispg,:) = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg)+poro*sat*khgasx(ispg,:)*1d3*tora*dgasa(ispg)
        dgasi(ispg) = ucv*1d3*dgasg(ispg) 
        
        agas(ispg,:)= ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khgas(ispg,:)*1d3
        agasx(ispg,:)= ucv*poro*(1.0d0-sat)*1d3+poro*sat*khgasx(ispg,:)*1d3
        
        do ispa = 1,nsp_aq 
            dkhgas_dmaq(ispg,ispa,:) = dkhgas_dpro(ispg,:)*dprodmaq(ispa,:)
            ddgas_dmaq(ispg,ispa,:) = poro*sat*dkhgas_dmaq(ispg,ispa,:)*1d3*tora*dgasa(ispg)
            dagas_dmaq(ispg,ispa,:) =  poro*sat*dkhgas_dmaq(ispg,ispa,:)*1d3
        enddo 
        
        do ispg2 = 1,nsp_gas 
            dkhgas_dmgas(ispg,ispg2,:) = dkhgas_dpro(ispg,:)*dprodmgas(ispg2,:)
            ddgas_dmgas(ispg,ispg2,:) = poro*sat*dkhgas_dmgas(ispg,ispg2,:)*1d3*tora*dgasa(ispg)
            dagas_dmgas(ispg,ispg2,:) =  poro*sat*dkhgas_dmgas(ispg,ispg2,:)*1d3
        enddo 
        
        do isps = 1, nsp_sld
            rxngas(ispg,:) =  rxngas(ispg,:) + (&
                & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
            drxngas_dmsld(ispg,isps,:) =  drxngas_dmsld(ispg,isps,:) + (&
                & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
            do ispg2 = 1,nsp_gas
                drxngas_dmgas(ispg,ispg2,:) =  drxngas_dmgas(ispg,ispg2,:) + (&
                    & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg2,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
            do ispa = 1,nsp_aq
                drxngas_dmaq(ispg,ispa,:) =  drxngas_dmaq(ispg,ispa,:) + ( &
                    & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
        enddo 
    enddo 
    
    ! print *,drxngas_dmaq(findloc(chrgas,'pco2',dim=1),findloc(chraq,'ca',dim=1),:)
    
    do iz = 1, nz
        
        do ispg = 1, nsp_gas
        
            row = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg
            
            pco2n_tmp = mgasx(ispg,max(1,iz-1))
            khco2n_tmp = khgasx(ispg,max(1,iz-1))
            edifn_tmp = dgas(ispg,max(1,iz-1))
            if (iz == 1) then 
                pco2n_tmp = mgasi(ispg)
                khco2n_tmp = khgasi(ispg)
                edifn_tmp = dgasi(ispg)
            endif 

            amx3(row,row) = ( &
                & (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz)) &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz)*dt &
                & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg,iz))*dt &
                & -drxngas_dmgas(ispg,ispg,iz)*dt &
                & ) &
                & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
            
            ymx3(row) = ( &
                & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz)) &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*dt &
                ! & -resp(iz)*dt &
                & -sum(stgas_ext(:,ispg)*rxnext(:,iz))*dt &
                & -rxngas(ispg,iz)*dt &
                & -mgassupp(ispg,iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg)*sw_red)
            
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                        & + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 
                    amx3(row,col+nsp3) = ( &
                            & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                            & ) &
                            & *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                enddo 
            
            endif 
            
            if (iz/=1) then 
                amx3(row,row-nsp3) = ( &
                    & -(- 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                    & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                    &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-khgasx(ispg,iz-1)*1d0)/dz(iz)*dt &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg,iz-1)*mgasx(ispg,iz-1))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,mgasx(ispg,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 

                    amx3(row,col-nsp3) = ( &
                        & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                        &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                enddo 
            endif 
            
            do isps = 1,nsp_sld
                col = nsp3*(iz-1) + isps 
                amx3(row,col) = ( &
                    & -drxngas_dmsld(ispg,isps,iz)*dt &
                    & -sum(stgas_ext(:,ispg)*drxnext_dmsld(:,isps,iz))*dt &
                    & ) &
                    & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
            enddo 
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1) + nsp_sld + ispa 
                amx3(row,col) = ( &
                    & (dagas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz)) &
                    & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & - 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/dz(iz)*dt &
                    & -drxngas_dmaq(ispg,ispa,iz)*dt &
                    & -sum(stgas_ext(:,ispg)*drxnext_dmaq(:,ispa,iz))*dt &
                    & ) &
                    & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
            enddo 
            
            do ispg2 = 1, nsp_gas
                if (ispg == ispg2) cycle
                col = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg2
                amx3(row,col) = ( &
                    & (dagas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz)) &
                    & -( 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & - 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,iz)-pco2n_tmp) &
                    &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/dz(iz)*dt &
                    & -drxngas_dmgas(ispg,ispg2,iz)*dt &
                    & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg2,iz))*dt &
                    & ) &
                    & *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
            enddo 
            
            flx_gas(ispg,itflx,iz) = ( &
                & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt &
                & )         
            flx_gas(ispg,idif,iz) = ( &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & )
            flx_gas(ispg,iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) &
                & )
            flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
            flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz)
            flx_gas(ispg,irxn_sld(:),iz) = &
                & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) 
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
            
            if (any(isnan(flx_gas(ispg,:,iz)))) then
                print *,flx_gas(ispg,:,iz)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/alpha(iz)
            ! ymx3(row) = ymx3(row)/alpha(iz)
        enddo 

    end do 
    
    ! amx3 = amx3/maxval(ksld)
    ! ymx3 = ymx3/maxval(ksld)
    
    ymx3=-1.0d0*ymx3

    if (any(isnan(amx3)).or.any(isnan(ymx3)).or.any(amx3>infinity).or.any(ymx3>infinity)) then 
    ! if (.true.) then 
        print*,'error in mtx'
        print*,'any(isnan(amx3)),any(isnan(ymx3))'
        print*,any(isnan(amx3)),any(isnan(ymx3))

        ! if (any(isnan(ymx3))) then 
            ! do ie = 1,nsp3*(nz)
                ! if (isnan(ymx3(ie))) then 
                    ! print*,'NAN is here...',ie
                ! endif
            ! enddo
        ! endif


        ! if (any(isnan(amx3))) then 
            ! do ie = 1,nsp3*(nz)
                ! do ie2 = 1,nsp3*(nz)
                    ! if (isnan(amx3(ie,ie2))) then 
                        ! print*,'NAN is here...',ie,ie2
                    ! endif
                ! enddo
            ! enddo
        ! endif
        
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
        do isps = 1, nsp_sld
            row = isps + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chrsld(isps)
                stop
            endif

            if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                msldx(isps,iz) = msldx(isps,iz)*1.5d0
            else if (ymx3(row) < -threshold) then 
                msldx(isps,iz) = msldx(isps,iz)*0.50d0
            else   
                msldx(isps,iz) = msldx(isps,iz)*exp(ymx3(row))
            endif
        enddo 
        
        do ispa = 1, nsp_aq
            row = ispa + nsp_sld + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chraq(ispa)
                stop
            endif

            if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                maqx(ispa,iz) = maqx(ispa,iz)*1.5d0
            else if (ymx3(row) < -threshold) then 
                maqx(ispa,iz) = maqx(ispa,iz)*0.50d0
            else   
                maqx(ispa,iz) = maqx(ispa,iz)*exp(ymx3(row))
            endif
        enddo 
        
        do ispg = 1, nsp_gas
            row = ispg + nsp_aq + nsp_sld + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chrgas(ispg)
                stop
            endif

            if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                mgasx(ispg,iz) = mgasx(ispg,iz)*1.5d0
            else if (ymx3(row) < -threshold) then 
                mgasx(ispg,iz) = mgasx(ispg,iz)*0.50d0
            else   
                mgasx(ispg,iz) = mgasx(ispg,iz)*exp(ymx3(row))
            endif
        enddo 
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(msldx)) .or. any(isnan(maqx)).or. any(isnan(mgasx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(msldx)),any(isnan(maqx)),any(isnan(mgasx))
        
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
    ! call calc_pH_v4( &
        ! & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        ! & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        ! & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v5( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
        & ,print_cb,print_loc,z &! input 
        & ,prox,ph_error &! output
        & ) 

    if (display) then 
        print *, 'silicate_dis error',error,info,iter,dt
    endif      
    iter = iter + 1 

    ! if (iter > iter_Max ) then
    if (iter > iter_Max .or. error > infinity) then
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
        
        exit 
    end if
    
#ifdef dispiter
    print *
    print *,'-=-=-=-=-=-= Aq species flx -=-=-=-=-=-=-='
    do ispa = 1, nsp_aq
        print *, trim(adjustl(chraq(ispa))), (sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
    print *,'-=-=-=-=-=-= Sld species flx -=-=-=-=-=-=-='
    do isps = 1, nsp_sld
        print *, trim(adjustl(chrsld(isps))), (sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
    do isps = 1, nsp_sld
        print *, 'omega_'//trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Gas species flx -=-=-=-=-=-=-='
    do ispg = 1, nsp_gas
        print *, trim(adjustl(chrgas(ispg))), (sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
#endif     

enddo

! just addint flx calculation at the end 

    
flx_sld = 0d0
flx_aq = 0d0
flx_gas = 0d0

! pH calculation and its derivative wrt aq and gas species

! call calc_pH_v4( &
    ! & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
    ! & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
    ! & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
    ! & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
    ! & ,prox &! output
    ! & ) 
call calc_pH_v5( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c &! input
    & ,print_cb,print_loc,z &! input 
    & ,prox,ph_error &! output
    & ) 

! saturation state calc. and their derivatives wrt aq and gas species

omega = 0d0

do isps =1, nsp_sld
    ! call calc_omega_v3( &
        ! & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        ! & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
        ! & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
        ! & ,prox,chrsld(isps) &! input 
        ! & ,omega(isps,:) &! output
        ! & )
    dummy = 0d0
    dummy2 = 0d0
    call calc_omega_dev( &
        & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all &! input
        & ,prox,chrsld(isps),'pro  ' &! input 
        & ,dummy,dummy2 &! output
        & )
    omega(isps,:) = dummy
enddo 

nonprec = 1d0 ! primary minerals only dissolve
if (cplprec)then
    do isps = 1, nsp_sld
        if (any(chrsld_2 == chrsld(isps))) then  
            nonprec(isps,:) = 0d0 ! allowing precipitation for secondary phases
        endif 
    enddo
endif 

! adding reactions that are not based on dis/prec of minerals
rxnext = 0d0

do irxn=1,nrxn_ext
    ! call calc_rxn_ext_v2( &
        ! & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
        ! & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
        ! & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
        ! & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
        ! & ,chrrxn_ext(irxn) &! input 
        ! & ,rxnext(irxn,:) &! output
        ! & )
    dummy = 0d0
    dummy2 = 0d0
    call calc_rxn_ext_dev( &
        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain &!input
        & ,chrrxn_ext(irxn),'pro  ' &! input 
        & ,dummy,dummy2 &! output
        & )
    rxnext(irxn,:) = dummy
enddo 


do iz = 1, nz  !================================
    
    do isps = 1, nsp_sld
        
        k_tmp = ksld(isps,iz)
        mv_tmp = mv(isps)
        omega_tmp = omega(isps,iz)
        omega_tmp_th = omega_tmp*nonprec(isps,iz)
        m_tmp = msldx(isps,iz)
        mth_tmp = msldth(isps) 
        mi_tmp = msldi(isps)
        mp_tmp = msldx(isps,min(nz,iz+1))
        msupp_tmp = msldsupp(isps,iz)
        rxn_ext_tmp = sum(stsld_ext(:,isps)*rxnext(:,iz))
        mprev_tmp = msld(isps,iz)
        
        if (iz==nz) mp_tmp = mi_tmp
        
        ! diffusion terms are filled with transition matrices 
        if (turbo2(isps).or.labs(isps)) then
            do iiz = 1, nz
                if (trans(iiz,iz,isps)==0d0) cycle
                    
                flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                    & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                    & )
            enddo
        else
            do iiz = 1, nz
                if (trans(iiz,iz,isps)==0d0) cycle
                    
                flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                    ! & - trans(iiz,iz,isps)/dz(iz)*msldx(iiz,isps)/(1d0-poro(iiz)) &
                    & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                    & )
            enddo
        endif
        
        flx_sld(isps,itflx,iz) = ( &
            & (m_tmp-mprev_tmp)/dt &
            & )
        flx_sld(isps,iadv,iz) = ( &
            & -w*(mp_tmp-m_tmp)/dz(iz)  &
            & )
        flx_sld(isps,irxn_sld(isps),iz) = ( &
            & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
            & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
            & )
        flx_sld(isps,irain,iz) = (&
            & - msupp_tmp  &
            & )
        flx_sld(isps,irxn_ext(:),iz) = (&
                & - stsld_ext(:,isps)*rxnext(:,iz)  &
                & )
        flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz))
        if (isnan(flx_sld(isps,ires,iz))) then 
            print *,chrsld(isps),iz,(flx_sld(isps,iflx,iz),iflx=1,nflx)
        endif   
    enddo 
end do  !================================

do iz = 1, nz
    
    do ispa = 1, nsp_aq
        
        d_tmp = daq(ispa)
        caq_tmp = maqx(ispa,iz)
        caq_tmp_prev = maq(ispa,iz)
        caq_tmp_p = maqx(ispa,min(nz,iz+1))
        caq_tmp_n = maqx(ispa,max(1,iz-1))
        caqth_tmp = maqth(ispa)
        caqi_tmp = maqi(ispa)
        caqsupp_tmp = maqsupp(ispa,iz)
        rxn_ext_tmp = sum(staq_ext(:,ispa)*rxnext(:,iz))
        rxn_tmp = 0d0
        drxndisp_tmp = 0d0
        do isps = 1,nsp_sld
            ! if (staq(isps,ispa)==0d0) cycle
            rxn_tmp = rxn_tmp + ( &
                & staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & ) 
            drxndisp_tmp = drxndisp_tmp + ( &
                & staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,ispa,iz)) &
                & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) &
                & ) 
        enddo
        
        if (iz==1) caq_tmp_n = caqi_tmp
            
        edif_tmp = 1d3*poro(iz)*sat(iz)*tora(iz)*d_tmp
        edif_tmp_p = 1d3*poro(min(iz+1,nz))*sat(min(iz+1,nz))*tora(min(iz+1,nz))*d_tmp
        edif_tmp_n = 1d3*poro(max(iz-1,1))*sat(max(iz-1,1))*tora(max(iz-1,1))*d_tmp
                
        flx_aq(ispa,itflx,iz) = (&
            & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
            & ) 
        flx_aq(ispa,iadv,iz) = (&
            & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
            & ) 
        flx_aq(ispa,idif,iz) = (&
            & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
            & ) 
        flx_aq(ispa,irxn_sld(:),iz) = (& 
            & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
            & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
            & ) 
        flx_aq(ispa,irain,iz) = (&
            & - caqsupp_tmp &
            & ) 
        flx_aq(ispa,irxn_ext(:),iz) = (&
            & - staq_ext(:,ispa)*rxnext(:,iz) &
            & ) 
        flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz))
        if (isnan(flx_aq(ispa,ires,iz))) then 
            print *,chraq(ispa),iz,(flx_aq(ispa,iflx,iz),iflx=1,nflx)
        endif 
    
    enddo 
    
end do  ! ==============================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2 & pO2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

khgas = 0d0
khgasx = 0d0

dgas = 0d0

agas = 0d0
agasx = 0d0

rxngas = 0d0

do ispg = 1, nsp_gas
    
    select case (trim(adjustl(chrgas(ispg))))
        case('pco2')
            khgas(ispg,:) = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
            khgasx(ispg,:) = kco2*(1d0+k1/prox + k1*k2/prox/prox)
        case('po2')
            khgas(ispg,:) = kho ! previous value; should not change through iterations 
            khgasx(ispg,:) = kho
    endselect 
    
    dgas(ispg,:) = ucv*poro*(1.0d0-sat)*1d3*torg*dgasg(ispg)+poro*sat*khgasx(ispg,:)*1d3*tora*dgasa(ispg)
    dgasi(ispg) = ucv*1d3*dgasg(ispg) 
    
    agas(ispg,:)= ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khgas(ispg,:)*1d3
    agasx(ispg,:)= ucv*poro*(1.0d0-sat)*1d3+poro*sat*khgasx(ispg,:)*1d3
    
    do isps = 1, nsp_sld
        rxngas(ispg,:) =  rxngas(ispg,:) + (&
            & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
            & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
            & )
    enddo 
enddo 

! print *,drxngas_dmaq(findloc(chrgas,'pco2',dim=1),findloc(chraq,'ca',dim=1),:)

do iz = 1, nz
    
    do ispg = 1, nsp_gas
        
        pco2n_tmp = mgasx(ispg,max(1,iz-1))
        khco2n_tmp = khgasx(ispg,max(1,iz-1))
        edifn_tmp = dgas(ispg,max(1,iz-1))
        if (iz == 1) then 
            pco2n_tmp = mgasi(ispg)
            khco2n_tmp = khgasi(ispg)
            edifn_tmp = dgasi(ispg)
        endif 
        
        flx_gas(ispg,itflx,iz) = ( &
            & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt &
            & )         
        flx_gas(ispg,idif,iz) = ( &
            & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & )
        flx_gas(ispg,iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) &
            & )
        flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
        flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz)
        flx_gas(ispg,irxn_sld(:),iz) = &
            & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
            & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) 
        flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
        
        if (any(isnan(flx_gas(ispg,:,iz)))) then
            print *,flx_gas(ispg,:,iz)
        endif 
    enddo 

end do 
    

endsubroutine alsilicate_aq_gas_1D_v3_1

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

#ifdef no_intr_fincloc
!ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
function findloc(chrlist_in,chrspecific,dim)
implicit none
character(*),intent(in)::chrlist_in(:),chrspecific
integer,intent(in)::dim
integer findloc,i

findloc = 0
do i=1, size(chrlist_in,dim=dim)
    if (trim(adjustl(chrspecific)) == trim(adjustl(chrlist_in(i)))) then
        findloc = i
        return
    endif 

enddo 

endfunction findloc
!ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
#endif 

endprogram weathering