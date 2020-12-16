program weathering

implicit none 

integer nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nz
character(5),dimension(:),allocatable::chraq,chrsld,chrgas,chrrxn_ext 
character(500) sim_name,runname_save
real(kind=8) ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain,zml_ref,tc
integer count_dtunchanged_Max

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

sim_name = 'chkchk'

call get_bsdvalues( &
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
    & ,count_dtunchanged_Max,tc &
    & )
    
call weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
    & ,count_dtunchanged_Max,tc &
    & )

contains 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
subroutine weathering_main( &
    & nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain  &! input
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save &! input
    & ,count_dtunchanged_Max,tc &! input 
    & )

implicit none

!-----------------------------

real(kind=8),intent(in) :: ztot != 3.0d0 ! m
real(kind=8),intent(in) :: ttot  ! yr
! real(kind=8) dz
integer,intent(in) :: nz != 30 
real(kind=8) z(nz),dz(nz)
real(kind=8),intent(in) :: tc != 15.0d0 ! deg celsius
real(kind=8) dt  ! yr 
integer, parameter :: nt = 50000000
real(kind=8) time

real(kind=8),parameter :: rg = 8.3d-3   ! kJ mol^-1 K^-1
real(kind=8),parameter :: rg2 = 8.2d-2  ! L mol^-1 atm K^-1

real(kind=8),parameter :: tempk_0 = 273d0
real(kind=8),parameter :: sec2yr = 60d0*60d0*24d0*365d0

real(kind=8) pco2i

real(kind=8),parameter :: fr_an_la = 0.6d0 ! Anorthite fraction for Labradorite (Beerling et al., 2020); 0.5 - 0.7
real(kind=8),parameter :: fr_an_ab = 0.0d0 ! Anorthite fraction for albite (Beerling et al., 2020); 0.0 - 0.1
real(kind=8),parameter :: fr_an_an = 1.0d0 ! Anorthite fraction for anorthite (Beerling et al., 2020); 0.9 - 1.0

real(kind=8),parameter :: mvka = 99.52d0 ! cm3/mol; molar volume of kaolinite; Robie et al. 1978
real(kind=8),parameter :: mvfo = 43.79d0 ! cm3/mol; molar volume of Fo; Robie et al. 1978
real(kind=8),parameter :: mvab_0 = 100.07d0 ! cm3/mol; molar volume of Ab(NaAlSi3O8); Robie et al. 1978 
real(kind=8),parameter :: mvan_0 = 100.79d0 ! cm3/mol; molar volume of An (CaAl2Si2O8); Robie et al. 1978
real(kind=8),parameter :: mvab = fr_an_ab*mvan_0 + (1d0-fr_an_ab)*mvab_0 ! cm3/mol; molar volume of albite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
real(kind=8),parameter :: mvan = fr_an_an*mvan_0 + (1d0-fr_an_an)*mvab_0 ! cm3/mol; molar volume of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
real(kind=8),parameter :: mvcc = 36.934d0 ! cm3/mol; molar volume of Cc (CaCO3); Robie et al. 1978
real(kind=8),parameter :: mvpy = 23.94d0 ! cm3/mol; molar volume of Pyrite (FeS2); Robie et al. 1978
real(kind=8),parameter :: mvgb = 31.956d0 ! cm3/mol; molar volume of Gibsite (Al(OH)3); Robie et al. 1978
real(kind=8),parameter :: mvct = 108.5d0 ! cm3/mol; molar volume of Chrysotile (Mg3Si2O5(OH)4); Robie et al. 1978
real(kind=8),parameter :: mvfa = 46.39d0 ! cm3/mol; molar volume of Fayalite (Fe2SiO4); Robie et al. 1978
real(kind=8),parameter :: mvgt = 20.82d0 ! cm3/mol; molar volume of Goethite (FeO(OH)); Robie et al. 1978
real(kind=8),parameter :: mvcabd = 129.77d0 ! cm3/mol; molar volume of Ca-beidellite (Ca(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
real(kind=8),parameter :: mvdp = 66.09d0 ! cm3/mol; molar volume of Diopside (MgCaSi2O6);  Robie et al. 1978
real(kind=8),parameter :: mvhb = 248.09d0/3.55d0 ! cm3/mol; molar volume of Hedenbergite (FeCaSi2O6); from a webpage
real(kind=8),parameter :: mvkfs = 108.72d0 ! cm3/mol; molar volume of K-feldspar (KAlSi3O8); Robie et al. 1978
real(kind=8),parameter :: mvom = 5.3d0 ! cm3/mol; molar volume of OM (CH2O); Lasaga and Ohmoto 2002
real(kind=8),parameter :: mvomb = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8),parameter :: mvg1 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8),parameter :: mvg2 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8),parameter :: mvg3 = 5.3d0 ! cm3/mol; assumed to be same as mvom
real(kind=8),parameter :: mvamsi = 25.739d0 ! cm3/mol; molar volume of amorphous silica taken as cristobalite (SiO2); Robie et al. 1978
real(kind=8),parameter :: mvarg = 34.15d0 ! cm3/mol; molar volume of aragonite; Robie et al. 1978
real(kind=8),parameter :: mvdlm = 64.34d0 ! cm3/mol; molar volume of dolomite; Robie et al. 1978
real(kind=8),parameter :: mvhm = 30.274d0 ! cm3/mol; molar volume of hematite; Robie et al. 1978
real(kind=8),parameter :: mvill = 139.35d0 ! cm3/mol; molar volume of illite (K0.6Mg0.25Al2.3Si3.5O10(OH)2); Wolery and Jove-Colon 2004
real(kind=8),parameter :: mvanl = 97.49d0 ! cm3/mol; molar volume of analcime (NaAlSi2O6*H2O); Robie et al. 1978
real(kind=8),parameter :: mvnph = 54.16d0 ! cm3/mol; molar volume of nepheline (NaAlSiO4); Robie et al. 1978
real(kind=8),parameter :: mvqtz = 22.688d0 ! cm3/mol; molar volume of quartz (SiO2); Robie et al. 1978
real(kind=8),parameter :: mvgps = 74.69d0 ! cm3/mol; molar volume of gypsum (CaSO4*2H2O); Robie et al. 1978
real(kind=8),parameter :: mvtm = 272.92d0 ! cm3/mol; molar volume of tremolite (Ca2Mg5(Si8O22)(OH)2); Robie et al. 1978
real(kind=8),parameter :: mvla = fr_an_la*mvan_0 + (1d0-fr_an_la)*mvab_0 ! cm3/mol; molar volume of Labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing

real(kind=8),parameter :: mwtka = 258.162d0 ! g/mol; formula weight of Ka; Robie et al. 1978
real(kind=8),parameter :: mwtfo = 140.694d0 ! g/mol; formula weight of Fo; Robie et al. 1978
real(kind=8),parameter :: mwtab_0 = 262.225d0 ! g/mol; formula weight of Ab; Robie et al. 1978
real(kind=8),parameter :: mwtan_0 = 278.311d0 ! g/mol; formula weight of An; Robie et al. 1978
real(kind=8),parameter :: mwtab = fr_an_ab*mwtan_0 + (1d0-fr_an_ab)*mwtab_0 ! g/mol; formula weight of albte (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
real(kind=8),parameter :: mwtan = fr_an_an*mwtan_0 + (1d0-fr_an_an)*mwtab_0 ! g/mol; formula weight of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
real(kind=8),parameter :: mwtcc = 100.089d0 ! g/mol; formula weight of Cc; Robie et al. 1978
real(kind=8),parameter :: mwtpy = 119.967d0 ! g/mol; formula weight of Py; Robie et al. 1978
real(kind=8),parameter :: mwtgb = 78.004d0 ! g/mol; formula weight of Gb; Robie et al. 1978
real(kind=8),parameter :: mwtct = 277.113d0 ! g/mol; formula weight of Ct; Robie et al. 1978
real(kind=8),parameter :: mwtfa = 203.778d0 ! g/mol; formula weight of Fa; Robie et al. 1978
real(kind=8),parameter :: mwtgt = 88.854d0 ! g/mol; formula weight of Gt; Robie et al. 1978
real(kind=8),parameter :: mwtcabd = 366.6252667d0 ! g/mol; formula weight of Cabd calculated from atmoic weight
real(kind=8),parameter :: mwtdp = 216.553d0 ! g/mol;  Robie et al. 1978
real(kind=8),parameter :: mwthb = 248.09d0 ! g/mol; from a webpage
real(kind=8),parameter :: mwtkfs = 278.33d0 ! g/mol; formula weight of Kfs; Robie et al. 1978
real(kind=8),parameter :: mwtom = 30d0 ! g/mol; formula weight of CH2O
real(kind=8),parameter :: mwtomb = 30d0 ! g/mol; formula weight of CH2O
real(kind=8),parameter :: mwtg1 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8),parameter :: mwtg2 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8),parameter :: mwtg3 = 30d0 ! g/mol; formula weight of CH2O
real(kind=8),parameter :: mwtamsi = 60.085d0 ! g/mol; formula weight of amorphous silica
real(kind=8),parameter :: mwtarg = 100.089d0 ! g/mol; formula weight of aragonite
real(kind=8),parameter :: mwtdlm = 184.403d0 ! g/mol; formula weight of dolomite
real(kind=8),parameter :: mwthm = 159.692d0 ! g/mol; formula weight of hematite
real(kind=8),parameter :: mwtill = 383.90053d0 ! g/mol; formula weight of Ill calculated from atmoic weight
real(kind=8),parameter :: mwtanl = 220.155d0 ! g/mol; formula weight of analcime
real(kind=8),parameter :: mwtnph = 142.055d0 ! g/mol; formula weight of nepheline
real(kind=8),parameter :: mwtqtz = 60.085d0 ! g/mol; formula weight of quartz
real(kind=8),parameter :: mwtgps = 172.168d0 ! g/mol; formula weight of gypsum
real(kind=8),parameter :: mwttm = 812.374d0 ! g/mol; formula weight of tremolite
real(kind=8),parameter :: mwtla = fr_an_la*mwtan + (1d0-fr_an_la)*mwtab ! g/mol; formula weight of Labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing

real(kind=8),parameter :: rho_grain = 2.7d0 ! g/cm3 as soil grain density 

! real(kind=8)::plant_rain = 1.4d-3 ! g C/g soil/yr; converted from 1.6d-4 mg C / g soil /hr from Georgiou et al. 2017 ! 
real(kind=8),intent(in)::plant_rain != 1d2 ! 1 t/ha/yr; approximate values from Vanveen et al. 1991 ! 
! real(kind=8)::plant_rain = 0.1d2 ! 

real(kind=8)::zsupp_plant = 0.3d0 !  e-folding decrease

! real(kind=8)::rainpowder = 40d2 !  g/m2/yr corresponding to 40 t/ha/yr (40x1e3x1e3/1e4)
! real(kind=8)::rainpowder = 0.5d2 !  g/m2/yr corresponding to 0.5 t/ha/yr (0.5x1e3x1e3/1e4)
real(kind=8),intent(in)::rainpowder != 30d2 !  g/m2/yr 
! real(kind=8)::rainpowder = 10d2 !  g/m2/yr corresponding to 10 t/ha/yr (0.5x1e3x1e3/1e4)


real(kind=8),intent(in)::zsupp != 0.3d0 !  e-folding decrease

real(kind=8) sat(nz), poro(nz), torg(nz), tora(nz)

! real(kind=8) :: poroi = 0.1d0 !*** default
real(kind=8),intent(in) :: poroi != 0.5d0

real(kind=8),intent(in) :: satup != 0.10d0

! real(kind=8) :: zsat = 30d0  ! water table depth [m] ** default 
real(kind=8),intent(in) :: zsat != 5d0  ! water table depth [m] 
! real(kind=8) :: zsat = 15d0

real(kind=8),intent(in) :: w != 5.0d-5 ! m yr^-1, uplift rate ** default 
! real(kind=8), parameter :: w = 1.0d-4 ! m yr^-1, uplift rate


! real(kind=8) :: qin = 1d-1 ! m yr^-1, advection (m3 water / m2 profile / yr)  ** default
real(kind=8),intent(in) :: qin != 10d-1 ! m yr^-1 
! real(kind=8) :: qin = 0.1d-1 ! m yr^-1 
real(kind=8) v(nz)

! real(kind=8) :: hr = 1d5 ! m^2 m^-3, reciprocal of hydraulic radius  ** default 
! real(kind=8) :: hr = 1d4 ! m^2 m^-3, reciprocal of hydraulic radius
real(kind=8) :: hrii = 1d5

! real(kind=8) :: p80 = 10d-6 ! m (**default?)
real(kind=8),intent(in) :: p80 != 1d-6 ! m 

real(kind=8) ssa_cmn,mvab_save,mvan_save,mvcc_save,mvfo_save,mvka_save,mvgb_save
real(kind=8),dimension(nz):: pro,prox,poroprev,hr,rough,hri

real(kind=8) kho,ucv,kco2,k1,kw,k2,khco2i

integer iz,it,ispa,ispg,isps,irxn,ispa2,ispg2,isps2,ico2

real(kind=8) error 
real(kind=8) :: tol = 1d-6

! integer, parameter :: nrec = 22
integer, parameter :: nrec = 20
real(kind=8) rectime(nrec)
character(3) chr
character(256) runname,workdir, chrz(3), chrq(3),base,fname, chrrain
character(500),intent(in):: runname_save
integer irec, iter


integer  iflx
! real(kind=8) :: maxdt = 10d0
real(kind=8) :: maxdt = 0.2d0 ! for basalt exp?

logical :: pre_calc = .false.
! logical :: pre_calc = .true.

logical :: read_data = .false.
! logical :: read_data = .true.

! logical :: incld_rough = .false.
logical :: incld_rough = .true.

! logical :: cplprec = .false.
logical :: cplprec = .true.

logical :: rain_wave = .false.
! logical :: rain_wave = .true.

logical :: al_inhibit = .false.
! logical :: al_inhibit = .true.

logical :: timestep_fixed = .false.
! logical :: timestep_fixed = .true.

! logical :: no_biot = .false.
logical :: no_biot = .true.

logical :: biot_fick = .false.
! logical :: biot_fick = .true.

logical :: biot_turbo2 = .false.
! logical :: biot_turbo2 = .true.

logical :: biot_labs = .false.
! logical :: biot_labs = .true.

logical :: biot_till = .false.
! logical :: biot_till = .true.

! logical :: display = .false.
logical :: display = .true.

! logical :: regular_grid = .false.
logical :: regular_grid = .true.

! logical :: method_precalc = .false.
logical :: method_precalc = .true.

! logical :: sld_enforce = .false.
logical :: sld_enforce = .true.

logical :: poroevol = .false.

logical :: surfevol1 = .false.
logical :: surfevol2 = .false.

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

real(kind=8) beta 

logical :: flgback = .false.
logical :: flgreducedt = .false.
logical :: flgreducedt_prev = .false.

real(kind=8) time_start, time_fin, progress_rate, progress_rate_prev
integer count_dtunchanged  
integer,intent(in):: count_dtunchanged_Max  

integer,intent(in)::nsp_sld != 5
integer,parameter::nsp_sld_2 = 14
integer,parameter::nsp_sld_all = 30
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
integer,parameter::nrxn_ext_all = 5
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
real(kind=8),dimension(nsp_sld,nz)::msldx,msld,ksld,omega,msldsupp,nonprec,rxnsld
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
real(kind=8),dimension(nsp_aq_all,2)::keqaq_s
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

real(kind=8),dimension(4,5 + nrxn_ext + nsp_sld,nz)::flx_co2sp
character(5),dimension(6)::chrco2sp

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

integer,dimension(nsp_aq)::iaqflx
integer,dimension(nsp_gas)::igasflx
integer,dimension(nsp_sld)::isldflx
integer,dimension(6)::ico2flx

integer,parameter::ibasaltrain = 15
integer isldprof,isldprof2,iaqprof,igasprof,isldsat,ibsd,irate 

logical,dimension(nsp_sld)::turbo2,labs,nonlocal,nobio,fick,till
real(kind=8),dimension(nz,nz,nsp_sld)::trans
real(kind=8),intent(in)::zml_ref
real(kind=8) dbl_ref
integer izml
integer :: nz_disp = 10

real(kind=8),dimension(nz)::so4f

real(kind=8) dt_prev

logical print_cb,ph_error,save_trans
character(500) print_loc
character(500),intent(in):: sim_name

real(kind=8) def_dust,def_rain,def_pr,def_OM_frc
character(5),dimension(5 + nrxn_ext + nsp_sld)::chrflx
character(20) chrfmt

integer::itflx,iadv,idif,irain,ires
data itflx,iadv,idif,irain/1,2,3,4/

integer,dimension(nsp_sld)::irxn_sld 
integer,dimension(nrxn_ext)::irxn_ext

integer nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save 
character(5),dimension(:),allocatable::chraq_save,chrsld_save,chrgas_save,chrrxn_ext_save 
real(kind=8),dimension(:,:),allocatable::msld_save,mgas_save,maq_save

integer t1,t2,t_rate,t_max,diff
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
irate = ibasaltrain + nsp_sld + nsp_gas + nsp_aq + 7

nflx = 5 + nrxn_ext + nsp_sld

do isps=1,nsp_sld
    irxn_sld(isps) = 4+isps
enddo 

do irxn=1,nrxn_ext
    irxn_ext(irxn) = 4+nsp_sld+irxn
enddo 

ires = nflx

chrflx(1:4) = (/'tflx ','adv  ','dif  ','rain '/)
if (nsp_sld > 0) chrflx(irxn_sld(:)) = chrsld
if (nrxn_ext > 0) chrflx(irxn_ext(:)) = chrrxn_ext
chrflx(nflx) = 'res  '

! print *,chrflx,irxn_sld,chrsld

! pause

! define all species and rxns definable in the model 
! note that rxns here exclude diss(/prec) of mineral 
! which are automatically included when associated mineral is chosen

chrsld_all = (/'fo   ','ab   ','an   ','cc   ','ka   ','gb   ','py   ','ct   ','fa   ','gt   ','cabd ' &
    & ,'dp   ','hb   ','kfs  ','om   ','omb  ','amsi ','arg  ','dlm  ','hm   ','ill  ','anl  ','nph  ' &
    & ,'qtz  ','gps  ','tm   ','la   ' &
    & ,'g1   ','g2   ','g3   '/)
chraq_all = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    '/)
chrgas_all = (/'pco2','po2 '/)
chrrxn_ext_all = (/'resp ','fe2o2','omomb','ombto','pyfe3'/)

! define the species and rxns explicitly simulated in the model in a fully coupled way
! should be chosen from definable species & rxn lists above 

! chrsld = (/'fo   ','ab   ','an   ','cc   ','ka   '/)
! chraq = (/'mg   ','si   ','na   ','ca   ','al   '/)
! chrgas = (/'pco2 ','po2  '/)
! chrrxn_ext = (/'resp '/)


! define solid species which can precipitate
! in default, all minerals only dissolve 
! should be chosen from the chrsld list
chrsld_2 = (/'cc   ','ka   ','gb   ','ct   ','gt   ','cabd ','amsi ','hm   ','ill  ','anl  ','gps  '  &
    ,'arg  ','dlm  ','qtz  '/) 

! below are species which are sensitive to pH 
chraq_ph = (/'mg   ','si   ','na   ','ca   ','al   ','fe2  ','fe3  ','so4  ','k    '/)
chrgas_ph = (/'pco2 '/)

chrco2sp = (/'co2g ','co2aq','hco3 ','co3  ','DIC  ','ALK  '/)

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

mv_all = (/mvfo,mvab,mvan,mvcc,mvka,mvgb,mvpy,mvct,mvfa,mvgt,mvcabd,mvdp,mvhb,mvkfs,mvom,mvomb,mvamsi &
    & ,mvarg,mvdlm,mvhm,mvill,mvanl,mvnph,mvqtz,mvgps,mvtm,mvla &
    & ,mvg1,mvg2,mvg3/)
mwt_all = (/mwtfo,mwtab,mwtan,mwtcc,mwtka,mwtgb,mwtpy,mwtct,mwtfa,mwtgt,mwtcabd,mwtdp,mwthb,mwtkfs,mwtom,mwtomb,mwtamsi &
    & ,mwtarg,mwtdlm,mwthm,mwtill,mwtanl,mwtnph,mwtqtz,mwtgps,mwttm,mwtla &
    & ,mwtg1,mwtg2,mwtg3/)

do isps = 1, nsp_sld 
    mv(isps) = mv_all(findloc(chrsld_all,chrsld(isps),dim=1))
    mwt(isps) = mwt_all(findloc(chrsld_all,chrsld(isps),dim=1))
enddo 

! maqi_all = 0d0
    
def_rain = 1d-20
def_pr = 1d-20
    
call get_rainwater( &
    & nsp_aq_all,chraq_all,def_rain &! input
    & ,maqi_all &! output
    & )
    
call get_parentrock( &
    & nsp_sld_all,chrsld_all,def_pr &! input
    & ,msldi_all &! output
    & )

msldi_all = msldi_all/mwt_all*rho_grain*1d6

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
mgasth_all = 1d-200
maqth_all = 1d-200
msldth_all = 1d-200


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
! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
! staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! Analcime; NaAlSi2O6*H2O
staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'anl',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! Nepheline; NaAlSiO4
staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'nph',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! K-feldspar; KAlSi3O8
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'k',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0
staq_all(findloc(chrsld_all,'kfs',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0
! Anothite; CaAl2Si2O8
! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
! staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'al',dim=1)) = 2d0
! Albite; CaxNa(1-x)Al(1+x)Si(3-x)O8
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_ab
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_ab
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_ab
staq_all(findloc(chrsld_all,'ab',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_ab
! Anothite; CaxNa(1-x)Al(1+x)Si(3-x)O8
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_an
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_an
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_an
staq_all(findloc(chrsld_all,'an',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_an
! Labradorite; CaxNa(1-x)Al(1+x)Si(3-x)O8
staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'ca',dim=1)) = fr_an_la
staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'na',dim=1)) = 1d0 - fr_an_la
staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'al',dim=1)) = 1d0 + fr_an_la
staq_all(findloc(chrsld_all,'la',dim=1), findloc(chraq_all,'si',dim=1)) = 3d0 - fr_an_la
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
! Hematite; Fe2O3
staq_all(findloc(chrsld_all,'hm',dim=1), findloc(chraq_all,'fe3',dim=1)) = 2d0
! Ca-beidellite; Ca(1/6)Al(7/3)Si(11/3)O10(OH)2
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0/6d0
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'al',dim=1)) = 7d0/3d0
staq_all(findloc(chrsld_all,'cabd',dim=1), findloc(chraq_all,'si',dim=1)) = 11d0/3d0
! Illite; K0.6Mg0.25Al2.3Si3.5O10(OH)2
staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'k',dim=1)) = 0.6d0
staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'mg',dim=1)) = 0.25d0
staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'al',dim=1)) = 2.3d0
staq_all(findloc(chrsld_all,'ill',dim=1), findloc(chraq_all,'si',dim=1)) = 3.5d0
! Diopside (MgCaSi2O6)
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'dp',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
! Hedenbergite (FeCaSi2O6)
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'hb',dim=1), findloc(chraq_all,'si',dim=1)) = 2d0
! Tremolite (Ca2Mg5(Si8O22)(OH)2)
staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'ca',dim=1)) = 2d0
staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'mg',dim=1)) = 5d0
staq_all(findloc(chrsld_all,'tm',dim=1), findloc(chraq_all,'si',dim=1)) = 8d0
! Amorphous silica; SiO2
staq_all(findloc(chrsld_all,'amsi',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
! Quartz; SiO2
staq_all(findloc(chrsld_all,'qtz',dim=1), findloc(chraq_all,'si',dim=1)) = 1d0
! Aragonite (CaCO3)
staq_all(findloc(chrsld_all,'arg',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'arg',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
! Dolomite (CaMg(CO3)2)
staq_all(findloc(chrsld_all,'dlm',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'dlm',dim=1), findloc(chraq_all,'mg',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'dlm',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 2d0
! Gypsum; CaSO4*2H2O
staq_all(findloc(chrsld_all,'gps',dim=1), findloc(chraq_all,'ca',dim=1)) = 1d0
staq_all(findloc(chrsld_all,'gps',dim=1), findloc(chraq_all,'so4',dim=1)) = 1d0
! OMs; CH2O
stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'g1',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'g2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0

stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 1d0
stgas_all(findloc(chrsld_all,'g3',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0
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
stgas_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1), findloc(chrgas_all,'po2',dim=1)) = -1d0/4d0
! SOC assimilation by microbes 
stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'om',dim=1)) = -1d0
stsld_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrsld_all,'omb',dim=1)) = 0.31d0
stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'pco2',dim=1)) = 0.69d0
stgas_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1), findloc(chrgas_all,'po2',dim=1)) = -0.69d0
! turnover of microbes 
stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'om',dim=1)) = 1d0
stsld_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1), findloc(chrsld_all,'omb',dim=1)) = -1d0
! pyrite oxidation by fe3
stsld_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chrsld_all,'py',dim=1)) = -1d0
staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe3',dim=1)) = -14d0
staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe2',dim=1)) = 15d0
staq_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'so4',dim=1)) = 2d0

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
! pyrite oxidation by fe3
stsld_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chrsld_all,'py',dim=1)) = 1d0
staq_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe2',dim=1)) = 1d0
staq_dext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1), findloc(chraq_all,'fe3',dim=1)) = 1d0

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

! rfrc_sld_plant_all = rfrc_sld_plant_all/mwt_all
! rfrc_sld_plant_all = rfrc_sld_plant_all/sum(rfrc_sld_plant_all)



do isps = 1, nsp_sld 
    rfrc_sld(isps) = rfrc_sld_all(findloc(chrsld_all,chrsld(isps),dim=1))
    rfrc_sld_plant(isps) = rfrc_sld_plant_all(findloc(chrsld_all,chrsld(isps),dim=1))
enddo


call get_switches( &
    & no_biot,biot_fick,biot_turbo2,biot_labs,biot_till,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc,regular_grid,sld_enforce &! inout
    & ,poroevol,surfevol1,surfevol2 &!
    & )

if (sld_enforce) nsp3 = nsp_aq + nsp_gas ! excluding solid phases

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
 
if (poroevol) then      
    base = trim(adjustl(base))//'_pevol'
endif 
if (surfevol1)then 
    base = trim(adjustl(base))//'_sevol1'
elseif (surfevol2) then 
    base = trim(adjustl(base))//'_sevol2'
#if defined(surfssa)
    base = trim(adjustl(base))//'_ssa'
#endif 
endif 

if (.not. regular_grid) then 
    base = trim(adjustl(base))//'_irr'
endif 

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

do ico2 = 1, 6
    ico2flx(ico2) = ibasaltrain + nsp_sld + nsp_aq + nsp_gas + ico2
enddo 

! print*,workdir
! print*,runname
! pause

call system ('mkdir -p '//trim(adjustl(workdir))//trim(adjustl(runname)))

call system ('cp gases.in solutes.in slds.in extrxns.in '//trim(adjustl(workdir))//trim(adjustl(runname)))

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

do ico2 = 1,6
    open(ico2flx(ico2), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
        & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', status='replace')
    write(ico2flx(ico2),trim(adjustl(chrfmt))) 'time',(chrflx(iflx),iflx=1,nflx)
    close(ico2flx(ico2))
enddo 

open(ibasaltrain, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'rain.txt', &
    & status='replace')
close(ibasaltrain)



!!!  MAKING GRID !!!!!!!!!!!!!!!!! 
beta = 1.00000000005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
beta = 1.00005d0  ! a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
call makegrid(beta,nz,ztot,dz,z,regular_grid)


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

rough = 1d0
if (incld_rough) then 
    rough = 10d0**(3.3d0)*p80**0.33d0 ! from Navarre-Sitchler and Brantley (2007)
endif 

ssa_cmn = -4.4528d0*log10(p80*1d6) + 11.578d0 ! m2/g

hrii = 1d0/p80
hri = hrii

hr = hri*rough
v = qin/poroi/sat
poro = poroi
torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)

! #ifdef surfssa
! hri = ssa_cmn*1d6/poro
! mvab_save = mvab
! mvan_save = mvan
! mvcc_save = mvcc
! mvfo_save = mvfo
! mvka_save = mvka
! mvab = mwtab 
! mvan = mwtan 
! mvcc = mwtcc 
! mvfo = mwtfo 
! mvka = mwtka 
! #endif 

dt = maxdt

dt = 1d-20 ! for basalt exp?

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
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro,poro,hr &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 

print_cb = .false. 
print_loc = './ph.txt'

call calc_pH_v7_2( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
    & ,print_cb,print_loc,z &! input 
    & ,pro,ph_error,so4f &! output
    & ) 

poroprev = poro

!  --------- read -----
if (read_data) then 
    ! runname_save = 'test_cpl_rain-0.40E+04_pevol_sevol1_q-0.10E-01_zsat-5' ! specifiy the file where restart data is stored 
    ! runname_save = runname  ! the working folder has the restart data 
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'prof_sld-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_sld-restart.txt')
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'prof_aq-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_aq-restart.txt')
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'prof_gas-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_gas-restart.txt')
    call system('cp '//trim(adjustl(workdir))//trim(adjustl(runname_save))//'/'//'bsd-save.txt '  &
        & //trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'bsd-restart.txt')
        
    call get_saved_variables_num( &
        & workdir,runname_save &! input
        & ,nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save &! output
        & )
    
    allocate(chraq_save(nsp_aq_save),chrsld_save(nsp_sld_save),chrgas_save(nsp_gas_save),chrrxn_ext_save(nrxn_ext_save))
    allocate(maq_save(nsp_aq_save,nz),msld_save(nsp_sld_save,nz),mgas_save(nsp_gas_save,nz))
        
    
    call get_saved_variables( &
        & workdir,runname_save &! input
        & ,nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save &! input
        & ,chraq_save,chrgas_save,chrsld_save,chrrxn_ext_save &! output
        & )
    
        
    open (isldprof, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_sld-restart.txt',  &
        & status ='old',action='read')
    open (iaqprof, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_aq-restart.txt',  &
        & status ='old',action='read')
    open (igasprof, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'prof_gas-restart.txt',  &
        & status ='old',action='read')
    open (ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'//'bsd-restart.txt',  &
        & status ='old',action='read')
    
    read (isldprof,'()')
    read (iaqprof,'()')
    read (igasprof,'()')
    read (ibsd,'()')
    
    do iz = 1, Nz
        read (isldprof,*) z(iz),(msld_save(isps,iz),isps=1,nsp_sld_save),time
        read (iaqprof,*) z(iz),(maq_save(ispa,iz),ispa=1,nsp_aq_save),pro(iz),time
        read (igasprof,*) z(iz),(mgas_save(ispg,iz),ispg=1,nsp_gas_save),time
        read (ibsd,*) z(iz),poro(iz),sat(iz),v(iz),hr(iz),time
    enddo 
    close(isldprof)

    pro = 10d0**(-pro) ! read data is -log10 (pro)
    time = 0d0
    
    torg = poro**(3.4d0-2.0d0)*(1.0d0-sat)**(3.4d0-1.0d0)
    tora = poro**(3.4d0-2.0d0)*(sat)**(3.4d0-1.0d0)
    
    do isps = 1,nsp_sld_save
        if (any(chrsld == chrsld_save(isps))) then 
            msld(findloc(chrsld,chrsld_save(isps),dim=1),:) = msld_save(isps,:)
        elseif (any(chrsld_cnst == chrsld_save(isps))) then
            msldc(findloc(chrsld_cnst,chrsld_save(isps),dim=1),:) = msld_save(isps,:)
        else 
            print *,'error in re-assignment of sld conc.'
        endif 
    enddo 
    
    do ispa = 1,nsp_aq_save
        if (any(chraq == chraq_save(ispa))) then 
            maq(findloc(chraq,chraq_save(ispa),dim=1),:) = maq_save(ispa,:)
        elseif (any(chraq_cnst == chraq_save(ispa))) then
            maqc(findloc(chraq_cnst,chraq_save(ispa),dim=1),:) = maq_save(ispa,:)
        else 
            print *,'error in re-assignment of aq conc.'
        endif 
    enddo 
    
    do ispg = 1,nsp_gas_save
        if (any(chrgas == chrgas_save(ispg))) then 
            mgas(findloc(chrgas,chrgas_save(ispg),dim=1),:) = mgas_save(ispg,:)
        elseif (any(chrgas_cnst == chrgas_save(ispg))) then
            mgasc(findloc(chrgas_cnst,chrgas_save(ispg),dim=1),:) = mgas_save(ispg,:)
        else 
            print *,'error in re-assignment of gas conc.'
        endif 
    enddo 
        
    if (display) then
        write(chrfmt,'(i0)') nz_disp
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        
        print *
        print *,' [concs] '
        print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
        if (nsp_aq>0) then 
            print *,' < aq species >'
            do ispa = 1, nsp_aq
                print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_sld>0) then 
            print *,' < sld species >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_gas>0) then 
            print *,' < gas species >'
            do ispg = 1, nsp_gas
                print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
    endif      
endif
    
call coefs_v2( &
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro,poro,hr &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
    
! zml_ref = 1.5d0 ! mixed layer depth [m]
! zml_ref = 0.5d0 ! mixed layer depth [m]
dbl_ref = 0d0  
labs = .false.
turbo2 = .false.
nobio = .false.
till = .false.
    
if (no_biot) nobio = .true.
if (biot_turbo2) turbo2 = .true.
if (biot_fick) fick = .true.
if (biot_labs) labs = .true.
if (biot_till) till = .true.

save_trans = .true.
call make_transmx(  &
    & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  &! input
    & ,trans,nonlocal,izml  &! output 
    & )
    
! --------- loop -----
print *, 'about to start time loop'
it = 0
irec = 0

count_dtunchanged = 0

!! @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@

do while (it<nt)
    ! call cpu_time(time_start)
    call system_clock(t1)
    
    if (display) then 
        print *
        print *, '-----------------------------------------'
        print '(i11,a)', it,': time iteration'
        print '(E11.3,a)',time,': time [yr]' 
        print *
    endif
    dt_prev = dt
    
    if (time>rectime(nrec)) exit
    
    if (it == 0) then 
        maxdt = 0.2d0
    endif 
    
    if (timestep_fixed) then 
        maxdt = 0.2d0
        ! maxdt = 0.02d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.005d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.002d0 ! working with p80 = 10 um
        ! maxdt = 0.001d0 ! when calcite is included smaller time step must be assumed 
        ! maxdt = 0.0005d0 ! working with p80 = 1 um
        
        ! if (time<1d-2) then  
            ! maxdt = 1d-6 
        ! elseif (time>=1d-2 .and. time<1d-1) then 
        ! if (time<1d-3) then  
            ! maxdt = 1d-7 
        ! elseif (time>=1d-3 .and. time<1d-2) then  
            ! maxdt = 1d-6 
        ! elseif (time>=1d-2 .and. time<1d-1) then  
            ! maxdt = 1d-5 
        ! elseif (time>=1d-1 .and. time<1d0) then  
        if ( time<1d0) then  
            maxdt = 1d-4 
        elseif (time>=1d0 .and. time<1d1) then 
            maxdt = 1d-3 
        elseif (time>=1d1 .and. time<1d2) then 
            maxdt = 1d-2  
        ! elseif (time>=1d2 .and. time<1d3) then 
            ! maxdt = 1d-1 
        ! elseif (time>=1d3 .and. time<1d4) then 
            ! maxdt = 1d0 
        ! elseif (time>=1d4 .and. time<1d5) then 
            ! maxdt = 1d1 
        ! elseif (time>=1d5 ) then 
            ! maxdt = 1d2 
        endif 
        
        ! maxdt = maxdt * 1d-1
    endif 
    
    ! count_dtunchanged_Max = 1000
    ! count_dtunchanged_Max = 10
    ! if (sld_enforce) count_dtunchanged_Max = 10
    ! if (dt<1d-5) then 
        ! count_dtunchanged_Max = 10
    ! elseif (dt>=1d-5 .and. dt<1d0) then
        ! count_dtunchanged_Max = 100
    ! elseif (dt>=1d0 ) then 
        ! count_dtunchanged_Max = 1000
    ! endif 

    ! -------- modifying dt --------

    !        if ((iter <= 10).and.(dt<1d1)) then
    if (dt<maxdt) then
        ! dt = dt*1.01d0
        dt = dt*10d0
        if (dt>maxdt) dt = maxdt
    endif
    ! if (iter > 300) then
        ! dt = dt/10d0
    ! end if
    
    ! if (dt/=dt_prev) pre_calc = .true.

    ! incase temperature&ph change
        
    call coefs_v2( &
        & nz,rg,rg2,tc,sec2yr,tempk_0,pro,poro,hr &! input
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
        & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
        & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all &!input
        & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! output
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
    pco2i = mgasi_all(findloc(chrgas_all,'pco2',dim=1))
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
    save_trans = .false.
    call make_transmx(  &
        & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  &! input
        & ,trans,nonlocal,izml  &! output 
        & )


    error = 1d4
    ! iter=0

100 continue

    mgasx = mgas
    msldx = msld
    maqx = maq
    
    prox = pro  


    !  raining dust & OM 
    maqsupp = 0d0
    mgassupp = 0d0
    do isps = 1, nsp_sld
        if (no_biot) then 
            msldsupp(isps,:) = rainpowder*rfrc_sld(isps)*exp(-z/zsupp)/zsupp
        else 
            msldsupp(isps,1) = rainpowder*rfrc_sld(isps)/dz(1)
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
                & + plant_rain/12d0*rfrc_sld_plant(isps)/dz(1) ! when plant_rain is in g_C/m2/yr
        endif 
    enddo 
    
    if (sld_enforce) then 
        if (any(chrgas=='pco2')) then 
            ! mgassupp(findloc(chrgas,'pco2',dim=1),:) = plant_rain/12d0*exp(-z/zsupp_plant)/zsupp_plant
            mgassupp(findloc(chrgas,'pco2',dim=1),:) = plant_rain/12d0/ztot
        endif 
    endif 
    
    
    ! if (it==0) pre_calc = .true.
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

    ! if ((.not.read_data) .and. it == 0 .and. iter == 0) then 
        ! do ispa = 1, nsp_aq
            ! if (chraq(ispa)/='so4') then
                ! maqx(ispa,1:) = 1d2
            ! endif 
        ! enddo
    ! endif 

    call alsilicate_aq_gas_1D_v3_1( &
        ! new input 
        & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
        & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
        & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
        & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
        & ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext &
        & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst &
        & ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all &
        & ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s &
        & ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &
        & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all &
        & ,turbo2,labs,trans,method_precalc,display,chrflx,sld_enforce &! input
        !  old inputs
        & ,hr,poro,z,dz,w,sat,pro,poroprev,tora,v,tol,it,nflx,kw & 
        & ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0  &
        ! old inout
        & ,dt,flgback &    
        ! output 
        & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp & 
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
        ! pre_calc = .true.
        dt = dt/1d1
        go to 100
    endif    
    
    if (poroevol) then 
        poroprev = poro
! #ifdef surfssa
        ! mvab = mvab_save 
        ! mvan = mvan_save 
        ! mvcc = mvcc_save 
        ! mvfo = mvfo_save 
        ! mvka = mvka_save 
! #endif 
        ! poro = poroi + (mabi-mabx)*(mvab)*1d-6  &
            ! & +(mfoi-mfox)*(mvfo)*1d-6 &
            ! & +(mani-manx)*(mvan)*1d-6 &
            ! & +(mcci-mccx)*(mvcc)*1d-6 &
            ! & +(mkai-mkax)*(mvka)*1d-6 
        poro = poroi
        do isps=1,nsp_sld
            poro = poro + (msldi(isps)-msldx(isps,:))*mv(isps)*1d-6
        enddo 
! #ifdef surfssa
        ! mvab = mwtab 
        ! mvan = mwtan 
        ! mvcc = mwtcc 
        ! mvfo = mwtfo 
        ! mvka = mwtka 
! #endif 
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
        hr = hri*rough
        if (surfevol1 ) then 
            hr = hri*rough*((1d0-poro)/(1d0-poroi))**(2d0/3d0)
        endif 
        if (surfevol2 ) then 
            hr = hri*rough*(poro/poroi)**(2d0/3d0)  ! SA increases with porosity 
        endif 
        if (display) then 
            write(chrfmt,'(i0)') nz_disp
            chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
            print *
            print *,' [porosity & surface area]'
            print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
            print trim(adjustl(chrfmt)),'poro',(poro(iz),iz=1,nz, nz/nz_disp)
            print trim(adjustl(chrfmt)),'SA',(hr(iz),iz=1,nz, nz/nz_disp)
            print *
        endif 
    endif  


    if (display) then 
        write(chrfmt,'(i0)') nz_disp
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        
        print *
        print *,' [concs] '
        print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
        if (nsp_aq>0) then 
            print *,' < aq species >'
            do ispa = 1, nsp_aq
                print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_sld>0) then 
            print *,' < sld species >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_gas>0) then 
            print *,' < gas species >'
            do ispg = 1, nsp_gas
                print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        
        print *
        print *,' [saturation & pH] '
        if (nsp_sld>0) then 
            print *,' < sld species omega >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        print *,' < pH >'
        print trim(adjustl(chrfmt)), 'ph', (-log10(prox(iz)),iz=1,nz, nz/nz_disp)
        
        
        
        write(chrfmt,'(i0)') nflx
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'
        
        print *
        print *,' [fluxes] '
        print trim(adjustl(chrfmt)),' ',(chrflx(iflx),iflx=1,nflx)
        
        write(chrfmt,'(i0)') nflx
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        if (nsp_aq>0) then 
            print *,' < aq species >'
            do ispa = 1, nsp_aq
                print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        if (nsp_sld>0) then 
            print *,' < sld species >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        if (nsp_gas>0) then 
            print *,' < gas species >'
            do ispg = 1, nsp_gas
                print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
            enddo 
        endif 
        
    endif 

    ! stop
    
    mgas = mgasx
    maq = maqx
    msld = msldx
    
    pro = prox

    if (time > savetime) then 
        
        open(isldprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_sld-save.txt', status='replace')
        open(igasprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_gas-save.txt', status='replace')
        open(iaqprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_aq-save.txt', status='replace')
        open(ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
            & //'bsd-save.txt', status='replace')
            
        write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
        write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
        write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' time '

        do iz = 1, Nz
            write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
            write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
            write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
            write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
        end do

        close(isldprof)
        close(iaqprof)
        close(igasprof)
        close(ibsd)
        
        savetime = savetime + dsavetime
    endif 

    if (time>=rectime(irec+1)) then
        write(chr,'(i3.3)') irec+1
        
        
        print_cb = .true. 
        print_loc = trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'chrge_balance-'//chr//'.txt'

        call calc_pH_v7_2( &
            & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
            & ,print_cb,print_loc,z &! input 
            & ,prox,ph_error,so4f &! output
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
        open(irate, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
            & //'rate-'//chr//'.txt', status='replace')
            
        write(chrfmt,'(i0)') nsp_sld+2
        chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
        write(isldprof,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
        write(isldprof2,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
        write(isldsat,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),'time'
        write(chrfmt,'(i0)') nsp_aq+3
        chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
        write(iaqprof,trim(adjustl(chrfmt))) 'z',(chraq(isps),isps=1,nsp_aq),'ph','time'
        write(chrfmt,'(i0)') nsp_gas+2
        chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
        write(igasprof,trim(adjustl(chrfmt))) 'z',(chrgas(isps),isps=1,nsp_gas),'time'
        write(chrfmt,'(i0)') 6
        chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a10))'
        write(ibsd,trim(adjustl(chrfmt))) 'z','poro', 'sat', 'v[m/yr]', 'm2/m3' ,'time'
        write(chrfmt,'(i0)') 2 + nsp_sld + nrxn_ext
        chrfmt = '('//trim(adjustl(chrfmt))//'(1x,a5))'
        write(irate,trim(adjustl(chrfmt))) 'z',(chrsld(isps),isps=1,nsp_sld),(chrrxn_ext(irxn),irxn=1,nrxn_ext),'time'

        do iz = 1, Nz
            write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
            write(isldprof2,*) z(iz),(msldx(isps,iz)*mwt(isps)*1d2/(rho_grain*1d6),isps = 1, nsp_sld),time
            write(isldsat,*) z(iz),(omega(isps,iz),isps = 1, nsp_sld),time
            write(igasprof,*) z(iz),(mgasx(ispg,iz),ispg = 1, nsp_gas),time
            write(iaqprof,*) z(iz),(maqx(ispa,iz),ispa = 1, nsp_aq),-log10(prox(iz)),time
            write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
            write(irate,*) z(iz), (rxnsld(isps,iz),isps=1,nsp_sld),(rxnext(irxn,iz),irxn=1,nrxn_ext), time 
        end do
        irec=irec+1

        close(isldprof)
        close(isldprof2)
        close(isldsat)
        close(iaqprof)
        close(igasprof)
        close(ibsd)
        close(irate)
        
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
        
        do ico2=1,6 
            open(ico2flx(ico2), file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
                & //'flx_co2sp-'//trim(adjustl(chrco2sp(ico2)))//'.txt', action='write',status='old',access='append')
            if (ico2 .le. 4) then 
                write(ico2flx(ico2),*) time,(sum(flx_co2sp(ico2,iflx,:)*dz(:)),iflx=1,nflx)
            elseif (ico2 .eq. 5) then 
                write(ico2flx(ico2),*) time &
                    & ,(sum(flx_co2sp(2,iflx,:)*dz(:))+sum(flx_co2sp(3,iflx,:)*dz(:))+sum(flx_co2sp(4,iflx,:)*dz(:)) &
                    & ,iflx=1,nflx)
            elseif (ico2 .eq. 6) then 
                write(ico2flx(ico2),*) time &
                    & ,(sum(flx_co2sp(3,iflx,:)*dz(:))+2d0*sum(flx_co2sp(4,iflx,:)*dz(:)) &
                    & ,iflx=1,nflx)
            endif 
            close(ico2flx(ico2))
        enddo 
        
        open(isldprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_sld-save.txt', status='replace')
        open(igasprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_gas-save.txt', status='replace')
        open(iaqprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
            & //'prof_aq-save.txt', status='replace')
        open(ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
            & //'bsd-save.txt', status='replace')
            
        write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
        write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
        write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
        write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' time '

        do iz = 1, Nz
            write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
            write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
            write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
            write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
        end do

        close(isldprof)
        close(iaqprof)
        close(igasprof)
        close(ibsd)
        
    end if

    it = it + 1
    time = time + dt
    count_dtunchanged = count_dtunchanged + 1
    
    progress_rate_prev = progress_rate
    
    ! call cpu_time(time_fin)
    call system_clock(t2,t_rate,t_max)
    if ( t2 < t1 ) then
        diff = (t_max - t1) + t2 + 1
    else
        diff = t2 - t1
    endif
    
    ! progress_rate = dt/(time_fin-time_start)*sec2yr ! (model yr)/(computer yr)
    ! progress_rate = (time_fin-time_start) ! (computer sec)
    progress_rate = diff/dble(t_rate) ! (computer sec)
    
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
    
    ! if (progress_rate ==0d0 .or. progress_rate_prev ==0d0) maxdt = 1d2
    ! print *,progress_rate,progress_rate_prev,maxdt,time_fin,time_start
    if (isnan(maxdt).or.maxdt ==0d0) then 
    ! if (.true.) then 
        print *
        print *, 'maxdt is nan or zero',progress_rate,progress_rate_prev,maxdt,time_fin,time_start
        stop
    endif 
    
    if (display) then 
        print *
        print '(E11.3,a)',progress_rate,': computation time per iteration [sec]'
        print '(E11.3,a)',maxdt, ': maxdt [yr]'
        print '(i11,a)',count_dtunchanged,': count_dtunchanged'
        print *, '-----------------------------------------'
        print *
    endif 
    
end do
        
        
open(isldprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
    & //'prof_sld-save.txt', status='replace')
open(igasprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
    & //'prof_gas-save.txt', status='replace')
open(iaqprof,file=trim(adjustl(workdir))//trim(adjustl(runname))//'/' &
    & //'prof_aq-save.txt', status='replace')
open(ibsd, file=trim(adjustl(workdir))//trim(adjustl(runname))//'/'  &
    & //'bsd-save.txt', status='replace')
            
write(isldprof,*) ' z ',(chrsld(isps),isps=1,nsp_sld),' time '
write(iaqprof,*) ' z ',(chraq(isps),isps=1,nsp_aq),' ph ',' time '
write(igasprof,*) ' z ',(chrgas(isps),isps=1,nsp_gas),' time '
write(ibsd,*) ' z ',' poro ', ' sat ', ' v[m/yr] ', ' m2/m3 ' ,' time '

do iz = 1, Nz
    write(isldprof,*) z(iz),(msldx(isps,iz),isps = 1, nsp_sld),time
    write(igasprof,*) z(iz),(mgasx(isps,iz),isps = 1, nsp_gas),time
    write(iaqprof,*) z(iz),(maqx(isps,iz),isps = 1, nsp_aq),-log10(prox(iz)),time
    write(ibsd,*) z(iz), poro(iz),sat(iz),v(iz),hr(iz),time
end do

close(isldprof)
close(iaqprof)
close(igasprof)
close(ibsd)

call system ('cp gases.in '//trim(adjustl(workdir))//trim(adjustl(runname))//'/gases.save')
call system ('cp solutes.in '//trim(adjustl(workdir))//trim(adjustl(runname))//'/solutes.save')
call system ('cp slds.in '//trim(adjustl(workdir))//trim(adjustl(runname))//'/slds.save')
call system ('cp extrxns.in '//trim(adjustl(workdir))//trim(adjustl(runname))//'/extrxns.save')

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

subroutine get_saved_variables_num( &
    & workdir,runname_save &! input 
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! output
    & )
implicit none

integer,intent(out):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext
character(256),intent(in):: workdir,runname_save
character(500) file_name

file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/slds.save'
call Console4(file_name,nsp_sld)
file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/solutes.save'
call Console4(file_name,nsp_aq)
file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/gases.save'
call Console4(file_name,nsp_gas)
file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/extrxns.save'
call Console4(file_name,nrxn_ext)

nsp_sld = nsp_sld - 1
nsp_aq = nsp_aq - 1
nsp_gas = nsp_gas - 1
nrxn_ext = nrxn_ext - 1

endsubroutine get_saved_variables_num

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_saved_variables( &
    & workdir,runname_save &! input 
    & ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext &! input
    & ,chraq,chrgas,chrsld,chrrxn_ext &! output
    & )
implicit none

integer,intent(in):: nsp_sld,nsp_aq,nsp_gas,nrxn_ext
character(5),dimension(nsp_sld),intent(out)::chrsld 
character(5),dimension(nsp_aq),intent(out)::chraq 
character(5),dimension(nsp_gas),intent(out)::chrgas 
character(5),dimension(nrxn_ext),intent(out)::chrrxn_ext 
character(256),intent(in):: workdir,runname_save

character(500) file_name
integer ispa,ispg,isps,irxn

if (nsp_aq>=1) then 
    file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/solutes.save'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do ispa =1,nsp_aq
        read(50,*) chraq(ispa) 
    enddo 
    close(50)
endif 

if (nsp_sld>=1) then 
    file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/slds.save'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do isps =1,nsp_sld
        read(50,*) chrsld(isps) 
    enddo  
    close(50)
endif 

if (nsp_gas>=1) then 
    file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/gases.save'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do ispg =1,nsp_gas
        read(50,*) chrgas(ispg) 
    enddo 
    close(50)
endif 

if (nrxn_ext>=1) then 
    file_name = trim(adjustl(workdir))//trim(adjustl(runname_save))//'/extrxns.save'
    open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
    read(50,'()')
    do irxn =1,nrxn_ext
        read(50,*) chrrxn_ext(irxn) 
    enddo 
    close(50)
endif 

endsubroutine get_saved_variables

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine get_bsdvalues( &
    & nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save &! output
    & ,count_dtunchanged_Max,tc &
    & )
implicit none

integer,intent(out):: nz
real(kind=8),intent(out)::ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,plant_rain,zml_ref,tc
character(500),intent(out)::sim_name,runname_save
integer,intent(out)::count_dtunchanged_Max

character(500) file_name

file_name = './frame.in'
open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')
read(50,*) ztot
read(50,*) nz
read(50,*) ttot
read(50,*) tc
read(50,*) rainpowder
read(50,*) plant_rain
read(50,*) zsupp
read(50,*) poroi
read(50,*) satup
read(50,*) zsat
read(50,*) zml_ref
read(50,*) w
read(50,*) qin
read(50,*) p80
read(50,*) count_dtunchanged_Max
read(50,*) runname_save
read(50,'()')
read(50,*) sim_name
close(50)

print*,'nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain,runname_save,count_dtunchanged_Max'
print*,nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,w,qin,p80,sim_name,plant_rain,runname_save,count_dtunchanged_Max

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
    & no_biot,biot_fick,biot_turbo2,biot_labs,biot_till,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc,regular_grid,sld_enforce &! inout
    & ,poroevol,surfevol1,surfevol2 &! inout
    & )
implicit none

character(100) chr_tmp
logical,intent(inout):: no_biot,biot_turbo2,biot_labs,display,read_data,incld_rough &
    & ,al_inhibit,timestep_fixed,method_precalc,biot_fick,biot_till,regular_grid,sld_enforce &
    & ,poroevol,surfevol1,surfevol2

character(500) file_name
integer i,n_tmp

file_name = './switches.in'

open(50,file=trim(adjustl(file_name)),status = 'old',action='read')
read(50,'()')

read(50,*) no_biot,chr_tmp
read(50,*) biot_fick,chr_tmp
read(50,*) biot_turbo2,chr_tmp
read(50,*) biot_labs,chr_tmp
read(50,*) biot_till,chr_tmp
read(50,*) display,chr_tmp
read(50,*) read_data,chr_tmp
read(50,*) incld_rough,chr_tmp
read(50,*) al_inhibit,chr_tmp
read(50,*) timestep_fixed,chr_tmp
read(50,*) method_precalc,chr_tmp
read(50,*) regular_grid,chr_tmp
read(50,*) sld_enforce,chr_tmp
read(50,*) poroevol,chr_tmp
read(50,*) surfevol1,chr_tmp
read(50,*) surfevol2,chr_tmp

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

subroutine makegrid(beta,nz,ztot,dz,z,regular_grid)  !  making grid, after Hoffmann & Chiang, 2000
implicit none
integer(kind=4),intent(in) :: nz
logical,intent(in)::regular_grid
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

if (regular_grid) then 
    dz = ztot/nz  ! when implementing regular grid
endif 

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
    & nz,rg,rg2,tc,sec2yr,tempk_0,pro,poro,hr &! input
    & ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all &! input
    & ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all &! input
    & ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all &!input
    & ,ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! output
    & ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,rg2,tc,sec2yr,tempk_0
real(kind=8),dimension(nz),intent(in)::pro,poro,hr
real(kind=8),dimension(nz)::oh,po2,kin,dkin_dmsp
real(kind=8) kho,po2th,mv_tmp,therm,ss_x
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
real(kind=8),dimension(nsp_aq_all,2),intent(out)::keqaq_s
real(kind=8),dimension(nsp_sld_all,nz),intent(out)::ksld_all
real(kind=8),dimension(nsp_sld_all),intent(in)::mv_all
real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all
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

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

integer isps

! real(kind=8)::thon = 1d0
real(kind=8)::thon = -1d100
character(5) mineral

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
keqaq_s = 0d0

! SO4-2 + H+ = HSO4- 
keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4) = &
    & k_arrhenius(10d0**(1.988d0),25d0+tempk_0,tc+tempk_0,3.85d0*cal2j,rg) ! from PHREEQC.DAT


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
! Al+3 + SO4-2 = AlSO4+
keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4) = &
    & k_arrhenius(10d0**(3.5d0),25d0+tempk_0,tc+tempk_0,2.29d0*cal2j,rg) ! from PHREEQC.DAT 
! Al+3 + 2SO4-2 = Al(SO4)2-
keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42) = &
    & k_arrhenius(10d0**(5.0d0),25d0+tempk_0,tc+tempk_0,3.11d0*cal2j,rg) ! from PHREEQC.DAT 


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
! Mg+2 + SO4-2 = MgSO4
keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4) = & 
    & k_arrhenius(10d0**(2.37d0),25d0+tempk_0,tc+tempk_0, 4.550d0*cal2j,rg) ! from PHREEQC.DAT 

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
! Ca+2 + SO4-2 = CaSO4
keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4) = &
    & k_arrhenius(10d0**(2.25d0),25d0+tempk_0,tc+tempk_0,1.325d0*cal2j,rg) ! from PHREEQC.DAT 

    
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
! Fe+2 + SO4-2 = FeSO4
keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4) = &
    & k_arrhenius(10d0**(2.25d0),25d0+tempk_0,tc+tempk_0,3.230d0*cal2j,rg) ! from PHREEQC.DAT 


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
! Fe+3 + SO4-2 = FeSO4+
keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4) = &
    & k_arrhenius(10d0**(4.04d0),25d0+tempk_0,tc+tempk_0,3.91d0*cal2j,rg) ! from PHREEQC.DAT 
! Fe+3 + 2 SO4-2 = Fe(SO4)2-
keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42) = &
    & k_arrhenius(10d0**(5.38d0),25d0+tempk_0,tc+tempk_0,4.60d0*cal2j,rg) ! from PHREEQC.DAT 



! Na+ + CO3-2 = NaCO3-
keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3) = & 
    & k_arrhenius(10d0**(1.27d0),25d0+tempk_0,tc+tempk_0, 8.91d0*cal2j,rg) ! from PHREEQC.DAT 
! Na+ + H + CO3- = NaHCO3
keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3) = & 
    & k_arrhenius(10d0**(-0.25d0),25d0+tempk_0,tc+tempk_0, -1d0*cal2j,rg) &! from PHREEQC.DAT for Na+ + HCO3- = NaHCO3
    & /keqgas_h(findloc(chrgas_all,'pco2',dim=1),ieqgas_h2)  ! HCO3- = CO32- + H+
! Na+ + SO4-2 = NaSO4-
keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4) = & 
    & k_arrhenius(10d0**(0.7d0),25d0+tempk_0,tc+tempk_0, 1.120d0*cal2j,rg) ! from PHREEQC.DAT 



! K+ + SO4-2 = KSO4-
keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4) = & 
    & k_arrhenius(10d0**(0.85d0),25d0+tempk_0,tc+tempk_0, 2.250d0*cal2j,rg) ! from PHREEQC.DAT 


! keqaq_s = 0d0


!!! ----------- Solid phases ------------------------!!
ksld_all = 0d0 
keqsld_all = 0d0

do isps = 1, nsp_sld_all
    mv_tmp = mv_all(isps)
    mineral = chrsld_all(isps)
    
    call sld_kin( &
        & nz,rg,tc,sec2yr,tempk_0,pro,poro,hr,kw,kho,mv_tmp &! input
        & ,mineral,'xxxxx' &! input 
        & ,kin,dkin_dmsp &! output
        & ) 
    ksld_all(isps,:) = kin
    
    ! check for solid solution 
    select case (trim(adjustl(mineral))) 
        case('la','ab','an')
            ss_x = staq_all(isps, findloc(chraq_all,'ca',dim=1))
        case default 
            ss_x = 0d0 ! non-zero if it is a solid solution 
    endselect 
    
    call sld_therm( &
        & rg,tc,tempk_0,ss_x &! input
        & ,mineral &! input
        & ,therm &! output
        & ) 
    keqsld_all(isps) = therm
enddo


!--------- other reactions -------------! 
krxn1_ext_all = 0d0
krxn2_ext_all = 0d0

krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) = &
    & 1d0
    ! & max(8.0d13*60.0d0*24.0d0*365.0d0*(kw/pro)**2.0d0, 1d-7*60.0d0*24.0d0*365.0d0) &   
    ! mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)excluding the term (c*po2)
    ! & *merge(0d0,1d0,po2<po2th*thon)
     
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




krxn1_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1),:) = & 
    & 10.0d0**(-6.07d0)*60.0d0*60.0d0*24.0d0*365.0d0  !! excluding the term (fe3**0.93/fe2**0.40)  


endsubroutine coefs_v2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine sld_kin( &
    & nz,rg,tc,sec2yr,tempk_0,prox,poro,hr,kw,kho,mv_tmp &! input
    & ,mineral,dev_sp &! input 
    & ,kin,dkin_dmsp &! output
    & ) 
implicit none

integer,intent(in)::nz
real(kind=8),intent(in)::rg,tc,sec2yr,tempk_0,mv_tmp,kw,kho
real(kind=8),dimension(nz),intent(in)::prox,poro,hr

real(kind=8) :: cal2j = 4.184d0 

character(5),intent(in)::mineral,dev_sp
real(kind=8),dimension(nz),intent(out)::kin,dkin_dmsp
real(kind=8) mh,moh,kinn_ref,kinh_ref,kinoh_ref,ean,eah,eaoh,tc_ref

! real(kind=8) k_arrhenius


kin = 0d0
dkin_dmsp = 0d0

select case(trim(adjustl(mineral)))
    case('ka')
        mh = 0.777d0
        moh = -0.472d0
        kinn_ref = 10d0**(-13.18d0)*sec2yr
        kinh_ref = 10d0**(-11.31d0)*sec2yr
        kinoh_ref = 10d0**(-17.05d0)*sec2yr
        ean = 22.2d0
        eah = 65.9d0
        eaoh = 17.9d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
    
    case('ab')
        mh = 0.457d0
        moh = -0.572d0
        kinn_ref = 10d0**(-12.56d0)*sec2yr
        kinh_ref = 10d0**(-10.16d0)*sec2yr
        kinoh_ref = 10d0**(-15.6d0)*sec2yr
        ean = 69.8d0
        eah = 65d0
        eaoh = 71d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('kfs')
        mh = 0.5d0
        moh = -0.823d0
        kinn_ref = 10d0**(-12.41d0)*sec2yr
        kinh_ref = 10d0**(-10.06d0)*sec2yr
        kinoh_ref = 10d0**(-9.68d0)*sec2yr*kw**(-moh)
        ean = 9.08*cal2j
        eah = 12.4d0*cal2j
        eaoh = 22.5d0*cal2j
        tc_ref = 25d0
        ! from Brantley et al 2008
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('fo')
        mh = 0.47d0
        moh = 0d0
        kinn_ref = 10d0**(-10.64d0)*sec2yr
        kinh_ref = 10d0**(-6.85d0)*sec2yr
        kinoh_ref = 0d0
        ean = 79d0
        eah = 67.2d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('fa')
        mh = 1d0
        moh = 0d0
        kinn_ref = 10d0**(-12.80d0)*sec2yr
        kinh_ref = 10d0**(-4.80d0)*sec2yr
        kinoh_ref = 0d0
        ean = 94.4d0
        eah = 94.4d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('an')
        mh = 1.411d0
        moh = 0d0
        kinn_ref = 10d0**(-9.12d0)*sec2yr
        kinh_ref = 10d0**(-3.5d0)*sec2yr
        kinoh_ref = 0d0
        ean = 17.8d0
        eah = 16.6d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('la')
        mh = 0.626d0
        moh = 0d0
        kinn_ref = 10d0**(-10.91d0)*sec2yr
        kinh_ref = 10d0**(-7.87d0)*sec2yr
        kinoh_ref = 0d0
        ean = 45.2d0
        eah = 42.1d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('cc')
        mh = 1d0
        moh = 0d0
        kinn_ref = 10d0**(-5.81d0)*sec2yr
        kinh_ref = 10d0**(-0.3d0)*sec2yr
        kinoh_ref = 0d0
        ean = 23.5d0
        eah = 14.4d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('arg')
        ! assumed to be the same as those for cc
        mh = 1d0
        moh = 0d0
        kinn_ref = 10d0**(-5.81d0)*sec2yr
        kinh_ref = 10d0**(-0.3d0)*sec2yr
        kinoh_ref = 0d0
        ean = 23.5d0
        eah = 14.4d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('dlm') ! for disordered dolomite
        mh = 0.500d0
        moh = 0d0
        kinn_ref = 10d0**(-7.53d0)*sec2yr
        kinh_ref = 10d0**(-3.19d0)*sec2yr
        kinoh_ref = 0d0
        ean = 52.2d0
        eah = 36.1d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('gb')
        mh = 0.992d0
        moh = -0.784d0
        kinn_ref = 10d0**(-11.50d0)*sec2yr
        kinh_ref = 10d0**(-7.65d0)*sec2yr
        kinoh_ref = 10d0**(-16.65d0)*sec2yr
        ean = 61.2d0
        eah = 47.5d0
        eaoh = 80.1d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('amsi')
        mh = 0d0
        moh = 0d0
        kinn_ref = 10d0**(-12.23d0)*sec2yr
        kinh_ref = 0d0
        kinoh_ref = 0d0
        ean = 74.5d0
        eah = 0d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        dkin_dmsp = 0d0

    case('qtz')
        mh = 0d0
        moh = 0d0
        kinn_ref = 10d0**(-13.40d0)*sec2yr
        kinh_ref = 0d0
        kinoh_ref = 0d0
        ean = 90.9d0
        eah = 0d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        dkin_dmsp = 0d0

    case('gt')
        mh = 0d0
        moh = 0d0
        kinn_ref = 10d0**(-7.94d0)*sec2yr
        kinh_ref = 0d0
        kinoh_ref = 0d0
        ean = 86.5d0
        eah = 0d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        dkin_dmsp = 0d0

    case('hm')
        mh = 1d0
        moh = 0d0
        kinn_ref = 10d0**(-14.60d0)*sec2yr
        kinh_ref = 10d0**(-9.39d0)*sec2yr
        kinoh_ref = 0d0
        ean = 66.2d0
        eah = 66.2d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        dkin_dmsp = 0d0

    case('ct')
        mh = 0d0
        moh = -0.23d0
        kinn_ref = 10d0**(-12d0)*sec2yr
        kinh_ref = 0d0
        kinoh_ref = 10d0**(-13.58d0)*sec2yr
        ean = 73.5d0
        eah = 0d0
        eaoh = 73.5d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('cabd','ill') ! illite kinetics is assumed to be the same as smectite (Bibi et al., 2011)
        mh = 0.34d0
        moh = -0.4d0
        kinn_ref = 10d0**(-12.78d0)*sec2yr
        kinh_ref = 10d0**(-10.98d0)*sec2yr
        kinoh_ref = 10d0**(-16.52d0)*sec2yr
        ean = 35d0
        eah = 23.6d0
        eaoh = 58.9d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('nph','anl') ! analcime kinetics is assumed to be the same as nepherine (cf. Ragnarsdottir, GCA, 1993)
        mh = 1.130d0
        moh = -0.200d0
        kinn_ref = 10d0**(-8.56d0)*sec2yr
        kinh_ref = 10d0**(-2.73d0)*sec2yr
        kinoh_ref = 10d0**(-10.76d0)*sec2yr
        ean = 65.4d0
        eah = 62.9d0
        eaoh = 37.8d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 

    case('dp')
        mh = 0.71d0
        moh = 0d0
        kinn_ref = 10d0**(-11.11d0)*sec2yr
        kinh_ref = 10d0**(-6.36d0)*sec2yr
        kinoh_ref = 0d0
        ean = 50.6d0
        eah = 96.1d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
    
    case('hb')
        mh = 0.70d0
        moh = 0d0
        kinn_ref = 10d0**(-11.97d0)*sec2yr
        kinh_ref = 10d0**(-6.82d0)*sec2yr
        kinoh_ref = 0d0
        ean = 78.0d0
        eah = 78.0d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
    
    case('tm')
        mh = 0.70d0
        moh = 0d0
        kinn_ref = 10d0**(-10.60d0)*sec2yr
        kinh_ref = 10d0**(-8.40d0)*sec2yr
        kinoh_ref = 0d0
        ean = 94.4d0
        eah = 18.9d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
    
    case('gps')
        mh = 0d0
        moh = 0d0
        kinn_ref = 10d0**(-2.79d0)*sec2yr
        kinh_ref = 0d0
        kinoh_ref = 0d0
        ean = 0d0
        eah = 0d0
        eaoh = 0d0
        tc_ref = 25d0
        ! from Palandri and Kharaka, 2004
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
        
    case('py')
        mh = 0d0
        moh = -0.11d0
        kinn_ref = 0d0
        kinh_ref = 0d0
        kinoh_ref = 10.0d0**(-8.19d0)*sec2yr*kho**0.5d0
        ean = 0d0
        eah = 0d0
        eaoh = 57d0
        tc_ref = 15d0
        ! from Williamson and Rimstidt (1994)
        kin = ( & 
            & k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) &
            & + prox**mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
            & + prox**moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
            & ) 
        select case(trim(adjustl(dev_sp)))
            case('pro')
                dkin_dmsp = ( & 
                    & + mh*prox**(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) &
                    & + moh*prox**(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) &
                    & ) 
            case default 
                dkin_dmsp = 0d0
        endselect 
        
    case('g1')
        kin = ( &
            & 1d0/(poro*hr*mv_tmp*1d-6) &
            & *1d0/1d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
            & )
        dkin_dmsp = 0d0
        
    case('g2')
        kin = ( &
            & 1d0/(poro*hr*mv_tmp*1d-6) &
            & *1d0/8d0 &! mol m^-2 yr^-1, just a value assumed; turnover time of 8 year as in Chen et al. (2010, AFM) 
            & )
        dkin_dmsp = 0d0
        
    case('g3')
        kin = ( &
            & 1d0/(poro*hr*mv_tmp*1d-6) &
            & *1d0/1d3 &! mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1k year  
            & )
        dkin_dmsp = 0d0
        
    case default 
        kin =0d0
        dkin_dmsp = 0d0

endselect  


endsubroutine sld_kin

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine sld_therm( &
    & rg,tc,tempk_0,ss_x &! input
    & ,mineral &! input
    & ,therm &! output
    & ) 
implicit none

real(kind=8),intent(in)::rg,tc,tempk_0,ss_x
real(kind=8) :: cal2j = 4.184d0 
real(kind=8),intent(out):: therm
character(5),intent(in):: mineral
real(kind=8) tc_ref,ha,therm_ref,delG
real(kind=8) tc_ref_1,ha_1,therm_ref_1,therm_1,delG_1
real(kind=8) tc_ref_2,ha_2,therm_ref_2,therm_2,delG_2

! real(kind=8) k_arrhenius

therm = 0d0

select case(trim(adjustl(mineral))) 
    case('ka') 
        ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
        therm_ref = 10d0**(7.435d0)
        ha = -35.3d0*cal2j
        tc_ref = 25d0
        ! from PHREEQC.DAT 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    ! case('ab')
        ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
        ! therm_ref = 10d0**3.412182823d0
        ! ha = -54.15042876d0
        ! tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('kfs')
        ! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
        therm_ref = 10d0**0.227294204d0
        ha = -26.30862098d0
        tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('anl')
        ! NaAlSi2O6*H2O  + 5 H2O  = Na+  + Al(OH)4-  + 2 Si(OH)4(aq)
        therm_ref = 10d0**(-16.06d0)
        ha = 101d0
        tc_ref = 25d0
        ! from Wilkin and Barnes 1998
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('nph')
        ! Nepheline  + 4 H+  = 2 H2O  + SiO2(aq)  + Al+++  + Na+
        therm_ref = 10d0**(14.93646757d0)
        ha = -130.8197467d0
        tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('fo')
        ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O
        therm_ref = 10d0**29.41364324d0
        ha = -208.5932252d0
        tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('fa')
        ! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O
        therm_ref = 10d0**19.98781342d0
        ha = -153.7676621d0
        tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    ! case('an')
        ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
        ! therm_ref = 10d0**28.8615308d0
        ! ha = -292.8769275d0
        ! tc_ref = 15d0
        ! from Kanzaki and Murakami 2018
        ! therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('cc')
        ! CaCO3 = Ca2+ + CO32-
        therm_ref = 10d0**(-8.43d0)
        ha = -8.028943471d0
        tc_ref = 15d0
        ! from Kanzaki and Murakami 2015
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('arg')
        ! CaCO3 = Ca2+ + CO32-
        therm_ref = 10d0**(-8.3d0)
        ha = -12d0
        tc_ref = 25d0
        ! from minteq.v4
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('dlm') ! disordered
        ! CaMg(CO3)2 = Ca+2 + Mg+2 + 2CO3-2
        therm_ref = 10d0**(-16.54d0)
        ha = -46.4d0
        tc_ref = 25d0
        ! from minteq.v4
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('gb')
        ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O
        therm_ref = 10d0**(8.11d0)
        ha = -22.80d0*cal2j
        tc_ref = 25d0
        ! from PHREEQC.DAT 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('amsi')
        ! SiO2 + 2 H2O = H4SiO4
        therm_ref = 10d0**(-2.71d0)
        ha = 3.340d0*cal2j
        tc_ref = 25d0
        ! from PHREEQC.DAT 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('qtz')
        ! SiO2 + 2H2O = H4SiO4
        therm_ref = 10d0**(-4d0)
        ha = 22.36d0
        tc_ref = 25d0
        ! from minteq.v4 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('gt')
        ! Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
        therm_ref = 10d0**(0.5345d0)
        ha = -61.53703d0
        tc_ref = 25d0
        ! from Sugimori et al. 2012 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('hm')
        ! Fe2O3 + 6H+ = 2Fe+3 + 3H2O
        therm_ref = 10d0**(-1.418d0)
        ha = -128.987d0
        tc_ref = 25d0
        ! from minteq.v4
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('ct')
        ! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
        therm_ref = 10d0**(32.2d0)
        ha = -46.800d0*cal2j
        tc_ref = 25d0
        ! from PHREEQC.DAT 
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('cabd')
        ! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
        therm_ref = 10d0**(7.269946518d0)
        ha = -157.0186168d0
        tc_ref = 15d0
        ! from Kanzaki & Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('ill')
        ! Illite  + 8 H+  = 5 H2O  + .6 K+  + .25 Mg++  + 2.3 Al+++  + 3.5 SiO2(aq)
        therm_ref = 10d0**(10.8063184d0)
        ha = -166.39733d0
        tc_ref = 15d0
        ! from Kanzaki & Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('dp')
        ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
        therm_ref = 10d0**(21.79853309d0)
        ha = -138.6020832d0
        tc_ref = 15d0
        ! from Kanzaki & Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('hb')
        ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
        therm_ref = 10d0**(20.20981116d0)
        ha = -128.5d0
        tc_ref = 15d0
        ! from Kanzaki & Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('tm')
        ! Tremolite  + 14 H+  = 8 H2O  + 8 SiO2(aq)  + 2 Ca++  + 5 Mg++
        therm_ref = 10d0**(61.6715d0)
        ha = -429.0d0
        tc_ref = 15d0
        ! from Kanzaki & Murakami 2018
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('gps')
        ! CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O
        therm_ref = 10d0**(-4.61d0)
        ha = 1d0
        tc_ref = 25d0
        ! from minteq.v4
        therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg)
    case('la','ab','an')
        ! CaxNa(1-x)Al(1+x)Si(3-x)O8 + (4x + 4) = xCa+2 + (1-x)Na+ + (1+x)Al+++ + (3-x)SiO2(aq) 
        ! obtaining Anorthite 
        therm_ref_1 = 10d0**28.8615308d0
        ha_1 = -292.8769275d0
        tc_ref_1 = 15d0
        ! from Kanzaki and Murakami 2018
        therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg) ! rg in kJ mol^-1 K^-1
        delG_1 = - rg*(tc+tempk_0)*log(therm_1) ! del-G = -RT ln K  now in kJ mol-1
        ! Then albite 
        therm_ref_2 = 10d0**3.412182823d0
        ha_2 = -54.15042876d0
        tc_ref_2 = 15d0
        ! from Kanzaki and Murakami 2018
        therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg)
        delG_2 = - rg*(tc+tempk_0)*log(therm_2) ! del-G = -RT ln K  now in kJ mol-1
        
        if (ss_x == 1d0) then 
            delG = delG_1 ! ideal anorthite
        elseif (ss_x == 0d0) then 
            delG = delG_2 ! ideal albite
        elseif (ss_x > 0d0 .and. ss_x < 1d0) then  ! solid solution 
            ! ideal(?) mixing (after Gislason and Arnorsson, 1993)
            delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x))
        endif 
        therm = exp(-delG/(rg*(tc+tempk_0)))
    case('g1')
        therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
    case('g2')
        therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
    case('g3')
        therm = 0.121d0 ! mo2 Michaelis, Davidson et al. (2012)
    case default 
        therm = 0d0
endselect 

endsubroutine sld_therm

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
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
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
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42
real(kind=8),dimension(nz)::nax,mgx,cax,so4x,pco2x,six,alx,fe2x,fe3x,kx,so4f
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
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

logical,intent(in)::print_cb
character(500),intent(in)::print_loc
logical so4_error

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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3  = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)

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
    so4x = maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = maqc(findloc(chraq_cnst,'so4',dim=1),:)
endif 

! netcat = nax + kx -2d0*so4x
netcat = kx -2d0*so4x

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


call calc_so4( &
    & nz,so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox &! input 
    & ,nsp_gas_all,nsp_aq_all,chraq_all,chrgas_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input 
    & ,so4f,so4_error &! output
    & )


! print*,'calc_pH'
do while (error > tol)
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        ! na
        & + nax*prox**2d0/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) &
        & -1d0*nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) &
        &       *k1naco3*k1*k2*kco2*pco2x  &
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
        ! na
        & + nax*prox*2d0/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) &
        & + nax*prox**2d0*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &       *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
        & -1d0*nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox)**2d0 &
        &       *k1naco3*k1*k2*kco2*pco2x  &
        &       *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
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
    write(88,*) ' z ',' h+ ',' oh- ',' na+ ',' naco3- ', ' nahco3 ', ' k+ ',' so42- ', 'hco3- ', ' co32- ' &
        & ,' h4sio4 ',' h3sio4- ',' h2sio42- ' &
        & ,' mg2+ ', ' mg(oh)+ ', ' mgco3 ', 'mghco3+ ', ' ca2+ ', ' ca(oh)+ ', ' caco3 ', ' cahco3+ ' &
        & ,' al3+ ', ' al(oh)2+ ', ' al(oh)2+ ', ' al(oh)3 ', ' al(oh)4- ' &
        & , ' fe22+ ', ' fe2(oh)+ ', ' fe2co3 ', ' fe2hco3+ ' &
        & ,' fe3+ ', ' fe3(oh)2+ ', ' fe3(oh)2+ ', ' fe3(oh)3 ', ' fe3(oh)4- ', ' total_charge ' 
    do iz=1,nz
        write(88,*) z(iz) &
        &, prox(iz) &
        & ,kw/prox(iz) &
        & ,nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & ,nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)     &
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
        & +(1d0)*nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        & +(-1d0)*nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*nax(iz)/(1d0+k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0+k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)) &
        &       *k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)      &
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

subroutine calc_pH_v6( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input
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
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42
real(kind=8),dimension(nz)::nax,mgx,cax,so4x,pco2x,six,alx,fe2x,fe3x,kx
real(kind=8),dimension(nz)::so4f,naf,kf,mgf,caf,fe2f,fe3f,alf,sif
real(kind=8),dimension(nz)::dso4f_dpro,dnaf_dpro,dkf_dpro,dmgf_dpro,dcaf_dpro,dfe2f_dpro,dfe3f_dpro,dalf_dpro,dsif_dpro
real(kind=8),dimension(nz),intent(in)::z
real(kind=8),dimension(nz),intent(inout)::prox
logical,intent(out)::ph_error

real(kind=8),dimension(nz)::df,f,netcat
real(kind=8) error,tol,dconc
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
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

logical,intent(in)::print_cb
character(500),intent(in)::print_loc
logical so4_error

error = 1d4
tol = 1d-6
dconc = 1d-14

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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3  = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)

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
    so4x = maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = maqc(findloc(chraq_cnst,'so4',dim=1),:)
endif 

! netcat = nax + kx -2d0*so4x
! netcat = kx -2d0*so4x
netcat = 0d0

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
    ! free SO42- (for simplicity only consider XSO4 complex where X is a cation)
    
    if (any(so4x==0d0)) then 
        so4f = so4x
        dso4f_dpro = 0d0
    else 
        ! print*,'so4x'
        ! print*,so4x
        call calc_so4( &
            & nz,so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox &! input 
            & ,nsp_gas_all,nsp_aq_all,chraq_all,chrgas_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input 
            & ,so4f,so4_error &! output
            & )
        call calc_so4( &
            & nz,so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox+dconc &! input 
            & ,nsp_gas_all,nsp_aq_all,chraq_all,chrgas_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input 
            & ,dso4f_dpro,so4_error &! output
            & )
        dso4f_dpro = (dso4f_dpro - so4f)/dconc 
    endif 
    
    kf = kx/(1d0+k1kso4*so4f)
    dkf_dpro = kx*(-1d0)/(1d0+k1kso4*so4f)**2d0*(k1kso4*dso4f_dpro)
    
    naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    dnaf_dpro = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
        & *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0+k1naso4*dso4f_dpro)
        
    caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    dcaf_dpro = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
        & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0 &
        &   +k1caso4*dso4f_dpro)
        
    mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    dmgf_dpro = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
        & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0 &
        &   +k1mgso4*dso4f_dpro)
        
    fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    dfe2f_dpro = fe2x*(-1d0) &
        & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
        & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0 &
        &   +k1fe2so4*dso4f_dpro)
        
    sif = six/(1d0+k1si/prox+k2si/prox**2d0)
    dsif_dpro = six*(-1d0)/(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0)
    
    alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)
    dalf_dpro = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
        & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0+k1also4*dso4f_dpro)
    
    fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)
    dfe3f_dpro = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
        & *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0+k1fe3so4*dso4f_dpro)
    
    f = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        ! so4
        & -2d0*so4f*prox**2d0 &
        & -1d0*so4f*prox*k1so4 &
        ! k
        & + kf*prox**2d0 &
        & - kf*prox**2d0*k1kso4*so4f &
        ! na
        & + naf*prox**2d0 &
        & -1d0*naf*k1naco3*k1*k2*kco2*pco2x  &
        & -1d0*naf*prox**2d0*k1naso4*so4f  &
        ! si 
        & - sif*k1si*prox &
        & - 2d0*sif*k2si &
        ! mg
        & + 2d0*mgf*prox**2d0 &
        & + mgf*k1mg*prox &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x*prox  &
        ! ca
        & + 2d0*caf*prox**2d0 &
        & + caf*k1ca*prox &
        & + caf*k1cahco3*k1*k2*kco2*pco2x*prox &
        ! al
        & + 3d0*alf*prox**2d0 &
        & + 2d0*alf*k1al*prox &
        & + alf*k2al &
        & - alf*k4al/prox**2d0 &
        & + alf*k1also4*so4f*prox**2d0 &
        ! fe2
        & + 2d0*fe2f*prox**2d0 &
        & + fe2f*k1fe2*prox &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        ! fe3
        & + 3d0*fe3f*prox**2d0 &
        & + 2d0*fe3f*k1fe3*prox &
        & + fe3f*k2fe3 &
        & - fe3f*k4fe3/prox**2d0 &
        & + fe3f*k1fe3so4*so4f*prox**2d0 
    df = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw) &
        !
        ! so4
        & -2d0*dso4f_dpro*prox**2d0 &
        & -2d0*so4f*prox*2d0 &
        & -1d0*dso4f_dpro*prox*k1so4 &
        & -1d0*so4f*k1so4 &
        !
        ! k
        & + dkf_dpro*prox**2d0 &
        & + kf*prox*2d0 &
        & - dkf_dpro*prox**2d0*k1kso4*so4f &
        & - kf*prox*2d0*k1kso4*so4f &
        & - kf*prox**2d0*k1kso4*dso4f_dpro &
        !
        ! na
        & + dnaf_dpro*prox**2d0 &
        & + naf*prox*2d0 &
        & -1d0*dnaf_dpro*k1naco3*k1*k2*kco2*pco2x  &
        & -1d0*dnaf_dpro*prox**2d0*k1naso4*so4f  &
        & -1d0*naf*prox*2d0*k1naso4*so4f  &
        & -1d0*naf*prox**2d0*k1naso4*dso4f_dpro  &
        !
        ! si
        !
        & - dsif_dpro*k1si*prox &
        & - sif*k1si*1d0 &
        & - 2d0*dsif_dpro*k2si &
        !
        ! mg
        !
        & + 2d0*dmgf_dpro*prox**2d0 &
        & + 2d0*mgf*prox*2d0 &
        & + dmgf_dpro*k1mg*prox &
        & + mgf*k1mg*1d0 &
        & + dmgf_dpro*k1mghco3*k1*k2*kco2*pco2x*prox  &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x*1d0  &
        ! 
        ! ca
        !
        & + 2d0*dcaf_dpro*prox**2d0 &
        & + 2d0*caf*prox*2d0 &
        & + dcaf_dpro*k1ca*prox &
        & + caf*k1ca*1d0 &
        & + dcaf_dpro*k1cahco3*k1*k2*kco2*pco2x*prox  &
        & + caf*k1cahco3*k1*k2*kco2*pco2x*1d0  &
        !
        ! al
        !
        & + 3d0*dalf_dpro*prox**2d0 &
        & + 3d0*alf*prox*2d0 &
        & + 2d0*dalf_dpro*k1al*prox &
        & + 2d0*alf*k1al*1d0 &
        & + dalf_dpro*k2al &
        & - dalf_dpro*k4al/prox**2d0 &
        & - alf*k4al*(-2d0)/prox**3d0 &
        & + dalf_dpro*k1also4*so4f*prox**2d0 &
        & + alf*k1also4*dso4f_dpro*prox**2d0 &
        & + alf*k1also4*so4f*prox*2d0 &
        !
        ! fe2
        ! 
        & + 2d0*dfe2f_dpro*prox**2d0 &
        & + 2d0*fe2f*prox*2d0 &
        & + dfe2f_dpro*k1fe2*prox &
        & + fe2f*k1fe2*1d0 &
        & + dfe2f_dpro*k1fe2hco3*k1*k2*kco2*pco2x*prox  &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x*1d0  &
        !
        ! fe3 
        !
        & + 3d0*dfe3f_dpro*prox**2d0 &
        & + 3d0*fe3f*prox*2d0 &
        & + 2d0*dfe3f_dpro*k1fe3*prox &
        & + 2d0*fe3f*k1fe3*1d0 &
        & + dfe3f_dpro*k2fe3 &
        & - dfe3f_dpro*k4fe3/prox**2d0 &
        & - fe3f*k4fe3*(-2d0)/prox**3d0 &
        & + dfe3f_dpro*k1fe3so4*so4f*prox**2d0 &
        & + fe3f*k1fe3so4*dso4f_dpro*prox**2d0 &
        & + fe3f*k1fe3so4*so4f*prox*2d0 
        
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

if (any(so4x==0d0)) then 
    so4f = so4x
else 
    call calc_so4( &
        & nz,so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox &! input 
        & ,nsp_gas_all,nsp_aq_all,chraq_all,chrgas_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input 
        & ,so4f,so4_error &! output
        & )
endif 
    
kf = kx/(1d0+k1kso4*so4f)

naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    
caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    
mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    
fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    
sif = six/(1d0+k1si/prox+k2si/prox**2d0)

alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)

fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)

if (print_cb) then 
    open(88,file = trim(adjustl(print_loc)),status='replace')
    write(88,*) ' z ',' h+ ',' oh- ',' na+ ',' naco3- ', ' nahco3 ', ' naso4- ', ' k+ ',' kso4- ' &
        & ,' so42- ', ' hso4- ', 'hco3- ', ' co32- ' &
        & ,' h4sio4 ',' h3sio4- ',' h2sio42- ' &
        & ,' mg2+ ', ' mg(oh)+ ', ' mgco3 ', 'mghco3+ ', 'mgso4 ' &
        & , ' ca2+ ', ' ca(oh)+ ', ' caco3 ', ' cahco3+ ', ' caso4 ' &
        & ,' al3+ ', ' al(oh)2+ ', ' al(oh)2+ ', ' al(oh)3 ', ' al(oh)4- ' , ' also4+ ' &
        & , ' fe22+ ', ' fe2(oh)+ ', ' fe2co3 ', ' fe2hco3+ ', ' fe2so4 ' &
        & ,' fe3+ ', ' fe3(oh)2+ ', ' fe3(oh)2+ ', ' fe3(oh)3 ', ' fe3(oh)4- ', ' fe3so4+ ' &
        & , ' total_charge ' 
    do iz=1,nz
        write(88,*) z(iz) &
        &, prox(iz) &
        & ,kw/prox(iz) &
        & ,naf(iz) &
        & ,naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)     &
        & ,naf(iz)*k1naso4*so4f(iz)     &
        & ,kf(iz) &
        & ,kf(iz)*k1kso4*so4f(iz) &
        & ,so4f(iz) &
        & ,so4f(iz)*k1so4/prox &
        & ,k1*kco2*pco2x(iz)/prox(iz) &
        & ,k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,sif(iz) &
        & ,sif(iz)*k1si/prox(iz) &
        & ,sif(iz)*k2si/prox(iz)**2d0 &
        & ,mgf(iz) &
        & ,mgf(iz)*k1mg/prox(iz)  &
        & ,mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & ,mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & ,mgf(iz)*k1mgso4*so4f(iz)  &
        & ,caf(iz) &
        & ,caf(iz)*k1ca/prox(iz) &
        & ,caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,caf(iz)*k1caso4*so4f(iz) &
        & ,alf(iz) &
        & ,alf(iz)*k1al/prox(iz) &
        & ,alf(iz)*k2al/prox(iz)**2d0 &
        & ,alf(iz)*k3al/prox(iz)**3d0 &
        & ,alf(iz)*k4al/prox(iz)**4d0 &
        & ,alf(iz)*k1also4*so4f(iz) &
        & ,fe2f(iz) &
        & ,fe2f(iz)*k1fe2/prox(iz) &
        & ,fe2f(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,fe2f(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,fe2f(iz)*k1fe2so4*so4f(iz) &
        & ,fe3f(iz) &
        & ,fe3f(iz)*k1fe3/prox(iz) &
        & ,fe3f(iz)*k2fe3/prox(iz)**2d0 &
        & ,fe3f(iz)*k3fe3/prox(iz)**3d0 &
        & ,fe3f(iz)*k4fe3/prox(iz)**4d0 &
        & ,fe3f(iz)*k1fe3so4*so4f(iz) &
        ! charge balance 
        & ,1d0*prox(iz) &
        & +(-1d0)*kw/prox(iz) &
        & +(1d0)*naf(iz) &
        & +(-1d0)*naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)      &
        & +(-1d0)*naf(iz)*k1naso4*so4f(iz)     &
        & +(1d0)*kf(iz) &
        & +(-1d0)*kf(iz)*k1kso4*so4f(iz) &
        & +(-2d0)*so4f(iz) &
        & +(-1d0)*so4f(iz)*k1so4/prox &
        & +(-1d0)*k1*kco2*pco2x(iz)/prox(iz) &
        & +(-2d0)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*sif(iz) &
        & +(-1d0)*sif(iz)*k1si/prox(iz) &
        & +(-2d0)*sif(iz)*k2si/prox(iz)**2d0 &
        & +(2d0)*mgf(iz) &
        & +(1d0)*mgf(iz)*k1mg/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & +(1d0)*mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgso4*so4f(iz)  &
        & +(2d0)*caf(iz) &
        & +(1d0)*caf(iz)*k1ca/prox(iz) &
        & +(0d0)*caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*caf(iz)*k1caso4*so4f(iz) &
        & +(3d0)*alx(iz) &
        & +(2d0)*alx(iz)*k1al/prox(iz) &
        & +(1d0)*alx(iz)*k2al/prox(iz)**2d0 &
        & +(0d0)*alx(iz)*k3al/prox(iz)**3d0 &
        & +(-1d0)*alx(iz)*k4al/prox(iz)**4d0 &
        & +(1d0)*alx(iz)*k1also4*so4f(iz) &
        & +(2d0)*fe2x(iz) &
        & +(1d0)*fe2x(iz)*k1fe2/prox(iz) &
        & +(0d0)*fe2x(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*fe2x(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*fe2x(iz)*k1fe2so4*so4f(iz) &
        & +(3d0)*fe3x(iz) &
        & +(2d0)*fe3x(iz)*k1fe3/prox(iz) &
        & +(1d0)*fe3x(iz)*k2fe3/prox(iz)**2d0 &
        & +(0d0)*fe3x(iz)*k3fe3/prox(iz)**3d0 &
        & +(-1d0)*fe3x(iz)*k4fe3/prox(iz)**4d0 &
        & +(-1d0)*fe3x(iz)*k1fe3so4*so4f(iz) 
    enddo 
    close(88)
endif 
            

endsubroutine calc_pH_v6

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH_v7( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
    & ,print_cb,print_loc,z &! input 
    & ,prox,ph_error,so4f &! output
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
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42,so4th
real(kind=8),dimension(nz)::nax,mgx,cax,so4x,pco2x,six,alx,fe2x,fe3x,kx
real(kind=8),dimension(nz)::naf,kf,mgf,caf,fe2f,fe3f,alf,sif
real(kind=8),dimension(nz)::dnaf_dpro,dkf_dpro,dmgf_dpro,dcaf_dpro,dfe2f_dpro,dfe3f_dpro,dalf_dpro,dsif_dpro
real(kind=8),dimension(nz)::dnaf_dso4f,dkf_dso4f,dmgf_dso4f,dcaf_dso4f,dfe2f_dso4f,dfe3f_dso4f,dalf_dso4f,dsif_dso4f
real(kind=8),dimension(nz),intent(in)::z
real(kind=8),dimension(nz),intent(inout)::prox,so4f
logical,intent(out)::ph_error

real(kind=8),dimension(nz)::df1,f1,netcat,f2,df2,df21,df12
real(kind=8) error,tol,dconc
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
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

logical,intent(in)::print_cb
character(500),intent(in)::print_loc
logical so4_error

real(kind=8),allocatable::amx(:,:),ymx(:)
integer,allocatable::ipiv(:)
integer info,nmx

external DGESV


error = 1d4
tol = 1d-6
dconc = 1d-14

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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3  = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)

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
    so4x = maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = maqc(findloc(chraq_cnst,'so4',dim=1),:)
endif 

! netcat = nax + kx -2d0*so4x
! netcat = kx -2d0*so4x
netcat = 0d0

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


so4th = maqth_all(findloc(chraq_all,'so4',dim=1))

so4f  = so4x 

nmx = 2*nz
! if (all(so4x==0d0)) then 
if (all(so4x<=so4th)) then 
    nmx = nz
endif  

allocate(amx(nmx,nmx),ymx(nmx),ipiv(nmx))

ph_error = .false.

! print*,'calc_pH'
do while (error > tol)
    ! free SO42- (for simplicity only consider XSO4 complex where X is a cation)
    
    kf = kx/(1d0+k1kso4*so4f)
    dkf_dso4f = kx*(-1d0)/(1d0+k1kso4*so4f)**2d0*(k1kso4)
    dkf_dpro = 0d0
    
    naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    dnaf_dpro = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
        & *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dnaf_dso4f = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
        & *(k1naso4)
        
    caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    dcaf_dpro = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
        & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dcaf_dso4f = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
        & *(k1caso4)
        
    mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    dmgf_dpro = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
        & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dmgf_dso4f = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
        & *(k1mgso4)
        
    fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    dfe2f_dpro = fe2x*(-1d0) &
        & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
        & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dfe2f_dso4f = fe2x*(-1d0) &
        & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
        & *(k1fe2so4)
        
    sif = six/(1d0+k1si/prox+k2si/prox**2d0)
    dsif_dpro = six*(-1d0)/(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0)
    dsif_dso4f = 0d0
    
    alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)
    dalf_dpro = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
        & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0)
    dalf_dso4f = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
        & *(k1also4)
    
    fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)
    dfe3f_dpro = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
        & *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0)
    dfe3f_dso4f = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
        & *(k1fe3so4)
    
    f1 = prox**5d0 + netcat*prox**4d0 - (k1*kco2*pco2x+kw)*prox**3d0 - 2d0*k2*k1*kco2*pco2x*prox**2d0  &
        ! so4
        & -2d0*so4f*prox**4d0 &
        & -1d0*so4f*prox**5d0*k1so4 &
        ! k
        & + kf*prox**4d0 &
        & - kf*prox**4d0*k1kso4*so4f &
        ! na
        & + naf*prox**4d0 &
        & -1d0*naf*k1naco3*k1*k2*kco2*pco2x*prox**2d0  &
        & -1d0*naf*prox**4d0*k1naso4*so4f  &
        ! si 
        & - sif*k1si*prox**3d0 &
        & - 2d0*sif*k2si*prox**2d0 &
        ! mg
        & + 2d0*mgf*prox**4d0 &
        & + mgf*k1mg*prox**3d0 &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x*prox**3d0  &
        ! ca
        & + 2d0*caf*prox**4d0 &
        & + caf*k1ca*prox**3d0 &
        & + caf*k1cahco3*k1*k2*kco2*pco2x*prox**3d0 &
        ! al
        & + 3d0*alf*prox**4d0 &
        & + 2d0*alf*k1al*prox**3d0 &
        & + alf*k2al*prox**2d0 &
        & - alf*k4al &
        & + alf*k1also4*so4f*prox**4d0 &
        ! fe2
        & + 2d0*fe2f*prox**4d0 &
        & + fe2f*k1fe2*prox**3d0 &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x*prox**3d0 &
        ! fe3
        & + 3d0*fe3f*prox**4d0 &
        & + 2d0*fe3f*k1fe3*prox**3d0 &
        & + fe3f*k2fe3*prox**2d0 &
        & - fe3f*k4fe3 &
        & + fe3f*k1fe3so4*so4f*prox**4d0 
        
    df1 = 5d0*prox**4d0 + 4d0*netcat*prox**3d0 - (k1*kco2*pco2x+kw)*3d0*prox**2d0 - 2d0*k2*k1*kco2*pco2x*prox*2d0 &
        !
        ! so4
        & -2d0*so4f*4d0*prox**3d0 &
        & -1d0*so4f*5d0*prox**4d0*k1so4 &
        !
        ! k
        & + dkf_dpro*prox**4d0 &
        & + kf*4d0*prox**3d0 &
        & - dkf_dpro*prox**4d0*k1kso4*so4f &
        & - kf*4d0*prox**3d0*k1kso4*so4f &
        !
        ! na
        & + dnaf_dpro*prox**4d0 &
        & + naf*4d0*prox**3d0 &
        & -1d0*dnaf_dpro*k1naco3*k1*k2*kco2*pco2x*prox**2d0  &
        & -1d0*naf*k1naco3*k1*k2*kco2*pco2x*prox*2d0  &
        & -1d0*dnaf_dpro*prox**4d0*k1naso4*so4f  &
        & -1d0*naf*4d0*prox**3d0*k1naso4*so4f  &
        !
        ! si
        !
        & - dsif_dpro*k1si*prox**3d0 &
        & - sif*k1si*3d0*prox**2d0 &
        & - 2d0*dsif_dpro*k2si*prox**2d0 &
        & - 2d0*sif*k2si*prox*2d0 &
        !
        ! mg
        !
        & + 2d0*dmgf_dpro*prox**4d0 &
        & + 2d0*mgf*4d0*prox**3d0 &
        & + dmgf_dpro*k1mg*prox**3d0 &
        & + mgf*k1mg*3d0*prox**2d0 &
        & + dmgf_dpro*k1mghco3*k1*k2*kco2*pco2x*prox**3d0  &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x*3d0*prox**2d0  &
        ! 
        ! ca
        !
        & + 2d0*dcaf_dpro*prox**4d0 &
        & + 2d0*caf*4d0*prox**3d0 &
        & + dcaf_dpro*k1ca*prox**3d0 &
        & + caf*k1ca*3d0*prox**2d0 &
        & + dcaf_dpro*k1cahco3*k1*k2*kco2*pco2x*prox**3d0  &
        & + caf*k1cahco3*k1*k2*kco2*pco2x*3d0*prox**2d0  &
        !
        ! al
        !
        & + 3d0*dalf_dpro*prox**4d0 &
        & + 3d0*alf*4d0*prox**3d0 &
        & + 2d0*dalf_dpro*k1al*prox**3d0 &
        & + 2d0*alf*k1al*3d0*prox**2d0 &
        & + dalf_dpro*k2al*prox**2d0 &
        & + alf*k2al*prox*2d0 &
        & - dalf_dpro*k4al &
        & + dalf_dpro*k1also4*so4f*prox**4d0 &
        & + alf*k1also4*so4f*4d0*prox**3d0 &
        !
        ! fe2
        ! 
        & + 2d0*dfe2f_dpro*prox**4d0 &
        & + 2d0*fe2f*4d0*prox**3d0 &
        & + dfe2f_dpro*k1fe2*prox**3d0 &
        & + fe2f*k1fe2*3d0*prox**2d0 &
        & + dfe2f_dpro*k1fe2hco3*k1*k2*kco2*pco2x*prox**3d0  &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x*3d0*prox**2d0  &
        !
        ! fe3 
        !
        & + 3d0*dfe3f_dpro*prox**4d0 &
        & + 3d0*fe3f*4d0*prox**3d0 &
        & + 2d0*dfe3f_dpro*k1fe3*prox**3d0 &
        & + 2d0*fe3f*k1fe3*3d0*prox**2d0 &
        & + dfe3f_dpro*k2fe3*prox**2d0 &
        & + fe3f*k2fe3*prox*2d0 &
        & - dfe3f_dpro*k4fe3 &
        & + dfe3f_dpro*k1fe3so4*so4f*prox**4d0 &
        & + fe3f*k1fe3so4*so4f*4d0*prox**3d0 
    
    df12 =   &
        ! so4
        & -2d0*1d0*prox**4d0 &
        & -1d0*1d0*prox**5d0*k1so4 &
        ! k
        & + dkf_dso4f*prox**4d0 &
        & - dkf_dso4f*prox**4d0*k1kso4*so4f &
        & - kf*prox**4d0*k1kso4*1d0 &
        ! na
        & + dnaf_dso4f*prox**4d0 &
        & -1d0*dnaf_dso4f*k1naco3*k1*k2*kco2*pco2x*prox**2d0  &
        & -1d0*dnaf_dso4f*prox**4d0*k1naso4*so4f  &
        & -1d0*naf*prox**4d0*k1naso4*1d0  &
        ! si 
        & - dsif_dso4f*k1si*prox**3d0 &
        & - 2d0*dsif_dso4f*k2si*prox**2d0 &
        ! mg
        & + 2d0*dmgf_dso4f*prox**4d0 &
        & + dmgf_dso4f*k1mg*prox**3d0 &
        & + dmgf_dso4f*k1mghco3*k1*k2*kco2*pco2x*prox**3d0  &
        ! ca
        & + 2d0*dcaf_dso4f*prox**4d0 &
        & + dcaf_dso4f*k1ca*prox**3d0 &
        & + dcaf_dso4f*k1cahco3*k1*k2*kco2*pco2x*prox**3d0 &
        ! al
        & + 3d0*dalf_dso4f*prox**4d0 &
        & + 2d0*dalf_dso4f*k1al*prox**3d0 &
        & + dalf_dso4f*k2al*prox**2d0 &
        & - dalf_dso4f*k4al &
        & + dalf_dso4f*k1also4*so4f*prox**4d0 &
        & + alf*k1also4*1d0*prox**4d0 &
        ! fe2
        & + 2d0*dfe2f_dso4f*prox**4d0 &
        & + dfe2f_dso4f*k1fe2*prox**3d0 &
        & + dfe2f_dso4f*k1fe2hco3*k1*k2*kco2*pco2x*prox**3d0 &
        ! fe3
        & + 3d0*dfe3f_dso4f*prox**4d0 &
        & + 2d0*dfe3f_dso4f*k1fe3*prox**3d0 &
        & + dfe3f_dso4f*k2fe3*prox**2d0 &
        & - dfe3f_dso4f*k4fe3 &
        & + dfe3f_dso4f*k1fe3so4*so4f*prox**4d0 &
        & + fe3f*k1fe3so4*1d0*prox**4d0 
        
        
        
    
    f2 = prox**4d0*so4x - prox**4d0*so4f*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )
        
    df2 =  - prox**4d0*1d0*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & ) &
        & - prox**4d0*so4f*(  &
        & +k1kso4*dkf_dso4f &
        & +k1naso4*dnaf_dso4f &
        & +k1caso4*dcaf_dso4f &
        & +k1mgso4*dmgf_dso4f &
        & +k1fe2so4*dfe2f_dso4f &
        & +k1also4*dalf_dso4f &
        & +k1fe3so4*dfe3f_dso4f &
        & )
        
    df21 = 4d0*prox**3d0*so4x- prox**4d0*so4f*( k1so4*1d0 &
        & +k1kso4*dkf_dpro &
        & +k1naso4*dnaf_dpro &
        & +k1caso4*dcaf_dpro &
        & +k1mgso4*dmgf_dpro &
        & +k1fe2so4*dfe2f_dpro &
        & +k1also4*dalf_dpro &
        & +k1fe3so4*dfe3f_dpro &
        & ) &
        & -4d0*prox**3d0*so4f*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )
        
        
        
    df1 = df1*prox
    df21 = df21*prox
    df2 = df2*so4f
    df12 = df12*so4f
    
    if (any(isnan(f1)).or.any(isnan(f2)).or.any(isnan(df1)).or.any(isnan(df2)) &
        & .or.any(isnan(df12)).or.any(isnan(df21))) then 
        print *,any(isnan(f1)),any(isnan(f2)),any(isnan(df1)),any(isnan(df2)) &
            & ,any(isnan(df12)),any(isnan(df21))
        print *,prox
        ! pause 
    endif 
    
    
    if (nmx/=nz) then 
        amx = 0d0
        ymx = 0d0
        
        ymx(1:nz) = f1
        ymx(nz+1:nmx) = f2
        
        do iz=1,nz
            amx(iz,iz)=df1(iz)
            amx(nz+iz,nz+iz)=df2(iz)
            amx(iz,nz+iz)=df12(iz)
            amx(nz+iz,iz)=df21(iz)
        enddo 
        ymx = -ymx
        
        call DGESV(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 
        
        prox = prox*exp( ymx(1:nz) )
        so4f = so4f*exp( ymx(nz+1:nmx) )
        
        error = maxval(abs(exp( ymx )-1d0))
        if (isnan(error)) then 
            error = 1d4
            ph_error = .true.
            exit 
        endif 
    else 
        prox = prox*exp( -f1/df1 )
        error = maxval(abs(exp( -f1/df1 )-1d0))
        if (isnan(error)) error = 1d4
    endif 
    
    iter = iter + 1
    if (iter > 1000) then 
        print *,'iteration exceeds 1000'
        ph_error = .true.
        exit
    endif 
enddo 
if (any(isnan(prox))) then     
    print *, (-log10(prox(iz)),iz=1,nz,nz/5)
    print*,'ph is nan'
    ph_error = .true.
    ! stop
endif 

    
kf = kx/(1d0+k1kso4*so4f)

naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    
caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    
mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    
fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    
sif = six/(1d0+k1si/prox+k2si/prox**2d0)

alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)

fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)

if (print_cb) then 
    open(88,file = trim(adjustl(print_loc)),status='replace')
    write(88,*) ' z ',' h+ ',' oh- ',' na+ ',' naco3- ', ' nahco3 ', ' naso4- ', ' k+ ',' kso4- ' &
        & ,' so42- ', ' hso4- ', 'hco3- ', ' co32- ' &
        & ,' h4sio4 ',' h3sio4- ',' h2sio42- ' &
        & ,' mg2+ ', ' mg(oh)+ ', ' mgco3 ', 'mghco3+ ', 'mgso4 ' &
        & , ' ca2+ ', ' ca(oh)+ ', ' caco3 ', ' cahco3+ ', ' caso4 ' &
        & ,' al3+ ', ' al(oh)2+ ', ' al(oh)2+ ', ' al(oh)3 ', ' al(oh)4- ' , ' also4+ ' &
        & , ' fe22+ ', ' fe2(oh)+ ', ' fe2co3 ', ' fe2hco3+ ', ' fe2so4 ' &
        & ,' fe3+ ', ' fe3(oh)2+ ', ' fe3(oh)2+ ', ' fe3(oh)3 ', ' fe3(oh)4- ', ' fe3so4+ ' &
        & , ' total_charge ' 
    do iz=1,nz
        write(88,*) z(iz) &
        & ,prox(iz) &
        & ,kw/prox(iz) &
        & ,naf(iz) &
        & ,naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)     &
        & ,naf(iz)*k1naso4*so4f(iz)     &
        & ,kf(iz) &
        & ,kf(iz)*k1kso4*so4f(iz) &
        & ,so4f(iz) &
        & ,so4f(iz)*k1so4*prox(iz) &
        & ,k1*kco2*pco2x(iz)/prox(iz) &
        & ,k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,sif(iz) &
        & ,sif(iz)*k1si/prox(iz) &
        & ,sif(iz)*k2si/prox(iz)**2d0 &
        & ,mgf(iz) &
        & ,mgf(iz)*k1mg/prox(iz)  &
        & ,mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & ,mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & ,mgf(iz)*k1mgso4*so4f(iz)  &
        & ,caf(iz) &
        & ,caf(iz)*k1ca/prox(iz) &
        & ,caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,caf(iz)*k1caso4*so4f(iz) &
        & ,alf(iz) &
        & ,alf(iz)*k1al/prox(iz) &
        & ,alf(iz)*k2al/prox(iz)**2d0 &
        & ,alf(iz)*k3al/prox(iz)**3d0 &
        & ,alf(iz)*k4al/prox(iz)**4d0 &
        & ,alf(iz)*k1also4*so4f(iz) &
        & ,fe2f(iz) &
        & ,fe2f(iz)*k1fe2/prox(iz) &
        & ,fe2f(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,fe2f(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,fe2f(iz)*k1fe2so4*so4f(iz) &
        & ,fe3f(iz) &
        & ,fe3f(iz)*k1fe3/prox(iz) &
        & ,fe3f(iz)*k2fe3/prox(iz)**2d0 &
        & ,fe3f(iz)*k3fe3/prox(iz)**3d0 &
        & ,fe3f(iz)*k4fe3/prox(iz)**4d0 &
        & ,fe3f(iz)*k1fe3so4*so4f(iz) &
        ! charge balance 
        & ,1d0*prox(iz) &
        & +(-1d0)*kw/prox(iz) &
        & +(1d0)*naf(iz) &
        & +(-1d0)*naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)      &
        & +(-1d0)*naf(iz)*k1naso4*so4f(iz)     &
        & +(1d0)*kf(iz) &
        & +(-1d0)*kf(iz)*k1kso4*so4f(iz) &
        & +(-2d0)*so4f(iz) &
        & +(-1d0)*so4f(iz)*k1so4*prox(iz) &
        & +(-1d0)*k1*kco2*pco2x(iz)/prox(iz) &
        & +(-2d0)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*sif(iz) &
        & +(-1d0)*sif(iz)*k1si/prox(iz) &
        & +(-2d0)*sif(iz)*k2si/prox(iz)**2d0 &
        & +(2d0)*mgf(iz) &
        & +(1d0)*mgf(iz)*k1mg/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & +(1d0)*mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgso4*so4f(iz)  &
        & +(2d0)*caf(iz) &
        & +(1d0)*caf(iz)*k1ca/prox(iz) &
        & +(0d0)*caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*caf(iz)*k1caso4*so4f(iz) &
        & +(3d0)*alf(iz) &
        & +(2d0)*alf(iz)*k1al/prox(iz) &
        & +(1d0)*alf(iz)*k2al/prox(iz)**2d0 &
        & +(0d0)*alf(iz)*k3al/prox(iz)**3d0 &
        & +(-1d0)*alf(iz)*k4al/prox(iz)**4d0 &
        & +(1d0)*alf(iz)*k1also4*so4f(iz) &
        & +(2d0)*fe2f(iz) &
        & +(1d0)*fe2f(iz)*k1fe2/prox(iz) &
        & +(0d0)*fe2f(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*fe2f(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*fe2f(iz)*k1fe2so4*so4f(iz) &
        & +(3d0)*fe3f(iz) &
        & +(2d0)*fe3f(iz)*k1fe3/prox(iz) &
        & +(1d0)*fe3f(iz)*k2fe3/prox(iz)**2d0 &
        & +(0d0)*fe3f(iz)*k3fe3/prox(iz)**3d0 &
        & +(-1d0)*fe3f(iz)*k4fe3/prox(iz)**4d0 &
        & +(1d0)*fe3f(iz)*k1fe3so4*so4f(iz) 
    enddo 
    close(88)
endif 
            

endsubroutine calc_pH_v7

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_pH_v7_2( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
    & ,print_cb,print_loc,z &! input 
    & ,prox,ph_error,so4f &! output
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
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42,so4th
real(kind=8),dimension(nz)::nax,mgx,cax,so4x,pco2x,six,alx,fe2x,fe3x,kx
real(kind=8),dimension(nz)::naf,kf,mgf,caf,fe2f,fe3f,alf,sif
real(kind=8),dimension(nz)::dnaf_dpro,dkf_dpro,dmgf_dpro,dcaf_dpro,dfe2f_dpro,dfe3f_dpro,dalf_dpro,dsif_dpro
real(kind=8),dimension(nz)::dnaf_dso4f,dkf_dso4f,dmgf_dso4f,dcaf_dso4f,dfe2f_dso4f,dfe3f_dso4f,dalf_dso4f,dsif_dso4f
real(kind=8),dimension(nz),intent(in)::z
real(kind=8),dimension(nz),intent(inout)::prox,so4f
logical,intent(out)::ph_error

real(kind=8),dimension(nz)::df1,f1,netcat,f2,df2,df21,df12
real(kind=8) error,tol,dconc
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
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

logical,intent(in)::print_cb
character(500),intent(in)::print_loc
logical so4_error

real(kind=8),allocatable::amx(:,:),ymx(:)
integer,allocatable::ipiv(:)
integer info,nmx

external DGESV


error = 1d4
tol = 1d-6
dconc = 1d-14

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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3  = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)

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
    so4x = maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = maqc(findloc(chraq_cnst,'so4',dim=1),:)
endif 

! netcat = nax + kx -2d0*so4x
! netcat = kx -2d0*so4x
netcat = 0d0

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


so4th = maqth_all(findloc(chraq_all,'so4',dim=1))

so4f  = so4x 

nmx = 2*nz
! if (all(so4x==0d0)) then 
if (all(so4x<=so4th)) then 
    nmx = nz
endif  

allocate(amx(nmx,nmx),ymx(nmx),ipiv(nmx))

ph_error = .false.

! print*,'calc_pH'
do while (error > tol)
    ! free SO42- (for simplicity only consider XSO4 complex where X is a cation)
    
    kf = kx/(1d0+k1kso4*so4f)
    dkf_dso4f = kx*(-1d0)/(1d0+k1kso4*so4f)**2d0*(k1kso4)
    dkf_dpro = 0d0
    
    naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    dnaf_dpro = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
        & *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dnaf_dso4f = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
        & *(k1naso4)
        
    caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    dcaf_dpro = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
        & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dcaf_dso4f = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
        & *(k1caso4)
        
    mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    dmgf_dpro = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
        & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dmgf_dso4f = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
        & *(k1mgso4)
        
    fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    dfe2f_dpro = fe2x*(-1d0) &
        & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
        & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
    dfe2f_dso4f = fe2x*(-1d0) &
        & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
        & *(k1fe2so4)
        
    sif = six/(1d0+k1si/prox+k2si/prox**2d0)
    dsif_dpro = six*(-1d0)/(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0)
    dsif_dso4f = 0d0
    
    alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)
    dalf_dpro = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
        & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0)
    dalf_dso4f = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
        & *(k1also4)
    
    fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)
    dfe3f_dpro = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
        & *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0)
    dfe3f_dso4f = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
        & *(k1fe3so4)
    
    f1 = prox**3d0 + netcat*prox**2d0 - (k1*kco2*pco2x+kw)*prox - 2d0*k2*k1*kco2*pco2x  &
        ! so4
        & -2d0*so4f*prox**2d0 &
        & -1d0*so4f*prox**3d0*k1so4 &
        ! k
        & + kf*prox**2d0 &
        & - kf*prox**2d0*k1kso4*so4f &
        ! na
        & + naf*prox**2d0 &
        & -1d0*naf*k1naco3*k1*k2*kco2*pco2x  &
        & -1d0*naf*prox**2d0*k1naso4*so4f  &
        ! si 
        & - sif*k1si*prox &
        & - 2d0*sif*k2si &
        ! mg
        & + 2d0*mgf*prox**2d0 &
        & + mgf*k1mg*prox &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x*prox  &
        ! ca
        & + 2d0*caf*prox**2d0 &
        & + caf*k1ca*prox &
        & + caf*k1cahco3*k1*k2*kco2*pco2x*prox &
        ! al
        & + 3d0*alf*prox**2d0 &
        & + 2d0*alf*k1al*prox &
        & + alf*k2al &
        & - alf*k4al/prox**2d0 &
        & + alf*k1also4*so4f*prox**2d0 &
        ! fe2
        & + 2d0*fe2f*prox**2d0 &
        & + fe2f*k1fe2*prox &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        ! fe3
        & + 3d0*fe3f*prox**2d0 &
        & + 2d0*fe3f*k1fe3*prox &
        & + fe3f*k2fe3 &
        & - fe3f*k4fe3/prox**2d0 &
        & + fe3f*k1fe3so4*so4f*prox**2d0 
        
    df1 = 3d0*prox**2d0 + 2d0*netcat*prox - (k1*kco2*pco2x+kw)*1d0  &
        ! so4
        & -2d0*so4f*prox*2d0 &
        & -1d0*so4f*3d0*prox**2d0*k1so4 &
        ! k
        & + dkf_dpro*prox**2d0 &
        & + kf*prox*2d0 &
        & - dkf_dpro*prox**2d0*k1kso4*so4f &
        & - kf*prox*2d0*k1kso4*so4f &
        ! na
        & + dnaf_dpro*prox**2d0 &
        & + naf*prox*2d0 &
        & -1d0*dnaf_dpro*k1naco3*k1*k2*kco2*pco2x  &
        & -1d0*dnaf_dpro*prox**2d0*k1naso4*so4f  &
        & -1d0*naf*prox*2d0*k1naso4*so4f  &
        ! si 
        & - dsif_dpro*k1si*prox &
        & - sif*k1si &
        & - 2d0*dsif_dpro*k2si &
        ! mg
        & + 2d0*dmgf_dpro*prox**2d0 &
        & + 2d0*mgf*prox*2d0 &
        & + dmgf_dpro*k1mg*prox &
        & + mgf*k1mg &
        & + dmgf_dpro*k1mghco3*k1*k2*kco2*pco2x*prox  &
        & + mgf*k1mghco3*k1*k2*kco2*pco2x  &
        ! ca
        & + 2d0*dcaf_dpro*prox**2d0 &
        & + 2d0*caf*prox*2d0 &
        & + dcaf_dpro*k1ca*prox &
        & + caf*k1ca &
        & + dcaf_dpro*k1cahco3*k1*k2*kco2*pco2x*prox &
        & + caf*k1cahco3*k1*k2*kco2*pco2x &
        ! al
        & + 3d0*dalf_dpro*prox**2d0 &
        & + 3d0*alf*prox*2d0 &
        & + 2d0*dalf_dpro*k1al*prox &
        & + 2d0*alf*k1al &
        & + dalf_dpro*k2al &
        & - dalf_dpro*k4al/prox**2d0 &
        & - alf*k4al*(-2d0)/prox**3d0 &
        & + dalf_dpro*k1also4*so4f*prox**2d0 &
        & + alf*k1also4*so4f*prox*2d0 &
        ! fe2
        & + 2d0*dfe2f_dpro*prox**2d0 &
        & + 2d0*fe2f*prox*2d0 &
        & + dfe2f_dpro*k1fe2*prox &
        & + fe2f*k1fe2 &
        & + dfe2f_dpro*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        & + fe2f*k1fe2hco3*k1*k2*kco2*pco2x &
        ! fe3
        & + 3d0*dfe3f_dpro*prox**2d0 &
        & + 3d0*fe3f*prox*2d0 &
        & + 2d0*dfe3f_dpro*k1fe3*prox &
        & + 2d0*fe3f*k1fe3 &
        & + dfe3f_dpro*k2fe3 &
        & - dfe3f_dpro*k4fe3/prox**2d0 &
        & - fe3f*k4fe3*(-2d0)/prox**3d0 &
        & + dfe3f_dpro*k1fe3so4*so4f*prox**2d0 &
        & + fe3f*k1fe3so4*so4f*prox*2d0 
    
    df12 =   &
        ! so4
        & -2d0*1d0*prox**2d0 &
        & -1d0*1d0*prox**3d0*k1so4 &
        ! k
        & + dkf_dso4f*prox**2d0 &
        & - dkf_dso4f*prox**2d0*k1kso4*so4f &
        & - kf*prox**2d0*k1kso4*1d0 &
        ! na
        & + dnaf_dso4f*prox**2d0 &
        & -1d0*dnaf_dso4f*k1naco3*k1*k2*kco2*pco2x  &
        & -1d0*dnaf_dso4f*prox**2d0*k1naso4*so4f  &
        & -1d0*naf*prox**2d0*k1naso4*1d0  &
        ! si 
        & - dsif_dso4f*k1si*prox &
        & - 2d0*dsif_dso4f*k2si &
        ! mg
        & + 2d0*dmgf_dso4f*prox**2d0 &
        & + dmgf_dso4f*k1mg*prox &
        & + dmgf_dso4f*k1mghco3*k1*k2*kco2*pco2x*prox  &
        ! ca
        & + 2d0*dcaf_dso4f*prox**2d0 &
        & + dcaf_dso4f*k1ca*prox &
        & + dcaf_dso4f*k1cahco3*k1*k2*kco2*pco2x*prox &
        ! al
        & + 3d0*dalf_dso4f*prox**2d0 &
        & + 2d0*dalf_dso4f*k1al*prox &
        & + dalf_dso4f*k2al &
        & - dalf_dso4f*k4al/prox**2d0 &
        & + dalf_dso4f*k1also4*so4f*prox**2d0 &
        & + alf*k1also4*1d0*prox**2d0 &
        ! fe2
        & + 2d0*dfe2f_dso4f*prox**2d0 &
        & + dfe2f_dso4f*k1fe2*prox &
        & + dfe2f_dso4f*k1fe2hco3*k1*k2*kco2*pco2x*prox &
        ! fe3
        & + 3d0*dfe3f_dso4f*prox**2d0 &
        & + 2d0*dfe3f_dso4f*k1fe3*prox &
        & + dfe3f_dso4f*k2fe3 &
        & - dfe3f_dso4f*k4fe3/prox**2d0 &
        & + dfe3f_dso4f*k1fe3so4*so4f*prox**2d0 &
        & + fe3f*k1fe3so4*1d0*prox**2d0 
        
        
        
    
    f2 = prox**2d0*so4x - prox**2d0*so4f*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )
        
    df2 =  - prox**2d0*1d0*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & ) &
        & - prox**2d0*so4f*(  &
        & +k1kso4*dkf_dso4f &
        & +k1naso4*dnaf_dso4f &
        & +k1caso4*dcaf_dso4f &
        & +k1mgso4*dmgf_dso4f &
        & +k1fe2so4*dfe2f_dso4f &
        & +k1also4*dalf_dso4f &
        & +k1fe3so4*dfe3f_dso4f &
        & )
        
    df21 = 2d0*prox*so4x- prox**2d0*so4f*( k1so4*1d0 &
        & +k1kso4*dkf_dpro &
        & +k1naso4*dnaf_dpro &
        & +k1caso4*dcaf_dpro &
        & +k1mgso4*dmgf_dpro &
        & +k1fe2so4*dfe2f_dpro &
        & +k1also4*dalf_dpro &
        & +k1fe3so4*dfe3f_dpro &
        & ) &
        & -2d0*prox*so4f*( 1d0+k1so4*prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )
        
        
        
    df1 = df1*prox
    df21 = df21*prox
    df2 = df2*so4f
    df12 = df12*so4f
    
    if (any(isnan(f1)).or.any(isnan(f2)).or.any(isnan(df1)).or.any(isnan(df2)) &
        & .or.any(isnan(df12)).or.any(isnan(df21))) then 
        print *,any(isnan(f1)),any(isnan(f2)),any(isnan(df1)),any(isnan(df2)) &
            & ,any(isnan(df12)),any(isnan(df21))
        print *,prox
        exit
        ! pause 
    endif 
    
    
    if (nmx/=nz) then 
        amx = 0d0
        ymx = 0d0
        
        ymx(1:nz) = f1
        ymx(nz+1:nmx) = f2
        
        do iz=1,nz
            amx(iz,iz)=df1(iz)
            amx(nz+iz,nz+iz)=df2(iz)
            amx(iz,nz+iz)=df12(iz)
            amx(nz+iz,iz)=df21(iz)
        enddo 
        ymx = -ymx
        
        call DGESV(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 
        
        prox = prox*exp( ymx(1:nz) )
        so4f = so4f*exp( ymx(nz+1:nmx) )
        
        error = maxval(abs(exp( ymx )-1d0))
        if (isnan(error)) then 
            error = 1d4
            ph_error = .true.
            exit 
        endif 
    else 
        prox = prox*exp( -f1/df1 )
        error = maxval(abs(exp( -f1/df1 )-1d0))
        if (isnan(error)) error = 1d4
    endif 
    
    iter = iter + 1
    
    if (iter > 3000) then 
        print *,'iteration exceeds 3000'
        ph_error = .true.
        return
    endif 
enddo 
if (any(isnan(prox))) then     
    print *, (-log10(prox(iz)),iz=1,nz,nz/5)
    print*,'ph is nan'
    ph_error = .true.
    ! stop
endif 

    
kf = kx/(1d0+k1kso4*so4f)

naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
    
caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
    
mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
    
fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
    
sif = six/(1d0+k1si/prox+k2si/prox**2d0)

alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)

fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)

if (print_cb) then 
    open(88,file = trim(adjustl(print_loc)),status='replace')
    write(88,*) ' z ',' h+ ',' oh- ',' na+ ',' naco3- ', ' nahco3 ', ' naso4- ', ' k+ ',' kso4- ' &
        & ,' so42- ', ' hso4- ', 'hco3- ', ' co32- ' &
        & ,' h4sio4 ',' h3sio4- ',' h2sio42- ' &
        & ,' mg2+ ', ' mg(oh)+ ', ' mgco3 ', 'mghco3+ ', 'mgso4 ' &
        & , ' ca2+ ', ' ca(oh)+ ', ' caco3 ', ' cahco3+ ', ' caso4 ' &
        & ,' al3+ ', ' al(oh)2+ ', ' al(oh)2+ ', ' al(oh)3 ', ' al(oh)4- ' , ' also4+ ' &
        & , ' fe22+ ', ' fe2(oh)+ ', ' fe2co3 ', ' fe2hco3+ ', ' fe2so4 ' &
        & ,' fe3+ ', ' fe3(oh)2+ ', ' fe3(oh)2+ ', ' fe3(oh)3 ', ' fe3(oh)4- ', ' fe3so4+ ' &
        & , ' total_charge ' 
    do iz=1,nz
        write(88,*) z(iz) &
        & ,prox(iz) &
        & ,kw/prox(iz) &
        & ,naf(iz) &
        & ,naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)     &
        & ,naf(iz)*k1naso4*so4f(iz)     &
        & ,kf(iz) &
        & ,kf(iz)*k1kso4*so4f(iz) &
        & ,so4f(iz) &
        & ,so4f(iz)*k1so4*prox(iz) &
        & ,k1*kco2*pco2x(iz)/prox(iz) &
        & ,k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,sif(iz) &
        & ,sif(iz)*k1si/prox(iz) &
        & ,sif(iz)*k2si/prox(iz)**2d0 &
        & ,mgf(iz) &
        & ,mgf(iz)*k1mg/prox(iz)  &
        & ,mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & ,mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & ,mgf(iz)*k1mgso4*so4f(iz)  &
        & ,caf(iz) &
        & ,caf(iz)*k1ca/prox(iz) &
        & ,caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,caf(iz)*k1caso4*so4f(iz) &
        & ,alf(iz) &
        & ,alf(iz)*k1al/prox(iz) &
        & ,alf(iz)*k2al/prox(iz)**2d0 &
        & ,alf(iz)*k3al/prox(iz)**3d0 &
        & ,alf(iz)*k4al/prox(iz)**4d0 &
        & ,alf(iz)*k1also4*so4f(iz) &
        & ,fe2f(iz) &
        & ,fe2f(iz)*k1fe2/prox(iz) &
        & ,fe2f(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & ,fe2f(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & ,fe2f(iz)*k1fe2so4*so4f(iz) &
        & ,fe3f(iz) &
        & ,fe3f(iz)*k1fe3/prox(iz) &
        & ,fe3f(iz)*k2fe3/prox(iz)**2d0 &
        & ,fe3f(iz)*k3fe3/prox(iz)**3d0 &
        & ,fe3f(iz)*k4fe3/prox(iz)**4d0 &
        & ,fe3f(iz)*k1fe3so4*so4f(iz) &
        ! charge balance 
        & ,1d0*prox(iz) &
        & +(-1d0)*kw/prox(iz) &
        & +(1d0)*naf(iz) &
        & +(-1d0)*naf(iz)*k1naco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*naf(iz)*k1nahco3*k1*k2*kco2*pco2x(iz)/prox(iz)      &
        & +(-1d0)*naf(iz)*k1naso4*so4f(iz)     &
        & +(1d0)*kf(iz) &
        & +(-1d0)*kf(iz)*k1kso4*so4f(iz) &
        & +(-2d0)*so4f(iz) &
        & +(-1d0)*so4f(iz)*k1so4*prox(iz) &
        & +(-1d0)*k1*kco2*pco2x(iz)/prox(iz) &
        & +(-2d0)*k2*k1*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(0d0)*sif(iz) &
        & +(-1d0)*sif(iz)*k1si/prox(iz) &
        & +(-2d0)*sif(iz)*k2si/prox(iz)**2d0 &
        & +(2d0)*mgf(iz) &
        & +(1d0)*mgf(iz)*k1mg/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0  &
        & +(1d0)*mgf(iz)*k1mghco3*k1*k2*kco2*pco2x(iz)/prox(iz)  &
        & +(0d0)*mgf(iz)*k1mgso4*so4f(iz)  &
        & +(2d0)*caf(iz) &
        & +(1d0)*caf(iz)*k1ca/prox(iz) &
        & +(0d0)*caf(iz)*k1caco3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*caf(iz)*k1cahco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*caf(iz)*k1caso4*so4f(iz) &
        & +(3d0)*alf(iz) &
        & +(2d0)*alf(iz)*k1al/prox(iz) &
        & +(1d0)*alf(iz)*k2al/prox(iz)**2d0 &
        & +(0d0)*alf(iz)*k3al/prox(iz)**3d0 &
        & +(-1d0)*alf(iz)*k4al/prox(iz)**4d0 &
        & +(1d0)*alf(iz)*k1also4*so4f(iz) &
        & +(2d0)*fe2f(iz) &
        & +(1d0)*fe2f(iz)*k1fe2/prox(iz) &
        & +(0d0)*fe2f(iz)*k1fe2co3*k1*k2*kco2*pco2x(iz)/prox(iz)**2d0 &
        & +(1d0)*fe2f(iz)*k1fe2hco3*k1*k2*kco2*pco2x(iz)/prox(iz) &
        & +(0d0)*fe2f(iz)*k1fe2so4*so4f(iz) &
        & +(3d0)*fe3f(iz) &
        & +(2d0)*fe3f(iz)*k1fe3/prox(iz) &
        & +(1d0)*fe3f(iz)*k2fe3/prox(iz)**2d0 &
        & +(0d0)*fe3f(iz)*k3fe3/prox(iz)**3d0 &
        & +(-1d0)*fe3f(iz)*k4fe3/prox(iz)**4d0 &
        & +(1d0)*fe3f(iz)*k1fe3so4*so4f(iz) 
    enddo 
    close(88)
endif 
            

endsubroutine calc_pH_v7_2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_so4( &
    & nz,so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox &! input 
    & ,nsp_gas_all,nsp_aq_all,chraq_all,chrgas_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s &! input 
    & ,so4f,so4_error &! output
    & )
implicit none 

integer,intent(in)::nz,nsp_gas_all,nsp_aq_all
real(kind=8),dimension(nz),intent(in)::so4x,nax,kx,mgx,cax,fe2x,alx,fe3x,pco2x,prox
real(kind=8),dimension(nz),intent(out)::so4f
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
character(5),dimension(nsp_aq_all),intent(in)::chraq_all
character(5),dimension(nsp_gas_all),intent(in)::chrgas_all
logical,intent(out)::so4_error

real(kind=8) kco2,k1,k2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,k1al,k2al,k3al,k4al &
    & ,k1fe2,k1fe2co3,k1fe2hco3,k1fe3,k2fe3,k3fe3,k4fe3,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42

real(kind=8),dimension(nz)::f,df
real(kind=8) error,tol 
integer iter

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/


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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3  = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)


error = 1d4
tol = 1d-6

iter = 0

so4f = so4x

do while (error > tol)
    f = so4x - so4f*( 1d0+k1so4/prox &
        & +k1kso4*kx/(1d0+k1kso4*so4f) &
        & +k1naso4*nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f) &
        & +k1caso4*cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f) &
        & +k1mgso4*mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f) &
        & +k1fe2so4*fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f) &
        & +k1also4*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f) &
        & +k1fe3so4*fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f) &
        & )
        
    df = - 1d0*( 1d0+k1so4/prox &
        & +k1kso4*kx/(1d0+k1kso4*so4f) &
        & +k1naso4*nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f) &
        & +k1caso4*cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f) &
        & +k1mgso4*mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f) &
        & +k1fe2so4*fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f) &
        & +k1also4*alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f) &
        & +k1fe3so4*fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f) &
        & ) &
        & - so4f*( &
        & +k1kso4*kx*(-1d0)/(1d0+k1kso4*so4f)**2d0*k1kso4 &
        & +k1naso4*nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0*k1naso4 &
        & +k1caso4*cax*(-1d0) &
        &       /(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0*k1caso4 &
        & +k1mgso4*mgx*(-1d0) &
        &       /(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0*k1mgso4 &
        & +k1fe2so4*fe2x*(-1d0) &
        &       /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0*k1fe2so4 &
        & +k1also4*alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0*k1also4 &
        & +k1fe3so4*fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0*k1fe3so4 &
        & )
    
    
    df = df*so4f
    so4f = so4f*exp( -f/df )
    error = maxval(abs(exp( -f/df )-1d0))
    if (any(isnan(f)).or.any(isnan(df))) then 
        ! print *,so4x
        ! print *,any(isnan(f))
        ! print *,f
        ! print *,any(isnan(df))
        ! print *,df
        so4_error = .true.
        exit
    endif 
    if (isnan(error)) error = 1d4
    
    ! print *,iter,error
    
    iter = iter + 1
enddo 

so4_error = .false.
if (any(isnan(so4f))) then     
    print *, so4f
    print*,'so4f is nan'
    print*,so4x
    so4_error = .true.
    stop
endif 

endsubroutine calc_so4

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
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all &! input
    & ,prox,mineral,sp_name &! input 
    & ,omega,domega_dmsp,omega_error &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka,keqgb,keqct,k1fe2,k1fe2co3,k1fe2hco3,keqfa,k1fe3,k2fe3,k3fe3,k4fe3,keqgt &
    & ,keqcabd,keqdp,keqhb,keqkfs,keqamsi,keqg1,keqg2,keqg3,po2th,k1naco3,k1nahco3
real(kind=8),dimension(nz),intent(in):: prox
real(kind=8),dimension(nz):: pco2x,cax,mgx,six,nax,alx,po2x,fe2x,fe3x,kx
real(kind=8),dimension(nz),intent(out):: domega_dmsp
real(kind=8),dimension(nz),intent(out):: omega
logical,intent(out)::omega_error
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
real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)


po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))


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
            & nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
            & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
            & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
            & /prox**4d0 &
            & /keqab &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox)**2d0 & 
                    & *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0) &
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    !
                    & +nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
                    & *alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0)**2d0 &
                    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    ! 
                    & +nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0*(-3d0)/(1d0+k1si/prox+k2si/prox**2d0)**4d0 &
                    & *(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0) &
                    & /prox**4d0 &
                    & /keqab &
                    ! 
                    & +nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqab &
                    & )
            case('na')
                domega_dmsp = ( & 
                    & 1d0/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
                    & *alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
                    & *1d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0) &
                    & *six**3d0/(1d0+k1si/prox+k2si/prox**2d0)**3d0 &
                    & /prox**4d0 &
                    & /keqab &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox) & 
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
        omega = 1d0 - po2x**0.5d0*merge(0d0,1d0,po2x<po2th)
        ! omega = merge(1d0,1d0 - po2x**0.5d0,po2x<po2th.or. isnan(po2x**0.5d0) .or. isnan(po2x**(-0.5d0)))
        ! omega = 0d0
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 0.5d0/po2x**(0.5d0)*merge(0d0,1d0,po2x<po2th) &
                    & )
                ! domega_dmsp = merge(0d0,-0.5d0*po2x**(-0.5d0),po2x<po2th.or. isnan(po2x**0.5d0) .or. isnan(po2x**(-0.5d0)))
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

omega_error = .false.
if (any(isnan(omega)) .or. any(isnan(domega_dmsp))) then 
    print *,'nan in calc_omega_dev',any(isnan(omega)),any(isnan(domega_dmsp)),mineral,sp_name
    omega_error = .true.
    ! stop
endif 

endsubroutine calc_omega_dev

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine calc_omega_dev_v2( &
    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
    & ,prox,mineral,sp_name &! input 
    & ,omega,domega_dmsp,omega_error &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &
    & ,k1al,k2al,k3al,k4al,keqka,keqgb,keqct,k1fe2,k1fe2co3,k1fe2hco3,keqfa,k1fe3,k2fe3,k3fe3,k4fe3,keqgt &
    & ,keqcabd,keqdp,keqhb,keqkfs,keqamsi,keqg1,keqg2,keqg3,po2th,k1naco3,k1nahco3,k1so4,k1kso4,k1naso4  &
    & ,k1caso4,k1mgso4,k1fe2so4,k1also4,k1also42,k1fe3so4,k1fe3so42,mo2g1,mo2g2,mo2g3,keq_tmp,ss_x
real(kind=8),dimension(nz),intent(in):: prox,so4f
real(kind=8),dimension(nz):: pco2x,cax,mgx,six,nax,alx,po2x,fe2x,fe3x,kx,so4x
real(kind=8),dimension(nz):: caf,mgf,sif,naf,alf,fe2f,fe3f,kf
real(kind=8),dimension(nz):: dcaf_dpro,dmgf_dpro,dsif_dpro,dnaf_dpro,dalf_dpro,dfe2f_dpro,dfe3f_dpro,dkf_dpro
real(kind=8),dimension(nz):: dcaf_dpco2,dmgf_dpco2,dsif_dpco2,dnaf_dpco2,dalf_dpco2,dfe2f_dpco2,dfe3f_dpco2,dkf_dpco2
real(kind=8),dimension(nz):: dcaf_dca,dmgf_dmg,dsif_dsi,dnaf_dna,dalf_dal,dfe2f_dfe2,dfe3f_dfe3,dkf_dk,dso4f_dso4,dso4f_dpro &
    & ,dso4f_dpco2,dkf_dso4f,dnaf_dso4f,dcaf_dso4f,dmgf_dso4f,dalf_dso4f,dfe2f_dso4f,dfe3f_dso4f,dsif_dso4f
real(kind=8),dimension(nz),intent(out):: domega_dmsp
real(kind=8),dimension(nz),intent(out):: omega
logical,intent(out)::omega_error
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
real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c,keqaq_s
real(kind=8),dimension(nsp_sld_all),intent(in)::keqsld_all
real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

! real(kind=8)::thon = 1d0
real(kind=8)::thon = -1d100

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
k1naco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_co3)
k1nahco3 = keqaq_c(findloc(chraq_all,'na',dim=1),ieqaq_hco3)

k1so4 = keqaq_s(findloc(chraq_all,'so4',dim=1),ieqaq_so4)
k1naso4 = keqaq_s(findloc(chraq_all,'na',dim=1),ieqaq_so4)
k1kso4 = keqaq_s(findloc(chraq_all,'k',dim=1),ieqaq_so4)
k1caso4 = keqaq_s(findloc(chraq_all,'ca',dim=1),ieqaq_so4)
k1mgso4 = keqaq_s(findloc(chraq_all,'mg',dim=1),ieqaq_so4)
k1fe2so4 = keqaq_s(findloc(chraq_all,'fe2',dim=1),ieqaq_so4)
k1also4 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so4)
k1fe3so4 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so4)
k1also42 = keqaq_s(findloc(chraq_all,'al',dim=1),ieqaq_so42)
k1fe3so42 = keqaq_s(findloc(chraq_all,'fe3',dim=1),ieqaq_so42)

mo2g1 = keqsld_all(findloc(chrsld_all,'g1',dim=1))
mo2g2 = keqsld_all(findloc(chrsld_all,'g2',dim=1))
mo2g3 = keqsld_all(findloc(chrsld_all,'g3',dim=1))

po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))


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
so4x = 0d0

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
if (any(chraq=='so4')) then 
    so4x = maqx(findloc(chraq,'so4',dim=1),:)
elseif (any(chraq_cnst=='so4')) then 
    so4x = maqc(findloc(chraq_cnst,'so4',dim=1),:)
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


    
kf = kx/(1d0+k1kso4*so4f)
dkf_dk = 1d0/(1d0+k1kso4*so4f)
dkf_dpro = 0d0
dkf_dpco2 = 0d0
dkf_dso4f = kx*(-1d0)/(1d0+k1kso4*so4f)**2d0 *k1kso4

naf = nax/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
dnaf_dna = 1d0/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)
dnaf_dpro = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
    & *(k1naco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1nahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
dnaf_dpco2 = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
    & *(k1naco3*k1*k2*kco2*1d0/prox**2d0+k1nahco3*k1*k2*kco2*1d0/prox)
dnaf_dso4f = nax*(-1d0)/(1d0+k1naco3*k1*k2*kco2*pco2x/prox**2d0+k1nahco3*k1*k2*kco2*pco2x/prox+k1naso4*so4f)**2d0 &
    & * k1naso4
    
caf = cax/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
dcaf_dca = 1d0/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)
dcaf_dpro = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
    & *(k1ca*(-1d0)/prox**2d0+k1caco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1cahco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
dcaf_dpco2 = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
    & *(k1caco3*k1*k2*kco2*1d0/prox**2d0+k1cahco3*k1*k2*kco2*1d0/prox)
dcaf_dso4f = cax*(-1d0)/(1d0+k1ca/prox+k1caco3*k1*k2*kco2*pco2x/prox**2d0+k1cahco3*k1*k2*kco2*pco2x/prox+k1caso4*so4f)**2d0 &
    & * k1caso4
    
mgf = mgx/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
dmgf_dmg = 1d0/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)
dmgf_dpro = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
    & *(k1mg*(-1d0)/prox**2d0+k1mgco3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1mghco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
dmgf_dpco2 = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
    & *(k1mgco3*k1*k2*kco2*1d0/prox**2d0+k1mghco3*k1*k2*kco2*1d0/prox)
dmgf_dso4f = mgx*(-1d0)/(1d0+k1mg/prox+k1mgco3*k1*k2*kco2*pco2x/prox**2d0+k1mghco3*k1*k2*kco2*pco2x/prox+k1mgso4*so4f)**2d0 &
    & * k1mgso4
    
fe2f = fe2x/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
dfe2f_dfe2 = 1d0/(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)
dfe2f_dpro = fe2x*(-1d0) &
    & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
    & *(k1fe2*(-1d0)/prox**2d0+k1fe2co3*k1*k2*kco2*pco2x*(-2d0)/prox**3d0+k1fe2hco3*k1*k2*kco2*pco2x*(-1d0)/prox**2d0)
dfe2f_dpco2 = fe2x*(-1d0) &
    & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
    & *(k1fe2co3*k1*k2*kco2*1d0/prox**2d0+k1fe2hco3*k1*k2*kco2*1d0/prox)
dfe2f_dso4f = fe2x*(-1d0) &
    & /(1d0+k1fe2/prox+k1fe2co3*k1*k2*kco2*pco2x/prox**2d0+k1fe2hco3*k1*k2*kco2*pco2x/prox+k1fe2so4*so4f)**2d0 &
    & * k1fe2so4
    
sif = six/(1d0+k1si/prox+k2si/prox**2d0)
dsif_dsi = 1d0/(1d0+k1si/prox+k2si/prox**2d0)
dsif_dpro = six*(-1d0)/(k1si*(-1d0)/prox**2d0+k2si*(-2d0)/prox**3d0)
dsif_dpco2 = 0d0
dsif_dso4f = 0d0

alf = alx/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)
dalf_dal = 1d0/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)
dalf_dpro = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 &
    & *(k1al*(-1d0)/prox**2d0+k2al*(-2d0)/prox**3d0+k3al*(-3d0)/prox**4d0+k4al*(-4d0)/prox**5d0)
dalf_dpco2 = 0d0
dalf_dso4f = alx*(-1d0)/(1d0+k1al/prox+k2al/prox**2d0+k3al/prox**3d0+k4al/prox**4d0+k1also4*so4f)**2d0 & 
    & * k1also4

fe3f = fe3x/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)
dfe3f_dfe3 = 1d0/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)
dfe3f_dpro = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
    & *(k1fe3*(-1d0)/prox**2d0+k2fe3*(-2d0)/prox**3d0+k3fe3*(-3d0)/prox**4d0+k4fe3*(-4d0)/prox**5d0)
dfe3f_dpco2 = 0d0
dfe3f_dso4f = fe3x*(-1d0)/(1d0+k1fe3/prox+k2fe3/prox**2d0+k3fe3/prox**3d0+k4fe3/prox**4d0+k1fe3so4*so4f)**2d0 &
    & * k1fe3so4

dso4f_dso4 = 1d0/( 1d0+k1so4/prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )

dso4f_dpro = so4x*(-1d0)/( 1d0+k1so4/prox &
        & +k1kso4*kf &
        & +k1naso4*naf &
        & +k1caso4*caf &
        & +k1mgso4*mgf &
        & +k1fe2so4*fe2f &
        & +k1also4*alf &
        & +k1fe3so4*fe3f &
        & )**2d0 &
        & * (k1so4*(-1d0)/prox**2d0 &
        & )
! dso4f_dpro = 0d0
dso4f_dpco2 = 0d0

select case(trim(adjustl(mineral)))
    case('fo')
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
        omega = ( & 
            & mgf**2d0 &
            & *sif &
            & /prox**4d0 &
            & /keqfo &
            & )
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & mgf*2d0*dmgf_dpro &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfo &
                    !
                    & +mgf**2d0 &
                    & *dsif_dpro &
                    & /prox**4d0 &
                    & /keqfo &
                    !
                    & +mgf**2d0 &
                    & *sif &
                    & *(-4d0)/prox**5d0 &
                    & /keqfo &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & mgf*2d0*dmgf_dso4f &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfo &
                    !
                    & +mgf**2d0 &
                    & *dsif_dso4f &
                    & /prox**4d0 &
                    & /keqfo &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & mgf*2d0*dmgf_dmg &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfo &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & mgf**2d0 &
                    & *dsif_dsi &
                    & /prox**4d0 &
                    & /keqfo &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & mgf*2d0*dmgf_dpco2 &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfo &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('fa')
    ! Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O 
        omega = ( & 
            & fe2f**2d0 &
            & *sif &
            & /prox**4d0 &
            & /keqfa &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & fe2f*2d0*dfe2f_dpro &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfa &
                    !
                    & +fe2f**2d0 &
                    & *dsif_dpro &
                    & /prox**4d0 &
                    & /keqfa &
                    !
                    & +fe2f**2d0 &
                    & *sif &
                    & *(-4d0)/prox**5d0 &
                    & /keqfa &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & fe2f*2d0*dfe2f_dso4f &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfa &
                    !
                    & +fe2f**2d0 &
                    & *dsif_dso4f &
                    & /prox**4d0 &
                    & /keqfa &
                    & )
            case('fe2')
                domega_dmsp = ( & 
                    & fe2f*2d0*dfe2f_dfe2 &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfa &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & fe2f**2d0 &
                    & *1d0*dsif_dsi &
                    & /prox**4d0 &
                    & /keqfa &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & fe2f*2d0*dfe2f_dpco2 &
                    & *sif &
                    & /prox**4d0 &
                    & /keqfa &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    ! case('ab')
    ! NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
        ! omega = ( & 
            ! & naf & 
            ! & *alf &
            ! & *sif**3d0 &
            ! & /prox**4d0 &
            ! & /keqab &
            ! & )
            
        ! select case(trim(adjustl(sp_name)))
            ! case('pro')
                ! domega_dmsp = ( & 
                    ! & 1d0*dnaf_dpro & 
                    ! & *alf &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    
                    ! & +naf & 
                    ! & *1d0*dalf_dpro &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    
                    ! & +naf & 
                    ! & *alf &
                    ! & *3d0*sif**2d0*dsif_dpro &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    
                    ! & +naf & 
                    ! & *alf &
                    ! & *sif**3d0 &
                    ! & *(-4d0)/prox**5d0 &
                    ! & /keqab &
                    ! & )
            ! case('so4f')
                ! domega_dmsp = ( & 
                    ! & 1d0*dnaf_dso4f & 
                    ! & *alf &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    
                    ! & +naf & 
                    ! & *1d0*dalf_dso4f &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    
                    ! & +naf & 
                    ! & *alf &
                    ! & *3d0*sif**2d0*dsif_dso4f &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    ! & )
            ! case('na')
                ! domega_dmsp = ( & 
                    ! & 1d0*dnaf_dna & 
                    ! & *alf &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    ! & )
            ! case('al')
                ! domega_dmsp = ( & 
                    ! & naf & 
                    ! & *1d0*dalf_dal &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    ! & )
            ! case('si')
                ! domega_dmsp = ( & 
                    ! & naf & 
                    ! & *alf &
                    ! & *3d0*sif**2d0*dsif_dsi &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    ! & )
            ! case('pco2')
                ! domega_dmsp = ( & 
                    ! & 1d0*dnaf_dpco2 & 
                    ! & *alf &
                    ! & *sif**3d0 &
                    ! & /prox**4d0 &
                    ! & /keqab &
                    ! & )
            ! case default 
                ! domega_dmsp = 0d0
        ! endselect 
        
        
    case('kfs')
    ! K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
        omega = ( & 
            & kf &
            & *alf &
            & *sif**3d0 &
            & /prox**4d0 &
            & /keqkfs   &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 1d0*dkf_dpro &
                    & *alf &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    ! 
                    & +kf &
                    & *1d0*dalf_dpro &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    ! 
                    & +kf &
                    & *alf &
                    & *3d0*sif**2d0*dsif_dpro &
                    & /prox**4d0 &
                    & /keqkfs   &
                    !
                    & +kf &
                    & *alf &
                    & *sif**3d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqkfs   &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 1d0*dkf_dso4f &
                    & *alf &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    ! 
                    & +kf &
                    & *1d0*dalf_dso4f &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    ! 
                    & +kf &
                    & *alf &
                    & *3d0*sif**2d0*dsif_dso4f &
                    & /prox**4d0 &
                    & /keqkfs   &
                    & )
            case('k')
                domega_dmsp = ( & 
                    & 1d0*dkf_dk &
                    & *alf &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & kf &
                    & *1d0*dalf_dal &
                    & *sif**3d0 &
                    & /prox**4d0 &
                    & /keqkfs   &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & kf &
                    & *alf &
                    & *3d0*sif**2d0*dsif_dsi &
                    & /prox**4d0 &
                    & /keqkfs   &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    ! case('an')
    ! CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
        ! omega = ( & 
            ! & caf &
            ! & *alf**2d0 &
            ! & *sif**2d0 &
            ! & /prox**8d0 &
            ! & /keqan &
            ! & )
            
        ! select case(trim(adjustl(sp_name)))
            ! case('pro')
                ! domega_dmsp = ( & 
                    ! & 1d0*dcaf_dpro &
                    ! & *alf**2d0 &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & + &
                    ! & caf &
                    ! & *alf*2d0*dalf_dpro &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & + &
                    ! & caf &
                    ! & *alf**2d0 &
                    ! & *sif*2d0*dsif_dpro &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & + &
                    ! & caf &
                    ! & *alf**2d0 &
                    ! & *sif**2d0 &
                    ! & *(-8d0)/prox**9d0 &
                    ! & /keqan &
                    ! & )
            ! case('so4f')
                ! domega_dmsp = ( & 
                    ! & 1d0*dcaf_dso4f &
                    ! & *alf**2d0 &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & + &
                    ! & caf &
                    ! & *alf*2d0*dalf_dso4f &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & + &
                    ! & caf &
                    ! & *alf**2d0 &
                    ! & *sif*2d0*dsif_dso4f &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & )
            ! case('pco2')
                ! domega_dmsp = ( & 
                    ! & 1d0*dcaf_dpco2 &
                    ! & *alf**2d0 &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & )
            ! case('ca')
                ! domega_dmsp = ( & 
                    ! & 1d0*dcaf_dca &
                    ! & *alf**2d0 &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & )
            ! case('al')
                ! domega_dmsp = ( & 
                    ! & caf &
                    ! & *alf*2d0*dalf_dal &
                    ! & *sif**2d0 &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & )
            ! case('si')
                ! domega_dmsp = ( & 
                    ! & caf &
                    ! & *alf**2d0 &
                    ! & *sif*2d0*dsif_dsi &
                    ! & /prox**8d0 &
                    ! & /keqan &
                    ! & )
            ! case default 
                ! domega_dmsp = 0d0
        ! endselect 
        
        
    case('la','ab','an')
        ! CaxNa(1-x)Al(1+x)Si(3-x)O8 + (4x + 4) = xCa+2 + (1-x)Na+ + (1+x)Al+++ + (3-x)SiO2(aq) 
        keq_tmp = keqsld_all(findloc(chrsld_all,mineral,dim=1))
        ss_x = staq_all(findloc(chrsld_all,mineral,dim=1), findloc(chraq_all,'ca',dim=1) )
        omega = ( & 
            & caf**ss_x &
            & *naf**(1d0-ss_x) &
            & *alf**(1d0+ss_x) &
            & *sif**(3d0-ss_x) &
            & /prox**(4d0 + 4d0*ss_x) &
            & /keq_tmp &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & ss_x*caf**(ss_x-1d0)*dcaf_dpro &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *(1d0-ss_x)*naf**(1d0-ss_x-1d0)*dnaf_dpro &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *(1d0+ss_x)*alf**(1d0+ss_x-1d0)*dalf_dpro &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *(3d0-ss_x)*sif**(3d0-ss_x-1d0)*dsif_dpro &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & *(-1d0*(4d0 + 4d0*ss_x))/prox**(4d0 + 4d0*ss_x+1d0) &
                    & /keq_tmp &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & ss_x*caf**(ss_x-1d0)*dcaf_dso4f &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *(1d0-ss_x)*naf**(1d0-ss_x-1d0)*dnaf_dso4f &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *(1d0+ss_x)*alf**(1d0+ss_x-1d0)*dalf_dso4f &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *(3d0-ss_x)*sif**(3d0-ss_x-1d0)*dsif_dso4f &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & ss_x*caf**(ss_x-1d0)*dcaf_dpco2 &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *(1d0-ss_x)*naf**(1d0-ss_x-1d0)*dnaf_dpco2 &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *(1d0+ss_x)*alf**(1d0+ss_x-1d0)*dalf_dpco2 &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & + &
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *(3d0-ss_x)*sif**(3d0-ss_x-1d0)*dsif_dpco2 &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & ss_x*caf**(ss_x-1d0)*dcaf_dca &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case('na')
                domega_dmsp = ( & 
                    & caf**ss_x &
                    & *(1d0-ss_x)*naf**(1d0-ss_x-1d0)*dnaf_dna &
                    & *alf**(1d0+ss_x) &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *(1d0+ss_x)*alf**(1d0+ss_x-1d0)*dalf_dal &
                    & *sif**(3d0-ss_x) &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & caf**ss_x &
                    & *naf**(1d0-ss_x) &
                    & *alf**(1d0+ss_x) &
                    & *(3d0-ss_x)*sif**(3d0-ss_x-1d0)*dsif_dsi &
                    & /prox**(4d0 + 4d0*ss_x) &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('cc')
        omega = ( &
            & caf &
            & *k1*k2*kco2*pco2x/(prox**2d0) &
            & /keqcc &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dpro &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    & + &
                    & caf &
                    & *k1*k2*kco2*pco2x*(-2d0)/(prox**3d0) &
                    & /keqcc &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dso4f &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dpco2 &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    & + &
                    & caf &
                    & *k1*k2*kco2*1d0/(prox**2d0) &
                    & /keqcc &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dca &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keqcc &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('arg')
        keq_tmp = keqsld_all(findloc(chrsld_all,'arg',dim=1))
        omega = ( &
            & caf &
            & *k1*k2*kco2*pco2x/(prox**2d0) &
            & /keq_tmp &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dpro &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keq_tmp &
                    & + &
                    & caf &
                    & *k1*k2*kco2*pco2x*(-2d0)/(prox**3d0) &
                    & /keq_tmp &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dso4f &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keq_tmp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dpco2 &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keq_tmp &
                    & + &
                    & caf &
                    & *k1*k2*kco2*1d0/(prox**2d0) &
                    & /keq_tmp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 1d0*dcaf_dca &
                    & *k1*k2*kco2*pco2x/(prox**2d0) &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('dlm')
        keq_tmp = keqsld_all(findloc(chrsld_all,'dlm',dim=1))
        omega = ( &
            & caf*mgf &
            & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
            & /keq_tmp &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & dcaf_dpro*mgf &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & + &
                    & caf*dmgf_dpro &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & + &
                    & caf*mgf &
                    & *2d0*(k1*k2*kco2*pco2x/(prox**2d0))*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0) &
                    & /keq_tmp &
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & dcaf_dso4f*mgf &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & + &
                    & caf*dmgf_dso4f &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dcaf_dpco2*mgf &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & + &
                    & caf*dmgf_dpco2 &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & + &
                    & caf*mgf &
                    & *2d0*(k1*k2*kco2*pco2x/(prox**2d0))*k1*k2*kco2*1d0/(prox**2d0) &
                    & /keq_tmp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & dcaf_dca*mgf &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & caf*dmgf_dmg &
                    & *(k1*k2*kco2*pco2x/(prox**2d0))**2d0 &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('ka')
    ! Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
        omega = ( &
            & alf**2d0 &
            & *sif**2d0 &
            & /prox**6d0 &
            & /keqka &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependences of al and si are from an case
                domega_dmsp = ( & 
                    & alf*2d0*dalf_dpro &
                    & *sif**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    & + &
                    & alf**2d0 &
                    & *sif*2d0*dsif_dpro &
                    & /prox**6d0 &
                    & /keqka &
                    & + &
                    & alf**2d0 &
                    & *sif**2d0 &
                    & *(-6d0)/prox**7d0 &
                    & /keqka &
                    ! 
                    & )
            case('so4f')
                ! ph dependences of al and si are from an case
                domega_dmsp = ( & 
                    & alf*2d0*dalf_dso4f &
                    & *sif**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    & + &
                    & alf**2d0 &
                    & *sif*2d0*dsif_dso4f &
                    & /prox**6d0 &
                    & /keqka &
                    ! 
                    & )
            case('al')
                domega_dmsp = ( & 
                    & alf*2d0*dalf_dal &
                    & *sif**2d0 &
                    & /prox**6d0 &
                    & /keqka &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & alf**2d0 &
                    & *sif*2d0*dsif_dsi &
                    & /prox**6d0 &
                    & /keqka &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('gb')
    ! Al(OH)3 + 3 H+ = Al+3 + 3 H2O 
        omega = ( &
            & alf &
            & /prox**3d0 &
            & /keqgb &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependence of al is from an case
                domega_dmsp = ( & 
                    & 1d0*dalf_dpro &
                    & /prox**3d0 &
                    & /keqgb &
                    & + &
                    & alf &
                    & *(-3d0)/prox**4d0 &
                    & /keqgb &
                    ! 
                    & )
            case('so4f')
                ! ph dependence of al is from an case
                domega_dmsp = ( & 
                    & 1d0*dalf_dso4f &
                    & /prox**3d0 &
                    & /keqgb &
                    ! 
                    & )
            case('al')
                domega_dmsp = ( & 
                    & 1d0*dalf_dal &
                    & /prox**3d0 &
                    & /keqgb &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('amsi')
    ! SiO2 + 2 H2O = H4SiO4
        omega = ( &
            & sif &
            & /keqamsi &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 1d0*dsif_dpro &
                    & /keqamsi &
                    ! 
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 1d0*dsif_dso4f &
                    & /keqamsi &
                    ! 
                    & )
            case('si')
                domega_dmsp = ( & 
                    & 1d0*dsif_dsi &
                    & /keqamsi &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
    case('qtz')
    !  SiO2 + 2H2O = H4SiO4
        keq_tmp = keqsld_all(findloc(chrsld_all,'qtz',dim=1))
        omega = ( &
            & sif &
            & /keq_tmp &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 1d0*dsif_dpro &
                    & /keq_tmp &
                    ! 
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 1d0*dsif_dso4f &
                    & /keq_tmp &
                    ! 
                    & )
            case('si')
                domega_dmsp = ( & 
                    & 1d0*dsif_dsi &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('gt')
    !  Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
        omega = (&
            & fe3f &
            & /prox**3d0 &
            & /keqgt &
            & )
            
        ! copied and pasted from gb case with replacing al and gb with fe3 and gt respectively
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & dfe3f_dpro &
                    & /prox**3d0 &
                    & /keqgt &
                    & + &
                    & fe3f &
                    & *(-3d0)/prox**4d0 &
                    & /keqgt &
                    ! 
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & dfe3f_dso4f &
                    & /prox**3d0 &
                    & /keqgt &
                    ! 
                    & )
            case('fe3')
                domega_dmsp = ( & 
                    & dfe3f_dfe3 &
                    & /prox**3d0 &
                    & /keqgt &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('hm')
    !  Fe2O3 + 6H+ = 2Fe+3 + 3H2O
        keq_tmp = keqsld_all(findloc(chrsld_all,'hm',dim=1))
        omega = (&
            & fe3f**2d0 &
            & /prox**6d0 &
            & /keq_tmp &
            & )
            
        ! copied and pasted from gb case with replacing al and gb with fe3 and gt respectively
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( & 
                    & 2d0*fe3f*dfe3f_dpro &
                    & /prox**6d0 &
                    & /keq_tmp &
                    & + &
                    & fe3f**2d0 &
                    & *(-6d0)/prox**7d0 &
                    & /keq_tmp &
                    ! 
                    & )
            case('so4f')
                domega_dmsp = ( & 
                    & 2d0*fe3f*dfe3f_dso4f &
                    & /prox**6d0 &
                    & /keq_tmp &
                    ! 
                    & )
            case('fe3')
                domega_dmsp = ( & 
                    & 2d0*fe3f*dfe3f_dfe3 &
                    & /prox**6d0 &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('ct')
    ! Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
        omega = ( &
            & mgf**3d0 &
            & *sif**2d0  &
            & /prox**6d0 &
            & /keqct &
            & )
            
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                ! ph dependence of mg and si are from above
                domega_dmsp = ( & 
                    & 3d0*mgf**2d0*dmgf_dpro &
                    & *sif**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & + &
                    & mgf**3d0 &
                    & *sif*2d0*dsif_dpro  &
                    & /prox**6d0 &
                    & /keqct &
                    & + &
                    & mgf**3d0 &
                    & *sif**2d0  &
                    & *(-6d0)/prox**7d0 &
                    & /keqct &
                    & )
            case('so4f')
                ! ph dependence of mg and si are from above
                domega_dmsp = ( & 
                    & 3d0*mgf**2d0*dmgf_dso4f &
                    & *sif**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & + &
                    & mgf**3d0 &
                    & *sif*2d0*dsif_dso4f  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case('pco2')
                ! pco2 dependence of mg is from fo case
                domega_dmsp = ( & 
                    & 3d0*mgf**2d0*dmgf_dpco2 &
                    & *sif**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & 3d0*mgf**2d0*dmgf_dmg &
                    & *sif**2d0  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & mgf**3d0 &
                    & *sif*2d0*dsif_dsi  &
                    & /prox**6d0 &
                    & /keqct &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('cabd')
    ! Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
        omega = ( &
            & caf**(1d0/6d0) &
            & *alf**(7d0/3d0) &
            & *sif**(11d0/3d0) &
            & /prox**(22d0/3d0) &
            & /keqcabd &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &
                    & (1d0/6d0)*caf**(1d0/6d0-1d0)*dcaf_dpro &
                    & *alf**(7d0/3d0) &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & + &
                    & caf**(1d0/6d0) &
                    & *(7d0/3d0)*alf**(7d0/3d0-1d0)*dalf_dpro &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & + &
                    & caf**(1d0/6d0) &
                    & *alf**(7d0/3d0) &
                    & *(11d0/3d0)*sif**(11d0/3d0-1d0)*dsif_dpro &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & + &
                    & caf**(1d0/6d0) &
                    & *alf**(7d0/3d0) &
                    & *sif**(11d0/3d0) &
                    & *(-22d0/3d0)/prox**(22d0/3d0+1d0) &
                    & /keqcabd &
                    & )
            case('so4f')
                domega_dmsp = ( &
                    & (1d0/6d0)*caf**(1d0/6d0-1d0)*dcaf_dso4f &
                    & *alf**(7d0/3d0) &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & + &
                    & caf**(1d0/6d0) &
                    & *(7d0/3d0)*alf**(7d0/3d0-1d0)*dalf_dso4f &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & + &
                    & caf**(1d0/6d0) &
                    & *alf**(7d0/3d0) &
                    & *(11d0/3d0)*sif**(11d0/3d0-1d0)*dsif_dso4f &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & (1d0/6d0)*caf**(1d0/6d0-1d0)*dcaf_dpco2 &
                    & *alf**(7d0/3d0) &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & (1d0/6d0)*caf**(1d0/6d0-1d0)*dcaf_dca &
                    & *alf**(7d0/3d0) &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & caf**(1d0/6d0) &
                    & *(7d0/3d0)*alf**(7d0/3d0-1d0)*dalf_dal &
                    & *sif**(11d0/3d0) &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & caf**(1d0/6d0) &
                    & *alf**(7d0/3d0) &
                    & *(11d0/3d0)*sif**(11d0/3d0-1d0)*dsif_dsi &
                    & /prox**(22d0/3d0) &
                    & /keqcabd &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('ill')
    ! Illite  + 8 H+  = 5 H2O  + .6 K+  + .25 Mg++  + 2.3 Al+++  + 3.5 SiO2(aq)
        keq_tmp = keqsld_all(findloc(chrsld_all,'ill',dim=1))
        omega = (  &
            & kf**(0.6d0) &
            & *mgf**(0.25d0) &
            & *alf**(2.3d0)  &
            & *sif**(3.5d0)   &
            & /prox**(8d0)    &
            & /keq_tmp      &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &
                    & 0.6d0*kf**(0.6d0-1d0)*dkf_dpro &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *0.25d0*mgf**(0.25d0-1d0)*dmgf_dpro &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *2.3d0*alf**(2.3d0-1d0)*dalf_dpro  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *3.5d0*sif**(3.5d0-1d0)*dsif_dpro   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & *(-8d0)/prox**(9d0)    &
                    & /keq_tmp      &
                    & )
            case('so4f')
                domega_dmsp = ( &
                    & 0.6d0*kf**(0.6d0-1d0)*dkf_dso4f &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *0.25d0*mgf**(0.25d0-1d0)*dmgf_dso4f &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *2.3d0*alf**(2.3d0-1d0)*dalf_dso4f  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *3.5d0*sif**(3.5d0-1d0)*dsif_dso4f   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & 0.6d0*kf**(0.6d0-1d0)*dkf_dpco2 &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *0.25d0*mgf**(0.25d0-1d0)*dmgf_dpco2 &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *2.3d0*alf**(2.3d0-1d0)*dalf_dpco2  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & + &
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *3.5d0*sif**(3.5d0-1d0)*dsif_dpco2   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case('k')
                domega_dmsp = ( & 
                    & 0.6d0*kf**(0.6d0-1d0)*dkf_dk &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & kf**(0.6d0) &
                    & *0.25d0*mgf**(0.25d0-1d0)*dmgf_dmg &
                    & *alf**(2.3d0)  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *2.3d0*alf**(2.3d0-1d0)*dalf_dal  &
                    & *sif**(3.5d0)   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & kf**(0.6d0) &
                    & *mgf**(0.25d0) &
                    & *alf**(2.3d0)  &
                    & *3.5d0*sif**(3.5d0-1d0)*dsif_dsi   &
                    & /prox**(8d0)    &
                    & /keq_tmp      &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('anl')
    ! NaAlSi2O6*H2O  + 5 H2O  = Na+  + Al(OH)4-  + 2 Si(OH)4(aq)
        keq_tmp = keqsld_all(findloc(chrsld_all,'anl',dim=1))
        omega = (  &
            & naf &
            & *alf*k4al/prox**4d0  &
            & *sif**2d0   &
            & /keq_tmp      &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &
                    & dnaf_dpro &
                    & *alf*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dpro*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf*k4al*(-4d0)/prox**5d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf*k4al/prox**4d0  &
                    & *2d0*sif*dsif_dpro   &
                    & /keq_tmp      &
                    & )
            case('so4f')
                domega_dmsp = ( &
                    & dnaf_dso4f &
                    & *alf*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dso4f*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf*k4al/prox**4d0  &
                    & *2d0*sif*dsif_dso4f   &
                    & /keq_tmp      &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dnaf_dpco2 &
                    & *alf*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dpco2*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf*k4al/prox**4d0  &
                    & *2d0*sif*dsif_dpco2   &
                    & /keq_tmp      &
                    & )
            case('na')
                domega_dmsp = ( & 
                    & dnaf_dna &
                    & *alf*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & naf &
                    & *dalf_dal*k4al/prox**4d0  &
                    & *sif**2d0   &
                    & /keq_tmp      &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & naf &
                    & *alf*k4al/prox**4d0  &
                    & *2d0*sif*dsif_dsi   &
                    & /keq_tmp      &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('nph')
    ! Nepheline  + 4 H+  = 2 H2O  + SiO2(aq)  + Al+++  + Na+
        keq_tmp = keqsld_all(findloc(chrsld_all,'nph',dim=1))
        omega = (  &
            & naf &
            & *alf  &
            & *sif   &
            & /prox**4d0  &
            & /keq_tmp      &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &
                    & dnaf_dpro &
                    & *alf  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dpro  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf  &
                    & *dsif_dpro   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf  &
                    & *sif   &
                    & *(-4d0)/prox**5d0 &
                    & /keq_tmp      &
                    & )
            case('so4f')
                domega_dmsp = ( &
                    & dnaf_dso4f &
                    & *alf  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dso4f  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf  &
                    & *dsif_dso4f   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dnaf_dpco2 &
                    & *alf  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *dalf_dpco2  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & + &
                    & naf &
                    & *alf  &
                    & *dsif_dpco2   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & )
            case('na')
                domega_dmsp = ( & 
                    & dnaf_dna &
                    & *alf  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & )
            case('al')
                domega_dmsp = ( & 
                    & naf &
                    & *dalf_dal  &
                    & *sif   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & naf &
                    & *alf  &
                    & *dsif_dsi   &
                    & /prox**4d0 &
                    & /keq_tmp      &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('dp')
    ! Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
        omega = ( &
            & caf &
            & *mgf &
            & *sif**2d0 &
            & /prox**4d0 &
            & /keqdp &
            & )
            
        ! dependences from an case
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & dcaf_dpro &
                    & *mgf &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *dmgf_dpro &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *mgf &
                    & *sif*2d0*dsif_dpro &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *mgf &
                    & *sif**2d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqdp &
                    & )
            case('so4f')
                domega_dmsp = ( &   
                    & dcaf_dso4f &
                    & *mgf &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *dmgf_dso4f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *mgf &
                    & *sif*2d0*dsif_dso4f &
                    & /prox**4d0 &
                    & /keqdp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dcaf_dpco2 &
                    & *mgf &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & + &
                    & caf &
                    & *dmgf_dpco2 &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & caf &
                    & *dmgf_dmg &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & dcaf_dca &
                    & *mgf &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqdp &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & caf &
                    & *mgf &
                    & *sif*2d0*dsif_dsi &
                    & /prox**4d0 &
                    & /keqdp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('hb')
    ! Hedenbergite  + 4 H+  = 2 H2O  + 2 SiO2(aq)  + Fe++  + Ca++
        omega = ( &
            & caf &
            & *fe2f &
            & *sif**2d0 &
            & /prox**4d0 &
            & /keqhb &
            & )
            
        ! copied and pasted from dp case with replacing mg and dp by fe2 and hb
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & dcaf_dpro &
                    & *fe2f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *dfe2f_dpro &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *fe2f &
                    & *sif*2d0*dsif_dpro &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *fe2f &
                    & *sif**2d0 &
                    & *(-4d0)/prox**5d0 &
                    & /keqhb &
                    & )
            case('so4f')
                domega_dmsp = ( &   
                    & dcaf_dso4f &
                    & *fe2f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *dfe2f_dso4f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *fe2f &
                    & *sif*2d0*dsif_dso4f &
                    & /prox**4d0 &
                    & /keqhb &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dcaf_dpco2 &
                    & *fe2f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & + &
                    & caf &
                    & *dfe2f_dpco2 &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & )
            case('fe2')
                domega_dmsp = ( & 
                    & caf &
                    & *dfe2f_dfe2 &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & dcaf_dca &
                    & *fe2f &
                    & *sif**2d0 &
                    & /prox**4d0 &
                    & /keqhb &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & caf &
                    & *fe2f &
                    & *sif*2d0*dsif_dsi &
                    & /prox**4d0 &
                    & /keqhb &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('tm')
    ! Tremolite  + 14 H+  = 8 H2O  + 8 SiO2(aq)  + 2 Ca++  + 5 Mg++
        keq_tmp = keqsld_all(findloc(chrsld_all,'tm',dim=1))
        omega = ( &
            & caf**2d0 &
            & *mgf**5d0 &
            & *sif**8d0 &
            & /prox**14d0 &
            & /keq_tmp &
            & )
            
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & 2d0*caf*dcaf_dpro &
                    & *mgf**5d0 &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *5d0*mgf**4d0*dmgf_dpro &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *mgf**5d0 &
                    & *8d0*sif**7d0*dsif_dpro &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *mgf**5d0 &
                    & *sif**8d0 &
                    & *(-14d0)/prox**15d0 &
                    & /keq_tmp &
                    & )
            case('so4f')
                domega_dmsp = ( &   
                    & 2d0*caf*dcaf_dso4f &
                    & *mgf**5d0 &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *5d0*mgf**4d0*dmgf_dso4f &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *mgf**5d0 &
                    & *8d0*sif**7d0*dsif_dso4f &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & 2d0*caf*dcaf_dpco2 &
                    & *mgf**5d0 &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & + &
                    & caf**2d0 &
                    & *5d0*mgf**4d0*dmgf_dpco2 &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & )
            case('mg')
                domega_dmsp = ( & 
                    & caf**2d0 &
                    & *5d0*mgf**4d0*dmgf_dmg &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & 2d0*caf*dcaf_dca &
                    & *mgf**5d0 &
                    & *sif**8d0 &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & )
            case('si')
                domega_dmsp = ( & 
                    & caf**2d0 &
                    & *mgf**5d0 &
                    & *8d0*sif**7d0*dsif_dsi &
                    & /prox**14d0 &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('gps')
    ! CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O
        keq_tmp = keqsld_all(findloc(chrsld_all,'gps',dim=1))
        omega = ( &
            & caf &
            & *so4f &
            & /keq_tmp &
            & )
            
        ! copied and pasted from dp case with replacing mg and dp by fe2 and hb
        select case(trim(adjustl(sp_name)))
            case('pro')
                domega_dmsp = ( &   
                    & dcaf_dpro &
                    & *so4f &
                    & /keq_tmp &
                    & + &
                    & caf &
                    & *dso4f_dpro &
                    & /keq_tmp &
                    & )
            case('so4f')
                domega_dmsp = ( &   
                    & dcaf_dso4f &
                    & *so4f &
                    & /keq_tmp &
                    & + &
                    & caf &
                    & *1d0 &
                    & /keq_tmp &
                    & )
            case('pco2')
                domega_dmsp = ( & 
                    & dcaf_dpco2 &
                    & *so4f &
                    & /keq_tmp &
                    & + &
                    & caf &
                    & *dso4f_dpco2 &
                    & /keq_tmp &
                    & )
            case('ca')
                domega_dmsp = ( & 
                    & dcaf_dca &
                    & *so4f &
                    & /keq_tmp &
                    & )
            case('so4')
                domega_dmsp = ( & 
                    & caf &
                    & *dso4f_dso4 &
                    & /keq_tmp &
                    & )
            case default 
                domega_dmsp = 0d0
        endselect 
        
        
    case('py')
    ! omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x**0.5d0
    ! i.e., 1.0 - omega_py = po2x**0.5 
        ! omega = 1d0 - po2x**0.5d0
        omega = 1d0 - po2x**0.5d0*merge(0d0,1d0,po2x<po2th*thon)
        ! omega = merge(1d0,1d0 - po2x**0.5d0,po2x<po2th.or. isnan(po2x**0.5d0) .or. isnan(po2x**(-0.5d0)))
        ! omega = 0d0
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 0.5d0/po2x**(0.5d0)*merge(0d0,1d0,po2x<po2th*thon) &
                    & )
                ! domega_dmsp = merge(0d0,-0.5d0*po2x**(-0.5d0),po2x<po2th.or. isnan(po2x**0.5d0) .or. isnan(po2x**(-0.5d0)))
                ! print *,omega
                ! print *
                ! print *,domega_dmsp
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
        ! omega = 1d0 
        ! domega_dmsp = 0d0
        ! print *,mineral,sp_name
        omega = 1d0 - po2x/(po2x+mo2g1)&
            & *merge(0d0,1d0,po2x < po2th*thon)
        
        ! omega = 1d0 - po2x/(po2x + keqg1)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 1d0/(po2x + mo2g1) &
                    & - po2x*(-1d0)/(po2x + mo2g1)**2d0 &
                    & ) &
                    & *merge(0d0,1d0,po2x < po2th*thon)
            case default 
                domega_dmsp = 0d0
        endselect 
        
        ! print*, omega
        
        
    case('g2')
    ! omega is defined so that kg2*poro*hr*mvg2*1d-6*mg2x*(1d0-omega_g2) = kg2*poro*hr*mvg2*1d-6*mg2x*po2x/(po2x+mo2)
    ! i.e., 1.0 - omega_g2 = po2x/(po2x+mo2) 
        omega = 1d0 - po2x/(po2x+mo2g2)&
            & *merge(0d0,1d0,po2x < po2th*thon)
        ! print *,mineral,sp_name
        
        
        ! omega = 1d0 
        ! domega_dmsp = 0d0
        
        ! omega = 1d0 - po2x/(po2x + keqg2)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 1d0/(po2x + mo2g2) &
                    & - po2x*(-1d0)/(po2x + mo2g2)**2d0 &
                    & ) &
                    & *merge(0d0,1d0,po2x < po2th*thon)
            case default 
                domega_dmsp = 0d0
        endselect 
        
        ! print*, omega
        
        
    case('g3')
    ! omega is defined so that kg3*poro*hr*mvg3*1d-6*mg3x*(1d0-omega_g3) = kg3*poro*hr*mvg3*1d-6*mg3x*po2x/(po2x+mo2)
    ! i.e., 1.0 - omega_g3 = po2x/(po2x+mo2) 
        omega = 1d0 - po2x/(po2x+mo2g3) &
            & *merge(0d0,1d0,po2x < po2th*thon)
        ! print *,mineral,sp_name
        ! omega = 1d0 
        ! domega_dmsp = 0d0
        
        ! omega = 1d0 - po2x/(po2x + keqg3)
        
        select case(trim(adjustl(sp_name)))
            case('po2')
                domega_dmsp = ( &
                    & - 1d0/(po2x + mo2g3) &
                    & - po2x*(-1d0)/(po2x + mo2g3)**2d0 &
                    & ) &
                    & *merge(0d0,1d0,po2x < po2th*thon)
            case default 
                domega_dmsp = 0d0
        endselect 
        
        ! print*, omega
        
        
    case default 
        ! print *,'non-specified'
        omega = 1d0
        domega_dmsp = 0d0
        ! print *,omega
        
endselect

omega_error = .false.
if (any(isnan(omega)) .or. any(isnan(domega_dmsp))) then 
    print *,'nan in calc_omega_dev_v2',any(isnan(omega)),any(isnan(domega_dmsp)),mineral,sp_name
    omega_error = .true.
    ! stop
endif 

endsubroutine calc_omega_dev_v2

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
    & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
    & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
    & ,rxn_name,sp_name &! input 
    & ,rxn_ext,drxnext_dmsp,rxnext_error &! output
    & )
implicit none
integer,intent(in)::nz
real(kind=8):: po2th,fe2th,mwtom,g1th,g2th,g3th,mvpy,fe3th
real(kind=8),dimension(nz):: po2x,vmax,mo2,fe2x,koxa,vmax2,mom2,komb,beta,omx,ombx &
    & ,mo2g1,mo2g2,mo2g3,kg1,kg2,kg3,g1x,g2x,g3x,pyx,fe3x,koxpy
real(kind=8),dimension(nz),intent(in):: poro,sat,hr,prox
real(kind=8),dimension(nz),intent(out):: drxnext_dmsp
real(kind=8),dimension(nz),intent(out):: rxn_ext
character(5),intent(in)::rxn_name,sp_name
logical,intent(out)::rxnext_error

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

integer,intent(in)::nsp_sld,nsp_sld_cnst,nsp_sld_all

character(5),dimension(nsp_sld),intent(in)::chrsld
character(5),dimension(nsp_sld_cnst),intent(in)::chrsld_cnst
character(5),dimension(nsp_sld_all),intent(in)::chrsld_all

real(kind=8),dimension(nsp_sld,nz),intent(in)::msldx
real(kind=8),dimension(nsp_sld_cnst,nz),intent(in)::msldc
real(kind=8),dimension(nsp_sld_all),intent(in)::msldth_all,mv_all

real(kind=8),intent(in)::rho_grain,kw

! real(kind=8):: thon = 1d0
real(kind=8):: thon = -1d100

vmax = krxn1_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)
mo2 = krxn2_ext_all(findloc(chrrxn_ext_all,'resp',dim=1),:)

po2th = mgasth_all(findloc(chrgas_all,'po2',dim=1))

koxa = krxn1_ext_all(findloc(chrrxn_ext_all,'fe2o2',dim=1),:) 

fe2th = maqth_all(findloc(chraq_all,'fe2',dim=1))

fe3th = maqth_all(findloc(chraq_all,'fe3',dim=1))

g1th = msldth_all(findloc(chrsld_all,'g1',dim=1))
g2th = msldth_all(findloc(chrsld_all,'g2',dim=1))
g3th = msldth_all(findloc(chrsld_all,'g3',dim=1))



vmax2 = krxn1_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)
mom2 = krxn2_ext_all(findloc(chrrxn_ext_all,'omomb',dim=1),:)

komb = krxn1_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)
beta = krxn2_ext_all(findloc(chrrxn_ext_all,'ombto',dim=1),:)

! kg1 = krxn1_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:)
! mo2g1 = krxn2_ext_all(findloc(chrrxn_ext_all,'g1dec',dim=1),:)

! kg2 = krxn1_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:)
! mo2g2 = krxn2_ext_all(findloc(chrrxn_ext_all,'g2dec',dim=1),:)

! kg3 = krxn1_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:)
! mo2g3 = krxn2_ext_all(findloc(chrrxn_ext_all,'g3dec',dim=1),:)

koxpy = krxn1_ext_all(findloc(chrrxn_ext_all,'pyfe3',dim=1),:)

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

fe3x = 0d0
if (any(chraq=='fe3')) then 
    fe3x = maqx(findloc(chraq,'fe3',dim=1),:)
elseif (any(chraq_cnst=='fe3')) then 
    fe3x = maqc(findloc(chraq_cnst,'fe3',dim=1),:)
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

pyx = 0d0
if (any(chrsld=='py')) then 
    pyx = msldx(findloc(chrsld,'py',dim=1),:)
elseif (any(chraq_cnst=='py')) then 
    pyx = msldc(findloc(chrsld_cnst,'py',dim=1),:)
endif 

mvpy = mv_all(findloc(chrsld_all,'py',dim=1))

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
        rxn_ext = ( &
            & poro*sat*1d3*fe2x*po2x &
            & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
            & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
            & )
        
        select case(trim(adjustl(sp_name)))
            case('pro')
                drxnext_dmsp = ( &
                    & poro*sat*1d3*fe2x*po2x &
                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)*2.0d0*(kw*(-1d0)/prox**2d0)) &
                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                    & )
            case('po2')
                drxnext_dmsp = ( &
                    & poro*sat*1d3*fe2x*1d0 &
                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                    & )
            case('fe2')
                drxnext_dmsp = ( &
                    & poro*sat*1d3*1d0*po2x &
                    & *(8.0d13*60.0d0*24.0d0*365.0d0*(kw/prox)**2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) &
                    & *merge(0d0,1d0,po2x < po2th*thon .or. fe2x < fe2th*thon) &
                    & )
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
    
    ! case('g1dec')
        ! rxn_ext = kg1*g1x*po2x/(po2x+mo2g1) &
            ! & *merge(0d0,1d0,po2x < po2th .or. g1x < g1th)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! drxnext_dmsp = (&
                    ! & kg1*g1x*1d0/(po2x+mo2g1) &
                    ! & +kg1*g1x*po2x*(-1d0)/(po2x+mo2g1)**2d0 &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g1x < g1th)
            ! case('g1')
                ! drxnext_dmsp = (&
                    ! & kg1*1d0*po2x/(po2x+mo2g1) &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g1x < g1th)
            ! case default
                ! drxnext_dmsp = 0d0
        ! endselect 
    
    ! case('g2dec')
        ! rxn_ext = kg2*g2x*po2x/(po2x+mo2g2)&
            ! & *merge(0d0,1d0,po2x < po2th .or. g2x < g2th)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! drxnext_dmsp = (&
                    ! & kg2*g2x*1d0/(po2x+mo2g2) &
                    ! & +kg2*g2x*po2x*(-1d0)/(po2x+mo2g2)**2d0 &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g2x < g2th)
            ! case('g2')
                ! drxnext_dmsp = (&
                    ! & kg2*1d0*po2x/(po2x+mo2g2) &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g2x < g2th)
            ! case default
                ! drxnext_dmsp = 0d0
        ! endselect 
    
    ! case('g3dec')
        ! rxn_ext = kg3*g3x*po2x/(po2x+mo2g3)&
            ! & *merge(0d0,1d0,po2x < po2th .or. g3x < g3th)
        
        ! select case(trim(adjustl(sp_name)))
            ! case('po2')
                ! drxnext_dmsp = (&
                    ! & kg3*g3x*1d0/(po2x+mo2g3) &
                    ! & +kg3*g3x*po2x*(-1d0)/(po2x+mo2g3)**2d0 &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g3x < g3th)
            ! case('g3')
                ! drxnext_dmsp = (&
                    ! & kg3*1d0*po2x/(po2x+mo2g3) &
                    ! & ) &
                    ! & *merge(0d0,1d0,po2x < po2th .or. g3x < g3th)
            ! case default
                ! drxnext_dmsp = 0d0
        ! endselect 
                
    
    case('pyfe3') 
        rxn_ext = ( &
            & koxpy*poro*hr*mvpy*pyx*fe3x**0.93d0*fe2x**(-0.40d0) &
            & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
            & )
        
        select case(trim(adjustl(sp_name)))
            case('py')
                drxnext_dmsp = (&
                    & koxpy*poro*hr*mvpy*1d0*fe3x**0.93d0*fe2x**(-0.40d0) &
                    & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                    & )
            case('fe3')
                drxnext_dmsp = (&
                    & koxpy*poro*hr*mvpy*pyx*(0.93d0)*fe3x**(0.93d0-1d0)*fe2x**(-0.40d0) &
                    & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                    & )
            case('fe2')
                drxnext_dmsp = (&
                    & koxpy*poro*hr*mvpy*pyx*fe3x**0.93d0*(-0.4d0)*fe2x**(-0.40d0-1d0) &
                    & *merge(0d0,1d0,fe3x<fe3th*thon .or. fe2x<fe2th*thon) &
                    & )
            case default
                drxnext_dmsp = 0d0
        endselect 
    
                
    case default 
        rxn_ext = 0d0
        
endselect

rxnext_error = .false.
if (any(isnan(rxn_ext)) .or. any(isnan(drxnext_dmsp))) then 
    print *,'nan in calc_rxn_ext_dev'
    rxnext_error = .true.
    ! stop
endif 

endsubroutine calc_rxn_ext_dev

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine make_transmx(  &
    & labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  &! input
    & ,trans,nonlocal,izml  &! output 
    & )
implicit none
integer,intent(in)::nsp_sld,nz
real(kind=8),intent(in)::dz(nz),poro(nz),z(nz),zml_ref,dbl_ref,tol
real(kind=8)::sporo(nz)
logical,intent(in)::labs(nsp_sld),turbo2(nsp_sld),nobio(nsp_sld),fick(nsp_sld),till(nsp_sld),save_trans
real(kind=8),intent(out)::trans(nz,nz,nsp_sld)
logical,intent(out)::nonlocal(nsp_sld)
integer,intent(out)::izml
integer iz,isp,iiz,izdbl
real(kind=8) :: translabs(nz,nz),dbio(nz),transdbio(nz,nz),transturbo2(nz,nz),transtill(nz,nz)
real(kind=8) :: zml(nsp_sld),probh,dbl
character(10) chr

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
    do iz = max(1,izdbl), izml
        if (iz==max(1,izdbl)) then 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        elseif (iz==izml) then 
            transdbio(iz,iz) = 0.5d0*(sporo(Iz)*dbio(iz)+sporo(Iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        else 
            transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        endif
    enddo
    ! do iz = max(1,izdbl), izml
        ! if (iz==max(1,izdbl)) then 
            ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! elseif (iz==izml) then 
            ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
        ! else 
            ! transdbio(iz,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  &
                ! + 0.5d0*(dbio(iz)+dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
            ! transdbio(iz-1,iz) = 0.5d0*(dbio(iz)+dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)))
            ! transdbio(iz+1,iz) = 0.5d0*(dbio(iz)+dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)))
        ! endif
    ! enddo
    
    ! Added; changes have been made here rather than in solving governing eqs.
    do iz=1,nz
        transdbio(:,iz) = transdbio(:,iz)/dz(iz)
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
    transturbo2 = 0d0
    probh = 0.0001d0
    do iz=1,izml 
        do iiz=1,izml
            if (iiz/=iz) then 
                transturbo2(iiz,iz) = probh!*dz(iz)/dz(iiz)
                transturbo2(iiz,iiz) = transturbo2(iiz,iiz) - transturbo2(iiz,iz)
            endif 
        enddo
    enddo
    
    ! trying inverse mixing 
    transtill = 0d0
    probh = 0.010d0
    probh = 0.10d0
    do iz=1,izml  ! when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
        ! transtill(iz,iz)=-probh*dz(iz)/dz(izml+1-iz) !*(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
        ! transtill(izml+1-iz,iz)=probh*dz(iz)/dz(izml+1-iz) ! *(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
        ! do iiz = izml+1-iz,iz+1,-1
            ! transtill(iiz,iz)= probh*dz(iz)/dz(iiz)*(iiz/real(izml+1-iz,kind=8))
        ! enddo 
        ! transtill(iz,iz) = -sum(transtill(:,iz))
        do iiz=1,izml
            if (iiz/=iz) then 
                if (iiz==iz-1 .or. iiz == iz+1 .or. iiz == izml + 1 - iz ) then 
                    transtill(iiz,iz)= probh !*dz(iz)/dz(iiz) 
                    ! transtill(iiz,iiz) = transtill(iiz,iiz) - transtill(iiz,iz)
                endif 
            endif 
        enddo 
        transtill(iz,iz) = -sum(transtill(:,iz))
    enddo
    

    ! if (turbo2(isp)) translabs = transturbo2   ! translabs temporarily used to represents nonlocal mixing 
    
    ! added 
    do iz =1,nz
        do iiz= 1,nz
            translabs(iiz,iz) = translabs(iiz,iz)/dz(iz)*dz(iiz)
            transturbo2(iiz,iz) = transturbo2(iiz,iz)/dz(iz)*dz(iiz)
            transtill(iiz,iz) = transtill(iiz,iz)/dz(iz)*dz(iiz)
        enddo 
    enddo 
    
    trans(:,:,isp) = 0d0 
    
    if (nobio(isp)) cycle
    
    if (fick(isp)) then 
        trans(:,:,isp) = trans(:,:,isp) + transdbio(:,:)
    endif 
    
    if (turbo2(isp)) then 
        trans(:,:,isp) = trans(:,:,isp) + transturbo2(:,:)
    endif 
    
    if (labs(isp)) then 
        trans(:,:,isp) = trans(:,:,isp) + translabs(:,:)
    endif 
    
    if (till(isp)) then 
        trans(:,:,isp) = trans(:,:,isp) + transtill(:,:)
    endif 
    
    ! if (any(abs(sum(trans(:,:,isp),dim=1))>tol)) then 
        ! print *, 'transition matrix can be non-conservative'
        ! write(chr,'(i3.3)') isp
        ! open(unit=88,file='./mtx-'//trim(adjustl(chr))//'.txt',action='write',status='replace')
        ! do iz=1,nz
            ! write(88,*) (trans(iz,iiz,isp),iiz=1,nz)  ! writing 
        ! enddo
        ! close(88)
        ! stop
    ! endif 
    
    if (save_trans) then 
        write(chr,'(i3.3)') isp
        open(unit=88,file='./mtx-'//trim(adjustl(chr))//'.txt',action='write',status='replace')
        do iz=1,nz
            write(88,*) (trans(iz,iiz,isp),iiz=1,nz)  ! writing 
        enddo
        close(88)
    endif 
    
    ! trans(:,:,isp) = transdbio(:,:)  !  firstly assume local mixing implemented by dbio 

    ! if (nonlocal(isp)) trans(:,:,isp) = translabs(:,:)  ! if nonlocal, replaced by either turbo2 mixing or labs mixing 
    ! if (nobio(isp)) trans(:,:,isp) = 0d0  ! if assuming no bioturbation, transition matrix is set at zero  
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
    & ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s &
    & ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &
    & ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all &
    & ,turbo2,labs,trans,method_precalc,display,chrflx,sld_enforce &! input
    !  old inputs
    & ,hr,poro,z,dz,w,sat,pro,poroprev,tora,v,tol,it,nflx,kw & 
    & ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0  &
    ! old inout
    & ,dt,flgback &    
    ! output 
    & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp & 
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,tol,kw,ucv,rho_grain,rg,tc,sec2yr,tempk_0
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,tora,v,poroprev,dz,torg,pro
real(kind=8),dimension(nz),intent(out)::prox
integer,intent(inout)::it
integer iter
logical,intent(in)::cplprec,method_precalc,display
logical,intent(inout)::flgback
real(kind=8),intent(in)::dt
real(kind=8) error

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
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,msldsupp 
real(kind=8),dimension(nsp_sld,nz),intent(inout)::ksld
real(kind=8),dimension(nz,nz,nsp_sld),intent(in)::trans
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx,omega,nonprec,rxnsld
real(kind=8),dimension(nsp_sld,nz)::domega_dpro,dmsld,dksld_dpro,drxnsld_dmsld,dksld_dso4f,domega_dso4f
real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq,dksld_dmaq,drxnsld_dmaq
real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas,dksld_dmgas,drxnsld_dmgas
real(kind=8),dimension(nsp_sld,nflx,nz),intent(out)::flx_sld
real(kind=8),dimension(nsp_aq,nz),intent(in)::maq,maqsupp
real(kind=8),dimension(nsp_aq,nz),intent(inout)::maqx 
real(kind=8),dimension(nsp_aq,nz)::dprodmaq,dmaq,dso4fdmaq 
real(kind=8),dimension(nsp_aq,nflx,nz),intent(out)::flx_aq
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,mgassupp
real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx 
real(kind=8),dimension(nsp_gas,nz)::khgasx,khgas,dgas,agasx,agas,rxngas,dkhgas_dpro,dprodmgas,dmgas,dso4fdmgas
real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
real(kind=8),dimension(nsp_gas,nflx,nz),intent(out)::flx_gas 
real(kind=8),dimension(nrxn_ext,nz),intent(inout)::rxnext
real(kind=8),dimension(nrxn_ext,nz)::drxnext_dpro
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
real(kind=8),dimension(nsp_sld_all,nsp_aq_all),intent(in)::staq_all
real(kind=8),dimension(nsp_gas_all,3),intent(in)::keqgas_h
real(kind=8),dimension(nsp_aq_all,4),intent(in)::keqaq_h
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_c
real(kind=8),dimension(nsp_aq_all,2),intent(in)::keqaq_s
real(kind=8),dimension(nsp_sld_all),intent(in)::keqsld_all,msldth_all,mv_all

character(5),dimension(nflx),intent(in)::chrflx

integer ieqgas_h0,ieqgas_h1,ieqgas_h2
data ieqgas_h0,ieqgas_h1,ieqgas_h2/1,2,3/

integer ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4
data ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4/1,2,3,4/

integer ieqaq_co3,ieqaq_hco3
data ieqaq_co3,ieqaq_hco3/1,2/

integer ieqaq_so4,ieqaq_so42
data ieqaq_so4,ieqaq_so42/1,2/

integer,intent(in)::nrxn_ext_all

character(5),dimension(nrxn_ext_all),intent(in)::chrrxn_ext_all

real(kind=8),dimension(nsp_gas_all),intent(in)::mgasth_all
real(kind=8),dimension(nsp_aq_all),intent(in)::maqth_all
real(kind=8),dimension(nrxn_ext_all,nz),intent(in)::krxn1_ext_all,krxn2_ext_all

real(kind=8),dimension(4,nflx,nz),intent(out)::flx_co2sp

integer iz,row,ie,ie2,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn,isps2,iiz
integer::itflx,iadv,idif,irain,ires
data itflx,iadv,idif,irain/1,2,3,4/

integer,dimension(nsp_sld)::irxn_sld 
integer,dimension(nrxn_ext)::irxn_ext 

real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,omega_tmp_th,rxn_ext_tmp &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,khco2n_tmp,pco2n_tmp,edifn_tmp,caqsupp_tmp,kco2,k1,k2,kho,sw_red &
    & ,flx_max,flx_max_max,pco2i,proi

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::fact = 1d-3
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 3d0
real(kind=8)::corr = 1.5d0
real(kind=8),dimension(nz)::dummy,dummy2,so4f,kin,dkin_dmsp  

logical print_cb,ph_error,omega_error,rxnext_error
character(500) print_loc
character(20) chrfmt

integer,parameter :: iter_max = 50
! integer,parameter :: iter_max = 300

integer :: nz_disp = 10

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz),emx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

logical::chkflx = .true.
logical::dt_norm = .true.
logical::kin_iter = .true.

! logical::sld_enforce = .false.
logical,intent(in)::sld_enforce != .true.

character(10) precstyle 
real(kind=8) msld_seed ,fact2
real(kind=8):: fact_tol = 1d-3
real(kind=8):: dt_th = 1d-3
real(kind=8) flx_tol != tol*fact_tol*(z(nz)+0.5d0*dz(nz))
integer solve_sld 

!-----------------------------------------------


precstyle = 'def'
! precstyle = 'full'
! precstyle = 'full_lim'
! precstyle = 'seed '

msld_seed = 1d-20

! flx_tol = tol*fact_tol*(z(nz)+0.5d0*dz(nz))
flx_tol = 1d-4

if (sld_enforce) then 
    solve_sld = 0
else
    solve_sld = 1
endif 

sw_red = 1d0
if (.not.method_precalc) sw_red = -1d100
sw_red = -1d100

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

    
nonprec = 1d0 ! primary minerals only dissolve
if (cplprec)then
    do isps = 1, nsp_sld
        if (any(chrsld_2 == chrsld(isps))) then  
            nonprec(isps,:) = 0d0 ! allowing precipitation for secondary phases
        endif 
    enddo
endif 

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

! print *, 'starting silciate calculation'

do while ((.not.isnan(error)).and.(error > tol*fact_tol))

    amx3=0.0d0
    ymx3=0.0d0 
    emx3=0.0d0 
    
    flx_sld = 0d0
    flx_aq = 0d0
    flx_gas = 0d0
    
    ! pH calculation and its derivative wrt aq and gas species
    
    call calc_pH_v7_2( &
        & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
        & ,print_cb,print_loc,z &! input 
        & ,prox,ph_error,so4f &! output
        & ) 
    if (ph_error) then 
        flgback = .true.
        return
    endif 
    
    dprodmaq = 0d0
    dso4fdmaq = 0d0
    do ispa=1,nsp_aq
        dmaq = 0d0
        dmgas = 0d0
        if (any (chraq_ph == chraq(ispa))) then 
            dmaq(ispa,:) = dconc
            dummy = 0d0
            dummy2 = 0d0
            call calc_pH_v7_2( &
                & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
                & ,maqx+dmaq,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
                & ,print_cb,print_loc,z &! input 
                & ,dummy,ph_error,dummy2 &! output
                & ) 
            dprodmaq(ispa,:) = (dummy - prox)/dconc
            dso4fdmaq(ispa,:) = (dummy2 - so4f)/dconc
        endif 
    enddo 
    
    dprodmgas = 0d0
    dso4fdmgas = 0d0
    do ispg=1,nsp_gas
        dmaq = 0d0
        dmgas = 0d0
        if (any (chrgas_ph == chrgas(ispg))) then 
            dmgas(ispg,:) = dconc
            dummy = 0d0
            dummy2 = 0d0
            call calc_pH_v7_2( &
                & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
                & ,maqx,maqc,mgasx+dmgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
                & ,print_cb,print_loc,z &! input 
                & ,dummy,ph_error,dummy2 &! output
                & ) 
            dprodmgas(ispg,:) = (dummy - prox)/dconc
            dso4fdmgas(ispg,:) = (dummy2 - so4f)/dconc
        endif 
    enddo
    
    ! recalculation of rate constants for mineral reactions
    if (kin_iter) then 
        ksld = 0d0
        dksld_dpro = 0d0
        dksld_dso4f = 0d0
        dksld_dmaq = 0d0
        dksld_dmgas = 0d0
        do isps =1,nsp_sld 
            call sld_kin( &
                & nz,rg,tc,sec2yr,tempk_0,prox,poro,hr,kw,kho,mv(isps) &! input
                & ,chrsld(isps),'pro  ' &! input 
                & ,kin,dkin_dmsp &! output
                & ) 
            ksld(isps,:) = kin
            dksld_dpro(isps,:) = dkin_dmsp
            
            do ispa = 1,nsp_aq
                if (any (chraq_ph == chraq(ispa)) .or. staq(isps,ispa)/=0d0 ) then 
                    call sld_kin( &
                        & nz,rg,tc,sec2yr,tempk_0,prox,poro,hr,kw,kho,mv(isps) &! input
                        & ,chrsld(isps),chraq(ispa) &! input 
                        & ,kin,dkin_dmsp &! output
                        & ) 
                    dksld_dmaq(isps,ispa,:) = dkin_dmsp + ( &
                        & dksld_dpro(isps,:)*dprodmaq(ispa,:) &
                        & +dksld_dso4f(isps,:)*dso4fdmaq(ispa,:) &
                        & )
                endif 
            enddo 
            
            do ispg = 1,nsp_gas
                if (any (chrgas_ph == chrgas(ispg)) .or. stgas(isps,ispg)/=0d0) then 
                    call sld_kin( &
                        & nz,rg,tc,sec2yr,tempk_0,prox,poro,hr,kw,kho,mv(isps) &! input
                        & ,chrsld(isps),chrgas(ispg) &! input 
                        & ,kin,dkin_dmsp &! output
                        & ) 
                    dksld_dmgas(isps,ispg,:) = dkin_dmsp + ( &
                        & dksld_dpro(isps,:)*dprodmgas(ispg,:) &
                        & +dksld_dso4f(isps,:)*dso4fdmgas(ispg,:) &
                        & )
                endif 
            enddo 
        
        enddo 
    else 
        dksld_dpro = 0d0
        dksld_dso4f = 0d0
        dksld_dmaq = 0d0
        dksld_dmgas = 0d0
    endif 
    ! saturation state calc. and their derivatives wrt aq and gas species
    
    ! print *,'ksld',ksld(findloc(chrsld,'gt',dim=1),:)
    
    omega = 0d0
    domega_dpro = 0d0
    domega_dso4f = 0d0
    domega_dmaq = 0d0
    domega_dmgas = 0d0
    
    do isps =1, nsp_sld
    
        dummy = 0d0
        dummy2 = 0d0
        call calc_omega_dev_v2( &
            & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
            & ,prox,chrsld(isps),'pro  ' &! input 
            & ,dummy,dummy2,omega_error &! output
            & )
        if (omega_error) then
            flgback = .true.
            return 
        endif 
        omega(isps,:) = dummy
        domega_dpro(isps,:) = dummy2
        
    
        dummy = 0d0
        dummy2 = 0d0
        call calc_omega_dev_v2( &
            & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
            & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
            & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
            & ,prox,chrsld(isps),'so4f ' &! input 
            & ,dummy,dummy2,omega_error &! output
            & )
        if (omega_error) then
            flgback = .true.
            return 
        endif 
        domega_dso4f(isps,:) = dummy2
        
        do ispa = 1, nsp_aq
            if (any (chraq_ph == chraq(ispa)) .or. staq(isps,ispa)/=0d0 ) then 
                
                dummy = 0d0
                dummy2 = 0d0
                call calc_omega_dev_v2( &
                    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
                    & ,prox,chrsld(isps),chraq(ispa) &! input 
                    & ,dummy,dummy2,omega_error &! output
                    & )
                if (omega_error) then
                    flgback = .true.
                    return 
                endif 
                
                domega_dmaq(isps,ispa,:) = dummy2 + ( &
                    & domega_dpro(isps,:)*dprodmaq(ispa,:) &
                    & +domega_dso4f(isps,:)*dso4fdmaq(ispa,:) &
                    & )
            endif 
        enddo
        do ispg = 1, nsp_gas
            if (any (chrgas_ph == chrgas(ispg)) .or. stgas(isps,ispg)/=0d0) then 
                
                dummy = 0d0
                dummy2 = 0d0
                call calc_omega_dev_v2( &
                    & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
                    & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
                    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
                    & ,prox,chrsld(isps),chrgas(ispg) &! input 
                    & ,dummy,dummy2,omega_error &! output
                    & )
                if (omega_error) then
                    flgback = .true.
                    return 
                endif 
                
                domega_dmgas(isps,ispg,:) = dummy2 + ( &
                    & domega_dpro(isps,:)*dprodmgas(ispg,:) &
                    & +domega_dso4f(isps,:)*dso4fdmgas(ispg,:) &
                    & )
            endif 
        enddo
    enddo 
    
    
    ! print *,'omega',omega(findloc(chrsld,'gt',dim=1),:)
    ! print *
    ! print *,'domega_dmaq',domega_dmaq(findloc(chrsld,'gt',dim=1),:,:)
    ! print *
    ! print *,'domega_dmgas',domega_dmgas(findloc(chrsld,'gt',dim=1),:,:)
    
    ! adding reactions that are not based on dis/prec of minerals
    rxnext = 0d0
    drxnext_dpro = 0d0
    drxnext_dmaq = 0d0
    drxnext_dmgas = 0d0
    drxnext_dmsld = 0d0
    
    do irxn=1,nrxn_ext
        dummy = 0d0
        dummy2 = 0d0
        call calc_rxn_ext_dev( &
            & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
            & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
            & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
            & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
            & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
            & ,chrrxn_ext(irxn),'pro  ' &! input 
            & ,dummy,dummy2,rxnext_error &! output
            & )
        if (rxnext_error) then
            flgback = .true.
            return 
        endif 
        rxnext(irxn,:) = dummy
        drxnext_dpro(irxn,:) = dummy2
        
        do ispg=1,nsp_gas
            if (stgas_dext(irxn,ispg)==0d0) cycle
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
                & ,chrrxn_ext(irxn),chrgas(ispg) &! input 
                & ,dummy,dummy2,rxnext_error &! output
                & )
            if (rxnext_error) then
                flgback = .true.
                return 
            endif 
            drxnext_dmgas(irxn,ispg,:) = dummy2 + drxnext_dpro(irxn,:)*dprodmgas(ispg,:)
        enddo 
        
        do ispa=1,nsp_aq
            if (staq_dext(irxn,ispa)==0d0) cycle
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
                & ,chrrxn_ext(irxn),chraq(ispa) &! input 
                & ,dummy,dummy2,rxnext_error &! output
                & )
            if (rxnext_error) then
                flgback = .true.
                return 
            endif 
            drxnext_dmaq(irxn,ispa,:) = dummy2 + drxnext_dpro(irxn,:)*dprodmaq(ispa,:)
        enddo 
        
        do isps=1,nsp_sld
            if (stsld_dext(irxn,isps)==0d0) cycle
            
            dummy = 0d0
            dummy2 = 0d0
            call calc_rxn_ext_dev( &
                & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
                & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
                & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
                & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
                & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
                & ,chrrxn_ext(irxn),chrsld(isps) &! input 
                & ,dummy,dummy2,rxnext_error &! output
                & )
            if (rxnext_error) then
                flgback = .true.
                return 
            endif 
            drxnext_dmsld(irxn,isps,:) = dummy2
        enddo 
    enddo 
    
    ! gas tansport
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
    enddo 
    
    ! sld phase reactions
    
    rxnsld = 0d0
    drxnsld_dmsld = 0d0
    drxnsld_dmaq = 0d0
    drxnsld_dmgas = 0d0
    
    call sld_rxn( &
        & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
        & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle &! input
        & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas &! input
        & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
        & ) 
    
    ! gas reactions 
    
    rxngas = 0d0
    drxngas_dmaq = 0d0
    drxngas_dmsld = 0d0
    drxngas_dmgas = 0d0
        
    do ispg = 1, nsp_gas
        do isps = 1, nsp_sld
            rxngas(ispg,:) =  rxngas(ispg,:) + (&
                ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & + stgas(isps,ispg)*rxnsld(isps,:) &
                & )
            drxngas_dmsld(ispg,isps,:) =  drxngas_dmsld(ispg,isps,:) + (&
                ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & + stgas(isps,ispg)*drxnsld_dmsld(isps,:) &
                & )
            do ispg2 = 1,nsp_gas
                drxngas_dmgas(ispg,ispg2,:) =  drxngas_dmgas(ispg,ispg2,:) + (&
                    ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg2,:)) &
                    ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    ! & + stgas(isps,ispg)*dksld_dmgas(isps,ispg2,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + stgas(isps,ispg)*drxnsld_dmgas(isps,ispg2,:) &
                    & )
            enddo 
            do ispa = 1,nsp_aq
                drxngas_dmaq(ispg,ispa,:) =  drxngas_dmaq(ispg,ispa,:) + ( &
                    ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                    ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    ! & + stgas(isps,ispg)*dksld_dmaq(isps,ispa,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + stgas(isps,ispg)*drxnsld_dmaq(isps,ispa,:) &
                    & )
            enddo 
        enddo 
    enddo 
            
    if (.not.sld_enforce) then 

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

                amx3(row,row) = ((1d0-poro(iz))/merge(1d0,dt,dt_norm)     &
                    & + w/dz(iz)*merge(dt,1d0,dt_norm)    &
                    & + drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                    & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                    & ) &
                    & * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red)

                ymx3(row) = ( &
                    & ((1d0-poro(iz))*m_tmp-(1d0-poroprev(iz))*mprev_tmp)/merge(1d0,dt,dt_norm) &
                    & -w*(mp_tmp-m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & + rxnsld(isps,iz)*merge(dt,1d0,dt_norm) &
                    & -msupp_tmp*merge(dt,1d0,dt_norm)  &
                    & -rxn_ext_tmp*merge(dt,1d0,dt_norm)  &
                    & ) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    
                if (iz/=nz) amx3(row,row+nsp3) = (-w/dz(iz))*merge(dt,1d0,dt_norm) *merge(1.0d0,mp_tmp,m_tmp<mth_tmp*sw_red)
                
                do ispa = 1, nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa
                    
                    amx3(row,col ) = ( &
                        & + drxnsld_dmaq(isps,ispa,iz)*merge(dt,1d0,dt_norm) &
                        & - sum(stsld_ext(:,isps)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *maqx(ispa,iz) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                enddo 
                
                do ispg = 1, nsp_gas 
                    col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                    amx3(row,col) = ( &
                        & + drxnsld_dmgas(isps,ispg,iz)*merge(dt,1d0,dt_norm) &
                        & - sum(stsld_ext(:,isps)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *mgasx(ispg,iz) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                enddo 
                
                do isps2 = 1,nsp_sld 
                    if (isps2 == isps) cycle
                    col = nsp3*(iz-1)+ isps2

                    amx3(row,col) = ( &
                        & - sum(stsld_ext(:,isps)*drxnext_dmsld(:,isps2,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *msldx(isps2,iz) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                enddo 
                
                ! diffusion terms are filled with transition matrices 
                ! if (turbo2(isps).or.labs(isps)) then
                    ! do iiz = 1, nz
                        ! col = nsp3*(iiz-1)+isps
                        ! if (trans(iiz,iz,isps)==0d0) cycle
                        ! amx3(row,col) = amx3(row,col) &
                            ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)
                        ! ymx3(row) = ymx3(row) &
                            ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz)
                            
                        ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                            ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                            ! & )
                    ! enddo
                ! else
                    ! do iiz = 1, nz
                        ! col = nsp3*(iiz-1)+isps
                        ! if (trans(iiz,iz,isps)==0d0) cycle
                            
                        ! amx3(row,col) = amx3(row,col) -trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                            ! & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        ! ymx3(row) = ymx3(row) - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                            ! & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                            
                        ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                            ! & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                            ! & )
                    ! enddo
                ! endif
                
                ! modifications with porosity and dz are made in make_trans subroutine
                do iiz = 1, nz
                    col = nsp3*(iiz-1)+isps
                    if (trans(iiz,iz,isps)==0d0) cycle
                        
                    amx3(row,col) = amx3(row,col) -trans(iiz,iz,isps)*msldx(isps,iiz)*merge(dt,1d0,dt_norm) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                    ymx3(row) = ymx3(row) - trans(iiz,iz,isps)*msldx(isps,iiz)*merge(dt,1d0,dt_norm) &
                        & *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red)
                        
                    flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        & - trans(iiz,iz,isps)*msldx(isps,iiz) &
                        & )
                enddo
                
                flx_sld(isps,itflx,iz) = ( &
                    ! & (m_tmp-mprev_tmp)/dt &
                    & ((1d0-poro(iz))*m_tmp-(1d0-poroprev(iz))*mprev_tmp)/dt &
                    & )
                flx_sld(isps,iadv,iz) = ( &
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_sld(isps,irxn_sld(isps),iz) = ( &
                    & + rxnsld(isps,iz) &
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
    
    endif 

    do iz = 1, nz
        
        do ispa = 1, nsp_aq

            row = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa
            
            d_tmp = daq(ispa)
            caq_tmp = maqx(ispa,iz)
            caq_tmp_prev = maq(ispa,iz)
            caq_tmp_p = maqx(ispa,min(nz,iz+1))
            caq_tmp_n = maqx(ispa,max(1,iz-1))
            caqth_tmp = maqth(ispa)
            caqi_tmp = maqi(ispa)
            caqsupp_tmp = maqsupp(ispa,iz) 
            rxn_ext_tmp = sum(staq_ext(:,ispa)*rxnext(:,iz))
            rxn_tmp = sum(staq(:,ispa)*rxnsld(:,iz))
            drxndisp_tmp = sum(staq(:,ispa)*drxnsld_dmaq(:,ispa,iz))
            
            if (iz==1) caq_tmp_n = caqi_tmp
                
            edif_tmp = 1d3*poro(iz)*sat(iz)*tora(iz)*d_tmp
            edif_tmp_p = 1d3*poro(min(iz+1,nz))*sat(min(iz+1,nz))*tora(min(iz+1,nz))*d_tmp
            edif_tmp_n = 1d3*poro(max(iz-1,1))*sat(max(iz-1,1))*tora(max(iz-1,1))*d_tmp

            amx3(row,row) = ( &
                & (poro(iz)*sat(iz)*1d3*1d0)/merge(1d0,dt,dt_norm)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*merge(0d0,-1d0,iz==nz)/( 0.5d0*(dz(iz)+dz(min(nz,iz+1))) ) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(1d0)/( 0.5d0*(dz(iz)+dz(max(1,iz-1))) ))/dz(iz) &
                & *merge(dt,1d0,dt_norm) &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                & -drxndisp_tmp*merge(dt,1d0,dt_norm) &
                & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp*sw_red)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/merge(1d0,dt,dt_norm)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                & *merge(dt,1d0,dt_norm) &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*merge(dt,1d0,dt_norm) &
                & - rxn_tmp*merge(dt,1d0,dt_norm) &
                & - caqsupp_tmp*merge(dt,1d0,dt_norm) &
                & - rxn_ext_tmp*merge(dt,1d0,dt_norm) &
                & ) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)

            if (iz/=1) then 
                amx3(row,row-nsp3) = ( &
                    & -(-0.5d0*(edif_tmp +edif_tmp_n)*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & *merge(dt,1d0,dt_norm) &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(-1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            endif 
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz) &
                    & *merge(dt,1d0,dt_norm) &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            endif 
            
            if (.not.sld_enforce) then 
                do isps = 1, nsp_sld
                    col = nsp3*(iz-1)+ isps
                    
                    amx3(row, col) = (     & 
                        ! & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz)) &
                        ! & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)*merge(dt,1d0,dt_norm)  &
                        & - staq(isps,ispa)*drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) &
                        & - sum(staq_ext(:,ispa)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *msldx(isps,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
                enddo 
            endif  
            
            do ispa2 = 1, nsp_aq
                col = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa2
                
                if (ispa2 == ispa) cycle
                
                amx3(row,col) = amx3(row,col) + (     & 
                    & - sum(staq(:,ispa)*drxnsld_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) &
                    & - sum(staq_ext(:,ispa)*drxnext_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *maqx(ispa2,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
            enddo 
            
            do ispg = 1, nsp_gas
                col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg
                
                amx3(row,col) = amx3(row,col) + (     & 
                    & - sum(staq(:,ispa)*drxnsld_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                    & - sum(staq_ext(:,ispa)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red)   ! commented out (is this necessary?)
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
                ! & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
                & - staq(:,ispa)*rxnsld(:,iz) &
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
            
            amx3(row,:) = amx3(row,:)*fact 
            ymx3(row) = ymx3(row)*fact 
        
        enddo 
        
    end do  ! ==============================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2 & pO2   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! print *,drxngas_dmaq(findloc(chrgas,'pco2',dim=1),findloc(chraq,'ca',dim=1),:)
    
    do iz = 1, nz
        
        do ispg = 1, nsp_gas
        
            row = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg
            
            pco2n_tmp = mgasx(ispg,max(1,iz-1))
            khco2n_tmp = khgasx(ispg,max(1,iz-1))
            edifn_tmp = dgas(ispg,max(1,iz-1))
            if (iz == 1) then 
                pco2n_tmp = mgasi(ispg)
                khco2n_tmp = khgasi(ispg)
                edifn_tmp = dgasi(ispg)
            endif 

            amx3(row,row) = ( &
                & (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  &
                & *merge(dt,1d0,dt_norm) &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz) *merge(dt,1d0,dt_norm) &
                & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) &
                & -drxngas_dmgas(ispg,ispg,iz)*merge(dt,1d0,dt_norm) &
                & ) &
                & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
            
            ymx3(row) = ( &
                & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/merge(1d0,dt,dt_norm) &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
                & *merge(dt,1d0,dt_norm) &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*merge(dt,1d0,dt_norm) &
                ! & -resp(iz) &
                & -sum(stgas_ext(:,ispg)*rxnext(:,iz))*merge(dt,1d0,dt_norm) &
                & -rxngas(ispg,iz)*merge(dt,1d0,dt_norm) &
                & -mgassupp(ispg,iz)*merge(dt,1d0,dt_norm) &
                & ) &
                & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg)*sw_red)
            
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                        & + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa 
                    amx3(row,col+nsp3) = ( &
                        & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                            
                enddo 
            
            endif 
            
            if (iz/=1) then 
                amx3(row,row-nsp3) = ( &
                    & -(- 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                    & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                    &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-khgasx(ispg,iz-1)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg,iz-1)*mgasx(ispg,iz-1))/dz(iz)*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *merge(0.0d0,mgasx(ispg,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa 

                    amx3(row,col-nsp3) = ( &
                        & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                        &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*merge(dt,1d0,dt_norm)  &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                            
                    if (row==2 .and. col==3) then
                        print*,row,col,'iz-1,dgas/daq'
                    endif 
                enddo 
            endif 
            
            if (.not.sld_enforce) then 
                do isps = 1,nsp_sld
                    col = nsp3*(iz-1) + isps 
                    amx3(row,col) = ( &
                        & -drxngas_dmsld(ispg,isps,iz)*merge(dt,1d0,dt_norm) &
                        & -sum(stgas_ext(:,ispg)*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) &
                        & ) &
                        & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                enddo 
            endif 
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa 
                amx3(row,col) = ( &
                    & (dagas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                    & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & - 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  &
                    & *merge(dt,1d0,dt_norm) &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                    & -drxngas_dmaq(ispg,ispa,iz)*merge(dt,1d0,dt_norm) &
                    & -sum(stgas_ext(:,ispg)*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
            enddo 
            
            do ispg2 = 1, nsp_gas
                if (ispg == ispg2) cycle
                col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg2
                amx3(row,col) = ( &
                    & (dagas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) &
                    & -( 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                    &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & - 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,iz)-pco2n_tmp) &
                    &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*merge(dt,1d0,dt_norm)  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) &
                    & -drxngas_dmgas(ispg,ispg2,iz)*merge(dt,1d0,dt_norm) &
                    & -sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg2,iz))*merge(dt,1d0,dt_norm) &
                    & ) &
                    & *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red)
                
            enddo 
            
            if (amx3(row,row)==0d0) then 
                print *,amx3(row,row),mgasx(ispg,iz)<mgasth(ispg)*sw_red,mgasx(ispg,iz) 
                print *, &
                & (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz)) &
                & ,-( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  &
                & ,+poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz) &
                & ,+poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz) &
                & ,-sum(stgas_ext(:,ispg)*drxnext_dmgas(:,ispg,iz)) &
                & ,-drxngas_dmgas(ispg,ispg,iz) 
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
            flx_gas(ispg,irxn_sld(:),iz) = ( &
                ! & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
                ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0)  &
                & - stgas(:,ispg)*rxnsld(:,iz) &
                & ) 
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
            
            if (any(isnan(flx_gas(ispg,:,iz)))) then
                print *,flx_gas(ispg,:,iz)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/alpha(iz)
            ! ymx3(row) = ymx3(row)/alpha(iz)
        enddo 

    end do 
    
    fact2= maxval(abs(amx3))
    
    amx3 = amx3/fact2
    ymx3 = ymx3/fact2
    
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
        if (.not.sld_enforce) then 
            do isps = 1, nsp_sld
                row = isps + nsp3*(iz-1)

                if (isnan(ymx3(row))) then 
                    print *,'nan at', iz,z(iz),chrsld(isps)
                    stop
                endif
                
                emx3(row) = (1d0-poro(iz))*msldx(isps,iz)*exp(ymx3(row)) -(1d0-poro(iz))*msldx(isps,iz)

                if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                    msldx(isps,iz) = msldx(isps,iz)*corr
                else if (ymx3(row) < -threshold) then 
                    msldx(isps,iz) = msldx(isps,iz)/corr
                else   
                    msldx(isps,iz) = msldx(isps,iz)*exp(ymx3(row))
                endif
                
                if ( msldx(isps,iz)<msldth(isps)) then ! too small trancate value and not be accounted for error 
                    msldx(isps,iz)=msldth(isps)
                    ymx3(row) = 0d0
                endif
            enddo 
        endif 
        
        do ispa = 1, nsp_aq
            row = ispa + nsp_sld*solve_sld + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chraq(ispa)
                stop
            endif
            
            emx3(row) = poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*exp(ymx3(row)) - poro(iz)*sat(iz)*1d3*maqx(ispa,iz)

            if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                maqx(ispa,iz) = maqx(ispa,iz)*corr
            else if (ymx3(row) < -threshold) then 
                maqx(ispa,iz) = maqx(ispa,iz)/corr
            else   
                maqx(ispa,iz) = maqx(ispa,iz)*exp(ymx3(row))
            endif
            
            if (maqx(ispa,iz)<maqth(ispa)) then ! too small trancate value and not be accounted for error 
                maqx(ispa,iz)=maqth(ispa)
                ymx3(row) = 0d0
            endif
        enddo 
        
        do ispg = 1, nsp_gas
            row = ispg + nsp_aq + nsp_sld*solve_sld + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chrgas(ispg)
                stop
            endif
            
            emx3(row) =agasx(ispg,iz)* mgasx(ispg,iz)*exp(ymx3(row)) - agasx(ispg,iz)*mgasx(ispg,iz) 

            if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
                mgasx(ispg,iz) = mgasx(ispg,iz)*corr
            else if (ymx3(row) < -threshold) then 
                mgasx(ispg,iz) = mgasx(ispg,iz)/corr
            else   
                mgasx(ispg,iz) = mgasx(ispg,iz)*exp(ymx3(row))
            endif
            
            if (mgasx(ispg,iz)<mgasth(ispg)) then ! too small trancate value and not be accounted for error 
                mgasx(ispg,iz)=mgasth(ispg)
                ymx3(row) = 0d0
            endif
        enddo 
        
    end do 

    if (fact_tol == 1d0) then 
        error = maxval(exp(abs(ymx3))) - 1.0d0
    else 
        error = maxval((abs(emx3)))
    endif 
    
    if (isnan(error)) error = 1d4

    if (isnan(error).or.info/=0 .or. any(isnan(msldx)) .or. any(isnan(maqx)).or. any(isnan(mgasx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(msldx)),any(isnan(maqx)),any(isnan(mgasx))'
        print*,isnan(error),info,any(isnan(msldx)),any(isnan(maqx)),any(isnan(mgasx))
        
        open(unit=11,file='amx.txt',status = 'replace')
        open(unit=12,file='ymx.txt',status = 'replace')
        do ie = 1,nsp3*(nz)
            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
            write(12,*) ymx3(ie)
        enddo 
        close(11)
        close(12)      
        
        ! dt = dt/10d0
        flgback = .true.
        ! pause
        exit
        
        
        ! stop
    endif

    if (display) then 
        print '(a,E11.3,a,i0,a,E11.3)', 'iteration error = ',error, ', iteration = ',iter,', time step [yr] = ',dt
    endif      
    iter = iter + 1 

    if (iter > iter_Max ) then
    ! if (iter > iter_Max .or. (method_precalc .and. error > infinity)) then
        ! dt = dt/1.01d0
        ! dt = dt/10d0
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
        write(chrfmt,'(i0)') nz_disp
        chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
        
        print *
        print *,' [concs] '
        print trim(adjustl(chrfmt)),'z',(z(iz),iz=1,nz,nz/nz_disp)
        if (nsp_aq>0) then 
            print *,' < aq species >'
            do ispa = 1, nsp_aq
                print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_sld>0) then 
            print *,' < sld species >'
            do isps = 1, nsp_sld
                print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        if (nsp_gas>0) then 
            print *,' < gas species >'
            do ispg = 1, nsp_gas
                print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/nz_disp)
            enddo 
        endif 
        print *
#endif     

enddo

! just addint flx calculation at the end 

    
flx_sld = 0d0
flx_aq = 0d0
flx_gas = 0d0

flx_co2sp = 0d0

! pH calculation and its derivative wrt aq and gas species

call calc_pH_v7_2( &
    & nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
    & ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all &!input
    & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all &! input
    & ,print_cb,print_loc,z &! input 
    & ,prox,ph_error,so4f &! output
    & ) 

    
! recalculation of rate constants for mineral reactions

if (kin_iter) then 

    ksld = 0d0

    do isps =1,nsp_sld 
        call sld_kin( &
            & nz,rg,tc,sec2yr,tempk_0,prox,poro,hr,kw,kho,mv(isps) &! input
            & ,chrsld(isps),'pro  ' &! input 
            & ,kin,dkin_dmsp &! output
            & ) 
        ksld(isps,:) = kin
    enddo 

endif 

! saturation state calc. and their derivatives wrt aq and gas species

omega = 0d0

do isps =1, nsp_sld
    dummy = 0d0
    dummy2 = 0d0
    call calc_omega_dev_v2( &
        & nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst &! input 
        & ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all &!input
        & ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,mgasth_all,keqaq_s,so4f,staq_all &! input
        & ,prox,chrsld(isps),'pro  ' &! input 
        & ,dummy,dummy2,omega_error &! output
        & )
    if (omega_error) then
        ! dt = dt/10d0
        flgback = .true.
        exit 
    endif 
    omega(isps,:) = dummy
enddo 

rxnsld = 0d0
    
call sld_rxn( &
    & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
    & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle &! input
    & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas &! input
    & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
    & ) 

! adding reactions that are not based on dis/prec of minerals
rxnext = 0d0

do irxn=1,nrxn_ext
    dummy = 0d0
    dummy2 = 0d0
    call calc_rxn_ext_dev( &
        & nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  &!input
        & ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst &! input
        & ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all &! input
        & ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw &!input
        & ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox &! input
        & ,chrrxn_ext(irxn),'pro  ' &! input 
        & ,dummy,dummy2,rxnext_error &! output
        & )
    if (rxnext_error) then
        flgback = .true.
        exit 
    endif 
    rxnext(irxn,:) = dummy
enddo 

if (.not.sld_enforce)then 
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
            ! if (turbo2(isps).or.labs(isps)) then
                ! do iiz = 1, nz
                    ! if (trans(iiz,iz,isps)==0d0) cycle
                        
                    ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        ! & - trans(iiz,iz,isps)/dz(iz)*dz(iiz)*msldx(isps,iiz) &
                        ! & )
                ! enddo
            ! else
                ! do iiz = 1, nz
                    ! if (trans(iiz,iz,isps)==0d0) cycle
                        
                    ! flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                        ! & - trans(iiz,iz,isps)/dz(iz)*msldx(isps,iiz) &
                        ! & )
                ! enddo
            ! endif
            
            do iiz = 1, nz
                if (trans(iiz,iz,isps)==0d0) cycle
                    
                flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( &
                    & - trans(iiz,iz,isps)*msldx(isps,iiz) &
                    & )
            enddo
            
            flx_sld(isps,itflx,iz) = ( &
                ! & (m_tmp-mprev_tmp)/dt &
                & ((1d0-poro(iz))*m_tmp-(1d0-poroprev(iz))*mprev_tmp)/dt &
                & )
            flx_sld(isps,iadv,iz) = ( &
                & -w*(mp_tmp-m_tmp)/dz(iz)  &
                & )
            ! flx_sld(isps,irxn_sld(isps),iz) = ( &
                ! & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                ! & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                ! & )
            flx_sld(isps,irxn_sld(isps),iz) = ( &
                & + rxnsld(isps,iz) &
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
endif 

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
        rxn_tmp = sum(staq(:,ispa)*rxnsld(:,iz))
        
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
            ! & -staq(:,ispa)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
            ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
            & - staq(:,ispa)*rxnsld(:,iz) &
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
            ! & stgas(isps,ispg)*ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
            ! & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
            & + stgas(isps,ispg)*rxnsld(isps,:) &
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
        flx_gas(ispg,irxn_sld(:),iz) = ( &
            ! & -stgas(:,ispg)*ksld(:,iz)*poro(iz)*hr(iz)*mv(:)*1d-6*msldx(:,iz)*(1d0-omega(:,iz)) &
            ! & *merge(0d0,1d0,1d0-omega(:,iz)*nonprec(:,iz) < 0d0) &
            & - stgas(:,ispg)*rxnsld(:,iz) &
            & )
        flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
        
        if (any(isnan(flx_gas(ispg,:,iz)))) then
            print *,flx_gas(ispg,:,iz)
        endif 
    enddo 
    
    if (any(chrgas=='pco2')) then 
        ispg = findloc(chrgas,'pco2',dim=1)
        pco2i = mgasi(ispg)
        proi = prox(max(1,iz-1))
        if (iz==1) proi = sqrt(kco2*k1*pco2i)
        
        ! gaseous CO2
        
        edifn_tmp = ucv*poro(max(1,iz-1))*(1.0d0-sat(max(1,iz-1)))*1d3*torg(max(1,iz-1))*dgasg(ispg)
        if (iz==1) edifn_tmp = dgasi(ispg)
        
        flx_co2sp(1,itflx,iz) = ( &
            & (ucv*poro(iz)*(1.0d0-sat(Iz))*1d3*mgasx(ispg,iz)-ucv*poroprev(iz)*(1.0d0-sat(Iz))*1d3*mgas(ispg,iz))/dt &
            & )  
        flx_co2sp(1,idif,iz) = ( &
            & -( 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) &
            &       +ucv*poro(min(nz,iz+1))*(1.0d0-sat(min(nz,iz+1)))*1d3*torg(min(nz,iz+1))*dgasg(ispg)) &
            &       *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) + edifn_tmp) &
            &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & ) 
        flx_co2sp(1,irxn_ext(:),iz) = -stgas_ext(:,ispg)*rxnext(:,iz)
        flx_co2sp(1,irain,iz) = - mgassupp(ispg,iz)
        flx_co2sp(1,irxn_sld(:),iz) = ( &
            & - stgas(:,ispg)*rxnsld(:,iz) &
            & )
            
        ! dissolved CO2
        
        edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*1d3*tora(max(1,iz-1))*dgasa(ispg)
        if (iz==1) edifn_tmp = 0d0
        
        flx_co2sp(2,itflx,iz) = ( &
            & (poro(iz)*sat(iz)*kco2*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*1d3*mgas(ispg,iz))/dt &
            & )  
        flx_co2sp(2,idif,iz) = ( &
            & -( 0.5d0*(poro(iz)*sat(iz)*kco2*1d3*tora(iz)*dgasa(ispg) &
            &       +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*1d3*tora(min(nz,iz+1))*dgasa(ispg)) &
            &       *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(poro(iz)*sat(iz)*kco2*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) &
            &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & ) 
        flx_co2sp(2,iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(kco2*mgasx(ispg,iz)- kco2*pco2n_tmp)/dz(iz) &
            & )
            
        ! HCO3-
        
        edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*k1/prox(max(1,iz-1))*1d3*tora(max(1,iz-1))*dgasa(ispg)
        if (iz==1) edifn_tmp = 0d0
        
        flx_co2sp(3,itflx,iz) = ( &
            & (poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*k1/pro(iz)*1d3*mgas(ispg,iz))/dt &
            & )  
        flx_co2sp(3,idif,iz) = ( &
            & -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*tora(iz)*dgasa(ispg) &
            &       +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*k1/prox(min(nz,iz+1))*1d3*tora(min(nz,iz+1))*dgasa(ispg)) &
            &       *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) &
            &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & ) 
        flx_co2sp(3,iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*( &
            &       kco2*k1/prox(iz)*mgasx(ispg,iz) &
            &       - kco2*k1/proi*pco2n_tmp)/dz(iz) &
            & )
            
        ! CO32-
        
        edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*k1*k2/prox(max(1,iz-1))**2d0*1d3*tora(max(1,iz-1))*dgasa(ispg)
        if (iz==1) edifn_tmp = 0d0
        
        flx_co2sp(4,itflx,iz) = (  &
            & (poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*mgasx(ispg,iz) &
            &       -poroprev(iz)*sat(iz)*kco2*k1*k2/pro(iz)**2d0*1d3*mgas(ispg,iz))/dt &
            & )  
        flx_co2sp(4,idif,iz) = ( &
            & -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*tora(iz)*dgasa(ispg) &
            &       +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*k1*k2/prox(min(nz,iz+1))**2d0*1d3*tora(min(nz,iz+1))*dgasa(ispg)) &
            &       *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)**2d0*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) &
            &       *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  &
            & ) 
        flx_co2sp(4,iadv,iz) = ( &
            & +poro(iz)*sat(iz)*v(iz)*1d3*( &
            &       kco2*k1*k2/prox(iz)**2d0*mgasx(ispg,iz) &
            &       - kco2*k1*k2/proi**2d0*pco2n_tmp)/dz(iz) &
            & )
            
            
        flx_co2sp(1,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
        flx_co2sp(2,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
        flx_co2sp(3,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
        flx_co2sp(4,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz))
    endif 
    
end do 

    
#ifdef dispiter
        
print *
print *,' [saturation & pH] '
if (nsp_sld>0) then 
    print *,' < sld species omega >'
    do isps = 1, nsp_sld
        print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/nz_disp)
    enddo 
endif 
print *,' < pH >'
print trim(adjustl(chrfmt)), 'ph', (-log10(prox(iz)),iz=1,nz, nz/nz_disp)
print *

write(chrfmt,'(i0)') nflx
chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,a11))'

print *
print *,' [fluxes] '
print trim(adjustl(chrfmt)),'time',(chrflx(iflx),iflx=1,nflx)

write(chrfmt,'(i0)') nflx
chrfmt = '(a5,'//trim(adjustl(chrfmt))//'(1x,E11.3))'
if (nsp_aq>0) then 
    print *,' < aq species >'
    do ispa = 1, nsp_aq
        print trim(adjustl(chrfmt)), trim(adjustl(chraq(ispa))), (sum(flx_aq(ispa,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
endif 
if (nsp_sld>0) then 
    print *,' < sld species >'
    do isps = 1, nsp_sld
        print trim(adjustl(chrfmt)), trim(adjustl(chrsld(isps))), (sum(flx_sld(isps,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
endif 
if (nsp_gas>0) then 
    print *,' < gas species >'
    do ispg = 1, nsp_gas
        print trim(adjustl(chrfmt)), trim(adjustl(chrgas(ispg))), (sum(flx_gas(ispg,iflx,:)*dz(:)),iflx=1,nflx)
    enddo 
endif 
print *
#endif     

if (chkflx .and. dt > dt_th) then 
    flx_max_max = 0d0
    do isps = 1, nsp_sld

        flx_max = 0d0
        do iflx = 1, nflx
            flx_max = max(flx_max,abs(sum(flx_sld(isps,iflx,:)*dz)))
        enddo 
        
        flx_max_max = max(flx_max_max,flx_max)
    enddo 

    do ispa = 1, nsp_aq

        flx_max = 0d0
        do iflx = 1, nflx
            flx_max = max(flx_max,abs(sum(flx_aq(ispa,iflx,:)*dz)))
        enddo 
        flx_max_max = max(flx_max_max,flx_max)
    enddo 

    do ispg = 1, nsp_gas

        flx_max = 0d0
        do iflx = 1, nflx
            flx_max = max(flx_max,abs(sum(flx_gas(ispg,iflx,:)*dz)))
        enddo 
        flx_max_max = max(flx_max_max,flx_max)
    enddo 
    
    if (.not.sld_enforce) then 
        do isps = 1, nsp_sld

            flx_max = 0d0
            do iflx = 1, nflx
                flx_max = max(flx_max,abs(sum(flx_sld(isps,iflx,:)*dz)))
            enddo 
            
            if (flx_max/flx_max_max > 1d-9 .and.  abs(sum(flx_sld(isps,ires,:)*dz))/flx_max > flx_tol ) then 
                print *, 'too large error in mass balance of sld phases'
                print *,chrsld(isps),abs(sum(flx_sld(isps,ires,:)*dz)),flx_max
                flgback = .true.
                return
            
                open(unit=11,file='amx.txt',status = 'replace')
                open(unit=12,file='ymx.txt',status = 'replace')
                do ie = 1,nsp3*(nz)
                    write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                    write(12,*) ymx3(ie)
                enddo 
                close(11)
                close(12)     
                
                stop
                ! dt = dt/10d0
            endif 
        enddo 
    endif 

    do ispa = 1, nsp_aq

        flx_max = 0d0
        do iflx = 1, nflx
            flx_max = max(flx_max,abs(sum(flx_aq(ispa,iflx,:)*dz)))
        enddo 
        
        if (flx_max/flx_max_max > 1d-9  .and. abs(sum(flx_aq(ispa,ires,:)*dz))/flx_max > flx_tol ) then 
            print *, 'too large error in mass balance of aq phases'
            print *,chraq(ispa),abs(sum(flx_aq(ispa,ires,:)*dz)),flx_max
            flgback = .true.
            return
        
            open(unit=11,file='amx.txt',status = 'replace')
            open(unit=12,file='ymx.txt',status = 'replace')
            do ie = 1,nsp3*(nz)
                write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                write(12,*) ymx3(ie)
            enddo 
            close(11)
            close(12)     
            
            stop
            ! dt = dt/10d0
        endif 
    enddo 

    do ispg = 1, nsp_gas

        flx_max = 0d0
        do iflx = 1, nflx
            flx_max = max(flx_max,abs(sum(flx_gas(ispg,iflx,:)*dz)))
        enddo 
        
        if (flx_max/flx_max_max > 1d-9  .and. abs(sum(flx_gas(ispg,ires,:)*dz))/flx_max > flx_tol ) then 
            print *, 'too large error in mass balance of gas phases'
            print *,chrgas(ispg),abs(sum(flx_gas(ispg,ires,:)*dz)),flx_max
            flgback = .true.
            return
        
            open(unit=11,file='amx.txt',status = 'replace')
            open(unit=12,file='ymx.txt',status = 'replace')
            do ie = 1,nsp3*(nz)
                write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
                write(12,*) ymx3(ie)
            enddo 
            close(11)
            close(12)     
            
            ! dt = dt/10d0
            stop
        endif 
    enddo 
endif 

endsubroutine alsilicate_aq_gas_1D_v3_1

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine sld_rxn( &
    & nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz &! input 
    & ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle &! input
    & ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas &! input
    & ,rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas &! output
    & ) 
implicit none 

integer,intent(in)::nz,nsp_sld,nsp_aq,nsp_gas
real(kind=8),intent(in)::msld_seed,dt
real(kind=8),dimension(nz),intent(in)::hr,poro,sat,dz
real(kind=8),dimension(nsp_sld),intent(in)::mv,msldth
real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
real(kind=8),dimension(nsp_aq),intent(in)::maqth
real(kind=8),dimension(nsp_aq,nz),intent(in)::maq
real(kind=8),dimension(nsp_gas),intent(in)::mgasth
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,agas
real(kind=8),dimension(nsp_sld,nz),intent(in)::ksld,omega,nonprec,msldx,msld
real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(in)::dksld_dmaq,domega_dmaq
real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(in)::dksld_dmgas,domega_dmgas
character(10),intent(in)::precstyle
real(kind=8),dimension(nsp_sld,nz),intent(out)::rxnsld,drxnsld_dmsld
real(kind=8),dimension(nsp_sld,nsp_aq,nz),intent(out)::drxnsld_dmaq
real(kind=8),dimension(nsp_sld,nsp_gas,nz),intent(out)::drxnsld_dmgas

integer ispa,isps,ispg,iz
real(kind=8),dimension(nsp_sld,nz)::maxdis,maxprec
    
rxnsld = 0d0
drxnsld_dmsld = 0d0
drxnsld_dmaq = 0d0
drxnsld_dmgas = 0d0

maxdis = 1d200
maxprec = -1d200
do isps=1,nsp_sld
    maxdis(isps,:) = min(maxdis(isps,:),(1d0-poro)*(msld(isps,:)-msldth(isps))/dt)
    do ispa = 1,nsp_aq
        if (staq(isps,ispa)<0d0) then 
            maxdis(isps,:) = min(maxdis(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
        endif 
    enddo 
    do ispg = 1,nsp_gas
        if (stgas(isps,ispg)<0d0) then 
            maxdis(isps,:) = min(maxdis(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
        endif 
    enddo 
    do ispa = 1,nsp_aq
        if (staq(isps,ispa)>0d0) then 
            maxprec(isps,:) = max(maxprec(isps,:), -1d0/staq(isps,ispa)*poro*sat*1d3*(maq(ispa,:)-maqth(ispa))/dt )
        endif 
    enddo 
    do ispg = 1,nsp_gas
        if (stgas(isps,ispg)>0d0) then 
            maxprec(isps,:) = max(maxprec(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:)*(mgas(ispg,:)-mgasth(ispg))/dt )
        endif 
    enddo 
enddo 

do isps = 1,nsp_sld
    select case(trim(adjustl(precstyle)))
    
        case ('full_lim') 
            
            do iz = 1,nz
                if (1d0-omega(isps,iz) > 0d0) then 
                    rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz))
                    if (rxnsld(isps,iz)> maxdis(isps,iz)) then 
                        rxnsld(isps,iz) = maxdis(isps,iz)
                        drxnsld_dmsld(isps,iz) = 0d0
                        drxnsld_dmaq(isps,:,iz) = 0d0
                        drxnsld_dmgas(isps,:,iz) = 0d0
                    else 
                        drxnsld_dmsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz))
                        drxnsld_dmaq(isps,:,iz) = ( &
                            & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,:,iz)) &
                            & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                            & )
                        drxnsld_dmgas(isps,:,iz) = ( &
                            & +ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmgas(isps,:,iz)) &
                            & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) &
                            & )
                    endif 
                elseif (1d0-omega(isps,iz) < 0d0) then 
                    if (nonprec(isps,iz)==1d0) then 
                        rxnsld(isps,iz) = 0d0
                        drxnsld_dmsld(isps,iz) = 0d0
                        drxnsld_dmaq(isps,:,iz) = 0d0
                        drxnsld_dmgas(isps,:,iz) = 0d0
                    elseif (nonprec(isps,iz)==0d0) then 
                        rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz))
                        if (rxnsld(isps,iz) < maxprec(isps,iz)) then 
                            rxnsld(isps,iz) = maxprec(isps,iz)
                            drxnsld_dmsld(isps,iz) = 0d0
                            drxnsld_dmaq(isps,:,iz) = 0d0
                            drxnsld_dmgas(isps,:,iz) = 0d0
                        else
                            rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz))
                            drxnsld_dmsld(isps,iz) = 0d0
                            drxnsld_dmaq(isps,:,iz) = ( &
                                & +ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmaq(isps,:,iz)) &
                                & +dksld_dmaq(isps,:,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                & )
                            drxnsld_dmgas(isps,:,iz) = ( &
                                & +ksld(isps,iz)*poro(iz)*hr(iz)*(-domega_dmgas(isps,:,iz)) &
                                & +dksld_dmgas(isps,:,iz)*poro(iz)*hr(iz)*(1d0-omega(isps,iz)) &
                                & )
                        endif 
                    endif 
                endif 
            enddo 
            
    
        case ('full') 
            rxnsld(isps,:) = ( &
                & + min(ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)),maxdis(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                &  + max(ksld(isps,:)*poro*hr*(1d0-omega(isps,:)),maxprec(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                & )
            
            drxnsld_dmsld(isps,:) = ( &
                & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                & )
                
            do ispa = 1,nsp_aq
                drxnsld_dmaq(isps,ispa,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                    & + dksld_dmaq(isps,ispa,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                    &  + ksld(isps,:)*poro*hr*(-domega_dmaq(isps,ispa,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                    &  + dksld_dmaq(isps,ispa,:)*poro*hr*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                    & )
            enddo 
                
            do ispg = 1,nsp_gas
                drxnsld_dmgas(isps,ispg,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                    & + dksld_dmgas(isps,ispg,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:) < 0d0) &
                    &  + ksld(isps,:)*poro*hr*(-domega_dmgas(isps,ispg,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                    &  + dksld_dmgas(isps,ispg,:)*poro*hr*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*(1d0-nonprec(isps,:)) > 0d0) &
                    & )
            enddo 
            
            
        case('seed')
            rxnsld(isps,:) = ( &
                & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
        
            drxnsld_dmsld(isps,:) = ( &
                & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
            
            do ispa = 1, nsp_aq
                drxnsld_dmaq(isps,ispa,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmaq(isps,ispa,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + dksld_dmaq(isps,ispa,:)*poro*hr*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
            
            do ispg = 1, nsp_gas
                drxnsld_dmgas(isps,ispg,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(-domega_dmgas(isps,ispg,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + dksld_dmgas(isps,ispg,:)*poro*hr*mv(isps)*1d-6*(msldx(isps,:)+msld_seed)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
        
        
        case default
            rxnsld(isps,:) = ( &
                & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
        
            drxnsld_dmsld(isps,:) = ( &
                & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*1d0*(1d0-omega(isps,:)) &
                & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                & )
            
            do ispa = 1, nsp_aq
                drxnsld_dmaq(isps,ispa,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmaq(isps,ispa,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + dksld_dmaq(isps,ispa,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
            
            do ispg = 1, nsp_gas
                drxnsld_dmgas(isps,ispg,:) = ( &
                    & + ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(-domega_dmgas(isps,ispg,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & + dksld_dmgas(isps,ispg,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*(1d0-omega(isps,:)) &
                    & *merge(0d0,1d0,1d0-omega(isps,:)*nonprec(isps,:) < 0d0) &
                    & )
            enddo 
            
            ! print *, 'max-rxnflx', isps, sum (ksld(isps,:)*poro*hr*mv(isps)*1d-6*msldx(isps,:)*dz)
    
    endselect
enddo   
   
endsubroutine sld_rxn

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

#ifdef no_intr_findloc
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