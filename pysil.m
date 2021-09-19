

function weathering
    cwd = pwd;

    [nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc] = get_variables_num
        
    [chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc] ... 
        = get_variables(nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc )
        

    [nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save ...
        ,count_dtunchanged_Max,tc] = get_bsdvalues

    weathering_main( ...
        nz,ztot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,ttot,plant_rain  ...% input
        ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save ...% input
        ,count_dtunchanged_Max,tc ...% input 
        ,nsld_kinspc,chrsld_kinspc,kin_sld_spc ...% input
        )

end


function weathering_main( ...
    nz,ztot,rainpowder,zsupp,poroi,satup0,zsat,zml_ref,w0,q0,p80,ttot,plant_rain  ...% input
    ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,chraq,chrgas,chrsld,chrrxn_ext,sim_name,runname_save ...% input
    ,count_dtunchanged_Max,tcin ...% input 
    ,nsld_kinspc_in,chrsld_kinspc_in,kin_sld_spc_in ...% input 
    )

    %-----------------------------

    nt = 50000000;

    rg = 8.3d-3   ;% kJ mol^-1 K^-1
    rg2 = 8.2d-2  ;% L mol^-1 atm K^-1

    tempk_0 = 273d0 ;
    sec2yr = 60d0*60d0*24d0*365d0 ;

    n2c_g1 = 0.1d0 ;% N to C ratio for OM-G1; Could be related to reactivity cf. Janssen 1996
    n2c_g2 = 0.1d0 ;
    n2c_g3 = 0.1d0 ;

    fr_an_ab = 0.0d0 ;% Anorthite fraction for albite (Beerling et al., 2020); 0.0 - 0.1
    fr_an_olg = 0.2d0 ;% Anorthite fraction for oligoclase (Beerling et al., 2020); 0.1 - 0.3
    fr_an_and = 0.4d0 ;% Anorthite fraction for andesine (Beerling et al., 2020); 0.3 - 0.5
    fr_an_la = 0.6d0 ;% Anorthite fraction for labradorite (Beerling et al., 2020); 0.5 - 0.7
    fr_an_by = 0.8d0 ;% Anorthite fraction for bytownite (Beerling et al., 2020); 0.7 - 0.9
    fr_an_an = 1.0d0 ;% Anorthite fraction for anorthite (Beerling et al., 2020); 0.9 - 1.0

    fr_hb_cpx = 0.5d0 ;% Hedenbergite fraction for clinopyroxene; 0.0 - 1.0
    fr_fer_opx = 0.5d0 ;% Ferrosilite fraction for orthopyroxene; 0.0 - 1.0
    fr_fer_agt = 0.5d0 ;% Ferrosilite (and Hedenbergite; or Fe/(Fe+Mg)) fraction for Augite; 0.0 - 1.0
    fr_opx_agt = 0.5d0 ;% OPX (or Ca/(Fe+Mg)) fraction for Augite; 0.0 - 1.0

    mvka = 99.52d0 ;% cm3/mol; molar volume of kaolinite; Robie et al. 1978
    mvfo = 43.79d0 ;% cm3/mol; molar volume of Fo; Robie et al. 1978
    mvab_0 = 100.07d0 ;% cm3/mol; molar volume of Ab(NaAlSi3O8); Robie et al. 1978 
    mvan_0 = 100.79d0 ;% cm3/mol; molar volume of An (CaAl2Si2O8); Robie et al. 1978
    mvab = fr_an_ab*mvan_0 + (1d0-fr_an_ab)*mvab_0 ;% cm3/mol; molar volume of albite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
    mvan = fr_an_an*mvan_0 + (1d0-fr_an_an)*mvab_0 ;% cm3/mol; molar volume of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
    mvby = fr_an_by*mvan_0 + (1d0-fr_an_by)*mvab_0 ;% cm3/mol; molar volume of bytownite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing 
    mvla = fr_an_la*mvan_0 + (1d0-fr_an_la)*mvab_0 ;% cm3/mol; molar volume of labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mvand = fr_an_and*mvan_0 + (1d0-fr_an_and)*mvab_0 ;% cm3/mol; molar volume of andesine (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mvolg = fr_an_olg*mvan_0 + (1d0-fr_an_olg)*mvab_0 ;% cm3/mol; molar volume of oligoclase (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mvcc = 36.934d0 ;% cm3/mol; molar volume of Cc (CaCO3); Robie et al. 1978
    mvpy = 23.94d0 ;% cm3/mol; molar volume of Pyrite (FeS2); Robie et al. 1978
    mvgb = 31.956d0 ;% cm3/mol; molar volume of Gibsite (Al(OH)3); Robie et al. 1978
    mvct = 108.5d0 ;% cm3/mol; molar volume of Chrysotile (Mg3Si2O5(OH)4); Robie et al. 1978
    mvfa = 46.39d0 ;% cm3/mol; molar volume of Fayalite (Fe2SiO4); Robie et al. 1978
    mvgt = 20.82d0 ;% cm3/mol; molar volume of Goethite (FeO(OH)); Robie et al. 1978
    mvcabd = 129.77d0 ;% cm3/mol; molar volume of Ca-beidellite (Ca(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
    mvkbd = 134.15d0 ;% cm3/mol; molar volume of K-beidellite (K(1/3)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
    mvnabd = 130.73d0 ;% cm3/mol; molar volume of Na-beidellite (Na(1/3)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
    mvmgbd = 128.73d0 ;% cm3/mol; molar volume of Mg-beidellite (Mg(1/6)Al(7/3)Si(11/3)O10(OH)2); Wolery and Jove-Colon 2004
    mvdp = 66.09d0 ;% cm3/mol; molar volume of Diopside (MgCaSi2O6);  Robie et al. 1978
    mvhb = 248.09d0/3.55d0 ;% cm3/mol; molar volume of Hedenbergite (FeCaSi2O6); from a webpage
    mvcpx = fr_hb_cpx*mvhb + (1d0-fr_hb_cpx)*mvdp  ;% cm3/mol; molar volume of clinopyroxene (FexMg(1-x)CaSi2O6); assuming simple ('ideal'?) mixing
    mvkfs = 108.72d0 ;% cm3/mol; molar volume of K-feldspar (KAlSi3O8); Robie et al. 1978
    mvom = 30d0/1.5d0 ;% cm3/mol; molar volume of OM (CH2O); calculated assuming 30 g/mol of molar weight and 1.2 g/cm3 of density (Mayer et al., 2004; Ruhlmann et al.,2006)
    mvomb = 30d0/1.5d0 ;% cm3/mol; assumed to be same as mvom
    mvg1 = 30d0/1.5d0 ;% cm3/mol; assumed to be same as mvom
    mvg2 = 30d0/1.5d0 ;% cm3/mol; assumed to be same as mvom
    mvg3 = 30d0/1.5d0 ;% cm3/mol; assumed to be same as mvom
    mvamsi = 25.739d0 ;% cm3/mol; molar volume of amorphous silica taken as cristobalite (SiO2); Robie et al. 1978
    mvarg = 34.15d0 ;% cm3/mol; molar volume of aragonite; Robie et al. 1978
    mvdlm = 64.34d0 ;% cm3/mol; molar volume of dolomite; Robie et al. 1978
    mvhm = 30.274d0 ;% cm3/mol; molar volume of hematite; Robie et al. 1978
    mvill = 139.35d0 ;% cm3/mol; molar volume of illite (K0.6Mg0.25Al2.3Si3.5O10(OH)2); Wolery and Jove-Colon 2004
    mvanl = 97.49d0 ;% cm3/mol; molar volume of analcime (NaAlSi2O6*H2O); Robie et al. 1978
    mvnph = 54.16d0 ;% cm3/mol; molar volume of nepheline (NaAlSiO4); Robie et al. 1978
    mvqtz = 22.688d0 ;% cm3/mol; molar volume of quartz (SiO2); Robie et al. 1978
    mvgps = 74.69d0 ;% cm3/mol; molar volume of gypsum (CaSO4*2H2O); Robie et al. 1978
    mvtm = 272.92d0 ;% cm3/mol; molar volume of tremolite (Ca2Mg5(Si8O22)(OH)2); Robie et al. 1978
    mven = 31.31d0 ;% cm3/mol; molar volume of enstatite (MgSiO3); Robie and Hemingway 1995
    mvfer = 33.00d0 ;% cm3/mol; molar volume of ferrosilite (FeSiO3); Robie and Hemingway 1995
    mvopx = fr_fer_opx*mvfer +(1d0-fr_fer_opx)*mven ;%  cm3/mol; molar volume of clinopyroxene (FexMg(1-x)SiO3); assuming simple ('ideal'?) mixing
    mvmscv = 140.71d0 ;% cm3/mol; molar volume of muscovite (KAl2(AlSi3O10)(OH)2); Robie et al. 1978
    mvplgp = 149.91d0 ;% cm3/mol; molar volume of phlogopite (KMg3(AlSi3O10)(OH)2); Robie et al. 1978
    mvantp = 274.00d0 ;% cm3/mol; molar volume of anthophyllite (Mg7Si8O22(OH)2); Robie and Bethke 1962
    mvagt = (fr_fer_agt*mvfer +(1d0-fr_fer_agt)*mven)*2d0*fr_opx_agt ... % (Fe2xyMg2(1-x)ySi2yO6y)
        + (fr_fer_agt*mvhb + (1d0-fr_fer_agt)*mvdp)*(1d0-fr_opx_agt) ;% (Fex(1-y)Mg(1-x)(1-y)Ca(1-y)Si2(1-y)O6(1-y))
        %  cm3/mol; molar volume of augite 
        % (Fe(2xy+x(1-y))Mg(2y-2xy+1+xy-x-y)Ca(1-y)Si2O6 = Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6)
        % ; assuming simple ('ideal'?) mixing
                                    
    mwtka = 258.162d0 ;% g/mol; formula weight of Ka; Robie et al. 1978
    mwtfo = 140.694d0 ;% g/mol; formula weight of Fo; Robie et al. 1978
    mwtab_0 = 262.225d0 ;% g/mol; formula weight of Ab; Robie et al. 1978
    mwtan_0 = 278.311d0 ;% g/mol; formula weight of An; Robie et al. 1978
    mwtab = fr_an_ab*mwtan_0 + (1d0-fr_an_ab)*mwtab_0 ;% g/mol; formula weight of albte (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtan = fr_an_an*mwtan_0 + (1d0-fr_an_an)*mwtab_0 ;% g/mol; formula weight of anorthite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtby = fr_an_by*mwtan_0 + (1d0-fr_an_by)*mwtab_0 ;% g/mol; formula weight of bytownite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtla = fr_an_la*mwtan_0 + (1d0-fr_an_la)*mwtab_0 ;% g/mol; formula weight of labradorite (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtand = fr_an_and*mwtan_0 + (1d0-fr_an_and)*mwtab_0 ;% g/mol; formula weight of andesine (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtolg = fr_an_olg*mwtan_0 + (1d0-fr_an_olg)*mwtab_0 ;% g/mol; formula weight of oligoclase (CaxNa(1-x)Al(1+x)Si(3-x)O8); assuming simple ('ideal'?) mixing
    mwtcc = 100.089d0 ;% g/mol; formula weight of Cc; Robie et al. 1978
    mwtpy = 119.967d0 ;% g/mol; formula weight of Py; Robie et al. 1978
    mwtgb = 78.004d0 ;% g/mol; formula weight of Gb; Robie et al. 1978
    mwtct = 277.113d0 ;% g/mol; formula weight of Ct; Robie et al. 1978
    mwtfa = 203.778d0 ;% g/mol; formula weight of Fa; Robie et al. 1978
    mwtgt = 88.854d0 ;% g/mol; formula weight of Gt; Robie et al. 1978
    mwtcabd = 366.6252667d0 ;% g/mol; formula weight of Cabd calculated from atmoic weight
    mwtkbd = 372.9783667d0 ;% g/mol; formula weight of Kbd calculated from atmoic weight
    mwtnabd = 367.6088333d0 ;% g/mol; formula weight of Nabd calculated from atmoic weight
    mwtmgbd = 363.9964333d0 ;% g/mol; formula weight of Mgbd calculated from atmoic weight
    mwtdp = 216.553d0 ;% g/mol;  Robie et al. 1978
    mwthb = 248.09d0 ;% g/mol; from a webpage
    mwtcpx = fr_hb_cpx*mwthb + (1d0-fr_hb_cpx)*mwtdp ;% g/mol; formula weight of clinopyroxene (FexMg(1-x)CaSi2O6); assuming simple ('ideal'?) mixing
    mwtkfs = 278.33d0 ;% g/mol; formula weight of Kfs; Robie et al. 1978
    mwtom = 30d0 ;% g/mol; formula weight of CH2O
    mwtomb = 30d0 ;% g/mol; formula weight of CH2O
    mwtg1 = 30d0 ;% g/mol; formula weight of CH2O
    mwtg2 = 30d0 ;% g/mol; formula weight of CH2O
    mwtg3 = 30d0 ;% g/mol; formula weight of CH2O
    mwtamsi = 60.085d0 ;% g/mol; formula weight of amorphous silica
    mwtarg = 100.089d0 ;% g/mol; formula weight of aragonite
    mwtdlm = 184.403d0 ;% g/mol; formula weight of dolomite
    mwthm = 159.692d0 ;% g/mol; formula weight of hematite
    mwtill = 383.90053d0 ;% g/mol; formula weight of Ill calculated from atmoic weight
    mwtanl = 220.155d0 ;% g/mol; formula weight of analcime
    mwtnph = 142.055d0 ;% g/mol; formula weight of nepheline
    mwtqtz = 60.085d0 ;% g/mol; formula weight of quartz
    mwtgps = 172.168d0 ;% g/mol; formula weight of gypsum
    mwttm = 812.374d0 ;% g/mol; formula weight of tremolite
    mwten = 100.389d0 ;% g/mol; formula weight of enstatite
    mwtfer = 131.931d0 ;% g/mol; formula weight of ferrosilite
    mwtopx = fr_fer_opx*mwtfer + (1d0 -fr_fer_opx)*mwten ;% g/mol; formula weight of clinopyroxene (FexMg(1-x)SiO3); assuming simple ('ideal'?) mixing
    mwtmscv = 398.311d0 ;% g/mol; formula weight of muscovite
    mwtplgp = 417.262d0 ;% g/mol; formula weight of phlogopite
    mwtantp = 780.976d0 ;% g/mol; formula weight of anthophyllite
    mwtagt = (fr_fer_agt*mwtfer +(1d0-fr_fer_agt)*mwten)*2d0*fr_opx_agt ...% (Fe2xyMg2(1-x)ySi2yO6y)
        + (fr_fer_agt*mwthb + (1d0-fr_fer_agt)*mwtdp)*(1d0-fr_opx_agt) ;% (Fex(1-y)Mg(1-x)(1-y)Ca(1-y)Si2(1-y)O6(1-y))
                                    %  g/mol; formula weight of augite 
                                    % (Fe(2xy+x(1-y))Mg(2y-2xy+1+xy-x-y)Ca(1-y)Si2O6 = Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6)
                                    % ; assuming simple ('ideal'?) mixing
                                    
    rho_grain = 2.7d0 ;% g/cm3 as soil grain density

    mvblk = mvka ;% for bulk soil assumed to be equal to kaolinite
    mwtblk = mwtka ;

    zsupp_plant = 0.3d0 ;%  e-folding decrease

    rough_c0 = 10d0^(3.3d0) ;
    rough_c1 = 0.33d0 ;

    tol = 1d-6 ;

    % nrec_prof = 22 ;
    nrec_prof = 20 ;
    nrec_flx = 60 ;

    % maxdt = 10d0 ;
    maxdt = 0.2d0 ;% for basalt exp?

    maxdt_max = 1d2  ;% default   
    % maxdt_max = 1d0   ;% when time step matters a reduced value might work 

    pre_calc = false ; 
    % pre_calc = true ; 

    read_data = false ; 
    % read_data = true ; 

    % incld_rough = false ; 
    incld_rough = true ; 

    % cplprec = false ; 
    cplprec = true ; 

    dust_wave = false ; 
    % dust_wave = true ; 

    al_inhibit = false ; 
    % al_inhibit = true ; 

    timestep_fixed = false ; 
    % timestep_fixed = true ; 

    % no_biot = false ; 
    no_biot = true ; 

    biot_fick = false ; 
    % biot_fick = true ; 

    biot_turbo2 = false ; 
    % biot_turbo2 = true ; 

    biot_labs = false ; 
    % biot_labs = true ; 

    biot_till = false ; 
    % biot_till = true ; 

    % display = false ; 
    display = true ; 

    % regular_grid = false ; 
    regular_grid = true ; 

    % method_precalc = false ; 
    method_precalc = true ; 

    % sld_enforce = false ; 
    sld_enforce = true ; 

    poroevol = false ; 

    surfevol1 = false ; 
    surfevol2 = false ; 

    % noncnstw = false ;% constant uplift rate 
    noncnstw = true  ;% varied with porosity

    display_lim = false ;% limiting display fluxes and concs. 
    % display_lim = true ; 

    dust_step = false ; 
    % dust_step = true ; 

    season = false ; 
    % season = true ; 

    step_tau = 0.1d0 ;% yr time duration during which dust is added
    tol_step_tau = 1d-6 ;% yr time duration during which dust is added

    wave_tau = 2d0 ;% yr periodic time for wave 
    dust_norm = 0d0;
    dust_norm_prev = 0d0;

    iwtype_cnst = 0;
    iwtype_pwcnst = 1;
    iwtype_spwcnst = 2;
    iwtype_flex = 3;

    imixtype_nobio = 0;
    imixtype_fick = 1;
    imixtype_turbo2 = 2;
    imixtype_till = 3;
    imixtype_labs = 4;

    rectime_prof = [1d1,3d1,1d2,3d2,1d3,3d3,1d4,3d4 ...
        ,1d5,2d5,3d5,4d5,5d5,6d5,7d5,8d5,9d5,1d6,1.1d6,1.2d6];
    savetime = 1d3;
    dsavetime = 1d3;

    flgback = false ; 
    flgreducedt = false ; 
    flgreducedt_prev = false ; 

    % #ifdef diss_only
    % integer,parameter::nsp_sld_2 = 0
    % #else
    nsp_sld_2 = 17;
    % #endif 
    nsp_sld_all = 44 ;
    nsp_aq_ph = 10;
    nsp_aq_all = 10;
    nsp_gas_ph = 2;
    nsp_gas_all = 4;
    nrxn_ext_all = 9;
    % an attempt to record psd
    nps = 50 ;% bins for particle size 
    % ps_min = 0.1d-6 ;% min particle size (0.1 um)
    ps_min = 10d-9 ;% min particle size (10 nm)
    % ps_min = 100d-9 ;% min particle size (100 nm)
    ps_max = 10d-3 ;% max particle size (10 mm)
    pi = 4d0*atan(1d0) ;% 
    % psd_th = 1d-12 ;% 
    % psd_th = 1d-3 ;% 
    psd_th = 1d0 ;% 
    tol_dvd = 1d-4 ;% 
    nps_rain_char = 4;
    ps_sigma_std = 1d0;
    % ps_sigma_std = 0.8d0 ;
    % ps_sigma_std = 0.2d0 ;
    nflx_psd = 6 ;
    % do_psd = false ; 
    do_psd = true ; 
    do_psd_norm = true ; 
    % do_psd_norm = false ; 
    % mineral specific psd
    do_psd_full = true ; 
    % do_psd_full = false ; 
    psd_lim_min = true ; 
    % psd_lim_min = false ; 
    % psd_vol_consv = true ; 
    psd_vol_consv = false ; 

    [ieqgas_h0,ieqgas_h1,ieqgas_h2] = deal(1,2,3);

    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);

    [ieqaq_co3,ieqaq_hco3] = deal(1,2);

    [ieqaq_so4,ieqaq_so42] = deal(1,2);

    idust = 15;
    nz_disp = 10;

    [itflx,iadv,idif,irain]=deal(1,2,3,4);

    % msldunit = 'sld ';
    msldunit = 'blk';
    %-------------------------
        
    tc = tcin;
    qin = q0;
    satup = satup0;

    nsp_sld_cnst = nsp_sld_all - nsp_sld;
    nsp_aq_cnst = nsp_aq_all - nsp_aq;
    nsp_gas_cnst = nsp_gas_all - nsp_gas;
    nsp3 = nsp_sld + nsp_aq + nsp_gas;

    isldprof = idust + nsp_sld + nsp_gas + nsp_aq + 1;
    isldprof2 = idust + nsp_sld + nsp_gas + nsp_aq + 2;
    isldprof3 = idust + nsp_sld + nsp_gas + nsp_aq + 3;
    iaqprof = idust + nsp_sld + nsp_gas + nsp_aq + 4;
    igasprof = idust + nsp_sld + nsp_gas + nsp_aq + 5;
    isldsat = idust + nsp_sld + nsp_gas + nsp_aq + 6;
    ibsd = idust + nsp_sld + nsp_gas + nsp_aq + 7;
    irate = idust + nsp_sld + nsp_gas + nsp_aq + 8;
    ipsd = idust + nsp_sld + nsp_gas + nsp_aq + 9;
    ipsdv = idust + nsp_sld + nsp_gas + nsp_aq + 10;
    ipsds = idust + nsp_sld + nsp_gas + nsp_aq + 11;
    ipsdflx = idust + nsp_sld + nsp_gas + nsp_aq + 12;
    isa = idust + nsp_sld + nsp_gas + nsp_aq + 13;

    nflx = 5 + nrxn_ext + nsp_sld;

    for isps=1:nsp_sld
        irxn_sld(isps) = 4+isps;
    end 

    for irxn=1:nrxn_ext
        irxn_ext(irxn) = 4+nsp_sld+irxn;
    end 

    ires = nflx;

    chrflx = string({'tflx'; 'adv'; 'dif' ;'rain'})
    if (nsp_sld > 0) ; chrflx(irxn_sld(:)) = chrsld; end
    if (nrxn_ext > 0) ; chrflx(irxn_ext(:)) = chrrxn_ext; end
    chrflx(nflx) = 'res'

    % print *,chrflx,irxn_sld,chrsld

    % pause

    % define all species and rxns definable in the model 
    % note that rxns here exclude diss(/prec) of mineral 
    % which are automatically included when associated mineral is chosen

    chrsld_all = string(...
        {'fo';'ab';'an';'cc';'ka';'gb';'py';'ct';'fa';'gt';'cabd' ...
        ;'dp';'hb';'kfs';'om';'omb';'amsi';'arg';'dlm';'hm';'ill';'anl';'nph' ...
        ;'qtz';'gps';'tm';'la';'by';'olg';'and';'cpx';'en';'fer';'opx';'kbd' ...
        ;'mgbd';'nabd';'mscv';'plgp';'antp';'agt' ...
        ;'g1';'g2';'g3'}...
        )
    chraq_all = string({'mg';'si';'na';'ca';'al';'fe2';'fe3';'so4';'k';'no3'})
    chrgas_all =string({'pco2';'po2';'pnh3';'pn2o'})
    chrrxn_ext_all = string({'resp';'fe2o2';'omomb';'ombto';'pyfe3';'amo2o';'g2n0';'g2n21';'g2n22'})

    % define the species and rxns explicitly simulated in the model in a fully coupled way
    % should be chosen from definable species & rxn lists above 

    % chrsld = (/'fo';'ab';'an';'cc';'ka'/)
    % chraq = (/'mg';'si';'na';'ca';'al'/)
    % chrgas = (/'pco2';'po2'/)
    % chrrxn_ext = (/'resp'/)

    % define solid species which can precipitate
    % in default, all minerals only dissolve 
    % should be chosen from the chrsld list
    % #ifdef diss_only
    % chrsld_2(:) = ''
    % #else
    chrsld_2 = string({'cc';'ka';'gb';'ct';'gt';'cabd';'amsi';'hm';'ill';'anl';'gps'  ...
        ;'arg';'dlm';'qtz';'mgbd';'nabd';'kbd'}) 
    % #endif 
    % below are species which are sensitive to pH 
    chraq_ph = string({'mg';'si';'na';'ca';'al';'fe2';'fe3';'so4';'k';'no3'})
    chrgas_ph = string({'pco2';'pnh3'})

    chrco2sp = string({'co2g';'co2aq';'hco3';'co3';'DIC';'ALK'})
    
    chraq_cnst = strings(nsp_aq_cnst,1);
    if (nsp_aq_cnst ~= 0)  
        for ispa = 1:nsp_aq_cnst
            for ispa2=1:nsp_aq_all
                if (~any(chraq==chraq_all(ispa2)) && ~any(chraq_cnst==chraq_all(ispa2)))  
                    chraq_cnst(ispa) = chraq_all(ispa2);
                    continue
                end 
            end
        end
        disp(chraq_cnst)
    end 
    
    chrgas_cnst = strings(nsp_gas_cnst,1);
    if (nsp_gas_cnst ~= 0)  
        for ispg = 1: nsp_gas_cnst
            for ispg2=1:nsp_gas_all
                if (~any(chrgas==chrgas_all(ispg2)) && ~any(chrgas_cnst==chrgas_all(ispg2)))  
                    chrgas_cnst(ispg) = chrgas_all(ispg2);
                    continue 
                end 
            end
        end 
        disp(chrgas_cnst)
    end 

    chrsld_cnst = strings(nsp_sld_cnst,1);
    if (nsp_sld_cnst ~= 0)  
        for isps = 1: nsp_sld_cnst
            for isps2=1:nsp_sld_all
                if (~any(chrsld==chrsld_all(isps2)) && ~any(chrsld_cnst==chrsld_all(isps2)))  
                    chrsld_cnst(isps) = chrsld_all(isps2);
                    continue 
                end 
            end
        end 
        disp(chrsld_cnst)
    end 

    % molar volume 

    mv_all = [mvfo,mvab,mvan,mvcc,mvka,mvgb,mvpy,mvct,mvfa,mvgt,mvcabd,mvdp,mvhb,mvkfs,mvom,mvomb,mvamsi ...
        ,mvarg,mvdlm,mvhm,mvill,mvanl,mvnph,mvqtz,mvgps,mvtm,mvla,mvby,mvolg,mvand,mvcpx,mven,mvfer,mvopx ...
        ,mvkbd,mvmgbd,mvnabd,mvmscv,mvplgp,mvantp,mvagt ...
        ,mvg1,mvg2,mvg3];
    mwt_all = [mwtfo,mwtab,mwtan,mwtcc,mwtka,mwtgb,mwtpy,mwtct,mwtfa,mwtgt,mwtcabd,mwtdp,mwthb,mwtkfs,mwtom,mwtomb,mwtamsi ...
        ,mwtarg,mwtdlm,mwthm,mwtill,mwtanl,mwtnph,mwtqtz,mwtgps,mwttm,mwtla,mwtby,mwtolg,mwtand,mwtcpx,mwten,mwtfer,mwtopx ...
        ,mwtkbd,mwtmgbd,mwtnabd,mwtmscv,mwtplgp,mwtantp,mwtagt ...
        ,mwtg1,mwtg2,mwtg3];


    for isps = 1: nsp_sld 
        mv(isps) = mv_all(find(chrsld_all==chrsld(isps)));
        mwt(isps) = mwt_all(find(chrsld_all==chrsld(isps)));
    end 

    % maqi_all = 0d0
        
    def_rain = 1d-20;
    def_pr = 1d-20;
    
    [maqi_all] = get_rainwater(nsp_aq_all,chraq_all,def_rain);
    
    [msldi_all] = get_parentrock(nsp_sld_all,chrsld_all,def_pr);

    % bulk soil concentration     
    mblki = 0d0;
    incld_blk = false;

    % adding the case where input wt% exceeds 100% 
    if ( sum(msldi_all) > 1d0)  
        fprintf ( 'parent rock comp. exceeds 100%% so rescale\n')
        msldi_all = msldi_all/sum(msldi_all);  % now the units are g/g
    % endif 
    % msldi_all = msldi_all/mwt_all*rho_grain*1d6 % converting g/g to mol/sld m3
    % msldi_all = (1d0 - poroi) * msldi_all       % mol/sld m3 to mol/bulk m3 

    % when input is less than 100wt% add bulk soil 
    elseif ( sum(msldi_all) < 1d0 )  
        fprintf ('parent rock comp. is less than 100% so add "bulk soil"\n')
        % sum(msldi_all)  + mblki = 1d0
        incld_blk = true;
        mblki = 1d0 - sum(msldi_all);
        if ( mblki < 0d0 ); mblki = 0d0; end
    end 

    rho_grain_calc = rho_grain;
    msldi_allx = msldi_all./mwt_all'*rho_grain_calc*1d6; %  converting g/g to mol/sld m3
    mblkix = mblki/mwtblk*rho_grain_calc*1d6; %  converting g/g to mol/sld m3
    % msldi_allx = msldi_allx/sum(msldi_allx*mv_all*1d-6)  %  try to make sure volume total must be 1 
    msldi_allx = msldi_allx/( sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 );  %  try to make sure volume total must be 1 (including bulk soil if any)
    mblkix = mblkix/( sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 );  %  try to make sure volume total must be 1 (including bulk soil if any)
    if (msldunit=='blk')  
        msldi_allx = msldi_allx*(1d0 - poroi) ; %  try to make sure volume total must be 1 - poroi
        mblkix = mblkix*(1d0 - poroi)  ;%  try to make sure volume total must be 1 - poroi
    end 
    % then the follwoing must be satisfied
    % (1d0 - poroi)*rho_grain_calc = msldi_all*mwt_all*1d-6
    % (1d0 - poroi) = msldi_all*mv_all*1d-6
    rho_error = 1d4;
    rho_tol = 1d-6;
    % poroi_calc = 1d0 - sum(msldi_allx*mv_all*1d-6);
    poroi_calc = 1d0 - ( sum( msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 ); % corrected for bulk soil if any 
    while (rho_error > rho_tol) 
        rho_grain_calcx = rho_grain_calc;
        
        % rho_grain_calc = sum(msldi_allx(:)*mwt_all(:)*1d-6) % /(1d0-poroi_calc)
        rho_grain_calc = sum(msldi_allx(:).*mwt_all(:)*1d-6)  + mblkix*mwtblk*1d-6 ; % /(1d0-poroi_calc) | corrected for bulk soil 
        if (msldunit=='blk') ; rho_grain_calc = rho_grain_calc / (1d0-poroi_calc); end
         
        msldi_allx = msldi_all./mwt_all'*rho_grain_calc*1d6 ;%  converting g/g to mol/sld m3
        mblkix = mblki/mwtblk*rho_grain_calc*1d6 ;%  converting g/g to mol/sld m3
        % msldi_allx = msldi_allx/sum(msldi_allx*mv_all*1d-6)  %  try to make sure volume total must be 1 
        msldi_allx = msldi_allx/( sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 ) ; %  try to make sure volume total must be 1 | corrected for bulk soil 
        if (msldunit=='blk')   
            msldi_allx = (1d0 - poroi) * msldi_allx  ;%  converting mol/sld m3 to mol/bulk m3
            mblkix = (1d0 - poroi) * mblkix  ;%  converting mol/sld m3 to mol/bulk m3
        end 
        
        % poroi_calc = 1d0 - sum(msldi_allx*mv_all*1d-6)
        poroi_calc = 1d0 - ( sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 ) ;% corrected for bulk soil
        
        rho_error = abs ((rho_grain_calc - rho_grain_calcx)/rho_grain_calc) ;
        
        fprintf('rho error = %7.6e\n', rho_error)
        
    end

    if (msldunit=='sld')  
        fprintf('ideally should be unity:%.d (ideal) vs %7.6e (calc.)\n',1d0,sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6)
        fprintf('grain rho as sums comparison: %7.6e vs %7.6e\n',sum(msldi_allx.*mwt_all'*1d-6) + mblkix*mwtblk*1d-6 ,rho_grain_calc)
    elseif (msldunit=='blk')  
        fprintf('ideally should be unity: %.d (ideal) vs %7.6e (calc.)\n',1d0,( sum(msldi_allx.*mv_all'*1d-6) + mblkix*mvblk*1d-6 )/(1d0 - poroi))
        fprintf('grain rho as sums comparison: %7.6e vs %7.6e\n', ( sum(msldi_allx.*mwt_all'*1d-6) + mblkix*mwtblk*1d-6 )/(1d0 - poroi),rho_grain_calc)
    end 
    % print *,mblkix
    % stop
    msldi_all = msldi_allx;
    mblki = mblkix;
    rho_grain = rho_grain_calc;
    
    [mgasi_all] = get_atm(nsp_gas_all,chrgas_all);

    % constant values are taken from the boundary values specified above 
    mgasc = zeros(nsp_gas_cnst,nz,'double');
    maqc = zeros(nsp_aq_cnst,nz,'double');
    msldc = zeros(nsp_sld_cnst,nz,'double');
    for ispg = 1:nsp_gas_cnst
        mgasc(ispg,:) = mgasi_all(find(chrgas_all==chrgas_cnst(ispg)));
    end 
    for ispa = 1:nsp_aq_cnst
        maqc(ispa,:) = maqi_all(find(chraq_all==chraq_cnst(ispa)));
    end 
    for isps = 1:nsp_sld_cnst
        msldc(isps,:) = msldi_all(find(chrsld_all==chrsld_cnst(isps)));
    end 
    
    % threshould values 
    mgasth_all = ones(nsp_gas_all,1,'double')*1d-200;
    maqth_all = ones(nsp_aq_all,1,'double')*1d-200;
    msldth_all = ones(nsp_sld_all,1,'double')*1d-200;
    
    % passing initial and threshold values to explcit variables 
    msldi = zeros(nsp_sld,1,'double');
    msldth = zeros(nsp_sld,1,'double');
    for isps = 1: nsp_sld    
        % disp( chrsld(isps))
        if (any(chrsld_all == chrsld(isps)))  
            msldi(isps) = msldi_all(find(chrsld_all==chrsld(isps)));
            msldth(isps) = msldth_all(find(chrsld_all==chrsld(isps)));
            % fprintf("%4.3e %4.3e\n" msldi(isps),msldi_all(find(chrsld_all==chrsld(isps))))
        end 
    end 
    maqi = zeros(nsp_aq,1,'double');
    maqth = zeros(nsp_aq,1,'double');
    for ispa = 1: nsp_aq    
        if (any(chraq_all == chraq(ispa)))  
            maqi(ispa) = maqi_all(find(chraq_all==chraq(ispa)));
            maqth(ispa) = maqth_all(find(chraq_all==chraq(ispa)));
        end 
    end 
    mgasi = zeros(nsp_gas,1,'double');
    mgasth = zeros(nsp_gas,1,'double');
    for ispg = 1: nsp_gas    
        if (any(chrgas_all == chrgas(ispg)))  
            mgasi(ispg) = mgasi_all(find(chrgas_all==chrgas(ispg)));
            mgasth(ispg) = mgasth_all(find(chrgas_all==chrgas(ispg)));
        end 
    end 

    % stoichiometry
    % mineral dissolution(/precipitation)
    staq_all = zeros(nsp_sld_all,nsp_aq_all,'double');
    stgas_all = zeros(nsp_sld_all,nsp_gas_all,'double');
    % Forsterite; Mg2SiO4
    staq_all(find(chrsld_all=='fo'), find(chraq_all=='mg')) = 2d0 ;
    staq_all(find(chrsld_all=='fo'), find(chraq_all=='si')) = 1d0 ;
    % Albite; NaAlSi3O8
    % staq_all(find(chrsld_all=='ab'), find(chraq_all=='na')) = 1d0 ;
    % staq_all(find(chrsld_all=='ab'), find(chraq_all=='si')) = 3d0 ;
    % staq_all(find(chrsld_all=='ab'), find(chraq_all=='al')) = 1d0 ;
    % Analcime; NaAlSi2O6*H2O
    staq_all(find(chrsld_all=='anl'), find(chraq_all=='na')) = 1d0 ;
    staq_all(find(chrsld_all=='anl'), find(chraq_all=='si')) = 2d0 ;
    staq_all(find(chrsld_all=='anl'), find(chraq_all=='al')) = 1d0 ;
    % Nepheline; NaAlSiO4
    staq_all(find(chrsld_all=='nph'), find(chraq_all=='na')) = 1d0 ;
    staq_all(find(chrsld_all=='nph'), find(chraq_all=='si')) = 1d0 ;
    staq_all(find(chrsld_all=='nph'), find(chraq_all=='al')) = 1d0 ;
    % K-feldspar; KAlSi3O8
    staq_all(find(chrsld_all=='kfs'), find(chraq_all=='k')) = 1d0 ;
    staq_all(find(chrsld_all=='kfs'), find(chraq_all=='si')) = 3d0 ;
    staq_all(find(chrsld_all=='kfs'), find(chraq_all=='al')) = 1d0 ;
    % Anothite; CaAl2Si2O8
    % staq_all(find(chrsld_all=='an'), find(chraq_all=='ca')) = 1d0 ;
    % staq_all(find(chrsld_all=='an'), find(chraq_all=='si')) = 2d0 ;
    % staq_all(find(chrsld_all=='an'), find(chraq_all=='al')) = 2d0 ;
    % Albite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='ab'), find(chraq_all=='ca')) = fr_an_ab;
    staq_all(find(chrsld_all=='ab'), find(chraq_all=='na')) = 1d0 - fr_an_ab;
    staq_all(find(chrsld_all=='ab'), find(chraq_all=='al')) = 1d0 + fr_an_ab;
    staq_all(find(chrsld_all=='ab'), find(chraq_all=='si')) = 3d0 - fr_an_ab;
    % Anothite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='an'), find(chraq_all=='ca')) = fr_an_an;
    staq_all(find(chrsld_all=='an'), find(chraq_all=='na')) = 1d0 - fr_an_an;
    staq_all(find(chrsld_all=='an'), find(chraq_all=='al')) = 1d0 + fr_an_an;
    staq_all(find(chrsld_all=='an'), find(chraq_all=='si')) = 3d0 - fr_an_an;
    % Labradorite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='la'), find(chraq_all=='ca')) = fr_an_la;
    staq_all(find(chrsld_all=='la'), find(chraq_all=='na')) = 1d0 - fr_an_la;
    staq_all(find(chrsld_all=='la'), find(chraq_all=='al')) = 1d0 + fr_an_la;
    staq_all(find(chrsld_all=='la'), find(chraq_all=='si')) = 3d0 - fr_an_la;
    % Andesine; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='and'), find(chraq_all=='ca')) = fr_an_and;
    staq_all(find(chrsld_all=='and'), find(chraq_all=='na')) = 1d0 - fr_an_and;
    staq_all(find(chrsld_all=='and'), find(chraq_all=='al')) = 1d0 + fr_an_and;
    staq_all(find(chrsld_all=='and'), find(chraq_all=='si')) = 3d0 - fr_an_and;
    % Oligoclase; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='olg'), find(chraq_all=='ca')) = fr_an_olg;
    staq_all(find(chrsld_all=='olg'), find(chraq_all=='na')) = 1d0 - fr_an_olg;
    staq_all(find(chrsld_all=='olg'), find(chraq_all=='al')) = 1d0 + fr_an_olg;
    staq_all(find(chrsld_all=='olg'), find(chraq_all=='si')) = 3d0 - fr_an_olg;
    % Bytownite; CaxNa(1-x)Al(1+x)Si(3-x)O8
    staq_all(find(chrsld_all=='by'), find(chraq_all=='ca')) = fr_an_by;
    staq_all(find(chrsld_all=='by'), find(chraq_all=='na')) = 1d0 - fr_an_by;
    staq_all(find(chrsld_all=='by'), find(chraq_all=='al')) = 1d0 + fr_an_by;
    staq_all(find(chrsld_all=='by'), find(chraq_all=='si')) = 3d0 - fr_an_by;
    % Calcite; CaCO3
    staq_all(find(chrsld_all=='cc'), find(chraq_all=='ca')) = 1d0 ;
    stgas_all(find(chrsld_all=='cc'), find(chrgas_all=='pco2')) = 1d0 ;
    % Kaolinite; Al2Si2O5(OH)4
    staq_all(find(chrsld_all=='ka'), find(chraq_all=='si')) = 2d0 ;
    staq_all(find(chrsld_all=='ka'), find(chraq_all=='al')) = 2d0 ;
    % Gibbsite; Al(OH)3
    staq_all(find(chrsld_all=='gb'), find(chraq_all=='al')) = 1d0 ;
    % Pyrite; FeS2
    staq_all(find(chrsld_all=='py'), find(chraq_all=='fe2')) = 1d0 ;
    staq_all(find(chrsld_all=='py'), find(chraq_all=='so4')) = 2d0 ;
    stgas_all(find(chrsld_all=='py'), find(chrgas_all=='po2')) = -7d0/2d0 ;
    % Chrysotile; Mg3Si2O5(OH)4
    staq_all(find(chrsld_all=='ct'), find(chraq_all=='si')) = 2d0 ;
    staq_all(find(chrsld_all=='ct'), find(chraq_all=='mg')) = 3d0 ;
    % Fayalite; Fe2SiO4
    staq_all(find(chrsld_all=='fa'), find(chraq_all=='si')) = 1d0 ;
    staq_all(find(chrsld_all=='fa'), find(chraq_all=='fe2')) = 2d0 ;
    % Goethite; FeO(OH)
    staq_all(find(chrsld_all=='gt'), find(chraq_all=='fe3')) = 1d0 ;
    % Hematite; Fe2O3
    staq_all(find(chrsld_all=='hm'), find(chraq_all=='fe3')) = 2d0 ;
    % Ca-beidellite; Ca(1/6)Al(7/3)Si(11/3)O10(OH)2
    staq_all(find(chrsld_all=='cabd'), find(chraq_all=='ca')) = 1d0/6d0 ;
    staq_all(find(chrsld_all=='cabd'), find(chraq_all=='al')) = 7d0/3d0 ;
    staq_all(find(chrsld_all=='cabd'), find(chraq_all=='si')) = 11d0/3d0 ;
    % Mg-beidellite; Mg(1/6)Al(7/3)Si(11/3)O10(OH)2
    staq_all(find(chrsld_all=='mgbd'), find(chraq_all=='mg')) = 1d0/6d0 ;
    staq_all(find(chrsld_all=='mgbd'), find(chraq_all=='al')) = 7d0/3d0 ;
    staq_all(find(chrsld_all=='mgbd'), find(chraq_all=='si')) = 11d0/3d0 ;
    % K-beidellite; K(1/3)Al(7/3)Si(11/3)O10(OH)2
    staq_all(find(chrsld_all=='kbd'), find(chraq_all=='k')) = 1d0/3d0 ;
    staq_all(find(chrsld_all=='kbd'), find(chraq_all=='al')) = 7d0/3d0 ;
    staq_all(find(chrsld_all=='kbd'), find(chraq_all=='si')) = 11d0/3d0 ;
    % Na-beidellite; Na(1/3)Al(7/3)Si(11/3)O10(OH)2
    staq_all(find(chrsld_all=='nabd'), find(chraq_all=='na')) = 1d0/3d0 ;
    staq_all(find(chrsld_all=='nabd'), find(chraq_all=='al')) = 7d0/3d0 ;
    staq_all(find(chrsld_all=='nabd'), find(chraq_all=='si')) = 11d0/3d0 ;
    % Illite; K0.6Mg0.25Al2.3Si3.5O10(OH)2
    staq_all(find(chrsld_all=='ill'), find(chraq_all=='k')) = 0.6d0 ;
    staq_all(find(chrsld_all=='ill'), find(chraq_all=='mg')) = 0.25d0 ;
    staq_all(find(chrsld_all=='ill'), find(chraq_all=='al')) = 2.3d0 ;
    staq_all(find(chrsld_all=='ill'), find(chraq_all=='si')) = 3.5d0 ;
    % Diopside (MgCaSi2O6)
    staq_all(find(chrsld_all=='dp'), find(chraq_all=='ca')) = 1d0 ;
    staq_all(find(chrsld_all=='dp'), find(chraq_all=='mg')) = 1d0 ;
    staq_all(find(chrsld_all=='dp'), find(chraq_all=='si')) = 2d0 ;
    % Hedenbergite (FeCaSi2O6)
    staq_all(find(chrsld_all=='hb'), find(chraq_all=='ca')) = 1d0 ;
    staq_all(find(chrsld_all=='hb'), find(chraq_all=='fe2')) = 1d0 ;
    staq_all(find(chrsld_all=='hb'), find(chraq_all=='si')) = 2d0 ;
    % Clinopyroxene (FexMg(1-x)CaSi2O6)
    staq_all(find(chrsld_all=='cpx'), find(chraq_all=='ca')) = 1d0 ;
    staq_all(find(chrsld_all=='cpx'), find(chraq_all=='fe2')) = fr_hb_cpx;
    staq_all(find(chrsld_all=='cpx'), find(chraq_all=='mg')) = 1d0 - fr_hb_cpx;
    staq_all(find(chrsld_all=='cpx'), find(chraq_all=='si')) = 2d0 ;
    % Enstatite (MgSiO3)
    staq_all(find(chrsld_all=='en'), find(chraq_all=='mg')) = 1d0 ;
    staq_all(find(chrsld_all=='en'), find(chraq_all=='si')) = 1d0 ;
    % Ferrosilite (FeSiO3)
    staq_all(find(chrsld_all=='fer'), find(chraq_all=='fe2')) = 1d0 ;
    staq_all(find(chrsld_all=='fer'), find(chraq_all=='si')) = 1d0 ;
    % Orthopyroxene (FexMg(1-x)SiO3)
    staq_all(find(chrsld_all=='opx'), find(chraq_all=='fe2')) = fr_fer_opx;
    staq_all(find(chrsld_all=='opx'), find(chraq_all=='mg')) = 1d0 - fr_fer_opx;
    staq_all(find(chrsld_all=='opx'), find(chraq_all=='si')) = 1d0 ;
    % Augite (Fe(xy+x)Mg(y-xy+1-x)Ca(1-y)Si2O6); x=fr_fer_agt ; y=fr_opx_agt 
    staq_all(find(chrsld_all=='agt'), find(chraq_all=='fe2')) = fr_fer_agt* (1d0 + fr_opx_agt);
    staq_all(find(chrsld_all=='agt'), find(chraq_all=='mg')) = (1d0 - fr_fer_agt )*(fr_opx_agt + 1d0);
    staq_all(find(chrsld_all=='agt'), find(chraq_all=='ca')) = 1d0 - fr_opx_agt;
    staq_all(find(chrsld_all=='agt'), find(chraq_all=='si')) = 2d0 ;
    % Tremolite (Ca2Mg5(Si8O22)(OH)2)
    staq_all(find(chrsld_all=='tm'), find(chraq_all=='ca')) = 2d0 ;
    staq_all(find(chrsld_all=='tm'), find(chraq_all=='mg')) = 5d0 ;
    staq_all(find(chrsld_all=='tm'), find(chraq_all=='si')) = 8d0 ;
    % Anthophyllite (Mg2Mg5(Si8O22)(OH)2)
    staq_all(find(chrsld_all=='antp'), find(chraq_all=='mg')) = 7d0 ;
    staq_all(find(chrsld_all=='antp'), find(chraq_all=='si')) = 8d0 ;
    % Muscovite; KAl2(AlSi3O10)(OH)2
    staq_all(find(chrsld_all=='mscv'), find(chraq_all=='k')) = 1d0 ;
    staq_all(find(chrsld_all=='mscv'), find(chraq_all=='al')) = 3d0 ;
    staq_all(find(chrsld_all=='mscv'), find(chraq_all=='si')) = 3d0 ;
    % Phlogopite; KMg3(AlSi3O10)(OH)2
    staq_all(find(chrsld_all=='plgp'), find(chraq_all=='k')) = 1d0 ;
    staq_all(find(chrsld_all=='plgp'), find(chraq_all=='mg')) = 3d0 ;
    staq_all(find(chrsld_all=='plgp'), find(chraq_all=='al')) = 1d0 ;
    staq_all(find(chrsld_all=='plgp'), find(chraq_all=='si')) = 3d0 ;
    % Amorphous silica; SiO2
    staq_all(find(chrsld_all=='amsi'), find(chraq_all=='si')) = 1d0 ;
    % Quartz; SiO2
    staq_all(find(chrsld_all=='qtz'), find(chraq_all=='si')) = 1d0 ;
    % Aragonite (CaCO3)
    staq_all(find(chrsld_all=='arg'), find(chraq_all=='ca')) = 1d0 ;
    stgas_all(find(chrsld_all=='arg'), find(chrgas_all=='pco2')) = 1d0 ;
    % Dolomite (CaMg(CO3)2)
    staq_all(find(chrsld_all=='dlm'), find(chraq_all=='ca')) = 1d0 ;
    staq_all(find(chrsld_all=='dlm'), find(chraq_all=='mg')) = 1d0 ;
    stgas_all(find(chrsld_all=='dlm'), find(chrgas_all=='pco2')) = 2d0 ;
    % Gypsum; CaSO4*2H2O
    staq_all(find(chrsld_all=='gps'), find(chraq_all=='ca')) = 1d0 ;
    staq_all(find(chrsld_all=='gps'), find(chraq_all=='so4')) = 1d0 ;
    % OMs; CH2O
    stgas_all(find(chrsld_all=='g1'), find(chrgas_all=='pco2')) = 1d0 ;
    stgas_all(find(chrsld_all=='g1'), find(chrgas_all=='po2')) = -1d0 ;
    stgas_all(find(chrsld_all=='g1'), find(chrgas_all=='pnh3')) = n2c_g1;

    stgas_all(find(chrsld_all=='g2'), find(chrgas_all=='pco2')) = 1d0 ;
    stgas_all(find(chrsld_all=='g2'), find(chrgas_all=='po2')) = -1d0 ;
    stgas_all(find(chrsld_all=='g2'), find(chrgas_all=='pnh3')) = n2c_g2;

    stgas_all(find(chrsld_all=='g3'), find(chrgas_all=='pco2')) = 1d0 ;
    stgas_all(find(chrsld_all=='g3'), find(chrgas_all=='po2')) = -1d0 ;
    stgas_all(find(chrsld_all=='g3'), find(chrgas_all=='pnh3')) = n2c_g3;
    % the above need to be modified to enable anoxic degradation 

    staq = zeros(nsp_sld,nsp_aq,'double');
    stgas = zeros(nsp_sld,nsp_gas,'double');

    for isps = 1: nsp_sld
        if (any(chrsld_all == chrsld(isps)))  
            for ispa = 1: nsp_aq 
                if (any(chraq_all == chraq(ispa)))  
                    staq(isps,ispa) = staq_all(find(chrsld_all==chrsld(isps)), find(chraq_all==chraq(ispa)));
                end 
            end 
            for ispg = 1: nsp_gas 
                if (any(chrgas_all == chrgas(ispg)))  
                    stgas(isps,ispg) = stgas_all(find(chrsld_all==chrsld(isps)), find(chrgas_all==chrgas(ispg)));
                end 
            end 
        end 
    end 
    
    % external reactions
    staq_ext_all = zeros(nrxn_ext_all,nsp_aq_all,'double');
    stgas_ext_all = zeros(nrxn_ext_all,nsp_gas_all,'double');
    stsld_ext_all = zeros(nrxn_ext_all,nsp_sld_all,'double');
    % respiration 
    stgas_ext_all(find(chrrxn_ext_all=='resp'), find(chrgas_all=='pco2')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='resp'), find(chrgas_all=='po2')) = -1d0;
    % fe2 oxidation 
    staq_ext_all(find(chrrxn_ext_all=='fe2o2'), find(chraq_all=='fe2')) = -1d0;
    staq_ext_all(find(chrrxn_ext_all=='fe2o2'), find(chraq_all=='fe3')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='fe2o2'), find(chrgas_all=='po2')) = -1d0/4d0;
    % SOC assimilation by microbes 
    stsld_ext_all(find(chrrxn_ext_all=='omomb'), find(chrsld_all=='om')) = -1d0;
    stsld_ext_all(find(chrrxn_ext_all=='omomb'), find(chrsld_all=='omb')) = 0.31d0;
    stgas_ext_all(find(chrrxn_ext_all=='omomb'), find(chrgas_all=='pco2')) = 0.69d0;
    stgas_ext_all(find(chrrxn_ext_all=='omomb'), find(chrgas_all=='po2')) = -0.69d0;
    % turnover of microbes 
    stsld_ext_all(find(chrrxn_ext_all=='ombto'), find(chrsld_all=='om')) = 1d0;
    stsld_ext_all(find(chrrxn_ext_all=='ombto'), find(chrsld_all=='omb')) = -1d0;
    % pyrite oxidation by fe3
    stsld_ext_all(find(chrrxn_ext_all=='pyfe3'), find(chrsld_all=='py')) = -1d0;
    staq_ext_all(find(chrrxn_ext_all=='pyfe3'), find(chraq_all=='fe3')) = -14d0;
    staq_ext_all(find(chrrxn_ext_all=='pyfe3'), find(chraq_all=='fe2')) = 15d0;
    staq_ext_all(find(chrrxn_ext_all=='pyfe3'), find(chraq_all=='so4')) = 2d0;
    % ammonia oxidation by O2 (NH4+ + 2O2 -> NO3- + H2O + 2 H+) 
    staq_ext_all(find(chrrxn_ext_all=='amo2o'), find(chraq_all=='no3')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='amo2o'), find(chrgas_all=='pnh3')) = -1d0;
    stgas_ext_all(find(chrrxn_ext_all=='amo2o'), find(chrgas_all=='po2')) = -2d0;
    % overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
    staq_ext_all(find(chrrxn_ext_all=='g2n0'), find(chraq_all=='no3')) = -4d0/5d0; % values relative to CH2O 
    stsld_ext_all(find(chrrxn_ext_all=='g2n0'), find(chrsld_all=='g2')) = -1d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n0'), find(chrgas_all=='pco2')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n0'), find(chrgas_all=='pnh3')) = n2c_g2;
    % stgas_ext_all(find(chrrxn_ext_all=='g2n0'), find(chrgas_all=='pn2')) = 2d0/5d0; % should be added after enabling pn2 
    % first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O) 
    staq_ext_all(find(chrrxn_ext_all=='g2n21'), find(chraq_all=='no3')) = -1d0;  % values relative to CH2O
    stsld_ext_all(find(chrrxn_ext_all=='g2n21'), find(chrsld_all=='g2')) = -1d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n21'), find(chrgas_all=='pco2')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n21'), find(chrgas_all=='pn2o')) = 0.5d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n21'), find(chrgas_all=='pnh3')) = n2c_g2;
    % 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O) 
    stsld_ext_all(find(chrrxn_ext_all=='g2n22'), find(chrsld_all=='g2')) = -1d0; % values relative to CH2O
    stgas_ext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='pco2')) = 1d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='pn2o')) = -2d0;
    stgas_ext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='pnh3')) = n2c_g2;
    % stgas_ext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='pn2')) = 2d0; % should be added after enabling pn2 
    
    % define 1 when a reaction is sensitive to a speces 
    stgas_dext_all = zeros(nrxn_ext_all,nsp_gas_all,'double');
    staq_dext_all = zeros(nrxn_ext_all,nsp_aq_all,'double');
    stsld_dext_all = zeros(nrxn_ext_all,nsp_sld_all,'double');
    % respiration 
    stgas_dext_all(find(chrrxn_ext_all=='resp'), find(chrgas_all=='po2')) = 1d0;
    % fe2 oxidation 
    stgas_dext_all(find(chrrxn_ext_all=='fe2o2'), find(chrgas_all=='po2')) = 1d0;
    stgas_dext_all(find(chrrxn_ext_all=='fe2o2'), find(chrgas_all=='pco2')) = 1d0;
    staq_dext_all(find(chrrxn_ext_all=='fe2o2'), find(chraq_all=='fe2')) = 1d0;
    % SOC assimilation by microbes 
    stsld_dext_all(find(chrrxn_ext_all=='omomb'), find(chrsld_all=='om')) = 1d0;
    stsld_dext_all(find(chrrxn_ext_all=='omomb'), find(chrsld_all=='omb')) = 1d0;
    % turnover of microbes 
    stsld_dext_all(find(chrrxn_ext_all=='ombto'), find(chrsld_all=='omb')) = 1d0;
    % pyrite oxidation by fe3
    stsld_dext_all(find(chrrxn_ext_all=='pyfe3'), find(chrsld_all=='py')) = 1d0;
    staq_dext_all(find(chrrxn_ext_all=='pyfe3'), find(chraq_all=='fe2')) = 1d0;
    staq_dext_all(find(chrrxn_ext_all=='pyfe3'), find(chraq_all=='fe3')) = 1d0;
    % ammonia oxidation by O2 (NH4+ + 2O2 -> NO3- + H2O + 2 H+) 
    stgas_dext_all(find(chrrxn_ext_all=='amo2o'), find(chrgas_all=='po2')) = 1d0;
    stgas_dext_all(find(chrrxn_ext_all=='amo2o'), find(chrgas_all=='pnh3')) = 1d0;
    % overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
    staq_dext_all(find(chrrxn_ext_all=='g2n0'), find(chraq_all=='no3')) = 1d0;
    stgas_dext_all(find(chrrxn_ext_all=='g2n0'), find(chrgas_all=='po2')) = 1d0;
    stsld_dext_all(find(chrrxn_ext_all=='g2n0'), find(chrsld_all=='g2')) = 1d0;
    % first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O) 
    staq_dext_all(find(chrrxn_ext_all=='g2n21'), find(chraq_all=='no3')) = 1d0;
    stgas_dext_all(find(chrrxn_ext_all=='g2n21'), find(chrgas_all=='po2')) = 1d0;
    stsld_dext_all(find(chrrxn_ext_all=='g2n21'), find(chrsld_all=='g2')) = 1d0;
    % 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O) 
    staq_dext_all(find(chrrxn_ext_all=='g2n22'), find(chraq_all=='no3')) = 1d0;
    % stgas_dext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='po2')) = 1d0;
    stgas_dext_all(find(chrrxn_ext_all=='g2n22'), find(chrgas_all=='pn2o')) = 1d0;
    stsld_dext_all(find(chrrxn_ext_all=='g2n22'), find(chrsld_all=='g2')) = 1d0;

    staq_ext = zeros(nrxn_ext,nsp_aq,'double');
    staq_dext = zeros(nrxn_ext,nsp_aq,'double');
    stgas_ext = zeros(nrxn_ext,nsp_gas,'double');
    stgas_dext = zeros(nrxn_ext,nsp_gas,'double');
    stsld_ext = zeros(nrxn_ext,nsp_sld,'double');
    stsld_dext = zeros(nrxn_ext,nsp_sld,'double');

    for irxn = 1: nrxn_ext
        if (any(chrrxn_ext_all == chrrxn_ext(irxn)))  
            for ispa = 1: nsp_aq 
                if (any(chraq_all == chraq(ispa)))  
                    staq_ext(irxn,ispa) = staq_ext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chraq_all==chraq(ispa)));
                    staq_dext(irxn,ispa) = staq_dext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chraq_all==chraq(ispa)));
                end 
            end 
            for ispg = 1: nsp_gas 
                if (any(chrgas_all == chrgas(ispg)))  
                    stgas_ext(irxn,ispg) = stgas_ext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chrgas_all==chrgas(ispg)));
                    stgas_dext(irxn,ispg) = stgas_dext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chrgas_all==chrgas(ispg)));
                end 
            end 
            for isps = 1: nsp_sld 
                if (any(chrsld_all == chrsld(isps)))  
                    stsld_ext(irxn,isps) = stsld_ext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chrsld_all==chrsld(isps)));
                    stsld_dext(irxn,isps) = stsld_dext_all(find(chrrxn_ext_all==chrrxn_ext(irxn)),find(chrsld_all==chrsld(isps)));
                end 
            end 
        end 
    end 
    
    def_dust = 0d0;
    [rfrc_sld_all] = get_dust(nsp_sld_all,chrsld_all,def_dust);
    rfrc_sld_all = rfrc_sld_all./mwt_all';
    
    def_OM_frc = 0d0;
    [rfrc_sld_plant_all] = get_OM_rain(nsp_sld_all,chrsld_all,def_OM_frc);
    
    rfrc_sld = zeros(nsp_sld,1,'double');
    rfrc_sld_plant = zeros(nsp_sld,1,'double');
    for isps = 1: nsp_sld 
        rfrc_sld(isps) = rfrc_sld_all(find(chrsld_all==chrsld(isps)));
        rfrc_sld_plant(isps) = rfrc_sld_plant_all(find(chrsld_all==chrsld(isps)));
    end
    
    [iwtype,imixtype,poroiter_in,display,display_lim_in,read_data,incld_rough ...
    ,al_inhibit,timestep_fixed,method_precalc,regular_grid,sld_enforce ...% inout
    ,poroevol,surfevol1,surfevol2,do_psd,lim_minsld_in,do_psd_full,season ...% inout
    ] = get_switches;


    no_biot = false;
    biot_turbo2 = false;
    biot_fick = false;
    biot_labs = false;
    biot_till = false;

    switch (imixtype)
        case(imixtype_nobio)
            no_biot = true;
        case(imixtype_fick)
            biot_fick = true;
        case(imixtype_turbo2)
            biot_turbo2 = true;
        case(imixtype_till)
            biot_till = true;
        case(imixtype_labs)
            biot_labs = true;
        otherwise 
            warning( '***| chosen number is not available for mixing styles (choose between 0 to 4)')
            warning( '***| thus choose default |---- > no mixing')
            no_biot = true;
    end

    fprintf('no_biot = %s\nFickian mix. = %s\nturbo2 = %s\ntilling = %s\nlabs-mix. = %s\n' ...
        ,string(no_biot),string(biot_fick),string(biot_turbo2),string(biot_till),string(biot_labs))
        
    
    switch(iwtype)
        case(iwtype_cnst)
            fprintf('const w : %d\n',iwtype)
        case(iwtype_flex)
            fprintf('w flex (cnst porosity profile) : %d\n',iwtype)
        case(iwtype_pwcnst)
            fprintf('w x porosity = cnst : %d\n',iwtype)
        case(iwtype_spwcnst)
            fprintf('w x (1 - porosity) = cnst : %d\n',iwtype)
        otherwise 
            warning( '***| chosen number is not available for advection styles (choose between 0 to 3)')
            warning( '***| thus choose default |---- > cnst w')
            iwtype = iwtype_cnst
    end


    if (poroiter_in)  
        fprintf('porosity iteration is ON\n')
    else 
        fprintf('porosity iteration is OFF\n')
    end 

    if (lim_minsld_in)  
        fprintf('limiting lowest mineral conc. is ON\n')
    else 
        fprintf('limiting lowest mineral conc. is OFF\n')
    end 

    if (do_psd_full); do_psd = true;end

    if (display_lim_in); display_lim = true; end

    if (sld_enforce); nsp3 = nsp_aq + nsp_gas; end % excluding solid phases

    % kinetic formulation type
    precstyle = strings(nsp_sld,1);
    precstyle(:) = 'def';
    % precstyle(:) = 'full';
    % precstyle(:) = 'full_lim';
    % precstyle(:) = 'seed ';
    % precstyle(:) = '2/3';
    % precstyle(:) = 'psd_full';

    for isps = 1: nsp_sld
        switch chrsld(isps)
            case{'g1','g2','g3'}
                precstyle(isps) = 'decay';
            otherwise
                precstyle(isps) = 'def';
                % precstyle(isps) = '2/3';
                % precstyle(isps) = '2/3noporo';
                % precstyle(isps) = 'psd_full';
        end
    end 

    while (rectime_prof(nrec_prof)>ttot) 
        rectime_prof = rectime_prof/10d0;
    end
    while (rectime_prof(nrec_prof)<ttot) 
        rectime_prof = rectime_prof*10d0;
    end

    rectime_flx = zeros(nrec_flx,1,'double');
    for irec_flx = 1:20
        rectime_flx(irec_flx) = irec_flx/20d0;
    end
    for irec_flx = 21:38
        rectime_flx(irec_flx) = rectime_flx(20) + (irec_flx-20)/20d0*10d0;
    end
    for irec_flx = 39:60
        rectime_flx(irec_flx) = rectime_flx(38) + (irec_flx-38)/20d0*100d0;
    end

    while (rectime_flx(nrec_flx)>ttot) 
        rectime_flx = rectime_flx/10d0;
    end
    while (rectime_flx(nrec_flx)<ttot) 
        rectime_flx = rectime_flx*10d0;
    end
    
    workdir=strcat('../pyweath_output/matlab/',sim_name,'/') 
    flxdir= strcat(workdir,'flx' )  
    profdir=strcat(workdir,'prof' )  
    
    % isldflx = zeros(nsp_sld,1);
    % for isps = 1: nsp_sld 
        % isldflx(isps) = idust + isps;
    % end 
        
    % iaqflx = zeros(nsp_aq,1);
    % for ispa = 1: nsp_aq 
        % iaqflx(ispa) = idust + nsp_sld  + ispa;
    % end 

    % igasflx = zeros(nsp_gas,1);
    % for ispg = 1: nsp_gas
        % igasflx(ispg) = idust + nsp_sld + nsp_aq + ispg;
    % end 

    % ico2flx = zeros(6,1);
    % for ico2 = 1: 6
        % ico2flx(ico2) = idust + nsp_sld + nsp_aq + nsp_gas + ico2;
    % end 
    
    if ~exist(flxdir, 'dir'); mkdir(flxdir); end
    if ~exist(profdir, 'dir'); mkdir(profdir); end

    fmt=[repmat('%s \t',1,nflx+1) '\n'];
    isldflx = strings(nsp_sld,1);
    for isps = 1:nsp_sld
        isldflx(isps) = sprintf('%s/flx_sld-%s.txt',flxdir,chrsld(isps));
        fid = fopen(isldflx(isps),'w');
        fprintf(fid,fmt, 'time',chrflx(1:nflx));
        fclose(fid);
    end 

    iaqflx = strings(nsp_aq,1);
    for ispa = 1:nsp_aq
        iaqflx(ispa) = sprintf('%s/flx_aq-%s.txt',flxdir,chraq(ispa));
        fid = fopen(iaqflx(ispa),'w');
        fprintf(fid,fmt, 'time',chrflx(1:nflx));
        fclose(fid);
    end 

    igasflx = strings(nsp_gas,1);
    for ispg = 1:nsp_gas
        igasflx(ispg) = sprintf('%s/flx_gas-%s.txt',flxdir,chrgas(ispg));
        fid = fopen(igasflx(ispg),'w');
        fprintf(fid,fmt, 'time',chrflx(1:nflx));
        fclose(fid);
    end 

    ico2flx = strings(6,1);
    for ico2 = 1:6
        ico2flx(ico2) = sprintf('%s/flx_co2sp-%s.txt',flxdir,chrco2sp(ico2));
        fid = fopen(ico2flx(ico2),'w');
        fprintf(fid,fmt, 'time',chrflx(1:nflx));
        fclose(fid);
    end 

    idust = sprintf('%s/dust.txt',flxdir);
    fid = fopen(idust,'w');
    fprintf(fid,'%s\t%s\n', 'time', 'dust(relative_to_average)');
    fclose(fid);

    climate = zeros(3,1,'logical');
    climate(:) = false;
    if (season); climate(:) = true; end


    idust = sprintf('%s/climate.txt',flxdir);
    fid = fopen(idust,'w');
    fprintf(fid, '%s\t%s\t%s\t%s\n', 'time', 'T(oC)', 'q(m/yr)', 'Wet(-)');
    fclose(fid);
    
    dust_norm_prev = 0d0; dust_norm = 0d0;

    clim_file = string({'T_temp.in','q_temp.in','Wet_temp.in'});

    dct = zeros(3,1,'double');
    ctau = zeros(3,1,'double');

    for iclim = 1:3
        if (climate(iclim))  
            [nclim(iclim)] = get_clim_num(clim_file(iclim));
            switch (iclim) 
                case(1)
                    clim_T = zeros(2,nclim(iclim),'double');
                    idust = sprintf('%s/%s',workdir,clim_file(iclim));
                    fid = fopen(idust, 'r');
                    tline = fgetl(fid);
                    cnt = 0;
                    while ischar(tline)
                        a = strsplit(tline);
                        tline = fgetl(fid);
                        if cnt ==0
                            cnt = cnt + 1;
                            continue
                        end
                        clim_T(1,cnt) = (str2num(char(a(1))));
                        clim_T(2,cnt) = (str2num(char(a(2))));
                        cnt = cnt + 1;
                    end
                    fclose(fid);
                    for ict = 1: nclim(iclim)
                        fprintf('%4.3e %4.3e\n',clim_T(:,ict));
                    end 
                    dct(iclim) = clim_T(1,2) - clim_T(1,1);
                    ctau(iclim) = clim_T(1,nclim(iclim)) + dct(iclim);
                case(2)
                    clim_q = zeros(2,nclim(iclim),'double')
                    idust = sprintf('%s/%s',workdir,clim_file(iclim));
                    fid = fopen(idust,'r');
                    tline = fgetl(fid);
                    cnt = 0;
                    while ischar(tline)
                        a = strsplit(tline);
                        tline = fgetl(fid);
                        if cnt ==0
                            cnt = cnt + 1;
                            continue
                        end
                        clim_q(1,cnt) = (str2num(char(a(1))));
                        clim_q(2,cnt) = (str2num(char(a(2))));
                        cnt = cnt + 1;
                    end
                    fclose(fid);
                    % converting mm/month to m/yr 
                    clim_q(2,:) = clim_q(2,:)*12d0/1d3;
                    for ict = 1: nclim(iclim)
                        fprintf('%4.3e %4.3e\n',clim_q(:,ict));
                    end 
                    dct(iclim) = clim_T(1,2) - clim_T(1,1);
                    ctau(iclim) = clim_T(1,nclim(iclim)) + dct(iclim);
                case(3)
                    clim_sat = zeros(2,nclim(iclim),'double')
                    idust = sprintf('%s/%s',workdir,clim_file(iclim));
                    fid = fopen(idust,'r');
                    tline = fgetl(fid);
                    cnt = 0;
                    while ischar(tline)
                        a = strsplit(tline);
                        tline = fgetl(fid);
                        if cnt ==0
                            cnt = cnt + 1;
                            continue
                        end
                        clim_sat(1,cnt) = (str2num(char(a(1))));
                        clim_sat(2,cnt) = (str2num(char(a(2))));
                        cnt = cnt + 1;
                    end
                    fclose(fid);
                    % converting mm/m to m/m  
                    clim_sat(2,:) = clim_sat(2,:)*1d0/1d3;
                    for ict = 1: nclim(iclim)
                        fprintf('%4.3e %4.3e\n',clim_q(:,ict));
                    end 
                    dct(iclim) = clim_T(1,2) - clim_T(1,1);
                    ctau(iclim) = clim_T(1,nclim(iclim)) + dct(iclim);
                otherwise
                    warning( 'error in obtaining climate')
            end 
        end 
    end 

    %%%  MAKING GRID %%%%%%%%%%%%%%%%% 
    beta = 1.00000000005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
    beta = 1.00005d0;  % a parameter to make a grid; closer to 1, grid space is more concentrated around the sediment-water interface (SWI)
    [dz,z] = makegrid(beta,nz,ztot,regular_grid);
    
    sat = zeros(nz,1,'double');
    sat(:) = min(1.0d0,(1d0-satup)*z(:)/zsat + satup);
    
    
    [nsld_sa] = get_sa_num;
    
    [hrii,chrsld_sa] =  get_sa(nsp_sld,chrsld,p80,nsld_sa);

    hri = zeros(nsp_sld,nz,'double');
    for isps = 1: nsp_sld
        hri(isps,:) = 1d0/hrii(isps);
    end
    
    rough = ones(nsp_sld,nz,'double');
    if (incld_rough)  
        % rough = 10d0**(3.3d0)*p80**0.33d0 ! from Navarre-Sitchler and Brantley (2007)
        for isps=1:nsp_sld
            rough(isps,:) = rough_c0*(1d0./hri(isps,:)).^rough_c1; % from Navarre-Sitchler and Brantley (2007)
        end 
    end

    hr = zeros(nsp_sld,nz,'double');
    for isps=1:nsp_sld
        hr(isps,:) = hri(isps,:).*rough(isps,:);
    end 
    
    hrb = zeros(nz,1,'double');
    ssab = zeros(nz,1,'double');

    v = zeros(nz,1,'double');
    w = zeros(nz,1,'double');
    poro = zeros(nz,1,'double');
    torg = zeros(nz,1,'double');
    tora = zeros(nz,1,'double');

    v(:) = qin/poroi./sat(:);
    poro(:) = poroi;
    torg(:) = poro(:).^(3.4d0-2.0d0).*(1.0d0-sat(:)).^(3.4d0-1.0d0);
    tora(:) = poro(:).^(3.4d0-2.0d0).*(sat(:)).^(3.4d0-1.0d0);

    w_btm = w0;
    w(:) = w_btm;

    % ------------ determine calculation scheme for advection (from IMP code)
    [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
    

    % attempting to do psd 
    if (do_psd)  
        ps = zeros(nps,1,'double');
        dps = zeros(nps,1,'double');
        for ips = 1: nps
            ps(ips) = log10(ps_min) + (ips - 1d0)*(log10(ps_max) - log10(ps_min))/(nps - 1d0);
        end 
        dps(:) = ps(2) - ps(1);
        
        if (do_psd_full)  % do psd for every mienral
            ipsd = sprintf('%s/psd_pr.txt',profdir);
            fmt=['%s\t' repmat('%7.6f \t',1,nps) '%s \n'];
            fid = fopen(ipsd,'w');
            fprintf(fid,fmt,'sldsp\log10(radius)', ps(1:nps), 'time');
            fclose(fid);
            fmt=['%s\t' repmat('%7.6f \t',1,nps+1) ' \n'];
            mpsd_pr = zeros(nsp_sld,nps,'double');
            for isps = 1: nsp_sld
                volsld = msldi(isps)*mv(isps)*1d-6;
                [psd_pr] = calc_psd_pr( ...
                    nps ...% input
                    ,pi,hrii(isps),ps_sigma_std,poroi,volsld,tol ...% input
                    ,ps,dps ...% input
                    ,msldunit ...% input
                    );
                mpsd_pr(isps,:) = psd_pr(:);
                fid = fopen(ipsd,'a');
                fprintf(fid,fmt, chrsld(isps),psd_pr(1:nps), 0d0);
                fclose(fid);
            end 

            % initially particle is distributed as in parent rock 
            psd = zeros(nps,nz,'double');
            mpsd = zeros(nsp_sld,nps,nz,'double');
            ssa = zeros(nsp_sld,nz,'double');
            ssav = zeros(nsp_sld,nz,'double');
            ssv = zeros(nsp_sld,nz,'double');
            for isps=1:nsp_sld
                for iz = 1: nz
                    mpsd(isps,:,iz) = mpsd_pr(isps,:); 
                end 
            
                if (~incld_rough)  
                    for iz=1:nz
                        ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0.*mpsd(isps,:,iz)'.*dps(:) );
                        ssav(isps,iz) = sum( 3d0/(10d0.^ps(:)).*mpsd(isps,:,iz)'.*dps(:) );
                    end 
                else 
                    for iz=1:nz
                        ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0 ...
                            *rough_c0.*(10d0.^ps(:)).^rough_c1.*mpsd(isps,:,iz)'.*dps(:) );
                        ssav(isps,iz) = sum( 3d0./(10d0.^ps(:)) ... 
                            *rough_c0.*(10d0.^ps(:)).^rough_c1.*mpsd(isps,:,iz)'.*dps(:) );
                    end 
                end 
                for iz=1:nz
                    ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*mpsd(isps,:,iz)'.*dps(:));
                end 
            end 
            for iz=1:nz
                hr(:,iz) = ssa(:,iz)./ssv(:,iz)/poro(iz);
            end 
        else % do psd only for bulk 
            volsld = sum(msldi.*mv'*1d-6) + mblki*mvblk*1d-6;
            [psd_pr] = calc_psd_pr( ...
                nps ...% input
                ,pi,p80,ps_sigma_std,poroi,volsld,tol ...% input
                ,ps,dps ...% input
                ,msldunit ...% input
                );
            ipsd = sprintf('%s/psd_pr.txt',profdir);
            fmt=['%s\t' repmat('%7.6f \t',1,nps) '%s \n'];
            fid = fopen(ipsd,'w');
            fprintf(fid,fmt,'depth\log10(radius)', ps(1:nps), 'time');
            fmt=[repmat('%7.6f \t',1,nps+2) ' \n'];
            fprintf(fid,fmt, ztot,psd_pr(1:nps), 0d0);
            fclose(fid);

            % initially particle is distributed as in parent rock 
            psd = zeros(nps,nz,'double');
            mpsd = zeros(nsp_sld,nps,nz,'double');
            ssa = zeros(nsp_sld,nz,'double');
            ssav = zeros(nsp_sld,nz,'double');
            ssv = zeros(nsp_sld,nz,'double');
            
            for iz = 1: nz
                psd(:,iz) = psd_pr(:); 
            end 

            % dM = M * [psd*dps*S(r)] * k *dt 
            % so hr = sum (psd(:)*dps(:)*S(:) ) where S in units m2/m3 and simplest way 1/r 
            % in this case hr = sum(  psd(:)*dps(:)*1d0/(10d0^(-ps(:))) )
            if (~incld_rough)  
                for iz=1:nz
                    ssa(:,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0.*psd(:,iz).*dps(:) );
                    ssav(:,iz) = sum( 3d0/(10d0.^ps(:)).*psd(:,iz).*dps(:) );
                end 
            else 
                for iz=1:nz
                    ssa(:,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0 *rough_c0.*(10d0.^ps(:)).^rough_c1.*psd(:,iz).*dps(:) );
                    ssav(:,iz) = sum( 3d0./(10d0.^ps(:)) *rough_c0.*(10d0.^ps(:)).^rough_c1.*psd(:,iz).*dps(:) );
                end 
            end 
            for iz=1:nz
                ssv(:,iz) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd(:,iz).*dps(:) );
            end 
            % hr = ssa *(1-poro)/poro % converting m2/sld-m3 to m2/pore-m3
            % hr = ssa 
            for iz=1:nz
                hr(:,iz) = ssa(:,iz)/poro(iz); % so that poro * hr * mv * msld becomes porosity independent
            end 
        end 
    else 
        psd = zeros(nps,nz,'double');
        mpsd = zeros(nsp_sld,nps,nz,'double');
        ssa = zeros(nsp_sld,nz,'double');
        ssav = zeros(nsp_sld,nz,'double');
        ssv = zeros(nsp_sld,nz,'double');
        ps = zeros(nps,1,'double');
        dps = zeros(nps,1,'double');
    end 


    minsld = ones(nsp_sld,1,'double')*1d-20;
    % when minimum psd is limited, msld min is limited
    if (do_psd && do_psd_full && psd_lim_min && psd_vol_consv)  
        lim_minsld_in = true;
        for isps=1:nsp_sld
            minsld(isps) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0*psd_th*dps(:)) /( mv(isps)*1d-6 ) %  m3/m3 / (m3 /mol ) =  mol/m3
        end 
    end


    dt = maxdt;

    dt = 1d-20; % for basalt exp?

    maq = zeros(nsp_aq,nz,'double');
    mgas = zeros(nsp_gas,nz,'double');
    msld = zeros(nsp_sld,nz,'double');
    for ispa = 1: nsp_aq
        maq(ispa,:)=maqi(ispa);
    end 
    for ispg = 1: nsp_gas
        mgas(ispg,:)=mgasi(ispg);
    end 
    for isps = 1: nsp_sld
        msld(isps,:) = msldi(isps);
    end 

    rho_grain_z=zeros(nz,1,'double');sldvolfrac=zeros(nz,1,'double');
    mblk = zeros(nz,1,'double');
    mblk(:) = mblki;
    

    omega = zeros(nsp_sld,nz,'double');

    pro = ones(nz,1,'double')*1d-5;
   
    nsld_kinspc = nsld_kinspc_in;
    chrsld_kinspc = strings(nsld_kinspc,1);
    kin_sld_spc = zeros(nsld_kinspc,1,'double');
    chrsld_kinspc(:) = chrsld_kinspc_in(:);
    kin_sld_spc(:) = kin_sld_spc_in(:);


    [...
        ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 ...! output
        ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all ...! output
        ] = coefs_v2( ...
        nz,rg,rg2,tc,sec2yr,tempk_0,pro ...! input
        ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all ...! input
        ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all ...! input
        ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all ...!input
        );


    print_cb = false; 
    print_loc = './ph.txt';

    so4f = zeros(nz,1,'double');
    [ ...
        dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
        ,pro,ph_error,so4f,ph_iter ...% output
        ] = calc_pH_v7_3( ...
        nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
        ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
        ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
        ,print_cb,print_loc,z ...% input 
        ,pro,so4f ...% inout
        ); 
    % fprintf(strcat('showing pH here:\n',[repmat('%7.6e\t',1,5)],'\n') ... 
        % ,-log10(pro(1:nz/5:nz)));
    
    so4fprev = so4f;
    proi = pro(1);
    
    poroprev = poro;
    
    
    %  --------- read -----
    if (read_data)  
        loc_runname_save = strcat('../pyweath_output/',runname_save,'prof');
        if (runname_save == 'self'); loc_runname_save = profdir; end
        
        copyfile( strcat(loc_runname_save,'/prof_sld-save.txt'), strcat(profdir,'/prof_sld-restart.txt') );
        copyfile( strcat(loc_runname_save,'/prof_aq-save.txt'), strcat(profdir,'/prof_aq-restart.txt') );
        copyfile( strcat(loc_runname_save,'/prof_gas-save.txt'), strcat(profdir,'/prof_gas-restart.txt') );
        copyfile( strcat(loc_runname_save,'/bsd-save.txt'), strcat(profdir,'/bsd-restart.txt') );
        copyfile( strcat(loc_runname_save,'/psd-save.txt'), strcat(profdir,'/psd-restart.txt') );
        copyfile( strcat(loc_runname_save,'/sa-save.txt'), strcat(profdir,'/sa-restart.txt') );
            
        dir_empty = '';
        [nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save,nsld_kinspc_save,nsld_sa_save] =...
            get_saved_variables_num( dir_empty,loc_runname_save ); % here 'dir_empty' was used instead of 'workdir' in Fortran code
        
        chraq_save = strings(nsp_aq_save,1); chrsld_save = strings(nsp_sld_save,1);
        chrgas_save = strings(nsp_gas_save,1); chrrxn_ext_save = strings(nrxn_ext_save,1);
        maq_save = zeros(nsp_aq_save,nz,'double'); msld_save = zeros(nsp_sld_save,nz,'double'); mgas_save = zeros(nsp_gas_save,nz,'double');
        chrsld_kinspc_save = strings(nsld_kinspc_save,1); kin_sldspc_save = zeros(nsld_kinspc_save,1,'double');
        chrsld_sa_save = strings(nsld_sa_save,1); hrii_save = zeros(nsld_sa_save,1,'double');
        hr_save = zeros(nsp_sld_save,nz,'double'); mpsd_save = zeros(nsp_sld_save,nps,nz,'double');
        
        [ ... 
            chraq_save,chrgas_save,chrsld_save,chrrxn_ext_save,chrsld_kinspc_save,kin_sldspc_save ...% output
            ,chrsld_sa_save,hrii_save ...% output 
            ] = get_saved_variables( ...
            dir_empty,loc_runname_save ...% input  % here 'dir_empty' was used instead of 'workdir' in Fortran code
            ,nsp_aq_save,nsp_sld_save,nsp_gas_save,nrxn_ext_save,nsld_kinspc_save,nsld_sa_save ...% input
            );
        
        isldprof = fopen(strcat(profdir,'/prof_sld-restart.txt'),'r');
        iaqprof = fopen(strcat(profdir,'/prof_aq-restart.txt'),'r');
        igasprof = fopen(strcat(profdir,'/prof_gas-restart.txt'),'r');
        ibsd = fopen(strcat(profdir,'/bsd-restart.txt'),'r');
        ipsd = fopen(strcat(profdir,'/psd-restart.txt'),'r');
        isa = fopen(strcat(profdir,'/sa-restart.txt'),'r');
        
        tline_1 = fgetl(isldprof);
        tline_2 = fgetl(iaqprof);
        tline_3 = fgetl(igasprof);
        tline_4 = fgetl(ibsd);
        tline_5 = fgetl(ipsd);
        tline_6 = fgetl(isa);
        cnt = 0;
        while ischar(tline_1)
            a = strsplit(tline_1);
            b = strsplit(tline_2);
            c = strsplit(tline_3);
            d = strsplit(tline_4);
            e = strsplit(tline_5);
            f = strsplit(tline_6);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            msld_save(1:nsp_sld_save,cnt) = str2num(char(a(2:nsp_sld_save+1)));
            maq_save(1:nsp_aq_save,cnt) = str2num(char(b(2:nsp_aq_save+1)));
            mgas_save(1:nsp_gas_save,cnt) = str2num(char(c(2:nsp_gas_save+1)));
            poro(iz) = str2num(char(d(2))); sat(iz)= str2num(char(d(3))); v(iz)= str2num(char(d(4))); 
            hrb(iz)= str2num(char(d(5))); w(iz)= str2num(char(d(6))); sldvolfrac(iz)= str2num(char(d(7))); 
            rho_grain_z(iz)= str2num(char(d(8))); mblk(iz)= str2num(char(d(9))); 
            psd(1:nps,cnt) = str2num(char(e(2:nps+1)));
            hr_save(1:nsp_sld_save,cnt) = str2num(char(f(2:nsp_sld_save+1)));
            
            ucvsld1 = 1d0;
            if (msldunit == 'blk'); ucvsld1 = 1d0 - poro(iz); end
            
            mblk(iz) = mblk(iz)/ ( mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6) );
            
            cnt = cnt + 1;
                
        end
        
        fclose(isldprof);
        fclose(iaqprof);
        fclose(igasprof);
        fclose(ibsd);
        fclose(ipsd);
        fclose(isa);

        pro(:) = 10d0.^(-pro(:)); % read data is -log10 (pro)
        
        torg(:) = poro(:).^(3.4d0-2.0d0).*(1.0d0-sat(:)).^(3.4d0-1.0d0);
        tora(:) = poro(:).^(3.4d0-2.0d0).*(sat(:)).^(3.4d0-1.0d0);
        
        for isps = 1:nsp_sld_save
            if (any(chrsld == chrsld_save(isps)))  
                msld(find(chrsld==chrsld_save(isps)),:) = msld_save(isps,:);
            elseif (any(chrsld_cnst == chrsld_save(isps))) 
                msldc(find(chrsld_cnst==chrsld_save(isps)),:) = msld_save(isps,:);
            else 
                error ('error in re-assignment of sld conc.');
            end 
        end 
        
        for ispa = 1:nsp_aq_save
            if (any(chraq == chraq_save(ispa)))  
                maq(find(chraq==chraq_save(ispa)),:) = maq_save(ispa,:);
            elseif (any(chraq_cnst == chraq_save(ispa))) 
                maqc(find(chraq_cnst==chraq_save(ispa)),:) = maq_save(ispa,:);
            else 
                error('error in re-assignment of aq conc.');
            end 
        end 
        
        for ispg = 1:nsp_gas_save
            if (any(chrgas == chrgas_save(ispg)))  
                mgas(find(chrgas==chrgas_save(ispg)),:) = mgas_save(ispg,:);
            elseif (any(chrgas_cnst == chrgas_save(ispg))) 
                mgasc(find(chrgas_cnst==chrgas_save(ispg)),:) = mgas_save(ispg,:);
            else 
                error('error in re-assignment of gas conc.');
            end 
        end 
        
        % counting sld species whose values are to be specificed in kinspc.save and not so yet when reading from kinspc.in
        nsld_kinspc_add = 0;
        for isps_kinspc = 1:nsld_kinspc_save
            if (any(chrsld_kinspc_in == chrsld_kinspc_save(isps_kinspc)))  % already specified 
                continue
            else 
                nsld_kinspc_add = nsld_kinspc_add + 1;
            end 
        end 
        
        if (nsld_kinspc_add > 0)  
            % re-define sld species number whose rate const. is specified 
            nsld_kinspc = nsld_kinspc + nsld_kinspc_add
            % allocate 
            kin_sld_spc = zeros(nsld_kinspc,1,'double'); chrsld_kinspc = strings(nsld_kinspc,1);
            % saving already specified consts. 
            chrsld_kinspc(1:nsld_kinspc_in) = chrsld_kinspc_in(1:nsld_kinspc_in);
            kin_sld_spc(1:nsld_kinspc_in) = kin_sld_spc_in(1:nsld_kinspc_in);
            % adding previously specified rate const. 
            nsld_kinspc_add = 0
            for isps_kinspc = 1:nsld_kinspc_save
                if (any(chrsld_kinspc_in == chrsld_kinspc_save(isps_kinspc)))  
                    continue
                else 
                    nsld_kinspc_add = nsld_kinspc_add + 1;
                    chrsld_kinspc(nsld_kinspc_in + nsld_kinspc_add) = chrsld_kinspc_save(isps_kinspc);
                    kin_sld_spc(nsld_kinspc_in + nsld_kinspc_add) = kin_sldspc_save(isps_kinspc);
                end 
            end 
        end 
        
        % overloading sa if saved 
        
        for isps = 1:nsp_sld_save
            if (any(chrsld == chrsld_save(isps)))  
                hr(find(chrsld==chrsld_save(isps)),:) = hr_save(isps,:);
            end 
        end     
        
        if (nsld_sa_save > 0)  
            for isps_sa =1: nsld_sa_save 
                if (any(chrsld_sa == chrsld_sa_save(isps_sa)))  % if SA of a sld sp. is already specified, save data does not overload 
                    continue 
                else  % if SA is not specified some species but was specified in the previous run, saved data is loaded 
                    if (any(chrsld == chrsld_sa_save(isps_sa))) 
                        hrii(find(chrsld==chrsld_sa_save(isps_sa))) = hrii_save(isps_sa);
                    end 
                end 
            end 
        end 
        % updating SA parameters 
        for isps = 1: nsp_sld
            hri(isps,:) = 1d0/hrii(isps);
        end

        rough(:,:) = 1d0;
        if (incld_rough)  
            % rough = 10d0^(3.3d0)*p80^0.33d0 % from Navarre-Sitchler and Brantley (2007)
            for isps=1:nsp_sld
                rough(isps,:) = rough_c0*(1d0./hri(isps,:)).^rough_c1; % from Navarre-Sitchler and Brantley (2007)
            end 
        end 
        
        if (do_psd_full) 
            % updating parentrock psd if hrii has been loaded from a previous run 
            if (nsld_sa_save > 0) 
                ipsd = sprintf('%s/psd_pr.txt',profdir);
                fid = fopen(ipsd,'w');
                fmt=['%s\t' repmat('%7.6f \t',1,nps) '%s \n'];
                fprintf(fid,fmt,'sldsp\log10(radius)', ps(1:nps), 'time');
                for isps = 1: nsp_sld
                    volsld = msldi(isps)*mv(isps)*1d-6;
                    [psd_pr] = calc_psd_pr( ...
                        nps ...% input
                        ,pi,hrii(isps),ps_sigma_std,poroi,volsld,tol ...% input
                        ,ps,dps ...% input
                        ,msldunit ...% input
                        )
                    mpsd_pr(isps,:) = psd_pr(:)';
                    fprintf(fid,fmt, chrsld(isps),psd_pr(1:nps), 0d0);
                end 
                fclose(fid);
                for isps=1:nsp_sld
                    for iz = 1: nz
                        mpsd(isps,:,iz) = mpsd_pr(isps,:); 
                    end 
                end 
            end 
            
            for isps = 1: nsp_sld_save % loading psds 
            
                copyfile( strcat(loc_runname_save,'/psd_',chrsld_save(isps),'-save.txt') ...
                    , strcat(profdir,'/psd_',chrsld_save(isps),'-restart.txt') );
                ipsd = strcat(profdir,'/psd_',chrsld_save(isps),'-restart.txt');
                fid = fopen(ipsd,'r');
                
                tline = fgetl(fid);
                cnt = 0;
                while ischar(tline)
                    a = strsplit(tline);
                    
                    if cnt ==0
                        cnt = cnt + 1;
                        continue
                    end
                    psd(1:nps,cnt) = str2num(char(a(2:nps+1)));
                    
                    cnt = cnt + 1;     
                end
                fclose(fid);
                
                mpsd_save(isps,:,:) = psd(:,:);
                if (any(chrsld == chrsld_save(isps)))  
                    mpsd(find(chrsld==chrsld_save(isps)),:,:) = mpsd_save(isps,:,:);
                end 
            end
            
            for isps=1:nsp_sld % SA properties calc with updated PSDs
                if (~incld_rough)  
                    for iz=1:nz
                        ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0.*mpsd(isps,:,iz)'.*dps(:));
                        ssav(isps,iz) = sum( 3d0/(10d0.^ps(:)).*mpsd(isps,:,iz)'.*dps(:));
                    end 
                else 
                    for iz=1:nz
                        ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0 ...
                            *rough_c0*(10d0.^ps(:)).^rough_c1.*mpsd(isps,:,iz)'.*dps(:));
                        ssav(isps,iz) = sum( 3d0/(10d0.^ps(:)) ... 
                            *rough_c0*(10d0.^ps(:))^rough_c1.*mpsd(isps,:,iz)'.*dps(:));
                    end 
                end 
                for iz=1:nz
                    ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*mpsd(isps,:,iz)'.*dps(:));
                end 
            end 
            
        end 
        
        
        % just to obtain so4f 
        print_cb = false ;
        print_loc = './ph.txt';
        
        [ ...
            dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
            ,prox,ph_error,so4f,ph_iter ...% output
            ] = calc_pH_v7_3( ...
            nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
            ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
            ,maq,maqc,mgas,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
            ,print_cb,print_loc,z ...% input 
            ,prox,so4f ...% output
            ) 
        so4fprev = so4f;
        
        time = 0d0;
            
        if (display)         
            fmt=['%s\t' repmat('%7.6E \t',1,nz_disp) '\n'];
            fprintf('\n');
            fprintf('[concs]\n');
            fprintf(fmt,'z', z(1:nz/nz_disp:nz));
            if (nsp_aq>0)  
                fprintf(' < aq species >\n');
                for ispa = 1: nsp_aq
                    fprintf(fmt, chraq(ispa), maq(ispa,1:nz/nz_disp:nz));
                end 
            end 
            if (nsp_sld>0)  
                fprintf(' < sld species >\n');
                for isps = 1: nsp_sld
                    fprintf( fmt, chrsld(isps), msld(isps,1:nz/nz_disp:nz));
                end 
            end 
            if (nsp_gas>0)  
                fprintf(' < gas species >\n');
                for ispg = 1: nsp_gas
                    fprintf( fmt, chrgas(ispg), mgas(ispg,1:nz/nz_disp:nz));
                end 
            end 
        end      
    end
    
    
    [ ...
        ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 ...% output
        ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all ...% output
        ] = coefs_v2( ...
        nz,rg,rg2,tc,sec2yr,tempk_0,pro ...% input
        ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all ...% input
        ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all ...% input
        ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all ...%input
        ) ; 
        
        
    dbl_ref = 0d0;  
    labs = zeros(nsp_sld,1,'logical');
    turbo2 = zeros(nsp_sld,1,'logical');
    nobio = zeros(nsp_sld,1,'logical');
    till = zeros(nsp_sld,1,'logical');
    fick = zeros(nsp_sld,1,'logical');

    if (no_biot); nobio(:) = true; end
    if (biot_turbo2); turbo2(:) = true; end
    if (biot_fick); fick(:) = true; end
    if (biot_labs); labs(:) = true; end
    if (biot_till); till(:) = true; end
        
    save_trans = true; 
    save_trans = false;
    [trans,nonlocal,izml] =  make_transmx(  ...
        labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  ...! input
        );
           
    % --------- loop -----
    fprintf ('about to start time loop\n');
    time = 0;
    it = 0;
    irec_prof = 0;
    irec_flx = 0;

    ict = 0;
    ict_prev = zeros(3,1);
    ict_change = zeros(3,1,'logical');

    count_dtunchanged = 0;
    progress_rate = 0;
    

    %% @@@@@@@@@@@@@@@   start of time integration  @@@@@@@@@@@@@@@@@@@@@@
        
    while (it<nt)
        t1 = cputime;
        
        if (display)  
            fprintf('\n');
            fprintf('-----------------------------------------\n');
            fprintf('%d%s\n', it,"'s time integral");
            fprintf('time =\t%7.6E\t[yr]\n',time); 
            fprintf('\n');
        end
        dt_prev = dt;
        
        if (time>rectime_prof(nrec_prof)); break; end
        
        if (it == 0)  
            maxdt = 0.2d0;
        end 

        if (timestep_fixed)  
            maxdt = 0.2d0;
            if ( time<1d0)   
                maxdt = 1d-4 ;
            elseif (time>=1d0 && time<1d1)  
                maxdt = 1d-3 ;
            elseif (time>=1d1 && time<1d2)  
                maxdt = 1d-2  ;
            end 
            
        end 
        
        if (dt<maxdt) 
            dt = dt*10d0;
            if (dt>maxdt); dt = maxdt; end
        end
        
        % if climate is changing in the model 
        if (any(climate))  
            ict_change(:) = false;
            for iclim = 1:3
                if (climate(iclim)) 
                    switch (iclim)
                        case(1)
                            if (dt > dct(iclim)/10d0); dt = dct(iclim)/10d0; end
                            for ict = 1: nclim(iclim)
                                fprintf('%7.6f\t%7.6f\t%7.6f\n', clim_T(1,ict),mod(time,ctau(iclim)),clim_T(1,ict) + dct(iclim));
                                if ( ...
                                    clim_T(1,ict) <= mod(time,ctau(iclim)) ... 
                                    && clim_T(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) ...
                                    )  
                                    % if (  ...
                                        % mod(time,ctau(iclim)) + dt - clim_T(1,ict) + dct(iclim) ...
                                        % > ctau(iclim) * tol_step_tau ...
                                        % )  
                                        % dt = clim_T(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    % endif 
                                    fprinf('%d\n',ict); 
                                    if (ict ~= ict_prev(iclim)); ict_change(iclim) = true; end
                                    ict_prev(iclim) = ict;
                                    break
                                end 
                            end 
                            if (ict ~= nclim(iclim))  
                                tc = ( clim_T(2,ict+1) - clim_T(2,ict) ) /( clim_T(1,ict+1) - clim_T(1,ict) ) ...
                                    * ( mod(time,ctau(iclim)) - clim_T(1,ict) ) + clim_T(2,ict);
                            elseif (ict == nclim(iclim))  
                                tc = ( clim_T(2,1) - clim_T(2,ict) ) /( dct(iclim)  ) ...
                                    * ( mod(time,ctau(iclim)) - clim_T(1,ict) ) + clim_T(2,ict);
                            end 
                            
                        case(2)
                            if (dt > dct(iclim)/10d0); dt = dct(iclim)/10d0;end
                            for ict = 1: nclim(iclim)
                                fprintf('%7.6f\t%7.6f\t%7.6f\n', clim_q(1,ict),mod(time,ctau(iclim)),clim_q(1,ict) + dct(iclim) );
                                if ( ...
                                    clim_q(1,ict) <= mod(time,ctau(iclim)) ... 
                                    && clim_q(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) ...
                                    )  
                                    % if (  ...
                                        % mod(time,ctau(iclim)) + dt - clim_q(1,ict) + dct(iclim) ...
                                        % > ctau(iclim) * tol_step_tau ...
                                        % )  
                                        % dt = clim_q(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    % end 
                                    fprinf('%d\n',ict); 
                                    if (ict ~= ict_prev(iclim)); ict_change(iclim) = true; end
                                    ict_prev(iclim) = ict;
                                    break
                                end 
                            end 
                            if (ict ~= nclim(iclim))  
                                qin = ( clim_q(2,ict+1) - clim_q(2,ict) ) /( clim_q(1,ict+1) - clim_q(1,ict) ) ...
                                    * ( mod(time,ctau(iclim)) - clim_q(1,ict) ) + clim_q(2,ict);
                            elseif (ict == nclim(iclim))  
                                qin = ( clim_q(2,1) - clim_q(2,ict) ) /( dct(iclim)  ) ...
                                    * ( mod(time,ctau(iclim)) - clim_q(1,ict) ) + clim_q(2,ict);
                            end 
                            
                        case(3)
                            if (dt > dct(iclim)/10d0); dt = dct(iclim)/10d0; end
                            for ict = 1: nclim(iclim)
                                fprintf('%7.6f\t%7.6f\t%7.6f\n', clim_sat(1,ict),mod(time,ctau(iclim)),clim_sat(1,ict) + dct(iclim) );
                                if ( ...
                                    clim_sat(1,ict) <= mod(time,ctau(iclim)) ... 
                                    && clim_sat(1,ict) + dct(iclim) >= mod(time,ctau(iclim)) ...
                                    )  
                                    % if (  ...
                                        % mod(time,ctau(iclim)) + dt - clim_sat(1,ict) + dct(iclim) ...
                                        % > ctau(iclim) * tol_step_tau ...
                                        % )  
                                        % dt = clim_sat(1,ict) + dct(iclim) - mod(time,ctau(iclim))
                                    % end 
                                    fprinf('%d\n',ict); 
                                    if (ict ~= ict_prev(iclim)); ict_change(iclim) = true; end
                                    ict_prev(iclim) = ict;
                                    break
                                end 
                            end 
                            if (ict ~= nclim(iclim))  
                                satup = ( clim_sat(2,ict+1) - clim_sat(2,ict) ) /( clim_sat(1,ict+1) - clim_sat(1,ict) ) ...
                                    * ( mod(time,ctau(iclim)) - clim_sat(1,ict) ) + clim_sat(2,ict);
                            elseif (ict == nclim(iclim))  
                                satup = ( clim_sat(2,1) - clim_sat(2,ict) ) /( dct(iclim)  ) ...
                                    * ( mod(time,ctau(iclim)) - clim_sat(1,ict) ) + clim_sat(2,ict);
                            end 
                    end
                end 
            end 
            if (dt >= minval(dct)/10d0 || any (ict_change) )  
                idust = sprintf('%s/climate.txt',flxdir);
                fid = fopen(idust,'a');
                fmt=[repmat('%7.6E \t',1,4) '\n'];
                fprintf(fid,fmt, time,tc,qin,satup);
                fclose(fid);
            end 
            
            sat(:) = min(1.0d0, 1d0-(1d0-satup)*(1d0-z(:)/zsat).^2d0);
            for iz=1:nz
                if (z(iz)>=zsat); sat(iz)=1d0; end
            end 
            v(:) = qin/poroi./sat(:);
            torg(:) = poro(:).^(3.4d0-2.0d0).*(1.0d0-sat(:)).^(3.4d0-1.0d0);
            tora(:) = poro(:).^(3.4d0-2.0d0).*(sat(:)).^(3.4d0-1.0d0);
            
        end 
    
        [ ...
            ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 ...% output
            ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all ...% output
            ] = coefs_v2( ...
            nz,rg,rg2,tc,sec2yr,tempk_0,pro ...% input
            ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all ...% input
            ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all ...% input
            ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all ...%input
            ); 
        
        ksld = zeros(nsp_sld,nz,'double');
        for isps = 1: nsp_sld
            ksld(isps,:) = ksld_all(find(chrsld_all==chrsld(isps)),:);
            % print *,chrsld(isps),ksld(isps,:)
        end
        
        daq = zeros(nsp_aq,1,'double');
        for ispa = 1: nsp_aq 
            daq(ispa) = daq_all(find(chraq_all==chraq(ispa)));
        end 
        
        dgasa = zeros(nsp_gas,1,'double');
        dgasg = zeros(nsp_gas,1,'double');
        for ispg = 1: nsp_gas 
            dgasa(ispg) = dgasa_all(find(chrgas_all==chrgas(ispg)));
            dgasg(ispg) = dgasg_all(find(chrgas_all==chrgas(ispg)));
        end 
        kho = keqgas_h(find(chrgas_all=='po2'),ieqgas_h0);
        kco2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0);
        knh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h0);
        kn2o = keqgas_h(find(chrgas_all=='pn2o'),ieqgas_h0);
        k1 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);
        k2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2);
        k1nh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h1);
        pco2i = mgasi_all(find(chrgas_all=='pco2'));
        pnh3i = mgasi_all(find(chrgas_all=='pnh3'));
        khco2i = kco2*(1d0+k1/proi + k1*k2/proi/proi);
        khnh3i = knh3*(1d0+proi/k1nh3);
        
        khgasi = zeros(nsp_gas,1,'double');
        for ispg = 1: nsp_gas 
            switch (chrgas(ispg))
                case('pco2')
                    khgasi(ispg) = khco2i;
                case('po2')
                    khgasi(ispg) = kho;
                case('pnh3')  
                    khgasi(ispg) = khnh3i;
                case('pn2o')  
                    khgasi(ispg) = kn2o;
            end
        end
        
        % kinetic inhibition 
        if (al_inhibit)  
            if (any(chraq == 'al'))  
                for isps = 1: nsp_sld
                    if (staq(isps,find(chraq=='al')) ~= 0d0)  
                        ksld(isps,:) = ksld(isps,:) ...
                            *10d0^(-4.84d0)./(10d0^(-4.84d0)+maq(find(chraq=='al'),:)); 
                    end 
                end 
            end 
        end 
        
        % nobio = true
        save_trans = false;
        [trans,nonlocal,izml] = make_transmx(  ...
            labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  ...% input
            );

        error = 1d4;
        
        flg_100 = true; % flag raised so that at least one loop is gone through 
% 100 continue
                            while (flg_100) % ****************************trying to do go to stuff in matlab
        flg_100 = false; 

        mgasx = mgas;
        msldx = msld;
        maqx = maq;
        
        prox = pro;  
        
        so4f = so4fprev;
        
        poroprev = poro;
        hrprev = hr;
        vprev = v;
        torgprev = torg;
        toraprev = tora;
        wprev = w; 
        
        mblkx = mblk;
        
        % whether or not you are using psd
        psd_old = psd;
        mpsd_old = mpsd;
        psd_error_flg = false;
        
        %  raining dust and OM 
        maqsupp = zeros(nsp_aq,nz,'double');
        mgassupp = zeros(nsp_gas,nz,'double');
        msldsupp = zeros(nsp_sld,nz,'double');
        for isps = 1: nsp_sld
            if (no_biot)  
                msldsupp(isps,:) = rainpowder*rfrc_sld(isps)*exp(-z(:)'/zsupp)/zsupp;
            else 
                msldsupp(isps,1) = rainpowder*rfrc_sld(isps)/dz(1);
            end
        end 
        
        
        % dust options check 
        if (dust_wave && dust_step)  
            error('CAUTION: options of dust_wave and dust_step are both ON');
        end
        
        % if defined wave function is imposed on dust 
        if (dust_wave)  
            % added for MATLAB because I do not know an equivalent to merge in fortran 
            if nearest(time/wave_tau)==floor(time/wave_tau)
                dust_norm_tmp = 2d0;
            else
                dust_norm_tmp = 0d0;
            end
            
            for isps = 1: nsp_sld
                if (no_biot)  
                    msldsupp(isps,:) = msldsupp(isps,:)*dust_norm_tmp;
                    % msldsupp(isps,:) = msldsupp(isps,:)*merge(2d0,0d0,nearest(time/wave_tau)==floor(time/wave_tau))
                else 
                    msldsupp(isps,1) = msldsupp(isps,1)*dust_norm_tmp;
                    % msldsupp(isps,1) = msldsupp(isps,1)*merge(2d0,0d0,nearest(time/wave_tau)==floor(time/wave_tau))
                end
            end 
            % if (time==0d0 || dust_norm ~= merge(2d0,0d0,nearest(time/wave_tau)==floor(time/wave_tau))) 
            if (time==0d0 || dust_norm ~= dust_norm_tmp ) 
                idust = sprintf('%s/dust.txt',flxdir);
                fid = fopen(idust,'a');
                fmt=[repmat('%7.6E \t',1,2) '\n'];
                fprintf(fid,fmt,time-dt,dust_norm);
                fprintf(fid,fmt,time-dt,dust_norm_tmp);
                dust_norm = dust_norm_tmp;
                fclose(fid);
            end
        end
        
        % non continueous
        if (dust_step)  
            
            dust_norm_prev = dust_norm;
            
            if (dt > step_tau)  
                dt = step_tau;
                % go to 100
            end
            
            if (time - floor(time) < step_tau && time + dt - floor(time) >= step_tau )  
                if ( abs (step_tau - ( time - floor(time) ) ) > step_tau * tol_step_tau )  
                    dt = step_tau - ( time - floor(time) );
                    % go to 100
                end
            end

            % an attempt to use fickian mixing when applying rock powder
            % when not applying the powder use the chosen mixing (can be fickian)
            labs = false;
            turbo2 = false;
            nobio = false;
            till = false;
            fick = false;   
            if (time - floor(time) >= step_tau)  
                % print *, 'no dust time', time 
                msldsupp(:,:) = 0d0;
                dust_norm = 0d0;
                % only implement chosen mixing 
                if (no_biot); nobio(:) = true; end
                if (biot_turbo2); turbo2(:) = true;end
                if (biot_fick); fick(:) = true;end
                if (biot_labs); labs(:) = true;end
                if (biot_till); till(:) = true;end
                % OM is mixed in fickian
                if (any(chrsld == 'g1'))  
                    fick(find(chrsld=='g1')) = true;
                    nobio(find(chrsld=='g1')) = false;
                    labs(find(chrsld=='g1')) = false;
                    turbo2(find(chrsld=='g1')) = false;
                    till(find(chrsld=='g1')) = false;
                end
                if (any(chrsld == 'g2'))  
                    fick(find(chrsld=='g2')) = true;
                    nobio(find(chrsld=='g2')) = false;
                    labs(find(chrsld=='g2')) = false;
                    turbo2(find(chrsld=='g2')) = false;
                    till(find(chrsld=='g2')) = false;
                end
                if (any(chrsld == 'g3'))  
                    fick(find(chrsld=='g3')) = true;
                    nobio(find(chrsld=='g3')) = false;
                    labs(find(chrsld=='g3')) = false;
                    turbo2(find(chrsld=='g3')) = false;
                    till(find(chrsld=='g3')) = false;
                end
            else 
                % print *, 'dust time %%', time 
                msldsupp(:,:) = msldsupp(:,:)/step_tau;
                dust_norm = 1d0/step_tau;
                % only implement fickian mixing 
                fick(:) = true; 
            end
            % mixing reload
            save_trans = false;
            [trans,nonlocal,izml] = make_transmx(  ...
                labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  ...% input
                  );
            
            if ( dust_norm ~= dust_norm_prev ) 
                idust = sprintf('%s/dust.txt',flxdir);
                fid = fopen(idust,'a');
                fmt=[repmat('%7.6E \t',1,2) '\n'];
                fprintf(fid,fmt,time-dt,dust_norm_prev);
                fprintf(fid,fmt,time-dt,dust_norm);
                fclose(fid);
            end
            
        end
        
        % overload with OM rain 
        for isps = 1: nsp_sld
            if (no_biot)  
                msldsupp(isps,:) = msldsupp(isps,:) ...
                    + plant_rain/12d0*rfrc_sld_plant(isps)*exp(-z(:)'/zsupp_plant)/zsupp_plant; % when plant_
            else 
                msldsupp(isps,1) = msldsupp(isps,1) ...
                    + plant_rain/12d0*rfrc_sld_plant(isps)/dz(1); % when plant_rain is in g_C/m2/yr
            end
        end 
        % when enforcing solid states without previous OM spin-up
        if (sld_enforce && (~read_data))  
            if (any(chrgas=='pco2'))  
                % mgassupp(find(chrgas=='pco2'),:) = plant_rain/12d0*exp(-z/zsupp_plant)/zsupp_plant
                mgassupp(find(chrgas=='pco2'),:) = plant_rain/12d0/ztot;
            end
        end
        
        if (any(isnan(msldsupp)))  
            error('error in dust');
        end
    
    
        % do PSD for raining dust & OM 
        if (do_psd)  
            psu_rain_list = zeros(nps_rain_char,1,'double');
            pssigma_rain_list = zeros(nps_rain_char,1,'double');
            psu_rain_list(:) = [log10(5d-6), log10(20d-6),  log10(50d-6), log10(70d-6)];
            pssigma_rain_list(:) = [0.2d0, 0.2d0,  0.2d0, 0.2d0 ];

            ipsd = strcat(profdir,'/psd_rain.txt');
            fid = fopen(ipsd,'w');
            fmt=['%s \t' repmat('%7.6E \t',1,nps) '%s\n'];
            if (do_psd_full)  
                fprintf(fid,fmt, 'sldsp\log10(radius)', ps(1:nps), 'time');
            else
                fprintf(fid,fmt, 'depth\log10(radius)', ps(1:nps), 'time');
            end 
            
            for iz = 1:nz
            
                % rained particle distribution 
                if (read_data)  
                    psd_rain(:,iz) = 0d0;
                    for ips = 1: nps_rain_char
                        psu_rain = psu_rain_list(ips);
                        pssigma_rain = pssigma_rain_list(ips);
                        psd_rain(:,iz) = psd_rain(:,iz) ...
                            + 1d0/pssigma_rain/sqrt(2d0*pi)*exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain ).^2d0 );
                    end 
                else
                    psu_rain = log10(p80);
                    pssigma_rain = 1d0;
                    pssigma_rain = ps_sigma_std;
                    psd_rain(:,iz) = 1d0/pssigma_rain/sqrt(2d0*pi)*exp( -0.5d0*( (ps(:) - psu_rain)/pssigma_rain ).^2d0 );
                end 
                
                psd_rain(:,iz) = psd_rain(:,iz)/sum(psd_rain(:,iz).*dps(:)); 
                
                psd_rain_tmp = psd_rain(:,iz);

                % balance for volumes
                % sum(msldsupp*mv*1d-6) *dt (m3/m3) must be equal to sum( 4/3(pi)r3 * psd_rain * dps) 
                % where psd is number / bulk m3 / log r
                if (do_psd_full)  
                    fmt=['%s \t' repmat('%7.6E \t',1,nps+1) '\n'];
                    for isps = 1: nsp_sld
                        volsld = msldsupp(isps,iz)*mv(isps)*1d-6;
                        psd_rain(:,iz) = psd_rain_tmp(:) * volsld *dt ...
                            /sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain_tmp(:).*dps(:));
                            
                        if ( abs( ( volsld *dt ...
                            - sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain(:,iz).*dps(:))) ...
                            / ( volsld*dt ...
                            ) ) > tol)  
                            error( '%7.6E\t%7.6E\t%7.6E\n', iz, volsld*dt ...
                                ,sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain(:,iz).*dps(:)) );
                        end 
                        mpsd_rain(isps,:,iz) = psd_rain(:,iz)';
                
                        if (iz==1); fprintf(fid,fmt, chrsld(isps),mpsd_rain(isps,1:nps,iz), time);end
                    end 
                else 
                    fmt=[repmat('%7.6E \t',1,nps+2) '\n'];
                    volsld = sum(msldsupp(:,iz).*mv(:)*1d-6);
                    psd_rain(:,iz) = psd_rain_tmp(:) * volsld *dt ...
                        /sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain_tmp(:).*dps(:) ) ;
                        
                    if ( abs( ( volsld *dt ...
                        - sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain(:,iz).*dps(:))) ...
                        / ( volsld*dt ...
                        ) ) > tol)  
                        error( '%7.6E\t%7.6E\t%7.6E\n', iz, volsld*dt ...
                            ,sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_rain(:,iz).*dps(:)) );
                    end 
                
                    fprintf(fid,fmt, z(iz),psd_rain(1:nps,iz), time );
                end 
            end 
            
            fclose(fid);
        end 
        
        
        if (method_precalc); pre_calc = true; end
        
        if (pre_calc)  
            % to be implemented ...? Possibly no need but could be useful to implementing some time-forward methods
            pre_calc = false;
            if (any(isnan(mgasx))||any(isnan(msldx))||any(isnan(maqx)))  
                error( 'error in precalc');
            end
        end 




        % before start of porosity iteration
        poro_iter = 0;
        poro_error = 1d4;
        poro_tol = 1d-6;
        poro_iter_max = 50;
        
        poro_tol = 1d-9;
        
        
        while (poro_error > poro_tol) % start of porosity iteration 

            porox = poro;
            wx = w;
        
            [ ... 
                flgback,w ...    % output 
                ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp,so4f ... 
                ] = alsilicate_aq_gas_1D_v3_1( ...
                nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext ...
                ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  ...
                ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi ...
                ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp ...
                ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext ...
                ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst ...
                ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all ...
                ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s,keqaq_no3,keqaq_nh3 ...
                ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...
                ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all,stgas_all ...
                ,turbo2,labs,trans,method_precalc,display,chrflx,sld_enforce ...% input
                ,nsld_kinspc,chrsld_kinspc,kin_sld_spc ...% input
                ,precstyle ...% in 
                ,hr,poro,z,dz,w_btm,sat,pro,poroprev,tora,v,tol,it,nflx,kw,so4fprev ... %  old inputs
                ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0,proi,poroi,up,dwn,cnr,adf,msldunit  ...
                ,dt,flgback,w ...    % old inout
                ,msldx,omega,maqx,mgasx ... % inout 
                );
            
            save_trans = false;
            [trans,nonlocal,izml] = make_transmx(  ...
                labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  ...% input
                );
                
            if ( incld_blk )  
                for iz=1:nz
                    mblkx(iz) = 1d0 - poro(iz) - sum(msldx(:,iz).*mv(:)*1d-6);
                    mblkx(iz) = mblkx(iz)/(mvblk*1d-6);
                end 
            else 
                mblkx = 0d0;
            end 

            if (flgback)  
                flgback = false ;
                flgreducedt = true;
                dt = dt/1d1;
                psd = psd_old;
                mpsd = mpsd_old;
                poro = poroprev;
                torg = torgprev;
                tora = toraprev;
                v = vprev;
                hr = hrprev;
                w = wprev;
                [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                % go to 100
                flg_100 = true;
                break % escape from porosity loop
            end    
            
            if (poroevol)  
                if (iwtype == iwtype_flex)  
                    poro(:) = poroi;
                    % not constant but calculated as defined (only applicable when unit of msld(x) is mol per bulk soil)
                    for iz = 1:nz
                        poro(iz) = 1d0 - sum(msldx(:,iz).*mv(:)*1d-6);
                    end 
                else 
                    [poro] = calc_poro( ...
                        nz,nsp_sld,nflx,idif,irain ...% in
                        ,flx_sld,mv,poroprev,w,poroi,w_btm,dz,tol,dt ...% in
                        ,poro ...% inout
                        );
                end 
                
                if (any(poro<0d0))  
                    error('negative porosity: stop');
                    
                    flgback = false; 
                    flgreducedt = true;
                    dt = dt/1d1;
                    psd = psd_old;
                    mpsd = mpsd_old;
                    poro = poroprev;
                    torg = torgprev;
                    tora = toraprev;
                    v = vprev;
                    hr = hrprev;
                    w = wprev;
                    [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                    % go to 100
                    flg_100 = true;
                    break % escape from porosity loop
                end 
                if (any(poro>1d0))  
                    error('porosity exceeds 1: stop');
                    
                    flgback = false; 
                    flgreducedt = true;
                    dt = dt/1d1;
                    psd = psd_old;
                    mpsd = mpsd_old;
                    poro = poroprev;
                    torg = torgprev;
                    tora = toraprev;
                    v = vprev;
                    hr = hrprev;
                    w = wprev;
                    [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                    % go to 100
                    flg_100 = true;
                    break % escape from porosity loop
                end 
                v(:) = qin./poro(:)./sat(:);
                torg(:) = poro(:).^(3.4d0-2.0d0).*(1.0d0-sat(:)).^(3.4d0-1.0d0);
                tora(:) = poro(:).^(3.4d0-2.0d0).*(sat(:)).^(3.4d0-1.0d0);
                
                for isps=1:nsp_sld
                    hr(isps,:) = hri(isps,:).*rough(isps,:);
                    if (surfevol1 )  
                        hr(isps,:) = hri(isps,:).*rough(isps,:).*((1d0-poro(:)')/(1d0-poroi)).^(2d0/3d0);
                    end 
                    if (surfevol2 )  
                        hr(isps,:) = hri(isps,:).*rough(isps,:).*(poro(:)/poroi).^(2d0/3d0);  % SA increases with porosity 
                    end 
                end 
                % if doing psd SA is calculated reflecting psd
                if (do_psd)  
                    for isps=1:nsp_sld
                        if (do_psd_full)  
                            hr(isps,:) = ssa(isps,:)./ssv(isps,:)./poro(:)'; 
                        else
                            hr(isps,:) = ssa(isps,:)./poro(:)'; % so that poro * hr * mv * msld becomes porosity independent 
                        end 
                    end 
                end 
                
                if (display && (~ display_lim))  
                    fmt=['%-6s' repmat('%11.3E',1,nz_disp) '\n'];
                    fprintf ('\n');
                    fprintf (' [porosity & surface area]\n');
                    fprintf (fmt,'z',z(1:nz/nz_disp:nz) );
                    fprintf (fmt,'poro',poro(1:nz/nz_disp:nz) );
                    fprintf (fmt,'SA',hrb(1:nz/nz_disp:nz) );
                    fprintf ('\n');
                end 
            end  
            
            if (poroiter_in)  
                if (iwtype == iwtype_flex)  
                    % poro_error = max ( abs (( w(:) - wx(:) )./wx(:) ) );
                    poro_error = max ( abs ( w(:) - wx(:) ) );
                else
                    % poro_error = max ( abs (( poro(:) - porox(:) )./porox(:) ) );
                    poro_error = max ( abs ( poro(:) - porox(:) ) );
                end 
                
                fprintf( "%d's porosity iteration with an error of \t%7.6E\n ",poro_iter,poro_error);
                poro_iter = poro_iter + 1;
                            
                if (poro_iter > poro_iter_max)  
                    fprintf( 'too much porosity iteration but does not converge within assumed threshold\n');
                    fprintf( 'reducing dt and move back\n');
                    flgback = false; 
                    flgreducedt = true;
                    psd = psd_old;
                    mpsd = mpsd_old;
                    poro = poroprev;
                    torg = torgprev;
                    tora = toraprev;
                    v = vprev;
                    hr = hrprev;
                    w = wprev;
                    [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                    dt = dt/1d1;
                    % go to 100
                    flg_100 = true;
                    break % escape from porosity loop
                end    
            else 
                poro_error = 0d0;
                break
            end 
        
        end % porosity iteration end
        
        if flg_100; continue; end %  added to enable 'go to 100'
        
        
        % attempt to do psd
        if (do_psd) 
            
            if (display) 
                fprintf('\n');
                fprintf('-- doing PSD\n');
            end 
            
            dmpsd = zeros(nsp_sld,nps,nz,'double');
            DV = zeros(nz,1,'double');
            dpsd = zeros(nps,nz,'double'); 
            flx_mpsd = zeros(nsp_sld,nps,nflx_psd,nz,'double');
            
            if (do_psd_full) 
                
                for isps = 1: nsp_sld
                
                    if ( precstyle(isps) == 'decay'); continue; end % solid species not related to SA
                
                    DV(:) = flx_sld(isps, 4 + isps,:)*mv(isps)*1d-6*dt;  
                
                    dpsd(:,:) = 0d0;
                    psd(:,:) = mpsd(isps,:,:);
                
                    [dpsd,psd_error_flg] = psd_diss( ...
                        nz,nps ...% in
                        ,z,DV,dt,pi,tol_dvd,poro ...% in 
                        ,incld_rough,rough_c0,rough_c1 ...% in
                        ,psd,ps,dps,ps_min,ps_max ...% in 
                        ,chrsld(isps) ...% in 
                        ,dpsd,psd_error_flg ...% inout
                        );
                        
                    if (psd_error_flg) 
                        flgback = false;
                        flgreducedt = true;
                        psd = psd_old;
                        mpsd = mpsd_old;
                        poro = poroprev;
                        torg = torgprev;
                        tora = toraprev;
                        v = vprev;
                        hr = hrprev;
                        w = wprev;
                        [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                        dt = dt/1d1;
                        % go to 100
                        flg_100 = true;
                        continue
                    end 
                    
                    dmpsd(isps,:,:) = dpsd(:,:);
                
                end 
        
            else 

                DV(:) = 0d0;
                for isps = 1:nsp_sld 
                    for iz=1:nz
                        % DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6*dt/(1d0 - poro(iz))  
                        DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6*dt;  
                    end 
                end 

                dpsd(:,:) = 0d0;
                
                [dpsd,psd_error_flg] = psd_diss( ...
                    nz,nps ...% in
                    ,z,DV,dt,pi,tol,poro ...% in 
                    ,incld_rough,rough_c0,rough_c1 ...% in
                    ,psd,ps,dps,ps_min,ps_max ...% in 
                    ,'blk' ...% in 
                    ,dpsd,psd_error_flg ...% inout
                    );
                    
                if (psd_error_flg) 
                    flgback = false;
                    flgreducedt = true;
                    psd = psd_old;
                    mpsd = mpsd_old;
                    poro = poroprev;
                    torg = torgprev;
                    tora = toraprev;
                    v = vprev;
                    hr = hrprev;
                    w = wprev;
                    [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                    dt = dt/1d1;
                    % go to 100
                    flg_100 = true;
                    continue
                end 
            
            end 
            
            if (do_psd_norm) 
            
                if (do_psd_full) 
                    
                    flx_max_max = 0d0;
                    
                    for isps=1:nsp_sld
                        
                        if ( precstyle(isps) == 'decay' ) % solid species not related to SA                    
                            flx_mpsd(isps,:,:,:) = 0d0;
                            for iz=1:nz
                                mpsdx(isps,:,iz) = mpsd_pr(isps,:);
                            end 
                            continue 
                        end 
                    
                        for ips=1:nps
                            psd_norm_fact(ips) = max(mpsd(isps,ips,:));
                            
                            psd_norm(ips,:) = mpsd(isps,ips,:) / psd_norm_fact(ips);
                            psd_pr_norm(ips) = mpsd_pr(isps,ips) / psd_norm_fact(ips);
                            dpsd_norm(ips,:) = dmpsd(isps,ips,:) / psd_norm_fact(ips);
                            psd_rain_norm(ips,:) = mpsd_rain(isps,ips,:) / psd_norm_fact(ips);
                        end 
                        
                        [ ...
                            flgback,flx_max_max ...% inout
                            ,psdx_norm,flx_psd_norm ...% out
                            ] = psd_implicit_all_v2( ...
                            nz,nsp_sld,nps,nflx_psd ...% in
                            ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev ...% in 
                            ,incld_rough,rough_c0,rough_c1 ...% in
                            ,trans ...% in
                            ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm ...% in  
                            ,chrsld(isps) ...% in 
                            ,flgback,flx_max_max ...% inout
                            );
                        
                        if (flgback) 
                            flgback = false;
                            flgreducedt = true;
                            psd = psd_old;
                            mpsd = mpsd_old;
                            poro = poroprev;
                            torg = torgprev;
                            tora = toraprev;
                            v = vprev;
                            hr = hrprev;
                            w = wprev;
                            [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                            dt = dt/1d1;
                            % go to 100
                            flg_100 = true;
                            continue
                        end 
                            
                        for ips=1:nps
                            mpsdx(isps,ips,:) = psdx_norm(ips,:)*psd_norm_fact(ips);
                            flx_mpsd(isps,ips,:,:) = flx_psd_norm(ips,:,:)*psd_norm_fact(ips);
                        end 
                    
                    end 
                
                else 
                    
                    for ips=1:nps
                        psd_norm_fact(ips) = max(psd(ips,:));
                        
                        psd_norm(ips,:) = psd(ips,:) / psd_norm_fact(ips);
                        psd_pr_norm(ips) = psd_pr(ips) / psd_norm_fact(ips);
                        dpsd_norm(ips,:) = dpsd(ips,:) / psd_norm_fact(ips);
                        psd_rain_norm(ips,:) = psd_rain(ips,:) / psd_norm_fact(ips);
                    end 
                    
                    flx_max_max = 0d0;
                    
                    [...
                        flgback,flx_max_max ...% inout
                        ,psdx_norm,flx_psd_norm ...% out
                        ] = psd_implicit_all_v2( ...
                        nz,nsp_sld,nps,nflx_psd ...% in
                        ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev ...% in 
                        ,incld_rough,rough_c0,rough_c1 ...% in
                        ,trans ...% in
                        ,psd_norm,psd_pr_norm,ps,dps,dpsd_norm,psd_rain_norm ...% in    
                        ,'blk' ...% in 
                        ,flgback,flx_max_max ...% inout
                        );
                    
                    if (flgback) 
                        flgback = false;
                        flgreducedt = true;
                        psd = psd_old;
                        mpsd = mpsd_old;
                        poro = poroprev;
                        torg = torgprev;
                        tora = toraprev;
                        v = vprev;
                        hr = hrprev;
                        w = wprev;
                        [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                        dt = dt/1d1;
                        % go to 100
                        flg_100 = true;
                        continue
                    end 
                        
                    for ips=1:nps
                        psdx(ips,:) = psdx_norm(ips,:)*psd_norm_fact(ips);
                        flx_psd(ips,:,:) = flx_psd_norm(ips,:,:)*psd_norm_fact(ips);
                    end 
                
                end 
            else
                
                if (do_psd_full) 
                
                    flx_max_max = 0d0;
                
                    for isps = 1: nsp_sld
                        
                        if ( precstyle(isps) == 'decay' ) % solid species not related to SA                    
                            flx_mpsd(isps,:,:,:) = 0d0;
                            for iz=1:nz
                                mpsdx(isps,:,iz) = mpsd_pr(isps,:);
                            end 
                            continue 
                        end 
            
                        if (display) 
                            fprintf('\n');
                            fprintf('<%s>\n',chrsld(isps));
                        end 
                        
                        psd(:,:) = mpsd(isps,:,:);
                        psd_pr(:) = mpsd_pr(isps,:);
                        dpsd(:,:) = dmpsd(isps,:,:);
                        psd_rain(:,:) = mpsd_rain(isps,:,:);
                
                        [...
                            flgback,flx_max_max ...% inout
                            ,psdx,flx_psd ...% out
                            ] = psd_implicit_all_v2( ...
                            nz,nsp_sld,nps,nflx_psd ...% in
                            ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev ...% in 
                            ,incld_rough,rough_c0,rough_c1 ...% in
                            ,trans ...% in
                            ,psd,psd_pr,ps,dps,dpsd,psd_rain ...% in    
                            ,chrsld(isps) ...% in 
                            ,flgback,flx_max_max ...% inout
                            );
                
                        if (flgback) 
                            flgback = false;
                            flgreducedt = true;
                            psd = psd_old;
                            mpsd = mpsd_old;
                            poro = poroprev;
                            torg = torgprev;
                            tora = toraprev;
                            v = vprev;
                            hr = hrprev;
                            w = wprev;
                            [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                            dt = dt/1d1;
                            % go to 100
                            flg_100 = true;
                            continue
                        end 
                        
                        mpsdx(isps,:,:) = psdx(:,:);
                        flx_mpsd(isps,:,:,:) = flx_psd(:,:,:);
                    end
                
                else 
                
                    flx_max_max = 0d0;
                    
                    [...
                        flgback,flx_max_max ...% inout
                        ,psdx,flx_psd ...% out
                        ] = psd_implicit_all_v2( ...
                        nz,nsp_sld,nps,nflx_psd ...% in
                        ,z,dz,dt,pi,tol,w_btm,w,poro,poroi,poroprev ...% in 
                        ,incld_rough,rough_c0,rough_c1 ...% in
                        ,trans ...% in
                        ,psd,psd_pr,ps,dps,dpsd,psd_rain ...% in    
                        ,'blk' ...% in 
                        ,flgback,flx_max_max ...% inout
                        );
                
                    if (flgback) 
                        flgback = false;
                        flgreducedt = true;
                        psd = psd_old;
                        mpsd = mpsd_old;
                        poro = poroprev;
                        torg = torgprev;
                        tora = toraprev;
                        v = vprev;
                        hr = hrprev;
                        w = wprev;
                        [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                        dt = dt/1d1;
                        % go to 100
                        flg_100 = true;
                        continue
                    end    
                
                end 
                
            end 
            
            if (do_psd_full) 
            
                mpsd = mpsdx; 

                if (any(isnan(mpsd),'all')) 
                    fprintf( 'nan in mpsd\n');
                    error('stop');
                end 
                if (any(mpsd<0d0,'all')) 
                    error_psd = 0d0;
                    for isps = 1: nsp_sld
                        for iz = 1: nz
                            for ips=1:nps
                                if (mpsd(isps,ips,iz)<0d0) 
                                    error_psd = min(error_psd,mpsd(isps,ips,iz));
                                    mpsd(isps,ips,iz) = 0d0;
                                end 
                            end 
                        end 
                    end 
                    if (abs(error_psd/max(mpsd,[],'all')) > tol) 
                        fprintf( 'negative mpsd\n');
                        flgback = false;
                        flgreducedt = true;
                        psd = psd_old;
                        mpsd = mpsd_old;
                        poro = poroprev;
                        torg = torgprev;
                        tora = toraprev;
                        v = vprev;
                        hr = hrprev;
                        w = wprev;
                        [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                        dt = dt/1d1;
                        % go to 100
                        flg_100 = true;
                        continue
                    end 
                end 
                
                % trancating small psd 
                if (psd_lim_min) 
                    mpsd (mpsd < psd_th)  = psd_th;
                end
            
            else 
            
                psd = psdx;

                if (any(isnan(psd),'all')) 
                    fprintf( 'nan in psd\n');
                    error('stop');
                end 
                if (any(psd<0d0,'all')) 
                    error_psd = 0d0;
                    for iz = 1: nz
                        for ips=1:nps
                            if (psd(ips,iz)<0d0) 
                                error_psd = min(error_psd,psd(ips,iz));
                                psd(ips,iz) = 0d0;
                            end 
                        end 
                    end 
                    if (abs(error_psd/max(psd,[],'all')) > tol) 
                        fprintf('negative psd\n');
                        flgback = false;
                        flgreducedt = true;
                        psd = psd_old;
                        mpsd = mpsd_old;
                        poro = poroprev;
                        torg = torgprev;
                        tora = toraprev;
                        v = vprev;
                        hr = hrprev;
                        w = wprev;
                        [up,dwn,cnr,adf] = calcupwindscheme(w,nz);
                        dt = dt/1d1;
                        % go to 100
                        flg_100 = true;
                        continue
                    end 
                end 
                
                % trancating small psd 
                if (psd_lim_min) 
                    psd (psd < psd_th)  = psd_th;
                end 
            end 
            
            if (display) 
                fprintf( '-- ending PSD\n');
                fprintf('\n');
            end 
        
        else 
            psd(:,:) = 0d0;
            mpsd(:,:,:) = 0d0;
        end 
        % ************************************* if flg_100 is raised within while loop, all is returned to the start of while loop 
                            end % ********************************* go to stuff loop end 
        
        % if doing psd SA is calculated reflecting psd
        if (do_psd)  
            if (do_psd_full)   
                for isps=1:nsp_sld
                    if (~ incld_rough)  
                        for iz=1:nz
                            ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0.*mpsd(isps,:,iz)'.*dps(:));
                            ssav(isps,iz) = sum( 3d0/(10d0.^ps(:)).*mpsd(isps,:,iz)'.*dps(:));
                        end 
                    else 
                        for iz=1:nz
                            ssa(isps,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0  ...
                                *rough_c0.*(10d0.^ps(:)).^rough_c1.*mpsd(isps,:,iz)'.*dps(:));
                            ssav(isps,iz) = sum( 3d0./(10d0.^ps(:))  ...
                                *rough_c0.*(10d0.^ps(:)).^rough_c1.*mpsd(isps,:,iz)'.*dps(:));
                        end 
                    end
                    for iz=1:nz
                        ssv(isps,iz) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*mpsd(isps,:,iz)'.*dps(:));
                    end 
                end 
            else 
                if (~ incld_rough)  
                    for iz=1:nz
                        ssa(:,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0.*psd(:,iz).*dps(:));
                        ssav(:,iz) = sum( 3d0/(10d0.^ps(:)).*psd(:,iz).*dps(:));
                    end 
                else 
                    for iz=1:nz
                        ssa(:,iz) = sum( 4d0*pi*(10d0.^ps(:)).^2d0*rough_c0.*(10d0.^ps(:)).^rough_c1.*psd(:,iz).*dps(:));
                        ssav(:,iz) = sum( 3d0./(10d0.^ps(:))*rough_c0.*(10d0.^ps(:)).^rough_c1.*psd(:,iz).*dps(:));
                    end 
                end
                for iz=1:nz
                    ssv(:,iz) = sum( 4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd(:,iz).*dps(:));
                end 
            end 
            for isps=1:nsp_sld
                if (do_psd_full)  
                    hr(isps,:) = ssa(isps,:)./ssv(isps,:)./poro(:)';
                else 
                    hr(isps,:) = ssa(isps,:)./poro(:)'; % so that poro * hr * mv * msld becomes porosity independent
                end 
            end
        end 
        
        if (any(poro < 1d-10))  
            error ('***| too small porosity: going to end sim as likely ending up crogging ');
            error ( '***| ... and no more reasonable simulation ...  ');
        end 

        if (display  && (~ display_lim))  
            fmt=['%-6s' repmat('%11.3E',1,nz_disp) '\n'];
            
            fprintf ('\n');
            fprintf (' [concs] \n');
            fprintf (fmt,'z',z(1:nz/nz_disp:nz) );
            if (nsp_aq>0)  
                fprintf(' < aq species >\n');
                for ispa = 1: nsp_aq
                    fprintf (fmt,chraq(ispa),maqx(ispa,1:nz/nz_disp:nz) );
                end 
            end 
            if (nsp_sld>0)  
                fprintf (' < sld species >\n');
                for isps = 1: nsp_sld
                    fprintf (fmt,chrsld(isps),msldx(isps,1:nz/nz_disp:nz) );
                end 
            end 
            if (nsp_gas>0)  
                fprintf (' < gas species >\n');
                for ispg = 1, nsp_gas
                    fprintf (fmt,chrgas(ispg),mgasx(ispg,1:nz/nz_disp:nz) );
                end 
            end 
            
            fprintf ('\n');
            fprintf (' [saturation & pH] \n');
            if (nsp_sld>0)  
                fprintf (' < sld species omega > \n');
                for isps = 1: nsp_sld
                    fprintf (fmt,chrsld(isps),omega(isps,1:nz/nz_disp:nz) );
                end 
            end 
            fprintf (' < pH >\n');
            fprintf (fmt,'ph',-log10(prox(1:nz/nz_disp:nz)) );
            
            
            fmt=['%-6s' repmat('%11s',1,nflx) '\n'];
            
            fprintf ('\n');
            fprintf (' [fluxes] \n');
            fprintf (fmt, ' ',chrflx(1:nflx));
            
            if (nsp_aq>0)  
                fprintf (' < aq species >\n');
                for ispa = 1: nsp_aq
                    fprintf ('%-6s',chraq(ispa));
                    for iflx=1:nflx
                        fprintf ('%11.3E',sum(squeeze(flx_aq(ispa,iflx,:)).*dz(:)));
                    end 
                    fprintf ('\n');
                end 
            end 
            if (nsp_sld>0)  
                fprintf (' < sld species >\n');
                for isps = 1: nsp_sld
                    fprintf ('%-6s',chrsld(isps));
                    for iflx=1:nflx
                        fprintf ('%11.3E',sum(squeeze(flx_sld(isps,iflx,:)).*dz(:)));
                    end 
                    fprintf ('\n');
                end 
            end 
            if (nsp_gas>0)  
                fprintf (' < gas species >\n');
                for ispg = 1, nsp_gas
                    fprintf ('%-11s',chrgas(ispg));
                    for iflx=1:nflx
                        fprintf ('%11.3E\t',sum(squeeze(flx_gas(ispg,iflx,:)).*dz(:)));
                    end 
                    fprintf ('\n');
                end 
            end 
            
            if (do_psd)  
                if (do_psd_full) 
                    fprintf ('\n');
                    fprintf (' [fluxes -- PSD]\n');
                    fmt =['%-6s%-6s' repmat('%11s',1,nflx_psd) '\n'];
                    fprintf (fmt,'sld','rad','tflx','adv','dif','rain','rxn','res');
                    for isps = 1: nsp_sld
                        for ips = 1: nps
                            fprintf('%-6s%-6f',chrsld(isps), ps(ips));
                            for iflx = 1:nflx_psd
                                fprintf('%11.3E', sum(squeeze(flx_mpsd(isps,ips,iflx,:)).*dz(:))); 
                            end 
                            fprintf ('\n');
                        end 
                        fprintf ('\n');
                    end 
                else 
                    fprintf ('\n');
                    fprintf (' [fluxes -- PSD]\n');
                    fmt =['%-6s' repmat('%11s',1,nflx_psd) '\n'];
                    fprintf (fmt,'rad','tflx','adv','dif','rain','rxn','res');
                    for ips = 1: nps
                        fprintf('%-6f',ps(ips));
                        for iflx = 1:nflx_psd
                            fprintf('%11.3E', sum(squeeze(flx_psd(ips,iflx,:)).*dz(:))); 
                        end 
                        fprintf ('\n');
                    end 
                end 
            end
            
            if (display_lim_in); display_lim = true; end 
            
        end 

        if (lim_minsld_in)  
            for isps = 1: nsp_sld
                for iz=1:nz
                    if (msldx(isps,iz) < minsld(isps)); msldx(isps,iz) = minsld(isps); end 
                end 
            end 
        end 
        
        mgas = mgasx;
        maq = maqx;
        msld = msldx;
        
        pro = prox;
        so4fprev = so4f;
        
        mblk = mblkx;
        
        for iz = 1: nz
            % accounting for blk soil
            rho_grain_z(iz) = sum(msldx(:,iz).*mwt(:)*1d-6) + mblkx(iz)*mwtblk*1d-6;
            sldvolfrac(iz) = sum(msldx(:,iz).*mv(:)*1d-6) + mblkx(iz)*mvblk*1d-6;
            
            if (msldunit=='blk')  
                rho_grain_z(iz) = rho_grain_z(iz) / ( 1d0 - poro(iz) );
                sldvolfrac(iz) = sldvolfrac(iz) / ( 1d0 - poro(iz) );
            end 
        end 
        
        % calculating volume weighted average surface area 
        for iz=1:nz
            hrb(iz) = sum( hr(:,iz).* msldx(:,iz).*mv(:)*1d-6) / sum(msldx(:,iz).*mv(:)*1d-6);
            ssab(iz) = sum( ssa(:,iz).* msldx(:,iz).*mv(:)*1d-6) / sum(msldx(:,iz).*mv(:)*1d-6);
        end 
        
        
        if (time > savetime)  
            
            isldprof = fopen (strcat(profdir,'/prof_sld-save.txt'), 'w');
            igasprof = fopen (strcat(profdir,'/prof_gas-save.txt'), 'w');
            iaqprof = fopen (strcat(profdir,'/prof_aq-save.txt'), 'w');
            ibsd = fopen (strcat(profdir,'/bsd-save.txt'), 'w');
            isa = fopen (strcat(profdir,'/sa-save.txt'), 'w');
            
            fprintf(isldprof,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'],'z',chrsld(1:nsp_sld),'time');
            fprintf(iaqprof,[repmat('%s\t',1,max(1,nsp_aq)+3) '\n'],'z',chraq(1:nsp_aq),'ph','time');
            fprintf(igasprof,[repmat('%s\t',1,max(1,nsp_gas)+2) '\n'], 'z',chrgas(1:nsp_gas),'time');
            fprintf(ibsd,[repmat('%s\t',1, 10) '\n'], 'z','poro', 'sat', 'v[m/yr]', 'm2/m3' , 'w[m/yr]'  ...
                , 'vol[m3/m3]','dens[g/cm3]', 'blk[wt%]','time');
            fprintf(isa,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'],'z',chrsld(1:nsp_sld),'time');

            for iz = 1: nz
                ucvsld1 = 1d0;
                if (msldunit == 'blk'); ucvsld1 = 1d0 - poro(iz); end
                fprintf(isldprof,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz),msldx(1:nsp_sld,iz),time);
                fprintf(igasprof,[repmat('%e\t',1,max(1,nsp_gas)+2) '\n'], z(iz),mgasx(1:nsp_gas,iz),time);
                fprintf(iaqprof,[repmat('%e\t',1,max(1,nsp_aq)+3) '\n'],z(iz),maqx(1:nsp_aq,iz),-log10(prox(iz)),time);
                fprintf(ibsd,[repmat('%e\t',1, 10) '\n'],z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz) ...
                    ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6), time);
                fprintf(isa,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz),hr(1:nsp_sld,iz),time);
            end 

            fclose(isldprof);
            fclose(iaqprof);
            fclose(igasprof);
            fclose(ibsd);
            fclose(isa);
            
            if (do_psd)  
                if (do_psd_full)  
                    for isps=1:nsp_sld
                        ipsd = fopen (strcat(profdir,'/psd_',chrsld(isps),'-save.txt'), 'w');
                        fprintf(ipsd,['%s\t' repmat('%f\t',1,nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                        for iz = 1: nz
                            fprintf(ipsd,[repmat('%e\t',1,nps+2) '\n'], z(iz), mpsd(isps,1:nps,iz), time );
                        end 
                        fclose(ipsd);
                    end 
                else 
                    ipsd = fopen (strcat(profdir,'/psd-save.txt'), 'w');
                    fprintf(ipsd,['%s\t' repmat('%f\t',1,nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                    for iz = 1: nz
                        fprintf(ipsd,[repmat('%e\t',1,nps+2) '\n'], z(iz), psd(1:nps,iz), time );
                    end 
                    fclose(ipsd);
                end 
            end 
            
            savetime = savetime + dsavetime;
            
        end 
        
        
        
        flx_recorded = false;
        
        if (time>=rectime_prof(irec_prof+1)) 
            chr = sprintf('%3.3d',irec_prof+1);
            
            
            print_cb = true; 
            print_loc = sprintf('%s/chrge_balance-%3.3d.txt',profdir,irec_prof+1);
            
            [ ...
                dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
                ,prox,ph_error,so4f,ph_iter ...% output
                ] = calc_pH_v7_3( ...
                nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
                ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
                ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
                ,print_cb,print_loc,z ...% input 
                ,prox,so4f ...% output
                ); 
            
            isldprof = fopen(strcat(profdir,'/prof_sld-',chr,'.txt'), 'w');
            isldprof2 = fopen(strcat(profdir,'/prof_sld(wt%)-',chr,'.txt'), 'w');
            isldprof3 = fopen(strcat(profdir,'/prof_sld(v%)-',chr,'.txt'), 'w');
            isldsat = fopen(strcat(profdir,'/sat_sld-',chr,'.txt'), 'w');
            igasprof = fopen(strcat(profdir,'/prof_gas-',chr,'.txt'), 'w');
            iaqprof = fopen(strcat(profdir,'/prof_aq-',chr,'.txt'), 'w');
            ibsd = fopen(strcat(profdir,'/bsd-',chr,'.txt'), 'w');
            irate = fopen(strcat(profdir,'/rate-',chr,'.txt'), 'w');
            isa = fopen(strcat(profdir,'/sa-',chr,'.txt'), 'w');
                
            fprintf(isldprof,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'], 'z',chrsld(1:nsp_sld),'time');
            fprintf(isldprof2,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'], 'z',chrsld(1:nsp_sld),'time');
            fprintf(isldprof3,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'], 'z',chrsld(1:nsp_sld),'time');
            fprintf(isldsat,[repmat('%s\t',1,max(1,nsp_sld)+2) '\n'],'z',chrsld(1:nsp_sld),'time');
            fprintf(iaqprof,[repmat('%s\t',1,max(1,nsp_aq)+3) '\n'],'z',chraq(1:nsp_aq),'ph','time');
            fprintf(igasprof,[repmat('%s\t',1,max(1,nsp_gas)+2) '\n'], 'z',chrgas(1:nsp_gas),'time');
            fprintf(ibsd,[repmat('%s\t',1,10) '\n'], 'z','poro', 'sat', 'v[m/yr]', 'm2/m3' , 'w[m/yr]' ...
                , 'vol[m3/m3]','dens[g/cm3]','blk[wt%]','time');
            fprintf(irate,[repmat('%s\t',1, 2 + max(1,nsp_sld) + max(1,nrxn_ext)) '\n'], 'z', chrsld(1:nsp_sld),chrrxn_ext(1:nrxn_ext),'time');
            fprintf(isa,[repmat('%s\t',1, max(1,nsp_sld) +2) '\n'], 'z',chrsld(1:nsp_sld),'time');

            for iz = 1: nz
                ucvsld1 = 1d0;
                if (msldunit == 'blk'); ucvsld1 = 1d0 - poro(iz); end
                
                fprintf(isldprof,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz),msldx(1:nsp_sld,iz),time);
                fprintf(isldprof2,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz) ...
                    ,msldx(1:nsp_sld,iz).*mwt(1:nsp_sld)'*1d2/ucvsld1/(rho_grain_z(iz)*1d6),time);
                fprintf(isldprof3,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz) ...
                    ,msldx(1:nsp_sld,iz).*mv(1:nsp_sld)'/ucvsld1*1d-6*1d2,time);
                fprintf(isldsat,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz),omega(1:nsp_sld,iz),time);
                fprintf(igasprof,[repmat('%e\t',1,max(1,nsp_gas)+2) '\n'], z(iz),mgasx(1:nsp_gas,iz),time);
                fprintf(iaqprof,[repmat('%e\t',1,max(1,nsp_aq)+3) '\n'], z(iz),maqx(1:nsp_aq,iz),-log10(prox(iz)),time);
                fprintf(ibsd,[repmat('%e\t',1, 10) '\n'], z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz)  ...
                    ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6),time);
                fprintf(irate,[repmat('%e\t',1,max(1,nsp_sld)+max(1,nrxn_ext)+2) '\n'], z(iz),rxnsld(1:nsp_sld,iz),rxnext(1:nrxn_ext,iz),time );
                fprintf(isa,[repmat('%e\t',1,max(1,nsp_sld)+2) '\n'], z(iz),hr(1:nsp_sld,iz),time);
            end

            fclose(isldprof);
            fclose(isldprof2);
            fclose(isldprof3);
            fclose(isldsat);
            fclose(iaqprof);
            fclose(igasprof);
            fclose(ibsd);
            fclose(irate);
            fclose(isa);
            
            if (do_psd)  
                
                if (do_psd_full)  
                    
                    for isps = 1: nsp_sld 
                        ipsd = fopen(strcat(profdir,'/psd_',chrsld(isps),'-',chr,'.txt'), 'w');
                        ipsdv = fopen(strcat(profdir,'/psd_',chrsld(isps),'(v%)-',chr,'.txt'), 'w');
                        ipsds = fopen(strcat(profdir,'/psd_',chrsld(isps),'(SA%)-',chr,'.txt'), 'w');
                        ipsdflx = fopen(strcat(flxdir,'/flx_psd_',chrsld(isps),'-',chr,'.txt'), 'w');

                        fprintf(ipsd,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                        fprintf(ipsdv,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                        fprintf(ipsds,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                        fprintf(ipsdflx,[repmat('%s\t',1, nflx_psd + 2) '\n'], 'time','log10(r[m])\flx','tflx','adv','dif','rain','rxn','res');
                        
                        for iz = 1: nz
                            ucvsld2 = 1d0 - poro(iz);
                            if (msldunit == 'blk'); ucvsld2 = 1d0; end
                            
                            fprintf(ipsd,[repmat('%e\t',1, nps + 2) '\n'], z(iz),mpsd(isps,1:nps,iz), time );
                            fprintf(ipsdv,[repmat('%e\t',1, nps + 2) '\n'], z(iz), 4d0/3d0*pi*(10d0.^ps(1:nps)).^3d0.*mpsd(isps,1:nps,iz)'.*dps(1:nps) ...
                                / (  msld(isps,iz)*mv(isps)*1d-6  ) * 1d2 ...
                                /ucvsld2  ...
                                , time );
                            if (~incld_rough)  
                                fprintf(ipsds,[repmat('%e\t',1, nps + 2) '\n'], z(iz), 4d0*pi*(10d0.^ps(1:nps)).^2d0.*mpsd(isps,1:nps,iz)'.*dps(1:nps) ...
                                    / ssa(isps,iz)  * 1d2 ...
                                    , time );
                            else
                                fprintf(ipsds,[repmat('%e\t',1, nps + 2) '\n'], z(iz), 4d0*pi*(10d0.^ps(1:nps)).^2d0 ...
                                    *rough_c0.*(10d0.^ps(1:nps)).^rough_c1.*mpsd(isps,1:nps,iz)'.*dps(1:nps) ...
                                    / ssa(isps,iz)  * 1d2 ...
                                    ,time );
                            end 
                        end
                        
                        for ips=1:nps
                            fprintf(ipsdflx,'%e\t%e\t', time,ps(ips));
                            for iflx=1:nflx_psd
                                fprintf(ipsdflx,'%e\t', sum(squeeze(flx_mpsd(isps,ips,iflx,:)).*dz(:)) );
                            end 
                            fprintf(ipsdflx,'\n');
                        end 
                        
                        fclose(ipsd);
                        fclose(ipsdv);
                        fclose(ipsds);
                        fclose(ipsdflx);
                    end 
                else 
                    ipsd = fopen(strcat(profdir,'/psd-',chr,'.txt'), 'w');
                    ipsdv = fopen(strcat(profdir,'/psd(v%)-',chr,'.txt'), 'w');
                    ipsds = fopen(strcat(profdir,'/psd(SA%)-',chr,'.txt'), 'w');
                    ipsdflx = fopen(strcat(flxdir,'/flx_psd-',chr,'.txt'), 'w');

                    fprintf(ipsd,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                    fprintf(ipsdv,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                    fprintf(ipsds,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',ps(1:nps),'time');
                    fprintf(ipsdflx,[repmat('%s\t',1, nflx_psd + 2) '\n'], 'time','log10(r[m])\flx','tflx','adv','dif','rain','rxn','res');
                    
                    for iz = 1: nz
                        ucvsld2 = 1d0 - poro(iz);
                        if (msldunit == 'blk'); ucvsld2 = 1d0; end 
                        
                        fprintf(ipsd,[repmat('%e\t',1, nps+2) '\n'], z(iz), (psd(1:nps,iz)), time );
                        fprintf(ipsdv,[repmat('%e\t',1, nps+2) '\n'], z(iz), (4d0/3d0*pi*(10d0.^ps(1:nps)).^3d0.*psd(1:nps,iz).*dps(1:nps) ...
                            / ( sum( msld(:,iz).*mv(:)*1d-6) + mblk(iz)*mvblk*1d-6 ) * 1d2 ...
                            /ucvsld2  ...
                            ), time );
                        if (~incld_rough)  
                            fprintf(ipsds,[repmat('%e\t',1, nps+2) '\n'], z(iz), (4d0*pi*(10d0.^ps(1:nps)).^2d0.*psd(1:nps,iz).*dps(1:nps) ...
                                / ssab(iz)  * 1d2 ...
                                ), time );
                        else
                            fprintf(ipsds,[repmat('%e\t',1, nps+2) '\n'], z(iz), (4d0*pi*(10d0.^ps(1:nps)).^2d0*rough_c0.*(10d0.^ps(1:nps)).^rough_c1.*psd(1:nps,iz).*dps(1:nps) ...
                                / ssab(iz)  * 1d2 ...
                                ), time );
                        end 
                    end
                    
                    for ips=1:nps
                        fprintf(ipsdflx,'%e\t%e\t', time,ps(ips));
                        for iflx=1:nflx_psd
                            fprintf(ipsdflx,'%e\t', sum(squeeze(flx_psd(ips,iflx,:)).*dz(:)) );
                        end 
                        fprintf(ipsdflx,'\n');
                    end 
                    
                    fclose(ipsd);
                    fclose(ipsdv);
                    fclose(ipsds);
                    fclose(ipsdflx);
                
                end 
            
            end 
            
            irec_prof=irec_prof+1;
            
            for isps=1:nsp_sld 
                fid = fopen(isldflx(isps),'a');
                fprintf(fid, '%e\t', time);
                for iflx=1:nflx
                    fprintf(fid,'%e\t', sum(squeeze(flx_sld(isps,iflx,:)).*dz(:)) );
                end 
                fprintf(fid,'\n');
                fclose(fid);
            end 
            
            for ispa=1:nsp_aq 
                fid = fopen(iaqflx(ispa),'a');
                fprintf(fid, '%e\t', time);
                for iflx=1:nflx
                    fprintf(fid,'%e\t', sum(squeeze(flx_aq(ispa,iflx,:)).*dz(:)) );
                end 
                fprintf(fid,'\n');
                fclose(fid);
            end 
            
            for ispg=1:nsp_gas 
                fid = fopen(igasflx(ispg),'a');
                fprintf(fid, '%e\t', time);
                for iflx=1:nflx
                    fprintf(fid,'%e\t', sum(squeeze(flx_gas(ispg,iflx,:)).*dz(:)) );
                end 
                fprintf(fid,'\n');
                fclose(fid);
            end 
            
            for ico2=1:6 
                fid = fopen(ico2flx(ico2),'a');
                fprintf(fid, '%e\t', time);
                if ico2 <= 4
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(ico2,iflx,:)).*dz(:)) );
                    end 
                elseif ico2 == 5
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(2,iflx,:)).*dz(:)) + sum(squeeze(flx_co2sp(3,iflx,:)).*dz(:)) ...
                            +  sum(squeeze(flx_co2sp(4,iflx,:)).*dz(:)) );
                    end 
                elseif ico2 == 6
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(3,iflx,:)).*dz(:)) + 2d0*sum(squeeze(flx_co2sp(4,iflx,:)).*dz(:)) );
                    end 
                end 
                fprintf(fid,'\n');
                fclose(fid);
            end 
            
            flx_recorded = true;
            
            isldprof = fopen(strcat(profdir,'/prof_sld-save.txt'), 'w');
            igasprof = fopen(strcat(profdir,'/prof_gas-save.txt'), 'w');
            iaqprof = fopen(strcat(profdir,'/prof_aq-save.txt'), 'w');
            ibsd = fopen(strcat(profdir,'/bsd-save.txt'), 'w');
            ipsd = fopen(strcat(profdir,'/psd-save.txt'), 'w');
            isa = fopen(strcat(profdir,'/sa-save.txt'), 'w');
                
            fprintf(isldprof,[repmat('%s\t',1, max(1,nsp_sld)+2) '\n'], 'z',(chrsld(1:nsp_sld)),'time');
            fprintf(iaqprof,[repmat('%s\t',1, max(1,nsp_aq)+3) '\n'], 'z',(chraq(1:nsp_aq)),'ph','time');
            fprintf(igasprof,[repmat('%s\t',1, max(1,nsp_gas)+2) '\n'], 'z',(chrgas(1:nsp_gas)),'time');
            fprintf(ibsd,[repmat('%s\t',1, 10) '\n'], 'z','poro', 'sat', 'v[m/yr]', 'm2/m3' ,'w[m/yr]'  ...
                , 'vol[m3/m3]','dens[g/cm3]', 'blk[wt%]', 'time');
            fprintf(ipsd,['%s\t' repmat('%f\t',1, nps) '%s\n'], 'z[m]\log10(r[m])',(ps(1:nps)),'time');
            fprintf(isa,[repmat('%s\t',1, max(1,nsp_sld)+2) '\n'], 'z',(chrsld(1:nsp_sld)),'time');

            for iz = 1: nz
                ucvsld1 = 1d0;
                if (msldunit == 'blk'); ucvsld1 = 1d0 - poro(iz); end 
                fprintf(isldprof,[repmat('%e\t',1, max(1,nsp_sld)+2) '\n'], z(iz),(msldx(1:nsp_sld,iz)),time);
                fprintf(igasprof,[repmat('%e\t',1, max(1,nsp_gas)+2) '\n'], z(iz),(mgasx(1:nsp_gas,iz)),time);
                fprintf(iaqprof,[repmat('%e\t',1, max(1,nsp_aq)+3) '\n'], z(iz),(maqx(1:nsp_aq,iz)),-log10(prox(iz)),time);
                fprintf(ibsd,[repmat('%e\t',1, 10) '\n'], z(iz), poro(iz),sat(iz),v(iz),hrb(iz),w(iz),sldvolfrac(iz),rho_grain_z(iz) ...
                    ,mblkx(iz)*mwtblk*1d2/ucvsld1/(rho_grain_z(iz)*1d6),time );
                fprintf(ipsd,[repmat('%e\t',1, nps+2) '\n'], z(iz), (psd(1:nps,iz)), time );
                fprintf(isa,[repmat('%e\t',1, max(1,nsp_sld)+2) '\n'], z(iz),(hr(1:nsp_sld,iz)),time );
            end

            fclose(isldprof);
            fclose(iaqprof);
            fclose(igasprof);
            fclose(ibsd);
            fclose(ipsd);
            fclose(isa);
            
            if (display_lim_in); display_lim = false; end 
            
        end    
        
        
        % saving flx when climate is changed within model 
        if ( (any(climate) && any (ict_change))  ... 
            ||(time>=rectime_flx(irec_flx+1)) )  
            
            if (time>=rectime_flx(irec_flx+1))  
                irec_flx = irec_flx + 1
            end 
            
            if (~ flx_recorded)  
                for isps=1:nsp_sld 
                    fid = fopen(isldflx(isps),'a');
                    fprintf(fid, '%e\t', time);
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_sld(isps,iflx,:)).*dz(:)) );
                    end 
                    fprintf(fid,'\n');
                    fclose(fid);
                end 
                
                for ispa=1:nsp_aq 
                    fid = fopen(iaqflx(ispa),'a');
                    fprintf(fid, '%e\t', time);
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_aq(ispa,iflx,:)).*dz(:)) );
                    end 
                    fprintf(fid,'\n');
                    fclose(fid);
                end 
                
                for ispg=1:nsp_gas 
                    fid = fopen(igasflx(ispg),'a');
                    fprintf(fid, '%e\t', time);
                    for iflx=1:nflx
                        fprintf(fid,'%e\t', sum(squeeze(flx_gas(ispg,iflx,:)).*dz(:)) );
                    end 
                    fprintf(fid,'\n');
                    fclose(fid);
                end 
                
                for ico2=1:6 
                    fid = fopen(ico2flx(ico2),'a');
                    fprintf(fid, '%e\t', time);
                    if ico2 <= 4
                        for iflx=1:nflx
                            fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(ico2,iflx,:)).*dz(:)) );
                        end 
                    elseif ico2 == 5
                        for iflx=1:nflx
                            fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(2,iflx,:)).*dz(:)) + sum(squeeze(flx_co2sp(3,iflx,:)).*dz(:)) ...
                                +  sum(squeeze(flx_co2sp(4,iflx,:)).*dz(:)) );
                        end 
                    elseif ico2 == 6
                        for iflx=1:nflx
                            fprintf(fid,'%e\t', sum(squeeze(flx_co2sp(3,iflx,:)).*dz(:)) + 2d0*sum(squeeze(flx_co2sp(4,iflx,:)).*dz(:)) );
                        end 
                    end 
                    fprintf(fid,'\n');
                    fclose(fid);
                end 
                
            end 
            
        end 
        
        % break
        
        
        
        
        
        
        
        it = it + 1;
        time = time + dt;
        count_dtunchanged = count_dtunchanged + 1;
        
        progress_rate_prev = progress_rate;
        
        % call cpu_time(time_fin)
        % call system_clock(t2,t_rate,t_max)
        t2 = cputime;
        % if ( t2 < t1 ) 
            % diff = (t_max - t1) + t2 + 1
        % else
            % diff = t2 - t1
        % end
        
        % progress_rate = dt/(time_fin-time_start)*sec2yr % (model yr)/(computer yr)
        % progress_rate = (time_fin-time_start) % (computer sec)
        % progress_rate = diff/dble(t_rate) % (computer sec)
        progress_rate = double(t2-t1); % (computer sec)
        progress_rate = max(progress_rate,1d-2);
        
        if (~timestep_fixed)  
            if (it~=1)  
                if (flgreducedt)  
                    maxdt = maxdt/10d0;
                    flgreducedt = false;
                    count_dtunchanged = 0;
                else
                    maxdt = maxdt* (progress_rate/progress_rate_prev)^(-0.33d0);
                    if (maxdt > maxdt_max); maxdt = maxdt_max; end
                    if (dt < maxdt); count_dtunchanged = 0; end
                end 
                
                if (count_dtunchanged > count_dtunchanged_Max)  
                    maxdt = maxdt*10d0;
                    count_dtunchanged = 0;
                end 
            end 
        end 
        
        % if (progress_rate ==0d0 || progress_rate_prev ==0d0) maxdt = 1d2
        % print *,progress_rate,progress_rate_prev,maxdt,time_fin,time_start
        if (isnan(maxdt)||maxdt ==0d0)  
            error( 'maxdt is nan or zero: progress_rate = \t%7.6E\n',progress_rate)
        end 
        
        if (display  && (~ display_lim))  
            fprintf('\n');
            fprintf('computation time per iteration = \t%7.6E\t [sec]\n',progress_rate)
            fprintf('maxdt = \t%7.6E\t [yr]\n',maxdt)
            fprintf('count_dtunchanged = \t%d\n',count_dtunchanged)
            fprintf('-----------------------------------------\n');
            fprintf('\n');
        end 

    end

end


function [nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc] = get_variables_num
    
    file_name = './slds.in';
    A = importdata(file_name);
    nsp_sld = size(A,1);
    % chrsld = A(2:end);
    % call Console4(file_name,nsp_sld)
    file_name = './solutes.in';
    A = importdata(file_name);
    nsp_aq = size(A,1);
    % chraq = A(2:end);
    % call Console4(file_name,nsp_aq)
    file_name = './gases.in';
    A = importdata(file_name);
    nsp_gas = size(A,1);
    % chrgas = A(2:end);
    % call Console4(file_name,nsp_gas)
    file_name = './extrxns.in';
    A = importdata(file_name);
    nrxn_ext = size(A,1);
    % chrrxn_ext = A(2:end);
    % call Console4(file_name,nrxn_ext)
    file_name = './kinspc.in';
    A = importdata(file_name);
    nsld_kinspc = size(A,1);
    % chrsld_kinspc = A(2:end);
    % call Console4(file_name,nsld_kinspc)

    nsp_sld = nsp_sld - 1;
    nsp_aq = nsp_aq - 1;
    nsp_gas = nsp_gas - 1;
    nrxn_ext = nrxn_ext - 1;
    nsld_kinspc = nsld_kinspc - 1;

end


function [chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc] ... 
    = get_variables(nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc )
        
    if (nsp_aq>=1)
        file_name = './solutes.in';
        A = importdata(file_name);
        chraq = A(2:end);
        chraq = string(chraq);
    else 
        chraq = "";
    end 

    if (nsp_sld>=1)  
        file_name = './slds.in';
        A = importdata(file_name);
        chrsld = A(2:end);
        chrsld = string(chrsld);
    else
        chrsld = "";
    end 

    if (nsp_gas>=1)
        file_name = './gases.in';
        A = importdata(file_name);
        chrgas = A(2:end);
        chrgas = string(chrgas);
    else
        chrgas = "";
    end 

    if (nrxn_ext>=1)
        file_name = './extrxns.in';
        A = importdata(file_name);
        chrrxn_ext = A(2:end);
        chrrxn_ext = string(chrrxn_ext);
    else
        chrrxn_ext = "";
    end 

    if (nsld_kinspc>=1)
        file_name = './kinspc.in';
        A = importdata(file_name);
        chrsld_kinspc = A(2:end,1);
        chrsld_kinspc = string(chrsld_kinspc);
        kin_sld_spc = A(2:end,2);
    else
        chrsld_kinspc = "";
        kin_sld_spc = nan;
    end 

end


function [nz,ztot,ttot,rainpowder,zsupp,poroi,satup,zsat,zml_ref,w,qin,p80,sim_name,plant_rain,runname_save ...
    ,count_dtunchanged_Max,tc] = get_bsdvalues
    
    file_name = './frame.in';

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        if cnt==1
            ztot = str2num(char(a(1)));
        elseif cnt==2
            nz = int32(str2num(char(a(1))));
        elseif cnt==3
            ttot = (str2num(char(a(1))));
        elseif cnt==4
            tc = (str2num(char(a(1))));
        elseif cnt==5
            rainpowder = (str2num(char(a(1))));
        elseif cnt==6
            plant_rain = (str2num(char(a(1))));
        elseif cnt==7
            zsupp = (str2num(char(a(1))));
        elseif cnt==8
            poroi = (str2num(char(a(1))));
        elseif cnt==9
            satup = (str2num(char(a(1))));
        elseif cnt==10
            zsat = (str2num(char(a(1))));
        elseif cnt==11
            zml_ref = (str2num(char(a(1))));
        elseif cnt==12
            w = (str2num(char(a(1))));
        elseif cnt==13
            qin = (str2num(char(a(1))));
        elseif cnt==14
            p80 = (str2num(char(a(1))));
        elseif cnt==15
            count_dtunchanged_Max = int32(str2num(char(a(1))));
        elseif cnt==16
            runname_save = char(a(1));
        elseif cnt==18
            sim_name = char(a(1));
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);

end


function [rain_all] = get_rainwater(nsp_aq_all,chraq_all,def_rain)

    file_name = './rain.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end

    % in default 
    rain_all = zeros(nsp_aq_all,1,'double');
    rain_all(:) = def_rain;

    if (n_tmp <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        chr_tmp = (char(a(2)));
        val_tmp = str2num(char(a(2)));
        if (any(chraq_all == chr_tmp))  
            rain_all(find(chraq_all==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);
end



function [parentrock_frct_all] = get_parentrock(nsp_sld_all,chrsld_all,def_pr)

    file_name = './parentrock.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end


    % in default 
    parentrock_frct_all = zeros(nsp_sld_all,1,'double');
    parentrock_frct_all(:) = def_pr;

    if (n_tmp <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        % disp(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        chr_tmp = (char(a(1)));
        val_tmp = str2num(char(a(2)));
        if (any(chrsld_all == chr_tmp))  
            parentrock_frct_all(find(chrsld_all==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);
end


function [atm_all] = get_atm(nsp_gas_all,chrgas_all)

    file_name = './atm.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end


    % in default 
    atm_all = zeros(nsp_gas_all,1,'double');
    atm_all(find(chrgas_all=='po2')) = 0.21d0;
    atm_all(find(chrgas_all=='pco2')) = 10d0^(-3.5d0);
    atm_all(find(chrgas_all=='pnh3')) = 1d-9;
    atm_all(find(chrgas_all=='pn2o')) = 270d-9;

    if (n_tmp <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        % disp(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        chr_tmp = (char(a(1)));
        val_tmp = str2num(char(a(2)));
        if (any(chrgas_all == chr_tmp))  
            atm_all(find(chrgas_all==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);
end


function [dust_frct_all] = get_dust(nsp_sld_all,chrsld_all,def_dust)

    file_name = './dust.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end


    % in default 
    dust_frct_all = zeros(nsp_sld_all,1,'double');
    dust_frct_all(:) = def_dust;

    if (n_tmp <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        % disp(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        chr_tmp = (char(a(1)));
        val_tmp = str2num(char(a(2)));
        if (any(chrsld_all == chr_tmp))  
            dust_frct_all(find(chrsld_all==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);
end


function [OM_frct_all] = get_OM_rain(nsp_sld_all,chrsld_all,def_OM_frc)

    file_name = './OM_rain.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end

    % in default 
    OM_frct_all = zeros(nsp_sld_all,1,'double');
    OM_frct_all(:) = def_OM_frc;

    if (n_tmp <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        % disp(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        chr_tmp = (char(a(1)));
        val_tmp = str2num(char(a(2)));
        if (any(chrsld_all == chr_tmp))  
            OM_frct_all(find(chrsld_all==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);

end


function [iwtype,imixtype,poroiter_in,display,display_lim_in,read_data,incld_rough ...
    ,al_inhibit,timestep_fixed,method_precalc,regular_grid,sld_enforce ...% inout
    ,poroevol,surfevol1,surfevol2,do_psd,lim_minsld_in,do_psd_full,season ...% inout
    ] = get_switches
    
    file_name = './switches.in';
    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        if cnt==1
            iwtype = int32(str2num(char(a(1))));
        elseif cnt==2
            imixtype = int32(str2num(char(a(1))));
        elseif cnt==3
            poroiter_in = (str2num(char(a(1))));
        elseif cnt==4
            lim_minsld_in = (str2num(char(a(1))));
        elseif cnt==5
            display = (str2num(char(a(1))));
        elseif cnt==6
            display_lim_in = (str2num(char(a(1))));
        elseif cnt==7
            read_data = (str2num(char(a(1))));
        elseif cnt==8
            incld_rough = (str2num(char(a(1))));
        elseif cnt==9
            al_inhibit = (str2num(char(a(1))));
        elseif cnt==10
            timestep_fixed = (str2num(char(a(1))));
        elseif cnt==11
            method_precalc = (str2num(char(a(1))));
        elseif cnt==12
            regular_grid = (str2num(char(a(1))));
        elseif cnt==13
            sld_enforce = (str2num(char(a(1))));
        elseif cnt==14
            poroevol = (str2num(char(a(1))));
        elseif cnt==15
            surfevol1 = (str2num(char(a(1))));
        elseif cnt==16
            surfevol2 = (str2num(char(a(1))));
        elseif cnt==17
            do_psd = (str2num(char(a(1))));
        elseif cnt==18
            do_psd_full = (str2num(char(a(1))));
        elseif cnt==19
            season = (str2num(char(a(1))));
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);

end


function [n_file] = get_clim_num(file_in)

    file_name = file_in; 
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end

end


function [dz,z] = makegrid(beta,nz,ztot,regular_grid)  %  making grid, after Hoffmann and Chiang, 2000

    z = zeros(nz,1,'double');
    dz = zeros(nz,1,'double');
    for iz = 1: nz 
        z(iz) = double(iz)*ztot/double(nz);  % regular grid 
        if (iz==1) 
            dz(iz) = ztot*log((beta+(z(iz)/ztot)^2d0)/(beta-(z(iz)/ztot)^2d0))/log((beta+1d0)/(beta-1d0));
        end
        if (iz~=1)  
            dz(iz) = ztot*log((beta+(z(iz)/ztot)^2d0)/(beta-(z(iz)/ztot)^2d0))/log((beta+1d0)/(beta-1d0)) - sum(dz(1:iz-1));
        end
    end

    if (regular_grid)  
        dz(:) = ztot/double(nz);  % when implementing regular grid
    end 

    for iz=1:nz  % depth is defined at the middle of individual layers 
        if (iz==1); z(iz)=dz(iz)*0.5d0 ; end
        if (iz~=1); z(iz) = z(iz-1)+dz(iz-1)*0.5d0 + 0.5d0*dz(iz); end
    end

end


function [n_tmp] = get_sa_num

    file_name = './sa.in';
    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end
end


function [hrii,chrsld_sa_dum] = get_sa(nsp_sld,chrsld,def_hr,nsld_sa)

    file_name = './sa.in';
    % in default 
    chrsld_sa_dum = strings(nsld_sa,1);
    hrii = zeros(nsp_sld,1,'double');
    hrii(:) = def_hr;

    if (nsld_sa <= 0); return; end

    fid = fopen(file_name);
    tline = fgetl(fid);
    cnt = 0;
    while ischar(tline)
        % disp(tline)
        a = strsplit(tline);
        tline = fgetl(fid);
        
        if cnt ==0
            cnt = cnt + 1;
            continue
        end
        
        chr_tmp = (char(a(1)));
        val_tmp = str2num(char(a(2)));
        chrsld_sa_dum(cnt) = chr_tmp; 
        
        if (any(chrsld == chr_tmp))  
            hrii(find(chrsld==chr_tmp)) = val_tmp;
        end
        
        cnt = cnt + 1;
            
    end
    fclose(fid);
end


function [up,dwn,cnr,adf] = calcupwindscheme(w,nz)

    % copied and pasted from iMP code and modified for weathering (to be implemented)

    % ------------ determine variables to realize advection 
    %  upwind scheme 
    %  up  ---- burial advection at grid i = sporo(i)*w(i)*(some conc. at i) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
    %  dwn ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i)*w(i)*(some conc. at i) 
    %  cnr ---- burial advection at grid i = sporo(i+1)*w(i+1)*(some conc. at i+1) - sporo(i-1)*w(i-1)*(some conc. at i - 1) 
    %  when burial rate is positive, scheme need to choose up, i.e., up = 1.  
    %  when burial rate is negative, scheme need to choose dwn, i.e., dwn = 1.  
    %  where burial change from positive to negative or vice versa, scheme chooses cnr, i.e., cnr = 1. for the mass balance sake 

    up = zeros(nz,1,'double');
    dwn= zeros(nz,1,'double');
    cnr = zeros(nz,1,'double');
    adf= ones(nz,1,'double');
    for iz=1:nz 
        if (iz==1)  
            if (w(iz)>=0d0 && w(iz+1)>=0d0)   % positive burial 
                up(iz) = 1;
            elseif (w(iz)<=0d0 && w(iz+1)<=0d0)   % negative burial 
                dwn(iz) = 1;
            else   %  where burial sign changes  
                if (~(w(iz)*w(iz+1) <=0d0))  
                    error('error')
                end
                cnr(iz) = 1;
            end
        elseif (iz==nz)  
            if (w(iz)>=0d0 && w(iz-1)>=0d0) 
                up(iz) = 1;
            elseif (w(iz)<=0d0 && w(iz-1)<=0d0) 
                dwn(iz) = 1;
            else 
                if (~(w(iz)*w(iz-1) <=0d0))  
                    error('error')
                end
                cnr(iz) = 1;
            end
        else 
            % if iz-1 and iz+1 have the same sign,  it can be assigned either as up or dwn
            % else cnr whose neighbor has a different sign 
            if (w(iz) >=0d0)  
                if (w(iz+1)>=0d0 && w(iz-1)>=0d0) 
                    up(iz) = 1;
                else
                    cnr(iz) = 1;
                end
            else  
                if (w(iz+1)<=0d0 && w(iz-1)<=0d0) 
                    dwn(iz) = 1;
                else
                    cnr(iz) = 1;
                end
            end
        end
    end        

    if (sum(up(:)+dwn(:)+cnr(:))~=nz) 
        error('error %4.3e %4.3e %4.3e\n',sum(up),sum(dwn),sum(cnr))
    end

    % try to make sure mass balance where advection direction changes 
    % 
    % case (i)
    %       :           w         direction    
    %     iz - 2        +             ^              w(iz-1) - w(iz-2)
    %     iz - 1        +             ^              w(iz  ) - w(iz-1)
    %     iz            +             ^            a[w(iz+1) - w(iz  )] + b[w(iz+1) - w(iz-1)]  
    %     iz + 1        -             v            c[w(iz+1) - w(iz  )] + d[w(iz+2) - w(iz  )]
    %     iz + 2        -             v              w(iz+2) - w(iz+1)
    %     iz + 3        -             v              w(iz+3) - w(iz+2) 
    % layers [iz] & [iz+1] must yield [w(iz+1) - w(iz  )]
    % and calculated as (a+b+c)w(iz+1) - (a+c+d)w(iz  ) - b w(iz-1) + d w(iz+1) 
    % thus b = d = 0 and   a + b + c = 1 and a + c + d = 1
    % a and c can be arbitrary as long as satisfying a + c = 1 
    % --------------------------------------------------------------------------------------------
    % case (ii)
    %       :           w         direction    
    %     iz - 2        -             v              w(iz-2) - w(iz-3)
    %     iz - 1        -             v              w(iz-1) - w(iz-2)
    %     iz            -             v            a[w(iz  ) - w(iz-1)] + b[w(iz+1) - w(iz-1)]  
    %     iz + 1        +             ^            c[w(iz+2) - w(iz+1)] + d[w(iz+2) - w(iz  )]
    %     iz + 2        +             ^              w(iz+3) - w(iz+2)
    %     iz + 3        +             ^              w(iz+4) - w(iz+3) 
    % layers [iz] & [iz+1] must yield [w(iz+2) - w(iz-1)]
    % and calculated as (c+d)w(iz+2) - (a+b)w(iz-1) + (a-d)w(iz  ) + (b-c)w(iz+1) 
    % thus c + d = 1, a + b = 1, a - d = 0, and b - c = 0
    % these can be satisfied by b = c = 1 - a and d = a and a can be arbitrary as long as 0 <= a <= 1
    cnr_save = cnr; 
    for iz=1:nz-1
        if (cnr_save(iz)==1 && cnr_save(iz+1)==1)  
        % if (cnr(iz)==1 && cnr(iz+1)==1)  
            if (w(iz) < 0d0 && w(iz+1) >= 0d0) 
                corrf = 5d0;  %  This assignment of central advection term helps conversion especially when assuming turbo2 mixing 
                cnr(iz+1)=abs(w(iz)^corrf)/(abs(w(iz+1)^corrf)+abs(w(iz)^corrf));
                cnr(iz)=abs(w(iz+1)^corrf)/(abs(w(iz+1)^corrf)+abs(w(iz)^corrf));
                dwn(iz+1)=1d0-cnr(iz+1);
                up(iz)=1d0-cnr(iz);
            end 
        end 
        if (cnr_save(iz)==1 && cnr_save(iz+1)==1)  
        % if (cnr(iz)==1 && cnr(iz+1)==1)  
            if (w(iz)>= 0d0 && w(iz+1) < 0d0) 
                cnr(iz+1)=0;
                cnr(iz)=0;
                up(iz+1)=1;
                dwn(iz)=1;
                adf(iz)=abs(w(iz+1))/(abs(w(iz+1))+abs(w(iz)));
                adf(iz+1)=abs(w(iz))/(abs(w(iz+1))+abs(w(iz)));
            end 
        end 
    end       

end


function [psd_pr] = calc_psd_pr( ...
    nps ...% input
    ,pi,p80,ps_sigma_std,poroi,volsld,tol ...% input
    ,ps,dps ...% input
    ,msldunit ...% input
    )


    psu_pr = log10(p80);
    pssigma_pr = ps_sigma_std;

    % calculate parent rock particle size distribution 
    psd_pr = zeros(nps,1,'double');
    psd_pr(:) = 1d0/pssigma_pr/sqrt(2d0*pi)*exp( -0.5d0*( (ps(:) - psu_pr)/pssigma_pr ).^2d0 );

    % to ensure sum is 1
    % print *, sum(psd_pr*dps)
    psd_pr(:) = psd_pr(:)/sum(psd_pr(:).*dps(:));  
    % print *, sum(psd_pr*dps)
    % stop

    % balance for volumes
    % sum(msldi*mv*1d-6) (m3/m3) must be equal to sum( 4/3(pi)r3 * psd_pr * dps) 
    % where psd is number / bulk m3 / log r 
    % (if msld is defined as mol/sld m3 then msldi needs to be multiplied by (1 - poro)
    % volsld = sum(msldi*mv*1d-6) + mblki*mvblk*1d-6 % bulk case
    psd_pr(:) = psd_pr(:)*( volsld )  / sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_pr(:).*dps(:));

    if ( msldunit == 'sld')  
        psd_pr = psd_pr*(1d0-poroi)

        if ( abs( ( ( volsld ) * (1d0-poroi) -sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_pr(:).*dps(:))) ...
            / ( ( volsld ) * (1d0-poroi) )  ) > tol)  
            error( '%4.3e %4.3e \n',( volsld )* (1d0-poroi) , sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_pr(:).*dps(:)) );
        end
    elseif ( msldunit == 'blk')  

        if ( abs( ( ( volsld ) -sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_pr(:).*dps(:))) / ( ( volsld ) )  ) > tol)  
            error( '%4.3e %4.3e \n',( volsld ) , sum(4d0/3d0*pi*(10d0.^ps(:)).^3d0.*psd_pr(:).*dps(:)) );
        end
    end

end


function res = k_arrhenius(kref,tempkref,tempk,eapp,rg)
    res = kref*exp(-eapp/rg*(1d0/tempk-1d0/tempkref));
end


function [mgasx_loc] = get_mgasx_all( ...
    nz,nsp_gas_all,nsp_gas,nsp_gas_cnst ...
    ,chrgas,chrgas_all,chrgas_cnst ...
    ,mgasx,mgasc ...
    )

    mgasx_loc=zeros(nsp_gas_all,nz,'double');

    for ispg = 1: nsp_gas_all
        if (any(chrgas==chrgas_all(ispg)))  
            mgasx_loc(ispg,:) =  mgasx(find(chrgas==chrgas_all(ispg)),:);
        elseif (any(chrgas_cnst==chrgas_all(ispg)))  
            mgasx_loc(ispg,:) =  mgasc(find(chrgas_cnst==chrgas_all(ispg)),:);
        end 
    end
    
end


function [kin,dkin_dmsp] = sld_kin( ...
    nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv_tmp ...% input
    ,nsp_gas_all,chrgas_all,mgas_loc ...% input
    ,mineral,dev_sp ...% input 
    ) 

    cal2j = 4.184d0 ;

    % initialize output variables
    kin = zeros(nz,1,'double');
    dkin_dmsp = zeros(nz,1,'double');

    % local variables
    pco2 = zeros(nz,1,'double');
    mco2 = 0;kinco2_ref=0;eaco2=0;
    mh=0;moh=0;kinn_ref=0;kinh_ref=0;kinoh_ref=0;ean=0;eah=0;eaoh=0;tc_ref=0;

    switch(mineral)
        case('ka')
            mh = 0.777d0;
            moh = -0.472d0;
            kinn_ref = 10d0^(-13.18d0)*sec2yr;
            kinh_ref = 10d0^(-11.31d0)*sec2yr;
            kinoh_ref = 10d0^(-17.05d0)*sec2yr;
            ean = 22.2d0;
            eah = 65.9d0;
            eaoh = 17.9d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case('ab')
            mh = 0.457d0;
            moh = -0.572d0;
            kinn_ref = 10d0^(-12.56d0)*sec2yr;
            kinh_ref = 10d0^(-10.16d0)*sec2yr;
            kinoh_ref = 10d0^(-15.6d0)*sec2yr;
            ean = 69.8d0;
            eah = 65d0;
            eaoh = 71d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('kfs')
            mh = 0.5d0;
            moh = -0.823d0;
            kinn_ref = 10d0^(-12.41d0)*sec2yr;
            kinh_ref = 10d0^(-10.06d0)*sec2yr;
            kinoh_ref = 10d0^(-9.68d0)*sec2yr*kw^(-moh);
            ean = 9.08*cal2j;
            eah = 12.4d0*cal2j;
            eaoh = 22.5d0*cal2j;
            tc_ref = 25d0;
            % from Brantley et al 2008
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('fo')
            mh = 0.47d0;
            moh = 0d0;
            kinn_ref = 10d0^(-10.64d0)*sec2yr;
            kinh_ref = 10d0^(-6.85d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 79d0;
            eah = 67.2d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('fa')
            mh = 1d0;
            moh = 0d0;
            kinn_ref = 10d0^(-12.80d0)*sec2yr;
            kinh_ref = 10d0^(-4.80d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 94.4d0;
            eah = 94.4d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('an')
            mh = 1.411d0;
            moh = 0d0;
            kinn_ref = 10d0^(-9.12d0)*sec2yr;
            kinh_ref = 10d0^(-3.5d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 17.8d0;
            eah = 16.6d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('la')
            mh = 0.626d0;
            moh = 0d0;
            kinn_ref = 10d0^(-10.91d0)*sec2yr;
            kinh_ref = 10d0^(-7.87d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 45.2d0;
            eah = 42.1d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('and')
            mh = 0.541d0;
            moh = 0d0;
            kinn_ref = 10d0^(-11.47d0)*sec2yr;
            kinh_ref = 10d0^(-8.88d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 57.4d0;
            eah = 53.5d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('olg')
            mh = 0.457d0;
            moh = 0d0;
            kinn_ref = 10d0^(-11.84d0)*sec2yr;
            kinh_ref = 10d0^(-9.67d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 69.8d0;
            eah = 65.0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('by')
            mh = 1.018d0;
            moh = 0d0;
            kinn_ref = 10d0^(-9.82d0)*sec2yr;
            kinh_ref = 10d0^(-5.85d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 31.5d0;
            eah = 29.3d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('cc')
            mh = 1d0;
            moh = 0d0;
            kinn_ref = 10d0^(-5.81d0)*sec2yr;
            kinh_ref = 10d0^(-0.3d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 23.5d0;
            eah = 14.4d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % adding co2 mechanism
            mco2 = 1d0;
            kinco2_ref = 10d0^(-3.48d0)*sec2yr;
            eaco2 = 35.4d0;
            pco2(:) = mgas_loc(find(chrgas_all=='pco2'),:);
            % from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                + pco2(:).^mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                case('pco2')
                    dkin_dmsp(:) = ( ... 
                        + mco2*pco2(:).^(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('arg')
            % assumed to be the same as those for cc
            mh = 1d0;
            moh = 0d0;
            kinn_ref = 10d0^(-5.81d0)*sec2yr;
            kinh_ref = 10d0^(-0.3d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 23.5d0;
            eah = 14.4d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % adding co2 mechanism
            mco2 = 1d0;
            kinco2_ref = 10d0^(-3.48d0)*sec2yr;
            eaco2 = 35.4d0;
            pco2(:) = mgas_loc(find(chrgas_all=='pco2'),:);
            % from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                + pco2(:).^mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                case('pco2')
                    dkin_dmsp(:) = ( ... 
                        + mco2*pco2(:).^(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('dlm') % for disordered dolomite
            mh = 0.500d0;
            moh = 0d0;
            kinn_ref = 10d0^(-7.53d0)*sec2yr;
            kinh_ref = 10d0^(-3.19d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 52.2d0;
            eah = 36.1d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % adding co2 mechanism
            mco2 = 0.5d0;
            kinco2_ref = 10d0^(-5.11d0)*sec2yr;
            eaco2 = 34.8d0;
            pco2(:) = mgas_loc(find(chrgas_all=='pco2'),:);
            % from Palandri and Kharaka, 2004 (excluding carbonate mechanism)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                + pco2(:).^mco2*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                case('pco2')
                    dkin_dmsp(:) = ( ... 
                        + mco2*pco2(:).^(mco2-1d0)*k_arrhenius(kinco2_ref,tc_ref+tempk_0,tc+tempk_0,eaco2,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('gb')
            mh = 0.992d0;
            moh = -0.784d0;
            kinn_ref = 10d0^(-11.50d0)*sec2yr;
            kinh_ref = 10d0^(-7.65d0)*sec2yr;
            kinoh_ref = 10d0^(-16.65d0)*sec2yr;
            ean = 61.2d0;
            eah = 47.5d0;
            eaoh = 80.1d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('amsi')
            mh = 0d0;
            moh = 0d0;
            kinn_ref = 10d0^(-12.23d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 0d0;
            ean = 74.5d0;
            eah = 0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            dkin_dmsp(:) = 0d0;

        case('qtz')
            mh = 0d0;
            moh = 0d0;
            kinn_ref = 10d0^(-13.40d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 0d0;
            ean = 90.9d0;
            eah = 0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            dkin_dmsp(:) = 0d0;

        case('gt')
            mh = 0d0;
            moh = 0d0;
            kinn_ref = 10d0^(-7.94d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 0d0;
            ean = 86.5d0;
            eah = 0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            dkin_dmsp(:) = 0d0;

        case('hm')
            mh = 1d0;
            moh = 0d0;
            kinn_ref = 10d0^(-14.60d0)*sec2yr;
            kinh_ref = 10d0^(-9.39d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 66.2d0;
            eah = 66.2d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            dkin_dmsp(:) = 0d0;

        case('ct')
            mh = 0d0;
            moh = -0.23d0;
            kinn_ref = 10d0^(-12d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 10d0^(-13.58d0)*sec2yr;
            ean = 73.5d0;
            eah = 0d0;
            eaoh = 73.5d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('mscv')
            mh = 0.370d0;
            moh = -0.22d0;
            kinn_ref = 10d0^(-13.55d0)*sec2yr;
            kinh_ref = 10d0^(-11.85d0)*sec2yr;
            kinoh_ref = 10d0^(-13.55d0)*sec2yr;
            ean = 22d0;
            eah = 22d0;
            eaoh = 22d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('plgp')
            mh = 0d0;
            moh = 0d0;
            kinn_ref = 10d0^(-12.4d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 0d0;
            ean = 29d0;
            eah = 0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case{'cabd','ill','kbd','nabd','mgbd'} % illite kinetics is assumed to be the same as smectite (Bibi et al., 2011)
            mh = 0.34d0;
            moh = -0.4d0;
            kinn_ref = 10d0^(-12.78d0)*sec2yr;
            kinh_ref = 10d0^(-10.98d0)*sec2yr;
            kinoh_ref = 10d0^(-16.52d0)*sec2yr;
            ean = 35d0;
            eah = 23.6d0;
            eaoh = 58.9d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case{'nph','anl'} % analcime kinetics is assumed to be the same as nepherine (cf. Ragnarsdottir, GCA, 1993)
            mh = 1.130d0;
            moh = -0.200d0;
            kinn_ref = 10d0^(-8.56d0)*sec2yr;
            kinh_ref = 10d0^(-2.73d0)*sec2yr;
            kinoh_ref = 10d0^(-10.76d0)*sec2yr;
            ean = 65.4d0;
            eah = 62.9d0;
            eaoh = 37.8d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 

        case('dp')
            mh = 0.71d0;
            moh = 0d0;
            kinn_ref = 10d0^(-11.11d0)*sec2yr;
            kinh_ref = 10d0^(-6.36d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 50.6d0;
            eah = 96.1d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case{'hb','cpx','agt'}
            mh = 0.70d0;
            moh = 0d0;
            kinn_ref = 10d0^(-11.97d0)*sec2yr;
            kinh_ref = 10d0^(-6.82d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 78.0d0;
            eah = 78.0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % for augite from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case{'en','opx','fer'}
            mh = 0.60d0;
            moh = 0d0;
            kinn_ref = 10d0^(-12.72d0)*sec2yr;
            kinh_ref = 10d0^(-9.02d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 80.0d0;
            eah = 80.0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % for enstatite from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case('tm')
            mh = 0.70d0;
            moh = 0d0;
            kinn_ref = 10d0^(-10.60d0)*sec2yr;
            kinh_ref = 10d0^(-8.40d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 94.4d0;
            eah = 18.9d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case('antp')
            mh = 0.440d0;
            moh = 0d0;
            kinn_ref = 10d0^(-14.24d0)*sec2yr;
            kinh_ref = 10d0^(-11.94d0)*sec2yr;
            kinoh_ref = 0d0;
            ean = 51.0d0;
            eah = 51.0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
        
        case('gps')
            mh = 0d0;
            moh = 0d0;
            kinn_ref = 10d0^(-2.79d0)*sec2yr;
            kinh_ref = 0d0;
            kinoh_ref = 0d0;
            ean = 0d0;
            eah = 0d0;
            eaoh = 0d0;
            tc_ref = 25d0;
            % from Palandri and Kharaka, 2004
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
            
        case('py')
            mh = 0d0;
            moh = -0.11d0;
            kinn_ref = 0d0;
            kinh_ref = 0d0;
            kinoh_ref = 10.0d0^(-8.19d0)*sec2yr*kho^0.5d0;
            ean = 0d0;
            eah = 0d0;
            eaoh = 57d0;
            tc_ref = 15d0;
            % from Williamson and Rimstidt (1994)
            kin(:) = ( ... 
                k_arrhenius(kinn_ref,tc_ref+tempk_0,tc+tempk_0,ean,rg) ...
                + prox(:).^mh*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                + prox(:).^moh*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                ); 
            switch(dev_sp)
                case('pro')
                    dkin_dmsp(:) = ( ... 
                        + mh*prox(:).^(mh-1d0)*k_arrhenius(kinh_ref,tc_ref+tempk_0,tc+tempk_0,eah,rg) ...
                        + moh*prox(:).^(moh-1d0)*k_arrhenius(kinoh_ref,tc_ref+tempk_0,tc+tempk_0,eaoh,rg) ...
                        ); 
                otherwise 
                    dkin_dmsp(:) = 0d0;
            end 
            
        case('g1')
            kin(:) = ( ...
                1d0/1d0 ...% mol m^-2 yr^-1, just a value assumed; turnover time of 1 year as in Chen et al. (2010, AFM) 
                );
            dkin_dmsp(:) = 0d0;
            
        case('g2')
            kin(:) = ( ...
                1d0/8d0 ...% mol m^-2 yr^-1, just a value assumed; turnover time of 8 year as in Chen et al. (2010, AFM) 
                );
            dkin_dmsp(:) = 0d0;
            
        case('g3')
            kin(:) = ( ...
                1d0/1d3 ...% mol m^-2 yr^-1, just a value assumed; picked up to represent turnover time of 1k year  
                );
            dkin_dmsp(:) = 0d0;
            
        otherwise 
            kin(:) =0d0;
            dkin_dmsp(:) = 0d0;

    end  

end


function [therm] = sld_therm( ...
    rg,tc,tempk_0,ss_x,ss_y ...% input
    ,mineral ...% input
    ) 

    cal2j = 4.184d0;

    % initialize output variables
    therm = 0d0;

    % local variables
    tc_ref=0;ha=0;therm_ref=0;delG=0;
    tc_ref_1=0;ha_1=0;therm_ref_1=0;therm_1=0;delG_1=0;
    tc_ref_2=0;ha_2=0;therm_ref_2=0;therm_2=0;delG_2=0;
    tc_ref_3=0;ha_3=0;therm_ref_3=0;therm_3=0;delG_3=0;
    tc_ref_4=0;ha_4=0;therm_ref_4=0;therm_4=0;delG_4=0;
    tc_ref_5=0;ha_5=0;therm_ref_5=0;therm_5=0;delG_5=0;
    tc_ref_6=0;ha_6=0;therm_ref_6=0;therm_6=0;delG_6=0;

    switch (mineral) 
        case('ka') 
            % Al2Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 2 Al+3 
            therm_ref = 10d0^(7.435d0);
            ha = -35.3d0*cal2j;
            tc_ref = 25d0;
            % from PHREEQC.DAT 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        % case('ab')
            % NaAlSi3O8 + 4 H+ = Na+ + Al3+ + 3SiO2 + 2H2O
            % therm_ref = 10d0^3.412182823d0
            % ha = -54.15042876d0
            % tc_ref = 15d0
            % from Kanzaki and Murakami 2018
            % therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('kfs')
            % K-feldspar  + 4 H+  = 2 H2O  + K+  + Al+++  + 3 SiO2(aq)
            therm_ref = 10d0^0.227294204d0;
            ha = -26.30862098d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('anl')
            % NaAlSi2O6*H2O  + 5 H2O  = Na+  + Al(OH)4-  + 2 Si(OH)4(aq)
            therm_ref = 10d0^(-16.06d0);
            ha = 101d0;
            tc_ref = 25d0;
            % from Wilkin and Barnes 1998
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('nph')
            % Nepheline  + 4 H+  = 2 H2O  + SiO2(aq)  + Al+++  + Na+
            therm_ref = 10d0^(14.93646757d0);
            ha = -130.8197467d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('fo')
            % Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O
            therm_ref = 10d0^29.41364324d0;
            ha = -208.5932252d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('fa')
            % Fa + 4H+ = 2Fe2+ + SiO2(aq) + 2H2O
            therm_ref = 10d0^19.98781342d0;
            ha = -153.7676621d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        % case('an')
            % CaAl2Si2O8 + 8H+ = Ca2+ + 2 Al3+ + 2SiO2 + 4H2O
            % therm_ref = 10d0^28.8615308d0;
            % ha = -292.8769275d0;
            % tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            % therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('cc')
            % CaCO3 = Ca2+ + CO32-
            therm_ref = 10d0^(-8.43d0);
            ha = -8.028943471d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2015
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('arg')
            % CaCO3 = Ca2+ + CO32-
            therm_ref = 10d0^(-8.3d0);
            ha = -12d0;
            tc_ref = 25d0;
            % from minteq.v4
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('dlm') % disordered
            % CaMg(CO3)2 = Ca+2 + Mg+2 + 2CO3-2
            therm_ref = 10d0^(-16.54d0);
            ha = -46.4d0;
            tc_ref = 25d0;
            % from minteq.v4
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('gb')
            % Al(OH)3 + 3 H+ = Al+3 + 3 H2O
            therm_ref = 10d0^(8.11d0);
            ha = -22.80d0*cal2j;
            tc_ref = 25d0;
            % from PHREEQC.DAT 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('amsi')
            % SiO2 + 2 H2O = H4SiO4
            therm_ref = 10d0^(-2.71d0);
            ha = 3.340d0*cal2j;
            tc_ref = 25d0;
            % from PHREEQC.DAT 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('qtz')
            % SiO2 + 2H2O = H4SiO4
            therm_ref = 10d0^(-4d0);
            ha = 22.36d0;
            tc_ref = 25d0;
            % from minteq.v4 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('gt')
            % Fe(OH)3 + 3 H+ = Fe+3 + 2 H2O
            therm_ref = 10d0^(0.5345d0);
            ha = -61.53703d0;
            tc_ref = 25d0;
            % from Sugimori et al. 2012 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('hm')
            % Fe2O3 + 6H+ = 2Fe+3 + 3H2O
            therm_ref = 10d0^(-1.418d0);
            ha = -128.987d0;
            tc_ref = 25d0;
            % from minteq.v4
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('ct')
            % Mg3Si2O5(OH)4 + 6 H+ = H2O + 2 H4SiO4 + 3 Mg+2
            therm_ref = 10d0^(32.2d0);
            ha = -46.800d0*cal2j;
            tc_ref = 25d0;
            % from PHREEQC.DAT 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('mscv')
            % KAl2(AlSi3O10)(OH)2 + 10 H+  = 6 H2O  + 3 SiO2(aq)  + K+  + 3 Al+++
            therm_ref = 10d0^(15.97690572d0);
            ha = -230.7845245d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('plgp')
            % KMg3(AlSi3O10)(OH)2 + 10 H+  = 6 H2O  + 3 SiO2(aq)  + Al+++  + K+  + 3 Mg++
            therm_ref = 10d0^(40.12256823d0);
            ha = -312.7817497d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018 
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('cabd')
            % Beidellit-Ca  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Ca++
            therm_ref = 10d0^(7.269946518d0);
            ha = -157.0186168d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('mgbd')
            % Beidellit-Mg  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .165 Mg++
            therm_ref = 10d0^(7.270517113d0);
            ha = -160.1864268d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('nabd')
            % Beidellit-Na  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .33 Na+
            therm_ref = 10d0^(7.288837383d0);
            ha = -150.7328834d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('kbd')
            % Beidellit-K  + 7.32 H+  = 4.66 H2O  + 2.33 Al+++  + 3.67 SiO2(aq)  + .33 K+
            therm_ref = 10d0^(6.928086412d0);
            ha = -145.6776905d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('ill')
            % Illite  + 8 H+  = 5 H2O  + .6 K+  + .25 Mg++  + 2.3 Al+++  + 3.5 SiO2(aq)
            therm_ref = 10d0^(10.8063184d0);
            ha = -166.39733d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        % case('dp')
            % Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
            % therm_ref = 10d0^(21.79853309d0);
            % ha = -138.6020832d0;
            % tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            % therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        % case('hb')
            % Diopside  + 4 H+  = Ca++  + 2 H2O  + Mg++  + 2 SiO2(aq)
            % therm_ref = 10d0^(20.20981116d0);
            % ha = -128.5d0;
            % tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            % therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('tm')
            % Tremolite  + 14 H+  = 8 H2O  + 8 SiO2(aq)  + 2 Ca++  + 5 Mg++
            therm_ref = 10d0^(61.6715d0);
            ha = -429.0d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('antp')
            % Anthophyllite (Mg2Mg5(Si8O22)(OH)2) + 14 H+  = 8 H2O  + 7 Mg++  + 8 SiO2(aq)
            therm_ref = 10d0^(70.83527792d0);
            ha = -508.6621624d0;
            tc_ref = 15d0;
            % from Kanzaki and Murakami 2018
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case('gps')
            % CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O
            therm_ref = 10d0^(-4.61d0);
            ha = 1d0;
            tc_ref = 25d0;
            % from minteq.v4
            therm = k_arrhenius(therm_ref,tc_ref+tempk_0,tc+tempk_0,ha,rg);
        case{'la','ab','an','by','olg','and'}
            % CaxNa(1-x)Al(1+x)Si(3-x)O8 + (4x + 4) = xCa+2 + (1-x)Na+ + (1+x)Al+++ + (3-x)SiO2(aq) 
            % obtaining Anorthite 
            therm_ref_1 = 10d0^28.8615308d0;
            ha_1 = -292.8769275d0;
            tc_ref_1 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg); % rg in kJ mol^-1 K^-1
            delG_1 = - rg*(tc+tempk_0)*log(therm_1); % del-G = -RT ln K  now in kJ mol-1
            % Then albite 
            therm_ref_2 = 10d0^3.412182823d0;
            ha_2 = -54.15042876d0;
            tc_ref_2 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg);
            delG_2 = - rg*(tc+tempk_0)*log(therm_2); % del-G = -RT ln K  now in kJ mol-1
            
            if (ss_x == 1d0)  
                delG = delG_1; % ideal anorthite
            elseif (ss_x == 0d0)  
                delG = delG_2; % ideal albite
            elseif (ss_x > 0d0 && ss_x < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x));
            end 
            therm = exp(-delG/(rg*(tc+tempk_0)));
        case{'cpx','hb','dp'}
            % FexMg(1-x)CaSi2O6 + 4 H+  = Ca++  + 2 H2O  + xFe++ + (1-x)Mg++  + 2 SiO2(aq)
            % obtaining hedenbergite 
            therm_ref_1 = 10d0^(20.20981116d0);
            ha_1 = -128.5d0;
            tc_ref_1 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg); % rg in kJ mol^-1 K^-1
            delG_1 = - rg*(tc+tempk_0)*log(therm_1); % del-G = -RT ln K  now in kJ mol-1
            % Then diopside 
            therm_ref_2 = 10d0^(21.79853309d0);
            ha_2 = -138.6020832d0;
            tc_ref_2 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg);
            delG_2 = - rg*(tc+tempk_0)*log(therm_2); % del-G = -RT ln K  now in kJ mol-1
            
            if (ss_x == 1d0)  
                delG = delG_1; % ideal hedenbergite
            elseif (ss_x == 0d0)  
                delG = delG_2; % ideal diopside
            elseif (ss_x > 0d0 && ss_x < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x));
            end 
            therm = exp(-delG/(rg*(tc+tempk_0)));
        case{'opx','en','fer'}
            % FexMg(1-x)SiO3 + 2 H+  = xFe++ + (1-x)Mg++  +  SiO2(aq)
            % obtaining ferrosilite
            therm_ref_1 = 10d0^(7.777162795d0);
            ha_1 = -60.08612326d0;
            tc_ref_1 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg); % rg in kJ mol^-1 K^-1
            delG_1 = - rg*(tc+tempk_0)*log(therm_1); % del-G = -RT ln K  now in kJ mol-1
            % Then enstatite 
            therm_ref_2 = 10d0^(11.99060855d0);
            ha_2 = -85.8218778d0;
            tc_ref_2 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg);
            delG_2 = - rg*(tc+tempk_0)*log(therm_2); % del-G = -RT ln K  now in kJ mol-1
            
            if (ss_x == 1d0)  
                delG = delG_1; % ideal ferrosilite
            elseif (ss_x == 0d0)  
                delG = delG_2; % ideal enstatite
            elseif (ss_x > 0d0 && ss_x < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x));
            end 
            therm = exp(-delG/(rg*(tc+tempk_0)));
        case('agt')
            % obtaining opx (FexMg(1-x)SiO3 + 2 H+  = xFe++ + (1-x)Mg++  +  SiO2(aq))
            % obtaining ferrosilite
            therm_ref_1 = 10d0^(7.777162795d0);
            ha_1 = -60.08612326d0;
            tc_ref_1 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_1 = k_arrhenius(therm_ref_1,tc_ref_1+tempk_0,tc+tempk_0,ha_1,rg); % rg in kJ mol^-1 K^-1
            % converting to the formula Fe2Si2O6 + 4 H+  = 2Fe++ + 2SiO2(aq)
            therm_1 = therm_1^2d0;
            delG_1 = - rg*(tc+tempk_0)*log(therm_1); % del-G = -RT ln K  now in kJ mol-1
            
            % Then enstatite 
            therm_ref_2 = 10d0^(11.99060855d0);
            ha_2 = -85.8218778d0;
            tc_ref_2 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_2 = k_arrhenius(therm_ref_2,tc_ref_2+tempk_0,tc+tempk_0,ha_2,rg);
            % converting to the formula Mg2Si2O6 + 4 H+  = 2Mg++ + 2SiO2(aq)
            therm_2 = therm_2^2d0;
            delG_2 = - rg*(tc+tempk_0)*log(therm_2); % del-G = -RT ln K  now in kJ mol-1
            
            if (ss_x == 1d0)  
                delG_3 = delG_1; % ideal ferrosilite
            elseif (ss_x == 0d0)  
                delG_3 = delG_2; % ideal enstatite
            elseif (ss_x > 0d0 && ss_x < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG_3 = ss_x*delG_1 + (1d0-ss_x)*delG_2 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x));
            end 
            therm_3 = exp(-delG_3/(rg*(tc+tempk_0)));
            
            % obtaining cpx (FexMg(1-x)CaSi2O6 + 4 H+  = Ca++  + 2 H2O  + xFe++ + (1-x)Mg++  + 2 SiO2(aq))
            % obtaining hedenbergite 
            therm_ref_4 = 10d0^(20.20981116d0);
            ha_4 = -128.5d0;
            tc_ref_4 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_4 = k_arrhenius(therm_ref_4,tc_ref_4+tempk_0,tc+tempk_0,ha_4,rg); % rg in kJ mol^-1 K^-1
            delG_4 = - rg*(tc+tempk_0)*log(therm_4); % del-G = -RT ln K  now in kJ mol-1
            % Then diopside 
            therm_ref_5 = 10d0^(21.79853309d0);
            ha_5 = -138.6020832d0;
            tc_ref_5 = 15d0;
            % from Kanzaki and Murakami 2018
            therm_5 = k_arrhenius(therm_ref_5,tc_ref_5+tempk_0,tc+tempk_0,ha_5,rg);
            delG_5 = - rg*(tc+tempk_0)*log(therm_5); % del-G = -RT ln K  now in kJ mol-1
            
            if (ss_x == 1d0)  
                delG_6 = delG_4; % ideal hedenbergite
            elseif (ss_x == 0d0)  
                delG_6 = delG_5; % ideal diopside
            elseif (ss_x > 0d0 && ss_x < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG_6 = ss_x*delG_4 + (1d0-ss_x)*delG_5 + rg*(tc+tempk_0)*(ss_x*log(ss_x)+(1d0-ss_x)*log(1d0-ss_x));
            end 
            therm_6 = exp(-delG_6/(rg*(tc+tempk_0)));
            
            % finally mixing opx and cpx 
            if (ss_y == 1d0)  
                delG = delG_3; % ideal opx
            elseif (ss_y == 0d0)  
                delG = delG_6; % ideal cpx
            elseif (ss_y > 0d0 && ss_y < 1d0)   % solid solution 
                % ideal(?) mixing (after Gislason and Arnorsson, 1993)
                delG = ss_y*delG_3 + (1d0-ss_y)*delG_6 + rg*(tc+tempk_0)*(ss_y*log(ss_y)+(1d0-ss_y)*log(1d0-ss_y));
            end 
            therm = exp(-delG/(rg*(tc+tempk_0)));
            
        case('g1')
            therm = 0.121d0; % mo2 Michaelis, Davidson et al. (2012)
        case('g2')
            therm = 0.121d0; % mo2 Michaelis, Davidson et al. (2012)
            % therm = 0.121d-1 % mo2 Michaelis, Davidson et al. (2012) x 0.1
            % therm = 0.121d-2 % mo2 Michaelis, Davidson et al. (2012) x 0.01
            % therm = 0.121d-3 % mo2 Michaelis, Davidson et al. (2012) x 0.001
            % therm = 0.121d-6 % mo2 Michaelis, Davidson et al. (2012) x 1e-6
        case('g3')
            therm = 0.121d0; % mo2 Michaelis, Davidson et al. (2012)
        otherwise
            therm = 0d0;
    end 

end


function [ ...
    ucv,kw,daq_all,dgasa_all,dgasg_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3,keqaq_nh3 ...% output
    ,ksld_all,keqsld_all,krxn1_ext_all,krxn2_ext_all ...% output
    ] = coefs_v2( ...
    nz,rg,rg2,tc,sec2yr,tempk_0,pro ...% input
    ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nrxn_ext_all ...% input
    ,chraq_all,chrgas_all,chrsld_all,chrrxn_ext_all ...% input
    ,nsp_gas,nsp_gas_cnst,chrgas,chrgas_cnst,mgas,mgasc,mgasth_all,mv_all,staq_all ...%input
    ) 

    cal2j = 4.184d0; 
    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);
    [ieqaq_co3,ieqaq_hco3]=deal(1,2);
    [ieqaq_so4,ieqaq_so42]=deal(1,2);
    [ieqaq_no3,ieqaq_no32]=deal(1,2);
    [ieqaq_nh3,ieqaq_nh32]=deal(1,2);
    thon = -1d100;

    % initialize output variables
    ucv = 0;kw = 0;
    daq_all = zeros(nsp_aq_all,1,'double');
    dgasa_all = zeros(nsp_gas_all,1,'double');
    dgasg_all = zeros(nsp_gas_all,1,'double');
    keqgas_h = zeros(nsp_gas_all,3,'double');
    keqaq_h = zeros(nsp_aq_all,4,'double');
    keqaq_c = zeros(nsp_aq_all,2,'double');
    keqaq_s = zeros(nsp_aq_all,2,'double');
    keqaq_no3 = zeros(nsp_aq_all,2,'double');
    keqaq_nh3 = zeros(nsp_aq_all,2,'double');
    ksld_all = zeros(nsp_sld_all,nz,'double');
    keqsld_all = zeros(nsp_sld_all,1,'double');
    krxn1_ext_all = zeros(nrxn_ext_all,nz,'double');
    krxn2_ext_all = zeros(nrxn_ext_all,nz,'double');

    % loca variables (int/float)
    oh=zeros(nz,1,'double');po2=zeros(nz,1,'double');kin=zeros(nz,1,'double');dkin_dmsp=zeros(nz,1,'double');
    kho=0;po2th=0;mv_tmp=0;therm=0;ss_x=0;ss_y=0;
    mgas_loc = zeros(nsp_gas_all,nz,'double');

    % --- start getting coeffs

    ucv = 1.0d0/(rg2*(tempk_0+tc));

    % Aq species diffusion from Li and Gregory 1974 except for Si which is based on Rebreanu et al. 2008
    daq_all(find(chraq_all=='fe2'))= k_arrhenius(1.7016d-2    , 15d0+tempk_0, tc+tempk_0, 19.615251d0, rg);
    daq_all(find(chraq_all=='fe3'))= k_arrhenius(1.5664d-2    , 15d0+tempk_0, tc+tempk_0, 14.33659d0 , rg);
    daq_all(find(chraq_all=='so4'))= k_arrhenius(2.54d-2      , 15d0+tempk_0, tc+tempk_0, 20.67364d0 , rg);
    daq_all(find(chraq_all=='no3'))= k_arrhenius(4.6770059d-2 , 15d0+tempk_0, tc+tempk_0, 18.00685d0 , rg);
    daq_all(find(chraq_all=='na')) = k_arrhenius(3.19d-2      , 15d0+tempk_0, tc+tempk_0, 20.58566d0 , rg);
    daq_all(find(chraq_all=='k'))  = k_arrhenius(4.8022699d-2 , 15d0+tempk_0, tc+tempk_0, 18.71816d0 , rg);
    daq_all(find(chraq_all=='mg')) = k_arrhenius(1.7218079d-2 , 15d0+tempk_0, tc+tempk_0, 18.51979d0 , rg);
    daq_all(find(chraq_all=='si')) = k_arrhenius(2.682396d-2  , 15d0+tempk_0, tc+tempk_0, 22.71378d0 , rg);
    daq_all(find(chraq_all=='ca')) = k_arrhenius(1.9023312d-2 , 15d0+tempk_0, tc+tempk_0, 20.219661d0, rg);
    daq_all(find(chraq_all=='al')) = k_arrhenius(1.1656226d-2 , 15d0+tempk_0, tc+tempk_0, 21.27788d0 , rg);

    % values used in Kanzaki and Murakami 2016 for oxygen 
    dgasa_all(find(chrgas_all=='po2')) = k_arrhenius(5.49d-2 , 15d0+tempk_0, tc+tempk_0, 20.07d0 , rg);
    dgasg_all(find(chrgas_all=='po2')) = k_arrhenius(6.09d2  , 15d0+tempk_0, tc+tempk_0, 4.18d0  , rg);

    % assuming a value of 0.14 cm2/sec (e.g., Pritchard and Currie, 1982) and O2 gas activation energy for CO2 gas 
    % and CO32- diffusion from Li and Greogy 1974 for aq CO2 
    dgasa_all(find(chrgas_all=='pco2')) = k_arrhenius(2.2459852d-2, 15d0+tempk_0, tc+tempk_0, 21.00564d0, rg);
    dgasg_all(find(chrgas_all=='pco2')) = k_arrhenius(441.504d0   , 15d0+tempk_0, tc+tempk_0, 4.18d0    , rg);

    % NH4+ diffusion for aqueous diffusion from Schulz and Zabel 2005
    % NH3 diffusion in air from Massman 1998
    dgasa_all(find(chrgas_all=='pnh3')) = k_arrhenius(4.64d-02    , 15d0+tempk_0, tc+tempk_0, 19.15308d0, rg);
    dgasg_all(find(chrgas_all=='pnh3')) = 0.1978d0*((tc+tempk_0)/(0d0+tempk_0))^1.81d0 * sec2yr *1d-4; % sec2yr*1d-4 converting cm2 to m2 and sec-1 to yr-1

    % assuming the same diffusion as CO2 diffusion (e.g., Pritchard and Currie, 1982) for gaseous N2O 
    % N2O(aq) diffusion from Schulz and Zabel 2005
    dgasa_all(find(chrgas_all=='pn2o')) = k_arrhenius(4.89d-02    , 15d0+tempk_0, tc+tempk_0, 20.33417d0, rg);
    dgasg_all(find(chrgas_all=='pn2o')) = k_arrhenius(441.504d0   , 15d0+tempk_0, tc+tempk_0, 4.18d0    , rg);

    kw = -14.93d0+0.04188d0*tc-0.0001974d0*tc^2d0+0.000000555d0*tc^3d0-0.0000000007581d0*tc^4d0 ; % Murakami et al. 2011
    kw = k_arrhenius(10d0^(-14.35d0), tempk_0+15.0d0, tempk_0+tc, 58.736742d0, rg); % from Kanzaki and Murakami 2015

    oh(:) = kw./pro(:);


    % kho = k_arrhenius(10.0d0^(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg)
    keqgas_h(find(chrgas_all=='po2'),ieqgas_h0) = ...
        k_arrhenius(10d0^(-2.89d0), tempk_0+25.0d0, tempk_0+tc, -13.2d0, rg);
    kho = keqgas_h(find(chrgas_all=='po2'),ieqgas_h0);

    keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0) = ...
        k_arrhenius(10d0^(-1.34d0), tempk_0+15.0d0, tempk_0+tc, -21.33183d0, rg); % from Kanzaki and Murakami 2015
    keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1) = ...
        k_arrhenius(10d0^(-6.42d0), tempk_0+15.0d0, tempk_0+tc, 11.94453d0, rg); % from Kanzaki and Murakami 2015
    keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2) = ...
        k_arrhenius(10d0^(-10.43d0), tempk_0+15.0d0, tempk_0+tc, 17.00089d0, rg); % from Kanzaki and Murakami 2015

    keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h0) = ...
        k_arrhenius(10d0^(1.770d0), tempk_0+25.0d0, tempk_0+tc, -8.170d0*cal2j, rg); % from WATEQ4F.DAT 
    keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h1) = ...
        k_arrhenius(10d0^(-9.252d0), tempk_0+25.0d0, tempk_0+tc, 12.48d0*cal2j, rg); % from WATEQ4F.DAT (NH4+ = NH3 + H+)

    keqgas_h(find(chrgas_all=='pn2o'),ieqgas_h0) = ...
        k_arrhenius(0.033928709d0, tempk_0+15.0d0, tempk_0+tc, -22.21661d0, rg); % % N2O solubility from Weiss and Price 1980 MC assuming 0 salinity


    % SO4-2 + H+ = HSO4- 
    % keqaq_s(find(chraq_all=='so4'),ieqaq_so4) = ...
        % k_arrhenius(10d0^(1.988d0),25d0+tempk_0,tc+tempk_0,3.85d0*cal2j,rg); % from PHREEQC.DAT
    keqaq_h(find(chraq_all=='so4'),ieqaq_h1) = ...
        k_arrhenius(10d0^(1.988d0),25d0+tempk_0,tc+tempk_0,3.85d0*cal2j,rg); % from PHREEQC.DAT
    % SO4-2 + NH4+ = NH4SO4-
    % keqaq_nh3(find(chraq_all=='so4'),ieqaq_nh3) = ...
        % k_arrhenius(10d0^(1.03d0),25d0+tempk_0,tc+tempk_0,0d0,rg); % from MINTEQV4.DAT 

    % H+ + NO3- = HNO3 
    % keqaq_no3(find(chraq_all=='no3'),ieqaq_no3) = 1d0/35.5d0 ;% from Levanov et al. 2017 
    % keqaq_no3(find(chraq_all=='no3'),ieqaq_no3) = 1d0/(10d0^1.3d0) ;% from Maggi et al. 2007 
    keqaq_h(find(chraq_all=='no3'),ieqaq_h1) = 1d0/35.5d0 ;% from Levanov et al. 2017 
    keqaq_h(find(chraq_all=='no3'),ieqaq_h1) = 1d0/(10d0^1.3d0) ;% from Maggi et al. 2007 
    % (temperature dependence is assumed to be 0) 

    % Al3+ + H2O = Al(OH)2+ + H+
    keqaq_h(find(chraq_all=='al'),ieqaq_h1) = ...
        k_arrhenius(10d0^(-5d0),25d0+tempk_0,tc+tempk_0,11.49d0*cal2j,rg); % from PHREEQC.DAT 
    % Al3+ + 2H2O = Al(OH)2+ + 2H+
    keqaq_h(find(chraq_all=='al'),ieqaq_h2) = ...
        k_arrhenius(10d0^(-10.1d0),25d0+tempk_0,tc+tempk_0,26.90d0*cal2j,rg) ;% from PHREEQC.DAT 
    % Al3+ + 3H2O = Al(OH)3 + 3H+
    keqaq_h(find(chraq_all=='al'),ieqaq_h3) = ...
        k_arrhenius(10d0^(-16.9d0),25d0+tempk_0,tc+tempk_0,39.89d0*cal2j,rg); % from PHREEQC.DAT 
    % Al3+ + 4H2O = Al(OH)4- + 4H+
    keqaq_h(find(chraq_all=='al'),ieqaq_h4) = ...
        k_arrhenius(10d0^(-22.7d0),25d0+tempk_0,tc+tempk_0,42.30d0*cal2j,rg); % from PHREEQC.DAT 
    % Al+3 + SO4-2 = AlSO4+
    keqaq_s(find(chraq_all=='al'),ieqaq_so4) = ...
        k_arrhenius(10d0^(3.5d0),25d0+tempk_0,tc+tempk_0,2.29d0*cal2j,rg); % from PHREEQC.DAT 
    % Al+3 + 2SO4-2 = Al(SO4)2-
    % ignoring for now
    % keqaq_s(find(chraq_all=='al'),ieqaq_so42) = ...
        % k_arrhenius(10d0^(5.0d0),25d0+tempk_0,tc+tempk_0,3.11d0*cal2j,rg) % from PHREEQC.DAT 

    % H4SiO4 = H3SiO4- + H+
    keqaq_h(find(chraq_all=='si'),ieqaq_h1) = ...
        k_arrhenius(10d0^(-9.83d0),25d0+tempk_0,tc+tempk_0,6.12d0*cal2j,rg); % from PHREEQC.DAT 
    % H4SiO4 = H2SiO4-2 + 2 H+
    keqaq_h(find(chraq_all=='si'),ieqaq_h2) = ...
        k_arrhenius(10d0^(-23d0),25d0+tempk_0,tc+tempk_0,17.6d0*cal2j,rg) ;% from PHREEQC.DAT 


    % Mg2+ + H2O = Mg(OH)+ + H+
    keqaq_h(find(chraq_all=='mg'),ieqaq_h1) = ...
        k_arrhenius(10d0^(-11.44d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg) ;% from PHREEQC.DAT 
    % Mg2+ + CO32- = MgCO3 
    keqaq_c(find(chraq_all=='mg'),ieqaq_co3) = ...
        k_arrhenius(10d0^(2.98d0),25d0+tempk_0,tc+tempk_0,2.713d0*cal2j,rg); % from PHREEQC.DAT 
    % Mg2+ + H+ + CO32- = MgHCO3
    keqaq_c(find(chraq_all=='mg'),ieqaq_hco3) = ... 
        k_arrhenius(10d0^(11.399d0),25d0+tempk_0,tc+tempk_0,-2.771d0*cal2j,rg); % from PHREEQC.DAT 
    % Mg+2 + SO4-2 = MgSO4
    keqaq_s(find(chraq_all=='mg'),ieqaq_so4) = ... 
        k_arrhenius(10d0^(2.37d0),25d0+tempk_0,tc+tempk_0, 4.550d0*cal2j,rg); % from PHREEQC.DAT 

    % Ca2+ + H2O = Ca(OH)+ + H+
    keqaq_h(find(chraq_all=='ca'),ieqaq_h1) =  ...
        k_arrhenius(10d0^(-12.78d0),25d0+tempk_0,tc+tempk_0,15.952d0*cal2j,rg); % from PHREEQC.DAT 
    % (No delta_h is reported so used the same value for Mg)
    % Ca2+ + CO32- = CaCO3 
    keqaq_c(find(chraq_all=='ca'),ieqaq_co3) = ...
        k_arrhenius(10d0^(3.224d0),25d0+tempk_0,tc+tempk_0,3.545d0*cal2j,rg); % from PHREEQC.DAT 
    % Ca2+ + H+ + CO32- = CaHCO3
    keqaq_c(find(chraq_all=='ca'),ieqaq_hco3) = ...
        k_arrhenius(10d0^(11.435d0),25d0+tempk_0,tc+tempk_0,-0.871d0*cal2j,rg); % from PHREEQC.DAT 
    % Ca+2 + SO4-2 = CaSO4
    keqaq_s(find(chraq_all=='ca'),ieqaq_so4) = ...
        k_arrhenius(10d0^(2.25d0),25d0+tempk_0,tc+tempk_0,1.325d0*cal2j,rg); % from PHREEQC.DAT 
    % Ca+2 + NO3- = CaNO3+
    % keqaq_no3(find(chraq_all=='ca'),ieqaq_no3) = ...
        % k_arrhenius(10d0^(0.5d0),25d0+tempk_0,tc+tempk_0,-5.4d0,rg); % from MINTEQV4.DAT           
    % Ca+2 + NH4+ = CaNH3+2 + H+
    % keqaq_nh3(find(chraq_all=='ca'),ieqaq_nh3) = ...
        % k_arrhenius(10d0^(-9.144d0),25d0+tempk_0,tc+tempk_0,0d0,rg); % from MINTEQV4.DAT 
    % Ca+2 + 2NH4+ = Ca(NH3)2+2 + 2H+
    % ignoring for now
    % keqaq_nh3(find(chraq_all=='ca'),ieqaq_nh32) = ...
        % k_arrhenius(10d0^(-18.788d0),25d0+tempk_0,tc+tempk_0,0d0,rg); % from MINTEQV4.DAT 

        
    % Fe2+ + H2O = Fe(OH)+ + H+
    keqaq_h(find(chraq_all=='fe2'),ieqaq_h1) = ...
        k_arrhenius(10d0^(-9.51d0),25d0+tempk_0,tc+tempk_0, 40.3d0,rg); % from Kanzaki and Murakami 2016
    % Fe2+ + CO32- = FeCO3 
    keqaq_c(find(chraq_all=='fe2'),ieqaq_co3) = ...
        k_arrhenius(10d0^(5.69d0),25d0+tempk_0,tc+tempk_0, -45.6d0,rg); % from Kanzaki and Murakami 2016
    % Fe2+ + H+ + CO32- = FeHCO3
    keqaq_c(find(chraq_all=='fe2'),ieqaq_hco3) = ...
        k_arrhenius(10d0^(1.47d0),25d0+tempk_0,tc+tempk_0, -18d0,rg) ...% from Kanzaki and Murakami 2016 
        /keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2); 
    % Fe+2 + SO4-2 = FeSO4
    keqaq_s(find(chraq_all=='fe2'),ieqaq_so4) = ...
        k_arrhenius(10d0^(2.25d0),25d0+tempk_0,tc+tempk_0,3.230d0*cal2j,rg); % from PHREEQC.DAT 


    % Fe3+ + H2O = Fe(OH)2+ + H+
    keqaq_h(find(chraq_all=='fe3'),ieqaq_h1) = ...
        k_arrhenius(10d0^(-2.19d0),25d0+tempk_0,tc+tempk_0,10.4d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe3+ + 2H2O = Fe(OH)2+ + 2H+
    keqaq_h(find(chraq_all=='fe3'),ieqaq_h2) = ...
        k_arrhenius(10d0^(-5.67d0),25d0+tempk_0,tc+tempk_0,17.1d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe3+ + 3H2O = Fe(OH)3 + 3H+
    keqaq_h(find(chraq_all=='fe3'),ieqaq_h3) = ...
        k_arrhenius(10d0^(-12.56d0),25d0+tempk_0,tc+tempk_0,24.8d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe3+ + 4H2O = Fe(OH)4- + 4H+
    keqaq_h(find(chraq_all=='fe3'),ieqaq_h4) = ...
        k_arrhenius(10d0^(-21.6d0),25d0+tempk_0,tc+tempk_0,31.9d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe+3 + SO4-2 = FeSO4+
    keqaq_s(find(chraq_all=='fe3'),ieqaq_so4) = ...
        k_arrhenius(10d0^(4.04d0),25d0+tempk_0,tc+tempk_0,3.91d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe+3 + 2 SO4-2 = Fe(SO4)2-
    % ignoring for now
    % keqaq_s(find(chraq_all=='fe3'),ieqaq_so42) = ...
        % k_arrhenius(10d0^(5.38d0),25d0+tempk_0,tc+tempk_0,4.60d0*cal2j,rg); % from PHREEQC.DAT 
    % Fe+3 + NO3- = FeNO3+2
    % keqaq_no3(find(chraq_all=='fe3'),ieqaq_no3) = ...
        % k_arrhenius(10d0^(1d0),25d0+tempk_0,tc+tempk_0,-37d0,rg); % from MINTEQV4.DAT 



    % Na+ + CO3-2 = NaCO3-
    keqaq_c(find(chraq_all=='na'),ieqaq_co3) = ... 
        k_arrhenius(10d0^(1.27d0),25d0+tempk_0,tc+tempk_0, 8.91d0*cal2j,rg); % from PHREEQC.DAT 
    % Na+ + H + CO3- = NaHCO3
    keqaq_c(find(chraq_all=='na'),ieqaq_hco3) = ... 
        k_arrhenius(10d0^(-0.25d0),25d0+tempk_0,tc+tempk_0, -1d0*cal2j,rg) ...% from PHREEQC.DAT for Na+ + HCO3- = NaHCO3
        /keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2);  % HCO3- = CO32- + H+
    % Na+ + SO4-2 = NaSO4-
    keqaq_s(find(chraq_all=='na'),ieqaq_so4) = ... 
        k_arrhenius(10d0^(0.7d0),25d0+tempk_0,tc+tempk_0, 1.120d0*cal2j,rg); % from PHREEQC.DAT 



    % K+ + SO4-2 = KSO4-
    keqaq_s(find(chraq_all=='k'),ieqaq_so4) = ... 
        k_arrhenius(10d0^(0.85d0),25d0+tempk_0,tc+tempk_0, 2.250d0*cal2j,rg); % from PHREEQC.DAT 


    % keqaq_s = 0d0


    %%% ----------- Solid phases ------------------------%%

    [mgas_loc] = get_mgasx_all( ...
        nz,nsp_gas_all,nsp_gas,nsp_gas_cnst ...
        ,chrgas,chrgas_all,chrgas_cnst ...
        ,mgas,mgasc ...
        );

    for isps = 1: nsp_sld_all
        mv_tmp = mv_all(isps);
        mineral = chrsld_all(isps);
        
        [kin,dkin_dmsp] = sld_kin( ...
            nz,rg,tc,sec2yr,tempk_0,pro,kw,kho,mv_tmp ...% input
            ,nsp_gas_all,chrgas_all,mgas_loc ...% input
            ,mineral,'xxxxx' ...% input 
            ); 
        ksld_all(isps,:) = kin(:);
        
        % check for solid solution 
        switch (mineral) 
            case{'la','ab','an','by','olg','and'}
                ss_x = staq_all(isps, find(chraq_all=='ca'));
                ss_y = 0d0; % non-zero if it is a solid solution 
            case{'cpx','hb','dp'} 
                ss_x = staq_all(isps, find(chraq_all=='fe2'));
                ss_y = 0d0; % non-zero if it is a solid solution 
            case{'opx','en','fer'} 
                ss_x = staq_all(isps, find(chraq_all=='fe2'));
                ss_y = 0d0; % non-zero if it is a solid solution 
            case('agt') 
                ss_y = 1d0 - staq_all(isps, find(chraq_all=='ca'));
                ss_x = staq_all(isps, find(chraq_all=='fe2'))/(1d0+ ss_y ); % non-zero if it is a solid solution 
            otherwise 
                ss_x = 0d0; % non-zero if it is a solid solution 
                ss_y = 0d0; % non-zero if it is a solid solution 
        end 
        
        [therm] = sld_therm( ...
            rg,tc,tempk_0,ss_x,ss_y ...% input
            ,mineral ...% input
            );
        
        % correction of thermodynamic data wrt primary species
        switch (mineral) 
            case('anl')
                % replacing Al(OH)42- with Al+++ as primary Al species
                therm = therm/keqaq_h(find(chraq_all=='al'),ieqaq_h4);
            otherwise 
                % do nothing
        end 
        
        keqsld_all(isps) = therm;
    end


    %--------- other reactions -------------% 

    krxn1_ext_all(find(chrrxn_ext_all=='fe2o2'),:) = ...
        1d0;
        % max(8.0d13*60.0d0*24.0d0*365.0d0*(kw/pro)^2.0d0, 1d-7*60.0d0*24.0d0*365.0d0) ...   
        % mol L^-1 yr^-1 (25 deg C), Singer and Stumm (1970)excluding the term (c*po2)
        % *merge(0d0,1d0,po2<po2th*thon)
         
    krxn1_ext_all(find(chrrxn_ext_all=='resp'),:) = 0.71d0; % vmax mol m^-3, yr^-1, max soil respiration, Wood et al. (1993)
    krxn1_ext_all(find(chrrxn_ext_all=='resp'),:) = ...
        krxn1_ext_all(find(chrrxn_ext_all=='resp'),:); %*1d1 % reducing a bit to be fitted with modern soil pco2

    krxn2_ext_all(find(chrrxn_ext_all=='resp'),:) = 0.121d0; % mo2 Michaelis, Davidson et al. (2012)
         
         
    krxn1_ext_all(find(chrrxn_ext_all=='omomb'),:) = 0.01d0*24d0*365d0; % mg C mg-1 MBC yr-1
    % converted from 0.01 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

    krxn2_ext_all(find(chrrxn_ext_all=='omomb'),:) = 250d0; % mg C g-1 soil  Georgiou et al. (2017)
         
         
    krxn1_ext_all(find(chrrxn_ext_all=='ombto'),:) = 0.00028d0*24d0*365d0; % mg C mg-1 MBC yr-1
    % converted from 0.00028 mg C mg-1 MBC hr-1 Georgiou et al. (2017)

    krxn2_ext_all(find(chrrxn_ext_all=='ombto'),:) = 2d0; % beta value Georgiou et al. (2017)


    krxn1_ext_all(find(chrrxn_ext_all=='pyfe3'),:) = ... 
        10.0d0^(-6.07d0)*60.0d0*60.0d0*24.0d0*365.0d0;  %% excluding the term (fe3^0.93/fe2^0.40)  

end


function [maqx_loc,mgasx_loc] = get_maqgasx_all( ...
    nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst ...
    ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst ...
    ,maqx,mgasx,maqc,mgasc ...
    )

    % initialize output variables
    maqx_loc = zeros(nsp_aq_all,nz,'double');
    mgasx_loc = zeros(nsp_gas_all,nz,'double');

    for ispa = 1: nsp_aq_all
        if (any(chraq==chraq_all(ispa)))  
            maqx_loc(ispa,:) =  maqx(find(chraq==chraq_all(ispa)),:);
        elseif (any(chraq_cnst==chraq_all(ispa)))  
            maqx_loc(ispa,:) =  maqc(find(chraq_cnst==chraq_all(ispa)),:);
        end 
    end 

    for ispg = 1: nsp_gas_all
        if (any(chrgas==chrgas_all(ispg)))  
            mgasx_loc(ispg,:) =  mgasx(find(chrgas==chrgas_all(ispg)),:);
        elseif (any(chrgas_cnst==chrgas_all(ispg)))  
            mgasx_loc(ispg,:) =  mgasc(find(chrgas_cnst==chrgas_all(ispg)),:);
        end 
    end 

end


function [base_charge] = get_base_charge( ...
    nsp_aq_all ... 
    ,chraq_all ... 
    )

    % initialize output variables
    base_charge = zeros(nsp_aq_all,1,'double');

    for ispa = 1: nsp_aq_all
        switch chraq_all(ispa)
            case('so4')
                base_charge(ispa) = -2d0;
            case('no3')
                base_charge(ispa) = -1d0;
            case('si')
                base_charge(ispa) = 0d0;
            case{'na','k'}
                base_charge(ispa) = 1d0;
            case{'fe2','mg','ca'}
                base_charge(ispa) = 2d0;
            case{'fe3','al'}
                base_charge(ispa) = 3d0;
            otherwise 
                error ('error in charge assignment')
        end 
    end
    
end


function [ ...
    dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
    ,maqf_loc  ...% output
    ] = get_maqf_all( ...
    nz,nsp_aq_all,nsp_gas_all ...
    ,chraq_all,chrgas_all ...
    ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
    ,mgasx_loc,maqx_loc,prox,so4f ...
    )
        
    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);

    % output initilize
    maqf_loc = zeros(nsp_aq_all,nz,'double');

    dmaqf_dpro = zeros(nsp_aq_all,nz,'double');
    dmaqf_dso4f =zeros(nsp_aq_all,nz,'double');
    dmaqf_dmaq = zeros(nsp_aq_all,nz,'double');
    dmaqf_dpco2 = zeros(nsp_aq_all,nz,'double');

    % local initilize
    kco2=0;k1=0;k2=0;k1no3=0;rspa_h=0;rspa_s=0;
    pco2x=zeros(nz,1,'double');

    % start 
    kco2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0);
    k1 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);
    k2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2);

    pco2x(:) = mgasx_loc(find(chrgas_all=='pco2'),:);

    k1no3 = keqaq_h(find(chraq_all=='no3'),ieqaq_h1);

    for ispa = 1: nsp_aq_all
        % annions
        if (chraq_all(ispa)=='no3' || chraq_all(ispa)=='so4')  
            switch(chraq_all(ispa))
                case('no3')
                    maqf_loc(ispa,:) = 1d0;
                    if (k1no3 > 0d0)  % currently NO3 complex with cations are ignored and thus no3f can be calculated analytically
                        maqf_loc(ispa,:) = maqf_loc(ispa,:) + k1no3*prox(:)';
                        dmaqf_dpro(ispa,:) = dmaqf_dpro(ispa,:) + k1no3;
                    end
                    dmaqf_dpro(ispa,:) = maqx_loc(ispa,:)*(-1d0)./maqf_loc(ispa,:).^2d0.*dmaqf_dpro(ispa,:);
                    dmaqf_dmaq(ispa,:) = 1d0./maqf_loc(ispa,:);
                    maqf_loc(ispa,:) = maqx_loc(ispa,:)./maqf_loc(ispa,:);
                case('so4') % currently SO4 complex with cations are included so that so4f is numerically calculated with pH (or charge balance)
                    maqf_loc(ispa,:) = so4f(:)';
                    dmaqf_dso4f(ispa,:) = 1d0;
            end
        % cations
        else 
            maqf_loc(ispa,:) = 1d0;
            % account for hydrolysis speces
            for ispa_h = 1:4
                rspa_h = double(ispa_h);
                if ( keqaq_h(ispa,ispa_h) > 0d0)  
                    maqf_loc(ispa,:) = maqf_loc(ispa,:) + keqaq_h(ispa,ispa_h)./prox(:)'.^rspa_h;
                    dmaqf_dpro(ispa,:) = dmaqf_dpro(ispa,:) + keqaq_h(ispa,ispa_h)*(-rspa_h)./prox(:)'.^(1d0+rspa_h);
                end 
            end 
            % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
            for ispa_c = 1:2
                if ( keqaq_c(ispa,ispa_c) > 0d0)  
                    if (ispa_c == 1)  % with CO3--
                        maqf_loc(ispa,:) = maqf_loc(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*pco2x(:)'./prox(:)'.^2d0;
                        dmaqf_dpro(ispa,:) = dmaqf_dpro(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*pco2x(:)'*(-2d0)./prox(:)'.^3d0;
                        dmaqf_dpco2(ispa,:) = dmaqf_dpco2(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*1d0./prox(:)'.^2d0;
                    elseif (ispa_c == 2)  % with HCO3- ( CO32- + H+)
                        maqf_loc(ispa,:) = maqf_loc(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*pco2x(:)'./prox(:)';
                        dmaqf_dpro(ispa,:) = dmaqf_dpro(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*pco2x(:)'*(-1d0)./prox(:)'.^2d0;
                        dmaqf_dpco2(ispa,:) = dmaqf_dpco2(ispa,:) + keqaq_c(ispa,ispa_c)*k1*k2*kco2*1d0./prox(:)';
                    end 
                end 
            end 
            % account for complexation with free SO4
            for ispa_s = 1:2
                rspa_s = double(ispa_s);
                if ( keqaq_s(ispa,ispa_s) > 0d0)  
                    maqf_loc(ispa,:) = maqf_loc(ispa,:) + keqaq_s(ispa,ispa_s)*so4f(:)'.^rspa_s;
                    dmaqf_dso4f(ispa,:) = dmaqf_dso4f(ispa,:) + keqaq_s(ispa,ispa_s)*rspa_s*so4f(:)'.^(rspa_s-1d0);
                end 
            end 
            % currently NO3 complexation with cations are ignored
            dmaqf_dpro(ispa,:) = maqx_loc(ispa,:)*(-1d0)./maqf_loc(ispa,:).^2d0.*dmaqf_dpro(ispa,:);
            dmaqf_dso4f(ispa,:) = maqx_loc(ispa,:)*(-1d0)/maqf_loc(ispa,:).^2d0.*dmaqf_dso4f(ispa,:);
            dmaqf_dpco2(ispa,:) = maqx_loc(ispa,:)*(-1d0)/maqf_loc(ispa,:).^2d0.*dmaqf_dpco2(ispa,:);
            dmaqf_dmaq(ispa,:) = 1d0./maqf_loc(ispa,:);
            maqf_loc(ispa,:) = maqx_loc(ispa,:)./maqf_loc(ispa,:);
        end 
    end     

end


function [ ...
    f1,df1,df12,df1dmaq,df1dmgas ...%output
    ,f2,df2,df21,df2dmaq,df2dmgas ...%output
    ] ...
    = calc_charge_so4_balance( ...
    nz,nsp_aq_all,nsp_gas_all ...
    ,chraq_all,chrgas_all ...
    ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
    ,base_charge ...
    ,mgasx_loc,maqf_loc ...
    ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
    ,z,prox,so4f,so4x ...
    ,print_loc,print_res,ph_add_order ...
    )

    [ieqgas_h0,ieqgas_h1,ieqgas_h2] = deal(1,2,3);

    % initialize output variables
    f1=zeros(nz,1,'double');df1=zeros(nz,1,'double');df12=zeros(nz,1,'double');f2=zeros(nz,1,'double');df2=zeros(nz,1,'double');df21=zeros(nz,1,'double');
    df1dmaq=zeros(nsp_aq_all,nz,'double');df2dmaq=zeros(nsp_aq_all,nz,'double');
    df1dmgas=zeros(nsp_gas_all,nz,'double');df2dmgas=zeros(nsp_gas_all,nz,'double');

    % initilize local variables
    ispa=0;ispa_h=0;ispa_c=0;ispa_s=0;iz=0;ipco2=0;ipnh3=0;
    kco2=0;k1=0;k2=0;knh3=0;k1nh3=0;rspa_h=0;rspa_s=0;ss_add=0;
    pco2x=zeros(nz,1,'double');pnh3x=zeros(nz,1,'double');
    f1_chk=zeros(nz,1,'double');

    % start
    if (print_res) 
        fid = fopen(print_loc,'w'); 
    end

    ipco2 = find(chrgas_all=='pco2');
    ipnh3 = find(chrgas_all=='pnh3');

    kco2 = keqgas_h(ipco2,ieqgas_h0);
    k1 = keqgas_h(ipco2,ieqgas_h1);
    k2 = keqgas_h(ipco2,ieqgas_h2);

    pco2x(:) = mgasx_loc(ipco2,:);

    knh3 = keqgas_h(ipnh3,ieqgas_h0);
    k1nh3 = keqgas_h(ipnh3,ieqgas_h1);

    pnh3x(:) = mgasx_loc(ipnh3,:);

    ss_add = ph_add_order;

    f1(:) = f1(:) + prox(:).^(ss_add+1d0) - kw*prox(:).^(ss_add-1d0);
    df1(:) = df1(:) + (ss_add+1d0)*prox(:).^ss_add - kw*(ss_add-1d0)*prox(:).^(ss_add-2d0);
    if (print_res); fprintf(fid,'%s\t%s\t%s\t', 'z','h', 'oh');  end

    % adding charges coming from aq species in eq with gases
    % pCO2
    f1(:) = f1(:)  -  k1*kco2*pco2x(:).*prox(:).^(ss_add-1d0)  -  2d0*k2*k1*kco2*pco2x(:).*prox(:).^(ss_add-2d0);
    df1(:) = df1(:)  -  k1*kco2*pco2x(:)*(ss_add-1d0).*prox(:).^(ss_add-2d0)  -  2d0*k2*k1*kco2*pco2x(:)*(ss_add-2d0).*prox(:).^(ss_add-3d0);
    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) -  k1*kco2*1d0*prox(:)'.^(ss_add-1d0)  -  2d0*k2*k1*kco2*1d0*prox(:)'.^(ss_add-2d0);
    if (print_res); fprintf(fid,'%s\t%s\t', 'hco3','co3');  end    
    % pNH3
    f1(:) = f1(:)  +  pnh3x(:)*knh3/k1nh3.*prox(:).^(ss_add+1d0);
    df1(:) = df1(:)  +  pnh3x(:)*knh3/k1nh3*(ss_add+1d0).*prox(:).^ss_add;
    df1dmgas(ipnh3,:) = df1dmgas(ipnh3,:)  +  1d0*knh3/k1nh3.*prox(:)'.^(ss_add+1d0);
    if (print_res) ; fprintf(fid,'%s\t', 'nh4');  end  

    %### SO4 mass balance ###

    f2(:) = so4x(:).*prox(:).^ss_add - so4f(:).*prox(:).^ss_add;
    df2(:) = - 1d0*prox(:).^ss_add;
    df21(:) = so4x(:)*ss_add.*prox(:).^(ss_add-1d0) - so4f(:)*ss_add.*prox(:).^(ss_add-1d0);
    df2dmaq(find(chraq_all=='so4'),:) =  1d0*prox(:).^ss_add;
    %### SO4 mass balance ###

    for ispa = 1: nsp_aq_all
        
        f1(:) = f1(:) + base_charge(ispa)*maqf_loc(ispa,:)'.*prox(:).^(ss_add);
        df1(:) = df1(:) + ( ...
            + base_charge(ispa)*dmaqf_dpro(ispa,:)'.*prox(:).^(ss_add)  ...
            + base_charge(ispa)*maqf_loc(ispa,:)'*(ss_add).*prox(:).^(ss_add-1d0)  ...
            );
        df12(:) = df12(:) + base_charge(ispa)*dmaqf_dso4f(ispa,:)'.*prox(:).^(ss_add);
        df1dmaq(ispa,:) = df1dmaq(ispa,:) + base_charge(ispa)*dmaqf_dmaq(ispa,:).*prox(:)'.^(ss_add); 
        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + base_charge(ispa)*dmaqf_dpco2(ispa,:) .*prox(:)'.^(ss_add);
        if (print_res) ; fprintf(fid,'%s\t', chraq_all(ispa));  end  
        
        % annions
        if (chraq_all(ispa)=='no3' || chraq_all(ispa)=='so4')  
            
            % account for speces associated with H+
            for ispa_h = 1:4
                if ( keqaq_h(ispa,ispa_h) > 0d0)  
                    rspa_h = double(ispa_h);
                    f1(:) = f1(:) + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'.*prox(:).^(rspa_h+ss_add);
                    df1(:) = df1(:) + ( ... 
                        + (base_charge(ispa) + rspa_h) ...
                               *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'*(rspa_h+ss_add).*prox(:).^(rspa_h+ss_add-1d0) ...
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpro(ispa,:)'.*prox(:).^(rspa_h+ss_add) ...
                        );
                    df12(:) = df12(:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dso4f(ispa,:)'.*prox(:).^(rspa_h+ss_add) ...
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dmaq(ispa,:).*prox(:)'.^(rspa_h+ss_add) ...
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpco2(ispa,:).*prox(:)'.^(rspa_h+ss_add) ...
                        );
                    if (print_res) 
                        fprintf(fid,'%s\t', strcat('h',string(ispa_h),chraq_all(ispa)));
                    end 
                    
                    % ### SO4 mass balance
                    % account for SO4 association with H+
                    if ( chraq_all(ispa)=='so4')  
                        f2(:) = f2(:) - keqaq_h(ispa,ispa_h)*so4f(:).*prox(:).^(rspa_h+ss_add);
                        df2(:) = df2(:) - keqaq_h(ispa,ispa_h)*1d0*prox(:).^(rspa_h+ss_add);
                        df21(:) = df21(:) - keqaq_h(ispa,ispa_h)*so4f(:)*(rspa_h+ss_add).*prox(:).^(rspa_h+ss_add-1d0);
                    end 
                    %### SO4 mass balance ###
                    
                end 
            end 
        % cations
        else 
            % account for hydrolysis speces
            for ispa_h = 1:4
                if ( keqaq_h(ispa,ispa_h) > 0d0)  
                    rspa_h = double(ispa_h);
                    f1(:) = f1(:) + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'.*prox(:).^(ss_add-rspa_h);
                    df1(:) = df1(:) + ( ...
                        + (base_charge(ispa) - rspa_h) ...
                              *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'*(ss_add-rspa_h).*prox(:).^(ss_add-rspa_h-1d0) ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpro(ispa,:)'.*prox(:).^(ss_add-rspa_h) ...
                        );
                    df12(:) = df12(:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dso4f(ispa,:)'.*prox(:).^(ss_add-rspa_h) ...
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dmaq(ispa,:).*prox(:)'.^(ss_add-rspa_h) ...
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpco2(ispa,:).*prox(:)'.^(ss_add-rspa_h) ...
                        );
                    if (print_res)  
                        fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(oh)',string(ispa_h)));
                    end 
                end 
            end 
            % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
            for ispa_c = 1:2
                if ( keqaq_c(ispa,ispa_c) > 0d0)  
                    if (ispa_c == 1)  % with CO3--
                        f1(:) = f1(:) + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0);
                        df1(:) = df1(:) + ( ... 
                            + (base_charge(ispa)-2d0) ...
                                  *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:)*(ss_add-2d0).*prox(:).^(ss_add-3d0) ...
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpro(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0) ...
                            );
                        df12(:) = df12(:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dso4f(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0) ...
                            );
                        df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dmaq(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-2d0) ...
                            );
                        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpco2(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-2d0) ...
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0.*prox(:)'.^(ss_add-2d0) ...
                            );
                        if (print_res) 
                            fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(co3)'));
                        end 
                    elseif (ispa_c == 2)  % with HCO3-
                        f1(:) = f1(:) + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0);
                        df1(:) = df1(:) + ( ... 
                            + (base_charge(ispa)-1d0) ...
                                  *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:)*(ss_add-1d0).*prox(:).^(ss_add-2d0) ...
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpro(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0) ...
                            );
                        df12(:) = df12(:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dso4f(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0) ...
                            );
                        df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dmaq(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-1d0) ...
                            );
                        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpco2(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-1d0) ...
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0.*prox(:)'.^(ss_add-1d0) ...
                            );
                        if (print_res) 
                            fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(hco3)'));
                        end 
                    end 
                end 
            end 
            % account for complexation with free SO4
            for ispa_s = 1:2
                if ( keqaq_s(ispa,ispa_s) > 0d0)  
                    rspa_s = double(ispa_s);
                    f1(:) = f1(:) + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add;
                    df1(:) = df1(:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dpro(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ... 
                        + (base_charge(ispa)-2d0*rspa_s) ...
                              *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s*ss_add.*prox(:).^(ss_add-1d0) ... 
                        );
                    df12(:) = df12(:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s) ...
                              *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'*rspa_s.*so4f(:).^(rspa_s-1d0).*prox(:).^ss_add ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dso4f(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ... 
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dmaq(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ... 
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dpco2(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ... 
                        );
                    if (print_res)  
                        fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(so4)',string(ispa_s)));
                    end 
                    % ### SO4 mass balance
                    % account for complexation with free SO4
                    f2(:) = f2(:) - rspa_s*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add;
                    df2(:) = df2(:) - ( ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'*rspa_s.*so4f(:).^(rspa_s-1d0).*prox(:).^ss_add ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dso4f(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ...
                        );
                    df21(:) = df21(:) - ( ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dpro(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s*ss_add.*prox(:).^(ss_add-1d0) ...
                        );
                    df2dmaq(ispa,:) = df2dmaq(ispa,:) - ( ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dmaq(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ...
                        );
                    df2dmgas(ipco2,:) = df2dmgas(ipco2,:) - ( ...
                        + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dpco2(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ...
                        );
                    %### SO4 mass balance ###
                        
                end 
            end 
            % currently NO3 complexation with cations are ignored
        end 
    end     

    if (print_res) ; fprintf(fid,'%s\n', 'tot_charge');  end  

    f1_chk(:) = 0d0;
    ss_add = 0d0;
    if (print_res) 
        for iz = 1: nz
            f1_chk(iz) = f1_chk(iz) + prox(iz)^(ss_add+1d0) - kw*prox(iz)^(ss_add-1d0);
            fprintf(fid,[repmat('%7.6e\t',1,3)], z(iz),prox(iz), kw/prox(iz)); 

            % adding charges coming from aq species in eq with gases
            % pCO2
            f1_chk(iz) = f1_chk(iz)  -  k1*kco2*pco2x(iz)*prox(iz)^(ss_add-1d0)  -  2d0*k2*k1*kco2*pco2x(iz)*prox(iz)^(ss_add-2d0);
            fprintf(fid,[repmat('%7.6e\t',1,2)], k1*kco2*pco2x(iz)/prox(iz),  k2*k1*kco2*pco2x(iz)/prox(iz)^2d0); 
            
            % pNH3
            f1_chk(iz) = f1_chk(iz)  +  pnh3x(iz)*knh3/k1nh3*prox(iz)^(ss_add+1d0);
            fprintf(fid,[repmat('%7.6e\t',1,1)], pnh3x(iz)*knh3/k1nh3*prox(iz)); 

            for ispa = 1: nsp_aq_all
                
                f1_chk(iz) = f1_chk(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*prox(iz)^(ss_add);
                fprintf(fid,[repmat('%7.6e\t',1,1)], maqf_loc(ispa,iz) ); 
                
                % annions
                if (chraq_all(ispa)=='no3' || chraq_all(ispa)=='so4')  
                    
                    % account for speces associated with H+
                    for ispa_h = 1:4
                        if ( keqaq_h(ispa,ispa_h) > 0d0)  
                            rspa_h = double(ispa_h);
                            f1_chk(iz) = f1_chk(iz) ...
                                + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^(rspa_h+ss_add);
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^rspa_h );
                        end 
                    end 
                % cations
                else 
                    % account for hydrolysis speces
                    for ispa_h = 1:4
                        if ( keqaq_h(ispa,ispa_h) > 0d0)  
                            rspa_h = double(ispa_h);
                            f1_chk(iz) = f1_chk(iz) ...
                                + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^(ss_add-rspa_h);
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)^rspa_h );
                        end 
                    end 
                    % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                    for ispa_c = 1:2
                        if ( keqaq_c(ispa,ispa_c) > 0d0)  
                            if (ispa_c == 1)  % with CO3--
                                f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-2d0) ...
                                    *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)^(ss_add-2d0);
                                fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                    ,keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz)^2d0 );
                                
                            elseif (ispa_c == 2)  % with HCO3-
                                f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-1d0) ...
                                    *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)^(ss_add-1d0);
                                fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                    ,keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz) );
                            end 
                        end 
                    end 
                    % account for complexation with free SO4
                    for ispa_s = 1:2
                        if ( keqaq_s(ispa,ispa_s) > 0d0)  
                            rspa_s = double(ispa_s);
                            f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)-2d0*rspa_s) ...
                                *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)^rspa_s*prox(iz)^ss_add;
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)^rspa_s );
                        end 
                    end 
                    % currently NO3 complexation with cations is ignored
                end 
            end     
            fprintf(fid,[repmat('%7.6e\n',1,1)],f1_chk(iz) ); 
        end 
        fclose(fid);
    end 

end


function [ ...
    df1,df12,df1dmaq,df1dmgas ...%output
    ,f1 ...% output
    ] = calc_charge( ...
    nz,nsp_aq_all,nsp_gas_all ...
    ,chraq_all,chrgas_all ...
    ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
    ,mgasx_loc,maqf_loc ...
    ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
    ,z,prox,so4f ...
    ,print_loc,print_res,ph_add_order ...
    )

    % fixed integers
    [ieqgas_h0,ieqgas_h1,ieqgas_h2] = deal(1,2,3);

    % output variable initilization
    f1=zeros(nz,1,'double');df1=zeros(nz,1,'double');df12=zeros(nz,1,'double');
    df1dmaq=zeros(nsp_aq_all,nz,'double');
    df1dmgas=zeros(nsp_gas_all,nz,'double');

    % local variables initilization
    ispa=0;ispa_h=0;ispa_c=0;ispa_s=0;iz=0;ipco2=0;ipnh3=0;
    kco2=0;k1=0;k2=0;knh3=0;k1nh3=0;rspa_h=0;rspa_s=0;ss_add=0;
    pco2x=zeros(nz,1,'double');pnh3x=zeros(nz,1,'double');
    f1_chk=zeros(nz,1,'double');
    base_charge=zeros(nsp_aq_all,1,'double');

    % -- start --- %


    if (print_res) 
        fid = fopen(print_loc,'w');
    end

    ipco2 = find(chrgas_all=='pco2');
    ipnh3 = find(chrgas_all=='pnh3');

    kco2 = keqgas_h(ipco2,ieqgas_h0);
    k1 = keqgas_h(ipco2,ieqgas_h1);
    k2 = keqgas_h(ipco2,ieqgas_h2);

    pco2x(:) = mgasx_loc(ipco2,:);

    knh3 = keqgas_h(ipnh3,ieqgas_h0);
    k1nh3 = keqgas_h(ipnh3,ieqgas_h1);

    pnh3x(:) = mgasx_loc(ipnh3,:);

    ss_add = ph_add_order;

    for ispa = 1: nsp_aq_all
        switch (chraq_all(ispa))
            case('so4')
                base_charge(ispa) = -2d0;
            case('no3')
                base_charge(ispa) = -1d0;
            case('si')
                base_charge(ispa) = 0d0;
            case{'na','k'}
                base_charge(ispa) = 1d0;
            case{'fe2','mg','ca'}
                base_charge(ispa) = 2d0;
            case{'fe3','al'}
                base_charge(ispa) = 3d0;
            otherwise
                error('error in charge assignment');
        end 
    end

    f1(:) = f1(:) + prox(:).^(ss_add+1d0) - kw*prox(:).^(ss_add-1d0);
    df1(:) = df1(:) + (ss_add+1d0)*prox(:).^ss_add - kw*(ss_add-1d0)*prox(:).^(ss_add-2d0);
    if (print_res); fprintf(fid,'%s\t%s\t%s\t', 'z','h', 'oh'); end

    % adding charges coming from aq species in eq with gases
    % pCO2
    f1(:) = f1(:)  -  k1*kco2*pco2x(:).*prox(:).^(ss_add-1d0)  -  2d0*k2*k1*kco2*pco2x(:).*prox(:).^(ss_add-2d0);
    df1(:) = df1(:)  -  k1*kco2*pco2x(:)*(ss_add-1d0).*prox(:).^(ss_add-2d0)  -  2d0*k2*k1*kco2*pco2x(:)*(ss_add-2d0).*prox(:).^(ss_add-3d0);
    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) -  k1*kco2*1d0*prox(:)'.^(ss_add-1d0)  -  2d0*k2*k1*kco2*1d0*prox(:)'.^(ss_add-2d0);
    if (print_res); fprintf(fid,'%s\t%s\t', 'hco3','co3'); end 
    % pNH3
    f1(:) = f1(:)  +  pnh3x(:)*knh3/k1nh3.*prox(:).^(ss_add+1d0);
    df1(:) = df1(:)  +  pnh3x(:)*knh3/k1nh3*(ss_add+1d0).*prox(:).^ss_add;
    df1dmgas(ipnh3,:) = df1dmgas(ipnh3,:)  +  1d0*knh3/k1nh3*prox(:)'.^(ss_add+1d0);
    if (print_res); fprintf(fid,'%s\t', 'nh4'); end  

    for ispa = 1: nsp_aq_all
        
        f1(:) = f1(:) + base_charge(ispa)*maqf_loc(ispa,:)'.*prox(:).^(ss_add);
        df1(:) = df1(:) + ( ...
            + base_charge(ispa)*dmaqf_dpro(ispa,:)'.*prox(:).^(ss_add)  ...
            + base_charge(ispa)*maqf_loc(ispa,:)'*(ss_add).*prox(:).^(ss_add-1d0)  ...
            );
        df12(:) = df12(:) + base_charge(ispa)*dmaqf_dso4f(ispa,:)'.*prox(:).^(ss_add); 
        df1dmaq(ispa,:) = df1dmaq(ispa,:) + base_charge(ispa)*dmaqf_dmaq(ispa,:).*prox(:)'.^(ss_add); 
        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + base_charge(ispa)*dmaqf_dpco2(ispa,:) .*prox(:)'.^(ss_add);
        if (print_res); fprintf(fid,'%s\t', chraq_all(ispa)); end  
        
        % annions
        if (chraq_all(ispa)=='no3' || chraq_all(ispa)=='so4')  
            
            % account for speces associated with H+
            for ispa_h = 1:4
                if ( keqaq_h(ispa,ispa_h) > 0d0)  
                    rspa_h = double(ispa_h);
                    f1(:) = f1(:) + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'.*prox(:).^(rspa_h+ss_add);
                    df1(:) = df1(:) + ( ... 
                        + (base_charge(ispa) + rspa_h) ...
                               *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'*(rspa_h+ss_add).*prox(:).^(rspa_h+ss_add-1d0) ...
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpro(ispa,:)'.*prox(:).^(rspa_h+ss_add) ...
                        );
                    df12(:) = df12(:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dso4f(ispa,:)'.*prox(:).^(rspa_h+ss_add) ...
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dmaq(ispa,:).*prox(:)'.^(rspa_h+ss_add) ...
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + (... 
                        + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpco2(ispa,:).*prox(:)'.^(rspa_h+ss_add) ...
                        );
                    if (print_res)  
                        fprintf(fid,'%s\t', strcat('h',string(ispa_h),chraq_all(ispa)));
                    end 
                end 
            end 
        % cations
        else 
            % account for hydrolysis speces
            for ispa_h = 1:4
                if ( keqaq_h(ispa,ispa_h) > 0d0)  
                    rspa_h = double(ispa_h);
                    f1(:) = f1(:) + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'.*prox(:).^(ss_add-rspa_h);
                    df1(:) = df1(:) + ( ...
                        + (base_charge(ispa) - rspa_h) ...
                              *keqaq_h(ispa,ispa_h)*maqf_loc(ispa,:)'*(ss_add-rspa_h).*prox(:).^(ss_add-rspa_h-1d0) ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpro(ispa,:)'.*prox(:).^(ss_add-rspa_h) ...
                        );
                    df12(:) = df12(:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dso4f(ispa,:)'.*prox(:).^(ss_add-rspa_h) ...
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dmaq(ispa,:).*prox(:)'.^(ss_add-rspa_h) ...
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ...
                        + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*dmaqf_dpco2(ispa,:).*prox(:)'.^(ss_add-rspa_h) ...
                        );
                    if (print_res)  
                        fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(oh)',string(ispa_h)));
                    end 
                end 
            end 
            % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
            for ispa_c = 1:2
                if ( keqaq_c(ispa,ispa_c) > 0d0)  
                    if (ispa_c == 1)  % with CO3--
                        f1(:) = f1(:) + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0);
                        df1(:) = df1(:) + ( ... 
                            + (base_charge(ispa)-2d0) ...
                                  *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:)*(ss_add-2d0).*prox(:).^(ss_add-3d0) ...
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpro(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0) ...
                            );
                        df12(:) = df12(:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dso4f(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-2d0) ...
                            );
                        df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dmaq(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-2d0) ...
                            );
                        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpco2(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-2d0) ...
                            + (base_charge(ispa)-2d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0.*prox(:)'.^(ss_add-2d0) ...
                            );
                        if (print_res) 
                            fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(co3)'));
                        end 
                    elseif (ispa_c == 2)  % with HCO3-
                        f1(:) = f1(:) + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0);
                        df1(:) = df1(:) + ( ... 
                            + (base_charge(ispa)-1d0) ...
                                  *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)'*k1*k2*kco2.*pco2x(:)*(ss_add-1d0).*prox(:).^(ss_add-2d0) ...
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpro(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0) ...
                            );
                        df12(:) = df12(:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dso4f(ispa,:)'*k1*k2*kco2.*pco2x(:).*prox(:).^(ss_add-1d0) ...
                            );
                        df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dmaq(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-1d0) ...
                            );
                        df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*dmaqf_dpco2(ispa,:)*k1*k2*kco2.*pco2x(:)'.*prox(:)'.^(ss_add-1d0) ...
                            + (base_charge(ispa)-1d0)*keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*1d0.*prox(:)'.^(ss_add-1d0) ...
                            );
                        if (print_res) 
                            fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(hco3)'));
                        end 
                    end 
                end 
            end 
            % account for complexation with free SO4
            for ispa_s = 1:2
                if ( keqaq_s(ispa,ispa_s) > 0d0)  
                    rspa_s = double(ispa_s);
                    f1(:) = f1(:) + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add;
                    df1(:) = df1(:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dpro(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ... 
                        + (base_charge(ispa)-2d0*rspa_s) ...
                              *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^rspa_s*ss_add.*prox(:).^(ss_add-1d0) ... 
                        );
                    df12(:) = df12(:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s) ...
                              *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'*rspa_s.*so4f(:).^(rspa_s-1d0).*prox(:).^ss_add ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dso4f(ispa,:)'.*so4f(:).^rspa_s.*prox(:).^ss_add ... 
                        );
                    df1dmaq(ispa,:) = df1dmaq(ispa,:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dmaq(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ... 
                        );
                    df1dmgas(ipco2,:) = df1dmgas(ipco2,:) + ( ... 
                        + (base_charge(ispa)-2d0*rspa_s)*keqaq_s(ispa,ispa_s)*dmaqf_dpco2(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ... 
                        );
                    if (print_res)  
                        fprintf(fid,'%s\t', strcat(chraq_all(ispa),'(so4)',string(ispa_s)));
                    end 
                end 
            end 
            % currently NO3 complexation with cations are ignored
        end 
    end     

    if (print_res) ; fprintf(fid,'%s\n', 'tot_charge'); end  

    f1_chk(:) = 0d0;
    ss_add = 0d0;
    if (print_res) 
        if (any(isnan(prox)) || any(prox <= 0d0)) 
            error (strcat('H+ conc is nan or <=0: showing H+ conc here: ',[repmat('%7.6e\t',1,nz)],'\n') ,prox(1:nz) ); 
        end 
        for iz = 1: nz
            f1_chk(iz) = f1_chk(iz) + prox(iz)^(ss_add+1d0) - kw*prox(iz)^(ss_add-1d0);
            fprintf(fid,[repmat('%7.6e\t',1,3)], z(iz),prox(iz), kw/prox(iz)); 

            % adding charges coming from aq species in eq with gases
            % pCO2
            f1_chk(iz) = f1_chk(iz)  -  k1*kco2*pco2x(iz)*prox(iz)^(ss_add-1d0)  -  2d0*k2*k1*kco2*pco2x(iz)*prox(iz)^(ss_add-2d0);
            fprintf(fid,[repmat('%7.6e\t',1,2)], k1*kco2*pco2x(iz)/prox(iz),  k2*k1*kco2*pco2x(iz)/prox(iz)^2d0); 
            
            % pNH3
            f1_chk(iz) = f1_chk(iz)  +  pnh3x(iz)*knh3/k1nh3*prox(iz)^(ss_add+1d0);
            fprintf(fid,[repmat('%7.6e\t',1,1)], pnh3x(iz)*knh3/k1nh3*prox(iz)); 
            
            for ispa = 1: nsp_aq_all
                
                f1_chk(iz) = f1_chk(iz) + base_charge(ispa)*maqf_loc(ispa,iz)*prox(iz)^(ss_add);
                fprintf(fid,[repmat('%7.6e\t',1,1)], maqf_loc(ispa,iz) ); 
                
                % annions
                if (chraq_all(ispa)=='no3' || chraq_all(ispa)=='so4')  
                    
                    % account for speces associated with H+
                    for ispa_h = 1:4
                        if ( keqaq_h(ispa,ispa_h) > 0d0)  
                            rspa_h = double(ispa_h);
                            f1_chk(iz) = f1_chk(iz) ...
                                + (base_charge(ispa) + rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^(rspa_h+ss_add);
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^rspa_h );
                        end 
                    end 
                % cations
                else 
                    % account for hydrolysis speces
                    for ispa_h = 1:4
                        if ( keqaq_h(ispa,ispa_h) > 0d0)  
                            rspa_h = double(ispa_h);
                            f1_chk(iz) = f1_chk(iz) ...
                                + (base_charge(ispa) - rspa_h)*keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)*prox(iz)^(ss_add-rspa_h);
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_h(ispa,ispa_h)*maqf_loc(ispa,iz)/prox(iz)^rspa_h );
                        end 
                    end 
                    % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                    for ispa_c = 1:2
                        if ( keqaq_c(ispa,ispa_c) > 0d0)  
                            if (ispa_c == 1)  % with CO3--
                                f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-2d0) ...
                                    *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)^(ss_add-2d0);
                                fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                    ,keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz)^2d0 );
                            elseif (ispa_c == 2)  % with HCO3-
                                f1_chk(iz) = f1_chk(iz) + (base_charge(ispa)-1d0) ...
                                    *keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)*prox(iz)^(ss_add-1d0);
                                fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                    ,keqaq_c(ispa,ispa_c)*maqf_loc(ispa,iz)*k1*k2*kco2*pco2x(iz)/prox(iz) );
                            end 
                        end 
                    end 
                    % account for complexation with free SO4
                    for ispa_s = 1:2
                        if ( keqaq_s(ispa,ispa_s) > 0d0)  
                            rspa_s = double(ispa_s);
                            f1_chk(iz) = f1_chk(iz)  + (base_charge(ispa)-2d0*rspa_s) ...
                                *keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)^rspa_s*prox(iz)^ss_add;
                            fprintf(fid,[repmat('%7.6e\t',1,1)] ...
                                ,keqaq_s(ispa,ispa_s)*maqf_loc(ispa,iz)*so4f(iz)^rspa_s );
                        end 
                    end 
                    % currently NO3 complexation with cations is ignored
                end 
            end     
            fprintf(fid,[repmat('%7.6e\n',1,1)],f1_chk(iz) ); 
        end 
        fclose(fid);
    end 

end


function [ ...
     df2,df21,df2dmaq,df2dmgas ...% output
     ,f2 ...% output
     ] = calc_so4_balance( ...
     nz,nsp_aq_all,nsp_gas_all ...
     ,chraq_all,chrgas_all ...
     ,keqaq_h,keqaq_s  ...
     ,maqf_loc ...
     ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
     ,prox,so4f,so4x ...
     ,ph_add_order ...
     )
    % f2 = prox^2d0*so4x - prox^2d0*so4f*( 1d0+k1so4*prox ...
        %  +k1kso4*kf ...
        %  +k1naso4*naf ...
        %  +k1caso4*caf ...
        %  +k1mgso4*mgf ...
        %  +k1fe2so4*fe2f ...
        %  +k1also4*alf ...
        %  +k1fe3so4*fe3f ...
        %  )

    % initilization of output variables
    f2=zeros(nz,1,'double');df2=zeros(nz,1,'double');df21=zeros(nz,1,'double');
    df2dmaq=zeros(nsp_aq_all,nz,'double');
    df2dmgas=zeros(nsp_gas_all,nz,'double');

    % initilization of local variables
    ispa=0;ispa_h=0;ispa_s=0;ipco2=0;
    rspa_h=0;rspa_s=0;ss_add=0;

    % -- start --
    ipco2 = find(chrgas_all=='pco2');

    ss_add = ph_add_order;

    f2(:) = so4x(:).*prox(:).^ss_add - so4f(:).*prox(:).^ss_add;
    df2(:) = - 1d0*prox(:).^ss_add;
    df21(:) = so4x(:)*ss_add.*prox(:).^(ss_add-1d0) - so4f(:)*ss_add.*prox(:).^(ss_add-1d0);
    df2dmaq(find(chraq_all=='so4'),:) =  1d0*prox(:)'.^ss_add ;

    f2(:) = 1d0;
    df2(:) = 0d0; 
    df21(:) = 0d0;

    for ispa = 1: nsp_aq_all
        
        switch (chraq_all(ispa))
            % annions
            case('so4')  
            
                % account for SO4 association with H+
                for ispa_h = 1:4
                    if ( keqaq_h(ispa,ispa_h) > 0d0)  
                        rspa_h = double(ispa_h);
                        f2(:) = f2(:) + keqaq_h(ispa,ispa_h)*prox(:).^(rspa_h);
                        df21(:) = df21(:) + keqaq_h(ispa,ispa_h)*rspa_h*prox(:).^(rspa_h-1d0);
                    end 
                end 
            
            case('no3')
                % do nothing because it is not associated with so4
            
            % cations
            case{'k','na','si','mg','ca','fe2','fe3','al'}  
                % account for complexation with free SO4
                for ispa_s = 1:2
                    if ( keqaq_s(ispa,ispa_s) > 0d0) 
                        rspa_s = double(ispa_s);
                        f2(:) = f2(:) + rspa_s*keqaq_s(ispa,ispa_s)*maqf_loc(ispa,:)'.*so4f(:).^(rspa_s-1d0);
                        df2(:) = df2(:) + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dso4f(ispa,:)'.*so4f(:).^(rspa_s-1d0);
                        df21(:) = df21(:) + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dpro(ispa,:)'.*so4f(:).^(rspa_s-1d0);
                        df2dmaq(ispa,:) = df2dmaq(ispa,:) - ( ...
                             + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dmaq(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ...
                             );
                        df2dmgas(ipco2,:) = df2dmgas(ipco2,:) - ( ...
                             + rspa_s*keqaq_s(ispa,ispa_s)*dmaqf_dpco2(ispa,:).*so4f(:)'.^rspa_s.*prox(:)'.^ss_add ...
                             );
                    end 
                end 
            
            otherwise
                
                error('** error: you should not come here @ calc_so4_balance')
            
        end
    end     

    df2(:) = - 1d0*prox(:).^ss_add.*f2(:) - so4f(:).*prox(:).^ss_add.*df2(:);
    df21(:) = so4x(:)*ss_add.*prox(:).^(ss_add-1d0) - so4f(:)*ss_add.*prox(:).^(ss_add-1d0).*f2(:) - so4f(:).*prox(:).^ss_add.*df21(:);
    f2(:) = so4x(:).*prox(:).^ss_add - so4f(:).*prox(:).^ss_add.*f2(:);

end


function [ ...
    dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
    ,prox,ph_error,so4f,ph_iter ...% output
    ] = calc_pH_v7_3( ...
    nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
    ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
    ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
    ,print_cb,print_loc,z ...% input 
    ,prox,so4f ...% inout
    ) 
    % solving charge balance and so4 species balance

    % initialize output variables
    ph_error = false;
    ph_iter=0;
    dprodmaq_all=zeros(nsp_aq_all,nz,'double');dso4fdmaq_all=zeros(nsp_aq_all,nz,'double');
    dprodmgas_all=zeros(nsp_gas_all,nz,'double');dso4fdmgas_all=zeros(nsp_gas_all,nz,'double');

    % local variables
    so4th=0;
    so4x=zeros(nz,1,'double');

    df1=zeros(nz,1,'double');f1=zeros(nz,1,'double');f2=zeros(nz,1,'double');
    df2=zeros(nz,1,'double');df21=zeros(nz,1,'double');df12=zeros(nz,1,'double');
    error=0;tol=0;dconc=0;ph_add_order =0;
    iter=0;iz=0;ispa=0;ispg=0;

    base_charge=zeros(nsp_aq_all,1,'double');
    maqx_loc=zeros(nsp_aq_all,nz,'double');maqf_loc=zeros(nsp_aq_all,nz,'double');
    dmaqf_dpro=zeros(nsp_aq_all,nz,'double');dmaqf_dso4f=zeros(nsp_aq_all,nz,'double');
    dmaqf_dmaq=zeros(nsp_aq_all,nz,'double');dmaqf_dpco2=zeros(nsp_aq_all,nz,'double');
    mgasx_loc=zeros(nsp_gas_all,nz,'double');
    df1dmaq=zeros(nsp_aq_all,nz,'double');df2dmaq=zeros(nsp_aq_all,nz,'double');
    df1dmgas=zeros(nsp_gas_all,nz,'double');df2dmgas=zeros(nsp_gas_all,nz,'double');

    dmaq=zeros(nsp_aq_all,nz,'double');
    dmgas=zeros(nsp_gas_all,nz,'double');
    df1_dum=zeros(nz,1,'double');f1_dum=zeros(nz,1,'double');f2_dum=zeros(nz,1,'double');
    df2_dum=zeros(nz,1,'double');df21_dum=zeros(nz,1,'double');df12_dum=zeros(nz,1,'double');
    df1dmaq_dum=zeros(nsp_aq_all,nz,'double');df2dmaq_dum=zeros(nsp_aq_all,nz,'double');
    df1dmgas_dum=zeros(nsp_gas_all,nz,'double');df2dmgas_dum=zeros(nsp_gas_all,nz,'double');

    so4_error=false;print_res=false;


    % start doing something

    error = 1d4;
    tol = 1d-6;
    dconc = 1d-9;
    ph_add_order = 2d0;

    iter = 0;

    if (any(isnan(maqx(:))) || any(isnan(maqc(:))))  
        error ('nan in input aqueosu species')
    end 

    [maqx_loc,mgasx_loc] = get_maqgasx_all( ...
        nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst ...
        ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst ...
        ,maqx,mgasx,maqc,mgasc ...
        );
        
    [base_charge] = get_base_charge( ...
        nsp_aq_all ... 
        ,chraq_all ... 
        );

    so4x(:) = maqx_loc(find(chraq_all=='so4'),:);

    so4th = maqth_all(find(chraq_all=='so4'));

    % so4f  = so4x 

    nmx = 2*nz;
    % if (all(so4x==0d0))  
    if (all(so4x<=so4th))  
        nmx = nz;
    end  

    amx = zeros(nmx,nmx,'double');
    ymx = zeros(nmx,1,'double');
    xmx = zeros(nmx,1,'double');
    rmx = 0;


    if (~ print_cb) 
        % obtaining ph and so4f from scratch
      
        so4f(:) = so4x(:);
        prox(:) = 1d0; 
        while (error > tol)
            % free SO42- (for simplicity only consider XSO4 complex where X is a cation)
            
            [ ...
                dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
                ,maqf_loc  ...% output
                ] = get_maqf_all( ...
                nz,nsp_aq_all,nsp_gas_all ...
                ,chraq_all,chrgas_all ...
                ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
                ,mgasx_loc,maqx_loc,prox,so4f ...
                );
            
            [ ...
                f1,df1,df12,df1dmaq,df1dmgas ...%output
                ,f2,df2,df21,df2dmaq,df2dmgas ...%output
                ] = calc_charge_so4_balance( ...
                nz,nsp_aq_all,nsp_gas_all ...
                ,chraq_all,chrgas_all ...
                ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
                ,base_charge ...
                ,mgasx_loc,maqf_loc ...
                ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
                ,z,prox,so4f,so4x ...
                ,print_loc,print_res,ph_add_order ...
                );
            
            df1(:) = df1(:).*prox(:);
            df21(:) = df21(:).*prox(:);
            df2(:) = df2(:).*so4f(:);
            df12(:) = df12(:).*so4f(:);
            
            if (any(isnan(f1))||any(isnan(f2))||any(isnan(df1))||any(isnan(df2)) ...
                ||any(isnan(df12))||any(isnan(df21)))  
                warning(strcat('found nan during the course of ph calc:\n' ...
                    ,'NAN in f1 = %s\t' ,'NAN in f2 = %s\t','NAN in df1 = %s\n' ...
                    ,'NAN in df2 = %s\t','NAN in df12 = %s\t','NAN in df21 = %s\n') ...
                    ,string(any(isnan(f1))),string(any(isnan(f2))),string(any(isnan(df1))) ...
                    ,string(any(isnan(df2))),string(any(isnan(df12))),string(any(isnan(df21))) ...
                    );
                fprintf( strcat('... displaying ph ...',[repmat('%7.6e\t',1,nz)],'\n') , prox(:) );
                ph_error = true;
                return; 
            end 
            
            
            if (nmx~=nz)  
                amx(:,:) = 0d0;
                ymx(:) = 0d0;
                
                ymx(1:nz) = f1(:);
                ymx(nz+1:nmx) = f2(:);
                
                for iz=1:nz
                    amx(iz,iz)=df1(iz);
                    amx(nz+iz,nz+iz)=df2(iz);
                    amx(iz,nz+iz)=df12(iz);
                    amx(nz+iz,iz)=df21(iz);
                end 
                ymx(:) = -ymx(:);
                
                % call DGESV(nmx,int(1),amx,nmx,ipiv,ymx,nmx,info) 
                [xmx,rmx] = linsolve(amx,ymx);
                ymx = xmx;
                
                prox(:) = prox(:).*exp( ymx(1:nz) );
                so4f(:) = so4f(:).*exp( ymx(nz+1:nmx) );
                
                error = max(abs(exp( ymx )-1d0));
                if (isnan(error))  
                    error = 1d4;
                    ph_error = true;
                    return; 
                end 
            else 
                prox(:) = prox(:).*exp( -f1(:)./df1(:) );
                error = max(abs(exp( -f1(:)./df1(:) )-1d0));
                if (isnan(error)); error = 1d4; end
            end 
            
            iter = iter + 1;
            
            % fprintf('%d\t%7.6e\n',iter,error);
            
            if (iter > 3000)  
                warning ('iteration exceeds 3000');
                ph_error = true;
                return
            end 
        end  
    end 

    ph_iter = iter;
    
    % pause;

    if (any(isnan(prox)) || any(prox<=0d0))      
        warning(strcat('ph is nan or <= zero: showing pH here: ',[repmat('%7.6e\t',1,5)],'\n') ... 
            ,-log10(prox(1:nz/5:nz)));
        ph_error = true;
    end 

    [...
        dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
        ,maqf_loc  ...% output
        ] = get_maqf_all( ...
        nz,nsp_aq_all,nsp_gas_all ...
        ,chraq_all,chrgas_all ...
        ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
        ,mgasx_loc,maqx_loc,prox,so4f ...
        );

    if (print_cb); print_res = true; end

    [ ...
        df1,df12,df1dmaq,df1dmgas ...%output
        ,f1 ...% output
        ] = calc_charge( ...
        nz,nsp_aq_all,nsp_gas_all ...
        ,chraq_all,chrgas_all ...
        ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
        ,mgasx_loc,maqf_loc ...
        ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
        ,z,prox,so4f ...
        ,print_loc,print_res,ph_add_order ...
        );

        
    [ ...
        df2,df21,df2dmaq,df2dmgas ...% output
        ,f2 ...% output
        ] = calc_so4_balance( ...
        nz,nsp_aq_all,nsp_gas_all ...
        ,chraq_all,chrgas_all ...
        ,keqaq_h,keqaq_s  ...
        ,maqf_loc ...
        ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
        ,prox,so4f,so4x ...
        ,ph_add_order ...
        );    
        
    for ispa=1:nsp_aq_all
        dmaq(:,:) = 0d0;
        dmaq(ispa,:) = dconc; %*maqx_loc(ispa,:)
        [ ...
            dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
            ,maqf_loc  ...% output
            ] = get_maqf_all( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
            ,mgasx_loc,maqx_loc+dmaq,prox,so4f ...
            );
            
        [ ...
            df1_dum,df12_dum,df1dmaq_dum,df1dmgas_dum ...%output
            ,f1_dum ...% output
            ] = calc_charge( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
            ,mgasx_loc,maqf_loc ...
            ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
            ,z,prox,so4f ...
            ,print_loc,print_res,ph_add_order ...
            );
        
        [ ...
            df2_dum,df21_dum,df2dmaq_dum,df2dmgas_dum ...% output
            ,f2_dum ...% output
            ] = calc_so4_balance( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,keqaq_h,keqaq_s  ...
            ,maqf_loc ...
            ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
            ,prox,so4f,so4x ...
            ,ph_add_order ...
            );   
        dprodmaq_all(ispa,:) = -f1_dum(:)'./df1_dum(:)'/dconc; %/maqx_loc(ispa,:) 
        dso4fdmaq_all(ispa,:) = -f2_dum(:)'./df2_dum(:)'/dconc;%/maqx_loc(ispa,:) 
    end 
        
    for ispg=1:nsp_gas_all
        dmgas(:,:) = 0d0;
        dmgas(ispg,:) = dconc; %*mgasx_loc(ispg,:)
        [ ...
            dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
            ,maqf_loc  ...% output
            ] = get_maqf_all( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
            ,mgasx_loc+dmgas,maqx_loc,prox,so4f ...
            );
            
        [ ...
            df1_dum,df12_dum,df1dmaq_dum,df1dmgas_dum ...%output
            ,f1_dum ...% output
            ] = calc_charge( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,kw,keqgas_h,keqaq_h,keqaq_c,keqaq_s  ...
            ,mgasx_loc+dmgas,maqf_loc ...
            ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
            ,z,prox,so4f ...
            ,print_loc,print_res,ph_add_order ...
            );
        
        [ ... 
            df2_dum,df21_dum,df2dmaq_dum,df2dmgas_dum ...% output
            ,f2_dum ...% output
            ] = calc_so4_balance( ...
            nz,nsp_aq_all,nsp_gas_all ...
            ,chraq_all,chrgas_all ...
            ,keqaq_h,keqaq_s  ...
            ,maqf_loc ...
            ,dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...
            ,prox,so4f,so4x ...
            ,ph_add_order ...
            );   
        dprodmgas_all(ispg,:) = -f1_dum(:)'./df1_dum(:)'/dconc;%/mgasx_loc(ispg,:)  
        dso4fdmgas_all(ispg,:) = -f2_dum(:)'./df2_dum(:)'/dconc; %/mgasx_loc(ispg,:)
    end 
    % 
    % solving two equations analytically:
    % df1/dph * dph/dmsp + df1/dso4f * dso4f/dmsp + df1/dmsp = 0 
    % df2/dph * dph/dmsp + df2/dso4f * dso4f/dmsp + df2/dmsp = 0 
    % using the variables in this subroutine and defining x = dph/dmsp and y = dso4f/dmsp
    % df1 * x + df12 * y + df1dmsp = 0
    % df21 * x + df2 * y + df2dmsp = 0
    % 
    for ispa = 1: nsp_aq_all
        dprodmaq_all(ispa,:) = - (df2(:)'.*df1dmaq(ispa,:) - df12(:)'.*df2dmaq(ispa,:))./(df2(:)'.*df1(:)' - df12(:)'.*df21(:)');   
        dso4fdmaq_all(ispa,:) = - ( df21(:)'.*df1dmaq(ispa,:) - df1(:)'.*df2dmaq(ispa,:) )./(df21(:)'.*df12(:)' - df1(:)'.*df2(:)' ); 
    end 

    for ispg = 1: nsp_gas_all
        dprodmgas_all(ispg,:) = - (df2(:)'.*df1dmgas(ispg,:) - df12(:)'.*df2dmgas(ispg,:) )./(df2(:)'.*df1(:)' -df12(:)'.*df21(:)');
        dso4fdmgas_all(ispg,:) = - ( df21(:)'.*df1dmgas(ispg,:) - df1(:)'.*df2dmgas(ispg,:) )./(df21(:)'.*df12(:)' - df1(:)'.*df2(:)' );  
    end 

end


function [n_tmp] = Console4(file_name)

    A = importdata(file_name);
    if isa(A,'struct')
        n_tmp = size(A.data,1);
    else
        n_tmp = size(A,1);
        n_tmp = n_tmp - 1;
    end
end 

function [nsp_sld,nsp_aq,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_save] ...
    = get_saved_variables_num(workdir,runname_save)

    file_name = strcat (workdir,runname_save,'/slds.save');
    [nsp_sld] = Console4(file_name);
    file_name = strcat(workdir,runname_save,'/solutes.save');
    [nsp_aq] = Console4(file_name);
    file_name = strcat(workdir,runname_save,'/gases.save');
    [nsp_gas] = Console4(file_name);
    file_name = strcat(workdir,runname_save,'/extrxns.save');
    [nrxn_ext] = Console4(file_name);
    file_name = strcat(workdir,runname_save,'/kinspc.save');
    [nsld_kinspc] = Console4(file_name);
    file_name = strcat(workdir,runname_save,'/sa.save');
    [nsld_sa_save] = Console4(file_name);

end


function [ ...
    chraq,chrgas,chrsld,chrrxn_ext,chrsld_kinspc,kin_sld_spc ...% output
    ,chrsld_sa_dum,hrii_dum ...% output 
    ] = get_saved_variables( ...
    workdir,runname_save ...% input 
    ,nsp_aq,nsp_sld,nsp_gas,nrxn_ext,nsld_kinspc,nsld_sa_dum ...% input
    )

    if (nsp_aq>=1)  
        file_name = srcat(workdir,runname_save,'/solutes.save');
        chraq = strings(nsp_aq,1);
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chraq(cnt) = chr_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end

    if (nsp_sld>=1)  
        file_name = srcat(workdir,runname_save,'/slds.save');
        chrsld = strings(nsp_sld,1);
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chrsld(cnt) = chr_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end

    if (nsp_gas>=1)  
        file_name = srcat(workdir,runname_save,'/gases.save');
        chrgas = strings(nsp_gas,1);
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chrgas(cnt) = chr_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end

    if (nrxn_ext>=1)  
        file_name = srcat(workdir,runname_save,'/extrxns.save');
        chrrxn_ext = strings(nrxn_ext,1);
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chrrxn_ext(cnt) = chr_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end

    if (nsld_kinspc>=1)  
        file_name = srcat(workdir,runname_save,'/kinspc.save');
        chrsld_kinspc = strings(nsld_kinspc,1);
        kin_sld_spc = zeros(nsld_kinspc,1,'double');
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chrsld_kinspc(cnt) = chr_tmp; 
            
            val_tmp = str2num(char(a(2)));
            kin_sld_spc(cnt) = val_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end

    if (nsld_sa_dum>=1)  
        file_name = srcat(workdir,runname_save,'/sa.save');
        chrsld_sa_dum = strings(nsld_sa_dum,1);
        hrii_dum = zeros(nsld_sa_dum,1,'double');
        tline = fgetl(fid);
        cnt = 0;
        while ischar(tline)
            % disp(tline)
            a = strsplit(tline);
            tline = fgetl(fid);
            
            if cnt ==0
                cnt = cnt + 1;
                continue
            end
            
            chr_tmp = (char(a(1)));
            chrsld_sa_dum(cnt) = chr_tmp; 
            
            val_tmp = str2num(char(a(2)));
            hrii_dum(cnt) = val_tmp; 
            
            cnt = cnt + 1;
                
        end
        fclose(fid);
    end
end


function [trans,nonlocal,izml] = make_transmx(  ...
    labs,nsp_sld,turbo2,nobio,dz,poro,nz,z,zml_ref,dbl_ref,fick,till,tol,save_trans  ...% input
    )
    
    % initilization of output variables
    trans = zeros(nz,nz,nsp_sld,'double');
    nonlocal = zeros(nsp_sld,1,'logical');
    izml = 0;

    % initilization of local variables
    sporo = zeros(nz,1,'double');dbio = zeros(nz,1,'double'); 
    iz=0;isp=0;iiz=0;izdbl=0;
    translabs = zeros(nz,nz,'double'); transdbio = zeros(nz,nz,'double'); 
    transturbo2 = zeros(nz,nz,'double');transtill = zeros(nz,nz,'double');
    zml = zeros(nsp_sld,1,'double');probh=0;dbl=0;

    % --- start ---% 

    sporo(:) = 1d0 - poro(:);
    trans(:,:,:) = 0d0;
    %~~~~~~~~~~~~ loading transition matrix from LABS ~~~~~~~~~~~~~~~~~~~~~~~~
    if (any(labs)) 
        translabs(:,:) = 0d0;
        % to be implemented ...?
    end

    if (true)   % devided by the time duration when transition matrices are created in LABS and weakening by a factor
        translabs = translabs *365.25d0/10d0*1d0/10d0;
    end
    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    zml(:)=zml_ref; % mixed layer depth assumed to be a reference value at first 

    dbl = dbl_ref;

    nonlocal(:) = false; % initial assumption 
    for isp=1:nsp_sld
        if (turbo2(isp) || labs(isp)); nonlocal(isp)=true; end % if mixing is made by turbo2 or labs,  nonlocal 
        
        dbio(:)=0d0;
        izdbl=0;
        for iz = 1: nz
            if (z(iz) <= dbl)  
                dbio(iz) = 0d0;
                izdbl = iz;
            elseif (dbl < z(iz) && z(iz) <=zml(isp)) 
                % dbio(iz) =  0.15d-4 ;  %  within mixed layer 150 cm2/kyr (Emerson, 1985) 
                dbio(iz) =  2d-4;   %  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                % dbio(iz) =  2d-4*exp(z(iz)/0.1d0);   %  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                % dbio(iz) =  2d-7*exp(z(iz)/1d0);   %  within mixed layer ~5-6e-7 m2/day (Astete et al., 2016) 
                % dbio(iz) =  2d-10;   %  just a small value 
                % dbio(iz) =  2d-3;   %  just a value changed 
                izml = iz;   % determine grid of bottom of mixed layer 
            else
                dbio(iz) =  0d0; % no biodiffusion in deeper depths 
            end
        end

        transdbio(:,:) = 0d0;   % transition matrix to realize Fickian mixing with biodiffusion coefficient dbio which is defined just above 
        for iz = max(1,izdbl): izml
            if (iz==max(1,izdbl))  
                transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
            elseif (iz==izml)  
                transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
                transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
            else 
                transdbio(iz,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz-1)))  ...
                    + 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(-1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
                transdbio(iz-1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz-1)*dbio(iz-1))*(1d0)/(0.5d0*(dz(iz)+dz(iz-1)));
                transdbio(iz+1,iz) = 0.5d0*(sporo(iz)*dbio(iz)+sporo(iz+1)*dbio(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(iz+1)));
            end
        end
        
        % Added; changes have been made here rather than in solving governing eqs.
        for iz=1:nz
            transdbio(:,iz) = transdbio(:,iz)/dz(iz);
        end 

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
        % transition matrix for random mixing 
        transturbo2(:,:) = 0d0;
        % ending up in upward mixing 
        probh = 0.0010d0;
        transturbo2(max(1,izdbl):izml,max(1,izdbl):izml) = probh;  % arbitrary assumed probability 
        for iz=1:izml  % when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
           transturbo2(iz,iz)=-probh*(izml-max(1,izdbl));  
        end
        % trying real homogeneous 
        transturbo2(:,:) = 0d0;
        probh = 0.001d0; % def used in IMP
        % probh = 0.01d0; % strong mixing
        probh = 0.1d0; % strong mixing
        % probh = 0.5d0; % strong mixing
        % probh = 0.0005d0; % just testing smaller mixing (used for tuning)
        % probh = 0.0001d0; % just testing smaller mixing for PSDs
        for iz=1:izml 
            for iiz=1:izml
                if (iiz~=iz)  
                    transturbo2(iiz,iz) = probh;%*dz(iz)/dz(iiz)
                    transturbo2(iiz,iiz) = transturbo2(iiz,iiz) - transturbo2(iiz,iz);
                end 
            end
        end
        
        % trying inverse mixing 
        transtill(:,:) = 0d0;
        probh = 0.010d0;
        probh = 0.10d0;
        for iz=1:izml  % when i = j, transition matrix contains probabilities with which particles are moved from other layers of sediment   
            % transtill(iz,iz)=-probh*dz(iz)/dz(izml+1-iz) %*(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
            % transtill(izml+1-iz,iz)=probh*dz(iz)/dz(izml+1-iz) % *(iz - izml*0.5d0)**2d0/(izml**2d0*0.25d0)
            % for iiz = izml+1-iz,iz+1,-1
                % transtill(iiz,iz)= probh*dz(iz)/dz(iiz)*(iiz/real(izml+1-iz,kind=8))
            % end 
            % transtill(iz,iz) = -sum(transtill(:,iz))
            for iiz=1:izml
                if (iiz~=iz)  
                    if (iiz==iz-1 || iiz == iz+1 || iiz == izml + 1 - iz )  
                        transtill(iiz,iz)= probh; %*dz(iz)/dz(iiz) 
                        % transtill(iiz,iiz) = transtill(iiz,iiz) - transtill(iiz,iz)
                    end 
                end 
            end 
            transtill(iz,iz) = -sum(transtill(:,iz));
        end
        

        % if (turbo2(isp)) translabs = transturbo2   % translabs temporarily used to represents nonlocal mixing 
        
        % added 
        for iz =1:nz
            for iiz= 1:nz
                translabs(iiz,iz) = translabs(iiz,iz)/dz(iz)*dz(iiz);
                transturbo2(iiz,iz) = transturbo2(iiz,iz)/dz(iz)*dz(iiz);
                transtill(iiz,iz) = transtill(iiz,iz)/dz(iz)*dz(iiz);
            end 
        end 
        
        trans(:,:,isp) = 0d0; 
        
        if (nobio(isp)); continue; end
        
        if (fick(isp))  
            trans(:,:,isp) = trans(:,:,isp) + transdbio(:,:);
        end 
        
        if (turbo2(isp))  
            trans(:,:,isp) = trans(:,:,isp) + transturbo2(:,:);
        end 
        
        if (labs(isp))  
            trans(:,:,isp) = trans(:,:,isp) + translabs(:,:);
        end 
        
        if (till(isp))  
            trans(:,:,isp) = trans(:,:,isp) + transtill(:,:);
        end 
        
        if (save_trans)  
            fid = fopen(strcat('./mtx-',string(isp),'.txt'),'w');
            fmt=[repmat('%7.6E \t',1,nz) '\n'];
            for iz=1:nz
                fprintf(fid,fmt,trans(iz,1:nz,isp));  % writing 
            end
            fclose(fid);
        end 
    end
    % even when all are local Fickian mixing, mixing treatment must be the same as in case of nonlocal 
    % if mixing intensity and depths are different between different species  
    if (all(~nonlocal))   
        for isp=1:nsp_sld-1
            if (any(trans(:,:,isp+1)~=trans(:,:,isp),'all')); nonlocal(:)=true; end
        end
    end 

end


function res = merge(a,b,c)

    % res = a.*c + b;
    
    res = a.*c + b.*(~c);
    
    % if size(c,1)==1
        % if c
            % res = a;
        % else 
            % res = b;
        % end 
    % else
        % res = zeros(size(c,1),1)
        % for i =1:size(c,1)
            % if c(i) 
                % res(i)=a;
            % else
                % res(i)=b;
            % end
        % end 
    % end 
end


function [ ...
    domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dso4f_loc ...% output
    ,omega,omega_error ...% output
    ] = calc_omega_v4( ...
    nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ... 
    ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all ...
    ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,so4f ...
    ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
    ,staq_all,stgas_all ...
    ,mineral ...
    )

    % output variables initialization 
    omega = zeros(nz,1,'double');
    omega_error = false;
    domega_dpro_loc = zeros(nz,1,'double');domega_dso4f_loc = zeros(nz,1,'double');
    domega_dmgas_all = zeros(nsp_gas_all,nz,'double');
    domega_dmaq_all = zeros(nsp_aq_all,nz,'double');

    % local variables initialization 
    k1=0;k2=0;kco2=0;po2th=0;mo2g1=0;mo2g2=0;mo2g3=0;keq_tmp=0;ss_x=0;ss_pro=0;ss_pco2=0;mo2_tmp=0;
    pco2x=zeros(nz,1,'double');po2x=zeros(nz,1,'double');
    maqx_loc=zeros(nsp_aq_all,nz,'double');maqf_loc=zeros(nsp_aq_all,nz,'double');
    dmaqf_dpro=zeros(nsp_aq_all,nz,'double');dmaqf_dso4f=zeros(nsp_aq_all,nz,'double');
    dmaqf_dmaq=zeros(nsp_aq_all,nz,'double');dmaqf_dpco2=zeros(nsp_aq_all,nz,'double');
    mgasx_loc=zeros(nsp_gas_all,nz,'double');

    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);
    [ieqaq_co3,ieqaq_hco3]=deal(1,2);
    [ieqaq_so4,ieqaq_so42]=deal(1,2);

    ispa=0;ipco2=0;ipo2=0;
    % thon = 1d0;
    thon = -1d100;

    % ---- start ----

    mo2g1 = keqsld_all(find(chrsld_all=='g1'));
    mo2g2 = keqsld_all(find(chrsld_all=='g2'));
    mo2g3 = keqsld_all(find(chrsld_all=='g3'));

    po2th = mgasth_all(find(chrgas_all=='po2'));

    kco2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0);
    k1 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);
    k2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2);

    [maqx_loc,mgasx_loc] = get_maqgasx_all( ...
        nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst ...
        ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst ...
        ,maqx,mgasx,maqc,mgasc ...
        );

    [ ...
        dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
        ,maqf_loc  ...% output
        ] = get_maqf_all( ...
        nz,nsp_aq_all,nsp_gas_all ...
        ,chraq_all,chrgas_all ...
        ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
        ,mgasx_loc,maqx_loc,prox,so4f ...
        );


    pco2x(:) = mgasx_loc(find(chrgas_all=='pco2'),:);
    po2x(:) = mgasx_loc(find(chrgas_all=='po2'),:);

    ipco2 = find(chrgas_all=='pco2');
    ipo2 = find(chrgas_all=='po2');

    domega_dmaq_all(:,:) =0d0;
    domega_dmgas_all(:,:) =0d0;
    domega_dso4f_loc(:) =0d0;
    domega_dpro_loc(:) =0d0;

    switch (mineral)

        % case default % (almino)silicates & oxides
        case { ...
            'fo','ab','an','ka','gb','ct','fa','gt','cabd','dp','hb','kfs','amsi','hm','ill','anl','nph' ...
            ,'qtz','tm','la','by','olg','and','cpx','en','fer','opx','mgbd','kbd','nabd','mscv','plgp','antp' ...
            ,'agt' ...
            }  % (almino)silicates & oxides
            keq_tmp = keqsld_all(find(chrsld_all==mineral));
            omega(:) = 1d0;
            ss_pro = 0d0;
            for ispa = 1:nsp_aq_all
                if (staq_all(find(chrsld_all==mineral),ispa) > 0d0)  

                    omega(:) = omega(:).*maqf_loc(ispa,:)'.^staq_all(find(chrsld_all==mineral),ispa);
                    
                    % derivatives are first given as d(log omega)/dc 
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dmaq(ispa,:) ...
                        );
                    domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dpco2(ispa,:) ...
                        );
                    domega_dpro_loc(:) = domega_dpro_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dpro(ispa,:)' ...
                        );
                    domega_dso4f_loc(:) = domega_dso4f_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dso4f(ispa,:)' ...
                        );

                    switch (chraq_all(ispa))
                        case{'na','k'}
                            ss_pro = ss_pro + staq_all(find(chrsld_all==mineral),ispa);
                        case{'fe2','ca','mg'}
                            ss_pro = ss_pro + 2d0*staq_all(find(chrsld_all==mineral),ispa);
                        case{'fe3','al'}
                            ss_pro = ss_pro + 3d0*staq_all(find(chrsld_all==mineral),ispa);
                    end
                end 
            end 
            
            if (ss_pro > 0d0)  
                omega(:) = omega(:) ./ prox(:).^ss_pro;

                % derivatives are first given as d(log omega)/dc 
                domega_dpro_loc(:) = domega_dpro_loc(:) - ss_pro./prox(:); 
            end 
            
            if (keq_tmp > 0d0)  
                omega(:) = omega(:) / keq_tmp;
            end         
            
            % derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
            for ispa=1:nsp_aq_all
                domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:).*omega(:)';
            end 
            domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:).*omega(:)';
            domega_dpro_loc(:) = domega_dpro_loc(:).*omega(:);
            domega_dso4f_loc(:) = domega_dso4f_loc(:).*omega(:);
            
        case{'cc','arg','dlm'} % carbonates
            keq_tmp = keqsld_all(find(chrsld_all==mineral));
            ss_pco2 = stgas_all(find(chrsld_all==mineral),find(chrgas_all=='pco2'));
            omega(:) = 1d0;
            
            for ispa = 1:nsp_aq_all
                if (staq_all(find(chrsld_all==mineral),ispa) > 0d0) ;
                    omega(:) = omega(:).*maqf_loc(ispa,:)'.^staq_all(find(chrsld_all==mineral),ispa);
                    
                    % derivatives are first given as d(log omega)/dc 
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dmaq(ispa,:) ...
                        );
                    domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dpco2(ispa,:) ...
                        );
                    domega_dpro_loc(:) = domega_dpro_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dpro(ispa,:)' ...
                        );
                    domega_dso4f_loc(:) = domega_dso4f_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dso4f(ispa,:)' ...
                        );
                end 
            end 
            
            if (ss_pco2 > 0d0) 
                omega(:) = omega(:).*(k1*k2*kco2*pco2x(:)./(prox(:).^2d0)).^ss_pco2;
                
                % derivatives are first given as d(log omega)/dc 
                domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:) + ss_pco2./pco2x(:)'; 
                domega_dpro_loc(:) = domega_dpro_loc(:) - 2d0*ss_pco2./prox(:);
            end 
            
            if (keq_tmp > 0d0)  
                omega(:) = omega(:) / keq_tmp;
            end     
            
            % derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
            for ispa=1:nsp_aq_all
                domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:).*omega(:)';
            end 
            domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:).*omega(:)';
            domega_dpro_loc(:) = domega_dpro_loc(:).*omega(:);
            domega_dso4f_loc(:) = domega_dso4f_loc(:).*omega(:);
            
        case('gps') % sulfates
        % CaSO4*2H2O = Ca+2 + SO4-2 + 2H2O
            keq_tmp = keqsld_all(find(chrsld_all==mineral));
            omega(:) = 1d0;
            
            for ispa = 1:nsp_aq_all
                if (staq_all(find(chrsld_all,mineral),ispa) > 0d0)  
                    omega(:) = omega(:).*maqf_loc(ispa,:)'.^staq_all(find(chrsld_all==mineral),ispa);
                    
                    % derivatives are first given as d(log omega)/dc 
                    domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dmaq(ispa,:) ...
                        );
                    domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:).*dmaqf_dpco2(ispa,:) ...
                        ); 
                    domega_dpro_loc(:) = domega_dpro_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dpro(ispa,:)' ...
                        );
                    domega_dso4f_loc(:) = domega_dso4f_loc(:) + ( ...
                        + staq_all(find(chrsld_all==mineral),ispa)./maqf_loc(ispa,:)'.*dmaqf_dso4f(ispa,:)' ...
                        );
                end 
            end 
            
            if (keq_tmp > 0d0)  
                omega(:) = omega(:) / keq_tmp;
            end     
            
            % derivatives are now d(omega)/dc ( = d(omega)/d(log omega) * d(log omega)/dc = omega * d(log omega)/dc)
            for ispa=1:nsp_aq_all
                domega_dmaq_all(ispa,:) = domega_dmaq_all(ispa,:).*omega(:)';
            end 
            domega_dmgas_all(ipco2,:) = domega_dmgas_all(ipco2,:).*omega(:)';
            domega_dpro_loc(:) = domega_dpro_loc(:).*omega(:);
            domega_dso4f_loc(:) = domega_dso4f_loc(:).*omega(:);
            
        %%% other minerals that are assumed not to be controlled by distance from equilibrium i.e. omega
        
        case('py') % sulfides (assumed to be totally controlled by kinetics)
        % omega is defined so that kpy*poro*hr*mvpy*1d-6*mpyx*(1d0-omega_py) = kpy*poro*hr*mvpy*1d-6*mpyx*po2x^0.5d0
        % i.e., 1.0 - omega_py = po2x^0.5 
            % omega = 1d0 - po2x^0.5d0
            omega(:) = 1d0 - po2x(:).^0.5d0.*merge(0d0,1d0,po2x(:)<po2th*thon);
            domega_dmgas_all(ipo2,:) = - 0.5d0*po2x(:)'.^(-0.5d0).*merge(0d0,1d0,po2x(:)<po2th*thon)';
            
        case{'om','omb'}
            omega(:) = 1d0; % these are not used  
            
        case{'g1','g2','g3'}
        % omega is defined so that kg1*poro*hr*mvg1*1d-6*mg1x*(1d0-omega_g1) = kg1*poro*hr*mvg1*1d-6*mg1x*po2x/(po2x+mo2)
        % i.e., 1.0 - omega_g1 = po2x/(po2x+mo2) 
            if (mineral == 'g1'); mo2_tmp = mo2g1; end
            if (mineral == 'g2'); mo2_tmp = mo2g2; end
            if (mineral == 'g3'); mo2_tmp = mo2g3; end
            omega(:) = 1d0 - po2x(:)./(po2x(:)+mo2_tmp).*merge(0d0,1d0,po2x(:) < po2th*thon);
            domega_dmgas_all(ipo2,:) = ( ...
                - 1d0./(po2x(:)'+mo2_tmp).*merge(0d0,1d0,po2x(:) < po2th*thon)' ...
                - po2x(:)'*(-1d0)./(po2x(:)'+mo2_tmp).^2d0.*merge(0d0,1d0,po2x(:) < po2th*thon)' ...
                );
            
        otherwise
            % this should not be selected
            omega(:) = 1d0;
            error( '*** CAUTION: mineral (%s) saturation state is not defined --- > stop\n',mineral);
            
    end

    omega_error = false;
    if (any(isnan(omega)))  
        warning ('nan in calc_omega_v4 for %s ----> continue \n',mineral);
        omega_error = true;
    end 

end


function [ ... 
    rxn_ext,drxnext_dmsp,rxnext_error ...% output
    ] = calc_rxn_ext_dev_2( ...
    nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
    ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
    ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
    ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
    ,rg,tempk_0,tc ...%input
    ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
    ,rxn_name,sp_name ...% input 
    )

    % output variables initialization 
    rxnext_error = false;
    drxnext_dmsp =zeros(nz,1,'double');
    rxn_ext =zeros(nz,1,'double');

    % local variables initialization 
    po2th=0;fe2th=0;mwtom=0;g1th=0;g2th=0;g3th=0;mvpy=0;fe3th=0;knh3=0;k1nh3=0;ko2=0;v_tmp=0;km_tmp1=0;
    km_tmp2=0;km_tmp3=0;kn2o=0;k1fe2=0;k1fe2co3=0;k1fe2hco3=0;k1fe2so4=0;kco2=0;k1=0;k2=0;
    po2x=zeros(nz,1,'double');vmax=zeros(nz,1,'double');mo2=zeros(nz,1,'double');fe2x=zeros(nz,1,'double');
    koxa=zeros(nz,1,'double');vmax2=zeros(nz,1,'double');mom2=zeros(nz,1,'double');komb=zeros(nz,1,'double');
    beta=zeros(nz,1,'double');omx=zeros(nz,1,'double');ombx=zeros(nz,1,'double');mo2g1=zeros(nz,1,'double');
    mo2g2=zeros(nz,1,'double');mo2g3=zeros(nz,1,'double');kg1=zeros(nz,1,'double');kg2=zeros(nz,1,'double');
    kg3=zeros(nz,1,'double');g1x=zeros(nz,1,'double');g2x=zeros(nz,1,'double');g3x=zeros(nz,1,'double');
    pyx=zeros(nz,1,'double');fe3x=zeros(nz,1,'double');koxpy=zeros(nz,1,'double');pnh3x=zeros(nz,1,'double');
    nh4x=zeros(nz,1,'double');dnh4_dpro=zeros(nz,1,'double');dnh4_dpnh3=zeros(nz,1,'double');no3x=zeros(nz,1,'double');
    pn2ox=zeros(nz,1,'double');dv_dph_tmp=zeros(nz,1,'double');fe2f=zeros(nz,1,'double');dfe2f_dfe2=zeros(nz,1,'double');
    dfe2f_dpco2=zeros(nz,1,'double');dfe2f_dpro=zeros(nz,1,'double');dfe2f_dso4f=zeros(nz,1,'double');pco2x=zeros(nz,1,'double');
    hrpy=zeros(nz,1,'double');

    % thon = 1d0;
    thon = -1d100;

    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);
    [ieqaq_co3,ieqaq_hco3]=deal(1,2);
    [ieqaq_so4,ieqaq_so42]=deal(1,2);

    % --- start ----

    vmax(:) = krxn1_ext_all(find(chrrxn_ext_all=='resp'),:);
    mo2(:) = krxn2_ext_all(find(chrrxn_ext_all=='resp'),:);

    po2th = mgasth_all(find(chrgas_all=='po2'));

    koxa(:) = krxn1_ext_all(find(chrrxn_ext_all=='fe2o2'),:); 

    fe2th = maqth_all(find(chraq_all=='fe2'));

    fe3th = maqth_all(find(chraq_all=='fe3'));

    g1th = msldth_all(find(chrsld_all=='g1'));
    g2th = msldth_all(find(chrsld_all=='g2'));
    g3th = msldth_all(find(chrsld_all=='g3'));

    ko2 = keqgas_h(find(chrgas_all=='po2'),ieqgas_h0);

    kco2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0);
    k1 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);
    k2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);

    knh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h0);
    k1nh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h1);

    kn2o = keqgas_h(find(chrgas_all=='pn2o'),ieqgas_h0);

    k1fe2 = keqaq_h(find(chraq_all=='fe2'),ieqaq_h1);
    k1fe2co3 = keqaq_c(find(chraq_all=='fe2'),ieqaq_co3);
    k1fe2hco3  = keqaq_c(find(chraq_all=='fe2'),ieqaq_hco3);
    k1fe2so4 = keqaq_s(find(chraq_all=='fe2'),ieqaq_so4);

    vmax2(:) = krxn1_ext_all(find(chrrxn_ext_all=='omomb'),:);
    mom2(:) = krxn2_ext_all(find(chrrxn_ext_all=='omomb'),:);

    komb(:) = krxn1_ext_all(find(chrrxn_ext_all=='ombto'),:);
    beta(:) = krxn2_ext_all(find(chrrxn_ext_all=='ombto'),:);

    koxpy(:) = krxn1_ext_all(find(chrrxn_ext_all=='pyfe3'),:);

    po2x(:) = 0d0;
    if (any(chrgas=='po2'))  
        po2x(:) = mgasx(find(chrgas=='po2'),:);
    elseif (any(chrgas_cnst=='po2'))  
        po2x(:) = mgasc(find(chrgas_cnst=='po2'),:);
    end 

    pco2x(:) = 0d0;
    if (any(chrgas=='pco2'))  
        pco2x(:) = mgasx(find(chrgas=='pco2'),:);
    elseif (any(chrgas_cnst=='pco2'))  
        pco2x(:) = mgasc(find(chrgas_cnst=='pco2'),:);
    end 

    pnh3x(:) = 0d0;
    if (any(chrgas=='pnh3'))  
        pnh3x(:) = mgasx(find(chrgas=='pnh3'),:);
    elseif (any(chrgas_cnst=='pnh3'))  
        pnh3x(:) = mgasc(find(chrgas_cnst=='pnh3'),:);
    end 
    nh4x(:) = pnh3x(:)*knh3.*prox(:)/k1nh3;
    dnh4_dpro(:) = pnh3x(:)*knh3*1d0/k1nh3;
    dnh4_dpnh3(:) = 1d0*knh3*prox(:)/k1nh3;

    pn2ox(:) = 0d0;
    if (any(chrgas=='pn2o'))  
        pn2ox(:) = mgasx(find(chrgas=='pn2o'),:);
    elseif (any(chrgas_cnst=='pn2o'))  
        pn2ox(:) = mgasc(find(chrgas_cnst=='pn2o'),:);
    end 

    fe2x(:) = 0d0;
    if (any(chraq=='fe2'))  
        fe2x(:) = maqx(find(chraq=='fe2'),:);
    elseif (any(chraq_cnst=='fe2'))  
        fe2x(:) = maqc(find(chraq_cnst=='fe2'),:);
    end 
        
    fe2f(:) = fe2x(:)./(1d0+k1fe2./prox(:)+k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:)+k1fe2so4*so4f(:));
    dfe2f_dfe2(:) = 1d0./(1d0+k1fe2./prox(:)+k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:)+k1fe2so4*so4f(:));
    dfe2f_dpro(:) = fe2x(:)*(-1d0) ...
        ./(1d0+k1fe2./prox(:)+k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:)+k1fe2so4*so4f(:)).^2d0 ...
        *(k1fe2*(-1d0)./prox(:).^2d0+k1fe2co3*k1*k2*kco2*pco2x(:)*(-2d0)./prox(:).^3d0+k1fe2hco3*k1*k2*kco2*pco2x(:)*(-1d0)./prox(:).^2d0);
    dfe2f_dpco2(:) = fe2x(:)*(-1d0) ...
        /(1d0+k1fe2./prox(:)+k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:)+k1fe2so4*so4f(:)).^2d0 ...
        *(k1fe2co3*k1*k2*kco2*1d0./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*1d0./prox(:));
    dfe2f_dso4f(:) = fe2x(:)*(-1d0) ...
        ./(1d0+k1fe2./prox(:)+k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0+k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:)+k1fe2so4*so4f(:)).^2d0 ...
        * k1fe2so4;

    fe3x(:) = 0d0;
    if (any(chraq=='fe3'))  
        fe3x(:) = maqx(find(chraq=='fe3'),:);
    elseif (any(chraq_cnst=='fe3'))  
        fe3x(:) = maqc(find(chraq_cnst=='fe3'),:);
    end 

    no3x(:) = 0d0;
    if (any(chraq=='no3'))  
        no3x(:) = maqx(find(chraq=='no3'),:);
    elseif (any(chraq_cnst=='no3'))  
        no3x(:) = maqc(find(chraq_cnst=='no3'),:);
    end 

    omx(:) = 0d0;
    if (any(chrsld=='om'))  
        omx(:) = msldx(find(chrsld=='om'),:);
    elseif (any(chraq_cnst=='om'))  
        omx(:) = msldc(find(chrsld_cnst=='om'),:);
    end 

    ombx(:) = 0d0;
    if (any(chrsld=='omb'))  
        ombx(:) = msldx(find(chrsld=='omb'),:);
    elseif (any(chraq_cnst=='omb'))  
        ombx(:) = msldc(find(chrsld_cnst=='omb'),:);
    end 

    g1x(:) = 0d0;
    if (any(chrsld=='g1'))  
        g1x(:) = msldx(find(chrsld=='g1'),:);
    elseif (any(chraq_cnst=='g1'))  
        g1x(:) = msldc(find(chrsld_cnst=='g1'),:);
    end 

    g2x(:) = 0d0;
    if (any(chrsld=='g2'))  
        g2x(:) = msldx(find(chrsld=='g2'),:);
    elseif (any(chraq_cnst=='g2'))  
        g2x(:) = msldc(find(chrsld_cnst=='g2'),:);
    end 

    g3x(:) = 0d0;
    if (any(chrsld=='g3'))  
        g3x(:) = msldx(find(chrsld=='g3'),:);
    elseif (any(chraq_cnst=='g3'))  
        g3x(:) = msldc(find(chrsld_cnst=='g3'),:);
    end 

    pyx(:) = 0d0;
    if (any(chrsld=='py'))  
        pyx(:) = msldx(find(chrsld=='py'),:);
    elseif (any(chraq_cnst=='py'))  
        pyx(:) = msldc(find(chrsld_cnst=='py'),:);
    end 

    mvpy = mv_all(find(chrsld_all=='py'));

    hrpy(:) = 0d0;
    if (any(chrsld=='py'))  
        hrpy(:) = hr(find(chrsld=='py'),:);
    end 

    switch rxn_name

        case('resp')
            rxn_ext(:) = vmax(:).*po2x(:)./(po2x(:)+mo2(:));
            
            switch sp_name
                case('po2')
                    drxnext_dmsp(:) = (...
                        vmax(:)*1d0./(po2x(:)+mo2(:)) ...
                        +vmax(:).*po2x(:)*(-1d0)./(po2x(:)+mo2(:)).^2d0 ...
                        );
                otherwise
                    drxnext_dmsp(:) = 0d0;
            end 
            
        case('fe2o2')
            % scheme = 'full'; % reflecting individual rate consts for different Fe2+ species (after Kanzaki and Murakami 2016)
            scheme = 'default'; % as a function of pH and pO2
            
            switch scheme
                case('full')
                    rxn_ext(:) = ( ...
                        + poro(:).*sat(:)*1d3.*fe2f(:).*( ...
                        + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                        + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                        + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)/prox(:)^2d0 ...
                        + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                        ).*po2x(:) ...
                        );
                    
                    switch sp_name
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                + poro(:).*sat(:)*1d3.*dfe2f_dpro(:).*( ...
                                + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                                ).*po2x(:) ...
                                + poro(:).*sat(:)*1d3.*fe2f(:).*( ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2*(-1d0)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg) ...
                                      *k1fe2co3*k1*k2*kco2*pco2x(:)*(-2d0)./prox(:).^3d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg) ...
                                      *k1fe2hco3*k1*k2*kco2*pco2x(:)*(-1d0)./prox(:).^2d0 ...
                                ).*po2x(:) ...
                                );
                        case('so4f')
                            drxnext_dmsp(:) = ( ...
                                + poro(:).*sat(:)*1d3.*dfe2f_dso4f(:).*( ...
                                + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                                ).*po2x(:) ...
                                );
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                + poro(:).*sat(:)*1d3.*fe2f(:).*( ...
                                + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                                )*1d0 ...
                                );
                        case('pco2')
                            drxnext_dmsp(:) = ( ...
                                + poro(:).*sat(:)*1d3.*dfe2f_dpco2(:).*( ...
                                + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                                ).*po2x(:) ...
                                + poro(:).*sat(:)*1d3.*fe2f(:).*( ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*1d0./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*1d0./prox(:) ...
                                ).*po2x(:) ...
                                );
                        case('fe2')
                            drxnext_dmsp(:) = ( ...
                                + poro(:).*sat(:)*1d3.*dfe2f_dfe2(:).*( ...
                                + k_arrhenius(10d0^(1.46d0),25d0+tempk_0,tc+tempk_0,46d0,rg) ...
                                + k_arrhenius(10d0^(8.34d0),25d0+tempk_0,tc+tempk_0,21.6d0,rg)*k1fe2./prox(:) ...
                                + k_arrhenius(10d0^(6.27d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2co3*k1*k2*kco2*pco2x(:)./prox(:).^2d0 ...
                                + k_arrhenius(10d0^(5.12d0),25d0+tempk_0,tc+tempk_0,29d0,rg)*k1fe2hco3*k1*k2*kco2*pco2x(:)./prox(:) ...
                                ).*po2x(:) ...
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end
                    
                otherwise
                    rxn_ext(:) = ( ...
                        poro(:).*sat(:)*1d3.*fe2x(:).*po2x(:) ...
                        .*(8.0d13*60.0d0*24.0d0*365.0d0*(kw./prox(:)).^2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) ...
                        .*merge(0d0,1d0,po2x(:) < po2th*thon | fe2x(:) < fe2th*thon) ...
                        );
                    
                    switch sp_name
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                poro(:).*sat(:)*1d3.*fe2x(:).*po2x(:) ...
                                .*(8.0d13*60.0d0*24.0d0*365.0d0*(kw./prox(:))*2.0d0*(kw*(-1d0)./prox(:).^2d0)) ...
                                .*merge(0d0,1d0,po2x(:) < po2th*thon | fe2x(:) < fe2th*thon) ...
                                );
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                poro(:).*sat(:)*1d3.*fe2x(:)*1d0 ...
                                .*(8.0d13*60.0d0*24.0d0*365.0d0*(kw./prox(:)).^2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) ...
                                .*merge(0d0,1d0,po2x(:) < po2th*thon | fe2x(:) < fe2th*thon) ...
                                );
                        case('fe2')
                            drxnext_dmsp(:) = ( ...
                                poro(:).*sat(:)*1d3*1d0.*po2x(:) ...
                                .*(8.0d13*60.0d0*24.0d0*365.0d0*(kw./prox(:)).^2.0d0 + 1d-7*60.0d0*24.0d0*365.0d0) ...
                                .*merge(0d0,1d0,po2x(:) < po2th*thon | fe2x(:) < fe2th*thon) ...
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end
            end 
        
        case('omomb')
            rxn_ext(:) = vmax2(:) ... % mg C / soil g /yr
                .*omx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...% mol/m3 converted to mg C/ soil g
                .*ombx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...
                ./(mom2(:) + (omx(:).*(1d0-poro(:))*rho_grain)*1d6*12d0*1d3) ...
                *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6); % converting mg_C/soil_g to mol_C/soil_m3
            
            switch sp_name
                case('om')
                    drxnext_dmsp(:) = ( ...
                        vmax2(:) ... % mg C / soil g /yr
                        *1d0.*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...% mol/m3 converted to mg C/ soil g
                        .*ombx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...
                        ./(mom2(:) + (omx(:).*(1d0-poro(:))*rho_grain)*1d6*12d0*1d3) ...
                        *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6) ...% converting mg_C/soil_g to mol_C/soil_m3
                        + vmax2(:) ... % mg C / soil g /yr
                        .*omx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...% mol/m3 converted to mg C/ soil g
                        .*ombx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...
                        *(-1d0)./(mom2(:) + (omx(:).*(1d0-poro(:))*rho_grain)*1d6*12d0*1d3).^2d0 ...
                        .*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...
                        *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6) ...% converting mg_C/soil_g to mol_C/soil_m3
                        ); 
                case('omb')
                    drxnext_dmsp(:) = ( ...
                        vmax2(:) ... % mg C / soil g /yr
                        .*omx(:).*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...% mol/m3 converted to mg C/ soil g
                        *1d0.*(1d0-poro(:))*rho_grain*1d6*12d0*1d3 ...
                        ./(mom2(:) + (omx(:).*(1d0-poro(:))*rho_grain)*1d6*12d0*1d3) ...
                        *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6) ...% converting mg_C/soil_g to mol_C/soil_m3
                        ); 
                otherwise
                    drxnext_dmsp(:) = 0d0;
            end
        
        case('ombto')
            rxn_ext(:) = komb(:).*(ombx(:).*(1d0-poro(:))*rho_grain*rho_grain*1d6).^beta(:) ...
                *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6); % converting mg_C/soil_g to mol_C/soil_m3
            
            switch sp_name
                case('omb')
                    drxnext_dmsp(:) = ( ...
                        komb(:).*beta(:).*(ombx(:).*(1d0-poro(:))*rho_grain*rho_grain*1d6).^(beta(:)-1d0) ...
                        .*(1d0-poro(:))*rho_grain*rho_grain*1d6 ...
                        *1d-3/12d0./((1d0-poro(:))*rho_grain*1d6) ...% converting mg_C/soil_g to mol_C/soil_m3
                        ); 
        
                otherwise
                    drxnext_dmsp(:) = 0d0;
                    
            end 
                    
        
        case('pyfe3') 
            rxn_ext(:) = ( ...
                koxpy(:).*poro(:).*hrpy(:)*mvpy.*pyx(:).*fe3x(:).^0.93d0.*fe2x(:).^(-0.40d0) ...
                .*merge(0d0,1d0,fe3x(:)<fe3th*thon | fe2x(:)<fe2th*thon) ...
                ./(1d0 - poro(:)) ...
                );
            
            switch sp_name
                case('py')
                    drxnext_dmsp(:) = (...
                        koxpy(:).*poro(:).*hrpy(:)*mvpy*1d0.*fe3x(:).^0.93d0.*fe2x(:).^(-0.40d0) ...
                        .*merge(0d0,1d0,fe3x(:)<fe3th*thon | fe2x(:)<fe2th*thon) ...
                        );
                case('fe3')
                    drxnext_dmsp(:) = (...
                        koxpy(:).*poro(:).*hrpy(:)*mvpy.*pyx(:)*(0.93d0).*fe3x(:).^(0.93d0-1d0).*fe2x(:).^(-0.40d0) ...
                        .*merge(0d0,1d0,fe3x(:)<fe3th*thon | fe2x(:)<fe2th*thon) ...
                        )
                case('fe2')
                    drxnext_dmsp(:) = (...
                        koxpy(:).*poro(:).*hrpy(:)*mvpy.*pyx(:).*fe3x(:).^0.93d0*(-0.4d0).*fe2x(:).^(-0.40d0-1d0) ...
                        .*merge(0d0,1d0,fe3x(:)<fe3th*thon | fe2x(:)<fe2th*thon) ...
                        )
                otherwise
                    drxnext_dmsp(:) = 0d0;
            end 
            
        case('amo2o')
            scheme = 'maggi08'; % Maggi et al. (2008) wihtout baterial, pH and water saturation functions
            % scheme = 'Fennel'; % from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 with a correction 
            % scheme = 'FennelOLD'; % from biogem_box_geochem.f90 in GENIE model referring to Fennel et al. 2005 without a correction 
            % scheme = 'Ozaki'; % from biogem_box_geochem.f90 in GENIE model referring to Ozaki et al. [EPSL ... ?]
            
            switch scheme
                case('maggi08')
                    v_tmp = 9.53d-6*60d0*60d0*24d0*365d0; % (~300 /yr)
                    % v_tmp = v_tmp/100d0; % (~3 /yr; default value produces too much nitrate (pH goes down to ~1)
                    km_tmp1 = 14d-5;
                    km_tmp2 = 2.41d-5;
                    rxn_ext(:) = ( ...
                        v_tmp ...
                        *nh4x(:)./(nh4x(:) + km_tmp1 ) ...
                        .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2 ) ...
                        .*min(2d0*sat(:),1d0) ...
                        .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                        );
                    
                    dv_dph_tmp(:) = 0d0;
                    for iz = 1:nz
                        if (3d0 < -log10(prox(iz)) && -log10(prox(iz)) < 7d0)
                            dv_dph_tmp(iz) = 0.25d0;
                        elseif (7d0 < -log10(prox(iz)) && -log10(prox(iz)) < 11d0)
                            dv_dph_tmp(iz) = -0.25d0;
                        elseif (-log10(prox(iz)) == 7d0)
                            dv_dph_tmp(iz) = 0d0;
                        elseif (-log10(prox(iz)) == 3d0)
                            dv_dph_tmp(iz) = 0.125d0;
                        elseif (-log10(prox(iz)) == 11d0)
                            dv_dph_tmp(iz) = -0.125d0;
                        else
                            dv_dph_tmp(iz) = 0d0;
                        end
                    end
                    
                    % when using modified version using normal distribution with sigma = 1 
                    dv_dph_tmp(:) = exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...
                        *-0.5d0*2d0.*((-log10(prox(:))-7d0)/1d0)  ...
                        *(-1d0) ...
                        *1d0/log(10d0)./prox(:);
                    
                    switch sp_name
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *nh4x(:)./(nh4x(:) + km_tmp1 ) ...
                                .*( ...
                                1d0*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                + po2x(:)*ko2*(-1d0)./(po2x(:)*ko2 + km_tmp2).^2d0 * ko2 ...
                                ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('pnh3')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * ( ... 
                                dnh4_dpnh3(:)./(nh4x(:) + km_tmp1 ) ...
                                + nh4x(:)*(-1d0)./(nh4x(:) + km_tmp1 ).^2d0 .* dnh4_dpnh3(:) ...
                                ) ...
                                .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * ( ... 
                                dnh4_dpro(:)./(nh4x(:) + km_tmp1 ) ...
                                + nh4x(:)*(-1d0)./(nh4x(:) + km_tmp1 ).^2d0 .* dnh4_dpro(:) ...
                                ) ...
                                .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                + ...
                                v_tmp ...
                                *nh4x(:)./(nh4x(:) + km_tmp1 ) ...
                                .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*dv_dph_tmp(:) ...
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end 
                    
                case{'Fennel','FennelOLD'}
                    if (scheme== 'Fennel'); v_tmp = 6.0d0; end % /yr
                    if (scheme== 'FennelOLD'); v_tmp = 0.16667d0; end % /yr
                    km_tmp2 = 2.0D-05;
                    rxn_ext(:) = ( ...
                        v_tmp ...
                        *nh4x(:) ...
                        .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2 ) ...
                        );
                    
                    switch sp_name
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *nh4x(:) ...
                                .*( ...
                                1d0*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                + po2x(:)*ko2*(-1d0)./(po2x(:)*ko2 + km_tmp2).^2d0 * ko2 ...
                                ) ...
                                );
                        case('pnh3')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *dnh4_dpnh3(:) ...
                                .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                );
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * dnh4_dpro(:) ...
                                .*po2x(:)*ko2./(po2x(:)*ko2 + km_tmp2) ...
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end 
                    
                case('Ozaki')
                    v_tmp = 18250.0d0/1027.649d0; % /yr
                    rxn_ext(:) = ( ...
                        v_tmp ...
                        *nh4x(:) ...
                        .*po2x(:)*ko2 ...
                        );
                    
                    switch sp_name
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *nh4x(:) ...
                                * 1d0*ko2 ...
                                );
                        case('pnh3')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *dnh4_dpnh3(:) ...
                                .*po2x(:)*ko2 ...
                                );
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * dnh4_dpro(:) ...
                                .*po2x(:)*ko2 ...
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end 
                    
                end
                    
        case{'g2n0','g2n21'} 
            % overall denitrification (4 NO3-  +  5 CH2O  +  4 H+  ->  2 N2  +  5 CO2  +  7 H2O) 
            % first of 2 step denitrification (2 NO3-  +  2 CH2O  +  2 H+  ->  N2O  +  2 CO2  +  3 H2O)   
            % (assuming that oxidation by N2O governs overall denitrification)
            scheme = 'maggi08'; % Maggi et al. (2008) wihtout baterial, pH and water saturation functions; vmax from oxidation by N2O (rate-limiting)
            
            switch scheme
                case('maggi08')
                    v_tmp = 1.23d-7*60d0*60d0*24d0*365d0;
                    km_tmp1 = 10d-5 * 1d6; % mol L-1 converted to mol m-3
                    km_tmp2 = 11.3d-5;
                    km_tmp3 = 2.52d-5;
                    rxn_ext(:) = ( ...
                        v_tmp ...
                        *g2x(:)./(g2x(:) + km_tmp1 ) ...
                        .*no3x(:)./(no3x(:) + km_tmp2 ) ...
                        *km_tmp3./(po2x(:)*ko2 + km_tmp3 ) ...
                        .*min(2d0*sat(:),1d0) ...
                        .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                        );
                    
                    % when using modified version using normal distribution with sigma = 1 
                    dv_dph_tmp(:) = exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...
                        *-0.5d0*2d0.*((-log10(prox(:))-7d0)/1d0)  ...
                        *(-1d0) ...
                        *1d0/log(10d0)./prox(:);
                        
                    switch sp_name
                        case('g2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * ( ... 
                                1d0./(g2x(:) + km_tmp1 ) ...
                                + g2x(:)*(-1d0)./(g2x(:) + km_tmp1 ).^2d0 * 1d0 ...
                                ) ...
                                .*no3x(:)./(no3x(:) + km_tmp2 ) ...
                                *km_tmp3./(po2x(:)*ko2 + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('no3')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                .* ( ... 
                                1d0./(no3x(:) + km_tmp2 ) ...
                                + no3x(:)*(-1d0)./(no3x(:) + km_tmp2 ).^2d0 * 1d0 ...
                                ) ...
                                *km_tmp3./(po2x(:)*ko2 + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('po2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                .*no3x(:)./(no3x(:) + km_tmp2 ) ...
                                *km_tmp3*(-1d0)./(po2x(:)*ko2 + km_tmp3 ).^2d0 * ko2 ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                .*no3x(:)./(no3x(:) + km_tmp2 ) ...
                                *km_tmp3./(po2x(:)*ko2 + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*dv_dph_tmp(:) ...% modified version using normal distribution with sigma = 1 
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end 
                    
            end 
                    
        case('g2n22') % 2nd of 2 step denitrification (2 N2O  +  CH2O  ->  2 N2  +  CO2  +  H2O)  
            scheme = 'maggi08'; % Maggi et al. (2008) wihtout baterial, pH and water saturation functions
            
            switch scheme
                case('maggi08')
                    v_tmp = 1.23d-7*60d0*60d0*24d0*365d0;
                    km_tmp1 = 10d-5 * 1d6; % mol L-1 converted to mol m-3
                    km_tmp2 = 11.3d-5;
                    km_tmp3 = 2.52d-5;
                    rxn_ext(:) = ( ...
                        v_tmp ...
                        *g2x(:)./(g2x(:) + km_tmp1 ) ...
                        *kn2o.*pn2ox(:)./(kn2o*pn2ox(:) + km_tmp2 ) ...
                        *km_tmp3./(no3x(:) + km_tmp3 ) ...
                        .*min(2d0*sat(:),1d0) ...
                        .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                        );
                    
                    % when using modified version using normal distribution with sigma = 1 
                    dv_dph_tmp(:) = exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...
                        *-0.5d0*2d0.*((-log10(prox(:))-7d0)/1d0)  ...
                        *(-1d0) ...
                        *1d0/log(10d0)./prox(:);
                        
                    switch sp_name
                        case('g2')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                * ( ... 
                                1d0./(g2x(:) + km_tmp1 ) ...
                                + g2x(:)*(-1d0)./(g2x(:) + km_tmp1 ).^2d0 * 1d0 ...
                                ) ...
                                *kn2o.*pn2ox(:)./(kn2o*pn2ox(:) + km_tmp2 ) ...
                                *km_tmp3./(no3x(:) + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('pn2o')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                .* ( ... 
                                kn2o./(kn2o*pn2ox(:) + km_tmp2 ) ...
                                + kn2o*pn2ox(:)*(-1d0)./(kn2o*pn2ox(:) + km_tmp2 ).^2d0 * kn2o ...
                                ) ...
                                *km_tmp3./(no3x(:) + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('no3')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                *kn2o.*pn2ox(:)./(kn2o*pn2ox(:) + km_tmp2 ) ...
                                *km_tmp3*(-1d0)./(no3x(:) + km_tmp3 ).^2d0 * 1d0 ...
                                .*min(2d0*sat(:),1d0) ...
                                .*exp(-0.5d0*((-log10(prox(:))-7d0)/1d0).^2d0) ...% modified version using normal distribution with sigma = 1 
                                );
                        case('pro')
                            drxnext_dmsp(:) = ( ...
                                v_tmp ...
                                *g2x(:)./(g2x(:) + km_tmp1 ) ...
                                *kn2o.*pn2ox(:)./(kn2o*pn2ox(:) + km_tmp2 ) ...
                                *km_tmp3./(no3x(:) + km_tmp3 ) ...
                                .*min(2d0*sat(:),1d0) ...
                                .*dv_dph_tmp(:) ...% modified version using normal distribution with sigma = 1 
                                );
                        otherwise
                            drxnext_dmsp(:) = 0d0;
                    end 
                    
            end 
            
            
        otherwise 
            rxn_ext(:) = 0d0;
            drxnext_dmsp(:) = 0d0;
            
    end

    rxnext_error = false;
    if (any(isnan(rxn_ext)) || any(isnan(drxnext_dmsp)))  
        warning('nan in calc_rxn_ext_dev_2');
        rxnext_error = true;
    end 

end 


function [ ... 
    khgas,khgasx,dkhgas_dpro,dkhgas_dso4f,dkhgas_dmaq,dkhgas_dmgas ...%output
    ] = calc_khgas_all( ...
    nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst ...
    ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas ...
    ,maq,mgas,maqx,mgasx,maqc,mgasc ...
    ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3  ...
    ,pro,prox,so4fprev,so4f ...
    )
    
    % output 
    khgas=zeros(nsp_gas_all,nz,'double');khgasx=zeros(nsp_gas_all,nz,'double');
    dkhgas_dpro=zeros(nsp_gas_all,nz,'double');dkhgas_dso4f=zeros(nsp_gas_all,nz,'double');
    dkhgas_dmgas=zeros(nsp_gas_all,nsp_gas_all,nz,'double');
    dkhgas_dmaq=zeros(nsp_gas_all,nsp_aq_all,nz,'double');

    % local 
    maqx_loc=zeros(nsp_aq_all,nz,'double');maq_loc=zeros(nsp_aq_all,nz,'double');
    maqf_loc=zeros(nsp_aq_all,nz,'double');maqf_loc_prev=zeros(nsp_aq_all,nz,'double');
    mgasx_loc=zeros(nsp_gas_all,nz,'double');mgas_loc=zeros(nsp_gas_all,nz,'double');
    dmaqf_dpro=zeros(nsp_aq_all,nz,'double');dmaqf_dso4f=zeros(nsp_aq_all,nz,'double');
    dmaqf_dmaq=zeros(nsp_aq_all,nz,'double');dmaqf_dpco2=zeros(nsp_aq_all,nz,'double');

    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    ispg=0;ispa=0;ispa_c=0;ipco2=0;ipnh3=0;io2=0;in2o=0;
    kco2=0;k1=0;k2=0;knh3=0;k1nh3=0;kho=0;kn2o=0;

    % --- start ----


    ipco2 = find(chrgas_all=='pco2');
    ipnh3 = find(chrgas_all=='pnh3');
    io2 = find(chrgas_all=='po2');
    in2o = find(chrgas_all=='pn2o');

    kco2 = keqgas_h(ipco2,ieqgas_h0);
    k1 = keqgas_h(ipco2,ieqgas_h1);
    k2 = keqgas_h(ipco2,ieqgas_h2);

    knh3 = keqgas_h(ipnh3,ieqgas_h0);
    k1nh3 = keqgas_h(ipnh3,ieqgas_h1);

    kho = keqgas_h(io2,ieqgas_h0);

    kn2o = keqgas_h(in2o,ieqgas_h0);

    khgas(:,:) = 0d0;
    khgasx(:,:) = 0d0;

    dkhgas_dpro(:,:) = 0d0;
    dkhgas_dso4f(:,:) = 0d0;
    dkhgas_dmgas(:,:,:) = 0d0;
    dkhgas_dmaq(:,:,:) = 0d0;

    for ispg = 1: nsp_gas_all
        switch (chrgas_all(ispg))
            case('pco2')
                khgas(ispg,:) = kco2*(1d0+k1./pro(:)' + k1*k2./pro(:)'./pro(:)'); % previous value; should not change through iterations 
                khgasx(ispg,:) = kco2*(1d0+k1./prox(:)' + k1*k2./prox(:)'./prox(:)');
                
                dkhgas_dpro(ispg,:) = kco2*(k1*(-1d0)./prox(:)'.^2d0 + k1*k2*(-2d0)./prox(:)'.^3d0);
                
                % obtain previous data 
                [maq_loc,mgas_loc] = get_maqgasx_all( ...
                    nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst ...
                    ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst ...
                    ,maq,mgas,maqc,mgasc ...
                    );
                [ ...
                    dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
                    ,maqf_loc_prev  ...% output
                    ] = get_maqf_all( ...
                    nz,nsp_aq_all,nsp_gas_all ...
                    ,chraq_all,chrgas_all ...
                    ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
                    ,mgas_loc,maq_loc,pro,so4fprev ...
                    );
                    
                [maqx_loc,mgasx_loc] = get_maqgasx_all( ...
                    nz,nsp_aq_all,nsp_gas_all,nsp_aq,nsp_gas,nsp_aq_cnst,nsp_gas_cnst ...
                    ,chraq,chraq_all,chraq_cnst,chrgas,chrgas_all,chrgas_cnst ...
                    ,maqx,mgasx,maqc,mgasc ...
                    );
                % getting free maq
                [ ... 
                    dmaqf_dpro,dmaqf_dso4f,dmaqf_dmaq,dmaqf_dpco2 ...% output
                    ,maqf_loc  ...% output
                    ] = get_maqf_all( ...
                    nz,nsp_aq_all,nsp_gas_all ...
                    ,chraq_all,chrgas_all ...
                    ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
                    ,mgasx_loc,maqx_loc,prox,so4f ...
                    );
                    
                % account for species associated with CO3-- (ispa_c =1) and HCO3- (ispa_c =2)
                for ispa = 1: nsp_aq_all
                    for ispa_c = 1:2
                        if ( keqaq_c(ispa,ispa_c) > 0d0)  
                            if (ispa_c == 1)  % with CO3--
                                khgas(ispg,:) = khgas(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc_prev(ispa,:)*k1*k2*kco2.*pro(:)'.^(-2d0) ...
                                    );
                                khgasx(ispg,:) = khgasx(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2.*prox(:)'.^(-2d0) ...
                                    );
                                dkhgas_dpro(ispg,:) = dkhgas_dpro(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*(-2d0).*prox(:)'.^(-3d0) ...
                                    + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-2d0) ...
                                    .*dmaqf_dpro(ispa,:) ...
                                    );
                                dkhgas_dso4f(ispg,:) = dkhgas_dso4f(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-2d0) ...
                                    .*dmaqf_dso4f(ispa,:) ...
                                    );
                                % dkhgas_dmgas(ispg,ipco2,:) = dkhgas_dmgas(ispg,ipco2,:) + ( ...
                                    % + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-2d0) ...
                                    % .*dmaqf_dpco2(ispa,:) ...
                                    % );
                                % dkhgas_dmaq(ispg,ispa,:) = dkhgas_dmaq(ispg,ispa,:) + ( ...
                                    % + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-2d0) ...
                                    % .*dmaqf_dmaq(ispa,:) ...
                                    % );
                                % MATLAB sucks 
                                for iz=1:nz
                                    dkhgas_dmgas(ispg,ipco2,iz) = dkhgas_dmgas(ispg,ipco2,iz) + ( ...
                                        + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(iz)^(-2d0) ...
                                        *dmaqf_dpco2(ispa,iz) ...
                                        );
                                    dkhgas_dmaq(ispg,ispa,iz) = dkhgas_dmaq(ispg,ispa,iz) + ( ...
                                        + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(iz)^(-2d0) ...
                                        *dmaqf_dmaq(ispa,iz) ...
                                        );
                                end 
                            elseif (ispa_c == 2)  % with HCO3-
                                khgas(ispg,:) = khgas(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc_prev(ispa,:)*k1*k2*kco2.*pro(:)'.^(-1d0) ... 
                                    );
                                khgasx(ispg,:) = khgasx(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2.*prox(:)'.^(-1d0) ... 
                                    );
                                dkhgas_dpro(ispg,:) = dkhgas_dpro(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*maqf_loc(ispa,:)*k1*k2*kco2*(-1d0).*prox(:)'.^(-2d0) ... 
                                    + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-1d0) ... 
                                    .*dmaqf_dpro(ispa,:) ...
                                    );
                                dkhgas_dso4f(ispg,:) = dkhgas_dso4f(ispg,:) + ( ...
                                    + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-1d0) ... 
                                    .*dmaqf_dso4f(ispa,:) ...
                                    );
                                % dkhgas_dmgas(ispg,ipco2,:) = dkhgas_dmgas(ispg,ipco2,:) + ( ...
                                    % + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-1d0) ... 
                                    % .*dmaqf_dpco2(ispa,:) ...
                                    % );
                                % dkhgas_dmaq(ispg,ispa,:) = dkhgas_dmaq(ispg,ispa,:) + ( ...
                                    % + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(:)'.^(-1d0) ... 
                                    % .*dmaqf_dmaq(ispa,:) ...
                                    % );
                                % MATLAB sucks
                                for iz=1:nz
                                    dkhgas_dmgas(ispg,ipco2,iz) = dkhgas_dmgas(ispg,ipco2,iz) + ( ...
                                        + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(iz)^(-1d0) ... 
                                        *dmaqf_dpco2(ispa,iz) ...
                                        );
                                    dkhgas_dmaq(ispg,ispa,iz) = dkhgas_dmaq(ispg,ispa,iz) + ( ...
                                        + keqaq_c(ispa,ispa_c)*k1*k2*kco2*prox(iz)^(-1d0) ... 
                                        *dmaqf_dmaq(ispa,iz) ...
                                        );
                                end 
                                
                            end 
                        end 
                    end 
                end 
                
            case('po2')
                khgas(ispg,:) = kho; % previous value; should not change through iterations 
                khgasx(ispg,:) = kho;

            case('pnh3')
                khgas(ispg,:) = knh3*(1d0+pro(:)/k1nh3); % previous value; should not change through iterations 
                khgasx(ispg,:) = knh3*(1d0+prox(:)/k1nh3);

                dkhgas_dpro(ispg,:) = knh3*(1d0/k1nh3);
            case('pn2o')
                khgas(ispg,:) = kn2o; % previous value; should not change through iterations 
                khgasx(ispg,:) = kn2o;
        end 

    end 
    
end


function [ ... 
    rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas ...% output
    ] = sld_rxn( ...
    nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz ...% input 
    ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle ...% input
    ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas ...% input
    ) 
    
    % output 
    rxnsld=zeros(nsp_sld,nz,'double');drxnsld_dmsld=zeros(nsp_sld,nz,'double');
    drxnsld_dmaq=zeros(nsp_sld,nsp_aq,nz,'double');
    drxnsld_dmgas=zeros(nsp_sld,nsp_gas,nz,'double');

    % local
    ispa=0;isps=0;ispg=0;iz=0;
    maxdis=zeros(nsp_sld,nz,'double');maxprec=zeros(nsp_sld,nz,'double');

    auth_th = 1d2;

    % -- start --
        
    rxnsld(:,:) = 0d0;
    drxnsld_dmsld(:,:) = 0d0;
    drxnsld_dmaq(:,:,:) = 0d0;
    drxnsld_dmgas(:,:,:) = 0d0;

    maxdis(:,:) = 1d200;
    maxprec(:,:) = -1d200;

    for isps=1:nsp_sld
        maxdis(isps,:) = min(maxdis(isps,:),(msld(isps,:)-msldth(isps))/dt);
        for ispa = 1:nsp_aq
            if (staq(isps,ispa)<0d0)  
                maxdis(isps,:) = min(maxdis(isps,:), -1d0/staq(isps,ispa)*poro(:)'.*sat(:)'*1d3.*(maq(ispa,:)-maqth(ispa))/dt );
            end 
        end 
        for ispg = 1:nsp_gas
            if (stgas(isps,ispg)<0d0)  
                maxdis(isps,:) = min(maxdis(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:).*(mgas(ispg,:)-mgasth(ispg))/dt );
            end 
        end 
        for ispa = 1:nsp_aq
            if (staq(isps,ispa)>0d0)  
                maxprec(isps,:) = max(maxprec(isps,:), -1d0/staq(isps,ispa)*poro(:)'.*sat(:)'*1d3.*(maq(ispa,:)-maqth(ispa))/dt );
            end 
        end 
        for ispg = 1:nsp_gas
            if (stgas(isps,ispg)>0d0)  
                maxprec(isps,:) = max(maxprec(isps,:), -1d0/stgas(isps,ispg)*agas(ispg,:).*(mgas(ispg,:)-mgasth(ispg))/dt );
            end 
        end 
    end 

    for isps = 1:nsp_sld
        switch (precstyle(isps))
        
            case ('full_lim') 
                
                for iz = 1:nz
                    if (1d0-omega(isps,iz) > 0d0)  
                        rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) ;
                        if (rxnsld(isps,iz)> maxdis(isps,iz))  
                            rxnsld(isps,iz) = maxdis(isps,iz);
                            drxnsld_dmsld(isps,iz) = 0d0;
                            drxnsld_dmaq(isps,:,iz) = 0d0;
                            drxnsld_dmgas(isps,:,iz) = 0d0;
                        else 
                            drxnsld_dmsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz)) ;
                            drxnsld_dmaq(isps,:,iz) = ( ...
                                +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmaq(isps,:,iz)) ...
                                +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                                );
                            drxnsld_dmgas(isps,:,iz) = ( ...
                                +ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(-domega_dmgas(isps,:,iz)) ...
                                +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                                );
                        end 
                    elseif (1d0-omega(isps,iz) < 0d0)  
                        if (nonprec(isps,iz)==1d0)  
                            rxnsld(isps,iz) = 0d0;
                            drxnsld_dmsld(isps,iz) = 0d0;
                            drxnsld_dmaq(isps,:,iz) = 0d0;
                            drxnsld_dmgas(isps,:,iz) = 0d0;
                        elseif (nonprec(isps,iz)==0d0)  
                            rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz));
                            if (rxnsld(isps,iz) < maxprec(isps,iz))  
                                rxnsld(isps,iz) = maxprec(isps,iz);
                                drxnsld_dmsld(isps,iz) = 0d0;
                                drxnsld_dmaq(isps,:,iz) = 0d0;
                                drxnsld_dmgas(isps,:,iz) = 0d0;
                            else
                                rxnsld(isps,iz) = ksld(isps,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz));
                                drxnsld_dmsld(isps,iz) = 0d0;
                                drxnsld_dmaq(isps,:,iz) = ( ...
                                    +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmaq(isps,:,iz)) ...
                                    +dksld_dmaq(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) ...
                                    );
                                drxnsld_dmgas(isps,:,iz) = ( ...
                                    +ksld(isps,iz)*poro(iz)*hr(isps,iz)*(-domega_dmgas(isps,:,iz)) ...
                                    +dksld_dmgas(isps,:,iz)*poro(iz)*hr(isps,iz)*(1d0-omega(isps,iz)) ...
                                    );
                            end 
                        end 
                    end 
                end 
                
        
            case ('full') 
                rxnsld(isps,:) = ( ...
                    + min(ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*msldx(isps,:).*(1d0-omega(isps,:))./(1d0-poro(:)),maxdis(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                     + max(ksld(isps,:).*poro(:)'.*hr(isps,:).*(1d0-omega(isps,:)),maxprec(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*(1d0-nonprec(isps,:)) > 0d0) ...
                    );
                
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6*1d0.*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                    );
                    
                for ispa = 1:nsp_aq
                    drxnsld_dmaq(isps,ispa,:) = ( ...
                        + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*msldx(isps,:).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                        + dksld_dmaq(isps,ispa,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                         + ksld(isps,:).*poro(:)'.*hr(isps,:).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*(1d0-nonprec(isps,:)) > 0d0) ...
                         + dksld_dmaq(isps,ispa,:).*poro(:)'.*hr(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*(1d0-nonprec(isps,:)) > 0d0) ...
                        );
                end 
                    
                for ispg = 1:nsp_gas
                    drxnsld_dmgas(isps,ispg,:) = ( ...
                        + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*msldx(isps,:).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                        ./(1d0-poro(:)') ...
                        + dksld_dmgas(isps,ispg,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:) < 0d0) ...
                        ./(1d0-poro(:)) ...
                         + ksld(isps,:).*poro(:)'.*hr(isps,:).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*(1d0-nonprec(isps,:)) > 0d0) ...
                         + dksld_dmgas(isps,ispg,:).*poro(:)'.*hr(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*(1d0-nonprec(isps,:)) > 0d0) ...
                        );
                end 
                
                
            case('seed')
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*(msldx(isps,:)+msld_seed).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
            
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6*1d0.*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    drxnsld_dmaq(isps,ispa,:) = ( ...
                        + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*(msldx(isps,:)+msld_seed).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmaq(isps,ispa,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*(msldx(isps,:)+msld_seed).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
                
                for ispg = 1: nsp_gas
                    drxnsld_dmgas(isps,ispg,:) = ( ...
                        + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*(msldx(isps,:)+msld_seed).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmgas(isps,ispg,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6.*(msldx(isps,:)+msld_seed).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
            
            
            case('decay')
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*msldx(isps,:).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
            
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:)*1d0.*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    % drxnsld_dmaq(isps,ispa,:) = ( ...
                        % + ksld(isps,:).*msldx(isps,:).*(-domega_dmaq(isps,ispa,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % + dksld_dmaq(isps,ispa,:).*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % );
                    % MATLAB sucks
                    for iz=1:nz
                        drxnsld_dmaq(isps,ispa,iz) = ( ...
                            + ksld(isps,iz)*msldx(isps,iz)*(-domega_dmaq(isps,ispa,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            + dksld_dmaq(isps,ispa,iz)*msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            );
                    end 
                end 
                
                for ispg = 1: nsp_gas
                    % drxnsld_dmgas(isps,ispg,:) = ( ...
                        % + ksld(isps,:).*msldx(isps,:).*(-domega_dmgas(isps,ispg,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % + dksld_dmgas(isps,ispg,:).*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % );
                    % MATLAB sucks
                    for iz=1:nz
                        drxnsld_dmgas(isps,ispg,iz) = ( ...
                            + ksld(isps,iz)*msldx(isps,iz)*(-domega_dmgas(isps,ispg,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            + dksld_dmgas(isps,ispg,iz)*msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            );
                    end 
                end 
            
            
            case('2/3noporo')
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
            
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:).*hr(isps,:)*(mv(isps)*1d-6)^(2d0/3d0) ...
                          *(2d0/3d0).*msldx(isps,:).^(-1d0/3d0).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    drxnsld_dmaq(isps,ispa,:) = ( ...
                        + ksld(isps,:).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmaq(isps,ispa,:).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
                
                for ispg = 1: nsp_gas
                    drxnsld_dmgas(isps,ispg,:) = ( ...
                        + ksld(isps,:).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmgas(isps,ispg,:).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0)*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
            
            
            case('2/3')
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.^(2d0/3d0).*hr(isps,:).*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
            
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.^(2d0/3d0).*hr(isps,:)*(mv(isps)*1d-6)^(2d0/3d0) ...
                          *(2d0/3d0).*msldx(isps,:).^(-1d0/3d0).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    drxnsld_dmaq(isps,ispa,:) = ( ...
                        + ksld(isps,:).*poro(:)'.^(2d0/3d0).*hr(isps,:) ...
                        .*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmaq(isps,ispa,:).*poro(:)'.^(2d0/3d0).*hr(isps,:) ...
                        .*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
                
                for ispg = 1: nsp_gas
                    drxnsld_dmgas(isps,ispg,:) = ( ...
                        + ksld(isps,:).*poro(:)'.^(2d0/3d0).*hr(isps,:) ...
                        .*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmgas(isps,ispg,:).*poro(:)'.^(2d0/3d0).*hr(isps,:) ...
                        .*(mv(isps)*1d-6*msldx(isps,:)).^(2d0/3d0).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
            
            
            case('psd_full')
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*hr(isps,:).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    drxnsld_dmaq(isps,ispa,:) = ( ...
                        + ksld(isps,:).*hr(isps,:).*(-domega_dmaq(isps,ispa,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmaq(isps,ispa,:).*hr(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
                
                for ispg = 1: nsp_gas
                    drxnsld_dmgas(isps,ispg,:) = ( ...
                        + ksld(isps,:).*hr(isps,:).*(-domega_dmgas(isps,ispg,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        + dksld_dmgas(isps,ispg,:).*hr(isps,:).*(1d0-omega(isps,:)) ...
                        .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        );
                end 
            
            
            otherwise
                rxnsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6 .*msldx(isps,:).*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
            
                drxnsld_dmsld(isps,:) = ( ...
                    + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6*1d0.*(1d0-omega(isps,:)) ...
                    .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                    );
                
                for ispa = 1: nsp_aq
                    % drxnsld_dmaq(isps,ispa,:) = ( ...
                        % + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6 .*msldx(isps,:).*(-domega_dmaq(isps,ispa,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % + dksld_dmaq(isps,ispa,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6 .*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % );
                    % MATLAB sucks
                    for iz=1:nz
                        drxnsld_dmaq(isps,ispa,iz) = ( ...
                            + ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6 *msldx(isps,iz)*(-domega_dmaq(isps,ispa,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            + dksld_dmaq(isps,ispa,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6 *msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            );
                    end 
                end 
                
                for ispg = 1: nsp_gas
                    % drxnsld_dmgas(isps,ispg,:) = ( ...
                        % + ksld(isps,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6 .*msldx(isps,:).*(-domega_dmgas(isps,ispg,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % + dksld_dmgas(isps,ispg,:).*poro(:)'.*hr(isps,:)*mv(isps)*1d-6 .*msldx(isps,:).*(1d0-omega(isps,:)) ...
                        % .*merge(0d0,1d0,1d0-omega(isps,:).*nonprec(isps,:) < 0d0) ...
                        % );
                    % matlab sucks
                    for iz=1:nz
                        drxnsld_dmgas(isps,ispg,iz) = ( ...
                            + ksld(isps,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6 *msldx(isps,iz)*(-domega_dmgas(isps,ispg,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            + dksld_dmgas(isps,ispg,iz)*poro(iz)*hr(isps,iz)*mv(isps)*1d-6*msldx(isps,iz)*(1d0-omega(isps,iz)) ...
                            *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0) ...
                            );
                    end 
                end 
                
                % fprintf('%7.6E\t%7.6E\n',rxnsld(isps,1),drxnsld_dmsld(isps,1));
                % fprintf('%7.6E\t%7.6E\t%7.6E\n',drxnsld_dmaq(isps,:,1));
                % fprintf('%7.6E\t%7.6E\t%7.6E\t%7.6E\t%7.6E\t%7.6E\t%7.6E\n' ...
                    % ,+ ksld(isps,1),poro(1),hr(isps,1),mv(isps),msldx(isps,1),(1d0-omega(isps,1)) ...
                    % , merge(0d0,1d0,1d0-omega(isps,1)*nonprec(isps,1) < 0d0) ...
                    % );
                % error('')
        end
    end   
   
end


function [ ... 
    flgback,w ...    % inout
    ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,nonprec,rxnsld,flx_co2sp,so4f ... % inout/out
    ] = alsilicate_aq_gas_1D_v3_1( ...
    nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext ...
    ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  ...
    ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi ...
    ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp ...
    ,stgas_ext,stgas_dext,staq_ext,stsld_ext,staq_dext,stsld_dext ...
    ,nsp_aq_all,nsp_gas_all,nsp_sld_all,nsp_aq_cnst,nsp_gas_cnst ...
    ,chraq_cnst,chraq_all,chrgas_cnst,chrgas_all,chrsld_all ...
    ,maqc,mgasc,keqgas_h,keqaq_h,keqaq_c,keqsld_all,keqaq_s,keqaq_no3,keqaq_nh3 ...
    ,nrxn_ext_all,chrrxn_ext_all,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...
    ,nsp_sld_cnst,chrsld_cnst,msldc,rho_grain,msldth_all,mv_all,staq_all,stgas_all ...
    ,turbo2,labs,trans,method_precalc,display,chrflx,sld_enforce ...% input
    ,nsld_kinspc,chrsld_kinspc,kin_sld_spc ...% input
    ,precstyle ...% input
    ,hr,poro,z,dz,w_btm,sat,pro,poroprev,tora,v,tol,it,nflx,kw,so4fprev ... %  old inputs
    ,ucv,torg,cplprec,rg,tc,sec2yr,tempk_0,proi,poroi,up,dwn,cnr,adf,msldunit  ...
    ,dt,flgback,w ...    % old inout 
    ,msldx,omega,maqx,mgasx ... % inout
    )

    % initiliza output variables
    prox=zeros(nz,1,'double');so4f=zeros(nz,1,'double');
    flx_aq=zeros(nsp_aq,nflx,nz,'double');
    flx_gas =zeros(nsp_gas,nflx,nz,'double');
    flx_sld=zeros(nsp_sld,nflx,nz,'double');
    flx_co2sp =zeros(4,nflx,nz,'double');
    rxnext=zeros(nrxn_ext,nz,'double');
    nonprec=zeros(nsp_sld,nz,'double');rxnsld=zeros(nsp_sld,nz,'double');


    % initiliza local variables
    iter=0;
    error=0;
    dgasi=zeros(nsp_gas,1,'double');
    domega_dpro=zeros(nsp_sld,nz,'double');dmsld=zeros(nsp_sld,nz,'double');dksld_dpro=zeros(nsp_sld,nz,'double');
    drxnsld_dmsld=zeros(nsp_sld,nz,'double');dksld_dso4f=zeros(nsp_sld,nz,'double');domega_dso4f=zeros(nsp_sld,nz,'double');
    domega_dmaq=zeros(nsp_sld,nsp_aq,nz,'double');dksld_dmaq=zeros(nsp_sld,nsp_aq,nz,'double');drxnsld_dmaq=zeros(nsp_sld,nsp_aq,nz,'double');
    domega_dmgas=zeros(nsp_sld,nsp_gas,nz,'double');dksld_dmgas=zeros(nsp_sld,nsp_gas,nz,'double');drxnsld_dmgas=zeros(nsp_sld,nsp_gas,nz,'double');
    dprodmaq=zeros(nsp_aq,nz,'double');dmaq=zeros(nsp_aq,nz,'double');dso4fdmaq =zeros(nsp_aq,nz,'double');
    khgasx=zeros(nsp_gas,nz,'double');khgas=zeros(nsp_gas,nz,'double');dgas=zeros(nsp_gas,nz,'double');
    agasx=zeros(nsp_gas,nz,'double');agas=zeros(nsp_gas,nz,'double');rxngas=zeros(nsp_gas,nz,'double');
    dkhgas_dpro=zeros(nsp_gas,nz,'double');dprodmgas=zeros(nsp_gas,nz,'double');dmgas=zeros(nsp_gas,nz,'double');
    dso4fdmgas=zeros(nsp_gas,nz,'double');dkhgas_dso4f=zeros(nsp_gas,nz,'double');
    dkhgas_dmaq=zeros(nsp_gas,nsp_aq,nz,'double');ddgas_dmaq=zeros(nsp_gas,nsp_aq,nz,'double');
    dagas_dmaq=zeros(nsp_gas,nsp_aq,nz,'double');drxngas_dmaq=zeros(nsp_gas,nsp_aq,nz,'double'); 
    drxngas_dmsld=zeros(nsp_gas,nsp_sld,nz,'double');
    dkhgas_dmgas=zeros(nsp_gas,nsp_gas,nz,'double');ddgas_dmgas=zeros(nsp_gas,nsp_gas,nz,'double');
    dagas_dmgas=zeros(nsp_gas,nsp_gas,nz,'double');drxngas_dmgas=zeros(nsp_gas,nsp_gas,nz,'double'); 
    drxnext_dpro=zeros(nrxn_ext,nz,'double');drxnext_dso4f=zeros(nrxn_ext,nz,'double');
    drxnext_dmgas=zeros(nrxn_ext,nsp_gas,nz,'double');
    drxnext_dmaq=zeros(nrxn_ext,nsp_aq,nz,'double');
    drxnext_dmsld=zeros(nrxn_ext,nsp_sld,nz,'double');

    dprodmaq_all=zeros(nsp_aq_all,nz,'double');dso4fdmaq_all=zeros(nsp_aq_all,nz,'double');
    dprodmgas_all=zeros(nsp_gas_all,nz,'double');dso4fdmgas_all=zeros(nsp_gas_all,nz,'double');

    domega_dpro_loc=zeros(nz,1,'double');domega_dso4f_loc=zeros(nz,1,'double');
    domega_dmgas_all=zeros(nsp_gas_all,nz,'double');
    domega_dmaq_all=zeros(nsp_aq_all,nz,'double');

    mgasx_loc=zeros(nsp_gas_all,nz,'double');
    khgas_all=zeros(nsp_gas_all,nz,'double');khgasx_all=zeros(nsp_gas_all,nz,'double');
    dkhgas_dpro_all=zeros(nsp_gas_all,nz,'double');dkhgas_dso4f_all=zeros(nsp_gas_all,nz,'double');
    dkhgas_dmaq_all=zeros(nsp_gas_all,nsp_aq_all,nz,'double');
    dkhgas_dmgas_all=zeros(nsp_gas_all,nsp_gas_all,nz,'double');

    [ieqgas_h0,ieqgas_h1,ieqgas_h2]=deal(1,2,3);
    [ieqaq_h1,ieqaq_h2,ieqaq_h3,ieqaq_h4]=deal(1,2,3,4);
    [ieqaq_co3,ieqaq_hco3]=deal(1,2);
    [ieqaq_so4,ieqaq_so42]=deal(1,2);

    iz=0;row=0;ie=0;ie2=0;iflx=0;isps=0;ispa=0;ispg=0;ispa2=0;ispg2=0;col=0;irxn=0;isps2=0;iiz=0;isps_kinspc=0;row_w=0;col_w=0;
    itflx=0;iadv=0;idif=0;irain=0;ires=0;
    ph_iter=0;ph_iter2=0;
    [itflx,iadv,idif,irain]=deal(1,2,3,4);

    irxn_sld= zeros( nsp_sld,1,'int32');
    irxn_ext =zeros( nrxn_ext,1,'int32');

    d_tmp=0;caq_tmp=0;caq_tmp_p=0;caq_tmp_n=0;caqth_tmp=0;caqi_tmp=0;rxn_tmp=0;caq_tmp_prev=0;drxndisp_tmp=0; 
    k_tmp=0;mv_tmp=0;omega_tmp=0;m_tmp=0;mth_tmp=0;mi_tmp=0;mp_tmp=0;msupp_tmp=0;mprev_tmp=0;omega_tmp_th=0;rxn_ext_tmp=0;
    edif_tmp=0;edif_tmp_n=0;edif_tmp_p=0;khco2n_tmp=0;pco2n_tmp=0;edifn_tmp=0;caqsupp_tmp=0;kco2=0;k1=0;k2=0;kho=0;sw_red=0;
    flx_max=0;flx_max_max=0;proi_tmp=0;knh3=0;k1nh3=0;kn2o=0;wp_tmp=0;w_tmp=0;sporo_tmp=0;sporop_tmp=0;sporoprev_tmp=0;
    mn_tmp=0;wn_tmp=0;sporon_tmp=0;

    infinity = Inf;
    fact = 1d-3;
    dconc = 1d-14;
    maxfact = 1d200;
    % threshold = log(maxfact);
    threshold = 10d0;
    % threshold = 3d0;
    % corr = 1.5d0;
    corr = exp(threshold);

    dummy=zeros(nz,1,'double');dummy2=zeros(nz,1,'double');dummy3=zeros(nz,1,'double');kin=zeros(nz,1,'double');
    dkin_dmsp=zeros(nz,1,'double');dumtest=zeros(nz,1,'double');sporo=zeros(nz,1,'double');

    print_cb=false;ph_error=false;omega_error=false;rxnext_error=false;

    iter_max = 50;
    % iter_max = 300;

    nz_disp = 10;

    amx3=zeros(nsp3*nz,nsp3*nz,'double'); ymx3=zeros(nsp3*nz,1,'double');emx3=zeros(nsp3*nz,1,'double');
    xmx3=zeros(nsp3*nz,1,'double');
    rmx3=0;
    ipiv3=zeros(nsp3*nz,1,'int32');
    info =0;

    chkflx = true;
    dt_norm = true;
    kin_iter = true;
    new_gassol = true;
    % new_gassol = false;

    % sld_enforce = false;
    msld_seed=0;fact2=0;
    fact_tol = 1d-3;
    % fact_tol = 1d-4;
    dt_th = 1d-6;
    flx_tol = 1d-4;%= tol*fact_tol*(z(nz)+0.5d0*dz(nz))
    % flx_tol = 1d-3 ;% desparate to make things converge 
    % flx_max_tol = 1d-9 ;%= tol*fact_tol*(z(nz)+0.5d0*dz(nz)) % working for most cases but not when spinup with N cycles
    flx_max_tol = 1d-6 ;%= tol*fact_tol*(z(nz)+0.5d0*dz(nz)) 
    solve_sld =0;

    %-----------------------------------------------
        
    % dlmwrite('msld_matlab.txt',msldx);
    % dlmwrite('maq_matlab.txt',maqx);
    
    % error('printing_stuff');
    % pause;

    msld_seed = 1d-20;

    if (sld_enforce)  
        solve_sld = 0;
    else
        solve_sld = 1;
    end 

    sw_red = 1d0;
    if (~method_precalc); sw_red = -1d100; end
    sw_red = -1d100;

    for isps=1:nsp_sld
        irxn_sld(isps) = 4+isps;
    end 

    for irxn=1:nrxn_ext
        irxn_ext(irxn) = 4+nsp_sld+irxn;
    end 

    ires = nflx;

    print_cb = false; 
    print_loc = './ph.txt';

    kco2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h0);
    k1 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h1);
    k2 = keqgas_h(find(chrgas_all=='pco2'),ieqgas_h2);

    kho = keqgas_h(find(chrgas_all=='po2'),ieqgas_h0);

    knh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h0);
    k1nh3 = keqgas_h(find(chrgas_all=='pnh3'),ieqgas_h1);

    kn2o = keqgas_h(find(chrgas_all=='pn2o'),ieqgas_h0);

    sporo(:) = 1d0 - poro(:);
    if (msldunit=='blk'); sporo(:) = 1d0; end
        
    nonprec(:,:) = 1d0; % primary minerals only dissolve
    if (cplprec)
        for isps = 1: nsp_sld
            if (any(chrsld_2 == chrsld(isps)))   
                nonprec(isps,:) = 0d0; % allowing precipitation for secondary phases
            end 
        end
    end 

    if (any(isnan(tora))) 
        error([repmat('%7.6E \t',1,nz) '\n'], tora(1:nz));
    end 

    dummy(:) = 0d0;
    dummy2(:) = 0d0;

    error = 1d4;
    iter = 0;

    while ((~isnan(error))&&(error > tol*fact_tol))

        amx3(:,:)=0.0d0;
        ymx3(:)=0.0d0;
        emx3(:)=0.0d0; 
        xmx3(:)=0.0d0; 
        rmx3=0.0d0; 
        
        flx_sld(:,:,:) = 0d0;
        flx_aq(:,:,:) = 0d0;
        flx_gas(:,:,:) = 0d0;
        
        % pH calculation and its derivative wrt aq and gas species
        
        [ ...
            dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
            ,prox,ph_error,so4f,ph_iter ...% output
            ] = calc_pH_v7_3( ...
            nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
            ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
            ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
            ,print_cb,print_loc,z ...% input 
            ,prox,so4f ...% inout
            ); 

        if (ph_error)  
            flgback = true;
            return
        end 
        
        dprodmaq(:,:) = 0d0;
        dso4fdmaq(:,:) = 0d0;
        for ispa=1:nsp_aq
            if (any (chraq_ph == chraq(ispa)))  
                dprodmaq(ispa,:)=dprodmaq_all(find(chraq_all==chraq(ispa)),:);
                dso4fdmaq(ispa,:)=dso4fdmaq_all(find(chraq_all==chraq(ispa)),:);
            end 
        end 
        
        % dlmwrite('dprodmaq.txt',dprodmaq);
        % dlmwrite('dso4fdmaq.txt',dso4fdmaq);
        % error('printing_stuff');
        % pause;
        
        dprodmgas(:,:) = 0d0;
        dso4fdmgas(:,:) = 0d0;
        for ispg=1:nsp_gas
            if (any (chrgas_ph == chrgas(ispg)))  
                dprodmgas(ispg,:)=dprodmgas_all(find(chrgas_all==chrgas(ispg)),:);
                dso4fdmGas(ispg,:)=dso4fdmgas_all(find(chrgas_all==chrgas(ispg)),:);
            end 
        end 
        
        % recalculation of rate constants for mineral reactions
        if (kin_iter)  
            ksld(:,:) = 0d0;
            dksld_dpro(:,:) = 0d0;
            dksld_dso4f(:,:) = 0d0;
            dksld_dmaq(:,:,:) = 0d0;
            dksld_dmgas(:,:,:) = 0d0;
            
            [mgasx_loc] = get_mgasx_all( ...
                nz,nsp_gas_all,nsp_gas,nsp_gas_cnst ...
                ,chrgas,chrgas_all,chrgas_cnst ...
                ,mgasx,mgasc ...
                );
            
            for isps =1:nsp_sld 
                [kin,dkin_dmsp] = sld_kin( ...
                    nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) ...% input
                    ,nsp_gas_all,chrgas_all,mgasx_loc ...% input
                    ,chrsld(isps),'pro' ...% input 
                    ); 
                ksld(isps,:) = kin(:);
                dksld_dpro(isps,:) = dkin_dmsp(:);
                
                for ispa = 1:nsp_aq
                    if (any (chraq_ph == chraq(ispa)) || staq(isps,ispa)~=0d0 )  
                        [kin,dkin_dmsp] = sld_kin( ...
                            nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) ...% input
                            ,nsp_gas_all,chrgas_all,mgasx_loc ...% input
                            ,chrsld(isps),chraq(ispa) ...% input 
                            ); 
                        % dksld_dmaq(isps,ispa,:) = dkin_dmsp(:)' + ( ...
                            % dksld_dpro(isps,:).*dprodmaq(ispa,:) ...
                            % +dksld_dso4f(isps,:).*dso4fdmaq(ispa,:) ...
                            % );
                        for iz=1:nz
                            dksld_dmaq(isps,ispa,iz) = dkin_dmsp(iz) + ( ...
                                dksld_dpro(isps,iz)*dprodmaq(ispa,iz) ...
                                +dksld_dso4f(isps,iz)*dso4fdmaq(ispa,iz) ...
                                );
                        end 
                    end 
                end 
                
                for ispg = 1:nsp_gas
                    if (any (chrgas_ph == chrgas(ispg)) || stgas(isps,ispg)~=0d0)  
                        [kin,dkin_dmsp] = sld_kin( ...
                            nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) ...% input
                            ,nsp_gas_all,chrgas_all,mgasx_loc ...% input
                            ,chrsld(isps),chrgas(ispg) ...% input 
                            ) ;
                        % dksld_dmgas(isps,ispg,:) = dkin_dmsp(:)' + ( ...
                            % dksld_dpro(isps,:).*dprodmgas(ispg,:) ...
                            % +dksld_dso4f(isps,:).*dso4fdmgas(ispg,:) ...
                            % );
                        for iz=1:nz
                            dksld_dmgas(isps,ispg,iz) = dkin_dmsp(iz) + ( ...
                                dksld_dpro(isps,iz)*dprodmgas(ispg,iz) ...
                                +dksld_dso4f(isps,iz)*dso4fdmgas(ispg,iz) ...
                                );
                        end 
                    end 
                end 
            
            end 
        
            % dlmwrite('ksld.txt',ksld);
            % dlmwrite('dksld_dpro.txt',dksld_dpro);
            % error('printing_stuff');
            % pause;
        
        else 
            dksld_dpro(:,:) = 0d0;
            dksld_dso4f(:,:) = 0d0;
            dksld_dmaq(:,:,:) = 0d0;
            dksld_dmgas(:,:,:) = 0d0;
        end 
        
        % if kin const. is specified in input file 
        if (nsld_kinspc > 0)  
            for isps_kinspc=1:nsld_kinspc    
                if ( any( chrsld == chrsld_kinspc(isps_kinspc)))  
                    switch (chrsld_kinspc(isps_kinspc))
                        case{'g1','g2','g3'} % for OMs, turn over year needs to be provided [yr]
                            ksld(find(chrsld,chrsld_kinspc(isps_kinspc)),:) = ( ...                   
                                1d0/kin_sld_spc(isps_kinspc) ...
                                ) ;
                            dksld_dpro(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                            dksld_dso4f(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                            dksld_dmaq(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                            dksld_dmgas(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                        otherwise % otherwise, usual rate constant [mol/m2/yr]
                            ksld(find(chrsld,chrsld_kinspc(isps_kinspc)),:) = ( ...                            
                                kin_sld_spc(isps_kinspc) ...
                                ) ;
                            dksld_dpro(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                            dksld_dso4f(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                            dksld_dmaq(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                            dksld_dmgas(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                    end
                end 
            end 
        end 
                        
        
        % saturation state calc. and their derivatives wrt aq and gas species
        
        omega(:,:) = 0d0;
        domega_dpro(:,:) = 0d0;
        domega_dso4f(:,:) = 0d0;
        domega_dmaq(:,:,:) = 0d0;
        domega_dmgas(:,:,:) = 0d0;
        
        for isps =1: nsp_sld
            
            dummy(:) = 0d0;
            domega_dpro_loc(:) = 0d0;
            domega_dso4f_loc(:) = 0d0;
            [ ... 
                domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dso4f_loc ...% output
                ,dummy,omega_error ...% output
                ] = calc_omega_v4( ...
                nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ... 
                ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all ...
                ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,so4f ...
                ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
                ,staq_all,stgas_all ...
                ,chrsld(isps) ...
                );
            if (omega_error) 
                flgback = true;
                return 
            end 
            omega(isps,:) = dummy(:);
            domega_dpro(isps,:) = domega_dpro_loc(:);
            domega_dso4f(isps,:) = domega_dso4f_loc(:);
            
            for ispa = 1: nsp_aq
                if (any (chraq_ph == chraq(ispa)) || staq(isps,ispa)~=0d0 )  
                
                    % domega_dmaq(isps,ispa,:) = domega_dmaq_all(find(chraq_all==chraq(ispa)),:)+ ( ...
                        % domega_dpro(isps,:).*dprodmaq(ispa,:) ...
                        % +domega_dso4f(isps,:).*dso4fdmaq(ispa,:) ...
                        % );
                    for iz=1:nz
                        domega_dmaq(isps,ispa,iz) = domega_dmaq_all(find(chraq_all==chraq(ispa)),iz)+ ( ...
                            domega_dpro(isps,iz)*dprodmaq(ispa,iz) ...
                            +domega_dso4f(isps,iz)*dso4fdmaq(ispa,iz) ...
                            );
                    end 
                end 
            end
            
            for ispg = 1: nsp_gas
                if (any (chrgas_ph == chrgas(ispg)) || stgas(isps,ispg)~=0d0)  
                
                    % domega_dmgas(isps,ispg,:) = domega_dmgas_all(find(chrgas_all==chrgas(ispg)),:)+ ( ...
                        % domega_dpro(isps,:).*dprodmgas(ispg,:) ...
                        % +domega_dso4f(isps,:).*dso4fdmgas(ispg,:) ...
                        % );
                    for iz=1:nz
                        domega_dmgas(isps,ispg,iz) = domega_dmgas_all(find(chrgas_all==chrgas(ispg)),iz)+ ( ...
                            domega_dpro(isps,iz)*dprodmgas(ispg,iz) ...
                            +domega_dso4f(isps,iz)*dso4fdmgas(ispg,iz) ...
                            );
                    end 
                end 
            end
        end 
        
        % dlmwrite('omega.txt',omega');
        % dlmwrite('domega_dpro.txt',domega_dpro');
        % error('printing_stuff');
        % pause;
        
        % adding reactions that are not based on dis/prec of minerals
        
        rxnext(:,:) = 0d0;
        drxnext_dpro(:,:) = 0d0;
        drxnext_dso4f(:,:) = 0d0;
        drxnext_dmaq(:,:,:) = 0d0;
        drxnext_dmgas(:,:,:) = 0d0;
        drxnext_dmsld(:,:,:) = 0d0;
        
        for irxn=1:nrxn_ext
            dummy(:) = 0d0;
            dummy2(:) = 0d0;
            [ ... 
                dummy,dummy2,rxnext_error ...% output
                ] = calc_rxn_ext_dev_2( ...
                nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
                ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
                ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
                ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
                ,rg,tempk_0,tc ...%input
                ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
                ,chrrxn_ext(irxn),'pro' ...% input 
                );
            if (rxnext_error) 
                flgback = true;
                return 
            end 
            rxnext(irxn,:) = dummy(:)';
            drxnext_dpro(irxn,:) = dummy2(:)';
            
            dummy(:) = 0d0;
            dummy2(:) = 0d0;
            [ ... 
                dummy,dummy2,rxnext_error ...% output
                ] = calc_rxn_ext_dev_2( ...
                nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
                ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
                ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
                ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
                ,rg,tempk_0,tc ...%input
                ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
                ,chrrxn_ext(irxn),'so4f ' ...% input 
                );
            if (rxnext_error) 
                flgback = true;
                return 
            end 
            drxnext_dso4f(irxn,:) = dummy2(:)';
            
            for ispg=1:nsp_gas
                if (stgas_dext(irxn,ispg)==0d0); continue; end
                
                dummy(:) = 0d0;
                dummy2(:) = 0d0;
                [ ... 
                    dummy,dummy2,rxnext_error ...% output
                    ] = calc_rxn_ext_dev_2( ...
                    nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
                    ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
                    ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
                    ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
                    ,rg,tempk_0,tc ...%input
                    ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
                    ,chrrxn_ext(irxn),chrgas(ispg) ...% input 
                    );
                if (rxnext_error) 
                    flgback = true;
                    return 
                end 
                % drxnext_dmgas(irxn,ispg,:) = dummy2(:)' + (...
                    % + drxnext_dpro(irxn,:).*dprodmgas(ispg,:) ...
                    % + drxnext_dso4f(irxn,:).*dso4fdmgas(ispg,:) ...
                    % );
                for iz=1:nz
                    drxnext_dmgas(irxn,ispg,iz) = dummy2(iz) + (...
                        + drxnext_dpro(irxn,iz)*dprodmgas(ispg,iz) ...
                        + drxnext_dso4f(irxn,iz)*dso4fdmgas(ispg,iz) ...
                        );
                end 
            end 
            
            for ispa=1:nsp_aq
                if (staq_dext(irxn,ispa)==0d0); continue; end
                
                dummy(:) = 0d0;
                dummy2(:) = 0d0;
                [ ...
                    ,dummy,dummy2,rxnext_error ...% output
                    ] = calc_rxn_ext_dev_2( ...
                    nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
                    ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
                    ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
                    ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
                    ,rg,tempk_0,tc ...%input
                    ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
                    ,chrrxn_ext(irxn),chraq(ispa) ...% input 
                    );
                if (rxnext_error) 
                    flgback = true;
                    return 
                end 
                % drxnext_dmaq(irxn,ispa,:) = dummy2(:)' + ( ...
                    % + drxnext_dpro(irxn,:)*dprodmaq(ispa,:) ...
                    % + drxnext_dso4f(irxn,:)*dso4fdmaq(ispa,:) ...
                    % );
                for iz=1:nz
                    drxnext_dmaq(irxn,ispa,iz) = dummy2(iz) + ( ...
                        + drxnext_dpro(irxn,iz)*dprodmaq(ispa,iz) ...
                        + drxnext_dso4f(irxn,iz)*dso4fdmaq(ispa,iz) ...
                        );
                end 
            end 
            
            for isps=1:nsp_sld
                if (stsld_dext(irxn,isps)==0d0); continue; end
                
                dummy(:) = 0d0;
                dummy2(:) = 0d0;
                [ ... 
                    dummy,dummy2,rxnext_error ...% output
                    ] = calc_rxn_ext_dev_2( ...
                    nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
                    ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
                    ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
                    ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
                    ,rg,tempk_0,tc ...%input
                    ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
                    ,chrrxn_ext(irxn),chrsld(isps) ...% input 
                    );
                if (rxnext_error) 
                    flgback = true;
                    return 
                end 
                % drxnext_dmsld(irxn,isps,:) = dummy2(:)';
                for iz=1:nz
                    drxnext_dmsld(irxn,isps,iz) = dummy2(iz);
                end 
            end 
        end 
        
        % gas tansport
        khgas(:,:) = 0d0;
        khgasx(:,:) = 0d0;
        dkhgas_dmaq(:,:,:) = 0d0;
        dkhgas_dmgas(:,:,:) = 0d0;
        % added
        dkhgas_dpro(:,:) = 0d0;
        dkhgas_dso4f(:,:) = 0d0;
        
        if (new_gassol)  
            [ ...
                khgas_all,khgasx_all,dkhgas_dpro_all,dkhgas_dso4f_all,dkhgas_dmaq_all,dkhgas_dmgas_all ...%output
                ] = calc_khgas_all( ...
                nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst ...
                ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas ...
                ,maq,mgas,maqx,mgasx,maqc,mgasc ...
                ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3  ...
                ,pro,prox,so4fprev,so4f ...
                );
            
            for ispg=1:nsp_gas
                khgas(ispg,:)=khgas_all(find(chrgas_all==chrgas(ispg)),:);
                khgasx(ispg,:)=khgasx_all(find(chrgas_all==chrgas(ispg)),:);
                dkhgas_dpro(ispg,:)=dkhgas_dpro_all(find(chrgas_all==chrgas(ispg)),:);
                dkhgas_dso4f(ispg,:)=dkhgas_dso4f_all(find(chrgas_all==chrgas(ispg)),:);
                for ispa=1:nsp_aq
                    % dkhgas_dmaq(ispg,ispa,:)= ...
                        % dkhgas_dmaq_all(find(chrgas_all==chrgas(ispg)),find(chraq_all==chraq(ispa)),:) ...
                        % + dkhgas_dpro(ispg,:).*dprodmaq(ispa,:) ...
                        % + dkhgas_dso4f(ispg,:).*dso4fdmaq(ispa,:);
                    % MATLAB sucks
                    for iz=1:nz
                        dkhgas_dmaq(ispg,ispa,iz)= ...
                            dkhgas_dmaq_all(find(chrgas_all==chrgas(ispg)),find(chraq_all==chraq(ispa)),iz) ...
                            + dkhgas_dpro(ispg,iz)*dprodmaq(ispa,iz) ...
                            + dkhgas_dso4f(ispg,iz)*dso4fdmaq(ispa,iz);
                    end 
                end 
                for ispg2=1:nsp_gas
                    % dkhgas_dmgas(ispg,ispg2,:)= ...
                        % dkhgas_dmgas_all(find(chrgas_all==chrgas(ispg)),find(chrgas_all==chrgas(ispg2)),:) ...
                        % + dkhgas_dpro(ispg,:).*dprodmgas(ispg2,:) ...
                        % + dkhgas_dso4f(ispg,:).*dso4fdmgas(ispg2,:);
                    % MATLAB sucks
                    for iz=1:nz
                        dkhgas_dmgas(ispg,ispg2,iz)= ...
                            dkhgas_dmgas_all(find(chrgas_all==chrgas(ispg)),find(chrgas_all==chrgas(ispg2)),iz) ...
                            + dkhgas_dpro(ispg,iz)*dprodmgas(ispg2,iz) ...
                            + dkhgas_dso4f(ispg,iz)*dso4fdmgas(ispg2,iz);
                    end 
                end 
            end 
        end
        
        dgas(:,:) = 0d0;
        ddgas_dmaq(:,:,:) = 0d0;
        ddgas_dmgas(:,:,:) = 0d0;
        
        agas(:,:) = 0d0;
        agasx(:,:) = 0d0;
        dagas_dmaq(:,:,:) = 0d0;
        dagas_dmgas(:,:,:) = 0d0;
        
        for ispg = 1: nsp_gas
            
            if (~ new_gassol)  % old way to calc solubility (to be removed?)
                switch (chrgas(ispg))
                    case('pco2')
                        khgas(ispg,:) = kco2*(1d0+k1./pro(:)' + k1*k2./pro(:)'./pro(:)'); % previous value; should not change through iterations 
                        khgasx(ispg,:) = kco2*(1d0+k1./prox(:)' + k1*k2./prox(:)'./prox(:)');
                
                        dkhgas_dpro(ispg,:) = kco2*(k1*(-1d0)./prox(:)'.^2d0 + k1*k2*(-2d0)./prox(:)'.^3d0);
                    case('po2')
                        khgas(ispg,:) = kho; % previous value; should not change through iterations 
                        khgasx(ispg,:) = kho;
                
                        dkhgas_dpro(ispg,:) = 0d0;
                    case('pnh3')
                        khgas(ispg,:) = knh3*(1d0+pro(:)/k1nh3); % previous value; should not change through iterations 
                        khgasx(ispg,:) = knh3*(1d0+prox(:)/k1nh3);
                
                        dkhgas_dpro(ispg,:) = knh3*(1d0/k1nh3);
                    case('pn2o')
                        khgas(ispg,:) = kn2o; % previous value; should not change through iterations 
                        khgasx(ispg,:) = kn2o;
                
                        dkhgas_dpro(ispg,:) = 0d0;
                end
            end 
            
            dgas(ispg,:) = ucv*poro(:).*(1.0d0-sat(:))*1d3.*torg(:)*dgasg(ispg)+poro(:).*sat(:).*khgasx(ispg,:)'*1d3.*tora(:)*dgasa(ispg);
            dgasi(ispg) = ucv*1d3*dgasg(ispg) ;
            
            agas(ispg,:)= ucv*poroprev(:).*(1.0d0-sat(:))*1d3+poroprev(:).*sat(:).*khgas(ispg,:)'*1d3;
            agasx(ispg,:)= ucv*poro(:).*(1.0d0-sat(:))*1d3+poro(:).*sat(:).*khgasx(ispg,:)'*1d3;
            
            for ispa = 1:nsp_aq 
                % if (~ new_gassol); dkhgas_dmaq(ispg,ispa,:) = dkhgas_dpro(ispg,:).*dprodmaq(ispa,:); end % old way to calc solubility (to be removed?)
                % ddgas_dmaq(ispg,ispa,:) = poro(:)'.*sat(:)'.*dkhgas_dmaq(ispg,ispa,:)*1d3.*tora(:)'*dgasa(ispg);
                % dagas_dmaq(ispg,ispa,:) =  poro(:)'.*sat(:)'.*dkhgas_dmaq(ispg,ispa,:)*1d3;
                % MATLAB sucks
                for iz=1:nz
                    if (~ new_gassol); dkhgas_dmaq(ispg,ispa,iz) = dkhgas_dpro(ispg,iz)*dprodmaq(ispa,iz); end % old way to calc solubility (to be removed?)
                    ddgas_dmaq(ispg,ispa,iz) = poro(iz)*sat(iz)*dkhgas_dmaq(ispg,ispa,iz)*1d3*tora(iz)*dgasa(ispg);
                    dagas_dmaq(ispg,ispa,iz) =  poro(iz)*sat(iz)*dkhgas_dmaq(ispg,ispa,iz)*1d3;
                end 
            end 
            
            for ispg2 = 1:nsp_gas 
                % if (~ new_gassol); dkhgas_dmgas(ispg,ispg2,:) = dkhgas_dpro(ispg,:).*dprodmgas(ispg2,:); end % old way to calc solubility (to be removed?)
                % ddgas_dmgas(ispg,ispg2,:) = poro(:)'.*sat(:)'.*dkhgas_dmgas(ispg,ispg2,:)*1d3.*tora(:)'*dgasa(ispg);
                % dagas_dmgas(ispg,ispg2,:) =  poro(:)'.*sat(:)'.*dkhgas_dmgas(ispg,ispg2,:)*1d3;
                % MATLAB sucks
                for iz=1:nz
                    if (~ new_gassol); dkhgas_dmgas(ispg,ispg2,iz) = dkhgas_dpro(ispg,iz)*dprodmgas(ispg2,iz); end % old way to calc solubility (to be removed?)
                    ddgas_dmgas(ispg,ispg2,iz) = poro(iz)*sat(iz)*dkhgas_dmgas(ispg,ispg2,iz)*1d3*tora(iz)*dgasa(ispg);
                    dagas_dmgas(ispg,ispg2,iz) =  poro(iz)*sat(iz)*dkhgas_dmgas(ispg,ispg2,iz)*1d3;
                end 
            end 
        end 
        
        % sld phase reactions
        
        rxnsld(:,:) = 0d0;
        drxnsld_dmsld(:,:) = 0d0;
        drxnsld_dmaq(:,:,:) = 0d0;
        drxnsld_dmgas(:,:,:) = 0d0;
        
        [ ... 
            rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas ...% output
            ] = sld_rxn( ...
            nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz ...% input 
            ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle ...% input
            ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas ...% input
            ); 
        
        % dlmwrite('rxnsld.txt',rxnsld');
        % dlmwrite('drxnsld_dmsld.txt',drxnsld_dmsld');
        % error('printing_stuff');
        % pause;
        
        % gas reactions 
        
        rxngas(:,:) = 0d0;
        drxngas_dmaq(:,:,:) = 0d0;
        drxngas_dmsld(:,:,:) = 0d0;
        drxngas_dmgas(:,:,:) = 0d0;
            
        for ispg = 1: nsp_gas
            for isps = 1: nsp_sld
                rxngas(ispg,:) =  rxngas(ispg,:) + (...
                    + stgas(isps,ispg)*rxnsld(isps,:) ...
                    );
                % drxngas_dmsld(ispg,isps,:) =  drxngas_dmsld(ispg,isps,:) + (...
                    % + stgas(isps,ispg)*drxnsld_dmsld(isps,:) ...
                    % );
                % MATLAB sucks
                for iz=1:nz
                    drxngas_dmsld(ispg,isps,iz) =  drxngas_dmsld(ispg,isps,iz) + (...
                        + stgas(isps,ispg)*drxnsld_dmsld(isps,iz) ...
                        );
                end 
                % change done
                for ispg2 = 1:nsp_gas
                    drxngas_dmgas(ispg,ispg2,:) =  drxngas_dmgas(ispg,ispg2,:) + (...
                        + stgas(isps,ispg)*drxnsld_dmgas(isps,ispg2,:) ...
                        );
                end 
                for ispa = 1:nsp_aq
                    drxngas_dmaq(ispg,ispa,:) =  drxngas_dmaq(ispg,ispa,:) + ( ...
                        + stgas(isps,ispg)*drxnsld_dmaq(isps,ispa,:) ...
                        );
                end 
            end 
        end 
                
        if (~sld_enforce)  

            for iz = 1: nz  %================================
                
                for isps = 1: nsp_sld
                
                    row = nsp3*(iz-1)+isps;
                    
                    k_tmp = ksld(isps,iz);
                    mv_tmp = mv(isps);
                    omega_tmp = omega(isps,iz);
                    omega_tmp_th = omega_tmp*nonprec(isps,iz);
                    m_tmp = msldx(isps,iz) ;
                    mth_tmp = msldth(isps) ;
                    mi_tmp = msldi(isps);
                    mp_tmp = msldx(isps,min(nz,iz+1));
                    msupp_tmp = msldsupp(isps,iz) ;
                    rxn_ext_tmp = sum(stsld_ext(:,isps).*rxnext(:,iz));
                    mprev_tmp = msld(isps,iz)  ;
                    w_tmp = w(iz) ;
                    wp_tmp = w(min(nz,iz+1)) ;
                    sporo_tmp = 1d0-poro(iz);
                    sporop_tmp = 1d0-poro(min(nz,iz+1)) ;
                    sporoprev_tmp = 1d0-poroprev(iz);
                    mn_tmp = msldx(isps,max(1,iz-1));
                    wn_tmp = w(max(1,iz-1));
                    sporon_tmp = 1d0-poro(max(1,iz-1));
                    
                    if (iz==1)  
                        mn_tmp = 0d0;
                        wn_tmp = 0d0;
                        sporon_tmp = 0d0;
                    end 
                    
                    if (iz==nz)  
                        mp_tmp = mi_tmp;
                        wp_tmp = w_btm ;
                        sporop_tmp = 1d0- poroi;
                    end 
                    
                    if (msldunit == 'blk')  
                        sporo_tmp = 1d0;
                        sporop_tmp = 1d0;
                        sporon_tmp = 1d0;
                        sporoprev_tmp = 1d0;
                    end 

                    amx3(row,row) = ( ...
                        1d0 *  sporo_tmp /merge(1d0,dt,dt_norm)     ...
                        + sporo_tmp*w_tmp/dz(iz)*merge(dt,1d0,dt_norm)    ...
                        + drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) ...
                        - sum(stsld_ext(:,isps).*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) ...
                        ) ...
                        * merge(1.0d0,m_tmp,m_tmp<mth_tmp*sw_red);

                    ymx3(row) = ( ...
                        ( sporo_tmp*m_tmp - sporoprev_tmp*mprev_tmp )/merge(1d0,dt,dt_norm) ...
                        - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  ...
                        + rxnsld(isps,iz)*merge(dt,1d0,dt_norm) ...
                        -msupp_tmp*merge(dt,1d0,dt_norm)  ...
                        -rxn_ext_tmp*merge(dt,1d0,dt_norm)  ...
                        ) ...
                        *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                    
                    % disp(ymx3(row));
                    % fprintf('%7.6E\t%7.6E\t%7.6E\t%7.6E\t%7.6E\n' ...
                        % ,( sporo_tmp*m_tmp - sporoprev_tmp*mprev_tmp )/merge(1d0,dt,dt_norm) ...
                        % ,- ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)*merge(dt,1d0,dt_norm)  ...
                        % ,+ rxnsld(isps,iz)*merge(dt,1d0,dt_norm) ...
                        % ,-msupp_tmp*merge(dt,1d0,dt_norm)  ...
                        % ,-rxn_ext_tmp*merge(dt,1d0,dt_norm)  ...
                        % );
                    % fprintf('%7.6E\t%7.6E\n' ...
                        % ,+ rxnsld(isps,iz),merge(dt,1d0,dt_norm) ...
                        % );
                    % if row == 2; error('chk'); end
                    % pause;
                        
                    if (iz~=nz) 
                        amx3(row,row+nsp3) = ( ...
                            (- sporop_tmp*wp_tmp/dz(iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *merge(1.0d0,mp_tmp,m_tmp<mth_tmp*sw_red);
                    end
                    
                    for ispa = 1: nsp_aq
                        col = nsp3*(iz-1) + nsp_sld + ispa;
                        
                        amx3(row,col ) = ( ...
                            + drxnsld_dmaq(isps,ispa,iz)*merge(dt,1d0,dt_norm) ...
                            - sum(stsld_ext(:,isps).*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *maqx(ispa,iz) ...
                            *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                    end 
                    
                    for ispg = 1: nsp_gas 
                        col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg;

                        amx3(row,col) = ( ...
                            + drxnsld_dmgas(isps,ispg,iz)*merge(dt,1d0,dt_norm) ...
                            - sum(stsld_ext(:,isps).*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *mgasx(ispg,iz) ...
                            *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                    end 
                    
                    for isps2 = 1:nsp_sld 
                        if (isps2 == isps); continue; end
                        col = nsp3*(iz-1)+ isps2;

                        amx3(row,col) = ( ...
                            - sum(stsld_ext(:,isps).*drxnext_dmsld(:,isps2,iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *msldx(isps2,iz) ...
                            *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                    end 
                    
                    % modifications with porosity and dz are made in make_trans subroutine
                    for iiz = 1: nz
                        col = nsp3*(iiz-1)+isps;
                        if (trans(iiz,iz,isps)==0d0); continue; end
                            
                        amx3(row,col) = amx3(row,col) - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) ...
                            *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                        ymx3(row) = ymx3(row) - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz)* merge(dt,1d0,dt_norm) ...
                            *merge(0.0d0,1d0,m_tmp<mth_tmp*sw_red);
                            
                        flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( ...
                            - trans(iiz,iz,isps)*msldx(isps,iiz)* sporo(iiz) ...
                            );
                    end
                    
                    % disp(ymx3(row));
                    % error('chk');
                    % pause;
                    
                    flx_sld(isps,itflx,iz) = ( ...
                        ( sporo_tmp*m_tmp- sporoprev_tmp*mprev_tmp)/dt ...
                        );
                    flx_sld(isps,iadv,iz) = ( ...
                        - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  ...
                        );
                    flx_sld(isps,irxn_sld(isps),iz) = ( ...
                        + rxnsld(isps,iz) ...
                        );
                    flx_sld(isps,irain,iz) = (...
                        - msupp_tmp  ...
                        );
                    flx_sld(isps,irxn_ext(:),iz) = (...
                        - stsld_ext(:,isps).*rxnext(:,iz)  ...
                        );
                    flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz));
                    if (isnan(flx_sld(isps,ires,iz)))  
                        warning( ['nan in flx_sld %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chrsld(isps),iz,flx_sld(isps,1:nflx,iz) );
                    end 
                end 
            end  %================================
        
        end 
        

        for iz = 1: nz
            
            for ispa = 1: nsp_aq

                row = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa;
                
                d_tmp = daq(ispa);
                caq_tmp = maqx(ispa,iz);
                caq_tmp_prev = maq(ispa,iz);
                caq_tmp_p = maqx(ispa,min(nz,iz+1));
                caq_tmp_n = maqx(ispa,max(1,iz-1));
                caqth_tmp = maqth(ispa);
                caqi_tmp = maqi(ispa);
                caqsupp_tmp = maqsupp(ispa,iz) ;
                rxn_ext_tmp = sum(staq_ext(:,ispa).*rxnext(:,iz));
                rxn_tmp = sum(staq(:,ispa).*rxnsld(:,iz));
                drxndisp_tmp = sum(staq(:,ispa).*drxnsld_dmaq(:,ispa,iz));
                
                if (iz==1); caq_tmp_n = caqi_tmp; end
                    
                edif_tmp = 1d3*poro(iz)*sat(iz)*tora(iz)*d_tmp;
                edif_tmp_p = 1d3*poro(min(iz+1,nz))*sat(min(iz+1,nz))*tora(min(iz+1,nz))*d_tmp;
                edif_tmp_n = 1d3*poro(max(iz-1,1))*sat(max(iz-1,1))*tora(max(iz-1,1))*d_tmp;

                amx3(row,row) = ( ...
                    (poro(iz)*sat(iz)*1d3*1d0)/merge(1d0,dt,dt_norm)  ...
                    -(0.5d0*(edif_tmp +edif_tmp_p)*merge(0d0,-1d0,iz==nz)/( 0.5d0*(dz(iz)+dz(min(nz,iz+1))) ) ...
                    -0.5d0*(edif_tmp +edif_tmp_n)*(1d0)/( 0.5d0*(dz(iz)+dz(max(1,iz-1))) ))/dz(iz) ...
                    *merge(dt,1d0,dt_norm) ...
                    + poro(iz)*sat(iz)*1d3*v(iz)*(1d0)/dz(iz)*merge(dt,1d0,dt_norm) ...
                    -drxndisp_tmp*merge(dt,1d0,dt_norm) ...
                    - sum(staq_ext(:,ispa).*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) ...
                    ) ...
                    *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp*sw_red);

                ymx3(row) = ( ...
                    (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/merge(1d0,dt,dt_norm)  ...
                    -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) ...
                    *merge(dt,1d0,dt_norm) ...
                    + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*merge(dt,1d0,dt_norm) ...
                    - rxn_tmp*merge(dt,1d0,dt_norm) ...
                    - caqsupp_tmp*merge(dt,1d0,dt_norm) ...
                    - rxn_ext_tmp*merge(dt,1d0,dt_norm) ...
                    ) ...
                    *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)

                if (iz~=1)  
                    amx3(row,row-nsp3) = ( ...
                        -(-0.5d0*(edif_tmp +edif_tmp_n)*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) ...
                        *merge(dt,1d0,dt_norm) ...
                        + poro(iz)*sat(iz)*1d3*v(iz)*(-1d0)/dz(iz)*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *caq_tmp_n ...
                        *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)
                end 
                
                if (iz~=nz)  
                    amx3(row,row+nsp3) = ( ...
                        -(0.5d0*(edif_tmp +edif_tmp_p)*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz) ...
                        *merge(dt,1d0,dt_norm) ...
                        ) ...
                        *caq_tmp_p ...
                        *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)
                end 
                
                if (~sld_enforce)  
                    for isps = 1: nsp_sld
                        col = nsp3*(iz-1)+ isps;
                        
                        amx3(row, col) = (     ... 
                            - staq(isps,ispa)*drxnsld_dmsld(isps,iz)*merge(dt,1d0,dt_norm) ...
                            - sum(staq_ext(:,ispa).*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *msldx(isps,iz) ...
                            *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)
                    end 
                end  
                
                for ispa2 = 1: nsp_aq
                    col = nsp3*(iz-1)+ nsp_sld*solve_sld + ispa2;
                    
                    if (ispa2 == ispa); continue; end
                    
                    amx3(row,col) = amx3(row,col) + (     ... 
                        - sum(staq(:,ispa).*drxnsld_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) ...
                        - sum(staq_ext(:,ispa).*drxnext_dmaq(:,ispa2,iz))*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *maqx(ispa2,iz) ...
                        *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)
                end 
                
                for ispg = 1: nsp_gas
                    col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg;
                    
                    amx3(row,col) = amx3(row,col) + (     ... 
                        - sum(staq(:,ispa).*drxnsld_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) ...
                        - sum(staq_ext(:,ispa).*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *mgasx(ispg,iz) ...
                        *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp*sw_red);   % commented out (is this necessary?)
                end 
                        
                flx_aq(ispa,itflx,iz) = (...
                    (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  ...
                    ); 
                flx_aq(ispa,iadv,iz) = (...
                    + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) ...
                    );
                flx_aq(ispa,idif,iz) = (...
                    -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) ...
                    );
                flx_aq(ispa,irxn_sld(:),iz) = (... 
                    - staq(:,ispa).*rxnsld(:,iz) ...
                    );
                flx_aq(ispa,irain,iz) = (...
                    - caqsupp_tmp ...
                    );
                flx_aq(ispa,irxn_ext(:),iz) = (...
                    - staq_ext(:,ispa).*rxnext(:,iz) ...
                    );
                flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz));
                if (isnan(flx_aq(ispa,ires,iz)))  
                    warning( ['nan in flx_aq %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chraq(ispa),iz,flx_aq(ispa,1:nflx,iz) );
                end 
                
                amx3(row,:) = amx3(row,:)*fact;
                ymx3(row) = ymx3(row)*fact;
            
            end 
            
        end  % ==============================
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    pCO2 pO2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        for iz = 1: nz
            
            for ispg = 1: nsp_gas
            
                row = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg;
                
                pco2n_tmp = mgasx(ispg,max(1,iz-1));
                khco2n_tmp = khgasx(ispg,max(1,iz-1));
                edifn_tmp = dgas(ispg,max(1,iz-1));
                if (iz == 1)  
                    pco2n_tmp = mgasi(ispg);
                    khco2n_tmp = khgasi(ispg);
                    edifn_tmp = dgasi(ispg);
                end 

                amx3(row,row) = ( ...
                    (agasx(ispg,iz) + dagas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) ...
                    -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    +0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) ...
                    - 0.5d0*(ddgas_dmgas(ispg,ispg,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  ...
                    *merge(dt,1d0,dt_norm) ...
                    +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) ...
                    +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg,iz)*mgasx(ispg,iz))/dz(iz) *merge(dt,1d0,dt_norm) ...
                    -sum(stgas_ext(:,ispg).*drxnext_dmgas(:,ispg,iz))*merge(dt,1d0,dt_norm) ...
                    -drxngas_dmgas(ispg,ispg,iz)*merge(dt,1d0,dt_norm) ...
                    ) ...
                    *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                
                ymx3(row) = ( ...
                    (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/merge(1d0,dt,dt_norm) ...
                    -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                          /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                    *merge(dt,1d0,dt_norm) ...
                    +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*merge(dt,1d0,dt_norm) ...
                    -sum(stgas_ext(:,ispg).*rxnext(:,iz))*merge(dt,1d0,dt_norm) ...
                    -rxngas(ispg,iz)*merge(dt,1d0,dt_norm) ...
                    -mgassupp(ispg,iz)*merge(dt,1d0,dt_norm) ...
                    ) ...
                    *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                
                
                if (iz~=nz)  
                    amx3(row,row+nsp3) = ( ...
                            -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                            + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) ...
                                  /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                    
                    for ispa = 1:nsp_aq
                        col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa;
                        amx3(row,col+nsp3) = ( ...
                            -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) ...
                                  /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                                
                    end 
                
                end 
                
                if (iz~=1)  
                    amx3(row,row-nsp3) = ( ...
                        -(- 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) ...
                        - 0.5d0*(ddgas_dmgas(ispg,ispg,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) ...
                              /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*merge(dt,1d0,dt_norm)  ...
                        +poro(iz)*sat(iz)*v(iz)*1d3*(-khgasx(ispg,iz-1)*1d0)/dz(iz)*merge(dt,1d0,dt_norm) ...
                        +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmgas(ispg,ispg,iz-1)*mgasx(ispg,iz-1))/dz(iz)*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *merge(0.0d0,mgasx(ispg,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                    
                    for ispa = 1:nsp_aq
                        col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa;

                        amx3(row,col-nsp3) = ( ...
                            -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) ...
                                  /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*merge(dt,1d0,dt_norm)  ...
                            +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                               
                    end 
                end 
                
                if (~sld_enforce)  
                    for isps = 1:nsp_sld
                        col = nsp3*(iz-1) + isps;
                        amx3(row,col) = ( ...
                            -drxngas_dmsld(ispg,isps,iz)*merge(dt,1d0,dt_norm) ...
                            -sum(stgas_ext(:,ispg).*drxnext_dmsld(:,isps,iz))*merge(dt,1d0,dt_norm) ...
                            ) ...
                            *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                    end 
                end 
                
                for ispa = 1: nsp_aq
                    col = nsp3*(iz-1) + nsp_sld*solve_sld + ispa; 
                    amx3(row,col) = ( ...
                        (dagas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) ...
                        -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                              /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                        - 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)  ...
                        *merge(dt,1d0,dt_norm) ...
                        +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) ...
                        -drxngas_dmaq(ispg,ispa,iz)*merge(dt,1d0,dt_norm) ...
                        -sum(stgas_ext(:,ispg).*drxnext_dmaq(:,ispa,iz))*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                    
                end 
                
                for ispg2 = 1: nsp_gas
                    if (ispg == ispg2); continue; end
                    col = nsp3*(iz-1) + nsp_sld*solve_sld + nsp_aq + ispg2;
                    amx3(row,col) = ( ...
                        (dagas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/merge(1d0,dt,dt_norm) ...
                        -( 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                              /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                        - 0.5d0*(ddgas_dmgas(ispg,ispg2,iz))*(mgasx(ispg,iz)-pco2n_tmp) ...
                              /(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*merge(dt,1d0,dt_norm)  ...
                        +poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmgas(ispg,ispg2,iz)*mgasx(ispg,iz))/dz(iz)*merge(dt,1d0,dt_norm) ...
                        -drxngas_dmgas(ispg,ispg2,iz)*merge(dt,1d0,dt_norm) ...
                        -sum(stgas_ext(:,ispg).*drxnext_dmgas(:,ispg2,iz))*merge(dt,1d0,dt_norm) ...
                        ) ...
                        *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg)*sw_red);
                    
                end 
                
                flx_gas(ispg,itflx,iz) = ( ...
                    (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt ...
                    );         
                flx_gas(ispg,idif,iz) = ( ...
                    -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                          /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                    - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                    );
                flx_gas(ispg,iadv,iz) = ( ...
                    +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) ...
                    );
                flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg).*rxnext(:,iz);
                flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz);
                flx_gas(ispg,irxn_sld(:),iz) = ( ...
                    - stgas(:,ispg).*rxnsld(:,iz) ...
                    );
                flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz));
                
                if (any(isnan(flx_gas(ispg,:,iz)))) 
                    warning( ['nan in flx_gas %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chrgas(ispg),iz,flx_gas(ispg,1:nflx,iz) );
                end 
            end 

        end 
        
        fact2= max(abs(amx3),[],'all');
        
        amx3 = amx3/fact2;
        ymx3 = ymx3/fact2;
        
        ymx3=-1.0d0*ymx3;

        if (any(isnan(amx3),'all')||any(isnan(ymx3))||any(amx3>infinity,'all')||any(ymx3>infinity))  
            warning('error in mtx:\ninsnan(amx) = %s\ninsnan(ymx) =  %s\n' ...
                ,string(any(isnan(amx3),'all')),string(any(isnan(ymx3))) );

            if (any(isnan(ymx3)))  
                for ie = 1:nsp3*(nz)
                    if (isnan(ymx3(ie)))  
                        warning('NAN in ymx is here...%d\n',ie);
                    end
                end
            end


            if (any(isnan(amx3),'all'))  
                for ie = 1:nsp3*(nz)
                    for ie2 = 1:nsp3*(nz)
                        if (isnan(amx3(ie,ie2)))  
                            warning('NAN in amx is here...%d\t%d\n',ie,ie2);
                        end
                    end
                end
            end
            
            flgback = true;
            return
        end
        
        % dlmwrite('amx3_matlab.txt',amx3);
        % dlmwrite('ymx3_matlab.txt',ymx3);
        % error('printing stuff');
        % pause;
        
        % call DGESV(nsp3*(Nz),int(1),amx3,nsp3*(Nz),IPIV3,ymx3,nsp3*(Nz),INFO) 
        [xmx3,rmx3] = linsolve(amx3,ymx3);
        ymx3 = xmx3;

        if (any(isnan(ymx3))) 
            warning('error in soultion');
            flgback = true;
            return
        end

        for iz = 1: nz
            if (~sld_enforce)  
                for isps = 1: nsp_sld
                    row = isps + nsp3*(iz-1);

                    if (isnan(ymx3(row)))  
                        error('nan at %d\t (z = %7.6f) for species %s\n', iz,z(iz),chrsld(isps));
                    end
                    
                    emx3(row) = msldx(isps,iz)*exp(ymx3(row)) - msldx(isps,iz);

                    if ((~isnan(ymx3(row))) && ymx3(row) >threshold)  
                        msldx(isps,iz) = msldx(isps,iz)*corr;
                    elseif (ymx3(row) < -threshold)  
                        msldx(isps,iz) = msldx(isps,iz)/corr;
                    else   
                        msldx(isps,iz) = msldx(isps,iz)*exp(ymx3(row));
                    end
                    
                    if ( msldx(isps,iz)<msldth(isps))  % too small trancate value and not be accounted for error 
                        msldx(isps,iz)=msldth(isps);
                        ymx3(row) = 0d0;
                    end
                end 
            end 
            
            for ispa = 1: nsp_aq
                row = ispa + nsp_sld*solve_sld + nsp3*(iz-1);

                if (isnan(ymx3(row)))  
                    error('nan at %d\t (z = %7.6f) for species %s\n', iz,z(iz),chraq(ispa));
                end
                
                emx3(row) = poro(iz)*sat(iz)*1d3*maqx(ispa,iz)*exp(ymx3(row)) - poro(iz)*sat(iz)*1d3*maqx(ispa,iz);

                if ((~isnan(ymx3(row))) && ymx3(row) >threshold)  
                    maqx(ispa,iz) = maqx(ispa,iz)*corr;
                elseif (ymx3(row) < -threshold)  
                    maqx(ispa,iz) = maqx(ispa,iz)/corr;
                else   
                    maqx(ispa,iz) = maqx(ispa,iz)*exp(ymx3(row));
                end
                
                if (maqx(ispa,iz)<maqth(ispa))  % too small trancate value and not be accounted for error 
                    maqx(ispa,iz)=maqth(ispa);
                    ymx3(row) = 0d0;
                end
            end 
            
            for ispg = 1: nsp_gas
                row = ispg + nsp_aq + nsp_sld*solve_sld + nsp3*(iz-1);

                if (isnan(ymx3(row)))  
                    error('nan at %d\t (z = %7.6f) for species %s\n', iz,z(iz),chrgas(ispg));
                end
                
                emx3(row) =agasx(ispg,iz)* mgasx(ispg,iz)*exp(ymx3(row)) - agasx(ispg,iz)*mgasx(ispg,iz);

                if ((~isnan(ymx3(row))) && ymx3(row) >threshold)  
                    mgasx(ispg,iz) = mgasx(ispg,iz)*corr;
                elseif (ymx3(row) < -threshold)  
                    mgasx(ispg,iz) = mgasx(ispg,iz)/corr;
                else   
                    mgasx(ispg,iz) = mgasx(ispg,iz)*exp(ymx3(row));
                end
                
                if (mgasx(ispg,iz)<mgasth(ispg))  % too small trancate value and not be accounted for error 
                    mgasx(ispg,iz)=mgasth(ispg);
                    ymx3(row) = 0d0;
                end
            end 

        end 

        if (fact_tol == 1d0)  
            error = max(exp(abs(ymx3))) - 1.0d0;
        else 
            error = max((abs(emx3)));
        end 
        
        if (isnan(error)); error = 1d4; end

        if (isnan(error)|| any(isnan(msldx),'all') || any(isnan(maqx),'all')|| any(isnan(mgasx),'all'))  
            error = 1d3;
            warning('error is NaN; values are returned to those before iteration with reducing dt\n')
            warning('Nan in error?\t%s\nNan in msldx?\t%s\nNan in maqx?\t%s\nNan in mgasx?\t%s\n' ...
                ,string(isnan(error)),string(any(isnan(msldx),'all')) ...
                ,string(any(isnan(maqx),'all')),string(any(isnan(mgasx),'all')) );
            flgback = true;
            return
        end

        if (display)  
            fprintf("error in %d's iteration = %7.6E\t with time step = %7.6E\t[yr]\n",iter,error,dt);
        end      
        iter = iter + 1; 

        if (iter > iter_max ) 
            if (dt==0d0)  
                error('dt==0d0; stop');
            end 
            flgback = true;
            return
        end 
    end

    % just adding flx calculation at the end 
     
    flx_sld(:,:,:) = 0d0;
    flx_aq(:,:,:) = 0d0;
    flx_gas(:,:,:) = 0d0;

    flx_co2sp(:,:,:) = 0d0;

    % pH calculation and its derivative wrt aq and gas species

    [ ...
        dprodmaq_all,dprodmgas_all,dso4fdmaq_all,dso4fdmgas_all ...% output
        ,prox,ph_error,so4f,ph_iter ...% output
        ] = calc_pH_v7_3( ...
        nz,kw,nsp_aq,nsp_gas,nsp_aq_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ...% input 
        ,chraq,chraq_cnst,chraq_all,chrgas,chrgas_cnst,chrgas_all ...%input
        ,maqx,maqc,mgasx,mgasc,keqgas_h,keqaq_h,keqaq_c,keqaq_s,maqth_all,keqaq_no3,keqaq_nh3 ...% input
        ,print_cb,print_loc,z ...% input 
        ,prox,so4f ...% output
        ); 
        
    % recalculation of rate constants for mineral reactions

    if (kin_iter)  

        ksld(:,:) = 0d0;
            
        [mgasx_loc] = get_mgasx_all( ...
            nz,nsp_gas_all,nsp_gas,nsp_gas_cnst ...
            ,chrgas,chrgas_all,chrgas_cnst ...
            ,mgasx,mgasc ...
            );

        for isps =1:nsp_sld 
            [kin,dkin_dmsp] = sld_kin( ...
                nz,rg,tc,sec2yr,tempk_0,prox,kw,kho,mv(isps) ...% input
                ,nsp_gas_all,chrgas_all,mgasx_loc ...% input
                ,chrsld(isps),'pro' ...% input 
                ); 
            ksld(isps,:) = kin(:);
        end 

    end 
        
    % if kin const. is specified in input file 
    if (nsld_kinspc > 0)  
        for isps_kinspc=1:nsld_kinspc    
            if ( any( chrsld == chrsld_kinspc(isps_kinspc)))  
                switch (chrsld_kinspc(isps_kinspc))
                    case{'g1','g2','g3'} % for OMs, turn over year needs to be provided [yr]
                        ksld(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = ( ...                   
                            1d0/kin_sld_spc(isps_kinspc) ...
                            ); 
                        dksld_dpro(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                        dksld_dso4f(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                        dksld_dmaq(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                        dksld_dmgas(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                    otherwise % otherwise, usual rate constant [mol/m2/yr]
                        ksld(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = ( ...                            
                            kin_sld_spc(isps_kinspc) ...
                            ); 
                        dksld_dpro(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                        dksld_dso4f(find(chrsld==chrsld_kinspc(isps_kinspc)),:) = 0d0;
                        dksld_dmaq(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                        dksld_dmgas(find(chrsld==chrsld_kinspc(isps_kinspc)),:,:) = 0d0;
                end
            end 
        end 
    end 

    % saturation state calc. and their derivatives wrt aq and gas species

    omega(:,:) = 0d0;

    for isps =1: nsp_sld
        dummy(:) = 0d0;
        [ ...
            domega_dmaq_all,domega_dmgas_all,domega_dpro_loc,domega_dso4f_loc ...% output
            ,dummy,omega_error ...% output
            ] = calc_omega_v4( ...
            nz,nsp_aq,nsp_gas,nsp_aq_all,nsp_sld_all,nsp_gas_all,nsp_aq_cnst,nsp_gas_cnst ... 
            ,chraq,chraq_cnst,chraq_all,chrsld_all,chrgas,chrgas_cnst,chrgas_all ...
            ,maqx,maqc,mgasx,mgasc,mgasth_all,prox,so4f ...
            ,keqsld_all,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3 ...
            ,staq_all,stgas_all ...
            ,chrsld(isps) ...
            );
        omega(isps,:) = dummy(:);
    end 

    rxnsld(:,:) = 0d0;
        
    [ ...
        rxnsld,drxnsld_dmsld,drxnsld_dmaq,drxnsld_dmgas ...% output
        ] = sld_rxn( ...
        nz,nsp_sld,nsp_aq,nsp_gas,msld_seed,hr,poro,mv,ksld,omega,nonprec,msldx,dz ...% input 
        ,dksld_dmaq,domega_dmaq,dksld_dmgas,domega_dmgas,precstyle ...% input
        ,msld,msldth,dt,sat,maq,maqth,agas,mgas,mgasth,staq,stgas ...% input
        );

    % adding reactions that are not based on dis/prec of minerals
    rxnext(:,:) = 0d0;

    for irxn=1:nrxn_ext
        dummy(:) = 0d0;
        dummy2(:) = 0d0;
        [ ... 
            dummy,dummy2,rxnext_error ...% output
            ] = calc_rxn_ext_dev_2( ...
            nz,nrxn_ext_all,nsp_gas_all,nsp_aq_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst  ...%input
            ,chrrxn_ext_all,chrgas,chrgas_all,chrgas_cnst,chraq,chraq_all,chraq_cnst ...% input
            ,poro,sat,maqx,maqc,mgasx,mgasc,mgasth_all,maqth_all,krxn1_ext_all,krxn2_ext_all ...% input
            ,nsp_sld,nsp_sld_cnst,chrsld,chrsld_cnst,msldx,msldc,rho_grain,kw ...%input
            ,rg,tempk_0,tc ...%input
            ,nsp_sld_all,chrsld_all,msldth_all,mv_all,hr,prox,keqgas_h,keqaq_h,keqaq_c,keqaq_s,so4f ...% input
            ,chrrxn_ext(irxn),'pro' ...% input 
            );
        if (rxnext_error) 
            flgback = true;
            return 
        end 
        rxnext(irxn,:) = dummy(:);
    end 

    if (~sld_enforce) 
        for iz = 1: nz  %================================
            
            for isps = 1: nsp_sld
                
                k_tmp = ksld(isps,iz);
                mv_tmp = mv(isps);
                omega_tmp = omega(isps,iz);
                omega_tmp_th = omega_tmp*nonprec(isps,iz);
                m_tmp = msldx(isps,iz);
                mth_tmp = msldth(isps);
                mi_tmp = msldi(isps);
                mp_tmp = msldx(isps,min(nz,iz+1));
                msupp_tmp = msldsupp(isps,iz);
                rxn_ext_tmp = sum(stsld_ext(:,isps).*rxnext(:,iz));
                mprev_tmp = msld(isps,iz);
                w_tmp = w(iz);
                wp_tmp = w(min(nz,iz+1));
                sporo_tmp = 1d0-poro(iz);
                sporop_tmp = 1d0-poro(min(nz,iz+1));
                sporoprev_tmp = 1d0-poroprev(iz);
                mn_tmp = msldx(isps,max(1,iz-1));
                wn_tmp = w(max(1,iz-1));
                sporon_tmp = 1d0-poro(max(1,iz-1));
                
                if (iz==1)  
                    mn_tmp = 0d0;
                    wn_tmp = 0d0;
                    sporon_tmp = 0d0;
                end 
                
                if (iz==nz)  
                    mp_tmp = mi_tmp;
                    wp_tmp = w_btm;
                    sporop_tmp = 1d0- poroi;
                end 
                
                if (msldunit == 'blk')  
                    sporo_tmp = 1d0;
                    sporop_tmp = 1d0;
                    sporon_tmp = 1d0;
                    sporoprev_tmp = 1d0;
                end 
                
                for iiz = 1: nz
                    if (trans(iiz,iz,isps)==0d0); continue; end
                        
                    flx_sld(isps,idif,iz) = flx_sld(isps,idif,iz) + ( ...
                        - trans(iiz,iz,isps)*msldx(isps,iiz) * sporo(iiz) ...
                        );
                end
                
                flx_sld(isps,itflx,iz) = ( ...
                    (sporo_tmp*m_tmp - sporoprev_tmp*mprev_tmp)/dt ...
                    );
                flx_sld(isps,iadv,iz) = ( ...
                    - ( sporop_tmp*wp_tmp*mp_tmp - sporo_tmp*w_tmp* m_tmp)/dz(iz)  ...
                    );
                flx_sld(isps,irxn_sld(isps),iz) = ( ...
                    + rxnsld(isps,iz) ...
                    );
                flx_sld(isps,irain,iz) = (...
                    - msupp_tmp  ...
                    );
                flx_sld(isps,irxn_ext(:),iz) = (...
                        - stsld_ext(:,isps).*rxnext(:,iz)  ...
                        );
                flx_sld(isps,ires,iz) = sum(flx_sld(isps,:,iz));
                if (isnan(flx_sld(isps,ires,iz)))  
                    warning( ['nan in flx_sld %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chrsld(isps),iz,flx_sld(isps,1:nflx,iz) );
                end   
                
            end 
        end  %================================
    end 

    for iz = 1: nz
        
        for ispa = 1: nsp_aq
            
            d_tmp = daq(ispa);
            caq_tmp = maqx(ispa,iz);
            caq_tmp_prev = maq(ispa,iz);
            caq_tmp_p = maqx(ispa,min(nz,iz+1));
            caq_tmp_n = maqx(ispa,max(1,iz-1));
            caqth_tmp = maqth(ispa);
            caqi_tmp = maqi(ispa);
            caqsupp_tmp = maqsupp(ispa,iz);
            rxn_ext_tmp = sum(staq_ext(:,ispa).*rxnext(:,iz));
            rxn_tmp = sum(staq(:,ispa).*rxnsld(:,iz));
            
            if (iz==1); caq_tmp_n = caqi_tmp; end
                
            edif_tmp = 1d3*poro(iz)*sat(iz)*tora(iz)*d_tmp;
            edif_tmp_p = 1d3*poro(min(iz+1,nz))*sat(min(iz+1,nz))*tora(min(iz+1,nz))*d_tmp;
            edif_tmp_n = 1d3*poro(max(iz-1,1))*sat(max(iz-1,1))*tora(max(iz-1,1))*d_tmp;
                    
            flx_aq(ispa,itflx,iz) = (...
                (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  ...
                );
            flx_aq(ispa,iadv,iz) = (...
                + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) ...
                );
            flx_aq(ispa,idif,iz) = (...
                -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) ...
                );
            flx_aq(ispa,irxn_sld(:),iz) = (... 
                - staq(:,ispa).*rxnsld(:,iz) ...
                ); 
            flx_aq(ispa,irain,iz) = (...
                - caqsupp_tmp ...
                ); 
            flx_aq(ispa,irxn_ext(:),iz) = (...
                - staq_ext(:,ispa).*rxnext(:,iz) ...
                ); 
            flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz));
            if (isnan(flx_aq(ispa,ires,iz)))  
                warning( ['nan in flx_aq %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chraq(ispa),iz,flx_aq(ispa,1:nflx,iz) );
            end 
        
        end 
        
    end  % ==============================

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    pCO2 pO2   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    khgas(:,:) = 0d0;
    khgasx(:,:) = 0d0;
    % added
    if (new_gassol)  
        [ ...
            khgas_all,khgasx_all,dkhgas_dpro_all,dkhgas_dso4f_all,dkhgas_dmaq_all,dkhgas_dmgas_all ...%output
            ] = calc_khgas_all( ...
            nz,nsp_aq_all,nsp_gas_all,nsp_gas,nsp_aq,nsp_aq_cnst,nsp_gas_cnst ...
            ,chraq_all,chrgas_all,chraq_cnst,chrgas_cnst,chraq,chrgas ...
            ,maq,mgas,maqx,mgasx,maqc,mgasc ...
            ,keqgas_h,keqaq_h,keqaq_c,keqaq_s,keqaq_no3  ...
            ,pro,prox,so4fprev,so4f ...
            );
            
        for ispg=1:nsp_gas
            khgas(ispg,:)=khgas_all(find(chrgas_all==chrgas(ispg)),:);
            khgasx(ispg,:)=khgasx_all(find(chrgas_all==chrgas(ispg)),:);
        end 
    end 

    dgas(:,:) = 0d0;

    agas(:,:) = 0d0;
    agasx(:,:) = 0d0;

    rxngas(:,:) = 0d0;

    for ispg = 1: nsp_gas
        
        if (~new_gassol)  % to be removed?
            switch (chrgas(ispg))
                case('pco2')
                    khgas(ispg,:) = kco2*(1d0+k1./pro(:)' + k1*k2./pro(:)'./pro(:)'); % previous value; should not change through iterations 
                    khgasx(ispg,:) = kco2*(1d0+k1./prox(:)' + k1*k2./prox(:)'./prox(:)');
                case('po2')
                    khgas(ispg,:) = kho; % previous value; should not change through iterations 
                    khgasx(ispg,:) = kho;
                case('pnh3')
                    khgas(ispg,:) = knh3*(1d0+pro(:)'/k1nh3); % previous value; should not change through iterations 
                    khgasx(ispg,:) = knh3*(1d0+prox(:)'/k1nh3);
                case('pn2o')
                    khgas(ispg,:) = kn2o; % previous value; should not change through iterations 
                    khgasx(ispg,:) = kn2o;
            end
        end 
        
        dgas(ispg,:) = ucv*poro(:)'.*(1.0d0-sat(:)')*1d3.*torg(:)'*dgasg(ispg)+poro(:)'.*sat(:)'.*khgasx(ispg,:)*1d3.*tora(:)'*dgasa(ispg);
        dgasi(ispg) = ucv*1d3*dgasg(ispg); 
        
        agas(ispg,:)= ucv*poroprev(:)'.*(1.0d0-sat(:)')*1d3+poroprev(:)'.*sat(:)'.*khgas(ispg,:)*1d3;
        agasx(ispg,:)= ucv*poro(:)'.*(1.0d0-sat(:)')*1d3+poro(:)'.*sat(:)'.*khgasx(ispg,:)*1d3;
        
        for isps = 1: nsp_sld
            rxngas(ispg,:) =  rxngas(ispg,:) + (...
                + stgas(isps,ispg)*rxnsld(isps,:) ...
                );
        end 
    end 

    for iz = 1: nz
        
        for ispg = 1: nsp_gas
            
            pco2n_tmp = mgasx(ispg,max(1,iz-1));
            khco2n_tmp = khgasx(ispg,max(1,iz-1));
            edifn_tmp = dgas(ispg,max(1,iz-1));
            if (iz == 1)  
                pco2n_tmp = mgasi(ispg);
                khco2n_tmp = khgasi(ispg);
                edifn_tmp = dgasi(ispg);
            end 
            
            flx_gas(ispg,itflx,iz) = ( ...
                (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz))/dt ...
                );    
            flx_gas(ispg,idif,iz) = ( ...
                -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                      /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                );
            flx_gas(ispg,iadv,iz) = ( ...
                +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz) ...
                );
            flx_gas(ispg,irxn_ext(:),iz) = -stgas_ext(:,ispg).*rxnext(:,iz);
            flx_gas(ispg,irain,iz) = - mgassupp(ispg,iz);
            flx_gas(ispg,irxn_sld(:),iz) = ( ...
                - stgas(:,ispg).*rxnsld(:,iz) ...
                );
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz));
            
            if (any(isnan(flx_gas(ispg,:,iz)))) 
                warning( ['nan in flx_gas %s\t %d\t' repmat('%7.6E \t',1,nflx) '\n'], chrgas(ispg),iz,flx_gas(ispg,1:nflx,iz) );
            end 
        end 
        
        if (any(chrgas=='pco2'))  
            ispg = find(chrgas=='pco2');
            
            pco2n_tmp = mgasx(ispg,max(1,iz-1));
            proi_tmp = prox(max(1,iz-1));
            if (iz == 1)  
                pco2n_tmp = mgasi(ispg);
                proi_tmp = proi;
            end 
            
            % gaseous CO2
            
            edifn_tmp = ucv*poro(max(1,iz-1))*(1.0d0-sat(max(1,iz-1)))*1d3*torg(max(1,iz-1))*dgasg(ispg);
            if (iz==1); edifn_tmp = dgasi(ispg); end
            
            flx_co2sp(1,itflx,iz) = ( ...
                (ucv*poro(iz)*(1.0d0-sat(Iz))*1d3*mgasx(ispg,iz)-ucv*poroprev(iz)*(1.0d0-sat(Iz))*1d3*mgas(ispg,iz))/dt ...
                );  
            flx_co2sp(1,idif,iz) = ( ...
                -( 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) ...
                      +ucv*poro(min(nz,iz+1))*(1.0d0-sat(min(nz,iz+1)))*1d3*torg(min(nz,iz+1))*dgasg(ispg)) ...
                      *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                      /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                - 0.5d0*(ucv*poro(iz)*(1.0d0-sat(iz))*1d3*torg(Iz)*dgasg(ispg) + edifn_tmp) ...
                      *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                ); 
            flx_co2sp(1,irxn_ext(:),iz) = -stgas_ext(:,ispg).*rxnext(:,iz);
            flx_co2sp(1,irain,iz) = - mgassupp(ispg,iz);
            flx_co2sp(1,irxn_sld(:),iz) = ( ...
                - stgas(:,ispg).*rxnsld(:,iz) ...
                );
                
            % dissolved CO2
            
            edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*1d3*tora(max(1,iz-1))*dgasa(ispg);
            if (iz==1); edifn_tmp = 0d0; end
            
            flx_co2sp(2,itflx,iz) = ( ...
                (poro(iz)*sat(iz)*kco2*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*1d3*mgas(ispg,iz))/dt ...
                );
            flx_co2sp(2,idif,iz) = ( ...
                -( 0.5d0*(poro(iz)*sat(iz)*kco2*1d3*tora(iz)*dgasa(ispg) ...
                      +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*1d3*tora(min(nz,iz+1))*dgasa(ispg)) ...
                      *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                      /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                - 0.5d0*(poro(iz)*sat(iz)*kco2*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) ...
                      *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                );
            flx_co2sp(2,iadv,iz) = ( ...
                +poro(iz)*sat(iz)*v(iz)*1d3*(kco2*mgasx(ispg,iz)- kco2*pco2n_tmp)/dz(iz) ...
                );
                
            % HCO3-
            
            edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*k1/prox(max(1,iz-1))*1d3*tora(max(1,iz-1))*dgasa(ispg);
            if (iz==1); edifn_tmp = 0d0; end
            
            flx_co2sp(3,itflx,iz) = ( ...
                (poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*mgasx(ispg,iz)-poroprev(iz)*sat(iz)*kco2*k1/pro(iz)*1d3*mgas(ispg,iz))/dt ...
                );  
            flx_co2sp(3,idif,iz) = ( ...
                -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*tora(iz)*dgasa(ispg) ...
                      +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*k1/prox(min(nz,iz+1))*1d3*tora(min(nz,iz+1))*dgasa(ispg)) ...
                      *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                      /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                - 0.5d0*(poro(iz)*sat(iz)*kco2*k1/prox(iz)*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) ...
                      *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                );
            flx_co2sp(3,iadv,iz) = ( ...
                +poro(iz)*sat(iz)*v(iz)*1d3*( ...
                      kco2*k1/prox(iz)*mgasx(ispg,iz) ...
                      - kco2*k1/proi_tmp*pco2n_tmp)/dz(iz) ...
                );
                
            % CO32-
            
            edifn_tmp = poro(max(1,iz-1))*sat(max(1,iz-1))*kco2*k1*k2/prox(max(1,iz-1))^2d0*1d3*tora(max(1,iz-1))*dgasa(ispg);
            if (iz==1); edifn_tmp = 0d0; end
            
            flx_co2sp(4,itflx,iz) = (  ...
                (poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)^2d0*1d3*mgasx(ispg,iz) ...
                      -poroprev(iz)*sat(iz)*kco2*k1*k2/pro(iz)^2d0*1d3*mgas(ispg,iz))/dt ...
                );  
            flx_co2sp(4,idif,iz) = ( ...
                -( 0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)^2d0*1d3*tora(iz)*dgasa(ispg) ...
                      +poro(min(nz,iz+1))*sat(min(nz,iz+1))*kco2*k1*k2/prox(min(nz,iz+1))^2d0*1d3*tora(min(nz,iz+1))*dgasa(ispg)) ...
                      *(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) ...
                      /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) ...
                - 0.5d0*(poro(iz)*sat(iz)*kco2*k1*k2/prox(iz)^2d0*1d3*tora(iz)*dgasa(ispg) + edifn_tmp) ...
                      *(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)  ...
                ); 
            flx_co2sp(4,iadv,iz) = ( ...
                +poro(iz)*sat(iz)*v(iz)*1d3*( ...
                      kco2*k1*k2/prox(iz)^2d0*mgasx(ispg,iz) ...
                      - kco2*k1*k2/proi_tmp^2d0*pco2n_tmp)/dz(iz) ...
                );
                
                
            flx_co2sp(1,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz),'all');
            flx_co2sp(2,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz),'all');
            flx_co2sp(3,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz),'all');
            flx_co2sp(4,ires,iz) = sum(flx_co2sp(:,1:nflx-1,iz),'all');
        end 
        
    end 

    if (chkflx && dt > dt_th)  
        flx_max_max = 0d0;
        for isps = 1: nsp_sld

            flx_max = 0d0;
            for iflx = 1: nflx
                flx_max = max(flx_max,abs(sum(squeeze(flx_sld(isps,iflx,:)).*dz(:))));
            end 
            
            flx_max_max = max(flx_max_max,flx_max);
        end 

        for ispa = 1: nsp_aq

            flx_max = 0d0;
            for iflx = 1: nflx
                flx_max = max(flx_max,abs(sum(squeeze(flx_aq(ispa,iflx,:)).*dz(:))));
            end 
            flx_max_max = max(flx_max_max,flx_max);
        end 

        for ispg = 1: nsp_gas

            flx_max = 0d0;
            for iflx = 1: nflx
                flx_max = max(flx_max,abs(sum(squeeze(flx_gas(ispg,iflx,:)).*dz(:))));
            end 
            flx_max_max = max(flx_max_max,flx_max);
        end 
        
        if (~sld_enforce)  
            for isps = 1: nsp_sld

                flx_max = 0d0;
                for iflx = 1: nflx
                    flx_max = max(flx_max,abs(sum(squeeze(flx_sld(isps,iflx,:)).*dz(:))));
                end 
                
                if (flx_max/flx_max_max > flx_max_tol &&  abs(sum(squeeze(flx_sld(isps,ires,:)).*dz(:)))/flx_max > flx_tol )  
                    warning('too large error in mass balance of sld phases\n');
                    warning('sp = %s\t| flx that raised the flag = %7.6E\t| flx_max = %7.6E\t| flx_tol =%7.6E\n' ...
                        ,chrsld(isps),abs(sum(squeeze(flx_sld(isps,ires,:)).*dz(:))),flx_max,flx_tol );
                    flgback = true;
                    return
                end 
            end 
        end 

        for ispa = 1: nsp_aq

            flx_max = 0d0;
            for iflx = 1: nflx
                flx_max = max(flx_max,abs(sum(squeeze(flx_aq(ispa,iflx,:)).*dz(:))));
            end 
            
            if (flx_max/flx_max_max > flx_max_tol  && abs(sum(squeeze(flx_aq(ispa,ires,:)).*dz(:)))/flx_max > flx_tol )  
                warning('too large error in mass balance of aq phases\n');
                warning('sp = %s\t| flx that raised the flag = %7.6E\t| flx_max = %7.6E\t| flx_tol =%7.6E\n' ...
                    ,chraq(ispa),abs(sum(squeeze(flx_aq(ispa,ires,:)).*dz(:))),flx_max,flx_tol );
                flgback = true;
                return
            end 
        end 

        for ispg = 1: nsp_gas

            flx_max = 0d0;
            for iflx = 1: nflx
                flx_max = max(flx_max,abs(sum(squeeze(flx_gas(ispg,iflx,:)).*dz(:))));
            end 
            
            if (flx_max/flx_max_max > flx_max_tol  && abs(sum(squeeze(flx_gas(ispg,ires,:)).*dz(:)))/flx_max > flx_tol )  
                warning('too large error in mass balance of gas phases\n');
                warning('sp = %s\t| flx that raised the flag = %7.6E\t| flx_max = %7.6E\t| flx_tol =%7.6E\n' ...
                    ,chrgas(ispg),abs(sum(squeeze(flx_gas(ispg,ires,:)).*dz(:))),flx_max,flx_tol );
                flgback = true;
                return
            end 
        end 
    end 

end


function [poro] = calc_poro( ...
    nz,nsp_sld,nflx,idif,irain ...% in
    ,flx_sld,mv,poroprev,w,poroi,w_btm,dz,tol,dt ...% in
    ,poro ...% inout
    )
    
    % local 
    DV=zeros(nz,1,'double');resi_poro=zeros(nz,1,'double');
    iz=0;isps=0;row=0;ie=0;ie2=0;
    w_tmp=0;wp_tmp=0;sporo_tmp=0;sporop_tmp=0;sporoprev_tmp=0;

    infinity = Inf;

    amx3 =zeros(nz,nz,'double');
    ymx3=zeros(nz,1,'double');
    xmx3=zeros(nz,1,'double');
    rmx = 0;
    % 
    % attempt to solve porosity under any kind of porosity - uplift rate relationship 
    % based on equation:
    % d(1-poro)/dt = d(1-poro)*w/dz - mv*1d-6*sum( flx_sld(mixing, dust, rxns) ) 

        
    ymx3(:) = 0d0;
    xmx3(:) = 0d0;
    amx3(:,:) = 0d0;
    DV(:) = 0d0;

    for iz=1:nz
        for isps = 1:nsp_sld 
            DV(iz) = DV(iz) + ( flx_sld(isps, 4 + isps,iz) + flx_sld(isps, idif ,iz) + flx_sld(isps, irain ,iz) ) ...
                *mv(isps)*1d-6 ;
        end 
        
        row = iz;
        
        w_tmp = w(iz);
        wp_tmp = w(min(nz,iz+1));
        sporo_tmp = 1d0-poro(iz);
        sporop_tmp = 1d0-poro(min(nz,iz+1)) ;
        sporoprev_tmp = 1d0-poroprev(iz);
                
        if (iz==nz)  
            wp_tmp = w_btm;
            sporop_tmp = 1d0 - poroi;
        end 
        
        if (iz~=nz)  
        
            ymx3(row) = ( ...
                + (1d0 - sporoprev_tmp)/dt    ...
                - ( 1d0*wp_tmp - 1d0*w_tmp)/dz(iz)  ...
                + DV(iz) ...
                );
                
            amx3(row,row) = ( ...
                + (-1d0 )/dt    ...
                - ( - (-1d0)*w_tmp)/dz(iz)  ...
                );
                
            amx3(row,row+1) = ( ...
                - ( -1d0*wp_tmp )/dz(iz)  ...
                );
            
        else 
        
            ymx3(row) = ( ...
                + (1d0 - sporoprev_tmp)/dt    ...
                - ( sporop_tmp*wp_tmp - 1d0*w_tmp)/dz(iz)  ...
                + DV(iz) ...
                );
                
            amx3(row,row) = ( ...
                + (-1d0 )/dt    ...
                - (- (-1d0)*w_tmp)/dz(iz)  ...
                );
        
        
        end 
        
    end 
        
    ymx3=-1.0d0*ymx3;

    if (any(isnan(amx3),'all') || any(isnan(ymx3)) || any(amx3>infinity,'all') || any(ymx3>infinity))  
        warning('porocalc: error in mtx');

        if (any(isnan(ymx3)))  
            for iz = 1: nz
                if (isnan(ymx3(iz)))  
                    warning('porocalc: NAN is here...%d\n',iz);
                end
            end 
        end


        if (any(isnan(amx3),'all'))  
            for ie = 1:(nz)
                for ie2 = 1:(nz)
                    if (isnan(amx3(ie,ie2)))  
                        warning('porocalc: NAN is here...%d\t%d\n',ie,ie2);
                    end
                end
            end
        end 
        error('stop');
        
    end

    % call DGESV(Nz,int(1),amx3,Nz,IPIV3,ymx3,Nz,INFO) 
    [xmx3,rmx3] = linsolve(amx3,ymx3);
    ymx3 = xmx3;

    poro = ymx3;
   
end


function [dpsd,psd_error_flg] = psd_diss( ...
    nz,nps ...% in
    ,z,DV,dt,pi,tol,poro ...% in 
    ,incld_rough,rough_c0,rough_c1 ...% in
    ,psd,ps,dps,ps_min,ps_max ...% in 
    ,chrsp ...% in 
    ,dpsd,psd_error_flg ...% inout
    )
    
    % local 
    dVd=zeros(nps,nz,'double');psd_old=zeros(nps,nz,'double');psd_new=zeros(nps,nz,'double');dpsd_tmp=zeros(nps,nz,'double');
    psd_tmp=zeros(nps,1,'double');dvd_tmp=zeros(nps,1,'double');
    ps_new=0;ps_newp=0;dvd_res=0;
    ips=0;iips=0;ips_new=0;iz=0;isps=0;
    safe_mode = false;
    % safe_mode = true;

    % attempt to do psd ( defined with particle number / bulk m3 / log (r) )
    % assumptions: 
    % 1. particle numbers are only affected by transport (including raining/dusting) 
    % 2. dissolution does not change particle numbers: it only affect particle distribution 
    % unless particle is the minimum radius. in this case particle can be lost via dissolution 
    % 3. when a mineral precipitates, it is assumed to increase particle radius?
    % e.g., when a 1 um of particle is dissolved by X m3, its radius is changed and this particle is put into a different bin of (smaller) radius 

    % sum of volume change of minerals at iz is DV = sum(flx_sld(5:5+nsp_sld,iz)*mv(:)*1d-6)*dt (m3 / m3) 
    % this must be distributed to different particle size bins (dV(r)) in proportion to psd * (4*pi*r^2)
    % dV(r) = DV/(psd*4*pi*r^2) where dV is m3 / bulk m3 / log(r) 
    % has to modify so that sum( dV * dps ) = DV 
    % new psd is obtained by dV(r) = psd(r)*( 4/3 * pi * r^3 - 4/3 * pi * r'^3 ) where r' is the new radius as a result of dissolution
    % if r' is exactly one of ps value (ps(ips) == r'), then  psd(r') = psd(r)
    % else: 
    % first find the closest r* value which is one of ps values.
    % then DV(r) = psd(r)* 4/3 * pi * r^3 - psd(r*)* 4/3 * pi * r*^3 
    % i.e., psd(r*) = [ psd(r)* 4/3 * pi * r^3 - DV(r)]  /( 4/3 * pi * r*^3)
    %               = [ psd(r)* 4/3 * pi * r^3 - psd(r)*( 4/3 * pi * r^3 - 4/3 * pi * r'^3 ) ] /( 4/3 * pi * r*^3)
    %               = psd(r) * 4/3 * pi * r'^3 /( 4/3 * pi * r*^3)
    %               = psd(r) * (r'/r*)^3
    % in this way volume is conservative? 
    % check: sum( psd(r) * 4/3 * pi * r^3 * dps) - sum( psd(r') * 4/3 * pi * r'^3 * dps) = DV 

    dpsd_tmp(:,:) = 0d0;
    psd_old = psd;
    psd_new = psd;
    for iz=1:nz
        
        % correct one?
        dVd(:,:) = 0d0;
        if (~incld_rough)  
            dVd(:,iz) = ( psd (:,iz) .* (10d0.^ps(:)).^2d0 );
        else
            dVd(:,iz) = ( psd (:,iz) .* (10d0.^ps(:)).^2d0 *rough_c0.*(10d0.^ps(:)).^rough_c1);
        end 
        
        if (all(dVd == 0d0))  
            fprintf('%s\t%s\n', chrsp,'all dissolved loc1?');
            psd_error_flg = true;
            return
        end 
        
        % scale with DV
        dVd(:,iz) = dVd(:,iz)*DV(iz)/sum(dVd(:,iz) .* dps(:));
        
        if ( abs( (sum(dVd(:,iz) .* dps(:)) - DV(iz))/DV(iz)) > tol ) 
            error('%s\t%s\t%7.6e\n', chrsp,' vol. balance failed somehow ',abs( (sum(dVd(:,iz) * dps(:)) - DV(iz))/DV(iz)));
            error('%d\t%7.6e\t%7.6e\n', iz, sum(dVd(:,iz) .* dps(:)), DV(iz));
        end 
        
        if (any(isnan(dVd(:,iz))) )  
            error('%s\t%s\n',chrsp, 'nan in dVd loc 1');
        end 
        
        for ips = 1: nps
            
            if ( psd(ips,iz) == 0d0); continue; end 
            
            if ( ips == 1 && dVd(ips,iz) > 0d0 )  
                % this is the minimum size dealt within the model 
                % so if dissolved (dVd > 0), particle number must reduce 
                % (revised particle volumes) = (initial particle volumes) - (volume change) 
                % psd'(ips,iz) * 4d0/3d0*pi*(10d0^ps(ips))^3d0 =  psd(ips,iz) * 4d0/3d0*pi*(10d0^ps(ips))^3d0 - dVd(ips,iz) 
                % [ psd'(ips,iz) - psd(ips,iz) ] * 4d0/3d0*pi*(10d0^ps(ips))^3d0 = - dVd(ips,iz) 
                if ( ~ safe_mode )   % for not care about producing negative particle number
                    dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0); 
                elseif ( safe_mode )    % never producing negative particle number
                    if ( dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0) < psd(ips,iz) )  % when dissolution does not consume existing particles 
                        dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0); 
                    else % when dissolution exceeds potential consumption of existing particles 
                        % dvd_res is defined as residual volume to be dissolved 
                        dvd_res = dVd(ips,iz) - psd(ips,iz)*(4d0/3d0*pi*(10d0^ps(ips))^3d0);  % residual 
                        
                        % distributing the volume to whole radius 
                        
                        % correct one?
                        dVd_tmp(:) = 0d0;
                        if (~incld_rough)  
                            dVd_tmp(ips+1:end) = ( psd (ips+1:end,iz) .* (10d0.^ps(ips+1:end)).^2d0 );
                        else
                            dVd_tmp(ips+1:end) = ( psd (ips+1:end,iz) .* (10d0.^ps(ips+1:end)).^2d0 *rough_c0.*(10d0.^ps(ips+1:end)).^rough_c1 );
                        end 
                        
                        if (all(dVd_tmp == 0d0))  
                            fprintf('%s\t%s\t%d\n',chrsp,'all dissolved loc2?',ips);
                            psd_error_flg = true;
                            return
                        end 
                        
                        % scale with dvd_res*dps
                        dVd_tmp(ips+1:end) = dVd_tmp(ips+1:end)*dvd_res*dps(ips)./sum(dVd_tmp(ips+1:end) .* dps(ips+1:end));
                        
                        % dVd(ips+1:,iz) = dVd(ips+1:,iz) + dvd_res/(nps - ips)
                        dVd(ips+1:end,iz) = dVd(ips+1:end,iz) + dVd_tmp(ips+1:end);
                        dVd(ips,iz) = dVd(ips,iz) - dvd_res;
                        
                        if ( abs( (sum(dVd(:,iz) .* dps(:)) - DV(iz))/DV(iz)) > tol ) 
                            error('%s\t%s\t%7.6e\n',chrsp, ' vol. balance failed somehow loc2 ',abs( (sum(dVd(:,iz) .* dps(:)) - DV(iz))/DV(iz)) );
                            error('%d\t%7.6e\t%7.6e\n', iz, sum(dVd(:,iz) * dps(:)), DV(iz) );
                        end 
                        
                        if (any(isnan(dVd(:,iz))) )  
                            error('%s\t%s\n', chrsp,'nan in dVd loc 2');
                        end 
                        
                        % if (dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0) > psd(ips,iz))  
                            % print *, 'error: stop',psd(ips,iz),dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0)
                            % stop
                        % end 
                        
                        % dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0) 
                        dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - psd(ips,iz); 
                    end 
                end 
            
            elseif ( ips == nps && dVd(ips,iz) < 0d0 )  
                % this is the max size dealt within the model 
                % so if precipirated (dVd < 0), particle number must increase  
                % (revised particle volumes) = (initial particle volumes) - (volume change) 
                % psd'(ips,iz) * 4d0/3d0*pi*(10d0^ps(ips))^3d0 =  psd(ips,iz) * 4d0/3d0*pi*(10d0^ps(ips))^3d0 - dVd(ips,iz) 
                % [ psd'(ips,iz) - psd(ips,iz) ] * 4d0/3d0*pi*(10d0^ps(ips))^3d0 = - dVd(ips,iz) 
                dpsd_tmp(ips,iz) = dpsd_tmp(ips,iz) - dVd(ips,iz)/(4d0/3d0*pi*(10d0^ps(ips))^3d0);
            
            else 
                % new r*^3 after dissolution/precipitation 
                ps_new =  ( 4d0/3d0*pi*(10d0^ps(ips))^3d0 - dVd(ips,iz) /psd(ips,iz) )/(4d0/3d0*pi);
                
                if (ps_new <= 0d0)   % too much dissolution so removing all particles cannot explain dVd at a given bin ips
                    ps_new = ps_min;
                    ps_newp = 10d0^ps(1);
                    dpsd_tmp(1,iz) =  dpsd_tmp(1,iz) + psd(ips,iz);
                    dpsd_tmp(ips,iz) =  dpsd_tmp(ips,iz) - psd(ips,iz);
                    
                    dvd_res = dVd(ips,iz) - psd(ips,iz)*(4d0/3d0*pi*(10d0^ps(ips))^3d0);  % residual
                    
                    % distributing the volume to whole radius 
                    
                    % correct one?
                    dVd_tmp(:) = 0d0;
                    if (~incld_rough)  
                        % dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0^ps(ips+1:))^2d0 )
                        for iips=1:nps
                            if (iips == ips); continue; end
                            dVd_tmp(iips) = ( psd (iips,iz) * (10d0^ps(iips))^2d0 );
                        end 
                    else
                        % dVd_tmp(ips+1:) = ( psd (ips+1:,iz) * (10d0^ps(ips+1:))^2d0 * rough_c0*(10d0^ps(ips+1:))^rough_c1 )
                        for iips = 1:nps
                            if (iips == ips); continue; end
                            dVd_tmp(iips) = ( psd (iips,iz) * (10d0^ps(iips))^2d0 * rough_c0*(10d0^ps(iips))^rough_c1 );
                        end 
                    end 
                    
                    if (all(dVd_tmp == 0d0))  
                        fprintf('%s\t%s\t%d\n',chrsp,'all dissolved loc3?',ips);
                        psd_error_flg = true;
                        return
                    end 
                    
                    % dVd_tmp(ips+1:) = dVd_tmp(ips+1:)*dvd_res*dps(ips)/sum(dVd_tmp(ips+1:) * dps(ips+1:))
                    dVd_tmp(:) = dVd_tmp(:)*dvd_res*dps(ips)/sum(dVd_tmp(:) .* dps(:));
                    dVd_tmp(ips) = - dvd_res;
                    
                    % dVd(ips+1:,iz) = dVd(ips+1:,iz) + dVd_tmp(ips+1:)
                    % dVd(ips,iz) = dVd(ips,iz) - dvd_res
                    dVd(:,iz) = dVd(:,iz) + dVd_tmp(:);
                    
                    if ( abs( (sum(dVd(:,iz) .* dps(:)) - DV(iz))/DV(iz)) > tol ) 
                        error('%s\t%s\t%7.6e\n',chrsp, ' vol. balance failed somehow loc4 ',abs( (sum(dVd(:,iz) .* dps(:)) - DV(iz))/DV(iz)) );
                        error('going to stop\n');
                        error('%d\t%7.6e\t%7.6e\n', iz, sum(dVd(:,iz) * dps(:)), DV(iz));
                    end 
                    
                    if (any(isnan(dVd(:,iz))) )  
                        error('%s\t%s\n',chrsp, 'nan in dVd loc 4');
                    end 
                    
                else 
                    % new r*
                    ps_new =  ps_new^(1d0/3d0); 
                    if (ps_new <= ps_min)  
                        ips_new = 1;
                    elseif (ps_new >= ps_max)  
                        ips_new = nps;
                    else 
                        for iips = 1: nps -1
                            if ( ( ps_new - 10d0^ps(iips) ) *  ( ps_new - 10d0^ps(iips+1) ) <= 0d0 )  
                                if ( log10(ps_new) <= 0.5d0*( ps(iips) + ps(iips+1) ) )  
                                    ips_new = iips;
                                else 
                                    ips_new = iips + 1; 
                                end 
                                break 
                            end 
                        end 
                    end 
                    ps_newp = 10d0^ps(ips_new);  % closest binned particle radius to r*
                    dpsd_tmp(ips_new,iz) = dpsd_tmp(ips_new,iz) + psd(ips,iz)*(ps_new/ps_newp)^3d0;
                    dpsd_tmp(ips,iz) =  dpsd_tmp(ips,iz) - psd(ips,iz);
                end 
            
            end 
        end 
    end 

    if (any(isnan(psd),'all'))  
        error('%s\t%s\n',chrsp, 'nan in psd');
    end 
    if (any(psd<0d0,'all'))  
        error('%s\t%s\n' ,chrsp, 'negative psd');
    end 
    if (any(isnan(dVd),'all'))  
        fprintf('%s\t%s\n', chrsp, 'nan in dVd'); 
        for iz = 1: nz
            for ips=1:nps
                if (isnan(dVd(ips,iz)))  
                    fprintf('%s\t%s\t%d\t%d\t%7.6e\t%7.6e\n',chrsp, 'ips,iz,dVd,psd',ips,iz,dVd(ips,iz),psd(ips,iz));
                end 
            end 
        end 
        error('stop');
    end 

    psd_new(:,:) = psd_new(:,:) + dpsd_tmp(:,:);
    for iz = 1: nz
        % if ( abs(DV(iz)) > tol  ...
            % && abs ( ( sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:)) ...
            % - sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:)) - DV(iz) ) / DV(iz) ) > tol )   
        if ( abs(DV(iz)) > tol  ...
            && abs ( ( - sum( dpsd_tmp(:,iz) * 4d0/3d0 * pi .* (10d0.^ps(:)).^3d0 .* dps(:)) ...
             - DV(iz) ) / DV(iz) ) > tol )   
            fprintf( '%s\t%s\t%7.6e\n', chrsp,'checking the vol. balance and failed ... ' ...
                , abs ( ( - sum( dpsd_tmp(:,iz) * 4d0/3d0 * pi .* (10d0.^ps(:)).^3d0 .* dps(:)) ...
                 - DV(iz) ) / DV(iz) ) );
            fprintf( '%s\t%7.6e\t%7.6e\t%7.6e\t%7.6e\n', iz, sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:)) ...
                ,sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:)) ... 
                ,sum( psd_new(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:)) ...
                  -sum( psd_old(:,iz) * 4d0/3d0 * pi * (10d0^ps(:))^3d0 * dps(:))  ...
                ,DV(iz) );
            psd_error_flg = true;
        end 
    end 

    if (any(isnan(dpsd_tmp),'all'))  
        fprintf('%s\t%s\n',chrsp, 'nan in dpsd _rxn' );
        for iz = 1: nz
            for ips=1:nps
                if (isnan(dpsd_tmp(ips,iz)))  
                    fprintf('%s\t%d\t%d\t%7.6e\t%7.6e\n', 'ips,iz,dpsd_tmp,psd',ips,iz,dpsd_tmp(ips,iz),psd(ips,iz) );
                end 
            end 
        end 
        error('stop');
    end 

    dpsd(:,:) = dpsd(:,:) + dpsd_tmp(:,:);
   
end 


function [ ...
    flgback,flx_max_max ...% inout
    ,psdx,flx_psd ...% out
    ] = psd_implicit_all_v2( ...
    nz,nsp_sld,nps,nflx_psd ...% in
    ,z,dz,dt,pi,tol,w0,w,poro,poroi,poroprev ...% in 
    ,incld_rough,rough_c0,rough_c1 ...% in
    ,trans ...% in
    ,psd,psd_pr,ps,dps,dpsd,psd_rain ...% in   
    ,chrsp ...% in 
    ,flgback,flx_max_max ...% inout
    )

    % output 
    psdx=zeros(nps,nz,'double');
    flx_psd=zeros(nps,nflx_psd,nz,'double'); % itflx,iadv,idif,irain,irxn,ires
    % local 
    psd_old=zeros(nps,nz,'double');dpsd_tmp=zeros(nps,nz,'double');
    DV=zeros(nz,1,'double');kpsd=zeros(nz,1,'double');sporo=zeros(nz,1,'double');
    fact_tol=zeros(nps,1,'double'); 
    iz=0;isps=0;ips=0;iiz=0;row=0;col=0;ie=0;ie2=0;iips=0;
    iter=0;iflx=0;
    error=0;fact=0;flx_max=0; % ,flx_max_max
    vol=0;surf=0;m_tmp=0;mp_tmp=0;mi_tmp=0;mprev_tmp=0;rxn_tmp=0;drxn_tmp=0;w_tmp=0;wp_tmp=0;trans_tmp=0;
    msupp_tmp=0;sporo_tmp=0; sporop_tmp=0;sporoprev_tmp=0;dtinv=0;dzinv=0;

    chrflx_psd=strings(nflx_psd,1);

    dt_norm = true;
    infinity = Inf;
    threshold = 20d0;
    % threshold = 3d0;
    % corr = 1.5d0;
    corr = exp(threshold);
    iter_max = 50;
    % nflx_psd = 6;
    flx_tol = 1d-3;
    % flx_tol = 1d-4;
    flx_max_tol = 1d-6;
    % flx_max_tol = 1d-5;
    dt_th = 1d-6;
    [itflx_psd,iadv_psd,idif_psd,irain_psd,irxn_psd,ires_psd]=deal(1,2,3,4,5,6);
    % chkflx = false;
    chkflx = true;

    amx3=zeros(nz,nz,'double');ymx3=zeros(nz,1,'double');emx3=zeros(nps,1,'double');emx3_loc=zeros(nz,1,'double');
    xmx3=zeros(nz,1,'double');rmx3=0;
            
    % for iz=1:nz
        % hr(iz) = sum( 4d0*pi*(10d0^ps(:))^2d0*psd(:,iz)*dps(:) )
    % end 
    % for iz=1:nz
        % hr(iz) = sum( 4d0*pi*(10d0^ps(:))^2d0*rough_c0*(10d0^ps(:))^rough_c1*psd(:,iz)*dps(:) )
    % end 
    % 
    % attempt to for psd ( defined with particle number / bulk m3 / log (r) )
    % solve as for msld 
    % one of particle size equations are used to give massbalance constraint  


    chrflx_psd = string({'tflx';'adv';'dif';'rain';'rxn';'res'});

    dtinv = 1d0/dt;

    sporo(:) = 1d0 - poro(:);
    sporo(:) = 1d0; 

    psdx = psd;

    kpsd(:) = 0d0;
    % for iz = 1: nz
        % for isps = 1:nsp_sld 
            % DV(iz) = DV(iz) + flx_sld(isps, 4 + isps,iz)*mv(isps)*1d-6
        % end 
        % kpsd(iz) = DV(iz) / hr(iz) %/ sum( msldx(:,iz) * mv(:) * 1d-6 )
    % end 


    error = 1d4; 
    iter = 1;
    emx3(:) = error;

    for ips =1:nps
        fact_tol(ips) = max(psd(ips,:)) * 1d-12;
        % fact_tol(ips) = max(psd(ips,:)) * 1d-15;
    end 
    % fact_tol(:) = 1d0;

    % while (error > tol*fact_tol) 
    while (error > 1d0) 
    % while ( any (emx3 > fact_tol )  ) 
        
        
        for ips = 1: nps
        
            if (emx3(ips) <= fact_tol(ips)); continue; end

            amx3(:,:) = 0d0;
            ymx3(:) = 0d0;
            xmx3(:) = 0d0;
            rmx3 = 0d0;
            
            for iz = 1: nz
            
                row =  iz;
                
                vol  = 4d0/3d0*pi*(10d0^ps(ips))^3d0;
                surf = 4d0*pi*(10d0^ps(ips))^2d0;
                        
                m_tmp = vol * psdx(ips,iz) * dps(ips);
                mprev_tmp = vol * psd(ips,iz) * dps(ips);        
                % rxn_tmp = vol * psdx(ips,iz)*dps(ips) ...
                    % * surf * psdx(ips,iz)*dps(ips) * kpsd(iz); 
                % rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz); 
                % rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt ;
                rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) * dtinv;
                
                % msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt;
                msupp_tmp = vol * psd_rain(ips,iz) * dps(ips) * dtinv;
                
                mi_tmp = vol * psd_pr(ips) * dps(ips);
                mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips);
                
                dzinv = 1d0/dz(iz);
                
                % drxn_tmp = ... 
                    % vol * 1d0 * dps(ips) ...
                    % * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) ...
                    % + vol * psdx(ips,iz) * dps(ips) ...
                    % * su;
                % drxn_tmp = surf * 1d0 * dps(ips) * kpsd(iz) ;
                drxn_tmp = 0d0; 
                
                w_tmp = w(iz);
                wp_tmp = w(min(nz,iz+1));

                sporo_tmp = 1d0-poro(iz);
                sporop_tmp = 1d0-poro(min(nz,iz+1));
                sporoprev_tmp = 1d0-poroprev(iz);
                
                if (iz==nz)  
                    mp_tmp = mi_tmp;
                    wp_tmp = w0;
                    sporop_tmp = 1d0- poroi;
                end 
                
                sporo_tmp = 1d0;
                sporop_tmp = 1d0;
                sporoprev_tmp = 1d0;

                amx3(row,row) = ( ...
                    sporo_tmp * vol * dps(iz) * 1d0  *  merge(1d0,dtinv,dt_norm)     ...
                    + sporo_tmp * vol * dps(iz) * w_tmp * dzinv  *merge(dt,1d0,dt_norm)    ...
                    + sporo_tmp * drxn_tmp * merge(dt,1d0,dt_norm) ...
                    ) ...
                    * psdx(ips,iz);

                ymx3(row) = ( ...
                    ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * merge(1d0,dtinv,dt_norm) ...
                    -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv * merge(dt,1d0,dt_norm)  ...
                    + sporo_tmp* rxn_tmp * merge(dt,1d0,dt_norm) ...
                    - sporo_tmp* msupp_tmp * merge(dt,1d0,dt_norm) ...
                    ) ...
                    * 1d0;
                            
                if (iz~=nz)  
                    col = iz + 1;
                    amx3(row,col) = ...
                        (- sporop_tmp * vol * dps(iz) * wp_tmp * dzinv) * merge(dt,1d0,dt_norm) * psdx(ips,min(iz+1,nz));
                end 
                
                for iiz = 1: nz
                    col = iiz;
                    trans_tmp = sum(trans(iiz,iz,:))/nsp_sld;
                    if (trans_tmp == 0d0); continue; end
                        
                    amx3(row,col) = amx3(row,col) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm);
                    ymx3(row) = ymx3(row) - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips)  * merge(dt,1d0,dt_norm);
                end
                
            end
        
            ymx3=-1.0d0*ymx3;

            if (any(isnan(amx3),'all')||any(isnan(ymx3))||any(amx3>infinity,'all')||any(ymx3>infinity))  
            % if (true)  
                fprintf('PSD--%s: error in mtx\n',chrsp);
                fprintf('PSD--%s: any(isnan(amx3)),any(isnan(ymx3))\n',chrsp);
                fprintf(string(any(isnan(amx3),'all')),string(any(isnan(ymx3))) );

                if (any(isnan(ymx3)))  
                    for iz = 1: nz
                        if (isnan(ymx3(iz)))  
                            fprintf('%s\t%d\t%d\n', 'NAN is here...',ips,iz);
                        end
                    end 
                end


                if (any(isnan(amx3),'all'))  
                    for ie = 1:(nz)
                        for ie2 = 1:(nz)
                            if (isnan(amx3(ie,ie2)))  
                                fprintf('%s\t%d\t%d\t%d\n', 'PSD: NAN is here...',ips,ie,ie2);
                            end
                        end
                    end
                end
                error('stop');
            end
        
            % call DGESV(Nz,int(1),amx3,Nz,IPIV3,ymx3,Nz,INFO) 
            [xmx3,rmx3] = linsolve(amx3,ymx3);
            ymx3 = xmx3;
        
            if (any(isnan(ymx3))) 
                fprintf( 'PSD--%s: error in soultion\n',chrsp);
                flgback = true;
                return
            end
        
            for iz = 1: nz
                
                row =  iz;

                if (isnan(ymx3(row)))  
                    fprintf('PSD--%s: nan at\t',chrsp); 
                    fprintf('%d\t%d\n',iz,ips);
                    error('stop');
                end
                
                emx3_loc(row) = dps(ips)*psdx(ips,iz)*exp(ymx3(row)) - dps(ips)*psdx(ips,iz);
                
                if ((~isnan(ymx3(row)))&&ymx3(row) >threshold)  
                    psdx(ips,iz) = psdx(ips,iz)*corr;
                elseif (ymx3(row) < -threshold)  
                    psdx(ips,iz) = psdx(ips,iz)/corr;
                else   
                    psdx(ips,iz) = psdx(ips,iz)*exp(ymx3(row));
                end
            end 

            if (all(fact_tol == 1d0))  
                emx3(ips) = max(exp(abs(ymx3))) - 1.0d0;
            else 
                emx3(ips) = max(abs(emx3_loc));
            end 
        
        end 
        
        error = max(emx3(:)./fact_tol(:));

        if ( isnan(error) || any(isnan(psdx),'all') )  
            error = 1d3;
            fprintf('PSD--%s: !! error is NaN; values are returned to those before iteration with reducing dt\n',chrsp);
            fprintf('PSD--%s: isnan(error), info~=0,any(isnan(pdsx))\n',chrsp);
            fprintf('%s\t%s\n', string(isnan(error)), string(any(isnan(psdx),'all')) );
            
            flgback = true;
            error('stop')
        end

        fprintf( 'PSD--%s: iteration error = %7.6e, iteration = %d, time step [yr] = %7.6e\n',chrsp,error, iter,dt);
        iter = iter + 1; 
        
        if (iter > iter_max ) 
            if (dt==0d0)  
                fprintf(chrsp,'dt==0d0\n');
                error('stop');
            end 
            flgback = true;
            
            return
        end 
        
        if (flgback); return; end

    end 

    % calculating flux 
    flx_psd(:,:,:) = 0d0;
    for ips = 1: nps
        
        for iz = 1: nz
        
            row =  iz;
            
            vol  = 4d0/3d0*pi*(10d0^ps(ips))^3d0;
            surf = 4d0*pi*(10d0^ps(ips))^2d0;
                    
            m_tmp = vol * psdx(ips,iz) * dps(ips);
            mprev_tmp = vol * psd(ips,iz) * dps(ips);        
            % rxn_tmp = vol * psdx(ips,iz)*dps(ips) ...
                % * surf * psdx(ips,iz)*dps(ips) * kpsd(iz); 
            % rxn_tmp = surf * psdx(ips,iz)*dps(ips) * kpsd(iz); 
            rxn_tmp =  - vol * dpsd(ips,iz) * dps(ips) / dt ;
            
            msupp_tmp = vol * psd_rain(ips,iz) * dps(ips)  / dt;
            
            mi_tmp = vol * psd_pr(ips) * dps(ips);
            mp_tmp = vol * psdx(ips,min(iz+1,nz)) * dps(ips);
            
            dzinv = 1d0/dz(iz);
            
            % drxn_tmp = ... 
                % vol * 1d0 * dps(ips) ...
                % * surf * psdx(ips,iz) * dps(ips) * kpsd(iz) ...
                % + vol * psdx(ips,iz) * dps(ips) ...
                % * su;
            % drxn_tmp = surf * 1d0 * dps(ips) * kpsd(iz) ;
            drxn_tmp = 0d0;
            
            w_tmp = w(iz);
            wp_tmp = w(min(nz,iz+1));

            sporo_tmp = 1d0-poro(iz);
            sporop_tmp = 1d0-poro(min(nz,iz+1));
            sporoprev_tmp = 1d0-poroprev(iz);
            
            if (iz==nz)  
                mp_tmp = mi_tmp;
                wp_tmp = w0;
                sporop_tmp = 1d0- poroi;
            end 
            
            sporo_tmp = 1d0;
            sporop_tmp = 1d0;
            sporoprev_tmp = 1d0;
            
            flx_psd(ips,itflx_psd,iz) = ( ...
                ( sporo_tmp * m_tmp - sporoprev_tmp*mprev_tmp ) * dtinv  ...
                );
            flx_psd(ips,iadv_psd,iz) = ( ...
                -( sporop_tmp * wp_tmp * mp_tmp - sporo_tmp * w_tmp * m_tmp ) * dzinv ...
                );
            flx_psd(ips,irxn_psd,iz) = ( ...
                + sporo_tmp* rxn_tmp  ...
                );
            flx_psd(ips,irain_psd,iz) = ( ...
                - sporo_tmp* msupp_tmp  ...
                );
            
            for iiz = 1: nz  
                trans_tmp = sum(trans(iiz,iz,:))/nsp_sld;
                if (trans_tmp == 0d0); continue; end
                
                flx_psd(ips,idif_psd,iz) = flx_psd(ips,idif_psd,iz) + ( ...
                    - trans_tmp * sporo(iiz) * vol * psdx(ips,iiz) * dps(ips) ...
                    );
            end
            
            
            flx_psd(ips,ires_psd,iz) = sum(flx_psd(ips,:,iz));
            
        end
    end
       

            
    if ( chkflx && dt > dt_th)
        for ips = 1: nps
            flx_max = 0d0;
            for iflx=1:nflx_psd 
                flx_max= max( flx_max, abs( sum(squeeze(flx_psd(ips,iflx,:)).*dz(:)) ) );
            end 
            flx_max_max = max( flx_max_max, flx_max);
        end 
        for ips = 1: nps
        
            if ( flx_max > flx_max_max*flx_max_tol && abs( sum(squeeze(flx_psd(ips,ires_psd,:)).*dz(:)) ) > flx_max * flx_tol )  
                
                fprintf('%s too large error in PSD flx?\n',chrsp);
                fprintf('flx_max, flx_max_max,tol = %7.6e, %7.6e, %7.6e\n', flx_max,flx_max_max,flx_max_tol);
                fprintf('res, max = %7.6e, %7.6e\n', abs( sum(squeeze(flx_psd(ips,ires_psd,:)).*dz(:)) ), flx_max);
                fprintf('res/max, target = %7.6e, %7.6e\n', abs( sum(squeeze(flx_psd(ips,ires_psd,:)).*dz(:)) )/flx_max, flx_tol);
                
                flgback = true;
            
            end 
            
        end 
    end  
   
end
