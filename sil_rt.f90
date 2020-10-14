

subroutine silicate_dis_co2_1D( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,resp,poroprev,daqc,dgasc  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,pco2,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,khco2i,pco2i,pco2th,torg,ucv &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,pco2x,flx_co2,khco2x &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::dz,w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,daqc,dgasc,khco2i,pco2i,pco2th,ucv
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,pco2,man,kan,mansupp,mcc,kcc,mccsupp,resp,poroprev,torg
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx,pco2x,khco2x 
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,flx_co2
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 9
integer iz,row,nmx,ie,ie2,isp,iflx,col
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg,domega_cc_dco2,domega_ab_dco2,domega_an_dco2,domega_fo_dco2,dprodco2,preccc,dkhco2_dna &
    & ,dkhco2_dmg,dkhco2_dca,dkhco2_dco2,dkhco2_dsi,dpreccc_dco2,dpreccc_dna,dpreccc_dca,dpreccc_dmg,dpreccc_dsi,khco2 &
    & ,edif,alpha,alphaprev,dalpha_dco2,dalpha_dca,dalpha_dna,dalpha_dmg,dalpha_dsi,dedif_dco2,dedif_dca,dedif_dna,dedif_dmg &
    & ,dedif_dsi,dpreccc_dcc 
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc,edifi,dalpha_disp &
    & ,dkhco2_disp,dkhco2_disp_n,dkhco2_disp_p,dedif_disp,dedif_disp_n,dedif_disp_p,dpreccc_disp
real(kind=8)::k1_fo = 10d0**(-6.85d0), E1_fo = 51.7d0, n1_fo = 0.5d0, k2_fo = 10d0**(-12.41d0),E2_fo = 38d0 &
    & ,k3_fo = 10d0**(-21.2d0),E3_fo = 94.1d0,n3_fo = -0.82d0  &
    & ,k1_ab = 10d0**(-10.16d0), E1_ab = 65d0, n1_ab = 0.457d0, k2_ab = 10d0**(-12.56d0),E2_ab = 69.8d8 &
    & ,k3_ab = 10d0**(-15.60d0),E3_ab = 71d0, n3_ab = -0.572d0 &
    & ,k1_an = 10d0**(-3.5d0),E1_an = 16.6d0,n1_an = 1.411d0,k2_an = 10d0**(-9.12d0), E2_an = 17.8d0 

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8),parameter::infinity = huge(0d0)
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
    flx_co2 = 0d0

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
    call calc_pH( &
        & nz,nax+2d0*(cx+cax+mgx-so4x)+3d0*c2x,pco2x+dconc,kw,kco2,k1,k2 &! input 
        & ,dprodco2 &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodco2 = (dprodco2-prox)/dconc
    dprodca = dprodmg
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    omega_fo = mgx**2d0*six/(prox**4d0)/keqfo
    domega_fo_dmg = 2d0*mgx*six/(prox**4d0)/keqfo + mgx**2d0*six*(-4d0)/(prox**5d0)*dprodmg/keqfo
    domega_fo_dsi = mgx**2d0/(prox**4d0)/keqfo
    domega_fo_dna = mgx**2d0*six*(-4d0)/(prox**5d0)*dprodna/keqfo
    domega_fo_dca = mgx**2d0*six*(-4d0)/(prox**5d0)*dprodca/keqfo
    domega_fo_dco2 = mgx**2d0*six*(-4d0)/(prox**5d0)*dprodco2/keqfo
    
    ! omega_fo = mg**2d0*si/(pro**4d0)/keqfo
    ! domega_fo_dmg = 0d0
    ! domega_fo_dsi = 0d0
    
    ! print *,omega_fo
    ! print *,domega_fo_dmg
    ! print *,domega_fo_dsi
    ! stop
    
    ! Ab + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)
    omega_ab = nax*six**2d0/prox/keqab
    domega_ab_dna = six**2d0/prox/keqab + nax*six**2d0*(-1d0)/(prox**2d0)/keqab*dprodna
    domega_ab_dsi = nax*(2d0)*six/prox/keqab
    domega_ab_dmg = nax*six**2d0*(-1d0)/(prox**2d0)/keqab*dprodmg
    domega_ab_dca = nax*six**2d0*(-1d0)/(prox**2d0)/keqab*dprodca
    domega_ab_dco2 = nax*six**2d0*(-1d0)/(prox**2d0)/keqab*dprodco2
    
    ! An + 2H+ + H2O = kaolinite + Ca2+ 
    omega_an = cax/(prox**2d0)/keqan
    domega_an_dca = 1d0/(prox**2d0)/keqan + cax*(-2d0)/(prox**3d0)*dprodca/keqan
    domega_an_dna = cax*(-2d0)/(prox**3d0)*dprodna/keqan
    domega_an_dmg = cax*(-2d0)/(prox**3d0)*dprodmg/keqan
    domega_an_dco2 = cax*(-2d0)/(prox**3d0)*dprodco2/keqan
    domega_an_dsi = 0d0
    
    ! Cc = Ca2+ + CO32- 
    omega_cc = cax*k1*k2*kco2*pco2x/(prox**2d0)/keqcc
    domega_cc_dca = k1*k2*kco2*pco2x/(prox**2d0)/keqcc + cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodca/keqcc
    domega_cc_dmg = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodmg/keqcc
    domega_cc_dna = cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodna/keqcc
    domega_cc_dco2 = cax*k1*k2*kco2/(prox**2d0)/keqcc + cax*k1*k2*kco2*pco2x*(-2d0)/(prox**3d0)*dprodco2/keqcc
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

                amx3(row,row) = (1.0d0  &
                    & + w/dz*dt  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp)*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp) &
                    & -w*(mi_tmp-m_tmp)/dz*dt &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp)*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & -msupp_tmp*dt  &
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

                amx3(row,row) = (1.0d0     &
                    & + w/dz*dt    &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*(1d0-omega_tmp)*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & * merge(1.0d0,m_tmp,m_tmp<mth_tmp)

                amx3(row,row+nsp3) = (-w/dz)*dt *merge(1.0d0,mp_tmp,m_tmp<mth_tmp)

                ymx3(row) = ( &
                    & (m_tmp-mprev_tmp) &
                    & -w*(mp_tmp-m_tmp)/dz*dt  &
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp)*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & -msupp_tmp*dt  &
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

                amx3(row,row + 8 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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

            if (.not.((iz == 1).or.(iz==nz))) then

                amx3(row,row) = ( &
                    & poro(iz)*sat(iz)*1d3  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(-2d0)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)/dz*dt  &
                    & -drxndisp_tmp*dt &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0))*dt &
                    & - poro(iz)*sat(iz)*1d3*v(iz)/dz*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                amx3(row,row+nsp3) = ( &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(1d0)/(dz**2d0))*dt &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz*dt &
                    & - rxn_tmp*dt &
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
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caq_tmp_n-2d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == 1) then

                amx3(row,row) = ( &
                    & poro(iz)*sat(iz)*1d3  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(-2d0)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)/dz*dt &
                    & - drxndisp_tmp*dt &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row+nsp3) = ( &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(1d0)/(dz**2d0))*dt &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caqi_tmp)/dz*dt &
                    & - rxn_tmp*dt &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caqi_tmp)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_p+caqi_tmp-2d0*caq_tmp)/(dz**2d0)) &
                        & ) 
                endif 

            else if (iz == nz) then

                amx3(row,row) = ( &
                    & poro(iz)*sat(iz)*1d3  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(-1d0)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(1d0)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)/dz*dt &
                    & - drxndisp_tmp*dt &
                    & ) &
                    & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

                amx3(row,row-nsp3) = ( &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(1d0)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1))*(-1d0)/(dz**2d0))*dt &
                    & - v(iz)/dz*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)

                ymx3(row) = ( &
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                    & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                    & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                    & *(caq_tmp-caq_tmp_n)/(dz**2d0))*dt &
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz*dt &
                    & - rxn_tmp*dt &
                    & ) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                
                if (isp==1) then 
                    flx_mg(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_mg(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==2) then 
                    flx_si(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_si(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==3) then 
                    flx_na(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_na(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                elseif (isp==4) then 
                    flx_ca(iadv,iz) = (&
                        & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz &
                        & ) 
                    flx_ca(idif,iz) = (&
                        & +(-poro(iz)*sat(iz)*1d3*d_tmp*tora(iz)*(caq_tmp_n-1d0*caq_tmp)/(dz**2d0) &
                        & -d_tmp*1d3*(poro(iz)*sat(iz)*tora(iz)-poro(iz-1)*sat(iz-1)*tora(iz-1)) &
                        & *(caq_tmp-caq_tmp_n)/(dz**2d0)) &
                        & ) 
                endif 


            end if 
            
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
            
                amx3(row,row  + 4) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
            
                amx3(row,row  + 3) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
            
                amx3(row,row  + 2) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
            
                amx3(row,row  + 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
    khco2x = kco2*(1d0+k1/prox + k1*k2/prox/prox)
    preccc = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(1d0-omega_cc) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dkhco2_dco2 = kco2*(k1*(-1d0)/(prox**2d0) + k1*k2*(-2d0)/(prox**3d0))*dprodco2
    dkhco2_dna = kco2*(k1*(-1d0)/(prox**2d0) + k1*k2*(-2d0)/(prox**3d0))*dprodna
    dkhco2_dca = kco2*(k1*(-1d0)/(prox**2d0) + k1*k2*(-2d0)/(prox**3d0))*dprodca
    dkhco2_dmg = kco2*(k1*(-1d0)/(prox**2d0) + k1*k2*(-2d0)/(prox**3d0))*dprodmg
    dkhco2_dsi = 0d0
    
    dpreccc_dco2 = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(-domega_cc_dco2) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dpreccc_dmg = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(-domega_cc_dmg) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dpreccc_dna = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(-domega_cc_dna) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dpreccc_dca = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(-domega_cc_dca) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dpreccc_dsi = kcc*poro(iz)*hr(iz)*mvcc*1d-6*mccx*(-domega_cc_dsi) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    dpreccc_dcc = kcc*poro(iz)*hr(iz)*mvcc*1d-6*(1d0-omega_cc) &
                    & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0) 
    
    ! preccc=0d0
    ! dpreccc_dca=0d0
    ! dpreccc_dcc=0d0
    ! dpreccc_dco2=0d0
    ! dpreccc_dmg=0d0
    ! dpreccc_dsi=0d0
    ! dpreccc_dna=0d0
    
    ! dkhco2_dco2 = 0d0
    ! dkhco2_dna = 0d0
    ! dkhco2_dca = 0d0
    ! dkhco2_dmg = 0d0
    ! dkhco2_dsi = 0d0
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc+poro*sat*khco2x*1d3*tora*daqc
    edifi = edif(1)
    edifi = ucv*1d3*dgasc 
    
    dedif_dco2 = poro*sat*dkhco2_dco2*1d3*tora*daqc
    dedif_dca = poro*sat*dkhco2_dca*1d3*tora*daqc
    dedif_dmg = poro*sat*dkhco2_dmg*1d3*tora*daqc
    dedif_dna = poro*sat*dkhco2_dna*1d3*tora*daqc
    dedif_dsi = poro*sat*dkhco2_dsi*1d3*tora*daqc
    
    alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2x*1d3
    alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khco2*1d3
    
    dalpha_dco2 = poro*sat*dkhco2_dco2*1d3
    dalpha_dca = poro*sat*dkhco2_dca*1d3
    dalpha_dna = poro*sat*dkhco2_dna*1d3
    dalpha_dmg = poro*sat*dkhco2_dmg*1d3
    dalpha_dsi = poro*sat*dkhco2_dsi*1d3
    
    
    
    do iz = 1, nz

        row = nsp3*(iz-1) + 9

        if (iz == 1) then

            amx3(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/dz &
                & +0.5d0*(dedif_dco2(iz))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edifi)*(1d0)/dz &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2i)/dz )/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx3(row,row+nsp3) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/dz &
                & + 0.5d0*(dedif_dco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz)/dz*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            ymx3(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edifi)*(pco2x(iz)-pco2i)/dz )/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2i*pco2i)/dz *dt&
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            do isp=1,4
                if (isp==1) then    
                    dedif_disp = dalpha_dmg(iz)
                    dedif_disp_p = dalpha_dmg(iz+1)
                    dkhco2_disp = dkhco2_dmg(iz)
                    dkhco2_disp_p = dkhco2_dmg(iz+1)
                    dalpha_disp = dalpha_dmg(iz)
                    caq_tmp = mgx(iz)
                    caq_tmp_p = mgx(iz+1)
                    dpreccc_disp = dpreccc_dmg(iz)
                elseif (isp==2) then 
                    dedif_disp = dalpha_dsi(iz)
                    dedif_disp_p = dalpha_dsi(iz+1)
                    dkhco2_disp = dkhco2_dsi(iz)
                    dkhco2_disp_p = dkhco2_dsi(iz+1)
                    dalpha_disp = dalpha_dsi(iz)
                    caq_tmp = six(iz)
                    caq_tmp_p = six(iz+1)
                    dpreccc_disp = dpreccc_dsi(iz)
                elseif (isp==3) then 
                    dedif_disp = dalpha_dna(iz)
                    dedif_disp_p = dalpha_dna(iz+1)
                    dkhco2_disp = dkhco2_dna(iz)
                    dkhco2_disp_p = dkhco2_dna(iz+1)
                    dalpha_disp = dalpha_dna(iz)
                    caq_tmp = nax(iz)
                    caq_tmp_p = nax(iz+1)
                    dpreccc_disp = dpreccc_dna(iz)
                elseif (isp==4) then 
                    dedif_disp = dalpha_dca(iz)
                    dedif_disp_p = dalpha_dca(iz+1)
                    dkhco2_disp = dkhco2_dca(iz)
                    dkhco2_disp_p = dkhco2_dca(iz+1)
                    dalpha_disp = dalpha_dca(iz)
                    caq_tmp = cax(iz)
                    caq_tmp_p = cax(iz+1)
                    dpreccc_disp = dpreccc_dca(iz)
                endif 
                
                col = (iz-1)*nsp3 + 4 + isp 
                
                amx3(row,col) = ( &
                    & (dalpha_disp*pco2x(iz)) &
                    & -( 0.5d0*(dedif_disp)*(pco2x(iz+1)-pco2x(iz))/dz &
                    & - 0.5d0*(dedif_disp)*(pco2x(iz)-pco2i)/dz )/dz*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_disp*pco2x(iz))/dz*dt &
                    & -dpreccc_disp*dt &
                    & ) &
                    & *caq_tmp &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
                    
                amx3(row,col+nsp3) = ( &
                    & -( 0.5d0*(dedif_disp_p)*(pco2x(iz+1)-pco2x(iz))/dz )/dz*dt  &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
            enddo
            
            flx_co2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edifi)*(pco2x(iz)-pco2i)/dz )/dz  &
                & )
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2i*pco2i)/dz &
                & )

        else if (iz == nz) then

            amx3(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)) &
                & -( - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/dz &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2x(iz-1))/dz)/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx3(row,row-nsp3) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/dz &
                & - 0.5d0*(dedif_dco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz)/dz *dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dco2(iz-1)*pco2x(iz-1))/dz*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx3(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            do isp=1,4
                if (isp==1) then    
                    dedif_disp = dalpha_dmg(iz)
                    dedif_disp_n = dalpha_dmg(iz-1)
                    dkhco2_disp = dkhco2_dmg(iz)
                    dkhco2_disp_n = dkhco2_dmg(iz-1)
                    dalpha_disp = dalpha_dmg(iz)
                    caq_tmp = mgx(iz)
                    caq_tmp_n = mgx(iz-1)
                    dpreccc_disp = dpreccc_dmg(iz)
                elseif (isp==2) then 
                    dedif_disp = dalpha_dsi(iz)
                    dedif_disp_n = dalpha_dsi(iz-1)
                    dkhco2_disp = dkhco2_dsi(iz)
                    dkhco2_disp_n = dkhco2_dsi(iz-1)
                    dalpha_disp = dalpha_dsi(iz)
                    caq_tmp = six(iz)
                    caq_tmp_n = six(iz-1)
                    dpreccc_disp = dpreccc_dsi(iz)
                elseif (isp==3) then 
                    dedif_disp = dalpha_dna(iz)
                    dedif_disp_n = dalpha_dna(iz-1)
                    dkhco2_disp = dkhco2_dna(iz)
                    dkhco2_disp_n = dkhco2_dna(iz-1)
                    dalpha_disp = dalpha_dna(iz)
                    caq_tmp = nax(iz)
                    caq_tmp_n = nax(iz-1)
                    dpreccc_disp = dpreccc_dna(iz)
                elseif (isp==4) then 
                    dedif_disp = dalpha_dca(iz)
                    dedif_disp_n = dalpha_dca(iz-1)
                    dkhco2_disp = dkhco2_dca(iz)
                    dkhco2_disp_n = dkhco2_dca(iz-1)
                    dalpha_disp = dalpha_dca(iz)
                    caq_tmp = cax(iz)
                    caq_tmp_n = cax(iz-1)
                    dpreccc_disp = dpreccc_dca(iz)
                endif 
                col = (iz-1)*nsp3 + 4 + isp
                amx3(row,col) = ( &
                    & (dalpha_disp*pco2x(iz)) &
                    & -(- 0.5d0*(dedif_disp)*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_disp*pco2x(iz))/dz*dt &
                    & -dpreccc_disp*dt &
                    & ) &
                    & *caq_tmp &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
                    
                amx3(row,col-nsp3) = ( &
                    & -( - 0.5d0*(dedif_disp_n)*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_disp_n*pco2x(iz-1))/dz*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
            enddo
            
            flx_co2(idif,iz) = ( &
                & -( 0d0 &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz  &
                & ) 
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz &
                & )

        else

            amx3(row,row) = ( &
                & (alpha(iz) + dalpha_dco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(-1d0)/dz &
                & +0.5d0*(dedif_dco2(iz))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(1d0)/dz  &
                & - 0.5d0*(dedif_dco2(iz))*(pco2x(iz)-pco2x(iz-1))/dz )/dz *dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dco2(iz)*pco2x(iz))/dz*dt &
                & -dpreccc_dco2(iz)*dt &
                & ) &
                & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

            amx3(row,row+nsp3) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/dz &
                & + 0.5d0*(dedif_dco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz )/dz*dt  &
                & ) &
                & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/dz &
                & - 0.5d0*(dedif_dco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dco2(iz-1)*pco2x(iz-1))/dz*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            ymx3(row) = ( &
                & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz*dt &
                & -resp(iz)*dt &
                & -preccc(iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
            
            do isp=1,4
                if (isp==1) then    
                    dedif_disp = dalpha_dmg(iz)
                    dedif_disp_p = dalpha_dmg(iz+1)
                    dedif_disp_n = dalpha_dmg(iz-1)
                    dkhco2_disp = dkhco2_dmg(iz)
                    dkhco2_disp_p = dkhco2_dmg(iz+1)
                    dkhco2_disp_n = dkhco2_dmg(iz-1)
                    dalpha_disp = dalpha_dmg(iz)
                    caq_tmp = mgx(iz)
                    caq_tmp_p = mgx(iz+1)
                    caq_tmp_n = mgx(iz-1)
                    dpreccc_disp = dpreccc_dmg(iz)
                elseif (isp==2) then 
                    dedif_disp = dalpha_dsi(iz)
                    dedif_disp_p = dalpha_dsi(iz+1)
                    dedif_disp_n = dalpha_dsi(iz-1)
                    dkhco2_disp = dkhco2_dsi(iz)
                    dkhco2_disp_p = dkhco2_dsi(iz+1)
                    dkhco2_disp_n = dkhco2_dsi(iz-1)
                    dalpha_disp = dalpha_dsi(iz)
                    caq_tmp = six(iz)
                    caq_tmp_p = six(iz+1)
                    caq_tmp_n = six(iz-1)
                    dpreccc_disp = dpreccc_dsi(iz)
                elseif (isp==3) then 
                    dedif_disp = dalpha_dna(iz)
                    dedif_disp_p = dalpha_dna(iz+1)
                    dedif_disp_n = dalpha_dna(iz-1)
                    dkhco2_disp = dkhco2_dna(iz)
                    dkhco2_disp_p = dkhco2_dna(iz+1)
                    dkhco2_disp_n = dkhco2_dna(iz-1)
                    dalpha_disp = dalpha_dna(iz)
                    caq_tmp = nax(iz)
                    caq_tmp_p = nax(iz+1)
                    caq_tmp_n = nax(iz-1)
                    dpreccc_disp = dpreccc_dna(iz)
                elseif (isp==4) then 
                    dedif_disp = dalpha_dca(iz)
                    dedif_disp_p = dalpha_dca(iz+1)
                    dedif_disp_n = dalpha_dca(iz-1)
                    dkhco2_disp = dkhco2_dca(iz)
                    dkhco2_disp_p = dkhco2_dca(iz+1)
                    dkhco2_disp_n = dkhco2_dca(iz-1)
                    dalpha_disp = dalpha_dca(iz)
                    caq_tmp = cax(iz)
                    caq_tmp_p = cax(iz+1)
                    caq_tmp_n = cax(iz-1)
                    dpreccc_disp = dpreccc_dca(iz)
                endif 
                col = (iz-1)*nsp3 + 4 + isp
                amx3(row,col) = ( &
                    & (dalpha_disp*pco2x(iz)) &
                    & -( 0.5d0*(dedif_disp)*(pco2x(iz+1)-pco2x(iz))/dz &
                    & - 0.5d0*(dedif_disp)*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt & 
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_disp*pco2x(iz))/dz*dt &
                    & -dpreccc_disp*dt &
                    & ) &
                    & *caq_tmp &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
                    
                amx3(row,col-nsp3) = ( &
                    & -( - 0.5d0*(dedif_disp_n)*(pco2x(iz)-pco2x(iz-1))/dz )/dz*dt  &
                    & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_disp_n*pco2x(iz-1))/dz*dt &
                    & ) &
                    & *caq_tmp_n &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
                    
                amx3(row,col+nsp3) = ( &
                    & -( 0.5d0*(dedif_disp_p)*(pco2x(iz+1)-pco2x(iz))/dz )/dz*dt  &
                    & ) &
                    & *caq_tmp_p &
                    & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
            enddo
            
            flx_co2(idif,iz) = ( &
                & -( 0.5d0*(edif(iz)+edif(iz+1))*(pco2x(iz+1)-pco2x(iz))/dz &
                & - 0.5d0*(edif(iz)+edif(iz-1))*(pco2x(iz)-pco2x(iz-1))/dz )/dz  &
                & ) 
            flx_co2(iadv,iz) = ( &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2x(iz-1)*pco2x(iz-1))/dz &
                & )

        end if 
        
        col = (iz-1)*nsp3 + 4 ! calcite
        amx3(row,col) = ( &
            & -dpreccc_dcc(iz)*dt &
            & ) &
            & *mccx(iz) &
            & *merge(0.0d0,1d0,pco2x(iz)<pco2th)
        
        flx_co2(itflx,iz) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz))/dt &
            & ) 
        flx_co2(irxn_fo,iz) = -resp(iz)
        flx_co2(irain,iz) = -preccc(iz)
        flx_co2(ires,iz) = sum(flx_co2(:,iz))
        
        if (any(isnan(flx_co2(:,iz)))) then
            print *,flx_co2(:,iz)
        endif 

    end do 
    
    ymx3=-1.0d0*ymx3

    if (any(isnan(amx3)).or.any(isnan(ymx3)).or. any(amx3>infinity).or. any(ymx3>infinity) ) then 
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
        print*,'sil_dis: error in soultion'
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
        
        row = 9 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'pco2'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                pco2x(iz) = 0.1d0*pco2th
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            pco2x(iz) = pco2x(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            pco2x(iz) = pco2x(iz)*0.50d0
        else
            pco2x(iz) = pco2x(iz)*exp(ymx3(row))
        endif
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)).or. any(isnan(cax)).or. any(isnan(pco2x)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx)).or. any(isnan(manx)).or. any(isnan(mccx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx)) &
            & ,any(isnan(cax)),any(isnan(manx)),any(isnan(mccx)),any(isnan(pco2x))
        
        open(unit=11,file='amx.txt',status = 'replace')
        open(unit=12,file='ymx.txt',status = 'replace')
        do ie = 1,nsp3*(nz)
            write(11,*) (amx3(ie,ie2),ie2 = 1,nsp3*nz)
            write(12,*) ymx3(ie)
        enddo 
        close(11)
        close(12)       
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
        pco2x = pco2
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
        
        row = 9 + nsp3*(iz-1)

        if (pco2x(iz) < 0.0d0) then
            pco2x(iz) = pco2x(iz)/exp(ymx3(row))*0.5d0
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

endsubroutine silicate_dis_co2_1D



!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine silicate_dis_co2_1D_v2( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
    & ,pco2th,resp &! intput
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2 &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,k1si,k2si,keqcca,authig,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3,pco2i,khco2i,ucv,dgasc,daqc,pco2th
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz,kcca,pco2,torg,resp
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx,pco2x
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc,omega_cca
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,flx_co2
integer,intent(inout)::iter,it
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 9
integer iz,row,nmx,ie,ie2,isp,iflx
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg,dprodsi,domega_fo_dpro,domega_ab_dpro,domega_an_dpro,domega_cc_dpro &
    & ,domega_cca_dmg,domega_cca_dca,domega_cca_dna,domega_cca_dsi,domega_cca_dpro,dprodpco2,domega_fo_dpco2,domega_ab_dpco2 &
    & ,domega_an_dpco2,domega_cc_dpco2,domega_cca_dpco2,dkhco2_dpro,dkhco2_dpco2,dpreccc_dpco2,preccc,khco2,khco2x,alpha,alphaprev &
    & ,dalpha_dpco2,edif,dedif_dpco2,dalpha_dna,dalpha_dsi,dalpha_dmg,dalpha_dca,dedif_dca,dedif_dmg,dedif_dsi,dedif_dna &
    & ,dkhco2_dca,dkhco2_dna,dkhco2_dsi,dkhco2_dmg,dpreccc_dca,dpreccc_dna,dpreccc_dsi,dpreccc_dmg,dpreccc_dmcc
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,st_cca,edifi,khco2n_tmp,pco2n_tmp,edifn_tmp
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
real(kind=8),dimension(nz)::disonly ! for cc [1---yes, 0---no]
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

    call calc_pH_v3( &
        & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,prox &! output
        & ) 
else 

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
    flx_co2 = 0d0

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
    call calc_pH_v3( &
        & nz,nax,mgx,cax,so4x,pco2x+dconc,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,dprodpco2 &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodsi = (dprodsi-prox)/dconc
    dprodca = (dprodca-prox)/dconc
    dprodpco2 = (dprodpco2-prox)/dconc
    
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
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x+dconc,cax,mgx,six,nax,prox,'fo' &! input 
        & ,domega_fo_dpco2 &! output
        & )
    domega_fo_dmg = (domega_fo_dmg-omega_fo)/dconc
    domega_fo_dsi = (domega_fo_dsi-omega_fo)/dconc
    domega_fo_dna = (domega_fo_dna-omega_fo)/dconc
    domega_fo_dca = (domega_fo_dca-omega_fo)/dconc
    domega_fo_dpro = (domega_fo_dpro-omega_fo)/dconc
    domega_fo_dpco2 = (domega_fo_dpco2-omega_fo)/dconc
    
    domega_fo_dmg = domega_fo_dmg + domega_fo_dpro*dprodmg
    domega_fo_dca = domega_fo_dca + domega_fo_dpro*dprodca
    domega_fo_dna = domega_fo_dna + domega_fo_dpro*dprodna
    domega_fo_dsi = domega_fo_dsi + domega_fo_dpro*dprodsi
    domega_fo_dpco2 = domega_fo_dpco2 + domega_fo_dpro*dprodpco2
    
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
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x+dconc,cax,mgx,six,nax,prox,'ab' &! input 
        & ,domega_ab_dpco2 &! output
        & )
    domega_ab_dmg = (domega_ab_dmg-omega_ab)/dconc
    domega_ab_dsi = (domega_ab_dsi-omega_ab)/dconc
    domega_ab_dna = (domega_ab_dna-omega_ab)/dconc
    domega_ab_dca = (domega_ab_dca-omega_ab)/dconc
    domega_ab_dpro = (domega_ab_dpro-omega_ab)/dconc
    domega_ab_dpco2 = (domega_ab_dpco2-omega_ab)/dconc
    
    domega_ab_dmg = domega_ab_dmg + domega_ab_dpro*dprodmg
    domega_ab_dca = domega_ab_dca + domega_ab_dpro*dprodca
    domega_ab_dna = domega_ab_dna + domega_ab_dpro*dprodna
    domega_ab_dsi = domega_ab_dsi + domega_ab_dpro*dprodsi
    domega_ab_dpco2 = domega_ab_dpco2 + domega_ab_dpro*dprodpco2
    
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
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,'an' &! input 
        & ,domega_an_dpco2 &! output
        & )
    domega_an_dmg = (domega_an_dmg-omega_an)/dconc
    domega_an_dsi = (domega_an_dsi-omega_an)/dconc
    domega_an_dna = (domega_an_dna-omega_an)/dconc
    domega_an_dca = (domega_an_dca-omega_an)/dconc
    domega_an_dpro = (domega_an_dpro-omega_an)/dconc
    domega_an_dpco2 = (domega_an_dpco2-omega_an)/dconc
    
    domega_an_dmg = domega_an_dmg + domega_an_dpro*dprodmg
    domega_an_dca = domega_an_dca + domega_an_dpro*dprodca
    domega_an_dna = domega_an_dna + domega_an_dpro*dprodna
    domega_an_dsi = domega_an_dsi + domega_an_dpro*dprodsi
    domega_an_dpco2 = domega_an_dpco2 + domega_an_dpro*dprodpco2
    
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
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,'cc' &! input 
        & ,domega_cc_dpco2 &! output
        & )
    domega_cc_dmg = (domega_cc_dmg-omega_cc)/dconc
    domega_cc_dsi = (domega_cc_dsi-omega_cc)/dconc
    domega_cc_dna = (domega_cc_dna-omega_cc)/dconc
    domega_cc_dca = (domega_cc_dca-omega_cc)/dconc
    domega_cc_dpro = (domega_cc_dpro-omega_cc)/dconc
    domega_cc_dpco2 = (domega_cc_dpco2-omega_cc)/dconc
    
    domega_cc_dmg = domega_cc_dmg + domega_cc_dpro*dprodmg
    domega_cc_dca = domega_cc_dca + domega_cc_dpro*dprodca
    domega_cc_dna = domega_cc_dna + domega_cc_dpro*dprodna
    domega_cc_dsi = domega_cc_dsi + domega_cc_dpro*dprodsi
    domega_cc_dpco2 = domega_cc_dpco2 + domega_cc_dpro*dprodpco2
    
    disonly = 0d0
    
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
    call calc_omega( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,'cc' &! input 
        & ,domega_cca_dpco2 &! output
        & )
    domega_cca_dmg = (domega_cca_dmg-omega_cca)/dconc
    domega_cca_dsi = (domega_cca_dsi-omega_cca)/dconc
    domega_cca_dna = (domega_cca_dna-omega_cca)/dconc
    domega_cca_dca = (domega_cca_dca-omega_cca)/dconc
    domega_cca_dpro = (domega_cca_dpro-omega_cca)/dconc
    domega_cca_dpco2 = (domega_cca_dpco2-omega_cca)/dconc
    
    domega_cca_dmg = domega_cca_dmg + domega_cca_dpro*dprodmg
    domega_cca_dca = domega_cca_dca + domega_cca_dpro*dprodca
    domega_cca_dna = domega_cca_dna + domega_cca_dpro*dprodna
    domega_cca_dsi = domega_cca_dsi + domega_cca_dpro*dprodsi
    domega_cca_dpco2 = domega_cca_dpco2 + domega_cca_dpro*dprodpco2

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
                disonly(iz) = merge(1d0,0d0,pco2x(iz)<pco2th)
                omega_tmp_th = omega_tmp*disonly(iz)
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

                amx3(row,row + 8 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dpco2(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0) &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(1d0-omega_cca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) 
                drxndisp_tmp = ( &
                    & st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0) &
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
            
                amx3(row,row  + 4) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
            
                amx3(row,row  + 3) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
            
                amx3(row,row  + 2) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dmg(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dsi(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dna(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dpco2(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
    khco2x = kco2*(1d0+k1/prox + k1*k2/prox/prox)
    
    disonly = merge(1d0,0d0,pco2x<pco2th)
    preccc = kcc*poro*hr*mvcc*1d-6*mccx*(1d0-omega_cc) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
        
    dpreccc_dpco2=0d0
    dkhco2_dpco2 = 0d0
        
    dpreccc_dpco2=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dpco2) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dna=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dna) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dmg=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dmg) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dsi=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dsi) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dca=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dca) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dmcc=kcc*poro*hr*mvcc*1d-6*(1d0-omega_cc) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dkhco2_dpro = kco2*(k1*(-1d0)/prox**2d0 + k1*k2*(-2d0)/prox**3d0)
    dkhco2_dpco2 = dkhco2_dpro*dprodpco2
    dkhco2_dna = dkhco2_dpro*dprodna
    dkhco2_dmg = dkhco2_dpro*dprodmg
    dkhco2_dsi = dkhco2_dpro*dprodsi
    dkhco2_dca = dkhco2_dpro*dprodca
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc+poro*sat*khco2x*1d3*tora*daqc
    edifi = edif(1)
    edifi = ucv*1d3*dgasc 
    
    dedif_dpco2 = poro*sat*dkhco2_dpco2*1d3*tora*daqc
    dedif_dna = poro*sat*dkhco2_dna*1d3*tora*daqc
    dedif_dmg = poro*sat*dkhco2_dmg*1d3*tora*daqc
    dedif_dsi = poro*sat*dkhco2_dsi*1d3*tora*daqc
    dedif_dca = poro*sat*dkhco2_dca*1d3*tora*daqc
    
    alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2x*1d3
    alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khco2*1d3
    
    dalpha_dpco2 = poro*sat*dkhco2_dpco2*1d3
    dalpha_dna = poro*sat*dkhco2_dna*1d3
    dalpha_dmg = poro*sat*dkhco2_dmg*1d3
    dalpha_dsi = poro*sat*dkhco2_dsi*1d3
    dalpha_dca = poro*sat*dkhco2_dca*1d3
    
    if (any(isnan(alpha)).or.any(isnan(alphaprev))) then 
        print *, alpha 
        print *, alphaprev
        print *, poroprev
        print *, khco2
        print *, pro
        stop
    endif 
    
    do iz = 1, nz

        row = nsp3*(iz-1) + 9
        
        pco2n_tmp = pco2x(max(1,iz-1))
        khco2n_tmp = khco2x(max(1,iz-1))
        edifn_tmp = edif(max(1,iz-1))
        if (iz == 1) then 
            pco2n_tmp = pco2i
            khco2n_tmp = khco2i
            edifn_tmp = edifi
        endif 

        amx3(row,row) = ( &
            & (alpha(iz) + dalpha_dpco2(iz)*pco2x(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & +0.5d0*(dedif_dpco2(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
            & - 0.5d0*(dedif_dpco2(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dpco2(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dpco2(iz)*dt &
            & ) &
            & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

        ymx3(row) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & ) &
            & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
        
        
        if (iz/=nz) then 
        
            amx3(row,row+nsp3) = ( &
                    & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & + 0.5d0*(dedif_dpco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-4) = ( &
                    & -( 0.5d0*(dedif_dmg(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,mgx(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-3) = ( &
                    & -( 0.5d0*(dedif_dsi(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,six(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-2) = ( &
                    & -( 0.5d0*(dedif_dna(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,nax(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-1) = ( &
                    & -( 0.5d0*(dedif_dca(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,cax(iz+1),pco2x(iz)<pco2th)
        
        endif 
        
        if (iz/=1) then 

            amx3(row,row-nsp3) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dpco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dpco2(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-4) = ( &
                & -(- 0.5d0*(dedif_dmg(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dmg(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,mgx(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-3) = ( &
                & -(- 0.5d0*(dedif_dsi(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dsi(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,six(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-2) = ( &
                & -(- 0.5d0*(dedif_dna(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dna(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,nax(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-1) = ( &
                & -(- 0.5d0*(dedif_dca(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dca(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,cax(iz-1),pco2x(iz)<pco2th)
        endif 
        
        amx3(row,row-5) = ( &
            & -dpreccc_dmcc(iz)*dt &
            & ) &
            & *merge(1.0d0,mccx(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-4) = ( &
            & (dalpha_dmg(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dmg(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dmg(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dmg(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dmg(iz)*dt &
            & ) &
            & *merge(1.0d0,mgx(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-3) = ( &
            & (dalpha_dsi(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dsi(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dsi(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dsi(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dsi(iz)*dt &
            & ) &
            & *merge(1.0d0,six(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-2) = ( &
            & (dalpha_dna(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dna(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dna(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dna(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dna(iz)*dt &
            & ) &
            & *merge(1.0d0,nax(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-1) = ( &
            & (dalpha_dca(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dca(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dca(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dca(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dca(iz)*dt &
            & ) &
            & *merge(1.0d0,cax(iz),pco2x(iz)<pco2th)
        
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
        flx_co2(irxn_fo,iz) = -resp(iz)
        flx_co2(irain,iz) = -preccc(iz)
        flx_co2(ires,iz) = sum(flx_co2(:,iz))
        
        if (any(isnan(flx_co2(:,iz)))) then
            print *,flx_co2(:,iz)
        endif 

    end do 
    
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
        
        row = 9 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'pco2'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                pco2x(iz) = 0.1d0*pco2th
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            pco2x(iz) = pco2x(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            pco2x(iz) = pco2x(iz)*0.50d0
        else
            pco2x(iz) = pco2x(iz)*exp(ymx3(row))
        endif
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)).or. any(isnan(cax)).or. any(isnan(pco2x)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx)).or. any(isnan(manx)).or. any(isnan(mccx))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx)) &
            & ,any(isnan(cax)),any(isnan(manx)),any(isnan(mccx)),any(isnan(pco2x))
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
        pco2x = pco2
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
        
        row = 9 + nsp3*(iz-1)

        if (pco2x(iz) < 0.0d0) then
            pco2x(iz) = pco2x(iz)/exp(ymx3(row))*0.5d0
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
    print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An, pCO2 -=-=-=-=-=-=-='
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
enddo

endsubroutine silicate_dis_co2_1D_v2













!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine alsilicate_dis_co2_1D_v2( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
    & ,pco2th,resp,k1al,k2al,k3al,k4al,keqka,kka,al,mvka,ali,mkai,alth,mkath,dal,mkasupp,mka &! intput
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,cplprec &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2,alx,flx_ka,flx_al,omega_ka,mkax &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,k1si,k2si,keqcca,authig,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3,pco2i,khco2i,ucv,dgasc,daqc,pco2th,k1al,k2al,k3al,k4al,keqka,mvka,ali,mkai,alth,mkath,dal
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz,kcca,pco2,torg,resp,kka,al,mkasupp,mka &
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx,pco2x,alx,mkax
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc,omega_cca,omega_ka
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,flx_co2,flx_ka,flx_al
integer,intent(inout)::iter,it
logical,intent(in)::cplprec
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp3 = 11
integer iz,row,nmx,ie,ie2,isp,iflx
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg,dprodsi,domega_fo_dpro,domega_ab_dpro,domega_an_dpro,domega_cc_dpro &
    & ,domega_cca_dmg,domega_cca_dca,domega_cca_dna,domega_cca_dsi,domega_cca_dpro,dprodpco2,domega_fo_dpco2,domega_ab_dpco2 &
    & ,domega_an_dpco2,domega_cc_dpco2,domega_cca_dpco2,dkhco2_dpro,dkhco2_dpco2,dpreccc_dpco2,preccc,khco2,khco2x,alpha,alphaprev &
    & ,dalpha_dpco2,edif,dedif_dpco2,dalpha_dna,dalpha_dsi,dalpha_dmg,dalpha_dca,dedif_dca,dedif_dmg,dedif_dsi,dedif_dna &
    & ,dkhco2_dca,dkhco2_dna,dkhco2_dsi,dkhco2_dmg,dpreccc_dca,dpreccc_dna,dpreccc_dsi,dpreccc_dmg,dpreccc_dmcc,dprodall &
    & ,domega_fo_dal,domega_ka_dal,domega_ka_dsi,domega_ka_dpco2,domega_ka_dpro,domega_ka_dna,domega_ka_dca,domega_ka_dmg  &
    & ,domega_ab_dal,domega_an_dal,domega_cc_dal,domega_cca_dal,dalpha_dal,dkhco2_dal,dedif_dal,dpreccc_dal,dprodal
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,st_cca,edifi,khco2n_tmp,pco2n_tmp,edifn_tmp,st_ka,caqsupp_tmp
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
real(kind=8),dimension(nz)::disonly ! for cc [1---yes, 0---no]
real(kind=8),dimension(nz)::disonly_ka ! for ka [1---yes, 0---no]
! real(kind=8)::disonly = 1d0 ! for cc 
! real(kind=8)::authig = 0d0 ! for cc whether to include authigenesis [1---yes, 0---no]
! real(kind=8)::authig = 1d0 ! 

integer,parameter :: iter_max = 300

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

    ! call calc_pH_v3( &
        ! & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
else 

    ! call calc_pH_v3( &
        ! & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
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
    flx_al = 0d0
    flx_cc = 0d0
    flx_ka = 0d0
    flx_co2 = 0d0

    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax,mgx+dconc,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,dprodmg &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax+dconc,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,alx,k1al,k2al,k3al,k4al &! input 
        & ,dprodna &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax+dconc,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,dprodca &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six+dconc,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,dprodsi &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx+dconc,k1al,k2al,k3al,k4al &! input
        & ,dprodal &! output
        & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x+dconc,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,dprodpco2 &! output
        & ) 
    dprodna = (dprodna-prox)/dconc
    dprodmg = (dprodmg-prox)/dconc
    dprodsi = (dprodsi-prox)/dconc
    dprodca = (dprodca-prox)/dconc
    dprodal = (dprodal-prox)/dconc
    dprodpco2 = (dprodpco2-prox)/dconc
    
    ! Fo + 4H+ = 2Mg2+ + SiO2(aq) + 2H2O 
    
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,omega_fo &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &!input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'fo' &! input 
        & ,domega_fo_dpco2 &! output
        & )
    domega_fo_dmg = (domega_fo_dmg-omega_fo)/dconc
    domega_fo_dsi = (domega_fo_dsi-omega_fo)/dconc
    domega_fo_dna = (domega_fo_dna-omega_fo)/dconc
    domega_fo_dca = (domega_fo_dca-omega_fo)/dconc
    domega_fo_dal = (domega_fo_dal-omega_fo)/dconc
    domega_fo_dpro = (domega_fo_dpro-omega_fo)/dconc
    domega_fo_dpco2 = (domega_fo_dpco2-omega_fo)/dconc
    
    domega_fo_dmg = domega_fo_dmg + domega_fo_dpro*dprodmg
    domega_fo_dca = domega_fo_dca + domega_fo_dpro*dprodca
    domega_fo_dna = domega_fo_dna + domega_fo_dpro*dprodna
    domega_fo_dsi = domega_fo_dsi + domega_fo_dpro*dprodsi
    domega_fo_dal = domega_fo_dal + domega_fo_dpro*dprodal
    domega_fo_dpco2 = domega_fo_dpco2 + domega_fo_dpro*dprodpco2
    
    ! Ab + H+ + 0.5H2O --> 0.5 kaolinite + Na+ + 2 SiO2(aq)
    
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,omega_ab &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ab' &! input 
        & ,domega_ab_dpco2 &! output
        & )
    domega_ab_dmg = (domega_ab_dmg-omega_ab)/dconc
    domega_ab_dsi = (domega_ab_dsi-omega_ab)/dconc
    domega_ab_dna = (domega_ab_dna-omega_ab)/dconc
    domega_ab_dca = (domega_ab_dca-omega_ab)/dconc
    domega_ab_dal = (domega_ab_dal-omega_ab)/dconc
    domega_ab_dpro = (domega_ab_dpro-omega_ab)/dconc
    domega_ab_dpco2 = (domega_ab_dpco2-omega_ab)/dconc
    
    domega_ab_dmg = domega_ab_dmg + domega_ab_dpro*dprodmg
    domega_ab_dca = domega_ab_dca + domega_ab_dpro*dprodca
    domega_ab_dna = domega_ab_dna + domega_ab_dpro*dprodna
    domega_ab_dsi = domega_ab_dsi + domega_ab_dpro*dprodsi
    domega_ab_dal = domega_ab_dal + domega_ab_dpro*dprodal
    domega_ab_dpco2 = domega_ab_dpco2 + domega_ab_dpro*dprodpco2
    
    ! An + 2H+ + H2O = kaolinite + Ca2+ 
    
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,omega_an &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'an' &! input 
        & ,domega_an_dpco2 &! output
        & )
    domega_an_dmg = (domega_an_dmg-omega_an)/dconc
    domega_an_dsi = (domega_an_dsi-omega_an)/dconc
    domega_an_dna = (domega_an_dna-omega_an)/dconc
    domega_an_dca = (domega_an_dca-omega_an)/dconc
    domega_an_dal = (domega_an_dal-omega_an)/dconc
    domega_an_dpro = (domega_an_dpro-omega_an)/dconc
    domega_an_dpco2 = (domega_an_dpco2-omega_an)/dconc
    
    domega_an_dmg = domega_an_dmg + domega_an_dpro*dprodmg
    domega_an_dca = domega_an_dca + domega_an_dpro*dprodca
    domega_an_dna = domega_an_dna + domega_an_dpro*dprodna
    domega_an_dsi = domega_an_dsi + domega_an_dpro*dprodsi
    domega_an_dal = domega_an_dal + domega_an_dpro*dprodal
    domega_an_dpco2 = domega_an_dpco2 + domega_an_dpro*dprodpco2
    
    ! Kaolinite + 6H+ = 2Al3+ + 2SiO2 
    
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,omega_ka &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'ka' &! input 
        & ,domega_ka_dpco2 &! output
        & )
    domega_ka_dmg = (domega_ka_dmg-omega_ka)/dconc
    domega_ka_dsi = (domega_ka_dsi-omega_ka)/dconc
    domega_ka_dna = (domega_ka_dna-omega_ka)/dconc
    domega_ka_dca = (domega_ka_dca-omega_ka)/dconc
    domega_ka_dal = (domega_ka_dal-omega_ka)/dconc
    domega_ka_dpro = (domega_ka_dpro-omega_ka)/dconc
    domega_ka_dpco2 = (domega_ka_dpco2-omega_ka)/dconc
    
    domega_ka_dmg = domega_ka_dmg + domega_ka_dpro*dprodmg
    domega_ka_dca = domega_ka_dca + domega_ka_dpro*dprodca
    domega_ka_dna = domega_ka_dna + domega_ka_dpro*dprodna
    domega_ka_dsi = domega_ka_dsi + domega_ka_dpro*dprodsi
    domega_ka_dal = domega_ka_dal + domega_ka_dpro*dprodal
    domega_ka_dpco2 = domega_ka_dpco2 + domega_ka_dpro*dprodpco2
    
    ! Cc = Ca2+ + CO32- 
    
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,omega_cc &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cc_dpco2 &! output
        & )
    domega_cc_dmg = (domega_cc_dmg-omega_cc)/dconc
    domega_cc_dsi = (domega_cc_dsi-omega_cc)/dconc
    domega_cc_dna = (domega_cc_dna-omega_cc)/dconc
    domega_cc_dca = (domega_cc_dca-omega_cc)/dconc
    domega_cc_dal = (domega_cc_dal-omega_cc)/dconc
    domega_cc_dpro = (domega_cc_dpro-omega_cc)/dconc
    domega_cc_dpco2 = (domega_cc_dpco2-omega_cc)/dconc
    
    domega_cc_dmg = domega_cc_dmg + domega_cc_dpro*dprodmg
    domega_cc_dca = domega_cc_dca + domega_cc_dpro*dprodca
    domega_cc_dna = domega_cc_dna + domega_cc_dpro*dprodna
    domega_cc_dsi = domega_cc_dsi + domega_cc_dpro*dprodsi
    domega_cc_dal = domega_cc_dal + domega_cc_dpro*dprodal
    domega_cc_dpco2 = domega_cc_dpco2 + domega_cc_dpro*dprodpco2
    
    disonly = 0d0
    disonly_ka = 0d0
    if (.not.cplprec)then
        disonly = 1d0
        disonly_ka = 1d0
    endif 
    
    ! authigenesis (prety much tentative parameterization)
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,omega_cca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx+dconc,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dmg &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six+dconc,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dsi &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,pco2x,cax,mgx,six,nax+dconc,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dna &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax+dconc,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dca &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox,alx+dconc,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dal &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x,cax,mgx,six,nax,prox+dconc,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dpro &! output
        & )
    call calc_omega_v2( &
        & nz,keqfo,keqab,keqan,keqcca,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input
        & ,pco2x+dconc,cax,mgx,six,nax,prox,alx,k1al,k2al,k3al,k4al,keqka,'cc' &! input 
        & ,domega_cca_dpco2 &! output
        & )
    domega_cca_dmg = (domega_cca_dmg-omega_cca)/dconc
    domega_cca_dsi = (domega_cca_dsi-omega_cca)/dconc
    domega_cca_dna = (domega_cca_dna-omega_cca)/dconc
    domega_cca_dca = (domega_cca_dca-omega_cca)/dconc
    domega_cca_dal = (domega_cca_dal-omega_cca)/dconc
    domega_cca_dpro = (domega_cca_dpro-omega_cca)/dconc
    domega_cca_dpco2 = (domega_cca_dpco2-omega_cca)/dconc
    
    domega_cca_dmg = domega_cca_dmg + domega_cca_dpro*dprodmg
    domega_cca_dca = domega_cca_dca + domega_cca_dpro*dprodca
    domega_cca_dna = domega_cca_dna + domega_cca_dpro*dprodna
    domega_cca_dsi = domega_cca_dsi + domega_cca_dpro*dprodsi
    domega_cca_dal = domega_cca_dal + domega_cca_dpro*dprodal
    domega_cca_dpco2 = domega_cca_dpco2 + domega_cca_dpro*dprodpco2

    do iz = 1, nz  !================================
        
        do isp = 1, 5
        
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
                ! disonly(iz) = merge(1d0,0d0,pco2x(iz)<pco2th)
                omega_tmp_th = omega_tmp*disonly(iz)
                m_tmp = mccx(iz)
                mth_tmp = mccth 
                mi_tmp = mcci
                mp_tmp = mccx(min(nz,iz+1))
                msupp_tmp = mccsupp(iz)
                mprev_tmp = mcc(iz)
            elseif (isp==5)then ! Ka
                k_tmp = kka(iz)
                mv_tmp = mvka
                omega_tmp = omega_ka(iz)
                ! disonly_ka(iz) = merge(1d0,0d0,alx(iz)<alth)
                omega_tmp_th = omega_tmp*disonly_ka(iz)
                m_tmp = mkax(iz)
                mth_tmp = mkath 
                mi_tmp = mkai
                mp_tmp = mkax(min(nz,iz+1))
                msupp_tmp = mkasupp(iz)
                mprev_tmp = mka(iz)
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
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 8 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 9 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)

                amx3(row,row + 10 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 8 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 9 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 8 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dmg(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dsi(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dna(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dca(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dal(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 7 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_cc_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & + kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dpco2(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                
            elseif (isp==5) then 
                amx3(row,row + 1 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 2 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 3 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 4 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 5 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                amx3(row,row + 6 ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_ka_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *pco2x(iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
                    
                flx_ka(itflx,iz) = (&
                    & (m_tmp-mprev_tmp)/dt &
                    & )
                flx_ka(iadv,iz) = (&
                    & -w*(mp_tmp-m_tmp)/dz(iz)  &
                    & )
                flx_ka(irxn_fo,iz) = (&
                    & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & )
                flx_ka(irain,iz) = (&
                    & -msupp_tmp  &
                    & )
                flx_ka(ires,iz) = sum(flx_ka(:,iz))
                if (isnan(flx_ka(ires,iz))) then 
                    print *,'ka',iz,(flx_ka(iflx,iz),iflx=1,nflx)
                endif 
            endif 
        enddo 
    end do  !================================

    do iz = 1, nz
        
        do isp = 1, 5

            row = nsp3*(iz-1)+5 + isp
            
            if (isp==1) then ! mg 
                d_tmp = dmg
                caq_tmp = mgx(iz)
                caq_tmp_prev = mg(iz)
                caq_tmp_p = mgx(min(nz,iz+1))
                caq_tmp_n = mgx(max(1,iz-1))
                caqth_tmp = mgth
                caqi_tmp = mgi
                caqsupp_tmp = mgsupp(iz)
                st_fo = 2d0
                st_ab = 0d0
                st_an = 0d0
                st_cc = 0d0
                st_cca = 0d0
                st_ka = 0d0
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
                caqsupp_tmp = sisupp(iz)
                st_fo = 1d0
                st_ab = 3d0
                st_an = 2d0
                st_cc = 0d0
                st_cca = 0d0
                st_ka = 2d0
                rxn_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(1d0-omega_fo(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
                    & + st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(1d0-omega_ka(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0) 
                drxndisp_tmp = st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0) &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
                    & + st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dsi(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0) 
            elseif (isp==3) then  ! na
                d_tmp = dna
                caq_tmp = nax(iz)
                caq_tmp_prev = na(iz)
                caq_tmp_p = nax(min(nz,iz+1))
                caq_tmp_n = nax(max(1,iz-1))
                caqth_tmp = nath
                caqi_tmp = nai
                caqsupp_tmp = nasupp(iz)
                st_fo = 0d0
                st_ab = 1d0
                st_an = 0d0
                st_cc = 0d0
                st_cca = 0d0
                st_ka = 0d0
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
                caqsupp_tmp = casupp(iz)
                st_fo = 0d0
                st_ab = 0d0
                st_an = 1d0
                st_cc = 1d0
                st_cca = 1d0
                st_ka = 0d0
                rxn_tmp = st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(1d0-omega_cc(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0) &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(1d0-omega_cca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) 
                drxndisp_tmp = ( &
                    & st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dca(iz)) &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0) &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dca(iz))*authig &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & )
            elseif (isp==5) then  ! al
                d_tmp = dal
                caq_tmp = alx(iz)
                caq_tmp_prev = al(iz)
                caq_tmp_p = alx(min(nz,iz+1))
                caq_tmp_n = alx(max(1,iz-1))
                caqth_tmp = alth
                caqi_tmp = ali
                caqsupp_tmp = alsupp(iz)
                st_fo = 0d0
                st_ab = 1d0
                st_an = 2d0
                st_cc = 0d0
                st_cca = 0d0
                st_ka = 2d0
                rxn_tmp = &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(1d0-omega_ab(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
                    & + st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(1d0-omega_an(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(1d0-omega_ka(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0) 
                drxndisp_tmp =  &
                    & + st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dal(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0) &
                    & + st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dal(iz)) &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0) &
                    & + st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dal(iz)) &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0) 
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
                & - caqsupp_tmp*dt &
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
            
            amx3(row,row  - isp - 4) = (     & 
                & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*1d0*(1d0-omega_fo(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                & ) &
                & *mfox(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -3) = (     & 
                & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*1d0*(1d0-omega_ab(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                & ) &
                & *mabx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -2) = (     & 
                & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*1d0*(1d0-omega_an(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                & ) &
                & *manx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp -1) = (     & 
                & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*1d0*(1d0-omega_cc(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                & ) &
                & *mccx(iz) &
                & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
            amx3(row,row  - isp ) = (     & 
                & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*1d0*(1d0-omega_ka(iz))*dt &
                & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                & ) &
                & *mkax(iz) &
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
            
                amx3(row,row  + 4) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 5) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
                flx_mg(irain,iz) = (&
                    & - caqsupp_tmp &
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
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 3) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 4) = (     & 
                    & - st_fo*kfo(iz)*poro(iz)*hr(iz)*mvfo*1d-6*mfox(iz)*(-domega_fo_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_fo(iz) < 0d0)  &
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
                flx_si(irain,iz) = (&
                    & - caqsupp_tmp &
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
            
                amx3(row,row  + 2) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 3) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
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
                flx_na(irain,iz) = (&
                    & - caqsupp_tmp &
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
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dmg(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dsi(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dna(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dal(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dal(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *alx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 2) = (     & 
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_cc*kcc(iz)*poro(iz)*hr(iz)*mvcc*1d-6*mccx(iz)*(-domega_cc_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_cc(iz)*disonly(iz) < 0d0)  &
                    & + st_cca*kcca(iz)*poro(iz)*sat(iz)*(-domega_cca_dpco2(iz))*authig*dt &
                    & *merge(0d0,1d0,1d0-omega_cca(iz) > 0d0) &
                    & ) &
                    & *pco2x(iz) &
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
                flx_ca(irain,iz) = (&
                    & - caqsupp_tmp &
                    & ) 
                flx_ca(ires,iz) = sum(flx_ca(:,iz))
                if (isnan(flx_ca(ires,iz))) then 
                    print *,'ca',iz,(flx_ca(iflx,iz),iflx=1,nflx)
                endif 
            elseif (isp==5) then 
            
                amx3(row,row  - 4) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dmg(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *mgx(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 3) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dsi(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *six(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 2) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dna(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *nax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  - 1) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dca(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *cax(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            
                amx3(row,row  + 1) = (     & 
                    & - st_ab*kab(iz)*poro(iz)*hr(iz)*mvab*1d-6*mabx(iz)*(-domega_ab_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ab(iz) < 0d0)  &
                    & - st_an*kan(iz)*poro(iz)*hr(iz)*mvan*1d-6*manx(iz)*(-domega_an_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_an(iz) < 0d0)  &
                    & - st_ka*kka(iz)*poro(iz)*hr(iz)*mvka*1d-6*mkax(iz)*(-domega_ka_dpco2(iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_ka(iz)*disonly_ka(iz) < 0d0)  &
                    & ) &
                    & *pco2x(iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                    
                flx_al(itflx,iz) = (&
                    & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)/dt  &
                    & ) 
                flx_al(iadv,iz) = (&
                    & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz) &
                    & ) 
                flx_al(idif,iz) = (&
                    & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz) &
                    & ) 
                flx_al(irxn_fo,iz) = (&
                    & - rxn_tmp &
                    & ) 
                flx_al(irain,iz) = (&
                    & - caqsupp_tmp &
                    & ) 
                flx_al(ires,iz) = sum(flx_al(:,iz))
                if (isnan(flx_al(ires,iz))) then 
                    print *,'si',iz,(flx_al(iflx,iz),iflx=1,nflx)
                endif 
            endif 
            
            ! amx3(row,:) = amx3(row,:)/(poro(iz)*sat(iz)*1d3)
            ! ymx3(row) = ymx3(row)/(poro(iz)*sat(iz)*1d3)
        
        enddo 
        
    end do  ! ==============================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    khco2 = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
    khco2x = kco2*(1d0+k1/prox + k1*k2/prox/prox)
    
    disonly = merge(1d0,0d0,pco2x<pco2th)
    preccc = kcc*poro*hr*mvcc*1d-6*mccx*(1d0-omega_cc) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
        
    dpreccc_dpco2=0d0
    dkhco2_dpco2 = 0d0
        
    dpreccc_dpco2=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dpco2) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dna=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dna) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dmg=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dmg) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dsi=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dsi) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dca=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dca) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dal=kcc*poro*hr*mvcc*1d-6*mccx*(-domega_cc_dal) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dpreccc_dmcc=kcc*poro*hr*mvcc*1d-6*(1d0-omega_cc) &
        & *merge(0d0,1d0,1d0-omega_cc*disonly < 0d0)
    dkhco2_dpro = kco2*(k1*(-1d0)/prox**2d0 + k1*k2*(-2d0)/prox**3d0)
    dkhco2_dpco2 = dkhco2_dpro*dprodpco2
    dkhco2_dna = dkhco2_dpro*dprodna
    dkhco2_dmg = dkhco2_dpro*dprodmg
    dkhco2_dsi = dkhco2_dpro*dprodsi
    dkhco2_dca = dkhco2_dpro*dprodca
    dkhco2_dal = dkhco2_dpro*dprodal
    
    edif = ucv*poro*(1.0d0-sat)*1d3*torg*dgasc+poro*sat*khco2x*1d3*tora*daqc
    edifi = edif(1)
    edifi = ucv*1d3*dgasc 
    
    dedif_dpco2 = poro*sat*dkhco2_dpco2*1d3*tora*daqc
    dedif_dna = poro*sat*dkhco2_dna*1d3*tora*daqc
    dedif_dmg = poro*sat*dkhco2_dmg*1d3*tora*daqc
    dedif_dsi = poro*sat*dkhco2_dsi*1d3*tora*daqc
    dedif_dca = poro*sat*dkhco2_dca*1d3*tora*daqc
    dedif_dal = poro*sat*dkhco2_dal*1d3*tora*daqc
    
    alpha = ucv*poro*(1.0d0-sat)*1d3+poro*sat*khco2x*1d3
    alphaprev = ucv*poroprev*(1.0d0-sat)*1d3+poroprev*sat*khco2*1d3
    
    dalpha_dpco2 = poro*sat*dkhco2_dpco2*1d3
    dalpha_dna = poro*sat*dkhco2_dna*1d3
    dalpha_dmg = poro*sat*dkhco2_dmg*1d3
    dalpha_dsi = poro*sat*dkhco2_dsi*1d3
    dalpha_dca = poro*sat*dkhco2_dca*1d3
    dalpha_dal = poro*sat*dkhco2_dal*1d3
    
    if (any(isnan(alpha)).or.any(isnan(alphaprev))) then 
        print *, alpha 
        print *, alphaprev
        print *, poroprev
        print *, khco2
        print *, pro
        stop
    endif 
    
    do iz = 1, nz

        row = nsp3*(iz-1) + 11
        
        pco2n_tmp = pco2x(max(1,iz-1))
        khco2n_tmp = khco2x(max(1,iz-1))
        edifn_tmp = edif(max(1,iz-1))
        if (iz == 1) then 
            pco2n_tmp = pco2i
            khco2n_tmp = khco2i
            edifn_tmp = edifi
        endif 

        amx3(row,row) = ( &
            & (alpha(iz) + dalpha_dpco2(iz)*pco2x(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*merge(0d0,-1d0,iz==nz)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & +0.5d0*(dedif_dpco2(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
            & - 0.5d0*(dedif_dpco2(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*1d0)/dz(iz)*dt &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dpco2(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dpco2(iz)*dt &
            & ) &
            & *merge(1.0d0,pco2x(iz),pco2x(iz)<pco2th)

        ymx3(row) = ( &
            & (alpha(iz)*pco2x(iz)-alphaprev(iz)*pco2(iz)) &
            & -( 0.5d0*(edif(iz)+edif(min(nz,iz+1)))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(edif(iz)+edifn_tmp)*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(khco2x(iz)*pco2x(iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*dt &
            & -resp(iz)*dt &
            & -preccc(iz)*dt &
            & -pco2supp(iz)*dt &
            & ) &
            & *merge(0.0d0,1.0d0,pco2x(iz)<pco2th)
        
        
        if (iz/=nz) then 
        
            amx3(row,row+nsp3) = ( &
                    & -( 0.5d0*(edif(iz)+edif(iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    & + 0.5d0*(dedif_dpco2(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,pco2x(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-5) = ( &
                    & -( 0.5d0*(dedif_dmg(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,mgx(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-4) = ( &
                    & -( 0.5d0*(dedif_dsi(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,six(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-3) = ( &
                    & -( 0.5d0*(dedif_dna(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,nax(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-2) = ( &
                    & -( 0.5d0*(dedif_dca(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,cax(iz+1),pco2x(iz)<pco2th)
        
            amx3(row,row+nsp3-1) = ( &
                    & -( 0.5d0*(dedif_dal(iz+1))*(pco2x(iz+1)-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                    & ) &
                    & *merge(0.0d0,alx(iz+1),pco2x(iz)<pco2th)
        
        endif 
        
        if (iz/=1) then 

            amx3(row,row-nsp3) = ( &
                & -(- 0.5d0*(edif(iz)+edif(iz-1))*(-1d0)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) &
                & - 0.5d0*(dedif_dpco2(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-khco2x(iz-1)*1d0)/dz(iz)*dt &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dpco2(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,pco2x(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-5) = ( &
                & -(- 0.5d0*(dedif_dmg(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dmg(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,mgx(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-4) = ( &
                & -(- 0.5d0*(dedif_dsi(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dsi(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,six(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-3) = ( &
                & -(- 0.5d0*(dedif_dna(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dna(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,nax(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-2) = ( &
                & -(- 0.5d0*(dedif_dca(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dca(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,cax(iz-1),pco2x(iz)<pco2th)

            amx3(row,row-nsp3-1) = ( &
                & -(- 0.5d0*(dedif_dal(iz-1))*(pco2x(iz)-pco2x(iz-1))/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhco2_dal(iz-1)*pco2x(iz-1))/dz(iz)*dt &
                & ) &
                & *merge(0.0d0,alx(iz-1),pco2x(iz)<pco2th)
        endif 
        
        amx3(row,row-7) = ( &
            & -dpreccc_dmcc(iz)*dt &
            & ) &
            & *merge(1.0d0,mccx(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-5) = ( &
            & (dalpha_dmg(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dmg(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dmg(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dmg(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dmg(iz)*dt &
            & ) &
            & *merge(1.0d0,mgx(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-4) = ( &
            & (dalpha_dsi(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dsi(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dsi(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dsi(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dsi(iz)*dt &
            & ) &
            & *merge(1.0d0,six(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-3) = ( &
            & (dalpha_dna(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dna(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dna(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dna(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dna(iz)*dt &
            & ) &
            & *merge(1.0d0,nax(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-2) = ( &
            & (dalpha_dca(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dca(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dca(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dca(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dca(iz)*dt &
            & ) &
            & *merge(1.0d0,cax(iz),pco2x(iz)<pco2th)
        
        amx3(row,row-1) = ( &
            & (dalpha_dal(iz)*pco2x(iz)) &
            & -( 0.5d0*(dedif_dal(iz))*(pco2x(min(nz,iz+1))-pco2x(iz))/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
            & - 0.5d0*(dedif_dal(iz))*(pco2x(iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
            & +poro(iz)*sat(iz)*v(iz)*1d3*(dkhco2_dal(iz)*pco2x(iz))/dz(iz)*dt &
            & -dpreccc_dal(iz)*dt &
            & ) &
            & *merge(1.0d0,alx(iz),pco2x(iz)<pco2th)
        
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
        flx_co2(irxn_fo,iz) = -resp(iz)
        flx_co2(irain,iz) = -preccc(iz) - pco2supp(iz)
        flx_co2(ires,iz) = sum(flx_co2(:,iz))
        
        if (any(isnan(flx_co2(:,iz)))) then
            print *,flx_co2(:,iz)
        endif 
        
        ! amx3(row,:) = amx3(row,:)/alpha(iz)
        ! ymx3(row) = ymx3(row)/alpha(iz)

    end do 
    
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
            print *,'nan at', iz,z(iz),'Ka'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            mkax(iz) = mkax(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            mkax(iz) = mkax(iz)*0.50d0
        else   
            mkax(iz) = mkax(iz)*exp(ymx3(row))
        endif
        
        row = 6 + nsp3*(iz-1)

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
        
        row = 7 + nsp3*(iz-1)

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
        
        row = 8 + nsp3*(iz-1)

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
        
        row = 9 + nsp3*(iz-1)

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
        
        row = 10 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'al'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                alx(iz) = 0.1d0*alth
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            alx(iz) = alx(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            alx(iz) = alx(iz)*0.50d0
        else
            alx(iz) = alx(iz)*exp(ymx3(row))
        endif
        
        row = 11 + nsp3*(iz-1)

        if (isnan(ymx3(row))) then 
            print *,'nan at', iz,z(iz),'pco2'
            if (z(iz)<zrxn) then 
                ymx3(row)=0d0
                pco2x(iz) = 0.1d0*pco2th
            endif
        endif

        if ((.not.isnan(ymx3(row))).and.ymx3(row) >threshold) then 
            pco2x(iz) = pco2x(iz)*1.5d0
        else if (ymx3(row) < -threshold) then 
            pco2x(iz) = pco2x(iz)*0.50d0
        else
            pco2x(iz) = pco2x(iz)*exp(ymx3(row))
        endif
        
    end do 

    error = maxval(exp(abs(ymx3))) - 1.0d0

    if (isnan(error).or.info/=0 .or. any(isnan(mgx)) .or. any(isnan(six)).or. any(isnan(cax)).or. any(isnan(pco2x)) &
        & .or. any(isnan(nax)) .or. any(isnan(mfox)).or. any(isnan(mabx)).or. any(isnan(manx)).or. any(isnan(mccx)) &
        & .or. any(isnan(alx)) .or. any(isnan(mkax))) then 
        error = 1d3
        print *, '!! error is NaN; values are returned to those before iteration with reducing dt'
        print*, 'isnan(error), info/=0,any(isnan(mgx)),any(isnan(mfox)))'
        print*,isnan(error),info/=0,any(isnan(mgx)),any(isnan(six)),any(isnan(mfox)),any(isnan(nax)),any(isnan(mabx)) &
            & ,any(isnan(cax)),any(isnan(manx)),any(isnan(mccx)),any(isnan(pco2x)),any(isnan(mkax)),any(isnan(alx))
        stop
        mgx = mg
        six = si
        nax = na
        cax = ca
        alx = al
        mfox = mfo
        mabx = mab
        manx = man
        mccx = mcc
        mkax = mka
        prox = pro
        pco2x = pco2
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
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
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

        if (mkax(iz) < 0.0d0) then
            mkax(iz) = mkax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 6 + nsp3*(iz-1)

        if (mgx(iz) < 0.0d0) then
            mgx(iz) = mgx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 7 + nsp3*(iz-1)

        if (six(iz) < 0.0d0) then
            six(iz) = six(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 8 + nsp3*(iz-1)

        if (nax(iz) < 0.0d0) then
            nax(iz) = nax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 9 + nsp3*(iz-1)

        if (cax(iz) < 0.0d0) then
            cax(iz) = cax(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 10 + nsp3*(iz-1)

        if (alx(iz) < 0.0d0) then
            alx(iz) = alx(iz)/exp(ymx3(row))*0.5d0
            error = 1.0d0
        end if
        
        row = 11 + nsp3*(iz-1)

        if (pco2x(iz) < 0.0d0) then
            pco2x(iz) = pco2x(iz)/exp(ymx3(row))*0.5d0
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
    print *,'-=-=-=-=-=-= Mg, Si, Na, Ca, Fo, Ab, An, pCO2 -=-=-=-=-=-=-='
    print *, 'mg:', (mgx(iz),iz=1,nz, nz/5)
    print *, 'si:', (six(iz),iz=1,nz, nz/5)
    print *, 'na:', (nax(iz),iz=1,nz, nz/5)
    print *, 'ca:', (cax(iz),iz=1,nz, nz/5)
    print *, 'al:', (alx(iz),iz=1,nz, nz/5)
    print *, 'fo:', (mfox(iz),iz=1,nz, nz/5)
    print *, 'ab:', (mabx(iz),iz=1,nz, nz/5)
    print *, 'an:', (manx(iz),iz=1,nz, nz/5)
    print *, 'ka:', (mkax(iz),iz=1,nz, nz/5)
    print *, 'cc:', (mccx(iz),iz=1,nz, nz/5)
    print *, 'omega_fo:', (omega_fo(iz),iz=1,nz, nz/5)
    print *, 'omega_ab:', (omega_ab(iz),iz=1,nz, nz/5)
    print *, 'omega_an:', (omega_an(iz),iz=1,nz, nz/5)
    print *, 'omega_ka:', (omega_ka(iz),iz=1,nz, nz/5)
    print *, 'omega_cc:', (omega_cc(iz),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
    print *,'-=-=-=-=-=-= CO2 -=-=-=-=-=-=-='
    print *, 'co2:', (pco2x(iz),iz=1,nz, nz/5)
    print *
#endif     
enddo

endsubroutine alsilicate_dis_co2_1D_v2














!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine alsilicate_dis_co2_1D_v3( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
    & ,pco2th,resp,k1al,k2al,k3al,k4al,keqka,kka,al,mvka,ali,mkai,alth,mkath,dal,mkasupp,mka &! intput
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,cplprec &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2,alx,flx_ka,flx_al,omega_ka,mkax &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,k1si,k2si,keqcca,authig,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3,pco2i,khco2i,ucv,dgasc,daqc,pco2th,k1al,k2al,k3al,k4al,keqka,mvka,ali,mkai,alth,mkath,dal
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz,kcca,pco2,torg,resp,kka,al,mkasupp,mka &
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx,pco2x,alx,mkax
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc,omega_cca,omega_ka
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,flx_co2,flx_ka,flx_al
integer,intent(inout)::iter,it
logical,intent(in)::cplprec
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp_sld = 5
integer,parameter::nsp_sld_2 = 2
integer,parameter::nsp_aq = 5
integer,parameter::nsp_aq_cnst = 4
integer,parameter::nsp_aq_ph = 5
integer,parameter::nsp_gas_ph = 1
integer,parameter::nsp_aq_all = 9
integer,parameter::nsp_gas = 1
integer,parameter::nsp3 = nsp_sld + nsp_aq + nsp_gas
character(5),dimension(nsp_sld)::chrsld
character(5),dimension(nsp_sld_2)::chrsld_2
character(5),dimension(nsp_aq)::chraq
character(5),dimension(nsp_aq_cnst)::chraq_cnst
character(5),dimension(nsp_aq_ph)::chraq_ph
character(5),dimension(nsp_gas_ph)::chrgas_ph
character(5),dimension(nsp_aq_all)::chraq_all
character(5),dimension(nsp_gas)::chrgas
real(kind=8),dimension(nsp_sld)::msldi,msldth,keq,mv
real(kind=8),dimension(nsp_aq)::maqi,maqth,daq,dmaq 
real(kind=8),dimension(nsp_gas)::mgasi,mgasth,dgasa,dgasg,dmgas,khgasi,dgasi
real(kind=8),dimension(nsp_sld,nsp_aq)::staq
real(kind=8),dimension(nsp_sld,nsp_gas)::stgas
real(kind=8),dimension(nsp_sld,nz)::msldx,msld,ksld,omega,nonprec,rxnsld,msldsupp,domega_dpro 
real(kind=8),dimension(nsp_sld,nsp3,nz)::domega_disp,drxnsld_disp 
real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq,drxnsld_dmaq 
real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas,drxnsld_dmgas 
real(kind=8),dimension(nsp_sld,nflx,nz)::flx_sld
real(kind=8),dimension(nsp_aq_cnst,nz)::maqcnst
real(kind=8),dimension(nsp_aq,nz)::maqx,maq,rxnaq,maqsupp,drxnaq_dpro,dprodmaq 
real(kind=8),dimension(nsp_aq,nsp3,nz)::drxnaq_disp
real(kind=8),dimension(nsp_aq,nflx,nz)::flx_aq
real(kind=8),dimension(nsp_gas,nz)::mgasx,mgas,khgasx,khgas,dgas,agasx,agas,rxngas,mgassupp,dkhgas_dpro,ddgas_dpro,dagas_dpro &
    & ,dprodmgas
real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
real(kind=8),dimension(nsp_gas,nflx,nz)::flx_gas 
integer iz,row,nmx,ie,ie2,isp,iflx,isps,ispa,ispg,ispa2,ispg2,col
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg,dprodsi,domega_fo_dpro,domega_ab_dpro,domega_an_dpro,domega_cc_dpro &
    & ,domega_cca_dmg,domega_cca_dca,domega_cca_dna,domega_cca_dsi,domega_cca_dpro,dprodpco2,domega_fo_dpco2,domega_ab_dpco2 &
    & ,domega_an_dpco2,domega_cc_dpco2,domega_cca_dpco2,dkhco2_dpro,dkhco2_dpco2,dpreccc_dpco2,preccc,khco2,khco2x,alpha,alphaprev &
    & ,dalpha_dpco2,edif,dedif_dpco2,dalpha_dna,dalpha_dsi,dalpha_dmg,dalpha_dca,dedif_dca,dedif_dmg,dedif_dsi,dedif_dna &
    & ,dkhco2_dca,dkhco2_dna,dkhco2_dsi,dkhco2_dmg,dpreccc_dca,dpreccc_dna,dpreccc_dsi,dpreccc_dmg,dpreccc_dmcc,dprodall &
    & ,domega_fo_dal,domega_ka_dal,domega_ka_dsi,domega_ka_dpco2,domega_ka_dpro,domega_ka_dna,domega_ka_dca,domega_ka_dmg  &
    & ,domega_ab_dal,domega_an_dal,domega_cc_dal,domega_cca_dal,dalpha_dal,dkhco2_dal,dedif_dal,dpreccc_dal,dprodal
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,st_cca,edifi,khco2n_tmp,pco2n_tmp,edifn_tmp,st_ka,caqsupp_tmp

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 100d0
real(kind=8),dimension(nz)::disonly ! for cc [1---yes, 0---no]
real(kind=8),dimension(nz)::disonly_ka ! for ka [1---yes, 0---no]
! real(kind=8)::disonly = 1d0 ! for cc 
! real(kind=8)::authig = 0d0 ! for cc whether to include authigenesis [1---yes, 0---no]
! real(kind=8)::authig = 1d0 ! 

integer,parameter :: iter_max = 300

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! passing major variables to more compact matrices 

chrsld = (/'fo','ab','an','cc','ka'/)
chraq = (/'mg','si','na','ca','al'/)
chrgas = (/'pco2'/)

! previous values
msld(findloc(chrsld,'fo',dim=1),:)=mfo(:)
msld(findloc(chrsld,'ab',dim=1),:)=mab(:)
msld(findloc(chrsld,'an',dim=1),:)=man(:)
msld(findloc(chrsld,'cc',dim=1),:)=mcc(:)
msld(findloc(chrsld,'ka',dim=1),:)=mka(:)

maq(findloc(chraq,'mg',dim=1),:)=mg(:)
maq(findloc(chraq,'si',dim=1),:)=si(:)
maq(findloc(chraq,'na',dim=1),:)=na(:)
maq(findloc(chraq,'ca',dim=1),:)=ca(:)
maq(findloc(chraq,'al',dim=1),:)=al(:)

mgas(findloc(chrgas,'pco2',dim=1),:)=pco2(:)

! output
msldx(findloc(chrsld,'fo',dim=1),:)=mfox(:)
msldx(findloc(chrsld,'ab',dim=1),:)=mabx(:)
msldx(findloc(chrsld,'an',dim=1),:)=manx(:)
msldx(findloc(chrsld,'cc',dim=1),:)=mccx(:)
msldx(findloc(chrsld,'ka',dim=1),:)=mkax(:)

maqx(findloc(chraq,'mg',dim=1),:)=mgx(:)
maqx(findloc(chraq,'si',dim=1),:)=six(:)
maqx(findloc(chraq,'na',dim=1),:)=nax(:)
maqx(findloc(chraq,'ca',dim=1),:)=cax(:)
maqx(findloc(chraq,'al',dim=1),:)=alx(:)

mgasx(findloc(chrgas,'pco2',dim=1),:)=pco2x(:)

! constants
msldi = (/mfoi,mabi,mani,mcci,mkai/)
msldth = (/mfoth,mabth,manth,mccth,mkath/)
maqi = (/mgi,sii,nai,cai,ali/)
maqth = (/mgth,sith,nath,cath,alth/)
mgasth = (/pco2th/)
mgasi = (/pco2i/)

ksld(findloc(chrsld,'fo',dim=1),:) = kfo
ksld(findloc(chrsld,'ab',dim=1),:) = kab
ksld(findloc(chrsld,'an',dim=1),:) = kan
ksld(findloc(chrsld,'cc',dim=1),:) = kcc
ksld(findloc(chrsld,'ka',dim=1),:) = kka

msldsupp(findloc(chrsld,'fo',dim=1),:) = mfosupp
msldsupp(findloc(chrsld,'ab',dim=1),:) = mabsupp
msldsupp(findloc(chrsld,'an',dim=1),:) = mansupp
msldsupp(findloc(chrsld,'cc',dim=1),:) = mccsupp
msldsupp(findloc(chrsld,'ka',dim=1),:) = mkasupp

maqsupp(findloc(chraq,'mg',dim=1),:) = mgsupp
maqsupp(findloc(chraq,'si',dim=1),:) = sisupp
maqsupp(findloc(chraq,'na',dim=1),:) = nasupp
maqsupp(findloc(chraq,'ca',dim=1),:) = casupp
maqsupp(findloc(chraq,'al',dim=1),:) = alsupp

mgassupp(findloc(chrgas,'pco2',dim=1),:) = pco2supp

daq = (/dmg,dsi,dna,dca,dal/)

dgasg = (/dgasc/)
dgasa = (/daqc/)

khgasi = (/khco2i/)

mv = (/mvfo,mvab,mvan,mvcc,mvka/)

chraq_cnst = (/'so4  ','fe2  ','fe3  ','k    '/)
maqcnst = 0d0

chrsld_2 = (/'cc','ka'/)

chraq_ph = chraq
chrgas_ph = chrgas

staq = 0d0
stgas = 0d0
! Forsterite; Mg2SiO4
staq(findloc(chrsld,'fo',dim=1), findloc(chraq,'mg',dim=1)) = 2d0
staq(findloc(chrsld,'fo',dim=1), findloc(chraq,'si',dim=1)) = 1d0
! Albite; NaAlSi3O8
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'na',dim=1)) = 1d0
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'si',dim=1)) = 3d0
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'al',dim=1)) = 1d0
! Anothite; CaAl2Si2O8
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'ca',dim=1)) = 1d0
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'si',dim=1)) = 2d0
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'al',dim=1)) = 2d0
! Calcite; CaCO3
staq(findloc(chrsld,'cc',dim=1), findloc(chraq,'ca',dim=1)) = 1d0
stgas(findloc(chrsld,'cc',dim=1), findloc(chrgas,'pco2',dim=1)) = 1d0
! Kaolinite; Al2Si2O5(OH)4
staq(findloc(chrsld,'ka',dim=1), findloc(chraq,'si',dim=1)) = 2d0
staq(findloc(chrsld,'ka',dim=1), findloc(chraq,'al',dim=1)) = 2d0

! print *, staq
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

    ! call calc_pH_v3( &
        ! & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
else 

    ! call calc_pH_v3( &
        ! & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,alx,k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
endif

! print *, 'starting silciate calculation'

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_sld = 0d0
    flx_aq = 0d0
    flx_gas = 0d0
    
    ! pH calculation and its derivative wrt aq and gas species
    
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
    do ispa=1,nsp_aq
        dmaq = 0d0
        dmgas = 0d0
        if (any (chraq_ph == chraq(ispa))) then 
            dmaq(ispa) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmaq(ispa,:) &! output
                & ) 
            dprodmaq(ispa,:) = (dprodmaq(ispa,:) - prox(:))/dconc
        endif 
    enddo 
    
    do ispg=1,nsp_gas
        dmaq = 0d0
        dmgas = 0d0
        if (any (chrgas_ph == chrgas(ispg))) then 
            dmgas(ispg) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmgas(ispg,:) &! output
                & ) 
            dprodmgas(ispg,:) = (dprodmgas(ispg,:) - prox(:))/dconc
        endif 
    enddo
    
    ! saturation state calc. and their derivatives wrt aq and gas species
    
    omega = 0d0
    domega_dpro = 0d0
    domega_dmaq = 0d0
    domega_dmgas = 0d0
    
    do isps =1, nsp_sld
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:) &
            & ,maqx(findloc(chraq,'ca',dim=1),:) &
            & ,maqx(findloc(chraq,'mg',dim=1),:) &
            & ,maqx(findloc(chraq,'si',dim=1),:) &
            & ,maqx(findloc(chraq,'na',dim=1),:) &
            & ,prox &
            & ,maqx(findloc(chraq,'al',dim=1),:) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,omega(isps,:) &! output
            & )
        dmaq = 0d0
        dmgas = 0d0
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
            & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
            & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
            & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
            & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
            & ,prox+dconc &
            & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,domega_dpro(isps,:) &! output
            & )
        domega_dpro(isps,:) = (domega_dpro(isps,:)-omega(isps,:))/dconc
        do ispa = 1, nsp_aq
            dmaq = 0d0
            dmgas = 0d0
            if (any (chraq_ph == chraq(ispa))) then 
                dmaq(ispa) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmaq(isps,ispa,:) &! output
                    & )
                domega_dmaq(isps,ispa,:) = (domega_dmaq(isps,ispa,:)-omega(isps,:))/dconc
                domega_dmaq(isps,ispa,:) = domega_dmaq(isps,ispa,:) + domega_dpro(isps,:)*dprodmaq(ispa,:)
            endif 
        enddo
        do ispg = 1, nsp_gas
            dmaq = 0d0
            dmgas = 0d0
            if (any (chrgas_ph == chrgas(ispg))) then 
                dmgas(ispg) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmgas(isps,ispg,:) &! output
                    & )
                domega_dmgas(isps,ispg,:) = (domega_dmgas(isps,ispg,:)-omega(isps,:))/dconc
                domega_dmgas(isps,ispg,:) = domega_dmgas(isps,ispg,:) + domega_dpro(isps,:)*dprodmgas(ispg,:)
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
            mprev_tmp = msld(isps,iz)
            
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
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1)+nsp_sld + ispa
                
                amx3(row,col ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmaq(isps,ispa,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *maqx(ispa,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            do ispg = 1, nsp_gas 
                col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                amx3(row,col) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmgas(isps,ispg,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            flx_sld(isps,itflx,iz) = ( &
                & (m_tmp-mprev_tmp)/dt &
                & )
            flx_sld(isps,iadv,iz) = ( &
                & -w*(mp_tmp-m_tmp)/dz(iz)  &
                & )
            flx_sld(isps,irxn_fo,iz) = ( &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & )
            flx_sld(isps,irain,iz) = (&
                & - msupp_tmp  &
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
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*dt &
                & - rxn_tmp*dt &
                & - caqsupp_tmp*dt &
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
            
            do isps = 1, nsp_sld
                col = nsp3*(iz-1)+ isps
                
                amx3(row, col) = (     & 
                    & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                    & ) &
                    & *msldx(isps,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            enddo 
            
            do ispa2 = 1, nsp_aq
                col = nsp3*(iz-1)+ nsp_sld + ispa2
                
                if (ispa2 == ispa) cycle
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmaq(isps,ispa2,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *maqx(ispa2,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                enddo 
            enddo 
            
            do ispg = 1, nsp_gas
                col = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmgas(isps,ispg,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *mgasx(ispg,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
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
            flx_aq(ispa,irxn_fo,iz) = (&
                & - rxn_tmp &
                & ) 
            flx_aq(ispa,irain,iz) = (&
                & - caqsupp_tmp &
                & ) 
            flx_aq(ispa,ires,iz) = sum(flx_aq(ispa,:,iz))
            if (isnan(flx_aq(ispa,ires,iz))) then 
                print *,chraq(ispa),iz,(flx_aq(ispa,iflx,iz),iflx=1,nflx)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/(poro(iz)*sat(iz)*1d3)
            ! ymx3(row) = ymx3(row)/(poro(iz)*sat(iz)*1d3)
        
        enddo 
        
    end do  ! ==============================
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    pCO2    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
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
    
        khgas(ispg,:) = kco2*(1d0+k1/pro + k1*k2/pro/pro) ! previous value; should not change through iterations 
        khgasx(ispg,:) = kco2*(1d0+k1/prox + k1*k2/prox/prox)
        
        dkhgas_dpro(ispg,:) = kco2*(k1*(-1d0)/prox**2d0 + k1*k2*(-2d0)/prox**3d0)
        
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
                & -drxngas_dmgas(ispg,ispg,iz)*dt &
                & ) &
                & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg))

            ymx3(row) = ( &
                & (agasx(ispg,iz)*mgasx(ispg,iz)-agas(ispg,iz)*mgas(ispg,iz)) &
                & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,min(nz,iz+1)))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & - 0.5d0*(dgas(ispg,iz)+edifn_tmp)*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                & +poro(iz)*sat(iz)*v(iz)*1d3*(khgasx(ispg,iz)*mgasx(ispg,iz)-khco2n_tmp*pco2n_tmp)/dz(iz)*dt &
                & -resp(iz)*dt &
                & -rxngas(ispg,iz)*dt &
                & -mgassupp(ispg,iz)*dt &
                & ) &
                & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg))
            
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                        & + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg))
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 
                    amx3(row,col+nsp3) = ( &
                            & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                            & ) &
                            & *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg))
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
                    & *merge(0.0d0,mgasx(ispg,iz-1),pco2x(iz)<pco2th)
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 

                    amx3(row,col-nsp3) = ( &
                        & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                        &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg))
                enddo 
            endif 
            
            do isps = 1,nsp_sld
                col = nsp3*(iz-1) + isps 
                amx3(row,col) = ( &
                    & -drxngas_dmsld(ispg,isps,iz)*dt &
                    & ) &
                    & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg))
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
                    & ) &
                    & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg))
                    
                ! if (trim(adjustl(chraq(ispa)))=='ca') then 
                    ! print *, &
                    ! & (dagas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz)) &
                    ! & ,-( 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,min(nz,iz+1))-mgasx(ispg,iz)) &
                    ! &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                    ! & - 0.5d0*(ddgas_dmaq(ispg,ispa,iz))*(mgasx(ispg,iz)-pco2n_tmp)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))) )/dz(iz)*dt  &
                    ! & ,+poro(iz)*sat(iz)*v(iz)*1d3*(dkhgas_dmaq(ispg,ispa,iz)*mgasx(ispg,iz))/dz(iz)*dt &
                    ! & ,-drxngas_dmaq(ispg,ispa,iz)*dt 
                ! endif 
                
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
            flx_gas(ispg,irxn_fo,iz) = -resp(iz)
            flx_gas(ispg,irain,iz) = -rxngas(ispg,iz) - mgassupp(ispg,iz)
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
            
            if (any(isnan(flx_gas(ispg,:,iz)))) then
                print *,flx_gas(ispg,:,iz)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/alpha(iz)
            ! ymx3(row) = ymx3(row)/alpha(iz)
        enddo 

    end do 
    
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
        do isps = 1, nsp_sld
            row = isps + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chrsld(isps)
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 

    co2 = kco2*mgasx(findloc(chrgas,'pco2',dim=1),:)
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

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
    print *
    print *,'-=-=-=-=-=-= Aq species -=-=-=-=-=-=-='
    do ispa = 1, nsp_aq
        print *, trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Sld species -=-=-=-=-=-=-='
    do isps = 1, nsp_sld
        print *, trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/5)
    enddo 
    do isps = 1, nsp_sld
        print *, 'omega_'//trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Gas species -=-=-=-=-=-=-='
    do ispg = 1, nsp_gas
        print *, trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
#endif     
enddo

mfox(:)=msldx(findloc(chrsld,'fo',dim=1),:)
mabx(:)=msldx(findloc(chrsld,'ab',dim=1),:)
manx(:)=msldx(findloc(chrsld,'an',dim=1),:)
mccx(:)=msldx(findloc(chrsld,'cc',dim=1),:)
mkax(:)=msldx(findloc(chrsld,'ka',dim=1),:)

mgx(:)=maqx(findloc(chraq,'mg',dim=1),:)
six(:)=maqx(findloc(chraq,'si',dim=1),:)
nax(:)=maqx(findloc(chraq,'na',dim=1),:)
cax(:)=maqx(findloc(chraq,'ca',dim=1),:)
alx(:)=maqx(findloc(chraq,'al',dim=1),:)

pco2x(:)=mgasx(findloc(chrgas,'pco2',dim=1),:)

omega_fo =omega(findloc(chrsld,'fo',dim=1),:)
omega_ab =omega(findloc(chrsld,'ab',dim=1),:)
omega_an =omega(findloc(chrsld,'an',dim=1),:)
omega_ka =omega(findloc(chrsld,'ka',dim=1),:)
omega_cc =omega(findloc(chrsld,'cc',dim=1),:)

flx_fo =flx_sld(findloc(chrsld,'fo',dim=1),:,:)
flx_ab =flx_sld(findloc(chrsld,'ab',dim=1),:,:)
flx_an =flx_sld(findloc(chrsld,'an',dim=1),:,:)
flx_ka =flx_sld(findloc(chrsld,'ka',dim=1),:,:)
flx_cc =flx_sld(findloc(chrsld,'cc',dim=1),:,:)

flx_mg=flx_aq(findloc(chraq,'mg',dim=1),:,:)
flx_si=flx_aq(findloc(chraq,'si',dim=1),:,:)
flx_na=flx_aq(findloc(chraq,'na',dim=1),:,:)
flx_ca=flx_aq(findloc(chraq,'ca',dim=1),:,:)
flx_al=flx_aq(findloc(chraq,'al',dim=1),:,:)

flx_co2=flx_gas(findloc(chrgas,'pco2',dim=1),:,:)


endsubroutine alsilicate_dis_co2_1D_v3














!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine alsilicate_aq_gas_1D( &
    & nz,mfo,mab,man,mcc,na,mg,si,ca,hr,poro,z,dz,w,kfo,kab,kan,kcc,keqfo,keqab,keqan,keqcc,mfoth,mabth,manth,mccth &! input
    & ,dmg,dsi,dna,dca,sat,dporodta,pro,mfoi,mabi,mani,mcci,mfosupp,mabsupp,mansupp,mccsupp,poroprev  &! input
    & ,kco2,k1,k2,mgth,sith,nath,cath,tora,v,tol,zrxn,it,cx,c2x,so4x,mgi,sii,nai,cai,mvfo,mvab,mvan,mvcc,nflx,kw &! input
    & ,k1si,k2si,kcca,keqcca,authig,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,pco2,pco2i,khco2i,ucv,torg,dgasc,daqc &! input 
    & ,pco2th,k1al,k2al,k3al,k4al,keqka,kka,al,mvka,ali,mkai,alth,mkath,dal,mkasupp,mka &! intput
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,cplprec,po2,po2th,po2i,kho,po2supp,dgaso,daqo,vmax,mo2 &! input
    & ,iter,error,dt,flgback &! inout
    & ,mgx,six,nax,cax,prox,co2,hco3,co3,dic,mfox,mabx,manx,mccx,omega_fo,omega_ab,omega_an,omega_cc,flx_fo,flx_mg &! output
    & ,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,omega_cca,pco2x,flx_co2,alx,flx_ka,flx_al,omega_ka,mkax,po2x,flx_o2,resp &! output
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,mfoth,tol,zrxn,dmg,dsi,mgth,sith,kco2,k1,k2,mfoi,mgi,sii,keqfo,mvfo,keqab,mabth,dna,mabi,nath &
    & ,nai,mvab,kw,keqan,manth,dca,mani,cath,cai,mvan,keqcc,mccth,mcci,mvcc,k1si,k2si,keqcca,authig,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3,pco2i,khco2i,ucv,dgasc,daqc,pco2th,k1al,k2al,k3al,k4al,keqka,mvka,ali,mkai,alth,mkath,dal &
    & ,po2th,po2i,kho,dgaso,daqo,vmax,mo2
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,dporodta,tora,v,mfo,kfo,cx,c2x,so4x,ca,pro,mfosupp,mg,si,mab,na,kab,mabsupp &
    & ,man,kan,mansupp,mcc,kcc,mccsupp,poroprev,dz,kcca,pco2,torg,kka,al,mkasupp,mka &
    & ,alsupp,sisupp,casupp,pco2supp,mgsupp,nasupp,po2,po2supp
real(kind=8),dimension(nz),intent(inout)::mgx,six,mfox,nax,mabx,cax,manx,mccx,pco2x,alx,mkax,po2x
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic,omega_fo,omega_ab,omega_an,omega_cc,omega_cca,omega_ka,resp
real(kind=8),dimension(nflx,nz),intent(out)::flx_fo,flx_mg,flx_si,flx_ab,flx_na,flx_an,flx_ca,flx_cc,flx_co2,flx_ka,flx_al &
    & ,flx_o2
integer,intent(inout)::iter,it
logical,intent(in)::cplprec
logical,intent(inout)::flgback
real(kind=8),intent(inout)::error,dt

integer,parameter::nsp_sld = 5
integer,parameter::nsp_sld_2 = 2
integer,parameter::nsp_aq = 5
integer,parameter::nsp_aq_cnst = 4
integer,parameter::nsp_aq_ph = 5
integer,parameter::nsp_gas_ph = 1
integer,parameter::nsp_aq_all = 9
integer,parameter::nsp_gas = 2
integer,parameter::nsp3 = nsp_sld + nsp_aq + nsp_gas
integer,parameter::nrxn_ext = 1
character(5),dimension(nsp_sld)::chrsld
character(5),dimension(nsp_sld_2)::chrsld_2
character(5),dimension(nsp_aq)::chraq
character(5),dimension(nsp_aq_cnst)::chraq_cnst
character(5),dimension(nsp_aq_ph)::chraq_ph
character(5),dimension(nsp_gas_ph)::chrgas_ph
character(5),dimension(nsp_aq_all)::chraq_all
character(5),dimension(nsp_gas)::chrgas
character(5),dimension(nrxn_ext)::chrrxn_ext
real(kind=8),dimension(nsp_sld)::msldi,msldth,keq,mv
real(kind=8),dimension(nsp_aq)::maqi,maqth,daq,dmaq 
real(kind=8),dimension(nsp_gas)::mgasi,mgasth,dgasa,dgasg,dmgas,khgasi,dgasi
real(kind=8),dimension(nsp_sld,nsp_aq)::staq
real(kind=8),dimension(nsp_sld,nsp_gas)::stgas
real(kind=8),dimension(nsp_sld,nz)::msldx,msld,ksld,omega,nonprec,rxnsld,msldsupp,domega_dpro 
real(kind=8),dimension(nsp_sld,nsp3,nz)::domega_disp,drxnsld_disp 
real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq,drxnsld_dmaq 
real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas,drxnsld_dmgas 
real(kind=8),dimension(nsp_sld,nflx,nz)::flx_sld
real(kind=8),dimension(nsp_aq_cnst,nz)::maqcnst
real(kind=8),dimension(nsp_aq,nz)::maqx,maq,rxnaq,maqsupp,drxnaq_dpro,dprodmaq 
real(kind=8),dimension(nsp_aq,nsp3,nz)::drxnaq_disp
real(kind=8),dimension(nsp_aq,nflx,nz)::flx_aq
real(kind=8),dimension(nsp_gas,nz)::mgasx,mgas,khgasx,khgas,dgas,agasx,agas,rxngas,mgassupp,dkhgas_dpro,ddgas_dpro,dagas_dpro &
    & ,dprodmgas
real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
real(kind=8),dimension(nsp_gas,nflx,nz)::flx_gas 
real(kind=8),dimension(nrxn_ext,nz)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_gas)::stgas_ext
real(kind=8),dimension(nrxn_ext,nsp_gas)::stgas_dext
real(kind=8),dimension(nrxn_ext,nsp_aq)::staq_ext
real(kind=8),dimension(nrxn_ext,nsp_sld)::stsld_ext
real(kind=8),dimension(nrxn_ext,nsp_gas,nz)::drxnext_dmgas
real(kind=8),dimension(nrxn_ext,nsp_aq,nz)::drxnext_dmaq
real(kind=8),dimension(nrxn_ext,nsp_sld,nz)::drxnext_dmsld

integer iz,row,nmx,ie,ie2,isp,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8),dimension(nz)::dprodna,dprodmg,domega_fo_dmg,domega_fo_dsi,domega_ab_dsi,domega_ab_dna,domega_ab_dmg,domega_fo_dna &
    & ,domega_ab_dca,domega_fo_dca,domega_an_dsi,domega_an_dna,domega_an_dmg,domega_an_dca,dprodca,domega_cc_dca,domega_cc_dna &
    & ,domega_cc_dsi,domega_cc_dmg,dprodsi,domega_fo_dpro,domega_ab_dpro,domega_an_dpro,domega_cc_dpro &
    & ,domega_cca_dmg,domega_cca_dca,domega_cca_dna,domega_cca_dsi,domega_cca_dpro,dprodpco2,domega_fo_dpco2,domega_ab_dpco2 &
    & ,domega_an_dpco2,domega_cc_dpco2,domega_cca_dpco2,dkhco2_dpro,dkhco2_dpco2,dpreccc_dpco2,preccc,khco2,khco2x,alpha,alphaprev &
    & ,dalpha_dpco2,edif,dedif_dpco2,dalpha_dna,dalpha_dsi,dalpha_dmg,dalpha_dca,dedif_dca,dedif_dmg,dedif_dsi,dedif_dna &
    & ,dkhco2_dca,dkhco2_dna,dkhco2_dsi,dkhco2_dmg,dpreccc_dca,dpreccc_dna,dpreccc_dsi,dpreccc_dmg,dpreccc_dmcc,dprodall &
    & ,domega_fo_dal,domega_ka_dal,domega_ka_dsi,domega_ka_dpco2,domega_ka_dpro,domega_ka_dna,domega_ka_dca,domega_ka_dmg  &
    & ,domega_ab_dal,domega_an_dal,domega_cc_dal,domega_cca_dal,dalpha_dal,dkhco2_dal,dedif_dal,dpreccc_dal,dprodal
real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp,st_fo,st_ab &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,st_an,omega_tmp_th,st_cc &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,st_cca,edifi,khco2n_tmp,pco2n_tmp,edifn_tmp,st_ka,caqsupp_tmp

real(kind=8),parameter::sec2yr = 60d0*60d0*60d0*24d0*365d0
real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 100d0
real(kind=8),dimension(nz)::disonly ! for cc [1---yes, 0---no]
real(kind=8),dimension(nz)::disonly_ka ! for ka [1---yes, 0---no]
real(kind=8),dimension(nz)::zeros  
! real(kind=8)::disonly = 1d0 ! for cc 
! real(kind=8)::authig = 0d0 ! for cc whether to include authigenesis [1---yes, 0---no]
! real(kind=8)::authig = 1d0 ! 

integer,parameter :: iter_max = 100

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV

! passing major variables to more compact matrices 

chrsld = (/'fo','ab','an','cc','ka'/)
chraq = (/'mg','si','na','ca','al'/)
chrgas = (/'pco2','po2 '/)

! previous values
msld(findloc(chrsld,'fo',dim=1),:)=mfo(:)
msld(findloc(chrsld,'ab',dim=1),:)=mab(:)
msld(findloc(chrsld,'an',dim=1),:)=man(:)
msld(findloc(chrsld,'cc',dim=1),:)=mcc(:)
msld(findloc(chrsld,'ka',dim=1),:)=mka(:)

maq(findloc(chraq,'mg',dim=1),:)=mg(:)
maq(findloc(chraq,'si',dim=1),:)=si(:)
maq(findloc(chraq,'na',dim=1),:)=na(:)
maq(findloc(chraq,'ca',dim=1),:)=ca(:)
maq(findloc(chraq,'al',dim=1),:)=al(:)

mgas(findloc(chrgas,'pco2',dim=1),:)=pco2(:)
mgas(findloc(chrgas,'po2',dim=1),:)=po2(:)

! output
msldx(findloc(chrsld,'fo',dim=1),:)=mfox(:)
msldx(findloc(chrsld,'ab',dim=1),:)=mabx(:)
msldx(findloc(chrsld,'an',dim=1),:)=manx(:)
msldx(findloc(chrsld,'cc',dim=1),:)=mccx(:)
msldx(findloc(chrsld,'ka',dim=1),:)=mkax(:)

maqx(findloc(chraq,'mg',dim=1),:)=mgx(:)
maqx(findloc(chraq,'si',dim=1),:)=six(:)
maqx(findloc(chraq,'na',dim=1),:)=nax(:)
maqx(findloc(chraq,'ca',dim=1),:)=cax(:)
maqx(findloc(chraq,'al',dim=1),:)=alx(:)

mgasx(findloc(chrgas,'pco2',dim=1),:)=pco2x(:)
mgasx(findloc(chrgas,'po2',dim=1),:)=po2x(:)

! constants
msldi = (/mfoi,mabi,mani,mcci,mkai/)
msldth = (/mfoth,mabth,manth,mccth,mkath/)
maqi = (/mgi,sii,nai,cai,ali/)
maqth = (/mgth,sith,nath,cath,alth/)
mgasth = (/pco2th,po2th/)
mgasi = (/pco2i,po2i/)

ksld(findloc(chrsld,'fo',dim=1),:) = kfo
ksld(findloc(chrsld,'ab',dim=1),:) = kab
ksld(findloc(chrsld,'an',dim=1),:) = kan
ksld(findloc(chrsld,'cc',dim=1),:) = kcc
ksld(findloc(chrsld,'ka',dim=1),:) = kka

msldsupp(findloc(chrsld,'fo',dim=1),:) = mfosupp
msldsupp(findloc(chrsld,'ab',dim=1),:) = mabsupp
msldsupp(findloc(chrsld,'an',dim=1),:) = mansupp
msldsupp(findloc(chrsld,'cc',dim=1),:) = mccsupp
msldsupp(findloc(chrsld,'ka',dim=1),:) = mkasupp

maqsupp(findloc(chraq,'mg',dim=1),:) = mgsupp
maqsupp(findloc(chraq,'si',dim=1),:) = sisupp
maqsupp(findloc(chraq,'na',dim=1),:) = nasupp
maqsupp(findloc(chraq,'ca',dim=1),:) = casupp
maqsupp(findloc(chraq,'al',dim=1),:) = alsupp

mgassupp(findloc(chrgas,'pco2',dim=1),:) = pco2supp
mgassupp(findloc(chrgas,'po2',dim=1),:) = po2supp

daq = (/dmg,dsi,dna,dca,dal/)

dgasg = (/dgasc,dgaso/)
dgasa = (/daqc,daqo/)

khgasi = (/khco2i,kho/)

mv = (/mvfo,mvab,mvan,mvcc,mvka/)

chraq_cnst = (/'so4  ','fe2  ','fe3  ','k    '/)
maqcnst = 0d0

chrsld_2 = (/'cc','ka'/)

chraq_ph = (/'mg','si','na','ca','al'/)
chrgas_ph = (/'pco2'/)

chrrxn_ext = (/'resp'/)

staq = 0d0
stgas = 0d0
! Forsterite; Mg2SiO4
staq(findloc(chrsld,'fo',dim=1), findloc(chraq,'mg',dim=1)) = 2d0
staq(findloc(chrsld,'fo',dim=1), findloc(chraq,'si',dim=1)) = 1d0
! Albite; NaAlSi3O8
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'na',dim=1)) = 1d0
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'si',dim=1)) = 3d0
staq(findloc(chrsld,'ab',dim=1), findloc(chraq,'al',dim=1)) = 1d0
! Anothite; CaAl2Si2O8
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'ca',dim=1)) = 1d0
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'si',dim=1)) = 2d0
staq(findloc(chrsld,'an',dim=1), findloc(chraq,'al',dim=1)) = 2d0
! Calcite; CaCO3
staq(findloc(chrsld,'cc',dim=1), findloc(chraq,'ca',dim=1)) = 1d0
stgas(findloc(chrsld,'cc',dim=1), findloc(chrgas,'pco2',dim=1)) = 1d0
! Kaolinite; Al2Si2O5(OH)4
staq(findloc(chrsld,'ka',dim=1), findloc(chraq,'si',dim=1)) = 2d0
staq(findloc(chrsld,'ka',dim=1), findloc(chraq,'al',dim=1)) = 2d0


staq_ext = 0d0
stgas_ext = 0d0
stsld_ext = 0d0
! respiration 
stgas_ext(findloc(chrrxn_ext,'resp',dim=1), findloc(chrgas,'pco2',dim=1)) = 1d0
stgas_ext(findloc(chrrxn_ext,'resp',dim=1), findloc(chrgas,'po2',dim=1)) = -1d0

stgas_dext = 0d0
stgas_dext(findloc(chrrxn_ext,'resp',dim=1), findloc(chrgas,'po2',dim=1)) = 1d0

! print *, staq
! print *, mgx
! print *, six
! print *, mfosupp
! stop

if (any(isnan(tora)))then 
    print*,tora
endif 

zeros = 0d0

error = 1d4
iter = 0

if (it ==0) then 

    ! call calc_pH_v3( &
        ! & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,zeros,zeros,zeros,zeros,mgasx(findloc(chrgas,'pco2',dim=1),:) &
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:) &
        & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
else 

    ! call calc_pH_v3( &
        ! & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:) &
        & ,maqx(findloc(chraq,'ca',dim=1),:),zeros,mgasx(findloc(chrgas,'pco2',dim=1),:) &
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:) &
        & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
endif

! print *, 'starting silciate calculation'

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_sld = 0d0
    flx_aq = 0d0
    flx_gas = 0d0
    
    ! pH calculation and its derivative wrt aq and gas species
    
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
    do ispa=1,nsp_aq
        dmaq = 0d0
        dmgas = 0d0
        if (any (chraq_ph == chraq(ispa))) then 
            dmaq(ispa) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmaq(ispa,:) &! output
                & ) 
            dprodmaq(ispa,:) = (dprodmaq(ispa,:) - prox(:))/dconc
        endif 
    enddo 
    
    dprodmgas = 0d0
    do ispg=1,nsp_gas
        dmaq = 0d0
        dmgas = 0d0
        if (any (chrgas_ph == chrgas(ispg))) then 
            dmgas(ispg) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmgas(ispg,:) &! output
                & ) 
            dprodmgas(ispg,:) = (dprodmgas(ispg,:) - prox(:))/dconc
        endif 
    enddo
    
    ! saturation state calc. and their derivatives wrt aq and gas species
    
    omega = 0d0
    domega_dpro = 0d0
    domega_dmaq = 0d0
    domega_dmgas = 0d0
    
    do isps =1, nsp_sld
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:) &
            & ,maqx(findloc(chraq,'ca',dim=1),:) &
            & ,maqx(findloc(chraq,'mg',dim=1),:) &
            & ,maqx(findloc(chraq,'si',dim=1),:) &
            & ,maqx(findloc(chraq,'na',dim=1),:) &
            & ,prox &
            & ,maqx(findloc(chraq,'al',dim=1),:) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,omega(isps,:) &! output
            & )
        dmaq = 0d0
        dmgas = 0d0
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
            & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
            & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
            & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
            & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
            & ,prox+dconc &
            & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,domega_dpro(isps,:) &! output
            & )
        domega_dpro(isps,:) = (domega_dpro(isps,:)-omega(isps,:))/dconc
        do ispa = 1, nsp_aq
            dmaq = 0d0
            dmgas = 0d0
            if (any (chraq_ph == chraq(ispa))) then 
                dmaq(ispa) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmaq(isps,ispa,:) &! output
                    & )
                domega_dmaq(isps,ispa,:) = (domega_dmaq(isps,ispa,:)-omega(isps,:))/dconc
                domega_dmaq(isps,ispa,:) = domega_dmaq(isps,ispa,:) + domega_dpro(isps,:)*dprodmaq(ispa,:)
            endif 
        enddo
        do ispg = 1, nsp_gas
            dmaq = 0d0
            dmgas = 0d0
            if (any (chrgas_ph == chrgas(ispg))) then 
                dmgas(ispg) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmgas(isps,ispg,:) &! output
                    & )
                domega_dmgas(isps,ispg,:) = (domega_dmgas(isps,ispg,:)-omega(isps,:))/dconc
                domega_dmgas(isps,ispg,:) = domega_dmgas(isps,ispg,:) + domega_dpro(isps,:)*dprodmgas(ispg,:)
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
        call calc_rxn_ext( &
            & nz,vmax,mo2,po2th,po2x,chrrxn_ext(irxn) &! input 
            & ,rxnext(irxn,:) &! output
            & )
        dmgas = 0d0
        do ispg=1,nsp_gas
            if (stgas_dext(irxn,ispg)==0d0) cycle
            dmgas(ispg) = dconc
            call calc_rxn_ext( &
                & nz,vmax,mo2,po2th,mgasx(findloc(chrgas,'po2',dim=1),:)+dmgas(findloc(chrgas,'po2',dim=1)) &
                & ,chrrxn_ext(irxn) &! input 
                & ,drxnext_dmgas(irxn,ispg,:)  &! output
                & )
            drxnext_dmgas(irxn,ispg,:) = (drxnext_dmgas(irxn,ispg,:) -rxnext(irxn,:))/dconc
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
            mprev_tmp = msld(isps,iz)
            
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
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1)+nsp_sld + ispa
                
                amx3(row,col ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmaq(isps,ispa,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *maqx(ispa,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            do ispg = 1, nsp_gas 
                col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                amx3(row,col) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmgas(isps,ispg,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            flx_sld(isps,itflx,iz) = ( &
                & (m_tmp-mprev_tmp)/dt &
                & )
            flx_sld(isps,iadv,iz) = ( &
                & -w*(mp_tmp-m_tmp)/dz(iz)  &
                & )
            flx_sld(isps,irxn_fo,iz) = ( &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & )
            flx_sld(isps,irain,iz) = (&
                & - msupp_tmp  &
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
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*dt &
                & - rxn_tmp*dt &
                & - caqsupp_tmp*dt &
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
            
            do isps = 1, nsp_sld
                col = nsp3*(iz-1)+ isps
                
                amx3(row, col) = (     & 
                    & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                    & ) &
                    & *msldx(isps,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            enddo 
            
            do ispa2 = 1, nsp_aq
                col = nsp3*(iz-1)+ nsp_sld + ispa2
                
                if (ispa2 == ispa) cycle
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmaq(isps,ispa2,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *maqx(ispa2,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                enddo 
            enddo 
            
            do ispg = 1, nsp_gas
                col = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmgas(isps,ispg,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *mgasx(ispg,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
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
            flx_aq(ispa,irxn_fo,iz) = (&
                & - rxn_tmp &
                & ) 
            flx_aq(ispa,irain,iz) = (&
                & - caqsupp_tmp &
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
                & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg))
            
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
                & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg))
            
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                        & + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg))
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 
                    amx3(row,col+nsp3) = ( &
                            & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                            & ) &
                            & *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg))
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
                    & *merge(0.0d0,mgasx(ispg,iz-1),mgasx(ispg,iz)<mgasth(ispg))
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 

                    amx3(row,col-nsp3) = ( &
                        & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                        &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg))
                enddo 
            endif 
            
            do isps = 1,nsp_sld
                col = nsp3*(iz-1) + isps 
                amx3(row,col) = ( &
                    & -drxngas_dmsld(ispg,isps,iz)*dt &
                    & ) &
                    & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg))
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
                    & ) &
                    & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg))
                
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
                    & *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg))
                
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
            ! flx_gas(ispg,irxn_fo,iz) = -resp(iz)
            flx_gas(ispg,irxn_fo,iz) = -sum(stgas_ext(:,ispg)*rxnext(:,iz))
            flx_gas(ispg,irain,iz) = -rxngas(ispg,iz) - mgassupp(ispg,iz)
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
            
            if (any(isnan(flx_gas(ispg,:,iz)))) then
                print *,flx_gas(ispg,:,iz)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/alpha(iz)
            ! ymx3(row) = ymx3(row)/alpha(iz)
        enddo 

    end do 
    
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
        do isps = 1, nsp_sld
            row = isps + nsp3*(iz-1)

            if (isnan(ymx3(row))) then 
                print *,'nan at', iz,z(iz),chrsld(isps)
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
                if (z(iz)<zrxn) then 
                    ymx3(row)=0d0
                endif
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
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 

    co2 = kco2*mgasx(findloc(chrgas,'pco2',dim=1),:)
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

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
    print *
    print *,'-=-=-=-=-=-= Aq species -=-=-=-=-=-=-='
    do ispa = 1, nsp_aq
        print *, trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Sld species -=-=-=-=-=-=-='
    do isps = 1, nsp_sld
        print *, trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/5)
    enddo 
    do isps = 1, nsp_sld
        print *, 'omega_'//trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Gas species -=-=-=-=-=-=-='
    do ispg = 1, nsp_gas
        print *, trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
#endif     
enddo

mfox(:)=msldx(findloc(chrsld,'fo',dim=1),:)
mabx(:)=msldx(findloc(chrsld,'ab',dim=1),:)
manx(:)=msldx(findloc(chrsld,'an',dim=1),:)
mccx(:)=msldx(findloc(chrsld,'cc',dim=1),:)
mkax(:)=msldx(findloc(chrsld,'ka',dim=1),:)

mgx(:)=maqx(findloc(chraq,'mg',dim=1),:)
six(:)=maqx(findloc(chraq,'si',dim=1),:)
nax(:)=maqx(findloc(chraq,'na',dim=1),:)
cax(:)=maqx(findloc(chraq,'ca',dim=1),:)
alx(:)=maqx(findloc(chraq,'al',dim=1),:)

pco2x(:)=mgasx(findloc(chrgas,'pco2',dim=1),:)
po2x(:)=mgasx(findloc(chrgas,'po2',dim=1),:)

omega_fo =omega(findloc(chrsld,'fo',dim=1),:)
omega_ab =omega(findloc(chrsld,'ab',dim=1),:)
omega_an =omega(findloc(chrsld,'an',dim=1),:)
omega_ka =omega(findloc(chrsld,'ka',dim=1),:)
omega_cc =omega(findloc(chrsld,'cc',dim=1),:)

flx_fo =flx_sld(findloc(chrsld,'fo',dim=1),:,:)
flx_ab =flx_sld(findloc(chrsld,'ab',dim=1),:,:)
flx_an =flx_sld(findloc(chrsld,'an',dim=1),:,:)
flx_ka =flx_sld(findloc(chrsld,'ka',dim=1),:,:)
flx_cc =flx_sld(findloc(chrsld,'cc',dim=1),:,:)

flx_mg=flx_aq(findloc(chraq,'mg',dim=1),:,:)
flx_si=flx_aq(findloc(chraq,'si',dim=1),:,:)
flx_na=flx_aq(findloc(chraq,'na',dim=1),:,:)
flx_ca=flx_aq(findloc(chraq,'ca',dim=1),:,:)
flx_al=flx_aq(findloc(chraq,'al',dim=1),:,:)

flx_co2=flx_gas(findloc(chrgas,'pco2',dim=1),:,:)
flx_o2=flx_gas(findloc(chrgas,'po2',dim=1),:,:)

resp = rxnext(findloc(chrrxn_ext,'resp',dim=1),:)

endsubroutine alsilicate_aq_gas_1D















!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

subroutine alsilicate_aq_gas_1D_v2( &
    ! new input 
    & nz,nsp_sld,nsp_sld_2,nsp_aq,nsp_aq_ph,nsp_gas_ph,nsp_gas,nsp3,nrxn_ext &
    & ,chrsld,chrsld_2,chraq,chraq_ph,chrgas_ph,chrgas,chrrxn_ext  &
    & ,msldi,msldth,mv,maqi,maqth,daq,mgasi,mgasth,dgasa,dgasg,khgasi &
    & ,staq,stgas,msld,ksld,msldsupp,maq,maqsupp,mgas,mgassupp &
    & ,stgas_ext,stgas_dext,staq_ext,stsld_ext &
    !  old inputs
    & ,hr,poro,z,dz,w,keqfo,keqab,keqan,keqcc,sat,pro,poroprev &
    & ,kco2,k1,k2,tora,v,tol,it,so4x,nflx,kw,k1si,k2si & 
    & ,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3,ucv,torg &
    & ,k1al,k2al,k3al,k4al,keqka,cplprec,kho,vmax,mo2  &
    ! old inout
    & ,iter,error,dt,flgback &    
    ! output 
    & ,msldx,omega,flx_sld,maqx,flx_aq,mgasx,flx_gas,rxnext,prox,co2,hco3,co3,dic & 
    & )
    
implicit none 

integer,intent(in)::nz,nflx
real(kind=8),intent(in)::w,tol,kco2,k1,k2,keqfo,keqab,kw,keqan,keqcc,k1si,k2si,k1mg,k1mgco3,k1mghco3 &
    & ,k1ca,k1caco3,k1cahco3,ucv,k1al,k2al,k3al,k4al,keqka,kho,vmax,mo2
real(kind=8),dimension(nz),intent(in)::hr,poro,z,sat,tora,v,so4x,poroprev,dz,torg,pro
real(kind=8),dimension(nz),intent(out)::prox,co2,hco3,co3,dic
integer,intent(inout)::iter,it
logical,intent(in)::cplprec
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
real(kind=8),dimension(nsp_aq)::dmaq 
real(kind=8),dimension(nsp_gas),intent(in)::mgasi,mgasth,dgasa,dgasg,khgasi
real(kind=8),dimension(nsp_gas)::dmgas,dgasi
real(kind=8),dimension(nsp_sld,nsp_aq),intent(in)::staq
real(kind=8),dimension(nsp_sld,nsp_gas),intent(in)::stgas
real(kind=8),dimension(nsp_sld,nz),intent(in)::msld,ksld,msldsupp 
real(kind=8),dimension(nsp_sld,nz),intent(inout)::msldx,omega
real(kind=8),dimension(nsp_sld,nz)::nonprec,domega_dpro 
real(kind=8),dimension(nsp_sld,nsp_aq,nz)::domega_dmaq
real(kind=8),dimension(nsp_sld,nsp_gas,nz)::domega_dmgas
real(kind=8),dimension(nsp_sld,nflx,nz),intent(out)::flx_sld
real(kind=8),dimension(nsp_aq,nz),intent(in)::maq,maqsupp
real(kind=8),dimension(nsp_aq,nz),intent(inout)::maqx 
real(kind=8),dimension(nsp_aq,nz)::dprodmaq 
real(kind=8),dimension(nsp_aq,nflx,nz),intent(out)::flx_aq
real(kind=8),dimension(nsp_gas,nz),intent(in)::mgas,mgassupp
real(kind=8),dimension(nsp_gas,nz),intent(inout)::mgasx 
real(kind=8),dimension(nsp_gas,nz)::khgasx,khgas,dgas,agasx,agas,rxngas,dkhgas_dpro,dprodmgas
real(kind=8),dimension(nsp_gas,nsp_aq,nz)::dkhgas_dmaq,ddgas_dmaq,dagas_dmaq,drxngas_dmaq 
real(kind=8),dimension(nsp_gas,nsp_sld,nz)::drxngas_dmsld 
real(kind=8),dimension(nsp_gas,nsp_gas,nz)::dkhgas_dmgas,ddgas_dmgas,dagas_dmgas,drxngas_dmgas 
real(kind=8),dimension(nsp_gas,nflx,nz),intent(out)::flx_gas 
real(kind=8),dimension(nrxn_ext,nz),intent(inout)::rxnext
real(kind=8),dimension(nrxn_ext,nsp_gas),intent(in)::stgas_ext,stgas_dext
real(kind=8),dimension(nrxn_ext,nsp_aq),intent(in)::staq_ext
real(kind=8),dimension(nrxn_ext,nsp_sld),intent(in)::stsld_ext
real(kind=8),dimension(nrxn_ext,nsp_gas,nz)::drxnext_dmgas
real(kind=8),dimension(nrxn_ext,nsp_aq,nz)::drxnext_dmaq
real(kind=8),dimension(nrxn_ext,nsp_sld,nz)::drxnext_dmsld

integer iz,row,ie,ie2,iflx,isps,ispa,ispg,ispa2,ispg2,col,irxn
integer::itflx,iadv,idif,irxn_fo,irain,ires
data itflx,iadv,idif,irxn_fo,irain,ires/1,2,3,4,5,6/

real(kind=8) d_tmp,caq_tmp,caq_tmp_p,caq_tmp_n,caqth_tmp,caqi_tmp,rxn_tmp,caq_tmp_prev,drxndisp_tmp &
    & ,k_tmp,mv_tmp,omega_tmp,m_tmp,mth_tmp,mi_tmp,mp_tmp,msupp_tmp,mprev_tmp,omega_tmp_th &
    & ,edif_tmp,edif_tmp_n,edif_tmp_p,khco2n_tmp,pco2n_tmp,edifn_tmp,caqsupp_tmp

real(kind=8),parameter::infinity = huge(0d0)
real(kind=8)::dconc = 1d-14
real(kind=8)::threshold = 10d0
! real(kind=8)::threshold = 100d0
real(kind=8),dimension(nz)::zeros  

integer,parameter :: iter_max = 100

real(kind=8) amx3(nsp3*nz,nsp3*nz),ymx3(nsp3*nz)
integer ipiv3(nsp3*nz)
integer info

external DGESV



! print *, staq
! print *, mgx
! print *, six
! print *, mfosupp
! stop

if (any(isnan(tora)))then 
    print*,tora
endif 

zeros = 0d0

error = 1d4
iter = 0

if (it ==0) then 

    ! call calc_pH_v3( &
        ! & nz,cx,cx,cx,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,zeros,zeros,zeros,zeros,mgasx(findloc(chrgas,'pco2',dim=1),:) &
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:) &
        & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
else 

    ! call calc_pH_v3( &
        ! & nz,nax,mgx,cax,so4x,pco2x,kw,kco2,k1,k2,six,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        ! & ,prox &! output
        ! & ) 
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:) &
        & ,maqx(findloc(chraq,'ca',dim=1),:),zeros,mgasx(findloc(chrgas,'pco2',dim=1),:) &
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:) &
        & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
endif

! print *, 'starting silciate calculation'

do while ((.not.isnan(error)).and.(error > tol))

    amx3=0.0d0
    ymx3=0.0d0 
    
    flx_sld = 0d0
    flx_aq = 0d0
    flx_gas = 0d0
    
    ! pH calculation and its derivative wrt aq and gas species
    
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 
    do ispa=1,nsp_aq
        dmaq = 0d0
        dmgas = 0d0
        if (any (chraq_ph == chraq(ispa))) then 
            dmaq(ispa) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmaq(ispa,:) &! output
                & ) 
            dprodmaq(ispa,:) = (dprodmaq(ispa,:) - prox(:))/dconc
        endif 
    enddo 
    
    dprodmgas = 0d0
    do ispg=1,nsp_gas
        dmaq = 0d0
        dmgas = 0d0
        if (any (chrgas_ph == chrgas(ispg))) then 
            dmgas(ispg) = dconc
            call calc_pH_v4( &
                & nz,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) & 
                & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &! input
                & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &! input
                & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                & ,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
                & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)),k1al,k2al,k3al,k4al &! input
                & ,dprodmgas(ispg,:) &! output
                & ) 
            dprodmgas(ispg,:) = (dprodmgas(ispg,:) - prox(:))/dconc
        endif 
    enddo
    
    ! saturation state calc. and their derivatives wrt aq and gas species
    
    omega = 0d0
    domega_dpro = 0d0
    domega_dmaq = 0d0
    domega_dmgas = 0d0
    
    do isps =1, nsp_sld
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:) &
            & ,maqx(findloc(chraq,'ca',dim=1),:) &
            & ,maqx(findloc(chraq,'mg',dim=1),:) &
            & ,maqx(findloc(chraq,'si',dim=1),:) &
            & ,maqx(findloc(chraq,'na',dim=1),:) &
            & ,prox &
            & ,maqx(findloc(chraq,'al',dim=1),:) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,omega(isps,:) &! output
            & )
        dmaq = 0d0
        dmgas = 0d0
        call calc_omega_v2( &
            & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
            & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
            & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
            & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
            & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
            & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
            & ,prox+dconc &
            & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
            & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
            & ,domega_dpro(isps,:) &! output
            & )
        domega_dpro(isps,:) = (domega_dpro(isps,:)-omega(isps,:))/dconc
        do ispa = 1, nsp_aq
            dmaq = 0d0
            dmgas = 0d0
            if (any (chraq_ph == chraq(ispa))) then 
                dmaq(ispa) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmaq(isps,ispa,:) &! output
                    & )
                domega_dmaq(isps,ispa,:) = (domega_dmaq(isps,ispa,:)-omega(isps,:))/dconc
                domega_dmaq(isps,ispa,:) = domega_dmaq(isps,ispa,:) + domega_dpro(isps,:)*dprodmaq(ispa,:)
            endif 
        enddo
        do ispg = 1, nsp_gas
            dmaq = 0d0
            dmgas = 0d0
            if (any (chrgas_ph == chrgas(ispg))) then 
                dmgas(ispg) = dconc
                call calc_omega_v2( &
                    & nz,keqfo,keqab,keqan,keqcc,k1,k2,kco2,k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! inpuy
                    & ,mgasx(findloc(chrgas,'pco2',dim=1),:)+dmgas(findloc(chrgas,'pco2',dim=1)) &
                    & ,maqx(findloc(chraq,'ca',dim=1),:)+dmaq(findloc(chraq,'ca',dim=1)) &
                    & ,maqx(findloc(chraq,'mg',dim=1),:)+dmaq(findloc(chraq,'mg',dim=1)) &
                    & ,maqx(findloc(chraq,'si',dim=1),:)+dmaq(findloc(chraq,'si',dim=1)) &
                    & ,maqx(findloc(chraq,'na',dim=1),:)+dmaq(findloc(chraq,'na',dim=1)) &
                    & ,prox &
                    & ,maqx(findloc(chraq,'al',dim=1),:)+dmaq(findloc(chraq,'al',dim=1)) &
                    & ,k1al,k2al,k3al,k4al,keqka,chrsld(isps) &! input 
                    & ,domega_dmgas(isps,ispg,:) &! output
                    & )
                domega_dmgas(isps,ispg,:) = (domega_dmgas(isps,ispg,:)-omega(isps,:))/dconc
                domega_dmgas(isps,ispg,:) = domega_dmgas(isps,ispg,:) + domega_dpro(isps,:)*dprodmgas(ispg,:)
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
        call calc_rxn_ext( &
            & nz,vmax,mo2,mgasth(findloc(chrgas,'po2',dim=1)),mgasx(findloc(chrgas,'po2',dim=1),:) &
            & ,chrrxn_ext(irxn) &! input 
            & ,rxnext(irxn,:) &! output
            & )
        dmgas = 0d0
        do ispg=1,nsp_gas
            if (stgas_dext(irxn,ispg)==0d0) cycle
            dmgas(ispg) = dconc
            call calc_rxn_ext( &
                & nz,vmax,mo2,mgasth(findloc(chrgas,'po2',dim=1)) &
                & ,mgasx(findloc(chrgas,'po2',dim=1),:)+dmgas(findloc(chrgas,'po2',dim=1)) &
                & ,chrrxn_ext(irxn) &! input 
                & ,drxnext_dmgas(irxn,ispg,:)  &! output
                & )
            drxnext_dmgas(irxn,ispg,:) = (drxnext_dmgas(irxn,ispg,:) -rxnext(irxn,:))/dconc
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
            mprev_tmp = msld(isps,iz)
            
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
            
            do ispa = 1, nsp_aq
                col = nsp3*(iz-1)+nsp_sld + ispa
                
                amx3(row,col ) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmaq(isps,ispa,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *maqx(ispa,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            do ispg = 1, nsp_gas 
                col = nsp3*(iz-1)+nsp_sld + nsp_aq + ispg

                amx3(row,col) = ( &
                    & k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(-domega_dmgas(isps,ispg,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                    & ) &
                    & *mgasx(ispg,iz) &
                    & *merge(0.0d0,1d0,m_tmp<mth_tmp)
            enddo 
            
            flx_sld(isps,itflx,iz) = ( &
                & (m_tmp-mprev_tmp)/dt &
                & )
            flx_sld(isps,iadv,iz) = ( &
                & -w*(mp_tmp-m_tmp)/dz(iz)  &
                & )
            flx_sld(isps,irxn_fo,iz) = ( &
                & + k_tmp*poro(iz)*hr(iz)*mv_tmp*1d-6*m_tmp*(1d0-omega_tmp) &
                & *merge(0d0,1d0,1d0-omega_tmp_th < 0d0) &
                & )
            flx_sld(isps,irain,iz) = (&
                & - msupp_tmp  &
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
                & ) &
                & *merge(1.0d0,caq_tmp,caq_tmp<caqth_tmp)

            ymx3(row) = ( &
                & (poro(iz)*sat(iz)*1d3*caq_tmp-poroprev(iz)*sat(iz)*1d3*caq_tmp_prev)  &
                & -(0.5d0*(edif_tmp +edif_tmp_p)*(caq_tmp_p-caq_tmp)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                & -0.5d0*(edif_tmp +edif_tmp_n)*(caq_tmp-caq_tmp_n)/(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt &
                & + poro(iz)*sat(iz)*1d3*v(iz)*(caq_tmp-caq_tmp_n)/dz(iz)*dt &
                & - rxn_tmp*dt &
                & - caqsupp_tmp*dt &
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
            
            do isps = 1, nsp_sld
                col = nsp3*(iz-1)+ isps
                
                amx3(row, col) = (     & 
                    & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*1d0*(1d0-omega(isps,iz))*dt &
                    & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                    & ) &
                    & *msldx(isps,iz) &
                    & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
            enddo 
            
            do ispa2 = 1, nsp_aq
                col = nsp3*(iz-1)+ nsp_sld + ispa2
                
                if (ispa2 == ispa) cycle
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmaq(isps,ispa2,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *maqx(ispa2,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
                enddo 
            enddo 
            
            do ispg = 1, nsp_gas
                col = nsp3*(iz-1) + nsp_sld + nsp_aq + ispg
                
                do isps = 1, nsp_sld
                    amx3(row,col) = amx3(row,col) + (     & 
                        & - staq(isps,ispa)*ksld(isps,iz)*poro(iz)*hr(iz)*mv(isps)*1d-6*msldx(isps,iz) &
                        & *(-domega_dmgas(isps,ispg,iz))*dt &
                        & *merge(0d0,1d0,1d0-omega(isps,iz)*nonprec(isps,iz) < 0d0)  &
                        & ) &
                        & *mgasx(ispg,iz) &
                        & *merge(0.0d0,1.0d0,caq_tmp<caqth_tmp)   ! commented out (is this necessary?)
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
            flx_aq(ispa,irxn_fo,iz) = (&
                & - rxn_tmp &
                & ) 
            flx_aq(ispa,irain,iz) = (&
                & - caqsupp_tmp &
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
                & *merge(1.0d0,mgasx(ispg,iz),mgasx(ispg,iz)<mgasth(ispg))
            
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
                & *merge(0.0d0,1.0d0,mgasx(ispg,iz)<mgasth(ispg))
            
            
            if (iz/=nz) then 
                amx3(row,row+nsp3) = ( &
                        & -( 0.5d0*(dgas(ispg,iz)+dgas(ispg,iz+1))*(1d0)/(0.5d0*(dz(iz)+dz(min(nz,iz+1)))) &
                        & + 0.5d0*(ddgas_dmgas(ispg,ispg,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                        &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,mgasx(ispg,iz+1),mgasx(ispg,iz)<mgasth(ispg))
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 
                    amx3(row,col+nsp3) = ( &
                            & -( 0.5d0*(ddgas_dmaq(ispg,ispa,iz+1))*(mgasx(ispg,iz+1)-mgasx(ispg,iz)) &
                            &       /(0.5d0*(dz(iz)+dz(min(nz,iz+1)))))/dz(iz)*dt &
                            & ) &
                            & *merge(0.0d0,maqx(ispa,iz+1),mgasx(ispg,iz)<mgasth(ispg))
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
                    & *merge(0.0d0,mgasx(ispg,iz-1),mgasx(ispg,iz)<mgasth(ispg))
                
                do ispa = 1,nsp_aq
                    col = nsp3*(iz-1) + nsp_sld + ispa 

                    amx3(row,col-nsp3) = ( &
                        & -(- 0.5d0*(ddgas_dmaq(ispg,ispa,iz-1))*(mgasx(ispg,iz)-mgasx(ispg,iz-1)) &
                        &       /(0.5d0*(dz(iz)+dz(max(1,iz-1)))))/dz(iz)*dt  &
                        & +poro(iz)*sat(iz)*v(iz)*1d3*(-dkhgas_dmaq(ispg,ispa,iz-1)*mgasx(ispg,iz-1))/dz(iz)*dt &
                        & ) &
                        & *merge(0.0d0,maqx(ispa,iz-1),mgasx(ispg,iz)<mgasth(ispg))
                enddo 
            endif 
            
            do isps = 1,nsp_sld
                col = nsp3*(iz-1) + isps 
                amx3(row,col) = ( &
                    & -drxngas_dmsld(ispg,isps,iz)*dt &
                    & ) &
                    & *merge(1.0d0,msldx(isps,iz),mgasx(ispg,iz)<mgasth(ispg))
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
                    & ) &
                    & *merge(1.0d0,maqx(ispa,iz),mgasx(ispg,iz)<mgasth(ispg))
                
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
                    & *merge(1.0d0,mgasx(ispg2,iz),mgasx(ispg,iz)<mgasth(ispg))
                
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
            ! flx_gas(ispg,irxn_fo,iz) = -resp(iz)
            flx_gas(ispg,irxn_fo,iz) = -sum(stgas_ext(:,ispg)*rxnext(:,iz))
            flx_gas(ispg,irain,iz) = -rxngas(ispg,iz) - mgassupp(ispg,iz)
            flx_gas(ispg,ires,iz) = sum(flx_gas(ispg,:,iz))
            
            if (any(isnan(flx_gas(ispg,:,iz)))) then
                print *,flx_gas(ispg,:,iz)
            endif 
            
            ! amx3(row,:) = amx3(row,:)/alpha(iz)
            ! ymx3(row) = ymx3(row)/alpha(iz)
        enddo 

    end do 
    
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
    call calc_pH_v4( &
        & nz,maqx(findloc(chraq,'na',dim=1),:),maqx(findloc(chraq,'mg',dim=1),:),maqx(findloc(chraq,'ca',dim=1),:) &! input
        & ,so4x,mgasx(findloc(chrgas,'pco2',dim=1),:) &! input
        & ,kw,kco2,k1,k2,maqx(findloc(chraq,'si',dim=1),:),k1si,k2si,k1mg,k1mgco3,k1mghco3,k1ca,k1caco3,k1cahco3 &! input 
        & ,maqx(findloc(chraq,'al',dim=1),:),k1al,k2al,k3al,k4al &! input
        & ,prox &! output
        & ) 

    co2 = kco2*mgasx(findloc(chrgas,'pco2',dim=1),:)
    hco3 = k1*co2/prox
    co3 = k2*hco3/prox
    dic = co2 + hco3 + co3

#ifdef display      
    print *, 'silicate_dis error',error,info,iter,dt
#endif      
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
    print *,'-=-=-=-=-=-= Aq species -=-=-=-=-=-=-='
    do ispa = 1, nsp_aq
        print *, trim(adjustl(chraq(ispa))), (maqx(ispa,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Sld species -=-=-=-=-=-=-='
    do isps = 1, nsp_sld
        print *, trim(adjustl(chrsld(isps))), (msldx(isps,iz),iz=1,nz, nz/5)
    enddo 
    do isps = 1, nsp_sld
        print *, 'omega_'//trim(adjustl(chrsld(isps))), (omega(isps,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= Gas species -=-=-=-=-=-=-='
    do ispg = 1, nsp_gas
        print *, trim(adjustl(chrgas(ispg))), (mgasx(ispg,iz),iz=1,nz, nz/5)
    enddo 
    print *,'-=-=-=-=-=-= pH -=-=-=-=-=-=-='
    print *, 'ph:', (-log10(prox(iz)),iz=1,nz, nz/5)
    print *
#endif     
enddo

endsubroutine alsilicate_aq_gas_1D_v2

