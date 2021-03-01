module patmo_ode
contains
  subroutine fex(neq,tt,nin,dy)
    use patmo_commons
    use patmo_constants
    use patmo_parameters
    use patmo_utils
    implicit none
    integer,intent(in)::neq
    real*8,intent(in)::tt,nin(neqAll)
    real*8,intent(out)::dy(neqAll)
    real*8::d_hp(cellsNumber,speciesNumber)
    real*8::d_hm(cellsNumber,speciesNumber)
    real*8::k_hp(cellsNumber)
    real*8::k_hm(cellsNumber)
    real*8::dzz_hp(cellsNumber),dzz_hm(cellsNumber)
    real*8::kzz_hp(cellsNumber),kzz_hm(cellsNumber)
    real*8::prem(cellsNumber)
    real*8::n(cellsNumber,speciesNumber)
    real*8::dn(cellsNumber,speciesNumber)
    real*8::Tgas(cellsNumber)
    real*8::n_p(cellsNumber,speciesNumber)
    real*8::n_m(cellsNumber,speciesNumber)
    real*8::m(speciesNumber),ngas(cellsNumber)
    real*8::ngas_hp(cellsNumber),ngas_hm(cellsNumber)
    real*8::ngas_p(cellsNumber),ngas_m(cellsNumber)
    real*8::Tgas_hp(cellsNumber),Tgas_hm(cellsNumber)
    real*8::Tgas_p(cellsNumber),Tgas_m(cellsNumber)
    real*8::ngas_hpp(cellsNumber)
    real*8::ngas_hmm(cellsNumber)
    real*8::ngas_hpz(cellsNumber)
    real*8::ngas_hmz(cellsNumber)
    real*8::therm_hp(cellsNumber)
    real*8::therm_hm(cellsNumber)
    real*8::dzzh_hp(cellsNumber)
    real*8::dzzh_hm(cellsNumber)
    real*8::iTgas_hp(cellsNumber)
    real*8::iTgas_hm(cellsNumber)
    integer::i,j

    !get mass of individual species
    m(:) = getSpeciesMass()

    !roll chemistry
    do i=1,speciesNumber
      n(:,i) = nin((i-1)*cellsNumber+1:(i*cellsNumber))
    end do

    !local copy of Tgas
    Tgas(:) = nin((positionTgas-1)*cellsNumber+1:(positionTgas*cellsNumber))
    ngas(:) = nTotAll(:)

    !forward grid points
    do j=1,cellsNumber-1
      dzz_hp(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j+1))
      kzz_hp(j) = .5d0*(eddyKzz(j)+eddyKzz(j+1))
      Tgas_hp(j) = .5d0*(Tgas(j)+Tgas(j+1))
      Tgas_p(j) = Tgas(j+1)
      ngas_p(j) = ngas(j+1)
      ngas_hp(j) = .5d0*(ngas(j)+ngas(j+1))
      n_p(j,:) = n(j+1,:)
    end do

    !forward grid points: boundary conditions
    dzz_hp(cellsNumber) = 0d0
    kzz_hp(cellsNumber) = 0d0
    Tgas_hp(cellsNumber) = Tgas_hp(cellsNumber-1)
    Tgas_p(cellsNumber) = Tgas_p(cellsNumber-1)
    ngas_p(cellsNumber) = ngas_p(cellsNumber-1)
    ngas_hp(cellsNumber) = ngas_hp(cellsNumber-1)
    n_p(cellsNumber,:) = n_p(cellsNumber-1,:)

    !bakcward grid points
    do j=2,cellsNumber
      dzz_hm(j) = .5d0*(diffusionDzz(j)+diffusionDzz(j-1))
      kzz_hm(j) = .5d0*(eddyKzz(j)+eddyKzz(j-1))
      Tgas_hm(j) = .5d0*(Tgas(j)+Tgas(j-1))
      Tgas_m(j) = Tgas(j-1)
      ngas_m(j) = ngas(j-1)
      ngas_hm(j) = .5d0*(ngas(j)+ngas(j-1))
      n_m(j,:) = n(j-1,:)
    end do

    !backward grid points: boundary conditions
    dzz_hm(1) = 0d0
    kzz_hm(1) = 0d0
    Tgas_hm(1) = Tgas_hm(2)
    Tgas_m(1) = Tgas_m(2)
    ngas_m(1) = ngas_m(2)
    ngas_hm(1) = ngas_hm(2)
    n_m(1,:) = n_m(2,:)

    !eqn.24 of Rimmer+Helling (2015), http://arxiv.org/abs/1510.07052
    therm_hp(:) = thermalDiffusionFactor/Tgas_hp(:)*(Tgas_p(:)-Tgas(:))
    therm_hm(:) = thermalDiffusionFactor/Tgas_hm(:)*(Tgas_m(:)-Tgas(:))
    dzzh_hp(:) = 0.5d0*dzz_hp(:)*idh2(:)
    dzzh_hm(:) = 0.5d0*dzz_hm(:)*idh2(:)
    iTgas_hp(:) = 1d0/Tgas_hp(:)
    iTgas_hm(:) = 1d0/Tgas_hm(:)
    do i=1,speciesNumber
      prem(:) = (meanMolecularMass-m(i))*gravity/kboltzmann*gridSpace(:)
      d_hp(:,i) =  dzzh_hp(:) &
          * (prem(:)*iTgas_hp(:) &
          - therm_hp(:))
      d_hm(:,i) = dzzh_hm(:) &
          * (prem(:)*iTgas_hm(:) &
          - therm_hm(:))
    end do

    k_hp(:) = (kzz_hp(:)+dzz_hp(:))*idh2(:)
    k_hm(:) = (kzz_hm(:)+dzz_hm(:))*idh2(:)

    dn(:,:) = 0d0

    dn(:,patmo_idx_CS2) = &
        - krate(:,1)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        + krate(:,6)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,7)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        + krate(:,10)*n(:,patmo_idx_SCSOH) &
        - krate(:,17)*n(:,patmo_idx_CS2) &
        + krate(:,18)*n(:,patmo_idx_SCSOH) &
        - krate(:,23)*n(:,patmo_idx_CS2) &
        - krate(:,24)*n(:,patmo_idx_CS2) &
        - krate(:,27)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_CS2E) = &
        - krate(:,6)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,7)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
        - krate(:,8)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,17)*n(:,patmo_idx_CS2) &
        + krate(:,23)*n(:,patmo_idx_CS2) &
        + krate(:,24)*n(:,patmo_idx_CS2) &
        + krate(:,25)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_SO) = &
        + krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,11)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        - krate(:,12)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        - krate(:,13)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        + krate(:,15)*n(:,patmo_idx_SO2) &
        - krate(:,26)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
        + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        + krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        + krate(:,30)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H)

    dn(:,patmo_idx_O) = 0d0
    !    + krate(:,3)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,12)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    + krate(:,14)*n(:,patmo_idx_O3) &
    !    + krate(:,15)*n(:,patmo_idx_SO2) &
    !    + krate(:,16)*n(:,patmo_idx_O2) &
    !    + krate(:,16)*n(:,patmo_idx_O2) &
    !    - krate(:,20)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
    !    - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)

    dn(:,patmo_idx_SO2) = &
        + krate(:,5)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        + krate(:,8)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        + krate(:,11)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
        + krate(:,12)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
        + krate(:,13)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
        - krate(:,15)*n(:,patmo_idx_SO2) &
        - krate(:,22)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
        - krate(:,25)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
        - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
        - krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O) &
        - krate(:,30)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H)

    dn(:,patmo_idx_HO2) = 0d0
    !    + krate(:,5)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    - krate(:,22)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_SCSOH) = &
        + krate(:,1)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
        - krate(:,2)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,10)*n(:,patmo_idx_SCSOH) &
        - krate(:,18)*n(:,patmo_idx_SCSOH) &
        + krate(:,19)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        + krate(:,27)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH)

    dn(:,patmo_idx_OH) = 0d0
    !    - krate(:,1)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
    !    + krate(:,10)*n(:,patmo_idx_SCSOH) &
    !    - krate(:,13)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    + krate(:,18)*n(:,patmo_idx_SCSOH) &
    !    - krate(:,27)*n(:,patmo_idx_CS2)*n(:,patmo_idx_OH) &
    !    + krate(:,30)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H)

    dn(:,patmo_idx_COS) = &
        + krate(:,2)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        + krate(:,3)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,4)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        - krate(:,19)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        - krate(:,20)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        - krate(:,21)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_HSO2) = &
        + krate(:,2)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
        - krate(:,5)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
        - krate(:,19)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
        + krate(:,22)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2)

    dn(:,patmo_idx_CO) = 0d0
    !    + krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    - krate(:,26)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_CS) = &
        - krate(:,3)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        - krate(:,4)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
        + krate(:,8)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
        - krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
        + krate(:,20)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
        + krate(:,21)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
        - krate(:,25)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
        + krate(:,26)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO)

    dn(:,patmo_idx_O2) = 0d0
    !    - krate(:,2)*n(:,patmo_idx_SCSOH)*n(:,patmo_idx_O2) &
    !    - krate(:,3)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,4)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
    !    - krate(:,5)*n(:,patmo_idx_HSO2)*n(:,patmo_idx_O2) &
    !    - krate(:,6)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
    !    - krate(:,8)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_O2) &
    !    - krate(:,9)*n(:,patmo_idx_CS)*n(:,patmo_idx_O2) &
    !    + krate(:,11)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
    !    - krate(:,12)*n(:,patmo_idx_SO)*n(:,patmo_idx_O2) &
    !    + krate(:,14)*n(:,patmo_idx_O3) &
    !    - krate(:,16)*n(:,patmo_idx_O2) &
    !    + krate(:,19)*n(:,patmo_idx_COS)*n(:,patmo_idx_HSO2) &
    !    + krate(:,20)*n(:,patmo_idx_COS)*n(:,patmo_idx_O) &
    !    - krate(:,21)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
    !    + krate(:,22)*n(:,patmo_idx_HO2)*n(:,patmo_idx_SO2) &
    !    + krate(:,23)*n(:,patmo_idx_CS2) &
    !    + krate(:,25)*n(:,patmo_idx_CS)*n(:,patmo_idx_SO2) &
    !    + krate(:,26)*n(:,patmo_idx_CO)*n(:,patmo_idx_SO) &
    !    - krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2) &
    !    + krate(:,29)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O)

    dn(:,patmo_idx_H) = 0d0
    !    + krate(:,13)*n(:,patmo_idx_SO)*n(:,patmo_idx_OH) &
    !    - krate(:,30)*n(:,patmo_idx_SO2)*n(:,patmo_idx_H)

    dn(:,patmo_idx_O3) = 0d0
    !    - krate(:,4)*n(:,patmo_idx_CS)*n(:,patmo_idx_O3) &
    !    - krate(:,11)*n(:,patmo_idx_SO)*n(:,patmo_idx_O3) &
    !    - krate(:,14)*n(:,patmo_idx_O3) &
    !    + krate(:,21)*n(:,patmo_idx_COS)*n(:,patmo_idx_O2) &
    !    + krate(:,28)*n(:,patmo_idx_SO2)*n(:,patmo_idx_O2)

    dn(:,patmo_idx_N2) = 0d0
    !    - krate(:,7)*n(:,patmo_idx_CS2E)*n(:,patmo_idx_N2) &
    !    + krate(:,24)*n(:,patmo_idx_CS2)

    ngas_hpp(:) = ngas_hp(:)/ngas_p(:)
    ngas_hpz(:) = ngas_hp(:)/ngas(:)
    ngas_hmm(:) = ngas_hm(:)/ngas_m(:)
    ngas_hmz(:) = ngas_hm(:)/ngas(:)

!!This is the calculation code for transport of all species
!    do i=1,chemSpeciesNumber
!      dn(:,i) = dn(:,i) &
!          + (k_hp(:)-d_hp(:,i)) * ngas_hpp(:) * n_p(:,i) &
!          - ((k_hp(:)+d_hp(:,i)) * ngas_hpz(:) &
!          + (k_hm(:)-d_hm(:,i)) * ngas_hmz(:)) * n(:,i) &
!          + (k_hm(:)+d_hm(:,i)) * ngas_hmm(:) * n_m(:,i)
!    end do

!Notice that transport of sulfur compounds only 
       do i=1,60 !COS
        dn(i,1) = dn(i,1) &
             + (k_hp(i)-d_hp(i,1)) * ngas_hpp(i) * n_p(i,1) &
             - ((k_hp(i)+d_hp(i,1)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,1)) * ngas_hmz(i)) * n(i,1) &
            + (k_hm(i)+d_hm(i,1)) * ngas_hmm(i) * n_m(i,1)
       end do    
      
       do i=1,60 !SCSOH
        dn(i,3) = dn(i,3) &
             + (k_hp(i)-d_hp(i,3)) * ngas_hpp(i) * n_p(i,3) &
             - ((k_hp(i)+d_hp(i,3)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,3)) * ngas_hmz(i)) * n(i,3) &
            + (k_hm(i)+d_hm(i,3)) * ngas_hmm(i) * n_m(i,3)
       end do       
      
       do i=1,60 !SO2
        dn(i,7) = dn(i,7) &
             + (k_hp(i)-d_hp(i,7)) * ngas_hpp(i) * n_p(i,7) &
             - ((k_hp(i)+d_hp(i,7)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,7)) * ngas_hmz(i)) * n(i,7) &
            + (k_hm(i)+d_hm(i,7)) * ngas_hmm(i) * n_m(i,7)
       end do        
 
       do i=1,60 !CS2E
        dn(i,9) = dn(i,9) &
             + (k_hp(i)-d_hp(i,9)) * ngas_hpp(i) * n_p(i,9) &
             - ((k_hp(i)+d_hp(i,9)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,9)) * ngas_hmz(i)) * n(i,9) &
            + (k_hm(i)+d_hm(i,9)) * ngas_hmm(i) * n_m(i,9)
       end do   
      
       do i=11,12 !CS, SO
        do j=1,60
         dn(j,i) = dn(j,i) &
             + (k_hp(j)-d_hp(j,i)) * ngas_hpp(j) * n_p(j,i) &
             - ((k_hp(j)+d_hp(j,i)) * ngas_hpz(j) &
             + (k_hm(j)-d_hm(j,i)) * ngas_hmz(j)) * n(j,i) &
             + (k_hm(j)+d_hm(j,i)) * ngas_hmm(j) * n_m(j,i)
        end do
      end do
       
      do i=1,60 !HSO2
        dn(i,14) = dn(i,14) &
             + (k_hp(i)-d_hp(i,14)) * ngas_hpp(i) * n_p(i,14) &
             - ((k_hp(i)+d_hp(i,14)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,14)) * ngas_hmz(i)) * n(i,14) &
            + (k_hm(i)+d_hm(i,14)) * ngas_hmm(i) * n_m(i,14)
       end do    
  
      do i=1,60 !CS2
        dn(i,16) = dn(i,16) &
             + (k_hp(i)-d_hp(i,16)) * ngas_hpp(i) * n_p(i,16) &
             - ((k_hp(i)+d_hp(i,16)) * ngas_hpz(i) &
             + (k_hm(i)-d_hm(i,16)) * ngas_hmz(i)) * n(i,16) &
            + (k_hm(i)+d_hm(i,16)) * ngas_hmm(i) * n_m(i,16)
       end do  

   !Calculation of troposphere-stratosphere calculation (not sure if really works!	
     !to stratosphere
       ! write(43,*)  (k_hm(13)+d_hm(13,1)) * ngas_hmm(13) * n_m(13,1) 
   
     !to troposphere 
         ! write(44,*)  (k_hm(13)-d_hm(13,1)) * ngas_hmz(13)* n(13,1) 	
   !!!!!End of transport section
   
     !emission
     dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) + 1.03d7/1d5           !OCS 1.3 Tg/y
      dn(1,patmo_idx_CS2) = dn(1,patmo_idx_CS2) + 2.20d7/1d5           !CS2 124-250 Khalil
  !   dn(1,patmo_idx_H2S) = dn(1,patmo_idx_H2S) + 1.40d8/1d5           !H2S 903.5-1252.6
     dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) + 7.50d8/1d5           !SO2 7196-7715
  !   dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) + 2.50d8/1d5   !DMS
   
   
     !dry deposition
     dn(1,patmo_idx_COS) = dn(1,patmo_idx_COS) - 1.585d-8*n(1,1)             !OCS  lifetime 2.0y
       dn(1,patmo_idx_SO2) = dn(1,patmo_idx_SO2) - 1.355d-7*n(1,7)            !SO2  lifetime 5.0d
  !     dn(1,patmo_idx_CH3SCH3) = dn(1,patmo_idx_CH3SCH3) - 7.940d-6*n(1,33)    !DMS  lifetime 0.5d
   
     !aerosol formation
  !    do i=13,34
  !       if (va(i) <= n(i,26) .and. pa(i) >= n(i,26)) then
  !      dn(i,patmo_idx_H2SO4) = dn(i,patmo_idx_H2SO4)- (n(i,26)-va(i))
  !        dn(i,patmo_idx_SO4) = dn(i,patmo_idx_SO4)+ (n(i,26)-va(i))
  !       end if
  !    end do	
   
     !gravity settling SO4 Aerosol (JAM-Kasten-1968,r=0.3)
  !    do j=60,2,-1
  !         dn(j,patmo_idx_SO4) = dn(j,patmo_idx_SO4)-gd(j)*n(j,patmo_idx_SO4)
  !       dn(j-1,patmo_idx_SO4) = dn(j-1,patmo_idx_SO4)+gd(j)*n(j,patmo_idx_SO4)
  !       end do  
  !      dn(1,19) = dn(1,19)-gd(1)*n(1,19)
   
     !wet deposition
     do j=12,2,-1
        do i=1,chemSpeciesNumber
            dn(j,i) = dn(j,i)-wetdep(j,i)*n(j,i)
          dn(j-1,i) = dn(j-1,i)+wetdep(j,i)*n(j,i)
        end do   												
       end do
       do i=1,chemSpeciesNumber
         dn(1,i) = dn(1,i)-wetdep(1,i)*n(1,i)
       end do  
   
    !DMS → SO2 (96%)
     ! do i=1,60
     !     dn(i,patmo_idx_CH3SCH3) = dn(i,patmo_idx_CH3SCH3)- (n(i,30)*0.96d0)
     !     dn(i,patmo_idx_SO2) = dn(i,patmo_idx_SO2)+ (n(i,30)*0.96d0)
     ! end do
   ! end if

    !unroll chemistry
    dy(:) = 0d0
    do i=1,speciesNumber
      dy((i-1)*cellsNumber+1:(i*cellsNumber)) = dn(:,i)
    end do

  end subroutine fex
end module patmo_ode
