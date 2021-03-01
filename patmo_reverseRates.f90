module patmo_reverseRates
contains

  !compute reverse rates using thermochemistry polynomials
  subroutine computeReverseRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(:)
    real*8::Tgas(cellsNumber)
    real*8::lnTgas(cellsNumber)
    real*8::Tgas2(cellsNumber)
    real*8::Tgas3(cellsNumber)
    real*8::Tgas4(cellsNumber)
    real*8::invTgas(cellsNumber)
    real*8::ntot(cellsNumber)
    integer::i

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)

    !extrapolate lower and upper limits
    do i=1,cellsNumber
      Tgas(i) = max(inTgas(i),2d2)
      Tgas(i) = min(Tgas(i),5d3)
    end do

    lnTgas(:) = log(Tgas(:))
    Tgas2(:) = Tgas(:)**2
    Tgas3(:) = Tgas(:)**3
    Tgas4(:) = Tgas(:)**4
    invTgas(:) = 1d0/Tgas(:)

    !SCSOH -> CS2 + OH
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,18) = krate(i,1)*exp(8.837417d-1*(lnTgas(i)-1d0) &
!            + 6.639121d-3*Tgas(i) &
!            - 1.182399d-5*Tgas2(i) &
!            + 6.691894d-9*Tgas3(i) &
!            - 1.506943d-12*Tgas4(i) &
!            - 2.240521d4*invTgas(i) &
!            + 9.828921d0)*(1.3806488d-22*Tgas(i))**(-1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,18) = krate(i,1)*exp(2.321241d0*(lnTgas(i)-1d0) &
!            - 6.394723d-3*Tgas(i) &
!            + 7.538362d-7*Tgas2(i) &
!            - 5.96604d-11*Tgas3(i) &
!            + 2.152074d-15*Tgas4(i) &
!            - 2.306419d4*invTgas(i) &
!            + 7.693985d0)*(1.3806488d-22*Tgas(i))**(-1)
!      else
        krate(i,18) = 0d0
!      end if
    end do

    !COS + HSO2 -> SCSOH + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,19) = krate(i,2)*exp(3.553689d0*(lnTgas(i)-1d0) &
            - 1.350409d-2*Tgas(i) &
            + 1.429369d-5*Tgas2(i) &
            - 7.152137d-9*Tgas3(i) &
            + 1.503259d-12*Tgas4(i) &
            - 4.066505d4*invTgas(i) &
            - 1.808164d1)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,19) = krate(i,2)*exp(-1.642933d0*(lnTgas(i)-1d0) &
            + 5.788839d-3*Tgas(i) &
            - 6.630917d-7*Tgas2(i) &
            + 5.234627d-11*Tgas3(i) &
            - 1.859325d-15*Tgas4(i) &
            - 4.088933d4*invTgas(i) &
            + 2.80749d0)
      else
        krate(i,19) = 0d0
      end if
    end do

    !COS + O -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,20) = krate(i,3)*exp(2.573447d0*(lnTgas(i)-1d0) &
            - 9.982074d-3*Tgas(i) &
            + 7.16588d-6*Tgas2(i) &
            - 3.355992d-9*Tgas3(i) &
            + 6.904257d-13*Tgas4(i) &
            - 2.038875d4*invTgas(i) &
            - 7.526717d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,20) = krate(i,3)*exp(-4.876406d-1*(lnTgas(i)-1d0) &
            - 3.447248d-4*Tgas(i) &
            + 6.608958d-8*Tgas2(i) &
            - 7.10943d-12*Tgas3(i) &
            + 3.52615d-16*Tgas4(i) &
            - 2.072572d4*invTgas(i) &
            + 5.698036d0)
      else
        krate(i,20) = 0d0
      end if
    end do

    !COS + O2 -> CS + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,21) = krate(i,4)*exp(1.584184d0*(lnTgas(i)-1d0) &
            - 7.598104d-3*Tgas(i) &
            + 7.298724d-6*Tgas2(i) &
            - 4.114045d-9*Tgas3(i) &
            + 9.597224d-13*Tgas4(i) &
            - 6.75034d4*invTgas(i) &
            - 4.507659d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,21) = krate(i,4)*exp(7.064366d0*(lnTgas(i)-1d0) &
            - 6.980988d-3*Tgas(i) &
            + 1.443677d-6*Tgas2(i) &
            - 1.577886d-10*Tgas3(i) &
            + 6.762342d-15*Tgas4(i) &
            - 6.505927d4*invTgas(i) &
            - 3.709273d1)
      else
        krate(i,21) = 0d0
      end if
    end do

    !HO2 + SO2 -> HSO2 + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,22) = krate(i,5)*exp(-4.568213d-1*(lnTgas(i)-1d0) &
            + 4.389574d-3*Tgas(i) &
            - 3.938993d-6*Tgas2(i) &
            + 2.014456d-9*Tgas3(i) &
            - 4.35746d-13*Tgas4(i) &
            - 5.791471d3*invTgas(i) &
            + 2.080435d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.5d3) then
        krate(i,22) = krate(i,5)*exp(5.001109d-1*(lnTgas(i)-1d0) &
            - 1.699002d-4*Tgas(i) &
            - 6.723714d-9*Tgas2(i) &
            + 2.481935d-12*Tgas3(i) &
            - 1.41653d-16*Tgas4(i) &
            - 5.822978d3*invTgas(i) &
            - 1.337804d0)
      else
        krate(i,22) = 0d0
      end if
    end do

    !CS2 -> CS2E + O2
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,23) = krate(i,6)*exp(3.782456d0*(lnTgas(i)-1d0) &
!            - 1.498367d-3*Tgas(i) &
!            + 1.641217d-6*Tgas2(i) &
!            - 8.067746d-10*Tgas3(i) &
!            + 1.621864d-13*Tgas4(i) &
!            + 1.063944d3*invTgas(i) &
!            + 3.657676d0)*(1.3806488d-22*Tgas(i))**(-1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,23) = krate(i,6)*exp(3.660961d0*(lnTgas(i)-1d0) &
!            + 3.281829d-4*Tgas(i) &
!            - 2.352494d-8*Tgas2(i) &
!            + 1.714983d-12*Tgas3(i) &
!            - 6.495672d-17*Tgas4(i) &
!            + 1.215977d3*invTgas(i) &
!            + 3.415363d0)*(1.3806488d-22*Tgas(i))**(-1)
!      else
        krate(i,23) = 0d0
!      end if
    end do

    !CS2 -> CS2E + N2
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,24) = krate(i,7)*exp(3.531005d0*(lnTgas(i)-1d0) &
!            - 6.183049d-5*Tgas(i) &
!            - 8.383324d-8*Tgas2(i) &
!            + 2.029422d-10*Tgas3(i) &
!            - 7.044062d-14*Tgas4(i) &
!            + 1.046976d3*invTgas(i) &
!            + 2.96747d0)*(1.3806488d-22*Tgas(i))**(-1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,24) = krate(i,7)*exp(2.952576d0*(lnTgas(i)-1d0) &
!            + 6.984502d-4*Tgas(i) &
!            - 8.210527d-8*Tgas2(i) &
!            + 6.550085d-12*Tgas3(i) &
!            - 2.303776d-16*Tgas4(i) &
!            + 9.239487d2*invTgas(i) &
!            + 5.871888d0)*(1.3806488d-22*Tgas(i))**(-1)
!      else
        krate(i,24) = 0d0
!      end if
    end do

    !CS + SO2 -> CS2E + O2
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,25) = krate(i,8)*exp(-1.451291d0*(lnTgas(i)-1d0) &
!            + 7.972313d-3*Tgas(i) &
!            - 6.985411d-6*Tgas2(i) &
!            + 3.721109d-9*Tgas3(i) &
!            - 8.390358d-13*Tgas4(i) &
!            - 1.624584d4*invTgas(i) &
!            + 3.123155d0)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,25) = krate(i,8)*exp(4.561796d-1*(lnTgas(i)-1d0) &
!            - 3.051946d-5*Tgas(i) &
!            + 9.916702d-9*Tgas2(i) &
!            - 1.934438d-14*Tgas3(i) &
!            - 4.758353d-17*Tgas4(i) &
!            - 1.615798d4*invTgas(i) &
!            - 4.34393d0)
!      else
        krate(i,25) = 0d0
!      end if
    end do

    !CO + SO -> CS + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,26) = krate(i,9)*exp(3.155756d-1*(lnTgas(i)-1d0) &
            - 1.58134d-3*Tgas(i) &
            + 1.611175d-6*Tgas2(i) &
            - 8.785323d-10*Tgas3(i) &
            + 1.937104d-13*Tgas4(i) &
            - 4.620286d4*invTgas(i) &
            - 1.667224d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,26) = krate(i,9)*exp(4.131292d-1*(lnTgas(i)-1d0) &
            - 1.708392d-4*Tgas(i) &
            + 1.567545d-8*Tgas2(i) &
            - 1.309359d-12*Tgas3(i) &
            + 9.876752d-17*Tgas4(i) &
            - 4.602858d4*invTgas(i) &
            - 2.916436d0)
      else
        krate(i,26) = 0d0
      end if
    end do

    !CS2 + OH -> SCSOH
    do i=1,cellsNumber
!      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
!        krate(i,27) = krate(i,10)*exp(-8.837417d-1*(lnTgas(i)-1d0) &
!            - 6.639121d-3*Tgas(i) &
!            + 1.182399d-5*Tgas2(i) &
!            - 6.691894d-9*Tgas3(i) &
!            + 1.506943d-12*Tgas4(i) &
!            + 2.240521d4*invTgas(i) &
!            - 9.828921d0)*(1.3806488d-22*Tgas(i))**(1)
!      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
!        krate(i,27) = krate(i,10)*exp(-2.321241d0*(lnTgas(i)-1d0) &
!            + 6.394723d-3*Tgas(i) &
!            - 7.538362d-7*Tgas2(i) &
!            + 5.96604d-11*Tgas3(i) &
!            - 2.152074d-15*Tgas4(i) &
!            + 2.306419d4*invTgas(i) &
!            - 7.693985d0)*(1.3806488d-22*Tgas(i))**(1)
!      else
        krate(i,27) = 0d0
!      end if
    end do

    !SO2 + O2 -> SO + O3
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,28) = krate(i,11)*exp(-4.312865d-1*(lnTgas(i)-1d0) &
            + 2.22883d-4*Tgas(i) &
            + 1.19644d-6*Tgas2(i) &
            - 1.100242d-9*Tgas3(i) &
            + 3.180969d-13*Tgas4(i) &
            - 5.339333d4*invTgas(i) &
            + 3.021177d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,28) = krate(i,11)*exp(7.254038d0*(lnTgas(i)-1d0) &
            - 6.945426d-3*Tgas(i) &
            + 1.461383d-6*Tgas2(i) &
            - 1.595621d-10*Tgas3(i) &
            + 6.770763d-15*Tgas4(i) &
            - 5.076969d4*invTgas(i) &
            - 3.873146d1)
      else
        krate(i,28) = 0d0
      end if
    end do

    !SO2 + O -> SO + O2
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,29) = krate(i,12)*exp(5.579769d-1*(lnTgas(i)-1d0) &
            - 2.161087d-3*Tgas(i) &
            + 1.063596d-6*Tgas2(i) &
            - 3.421897d-10*Tgas3(i) &
            + 4.880018d-14*Tgas4(i) &
            - 6.278683d3*invTgas(i) &
            + 2.11912d-3)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,29) = krate(i,12)*exp(-2.979689d-1*(lnTgas(i)-1d0) &
            - 3.091634d-4*Tgas(i) &
            + 8.379577d-8*Tgas2(i) &
            - 8.882901d-12*Tgas3(i) &
            + 3.610358d-16*Tgas4(i) &
            - 6.436141d3*invTgas(i) &
            + 4.059304d0)
      else
        krate(i,29) = 0d0
      end if
    end do

    !SO2 + H -> SO + OH
    do i=1,cellsNumber
      if(Tgas(i)<1d3.and.Tgas(i).ge.2d2) then
        krate(i,30) = krate(i,13)*exp(1.435772d0*(lnTgas(i)-1d0) &
            - 3.502913d-3*Tgas(i) &
            + 1.298996d-6*Tgas2(i) &
            - 3.693508d-10*Tgas3(i) &
            + 6.04065d-14*Tgas4(i) &
            - 1.436012d4*invTgas(i) &
            - 1.260939d0)
      elseif(Tgas(i).ge.1d3.and.Tgas(i).le.6d3) then
        krate(i,30) = krate(i,13)*exp(-1.076762d0*(lnTgas(i)-1d0) &
            - 9.729794d-5*Tgas(i) &
            + 5.762229d-8*Tgas2(i) &
            - 6.67916d-12*Tgas3(i) &
            + 2.808699d-16*Tgas4(i) &
            - 1.510228d4*invTgas(i) &
            + 1.185787d1)
      else
        krate(i,30) = 0d0
      end if
    end do

  end subroutine computeReverseRates

end module patmo_reverseRates
