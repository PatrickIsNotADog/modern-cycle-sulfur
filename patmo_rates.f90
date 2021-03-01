module patmo_rates
contains

  !***************
  subroutine computeRates(inTgas)
    use patmo_commons
    use patmo_parameters
    implicit none
    real*8,intent(in)::inTgas(cellsNumber)
    real*8::Tgas,T,invT,ntot(cellsNumber)
    integer::icell

    !total density per layer
    ntot(:) = sum(nAll(:,1:chemSpeciesNumber),2)
    open(59,file="H2SO4_photorate.txt",status="old")
    open(60,file="SO2_O_rate.txt",status="old")
    open(61,file="SO2_OH_rate.txt",status="old")
    open(62,file="CS2_OH_rate.txt",status="old")
    open(63,file="DMS_OH_rate.txt",status="old")

    !loop on cells
    do icell=1,cellsNumber
      Tgas = inTgas(icell)
      T = Tgas
      invT = 1d0/Tgas
      !CS2 + OH -> SCSOH
      krate(icell,1) = 8.00d-12

      !SCSOH + O2 -> COS + HSO2
      krate(icell,2) = 2.80d-11

      !CS + O2 -> COS + O
      krate(icell,3) = 2.9d-19

      !CS + O3 -> COS + O2
      krate(icell,4) = 3.0d-16

      !HSO2 + O2 -> HO2 + SO2
      krate(icell,5) = 3.01d-13

      !CS2E + O2 -> CS2
      krate(icell,6) = 2.5d-11

      !CS2E + N2 -> CS2
      krate(icell,7) = 2.5d-11

      !CS2E + O2 -> CS + SO2
      krate(icell,8) = 1.25d-12

      !CS + O2 -> CO + SO
      krate(icell,9) = 3.01d-18

      !SCSOH -> CS2 + OH
      krate(icell,10) = 6.16d+3

      !SO + O3 -> SO2 + O2
      krate(icell,11) = 4.50d-12*exp(-1170/T)

      !SO + O2 -> SO2 + O
      krate(icell,12) = 1.60d-13*exp(-2282/T)

      !SO + OH -> SO2 + H
      krate(icell,13) = 2.70d-11*exp(335/T)

    end do

    close(59)
    close(60)
    close(61)
    close(62)
    close(63)

  end subroutine computeRates

end module patmo_rates
