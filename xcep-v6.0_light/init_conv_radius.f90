  ! 
  !============ covalent radius of atoms ============
  ! 
  subroutine init_conv_radius( melement, cov_radius)
    implicit none 
    integer :: melement
    real(8) :: cov_radius(melement)
    !
    ! numbers from wiki 
    ! http://pubs.rsc.org/en/Content/ArticleLanding/2008/DT/b801115j#!divAbstract
    ! Dalton Trans., 2008, 2832-2838
    !
     cov_radius(1) = 0.31d0
     cov_radius(2) = 0.28d0
     cov_radius(3) = 1.28d0
     cov_radius(4) = 0.96d0
     cov_radius(5) = 0.84d0
     cov_radius(6) = 0.76d0  
     cov_radius(7) = 0.71d0
     cov_radius(8) = 0.66d0
     cov_radius(9) = 0.57d0
     cov_radius(10) = 0.58d0
     cov_radius(11) = 1.66d0 
     cov_radius(12) = 1.41d0
     cov_radius(13) = 1.21d0 
     cov_radius(14) = 1.11d0
     cov_radius(15) = 1.07d0
     cov_radius(16) = 1.05d0 
     cov_radius(17) = 1.02d0
     cov_radius(18) = 1.06d0 
     cov_radius(19) = 2.03d0
     cov_radius(20) = 1.76d0 
     cov_radius(21) = 1.70d0 
     cov_radius(22) = 1.60d0
     cov_radius(23) = 1.53d0 
     cov_radius(24) = 1.39d0
     cov_radius(25) = 1.50d0  ! l.s
     cov_radius(26) = 1.32d0  ! l.s
     cov_radius(27) = 1.26d0  ! l.s
     cov_radius(28) = 1.24d0 
     cov_radius(29) = 1.32d0
     cov_radius(30) = 1.22d0
     cov_radius(31) = 1.22d0 
     cov_radius(32) = 1.20d0 
     cov_radius(33) = 1.19d0
     cov_radius(34) = 1.20d0
     cov_radius(35) = 1.20d0
     cov_radius(36) = 1.16d0
     cov_radius(37) = 2.20d0
     cov_radius(38) = 1.95d0
     cov_radius(39) = 1.90d0
     cov_radius(40) = 1.75d0
     cov_radius(41) = 1.64d0
     cov_radius(42) = 1.54d0
     cov_radius(43) = 1.47d0
     cov_radius(44) = 1.46d0
     cov_radius(45) = 1.42d0
     cov_radius(46) = 1.39d0
     cov_radius(47) = 1.45d0
     cov_radius(48) = 1.44d0
     cov_radius(49) = 1.42d0
     cov_radius(50) = 1.39d0
     cov_radius(51) = 1.39d0
     cov_radius(52) = 1.38d0
     cov_radius(53) = 1.39d0
     cov_radius(54) = 1.40d0
     cov_radius(55) = 2.44d0
     cov_radius(56) = 2.15d0
     cov_radius(57) = 2.07d0
     cov_radius(58) = 2.04d0
     cov_radius(59) = 2.03d0
     cov_radius(60) = 2.01d0
     cov_radius(61) = 1.99d0
     cov_radius(62) = 1.98d0
     cov_radius(63) = 1.98d0
     cov_radius(64) = 1.96d0
     cov_radius(65) = 1.94d0
     cov_radius(66) = 1.92d0
     cov_radius(67) = 1.92d0
     cov_radius(68) = 1.89d0
     cov_radius(69) = 1.90d0
     cov_radius(70) = 1.87d0
     cov_radius(71) = 1.87d0
     cov_radius(72) = 1.75d0
     cov_radius(73) = 1.70d0
     cov_radius(74) = 1.62d0
     cov_radius(75) = 1.51d0
     cov_radius(76) = 1.44d0
     cov_radius(77) = 1.41d0
     cov_radius(78) = 1.36d0
     cov_radius(79) = 1.36d0
     cov_radius(80) = 1.32d0
     cov_radius(81) = 1.45d0
     cov_radius(82) = 1.46d0
     cov_radius(83) = 1.48d0
     cov_radius(84) = 1.40d0
     cov_radius(85) = 1.50d0
     cov_radius(86) = 1.50d0
     cov_radius(87) = 2.60d0
     cov_radius(88) = 2.21d0
     cov_radius(89) = 2.15d0
     cov_radius(90) = 2.06d0
     cov_radius(91) = 2.00d0
     cov_radius(92) = 1.96d0
     cov_radius(93) = 1.90d0
     cov_radius(94) = 1.87d0
     cov_radius(95) = 1.80d0
     cov_radius(96) = 1.69d0


  end subroutine init_conv_radius

