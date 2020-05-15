!!****m* ABINIT/interfaces_21_psiesta_noabirule
!! NAME
!! interfaces_21_psiesta_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/21_psiesta_noabirule
!!
!! COPYRIGHT
!! Copyright (C) 2010-2014 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_21_psiesta_noabirule

 implicit none

interface
 SUBROUTINE EGOFV(H,S,NR,E,G,Y,L,Z,A,B,RMAX,NPRIN,NNODE,DGDR)
  use defs_basis
  implicit none
  integer :: L
  integer :: NNODE
  integer :: NPRIN
  integer :: NR
  real(dp) :: A
  real(dp) :: B
  real(dp) :: DGDR
  real(dp) :: E
  real(dp) :: RMAX
  real(dp) :: Z
  real(dp) :: G(NR)
  real(dp) :: H(NR)
  real(dp) :: S(NR)
  real(dp) :: Y(NR)
 end subroutine EGOFV
end interface

interface
 SUBROUTINE YOFE( E, DE, DGDR, RMAX, H, S, Y, NR, L,&  
  &  NCOR, NNODE, Z, A, B )
  use defs_basis
  implicit none
  integer,           intent(in) :: L
  integer,           intent(in) :: NCOR
  integer,          intent(out) :: NNODE
  integer,           intent(in) :: NR
  real(dp)        ,  intent(in) :: A
  real(dp)        ,  intent(in) :: B
  real(dp)        , intent(out) :: DE
  real(dp)        ,  intent(in) :: DGDR
  real(dp)        ,  intent(in) :: E
  real(dp)        ,  intent(in) :: RMAX
  real(dp)        ,  intent(in) :: Z
  real(dp)        ,  intent(in) :: H(NR)
  real(dp)        ,  intent(in) :: S(NR)
  real(dp)        , intent(out) :: Y(NR)
 end subroutine YOFE
end interface

interface
 SUBROUTINE NRMLZG(G,S,N)
  use defs_basis
  implicit none
  integer :: N
  real(dp) :: G(N)
  real(dp) :: S(N)
 end subroutine NRMLZG
end interface

interface
 SUBROUTINE BCORGN(E,H,S,L,ZDR,Y2)
  use defs_basis
  implicit none
  integer :: L
  real(dp) :: E
  real(dp) :: Y2
  real(dp) :: ZDR
  real(dp) :: H(3)
  real(dp) :: S(3)
 end subroutine BCORGN
end interface

interface
 SUBROUTINE BCRMAX(E,DGDR,RMAX,H,S,N,YN,A,B)
  use defs_basis
  implicit none
  integer :: N
  real(dp) :: A
  real(dp) :: B
  real(dp) :: DGDR
  real(dp) :: E
  real(dp) :: RMAX
  real(dp) :: YN
  real(dp) :: H(N+1)
  real(dp) :: S(N+1)
 end subroutine BCRMAX
end interface

interface
 SUBROUTINE NUMIN(E,H,S,Y,NR,NNODE,YN,G,GSG,DY,KNK)
  use defs_basis
  implicit none
  integer :: KNK
  integer :: NNODE
  integer :: NR
  real(dp) :: DY
  real(dp) :: E
  real(dp) :: G
  real(dp) :: GSG
  real(dp) :: YN
  real(dp) :: H(NR)
  real(dp) :: S(NR)
  real(dp) :: Y(NR)
 end subroutine NUMIN
end interface

interface
 SUBROUTINE NUMOUT( E, H, S, Y, NCOR, KNK, NNODE, Y2, G, GSG, DY )
  use defs_basis
  implicit none
  integer,        intent(inout) :: KNK
  integer,           intent(in) :: NCOR
  integer,          intent(out) :: NNODE
  real(dp)        , intent(out) :: DY
  real(dp)        ,  intent(in) :: E
  real(dp)        , intent(out) :: G
  real(dp)        , intent(out) :: GSG
  real(dp)        ,  intent(in) :: Y2
  real(dp)        ,  intent(in) :: H(KNK)
  real(dp)        ,  intent(in) :: S(KNK)
  real(dp)        , intent(out) :: Y(KNK)
 end subroutine NUMOUT
end interface

interface
 SUBROUTINE VHRTRE(R2RHO,V,R,DRDI,SRDRDI,NR,A)
  use defs_basis
  implicit none
  integer :: NR
  real(dp) :: A
  real(dp) :: DRDI(NR)
  real(dp) :: R(NR)
  real(dp) :: R2RHO(NR)
  real(dp) :: SRDRDI(NR)
  real(dp) :: V(NR)
 end subroutine VHRTRE
end interface

interface
 SUBROUTINE INTEGRATOR(F,S,NP,VAL)
  use defs_basis
  implicit none
  integer :: NP
  real(dp) :: VAL
  real(dp) :: F(NP)
  real(dp) :: S(NP)
 end subroutine INTEGRATOR
end interface

interface
 subroutine atomxc( FUNCTL, AUTHOR, IREL, NR, MAXR, RMESH,&  
  &  nspin, Dens,&  
  &  EX, EC, DX, DC, VXC )
  use defs_basis
  implicit none
  integer,   intent(in) :: IREL
  integer,   intent(in) :: MAXR
  integer,   intent(in) :: NR
  integer,   intent(in) :: nspin
  character(len=*), intent(in) :: AUTHOR
  real(dp),  intent(out) :: DC
  real(dp),  intent(out) :: DX
  real(dp),  intent(out) :: EC
  real(dp),  intent(out) :: EX
  character(len=*), intent(in) :: FUNCTL
  real(dp),  intent(in) :: Dens(MAXR,nspin)
  real(dp),  intent(in) :: RMESH(MAXR)
  real(dp),  intent(out) :: VXC(MAXR,nspin)
 end subroutine atomxc
end interface

interface
 subroutine exchng( IREL, NSP, DS, EX, VX )
  use defs_basis
  implicit none
  integer, intent(in) :: irel
  integer, intent(in) :: nsp
  real(dp), intent(out) :: EX
  real(dp), intent(in) :: DS(NSP)
  real(dp), intent(out) :: VX(NSP)
 end subroutine exchng
end interface

interface
 SUBROUTINE GGAXC( AUTHOR, IREL, nspin, D, GD,&  
  &  EPSX, EPSC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  character*(*) :: AUTHOR
  real(dp) :: EPSC
  real(dp) :: EPSX
  real(dp) :: D(nspin)
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: GD(3,nspin)
 end subroutine GGAXC
end interface

interface
 SUBROUTINE LDAXC( AUTHOR, IREL, nspin, D, EPSX, EPSC, VX, VC,&  
  &  DVXDN, DVCDN )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  character*(*) :: AUTHOR
  real(dp) :: EPSC
  real(dp) :: EPSX
  real(dp) :: D(nspin)
  real(dp) :: DVCDN(nspin,nspin)
  real(dp) :: DVXDN(nspin,nspin)
  real(dp) :: VC(nspin)
  real(dp) :: VX(nspin)
 end subroutine LDAXC
end interface

interface
 SUBROUTINE PBEXC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine PBEXC
end interface

interface
 SUBROUTINE REVPBEXC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine REVPBEXC
end interface

interface
 SUBROUTINE PW91XC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine PW91XC
end interface

interface
 SUBROUTINE PW92C( nspin, Dens, EC, VC )
  use defs_basis
  implicit none
  integer :: nspin
  real(dp) :: EC
  real(dp) :: Dens(nspin)
  real(dp) :: VC(nspin)
 end subroutine PW92C
end interface

interface
 SUBROUTINE PW92XC( IREL, nspin, Dens, EPSX, EPSC, VX, VC )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EPSC
  real(dp) :: EPSX
  real(dp) :: Dens(nspin)
  real(dp) :: VC(nspin)
  real(dp) :: VX(nspin)
 end subroutine PW92XC
end interface

interface
 SUBROUTINE PZXC( IREL, NSP, DS, EX, EC, VX, VC, DVXDN, DVCDN )
  use defs_basis
  implicit none
  integer :: irel
  integer :: nsp
  real(dp) :: ec
  real(dp) :: ex
  real(dp) :: DS(NSP)
  real(dp) :: DVCDN(NSP,NSP)
  real(dp) :: DVXDN(NSP,NSP)
  real(dp) :: VC(NSP)
  real(dp) :: VX(NSP)
 end subroutine PZXC
end interface

interface
 subroutine blypxc(nspin,dens,gdens,EX,EC,&  
  &  dEXdd,dECdd,dEXdgd,dECdgd) 
  use defs_basis
  implicit none
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: dECdd(nspin)
  real(dp) :: dECdgd(3,nspin)
  real(dp) :: dEXdd(nspin)
  real(dp) :: dEXdgd(3,nspin)
  real(dp) :: dens(nspin)
  real(dp) :: gdens(3,nspin)
 end subroutine blypxc
end interface

interface
 SUBROUTINE RPBEXC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine RPBEXC
end interface

interface
 SUBROUTINE WCXC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine WCXC
end interface

interface
 SUBROUTINE PBESOLXC( IREL, nspin, Dens, GDens,&  
  &  EX, EC, DEXDD, DECDD, DEXDGD, DECDGD )
  use defs_basis
  implicit none
  integer :: IREL
  integer :: nspin
  real(dp) :: EC
  real(dp) :: EX
  real(dp) :: DECDD(nspin)
  real(dp) :: DECDGD(3,nspin)
  real(dp) :: DEXDD(nspin)
  real(dp) :: DEXDGD(3,nspin)
  real(dp) :: Dens(nspin)
  real(dp) :: GDens(3,nspin)
 end subroutine PBESOLXC
end interface

end module interfaces_21_psiesta_noabirule
!!***
