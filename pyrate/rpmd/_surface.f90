!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   PyRate - Python tools for computing chemical reaction rates
!
!   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
!                         Yury V. Suleimanov (ysuleyma@mit.edu)
!                         William H. Green (whgreen@mit.edu)
!
!   Permission is hereby granted, free of charge, to any person obtaining a 
!   copy of this software and associated documentation files (the "Software"), 
!   to deal in the Software without restriction, including without limitation
!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!   and/or sell copies of the Software, and to permit persons to whom the 
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be included in
!   all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
!   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
!   DEALINGS IN THE SOFTWARE. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! Return the value of the dividing surface function for each of the
! equivalent transition states that define the dividing surface.
! Parameters:
!   position - A 3 x Natoms array of atomic positions
!   Natoms - The number of atoms
!   Nts - The number of equivalent transition states that define the dividing surface
!   forming_bonds - An array listing the pairs of indices of each forming bond in each transition state
!   forming_bond_lengths - An array listing the lengths of each forming bond in each transition state
!   number_of_forming_bonds - The number of bonds being formed by the reaction
!   breaking_bonds - An array listing the pairs of indices of each breaking bond in each transition state
!   breaking_bond_lengths - An array listing the lengths of each breaking bond in each transition state
!   number_of_breaking_bonds - The number of bonds being broken by the reaction
! Returns:
!   values - The value of the dividing surface function for each transition state
subroutine transition_state_evaluate_all(position, Natoms, Nts, &
    forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
    breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, values)

    implicit none
    integer, intent(in) :: Natoms, Nts
    integer, intent(in) :: number_of_forming_bonds, number_of_breaking_bonds
    double precision, intent(in) :: position(3,Natoms)
    integer, intent(in) :: forming_bonds(Nts,number_of_forming_bonds,2)
    integer, intent(in) :: breaking_bonds(Nts,number_of_breaking_bonds,2)
    double precision, intent(in) :: forming_bond_lengths(Nts,number_of_forming_bonds)
    double precision, intent(in) :: breaking_bond_lengths(Nts,number_of_breaking_bonds)
    double precision, dimension(Nts), intent(out) :: values

    integer :: m, n, atom1, atom2
    double precision :: Rx, Ry, Rz, R

    values(:) = 0.0

    do n = 1, Nts

        do m = 1, number_of_forming_bonds
            atom1 = forming_bonds(n,m,1)
            atom2 = forming_bonds(n,m,2)
            Rx = position(1,atom1) - position(1,atom2)
            Ry = position(2,atom1) - position(2,atom2)
            Rz = position(3,atom1) - position(3,atom2)
            R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            values(n) = values(n) + (forming_bond_lengths(n,m) - R)
        end do

        do m = 1, number_of_breaking_bonds
            atom1 = breaking_bonds(n,m,1)
            atom2 = breaking_bonds(n,m,2)
            Rx = position(1,atom1) - position(1,atom2)
            Ry = position(2,atom1) - position(2,atom2)
            Rz = position(3,atom1) - position(3,atom2)
            R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            values(n) = values(n) - (breaking_bond_lengths(n,m) - R)
        end do

    end do

end subroutine transition_state_evaluate_all

! Return the value of the transition state dividing surface function. This is
! the maximum of the individual values for each equivalent transition state,
! as determined by transition_state_evaluate_all().
! Parameters:
!   position - A 3 x Natoms array of atomic positions
!   Natoms - The number of atoms
!   Nts - The number of equivalent transition states that define the dividing surface
!   forming_bonds - An array listing the pairs of indices of each forming bond in each transition state
!   forming_bond_lengths - An array listing the lengths of each forming bond in each transition state
!   number_of_forming_bonds - The number of bonds being formed by the reaction
!   breaking_bonds - An array listing the pairs of indices of each breaking bond in each transition state
!   breaking_bond_lengths - An array listing the lengths of each breaking bond in each transition state
!   number_of_breaking_bonds - The number of bonds being broken by the reaction
! Returns:
!   s1 - The value of the dividing surface function
subroutine transition_state_value(position, Natoms, Nts, &
    forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
    breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, s1)

    implicit none
    integer, intent(in) :: Natoms, Nts
    integer, intent(in) :: number_of_forming_bonds, number_of_breaking_bonds
    double precision, intent(in) :: position(3,Natoms)
    integer, intent(in) :: forming_bonds(Nts,number_of_forming_bonds,2)
    integer, intent(in) :: breaking_bonds(Nts,number_of_breaking_bonds,2)
    double precision, intent(in) :: forming_bond_lengths(Nts,number_of_forming_bonds)
    double precision, intent(in) :: breaking_bond_lengths(Nts,number_of_breaking_bonds)
    double precision, intent(out) :: s1

    double precision :: values(Nts)

    call transition_state_evaluate_all(position, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, values)

    s1 = maxval(values)

end subroutine transition_state_value

! Return the value of the gradient of the transition state dividing surface
! function.
! Parameters:
!   position - A 3 x Natoms array of atomic positions
!   Natoms - The number of atoms
!   Nts - The number of equivalent transition states that define the dividing surface
!   forming_bonds - An array listing the pairs of indices of each forming bond in each transition state
!   forming_bond_lengths - An array listing the lengths of each forming bond in each transition state
!   number_of_forming_bonds - The number of bonds being formed by the reaction
!   breaking_bonds - An array listing the pairs of indices of each breaking bond in each transition state
!   breaking_bond_lengths - An array listing the lengths of each breaking bond in each transition state
!   number_of_breaking_bonds - The number of bonds being broken by the reaction
! Returns:
!   ds1 - The gradient of the dividing surface function
subroutine transition_state_gradient(position, Natoms, Nts, &
    forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
    breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, ds1)

    implicit none
    integer, intent(in) :: Natoms, Nts
    integer, intent(in) :: number_of_forming_bonds, number_of_breaking_bonds
    double precision, intent(in) :: position(3,Natoms)
    integer, intent(in) :: forming_bonds(Nts,number_of_forming_bonds,2)
    integer, intent(in) :: breaking_bonds(Nts,number_of_breaking_bonds,2)
    double precision, intent(in) :: forming_bond_lengths(Nts,number_of_forming_bonds)
    double precision, intent(in) :: breaking_bond_lengths(Nts,number_of_breaking_bonds)
    double precision, intent(out) :: ds1(3,Natoms)

    double precision :: values(Nts)
    integer :: m, n, atom1, atom2
    double precision :: Rx, Ry, Rz, Rinv

    ds1(:,:) = 0.0

    call transition_state_evaluate_all(position, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, values)

    n = maxloc(values, 1)

    do m = 1, number_of_forming_bonds
        atom1 = forming_bonds(n,m,1)
        atom2 = forming_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        ds1(1,atom1) = ds1(1,atom1) - Rx * Rinv
        ds1(2,atom1) = ds1(2,atom1) - Ry * Rinv
        ds1(3,atom1) = ds1(3,atom1) - Rz * Rinv
        ds1(1,atom2) = ds1(1,atom2) + Rx * Rinv
        ds1(2,atom2) = ds1(2,atom2) + Ry * Rinv
        ds1(3,atom2) = ds1(3,atom2) + Rz * Rinv
    end do

    do m = 1, number_of_breaking_bonds
        atom1 = breaking_bonds(n,m,1)
        atom2 = breaking_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        ds1(1,atom1) = ds1(1,atom1) + Rx * Rinv
        ds1(2,atom1) = ds1(2,atom1) + Ry * Rinv
        ds1(3,atom1) = ds1(3,atom1) + Rz * Rinv
        ds1(1,atom2) = ds1(1,atom2) - Rx * Rinv
        ds1(2,atom2) = ds1(2,atom2) - Ry * Rinv
        ds1(3,atom2) = ds1(3,atom2) - Rz * Rinv
    end do

end subroutine transition_state_gradient

! Return the value of the Hessian of the transition state dividing surface
! function.
! Parameters:
!   position - A 3 x Natoms array of atomic positions
!   Natoms - The number of atoms
!   Nts - The number of equivalent transition states that define the dividing surface
!   forming_bonds - An array listing the pairs of indices of each forming bond in each transition state
!   forming_bond_lengths - An array listing the lengths of each forming bond in each transition state
!   number_of_forming_bonds - The number of bonds being formed by the reaction
!   breaking_bonds - An array listing the pairs of indices of each breaking bond in each transition state
!   breaking_bond_lengths - An array listing the lengths of each breaking bond in each transition state
!   number_of_breaking_bonds - The number of bonds being broken by the reaction
! Returns:
!   d2s1 - The Hessian of the dividing surface function
subroutine transition_state_hessian(position, Natoms, Nts, &
    forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
    breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, d2s1)

    implicit none
    integer, intent(in) :: Natoms, Nts
    integer, intent(in) :: number_of_forming_bonds, number_of_breaking_bonds
    double precision, intent(in) :: position(3,Natoms)
    integer, intent(in) :: forming_bonds(Nts,number_of_forming_bonds,2)
    integer, intent(in) :: breaking_bonds(Nts,number_of_breaking_bonds,2)
    double precision, intent(in) :: forming_bond_lengths(Nts,number_of_forming_bonds)
    double precision, intent(in) :: breaking_bond_lengths(Nts,number_of_breaking_bonds)
    double precision, intent(out) :: d2s1(3,Natoms,3,Natoms)

    double precision :: values(Nts)
    integer :: m, n, atom1, atom2
    double precision :: Rx, Ry, Rz, Rinv
    double precision :: dxx, dyy, dzz, dxy, dxz, dyz

    d2s1(:,:,:,:) = 0.0

    call transition_state_evaluate_all(position, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, values)

    n = maxloc(values, 1)

    do m = 1, number_of_forming_bonds
        atom1 = forming_bonds(n,m,1)
        atom2 = forming_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        dxx = -(Ry * Ry + Rz * Rz) * (Rinv * Rinv * Rinv)
        dyy = -(Rz * Rz + Rx * Rx) * (Rinv * Rinv * Rinv)
        dzz = -(Rx * Rx + Ry * Ry) * (Rinv * Rinv * Rinv)
        dxy = Rx * Ry * (Rinv * Rinv * Rinv)
        dxz = Rx * Rz * (Rinv * Rinv * Rinv)
        dyz = Ry * Rz * (Rinv * Rinv * Rinv)

        d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx
        d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy
        d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz
        d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy
        d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy
        d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz
        d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz
        d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz
        d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz

        d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx
        d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy
        d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz
        d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy
        d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy
        d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz
        d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz
        d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz
        d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz

        d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx
        d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy
        d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz
        d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy
        d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy
        d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz
        d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz
        d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz
        d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz

        d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx
        d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy
        d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz
        d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy
        d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy
        d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz
        d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz
        d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz
        d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz
    end do

    do m = 1, number_of_breaking_bonds
        atom1 = breaking_bonds(n,m,1)
        atom2 = breaking_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        dxx = (Ry * Ry + Rz * Rz) * (Rinv * Rinv * Rinv)
        dyy = (Rz * Rz + Rx * Rx) * (Rinv * Rinv * Rinv)
        dzz = (Rx * Rx + Ry * Ry) * (Rinv * Rinv * Rinv)
        dxy = -Rx * Ry * (Rinv * Rinv * Rinv)
        dxz = -Rx * Rz * (Rinv * Rinv * Rinv)
        dyz = -Ry * Rz * (Rinv * Rinv * Rinv)

        d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx
        d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy
        d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz
        d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy
        d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy
        d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz
        d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz
        d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz
        d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz

        d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx
        d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy
        d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz
        d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy
        d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy
        d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz
        d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz
        d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz
        d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz

        d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx
        d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy
        d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz
        d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy
        d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy
        d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz
        d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz
        d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz
        d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz

        d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx
        d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy
        d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz
        d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy
        d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy
        d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz
        d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz
        d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz
        d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz
    end do

end subroutine transition_state_hessian
