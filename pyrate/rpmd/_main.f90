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

! Compute the value, gradient, and Hessian of the reaction coordinate.
! Parameters:
!   centroid - The centroid of each atom
!   Natoms - The number of atoms in the molecular system
!   Rinf - The distance at which the reactant interaction becomes negligible
!   massfrac - The mass fraction of each atom
!   reactant1_atoms - An array of indices for each atom in the first reactant
!   Nreactant1_atoms - The number of atoms in the first reactant
!   reactant2_atoms - An array of indices for each atom in the second reactant
!   Nreactant2_atoms - The number of atoms in the second reactant
!   Nts - The number of equivalent transition states that define the dividing surface
!   forming_bonds - An array listing the pairs of indices of each forming bond in each transition state
!   forming_bond_lengths - An array listing the lengths of each forming bond in each transition state
!   number_of_forming_bonds - The number of bonds being formed by the reaction
!   breaking_bonds - An array listing the pairs of indices of each breaking bond in each transition state
!   breaking_bond_lengths - An array listing the lengths of each breaking bond in each transition state
!   number_of_breaking_bonds - The number of bonds being broken by the reaction
!   xi_current - The maximum of the reaction coordinate at the current temperature
!   mode - 1 for umbrella integration, 2 for recrossing factor
! Returns:
!   xi - The value of the reaction coordinate
!   dxi - The gradient of the reaction coordinate
!   d2xi - The Hessian of the reaction coordinate
subroutine get_reaction_coordinate(centroid, xi, dxi, d2xi, Natoms, &
  Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
  Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
  breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
  xi_current, mode)

    integer, intent(in) :: Natoms
    double precision, intent(in) :: centroid(3,Natoms)
    integer, intent(in) :: mode
    double precision, intent(in) :: xi_current
    ! Reactant dividing surface
    integer, intent(in) :: Nreactant1_atoms, Nreactant2_atoms
    double precision, intent(in) :: massfrac(Natoms)
    integer, intent(in) :: reactant1_atoms(Nreactant1_atoms), reactant2_atoms(Nreactant2_atoms)
    double precision, intent(in) :: Rinf
    ! Transition state dividing surface
    integer, intent(in) :: Nts
    integer, intent(in) :: number_of_forming_bonds, number_of_breaking_bonds
    integer, intent(in) :: forming_bonds(Nts,number_of_forming_bonds,2)
    integer, intent(in) :: breaking_bonds(Nts,number_of_breaking_bonds,2)
    double precision, intent(in) :: forming_bond_lengths(Nts,number_of_forming_bonds)
    double precision, intent(in) :: breaking_bond_lengths(Nts,number_of_breaking_bonds)
    ! Result
    double precision, intent(out) :: xi, dxi(3,Natoms), d2xi(3,Natoms,3,Natoms)

    double precision :: s0, ds0(3,Natoms), d2s0(3,Natoms,3,Natoms)
    double precision :: s1, ds1(3,Natoms), d2s1(3,Natoms,3,Natoms)
    integer :: i1, i2, j1, j2

    xi = 0.0d0
    dxi(:,:) = 0.0d0
    d2xi(:,:,:,:) = 0.0d0

    ! Evaluate reactants dividing surface
    call reactants_value(centroid, Natoms, Rinf, massfrac, &
        reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
        s0)
    ! Evaluate reactants dividing surface gradient
    call reactants_gradient(centroid, Natoms, massfrac, &
        reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
        ds0)
    ! Evaluate reactants dividing surface Hessian
    call reactants_hessian(centroid, Natoms, massfrac, &
        reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
        d2s0)

    ! Evaluate transition state dividing surface
    call transition_state_value(centroid, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, s1)
    ! Evaluate transition state dividing surface gradient
    call transition_state_gradient(centroid, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, ds1)
    ! Evaluate transition state dividing surface Hessian
    call transition_state_hessian(centroid, Natoms, Nts, &
        forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, d2s1)

    ! Compute reaction coordinate and its gradient and Hessian
    ! The functional form is different depending on whether we are doing
    ! umbrella integration or a recrossing factor calculation
    if (mode .eq. 1) then
        ! Umbrella integration
        xi = s0 / (s0 - s1)
        dxi = (s0 * ds1 - s1 * ds0) / ((s0 - s1) * (s0 - s1))
        do i1 = 1, 3
            do j1 = 1, Natoms
                do i2 = 1, 3
                    do j2 = 1, Natoms
                        d2xi(i1,j1,i2,j2) = ((s0 * d2s1(i1,j1,i2,j2) + ds0(i2,j2) * ds1(i1,j1) &
                            - ds1(i2,j2) * ds0(i1,j1) - s1 * d2s0(i1,j1,i2,j2)) * (s0 - s1) &
                            - 2.0d0 * (s0 * ds1(i1,j1) - s1 * ds0(i1,j1)) &
                            * (ds0(i2,j2) - ds1(i2,j2))) &
                            / ((s0 - s1) * (s0 - s1) * (s0 - s1))
                    end do
                end do
            end do
        end do
    elseif (mode .eq. 2) then
        ! Recrossing factor
        xi = xi_current * s1 + (1 - xi_current) * s0
        dxi = xi_current * ds1 + (1 - xi_current) * ds0
        ! Don't need Hessian for recrossing factor, so don't compute it
    end if

end subroutine get_reaction_coordinate

! Compute the total energy of all ring polymers in the RPMD system.
! Parameters:
!   q - The position of each bead in each atom
!   mass - The mass of each atom
!   beta - The inverse temperature of interest
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   Ering - The total energy of all ring polymers
subroutine get_ring_polymer_energy(q, mass, beta, Natoms, Nbeads, Ering)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: q(3,Natoms,Nbeads), mass(Natoms)
    double precision, intent(in) :: beta
    double precision, intent(out) :: Ering

    double precision :: wn, dx, dy, dz
    integer :: j, k

    Ering = 0.0d0
    wn = Nbeads / beta
    do j = 1, Natoms
        dx = q(1,j,1) - q(1,j,Nbeads)
        dy = q(2,j,1) - q(2,j,Nbeads)
        dz = q(3,j,1) - q(3,j,Nbeads)
        Ering = Ering + 0.5d0 * mass(j) * wn * wn * (dx * dx + dy * dy + dz * dz)
        do k = 2, Nbeads
            dx = q(1,j,k-1) - q(1,j,k)
            dy = q(2,j,k-1) - q(2,j,k)
            dz = q(3,j,k-1) - q(3,j,k)
            Ering = Ering + 0.5d0 * mass(j) * wn * wn * (dx * dx + dy * dy + dz * dz)
        end do
    end do

end subroutine get_ring_polymer_energy

! Compute the kinetic energy of the molecular system.
! Parameters:
!   p - The momentum of each bead in each atom
!   mass - The mass of each atom
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   Ek - The kinetic energy of the system
subroutine get_kinetic_energy(p, mass, Natoms, Nbeads, Ek)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: p(3,Natoms,Nbeads), mass(Natoms)
    double precision, intent(out) :: Ek

    integer :: i, j, k

    Ek = 0.0d0
    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                Ek = Ek + 0.5d0 * p(i,j,k) * p(i,j,k) / mass(j)
            end do
        end do
    end do

end subroutine get_kinetic_energy

! Compute the center of mass position of the molecular system.
! Parameters:
!   q - The position of each bead in each atom
!   mass - The mass of each atom
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   cm - The center of mass of the system
subroutine get_center_of_mass(q, mass, Natoms, Nbeads, cm)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: q(3,Natoms,Nbeads), mass(Natoms)
    double precision, intent(out) :: cm(3)

    double precision :: total_mass
    integer :: i, j, k

    cm(:) = 0.0d0
    total_mass = sum(mass)
    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                cm(i) = cm(i) + q(i,j,k) * mass(j)
            end do
        end do
        cm(i) = cm(i) / total_mass
    end do

end subroutine get_center_of_mass

! Compute the centroid position of each atom in the molecular system.
! Parameters:
!   q - The position of each bead in each atom
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   centroid - The centroid of each atom
subroutine get_centroid(q, Natoms, Nbeads, centroid)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: q(3,Natoms,Nbeads)
    double precision, intent(out) :: centroid(3,Natoms)

    integer :: i, j, k

    centroid(:,:) = 0.0d0
    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                centroid(i,j) = centroid(i,j) + q(i,j,k)
            end do
            centroid(i,j) = centroid(i,j) / Nbeads
        end do
    end do

end subroutine get_centroid

! Compute the radius of gyration of each atom in the molecular system. This is
! a useful quantity to check while debugging.
! Parameters:
!   q - The position of each bead in each atom
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   R - The radius of gyration of each atom
subroutine get_radius_of_gyration(q, Natoms, Nbeads, R)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: q(3,Natoms,Nbeads)
    double precision, intent(out) :: R(Natoms)

    double precision :: centroid(3,Natoms), dx
    integer :: i, j, k

    call get_centroid(q, Natoms, Nbeads, centroid)

    R(:) = 0.0d0
    do j = 1, Natoms
        do i = 1, 3
            do k = 1, Nbeads
                dx = q(i,j,k) - centroid(i,j)
                R(j) = R(j) + dx * dx
            end do
        end do
        R(j) = sqrt(R(j) / Nbeads)
    end do

end subroutine get_radius_of_gyration

! Compute a pseudo-random number uniformly distributed in [0,1].
! Returns:
!   rn - The pseudo-random number
subroutine random(rn)
    implicit none
    double precision, intent(out) :: rn
    call random_number(rn)
end subroutine random

! Compute a pseudo-random number weighted by a standard normal distribution.
! The Marsaglia polar method is used to convert from a uniform distribution to
! a normal distribution.
! Returns:
!   rn - The pseudo-random number
subroutine randomn(rn)

    implicit none
    double precision, intent(out) :: rn
    integer :: iset
    double precision :: gset, u, v, S, fac

    save iset,gset

    data iset/0/

    ! The Marsaglia polar method
    if (iset .eq. 0) then
        S = 1.d0
        do while (S .ge. 1.d0 .or. S .eq. 0.d0)
            call random_number(u)
            call random_number(v)
            u = 2.d0 * u - 1.d0
            v = 2.d0 * v - 1.d0
            S = u * u + v * v
        end do
        fac = sqrt(-2 * log(S) / S)
        gset = u * fac
        rn = v * fac
        iset = 1
    else
        rn = gset
        iset = 0
    end if

end subroutine randomn

! Compute the real fast Fourier transform of the given array of data.
! Parameters:
!   x - The array of data to transform, in half-complex form
!   N - The length of the array of data
! Returns:
!   x - The transformed array of data
subroutine rfft(x,N)

    implicit none
    integer, intent(in) :: N
    double precision, intent(inout) :: x(N)
    integer*8 plan

    integer Np, Nmax
    parameter (Nmax = 512)
    double precision copy(Nmax), factor

    data Np /0/
    save copy, factor, plan, Np

    if (N .ne. Np) then
        call dfftw_plan_r2r_1d(plan,N,copy,copy,0,64)
        factor = dsqrt(1.d0/N)
        Np = N
    end if

    copy(1:N) = x
    call dfftw_execute(plan)
    x = factor * copy(1:N)

end subroutine rfft

! Compute the inverse real fast Fourier transform of the given array of data.
! Parameters:
!   x - The array of data to transform, in half-complex form
!   N - The length of the array of data
! Returns:
!   x - The transformed array of data
subroutine irfft(x,N)

    implicit none
    integer, intent(in) :: N
    double precision, intent(inout) :: x(N)
    integer*8 plan

    integer Np, Nmax
    parameter (Nmax = 512)
    double precision copy(Nmax), factor

    data Np /0/
    save copy, factor, plan, Np

    if (N .ne. Np) then
        call dfftw_plan_r2r_1d(plan,N,copy,copy,1,64)
        factor = dsqrt(1.d0/N)
        Np = N
    end if

    copy(1:N) = x
    call dfftw_execute(plan)
    x = factor * copy(1:N)

end subroutine irfft
