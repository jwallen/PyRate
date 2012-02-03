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

! Evolve the RPMD system for one time step, optionally applying a constraint
! that the system remain on the transition state dividing surface.
! Parameters:
!   p - The momentum of each bead in each atom
!   q - The position of each bead in each atom
!   V - The potential of each bead
!   dVdq - The force exerted on each bead in each atom
!   xi - The value of the reaction coordinate
!   dxi - The gradient of the reaction coordinate
!   d2xi - The Hessian of the reaction coordinate
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
!   constrain - 1 to constrain the system to the dividing surface, 0 otherwise
!   beta - The inverse temperature at which to compute the transmission coefficient
!   dt - The time step to use in the simulation
!   mass - The mass of each atom in the molecular system
!   potential - A function that evaluates the potential and force for a given position
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
!   p - The updated momentum of each bead in each atom
!   q - The updated position of each bead in each atom
!   V - The updated potential of each bead
!   dVdq - The updated force exerted on each bead in each atom
!   xi - The updated value of the reaction coordinate
!   dxi - The updated gradient of the reaction coordinate
!   d2xi - The updated Hessian of the reaction coordinate
subroutine evolve(p, q, V, dVdq, xi, dxi, d2xi, Natoms, Nbeads, constrain, &
  beta, dt, mass, potential, &
  Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
  Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
  breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
  xi_current, mode)

    implicit none
    external potential
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(inout) :: p(3,Natoms,Nbeads), q(3,Natoms,Nbeads)
    double precision, intent(inout) :: V(Nbeads), dVdq(3,Natoms,Nbeads)
    double precision, intent(inout) :: xi, dxi(3,Natoms), d2xi(3,Natoms,3,Natoms)
    integer, intent(in) :: constrain
    double precision, intent(in) :: beta, dt, xi_current, mass(Natoms)
    integer, intent(in) :: mode
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

    double precision :: centroid(3,Natoms)
    integer :: i, j

    ! Update momentum (half time step)
    p = p - 0.5d0 * dt * dVdq

    ! Update position (full time step)
    if (Nbeads .eq. 1) then
        do i = 1, 3
            do j = 1, Natoms
                q(i,j,1) = q(i,j,1) + p(i,j,1) * dt / mass(j)
            end do
        end do
    else
        call evolve_free_ring_polymer(p, q, beta, dt, mass, Natoms, Nbeads)
    end if

    ! If constrain is on, the evolution will be constrained to the
    ! transition state dividing surface
    if (constrain .eq. 1) call constrain_to_dividing_surface(p, q, dxi, mass, &
        dt, Natoms, Nbeads, &
        Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
        Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
        xi_current, mode)

    ! Update reaction coordinate and its gradient and Hessian
    call get_centroid(q, Natoms, Nbeads, centroid)
    call get_reaction_coordinate(centroid, xi, dxi, d2xi, Natoms, &
        Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
        Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
        breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
        xi_current, mode)

    ! Update potential using new position
    call get_potential(q, xi, dxi, d2xi, V, dVdq, Natoms, Nbeads, potential)

    ! Update momentum (half time step) using new potential
    p = p - 0.5d0 * dt * dVdq

    ! Constrain momentum again
    if (constrain .eq. 1) call constrain_momentum_to_dividing_surface(p, mass, dxi, Natoms, Nbeads)

end subroutine evolve

! Evolve the position and momentum in time according to the free ring polymer
! term in the Hamiltonian. This is most efficiently done in normal mode space
! instead of Cartesian space; fast Fourier transforms are employed to move to
! and from normal mode space.
! Parameters:
!   p - The momentum of each bead in each atom
!   q - The position of each bead in each atom
!   beta - The inverse temperature at which to compute the transmission coefficient
!   dt - The time step to use in the simulation
!   mass - The mass of each atom in the molecular system
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   p - The updated momentum of each bead in each atom
!   q - The updated position of each bead in each atom
subroutine evolve_free_ring_polymer(p, q, beta, dt, mass, Natoms, Nbeads)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(inout) :: p(3,Natoms,Nbeads), q(3,Natoms,Nbeads)
    double precision, intent(in) :: beta, dt, mass(Natoms)

    double precision :: poly(4,Nbeads)
    double precision :: pi, beta_n, twown, pi_n, wk, wt, wm, cos_wt, sin_wt, p_new
    integer :: i, j, k

    pi = dacos(-1.0d0)

    ! Transform to normal mode space
    do i = 1, 3
        do j = 1, Natoms
            call rfft(p(i,j,:), Nbeads)
            call rfft(q(i,j,:), Nbeads)
        end do
    end do

    do j = 1, Natoms

        poly(1,1) = 1.0d0
        poly(2,1) = 0.0d0
        poly(3,1) = dt / mass(j)
        poly(4,1) = 1.0d0

        if (Nbeads .gt. 1) then
            beta_n = beta / Nbeads
            twown = 2.0d0 / beta_n
            pi_n = pi / Nbeads
            do k = 1, Nbeads / 2
                wk = twown * dsin(k * pi_n)
                wt = wk * dt
                wm = wk * mass(j)
                cos_wt = dcos(wt)
                sin_wt = dsin(wt)
                poly(1,k+1) = cos_wt
                poly(2,k+1) = -wm*sin_wt
                poly(3,k+1) = sin_wt/wm
                poly(4,k+1) = cos_wt
            end do
            do k = 1, (Nbeads - 1) / 2
                poly(1,Nbeads-k+1) = poly(1,k+1)
                poly(2,Nbeads-k+1) = poly(2,k+1)
                poly(3,Nbeads-k+1) = poly(3,k+1)
                poly(4,Nbeads-k+1) = poly(4,k+1)
            end do
        end if

        do k = 1, Nbeads
            do i = 1, 3
                p_new = p(i,j,k) * poly(1,k) + q(i,j,k) * poly(2,k)
                q(i,j,k) = p(i,j,k) * poly(3,k) + q(i,j,k) * poly(4,k)
                p(i,j,k) = p_new
            end do
        end do

    end do

    ! Transform back to Cartesian space
    do i = 1, 3
        do j = 1, Natoms
            call irfft(p(i,j,:), Nbeads)
            call irfft(q(i,j,:), Nbeads)
        end do
    end do

end subroutine evolve_free_ring_polymer

! Constrain the position and the momentum to the dividing surface, using the
! SHAKE/RATTLE algorithm.
! Parameters:
!   p - The momentum of each bead in each atom
!   q - The position of each bead in each atom
!   dxi - The gradient of the reaction coordinate
!   mass - The mass of each atom
!   dt - The time step to use in the simulation
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
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
!   p - The constrained momentum of each bead in each atom
!   q - The constrained position of each bead in each atom
subroutine constrain_to_dividing_surface(p, q, dxi, mass, &
  dt, Natoms, Nbeads, &
  Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
  Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
  breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
  xi_current, mode)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(inout) :: p(3,Natoms,Nbeads), q(3,Natoms,Nbeads)
    integer, intent(in) :: mode
    double precision, intent(in) :: xi_current, dt, mass(Natoms)
    double precision, intent(in) :: dxi(3,Natoms)
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

    double precision :: centroid(3,Natoms), qctemp(3,Natoms)
    integer :: i, j, k, maxiter, iter
    double precision :: xi_new, dxi_new(3,Natoms), d2xi_new(3,Natoms,3,Natoms)
    double precision :: mult, sigma, dsigma, dx, coeff

    call get_centroid(q, Natoms, Nbeads, centroid)

    ! The Lagrange multiplier for the constraint
    mult = 0.0d0

    qctemp(:,:) = 0.0d0

    maxiter = 100
    do iter = 1, maxiter

        coeff = mult * dt * dt / Nbeads

        do i = 1, 3
            do j = 1, Natoms
                qctemp(i,j) = centroid(i,j) + coeff * dxi(i,j) / mass(j)
            end do
        end do

        call get_reaction_coordinate(qctemp, xi_new, dxi_new, d2xi_new, Natoms, &
            Rinf, massfrac, reactant1_atoms, Nreactant1_atoms, reactant2_atoms, Nreactant2_atoms, &
            Nts, forming_bonds, forming_bond_lengths, number_of_forming_bonds, &
            breaking_bonds, breaking_bond_lengths, number_of_breaking_bonds, &
            xi_current, mode)

        sigma = xi_new
        dsigma = 0.0d0
        do i = 1, 3
            do j = 1, Natoms
                dsigma = dsigma + dxi_new(i,j) * dt * dt * dxi(i,j) / (mass(j) * Nbeads)
            end do
        end do

        dx = sigma / dsigma
        mult = mult - dx
        if (dabs(dx) .lt. 1.0d-9 .or. dabs(sigma) .lt. 1.0d-13) exit

        if (iter .eq. maxiter) write (*,fmt='(A)') 'SHAKE exceeded maximum number of iterations.'

    end do

    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                q(i,j,k) = q(i,j,k) + coeff / mass(j) * dxi(i,j)
                p(i,j,k) = p(i,j,k) + mult * dt / Nbeads * dxi(i,j)
            end do
        end do
    end do

end subroutine constrain_to_dividing_surface

! Constrain the momentum to the reaction coordinate, to ensure that the time
! derivative of the dividing surface is zero.
! Parameters:
!   p - The momentum of each bead in each atom
!   mass - The mass of each atom
!   dxi - The gradient of the reaction coordinate
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
! Returns:
!   p - The constrained momentum of each bead in each atom
subroutine constrain_momentum_to_dividing_surface(p, mass, dxi, Natoms, Nbeads)

    implicit none
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: dxi(3,Natoms), mass(Natoms)
    double precision, intent(inout) :: p(3,Natoms,Nbeads)

    double precision :: coeff1, coeff2, lambda
    integer :: i, j, k

    coeff1 = 0.0d0
    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                coeff1 = coeff1 + dxi(i,j) * p(i,j,k) / mass(j)
            end do
        end do
    end do

    coeff2 = 0.0d0
    do i = 1, 3
        do j = 1, Natoms
            coeff2 = coeff2 + dxi(i,j) * dxi(i,j) / mass(j)
        end do
    end do

    lambda = -coeff1 / coeff2 / Nbeads
    do i = 1, 3
        do j = 1, Natoms
            do k = 1, Nbeads
                p(i,j,k) = p(i,j,k) + lambda * dxi(i,j)
            end do
        end do
    end do

    ! DEBUG: Check that constraint is correct: coeff1 should now evaluate to
    ! zero within numerical precision
    !coeff1 = 0.0d0
    !do i = 1, 3
    !    do j = 1, Natoms
    !        do k = 1, Nbeads
    !            coeff1 = coeff1 + dxi(i,j) * p(i,j,k) / mass(j)
    !        end do
    !    end do
    !end do

end subroutine constrain_momentum_to_dividing_surface

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

! Compute the potential for the given position.
! Parameters:
!   q - The momentum of each bead in each atom
!   xi - The value of the reaction coordinate
!   dxi - The gradient of the reaction coordinate
!   d2xi - The Hessian of the reaction coordinate
!   Natoms - The number of atoms in the molecular system
!   Nbeads - The number of beads to use per atom
!   potential - A function that evaluates the potential and force for a given position
! Returns:
!   V - The potential of each bead
!   dVdq - The force exerted on each bead in each atom
subroutine get_potential(q, xi, dxi, d2xi, V, dVdq, Natoms, Nbeads, potential)

    implicit none
    external potential
    integer, intent(in) :: Natoms, Nbeads
    double precision, intent(in) :: q(3,Natoms,Nbeads)
    double precision, intent(in) :: xi, dxi(3,Natoms), d2xi(3,Natoms,3,Natoms)
    double precision, intent(out) :: V(Nbeads), dVdq(3,Natoms,Nbeads)

    call potential(q, V, dVdq, Natoms, Nbeads)

end subroutine get_potential

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
