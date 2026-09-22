!###############################################################################
!#                                                                             #
!# glm_sed_heat.F90                                                           #
!#                                                                             #
!# The dynamic (sed_heat_model = 2) sediment/soil heat-conduction solver.     #
!#                                                                             #
!# This is a self-contained 1-D finite-volume heat-conduction model with no   #
!# dependency on the AED water-quality library: it is a physical-boundary     #
!# process (bed <-> water heat exchange), not biogeochemistry, so it is       #
!# always compiled into GLM core and works whether or not a WQ module         #
!# (aed/api) is active. Ported from libaed-water's aed_util.F90               #
!# (SoilCond/SoilCap/SoilTempFV/InitialTempFV), which is only reachable when  #
!# AED is linked in.                                                          #
!#                                                                             #
!# Developed by :                                                             #
!#     AquaticEcoDynamics (AED) Group                                         #
!#     School of Agriculture and Environment                                  #
!#     The University of Western Australia                                    #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!# Copyright 2013-2026 : The University of Western Australia                  #
!#                                                                             #
!#  This file is part of GLM (General Lake Model)                             #
!#                                                                             #
!#  GLM is free software: you can redistribute it and/or modify               #
!#  it under the terms of the GNU General Public License as published by      #
!#  the Free Software Foundation, either version 3 of the License, or         #
!#  (at your option) any later version.                                       #
!#                                                                             #
!#  GLM is distributed in the hope that it will be useful,                    #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#  GNU General Public License for more details.                              #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License         #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!###############################################################################

#undef MISVAL

#include "glm.h"

MODULE glm_sed_heat

   USE ISO_C_BINDING
   USE glm_types

   IMPLICIT NONE

   PRIVATE

   PUBLIC soil_props_t, SoilTempFV, InitialTempFV

   AED_REAL,PARAMETER :: zero_ = 0.0, one_ = 1.0

   !# Soil thermal properties for the cell-centred soil-temperature solver
   !# (mirrors the thermal fields of intertidal_soil.SoilParams).
   TYPE :: soil_props_t
      AED_REAL :: k_mineral, k_water
      AED_REAL :: c_mineral, c_water, c_air
      AED_REAL :: bulk_density, mineral_density, porosity
   END TYPE soil_props_t

CONTAINS

!===============================================================================
! Faithful port of intertidal-soil-py (intertidal_soil/temperature.py
! soil_temp + initial_temp; conductivity/heat-capacity from soil_params.py).
! Cell-centred finite volume, fully-implicit (backward Euler), Johansen
! property-based thermal conductivity, harmonic-mean interface conductivities,
! Dirichlet top (surfTemp) and bottom (deepTemp) boundary conditions.
!===============================================================================
PURE FUNCTION SoilCond(theta, p) RESULT(lam)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: theta
   TYPE(soil_props_t),INTENT(in) :: p
!LOCALS
   AED_REAL :: lam, Sr, Ke, k_dry, k_sat
!
!-------------------------------------------------------------------------------
!BEGIN
   Sr = MIN(MAX(theta/p%porosity, zero_), one_)
   Ke = MIN(MAX(LOG10(MAX(Sr, 1e-10)) + one_, zero_), one_)
   k_dry = (0.135*p%bulk_density*1e3 + 64.7) &
         / (p%mineral_density*1e3 - 0.947*p%bulk_density*1e3)
   k_sat = p%k_mineral**(one_ - p%porosity) * p%k_water**p%porosity
   lam = k_dry + (k_sat - k_dry)*Ke
END FUNCTION SoilCond
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
PURE FUNCTION SoilCap(theta, p) RESULT(cap)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,INTENT(in) :: theta
   TYPE(soil_props_t),INTENT(in) :: p
!LOCALS
   AED_REAL :: cap
!
!-------------------------------------------------------------------------------
!BEGIN
   cap = (one_ - p%porosity)*p%c_mineral + theta*p%c_water + (p%porosity - theta)*p%c_air
END FUNCTION SoilCap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
SUBROUTINE SoilTempFV(n, d, vwc, surfTemp, deepTemp, dt, props, temp, heatflux)
!-------------------------------------------------------------------------------
!  n         : number of interior cells (= n_sed_layers - 2)
!  d(0:n+1)  : node depths (m); d(0)=0 (surface), d(n+1)=column depth
!  vwc(1:n)  : volumetric water content per cell
!  surfTemp  : Dirichlet surface temperature (degC)
!  deepTemp  : Dirichlet deep-boundary temperature (degC)
!  dt        : time-step (s)
!  temp(1:n) : (inout) cell-centred temperature profile (degC)
!  heatflux  : (out) surface heat flux W/m2 (positive = into soil)
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: n
   AED_REAL,INTENT(in) :: d(0:n+1), vwc(n), surfTemp, deepTemp, dt
   TYPE(soil_props_t),INTENT(in) :: props
   AED_REAL,INTENT(inout) :: temp(n)
   AED_REAL,INTENT(out)   :: heatflux
!LOCALS
   INTEGER  :: i
   AED_REAL :: dz(n), zc(n), lam(n), cap(n)
   AED_REAL :: dzc(0:n), ki(0:n)
   AED_REAL :: a(n), b(n), c(n), rhs(n), cp(n), dp(n), piv, fa, fb, capv
!
!-------------------------------------------------------------------------------
!BEGIN
   zc(1) = zero_ ; lam(1) = zero_
   b(1) = zero_ ; c(1) = zero_ ; rhs(1) = zero_

   DO i = 1, n
      dz(i)  = d(i+1) - d(i)
      zc(i)  = 0.5*(d(i) + d(i+1))
      lam(i) = SoilCond(vwc(i), props)
      cap(i) = SoilCap(vwc(i), props)
   ENDDO
   !# interface conductivities (0 = surface->cell1 .. n = cellN->deep)
   dzc(0) = zc(1)            ;  ki(0) = lam(1)
   DO i = 1, n-1
      dzc(i) = zc(i+1) - zc(i)
      ki(i)  = 2.0*lam(i+1)*lam(i) / (lam(i+1) + lam(i) + 1e-30)
   ENDDO
   dzc(n) = d(n+1) - zc(n)   ;  ki(n) = lam(n)
   !# assemble tridiagonal system (Dirichlet top & bottom)
   DO i = 1, n
      capv = cap(i)*dz(i)/dt
      fb   = ki(i)/dzc(i)         !# conductance to cell below / deep boundary
      fa   = ki(i-1)/dzc(i-1)     !# conductance from cell above / surface boundary
      b(i)   = capv + fa + fb
      rhs(i) = capv*temp(i)
      IF (i > 1) a(i) = -fa
      IF (i < n) c(i) = -fb
      IF (i == 1) rhs(i) = rhs(i) + fa*surfTemp
      IF (i == n) rhs(i) = rhs(i) + fb*deepTemp
   ENDDO
   !# Thomas tridiagonal solve
   cp(1) = c(1)/b(1)  ;  dp(1) = rhs(1)/b(1)
   DO i = 2, n
      piv = b(i) - a(i)*cp(i-1)
      IF (i < n) cp(i) = c(i)/piv
      dp(i) = (rhs(i) - a(i)*dp(i-1))/piv
   ENDDO
   temp(n) = dp(n)
   DO i = n-1, 1, -1
      temp(i) = dp(i) - cp(i)*temp(i+1)
   ENDDO
   heatflux = ki(0)/dzc(0) * (surfTemp - temp(1))
END SUBROUTINE SoilTempFV
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!===============================================================================
SUBROUTINE InitialTempFV(n, d, vwc_init, surfTemp, deepTemp, spinup_days, dt, props, temp)
!-------------------------------------------------------------------------------
! Spin up the cell-centred profile to quasi-equilibrium (port of initial_temp).
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)  :: n
   AED_REAL,INTENT(in) :: d(0:n+1), vwc_init, surfTemp, deepTemp, spinup_days, dt
   TYPE(soil_props_t),INTENT(in) :: props
   AED_REAL,INTENT(out) :: temp(n)
!LOCALS
   INTEGER  :: i, s, nsteps
   AED_REAL :: vwc(n), hf
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i = 1, n
      temp(i) = surfTemp + (deepTemp - surfTemp)*REAL(i-1)/REAL(MAX(n-1,1))
      vwc(i)  = vwc_init
   ENDDO
   nsteps = INT(spinup_days*86400.0/dt)
   DO s = 1, nsteps
      CALL SoilTempFV(n, d, vwc, surfTemp, deepTemp, dt, props, temp, hf)
   ENDDO
END SUBROUTINE InitialTempFV
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTempFV(n,depth,vwc_init,surfTemp,deepTemp,spinup_days,dt,tNew) &
                                                  BIND(C, name="InitialTempFV")
!-------------------------------------------------------------------------------
!# C-callable spin-up for the cell-centred soil-temperature model.  Reads the
!# soil thermal properties from the (Fortran-owned) glm_types globals and
!# returns the n interior-cell temperatures in tNew(1:n).
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,intent(in)   :: n
   AED_REAL,intent(in)  :: depth(0:n+1), vwc_init, surfTemp, deepTemp, spinup_days, dt
   AED_REAL,intent(out) :: tNew(n)
!LOCALS
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL InitialTempFV(n, depth, vwc_init, surfTemp, deepTemp, spinup_days, dt, props, tNew)
END SUBROUTINE GInitialTempFV
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ZSoilTemp(izone) BIND(C, name="zZSoilTemp")
!-------------------------------------------------------------------------------
!ARGUMENTS
   !# izone is a C pointer passed BY VALUE from glm_surface.c (ZSoilTemp(&theZones[z]));
   !# it MUST carry the VALUE attribute or Fortran treats it as pass-by-reference and
   !# dereferences the struct's first word as the address.
   TYPE(C_PTR),VALUE :: izone
!LOCALS
   TYPE(ZoneType),POINTER :: zone
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
   INTEGER :: nsl, n, kk
   AED_REAL :: hf, deepTemp
   AED_REAL,DIMENSION(:),ALLOCATABLE :: sdepth, cellvwc, cellt
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(izone, zone);
   nsl = zone%n_sed_layers
   n   = nsl - 2
   CALL C_F_POINTER(zone%c_layers, layers, [nsl]);
   !# Cell-centred layout: node depths d(0:n+1) = layers(1:nsl)%depth; the n
   !# interior cells (the prognostic state) live in slots 2..nsl-1; layers(1) is
   !# the surface boundary node and layers(nsl) the deep boundary node.
   ALLOCATE(sdepth(0:n+1), cellvwc(n), cellt(n))
   DO kk = 1, nsl
      sdepth(kk-1) = layers(kk)%depth
   ENDDO
   DO kk = 1, n
      cellvwc(kk) = layers(kk+1)%vwc
      cellt(kk)   = layers(kk+1)%temp
   ENDDO
   deepTemp = layers(nsl)%temp
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL SoilTempFV(n, sdepth, cellvwc, zone%ztemp, deepTemp, soil_dt, props, cellt, hf)
   zone%heatflux = hf
   DO kk = 1, n
      layers(kk+1)%temp = cellt(kk)
   ENDDO
   layers(1)%temp = zone%ztemp
   DEALLOCATE(sdepth, cellvwc, cellt)
END SUBROUTINE ZSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_sed_heat
