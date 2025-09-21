#
# EZmock.pyx: this file is part of the EZmock python package.
#
# EZmock: Effective Zel'dovich approximation mock generator.
#
# Github repository:
#       https://github.com/cheng-zhao/EZmock
#
# Copyright (c) 2025 Yunyi Tang <yunyi-tang@outlook.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import numpy as np
from numpy.typing import NDArray
from typing import Optional, Literal, Tuple

class EZmock:
    def __init__(
        self,
        Lbox: float,
        Ngrid: float,
        rng: Literal[0, 1] = 1,
        seed: int = 1,
        nthreads: int = 0
    ) -> None:
        '''
        Initialize the EZmock instance.

        Parameters
        ----------
        Lbox: float
            Side length of the cubic simulation box.
        Ngrid: int
            Number of grid cells per box side for the density field.
        rng: Literal[0, 1], optional
            Random number generation algorithm (0: MRG32K3A, 1: MT19937).
            Default: 1.
        seed: int, optional
            Random seed.
            Default: 1.
        nthread: int, optional
            Number of OpenMP threads to be used.
            Default: 0 (maximum number of threads available).
        '''
        ...

    def set_growth_params(
        self,
        growth_pk: float,
        vel_fac: float,
    ) -> None:
        '''
        Set the structure growth parameters directly.

        Parameters
        ----------
        growth_pk: float
            Factor that normalizes the input linear power spectrum at the
            desired redshift of the output catalog: ( D(z_out) / D(z_pk) )**2,
            where D indicates the linear growth factor.
        vel_fac: float
            Factor for converting Lagrangian displacements into physical peculiar
            velocities: f(z_out) * H(z_out) * a(z_out) / h, where f is the linear
            growth rate, H the Hubble parameter, a the scale factor, and h the
            dimensionless Hubble parameter.
        '''
        ...

    def eval_growth_params(
        self,
        z_out: float,
        Omega_m: float,
        z_pk: float = 0,
        Omega_nu: float = 0,
        w: float = -1
    ) -> None:
        '''
        Evaluate the structure growth parameters with cosmological parameters.

        Parameters
        ----------
        z_out: float
            Redshift of the output catalog.
        Omega_m: float
            Matter (without neutrino) density parameter at z = 0.
        z_pk: float, optional
            Redshift of the input linear power spectrum.
            Default: 0.
        Omega_nu: float, optional
            Neutrino density parameter at z = 0.
            Default: 0.
        w: float, optional
            Dark energy equation of state.
            Default: -1.
        '''
        ...

    def setup_linear_pk(
        self,
        k: NDArray,
        Plin: NDArray,
        Pnw: Optional[NDArray] = None,
        BAO_enhance: float = 0,
        interp_log: bool = False
    ) -> None:
        '''
        Setup the linear power spectrum for generating the initial condition.

        Parameters
        ----------
        k: float numpy array
            An ascending array of wavenumbers for the input power spectra.
        Plin: float numpy array
            The linear power spectrum.
        Pnw: float numpy array, optional
            The linear power spectrum without BAO wiggles (non-wiggle).
            It is only necessary if `BAO_enhance` != 0.
            Default: None.
        BAO_enhance: float, optional
            The BAO enhancement parameter (positive: enhance BAO;
            negative: damp BAO; 0: no effect).
            Default: 0.
        interp_log: boolean, optional
            True for interpolating the power spectra in log scale; False for
            interpolating in linear scale.
            Default: False.

        Prerequisite
        ------------
        Function `set_growth_params` or `eval_growth_params`.
        '''
        ...


    def create_dens_field_from_ic(
        self,
        fixamp: bool = False,
        iphase: bool = False
    ) -> None:
        '''
        Generate the EZmock density field from a fresh initial condition based on
        the input linear power spectra.

        Parameters
        ----------
        fixamp: boolean, optional
            True for fixing the amplitude of the initial condition.
            Default: False.
        iphase: boolean, optional
            True for inverting the phases of the initial condition.
            Default: False.

        Prerequisite
        ------------
        Function `setup_linear_pk`.
        '''
        ...

    def create_dens_field_from_wn(
        self,
        delta: NDArray,
    ) -> None:
        '''
        Generate the EZmock density field given a white noise field and the
        input linear power spectra.

        Parameters
        ----------
        delta: float numpy array
            The configuration-space white noise field, with dimension `Ngrid`**3,
            and C-order indices.

        Prerequisite
        ------------
        Function `setup_linear_pk`.
        '''
        ...

    def create_dens_field_from_disp(
        self,
        dx: NDArray,
        dy: NDArray,
        dz: NDArray,
        deepcopy: bool = False
        ) -> None:
        '''
        Generate the EZmock density field given the displacement fields on
        different directions.

        Parameters
        ----------
        dx, dy, dz: float numpy array
            The displacement fields on the x, y, and z directions, respectively,
            with dimensions being `Ngrid`**3, and indices in C order.

        deepcopy: boolean, optional
            True for maintaining a copy of the displacements by the EZmock library,
            otherwise only the references are passed.
            In any case it is the user's responsibility to release memory for the
            displacement field arrays.
            Default: False.
        '''
        ...

    def populate_tracer(
        self,
        rho_c: float,
        rho_exp: float,
        pdf_base: float,
        sigma_v: float,
        ntracer: int,
        att_part: bool = False
    ) -> Tuple[NDArray, NDArray, NDArray, NDArray, NDArray, NDArray]:
        '''
        Generate the EZmock tracer catalog based on the density field and
        effective tracer bias model.

        Parameters
        ----------
        rho_c: float
            Critical dark matter density as the threshold of tracer formation.
            See Eq. (16) of the reference. Allowed range: [0, infty).
        rho_exp: float
            Exponential cut-off of the effective bias model.
            See Eq. (16) of the reference. Allowed range: (0, infty).
        pdf_base: float
            Base number of the power-law tracer probability distribution function.
            See Eq. (18) of the reference. Allowed range: (0, 1).
        sigma_v: float
            Standard deviation for random local peculiar motions.
            See Eq. (24) of the reference. Allowed range: [0, infty).
        ntracer: int
            The expected number of tracers to be generated.
        att_part: boolean, optional
            True for attaching tracers to dark matter partciles whenever possible.
            Default: False.

        Return
        ------
        x, y, z, vx, vy, vz: float numpy arrays
            Arrays for the coordinates and velocities of the tracers.

        Prerequisite
        ------------
        Function `create_dens_field_from_ic` or `create_dens_field_from_wn`
        or `create_dens_field_from_disp`.

        Reference
        ---------
        https://arxiv.org/abs/2007.08997
        '''
        ...

    def populate_tracer_to_file(
        self,
        rho_c: float,
        rho_exp: float,
        pdf_base: float,
        sigma_v: float,
        ntracer: float,
        fname: str,
        att_part: bool = False,
        rsd_fac: float = 0,
        header: bool = True
    ) -> None:
        '''
        Generate the EZmock tracer catalog based on the density field and
        effective tracer bias model, then save the catalog to an ASCII file.

        Parameters
        ----------
        rho_c: float
            Critical dark matter density as the threshold of tracer formation.
            See Eq. (16) of the reference. Allowed range: [0, infty).
        rho_exp: float
            Exponential cut-off of the effective bias model.
            See Eq. (16) of the reference. Allowed range: (0, infty).
        pdf_base: float
            Base number of the power-law tracer probability distribution function.
            See Eq. (18) of the reference. Allowed range: (0, 1).
        sigma_v: float
            Standard deviation for random local peculiar motions.
            See Eq. (24) of the reference. Allowed range: [0, infty).
        ntracer: int
            The expected number of tracers to be generated.
        fname: str
            Filename of the output tracer catalog.
        att_part: boolean, optional
            True for attaching tracers to dark matter partciles whenever possible.
            Default: False.
        rsd_fac: double, optional
            Positive for multiplying this factor to the z-velocity for
            redshift-space z coordinate, then write the redshift-space coordinates;
            otherwise, write the real-space coordinates and velocities.
            Default: 0.
        header: boolean, optional
            True for writing the header to the output file.
            Default: True.

        Prerequisite
        ------------
        Function `create_dens_field_from_ic` or `create_dens_field_from_wn`
        or `create_dens_field_from_disp`.

        Reference
        ---------
        https://arxiv.org/abs/2007.08997
        '''
        ...
