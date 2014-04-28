Numerical Method
================

During each time-integration step, the code calculate the fluxes :math:`\vec{F'}_{i+1/2,j}` at every ":math:`i`" half-point locations, and :math:`\vec{G'}_{i,j+1/2}` at every ":math:`j`" half-point locations. In order to obtain the flux terms properly treated with consideration of characteristics of wave propagation, MUSCL differencing should first be used to extrapolate the state vectors to every half point locations. After then AUSMPW+ scheme applies to those points for solving the inviscid flux terms. In addition, to evaluate the viscous flux terms, shear stress and heat flux terms should then be calculated at every half-points.

Flux vector evaluation
----------------------

Here, the description and formulation of AUSMPW+ scheme are not repeated. The viscous flux vectors in the generalized coordinates can then be evaluated by solving the following forms

.. math::

   \vec{F}_{V}' = \frac{1}{J}\left ( \xi_{x} \vec{F}_{V} + \xi_{y} \vec{G}_{V} \right )

   \vec{G}_{V}' = \frac{1}{J}\left ( \eta_{x} \vec{F}_{V} + \eta_{y} \vec{G}_{V} \right )

The above froms can be rearranged to the following form composed of shear stress and heat flux terms:

.. math::

   \vec{F}_{V}' \: \text{or} \: \vec{G}_{V}' = \frac{1}{J}\begin{bmatrix} 0 \\ m_{x} \tau_{xx} + m_{y} \tau_{xy}\\ m_{x} \tau_{xy} + m_{y} \tau_{yy}\\ m_{x}\left ( u \tau_{xx} + v \tau_{xy} - q_{x} \right ) + m_{y} \left ( u \tau_{xy} + v \tau_{yy} - q_{y} \right ) \end{bmatrix}

where,

.. math::

   m_{x} = \xi_{x} \:\:\: \text{and} \:\:\: m_{y} = \xi_{y} \:\:\: \text{for} \:\:\: \vec{F}_{V}'

   m_{x} = \eta_{x} \:\:\: \text{and} \:\:\: m_{y} = \eta_{y} \:\:\: \text{for} \:\:\: \vec{G}_{V}'

The nondimensional form of shear stress and heat flux terms are given by:

.. math::

   \tau_{xx} = \frac{2\mu}{3\text{Re}_{L}} \left ( 2 \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y} \right )

   \tau_{yy} = \frac{2\mu}{3\text{Re}_{L}} \left ( 2 \frac{\partial v}{\partial y} - \frac{\partial u}{\partial x} \right )

   \tau_{xy} = \frac{\mu}{\text{Re}_{L}} \left ( 2 \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right )

   q_{x} = - \frac{\mu}{(\gamma - 1)M_{\infty}^{2} \text{Re}_{L} \text{Pr}} \frac{\partial T}{\partial x}

   q_{y} = - \frac{\mu}{(\gamma - 1)M_{\infty}^{2} \text{Re}_{L} \text{Pr}} \frac{\partial T}{\partial y}

where :math:`M_{\infty}` is the free stream Mach number,

.. math::

   M_{\infty} = \frac{V_{\infty}}{\sqrt{\gamma R T_{\infty}}}

and the perfect gas equations of state become

.. math::

   p = (\gamma - 1) \rho e

   T = \frac{\rho M_{\infty}^{2} p}{\rho}

Note that the above equations should be evaluated by non-dimensional forms.

The coefficients of viscosity and thermal conductivity can be related to the thermodynamic variables according to the gas kinetic theory. Based on this theory, Sutherland's formulas can be applied to evaluate the mass and thermal diffusivity by solving:

.. math::

   \mu = C_{1} \frac{T^{3/2}}{T + C_{2}} \:\: , \:\:\: k = C_{3} \frac{T^{3/2}}{T + C_{4}}

where :math:`C_{1}` to :math:`C_{4}` are constants for a given gas. For air at moderate temperatures, :math:`C_{1} = 1.458 \times 10^{-6}`, :math:`C_{2} = 110.4`, :math:`C_{3} = 2.495 \times 10^{-3}` and :math:`C_{4} = 194` in SI units. Note that the temperature :math:`T` here must be here assumed to be dimenional variable in SI unit.

In this project, the Prandtl number is assumed to variable based on the following definition:

.. math::

   \text{Pr} = \frac{c_{p} \mu}{k}


Initial Conditions
------------------

At the beginning of simulation, the 2DNS code sets the initial condition. After then the code set the boundary conditions at every time step. The initial conditions at all grid points is set on the basis of following pre-specified flow quantities in nondimensional forms:

.. math::
   M = 2.0, \;\;\; \rho = 1.0,\; \;\; \; u = 1.0,\; \;\; \;v = 1.0, \; \;\; \; \gamma = 1.4,\; \;\; \; p = \frac{1}{\gamma M^{2}}, \;\;\; T = 1.0

The reference free stream conditions used to nondimensionalize the flow variables are given by:

.. math::

   \rho_{\infty} = 0.01 [\text{kg}/\text{m}^{3}], \;\;\; T_{\infty} = 300 [\text{K}], \;\;\; V_{\infty} = 694.44 [\text{m}/\text{sec}], \;\;\; L_{\infty} = 1.0 [\text{m}]


Boundary Conditions
-------------------

The flow is assumed to be coming in and blowing out at both inlet and outlet under a supersonic condition. Thus the following boundary conditions can be suitable:

Inflow: :math:`\vec{u}_{1,j}^{n}` = fixed at initial conditions at every time step

Outflow: :math:`\vec{u}_{imax,j}^{n} = \vec{u}_{imax - 1,j}^{n}` (1st order extrapolation for all n)

Inviscid wall (top): No velocity in the :math:`\eta` direction. The 2DNS code uses a 2nd order extrapolation that is described in the 3rd computer project assignment.

**Bottom wall**: 

In this project, two different wall boundary conditions are employed based on the treatment of wall temperature: adiabatic wall and isothermal wall boundaries. For the adiabatic wall boundary condition, the heat flux normal to the surface is enforced to be zero by:

.. math::

   T_{i,1} = T_{i,2}

For the isothermal wall boundary, the pre-specified wall temperature is applied to every :math:`j = 1` node points. In this project, 300 K is applied to the wall temperature as isothermal boundary condition.

In addition to the wall temperature BC, viscous wall boundary should be taken account. This can be made by assuming no-slip wall boundary for the moderate gas pressure. In this approach, the velocity right at the wall is set to zero. Then the surface pressure is calculated from the assumption that the normal component of the momentum equation is zero. This can be implemented by resolving the following relation for a non-orthogonal grid with no-slip condition:

.. math::

   \left ( x_{\xi}^{2} + y_{\xi}^{2} \right )\frac{\partial p}{\partial \eta} = \left ( x_{\xi} x_{\eta} + y_{\xi} y_{\eta} \right ) \frac{\partial p}{\partial \xi}

From the pressure and temperature resolved above, the density is then enforced by solving the gas equations of state.


Convergence Log (RMS error)
---------------------------

In order to see the convergence history, the 2DNS code calculates the following RMS error at every time step:

.. math::
   \text{RMS}^{n} = \sqrt{\frac{1}{N}\sum_{m=1}^{4} \sum_{i=1}^{imax} \sum_{j=1}^{jmax} \left [ \left ( \vec{U}_{i,j}^{n+1} - \vec{U}_{i,j}^{n} \right )^{2} \right ]}

This error log was used for checking convergence history only not for the termination of program. Every cases in this project was run by 40,000 iterations.
