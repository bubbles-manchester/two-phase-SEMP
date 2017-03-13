program SEMP
!$$ A program for solving the governing equations of a multiphase compressible 
!$$ viscoelastic fluid using the spectral element and marker particle methods.
!$$ Authors: S J Lind and P C Bollada
 
!$$ Important references are:
!$$ Lind and Phillips (2011) Int. J. Numer. Meth. Fluids, DOI: 10.1002/fld.2737
!$$ Lind (2010),A Numerical study of the effect of viscoelasticity on
!$$ cavitation and bubble dynamics, PhD Thesis, Cardiff University.
!$$ Bollada and Phillips(2008) Rheol. Acta, 47,719-725

!$$ Note: Eqn numbers mentioned herein refer to those in Lind and Phillips(2011)

implicit none


!$$$$$$$$$$$$$$$$$ ESSENTIAL PARAMETER ALLOCATION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

             integer    :: N, alphamax, betamax, bmax, imax, Gmax, F, alpha, beta, Nonz,a

            parameter     (N = 8, alphamax = 8, betamax = 8, bmax = 30, imax = 100,&
                           Gmax = (N*(alphamax + 1) - 1)*(N*(betamax + 1) - 1)*2, F = 50)

!$$ N: The order of Legendre polynomial
!$$ alphamax: Max number of spectral elements in x-direction
!$$ betamax: Max number of spectral elements in y-direction
!$$ alpha: x-direction label for an element
!$$ beta: y-direction label for element
!$$ imax: Max number of iterations in updating time step
!$$ Gmax: Dimension of Global stiffness matrix, G(Gmax,Gmax)
!$$ F: Max number of time steps

            integer     :: time
      double precision  :: dt, time2
      parameter           (dt = 5d-3)
      
!$$ dt: time step size
!$$ time: the time step counter (integer)
!$$ time2: the time

             integer    :: nnmax, MiCG5, NiCG5, status, flag, i, j, k, kk,b, rowIndex(Gmax+1)
             integer*8  pt(64)
 
            integer, allocatable  :: iCG5(:,:,:), columns(:)
            integer, allocatable  :: mmax(:)
   double precision, allocatable  :: values(:)


!$$ nnmax, MiCG5, NiCG5 concern the calculation of parts of the stiffness matrix (more details to follow).
!$$ status and flag are code "switches" used in time step iteration.
!$$ i,j,k,kk,b are frequently used indices.


!$$$$$$     PHYSICAL PARAMETER ALLOCATION $$$$$$$$$$$$$$$$$$$$$$$$$$$

      double precision  :: heat0, c2, nu
      parameter           (heat0 = 0d0, c2 = 1d2, nu = 0d0)
      
 !$$ heat0: heat constant (not currently used)
 !$$ c2: speed of sound squared 
 !$$ nu: 2nd viscosity coefficent (not currently used)

      double precision  :: mu_s0(2), mu_p0(2), lambda10(2), rho0(2)

 !$$ mu_s0: array containing initial solvent viscosities (1 for bubble, 2 for ambient fluid)
 !$$ mu_p0: Initial polymeric viscosity (1 for bubble, 2 for fluid)
 !$$ lambda10: Initial relaxation time (1 for bubble, 2 for fluid)
 !$$ rho0: Initial density (1 for bubble, 2 for fluid)

      double precision  :: xc, yc, rad
      parameter           (xc = 5d0, yc = 1.1d0, rad = 1d0) 
 
 !$$ speed: boundary wall speed (used for Couette flow test case)  
 !$$ xc, yc: initial coordinates of bubble centre
 !$$ rad: initial bubble radius

 !$$$$$$      MESH PARAMETERS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 !$$ The mesh consists of a refined rectangular centre around the bubble
 !$$ with elements fanning outwards to the domain boundary (see Lind and Phillips(2011), Fig. 4)

      integer           :: reg_grid, no_outer, mesh_info  

      parameter            (reg_grid = 0, no_outer = 1, mesh_info = 0)

      double precision  :: Len, Height, rl, rh, LL, HH, fl, fh, eps  
 
      parameter            (Len = 10d0, Height = 10d0, rl = 2d0, rh = 3d0, eps = 1d-11)

 !$$ reg_grid: A switch. If reg_grid=1, mesh forms regular rectangular arrangment instead.
 !$$ no_outer: Number of sets of outer fanning elements
 !$$ mesh_info: A switch: If mesh_info=1, all stiffness matrix calculations to do with mesh are done,
 !$$            so read them in from a file in this directory. If no such files exist, set mesh_info=0 
 !$$ rl, rh: ratio of length and height of inner mesh to outer (redundant if reg_grid=1)
 !$$ Len : Length of whole domain.
 !$$ Height: Height of whole domain.
 !$$ eps: A very small number that becomes useful at various points. 

!$$$$$$$$  MARKER PARTICLE PARAMETERS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      integer                    :: Nx, Ny, Np, TNP, Npx, Npy, Np1, N_jet

      parameter                    (Nx = 125, Ny = 125, Np = 4, N_jet = 1)

   double precision              :: dx_c, dy_c, dx_p, dy_p

   double precision, allocatable :: xp(:), yp(:), Cp(:,:), Upn(:,:)  

   double precision              :: xp_jet(N_jet), yp_jet(N_jet), U_jet(2)

!$$ TNP: Total Number of Marker Particles 
!$$ Nx: No of particle cells in x direction
!$$ Ny: No of particle cells in y direction
!$$ Np: No of particles per cell is Np**2
!$$ Npx: No of particles in x-direction
!$$ Npy: No of particles in y-direction
!$$ dx_c: Thickness of marker part. cell in x-direction
!$$ dy_c: Thickness of MP cell in y-direction
!$$ dx_p: intial distance between particles x-direction
!$$ dy_p: intial distance between particles y-direction
!$$ xp,yp: Position of particle labelled (:)
!$$ Cp: Colour of particle labelled (:)
!$$ Upn: Previous velocity of particle labelled (:)

!$$$$$$      FILE NAMES  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$ These are character strings for use when writing various files

      character(len=30) :: filename, filename2
      character(len=5)  :: al_name, be_name, N_name, nout_n, rl_n, rh_n,&
                           lamb_name, mu_p_name, mu_s_name, len_name,&
                           rhob_name, h_name, dt_name, snd_name
      
      parameter (filename ='solution.dat', filename2 ='mesh_test.dat' )


!$$$$$$   VARIABLE DECLARATION  $$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      double precision  & 
     
      TG4(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !Global solvent stress tensor
      gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !Global velocity gradient
      graduN(0:N*(alphamax+1),0:N*(betamax+1),2,2),& ! Global vel grad at previous time step
      divu(0:N*(alphamax+1),0:N*(betamax+1)),&	!Global velocity divergence
      U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),& !Local velocity
      Un(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),& !Local velocity at previous time step x^n
      W(0:N),P(0:50),L(0:N),L1(0:N),& !Weights, GLL points, Legendre polynomial and derivative
      h1(0:N,0:N),h5(0:N,0:N,0:N,0:N,2),& !Lagrange interpolant derivative and procuct of derivative
      J4(0:alphamax,0:betamax,0:N,0:N),& !Jacobian from physical to parent element
      Z6(0:alphamax,0:betamax,0:N,0:N,2,2),& !Z matrix (Eqn 38)
      B3_Left(1:N*(betamax + 1) - 1, 3),&  !velocity boundary conditions - one for each side of domain
      B3_Right(1:N*(betamax + 1) - 1, 3),&
      B3_Top(0:N*(alphamax + 1), 3),&
      B3_Bot(0:N*(alphamax + 1), 3),&
      xq(0:alphamax,0:betamax,5),&
      yq(0:alphamax,0:betamax,5),& !Vertices of quadrilateral spectral elements - xq(5) = xq(1), for ease of coding
      divu_max,& !max values of some tensors
      TG4_max,&
      CG5_max,&
      Pstr_max,&
      ZG3_max,&
      mu_p_sclr, mu_sclr, time3,&
      Clij(0:alphamax,0:betamax,0:N,0:N,2),& !Local colour function
      PstressL(0:alphamax,0:betamax,0:N,0:N,2,2),&!Local polymeric stress
      PstressnL(0:alphamax,0:betamax,0:N,0:N,2,2),&!Local polymeric stress at previous t-step, x^n
      Norm_Pstress(0:N*(alphamax+1),0:N*(betamax+1)),&!Normal polymeric stress
      Shear_Pstress(0:N*(alphamax+1),0:N*(betamax+1)),&!Shear polymeric stress
      Norm_TG4(0:N*(alphamax+1),0:N*(betamax+1)),&!Normal extra stress
      Shear_TG4(0:N*(alphamax+1),0:N*(betamax+1)),&!Shear extra stress
      speed(0:N*(alphamax+1),0:N*(betamax+1)),&!velocity magnitude 
      x_bub_c,y_bub_c,Q_bub_c,&!Coordinates (not important)
      De_BEM,Re_BEM,mu_BEM,& !Equivalent Re, De, mu using BEM definitions 
      Mass, Mass_error !Errors in particle mass when doing test problems
                       

      double precision  &
   
      UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),& !Global velocity and density
      AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),& !Global tensor component of stiffness matrix (Eqn 36)
      CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&!Global tensor comp. of stiff matrix (Eqn 37)
      TG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&!
      TG7(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&!Various comp. of stiff matrix (Eqn 44)
      TG8(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&!
      HG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&!
      Qn(0:N*(alphamax+1),0:N*(betamax+1)),& ! Previous global density
      fG3(0:N*(alphamax+1),0:N*(betamax+1),3),&!Part of known vector on RHS 
      YG3(0:N*(alphamax+1),0:N*(betamax+1),2),& !Part of known vector on RHS
      VG3(N*(alphamax+1)-1, N*(betamax+1)-1,2),&!Known vector on RHS (Eqn 47)
      mu_t(0:N*(alphamax+1),0:N*(betamax+1)),& !"Total viscosity" (not a viscosity but a sum of various fluid parameters)
      mu_s(0:N*(alphamax+1),0:N*(betamax+1)),& !solvent viscosity
      lambda1(0:N*(alphamax+1),0:N*(betamax+1)),& !relaxation time
      Cij(0:N*(alphamax+1),0:N*(betamax+1),2),& !Global interpolated colour function
      xg(0:N*(alphamax+1),0:N*(betamax+1)),& !Global GLL points coordinates (xg,yg)
      yg(0:N*(alphamax+1),0:N*(betamax+1)),& 
      mu_p(0:N*(alphamax+1),0:N*(betamax+1)),& !polymeric viscosity
      Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !Global polymeric stress
      Bn(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !Part of polymeric stress
      Pstressn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!Previous global polymeric stress
      Extra_Stress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!Global extra_stress=solvent+polymeric stress 
      DG3(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !Part of polymeric stress
      ZG3(0:N*(alphamax+1),0:N*(betamax+1),2),&!part of polymeric stress
      F_def(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!part of polymeric stress
      F_defn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!part of polymeric stress
      S_convn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!part of polymeric stress
      S_conv(0:N*(alphamax+1),0:N*(betamax+1),2,2),&!part of polymeric stress
      GammaN(0:N*(alphamax+1),0:N*(betamax+1),2,2),& !component matrices of stiffness matrix (Eqn 44)
      GammaP(0:N*(alphamax+1),0:N*(betamax+1),2,2,2,2),&
      S(0:3),&!Array used for summing velocity components
      xi,zeta,&!coordinates on parent element
      xx,yy,&!temp x,y coordinates of particles
      num_u,&!random variable
      vol,&!volume of bubble
      rad_av,&!average radius of bubble
      cputime1, cputime2, cputime3, cputime4, cputime5, cputime6!floating point numbers for checking cpu time

   double precision :: MGlobal(Gmax,Gmax), Vglobal(Gmax) !Global stiffness matrix and RHS vector
       integer      :: ipiv(Gmax) !Integer array for matrix computation using PARDISO.

 !!!$$$$$$$$$$ INITIALISE VARIABLES $$$$$$$$$$$$$$$$$

      
      TG4 = 0d0
    gradu = 0d0
   graduN = 0d0
     divu = 0d0	  
        U = 0d0
       Un = 0d0
        W = 0d0
        P = 0d0
        L = 0d0
       L1 = 0d0
       h1 = 0d0 
       h5 = 0d0
       J4 = 0d0
       Z6 = 0d0
  B3_Left = 0d0
 B3_Right = 0d0
   B3_Top = 0d0
   B3_Bot = 0d0
       xq = 0d0
       yq = 0d0 
      UG3 = 0d0
      DG3 = 0d0
      AG4 = 0d0
      CG5 = 0d0
      TG6 = 0d0
      TG7 = 0d0
      TG8 = 0d0
      HG6 = 0d0
       Qn = 0d0
      fG3 = 0d0
      YG3 = 0d0
      ZG3 = 0d0
      VG3 = 0d0
  MGlobal = 0d0
  Vglobal = 0d0
    F_def = 0d0
   F_defn = 0d0
  S_convn = 0d0
   S_conv = 0d0
       Bn = 0d0
 Pstressn = 0d0
  Pstress = 0d0
  GammaN = 0d0
  GammaP = 0d0
       xg = 0d0
       yg = 0d0
      Cij = 0d0
      Clij = 0d0
 PstressnL = 0d0
    
 
     ipiv = 0

!!!$$$ Set material parameters...$$$$$$

   mu_s0(1) = 1d-5    !$ Solvent viscosity. Label 1 refers to bubble
   mu_s0(2) = 1d0     !$ Label 2 refers to ambient fluid

   mu_p0(1) = 0d0     !$ polymer viscosity
   mu_p0(2) = 1d0

   lambda10(1) = 0d0  !$ relaxation time
   lambda10(2) = 1d0

   rho0(1) = 0d0      !$ initial density
   rho0(2) = log(4d0)

    xp_jet = 5d0      !$ initial coordinates of liquid jet
    yp_jet = yc + rad
     U_jet = 0d0      !$ initial velocity of liquid jet

!!!$$$$$$$$$$$ PROGRAM START $$$$$$$$$$$$$$$$$$$$
 !!!    open(*, file='error_messages_var_mu.txt')
 !!!    open(*, file='run_notes_var_mu_static_link.txt')
     write(*,*) 'program started'
     write(*,*) 'material variables'
     write(*,*) 'mus0 1,2=',mu_s0(1),mu_s0(2)
     write(*,*) 'mup0 1,2=',mu_p0(1),mu_p0(2)
     write(*,*) 'lambda10 1,2=',lambda10(1),lambda10(2)
     write(*,*) 'rho0 1,2=',rho0(1),rho0(2) 
     write(*,*) 'c2=',c2   
     write(*,*) 'height=',yc
     write(*,*) 'time step, dt=',dt 

!$$ Write parameters to character strings to use in filenames

      write(al_name,100) alphamax
      write(be_name,100) betamax
      write(N_name,100) N
      write(nout_n,100) no_outer
      write(rl_n,101) rl
      write(rh_n,101) rh
      write(lamb_name,104) lambda10(2)
      write(mu_p_name,101) mu_p0(2)
      write(mu_s_name,101) mu_s0(2)
      write(rhob_name,101) rho0(2)
      write(h_name,101) yc
      write(dt_name,101) dt
      write(snd_name,103) c2
!!      write(len_name,102) Len
    
      100 format(I3)
      101 format(f5.3)
      102 format(f4.1)
      103 format(f5.0)
      104 format(f5.1)

 !$$ Files containing mesh/stiffness matrix data.   
 !$$ (constant for a given mesh - saves working them out again for future runs)
       open(11, file = 'AG4_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(12, file = 'CG5_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(13, file = 'NiCG5_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(14, file = 'MiCG5_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(15, file = 'nnmax_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(16, file = 'mmax_al='//al_name//'be='//be_name//'N='//N_name//'.dat')
       open(17, file = 'iCG5_al='//al_name//'be='//be_name//'N='//N_name//'.dat')

 !!$ Solution files

   open(998, file='colour_imp_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
   open(99, file='soln_imp_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                 '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
   open(888, file='Pstress_imp_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
   open(889, file='extra_stress_imp_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
   open(892, file='gradu_EF_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
  
   open(991, file='mu_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')

   open(32,file='decomposed_stress_smth_den_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
  
   open(33,file='2BEM_parameters_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')

   open(77,file='2Bubble_mass_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')

    open(88,file='Bubble_vol_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')

    open(777,file='jet_vel_pos_c2='//snd_name//'_lam='//lamb_name//'_mu_p='//mu_p_name//&
                  '_mu_s='//mu_s_name//'_rhob'//rhob_name//'_h='//h_name//'_dt='//dt_name//'N='//N_name//'.dat')
     write(777,*) 'VARIABLES = "time","x jet","y jet","u jet","v jet"'



!!$ Calculating GLL points, interpolants etc...

      call jacobl(0d0,0d0,N,P) !Calculate GLL points
      call assignL(N,P,L) ! Calc Legendre polys
      call calcW(N,L,W) !Calculate weights
      call calcL1(N,L1) ! Calc derivative of Legendre polys
      call calch1(N,L,L1,P,h1) !Calculate Lagrange interpolant derivative

 !!$ Calculate mesh - find element vertices (xq, yq)
       LL = Len
       HH = Height
       fl = rl
       fh = rh

      call mesh_generator(LL, HH, alphamax, betamax, fl, fh, no_outer, reg_grid, xq, yq) 

 !!$ Calc Jacobian
      call calcJ4(N, P, alphamax, betamax, xq, yq, J4)

 !!$ Calc Z matrix (Eqn 38)
      call calcZ6(N, P, alphamax, betamax, xq, yq, Z6)
   
 !!$ Calc product of interpolant derivative and interpolant (in Eqn 37)
      call calch5(N,h1,h5)
       
 !!$ If mesh matrices have not been calculated previously do so... 

     if(mesh_info.eq.0)then

     ! Calc matrix A (Eqn 36)
      call calcAG4(N,alphamax,betamax,W,J4,AG4)

     ! Calc matrix C (Eqn 37)
      call calcCG5(N,W,h5,alphamax,betamax,Z6,CG5)

      call precalciCG5(N,alphamax,betamax,CG5,NiCG5,MiCG5, eps)
 
      allocate (iCG5(NiCG5, MiCG5, 2), stat = status)
      allocate (mmax(NiCG5), stat = status)
 
      call calciCG5(N, alphamax, betamax, CG5, nnmax, mmax, NiCG5, MiCG5, iCG5, eps)
     
!      call write_mesh_details(N, alphamax, betamax, AG4, CG5, NiCG5, MiCG5, nnmax, mmax, iCG5)

      else 

 !!$ ...else read in previously written files for previous calculated mesh

      call read_mesh_details_a(N, alphamax, betamax, AG4, CG5, NiCG5, MiCG5, nnmax)

      allocate (iCG5(NiCG5, MiCG5, 2), stat = status)
      allocate (mmax(NiCG5), stat = status)
 
      call read_mesh_details_b(NiCG5, MiCG5, mmax, iCG5)

      write(*,*) 'have read mesh details'

      endif

       time = 1
      time2 = 0d0 ! Initialise time variables
      time3 = dt


     !$ Calculate total number of particles (TNP), given particle cell width and number of 
     !$ particles per cell
      call calc_no_particles(TNP, Npx, Npy, Nx, Ny, Np, dx_p, dy_p, dx_c, dy_c, Len, Height)

      write(*,*) 'TNP=', TNP

     !$ Allocate particle arrays according to total particle number
     !$ (x position, y pos, colour, velocity, respectively)
      allocate (xp(TNP), yp(TNP), Cp(TNP,2), Upn(TNP,2))

     !$ Calc (xg,yg): global coordinates of GLL points
      call initUG3(N,alphamax,betamax,UG3,U,imax,xg,yg,xq,yq,P)   

     !$ Initialise postion of marker particles
      call calc_ptcl_pstn(TNP, Np1, Npx, Npy, xp, yp, Cp, xc, yc, rad, dx_p, dy_p)  
 
     !$ Initialise colour of particles
      call calc_colour(N, TNP, alphamax, betamax, Cp, Cij, dx_c, dy_c, xp, yp, xg, yg)

     !$ Send global colour array to local array
      call calc_local_colour(N, alphamax, betamax, Cij, Clij) 

     !$ Initialise local velocity and density
      call initU(N,alphamax,betamax,U,imax,heat0, xq, yq, P, xc, yc, rho0, rad, Clij)

     !$ Initialise global velocity and density field
      call initUG3(N,alphamax,betamax,UG3,U,imax,xg,yg,xq,yq,P)

     !$ Calculate initial bubble mass (for validation purposes)
      call calc_bubble_mass(N,imax,alphamax,betamax,W,J4,U,Mass,Mass_error,rad,Clij,rho0(1),time)

     !$ Write solutions (velocity, density, stresss) to output files
     call writepoints_grid(alphamax, betamax, P, N, U, J4, Un, xq, yq, time2, imax, i-1,PstressnL)
                         
      write(*,*) 'writepoints grid' 
     !$ Write colour to output file
     call writepoints_colour(N, alphamax, betamax, Cij, xg, yg, time2)

      write(*,*) 'writepoints_colour'
 
     !$ Note cputime up to this point  
     call cpu_time(cputime3)


  do while (time.le.F) !!$$ BEGIN THE TIME LOOP!!!

  !$ Write mass and mass_error to file (for test cases)
    write(77,*) time2, Mass, Mass_error

  !$ Calc bubble vol (notable physical quantity)
     call calc_bubble_vol(N,alphamax,betamax,W,J4,Clij,vol,rad_av)

    !$ Write bubble volume and average radius to output file
     write(88,*) time2, vol, rad_av

   

     call cpu_time(cputime1)

    !$ Calc bubble centroid and density at bubble centroid
     call calc_bubble_xy(N,imax,alphamax,betamax,TNP,Np1,P,xp,yp,xq,yq,U,Cp,&
                        x_bub_c,y_bub_c,Q_bub_c)

    !$ Calc equivalent Re, De numbers using BEM definitions (see Lind 2010)
     call calc_BEM_parameters(c2,rho0,Q_bub_c,lambda10(2),mu_p0(2),De_BEM,Re_BEM,mu_BEM)
 
    !$ Write Re, De to file to show variation in time
     write(33,*) time2,De_BEM,Re_BEM

 
   !$$ Calculate material variables at each GLL point

     call calc_mu_s(N, TNP, alphamax, betamax, Cij, mu_s, mu_s0) !solvent viscosity
     call calc_mu_p(N, TNP, alphamax, betamax, Cij, mu_p, mu_p0) !polymeric viscosity
     call calc_lambda(N, TNP, alphamax, betamax, Cij, lambda1, lambda10)!relaxation time
     call calc_mu_t(N, alphamax, betamax, mu_s, mu_p, lambda1, dt, mu_t)!total viscosity

   !$$ Calc material variable arrays (GammaN, GammaP (Eqns 45,46)) for stiffness matrix
     call calc_Gamma(N, alphamax, betamax, GammaN, GammaP, mu_t, lambda1, Pstress,dt)
   
     call cpu_time(cputime2)

    write(*,*) 'calcing material variables',(cputime2-cputime1)

   !$$ Calculate components of stiffness matrix

      call cpu_time(cputime1)

       call calcTGG(N,alphamax,betamax,nnmax,mmax,NiCG5,MiCG5,iCG5,AG4,CG5,TG6,TG7,TG8,GammaN,GammaP)

      call cpu_time(cputime2)

     write(*,*) 'time in TGG',(cputime2-cputime1)

    !$$ Begin assembly of stiffness matrix 
      call adjHG6(N,alphamax,betamax,c2,dt,TG6,TG7,TG8,HG6)
 
      call cpu_time(cputime1)

     write(*,*)  'time in adjHG6',(cputime1-cputime2)

   !!$ Complete assembly of global stiffness matrix

      call calcglobalmatrix(N,alphamax,betamax,Gmax,HG6,Mglobal,Nonz)
    
     call cpu_time(cputime2)

     write(*,*) 'time in calc global matrix',(cputime2-cputime1)

     allocate (values(Nonz),columns(Nonz))

     call cpu_time(cputime1)

    !$$ Recast stiffness matrix in a form suitable for PARDISO

    call calcStorageArrays(Gmax,Nonz,Mglobal,values,columns,rowIndex)

    call cpu_time(cputime2)

    write(*,*) 'time in calc storage array',(cputime2-cputime1)

    !$$ Factorise Matrix using PARDISO

    call Pardiso_Fact(pt, Nonz, Gmax, values,columns,rowIndex)

!$$$$$ Alternatively can use NAG library routine below
!$$$$$ Permits greater portability as PARDISO only available to
!$$$$$ me on MERLIN

    !!$ Factorise stiffness matrix using Chloesky factorisation
    !!$ (This is where the code spends most time)

    !!$  call CholeskyFact(Gmax,MGlobal,ipiv) !Uses NAG libs

      call cpu_time(cputime1)

     write(*,*) 'time in pardiso fact',(cputime1-cputime2)

     !!$ Initialise boundary conditions
       call calcB3(N,alphamax, betamax, B3_Left, B3_Right, B3_Top, B3_Bot)

         i = 0 !Counter for iteration
      flag = 2 !switch to exit iteration once converged

      !!$ A routine that rewrites a global tensor as its local equivalent
      !!$ (In this case - Pstress -> PstressL)
       call global_2_local_tensor(N,alphamax,betamax,Pstress,PstressL)

       call cpu_time(cputime2)

      write(*,*) 'time for boundary condit. and local tensor calc',(cputime2-cputime1)

       call cpu_time(cputime3)

   
    !!$ ITERATIVE PROCEDURE (AS IN SEC. 3.1.1) TO CALC VELOCITY
       do while(flag.ge.1.and.i.le.imax)

         
      !!$ Calc previous velocity at x^n
         call calcUn(alphamax,betamax,N, xq, yq, P, dt, i, U, Un, Qn, flag, imax, speed, Len, Height,&
                     PstressL,Pstressn,PstressnL,eps,time)
  
           write(*,*) 'calc un'  
        !!$ Calc global matrix associated with previous velocity
         call calcfG3(N,i,alphamax,betamax,Un,W,J4,fG3,imax)   

             write(*,*) 'calc fg3'  

         !!$ Calc RHS vector of "knowns" (Eqn 47)
         call calcVG3(imax,N,dt,c2,alphamax,betamax,B3_Left, B3_Right, B3_Top, B3_Bot,&
                      CG5,Qn,fG3,HG6,UG3,YG3,ZG3,VG3,Pstressn,Bn,lambda1)
       
             write(*,*) 'calcvg3'  

         !!$ Assemble global RHS vector (Eqn 47)
         call calcglobalvector(N,alphamax,betamax,Gmax,VG3,HG6,VGlobal)

             write(*,*) 'calc global vector'  

         call cpu_time(cputime1)

      
         !!$ Solve linear system (Eqn 43) using PARDISO
          call Pardiso_Solve(pt, Nonz, Gmax,values,columns,rowIndex,VGlobal)

        !$$$$$$ Alternatively can use NAG routine below (greater portability)

         !!$ Solve linear system Mu=V for global velocity (using lapack routine)
         !!$ call solveSym4U(Gmax,MGlobal,VGlobal,ipiv)

        !$$$$$$ End of NAG routine call

 
         call cpu_time(cputime2)

           write(*,*)'time in solveSym4u',(cputime2-cputime1)  

         !!$ Update the new found global velocity
         call updateUG3(N,Gmax,i+1,dt,alphamax,betamax,VGlobal,&
                        Qn,CG5,AG4,UG3,HG6,imax,1)

               write(*,*) 'calc updateUG3'  
        
         !!$ Update the local velocity
         call updateU(N,i+1,alphamax,betamax,UG3,U,imax,1)

                write(*,*) 'calc updateU'  


         i = i + 1 !Increase iteration counter

       enddo !!$ END OF ITERATION (EITHER CONVERGED - flag=1, or NOT - i=imax)

       if(i.ge.imax-1)then
        write(*,*)'iteration for Un exceeded'
        stop
       endif

       call cpu_time(cputime4)

       write(*,*) 'time in iteration',(cputime4-cputime3)
     

     !$ Out of iteration: now update global density    
      call updateUG3(N,Gmax,0,dt,alphamax,betamax,VGlobal,&
                     Qn,CG5,AG4,UG3,HG6,imax,0)
     !$ Update local density
      call updateU(N,0,alphamax,betamax,UG3,U,imax,0)

       write(*,*) 'done update U'

    !$ Calc solvent stress TG4
      call calcVG2X(0,imax,N,alphamax,betamax,mu_s,nu,CG5,UG3,AG4,TG4,&
                    gradu,graduN,divu)

      write(*,*) 'calc VG2X'

     !$ Calc polymeric stress, Pstress
     
      write(*,*) 'oldPstress',maxval(Pstress),minval(Pstress)

      call calcOldBstress3(N,alphamax,betamax,dt,gradu,mu_p,lambda1,Pstress,Pstressn)

      write(*,*) 'newPstress',maxval(Pstress),minval(Pstress)

     !$ Decompose into shear and normal stress components as in Bollada and Phillips(2008)
     !$ (Not for any direct purpose, other than interest)
      call calcStressDecomp(imax,N,alphamax,betamax,TG4,Pstress,UG3,&
                            Norm_Pstress,Shear_Pstress,Norm_TG4,Shear_TG4)

     
     !$ Calc matrix YG3 (part of RHS vector VG3 (Eqn 47))
      call calcYG3(N,imax,alphamax,betamax,UG3,CG5,TG4,Pstress,Extra_Stress,YG3)

          write(*,*) 'calc yg3'

      call cpu_time(cputime1)

     !$ With new velocity, update postion of all marker particles in Lagrangian fashion
     !$ (using simple Euler scheme)
      call update_ptcl(TNP, alphamax, betamax, N, imax, P, xq, yq, Len, Height, xp, yp, U, U_jet, dt, eps,time,0)

     !$ Output position and velocity of liquid jet.
       write(777,*) time2, xp_jet(1), yp_jet(1), U_jet(1), U_jet(2)
    
     !$ Calc velocity of "jet" particle and also update
      call update_ptcl(N_jet, alphamax, betamax, N, imax, P, xq, yq, Len, Height, xp_jet, yp_jet, U, U_jet, dt, eps,time,1)

    

      call cpu_time(cputime2)

            write(*,*) 'time to  update ptcl',(cputime2-cputime1)    

     !$ Calc new interpolated colour function for new particle positions
      call calc_colour(N, TNP, alphamax, betamax, Cp, Cij, dx_c, dy_c, xp, yp, xg, yg)

     !$ Calc local version of global colour array
      call calc_local_colour(N, alphamax, betamax, Cij, Clij) 

     !$ Calc bubble mass (for validatory test cases)
      call calc_bubble_mass(N,imax,alphamax,betamax,W,J4,U,Mass,Mass_error,rad,Clij,rho0(1),time)

       write(*,*) 'time step',time,time2

      !$ Update time step counters
       time = time + 1
      time2 = time2 + dt
      time3 = time3 + dt


     !$ Write solution to various output files
      if((mod((time-1),40).eq.0).or.(time.le.100))then 
        call writepoints_grid(alphamax, betamax, P, N, U, J4, Un, xq, yq, time2, imax, i-1, PstressnL)
        call writepoints_colour(N, alphamax, betamax, Cij, xg, yg, time2)
        call writepoints_tensors(N, alphamax, betamax, Extra_Stress, xg, yg, time2, 889)    
      endif

     !$ Deallocate memory used by PARDISO
    call Pardiso_Release_Memory(pt,Gmax)

    deallocate(values)
    deallocate(columns)

    call cpu_time(cputime6)

    write(*,*) 'time for time step',(cputime6-cputime5)


  
  enddo !$$ END OF TIME LOOP

      call cpu_time(cputime4)

     write(*,*) 'time for run',(cputime4-cputime3)

   !$$ Close various output files
        close(998)
        close(99)
        close(888)
        close(889)
        close(890)
        close(891)
       
     !$$ deallocate various arrays   
        deallocate(iCG5,stat = status)
	deallocate(mmax,stat = status)

 
end program SEMP
      

   !$ Calc Legendre polynomials
      double precision recursive function calcL(N,x) result(answer)
      implicit none
      integer N
      double precision x
      if (N.eq.0) then
         answer=1d0
      else if (N.eq.1) then 
         answer=x
      else
      answer=1d0/N*((2d0*N-1)*x*calcL(N-1,x)-(N-1)*calcL(N-2,x))
      endif
      return
      end function calcL


      subroutine assignL(N,P,L)
      implicit none
      integer i,N
      double precision L(0:N),P(0:N)
      double precision calcL
      do i=0,N
         L(i)=calcL(N,P(i))
      enddo
      end               

   !$ Calc weights for quadrature rule
      subroutine calcW(N,L,W)
      implicit none
      integer i,N      
      double precision L(0:N),W(0:N)
      do i=0,N
         W(i)=2d0/N/(N+1d0)/L(i)/L(i)
      enddo
      end      
 
    !$ Calc derivative of Legendre poly
      subroutine calcL1(N,L1)
      implicit none
      integer i,F,N
      double precision L1(0:N)
      F=1
      if(N.eq.2*int(N/2))F=-1
      L1(N)=N*(N+1d0)/2d0
      L1(0)=F*L1(N)
      do i=1,N-1
         L1(i)=0d0
      enddo
      end

    !$ Lagrange interpolant h (at interpol. points)
      double precision function h(i,j)
      implicit none
      integer i,j
      if(i.eq.j) then
         h=1d0
      else
         h=0d0
      endif
      end
      
     !$ Derivative of interpolant
      subroutine calch1(N,L,L1,P,h1)
      implicit none
      integer i,j,N
      double precision h1(0:N,0:N),L(0:N),L1(0:N),P(0:N)
      do i=0,N
         do j=0,N
            if(i.eq.j)then
               h1(i,j)=0.5d0*L1(i)/L(i)
            else
               h1(i,j)=L(j)/(L(i)*(P(j)-P(i)))
            end if
         enddo
      enddo
      end



     !$ Product of interpolants (see Eqn 37)
      subroutine calch5(N,h1,h5)
      implicit none
      integer i,j,k,l,N
      double precision h,h1(0:N,0:N)
      double precision h5(0:N,0:N,0:N,0:N,2)
      do i=0,N
         do j=0,N
            do k=0,N
               do l=0,N
                  h5(k,i,l,j,1)=h1(k,i)*h(l,j)
                  h5(k,i,l,j,2)=h(k,i)*h1(l,j)
               enddo
            enddo
         enddo
      enddo
      end


    !$ Initialise local velocity and density
     subroutine initU(N,alphamax,betamax,U,imax,heat0, xq, yq, P, xc, yc, rho0, rad, Clij)
      implicit none
      integer          :: al,be,i,j,alphamax,betamax,N,imax, n_crit
      double precision :: heat0,&
                          U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                          P(0:N), Len, Height, xc, yc, xx, yy, rho0(2),&
                          rad, x, y, rr, xq(0:alphamax,0:betamax,5),&
                          yq(0:alphamax,0:betamax,5),&
                          Clij(0:alphamax,0:betamax,0:N,0:N,2),aa,bb
                        
    !$ Density given by U(:::0::) entry   
    !$ Horizontal velocity given by U(:::1::)
    !$ Vertical velocity given by U(:::2::)
    !$ Temperature in U(:::3::) (not used here)

    aa = rad - 0.1d0
    bb = rad + 0.1d0

    do al = 0, alphamax
       do be = 0, betamax
         do i = 0, N
           do j = 0, N

              xx = x(al, alphamax, be, betamax, P(i), P(j), xq)
              yy = y(al, alphamax, be, betamax, P(i), P(j), yq)

              rr = dsqrt((xx - xc)**2 + (yy - yc)**2)
             
            !$Initial density
              if(rr.lt.aa)then

                U(0,al,be,0,i,j) = 0d0

              elseif((rr.le.bb).and.(rr.ge.aa))then

                U(0,al,be,0,i,j) = (rho0(2)/(bb-aa))*(rr - aa)
             
              else
       
                U(0,al,be,0,i,j) = rho0(2)

              endif
 
             !$ Initial velocities 
              U(0,al,be,1,i,j) = 0d0!!!!speed*(y(al, alphamax, be, betamax, P(i), P(j), yq))/Height
              U(0,al,be,2,i,j) = 0d0            
             
             !$ Initial temp
              U(0,al,be,3,i,j) = heat0
           
             enddo
            enddo
         enddo
      enddo
      endsubroutine


    !$ Initialise global velocity and density  
      subroutine initUG3(N,alphamax,betamax,UG3,U,imax,xg,yg,xq,yq,P)
      implicit none
      integer           :: i,j,a,alphamax,betamax,N,Nhat,Mhat,imax,al,be,k,l
      double precision  :: heat0,&
                           UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                           U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),P(0:N),&
                           xg(0:N*(alphamax+1),0:N*(betamax+1)),&
                           yg(0:N*(alphamax+1),0:N*(betamax+1)),&
                           xq(0:alphamax,0:betamax,5),&
                           yq(0:alphamax,0:betamax,5),x,y
                          
     
      Nhat = N*(alphamax+1)
      Mhat = N*(betamax+1)


    do a = 0,3
      do i = 0, Nhat !Loop over all GLL points in x
         do j = 0, Mhat !Loop over all GLL points in y


  !$ Mapping from global indices i,j, to local coordinates, al, be, k,l
              al = (i-1)/N   !unfortunately this does not work for N=1
              be = (j-1)/N
               k = mod(i+al,N+1)
               l = mod(j+be,N+1)
         
    !$ Global positions of mesh nodes (xg,yg) from local.
          if(a.eq.0)then
            xg(i,j) =  x(al, alphamax, be, betamax, P(k), P(l), xq)
            yg(i,j) =  y(al, alphamax, be, betamax, P(k), P(l), yq)
         endif

      !$ Set global velocity from local
 
             UG3(0,a,i,j) = U(0,al,be,a,k,l)

       
           enddo
         enddo
        enddo

      
      return
      end
            

!$$$$$$ NEED

      subroutine jacobl(alpha,beta,N,xjac)
!
!   computes the Gauss-Lobatto collocation points for the jacobi polynomials
!
!   n:      degree of approximation
!   alpha:  parameter in jacobi weight
!   beta:   parameter in jacobi weight
!   
!   xjac:   output array with the Gauss-Lobatto roots
!           they are ordered from the largest (+1.0) to the smallest (-1.0)???
!     smallest to largest I believe!
!
      implicit none
      integer N
      double precision alpha,beta,xjac(0:50)
!
      integer np,nh,npp,i,j,jm,k,kstop
      double precision pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,&
                       pnp1m,pdnp1m,pnm,pnm1m,pdnm,det,&
                       pnm1,cs,x,pnp1,pdnp1,pn,pdn,&
                       rp,rm,ag,bg,dth,cd,sd,ss,poly,pder,&
                       cssave,delx,epsg,recsum,hulpar(0:64),pi
      double precision alp,bet,rv
      common /jacpar/ alp,bet,rv
      data kstop/10/
      data epsg/1.0d-25/
!
      pi = 4d0*datan(1d0)
      alp = alpha
      bet = beta
      rv = 1 + alp
      np = n+1
!
!   compute the parameters in the polynomial whose roots are desired
!
      call jacobf(np,pnp1p,pdnp1p,pnp,pdnp,pnm1p,pdnm1,1d0)
      call jacobf(np,pnp1m,pdnp1m,pnm,pdnm,pnm1m,pdnm1,-1d0)
      det = pnp*pnm1m-pnm*pnm1p
      rp = -pnp1p
      rm = -pnp1m
      ag = (rp*pnm1m-rm*pnm1p)/det
      bg = (rm*pnp-rp*pnm)/det
!
      xjac(1) = 1d0
      nh = n
!
!   set-up recursion relation for initial guess for the roots
!
      dth = pi/(2*n+1)
      cd = cos(2d0*dth)
      sd = sin(2d0*dth)
      cs = cos(dth)
      ss = sin(dth)
!
!   compute the first half of the roots by polynomial deflation
!
      do 39 j = 2,nh
        x = cs
        do 29 k=1,kstop
          call jacobf(np,pnp1,pdnp1,pn,pdn,pnm1,pdnm1,x)
          poly = pnp1+ag*pn+bg*pnm1
          pder = pdnp1+ag*pdn+bg*pdnm1
          recsum = 0d0
          jm = j-1
          do 27 i=1,jm
            recsum = recsum + 1d0/(x-xjac(i))
   27     continue
          delx = -poly/(pder-recsum*poly)
          x = x + delx
          if (abs(delx).lt.epsg) goto 30
   29   continue
   30   continue
        xjac(j) = x
        cssave = cs*cd-ss*sd
        ss = cs*sd+ss*cd
        cs = cssave
   39 continue
      xjac(np) = -1d0
      npp = n+2
!
!   use symmetry for second half of the roots
! 
      do 112 i=1,N+1
        hulpar(N-i+1) = xjac(i)
  112 continue
      do 113 i=0,N
        xjac(i) =hulpar(i)
  113 continue
      if (N.eq.2*int(N/2.))xjac(N/2)=0. 
      return
      end
!

!$$$$$$ NEED
      subroutine jacobf(n,poly,pder,polym1,pderm1,polym2,pderm2,x)
!
!   computes the jacobi polynomial (poly) and its derivative
!   (pder) of degree n at x
!
      implicit double precision(a-h,o-z)
      common /jacpar/alp,bet,rv
      apb = alp + bet
      poly = 1d0
      pder = 0d0
      if (n.eq.0) return
      polylst = poly
      pderlst = pder
      poly = 5d-1*(1d0+bet)*(x-1d0) + 5d-1*(1d0+alp)*(x+1d0)
      pder = 5d-1*(2d0+apb)
      if (n.eq.1) return
      do 19 k=2,n
        a1 = 2d0*k*(k+apb)*(2d0*k+apb-2d0)
        a2 = (2d0*k+apb-1d0)*(alp**2-bet**2)
        b3 = (2d0*k+apb-2d0)
        a3 = b3*(b3+1d0)*(b3+2d0)
        a4 = 2d0*(k+alp-1d0)*(k+bet-1d0)*(2d0*k+apb)
        polyn = ((a2+a3*x)*poly-a4*polylst)/a1
        pdern = ((a2+a3*x)*pder-a4*pderlst+a3*poly)/a1
        psave = polylst
        pdsave = pderlst
        polylst = poly
        poly = polyn
        pderlst = pder
        pder = pdern
   19 continue
      polym1 = polylst
      pderm1 = pderlst
      polym2 = psave
      pderm2 = pdsave
      return
      end

    !$ Calc Jacobian for quadrilateral elements 
   
      subroutine calcJ4(N, P, alphamax, betamax, xq, yq, J4)
      implicit none
      
      integer           :: N,alphamax,betamax,a,b,i,j
      double precision  :: xq(0:alphamax,0:betamax,5), yq(0:alphamax,0:betamax,5), P(0:N) 
      double precision  :: J4(0:alphamax,0:betamax,0:N,0:N), a11, a12, a21, a22
      
      do a=0,alphamax ! Loop over elements in x
         do b=0,betamax ! Loop over elements in y
            do i=0,N     ! Loop over GLL points
               do j=0,N
                  
                  a11 = (xq(a,b,2) - xq(a,b,1))*(1d0 - P(j)) + (xq(a,b,3) - xq(a,b,4))*(1d0 + P(j))
                  a12 = (yq(a,b,2) - yq(a,b,1))*(1d0 - P(j)) + (yq(a,b,3) - yq(a,b,4))*(1d0 + P(j))
                  a21 = (xq(a,b,4) - xq(a,b,1))*(1d0 - P(i)) + (xq(a,b,3) - xq(a,b,2))*(1d0 + P(i))
                  a22 = (yq(a,b,4) - yq(a,b,1))*(1d0 - P(i)) + (yq(a,b,3) - yq(a,b,2))*(1d0 + P(i))

                !Jacobian
                  J4(a,b,i,j) = (1d0/16d0)*(a11*a22 - a21*a12)                    


               enddo
            enddo
         enddo
      enddo
      end 
      

     !$ Calc Z matrix (see Eqn 38)
      subroutine calcZ6(N, P, alphamax, betamax, xq, yq, Z6)
      implicit none
      integer          :: a,b,i,j,alphamax,betamax,N
      double precision :: Z6(0:alphamax,0:betamax,0:N,0:N,2,2), P(0:N),x,y
      double precision :: xq(0:alphamax,0:betamax,5), yq(0:alphamax,0:betamax,5)
      
     do a = 0, alphamax
         do b = 0, betamax !Loop over elements
            do i = 0, N   
               do j = 0, N  !Loop over GLL points

            Z6(a,b,i,j,1,1) = 0.25d0*((yq(a,b,4) - yq(a,b,1))*(1d0 - P(i)) + (yq(a,b,3) - yq(a,b,2))*(1d0 + P(i)))
            Z6(a,b,i,j,2,1) = -0.25d0*((yq(a,b,2) - yq(a,b,1))*(1d0 - P(j)) + (yq(a,b,3) - yq(a,b,4))*(1d0 + P(j)))
            Z6(a,b,i,j,1,2) = -0.25d0*((xq(a,b,4) - xq(a,b,1))*(1d0 - P(i)) + (xq(a,b,3) - xq(a,b,2))*(1d0 + P(i)))
            Z6(a,b,i,j,2,2) = 0.25d0*((xq(a,b,2) - xq(a,b,1))*(1d0 - P(j)) + (xq(a,b,3) - xq(a,b,4))*(1d0 + P(j)))

              enddo
            enddo
         enddo
      enddo

    endsubroutine 
      


      double precision function x(al, alphamax, be, betamax, xi, zeta, xq)
      implicit none
      integer          ::  al, be, alphamax, betamax
      double precision ::  xi, zeta, xq(0:alphamax, 0:betamax, 5),&
                           N1, N2, N3, N4
      
      !!$ Given an element (al,be) and postion within (xi,zeta), function returns
      !!$ the x-coordinate
      !!$ Vertices of quadrilateral given by xq. Note xq(1)=xq(5).

      N1 = 0.25d0*(1d0 - xi)*(1d0 - zeta)
      N2 = 0.25d0*(1d0 + xi)*(1d0 - zeta)
      N3 = 0.25d0*(1d0 + xi)*(1d0 + zeta)
      N4 = 0.25d0*(1d0 - xi)*(1d0 + zeta)

      x = N1*xq(al,be,1) + N2*xq(al,be,2) + N3*xq(al,be,3) + N4*xq(al,be,4)

      return
      end
      

     
      double precision function y(al, alphamax, be, betamax, xi, zeta, yq)
      implicit none
      integer          ::  al, be, alphamax, betamax
      double precision ::  xi, zeta, yq(0:alphamax, 0:betamax, 5),&
                           N1, N2, N3, N4
      
      !!$ Given an element (al,be) and postion within (xi,zeta), function returns
      !!$ the y-coordinate
      !!$ Vertices of quadrilateral given by xq
      
      N1 = 0.25d0*(1d0 - xi)*(1d0 - zeta)
      N2 = 0.25d0*(1d0 + xi)*(1d0 - zeta)
      N3 = 0.25d0*(1d0 + xi)*(1d0 + zeta)
      N4 = 0.25d0*(1d0 - xi)*(1d0 + zeta)

      y = N1*yq(al,be,1) + N2*yq(al,be,2) + N3*yq(al,be,3) + N4*yq(al,be,4)

      return
      end


     !$ Lagrange interpolant (not just at GLL points) 
      double precision function hbasis(i,xz,N,P,eps)
      implicit none
      integer          :: N,i,j
      double precision :: P(0:N),S,xz,eps

!!$      if(xz.lt.-(1d0+eps))then!!.or.(xz.gt.(1d0+eps)))then
!!$       !!  write(*,*)'h out of range',xz
!!$       xz = -1d0
!!$      elseif(xz.gt.(1d0+eps))then
!!$       xz = 1d0
!!$      endif
  
    
      S=1d0
      do j=0,N
         if(P(i).ne.P(j))S=S*(xz-P(j))/(P(i)-P(j))
      enddo
      hbasis=S

      return
      end
      

      
     !$ Write solution to output files
      subroutine writepoints_grid(alphamax, betamax,&
      P, N, U, J4, Un, xq, yq, time2, imax, ii, PstressL)
      implicit none
      integer           :: alphamax,betamax,N,al,be,i,j,imax,l,k,ii
      double precision  :: x, y, P(0:N), U(0:imax, 0:alphamax, 0:betamax, 0:3, 0:N, 0:N),&
                           J4(0:alphamax,0:betamax,0:N,0:N), xq(0:alphamax,0:betamax,5),&
                           yq(0:alphamax,0:betamax,5), time2, Un(0:imax, 0:alphamax, 0:betamax, 0:3, 0:N, 0:N),&
                           PstressL(0:alphamax,0:betamax,0:N,0:N,2,2)
                            
                           
    !!$ Tecplot output file syntax           
      write(99,*) 'VARIABLES = "x","y","u","v","Q","P11","P12","P22"'
      write(99,*) 'ZONE T="time:',time2,'" I =', ((betamax + 1)*(N + 1)),', J =', ((alphamax + 1)*(N + 1)),',F=POINT' 

      l = 0
   
      do al = 0, alphamax   
         do i = 0, N
             k = 0  
            do be = 0, betamax
               do j = 0, N
                 
                  
      write(99,*) x(al, alphamax, be, betamax, P(i), P(j), xq), y(al, alphamax, be, betamax, P(i), P(j), yq),&
                  U(0,al,be,1,i,j), U(0,al,be,2,i,j),U(0,al,be,0,i,j),&
                  PstressL(al,be,i,j,1,1),PstressL(al,be,i,j,1,2),PstressL(al,be,i,j,2,2)

                      
                  k = k + 1             
               enddo
             enddo
         l = l + 1
         enddo
      enddo
      end 
     
     !$ Output solutions using global variables     
      subroutine writepoints_global(alphamax, betamax, P, N, UG3, J4, imax, time2, ii)
      implicit none
      integer           :: alphamax,betamax,N,al,be,i,j,imax,l,k, Nhat, Mhat, ii
      double precision  :: x, y, P(0:N), U(0:imax, 0:alphamax, 0:betamax, 0:3, 0:N, 0:N),&
                           J4(0:alphamax,0:betamax,0:N,0:N), Len, Height,&
                           UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),time2
                           
   
       
      
      write(888,*) 'VARIABLES = "u", "v"'
      write(888,*) 'ZONE T="',time2,'" I =', ((betamax + 1)*N + 1),', J =', ((alphamax + 1)*N + 1),',F=POINT' 

      Nhat = N*(alphamax + 1)
      Mhat = N*(betamax + 1)
   
      do i = 0, Nhat
       do j = 0, Mhat
                 
        write(888,*) UG3(ii,1,i,j), UG3(ii,2,i,j)   
                 
       enddo
      enddo
         

      end 
                 


     !$ Calc gloabl version of matrix A (Eqn 36)
      subroutine calcAG4(N,alphamax,betamax,W,J4,AG4)
      implicit none
      integer           :: i,j,k,l,alphamax,betamax,al,be,N,Mhat,Nhat
      double precision  :: W(0:N),&
                           J4(0:alphamax,0:betamax,0:N,0:N),S,AL4,&
                           AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
                           comp_time
      
      Mhat = N*(betamax+1)
      Nhat = N*(alphamax+1)

      do i = 0, Nhat
      do j = 0,Mhat
      do k = 0,Nhat
      do l = 0,Mhat
      S = 0d0
      do al = 0, alphamax
         do be = 0, betamax
           
            S = S + AL4(al,be,i-al*N,j-be*N,k-al*N,l-be*N,N,alphamax,betamax,W,J4)
        
         enddo
       enddo
      AG4(i,j,k,l) = S

      enddo
      enddo
      enddo
      enddo
      return
      end 
    

      !Calc (local) matrix A (Eqn 36)
  double precision function & 
      AL4(al,be,i,j,k,l,N,alphamax,betamax,W,J4)
      implicit none
      integer          :: i,j,k,l,alphamax,betamax,al,be,N
      double precision :: W(0:N),h,&
                          J4(0:alphamax,0:betamax,0:N,0:N)

      if((i.lt.0).or.(j.lt.0).or.(k.lt.0).or.(l.lt.0).or.&
        (i.gt.N).or.(j.gt.N).or.(k.gt.N).or.(l.gt.N))then
         AL4 = 0d0
      else
         AL4 = J4(al,be,i,j)*W(i)*W(j)*h(k,i)*h(l,j)
      endif
      end 


     !$ Calc gloabl C matrix (see Eqn 37)
      subroutine &    
      calcCG5(N,W,h5,alphamax,betamax,Z6,CG5)
      implicit none
      integer           :: al,be,i,j,k,l,b,N,alphamax,betamax
      double precision  :: W(0:N),h5(0:N,0:N,0:N,0:N,2),&
                           Z6(0:alphamax,0:betamax,0:N,0:N,2,2),S,CL5,h, &
                           CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&
                           comp_time
     

      

      do i = 0,N*(alphamax+1)
      do j = 0,N*(betamax+1)
      do k = 0,N*(alphamax+1)
      do l = 0,N*(betamax+1)

      do b = 1, 2
      S = 0d0
      do al = 0,alphamax
         do be = 0,betamax

            S = S + CL5(al,be,i-al*N,j-be*N,k-al*N,l-be*N,&
                        b,N,W,h5,alphamax,betamax,Z6)
         enddo
        enddo
         CG5(i,j,k,l,b) = S
  
        

      enddo
      enddo
      enddo
      enddo
      enddo

      end 

      !$ Calc (local) C matrix (see Eqn 37)
       double precision function & 
      C5(al,be,i,j,k,l,b,N,W,h5,alphamax,betamax,Z6)
      implicit none
      integer           :: al,be,i,j,k,l,b,a,N,alphamax,betamax
      double precision  :: W(0:N),h5(0:N,0:N,0:N,0:N,2),&
                           Z6(0:alphamax,0:betamax,0:N,0:N,2,2),S    
      S = 0d0
      do a = 1, 2
          S = S + Z6(al,be,i,j,a,b)*W(i)*W(j)*h5(k,i,l,j,a)
      enddo
      C5 = S
      return
      end

      double precision function & 
      CL5(al,be,i,j,k,l,b,N,W,h5,alphamax,betamax,Z6)
      implicit none
      integer           :: al,be,i,j,k,l,b,N,alphamax,betamax
      double precision  :: W(0:N),h5(0:N,0:N,0:N,0:N,2), &
                           Z6(0:alphamax,0:betamax,0:N,0:N,2,2),C5    

      if((i.lt.0).or.(j.lt.0).or.(k.lt.0).or.(l.lt.0).or. &
        (i.gt.N).or.(j.gt.N).or.(k.gt.N).or.(l.gt.N))then
         CL5 = 0d0
      else
         CL5 = C5(al,be,i,j,k,l,b,N,W,h5,alphamax,betamax,Z6)
      endif
      return
      end
   
  !$ Pre-calculate iCG5, an array that takes advantage of sparsity of
  !$ matrix CG5, for efficiency purposes
      subroutine precalciCG5(N,alphamax,betamax,CG5,nnmax,mmax,eps)
      implicit none
      integer            :: i,j,k,l,p,q,nn,mm,N,alphamax,betamax,nnmax,mmax
      double precision   :: &
                            CG5(0:N*(alphamax+1), 0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&
                            x1,x2,y1,y2, eps, comp_time
      
    nn = 1
    mm = 1
    mmax = 0

      do l=0,N*(betamax+1)
      do k=0,N*(alphamax+1)
      do j=0,N*(betamax+1)
      do i=0,N*(alphamax+1)

   

         mm = 4
         do q = 0,N*(betamax+1)
         do p = 0,N*(alphamax+1)
            x1 = CG5(p,q,i,j,1)
            x2 = CG5(p,q,i,j,2)
            y1 = CG5(p,q,k,l,1)
            y2 = CG5(p,q,k,l,2)

           if((x1.ne.0d0.or.x2.ne.0d0).and.(y1.ne.0d0.or.y2.ne.0d0)) then

             mm = mm + 1

            endif
         enddo
         enddo
         if(mm.gt.mmax)mmax=mm
         if(mm.gt.4)nn=nn+1
      enddo
      enddo
      enddo
      enddo

        if(mm.gt.4)then
            nnmax=nn
        else
            nnmax=nn-1
        endif
      end

  !$ Calculate iCG5, an array that takes advantage of sparsity of
  !$ matrix CG5, for efficiency purposes

      subroutine calciCG5(N,alphamax,betamax,CG5,nnmax,mmax,NiCG5,MiCG5,iCG5, eps)
      implicit none
      integer           :: i,j,k,l,p,q,nn,mm,N,alphamax,betamax,nnmax
      integer           :: NiCG5,MiCG5
      integer           :: iCG5(NiCG5,MiCG5,2), mmax(NiCG5)
      double precision  :: CG5(0:N*(alphamax+1),0:N*(betamax+1),&
                           0:N*(alphamax+1),0:N*(betamax+1),2),x1,x2,y1,y2, eps,&
                           comp_time
      
     !!$ determines which components of CG5 are non-zero, allowing quicker construction of TG matrices.
     !!$ See subroutine calcTGG      

    mmax=0
    iCG5=0
      nn=1


      do l=0,N*(betamax+1)
      do k=0,N*(alphamax+1)
      do j=0,N*(betamax+1)
      do i=0,N*(alphamax+1)
    

        iCG5(nn,1,1)=i
        iCG5(nn,2,1)=j 
        iCG5(nn,3,1)=k 
        iCG5(nn,4,1)=l    

       
        mm=4
         do q=0,N*(betamax+1)
         do p=0,N*(alphamax+1)
            x1=CG5(p,q,i,j,1)
            x2=CG5(p,q,i,j,2)
            y1=CG5(p,q,k,l,1)
            y2=CG5(p,q,k,l,2)

              if((x1.ne.0d0.or.x2.ne.0d0).and.(y1.ne.0d0.or.y2.ne.0d0))then
  
              mm=mm+1
               iCG5(nn,mm,1)=p                                                              
               iCG5(nn,mm,2)=q
         
              endif
         enddo
         enddo
         mmax(nn)=mm
         if(mm.gt.4)nn=nn+1
      enddo
      enddo
      enddo
      enddo

        if(mm.gt.4)then
            nnmax=nn
        else
            nnmax=nn-1
        endif
      endsubroutine      

      !$ Calculate tensors for use in the stiffness matrix

      subroutine &
      calcTGG(N,alphamax,betamax,nnmax,mmax,NiCG5,MiCG5,iCG5,AG4,CG5,TG6,TG7,TG8,GammaN,GammaP)
      implicit none
      integer           :: i,j,k,l,a,b,c,d,p,q,N,alphamax,betamax,nn,mm
      integer           :: nnmax
    integer, intent(in) :: NiCG5,MiCG5
      integer           :: iCG5(NiCG5,MiCG5,2),mmax(NiCG5)     
      double precision   &
      
      CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),h,&
      AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
      TG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      TG7(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      TG8(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&    
      S1,S2,S3,S4,&
      GammaN(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      GammaP(0:N*(alphamax+1),0:N*(betamax+1),2,2,2,2)
      

      do nn = 1, nnmax

         i = iCG5(nn,1,1)
         j = iCG5(nn,2,1) 
         k = iCG5(nn,3,1) 
         l = iCG5(nn,4,1) 

         do a = 1,2
         do b = 1,2
   
            S1 = 0d0
            S2 = 0d0
            S3 = 0d0 
            S4 = 0d0  
 
          if(mmax(nn).gt.4)then

           do mm = 5, mmax(nn)
               p = iCG5(nn,mm,1)                                                              
               q = iCG5(nn,mm,2)

             
              do c = 1, 2
                     S1 = S1 + GammaN(p,q,a,b)*CG5(p,q,i,j,c)/AG4(p,q,p,q)*CG5(p,q,k,l,c)
                     S2 = S2 + GammaN(p,q,a,c)*CG5(p,q,i,j,b)/AG4(p,q,p,q)*CG5(p,q,k,l,c)

                do d = 1, 2

                     S3 = S3 + GammaP(p,q,a,b,c,d)*CG5(p,q,i,j,d)/AG4(p,q,p,q)*CG5(p,q,k,l,c)

                enddo
 
             enddo                           
            
                   S4 = S4 + CG5(p,q,i,j,a)/AG4(p,q,p,q)*CG5(p,q,k,l,b)
            enddo

          endif
          ! Each of these tensors corresponds to the different terms in the 
          ! stiffness matrix (see Eqn 44)
            TG6(i,j,k,l,a,b) = S1 + S2 + S3 ! Middle terms on RHS of Eqn 44
            TG7(i,j,k,l,a,b) = AG4(i,j,k,l)*h(a,b) !This is 1st term on RHS of Eqn 44
            TG8(i,j,k,l,a,b) = S4 ! This is last term on RHS of Eqn 44

           enddo
         enddo         
      enddo

      endsubroutine

!!$ Combines TG matrices before insertion into final global stiffness matrix
  subroutine  adjHG6(N,alphamax,betamax,c2,dt,TG6,TG7,TG8,HG6)

      implicit none
      integer            i,j,k,l,a,b,N,alphamax,betamax
      double precision   &
      HG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      TG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      TG7(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      TG8(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      c2,dt          

      do i = 0,N*(alphamax+1)
      do j = 0,N*(betamax+1)
      do k = 0,N*(alphamax+1)
      do l = 0,N*(betamax+1)
         do a=1,2
         do b=1,2
          !Tensor components of stiffness matrix (see Eqn 44)    
               HG6(i,j,k,l,a,b) = dt*TG6(i,j,k,l,a,b) + TG7(i,j,k,l,a,b) + &         !!!
                                  dt*(c2*dt)*TG8(i,j,k,l,a,b)   !nu is treated as zero-probably insignificant
          
          enddo
         enddo
      enddo
      enddo
      enddo
      enddo
      endsubroutine

     !!$ Assemble final gloabl stiffness matrix
      subroutine calcglobalmatrix(N,alphamax,betamax,Gmax,HG6X,Mglobal,Nonz)
      implicit none
      integer            i,j,k,l,a,b,N,p,q,alphamax,betamax,Nhat,Mhat,Gmax,Nonz
      double precision   &
      HG6X(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
      MGlobal(Gmax,Gmax),scal1,scal2

      !!$ Assemble global matrix

      Nonz = 0 

      Nhat = N*(alphamax+1)-1
      Mhat = N*(betamax+1)-1
      do a = 1, 2
      do j = 1, Mhat
      do i = 1, Nhat
         do b=1,2
         do l=1, Mhat      
         do k=1, Nhat
        
            p=i+(j-1)*Nhat+(a-1)*Nhat*Mhat
            q=k+(l-1)*Nhat+(b-1)*Nhat*Mhat
          !  MGlobal(p,q)=(HG6X(i,j,k,l,a,b)+HG6X(k,l,i,j,b,a))/2d0/dsqrt(HG6X(k,l,k,l,b,b)*HG6X(i,j,i,j,a,a))
          
            scal1=abs(HG6X(k,l,k,l,b,b))
            scal2=abs(HG6X(i,j,i,j,a,a))

            MGlobal(p,q)=(HG6X(i,j,k,l,a,b))/dsqrt(scal1*scal2)!!dsqrt(HG6X(k,l,k,l,b,b)*HG6X(i,j,i,j,a,a))
            !Scale by diagonal elements (Jacobi preconditioning)

        if(abs(MGlobal(p,q)).gt.1d-15)then
          Nonz = Nonz + 1
         endif
         
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo


endsubroutine

 !Recast stiffness matrix into 1D arrays for use in PARDISO
 subroutine calcStorageArrays(Gmax,Nonz,Mglobal,values,columns,rowIndex)
 implicit none


 integer :: Nonz, first, i, Gmax, p, q, rowIndex(Gmax+1),columns(Nonz)
 double precision :: values(Nonz),MGlobal(Gmax,Gmax)
 

 i = 1

 do p = 1, Gmax

   first = 0

  do q = 1, Gmax

   if(abs(MGlobal(p,q)).gt.1d-15)then  

    values(i) = MGlobal(p,q) 
   columns(i) = q

    if(first.eq.0)then
      rowIndex(p) = i
            first = 1
    endif
 
    i = i + 1

   endif

  enddo

enddo 

 rowIndex(Gmax+1) = Nonz + 1


endsubroutine

! Subroutine which calls external PARDISO solvers
subroutine Pardiso_Fact(pt, Nonz, Gmax, values,columns,rowIndex)
implicit none

 integer*8  pt(64)
 integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
 integer :: solver, idum, Gmax, Nonz
 integer :: iparm(64), rowIndex(Gmax+1),columns(Nonz)
 
 double precision :: values(Nonz)
 double precision :: dparm(64), ddum
    
      pt = 0 !internal pointer - MUST NOT BE CHANGED FROM HERE ON IN
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
  maxfct = 1 ! max number of matrices (i think)

   mtype = 11 !real unsymmetric matrix 
 ! mtype = 2 for real symmetric pos def matrix

   phase = 12 !perform analysis and numerical factorisation

  msglvl = 0 !print no (=0) statistical information on screen 

  iparm(1) = 0 !use default solver settings
  iparm(3) = 1

  idum = 0 !dummy variables which aren't used for this particular phase=12
  ddum = 0

   call pardiso (pt, maxfct, mnum, mtype, phase, Gmax, values, rowIndex, columns,& 
                 idum, nrhs, iparm, msglvl, ddum, ddum, error) 


   write(*,*) 'error from pardiso',error

 endsubroutine

 !!$ Lapack (NAG) routine for Cholesky factorisation
subroutine CholeskyFact(M,G2,ipiv)
      implicit none
      integer           M,IFAIL,ipiv(M),info
      double precision  G2(M,M)

      IFAIL=0
    !!  call F07FDF('U',M,G2,M,IFAIL)!cholesky
     !!   call DPOTRF('U',M,G2,M,info)
     !!$ Lapack routine for Cholesky factorisation

      call DGETRF(M, M, G2, M, ipiv, info)
     !!$ Lapack routine for general LU factorisation   

         write(*,*)'info is',info   
      endsubroutine

   !$$ Set velocity boundary conditions on wall
   subroutine calcB3(N,alphamax, betamax, B3_Left, B3_Right, B3_Top, B3_Bot)

   implicit none
      integer           :: N, alphamax,betamax
      double precision  :: B3_Left(1:N*(betamax + 1) - 1, 3),&
                           B3_Right(1:N*(betamax + 1) - 1, 3),&
                           B3_Top(0:N*(alphamax + 1), 3),&
                           B3_Bot(0:N*(alphamax + 1), 3)
   
!!$ BOUNDARY CONDITIONS ON WALL (NO SLIP)                      
                          
          B3_Left = 0d0
         B3_Right = 0d0
           B3_Top = 0d0
           B3_Bot = 0d0
     

      return
      endsubroutine

   !$ Calc velocity Un at previous time and particle position x^n
   !$ (see semi Lagrangian scheme as described in Lind and Phillips 2011)
       subroutine calcUn(alphamax, betamax, N, xq, yq, P,dt,it,&
                         U,Un,Qn,flag,imax,speed, Len, Height, PstressL, Pstressn, PstressnL,eps,time)
       implicit none
       integer          :: al,be,i,j,ii,jj,alphamax,betamax,N,alpha,beta,&
                           a,it,flag,imax, k, l,b,time

      double precision  :: x,y,dt,Height,Len,delta,&
                           U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                           Un(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                           RR(0:alphamax,0:betamax,0:N,0:N),&
                           S(0:3),P(0:N),xx,yy,xi,zeta,&
                           Qn(0:N*(alphamax+1),0:N*(betamax+1)),&
                           speed, xq(0:alphamax, 0:betamax,5),&
                           yq(0:alphamax, 0:betamax,5), eps,&
                           Pxy(2,2),&
                           PstressL(0:alphamax,0:betamax,0:N,0:N,2,2),&
                           PstressnL(0:alphamax,0:betamax,0:N,0:N,2,2),&
                           Pstressn(0:N*(alphamax+1),0:N*(betamax+1),2,2)
      


     RR = 0d0
 
      do al = 0, alphamax
      do be = 0, betamax
         do i = 0, N
         do j = 0, N
          
!!$ Calculate previous postion of fluid particle      
!!$ change .eq. to .ge. to use implicit Euler for x^n            
            if(it.eq.0)then 
               xx = x(al, alphamax, be, betamax, P(i), P(j), xq) - dt*U(it,al,be,1,i,j)
               yy = y(al, alphamax, be, betamax, P(i), P(j), yq) - dt*(U(it,al,be,2,i,j))!+1d-1)
            else
!*     midpoint rule for x^n
               xx = x(al, alphamax, be, betamax, P(i), P(j), xq) - 0.5d0*dt*(U(it,al,be,1,i,j) + Un(it-1,al,be,1,i,j))
               yy = y(al, alphamax, be, betamax, P(i), P(j), yq) - 0.5d0*dt*(U(it,al,be,2,i,j) + Un(it-1,al,be,2,i,j))!+2d-1)
            endif 

 
!!$ If the previous point outside the domain, the set it to be on the domain

        if((xx.lt.eps).or.(xx.gt.Len).or.(yy.lt.eps).or.(yy.gt.Height))then

              if(xx.le.eps)then
              
               xx = 0d0              
           
              endif

              if(xx.ge.Len)then
            
               xx = Len 
       
              endif

               if(yy.le.eps)then
              
               yy = 0d0             
                
               endif

                if(yy.ge.Height)then
               
                yy = Height             

                endif

           endif

     !!$ For a given (x,y),calculate the containing element (alpha,beta) and position (xi,zeta) within that element
             call x2xi(alphamax,betamax,alpha,beta, xx, yy, xq, yq, xi, zeta, eps)
            
     !!$ Given this postion - calculate the (local) velocity Un
             call calcUxy(U,alpha,beta,xi,zeta,0,N,P,alphamax,betamax,S,imax,0d0)!also includes calculation for Qn
           
             do a = 0,3
               Un(it,al,be,a,i,j) = S(a)
             enddo

             Pxy = 0d0

         !!$ Calculate the value of the polymeric stress at x^n 
             call calcPxy(PstressL,alpha,beta,xi,zeta,N,P,alphamax,betamax,Pxy,0d0)

           do a = 1, 2
            do b = 1, 2
               PstressnL(al,be,i,j,a,b) = Pxy(a,b)
            enddo
           enddo

     !!$ Check error RR between subsequent iterations for Un.
             if(it.gt.0)RR(al,be,i,j) = abs(Un(it,al,be,2,i,j) - Un(it-1,al,be,2,i,j))&
                                      + abs(Un(it,al,be,1,i,j) - Un(it-1,al,be,1,i,j))

 
          enddo
         enddo
      enddo
      enddo

 
      !$$ Update global density and polymeric stress for use in 
      !$$ subsequent iterations
         do i = 0, N*(alphamax + 1)
           do j = 0, N*(betamax + 1) 
               
               al = (i-1)/N !Unfortunately this does not work for  N=1
               be = (j-1)/N
                k = mod(i+al,N+1)
                l = mod(j+be,N+1)
             
              Qn(i,j) = Un(it,al,be,0,k,l)

              do a = 1, 2
               do b = 1, 2
                 Pstressn(i,j,a,b) =  PstressnL(al,be,k,l,a,b)
               enddo
              enddo
             
            enddo
          enddo

 !$ Check if convergence criterion for iteration satisfied

     if(it.gt.0)then
         delta = maxval(RR) !abs(Un(it,0,0,2,2,0)-Un(it-1,0,0,2,2,0))
         write(*,*) 'maxval RR',delta,maxloc(RR)
         if(delta.lt.1d-5.or.flag.eq.1)flag = flag-1 !convergence ends iteration  
      endif      
      
      return
      endsubroutine

     !$$ Given an (x,y), find the surrounding element and local coordinates (xi,zeta)
      subroutine x2xi(alphamax,betamax,alpha,beta, xx, yy, xq, yq, xi, zeta, eps)
      implicit none

      integer           :: alphamax, betamax, al, be, alpha, beta, alpha1, beta1,i

      double precision  :: xi, zeta, xx, yy, xq(0:alphamax, 0:betamax,5),& 
                           yq(0:alphamax, 0:betamax,5), A, Area, vx(4), vy(4), eps


      alpha1 = -1
        beta1 = -1

!!$ Finds the element inwhich a point (x,y) resides

alphaloop: do al = 0, alphamax !Loop over all elements
            do be = 0, betamax
      iloop: do i = 1, 4 !Loop over vertices
  
               A = Area(xq(al,be,i+1),xq(al,be,i),yq(al,be,i+1),yq(al,be,i),xx,yy) 
 !Determines if point is on left (A>0) of an edge of element (see Computational Geometry in C. J. Roukre)

                if(A.lt.-1d-14)exit iloop
                if(i.eq.4)then
                 alpha1 = al            !!$ If point is left of all sides of element, then
                  beta1 = be            !!$ point is within that element
                  exit alphaloop
                 endif

              enddo iloop
            enddo
           enddo alphaloop
   
          if(alpha1.lt.0.or.beta1.lt.0)then
            write(*,*) 'could not find point within elements',xx,yy
              stop
           endif

            alpha = alpha1  !!$ Element labels which contain point (x,y)
             beta = beta1

        !! Found elements alpha, beta, now find xi, zeta within element
          
          do i = 1, 4
            vx(i) = xq(alpha1,beta1,i)
            vy(i) = yq(alpha1,beta1,i)
           enddo

        !!$ Inverts the transfinite map to find (xi, zeta), given (x, y)
       call quad2squ(al, be, xx, yy, vx, vy, xi, zeta, eps)


    !!$ If cant find a (xi,zeta) - error
      if(xi.lt.-(1d0+eps).or.xi.gt.(1d0+eps))then
         write(*,*)'error xi=',xi,'al=',alpha,'be=',beta
       stop
      endif
      if(zeta.lt.-(1d0+eps).or.zeta.gt.(1d0+eps))then
      write(*,*)'error zeta=',zeta,'al=',alpha,'be=',beta          
       stop      
      endif

  return
  endsubroutine

     !$$ Function used in determining the element in which 
     !$$ a particle resides.
      double precision function Area(x2,x1,y2,y1,xx,yy)
      implicit none
      double precision :: x2,x1,y2,y1,xx,yy

      Area = (x2 - x1)*(yy - y1) - (xx - x1)*(y2 - y1)

      return
      end

     !$$ Inverts transfinite map to determine a (xi,zeta), given an (x,y)
      subroutine quad2squ(al, be, xx, yy, xxq, yyq, xi, zeta, eps)
      implicit none

      integer          :: al, be
      double precision :: a,b,c,d,e,f,g,h,aa,bb,cc,dd,ee,ff,&
                          x_pos, x_neg, xi, zeta, xxq(4), yyq(4),&
                          xx, yy, y_pos, y_neg, eps


        a = (xxq(1) + xxq(2) + xxq(3) + xxq(4))/4d0
        b = (-xxq(1) + xxq(2) + xxq(3) -xxq(4))/4d0
        c = (-xxq(1) - xxq(2) + xxq(3) + xxq(4))/4d0
        d = (xxq(1) - xxq(2) + xxq(3) - xxq(4))/4d0

        e = (yyq(1) + yyq(2) + yyq(3) + yyq(4))/4d0
        f = (-yyq(1) + yyq(2) + yyq(3) - yyq(4))/4d0
        g = (-yyq(1) - yyq(2) + yyq(3) + yyq(4))/4d0
        h = (yyq(1) - yyq(2) + yyq(3) - yyq(4))/4d0

       aa = -f*d + h*b
       bb = (yy - e)*d + (a - xx)*h - c*f + g*b 
       cc = (yy - e)*c + (a - xx)*g

       dd = h*c - g*d
       ee = (yy - e)*d + (a - xx)*h + c*f - g*b
       ff = (yy - e)*b + (a - xx)*f
     
     !!$ solves quadratic equation seen in the inverse map using
     !!$ simple formula
      call quadratic_eqn_solver(aa, bb, cc, x_pos, x_neg, eps)
      call quadratic_eqn_solver(dd, ee, ff, y_pos, y_neg, eps)

    
     !!$ Of the solutions, determine which "make sense"
     !!$ (xi,zeta) should both both be between [-1,1]

      call root_decider(x_pos, x_neg, xi, eps)
      call root_decider(y_pos, y_neg, zeta, eps)
     
     !!$ This does always work as the map IS bijective over each
     !!$ element.

      return
      endsubroutine

     !!$ Solves quadratic eqn that appears in inverse transfinite map
      subroutine quadratic_eqn_solver(a1, b1, c1, s_pos, s_neg, eps)
      implicit none

      double precision :: a1, b1, c1, s_pos, s_neg, determinant, eps, abs_a1,abs_b1

              

      determinant = b1**2d0 - 4d0*a1*c1

           abs_a1 = dabs(a1)
           abs_b1 = dabs(b1)

      if((abs_a1.le.eps))then

      s_pos = -c1/b1
      s_neg = 2d0

      elseif(determinant.lt.0d0)then

      write(*,*) 'xi, zeta, complex roots'
      stop

      else

      s_pos = (-b1 + dsqrt(determinant))/(2d0*a1)
      s_neg = (-b1 - dsqrt(determinant))/(2d0*a1)

     
      endif
      return

      

      endsubroutine

     !$$ Decides which solutions of quadratic equations 
     !$$ are the "correct" ones.
      subroutine root_decider(s_pos, s_neg, soln, eps)
      implicit none

      double precision :: s_pos, s_neg, soln, eps

     

       if((((s_pos.gt.(1d0+eps))).or.(s_pos.lt.(-(1d0+eps)))).and.(((s_neg.gt.(1d0+eps))).or.(s_neg.lt.(-(1d0+eps)))))then
    
        write(*,*) 'xi, zeta out of limits',s_pos,s_neg
        stop

       elseif(((s_neg.gt.(1d0+eps))).or.(s_neg.lt.(-(1d0+eps))))then
        
       soln = s_pos

       elseif(((s_pos.gt.(1d0+eps))).or.(s_pos.lt.(-(1d0+eps))))then

       soln = s_neg

       else

       write(*,*) 'multiple xi, zeta within parent element'
       write(*,*) s_pos,s_neg
       stop

       endif

       return
       endsubroutine
     
     !!$ Determines interpolated velocity from sum of basis functions 
     !!$ (See Eqn 32)
      subroutine calcUxy(U,al,be,xi,zeta,it,N,P,alphamax,betamax,Uxy,imax, eps)
      implicit none
      integer           :: it,al,be,N,a,i,j,alphamax,betamax,imax
      double precision  :: S,P(0:N),hbasis,xi,zeta,Uxy(0:3), eps,&
                           U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N)

      !!$ determine interpolated velocity from sum of basis functions 

      do a = 0,3
         S = 0d0
         do i = 0,N
         do j = 0,N
            S = S + U(it,al,be,a,i,j)*&
                    hbasis(i,xi,N,P,eps)*hbasis(j,zeta,N,P,eps)
         enddo
         enddo
         Uxy(a)=S
       enddo

       endsubroutine

     !!$ Determine interpolated stress from sum of basis functions 
     !!$ (See Eqn 33)
       subroutine calcPxy(PstressL,al,be,xi,zeta,N,P,alphamax,betamax,Pxy,eps)
       implicit none
       integer           :: it,al,be,N,a,b,i,j,alphamax,betamax
       double precision  :: S,P(0:N),hbasis,xi,zeta,Pxy(2,2),eps,&
                            PstressL(0:alphamax,0:betamax,0:N,0:N,2,2)

     

      do a = 1,2
       do b = 1,2
         S = 0d0
         do i = 0,N
         do j = 0,N
            S = S + PstressL(al,be,i,j,a,b)*&
                    hbasis(i,xi,N,P,eps)*hbasis(j,zeta,N,P,eps)
         enddo
         enddo
         Pxy(a,b)=S
       enddo
       enddo

       
       endsubroutine
   

      !$$ Construct global matrix containing Un for use
      !$$ in solution of gloabl linear system
      subroutine &     
      calcfG3(N,it,alphamax,betamax,Un,W,J4,fG3,imax)
      implicit none
      integer           :: al,be,k,l,b,N,alphamax,betamax,it,imax
      double precision  :: W(0:N),J4(0:alphamax,0:betamax,0:N,0:N),&
                           Un(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),fL3,S,&
                           fG3(0:N*(alphamax+1),0:N*(betamax+1),3)
     
      do k = 0, N*(alphamax + 1)
      do l = 0, N*(betamax + 1)
         do b = 1, 3 !3rd component is temperature               
            S = 0d0
            
            do al = 0, alphamax
            do be = 0, betamax
            
              S = S + fL3(al, be, k - al*N, l - be*N, b,N, it, alphamax, betamax, W, Un, J4, imax)
     
            enddo
            enddo
            fG3(k,l,b) = S
             
          

         enddo
      enddo
      enddo
      end

       double precision function & 
       fL3(al,be,k,l,b,N,it,alphamax,betamax,W,Un,J4,imax)
       implicit none
       integer           :: al,be,k,l,b,N,alphamax,betamax,it,imax
       double precision  :: W(0:N),J4(0:alphamax,0:betamax,0:N,0:N),&
                            Un(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N)
       
        if((k.lt.0).or.(l.lt.0).or. &
          (k.gt.N).or.(l.gt.N))then         
         fL3 = 0d0
      else
         fL3 = J4(al,be,k,l)*W(k)*W(l)*Un(it,al,be,b,k,l)
      endif
      return
      end
          

           !$$ Global vector of knowns (Eqn 47) on RHS of linear system (Eqn 43) 
              subroutine &    
              calcVG3(imax,N,dt,c2,alphamax,betamax,B3_Left, B3_Right, B3_Top, B3_Bot,&
                      CG5,Qn,fG3,HG6,UG3,YG3,ZG3,VG3,Pstressn,Bn,lambda1)

              !!$ Vector of knowns on RHS of Mu=V (Eqn 43)

              implicit none
              integer          :: imax,k,l,a,b,N,alphamax,betamax,i,j,Nhat,Mhat
                           
              double precision :: S1, S2, dt,c2,C1,C3,&
                                  B3_Left(1:N*(betamax + 1) - 1, 3),&
                                  B3_Right(1:N*(betamax + 1) - 1, 3),&
                                  B3_Top(0:N*(alphamax + 1), 3),&
                                  B3_Bot(0:N*(alphamax + 1), 3),&
                                  fG3(0:N*(alphamax+1),0:N*(betamax+1),3),&   !test fG3 has infact 3 components (not 2) 
                                  CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&
                                  VG3(N*(alphamax+1)-1,N*(betamax+1)-1,2),&
                                  YG3(0:N*(alphamax+1),0:N*(betamax+1),2),&! grad Q.T term
                                  ZG3(0:N*(alphamax+1),0:N*(betamax+1),2),&
                                  Qn(0:N*(alphamax+1),0:N*(betamax+1)),&
                                  HG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                                  UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                                  Pstressn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                                  Bn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                                  lambda1(0:N*(alphamax+1),0:N*(betamax+1))
     

          
        Nhat = N*(alphamax + 1)
        Mhat = N*(betamax + 1)
       

!$$ Known velocities on the boundary

       do l = 0, Nhat
 
          do i = 0, imax 
            do b = 1, 3 
                          
                    UG3(i,b,l,0) = B3_Bot(l,b)
                 UG3(i,b,l,Mhat) = B3_Top(l,b)

            enddo 
                
            enddo
      enddo
      
       do l = 1, Mhat - 1
          do i = 0, imax 
            do b = 1, 3 
                            
                  UG3(i,b,0,l) = B3_Left(l,b)
               UG3(i,b,Nhat,l) = B3_Right(l,b)

            enddo 
                 
             enddo
      enddo

      do k = 1, N*(alphamax + 1) - 1
        do l = 1, N*(betamax + 1) - 1


          do b = 1, 2

             S1 = 0d0
             S2 = 0d0

            do j = 0, Mhat
               do i = 0, Nhat
                   
                 S1 = S1 + c2*dt*Qn(i,j)*CG5(i,j,k,l,b)

                 C1 = (lambda1(i,j)/(1d0 + lambda1(i,j)/dt))
            
               do a = 1, 2
                 S2 = S2 + C1*Pstressn(i,j,a,b)*CG5(i,j,k,l,a)
               enddo
                
               enddo
             enddo
 
    !$ Insert known boundary conditions here
               do j = 1, Mhat - 1  
                  do a = 1, 2                            
                   S1 = S1 - B3_Left(j,a)*HG6(0,j,k,l,a,b)  
                   S1 = S1 - B3_Right(j,a)*HG6(Nhat,j,k,l,a,b)
               enddo
             enddo
   
             do i = 0, Nhat
                  do a = 1, 2
                    S1 = S1 - B3_Bot(i,a)*HG6(i,0,k,l,a,b)        
                    S1 = S1 - B3_Top(i,a)*HG6(i,Mhat,k,l,a,b)  
                   enddo
               enddo
           
          VG3(k,l,b) = S1 - S2 + fG3(k,l,b) + dt*YG3(k,l,b)

            enddo
       enddo
       enddo
       endsubroutine


      !$$ Calc tensor which appears on RHS of global linear system (Eqn 43)
      !$$ See last term of Eqn 47.
      subroutine & !Evaluates grad Q.T
      calcYG3(N,imax,alphamax,betamax,UG3,CG5,TG4,Pstress,Extra_Stress,YG3)
      implicit none
      integer          ::  i,j,a,b,c,p,q,N,alphamax,betamax,imax
      double precision :: &
                          CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),h,&
                          TG4(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                          QG3(0:N*(alphamax+1),0:N*(betamax+1),2),&
                          YG3(0:N*(alphamax+1),0:N*(betamax+1),2),&
                          UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                          S1, Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                          Extra_Stress(0:N*(alphamax+1),0:N*(betamax+1),2,2)

   
      do i = 0, N*(alphamax+1)
      do j = 0, N*(betamax+1)
         
        do a = 1, 2

            S1 = 0d0
            do p = 0,N*(alphamax+1)
            do q = 0,N*(betamax+1)

                  S1 = S1 + UG3(0,0,p,q)*CG5(i,j,p,q,a)

            enddo
            enddo

            QG3(i,j,a) = S1

         enddo
      enddo
      enddo

!    YG3 represents the (grad q).T term

      do i = 0,N*(alphamax+1)
      do j = 0,N*(betamax+1)
 
         do a = 1, 2

            S1 = 0d0

            do b = 1, 2

              Extra_Stress(i,j,a,b) = TG4(i,j,a,b) + Pstress(i,j,a,b)

              !S1 = S1 + QG3(i,j,b)*(TG4(i,j,a,b) + Pstress(i,j,a,b))
              S1 = S1 + QG3(i,j,b)*(Extra_Stress(i,j,a,b))
            enddo

            

            if(i.eq.0.or.i.eq.N*(alphamax+1).or.j.eq.0.or.j.eq.N*(betamax + 1))then
       
                  YG3(i,j,a) = 0d0
            else
                  YG3(i,j,a) = S1
            endif

           

         enddo
      enddo
      enddo      
      end

    !$$ Assemble global vector on RHS on Eqn 43 for use in solvers
      subroutine calcglobalvector(N,alphamax,betamax,Gmax,VG3,HG6,VGlobal)
      implicit none
      integer           :: i,k,l,b,N,alphamax,betamax,Nhat,Mhat,Gmax
      double precision  :: &
                          VG3(N*(alphamax+1)-1,N*(betamax+1)-1,2),&
                          HG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                          VGlobal(Gmax),scal1
     
      Nhat = N*(alphamax+1)-1
      Mhat = N*(betamax+1)-1
      
      do b = 1, 2
         do k = 1, Nhat
            do l = 1, Mhat

               i = k+(l-1)*Nhat+(b-1)*Nhat*Mhat

                scal1=abs(HG6(k,l,k,l,b,b))
               !$$ Scaling factor

               VGlobal(i) = VG3(k,l,b)/dsqrt(scal1)!!dsqrt(HG6(k,l,k,l,b,b))
            
           enddo
         enddo
      enddo

      endsubroutine

     !$$ Solve Linear system (Eqn 43) using PARDISO solver

      subroutine Pardiso_Solve(pt, Nonz, Gmax, values,columns,rowIndex,VGlobal)
      implicit none

       integer*8  pt(64)
       integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
       integer :: solver, idum, Gmax, Nonz
       integer :: iparm(64), rowIndex(Gmax+1),columns(Nonz)
 
      double precision :: values(Nonz),VGlobal(Gmax),soln(Gmax)
      double precision :: dparm(64), ddum
    
      nrhs = 1 ! number of rhs vectors
      mnum = 1 ! number of matrices
    maxfct = 1 ! max number of matrices (i think)

     mtype = 11 !real unsymmetric matrix 
 ! mtype = 2 for real symmetric pos def matrix

     phase = 33 !solves equation

    msglvl = 0 !print no (=0) statistical information on screen 

  iparm(1) = 0 !use default solver settings
  iparm(3) = 1
  iparm(6) = 1 !write solution onto vglobal

  idum = 0 !dummy variables which aren't used for this particular phase=12
  ddum = 0


   call pardiso (pt, maxfct, mnum, mtype, phase, Gmax, values, rowIndex, columns,& 
                 idum, nrhs, iparm, msglvl, VGlobal, soln, error) 

     write(*,*) 'error from pardisosolve',error

  endsubroutine

    !$$ Alternatively, solve factorised linear system
    !$$ using Lapack/NAG rountines

      subroutine solveSym4U(M,G2,V1,ipiv)
      implicit none
      integer          :: M,IFAIL,ipiv(M),info
      double precision :: G2(M,M),V1(M)

      IFAIL=0

     !! call F07FEF('U',M,1,G2,M,V1,M,IFAIL)
    !!   call DPOTRS('U',M,1,G2,M,V1,M,IFAIL)

     !! Lapack routine for solving linear system with 
     !! with matrix having been previously factorised by Chloseky
      
     call dgetrs('N',M,1,G2,M,ipiv,V1,M,info)

     !! Lapack routine for solving general LU factorised linear system

     write(*,*) 'info in solvesym4u is',info

      endsubroutine     


       !$$ Update global velocity with solution from linear system

       subroutine & 
       updateUG3(N,Gmax,it,dt,alphamax,betamax,solution,Qn,&
                 CG5,AG4,UG3,HG6,imax,F)


     !!$ Set solutions as UG3 (global velocity)

      implicit none
      integer           :: N, Gmax, alphamax, betamax, it, imax, F
      integer           :: a, i, j, k, l, Nhat, Mhat, alpha   
      double precision  :: &
                           solution(Gmax),S,dt,&
                           AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
                           HG6(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                           CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&
                           UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                           Qn(0:N*(alphamax+1),0:N*(betamax+1)),scal2

      Nhat = N*(alphamax + 1)
      Mhat = N*(betamax + 1)
      
      alpha = 0

      do a = 1, 2
         do j = 1, Mhat - 1   
            do i = 1, Nhat - 1

               alpha = alpha + 1

               scal2=abs(HG6(i,j,i,j,a,a))

               UG3(it,a,i,j) = solution(alpha)/dsqrt(scal2)!!sqrt(HG6(i,j,i,j,a,a))!updates velocity

            enddo
         enddo
      enddo 
    !$$ This rountine is also called within the Lagrangian iteration
    if(F.eq.0) then !If out of iteratation then update density (Q)

      do k = 0, Nhat
         do l = 0,Mhat

            S = 0d0

            do j = 0, Mhat
               do i = 0, Nhat

                  do a = 1, 2

                     S = S - UG3(it,a,i,j)*CG5(k,l,i,j,a)/AG4(k,l,k,l)*dt

                  enddo
               enddo
            enddo
            UG3(0,0,k,l) = S + Qn(k,l) !***********************
        
           

        enddo
      enddo
      endif
          
   end

    !$$ Update the local velocity using the global
    subroutine &
      updateU(N,it,alphamax,betamax,UG3,U,imax,F)
      implicit none
      integer          :: N,alphamax,betamax,it,imax
      integer          :: a,i,j,k,l,al,be,Nhat,Mhat,alpha,F   
      double precision :: &
                          U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                          UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1))

    !! update the local velocity from the global


       Nhat = N*(alphamax+1)
       Mhat = N*(betamax+1)
      alpha = 0

     
      
      do a = F, 3             !update from F=0 if out of iteration else from F=1
         do l = 0, Mhat       !(i.e. include density in update if out of iteration)
            do k = 0, Nhat
               
               al = (k-1)/N   !unfortunately this does not work for N=1
               be = (l-1)/N
               i = mod(k+al,N+1)
               j = mod(l+be,N+1)

               U(it,al,be,a,i,j) = UG3(it,a,k,l)
            
            enddo
         enddo
      enddo

     
     do a = F,3  !update from F=0 if out of iteration else from F=1
        
       do al = 1, alphamax
            do be = 0, betamax
               do j = 0 , N 
                  U(it,al,be,a,0,j) = U(it,al-1,be,a,N,j) !resets x-direction duplicates
               enddo                                      !on element boundaries
            enddo
             
         enddo

         do be = 1, betamax
            do al = 0, alphamax
               do i = 0, N
                  U(it,al,be,a,i,0) = U(it,al,be-1,a,i,N)!resets y-direction duplicates
               enddo                                     !on element boundaries
            enddo
               
         enddo
      enddo 

  endsubroutine      

          
!$$ This routine calcs div u, grad u and the solvent part of the stress (TG4)
subroutine calcVG2X(it,imax,N,alphamax,betamax,mu_s,nu,CG5,UG3,AG4,TG4,&
                    gradu,graduN,divu)


 
      implicit none
      integer             :: i,j,k,l,a,b,c,alphamax,betamax,N,imax,it
        double precision  :: S1,S2,nu,&
                             TG4(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                             gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                             graduN(0:N*(alphamax+1),0:N*(betamax+1),2,2),&     
                             divu(0:N*(alphamax+1),0:N*(betamax+1)),h,&
                             CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&
                             AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
                             UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),& 
                             Usmth(0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                             mu_s(0:N*(alphamax+1),0:N*(betamax+1))                                                     !,pressure
  
     
!************************************
!*  calculate kinematic stress TG4  *
!************************************    
      
 
      graduN = gradu

      gradu = 0d0
       divu = 0d0

      do k = 0, N*(alphamax+1)
      do l = 0, N*(betamax+1)
         
         do a = 1, 2
         do b = 1, 2

            S1 = 0d0
            S2 = 0d0

            do i = 0, N*(alphamax+1)
            do j = 0, N*(betamax+1)

               S1 = S1 + UG3(it,a,i,j)*CG5(k,l,i,j,b)

            do c = 1, 2

               S2 = S2 + UG3(it,c,i,j)*CG5(k,l,i,j,c)
            
            enddo
            enddo
            enddo

            gradu(k,l,a,b) = S1/AG4(k,l,k,l) !Gradient of velocity
             
                 divu(k,l) = S2/AG4(k,l,k,l) !Divergence of velocity
         
            enddo
         enddo
          
          !Solvent part of extra stress
         do a = 1, 2
         do b = 1, 2
      
            TG4(k,l,a,b) = mu_s(k,l)*(gradu(k,l,a,b) + gradu(k,l,b,a)) + nu*divu(k,l)*h(a,b)

         enddo
         enddo
      
 
      enddo
      enddo
      
      endsubroutine

     

!$ Generates vertices (xq, yq) of elements for spectral element mesh
subroutine mesh_generator(L, H, alphamax, betamax, fl, fh, no_outer, switch, xq, yq) 
implicit none



integer                      :: alpha, beta, alphamax, betamax, switch, no_outer,&
                                deltamax, deltamin, gamamax, gamamin, i, no_inner_x, no_inner_y


double precision             :: L, H, Lo, dL, Ho, Ox, Oy, vertices(0:alphamax, 0:betamax, 2, 4), fl, fh,&
                                re_no_inner_x, re_no_inner_y, re_alphamax, re_betamax,&
                                xq(0:alphamax, 0:betamax, 5), yq(0:alphamax, 0:betamax, 5)


!$$ This routine assumes an increasing level of refinemnt within a number of "boxes"
!$$ around the bubble near the wall

!$$ Innermost box is a regular rectangular array of equally-sized elements.
!$$ Outer boxes contain quadrilateral elements which fan outwards from
!$$ innermost box to next outer box.

 re_alphamax = real(alphamax) !float point versions of the number of elements 
  re_betamax = real(betamax)

write(*,*) 'al,be=',re_alphamax, re_betamax


if(switch.eq.1)then
      fl = (re_alphamax + 1d0)/(re_alphamax - 1d0)  ! Ratio of lengths of inner and outer boxes fl > 1
      fh = (re_betamax + 1d0)/(re_betamax)          ! Ratio of heights of inner and outer boxes

 write(*,*) 'fl,fh=',fl, fh

endif
     

      Lo = L/fl !Length of inner box
      Ho = H/fh !Height of inner box

      dL = (L - Lo)/2d0 !Length of gap either side of inner or outer box 
      
      Ox = 0d0  !Coordinates of lower left hand corner of outer box
      Oy = 0d0  

  gamamin = 0          !Label of element of outer box
  gamamax = alphamax

 deltamin = 0          !Label of element of outer box
 deltamax = betamax

 !Build mesh from bottom left corner of outermost box
 !(the domain boundary)
 do i = 1, no_outer

 !Calulates coordinates of vertices of elements in outer mesh (fanning elements).
 call outer_mesh(L, Lo, H, Ho, Ox, Oy, gamamin, gamamax, deltamin, deltamax, alphamax, betamax, vertices)


 !Reset inner box as next outer box before returning to
 !outer_mesh subroutine again.
 !For all cases considered in Lind and Phillips (2011), the 
 !number of outer meshes is one. 

  L = Lo
  H = Ho

 gamamin = gamamin + 1
 gamamax = gamamax - 1
deltamax = deltamax - 1

  no_inner_x = gamamax - gamamin + 1
  no_inner_y = deltamax - deltamin + 1
  
 re_no_inner_x = real(no_inner_x)
 re_no_inner_y = real(no_inner_y)

if(switch.eq.1)then
  fl =  re_no_inner_x/(re_no_inner_x - 2d0)
  fh =  re_no_inner_y/(re_no_inner_y - 1d0)
endif

 Lo = L/fl
 Ho = H/fh

 Ox = Ox + dL 
 dL = (L - Lo)/2d0

 enddo

 !Calulates coordinates of vertices of elements in innermost mesh
 !(regular rectangular arrangement)
 call inner_mesh(L, H, Ox, Oy, gamamin, gamamax, deltamin, deltamax, alphamax, betamax, vertices)

    

      do alpha = 0, alphamax
        do beta = 0, betamax
          do i = 1, 4
           
            xq(alpha,beta,i) = vertices(alpha,beta,1,i)
            yq(alpha,beta,i) = vertices(alpha,beta,2,i)
           
          enddo
          
            xq(alpha,beta,5) = vertices(alpha,beta,1,1)
            yq(alpha,beta,5) = vertices(alpha,beta,2,1)
          
          enddo
       enddo

      
endsubroutine


! Calculates vertices of elements in inner refined mesh - a regular rectangular arrangement
subroutine inner_mesh(L, H, Ox, Oy, gamamin, gamamax, deltamin, deltamax, alphamax, betamax, vertices)
implicit none

  integer                      :: alpha, beta,  deltamax, deltamin, gamamax, gamamin,&
                                  i ,j, no_inner_x, no_inner_y, betamax, alphamax


  double precision             :: L, H, Ox, Oy, vertices(0:alphamax, 0:betamax, 2, 4),&
                                  ii, jj, re_no_inner_x, re_no_inner_y

       no_inner_x = gamamax - gamamin + 1
       no_inner_y = deltamax - deltamin + 1

    re_no_inner_x = real(no_inner_x)
    re_no_inner_y = real(no_inner_y)  

  i = 0
 do alpha = gamamin, gamamax
   j = 0
   do beta = deltamin, deltamax
    
    ii = real(i)
    jj = real(j)
    
     

     vertices(alpha, beta, 1, 1) = Ox + ii*L/(re_no_inner_x)
     vertices(alpha, beta, 2, 1) = Oy + jj*H/(re_no_inner_y)

     vertices(alpha, beta, 1, 2) = Ox + (ii + 1d0)*L/(re_no_inner_x)
     vertices(alpha, beta, 2, 2) = Oy + jj*H/(re_no_inner_y)

     vertices(alpha, beta, 1, 3) = Ox + (ii + 1d0)*L/(re_no_inner_x)
     vertices(alpha, beta, 2, 3) = Oy + (jj + 1d0)*H/(re_no_inner_y)

     vertices(alpha, beta, 1, 4) = Ox + ii*L/(re_no_inner_x)
     vertices(alpha, beta, 2, 4) = Oy + (jj + 1d0)*H/(re_no_inner_y)

     j = j + 1

   enddo
   i = i + 1
   enddo

return
end subroutine   

! Calculates vertices of elements in outer mesh.
! Elements fan outwards (increasing in size) from a specified inner box
subroutine outer_mesh(L, Lo, H, Ho, Ox, Oy, gamamin, gamamax, deltamin, deltamax, alphamax, betamax, vertices)

implicit none

integer                      :: alpha, beta, deltamax, deltamin, gamamax, gamamin,&
                                alphamax, betamax, i, j

double precision             :: L, Lo, dL, H, Ho, Ox, Oy, vertices(0:alphamax, 0:betamax, 2, 4),&
                                ii, jj


dL = (L - Lo)/2d0

j = 0

do beta = deltamin, deltamax

      jj = real(j)

  vertices(gamamin, beta, 1, 1) = Ox 
  vertices(gamamin, beta, 2, 1) = Oy + jj*H/(deltamax + 1d0)
  
  vertices(gamamin, beta, 1, 2) = Ox + dL
  vertices(gamamin, beta, 2, 2) = Oy + jj*Ho/deltamax
  
  vertices(gamamin, beta, 1, 3) = Ox + dL
  vertices(gamamin, beta, 2, 3) = Oy + (jj + 1d0)*Ho/deltamax
  
  vertices(gamamin, beta, 1, 4) = Ox
  vertices(gamamin, beta, 2, 4) = Oy + (jj + 1d0)*H/(deltamax + 1d0)

 if(beta.eq.deltamax)then

   vertices(gamamin, beta, 1, 3) = Ox + L/(gamamax - gamamin + 1d0)
   vertices(gamamin, beta, 2, 3) = Oy + H
   
 endif

  vertices(gamamax, beta, 1, 1) = Ox + dL + Lo
  vertices(gamamax, beta, 2, 1) = Oy + jj*Ho/(deltamax)
  
  vertices(gamamax, beta, 1, 2) = Ox + L
  vertices(gamamax, beta, 2, 2) = Oy + jj*H/(deltamax + 1d0)
  
  vertices(gamamax, beta, 1, 3) = Ox + L
  vertices(gamamax, beta, 2, 3) = Oy + (jj + 1d0)*H/(deltamax + 1d0)
  
  vertices(gamamax, beta, 1, 4) = Ox + dL + Lo
  vertices(gamamax, beta, 2, 4) = Oy + (jj + 1d0)*Ho/(deltamax)

 if(beta.eq.deltamax)then

   vertices(gamamax, beta, 1, 4) = Ox + (gamamax - gamamin)*L/(gamamax - gamamin + 1d0)
   vertices(gamamax, beta, 2, 4) = Oy + H
    
 endif

   j = j + 1

 enddo

   i = 0

 do alpha = gamamin, gamamax

   ii = real(i)
   
  vertices(alpha, deltamax, 1, 1) = Ox + dL + (ii - 1d0)*Lo/(gamamax - gamamin - 1d0)
  vertices(alpha, deltamax, 2, 1) = Oy + Ho
  
  vertices(alpha, deltamax, 1, 2) = Ox + dL + ii*Lo/(gamamax - gamamin - 1d0)
  vertices(alpha, deltamax, 2, 2) = Oy + Ho
  
  vertices(alpha, deltamax, 1, 3) = Ox + (ii + 1d0)*L/(gamamax - gamamin + 1)
  vertices(alpha, deltamax, 2, 3) = Oy + H
  
  vertices(alpha, deltamax, 1, 4) = Ox + ii*L/(gamamax - gamamin + 1)
  vertices(alpha, deltamax, 2, 4) = Oy + H

  if(alpha.eq.gamamin)then

  vertices(alpha, deltamax, 1, 1) = Ox 
  vertices(alpha, deltamax, 2, 1) = Oy + deltamax*H/(deltamax + 1d0)

  elseif(alpha.eq.gamamax)then

  vertices(alpha, deltamax, 1, 2) = Ox + L
  vertices(alpha, deltamax, 2, 2) = Oy + deltamax*H/(deltamax + 1d0)

  endif

  i = i + 1
  
  enddo  

 return
 end subroutine


!$$ Write commonly used tensors (fixed for a particular mesh) to a file so they only
!$$ have to be read in and not recalculated when performing a new run on same mesh.
!$$ This is slightly faster than recalculating them from scratch for every run.
 subroutine write_mesh_details(N, alphamax, betamax, AG4, CG5, NiCG5, MiCG5, nnmax, mmax, iCG5)
  implicit none

    integer       :: nnmax, MiCG5, NiCG5, i, j, k, l, N, alphamax, betamax, b 
 
    integer       :: iCG5(NiCG5,MiCG5,2), mmax(NiCG5)

 double precision :: AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
                     CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2)
  



     
 do i = 0, N*(alphamax+1)
  do j = 0, N*(betamax+1)
   do k = 0, N*(alphamax+1)
    do l = 0, N*(betamax+1)
 
       write(11,*) AG4(i,j,k,l)

      do b = 1,2

       write(12,*) CG5(i,j,k,l,b)

       enddo
      enddo
    enddo
   enddo
  enddo

   close(11)
   close(12)

  write(13,*) NiCG5
  write(14,*) MiCG5
  write(15,*) nnmax

  close(13)
  close(14)
  close(15)

  do i = 1, NiCG5
   
    write(16,*) mmax(i)

    do j = 1, MiCG5
      do b = 1,2
 
      write(17,*) iCG5(i,j,b)

       enddo
    enddo
   enddo

   close(16)
   close(17)
    
  end subroutine
    

  !$$ Read in mesh details/tensors from previously written mesh files
  subroutine read_mesh_details_a(N, alphamax, betamax, AG4, CG5, NiCG5, MiCG5, nnmax)
  implicit none

          integer       :: nnmax, MiCG5, NiCG5, i, j, k, l, N, alphamax, betamax, b 


       double precision :: AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&
                           CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2)
  

     


do i = 0, N*(alphamax+1)
  do j = 0, N*(betamax+1)
   do k = 0, N*(alphamax+1)
    do l = 0, N*(betamax+1)
 
       read(11,*) AG4(i,j,k,l)

      do b = 1,2

       read(12,*) CG5(i,j,k,l,b)

       enddo
      enddo
    enddo
   enddo
  enddo

  close(11)
  close(12)

  read(13,*) NiCG5
  read(14,*) MiCG5
  read(15,*) nnmax

 close(13)
 close(14)
 close(15) 

  return
  end subroutine
         
  subroutine read_mesh_details_b(NiCG5, MiCG5, mmax, iCG5)
  implicit none

  integer      :: NiCG5, MiCG5, i, j, b
  integer      :: mmax(NiCG5), iCG5(NiCG5, MiCG5, 2)

   do i = 1, NiCG5
   
    read(16,*) mmax(i)

    do j = 1, MiCG5
      do b = 1,2
 
      read(17,*) iCG5(i,j,b)

       enddo
    enddo

   enddo

  close(16)
  close(17)

  return
  end subroutine


!$$ Given the number of particle cells in x and y direction (Nx, Ny resp), 
!$$ and the number of particles per cell (Np^2), Calculate total number
!$$ of particles

subroutine calc_no_particles(TNP, Npx, Npy, Nx, Ny, Np, dx_p, dy_p, dx_c, dy_c, LL, HH)
implicit none

!$ calculating numbers associted with marker particles

integer           :: TNP, Npx, Npy, Nx, Ny, Np

double precision  :: dx_p, dy_p, dx_c, dy_c, LL, HH, re_Nx, re_Ny

      re_Nx = real(Nx) !Number of cells in x
      re_Ny = real(Ny) !Number of cells in y
     
  
      Npx = (Np - 1)*Nx + 1  !! No of particles in x direction
      Npy = (Np - 1)*Ny + 1  !! No of particles in y

      TNP = Npx*Npy !! Total number of particles

     dx_p = LL/(Npx - 1d0) !! Spacing between particles on reg grid in x
     dy_p = HH/(Npy - 1d0) !! Spacing between particles on reg grid in y

     dx_c = LL/(re_Nx) !! Spacing btw cell edges (width of cells)
     dy_c = HH/(re_Ny) !! Spacing btw cell edges (height of cells)


return
endsubroutine


!$$ Counter over all particles and initialise their colour and their position
subroutine calc_ptcl_pstn(TNP, Np1, Npx, Npy, xp, yp, Cp, xc, yc, rad, dx_p, dy_p)
implicit none

   integer       :: Npx, Npy, i, j, l, TNP, Np1

double precision :: xp(TNP), yp(TNP), Cp(TNP,2), xc, yc, rad, rr, dx_p, dy_p

Np1 = 0
  l = 1

!$ Initial particle position (xp,yp)

 do i = 1, Npx
  do j = 1, Npx

    xp(l) = (i - 1d0)*dx_p
    yp(l) = (j - 1d0)*dy_p

 rr = (xp(l) - xc)**2 + (yp(l) - yc)**2

 if(rr.le.rad**2)then !$$ If particle within bubble

 !$ Set colour of each particle
 !$ Cp(1)=1 if bubble 
 !$ Cp(1)=0 otherwise

   Cp(l,1) = 1d0 !$ Discrete colour of particle l. 
   Cp(l,2) = 0d0

  Np1 = Np1 + 1

else

   Cp(l,1) = 0d0
   Cp(l,2) = 1d0

endif

  l = l + 1

 enddo
enddo

return
endsubroutine


  !$$ Output to file the particle data (testing purposes)
   subroutine writepoints_ptcl(TNP, Np1, Npx, Npy, xp, yp, Cp, time2)
      implicit none
      integer           :: TNP, i, Npx, Npy, Np1
      double precision  :: xp(TNP), yp(TNP), Cp(TNP,2), time2
                            
                           
      write(999,*) 'VARIABLES = "x", "y"'
      write(999,*) 'ZONE T="time:',time2,'" I =',Np1,', J =',1,',F=POINT' 

      do i = 1, TNP  
                 
      if(Cp(i,1).eq.1d0)then
        write(999,*) xp(i),yp(i)
      endif

      enddo
      
       endsubroutine

 !$$ Given a velocity, update marker particle positions in time using simple Euler stepping
 subroutine update_ptcl(TNP, alphamax, betamax, N, imax, P, xq, yq, Len, Height, xp, yp, U, U_jet, dt, eps, time, jet_flag)
 implicit none

        integer  :: TNP, alphamax, betamax, k, N, alpha, beta, a, imax,time, jet_flag

double precision :: xx, yy, xx_new, yy_new, xp(TNP), yp(TNP), xi, zeta, eps,&
                    S(0:3),P(0:N), xq(0:alphamax, 0:betamax,5), yq(0:alphamax, 0:betamax,5),&
                    U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N), dt, Up(TNP,2), U_jet(2), Len, Height

 do k = 1, TNP !$ Counter over all particles

   xx = xp(k)
   yy = yp(k)

   !calculate local coordinates (xi,zeta) on element
     call x2xi(alphamax,betamax,alpha,beta, xx, yy, xq, yq, xi, zeta, eps) 
   !calculate velocity using Eqn 32     
     call calcUxy(U,alpha,beta,xi,zeta,0,N,P,alphamax,betamax,S,imax, eps)
           
       do a = 1,2
          Up(k,a)=S(a)

       if(jet_flag.eq.1)then
          U_jet(a) = S(a)
       endif

       enddo

    !$$ Update particle positions 
!!  if(time.eq.1)then
        xx_new = xx + dt*Up(k,1)
        yy_new = yy + dt*Up(k,2)
!!  else
!!        xx_new = xx + dt*(1.5d0*Up(k,1) - 0.5d0*Upn(k,1))  !!2nd order adams-basforth
!!        yy_new = yy + dt*(1.5d0*Up(k,2) - 0.5d0*Upn(k,2))  
!!  endif

!$$ If particles leave domain, reflect inwards
!$$ by amount they leave.

if((xx_new.lt.eps).or.(xx_new.gt.Len).or.&
(yy_new.lt.eps).or.(yy_new.gt.Height))then

            if(xx_new.le.eps)then
              
               xx_new = eps + (eps - xx_new)              
           
              endif

              if(xx_new.ge.Len)then
            
               xx_new = Len - (xx_new - Len)
       
              endif

               if(yy_new.le.eps)then
              
               yy_new = eps + (eps - yy_new)              
                
               endif

                if(yy_new.ge.Height)then
               
                yy_new = Height - (yy_new - Height)             

                endif

           endif

      xp(k) = xx_new
      yp(k) = yy_new

enddo

! do k = 1, TNP
!  do a = 1,2
!   Upn(k,a) = Up(k,a)
!  enddo
! enddo



return
endsubroutine

!$$ Calculate the interpolated colour function Cij from discrete colour function Cp
!$$ (See Eqn 51)

subroutine calc_colour(N, TNP, alphamax, betamax, Cp, Cij, dx_c, dy_c, xp, yp, xg, yg)
implicit none

         integer :: i, j, p, TNP, alphamax, betamax, N
double precision :: Sum1, Sum2, S, dx_c, dy_c, Cp(TNP,2),& 
                    Cij(0:N*(alphamax + 1),0:N*(betamax + 1),2),&
                    xp(TNP), yp(TNP),&
                    xg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    yg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    xx, yy, xxg, yyg
                   
!!!write(*,*) 'x,dy',dx_c,dy_c                   
 
!!$ calculate colour at each GLL point through averaging of all coloured particles
!!$ within square of area 4*dx_c*dy_c around node

do i = 0, N*(alphamax + 1)
do j = 0, N*(betamax + 1)

 Sum1 = 0d0
 Sum2 = 0d0

! If on boundary, then colour is fixed (reasonable as no slip bc)
if(i.eq.0.or.i.eq.(N*(alphamax + 1)).or.j.eq.0.or.j.eq.(N*(betamax + 1)))then

Cij(i,j,1) = 0d0
Cij(i,j,2) = 1d0

else

do p = 1, TNP

  xx = xp(p)
  yy = yp(p)
 xxg = xg(i,j)
 yyg = yg(i,j) 

!$$ Sum bilinear interpolating function, S.
  Sum1 = Sum1 + S(xx,xxg,yy,yyg,dx_c,dy_c)*Cp(p,1)
 
  Sum2 = Sum2 + S(xx,xxg,yy,yyg,dx_c,dy_c)

enddo

!!!write(*,*)'sum',i,j,Sum1,Sum2

if(Sum2.eq.0d0)then
write(*,*) 'no marker particles around node',i,j
stop
endif

Cij(i,j,1) = Sum1/Sum2

Cij(i,j,2) = 1d0 - Cij(i,j,1)

endif

enddo
enddo



return
endsubroutine

!$$ Output interpolated colour data to file
subroutine writepoints_colour(N, alphamax, betamax, Cij, xg, yg, time2)
implicit none

         integer :: N, alphamax, betamax, i, j

double precision :: Cij(0:N*(alphamax + 1),0:N*(betamax + 1),2),&
                    xg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    yg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    time2
            
      write(998,*) 'VARIABLES = "x", "y","C"'
      write(998,*) 'ZONE T="',time2,'" I =', ((betamax + 1)*N + 1),', J =', ((alphamax + 1)*N + 1),',F=POINT' 
        
 do i = 0, N*(alphamax + 1)
  do j = 0, N*(betamax + 1)
    write(998,*) xg(i,j),yg(i,j),Cij(i,j,1)
  enddo
 enddo
    
endsubroutine 

!$$ Output viscosity to file (testing purposes)
subroutine writepoints_mu(N, alphamax, betamax, mu, xg, yg, time2)
implicit none

         integer :: N, alphamax, betamax, i, j

double precision :: mu(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    xg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    yg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    time2
            
      write(991,*) 'VARIABLES = "x", "y","mu"'
      write(991,*) 'ZONE T="',time2,'" I =', ((betamax + 1)*N + 1),', J =', ((alphamax + 1)*N + 1),',F=POINT' 
        
 do i = 0, N*(alphamax + 1)
  do j = 0, N*(betamax + 1)
    write(991,*) xg(i,j),yg(i,j),mu(i,j)
  enddo
 enddo
    
endsubroutine 

!$$ General output routine for writing tensors of same (common) indexing structure.
subroutine writepoints_tensors(N, alphamax, betamax, Tensor, xg, yg, time2, a)
implicit none

         integer :: N, alphamax, betamax, i, j, a

double precision :: Tensor(0:N*(alphamax + 1),0:N*(betamax + 1),2,2),&
                    xg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    yg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    time2

      write(a,*) 'VARIABLES = "x", "y","T11","T12","T21","T22"'
      write(a,*) 'ZONE T="',time2,'" I =', ((betamax + 1)*N + 1),', J =', ((alphamax + 1)*N + 1),',F=POINT' 
        
 do i = 0, N*(alphamax + 1)
  do j = 0, N*(betamax + 1)
    write(a,*) xg(i,j),yg(i,j),Tensor(i,j,1,1),Tensor(i,j,1,2),Tensor(i,j,2,1),Tensor(i,j,2,2)
  enddo
 enddo


    
endsubroutine 

!$$ Output gloabl velocity data to file.
subroutine writepoints_UG(N, alphamax, betamax, UG, xg, yg, time2, a)
implicit none

         integer :: N, alphamax, betamax, i, j, a

double precision :: UG(0:3,0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    xg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    yg(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    time2

      write(a,*) 'VARIABLES = "x", "y","U1","U2"'
      write(a,*) 'ZONE T="',time2,'" I =', ((betamax + 1)*N + 1),', J =', ((alphamax + 1)*N + 1),',F=POINT' 
        
 do i = 0, N*(alphamax + 1)
  do j = 0, N*(betamax + 1)
    write(a,*) xg(i,j),yg(i,j),UG(1,i,j),UG(2,i,j)
  enddo
 enddo


    
endsubroutine 


!$$ Calculate interpolated  solvent viscosity using colour function (Eqn 50)
subroutine calc_mu_s(N, TNP, alphamax, betamax, Cij, mu_s, mu_s0)
implicit none


         integer :: TNP, alphamax, betamax,N, i, j

double precision :: mu_s(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    Cij(0:N*(alphamax + 1),0:N*(betamax + 1),2),&
                    mu_s0(2)

  do i = 0, N*(alphamax + 1)
   do j = 0, N*(betamax + 1)

     mu_s(i,j) = mu_s0(1)*Cij(i,j,1) + mu_s0(2)*Cij(i,j,2)

   enddo
  enddo



return
endsubroutine

!$$ Calculate interpolated polymeric viscosity using colour function (Eqn 50)
subroutine calc_mu_p(N, TNP, alphamax, betamax, Cij, mu_p, mu_p0)
implicit none


         integer :: TNP, alphamax, betamax,N, i, j

double precision :: mu_p(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    Cij(0:N*(alphamax + 1),0:N*(betamax + 1),2),&
                    mu_p0(2)

  do i = 0, N*(alphamax + 1)
   do j = 0, N*(betamax + 1)
  
     mu_p(i,j) = mu_p0(1)*Cij(i,j,1) + mu_p0(2)*Cij(i,j,2)

   enddo
  enddo



return
endsubroutine

!$$ Calculate interpolated relaxation time using colour function (Eqn 50)
subroutine calc_lambda(N, TNP, alphamax, betamax, Cij, lambda1, lambda10)
implicit none


         integer :: TNP, alphamax, betamax,N, i, j

double precision :: lambda1(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    Cij(0:N*(alphamax + 1),0:N*(betamax + 1),2),&
                    lambda10(2)

  do i = 0, N*(alphamax + 1)
   do j = 0, N*(betamax + 1)
  
      lambda1(i,j) = lambda10(1)*Cij(i,j,1) + lambda10(2)*Cij(i,j,2)
   
    enddo
  enddo




return
endsubroutine

!$$ Calculate interpolated "total viscosity" (not a viscosity but a material parameter made up of a number 
!$$ of others).
subroutine calc_mu_t(N, alphamax, betamax, mu_s, mu_p, lambda1, dt, mu_t)
implicit none

         integer :: N, alphamax, betamax, i, j
double precision :: mu_t(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    mu_s(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    mu_p(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    lambda1(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    C1,dt

 do i = 0, N*(alphamax + 1)
   do j = 0, N*(betamax + 1)

            C1 = 1d0 + lambda1(i,j)/dt

     mu_t(i,j) = mu_s(i,j) + mu_p(i,j)/C1

   enddo
 enddo


return
endsubroutine

!$$ Calc Gamma tensors. Tensors which contain material and polymeric stress data which 
!$$ are used in the stiffness matrix Eqn 44.
subroutine calc_Gamma(N, alphamax, betamax, GammaN, GammaP, mu_t, lambda1, Pstress,dt)
implicit none

         integer :: N, alphamax, betamax, i, j, a, b,c,d
double precision :: mu_t(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    lambda1(0:N*(alphamax + 1),0:N*(betamax + 1)),&
                    GammaN(0:N*(alphamax + 1),0:N*(betamax + 1),2,2),&
                    GammaP(0:N*(alphamax + 1),0:N*(betamax + 1),2,2,2,2),&
                    Pstress(0:N*(alphamax + 1),0:N*(betamax + 1),2,2),&
                    C1,dt,h

 
   do i = 0, N*(alphamax + 1)
    do j = 0, N*(betamax + 1)

                 C1 = 1d0 + lambda1(i,j)/dt

       do a = 1, 2
        do b = 1, 2

          GammaN(i,j,a,b) = mu_t(i,j)*h(a,b) !!+ Pstress(i,j,a,b)*(lambda1(i,j)/C1)

         do c = 1, 2
          do d = 1, 2

           GammaP(i,j,a,b,c,d) = (h(a,b)*Pstress(i,j,c,d) + h(a,c)*Pstress(i,j,b,d))*(lambda1(i,j)/C1)
        
         enddo
        enddo

       enddo
      enddo

  enddo 
 enddo

return
endsubroutine


 !$$ Bilinear interpolation function S (see Eqn 51)
  double precision function & 
   S(xx, xxi, yy, yyi, dx_c, dy_c)
   implicit none
  
!!$ Interpolation function

   double precision  :: xx, xxi, yy, yyi, dx_c, dy_c, a, b    

    a = abs((xx - xxi)/dx_c)  
    b = abs((yy - yyi)/dy_c)

    if(a.gt.0d0.and.a.lt.1d0.and.b.gt.0d0.and.b.lt.1d0)then

    S = (1d0 - a)*(1d0 - b)

    else

    S = 0d0        

    endif

    return
    end


   
                        
  !$$ Now redunant rountine for calculating polymeric stress
  !$$ (Oldroyd B /UCM). Not used currently
     subroutine calcOldBstress(N,alphamax,betamax,imax,dt,gradu,&
                               mu_p,lambda1,Pstress,Pstress2,CG5,UG3,AG4)

     
      implicit none
      
      integer ::  N,alphamax,betamax,imax,i,j,k,l,a,b,c,d,refine_dt,&
                  IFAIL,IPIV(3)

     double precision :: DG3(4,0:N*(alphamax+1),0:N*(betamax+1)),& !gammadot
                         UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),& 
                         AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&          
                         gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),h,dt,&
                         M(4,4),V(3),MM(3,3),S(2,2),&
                         Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         Pstress2(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&     
                         lambda1,mu_p,ddt



     ddt = dt

      IFAIL = 0      

          S = 0d0

     do l = 0, N*(betamax + 1) 
     do k = 0, N*(alphamax + 1)      

               do d=1,2
               do c=1,2
               do b=1,2
               do a=1,2

             M(a + 2*(b - 1), c + 2*(d - 1))=h(a,c)*h(b,d)&
             *(1d0 + ddt/lambda1) - 0.5d0*ddt*&
              (&
                gradu(k,l,a,c)*h(b,d) + gradu(k,l,b,d)*h(a,c) +&
                gradu(k,l,a,d)*h(c,b) + gradu(k,l,b,c)*h(d,a) &
              )

               enddo
               enddo
                  DG3(c + 2*(d-1),k,l) = h(c,d)
               enddo
               enddo               

                  do a = 1, 2 !convection term u.grad(T)
                   do b = 1, 2
 
                     S(a,b) = 0d0
 
                     do j = 0, N*(betamax+1)
                     do i = 0, N*(alphamax+1)
                     do c = 1, 2
                        S(a,b) = S(a,b) - UG3(0,c,k,l)*Pstress(i,j,a,b)*CG5(k,l,i,j,c)
                     enddo
                     enddo
                     enddo
 
                  enddo
                  enddo

!!$               write(*,*) 'before set v'
!!$
!!$               write(*,*) 'P',Pstress(k,l,1,2),Pstress(k,l,1,2),Pstress(k,l,1,4)
!!$               write(*,*) 'mup',mu_p(k,l)
!!$               write(*,*) 'Dg3',DG3(1,k,l),DG3(2,k,l),DG3(4,k,l)
!!$               write(*,*) 's',S(1,1),S(1,2),S(2,2)

                  V(1) = Pstress(k,l,1,1) + ddt*((mu_p/lambda1)*DG3(1,k,l) + S(1,1)/AG4(k,l,k,l))
                  V(2) = Pstress(k,l,1,2) + ddt*((mu_p/lambda1)*DG3(2,k,l) + S(1,2)/AG4(k,l,k,l))
                  V(3) = Pstress(k,l,2,2) + ddt*((mu_p/lambda1)*DG3(4,k,l) + S(2,2)/AG4(k,l,k,l))! 
              
             !!  write(*,*) 'after set V'
                                      !
               MM(1,1) = M(1,1)
               MM(2,1) = M(2,1)
               MM(3,1) = M(4,1)
               MM(1,2) = M(1,2) + M(1,3)
               MM(2,2) = M(2,2) + M(2,3)   !contracts M(4 x 4) to MM(3 x 3)
               MM(3,2) = M(4,2) + M(4,3)
               MM(1,3) = M(1,4)
               MM(2,3) = M(2,4)
               MM(3,3) = M(4,4)   !DGESV alters MM  

         !!!     write(*,*) 'before DGESV'                                     
                                                                              !solve for stress
               call DGESV(3,1,MM,3,IPIV,V,3,IFAIL)                            
                                       
         !!!     write(*,*) 'after DGESV'
                                      
                  Pstress(k,l,1,1) = V(1) 
                  Pstress(k,l,1,2) = V(2) 
                  Pstress(k,l,2,1) = V(2) 
                  Pstress(k,l,2,2) = V(3)     

                  Pstress2(k,l,1,1) = Pstress(k,l,1,1) - (mu_p/lambda1)
                  Pstress2(k,l,1,2) = Pstress(k,l,1,2)
                  Pstress2(k,l,2,1) = Pstress(k,l,2,1)
                  Pstress2(k,l,2,2) = Pstress(k,l,2,2) - (mu_p/lambda1)               

      enddo !k
      enddo !l

     end

  !$$ Another redundant routine for calculating polymeric stress
   subroutine calcOldBstress2(N,alphamax,betamax,imax,dt,gradu,&
                              mu_p,lambda1,Pstress,Bn,DG3)

     
      implicit none
      
      integer ::  N,alphamax,betamax,imax,i,j,k,l,a,b,c,d
                  

     double precision :: gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),h,dt,&
                         Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         lambda1(0:N*(alphamax+1),0:N*(betamax+1)),&
                         mu_p(0:N*(alphamax+1),0:N*(betamax+1)),&
                         Bn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         DG3(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         ct(0:N*(alphamax+1),0:N*(betamax+1))

      
    !! ct = 1d0 + lambda1/dt
       
      

     do l = 0, N*(betamax + 1) 
     do k = 0, N*(alphamax + 1)   
   
       ct(k,l) = 1d0 + 3d0*lambda1(k,l)/(2d0*dt)
                       
            do b = 1, 2   !!rate of Deformation
             do a = 1, 2
                 
               DG3(k,l,a,b) = mu_p(k,l)*(gradu(k,l,a,b) + gradu(k,l,b,a))

              enddo
             enddo            

               
               Pstress(k,l,1,1) = DG3(k,l,1,1)/ct(k,l) + Bn(k,l,1,1)
               Pstress(k,l,1,2) = DG3(k,l,1,2)/ct(k,l) + Bn(k,l,1,2)
               Pstress(k,l,2,1) = DG3(k,l,2,1)/ct(k,l) + Bn(k,l,2,1)
               Pstress(k,l,2,2) = DG3(k,l,2,2)/ct(k,l) + Bn(k,l,2,2)


         enddo !k
         enddo !l

   endsubroutine

 !$$ The currently used subroutine for calculating the polymeric stress, 
 !$$ using Oldroyd B/ UCM. 

 subroutine calcOldBstress3(N,alphamax,betamax,dt,gradu,mu_p,lambda1,Pstress,Pstressn)

     
      implicit none
      
      integer ::  N,alphamax,betamax,imax,i,j,k,l,a,b,c,d,refine_dt,&
                  IFAIL,IPIV(3)

     double precision :: DG3(2,2),& !gammadot
                         gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),h,dt,&
                         M(4,4),V(3),MM(3,3),&
                         Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         Pstressn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                         lambda1(0:N*(alphamax+1),0:N*(betamax+1)),&
                         mu_p(0:N*(alphamax+1),0:N*(betamax+1)),ddt
                                          



     ddt = dt

      IFAIL = 0      

    do l = 0, N*(betamax + 1) 
     do k = 0, N*(alphamax + 1)      

      !$$ Construct matrix of deformation terms to solve 
      !$$ for Pstress

       do d=1,2
        do c=1,2
         do b=1,2
          do a=1,2

           M(a + 2*(b - 1), c + 2*(d - 1))=h(a,c)*h(b,d)&
           *(1d0 + lambda1(k,l)/dt) - 0.5d0*lambda1(k,l)*&
              (&
                gradu(k,l,a,c)*h(b,d) + gradu(k,l,b,d)*h(a,c) +&
                gradu(k,l,a,d)*h(c,b) + gradu(k,l,b,c)*h(d,a) &
              )

               enddo
               enddo
               
               enddo
               enddo                   


     do a = 1,2
       do b = 1,2 !$$ Rate of deformation + Poly stress at previous particle time/location
            DG3(a,b) = mu_p(k,l)*(gradu(k,l,a,b)+gradu(k,l,b,a)) + (lambda1(k,l)/dt)*Pstressn(k,l,a,b)
       enddo
     enddo               

                  V(1) = DG3(1,1)
                  V(2) = DG3(1,2)
                  V(3) = DG3(2,2)
              

               MM(1,1) = M(1,1)
               MM(2,1) = M(2,1)
               MM(3,1) = M(4,1)
               MM(1,2) = M(1,2) + M(1,3)
               MM(2,2) = M(2,2) + M(2,3)   !contracts M(4 x 4) to MM(3 x 3)
               MM(3,2) = M(4,2) + M(4,3)
               MM(1,3) = M(1,4)
               MM(2,3) = M(2,4)
               MM(3,3) = M(4,4)   !DGESV alters MM  

              IPIV=0

         !$$ Invert matrix of deformation terms
                               
         call DGESV(3,1,MM,3,IPIV,V,3,IFAIL)                            
                                       
        if(Ifail.ne.0)then
         write(*,*) 'ifail in oldbstress',IFAIL
        endif                     
                      
         !$$ And find Pstress (polymeric stress)

                  Pstress(k,l,1,1) = V(1) 
                  Pstress(k,l,1,2) = V(2) 
                  Pstress(k,l,2,1) = V(2) 
                  Pstress(k,l,2,2) = V(3)     

               
      enddo !k
      enddo !l

   endsubroutine


 !$$ Redunant routine from previous version of code (not used currently)
  subroutine calcZG3(N,imax,alphamax,betamax,UG3,CG5,ZG3,mu_p,lambda1)
  implicit none
 
    integer          ::  i,j,a,b,c,p,q,N,alphamax,betamax,imax
      double precision :: &
                          CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),h,&
                          TG4(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                          QG3(0:N*(alphamax+1),0:N*(betamax+1),2),&
                          ZG3(0:N*(alphamax+1),0:N*(betamax+1),2),&
                          UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                          S1,&
                          mu_p(0:N*(alphamax+1),0:N*(betamax+1)), lambda1

   
      do i = 0, N*(alphamax+1)
      do j = 0, N*(betamax+1)
         
        do a = 1, 2

            S1 = 0d0

            do p = 0,N*(alphamax+1)
            do q = 0,N*(betamax+1)

                  S1 = S1 + UG3(0,0,p,q)*CG5(i,j,p,q,a)

            enddo
            enddo

              QG3(i,j,a) = S1

         enddo
      enddo
      enddo

!    YG3 represents the grad q.T term

      do i = 0,N*(alphamax+1)
      do j = 0,N*(betamax+1)
 
         do a = 1, 2

           S1 = QG3(i,j,a)*mu_p(i,j)/lambda1

            if(i.eq.0.or.i.eq.N*(alphamax+1).or.j.eq.0.or.j.eq.N*(betamax + 1))then
              ZG3(i,j,a) = 0d0
            else
              ZG3(i,j,a) = S1
            endif

        enddo
      enddo
      enddo      
      end
  

 !$$ Routine which calculates analytical solution for viscoelastic Poiseuille flow,
 !$$ and compares with numerical solution (not used for bubble problems)  
subroutine cf_analytic(N,alphamax,betamax,imax,U,time,mu,mu_p,lambda1,Height,Len,P,xq,yq)
implicit none

integer          :: i,j, N, alphamax, betamax, imax, n_crit, alpha, beta
double precision :: U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N), yy, xx, crit,& 
                    anal_u, ns, sum, abs_dsum, A, B, nns, pi,&
                    time, mu, Height, dsum, num_u, error,&
                    alpN, betN, S1, S2, GN,&
                    lambda1, mu_p, bet2, aa,&
                    Len, eps,&
                    S(0:3), xi, zeta, P(0:N),&
                    xq(0:alphamax,0:betamax,5),&
                    yq(0:alphamax,0:betamax,5)
                    
                          


         crit = 1d-6
      
           pi = 4d0*datan(1d0)

           S1 = lambda1*(mu + mu_p)
           S2 = lambda1*mu

           yy = 0.5d0

           ns = 1d0
          sum = 0d0
     abs_dsum = 1d0
       n_crit = 0

      do 

       nns = (2d0*ns - 1d0)*pi

      alpN = 1d0 + S2*nns**2 

      bet2 = (1d0 + S2*nns**2)**2 - 4d0*S1*(nns**2)

       if(bet2.lt.0d0)then

      bet2 = -bet2

      betN = dsqrt(bet2)

        GN = dcos(betN*time/(2d0*S1))+((1d0 + (nns**2)*(S2 - 2d0*S1))/betN)*dsin(betN*time/(2d0*S1))       

      else

      betN = dsqrt(bet2)

        aa = betN*time/(2d0*S1)
      
        GN = dcosh(aa) + ((1d0 + (nns**2)*(S2 - 2d0*S1))/betN)*dsinh(aa)

      endif

         A = nns*yy/Height

         B = -alpN*time/(2d0*S1*Height**2)

      dsum = (1d0/nns**3)*sin(A)*exp(B)*GN

       sum = sum  + dsum

  abs_dsum = abs(dsum)

      if(abs_dsum.lt.crit)then
        
        n_crit = n_crit + 1

      else

        n_crit = 0

      endif

         ns = ns + 1d0

       if(abs_dsum.lt.crit.and.n_crit.gt.3)exit

      enddo

       anal_u = -4d0*yy*(yy - 1d0) - 32d0*sum    

           xx = Len/2d0

   call x2xi(alphamax,betamax,alpha,beta, xx, yy, xq, yq, xi, zeta, eps)

   call calcUxy(U,alpha,beta,xi,zeta,0,N,P,alphamax,betamax,S,imax, eps)!also includes calculation for Qn
           
        num_u = S(1) 

    error = abs(anal_u - num_u)

  write(888,*) time, anal_u, num_u, error

return
endsubroutine


!$$ Routine which calculates velocity boundary conditions for Newtonian Poiseuille flow
!$$ Not used for bubble problems
  subroutine calcB3_Newt(N,alphamax, betamax, P, heat0, speed, B3_Left, B3_Right, B3_Top, B3_Bot,  Height, yq, time, mu)
      implicit none
      integer           :: al,be,i,j,l,N, alphamax,betamax, it
      double precision  :: B3_Left(1:N*(betamax + 1) - 1, 3),&
                           B3_Right(1:N*(betamax + 1) - 1, 3),&
                           B3_Top(0:N*(alphamax + 1), 3),&
                           B3_Bot(0:N*(alphamax + 1), 3),&
                           P(0:N), heat0, speed, Height, y,&
                           yq(0:alphamax,0:betamax,5),&
                           mu, yy, nns, sum, dsum, abs_dsum,&
                           time, A, B, pi, crit, ns
                          
 
        !!BOUNDARY CONDITIONS FOR Couette Flow

!! BOUNDARY CONDITIONS ALONG X AXIS

     !!   it = 1d15
      crit = 1d-10
        pi = 4d0*datan(1d0)


        B3_Bot(0, 1) = 0d0            
        B3_Top(0, 1) = 0d0
        
        B3_Bot(0, 2) = 0d0            
        B3_Top(0, 2) = 0d0            
        
        B3_Bot(0, 3) = heat0           
        B3_Top(0, 3) = heat0 

i = 1      
      do al = 0, alphamax
        do  l = 1, N
        
        B3_Bot(i, 1) = 0d0            
        B3_Top(i, 1) = 0d0
        
        B3_Bot(i, 2) = 0d0            
        B3_Top(i, 2) = 0d0            
        
        B3_Bot(i, 3) = heat0           
        B3_Top(i, 3) = heat0 
         
      i = i + 1                        
         enddo
        enddo

   
!!! BOUNDARY CONDITIONS ALONG Y AXIS

 j = 1       
      do be = 0, betamax - 1
        do l = 1, N

           yy = y(0, alphamax, be, betamax, P(0), P(l), yq)

           ns = 1d0
          sum = 0d0
     abs_dsum = 1d0

      do while(abs_dsum.gt.crit)

       nns = (2d0*ns - 1d0)*pi

         A = nns*yy/Height

         B = -(nns**2)*mu*time/(Height**2)

      dsum = (1d0/nns**3)*sin(A)*exp(B)

       sum = sum  + dsum

  abs_dsum = abs(dsum)

         ns = ns + 1d0
       

       enddo

         B3_Left(j, 1) = -4d0*yy*(yy - 1d0) - 32d0*sum      
        B3_Right(j, 1) = -4d0*yy*(yy - 1d0) - 32d0*sum      

    !!     write(*,*) B3_Left(j,1), yy, sum

         B3_Left(j, 2) = 0d0        
        B3_Right(j, 2) = 0d0

         B3_Left(j, 3) = heat0       
        B3_Right(j, 3) = heat0   
     
      j = j + 1   
        enddo
       enddo

    be = betamax

  do l = 1, N-1

            yy = y(0, alphamax, be, betamax, P(0), P(l), yq)

           ns = 1d0
          sum = 0d0
     abs_dsum = 1d0

      do while(abs_dsum.gt.crit)

       nns = (2d0*ns - 1d0)*pi

         A = nns*yy/Height

         B = -(nns**2)*mu*time/(Height**2)

      dsum = (1d0/nns**3)*sin(A)*exp(B)

       sum = sum  + dsum

  abs_dsum = abs(dsum)

         ns = ns + 1d0
     

       enddo

         B3_Left(j, 1) = -4d0*yy*(yy - 1d0) - 32d0*sum      
        B3_Right(j, 1) = -4d0*yy*(yy - 1d0) - 32d0*sum  
   
    !!     write(*,*) B3_Left(j,1), yy, sum

         B3_Left(j, 2) = 0d0        
        B3_Right(j, 2) = 0d0

         B3_Left(j, 3) = heat0       
        B3_Right(j, 3) = heat0   

    j = j + 1

   enddo
 

      return
      end

!$$ Redudant rountine that worked out tensor Bn for use in calculation of polymeric stress
!$$ This Bn is that used in explicit time scheme A found in Lind (2010)
 subroutine calcBn(N,alphamax,betamax,imax,dt,gradu,Pstress,Pstressn,F,Fn,S,Sn,CG5,UG3,AG4,Bn,lambda1,Usmth)
  implicit none
 
          integer :: N, alphamax, betamax, imax, a, b, c, k, l, i, j
 double precision :: UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),& 
                     AG4(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1)),&          
                     gradu(0:N*(alphamax+1),0:N*(betamax+1),2,2),h,dt,&
                     S(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Sn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     CG5(0:N*(alphamax+1),0:N*(betamax+1),0:N*(alphamax+1),0:N*(betamax+1),2),&     
                     lambda1(0:N*(alphamax+1),0:N*(betamax+1)),&
                     ct(0:N*(alphamax+1),0:N*(betamax+1)),&
                     S1, S2, F(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Bn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Fn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Pstressn(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     graduT(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     F1(2,2), F2(2,2),&
                     Usmth(0:3,0:N*(alphamax+1),0:N*(betamax+1))
                       


   !! ct = 1d0 + lambda1/dt

!!$  do l = 0, N*(betamax + 1)
!!$   do k = 0, N*(alphamax + 1)
!!$
!!$    graduT(k,l,1,1) = gradu(k,l,1,1)
!!$    graduT(k,l,2,1) = gradu(k,l,1,2)
!!$    graduT(k,l,1,2) = gradu(k,l,2,1)
!!$    graduT(k,l,2,2) = gradu(k,l,2,2)
!!$
!!$   enddo
!!$  enddo

  do j = 0, N*(betamax + 1) 
     do i = 0, N*(alphamax + 1)   

      ct(i,j) = 1d0/(1d0 + 3d0*lambda1(i,j)/(2d0*dt)) 
 
     enddo
   enddo

  do l = 0, N*(betamax + 1) 
     do k = 0, N*(alphamax + 1)   
   
           do b = 1, 2   !! Deformation terms (gradu)P + P(gradu)^T
             do a = 1, 2
                 
               S1 = 0d0
               S2 = 0d0

              do c = 1, 2

                 S1 = S1 + gradu(k,l,a,c)*Pstress(k,l,c,b) 
                 S2 = S2 + Pstress(k,l,a,c)*gradu(k,l,b,c)

               enddo 

            F(k,l,a,b) = S1 + S2
           
              enddo
             enddo            

!!$            do a = 1, 2
!!$             do b = 1, 2
!!$               F(k,l,a,b) = S1 + S2
!!$             enddo
!!$            enddo

                 do a = 1, 2 !convection term u.grad(P)
                   do b = 1, 2
 
                     S(k,l,a,b) = 0d0
 
                     do j = 0, N*(betamax+1)
                     do i = 0, N*(alphamax+1)
                     do c = 1, 2
                    !    S(k,l,a,b) = S(k,l,a,b) - UG3(0,c,k,l)*Pstress(i,j,a,b)*CG5(k,l,i,j,c)

                         S(k,l,a,b) = S(k,l,a,b) - Usmth(c,k,l)*Pstress(i,j,a,b)*CG5(k,l,i,j,c)
                     enddo
                     enddo
                     enddo

                    if(AG4(k,l,k,l).eq.0d0)then
                     write(*,*) 'AG4 is zero'
                    stop 
                    endif

                   S(k,l,a,b) = S(k,l,a,b)/AG4(k,l,k,l)
                    
                   enddo
                  enddo

                  Bn(k,l,1,1) = F(k,l,1,1)
                  Bn(k,l,1,2) = F(k,l,1,2)
                  Bn(k,l,2,1) = F(k,l,2,1)
                  Bn(k,l,2,2) = F(k,l,2,2)
                

  !!! if(time.eq.1)then

!!$               Bn(k,l,1,1) = (lambda1/ct)*(Pstress(k,l,1,1)/dt + S(k,l,1,1) + F(k,l,1,1))
!!$               Bn(k,l,1,2) = (lambda1/ct)*(Pstress(k,l,1,2)/dt + S(k,l,1,2) + F(k,l,1,2))
!!$               Bn(k,l,2,1) = (lambda1/ct)*(Pstress(k,l,2,1)/dt + S(k,l,2,1) + F(k,l,2,1))
!!$               Bn(k,l,2,2) = (lambda1/ct)*(Pstress(k,l,2,2)/dt + S(k,l,2,2) + F(k,l,2,2))

  !!  else

!!$               Bn(k,l,1,1) = ct(k,l)*((lambda1(k,l)/(2d0*dt))*(4d0*Pstress(k,l,1,1) - Pstressn(k,l,1,1))&
!!$                           + lambda1(k,l)*(2d0*S(k,l,1,1) - Sn(k,l,1,1)) + lambda1(k,l)*F(k,l,1,1))!!(2d0*F(k,l,1,1) - Fn(k,l,1,1)))
!!$              
!!$               Bn(k,l,1,2) = ct(k,l)*((lambda1(k,l)/(2d0*dt))*(4d0*Pstress(k,l,1,2) - Pstressn(k,l,1,2))&
!!$                           + lambda1(k,l)*(2d0*S(k,l,1,2) - Sn(k,l,1,2)) + lambda1(k,l)*F(k,l,1,2))!!(2d0*F(k,l,1,2) - Fn(k,l,1,2)))
!!$              
!!$               Bn(k,l,2,1) = ct(k,l)*((lambda1(k,l)/(2d0*dt))*(4d0*Pstress(k,l,2,1) - Pstressn(k,l,2,1))&
!!$                           + lambda1(k,l)*(2d0*S(k,l,2,1) - Sn(k,l,2,1)) + lambda1(k,l)*F(k,l,2,1))!!(2d0*F(k,l,2,1) - Fn(k,l,2,1)))
!!$
!!$               Bn(k,l,2,2) = ct(k,l)*((lambda1(k,l)/(2d0*dt))*(4d0*Pstress(k,l,2,2) - Pstressn(k,l,2,2))&
!!$                           + lambda1(k,l)*(2d0*S(k,l,2,2) - Sn(k,l,2,2)) + lambda1(k,l)*F(k,l,2,2))!!(2d0*F(k,l,2,2) - Fn(k,l,2,2)))


!!$ Below does not including varying lambda. Needs to go in calcVG3.

!!$               Bn(k,l,1,1) = (1d0/(2d0*dt))*(4d0*Pstress(k,l,1,1) - Pstressn(k,l,1,1))&
!!$                           + 2d0*S(k,l,1,1) - Sn(k,l,1,1) + F(k,l,1,1)!!(2d0*F(k,l,1,1) - Fn(k,l,1,1)))
!!$              
!!$               Bn(k,l,1,2) = (1d0/(2d0*dt))*(4d0*Pstress(k,l,1,2) - Pstressn(k,l,1,2))&
!!$                           + 2d0*S(k,l,1,2) - Sn(k,l,1,2) + F(k,l,1,2)!!(2d0*F(k,l,1,2) - Fn(k,l,1,2)))
!!$              
!!$               Bn(k,l,2,1) = (1d0/(2d0*dt))*(4d0*Pstress(k,l,2,1) - Pstressn(k,l,2,1))&
!!$                           + 2d0*S(k,l,2,1) - Sn(k,l,2,1) + F(k,l,2,1)!!(2d0*F(k,l,2,1) - Fn(k,l,2,1)))
!!$
!!$               Bn(k,l,2,2) = (1d0/(2d0*dt))*(4d0*Pstress(k,l,2,2) - Pstressn(k,l,2,2))&
!!$                           + 2d0*S(k,l,2,2) - Sn(k,l,2,2) + F(k,l,2,2)!!(2d0*F(k,l,2,2) - Fn(k,l,2,2)))





   !!  endif


    enddo
    enddo

       
     do l = 0, N*(betamax + 1) 
     do k = 0, N*(alphamax + 1)  

          do a = 1, 2
           do b = 1, 2

             Pstressn(k,l,a,b) = Pstress(k,l,a,b)
                   Fn(k,l,a,b) = F(k,l,a,b)
                   Sn(k,l,a,b) = S(k,l,a,b)

            enddo
           enddo

      enddo
      enddo

return
endsubroutine


!$$ Form a local (element-wise) colour function from global version
subroutine calc_local_colour(N, alphamax, betamax, Cij, Clij) 
implicit none

         integer :: i, j, k , l, N, al ,be, alphamax, betamax, a 
double precision :: Clij(0:alphamax,0:betamax,0:N,0:N,2), Cij(0:N*(alphamax+1),0:N*(betamax+1),2)

     do a = 1, 2             
         do l = 0, N*(betamax+1)
            do k = 0, N*(alphamax+1)
               
               al = (k-1)/N   !unfortunately this does not work for N=1
               be = (l-1)/N
               i = mod(k+al,N+1)
               j = mod(l+be,N+1)

           !!    U(it,al,be,a,i,j) = UG3(it,a,k,l)
               Clij(al,be,i,j,a) = Cij(k,l,a)

            enddo
         enddo
      enddo

     
 do a = 1,2  !update from F=0 if out of iteration else from F=1
        
       do al = 1, alphamax
            do be = 0, betamax
               do j = 0 , N 
                !!  U(it,al,be,a,0,j) = U(it,al-1,be,a,N,j) !resets x-direction duplicates
                  Clij(al,be,0,j,a) = Clij(al-1,be,N,j,a)    
               enddo
            enddo
             
         enddo

         do be = 1, betamax
            do al = 0, alphamax
               do i = 0, N
               !   U(it,al,be,a,i,0) = U(it,al,be-1,a,i,N)
                  Clij(al,be,i,0,a) = Clij(al,be-1,i,N,a) 
               enddo
            enddo
               
         enddo
      enddo 

endsubroutine


!$$ Experimental routine that looked at smoothing velocity (not used - and
!$$ probably never should be!)
subroutine smoothU(imax,N,alphamax,betamax,xg,yg,xq,yq,P,W,J4,U,UG3,Usmth,Len,Height)
implicit none

          integer :: i, j, k, l, N, alphamax, betamax,  a, imax, Mhat, Nhat,&
                     al, be, alpha, beta, almax, almin, bemax, bemin, ii, jj
 
 double precision :: Int, x0, y0, x1, y1, x2, y2, dx, dy,x,y,&
                     D, re,&
                     UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                     Usmth(0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                     U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                     xg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     yg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     xq(0:alphamax,0:betamax,5),&
                     yq(0:alphamax,0:betamax,5),& 
                     W(0:N),P(0:N),&
                     J4(0:alphamax,0:betamax,0:N,0:N),&
                     Kernel, Len, Height, S


   re = 5d0   ! Width parameter in Kernel (bigger re -> thinner kernel (dirac delta func))
    D = 0.1d0  ! Support of Kernel (domain for which Kernel.ne.0) 


 Mhat = N*(alphamax + 1)
 Nhat = N*(betamax + 1)

do k = 0, Mhat  ! Loop over all points in domain
 do l = 0, Nhat

  do a = 0, 2

  Usmth(a,k,l) = UG3(0,a,k,l)

  enddo

     alpha = (k-1)/N   !unfortunately this does not work for N=1
     beta  = (l-1)/N
   
if(alpha.ge.1.and.alpha.le.(alphamax-1).and.beta.ge.0.and.beta.le.(betamax-1))then

   !! ii = mod(k+alpha,N+1)  !(alpha,beta,i,j) local coordinates of point (k,l) to be smoothed
      !! jj = mod(l+beta,N+1)

    x0 = xg(k,l)!!x(alpha, alphamax, beta, betamax, P(ii), P(jj), xq) !(x,y) coordinates of smoothed point
    y0 = yg(k,l)!!y(alpha, alphamax, beta, betamax, P(ii), P(jj), xq)

!! if(((x0 + D).lt.Len).and.((x0 - D).gt.0d0).and.&
!!    ((y0 - D).gt.0d0).and.((y0 + D).lt.Height))then

  if(alpha.eq.0)then
    almin = 0
    almax = alpha + 1
   elseif(alpha.eq.alphamax)then
    almin = alpha - 1
    almax = alphamax
   else
    almin = alpha - 1
    almax = alpha + 1
  endif

  if(beta.eq.0)then   !! If on edge elements, then dont integrate beyond domain
    bemin = 0
    bemax = beta + 1
   elseif(beta.eq.betamax)then
    bemin = beta - 1
    bemax = betamax
   else
    bemin = beta - 1
    bemax = beta + 1
  endif

 do a = 0,2 ! Loop over velocity components (start from 0 to include density)

    S = 0d0

  do al = almin, almax !Loop over surrounding elements
   do be = bemin, bemax

     do i = 0, N   !Integrate over element (al,be)
      do j = 0, N 

       x1 = x(al, alphamax, be, betamax, P(i), P(j), xq)
       y1 = y(al, alphamax, be, betamax, P(i), P(j), yq)

        S = S + Kernel(x1,x0,y1,y0,re,D)*U(0,al,be,a,i,j)*W(i)*W(j)*J4(al,be,i,j)

      enddo !(ij)loop
     enddo

    enddo  !(al,be)loop
   enddo

    UG3(0,a,k,l) = S
    Usmth(a,k,l) = S

 enddo !(a) loop

!! endif

endif

 enddo !(k,l) loop
enddo
 

endsubroutine

!$$ Experimental routine that looked at smoothing density (not used - and
!$$ probably never should be!)
subroutine smooth_density(imax,N,alphamax,betamax,xg,yg,xq,yq,P,W,J4,U,Usmth,Len,Height)
implicit none

          integer :: i, j, k, l, N, alphamax, betamax,  a, imax, Mhat, Nhat,&
                     al, be, alpha, beta, almax, almin, bemax, bemin, ii, jj
 
 double precision :: Int, x0, y0, x1, y1, x2, y2, dx, dy,x,y,&
                     D, re,&
                     Usmth(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                     U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                     xg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     yg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     xq(0:alphamax,0:betamax,5),&
                     yq(0:alphamax,0:betamax,5),& 
                     W(0:N),P(0:N),&
                     J4(0:alphamax,0:betamax,0:N,0:N),&
                     Kernel, Len, Height, S


   re = 20d0   ! Width parameter in Kernel (bigger re -> thinner kernel (dirac delta func))
    D = 0.5d0  ! Support of Kernel (domain for which Kernel.ne.0) 


 Mhat = N*(alphamax + 1)
 Nhat = N*(betamax + 1)

 Usmth = U

do alpha = 0, alphamax  ! Loop over all points in domain
 do beta = 0, betamax
  do k = 0, N
   do l = 0, N

   
   
if(alpha.ge.1.and.alpha.le.(alphamax-1).and.beta.ge.0.and.beta.le.(betamax-1))then

   !! ii = mod(k+alpha,N+1)  !(alpha,beta,i,j) local coordinates of point (k,l) to be smoothed
      !! jj = mod(l+beta,N+1)

    x0 = x(alpha, alphamax, beta, betamax, P(k), P(l), xq) !(x,y) coordinates of smoothed point
    y0 = y(alpha, alphamax, beta, betamax, P(k), P(l), xq)

!! if(((x0 + D).lt.Len).and.((x0 - D).gt.0d0).and.&
!!    ((y0 - D).gt.0d0).and.((y0 + D).lt.Height))then

  if(alpha.eq.0)then
    almin = 0
    almax = alpha + 1
   elseif(alpha.eq.alphamax)then
    almin = alpha - 1
    almax = alphamax
   else
    almin = alpha - 1
    almax = alpha + 1
  endif

  if(beta.eq.0)then   !! If on edge elements, then dont integrate beyond domain
    bemin = 0
    bemax = beta + 1
   elseif(beta.eq.betamax)then
    bemin = beta - 1
    bemax = betamax
   else
    bemin = beta - 1
    bemax = beta + 1
  endif

 do a = 0,2 ! Loop over velocity components (start from 0 to include density)

    S = 0d0

  do al = almin, almax !Loop over surrounding elements
   do be = bemin, bemax

     do i = 0, N   !Integrate over element (al,be)
      do j = 0, N 

       x1 = x(al, alphamax, be, betamax, P(i), P(j), xq)
       y1 = y(al, alphamax, be, betamax, P(i), P(j), yq)

        S = S + Kernel(x1,x0,y1,y0,re,D)*U(0,al,be,a,i,j)*W(i)*W(j)*J4(al,be,i,j)

      enddo !(ij)loop
     enddo

    enddo  !(al,be)loop
   enddo

    Usmth(0,alpha,beta,a,k,l) = S

 enddo !(a) loop

!! endif

endif


   enddo !(k,l) loop
  enddo
 enddo
enddo 

endsubroutine


!$$ Experimental routine that looked at smoothing stresses (not used - and
!$$ probably never should be!)
subroutine smooth_tensor(imax,N,alphamax,betamax,xg,yg,xq,yq,P,W,J4,TG,Len,Height)
implicit none

          integer :: i, j, k, l, N, alphamax, betamax,  a, b, imax, Mhat, Nhat,&
                     al, be, alpha, beta, almax, almin, bemax, bemin, ii, jj
 
 double precision :: Int, x0, y0, x1, y1, x2, y2, dx, dy,x,y,&
                     D, re,&
                     TG(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     Tsmth(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                     T_loc(0:alphamax,0:betamax,0:N,0:N,2,2),&
                     xg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     yg(0:N*(alphamax+1),0:N*(betamax+1)),&
                     xq(0:alphamax,0:betamax,5),&
                     yq(0:alphamax,0:betamax,5),& 
                     W(0:N),P(0:N),&
                     J4(0:alphamax,0:betamax,0:N,0:N),&
                     Kernel, Len, Height, S


   re = 20d0   ! Width parameter in Kernel (bigger re -> thinner kernel (dirac delta func))
    D = 0.5d0  ! Support of Kernel (domain for which Kernel.ne.0) 


 Mhat = N*(alphamax + 1)
 Nhat = N*(betamax + 1)

    do a = 1, 2  
     do b = 1, 2           
         do l = 0, N*(betamax+1)
            do k = 0, N*(alphamax+1)
               
               al = (k-1)/N   !unfortunately this does not work for N=1
               be = (l-1)/N
               i = mod(k+al,N+1)
               j = mod(l+be,N+1)

           !!    U(it,al,be,a,i,j) = UG3(it,a,k,l)
               T_loc(al,be,i,j,a,b) = TG(k,l,a,b)
                     Tsmth(k,l,a,b) = TG(k,l,a,b)

            enddo
         enddo
      enddo
    enddo
     
 do a = 1,2  !update from F=0 if out of iteration else from F=1
  do b = 1, 2   
     
       do al = 1, alphamax
            do be = 0, betamax
               do j = 0 , N 
                !!  U(it,al,be,a,0,j) = U(it,al-1,be,a,N,j) !resets x-direction duplicates
                  T_loc(al,be,0,j,a,b) = T_loc(al-1,be,N,j,a,b)    
               enddo
            enddo
             
         enddo

         do be = 1, betamax
            do al = 0, alphamax
               do i = 0, N
               !   U(it,al,be,a,i,0) = U(it,al,be-1,a,i,N)
                  T_loc(al,be,i,0,a,b) = T_loc(al,be-1,i,N,a,b) 
               enddo
            enddo
               
         enddo
      enddo 
   enddo

do k = 0, Mhat  ! Loop over all points in domain
 do l = 0, Nhat


     alpha = (k-1)/N   !unfortunately this does not work for N=1
     beta  = (l-1)/N
   
if(beta.eq.0.or.(beta.le.(betamax-1).and.alpha.ge.1.and.alpha.le.(alphamax-1)))then

!!if(alpha.ge.1.and.alpha.le.(alphamax-1).and.beta.ge.0.and.beta.le.(betamax-1))then

   !! ii = mod(k+alpha,N+1)  !(alpha,beta,i,j) local coordinates of point (k,l) to be smoothed
      !! jj = mod(l+beta,N+1)

    x0 = xg(k,l)!!x(alpha, alphamax, beta, betamax, P(ii), P(jj), xq) !(x,y) coordinates of smoothed point
    y0 = yg(k,l)!!y(alpha, alphamax, beta, betamax, P(ii), P(jj), xq)

!! if(((x0 + D).lt.Len).and.((x0 - D).gt.0d0).and.&
!!    ((y0 - D).gt.0d0).and.((y0 + D).lt.Height))then

  if(alpha.eq.0)then
    almin = 0
    almax = alpha + 1
   elseif(alpha.eq.alphamax)then
    almin = alpha - 1
    almax = alphamax
   else
    almin = alpha - 1
    almax = alpha + 1
  endif

  if(beta.eq.0)then   !! If on edge elements, then dont integrate beyond domain
    bemin = 0
    bemax = beta + 1
   elseif(beta.eq.betamax)then
    bemin = beta - 1
    bemax = betamax
   else
    bemin = beta - 1
    bemax = beta + 1
  endif

 do a = 1,2 ! Loop over velocity components (start from 0 to include density)
  do b = 1,2

    S = 0d0

  do al = almin, almax !Loop over surrounding elements
   do be = bemin, bemax

     do i = 0, N   !Integrate over element (al,be)
      do j = 0, N 

       x1 = x(al, alphamax, be, betamax, P(i), P(j), xq)
       y1 = y(al, alphamax, be, betamax, P(i), P(j), yq)

        S = S + Kernel(x1,x0,y1,y0,re,D)*T_loc(al,be,i,j,a,b)*W(i)*W(j)*J4(al,be,i,j)

      enddo !(ij)loop
     enddo

    enddo  !(al,be)loop
   enddo

   Tsmth(k,l,a,b) = S

 enddo !b loop
 enddo !(a) loop

 endif

!!endif

  enddo !(k,l) loop
enddo

  TG = Tsmth 
 
endsubroutine


!$$ Gaussian kernel function
!$$ For use in above smoothing routines
double precision function Kernel(x1,x0,y1,y0,re,D)
implicit none

double precision :: x1, y1, x0, y0, re, pi, D, rsqd 

    pi = 4d0*datan(1d0)

  rsqd = ((x1 - x0)**2 + (y1 - y0)**2)
 
if(rsqd.lt.D**2)then 

 Kernel = (re/pi)*exp(-re*rsqd) 

else

 Kernel = 0d0

endif

endfunction


!$$ Routine that writes global tensor data to a local version
subroutine global_2_local_tensor(N,alphamax,betamax,TG,TL)
implicit none

         integer :: N, alphamax, betamax, i,j,k,l,a,b,al,be

double precision :: TG(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                    TL(0:alphamax,0:betamax,0:N,0:N,2,2)

          do l = 0, N*(betamax+1)
            do k = 0, N*(alphamax+1)
               
               al = (k-1)/N   !unfortunately this does not work for N=1
               be = (l-1)/N
               i = mod(k+al,N+1)
               j = mod(l+be,N+1)

           !!    U(it,al,be,a,i,j) = UG3(it,a,k,l)

              do a = 1, 2
               do b = 1, 2
              
                  TL(al,be,i,j,a,b) = TG(k,l,a,b)
                   
               enddo
              enddo   

            enddo
         enddo
    
     
 do a = 1,2  !update from F=0 if out of iteration else from F=1
  do b = 1, 2   
     
       do al = 1, alphamax
            do be = 0, betamax
               do j = 0 , N 
                !!  U(it,al,be,a,0,j) = U(it,al-1,be,a,N,j) !resets x-direction duplicates
                  TL(al,be,0,j,a,b) = TL(al-1,be,N,j,a,b)    
               enddo
            enddo
        enddo

         do be = 1, betamax
            do al = 0, alphamax
               do i = 0, N
               !   U(it,al,be,a,i,0) = U(it,al,be-1,a,i,N)
                  TL(al,be,i,0,a,b) = TL(al,be-1,i,N,a,b) 
               enddo
            enddo
         enddo
      
     enddo 
   enddo

endsubroutine


 !$$ Decomposes stress into a normal and shear component to the flow as 
 !$$ suggested by Bollada and Phillips (2008).
 !$$ Plays no critical role - just for interest
 subroutine calcStressDecomp(imax,N,alphamax,betamax,TG4,Pstress,UG3,&
                             Norm_Pstress,Shear_Pstress,Norm_TG4,Shear_TG4)

 implicit none

  integer :: imax, N, alphamax, betamax, k , l,a, b
  double precision ::  TG4(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                       Pstress(0:N*(alphamax+1),0:N*(betamax+1),2,2),&
                       UG3(0:imax,0:3,0:N*(alphamax+1),0:N*(betamax+1)),&
                       Norm_Pstress(0:N*(alphamax+1),0:N*(betamax+1)),&
                       Shear_Pstress(0:N*(alphamax+1),0:N*(betamax+1)),&
                       Norm_TG4(0:N*(alphamax+1),0:N*(betamax+1)),&
                       Shear_TG4(0:N*(alphamax+1),0:N*(betamax+1)),&
                       U(2),V(2), modU, modV, crossp
                       

  !! V is perpendicular to U - ie. U.V=0

   do k = 0, N*(alphamax+1)
    do l = 0, N*(betamax+1)

     !Calculates U 

     if(k.eq.0)then

      U(1) = 1d0 
      U(2) = 0d0

    elseif(l.eq.0)then

      U(1) = 0d0
      U(2) = 1d0

    elseif(k.eq.(N*(alphamax+1)))then

      U(1) = -1d0
      U(2) = 0d0

     elseif(l.eq.(N*(betamax+1)))then

      U(1) = 0d0
      U(2) = -1d0

     else

      U(1) = UG3(0,1,k,l)     
      U(2) = UG3(0,2,k,l)

      modU = dsqrt(U(1)**2 + U(2)**2)

     if(modU.eq.0d0)then

         U = 0d0
      
     else

         U = U/modU

      endif

     endif

      V(1) = -U(2)
      V(2) = U(1)

      modV = dsqrt(V(1)**2 + V(2)**2)

    if(modV.eq.0d0)then
      
        V = 0d0

    else

        V = V/modV

    endif

     ! Calculates V

   Norm_Pstress(k,l) = 0d0
  Shear_Pstress(k,l) = 0d0
       Norm_TG4(k,l) = 0d0
      Shear_TG4(k,l) = 0d0

  do a = 1,2
   do b = 1,2

    Norm_Pstress(k,l) = Norm_Pstress(k,l) + U(a)*Pstress(k,l,a,b)*U(b)
   Shear_Pstress(k,l) = Shear_Pstress(k,l) + V(a)*Pstress(k,l,a,b)*U(b)

        Norm_TG4(k,l) = Norm_TG4(k,l) + U(a)*TG4(k,l,a,b)*U(b)
       Shear_TG4(k,l) = Shear_TG4(k,l) + V(a)*TG4(k,l,a,b)*U(b)

   enddo
  enddo


    enddo
   enddo


 endsubroutine

!$$ Very crude means of calculating centre of bubble and the 
!$$ density at that point. 

subroutine calc_bubble_xy(N,imax,alphamax,betamax,TNP,Np1,P,xp,yp,xq,yq,U,Cp,&
                          x_bub_c,y_bub_c,Q_bub_c)
implicit none

         integer :: l, pp, TNP, Np1,N,alphamax,betamax,alpha,beta,imax
double precision :: xp(TNP), yp(TNP),&
                    x_bub(Np1), y_bub(Np1), Cp(TNP,2),&
                    x_bub_min, x_bub_max, y_bub_min, y_bub_max,&
                    U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                    xi,zeta,&
                    xq(0:alphamax, 0:betamax,5),&
                    yq(0:alphamax, 0:betamax,5),eps,&
                    x_bub_c, y_bub_c, Q_bub_c,&  
                    S(0:3),P(0:N)
                    
     eps = 1d-12            
 x_bub_c = 0d0
 y_bub_c = 0d0
  
       l = 1

do pp = 1, TNP

 
  if(Cp(pp,1).eq.1d0)then

   x_bub(l) = xp(pp)
   y_bub(l) = yp(pp) 

   l = l + 1

  endif

enddo

 x_bub_min = minval(x_bub)
 x_bub_max = maxval(x_bub)

 y_bub_min = minval(y_bub)
 y_bub_max = maxval(y_bub)

 x_bub_c = x_bub_min + (x_bub_max - x_bub_min)/2d0 ! x pos of bubble centre
 y_bub_c = y_bub_min + (y_bub_max - y_bub_min)/2d0 ! y pos of bubble centre

 
  !!$ For a given (x,y), calculate the surrounding element (alpha,beta) and position (xi,zeta) within that element
    call x2xi(alphamax,betamax,alpha,beta, x_bub_c, y_bub_c, xq, yq, xi, zeta, eps)
            
     !!$ Given this postion - calculate the velocity Un
    call calcUxy(U,alpha,beta,xi,zeta,0,N,P,alphamax,betamax,S,imax,0d0)!also includes calculation for Qn
           
       Q_bub_c = S(0)
        
Write(*,*) 'bubble centre is',x_bub_c,y_bub_c
write(*,*) 'density at this point=',Q_bub_c

return
endsubroutine

!$$ Given the present values of Reynolds and Deborah number, this 
!$$ routine calculates the equivalent values using the boundary element (BEM) 
!$$ definitions (see Lind 2010)
!$$ Plays no critical role - allows for comparisons.

 subroutine calc_BEM_parameters(c2,rho0,Q_bub_c,lambda,mu_SEM,De_BEM,Re_BEM,mu_BEM)
 implicit none

          
 double precision :: c2, c, rhof, rhob, rho0(2), Q_bub_c,&
                     drho, lambda, De_BEM, Re_BEM, mu_BEM, mu_SEM

 
    c = dsqrt(c2)

 rhof = exp(rho0(2))

 rhob = exp(Q_bub_c) 

 drho = rhof - rhob

 if(drho.eq.0d0)then

  De_BEM = 0d0
  Re_BEM = 0d0

 else

     De_BEM = lambda*c*dsqrt(abs(drho)/rhof) !Deborah number (BEM defn)

     mu_BEM = rhof*mu_SEM/(c*dsqrt(abs(drho)*rhof)) !Normalised viscosity (BEM defn)

     Re_BEM = 1d0/mu_BEM ! Reynolds number (BEM defn)

  endif

 
 return
 endsubroutine

!$$ Calcs bubble mass (for testing purposes)
subroutine calc_bubble_mass(N,imax,alphamax,betamax,W,J4,U,Mass,Mass_error,rad,Clij,rho_b0,time)
implicit none

integer          :: alphamax, betamax, imax, al, be, i, j, N, time

double precision :: W(0:N),&
                    J4(0:alphamax,0:betamax,0:N,0:N),&
                    U(0:imax,0:alphamax,0:betamax,0:3,0:N,0:N),&
                    Mass, rad, Mass_error, Mass_anal, pi,&
                    Clij(0:alphamax,0:betamax,0:N,0:N,2),&
                    rho_b, rho_b0, rho_l0, bub
Mass = 0d0
  pi = 4d0*datan(1d0) 

! Analytical (initial) mass 
 Mass_anal = exp(rho_b0)*pi*rad**2

 write(*,*) 'Mass_anal',Mass_anal
 
 do al = 0, alphamax
  do be = 0, betamax
   do i = 0, N
    do j = 0, N

       bub = Clij(al,be,i,j,1)

     if(bub.gt.0.5d0)then

     rho_b = exp(U(0,al,be,0,i,j))

      Mass = Mass + rho_b*J4(al,be,i,j)*W(i)*W(j)
  
     endif

    enddo
   enddo
  enddo
 enddo

!$$ Compares initial mass with subsequent masses 
!$$ calculated using the bubble density and volume.
 Mass_error = abs(Mass - Mass_anal)/Mass_anal

return
endsubroutine

!$$ Routine calculates bubble volume and subsequently
!$$ the (equivalent) bubble radius.
subroutine calc_bubble_vol(N,alphamax,betamax,W,J4,Clij,vol,rad_av)
implicit none

integer          :: alphamax, betamax, imax, al, be, i, j, N

double precision :: W(0:N),&
                    J4(0:alphamax,0:betamax,0:N,0:N),&
                    Clij(0:alphamax,0:betamax,0:N,0:N,2), vol, rad_av, pi

   vol = 0d0
 
    pi = 4d0*datan(1d0)
 
 do al = 0, alphamax
  do be = 0, betamax
   do i = 0, N
    do j = 0, N

      vol = vol + Clij(al,be,i,j,1)*J4(al,be,i,j)*W(i)*W(j)
  
     enddo
   enddo
  enddo
 enddo

   rad_av = dsqrt(vol/pi)

return
endsubroutine

!$$ Necessary function call that releases the memory 
!$$ held by PARDISO rounties

subroutine Pardiso_Release_Memory(pt,Gmax)
implicit none

 integer*8  pt(64)
 integer :: maxfct, mnum, mtype, phase, nrhs, error, msglvl, i
 integer :: solver, idum, Gmax, Nonz
 integer :: iparm(64)
 
 double precision :: dparm(64), ddum
    
  
    nrhs = 1 ! number of rhs vectors
    mnum = 1 ! number of matrices
  maxfct = 1 ! max number of matrices (i think)

   mtype = 11 !real unsymmetric matrix 
 ! mtype = 2 for real symmetric pos def matrix

   phase = -1 !release memory

  msglvl = 0 !print no (=0) statistical information on screen 

  iparm(1) = 0 !use default solver settings
  iparm(3) = 1

  idum = 0 !dummy variables which aren't used for this particular phase
  ddum = 0

   call pardiso (pt, maxfct, mnum, mtype, phase, Gmax, ddum, idum, idum,& 
                 idum, nrhs, iparm, msglvl, ddum, ddum, error) 


   write(*,*) 'error from pardiso',error

 endsubroutine

!$$ That's all folks!
