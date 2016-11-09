module o_param
  implicit none
  save
  !    
  ! *** Fixed parameters ***
  real(kind=8), parameter  	:: pi=3.141592653589793, rad=pi/180.0  
  real(kind=8), parameter  	:: omega=2.0*pi/(24.0*60.0*60.0)
  real(kind=8), parameter  	:: g=9.81                       ![m/s^2]
  real(kind=8), parameter  	:: r_earth=6.3675e6             ![m]
  real(kind=8), parameter  	:: rho0=1030.                   ![kg/m^3]
  real(kind=8), parameter 	:: rho0r=1.0/rho0 
  real(kind=8), parameter  	:: vcpw=4.2e6                   ![J/m^3/K]volum. heat cap. of water
  real(kind=8), parameter	:: small=1.0e-8                 !small value

  ! *** mixing and friction setting ***
  real             	        :: Ah0=6000.            	!Lapl. hori. visc, [m^2/s]
  real             	        :: Ahb0=2.7e13             	!Bihar. hori. visc,[m^4/s] <8.0e12/1degree
  real             	        :: Kh0=600.                     !lateral diff 
  !
  logical                  	:: biharmonic_visc=.false.      !bihar. true will turn off Lapl.
  logical                  	:: smagorinsky_visc=.false. 
  !
  logical                       :: increase_equ_zonal_visc=.true.!increase zonal viscosity at equator
  real                          :: fac_visc_increase=3.0         !increase factor
  !
  logical                  	:: scale_mixing_h=.true.   	!scale hor. diff/visc. coef. by resol.
  integer			:: scale_mixing_type=2		!1:scale by \delt x^2; 2:scale by \delta x
  real(kind=8)             	:: scalevol=5.e9          	!associated with 100km resolution
  !
  logical	  	   	:: Redi_GM=.true.         	!flag for Redi/GM scheme
  logical                       :: ODM95=.true.                 !taper fcn ODM95
  logical                       :: LDD97=.true.                 !taper fcn LDD97
  real(kind=8)		   	:: ratio_K_GM=1.0          	!ratio of K_GM to Kh
  real(kind=8)                  :: Slope_critical=4.0e-3        !upper slope limit for applying Redi/GM
  integer                       :: nslope_version=1             !1-neutral slope over prism,2-over tetrahedra
  !
  real             	        :: Av0=1.e-4      	        !background (or internal wave) vert. mixing
  real              	        :: Kv0=1.e-5                    !m^2/s
  real                          :: visc_conv_limit=0.1          !visc due to convective instability
  real                          :: diff_conv_limit=0.1          !diff due to convective instability
  !
  character(5)                 	:: mix_scheme='KPP'		!'KPP','PP', 'MY2p5', 'no'
  !
  real                          :: visc_sh_limit=5.0e-3         !for kpp,max visc due to shear instability
  real                          :: diff_sh_limit=5.0e-3         !for kpp,max diff due to shear instability
  logical                       :: double_diffusion=.true.      !for KPP,dd switch
  logical                       :: smooth_blmc=.true.           !for KPP,hori. smooth of b.l. mixing coeff.
  !
  real             	        :: PP_max_mix_coeff=5.0e-3     	!for PP, max Kv/Av 
  real                          :: wndmix=1.0e-3                !for PP, to simulate missing high frequency wind
  logical                       :: allow_convect_global=.true.  !for PP, convection for global or only NH
  !
  logical                  	:: add_TB04_to_PP=.false.   	!for PP, TB04 switch
  real(kind=8)             	:: modiff=0.01                  !for PP, vert. mix. coeff. for TB04
  !
  logical			:: tidal_mixing=.false.		!switch for tidal mixing
  logical			:: use_drag_dissipation=.true. 	!barotropic
  logical			:: use_wave_dissipation=.false.	!baroclinic
  logical			:: read_tide_speed=.true.	!read tide speed or use default
  real(kind=8)                  :: default_tide_speed=0.01      !(m/s)
  real(kind=8)                  :: max_drag_diffusivity=5.e-3   !m2/s
  real(kind=8)                  :: max_wave_diffusivity=5.e-3   !m2/s
  character(2) 	                :: Tmix_tidalconstituent='M2'   !which tidal constituent 
  character(15)	                :: Tmix_tidalmodelname='tpxo71' !source model name
  !
  real(kind=8)             	:: C_d=0.0025               	!Bottom fri. coeff.

  namelist /viscdiff/ Ah0, Ahb0, Kh0, Av0, Kv0, &
       biharmonic_visc, smagorinsky_visc, &
       increase_equ_zonal_visc, fac_visc_increase, &
       scale_mixing_h, scale_mixing_type, scalevol, &
       Redi_GM, ODM95, LDD97, ratio_K_GM, Slope_critical, nslope_version, &
       mix_scheme, visc_sh_limit, diff_sh_limit, visc_conv_limit, diff_conv_limit, &
       double_diffusion, smooth_blmc, PP_max_mix_coeff, wndmix, &
       allow_convect_global, add_TB04_to_PP, modiff, &
       tidal_mixing, use_drag_dissipation, &
       use_wave_dissipation, read_tide_speed, &
       default_tide_speed, max_drag_diffusivity, max_wave_diffusivity, &
       Tmix_tidalconstituent, Tmix_tidalmodelname, &
       C_d

  ! *** surface and open boundary setting ***
  logical                  	:: ts_surfbd=.true.     
  !
  logical                       :: ref_sss_local=.false.	!virtual salt flux using local SSS or ref_sss
  real(kind=8)             	:: ref_sss=34.7			!ref. value for cal. virtual salt flux
  !
  real(kind=8)             	:: restore_s_surf=10./(180.*86400.)	! m/s
  real(kind=8)			:: restore_t_surf=0.0
  !
  logical                       :: balance_salt_water=.true.    !balance virtual salt or water flux or not
  !
  logical                       :: buffer_zone=.false.
  real(kind=8)             	:: restore_ts_buff= 1./(86400.*5.)     	! timescale for buffer zone [1/s]

  namelist /boundary/ ts_surfbd, ref_sss_local, ref_sss, restore_s_surf, &
       restore_t_surf, balance_salt_water, buffer_zone, restore_ts_buff

  ! *** numerical schemes
  real(KIND=8)             	:: gamma_stab=0.99          	!stabilization for ssh
  real(KIND=8)             	:: gamma_stab_nh=0.5       	!stabilization for nh-pressure
  real(kind=8)             	:: gamma_fct=0.4           	!param for tracer FCT scheme
  real(kind=8)             	:: alpha_AB=1.55            	!when Adams-Bashforth Coriolis
  real(kind=8)             	:: alpha_trapez=0.55        	!when semi-implicit Coriolis
  real(kind=8)			:: theta_ssh=0.5		!semi-implicit ssh when semiimpl sheme	
  real(kind=8)			:: theta_vel=0.5		!semi-implicit baro. vel
  !
  logical                       :: use_vertvisc_impl=.true.     !if implicit vertical viscosity,keep true
  logical                       :: use_vertdiff_impl=.true.     !if implicit vertical diff., keep ture
  logical                       :: use_cori_semi=.false.        !if semiimplicit coriolis force
  !
  logical                  	:: lump_uv_matrix=.true.   	!for mass vel. matrix case, keep true!
  logical                  	:: lump_ts_matrix=.true.   	!for mass T/S matrix case, keep true!
  integer                  	:: num_iter_solve=3        	!iteration # for mass matrix case

  namelist /oce_scheme/ gamma_stab, gamma_stab_nh, gamma_fct, alpha_AB, alpha_trapez, &
       theta_ssh, theta_vel, use_cori_semi

  ! *** density and pressure force ***
  logical                  	:: density_linear=.false.
  logical                 	:: use_ref_density=.true.     

  namelist /denspress/ density_linear, use_ref_density

  ! *** parameters for nonlinear free surface cases ***
  real(kind=8)             	:: max_ice_loading=5.0		!m, maximal pressure from ice felt by the ocean

  namelist /param_freesurf/ max_ice_loading

  ! *** tide configuration ***
  integer     	   	        :: nmbr_tidal_cons=4
  character(20)	                :: tidal_constituent='M2S2K1O1' !M2 S2 N2 K2 K1 O1 P1 Q1
  character(15)			:: tidemodelname='tpxo71'
  character(10)                 :: tide_opbnd_type='Flather' 	!ssh, vel, Flather 
  real(kind=8)                  :: tide_amplify_coeff=1.0       !amplify tidal amplitude by this factor  

  namelist /tide_obc/ nmbr_tidal_cons, tidal_constituent, tidemodelname, tide_opbnd_type, tide_amplify_coeff
  
  ! *** passive tracer ***
  logical                       :: use_passive_tracer=.false.    !use passive tracer
  integer                       :: num_passive_tracer=1          !only 1 before update
  integer                       :: ptr_start_year=1948           !when to start having ptr
  logical                       :: passive_tracer_flux=.false.   !ptr enters from sfc flux
  logical                       :: passive_tracer_restore=.true. !add ptr from restoring 
  logical                       :: ptr_restore_in_volume=.true.  !restoring in the 3D region
  real(kind=8)                  :: ptr_background_value=0.0      !ptr init. background value
  real(kind=8)                  :: ptr_restore_value=1.0         !restore value for ptr
 
  namelist /passive_tracer/ use_passive_tracer, num_passive_tracer, &
       ptr_start_year, passive_tracer_flux, passive_tracer_restore, &
       ptr_restore_in_volume, ptr_background_value, ptr_restore_value
  
  ! *** age tracer ***
  logical                       :: use_age_tracer=.false.        !use age tracer
  integer                       :: num_age_tracer=1              !only 1 before update
  integer                       :: age_tracer_start_year=1948    !when to start having age tr.
  logical                       :: zero_age_at_release=.true.    !keep zero age in rel. zone
  logical                       :: age_release_in_volume=.false. !age tr. release in 3D regions
  real(kind=8)                  :: age_tracer_restore_time=864000.!restore time scale in rel. zone

  namelist /age_tracer/ use_age_tracer, age_release_in_volume, zero_age_at_release, &
       age_tracer_restore_time, num_age_tracer, age_tracer_start_year       

  ! *** others ***
  integer                       :: num_tracer

end module o_param
  !  
  !no user specification beyond this line!
  !---------------------------------------------------------------------------------------
  !
module o_array
  implicit none
  save
  
  character(4), allocatable, dimension(:)         :: prog_tracer_name

  real(kind=8), allocatable, dimension(:)         :: coriolis_param_elem2D
  real(kind=8), allocatable, dimension(:)         :: coriolis_param_nod2D  
  real(kind=8), allocatable, dimension(:)    	  :: stress_x
  real(kind=8), allocatable, dimension(:)    	  :: stress_y
  real(kind=8), allocatable, dimension(:,:)  	  :: stress_x_t
  real(kind=8), allocatable, dimension(:,:)  	  :: stress_y_t
  real(kind=8), allocatable, dimension(:)    	  :: Tsurf, Ssurf
  real(kind=8), allocatable, dimension(:)  	  :: heat_flux
  real(kind=8), allocatable, dimension(:,:)       :: heat_flux_t  
  real(kind=8), allocatable, dimension(:)         :: water_flux
  real(kind=8), allocatable, dimension(:,:)       :: water_flux_t         

  real(kind=8), allocatable, target, dimension(:) :: uv_rhs, uf, uf0, duf 
  real(kind=8), allocatable, target, dimension(:) :: ssh, ssh0, dssh, ssh_rhs
#ifndef use_non_hydrostatic
  real(kind=8), allocatable, target, dimension(:) :: wrhs, w  
#else
  real(kind=8), allocatable, target, dimension(:) :: nhp_rhs, nhp, nhp0 
#endif
  real(kind=8), allocatable, target, dimension(:,:) :: tracer_rhs, tracer, dtracer
  
  real(kind=8), allocatable, target, dimension(:,:) :: ts_sfc_force
  real(kind=8), allocatable, target, dimension(:,:) :: uv_sfc_force, uv_bott_force

  real(kind=8), allocatable, dimension(:)         :: virtual_salt, relax_salt
  real(kind=8), allocatable, dimension(:)         :: real_salt_flux

  real(kind=8), allocatable, target, dimension(:) :: ucori, vcori, ucori_back, vcori_back
 
  real(kind=8), allocatable, dimension(:)         :: density_ref, density_insitu  
  real(kind=8), allocatable, dimension(:)         :: hpressure
  real(kind=8), allocatable, dimension(:,:,:)     :: PGF

  ! buffer zone restoring 
  real(kind=8), allocatable, target, dimension(:,:) :: tracer0
  real(kind=8), allocatable, dimension(:)         :: tracer_restore_coeff

  ! restoring open boundary
  real(kind=8), allocatable, dimension(:)         :: opbnd_ssh_rhs

  ! tidal open boundary
  real(kind=8), allocatable, dimension(:)         :: opbnd_u_tide, opbnd_v_tide
  real(kind=8), allocatable, dimension(:)         :: opbnd_z_tide, opbnd_z0_tide
  real(kind=8), allocatable, dimension(:)    	  :: tide_period_coeff
  real(kind=8), allocatable, dimension(:,:)	  :: tide_u_amp, tide_u_pha
  real(kind=8), allocatable, dimension(:,:)	  :: tide_v_amp, tide_v_pha
  real(kind=8), allocatable, dimension(:,:)	  :: tide_z_amp, tide_z_pha
  real(kind=8), allocatable, dimension(:)         :: opbnd_dep

  ! vertical mixing 
  real, allocatable, dimension(:) 	  :: Av
  real, allocatable, dimension(:,:)  :: Kv 

  !FCT advection scheme
  real(kind=8), allocatable, dimension(:,:)       :: tral   
  real(kind=8), allocatable, dimension(:,:)       :: trafluxes 
  real(kind=8), allocatable, dimension(:)         :: pplus, pminus

  !Redi/GM flag
  real(kind=8), allocatable, dimension(:,:,:)  	  :: neutral_slope 
  real(kind=8), allocatable, dimension(:,:)  	  :: neutral_slope_elem

end module o_array
!
!----------------------------------------------------------------------
!
module o_solver
  implicit none
  save
  !solver index
  integer                  	:: solve_u=1
  integer       	 	:: solve_v=2
  integer                  	:: solve_w=6
  integer                  	:: solve_ssh=7
  integer                  	:: solve_nhp=8
  integer                  	:: solve_tra=10
  logical                  	:: iter_first=.true.
  logical                  	:: iteruv_first=.true.
end module o_solver
!
!----------------------------------------------------------------------
!
module g_parfe
  implicit none
  save 
  !
#ifdef PETSC
#include "finclude/petsc.h"
#else
  include 'mpif.h'
#endif

  ! communication part
  type com_struct
     integer    :: rPEnum
     integer, dimension(:), pointer :: rPE
     integer, dimension(:), pointer :: rptr
     integer, dimension(:), pointer :: rlist
     integer    :: sPEnum
     integer, dimension(:), pointer :: sPE
     integer, dimension(:), pointer :: sptr
     integer, dimension(:), pointer :: slist
  end type com_struct
  type(com_struct)                      	:: com_nod2D, com_nod3D
  type com_array
     real(kind=8), dimension(:), pointer :: array
  end  type com_array
  type(com_array), allocatable          	:: s_buff_2d(:), r_buff_2d(:)
  type(com_array), allocatable          	:: s_buff_3d(:), r_buff_3d(:)

  ! general MPI part
  integer                               	:: MPIERR
  integer           				:: npes
  integer        				:: mype
  integer, allocatable, dimension(:)  		:: part2D, part3D       

  ! Mesh partition
  integer                             		:: myDim_nod2D, eDim_nod2D, ToDim_nod2D
  integer, allocatable, dimension(:)  		:: myList_nod2D
  integer                             		:: myDim_nod3D, eDim_nod3D, ToDim_nod3D
  integer, allocatable, dimension(:)  		:: myList_nod3D
  integer                             		:: myDim_elem2D
  integer, allocatable, dimension(:)  		:: myList_elem2D
  integer                             		:: myDim_elem3D  
  integer, allocatable, dimension(:)  		:: myList_elem3D
end module g_PARFE
!
!----------------------------------------------------------------------
!
module o_DATA_TYPES
  implicit none
  save
  !
  type sparse_matrix 
     integer :: nza
     integer :: dim
     real(kind=8), pointer, dimension(:)      :: values
     integer(KIND=4), pointer,   dimension(:) :: colind
     integer(KIND=4), pointer,   dimension(:) :: rowptr
  end type sparse_matrix
  !
  type addresstype
     integer                                :: nmb
     integer(KIND=4), dimension(:), pointer :: addresses
  end type addresstype
  !
end module o_DATA_TYPES
!
!----------------------------------------------------------------------
!
module o_MATRICES
  ! 
  use o_DATA_TYPES
  implicit none
  save
  !
  type(sparse_matrix)                             :: uvstiff, sshstiff
#ifndef use_non_hydrostatic  
  real(kind=8), allocatable, dimension(:,:,:)     :: wpot_matrix
#else
  type(sparse_matrix)                             :: nhpstiff
#endif
  type(sparse_matrix)                             :: tsstiff
  !
  real(kind=8), allocatable, dimension(:)         :: uv_lump 
  real(kind=8), allocatable, dimension(:)         :: uv_lump_prev 
  real(kind=8), allocatable, dimension(:)         :: ts_lump
end module o_MATRICES
!
!--------------------------------------------------------------------
!
module o_mesh
  !
  use o_DATA_TYPES
  implicit none
  save
  !
  integer, allocatable, dimension(:)           :: mapping
  integer, allocatable, dimension(:)           :: col_pos
  !
  integer                                      :: nod2D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod2D  
  integer, allocatable, dimension(:)           :: index_nod2D  
  integer                                      :: nod3D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod3D  
  integer(KIND=4), allocatable, dimension(:)   :: index_nod3D  
  real(kind=8), allocatable, dimension(:)      :: cos_elem2D 
  real(kind=8), allocatable, dimension(:)      :: geolat
  !
  type(addresstype), allocatable, dimension(:) :: nod_in_elem3D     
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D     
  type(addresstype), allocatable, dimension(:) :: nod_in_opbnd_tri
  type(addresstype), allocatable, dimension(:) :: nghbr_nod3D
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D
  type(addresstype), allocatable, dimension(:) :: nghbr_elem3D
  type(addresstype), allocatable, dimension(:) :: nghbr_elem2D
  !
  integer                                      :: max_num_layers
  integer, allocatable, dimension(:)           :: num_layers_below_nod2D
  integer(KIND=4), allocatable, dimension(:,:) :: nod3D_below_nod2D  
  integer(KIND=4), allocatable, dimension(:)   :: nod2D_corresp_to_nod3D 
  !
  integer(KIND=4), allocatable, dimension(:)   :: bt_nds
  !
  integer                                      :: nmbr_opbnd_n2D, nmbr_opbnd_t2D
  integer                                      :: nmbr_opbnd_n3D, nmbr_opbnd_tri
  integer                                      :: nmbr_opbnd_edg
  integer, allocatable, dimension(:)           :: opbnd_n2D, opbnd_n3D
  integer, allocatable, dimension(:)           :: mapping_opbnd_n2d
  integer, allocatable, dimension(:,:)         :: opbnd_tri, opbnd_edg
  real(kind=8), allocatable, dimension(:,:)    :: opbnd_nv, opbnd_edg_nv 

  ! for cases using cavity
  integer, allocatable, dimension(:)           :: cavity_flag_nod2d

  ! for sigma or hybrid grids
  integer, allocatable, dimension(:)           :: grid_type_elem2d
  integer, allocatable, dimension(:,:,:)       :: dens_interp_nodes
  integer, allocatable, dimension(:)           :: elem3d_layer
  real(kind=8), allocatable, dimension(:,:,:)  :: grid_slope
end module o_mesh
!
!----------------------------------------------------------------------------
!
module o_elements
  implicit none
  save
  integer                                      :: elem2D
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nodes 
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nghbrs 
  integer                                      :: elem3D
  integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nodes 
  integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nghbrs  
  integer(KIND=4), allocatable, dimension(:)   :: elem2D_corresp_to_elem3D 
  !
  real(kind=8)                                 :: Vol2D, Vol3D   
  real(kind=8)                                 :: sProd_2Di, sProd_3Di 
  real(kind=8), allocatable, dimension(:,:)    :: sProd_2Dij, sProd_3Dij 
  real(kind=8), allocatable, dimension(:,:)    :: derivative_stdbafu_x_2D 
  real(kind=8), allocatable, dimension(:,:)    :: derivative_stdbafu_x_3D  
  real(kind=8), allocatable, dimension(:,:)    :: bafux_3d, bafuy_3d, bafuz_3d
  real(kind=8), allocatable, dimension(:)      :: voltetra
  real(kind=8), allocatable, dimension(:,:)    :: bafux_2d, bafuy_2d
  real(kind=8), allocatable, dimension(:)      :: voltriangle
  !
#ifdef use_fullfreesurf
  integer(kind=4), allocatable, dimension(:)   :: map_elem
  real(kind=8), allocatable, dimension(:)      :: voltetra_new
  real(kind=8), allocatable, dimension(:,:)    :: bafux_3d_new, bafuy_3d_new, bafuz_3d_new
#endif
  !
  real(kind=8)                                 :: ocean_area
  real(kind=8),allocatable, dimension(:)       :: cluster_area_2D
  
end module o_elements
!
module i_dyn_parms
  implicit none
  save

  ! *** ice_stress ***
  real(kind=8)             :: Pstar = 30000.           ! [N/m^2] 15000. 27500. 30000.
  real(kind=8)             :: ellipse =2.              !
  real(kind=8)             :: c_pressure =20.0         !
  real(kind=8)             :: delta_min=1.0e-11        ! [s^(-1)]

  namelist /ice_stress/ Pstar, ellipse, c_pressure, delta_min

  ! *** friction setting ***
  real*8                   :: Cd_oce_ice = 5.0e-3      ! drag coefficient ocean-ice 3.e-3 5.5e-3
  real*8                   :: Kh_ice=0.0               ! numerical ice/snow diffusivity

  namelist /ice_fric/ Cd_oce_ice, Kh_ice

  ! *** rheology ***
  logical                  :: EVP_rheology=.true.
  integer                  :: evp_rheol_steps=120      ! substeps for EVP rheology
  integer                  :: evp_Tdamp_ratio=3        ! ratio dt/T_damp
  integer                  :: vp_rheol_steps=500       ! substeps for VP rheology

  namelist /ice_rheology/  EVP_rheology, evp_rheol_steps, evp_Tdamp_ratio, vp_rheol_steps

  ! *** numerical scheme ***
  real(kind=8)             :: ice_gamma_fct=0.3
  logical                  :: lump_ice_matrix=.true.   !for mass ice matrix case
  integer                  :: num_iter_solve_ice=3 

  namelist /ice_scheme/ ice_gamma_fct

  ! *** others ***
  real(kind=8)             :: Tevp_inv

end module i_dyn_parms

!=======================================================================

module i_therm_parms
  implicit none
  real*8, parameter  :: rhoair=  1.3      ! Air density ! AOMIP
  real*8, parameter  :: rhowat= 1025.     ! Water density
  real*8, parameter  :: rhoice=  910.     ! Ice density
  real*8, parameter  :: rhosno=  290.     ! Snow density

  real*8, parameter  :: cpair=1005.       ! Specific heat of air [J/(kg * K)] / 1004 
  real*8, parameter  :: cc=rhowat*4190.0  ! Volumetr. heat cap. of water [J/m**3/K](cc = rhowat*cp_water)
  real*8, parameter  :: cl=rhoice*3.34e5  ! Volumetr. latent heat of ice fusion [J/m**3](cl=rhoice*Lf) 
  real*8, parameter  :: clhw=2.501e6      ! Specific latent heat [J/kg]: water	-> water vapor
  real*8, parameter  :: clhi=2.835e6      !                              sea ice-> water vapor

  real*8, parameter  :: tmelt=273.15      ! 0 deg C expressed in K 
  real*8, parameter  :: boltzmann=5.67E-8 ! S. Boltzmann const.*longw. emissivity

  real*8, parameter  :: con   = 2.1656    ! Thermal conductivities: ice; W/m/K
  real*8, parameter  :: consn = 0.31	  !                         snow

  real*8, parameter  :: hmin= 0.05        ! Cut-off ice thickness     !!
  real*8, parameter  :: Armin=0.05        ! Minimum ice concentration !!

  integer,parameter  :: iclasses=7        ! Number of ice thickness gradations for ice growth calcs.

  real*8    :: Sice = 5.0        ! Ice salinity 3.2--5.0 ppt.
  real*8    :: h0=1.0	         ! Lead closing parameter [m] ! 0.5

  real*8    :: emiss_ice=0.98    ! Emissivity of Snow/Ice
  real*8    :: emiss_wat=0.97    ! of open water  

  real*8    :: albsn=	0.81 	 ! Albedo: frozen snow
  real*8    :: albsnm=  0.77   	 !         melting snow
  real*8    :: albi=	0.70 	 !         frozen ice
  real*8    :: albim=	0.68 	 !         melting ice
  real*8    :: albw=	0.1	 !         open water

  namelist /ice_therm/ Sice, h0, emiss_ice, emiss_wat, albsn, albsnm, albi, albim, albw

end module i_therm_parms

!=======================================================================

module i_array
  use o_data_types
  implicit none
  save
  real(kind=8), allocatable, dimension(:)         :: u_ice, v_ice
  real(kind=8), allocatable, dimension(:)         :: m_ice, a_ice, m_snow  
  real(kind=8), allocatable, dimension(:)         :: rhs_u, rhs_v
  real(kind=8), allocatable, dimension(:)         :: rhs_m, rhs_a, rhs_ms
  real(kind=8), allocatable, dimension(:)         :: dm_ice, da_ice, dm_snow
 
  real(kind=8), allocatable, dimension(:)         :: t_skin

  type(sparse_matrix)                             :: icestiff
  real(kind=8), allocatable, dimension(:)         :: ice_lump

  real(kind=8), allocatable, dimension(:)         :: m_icel, a_icel, m_snowl
  real(kind=8), allocatable, dimension(:)         :: icepplus, icepminus 
  real(kind=8), allocatable, dimension(:,:)       :: icefluxes

  real(kind=8), allocatable, dimension(:)         :: sigma11, sigma12, sigma22  

  real(kind=8), allocatable, dimension(:)         :: u_w, v_w
  real(kind=8), allocatable, dimension(:)         :: fresh_wa_flux
  real(kind=8), allocatable, dimension(:)         :: net_heat_flux
  real(kind=8), allocatable, dimension(:)         :: S_oc_array, T_oc_array
  real(kind=8), allocatable, dimension(:)         :: elevation

  real(kind=8), allocatable, dimension(:)         :: stress_iceoce_x         
  real(kind=8), allocatable, dimension(:)         :: stress_iceoce_y
  real(kind=8), allocatable, dimension(:)         :: stress_atmice_x         
  real(kind=8), allocatable, dimension(:)         :: stress_atmice_y
  real(kind=8), allocatable, dimension(:)         :: stress_atmoce_x         
  real(kind=8), allocatable, dimension(:)         :: stress_atmoce_y

end module i_array

!=====================================================================

module i_solver
implicit none
save
  integer                  	:: solve_m_ice=101
  integer       	 	:: solve_a_ice=102
  integer			:: solve_m_snow=103
end module i_solver
module g_config
  implicit none
  save

  ! *** Modelname ***
  character(5)             	:: runid='test1'                ! a model/setup name

  namelist /modelname/ runid

  ! *** time step ***
  integer                  	:: step_per_day=12           	!number of steps per day
  integer                  	:: run_length=1	                !run length
  character                     :: run_length_unit='y'          !unit: y, d, s

  namelist /timestep/ step_per_day, run_length, run_length_unit

  ! *** Paths for all in and out ***
  character(100)                :: MeshPath='./mesh/'
  character(100)                :: OpbndPath='./opbnd/'
  character(100)                :: ClimateDataPath='./hydrography/'
  character(100)                :: ForcingDataPath='./forcing/'
  character(100)                :: TideForcingPath='./tide_forcing/'
  character(100)                :: ResultPath='./result/'

  namelist /paths/  MeshPath, OpbndPath, ClimateDataPath, ForcingDataPath, &
       TideForcingPath, ResultPath

  ! *** ocean climatology data name ***
  character(100)                :: OceClimaDataName='annual_woa01_ts.out'
  logical                       :: use_prepared_init_ice=.false.     !how to initial. ice at the beginning 

  namelist /initialization/ OceClimaDataName, use_prepared_init_ice

  ! *** in out ***
  character*4               	:: restartflag='last'  	             !restart from which saved record,'#','last'
  integer                       :: output_length=1                   !valid for d,h,s
  character                	:: output_length_unit='m'      	     !output period: y, m, d, h, s 
  integer                       :: logfile_outfreq=1                 !in logfile info. output frequency, # steps

  namelist /inout/ restartflag, output_length, output_length_unit, logfile_outfreq

  ! *** mesh ***
  integer                       :: grid_type=1              	! z-level, 2 sigma, 3 sigma + z-level

  namelist /mesh_def/ grid_type

  ! *** model geometry
  logical                  	:: cartesian=.false.
  logical                  	:: fplane=.false.
  logical                  	:: betaplane=.false.
  real(kind=8)             	:: f_fplane=-1.4e-4        	![1/s]
  real(kind=8)             	:: beta_betaplane=2.0e-11  	![1/s/m]
  real(kind=8)             	:: domain_length=360.    	![degree]
  !
  logical                  	:: rotated_grid=.true.    	!option only valid for coupled model case now
  real(kind=8)             	:: alphaEuler=-30. 		![degree] Euler angles, convention:
  real(kind=8)             	:: betaEuler=-90.  		![degree] first around z, then around new x,
  real(kind=8)			:: gammaEuler=-90.		![degree] then around new z.

  namelist /geometry/  cartesian, fplane, betaplane, f_fplane, beta_betaplane, &
       domain_length, rotated_grid, alphaEuler, betaEuler, gammaEuler

  ! *** fleap_year ***
  logical                       :: include_fleapyear=.false.
  
  namelist /calendar/ include_fleapyear
  
   ! *** machine ***
  integer                       :: system=1                     ! XD1 2(byte), HLRN 1(word)
  
  namelist /machine/ system


  ! *** others ***
  real(kind=8)             	:: dt, dt_inv
  integer                  	:: istep, nsteps
  integer                       :: save_count
  logical                       :: r_restart

end module g_config
module g_clock
  !combining RT and Lars version
  !
  use g_config
  implicit none
  save
  real(kind=8)             :: timeold, timenew     !time in a day, unit: sec
  integer                  :: dayold, daynew       !day in a year
  integer                  :: yearold, yearnew     !year before and after time step
  integer                  :: month, day_in_month  !month and day in a month
  integer                  :: fleapyear            !1 fleapyear, 0 not 
  integer                  :: ndpyr                !number of days in yearnew 
  integer                  :: num_day_in_month(0:1,12)
  character(4)             :: cyearold, cyearnew   !year as character string      

  data num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/


contains
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock

    implicit none
    integer         :: i
    real(kind=8)    :: aux1, aux2
    !
    timeold=timenew 
    dayold=daynew
    yearold=yearnew

    ! update time
    timenew=timenew+dt          

    ! update day
    if (timenew>86400.) then  !assumed that time step is less than one day!
       daynew=daynew+1
       timenew=timenew-86400.
    endif

    ! update year
    if (daynew>ndpyr) then
       daynew=1
       yearnew=yearnew+1
       call check_fleapyr(yearnew, fleapyear)
       ndpyr=365+fleapyear
       write(cyearold,'(i4)') yearold
       write(cyearnew,'(i4)') yearnew
    endif

    ! find month and dayinmonth at new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

  end subroutine clock
  !
  !--------------------------------------------------------------------------------
  !
  subroutine clock_init
    use o_param
    use g_parfe
    implicit none
    integer         :: i, daystart, yearstart
    real(kind=8)    :: aux1, aux2, timestart

    ! the model inialized at
    timestart=timenew
    daystart=daynew
    yearstart=yearnew

    ! init clock for this run
    open(99,file=trim(ResultPath)//runid//'.clock',status='old')
    read(99,*) timeold, dayold, yearold
    read(99,*) timenew, daynew, yearnew
    close(99)
    if(daynew==0) daynew=1

    ! check if this is a restart or not
    if(yearnew==yearstart .and. daynew==daystart .and. timenew==timestart) then
       r_restart=.false.
       yearold=yearnew-1 !required for checking if create new output files
    else
       r_restart=.true.
    end if

    ! year as character string 
    write(cyearold,'(i4)') yearold
    write(cyearnew,'(i4)') yearnew

    ! if restart model at beginning of a day, set timenew to be zero
    if (timenew==86400.) then  
       timenew=0.0
       daynew=daynew+1
    endif

    ! set timeold to be timenew, ready for initializing forcing fields,
    ! yearold should not be updated here, which is requird to open input files.
    timeold=timenew 
    dayold=daynew

    ! check fleap year
    call check_fleapyr(yearnew, fleapyear)
    ndpyr=365+fleapyear

    ! find month and dayinmonth at the new time step
    aux1=0
    do i=1,12
       aux2=aux1+num_day_in_month(fleapyear,i)
       if(daynew>aux1 .and. daynew<=aux2) then
          month=i
          day_in_month=daynew-aux1
          exit
       end if
       aux1=aux2
    end do

    if(mype==0) then
       write(*,*)'clock initialized at time ', real(timenew,4), daynew, yearnew
       if(r_restart) then
          write(*,*) 'THIS RUN IS A RESTART RUN!'
       end if
    end if

  end subroutine clock_init
  !
  !-------------------------------------------------------------------------------
  !
  subroutine clock_finish
    use o_param
    implicit none
    !
    if ((daynew==ndpyr) .and. (timenew==86400.)) then
       timenew=0.0
       daynew=1
       yearnew=yearold+1
    endif

    open(99,file=trim(ResultPath)//runid//'.clock',status='unknown')
    write(99,*) timeold, dayold, yearold
    write(99,*) timenew, daynew, yearnew
    close(99)
  end subroutine clock_finish
  !
  !----------------------------------------------------------------------------
  !
  subroutine clock_newyear
    implicit none
    !
    if ((daynew>=ndpyr).and.(timenew==86400.)) then
       timenew=0.0
       daynew=1
       yearnew=yearold+1
       write(cyearnew,'(i4)') yearnew
    endif
  end subroutine clock_newyear
  !
  !----------------------------------------------------------------------------
  !
  subroutine check_fleapyr(year, flag)
    use O_param
    implicit none
    integer, intent(in) :: year      
    integer, intent(out):: flag

    flag=0

    if(.not.include_fleapyear) return

    if ((mod(year,4)==0.and.mod(year,100)/=0) .or. mod(year,400)==0) then
       flag=1
    endif
  end subroutine check_fleapyr
  !
  !----------------------------------------------------------------------------
  !
end module g_clock
module g_rotate_grid
  use g_config
  use g_parfe

  implicit none
  save
  real(kind=8)        :: rotate_matrix(3,3)

contains


  !----------------------------------------------------------------
  !
  subroutine calculate_rotate_matrix
    ! A, B, G [radian] are Euler angles.
    ! The convention (A, B, G) used here is: the first rotation is by an
    ! angle A around z-axis, the second is by an angle B about the new 
    ! x-axis, and the third is by an angle G about the new z-axis.   
    implicit none

    real(kind=8)      :: al, be, ga

    al=alphaEuler
    be=betaEuler
    ga=gammaEuler

    ! rotation matrix
    rotate_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
    rotate_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
    rotate_matrix(1,3)=sin(ga)*sin(be)
    rotate_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
    rotate_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
    rotate_matrix(2,3)=cos(ga)*sin(be)
    rotate_matrix(3,1)=sin(be)*sin(al) 
    rotate_matrix(3,2)=-sin(be)*cos(al)  
    rotate_matrix(3,3)=cos(be)

    if(mype==0) write(*,*) 'rotation matrix for rotated model grids prepared'
  end subroutine calculate_rotate_matrix
  !
  !----------------------------------------------------------------
  !
  subroutine r2g(lon, lat, rlon, rlat)
    ! Convert the rotated coordinates to geographical coordinates  
    ! lon, lat		:: [radian] geographical coordinates
    ! rlon, rlat	:: [radian] rotated coordinates
    !
    implicit none
    real(kind=8), intent(out)      :: lon, lat
    real(kind=8), intent(in)       :: rlon, rlat
    real(kind=8)                   :: xr, yr, zr, xg, yg, zg
    !
    ! Rotated Cartesian coordinates:
    xr=cos(rlat)*cos(rlon)
    yr=cos(rlat)*sin(rlon)
    zr=sin(rlat)

    ! Geographical Cartesian coordinates:
    xg=rotate_matrix(1,1)*xr + rotate_matrix(2,1)*yr + rotate_matrix(3,1)*zr
    yg=rotate_matrix(1,2)*xr + rotate_matrix(2,2)*yr + rotate_matrix(3,2)*zr  
    zg=rotate_matrix(1,3)*xr + rotate_matrix(2,3)*yr + rotate_matrix(3,3)*zr  

    ! Geographical coordinates:
    lat=asin(zg)
    if(yg==0. .and. xg==0.) then
       lon=0.0     ! exactly at the poles
    else
       lon=atan2(yg,xg)
    end if
  end subroutine r2g
  !
  !----------------------------------------------------------------
  !
  subroutine g2r(lon, lat, rlon, rlat)
    ! Convert the geographical coordinates to rotated coordinates  
    ! lon, lat		:: [radian] geographical coordinates
    ! rlon, rlat	:: [radian] rotated coordinates
    !
    implicit none
    real(kind=8), intent(in)       :: lon, lat
    real(kind=8), intent(out)      :: rlon, rlat
    real(kind=8)                   :: xr, yr, zr, xg, yg, zg
    !
    ! geographical Cartesian coordinates:
    xg=cos(lat)*cos(lon)
    yg=cos(lat)*sin(lon)
    zg=sin(lat)

    ! rotated Cartesian coordinates:
    xr=rotate_matrix(1,1)*xg + rotate_matrix(1,2)*yg + rotate_matrix(1,3)*zg
    yr=rotate_matrix(2,1)*xg + rotate_matrix(2,2)*yg + rotate_matrix(2,3)*zg  
    zr=rotate_matrix(3,1)*xg + rotate_matrix(3,2)*yg + rotate_matrix(3,3)*zg  

    ! rotated coordinates:
    rlat=asin(zr)
    if(yr==0. .and. xr==0.) then
       rlon=0.0     ! exactly at the poles
    else
       rlon=atan2(yr,xr)
    end if
  end subroutine g2r
  !
  !--------------------------------------------------------------------
  !
  subroutine vector_g2r(tlon, tlat, lon, lat, flag_coord)
    ! rotate 2d vector (tlon, tlat) to be in the rotated coordinates
    ! tlon, tlat (in)	:: lon & lat components of a vector in geographical coordinates 
    !            (out)	:: lon & lat components of the vector in rotated coordinates              
    ! lon, lat	        :: [radian] coordinates
    ! flag_coord        :: 1, (lon,lat) is the geographical coord.; else, rotated coord.
    !
    implicit none
    integer, intent(in)           :: flag_coord
    real(kind=8), intent(inout)   :: tlon, tlat
    real(kind=8), intent(in)      :: lon, lat
    real(kind=8)                  :: rlon, rlat, glon, glat
    real(kind=8)		  :: txg, tyg, tzg, txr, tyr, tzr
    !
    ! geographical coordinate
    if(flag_coord==1) then  ! input is geographical coordinates
       glon=lon
       glat=lat
       call g2r(glon,glat,rlon,rlat)
    else                    ! input is rotated coordinates 
       rlon=lon
       rlat=lat
       call r2g(glon,glat,rlon,rlat)
    end if
    !
    ! vector in Cartesian
    txg=-tlat*sin(glat)*cos(glon)-tlon*sin(glon)
    tyg=-tlat*sin(glat)*sin(glon)+tlon*cos(glon)
    tzg=tlat*cos(glat)
    !
    ! vector in rotated Cartesian
    txr=rotate_matrix(1,1)*txg + rotate_matrix(1,2)*tyg + rotate_matrix(1,3)*tzg 
    tyr=rotate_matrix(2,1)*txg + rotate_matrix(2,2)*tyg + rotate_matrix(2,3)*tzg 
    tzr=rotate_matrix(3,1)*txg + rotate_matrix(3,2)*tyg + rotate_matrix(3,3)*tzg 
    !
    ! vector in rotated coordinate
    tlat=-sin(rlat)*cos(rlon)*txr - sin(rlat)*sin(rlon)*tyr + cos(rlat)*tzr
    tlon=-sin(rlon)*txr + cos(rlon)*tyr

  end subroutine vector_g2r
  !
end module g_rotate_grid
module g_read_CORE_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid

contains 

  subroutine read_CORE_NetCDF(file,vari,itime,ncdata)
    use g_config
    use g_parfe
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer, dimension(3)           	:: istart, icount
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    real(kind=8), dimension (nci,ncj)   :: ncdata
    character(15)                       :: vari
    character(80)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    io=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    io=nf_close(ncid)

    return
  end subroutine read_CORE_NetCDF

end module g_read_CORE_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_NCEP_NetCDF
  ! read NetCDF data on T62 NCEP/NCAR grid

contains 

  subroutine read_NCEP_NetCDF(file,vari,itime,ncdata)
    use g_config
    use g_parfe
    implicit none

#include "netcdf.inc" 

    integer, parameter             	:: nci=192, ncj=94  ! T62 grid dimension
    integer				:: ncid, varid, io
    integer                             :: i, j, itime
    integer, dimension(3)           	:: istart, icount
    integer(kind=2), dimension(nci,ncj) :: iuw 
    real(kind=4)                        :: xscale, xoff, miss
    real(kind=8), dimension(nci,ncj)    :: ncdata
    character(15)                       :: vari
    character(80)                      	:: file

    istart = (/1,1,itime/)
    icount= (/nci,ncj,1/)

    ! open netcdf input file
    io=nf_open(file, nf_nowrite, ncid)

    if (io.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ FORCING FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop 'Fatal error in open_netcdf'
    endif

    ! get handle for data variable    
    io=nf_inq_varid(ncid, vari, varid)

    ! get att
    io=nf_get_att_real(ncid,varid,'scale_factor',xscale)
    io= nf_get_att_real(ncid,varid,'add_offset',xoff)
    ! io= nf_get_att_int(ncid,varid,'missing_value',miss)

    ! get variable
    io=nf_get_vara_int2(ncid,varid,istart,icount,iuw)

    ! close file 
    io=nf_close(ncid)

    ! ncdata
    do j=1,ncj
       do i=1,nci
          ncdata(i,j)=real(iuw(i,j))*xscale+xoff
       enddo
    end do

    return
  end subroutine read_NCEP_NetCDF
  !
  !--------------------------------------------------------------
  !
  subroutine upside_down(array,ni,nj)
    implicit none

    integer                         :: i, j, ji
    integer                         :: ni, nj
    real(kind=8), dimension(ni,nj)  :: array
    real(kind=8), dimension(nj)     :: ywrk

    do i=1,ni  
       do j=1,nj
          ywrk(j)=array(i,j)
       end do
       do j=1,nj
          ji =nj-j+1
          array(i,ji)=ywrk(j)            
       enddo
    enddo

    return
  end subroutine upside_down

end module g_read_NCEP_NetCDF
!
!------------------------------------------------------------------------------------
!
module g_read_other_NetCDF
  ! Read global data with NetCDF format and interpolate to the model grid.
  ! Having different versions for different cases 

contains

  subroutine read_other_NetCDF(file, vari, itime, model_2Darray, check_dummy)
    ! Read 2D data and interpolate to the model grid.
    ! Currently used for reading runoff and SSS.
    ! First missing value will be replaced on the raw regular grid;
    ! Then calling interp_2d_field to do interpolation.

    ! The check_dummy part should be carefully modified in new applications!
    ! if check_dummy=.true., replace missing value with a meaningful value nearby
    ! if check_dummy=.false., replace missing value with 0.0

    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parfe
    implicit none

#include "netcdf.inc" 

    integer			:: i, j, ii, jj, k, n, num, flag, cnt
    integer			:: itime, latlen, lonlen
    integer			:: status, ncid, varid
    integer			:: lonid, latid
    integer			:: istart(3), icount(3)
    real(kind=8)		:: x, y, miss, aux
    real(kind=8), allocatable	:: lon(:), lat(:)
    real(kind=8), allocatable	:: ncdata(:,:), ncdata_temp(:,:)
    real(kind=8), allocatable	:: temp_x(:), temp_y(:)
    real(kind=8)		:: model_2Darray(myDim_nod2d+eDim_nod2D)   
    character(15)		:: vari
    character(80)              	:: file
    logical                     :: check_dummy

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    allocate(lat(latlen))
    status=nf_inq_varid(ncid, 'lat', varid)
    status=nf_get_vara_double(ncid,varid,1,latlen,lat)

    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
    allocate(lon(lonlen))
    status=nf_inq_varid(ncid, 'lon', varid)
    status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
    ! make sure range 0. - 360.
    do n=1,lonlen
       if(lon(n)<0.0) then
          lon(n)=lon(n)+360.
       end if
    end do

    ! data
    allocate(ncdata(lonlen,latlen), ncdata_temp(lonlen,latlen))
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1,1,itime/)
    icount= (/lonlen,latlen,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
    status= nf_get_att_double(ncid,varid,'missing_value',miss)
    !write(*,*)'miss', miss
    !write(*,*)'raw',minval(ncdata),maxval(ncdata)
    ncdata_temp=ncdata
    do i=1,lonlen
       do j=1,latlen
          if(ncdata(i,j)==miss .or. ncdata(i,j)==-99.0_8) then  !!
             if(check_dummy) then
                aux=0.0
                cnt=0
                do k=1,30
                   do ii=max(1,i-k),min(lonlen,i+k)
                      do jj=max(1,j-k),min(latlen,j+k)
                         if(ncdata_temp(ii,jj)/=miss .and. ncdata_temp(ii,jj)/=-99.0_8) then  !!
                            aux=aux+ncdata_temp(ii,jj)
                            cnt=cnt+1                         
                         end if
                      end do	!ii
                   end do	!jj
                   if(cnt>0) then
                      ncdata(i,j)=aux/cnt
                      exit
                   end if
                end do  	!k    
             else
                ncdata(i,j)=0.0
             end if
          end if
       end do
    end do
    !write(*,*) 'post',minval(ncdata), maxval(ncdata)

    ! close file
    status=nf_close(ncid)
    ! model grid coordinates
    num=myDim_nod2d+eDim_nod2d
    allocate(temp_x(num), temp_y(num))  
    do n=1, num                        
       if(rotated_grid) then
          call r2g(x, y, coord_nod2d(1,n), coord_nod2d(2,n))
          temp_x(n)=x/rad   ! change unit to degree  
          temp_y(n)=y/rad                             
       else
          temp_x(n)=coord_nod2d(1,n)/rad              
          temp_y(n)=coord_nod2d(2,n)/rad             
       end if
       ! change lon range to [0 360]
       if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
    end do

    ! interpolation
    flag=0
    call interp_2d_field(lonlen, latlen, lon, lat, ncdata, num, temp_x, temp_y, & 
         model_2Darray, flag) 
    deallocate(temp_y, temp_x, ncdata_temp, ncdata, lon, lat)

  end subroutine read_other_NetCDF
  !
  !------------------------------------------------------------------------------------
  !
  subroutine read_surf_hydrography_NetCDF(file, vari, itime, model_2Darray)
    ! Read WOA (NetCDF) surface T/S and interpolate to the model grid.
    ! Currently used for surface restoring in case of ocean-alone models
    ! Calling interp_2d_field_v2 to do interpolation, which also treats the dummy value.

    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parfe
    implicit none

#include "netcdf.inc" 

    integer			:: i, j,  n, num
    integer			:: itime, latlen, lonlen
    integer			:: status, ncid, varid
    integer			:: lonid, latid
    integer			:: istart(4), icount(4)
    real(kind=8)		:: x, y, miss
    real(kind=8), allocatable	:: lon(:), lat(:)
    real(kind=8), allocatable	:: ncdata(:,:)
    real(kind=8), allocatable	:: temp_x(:), temp_y(:)
    real(kind=8)		:: model_2Darray(myDim_nod2d+eDim_nod2D)   
    character(15)		:: vari
    character(80)              	:: file
    logical                     :: check_dummy

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! lat
    status=nf_inq_dimid(ncid, 'lat', latid)
    status=nf_inq_dimlen(ncid, latid, latlen)
    allocate(lat(latlen))
    status=nf_inq_varid(ncid, 'lat', varid)
    status=nf_get_vara_double(ncid,varid,1,latlen,lat)

    ! lon
    status=nf_inq_dimid(ncid, 'lon', lonid)
    status=nf_inq_dimlen(ncid, lonid, lonlen)
    allocate(lon(lonlen))
    status=nf_inq_varid(ncid, 'lon', varid)
    status=nf_get_vara_double(ncid,varid,1,lonlen,lon)
    ! make sure range 0. - 360.
    do n=1,lonlen
       if(lon(n)<0.0) then
          lon(n)=lon(n)+360.
       end if
    end do

    ! data
    allocate(ncdata(lonlen,latlen))
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1,1,1,itime/)
    icount= (/lonlen,latlen,1,1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)

    ! missing value
    status= nf_get_att_double(ncid,varid,'missing_value',miss)
    !write(*,*)'miss', miss
    !write(*,*)'raw',minval(ncdata),maxval(ncdata)

    ! close file
    status=nf_close(ncid)

    ! the next step is to interpolate data to model grids

    ! model grid coordinates
    num=myDim_nod2d+eDim_nod2d
    allocate(temp_x(num), temp_y(num))  
    do n=1, num                        
       if(rotated_grid) then
          call r2g(x, y, coord_nod2d(1,n), coord_nod2d(2,n))
          temp_x(n)=x/rad   ! change unit to degree  
          temp_y(n)=y/rad                             
       else
          temp_x(n)=coord_nod2d(1,n)/rad              
          temp_y(n)=coord_nod2d(2,n)/rad             
       end if
       ! change lon range to [0 360]
       if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
    end do

    ! interpolation
    call interp_2d_field_v2(lonlen, latlen, lon, lat, ncdata, miss, &
         num, temp_x, temp_y, model_2Darray) 

    deallocate(temp_y, temp_x, ncdata, lon, lat)

  end subroutine read_surf_hydrography_NetCDF
  !
  !------------------------------------------------------------------------------------
  !
  subroutine read_2ddata_on_grid_NetCDF(file, vari, itime, model_2Darray)  
    ! read 2D data which are already prepared on the model grid
    use g_config
    use o_param
    use o_mesh
    use g_rotate_grid
    use g_parfe
    implicit none

#include "netcdf.inc" 

    integer			:: n, i
    integer			:: itime
    integer			:: status, ncid, varid
    integer			:: istart(2), icount(2)
    real(kind=8)           	:: ncdata(nod2D)
    real(kind=8), intent(out)	:: model_2Darray(ToDim_nod2d)
    character(80), intent(in) 	:: file
    character(15), intent(in)  :: vari

    ! open file
    status=nf_open(file, nf_nowrite, ncid)
    if (status.ne.nf_noerr)then
       print*,'ERROR: CANNOT READ runoff FILE CORRECTLY !!!!!'
       print*,'Error in opening netcdf file'//file
       stop
    endif

    ! get variables
    status=nf_inq_varid(ncid, vari, varid)
    istart = (/1, itime/)
    icount= (/nod2D, 1/)
    status=nf_get_vara_double(ncid,varid,istart,icount,ncdata)
    status=nf_close(ncid)
      
    model_2Darray=ncdata(myList_nod2D)     

  end subroutine read_2ddata_on_grid_NetCDF
  
end module g_read_other_NetCDF


module g_forcing_param
  implicit none
  save   

  ! *** exchange coefficients ***
  real*8    :: Ce_atm_oce=1.75e-3 ! exchange coeff. of latent heat over open water
  real*8    :: Ch_atm_oce=1.75e-3 ! exchange coeff. of sensible heat over open water
  real*8    :: Cd_atm_oce=1.0e-3  ! drag coefficient between atmosphere and water

  real*8    :: Ce_atm_ice=1.75e-3 ! exchange coeff. of latent heat over ice
  real*8    :: Ch_atm_ice=1.75e-3 ! exchange coeff. of sensible heat over ice
  real*8    :: Cd_atm_ice=1.32e-3 ! drag coefficient between atmosphere and ice 

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice


  ! *** forcing source and type ***
  character(10)                 :: wind_data_source='CORE2'
  character(10)                 :: rad_data_source='CORE2'
  character(10)                 :: precip_data_source='CORE2'
  character(10)                 :: runoff_data_source='CORE2'
  character(10)                 :: sss_data_source='CORE2'
  integer                       :: wind_ttp_ind=1
  integer                       :: rad_ttp_ind=2
  integer                       :: precip_ttp_ind=3
  integer                       :: runoff_ttp_ind=0
  integer                       :: sss_ttp_ind=4


  namelist /forcing_source/ wind_data_source, rad_data_source, precip_data_source, &
       runoff_data_source, sss_data_source, wind_ttp_ind, rad_ttp_ind, precip_ttp_ind, &
       runoff_ttp_ind, sss_ttp_ind

  ! *** coefficients in bulk formulae ***
  logical                       :: AOMIP_drag_coeff=.false.
  logical                       :: ncar_bulk_formulae=.false.

  namelist /forcing_bulk/ AOMIP_drag_coeff, ncar_bulk_formulae

  ! *** add land ice melt water ***
  logical                       :: use_landice_water=.false.
  integer                       :: landice_start_mon=1
  integer                       :: landice_end_mon=12

  namelist /land_ice/ use_landice_water, landice_start_mon, landice_end_mon

end module g_forcing_param
!
!----------------------------------------------------------------------------
!
module g_forcing_arrays
  implicit none
  save    

  ! forcing arrays
  real(kind=8), allocatable, dimension(:)         :: u_wind, v_wind 
  real(kind=8), allocatable, dimension(:)         :: Tair, shum
  real(kind=8), allocatable, dimension(:)         :: shortwave, longwave
  real(kind=8), allocatable, dimension(:)         :: prec_rain, prec_snow
  real(kind=8), allocatable, dimension(:)         :: runoff, evaporation
  
  real(kind=8), allocatable, dimension(:)         :: runoff_landice
  real(kind=8)                                    :: landice_season(12)

  ! shortwave penetration
  real(kind=8), allocatable, dimension(:)         :: chl, sw_3d

  real(kind=8), allocatable, dimension(:)         :: thdgr, thdgrsn, flice
  real(kind=8), allocatable, dimension(:)         :: olat_heat, osen_heat, olwout

  ! drag coefficient Cd_atm_oce and transfer coefficients for evaporation
  ! Ce_atm_oce and sensible heat Ch_atm_oce between atmosphere and ocean
  real(kind=8), allocatable, dimension(:)	  :: Cd_atm_oce_arr
  real(kind=8), allocatable, dimension(:)	  :: Ch_atm_oce_arr
  real(kind=8), allocatable, dimension(:)	  :: Ce_atm_oce_arr

  ! drag coefficient Cd_atm_oce between atmosphere and ice
  real(kind=8), allocatable, dimension(:)	  :: Cd_atm_ice_arr

end module g_forcing_arrays
!
!----------------------------------------------------------------------------
!
module g_forcing_index
  ! I uncommonted the part required for interpolation in time
  ! Interpolation is not supposed to be used
  ! 12,Oct. 2010
  
  implicit none
  save

  ! arrays for temporal interpolation
  integer                                         :: update_forcing_flag(4)
  integer                                         :: forcing_rec(4)
!!$  real(kind=8)                                    :: interp_coef(4)

contains

  subroutine forcing_index
    use g_clock
    implicit none

    real(kind=8)          :: sixhour_sec
    real(kind=8)          :: oneday_sec
    real(kind=8)          :: modtimeold

    data sixhour_sec /21600.0/, oneday_sec /86400.0/

    modtimeold=mod(timeold,oneday_sec)

    ! if update forcing or not
    update_forcing_flag=0
    if(mod(timeold, sixhour_sec)==0.0) update_forcing_flag(1)=1
    if(modtimeold==0.0) update_forcing_flag(2)=1
    if(day_in_month==1 .and. modtimeold==0.0) update_forcing_flag(3:4)=1

    ! which record should be used as the first one in interpolation
    forcing_rec(1) = 1+int(modtimeold/sixhour_sec)+4*(daynew-1)
    forcing_rec(2) = daynew
    forcing_rec(3) = month
    forcing_rec(4) = month

!!$    ! interpolation coefficients
!!$    interp_coef(1)=mod(timeold, sixhour_sec)/sixhour_sec
!!$    interp_coef(2)=modtimeold/oneday_sec
!!$    interp_coef(3)=(day_in_month-1.0+modtimeold/oneday_sec) &
!!$         /real(num_day_in_month(fleapyear,month))
!!$    interp_coef(4)=interp_coef(3)
!!$
!!$    if(any(interp_coef>1.) .or. any(interp_coef<0.)) then
!!$       write(*,*) 'error in interp_coef'
!!$       stop
!!$    end if

  end subroutine forcing_index

end module g_forcing_index
!
!----------------------------------------------------------------------------
!
module g_forcing_interp
  !This module prepare the weights for interpolating 
  !forcing data Tair, humidity, wind velocities,
  !precipitation, radiations, etc.
  !Based on assumption forcing data are on the T62 NCEP/NCAR grid

  integer, allocatable      :: lint_ind(:,:,:)
  real(kind=8), allocatable :: lint_weight(:,:)

contains

  !------------------------------------------------------------------
  subroutine  init_forcing_interp

    ! routine to calculate neighbours and weights for linear interpolation
    !
    ! required information
    !    xmod(nmp)  longitudes of model point on geographical grid in degree (no order assumed)
    !    ymod(nmp)  latitudes of model points on geographical grid in degree (no order assumed)
    !         nmp   number of model points where data are needed
    !    cx(ni)     longitudes of data points on regular geographical grid
    !               by now must be in range[0,..360] in ascending order
    !    cy(nj)     latitudes of data points on regular geographical grid 
    !               by now must be in range [-90,..90] in ascending order
    !
    ! OUTPUT
    !    lint_ind(4,2,nmp)   /i,j/ adress for four neighbors of each model node
    !    lint_weight(4,nmp)  interpolation weights

    use o_mesh
    use o_param
    use g_config
    use g_rotate_grid
    use g_parfe
    implicit none

    integer, parameter 	:: ni=192, nj=94  ! NCEP and CORE are on the same grid.
    integer     	:: i, ii, j, n, row, n2
    real(kind=8)      	:: rlon, rlat, aux
    real(kind=8)      	:: cx(ni), cy(nj)
    real(kind=8)      	:: xmod(myDim_nod2D+eDim_nod2D), ymod(myDim_nod2D+eDim_nod2D)
    real(kind=8)        :: wt(4)

    n2=myDim_nod2D+eDim_nod2D     

    allocate(lint_ind(4,2,n2))   
    allocate(lint_weight(4,n2))  

    ! NCEP/CORE latitude
    data cy /-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
         -77.1394, -75.2351, -73.3307, -71.4262, -69.5217, -67.6171, &  
         -65.7125, -63.8079, -61.9033, -59.9986, -58.0939, -56.1893, &  
         -54.2846, -52.3799, -50.4752, -48.5705, -46.6658, -44.7611,&  
         -42.8564, -40.9517, -39.0470, -37.1422, -35.2375, -33.3328, &  
         -31.4281, -29.5234, -27.6186, -25.7139, -23.8092, -21.9044, &  
         -19.9997, -18.0950, -16.1902, -14.2855, -12.3808, -10.47604, &  
         -8.57131, -6.66657, -4.76184, -2.8571, -0.952368, 0.952368, &  
         2.8571, 4.76184, 6.66657, 8.57131, 10.47604, 12.3808, &  
         14.2855, 16.1902, 18.095, 19.9997, 21.9044, 23.8092, &  
         25.7139, 27.6186, 29.5234, 31.4281, 33.3328, 35.2375,&  
         37.1422, 39.047,  40.9517, 42.8564, 44.7611, 46.6658,&  
         48.5705, 50.4752, 52.3799, 54.2846, 56.1893, 58.0939,&  
         59.9986, 61.9033, 63.8079, 65.7125, 67.6171, 69.5217, &  
         71.4262, 73.3307, 75.2351, 77.1394, 79.0435, 80.9473, &  
         82.8508, 84.7532, 86.6531, 88.542 /

    ! NCEP/CORE longitude
    cx(1)=0.0
    do i=2,ni
       cx(i)=cx(i-1)+1.875
    enddo

    ! some checks, range of cx and cy
    if(cx(ni)-cx(1).gt.360.) then
       write(*,*) 'gauss_init: x-interval gt 360'
       call abort
    endif
    if(cy(nj)-cy(1).gt.180.) then
       write(*,*) 'gauss_init: y-interval gt 180'
       call abort
    endif
    if(cx(1).ge.360.0) then
       aux=int(cx(1)/360.)*360.
       do i=1,ni
          cx(i)=cx(i)-aux
       enddo
    endif
    if(cx(ni).gt.360.) call abort

    ! in the following we need cx and cy in unit radian
    cx=cx*rad
    cy=cy*rad

    ! model grid coordinate (in radian, between 0 and 2*pi)
    do row=1,n2                    
       if(rotated_grid) then
          rlon=coord_nod2D(1,row)
          rlat=coord_nod2D(2,row)
          call r2g(xmod(row), ymod(row), rlon, rlat)
          if (xmod(row)<0.0) xmod(row)=2*pi+xmod(row)
       else
          xmod(row)=coord_nod2D(1,row)
          ymod(row)=coord_nod2D(2,row)
          if (xmod(row)<0.0) xmod(row)=2*pi+xmod(row)	
       endif
    enddo

    ! linear interpolation: nodes and weight
    do row=1,n2
       if(xmod(row)<cx(ni)) then
          do i=1,ni
             if(xmod(row)<cx(i)) then
                lint_ind(1,1,row)=i-1
                lint_ind(2,1,row)=i-1
                lint_ind(3,1,row)=i
                lint_ind(4,1,row)=i
                aux=(cx(i)-xmod(row))/(cx(i)-cx(i-1))
                exit
             end if
          end do
       else
          lint_ind(1,1,row)=ni
          lint_ind(2,1,row)=ni
          lint_ind(3,1,row)=1
          lint_ind(4,1,row)=1
          aux=(360.0_8-xmod(row))/(360.0_8-cx(ni))
       end if
       wt(1)=aux
       wt(2)=aux
       wt(3)=1.0_8-aux
       wt(4)=1.0_8-aux

       if(ymod(row)<cy(nj)) then
          do j=1,nj
             if(ymod(row)<cy(j)) then
                lint_ind(1,2,row)=j-1
                lint_ind(2,2,row)=j
                lint_ind(3,2,row)=j
                lint_ind(4,2,row)=j-1
                aux=(cy(j)-ymod(row))/(cy(j)-cy(j-1))
                exit
             end if
          end do
       else
          lint_ind(1,2,row)=nj
          lint_ind(2,2,row)=nj
          lint_ind(3,2,row)=nj
          lint_ind(4,2,row)=nj
          aux=1.0_8
       end if
       lint_weight(1,row)=wt(1)*aux
       lint_weight(2,row)=wt(2)*(1.0_8-aux)
       lint_weight(3,row)=wt(3)*(1.0_8-aux)
       lint_weight(4,row)=wt(4)*aux
    end do

    if(mype==0)  write(*,*) 'weights for interpolating forcing / 2D fields computed'
    return     
  end subroutine init_forcing_interp


  !------------------------------------------------------------------------
  subroutine forcing_linear_ip(zd,idim,jdim,ind,weights,zi,nmpt)
    !  this subroutine interpolates data using prepared weights- 
    !  see subroutine init_linear_forcing_interp
    !
    !  INPUT
    !        zd(idim,jdim)         available data set
    !        nmpt                  number of model positions, where data are wanted
    !        indx(4,2,nmpt)        i,j index of four neighbors of each node
    !        weights(4,nmpt)       interpolation weights
    !        
    !  OUTPUT
    !        zi(nmpt)              array of interpolated values for model points

    use g_parfe
    implicit none                                             

    integer      :: idim, jdim, nmpt                          
    integer      :: ind(4,2,nmpt)
    integer      :: i, n                      
    real(kind=8) :: zd(idim,jdim), zi(nmpt)
    real(kind=8) :: weights(4,nmpt)
    real(kind=8) :: fx

    do n=1,nmpt
       fx=0.0
       do i=1,4
          fx=fx+weights(i,n)*zd(ind(i,1,n),ind(i,2,n))
       enddo
       zi(n)=fx
    enddo

    return
  end subroutine forcing_linear_ip

end module g_forcing_interp
module g_diag
  ! diagnose flag
  implicit none
  save

  ! run state
  logical                                  :: check_run_state=.true.       ! use salinity to check blowup

  ! ocean
  logical                                  :: diag_oce=.true.
  logical                                  :: diag_oce_mix_layer=.true.
  logical                                  :: diag_oce_transp=.true.
  logical                                  :: diag_oce_GM_vel=.true.
  logical                                  :: diag_oce_SGS_transp=.true.
  
  ! ice
  logical                                  :: diag_ice=.true.
  ! forcing
  logical                                  :: diag_forcing=.true.

  namelist /diag_flag/ check_run_state, &
       diag_oce, diag_oce_mix_layer, diag_oce_transp, &
       diag_oce_GM_vel, diag_oce_SGS_transp, diag_ice, diag_forcing

end module g_diag
!
!------------------------------------------------------------------------------
!
module g_meanarrays
  ! mean and (mean)diagnose arrays
  implicit none
  save

  ! counter
  integer                                  :: meancounter

  ! mean of prediction variables

  ! ocean
  real(kind=8), allocatable, dimension(:,:):: tracermean
  real(kind=8), allocatable, dimension(:)  :: ufmean, sshmean
#ifndef use_non_hydrostatic
  real(kind=8), allocatable, dimension(:)  :: wrhsmean
#endif

  ! ice
  real(kind=8), allocatable, dimension(:)  :: a_ice_mean, m_ice_mean, m_snow_mean
  real(kind=8), allocatable, dimension(:)  :: u_ice_mean, v_ice_mean


  ! (mean) diagnose variables

  ! ocean
  real(kind=8), allocatable, dimension(:)  :: uTFmean, vTFmean
  real(kind=8), allocatable, dimension(:)  :: uSFmean, vSFmean
  real(kind=8), allocatable, dimension(:)  :: sgs_u, sgs_v
  real(kind=8), allocatable, dimension(:)  :: sgs_ut, sgs_vt
  real(kind=8), allocatable, dimension(:)  :: sgs_us, sgs_vs
  real(kind=8), allocatable, dimension(:)  :: mixlay_dep_mean

  ! ice
  real(kind=8), allocatable, dimension(:)  :: thdgr_mean, thdgrsn_mean
  real(kind=8), allocatable, dimension(:)  :: uhice_mean, vhice_mean
  real(kind=8), allocatable, dimension(:)  :: uhsnow_mean, vhsnow_mean
  real(kind=8), allocatable, dimension(:)  :: flice_mean

  ! forcing
  real(kind=8), allocatable, dimension(:)  :: tair_mean, shum_mean
  real(kind=8), allocatable, dimension(:)  :: uwind_mean, vwind_mean
  real(kind=8), allocatable, dimension(:)  :: rain_mean, snow_mean
  real(kind=8), allocatable, dimension(:)  :: runoff_mean, evap_mean
  real(kind=8), allocatable, dimension(:)  :: lwrd_mean, swrd_mean
  real(kind=8), allocatable, dimension(:)  :: qnet_mean, wnet_mean
  real(kind=8), allocatable, dimension(:)  :: olat_mean, osen_mean
  real(kind=8), allocatable, dimension(:)  :: olwout_mean
  real(kind=8), allocatable, dimension(:)  :: virtual_salt_mean, relax_salt_mean
  real(kind=8), allocatable, dimension(:)  :: stress_x_mean, stress_y_mean

end module g_meanarrays
!
!----------------------------------------------------------------------------

module o_mixing_kpp_mod
  ! Original numerical algorithm by Bill Large at NCAR, 1994
  ! Modified from MOM4 to FESOM by Qiang, Feb. 2011
  ! Equation numbers in the code refer to the Large etal paper. 

  use o_MESH
  use o_ELEMENTS
  use g_config
  use o_PARAM
  use o_array
  use g_forcing_arrays
  use g_PARFE

  implicit none

  private

  public oce_mixing_kpp_init
  public oce_mixing_kpp

  public hbl, ghats

  private bldepth
  private wscale
  private ri_iwmix
  private ddmix
  private blmix_kpp
  private enhance
  private cal_nodal_alpha_beta


  real, dimension(:), allocatable      :: bfsfc    ! surface buoyancy forcing    (m^2/s^3)
  real, dimension(:), allocatable      :: caseA    ! = 1 in case A; =0 in case B
  real, dimension(:), allocatable      :: stable   ! = 1 in stable forcing; =0 in unstable
  real, dimension(:,:), allocatable    :: dkm1     ! boundary layer diff_cbt at kbl-1 level
  real, dimension(:,:), allocatable    :: blmc     ! boundary layer mixing coefficients
  real, dimension(:), allocatable      :: Talpha   ! -d(rho)/ d(pot.temperature)  (kg/m^3/C)
  real, dimension(:), allocatable      :: Sbeta    ! d(rho)/ d(salinity)       (kg/m^3/PSU)
  real, dimension(:), allocatable      :: ustar    ! surface friction velocity       (m/s)
  real, dimension(:), allocatable      :: Bo       ! surface turb buoy. forcing  (m^2/s^3)
  real, dimension(:), allocatable      :: dbloc    ! local delta buoy at interfaces(m/s^2)
  real, dimension(:), allocatable      :: dbsfc    ! delta buoy w/ respect to sfc  (m/s^2)
  real, dimension(:), allocatable      :: dVsq     ! (velocity shear re sfc)^2   (m/s)^2

  integer, dimension(:), allocatable   :: kbl      ! index of first grid level below hbl

  real, dimension(:), allocatable      :: ghats    ! nonlocal transport (s/m^2)
  real, dimension(:), allocatable      :: hbl      ! boundary layer depth

  real, parameter :: epsln   = 1.0e-20 ! a small value

  real, parameter :: epsilon = 0.1
  real, parameter :: vonk    = 0.4
  real, parameter :: conc1   = 5.0

  real, parameter :: zmin    = -4.e-7  ! m3/s3 limit for lookup table of wm and ws
  real, parameter :: zmax    = 0.0     ! m3/s3 limit for lookup table of wm and ws
  real, parameter :: umin    = 0.0     ! m/s limit for lookup table of wm and ws
  real, parameter :: umax    = 0.04    ! m/s limit for lookup table of wm and ws

  real :: Ricr   = 0.3  ! critical bulk Richardson Number
  real :: concv  = 1.6  ! constant for pure convection (eqn. 23) (Large 1.5-1.6; MOM default 1.8)

  real :: cg            ! non-dimensional coefficient for counter-gradient term
  real :: Vtc           ! non-dimensional coefficient for velocity 
  ! scale of turbulant velocity shear        
  ! (=function of concv,concs,epsilon,vonk,Ricr)

  real :: deltaz        ! delta zehat in table
  real :: deltau        ! delta ustar in table

  integer, parameter :: nni = 890         ! number of values for zehat in the look up table
  integer, parameter :: nnj = 480         ! number of values for ustar in the look up table
  real, dimension(0:nni+1,0:nnj+1) :: wmt ! lookup table for wm, the turbulent velocity scale for momentum
  real, dimension(0:nni+1,0:nnj+1) :: wst ! lookup table for ws, the turbulent velocity scale scalars

  ! namelist parameters in the future: 
  ! concv, Ricr

contains


  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp_init">
  !
  ! Initialization for the KPP vertical mixing scheme
  !
  !     output:
  !       Vtc = non-dimensional constant used in calc. bulk Ri              
  !       cg  = constant used in calc.nonlocal transport term                
  !       wmt = turbulent velocity scale for momentum                         
  !       wst = turbulent velocity scale for scaler                          

  subroutine oce_mixing_kpp_init

    implicit none

    real, parameter :: cstar  = 10.0   ! proportionality coefficient for nonlocal transport
    real, parameter :: conam  = 1.257
    real, parameter :: concm  = 8.380 
    real, parameter :: conc2  = 16.0
    real, parameter :: zetam  = -0.2
    real, parameter :: conas  = -28.86
    real, parameter :: concs  = 98.96
    real, parameter :: conc3  = 16.0
    real, parameter :: zetas  = -1.0

    real :: zehat ! = zeta * ustar**3
    real :: zeta  ! = stability parameter d/L
    real :: usta

    integer :: i, j

    allocate (ghats(ToDim_nod3d))
    allocate (hbl(ToDim_nod2d))

    allocate (bfsfc(ToDim_nod2d))         ! surface buoyancy forcing    (m^2/s^3)
    allocate (caseA(ToDim_nod2d))         ! = 1 in case A; =0 in case B
    allocate (stable(ToDim_nod2d))        ! = 1 in stable forcing; =0 in unstable
    allocate (dkm1(ToDim_nod2d,3))        ! boundary layer diff at kbl-1 level
    allocate (blmc(ToDim_nod3d,3))        ! boundary layer mixing coefficients
    allocate (talpha(Todim_nod3d))        ! -d(rho)/ d(pot.temperature)  (g/m^3/C)
    allocate (sbeta(ToDim_nod3d))         ! d(rho)/ d(salinity)       (g/m^3/PSU)
    allocate (ustar(ToDim_nod2d))         ! surface friction velocity       (m/s)
    allocate (Bo(ToDim_nod2d))            ! surface turb buoy. forcing  (m^2/s^3)
    allocate (dbloc(ToDim_nod3d))         ! local delta buoy at interfaces(m/s^2)
    allocate (dbsfc(ToDim_nod3d))         ! delta buoy w/ respect to sfc  (m/s^2)
    allocate (dVsq(ToDim_nod3d))          ! (velocity shear re sfc)^2   (m/s)^2
    allocate (kbl(ToDim_nod2d))           ! index of first grid level below hbl

    ghats       = 0.0
    hbl         = 0.0

    !-----------------------------------------------------------------------
    !     initialize some constants for kmix subroutines, and initialize
    !     for kmix subroutine "wscale" the 2D-lookup table for wm and ws
    !     as functions of ustar and zetahat (=vonk*sigma*hbl*bfsfc).
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! define some non-dimensional constants  (recall epsilon=0.1)
    !-----------------------------------------------------------------------

    !     Vtc used in eqn. 23
    Vtc     = concv * sqrt(0.2/concs/epsilon) / vonk**2 / Ricr

    !     cg = cs in eqn. 20
    cg      = cstar * vonk * (concs * vonk * epsilon)**(1./3.)

    !-----------------------------------------------------------------------
    ! construct the wm and ws lookup tables (eqn. 13 & B1)
    !-----------------------------------------------------------------------

    deltaz = (zmax-zmin)/(nni+1) 
    deltau = (umax-umin)/(nnj+1)

    do i=0,nni+1
       zehat = deltaz*(i) + zmin
       do j=0,nnj+1
          usta = deltau*(j) + umin
          zeta = zehat/(usta**3+epsln)

          if(zehat >= 0.) then
             wmt(i,j) = vonk*usta/(1.+conc1*zeta)
             wst(i,j) = wmt(i,j)
          else
             if(zeta > zetam) then
                wmt(i,j) = vonk* usta * (1.-conc2*zeta)**(1./4.)
             else
                wmt(i,j) = vonk* (conam*usta**3-concm*zehat)**(1./3.)
             endif
             if(zeta > zetas) then
                wst(i,j) = vonk* usta * (1.-conc3*zeta)**(1./2.)
             else
                wst(i,j) = vonk* (conas*usta**3-concs*zehat)**(1./3.)
             endif
          endif
       enddo
    enddo

  end subroutine oce_mixing_kpp_init


  !#######################################################################
  ! <SUBROUTINE NAME="oce_mixing_kpp">
  !
  ! This subroutine computes the vertical diffusivity and viscosity according
  ! to the KPP scheme of Large etal.  In brief, the scheme does the 
  ! following:
  !
  ! --Compute interior mixing everywhere:                               
  !   interior mixing gets computed due to constant internal wave 
  !   background activity ("visc_cbu_iw" and "diff_cbt_iw").
  !   Mixing is enhanced in places of static instability (local Ri < 0).
  !   Additionally, mixing can be enhanced by contribution from shear 
  !   instability which is a function of the local Ri.
  !
  ! --Double diffusion:
  !   Interior mixing can be enhanced by double diffusion due to salt
  !   fingering and diffusive convection.
  !
  ! --Boundary layer:
  !
  !   (A) Boundary layer depth:
  !       at every gridpoint the depth of the oceanic boundary layer 
  !       ("hbl") gets computed by evaluating bulk richardson numbers.
  !
  !   (B) Boundary layer mixing:
  !       within the boundary layer, above hbl, vertical mixing is 
  !       determined by turbulent surface fluxes, and interior mixing at
  !       the lower boundary, i.e. at hbl.
  !
  ! outputs
  !
  !  hbl   = boundary layer depth (meters)
  !  ghats = nonlocal transport coefficient (s/m^2)
  !  viscA = viscosity coefficient  (m^2/s) 
  !  diffK = diffusion coefficient (m^2/s) 
  !
  subroutine oce_mixing_kpp(viscA, diffK)

    implicit none

    integer                             :: i, k, row, kn
    integer                             :: num, nodlo, nodup, ns, lay, lay_mi
    real                                :: smftu, smftv, aux
    real                                :: dens_up, minmix
    real                                :: usurf, vsurf, tsurf, ssurf
    real, dimension(:), intent(inout)   :: viscA !for momentum
    real, dimension(:,:), intent(inout) :: diffK !for T and S

    minmix=2.0e-3

    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)+1

       !-----------------------------------------------------------------------
       !     compute squared vertical difference of velocity
       !-----------------------------------------------------------------------

       usurf=0.5*(uf(nod3d_below_nod2d(1,i))+uf(nod3d_below_nod2d(2,i)))
       vsurf=0.5*(uf(nod3d_below_nod2d(1,i)+ToDim_nod3d) &
            +uf(nod3d_below_nod2d(2,i)+ToDim_nod3d))
       dVsq(nod3d_below_nod2d(1,i))=0.0
       do k=2, kn
          row=nod3d_below_nod2d(k,i)
          dVsq(row)=(usurf-uf(row))**2+(vsurf-uf(ToDim_nod3d+row))**2
       end do

       !-----------------------------------------------------------------------
       !     compute dbloc and dbsfc
       !     local buoyancy gradient at km interfaces:                  (m/s2)
       !           dbloc = g/rho{k+1,k+1} * [ drho{k,k+1}-drho{k+1,k+1} ]
       !     buoyancy difference with respect to "zref",i.e. the surface(m/s2)
       !           dbsfc = g * [ drho{1,k}/rho{1,k} - drho{k,k}/rho{k,k} ]
       !----------------------------------------------------------------------- 

       ns=nod3d_below_nod2d(1,i)
       tsurf=0.5*(tracer(ns,1)+tracer(nod3d_below_nod2d(2,i),1))
       ssurf=0.5*(tracer(ns,2)+tracer(nod3d_below_nod2d(2,i),2))
       dbsfc(ns)=0.0
       do k=2, kn
          nodup=nod3d_below_nod2D(k-1,i)
          nodlo=nod3D_below_nod2D(k,i)        
          call fcn_density(tracer(nodup,1),tracer(nodup,2),coord_nod3d(3,nodlo),dens_up)
          dbloc(nodup)=-g * (dens_up-density_insitu(nodlo))/density_insitu(nodlo)
          call fcn_density(tsurf,ssurf,coord_nod3d(3,nodlo),dens_up)
          dbsfc(nodlo)=-g * (dens_up-density_insitu(nodlo))/density_insitu(nodlo)
       end do
       dbloc(nodlo)=dbloc(nodup)
    end do

    !-----------------------------------------------------------------------
    !     compute thermal and haline expansion coefficients (without factor of rho).
    !     thermal expansion coefficient without 1/rho factor       (kg/m3/C)
    !           talpha= -d(rho{k,k})/d(T(k))
    !     salt expansion coefficient without 1/rho factor        (kg/m3/PSU)
    !           sbeta = d(rho{k,k})/d(S(k))
    !-----------------------------------------------------------------------

    call cal_nodal_alpha_beta

    !-----------------------------------------------------------------------
    ! friction velocity, turbulent sfc buoyancy forcing
    ! ustar = sqrt( sqrt( stress_x^2 + stress_y^2 ) / rho ) (m/s)
    ! bo =  -g * ( Talpha*heat_flux/vcpw + Sbeta * salinity*water_flux ) (m^2/s^3)
    !-----------------------------------------------------------------------

    do i=1,ToDim_nod2d
       row=nod3d_below_nod2d(1,i)
       ustar(i) = sqrt( sqrt(stress_x(i)**2 + stress_y(i)**2)*rho0r)
       Bo(i)    = -g * (Talpha(row) * heat_flux(i)/vcpw  &   !heat_flux & water_flux: positive up
            + Sbeta(row) * water_flux(i) * tracer(row,2))                                    
    enddo

    !-----------------------------------------------------------------------
    !     compute interior mixing coefficients everywhere, due to constant 
    !     internal wave activity, static instability, and local shear 
    !     instability.
    !-----------------------------------------------------------------------
    call ri_iwmix(viscA, diffK)

    !-----------------------------------------------------------------------
    !     add double diffusion
    !-----------------------------------------------------------------------

    if (double_diffusion) then
       call ddmix(diffK)
    endif

    !-----------------------------------------------------------------------
    !     boundary layer mixing coefficients: diagnose new b.l. depth
    !-----------------------------------------------------------------------

    call bldepth

    !-----------------------------------------------------------------------
    !     boundary layer diffusivities
    !-----------------------------------------------------------------------

    call blmix_kpp(viscA, diffK)

    !-----------------------------------------------------------------------
    !     enhance diffusivity at interface kbl - 1
    !-----------------------------------------------------------------------

    call enhance(viscA, diffK) 

    !-----------------------------------------------------------------------
    !     combine interior and b.l. coefficients and nonlocal term;
    !     add horizontal smoothing
    !-----------------------------------------------------------------------

    ! first do horizontal smoothing
    ! In structured mesh case: smooth diffusivities with a 5 point filter
    ! to remove 2 delta x noise
    ! In FEOM case: do smoothing over the cluster over each grid layer

    if (smooth_blmc) then 
       ! use dVsq and bfsfc to save memory
       dVsq=blmc(:,2)
       dbsfc=blmc(:,3)
       do i=1, myDim_nod2d  ! only over myDim_nod2d
          num=nghbr_nod2d(i)%nmb
          lay=minval(num_layers_below_nod2d(nghbr_nod2d(i)%addresses))+1
          lay_mi=maxval(kbl(nghbr_nod2d(i)%addresses))
          lay=min(lay_mi,lay)
          do k=1, lay
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),2))
             aux=aux+blmc(nod3d_below_nod2d(k,i),2)*(num-2)
             dVsq(nod3d_below_nod2d(k,i))=aux/2./(num-1)
             aux=sum(blmc(nod3d_below_nod2d(k,nghbr_nod2d(i)%addresses),3))
             aux=aux+blmc(nod3d_below_nod2d(k,i),3)*(num-2)
             dbsfc(nod3d_below_nod2d(k,i))=aux/2./(num-1)                   
          end do
       end do
       blmc(:,2)=dVsq
       blmc(:,3)=dbsfc
    end if

    ! then combine blmc and viscA/diffK

    do i=1, myDim_nod2d
       do k=1,num_layers_below_nod2d(i)+1
          if (k < kbl(i)) then
             !viscA(nod3d_below_nod2d(k,i))   = &
             !     max(Av0, blmc(nod3d_below_nod2d(k,i),1))
             !diffK(nod3d_below_nod2d(k,i),1) = &
             !     max(Kv0, blmc(nod3d_below_nod2d(k,i),2))
             !diffK(nod3d_below_nod2d(k,i),2) = &
             !     max(Kv0, blmc(nod3d_below_nod2d(k,i),3))
             viscA(nod3d_below_nod2d(k,i))   = &
                  max(viscA(nod3d_below_nod2d(k,i)), blmc(nod3d_below_nod2d(k,i),1))
             diffK(nod3d_below_nod2d(k,i),1) = &
                  max(diffK(nod3d_below_nod2d(k,i),1), blmc(nod3d_below_nod2d(k,i),2))
             diffK(nod3d_below_nod2d(k,i),2) = &
                  max(diffK(nod3d_below_nod2d(k,i),2), blmc(nod3d_below_nod2d(k,i),3))
          else
             ghats(nod3d_below_nod2d(k,i))=0.0
          endif
       enddo
       !do k=1, min(3,num_layers_below_nod2d(i)+1)
       !	  viscA(nod3d_below_nod2d(k,i)) = &
       !	max(minmix,viscA(nod3d_below_nod2d(k,i)))
       !enddo
    enddo

    ! Set the mixing coeff. in the first layer above some limiting value
    ! this is very helpful to avoid huge surface velocity when vertical
    ! viscosity is very small derived from the KPP scheme.
    ! I strongly recommend this trick, at least in the current FESOM version.	
    !minmix=1.0e-3
    !where(viscA(1:ToDim_nod2D)<minmix) 
    !   viscA(1:ToDim_nod2d)=minmix 
    !end where

    !-----------------------------------------------------------------------
    ! add non-local contribution will be added to tracer_rhs directly

  end subroutine oce_mixing_kpp


  !#######################################################################
  ! <SUBROUTINE NAME="bldepth">
  !
  ! The oceanic planetray boundary layer depth, hbl, is determined as
  ! the shallowest depth where the bulk richardson number is
  ! equal to the critical value, Ricr.
  !
  ! Bulk Richardson numbers are evaluated by computing velocity and
  ! buoyancy differences between values at level k and surface
  ! reference values.
  !
  ! In this configuration, the reference values are equal to the
  ! mean values between the first two levels.
  ! When using a very fine vertical grid, these values should be 
  ! computed as the vertical average of velocity and buoyancy from 
  ! the surface down to epsilon*zcoord(k).
  !
  ! When the bulk richardson number at k exceeds Ricr, hbl is
  ! linearly interpolated between grid levels zcoord(k) and zcoord(k-1).
  !
  ! The water column and the surface forcing are diagnosed for 
  ! stable/ustable forcing conditions, and where hbl is relative 
  ! to grid points (caseA), so that conditional branches can be 
  ! avoided in later subroutines.
  !
  !
  !  input
  !      real dbloc(t2d,nk)   = local delta buoyancy         (m/s^2)        
  !      real dbsfc(t2d,nk)   = delta buoyancy w/ respect to sfc (m/s)^2   
  !      real ustar(t2d)      = surface friction velocity     (m/s)      
  !      real dVsq(t3d)       = (velocity shear re sfc)^2   (m/s)^2
  !      real Bo(t2d)         = surface turbulent buoyancy forcing(m^2/s^3)  
  !      real sw_3d(t3d)      = radiative buoyancy forcing (C m/s)          
  !      real f_c(t2d)        = Coriolis parameter            (1/s)            
  !
  !  output
  !      real hbl(t2d)        ! boundary layer depth              (m)      
  !      real bfsfc(t2d)      ! Bo+radiation absorbed (m^2/s^3)      
  !      real stable(t2d)     ! =1 in stable forcing; =0 unstable          
  !      real caseA(t2d)      ! =1 in case A, =0 in case B                
  !      integer kbl(t2d)     ! index of first grid level below hbl        
  !
  subroutine bldepth

    real            :: Ritop, bvsq, Vtsq
    real            :: hekman, hmonob, hlimit
    real            :: Rib_km1, Rib_k, coeff_sw, zk, zkm1
    real            :: sigma,zehat, wm, ws, dzup, dzlo
    integer         :: i, k, nk, ns, nod, nodup, nodlo

    real, parameter :: cekman = 0.7  ! constant for Ekman depth
    real, parameter :: cmonob = 1.0  ! constant for Monin-Obukhov depth

    ! initialize hbl and kbl to bottomed out values
    kbl=num_layers_below_nod2d(:)+1
    hbl=abs(coord_nod3d(3,bt_nds))

    do i=1,ToDim_nod2d
#ifdef use_sw_pene      
       ns=nod3d_below_nod2d(1,i)
       coeff_sw=g*Talpha(ns)
#endif
       Rib_km1=0.0
       nk=num_layers_below_nod2d(i)+1
       bfsfc(i)=Bo(i) 

       do k=2,nk
          nod=nod3d_below_nod2d(k,i)
          nodup=nod3d_below_nod2d(k-1,i)
          zk=-coord_nod3d(3,nod)     
          zkm1=-coord_nod3d(3,nodup)

          ! bfsfc = Bo + sw contricution
#ifdef use_sw_pene
          bfsfc(i)=Bo(i) + coeff_sw*(sw_3d(ns)-sw_3d(nod))  ! sw_3d: [C m/s], positive downward
#endif
          stable(i) = 0.5 + sign( 0.5, bfsfc(i))
          sigma  = stable(i) + (1.0-stable(i)) * epsilon

          !-----------------------------------------------------------------------
          ! compute velocity scales at sigma, for z=abs(coord_nod3d(3,nod)):
          !-----------------------------------------------------------------------

          zehat= vonk * sigma * zk * bfsfc(i)
          call wscale(zehat, ustar(i), wm, ws)

          !-----------------------------------------------------------------------
          ! compute the turbulent shear contribution to Rib
          ! eqn. (23)
          !-----------------------------------------------------------------------

          dzup=zk-zkm1
          if(k<nk) then
             dzlo=-coord_nod3d(3,nod3d_below_nod2d(k+1,i))-zk
             bvsq =0.5*(dbloc(nodup) / dzup + dbloc(nod) / dzlo)
          else
             bvsq =dbloc(nodup)/dzup
          end if
          Vtsq = zk * ws * sqrt(abs(bvsq)) * Vtc

          !-----------------------------------------------------------------------
          !  compute bulk Richardson number at new level
          !  eqn. (21)
          !-----------------------------------------------------------------------

          Ritop = zk*dbsfc(nod)
          Rib_k = Ritop / (dVsq(nod) + Vtsq + epsln )

          if(Rib_k > Ricr) then
             ! linearly interpolate to find hbl where Rib = Ricr
             hbl(i) = zkm1 + dzup*(Ricr-Rib_km1)/(Rib_k-Rib_km1+epsln)
             kbl(i) = k
             exit
          else
             Rib_km1=Rib_k
          endif
       end do

       !-----------------------------------------------------------------------
       ! find stability and buoyancy forcing for boundary layer
       !-----------------------------------------------------------------------
#ifdef use_sw_pene
       ! Linear interpolation of sw_3d to depth of hbl
       bfsfc(i) =Bo(i) +  & 
            coeff_sw*(sw_3d(ns)-(sw_3d(nodup)+(sw_3d(nod)-sw_3d(nodup))*(hbl(i)-zkm1)/dzup))
       stable(i)=0.5 + sign( 0.5, bfsfc(i) )
       bfsfc(i) =bfsfc(i) + stable(i) * epsln  ! ensures bfsfc never=0
#endif

       !-----------------------------------------------------------------------
       !        check hbl limits for hekman or hmonob
       !        eqn. (24)
       !-----------------------------------------------------------------------
       if (bfsfc(i) > 0.0) then
          hekman = cekman * ustar(i) / max(abs(coriolis_param_nod2d(i)), epsln)
          hmonob = cmonob * ustar(i)*ustar(i)*ustar(i)     &
               /vonk / (bfsfc(i)+epsln) 
          hlimit = stable(i) * AMIN1(hekman,hmonob)
          hbl(i) = AMIN1(hbl(i), hlimit)
          hbl(i) = max(hbl(i), -coord_nod3d(3,nod3d_below_nod2d(2,i))) 
       endif

       !-----------------------------------------------------------------------
       !     find new kbl 
       !-----------------------------------------------------------------------
       kbl(i)=nk
       do k=2,nk
          if(-coord_nod3d(3,nod3d_below_nod2d(k,i)) > hbl(i)) then
             kbl(i) = k
             exit
          endif
       enddo

       !-----------------------------------------------------------------------
       !     find stability and buoyancy forcing for final hbl values
       !-----------------------------------------------------------------------      
#ifdef use_sw_pene
       ! Linear interpolation of sw_3d to depth of hbl
       nod=nod3d_below_nod2d(kbl(i),i)
       nodup=nod3d_below_nod2d(kbl(i)-1,i)
       bfsfc(i) =Bo(i) +  & 
            coeff_sw*(sw_3d(ns)-(sw_3d(nodup)+(sw_3d(nod)-sw_3d(nodup)) &
            *(hbl(i)+coord_nod3d(3,nodup)) &
            /(coord_nod3d(3,nodup)-coord_nod3d(3,nod))))
       stable(i)=0.5 + sign( 0.5, bfsfc(i) )
       bfsfc(i) =bfsfc(i) + stable(i) * epsln 
#endif

       !-----------------------------------------------------------------------
       !     determine caseA and caseB
       !     (if hbl is below (deeper than) the mid point of level kbl
       !     then caseA=0  else caseA=1)
       !-----------------------------------------------------------------------

       nod=nod3d_below_nod2d(kbl(i),i)
       nodup=nod3d_below_nod2d(kbl(i)-1,i)
       dzup=coord_nod3d(3,nodup)-coord_nod3d(3,nod)
       caseA(i)  = 0.5 + sign( 0.5, -coord_nod3d(3,nod)-0.5*dzup-hbl(i))

    enddo  !i

  end subroutine bldepth



  !#######################################################################
  ! <SUBROUTINE NAME="wscale">
  !
  ! Compute turbulent velocity scales.
  ! Use a 2D-lookup table for wm and ws as functions of ustar and
  ! zetahat (=vonk*sigma*hbl*bfsfc=zeta*ustar**3).
  !
  ! Note: the lookup table is only used for unstable conditions
  ! (zehat <= 0), in the stable domain wm (=ws) gets computed
  ! directly.
  !
  !  input
  !      real zehat    zeta*us**3
  !      real us       utar(nod), surface friction velocity    (m/s)    
  !  output                                                             
  !      real wm, ws   turbulent velocity scales 
  !
  subroutine wscale(zehat, us, wm, ws)

    implicit none

    real, intent(in)     :: zehat, us
    real, intent(out)    :: wm, ws 
    real                 :: zdiff, udiff, zfrac, ufrac, fzfrac
    real                 :: wam, wbm, was, wbs, u3
    integer              :: iz, izp1, ju, jup1

    !-----------------------------------------------------------------------
    !     use lookup table for zehat < zmax only; otherwise use
    !     stable formulae
    !-----------------------------------------------------------------------

    ! zehat = vonk * sigma * abs(depth) * bfsfc

    if (zehat <= zmax) then

       zdiff  = zehat-zmin
       iz = int( zdiff/deltaz )
       iz = min( iz , nni )
       iz = max( iz , 0  )
       izp1=iz+1

       udiff  = us-umin

       ju = int( min(udiff/deltau,float(nnj)))
       ju = max( ju , 0  )
       jup1=ju+1

       zfrac = zdiff/deltaz - float(iz)
       ufrac = udiff/deltau - float(ju)

       fzfrac= 1.-zfrac
       wam   = (fzfrac)  * wmt(iz,jup1) + zfrac*wmt(izp1,jup1)
       wbm   = (fzfrac)  * wmt(iz,ju  ) + zfrac*wmt(izp1,ju  )
       wm = (1.-ufrac)* wbm          + ufrac*wam

       was   = (fzfrac)  * wst(iz,jup1) + zfrac*wst(izp1,jup1)
       wbs   = (fzfrac)  * wst(iz,ju  ) + zfrac*wst(izp1,ju  )
       ws = (1.-ufrac)* wbs          + ufrac*was
    else
       u3    = us*us*us
       wm = vonk * us * u3 / ( u3 + conc1*zehat + epsln )
       ws = wm
    endif
  end subroutine wscale



  !#######################################################################
  ! <SUBROUTINE NAME="ri_iwmix">
  !
  ! Compute interior viscosity and diffusivity due 
  ! to shear instability (dependent on a local richardson number),
  ! to background internal wave activity, and 
  ! to static instability (local richardson number < 0).                    
  ! outputs:
  !    visc = viscosity coefficient (m**2/s)       
  !    diff = diffusion coefficient (m**2/s)     
  !
  subroutine ri_iwmix(viscA, diffK)

    implicit none

    integer           :: i, row, k, kn, nodup, nodlo, mr
    real, parameter   :: Riinfty = 0.8  ! local Richardson Number limit for shear instability
    real              :: dz, ri_prev, tmp, fx, dens_up
    real              :: Rigg, ratio, frit, fcont
    real, dimension(:), intent(out)   :: viscA
    real, dimension(:,:), intent(out) :: diffK

    logical :: smooth_richardson_number = .false.
    integer :: num_smoothings = 1 ! for vertical smoothing of Richardson number

    fx = -g*rho0r

    ! compute Richardson number and store it as viscA to save memory
    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)
       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)        
          call fcn_density(tracer(nodup,1),tracer(nodup,2),coord_nod3d(3,nodlo),dens_up)
          dz=coord_nod3d(3,nodup)-coord_nod3d(3,nodlo)
          viscA(nodup)=fx*dz*(dens_up-density_insitu(nodlo)) / &
               ((uf(nodup)-uf(nodlo))**2+(uf(ToDim_nod3d+nodup)-uf(ToDim_nod3d+nodlo))**2+epsln)
       end do

       ! smooth Richardson number in the vertical using a 1-2-1 filter
       if(smooth_richardson_number .and. kn>2) then
          do mr=1,num_smoothings
             ri_prev=0.25*viscA(nod3d_below_nod2d(1,i))
             do k=2,kn-1
                nodup=nod3d_below_nod2d(k,i)
                nodlo=nod3d_below_nod2d(k+1,i)
                tmp=viscA(nodup)
                viscA(nodup)=ri_prev+0.5*viscA(nodup)+0.25*viscA(nodlo)
                ri_prev=0.25*tmp
             end do
          end do
       end if
    end do

    ! compute viscA and diffK
    do row=1, ToDim_nod3D

       !-----------------------------------------------------------------------
       !           evaluate function of Ri# for shear instability eqn. (28b&c)
       !-----------------------------------------------------------------------

       Rigg  = AMAX1( viscA(row) , 0.0)
       ratio = AMIN1( Rigg/Riinfty , 1.0 )
       frit  = (1. - ratio*ratio)
       frit  = frit*frit*frit

       !-----------------------------------------------------------------------
       !           evaluate function of Ri# for convection eqn. (28a)
       !-----------------------------------------------------------------------

       fcont  = 0.5 * ( abs(viscA(row)) - viscA(row) )
       fcont  = fcont / (AMAX1(fcont,epsln))

       !-----------------------------------------------------------------------
       !           add mixing due to internal wave activity (equ 29),
       !           static instability and shear instability 
       !-----------------------------------------------------------------------

       viscA(row)     = Av0 + fcont*visc_conv_limit + visc_sh_limit*frit 
       diffK(row,1)   = Kv0 + fcont*diff_conv_limit + diff_sh_limit*frit 
       diffK(row,2)   = diffK(row,1)         

    enddo

  end subroutine ri_iwmix


  !#######################################################################
  ! <SUBROUTINE NAME="ddmix">
  !
  ! Rrho dependent interior flux parameterization.
  ! Add double-diffusion diffusivities to Ri-mix values at blending
  ! interface and below. 
  ! Salt fingering code modified july 2003
  ! by stephen.griffies@noaa.gov based on NCAR CCSM2.x
  !
  ! output: update diffu
  !
  subroutine ddmix(diffK)

    implicit none

    real, parameter :: Rrho0               = 1.9    ! limit for double diffusive density ratio
    real, parameter :: dsfmax              = 1.e-4  ! (m^2/s) max diffusivity in case of salt fingering
    real, parameter :: viscosity_molecular = 1.5e-6 ! (m^2/s)

    integer :: i, k, nodlo, nodup
    real    :: alphaDT, betaDS
    real    :: diffdd, Rrho, prandtl
    real, dimension(:,:), intent(inout) :: diffK

    do i=1,ToDim_nod2d
       do k=1, num_layers_below_nod2d(i)

          nodup=nod3d_below_nod2d(k,i)
          nodlo=nod3d_below_nod2d(k+1,i)

          alphaDT=0.5 * (Talpha(nodup)+Talpha(nodlo)) &
               * (tracer(nodup,1)-tracer(nodlo,1))
          betaDS =0.5 * (Sbeta(nodup)+Sbeta(nodlo)) &
               * (tracer(nodup,2)-tracer(nodlo,2))

          if (alphaDT > betaDS .and. betaDS > 0.) then

             ! salt fingering case,  eqn. (31)

             Rrho   = min(alphaDT/betaDS, Rrho0)

             !        diffdd = dsfmax*(1.0-((Rrho-1)/(Rrho0-1))**2)**pexp2  ! (very old code) 
             !        diffdd = 1.0-((Rrho-1)/(Rrho0-1))**2                  ! (less old code)
             diffdd = 1.0-((Rrho-1.0)/(Rrho0-1.0))                    ! (new code)
             diffdd = dsfmax*diffdd*diffdd*diffdd

             diffK(nodup,1) = diffK(nodup,1) + 0.7*diffdd !for temperature
             diffK(nodup,2) = diffK(nodup,2) + diffdd !for salinity

          else if ( alphaDT < 0. .and. alphaDT > betaDS ) then

             ! diffusive convection eqn. (32)

             Rrho    = alphaDT/betaDS 
             diffdd  = viscosity_molecular*0.909*exp(4.6*exp(-0.54*(1/Rrho-1.)))  

             ! eqn. (34)
             prandtl = 0.15*Rrho
             if (Rrho > 0.5) prandtl = (1.85-0.85/Rrho)*Rrho

             diffK(nodup,1) = diffK(nodup,1) + diffdd  !for temperature
             diffK(nodup,2) = diffK(nodup,2) + prandtl*diffdd  !for salinity

          endif

       enddo
    enddo

  end subroutine ddmix



  !#######################################################################
  ! <SUBROUTINE NAME="blmix_kpp">
  !
  ! Mixing coefficients within boundary layer depend on surface
  ! forcing and the magnitude and gradient of interior mixing below
  ! the boundary layer ("matching").
  !
  !     inputs:
  !
  !      real ustar(2d)    ! surface friction velocity         (m/s) 
  !      real bfsfc(2d)    ! surface buoyancy forcing     (m^2/s^3)  
  !      real hbl(2d)      ! boundary layer depth              (m)   
  !      real stable(2d)   ! = 1 in stable forcing                   
  !      real caseA(2d)    ! = 1 in case A                           
  !      integer kbl(2d)   ! index of first grid level below hbl
  !
  !     outputs:
  !
  !      real dkm1(2d,3) = boundary layer diff_cbt at kbl-1 level 
  !      real blmc(3d,3) = boundary layer mixing coeff.(m**2/s)   
  !      real ghats(3d)  = nonlocal scalar transport              
  !
  subroutine blmix_kpp(viscA,diffK)

    implicit none

    integer  :: i, nl, kn, knm1, knp1, ki, row
    integer  :: nod_col(1:max_num_layers)
    real     :: delhat, R, dvdzup, dvdzdn, wm, ws, sigma
    real     :: viscp, difsp, diftp, visch, difsh, difth, f1
    real     :: sig, a1, a2, a3, Gm, Gs, Gt
    real     :: zehat
    real     :: gat1m, gat1t, gat1s, dat1m, dat1s, dat1t
    real     :: z_col(1:max_num_layers), z_col_mid(1:max_num_layers)
    real     :: dthick(1:max_num_layers), diff_col(1:max_num_layers,3)
    real, dimension(:), intent(inout)   :: viscA
    real, dimension(:,:), intent(inout) :: diffK

    blmc = 0.0

    do i=1,ToDim_nod2d

       nl=num_layers_below_nod2d(i)+1

       if(nl<3) cycle  ! a temporary solution

       nod_col(1:nl)=nod3d_below_nod2d(1:nl,i)
       z_col(1:nl)=-coord_nod3d(3,nod_col(1:nl))
       z_col_mid(1:nl-1)=z_col(1:nl-1)+0.5*(z_col(2:nl)-z_col(1:nl-1))
       dthick(2:nl-1)=0.5*(z_col(3:nl)-z_col(1:nl-2))
       dthick(1)=dthick(2)
       dthick(nl)=dthick(nl-1)
       diff_col(1:nl-1,1)=viscA(nod_col(1:nl-1))
       diff_col(1:nl-1,2:3)=diffK(nod_col(1:nl-1),:)
       diff_col(nl,:)=diff_col(nl-1,:)

       !-----------------------------------------------------------------------
       !       compute velocity scales at hbl. (recall epsilon=0.1)
       !-----------------------------------------------------------------------
       sigma = stable(i) * 1.0 + (1.-stable(i)) * epsilon
       zehat= vonk * sigma * hbl(i) * bfsfc(i)
       call wscale(zehat, ustar(i), wm, ws)

       kn = ifix(caseA(i)+epsln) *(kbl(i) -1) +   &
            (1-ifix(caseA(i)+epsln)) * kbl(i)
       knm1 = max(kn-1,1)
       knp1 = min(kn+1,nl)

       !-----------------------------------------------------------------------
       !         find the interior viscosities and derivatives at hbl(i) 
       !         eqn. (18)
       !-----------------------------------------------------------------------

       delhat = z_col_mid(kn)-hbl(i)
       R      = 1.0 - delhat / dthick(kn)

       dvdzup = (diff_col(knm1,1) - diff_col(kn,1))/dthick(kn)
       dvdzdn = (diff_col(kn,1) - diff_col(knp1,1))/dthick(knp1)
       viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       dvdzup = (diff_col(knm1,3) - diff_col(kn,3))/dthick(kn)
       dvdzdn = (diff_col(kn,3) - diff_col(knp1,3))/dthick(knp1)     
       difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       dvdzup = (diff_col(knm1,2) - diff_col(kn,2))/dthick(kn)
       dvdzdn = (diff_col(kn,2) - diff_col(knp1,2))/dthick(knp1)
       diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+         &
            R  * (dvdzdn + abs(dvdzdn)) )

       visch  = diff_col(kn,1) + viscp * delhat
       difsh  = diff_col(kn,3) + difsp * delhat
       difth  = diff_col(kn,2) + diftp * delhat

       f1 = stable(i) * conc1 * bfsfc(i) / (ustar(i)**4+epsln) 

       gat1m = visch / (hbl(i)+epsln) / (wm+epsln)
       dat1m = -viscp / (wm+epsln) + f1 * visch
       dat1m = min(dat1m,0.) 

       gat1s = difsh  / (hbl(i)+epsln) / (ws+epsln)
       dat1s = -difsp / (ws+epsln) + f1 * difsh 
       dat1s = min(dat1s,0.) 

       gat1t = difth /  (hbl(i)+epsln) / (ws+epsln)
       dat1t = -diftp / (ws+epsln) + f1 * difth 
       dat1t = min(dat1t,0.) 

       do ki=1,nl-1

          if (ki >= kbl(i)) exit
          row=nod_col(ki)

          !-----------------------------------------------------------------------
          !         compute turbulent velocity scales on the interfaces
          !-----------------------------------------------------------------------

          sig   = z_col_mid(ki) / (hbl(i)+epsln)
          sigma = stable(i)*sig                          &
               + (1.-stable(i))*AMIN1(sig,epsilon)
          zehat= vonk * sigma * hbl(i) * bfsfc(i)
          call wscale(zehat, ustar(i), wm, ws)

          !-----------------------------------------------------------------------
          !         compute the dimensionless shape functions at the interfaces
          !         eqn. (11)
          !-----------------------------------------------------------------------

          a1 = sig - 2.
          a2 = 3.-2.*sig
          a3 = sig - 1.

          Gm = a1 + a2 * gat1m + a3 * dat1m
          Gs = a1 + a2 * gat1s + a3 * dat1s
          Gt = a1 + a2 * gat1t + a3 * dat1t

          !-----------------------------------------------------------------------
          !          compute boundary layer diffusivities at the interfaces
          !          eqn. (10)
          !-----------------------------------------------------------------------

          blmc(row,1) = hbl(i) * wm * sig * (1.+sig*Gm) 
          blmc(row,2) = hbl(i) * ws * sig * (1.+sig*Gt) 
          blmc(row,3) = hbl(i) * ws * sig * (1.+sig*Gs) 

          !-----------------------------------------------------------------------
          !          nonlocal transport term = ghats * <ws>o (eqn. 20)
          !-----------------------------------------------------------------------

          ghats(row) = (1.-stable(i)) * cg    &
               / (ws * hbl(i) + epsln)

       enddo

       !-----------------------------------------------------------------------
       !     find diffusivities at kbl-1 grid level 
       !-----------------------------------------------------------------------

       sig   =  z_col(kbl(i)-1)  / (hbl(i)+epsln)
       sigma =stable(i) * sig                  &
            + (1.-stable(i)) * min(sig,epsilon)
       zehat= vonk * sigma * hbl(i) * bfsfc(i)
       call wscale(zehat, ustar(i), wm, ws)

       a1= sig - 2.
       a2 = 3.-2.*sig
       a3 = sig - 1.
       Gm = a1 + a2 * gat1m + a3 * dat1m
       Gs = a1 + a2 * gat1s + a3 * dat1s
       Gt = a1 + a2 * gat1t + a3 * dat1t
       dkm1(i,1) = hbl(i) * wm * sig * (1. + sig * Gm)
       dkm1(i,2) = hbl(i) * ws * sig * (1. + sig * Gt)
       dkm1(i,3) = hbl(i) * ws * sig * (1. + sig * Gs)

    enddo

  end subroutine blmix_kpp


  !#######################################################################
  ! <SUBROUTINE NAME="enhance">
  !
  ! Enhance the diffusivity at the kbl-.5 interface
  !
  ! input
  !      integer kbl(n2)   =  grid above hbl                     
  !      real hbl(n2)      =  boundary layer depth (m)           
  !      real dkm1(n2,3)   =  bl diffusivity at kbl-1 grid level 
  !      real caseA(n2)    =  1 in caseA, = 0 in case B
  !
  ! input/output
  !      real ghats(ij_bounds,nk)  =  nonlocal transport     (s/m**2)
  !      modified ghats at kbl(i)-1 interface        
  ! output
  !      real blmc(n3,3) = enhanced boundary layer mixing coefficient
  !
  subroutine enhance(viscA, diffK)
    implicit none

    integer :: i, k, nodk, nodkbl
    real    :: delta, dkmp5, dstar
    real, dimension(:), intent(inout)   :: viscA
    real, dimension(:,:), intent(inout) :: diffK

    do i=1,ToDim_nod2d
       k = kbl(i) - 1
       nodk=nod3d_below_nod2d(k,i)
       nodkbl=nod3d_below_nod2d(kbl(i),i)

       delta = (hbl(i)+coord_nod3d(3,nodk))/(coord_nod3d(3,nodk)-coord_nod3d(3,nodkbl))

       ! momentum
       dkmp5 = caseA(i) * viscA(nodk)      &
            + (1.-caseA(i)) * blmc(nodk,1)
       dstar = (1.-delta)**2 * dkm1(i,1) + delta**2 * dkmp5
       blmc(nodk,1) = (1.-delta) * viscA(nodk)  &
            + delta * dstar

       ! temperature:
       dkmp5 = caseA(i) * diffK(nodk,1)  &
            + (1.-caseA(i)) * blmc(nodk,2)
       dstar = (1.-delta)**2 * dkm1(i,2) + delta**2 * dkmp5    
       blmc(nodk,2) = (1.-delta) * diffK(nodk,1)  &
            + delta * dstar

       ! salinity:   
       dkmp5 = caseA(i) * diffK(nodk,2)  &
            + (1.-caseA(i)) * blmc(nodk,3)
       dstar = (1.-delta)**2 * dkm1(i,3) + delta**2 * dkmp5
       blmc(nodk,3) = (1.-delta) * diffK(nodk,2)  &
            + delta * dstar

       ghats(nodk) = (1.-caseA(i)) * ghats(nodk)
    enddo

  end subroutine enhance
  !
  !----------------------------------------------------------------------------
  !
  subroutine cal_nodal_alpha_beta
    !   A function to calculate the thermal expansion coefficient
    !   and saline contraction coefficient. 
    !
    ! REFERENCE:
    !    McDougall, T.J. 1987.  Neutral Surfaces
    !    Journal of Physical Oceanography, vol 17, 1950-1964,
    !
    ! INPUT:
    !   tracer(:,1) = potential temperature [degree C (ITS-90)]
    !   tracer(:,2) = salinity              [psu      (PSS-78)]
    !   z           = pressure (or -depth)  [db]
    !
    ! OUTPUT:
    !   Talpha = Thermal expansion coeff (alpha) [degree_C.^-1]
    !   Sbeta  = Saline contraction coeff (beta) [psu.^-1]
    !
    ! Qiang Wang, 25,11,2004
    ! Modified to compute nodal values for KPP, Feb. 2011, Qiang
    !-----------------------------------------------------------------
    ! CHECK VALUE:
    !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
    !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
    !-----------------------------------------------------------------

    implicit none

    integer    :: i, ni
    real       :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2 
    real       :: a_over_b    

    ! we should ensure that we are evoluting potential temperature in the model.

    ! cycle over node
    do i = 1, ToDim_nod3d
       ! prepare values that will be used below
       t1 = tracer(i,1)*1.00024_8
       s1 = tracer(i,2)
       p1 = abs(coord_nod3D(3,i)) 

       t1_2 = t1*t1
       t1_3 = t1_2*t1
       t1_4 = t1_3*t1
       p1_2 = p1*p1
       p1_3 = p1_2*p1
       s35 = s1-35.0_8
       s35_2 = s35*s35

       ! calculate beta
       Sbeta(i) = 0.785567e-3 - 0.301985e-5*t1 &
            + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
            + s35*(-0.356603e-6 + 0.788212e-8*t1 &
            + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
            + s35_2*(0.515032e-8) & 
            + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
            + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
            + p1_3*(0.121551e-17)

       ! calaculate the thermal expansion / saline contraction ratio
       a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
            - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
            - 0.255019e-7*t1_4 &
            + s35*(0.378110e-2 - 0.846960e-4*t1 &
            - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
            + s35_2*(-0.678662e-5) &
            + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
            + p1_2*t1_2*(0.512857e-12) &
            - p1_3*(0.302285e-13)

       ! calculate alpha
       Talpha(i) = a_over_b*Sbeta(i)
    end do
  end subroutine cal_nodal_alpha_beta

end module o_mixing_KPP_mod

module o_mixing_pp_mod
  !Richardson number dependent Av and Kv following
  !Pacanowski and Philander, 1981

  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use i_array
  use g_config
  use g_PARFE

  implicit none

  private

  public oce_mixing_pp
 
contains

  subroutine oce_mixing_pp
    !Compute Richardson number dependent Av and Kv following
    !Pacanowski and Philander, 1981
    !Av = c * factor**2 + Av0,
    !Kv = Av * factor + Kv0, 
    !factor=1/(1+10Ri),  Ri=N**2/(dU/dz)**2 is the Richardson number
    !                    N is buoyancy frequency
    !                   dU/dz is vertical velocity shear
    !c is a tunable parameter (named mix_coeff_PP in the code)
    !Kv0 ~ 1.e-5, Av0 ~ 1.e-4

    implicit none
    !
    integer     	:: col,n,j,node_low,node_up, node, k, n3
    real   	:: density_up, drho, g_rhoinv
    real	        :: dz_inv, dzu, dzv, zm
    real	        :: velocity_shear, dzb, factor
    real          :: mixlength
    logical     :: mo_on 

    g_rhoinv=g*rho0r 
    n3=ToDim_nod3d

    do col=1,ToDim_nod2d

       n=num_layers_below_nod2D(col)

       do j=1,n       
          node=nod3D_below_nod2D(j,col)
          node_up = node
          node_low = nod3D_below_nod2D(j+1,col)

          dz_inv=1./(coord_nod3d(3,node_up)-coord_nod3d(3,node_low))

          call fcn_density(tracer(node_up,1),tracer(node_up,2), &
               coord_nod3d(3,node_low),density_up)

          ! density vert. grad.
          drho = density_up - density_insitu(node_low)

          if(drho>=0.0) then
             factor=1.0 
          else
             ! BV frequency squared
             dzb=drho*dz_inv*g_rhoinv  
             ! velocity shear
             dzu=(uf(node_up)- uf(node_low)) * dz_inv
             dzv=(uf(node_up+n3) - uf(node_low+n3)) * dz_inv 
             velocity_shear=dzu*dzu + dzv*dzv
             factor=velocity_shear/( velocity_shear - 5.0 * dzb)      
          end if

          Av(node)=PP_max_mix_coeff*(factor**2)
          Kv(node,1)=Av(node)*factor

          ! add background mixing
          Av(node)=Av(node)+Av0
          Kv(node,1)=Kv(node,1)+Kv0

          ! turn on convection if unstable 
          if(drho>0.0) then
             if(allow_convect_global) then
                Kv(node,1)=diff_conv_limit
                Av(node)=visc_conv_limit
             elseif(coord_nod2d(2,col)>0.0) then     
                Kv(node,1)=diff_conv_limit
                Av(node)=visc_conv_limit
             endif
          end if

          ! apply the tb04 scheme
#ifdef use_ice
#ifdef use_cavity
          mo_on=(add_tb04_to_PP .and. cavity_flag_nod2d(col)==0)
#else
          mo_on=add_tb04_to_PP
#endif
          if(mo_on) then
             call mo_length(water_flux(col),heat_flux(col), &
                  stress_x(col),stress_y(col), &
                  u_ice(col),v_ice(col),a_ice(col), &
                  dt, mixlength)
             zm=0.5*(coord_nod3d(3,node_up)+coord_nod3d(3,node_low))
             if(abs(zm)<=mixlength) then
                Av(node)=Av(node)+modiff
                Kv(node,1)=Kv(node,1)+modiff
             end if
          end if
#endif

       end do
    end do

    ! approximation for high freq. wind mixing near the surface
    ! (in case not in the forcing)
    where(Av(1:ToDim_nod2D)<wndmix) 
       Av(1:ToDim_nod2d)=wndmix 
    end where
    where(Kv(1:ToDim_nod2D,1)<wndmix) 
       Kv(1:ToDim_nod2d,1)=wndmix 
    end where

  end subroutine oce_mixing_PP
  !
  !---------------------------------------------------------------------------
  !
  subroutine mo_length(water_flux,heat_flux,stress_x,stress_y,  &
       u_ice,v_ice,a_ice,dt,mixlength)
    ! vertical mixing scheme of Timmermann and Beckmann, 2004.
    ! computes the mixing length derived from the Monin-Obukhov length
    ! Ralph Timmermann, 14.06.2006

    implicit none

    real*8              :: water_flux, heat_flux, stress_x, stress_y
    real*8              :: u_ice, v_ice, a_ice, uabs
    real*8              :: dt, ret, rtc, mixlength
    real*8              :: qfm, qtm, qw
    real*8              :: ustar,tau, obuk
    real(kind=8), parameter :: cosgam = 0.913632  ! cos(24.*3.14/180.)

    qfm            = water_flux * 34.       ! note that water_flux>0
    ![psu * m/s]   [m/s]   [psu]    ! f. upward fresh water flux

    qtm            = - 2.38e-7 * heat_flux  ! heat_flux>0 f. upward heat flux
    ![K * m/s]

    tau = sqrt(stress_x**2+stress_y**2)
    ustar = sqrt(tau/1030.)
    uabs = sqrt(u_ice**2+v_ice**2)

    qw = 1.25 * ustar**3 * (1-a_ice) + 0.005 * uabs**3 * cosgam * a_ice  !Eq. 8 of TB04

    call pmlktmo(qfm,qtm,qw,obuk)

    rtc=10.0*86400.0/dt     ! time constant of mixed layer retreat

    if (obuk.lt.mixlength) then
       ret=(obuk-mixlength)/rtc
       mixlength=mixlength+ret
    else
       mixlength=obuk
    endif

  end subroutine mo_length
  !
  !=========================================================================	    	    
  !
  subroutine pmlktmo(qfm,qtm,qw,obuk)
    ! gives the Monin-Obukhov length
    ! qtm  = Heat Flux into ML                                 [K m/s]
    ! qfm  = salinity flux into ML                             [psu m/s]
    ! qw   = production of turbulent kinetic energy 
    !-----------------------------------------------------------------------
    implicit none

    integer           :: iter
    real*8            :: qtm, qfm, qw, obuk
    real*8, parameter :: qhw   = 7.0                           ! [m]
    real*8, parameter :: betas = 0.0008
    real*8, parameter :: betat = 0.00004
    real*8            :: a1, f0, f1, ttmp, qrho

    qrho=betas*qfm-betat*qtm
    ttmp=60.

    if(qrho>0.) then
       ttmp=0.
    else
       do iter=1,5
          a1 = exp(-ttmp/qhw) 
          f0 = 2.* qw * a1 + 9.81 * qrho * ttmp
          f1 =-2.* qw * a1 / qhw + 9.81 * qrho
          ttmp = ttmp - f0 / f1
          ttmp = max(ttmp,10.) 
       enddo
    end if

    obuk=max(ttmp,10.)

  end subroutine pmlktmo

end module o_mixing_pp_mod


module o_mixing_my2p5_mod
  ! Mellor-Yamada 2.5 mixing scheme (Mellor and Yamada, 1982)
  ! The Cantha and Clayson (1994) stability functions are used.
  ! Constrains on turbulent scales following Galperin etal (1988) are enforced.
  ! Here we took a simple/approximate way to deal with advection and diffusion.
  !
  ! It was coded following the structure of POM with changes in detail. Further
  ! modification was then undertaken to use the Cantha and Clayson stability 
  ! functions. I am still not very sure about the treatment of the upper and 
  ! bottom boundaries. A check-through for it is required at a later stage.
  !
  ! The background values (Kv0 and Av0) should be added when operating vertical
  ! mixing for the momentum and tracer equations. They are traditionally set to
  ! Kv0~1e-5 and Av0~1e-4. Define them in the namelist file.
  !
  ! Qiang, 06.08.2010

  use o_mesh
  use o_elements
  use o_matrices
  use o_param
  use o_array 
  use g_config
  use g_parfe

  implicit none

  public oce_mixing_my2p5_init
  public oce_mixing_my2p5

  real, allocatable, dimension(:) 	  :: Kq
  real, allocatable, dimension(:) 	  :: q2, q2b, q2l, q2lb
  real, allocatable, dimension(:) 	  :: rhs_q2, rhs_q2l

contains

  subroutine oce_mixing_my2p5_init
    ! initialize the MY2p5 scheme
    implicit none

    integer	    :: m, row, k, kb, nodbot, nodk
    real(kind=8)    :: l, smallvalue

    data smallvalue /1.0e-8/

    allocate(Kq(ToDim_nod3d))
    allocate(q2(ToDim_nod3d), q2b(ToDim_nod3d))
    allocate(q2l(Todim_nod3d), q2lb(ToDim_nod3d))
    allocate(rhs_q2(ToDim_nod3d), rhs_q2l(ToDim_nod3d))
    do row=1,ToDim_nod2d
       kb=num_layers_below_nod2d(row)+1
       nodbot=bt_nds(row)
       l=0.1*(coord_nod3d(3,row)-coord_nod3d(3,nodbot))
       do k=1,kb
          nodk=nod3d_below_nod2d(k,row)
          q2b(nodk)=smallvalue
          q2lb(nodk)=l*q2b(nodk)
          Kv(nodk,1)=l*sqrt(q2b(nodk))
          Av(nodk)=Kv(nodk,1)
          Kq(nodk)=Kv(nodk,1)
          q2(nodk)=q2b(nodk)
          q2l(nodk)=q2lb(nodk)
       end do
    end do

  end subroutine oce_mixing_my2p5_init
  !
  !------------------------------------------------------------------------------
  !
  subroutine oce_mixing_my2p5
    ! update MY2.5

    implicit none

    integer   	:: m, row, n, i, j, elnodes(4), elem, elem2
    integer	:: k, ki, kbm1, kb, lay, n3, rowd, elem_type
    integer	:: nodbot, nodk, nodkm1, nodkp1
    integer	:: nod_col(max_num_layers)
    real(kind=8)        :: dt2, stress_sur, stress_bot, ustr, vstr
    real(kind=8)        :: aux, dudz, dvdz, v, vtr, mass
    real(kind=8)	:: dx(4), dy(4), dz(4), qdx, qdy, qdz
    real(kind=8)	:: um, vm, wm, u_el(4), v_el(4), q_el(4)
    real(kind=8)	:: Ah, d1, d2, S(3)
    real(kind=8)	:: rotate_coe, temp_coe, temp_coe2, temp_coe3
    real(kind=8)	:: zb, z0, z, tp, sp, p, kn, Gh
    real(kind=8)        :: dens_k, grad_rho, grad_rho_up
    real(kind=8)	:: a1, a2, b1, b2, c1, c2, c3
    real(kind=8)	:: e1, e2, l0, stf, smallvalue
    real(kind=8)	:: cbcnst, surfl, sef, smoth, inv4
    real(kind=8)        :: coef1, coef2, coef3, coef4, coef5
    real(kind=8)	:: a(max_num_layers), c(max_num_layers)
    real(kind=8)	:: ee(max_num_layers), gg(max_num_layers)
    real(kind=8)	:: rhs_q2_col(max_num_layers), rhs_q2l_col(max_num_layers)
    real(kind=8)        :: dz_col(max_num_layers), l(max_num_layers)
    real(kind=8)        :: boygr(max_num_layers), cc(max_num_layers)
    real(kind=8)        :: prod(max_num_layers), q_col(max_num_layers)
    real(kind=8)	:: dtef(max_num_layers)
    real(kind=8)	:: sh(max_num_layers), sm(max_num_layers), maxgrad
    real(kind=8)	:: kappa !von Karman's constant

    data kappa /0.4/
    data a1,b1,a2,b2,c1,c2,c3 /0.92, 16.6, 0.74, 10.1, 0.08, 0.7, 0.2/
    data e1/1.8/, e2/1.33/
    data cbcnst/100./, surfl/2.e5/, sef/1./

    data smoth/0.05/, inv4/0.25/
    data smallvalue /1.0e-8/

    dt2=2.0*dt

    coef1=a2*(1.0-6.0*a1/b1)
    coef2=3.0*a2*(6.0*a1+b2*(1.0-c3))
    coef3=a1*(1.0-3.0*c1-6.0*a1/b1)
    coef4=a1*9.0*(2.0*a1+a2*(1-c2))
    coef5=9.0*a1*a2

    n3=myDim_nod3D+eDim_nod3D              

    rhs_q2=0.0
    rhs_q2l=0.0

    ! assemble rhs for q2 and q2l (adv. and hori. diff.)
    ! use leapfrog for adv, Euler for diff.
    ! use lumped mass matrix.
    do elem=1, myDim_elem3d        
       elem2=elem2d_corresp_to_elem3d(elem)
       elem_type=grid_type_elem2d(elem2)
       elnodes=elem3D_nodes(:,elem)
       dx=bafux_3D(:,elem)
       dy=bafuy_3D(:,elem)
       dz=bafuz_3D(:,elem)
       v=voltetra(elem)

       ! mean velocity  
       u_el=uf(elnodes)
       v_el=uf(elnodes+n3)                       
       um=sum(u_el)*inv4
       vm=sum(v_el)*inv4
#ifndef use_non_hydrostatic
       wm=sum(w(elnodes)*dz)
#else
       wm=sum(uf(elnodes+2*n3))*inv4              
#endif

       ! horizontal viscosity 
       if(smagorinsky_visc) then
          vtr=voltriangle(elem2)
          d1=sum(u_el*dx)-sum(v_el*dy)
          d2=sum(u_el*dy)+sum(v_el*dx)
          Ah=sqrt(d1*d1+d2*d2)*vtr*4.0
          Ah=max(Ah, 40.0/1.0e8*vtr)              !limit from below: 20m^2/s on 10km
          if(Ah*dt/vtr > 0.05) Ah=0.05*vtr*dt_inv !limit from above
       else
          Ah=Ah0
          if(scale_mixing_h) Ah=Ah*(voltriangle(elem2)/scalevol)**(1.0/real(scale_mixing_type))
       endif

       ! rhs for q2
       ! adv
       q_el=q2(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       rhs_q2(elnodes)=rhs_q2(elnodes) - (um*qdx+vm*qdy+wm*qdz)*inv4*v
       ! diff 
       q_el=q2b(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       if(elem_type==1) then  !sigma diff.
          lay=elem3d_layer(elem)
          S(1)=grid_slope(1,lay,elem2)
          S(2)=grid_slope(2,lay,elem2)
          aux=S(1)**2+S(2)**2
          S(3)=sqrt(aux)
          rotate_coe=1.0/(1.0+aux) 

          !diagonal part1 (lateral)
          temp_coe=Ah*rotate_coe*v
          rhs_q2(elnodes)=rhs_q2(elnodes)-(dx*qdx*(1.0+S(2)*S(2))+ &
               dy*qdy*(1.0+S(1)*S(1)))*temp_coe
          !diagonal part2 (cross slope)
          temp_coe2=S(3)*S(3)*temp_coe
          rhs_q2(elnodes)=rhs_q2(elnodes)-dz*qdz*temp_coe2
          !off diagonal part1 (lateral) --> (1,3),(2,3)
          rhs_q2(elnodes)=rhs_q2(elnodes)-(S(1)*dx+S(2)*dy)*qdz*temp_coe
          !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
          temp_coe3=S(1)*S(2)*temp_coe
          rhs_q2(elnodes)=rhs_q2(elnodes) + (dx*qdy+dy*qdx)*temp_coe3    
          !off diagonal part2 (cross slope) --> (3,1),(3,2)
          rhs_q2(elnodes)=rhs_q2(elnodes)-(S(1)*qdx+S(2)*qdy)*dz*temp_coe
       else
          rhs_q2(elnodes)=rhs_q2(elnodes) - Ah*(dx*qdx+dy*qdy)*v
       endif

       ! rhs for q2l
       ! adv
       q_el=q2l(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       rhs_q2l(elnodes)=rhs_q2l(elnodes) - (um*qdx+vm*qdy+wm*qdz)*inv4*v
       ! diff 
       q_el=q2lb(elnodes)
       qdx=sum(dx*q_el)
       qdy=sum(dy*q_el)
       qdz=sum(dz*q_el)
       if(elem_type==1) then  !sigma diff.
          !diagonal part1 (lateral)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(dx*qdx*(1.0+S(2)*S(2))+ &
               dy*qdy*(1.0+S(1)*S(1)))*temp_coe
          !diagonal part2 (cross slope)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-dz*qdz*temp_coe2
          !off diagonal part1 (lateral) --> (1,3),(2,3)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(S(1)*dx+S(2)*dy)*qdz*temp_coe
          !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
          rhs_q2l(elnodes)=rhs_q2l(elnodes) + (dx*qdy+dy*qdx)*temp_coe3    
          !off diagonal part2 (cross slope) --> (3,1),(3,2)
          rhs_q2l(elnodes)=rhs_q2l(elnodes)-(S(1)*qdx+S(2)*qdy)*dz*temp_coe
       else
          rhs_q2l(elnodes)=rhs_q2l(elnodes) - Ah*(dx*qdx+dy*qdy)*v
       endif
    end do ! elem

    ! solve for rhs_q2 and rhs_q2l
    rhs_q2 = q2b + dt2*rhs_q2/uv_lump    
    rhs_q2l = q2lb + dt2*rhs_q2l/uv_lump

    ! add other terms in the q2 and q2l equations
    do row=1,myDim_nod2D                    

       kbm1=num_layers_below_nod2D(row)
       kb=kbm1+1

       ! nodes and dz in this colume
       nod_col(1)=nod3d_below_nod2d(1,row)
       do k=2, kb
          nod_col(k)=nod3d_below_nod2d(k,row)
          dz_col(k-1)=coord_nod3d(3,nod_col(k-1))-coord_nod3d(3,nod_col(k))
       end do

       ! surface and bottom depth and nodes
       rowd=nod_col(1) 
       z0=coord_nod3d(3,rowd)  
       nodbot=nod3d_below_nod2d(kb,row)
       zb=coord_nod3d(3,nodbot)

       !rhs_q2 and rhs_q2l in this colume
       rhs_q2_col(1:kb)=rhs_q2(nod_col(1:kb))
       rhs_q2l_col(1:kb)=rhs_q2l(nod_col(1:kb))

       !-----------------------------------------------------------------------
       !Surface and bottom boundary conditions

       ! amplitude of surface and bottom stress
       if(cavity_flag_nod2d(row)==0) then       
          stress_sur=sqrt(stress_x(row)**2+stress_y(row)**2)*rho0r
       else            
          aux=sqrt(uf(rowd)**2+uf(n3+rowd)**2)       
          ustr=-C_d*uf(rowd)*aux                      
          vstr=-C_d*uf(n3+rowd)*aux                 
          stress_sur=sqrt(ustr**2+vstr**2)
       end if
       aux=sqrt(uf(nodbot)**2+uf(n3+nodbot)**2)      
       ustr=-C_d*uf(nodbot)*aux
       vstr=-C_d*uf(n3+nodbot)*aux                   
       stress_bot=sqrt(ustr**2+vstr**2)        

       ee(1)=0.0
       aux=(16.6**(2.0/3.0))*sef
       if(cavity_flag_nod2d(row)==0) then      
          ! Wave breaking energy- a variant of Craig & Banner (1994)
          ! see Mellor and Blumberg, 2004.
          gg(1)=(15.8*cbcnst)**(2.0/3.0)*stress_sur 
          ! Surface length scale following Stacey (1999).
          l0=surfl*stress_sur/g
       else  !under the cavity
          gg(1)=stress_sur*aux   ! similar to ocean bottom
          l0=0.0
       end if
       rhs_q2_col(kb)=stress_bot*aux     


       !-----------------------------------------------------------------------
       ! The following section solves the equation:
       ! dti2*(kq*q2')' - q2*(2.*dti2*dtef+1.) = -q2b

       ! coefficients a, c
       do k=2, kbm1
          nodk=nod_col(k)
          nodkm1=nod_col(k-1)
          nodkp1=nod_col(k+1)
          aux=dz_col(k)+dz_col(k-1)
          a(k)=-dt2*(Kq(nodk)+Kq(nodkp1)+2.0*Kv0)/(dz_col(k)*aux)
          c(k)=-dt2*(Kq(nodkm1)+Kq(nodk)+2.0*Kv0)/(dz_col(k-1)*aux)
       end do

       ! the original way of computing buoyancy frequency in POM
       ! Calculate speed of sound:
       !do k=2,kbm1
       !   nodk=nod_col(k)
       !   tp=TF(nodk)
       !   sp=SF(nodk)
       !   ! pressure in units of decibars
       !   ! p=abs(coord_nod3d(3,nodk))*g*rho0*1.0e-4
       !   p=abs(coord_nod3d(3,nodk))
       !   cc(k)=1449.1+0.00821*p+4.55*tp-0.045*tp**2+1.34*(sp-35.0)
       !   cc(k)=cc(k)/sqrt((1.0-0.01642*p/cc(k))*(1.0-0.4*p/cc(k)**2))
       !end do
       !
       ! Calculate buoyancy gradient:
       !do k=2,kbm1
       !   nodk=nod_col(k)
       !   nodkm1=nod_col(k-1)
       !   nodkp1=nod_col(k+1)
       !   boygr(k)=g*0.5*((density_insitu(nodkm1)-density_insitu(nodk))/dz_col(k-1) + &
       !	(density_insitu(nodk)-density_insitu(nodkp1))/dz_col(k))*rho0r
       !   if(.not.density_linear) boygr(k)=boygr(k)+(g**2)/(cc(k)**2)       
       !end do

       ! my way of calculating buoyancy frequency:
       nodk=nod_col(1)
       nodkp1=nod_col(2)
       call fcn_density(tracer(nodk,1), tracer(nodk,2), coord_nod3d(3,nodkp1), dens_k)
       grad_rho_up=(dens_k-density_insitu(nodkp1))/dz_col(1)
       boygr(1)=g*rho0r*grad_rho_up
       do k=2,kbm1
          nodk=nod_col(k)
          nodkp1=nod_col(k+1)
          call fcn_density(tracer(nodk,1), tracer(nodk,2), coord_nod3d(3,nodkp1), dens_k)
          grad_rho=(dens_k-density_insitu(nodkp1))/dz_col(k)
          boygr(k)=g*rho0r*0.5*(grad_rho+grad_rho_up)
          grad_rho_up=grad_rho
       end do
       boygr(kb)=g*rho0r*grad_rho_up
       ! boygr(1) and boygr(kb) are only required for limiting l

       ! Calculate production of turbulent kinetic energy:
       do k=2,kbm1
          nodk=nod_col(k)
          nodkm1=nod_col(k-1)
          nodkp1=nod_col(k+1)
          dudz=(((uf(nodkm1)-uf(nodk))/dz_col(k-1))**2 + &
               ((uf(nodk)-uf(nodkp1))/dz_col(k))**2)*0.5
          dvdz=(((uf(n3+nodkm1)-uf(n3+nodk))/dz_col(k-1))**2 + &   
               ((uf(n3+nodk)-uf(n3+nodkp1))/dz_col(k))**2)*0.5      
          prod(k)=Av(nodk)*sef*(dudz+dvdz)  
          prod(k)=prod(k) + Kv(nodk,1)*boygr(k)       
       end do

       ! length scale at time step  n
       l(1)=kappa*l0
       do k=2,kbm1
          nodk=nod_col(k)
          q2b(nodk)=max(q2b(nodk), smallvalue)    !necessary to use q2b and q2lb
          q2lb(nodk)=max(q2lb(nodk),smallvalue)   !to derive l (not q2 and q2l)
          l(k)=q2lb(nodk)/q2b(nodk)           
          if((z0-coord_nod3d(3,nodk))<0.5*(z0-zb)) l(k)=max(l(k), l(1))  
       end do
       l(kb)=0.0

       !dissipation rate: dtef*q2
       do k=2, kbm1
          nodk=nod_col(k)
          dtef(k)=sqrt(q2(nodk))/(b1*l(k)+smallvalue)
       end do

       ! compute ee, gg
       ! ee(1) gg(1) was set boundary condition; ee(kb) gg(kb) is not required.		
       do k=2,kbm1
          gg(k)=1.0/(a(k)+c(k)*(1.0-ee(k-1))-2.0*dt2*dtef(k)-1.0)
          ee(k)=a(k)*gg(k)
          gg(k)=(-2.0*dt2*prod(k)-rhs_q2_col(k)+c(k)*gg(k-1))*gg(k)
       end do

       ! compute q2 at n+1
       ! rhs_q2_col(kb) was set boundary condition
       do k=1,kbm1
          ki=kb-k
          rhs_q2_col(ki)=ee(ki)*rhs_q2_col(ki+1)+gg(ki)
       end do

       !-----------------------------------------------------------------------
       ! The following section solves the equation:
       ! dti2(kq*q2l')' - q2l*(dti2*dtef+1.) = -q2lb

       !boundary 
       ee(1)=0.0  
       gg(1)=0.0
       rhs_q2l_col(kb)=0.0
       !ee(2)=0.0
       !gg(2)=kappa*(z0-coord_nod3d(3,nod_col(2)))*q2(nod_col(2))
       !rhs_q2l_col(kbm1)=kappa*(coord_nod3d(3,nod_col(kbm1))-zb)*q2(nod_col(kbm1))

       ! dissipation rate
       do k=2,kbm1 
          z=coord_nod3d(3,nod_col(k))
          dtef(k)=dtef(k)*(1.0+e2*((1.0/abs(z-z0)+1.0/abs(z-zb))*l(k)/kappa)**2)   
       end do

       ! compute ee, gg 
       ! ee(1) gg(1) was set boundary condition; ee(kb) gg(kb) is not required.	
       do k=2, kbm1 !3
          gg(k)=1.0/(a(k)+c(k)*(1.0-ee(k-1))-dt2*dtef(k)-1.0)
          ee(k)=a(k)*gg(k)
          gg(k)=(-dt2*prod(k)*l(k)*e1-rhs_q2l_col(k)+c(k)*gg(k-1))*gg(k)
       end do

       ! compute q2l at n+1
       ! rhs_q2l_col(kb) was set boundary condition	
       do k=1,kb-1  !kb-2        
          ki=kb-k
          rhs_q2l_col(ki)=ee(ki)*rhs_q2l_col(ki+1)+gg(ki)
       end do

       ! avoid negative values of q2 and q2l
       do k=1,kb
          if(rhs_q2_col(k)<smallvalue .or. rhs_q2l_col(k)<smallvalue) then
             rhs_q2_col(k)=smallvalue
             rhs_q2l_col(k)=smallvalue
          end if
       end do

       !-----------------------------------------------------------------------
       ! The following section solves for Av and Kv:

       ! new length scale
       l(2:kbm1)=rhs_q2l_col(2:kbm1)/rhs_q2_col(2:kbm1)
       ! constrain l according to Galperin etal. 1988
       q_col(1:kb)=sqrt(rhs_q2_col(1:kb))
       do k=1,kb
          if(boygr(k)<0.0) then
             l(k)=min(l(k), 0.53*q_col(k)/(sqrt(-boygr(k))+smallvalue))
          end if
       end do

       ! stability functions: Cantha and Clayson, 1994
       do k=1,kb
          Gh=l(k)**2/rhs_q2_col(k)*boygr(k)
          Gh=min(max(Gh,-0.28), 0.029)
          sh(k)=coef1/(1.0-coef2*Gh)
          sm(k)=coef3+sh(k)*coef4*Gh
          sm(k)=sm(k)/(1.0-coef5*Gh)
       end do

       ! POM note: There are 2 options for Kq which, unlike Av and Kv, 
       ! was not derived by Mellor and Yamada but was purely
       ! empirical based on neutral boundary layer data.
       ! The choice is whether or not it should be subject to
       ! the stability factor, sh. Generally, there is not a great
       ! difference in output. 
       ! Q:The first was mentioned/used by Mellor and Blumberg, 2004; 
       ! the second was used by Cantha and Clayson 1994)

       ! POM takes the mean of the current step and the last step for
       ! Av, Kv and Kq. I only take the new value as I did no see the 
       ! reason to do that.
       ! POM uses old l and q2 to compute Gh and kn; I decided to
       ! use the new l and q2 as I did not see the reason to do that.
       ! If we recognize what we did is wrong, modification is required here!

       do k=1,kb
          nodk=nod_col(k)
          kn=l(k)*q_col(k)
          !Kq(nodk)=kn*0.41*sh(k)
          Kq(nodk)=kn*0.2
          Av(nodk)=kn*sm(k)
          Kv(nodk,1)=kn*sh(k)
       end do

       !-----------------------------------------------------------------------
       ! Average Kv and Av to values at the mid-edges (between every two nodes)
       ! because currently the implicit vertical mixing equations require such values. 
       ! Storage convention: store a mean value to the upper node location of the edge
       do k=1,kbm1
          nodk=nod_col(k)
          Av(nodk)=0.5*(Av(nodk)+Av(nod_col(k+1)))
          Kv(nodk,1)=0.5*(Kv(nodk,1)+Kv(nod_col(k+1),1))
       end do

       !----------------------------------------------------------------------- 
       ! updating arrays q2 q2b q2l q2lb and filtering (Robert-Asselin time filter)
       do k=1,kb
          nodk=nod_col(k)
          q2(nodk)=q2(nodk)+0.5*smoth*(rhs_q2_col(k)+q2b(nodk)-2.0*q2(nodk))
          q2b(nodk)=q2(nodk)
          q2(nodk)=rhs_q2_col(k)
          q2l(nodk)=q2l(nodk)+0.5*smoth*(rhs_q2l_col(k)+q2lb(nodk)-2.0*q2l(nodk))
          q2lb(nodk)=q2l(nodk)
          q2l(nodk)=rhs_q2l_col(k)
       end do

    end do

    ! communicate partition neighbours
    call com_3d(Av)
    call com_3d(Kv(:,1))
    call com_3d(q2)
    call com_3d(q2b)
    call com_3d(q2l)
    call com_3d(q2lb)

  end subroutine oce_mixing_MY2p5

end module o_mixing_MY2p5_mod
module o_mixing_tidal_mod
  ! Note: currently only barotropic part is considered.
  ! 
  ! Compute vertical diffusivity and vertical viscosity
  ! deduced from barotropic and baroclinic tidal 
  ! dissipation. For the baroclinic dissipation, we follow
  ! Simmons etal, and for the barotropic dissipation we follow 
  ! Lee etal. Assume Prandtl number unity. 
  !
  ! Reference:
  !
  ! Simmons, Jayne, St. Laurent, and Weaver, 2004:
  ! Tidally driven mixing in a numerical model of the ocean 
  ! general circulation.  Ocean Modelling, vol. 6,
  ! pages 245-263.
  !
  ! Jayne and St. Laurent, 2001:
  ! Parameterizing tidal dissipation over rough topography.
  ! Geophysical Research Letters, vol. 28, pages 811-814.
  !
  ! Hyun-Chul Lee, A. Rosati, and M.J. Spelman, 2006: 
  ! Barotropic tidal mixing effects in a coupled climate model:
  ! ocean conditions in the northern Atlantic
  ! Ocean Modelling, vol 11, pages 464--477
  !
  ! Osborn, T.R., 1980: Estimates of the local rate of vertical diffusion 
  ! from dissipation measurements.  JPO, vol. 10, pages 83-89.
  !
  ! Munk and Anderson, 1948: Notes on a theory of the thermocline. 
  ! Journal of Marine Research, vol 3. pages 276-295.

  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_rotate_grid
  use g_parfe

  implicit none

  private

  public oce_mixing_tidal 
  public oce_mixing_tidal_init
  private vert_mix_drag 
  private compute_bvfreq

  real, dimension(:), allocatable :: tide_speed      ! tide speed (m/s) from barotropic tide model 
  real, dimension(:), allocatable :: bvfreq          ! buoyancy frequency (sec^-1) 

  real    :: munk_anderson_p          = 0.25     ! (dimensionless) from Munk and Anderson 
  real    :: munk_anderson_sigma      = 3.0      ! (dimensionless) from Munk and Anderson 

  real    :: von_karman               = 0.4
  real    :: bott_drag_cd             = 2.4e-3   ! (dimensionless) bottom drag from Lee etal

  real    :: background_diffusivity   = 0.1e-4   ! (m^2/sec)
  real    :: background_viscosity     = 0.1e-4   ! (m^2/sec)

  real    :: drhodz_min               = 1.e-10   ! (kg/m^4) minimum abs(drhodz) used to compute N^2 
  real    :: tide_speed_min           = 5.e-3    ! (m/s) below which set a mask=0 for drag mixing diffusivity 

  logical :: smooth_bf                = .true.   ! for smoothing N in vertical with 1-2-1
  integer :: num_121_smoothings       = 1        ! number of 1-2-1 smoothing 

contains


  !#######################################################################
  ! Initialization for the oce_mixing_tidal module.
  subroutine oce_mixing_tidal_init 

    integer					:: i, j, fid
    integer					:: num_lon_reg, num_lat_reg
    integer                                     :: p_flag
    real(kind=8)					:: lon, lat
    real(kind=8), allocatable, dimension(:)	:: lon_reg_4u, lat_reg_4u
    real(kind=8), allocatable, dimension(:)	:: lon_reg_4v, lat_reg_4v
    real(kind=8), allocatable, dimension(:)	:: xn2d, yn2d
    real, allocatable, dimension(:,:)		:: amp_reg
    real, allocatable, dimension(:)		:: tide_u_amp, tide_v_amp

    allocate(bvfreq(ToDim_nod3d))
    bvfreq = 0.0  

    allocate(tide_speed(ToDim_nod2D))
    tide_speed = default_tide_speed

    ! read tidal speed (m/s) from a tide model, such as the
    ! Global Inverse Solution TPX07.1 created by OSU.

    if(read_tide_speed) then 

       open(101,file=trim(TideForcingPath)//'lonlat_4U_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       read(101,*) num_lon_reg, num_lat_reg 
       allocate(lon_reg_4u(num_lon_reg), lat_reg_4u(num_lat_reg))
       do i=1,num_lon_reg
          read(101,*) lon_reg_4u(i)
       end do
       do i=1,num_lat_reg
          read(101,*) lat_reg_4u(i)
       end do
       close(101)
       open(102,file=trim(TideForcingPath)//'lonlat_4V_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       read(102,*) num_lon_reg, num_lat_reg 
       allocate(lon_reg_4v(num_lon_reg), lat_reg_4v(num_lat_reg))
       do i=1,num_lon_reg
          read(102,*) lon_reg_4v(i)
       end do
       do i=1,num_lat_reg
          read(102,*) lat_reg_4v(i)
       end do
       close(102)

       ! allocate arrays for global data on regular grids
       allocate(amp_reg(num_lon_reg, num_lat_reg)) 
       allocate(tide_u_amp(ToDim_nod2d))
       allocate(tide_v_amp(ToDim_nod2d))
       allocate(xn2d(ToDim_nod2d))
       allocate(yn2d(ToDim_nod2d))

       do i=1, ToDim_nod2d
          if(rotated_grid) then
             call r2g(lon, lat, coord_nod2d(1,i), coord_nod2d(2,i))
             xn2d(i)=lon/rad   ! change unit to degree
             yn2d(i)=lat/rad
          else
             xn2d(i)=coord_nod2d(1,i)/rad
             yn2d(i)=coord_nod2d(2,i)/rad
          end if
          ! change lon range to [0 360]
          if(xn2d(i)<0.) xn2d(i)=xn2d(i) + 360.0  
       end do

       fid=103
       open(fid,file=trim(TideForcingPath)//'U_'//Tmix_tidalconstituent//'_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       do i=1, num_lon_reg
          do j=1, num_lat_reg
             read(fid, *) amp_reg(i,j)         
          end do
       end do
       close(fid)
       p_flag=0
       call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, amp_reg, &
            ToDim_nod2d, xn2d, yn2d, tide_u_amp, p_flag)
       tide_u_amp=tide_u_amp/abs(coord_nod3d(3,bt_nds))  ! change transport (m^2/s) to depth mean velocity (m/s)

       open(fid,file=trim(TideForcingPath)//'V_'//Tmix_tidalconstituent//'_'//trim(Tmix_tidalmodelname)//'.dat', status='old')
       do i=1, num_lon_reg
          do j=1, num_lat_reg
             read(fid, *) amp_reg(i,j)         
          end do
       end do
       close(fid)
       p_flag=0
       call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, amp_reg, &
            ToDim_nod2d, xn2d, yn2d, tide_v_amp, p_flag)
       tide_v_amp=tide_v_amp/abs(coord_nod3d(3,bt_nds))  ! change transport (m^2/s) to depth mean velocity (m/s)

       tide_speed=sqrt(tide_u_amp*tide_u_amp + tide_v_amp*tide_v_amp)

       deallocate(amp_reg, tide_u_amp, tide_v_amp)
       deallocate(xn2d, yn2d)
       deallocate(lat_reg_4v,lon_reg_4v,lat_reg_4u,lon_reg_4u)

    else
       tide_speed=default_tide_speed
    endif

  end subroutine oce_mixing_tidal_init


  !#######################################################################
  !
  ! Compute vertical diffusivity and viscosity
  ! based on one or both of the following dissipation mechanisms:
  !
  ! 1. internal wave breaking as parameterized by Simmons etal.
  !
  ! 2. barotropic tides feeling the bottom drag, as parameterized by 
  !    Lee etal.  
  !
  ! Note: only the second part is coded. Part 1 need to be updated when required.

  subroutine oce_mixing_tidal(viscA, diffK)

    real, dimension(:,:), intent(inout) :: diffK
    real, dimension(:),   intent(inout) :: viscA

    call compute_bvfreq

    !if(use_wave_dissipation)  call vert_mix_wave(viscA, diffK)
    if(use_drag_dissipation)  call vert_mix_drag(viscA, diffK)

    diffK = diffK + background_diffusivity
    viscA = viscA + background_viscosity

  end subroutine oce_mixing_tidal


  !#####################################################################
  !
  ! Compute Brunt-Vaisala (buoyancy) frequency. 
  !
  subroutine compute_bvfreq

    real    :: dz, aux, bf_prev, tmp, fx, dens_up
    integer :: i, k, kn, mr, nodup, nodlo

    fx=g*rho0r
    do i=1,ToDim_nod2d
       kn=num_layers_below_nod2D(i)
       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)        
          call fcn_density(tracer(nodup,1),tracer(nodup,2),coord_nod3d(3,nodlo),dens_up)
          dz=coord_nod3d(3,nodup)-coord_nod3d(3,nodlo)
          aux=-(dens_up-density_insitu(nodlo)) / dz
	  bvfreq(nodup)=sqrt(fx*max(aux, drhodz_min))
       end do

       if(smooth_bf .and. kn>2) then
          do mr=1,num_121_smoothings
             bf_prev=0.25*bvfreq(nod3d_below_nod2d(1,i))
             do k=2,kn-1
                nodup=nod3d_below_nod2d(k,i)
                nodlo=nod3d_below_nod2d(k+1,i)
                tmp=bvfreq(nodup)
                bvfreq(nodup)=bf_prev+0.5*bvfreq(nodup)+0.25*bvfreq(nodlo)
                bf_prev=0.25*tmp
             end do
          end do
       end if
    end do

  end subroutine compute_bvfreq


  !#####################################################################
  !
  ! Computes tracer diffusivity based on the methods of Lee etal., 
  ! which consider the dissipation from barotropic tides
  ! rubbing against the ocean bottom. 
  !
  ! Assume a unit Prandtl number.
  !
  subroutine vert_mix_drag(viscA, diffK)

    integer :: i, k, kn, nodup, nodlo
    real    :: bottdep, speedr, fx, height, ri, diff_drag
    real, dimension(:,:), intent(inout) :: diffK
    real, dimension(:),   intent(inout) :: viscA

    fx=von_karman/sqrt(bott_drag_cd)

    do i=1,ToDim_nod2d
       if(tide_speed(i)<tide_speed_min) cycle
       kn=num_layers_below_nod2D(i)
       bottdep=-coord_nod3d(3,bt_nds(i))
       speedr=fx/tide_speed(i)

       do k=1, kn
          nodup=nod3d_below_nod2D(k,i)
          nodlo=nod3D_below_nod2D(k+1,i)
          height=bottdep + &
               0.5*(coord_nod3d(3,nodup)+coord_nod3d(3,nodlo))      
          ri = 2.0*(bvfreq(nodup)*height*speedr)**2

          ! compute drag induced diffusivity 
          ! (Lee etal equations (1), (2), and (3))
          diff_drag=max_drag_diffusivity &
               *(1.0 + munk_anderson_sigma*ri)**(-munk_anderson_p)
          diffK(nodup,:)=diffK(nodup,:)+diff_drag
          viscA(nodup)=viscA(nodup)+diff_drag
       enddo
    enddo

  end subroutine vert_mix_drag

end module o_mixing_tidal_mod
module o_age_tracer_mod
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_clock
  use g_parfe
  implicit none

  integer, allocatable, dimension(:)   :: index_age_tracer  
  integer, allocatable, dimension(:,:) :: age_tracer_loc_index

contains


  subroutine age_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cageind
    character(4)         :: tr_name
    character(100)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_age_tracer(num_age_tracer))
    do j=1, num_age_tracer
       write(cageind,'(i1)') j
       tr_name='age'//cageind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_age_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! set initial values
    do j=1, num_age_tracer
       tracer(:,index_age_tracer(j))=0.0
    end do

    !-------------------------------------------------------------- 
    ! restore time scale at the release region
    if(zero_age_at_release) age_tracer_restore_time=dt

    !--------------------------------------------------------------    	
    ! set age tracer location index: 1 at release, 0 otherwise

    allocate(age_tracer_loc_index(ToDim_nod3d,num_age_tracer))
    age_tracer_loc_index=0

    allocate(temp_arr2d(nod2d))
    temp_arr2d=0
    do n=1, ToDim_nod2D
       temp_arr2d(myList_nod2D(n))=n
    end do

    do j=1, num_age_tracer
       write(cageind,'(i1)') j
       tr_name='age'//cageind
       file_name=trim(meshpath)//'age_tracer_release_nodes_'//tr_name//'.out'
       fileID=160
       open(fileID, file=file_name)
       read(fileID,*) num_nod
       allocate(nodes_release(num_nod))
       read(fileID,*) nodes_release
       close(fileID)
       do n=1,num_nod
          n_loc=temp_arr2d(nodes_release(n))
          if(n_loc>0) then
             n_loc=nod3d_below_nod2d(1,n_loc)
             age_tracer_loc_index(n_loc,j)=1
          end if
       end do
       deallocate(nodes_release)
    end do

    deallocate(temp_arr2d)

    !--------------------------------------------------------------
    ! in case release in volume

    if(age_release_in_volume) then
       do i=1,ToDim_nod2d
          row=nod3d_below_nod2d(1,i)
          do j=1, num_age_tracer
             if(age_tracer_loc_index(row,j)==1) then
                do k=2,num_layers_below_nod2d(i)+1
                   n3=nod3d_below_nod2d(k,i)
                   age_tracer_loc_index(n3,j)=1
                end do
             end if
          end do
       end do
    end if

  end subroutine age_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_tendency

    integer      :: j, elem, elnodes(4), res_ind_elem(4)
    real(kind=8) :: inv20, vol, source_elem(4), sum_source
    real(kind=8) :: inside_val(4), outside_val

    inv20=1./20.
    outside_val=1.0/(86400.0*(365+fleapyear))

    do elem=1,myDim_elem3d              
       elnodes=elem3D_nodes(:,elem)
       vol=voltetra(elem)*inv20   
       do j=1, num_age_tracer
          inside_val=-tracer(elnodes,index_age_tracer(j))/age_tracer_restore_time
          res_ind_elem=age_tracer_loc_index(elnodes,j) ! 1-inside, 0-outside
          source_elem=inside_val*res_ind_elem + outside_val*(1.0-res_ind_elem)
          sum_source=sum(source_elem)
          tracer_rhs(elnodes,index_age_tracer(j))=tracer_rhs(elnodes,index_age_tracer(j)) &
               + (sum_source+source_elem(elnodes))*vol
       end do
    end do

  end subroutine age_tracer_tendency
  !
  !-------------------------------------------------------------------------
  !
  subroutine age_tracer_cutoff_restore

    integer   :: j, row

    do j=1, num_age_tracer
       do row=1,ToDim_nod3d
          tracer(row,index_age_tracer(j))= &
               max(tracer(row,index_age_tracer(j)),0.)
          if(zero_age_at_release .and. age_tracer_loc_index(row,j)==1) then
             tracer(row,index_age_tracer(j))=0.0
          end if
       end do
    end do

  end subroutine age_tracer_cutoff_restore

end module o_age_tracer_mod
module o_passive_tracer_mod
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer, allocatable, dimension(:)        :: index_passive_tracer
  integer, allocatable, dimension(:,:)      :: passive_tracer_loc_index
  real(kind=8), allocatable, dimension(:,:) :: ptr_sfc_force

contains


  subroutine passive_tracer_init

    integer              :: i, j, k, n, n3, row, fileID
    integer              :: n_loc, num_nod
    integer, allocatable :: temp_arr2d(:), nodes_release(:)
    character(1)         :: cptrind
    character(4)         :: tr_name
    character(100)       :: file_name

    !--------------------------------------------------------------  
    ! find index
    allocate(index_passive_tracer(num_passive_tracer))
    do j=1, num_passive_tracer
       write(cptrind,'(i1)') j
       tr_name='ptr'//cptrind
       do i=1, num_tracer
          if(prog_tracer_name(i) == tr_name) index_passive_tracer(j)=i
       end do
    end do

    !--------------------------------------------------------------  
    ! initial values
    do j=1, num_age_tracer
       tracer(:,index_passive_tracer(j))=ptr_background_value
    end do

    !--------------------------------------------------------------
    ! in case that p.tr is restored in a region
    if(passive_tracer_restore) then

       ! set passive tracer location index: 1 at release, 0 otherwise

       allocate(passive_tracer_loc_index(ToDim_nod3d,num_passive_tracer))
       passive_tracer_loc_index=0

       allocate(temp_arr2d(nod2d))
       temp_arr2d=0
       do n=1, ToDim_nod2D
          temp_arr2d(myList_nod2D(n))=n
       end do

       do j=1, num_passive_tracer
          write(cptrind,'(i1)') j
          tr_name='ptr'//cptrind
          file_name=trim(meshpath)//'passive_tracer_restore_nodes_'//tr_name//'.out'
          fileID=160
          open(fileID, file=file_name)
          read(fileID,*) num_nod
          allocate(nodes_release(num_nod))
          read(fileID,*) nodes_release
          close(fileID)
          do n=1,num_nod
             n_loc=temp_arr2d(nodes_release(n))
             if(n_loc>0) then
                n_loc=nod3d_below_nod2d(1,n_loc)
                passive_tracer_loc_index(n_loc,j)=1
                tracer(n_loc,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
          deallocate(nodes_release)
       end do

       deallocate(temp_arr2d)

       !--------------------------------------------------------------
       ! in case restore volume
       if(ptr_restore_in_volume) then
          do i=1,ToDim_nod2d
             row=nod3d_below_nod2d(1,i)
             do j=1, num_passive_tracer
                if(passive_tracer_loc_index(row,j)==1) then
                   do k=2,num_layers_below_nod2d(i)+1
                      n3=nod3d_below_nod2d(k,i)
                      passive_tracer_loc_index(n3,j)=1
                      tracer(n3,index_passive_tracer(j))=ptr_restore_value
                   end do
                end if
             end do
          end do
       end if

    end if

    !--------------------------------------------------------------
    ! in case that passive tracers enter through surface fluxes
    if(passive_tracer_flux) then
       allocate(ptr_sfc_force(ToDim_nod2d,num_passive_tracer))
    end if

  end subroutine passive_tracer_init
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_sfc_bc

    integer      :: elem, j, elnodes(3), elnodes2(3)
    real(kind=8) :: auxf, entries(3)

    if(passive_tracer_flux) then

       ptr_sfc_force=0.0

       do elem=1,myDim_elem2d             
          elnodes2=elem2D_nodes(:,elem)
          elnodes=nod3D_below_nod2D(1,elnodes2)   
          auxf=voltriangle(elem)/12.0_8    
          do j=1,num_passive_tracer
             entries=-auxf*(tracer(elnodes,2)+tracer(elnodes,index_passive_tracer(j))) &
                  * runoff_landice(elnodes2)*landice_season(month)
             ptr_sfc_force(elnodes2,j)=ptr_sfc_force(elnodes2,j)+sum(entries)+entries
          end do
       end do

    end if

  end subroutine ptr_sfc_bc
  !
  !-------------------------------------------------------------------------
  !
  subroutine ptr_cutoff_restore

    integer   :: j, row

    if(passive_tracer_restore) then
       do j=1, num_passive_tracer
          do row=1,ToDim_nod3d
             if(passive_tracer_loc_index(row,j)==1) then
                tracer(row,index_passive_tracer(j))=ptr_restore_value
             end if
          end do
       end do
    endif

  end subroutine ptr_cutoff_restore


end module o_passive_tracer_mod
!=============================================================================
!  Performs a time step of the ocean model
!=============================================================================

subroutine ocean_step
  use o_param
  use o_array
  use o_mixing_kpp_mod
  use o_mixing_pp_mod
  use o_mixing_my2p5_mod
  use o_mixing_tidal_mod
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use o_mesh
  use o_elements
  use o_solver
  use g_config
  use g_parfe

  implicit none
  integer      :: i, m, row, row2, row3, n3
  real(kind=8) :: t0,t1,t2,t3,t4,t5,t6,t7, t8,t9,t10, t11

  t0=MPI_Wtime() 

  n3=ToDim_nod3d

  do row=1,n3                           
     row2=row+n3                   
     uf0(row)=uf(row)                 ! uf0 & uf: u^n
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                     
     uf0(row3)=uf(row3)
#endif
  end do

  !-----------------------------------------------------------

  call compute_density

  if(grid_type/=2) then
     call compute_pressure            ! compute hpressure
  end if
  if(grid_type/=1) then
     call compute_pressure_force      ! compute pressure gradient
  end if
  t2=MPI_Wtime()    

  !-----------------------------------------------------------

  if(trim(mix_scheme)=='KPP') then
     call oce_mixing_kpp(Av, Kv)
  elseif(trim(mix_scheme)=='PP') then
     call oce_mixing_pp
  elseif(trim(mix_scheme)=='MY2p5') then
     call oce_mixing_MY2p5
  end if
  if(tidal_mixing) call oce_mixing_tidal(Av, Kv)
  t3=MPI_Wtime()

  !------------------------------------------------------------

  call velocity_rhs                 ! rhs for u*
  if(use_vertvisc_impl) call uv_sfc_bott_bc
  t4=MPI_Wtime()    

  if(lump_uv_matrix) then           ! solver for du*
     call uv_solve
  else
     call solve(solve_u)   
     iteruv_first=.false.
     call solve(solve_v)
     call com_3d(duf(1:n3))              
     call com_3d(duf(1+n3:2*n3))           
  endif

#ifdef use_non_hydrostatic 
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))          
#endif 
  do row=1,n3                            
     row2=row+n3                   
     uf(row)=uf(row)+duf(row)      ! uf: u*
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                         
     uf(row3)=uf(row3)+duf(row3)
#endif
  end do

  if(use_vertvisc_impl) then      ! apply implicit vertical viscosity
     call impl_vertvisc
  end if
  t5=MPI_Wtime() 

  !--------------------------------------------------------------

  ssh0=ssh                        ! ssh & ssh0: ssh^n  

#ifdef use_opbnd_tide
  call update_tidal_opbnd         ! update tidal open boundary
#endif

  call compute_ssh_rhs            ! ssh rhs

  call solve(solve_ssh)           ! solve dssh
  ssh=ssh0+dssh                   ! ssh: ssh^n+1
  call com_2D(ssh)                       

#ifdef use_non_hydrostatic
  nhp0=nhp
  call compute_nhp_rhs            ! nhp rhs
  call solve(solve_nhp)           ! solve nhp
  call com_3D(nhp)                       
#endif   
  t6=MPI_Wtime()      

  do row=1,n3                              
     row2=row+n3                        
     uf0(row)=uf(row)             ! uf0: u*     
     uf0(row2)=uf(row2) 
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf0(row3)=uf(row3)
#endif	
  end do

  !---------------------------------------------------------------

  call velocity_rhs_update        ! Update rhs: contribution from ssh/nhp
  t7=MPI_Wtime()

  if(lump_uv_matrix) then         ! solve for full du 
     call uv_solve
  else
     call solve(solve_u)                          
     call solve(solve_v)
     call com_3d(duf(1:n3))             
     call com_3d(duf(1+n3:2*n3))        
  endif
#ifdef use_non_hydrostatic
  call solve(solve_w)
  call com_3d(duf(1+2*n3:3*n3))         
#endif
  do row=1,n3                            
     row2=row+n3                         
     uf(row)=uf(row)+duf(row)     ! uf: u^n+1  
     uf(row2)=uf(row2)+duf(row2)
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uf(row3)=uf(row3)+duf(row3)
#endif	
  end do
  t8=MPI_Wtime()

  !---------------------------------------------------------

#ifndef use_non_hydrostatic
  call compute_vvel_rhs           ! vertical rhs
  call solve_wpot                 ! solve for w potential 
#endif
  t9=MPI_Wtime()  

  !----------------------------------------------------------

#ifdef use_fullfreesurf
  call update_mesh
  call update_matrices
#endif 
  t10=MPI_Wtime() 

  !-----------------------------------------------------------

#ifdef use_cavity
  call cavity_heat_water_fluxes
#endif

  !-----------------------------------------------------------

  if(Redi_GM) call compute_neutral_slope   ! calc. neutral slope

  !-----------------------------------------------------------

#ifdef use_tracer_gls
  call tsstiff_fill_gls          ! tracer matrix/rhs
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  do i=1,num_tracer
     call solve(solve_tra+i-1)   ! solve for tracer
     call com_3D(dtracer(:,i))          
  end do
  do i=1,num_tracer                     
     tracer(:,i)=tracer(:,i)+dtracer(:,i)
  end do

#else

#ifdef use_tracer_fct
  call tracer_rhs_tg
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  call fct_tracer_solve 
#else
  call tracer_rhs_tg
  call ts_sfc_bc
  if(use_passive_tracer) call ptr_sfc_bc
  if(use_age_tracer) call age_tracer_tendency
  t11=MPI_Wtime()
  if(lump_ts_matrix) then
     call tracer_solve
  else
     do i=1,num_tracer
        call solve(solve_tra+i-1)
	call com_3D(dtracer(:,i))      
     end do
  endif
  do i=1,num_tracer
     do row=1,n3                         
        tracer(row,i)=tracer(row,i)+dtracer(row,i)
     end do
  end do
#endif

  if(use_vertdiff_impl) then
     call impl_vertdiff        ! apply implicit vertical diff.
  end if
#endif

  if(use_passive_tracer) call ptr_cutoff_restore
  if(use_age_tracer) call age_tracer_cutoff_restore

  !--------------------------------------------------------------

  t1=MPI_Wtime()
  iter_first = .false.
  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then      
     write(*,*)
     write(*,*) 'ocean  took   ', t1-t0
     write(*,*) 'Dens/pressure ', t2-t0
     write(*,*) 'mixing scheme ', t3-t2
     write(*,*) 'v rhs         ', t4-t3
     write(*,*) 'solve v_star  ', t5-t4
     write(*,*) 'pressure_solve', t6-t5
     write(*,*) 'rhs_update    ', t7-t6
     write(*,*) 'solve full v  ', t8-t7
#ifndef use_non_hydrostatic
     write(*,*) 'vert_vel_solve', t9-t8
#endif
#ifdef use_fullfreesurf
     write(*,*) 'update mesh   ', t10-t9
#endif
     write(*,*) 'tra_assemble  ', t11-t10
     write(*,*) 'tra_solve     ', t1-t11
     write(*,*)
  endif

end subroutine ocean_step
!=================================================================
subroutine ocean_init
  !reads the initial state or the restart file for the ocean
  use o_param
  use o_array
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_clock
  use g_parfe
  implicit none
  
  integer       :: i, node

  ! ocean dynamic fields and active tracers

  ! read ocean status
  if (.not.r_restart) then
     uf=0.0
     ssh=0.0
     if(mype==0) write(*,*) 'read ocean T&S climate data ', trim(OceClimaDataName)
    call ocean_init_fets1
if(mype==0) write(*,*) 'ocean_init_fets: finished '
     
#ifdef use_cavity
     call init_cavity_ts
#endif
  else
     if(mype==0) write(*,*) 'read ocean restart file'
     call oce_input
     if(mix_scheme=='MY2p5') call read_MY_vara ! this shoulb be separated later!!
  end if
  ! init surface values,
  ! will be updated during iteration according to surface restoring spec.
! EDIT 
     if(mype==0) write(*,*) 'control T&S climate data ', trim(OceClimaDataName)
 
  do i=1, ToDim_nod2D     
     node=nod3D_below_nod2D(1,i)       
     Tsurf(i)=tracer(node,1)          
     Ssurf(i)=tracer(node,2)          
  end do
 
    
  ! ocean passive tracers
 
  if (use_passive_tracer) then
     if(mype==0) write(*,*) 'initialize passive tracers'
     call passive_tracer_init
     if(ptr_start_year<yearnew .or. &
          (ptr_start_year==yearnew .and. (daynew>1 .or. timenew>dt))) then
        if(mype==0) write(*,*) 'read passive tracer restart fields'
        call passive_tracer_input
     end if
  end if    
  
    
  ! ocean age tracers
 
  if (use_age_tracer) then
     if(mype==0) write(*,*) 'initialize age tracers'
     call age_tracer_init
     if(age_tracer_start_year<yearnew .or. &
          (age_tracer_start_year==yearnew .and. (daynew>1 .or. timenew>dt))) then
        call age_tracer_input
     end if
  end if


  ! backup all tracers

  tracer0=tracer

end subroutine ocean_init
!
!----------------------------------------------------------------------------
!
subroutine ocean_init_back
  !A backup for ocean model test cases
  !reads the initial state or the restart file for the ocean
  use o_param
  use o_array
  use g_parfe
  use o_mesh
  use g_config
  use g_clock
  implicit none
  !
  integer                     :: i, eof_flag, num
  real(kind=8)                :: x, y, z, h
  real(kind=8)                :: dtdz, Tb
  real(kind=8)                :: Bu, N, rs

  tracer=0.

  dtdz=0.
  tracer(:,1)=10.0
  dtdz=2.0/rho0/2.e-4/5000.
  tracer(:,1)=10.0+dtdz*coord_nod3D(3,:)

  do i=1,myDim_nod2d+eDim_nod2d         
     y=coord_nod2d(2,i)
     x=coord_nod2d(1,i)
     rs=sqrt((y)**2+(x+250.e3/r_earth)**2)
     if(rs<=25.e3/r_earth) then

        num=num_layers_below_nod2d(i)

        ! !!!!!!!!!!!!!!!!!!!!!! deleted
     end if
  end do

  !Bu=1.0
  !N=Bu*1.0e-4*25.0e3/4500.0
  !!vertical gradient of temperature
  !dtdz=(N**2)/g/2.e-4
  !!background T
  !tracer(:,1)=dtdz*coord_nod3D(3,:)

  do i=1, myDim_nod2D+eDim_nod2D    
     num=nod3D_below_nod2D(1,i)     
     Tsurf(i)=tracer(num,1)         
     Ssurf(i)=tracer(num,2)        
  end do
  tracer0=tracer

  if(mype==0) write(*,*) 'ocean initialization done'
end subroutine ocean_init_back
!
!----------------------------------------------------------------------------
!
subroutine ocean_init_fets1

  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_PARFE
  use g_clock
  use g_read_other_NetCDF

  implicit none
  !
  integer                       :: i, j, num, n2, nlay
  integer                       :: num_reg
  real(kind=8)                  :: x, y

  real(kind=8), dimension(47)   :: sa_ae, te_ae, x_reg_ae
  real(kind=8), dimension(1212) :: sa_tss, te_tss, x_reg_tss
  real(kind=8), dimension(452)  :: sa_bl, te_bl, x_reg_bl

  real(kind=8), dimension(5)    :: auxx_tss, auxy_tss
  real(kind=8), dimension(5)    :: auxx_bl, auxy_bl

  integer, allocatable          :: row(:)
  character*80                  :: filename
  logical                       :: flag, inside
 

  !==========================
  ! open T/S station data files
  ! data is from a cruise ( fall )  
  !==========================

  if (mype==0) write(*,*) 'reading aegean1.txt ...'
  
  filename=trim(ClimateDataPath)//'aegean1.txt'
  open(19, file=filename, status='old')

  num_reg=47

  ! read raw data and do interpolation
  do i=1,num_reg
    read(19,*) x_reg_ae(i), te_ae(i), sa_ae(i)
  end do
  close(19)

  do n2=1,myDim_nod2d+eDim_nod2D
     nlay=num_layers_below_nod2d(n2)+1
     allocate(row(nlay))
     row=nod3d_below_nod2d(1:nlay,n2)

     do i=1,num_reg
       do j=1,nlay
         if (ABS(x_reg_ae(i)-coord_nod3d(3,row(j)))< 0.01) then
          tracer(row(j),1)=te_ae(i)
          tracer(row(j),2)=sa_ae(i)
         end if
       end do
     end do

      do j=1,nlay
        if (x_reg_ae(num_reg) .lt. coord_nod3d(3,row(j))) then 
         tracer(row(j),1)=te_ae(num_reg)
         tracer(row(j),2)=sa_ae(num_reg)
        end if

        if (x_reg_ae(1) .gt. coord_nod3d(3,row(j))) then 
         tracer(row(j),1)=te_ae(1)
         tracer(row(j),2)=sa_ae(1)
        end if
      end do 

    deallocate(row)
  end do


  if (mype==0) write(*,*) 'reading tss_ts.txt ...'

  filename=trim(ClimateDataPath)//'tss_ts.txt'
  open(19,file=filename, status='old')

  num_reg=1212

  ! correct small area in aegean sea

  auxx_tss(1) = 26.50 ;  auxy_tss(1) = 40.35
  auxx_tss(2) = 26.70 ;  auxy_tss(2) = 40.48
  auxx_tss(3) = 26.90 ;  auxy_tss(3) = 40.60
  auxx_tss(4) = 26.90 ;  auxy_tss(4) = 40.80
  auxx_tss(5) = 26.50 ;  auxy_tss(5) = 40.80

  ! read raw data and do interpolation
  do i=1,num_reg
     read(19, *) x_reg_tss(i),te_tss(i),sa_tss(i)
  end do
  close(19)

  do n2=1,myDim_nod2d+eDim_nod2D
     x=coord_nod2d(1,n2)/rad             ! latitude  (radian)
     y=coord_nod2d(2,n2)/rad             ! longitude (radian)
     flag=inside(x,y,auxx_tss,auxy_tss,5)
     if(x>26.6 .AND. y>40.0) then
  !   if (flag==.true.) cycle
     if (flag) cycle
     nlay=num_layers_below_nod2d(n2)+1
     allocate(row(nlay))
     row=nod3d_below_nod2d(1:nlay,n2)

     do i=1,num_reg
       do j=1,nlay
         if (ABS(x_reg_tss(i)-coord_nod3d(3,row(j)))< 0.001) then
          tracer(row(j),1)=te_tss(i)
          tracer(row(j),2)=sa_tss(i)
         end if
       end do
     end do

      do j=1,nlay
        if (x_reg_tss(num_reg) .lt. coord_nod3d(3,row(j))) then 
         tracer(row(j),1)=te_tss(num_reg)
         tracer(row(j),2)=sa_tss(num_reg)
        end if

        if (x_reg_tss(1) .gt. coord_nod3d(3,row(j))) then 
         tracer(row(j),1)=te_tss(1)
         tracer(row(j),2)=sa_tss(1)
        end if
      end do

    deallocate(row)
   end if 
!   end if 
  end do

  if (mype==0) write(*,*) 'reading black1.txt ...'

  filename=trim(ClimateDataPath)//'black1.txt'
  open(19,file=filename, status='old')

  num_reg=452

  ! read raw data and do interpolation
  do i=1,num_reg
     read(19, *) x_reg_bl(i),te_bl(i),sa_bl(i)
  end do
  close(19)

  ! correct small area in the black sea

  auxx_bl(1) = 29.50 ;  auxy_bl(1) = 41.00
  auxx_bl(2) = 30.50 ;  auxy_bl(2) = 41.00
  auxx_bl(3) = 32.00 ;  auxy_bl(3) = 40.90
  auxx_bl(4) = 32.00 ;  auxy_bl(4) = 41.40
  auxx_bl(5) = 29.50 ;  auxy_bl(5) = 41.40

  do n2=1,myDim_nod2d+eDim_nod2D
   x=coord_nod2d(1,n2)/rad             ! latitude  (radian)
   y=coord_nod2d(2,n2)/rad             ! longitude (radian)
   flag=inside(x,y,auxx_bl,auxy_bl,5)

!   if (flag==.true.) then 
   if (flag) then 
!     write(*,*) n2, x, y
     nlay=num_layers_below_nod2d(n2)+1
     allocate(row(nlay))

     row=nod3d_below_nod2d(1:nlay,n2)

     do i=1,num_reg
      do j=1,nlay
       if (ABS(x_reg_bl(i)-coord_nod3d(3,row(j)))< 0.0001) then
        tracer(row(j),1)=te_bl(i)
        tracer(row(j),2)=sa_bl(i)
       end if
      end do
     end do

     do j=1,nlay
      if (x_reg_bl(num_reg) .lt. coord_nod3d(3,row(j))) then
       tracer(row(j),1)=te_bl(num_reg)
       tracer(row(j),2)=sa_bl(num_reg)
      end if

      if (x_reg_bl(1) .gt. coord_nod3d(3,row(j))) then
       tracer(row(j),1)=te_bl(1)
       tracer(row(j),2)=sa_bl(1)
      end if
     end do

     deallocate(row)

    end if 

    if(y>41.15) then
     nlay=num_layers_below_nod2d(n2)+1
     allocate(row(nlay))

     row=nod3d_below_nod2d(1:nlay,n2)

     do i=1,num_reg
      do j=1,nlay
       if (ABS(x_reg_bl(i)-coord_nod3d(3,row(j)))< 0.0001) then
        tracer(row(j),1)=te_bl(i)
        tracer(row(j),2)=sa_bl(i)
       end if
      end do
     end do

     do j=1,nlay
      if (x_reg_bl(num_reg) .lt. coord_nod3d(3,row(j))) then 
       tracer(row(j),1)=te_bl(num_reg)
       tracer(row(j),2)=sa_bl(num_reg)
      end if

      if (x_reg_bl(1) .gt. coord_nod3d(3,row(j))) then 
       tracer(row(j),1)=te_bl(1)
       tracer(row(j),2)=sa_bl(1)
      end if
     end do
  
     deallocate(row)
  end if
 end do

  do i=1, myDim_nod2D+eDim_nod2D
     num=nod3D_below_nod2D(1,i)
     Tsurf(i)=tracer(num,1)
     Ssurf(i)=tracer(num,2)
  end do
  tracer0=tracer

end subroutine ocean_init_fets1
!
!----------------------------------------------------------------------------
!
subroutine ocean_array_setup
  ! Sets up the arrays needed by the ocean part
  use o_param
  use o_array
  use o_mixing_kpp_mod
  use o_mixing_MY2p5_mod
  use o_mixing_tidal_mod
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  integer         :: size3D, size2D         

  size3d=ToDim_nod3D
  size2d=ToDim_nod2D

  ! density and pressure arrays
  allocate(density_insitu(size3D), density_ref(size3D)) 
  density_insitu=0.
  density_ref=0.
  if(grid_type/=2) then
     allocate(hpressure(size3D))
     hpressure=0.
  end if
  if(grid_type/=1) then
     allocate(PGF(2,max_num_layers-1,myDim_elem2D))
     PGF=0.
     call init_pressure_force
  end if

  ! Arrays used for ocean boundary forcing
  allocate(stress_x(size2d))
  allocate(stress_y(size2d)) 
  allocate(heat_flux(size2d)) 
  allocate(water_flux(size2d))
  allocate(Tsurf(size2d))
  allocate(Ssurf(size2d)) 
  allocate(ts_sfc_force(size2d, 2))
  allocate(uv_sfc_force(size2d, 2))
  allocate(uv_bott_force(size2d, 2))
  stress_x=0.
  stress_y=0.
  heat_flux=0.
  water_flux=0.
  ts_sfc_force=0.0
  uv_sfc_force=0.0
  uv_bott_force=0.0

  ! T, S fields, their increments and rhs
  allocate(tracer(size3d,num_tracer))     
  allocate(dtracer(size3d,num_tracer))   
  allocate(tracer_rhs(size3d,num_tracer)) 
  allocate(tracer0(size3d,num_tracer))
  tracer=0.0
  tracer0=0.0
  dtracer=0.0
  tracer_rhs=0.0

  ! u, v fields and their rhs 
#ifndef use_non_hydrostatic
  allocate(uf(2*size3D))             
  allocate(uf0(2*size3D))
  allocate(duf(2*size3D))           
  allocate(uv_rhs(2*size3D))
#else
  allocate(nhp(size3D), nhp0(size3D), nhp_rhs(size3D))   
  allocate(uf(3*size3D))          
  allocate(uf0(3*size3D))
  allocate(duf(3*size3D))           
  allocate(uv_rhs(3*size3D)) 
  nhp=0.
  nhp0=0.
#endif
  uf=0.
  uf0=0.
  duf=0.
  uv_rhs=0.

  !ssh
  allocate(ssh(size2D), ssh0(size2D), dssh(size2d), ssh_rhs(size2D)) 
  ssh=0.
  ssh0=0.
  dssh=0.
  uf=0.       ! ::OG::  not needed
  uf0=0.      ! ::OG::  not needed

  ! arrays for the AB2 coriolis case
  if(.not.use_cori_semi) then
     allocate(ucori(size3d), vcori(size3d))
     allocate(ucori_back(size3D), vcori_back(size3D))
     ucori=0.
     vcori=0.
     ucori_back=0.
     vcori_back=0.
  endif

  ! rhs of w equation and w-potential field
#ifndef use_non_hydrostatic
  allocate(wrhs(size3D),w(size3D))   
  wrhs=0.
  w=0.
#endif

  ! arrays for salt fluxes
  allocate(virtual_salt(size2d), relax_salt(size2d))
  virtual_salt=0.
  relax_salt=0.
#ifdef use_fullfreesurf
  allocate(real_salt_flux(size2d))
  real_salt_flux=0.
#endif  

  ! Redi/GM
  if (Redi_GM) then    
     if(nslope_version==1) then 
        allocate(neutral_slope(3,max_num_layers-1,myDim_elem2d))
        neutral_slope=0.0
     else
        allocate(neutral_slope_elem(3,myDim_elem3d))
        neutral_slope_elem=0.0
     end if
  end if

  ! vertical mixing 
  allocate(Av(ToDim_nod3d))
  Av=0.0
  if(trim(mix_scheme)=='KPP') then
     allocate(Kv(ToDim_nod3d,2))
     Kv=0.0
     call oce_mixing_kpp_init
  else
     allocate(Kv(ToDim_nod3d,1))
     Kv=0.0
  end if
  if(trim(mix_scheme)=='MY2p5') then
     call oce_mixing_MY2p5_init
  end if
  if(tidal_mixing) call oce_mixing_tidal_init 
 
  ! initialize the fct scheme
#ifdef use_tracer_fct    
  call fct_init   
#endif

  if(mype==0) write(*,*) 'Ocean arrays have been set up'
end subroutine ocean_array_setup
!
!----------------------------------------------------------------------------
!
subroutine set_coriolis_param
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer         :: i, j
  real(kind=8)    :: midlat, lon, lat, rlon, rlat
  
  ! nodal

  allocate(coriolis_param_nod2D(ToDim_nod2d))  
  coriolis_param_nod2D=0.0
  if (.not.rotated_grid) then
     coriolis_param_nod2D=2.0*omega*sin(coord_nod2D(2,:))
  else
     do i=1,ToDim_nod2d                        
        rlat=coord_nod2D(2,i)
        rlon=coord_nod2D(1,i)
        call r2g(lon, lat, rlon, rlat)
        coriolis_param_nod2D(i)=2.0*Omega*sin(lat)
     end do
  end if
  if(fplane) then
     coriolis_param_nod2D=f_fplane 
  end if
  if(betaplane) then    
     midlat=(maxval(coord_nod2d(2,:))+minval(coord_nod2d(2,:)))/2.0 
     coriolis_param_nod2D=f_fplane+beta_betaplane*(coord_nod2D(2,:)-midlat)*r_earth
  end if

  ! elementwise 

  allocate(coriolis_param_elem2D(myDim_elem2D))    
  coriolis_param_elem2D=0.0
  do i=1,myDim_elem2D                                            
     coriolis_param_elem2d(i)=sum(coriolis_param_nod2D(elem2D_nodes(:,i)))/3.0
  end do
end subroutine set_coriolis_param
!
!----------------------------------------------------------------------------
!
subroutine prepare_init_data
  ! this routine is only kept here for a backup, not used by the model
  ! read nc data and save to formatted data 
  use o_PARAM
  use g_config
  use g_parfe
  implicit none

#include "netcdf.inc" 

  integer			:: k, i, j
  integer,parameter             :: num_z=33
  integer			:: itime, latlen, lonlen
  integer			:: status, ncid, varid
  integer			:: lonid, latid
  integer			:: istart(3), icount(3)
  real(kind=8), allocatable	:: lon(:), lat(:)
  real(kind=8), allocatable	:: nc_temp(:,:,:), nc_salt(:,:,:)
  real(kind=8)                  :: dep(num_z)
  character                     :: cdep1*2, cdep2*3, cdep3*4
  character(15)			:: vari
  character(80)                	:: file  

  file=trim(ClimateDataPath)//'Winter_phc2.1_beta.dat.nc'

  ! open file
  status=nf_open(file, nf_nowrite, ncid)
  if (status.ne.nf_noerr)then
     print*,'ERROR: CANNOT READ init_data FILE CORRECTLY !!!!!'
     print*,'Error in opening netcdf file'//file
     call par_ex
     stop
  endif

  ! lat
  status=nf_inq_dimid(ncid, 'lat', latid)
  status=nf_inq_dimlen(ncid, latid, latlen)
  allocate(lat(latlen))
  status=nf_inq_varid(ncid, 'lat', varid)
  status=nf_get_vara_double(ncid,varid,1,latlen,lat)

  ! lon
  status=nf_inq_dimid(ncid, 'lon', lonid)
  status=nf_inq_dimlen(ncid, lonid, lonlen)
  allocate(lon(lonlen))
  status=nf_inq_varid(ncid, 'lon', varid)
  status=nf_get_vara_double(ncid,varid,1,lonlen,lon)

  ! depth
  dep=(/ 0., 10., 20., 30., 50., 75., 100., 125., 150., 200., &
       250., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., &
       1200., 1300., 1400., 1500., 1750., 2000., 2500., 3000., 3500., &
       4000., 4500., 5000., 5500. /)

  ! data
  allocate(nc_temp(lonlen,latlen,33))
  allocate(nc_salt(lonlen,latlen,33))

  do k=1,33
     if(k==1) then
        vari='MAM_Temp_00m'
     elseif(dep(k)<100.) then
        write(cdep1,'(i2)') int(dep((k)))
        vari='MAM_Temp_'//cdep1//'m'
     elseif(dep(k)<1000.) then
        write(cdep2,'(i3)') int(dep((k)))
        vari='MAM_Temp_'//cdep2//'m'
     else
        write(cdep3,'(i4)') int(dep((k)))
        vari='MAM_Temp_'//cdep3//'m'
     endif

     status=nf_inq_varid(ncid, vari, varid)
     if (status.ne.nf_noerr)then
        write(*,*) 'error by getting varid for temp'
        call abort
     end if
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,nc_temp(:,:,k))
     if (status.ne.nf_noerr)then
        write(*,*) 'error when reading temp'
        call abort
     end if

     if(k==1) then
        vari='MAM_Salt_00m'
     elseif(dep(k)<100.) then
        write(cdep1,'(i2)') int(dep((k)))
        vari='MAM_Salt_'//cdep1//'m'
     elseif(dep(k)<1000.) then
        write(cdep2,'(i3)') int(dep((k)))
        vari='MAM_Salt_'//cdep2//'m'
     else
        write(cdep3,'(i4)') int(dep((k)))
        vari='MAM_Salt_'//cdep3//'m'
     endif

     status=nf_inq_varid(ncid, vari, varid)
     if (status.ne.nf_noerr)then
        write(*,*) 'error by getting varid for salt'
        call abort
     end if
     istart = (/1,1,itime/)
     icount= (/lonlen,latlen,1/)
     status=nf_get_vara_double(ncid,varid,istart,icount,nc_salt(:,:,k))
     if (status.ne.nf_noerr)then
        write(*,*) 'error when reading salt'
        call abort
     end if
  end do

  ! close file
  status=nf_close(ncid)

  ! change dep to negative values
  dep=-dep

!!$  ! save to formated output
!!$  open(36,file='Winter_phc2.1_beta.dat')
!!$!  write(36,*) lonlen, latlen, num_z 
!!$!  write(36,'(1f8.2)') lon
!!$!  write(36,'(1f8.2)') lat
!!$!  write(36,'(1f8.2)') dep
!!$  do i=1, lonlen
!!$     do j=1, latlen
!!$        write(36,'(1f10.5)') nc_temp(i,j,1:num_z)         
!!$     end do
!!$  end do
!!$  do i=1, lonlen
!!$     do j=1, latlen
!!$        write(36,'(1f10.5)') nc_salt(i,j,1:num_z)         
!!$     end do
!!$  end do
!!$  close(36) 


  open (1,file='PHC_t.dat')
  open (2,file='PHC_s.dat')
  write(1, *) -dep
  write(2, *) -dep
  do i=1, lonlen
     do j=1, latlen              
        write(1, *) lon(i), lat(j), (nc_temp(i,j, k), k=1, num_z)
        write(2, *) lon(i), lat(j), (nc_salt(i,j, k), k=1, num_z)
     end do
  end do
  close(1)
  close(2)  


  deallocate(lon, lat, nc_temp, nc_salt)
end subroutine prepare_init_data
!
!-----------------------------------------------------------------------------------
!
!
      FUNCTION inside (Xo,Yo,Xb,Yb,Nb)
!
!svn $Id: inside.F 435 2010-01-02 15:55:10Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2010 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!      implicit none
      integer Nstep
      parameter (Nstep=128)
!
      logical inside
      integer Nb, crossings, i, inc, k, kk, nc
      integer indexx(Nstep)
      real*8 Xb(Nb+1), Yb(Nb+1), Xo, Yo, dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      do kk=0,Nb-1,Nstep
        nc=0
        do k=kk+1,MIN(kk+Nstep,Nb)
          if (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0).and.                     &
     &        (Xb(k).ne.Xb(k+1))) then
            nc=nc+1
            indexx(nc)=k
          end if
        end do
        do i=1,nc
          k=indexx(i)
          if (Xb(k).ne.Xb(k+1)) then
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            if ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) then
              crossings=1
              go to 10
            else if (((dx1.eq.0.0).and.(Yo.ge.Yb(k  ))).or.             &
     &              ((dx2.eq.0.0).and.(Yo.ge.Yb(k+1)))) then
              inc=1
            else if ((dx1*dx2.gt.0.0).and.                              &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0)) then
              inc=2
            end if
            if (Xb(k+1).gt.Xb(k)) then
              crossings=crossings+inc
            else
              crossings=crossings-inc
            end if
          end if
        end do
      end do
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  if (crossings.eq.0) then
        inside=.false.
      else
        inside=.true.
      end if 
      return
      end
!

!example routines for initializing and/or restoring the
!tracer fields in process studies (e.g., overflows)
!--------------------------------------------------------------------


subroutine init_ts_CavH1
  ! init T/S source for cavity test case of Hunter 1NN
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  tracer(:,1)=-1.9
  tracer(:,2)=34.4
  
  tracer0=tracer  
end subroutine init_ts_CavH1
!
!--------------------------------------------------------------------
!
subroutine init_tracers_FOtide
  ! init T/S source for FOtide
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  integer                     :: row, n, ind, fileID
  integer, allocatable        :: index_source(:), temp_arr2d(:)
  real(kind=8)                :: z
  real(kind=8)                :: zn(5), tn(5), sn(5)
  real(kind=8)                :: ts(4), ss(4), tint(4), sint(4)
  character*80                :: file_name

  data zn /0.0, -300.0, -550.0, -2000.0, -4000.0/
  data tn /-1.9, -1.75, 0.7, 0.0, -0.3/
  data sn /34.3, 34.42, 34.68, 34.66, 34.65/

  ts=(tn(1:4)-tn(2:5))/(zn(1:4)-zn(2:5))
  ss=(sn(1:4)-sn(2:5))/(zn(1:4)-zn(2:5))
  tint=tn(1:4)-ts*zn(1:4)
  sint=sn(1:4)-ss*zn(1:4)

  !background
  do row=1,myDim_nod3d+eDim_nod3d
     z=coord_nod3d(3,row)

     if(z>zn(2)) then
        tracer(row,1)=ts(1)*z+tint(1)
        tracer(row,2)=ss(1)*z+sint(1)        
     else if(z>zn(3)) then
        tracer(row,1)=ts(2)*z+tint(2)
        tracer(row,2)=ss(2)*z+sint(2)  
     else if(z>zn(4)) then
        tracer(row,1)=ts(3)*z+tint(3)
        tracer(row,2)=ss(3)*z+sint(3)    
     else if(z>zn(5)) then
        tracer(row,1)=ts(4)*z+tint(4)
        tracer(row,2)=ss(4)*z+sint(4)     
     else
        tracer(row,1)=tn(5)
        tracer(row,2)=sn(5)
     end if
     tracer(row,3)=0.0
  end do

  !source water (ISW)

  !read in info. specifying the initial source region 
  allocate(index_source(myDim_nod2D+eDim_nod2D))	 			

  file_name=trim(meshpath)//'index_source_region.out' 
  fileID=150
  open(fileID, file=file_name)

  allocate(temp_arr2d(nod2d))
  temp_arr2d=0
  do n=1, myDim_nod2D+eDim_nod2D
     temp_arr2d(myList_nod2D(n))=n
  end do

  do n=1,nod2D
     read(fileID,*) ind
     if (temp_arr2d(n)>0) then
        index_source(temp_arr2d(n))=ind
     end if
  end do
  close(fileID)

  ! set source values
  do n=1, myDim_nod3D+eDim_nod3D
     row=nod2d_corresp_to_nod3d(n)
     if(index_source(row)>0 .and. coord_nod3d(3,n)<=-400.0) then
        tracer(n,1)=-2.2
        tracer(n,2)=34.6
        tracer(n,3)=1.0
     end if
  end do

  tracer0=tracer  
 
  deallocate(index_source, temp_arr2d)

  if(mype==0) write(*,*) 'ambient and source water prescribed: FOtide'

end subroutine init_tracers_FOtide
!
!--------------------------------------------------------------------
!
subroutine restore_source_FOtide
  ! restore T/S source for FOtide setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
 
  integer                     :: row
  real(kind=8)                :: x, y

  ! restoring for tracers

  do row=1,myDim_nod3d+eDim_nod3d
     x=coord_nod3d(1,row)
     y=coord_nod3d(2,row)
 
     if(y<=-76.0 .and. x>-42.0 .and. x<-31.0) then   
        tracer(row,:)=tracer0(row,:)
     end if

  end do
  
end subroutine restore_source_FOtide
!
!=============================================================================
!
subroutine init_source_RS1
  ! init T/S source for RS1 setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none

  integer                     :: i, k, n, n2, num_hor_reg, num_ver_reg, num_sn
  integer, allocatable        :: source_nodes(:)
  real(kind=8)                :: x, y, z, rx, ry
  real(kind=8), allocatable   :: ver_reg(:), hor_reg(:), raw_t(:,:), raw_s(:,:)


  !tracer fields for restoring

  !passive tracer
  tracer(:,3)=0.0

  !source profile
  open(39,file=trim(MeshPath)//'RS_init_special.out', status='old')
  read(39,*) num_hor_reg, num_ver_reg
  allocate(hor_reg(num_hor_reg),ver_reg(num_ver_reg))
  read(39,*) hor_reg  
  read(39,*) ver_reg
  allocate(raw_t(num_hor_reg,num_ver_reg))
  allocate(raw_s(num_hor_reg,num_ver_reg))
  do i=1,num_hor_reg
     read(39, *) raw_t(i,1:num_ver_reg)
  end do
  do i=1,num_hor_reg
     read(39, *) raw_s(i,1:num_ver_reg)
  end do
  close(39)

  !2d nodes in the source region
  open(40,file=trim(MeshPath)//'source_2d_nodes.out', status='old')
  read(40,*) num_sn
  allocate(source_nodes(num_sn))
  read(40,*) source_nodes 
  close(40)

  !interpolate
  do n2=1,myDim_nod2d+eDim_nod2d
     if(any(source_nodes==myList_nod2d(n2))) then
	if(rotated_grid) then
           rx=coord_nod2d(1,n2)
           ry=coord_nod2d(2,n2)
           call r2g(x, y, rx, ry)
           x=x/rad   ! degree
           y=y/rad
        else
           x=coord_nod2d(1,n2)
           y=coord_nod2d(2,n2)
        end if
        do k=1,num_layers_below_nod2d(n2)+1
           n=nod3d_below_nod2d(k,n2)
           z=coord_nod3d(3,n)
           call interp_vertical_section(num_hor_reg, num_ver_reg, hor_reg, ver_reg, raw_t, &
                1, y, z, tracer(n,1))
           call interp_vertical_section(num_hor_reg, num_ver_reg, hor_reg, ver_reg, raw_s, &
                1, y, z, tracer(n,2))
           if(tracer(n,2)>=34.62) tracer(n,3)=1.0
        end do
     end if
  end do

  tracer0=tracer  !!!

  deallocate(ver_reg, hor_reg, raw_t, raw_s, source_nodes)

  if(mype==0) write(*,*) 'source water prescribed'

end subroutine init_source_RS1
!
!--------------------------------------------------------------------
!
subroutine restore_source_RS1
  ! restore T/S source for RS1 setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
 
  integer                     :: row
  real(kind=8)                :: x, y, rx, ry

  ! restoring for tracers

  do row=1,myDim_nod3d+eDim_nod3d
     
     if(rotated_grid) then
        rx=coord_nod3d(1,row)
        ry=coord_nod3d(2,row)
        call r2g(x, y, rx, ry)
        x=x/rad   ! degree
        y=y/rad
     else
        x=coord_nod3d(1,row)/rad
        y=coord_nod3d(2,row)/rad
     end if

     ! if((x>160. .and. x<170. .and. y<=-74.7) .or. (y<0.6897*x-193.3103 .and. y<-72.9)) then

     if(x>160. .and. x<171.2 .and. y<=-74.5) then   
        tracer(row,:)=tracer0(row,:)
     end if

  end do
  
end subroutine restore_source_RS1
subroutine fcn_density(t,s,z,rho)
  !
  ! - calculates insitu density as a function of potential temperature
  !   (t is relative to the surface)
  !   using the Jackett and McDougall equation of state (1992?)
  ! Qiang 02,07,2010:  Should this be updated (1995 or 2003)? The current 
  !   version is also different to the international equation of state 
  !   (Unesco 1983). What is the exact reference for this version then? 

  use o_PARAM
  implicit none

  real(kind=8), intent(IN)       :: t, s, z
  real(kind=8), intent(OUT)      :: rho                 
  real(kind=8)                   :: rhopot, bulk

  if(density_linear) then
     rho = -rho0 * 2.0e-4 * t
  else
     bulk = 19092.56 + t*(209.8925 				&
          - t*(3.041638 - t*(-1.852732e-3			&
          - t*(1.361629e-5))))				&
          + s*(104.4077 - t*(6.500517			&
          -  t*(.1553190 - t*(-2.326469e-4))))		&
          + sqrt(s**3)*(-5.587545				&
          + t*(0.7390729 - t*(1.909078e-2)))		&
          - z *(4.721788e-1 + t*(1.028859e-2		&
          + t*(-2.512549e-4 - t*(5.939910e-7))))		&
          - z*s*(-1.571896e-2				&
          - t*(2.598241e-4 + t*(-7.267926e-6)))		&
          - z*sqrt(s**3)					&
          *2.042967e-3 + z*z*(1.045941e-5			&
          - t*(5.782165e-10 - t*(1.296821e-7)))		&
          + z*z*s						&
          *(-2.595994e-7					&
          + t*(-1.248266e-9 + t*(-3.508914e-9)))

     rhopot = ( 999.842594					&
          + t*( 6.793952e-2			&
          + t*(-9.095290e-3			&
          + t*( 1.001685e-4			&
          + t*(-1.120083e-6			&
          + t*( 6.536332e-9)))))			&
          + s*( 0.824493				&
          + t *(-4.08990e-3			&
          + t *( 7.64380e-5			&
          + t *(-8.24670e-7			&
          + t *( 5.38750e-9)))))			&
          + sqrt(s**3)*(-5.72466e-3		&
          + t*( 1.02270e-4			&
          + t*(-1.65460e-6)))			&
          + 4.8314e-4*s**2)
     rho = rhopot / (1.0 + 0.1*z/bulk)
  end if
end subroutine fcn_density
!
!----------------------------------------------------------------------------
!
subroutine fcn_dens0(t,s,rho)
  ! - calculate density at ocean surface.
  ! - if input t is potential temperature, rho is then
  !   potential density (both referenced to surface).
  ! - This routine is used when diagnosing mixed layer thickness.
  !   Using this instead of fcn_density(t,s,0,rho) reduces computation load.

  implicit none
  real(kind=8), intent(IN)       :: t, s
  real(kind=8), intent(OUT)      :: rho                 
  real(kind=8)                    :: tn

  tn=t*1.00024

  rho = ( 999.842594			        & 
       + tn*( 6.793952e-2			&
       + tn*(-9.095290e-3			&
       + tn*( 1.001685e-4			&
       + tn*(-1.120083e-6			&
       + tn*( 6.536332e-9)))))			&
       + s*( 0.824493				&
       + tn *(-4.08990e-3			&
       + tn *( 7.64380e-5			&
       + tn *(-8.24670e-7			&
       + tn *( 5.38750e-9)))))			&
       + sqrt(s**3)*(-5.72466e-3		&
       + tn*( 1.02270e-4			&
       + tn*(-1.65460e-6)))			&
       + 4.8314e-4*s**2)

end subroutine fcn_dens0
!
!----------------------------------------------------------------------------
!
subroutine compute_ref_density
  use o_MESH
  use o_PARAM
  use o_array
  use g_PARFE
  implicit none
  !
  integer         :: n2, n, k
  real(kind=8)    :: T, S, z
  
  ! this should be modified later to get level mean t and s
  
  density_ref=0.0
  S=34.  
  T=4.0   

  do n2=1,myDim_nod2d+eDim_nod2d
     do k=1,num_layers_below_nod2d(n2)+1
        n=nod3d_below_nod2d(k,n2)
        z=min(coord_nod3D(3,n),0.)
        call fcn_density(T, S, z, density_ref(n))
     end do
  end do
  if(mype==0) write(*,*) 'Reference density computed.'
end subroutine compute_ref_density
!
!----------------------------------------------------------------------------
!
subroutine compute_density
  use o_MESH
  use o_PARAM
  use o_array
  use g_PARFE
  implicit none
  !
  integer         :: n2, n, k
  real(kind=8)    :: z
  !
  do n2=1,myDim_nod2d+eDim_nod2d
     do k=1,num_layers_below_nod2d(n2)+1
        n=nod3d_below_nod2d(k,n2)
        z=min(coord_nod3D(3,n), 0.0) 
        call fcn_density(tracer(n,1), tracer(n,2), z, density_insitu(n))
     end do
  end do
end subroutine compute_density
!
!----------------------------------------------------------------------------
!
subroutine compute_pressure
  ! Qiang, 16.12.2010: add the nonlinear free surface option
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use o_MATRICES
  use g_PARFE
#ifdef use_fullfreesurf
  use i_therm_parms
  use i_array
#endif
  implicit none
  !
  integer         :: i, n, node_lo, node_hi
  real(kind=8)    :: dens, denscor, z_up, z_lo, wd_ice_eff
  
  do n=1,myDim_nod2d+eDim_nod2d    
     node_hi = nod3D_below_nod2D(1,n)
     z_up=0.0 ! already re-checked, qiang, 15.03.2011
#ifdef use_fullfreesurf
     wd_ice_eff=(rhoice*m_ice(n)+rhosno*m_snow(n))*rho0r
     hpressure(node_hi)=g*min(wd_ice_eff,max_ice_loading)
#else
     hpressure(node_hi)=0.0
#endif
     do i=2, num_layers_below_nod2D(n)+1
        node_lo = nod3D_below_nod2D(i,n) 
        dens=0.5_8*(density_insitu(node_hi)-density_ref(node_hi) &
             +density_insitu(node_lo)-density_ref(node_lo))
        z_lo=coord_nod3d(3,node_lo)
        hpressure(node_lo)=hpressure(node_hi)+g*(z_up-z_lo)*dens*rho0r
        node_hi=node_lo
        z_up=z_lo
     end do
  end do

end subroutine compute_pressure
!
!----------------------------------------------------------------------------
!
subroutine init_pressure_force
  ! Initialization for PGF calculation on sigma grids
  ! -- assume that the difference in the number of layers between neighbour nodes
  ! is not more than one on sigma grids
  ! -- in case under cavities, assume that the cavity draft does not go deeper than
  ! the highest bottom depth under each 2d element 
  ! breaking these assumptions the model can not work!
  ! Qiang, 9,4,2010
  use o_mesh
  use o_elements
  use o_array
  use g_parfe
  implicit none

  integer       :: k, j, kk, elem, n2, n3, mloc(1), lay_el(3), lay_m
  integer       :: elnodes2(3), bottnodes(3), surfnodes(3), elnodes(6)
  real(kind=8)  :: zmean, bottdep(3), surfdep(3)

  allocate(dens_interp_nodes(4, max_num_layers-1, myDim_elem2D))

  !find nodes for the (4 points) cubic interpolation 
  do elem=1,myDim_elem2D                   
     if(grid_type_elem2d(elem)==0) cycle !not sigma grid part
     elnodes2=elem2d_nodes(:,elem)
     lay_el=num_layers_below_nod2d(elnodes2)+1
     lay_m=minval(lay_el)
     if(any(lay_el>lay_m)) lay_m=lay_m+1
     lay_el(1)=min(lay_el(1),lay_m)
     lay_el(2)=min(lay_el(2),lay_m)
     lay_el(3)=min(lay_el(3),lay_m)
     bottnodes(1)=nod3d_below_nod2d(lay_el(1),elnodes2(1))
     bottnodes(2)=nod3d_below_nod2d(lay_el(2),elnodes2(2))
     bottnodes(3)=nod3d_below_nod2d(lay_el(3),elnodes2(3))
     bottdep=coord_nod3d(3,bottnodes)
     surfnodes=nod3d_below_nod2d(1,elnodes2)
     surfdep=coord_nod3d(3,surfnodes)
     do k=1,lay_m-1
        elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
        if(k==lay_m-1) then
           elnodes(4:6)=bottnodes
        else
           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
        end if

        zmean=sum(coord_nod3d(3,elnodes))/6.0
        if(any(surfdep<zmean)) then
           mloc=minloc(surfdep)
           dens_interp_nodes(4, k, elem)=mloc(1)+3
           zmean=minval(surfdep)
        elseif(any(bottdep>zmean)) then
           mloc=maxloc(bottdep)
           dens_interp_nodes(4, k, elem)=mloc(1)    
           zmean=maxval(bottdep)
        else
           dens_interp_nodes(4, k, elem)=0        
        end if

        do j=1,3
           n2=elnodes2(j)
           do kk=1,lay_el(j)
              n3=nod3d_below_nod2d(kk,n2)
              if(coord_nod3d(3,n3)<=zmean) then
                 dens_interp_nodes(j, k, elem)=kk-1 
                 if(kk==1) dens_interp_nodes(j,k,elem)=1
                 exit
              end if
           end do ! kk   
        end do ! j
     end do ! k
  end do ! elem 

end subroutine init_pressure_force
!
!----------------------------------------------------------------------------
!
subroutine compute_pressure_force
  ! Calculating PGF on sigma grids
  ! -- assume that the difference in the number of layers between neighbour nodes
  ! is not more than one on sigma grids
  ! -- in case under cavities, assume that the cavity draft does not go deeper than
  ! the highest bottom depth under each 2d element 
  ! breaking these assumptions the model can not work!
  ! Qiang, 9,4,2010
  use o_param
  use o_mesh
  use o_elements
  use o_array
  use g_parfe
  implicit none

  integer            :: j, k, elem, bn, lay_el(3), lay_m
  integer            :: elnodes(6), elnodes2(3), bottnodes(3), surfnodes(3)
  integer            :: n2, ind, flag, s_nodes(4), s_ind(4)
  real(kind=8)       :: dx2d(3), dy2d(3), z, dz, bottdep(3), surfdep(3)
  real(kind=8)       :: p_grad(2), intz(2), intratio, intdens(2)
  real(kind=8)       :: rho(3), rho_grad(2)
  real(kind=8)       :: s_z(4), s_dens(4), s_H, aux1, aux2, aux(2)
  real(kind=8)       :: s_dup, s_dlo, a, b, c, d

  do elem=1, myDim_elem2D                  
     if(grid_type_elem2d(elem)==0) cycle !not sigma grid part
     elnodes2=elem2d_nodes(:,elem)
     dx2d=bafux_2d(:,elem)
     dy2d=bafuy_2d(:,elem)
     lay_el=num_layers_below_nod2d(elnodes2)+1
     lay_m=minval(lay_el)
     if(any(lay_el>lay_m)) lay_m=lay_m+1
     lay_el(1)=min(lay_el(1),lay_m)
     lay_el(2)=min(lay_el(2),lay_m)
     lay_el(3)=min(lay_el(3),lay_m)
     bottnodes(1)=nod3d_below_nod2d(lay_el(1),elnodes2(1))
     bottnodes(2)=nod3d_below_nod2d(lay_el(2),elnodes2(2))
     bottnodes(3)=nod3d_below_nod2d(lay_el(3),elnodes2(3))
     bottdep=coord_nod3d(3,bottnodes)
     surfnodes=nod3d_below_nod2d(1,elnodes2)
     surfdep=coord_nod3d(3,surfnodes)
     p_grad=0.0

     do k=1,lay_m-1
        elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
        if(k==lay_m-1) then
           elnodes(4:6)=bottnodes
        else
           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
        end if
        
        bn=dens_interp_nodes(4,k,elem)         
        if(bn==0) then
           z=sum(coord_nod3d(3,elnodes))/6.0
        elseif(bn<4) then
           z=bottdep(bn)
        else
           z=surfdep(bn-3)
        end if

        do j=1,3
           n2=elnodes2(j)
           ind=dens_interp_nodes(j,k,elem)    
           !cubic-spline interpolation
           s_ind=(/ind-1, ind, ind+1, ind+2/) 
           flag=0
           if(ind==1) then !assume the mesh has at lease 3 layers
              flag=-1
              s_ind(1)=1
           elseif(ind==lay_el(j)-1) then
              flag=1
              s_ind(4)=ind+1
           end if
           s_nodes=nod3d_below_nod2d(s_ind,n2)
           s_z=coord_nod3d(3,s_nodes)
           s_dens=density_insitu(s_nodes)-density_ref(s_nodes)
           s_H=s_z(3)-s_z(2)
           aux1=(s_dens(3)-s_dens(2))/s_H
           ! derivatives computed in a way to get monotonic profile
           if(flag==0) then
              aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
              s_dup=0.0
              if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
              aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
              s_dlo=0.0
              if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
           elseif(flag==1) then
              aux2=(s_dens(2)-s_dens(1))/(s_z(2)-s_z(1))
              s_dup=0.0
              if(aux1*aux2>0.)  s_dup=2.0*aux1*aux2/(aux1+aux2)
              s_dlo=1.5*aux1-0.5*s_dup
           else
              aux2=(s_dens(4)-s_dens(3))/(s_z(4)-s_z(3))
              s_dlo=0.0
              if(aux1*aux2>0.) s_dlo=2.0*aux1*aux2/(aux1+aux2)
              s_dup=1.5*aux1-0.5*s_dlo
           end if
           ! cubic polynomial coefficients
           a=s_dens(2)
           b=s_dup
           c=-(2.0*s_dup+s_dlo)/s_H + 3.0*(s_dens(3)-s_dens(2))/s_H**2
           d=(s_dup+s_dlo)/s_H**2 - 2.0*(s_dens(3)-s_dens(2))/s_H**3
           ! interploate
           dz=z-s_z(2)
           rho(j)=a+b*dz+c*dz**2+d*dz**3
        end do
        rho_grad(1)=sum(rho*dx2d) !elementwise density gradient
        rho_grad(2)=sum(rho*dy2d)

        !next, compute pressure gradient
        dz=(sum(coord_nod3d(3,elnodes(1:3)))-sum(coord_nod3d(3,elnodes(4:6))))/3.0
        aux=g*dz*rho_grad/rho0
        PGF(:,k,elem) = p_grad + aux*0.5    
        p_grad=p_grad + aux

     end do !k
  end do ! elem           

end subroutine compute_pressure_force
!
!----------------------------------------------------------------------------
!
function theta(s,t,p,pr)
  ! to compute local potential temperature at pr
  ! using bryden 1973 polynomial for adiabatic lapse rate
  ! and runge-kutta 4-th order integration algorithm.
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! fofonoff,n.,1977,deep-sea res.,24,489-491
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       reference prs   pr       decibars
  !       potential tmp.  theta    deg celsius
  ! checkvalue: theta= 36.89073 c,s=40 (ipss-78),t=40 deg c,
  ! p=10000 decibars,pr=0 decibars
  implicit none
  real*8 			:: theta, s, t, p, pr
  real*8 			:: h, xk, q
  real*8, external	        :: atg

  h = pr - p
  xk = h*atg(s,t,p)
  t = t + 0.5*xk
  q = xk
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  t = t + 0.29289322*(xk-q)
  q = 0.58578644*xk + 0.121320344*q
  xk = h*atg(s,t,p)
  t = t + 1.707106781*(xk-q)
  q = 3.414213562*xk - 4.121320344*q
  p = p + 0.5*h
  xk = h*atg(s,t,p)
  theta = t + (xk-2.0*q)/6.0
  return
end function theta
!
!-----------------------------------------------------------------
!
function atg(s,t,p)
  ! adiabatic temperature gradient deg c per decibar
  ! ref: bryden,h.,1973,deep-sea res.,20,401-408
  ! units:
  !       pressure        p        decibars
  !       temperature     t        deg celsius (ipts-68)
  !       salinity        s        (ipss-78)
  !       adiabatic       atg      deg. c/decibar
  ! checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
  ! t=40 deg c,p0=10000 decibars
  implicit none
  real*8  atg, s, t, p, ds

  ds = s - 35.0
  atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p   &
       +((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t        &
       +8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p             &
       +(-4.2393e-8*t+1.8932e-6)*ds                          &
       +((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5

  return
end function atg

subroutine ocean_mesh_setup
  use o_param
  use o_elements
  use o_mesh
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  if(rotated_grid) call calculate_rotate_matrix

  call read_mesh                            ! Mesh is read (routine in a separate file)

  call mesh_scaling                         ! long., lat. are transf. into rad

  call standard_element_definition_2D       ! Basis functions and scalar
  call standard_element_definition_3D        
  call basisfunctions_3D                    ! Local basis f. deriv. are computed
  call basisfunctions_2D                     

  if(mype==0) write(*,*) 'element volume and basis function derivatives prepared' 

  call check_mesh_quality_resolution

  call build_nghbr_arrays                   ! Builds arrays nod_in_elem2D

  call find_bottom_nodes

  if(grid_type/=1 .or. Redi_GM) call find_layer_elem3d

#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
  call check_mesh_quality_opbd
  call build_open_boundary_arrays 	    ! build OB arrays
#endif

#ifdef use_fullfreesurf
  call find_updating_element    
#endif

end subroutine ocean_mesh_setup
!
!---------------------------------------------------------------------------
!
subroutine check_mesh_quality_resolution
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_parfe
  implicit none

  real(kind=8)   :: max_area, min_area, g_max_area, g_min_area

  ! display minimal and maximal triangle area
  max_area=maxval(voltriangle)
  min_area=minval(voltriangle)
  g_max_area=0.0
  g_min_area=0.0
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE(max_area, g_max_area, &
       1, MPI_DOUBLE_PRECISION,MPI_MAX, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(min_area, g_min_area, &
       1, MPI_DOUBLE_PRECISION,MPI_MIN, &
       MPI_COMM_WORLD, MPIerr)
  if(mype==0) then
     write(*,*) 'The minimal triangle area is: ', g_min_area, ' m**2'
     write(*,*) 'The maximal triangle area is: ', g_max_area, ' m**2'    
  end if

end subroutine check_mesh_quality_resolution
!
!--------------------------------------------------------------
!
subroutine check_mesh_quality_opbd
  ! check open boundary nodes
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_parfe
  implicit none

  integer        :: elem, elnodes(3), ind(3), cnt

  cnt=0
  do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     ind=index_nod3d(nod3d_below_nod2d(1,elnodes))
     if(count(ind==10)==0 .and. any(ind==12)) then
        cnt=cnt+1
     end if
  end do
  if(cnt>0) then
     write(*,*) '-----------------------------------------------------------------------'
     write(*,*) 'Warning:'
     write(*,*) 'There are ', cnt, 'surface triangles without interior nodes at the'
     write(*,*) 'open boundary. Problems can occur when defining/searching open boundary'
     write(*,*) 'nodes. You are suggested to modify your mesh before continuing.'
     write(*,*) '-----------------------------------------------------------------------'
  end if

  ! note: more check could be updated here, e.g., check the quality of triangles (angle size).
end subroutine check_mesh_quality_opbd
!
!--------------------------------------------------------------
!
subroutine mesh_scaling  
  !
  ! Transforms degrees in rad in coord_nod2D(2,myDim_nod2D+eDim_nod2D)     
  ! Constructs the arrays cos_elem2D(myDim_elem2D)
  ! Constructs num_layers_below_nod2D, 
  ! and does transform to rad in coord_nod3D

  use g_config
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  use g_rotate_grid

  implicit none

  integer         :: i,j,ind2,ind3
  integer         :: n,  node, nodeup
  real(kind=8)    :: lon,lat,rlon,rlat

  ! =======================
  !  Lon and lat to radians
  ! =======================
  coord_nod2D(1,:)=coord_nod2D(1,:)*rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*rad

  ! =======================
  ! Mean cos on 2D elements
  ! This sets spherical geometry!
  ! =======================
  allocate(cos_elem2D(myDim_elem2D))
  do i=1, myDim_elem2D  
     cos_elem2D(i)=sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0
  end do
  if(cartesian) cos_elem2D=1.0  

  ! =======================
  ! number of layers 
  ! =======================
  allocate(num_layers_below_nod2D(myDim_nod2D+eDim_nod2D))
  num_layers_below_nod2D=-1
  do n=1,myDim_nod2D+eDim_nod2D
     do j=1,max_num_layers
        node=nod3D_below_nod2D(j,n)
        if (node > 0) then
           num_layers_below_nod2D(n)=num_layers_below_nod2D(n) + 1
        else
           exit
        end if
     end do
  end do
  ! ========================
  ! Lon and lat to radians for 3D nodes
  ! ========================

  do n=1,myDim_nod2D+eDim_nod2D
     !   coord_nod3D correction:
     do j=1,max_num_layers
        node=nod3D_below_nod2D(j,n)
        if (node < 1) exit
        coord_nod3D(1,node)=coord_nod2D(1,n)
        coord_nod3D(2,node)=coord_nod2D(2,n)
     end do
  end do

  ! setup geolat, which contains geographic latitude
  allocate(geolat(ToDim_nod3d))
  if(rotated_grid) then
     do i=1,ToDim_nod3d
        rlon=coord_nod3d(1,i)
        rlat=coord_nod3d(2,i)
        call r2g(lon,lat,rlon,rlat)
        geolat(i)=lat
     end do
  else
     geolat=coord_nod3d(2,:)
  end if

  if(mype==0) write(*,*) 'converting coord. from degree to rad. & defining num. layers done'
end subroutine mesh_scaling
!
!=====================================================================
!
subroutine standard_element_definition_2D
  !
  !    - BASISFUNCTIONS ON 2D STANDARD ELEMENTS 
  !         stdbafunc(1)=1-x-y  stdbafunc(2)=x  stdbafunc(3)=y
  !
  use o_ELEMENTS
  use o_mesh
  !
  implicit none
  !
  integer :: i,j
  !
  Vol2D =  1.0_8/2.0_8               ! Vol2D = < 1,1>
  !        
  sProd_2Di = 1.0_8/6.0_8            ! <stdbafunc(i),1.>
  !
  allocate(sProd_2Dij(3,3))
  sProd_2Dij=1.0_8/24.0_8            ! <stdbafunc(i),stdbafunc(j)>
  do j=1,3
     sProd_2Dij(j,j)=1.0_8/12.0_8 
  end do

  ! Scalar products are only required as sProd_2D/Vol2D:
  sProd_2Dij=sProd_2Dij/Vol2D
  sProd_2Di=sProd_2Di/Vol2D

  !   derivative_stdbafu_x(i,j) = d(Fi(j))/dx(i) on the standard element 
  allocate(derivative_stdbafu_x_2D(2,3))
  !
  derivative_stdbafu_x_2D= 0.
  derivative_stdbafu_x_2D(:,1)= -1.
  derivative_stdbafu_x_2D(1,2)= 1.
  derivative_stdbafu_x_2D(2,3)= 1.
  !
end subroutine standard_element_definition_2D
!
!--------------------------------------------------------------
!
subroutine standard_element_definition_3D
  !
  !    - BASISFUNCTIONS ON 3D STANDARD ELEMENTS (stdbafu)
  !      stdbafunc(1)=1-x-y-z, stdbafunc(2)=x, stdbafunc(3)=y, stdbafunc(4)=z     
  use o_ELEMENTS
  use o_MESH 
  !
  implicit none 
  !
  integer :: i,j
  !
  Vol3D = 1.0_8/6.0_8                  ! Vol3D = < 1.,1.>  
  !
  sProd_3Di = 1.0_8/24.0_8             ! sProd_3Di=<stdbafunc(i),1.>
  !
  allocate(sProd_3Dij(4,4))
  sProd_3Dij=1./120.                   ! sProd_3Dij=<stbafunc(i),stdbafunc(j)>
  do j=1,4
     sProd_3Dij(j,j)=1./60.
  end do


  ! Scalar products are only required as sProd_3D/Vol3D:
  sProd_3Dij=sProd_3Dij/Vol3D
  sProd_3Di=sProd_3Di/Vol3D

  !   derivative_stdbafu_x(j,i) = d(Fi(j))/dx(i) on the standard element
  allocate(derivative_stdbafu_x_3D(3,4))
  !
  do j=1,4
     do i=1,3
        derivative_stdbafu_x_3D(i,j)= 0.0
        if(j==1)   derivative_stdbafu_x_3D(i,j)=-1.
        if(i==j-1) derivative_stdbafu_x_3D(i,j)= 1.
     end do
  end do
  !
end subroutine standard_element_definition_3D
!
!=====================================================================
! 
subroutine basisfunctions_2D
  use o_ELEMENTS
  use o_MESH
  use g_PARFE
  implicit none
  !
  real(kind=8)                           :: DET2D
  real(kind=8), dimension(2,3)           :: derivative_locbafu_x_2D
  real(kind=8), dimension(2,2)           :: jacobian2D, jacobian2D_inv
  integer                                :: elem,i


  allocate(bafux_2d(3,myDim_elem2d), bafuy_2d(3,myDim_elem2d))
  allocate(voltriangle(myDim_elem2d))

  bafux_2d = 0.0
  bafuy_2d = 0.0
  voltriangle = 0.0

  do elem=1,myDim_elem2d
     call local_element_def_2D(elem, DET2D, derivative_locbafu_x_2D)
     do i=1,3
        bafux_2d(i,elem) = derivative_locbafu_x_2D(1,i)
        bafuy_2d(i,elem) = derivative_locbafu_x_2D(2,i)
     enddo
     voltriangle(elem) = abs(DET2D) * Vol2D
  enddo

end subroutine basisfunctions_2D
!
!=====================================================================
! 
subroutine basisfunctions_3D
  !  Derivatives of basis functions and volumes of tetrahedra are computed
  !  through transform from a standard element to a real one.
  !  Local Cartesian metrics is used
  use o_ELEMENTS
  use o_MESH
  use g_PARFE
  implicit none
  !
  real(kind=8), dimension(3,3)           :: jacobian3D
  real(kind=8), dimension(3,3)           :: jacobian3D_inv
  real(kind=8)                           :: DET3D
  real(kind=8), dimension(3,4)           :: derivative_locbafu_x_3D
  integer                                :: elem,i

  allocate(bafux_3d(4,myDim_elem3d), bafuy_3d(4,myDim_elem3d))
  allocate(bafuz_3d(4,myDim_elem3d),voltetra(myDim_elem3d))
  bafux_3d = 0.0
  bafuy_3d = 0.0
  bafuz_3d = 0.0
  voltetra = 0.0

  do elem=1,myDim_elem3d
     call local_element_def_3D(elem, DET3D, derivative_locbafu_x_3D)
     do i=1,4
        bafux_3d(i,elem) = derivative_locbafu_x_3D(1,i)
        bafuy_3d(i,elem) = derivative_locbafu_x_3D(2,i)
        bafuz_3d(i,elem) = derivative_locbafu_x_3D(3,i)
     enddo
     voltetra(elem) = abs(DET3D) * Vol3D
  enddo

end subroutine basisfunctions_3D
!
! =======================================================================
!
subroutine local_element_def_2D(element, DET, derivative_locbafu_x_2D)

  use o_ELEMENTS
  use o_MESH
  use o_param
  use g_config
  implicit none
  !
  integer, intent(IN)                        :: element
  real(kind=8), dimension(2,2)               :: jacobian2D
  real(kind=8), dimension(2,2)               :: jacobian2D_inv
  real(kind=8), intent(OUT)                  :: DET
  real(kind=8), dimension(2,3), intent(OUT)  :: derivative_locbafu_x_2D
  !
  real(kind=8), dimension(2,3)               :: local_cart
  real(kind=8), dimension(3,2)               :: der_transp
  integer                                    :: i, node
  real(kind=8)                               :: meancos
  !
  meancos=cos_elem2D(element)
  do i=1,3
     node=elem2D_nodes(i,element)
     !
     !  scaled cartesian coordinates
     !
     local_cart(1,i)=coord_nod2D(1,node) 
     local_cart(2,i)=r_earth * coord_nod2D(2,node) 
  end do
  !
  !  jacobian
  !
  do i=1,2
     jacobian2D(:,i)= local_cart(:,i+1)-local_cart(:,1)
     if (jacobian2D(1,i)> domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)-domain_length
     if (jacobian2D(1,i)<-domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)+domain_length
  end do
  jacobian2D(1,:)=jacobian2D(1,:)*meancos *r_earth
  !
  !  inverse of jacobian
  !
  call matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)
  !
  der_transp=matmul(transpose(derivative_stdbafu_x_2D), jacobian2D_inv)
  derivative_locbafu_x_2D=transpose(der_transp)
  !
  !
end subroutine local_element_def_2D
!
!============================================================================
!
subroutine local_element_def_3D(element, DET, derivative_locbafu_x_3D)
  !  Auxilliary routine to basisfunctions_3D
  use o_ELEMENTS
  use o_MESH
  use o_param
  use g_config
  !
  implicit none
  !
  integer, intent(IN)                                 :: element
  real(kind=8), dimension(3,3)                        :: jacobian3D
  real(kind=8), dimension(3,3)                        :: jacobian3D_inv
  real(kind=8), intent(OUT)                           :: DET
  real(kind=8), dimension(3,4), intent(OUT)           :: derivative_locbafu_x_3D
  !
  real(kind=8), dimension(3,4)                        :: local_cart
  real(kind=8), dimension(4,3)                        :: der_transp
  integer                                             :: i, node
  real(kind=8)                                        :: meancos 

  !
  !  scaled cartesian coordinates on elements
  !

  do i=1,4
     node = elem3D_nodes(i,element)
     local_cart(1,i) = coord_nod3D(1,node)
     local_cart(2,i) = r_earth * coord_nod3D(2,node)
     local_cart(3,i) = coord_nod3D(3,node)
  end do
  meancos=cos_elem2D(elem2D_corresp_to_elem3D(element))
  !
  !  transformation matrix Xl(k)=jacobian(i,k)*Xs(i)
  !
  do i=1,3
     jacobian3D(:,i) = local_cart(:,i+1)-local_cart(:,1)
     if (jacobian3D(1,i)> domain_length/2.0) jacobian3D(1,i)=jacobian3D(1,i)-domain_length
     if (jacobian3D(1,i)<-domain_length/2.0) jacobian3D(1,i)=jacobian3D(1,i)+domain_length
  end do

  jacobian3D(1,:)=jacobian3D(1,:)*meancos*r_earth
  !
  !  inverse jacobian
  !
  call matrix_inverse_3x3(jacobian3D, jacobian3D_inv, DET)
  !
  !  derivatives of local basis functions
  !
  der_transp=matmul(transpose(derivative_stdbafu_x_3D), jacobian3D_inv)
  derivative_locbafu_x_3D=transpose(der_transp)

end subroutine local_element_def_3D
!
!=======================================================================
!
subroutine  matrix_inverse_2x2 (A, AINV, DET)
  !
  !
  implicit none
  !
  real(kind=8), dimension(2,2), intent(IN)  :: A
  real(kind=8), dimension(2,2), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0 )  then
     do j=1,2
        write(*,*) (A(i,j),i=1,2)
     end do
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
!
!=======================================================================
!
subroutine matrix_inverse_3x3(A, AINV, DET)
  !
  implicit none
  !
  real(kind=8), dimension(3,3), intent(IN)  :: A
  real(kind=8), dimension(3,3), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  AINV(1,1) =  A(2,2)*A(3,3) - A(3,2)*A(2,3)
  AINV(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
  AINV(3,1) =  A(2,1)*A(3,2) - A(3,1)*A(2,2)
  AINV(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
  AINV(2,2) =  A(1,1)*A(3,3) - A(3,1)*A(1,3)
  AINV(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
  AINV(1,3) =  A(1,2)*A(2,3) - A(2,2)*A(1,3)
  AINV(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
  AINV(3,3) =  A(1,1)*A(2,2) - A(2,1)*A(1,2)
  DET = A(1,1)*AINV(1,1) + A(1,2)*AINV(2,1) + A(1,3)*AINV(3,1)
  !
  if ( DET .eq. 0.0 )  then
     do j=1,3
        write(*,*) (A(i,j),i=1,3)
     end do
     stop 'SINGULAR 3X3 MATRIX'
  else
     AINV = AINV/DET
  endif
  !
end subroutine matrix_inverse_3x3
!
! =======================================================================
!
subroutine build_nghbr_arrays
  !
  ! Assembles additional arrays which list for each node the elements 
  ! containing the node and node neighbours

  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none

  integer                            :: j,k,m,a,tr(3),tet(4), counter, el,ml(1)
  integer, allocatable, dimension(:) :: ind
  integer, dimension(100)            :: AUX=0

  !--------------- 2D mesh:
  ! Builds nod_in_elem2D
     
  allocate(ind(myDim_nod3D+eDim_nod3D))

  ind=0
  do j=1,myDim_elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(myDim_nod2D))
  nod_in_elem2D%nmb=ind(1:myDim_nod2D)    
  do j=1,myDim_nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,myDim_elem2D   
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
     do k=1,3
        if(tr(k)<=myDim_nod2D) then 
           nod_in_elem2D(tr(k))%addresses(ind(tr(k)))=j
        end if
     end do
  end do
  ! the list of elements is ordered, and no sorting is needed

  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(myDim_nod2D))
  ind=0
  do j=1, myDim_nod2D
     counter=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=counter
     allocate(nghbr_nod2D(j)%addresses(counter))

     ! we need to sort array aux(1:counter)
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod2D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

  !---------------3D mesh---------------------
  ! Construction of nod_in_elem3D
  ind=0
  do j=1,myDim_elem3D
     tet=elem3D_nodes(:,j)
     ind(tet)=ind(tet)+1
  end do
  allocate(nod_in_elem3D(myDim_nod3D))
  nod_in_elem3D%nmb=ind(1:myDim_nod3D)
  do j=1,myDim_nod3D   
     allocate(nod_in_elem3D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,myDim_elem3D   
     tet=elem3D_nodes(:,j)
     ind(tet)=ind(tet)+1
     do k=1,4
        if(tet(k)<=myDim_nod3D) then                         
           nod_in_elem3D(tet(k))%addresses(ind(tet(k)))=j     
        end if
     end do
  end do
  ! the list of elements is ordered, and no sorting is needed

  ! Builds nghbr_nod3D

  allocate(nghbr_nod3D(myDim_nod3D))

  ind=0
  do j=1, myDim_nod3D
     counter=0
     do m=1,nod_in_elem3D(j)%nmb
        el=nod_in_elem3D(j)%addresses(m)
        do k=1, 4
           a=elem3D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod3D(j)%nmb=counter
     allocate(nghbr_nod3D(j)%addresses(counter))

     ! we need to sort array(aux(1:counter))
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod3D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

  deallocate(ind)

  if(mype==0) write(*,*) 'node and element neighborhood prepared'
end subroutine build_nghbr_arrays
!
!==========================================================================
!
subroutine find_bottom_nodes
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  integer         :: j
  !
  allocate(bt_nds(myDim_nod2D+eDim_nod2D))
  do j=1, myDim_nod2D+eDim_nod2D
     bt_nds(j)=nod3D_below_nod2D(num_layers_below_nod2D(j)+1,j)
  end do
end subroutine find_bottom_nodes
!
!==========================================================================
!
subroutine find_cluster_area
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  use g_rotate_grid
  implicit none
  
  integer         :: i, elem, elnodes(3)
  real(kind=8)    :: inv3, vol

  inv3=1.0/3.0_8

  allocate(cluster_area_2d(ToDim_nod2D))
  cluster_area_2d=0.0

  do elem=1,myDim_elem2d
     elnodes=elem2d_nodes(:,elem)
     vol=voltriangle(elem)*inv3
     cluster_area_2d(elnodes)=cluster_area_2d(elnodes)+vol
  end do

  call com_2d(cluster_area_2d)

  vol=0.0
  do i=1,myDim_nod2d
     vol=vol+cluster_area_2d(i)
  end do

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  ocean_area=0.0
  call MPI_AllREDUCE(vol, ocean_area, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  
end subroutine find_cluster_area
!
!==========================================================================
!
subroutine find_layer_elem3d
  !find the layer number of 3d elements
  use o_MESH
  use o_ELEMENTS
  use g_parfe

  implicit none

  integer                   :: i, j, k, elem3, tetra_nodes(4)
  integer, allocatable      :: auxind(:)

  allocate(elem3d_layer(myDim_elem3D))
  allocate(auxind(myDim_nod3d+eDim_nod3d))	

  do i=1,myDim_nod2d+eDim_nod2d
     do k=1,max_num_layers
        j=nod3d_below_nod2d(k,i)  
        if(j<0) exit
        auxind(j)=k
     end do
  end do
  do elem3=1, myDim_elem3D
     tetra_nodes=elem3d_nodes(:,elem3)
     elem3d_layer(elem3)=minval(auxind(tetra_nodes))
  end do

  deallocate(auxind)
end subroutine find_layer_elem3d
!
!==========================================================================
!
subroutine build_open_boundary_arrays
  use o_MESH
  use o_ELEMENTS
  use o_param
  use o_array
  use g_config
  use g_PARFE
  implicit none
  integer         	:: j, k, m, cnt, col, nodes(4), nodes2(3), nodes23(4)
  integer         	:: tri(3,4), fourth_node(4), ind(4)
  integer               :: ind2(3), edg(2,3), third_node(3)
  real(kind=8)       	:: vec1(3), vec2(3), nvec(3), vec_dir(3)
  real(kind=8)          :: meancos, node_cart(3,4)
  real(kind=8)          :: vec2d(2), nvec2d(2), vec_dir2d(2), node_cart2(2,3)
  integer, allocatable, dimension(:,:) :: auxtri, auxedg
  integer, allocatable, dimension(:)   :: aux, aux2

  ! =============
  ! First, find the list of triangles on the open boundary
  ! =============

  cnt=0
  do j=1, myDim_elem3D        
     nodes=elem3D_nodes(:,j)
     nodes23=nod2D_corresp_to_nod3D(nodes)        
     ind=index_nod3D(nod3D_below_nod2D(1,nodes23))
     if (count(ind==12)==0) cycle
     ! form all possible triangles
     tri(:,1)=nodes(1:3)
     tri(:,2)=nodes(2:4)
     tri(:,3)=nodes((/1, 3, 4/))
     tri(:,4)=nodes((/1, 2, 4/))
     fourth_node(1)=nodes(4)
     fourth_node(2)=nodes(1)
     fourth_node(3)=nodes(2)
     fourth_node(4)=nodes(3)
     do k=1, 4
        !(1) if tri(:,k) belongs to open boundary it should project into
        ! a line at the surface (the triangle lies in the vertical plane).

        nodes2= nod2D_corresp_to_nod3D(tri(:,k))
	if((nodes2(1)==nodes2(2)).or.(nodes2(1)==nodes2(3)).or.(nodes2(2)==nodes2(3))) then

           !(2) there are no wet nodes and at least one node belongs
           !to the open boundary

           ind(1:3)=index_nod3D(nod3d_below_nod2D(1,nodes2))
           if((count(ind(1:3)/=10)==3).and.(count(ind(1:3)==12)>0)) then
              ! (3) the segment given by two different nodes from nodes2
              !    can be shared by 2 triangles (in the corners on bad
              !    meshes; this however implies that there are triangles
              !    formed by boundary nodes). We assume that such triangles are
              !    absent. 
              cnt=cnt+1
              exit
           end if
	end if
     end do
  end do

  allocate(auxtri(4,cnt))
  cnt=0	
  do j=1, myDim_elem3D          
     nodes=elem3D_nodes(:,j)
     nodes23=nod2D_corresp_to_nod3D(nodes)        
     ind=index_nod3D(nod3D_below_nod2D(1,nodes23))
     if (count(ind==12)==0) cycle
     ! form all possible triangles
     tri(:,1)=nodes(1:3)
     tri(:,2)=nodes(2:4)
     tri(:,3)=nodes((/1, 3, 4/))
     tri(:,4)=nodes((/1, 2, 4/))
     fourth_node(1)=nodes(4)
     fourth_node(2)=nodes(1)
     fourth_node(3)=nodes(2)
     fourth_node(4)=nodes(3)
     do k=1, 4
        !(1) if tri(:,k) belongs to open boundary it should project into
        !a line at the surface (the triangle lies in the vertical plane).

        nodes2= nod2D_corresp_to_nod3D(tri(:,k))
	if((nodes2(1)==nodes2(2)).or.(nodes2(1)==nodes2(3)).or.(nodes2(2)==nodes2(3))) then
           !(2) there are no wet nodes and at least one node belongs
           !to the open boundary

           ind(1:3)=index_nod3D(nod3d_below_nod2D(1,nodes2))
           if((count(ind(1:3)/=10)==3).and.(count(ind(1:3)==12)>0)) then
              ! (3) the segment given by two different nodes from nodes2
              !    can be shared by 2 triangles (in the corners on bad
              !    meshes; this however implies that there are triangles
              !    formed by boundary nodes). We assume that such triangles are
              !    absent. 
              cnt=cnt+1
              auxtri(1:3,cnt)=tri(:,k)
              auxtri(4,cnt)=fourth_node(k)
              exit
           end if
	end if
     end do
  end do

  nmbr_opbnd_tri=cnt
  !if(nmbr_opbnd_tri==0) then
  !   deallocate(auxtri)     
  !   return                           
  !end if
  allocate(opbnd_tri(nmbr_opbnd_tri,4))
  do j=1,nmbr_opbnd_tri
     opbnd_tri(j,:)=auxtri(:,j)
  end do
  deallocate(auxtri)  

  ! =============
  ! Second, compute outer normal to open boundary triangles and triangle areas
  ! =============
  allocate(opbnd_nv(nmbr_opbnd_tri,4))
  do j=1, nmbr_opbnd_tri
     ! (a) compute cartesian coordinates
     do m=1,4
        node_cart(:,m)=coord_nod3D(:,opbnd_tri(j,m))
     end do
     meancos=0.
     do m=1,3
        meancos=meancos+cos(node_cart(2,m))
     end do
     meancos=meancos/3.0

     ! (b) two vectors in the triangle plane and the vector to the fourth node
     vec1=node_cart(:,2)-node_cart(:,1)
     vec2=node_cart(:,3)-node_cart(:,1)
     vec_dir=node_cart(:,4)-node_cart(:,1)

     ! (c) check for cyclicity:
     if(vec1(1)>domain_length/2.0_8) vec1(1)=vec1(1)-domain_length
     if(vec1(1)<-domain_length/2.0_8) vec1(1)=vec1(1)+domain_length
     if(vec2(1)>domain_length/2.0_8) vec2(1)=vec2(1)-domain_length
     if(vec2(1)<-domain_length/2.0_8) vec2(1)=vec2(1)+domain_length
     if(vec_dir(1)>domain_length/2.0_8) vec_dir(1)=vec_dir(1)-domain_length
     if(vec_dir(1)<-domain_length/2.0_8) vec_dir(1)=vec_dir(1)+domain_length

     vec1(1:2)=vec1(1:2)*r_earth
     vec2(1:2)=vec2(1:2)*r_earth
     vec_dir(1:2)=vec_dir(1:2)*r_earth
     ! Local Cartesian coordinates:
     vec1(1)=vec1(1)*meancos
     vec2(1)=vec2(1)*meancos
     vec_dir(1)=vec_dir(1)*meancos

     ! (d) compute vector product of vec1 and vec2
     nvec(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
     nvec(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
     nvec(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
     opbnd_nv(j,4)=0.5*sqrt(sum(nvec*nvec))  ! area of surface tri in sca coord.
     ! Return to non-scaled nvec:
     nvec=nvec/sqrt(sum(nvec*nvec))

     ! (e) outer normal
     if (sum(nvec*vec_dir) > 0) nvec = -nvec
     opbnd_nv(j,1:3)=nvec
  end do

  ! ==================
  ! Third, find the list of 3d and 2d nodes on the open boundary
  ! ==================
  ! list of 3d nodes on the open boundary
  allocate(aux2(myDim_nod3D+eDim_nod3D))  
  aux2=0
  do j=1, nmbr_opbnd_tri
     do m=1,3
        if(aux2(opbnd_tri(j,m))==0) aux2(opbnd_tri(j,m))=1
     end do
  end do
  nmbr_opbnd_n3D=sum(aux2)
  allocate(aux(nmbr_opbnd_n3D))
  cnt=0
  do j=1, myDim_nod3D+eDim_nod3D           
     if(aux2(j)/=0) then
        cnt=cnt+1
        aux(cnt)=j
     end if
  end do
  allocate(opbnd_n3D(nmbr_opbnd_n3D))
  opbnd_n3D=aux

  ! number of 2D nodes that have index_nod3D=12 on opbnd
  cnt=0
  do j=1, myDim_nod2D+eDim_nod2D        
     m=nod3D_below_nod2D(1,j)           
     if (index_nod3D(m)==12) cnt=cnt+1
  end do
  nmbr_opbnd_n2D=cnt
  ! number of 2D nodes that have index_nod3D=11 on opbnd
  cnt=0
  do j=1,nmbr_opbnd_n3D
     col=aux(j)
     if (index_nod3D(col)==11) cnt=cnt+1
  end do
  nmbr_opbnd_t2D=cnt
  ! total 2D nodes on opbnd
  nmbr_opbnd_t2D=nmbr_opbnd_n2D+nmbr_opbnd_t2D

  allocate(opbnd_n2D(nmbr_opbnd_t2D))
  allocate(mapping_opbnd_n2d(myDim_nod2d+eDim_nod2D)) 
  cnt=0
  do j=1, myDim_nod2D+eDim_nod2D    
     m=nod3D_below_nod2D(1,j)       
     if(index_nod3D(m)==12) then
        cnt=cnt+1
        opbnd_n2D(cnt)=j
        mapping_opbnd_n2d(j)=cnt
     end if
  end do
  cnt=nmbr_opbnd_n2D
  do j=1,nmbr_opbnd_n3D
     col=aux(j)
     if (index_nod3D(col)==11) then
        col=nod2d_corresp_to_nod3d(col)
        cnt=cnt+1
        opbnd_n2D(cnt)=col
        mapping_opbnd_n2d(col)=cnt
     end if
  end do
  deallocate(aux, aux2)

  ! ================
  ! 4th, depth at the open boundary nodes
  ! ================
  allocate(opbnd_dep(nmbr_opbnd_t2d))
  do k=1, nmbr_opbnd_t2d
     m=opbnd_n2d(k)
     cnt=num_layers_below_nod2D(m)+1
     j=nod3d_below_nod2d(cnt,m)
     opbnd_dep(k)=abs(coord_nod3d(3,j))
  end do


#ifdef use_ice
  ! ================
  ! 5th, find the list of edges on the open boundary
  ! ================
  allocate(auxedg(3,myDim_elem2D))
  cnt=0
  do j=1, myDim_elem2D
     nodes2=elem2D_nodes(:,j)
     nodes(1:3)=nod3D_below_nod2D(1,nodes2)   
     ind2=index_nod3D(nodes(1:3))             
     if (count(ind2==12)==0) cycle
     ! form all possible edges
     edg(:,1)=nodes2(1:2)
     edg(:,2)=nodes2(2:3)
     edg(:,3)=nodes2((/3, 1/))

     third_node(1)=nodes2(3)
     third_node(2)=nodes2(1)
     third_node(3)=nodes2(2)

     do k=1, 3
        ind2(1:2)=index_nod3d(nod3D_below_nod2D(1,edg(:,k)))  
        if(count(ind2(1:2)/=10)==2) then
	   if(count(ind(1:2)==12)>0) then
              cnt=cnt+1
              auxedg(1:2,cnt)=edg(:,k)
              auxedg(3,cnt)=third_node(k) 
              exit
           end if
        end if
     end do
  end do
  nmbr_opbnd_edg=cnt
  allocate(opbnd_edg(nmbr_opbnd_edg,3))
  do j=1,nmbr_opbnd_edg
     opbnd_edg(j,:)=auxedg(:,j)
  end do
  deallocate(auxedg)  

  ! ==================
  ! 6th, compute outer normal to open boundary edges and edge length
  ! ==================
  allocate(opbnd_edg_nv(nmbr_opbnd_edg, 3))
  do j=1, nmbr_opbnd_edg
     ! (a) compute cartesian coordinates
     do m=1,3
        node_cart2(:,m)=coord_nod2D(:,opbnd_edg(j,m))
     end do
     meancos=0.
     do m=1,2
        meancos=meancos+cos(node_cart2(2,m))
     end do
     meancos=meancos/2.0

     ! (b) vectors of the edge and the vector to the third node
     vec2d=node_cart2(:,2)-node_cart2(:,1)
     vec_dir2d=node_cart2(:,3)-node_cart2(:,1)

     ! (c) check for cyclicity:
     if(vec2d(1)>domain_length/2.0_8) vec2d(1)=vec2d(1)-domain_length
     if(vec2d(1)<-domain_length/2.0_8) vec2d(1)=vec2d(1)+domain_length
     if(vec_dir2d(1)>domain_length/2.0_8) vec_dir2d(1)=vec_dir2d(1)-domain_length
     if(vec_dir2d(1)<-domain_length/2.0_8) vec_dir2d(1)=vec_dir2d(1)+domain_length

     vec2d(1:2)=vec2d(1:2)*r_earth
     vec_dir2d(1:2)=vec_dir2d(1:2)*r_earth
     ! Local Cartesian coordinates:
     vec2d(1)=vec2d(1)*meancos
     vec_dir2d(1)=vec_dir2d(1)*meancos

     ! (d) normal vector
     nvec2d(1)=vec2d(2)
     nvec2d(2)=vec2d(1)
     opbnd_edg_nv(j,3)=sqrt(nvec2d(1)**2+nvec2d(2)**2)
     nvec2d=nvec2d/opbnd_edg_nv(j,3)

     ! (e) outer normal
     if (sum(nvec2d*vec_dir2d) > 0.0) nvec2d = -nvec2d
     opbnd_edg_nv(j,1:2)=nvec2d
  end do
#endif

  if(mype==0) write(*,*) 'open boundary grid prepared'

end subroutine build_open_boundary_arrays
!
!=========================================================================
!
#ifdef use_fullfreesurf
subroutine find_updating_element
  use o_param
  use o_ELEMENTS
  use o_MESH
  use g_parfe
  implicit none
  !
  integer                                :: el, n, n_el, elnodes(4)
  !
  allocate(map_elem(myDim_elem3d))   
  map_elem=0
  n=0
  do el=1,myDim_elem3d              
     elnodes=elem3d_nodes(:,el)
     if(any(index_nod3D(elnodes)<=12)) then  
        n=n+1
        map_elem(el)=n
     end if
  end do
  !
  allocate(voltetra_new(n), bafux_3d_new(4,n))
  allocate(bafuy_3d_new(4,n), bafuz_3d_new(4,n))
  !
  do el=1,myDim_elem3D
     n_el=map_elem(el)
     if(n_el==0) cycle
     voltetra_new(n_el)=voltetra(el)
     bafux_3d_new(:,n_el)=bafux_3d(:,el)
     bafuy_3d_new(:,n_el)=bafuy_3d(:,el)
     bafuz_3d_new(:,n_el)=bafuz_3d(:,el)
  end do
  
  if(mype==0) write(*,*) 'index of first layer of elements recorded' 
end subroutine find_updating_element
#endif
!
!--------------------------------------------------------------------------
!
#ifdef use_fullfreesurf
subroutine update_mesh
  use o_param
  use o_array
  use o_mesh
  use o_ELEMENTS
  use g_parfe
  implicit none
  !
  integer                              :: m, el, n_el, i, row
  real(kind=8), dimension(3,3)         :: jacobian3D
  real(kind=8), dimension(3,3)         :: jacobian3D_inv
  real(kind=8), dimension(3,4)         :: derivative_locbafu_x_3D
  real(kind=8)                         :: DET3D
  
  ! update z coordiante
  do m=1,myDim_nod2d+eDim_nod2d     
     row=nod3d_below_nod2d(1,m)
     coord_nod3d(3,row)=ssh(m) 
  end do
  
  ! update new vol and derivatives
  do el=1, myDim_elem3d                
     n_el=map_elem(el)
     if(n_el==0) cycle
     bafux_3d(:,el)=bafux_3d_new(:,n_el)
     bafuy_3d(:,el)=bafuy_3d_new(:,n_el)
     bafuz_3d(:,el)=bafuz_3d_new(:,n_el)
     voltetra(el)=voltetra_new(n_el)
     call local_element_def_3D(el, DET3D, derivative_locbafu_x_3D)
     do i=1,4
        bafux_3d_new(i,n_el) = derivative_locbafu_x_3D(1,i)
        bafuy_3d_new(i,n_el) = derivative_locbafu_x_3D(2,i)
        bafuz_3d_new(i,n_el) = derivative_locbafu_x_3D(3,i)
     enddo
     voltetra_new(n_el) = abs(DET3D)*Vol3D
  enddo
end subroutine update_mesh
#endif
!--------------------------------------------------------------
! Reads mesh and communication information in a distributed way
! some extra info. (nodal flag for cavity, region type
! of 2d elements, sigma grid slope) is also read in here.

subroutine read_mesh
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use o_ARRAY
  use g_config
  use g_PARfe 
  implicit none

  integer        n, m, fileID, ind, nini, nend, n1, n2, n3, n4
  integer        vert_nodes(120)
  real(kind=8)   x, y, z
  character*10   mype_string
  character*80   file_name
  character*80   dist_mesh_dir

  write(mype_string,'(i4.4)') mype  
  dist_mesh_dir=trim(meshpath)//'dist/'

  !=======================
  ! rank partitioning vectors
  !=======================
  file_name=trim(dist_mesh_dir)//'rpart.out' 
  fileID=10+mype
  open(fileID, file=trim(file_name)) 
  allocate(part2D(npes+1), part3D(npes+1))

  read(fileID,*) n
  if (n.ne.npes) then
     write(*,*) 'current NPES does not coincide with that used for mesh pre-partition'
     call par_ex
     stop
  end if
  part2D(1)=1
  read(fileID,*) part2D(2:npes+1)
  do n=2, npes+1
     part2D(n)=part2D(n-1)+part2D(n)
  end do

  part3D(1)=1
  read(fileID,*) part3D(2:npes+1)
  do n=2, npes+1
     part3D(n)=part3D(n-1)+part3D(n)
  end do
  close(fileID)
  !write(*,*) 'rpart is read'

  !===========================
  ! Lists of nodes and elements 
  ! in global indexing. Not everything
  ! is needed
  !===========================

  file_name=trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
  fileID=10+mype  

  open(fileID, file=trim(file_name))
  read(fileID,*) n

  read(fileID,*) myDim_nod2D
  read(fileID,*) eDim_nod2D
  ToDim_nod2d=myDim_nod2d+eDim_nod2d
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
  read(fileID,*) myList_nod2D

  read(fileID,*) myDim_nod3D
  read(fileID,*) eDim_nod3D 	 
  ToDim_nod3d=myDim_nod3d+eDim_nod3d
  allocate(myList_nod3D(myDim_nod3D+eDim_nod3D)) 	 
  read(fileID,*) myList_nod3D

  read(fileID,*) myDim_elem2D
  allocate(myList_elem2D(myDim_elem2D))
  read(fileID,*) myList_elem2D

  read(fileID,*) myDim_elem3D
  allocate(myList_elem3D(myDim_elem3D))
  read(fileID,*) myList_elem3D ! m

  close(fileID)

  !==============================
  ! Allocate mapping array
  !==============================
  nod3D=part3D(npes+1)-1
  nod2D=part2D(npes+1)-1
  file_name=trim(meshpath)//'elem3d.out'
  open(fileID, file=file_name)
  read(fileID,*) elem3D
  close(fileID)
  allocate(mapping(elem3D))  
  mapping=0 
  !==============================
  ! It will be used for several purposes 
  ! and finally will be filled with a correct mapping 
  !==============================
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  !==============================
  ! read 2d node data
  !==============================

  allocate(coord_nod2D(2,myDim_nod2D+eDim_nod2D))
  allocate(index_nod2D(myDim_nod2D+eDim_nod2D))	 			

  file_name=trim(meshpath)//'nod2d.out' 
  open(fileID, file=file_name)
  read(fileID,*) n      ! nod2D, we know it already

  do n=1,nod2D
     read(fileID,*) m, x, y, ind
     if (mapping(n)>0) then
        coord_nod2D(1,mapping(n))=x
        coord_nod2D(2,mapping(n))=y
        index_nod2D(mapping(n))=ind
     end if
  end do
  mapping(1:nod2D)=0
  close(fileID)

  !==============================
  ! read 2d elem data
  !==============================
  file_name=trim(meshpath)//'elem2d.out' 
  open(fileID, file=file_name)

  allocate(elem2D_nodes(3, myDim_elem2D))
  do n=1, myDim_elem2D
     mapping(myList_elem2D(n))=n
  end do
  read(fileID,*) elem2d    
  do n=1,elem2D
     read(fileID,*) n1, n2, n3
     if (mapping(n)>0) then
        elem2D_nodes(1,mapping(n))=n1
        elem2D_nodes(2,mapping(n))=n2
        elem2D_nodes(3,mapping(n))=n3
     end if
  end do
  close(fileID)
  ! nodes in elem2d are in global numbering. convert to local:

  mapping(1:elem2D)=0
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do
  do n=1, myDim_elem2D
     do m=1,3
        n1=elem2D_nodes(m,n)	 
        elem2D_nodes(m,n)=mapping(n1)	 
     end do
  end do
  mapping(1:nod2D)=0

  !==============================
  ! Ice shelf variables
  !==============================
  allocate(cavity_flag_nod2d(myDim_nod2d+eDim_nod2d))
  cavity_flag_nod2d=0
#ifdef use_cavity
  file_name=trim(meshpath)//'cavity_flag_nod2d.out'
  open(fileID, file=file_name)
  do n=1, myDim_nod2D+eDim_nod2d
     mapping(myList_nod2D(n))=n
  end do
  do n=1,nod2D
     read(fileID,*) n1
     if(mapping(n)>0) cavity_flag_nod2d(mapping(n))=n1
  end do
  mapping(1:nod2D)=0
  close(fileID)
#endif

  !=============================
  ! Region type of 2d elements: z-level or sigma
  !=============================
  allocate(grid_type_elem2d(myDim_elem2d))
  if(grid_type==1) then
     grid_type_elem2d=0
  elseif(grid_type==2) then
     grid_type_elem2d=1
  else
     file_name=trim(meshpath)//'grid_type_elem2d.out'
     open(fileID, file=file_name)
     do n=1, myDim_elem2D
        mapping(myList_elem2D(n))=n
     end do
     do n=1,elem2D
        read(fileID,*) n1
	if(mapping(n)>0) grid_type_elem2d(mapping(n))=n1
     end do
     mapping(1:elem2D)=0
  end if
  close(fileID)

  !==============================
  ! read 3d node data
  !==============================

  allocate(coord_nod3D(3,myDim_nod3D+eDim_nod3D))
  allocate(index_nod3D(myDim_nod3D+eDim_nod3D))	 			
  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do

  file_name=trim(meshpath)//'nod3d.out' 
  open(fileID, file=file_name)
  read(fileID,*) n      ! nod3D, we know it already

  do n=1,nod3D
     read(fileID,*) m, x, y, z, ind
     if (mapping(n)>0) then
        coord_nod3D(1,mapping(n))=x
        coord_nod3D(2,mapping(n))=y
        coord_nod3D(3,mapping(n))=z
        index_nod3D(mapping(n))=ind
     end if
  end do

  mapping(1:nod3D)=0
  close(fileID)

  !==============================
  ! read 3d elem data
  !==============================

  file_name=trim(meshpath)//'elem3d.out' 
  open(fileID, file=file_name)

  allocate(elem3D_nodes(4, myDim_elem3D))
  do n=1, myDim_elem3D
     mapping(myList_elem3D(n))=n
  end do
  read(fileID,*) elem3d    
  do n=1,elem3D
     read(fileID,*) n1, n2, n3, n4
     if (mapping(n)>0) then
        elem3D_nodes(1,mapping(n))=n1
        elem3D_nodes(2,mapping(n))=n2
        elem3D_nodes(3,mapping(n))=n3
        elem3D_nodes(4,mapping(n))=n4
     end if
  end do

  ! nodes in elem3d are in natural numbering. convert to local:

  mapping(1:elem3D)=0
  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do
  do n=1, myDim_elem3D
     do m=1,4
        n1=elem3D_nodes(m,n)	 
        elem3D_nodes(m,n)=mapping(n1)	 
     end do
  end do
  mapping(1:nod3D)=0
  close(fileID)

  !==============================
  ! read aux. arrays 
  !==============================

  file_name=trim(meshpath)//'aux3d.out'  
  open(fileID, file=file_name)

  read(fileID, *) max_num_layers

  !=============================
  ! nod3D_below_nod2D
  !============================= 
  ! ATTENTION: the array is to be stored in slices of max_num_layers
  !            or as a column   
  allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D+eDim_nod2D))       
  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  do n=1, nod2D
     read(fileID, *) vert_nodes(1:max_num_layers)
     if (mapping(n)>0)  then
        nod3D_below_nod2D(:,mapping(n))=vert_nodes(1:max_num_layers)
     end if
  end do
  mapping(1:nod2D)=0

  do n=1, myDim_nod3D+eDim_nod3D
     mapping(myList_nod3D(n))=n
  end do

  do n=1, myDim_nod2D+eDim_nod2D
     do m=1, max_num_layers
        n1=nod3D_below_nod2D(m,n)
        if(n1>0)  then
           nod3D_below_nod2D(m,n)=mapping(n1)
        end if
     end do
  end do

  !=============================
  ! nod2D_corresp_to_nod3D
  !============================= 

  allocate(nod2D_corresp_to_nod3D(myDim_nod3D+eDim_nod3D)) 

  do n=1, nod3D
     read(fileID, *) m
     if (mapping(n)>0)  then
        nod2D_corresp_to_nod3D(mapping(n))=m
     end if
  end do
  mapping(1:nod3D)=0

  do n=1, myDim_nod2D+eDim_nod2D
     mapping(myList_nod2D(n))=n
  end do

  do n=1,myDim_nod3D+eDim_nod3D
     m=nod2D_corresp_to_nod3D(n)
     nod2D_corresp_to_nod3D(n)=mapping(m)
  end do
  mapping(1:nod2D)=0

  !=============================
  ! elem2D_corresp_to_elem3D
  !=============================

  allocate(elem2D_corresp_to_elem3D(myDim_elem3D)) 
  do n=1, myDim_elem3D
     mapping(myList_elem3D(n))=n
  end do
  do n=1,elem3D
     read(fileID,*) m
     if(mapping(n)>0) then
        elem2D_corresp_to_elem3D(mapping(n))=m
     end if
  end do
  mapping(1:elem3D)=0
  do n=1, myDim_elem2D
     mapping(myList_elem2D(n))=n
  end do
  do n=1,myDim_elem3D
     m=elem2D_corresp_to_elem3D(n)
     elem2D_corresp_to_elem3D(n)=mapping(m)
  end do
  mapping(1:elem2D)=0

  close(fileID)

  !correcting 3d nodal indices
#ifndef use_opbnd_restoring
#ifndef use_opbnd_tide
  do n=1,myDim_nod3d+eDim_nod3D
     if (index_nod3D(n)==12) index_nod3D(n)=11
     if (index_nod3D(n)==22) index_nod3D(n)=21
     if (index_nod3D(n)==32) index_nod3D(n)=31
  end do
#endif
#endif

  ! ==============================
  ! read sigma slope and define layer of elem.
  ! ==============================
  if(grid_type/=1) call read_grid_slope	
   

  if(mype==0) then
     write(*,*) 'mesh (according to pre-partition) is read in'	
     write(*,*) 'configured with nod2D=',nod2D,' nod3D=',nod3D
  end if

  ! ==============================
  ! Communication information
  ! ==============================
  file_name=trim(dist_mesh_dir)//'com_info'//trim(mype_string)//'.out'  
  fileID=10+mype  
  open(fileID, file=file_name)
  read(fileID,*)  n
  read(fileID,*) com_nod2D%rPEnum
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  read(fileID,*) com_nod2D%rPE
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1))
  read(fileID,*) com_nod2D%rptr
  allocate(com_nod2D%rlist(eDim_nod2D))
  read(fileID,*) com_nod2D%rlist

  read(fileID,*) com_nod2D%sPEnum
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  read(fileID,*) com_nod2D%sPE
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1))
  read(fileID,*) com_nod2D%sptr
  n=com_nod2D%sptr(com_nod2D%sPEnum+1)-1
  allocate(com_nod2D%slist(n))
  read(fileID,*) com_nod2D%slist

  read(fileID,*) com_nod3D%rPEnum
  allocate(com_nod3D%rPE(com_nod3D%rPEnum))
  read(fileID,*) com_nod3D%rPE
  allocate(com_nod3D%rptr(com_nod3D%rPEnum+1))
  read(fileID,*) com_nod3D%rptr
  allocate(com_nod3D%rlist(eDim_nod3D))
  read(fileID,*) com_nod3D%rlist

  read(fileID,*) com_nod3D%sPEnum
  allocate(com_nod3D%sPE(com_nod3D%sPEnum))
  read(fileID,*) com_nod3D%sPE
  allocate(com_nod3D%sptr(com_nod3D%sPEnum+1))
  read(fileID,*) com_nod3D%sptr
  n=com_nod3D%sptr(com_nod3D%sPEnum+1)-1
  allocate(com_nod3D%slist(n))
  read(fileID,*) com_nod3D%slist

  read(fileID, *) mapping(1:nod3D+nod2D)

  close(fileID)

  allocate(col_pos(myDim_nod3D+eDim_nod3D))
  col_pos=0

  if(mype==0) write(*,*) 'comm and mapping is read'
end subroutine  read_mesh
!
!=========================================================================
!
subroutine read_grid_slope  
  ! read sigma grid slope
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_config
  use g_parfe

  implicit none

  integer                   :: i, j
  real(kind=8)              :: temp(2)

  ! sigma grid slope
  allocate(grid_slope(2,max_num_layers-1,myDim_elem2d))
  grid_slope=0.0
  do i=1,myDim_elem2D
     mapping(myList_elem2D(i))=i
  end do
  open(13,file=trim(MeshPath)//'sigma_grid_slope_elem.out', status='old')
  do i=1,elem2d
     if(mapping(i)>0) then
        do j=1,max_num_layers-1
           read(13,*) temp(:)
           grid_slope(:,j,mapping(i))=temp(:)
        end do
     else
        do j=1,max_num_layers-1
           read(13,*) temp(:)
        end do
     end if
  end do
  close(13)
  mapping(1:elem2D)=0

end subroutine read_grid_slope
!
!=========================================================================
subroutine ocean_matrices_setup
  use o_param
  use o_mesh
  use g_parfe
  implicit none

  call set_coriolis_param      ! set f before assembling ssh matrix

  call sshstiff_matrix

#ifdef use_non_hydrostatic
  call nhpstiff_matrix  
#endif

  call tsstiff_construct

#ifndef use_tracer_gls
  call tsstiff_fill_tg
#endif

  call uvstiff_matrix

#ifndef use_non_hydrostatic
  call build_wpot_matrix
#endif

  if(mype==0) write(*,*) 'Ocean matrices have been set up'

end subroutine ocean_matrices_setup
!
!----------------------------------------------------------------------
!
subroutine uvstiff_matrix
  ! Stiffness matrix for 3D velocities (u and v entries are not coupled).
  ! mass matrix does not contain dt_inv 
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer                            :: k
  integer                            :: row, col, elnodes(4)
  integer                            :: i, q, elem, ind
  integer                            :: ipos, offset, is, ie
  real(kind=8)                       :: vol, inv20
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  uvstiff%dim=nod3D
  allocate(uvstiff%rowptr(myDim_nod3D+1))        
  uvstiff%rowptr(1)=1
  do k=1,myDim_nod3D                            
     uvstiff%rowptr(k+1)=uvstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  uvstiff%nza=uvstiff%rowptr(myDim_nod3D+1)-1    
  ! b)
  allocate(uvstiff%colind(uvstiff%nza))
  allocate(uvstiff%values(uvstiff%nza))
  uvstiff%values=0.0

  ! =================================            
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes), rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=uvstiff%nza
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  uvstiff%rowptr=uvstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  uvstiff%nza=sum(rpnza(1:npes))                  
  deallocate(rpnza, pnza)

  ! c)
  do k=1,myDim_nod3D                                 
     nini=uvstiff%rowptr(k)-uvstiff%rowptr(1)+1      
     nend=uvstiff%rowptr(k+1)-uvstiff%rowptr(1)      
     uvstiff%colind(nini:nend)= &                    
          nghbr_nod3D(k)%addresses
  end do

  ! ==================================
  ! addresses are now in local numbering. We need to make them global
  ! and then local contiguous
  ! (i) global natural: 
  do k=1,uvstiff%rowptr(myDim_nod3D+1)-uvstiff%rowptr(1)   
     uvstiff%colind(k)=myList_nod3D(uvstiff%colind(k))     
  end do
  ! (ii) global PE contiguous: 
  ! mapping(1:nod3D) is 3D mapping
  do k=1,uvstiff%rowptr(myDim_nod3D+1)-uvstiff%rowptr(1)   
     uvstiff%colind(k)=mapping(uvstiff%colind(k))          
  end do
  ! ==================================

  ! d) fill in
  allocate(uv_lump(myDim_nod3D+eDim_nod3D))           
  uv_lump=0.
  inv20=1.0_8/20.0_8

  do elem=1,myDim_elem3d                         
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*inv20
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                   
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=uvstiff%rowptr(row)-uvstiff%rowptr(1)  
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   uvstiff%values(ipos)=uvstiff%values(ipos)+vol
           if (q==i) then
              uvstiff%values(ipos)=uvstiff%values(ipos)+vol
           end if
        end do
        uv_lump(row)=uv_lump(row) + vol*5.0
     end do
  end do

  !boundary condition (no slip bc.)
  do i=1,myDim_nod3d                              
     ind = index_nod3D(i)
     if (ind==11 .or. ind==21 .or. ind==31) then
        is=uvstiff%rowptr(i)-uvstiff%rowptr(1)+1  
        ie=uvstiff%rowptr(i+1)-uvstiff%rowptr(1)   
        where(uvstiff%colind(is:ie)==mapping(i))   
           uvstiff%values(is:ie)=1.0
        elsewhere
           uvstiff%values(is:ie)=0.0
        end where
     end if
  end do

  call com_3D(uv_lump)  

end subroutine uvstiff_matrix
!
!========================================================================
!
subroutine sshstiff_matrix
  ! Stiffness matrix for ssh (surface pressure)
  use o_PARAM
  use o_array
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use g_config
  use G_parfe
  implicit none
  integer                           :: i, k, q, ipos, is, ie, row, col
  integer                           :: elem, elem2, elnodes2(3), mn(3)
  real(kind=8)                      :: vol, dx(3), dy(3), val, aux, tri_v(3)
  real(kind=8)                      :: cori_p, dparam,  beta, gamma
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  sshstiff%dim=nod2D
  allocate(sshstiff%rowptr(myDim_nod2D+1))         
  sshstiff%rowptr(1)=1              ! This will be updated as
  ! contiguous numbering is required 
  do k=1,myDim_nod2D                                
     sshstiff%rowptr(k+1)=sshstiff%rowptr(k)+nghbr_nod2D(k)%nmb
  end do
  sshstiff%nza=sshstiff%rowptr(myDim_nod2D+1)-1      

  ! b)
  allocate(sshstiff%colind(sshstiff%nza))
  allocate(sshstiff%values(sshstiff%nza))
  sshstiff%values=0.0

  ! =================================                  
  ! Now we need to exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza=0
  pnza(mype+1)=sshstiff%nza
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  sshstiff%rowptr=sshstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  sshstiff%nza=sum(rpnza(1:npes))                    
  deallocate(rpnza, pnza)

  ! c) FIND colind 
  do k=1,myDim_nod2D
     nini=sshstiff%rowptr(k)-sshstiff%rowptr(1)+1     
     nend=sshstiff%rowptr(k+1)-sshstiff%rowptr(1)     
     sshstiff%colind(nini:nend)= &
          nghbr_nod2D(k)%addresses
  end do
  ! ==================================
  ! colind are now in local numbering. We need to make them global
  ! and then global contiguous
  ! (i) global natural: 
  do k=1,sshstiff%rowptr(myDim_nod2D+1)-sshstiff%rowptr(1)  
     sshstiff%colind(k)=myList_nod2D(sshstiff%colind(k))   
  end do
  ! (ii) global PE contiguous: 
  ! mapping(1+nod3D:nod2D+nod3D) is 2D mapping
  do k=1,sshstiff%rowptr(myDim_nod2D+1)-sshstiff%rowptr(1)  
     sshstiff%colind(k)=mapping(nod3D+sshstiff%colind(k))   
  end do

  ! ==================================
  ! d) fill in

  ! g*dt*\int_{-H}^0(\nabla\phi_i * \nabla\phi_j) d\Omega
  do elem=1, myDim_elem3d                         
     elem2 = elem2D_corresp_to_elem3D(elem)
     elnodes2=elem2d_nodes(:,elem2)
     vol=voltetra(elem)
     dx=bafux_2d(:,elem2)
     dy=bafuy_2d(:,elem2)
     if(use_cori_semi) then
        !coriolis parameter
        cori_p=coriolis_param_elem2d(elem2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
#ifdef use_semiimplicit_scheme
     vol=vol*theta_ssh*theta_vel
#endif
     do i=1,3
        row=elnodes2(i)
	if (row>myDim_nod2D) cycle                  
        do q=1, nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q)) = q
        enddo
        do q=1,3
           col=elnodes2(q)
           if(use_cori_semi) then
              val=(dx(i)*dx(q)+dy(i)*dy(q))*vol*g*beta
              val=val+(dx(i)*dy(q)-dy(i)*dx(q))*vol*g*gamma
           else
              val = (dx(i)*dx(q) + dy(i)*dy(q))*vol*g*dt
           endif
           ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1) 
           sshstiff%values(ipos)=sshstiff%values(ipos)+val
        end do
     end do
  end do

  ! Mass matrix part
  do i=1, myDim_nod2d                                  
     is=sshstiff%rowptr(i)-sshstiff%rowptr(1)+1       
     ie=sshstiff%rowptr(i+1)-sshstiff%rowptr(1)       
     do q=1, nod_in_elem2D(i)%nmb
        elem2=nod_in_elem2D(i)%addresses(q)
        vol=voltriangle(elem2)
        elnodes2=elem2D_nodes(:, elem2)
        do k=1, 3
           col=elnodes2(k)
           aux=1./12.
           if (col==i) aux=1./6.
           row=mapping(nod3D+myList_nod2D(col))   
           where(sshstiff%colind(is:ie)==row)     
              sshstiff%values(is:ie)=sshstiff%values(is:ie)+aux*vol*dt_inv
           end where
        end do
     end do
  end do

  ! open boundary
#ifdef use_opbnd_tide
  if(nmbr_opbnd_tri==0) return                           
  if(trim(tide_opbnd_type)=='ssh') then
     do i=1,nmbr_opbnd_t2d
        row=opbnd_n2d(i)
        is=sshstiff%rowptr(row)-sshstiff%rowptr(1)+1     
        ie=sshstiff%rowptr(row+1)-sshstiff%rowptr(1)     
	q=mapping(nod3D+myList_nod2D(row))                 
        where(sshstiff%colind(is:ie)==q)                 
           sshstiff%values(is:ie)=1.0
        elsewhere
           sshstiff%values(is:ie)=0.0
        end where
     end do
  elseif(trim(tide_opbnd_type)=='Flather') then
     do elem2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(elem2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)    
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(elem2,4)
        tri_v=sqrt(g/opbnd_dep(mn)) 
        do i=1,3        
           row=elnodes2(i)
           if(row>myDim_nod2d) cycle
           is=sshstiff%rowptr(row)-sshstiff%rowptr(1)+1        
           ie=sshstiff%rowptr(row+1)-sshstiff%rowptr(1)         
           do k=1, 3
              col=elnodes2(k)
              aux=1./12.
              if (row==col) aux=1./6.
	      q=mapping(nod3D+myList_nod2D(col))           
              where(sshstiff%colind(is:ie)==q)               
                 sshstiff%values(is:ie)=sshstiff%values(is:ie)+aux*vol*tri_v(k)
              end where
           end do
        end do
     end do
  end if
#endif

end subroutine sshstiff_matrix
!
!============================================================================
!
#ifdef use_non_hydrostatic
subroutine nhpstiff_matrix
  ! Stiffness matrix for non-hydrostatic pressure
  use o_PARAM
  use o_array
  use o_MESH
  use o_ELEMENTS
  use o_MATRICES
  use g_config
  use g_parfe
  implicit none
  integer                           :: i, k, ipos, q
  integer                           :: row, col, elem, elnodes(4), elem2
  real(kind=8)                      :: vol, dx(4), dy(4), dz(4), val
  real(kind=8)                      :: cori_p, dparam, beta, gamma
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend
  !
  ! a)
  nhpstiff%dim=nod3D
  allocate(nhpstiff%rowptr(myDim_nod3D+1))          
  nhpstiff%rowptr(1)=1
  do k=1,myDim_nod3D                                 
     nhpstiff%rowptr(k+1)=nhpstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  nhpstiff%nza=nhpstiff%rowptr(myDim_nod3D+1)-1      
  ! b)
  allocate(nhpstiff%colind(nhpstiff%nza))
  allocate(nhpstiff%values(nhpstiff%nza))
  nhpstiff%values=0.0
  ! =================================                 
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes), rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=nhpstiff%nza
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  if (mype==0) then
     i=0
  else
     i=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  nhpstiff%rowptr=nhpstiff%rowptr+i       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  nhpstiff%nza=sum(rpnza(1:npes))
  deallocate(rpnza, pnza)                           

  ! c)
  do k=1,myDim_nod3D                                
     nini=nhpstiff%rowptr(k)-nhpstiff%rowptr(1)+1     
     nend=nhpstiff%rowptr(k+1)-nhpstiff%rowptr(1)     
     nhpstiff%colind(nini:nend)= &                    
          nghbr_nod3D(k)%addresses
  end do
  ! ==================================
  ! addresses are now in local numbering. We need to make them global
  ! and then local contiguous
  ! (i) global natural: 
  do k=1,nhpstiff%rowptr(myDim_nod3D+1)-nhpstiff%rowptr(1) 
     nhpstiff%colind(k)=myList_nod3D(nhpstiff%colind(k))   
  end do
  ! (ii) global PE contiguous: 
  ! mapping(1:nod3D) is 3D mapping
  do k=1,nhpstiff%rowptr(myDim_nod3D+1)-nhpstiff%rowptr(1) 
     nhpstiff%colind(k)=mapping(nhpstiff%colind(k))        
  end do
  ! ==================================


  ! d) fill in 
  do elem=1,myDim_elem3d                  
     elnodes=elem3d_nodes(:,elem)
     vol=voltetra(elem)
     dx=bafux_3d(:,elem)
     dy=bafuy_3d(:,elem)
     dz=bafuz_3d(:,elem)
     if(use_cori_semi) then
        elem2=elem2d_corresp_to_elem3d(elem)
        !coriolis parameter
        cori_p=coriolis_param_elem2d(elem2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
     do i=1,4
        row=elnodes(i)
	if (row>myDim_nod3D) cycle                        
        do q=1, nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q)) = q
        enddo
        do q=1,4
           col=elnodes(q)
           if(use_cori_semi) then
              val=(dx(i)*dx(q)+dy(i)*dy(q))*beta
              val=val+(dx(i)*dy(q)-dy(i)*dx(q))*gamma
              val=(val+dz(i)*dz(q)*dt)*vol
           else
              val=(dx(i)*dx(q)+dy(i)*dy(q)+dz(i)*dz(q))*vol*dt
           endif
           ipos=nhpstiff%rowptr(row)+col_pos(col)-nhpstiff%rowptr(1)  
           nhpstiff%values(ipos)=nhpstiff%values(ipos)+val
        end do
     end do
  end do
end subroutine nhpstiff_matrix
#endif
!
!============================================================================
!
#ifndef use_non_hydrostatic
subroutine build_wpot_matrix
  ! w potential matrix
  ! Taking into account the property of vertical dirivatives in tetrahedra:
  ! only a three-diagonal matrix is required.
  use o_param
  use o_mesh
  use o_elements
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer            :: n2, nlay, k, i, q, p
  integer            :: row, elem, elnodes(4)
  integer            :: nodup, nodlo
  real(kind=8)       :: vol
  real(kind=8)       :: a(max_num_layers), b(max_num_layers)
  real(kind=8)       :: c(max_num_layers)

  allocate(wpot_matrix(3,max_num_layers,myDim_nod2d))
  wpot_matrix=0.0

  do n2=1,myDim_nod2d
     nlay=num_layers_below_nod2d(n2)+1

     ! assemble three diagonal matrix
     do k=1,nlay
        a(k)=0.
        b(k)=0.
        c(k)=0.

        row=nod3d_below_nod2d(k,n2)
        if(k>1) nodup=nod3d_below_nod2d(k-1,n2)
        if(k<nlay) nodlo=nod3d_below_nod2d(k+1,n2)
        do i=1,nod_in_elem3D(row)%nmb
           elem=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,elem)  
           vol=voltetra(elem)

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    a(k)=a(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * vol
                    exit
                 end if
              end do
           end if

           !second entry
           b(k)=b(k) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * vol

           !third entry
           if(k<nlay) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    c(k)=c(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem)  * vol
                    exit
                 end if
              end do
           end if

        end do  !i
     end do  !k

     wpot_matrix(1,1:nlay,n2)=a(1:nlay)      
     wpot_matrix(2,1:nlay,n2)=b(1:nlay)    
     wpot_matrix(3,1:nlay,n2)=c(1:nlay)    

  end do  !m

end subroutine build_wpot_matrix
#endif
!
!=========================================================================
!
subroutine tsstiff_construct
  !T/S stiffness matrix
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use g_config
  use g_PARFE
  implicit none
  !
  integer                           :: k
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend

  ! a)
  tsstiff%dim=nod3D
  allocate(tsstiff%rowptr(myDim_nod3D+1))     
  tsstiff%rowptr(1)=1
  do k=1,myDim_nod3D                          
     tsstiff%rowptr(k+1)=tsstiff%rowptr(k)+nghbr_nod3D(k)%nmb
  end do
  tsstiff%nza=tsstiff%rowptr(myDim_nod3D+1)-1   

  ! b)
  allocate(tsstiff%colind(tsstiff%nza))
  allocate(tsstiff%values(tsstiff%nza))
  tsstiff%values=0.0
  ! =================================           
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=tsstiff%nza

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  if (mype==0) then
     k=0
  else
     k=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  tsstiff%rowptr=tsstiff%rowptr+k       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  tsstiff%nza=sum(rpnza(1:npes))              
  deallocate(rpnza,pnza)
  ! c)
  ! ===== FIND irowind ======
  do k=1,myDim_nod3D                            
     nini=tsstiff%rowptr(k)-tsstiff%rowptr(1)+1 
     nend=tsstiff%rowptr(k+1)-tsstiff%rowptr(1) 
     tsstiff%colind(nini:nend)= &                
          nghbr_nod3D(k)%addresses
  end do
  ! ==================================
  ! addresses are now in local numbering. We need to make them global
  ! and then local contiguous
  ! (i) global natural: 
  do k=1,tsstiff%rowptr(myDim_nod3D+1)-tsstiff%rowptr(1)
     tsstiff%colind(k)=myList_nod3D(tsstiff%colind(k))   
  end do
  ! (ii) global PE contiguous: 
  ! mapping(1:nod3D) is 3D mapping
  do k=1,tsstiff%rowptr(myDim_nod3D+1)-tsstiff%rowptr(1) 
     tsstiff%colind(k)=mapping(tsstiff%colind(k))        
  end do
  ! ==================================

end subroutine tsstiff_construct
!
!=========================================================================
!
#ifdef use_tracer_gls
subroutine tsstiff_fill_gls
  ! Fill in T/S stiffness matrix for GLS scheme
  ! 1) Fills in  advective, diffusive and time contributions
  ! 2) Computes RHSs for temperature and salinity
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_forcing_arrays
  use g_PARFE
  implicit none 

  integer                       :: k, i, q, elem, elem2, elem_type
  integer                       :: elnodes(4), is, ie, lay
  integer                       :: row, col, col2, ipos, offset, ntr
  real(kind=8)                  :: uvel(4), vvel(4), usum, vsum
  real(kind=8)                  :: u2d, v2d, um, vm, wm
  real(kind=8)                  :: dx(4), dy(4), dz(4), tra_elem(4,num_tracer)
  real(kind=8)                  :: auxf, advection, diffusion
  real(kind=8)                  :: vnabla_t, ept, ept_dt_inv
  real(kind=8)                  :: v, vtri, hh, hv, entries(4), entries_t(4)
  real(kind=8)                  :: Kh, Kv_el, r_coeff
  real(kind=8)                  :: inv2, inv4, inv5, inv20
  real(kind=8)                  :: fc, aux1, aux2, dif(4)
  real(kind=8)                  :: rotate_coe, temp_coe, temp_coe2
  real(kind=8)                  :: swr_conv
  !
  !variables used for Redi/GM scheme 
  real(kind=8)                  :: K_GM, Kh_diag, Kh_other
  real(kind=8)		        :: S(3), fcn1,fcn2, lambd, Z_mean, depth_scale
  real(kind=8)        		:: S_d, c_speed
  data S_d, c_speed /1.0e-3, 2.0/
  !
  real(kind=8)                  :: dparam, beta, gamma
#ifndef use_non_hydrostatic
  integer                       :: elnodes2(3)
  real(kind=8)                  :: vc(3)
#else
  integer                       :: elnodes23(4)
  real(kind=8)                  :: vc(4), wsum, wvel(4), w2d
#endif
#ifdef use_fullfreesurf
  integer                       :: n_el
  real(kind=8)                  :: wmove(4), wmove_sum
  real(kind=8)                  :: v_new, auxf_new, entries_rhs(4)
  logical                       :: flag_move=.false.
#endif

  inv2=0.5_8
  inv4=0.25_8
  inv5=0.2_8
  inv20=0.05_8
  
  tracer_rhs=0.0
  tsstiff%values(:)=0.0
 
  do elem=1, myDim_elem3d     
     elnodes=elem3D_nodes(:,elem)
     elem2 = elem2D_corresp_to_elem3D(elem)
     elem_type=grid_type_elem2d(elem2)
#ifndef use_non_hydrostatic
     elnodes2=elem2D_nodes(:,elem2)
#else
     elnodes23=nod2d_corresp_to_nod3d(elnodes)   
#endif
     fc=coriolis_param_elem2d(elem2)

     ! elementwise lateral(neutral) and vertical(dianeutral) diffusivity 
     Kh=Kh0
     if(scale_mixing_h) Kh=Kh*(voltriangle(elem2)/scalevol)**(1.0/real(scale_mixing_type))
     Kv_el=Kv0+sum(Kv(elnodes,1))/4.0

     !-------------------------------------------------------------------------------
     if (Redi_GM .and. elem_type==0) then 

        !GM diffusivity
        K_GM  = Kh*ratio_K_GM   
           
        !neutral slope
        if(nslope_version==1) then
           lay=elem3d_layer(elem)
           S = neutral_slope(:,lay,elem2)   ! S(1:3): Sx,Sy and |S|
        else
           S = neutral_slope_elem(:,elem) 
        end if

        !prepare for tapering
        !define 2 functions fcn1 and fcn2, which are required for tapering
        fcn1=1.0_8
        fcn2=1.0_8

        ! fcn1, hyperbolic tangent, used for steep slope region
        if(ODM95) fcn1 = 0.5_8*(1.0_8 + tanh((Slope_critical - S(3))/S_d))

        !we need to check if the element is near the surface 
        !If yes, then we need the tapering function fcn2, a sine function of depth.
        if(LDD97) then
           !the first baroclinic Rossby radius
           lambd = c_speed/abs(fc)

           !limit lambda [following Large et al(1997)] to handle singularity
           !at the equator
           if (lambd < 15000.) lambd = 15000.
           if (lambd > 100000.) lambd = 100000.

           !critical depth, above which sine tapering is necessary.
           depth_scale = 100.0 + lambd*S(3)
           !here 100.0m depth is assumed to be turbulent and does not need Redi/GM
           !a well defined turbulent layer should be updated later

           !the mean depth of this element
           Z_mean = abs(sum(coord_nod3D(3,elnodes)))/4.0

           !if in the surface layer we add tapering function f2
           if (Z_mean < depth_scale)  then
              fcn2 = 0.5*(1.0 + sin(pi*Z_mean/depth_scale - pi/2.0))
           end if
        end if
      
        ! apply tapering
        ! For steep slope region:
        ! a) no taper applied to diagonal piece of horizontal neutral operator
        ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
        !    horizontal operator and to diagonal and off-diagonal piece of vertical
        !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
        !    to a horizontal-vertical manner in regions of steep neutral slopes.
        ! c) Exponential taper applied to GM operator.
        ! For surface layer with small slope:
        ! a) sine taper applied to both neutral operator and GM operator, except the
        !    diagonal piece of the horizontal diffusion.
        ! In one word, here we use ldd97, but always keep the diagonal part of the 
        ! horizontal diffusion following the suggestion of Griffies (2004).
        
        ! diffusion part:
        Kh_diag = Kh
        Kh_other = Kh*fcn1*fcn2
        ! skewion part:
        K_GM = K_GM*fcn1*fcn2
        
     end if  	!Redi_GM:  tapered neutral diffusivity computed
     !------------------------------------------------------------------------------

     vtri=voltriangle(elem2)
     v=voltetra(elem)
     hh=sqrt(vtri)
     hv=v/vtri     

     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)

     uvel=uf0(elnodes)         !u*
     vvel=uf0(elnodes+myDim_nod3d+eDim_nod3D)   !v*      
     usum=sum(uvel)
     vsum=sum(vvel)

     tra_elem=tracer(elnodes,:)

     if(use_cori_semi) then
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
     endif

     ! shortwave penetration
#ifdef use_sw_pene
     swr_conv=sum(dz*sw_3d(elnodes))*v
#endif

#ifndef use_non_hydrostatic
#ifdef use_semiimplicit_scheme
     vc=theta_ssh*ssh(elnodes2)-(gamma_stab-1.0+theta_ssh)*ssh0(elnodes2)	
#else
     vc=ssh(elnodes2)-gamma_stab*ssh0(elnodes2)
#endif
     aux1=sum(vc*bafux_2d(:,elem2))*g
     aux2=sum(vc*bafuy_2d(:,elem2))*g
     if(use_cori_semi) then
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d	
     wm=sum(dz*w(elnodes))
#else
     !non-hydrostatic case
     wvel=uf0(elnodes+2*(myDim_nod3D+eDim_nod3D))          
     wsum=sum(wvel)
     vc=g*ssh(elnodes23)+nhp(elnodes) - &
          gamma_stab*g*ssh0(elnodes23)-gamma_stab_nh*nhp0(elnodes)
     aux1=sum(vc*dx)
     aux2=sum(vc*dy)
     if(use_cori_semi) then  
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     w2d=sum(vc*dz)*dt
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d
     wm=inv4*wsum-w2d
#endif

#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     flag_move=.false.
     v_new=v
     if(n_el/=0) then
        flag_move=.true.
        v_new=voltetra_new(n_el)
        wmove=0.
        do i=1, 4
           row=elnodes(i)
           if(row<=nod2d) wmove(i)=-(ssh(row)-ssh0(row))*dt_inv
        end do
        wmove_sum=sum(wmove)  
     end if
#endif

     ept=2.0*Kh/(hh*hh)+Kv_el/(hv*hv)+abs(um)/hh+abs(vm)/hh+abs(wm)/hv + 0.01*dt_inv 
#ifdef use_fullfreesurf
     ept=ept+abs(wmove_sum)/hv 
#endif
     ept=v/ept 
     ept_dt_inv=ept*dt_inv
     auxf=v*dt_inv*inv20
#ifdef use_fullfreesurf
     auxf_new=v_new*dt_inv*inv20
#endif

     !for along sigma mixing
     if (elem_type==1) then
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if


     entries=0.
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
        if(row>myDim_nod3D) cycle                      
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        vnabla_t=um*dx(i)+vm*dy(i)+wm*dz(i) 

#ifndef use_fullfreesurf  
	! linear freesurface
        do q=1,4          ! all columns 
           if(elem_type==1) then !sigma diff.
              diffusion= &
                   (dx(i)*dx(q)*(1.0+S(2)*S(2))+dy(i)*dy(q)* &
                   (1.0+S(1)*S(1)))*Kh*rotate_coe + &
                   dz(i)*dz(q)*(Kv_el+Kh*S(3)*S(3)*rotate_coe) + &  
                   (S(1)*dx(i)+S(2)*dy(i))*dz(q)*Kh*rotate_coe - &  
                   (dx(i)*dy(q)+dy(i)*dx(q))*S(1)*S(2)*Kh*rotate_coe + &
                   (S(1)*dx(q)+S(2)*dy(q))*dz(i)*Kh*rotate_coe
           else
	      if (Redi_GM) then
                 diffusion=Kh_diag*(dx(i)*dx(q)+dy(i)*dy(q)) &
                      + (Kh_other-K_GM)*(S(1)*dx(i)+S(2)*dy(i))*dz(q) &
                      + (Kv_el+Kh_other*S(3)*S(3))*dz(i)*dz(q) &    
                      + (Kh_other+K_GM)*(S(1)*dx(q)+S(2)*dy(q))*dz(i) 
              else
                 diffusion=Kh*(dx(i)*dx(q)+dy(i)*dy(q))+Kv_el*dz(i)*dz(q)
	      end if
           end if
#ifndef use_non_hydrostatic
           advection=((usum+uvel(i))*inv5-u2d)*dx(q) + &
                ((vsum+vvel(i))*inv5-v2d)*dy(q) + & 
                wm*dz(q) 
#else
           advection=((usum+uvel(i))*inv5-u2d)*dx(q) + &
                ((vsum+vvel(i))*inv5-v2d)*dy(q) + & 
                ((wsum+wvel(i))*inv5-w2d)*dz(q) 
#endif  
	   entries(q)=(diffusion+advection*inv4)*v              ! diff and adv

	   advection= um*dx(q) + vm*dy(q) + wm*dz(q) 
           entries(q)=entries(q)+ ept*vnabla_t*advection       	! stabilization, upwind 

           col=elnodes(q)
           do ntr=1,num_tracer
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)-entries(q)*tra_elem(q,ntr)   ! fill tracer_rhs
           end do

           entries_t(q)=auxf                            	! mass matrix
           entries_t(q)=entries_t(q)+inv4*vnabla_t*ept_dt_inv 	! stabilization, t
        end do

        entries_t(i)=entries_t(i)+auxf                    	! completes mass matrix 
        entries=inv2*entries+entries_t
        ! put the entries to the appropriate place
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)     
        do q=1,4        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+entries(q)
        end do

#else   
	! full free surface, mass matrix updated
        do q=1,4        ! all columns 
           if(elem_type==1) then !sigma diff.
              diffusion= &
                   (dx(i)*dx(q)*(1.0+S(2)*S(2))+dy(i)*dy(q)* &
                   (1.0+S(1)*S(1)))*Kh*rotate_coe + &
                   dz(i)*dz(q)*(Kv_el+Kh*S(3)*S(3)*rotate_coe) + &  
                   (S(1)*dx(i)+S(2)*dy(i))*dz(q)*Kh*rotate_coe - &  
                   (dx(i)*dy(q)+dy(i)*dx(q))*S(1)*S(2)*Kh*rotate_coe + &
                   (S(1)*dx(q)+S(2)*dy(q))*dz(i)*Kh*rotate_coe
           else
	      if (Redi_GM) then
                 diffusion=Kh_diag*(dx(i)*dx(q)+dy(i)*dy(q)) &
                      + (Kh_other-K_GM)*(S(1)*dx(i)+S(2)*dy(i))*dz(q) &
                      + (Kv_el+Kh_other*S(3)*S(3))*dz(i)*dz(q) &    
                      + (Kh_other+K_GM)*(S(1)*dx(q)+S(2)*dy(q))*dz(i) 
	      else
                 diffusion=Kh*(dx(i)*dx(q)+dy(i)*dy(q))+Kv_el*dz(i)*dz(q)
	      end if
           end if

#ifndef use_non_hydrostatic
           advection=(dx(i)*(usum+uvel(q))+dy(i)*(vsum+vvel(q)))*inv20 + &
		dz(i)*wm*inv4 - (dx(i)*u2d+dy(i)*v2d)*inv4 
           if(flag_move) advection=advection+dz(i)*(wmove_sum+wmove(q))*inv20
#else
           advection=(dx(i)*(usum+uvel(q))+dy(i)*(vsum+vvel(q))+ &
		dz(i)*(wsum+wvel(q)))*inv20 - &
		(dx(i)*u2d+dy(i)*v2d+dz(i)*w2d)*inv4 
           if(flag_move) advection=advection+dz(i)*(wmove_sum+wmove(q))*inv20
#endif  
	   entries(q)=(diffusion-advection)*v               	! diff and adv

	   advection= um*dx(q) + vm*dy(q) + wm*dz(q)
           if(flag_move) advection=advection+wmove_sum*inv4*dz(q)
           entries(q)=entries(q)+ ept*vnabla_t*advection     	! stabilization, upwind 

           entries_t(q) = auxf_new                            
           entries_rhs(q) = entries(q)+(auxf_new-auxf)     
           if(i==q) then
              entries_t(q)=entries_t(q)+auxf_new              	! mass matrix
              entries_rhs(q)=entries_rhs(q)+(auxf_new-auxf)  	! rhs entries
           endif

           col=elnodes(q)
           do ntr=1,num_tracer
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)-entries_rhs(q)*tra_elem(q,ntr)   ! fill tracer_rhs
           end do

           entries_t(q)=entries_t(q)+inv4*vnabla_t*ept_dt_inv 	! stabilization, t
        end do

        entries=inv2*entries+entries_t
        ! put the entries to the appropriate place
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)         
        do q=1,4        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+entries(q)
        end do

#endif

        ! other rhs terms

        ! in case of considering shortwave penetration into the ocean
#ifdef use_sw_pene
        tracer_rhs(row,1)=tracer_rhs(row,1)+swr_conv*inv4 
#endif

        if(buffer_zone) then
           aux1=tracer_restore_coeff(row)*v*inv20
           do ntr=1,num_tracer
              dif=tracer0(elnodes,ntr)-tra_elem(:,ntr)
              tracer_rhs(row,ntr)=tracer_rhs(row,ntr)+aux1*(sum(dif)+dif(i))
           end do
        end if

     end do  !rows to which this element can contribute 
  end do  ! 3d element

end subroutine tsstiff_fill_gls
#endif
!
!--------------------------------------------------------------------------
!
#ifndef use_tracer_gls
subroutine tsstiff_fill_tg
  ! Fill in T/S stiffness matrix for TG scheme
  ! mass matrix contains dt_inv
  use o_MATRICES
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                            :: m, i, q, row, col, ipos, offset
  integer                            :: elem, elnodes(4)
  real(kind=8)                       :: vol, inv20
  !
  inv20=0.05_8
  tsstiff%values=0.0
  allocate(ts_lump(myDim_nod3d+eDim_nod3D))      
  ts_lump=0.
  !
  !
  do elem=1,myDim_elem3d                                                 
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*dt_inv*inv20    ! contains dt_inv
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle               
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)   
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   tsstiff%values(ipos)=tsstiff%values(ipos)+vol
           if (q==i) then
              tsstiff%values(ipos)=tsstiff%values(ipos)+vol
           end if
        end do
        ts_lump(row)=ts_lump(row) + vol*5.0
     end do
  end do

  call com_3d(ts_lump)   !required for the FCT scheme
  
end subroutine tsstiff_fill_tg
#endif
!
!--------------------------------------------------------------------------
!
#ifdef use_fullfreesurf
subroutine update_matrices
  use o_param
  use o_array
  use o_ELEMENTS
  use o_MESH
  use o_matrices
  use g_config
  use g_parfe
  implicit none
  !
  integer          :: m, el, n_el, el2, elnodes(4), elnodes2(3), n2
  integer          :: i, q, p, row, k, col, ind, offset, ipos
  integer          :: nodup, nodlo
  real(kind=8)     :: vol, vol_new, dx2(3), dy2(3), dz(4), dz_new(4)
  real(kind=8)     :: aux, aux2, inv20, a, b, c
  real(kind=8)     :: cori_p, dparam,  beta, gamma
#ifdef use_non_hydrostatic
  real(kind=8)     :: dx(4), dy(4), dx_new(4), dy_new(4)
#endif

  inv20=0.05_8

  do el=1,myDim_elem3d                           
     n_el=map_elem(el)
     if(n_el==0) cycle 
     elnodes=elem3D_nodes(:,el)
     el2=elem2D_corresp_to_elem3D(el)
     elnodes2=elem2d_nodes(:,el2)
     vol=voltetra(el)
     vol_new=voltetra_new(n_el)
     dx2=bafux_2d(:,el2)
     dy2=bafuy_2d(:,el2)
     dz=bafuz_3d(:,el)
     dz_new=bafuz_3d_new(:,n_el)
#ifdef use_non_hydrostatic
     dx=bafux_3d(:,el)
     dy=bafuy_3d(:,el)
     dx_new=bafux_3d_new(:,n_el)
     dy_new=bafuy_3d_new(:,n_el)
#endif
     !
     !
     ! update uv matrix
     aux=(vol_new-vol)*inv20
     do i=1,4             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                      
        ind=index_nod3d(row)
        if(mod(ind,10)==1) cycle
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=uvstiff%rowptr(row)-uvstiff%rowptr(1)   
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   uvstiff%values(ipos)=uvstiff%values(ipos)+aux
           if (q==i) then
              uvstiff%values(ipos)=uvstiff%values(ipos)+aux
           end if
        end do
     end do
     !
     !
     if(biharmonic_visc .or. lump_uv_matrix) then
        aux=(vol_new-vol)*inv20
        do i=1,4  
           row=elnodes(i)
	   if(row>myDim_nod3D) cycle                  
           ind=index_nod3d(row)
           if(mod(ind,10)==1) cycle
           uv_lump(row)=uv_lump(row) + aux*5.0
        end do
     endif
     !
     !
     ! update ssh matrix
     if(use_cori_semi) then
        !coriolis parameter
        cori_p=coriolis_param_elem2d(el2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     endif
     do i=1,3
        row=elnodes2(i)
	if(row>myDim_nod2D) cycle                         
        do q=1, nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q)) = q
        enddo
        do q=1,3
           col=elnodes2(q)
           if(use_cori_semi) then
              aux=(dx2(i)*dx2(q)+dy2(i)*dy2(q))*(vol_new-vol)*g*beta
              aux=aux+(dx2(i)*dy2(q)-dy2(i)*dx2(q))*(vol_new-vol)*g*gamma
           else
              aux = (dx2(i)*dx2(q) + dy2(i)*dy2(q))*(vol_new-vol)*g*dt
           endif
           ipos=sshstiff%rowptr(row)+col_pos(col)-sshstiff%rowptr(1)   
           sshstiff%values(ipos)=sshstiff%values(ipos)+aux
        end do
     end do
     !
     !
#ifdef use_non_hydrostatic
     ! update nhp matrix
     do i=1,4
        row=elnodes(i)
	if(row>myDim_nod3D) cycle                          
        do q=1, nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q)) = q
        enddo
        do q=1,4
           col=elnodes(q)
           if(use_cori_semi) then
              aux=(dx_new(i)*dx_new(q)+dy_new(i)*dy_new(q))*beta
              aux=aux+(dx_new(i)*dy_new(q)-dy_new(i)*dx_new(q))*gamma
              aux=(aux+dz_new(i)*dz_new(q)*dt)*vol_new
              aux2=(dx(i)*dx(q)+dy(i)*dy(q))*beta
              aux2=aux2+(dx(i)*dy(q)-dy(i)*dx(q))*gamma
              aux=aux-(aux2+dz(i)*dz(q)*dt)*vol
           else
              aux=(dx_new(i)*dx_new(q)+dy_new(i)*dy_new(q)+ &
                   dz_new(i)*dz_new(q))*vol_new*dt
              aux=aux-(dx(i)*dx(q)+dy(i)*dy(q)+dz(i)*dz(q))*vol*dt
           endif
           ipos=nhpstiff%rowptr(row)+col_pos(col)-nhpstiff%rowptr(1)  
           nhpstiff%values(ipos)=nhpstiff%values(ipos)+aux
        end do
     end do
#endif
     !
     !
#ifndef use_tracer_gls
     ! update ts matrix 
     aux=(vol_new-vol)*dt_inv*inv20
     do i=1,4             ! all rows into which the element el could contribute
        row=elnodes(i)
        if(row>myDim_nod3D) cycle                             
        do q=1,nghbr_nod3D(row)%nmb
           col_pos(nghbr_nod3D(row)%addresses(q))=q
        end do
        offset=tsstiff%rowptr(row)-tsstiff%rowptr(1)            
        do q=1,4          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           tsstiff%values(ipos)=tsstiff%values(ipos)+aux
           if (q==i) then
              tsstiff%values(ipos)=tsstiff%values(ipos)+aux
           end if
        end do
        ts_lump(row)=ts_lump(row) + aux*5.0
     end do
#endif
  end do
  !
  !
#ifndef use_non_hydrostatic
  ! update w potential matrix
  
  do n2=1,myDim_nod2d
     do k=1,2   !only the first layer is moved
        a=0.
        b=0.
        c=0.
        row=nod3d_below_nod2d(k,n2)
        if(k==2) nodup=nod3d_below_nod2d(k-1,n2)
        if(k==1) nodlo=nod3d_below_nod2d(k+1,n2)

        do i=1,nod_in_elem3D(row)%nmb
           el=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,el)  
           n_el=map_elem(el)
           if(n_el==0) cycle
           vol=voltetra(el)
           vol_new=voltetra_new(n_el)
           dz=bafuz_3d(:,el)
           dz_new=bafuz_3d_new(:,n_el)

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do         

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    a=a + dz_new(p)*dz_new(q)*vol_new - dz(p)*dz(q)*vol
                    exit
                 end if
              end do
           end if
           
           !second entry
           b=b + dz_new(p)*dz_new(p)*vol_new - dz(p)*dz(p)*vol
        
           !third entry
           if(k<2) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    c=c + dz_new(p)*dz_new(q)*vol_new - dz(p)*dz(q)*vol
                    exit
                 end if
              end do
           end if
        end do  !i

        wpot_matrix(1,k,n2)= wpot_matrix(1,k,n2)+a
        wpot_matrix(2,k,n2)= wpot_matrix(2,k,n2)+b
        wpot_matrix(3,k,n2)= wpot_matrix(3,k,n2)+c

     end do  !k
  end do  !n2
#endif

end subroutine update_matrices
#endif
!
!------------------------------------------------------------------
!
subroutine uv_sfc_bott_bc
  ! momentum surface and bottom boundary conditions
  ! This routine is only required when using implicit vertical mixing scheme.
  ! In case of the explicit scheme, the assembling is done directly within
  ! the velocity_rhs routine. This is different from the tracer bc treatment.
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none

  integer      :: elem2, q, row, elnodes2(3), elnodes2_3d(3)
  real(kind=8) :: vol, usum, vsum, aux, uvel(3), vvel(3)
  real(kind=8) :: inv3, inv12

  inv3=1./3.
  inv12=1./12.

  uv_sfc_force=0.
  uv_bott_force=0.

  do elem2=1, myDim_elem2d    
     vol=voltriangle(elem2)
     elnodes2=elem2D_nodes(:,elem2)

     !surface             
     !Wind/ice stress, or stress from ice-cavity    
#ifdef use_cavity
     if(all(cavity_flag_nod2d(elnodes2)==0)) then   
#endif
        usum=sum(stress_x(elnodes2))
        vsum=sum(stress_y(elnodes2))
        aux=rho0r*inv12*vol
        do q=1,3
           row=elnodes2(q)
           uv_sfc_force(row,1)=uv_sfc_force(row,1) + aux*(usum+stress_x(row)) 
           uv_sfc_force(row,2)=uv_sfc_force(row,2) + aux*(vsum+stress_y(row)) 
        end do
#ifdef use_cavity
     else
        uvel=uf(nod3D_below_nod2D(1,elnodes2))      
        vvel=uf(nod3D_below_nod2D(1,elnodes2)+ToDim_nod3d)   
        usum=sum(uvel)
        vsum=sum(vvel)
        aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*vol
        do q=1,3
           row=elnodes2(q)
           uv_sfc_force(row,1)=uv_sfc_force(row,1) + aux*(usum+uvel(q)) 
           uv_sfc_force(row,2)=uv_sfc_force(row,2) + aux*(vsum+vvel(q)) 
        end do
     end if
#endif

     !bottom 
     !Bottom drag contribution  -C_d U|U| 
     elnodes2_3d=bt_nds(elnodes2)        
     uvel=uf(elnodes2_3d)
     vvel=uf(elnodes2_3d+ToDim_nod3D)                           
     usum=sum(uvel)
     vsum=sum(vvel)
     aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*vol
     do q=1,3
        row=elnodes2(q)
        uv_bott_force(row,1)=uv_bott_force(row,1) + aux*(usum+uvel(q)) 
        uv_bott_force(row,2)=uv_bott_force(row,2) + aux*(vsum+vvel(q))        
     end do
  end do

  do q=1, myDim_nod2d  !eDim_nod2d is not required 
     if(mod(index_nod3d(q),10)==1) then
     	uv_sfc_force(q,1)=0.
	uv_sfc_force(q,2)=0.
	uv_bott_force(q,1)=0.
	uv_bott_force(q,1)=0.
     endif
     if(mod(index_nod3d(bt_nds(q)),10)==1) then
        uv_bott_force(q,1)=0.
	uv_bott_force(q,2)=0.
     endif
  end do

end subroutine uv_sfc_bott_bc
!
!----------------------------------------------------------------------------
!
subroutine velocity_rhs
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer                        :: m, i, elem, elem2, elnodes(4), elnodes23(4)
  integer                        :: q, row, row2, elnodes2(3), ind, lay
  integer                        :: n3, elem_type
  real(kind=8)                   :: fc, aux, adv, coe, Ah, Av_el, d1, d2
  real(kind=8)                   :: nhp_el(4)
  real(kind=8)                   :: px_ssh, py_ssh, uvel(4), vvel(4)
  real(kind=8)                   :: usum, vsum, um, vm, wm
  real(kind=8)                   :: pxsum ,pysum
  real(kind=8)                   :: udx, udy, udz, vdx, vdy, vdz
  real(kind=8)                   :: dx(4), dy(4), dz(4), v, vtr, val
  real(kind=8)                   :: inv5, inv4, inv3, inv12, inv20
  real(kind=8)                   :: urhs_temp, vrhs_temp, cori_u_tmp, cori_v_tmp
  real(kind=8)                   :: S(3), rotate_coe, temp_coe, temp_coe2
  real(kind=8)                   :: dparam, beta, gamma
#ifdef use_non_hydrostatic
  integer                        :: row3
  real(kind=8)                   :: wvel(4), wsum, npz, npx, npy
  real(kind=8)                   :: wdx, wdy, wdz, wrhs_temp
#endif
#ifdef use_fullfreesurf
  integer                        :: n_el
  real(kind=8)                   :: coe_diff, uusum, uvsum, vvsum, wmove(4)
  real(kind=8)                   :: wmove_sum, uwmove_sum, vwmove_sum
  logical                        :: flag_move=.false.
#ifdef use_non_hydrostatic
  real(kind=8)                   :: wusum, wvsum, wwsum, wwmove_sum
#endif
#endif

  inv5=0.2_8
  inv4=0.25_8
  inv3=1.0_8/3.0_8
  inv12=1.0_8/12.0_8
  inv20=0.05_8
  n3=ToDim_nod3d
  do row=1, myDim_nod3d                   
     row2=row+n3                          
     uv_rhs(row)=0.
     uv_rhs(row2)=0.
#ifdef use_non_hydrostatic
     row3=row2+n3                       
     uv_rhs(row3)=0.
#endif
  enddo

  !3D contributions	 
  do elem=1, myDim_elem3d       
     elnodes=elem3D_nodes(:,elem)
     elem2 = elem2D_corresp_to_elem3D(elem)
     elnodes2=elem2D_nodes(:,elem2)
     elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)
     v=voltetra(elem)
     elem_type=grid_type_elem2d(elem2)

     !velocity
     uvel=uf(elnodes)
     vvel=uf(elnodes+n3)                  
     usum=sum(uvel)         
     vsum=sum(vvel)
     um=usum*inv4
     vm=vsum*inv4
#ifndef use_non_hydrostatic
     wm=sum(w(elnodes)*dz)
#else
     wvel=uf(elnodes+2*n3)                
     wsum=sum(wvel)
     wm=wsum*inv4
#endif

     !coriolis parameter
     fc=coriolis_param_elem2D(elem2)
     !terms for semi-implicit schemes
     if(use_cori_semi) then
        val=v
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
     else
        val=v*dt
     endif

     !hp
     if(elem_type==0) then
        pxsum=sum(dx*hpressure(elnodes))
        pysum=sum(dy*hpressure(elnodes)) 
     else 
        lay=elem3d_layer(elem) 
        pxsum=PGF(1,lay,elem2)                
        pysum=PGF(2,lay,elem2)
     end if

     !ssh/non-hydrostatic p
     nhp_el=g*ssh(elnodes23)
     px_ssh=sum(dx*nhp_el)
     py_ssh=sum(dy*nhp_el)
#ifdef use_non_hydrostatic
     nhp_el=nhp(elnodes)
     npx=sum(dx*nhp_el)
     npy=sum(dy*nhp_el)
     npz=sum(dz*nhp_el)
#endif

     !visc/adv
     udx=sum(dx*uvel)
     udy=sum(dy*uvel)
     udz=sum(dz*uvel)
     vdx=sum(dx*vvel)
     vdy=sum(dy*vvel)
     vdz=sum(dz*vvel)
#ifdef use_non_hydrostatic
     wdx=sum(dx*wvel)
     wdy=sum(dy*wvel)
     wdz=sum(dz*wvel)
#endif

     !viscosity
     vtr=voltriangle(elem2)
     if(smagorinsky_visc) then
	d1=udx-vdy
        d2=udy+vdx
        Ah=sqrt(d1*d1+d2*d2)*vtr*4.0
        Ah=max(Ah, 40.0/1.0e8*vtr)            !limit from below: 20m^2/s on 10km
	if(Ah*dt/vtr > 0.05) Ah=0.05*vtr/dt   !limit from above
     else
        Ah=Ah0
        if(scale_mixing_h) Ah=Ah*(vtr/scalevol)**(1.0/real(scale_mixing_type))
     endif
     if(use_vertvisc_impl) then
        Av_el=0.0
     else
        Av_el=sum(Av(elnodes))/4.0
     endif

     !required variables for free surface
#ifdef use_fullfreesurf
     uusum=sum(uvel*uvel)
     vvsum=sum(vvel*vvel)
     uvsum=sum(uvel*vvel)
#ifdef use_non_hydrostatic
     wusum=sum(wvel*uvel)
     wvsum=sum(wvel*vvel)
     wwsum=sum(wvel*wvel)
#endif
     n_el=map_elem(elem)
     flag_move=.false.
     if(n_el/=0) then
        flag_move=.true.
        coe_diff=inv20*(v-voltetra_new(n_el))
        wmove=0.
        do i=1, 4   
           if(index_nod3d(elnodes(i))>12) cycle
           row=nod2d_corresp_to_nod3d(elnodes(i))
           wmove(i)=-(ssh(row)-ssh0(row))*dt_inv 
        end do
        wmove_sum=sum(wmove)
        uwmove_sum=sum(uvel*wmove)
        vwmove_sum=sum(vvel*wmove)
#ifdef use_non_hydrostatic
        wwmove_sum=sum(wvel*wmove)
#endif
     end if
#endif

     !for along sigma mixing
     if (elem_type==1) then
        lay=elem3d_layer(elem)      
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux=S(1)**2+S(2)**2
        S(3)=sqrt(aux)
        rotate_coe=1.0/(1.0+aux) 
     end if

     coe=inv20*val
     do q=1,4
        row=elnodes(q)
        row2=row+n3                       
#ifdef use_non_hydrostatic
        row3=row2+n3                         
#endif
        aux=(um*dx(q)+vm*dy(q)+wm*dz(q))*dt*0.5_8*val
#ifdef use_fullfreesurf
        if(flag_move) then
           aux=aux+wmove_sum*inv4*dz(q)*dt*0.5_8*val  
        end if
#endif  
        !
#ifdef use_fullfreesurf
        if(flag_move) then
           uv_rhs(row)=uv_rhs(row)+(usum+uvel(q))*coe_diff
           uv_rhs(row2)=uv_rhs(row2)+(vsum+vvel(q))*coe_diff
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+(wsum+wvel(q))*coe_diff
#endif   
        endif
#endif
        !coriolis force
        cori_u_tmp=fc*(vsum+vvel(q))*coe
        cori_v_tmp=-fc*(usum+uvel(q))*coe 
        !stabilization - coriolis
        cori_u_tmp=cori_u_tmp+aux*fc*vm
        cori_v_tmp=cori_v_tmp-aux*fc*um	   

        !hydrostatic pressure
        urhs_temp=-pxsum*val*inv4
        vrhs_temp=-pysum*val*inv4
        !stabilization - hp
        urhs_temp=urhs_temp-aux*pxsum
        vrhs_temp=vrhs_temp-aux*pysum

        !ssh/non-hydrostatic pressure
#ifndef use_non_hydrostatic
        urhs_temp=urhs_temp-inv4*val*px_ssh*gamma_stab
        vrhs_temp=vrhs_temp-inv4*val*py_ssh*gamma_stab
        !stabilization - ssh
        urhs_temp=urhs_temp-aux*px_ssh
        vrhs_temp=vrhs_temp-aux*py_ssh
#else
        urhs_temp=urhs_temp-inv4*val*(px_ssh*gamma_stab+npx*gamma_stab_nh)
        vrhs_temp=vrhs_temp-inv4*val*(py_ssh*gamma_stab+npy*gamma_stab_nh)
        wrhs_temp=-inv4*val*npz*gamma_stab_nh
        !stabilization - ssh/nhp
        urhs_temp=urhs_temp-aux*(px_ssh+npx)
        vrhs_temp=vrhs_temp-aux*(py_ssh+npy)
        wrhs_temp=wrhs_temp-aux*npz
#endif

        !viscosity
        if(elem_type==1) then  !along-sigma viscosity
           !diagonal part1 (lateral)
	   temp_coe=Ah*rotate_coe*val
           urhs_temp=urhs_temp-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
           vrhs_temp=vrhs_temp-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
           !diagonal part2 (cross slope)
	   temp_coe2=Av_el*val+S(3)*S(3)*temp_coe
           urhs_temp=urhs_temp-dz(q)*udz*temp_coe2
           vrhs_temp=vrhs_temp-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp-dz(q)*wdz*temp_coe2
#endif
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           urhs_temp=urhs_temp-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
           vrhs_temp=vrhs_temp-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           urhs_temp=urhs_temp + (dx(q)*udy+dy(q)*udx)*temp_coe2
           vrhs_temp=vrhs_temp + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           urhs_temp=urhs_temp-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
           vrhs_temp=vrhs_temp-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
        else
           if(increase_equ_zonal_visc .and. abs(geolat(row))<5.*rad) then
              urhs_temp=urhs_temp-(fac_visc_increase*Ah*dx(q)*udx + &
                   Ah*dy(q)*udy + Av_el*dz(q)*udz)*val
              vrhs_temp=vrhs_temp-(fac_visc_increase*Ah*dx(q)*vdx + &
                   Ah*dy(q)*vdy + Av_el*dz(q)*vdz)*val
           else
              urhs_temp=urhs_temp-(Ah*dx(q)*udx + &
                   Ah*dy(q)*udy + Av_el*dz(q)*udz)*val
              vrhs_temp=vrhs_temp-(Ah*dx(q)*vdx + &
                   Ah*dy(q)*vdy + Av_el*dz(q)*vdz)*val
           end if
#ifdef use_non_hydrostatic
           wrhs_temp=wrhs_temp-(Ah*dx(q)*wdx + &
                Ah*dy(q)*wdy+Av_el*dz(q)*wdz)*val
#endif
        endif

        !advection
#ifndef use_non_hydrostatic
#ifndef use_fullfreesurf
        adv=(usum+uvel(q))*udx*inv5+(vsum+vvel(q))*udy*inv5+wm*udz
        urhs_temp=urhs_temp-adv*inv4*val
        adv=(usum+uvel(q))*vdx*inv5+(vsum+vvel(q))*vdy*inv5+wm*vdz 
        vrhs_temp=vrhs_temp-adv*inv4*val
#else
        adv=(dx(q)*(usum*usum+uusum)+dy(q)*(vsum*usum+uvsum))*inv20
        adv=adv+dz(q)*wm*usum*inv4
        if(flag_move) adv=adv+dz(q)*(wmove_sum*usum+uwmove_sum)*inv20
        urhs_temp=urhs_temp+adv*val
        adv=(dx(q)*(usum*vsum+uvsum)+dy(q)*(vsum*vsum+vvsum))*inv20
        adv=adv+dz(q)*wm*vsum*inv4
        if(flag_move) adv=adv+dz(q)*(wmove_sum*vsum+vwmove_sum)*inv20
        vrhs_temp=vrhs_temp+adv*val
#endif
        !stabilization - adv
        adv=um*udx+vm*udy+wm*udz
        urhs_temp=urhs_temp-aux*adv
        adv=um*vdx+vm*vdy+wm*vdz
        vrhs_temp=vrhs_temp-aux*adv
#ifdef use_fullfreesurf
        if(flag_move) then
           urhs_temp=urhs_temp-aux*wmove_sum*inv4*udz
           vrhs_temp=vrhs_temp-aux*wmove_sum*inv4*vdz
        end if
#endif 
#else
#ifndef use_fullfreesurf
        adv=(usum+uvel(q))*udx+(vsum+vvel(q))*udy+(wsum+wvel(q))*udz
        urhs_temp=urhs_temp-adv*coe
        adv=(usum+uvel(q))*vdx+(vsum+vvel(q))*vdy+(wsum+wvel(q))*vdz
        vrhs_temp=vrhs_temp-adv*coe
        adv=(usum+uvel(q))*wdx+(vsum+vvel(q))*wdy+(wsum+wvel(q))*wdz
        wrhs_temp=wrhs_temp-adv*coe
#else
        adv=dx(q)*(usum*usum+uusum)+dy(q)*(vsum*usum+uvsum) + &
             dz(q)*(wsum*usum+wusum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*usum+uwmove_sum)
        urhs_temp=urhs_temp+adv*coe
        adv=dx(q)*(usum*vsum+uvsum)+dy(q)*(vsum*vsum+vvsum) + &
             dz(q)*(wsum*vsum+wvsum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*vsum+vwmove_sum)
        vrhs_temp=vrhs_temp+adv*coe
        adv=dx(q)*(usum*wsum+wusum)+dy(q)*(vsum*wsum+wvsum) + &
             dz(q)*(wsum*wsum+wwsum)
        if(flag_move) adv=adv+dz(q)*(wmove_sum*wsum+wwmove_sum)
        wrhs_temp=wrhs_temp+adv*coe
#endif
        !stabilization - adv
        adv=um*udx+vm*udy+wm*udz
        urhs_temp=urhs_temp-aux*adv
        adv=um*vdx+vm*vdy+wm*vdz
        vrhs_temp=vrhs_temp-aux*adv
        adv=um*wdx+vm*wdy+wm*wdz
        wrhs_temp=wrhs_temp-aux*adv
#ifdef use_fullfreesurf
        if(flag_move) then
           urhs_temp=urhs_temp-aux*wmove_sum*inv4*udz
           vrhs_temp=vrhs_temp-aux*wmove_sum*inv4*vdz
           wrhs_temp=wrhs_temp-aux*wmove_sum*inv4*wdz
        end if
#endif 
#endif

        !put rhs temp values to the uv_rhs array
        if(use_cori_semi) then
           urhs_temp=urhs_temp+cori_u_tmp
           vrhs_temp=vrhs_temp+cori_v_tmp
           uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
           uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+wrhs_temp*dt
#endif
        else
           uv_rhs(row)=uv_rhs(row)+urhs_temp
           uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
           ucori(row)=ucori(row)+cori_u_tmp
           vcori(row)=vcori(row)+cori_v_tmp
#ifdef use_non_hydrostatic
           uv_rhs(row3)=uv_rhs(row3)+wrhs_temp
#endif
        endif
     end do
  end do


  !2D contributions: surface forcing and bottom stress
  if(.not.use_vertvisc_impl) then
     do elem2=1, myDim_elem2d    
        v=voltriangle(elem2)
        if(use_cori_semi) then
           val=v
           !coriolis parameter
           fc=coriolis_param_elem2D(elem2)
           !terms for semi-implicit schemes
           dparam=dt_inv**2 + alpha_trapez**2*fc**2
           beta=dt_inv/dparam
           gamma=fc*alpha_trapez/dparam
        else
           val=v*dt
        endif

        !surface             
        !Wind/ice stress contribution, or stress from ice-cavity    
        elnodes2=elem2D_nodes(:,elem2)   
#ifdef use_cavity
        if(all(cavity_flag_nod2d(elnodes2)==0)) then   
#endif
           usum=sum(stress_x(elnodes2))
           vsum=sum(stress_y(elnodes2))
           aux=rho0r*inv12*val
           do q=1,3
              row=elnodes2(q)
              urhs_temp=aux*(usum+stress_x(row)) 
              vrhs_temp=aux*(vsum+stress_y(row))
              row=nod3D_below_nod2D(1,elnodes2(q))       
              row2=row+n3                                 
              if(use_cori_semi) then
                 uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
              else
                 uv_rhs(row)=uv_rhs(row)+urhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
              endif
           end do
#ifdef use_cavity
        else
           uvel(1:3)=uf(nod3D_below_nod2D(1,elnodes2))      
           vvel(1:3)=uf(nod3D_below_nod2D(1,elnodes2)+n3)   
           usum=sum(uvel(1:3))
           vsum=sum(vvel(1:3))
           aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*val
           do q=1,3
              row=nod3D_below_nod2D(1,elnodes2(q))          
              row2=row+n3                                   
              urhs_temp=aux*(usum+uvel(q))
              vrhs_temp=aux*(vsum+vvel(q))
              if(use_cori_semi) then
                 uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
              else
                 uv_rhs(row)=uv_rhs(row)+urhs_temp
                 uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
              endif
           end do
        end if
#endif

        !bottom 
        !Bottom drag contribution  -C_d U|U| 
        elnodes2=bt_nds(elnodes2)        
        uvel(1:3)=uf(elnodes2)
        vvel(1:3)=uf(elnodes2+n3)                           
        usum=sum(uvel(1:3))
        vsum=sum(vvel(1:3))
        aux=-C_d*sqrt(usum**2+vsum**2)*inv3*inv12*val
        do q=1,3
           row=elnodes2(q)
           row2=row+n3                                     
           urhs_temp=aux*(usum+uvel(q))
           vrhs_temp=aux*(vsum+vvel(q))
           if(use_cori_semi) then
              uv_rhs(row)=uv_rhs(row)+beta*urhs_temp+gamma*vrhs_temp
              uv_rhs(row2)=uv_rhs(row2)+beta*vrhs_temp-gamma*urhs_temp
           else
              uv_rhs(row)=uv_rhs(row)+urhs_temp
              uv_rhs(row2)=uv_rhs(row2)+vrhs_temp
           endif
        end do
     end do
  end if


  if(biharmonic_visc) call biharmonic_viscosity


  !update A-B part (coriolis) and vertical boundary conditions
  do row=1, myDim_nod3d              
     row2=row+n3                           
     if(.not.use_cori_semi) then
        uv_rhs(row)=uv_rhs(row)+alpha_AB*ucori(row)+(1.0-alpha_AB)*ucori_back(row)
        uv_rhs(row2)=uv_rhs(row2)+alpha_AB*vcori(row)+(1.0-alpha_AB)*vcori_back(row)
        ucori_back(row)=ucori(row)
        vcori_back(row)=vcori(row)
        ucori(row)=0.
        vcori(row)=0.
     endif
     ! boundary conditions
     ind=index_nod3d(row)
     if((ind==11).or.(ind==21).or.(ind==31)) then
        uv_rhs(row)=0.0
        uv_rhs(row2)=0.0
#ifdef use_non_hydrostatic
        uv_rhs(row2+n3)=0.0              
#endif 
     end if
  end do
end subroutine velocity_rhs
!
!----------------------------------------------------------------------------
!
subroutine velocity_rhs_update
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer                        :: m, elem, el2, elnodes(4), elnodes23(4)
  integer                        :: row, row2, ind, n3
  real(kind=8)                   :: aux, vc(4), dx(4), dy(4), npx, npy
  real(kind=8)                   :: inv4
  real(kind=8)                   :: cori_p, dparam, beta, gamma
#ifdef use_non_hydrostatic
  real(kind=8)                   :: dz(4), npz
#endif
  !
  inv4=0.25_8
  n3=myDim_nod3D+eDim_nod3D          
  do row=1, myDim_nod3d                 
     row2=row+n3                     
     uv_rhs(row)=0.
     uv_rhs(row2)=0.
#ifdef use_non_hydrostatic
     uv_rhs(row2+n3)=0.               
#endif
  enddo
  !
  do elem=1, myDim_elem3d             
     elnodes=elem3D_nodes(:,elem)
     elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
#ifdef use_non_hydrostatic
     dz=bafuz_3D(:,elem)
#endif
     if(use_cori_semi) then
        aux=voltetra(elem)*inv4
        !coriolis parameter
        el2=elem2d_corresp_to_elem3d(elem)
        cori_p=coriolis_param_elem2d(el2)
        !terms for semi-implicit schemes
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
     else
        aux=voltetra(elem)*inv4*dt
     endif
     !ssh/non-hydrostatic pressure:
#ifndef use_non_hydrostatic
#ifdef use_semiimplicit_scheme
     vc=ssh0(elnodes23)*(gamma_stab-1.0+theta_ssh)-ssh(elnodes23)*theta_ssh
#else
     vc=ssh0(elnodes23)*gamma_stab-ssh(elnodes23)
#endif
     npx=g*sum(dx*vc)
     npy=g*sum(dy*vc)
#else
     vc=g*ssh0(elnodes23)*gamma_stab+nhp0(elnodes)*gamma_stab_nh - &
          (g*ssh(elnodes23)+nhp(elnodes))
     npx=sum(dx*vc)
     npy=sum(dy*vc)
     npz=sum(dz*vc) 
#endif
     !
     !assembling over nodes
     if(use_cori_semi) then
        uv_rhs(elnodes)=uv_rhs(elnodes)+(beta*npx+gamma*npy)*aux
        uv_rhs(elnodes+n3)=uv_rhs(elnodes+n3)+(beta*npy-gamma*npx)*aux 
#ifdef use_non_hydrostatic   
        uv_rhs(elnodes+2*n3)=uv_rhs(elnodes+2*n3)+npz*aux*dt            
#endif
     else
        uv_rhs(elnodes)=uv_rhs(elnodes)+npx*aux
        uv_rhs(elnodes+n3)=uv_rhs(elnodes+n3)+npy*aux                 
#ifdef use_non_hydrostatic   
        uv_rhs(elnodes+2*n3)=uv_rhs(elnodes+2*n3)+npz*aux              
#endif  
     endif
  end do
  !
  ! boundary conditions
  do row=1,myDim_nod3d                      
     row2=row+n3                           
     ind=index_nod3d(row)
     if((ind==11).or.(ind==21).or.(ind==31)) then
        uv_rhs(row)=0.0
        uv_rhs(row2)=0.0 
#ifdef use_non_hydrostatic
        uv_rhs(row2+n3)=0.0                
#endif 
     end if
  end do
end subroutine velocity_rhs_update
!
!----------------------------------------------------------------------------
!
subroutine impl_vertvisc
  ! apply separated implicit vertical viscosity
  ! lumped mass matrix is used
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none

  integer            :: m, n2, nlay, k, i, q, p, cnt
  integer            :: row, elem, elnodes(4), n_el
  integer            :: nodup, nodlo
  real(kind=8)       :: vol, vol_visc, inv4, Av_el
  real(kind=8)       :: a(max_num_layers), b(max_num_layers), c(max_num_layers)
  real(kind=8)       :: ru(max_num_layers), rv(max_num_layers)

  inv4=0.25_8

  do n2=1,myDim_nod2d  !only myDim_nod2d required

     if(mod(index_nod3d(n2),10)==1) cycle  !no-slip bc.	

     nlay=num_layers_below_nod2d(n2)+1

     ! assemble three diagonal matrix and rhs
     do k=1,nlay
        a(k)=0.
        b(k)=0.
        c(k)=0.
        ru(k)=0.
        rv(k)=0.

        row=nod3d_below_nod2d(k,n2)

        if(mod(index_nod3d(row),10)==1) then
           b(k)=1.0
           cycle !no-slip bc.
        end if

        if(k>1) nodup=nod3d_below_nod2d(k-1,n2)
        if(k<nlay) nodlo=nod3d_below_nod2d(k+1,n2)
        do i=1,nod_in_elem3D(row)%nmb
           elem=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,elem)
           vol_visc=voltetra(elem)
#ifdef use_fullfreesurf
           n_el=map_elem(elem)
           if(n_el==0) then
              vol=voltetra(elem)
           else
              vol=voltetra_new(n_el)
           end if
#else
           vol=voltetra(elem)
#endif

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do

           cnt=0
           Av_el=0.0

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    Av_el=Av(nodup)    !Av(nodup) is the Kv for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Av_el=Av_el+Av0
                    end if
                    a(k)=a(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Av_el * vol_visc
                    cnt=1
                    exit
                 end if
              end do
           end if

           !third entry
           if(k<nlay) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    Av_el=Av(row)    !Av(row) is the Kv for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Av_el=Av_el+Av0
                    end if
                    c(k)=c(k) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Av_el * vol_visc
                    cnt=1
                    exit
                 end if
              end do
           end if

           !second entry
           b(k)=b(k) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Av_el * vol_visc

           !second entry, mass matrix part
           b(k)=b(k) + vol*inv4*dt_inv
           !rhs
           ru(k)=ru(k) + uf(row)*vol*inv4*dt_inv
           rv(k)=rv(k) + uf(myDim_nod3d+eDim_nod3D+row)*vol*inv4*dt_inv 

        end do  !i
     end do  !k

     ! surface and bottom boundary conditions
     ru(1)=ru(1) + uv_sfc_force(n2,1)
     rv(1)=rv(1) + uv_sfc_force(n2,2)
     ru(nlay)=ru(nlay) + uv_bott_force(n2,1)
     rv(nlay)=rv(nlay) + uv_bott_force(n2,2)

     ! the sweep algorithm
     ! forward step
     ! prepare coefficients 
     c(1)=-c(1)/b(1)
     ru(1)=ru(1)/b(1)
     rv(1)=rv(1)/b(1)
     do k=2,nlay
        b(k)=1.0/(a(k)*c(k-1)+b(k))
        c(k)=-c(k)*b(k)
        ru(k)=(ru(k)-a(k)*ru(k-1))*b(k)
        rv(k)=(rv(k)-a(k)*rv(k-1))*b(k)
     end do
     ! backward sweep
     do k=nlay-1,1,-1
        ru(k)=ru(k)+c(k)*ru(k+1)
        rv(k)=rv(k)+c(k)*rv(k+1)
     end do

     ! update the velocity
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        uf(row)=ru(k)
        uf(row+myDim_nod3d+eDim_nod3D)=rv(k)   
     end do

  end do  !m

  call com_3D(uf(1:myDim_nod3d+eDim_nod3D))     
  call com_3D(uf(1+myDim_nod3d+eDim_nod3D:2*(myDim_nod3d+eDim_nod3D))) 

end subroutine impl_vertvisc
!
!============================================================================
!
subroutine compute_ssh_rhs
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer             :: elem, el2, elnodes(4), elnodes2(3)
  integer             :: q, row, m, mn(3), n3
  real(kind=8)        :: vol, dx(3), dy(3), us, vs
  real(kind=8)        :: inv4, inv12, aux, aux1, aux2
  real(kind=8)        :: u_el(4), v_el(4), tri_u(3), tri_v(3)
  real(kind=8)        :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer             :: n_el
#endif
  !
  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D       
  do row=1,myDim_nod2d            
     ssh_rhs(row)=0.
  enddo

  ! divergence contribution
  do elem=1,myDim_elem3d        
     el2=elem2d_corresp_to_elem3d(elem)
     elnodes2=elem2d_nodes(:,el2)
     elnodes=elem3D_nodes(:,elem)
     !elnodes23=nod2D_corresp_to_nod3D(elnodes)
     dx=bafux_2D(:,el2)
     dy=bafuy_2D(:,el2)
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        vol=voltetra(elem)
     else
        vol=voltetra_new(n_el)
     end if
#else
     vol=voltetra(elem)
#endif    
#ifdef use_semiimplicit_scheme
     u_el=uf(elnodes)*theta_vel+uf0(elnodes)*(1.0-theta_vel)
     v_el=uf(elnodes+n3)*theta_vel+uf0(elnodes+n3)*(1.0-theta_vel)   
#else
     u_el=uf(elnodes)
     v_el=uf(elnodes+n3)                                           
#endif
     !
     aux=g*(gamma_stab-1.0)
#ifdef use_semiimplicit_scheme
     aux=aux*theta_vel
#endif
     aux1=sum(dx*ssh0(elnodes2))*aux
     aux2=sum(dy*ssh0(elnodes2))*aux
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2d(el2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        us=sum(u_el)*inv4+beta*aux1+gamma*aux2
        vs=sum(v_el)*inv4+beta*aux2-gamma*aux1 
     else
        us=sum(u_el)*inv4+dt*aux1
        vs=sum(v_el)*inv4+dt*aux2 
     endif
     !
     do q=1,3
	row=elnodes2(q)
	ssh_rhs(row)=ssh_rhs(row)+(dx(q)*us+dy(q)*vs)*vol
     end do
  end do

#ifdef use_fullfreesurf
  ! P-E contribution
  do el2=1,myDim_elem2d                
     elnodes2=elem2D_nodes(:,el2)
     vs=sum(water_flux(elnodes2))  ! '+' is defined as upwards 
     vol=voltriangle(el2)*inv12
     do q=1, 3
        row=elnodes2(q)
        ssh_rhs(row)=ssh_rhs(row) - (vs+water_flux(row))*vol
     end do
  end do
#endif

  ! open boundary contribution
#ifdef use_opbnd_tide
  if(trim(tide_opbnd_type)=='ssh') then
     do q=1, nmbr_opbnd_t2d
        row=opbnd_n2d(q)
        ssh_rhs(row)=opbnd_z_tide(q)-opbnd_z0_tide(q)  ! increments
     end do
  else
     do el2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(el2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(el2,4)*inv12
	tri_u=opbnd_u_tide(mn)
 	tri_v=opbnd_v_tide(mn)
        tri_v=tri_u*opbnd_nv(el2,1)+tri_v*opbnd_nv(el2,2)
        if(trim(tide_opbnd_type)=='Flather') then
           tri_v=tri_v + sqrt(g/opbnd_dep(mn))*(ssh(elnodes2)-opbnd_z_tide(mn)) 
        end if
        vs=sum(tri_v)
        do q=1, 3
	   row=elnodes2(q)
           ssh_rhs(row)=ssh_rhs(row) - (vs+tri_v(q))*vol 
        end do
     end do
  end if
#endif
#ifdef use_opbnd_restoring
  do q=1, nmbr_opbnd_t2d
     row=opbnd_n2D(q)
     m=mapping_opbnd_n2d(row)
     ssh_rhs(row)=ssh_rhs(row)+opbnd_ssh_rhs(m)
  end do
#endif

end subroutine compute_ssh_rhs
!===================================================================	      

#ifdef use_non_hydrostatic
subroutine compute_nhp_rhs
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer             :: elem, elem2, q, row, m, n3
  integer             :: elnodes(4), elnodes23(4), elnodes2(3)
  real(kind=8)        :: vol, dx(4), dy(4), dz(4), us, vs, ws
  real(kind=8)        :: inv4, inv12, vc(4), aux1, aux2
  real(kind=8)        :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer             :: n_el
#endif
  !
  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1,myDim_nod3d              
     nhp_rhs(row)=0.
  enddo
  !
  do elem=1,myDim_elem3d          
     elnodes=elem3D_nodes(:,elem)
     elnodes23=nod2d_corresp_to_nod3d(elnodes)
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        dx=bafux_3D(:,elem)
        dy=bafuy_3D(:,elem)
        dz=bafuz_3D(:,elem)
        vol=voltetra(elem)
     else
        dx=bafux_3D_new(:,n_el)
        dy=bafuy_3D_new(:,n_el)
        dz=bafuz_3D_new(:,n_el)
        vol=voltetra_new(n_el)
     end if
#else
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     dz=bafuz_3D(:,elem)
     vol=voltetra(elem)
#endif    
     !
     vc=nhp(elnodes)*gamma_stab_nh+g*ssh0(elnodes23)*gamma_stab-g*ssh(elnodes23)
     aux1=sum(dx*vc)
     aux2=sum(dy*vc)
     if(use_cori_semi) then
        elem2=elem2d_corresp_to_elem3d(elem)
        cori_p=coriolis_param_elem2d(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        us=sum(uf(elnodes))*inv4 + beta*aux1+gamma*aux2
        vs=sum(uf(elnodes+n3))*inv4 + beta*aux2-gamma*aux1    
     else
        us=sum(uf(elnodes))*inv4 + dt*aux1
        vs=sum(uf(elnodes+n3))*inv4 + dt*aux2                   
     endif
     ws=sum(uf(elnodes+2*n3))*inv4 + dt*sum(dz*nhp(elnodes))*gamma_stab_nh 
     do q=1,4
	row=elnodes(q)
	nhp_rhs(row)=nhp_rhs(row)+(dx(q)*us+dy(q)*vs+dz(q)*ws)*vol
     end do
  end do
  !
  ! surface integral 
  do elem2=1, myDim_elem2d               
     vol=voltriangle(elem2)*inv12
     elnodes2=elem2D_nodes(:,elem2)
     us=sum(dssh(elnodes2))
     vs=sum(water_flux(elnodes2))
     do q=1,3
        row=elnodes2(q)
        nhp_rhs(row)  = nhp_rhs(row) - &
             vol*((us+dssh(row))*dt_inv+vs+water_flux(row))
     end do
  end do
end subroutine compute_nhp_rhs
#endif
!===================================================================
!
subroutine uv_solve
  use o_matrices
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                         :: n, m, cn, clo, clo2, location(100)
  integer                         :: row, row2, n3
  real(kind=8)                    :: u_rhs_new
#ifdef use_non_hydrostatic
  integer                         :: row3
  real(kind=8), allocatable       :: temp_array(:)
  allocate(temp_array(myDim_nod3d))
#endif
  !
  n3=myDim_nod3D+eDim_nod3D
  !the first approximation
  do row=1,myDim_nod3d                   
     row2=row+n3    	                  
     duf(row)=uv_rhs(row)/uv_lump(row)
     duf(row2)=uv_rhs(row2)/uv_lump(row)
#ifdef use_non_hydrostatic
     row3=row2+n3                         
     duf(row3)=uv_rhs(row3)/uv_lump(row)
#endif
  end do
  call com_3d(duf(1:n3))                   
  call com_3d(duf(1+n3:2*n3))             
#ifdef use_non_hydrostatic
  call com_3d(duf(1+2*n3:3*n3))           
#endif
  !

  !iterate

  do n=1, num_iter_solve-1
     do row=1,myDim_nod3d                 
 	row2=row+n3                     
        clo=uvstiff%rowptr(row)-uvstiff%rowptr(1)+1     
        clo2=uvstiff%rowptr(row+1)-uvstiff%rowptr(1)    
        cn=clo2-clo+1
        location(1:cn)=nghbr_nod3D(row)%addresses      
        u_rhs_new=uv_rhs(row) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)))
        dtracer(row,1)=duf(row)+u_rhs_new/uv_lump(row)
	u_rhs_new=uv_rhs(row2) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)+n3)) 
        dtracer(row,2)=duf(row2)+u_rhs_new/uv_lump(row)
#ifdef use_non_hydrostatic
        row3=row2+n3                                    
        u_rhs_new=uv_rhs(row3) - sum(uvstiff%values(clo:clo2)*duf(location(1:cn)+2*n3))
      	temp_array(m)=duf(row3)+u_rhs_new/uv_lump(row)
#endif
     end do
     do row=1,myDim_nod3d                  
	row2=row+n3                     
        duf(row)=dtracer(row,1)
        duf(row2)=dtracer(row,2)
#ifdef use_non_hydrostatic
        row3=row2+n3                     
        duf(row3)=temp_array(m)
#endif
     end do
     call com_3d(duf(1:n3))              
     call com_3d(duf(1+n3:2*n3))         
#ifdef use_non_hydrostatic
     call com_3d(duf(1+2*n3:3*n3))       
#endif

  end do
  !
#ifdef use_non_hydrostatic
  deallocate(temp_array)
#endif  
  !
end subroutine uv_solve
!
!========================================================================
!
#ifndef use_non_hydrostatic
subroutine compute_vvel_rhs
  !Compute RHS in vertical velocity equation
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                 :: i, m, n, q, elem, elem2, row, n3
  integer                 :: elnodes(4), elnodes2(3), mn(3)
  real(kind=8)            :: val, ux, vy, inv4, inv12, dx(4), dy(4)
  real(kind=8)            :: vc(3), aux1, aux2, tri_u(3), tri_v(3)
  real(kind=8)            :: cori_p, dparam, beta, gamma
#ifdef use_fullfreesurf
  integer                 :: n_el
#endif  

  inv4=0.25_8
  inv12=1.0_8/12.0_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1,myDim_nod3D              
     wrhs(row)=0.0
  enddo

  ! compute -\nabla\tilde\phi v* d\Omega   
  do elem=1, myDim_elem3d             
     elem2=elem2d_corresp_to_elem3d(elem)
     elnodes=elem3D_nodes(:,elem)
     elnodes2=elem2D_nodes(:,elem2)
     ux=inv4*sum(uf0(elnodes))
     vy=inv4*sum(uf0(elnodes+n3))    
#ifdef use_fullfreesurf
     n_el=map_elem(elem)
     if(n_el==0) then
        dx=bafux_3D(:,elem)
        dy=bafuy_3D(:,elem)
        val=voltetra(elem)
     else
        dx=bafux_3D_new(:,n_el)
        dy=bafuy_3D_new(:,n_el)
        val=voltetra_new(n_el)
     endif
#else
     dx=bafux_3D(:,elem)
     dy=bafuy_3D(:,elem)
     val=voltetra(elem)
#endif
     !
#ifdef use_semiimplicit_scheme
     vc=-g*(theta_ssh*ssh(elnodes2)-(gamma_stab-1.0+theta_ssh)*ssh0(elnodes2))
#else
     vc=-g*(ssh(elnodes2)-gamma_stab*ssh0(elnodes2))
#endif
     aux1=sum(vc*bafux_2D(:,elem2))
     aux2=sum(vc*bafuy_2D(:,elem2))
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2d(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        ux=ux+beta*aux1+gamma*aux2
        vy=vy+beta*aux2-gamma*aux1
     else
        ux=ux+dt*aux1
        vy=vy+dt*aux2
     endif
     !
     do q=1,4
        row=elnodes(q)
        wrhs(row) = wrhs(row)  - val*(dx(q)*ux+dy(q)*vy)
     end do
  end do

  ! sea surface boundary condition/integral
  do elem2=1, myDim_elem2d                   
     val=voltriangle(elem2)*inv12
     elnodes2=elem2D_nodes(:,elem2)
     ux=sum(dssh(elnodes2))
#ifdef use_fullfreesurf
     vy=sum(water_flux(elnodes2))
#endif
     do q=1,3
        row=elnodes2(q)
        n=nod3d_below_nod2d(1,row)
        wrhs(n) = wrhs(n) + &
             val*(ux+dssh(row))*dt_inv
#ifdef use_fullfreesurf
        wrhs(n) = wrhs(n) + val*(vy+water_flux(row))
#endif   
     end do
  end do


!!$  ! open boundary contribution
!!$#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
!!$  do elem2=1, nmbr_opbnd_tri
!!$     elnodes2=opbnd_tri(elem2,1:3)  !contains the 3d nodes on the opbnd.
!!$     val=opbnd_nv(elem2,4)*inv12
!!$     tri_u=uf(elnodes2)
!!$     tri_v=uf(elnodes2+nod3d)
!!$     tri_v=tri_u*opbnd_nv(elem2,1)+tri_v*opbnd_nv(elem2,2)
!!$     vy=sum(tri_v)
!!$     do q=1,3
!!$        row=elnodes2(q)
!!$        wrhs(row)=wrhs(row) + (vy+tri_v(q))*val 
!!$     end do
!!$  end do
!!$#endif


  ! Check solvability
  do row=1, myDim_nod2d     
     val=0.0
     do n=1,num_layers_below_nod2D(row)+1
        m = nod3D_below_nod2D(n,row)
        val=val+wrhs(m) 
     end do
     ! As a rule, val is a factor 10^(-4) smaller
     ! than wrhs(row), and is even smaller compared to 
     ! wrhs at nodes located deeper. Yet this accuracy
     ! is insufficient and thus solvability is enforced.
     val=val/real(num_layers_below_nod2D(row)+1)
     do n=1,num_layers_below_nod2D(row)+1
        m = nod3D_below_nod2D(n,row)
        wrhs(m)= wrhs(m)-val
     end do
  end do

end subroutine compute_vvel_rhs
#endif
!
!========================================================================
!
subroutine solve_wpot
  ! solve vertical velocity potential using sweep algorithm
  use o_param
  use o_array
  use o_mesh
  use o_matrices
  use g_config
  use g_parfe
  implicit none

  integer            :: n2, nlay, k, row
  real(kind=8)       :: a(max_num_layers), b(max_num_layers)
  real(kind=8)       :: c(max_num_layers), rw(max_num_layers)

  do n2=1,myDim_nod2d                      
     nlay=num_layers_below_nod2d(n2)+1

     ! matrix entries
     a(1:nlay)=wpot_matrix(1,1:nlay,n2) 
     b(1:nlay)=wpot_matrix(2,1:nlay,n2)  
     c(1:nlay)=wpot_matrix(3,1:nlay,n2)  

     !rhs
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        rw(k)=wrhs(row)
     end do

     ! due to increased order of operator in this equation, we should/can 
     ! specify the boundary (surface) value. Otherwise, the full matrix is 
     ! singular and not solvable.
     rw(1)=0.0

     if(nlay>2) then
        ! the sweep algorithm
        ! forward step
        ! prepare coefficients 
        c(2)=-c(2)/b(2)
        rw(2)=rw(2)/b(2)
        do k=3,nlay
           b(k)=a(k)*c(k-1)+b(k)
           c(k)=-c(k)/b(k)
           rw(k)=(rw(k)-a(k)*rw(k-1))/b(k)
        end do
        ! backward sweep
        do k=nlay-1,2,-1
           rw(k)=rw(k)+c(k)*rw(k+1)
        end do
     else
        rw(2)=rw(2)/b(2)
     end if

     ! put results to w array
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        w(row)=rw(k)
     end do

  end do  !m

  call com_3D(w)

end subroutine solve_wpot
!
!========================================================================
!
#ifndef use_non_hydrostatic
subroutine vvel_nodes
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use g_config
  use g_parfe
  implicit none
  !
  integer   :: col, elem, n, m, elnodes(4)
  real      :: vol, dz(4), elemw
  !
  !  Computes approximate vertical velocities at nodes as mean over neighbour
  !  elements

  do col=1,myDim_nod3d              
     wrhs(col)=0.
     vol=0.0
     do n=1,nod_in_elem3D(col)%nmb
        elem = nod_in_elem3D(col)%addresses(n)
        elnodes=elem3D_nodes(:,elem)
        dz=bafuz_3D(:,elem)
        vol=vol+voltetra(elem)
        elemw=sum(dz*w(elnodes))
        wrhs(col)=wrhs(col)+elemw*voltetra(elem) 
     end do
     wrhs(col)=wrhs(col)/vol
  end do
end subroutine vvel_nodes
#endif        
!============================================================================   
!
subroutine biharmonic_viscosity
  use o_mesh
  use o_elements
  use o_param 
  use o_array
  use o_matrices
  use g_config
  use g_parfe
  implicit none
  integer                      :: m, q, row, row2, elnodes(4)
  integer                      :: elem, elem2, ind, lay, n3, elem_type
  real(kind=8)                 :: dx(4), dy(4), dz(4), vol
  real(kind=8)                 :: vtr, Ah, d1, d2
  real(kind=8)                 :: u_el(4), v_el(4), inv4 
  real(kind=8)                 :: udx, udy, udz, aux1
  real(kind=8)                 :: vdx, vdy, vdz, aux2
  real(kind=8)                 :: S(3), rotate_coe, temp_coe, temp_coe2
#ifdef use_non_hydrostatic
  integer                      :: row3
  real(kind=8)                 :: w_el(4), wdx, wdy, wdz, aux3
#endif
  real(kind=8)                 :: cori_p, dparam, beta, gamma
  !
  inv4=0.25_8
  n3=myDim_nod3D+eDim_nod3D
  do row=1, myDim_nod3d        
     row2=row+n3                
     duf(row)=0.
     duf(row2)=0.
#ifdef use_non_hydrostatic
     row3=row2+n3               
     duf(row3)=0.
#endif
  enddo
  !
  ! ==== store laplacian of velocity, with terms due to
  !      differentiation of metrics ignored
  do elem=1,myDim_elem3d        
     elnodes=elem3D_nodes(:,elem)
     elem2=elem2d_corresp_to_elem3d(elem)
     elem_type=grid_type_elem2d(elem2)
     vol=voltetra(elem)
     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     u_el=uf(elnodes)
     v_el=uf(n3+elnodes)       
     udx=sum(dx*u_el)
     udy=sum(dy*u_el)
     vdx=sum(dx*v_el)
     vdy=sum(dy*v_el)
#ifdef use_non_hydrostatic
     w_el=uf(n3*2+elnodes)    
     wdx=sum(dx*w_el)
     wdy=sum(dy*w_el)
#endif
     if (elem_type==1) then
        dz=bafuz_3d(:,elem)
        udz=sum(dz*u_el)
        vdz=sum(dz*v_el)
#ifdef use_non_hydrostatic
        wdz=sum(dz*w_el)
#endif
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if
     !
     do q=1,4 
        row=elnodes(q)
        row2=row+n3          
	ind=index_nod3D(row)
        if((ind==11).or.(ind==21).or.(ind==31)) cycle
        !
        if(elem_type==1) then  !along-sigma viscosity
           !diagonal part1 (lateral)
	   temp_coe=rotate_coe*vol
           aux1=-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
           aux2=-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
           !diagonal part2 (cross slope)
	   temp_coe2=S(3)*S(3)*temp_coe
           aux1=aux1-dz(q)*udz*temp_coe2
           aux2=aux2-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3-dz(q)*wdz*temp_coe2
#endif
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           aux1=aux1-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
           aux2=aux2-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           aux1=aux1 + (dx(q)*udy+dy(q)*udx)*temp_coe2
           aux2=aux2 + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3 + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           aux1=aux1-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
           aux2=aux2-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
        else
           aux1=-(dx(q)*udx + dy(q)*udy)*vol
           aux2=-(dx(q)*vdx + dy(q)*vdy)*vol
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx + dy(q)*wdy)*vol
#endif
        endif
        !
        duf(row)=duf(row) + aux1
        duf(row2)=duf(row2) +aux2
#ifdef use_non_hydrostatic
        row3=row2+n3                         
        duf(row3)=duf(row3) + aux3
#endif
     end do
  end do
  !
  ! ==== Return to physical space ===============
  do row=1,myDim_nod3d                           
     row2=row+n3                                 
     ind=index_nod3D(row)
     if((ind==11).or.(ind==21).or.(ind==31)) cycle
     duf(row)=duf(row)/uv_lump(row)  
     duf(row2)=duf(row2)/uv_lump(row)
#ifdef use_non_hydrostatic
     row3=row2+n3                                 
     duf(row3)=duf(row3)/uv_lump(row)
#endif
  end do
  !
  ! ==== Communication ==========================
  call com_3d(duf(1:n3))                          
  call com_3d(duf(1+n3:2*n3))                    
#ifdef use_non_hydrostatic
  call com_3d(duf(1+2*n3:3*n3))                  
#endif
  !
  ! ==== Compute Laplacian of Laplacian =========
  do elem=1,myDim_elem3d                     
     elnodes=elem3D_nodes(:,elem)
     elem2=elem2d_corresp_to_elem3d(elem)
     elem_type=grid_type_elem2d(elem2)
     if(use_cori_semi) then
        cori_p=coriolis_param_elem2D(elem2)
        dparam=dt_inv**2 + alpha_trapez**2*cori_p**2
        beta=dt_inv/dparam
        gamma=cori_p*alpha_trapez/dparam
        vol=voltetra(elem)
     else
        vol=voltetra(elem)*dt
     endif
     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     !
     u_el=duf(elnodes)
     v_el=duf(elnodes+n3)                       
     udx=sum(dx*u_el)
     udy=sum(dy*u_el)
     vdx=sum(dx*v_el)
     vdy=sum(dy*v_el)
#ifdef use_non_hydrostatic
     w_el=duf(n3*2+elnodes)                    
     wdx=sum(dx*w_el)
     wdy=sum(dy*w_el)
#endif
     !
     !biharmonic viscosity 
     vtr=voltriangle(elem2)
     if(smagorinsky_visc) then
	d1=sum(dx*uf(elnodes))-sum(dy*uf(elnodes+n3))  
	d2=sum(dy*uf(elnodes))+sum(dx*uf(elnodes+n3))  
	Ah=sqrt(d1*d1+d2*d2)*vtr*4.0           !laplacian viscosity
        Ah=max(Ah, 40.0/1.0e8*vtr)             !limit from below: 20m^2/s on 10km
	if(Ah*dt/vtr > 0.02) Ah=0.02*vtr/dt    !limit from above (0.05)
        Ah=-Ah*vtr/8.0                         !biharmonic viscosity
     else
        Ah=-Ahb0   
        if(scale_mixing_h) Ah=Ah*(sqrt(vtr/scalevol))**3
     endif

     !special: test increasing boundary Ah, Qiang
     if(any(mod(index_nod3d(elnodes),10)==1)) then
!        Ah=Ah*3.0
     end if

     if (elem_type==1) then
        dz=bafuz_3d(:,elem)
        udz=sum(dz*u_el)
        vdz=sum(dz*v_el)
#ifdef use_non_hydrostatic
        wdz=sum(dz*w_el)
#endif
        lay=elem3d_layer(elem)
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if
     !
     !
     do q=1, 4 
        row=elnodes(q)
        row2=row+n3                          
        !
        if(elem_type==1) then  !along-sigma viscosity
           !diagonal part1 (lateral)
	   temp_coe=Ah*rotate_coe*vol
           aux1=-(dx(q)*udx*(1.0+S(2)*S(2))+ &
                dy(q)*udy*(1.0+S(1)*S(1)))*temp_coe
           aux2=-(dx(q)*vdx*(1.0+S(2)*S(2))+ &
                dy(q)*vdy*(1.0+S(1)*S(1)))*temp_coe
#ifdef use_non_hydrostatic
           aux3=-(dx(q)*wdx*(1.0+S(2)*S(2))+ &
                dy(q)*wdy*(1.0+S(1)*S(1)))*temp_coe
#endif
           !diagonal part2 (cross slope)
	   temp_coe2=S(3)*S(3)*temp_coe
           aux1=aux1-dz(q)*udz*temp_coe2
           aux2=aux2-dz(q)*vdz*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3-dz(q)*wdz*temp_coe2
#endif
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           aux1=aux1-(S(1)*dx(q)+S(2)*dy(q))*udz*temp_coe
           aux2=aux2-(S(1)*dx(q)+S(2)*dy(q))*vdz*temp_coe
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*dx(q)+S(2)*dy(q))*wdz*temp_coe
#endif
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           aux1=aux1 + (dx(q)*udy+dy(q)*udx)*temp_coe2
           aux2=aux2 + (dx(q)*vdy+dy(q)*vdx)*temp_coe2
#ifdef use_non_hydrostatic
           aux3=aux3 + (dx(q)*wdy+dy(q)*wdx)*temp_coe2
#endif
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           aux1=aux1-(S(1)*udx+S(2)*udy)*dz(q)*temp_coe
           aux2=aux2-(S(1)*vdx+S(2)*vdy)*dz(q)*temp_coe 
#ifdef use_non_hydrostatic
           aux3=aux3-(S(1)*wdx+S(2)*wdy)*dz(q)*temp_coe 
#endif
        else
           if(increase_equ_zonal_visc .and. abs(geolat(row))<5.*rad) then
              aux1=-Ah*(fac_visc_increase*dx(q)*udx + dy(q)*udy)*vol
              aux2=-Ah*(fac_visc_increase*dx(q)*vdx + dy(q)*vdy)*vol
           else
              aux1=-Ah*(dx(q)*udx + dy(q)*udy)*vol
              aux2=-Ah*(dx(q)*vdx + dy(q)*vdy)*vol 
           end if   
#ifdef use_non_hydrostatic
           aux3=-Ah*(dx(q)*wdx + dy(q)*wdy)*vol
#endif
        endif
        !
        if(use_cori_semi) then
           uv_rhs(row)=uv_rhs(row)+aux1*beta+aux2*gamma
           uv_rhs(row2)=uv_rhs(row2)+aux2*beta-aux1*gamma
#ifdef use_non_hydrostatic
           row3=row2+n3                            
           uv_rhs(row3)=uv_rhs(row3)+aux3*dt
#endif
        else
           uv_rhs(row)=uv_rhs(row)+aux1
           uv_rhs(row2)=uv_rhs(row2)+aux2
#ifdef use_non_hydrostatic
           row3=row2+n3                            
           uv_rhs(row3)=uv_rhs(row3)+aux3
#endif
        endif
     end do
  end do
  !
end subroutine biharmonic_viscosity
!
!--------------------------------------------------------------------------
subroutine ts_sfc_bc
  ! assemble the tracer surface boundary condition
  ! In both the implicit and explicit vertical mixing schemes cases,
  ! this routine should be called to set up the necessary arrays. 
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer         :: row, m, elnodes2(3), q, col, elnodes(3)
  real(kind=8)    :: auxf, entries(6), rsss

  if(.not.ts_surfbd) return

  ts_sfc_force=0.0

  ! sss restoring flux
  call cal_relax_salt

  ! virtual salt flux
#ifndef use_fullfreesurf
  rsss=ref_sss
  do row=1,ToDim_nod2d
     m=nod3D_below_nod2D(1,row)	    
     if(ref_sss_local) rsss=tracer(m,2)
     virtual_salt(row)=rsss*water_flux(row) 
  end do
#endif

  ! normalize the salt fluxes: remove the global mean
  if(balance_salt_water) call check_imb_salt_flux


  !apply fluxes and restoring to the RHS
  do row=1,myDim_elem2d             
     elnodes2=elem2D_nodes(:,row)
     elnodes=nod3D_below_nod2D(1,elnodes2)         
     auxf=voltriangle(row)/12.0_8    

#ifndef use_fullfreesurf   

#ifdef use_cavity
     if(all(cavity_flag_nod2d(elnodes2)==0)) then   
#endif   
        entries(1:3)=auxf*(restore_t_surf*(-tracer(elnodes,1)+Tsurf(elnodes2)) &
             - heat_flux(elnodes2)/vcpw)            
        entries(4:6)=auxf*(relax_salt(elnodes2)+virtual_salt(elnodes2))
#ifdef use_cavity
     else
        entries(1:3)=-auxf*heat_flux(elnodes2)/vcpw
        entries(4:6)=auxf*(relax_salt(elnodes2)+virtual_salt(elnodes2)) 
     end if
#endif

#else
     ! then full free surface

#ifdef use_cavity
     if(all(cavity_flag_nod2d(elnodes2)==0)) then   
#endif     
        entries(1:3)=auxf*(restore_t_surf*(-tracer(elnodes,1)+Tsurf(elnodes2)) &
             - heat_flux(elnodes2)/vcpw - tracer(elnodes,1)*water_flux(elnodes2))  !to update
        entries(4:6)=auxf*(relax_salt(elnodes2)+real_salt_flux(elnodes2))
#ifdef use_cavity
     else
        entries(1:3)=-auxf*(heat_flux(elnodes2)/vcpw + tracer(elnodes,1)*water_flux(elnodes2)) !to update
        entries(4:6)=auxf*relax_salt(elnodes2)  !only due to global normalization
     end if
#endif
#endif

     do q=1,3
        col=elnodes2(q) !elnodes2
        ts_sfc_force(col,1)=ts_sfc_force(col,1)+sum(entries(1:3))+entries(q)
        ts_sfc_force(col,2)=ts_sfc_force(col,2)+sum(entries(4:6))+entries(q+3) 
     end do
     if(.not.use_vertdiff_impl) then
        do q=1,3
           col=elnodes(q)                               
           tracer_rhs(col,1)=tracer_rhs(col,1)+sum(entries(1:3))+entries(q)
           tracer_rhs(col,2)=tracer_rhs(col,2)+sum(entries(4:6))+entries(q+3)
        end do
     end if
  end do

end subroutine ts_sfc_bc
!
!============================================================================
!
subroutine cal_relax_salt
  ! prepare restoring surface salt flux
  use o_mesh
  use o_param
  use o_array
  use i_array
  use g_config
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer         :: row, m
  real(kind=8)    :: aux

#ifdef use_ice
  do row=1,ToDim_nod2d
     !if(a_ice(row)>1.e-2) then
     !   !relax_salt(row)=0.0
     !	m=nod3D_below_nod2D(1,row)
     !	aux=Ssurf(row)-tracer(m,2)
     !   relax_salt(row)=restore_s_surf*aux/4.8667  !50m/4year
     !else
        m=nod3D_below_nod2D(1,row)	
        aux=Ssurf(row)-tracer(m,2)
     !   aux=sign(min(abs(aux),0.5), aux) 
        relax_salt(row)=restore_s_surf*aux
     !end if
  end do
#else
  do row=1,ToDim_nod2d
     m=nod3D_below_nod2D(1,row)	
     aux=Ssurf(row)-tracer(m,2)
     relax_salt(row)=restore_s_surf*aux
  end do
#endif  

end subroutine cal_relax_salt
!
!============================================================================
!
#ifndef use_tracer_gls
subroutine tracer_rhs_tg
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_forcing_arrays
  use g_PARFE
  implicit none 

  integer                      :: i, j, row, elnodes(4)
  integer                      :: elem, elem2, lay, elem_type
  real(kind=8)                 :: Kh, Kv_el, aux, adv, coe, fc
  real(kind=8)                 :: inv2, inv4, inv5, inv20
  real(kind=8)                 :: dx(4), dy(4), dz(4), vol
  real(kind=8)                 :: uvel(4), vvel(4), u2d, v2d
  real(kind=8)                 :: usum, vsum, um, vm, wm
  real(kind=8)                 :: tr_elem(4,num_tracer), entries(4)
  real(kind=8)                 :: aux1, aux2, dif(4)
  real(kind=8)                 :: rotate_coe, temp_coe, temp_coe2
  real(kind=8)                 :: swr_conv
  !
  !variables used for Redi/GM scheme
  real(kind=8)                 :: K_GM, Kh_diag, Kh_other
  real(kind=8)	               :: S(3), fcn1, fcn2, lambd, Z_mean, depth_scale
  real(kind=8)                 :: S_d, c_speed, turb_layer
  data S_d, c_speed /1.0e-3, 2.0/
  !
  real(kind=8)                 :: dparam, beta, gamma
#ifndef use_non_hydrostatic
  integer                      :: elnodes2(3) 
  real(kind=8)                 :: vc(3)
#else
  integer                      :: elnodes23(4)
  real(kind=8)                 :: vc(4), wsum, wvel(4), w2d
#endif
  !
#ifdef use_fullfreesurf
  integer                      :: n_el
  real(kind=8)                 :: coe_diff, wmove(4), wmove_sum
  logical                      :: flag_move=.false.
#endif

  inv2=0.5_8
  inv4=0.25_8
  inv5=0.2_8
  inv20=0.05_8

  do j=1,num_tracer
     do row=1, myDim_nod3d          
        tracer_rhs(row,j)=0.
     enddo
  enddo

  ! assembling tracer rhs
  do elem=1,myDim_elem3d              
     elnodes=elem3D_nodes(:,elem)
     elem2=elem2D_corresp_to_elem3D(elem)
     elem_type=grid_type_elem2d(elem2)
#ifndef use_non_hydrostatic
     elnodes2=elem2d_nodes(:,elem2)
#else
     elnodes23=nod2d_corresp_to_nod3d(elnodes)
#endif
     fc=coriolis_param_elem2d(elem2)

     !diffusivity
     Kh=Kh0
     if(scale_mixing_h) then
        Kh=Kh*(voltriangle(elem2)/scalevol)**(1.0/real(scale_mixing_type))
     end if
     if(use_vertdiff_impl) then
        Kv_el=0.0
     else
        Kv_el=sum(Kv(elnodes,1))/4.0 
     end if

     !-------------------------------------------------------------------------------
     if (Redi_GM .and. elem_type==0) then 

        !GM diffusivity
        K_GM  = Kh*ratio_K_GM   

        !neutral slope
        if(nslope_version==1) then
           lay=elem3d_layer(elem)
           S = neutral_slope(:,lay,elem2)   ! S(1:3): Sx,Sy and |S|
        else
           S = neutral_slope_elem(:,elem)
        end if
        !prepare for tapering
        !define 2 functions fcn1 and fcn2, which are required for tapering
        fcn1=1.0_8
        fcn2=1.0_8

        ! fcn1, hyperbolic tangent, used for steep slope region
        if(ODM95) fcn1 = 0.5_8*(1.0_8 + tanh((Slope_critical - S(3))/S_d))

        !we need to check if the element is near the surface in the transit layer
        !If yes, then we need the tapering function fcn2, a sine function of depth.
        if(LDD97) then
           !the first baroclinic Rossby radius
           lambd = c_speed/abs(fc)

           !limit lambda [following Large et al(1997)] to handle singularity
           !at the equator
           if (lambd < 15000.) lambd = 15000.
           if (lambd > 100000.) lambd = 100000.

           !critical depth, above which sine tapering is necessary.
           depth_scale = lambd*S(3)

           !the mean depth of this element
           Z_mean = abs(sum(coord_nod3D(3,elnodes)))/4.0

           turb_layer=100.0
           !here 100.0m depth is assumed to be turbulent and does not need Redi/GM
           !a well defined turbulent layer should be updated later

           !if in the turbulent layer or transit layer, calculate fcn2
           if(Z_mean<turb_layer) then   !turbulent layer
              fcn2=0.0
           elseif(Z_mean < turb_layer+depth_scale)  then   !transit layer
              fcn2 = 0.5*(1.0 + sin(pi*(Z_mean-turb_layer)/depth_scale - pi/2.0))
           end if
        end if

        ! apply tapering
        ! For steep slope region:
        ! a) no taper applied to diagonal piece of horizontal neutral operator
        ! b) hyperbolic tangent(exponential) taper applied to off-diagonal piece of
        !    horizontal operator and to diagonal and off-diagonal piece of vertical
        !    neutral diffusion operator. a)+b) means we transfer the tracer diffusion
        !    to a horizontal-vertical manner in regions of steep neutral slopes.
        ! c) Exponential taper applied to GM operator.
        ! For surface layer with small slope:
        ! a) sine taper applied to both neutral operator and GM operator, except the
        !    diagonal piece of the horizontal diffusion.
        ! In one word, here we use ldd97, but always keep the diagonal part of the 
        ! horizontal diffusion following the suggestion of Griffies (2004).

        ! diffusion part:
        Kh_diag = Kh
        Kh_other = Kh*fcn1*fcn2
        ! skewion part:
        K_GM = K_GM*fcn1*fcn2

	if(any(mod(index_nod3d(elnodes),10)==1)) K_GM=0.0 !b.c.

     end if  	!Redi_GM:  tapered neutral diffusivity computed
     !------------------------------------------------------------------------------

     !derivatives
     dx=bafux_3d(:, elem)
     dy=bafuy_3d(:, elem)
     dz=bafuz_3d(:, elem)
     vol=voltetra(elem)
     !velocity
     uvel=uf0(elnodes)         !u*
     vvel=uf0(elnodes+myDim_nod3d+eDim_nod3D)   !v*    
     usum=sum(uvel)
     vsum=sum(vvel)
#ifndef use_non_hydrostatic
#ifdef use_semiimplicit_scheme
     vc=theta_ssh*ssh(elnodes2)-(gamma_stab-1.0+theta_ssh)*ssh0(elnodes2)
#else
     vc=ssh(elnodes2)-gamma_stab*ssh0(elnodes2)
#endif
     aux1=sum(vc*bafux_2d(:,elem2))*g
     aux2=sum(vc*bafuy_2d(:,elem2))*g
     if(use_cori_semi) then
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d	
     wm=sum(dz*w(elnodes))
#else
     !non_hydrostatic case
     wvel=uf0(elnodes+2*(myDim_nod3d+eDim_nod3D))  
     wsum=sum(wvel)
     vc=g*ssh(elnodes23)+nhp(elnodes) - &
          gamma_stab*g*ssh0(elnodes23)-gamma_stab_nh*nhp0(elnodes)
     aux1=sum(vc*dx)
     aux2=sum(vc*dy)
     if(use_cori_semi) then
        dparam=dt_inv**2 + alpha_trapez**2*fc**2
        beta=dt_inv/dparam
        gamma=fc*alpha_trapez/dparam
        u2d=beta*aux1+gamma*aux2
        v2d=beta*aux2-gamma*aux1
     else
        u2d=aux1*dt
        v2d=aux2*dt	
     endif
     w2d=sum(vc*dz)*dt
     um=inv4*usum-u2d
     vm=inv4*vsum-v2d
     wm=inv4*wsum-w2d
#endif

     !element nodal tracer fields
     tr_elem=tracer(elnodes,:)
     ! shortwave penetration
#ifdef use_sw_pene
     swr_conv=sum(dz*sw_3d(elnodes))*vol
#endif

     ! used for free surface case
#ifdef use_fullfreesurf
     coe=inv20*vol*dt_inv
     n_el=map_elem(elem)
     flag_move=.false.
     if(n_el/=0) then
        flag_move=.true.
        coe_diff=coe-inv20*voltetra_new(n_el)*dt_inv
        wmove=0.
        do i=1, 4
           if(index_nod3d(elnodes(i))>12) cycle
           row=nod2D_corresp_to_nod3D(elnodes(i))     
           wmove(i)=-(ssh(row)-ssh0(row))*dt_inv
        end do
        wmove_sum=sum(wmove)
     end if
#endif
     !for along sigma mixing
     if (elem_type==1) then
        lay=elem3d_layer(elem)     
        S(1)=grid_slope(1,lay,elem2)
        S(2)=grid_slope(2,lay,elem2)
        aux1=S(1)**2+S(2)**2
        S(3)=sqrt(aux1)
        rotate_coe=1.0/(1.0+aux1) 
     end if


     !assembling over nodes
     do i=1,4 
        row=elnodes(i)
        entries=0.0

        ! part of mass matrix due to mesh movement
#ifdef use_fullfreesurf
        if(flag_move) then
           do j=1,num_tracer
              tracer_rhs(row,j)=tracer_rhs(row,j)+(sum(tr_elem(:,j))+tr_elem(i,j))*coe_diff          
           end do
        end if
#endif

        !diffusion
        if(elem_type==1) then  !along-sigma diffusion 
	   !app.: vertical diffusion is kept in the vertical 
           !diagonal part1 (lateral)
	   temp_coe=Kh*rotate_coe
           entries=entries-(temp_coe*dx(i)*(1.0+S(2)*s(2))*dx + &
                temp_coe*dy(i)*(1.0+S(1)*S(1))*dy)
           !diagonal part2 (cross slope)
	   temp_coe2=Kv_el+S(3)*S(3)*temp_coe
           entries=entries-temp_coe2*dz(i)*dz
           !off diagonal part1 (lateral) --> (1,3),(2,3)
           entries=entries-temp_coe*(S(1)*dx(i)+S(2)*dy(i))*dz
           !off diagonal part1 (lateral)--> (1,2),(2,1), '- -' -> '+'
	   temp_coe2=S(1)*S(2)*temp_coe
           entries=entries+(temp_coe2*dx(i)*dy + temp_coe2*dy(i)*dx)
           !off diagonal part2 (cross slope) --> (3,1),(3,2)
           entries=entries-dz(i)*temp_coe*(S(1)*dx+S(2)*dy)
        else
	   if(Redi_GM) then  ! Redi diffusion + GM
       !app: small slope
       !diagonal part1 (lateral)
              entries=entries - Kh_diag*dx(i)*dx - Kh_diag*dy(i)*dy
              !diagonal part2 (cross neutral)
	      aux=(Kv_el+Kh_other*S(3)*S(3))*dz(i)
              entries=entries-aux*dz
              !off diagonal part1 (lateral) 
              aux=(Kh_other-K_GM)*(S(1)*dx(i)+S(2)*dy(i))
              entries=entries-aux*dz
              !off diagonal part2 (cross neutral)
              aux=(Kh_other+K_GM)*dz(i)
              entries=entries-(aux*S(1)*dx+aux*S(2)*dy)
	   else  !horizontal diffusion
              entries=entries - Kh*dx(i)*dx - Kh*dy(i)*dy - Kv_el*dz(i)*dz
	   end if
        end if

        !advection
#ifndef use_non_hydrostatic
#ifndef use_fullfreesurf
        entries=entries - ((usum+uvel(i))*inv5-u2d)*inv4*dx - &
             ((vsum+vvel(i))*inv5-v2d)*inv4*dy - wm*inv4*dz
#else
        entries=entries+(dx(i)*(usum+uvel)+dy(i)*(vsum+vvel))*inv20
        entries=entries+dz(i)*wm*inv4
        entries=entries-(dx(i)*u2d+dy(i)*v2d)*inv4
        if(flag_move) entries=entries+dz(i)*(wmove_sum+wmove)*inv20
#endif
#else
#ifndef use_fullfreesurf
        entries=entries - ((usum+uvel(i))*inv5-u2d)*inv4*dx - &
             ((vsum+vvel(i))*inv5-v2d)*inv4*dy - &
             ((wsum+wvel(i))*inv5-w2d)*inv4*dz
#else
        entries=entries+(dx(i)*(usum+uvel)+dy(i)*(vsum+vvel) + &
             dz(i)*(wsum+wvel))*inv20
        entries=entries-(dx(i)*u2d+dy(i)*v2d+dz(i)*w2d)*inv4
        if(flag_move) entries=entries+dz(i)*(wmove_sum+wmove)*inv20
#endif
#endif

        !TG stabilization
        aux=(um*dx(i)+vm*dy(i)+wm*dz(i))*dt*inv2   
#ifdef use_fullfreesurf
        if(flag_move) aux=aux+wmove_sum*inv4*dz(i)*dt*inv2 
#endif
        entries=entries-aux*(um*dx+vm*dy+wm*dz)
#ifdef use_fullfreesurf
        if(flag_move) then
           entries=entries-aux*wmove_sum*inv4*dz
        end if
#endif     

        !sum up
        do j=1,num_tracer
           tracer_rhs(row,j)=tracer_rhs(row,j)+sum(entries*tr_elem(:,j))*vol
        end do

        ! in case of considering shortwave penetration into the ocean
#ifdef use_sw_pene
        tracer_rhs(row,1)=tracer_rhs(row,1)+swr_conv*inv4 
#endif

        ! local or near open boundary restoring
        if(buffer_zone) then
           aux=tracer_restore_coeff(row)*vol*inv20
           do j=1,num_tracer  !!!
              dif=tracer0(elnodes,j)-tr_elem(:,j)
              tracer_rhs(row,j)=tracer_rhs(row,j)+aux*(sum(dif)+dif(i))
           end do
        end if

     end do  ! i node

     ! Add elementwise SGS velocity and fluxes to temporary arrays for diagnostics
     ! currently only fluxes of T/S are saved. Modify the code to save more tracers
#ifdef allow_diag
     if(diag_oce) then
        if(elem_type==1) then ! sigma grid
           if(diag_oce_SGS_transp) then
              aux1=sum(dx*tr_elem(:,1))
              aux2=sum(dy*tr_elem(:,1))
              aux=sum(dz*tr_elem(:,1))
              sgs_ut(elem)=sgs_ut(elem)-temp_coe*aux1*(1.0+S(2)*S(2)) &
                   -temp_coe*S(1)*aux+temp_coe2*aux2
              sgs_vt(elem)=sgs_vt(elem)-temp_coe*aux2*(1.0+S(1)*S(1)) &
                   -temp_coe*S(2)*aux+temp_coe2*aux1
              aux1=sum(dx*tr_elem(:,2))
              aux2=sum(dy*tr_elem(:,2))
              aux=sum(dz*tr_elem(:,2))
              sgs_us(elem)=sgs_us(elem)-temp_coe*aux1*(1.0+S(2)*S(2)) &
                   -temp_coe*S(1)*aux+temp_coe2*aux2
              sgs_vs(elem)=sgs_vs(elem)-temp_coe*aux2*(1.0+S(1)*S(1)) &
                   -temp_coe*S(2)*aux+temp_coe2*aux1
           endif
        else
           if(Redi_GM) then
              if(diag_oce_GM_vel) then
                 ! GM velocity
                 sgs_u(elem)=sgs_u(elem)-K_GM*S(1)   
                 sgs_v(elem)=sgs_v(elem)-K_GM*S(2)  
              end if
              ! flux
              if(diag_oce_SGS_transp) then
                 aux=sum(dz*tr_elem(:,1))
                 sgs_ut(elem)=sgs_ut(elem)-Kh_diag*sum(dx*tr_elem(:,1))-(Kh_other-K_GM)*S(1)*aux
                 sgs_vt(elem)=sgs_vt(elem)-Kh_diag*sum(dy*tr_elem(:,1))-(Kh_other-K_GM)*S(2)*aux
                 aux=sum(dz*tr_elem(:,2))
                 sgs_us(elem)=sgs_us(elem)-Kh_diag*sum(dx*tr_elem(:,2))-(Kh_other-K_GM)*S(1)*aux
                 sgs_vs(elem)=sgs_vs(elem)-Kh_diag*sum(dy*tr_elem(:,2))-(Kh_other-K_GM)*S(2)*aux
              end if
           else
              if(diag_oce_SGS_transp) then
                 sgs_ut(elem)=sgs_ut(elem)-Kh*sum(dx*tr_elem(:,1))
                 sgs_vt(elem)=sgs_vt(elem)-Kh*sum(dy*tr_elem(:,1))
                 sgs_us(elem)=sgs_us(elem)-Kh*sum(dx*tr_elem(:,2))
                 sgs_vs(elem)=sgs_vs(elem)-Kh*sum(dy*tr_elem(:,2))
              end if
           end if
        end if
     end if
#endif

  end do  ! element

end subroutine tracer_rhs_tg
#endif
!
!----------------------------------------------------------------------------
!
subroutine impl_vertdiff
  ! apply implicit vertical diffusivity
  ! lumped mass matrix is used
  use o_param
  use o_array
  use o_mesh
  use o_elements
  use o_mixing_kpp_mod
  use o_passive_tracer_mod
  use g_config
  use g_parfe
  implicit none

  integer            :: n2, nlay, k, i, j, q, p, cnt
  integer            :: row, elem, elnodes(4), n_el
  integer            :: nodup, nodlo
  real(kind=8)       :: vol, vol_visc, inv4, Kv_el, Kv_el2
  real(kind=8)       :: rhs_nonloc(2)
  real(kind=8)       :: a(max_num_layers,2), b(max_num_layers,2), c(max_num_layers,2)
  real(kind=8)       :: trr(num_tracer, max_num_layers), tr_nod(num_tracer)

  inv4=0.25_8

  do n2=1,myDim_nod2d           
     nlay=num_layers_below_nod2d(n2)+1

     ! assemble three diagonal matrix and rhs
     trr=0.
     do k=1,nlay
        a(k,:)=0.
        b(k,:)=0.
        c(k,:)=0.

        row=nod3d_below_nod2d(k,n2)
        tr_nod=tracer(row,:)

        if(k>1) nodup=nod3d_below_nod2d(k-1,n2)
        if(k<nlay) nodlo=nod3d_below_nod2d(k+1,n2)
        do i=1,nod_in_elem3D(row)%nmb
           elem=nod_in_elem3D(row)%addresses(i)
           elnodes=elem3D_nodes(:,elem)
           vol_visc=voltetra(elem)
#ifdef use_fullfreesurf
           n_el=map_elem(elem)
           if(n_el==0) then
              vol=voltetra(elem)
           else
              vol=voltetra_new(n_el)
           end if
#else
           vol=voltetra(elem)
#endif

           do q=1,4
              if(elnodes(q)==row) then
                 p=q
                 exit
              end if
           end do

           cnt=0
           rhs_nonloc=0.0

           !first entry
           if(k>1) then
              do q=1,4
                 if(elnodes(q)==nodup) then
                    Kv_el=Kv(nodup,1)  !Kv(nodup) is the Kv for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Kv_el=Kv_el+Kv0
                    end if
                    a(k,1)=a(k,1) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Kv_el * vol_visc
                    !second entry
                    b(k,1)=b(k,1) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Kv_el * vol_visc
                    if(trim(mix_scheme)=='KPP') then
                       if(double_diffusion) then
                          Kv_el2=Kv(nodup,2)
                          a(k,2)=a(k,2) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Kv_el2 * vol_visc
                          !second entry
                          b(k,2)=b(k,2) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Kv_el2 * vol_visc
                       else
                          Kv_el2=Kv_el
                       end if
                       ! nonlocal transport to the rhs (only T and S currently)
                       rhs_nonloc(1)=bafuz_3d(p,elem)*vol_visc*min(Kv_el*ghats(nodup),1.)*heat_flux(n2)/vcpw
                       rhs_nonloc(2)=-bafuz_3d(p,elem)*vol_visc*min(Kv_el2*ghats(nodup),1.)*tracer(n2,2)*water_flux(n2)
                    end if

                    cnt=1
                    exit
                 end if
              end do
           end if

           !third entry
           if(k<nlay) then
              do q=1,4
                 if(elnodes(q)==nodlo) then
                    Kv_el=Kv(row,1)    !Kv(row) is the Kv for this segment
                    if(trim(mix_scheme)=='MY2p5') then
                       Kv_el=Kv_el+Kv0
                    end if
                    c(k,1)=c(k,1) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Kv_el * vol_visc
                    !second entry
                    b(k,1)=b(k,1) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Kv_el * vol_visc
                    if(trim(mix_scheme)=='KPP') then
                       if(double_diffusion) then
                          Kv_el2=Kv(row,2)
                          c(k,2)=c(k,2) + bafuz_3d(p,elem) * bafuz_3d(q,elem) * Kv_el2 * vol_visc
                          !second entry
                          b(k,2)=b(k,2) + bafuz_3d(p,elem) * bafuz_3d(p,elem) * Kv_el2 * vol_visc
                       else
                          Kv_el2=Kv_el
                       end if
                       ! nonlocal transport to the rhs (only T and S currently)
                       rhs_nonloc(1)=bafuz_3d(p,elem)*vol_visc*min(Kv_el*ghats(row),1.)*heat_flux(n2)/vcpw
                       rhs_nonloc(2)=-bafuz_3d(p,elem)*vol_visc*min(Kv_el2*ghats(row),1.)*tracer(n2,2)*water_flux(n2)
                    end if

                    cnt=1
                    exit
                 end if
              end do
           end if

           !second entry, mass matrix part
           b(k,:)=b(k,:) + vol*inv4*dt_inv

           !rhs
           do j=1,num_tracer
              trr(j,k)=trr(j,k)+tr_nod(j)*vol*inv4*dt_inv
           end do
           !add non-local term
           trr(1,k)=trr(1,k)+rhs_nonloc(1)
           trr(2,k)=trr(2,k)+rhs_nonloc(2)          

        end do  !i
     end do  !k

     !add surface forcing to rhs
     trr(1:2,1)=trr(1:2,1)+ts_sfc_force(n2,1:2)  !T&S
     if(use_passive_tracer .and. passive_tracer_flux) then
        do j=1, num_passive_tracer
           trr(index_passive_tracer(j),1)=trr(index_passive_tracer(j),1) + ptr_sfc_force(n2,j)
        end do
     end if

     ! the sweep algorithm
     if(trim(mix_scheme)=='KPP' .and. double_diffusion) then
        ! forward step
        ! prepare coefficients 
        c(1,:)=-c(1,:)/b(1,:)
        trr(1,1)=trr(1,1)/b(1,1)
        trr(2:,1)=trr(2:,1)/b(1,2)
        do k=2,nlay
           b(k,:)=1.0/(a(k,:)*c(k-1,:)+b(k,:))
           c(k,:)=-c(k,:)*b(k,:)
           trr(1,k)=(trr(1,k)-a(k,1)*trr(1,k-1))*b(k,1)
           trr(2:,k)=(trr(2:,k)-a(k,2)*trr(2:,k-1))*b(k,2)
        end do
        ! backward sweep
        do k=nlay-1,1,-1
           trr(1,k)=trr(1,k)+c(k,1)*trr(1,k+1)
           trr(2:,k)=trr(2:,k)+c(k,2)*trr(2:,k+1)
        end do
     else
        ! forward step
        ! prepare coefficients 
        c(1,1)=-c(1,1)/b(1,1)
        trr(:,1)=trr(:,1)/b(1,1)
        do k=2,nlay
           b(k,1)=1.0/(a(k,1)*c(k-1,1)+b(k,1))
           c(k,1)=-c(k,1)*b(k,1)
           trr(:,k)=(trr(:,k)-a(k,1)*trr(:,k-1))*b(k,1)
        end do
        ! backward sweep
        do k=nlay-1,1,-1
           trr(:,k)=trr(:,k)+c(k,1)*trr(:,k+1)
        end do
     end if

     ! update tracers
     do k=1,nlay
        row=nod3d_below_nod2d(k,n2)
        tracer(row,:)=trr(:,k)
     end do

  end do  !2d nodes

  do j=1,num_tracer
     call com_3D(tracer(:,j))
  end do

end subroutine impl_vertdiff
!
!----------------------------------------------------------------------------
!
subroutine tracer_solve
  use o_matrices
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                                 :: n, clo, clo2, cn, location(100), row, j
  real(kind=8)                            :: rhs_new
  real(kind=8), allocatable               :: auxarray(:,:)

  allocate(auxarray(myDim_nod3d,num_tracer))

  !the first approximation
  do j=1,num_tracer
     do row=1,myDim_nod3d              
        dtracer(row,j)=tracer_rhs(row,j)/ts_lump(row)
     end do
     call com_3D(dtracer(:,j))
  end do

  !iterate 
  do n=1,num_iter_solve-1                  
     do row=1,myDim_nod3d           
        clo=tsstiff%rowptr(row)-tsstiff%rowptr(1)+1  
        clo2=tsstiff%rowptr(row+1)-tsstiff%rowptr(1)  
        cn=clo2-clo+1
        location(1:cn)=nghbr_nod3D(row)%addresses    
        do j=1,num_tracer
           rhs_new=tracer_rhs(row,j) - sum(tsstiff%values(clo:clo2)*dtracer(location(1:cn),j))
           auxarray(row,j)=dtracer(row,j)+rhs_new/ts_lump(row)  
        end do
     end do

     do j=1,num_tracer 
        do row=1,myDim_nod3d            
           dtracer(row,j)=auxarray(row,j)  
        end do
        call com_3D(dtracer(:,j))
     end do
  end do

  deallocate(auxarray)
end subroutine tracer_solve
!
!----------------------------------------------------------------------------
!
subroutine compute_neutral_slope
  ! calculate neutral slopes 
  !
  ! OUTPUT:
  ! version 1: neutral_slope(3,layer,elem2d), over prism
  ! version 2: neutral_slope_elem(3, elem3d), elementwise neutral slope   
  ! 
  ! REFERENCE:
  !    McDougall, T.J. and  D.R. Jackett, 1988  
  !    On the helical nature of neutral surfaces
  !    Progress in Oceanography, vol 20, Pergamon, 153-183
  !    (about the definition of normal direction of neutral surface)
  !
  !    Griffies, S.M. et al. 1998
  !    Isoneutral Diffusion in a z-coordinate ocean model
  !    JPO, vol 28, 805-830
  !
  ! Qiang, 25,11,2004
  ! Qiang, revised on 02,07,2010
  !-------------------------------------------------------------------
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer :: i, k, nod, elem, lay_m, elnodes2(3), elnodes(6), lay_el(3)
  integer :: elnodes_tetra(4)
  real    :: rho(6), zm, dz, dx2d(3), dy2d(3), denom, ano_1, ano_2
  real    :: rhograd_x(max_num_layers), rhograd_y(max_num_layers), rhograd_z(max_num_layers)
  real    :: t, s, p, sw_alpha, sw_beta, g_t_x, g_t_y, g_t_z, g_s_x, g_s_y, g_s_z
  real    :: inv2, inv3, inv4, inv6

  inv2=0.5
  inv3=1.0/3.0
  inv4=0.25
  inv6=1.0/6.0

  if(nslope_version==1) then
     do elem=1, myDim_elem2D                  
        if(grid_type_elem2d(elem)==1) cycle  !only apply GM to z-level grids
        elnodes2=elem2d_nodes(:,elem)
        dx2d=bafux_2d(:,elem)
        dy2d=bafuy_2d(:,elem)
        lay_el=num_layers_below_nod2d(elnodes2)
        lay_m=minval(lay_el)

        do k=1,lay_m
           elnodes(1:3)=nod3d_below_nod2d(k,elnodes2)
           elnodes(4:6)=nod3d_below_nod2d(k+1,elnodes2)
           zm=sum(coord_nod3d(3,elnodes))*inv6
           dz=1.0/(coord_nod3d(3,elnodes(1))-coord_nod3d(3,elnodes(4)))
           do i=1,6
              nod=elnodes(i)
              call fcn_density(tracer(nod,1),tracer(nod,2),zm,rho(i))
           end do
           rhograd_x(k)=(sum(dx2d*rho(1:3))+sum(dx2d*rho(4:6)))*inv2
           rhograd_y(k)=(sum(dy2d*rho(1:3))+sum(dy2d*rho(4:6)))*inv2
           rhograd_z(k)=(sum(rho(1:3))-sum(rho(4:6)))*inv3*dz
        end do

        do k=1,lay_m
           if(lay_m>1) then
              if(k==1) then
                 denom=(2.0*rhograd_z(1)+rhograd_z(2))*inv3
              elseif(k<lay_m) then
                 denom=(rhograd_z(k-1)+2.0*rhograd_z(k)+rhograd_z(k+1))*inv4
              else
                 denom=(2.0*rhograd_z(lay_m)+rhograd_z(lay_m-1))*inv3
              end if
           else
              denom=rhograd_z(1)
           end if
           if(denom<0.0) then
              denom=denom-1.0e-20
              neutral_slope(1,k,elem)=-rhograd_x(k)/denom
              neutral_slope(2,k,elem)=-rhograd_y(k)/denom
              neutral_slope(3,k,elem) = sqrt(neutral_slope(1,k,elem)**2 + neutral_slope(2,k,elem)**2)
           else
              neutral_slope(:,k,elem)=0.0
           end if

        end do
     end do

  else !version 2

     uv_rhs(1:ToDim_nod3d)=0.   ! to save memory we use uv_rhs as a temporary array

     do elem = 1, myDim_elem3d
        elnodes_tetra=elem3D_nodes(:,elem)

        ! alpha and beta
        t=inv4*sum(tracer(elnodes_tetra,1))
        s=inv4*sum(tracer(elnodes_tetra,2))
        p=abs(inv4*sum(coord_nod3D(3,elnodes_tetra))) 
        call sw_alpha_beta(t,s,p,sw_alpha,sw_beta)

        g_t_z = sum(bafuz_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_z = sum(bafuz_3d(:,elem)*tracer(elnodes_tetra,2))
        denom = -sw_alpha*g_t_z + sw_beta*g_s_z

        if(denom>=0.0) then
           uv_rhs(elnodes_tetra)=1.
           neutral_slope_elem(:,elem)=0.0
           cycle
        end if

        denom=denom-1.e-20

        g_t_x = sum(bafux_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_x = sum(bafux_3d(:,elem)*tracer(elnodes_tetra,2))
        g_t_y = sum(bafuy_3d(:,elem)*tracer(elnodes_tetra,1))
        g_s_y = sum(bafuy_3d(:,elem)*tracer(elnodes_tetra,2))  

        ano_1 = -sw_alpha*g_t_x + sw_beta*g_s_x
        ano_2 = -sw_alpha*g_t_y + sw_beta*g_s_y

        neutral_slope_elem(1,elem) = -ano_1/denom
        neutral_slope_elem(2,elem) = -ano_2/denom
        neutral_slope_elem(3,elem) = &
             sqrt(neutral_slope_elem(1,elem)**2 + neutral_slope_elem(2,elem)**2)
     end do

     ! in case of static instability
     do i=1,myDim_nod3D
        if(uv_rhs(i)>0.5) then
           neutral_slope_elem(:,nod_in_elem3d(i)%addresses)=0.0
        end if
     end do

  end if
end subroutine compute_neutral_slope
!
!----------------------------------------------------------------------------
!
subroutine sw_alpha_beta(t1, s1, p1, alpha, beta)
  !   A function to calculate the thermal expansion coefficient
  !   and saline contraction coefficient.
  !
  ! INPUT: t1 (c), s1 (psu), p1 (db)
  !
  ! OUTPUT:
  !    alpha = Thermal expansion coeff (alpha) [degree_C.^-1]
  !    beta  = Saline contraction coeff (beta) [psu.^-1]
  !
  ! REFERENCE:
  !    McDougall, T.J. 1987.  Neutral Surfaces
  !    Journal of Physical Oceanography, vol 17, 1950-1964,
  !
  ! CHECK VALUE:
  !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
  !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
  !
  ! Qiang Wang, 25,11,2004
  !-----------------------------------------------------------------

  implicit none

  real :: t1, s1, p1, alpha, beta
  real :: t1_2,t1_3,t1_4,p1_2,p1_3,s35,s35_2 
  real :: a_over_b    

  t1_2 = t1*t1
  t1_3 = t1_2*t1
  t1_4 = t1_3*t1
  p1_2 = p1*p1
  p1_3 = p1_2*p1
  s35 = s1-35.0_8
  s35_2 = s35*s35

  ! calculate beta
  beta = 0.785567e-3 - 0.301985e-5*t1 &
       + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
       + s35*(-0.356603e-6 + 0.788212e-8*t1 &
       + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
       + s35_2*(0.515032e-8) & 
       + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
       + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
       + p1_3*(0.121551e-17)

  ! calaculate the thermal expansion / saline contraction ratio
  a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
       - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
       - 0.255019e-7*t1_4 &
       + s35*(0.378110e-2 - 0.846960e-4*t1 &
       - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
       + s35_2*(-0.678662e-5) &
       + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
       + p1_2*t1_2*(0.512857e-12) &
       - p1_3*(0.302285e-13)

  ! calculate alpha
  alpha = a_over_b*beta

end subroutine sw_alpha_beta
!
! This file collect subroutines implementing FE-FCT
! advection scheme by Loehner et al.
! There is a tunable paremeter gamma_fct in ts_solve_low_order and fem_fct.
! Increasing it leads to positivity preserving solution.
!
!----------------------------------------------------------------------------
!
subroutine fct_init
  use o_param
  use o_mesh
  use o_array
  use g_config
  use g_parfe
  !
  allocate(tral(myDim_nod3D+eDim_nod3D,num_tracer))       
  allocate(trafluxes(4,myDim_elem3D))
  allocate(pplus(myDim_nod3D+eDim_nod3D), pminus(myDim_nod3D+eDim_nod3D)) 
  tral=0.0

end subroutine fct_init
!
!----------------------------------------------------------------------------
!
subroutine fct_tracer_solve
  use o_solver
  use o_param
  use o_array
  use g_config
  use g_parfe

  integer     :: i

  ! high order solution
  if(lump_ts_matrix) then
     call tracer_solve
  else
     do i=1,num_tracer
        call solve(solve_tra+i-1)
     end do
  endif 

  ! low order solution
  call tra_solve_low_order

  ! fct
  do i=1,num_tracer
     call fem_fct(i)
  end do
  do i=1,num_tracer
     call com_3d(tracer(:,i))
  end do

end subroutine fct_tracer_solve
!
!----------------------------------------------------------------------------
!
subroutine tra_solve_low_order
  ! Low-order solution
  ! One adds diffusive contribution to the rhs. It is realized as
  ! difference between the consistent and lumped mass matrices
  ! acting on the field from the previous time step.   
  !
  use o_matrices
  use o_mesh
  use o_array
  use o_param
  use g_config
  use g_parfe
  implicit none
  
  integer      ::  j, row, clo, clo2, cn, location(100)
  
  do row=1,myDim_nod3D               
     clo=tsstiff%rowptr(row)-tsstiff%rowptr(1)+1  
     clo2=tsstiff%rowptr(row+1)-tsstiff%rowptr(1) 
     cn=clo2-clo+1
     location(1:cn)=nghbr_nod3D(row)%addresses    
     do j=1,num_tracer
        tral(row,j)=(tracer_rhs(row,j)+gamma_fct*sum(tsstiff%values(clo:clo2)* &
             tracer(location(1:cn),j)))/ts_lump(row) + &
             (1.-gamma_fct)*tracer(row,j)
     end do
   end do
  
   do j=1,num_tracer
      call com_3D(tral(:,j))     ! solution must be known to neighbours
   end do
end subroutine tra_solve_low_order
!
!----------------------------------------------------------------------------
!
subroutine fem_fct(tra_id)
  ! Flux corrected transport algorithm for tracer advection
  !
  ! It is based on Loehner et al. (Finite-element flux-corrected 
  ! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
  ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
  ! Turek. (kuzmin@math.uni-dortmund.de) 
  !
  ! Steps:
  !
  ! Construct a low-order solution by adding an artificial diffusion 
  ! to the rhs and lumping the mass matrix
  !
  ! Compute a high-order solution (run ts_solve)
  ! 
  ! Limit antidiffusive fluxes 
  ! 
  ! Update the low-order solution to the high order but using 
  ! the limited fluxes
  !
  use o_matrices
  use o_mesh
  use o_elements
  use o_array
  use o_param
  use g_config
  use g_parfe
  implicit none

  integer        :: tra_id
  integer        :: i, n, q, row
  integer        :: elem, elnodes(4)
  real(kind=8)   :: flux, ae, vol, aux(4), icoef(4,4), inv20
  real(kind=8), allocatable :: tmax(:), tmin(:) 

  allocate(tmax(myDim_nod3D), tmin(myDim_nod3D))
  inv20=1.0_8/20.0_8

  !==========================
  ! Compute elemental antidiffusive fluxes to nodes
  !==========================
  ! This is the most unpleasant part: it takes memory and time. 
  ! For every element we need its antidiffusive contribution to 
  ! each of its 4 nodes
  !
  ! Auxiliary elemental operator (lumped mass matrix - mass matrix)
  icoef=-inv20
  do n=1,4   
     icoef(n,n)=3.0*inv20
  end do
  ! antidiffusive fluxes 
  do elem=1, myDim_elem3D            
     elnodes=elem3D_nodes(:,elem)
     vol=voltetra(elem)*dt_inv   
     aux=gamma_fct*tracer(elnodes,tra_id) + dtracer(elnodes,tra_id)
     do q=1,4       
        trafluxes(q,elem)=sum(icoef(:,q)*aux)*vol/ts_lump(elnodes(q))
     end do
  end do

  !==========================   
  ! Screening the low-order solution
  !==========================
  ! TO BE ADDED IF FOUND NECESSARY
  ! Screening means comparing low-order solutions with the
  ! solution on the previous time step and using whichever 
  ! is greater/smaller in computations of max/min below
  !==========================
  ! Cluster min/max
  !==========================
  do row=1, myDim_nod3D               
     n=nghbr_nod3D(row)%nmb
     !-----------------------
     tmax(row)=maxval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id))  
     tmin(row)=minval(tral(nghbr_nod3D(row)%addresses(1:n),tra_id))  
     !-----------------------
     ! Admissible increments
     tmax(row)=tmax(row)-tral(row,tra_id)    
     tmin(row)=tmin(row)-tral(row,tra_id)   
  end do

  !=========================
  ! Sums of positive/negative fluxes to node row
  !=========================
  pplus=0.
  pminus=0.
  do elem=1, myDim_elem3D           
     elnodes=elem3D_nodes(:,elem)
     do q=1,4
        n=elnodes(q) 
        flux=trafluxes(q,elem)    
        if (flux>0.) then
           pplus(n)=pplus(n)+flux
        else
           pminus(n)=pminus(n)+flux	  
        end if
     end do
  end do

  !========================
  ! The least upper bound for the correction factors
  !========================
  do n=1,myDim_nod3D                
     flux=pplus(n)
     if (abs(flux)>0.) then    !!!!
        pplus(n)=min(1.0,tmax(n)/flux)  
     else
        pplus(n)=0.
     end if
     flux=pminus(n)
     if (abs(flux)>0.) then    !!!! 
        pminus(n)=min(1.0,tmin(n)/flux)  
     else
        pminus(n)=0.
     end if
  end do
  ! pminus and pplus are to be known to neighbouring PE
  call com_3D(pplus)
  call com_3D(pminus)	   

  !========================	 
  ! Limiting
  !========================	 
  do elem=1, myDim_elem3D           
     elnodes=elem3D_nodes(:,elem)
     ae=1.0
     do q=1,4
        n=elnodes(q)  
        flux=trafluxes(q,elem)            
        if(flux>=0.) ae=min(ae,pplus(n))
        if(flux<0.) ae=min(ae,pminus(n))
     end do
     trafluxes(:,elem)=ae*trafluxes(:,elem) 
  end do

  !==========================
  ! Update the solution
  !==========================
  do n=1,myDim_nod3D                 
     tracer(n,tra_id)=tral(n,tra_id)
  end do
  do elem=1, myDim_elem3D             
     elnodes=elem3D_nodes(:,elem)
     do q=1,4
        n=elnodes(q)  
        tracer(n,tra_id)=tracer(n,tra_id)+trafluxes(q,elem) 
     end do
  end do

  deallocate(tmin, tmax)
end subroutine fem_fct
!
!----------------------------------------------------------------------------
!
subroutine convect_adjust
  ! Increase vertical mixing in case of static instability
  ! It is considered (included) by PP scheme already.
  ! Only required for special cases
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use o_array
  use g_config
  use g_PARFE
  implicit none
  !
  integer     	:: col,n,j,node_low,node_up,node
  real(kind=8)	:: density_up, drho

  do col=1,myDim_nod2d+eDim_nod2d    

     n=num_layers_below_nod2D(col)

     do j=1,n       
        node=nod3D_below_nod2D(j,col)
        node_up = node
        node_low = nod3D_below_nod2D(j+1,col)

        call fcn_density(tracer(node_up,1),tracer(node_up,2), &
             coord_nod3d(3,node_low),density_up)

        drho = density_up - density_insitu(node_low)

        if(drho>0.0) then
           if(allow_convect_global) then
              Kv(node,1)=0.005
              Av(node)=0.1
           elseif(coord_nod2d(2,col)>0.0) then     
              Kv(node,1)=0.005     
              Av(node)=0.1
           endif
        else
           Kv(node,1)=0.0
           Av(node)=0.0
        end if

     end do
  end do
end subroutine convect_adjust

!example routines for initializing buffer zone restoring and 
!open boundary velocity restoring.
!----------------------------------------------------------------

subroutine init_restoring_bufferzone
  ! init tracer buffer zone for FO001 setup (a southern ocean setup)
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer                     :: row
  real(kind=8)                :: y, d, buffer_dist
  real(kind=8)                :: ymax_local, ymax_global

  buffer_dist=150.0e3  !in m, buffer zone scale

  ymax_local=maxval(coord_nod2d(2,:))

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE(ymax_local, ymax_global, &
       1, MPI_DOUBLE_PRECISION, MPI_MAX, &
       MPI_COMM_WORLD, MPIerr)

  allocate(tracer_restore_coeff(myDim_nod3d+eDim_nod3d))
  tracer_restore_coeff=0.0

  do row=1,myDim_nod3d+eDim_nod3d
     
     y=coord_nod3d(2,row)
     d=(ymax_global-y)*r_earth
     
     if(d<buffer_dist) then
        tracer_restore_coeff(row)=(buffer_dist-d)/buffer_dist*restore_ts_buff
     end if
  end do

  if(mype==0) write(*,*) 'restoring buffer zone ready: FO001'

end subroutine init_restoring_bufferzone
!
!--------------------------------------------------------------------
!
subroutine init_restoring_bufferzone_FOt
  ! init tracer buffer zone for FOt setup
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  
  integer                     :: row, n
  real(kind=8)                :: x, y, d, d_min, buffer_dist

  integer, allocatable        :: ind_op(:), ind_op_glo(:)
  real(kind=8), allocatable   :: x_op(:), y_op(:)
  real(kind=8), allocatable   :: x_op_glo(:), y_op_glo(:)

  buffer_dist=60.0e3  !in m, buffer zone scale

  allocate(ind_op(nod2d), x_op(nod2d), y_op(nod2d))
  ind_op=0
  x_op=0.0
  y_op=0.0
  allocate(ind_op_glo(nod2d), x_op_glo(nod2d), y_op_glo(nod2d))
  ind_op_glo=0
  x_op_glo=0.0
  y_op_glo=0.0

  do row=1,myDim_nod2d
     if(index_nod3d(nod3d_below_nod2d(1,row))==12) then
        n=myList_nod2d(row)
        ind_op(n)=1
        x_op(n)=coord_nod2d(1,row)
        y_op(n)=coord_nod2d(2,row)
     end if
  end do

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE(ind_op, ind_op_glo, &
       nod2d, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(x_op, x_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(y_op, y_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  deallocate(ind_op, x_op, y_op)

  allocate(tracer_restore_coeff(myDim_nod3d+eDim_nod3d))
  tracer_restore_coeff=0.0

  do row=1,myDim_nod3d+eDim_nod3d
     
     x=coord_nod3d(1,row)
     y=coord_nod3d(2,row)
     
     d_min=1000.0e3   !dist. in m
     do n=1,nod2d
        if(ind_op_glo(n)/=1) cycle
        call dist_on_earth(x, y, x_op_glo(n), y_op_glo(n), d)
        d_min=min(d_min, d)
     end do

     if(d_min<buffer_dist) then
        tracer_restore_coeff(row)=(buffer_dist-d_min)/buffer_dist*restore_ts_buff
     end if
  end do

  deallocate(ind_op_glo, x_op_glo, y_op_glo)

  if(mype==0) write(*,*) 'restoring buffer zone ready: FOt'

end subroutine init_restoring_bufferzone_FOt
!
!--------------------------------------------------------------------
!
subroutine init_restoring_bufferzone_Arctic
  ! init T/S buffer zone for Arctic setup
  ! T/S restored to PHC2.1 annual mean.
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer                     :: i, j, n
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  real(kind=8)                :: pp, pr, tt, ss, lon, lat
  real(kind=8)                :: rest_bound, rest_range
  real(kind=8), external      :: theta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data(:,:,:)
  real(kind=8), allocatable   :: temp_x(:), temp_y(:)


  !1) T S fields for restoring: T_0, S_0

  ! open global T/S data files
  !open(19,file=trim(ClimateDataPath)//'Winter_phc2.1_beta_ts.out', status='old')
  open(19,file=trim(ClimateDataPath)//'annual_phc_ts.out', status='old')
  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data(num_lon_reg,num_lat_reg,num_lay_reg))

  ! model grid coordinates
  allocate(temp_x(myDim_nod3d+eDim_nod3D), temp_y(myDim_nod3d+eDim_nod3D))
  do n=1, myDim_nod3d+eDim_nod3D        
     if(rotated_grid) then
        call r2g(lon, lat, coord_nod3d(1,n), coord_nod3d(2,n))
        temp_x(n)=lon/rad   ! change unit to degree
        temp_y(n)=lat/rad
     else
        temp_x(n)=coord_nod3d(1,n)/rad   
        temp_y(n)=coord_nod3d(2,n)/rad
     end if
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
  end do

  ! read raw data and do interpolation
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, nod3d, temp_x, temp_y, coord_nod3d(3,:), tracer0(:,1))

  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, nod3d, temp_x, temp_y, coord_nod3d(3,:), tracer0(:,2))

  close(19) 

  ! Convert in situ temperature into potential temperature
  pr=0.0_8
  do i=1,myDim_nod3d+eDim_nod3D                       
     tt=tracer0(i,1)
     ss=tracer0(i,2)
     pp=abs(coord_nod3D(3,i))
     tracer0(i,1)=theta(ss, tt, pp, pr)
  end do

  deallocate(temp_y, temp_x, raw_data, lay_reg, lat_reg, lon_reg)


  !2) where to apply restoring for T/S

  allocate(tracer_restore_coeff(myDim_nod3D+eDim_nod3D)) 
  tracer_restore_coeff=0.0
  ! east boundary
  rest_bound=maxval(coord_nod3d(1,:))
  rest_range=1.0*rad
  do i=1, myDim_nod3D+eDim_nod3D         
     if(coord_nod3D(1,i) >= rest_bound-rest_range) then
        tracer_restore_coeff(i)=1.0-abs(rest_bound - coord_nod3D(1,i))/rest_range
     end if
  end do
  ! west boundary
  rest_bound=minval(coord_nod3d(1,:))
  rest_range=3.0*rad
  do i=1, myDim_nod3D+eDim_nod3D        
     if(coord_nod3D(1,i) <= rest_bound+rest_range) then
        tracer_restore_coeff(i)=1.0-abs(rest_bound - coord_nod3D(1,i))/rest_range
     end if
  end do
  tracer_restore_coeff=tracer_restore_coeff*restore_ts_buff

  if(mype==0) write(*,*) 'restoring buffer zone ready'

end subroutine init_restoring_bufferzone_Arctic
!
!------------------------------------------------------------------------------
!
subroutine init_restoring_vel_Arctic
  ! init open boundary velocity for Arctic setup_Arctic
  ! In this version the barotropic velocity is specified (derived from
  ! streamfunction offline)
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none
  !
  integer                     :: i, el2, elnodes2(3)
  integer                     :: q, row, mn(3)
  real(kind=8)                :: vol, aux
  real(kind=8)                :: tri_u(3), tri_v(3)
  real(kind=8)                :: uu, vv
  real(kind=8), allocatable   :: opbnd_u2d(:), opbnd_v2d(:)
  integer, allocatable        :: map_loc(:)
  ! read barotropic velocity on open boundaries

  if(nmbr_opbnd_tri>0) then   !=======================    
     allocate(opbnd_ssh_rhs(nmbr_opbnd_t2d))
     allocate(map_loc(nod2d))                          
     map_loc=0                                         
     do i=1, nmbr_opbnd_n2D                                
        map_loc(myList_nod2D(opbnd_n2d(i)))=i        
     end do
     allocate(opbnd_u2d(nmbr_opbnd_t2D))
     allocate(opbnd_v2d(nmbr_opbnd_t2D))
     opbnd_u2d=0.
     opbnd_v2d=0.
     open(20,file=trim(OpbndPath)//'Arc_ob_int_vel.out', status='old')
     do i=1,nod2d
        read(20,*) uu, vv         !depth integrated velocity
        if(map_loc(i).ne.0) then                        
           opbnd_u2d(map_loc(i))=uu                    
           opbnd_v2d(map_loc(i))=vv                    
        end if
     end do
     close(20) 
     ! depth averaged velocity
     opbnd_u2d=opbnd_u2d/opbnd_dep
     opbnd_v2d=opbnd_v2d/opbnd_dep

     ! prepare ssh_rhs contribution from open boundary
     ! and correct the result to make sure IN = OUT
     opbnd_ssh_rhs=0.0
     do el2=1, nmbr_opbnd_tri
        elnodes2=opbnd_tri(el2,1:3)
        elnodes2=nod2d_corresp_to_nod3d(elnodes2)
        mn=mapping_opbnd_n2d(elnodes2)
        vol=opbnd_nv(el2,4)/12.0
        tri_u=opbnd_u2d(mn)
        tri_v=opbnd_v2d(mn)
        tri_v=tri_u*opbnd_nv(el2,1)+tri_v*opbnd_nv(el2,2)
        aux=sum(tri_v)
        do q=1, 3
           row=elnodes2(q)
           opbnd_ssh_rhs(mn)=opbnd_ssh_rhs(mn) - (aux+tri_v(q))*vol 
        end do
     end do
     deallocate(opbnd_v2d, opbnd_u2d, map_loc)
  end if                        !=======================    

  ! correction
  aux=0.0
  row=0
  if(nmbr_opbnd_tri>0) then                         
     ! Summation should go over my nodes:            
     do i=1, nmbr_opbnd_t2D     ! all 2D OB nodes    
        q=opbnd_n2D(i)                               
        if(q<=myDim_nod2D) then                      
           aux=aux+opbnd_ssh_rhs(i)                   
           row=row+1                                 !! .... 
        end if
     end do
  end if

  vol=0.0
  q=0
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( aux, vol, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE( row, q, 1, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  ! vol contains global sum of ssh_rhs, and q the global number 
  ! of open boundary nodes

  if(nmbr_opbnd_tri>0) then                         
     opbnd_ssh_rhs=opbnd_ssh_rhs-vol/real(q)
  end if

  if(mype==0) write(*,*) 'restoring open boundary velocity ready'

end subroutine init_restoring_vel_Arctic
subroutine init_cavity_ts
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_PARFE
  implicit none
  
  integer                     :: i, j, n, n2, n_bd, row, ind
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  integer, allocatable        :: ind_op(:), ind_op_glo(:), nearest_bd_nod(:)
  real(kind=8)                :: pp, pr, aux1,aux2, lon, lat, d, d_min
  real(kind=8)                :: x, y, temp_x, temp_y, temp_z
  real(kind=8), external      :: theta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data_T(:,:,:), raw_data_S(:,:,:)
  real(kind=8), allocatable   :: x_op(:), y_op(:), dep_op(:)
  real(kind=8), allocatable   :: x_op_glo(:), y_op_glo(:), dep_op_glo(:)

  ! define the cavity boundary line
  allocate(ind_op(nod2d), x_op(nod2d), y_op(nod2d), dep_op(nod2d))
  ind_op=0
  x_op=0.0
  y_op=0.0
  dep_op=0.0
  allocate(ind_op_glo(nod2d), x_op_glo(nod2d), y_op_glo(nod2d), dep_op_glo(nod2d))
  ind_op_glo=0
  x_op_glo=0.0
  y_op_glo=0.0
  dep_op_glo=0.0

  do row=1,myDim_nod2d  ! should not include eDim_nod2d 
     if(cavity_flag_nod2d(row)==0 .and. &
          any(cavity_flag_nod2d(nghbr_nod2d(row)%addresses)==1)) then 
        n=myList_nod2d(row)
        ind_op(n)=1
        x_op(n)=coord_nod2d(1,row)
        y_op(n)=coord_nod2d(2,row)
        if(rotated_grid) then
           call r2g(lon, lat, x_op(n), y_op(n))
           x_op(n)=lon
           y_op(n)=lat
        end if
        dep_op(n)=coord_nod3d(3,bt_nds(row))
     end if
  end do

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE(ind_op, ind_op_glo, &
       nod2d, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(x_op, x_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(y_op, y_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  call MPI_AllREDUCE(dep_op, dep_op_glo, &
       nod2d, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  deallocate(ind_op, x_op, y_op, dep_op)

  allocate(nearest_bd_nod(myDim_nod2d+eDim_nod2d))
  do row=1,myDim_nod2d+eDim_nod2d
     if(cavity_flag_nod2d(row)==0) cycle
     
     x=coord_nod2d(1,row)
     y=coord_nod2d(2,row)
     if(rotated_grid) then
        call r2g(lon, lat, x, y)
        x=lon
        y=lat
     end if
          
     d_min=3000.0e3   !dist. in m
     do n=1,nod2d
        if(ind_op_glo(n)/=1) cycle
        call dist_on_earth(x, y, x_op_glo(n), y_op_glo(n), d)
        if(d<d_min) then
           ind=n
           d_min=d
        end if
     end do
     nearest_bd_nod(row)=ind
  end do

  ! open global T/S dataset files
  open(19,file=trim(ClimateDataPath)//trim(OceClimaDataName), status='old')
  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data_T(num_lon_reg,num_lat_reg,num_lay_reg))
  allocate(raw_data_S(num_lon_reg,num_lat_reg,num_lay_reg))
  ! read raw data: T
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data_T(i,j,1:num_lay_reg)         
     end do
  end do
  ! read raw data: S
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data_S(i,j,1:num_lay_reg)         
     end do
  end do
  close(19) 

  ! interpolate to the nodes under cavity
  do row=1, myDim_nod3d+eDim_nod3D  
     n2=nod2d_corresp_to_nod3d(row)
     if(cavity_flag_nod2d(n2)==0) cycle
     n_bd=nearest_bd_nod(n2)
     temp_z=coord_nod3d(3,row)
     if(temp_z<dep_op_glo(n_bd)) temp_z=dep_op_glo(n_bd)
     temp_x=x_op_glo(n_bd)/rad
     temp_y=y_op_glo(n_bd)/rad
     ! change lon range to [0 360]
     if(temp_x<0.) temp_x=temp_x + 360.0  
     
     call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
          lon_reg, lat_reg, lay_reg, &
          raw_data_T, 1, temp_x, temp_y, temp_z, aux1)

     call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
          lon_reg, lat_reg, lay_reg, &
          raw_data_S, 1, temp_x, temp_y, temp_z, aux2)
     tracer(row,2)=aux2

     ! Convert in situ temperature into potential temperature
     pr=0.0_8
     pp=abs(temp_z)
     tracer(row,1)=theta(aux2, aux1, pp, pr)
  end do    

  deallocate(ind_op_glo, x_op_glo, y_op_glo, dep_op_glo, nearest_bd_nod)
  deallocate(raw_data_T, raw_data_S, lay_reg, lat_reg, lon_reg)
end subroutine init_cavity_ts
!
!----------------------------------------------------------------------------------------
!
subroutine cavity_heat_water_fluxes
  ! compute the heat and freshwater fluxes under ice cavity
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer        :: m, n, row
  real(kind=8)   :: gama, L
  real(kind=8)   :: c2, c3, c4, c5, c6
  real(kind=8)   :: t_i, s_i, p, t_fz

  ! parameter for computing heat and water fluxes
  gama = 1.0e-4     ! heat exchange velocity [m/s]
  L    = 334000.    ! water to ice latent heat [J/Kg], same as used by the ice model

  ! parameter for computing freezing temperature
  c3 = 1.710523e-3
  c4 = -2.154996e-4
  c5 = -0.0575
  c6 = -7.53e-4

  do n=1,myDim_nod2D+eDim_nod2D      
     if(cavity_flag_nod2d(n)==0) cycle   
     row=nod3d_below_nod2d(1,n)
     t_i = tracer(row,1)	
     s_i = tracer(row,2)
     t_fz = c3*(s_i**(3./2.)) + c4*(s_i**2) + c5*s_i + c6*abs(coord_nod3d(3,row))
     !
     heat_flux(n)=vcpw*gama*(t_i - t_fz)  !by Hunter2006 cpw=3974J/Kg is used
     water_flux(n) = -1.0*heat_flux(n)/(L*1000.0)  
  enddo

end subroutine cavity_heat_water_fluxes
!
!----------------------------------------------------------------------------------------
!
!==========================================================================
!
subroutine init_tidal_opbnd
  ! initialize tides: read constituents (amplitude and phase)
  ! The global tidal data should be on a grid with longitude range 0-360 degree.
  ! Assumption: the number of columes/rows of the grids for global tidal data
  !             are the same for all constituents and types (U,V, z), though 
  !             the exact lon/lat values can be (are) different for different
  ! 	        types.
  use o_param
  use o_array
  use o_mesh
  use g_config
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer					:: i, j, n, num, fid, nodmin,m(1)
  integer					:: num_lon_reg, num_lat_reg
  integer                                       :: p_flag
  character(2)    				:: cons_name
  real(kind=8)					:: lon, lat, dep
  real(kind=8)					:: x, y, x2, y2, d, dmin
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4u, lat_reg_4u
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4v, lat_reg_4v
  real(kind=8), allocatable, dimension(:)	:: lon_reg_4z, lat_reg_4z
  real(kind=8), allocatable, dimension(:,:)	:: amp_reg, pha_reg
  real(kind=8), allocatable, dimension(:)	:: lon_opbnd_n2d, lat_opbnd_n2d

  ! check
  num=len(trim(tidal_constituent))
  if(mod(num,2) /= 0 .or. num/2 /= nmbr_tidal_cons) then
     write(*,*) 'wrong specification of tidal constituents in O_module'
     call par_ex
     stop
  end if

  ! allocate arrays
  if(trim(tide_opbnd_type)=='Flather') then
     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
     allocate(opbnd_z_tide(nmbr_opbnd_t2D))
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     opbnd_z_tide=0.
  elseif(trim(tide_opbnd_type)=='ssh') then
     allocate(opbnd_z_tide(nmbr_opbnd_t2D))
     allocate(opbnd_z0_tide(nmbr_opbnd_t2D))
     opbnd_z_tide=0.  
     opbnd_z0_tide=0.
  else !'vel'
     allocate(opbnd_u_tide(nmbr_opbnd_t2D))
     allocate(opbnd_v_tide(nmbr_opbnd_t2D))
     opbnd_u_tide=0.
     opbnd_v_tide=0.
  end if

  ! allocate tidal amp and pha arrays
  allocate(tide_period_coeff(nmbr_tidal_cons))
  tide_period_coeff=0.
  if(trim(tide_opbnd_type)/='ssh') then
     allocate(tide_u_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_v_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_u_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_v_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     tide_u_amp=0.
     tide_v_amp=0.
     tide_u_pha=0.
     tide_v_pha=0.
  end if
  if(trim(tide_opbnd_type)/='vel') then
     allocate(tide_z_amp(nmbr_opbnd_t2d, nmbr_tidal_cons))
     allocate(tide_z_pha(nmbr_opbnd_t2d, nmbr_tidal_cons))
     tide_z_amp=0.
     tide_z_pha=0.
  end if

  ! read the regular grid (for velocity)
  if(trim(tide_opbnd_type)/='ssh') then
     open(101,file=trim(TideForcingPath)//'lonlat_4U_'//trim(tidemodelname)//'.dat', status='old')
     read(101,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4u(num_lon_reg), lat_reg_4u(num_lat_reg))
     do i=1,num_lon_reg
        read(101,*) lon_reg_4u(i)
     end do
     do i=1,num_lat_reg
        read(101,*) lat_reg_4u(i)
     end do
     close(101)
     open(102,file=trim(TideForcingPath)//'lonlat_4V_'//trim(tidemodelname)//'.dat', status='old')
     read(102,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4v(num_lon_reg), lat_reg_4v(num_lat_reg))
     do i=1,num_lon_reg
        read(102,*) lon_reg_4v(i)
     end do
     do i=1,num_lat_reg
        read(102,*) lat_reg_4v(i)
     end do
     close(102)
  end if

  ! read the regular grid (for ssh)
  if(trim(tide_opbnd_type)/='vel') then
     open(103,file=trim(TideForcingPath)//'lonlat_4z_'//trim(tidemodelname)//'.dat', status='old')
     read(103,*) num_lon_reg, num_lat_reg 
     allocate(lon_reg_4z(num_lon_reg), lat_reg_4z(num_lat_reg))
     do i=1,num_lon_reg
        read(103,*) lon_reg_4z(i)
     end do
     do i=1,num_lat_reg
        read(103,*) lat_reg_4z(i)
     end do
     close(103)
  end if

  ! allocate arrays for global data on regular grids
  allocate(amp_reg(num_lon_reg, num_lat_reg), pha_reg(num_lon_reg, num_lat_reg)) 

  ! open-boundary nodes coordinates
  allocate(lon_opbnd_n2d(nmbr_opbnd_n2d), lat_opbnd_n2d(nmbr_opbnd_n2d))
  do i=1, nmbr_opbnd_n2d
     n=opbnd_n2d(i)
     if(rotated_grid) then
        call r2g(lon, lat, coord_nod2d(1,n), coord_nod2d(2,n))
        lon_opbnd_n2d(i)=lon/rad   ! change unit to degree
        lat_opbnd_n2d(i)=lat/rad
     else
        lon_opbnd_n2d(i)=coord_nod2d(1,n)/rad
        lat_opbnd_n2d(i)=coord_nod2d(2,n)/rad
     end if
     ! change lon range to [0 360]
     if(lon_opbnd_n2d(i)<0.) lon_opbnd_n2d(i)=lon_opbnd_n2d(i) + 360.0  
  end do

  ! read and interpolate each constituent (for U,V,z and their phase) 
  do n=1,nmbr_tidal_cons
     cons_name=tidal_constituent(2*n-1:2*n)
     call tide_period(cons_name, tide_period_coeff(n))

     if(trim(tide_opbnd_type)/='ssh') then
        fid=103+n
        open(fid,file=trim(TideForcingPath)//'U_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_amp(1:nmbr_opbnd_n2d,n),p_flag)
        tide_u_amp(:,n)=tide_u_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4u, lat_reg_4u, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_u_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_u_pha(:,n)=tide_u_pha(:,n)*rad   ! change units of phase from degree to radian 


        open(fid,file=trim(TideForcingPath)//'V_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_amp(1:nmbr_opbnd_n2d,n),p_flag)
        tide_v_amp(:,n)=tide_v_amp(:,n)/opbnd_dep  ! change transport (m^2/s) to depth mean velocity (m/s)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4v, lat_reg_4v, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_v_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_v_pha(:,n)=tide_v_pha(:,n)*rad   ! change units of phase from degree to radian 

     end if


     if(trim(tide_opbnd_type)/='vel') then
        open(fid,file=trim(TideForcingPath)//'z_'//cons_name//'_'//trim(tidemodelname)//'.dat', status='old')
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) amp_reg(i,j)         
           end do
        end do
        do i=1, num_lon_reg
           do j=1, num_lat_reg
              read(fid, *) pha_reg(i,j)         
           end do
        end do
        close(fid)

        p_flag=0
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, amp_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_amp(1:nmbr_opbnd_n2d,n),p_flag)

        p_flag=1
        call interp_2d_field(num_lon_reg, num_lat_reg, lon_reg_4z, lat_reg_4z, pha_reg, &
             nmbr_opbnd_n2D, lon_opbnd_n2d, lat_opbnd_n2d, tide_z_pha(1:nmbr_opbnd_n2d,n),p_flag)
        tide_z_pha(:,n)=tide_z_pha(:,n)*rad   ! change units of phase from degree to radian 
     end if
  end do

  ! amplifying the magnitude for process studies
  if(trim(tide_opbnd_type)/='ssh') then
     tide_u_amp=tide_u_amp*tide_amplify_coeff
     tide_v_amp=tide_v_amp*tide_amplify_coeff
  end if
  if(trim(tide_opbnd_type)/='vel') then
     tide_z_amp=tide_z_amp*tide_amplify_coeff
  end if


  ! deallocate temporary arrays
  deallocate(lat_opbnd_n2d, lon_opbnd_n2d)
  deallocate(pha_reg, amp_reg)
  if(trim(tide_opbnd_type)/='ssh') deallocate(lat_reg_4v,lon_reg_4v,lat_reg_4u,lon_reg_4u)
  if(trim(tide_opbnd_type)/='vel') deallocate(lat_reg_4z,lon_reg_4z)

  if(mype==0) write(*,*) 'Tidal constituents ', trim(tidal_constituent), ' have been loaded for opbnd'
end subroutine init_tidal_opbnd
!
!==========================================================================
!
subroutine update_tidal_opbnd
  use o_param
  use o_array
  use o_mesh
  use g_config
  implicit none
  integer		:: n
  real(kind=8)		:: aux
  !
  aux=2.0*pi*istep*dt
  if(trim(tide_opbnd_type)=='Flather') then
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     opbnd_z_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
     end do
  elseif(trim(tide_opbnd_type)=='ssh') then
     opbnd_z0_tide=opbnd_z_tide
     opbnd_z_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_z_tide=opbnd_z_tide + tide_z_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_z_pha(:,n))
     end do
  else
     opbnd_u_tide=0.
     opbnd_v_tide=0.
     do n=1, nmbr_tidal_cons
        opbnd_u_tide=opbnd_u_tide + tide_u_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_u_pha(:,n))
        opbnd_v_tide=opbnd_v_tide + tide_v_amp(:,n)*cos(aux/tide_period_coeff(n)-tide_v_pha(:,n))
     end do
  end if
  !
end subroutine update_tidal_opbnd
!
!==========================================================================
!
subroutine tide_period(cons_name, period)
  use g_parfe
  use g_config
  implicit none
  real(kind=8)		:: period
  real(kind=8)          :: sec_p_h
  character(2) 		:: cons_name
  !  
  sec_p_h=3600.0_8
  if(cons_name=='M2') then
     period=12.420601*sec_p_h
  elseif(cons_name=='S2') then
     period=12.0*sec_p_h
  elseif(cons_name=='N2') then
     period=12.658348*sec_p_h 
  elseif(cons_name=='K2') then
     period=11.967235*sec_p_h  
  elseif(cons_name=='K1') then
     period=23.93447*sec_p_h 
  elseif(cons_name=='O1') then
     period=25.819342*sec_p_h
  elseif(cons_name=='P1') then
     period=24.06589*sec_p_h
  elseif(cons_name=='Q1') then
     period=26.868357*sec_p_h
  elseif(cons_name=='S1') then
     period=24.0*sec_p_h
  elseif(cons_name=='Mf') then
     period=13.661*24.0*sec_p_h
  elseif(cons_name=='Mm') then
     period=27.555*24.0*sec_p_h
  else
     write(*,*) 'Error! One or more tidal constituents you selected are not specified in the code.'
     call par_ex
     stop
  endif
end subroutine tide_period
!
subroutine cal_shortwave_rad
  ! Calculate shortwave radiation in the ocean according to chlorophyll climatology
  ! Assuming under ice no penetration. A decent way for ice region is to be discussed.
  ! This routine should be called after ice2oce coupling done if ice model is used.
  ! Ref.: Morel and Antoine 1994, Sweeney et al. 2005
  use o_mesh
  use o_param
  use o_array
  use g_forcing_arrays
  use g_parfe
  use g_config
  use i_therm_parms
  use i_array

  implicit none

  integer      :: m, n2, n3, k
  real(kind=8) :: swsurf, aux, z
  real(kind=8) :: c, c2, c3, c4, c5
  real(kind=8) :: v1, v2, sc1, sc2

  sw_3d=0.0

  do n2=1,myDim_nod2D+eDim_nod2D     
#ifdef use_cavity
     if(cavity_flag_nod2d(n2)==1) cycle   
#endif
#ifdef use_ice  
     if(a_ice(n2)>0.0) cycle !assume in ice region no penetration
#endif     

     ! shortwave rad.
     swsurf=(1.0-albw)*shortwave(n2)
     ! the visible part (300nm-750nm)
     swsurf=swsurf*0.54
     ! subtract visible sw rad. from heat_flux, which is '+' for upward
     heat_flux(n2)=heat_flux(n2)+swsurf 

     ! attenuation func. for vis. sw rad. according to Morel/Antoine param.
     ! the four parameters in the func.

     ! limit chl from below
     if(chl(n2)<0.02) chl(n2)=0.02

     c=log10(chl(n2))  
     c2=c*c
     c3=c2*c
     c4=c3*c
     c5=c4*c
     v1=0.008*c+0.132*c2+0.038*c3-0.017*c4-0.007*c5
     v2=0.679-v1
     v1=v1+0.321
     sc1=1.54-0.197*c+0.166*c2-0.252*c3-0.055*c4+0.042*c5
     sc2=7.925-6.644*c+3.662*c2-1.815*c3-0.218*c4+0.502*c5

     ! convert from heat flux [W/m2] to temperature flux [K m/s]
     swsurf=swsurf/vcpw

     ! vis. sw. rad. in the colume
     sw_3d(nod3d_below_nod2d(1,n2))=swsurf
     do k=2,num_layers_below_nod2d(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        z=coord_nod3d(3,n3)
        aux=(v1*exp(z/sc1)+v2*exp(z/sc2))
        sw_3d(n3)=swsurf*aux
        if(aux<0.00001) exit
     end do
  end do

end subroutine cal_shortwave_rad
subroutine check_imb_freshwater
  use o_param
  use o_mesh 
  use o_elements
  use o_array
  use g_config
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer      :: row
  real(kind=8) :: flux, corr

  flux=0.0
  do row=1,myDim_nod2D
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle   
#endif 
     flux=flux+(evaporation(row)+prec_rain(row)+ &
          prec_snow(row)+runoff(row))*cluster_area_2d(row)
  end do

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  corr=0.0
  call MPI_AllREDUCE( flux, corr, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  corr=corr/ocean_area
  water_flux=water_flux+corr  ! the + sign should be used here

end subroutine check_imb_freshwater
!
!-------------------------------------------------------------------------
!
subroutine check_imb_salt_flux
  use o_param
  use o_mesh 
  use o_elements
  use o_array
  use g_config
  use g_parfe
  implicit none

  integer      :: row
  real(kind=8) :: flux_rel, flux_vir, corr

  flux_rel=0.0
  flux_vir=0.0

  do row=1,myDim_nod2D
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle   
#endif 
     flux_rel=flux_rel+relax_salt(row)*cluster_area_2d(row)
#ifndef use_fullfreesurf
     flux_vir=flux_vir+virtual_salt(row)*cluster_area_2d(row)
#endif
  end do

  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  corr=0.0
  call MPI_AllREDUCE(flux_rel, corr, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  corr=corr/ocean_area
  relax_salt=relax_salt-corr

#ifndef use_fullfreesurf  
  corr=0.0
  call MPI_AllREDUCE(flux_vir, corr, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  corr=corr/ocean_area
  virtual_salt=virtual_salt-corr
#endif
  
end subroutine

subroutine landice_water_init
  ! init land ice melting rate
  use o_mesh
  use o_elements
  use o_param
  use o_array
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  integer                     :: n, i, j, num_reg, num_nod, num_days
  integer                     :: n_loc, fileID
  integer, allocatable        :: temp_arr2d(:), nodes_in_region(:)
  real(kind=8)                :: vol, vol_sum, aux
  real(kind=8), allocatable   :: totalflux(:)
  character*80                :: file_name
  character                   :: c_reg_ind

  ! read the number of region and the total yearly discharge in each region:
  file_name=trim(meshpath)//'landice_yearly_mass_loss.out' 
  fileID=160
  open(fileID, file=file_name)
  read(fileID,*) num_reg
  allocate(totalflux(num_reg))
  read(fileID,*) totalflux      !unit: Gt/year
  close(fileID)

  allocate(temp_arr2d(nod2d))
  temp_arr2d=0
  do n=1, myDim_nod2D  !note: eDim_nod2D should not be included in this case
     temp_arr2d(myList_nod2D(n))=n
  end do

  ! yearly mean runoff
  do i=1,num_reg
     
     !read in nodes in the region
     write(c_reg_ind,'(i1)') i   !assume that num_reg is less than 10
     file_name=trim(meshpath)//'landice_nodes_in_region_'//c_reg_ind//'.out' 
     fileID=160
     open(fileID, file=file_name)
     read(fileID,*) num_nod
     allocate(nodes_in_region(num_nod))
     read(fileID,*) nodes_in_region
     close(fileID)

     vol=0.0
     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           do j=1,nod_in_elem2d(n_loc)%nmb
              vol=vol+voltriangle(nod_in_elem2d(n_loc)%addresses(j))
           end do
        end if
     end do

     vol_sum=0.0
     call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
     call MPI_AllREDUCE(vol, vol_sum, 1, MPI_DOUBLE_PRECISION,MPI_SUM, &
          MPI_COMM_WORLD, MPIerr)
 
     vol_sum=vol_sum/3.0
     
     aux=totalflux(i)*1.0e9/real(ndpyr)/86400.0_8/vol_sum  !m/s

     do n=1,num_nod
        n_loc=temp_arr2d(nodes_in_region(n))
        if(n_loc>0) then
           runoff_landice(n_loc)=aux
        end if
     end do

     deallocate(nodes_in_region)
  end do

  call com_2D(runoff_landice)

  ! seasonality
  num_days=sum(num_day_in_month(fleapyear,landice_start_mon:landice_end_mon))
  aux=real(ndpyr)/real(num_days)
  landice_season=0.0
  landice_season(landice_start_mon:landice_end_mon)=aux
    
  deallocate(temp_arr2d, totalflux)

  if(mype==0) write(*,*) 'Land-ice melt water fluxes prepared.'

end subroutine landice_water_init
!
!----------------------------------------------------------------------------------
!
subroutine add_landice_water
  use o_array
  use g_forcing_arrays
  use g_clock
  use g_parfe
  implicit none

  water_flux=water_flux-runoff_landice*landice_season(month)
  !water_flux is positive for upward

end subroutine add_landice_water
subroutine solve(ident)

  ! All stiffness matrices are in the transposed form ?
  ! PIL_PMVALS --- each processor takes only owned rows
  ! PIL_MVALS  --- matrix should be present at each processor
  ! RHS  --- only owned rows are taken as default

  use g_PARFE
  use g_config
  use o_PARAM
  use o_array
  use o_MATRICES
  use o_MESH
  use o_solver
  implicit none

#include "petscf.h"
  !include "petscf.h"
  !include "pilutf.h"
  !include "hypref.h"

  integer,intent(IN)    :: ident
  integer               :: Pmode, n3 
  integer               :: maxiter, restart, lutype, fillin
  integer               :: COMM, myrows
  real(kind=8)          :: droptol, soltol, rinfo(20,20)
  save rinfo

  COMM=MPI_COMM_WORLD

  maxiter=2000
  restart=15
  fillin=3
  lutype=2
  droptol=1.e-6
  soltol=1.e-6
  n3=myDim_nod3D+eDim_nod3D
  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

#ifdef use_fullfreesurf     
  ! matrix is updated

  if (ident==solve_u) then  ! du* or du
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM &
          + PET_QUIET 
     if (.not. iteruv_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
     call PETSC_S(Pmode, &
          ident, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1:n3), duf(1:n3), rinfo(:,1), COMM) 
  endif

  if (ident==solve_v) then  ! dv* or dv
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3:n3*2), duf(1+n3:n3*2), rinfo(:,2), COMM) 
  endif

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3*2:n3*3), duf(1+n3*2:n3*3), rinfo(:,3), COMM) 
  endif
#endif

  if (ident==solve_ssh) then  !dssh
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_ILU + PET_PCBJ + PET_QUIET
     !if (.not. iter_first)   Pmode = Pmode - PET_STRUCT
     if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
     call PETSC_S(Pmode,ident,sshstiff%dim, sshstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*1.0e-2, &  
          soltol*1.0e-2,  &   ! Sol. Residuum reduction 
          part2D, sshstiff%rowptr, sshstiff%colind, sshstiff%values, &
          ssh_rhs, dssh, rinfo(:,4), COMM) 
  endif

#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
     Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_ILU + PET_PCBJ +PET_QUIET 
     !if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  + PET_NEWPC
     if (.not. iter_first)   Pmode = Pmode - PET_STRUCT  
     call PETSC_S(Pmode,ident,nhpstiff%dim, nhpstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*10., &  
          soltol*10.,  &  
          part3D, nhpstiff%rowptr, nhpstiff%colind, nhpstiff%values, &
          nhp_rhs, nhp, rinfo(:,5), COMM) 
  endif
#endif

  if(ident>=solve_tra) then
     if(ident==solve_tra) then  ! dTF
        Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM &
             + PET_QUIET
        if (.not. iter_first)  Pmode = Pmode - PET_STRUCT  + PET_NEWPC 
        call PETSC_S(Pmode, &
             ident,tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  ! droptol 
             soltol, &   ! Sol. Residuum reduction 
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,1), dtracer(:,1), rinfo(:,6), COMM) 
     else  !other tracers
        Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET  
        call PETSC_S(Pmode , &
             solve_tra, &
             tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &   ! Sol. Residuum reduction 
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,ident-solve_tra+1), dtracer(:,ident-solve_tra+1), rinfo(:,7), COMM)
     endif
  endif


  !-------------------------------------------------------------------------------------
#else   
  ! linear surface

  if (ident==solve_u) then  ! du* or du
     Pmode =PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
     if (iteruv_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_RCM
     call PETSC_S(Pmode, &
          ident, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1:n3), duf(1:n3), rinfo(:,1), COMM) 
  endif

  if (ident==solve_v) then  ! dv* or dv
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3:n3*2), duf(1+n3:n3*2), rinfo(:,2), COMM) 
  endif

#ifdef use_non_hydrostatic
  if (ident==solve_w) then  ! dw* or dw for non-hydrostatic case
     Pmode = PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET 
     call PETSC_S(Pmode, &
          solve_u, uvstiff%dim, uvstiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, & 
          soltol,  &   
          part3D, uvstiff%rowptr, uvstiff%colind, uvstiff%values, &
          uv_rhs(1+n3*2:n3*3), duf(1+n3*2:n3*3), rinfo(:,3), COMM) 
  endif
#endif

  if (ident==solve_ssh) then  ! dssh
     Pmode =PET_BLOCKP+PET_SOLVE+ PET_BICGSTAB+PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_PCBJ + PET_ILU
     call PETSC_S(Pmode,ident,sshstiff%dim, sshstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*1.0e-2, &  
          soltol*1.0e-2,  &   ! Sol. Residuum reduction 
          part2D, sshstiff%rowptr, sshstiff%colind, sshstiff%values, &
          ssh_rhs, dssh, rinfo(:,4), COMM) 
  endif

#ifdef use_non_hydrostatic  
  if (ident==solve_nhp) then
     Pmode =PET_BLOCKP+PET_SOLVE+ PET_BICGSTAB !+PET_QUIET 
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS+PET_PCBJ+PET_ILU
     call PETSC_S(Pmode,ident,nhpstiff%dim, nhpstiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol*10., &  
          soltol*10.,  &  
          part3D, nhpstiff%rowptr, nhpstiff%colind, nhpstiff%values, &
          nhp_rhs, nhp, rinfo(:,5), COMM) 
  endif
#endif

  if(ident>=solve_tra) then
     if(ident==solve_tra) then  ! dTF
#ifdef use_tracer_gls
        Pmode = PET_BLOCKP+PET_STRUCT + PET_PMVALS + PET_SOLVE + PET_BICGSTAB + PET_RCM + PET_QUIET
        if (.not. iter_first)  Pmode = Pmode - PET_STRUCT  + PET_NEWPC
#else
        Pmode =PET_BLOCKP+PET_SOLVE + PET_BICGSTAB + PET_QUIET
        if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS + PET_RCM
#endif
        call PETSC_S(Pmode, &
             ident,tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &   
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,1), dtracer(:,1), rinfo(:,2), COMM) 

     else  !other tracers

        Pmode =PET_BLOCKP+ PET_SOLVE + PET_BICGSTAB + PET_QUIET  
        call PETSC_S(Pmode , &
             solve_tra, &
             tsstiff%dim, tsstiff%nza,  myrows, &
             maxiter, & 
             restart, &  
             fillin,  &  
             droptol, &  
             soltol, &  
             part3D, tsstiff%rowptr, tsstiff%colind, tsstiff%values, &
             tracer_rhs(:,ident-solve_tra+1), dtracer(:,ident-solve_tra+1), rinfo(:,2), COMM)
     endif
  endif

#endif
end subroutine solve


!=======================================================================
subroutine ice2ocean
  !transmits the relevant fields from the ice to the ocean model
  !Ralph Timmermann, 15.9.2004

  use o_PARAM
  use o_array
  use o_mesh
  use i_array
  use g_parfe
  use i_dyn_parms
  use i_therm_parms
  use g_config
  implicit none

  integer          :: row
  real*8           :: aux

  do row=1,myDim_nod2d+eDim_nod2d            
#ifdef use_cavity
     if(cavity_flag_nod2d(row)==1) cycle
#endif   

     !heat flux:
     heat_flux(row) = -net_heat_flux(row)     ! W/m2  !heat_flux >0 f. upward 
     ! net_heat_flux > 0 f. downward
     heat_flux(row)=0.0
     !fresh water flux:
     water_flux(row) = -fresh_wa_flux(row)    ! m freshwater/s  !water_flux > 0 f. upward
     ! fresh_wa_flux > 0 f. downward
     water_flux(row)=0.0
     !momentum flux:
     aux=sqrt((u_ice(row)-u_w(row))**2+(v_ice(row)-v_w(row))**2)*rho0*Cd_oce_ice
     stress_iceoce_x(row) = aux * (u_ice(row)-u_w(row))
     stress_iceoce_y(row) = aux * (v_ice(row)-v_w(row))
     stress_x(row)=stress_iceoce_x(row)*a_ice(row) + stress_atmoce_x(row)*(1.-a_ice(row))
     stress_y(row)=stress_iceoce_y(row)*a_ice(row) + stress_atmoce_y(row)*(1.-a_ice(row))
     if(a_ice(row)<0.001) then
        stress_iceoce_x(row)=0.0
        stress_iceoce_y(row)=0.0
     end if
  end do

end subroutine ice2ocean
!=======================================================================


!=======================================================================
subroutine ocean2ice
  !transmits the relevant fields from the ocean to the ice model
  !Ralph Timmermann, 15.9.2004

  use o_param
  use o_array
  use i_array
  use o_MESH
  use g_parfe
  use g_config
  implicit none

  integer :: m, row

  ! the arrays in the ice model are renamed

  do row=1, myDim_nod2d+eDim_nod2d       
     m=nod3D_below_nod2D(1,row)       
     T_oc_array(row)=tracer(m,1)  
     S_oc_array(row)=tracer(m,2)  
     u_w(row) = uf(m)                        
     v_w(row) = uf(m+myDim_nod3d+eDim_nod3D)  
     elevation(row)= ssh(row)
  enddo
end subroutine ocean2ice
!=======================================================================

!ice module iteration
!==============================================================================

subroutine ice_step
  use o_mesh
  use o_param
  use i_therm_parms
  use i_dyn_parms
  use i_array
  use i_solver
  use g_parfe
  use g_config
  implicit none
  !
  integer        :: m, row
  real(kind=8)   :: t0, t1, t2, t3

  t0=MPI_Wtime()

  !1) Dynamic part

  ! to avoid zero ice velocities (to avoid the solver's sudden death)
  do row=1,myDim_nod2d+eDim_nod2d    
     if (m_ice(row) < 0.001) then
        u_ice(row)=small*myList_nod2D(row)   !! To reach correspondence with 
        v_ice(row)=small*myList_nod2D(row)   !! global memory code
     endif
  enddo

  ! to solve u,v
  call rheology      

#ifdef use_cavity
  call clean_cavity_vel
#endif    
  t1=MPI_Wtime()

  ! to solve m, a, ms due to advection
#ifdef use_ice_gls
  call icestiff_fill_gls
  call solveIce(solve_m_ice)       
  call solveIce(solve_a_ice)       
  call solveIce(solve_m_snow)      
  m_ice=m_ice+dm_ice
  a_ice=a_ice+da_ice
  m_snow=m_snow+dm_snow 
  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)  
#else

#ifdef use_ice_fct
  call ice_rhs_tg

  call fct_ice_solve

#else
  call ice_rhs_tg                 
  if(lump_ice_matrix) then
     call ice_solve
  else                            
     call solveIce(solve_m_ice)       
     call solveIce(solve_a_ice)       
     call solveIce(solve_m_snow)      
  end if
  m_ice=m_ice+dm_ice
  a_ice=a_ice+da_ice
  m_snow=m_snow+dm_snow
  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)  

#endif
#endif

  call cut_off

#ifdef use_cavity
  call clean_cavity_m
#endif
  t2=MPI_Wtime()   

  !2) Thermodynamic part (update m, a, ms)

  call thermodynamics
  t3=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 
     write(*,*) 'ice took      ', t3-t0
     write(*,*) 'ice uv        ', t1-t0
     write(*,*) 'ice mass      ', t2-t1
     write(*,*) 'ice thermo    ', t3-t2
  endif
end subroutine ice_step
!
!------------------------------------------------------------------------------
!
subroutine clean_cavity_vel
  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: m, row

  do row=1,myDim_nod2d+eDim_nod2d           
     if(cavity_flag_nod2d(row)==1) then  
        u_ice(row)=0.
        v_ice(row)=0.
     endif
  enddo
end subroutine clean_cavity_vel
!
!------------------------------------------------------------------------------
!
subroutine clean_cavity_m
  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: m, row

  do row=1,myDim_nod2d+eDim_nod2d            
     if(cavity_flag_nod2d(row)==1) then 
        m_ice(row)=0.
        a_ice(row)=0.
        m_snow(row)=0.
     endif
  enddo
end subroutine clean_cavity_m
!
!------------------------------------------------------------------------------
!
subroutine ice_init
  !sets inital values or reads restart file for ice model
  use i_ARRAY
  use o_MESH    
  use o_PARAM   
  use o_array        
  use g_parfe 
  use g_config
  implicit none
  !
  integer        :: i
  character*100  :: filename

  if (.not.r_restart) then
     m_ice =0.
     a_ice =0.
     u_ice =0.
     v_ice =0.
     m_snow=0.

     if(use_prepared_init_ice) then
        if(mype==0) write(*,*) 'initialize the sea ice with previous simulation results'
        call read_prepared_initial_ice
     else           
        if(mype==0) write(*,*) 'initialize the sea ice according to initial SST'
 	do i=1,ToDim_nod2d          
           if(cavity_flag_nod2d(i)==1) cycle              
           if (tracer(nod3D_below_nod2D(1,i),1) < -1.0) then  
              m_ice(i) = 1.0
              a_ice(i) = 0.9
              u_ice(i) = 0.
              v_ice(i) = 0.
              m_snow(i)= 0. 
           endif
        enddo
     end if

  else

     if(mype==0) write(*,*) 'read ice restart file'
     call ice_input     

  endif
end subroutine ice_init
!
!----------------------------------------------------------------------------
!
subroutine ice_array_setup
  !inializing sea ice model 
  use o_param
  use o_mesh
  use o_elements
  use i_array
  use g_forcing_arrays
  use g_rotate_grid
  use g_parfe
  implicit none
  !
  integer       :: i, n2
  real(kind=8)  :: midlat, lon, lat, rlon, rlat

  n2=ToDim_nod2D           

  ! Allocate memory for variables of ice model      
  allocate(u_ice(n2), v_ice(n2))
  allocate(m_ice(n2), a_ice(n2), m_snow(n2))
  allocate(dm_ice(n2), da_ice(n2), dm_snow(n2))
  allocate(sigma11(myDim_elem2D), sigma12(myDim_elem2D), sigma22(myDim_elem2D))
  allocate(rhs_m(n2), rhs_a(n2), rhs_u(n2), rhs_v(n2))
  allocate(rhs_ms(n2))
  allocate(t_skin(n2))
  rhs_m=0.
  rhs_ms=0.
  rhs_a=0.
  rhs_u=0.
  rhs_v=0.
  u_ice=0.
  v_ice=0.
  m_ice=0.
  a_ice=0.
  m_snow=0.
  dm_ice=0.
  da_ice=0.
  dm_snow=0.
  sigma11=0.
  sigma12=0.
  sigma22=0.
  t_skin=0.

  ! Allocate memory used for coupling
  allocate(S_oc_array(n2), T_oc_array(n2))
  allocate(fresh_wa_flux(n2), net_heat_flux(n2))
  allocate(stress_atmice_x(n2), stress_atmice_y(n2))    
  allocate(stress_atmoce_x(n2), stress_atmoce_y(n2))    
  allocate(stress_iceoce_x(n2), stress_iceoce_y(n2)) 
  allocate(u_w(n2), v_w(n2))
  allocate(elevation(n2))

#ifdef use_ice_fct
  call fct_ice_init
#endif

  if(mype==0) write(*,*) 'ice arrays have been set up'   

end subroutine ice_array_setup



!===================================================================
subroutine thermodynamics
  !
  ! For every surface node, this subroutine extracts the information
  ! needed for computation of thermodydnamics, calls the relevant
  ! subroutine, and returns the information to the vectors of prognostic
  ! variables.
  !
  ! Originally programmed by N. Yakovlev/S. Danilov.
  !
  ! Adjusted for upgraded model physics (NCEP forcing data; parameterization
  ! of ice-ocean heat flux considering friction velocity) by Ralph Timmermann.
  !
  ! Adjusted for general forcing data and NlFs, 13.01.2009
  !===================================================================
  use o_param
  use o_mesh
  use o_elements 
  use o_array,  only: real_salt_flux
  use i_therm_parms
  use i_dyn_parms
  use i_array
  use g_config
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe

  implicit none
  real*8  :: h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss,rsf
  real*8  :: ug,ustar,T_oc,S_oc,h_ml,t,ch,ce,ch_i,ce_i,fw,ehf,evap
  real*8  :: ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout
  integer :: i, j

  rsss=ref_sss

  do i=1, myDim_nod2d+eDim_nod2d   
#ifdef use_cavity
     if(cavity_flag_nod2d(i)==1) cycle
#endif
     h       = m_ice(i)
     hsn     = m_snow(i)
     A       = a_ice(i)
     fsh     = shortwave(i)
     flo     = longwave(i)
     Ta      = Tair(i)
     qa      = shum(i)  
     if(precip_data_source=='NCEP') then
        if(Ta>=0.0) then
           rain=prec_rain(i)
           snow=0.0
        else
           rain=0.0
           snow=prec_rain(i)
        endif
     else
        rain = prec_rain(i)
        snow = prec_snow(i)
     end if
     runo    = runoff(i)
     ug      = sqrt(u_wind(i)**2+v_wind(i)**2)
     ustar   = sqrt(Cd_oce_ice)*sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
     T_oc    = T_oc_array(i)      
     S_oc    = S_oc_array(i)
     if(ref_sss_local) rsss = S_oc
     t       = t_skin(i)   
     ch	     = Ch_atm_oce_arr(i)
     ce	     = Ce_atm_oce_arr(i)
     ch_i    = Ch_atm_ice
     ce_i    = Ce_atm_ice
     h_ml    = 10.       	         ! 10.0 or 30. used previously
     fw      = 0.
     ehf     = 0.

     call therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
          ug,ustar,T_oc,S_oc,h_ml,t,dt,ch,ce,ch_i,ce_i,fw,ehf,evap, &
          rsf, ithdgr, ithdgrsn, iflice, hflatow, hfsenow, hflwrdout)

     m_ice(i)         = h
     m_snow(i)        = hsn
     a_ice(i)         = A
     t_skin(i)        = t
     fresh_wa_flux(i) = fw
     net_heat_flux(i) = ehf
     evaporation(i)   = evap    !negative up
#ifdef use_fullfreesurf
     real_salt_flux(i)= rsf
#endif         

     thdgr(i)         = ithdgr
     thdgrsn(i)       = ithdgrsn
     flice(i)         = iflice
     olat_heat(i)     = hflatow
     osen_heat(i)     = hfsenow
     olwout(i)        = hflwrdout

  end do

end subroutine thermodynamics
!
!===================================================================
!
subroutine therm_ice(h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss, &
     ug,ustar,T_oc,S_oc,H_ML,t,dt,ch,ce,ch_i,ce_i,fw,ehf,evap, &
     rsf, dhgrowth, dhsngrowth, iflice, hflatow, hfsenow, hflwrdout)
  ! Ice Thermodynamic growth model     
  !
  ! Input parameters:
  !------------------
  ! h - ice mass [m]
  ! hsn - snow mass [m]
  ! A - ice compactness
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! Ta - air temperature
  ! qa - specific humidity
  ! rain - precipitation rain
  ! snow - precipitation snow
  ! runo - runoff
  ! ug - wind speed
  ! ustar - friction velocity
  ! T_oc, S_oc - ocean temperature and salinity beneath the ice (mixed layer)
  ! H_ML - mixed layer depth - should be specified.
  ! t - temperature of snow/ice top surface
  ! dt - time step [s]
  ! ch - transfer coefficient for sensible heat (for open ocean)
  ! ce - transfer coefficient for evaporation   (for open ocean)
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice)  

  ! Output parameters:
  !-------------------
  ! h - ice mass
  ! hsn - snow mass
  ! A - ice compactness
  ! t - temperature of snow/ice top surface
  ! fw - freshwater flux due to ice melting [m water/dt]
  ! ehf - net heat flux at the ocean surface [W/m2]        !RTnew

  use i_therm_parms
  implicit none

  integer k
  real*8  h,hsn,A,fsh,flo,Ta,qa,rain,snow,runo,rsss
  real*8  ug,ustar,T_oc,S_oc,H_ML,t,dt,ch,ce,ch_i,ce_i,fw,ehf
  real*8  dhgrowth,dhsngrowth,ahf,prec,subli,subli_i,rsf
  real*8  rhow,show,rhice,shice,sh,thick,thact
  real*8  rh,rA,qhst,sn,hsntmp,o2ihf,evap
  real*8  iflice,hflatow,hfsenow,hflwrdout
  real*8, external  :: TFrez  ! Sea water freeze temperature.

  ! Store ice thickness at start of growth routine
  dhgrowth=h  	  

  ! determine h(i,j)/a(i,j) = actual ice thickness.
  ! if snow layer is present, add hsn weighted with quotient
  ! of conductivities of ice and snow, according to 0-layer approach
  ! of Semtner (1976).   	    
  ! thickness at the ice covered part
  thick=hsn*(con/consn)/max(A,Armin)    ! Effective snow thickness
  thick=thick+h/max(A,Armin)            ! Effective total snow-ice thickness

  ! Growth rate for ice in open ocean
  rhow=0.0
  evap=0.0
  call obudget(qa,fsh,flo,T_oc,ug,ta,ch,ce,rhow,evap,hflatow,hfsenow,hflwrdout) 
  hflatow=hflatow*(1.0-A)
  hfsenow=hfsenow*(1.0-A)
  hflwrdout=hflwrdout*(1.0-A)

  ! add heat loss at open ocean due to melting snow fall
  rhow=rhow+snow*1000.0/rhoice
  ! dt and (1-A) will be multiplied afterwards

  ! growth rate of ice in ice covered part
  ! following Hibler 1984
  ! assuming ice thickness has an euqal, 7-level distribution from zero to two times h 
  rhice=0.0                      
  subli=0.0
  if (thick.gt.hmin) then
     do k=1,iclasses			  
        thact = (2*k-1)*thick/float(iclasses)  	! Thicknesses of actual ice class
        call budget(thact,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,shice,subli_i) 
        !Thick ice K-class growth rate
        rhice=rhice+shice/float(iclasses)      	! Add to average heat flux
        subli=subli+subli_i/float(iclasses)
     end do
  end if

  ! Convert growth rates [m ice/sec] into growth per time step DT.
  rhow=rhow*dt
  rhice=rhice*dt

  ! Multiply ice growth of open water and ice
  ! with the corresponding areal fractions of grid cell
  show =rhow*(1.0-A)
  shice=rhice*A
  sh   =show+shice

  ! Store atmospheric heat flux, average over grid cell [W/m**2]
  ahf=-cl*sh/dt   

  ! precipitation (into the ocean)
  prec=rain+runo+snow*(1.0-A)  	        ! m water/s

  ! snow fall above ice
  hsn=hsn+snow*dt*A*1000.0/rhosno	! Add snow fall to temporary snow thickness    !!!
  dhsngrowth=hsn   		        ! Store snow thickness after snow fall 

  ! evap
  evap=evap*(1.0-A)    			! m water/s
  subli=subli*A

  ! If there is atmospheric melting, first melt any snow that is present.
  ! Atmospheric heat flux available for melting
  ! sh = MINUS atm. heat flux / specific latent heat of sea ice
  ! Note: (sh<0) for melting, (sh>0) for freezing
  hsntmp= -min(sh,0.0)*rhoice/rhosno

  ! hsntmp is the decrease in snow thickness due to atmospheric melting
  ! [m/DT]. Do not melt more snow than available
  hsntmp=min(hsntmp,hsn)
  hsn=hsn-hsntmp  ! Update snow thickness after atmospheric snow melt

  ! Negative atmospheric heat flux left after melting of snow
  ! Note: (sh<0) and (hsntmp>0) for melting conditions
  ! hsntmp=0 for non-snow-melting conditions
  rh=sh+hsntmp*rhosno/rhoice
  h=max(h,0.)

  ! Compute heat content qhst of mixed layer - sea ice system
  !
  ! Total heat content is the sum of
  !	h	ice thickness after calculation of dynamic effects
  !	rh	change in ice thickness due to atmospheric forcing
  ! and heat available in mixed layer, with
  !	T_oc	temperature of ocean surface layer
  !	Tfrez	freezing point of sea water
  !	H_ML	thickness of uppermost layer
  !
  !RT:
  ! There are three possibilities to do this.
  ! 1.: Assume an instantaneous adjustment of mixed layer heat content.
  !     Any heat available is then instantaneously used to melt ice.
  !     (so-called ice-bath approach)
  !     This is what used to be used in the Lemke sea ice-mixed layer model.
  ! rh=rh-(T_oc-TFrez(S_oc))*H_ML*cc/cl
  ! qhst=h+rh 
  !
  ! 2.: Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference. For a first step 
  !     we can assume a constant exchange coefficient gamma_t:
  ! o2ihf= (T_oc-TFrez(S_oc))*gamma_t*cc*A     &
  !        +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A) ! [W/m2]
  ! rh=rh-o2ihf*dt/cl
  ! qhst=h+rh		                      	! [m]
  !
  ! 3.  Parameterize the ocean-to-ice heat flux (o2ihf)
  !     as a function of temperature difference and the
  !     friction velocity:
  o2ihf= (T_oc-TFrez(S_oc))*0.006*ustar*cc*A  &
       +(T_oc-Tfrez(S_oc))*H_ML/dt*cc*(1.0-A)  	! [W/m2]
  rh=rh-o2ihf*dt/cl
  qhst=h+rh		              		! [m]

  ! Melt snow if there is any ML heat content left (qhst<0).
  ! This may be the case if advection moves ice (with snow) to regions
  ! with a warm mixed layer.
  sn=hsn+min(qhst,0.)*rhoice/rhosno

  ! New temporary snow thickness must not be negative:
  sn=max(sn,0.)

  ! Update snow and ice depth
  hsn=sn
  h=max(qhst,0.)     
  if (h.lt.1E-6) h=0.        ! Avoid very small ice thicknesses

  ! heat and fresh water fluxes
  dhgrowth=h-dhgrowth        ! Change in ice thickness due to thermodynamic effects
  dhsngrowth=hsn-dhsngrowth  ! Change in snow thickness due to thermodynamic melting
  ! (without snow fall). This is a negative value (MINUS snow melt)

  dhgrowth=dhgrowth/dt       ! Conversion: 'per time step' -> 'per second'
  dhsngrowth=dhsngrowth/dt   ! Conversion: 'per time step' -> 'per second'

  ! (radiation+turbulent) + freezing(-melting) sea-ice&snow 
  ehf = ahf + cl*(dhgrowth+(rhosno/rhoice)*dhsngrowth) 
  ! (prec+runoff)+evap - freezing(+melting) ice&snow
#ifdef use_fullfreesurf
  fw= prec+evap - dhgrowth*rhoice/rhowat - dhsngrowth*rhosno/rhowat 
  rsf= -dhgrowth*rhoice/rhowat*Sice
#else
  fw= prec+evap - dhgrowth*rhoice/rhowat*(rsss-Sice)/rsss - dhsngrowth*rhosno/rhowat 
#endif

  ! Changes in compactnesses (equation 16 of Hibler 1979)
  rh=-min(h,-rh)   ! Make sure we do not try to melt more ice than is available
  rA= rhow - o2ihf*dt/cl !Qiang: it was -(T_oc-TFrez(S_oc))*H_ML*cc/cl, changed in June 2010
  A=A + 0.5*min(rh,0.)*A/max(h,hmin) + max(rA,0.)*(1.-A)/h0   
  !meaning:       melting                  freezing

  A=min(A,h*1.e6)     ! A -> 0 for h -> 0
  A=min(max(A,0.),1.) ! A >= 0, A <= 1

  ! Flooding (snow to ice conversion)
  iflice=h
  call flooding(h,hsn)     
  iflice=(h-iflice)/dt

  ! to maintain salt conservation for the current model version
  !(a way to avoid producing net salt from snow-type-ice) 
#ifdef use_fullfreesurf
  rsf=rsf-iflice*rhoice/rhowat*Sice
#else
  fw=fw+iflice*rhoice/rhowat*Sice/rsss
#endif   

  ! add sublimation to evap
  evap=evap+subli  !negative up, m/s

  return
end subroutine therm_ice
!
!=====================================================================================
!
subroutine budget (hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh,subli)
  ! Thick ice growth rate [m ice/sec]
  !
  ! INPUT:
  ! hice - actual ice thickness [m]
  ! hsn - snow thickness, used for albedo parameterization [m]
  ! t - temperature of snow/ice surface [C]
  ! ta - air temperature [C]
  ! qa - specific humidity [Kg/Kg]
  ! fsh - shortwave radiation [W/m**2]
  ! flo - longwave radiation  [W/m**2]
  ! ug - wind speed [m/sec]
  ! S_oc - ocean salinity for the temperature of the ice base calculation [ppt]
  ! ch_i - transfer coefficient for sensible heat (for ice)
  ! ce_i - transfer coefficient for evaporation   (for ice) 
  !
  ! OUTPUT: fh - growth rate
  !
  ! qiang: The formular for saturated humidity was modified according to Large/Yeager2004
  ! to allow direct comparison with the CORE results (Griffies et al. 2009). The new
  ! formular does not require sea level pressure.
  ! A similar change was also made for the obudget routine.
  ! It was found through experiments that the results are quite similar to that from the
  ! original code, and the simulated ice volume is only slightly larger after modification. 
  
  use i_therm_parms
  implicit none

  integer iter, imax      ! Number of iterations
  real*8  hice,hsn,t,ta,qa,fsh,flo,ug,S_oc,ch_i,ce_i,fh
  real*8  hfsen,hfrad,hflat,hftot,subli         
  real*8  alb             ! Albedo of sea ice
  real*8  q1, q2	  ! coefficients for saturated specific humidity
  real*8  A1,A2,A3,B,C, d1, d2, d3   
  real*8, external :: TFrez

  data q1 /11637800.0/, q2 /-5897.8/ 
  data imax /5/

  ! set albedo
  ! ice and snow, freezing and melting conditions are distinguished.
  if (t<0.0) then	        ! freezing condition    
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsn         	
     else              		!   no snow cover       
        alb=albi       	
     endif
  else			        ! melting condition     
     if (hsn.gt.0.0) then	!   snow cover present  
        alb=albsnm	    	
     else			!   no snow cover       
        alb=albim		
     endif
  endif

  d1=rhoair*cpair*Ch_i
  d2=rhoair*Ce_i
  d3=d2*clhi

  ! total incoming atmospheric heat flux
  A1=(1.0-alb)*fsh + flo + d1*ug*ta + d3*ug*qa   ! in LY2004 emiss is multiplied wiht flo

  ! NEWTON-RHAPSON TO GET TEMPERATURE AT THE TOP OF THE ICE LAYER

  do iter=1,imax

     B=q1/rhoair*exp(q2/(t+tmelt))		! (saturated) specific humidity over ice
     A2=-d1*ug*t-d3*ug*B &
          -emiss_ice*boltzmann*((t+tmelt)**4)	! sensible and latent heat,and outward radiation
     A3=-d3*ug*B*q2/((t+tmelt)**2)		! gradient coefficient for the latent heat part
     C=con/hice                     		! gradient coefficient for downward heat conductivity
     A3=A3+C+d1*ug & 			! gradient coefficient for sensible heat and radiation 
          +4.0*emiss_ice*boltzmann*((t+tmelt)**3)    
     C=C*(TFrez(S_oc)-t)       		! downward conductivity term

     t=t+(A1+A2+C)/A3 		        ! NEW ICE TEMPERATURE AS THE SUM OF ALL COMPONENTS
  end do

  t=min(0.0,t)

  ! heat fluxes [W/m**2]:

  hfrad= (1.0-alb)*fsh &	        ! absorbed short wave radiation
       +flo &           	        ! long wave radiation coming in  ! in LY2004 emiss is multiplied
       -emiss_ice*boltzmann*((t+tmelt)**4) 	! long wave radiation going out

  hfsen=d1*ug*(ta-t) 			! sensible heat 
  subli=d2*ug*(qa-B) 			! sublimation
  hflat=clhi*subli                     	! latent heat

  hftot=hfrad+hfsen+hflat               ! total heat

  fh= -hftot/cl                         ! growth rate [m ice/sec]
  !                                      	+: ML gains energy, ice melts
  !                                      	-: ML loses energy, ice grows
  subli=subli/rhowat                    ! negative upward

  return
end subroutine budget
!
!======================================================================================
!
subroutine obudget (qa,fsh,flo,t,ug,ta,ch,ce,fh,evap,hflatow,hfsenow,hflwrdout)  
  ! Ice growth rate for open ocean [m ice/sec]
  !
  ! INPUT:
  ! t - temperature of open water [C]
  ! fsh - shortwave radiation
  ! flo - longwave radiation
  ! ta - air temperature [C]
  ! qa  - specific humidity             
  ! ug - wind speed [m/sec]
  ! ch - transfer coefficient for sensible heat
  ! ce - transfer coefficient for evaporation
  !
  ! OUTPUT: fh - growth rate
  !         evap - evaporation

  use i_therm_parms
  implicit none

  real*8 qa,t,Ta,fsh,flo,ug,ch,ce,fh,evap
  real*8 hfsenow,hfradow,hflatow,hftotow,hflwrdout,b
  real*8 q1, q2 		! coefficients for saturated specific humidity

  data q1 /640380./, q2 /-5107.4/

  ! (saturated) surface specific humidity
  b=0.98*q1/rhoair*exp(q2/(t+tmelt)) 		! LY2004 

  ! heat fluxes [W/m**2]:

  hfradow= (1.0-albw)*fsh &	            	! absorbed short wave radiation
       +flo             	                ! long wave radiation coming in !put emiss/check
  hflwrdout=-emiss_wat*boltzmann*((t+tmelt)**4) 	! long wave radiation going out !in LY2004 emiss=1
  hfradow=hfradow+hflwrdout

  hfsenow=rhoair*cpair*ch*ug*(ta-t)             ! sensible heat 
  evap=rhoair*ce*ug*(qa-b)			! evaporation kg/m2/s
  hflatow=clhw*evap                       	! latent heat W/m2

  hftotow=hfradow+hfsenow+hflatow             	! total heat W/m2

  fh= -hftotow/cl                             	! growth rate [m ice/sec]
  !                                           	+: ML gains energy, ice melts
  !                                           	-: ML loses energy, ice grows
  evap=evap/rhowat 	 			! evaporation rate [m water/s],negative up !!!

  return
end subroutine obudget
!
!======================================================================================
!
subroutine flooding (h,hsn)
  use i_therm_parms

  real*8 h,hsn,hdraft,hflood

  hdraft=(rhosno*hsn+h*rhoice)/rhowat ! Archimedes: displaced water
  hflood=hdraft-min(hdraft,h)         ! Increase in mean ice thickness due to flooding
  h=h+hflood                          ! Add converted snow to ice volume
  hsn=hsn-hflood*rhoice/rhosno        ! Subtract snow from snow layer

  !RT   This is what all AWI sea ice models do, but
  !RT   I wonder whether it really is correct for the heat budget.
  !RT   I suggest we initially keep it to allow for a comparison with BRIOS results
  !RT   and rethink it at a later stage.

  return
end subroutine flooding
!
!======================================================================================
!
function TFrez(S)
  ! Nonlinear correlation for the water freezing temperature.
  ! Millero (1978) - UNESCO. Reference - See A. Gill, 1982.
  implicit none
  real*8 S, TFrez

  TFrez= -0.0575*S+1.7105e-3 *sqrt(S**3)-2.155e-4 *S*S

end function TFrez



subroutine ice_matrices_setup
  use g_parfe
  implicit none

  call icestiff_matrix
  call icestiff_fill_tg  ! fill mass matrix and get lumped mass matrix 

end subroutine ice_matrices_setup
!
!-------------------------------------------------------------------------
!
subroutine icestiff_matrix
  ! Sets the structure of the stiffness matrix for ice advection
  use i_array
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use g_PARFE
  use g_config
  implicit none
  !
  integer                           :: k
  integer, allocatable              :: pnza(:), rpnza(:) 
  integer                           :: nini, nend
  !
  ! a)
  icestiff%dim=nod2D
  allocate(icestiff%rowptr(myDim_nod2D+1))         
  icestiff%rowptr(1)=1
  do k=1,myDim_nod2D                               
     icestiff%rowptr(k+1)=icestiff%rowptr(k)+nghbr_nod2D(k)%nmb
  end do
  icestiff%nza=icestiff%rowptr(myDim_nod2D+1)-1    
  !

  ! b)
  allocate(icestiff%colind(icestiff%nza))
  allocate(icestiff%values(icestiff%nza))
  icestiff%values=0.0
  !
  ! =================================            
  ! Exchange between PE to know their 
  ! numbers of non-zero entries (nza):
  ! =================================     
  allocate(pnza(npes),rpnza(npes))
  pnza(1:npes)=0
  pnza(mype+1)=icestiff%nza
  call MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call MPI_AllREDUCE( pnza, rpnza, &
       npes, MPI_INTEGER,MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)

  if (mype==0) then
     k=0
  else
     k=sum(rpnza(1:mype))  ! This is offset for mype 
  end if
  icestiff%rowptr=icestiff%rowptr+k       

  ! =================================
  ! replace local nza with a global one
  ! =================================
  icestiff%nza=sum(rpnza(1:npes))             
  deallocate(rpnza,pnza)

  !c)
  do k=1,myDim_nod2D                                 
     nini=icestiff%rowptr(k)-icestiff%rowptr(1)+1    
     nend=icestiff%rowptr(k+1)-icestiff%rowptr(1)    
     icestiff%colind(nini:nend)= &
          nghbr_nod2D(k)%addresses
  end do

  ! ==================================
  ! addresses are now in local numbering. We need to make them global
  ! and then local contiguous
  ! (i) global natural: 
  do k=1,icestiff%rowptr(myDim_nod2D+1)-icestiff%rowptr(1)   
     icestiff%colind(k)=myList_nod2D(icestiff%colind(k))     
  end do
  ! (ii) global PE contiguous: 
  ! mapping(1+nod3D:nod2D+nod3D) is 2D mapping
  do k=1,icestiff%rowptr(myDim_nod2D+1)-icestiff%rowptr(1)    
     icestiff%colind(k)=mapping(nod3D+icestiff%colind(k))     
  end do

end subroutine icestiff_matrix
!
!----------------------------------------------------------------
!
subroutine icestiff_fill_tg
  ! Fill in ice stiffness(mass) matrix for TG/FCT scheme
  use o_MESH
  use o_ELEMENTS
  use o_PARAM
  use i_array
  use g_parfe
  use g_config
  implicit none
  !
  integer                            :: i, q, row, col, ipos, offset
  integer                            :: elem, elnodes(3)
  real(kind=8)                       :: vol, inv12
  !
  inv12=1.0_8/12.0_8
  icestiff%values=0.0
  allocate(ice_lump(myDim_nod2d+eDim_nod2D))              
  ice_lump=0.0_8
  !
  do elem=1,myDim_elem2D                                 

     elnodes=elem2D_nodes(:,elem)
     vol=voltriangle(elem)*dt_inv*inv12
     do i=1,3             ! all rows into which the element elem could contribute
        row=elnodes(i)
	if(row>myDim_nod2D) cycle                         
        do q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        end do
        offset=icestiff%rowptr(row)-icestiff%rowptr(1)    
        do q=1,3          ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
	   icestiff%values(ipos)=icestiff%values(ipos)+vol
           if (q==i) then
              icestiff%values(ipos)=icestiff%values(ipos)+vol
           end if
        end do
        ice_lump(row)=ice_lump(row) + vol*4.0

     end do
  end do
  call com_2D(ice_lump)

end subroutine icestiff_fill_tg
!
!------------------------------------------------------------------------
!
#ifndef use_ice_gls
subroutine ice_rhs_tg
  use o_mesh
  use o_elements
  use o_param
  use i_dyn_parms
  use i_array
  use g_PARFE
  use g_config
  implicit none 
  !
  integer                      :: m, i, row, elnodes(3), elem, edgnod(2)
  real(kind=8)                 :: dx(3), dy(3), vol, aux, adv
  real(kind=8)                 :: uvel(3), vvel(3), usum, vsum, um, vm
  real(kind=8)                 :: m_el(3), a_el(3), ms_el(3)
  real(kind=8)                 :: mdiffx, mdiffy, adiffx, adiffy 
  real(kind=8)                 :: msdiffx, msdiffy, m_sum, a_sum, ms_sum
  real(kind=8)                 :: um_sum, vm_sum, ua_sum, va_sum
  real(kind=8)                 :: ums_sum, vms_sum, div
  real(kind=8)                 :: inv2, inv3, inv6, inv12
  real(kind=8)                 :: nv(2), uedg(2), vedg(2)
  real(kind=8)                 :: m_edg(2), a_edg(2), ms_edg(2)
  real(kind=8)                 :: aux1, aux2, aux3

  inv2=0.5_8
  inv3=1.0_8/3.0_8
  inv12=1.0_8/12.0_8
  do row=1, myDim_nod2d          

     rhs_m(row)=0.
     rhs_a(row)=0.
     rhs_ms(row)=0.
  enddo
  !
  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)
     dx=bafux_2d(:, elem)
     dy=bafuy_2d(:, elem)
     vol=voltriangle(elem)
     uvel=u_ice(elnodes)
     vvel=v_ice(elnodes)
     usum=sum(uvel)
     vsum=sum(vvel)
     um=inv3*usum
     vm=inv3*vsum	
     m_el=m_ice(elnodes)
     a_el=a_ice(elnodes)
     ms_el=m_snow(elnodes)
     mdiffx=sum(dx*m_el)
     mdiffy=sum(dy*m_el)
     adiffx=sum(dx*a_el)
     adiffy=sum(dy*a_el)
     msdiffx=sum(dx*ms_el)
     msdiffy=sum(dy*ms_el)
     m_sum=sum(m_el)
     a_sum=sum(a_el)
     ms_sum=sum(ms_el)
     um_sum=sum(m_el*uvel)
     vm_sum=sum(m_el*vvel)
     ua_sum=sum(a_el*uvel)
     va_sum=sum(a_el*vvel)
     ums_sum=sum(ms_el*uvel)
     vms_sum=sum(ms_el*vvel)
     div=sum(dx*uvel+dy*vvel)
     !
     !assembling over nodes
     do i=1,3 
        row=elnodes(i)

        !diffusion
        rhs_m(row)=rhs_m(row)-Kh_ice*(mdiffx*dx(i)+mdiffy*dy(i))*vol
        rhs_a(row)=rhs_a(row)-Kh_ice*(adiffx*dx(i)+adiffy*dy(i))*vol
        rhs_ms(row)=rhs_ms(row)-Kh_ice*(msdiffx*dx(i)+msdiffy*dy(i))*vol

        !advection
        adv=(dx(i)*(usum*m_sum+um_sum)+dy(i)*(vsum*m_sum+vm_sum))*inv12
        rhs_m(row)=rhs_m(row)+adv*vol
        adv=(dx(i)*(usum*a_sum+ua_sum)+dy(i)*(vsum*a_sum+va_sum))*inv12
        rhs_a(row)=rhs_a(row)+adv*vol
        adv=(dx(i)*(usum*ms_sum+ums_sum)+dy(i)*(vsum*ms_sum+vms_sum))*inv12
        rhs_ms(row)=rhs_ms(row)+adv*vol

        !TG stabilization
        aux=(um*dx(i)+vm*dy(i))*dt*inv2*vol
        adv=um*mdiffx+vm*mdiffy+div*m_sum*inv3
        rhs_m(row)=rhs_m(row)-aux*adv
        adv=um*adiffx+vm*adiffy+div*a_sum*inv3
        rhs_a(row)=rhs_a(row)-aux*adv
        adv=um*msdiffx+vm*msdiffy+div*ms_sum*inv3
        rhs_ms(row)=rhs_ms(row)-aux*adv
     end do
  end do
  !
#ifdef use_opbnd_restoring 
  inv6=1.0_8/6.0_8
  do m=1,nmbr_opbnd_edg
     edgnod=opbnd_edg(m,1:2)
     nv=opbnd_edg_nv(m,1:2)
     vol=opbnd_edg_nv(m,3)
     uedg=u_ice(edgnod)
     vedg=v_ice(edgnod)
     vedg=uedg*nv(1)+vedg*nv(2)
     if(vedg(1)<0.) vedg(1)=0.     ! assume no ice is entering the model domain
     if(vedg(2)<0.) vedg(2)=0.

     m_edg=m_ice(edgnod)
     a_edg=a_ice(edgnod)
     ms_edg=m_snow(edgnod)

     aux1=(vedg(1)*m_edg(2)+vedg(2)*m_edg(1)+sum(vedg*m_edg))*inv12
     aux2=(vedg(1)*a_edg(2)+vedg(2)*a_edg(1)+sum(vedg*a_edg))*inv12
     aux3=(vedg(1)*ms_edg(2)+vedg(2)*ms_edg(1)+sum(vedg*ms_edg))*inv12

     do i=1,2
        row=edgnod(i)
        rhs_m(row)=rhs_m(row)-(aux1+vedg(i)*m_edg(i)*inv6)*vol
        rhs_a(row)=rhs_a(row)-(aux2+vedg(i)*a_edg(i)*inv6)*vol
        rhs_ms(row)=rhs_ms(row)-(aux3+vedg(i)*ms_edg(i)*inv6)*vol
     end do
  end do
#endif

end subroutine ice_rhs_tg
#endif
!
!=========================================================================
!
#ifdef use_ice_gls
subroutine icestiff_fill_gls
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use i_dyn_parms
  use i_array
  use g_PARFE
  use g_config
  implicit none

  integer                           :: m, i, q, ipos, offset, edgnod(2)
  integer                           :: col, col2, row, elem, elnodes(3)
  real(kind=8)                      :: dx(3), dy(3), entries(3), entries_t(3)
  real(kind=8)                      :: v, uel(3), vel(3), usum, vsum, um, vm
  real(kind=8)                      :: div, aux, hh, v_nabla_psi, adv, diff
  real(kind=8)                      :: ept, ept_invdt, inv12, inv6, inv3, inv2
  real(kind=8)                      :: nv(2), uedg(2), vedg(2)
  real(kind=8)                      :: m_edg(2), a_edg(2), ms_edg(2)
  real(kind=8)                      :: aux1, aux2, aux3

  inv12=1.0_8/12.0_8
  inv3=1.0_8/3.0_8
  inv2=0.5_8
  do row=1,myDim_nod2d             
     rhs_m(row)=0.
     rhs_a(row)=0.
     rhs_ms(row)=0.
     col=icestiff%rowptr(row)-icestiff%rowptr(1)+1   
     col2=icestiff%rowptr(row+1)-icestiff%rowptr(1)  
     icestiff%values(col:col2)=0.0
  enddo

  do elem=1, myDim_elem2d           
     elnodes=elem2D_nodes(:,elem)
     v=voltriangle(elem)
     dx=bafux_2D(:,elem)
     dy=bafuy_2D(:,elem)
     uel=u_ice(elnodes)
     vel=v_ice(elnodes)
     usum=sum(uel)
     vsum=sum(vel)
     um=usum*inv3
     vm=vsum*inv3
     div=sum(uel*dx + vel*dy)
     aux=v*dt_inv*inv12   
     hh=sqrt(v)
     ept=0.1*dt_inv + (abs(um)+abs(vm))/hh + 2.0*Kh_ice/(hh*hh)
     ept=v/ept
     ept_invdt=ept*dt_inv

     do i=1,3             ! all rows into which the element elem could contribute
        row=elnodes(i)
        if(row>myDim_nod2D) cycle     
        do q=1,nghbr_nod2D(row)%nmb
           col_pos(nghbr_nod2D(row)%addresses(q))=q
        end do

        v_nabla_psi=um*dx(i)+vm*dy(i) 

        entries=0.0
        do q=1,3          ! all columns 

           diff=Kh_ice*(dx(i)*dx(q)+dy(i)*dy(q))                ! pure numerical things
           adv=-(dx(i)*(usum+uel(q))+dy(i)*(vsum+vel(q)))*inv12        
	   entries(q)=(adv+diff)*v                              ! advection + numerical diff

	   adv= um*dx(q) + vm*dy(q) + div*inv3 
           entries(q)=entries(q)+ ept*v_nabla_psi*adv	        ! stabilization

           col=elnodes(q)
           rhs_m(row)=rhs_m(row)-entries(q)*m_ice(col)		
           rhs_a(row)=rhs_a(row)-entries(q)*a_ice(col)		
           rhs_ms(row)=rhs_ms(row)-entries(q)*m_snow(col)       ! fill rhs

           entries_t(q)=aux                             	! mass matrix

           entries_t(q)=entries_t(q)+ept_invdt*v_nabla_psi*inv3	! stabilization, t          

        end do

        entries_t(i)=entries_t(i)+aux                    	! completes mass matrix 
        entries=inv2*entries+entries_t       

        ! put the entries to the appropriate place
        offset=icestiff%rowptr(row)-icestiff%rowptr(1)    
        do q=1,3        ! all columns 
           col=elnodes(q)
           ipos = offset+col_pos(col)
           icestiff%values(ipos)=icestiff%values(ipos)+entries(q)
        end do
     end do
  end do

#ifdef use_opbnd_restoring 
  inv6=1.0_8/6.0_8
  do m=1,nmbr_opbnd_edg
     edgnod=opbnd_edg(m,1:2)
     nv=opbnd_edg_nv(m,1:2)
     v=opbnd_edg_nv(m,3)
     uedg=u_ice(edgnod)
     vedg=v_ice(edgnod)
     vedg=uedg*nv(1)+vedg*nv(2)
     if(vedg(1)<0.) vedg(1)=0.     ! assume no ice coming into the model domain
     if(vedg(2)<0.) vedg(2)=0.

     m_edg=m_ice(edgnod)
     a_edg=a_ice(edgnod)
     ms_edg=m_snow(edgnod)

     aux1=(vedg(1)*m_edg(2)+vedg(2)*m_edg(1)+sum(vedg*m_edg))*inv12
     aux2=(vedg(1)*a_edg(2)+vedg(2)*a_edg(1)+sum(vedg*a_edg))*inv12
     aux3=(vedg(1)*ms_edg(2)+vedg(2)*ms_edg(1)+sum(vedg*ms_edg))*inv12

     do i=1,2
        row=edgnod(i)
        rhs_m(row)=rhs_m(row)-(aux1+vedg(i)*m_edg(i)*inv6)*v
        rhs_a(row)=rhs_a(row)-(aux2+vedg(i)*a_edg(i)*inv6)*v
        rhs_ms(row)=rhs_ms(row)-(aux3+vedg(i)*ms_edg(i)*inv6)*v
     end do
  end do
#endif

end subroutine icestiff_fill_gls
#endif
!
!----------------------------------------------------------------------------
!
subroutine ice_solve
  use o_mesh
  use i_array
  use o_array
  use g_parfe
  use i_dyn_parms
  use g_config
  implicit none
  !
  integer                                 :: n, row, clo, clo2, cn, location(100)
  real(kind=8)                            :: rhs_new
  real(kind=8), allocatable               :: auxarray(:,:)

  allocate(auxarray(3,myDim_nod2d))

  !the first approximation
  do row=1,myDim_nod2D                  
     dm_ice(row)=rhs_m(row)/ice_lump(row)
     da_ice(row)=rhs_a(row)/ice_lump(row)
     dm_snow(row)=rhs_ms(row)/ice_lump(row)
  end do

  call com_2D(dm_ice)
  call com_2D(da_ice)
  call com_2D(dm_snow)

  !iterate 
  do n=1,num_iter_solve_ice-1
     do row=1,myDim_nod2D               
        clo=icestiff%rowptr(row)-icestiff%rowptr(1)+1  
        clo2=icestiff%rowptr(row+1)-icestiff%rowptr(1) 
        cn=clo2-clo+1
        location(1:cn)=nghbr_nod2D(row)%addresses      
        rhs_new=rhs_m(row) - sum(icestiff%values(clo:clo2)*dm_ice(location(1:cn)))
        auxarray(1,row)=dm_ice(row)+rhs_new/ice_lump(row)    
        rhs_new=rhs_a(row) - sum(icestiff%values(clo:clo2)*da_ice(location(1:cn)))
        auxarray(2,row)=da_ice(row)+rhs_new/ice_lump(row)     
        rhs_new=rhs_ms(row) - sum(icestiff%values(clo:clo2)*dm_snow(location(1:cn)))
        auxarray(3,row)=dm_snow(row)+rhs_new/ice_lump(row)  
     end do

     do row=1,myDim_nod2D             
        dm_ice(row)=auxarray(1,row)      
        da_ice(row)=auxarray(2,row)      
	dm_snow(row)=auxarray(3,row)    
     end do
     call com_2D(dm_ice)
     call com_2D(da_ice)
     call com_2D(dm_snow)
  end do

  deallocate(auxarray)
end subroutine ice_solve
subroutine fct_ice_init
  use o_mesh
  use i_array
  use g_parfe
  implicit none
  integer        :: n2

  n2=myDim_nod2D+eDim_nod2D   

  allocate(m_icel(n2), a_icel(n2), m_snowl(n2))  ! low-order solutions
  allocate(icefluxes(myDim_elem2D,3))
  allocate(icepplus(n2), icepminus(n2))

  m_icel=0.0
  a_icel=0.0 
  m_snowl=0.0

end subroutine fct_ice_init
!
!----------------------------------------------------------------------------
!
subroutine fct_ice_solve
  ! driving routine of the ice fct scheme
  use i_array
  use i_solver
  use i_dyn_parms
  use g_PARFE
  implicit none

  if(lump_ice_matrix) then
     call ice_solve    
  else                            
     call solveIce(solve_m_ice)       
     call solveIce(solve_a_ice)       
     call solveIce(solve_m_snow)
     call com_2D(dm_ice)
     call com_2D(da_ice)
     call com_2D(dm_snow)      
  end if

  call ice_solve_low_order

  call ice_fem_fct(1)    ! m_ice   
  call ice_fem_fct(2)    ! a_ice
  call ice_fem_fct(3)    ! m_snow

  call com_2D(m_ice)
  call com_2D(a_ice)
  call com_2D(m_snow)
end subroutine fct_ice_solve
!
!----------------------------------------------------------------------------
!
subroutine ice_solve_low_order
  ! Low-order solution
  ! It is assumed that m_ice, a_ice and m_snow from the previous time step 
  ! are known at 1:myDim_nod2D+eDim_nod2D.
  ! One adds diffusive contribution to the rhs. It is realized as
  ! difference between the consistent and lumped mass matrices
  ! acting on the field from the previous time step. The mass matrix on the 
  ! lhs is replaced with lumped one.   
  use o_mesh
  use i_array
  use i_dyn_parms
  use g_parfe
  implicit none

  integer      :: row, clo, clo2, cn, location(100)
  real(kind=8) :: gamma

  gamma=ice_gamma_fct     
  do row=1,myDim_nod2D             
     clo=icestiff%rowptr(row)-icestiff%rowptr(1)+1  
     clo2=icestiff%rowptr(row+1)-icestiff%rowptr(1) 
     cn=clo2-clo+1
     location(1:cn)=nghbr_nod2D(row)%addresses      
     m_icel(row)=(rhs_m(row)+gamma*sum(icestiff%values(clo:clo2)* &
          m_ice(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*m_ice(row)
     a_icel(row)=(rhs_a(row)+gamma*sum(icestiff%values(clo:clo2)* &
          a_ice(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*a_ice(row)
     m_snowl(row)=(rhs_ms(row)+gamma*sum(icestiff%values(clo:clo2)* &
          m_snow(location(1:cn))))/ice_lump(row) + &
          (1.-gamma)*m_snow(row)
  end do

  call com_2D(m_icel)
  call com_2D(a_icel)
  call com_2D(m_snowl)
end subroutine ice_solve_low_order
!
!----------------------------------------------------------------------------
!
subroutine ice_fem_fct(tr_array_id)
  ! Flux corrected transport algorithm for ice advection
  !
  ! It is based on Loehner et al. (Finite-element flux-corrected 
  ! transport (FEM-FCT) for the Euler and Navier-Stokes equation, 
  ! Int. J. Numer. Meth. Fluids, 7 (1987), 1093--1109) as described by Kuzmin and
  ! Turek. (kuzmin@math.uni-dortmund.de) 
  !
  use o_mesh
  use o_elements
  use i_array
  use i_dyn_parms
  use g_parfe
  use g_config
  implicit none

  integer        :: tr_array_id
  integer        :: n, q, elem, elnodes(3), row
  real(kind=8), allocatable, dimension(:) :: tmax, tmin 
  real(kind=8)   :: vol, flux, ae, gamma, inv12, icoef(3,3)

  inv12=1.0_8/12.0_8

  gamma=ice_gamma_fct         

  !==========================
  ! Compute elemental antidiffusive fluxes to nodes
  !==========================
  ! This is the most unpleasant part --- 
  ! it takes memory and time. For every element 
  ! we need its antidiffusive contribution to 
  ! each of its 3 nodes

  allocate(tmax(myDim_nod2D), tmin(myDim_nod2D))

  ! Auxiliary elemental operator (mass matrix- lumped mass matrix)
  icoef=1.0
  do n=1,3   ! three upper nodes
     icoef(n,n)=-2.0
  end do
  icoef=icoef*inv12

  do elem=1, myDim_elem2D     
     elnodes=elem2D_nodes(:,elem)
     vol=voltriangle(elem)*dt_inv

     if (tr_array_id==1) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_ice(elnodes) + & 
                dm_ice(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if

     if (tr_array_id==2) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*a_ice(elnodes) + &  
                da_ice(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if

     if (tr_array_id==3) then
        do q=1,3       
           icefluxes(elem,q)=-sum(icoef(:,q)*(gamma*m_snow(elnodes) + &  
                dm_snow(elnodes)))*vol/ice_lump(elnodes(q))
        end do
     end if
  end do

  !==========================   
  ! Screening the low-order solution
  !==========================
  ! TO BE ADDED IF FOUND NECESSARY
  ! Screening means comparing low-order solutions with the
  ! solution on the previous time step and using whichever 
  ! is greater/smaller in computations of max/min below

  !==========================
  ! Cluster min/max
  !==========================
  if (tr_array_id==1) then
     do row=1, myDim_nod2D       
        n=nghbr_nod2D(row)%nmb
        tmax(row)=maxval(m_icel(nghbr_nod2D(row)%addresses(1:n)))   
        tmin(row)=minval(m_icel(nghbr_nod2D(row)%addresses(1:n)))    
        ! Admissible increments
        tmax(row)=tmax(row)-m_icel(row)                             
        tmin(row)=tmin(row)-m_icel(row)                             
     end do
  end if

  if (tr_array_id==2) then
     do row=1, myDim_nod2D          
        n=nghbr_nod2D(row)%nmb
        tmax(row)=maxval(a_icel(nghbr_nod2D(row)%addresses(1:n)))  
        tmin(row)=minval(a_icel(nghbr_nod2D(row)%addresses(1:n)))  
        ! Admissible increments
        tmax(row)=tmax(row)-a_icel(row)                            
        tmin(row)=tmin(row)-a_icel(row)                           
     end do
  end if

  if (tr_array_id==3) then
     do row=1, myDim_nod2D       
        n=nghbr_nod2D(row)%nmb
        tmax(row)=maxval(m_snowl(nghbr_nod2D(row)%addresses(1:n))) 
        tmin(row)=minval(m_snowl(nghbr_nod2D(row)%addresses(1:n)))  
        ! Admissible increments
        tmax(row)=tmax(row)-m_snowl(row)                         
        tmin(row)=tmin(row)-m_snowl(row)                            
     end do
  end if


  !=========================
  ! Sums of positive/negative fluxes to node row
  !=========================

  icepplus=0.
  icepminus=0.
  do elem=1, myDim_elem2D          
     elnodes=elem2D_nodes(:,elem)
     do q=1,3
        n=elnodes(q) 
        flux=icefluxes(elem,q)     
        if (flux>0) then
           icepplus(n)=icepplus(n)+flux
        else
           icepminus(n)=icepminus(n)+flux	  
        end if
     end do
  end do

  !========================
  ! The least upper bound for the correction factors
  !========================
  do n=1,myDim_nod2D              
     flux=icepplus(n)
     if (abs(flux)>0) then
        icepplus(n)=min(1.0,tmax(n)/flux)     
     else
        icepplus(n)=0.
     end if

     flux=icepminus(n)
     if (abs(flux)>0) then
        icepminus(n)=min(1.0,tmin(n)/flux)      
     else
        icepminus(n)=0.
     end if
  end do
  ! pminus and pplus are to be known to neighbouting PE
  call com_2D(icepminus)
  call com_2D(icepplus) 

  !========================	 
  ! Limiting
  !========================	 
  do elem=1, myDim_elem2D                                          
     elnodes=elem2D_nodes(:,elem)
     ae=1.0
     do q=1,3
        n=elnodes(q)  
        flux=icefluxes(elem,q)     
        if(flux>=0.) ae=min(ae,icepplus(n))
        if(flux<0.) ae=min(ae,icepminus(n))
     end do
     icefluxes(elem,:)=ae*icefluxes(elem,:) 
     !if (ae.le.0.0) write (*,*) 'ae is too large', ae 
  end do


  !==========================
  ! Update the solution 
  !==========================
  if(tr_array_id==1) then
     do n=1,myDim_nod2D          
        m_ice(n)=m_icel(n)
     end do
     do elem=1, myDim_elem2D      
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           m_ice(n)=m_ice(n)+icefluxes(elem,q)   
        end do
     end do
  end if

  if(tr_array_id==2) then
     do n=1,myDim_nod2D           
        a_ice(n)=a_icel(n)
     end do
     do elem=1, myDim_elem2D      
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           a_ice(n)=a_ice(n)+icefluxes(elem,q)    
        end do
     end do
  end if

  if(tr_array_id==3) then
     do n=1,myDim_nod2D         
        m_snow(n)=m_snowl(n)
     end do
     do elem=1, myDim_elem2D     
        elnodes=elem2D_nodes(:,elem)
        do q=1,3
           n=elnodes(q)  
           m_snow(n)=m_snow(n)+icefluxes(elem,q)    
        end do
     end do
  end if

  deallocate(tmin, tmax)
end subroutine ice_fem_fct
subroutine stress_tensor
  !
  ! In contrast to the 'old' version of this subroutine, this version 
  ! does not differentiate 
  ! the metric terms. This is not quite exact, but it saves a lot
  ! of computing time. Sea ice dynamics are an estimate anyway,
  ! so we do not spoil beautiful physics by simplifying the
  ! approximation a bit.
  !
  ! idea and coding: Sergey Danilov
  ! reasoning and irony: Ralph Timmermann   March 3.06
  !===================================================================
  
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use o_mesh
  use o_elements 
  use o_param
  use g_parfe
  use g_config
  implicit none

  integer        :: i, elem, elnodes(3)
  real(kind=8)   :: dx(3), dy(3)
  real(kind=8)   :: eps11, eps12, eps22, eta, xi, pressure, delta, aa
  real(kind=8)   :: elcos(3), val3, meancos, usum, vsum, vale, xi_limit
  real(kind=8)   :: det1, det2, r1, r2, r3, si1, si2, dte

  Tevp_inv=dt_inv*real(evp_Tdamp_ratio)
  val3=1.0_8/3.0_8
  vale=1.0_8/(ellipse**2)

  if(EVP_rheology) then 
     dte=dt/(1.0_8*evp_rheol_steps)
     det1=1.0_8+0.5_8*Tevp_inv*dte
     det2=1.0_8+0.5_8*Tevp_inv*dte*ellipse**2    !RTSD corrected 8.3.2006
     det1=1.0_8/det1
     det2=1.0_8/det2
  else
     xi_limit=0.5*rhoice*9.e8*vp_rheol_steps/dt 
     ! Limit moduli to satisfy the CFL - type criterium in Euler
     ! forward time stepping;  9.0e8 is the squared mesh step ??
  end if


  do elem=1,myDim_elem2d                    
     elnodes=elem2D_nodes(:,elem)

     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0.0) cycle

     dx=bafux_2d(:,elem)
     dy=bafuy_2d(:,elem)     
     !     vsum=sum(v_ice(elnodes))
     !     usum=sum(u_ice(elnodes))
     !     elcos=cos(coord_nod2D(2,elnodes))
     !     meancos=cos_elem2D(elem)

     ! Deformation rate tensor on element elem:
     eps11=sum(dx*u_ice(elnodes))
     !     eps11=eps11+val3*vsum*sum(dy*elcos)/meancos
     eps22=sum(dy*v_ice(elnodes))
     eps12=0.5_8*sum(dy*u_ice(elnodes) + dx*v_ice(elnodes))
     !     eps12=eps12-0.5_8*val3*usum*sum(dy*elcos)/meancos  

     ! moduli:
     delta=(eps11**2+eps22**2)*(1.0_8+vale)+4.0_8*vale*eps12**2 + &
          2.0_8*eps11*eps22*(1.0_8-vale)
     delta=sqrt(delta)
    
     vsum=sum(m_ice(elnodes))*val3
     usum=sum(a_ice(elnodes))*val3

     pressure=pstar*vsum*exp(-c_pressure*(1.0_8-usum))
    
     if (.not.EVP_rheology) then   
        xi=0.5_8*pressure/(delta+delta_min)

        if(xi>xi_limit*vsum) then              ! Limiting moduli 
           xi=xi_limit*vsum
        end if

        eta=xi*vale
        pressure=0.5_8*pressure*delta/(delta+delta_min)

        ! stress arrays 
        sigma11(elem)=(xi-eta)*(eps11+eps22)-0.5_8*pressure
        sigma12(elem)=2.0_8*eta*eps12
        sigma22(elem)=2.0_8*eta*eps22+sigma11(elem)
        sigma11(elem)=2.0_8*eta*eps11+sigma11(elem) 
     else
        pressure=0.5_8*pressure*Tevp_inv
        pressure=pressure*delta
        delta=1.0_8/(delta+delta_min)
        pressure=pressure*delta

        r1=pressure*((eps11+eps22)*delta -1.0_8) 
        r2=pressure*(eps11-eps22)*delta
        r3=pressure*eps12*delta
        si1=sigma11(elem)+sigma22(elem)
        si2=sigma11(elem)-sigma22(elem)

        si1=det1*(si1+dte*r1)
        si2=det2*(si2+dte*r2)
        sigma12(elem)=det2*(sigma12(elem)+dte*r3)
        sigma11(elem)=0.5_8*(si1+si2)
        sigma22(elem)=0.5_8*(si1-si2)
     end if
  end do

end subroutine stress_tensor
!
!===================================================================
!
subroutine stress2rhs
  ! Sergey Danilov, March 2006
  ! Qiang, 16.12.2010: add the nonlinear free surface option
  use o_MESH
  use o_ELEMENTS
  use o_param
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use g_PARFE
  use g_clock
  use g_config
  implicit none
  !
  integer       :: i ,k, row, elem, elnodes(3)
  real(kind=8)  :: dx(3), dy(3), vol
  real(kind=8)  :: val3, val12, meancos, elcos(3), aa, aux(3)
  real(kind=8)  :: mass, cluster_area, elevation_elem(3)

  val3=1.0_8/3.0_8
  val12=1.0_8/12.0_8

  do row=1, myDim_nod2d 
     rhs_u(row)=0.0
     rhs_v(row)=0.0
     rhs_a(row)=0.0
     rhs_m(row)=0.0
  end do

  do elem=1,myDim_elem2d         
     elnodes=elem2D_nodes(:,elem)

     aa=product(m_ice(elnodes))*product(a_ice(elnodes))
     if (aa==0.0) cycle

     vol=voltriangle(elem)
     dx=bafux_2d(:,elem)
     dy=bafuy_2d(:,elem)     
     elevation_elem=elevation(elnodes)
#ifdef use_fullfreesurf
     aux=(rhoice*m_ice(elnodes)+rhosno*m_snow(elnodes))*rho0r
     do i=1,3
        aux(i)=min(aux(i),max_ice_loading)
     end do
     elevation_elem=elevation_elem+aux
     !rho0r is used here to be consistent with the computation of ocean hydrostatic pressure
#endif

     !elcos=cos(coord_nod2D(2,elnodes))
     !meancos=cos_elem2D(elem)

     do k=1,3
        row=elnodes(k)
        rhs_u(row)=rhs_u(row) - vol* &
             (sigma11(elem)*dx(k)+sigma12(elem)*dy(k))   ! ) &
        !+val3*sum(dy*elcos)/meancos))
        rhs_v(row)=rhs_v(row) - vol* &
             (sigma12(elem)*dx(k)+sigma22(elem)*dy(k))   ! ) &
        ! +sigma11(elem)*val3*sum(dy*elcos)/meancos)
     end do

     ! use rhs_m and rhs_a for storing the contribution from elevation:
     aa=g*val3*vol
     do k=1,3
        row=elnodes(k)
        rhs_a(row)=rhs_a(row)-aa*sum(dx*elevation_elem)	    
        rhs_m(row)=rhs_m(row)-aa*sum(dy*elevation_elem)
     end do
  end do

  do row=1, myDim_nod2d               
     cluster_area=ice_lump(row)*dt  !note: ice_lump contains dt_inv
     mass=cluster_area*(m_ice(row)*rhoice+m_snow(row)*rhosno)

     if (mass>1.0e-4) then
        rhs_u(row)=rhs_u(row)/mass + rhs_a(row)/cluster_area 
        rhs_v(row)=rhs_v(row)/mass + rhs_m(row)/cluster_area 
     else
        rhs_u(row)=0.  
        rhs_v(row)=0.
     end if
  end do

  ! Clean neighbours:
  do row=myDim_nod2d+1,eDim_nod2d+myDim_nod2d  
     rhs_u(row)=0. 
     rhs_v(row)=0.
  end do

end subroutine stress2rhs
!
!===================================================================
!
subroutine rheology
  use o_param
  use o_MESH
  use o_ELEMENTS
  use o_array, only: coriolis_param_nod2D
  use i_dyn_parms
  use i_therm_parms
  use i_array
  use g_parfe
  use g_config

  implicit none
  integer         :: steps, shortstep, i, j
  real(kind=8)    :: rdt, drag, det, fc
  real(kind=8)    :: thickness, inv_thickness, umod, rhsu, rhsv


  if(.not.EVP_rheology) then
     rdt=dt/(1.0*vp_rheol_steps)
     steps=vp_rheol_steps
  else
     rdt=dt/(1.0*evp_rheol_steps)
     steps=evp_rheol_steps
  end if

  do shortstep=1, steps 
     ! ===== Boundary condition      
     do i=1, myDim_nod2D+eDim_nod2D    
        j=nod3D_below_nod2D(1,i)      
        if(index_nod3D(j)==11) then     
           u_ice(i)=0.0
           v_ice(i)=0.0
        end if
     end do

     call stress_tensor
     call stress2rhs
     do i=1,myDim_nod2d                 
        j=nod3D_below_nod2d(1,i)        
        if (index_nod3D(j)==11) cycle       ! Skip boundary nodes
        if (a_ice(i) > Armin) then          ! If ice is present, update velocities
           thickness=(rhoice*m_ice(i)+rhosno*m_snow(i))/a_ice(i)
           thickness=max(thickness, 90.0)   ! Limit the weighted thickness if it is too small
           inv_thickness=1.0/thickness

           umod=sqrt((u_ice(i)-u_w(i))**2+(v_ice(i)-v_w(i))**2)
           drag=Cd_oce_ice*umod*rho0*inv_thickness

           !rhs for water stress, air stress, and rhs_u/v (internal stress + ssh)
           rhsu=u_ice(i)+drag*rdt*u_w(i)+rdt*(inv_thickness*stress_atmice_x(i)+rhs_u(i))
           rhsv=v_ice(i)+drag*rdt*v_w(i)+rdt*(inv_thickness*stress_atmice_y(i)+rhs_v(i))

           !solve (coriolis and water stress are treated implicitly)
           fc=coriolis_param_nod2D(i)
           det=(1.+drag*rdt)**2+(rdt*fc)**2
           det=1.0_8/det
           u_ice(i)=det*((1.0+drag*rdt)*rhsu+rdt*fc*rhsv)
           v_ice(i)=det*((1.0+drag*rdt)*rhsv-rdt*fc*rhsu)

           !elseif(a_ice(i) > 0.001)          ! Set ice velocity equal to water velocity
           !   u_ice(i)=u_w(i)             
           !   v_ice(i)=v_w(i)
        end if
     end do
     call com_2D(u_ice)
     call com_2D(v_ice)
  end do
end subroutine rheology
!
!===================================================================
!
subroutine cut_off
  use i_array
  use i_therm_parms
  implicit none
  !
  real(kind=8)   :: m_min, a_min 
  !
  a_min=1.0e-3
  m_min=1.0e-3

  !where ((m_ice< m_min).or.(a_ice<a_min))
  !   a_ice=0.
  !   m_ice=0.
  !end where

  where(a_ice>1.0_8)
     a_ice=1.0_8
  end where

  where(a_ice<0.0_8)
     a_ice=0.0_8
  end where
end subroutine cut_off
!
!===================================================================


subroutine solveIce(ident)
  use g_PARFE
  use o_mesh
  use o_solver
  USE i_array
  use i_solver
  implicit none

#include "petscf.h"
  !include "pilutf.h"
  !include "hypref.h"

  integer        :: ident
  INTEGER        :: Pmode, restart, maxiter, lutype, fillin
  INTEGER        :: COMM, myrows
  REAL(kind=8)   :: droptol, soltol, rinfo(20,20)
  save rinfo

  COMM=MPI_COMM_WORLD

  maxiter=2000
  restart=15
  fillin=2
  lutype=2
  droptol=1.e-6
  soltol=1.e-6

  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

  if (ident==solve_m_ice) then ! m_ice
#ifdef use_ice_gls
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB+PET_PMVALS+PET_RCM + PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT 
#else
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB + PET_QUIET !+PET_REPORT
     if (iter_first) Pmode=Pmode+PET_STRUCT+PET_PMVALS+PET_RCM
#endif
     call PETSC_S(Pmode,ident,icestiff%dim, icestiff%nza,  myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &  
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_m, dm_ice, rinfo(:,12), COMM)
  endif

  if (ident==solve_a_ice) then ! a_ice
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB +PET_QUIET !+PET_REPORT
     call PETSC_S(Pmode, &
          ident-1, &
          icestiff%dim, icestiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &   
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_a, da_ice,  rinfo(:,13), COMM)
  endif

  if (ident==solve_m_snow) then ! m_snow
     Pmode=PET_BLOCKP+PET_SOLVE+PET_BICGSTAB +PET_QUIET !+PET_REPORT
     call PETSC_S(Pmode, &
          ident-2, &
          icestiff%dim, icestiff%nza, myrows, &
          maxiter, & 
          restart, &  
          fillin,  &  
          droptol, &  
          soltol,  &  
          part2D, icestiff%rowptr, icestiff%colind, icestiff%values, &
          rhs_ms, dm_snow,  rinfo(:,14), COMM)
  endif
end subroutine solveIce

! ==================================================================
subroutine par_init      ! initializes MPI
  use g_PARFE
  implicit none

  integer :: i

  call MPI_INIT(i)
  call MPI_Comm_Size(MPI_COMM_WORLD,npes,i)
  call MPI_Comm_Rank(MPI_COMM_WORLD,mype,i)

  call PETSCInitialize(PETSC_NULL_CHARACTER,i)

  if(mype==0) write(*,*) 'MPI has been initialized'
end subroutine par_init
!
!============================================================================
!
subroutine par_ex       ! finalizes MPI
  use g_PARFE

  call  MPI_Barrier(MPI_COMM_WORLD,MPIerr)
  call  MPI_Finalize(MPIerr)

end subroutine par_ex
!
!============================================================================
!
subroutine set_par_support_ini
  ! used during preparing distributed memory
  use o_MESH
  use o_elements
  use o_matrices
  use g_PARFE 
  use o_MATRICES
  implicit none
  !
  integer   n, j, k, nini, nend, count, dim_array

  ! Define partitioning vector

  allocate(part2D(nod2D))
  allocate(part3D(nod3D))

  do n=0, npes-1
     nini=(nod2D/npes)*n+1
     nend=(nod2D/npes)*(n+1)
     if (n==npes-1) nend=nod2D
     part2D(nini:nend)=n
  end do
  call partit(sshstiff%dim,sshstiff%rowptr,sshstiff%colind,npes, part2D)

  ! propagate partitioning vector down:
  do n=1,nod2D
     do j=1,num_layers_below_nod2D(n)+1
        k = nod3D_below_nod2D(j,n)
        part3D(k)=part2D(n)
     end do
  end do

  call communication_nod
  call mymesh


  if(mype==0) write(*,*) 'Communication arrays are set up'   
end subroutine set_par_support_ini

!=======================================================================

subroutine set_par_support
  ! use during model run
  use o_MESH
  use g_PARFE 
  use o_MATRICES
  implicit none

  integer   n, cnt
  ! In the distributed memory version, most of the job is already done 
  ! at the initialization phase and is taken into account in read_mesh
  ! routine. Here only communication buffers are set. 

  ! Allocate communication buffers: 

  if (npes>1) then
     allocate(s_buff_2d(com_nod2D%sPEnum),r_buff_2d(com_nod2D%sPEnum))
     do n=1, com_nod2D%sPEnum
        cnt=com_nod2D%sptr(n+1) - com_nod2D%sptr(n)
        allocate(s_buff_2d(n)%array(cnt))
     end do
     do n=1, com_nod2D%rPEnum
        cnt=com_nod2D%rptr(n+1) - com_nod2D%rptr(n)
        allocate(r_buff_2d(n)%array(cnt))
     end do

     allocate(s_buff_3d(com_nod3D%sPEnum),r_buff_3d(com_nod3D%sPEnum))
     do n=1, com_nod3D%sPEnum
        cnt=com_nod3D%sptr(n+1) - com_nod3D%sptr(n)
        allocate(s_buff_3d(n)%array(cnt))
     end do
     do n=1, com_nod3D%rPEnum
        cnt=com_nod3D%rptr(n+1) - com_nod3D%rptr(n)
        allocate(r_buff_3d(n)%array(cnt))
     end do
  end if

  if(mype==0) write(*,*) 'Communication buffer arrays are set up' 
end subroutine set_par_support

!
!============================================================================
!
subroutine communication_nod
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  !
  integer n,np, nz, prank, elem, elnodes(3), epe(3), counter, nini, nend
  integer, allocatable :: aux(:,:), pnum(:,:)

  ! Assume we have 2D partitioning vector in part. Find communication
  ! rules
  allocate(aux(npes,nod2D))
  allocate(pnum(npes,npes))
  aux=0
  pnum=0
  do elem=1,elem2D
     elnodes=elem2D_nodes(:,elem)
     epe=part2D(elnodes)+1
     if(epe(1).ne.epe(2)) then
        aux(epe(1), elnodes(2))=1
        aux(epe(2), elnodes(1))=1
     end if

     if(epe(2).ne.epe(3)) then
        aux(epe(3), elnodes(2))=1
        aux(epe(2), elnodes(3))=1
     end if
     if(epe(1).ne.epe(3)) then
        aux(epe(1), elnodes(3))=1
        aux(epe(3), elnodes(1))=1
     end if
  end do

  do n=1, nod2D
     do np=1, npes
        if(aux(np,n).ne.0) then 
           pnum(np,part2D(n)+1)=pnum(np,part2D(n)+1)+1
        end if
     end do
  end do

  ! We know how many external nodes each PE needs
  ! This is the 'receive' list   
  ! com_nod2D for 2D nodes

  ! The number of external PE I receive information from
  com_nod2D%rPEnum=0
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        com_nod2D%rPEnum=com_nod2D%rPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%rPE(com_nod2D%rPEnum))
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rPE(counter)=n-1
     end if
  end do


  ! Ptr to list of external nodes ordered by external PE ranks

  counter=0
  allocate(com_nod2D%rptr(com_nod2D%rPEnum+1)) 
  com_nod2D%rptr(1)=1
  do n=1, npes
     if(pnum(mype+1,n).ne.0) then
        counter=counter+1
        com_nod2D%rptr(counter+1)=com_nod2D%rptr(counter)+ pnum(mype+1,n)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%rlist(com_nod2D%rptr(com_nod2D%rPEnum+1)-1)) 
  do np=1,com_nod2D%rPEnum
     prank=com_nod2D%rPE(np)
     do n=1, nod2D
        if((aux(mype+1,n)==1).and.(part2D(n)==prank)) then
           counter=counter+1
           com_nod2D%rlist(counter)=n
        end if
     end do
  end do

  ! Summary of this piece: mype receives
  ! information on external 2D nodes from
  ! comm_nod2D%rPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%rPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%rptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)
  ! Putting everything into structure takes many operations, but
  ! without the structure we will need to many names and arrays
  ! Do not forget that we need also send part, and we need analogous
  ! information for 3D nodes.

  ! SENDING PART
  com_nod2D%sPEnum=0
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        com_nod2D%sPEnum=com_nod2D%sPEnum+1
     end if
  end do

  ! Their ranks (PE numbers)

  counter=0
  allocate(com_nod2D%sPE(com_nod2D%sPEnum))
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sPE(counter)=n-1
     end if
  end do

  ! Ptr to list of external nodes ordered by external PE ranks
  counter=0
  allocate(com_nod2D%sptr(com_nod2D%sPEnum+1)) 
  com_nod2D%sptr(1)=1
  do n=1, npes
     if(pnum(n,mype+1).ne.0) then
        counter=counter+1
        com_nod2D%sptr(counter+1)=com_nod2D%sptr(counter)+ pnum(n,mype+1)
     end if
  end do

  ! List itself

  counter=0
  allocate(com_nod2D%slist(com_nod2D%sptr(com_nod2D%sPEnum+1)-1)) 
  do np=1,com_nod2D%sPEnum
     prank=com_nod2D%sPE(np)
     do n=1, nod2D
        if((aux(prank+1,n)==1).and.(part2D(n)==mype)) then
           counter=counter+1
           com_nod2D%slist(counter)=n
        end if
     end do
  end do

  ! mype sends its data to
  ! comm_nod2D%sPEnum external PEs
  ! Their ranks (numbers) are in array
  ! comm_nod2D%sPE(:)
  ! Pointers to external node numbers are in
  ! comm_nod2D%sptr(:)
  ! The node numbers are in 
  ! comm_nod2D%list(:)

  ! 3D part looks similar and simply requres to use larger arrays
  ! to store indices:
  com_nod3D%rPEnum=com_nod2D%rPEnum
  com_nod3D%sPEnum=com_nod2D%sPEnum
  allocate(com_nod3D%sPE(com_nod2D%sPEnum))
  allocate(com_nod3D%rPE(com_nod2D%rPEnum))
  com_nod3D%rPE=com_nod2D%rPE
  com_nod3D%sPE=com_nod2D%sPE
  allocate(com_nod3D%rptr(com_nod2D%rPEnum+1))
  allocate(com_nod3D%sptr(com_nod2D%sPEnum+1))
  com_nod3D%rptr(1)=1
  com_nod3D%sptr(1)=1

  ! Compute pointers and lists of communicating 3D nodes

  do np=1, com_nod3D%rPEnum
     counter=0
     nini=com_nod2D%rptr(np)
     nend=com_nod2D%rptr(np+1)-1
     do n=nini,nend
        counter=counter+num_layers_below_nod2D(com_nod2D%rlist(n))+1
     end do
     com_nod3D%rptr(np+1)=com_nod3D%rptr(np)+counter
  end do
  allocate(com_nod3D%rlist(com_nod3D%rptr(com_nod3D%rPEnum+1)-1)) 
  counter=0
  do np=1, com_nod3D%rPEnum
     nini=com_nod2D%rptr(np)
     nend=com_nod2D%rptr(np+1)-1
     do n=nini,nend
        do nz=1,num_layers_below_nod2D(com_nod2D%rlist(n))+1
           counter=counter+1
           com_nod3D%rlist(counter)=nod3D_below_nod2D(nz,com_nod2D%rlist(n))
        end do
     end do
  end do

  do np=1, com_nod3D%sPEnum
     counter=0
     nini=com_nod2D%sptr(np)
     nend=com_nod2D%sptr(np+1)-1
     do n=nini,nend
        counter=counter+num_layers_below_nod2D(com_nod2D%slist(n))+1
     end do
     com_nod3D%sptr(np+1)=com_nod3D%sptr(np)+counter
  end do
  allocate(com_nod3D%slist(com_nod3D%sptr(com_nod3D%sPEnum+1)-1)) 
  counter=0
  do np=1, com_nod3D%sPEnum
     nini=com_nod2D%sptr(np)
     nend=com_nod2D%sptr(np+1)-1
     do n=nini,nend
        do nz=1,num_layers_below_nod2D(com_nod2D%slist(n))+1
           counter=counter+1
           com_nod3D%slist(counter)=nod3D_below_nod2D(nz,com_nod2D%slist(n))
        end do
     end do
  end do
  ! For the case of distributed memory, numbering of nodes should be 
  ! made consistent before this procedure. Then transition between
  ! global--local is elementary

  deallocate(pnum,aux)
end subroutine communication_nod
!
!==========================================================================
!
subroutine mymesh
  use o_MESH
  use o_ELEMENTS
  use g_PARFE 
  implicit none
  !
  integer     n, counter, q

  !======= NODES 

  ! 2D nodes
  ! Owned nodes + external nodes which I need:

  counter=0
  do n=1, nod2D
     if (part2D(n)==mype) counter=counter+1
  end do
  myDim_nod2D=counter
  eDim_nod2D=com_nod2D%rptr(com_nod2D%rPEnum+1)-1   
  allocate(myList_nod2D(myDim_nod2D+eDim_nod2D))
  counter=0   
  do n=1, nod2D
     if (part2D(n)==mype) then
        counter=counter+1
        myList_nod2D(counter)=n
     end if
  end do
  myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)=&
       com_nod2D%rlist

  ! Summary:  	     
  ! myList_nod2D(myDim_nod2D+1:myDim_nod2D+eDim_nod2D)
  ! contains external nodes which mype needs;    
  ! myList_nod2D(1:myDim_nod2D) contains owned nodes

  ! 3D nodes 
  counter=0

  do n=1, nod3D
     if (part3D(n)==mype) counter=counter+1
  end do
  myDim_nod3D=counter
  eDim_nod3D=com_nod3D%rptr(com_nod3D%rPEnum+1)-1   
  allocate(myList_nod3D(myDim_nod3D+eDim_nod3D))
  counter=0   
  do n=1, nod3D
     if (part3D(n)==mype) then
        counter=counter+1
        myList_nod3D(counter)=n
     end if
  end do
  myList_nod3D(myDim_nod3D+1:myDim_nod3D+eDim_nod3D)=&
       com_nod3D%rlist

  ! Summary:  	     
  ! myList_nod3D(myDim_nod3D+1:myDim_nod3D+eDim_nod3D)
  ! contains external nodes which mype needs;    
  ! myList_nod3D(1:myDim_nod3D) contains owned nodes

  !======= ELEMENTS

  ! 2D elements 
  counter=0
  do n=1, elem2D
     do q=1,3
        if(part2D(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem2D=counter
  allocate(myList_elem2D(myDim_elem2D))

  counter=0
  do n=1, elem2D
     do q=1,3
        if(part2D(elem2D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem2D(counter)=n
           exit
        end if
     end do
  end do

  ! 3D elements
  counter=0
  do n=1, elem3D
     do q=1,4
        if(part3D(elem3D_nodes(q,n))==mype) then
           counter=counter+1
           exit
        end if
     end do
  end do
  myDim_elem3D=counter
  allocate(myList_elem3D(myDim_elem3D))

  counter=0
  do n=1, elem3D
     do q=1,4
        if(part3D(elem3D_nodes(q,n))==mype) then
           counter=counter+1
           myList_elem3D(counter)=n
           exit
        end if
     end do
  end do

  ! Summary: element lists are only needed for parallel assembling.
  ! For this reason, we do not need to distinguish between owned
  ! and shared (but visited) elements --- they should be in the list.
end subroutine mymesh
subroutine com_2D(arr2d)
  use o_MESH
  use i_ARRAY
  use g_PARFE 
  implicit none

  integer       :: sreq(npes)
  integer       :: rreq(npes)
  integer       :: sstat(MPI_STATUS_SIZE,npes)
  integer       :: rstat(MPI_STATUS_SIZE,npes)
  integer       :: n, sn, rn, dest, nini, nend
  integer       :: offset, count, source
  real(kind=8)	:: arr2d(nod2d)

  ! Put data to be communicated into send buffer 
  sn=com_nod2D%sPEnum
  rn=com_nod2D%rPEnum

  do n=1, sn
     nini=com_nod2D%sptr(n)
     nend=com_nod2D%sptr(n+1) - 1
     s_buff_2d(n)%array=arr2d(com_nod2D%slist(nini:nend))
  end do

  do n=1, sn
     dest=com_nod2D%sPE(n)
     nini=com_nod2D%sptr(n)
     count=com_nod2D%sptr(n+1) - nini

     call MPI_ISEND(s_buff_2d(n)%array, count, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)

     source=com_nod2D%rPE(n)
     nini=com_nod2D%rptr(n)
     count=com_nod2D%rptr(n+1) - nini

     call MPI_IRECV(r_buff_2d(n)%array, count, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do

  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(sn,rreq,rstat, MPIerr)

  ! Put received data to their destination
  do n=1, rn
     nini=com_nod2D%rptr(n)
     nend=com_nod2D%rptr(n+1) - 1
     count=com_nod2D%rptr(n+1) - nini
     arr2d(com_nod2D%rlist(nini:nend))=r_buff_2d(n)%array
  end do
end subroutine com_2D
!
!===================================================================
!
subroutine com_3D(arr3d)
  use o_MESH
  use o_ARRAY
  use g_PARFE 
  implicit none
  !
  integer  	:: sreq(npes)
  integer  	:: rreq(npes)
  integer  	:: sstat(MPI_STATUS_SIZE,npes)
  integer  	:: rstat(MPI_STATUS_SIZE,npes)
  integer  	:: n, sn, rn, dest, nini, nend
  integer	:: offset, count, source
  real(kind=8)	:: arr3d(nod3d)


  ! Put data to be communicated into send buffer 
  sn=com_nod3D%sPEnum
  rn=com_nod3D%rPEnum

  do n=1, sn
     nini=com_nod3D%sptr(n)
     nend=com_nod3D%sptr(n+1) - 1
     count=com_nod3D%sptr(n+1) - nini
     s_buff_3d(n)%array=arr3d(com_nod3D%slist(nini:nend))
  end do

  do n=1, sn
     dest=com_nod3D%sPE(n)
     nini=com_nod3D%sptr(n)
     count=com_nod3D%sptr(n+1) - nini

     call MPI_ISEND(s_buff_3d(n)%array, count, MPI_DOUBLE_PRECISION, dest, mype, & 
          MPI_COMM_WORLD, sreq(n), MPIerr)

     source=com_nod3D%rPE(n)
     nini=com_nod3D%rptr(n)
     count=com_nod3D%rptr(n+1) - nini

     call MPI_IRECV(r_buff_3d(n)%array, count, MPI_DOUBLE_PRECISION, source, &
          source, MPI_COMM_WORLD, rreq(n), MPIerr) 
  end do

  call MPI_WAITALL(sn,sreq,sstat, MPIerr)
  call MPI_WAITALL(sn,rreq,rstat, MPIerr)

  ! Put received data to their destination
  do n=1, rn
     nini=com_nod3D%rptr(n)
     nend=com_nod3D%rptr(n+1) - 1
     count=com_nod3D%rptr(n+1) - nini
     arr3d(com_nod3D%rlist(nini:nend))=r_buff_3d(n)%array
  end do
end subroutine com_3D
!
!===================================================================
!
subroutine broadcast3D(arr3D, arr3Dglobal)
  ! Makes nodal information available to all PE
  ! arr3d is any array like TF or SF of local size.
  ! arr3Dglobal is an array of nod3D size which 
  ! should be allocated before calling this routine.
  ! It will be filled with information on other PE in 
  ! natural numbering. The routine can be used to organize
  ! output in the same way as in global memory setup  
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer :: ireals
  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr3D(myDim_nod3D+eDim_nod3D)
  real(kind=8) ::  arr3Dglobal(nod3D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)

  if ( mype == 0 ) then
     if (npes>1) then
        arr3Dglobal(myList_nod3D(1:myDim_nod3D))=arr3D(1:myDim_nod3D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_WORLD, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_WORLD, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_WORLD, status, MPIerr )

        do i = 1, nTS
           arr3Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod3D), isendbuf(1:myDim_nod3D) )
     do n = 1, myDim_nod3D
        isendbuf(n) = myList_nod3D(n)
        sendbuf(n)  = arr3D(n)
     enddo
     call MPI_SEND( myDim_nod3D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod3D, MPI_INTEGER, 0, 1, &
          MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod3D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_WORLD, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr3Dglobal, nod3d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, MPIerr)

end subroutine broadcast3D
!
!===================================================================
!
subroutine broadcast2D(arr2D, arr2Dglobal)
  ! Makes nodal information available to all PE 
  ! As the preceeding routine, but for 2D arrays
  !===================================================================

  use g_PARFE
  use o_MESH
  use o_ELEMENTS

  implicit none

  integer :: ireals
  integer      ::  i, n, nTS, sender, status(MPI_STATUS_SIZE)
  integer, allocatable, dimension(:) ::  isendbuf, irecvbuf

  real(kind=8) ::  arr2D(myDim_nod2D+eDim_nod2D)
  real(kind=8) ::  arr2Dglobal(nod2D)
  real(kind=8), allocatable, dimension(:) ::  sendbuf, recvbuf

  call MPI_Barrier(MPI_COMM_WORLD, MPIERR)
  if ( mype == 0 ) then
     if (npes>1) then
        arr2Dglobal(myList_nod2D(1:myDim_nod2D))=arr2D(1:myDim_nod2D)
     end if
     do  n = 1, npes-1

        call MPI_RECV( nTS, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
             0, MPI_COMM_WORLD, status, MPIerr )
        sender = status(MPI_SOURCE)
        allocate( recvbuf(1:nTS), irecvbuf(1:nTS) )
        call MPI_RECV( irecvbuf(1), nTS, MPI_INTEGER, sender, &
             1, MPI_COMM_WORLD, status, MPIerr )
        call MPI_RECV( recvbuf(1), nTS, MPI_DOUBLE_PRECISION, sender, &
             2, MPI_COMM_WORLD, status, MPIerr )

        do i = 1, nTS
           arr2Dglobal(irecvbuf(i)) = recvbuf(i)
        enddo
        deallocate( recvbuf, irecvbuf )

     enddo

  else

     allocate( sendbuf(1:myDim_nod2D), isendbuf(1:myDim_nod2D) )
     do n = 1, myDim_nod2D
        isendbuf(n) = myList_nod2D(n)
        sendbuf(n)  = arr2D(n)
     enddo
     call MPI_SEND( myDim_nod2D, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( isendbuf(1), myDim_nod2D, MPI_INTEGER, 0, 1, &
          MPI_COMM_WORLD, MPIerr )
     call MPI_SEND( sendbuf(1), myDim_nod2D, MPI_DOUBLE_PRECISION, &
          0, 2, MPI_COMM_WORLD, MPIerr )
     deallocate( sendbuf, isendbuf )

  endif

  call MPI_BCAST( arr2Dglobal, nod2d, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, MPIerr)

end subroutine broadcast2D
!===================================================================
subroutine forcing_array_setup
  !inializing forcing fields 
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  use g_config
  implicit none

  integer    :: n2

  n2=myDim_nod2D+eDim_nod2D  

  ! Allocate memory for atmospheric forcing 
  allocate(shortwave(n2), longwave(n2))
  allocate(prec_rain(n2), prec_snow(n2))
  allocate(u_wind(n2), v_wind(n2))
  allocate(Tair(n2), shum(n2))
  allocate(runoff(n2), evaporation(n2))
  shortwave=0.
  longwave=0.
  prec_rain=0.
  prec_snow=0.
  u_wind=0.
  v_wind=0.
  Tair=0.
  shum=0.
  runoff=0.

  if(use_landice_water) then
     allocate(runoff_landice(n2))
     runoff_landice=0.0
  end if
 
  ! shortwave penetration
#ifdef use_sw_pene
  allocate(chl(n2))
  allocate(sw_3d(myDim_nod3d+eDim_nod3D))
  chl=0.0
#endif

  !for ice diagnose
#ifdef use_ice
  allocate(thdgr(n2), thdgrsn(n2), flice(n2))
  allocate(olat_heat(n2), osen_heat(n2), olwout(n2))
  thdgr=0.
  thdgrsn=0.
  flice=0.
  olat_heat=0.
  osen_heat=0.
  olwout=0.
#endif 

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  allocate(Ce_atm_oce_arr(n2))
  allocate(Ch_atm_oce_arr(n2))
  Cd_atm_oce_arr=Cd_atm_oce
  Ce_atm_oce_arr=Ce_atm_oce 
  Ch_atm_oce_arr=Ch_atm_oce
#ifdef use_ice
  allocate(Cd_atm_ice_arr(n2)) 
  Cd_atm_ice_arr=Cd_atm_ice   
#endif

  if(mype==0) write(*,*) 'forcing arrays have been set up'   

end subroutine forcing_array_setup
!
!----------------------------------------------------------------------
!
subroutine forcing_array_setup_OnlyOcean
  !inializing forcing fields for an ocean-alone case
  !currently only wind is applied.
  use o_param
  use o_mesh
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer    :: n2

  n2=myDim_nod2D+eDim_nod2D  

  ! Allocate memory for atmospheric forcing 
  allocate(u_wind(n2), v_wind(n2))
  u_wind=0.
  v_wind=0.

  ! drag coefficient and transfer coefficients for latent and sensible heat
  allocate(Cd_atm_oce_arr(n2))      
  Cd_atm_oce_arr=Cd_atm_oce

  if(mype==0) write(*,*) 'forcing arrays (for an ocean-alone case) have been set up'   

end subroutine forcing_array_setup_OnlyOcean
subroutine init_atm_forcing_OnlyOcean
  ! initialize the atmospheric forcing data for the ocean-alone model
  ! assume forcing data on T62 NCEP/NCAR grid
  use o_PARAM
  use o_MESH
  use o_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, i, k, n2
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:) 

  n2=myDim_nod2D+eDim_nod2D       

  ! predefinition/correction
  ! for the CORE case:
  if(wind_data_source=='CORE1' .or. wind_data_source=='CORE2') wind_ttp_ind=1

  if(mype==0) write(*,*) 'Forcing data which are constant in time are initialized'

end subroutine init_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------ 
!
subroutine update_atm_forcing_OnlyOcean
  ! update atmospheric forcing data for ocean-alone cases
  use o_PARAM
  use o_MESH
  use o_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_parfe
  use g_clock
  use g_config
  implicit none

  integer		:: i
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy
  real(kind=8)          :: rhoair_local
  real              	:: t1, t2

  data rhoair_local /1.3/

  t1=MPI_Wtime()  

  ! first, read forcing data
  call read_new_atm_forcing_OnlyOcean

  ! second, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call cal_wind_drag_coeff
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
     call ncar_ocean_fluxes_mode
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if

  ! third, compute wind stress
  do i=1,myDim_nod2d+eDim_nod2d     
     aux=sqrt(u_wind(i)**2+v_wind(i)**2)*rhoair_local 
     stress_x(i) = Cd_atm_oce_arr(i)*aux*u_wind(i)
     stress_y(i) = Cd_atm_oce_arr(i)*aux*v_wind(i)
  end do

  t2=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if

end subroutine update_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------
!
subroutine read_new_atm_forcing_OnlyOcean
  ! read the second record of atmospheric forcing data for the ocean alone cases 
  ! assume forcing data on T62 NCEP/NCAR grid
  use o_PARAM
  use o_MESH
  use o_array
   use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, m, i, k, n2
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:)       

  n2=myDim_nod2D+eDim_nod2D 

  !==========================================================================
  ! wind u and v
                
  if(wind_data_source=='CORE2') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind) 
       
        ! 10-m wind m/s ----------------------------------------

        filevari='u_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

     end if

  elseif(wind_data_source=='CORE1') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        itime=forcing_rec(wind_ttp_ind)
        
        ! 10-m wind m/s ----------------------------------------

        filevari='u_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

     end if
  endif

  !==========================================================================
  ! climatology of ssT/ssS
  ! in this case we use WOA05 (NetCDF files)

  if(update_forcing_flag(sss_ttp_ind)==1) then

     itime=forcing_rec(sss_ttp_ind)
     
     filevari='t0112an1'
     file=trim(ClimateDataPath)//trim(filevari)//'.nc'
     vari='t0112an1'
     call read_surf_hydrography_NetCDF(file, vari, itime, Tsurf)

     filevari='s0112an1'
     file=trim(ClimateDataPath)//trim(filevari)//'.nc'
     vari='s0112an1'
     call read_surf_hydrography_NetCDF(file, vari, itime, Ssurf)
  end if

end subroutine read_new_atm_forcing_OnlyOcean
!
!------------------------------------------------------------------------------
!
subroutine init_atm_forcing
  ! read in forcing data that are constant in time
  ! the time varying forcing fields will be read in by read_new_atm_forcing
  ! assume atmosphere forcing data on T62 NCEP/NCAR grid
  use o_PARAM
  use o_MESH
  use o_array
  use i_therm_parms
  use i_array
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, i, k, pr_n2d
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:) 

  pr_n2d=myDim_nod2D+eDim_nod2D       

  ! predefinition/correction
  ! for the CORE case:
  if(wind_data_source=='CORE1' .or. wind_data_source=='CORE2') wind_ttp_ind=1
  if(rad_data_source=='CORE1' .or. rad_data_source=='CORE2') rad_ttp_ind=2
  if(precip_data_source=='CORE1' .or. precip_data_source=='CORE2') precip_ttp_ind=3
  if(runoff_data_source=='CORE1' .or. runoff_data_source=='CORE2') runoff_ttp_ind=0
  if(sss_data_source=='CORE1' .or. sss_data_source=='CORE2') sss_ttp_ind=4


  !==========================================================================
  ! runoff    

  if(runoff_data_source=='CORE1' .or. runoff_data_source=='CORE2' ) then

     ! runoff in CORE is constant in time

     ! Warning: For a global mesh, conservative scheme is to be updated!!

     file=trim(ForcingDataPath)//trim(runoff_data_source)//'/runoff.nc'
     vari='Foxx_o_roff'
     check_dummy=.false.

     itime=1
     call read_other_NetCDF(file, vari, itime, runoff, check_dummy) 
     runoff=runoff/1000.0  ! Kg/s/m2 --> m/s
  end if
  

  !==========================================================================
  ! sss restoring

  if(restore_s_surf>0.) then
     if(sss_data_source=='AAOMIP') then

        ! taking the annual mean of PHC2 SSS

        file=trim(ForcingDataPath)//'CORE2'//'/PHC2_salx.nc'
        vari='SALT'
        check_dummy=.true.

        Ssurf=0.0

        allocate(aux(pr_n2d))
        do itime=1,12
           call read_other_NetCDF(file, vari, itime, aux, check_dummy) 
           Ssurf=Ssurf+aux
        end do
        Ssurf=Ssurf/12.0
        deallocate(aux)
     endif
  end if

  if(mype==0) write(*,*) 'Parts of forcing data (only constant in time fields) are read'
 
end subroutine init_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine update_atm_forcing
  ! update atmospheric forcing data
  use o_PARAM
  use o_MESH
  use o_array
  use i_array
  use i_dyn_parms
  use i_therm_parms
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_parfe
  use g_clock
  use g_config
  implicit none

  integer		:: i
  real(kind=8)		:: i_coef, aux
  real(kind=8)		:: dux, dvy
  real              	:: t1, t2

  t1=MPI_Wtime()  

  ! first, read forcing data
  call read_new_atm_forcing

  ! second, compute exchange coefficients
  ! 1) drag coefficient 
  if(AOMIP_drag_coeff) then
     call cal_wind_drag_coeff
  end if
  ! 2) drag coeff. and heat exchange coeff. over ocean in case using ncar formulae
  if(ncar_bulk_formulae) then
     call ncar_ocean_fluxes_mode
  elseif(AOMIP_drag_coeff) then
     cd_atm_oce_arr=cd_atm_ice_arr
  end if

  ! third, compute wind stress
  do i=1,myDim_nod2d+eDim_nod2d     
     dux=u_wind(i)-u_w(i) 
     dvy=v_wind(i)-v_w(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmoce_x(i) = Cd_atm_oce_arr(i)*aux*dux
     stress_atmoce_y(i) = Cd_atm_oce_arr(i)*aux*dvy
     dux=u_wind(i)-u_ice(i) 
     dvy=v_wind(i)-v_ice(i)
     aux=sqrt(dux**2+dvy**2)*rhoair
     stress_atmice_x(i) = Cd_atm_ice_arr(i)*aux*dux
     stress_atmice_y(i) = Cd_atm_ice_arr(i)*aux*dvy
  end do

  ! heat and fresh water fluxes are treated in i_therm and ice2ocean

  t2=MPI_Wtime()

  if (mod(istep,logfile_outfreq)==0 .and. mype==0) then
     write(*,*) 'update forcing data took', t2-t1
  end if

end subroutine update_atm_forcing
!
!------------------------------------------------------------------------------------
!
subroutine read_new_atm_forcing
  ! read the second record of atmospheric forcing data 
  ! assume forcing data on T62 NCEP/NCAR grid
  use o_PARAM
  use o_MESH
  use o_array
  use i_therm_parms
  use g_forcing_param
  use g_forcing_arrays
  use g_forcing_index
  use g_forcing_interp
  use g_read_CORE_NetCDF
  use g_read_NCEP_NetCDF
  use g_read_other_NetCDF
  use g_clock
  use g_parfe
  use g_config
  implicit none
  !
  integer, parameter        		:: nci=192, ncj=94 ! T62 grid
  integer                   		:: itime, m, i, k, n2
  integer                               :: readtype
  character(80)             		:: file
  character(15)             		:: vari, filevari
  character(4)				:: fileyear
  real(kind=8), dimension(nci,ncj)	:: array_nc, array_nc2
  real(kind=8), dimension(nod2D)    	:: array_fe
  logical                               :: check_dummy
  real(kind=8), allocatable             :: aux(:)       

  n2=myDim_nod2D+eDim_nod2D  

  !==========================================================================
  ! wind u and v, Tair, and shum              
  if(wind_data_source=='NCEP') then

     if(update_forcing_flag(wind_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind)

        ! three temporal types (6 hourly, daily and monthly) are possible 

        if(wind_ttp_ind==1) then ! 6 hourly data

           ! 10-m wind m/s ----------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'uwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'vwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'air.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)   
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'shum.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2)    

        elseif(wind_ttp_ind==2) then ! daily data      

           ! 10-m wind --------------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'uwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'vwnd.10m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)  
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'air.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)   
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'shum.2m.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

        elseif(wind_ttp_ind==3) then ! monthly data

           ! 10-m wind m/s ----------------------------------------

           vari='uwnd'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'uwnd10m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)

           vari='vwnd'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'vwnd10m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc2)
           call upside_down(array_nc2,nci,ncj)

           ! rotate wind
           if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

           ! interp wind to model grid
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
           call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

           ! 2-m temperature --------------------------------------

           vari='air'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'air2m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2)  
           Tair=Tair-tmelt  ! Kelvin --> degree Celcius

           ! 2 m specific humdity  Kg/Kg -------------------------

           vari='shum'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'shum2m.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

        end if

     end if

  elseif(wind_data_source=='CORE2') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(wind_ttp_ind) 

        ! 10-m wind m/s ----------------------------------------

        filevari='u_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2) 
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 
 !       u_wind=-4.0/1.4142
 !       v_wind=-4.0/1.4142
 
        u_wind=0.0
        v_wind=0.0
	
        ! 10-m temperature -------------------------------------

        filevari='t_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='T_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2) 
        Tair=Tair-tmelt  ! Kelvin --> degree celcium

        ! 10 m specific humdity  Kg/Kg -------------------------

        filevari='q_10.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='Q_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2) 

     end if

  elseif(wind_data_source=='CORE1') then

     ! in CORE 6-hourly wind is used 

     if(update_forcing_flag(wind_ttp_ind)==1) then

        itime=forcing_rec(wind_ttp_ind)

        ! 10-m wind m/s ----------------------------------------

        filevari='u_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='U_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)

        filevari='v_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='V_10_MOD'   
        call read_CORE_NetCDF(file, vari, itime, array_nc2)

        ! rotate wind
        if(rotated_grid) call rotate_T62_wind(array_nc, array_nc2)

        ! interp wind to model grid
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,u_wind,n2)   
        call forcing_linear_ip(array_nc2,nci,ncj,lint_ind,lint_weight,v_wind,n2) 

        ! 10-m temperature -------------------------------------

        filevari='t_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='T_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,Tair,n2) 
        Tair=Tair-tmelt  ! Kelvin --> Degree Celcius

        ! 10 m specific humdity  Kg/Kg -------------------------

        filevari='q_10'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='Q_10_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shum,n2)  

     end if

  endif


  !==========================================================================
  ! radiation 

  if(rad_data_source=='NCEP') then

     if(update_forcing_flag(rad_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)

        ! two temporal types (6 hourly, daily) are possible 

        if(rad_ttp_ind==1) then ! 6 hourly data

           ! short wave W/m2 --------------------------------------

           vari='dswrf'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'dswrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

           ! long wave W/m2 ---------------------------------------

           vari='dlwrf'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'dlwrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2)  

        elseif(rad_ttp_ind==2) then ! daily data

           ! short wave W/m2 --------------------------------------

           vari='dswrf'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'dswrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

           ! long wave W/m2 ---------------------------------------

           vari='dlwrf'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'dlwrf.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2) 

        end if

     end if

  elseif(rad_data_source=='CORE2') then

     ! in CORE daily radiation fluxes are used 

     if(update_forcing_flag(rad_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(rad_ttp_ind)

        ! short wave W/m2 --------------------------------------

        filevari='ncar_rad.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='SWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

        ! long wave W/m2 ---------------------------------------

        vari='LWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2) 

     end if

  elseif(rad_data_source=='CORE1') then

     ! in CORE daily radiation fluxes are used 

     if(update_forcing_flag(rad_ttp_ind)==1) then

        itime=forcing_rec(rad_ttp_ind)

        ! short wave W/m2 --------------------------------------

        filevari='ncar_rad'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='SWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,shortwave,n2) 

        ! long wave W/m2 ---------------------------------------

        vari='LWDN_MOD'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,longwave,n2)  

     end if

  end if


  !==========================================================================
  ! precipitation

  if(precip_data_source=='NCEP') then

     if(update_forcing_flag(precip_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)

        ! four temporal types (6 hourly, daily and monthly, monthly ltm) are possible 

        if(precip_ttp_ind==1) then ! 6 hourly data

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_6hourly/'//'prate.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)  
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==2) then ! daily data      

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_daily/'//'prate.sfc.gauss.'//fileyear//'.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==3) then ! monthly data 

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_monthly/'//'prate.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        elseif(precip_ttp_ind==4) then ! monthly ltm data 

           ! total precip mm/s ------------------------------------

           vari='prate'
           file=trim(ForcingDataPath)//'NCEP_mon_ltm/'//'prate.mon.mean.nc'
           call read_NCEP_NetCDF(file, vari, itime, array_nc)
           call upside_down(array_nc,nci,ncj)
           call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
           prec_rain=prec_rain/1000.  ! mm/s --> m/s

        end if

     end if

  elseif(precip_data_source=='CORE2') then

     ! in CORE monthly precipitation is used; 
     ! And rain and snow are separated.

     if(update_forcing_flag(precip_ttp_ind)==1) then

        fileyear=cyearnew
        itime=forcing_rec(precip_ttp_ind)

        ! rain mm/s --------------------------------------------

        filevari='ncar_precip.'
        file=trim(ForcingDataPath)//'CORE2/'//trim(filevari)//fileyear//'.nc'
        vari='RAIN'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2) 
        prec_rain=prec_rain/1000.  ! mm/s --> m/s

        ! snow mm/s --------------------------------------------

        vari='SNOW'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_snow,n2)  
        prec_snow=prec_snow/1000.  ! mm/s --> m/s

     end if

  elseif(precip_data_source=='CORE1') then

     ! in CORE monthly precipitation is used; 
     ! And rain and snow are separated.

     if(update_forcing_flag(precip_ttp_ind)==1) then

        itime=forcing_rec(precip_ttp_ind)

        ! rain mm/s --------------------------------------------

        filevari='ncar_precip'
        file=trim(ForcingDataPath)//'CORE1/'//trim(filevari)//'.nc'
        vari='RAIN'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_rain,n2)  
        prec_rain=prec_rain/1000.  ! mm/s --> m/s

        ! snow mm/s --------------------------------------------

        vari='SNOW'
        call read_CORE_NetCDF(file, vari, itime, array_nc)
        call forcing_linear_ip(array_nc,nci,ncj,lint_ind,lint_weight,prec_snow,n2)  
        prec_snow=prec_snow/1000.  ! mm/s --> m/s

     end if

  end if


  !==========================================================================
  ! runoff  

  if(runoff_data_source=='Dai09') then

     if(update_forcing_flag(runoff_ttp_ind)==1) then
        if(runoff_ttp_ind==4) then
           !climatology monthly mean

           itime=forcing_rec(runoff_ttp_ind)
           file=trim(MeshPath)//'runoff_on_grid//runoff_clim.nc' 
           vari='runoff'

           call read_2ddata_on_grid_NetCDF(file,vari,itime,runoff)

           !kg/m2/s -> m/s
           runoff=runoff/1000.

        elseif(runoff_ttp_ind==3) then
           !monthly data

           write(*,*) 'Monthly runoff need to be updated. Forced to stop.'
           call par_ex
           stop
        end if
     end if

  elseif(runoff_data_source=='AAOMIP') then

     ! runoff is monthly ltm in AOMIP/AAOMIP

     if(update_forcing_flag(runoff_ttp_ind)==1) then

        allocate(aux(nod2D))

        itime=forcing_rec(runoff_ttp_ind)

        file=trim(ForcingDataPath)//'AAOMIP'//'/river_runoff.dat'
        if(system==1) then
           readtype=2
        else
           readtype=8
        end if

        open(101,file=trim(file),form='unformatted', access='direct',recl=readtype*nod2d)
        read(101,rec=itime) aux                
        runoff=aux(myList_nod2D)        
        close(101)

        deallocate(aux)
     end if

  end if


  !==========================================================================
  ! sss restoring

  if(restore_s_surf>0.) then
     if(sss_data_source=='CORE1' .or. sss_data_source=='CORE2') then

        ! sss is monthly ltm in CORE cases

        if(update_forcing_flag(sss_ttp_ind)==1) then

           itime=forcing_rec(sss_ttp_ind)

           file=trim(ForcingDataPath)//trim(sss_data_source)//'/PHC2_salx.nc'
           vari='SALT'
           check_dummy=.true.
           call read_other_NetCDF(file, vari, itime, Ssurf, check_dummy)  

        end if

     end if
  end if


  !==========================================================================
  ! chlorophyll climatology

#ifdef use_sw_pene
  ! currently only one type supported
  ! perpetual monthly chl climatology prepared by Sweeney et al. 2005

  if(update_forcing_flag(4)==1) then

     allocate(aux(nod2D))

     itime=forcing_rec(4)

     if(system==1) then
        readtype=2
     else
        readtype=8
     end if

     file=trim(ForcingDataPath)//'chlorophyll'//'/Chl_Sweeney.dat'
     open(101,file=trim(file),form='unformatted', access='direct',recl=readtype*nod2d)
     read(101,rec=itime) aux       
     chl=aux(myList_nod2D)    
     close(101)

     deallocate(aux)
  end if

#endif

end subroutine read_new_atm_forcing
!
!---------------------------------------------------------------------------------------------------
!
subroutine rotate_T62_wind(xarray, yarray)
  ! rotate wind on T62 grid from geographical coord. to rotated coordinates.
  use o_param
  use g_config
  use g_rotate_grid
  implicit none

  integer, parameter 	:: ni=192, nj=94  ! NCEP and CORE are on the same grid.
  integer               :: i, j
  real(kind=8)      	:: cx(ni), cy(nj), xarray(ni,nj), yarray(ni,nj) 

  ! NCEP/CORE latitude
  cy=(/-88.542, -86.6531, -84.7532, -82.8508, -80.9473, -79.0435, &  
       -77.1394, -75.2351, -73.3307, -71.4262, -69.5217, -67.6171, &  
       -65.7125, -63.8079, -61.9033, -59.9986, -58.0939, -56.1893, &  
       -54.2846, -52.3799, -50.4752, -48.5705, -46.6658, -44.7611,&  
       -42.8564, -40.9517, -39.0470, -37.1422, -35.2375, -33.3328, &  
       -31.4281, -29.5234, -27.6186, -25.7139, -23.8092, -21.9044, &  
       -19.9997, -18.0950, -16.1902, -14.2855, -12.3808, -10.47604, &  
       -8.57131, -6.66657, -4.76184, -2.8571, -0.952368, 0.952368, &  
       2.8571, 4.76184, 6.66657, 8.57131, 10.47604, 12.3808, &  
       14.2855, 16.1902, 18.095, 19.9997, 21.9044, 23.8092, &  
       25.7139, 27.6186, 29.5234, 31.4281, 33.3328, 35.2375,&  
       37.1422, 39.047,  40.9517, 42.8564, 44.7611, 46.6658,&  
       48.5705, 50.4752, 52.3799, 54.2846, 56.1893, 58.0939,&  
       59.9986, 61.9033, 63.8079, 65.7125, 67.6171, 69.5217, &  
       71.4262, 73.3307, 75.2351, 77.1394, 79.0435, 80.9473, &  
       82.8508, 84.7532, 86.6531, 88.542 /)*rad

  ! NCEP/CORE longitude
  cx(1)=0.0
  do i=2,ni
     cx(i)=cx(i-1)+1.875*rad
  enddo

  !rotate wind
  !cx cy are in radian
  do i=1,ni
     do j=1,nj
        call vector_g2r(xarray(i,j), yarray(i,j), cx(i), cy(j), 1)
     end do
  end do
end subroutine rotate_T62_wind
!
!---------------------------------------------------------------------------------------------------
subroutine ncar_ocean_fluxes_mode 
  ! Compute drag coefficient and the transfer coefficients for evaporation and sensible heat
  ! according to LY2004
  ! In this routine we assume air temperature and humidity are at the same height as wind speed.
  ! Otherwise, the code should be modified.
  ! There is a parameter z, which sets the height of wind speed. For the CORE data, z=10.0
  !
  ! Code from CORE website is adopted.
  ! original note:
  ! Over-ocean fluxes following Large and Yeager (used in NCAR models)           
  ! Coded by Mike Winton (Michael.Winton@noaa.gov) in 2004
  ! A bug was found by Laurent Brodeau (brodeau@gmail.com) in 2007.
  ! Stephen.Griffies@noaa.gov updated the code with the bug fix. 

  use o_mesh
  use i_therm_parms
  use i_array
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer, parameter :: n_itts = 2
  integer            :: i, j, m
  real :: cd_n10, ce_n10, ch_n10, cd_n10_rt    ! neutral 10m drag coefficients
  real :: cd, ce, ch, cd_rt                    ! full drag coefficients @ z
  real :: zeta, x2, x, psi_m, psi_h, stab      ! stability parameters
  real :: t, ts, q, qs, u, u10, tv, xx, dux, dvy
  real :: tstar, qstar, ustar, bstar
  real, parameter :: grav = 9.80, vonkarm = 0.40
  real, parameter :: q1=640380., q2=-5107.4    ! for saturated surface specific humidity
  real, parameter :: z = 10.0

  do i=1,myDim_nod2d+eDim_nod2d       
 
     t=tair(i) + tmelt					      ! degree celcium to Kelvin
     ts=t_oc_array(i) + tmelt				      !
     q=shum(i)
     qs=0.98*q1/rhoair*exp(q2/ts) 			      ! L-Y eqn. 5 
     tv = t*(1.0+0.608*q)
     dux=u_wind(i)-u_w(i)
     dvy=v_wind(i)-v_w(i)
     u = max(sqrt(dux**2+dvy**2), 0.5)           	      ! 0.5 m/s floor on wind (undocumented NCAR)
     u10 = u                                                  ! first guess 10m wind

     cd_n10 = (2.7/u10+0.142+0.0764*u10)/1.0e3                ! L-Y eqn. 6a
     cd_n10_rt = sqrt(cd_n10) 
     ce_n10 = 34.6 *cd_n10_rt/1.0e3       		      ! L-Y eqn. 6b
     stab = 0.5 + sign(0.5,t-ts)
     ch_n10 = (18.0*stab+32.7*(1.0-stab))*cd_n10_rt/1.e3      ! L-Y eqn. 6c

     cd = cd_n10                                 	      ! first guess for exchange coeff's at z
     ch = ch_n10
     ce = ce_n10
     do j=1,n_itts                                            ! Monin-Obukhov iteration
        cd_rt = sqrt(cd)
        ustar    = cd_rt*u                                    ! L-Y eqn. 7a
        tstar    = (ch/cd_rt)*(t-ts)              	      ! L-Y eqn. 7b
        qstar    = (ce/cd_rt)*(q-qs)              	      ! L-Y eqn. 7c
        bstar    = grav*(tstar/tv+qstar/(q+1.0/0.608))
        zeta     = vonkarm*bstar*z/(ustar*ustar) 	      ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta )          ! undocumented NCAR
        x2 = sqrt(abs(1.-16.*zeta))                           ! L-Y eqn. 8b
        x2 = max(x2, 1.0)                                     ! undocumented NCAR
        x = sqrt(x2)

        if (zeta > 0.) then
           psi_m = -5.*zeta                                    ! L-Y eqn. 8c
           psi_h = -5.*zeta                                    ! L-Y eqn. 8c
        else
           psi_m = log((1.+2.*x+x2)*(1+x2)/8.)-2.*(atan(x)-atan(1.0))  ! L-Y eqn. 8d
           psi_h = 2.*log((1.+x2)/2.)                                  ! L-Y eqn. 8e
        end if

        u10 = u/(1.0+cd_n10_rt*(log(z/10.)-psi_m)/vonkarm)        ! L-Y eqn. 9 !why cd_n10_rt not cd_rt
        cd_n10 = (2.7/u10+0.142+0.0764*u10)/1.e3                  ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10) 
        ce_n10 = 34.6*cd_n10_rt/1.e3                              ! L-Y eqn. 6b again
        stab = 0.5 + sign(0.5,zeta)
        ch_n10 = (18.0*stab+32.7*(1.0-stab))*cd_n10_rt/1.e3       ! L-Y eqn. 6c again
        !z0 = 10*exp(-vonkarm/cd_n10_rt)                          ! diagnostic

        xx = (log(z/10.)-psi_m)/vonkarm
        cd = cd_n10/(1.0+cd_n10_rt*xx)**2             		  ! L-Y 10a
        xx = (log(z/10.)-psi_h)/vonkarm
        ch = ch_n10/(1.0+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10b (corrected code aug2007)
        ce = ce_n10/(1.0+ce_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10)     ! 10c (corrected code aug2007)
     end do

     cd_atm_oce_arr(i)=cd
     ch_atm_oce_arr(i)=ch
     ce_atm_oce_arr(i)=ce 
  end do

end subroutine ncar_ocean_fluxes_mode
!
!---------------------------------------------------------------------------------------------------
!
subroutine cal_wind_drag_coeff
  ! Compute wind-ice drag coefficient following AOMIP
  use o_mesh
  use i_array
  use g_forcing_arrays
  use g_parfe
  implicit none

  integer            :: i, m
  real(kind=8)       :: ws

  do i=1,myDim_nod2d+eDim_nod2d    
     ws=sqrt(u_wind(i)**2+v_wind(i)**2)
     cd_atm_ice_arr(i)=(1.1+0.04*ws)*1.0e-3
  end do

end subroutine cal_wind_drag_coeff
!
!---------------------------------------------------------------------------------------------------
subroutine oce_input
  ! read restart fields for ocean dynamics and active tracer variables
  use o_param
  use o_mesh
  use o_array
  use g_clock
  use g_config
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: ssh_varid, tra_varid(2)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: istart(2), icount(2), n3
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(nod2D), aux3(nod3D)) 
  n3=ToDim_nod3D           

  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.oce.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  status=nf_inq_varid(ncid, 'ssh', ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'u', u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'v', v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#ifdef use_non_hydrostatic
  status=nf_inq_varid(ncid, 'w', w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#else
  status=nf_inq_varid(ncid, 'wpot', wpot_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  status=nf_inq_varid(ncid, 'temp', tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'salt', tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  ! 2d fields
  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, ssh_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  ssh=aux2(myList_nod2D)         

  ! 3d fields
  istart=(/1,nrec/)
  icount=(/nod3d, 1/)

  status=nf_get_vara_double(ncid, u_varid, istart, icount, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1:n3)=aux3(myList_nod3D)     

  status=nf_get_vara_double(ncid, v_varid, istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1+n3:2*n3)=aux3(myList_nod3D)  

#ifdef use_non_hydrostatic
  status=nf_get_vara_double(ncid, w_varid, istart, icount, aux3) 
  uf(1+2*n3:3*n3)=aux3(myList_nod3D)
#else
  status=nf_get_vara_double(ncid, wpot_varid, istart, icount, aux3)
  w=aux3(myList_nod3D)             
#endif
  if (status .ne. nf_noerr) call handle_err(status)

  do j=1,2
     status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,j)=aux3(myList_nod3D)  
  end do

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! the next record to be saved
  save_count=nrec+1

  deallocate(aux3, aux2)   

end subroutine oce_input
!
!-------------------------------------------------------------------------
!
subroutine age_tracer_input
  use o_param
  use o_mesh
  use o_array
  use o_age_tracer_mod
  use g_clock
  use g_config
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: tra_varid(num_age_tracer)
  integer                   :: istart(2), icount(2), n3
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux3(:) 

  allocate(aux3(nod3D)) 
  n3=ToDim_nod3D           

  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.oce.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  do j=1,num_age_tracer 
     write(trind,'(i1)') j
     status=nf_inq_varid(ncid, 'age'//trind, tra_varid(j))
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  ! 3d age tracer fields
  istart=(/1,nrec/)
  icount=(/nod3d, 1/)

  do j=1,num_age_tracer
     status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,index_age_tracer(j))=aux3(myList_nod3D)  
  end do

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux3)   

end subroutine age_tracer_input
!
!-------------------------------------------------------------------------
!
subroutine passive_tracer_input
  use o_param
  use o_mesh
  use o_array
  use o_age_tracer_mod
  use o_passive_tracer_mod
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: tra_varid(num_age_tracer)
  integer                   :: istart(2), icount(2), n3
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux3(:) 

  allocate(aux3(nod3D)) 
  n3=ToDim_nod3D           

  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.oce.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  do j=1,num_passive_tracer 
     write(trind,'(i1)') j
     status=nf_inq_varid(ncid, 'ptr'//trind, tra_varid(j))
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  ! 3d age tracer fields
  istart=(/1,nrec/)
  icount=(/nod3d, 1/)

  do j=1,num_passive_tracer
     status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,index_passive_tracer(j))=aux3(myList_nod3D)  
  end do

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux3)   

end subroutine passive_tracer_input
!
!-------------------------------------------------------------------------
!
subroutine ice_input
  use o_mesh
  use i_array
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, nrec
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: istart(2), icount(2)
  character(100)            :: filename
  real(kind=8), allocatable :: aux2(:)

  allocate(aux2(nod2D))  

  ! open files
  filename=trim(ResultPath)//runid//'.'//cyearold//'.ice.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  status=nf_inq_varid(ncid, 'area', area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hice', hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'uice', uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'vice', vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, area_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  a_ice=aux2(myList_nod2D)     
  status=nf_get_vara_double(ncid, hice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, hsnow_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_snow=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, uice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  u_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, vice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  v_ice=aux2(myList_nod2D)       
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux2)   
end subroutine ice_input
!
!-------------------------------------------------------------------------
!
subroutine read_prepared_initial_ice
  use o_mesh
  use i_array
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, nrec
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: istart(2), icount(2)
  character(100)            :: filename
  real(kind=8), allocatable :: aux2(:)

  allocate(aux2(nod2D))  

  ! open files
  filename=trim(ResultPath)//runid//'.'//'initial_ice.nc'
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) then
     print*,'ERROR: CANNOT READ initial ice FILE CORRECTLY !'
     print*,'Error in opening netcdf file'//filename
     call par_ex 
     stop
  endif
  
  ! inquire variable id
  status=nf_inq_varid(ncid, 'area', area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hice', hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'uice', uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'vice', vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! which record to read
  if(restartflag=='last') then
     status = nf_inq_dimid(ncid, 'T', dimid_rec)
     if(status .ne. nf_noerr) call handle_err(status)
     status = nf_inq_dimlen(ncid, dimid_rec, nrec)
     if(status .ne. nf_noerr) call handle_err(status)
  else
     read(restartflag,'(i4)') nrec
  end if

  istart=(/1,nrec/)
  icount=(/nod2d, 1/)
  status=nf_get_vara_double(ncid, area_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  a_ice=aux2(myList_nod2D)     
  status=nf_get_vara_double(ncid, hice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, hsnow_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  m_snow=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, uice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  u_ice=aux2(myList_nod2D)      
  status=nf_get_vara_double(ncid, vice_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  v_ice=aux2(myList_nod2D)       
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux2)   
end subroutine read_prepared_initial_ice
!
!-----------------------------------------------------------------------------
!
subroutine read_init_ts
  use o_mesh
  use o_param
  use o_array
  use g_config
  use g_rotate_grid
  use g_PARFE
  implicit none
  !
  integer                     :: i, j, n
  integer                     :: num_lat_reg, num_lon_reg, num_lay_reg
  real(kind=8)                :: pp, pr, tt, ss, lon, lat
  real(kind=8), external      :: theta
  real(kind=8), allocatable   :: lon_reg(:), lat_reg(:), lay_reg(:)
  real(kind=8), allocatable   :: raw_data(:,:,:)
  real(kind=8), allocatable   :: temp_x(:), temp_y(:)

  ! open global T/S data files
  open(19,file=trim(ClimateDataPath)//trim(OceClimaDataName), status='old')

  ! read reg. grid
  read(19,*) num_lon_reg, num_lat_reg, num_lay_reg
  allocate(lon_reg(num_lon_reg))
  allocate(lat_reg(num_lat_reg))
  allocate(lay_reg(num_lay_reg))
  read(19,*) lon_reg
  read(19,*) lat_reg
  read(19,*) lay_reg
  allocate(raw_data(num_lon_reg,num_lat_reg,num_lay_reg))

  ! model grid coordinates
  allocate(temp_x(myDim_nod3d+eDim_nod3D), temp_y(myDim_nod3d+eDim_nod3D)) 
  do n=1, myDim_nod3d+eDim_nod3D                    
     if(rotated_grid) then
        call r2g(lon, lat, coord_nod3d(1,n), coord_nod3d(2,n))
        temp_x(n)=lon/rad   ! change unit to degree
        temp_y(n)=lat/rad
     else
        temp_x(n)=coord_nod3d(1,n)/rad   
        temp_y(n)=coord_nod3d(2,n)/rad
     end if
     ! change lon range to [0 360]
     if(temp_x(n)<0.) temp_x(n)=temp_x(n) + 360.0  
  end do

  ! read raw data and do interpolation
  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, myDim_nod3D+eDim_nod3D, temp_x, temp_y, coord_nod3d(3,:), tracer(:,1))

  do i=1, num_lon_reg
     do j=1, num_lat_reg
        read(19, *) raw_data(i,j,1:num_lay_reg)         
     end do
  end do
  call interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, lon_reg, lat_reg, lay_reg, &
       raw_data, myDim_nod3d+eDim_nod3D, temp_x, temp_y, coord_nod3d(3,:), tracer(:,2))

  close(19) 

  ! Convert in situ temperature into potential temperature
  pr=0.0_8
  do i=1,myDim_nod3d+eDim_nod3D    
     tt=tracer(i,1)
     ss=tracer(i,2)
     pp=abs(coord_nod3D(3,i))
     tracer(i,1)=theta(ss, tt, pp, pr)
  end do

  deallocate(temp_y, temp_x, raw_data, lay_reg, lat_reg, lon_reg)
end subroutine read_init_ts
!
!-----------------------------------------------------------------------------
!
subroutine read_MY_vara
  use o_MESH
  use o_array
  use o_mixing_my2p5_mod
  use g_config
  use g_PARFE
  implicit none

  real(kind=8), allocatable :: aux3(:) 

  allocate(aux3(nod3D))   

  open(35,file=trim(ResultPath)//runid//'_MY_restart.out', status='old')

  read(35,*) aux3
  Kv(:,1)=aux3(myList_nod3d)
  read(35,*) aux3
  Av=aux3(myList_nod3d)
  read(35,*) aux3
  Kq=aux3(myList_nod3d)
  read(35,*) aux3
  q2=aux3(myList_nod3d)
  read(35,*) aux3
  q2b=aux3(myList_nod3d)
  read(35,*) aux3
  q2l=aux3(myList_nod3d)
  read(35,*) aux3
  q2lb=aux3(myList_nod3d)

  close(35)

  deallocate(aux3)
end subroutine Read_MY_vara
subroutine init_output
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  character(100)            :: longname
  character(100)            :: filename
  character(1)              :: trind

  if(yearnew==yearold) return
  if (mype/=0) return

  write(*,*) 'initialize new output files'

  ! first, snapshots

  ! ocean

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'ssh', NF_DOUBLE, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_DOUBLE, 2, dimids, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_DOUBLE, 2, dimids, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_DOUBLE, 2, dimids, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  status = nf_def_var(ncid, 'wpot', NF_DOUBLE, 2, dimids, wpot_varid)
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  status = nf_def_var(ncid, 'temp', NF_DOUBLE, 2, dimids, tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_DOUBLE, 2, dimids, tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'ptr'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_passive_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'age'//trind, NF_DOUBLE, 2, dimids, &
             tra_varid(index_age_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
  longname='vertical velocity potential'
  status = nf_put_att_text(ncid, wpot_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, wpot_varid, 'units', 5, 'm.m/s')
  if (status .ne. nf_noerr) call handle_err(status)
#endif
  longname='potential temperature'
  status = nf_put_att_text(ncid, tra_varid(1), 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(1), 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='salinity'
  status = nf_put_att_text(ncid, tra_varid(2), 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(2), 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        longname='passive tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        !status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
        !'units', 3, 'NaN')
        !if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        longname='age tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'units', 4, 'Year')
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! ice

#ifdef use_ice
  filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'area', NF_DOUBLE, 2, dimids, area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hice', NF_DOUBLE, 2, dimids, hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hsnow', NF_DOUBLE, 2, dimids, hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'uice', NF_DOUBLE, 2, dimids, uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vice', NF_DOUBLE, 2, dimids, vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='ice concentration [0 to 1]'
  status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective ice thickness'
  status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective snow thickness'
  status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

#endif


  ! second, mean fields
#if defined(allow_calcmeans) || defined(allow_diag)
  call init_output_mean
#endif


  ! initialize the counter for saving results
  save_count=1

end subroutine init_output
!
!----------------------------------------------------------------------------
!
subroutine init_output_mean
  use o_param
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, dimid_rec, j
  integer                   :: dimid_2d, dimid_3d, dimids(2)
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid
  integer                   :: utemp_varid, vtemp_varid
  integer                   :: usalt_varid, vsalt_varid
  integer                   :: mixlay_varid, Kv_varid
  integer                   :: sgs_u_varid, sgs_v_varid, sgs_ut_varid
  integer                   :: sgs_vt_varid, sgs_us_varid, sgs_vs_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: thdgr_varid, thdgrsn_varid
  integer                   :: uhice_varid, vhice_varid
  integer                   :: uhsnow_varid, vhsnow_varid
  integer                   :: flice_varid, tair_varid
  integer                   :: shum_varid, uwind_varid, vwind_varid
  integer                   :: rain_varid, snow_varid, runoff_varid
  integer                   :: evap_varid, lwrd_varid, swrd_varid
  integer                   :: qnet_varid, wnet_varid
  integer                   :: olat_varid, osen_varid, olwout_varid
  integer                   :: virtual_salt_varid, relax_salt_varid
  integer                   :: stress_x_varid, stress_y_varid
  character(100)            :: longname
  character(100)            :: filename
  character(1)              :: trind

  if (mype/=0) return

#ifdef allow_calcmeans

  ! ocean

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.mean.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'ssh', NF_FLOAT, 2, dimids, ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 3D fields
  dimids(1) = dimid_3d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'u', NF_FLOAT, 2, dimids, u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_FLOAT, 2, dimids, v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_FLOAT, 2, dimids, w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'temp', NF_FLOAT, 2, dimids, tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_FLOAT, 2, dimids, tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'ptr'//trind, NF_FLOAT, 2, dimids, &
             tra_varid(index_passive_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        status = nf_def_var(ncid, 'age'//trind, NF_FLOAT, 2, dimids, &
             tra_varid(index_age_tracer(j)))
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='mean sea surface elevation'
  status = nf_put_att_text(ncid, ssh_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, ssh_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean zonal velocity'
  status = nf_put_att_text(ncid, u_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, u_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean meridional velocity'
  status = nf_put_att_text(ncid, v_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, v_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean vertical velocity'
  status = nf_put_att_text(ncid, w_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, w_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean potential temperature'
  status = nf_put_att_text(ncid, tra_varid(1), 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(1), 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='mean salinity'
  status = nf_put_att_text(ncid, tra_varid(2), 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, tra_varid(2), 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(trind,'(i1)') j
        longname='passive tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        !status = nf_put_att_text(ncid, tra_varid(index_passive_tracer(j)), &
        !'units', 3, 'NaN')
        !if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(trind,'(i1)') j
        longname='age tracer '//trind
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), &
             'description', len_trim(longname), longname) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, tra_varid(index_age_tracer(j)), & 
             'units', 4, 'Year')
        if (status .ne. nf_noerr) call handle_err(status)
     end do
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! ice
#ifdef use_ice

  filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.mean.nc'

  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the time and iteration variables
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_2d
  dimids(2) = dimid_rec

  status = nf_def_var(ncid, 'area', NF_FLOAT, 2, dimids, area_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hice', NF_FLOAT, 2, dimids, hice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'hsnow', NF_FLOAT, 2, dimids, hsnow_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'uice', NF_FLOAT, 2, dimids, uice_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'vice', NF_FLOAT, 2, dimids, vice_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='iteration_count'
  status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)

  longname='ice concentration [0 to 1]'
  status = nf_PUT_ATT_TEXT(ncid, area_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective ice thickness'
  status = nf_PUT_ATT_TEXT(ncid, hice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hice_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='effective snow thickness'
  status = nf_PUT_ATT_TEXT(ncid, hsnow_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, hsnow_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='zonal velocity'
  status = nf_PUT_ATT_TEXT(ncid, uice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, uice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='meridional velocity'
  status = nf_PUT_ATT_TEXT(ncid, vice_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, vice_varid, 'units', 3, 'm/s')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

#endif
#endif

  ! diagnose
#ifdef allow_diag

  ! coean
  if(diag_oce) then

     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'nodes_3d', nod3d, dimid_3d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     if(diag_oce_mix_layer) then
        status = nf_def_var(ncid, 'mixlay', NF_FLOAT, 2, dimids, mixlay_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! Define the netCDF variables for 3D fields
     dimids(1) = dimid_3d
     dimids(2) = dimid_rec

     if(diag_oce_transp) then
        status = nf_def_var(ncid, 'utemp', NF_FLOAT, 2, dimids, utemp_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vtemp', NF_FLOAT, 2, dimids, vtemp_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'usalt', NF_FLOAT, 2, dimids, usalt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'vsalt', NF_FLOAT, 2, dimids, vsalt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(Redi_GM .and. diag_oce_GM_vel) then
        status = nf_def_var(ncid, 'sgs_u', NF_FLOAT, 2, dimids, sgs_u_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_v', NF_FLOAT, 2, dimids, sgs_v_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        status = nf_def_var(ncid, 'sgs_ut', NF_FLOAT, 2, dimids, sgs_ut_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vt', NF_FLOAT, 2, dimids, sgs_vt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_us', NF_FLOAT, 2, dimids, sgs_us_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_def_var(ncid, 'sgs_vs', NF_FLOAT, 2, dimids, sgs_vs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     status = nf_def_var(ncid, 'Kv', NF_FLOAT, 2, dimids, Kv_varid)
     if (status .ne. nf_noerr) call handle_err(status)


     ! Assign long_name and units attributes to variables.

     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     if(diag_oce_mix_layer) then
        longname='mixed layer thickness'
        status = nf_put_att_text(ncid, mixlay_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, mixlay_varid, 'units', 1, 'm')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_transp) then
        longname='zonal advective flux of temperature (zonal velocity * temperature)'
        status = nf_put_att_text(ncid, utemp_varid, 'description', len_trim(longname), trim(longname)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, utemp_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='meridional advective flux of temperature (meridional velocity * temperature)'
        status = nf_put_att_text(ncid, vtemp_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vtemp_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='zonal advective flux of salinity (zonal velocity * salinity)'
        status = nf_put_att_text(ncid, usalt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, usalt_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='meridional advective flux of salinity (meridional velocity * salinity)'
        status = nf_put_att_text(ncid, vsalt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, vsalt_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(Redi_GM .and. diag_oce_GM_vel) then
        longname='SGS (GM) zonal velocity integrated from bottom (k*S_x)'
        status = nf_put_att_text(ncid, sgs_u_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_u_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS (GM) meridional velocity integrated from bottom (k*S_y)'
        status = nf_put_att_text(ncid, sgs_v_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_v_varid, 'units', 5, 'm^2/s')
        if (status .ne. nf_noerr) call handle_err(status)
     endif

     if(diag_oce_SGS_transp) then
        longname='SGS zonal temperature flux'
        status = nf_put_att_text(ncid, sgs_ut_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_ut_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional temperature flux'
        status = nf_put_att_text(ncid, sgs_vt_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vt_varid, 'units', 8, 'degC m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS zonal salinity flux'
        status = nf_put_att_text(ncid, sgs_us_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_us_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
        longname='SGS meridional salinity flux'
        status = nf_put_att_text(ncid, sgs_vs_varid, 'description', len_trim(longname), trim(longname))
        if (status .ne. nf_noerr) call handle_err(status)
        status = nf_put_att_text(ncid, sgs_vs_varid, 'units', 7, 'psu m/s')
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     longname='Instantaneous vertical diffusivity'
     status = nf_put_att_text(ncid, Kv_varid, 'description', len_trim(longname), trim(longname))
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, Kv_varid, 'units', 5, 'm^2/s')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  endif  !ocean

  ! ice

#ifdef use_ice
  if(diag_ice) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     status = nf_def_var(ncid, 'thdgr', NF_FLOAT, 2, dimids, thdgr_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'thdgrsn', NF_FLOAT, 2, dimids, thdgrsn_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uhice', NF_FLOAT, 2, dimids, uhice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vhice', NF_FLOAT, 2, dimids, vhice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uhsnow', NF_FLOAT, 2, dimids, uhsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vhsnow', NF_FLOAT, 2, dimids, vhsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'flice', NF_FLOAT, 2, dimids, flice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Assign long_name and units attributes to variables.
     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     longname='thermodynamic growth rate of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, thdgr_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, thdgr_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='melting rate of snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, thdgrsn_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, thdgrsn_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal advective flux of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, uhice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uhice_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional advective flux of eff. ice thickness'
     status = nf_PUT_ATT_TEXT(ncid, vhice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vhice_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal advective flux of eff. snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, uhsnow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uhsnow_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional advective flux of eff. snow thickness'
     status = nf_PUT_ATT_TEXT(ncid, vhsnow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vhsnow_varid, 'units', 5, 'm.m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='rate of flooding snow to ice'
     status = nf_PUT_ATT_TEXT(ncid, flice_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, flice_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  endif
#endif


  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.forcing.diag.nc'

     ! create a file
     status = nf_create(filename, nf_clobber, ncid)
     if (status.ne.nf_noerr) call handle_err(status)

     ! Define the dimensions
     status = nf_def_dim(ncid, 'nodes_2d', nod2d, dimid_2d)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_dim(ncid, 'T', NF_UNLIMITED, dimid_rec)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the time and iteration variables
     status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, dimid_rec, time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'iter', NF_INT, 1, dimid_rec, iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Define the netCDF variables for 2D fields.
     ! In Fortran, the unlimited dimension must come
     ! last on the list of dimids.
     dimids(1) = dimid_2d
     dimids(2) = dimid_rec

     status = nf_def_var(ncid, 'tair', NF_FLOAT, 2, dimids, tair_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'shum', NF_FLOAT, 2, dimids, shum_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'uwind', NF_FLOAT, 2, dimids, uwind_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'vwind', NF_FLOAT, 2, dimids, vwind_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'rain', NF_FLOAT, 2, dimids, rain_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'snow', NF_FLOAT, 2, dimids, snow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'runoff', NF_FLOAT, 2, dimids, runoff_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'evap', NF_FLOAT, 2, dimids, evap_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'lwrd', NF_FLOAT, 2, dimids, lwrd_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'swrd', NF_FLOAT, 2, dimids, swrd_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_def_var(ncid, 'qnet', NF_FLOAT, 2, dimids, qnet_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'olat', NF_FLOAT, 2, dimids, olat_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'osen', NF_FLOAT, 2, dimids, osen_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'olwout', NF_FLOAT, 2, dimids, olwout_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'wnet', NF_FLOAT, 2, dimids, wnet_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'virtual_salt', NF_FLOAT, 2, dimids, virtual_salt_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'relax_salt', NF_FLOAT, 2, dimids, relax_salt_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'stress_x', NF_FLOAT, 2, dimids, stress_x_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_def_var(ncid, 'stress_y', NF_FLOAT, 2, dimids, stress_y_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! Assign long_name and units attributes to variables.
     longname='model time'
     status = nf_PUT_ATT_TEXT(ncid, time_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, time_varid, 'units', 1, 's')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='iteration_count'
     status = nf_PUT_ATT_TEXT(ncid, iter_varid, 'long_name', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)

     longname='air temperature'
     status = nf_PUT_ATT_TEXT(ncid, tair_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, tair_varid, 'units', 4, 'degC')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='air specific humidity'
     status = nf_PUT_ATT_TEXT(ncid, shum_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, shum_varid, 'units', 5, 'kg/kg')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='zonal wind speed'
     status = nf_PUT_ATT_TEXT(ncid, uwind_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, uwind_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='meridional wind speed'
     status = nf_PUT_ATT_TEXT(ncid, vwind_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, vwind_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='precipitation rain'
     status = nf_PUT_ATT_TEXT(ncid, rain_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, rain_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='precipitation snow (in m/s water)'
     status = nf_PUT_ATT_TEXT(ncid, snow_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, snow_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='runoff'
     status = nf_PUT_ATT_TEXT(ncid, runoff_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, runoff_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='evaporation'
     status = nf_PUT_ATT_TEXT(ncid, evap_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, evap_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='atmosphere longwave radiation'
     status = nf_PUT_ATT_TEXT(ncid, lwrd_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, lwrd_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='atmosphere shortwave radiation'
     status = nf_PUT_ATT_TEXT(ncid, swrd_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, swrd_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)

     longname='net heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, qnet_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, qnet_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='latent heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, olat_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, olat_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='sensible heat flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, osen_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, osen_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='longwave radiation from ocean, downward positve'
     status = nf_PUT_ATT_TEXT(ncid, olwout_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, olwout_varid, 'units', 5, 'W/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='net freshwater flux to ocean, downward positive'
     status = nf_PUT_ATT_TEXT(ncid, wnet_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, wnet_varid, 'units', 3, 'm/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='virtual salt flux to ocean, >0 increase salinity'
     status = nf_PUT_ATT_TEXT(ncid, virtual_salt_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, virtual_salt_varid, 'units', 7, 'psu m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface salinity relaxation, >0 increase salinity'
     status = nf_PUT_ATT_TEXT(ncid, relax_salt_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, relax_salt_varid, 'units', 7, 'psu m/s')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface zonal wind stress'
     status = nf_PUT_ATT_TEXT(ncid, stress_x_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, stress_x_varid, 'units', 5, 'N/m^2')
     if (status .ne. nf_noerr) call handle_err(status)
     longname='ocean surface meridional wind stress'
     status = nf_PUT_ATT_TEXT(ncid, stress_y_varid, 'description', len_trim(longname), trim(longname)) 
     if (status .ne. nf_noerr) call handle_err(status)
     status = nf_put_att_text(ncid, stress_y_varid, 'units', 5, 'N/m^2')
     if (status .ne. nf_noerr) call handle_err(status)

     status = nf_enddef(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)

  endif
#endif

#endif

end subroutine init_output_mean
!
!--------------------------------------------------------------------------------------------
!
subroutine write_snapshots
  use o_array
  use o_mesh
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use i_array
  use g_config
  use g_clock
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: start(2), count(2), n3
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(nod2D), aux3(nod3D)) 
  n3=myDim_nod3D+eDim_nod3D             

  if (mype==0) then 

     sec_in_year=dt*istep

     ! ocean

     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#ifndef use_non_hydrostatic
     status=nf_inq_varid(ncid, 'wpot', wpot_varid)
     if (status .ne. nf_noerr) call handle_err(status)
#endif
     status=nf_inq_varid(ncid, 'temp', tra_varid(1))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'salt', tra_varid(2))
     if (status .ne. nf_noerr) call handle_err(status)

     if(use_passive_tracer) then
        do j=1,num_passive_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'ptr'//trind, tra_varid(index_passive_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     if(use_age_tracer) then
        do j=1,num_age_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'age'//trind, tra_varid(index_age_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)
  end if    !! mype==0     

  ! 2d fields
  call broadcast2D(ssh,aux2)  
  if(mype==0) then            
     start=(/1,save_count/)
     count=(/nod2d, 1/)
     status=nf_put_vara_double(ncid, ssh_varid, start, count, aux2) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  ! 3d fields
  call broadcast3D(uf(1:n3), aux3)   
  if (mype==0) then                  
     start=(/1,save_count/)
     count=(/nod3d, 1/)
     status=nf_put_vara_double(ncid, u_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  call broadcast3D(uf(1+n3:2*n3),aux3)  
  if(mype==0) then                      
     status=nf_put_vara_double(ncid, v_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

#ifdef use_non_hydrostatic
  call broadcast3D(uf(1+2*n3:3*n3), aux3)  
  if(mype==0) then                       
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#else
  call broadcast3D(wrhs,aux3)             
  if(mype==0) then                        
     status=nf_put_vara_double(ncid, w_varid, start, count, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast3D(w,aux3)             
  if(mype==0) then                     
     status=nf_put_vara_double(ncid, wpot_varid, start, count, aux3)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif
  do j=1,num_tracer
     call broadcast3D(tracer(:,j),aux3)    
     if(mype==0) then                      
        status=nf_put_vara_double(ncid, tra_varid(j), start, count, aux3) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! ice

#ifdef use_ice
  if(mype==0) then                     
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hice', hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'uice', uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vice', vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if     !! mype=0                      
  call broadcast2D(a_ice,aux2)                
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, area_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, hice_varid, start, count, aux2)   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_snow,aux2)              
  if(mype==0) then                            
     status=nf_put_vara_double(ncid, hsnow_varid, start, count, aux2)  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(u_ice,aux2)              
  if(mype==0) then                           
     status=nf_put_vara_double(ncid, uice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(v_ice,aux2)              
  if(mype==0) then                          
     status=nf_put_vara_double(ncid, vice_varid, start, count, aux2)    
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

  deallocate(aux3, aux2)
end subroutine write_snapshots
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part1
  ! write mean arrays and diagnose variables
  ! SGS parameterizations are saved by write_means_part2
  use o_mesh
  use o_array
  use o_passive_tracer_mod
  use o_age_tracer_mod
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j
  integer                   :: time_varid, iter_varid
  integer                   :: ssh_varid, tra_varid(num_tracer)
  integer                   :: u_varid, v_varid, w_varid
  integer                   :: utemp_varid, vtemp_varid
  integer                   :: usalt_varid, vsalt_varid
  integer                   :: mixlay_varid, Kv_varid
  integer                   :: area_varid, hice_varid, hsnow_varid
  integer                   :: uice_varid, vice_varid
  integer                   :: thdgr_varid, thdgrsn_varid
  integer                   :: uhice_varid, vhice_varid
  integer                   :: uhsnow_varid, vhsnow_varid
  integer                   :: flice_varid, tair_varid
  integer                   :: shum_varid, uwind_varid, vwind_varid
  integer                   :: rain_varid, snow_varid, runoff_varid
  integer                   :: evap_varid, lwrd_varid, swrd_varid
  integer                   :: qnet_varid, wnet_varid
  integer                   :: olat_varid, osen_varid, olwout_varid
  integer                   :: virtual_salt_varid, relax_salt_varid
  integer                   :: stress_x_varid, stress_y_varid
  integer                   :: start(2), count(2), n3
  real(kind=8)              :: sec_in_year
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  n3=myDim_nod3D+eDim_nod3D         
  allocate(aux2(nod2D), aux3(nod3D))  

  sec_in_year=dt*istep

#ifdef allow_calcmeans 
  if (mype==0) then 

     ! ocean

     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.mean.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'ssh', ssh_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'u', u_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'v', v_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'w', w_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'temp', tra_varid(1))
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'salt', tra_varid(2))
     if (status .ne. nf_noerr) call handle_err(status)

     if(use_passive_tracer) then
        do j=1,num_passive_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'ptr'//trind, tra_varid(index_passive_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     if(use_age_tracer) then
        do j=1,num_age_tracer
           write(trind,'(i1)') j
           status = nf_inq_varid(ncid, 'age'//trind, tra_varid(index_age_tracer(j)))
           if (status .ne. nf_noerr) call handle_err(status)
        end do
     end if

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if    !! mype=0        

  call broadcast2D(sshmean, aux2)       
  if(mype==0) then
     status=nf_put_vara_real(ncid, ssh_varid, start, count, real(aux2,4)) 
     if (status .ne. nf_noerr) call handle_err(status)    
     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)
  end if
  call broadcast3D(ufmean(1:n3),aux3)       
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, u_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast3D(ufmean(1+n3:2*n3),aux3)   
  if(mype==0) then                     
     status=nf_put_vara_real(ncid, v_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#ifdef use_non_hydrostatic
  call broadcast3D(ufmean(1+2*n3:3*n3),aux3)   
  if(mype==0) then                      
     status=nf_put_vara_real(ncid, w_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#else
  call broadcast3D(wrhsmean, aux3)      
  if(mype==0) then
     status=nf_put_vara_real(ncid, w_varid, start, count, real(aux3,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

  do j=1,num_tracer
     call broadcast3D(tracermean(:,j),aux3) 
     if(mype==0) then                    
        status=nf_put_vara_real(ncid, tra_varid(j), start, count, real(aux3,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end do

  if(mype==0) then                       
     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if


  ! ice
#ifdef use_ice
  if(mype==0) then                      
     ! open files
     filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.mean.nc'
     status = nf_open(filename, nf_write, ncid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! inquire variable id
     status=nf_inq_varid(ncid, 'time', time_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'iter', iter_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'area', area_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hice', hice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'hsnow', hsnow_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'uice', uice_varid)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_inq_varid(ncid, 'vice', vice_varid)
     if (status .ne. nf_noerr) call handle_err(status)

     ! write variables

     ! time and iteration
     status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
     if (status .ne. nf_noerr) call handle_err(status)
     status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
     if (status .ne. nf_noerr) call handle_err(status)

     ! 2d fields
     start=(/1,save_count/)
     count=(/nod2d, 1/)
  end if   !! mype=0                        

  call broadcast2D(a_ice_mean, aux2)        
  if(mype==0) then                         
     status=nf_put_vara_real(ncid, area_varid, start, count, real(aux2,4)) 
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_ice_mean, aux2)         
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, hice_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(m_snow_mean, aux2)        
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, hsnow_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(u_ice_mean, aux2)        
  if(mype==0) then                           
     status=nf_put_vara_real(ncid, uice_varid, start, count, real(aux2,4))  
     if (status .ne. nf_noerr) call handle_err(status)
  end if
  call broadcast2D(v_ice_mean, aux2)       
  if(mype==0) then                          
     status=nf_put_vara_real(ncid, vice_varid, start, count, real(aux2,4))   
     if (status .ne. nf_noerr) call handle_err(status)

     status=nf_close(ncid)
     if (status .ne. nf_noerr) call handle_err(status)
  end if
#endif

#endif
  ! mean arrays have been saved


#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(mype==0) then                          
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        if(diag_oce_mix_layer) then
           status=nf_inq_varid(ncid, 'mixlay', mixlay_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        endif

        if(diag_oce_transp) then
           status=nf_inq_varid(ncid, 'utemp', utemp_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vtemp', vtemp_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'usalt', usalt_varid)
           if (status .ne. nf_noerr) call handle_err(status)
           status=nf_inq_varid(ncid, 'vsalt', vsalt_varid)
           if (status .ne. nf_noerr) call handle_err(status)
        endif

        status=nf_inq_varid(ncid, 'Kv', Kv_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if    !! mype=0                        

     if(diag_oce_mix_layer) then
        call broadcast2D(mixlay_dep_mean,aux2)    
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, mixlay_varid, start, count, real(aux2,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     endif

     ! 3d fields
     start=(/1,save_count/)
     count=(/nod3d, 1/)

     ! transport
     if(diag_oce_transp) then
        ! the fields which are read
        call broadcast3D(uTFmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, utemp_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vTFmean,aux3)          
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vtemp_varid, start, count,real(aux3,4)) 
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(uSFmean,aux3)           
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, usalt_varid, start, count,real(aux3,4))  
           if (status .ne. nf_noerr) call handle_err(status)
        end if
        call broadcast3D(vSFmean,aux3)            
        if(mype==0) then                          
           status=nf_put_vara_real(ncid, vsalt_varid, start, count,real(aux3,4)) 
           if (status .ne. nf_noerr) call handle_err(status)
        end if
     endif

     ! Kv
     call broadcast3D(Kv,aux3)                  
     if(mype==0) then                          
        status=nf_put_vara_real(ncid, Kv_varid, start, count, real(aux3,4))    
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif  ! diag ocean

  ! ice
#ifdef use_ice
  if(diag_ice) then
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.ice.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'thdgr', thdgr_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'thdgrsn', thdgrsn_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uhice', uhice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vhice', vhice_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uhsnow', uhsnow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vhsnow', vhsnow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'flice', flice_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if    !! mype=0                        
     call broadcast2D(thdgr_mean,aux2)          
     if(mype==0) then                          
        status=nf_put_vara_real(ncid, thdgr_varid, start, count, real(aux2,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(thdgrsn_mean,aux2)        
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, thdgrsn_varid, start, count, real(aux2,4))
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(uhice_mean,aux2)          
     if(mype==0) then                            
        status=nf_put_vara_real(ncid, uhice_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(vhice_mean,aux2)          
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, vhice_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(uhsnow_mean,aux2)          
     if(mype==0) then                            
        status=nf_put_vara_real(ncid, uhsnow_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(vhsnow_mean,aux2)         
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, vhsnow_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
     end if
     call broadcast2D(flice_mean,aux2)          
     if(mype==0) then                           
        status=nf_put_vara_real(ncid, flice_varid, start, count, real(aux2,4)) 
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif
#endif


  ! forcing 
#ifdef use_ice
  if(diag_forcing) then
     ! open files
     if(mype==0) then 
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.forcing.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! inquire variable id
        status=nf_inq_varid(ncid, 'time', time_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'iter', iter_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'tair', tair_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'shum', shum_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'uwind', uwind_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'vwind', vwind_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'rain', rain_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'snow', snow_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'runoff', runoff_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'evap', evap_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'lwrd', lwrd_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'swrd', swrd_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'qnet', qnet_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'olat', olat_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'osen', osen_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'olwout', olwout_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'wnet', wnet_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'virtual_salt', virtual_salt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'relax_salt', relax_salt_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'stress_x', stress_x_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_inq_varid(ncid, 'stress_y', stress_y_varid)
        if (status .ne. nf_noerr) call handle_err(status)

        ! write variables

        ! time and iteration
        status=nf_put_vara_double(ncid, time_varid, save_count, 1, sec_in_year)
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_put_vara_int(ncid, iter_varid, save_count, 1, istep)
        if (status .ne. nf_noerr) call handle_err(status)

        ! 2d fields
        start=(/1,save_count/)
        count=(/nod2d, 1/)
     end if     ! mype=0                            
     call broadcast2D(tair_mean,aux2)                
     if(mype==0) then                                
        status=nf_put_vara_real(ncid, tair_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                           
     call broadcast2D(shum_mean,aux2)               
     if(mype==0) then                               
        status=nf_put_vara_real(ncid, shum_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                            
     call broadcast2D(uwind_mean,aux2)               
     if(mype==0) then                               
        status=nf_put_vara_real(ncid, uwind_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(vwind_mean,aux2)                 
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, vwind_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(rain_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, rain_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(snow_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, snow_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(runoff_mean,aux2)                
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, runoff_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(evap_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, evap_varid, start, count, real(aux2,4))     
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(lwrd_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, lwrd_varid, start, count, real(aux2,4))     
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(swrd_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, swrd_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(qnet_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, qnet_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(olat_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, olat_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(osen_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, osen_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(olwout_mean,aux2)                
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, olwout_varid, start, count, real(aux2,4))  
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(wnet_mean,aux2)                  
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, wnet_varid, start, count, real(aux2,4))    
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(virtual_salt_mean,aux2)          
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, virtual_salt_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(relax_salt_mean,aux2)            
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, relax_salt_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(stress_x_mean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, stress_x_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
     end if     ! mype=0                               
     call broadcast2D(stress_y_mean,aux2)              
     if(mype==0) then                                    
        status=nf_put_vara_real(ncid, stress_y_varid, start, count, real(aux2,4))   
        if (status .ne. nf_noerr) call handle_err(status)
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  endif
#endif

#endif

  deallocate(aux3, aux2)
end subroutine write_means_part1
!
!--------------------------------------------------------------------------------------------
!
subroutine write_means_part2
  ! write ocean SGS parameterizations
  use o_param
  use o_mesh
  use o_elements
  use o_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe
  use g_clock
  implicit none

#include "netcdf.inc" 

  integer                   :: m, elem, elnodes(4)
  integer                   :: status, ncid, sgs_varid
  integer                   :: start(2), count(2)
  real(kind=8)              :: array_3d(nod3d)
  character(100)            :: filename

  ! prepare cluster volume
  ! use wrhs as a temporary array
  wrhs=0.0
  do elem=1,myDim_elem3d                                          
     elnodes=elem3d_nodes(:,elem)
     wrhs(elnodes)=wrhs(elnodes)+voltetra(elem)
  end do

  ! 3d fields
  start=(/1,save_count/)
  count=(/nod3d, 1/)

  if(Redi_GM .and. diag_oce_GM_vel) then
     ! processing
     call process_elem2node(1,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_u', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(2,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_v', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

  end if

  if(diag_oce_SGS_transp) then
     ! processing
     call process_elem2node(3,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_ut', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(4,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vt', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(5,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_us', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if

     ! processing
     call process_elem2node(6,array_3d)
     if(mype==0) then
        ! open files
        filename=trim(ResultPath)//runid//'.'//cyearnew//'.oce.diag.nc'
        status = nf_open(filename, nf_write, ncid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! inquire variable id
        status=nf_inq_varid(ncid, 'sgs_vs', sgs_varid)
        if (status .ne. nf_noerr) call handle_err(status)
        ! write variables
        status=nf_put_vara_real(ncid, sgs_varid, start, count,real(array_3d,4))
        if (status .ne. nf_noerr) call handle_err(status)
        ! close
        status=nf_close(ncid)
        if (status .ne. nf_noerr) call handle_err(status)
     end if
  end if

end subroutine write_means_part2
!
!--------------------------------------------------------------------------------------------
!
subroutine output(directionflag)
  use g_config
  use g_clock
  use g_diag
  use g_PARFE
  implicit none

  logical :: do_output=.false.
  integer :: directionflag

  !check whether we want to do output
  if (output_length_unit.eq.'y') then
     call annual_output(do_output)
  else if (output_length_unit.eq.'m') then 
     call monthly_output(do_output) 
  else if (output_length_unit.eq.'d') then
     call daily_output(do_output)  
  else if (output_length_unit.eq.'h') then
     call hourly_output(do_output) 
  else if (output_length_unit.eq.'s') then
     call step_output(do_output) 
  else
     write(*,*) 'You did not specify a supported outputflag.'
     write(*,*) 'The program will stop to give you opportunity to do it.'
     call par_ex
     stop
  endif

  if (directionflag.eq.1) do_output=.true.  

  if (.not.do_output) return

  ! write results
  if(mype==0) write(*,*)'Do output (netCDF) ...'
  call write_snapshots

  ! write mean fields
#if defined(allow_calcmeans) || defined(allow_diag)
  call compute_means
  call write_means_part1
#ifdef allow_diag
  if(diag_oce .and. (diag_oce_SGS_transp .or. diag_oce_GM_vel)) call write_means_part2
#endif
  call clean_meanarrays
#endif

  save_count=save_count+1

end subroutine output
!
!--------------------------------------------------------------------------------------------
!
subroutine annual_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if ((daynew == ndpyr) .and. (timenew==86400.)) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine annual_output
!
!--------------------------------------------------------------------------------------------
!
subroutine monthly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (day_in_month==num_day_in_month(fleapyear,month) .and. &
       timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  end if

end subroutine monthly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine daily_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(daynew,output_length)==0 .and. timenew==86400.) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine daily_output
!
!--------------------------------------------------------------------------------------------
!
subroutine hourly_output(do_output)
  !decides whether it's time to do output
  use g_clock
  implicit none

  logical :: do_output

  if (mod(timenew, 3600.*output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine hourly_output
!
!--------------------------------------------------------------------------------------------
!
subroutine step_output(do_output)
  !decides whether it's time to do output
  use g_config
  implicit none

  logical :: do_output

  if (mod(istep, output_length)==0) then
     do_output=.true.
  else
     do_output=.false.
  endif

end subroutine step_output
!
!--------------------------------------------------------------------------------------------
!
subroutine handle_err(errcode)
  use g_parfe
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  call par_ex
  stop
end subroutine handle_err
!
!--------------------------------------------------------------------------------------------
!
subroutine oce_out
  use o_param
  use o_MESH
  use o_array
  use g_PARFE
  use g_config
  implicit none
  !
  integer                    :: i, j, n3
  real(kind=8), allocatable  :: aux2(:), aux3(:)

  n3=myDim_nod3D+eDim_nod3D              
  allocate(aux2(nod2D), aux3(nod3D))     

  call broadcast2D(ssh,aux2)
  call broadcast3D(uf(1:n3),aux3)
  if (mype==0) then      
     write(*,*) 'writing ocean results (ASCII)'
     open(35,file='ssh.out')
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
  end if
  call broadcast3D(uf(1+n3:2*n3), aux3)   
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1f9.5)') aux3(i)
        !write(35,'(1e10.3)') aux3(i)
     end do
     !
     do i=1,nod2d
        write(35,'(1f9.5)') aux2(i)
        !write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  if(mype==0) open(36,file='TS.out')
  do j=1,num_tracer
     call broadcast3D(tracer(:,j),aux3)
     if(mype==0) then
        do i=1,nod3D
           write(36,'(1f9.5)') aux3(i)
           !write(36,'(1e12.4)') aux3(i)
        end do
     end if
  end do
  if(mype==0) close(36)

#ifndef use_non_hydrostatic
  call broadcast3D(wrhs,aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#else
  call broadcast3D(uf(2*n3+1:3*n3), aux3)
  if(mype==0) then
     open(37,file='vvel.out') 
     do i=1,nod3D
        write(37,'(1e10.3)') aux3(i)
     end do
     close(37)
  end if
#endif

  call broadcast3D(Kv, aux3)
  if(mype==0) then 
     open(38,file='mix_coeff.out')
     do i=1,nod3D
        write(38,'(1e10.3)') aux3(i)
     end do
     close(38)
  end if

  deallocate(aux3, aux2)

end subroutine oce_out
!
!--------------------------------------------------------------------------------------------
!
subroutine ice_out
  use o_MESH
  use i_array
  use g_parfe

  implicit none
  integer :: i
  real(kind=8), allocatable   :: aux2(:) 

  allocate(aux2(nod2D))         

  call broadcast2D(m_ice, aux2) 
  if (mype==0) then 
     write(*,*) 'writing ice results (ASCII)'
     open(35,file='m_ice.out') 
     do i=1,nod2D
        !write(35,'(1f9.5)')  aux2(i)
        write(35,'(1e10.3)') aux2(i)
     end do
     close(35)
  end if

  call broadcast2D(a_ice, aux2)
  if(mype==0) then
     open(36,file='a_ice.out') 
     do i=1,nod2D
        !write(36, '(1f9.5)') aux2(i)
        write(36,'(1e10.3)') aux2(i)
     end do
     close(36)
  end if

  call broadcast2D(m_snow, aux2)
  if(mype==0) then 
     open(37,file='m_snow.out') 
     do i=1,nod2D
        !write(37, '(1f9.5)') aux2(i)
        write(37,'(1e10.3)') aux2(i)
     end do
     close(37)
  end if

  call broadcast2D(u_ice, aux2)
  if(mype==0) then
     open(38,file='u_ice.out')
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
  end if
  call broadcast2D(v_ice, aux2)
  if(mype==0) then
     do i=1,nod2D
        write(38, '(1f9.5)') aux2(i)
     end do
     close(38)
  end if

  call broadcast2D(net_heat_flux, aux2)
  if(mype==0) then
     open(39,file='heat_water_flux.out')
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
  end if
  call broadcast2D(fresh_wa_flux,aux2)
  if(mype==0) then
     do i=1,nod2D
        write(39, *) aux2(i)
     end do
     close(39)
  end if

  deallocate(aux2)

end subroutine ice_out
!
!--------------------------------------------------------------------------------------------
!
subroutine save_MY_vara
  use o_MESH
  use o_array
  use o_mixing_my2p5_mod
  use g_PARFE
  use g_config
  implicit none

  integer :: i
  real(kind=8), allocatable  :: aux3(:)

  allocate(aux3(nod3D))  

  if (mype==0) then
     write(*,*) 'writing MY2.5 variables (ascII)'
     Write(*,*) 'Notice: Currently MY variables are only saved once at the end of a run.'
     ! Later we may update to save MY in Netcdf formate when required.

     open(35,file=trim(ResultPath)//runid//'_MY_restart.out')
  end if

  call broadcast3D(Kv,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Av,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(Kq,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2b,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2l,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  call broadcast3D(q2lb,aux3)
  if(mype==0) then
     do i=1,nod3D
        write(35,'(1e11.4)') aux3(i)
     end do
  end if

  if(mype==0) close(35)

  deallocate(aux3)
end subroutine Save_MY_vara
subroutine init_meanarrays
  ! allocates and initializes fields used for mean arrays
  use o_mesh
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  use g_parfe

  implicit none
  integer       :: n2, n3

  n2=myDim_nod2D+eDim_nod3D        
  n3=myDim_nod3D+eDim_nod3D       

  ! for mean fields

#ifdef allow_calcmeans

  ! ocean part                    
  allocate(sshmean(n2))
#ifndef use_non_hydrostatic
  allocate(ufmean(2*n3)) 
  allocate(wrhsmean(n3))
#else
  allocate(ufmean(3*n3)) 
#endif
  allocate(tracermean(n3,num_tracer))

  ! ice part
#ifdef use_ice
  allocate(m_ice_mean(n2), a_ice_mean(n2), m_snow_mean(n2))
  allocate(u_ice_mean(n2), v_ice_mean(n2))
#endif

#endif

  ! for diagnostics

#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(diag_oce_transp) then
        allocate(uTFmean(n3),vTFmean(n3))
        allocate(uSFmean(n3),vSFmean(n3))
     end if
     if(Redi_GM .and. diag_oce_GM_vel) then
        allocate(sgs_u(myDim_elem3d), sgs_v(myDim_elem3d))
     end if
     if(diag_oce_SGS_transp) then
        allocate(sgs_ut(myDim_elem3d), sgs_vt(myDim_elem3d))
        allocate(sgs_us(myDim_elem3d), sgs_vs(myDim_elem3d))
     end if
     if(diag_oce_mix_layer) then
        allocate(mixlay_dep_mean(n2))
     end if
  end if

#ifdef use_ice
  ! ice
  if(diag_ice) then
     allocate(thdgr_mean(n2), thdgrsn_mean(n2))
     allocate(uhice_mean(n2), vhice_mean(n2))
     allocate(uhsnow_mean(n2), vhsnow_mean(n2))
     allocate(flice_mean(n2))
  end if

  ! forcing
  if(diag_forcing) then
     allocate(tair_mean(n2), shum_mean(n2))
     allocate(uwind_mean(n2), vwind_mean(n2))
     allocate(rain_mean(n2), snow_mean(n2))
     allocate(runoff_mean(n2), evap_mean(n2))
     allocate(lwrd_mean(n2), swrd_mean(n2))
     allocate(qnet_mean(n2), wnet_mean(n2))
     allocate(olat_mean(n2), osen_mean(n2))
     allocate(olwout_mean(n2))
     allocate(virtual_salt_mean(n2), relax_salt_mean(n2))
     allocate(stress_x_mean(n2), stress_y_mean(n2))
  end if
#endif

#endif

  call clean_meanarrays

  if(mype==0) write(*,*) 'Mean arrays have been set up'

  return
end subroutine init_meanarrays
!=============================================================================!


!=============================================================================!
subroutine add2meanarrays
  ! adds values to the mean-arrays
  use o_mesh
  use o_param
  use o_array
  use i_array
  use g_config
  use g_diag
  use g_meanarrays
  use g_forcing_arrays
  use g_parfe
  implicit none
  !
  integer       :: m, row, row2, row3, j

  meancounter=meancounter+1

#ifdef allow_calcmeans
#ifndef use_non_hydrostatic
  call vvel_nodes 
#endif
#endif

  do row=1,myDim_nod3d                      

     row2=row+myDim_nod3d+eDim_nod3D       

#ifdef allow_calcmeans

     ! ocean
     ufmean(row)        = ufmean(row) + uf(row)
     ufmean(row2)       = ufmean(row2) + uf(row2)
#ifdef use_non_hydrostatic
     row3=row2+myDim_nod3D+eDim_nod3D       
     ufmean(row3)       = ufmean(row3) + uf(row3)
#else
     wrhsmean(row)      = wrhsmean(row) + wrhs(row)
#endif
     do j=1,num_tracer
        tracermean(row,j)  = tracermean(row,j) + tracer(row,j)
     end do
#endif

#ifdef allow_diag
     if(diag_oce .and. diag_oce_transp) then
        uTFmean(row)       = uTFmean(row) + uf(row)*tracer(row,1)
        vTFmean(row)       = vTFmean(row) + uf(row2)*tracer(row,1)
        uSFmean(row)       = uSFmean(row) + uf(row)*tracer(row,2)
        vSFmean(row)       = vSFmean(row) + uf(row2)*tracer(row,2)
     endif
#endif
  end do

  !-----------------------------------------------------------

  do row=1,myDim_nod2d  

     ! ocean
#ifdef allow_calcmeans 
     sshmean(row)       = sshmean(row) + ssh(row)
#endif

     ! ice
#ifdef use_ice
#ifdef allow_calcmeans
     a_ice_mean(row)    = a_ice_mean(row) + a_ice(row)
     m_ice_mean(row)    = m_ice_mean(row) + m_ice(row)     
     m_snow_mean(row)   = m_snow_mean(row) + m_snow(row)
     u_ice_mean(row)    = u_ice_mean(row)  + u_ice(row)
     v_ice_mean(row)    = v_ice_mean(row)  + v_ice(row)
#endif
#ifdef allow_diag
     if(diag_ice) then
        thdgr_mean(row)    = thdgr_mean(row)  + thdgr(row)
        thdgrsn_mean(row)  = thdgrsn_mean(row) + thdgrsn(row)
        uhice_mean(row)    = uhice_mean(row) + u_ice(row)*m_ice(row)
        vhice_mean(row)    = vhice_mean(row) + v_ice(row)*m_ice(row)
        uhsnow_mean(row)   = uhsnow_mean(row) + u_ice(row)*m_snow(row)
        vhsnow_mean(row)   = vhsnow_mean(row) + v_ice(row)*m_snow(row)
        flice_mean(row)    = flice_mean(row) + flice(row)
     endif
#endif
#endif

     ! forcing
#ifdef allow_diag
#ifdef use_ice
     if(diag_forcing) then
        tair_mean(row)     = tair_mean(row) + Tair(row)
        shum_mean(row)     = shum_mean(row) + shum(row)
        uwind_mean(row)    = uwind_mean(row) + u_wind(row)
        vwind_mean(row)    = vwind_mean(row) + v_wind(row)   
        rain_mean(row)     = rain_mean(row) + prec_rain(row)
        snow_mean(row)     = snow_mean(row) + prec_snow(row)
        runoff_mean(row)   = runoff_mean(row) + runoff(row)
        evap_mean(row)     = evap_mean(row) + evaporation(row)
        lwrd_mean(row)     = lwrd_mean(row) + longwave(row)
        swrd_mean(row)     = swrd_mean(row) + shortwave(row)
        qnet_mean(row)     = qnet_mean(row) + net_heat_flux(row)
        wnet_mean(row)     = wnet_mean(row) + fresh_wa_flux(row)
        olat_mean(row)     = olat_mean(row) + olat_heat(row)
        osen_mean(row)     = osen_mean(row) + osen_heat(row)
        olwout_mean(row)   = olwout_mean(row) + olwout(row)
        virtual_salt_mean(row) = virtual_salt_mean(row) + virtual_salt(row)
        relax_salt_mean(row)   = relax_salt_mean(row) + relax_salt(row)
        stress_x_mean(row) = stress_x_mean(row) + stress_x(row)
        stress_y_mean(row) = stress_y_mean(row) + stress_y(row)
     endif
#endif
#endif

  enddo

  ! other diagnostics

  ! mixed layer thickness
#ifdef allow_diag
  if(diag_oce .and. diag_oce_mix_layer) then
     call compute_mixlay
  end if
#endif

  ! diagnose for SGS parameterization is done in the ts_rhs routine

  return
end subroutine add2meanarrays
!=============================================================================!


!=============================================================================!
subroutine compute_means
  !computes the mean values for output
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  real(kind=8)   :: cnt

  cnt=float(meancounter)

#ifdef allow_calcmeans

  ! ocean
  sshmean               = sshmean            /cnt
  ufmean                = ufmean             /cnt
#ifndef use_non_hydrostatic
  wrhsmean              = wrhsmean           /cnt
#endif
  tracermean            = tracermean         /cnt

  ! ice
#ifdef use_ice
  a_ice_mean            = a_ice_mean         /cnt
  m_ice_mean            = m_ice_mean         /cnt
  m_snow_mean           = m_snow_mean        /cnt
  u_ice_mean            = u_ice_mean         /cnt
  v_ice_mean            = v_ice_mean         /cnt
#endif

#endif

  ! diagnose
#ifdef allow_diag

  ! ocean
  if(diag_oce) then
     if(diag_oce_transp) then
        uTFmean               = uTFmean            /cnt
        vTFmean               = vTFmean            /cnt
        uSFmean               = uSFmean            /cnt
        vSFmean               = vSFmean            /cnt
     endif
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u                 = sgs_u              /cnt
        sgs_v                 = sgs_v              /cnt
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut                = sgs_ut             /cnt
        sgs_vt                = sgs_vt             /cnt
        sgs_us                = sgs_us             /cnt
        sgs_vs                = sgs_vs             /cnt
     endif
     if(diag_oce_mix_layer) then
        mixlay_dep_mean       = mixlay_dep_mean    /cnt
     endif
  endif

  ! ice
#ifdef use_ice
  if(diag_ice) then
     thdgr_mean            = thdgr_mean         /cnt
     thdgrsn_mean          = thdgrsn_mean       /cnt
     uhice_mean            = uhice_mean         /cnt
     vhice_mean            = vhice_mean         /cnt
     uhsnow_mean           = uhsnow_mean        /cnt
     vhsnow_mean           = vhsnow_mean        /cnt
     flice_mean            = flice_mean         /cnt
  endif
#endif

  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     tair_mean             = tair_mean          /cnt
     shum_mean             = shum_mean          /cnt
     uwind_mean            = uwind_mean         /cnt
     vwind_mean            = vwind_mean         /cnt
     rain_mean             = rain_mean          /cnt
     snow_mean             = snow_mean          /cnt
     runoff_mean           = runoff_mean        /cnt
     evap_mean             = evap_mean          /cnt
     lwrd_mean             = lwrd_mean          /cnt
     swrd_mean             = swrd_mean          /cnt
     qnet_mean             = qnet_mean          /cnt
     wnet_mean             = wnet_mean          /cnt
     olat_mean             = olat_mean          /cnt
     osen_mean             = osen_mean          /cnt
     olwout_mean           = olwout_mean        /cnt
     virtual_salt_mean     = virtual_salt_mean  /cnt
     relax_salt_mean       = relax_salt_mean    /cnt
     stress_x_mean         = stress_x_mean      /cnt
     stress_y_mean         = stress_y_mean      /cnt 
  endif
#endif

#endif

end subroutine compute_means
!=============================================================================!


!=============================================================================!
subroutine clean_meanarrays
  ! puts zeros into the mean-arrays
  use o_param
  use g_config
  use g_diag
  use g_meanarrays
  implicit none

  meancounter=0

#ifdef allow_calcmeans

  ! coean
  sshmean=0.
  ufmean=0.
#ifndef use_non_hydrostatic
  wrhsmean=0.
#endif
  tracermean=0.

  ! ice
#ifdef use_ice
  a_ice_mean=0.
  m_ice_mean=0.
  m_snow_mean=0.
  u_ice_mean=0.
  v_ice_mean=0.
  thdgr_mean=0.
  thdgrsn_mean=0.
  uhice_mean=0.
  vhice_mean=0.
  uhsnow_mean=0.
  vhsnow_mean=0.
  flice_mean=0.
#endif

#endif

  ! diagnose
#ifdef allow_diag

  if(diag_oce) then
     ! ocean
     if(diag_oce_transp) then
        uTFmean=0.
        vTFmean=0.
        uSFmean=0.
        vSFmean=0.
     endif
     if(Redi_GM .and. diag_oce_GM_vel) then
        sgs_u=0.
        sgs_v=0.
     endif
     if(diag_oce_SGS_transp) then
        sgs_ut=0.
        sgs_vt=0.
        sgs_us=0.
        sgs_vs=0.
     endif
     if(diag_oce_mix_layer) then
        mixlay_dep_mean=0.
     endif
  endif

  ! ice
#ifdef use_ice
  if(diag_ice) then
     thdgr_mean=0.
     thdgrsn_mean=0.
     uhice_mean=0.
     vhice_mean=0.
     uhsnow_mean=0.
     vhsnow_mean=0.
     flice_mean=0.
  endif
#endif

  ! forcing
#ifdef use_ice
  if(diag_forcing) then
     tair_mean=0.
     shum_mean=0.
     uwind_mean=0.
     vwind_mean=0.
     rain_mean=0.
     snow_mean=0.
     runoff_mean=0.
     evap_mean=0.
     lwrd_mean=0.
     swrd_mean=0.
     qnet_mean=0.
     wnet_mean=0.
     olat_mean=0.
     osen_mean=0.
     olwout_mean=0.
     virtual_salt_mean=0.
     relax_salt_mean=0.
     stress_x_mean=0.
     stress_y_mean=0.
  endif
#endif

#endif

  return
end subroutine clean_meanarrays
!=============================================================================!
subroutine interp_2d_field_v2(num_lon_reg, num_lat_reg, lon_reg, lat_reg, data_reg, missvalue, &
     num_mod, lon_mod, lat_mod, data_mod)
  !-------------------------------------------------------------------------------------
  ! A second version of 2D interpolation.
  ! This routine does 2d interpolation from a regular grid to specified nodes
  ! on the surface grid. The regular grid is assumed to be global. 
  ! In this version the missing value will be checked. If a model grid point is outside
  ! the ocean, the nearest value will be assigned.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! data_reg(num_lon_reg, num_lat_reg) 	data on the regular grid
  ! missvalue                           missing value in the raw data
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i, ii, jj, k, nod_find
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer, intent(in)         	:: num_lon_reg, num_lat_reg, num_mod
  real(kind=8) 			:: x, y, diff, d, dmin
  real(kind=8)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=8)                  :: data(2,2)
  real(kind=8)                  :: data_lo, data_up
  real(kind=8), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=8), intent(in)	:: data_reg(num_lon_reg, num_lat_reg), missvalue
  real(kind=8), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod)
  real(kind=8), intent(out)  	:: data_mod(num_mod)
  !
  if(lon_reg(1)<0.0 .or. lon_reg(num_lon_reg)>360.) then
     write(*,*) 'Error in 2D interpolation!'
     write(*,*) 'The regular grid is not in the proper longitude range.'
     call par_ex
     stop
  end if

  do n=1,num_mod
     x=lon_mod(n)
     y=lat_mod(n)
     ! find the surrounding rectangular box and get interpolation ratios
     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0-rt_lat1
     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
        rt_lon1=1.0-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     end if
     !
     data(1,1)=data_reg(ind_lon_l,ind_lat_l)
     data(1,2)=data_reg(ind_lon_l,ind_lat_h)
     data(2,1)=data_reg(ind_lon_h,ind_lat_l)
     data(2,2)=data_reg(ind_lon_h,ind_lat_h)
     !
     ! interpolate data
     if(any(data==missvalue)) then
        dmin=10000.0
        nod_find=0
        do k=1,5
           do ii=max(1,ind_lon_l-k+1),min(num_lon_reg,ind_lon_l+k)
              do jj=max(1,ind_lat_l-k+1),min(num_lat_reg,ind_lat_l+k)
                 if(data_reg(ii,jj)==missvalue) cycle
                 d=(x-lon_reg(ii))**2 + (y-lat_reg(jj))**2
                 if(d<dmin) then
                    nod_find=1
                    data_mod(n)=data_reg(ii,jj)
                    dmin=d
                 end if
              end do
           end do
           if(nod_find==1) exit
        end do
     else
        data_mod(n)=(data(1,1)*rt_lon1 + data(2,1)*rt_lon2)*rt_lat1 + &
             (data(1,2)*rt_lon1 + data(2,2)*rt_lon2)*rt_lat2
     end if
  end do
end subroutine interp_2d_field_v2
!
!---------------------------------------------------------------------------
!
subroutine interp_2d_field(num_lon_reg, num_lat_reg, lon_reg, lat_reg, data_reg, &
     num_mod, lon_mod, lat_mod, data_mod, phase_flag)
  !-------------------------------------------------------------------------------------
  ! This routine does 2d interpolation from a regular grid to specified nodes
  ! on the surface grid. The regular grid is assumed to be global.
  ! This version assumes that no dummy values exist in the raw data.
  ! The argument 'phase_flag' is used to set if a phase angle field is to be interpolated.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! data_reg(num_lon_reg, num_lat_reg) 	data on the regular grid
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! phase_flag                          1: interpolate phase angle (0-360°); 0: otherwise
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer, intent(in)         	:: num_lon_reg, num_lat_reg, num_mod
  integer, intent(in)          	:: phase_flag
  real(kind=8) 			:: x, y, diff
  real(kind=8)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=8)                  :: data_ll, data_lh, data_hl, data_hh
  real(kind=8)                  :: data_lo, data_up
  real(kind=8), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=8), intent(in)	:: data_reg(num_lon_reg, num_lat_reg)
  real(kind=8), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod)
  real(kind=8), intent(out)  	:: data_mod(num_mod)
  !
  if(lon_reg(1)<0.0 .or. lon_reg(num_lon_reg)>360.) then
     write(*,*) 'Error in 2D interpolation!'
     write(*,*) 'The regular grid is not in the proper longitude range.'
     call par_ex
     stop
  end if

  do n=1,num_mod
     x=lon_mod(n)
     y=lat_mod(n)
     ! find the surrounding rectangular box and get interpolation ratios
     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0-rt_lat1
     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
        rt_lon1=1.0-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     end if
     !
     data_ll=data_reg(ind_lon_l,ind_lat_l)
     data_lh=data_reg(ind_lon_l,ind_lat_h)
     data_hl=data_reg(ind_lon_h,ind_lat_l)
     data_hh=data_reg(ind_lon_h,ind_lat_h)
     !
     ! interpolate data
     if(phase_flag==1) then   ! interpolate phase value (0,360)
        if(abs(data_ll-data_hl)>180.) then
           if(data_ll<data_hl) then
              data_ll=data_ll+360.
           else
              data_hl=data_hl+360.
           end if
        end if
        if(abs(data_lh-data_hh)>180.) then
           if(data_lh<data_hh) then
              data_lh=data_lh+360.
           else
              data_hh=data_hh+360.
           end if
        end if
        data_lo = data_ll*rt_lon1 + data_hl*rt_lon2                
        data_up = data_lh*rt_lon1 + data_hh*rt_lon2
        if(abs(data_lo-data_up)>180.) then
           if(data_lo<data_up) then
              data_lo=data_lo+360.
           else
              data_up=data_up+360.
           end if
        end if
        data_mod(n)=data_lo*rt_lat1 + data_up*rt_lat2
        if(data_mod(n)>=360.) then
           data_mod(n)=mod(data_mod(n), 360.)
        end if
     else   ! other case
        data_mod(n)=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
             (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2
     end if
  end do
end subroutine interp_2d_field
!
!---------------------------------------------------------------------------
!
subroutine interp_3d_field(num_lon_reg, num_lat_reg, num_lay_reg, &
     lon_reg, lat_reg, lay_reg, data_reg, &
     num_mod, lon_mod, lat_mod, lay_mod, data_mod)
  !-------------------------------------------------------------------------------------
  ! This routine does 3d interpolation from a regular grid to specified nodes.
  ! The regular grid is assumed to be global.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_lon_reg				number of regular grid points in the longitude direction
  ! num_lat_reg				number of regular grid points in the latitude direction
  ! num_lay_reg				number of regular grid points in the vertical direction
  ! lon_reg(num_lon_reg)		longitude of the regular grid points
  ! lat_reg(num_lat_reg)        	latitude of the regular grid points
  ! lay_reg(num_lay_reg)                depth of the regular grid points
  ! data_reg(:,:,:)  	                data on the regular grid
  ! num_mod      			number of interpolation nodes
  ! lon_mod(num_mod)			longitude of interpolation nodes
  ! lat_mod(num_mod)			latitude of interpolation nodes
  ! lay_mod(num_mod)                    depth of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! Unit of lon_reg, lat_reg, lon_mod, lat_mod: degree
  ! Unit of lay_reg, lay_mod: m
  ! Order of lon_reg: monotonically increasing in the range [0 360]  
  !                   (range [180,360] is associated with longitude [-180, 0])
  ! Order of lat_reg: monotonically increasing in the range [-90 90]
  ! Order of lay_reg: monotonically decreasing from surface to bottom
  ! Make sure that 'lon_mod' is also in the [0, 360] range.
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i, flag
  integer			:: ind_lat_h, ind_lat_l, ind_lon_h, ind_lon_l
  integer                       :: ind_lay_h, ind_lay_l
  integer, intent(in)         	:: num_lon_reg, num_lat_reg, num_lay_reg
  integer, intent(in)          	:: num_mod
  real(kind=8) 			:: x, y, z, diff 
  real(kind=8)			:: rt_lat1, rt_lat2, rt_lon1, rt_lon2
  real(kind=8)                  :: rt_lay1, rt_lay2, v_dup, v_dlo
  real(kind=8)                  :: data_ll, data_lh, data_hl, data_hh
  real(kind=8)                  :: v_col(4), z_col(4), H, aux1, aux2
  real(kind=8)                  :: dz, a, b, c, d
  real(kind=8), intent(in)	:: lon_reg(num_lon_reg), lat_reg(num_lat_reg)
  real(kind=8), intent(in)      :: lay_reg(num_lay_reg)
  real(kind=8), intent(in)	:: data_reg(num_lon_reg, num_lat_reg, num_lay_reg)
  real(kind=8), intent(in)	:: lon_mod(num_mod), lat_mod(num_mod), lay_mod(num_mod)
  real(kind=8), intent(out)  	:: data_mod(num_mod)

  do n=1,num_mod
     x=lon_mod(n)
     y=lat_mod(n)
     z=lay_mod(n)

     ! find the surrounding box and get interpolation ratios

     ! 1) north-south direction
     if(y<lat_reg(1)) then
        ind_lat_h=2
        ind_lat_l=1
        y=lat_reg(1)
     elseif(y>lat_reg(num_lat_reg)) then
        ind_lat_h=num_lat_reg
        ind_lat_l=num_lat_reg-1
        y=lat_reg(num_lat_reg)
     else
        do i=2,num_lat_reg
           if(lat_reg(i)>=y) then
              ind_lat_h=i
              ind_lat_l=i-1
              exit
           end if
        end do
     end if
     diff=lat_reg(ind_lat_h)-lat_reg(ind_lat_l)
     rt_lat1=(lat_reg(ind_lat_h)-y)/diff
     rt_lat2=1.0-rt_lat1

     ! 2) east_west direction
     if(x<lon_reg(1)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     elseif(x>lon_reg(num_lon_reg)) then
        ind_lon_h=1
        ind_lon_l=num_lon_reg
        diff=lon_reg(ind_lon_h)+(360.-lon_reg(ind_lon_l))
        rt_lon2=(x-lon_reg(ind_lon_l))/diff
        rt_lon1=1.0-rt_lon2
     else
        do i=2,num_lon_reg
           if(lon_reg(i)>=x) then
              ind_lon_h=i
              ind_lon_l=i-1
              exit
           end if
        end do
        diff=lon_reg(ind_lon_h)-lon_reg(ind_lon_l)
        rt_lon1=(lon_reg(ind_lon_h)-x)/diff
        rt_lon2=1.0-rt_lon1
     end if

     ! 3) up-down direction
     if(z>lay_reg(1)) then
        ind_lay_h=1
        ind_lay_l=2
        z=lay_reg(1)
     elseif(z<lay_reg(num_lay_reg)) then
        ind_lay_h=num_lay_reg-1
        ind_lay_l=num_lay_reg
        z=lay_reg(num_lay_reg)
     else
        do i=2,num_lay_reg
           if(lay_reg(i)<=z) then
              ind_lay_h=i-1
              ind_lay_l=i
              exit
           end if
        end do
     end if
     diff=lay_reg(ind_lay_h)-lay_reg(ind_lay_l)
     rt_lay1=(z-lay_reg(ind_lay_l))/diff
     rt_lay2=1.0-rt_lay1

     data_ll=data_reg(ind_lon_l,ind_lat_l,ind_lay_h)
     data_lh=data_reg(ind_lon_l,ind_lat_h,ind_lay_h)
     data_hl=data_reg(ind_lon_h,ind_lat_l,ind_lay_h)
     data_hh=data_reg(ind_lon_h,ind_lat_h,ind_lay_h)    
     v_dup=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
          (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2

     data_ll=data_reg(ind_lon_l,ind_lat_l,ind_lay_l)
     data_lh=data_reg(ind_lon_l,ind_lat_h,ind_lay_l)
     data_hl=data_reg(ind_lon_h,ind_lat_l,ind_lay_l)
     data_hh=data_reg(ind_lon_h,ind_lat_h,ind_lay_l)    
     v_dlo=(data_ll*rt_lon1 + data_hl*rt_lon2)*rt_lat1 + &
          (data_lh*rt_lon1 + data_hh*rt_lon2)*rt_lat2

     data_mod(n)=v_dup*rt_lay1 + v_dlo*rt_lay2

  end do
end subroutine interp_3d_field
!
!---------------------------------------------------------------------------
!
subroutine interp_vertical_section(num_hor_reg, num_ver_reg, hor_reg, ver_reg, data_reg, &
	num_mod, hor_mod, ver_mod, data_mod)
  !-------------------------------------------------------------------------------------
  ! This routine does 2d interpolation from a regular grid to specified nodes
  ! in a vertical section.
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_hor_reg				number of regular grid points in the horizontal direction
  ! num_ver_reg				number of regular grid points in the vertical direction
  ! hor_reg(num_hor_reg)		horizontal coord. of the regular grid points
  ! ver_reg(num_ver_reg)        	vertical coord. of the regular grid points
  ! data_reg(num_hor_reg, num_ver_reg) 	data on the regular grid
  ! num_mod      			number of interpolation nodes
  ! hor_mod(num_mod)			horizontal coord. of interpolation nodes
  ! ver_mod(num_mod)			vertical coord. of interpolation nodes
  ! data_mod(num_mod)			output data on interpolation nodes
  ! Unit of hor_reg and hor_mod should be the same (in lon or lat unit)
  ! Unit of ver_reg and ver_mod should be the same (usually in m)
  ! Order of hor_reg: monotonically increasing
  ! Order of ver_reg: monotonically decreasing (from ocean top to bottom)
   !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i
  integer			:: ind_ver_h, ind_ver_l, ind_hor_h, ind_hor_l
  integer, intent(in)         	:: num_ver_reg, num_hor_reg, num_mod
  real(kind=8) 			:: x, z, diff
  real(kind=8)			:: rt_ver1, rt_ver2, rt_hor1, rt_hor2
  real(kind=8)                  :: data_ll, data_lh, data_hl, data_hh
  real(kind=8), intent(in)	:: ver_reg(num_ver_reg), hor_reg(num_hor_reg)
  real(kind=8), intent(in)	:: data_reg(num_hor_reg, num_ver_reg)
  real(kind=8), intent(in)	:: ver_mod(num_mod), hor_mod(num_mod)
  real(kind=8), intent(out)  	:: data_mod(num_mod)
  !
  do n=1,num_mod
        x=hor_mod(n)
        z=ver_mod(n)
        ! find the surrounding rectangular box and get interpolation ratios
        ! 1) vertical direction
        if(z>ver_reg(1)) then
           ind_ver_h=2
           ind_ver_l=1
           z=ver_reg(1)
        elseif(z<ver_reg(num_ver_reg)) then
           ind_ver_h=num_ver_reg
           ind_ver_l=num_ver_reg-1
           z=ver_reg(num_ver_reg)
        else
           do i=2,num_ver_reg
              if(ver_reg(i)<=z) then
                 ind_ver_h=i
                 ind_ver_l=i-1
                 exit
              end if
           end do
        end if
        diff=ver_reg(ind_ver_h)-ver_reg(ind_ver_l)
        rt_ver1=(ver_reg(ind_ver_h)-z)/diff
        rt_ver2=1.0-rt_ver1
        ! 2) east_west or north_south direction
        if(x<hor_reg(1)) then
           ind_hor_h=2
           ind_hor_l=1
	   x=hor_reg(1)
        elseif(x>hor_reg(num_hor_reg)) then
           ind_hor_h=num_hor_reg
           ind_hor_l=num_hor_reg-1
	   x=hor_reg(num_hor_reg)
        else
           do i=2,num_hor_reg
              if(hor_reg(i)>=x) then
                 ind_hor_h=i
                 ind_hor_l=i-1
                 exit
              end if
           end do
        end if
        diff=hor_reg(ind_hor_h)-hor_reg(ind_hor_l)
        rt_hor1=(hor_reg(ind_hor_h)-x)/diff
        rt_hor2=1.0-rt_hor1  
        !
        data_ll=data_reg(ind_hor_l,ind_ver_l)
        data_lh=data_reg(ind_hor_l,ind_ver_h)
        data_hl=data_reg(ind_hor_h,ind_ver_l)
        data_hh=data_reg(ind_hor_h,ind_ver_h)
        !
        ! interpolate data
        data_mod(n)=(data_ll*rt_hor1 + data_hl*rt_hor2)*rt_ver1 + &
           (data_lh*rt_hor1 + data_hh*rt_hor2)*rt_ver2     
  end do 

end subroutine interp_vertical_section
!
!---------------------------------------------------------------------------
!
subroutine interp_1d(num_reg, x_reg, data_reg, num_mod, x_mod, data_mod)
  ! not tested yet
  !-------------------------------------------------------------------------------------
  ! This routine does 1d interpolation from a regular grid to specified nodes
  ! The output is data_mod, others are input.
  ! Arguments:
  ! num_reg			number of regular grid points
  ! x_reg(num_hor_reg)	        coord. of the regular grid points
  ! data_reg(num_regg) 	        data on the regular grid
  ! num_mod      		number of interpolation nodes
  ! x_mod(num_mod)		coord. of interpolation nodes
  ! data_mod(num_mod)		output data on interpolation nodes
  ! Unit of x_reg and x_mod should be the same (in lon or lat unit)
  ! Order of x_reg: monotonically increasing
  !-------------------------------------------------------------------------------------
  implicit none
  integer             		:: n, i
  integer			:: ind_h, ind_l
  integer, intent(in)         	:: num_reg, num_mod
  real(kind=8) 			:: x, diff
  real(kind=8)			:: rt_1, rt_2
  real(kind=8), intent(in)	:: x_reg(num_reg)
  real(kind=8), intent(in)	:: data_reg(num_reg)
  real(kind=8), intent(in)	:: x_mod(num_mod)
  real(kind=8), intent(out)  	:: data_mod(num_mod)
  !
  do n=1,num_mod
     x=x_mod(n)
     ! find the two surrounding grid points and get interpolation ratios
     if(x<x_reg(1)) then
        ind_h=2
        ind_l=1
        x=x_reg(1)
     elseif(x>x_reg(num_reg)) then
        ind_h=num_reg
        ind_l=num_reg-1
        x=x_reg(num_reg)
     else
        do i=2,num_reg
           if(x_reg(i)>=x) then
              ind_h=i
              ind_l=i-1
              exit
           end if
        end do
     end if
     diff=x_reg(ind_h)-x_reg(ind_l)
     rt_1=(x_reg(ind_h)-x)/diff
     rt_2=1.0-rt_1  
     !
     ! interpolate data
     data_mod(n)=data_reg(ind_l)*rt_1 + data_reg(ind_h)*rt_2
  end do

end subroutine interp_1d
!
!---------------------------------------------------------------------------




subroutine setup_model
  implicit none

  call read_namelist              ! Read Namelists, this should be before clock_init
  call calculate_time_step      ! time step size
  call pre_handling_parameters  ! such as transform units 
  call define_prog_tracer

end subroutine setup_model
!
!-------------------------------------------------------------------
!
subroutine read_namelist
  ! Routine reads namelist files to overwrite default parameters.

  use g_config
  use o_param
  use i_dyn_parms
  use i_therm_parms
  use g_diag
  use g_forcing_param
  use g_parfe
  use g_clock, only: timenew, daynew, yearnew
  implicit none

  character(len=100)   :: nmlfile
  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (20,file=nmlfile)
  read (20,NML=modelname)
  read (20,NML=timestep)
  read (20,NML=clockinit) 
  read (20,NML=paths)
  read (20,NML=initialization)  
  read (20,NML=inout)
  read (20,NML=mesh_def)
  read (20,NML=geometry)
  read (20,NML=calendar)
  close (20)

  nmlfile ='namelist.oce'    ! name of ocean namelist file
  open (20,file=nmlfile)
  read (20,NML=viscdiff)
  read (20,NML=boundary)
  read (20,NML=oce_scheme)
  read (20,NML=denspress)
  read (20,NML=param_freesurf)
  read (20,NML=tide_obc)
  read (20,NML=passive_tracer)
  read (20,NML=age_tracer)
  close (20)

  nmlfile ='namelist.forcing'    ! name of forcing namelist file
  open (20,file=nmlfile)
  read (20,NML=forcing_exchange_coeff)
  read (20,NML=forcing_source)
  read (20,NML=forcing_bulk)
  read (20,NML=land_ice)
  close (20)

#ifdef use_ice
  nmlfile ='namelist.ice'    ! name of ice namelist file
  open (20,file=nmlfile)
  read (20,NML=ice_stress)
  read (20,NML=ice_fric)
  read (20,NML=ice_rheology)
  read (20,NML=ice_scheme)
  read (20,NML=ice_therm)
  close (20)
#endif

  nmlfile='namelist.diag'    ! name of diagnose namelist file
  open(20, file=nmlfile)
  read(20,NML=diag_flag)
  close(20)

  if(mype==0) write(*,*) 'Namelist files are read in'
end subroutine read_namelist
!
!-------------------------------------------------------------------
!
subroutine calculate_time_step
  use g_config
  use g_parfe
  implicit none

  ! compute dt and dt_inv
  dt=86400./float(step_per_day)
  dt_inv=1.0/dt    
  if(mype==0) write(*,*) 'time step size is set to ', real(dt,4), 'sec'

end subroutine calculate_time_step
!
!-------------------------------------------------------------------
!
subroutine pre_handling_parameters
  use o_param
  use g_config
  use g_parfe
  implicit none

  ! change parameter unit: degree to radian
  domain_length=domain_length*rad
  if(rotated_grid) then 
     alphaEuler=alphaEuler*rad 	
     betaEuler=betaEuler*rad
     gammaEuler=gammaEuler*rad
  end if

end subroutine pre_handling_parameters
!
!-------------------------------------------------------------------
!
subroutine define_prog_tracer
  use o_param
  use o_array, only : prog_tracer_name
  use g_parfe
  implicit none

  integer      :: j, num
  character(1) :: cageind
  character(4) :: tr_name

  ! allocate prog_tracer_name

  num_tracer=2  ! t and s
  if(use_passive_tracer) then
     num_tracer=num_tracer+num_passive_tracer
  end if
  if(use_age_tracer) then
     num_tracer=num_tracer+num_age_tracer
  end if

  allocate(prog_tracer_name(num_tracer))

  ! fill prog_tracer_name

  num=2  ! t and s

  prog_tracer_name(1)='temp'
  prog_tracer_name(2)='salt'

  if(use_passive_tracer) then
     do j=1,num_passive_tracer
        write(cageind,'(i1)') j
        tr_name='ptr'//cageind
	prog_tracer_name(num+j)=tr_name
     end do
     num=num+num_passive_tracer
  end if

  if(use_age_tracer) then
     do j=1,num_age_tracer
        write(cageind,'(i1)') j
        tr_name='age'//cageind
	prog_tracer_name(num+j)=tr_name
     end do
     num=num+num_age_tracer
  end if

  if(mype==0) write(*,*) 'Number of prognostic ocean tracers: ',num

end subroutine define_prog_tracer
!
!-------------------------------------------------------------------
!
subroutine get_run_steps
  use g_config
  use g_clock
  use g_parfe
  implicit none

  integer      :: i, temp_year, temp_mon, temp_fleapyear

  ! clock should have been inialized before calling this routine

  if(run_length_unit=='s') then
     nsteps=run_length
  elseif(run_length_unit=='d') then
     nsteps=step_per_day*run_length
  elseif(run_length_unit=='m') then
     nsteps=0
     temp_mon=month-1
     temp_year=yearnew
     temp_fleapyear=fleapyear
     do i=1,run_length
        temp_mon=temp_mon+1
        if(temp_mon>12) then
           temp_year=temp_year+1
           temp_mon=1
           call check_fleapyr(temp_year, temp_fleapyear)
        end if
        nsteps=nsteps+step_per_day*num_day_in_month(temp_fleapyear,temp_mon)
     end do
  elseif(run_length_unit=='y') then
     nsteps=0
     do i=1,run_length
        temp_year=yearnew+i-1
        call check_fleapyr(temp_year, temp_fleapyear)
        nsteps=nsteps+step_per_day*(365+temp_fleapyear)
     end do
  else
     write(*,*) 'Run length unit ', run_length_unit, ' is not defined.'
     write(*,*) 'Please check and update the code.'
     call par_ex
     stop
  end if

  if(mype==0) write(*,*) nsteps, ' steps to run for this (', runid, ') job submission'
end subroutine get_run_steps
subroutine compute_mixlay
  use o_MESH
  use o_PARAM
  use o_array
  use g_meanarrays
  use g_PARFE
  implicit none

  integer         :: m, n2, n3, k, row
  real(kind=8)    :: dens_surf, dens, md, bf, bf_up 
  real(kind=8)    :: buoyancy_crit, smallvalue

  smallvalue=1.0e-20
  buoyancy_crit=0.0003

  !mixed layer depth
  do n2=1,ToDim_nod2d
     row=nod3D_below_nod2D(1,n2)  
     call fcn_dens0(tracer(row,1),tracer(row,2),dens_surf)
     md=0.0
     bf_up=0.0
     do k=2,num_layers_below_nod2d(n2)+1
        n3=nod3d_below_nod2d(k,n2)
        call fcn_dens0(tracer(n3,1),tracer(n3,2),dens)
        bf=g*(dens-dens_surf)/dens
        if(bf>=buoyancy_crit) then
           md=md+(coord_nod3d(3,n3)-md)/(bf-bf_up+smallvalue)*(buoyancy_crit-bf_up)
           exit
        else
           md=coord_nod3d(3,n3)
           bf_up=bf
        end if
     end do
     mixlay_dep_mean(n2)=mixlay_dep_mean(n2) + abs(md)     
  end do


  !!mixed layer depth
  !do n2=1,ToDim_nod2d
  !   row=nod3D_below_nod2D(1,n2)           
  !   md=0.0
  !   bf_up=0.0
  !   do k=2,num_layers_below_nod2d(n2)+1
  !      n3=nod3d_below_nod2d(k,n2)
  !      call fcn_density(tracer(row,1),tracer(row,2),coord_nod3d(3,n3),dens_surf)
  !      bf=g*(density_insitu(n3)-dens_surf)/density_insitu(n3)
  !      if(bf>=buoyancy_crit) then
  !         md=md+(coord_nod3d(3,n3)-md)/(bf-bf_up+smallvalue)*(buoyancy_crit-bf_up)
  !         exit
  !      else
  !         md=coord_nod3d(3,n3)
  !         bf_up=bf
  !      end if
  !   end do
  !   mixlay_dep_mean(n2)=mixlay_dep_mean(n2) + abs(md)     
  !end do

end subroutine compute_mixlay
!
!--------------------------------------------------------------------------------------------
!
subroutine process_elem2node(id,array_3d)
  ! Average elementwise variables to nodes
  ! Used for output ocean SGS parameterization
  ! Called in write_means_part2
  use o_mesh
  use o_elements
  use o_array
  use g_meanarrays
  use g_parfe
  implicit none

  integer       :: id,  elem, elnodes(4), row
  real(kind=8)  :: array_3d(nod3D)
  real(kind=8), allocatable  :: aux(:)  

  allocate(aux(myDim_nod3D+eDim_nod3D)) 
  aux=0.0 

  if(id==1) then
     do elem=1,myDim_elem3d       
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_u(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d         
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d) 
  elseif(id==2) then
     do elem=1,myDim_elem3d          
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_v(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d           
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d) 
  elseif(id==3) then
     do elem=1,myDim_elem3d          

        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_ut(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d           
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume  
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==4) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_vt(elem)*voltetra(elem)
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume 
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==5) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_us(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume 
     end do
     call broadcast3d(aux,array_3d)  
  elseif(id==6) then
     do elem=1,myDim_elem3d         
        elnodes=elem3d_nodes(:,elem)
        aux(elnodes)=aux(elnodes)+sgs_vs(elem)*voltetra(elem) 
     end do
     do row=1,myDim_nod3d          
        aux(row)=aux(row)/wrhs(row)  !wrhs contains cluster volume
     end do
     call broadcast3d(aux,array_3d)  
  end if
  deallocate(aux)
end subroutine process_elem2node
subroutine config_remind_warning_info
  ! check configuration and options and provide reminding or warning
  use o_param
  use o_array          
  use o_solver
  use g_config
  use g_parfe
  use g_clock
  use g_forcing_param
  use g_forcing_index
  use g_forcing_arrays
  implicit none

  if(mype/=0) return

#ifdef use_opbnd_restoring
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Reminding:'
  write(*,*) 'It is configured to restore velocity at open boundarie.'
  write(*,*) 'Routine init_restoring_vel should be specifically prepared.'
  write(*,*) '---------------------------------------------------'
#endif

  if(buffer_zone) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'It is configured to restore tracers at open boundaries.'
     write(*,*) 'Routine init_restoring_bufferzone should be specifically'
     write(*,*) 'prepared.'
     write(*,*) '---------------------------------------------------'
  end if

#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
  if(.not.buffer_zone) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Warning:'
     write(*,*) 'It is recommended to restore tracers in a buffer zone near o.b.'
     write(*,*) '---------------------------------------------------'
  end if
#endif

#ifdef use_opbnd_tide
#ifndef use_semiimplicit_scheme
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Warning:'
  write(*,*) 'It is recommended to use the semi-implicit scheme when'
  write(*,*) 'simulating tides. To change, set it in the Makefile.'
  write(*,*) '---------------------------------------------------'
#endif
#endif

#ifdef use_ice
  if(restore_s_surf>0.0) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'It is specified to restore SSS. Check which climatology is'
     write(*,*) 'to be used (yearly or monthly, which source) and modify the'
     write(*,*) 'code in file gen_forcing_couple.F90 for your purpose.'
     write(*,*) 'Also check how restoring under ice is done in the code;'
     write(*,*) 'you might need to change the code for your application.'
     write(*,*) '---------------------------------------------------'
  end if
#endif

  if(balance_salt_water) then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'You specified to balance salt and freshwater fluxes.'
     write(*,*) 'The default is to do correction for every time step.'
     write(*,*) 'You might need your own correction strategy. Check'
     write(*,*) 'the code to be sure it does what you want.'
     write(*,*) '---------------------------------------------------'
   end if  

#ifdef use_sw_pene
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Reminding:'
  write(*,*) 'It is specified to consider shortwave penetration. The'
  write(*,*) 'chlorophyll climatology on meshes should have been prepared'
  write(*,*) 'offline (unformated). Automatic interpolation to model grids is'
  write(*,*) 'not supported in the code!'
  write(*,*) '---------------------------------------------------'
#endif

#ifndef use_ice
  write(*,*) '---------------------------------------------------'
  write(*,*) 'Warning:'
  write(*,*) 'You are running the ocean-alone model. The surface forcing routine'
  write(*,*) 'should be adjusted/checked to properly apply your particular surface'
  write(*,*) 'forcing to the ocean.'
  write(*,*) '---------------------------------------------------'
#endif

  if(wind_data_source=='NCEP' .and. ncar_bulk_formulae) then
      write(*,*) '---------------------------------------------------'
      write(*,*) 'Warning:'
      write(*,*) 'you are using NCEP 2m air temperature and humidity. The current'
      write(*,*) 'formulae for calculating drag and heat exchange coefficients'
      write(*,*) 'only support 10m data. If you plan to use these formulae, a small'
      write(*,*) 'update is required. Before doing it, turn off ncar_bulk_formulae!'
      write(*,*) '---------------------------------------------------'
   end if

  if(mix_scheme=='MY2p5') then
     write(*,*) '---------------------------------------------------'
     write(*,*) 'Reminding:'
     write(*,*) 'MY2.5 mixing scheme is to be used. Its variables will'
     write(*,*) 'be saved at the end of the run for the next restart job.'
     write(*,*) 'This file will be replaced in each restart run. So it'
     write(*,*) 'is not possible to restart from intermediate snapshots,'
     write(*,*) 'or from previous runs if the file is not backed up'
     write(*,*) 'manually. To get rid of this limit, code need update.'
     write(*,*) '---------------------------------------------------'
  end if

  if(use_passive_tracer) then
     if(ptr_start_year>yearnew) then
        write(*,*) '---------------------------------------------------'
        write(*,*) 'Warning:'
        write(*,*) 'You specify to use passive tracers not at the beginning of'
        write(*,*) 'this job. This is not supported. The model will start to'
        write(*,*) 'include the passive tracers from the beginning of this job.'
        write(*,*) 'If you do not want this, cancel this job, turn off'
        write(*,*) 'use_passive_tracer, run the model until the moment when'
        write(*,*) 'you want to have the passive tracers, and then re-start the'
        write(*,*) 'simulation with passive tracers used.'
        write(*,*) '---------------------------------------------------'
     end if
  end if

  if(use_age_tracer) then
     if(age_tracer_start_year>yearnew) then
        write(*,*) '---------------------------------------------------'
        write(*,*) 'Warning:'
        write(*,*) 'You specify to use age tracers not at the beginning of'
        write(*,*) 'this job. This is not supported. The model will start to'
        write(*,*) 'include the age tracers from the beginning of this job.'
        write(*,*) 'If you do not want this, cancel this job, turn off'
        write(*,*) 'use_age_tracer, run the model until the moment when'
        write(*,*) 'you want to have the age tracers, and then re-start the'
        write(*,*) 'simulation with age tracers used.'
        write(*,*) '---------------------------------------------------'
     end if
  end if
end subroutine config_remind_warning_info
!
!--------------------------------------------------------------------------
!
subroutine check_blowup
  !check if the model blows up and cancel the job if it is the case
  !Salinity is used as an indicator.
  !ALLREDUCE is used, so it slows down the code. 
  !A better way is requried!!! Qiang
  use o_array          
  use g_parfe
  implicit none

  integer     :: flag, g_flag

  flag=0
  if(any(tracer(:,2)<0.0)) then
     flag=1
  end if
   
  g_flag=0
  call MPI_AllREDUCE(flag, g_flag, 1, MPI_INTEGER, MPI_SUM, &
       MPI_COMM_WORLD, MPIerr)
  
  if(g_flag>0) then
     if(mype==0) then
        write(*,*) 'Negative salinity found. The model blows up!'
        write(*,*) 'The program will be forced to stop here.'
     end if
     
     call par_ex
     stop
  end if
end subroutine check_blowup
subroutine dist_on_earth(lon1, lat1, lon2, lat2, dist)
  ! distance on the earth between two points
  ! input: lon1 lat2 and lon2 lat2 in radian
  ! output: dist in m
  use o_param
  implicit none
  real(kind=8)  :: lon1, lat1, lon2, lat2, alpha, dist

  alpha=acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
  dist=r_earth*abs(alpha)

end subroutine dist_on_earth
program main
  !=============================================================================!
  !
  !                     Finite Element Ocean Model
  !
  !=============================================================================!
  !                      The main driving routine
  !=============================================================================!         

  use g_config
  use o_param
  use o_array          
  use o_solver
  use o_mixing_kpp_mod
  use g_parfe
  use g_clock
  use g_forcing_index
  use g_forcing_param
  use g_forcing_arrays
  use g_diag
  use g_forcing_interp
#ifdef use_ice
  use i_array
#endif
  implicit none


  ! MPI initialization
  call par_init                 ! initializes MPI
  if(mype==0) write(*, *) 'Running on ', npes, ' PEs'


  if(mype==0) write(*,*) '*************************************************************'
  ! read namelist, initialize clock, prepare basic configuration etc.

  call setup_model             ! setup basic config, do it before clock_init
  call clock_init               ! read the clock file
  call get_run_steps            ! compute total steps to run 
  call config_remind_warning_info


  if(mype==0) write(*,*) '*************************************************************'
  ! mesh and communication buffers

  call ocean_mesh_setup         ! setup the 2D/3D meshes
  call set_par_support
  call find_cluster_area        ! cluster area for 2D nodes

  if(mype==0) write(*,*) '*************************************************************'
  ! ocean: matrices, arrays, initialization, buffer zone, tide etc.

  call ocean_matrices_setup     ! Builds matrices and call partitioning
  call ocean_array_setup        ! allocate ocean arrays 
  call ocean_init               ! initialize the oce or read restart files

  if(use_ref_density) then
     call compute_ref_density   ! Fills in ocean reference density 
  endif

#ifdef use_opbnd_restoring
  call init_restoring_vel
#endif

  if(buffer_zone) then
     call init_restoring_bufferzone
  end if

#ifdef use_opbnd_tide
  call init_tidal_opbnd         ! initialize tidal ocean open boundary
#endif

#ifdef use_ice
  if(mype==0) write(*,*) '*************************************************************'
  ! ice: matrices, arrays, initialization

  call ice_matrices_setup       ! Build ice matrices
  call ice_array_setup          ! allocate ice arrays, setup ice adv matrix
  call ice_init                 ! initialize the ice or read restart files
#endif


  if(mype==0) write(*,*) '*************************************************************'
  ! forcing: arrays, initialization, interpolation preparation  

#ifdef use_ice
  call forcing_array_setup
  call init_forcing_interp      ! calculates the forcing interpolation weights
  call init_atm_forcing         ! initialize forcing fields
#else
#ifndef toy_ocean
  call forcing_array_setup_OnlyOcean
  call init_forcing_interp 
  call init_atm_forcing_OnlyOcean 
#endif 
#endif

  if(use_landice_water) call landice_water_init


  if(mype==0) write(*,*) '*************************************************************'
  ! init mean arrays, adjust clock and create output files if required  

#if defined(allow_calcmeans) || defined(allow_diag)
  call init_meanarrays          ! allocate arrays for mean fields
#endif
  call clock_newyear  		! check if it is a new year
  call init_output              ! create new output files


#ifdef use_fullfreesurf
  if(mype==0) write(*,*) '*************************************************************'
  ! updating mesh and matrices in case of full free surface setup 

  if(any(abs(ssh)>small)) then
     call update_mesh
     call update_matrices
     call update_mesh
  endif
#endif

  ! set some flags for solvers
  iter_first=.true.             ! iter_first & iteruv_first should be 'true' at start
  iteruv_first=.true.


  if(mype==0) then
     write(*,*) '*************************************************************'
     write(*,*) 'iteration starts ...'
     write(*,*) '*************************************************************'
     write(*,*)
  end if

  !   call output(1)

  ! preparation done
  !----------------------------------------------------------------------------
  ! start iteration

  do istep=1, nsteps   
     call clock
     call init_output

     call forcing_index 

#ifdef use_ice    
     call ocean2ice
     call update_atm_forcing
     call ice_step
     call ice2ocean
     if(use_landice_water) call add_landice_water
#else
#ifndef toy_ocean
     call update_atm_forcing_OnlyOcean
#endif
#endif

#ifdef use_fullfreesurf 
     if(balance_salt_water) call check_imb_freshwater
#endif

#ifdef use_sw_pene
     call cal_shortwave_rad
#endif

     call ocean_step

#if defined(allow_calcmeans) || defined(allow_diag)
     call add2meanarrays      
#endif

     ! save (NetCDF)
     call output(0) 

!!$     ! save (ascii, for an easy debugging during new setup test phase)
!!$     if(mod(istep,step_per_day*5)==0) then
!!$        call oce_out     
!!$#ifdef use_ice
!!$        call ice_out
!!$#endif
!!$     end if

     if(mod(istep,logfile_outfreq)==0) then	
        ! log file output (for debugging during new setup test phase)
        !write(*,*) 'uf max', maxval(abs(uf)), istep
        !write(*,*) 's min', minval(tracer(:,2)), istep

        if (mype==0) then
           write(*,*) 'Step', istep, '  day', daynew, '  year', yearnew
  	   write(*,*)
        endif
     end if

     if(check_run_state) call check_blowup	! check if the program blows up

  end do

  ! iteration done
  !--------------------------------------------------------------------
  ! some finishing-up routines follow


  !#ifndef use_ice
  ! call close_ocean_forcing	! required in special cases
  !#endif

  if(mix_scheme=='MY2p5') call save_MY_vara  ! save MY2.5 variables for next restart

  ! save (ascii, for an easy debugging during new setup test phase)
  !call oce_out
#ifdef use_ice
  !call ice_out
#endif

  call clock_finish		! save clock

  if(mype==0) then
     open(unit=50, file='goodfile')
     write(50,*)'go on'
     close(50)
  end if

  call par_ex 			! finalizes MPI
  write(*,*) 'Experiment '//runid//' successfully completed'

end program main


