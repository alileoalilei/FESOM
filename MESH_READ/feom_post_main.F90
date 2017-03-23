program feom_post_main
  use g_config
  use o_param
  use o_elements
  use o_mesh
  use o_array
  use o_read
  use g_clock
  use g_parfe
  use FEOM_FORC_DIAG 
  use FEOM_POST_DIAG     
  use FEOM_DART_DIAG
  use FEOM_OBS_DIAG 


implicit none

  real              :: t0, t1, t2, t3, t4, t5, t6, &
                       t7, t8, t9, t10, t11, t12,t13, &
                       t14,t15,t16, t17


  call read_namelist
  call cpu_time(t0); print*, 'time elapsed:',t0
  call read_elem
  call cpu_time(t1); print*, 'time elapsed:',t1
  call read_node
  call cpu_time(t2); print*, 'time elapsed:',t2
  call read_aux3
  call cpu_time(t3); print*, 'time elapsed:',t3
  call read_depth
  call cpu_time(t4); print*, 'time elapsed:',t4
  call ocean_array_setup
  call cpu_time(t5); print*, 'time elapsed:',t5
  call ocean_mesh_setup
  call cpu_time(t6); print*, 'time elapsed:',t6
  call find_cluster_area
  call cpu_time(t7); print*, 'time elapsed:',t7
  if ( tool.eq.1 ) then
  call basin_mean_evolution
  else if ( tool.eq.2 ) then
  call READ_THALWEG_FROM_NC
  else if ( tool.eq.3 ) then
  call READ_SECTION_FROM_NC
  else if ( tool.eq.4 ) then
  call marmara_mean_evolution
  else if ( tool.eq.5 ) then
  call CALC_SECTION_MONTHLY_MEAN 
  else if ( tool.eq.6 ) then
  call CALC_THALWEG_MONTHLY_MEAN 
  else if ( tool.eq.7 ) then
  call READ_ENSEMBLE_FROM_NC
  else if ( tool.eq.8 ) then
  call SYNTHETIC_FERRYBOX_FROM_NR
  else if ( tool.eq.9 ) then
  call READ_SECTION_FROM_INO
  else if ( tool.eq.10 ) then
  call READ_CTD_DATA
  else if ( tool.eq.11 ) then
  call PROFILE_FROM_NC
  else if ( tool.eq.12 ) then
  call velocity_at_the_exit
  else if ( tool.eq.13) then
  call READ_SECTION_FROM_INC
  else if ( tool.eq.14 ) then
  call dardanelles_for_MFS
  else if ( tool.eq.15 ) then
  call total_kinetic_energy
  else if ( tool.eq.16 ) then
  call surface_kinetic_energy
  call cpu_time(t16)
  else if ( tool.eq.17 ) then
  call CALC_SECTION_ANNUAL_MEAN 
  else if ( tool.eq.18 ) then
  call CALC_THALWEG_ANNUAL_MEAN 
  else if ( tool.eq.19 ) then
  call compute_vorticity
  else if ( tool.eq.20 ) then
  call forcing_array_setup
  call compute_wind_stress_curl
  else if ( tool.eq.21 ) then
  call compute_net_flux
  else if ( tool.eq.22 ) then
  call forcing_array_setup
  call compute_surface_buoyancy
  else if ( tool.eq.23 ) then
  call forcing_array_setup
  call compute_forcing_monthly_timeseries
  else if ( tool.eq.24 ) then
  call forcing_array_setup
  call compute_wind_work
  else if ( tool.eq.25 ) then
  call READ_SHIP_TRACK
  end if
end program feom_post_main
