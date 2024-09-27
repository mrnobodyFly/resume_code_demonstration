reset_timestep 0
variable is_converge equal 0
variable opt_i loop ${opt_loop_size}
dump dump_opt ${dump_group} custom ${dump_interval} ./dump_files/${dump_fname}.* id type xu yu zu v_group_label
label struct_opt_loop
  include ${lmp_t_dir}/e_min/fire_min.mod
  include ${lmp_t_dir}/e_min/cg_min.mod
  if "${fmax} < ${fmax_tol}" then &
     "variable is_converge equal 1" &
     "print '${step} ${fmax} ${fnorm} ${is_converge}' append ${min_result_fname}" &
     "jump SELF struct_opt_break"
  next opt_i
  jump SELF struct_opt_loop
  print "${step} ${fmax} ${fnorm} ${is_converge}" append ${min_result_fname}
label struct_opt_break
undump dump_opt
variable opt_i delete
