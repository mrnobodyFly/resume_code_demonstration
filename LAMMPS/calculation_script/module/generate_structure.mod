# gsc_method: scratch | data | dump | restart
# dump_fname, data_fname, restart_fname
# dumpfile_step

# call special init command before creating simulation box
shell ${lmp_t_dir}/../shell/util/check_special_init ${pot_prefix}/init.mod
include special_init.mod

if "${gsc_method} == scratch" then &
   "jump SELF sc_from_scratch" &
 elif "${gsc_method} == data" &
   "jump SELF sc_from_data" &
 elif "${gsc_method} == dump" &
   "jump SELF sc_from_dump" &
 elif "${gsc_method} == restart" &
   "jump SELF sc_from_restart"
#-----------------------------------------------------
# generate structure from scratch
#-----------------------------------------------------
label sc_from_scratch
# read lattice parameter
# include ${pot_prefix}/unit_cell/${atom_type}/${xtal_structure}/lat_param_${atom_type}_${xtal_structure}_static_0K.mod
# variable lat_a equal v_lat_${xtal_structure}_a_${atom_type}
# variable lat_k equal v_lat_${xtal_structure}_k_${atom_type}
if $(is_defined(variable,change_atom)) then &
   "variable atom_comp string ${atom_type}_${change_atom}" &
else &
   "variable atom_comp string ${atom_type}"
include ${pot_prefix}/unit_cell/${atom_comp}/${xtal_structure}/lat_param_${atom_comp}_${xtal_structure}_static_0K.mod
variable lat_a equal v_lat_${xtal_structure}_a_${atom_comp}
variable lat_k equal v_lat_${xtal_structure}_k_${atom_comp}
# orientation
include ${lmp_t_dir}/module/build_model/${xtal_structure}/general_orient/${orientation}.mod # specify orientation
include ${lmp_t_dir}/module/build_model/rotation_matrix.mod # calculate rotation matrix
# lattice
include ${lmp_t_dir}/module/build_model/${xtal_structure}/lattice_standard.mod # define lattice vectors, basis vectors
include ${lmp_t_dir}/module/build_model/${xtal_structure}/lattice_rotated.mod # create lattice
# box
include ${lmp_t_dir}/module/build_model/periodic_box_dimension.mod # define periodic box dimension
# region box prism &
# 0 ${supercell_lx} 0 ${supercell_ly} 0 ${pbox_lz} &
# 0 ${pbox_xz} ${pbox_yz} units box
variable box_lx equal ${supercell_lx}+20
variable box_ly equal ${supercell_ly}+20
region box prism &
0 ${box_ly} 0 ${box_ly} 0 ${pbox_lz} &
0 ${pbox_xz} ${pbox_yz} units box
create_box ${nt_atom} box
# potential
include ${pot_prefix}/potential.mod
# create atoms
# region atom_bk prism &
# 0 ${supercell_lx} 0 ${supercell_ly} 0 ${pbox_lz} &
# 0 ${pbox_xz} ${pbox_yz} units box
region atom_bk prism &
0 ${box_lx} 0 ${box_ly} 0 ${pbox_lz} &
0 ${pbox_xz} ${pbox_yz} units box
create_atoms 1 region atom_bk
# change atom type for alloys
include change_atom.mod # created by run.sh
jump SELF end_of_gsc

#-----------------------------------------------------
# generate structure from data file
#-----------------------------------------------------
label sc_from_data
read_data ${data_fname} # data_fname=../../data/data.relaxed_0Pa
include ${pot_prefix}/potential.mod
jump SELF end_of_gsc

#-----------------------------------------------------
# generate structure from dump file
#-----------------------------------------------------
label sc_from_dump
# read_dump must be used after simulation box is created
region box block 0 1 0 1 0 1 units box
create_box ${nt_atom} box
# change_box all triclinic
read_dump ${dump_fname} ${dumpfile_step} x y z box yes timestep no purge yes add keep format native
include ${pot_prefix}/potential.mod
jump SELF end_of_gsc

#-----------------------------------------------------
# generate structure from restart file
#-----------------------------------------------------
label sc_from_restart
read_restart ${restart_fname} # restart_fname=../../data/restart.relaxed_0Pa
include ${pot_prefix}/potential.mod
jump SELF end_of_gsc

label end_of_gsc

# run 0 # run 0 will discard all previously defined fix
