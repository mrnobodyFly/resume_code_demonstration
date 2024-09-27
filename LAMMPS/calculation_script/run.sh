#!/usr/bin/env bash
# dir: basis_composition/lattice_pctl_changeatom/temperature/size
. $HOME/template/shell/bash/config
. $HOME/template/shell/bash/string
. $HOME/template/shell/bash/lmp_pot_spec
#=============================================
#=============================================
# parsing variables, going down directories
#============================================
run() {
    if [ $# -eq 1 ]; then
        local root_path=$1
    else
        local root_path=./
    fi
    local cwd_parsing=$(pwd)
    for crack_calc in ${crack_calc_all[@]}; do
        local index_tok=($(str_tokenizer $crack_calc ':'))
        local crack_tok=($(str_tokenizer ${crack[${index_tok[0]}]} ':'))
        local size_tok=($(str_tokenizer ${size[${index_tok[1]}]} '_'))
        # local dirname_tok=($(str_tokenizer ${dir_name[${index_tok[2]}]} ':'))
        local k_load_dir=${dir_name[${index_tok[2]}]}
        local k_tok=($(str_tokenizer ${k[${index_tok[3]}]} ':'))
        # local restart_tok=($(str_tokenizer ${restart[${index_tok[4]}]} ':'))
        if [[ ${index_tok[4]} == x ]]; then
            local restart_tok=(0 restart.crack)
        else
            local restart_tok=($(str_tokenizer ${restart[${index_tok[4]}]} ':'))
        fi
        local ovito_plot_files=(${index_tok[@]:5})

        # crack_tok
        local atom_type=${crack_tok[0]}
        local xtal_structure_c1_c2_c3=${crack_tok[1]}
        local crack_type=${crack_tok[2]}
        local n=${#crack_tok[@]}
        if (( $n>=4 )); then
            local fc1_fc2=${crack_tok[3]}
            local fc1_fc2_tok=($(str_tokenizer $fc1_fc2 '_'))
            local frac_c1=${fc1_fc2_tok[0]}; local frac_c2=${fc1_fc2_tok[1]};
        else 
            local frac_c1=0; local frac_c2=0;
        fi
        if (( $n>=5 )); then
            local inclined_cancel=${crack_tok[4]}
        else
            local inclined_cancel=0
        fi

        # size_tok
        local supercell_lx=${size_tok[0]}
        local supercell_ly=${size_tok[1]}
        local supercell_nz=${size_tok[2]}

        # k_tok
        local kii0_ki0_kiii0=${k_tok[0]}
        local dkii_dki_dkiii=${k_tok[1]}
        local k_load_step=${k_tok[2]}

        # restart_tok
        local is_relaxed=${restart_tok[0]}
        local restart_fname=${restart_tok[1]}
        if [[ "${is_relaxed}" -eq 1  && "${restart_fname}" == restart.k_load_* ]]; then
            local kii0_ki0_kiii0=${restart_fname#restart.k_load_}
            echo "------------------------------------------------------"
            echo "Restart calculation using ${restart_fname}"
            echo "K0=${kii0_ki0_kiii0}, dK=${dkii_dki_dkiii}, k_load_step=${k_load_step}"
        fi

        # xtal_structure_c1_c2_c3
        local lat_c_tok=($(str_tokenizer ${xtal_structure_c1_c2_c3} '_'))
        local xtal_structure=${lat_c_tok[0]}
        local c1=${lat_c_tok[1]}
        local c2=${lat_c_tok[2]}
        local c3=${lat_c_tok[3]}

        # kii_ki_kiii tok
        local ki_tok=($(str_tokenizer $kii0_ki0_kiii0 '_'))
        local kii0=${ki_tok[0]}
        local ki0=${ki_tok[1]}
        local kiii0=${ki_tok[2]}

        local dki_tok=($(str_tokenizer $dkii_dki_dkiii '_'))
        local dkii0=${dki_tok[0]}
        local dki0=${dki_tok[1]} 
        local dkiii0=${dki_tok[2]}

        atom_type

        cd ${cwd_parsing}
    done
}
atom_type() {
	  local cwd_atom_type=$(pwd)

    local atom_type_dir=${atom_type}
    prepare_run_path ${atom_type_dir} 'yes'
    cd ${atom_type_dir}
    lattice_orient
    cd ${cwd_atom_type}

}
lattice_orient() {
    local cwd_lattice=$(pwd)
    # local lattice_dir=${xtal_structure_c1_c2_c3}
    if [ -z ${change_atom+x} ]; then
        local lattice_dir=${xtal_structure_c1_c2_c3}
        local atom_comp=${atom_type}
        local req_change_atom=false
    else
        local lattice_dir=${xtal_structure_c1_c2_c3}_${change_atom}
        local atom_comp=${atom_type}_${change_atom}
        local req_change_atom=true
    fi

    prepare_run_path ./${lattice_dir} 'yes'
    cd ./${lattice_dir}
    size

    cd ${cwd_lattice}
}
size() {
    local cwd_size=$(pwd)
    local size_dir=${supercell_lx}A_${supercell_ly}A_${supercell_nz}
    prepare_run_path ./${size_dir} 'yes'
    cd ./${size_dir}
    k_load

    cd ${cwd_size}
}
k_load() {
    local cwd_k_load=$(pwd)
    local k_load_dir=${crack_type}_${frac_c1}_${frac_c2}/${k_load_dir}
    prepare_run_path ./${k_load_dir} 'yes'
    cd ./${k_load_dir}

    local crack_type_full=$crack_type
    # parsing crack_type
    local crack_type_tok=($(str_tokenizer ${crack_type} '_'))
    local crack_type=${crack_type_tok[0]}
    local crack_nheight=0; local cancel_nheight=2; local plane_normal1=null; local plane_normal2=null;
    local n=${#crack_type_tok[@]}

    if (( $n>=2 )); then
        local crack_nheight=${crack_type_tok[1]}
    fi
    if (( $n>=3 )); then
        local cancel_nheight=${crack_type_tok[2]}
    fi
    if (( $n>=4 )); then
        local plane_normal1=${crack_type_tok[3]}
    fi
    if (( $n>=5 )); then
        local plane_normal2=${crack_type_tok[4]}
    fi

    # compute stroh tensor
    if [[ -f ./data/orientation.dat && ! -f ./data/stroh_tensor.dat ]]; then
        echo "Computing Stroh tensor"
        compute_stroh
    fi
    if [[ -f ./data/orientation.dat &&  ! -f ${cwd_size}/orientation.dat ]]; then
        rsync -avP ./data/orientation.dat ${cwd_size}/orientation.dat
        rsync -avP ./data/stroh_tensor.dat ${cwd_size}/stroh_tensor.dat
    fi

    mode

    cd ${cwd_k_load}
}

mode() {
    case "$calc_mode" in
        build)
            local lmp_infile=build_crack.lammps
            if [[ -f restart.crack ]]; then
                echo "Warning!!! restart.crack already exists !!!"
            fi
            ;;
        load)
            local lmp_infile=k_load.lammps
            ;;
        *)
            send_msg "Unknown calc_mode: $calc_mode" 1
            return
            ;;
    esac
    case "$mode" in
        prep)
            prepare_run_file
            ;;
        calc)
            calc
            ;;
        nonlinear)
            echo "Computing nonlinear displacement"
            echo $(pwd)
            nonlinear_displacement
            ;;
        ovito)
            echo $(pwd)
            echo "Plotting ovito figure"
            # ovito_plot
            # python3 ${root_path}/python/modeII_twin_growth.py
            # local tip_type=sharp
            # local tip_type=one_step
            local tip_type=blunt
            rsync -aP ${ovito_plot_files[@]} ~/Documents/article/bcc_crack/core_data/k_load/ovito_image/${atom_type}/${c2}_${c3}/${tip_type}/
            rsync -aP ./data/crack_info.mod dump.crack ~/Documents/article/bcc_crack/core_data/k_load/ovito_image/${atom_type}/${c2}_${c3}/${tip_type}/
            ;;
        view)
            echo $(pwd)
            echo "Ovito visualization"
            # ovito3.6 dump.k_load_*_${ki0}_0.0000 > /dev/null 2>&1
            # ovito3.6 dump.crack dump.k_load_0.0000_*_0.0000 > /dev/null 2>&1
            ovito3.8 dump.crack > /dev/null 2>&1
            ;;
        path)
            echo $(pwd)
            ;;
        python_load)
            local orient_fname="./data/orientation.dat"
            local cij_fname=${lmp_temp_dir}/potential/$(pot_element_path ${pot_element})/${pot_id}/elastic_tensor/${atom_type}/${xtal_structure}/elastic_c_${atom_type}_${xtal_structure}_static_0K.mod
            python3 ${root_path}/python/k_load_crack.py -k ${kii0} ${ki0} ${kiii0} -d ${dkii0} ${dki0} ${dkiii0} -n ${k_load_step} -f ./dump.crack -p ./data/crack_pos.mod -e ${cij_fname} -o ${orient_fname}
            ;;
        debug)
            debug
            ;;
        *)
            send_msg "Unknown mode: $mode" 1
            return
            ;;
    esac
}
ovito_plot() {
    if [[ -z "$ovito_plot_files" ]]; then
        local ovito_plot_files=dump.crack
    fi
    # prepare_run_path 'figure' 'yes'
    module load ovito/3.0.0_dev329_gcc_dp_znver2
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m custom
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m custom -o _blunt
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m custom -o _enlarge
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m ir
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m enlarge -o _blunt
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m enlarge -o _onestep
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m enlarge -o _angle
    # ovitos ${root_path}/python/ovito_visualization.py -f ${ovito_plot_files[@]} -m enlarge -o _plane

    # rsync -avP ./figure/* ~/Documents/article/bcc_crack/data/ki_simulation/ovito_image/${atom_type}/${c2}_${c3}/
    # rsync -avP ${ovito_plot_files[@]}  ~/Documents/article/ir_fcc_crack/data/ovito_image/${c2}_${c3}/
    # rsync -avP ./figure/* ~/Documents/article/ir_fcc_crack/data/ovito_image/${c2}_${c3}/
}
nonlinear_displacement() {
    for fname in ${ovito_plot_files[@]}; do
        echo $fname
        python3 ${root_path}/python/nonlinear_displacement.py -f $fname
        python3 ${root_path}/python/atomman_interpolate_contour_plot.py -f $fname -p e_tr
        python3 ${root_path}/python/atomman_interpolate_contour_plot.py -f $fname -p e_tt
    done
}
compute_stroh() {
    local orient_fname="./data/orientation.dat"
    # local cij_fname=${lmp_temp_dir}/potential/$(pot_element_path ${pot_element})/${pot_id}/elastic_tensor/${atom_type}/${xtal_structure}/elastic_c_${atom_type}_${xtal_structure}_static_0K.mod
    local cij_fname=${lmp_temp_dir}/potential/$(pot_element_path ${pot_element})/${pot_id}/elastic_tensor/${atom_comp}/${xtal_structure}/elastic_c_${atom_comp}_${xtal_structure}_static_0K.mod
    local stroh_tensor_fname="./data/stroh_tensor.dat"
    python3 ${root_path}/python/stroh_tensor.py -o ${orient_fname} -e ${cij_fname} -s ${stroh_tensor_fname}
    rsync -avP ./data/orientation.dat ../../orientation.dat
}
prepare_run_file() {
	  local pot_prefix=${lmp_temp_dir}/potential/${pot_element_path}/${pot_id}
	  if [[ "${req_change_atom}" == true ]]; then
        # change_atom.mod: parsing change_atom; instruct lamps to change atom type;
        # atomic_composition.mod: store change_atom
		    change_atom ${atom_type} ${change_atom} ${pot_prefix} ${xtal_structure} # generate change_atom.mod, atomic_composition.mod
	  else
		    if [ -f atomic_composition.mod ]; then
			      rm -f atomic_composition.mod
		    fi
		    echo "# no change_atom defined" > change_atom.mod
        local lat_param_file=${pot_prefix}/unit_cell/${atom_comp}/${xtal_structure}/lat_param_${atom_comp}_${xtal_structure}_static_0K.mod
	      check_lat_param_file ${lat_param_file} atomic_composition.mod # define lat_param_file_exist in atomic_composition.mod
	  fi

	  prepare_run_path 'data' 'yes'
	  prepare_run_path 'dump_files' 'yes'

    # local lmp_infile=in.lammps
    rsync -a ${root_path}/${lmp_infile} ${lmp_infile}
    rsync -a ${root_path}/init.mod init.mod
    rsync -a ${root_path}/module/ module
    rsync -a ${root_path}/crack_tip/ crack_tip

    local lmp_infiles=()
    lmp_infiles+=(init.mod)

    local lmp_variables=()
    lmp_variables=($(get_lmp_replace_variable init.mod))
    change_lmp_infile # replace REPLACE_VARNAME with ${varname} in lmp_infiles
    change_lmp_basis ${atom_type} module/generate_structure.mod # modify create_atom command

    local jobname=${atom_type}_${calc_name}
    generate_jobfile $jobname job.${jobname}.sh ${lmp_infile}
}
calc() {
    local jobname=${atom_type}_${calc_name}
    local job_fname=job.${jobname}.sh
    ${job_submit_cmd} --exclude=gauss ${job_fname}
}
debug() {
    echo $(pwd)
    module load ${module}
    mpirun lmp -in ${lmp_infile} -screen /dev/null
    module purge
}
#====================================================================

usage() {
    send_msg "usage: run.sh prep|calc spec_file" ${default_info_msg_level}
    exit 0
}
process_cmd() {
    if [ $# -lt 2 ]; then
        usage
    elif [ $# -lt 3 ]; then
        mode=$1
        spec_file=$2
        send_msg "---------------------------------------" 1
        send_msg "    RUN MODE: ${mode} using ${spec_file}" ${default_info_msg_level}
        send_msg "---------------------------------------" 1
    else
        mode=$1
        spec_file=$2
        calc_mode_cl=$3
        send_msg "---------------------------------------" 1
        send_msg "    RUN MODE: ${mode} using ${spec_file}" ${default_info_msg_level}
        send_msg "---------------------------------------" 1
    fi
}

process_cmd $@
. $spec_file
if [[ ! -z $calc_mode_cl ]]; then
    calc_mode=$calc_mode_cl
fi


root_path=$PWD
pot_element_path=$(pot_element_path ${pot_element})

send_msg "run_path: ${run_path}" 1
prepare_run_path ${run_path} 'yes'

run_path=$(readlink -f ${run_path})
send_msg "changing directory to ${run_path}" 1
cd ${run_path}

run $root_path

send_msg "changing directory to $root_path" 0
cd $root_path
send_msg "---------------------------------------"
send_msg "    RUN $mode COMPLETED."
send_msg "---------------------------------------"
