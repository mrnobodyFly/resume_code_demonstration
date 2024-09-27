#!/usr/bin/env bash
# clean all redundant files output by latex
# usage: clean_latex_output [-g] [-c] [-f fname] [-d maxdepth]
function clean_latex_output(){
    local is_git_repo=false; local remove_command="rm"
    local require_confirmation=false
    local fname="*"
    local maxdepth=1

    local optstring='gcf:d:'
    local option
    local OPTIND OPTARG # very important for repetitive call of this function
    while getopts ${optstring} option; do
        case ${option} in
            g)
                is_git_repo=true
                remove_command='git rm'
                ;;
            c)
                require_confirmation=true
                ;;
            f)
                fname=$OPTARG
                ;;
            d)
                maxdepth=$OPTARG
                ;;
            ?)
            echo "Unknown option; Usage: clean_latex_output [-g] [-c] [-f fname]"
            exit 1
            ;;
        esac
    done

    local file_suffix=(
        aux bbl bcf blg dvi fdb_latexmk fls lof log lot out out.ps run.xml synctex.gz toc
    )

    local junk_files=()
    local suffix
    for suffix in ${file_suffix[@]}; do
        junk_files=(${junk_files[@]} $(find ./ -maxdepth ${maxdepth} -type f -and -name "${fname}.${suffix}"))
    done

    if ${require_confirmation}; then
        echo "Following files will be removed by ${remove_command}:"
        echo ${junk_files[@]}
        echo -n "Proceed?(y/n): "
        read answer
        if [[ ${answer} != 'y' ]]; then
           return
        fi
    fi

    ${remove_command} ${junk_files[@]}
}

clean_latex_output $@
