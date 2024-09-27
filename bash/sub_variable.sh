#!/usr/bin/env bash
#------------------------------------------------
# usage substitute_variables [-q] fname [open] [close]
# substitute variables in fname using values of current bash variable
# if variable is not defined, comment the found variable line
#------------------------------------------------
function substitute_variables() {
    local warning=true
    local OPTSTRING=':q'
    local option
    local OPTIND OPTARG # very important for repetitive call of this function
    while getopts ${OPTSTRING} option; do
        case ${option} in
            q)
                warning=false
                ;;
            ?)
            echo "Unknown option; Usage: substitute_variables [-q] fname [open] [close]"
            exit
            ;;
        esac
    done
    shift $((OPTIND - 1))
    if [ $# -eq 1 ]; then
        local fname=$1
        local open='<|'
        local close='|>'
    elif [ $# -eq 3 ]; then
        local fname=$1
        local open=$2
        local close=$3
    else
        echo "Wrong number of arguments; Usage: substitute_variables [-q] fname [open] [close]"
        exit
    fi

    local variable_keys=($(grep -o "${open}.*${close}"  $fname | sort -u))
    for vk in ${variable_keys[@]}; do
        local vn=${vk#$open}; vn=${vn%$close}
        if [ -z ${!vn+x} ];  then # indirect expansion
            if ${warning}; then
                echo "variable ${vn} is not defined."
            fi
            sed -i "s@^variable.*${open}${vn}${close}@# &@g" $fname
        else
            sed -i "s@${open}${vn}${close}@${!vn}@g" $fname
        fi
    done
}
