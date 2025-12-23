#!/bin/bash
case $1 in
    # beginning with -, directly pass to native cut
    -*) /usr/bin/cut "$@";;

    # specify range, complement the -f and pass to native cut
    *-*)/usr/bin/cut -f"$@";;

    # other case of beginning with number, use awk
    [0-9]*)
        if [ $# -eq 1 ]; then
	    awk '{print $'${1//,/,$}'}' # ex. 1,3 -> $1,$3
	else
	    awk -F "$2" '{print $'${1//,/,$}'}'
	fi;;

    # beginning with non number, it's error.
    *)  echo "Usage: "$(basename $0)" N [SEP]"; exit 1;;
esac
