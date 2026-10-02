#!/usr/bin/env bash
set -e

usage() {
  cat <<EOF
usage: $0 options

Install or update an environment

OPTIONS:
   -h     Show this message
   -f     Environment YAML file
   -c     Command to install conda packages, either 'mamba' or 'conda' (default: conda)
   -q     Quiet installation
EOF
}

CONDA_CMD="conda"
QUIET="" # "-q"
ENVS_DIR=$(dirname $0)
EXECUTE=true

while getopts "hf:c:nq" OPTION; do
    case $OPTION in
        h) usage; exit 0;;
        f) FILE=$OPTARG;;
        c) CONDA_CMD=$OPTARG;;
        n) EXECUTE=false;;
        q) QUIET="-q";;
        ?) usage; exit;;
        :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
    esac
done

if [[ $FILE == "" ]]; then
    echo "Error: Must specify environment YAML file with -f" >&2
    usage

    exit 1
fi

ENV=$(basename $FILE | cut -d. -f1)
ALL_ENVS=$($CONDA_CMD env list)

if [[ $ALL_ENVS == *"$ENV"* ]]; then
    operation="update"
    if [[ $CONDA_CMD == "mamba" ]]; then
      operation="$operation -y"
    fi
else
    operation="create -y"
fi
echo "$operation $ENV from $FILE..."
if $EXECUTE; then
    $CONDA_CMD env $operation $QUIET --file $FILE
fi

# Dependencies that are not available from conda channels are installed by a post-deploy
# script next to the YAML file, which Snakemake also runs for environments created from YAML
POST_DEPLOY="${FILE%.yaml}.post-deploy.sh"
if [[ -f $POST_DEPLOY ]]; then
    echo "run $POST_DEPLOY in $ENV..."
    if $EXECUTE; then
        $CONDA_CMD run --live-stream -n $ENV bash $POST_DEPLOY
    fi
fi
