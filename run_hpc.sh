#! /bin/zsh

if [ "$#" -ne 2 ] || [ "$1" = "--help" ]; then
  echo "Run the hybberish utility."
  echo "Arguments: <experiment name> <node nr 0..7>"
  exit 1
fi

cd /project_antwerp/hybberish/

EXP_NAME=$1
NODE_NR=$2

julia --project=. -e"using Pkg; Pkg.instantiate();"

for i in $(seq 0 3)
do
    FILE_NR=$((NODE_NR * 4 + i))
    INPUT_FILE="./${EXPERIMENT_NAME}/input/node${NODE_NR}/input_${FILE_NR}.csv"
    OUTPUT_FILE="./${EXPERIMENT_NAME}/out/output_${FILE_NR}.csv"
    echo "Processing '${INPUT_FILE}' into '${OUTPUT_FILE}'."
    julia --project=. ./examples/trucktrailer.jl 0.01 10 "${INPUT_FILE}" "${OUTPUT_FILE}" &
done
