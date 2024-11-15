INPUT_PATH=$1
OUTPUT_PATH=$2
NUM='[+-]?[0-9]+(\.[0-9]+)?'
# Check script argument requirements.
if [ ! $# -eq 2 ]; then
    echo "Provided $# arguments. Expected 2 arguments: the input file containing clean flowpipes to process, and the output path."
    exit 1
fi

# Wrap the IntervalBox instances in a list
echo "get_clean_boxes() = [" > $OUTPUT_PATH

# Only keep lines starting in "x: " or "y: "
grep -E '(^x: )|(^y: )' $INPUT_PATH |
# Get rid of pesky newlines, which sed play well with
tr -d '\n\r' |
# Translate to julia IntervalBox objects
sed -E "s/x: \[($NUM), ($NUM)\]y: \[($NUM), ($NUM)\]/    IntervalBox\(\1..\3, \5..\7\),\n/g" >> $OUTPUT_PATH

# Wrap the IntervalBox instances in a list
echo "]" >> $OUTPUT_PATH
