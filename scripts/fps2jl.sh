INPUT_PATH=$1
OUTPUT_PATH=$2
REGEX_IN='x: (\[.+\])\ny: (\[.+\])'

# Check script argument requirements.
if [ ! $# -eq 2 ]; then
    echo "Provided $# arguments. Expected 2 arguments: the input file containing clean flowpipes to process, and the output path."
    exit 1
fi

# Wrap the IntervalBox instances in a list
echo "CLEAN_BOXES = [" > $OUTPUT_PATH

# Only keep lines starting in "x: " or "y: "
grep -E '(^x: )|(^y: )' $INPUT_PATH |
# Get rid of pesky newlines, which sed play well with
tr -d '\n\r' |
# Translate to julia IntervalBox objects
sed -E 's/x: (\[[0-9\., \+\-]+\])y: (\[[0-9\., \+\-]+\])/    IntervalBox\(\1, \2\),\n/g' >> $OUTPUT_PATH

# Wrap the IntervalBox instances in a list
echo "]" >> $OUTPUT_PATH
