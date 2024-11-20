INPUT_PATH=$1
OUTPUT_PATH=$2
NUM='[+-]?[0-9]+(\.[0-9]+)?'
JINT="\($NUM\)\.\.\($NUM\)" # Julia interval
# Check script argument requirements.
if [ ! $# -eq 2 ]; then
    echo "Provided $# arguments. Expected 2 arguments: the input file containing clean flowpipes to process, and the output path."
    exit 1
fi

# Wrap the IntervalBox instances in a list
echo "get_clean_boxes() = [" > $OUTPUT_PATH

# Only keep lines starting in "var: " or empty lines
grep -E '(^\w+: )|(^\s+$)' $INPUT_PATH |
# Reduce "\n\r" to "\n". Empty lines reduce to "\n".
sed -E 's/\r$//g' |
# sed does not allow multi-line matching, so replace
# newlines by whitespaces. Empty lines reduce to two
# whitespaces.
tr '\n' ' ' |
# Replace each "var: [lo, hi]" by "(lo)..(hi),"
sed -E "s/\w+: \[($NUM), ($NUM)\]/\(\1\)..\(\3\),/g" |
# Prune the last julia interval, and re-insert newlines
sed -E "s/(($JINT, )+)$JINT,  /\1\n/g" |
# Wrap the intervals in a julia IntervalBox object
sed -E "s/(.*)/    IntervalBox(\1),/g" >> $OUTPUT_PATH

# Wrap the IntervalBox instances in a list
echo "]" >> $OUTPUT_PATH
