#!/bin/bash
# Generate kcachegrind (callgrind) output for matrix IDs 27-33

EXECUTABLE="./fracessa/build/fracessa"
OUTPUT_DIR="./callgrind_reports"

# Matrix data for IDs 27-33
declare -A MATRICES
MATRICES[27]="18#8,7,5,5,5,5,7,8,-1"
MATRICES[28]="19#13,31,31,27,31,27,13,13,27"
MATRICES[29]="19#1,2,2,2,2,2,1,1,2"
MATRICES[30]="20#13,13,13,8,8,8,13,13,13,-1"
MATRICES[31]="21#15,15,7,15,15,7,7,15,7,15"
MATRICES[32]="22#2,2,4,4,5,5,4,4,2,2,-1"
MATRICES[33]="23#27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939"

# Check if valgrind is available
if ! command -v valgrind &> /dev/null; then
    echo "Error: valgrind is not installed. Please install it with: sudo apt-get install valgrind"
    exit 1
fi

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Generating kcachegrind (callgrind) output for matrix IDs 27-33..."
echo "Output directory: $OUTPUT_DIR"
echo ""

# Process each matrix ID
for ID in {27..33}; do
    MATRIX="${MATRICES[$ID]}"
    OUTPUT_FILE="$OUTPUT_DIR/callgrind.out.$ID"
    
    echo "Processing matrix ID $ID..."
    echo "  Matrix: $MATRIX"
    echo "  Output: $OUTPUT_FILE"
    
    # Run valgrind with callgrind
    valgrind --tool=callgrind \
             --callgrind-out-file="$OUTPUT_FILE" \
             --dump-instr=yes \
             --collect-jumps=yes \
             "$EXECUTABLE" -m "$ID" "$MATRIX" > /dev/null 2>&1
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Successfully generated $OUTPUT_FILE"
    else
        echo "  ✗ Failed to generate profile for ID $ID"
    fi
    echo ""
done

echo "Done! Callgrind output files are in $OUTPUT_DIR"
echo "Open them with: kcachegrind $OUTPUT_DIR/callgrind.out.XX"

