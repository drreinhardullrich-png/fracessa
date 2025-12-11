#!/bin/bash
# Run all matrices (except 34) with logging enabled

EXECUTABLE="./fracessa/build/fracessa"
OUTPUT_FILE="matrix_results.txt"

# Matrix data - excluding ID 34
# For circular symmetric (is_cs=true): use compact format (n/2 values) from "matrix" field
# For symmetric (is_cs=false): use upper triangular format (n*(n+1)/2 values) from "matrix" field
declare -A MATRICES
MATRICES[1]="2#0,1,0"  # symmetric, upper triangular
MATRICES[2]="2#3,3/2,4"  # symmetric, upper triangular
MATRICES[3]="3#4,13/2,1/2,5,11/2,3"  # symmetric, upper triangular
MATRICES[4]="3#1,2,1,1,1,1"  # symmetric, upper triangular
MATRICES[5]="3#0,0,1/2,0,1/2,0"  # symmetric, upper triangular
MATRICES[6]="3#0,2,-1,0,3,1/2"  # symmetric, upper triangular
MATRICES[7]="4#-1,1,0,0,-1,0,0,0,-1,0"  # symmetric, upper triangular
MATRICES[8]="5#1,3"  # circular symmetric
MATRICES[9]="5#1,0,2,2,2,1,2,2,2,1,0,0,1,0,1"  # symmetric, upper triangular
MATRICES[10]="6#1,-1,1"
MATRICES[11]="7#1,5,8"
MATRICES[12]="7#8,13,2"
MATRICES[13]="8#5,5,13,11"
MATRICES[14]="9#1,5,10,13"
MATRICES[15]="10#13,19,25,7,5"
MATRICES[16]="11#1,4,6,6,5"
MATRICES[17]="11#15,1/100,14,11/10,38/5"
MATRICES[18]="12#533,118,326,357,118,477"
MATRICES[19]="12#287,50,200,150,113,200"
MATRICES[20]="12#292,-1,293,-1,293,-1"
MATRICES[21]="12#16029525,4007380,9088537,12022143,2147548,16029525"
MATRICES[22]="13#7,18,18,10,7,10"
MATRICES[23]="14#4,3,2,2,3,4,-1"
MATRICES[24]="15#9,9,5,9,3,5,9"
MATRICES[25]="16#7,5,4,4,4,5,7,-1"
MATRICES[26]="17#1627,1039,1527,1075,1527,1039,1627,648"
MATRICES[27]="18#8,7,5,5,5,5,7,8,-1"
MATRICES[28]="19#13,31,31,27,31,27,13,13,27"
MATRICES[29]="19#1,2,2,2,2,2,1,1,2"
MATRICES[30]="20#13,13,13,8,8,8,13,13,13,-1"
MATRICES[31]="21#15,15,7,15,15,7,7,15,7,15"
MATRICES[32]="22#2,2,4,4,5,5,4,4,2,2,-1"
MATRICES[33]="23#27478,22664,10976,25676,18552,18552,25676,10976,22664,27478,17939"
MATRICES[35]="18#581,294,539,448,702,431,676,568,431"

# Expected ESS counts
declare -A EXPECTED_ESS
EXPECTED_ESS[1]=1
EXPECTED_ESS[2]=2
EXPECTED_ESS[3]=1
EXPECTED_ESS[4]=1
EXPECTED_ESS[5]=0
EXPECTED_ESS[6]=1
EXPECTED_ESS[7]=0
EXPECTED_ESS[8]=5
EXPECTED_ESS[9]=6
EXPECTED_ESS[10]=9
EXPECTED_ESS[11]=14
EXPECTED_ESS[12]=14
EXPECTED_ESS[13]=20
EXPECTED_ESS[14]=30
EXPECTED_ESS[15]=50
EXPECTED_ESS[16]=66
EXPECTED_ESS[17]=55
EXPECTED_ESS[18]=105
EXPECTED_ESS[19]=12
EXPECTED_ESS[20]=24
EXPECTED_ESS[21]=105
EXPECTED_ESS[22]=143
EXPECTED_ESS[23]=224
EXPECTED_ESS[24]=360
EXPECTED_ESS[25]=512
EXPECTED_ESS[26]=493
EXPECTED_ESS[27]=1152
EXPECTED_ESS[28]=1444
EXPECTED_ESS[29]=19
EXPECTED_ESS[30]=2560
EXPECTED_ESS[31]=4410
EXPECTED_ESS[32]=5632
EXPECTED_ESS[33]=2507
EXPECTED_ESS[35]=258

echo "Matrix ID | Expected ESS | Found ESS | Status" > "$OUTPUT_FILE"
echo "----------|--------------|-----------|--------" >> "$OUTPUT_FILE"

# Run each matrix
for ID in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 35; do
    MATRIX="${MATRICES[$ID]}"
    EXPECTED="${EXPECTED_ESS[$ID]}"
    
    echo "Running matrix ID $ID..."
    RESULT=$("$EXECUTABLE" --log -m "$ID" "$MATRIX" 2>&1 | tail -1)
    
    if [ -z "$RESULT" ]; then
        RESULT="ERROR"
    fi
    
    if [ "$RESULT" = "$EXPECTED" ]; then
        STATUS="✓ MATCH"
    else
        STATUS="✗ MISMATCH"
    fi
    
    echo "$ID | $EXPECTED | $RESULT | $STATUS" >> "$OUTPUT_FILE"
    echo "  Expected: $EXPECTED, Found: $RESULT"
done

echo ""
echo "Results saved to $OUTPUT_FILE"

