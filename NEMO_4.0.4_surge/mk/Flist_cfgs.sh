#!/bin/bash

echo -e "\n  ¤ Reference configurations with default sub-components (can be updated by a new set)"
cat ${MAIN_DIR}/cfgs/ref_cfgs.txt \
| awk '{printf "%-20s", $1}{$1 = ""; printf "%s\n", $0}'

echo -e "\n  ¤ Demonstrations cases (see https://github.com/sflavoni/NEMO-test-cases for more)"
cat ${MAIN_DIR}/tests/demo_cfgs.txt \
| awk '{printf "%-20s", $1}{$1 = ""; printf "%s\n", $0}'

echo -e "\n  ¤ Full scripted remote configurations (CPP file + EXP00 inputs + MY_SRC + ...)"
cat ${MAIN_DIR}/tests/rmt_cfgs.txt

echo -e "\n  ¤ Available sub-components ('OCE' is mandatory in any set)"
ls  ${MAIN_DIR}/src | awk -F/ '{ print $NF }' | column

echo
