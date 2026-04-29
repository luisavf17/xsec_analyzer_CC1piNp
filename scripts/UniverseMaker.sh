# Number of expected command-line arguments
num_expected=3

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: UniverseMaker.sh FPM_Config SEL_CONFIG OUTPUT_ROOT_FILE"
  exit 1
fi

FPM_Config=$1
BIN_CONFIG=$2
OUTPUT_ROOT_FILE=$3

if [ ! -f "$FPM_Config" ]; then
  echo "FPM_Config \"${FPM_Config}\" not found"
  exit 1
fi

if [ ! -f "$BIN_CONFIG" ]; then
  echo "BIN_CONFIG \"${BIN_CONFIG}\" not found"
  exit 1
fi

univmake ${FPM_Config} ${BIN_CONFIG} ${OUTPUT_ROOT_FILE} ${FPM_Config}
