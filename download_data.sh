#!/bin/bash
set -euo pipefail

roman_to_int() {
  case "$1" in
    i|I) echo 1 ;;
    ii|II) echo 2 ;;
    iii|III) echo 3 ;;
    iv|IV) echo 4 ;;
    v|V) echo 5 ;;
    vi|VI) echo 6 ;;
    vii|VII) echo 7 ;;
    viii|VIII) echo 8 ;;
    ix|IX) echo 9 ;;
    x|X) echo 10 ;;
    *) echo 0 ;;
  esac
}

DATA_DIR="mosden/data/testing"

# ENDF --------------------------------------------------------------------
ENDF_VERSION="VII.1"
ALLOWED_VERSIONS=("VIII.0" "VII.1")

if [[ ! " ${ALLOWED_VERSIONS[*]} " =~ " ${ENDF_VERSION} " ]]; then
    echo "Error: Invalid ENDF version '${ENDF_VERSION}'"
    echo "Allowed versions: ${ALLOWED_VERSIONS[*]}"
    exit 1
fi

mkdir -p "$DATA_DIR"

LOWERCASE_VERSION=$(echo "${ENDF_VERSION//./}" | tr '[:upper:]' '[:lower:]')

ROMAN_PART="${LOWERCASE_VERSION//[0-9]/}"
DIGIT_PART="${LOWERCASE_VERSION//[^0-9]/}"
INTEGER_VALUE=$(roman_to_int "$ROMAN_PART")
LOWERCASE_VERSION="${INTEGER_VALUE}${DIGIT_PART}"

ENDF_DIR="${DATA_DIR}/endfb${LOWERCASE_VERSION}"
NFY_DIR="${ENDF_DIR}"
mkdir -p "$NFY_DIR"

if [[ "${ENDF_VERSION}" == "VII.1" ]]; then
  SEPARATOR="-"
elif [[ "${ENDF_VERSION}" == "VIII.0" ]]; then
  SEPARATOR="_"
fi

NFY_ZIP_NAME="ENDF-B-${ENDF_VERSION}${SEPARATOR}nfy.zip"
NFY_URL="https://www.nndc.bnl.gov/endf-b7.1/zips/${NFY_ZIP_NAME}"

echo "Downloading NFY data for ENDF/B-${ENDF_VERSION}..."
TEMP_ZIP="${NFY_DIR}/${NFY_ZIP_NAME}"
wget --show-progress -O "$TEMP_ZIP" "$NFY_URL"
echo "Extracting NFY data..."
unzip "$TEMP_ZIP" -d "$NFY_DIR"
rm "$TEMP_ZIP"
echo "NFY data handled"

# /ENDF --------------------------------------------------------------------

# IAEA --------------------------------------------------------------------
IAEA_DIR="${DATA_DIR}/iaea"
IAEA_FILE="$IAEA_DIR/eval.csv"
URL="https://www-nds.iaea.org/relnsd/delayedn/eval.csv"
mkdir -p "$IAEA_DIR"

echo "Downloading IAEA delayed neutron data..."
wget -q --show-progress -O "$IAEA_FILE" "$URL"
echo "Saved to $IAEA_FILE"

# /IAEA --------------------------------------------------------------------

# OpenMC --------------------------------------------------------------------