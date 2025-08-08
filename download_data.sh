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

DATA_DIR="mosden/data/unprocessed"

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
IAEA_URL="https://www-nds.iaea.org/relnsd/delayedn/eval.csv"
mkdir -p "$IAEA_DIR"

echo "Downloading IAEA delayed neutron data..."
wget -q --show-progress -O "$IAEA_FILE" "$IAEA_URL"
echo "Saved to $IAEA_FILE"

# /IAEA --------------------------------------------------------------------

# OpenMC --------------------------------------------------------------------
OPENMC_DIR="${ENDF_DIR}/omcchain/"
mkdir -p "$OPENMC_DIR"
wget -q --show-progress -O "${OPENMC_DIR}chain_casl_pwr.xml" "https://anl.box.com/shared/static/3nvnasacm2b56716oh5hyndxdyauh5gs.xml"
wget -q --show-progress -O "${OPENMC_DIR}chain_casl_sfr.xml" "https://anl.box.com/shared/static/9fqbq87j0tx4m6vfl06pl4ccc0hwamg9.xml"
if [[ "${ENDF_VERSION}" == "VII.1" ]]; then
  wget -q --show-progress -O "${OPENMC_DIR}chain_endfb71_pwr.xml" "https://anl.box.com/shared/static/os1u896bwsbopurpgas72bi6aij2zzdc.xml"
  wget -q --show-progress -O "${OPENMC_DIR}chain_endfb71_sfr.xml" "https://anl.box.com/shared/static/9058zje1gm0ekd93hja542su50pccvj0.xml"
elif [[ "${ENDF_VERSION}" == "VIII.0" ]]; then
  wget -q --show-progress -O "${OPENMC_DIR}chain_endfb71_pwr.xml" "https://anl.box.com/shared/static/nyezmyuofd4eqt6wzd626lqth7wvpprr.xml"
  wget -q --show-progress -O "${OPENMC_DIR}chain_endfb71_sfr.xml" "https://anl.box.com/shared/static/x3kp739hr5upmeqpbwx9zk9ep04fnmtg.xml"
fi
