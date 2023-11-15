#!/bin/bash

# Parameters
version=35


# Main directory path
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"


# Parameter checks
if [ "${version}" == "0" ]
then
	echo "Please set a Pfam version"
fi


# Untar interaction benchmark
echo "Decompressing Pfam interaction database..."
find ${DIR}/db/interaction_database/ -name '*.tgz' -execdir tar -xzvf '{}' \; >/dev/null 2>&1
echo ""

# Create executable
mkdir bin/
python3_path=`which python3`
if [[ "${python3_path}" != "" ]]
then
	echo "#!${python3_path}" | cat - src/pfam2pdb.py > bin/pfam2pdb
	chmod 777 bin/pfam2pdb
fi

echo "Upgrading pip..."
pip install --upgrade pip >/dev/null 2>&1
echo ""

echo "Checking python3 packages..."
pip install -r python3_requirements.txt >/dev/null 2>&1

# Download Pfam version and unzip files
python3 ${DIR}/src/change_pfam_version.py ${version}

echo "Installation complete!"
