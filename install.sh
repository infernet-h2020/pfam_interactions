#!/bin/bash

# Parameters
version=32


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
find ${DIR}/db/interaction_database/ -name '*.tgz' -execdir tar -xzvf '{}' \;


# Download Pfam version and unzip files
python3 ${DIR}/src/change_pfam_version.py ${version}
find ${DIR}/db/external_resources/Pfam_${version}/database_files/ -name '*.gz' -execdir gunzip '{}' \;
