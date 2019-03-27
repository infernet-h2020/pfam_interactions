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


echo "Upgrading pip..."
pip install --upgrade pip 2&>1 > /dev/null
echo ""

echo "Checking python3 packages..."
pip install -r python3_requirements.txt 2&>1 > /dev/null

# Download Pfam version and unzip files
python3 ${DIR}/src/change_pfam_version.py ${version}

echo "Decompressing." 
echo "Warning: the decompressed files will take up these amounts of disk space:"
for fn in `find ${DIR}/db/external_resources/Pfam_${version}/database_files/ -name '*.gz'`
do
        gzip -dc ${fn} | wc -c | ${DIR}/src/byte_to_human.sh | sed 's/^[ \t]*//;s/[ \t]*$//' | awk -v fn=`basename ${fn}` 'NF>1{print fn, $0}'
done

while true
do
read -r -p "Proceed with decompression? [Y/n] " input

case $input in
    [yY][eE][sS]|[yY]|"")
 echo "Decompressing..."
 break
 ;;
    [nN][oO]|[nN])
 echo "Installation was interrupted."
 exit 1
       ;;
    *)
 echo "Invalid input..."
 ;;
esac
done

find ${DIR}/db/external_resources/Pfam_${version}/database_files/ -name '*.gz' -execdir gunzip '{}' \;

echo "Installation completed!"
