#!/bin/bash

DIR=${1}
version=${2}

echo ""
echo "Decompressing."
echo ""
echo "Warning: the decompressed files will take up these amounts of disk space:"
for fn in `find ${DIR}/db/external_resources/Pfam_${version}/database_files/ -name '*.gz'`
do
        gzip -dc ${fn} | wc -c | ${DIR}/src/byte_to_human.sh | sed 's/^[ \t]*//;s/[ \t]*$//' | awk -v fn=`basename ${fn}` 'NF>1{print fn, $0}'
done

while true
do
echo ""
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

echo "Done"
echo ""
