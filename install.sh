version=32

if [ "${version}" == "0" ]
then
	echo "Please set a Pfam version"
fi

find db/interaction_database/ -name '*.tgz' -execdir tar -xzvf '{}' \;

cd src/
python3 change_pfam_version.py ${version}

find db/external_resources/Pfam_${version}/database_files/ -name '*.gz' -execdir gunzip '{}' \;
