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
 break
       ;;
    *)
 echo "Invalid input..."
 ;;
esac
done
