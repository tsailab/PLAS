#!/bin/bash
a="ttslbthamd"
keep="no"
dir="./"
lastArg="${!#}"
if [[ -f $lastArg && "$lastArg" != "./graphbylength.sh" ]]; then
        echo "Input file: $lastArg"
else
        echo 'No input file! Use "./graphbylength.sh -h" for help'
        exit
fi

while getopts ":c:n:d:kh" opt; do
        case $opt in
        h)
                echo "#######################################"
                echo "#        Read instance counter        #"
                echo "#######################################"
                echo "# FASTA format input required         #"
                echo "# ./graphbylength.sh [OPTION] [FASTA] #"
                echo "# Options:                            #"
                echo "# -c: Cutoff. UPPER LIMIT value for   #"
                echo "#     density plot consideration      #"
                echo "# -k: Keep intermediate files         #"
                echo "#     for debugging purposes          #"
                echo "# -n: Change intermediate file name   #"
                echo "# -h: Help. Display this message      #"
                echo "#######################################"
                exit 0
                ;;
        c)
                echo "Setting cutoff to $OPTARG" >&1
                cutoff=$OPTARG
                ;;
		d)
				echo "Setting destination directory to $OPTARG" >&1
				dir=$OPTARG
				;;
        k)
                echo "Keeping intermediate files" >&1
                keep="yes"
                ;;
        n)
                echo "Changing intermediate filename to $OPTARG" >&1
                a=$OPTARG
                ;;
        \?)
                echo 'Invalid argument: -$OPTARG. Use "./graphbylength.sh -h" for argument listing' >&2
                exit 1
                ;;
        :)
                echo "Cutoff flag -$OPTARG requires a value." >&2
                exit 1
                ;;
        esac
done
mkdir $dir

if [ !-d $dir ]; then
	echo "Failed to create directory!"
	exit 0
fi

###Linearize file with AWK, sort by length, convert back to fasta format
echo "Linearizing file and sorting..."
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'  $lastArg  |\
awk -F '\t' '{printf("%d\t%s\n",length($2),$0);}' |\
sort -k1,1n | cut -f 2- | tr "\t" "\n" > "$dir"/Sorted_"$lastArg"

###Get length of each line, remove every other line (header lines), sort and count
echo "Getting line length..."
awk '{ print length($0); }' "$dir"/Sorted_"$lastArg"  > "$dir/$a".out

echo "Removing header lines..."
sed -n '1~2!p' < "$dir/$a".out > "$dir/$a".out.cut

echo "Sorting and counting by length..."
sort -n "$dir/$a".out.cut | uniq -c > "$dir/$a".out.cut.counted

###Re-sort, swap columns, change to csv
echo "Re-sorting file..."
sort -n "$dir/$a".out.cut.counted > "$dir/$a".out.cut.counted.sorted

echo "Swapping columns..."
awk '{ t = $1; $1 = $2; $2 = t; print; }' "$dir/$a".out.cut.counted.sorted > "$dir/$a".out.cut.counted.sorted.swapped

echo "Sorting one last time..."
sort -n "$dir/$a".out.cut.counted.sorted.swapped > "$dir/$a".out.cut.counted.sorted.swapped.sorted

echo "Adding headers..."
sed '1i length count' < "$dir/$a".out.cut.counted.sorted.swapped.sorted > "$dir/$a".out.cut.counted.sorted.swapped.sorted.header

echo "Converting to csv..."
tr ' ' ',' < "$dir/$a".out.cut.counted.sorted.swapped.sorted.header > Final.csv

if [ $keep == "no" ]; then
        echo "Removing intermediate files..."
        rm Sorted_"$lastArg"
        rm "$dir/$a".out*
fi
echo "All done!"

###Use R script to create density plot of length values
module load R
echo "Creating R graph..."
Rscript creategraph.R
echo 'Created plot in "Rplots.pdf"!'
