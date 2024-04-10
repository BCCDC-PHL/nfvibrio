if [ $# -ne 2 ]; then 
	echo "USAGE: bash make_sample_sheet.sh [DIRNAME] [OUTFILE]"
	exit 0 
fi


DIR=$1
OUTFILE=$2

echo "sample,fastq_1,fastq_2" > $OUTFILE

for f in $DIR/*R1*gz; do
	file1=`basename $f`
	name=`echo $file1 | cut -d\. -f1 | cut -d_ -f1` 
	file2=`echo "$file1" | sed 's/R1/R2/g'`
	
	echo $name 
	echo $file1

	echo "$name,$DIR/$file1,$DIR/$file2" >> $OUTFILE 
done 

