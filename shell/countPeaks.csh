#!/bin/csh -ef
set header =  "name,chr,start,end"
foreach b (`ls -1 *.bam`)
	set bn = $b:r:r
	set header = "$header","$bn"
end
echo $header
set samtools = /data/lowelab/software/samtools-1.9/samtools
foreach i (`cat MATAxLITC_merged.sorted.bed | tr '\t' '.'`)
	
	set name = $i:r:r:r:r:r:r:r:r:r:r:e
	set start = $i:r:r:r:r:r:r:r:r:r:r:r:r:e
	set end = $i:r:r:r:r:r:r:r:r:r:r:r:e
	set chr = $i:r:r:r:r:r:r:r:r:r:r:r:r:r
	set line = "$name","$chr","$start","$end"
	
	foreach bam (`ls -1 *.bam`)
		set n = $bam:r:r
		set count = `$samtools view $bam "$chr":"$start"-"$end" | wc -l`
		set line = "$line","$count"
	end
	echo $line
end
echo DONE	
