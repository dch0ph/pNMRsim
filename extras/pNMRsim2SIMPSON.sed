s/^proc  *{/proc main {} {\n  global par\n  set f [fsimpson]/
s/^pulseq  *{/proc pulseq {} {\n  global par/
s/\<verbose .*/verbose 1101/
s/\<addlb /faddlb \$f /
s/\([0-9]\)V[0-9]* /\1 /g
s/\([0-9]\)V[0-9]*$/\1/
s/\<ft/fft $f/
s/\<zerofill /fzerofill $f /
s/\<ft2d/fft $f 0 0 0 0/
s/\<save /fsave $f /
s/ -octant//
s/ -sphere//
s/ -hemisphere//
s/ -start//
s/ -middle//
s/ -end//
s/ -both//
/\#\#DELETENEXT/{N;d}
s/\#\#INSERT//
/\<truncate /d
/\<fit /d
/\<scalefirst /d
/\<histogram /d
/\<auto_opt /d
/\<tolerance /d
/\<precision /d
/\<eigenvalues /d
/\<sideband /d
/\<zerotolerance /d
/\<chebsyhev_iterations /d
/\<pseudohistogram /d
/\<cache_limit /d
/\<transients /d
/\<log_file /d
s/$(\(.*\))/$par(\1)/g
