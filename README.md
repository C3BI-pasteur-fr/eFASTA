# eFASTA v.1.0.170208ac
eFASTA: extracting nucleotide segments from FASTA file


## INSTALLATION ##
Use gcj compiler to compile the program from the code source. To create the binary eFASTA from the terminal, go to the directory containing the source code EFASTA.java and type:<br /> 
  make


## USAGE ##
Launch eFASTA without arguments to read the documentation.<br /> 

USAGE: eFASTA \<options\><br />
<br />
<br />
  where options are:<br />
<br />
-fasta, -f \<infile\><br />
(mandatory) FASTA-formatted nucleotide sequence file name<br />
<br />
-coord, -c \<pattern:start-end\><br />
(mandatory) DNA segment to extract from [infile] defined by the first FASTA header containing [pattern] and the region comprised  between nucleotide  indexes [start] and [end] (both inclusive); if [start] > [end], the extracted nucleotide segment is reverse-complemented (mandatory)<br />
<br />
-outname, -o \<basename\><br />
extracted nucleotide  segment is written  into the FASTA-formatted output file [basename].fna (default: "seq")<br />
<br />
-cds<br />
to indicate  that the  extracted  DNA  segment is  a CDS, leading to  the writing  of its  amino  acid  translation (standard  genetic  code)  into a  second FASTA-formatted output file [basename].faa<br />
<br />
-fcds<br />
same as  option -cds  but to search for the full CDS that includes the specified region,  i.e. first occuring codon START before index [start]  and first occuring codon STOP after index [end]
