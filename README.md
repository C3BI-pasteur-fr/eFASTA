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
   -fasta, -f \<infile\>             (mandatory) FASTA-formatted nucleotide sequence file name<br />
   -coord, -c \<pattern:start-end\>  (mandatory) DNA segment to extract  from [infile] defined<br />
                                   by the first FASTA  header  containing [pattern]  and the<br />
                                   region comprised  between nucleotide  indexes [start] and<br />
                                   [end] (both inclusive); if [start] > [end], the extracted<br />
                                   nucleotide segment is reverse-complemented (mandatory)<br />
   -outname, -o \<basename\>         extracted nucleotide  segment is written  into the FASTA-<br />
                                   formatted output file [basename].fna (default: "seq")<br />
   -cds                            to indicate  that the  extracted  DNA  segment is  a CDS,<br />
                                   leading to  the writing  of its  amino  acid  translation<br />
                                   (standard  genetic  code)  into a  second FASTA-formatted<br />
                                   output file [basename].faa<br />
   -fcds                           same as  option -cds  but to search for the full CDS that<br />
                                   includes the specified region,  i.e. first occuring codon<br />
                                   START before index [start]  and first occuring codon STOP<br />
                                   after index [end]<br />
