# eFASTA v.1.0.170208ac
eFASTA: extracting nucleotide segments from FASTA file


## INSTALLATION ##
Use gcj compiler to compile the program from the code source. To create the binary eFASTA from the terminal, go to the directory containing the source code EFASTA.java and type:
  make


## USAGE ##
Launch eFASTA without arguments to read the following documentation:

 USAGE:   eFASTA <options>

  where options are:

   -fasta, -f <infile>             (mandatory) FASTA-formatted nucleotide sequence file name
   -coord, -c <pattern:start-end>  (mandatory) DNA segment to extract  from [infile] defined
                                   by the first FASTA  header  containing [pattern]  and the
                                   region comprised  between nucleotide  indexes [start] and
                                   [end] (both inclusive); if [start] > [end], the extracted
                                   nucleotide segment is reverse-complemented (mandatory)
   -outname, -o <basename>         extracted nucleotide  segment is written  into the FASTA-
                                   formatted output file [basename].fna (default: "seq")
   -cds                            to indicate  that the  extracted  DNA  segment is  a CDS,
                                   leading to  the writing  of its  amino  acid  translation
                                   (standard  genetic  code)  into a  second FASTA-formatted
                                   output file [basename].faa
   -fcds                           same as  option -cds  but to search for the full CDS that
                                   includes the specified region,  i.e. first occuring codon
                                   START before index [start]  and first occuring codon STOP
                                   after index [end]
