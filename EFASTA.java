/*
  #################################################################################
  #                                                                               #
  # EFASTA: extracting nucleotide segments from FASTA file                        #
  # Copyright (C) 2017  Alexis Criscuolo                                          #
  #                                                                               #
  # EFASTA is free software;  you can redistribute it and/or modify it under      #
  # the terms  of the GNU  General Public License  as published by the Free       #
  # Software  Foundation;  either version 2  of the  License,  or  (at your       #
  # option) any later version.                                                    #
  #                                                                               #
  # EFASTA is distributed  in the hope that  it will be useful,  but WITHOUT      # 
  # ANY WARRANTY;  without even the  implied warranty of  MERCHANTABILITY or      #
  # FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for      #
  # more details.                                                                 #
  #                                                                               #
  # You should  have received  a copy  of the  GNU  General  Public  License      #
  # along   with  this   program;   if  not,  write  to  the  Free  Software      #
  # Foundation Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307            #
  #                                                                               # 
  # Contact:                                                                      #
  #  Alexis Criscuolo                                                             #
  #  HUB BioIT - C3BI                                                             #
  #  INSTITUT PASTEUR                                                             #
  #  25-28 rue du Dr Roux                                                         #
  #  75015 Paris (France)                                                         #
  #                                                                               #
  #  alexis.criscuolo@pasteur.fr                                                  #
  #                                                                               #
  #################################################################################
*/

import java.io.*;
import java.util.*;

public class EFASTA {

    // constants
    final static String VERSION = "1.0.170208ac";
    final static String[] TABLE_GENETIC_CODE_CODON =   {"AAA" , "AAC" , "AAG" , "AAT" ,  // A A -
							"ACA" , "ACC" , "ACG" , "ACT" ,  // . C -
							"AGA" , "AGC" , "AGG" , "AGT" ,  // . G -
							"ATA" , "ATC" , "ATG" , "ATT" ,  // . T -
							"CAA" , "CAC" , "CAG" , "CAT" ,  // C A -
							"CCA" , "CCC" , "CCG" , "CCT" ,  // . C -
							"CGA" , "CGC" , "CGG" , "CGT" ,  // . G -
							"CTA" , "CTC" , "CTG" , "CTT" ,  // . T -
							"GAA" , "GAC" , "GAG" , "GAT" ,  // G A -
							"GCA" , "GCC" , "GCG" , "GCT" ,  // . C -
							"GGA" , "GGC" , "GGG" , "GGT" ,  // . G -
							"GTA" , "GTC" , "GTG" , "GTT" ,  // . T -
							"TAA" , "TAC" , "TAG" , "TAT" ,  // T A -
							"TCA" , "TCC" , "TCG" , "TCT" ,  // . C -
							"TGA" , "TGC" , "TGG" , "TGT" ,  // . G -
							"TTA" , "TTC" , "TTG" , "TTT" }; // . T -
    final static char[] TABLE_GENETIC_CODE_AMINO_ACID = {'K'  ,  'N'  ,  'K'  ,  'N'  ,  // A A -
							 'T'  ,  'T'  ,  'T'  ,  'T'  ,  // . C -
							 'R'  ,  'S'  ,  'R'  ,  'S'  ,  // . G -
							 'I'  ,  'I'  ,  'M'  ,  'I'  ,  // . T -  ('M' = START)
							 'Q'  ,  'H'  ,  'Q'  ,  'H'  ,  // C A -
							 'P'  ,  'P'  ,  'P'  ,  'P'  ,  // . C -
							 'R'  ,  'R'  ,  'R'  ,  'R'  ,  // . G -
							 'L'  ,  'L'  ,  'L'  ,  'L'  ,  // . T -
							 'E'  ,  'D'  ,  'E'  ,  'D'  ,  // G A -
							 'A'  ,  'A'  ,  'A'  ,  'A'  ,  // . C -
							 'G'  ,  'G'  ,  'G'  ,  'G'  ,  // . G -
							 'V'  ,  'V'  ,  'V'  ,  'V'  ,  // . T -
							 '?'  ,  'Y'  ,  '?'  ,  'Y'  ,  // T A -  ('?' for STOP)
							 'S'  ,  'S'  ,  'S'  ,  'S'  ,  // . C -
							 '?'  ,  'C'  ,  'W'  ,  'C'  ,  // . G -  ('?' for STOP)
							 'L'  ,  'F'  ,  'L'  ,  'F'  }; // . T -

    // options
    static File dnaFile;
    static String fh;
    static int begin, end;
    static boolean isCDS, completeCDS;
    static String outmodel;
    
    // io
    static BufferedReader in;
    static BufferedWriter out;
    
    // variables
    static StringBuffer dnaSeq;
    static boolean perfect;
    static String dna, aa;
    
    // misc
    static int i, b, b1 ,b2, o, lgt;
    static boolean ok;
    static char c;
    static String line, codon;
    
    public static void main(String[] args) throws IOException {

	// ###########################
	// ##### documentation   #####
	// ###########################
	if ( args.length < 4 ) {
	    System.out.println("");
	    System.out.println(" eFASTA v." + VERSION);
	    System.out.println("");
	    System.out.println(" USAGE:   eFASTA <options>");
	    System.out.println("");
	    System.out.println("  where options are:");
	    System.out.println("");
	    System.out.println("   -fasta, -f <infile>             (mandatory) FASTA-formatted nucleotide sequence file name");
	    System.out.println("   -coord, -c <pattern:start-end>  (mandatory) DNA segment to extract  from [infile] defined");
	    System.out.println("                                   by the first FASTA  header  containing [pattern]  and the");
	    System.out.println("                                   region comprised  between nucleotide  indexes [start] and");
	    System.out.println("                                   [end] (both inclusive); if [start] > [end], the extracted");
	    System.out.println("                                   nucleotide segment is reverse-complemented (mandatory)");
	    System.out.println("   -outname, -o <basename>         extracted nucleotide  segment is written  into the FASTA-");
	    System.out.println("                                   formatted output file [basename].fna (default: \"seq\")");
	    System.out.println("   -cds                            to indicate  that the  extracted  DNA  segment is  a CDS,");
	    System.out.println("                                   leading to  the writing  of its  amino  acid  translation");
	    System.out.println("                                   (standard  genetic  code)  into a  second FASTA-formatted");
	    System.out.println("                                   output file [basename].faa");
	    System.out.println("   -fcds                           same as  option -cds  but to search for the full CDS that");
	    System.out.println("                                   includes the specified region,  i.e. first occuring codon");
	    System.out.println("                                   START before index [start]  and first occuring codon STOP");
	    System.out.println("                                   after index [end]");
	    System.out.println("");
	    System.out.println(" EXAMPLES:");
	    System.out.println("");
	    System.out.println("   eFASTA  -fasta Ecoli.fna  -coord NZ_AVRI01000134:1966-1328  -cds");
	    System.out.println("   eFASTA  -f Ecoli.fna  -c NZ_AVRI01000134:1966-1328  -fcds  -o gene ");
	    System.out.println("");
	    System.exit(0);
	}

	
	// ###########################
	// ##### reading options #####
	// ###########################
	dnaFile = new File("N.o.F.I.L.E");
	fh = "N.o.F.H";
	begin = 0;
	end = 1;
	isCDS = false;
	completeCDS = false;
	outmodel = "seq";
	o = -1;
	while ( ++o < args.length ) {
	    // dna fasta file
	    if ( args[o].equals("-fasta") || args[o].equals("-f") ) {
		dnaFile = new File( args[++o] );
		if ( ! dnaFile.exists() ) { System.out.println(" incorrect file name: " + dnaFile.toString() + " (option -f)"); System.exit(1); }
		continue;
	    }
	    // coordinates of the dna segement to output
	    if ( args[o].equals("-coord") || args[o].equals("-c") ) {
		line = args[++o];
		if ( line.endsWith(":") || line.endsWith("-") ) { System.out.println(" incorrect pattern: " + line + " (option -c)"); System.exit(1); }
		if ( (b=line.indexOf(":")) == -1 ) { System.out.println(" missing \":\": " + line + " (option -c)"); System.exit(1); }
		fh = line.substring(0 , b); line = line.substring(++b);
		if ( (b=line.indexOf("-")) == -1 ) { System.out.println(" missing \"-\": " + line + " (option -c)"); System.exit(1); }
		try { begin = Integer.parseInt( line.substring(0 , b) ); end = Integer.parseInt( line.substring(++b) ); }
		catch ( NumberFormatException e ) { System.out.println(" incorrect pattern: " + line + " (option -c)"); System.exit(1); }
		continue;
	    }
	    // is dna fragment a cds?
	    if ( args[o].equals("-cds") ) { isCDS = true; continue; }
	    // complete cds required?
	    if ( args[o].equals("-fcds") ) { completeCDS = true; continue; }
	    // outfile pattern
	    if ( args[o].equals("-outname") || args[o].equals("-o") ) { outmodel = args[++o]; continue; }
	} 
	if ( dnaFile.toString().equals("N.o.F.I.L.E") ) { System.out.println(" missing file name (option -f)"); System.exit(1); }


	// ##################################
	// ##### reading dna fasta file #####
	// ##################################
	in = new BufferedReader(new FileReader(dnaFile));
	line = "";
	while ( (! line.startsWith(">")) || (line.indexOf(fh) == -1) ) 
	    try { line = in.readLine().trim(); } catch ( NullPointerException e ) { in.close(); System.out.println(fh + " not found inside " + dnaFile.toString()); System.exit(1); }
	fh = line;
 	dnaSeq = new StringBuffer(""); lgt = ( begin < end ) ? end : begin; lgt += ( completeCDS ) ? 10000 : 0;
	while ( dnaSeq.length() < lgt ) {
	    try { if ( (line=in.readLine().trim()).startsWith(">") ) { in.close(); break; } } catch ( NullPointerException e ) { in.close(); break; }
	    dnaSeq = dnaSeq.append(line.toUpperCase());
	}
	//System.out.println(fh); System.out.println("   begin=" + begin + "   end=" + end);


	// ##################################
	// ##### getting DNA segment    #####
	// ##################################
	// no CDS
	if ( (! isCDS) && (! completeCDS) ) {
	    if ( begin < end ) { dna = dnaSeq.substring(--begin , end); ++begin; fh = fh + "::" + begin + "-" + end; }
	    else { dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; fh = fh + "::" + begin + "-" + end; }
	    out = new BufferedWriter(new FileWriter(new File(outmodel + ".fna"))); out.write(fh); out.newLine(); out.write(dna); out.newLine(); out.close();
	    System.exit(0);
	}
	// CDS expected
	if ( begin < end ) { 
	    dna = dnaSeq.substring(--begin , end); ++begin; aa = dna2aa(dna); 
	    // complete CDS required
	    if ( completeCDS ) {
		// looking for start codon
		while ( (begin > 3) && ((c=aa.charAt(0)) != 'M') && (c != '?') && (! (codon=dna.substring(0,3)).equals("GTG")) && (! codon.equals("TTG")) ) { 
		    begin -= 3; dna = dnaSeq.substring(--begin, end); ++begin; aa = dna2aa(dna); 
		}
		if ( aa.charAt(0) == '?' ) { begin += 3; dna = dnaSeq.substring(--begin, end); ++begin; aa = dna2aa(dna); }
		// looking for end codon
		lgt = dnaSeq.length(); lgt -= 2; while ( (end < lgt) && (lastChar(aa) != '?') ) { end += 3; dna = dnaSeq.substring(--begin, end); ++begin; aa = dna2aa(dna); } lgt += 2;
		if ( lastChar(aa) == '?' ) { end -= 3; dna = dnaSeq.substring(--begin, end); ++begin; aa = dna2aa(dna); }
	    }
	    fh = fh + "::" + begin + "-" + end;
	    out = new BufferedWriter(new FileWriter(new File(outmodel + ".fna"))); out.write(fh); out.newLine(); out.write(dna); out.newLine(); out.close();
	    out = new BufferedWriter(new FileWriter(new File(outmodel + ".faa"))); out.write(fh); out.newLine(); out.write(aa);  out.newLine(); out.close();
	    System.exit(0);
	}
	// CDS expected on reverse complement
	dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; aa = dna2aa(dna);
	
	// complete CDS required
	if ( completeCDS ) {
	    // looking for start codon
	    lgt = dnaSeq.length(); lgt -= 2; 
	    while ( (begin < lgt) && ((c=aa.charAt(0)) != 'M') && (c != '?') && (! (codon=dna.substring(0,3)).equals("GTG")) && (! codon.equals("TTG")) ) { 
		begin += 3; dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; aa = dna2aa(dna); 
	    }
	    if ( aa.charAt(0) == '?' ) { begin -= 3; dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; aa = dna2aa(dna); }
	    // looking for end codon
	    while ( (end > 3) && (lastChar(aa) != '?') ) { end -= 3; dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; aa = dna2aa(dna); }
	    if ( lastChar(aa) == '?' ) { end += 3; dna = dnaReverseComplement( dnaSeq.substring(--end, begin) ); ++end; aa = dna2aa(dna); }
	}
	fh = fh + "::" + begin + "-" + end;
	out = new BufferedWriter(new FileWriter(new File(outmodel + ".fna"))); out.write(fh); out.newLine(); out.write(dna); out.newLine(); out.close();
	out = new BufferedWriter(new FileWriter(new File(outmodel + ".faa"))); out.write(fh); out.newLine(); out.write(aa);  out.newLine(); out.close();
    }




    // ##### a method to translate DNA sequence String to amino acid String #####
    final static String dna2aa ( String sequence ) {
	int _b1, _b2, _m = sequence.length() / 3;
	StringBuffer _sb = new StringBuffer(sequence.substring(0 , _m));
	String _codon;
	_b1 = -1;
	while ( ++_b1 < _m ) {
	    _codon = sequence.substring(3*_b1 , 3*(++_b1));
	    _b2 = Arrays.binarySearch(TABLE_GENETIC_CODE_CODON , _codon);
	    if ( _b2 >= 0 ) _sb.setCharAt(--_b1 , TABLE_GENETIC_CODE_AMINO_ACID[_b2]);
	    else  _sb.setCharAt(--_b1 , 'X');
	}
	return _sb.toString().trim();
    }



    // ##### a method to translate DNA sequence String to its reverse complement #####
    final static String dnaReverseComplement ( String sequence ) {
	int __b1 = sequence.length(); int __b2 = -1;
	StringBuffer __sb = new StringBuffer(sequence);
	while ( --__b1 >= 0 ) 
	    switch ( sequence.charAt(__b1) ) {
	    case 'A': __sb.setCharAt(++__b2 , 'T'); break;
	    case 'C': __sb.setCharAt(++__b2 , 'G'); break;
	    case 'G': __sb.setCharAt(++__b2 , 'C'); break;
	    case 'T': __sb.setCharAt(++__b2 , 'A'); break;
	    default:  __sb.setCharAt(++__b2 , '?'); break;
	    }
	return __sb.toString().trim();
    }

    
    static char lastChar( String sequence ) {
	return sequence.charAt(sequence.length() - 1);
    }

}

