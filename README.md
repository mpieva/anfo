Where to start?
===
Get the prerequisites, compile, create documentation, test.  Also see the Github Wiki.

Prerequisites are...

 - Not too ancient G++; 4.0 and later work fine, 3.x might cause some
   difficulties.  Also make, a shell, and so on.

 - protobuf.  Find it at http://code.google.com/p/protobuf

 - popt.  You might already have it, else ask your admin or grab it at
   http://rpm5.org/files/popt

 - Doxygen, if you want the documentation.  Can be found at http://doxygen.org

 - zlib from zlib.net

 - libbz2 from bzip2.org

 - Elk from http://sam.zoy.org/elk/

If you have those, call './configure', 'make', then 'make doc' to create
the documentation, then 'make install'.


How to use it?
===
To test, you need a genome, an index, and some input in FastA or FastQ format:

 - Locate a genome in FASTA format (or any other format that has a converter;
   .2bit is fine, my own .dna isn't).  Run fa2dna on it (--help tells
   you how); you can pipe the input if you want.

 - On your new .dna file, run file-info.  It should spit out the meta data for
   that genome.  If not, you're hosed...

 - Build an index for the .dna file using dnaindex (--help tells you
   how, the defaults are probably fine).  Use only the basename of the genome,
   it will figure the extension out by itself.

 - You can run file-info on that, too.  It will tell you the word size
   you just set.  If not, you're hosed again...

 - Write a sensible anfo.cfg using the examples in example/ and
   config.proto as guideline.

 - You can run index-test or anfo now.  Index-test just does lookups in
   the index, and optionally (on simulated data) checks if a seed in the
   correct region was found.  

 - anfo can be run on any FASTA/FASTQ file now.

 - anfo-tool can operate on the output files (but it's cumbersome to
   use)

 - the Elk binding can do what anfo and anfo-tool do, is much more
   flexible, but undocumented.


What's stable and what isn't?
===
Pretty much stabilized are...

 - Genome handling and creation (files sequence.h, sequence.cc, index.h, index.cc, fa2dna.cc)

   I store genomes with four bits per nucleotide, which sounded like a
   fscking brilliant idea, since it's only half as wasteful as the text
   form and still allows ambiguity codes.  The headaches came with the
   implementation of the auxilliary DnaP class...  Anyway, the
   nucleotides A,C,T,G map to bits 0,1,2,3.  That order (you did notice
   T coming before G, didn't you?) also sounded like a fscking brilliant
   idea, but in realilty it doesn't matter and I keep mixing it up.
   Chromosomes are split into contigs at long stretches of Ns, contigs
   are separated by single gaps and the first and last ones terminate in
   a gap at either side.  This means you can start anywhere and safely
   run forward or backward until you hit a gap.

 - Index handling and creation (files index.h, index.cc, dnaindex.cc)

   The index is quite simple: oligos are mapped to integer offsets, a
   first level array contains a pointer to a second array for each of
   the possible oligos, and the second level array contains one long
   list of positions where these oligos were found.  Only the forward
   strand is indexed, but lookup is of course done for both strands.
   (This primitive thing seems such a waste when the much better
   FM-index family is just out of reach... to bad, for now.)

 - Handling of reads (files sequence.h and sequence.cc)

   Just a simple structure for sequences with quality scores and a
   FastA/FastQ/FourQ reader.  There's no support for mate pairs, I
   haven't even decided what to do about them.

 - File formats (genome, index, config file, output)

   Genome and index files are documented somewhere in the code.
   Metadata, config, complex output, etc. is encapsulated in protobuf
   messages.  That way I don't need to mess with parsers and pretty
   printers and the messages are extensible.  The two .proto files are
   more or less finalized, but it's easy to change them without breaking
   stuff, which also means you should think before changing them.

 - The aligner (align.h, align.cc)

   Completely rewritten now, it is much faster than the previous version
   and practically finished.  The general structure and the specific
   alignment mechanisms are somewhat separated from each other, making
   extensions like a special 454 aligner possible at least in principle.


Copyright

Anfo is (C) 2009 by Udo Stenzel <udo_stenzel@eva.mpg.de>

ANFO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


