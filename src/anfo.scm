;    Copyright 2009 Udo Stenzel
;    This file is part of ANFO
;
;    ANFO is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    Anfo is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with Anfo.  If not, see <http://www.gnu.org/licenses/>.

; Guile extension for Anfo
;
; This is supposed to provide high-level access to all sorts of
; post-processing facilities that used to live in anfo-tool, maybe even
; the complete Anfo functionality.  The hope is that a scripting
; language provides a powerful and more regular interface than an ad-hoc
; command line ever could.
;
; To use, start up guile in any suitable way (e.g. from the command line
; or as script interpreter; see the Guile manual) and issue
;
; (use-modules (anfo))

(define-module (anfo))
(use-modules (ice-9 optargs))

(load-extension "libanfo" "init_anfo_guile")

; setting verbosity.  Valid args: #f for silence, int for level (smaller
; values make more noise), one of [debug, info, notice, warning, error,
; critical] for level
(export verbosity)

(define*-public 
  (sort-by-pos #:key (mem 512) (handles 256) genome)
  (prim-sort-by-pos mem handles genome))

(define*-public 
  (sort-by-name #:key (mem 512) (handles 256))
  (prim-sort-by-name mem handles))

(define*-public
  (filter-by-score #:key (slope 7.5) (intercept 20) genome)
  (prim-filter-by-score sclope intercept genome))

(define*-public
  (filter-by-mapq mapq #:key genome)
  (prim-filter-by-mapq mapq genome))

(define*-public
  (ensure-hit #:key genome sequence)
  (prim-ensure-hit genome sequence))
 
(define*-public
  (delete-hit #:key genome sequence)
  (prim-delete-hit genome sequence))
 
(define*-public
  (edit-header #:key editor)
  (prim-edit-header editor))

(define*-public
  (rmdup #:key (slope 7.5) (intercept 20) (maxq 127))
  (prim-rmdup slope intercept maxq))

(export filter-by-length filter-by-qual ensure-multiplicity subsample)

; XXX awkward: region definitions shouldn't need to come from an effing file
(export regions-only not-regions)

; Different ways of merging streams.  All get a single list of streams
; as argument and nothing else.
(export concat-streams merge-streams join-streams mega-merge-streams)

; writes a SAM file to file/stdout, can be restricted to genome
(define*-public
  (write-sam fn #:key genome)
  (prim-write-sam fn genome))

; simple output formats: 3Aln, Glz, Fastq, native text
(export write-3aln write-glz write-text write-fastq)

; Write alignments in Fasta format.  Can be restricted to genome,
; context can be included.
(define*-public
  (write-fasta fn #:key genome context)
  (prim-write-fasta fn genome context))

; Write simple table to file or stdout, optionally restricted to genome.
(define*-public
  (write-table fn #:key genome)
  (prim-write-table fn genome))

; Write statistics to file.
; This is bit unwieldy; realistically we might want to put the
; statistics somewhere in Scheme-land.  But where?
(define*-public
  (write-stats fn #:key genome)
  (prim-write-stats fn genome))

; Creates a pseudo-assembly.  Hits to the given genome are used, else
; everything duct-tapes to the genome it hits best.  Contigs are names
; by the given name followed by a number.
(define*-public 
  (duct-tape #:key genome (name "contig"))
  (prim-duct-tape genome name))

; Reads a stream from a file.  All sorts of streams are supported,
; currently Anfo native files, FastA/FastQ, and the stream is
; transparently decompressed, which currently works for GZip and BZip2.
;
; Hints for decoding some formats are taken as keyword arguments,
; currently the details of the ill-specified FastQ format.  Hints that
; don't apply to the stream actually found are ignored.

(define*-public 
  (read-file fn #:key solexa-scale (origin 33))
  (prim-read-file fn solexa-scale origin))

; Writes an Anfo native file.  Default is best available compression,
; faster or no compression can be requested by setting the compression
; 'level' to something from 0 to 100.  No filename means to write to
; stdout.
(define*-public
  (write-anfo-file #:optional fn #:key (level 100))
  (prim-write-anfo-file fn level))
     
; Copies from one stream to another.  I'm not sure this is even useful;
; composing the two streams and having the result compute its status
; (presumably end-of-stream) should work, too.  Takes two streams as
; arguments.
(export transfer)

; Composes streams into a chain.  If it starts with an input filter, it
; results in an input filter.  If it ends with an output filter, the
; result is an output filter.  A full chain from input to output should
; behave like 'transfer' (XXX check?), only filters results in a filter,
; and everything else will misbehave.
(export chain)
