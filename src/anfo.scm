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
; (use modules (anfo))

(define-module (anfo))
(use-modules (ice-9 optargs))

(load-extension "libanfo" "init_anfo_guile")

; Stream* mk_sort_by_pos( const ParamBlock& p )
; { return new SortingStream<by_genome_coordinate>( parse_int( p.arg, 1024 ) * 1024 * 1024, 256, by_genome_coordinate(p.genome) ) ; }
; 
; Stream* mk_sort_by_name( const ParamBlock& p )
; { return new SortingStream<by_seqid>( parse_int( p.arg, 1024 ) * 1024 * 1024 ) ; }
; 
; Stream* mk_filter_by_length( const ParamBlock& p )
; { return new LengthFilter( parse_int( p.arg ) ) ; }
; 
; Stream* mk_filter_by_score( const ParamBlock& p )
; { return new ScoreFilter( p.slope, p.intercept, p.genome ) ; }
; 
; Stream* mk_filter_by_mapq( const ParamBlock& p )
; { return new MapqFilter( p.genome, parse_int( p.arg ) ) ; }
; 
; Stream* mk_filter_by_hit( const ParamBlock& p )
; { return new HitFilter( p.genome, p.arg ) ; }
; 
; Stream* mk_delete_hit( const ParamBlock& p )
; { return new IgnoreHit( p.genome, p.arg ) ; }
; 
; Stream* mk_filter_qual( const ParamBlock& p )
; { return new QualFilter( parse_int( p.arg ) ) ; }
; 
; Stream* mk_filter_multi( const ParamBlock& p )
; { return new MultiFilter( parse_int( p.arg, 2 ) ) ; }
; 
; Stream* mk_subsample( const ParamBlock& p )
; { return new Subsample( parse_float( p.arg ) ) ; }
; 
; Stream* mk_edit_header( const ParamBlock& p )
; { return new RepairHeaderStream( p.arg ) ; }
; 
; Stream* mk_rmdup( const ParamBlock& p )
; { return new RmdupStream( p.slope, p.intercept, parse_int( p.arg, 127 ) ) ; }
; 
; Stream* mk_regions_only( const ParamBlock& p )
; { return new InsideRegion( p.arg ) ; }
; 
; Stream* mk_not_regions( const ParamBlock& p )
; { return new OutsideRegion( p.arg ) ; }
; 
; StreamBundle* mk_merge( const ParamBlock& )
; { return new MergeStream() ; }
; 
; StreamBundle* mk_join( const ParamBlock& )
; { return new BestHitStream() ; }
; 
; StreamBundle* mk_mega_merge( const ParamBlock& )
; { return new MegaMergeStream() ; }
; 
; StreamBundle* mk_concat( const ParamBlock& )
; { return new ConcatStream() ; }
; 
; Stream* mk_output_text( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new TextWriter( 1 ) : new TextWriter( p.arg ) ; } 
; 
; Stream* mk_output_sam( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new SamWriter( cout.rdbuf(), p.genome ) : new SamWriter( p.arg, p.genome ) ; } 
; 
; Stream* mk_output_glz( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new GlzWriter( 1 ) : new GlzWriter( p.arg ) ; } 
; 
; Stream* mk_output_3aln( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new ThreeAlnWriter( std::cout.rdbuf() ) : new ThreeAlnWriter( p.arg ) ; } 
; 
; Stream* mk_output_fasta( const ParamBlock& p )
; {
; 	return is_stdout( p.arg ) ? new FastaAlnWriter( cout.rdbuf(), p.genome, p.context )
; 	                          : new FastaAlnWriter( p.arg, p.genome, p.context ) ; 
; } 
; 
; Stream* mk_output_fastq( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new FastqWriter( cout.rdbuf() ) : new FastqWriter( p.arg ) ; } 
; 
; Stream* mk_output_table( const ParamBlock& p )
; { return is_stdout( p.arg ) ? new TableWriter( cout.rdbuf(), p.genome ) : new TableWriter( p.arg, p.genome ) ; }
; 
; Stream* mk_duct_tape( const ParamBlock& p )
; { return new DuctTaper( p.genome, p.arg ? p.arg : "contig" ) ; }
; 
; Stream* mk_stats( const ParamBlock& p )
; { return new StatStream( p.arg, p.genome ) ; }



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

