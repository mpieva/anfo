// read ANFO files, sort and merge them
// Outline: keep a queue of opened files and a deque of results.  A new
// file (e.g. from command line) is added to the queue if it is sorted,
// else it is read and added to the deque.  If the deque becomes large,
// it is sorted and written to tempory storage, which is then added to
// the queue.  If the queue becomes long, it is merged into temporary
// storage, the new file is added to the queue.  Finally, all open files
// and the contents of the deque are merged and written out.

struct by_genome_coordinate {
    bool operator() ( const Result *a, const Result* b ) {
	if( a->has_best_to_genome() && !b->has_best_to_genome() ) return true ;
	if( !a->has_best_to_genome() ) return false ;
	const Hit& u = a->best_to_genome(), v = b->best_to_genome() ;
	if( u.sequence() < v.sequence() ) return true ;
	if( v.sequence() < u.sequence() ) return false ;
	if( u.start_pos() + max(u.aln_length(),0) < v.start_pos() + max(v.aln_length(),0) ) return true ;
	if( v.start_pos() + max(v.aln_length(),0) < u.start_pos() + max(u.aln_length(),0) ) return false ;
	return false ;
    }
} ;


    deque< const Result* > arr ;
    for( map<string,Result>::const_iterator b = buffer.begin(), e = buffer.end() ; b!=e ; ++b )
	arr.push_back( &b->second ) ;
    sort( arr.begin(), arr.end(), by_genome_coordinate() ) ;
    for( deque< const Result* >::const_iterator b = arr.begin(), e = arr.end() ; b!=e ; ++b )
	write_delimited_message( cos, 2, **b ) ;

