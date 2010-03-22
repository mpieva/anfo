{-# LANGUAGE FlexibleInstances, MultiParamTypeClasses, DeriveDataTypeable #-}
{- 
 - Construction of human + common ancestor (of human and chimp) genome.
 - Input is emf files from EPO (compara) output.
 -
 - What do we want:  find all known trees that contain at least human or
 - chimp.  From those, take all the human sequences (none, one or more)
 - and the ancestral sequence and construct a consensus using ambiguity
 - codes.  The state of indels is simply taken from the ancestor.
 - Create CIGAR strings to reconstruct the chimp and human sequences in
 - that alignment, store them along with coordinates.  If contigs are
 - adjacent in human coordinates, join them.  Keep track of covered
 - stretches, check if anything was lost, add it if so.
 -
 - We classify nodes in the phylogenetic tree as follows:  A node is an
 - 'uncommon ancestor' if it has chimp or human descendants, but not
 - both.  It is a 'common ancestor' if it contains both and it is
 - 'boring' if it contains neither.  A 'common ancestor' is essential
 - iff it has at least one 'uncommon' child, it is 'exposed' if it has
 - no ancestor that is itself 'essential'.  Subtrees that are rooted at
 - 'exposed essential common ancestors' go into the CA genome, and so do
 - trees where the root is an uncommon ancestor.
 -
 - For every such tree, the consensus sequence of all chimp and all
 - human sequences is stored along with cigar lines to reconstruct each
 - one of the inputs.  Boring subtrees are pruned.  A tree is 'useful'
 - if one of its children is a human sequence and the other contains
 - chimp sequences only. 
 -
 -}

import           Control.Monad.State hiding ( forM_ )
import           Data.Array.Unboxed
import           Data.Attoparsec.Char8
import           Data.Binary.Put
import           Data.Bits ( (.|.), shiftL )
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy.Char8 as B
import           Data.Char ( isSpace, toLower, chr )
import           Data.Foldable hiding ( foldr, concatMap )
import qualified Data.Map as M
import qualified Data.Sequence as Q
import           Data.Sequence ( (|>), ViewR(..) )
import qualified Data.Set as Set
import           Data.Word ( Word8 )
import           System.Environment ( getArgs )
import           System.IO
import           Text.ProtocolBuffers

import qualified Config.Contig
import qualified Config.Homolog
import qualified Config.Genome
import qualified Config.Sequence

data Label = Boring | Chimp | Human | Common | Essential deriving (Show, Eq)

data Tree = Empty
          | Node { n_label    :: {-# UNPACK #-} !Label
                 , n_chrom    :: {-# UNPACK #-} !S.ByteString
                 , n_start    :: {-# UNPACK #-} !Int32
                 , n_end      :: {-# UNPACK #-} !Int32
                 , n_strand   :: {-# UNPACK #-} !Bool
                 , n_left     :: {-# UNPACK #-} !Tree
                 , n_right    :: {-# UNPACK #-} !Tree
                 , n_sequence :: {-# UNPACK #-} !Int
                 } deriving Show


parse_epo_entry :: B.ByteString -> (Tree, [B.ByteString], B.ByteString)
parse_epo_entry s0 = case parse epo_meta s0 of
    (rest, Right phylo) -> case split_sep rest of (bdef, more) -> phylo `seq` (phylo, B.lines bdef, more)
    (rest, Left err)    -> error $ err ++ "near " ++ show (B.take 30 rest)
  where
    epo_meta = do skipSpace ; many comment ; contig

    split_sep s1 = case B.span (/= '/') s1 of
                    (l,r) | B.null r                       -> (l, r)
                          | B.pack "//\n" `B.isPrefixOf` r -> (l, B.drop 3 r)
                          | otherwise -> error "found '/' but no separator"
                   

lexeme :: Parser a -> Parser a
lexeme p = p <* skipWhile isSpace

word :: Parser B.ByteString
word = lexeme (takeTill isSpace)

comment :: Parser ()
comment = lexeme $ char '#' *> skipWhile (/= '\n')

contig :: Parser Tree
contig = do 
    seqdefs <- many1 (seqline <?> "sequence definition")
    phylo <- treeline seqdefs <?> "tree topology definition"
    lexeme $ try (string (B.pack "DATA") <?> "data keyword")
    return phylo

seqline :: Parser (B.ByteString, B.ByteString, B.ByteString, Int, Int, Bool)
seqline = do
    try $ lexeme $ string (B.pack "SEQ")
    family <- takeTill (== '_')
    char '_'
    species <- word
    chrom <- word
    start <- fromIntegral <$> lexeme (integer <?> "start coordinate")
    end <- fromIntegral <$> lexeme (integer <?> "end coordinate")
    strand <- (>=0) <$> lexeme (integer <?> "strand indicator")
    lexeme $ skipMany (notChar '\n')
    return (family, species, chrom, start, end, strand)

treeline :: [(B.ByteString, B.ByteString, B.ByteString, Int, Int, Bool)] -> Parser Tree
treeline seqs = try (lexeme $ string (B.pack "TREE")) *> node <* lexeme (char ';')
  where
    node :: Parser Tree
    node =     do char '(' <?> "opening parenthesis"
                  l <- node
                  char ',' <?> "comma"
                  r <- node
                  char ')' <?> "closing parenthesis" 
                  def <- short_seq_def
                  return $ find_seq seqs def (make_anc_node l r) 0
           <|> do def <- short_seq_def
                  return $ find_seq seqs def make_leaf_node 0

short_seq_def :: Parser (B.ByteString, B.ByteString, Int, Int, Bool)
short_seq_def = reorder <$> (takeTill (== '_') <* char '_' <?> "binary short name")
                        <*> (rest <* skipWhile (notInClass ",);"))
  where
    rest = try rest1 <|> (combine <$> (takeTill (== '_') <* char '_'<?> "sequence name") <*> rest)
    rest1 = (,,,) [] <$> (fromIntegral <$> (integer <* char '_') <?> "start position")
                     <*> (fromIntegral <$> (integer <* char '[') <?> "end position")
                     <*> ((const True <$> char '+' <|> const False <$> char '-') <?> "strand specifier")
    combine a (b,c,d,e) = (a:b,c,d,e)
    reorder a (bs,c,d,e) = (a, B.intercalate (B.pack "_") bs, c, d, e)
                
make_anc_node :: Tree -> Tree -> B.ByteString -> (Label, Tree, Tree)
make_anc_node l r _genome = (new_label (n_label l) (n_label r), l, r)

make_leaf_node :: B.ByteString -> (Label, Tree, Tree)
make_leaf_node genome = (l (B.unpack genome), Empty, Empty)
  where l "Hsap" = Human ; l "Ptro" = Chimp ; l _ = Boring

find_seq :: [(B.ByteString, B.ByteString, B.ByteString, Int, Int, Bool)]
         -> (B.ByteString, B.ByteString, Int, Int, Bool)
         -> (B.ByteString -> (Label, Tree, Tree)) -> Int -> Tree
find_seq ((fam,spec,chr1,beg1,end1,str1):ss) d@(genome,chrom,beg,end,str) make_children i
    | B.tail genome `B.isPrefixOf` spec && (B.head genome,chr1,beg1,end1,str1) == (B.head fam,chrom,beg,end,str)
        = Node lbl (shelve chrom) (fromIntegral beg) (fromIntegral end) str l r i
    | otherwise = find_seq ss d make_children (i+1)
  where (lbl, l, r) = make_children genome
find_seq [] (genome,chrom,_,_,_) _ _
    = error $ "sequence not found: " ++ B.unpack genome ++ "_" ++ B.unpack chrom

shelve :: B.ByteString -> S.ByteString
shelve s = case B.toChunks s of [] -> S.empty ; [x] -> S.copy x ; xs -> S.concat xs

new_label :: Label -> Label -> Label 
new_label Boring Essential = Common
new_label Boring x = x
new_label Essential Boring = Common
new_label x Boring = x
new_label Chimp Chimp = Chimp
new_label Chimp _ = Essential
new_label _ Chimp = Essential
new_label Human Human = Human
new_label Human _ = Essential
new_label _ Human = Essential
new_label _ _ = Common

                       
-- get either the essential common ancestors from a tree where the root
-- is a common ancestor, or the single non-boring leaf or the only
-- ancestor that has two non-boring children
interesting_trees :: Tree -> [Tree]
interesting_trees t@Node{ n_label = Human } = interesting_trees_lbl Human t
interesting_trees t@Node{ n_label = Chimp } = interesting_trees_lbl Chimp t
interesting_trees t = interesting_trees_common t

interesting_trees_lbl :: Label -> Tree -> [Tree]
interesting_trees_lbl l   Node{ n_label = l' } | l /= l' = []
interesting_trees_lbl _ t@Node{ n_left = Empty, n_right = Empty } = [t]
interesting_trees_lbl l t@Node{ n_left = tl, n_right = tr }
    | n_label tl == l && n_label tr == l = [t]
    | otherwise = interesting_trees_lbl l tl ++ interesting_trees_lbl l tr
interesting_trees_lbl _ _ = [] 

interesting_trees_common :: Tree -> [Tree]
interesting_trees_common t@Node{ n_label = Essential } = [t]
interesting_trees_common t@Node{ n_label = Common } = interesting_trees_common (n_left t) ++ interesting_trees_common (n_right t)
interesting_trees_common _ = []

data Op a = Delete {-# UNPACK #-} !Word32 | Insert a | Match {-# UNPACK #-} !Word32 | Replace a | InsN {-# UNPACK #-} !Word32 | ReplN {-# UNPACK #-} !Word32 | NoOp deriving Show

get_trees :: Label -> Tree -> [Tree]
get_trees _ Empty = []
get_trees s t@Node{ n_left = Empty, n_right = Empty, n_label = l } | l == s = [t]
get_trees s t = get_trees s (n_left t) ++ get_trees s (n_right t)


-- consensus sequence: common ancestor and all human or chimp seqs are included, structure is taken from ancestor, cigar lines are
-- created for all C and H sequences
consensus :: Tree -> [B.ByteString] -> [( Word8, [Op Word8] )]
consensus t = map $ blockcons (map fromIntegral . map n_sequence $ t : get_trees Human t)
                              (map fromIntegral . map n_sequence $ get_trees Human t ++ get_trees Chimp t)
  where
    blockcons :: [Int64] -> [Int64] -> B.ByteString -> ( Word8, [Op Word8] )
    blockcons cs os ss = ( k, cigs )
      where k | null cs                     = 15        -- N because nothing goes in at all
              | B.index ss (head cs) == '-' = 0         -- gap
              | otherwise                   = foldl' (.|.) 0 . map (to_bits' . B.index ss) $ cs
            cigs = map (mkcigar k . B.index ss) os

    to_bits x = case toLower x of 
                    '-' -> 0
                    'a' -> 1
                    'c' -> 2
                    'm' -> 3
                    't' -> 4
                    'w' -> 5
                    'y' -> 6
                    'h' -> 7
                    'g' -> 8
                    'r' -> 9
                    's' -> 10
                    'v' -> 11
                    'k' -> 12
                    'd' -> 13
                    'b' -> 14
                    'n' -> 15
                    _   -> error $ "does not code for nucleotide: " ++ show x
                    
    to_bits' 'N' = 0
    to_bits' 'n' = 0
    to_bits' x = to_bits x

    mkcigar 0 '-'                 = NoOp
    mkcigar 0  b                  = Insert $! to_bits b
    mkcigar _ '-'                 = Delete 1
    mkcigar a  b | a == to_bits b = Match 1
                 | True           = Replace $! to_bits b

main :: IO ()
main = getArgs              >>= \ (outf : infs) ->
       mapM B.readFile infs >>= 
       write_dna_file outf . concatMap parse_many

parse_many :: B.ByteString -> [ ( Tree, [B.ByteString] ) ]
parse_many = drop_dups Set.empty . do_parse
  where 
    do_parse inp | B.all isSpace inp = []
                 | otherwise         = case parse_epo_entry inp of { (tree, bases, rest) ->
                                       [ (t, bases ) | t <- interesting_trees tree ] ++ do_parse rest }

    drop_dups _seen [] = []
    drop_dups  seen (t:ts) | n `Set.member` seen =     drop_dups               seen  ts
                           | otherwise           = t : drop_dups (Set.insert n seen) ts
        where n = n_chrom (fst t)

    
type Cigar = Q.Seq (Op (Q.Seq Word8))

data OutputState = Out { out_p :: {-# UNPACK #-} !Word32
                       , out_byte :: {-# UNPACK #-} !Word8
                       , out_hdl :: {-# UNPACK #-} !Handle
                       , out_seqdefs :: ![Contig]
                       , out_missing_human :: !Missing
                       , out_missing_chimp :: !Missing }
type Out a = StateT OutputState IO a 

-- newtype CodedSeqdef = CodedSeqdef { unCodedSeqdef :: S.ByteString }

putNybble :: Word8 -> Out ()
putNybble w = do
    s@Out { out_p = p, out_byte = v, out_hdl = h } <- get
    when (p `mod` 2 == 1) $
        lift $ hPutChar h $ chr $ fromIntegral $ w .|. (v `shiftL` 4)
    put $! s { out_byte = w, out_p = succ p }

reduce_seqs :: [ ( Tree, [B.ByteString] ) ] -> Out ()
reduce_seqs [] = do
    putNybble 0
    p <- gets out_p 
    when (p `mod` 2 == 1) (putNybble 0)
    return ()
                         
reduce_seqs (x:xs) = do
    putNybble 0
    sd <- uncurry reduce_seq x
    let humans = get_trees Human (fst x)
        chimps = get_trees Chimp (fst x)
        mark'em = foldl' (\m t -> mark (min (n_start t) (n_end t) - 1) (max (n_start t) (n_end t)) (n_chrom t) m)

    st <- get
    put $! st { out_seqdefs = sd : out_seqdefs st
              , out_missing_human = mark'em (out_missing_human st) humans
              , out_missing_chimp = mark'em (out_missing_chimp st) chimps }
    reduce_seqs xs
                
data Homolog = Hom { hom_label :: {-# UNPACK #-} !Label
                   , hom_chrom :: {-# UNPACK #-} !S.ByteString
                   , hom_start :: {-# UNPACK #-} !Int32
                   , hom_cigar :: {-# UNPACK #-} !(UArray Int Word32) }

data Contig = Ctg { ctg_offset :: {-# UNPACK #-} !Word32
                  , ctg_length :: {-# UNPACK #-} !Word32
                  , ctg_name   :: {-# UNPACK #-} !S.ByteString
                  , ctg_homologs :: ![ Homolog ] }

reduce_seq :: Tree -> [B.ByteString] -> Out Contig
reduce_seq t bs = do
    lift $ hPutStrLn stderr $ ".o0O( " ++ S.unpack (n_chrom t) ++ " )"
    p <- gets out_p
    cigars <- r $ consensus t bs
    p' <- gets out_p

    let mkhom c t' = Hom (n_label t') (n_chrom t') start (listArray' (encodeCigar c))
            where start = if n_strand t' then n_start t' else - n_end t'
                  listArray' xs = listArray (1,length xs) xs

    return $! Ctg (p+1) (p'-p) (n_chrom t) $!
        zipWith' mkhom cigars $ get_trees Human t ++ get_trees Chimp t

  where
    r [] = return []
    r xs@( (_,ops) : _ ) = r' xs (replicate (length ops) Q.empty) 

    r' [              ] cs = return cs
    r' ( (x,ops) : xs ) cs = do putNybble x ; r' xs $! zipWith' app_op cs ops

        
app_op :: Cigar -> Op Word8 -> Cigar
app_op c o = case ( Q.viewr c, o ) of
    ( c' :> Match n,    Match m   ) -> (|>) c' $! Match $! n+m
    ( c' :> InsN n,     InsN m    ) -> (|>) c' $! InsN $! n+m
    ( c' :> InsN n,     Insert 15 ) -> (|>) c' $! InsN $! n+1
    ( c' :> Delete n,   Delete m  ) -> (|>) c' $! Delete $! n+m
    ( c' :> ReplN n,    ReplN m   ) -> (|>) c' $! ReplN $! n+m
    ( c' :> ReplN n,    Replace 15) -> (|>) c' $! ReplN $! n+1
    ( c' :> Replace xs, Replace x ) | x /= 15 -> (|>) c' $! Replace $! ((|>) xs $! x)
    ( c' :> Insert xs,  Insert x  ) | x /= 15 -> (|>) c' $! Insert $! ((|>) xs $! x)
    ( _,                _         ) -> app2 o
  where 
    app2 (Replace 15)= (|>) c $! ReplN 1
    app2 (Replace x) = (|>) c $! Replace $! (Q.singleton $! x)
    app2 (Insert 15) = (|>) c $! InsN 1
    app2 (Insert  x) = (|>) c $! Insert  $! (Q.singleton $! x)
    app2 (Match   n) = (|>) c $! Match n
    app2 (Delete  n) = (|>) c $! Delete n
    app2 (InsN    n) = (|>) c $! InsN n
    app2 (ReplN   n) = (|>) c $! ReplN n
    app2  NoOp       = c

encodeCigar :: Cigar -> [Word32]
encodeCigar c = foldMap encodeOp c []
  where
    -- encoding CIGAR as sequence of ints: low two bits are operation (0 == match, 1 == delete, 2 == insert), rest is either a
    -- nucleotide (for insert) or a number (otherwise); replacement becomes deletion + insertion
    encodeOp :: Op (Q.Seq Word8) -> [Word32] -> [Word32]
    encodeOp (Match n)   = (:) (n `shiftL` 2 .|. 0)
    encodeOp (Delete n)  = (:) (n `shiftL` 2 .|. 1)
    encodeOp (Insert xs) = foldMap (\x -> (:) (fromIntegral x `shiftL` 2 .|. 2)) xs
    encodeOp (InsN n)    = (:) (n `shiftL` 2 .|. 3)
    encodeOp (Replace xs) = encodeOp (Delete (fromIntegral (Q.length xs))) . encodeOp (Insert xs)
    encodeOp (ReplN n) = encodeOp (Delete n) . encodeOp (InsN n)
    encodeOp NoOp = id

{-# INLINE zipWith' #-}
zipWith' :: ( a -> b -> c ) -> [a] -> [b] -> [c]
zipWith' f (x:xs) (y:ys) = ((:) $! f x y) $! zipWith' f xs ys
zipWith' _ _ _ = []

write_dna_file :: FilePath -> [ ( Tree, [B.ByteString] ) ] -> IO ()
write_dna_file path aseqs = 
    withBinaryFile path WriteMode $ \h -> do
        B.hPut h $ B.pack "DNA1\xD\xE\xA\xD\xB\xE\xE\xF"
        hPutStrLn stderr "done writing header"

        Out { out_seqdefs = seqdefs, out_p = total, out_missing_human = missing_human, out_missing_chimp = missing_chimp }
            <- execStateT (reduce_seqs aseqs) $
                    Out { out_p = 24, out_seqdefs = [], out_byte = 0, out_hdl = h
                        , out_missing_human = human_tab, out_missing_chimp = chimp_tab }
        hPutStrLn stderr "done writing bases"

        {- 
        forM_ seqdefs $ \c -> do
            forM_ (ctg_homologs c) $ \hom -> do
                hPrint stderr . hom_chrom $ h
                hPrint stderr . uncurry (-) . bounds . hom_cigar $ h
            hPrint stderr $ ctg_name c
               
        hPutStrLn stderr "metainfo looks fine"
        -}

        let genome Chimp = "pt2"
            genome Human = "hg18"
            genome _     = error "oops, picked wrong tree"

            repack_homolog hom = defaultValue { 
                    Config.Homolog.genome = Just . Utf8 . B.pack . genome $ hom_label hom,
                    Config.Homolog.sequence = Just . Utf8 $ B.fromChunks [hom_chrom hom],
                    Config.Homolog.start = hom_start hom,
                    Config.Homolog.cigar = Q.fromList $ elems $ hom_cigar hom }

            repack_contig c = defaultValue {
                                Config.Sequence.contig = Q.singleton $ defaultValue {
                                  Config.Contig.offset      = ctg_offset c,
                                  Config.Contig.range_start = 0,
                                  Config.Contig.range_end   = ctg_length c - 1,
                                  Config.Contig.homolog     = Q.fromList $ map repack_homolog $ ctg_homologs c },
                                Config.Sequence.name = Utf8 (B.fromChunks [ctg_name c]) }

        index_start <- hTell h
        B.hPut h $ runPut $ messagePutM $ defaultValue {
                Config.Genome.name        = Just $ Utf8 $ B.pack "CA",
                Config.Genome.description = Just $ Utf8 $ B.pack "Human-Chimp common ancestor",
                Config.Genome.sequence    = Q.fromList $ map repack_contig $ reverse seqdefs,
                Config.Genome.total_size  = total }

        index_end <- hTell h
        hPutStrLn stderr "done writing metainfo"
        
        hSeek h AbsoluteSeek 4
        B.hPut h $ runPut $ do putWord32le (fromIntegral index_start)
                               putWord32le (fromIntegral $ index_end-index_start)
        hPutStrLn stderr "done writing pointers\n"
    
        let show_inner_map = M.foldWithKey (\a b k -> shows a . ('-':) . shows b . (", "++) . k) id
            show_outer_map = M.foldWithKey (\c m k -> ("  " ++) . (S.unpack c ++) . (": " ++) . show_inner_map m . ('\n':) . k) id

        putStrLn "Missing from human: "
        putStr $ show_outer_map missing_human "\n"

        putStrLn "Missing from chimp: "
        putStr $ show_outer_map missing_chimp "\n"


type Missing = M.Map S.ByteString (M.Map Int32 Int32)

human_tab :: Missing
human_tab = M.fromList [ "1" ==> 247249719, "1_random" ==> 1663265, "10" ==> 135374737, "10_random" ==> 113275, "11" ==> 134452384,
    "11_random" ==> 215294, "12" ==> 132349534, "13" ==> 114142980, "13_random" ==> 186858, "14" ==> 106368585, "15" ==> 100338915,
    "15_random" ==> 784346, "16" ==> 88827254, "16_random" ==> 105485, "17" ==> 78774742, "17_random" ==> 2617613, "18" ==>
    76117153, "18_random" ==> 4262, "19" ==> 63811651, "19_random" ==> 301858, "2" ==> 242951149, "2_random " ==> 185571, "20" ==>
    62435964, "21" ==> 46944323, "21_random" ==> 1679693, "22" ==> 49691432, "22_random" ==> 257318, "3" ==> 199501827, "3_random"
    ==> 749256, "4" ==> 191273063, "4_random" ==> 842648, "5" ==> 180857866, "5_random" ==> 143687, "6" ==> 170899992, "6_random"
    ==> 1875562, "7" ==> 158821424, "7_random" ==> 549659, "8" ==> 146274826, "8_random" ==> 943810, "9" ==> 140273252, "9_random"
    ==> 1146434, "M" ==> 16571, "X" ==> 154913754, "X_random" ==> 1719168, "Y" ==> 57772954 ]
  where c ==> l = (S.pack c, M.singleton 0 l)

chimp_tab :: Missing
chimp_tab = M.fromList [ "1" ==> 229974691, "1_random" ==> 9420409, "10" ==> 135001995, "10_random" ==> 8402541, "11" ==> 134204764,
    "11_random" ==> 8412303, "12" ==> 135371336, "12_random" ==> 2259969, "13" ==> 115868456, "13_random" ==> 9231045, "14" ==>
    107349158, "14_random" ==> 2108736, "15" ==> 100063422, "15_random" ==> 3087076, "16" ==> 90682376, "16_random" ==> 6370548,
    "17" ==> 83384210, "17_random" ==> 5078517, "18" ==> 77261746, "18_random" ==> 1841920, "19" ==> 64473437, "19_random" ==>
    2407237, "20" ==> 62293572, "20_random" ==> 1792361, "21" ==> 46489110, "22" ==> 50165558, "22_random" ==> 1182457, "2a" ==>
    114460064, "2a_random" ==> 3052259, "2b" ==> 248603653, "2b_random" ==> 2186977, "3" ==> 203962478, "3_random" ==> 3517036, "4"
    ==> 194897272, "4_random" ==> 6711082, "5" ==> 183994906, "5_random" ==> 3159943, "6" ==> 173908612, "6_random" ==> 9388360, "7"
    ==> 160261443, "7_random" ==> 7240870, "8" ==> 145085868, "8_random" ==> 7291619, "9" ==> 138509991, "9_random" ==> 7733331, "M"
    ==> 16554, "Un" ==> 58616431, "X" ==> 155361357, "X_random" ==> 3548706, "Y" ==> 23952694, "Y_random" ==> 772887 ] 
  where c ==> l = (S.pack c, M.singleton 0 l)

-- mark an interval: find all overlapping intervals, cut them back,
-- reassemble a finite map
-- Intervals are sorted by start coordinate.  Splitlookup is a good
-- start: candidates are the last interval in the left part, the direct
-- hit, and the minima of the right part.

mark :: Int32 -> Int32 -> S.ByteString -> Missing -> Missing
mark l r c  _ | l `seq` r `seq` c `seq` False = undefined
mark l r c m0 = case M.lookup c m0 of 
                    Just m -> (M.insert c $! mark1 m) m0 
                    Nothing -> m0
  where 
    mark1 m = mleft' `M.union` hit' `M.union` mright' mright
      where
        (mleft, mhit, mright) = M.splitLookup l m
        hit' = case mhit of
                    Nothing -> M.empty 
                    Just i -> markI (l,i)

        mleft' = case M.maxViewWithKey mleft of
                    Nothing -> M.empty
                    Just (i,mleft1) -> mleft1 `M.union` markI i

        mright' mr = case M.maxViewWithKey mr of
                        Nothing -> M.empty
                        Just ((u,_),_) | u >= r -> mr
                        Just (i,mr') -> markI i `M.union` mright' mr'

    markI (u,v) | l >= v || r <= u = M.singleton u v
                | l >  u && r <  v = M.singleton u l `M.union` M.singleton r v
                | l >  u           = M.singleton u l
                |           r <  v = M.singleton r v
                | otherwise        = M.empty


