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

import           Control.Monad.State
import           Data.Attoparsec.Char8
import           Data.Binary.Put
import           Data.Bits ( (.|.), shiftL )
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy.Char8 as B
import           Data.Char ( isSpace, toLower, chr )
import           Data.Foldable hiding ( foldr, concatMap )
import qualified Data.Sequence as Q
import           Data.Sequence ( (|>), (><), ViewR(..) )
import qualified Data.Set as Set
import           Data.Word ( Word8 )
import           System.Environment ( getArgs )
import           System.IO
import           Text.ProtocolBuffers
import           Text.ProtocolBuffers.WireMessage ( getFromBS )

import qualified Config.Contig
import qualified Config.Homolog
import qualified Config.Genome
import qualified Config.Sequence

data Label = Boring | Chimp | Human | Common | Essential deriving (Show, Eq)

data Tree = Empty
          | Node { n_label    :: {-# UNPACK #-} !Label
                 , n_chrom    :: {-# UNPACK #-} !S.ByteString
                 , n_start    :: {-# UNPACK #-} !Int32
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
        = Node lbl (shelve chrom) (fromIntegral $ if str then beg else -end) l r i
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

                       
interesting_trees :: Tree -> [Tree]
interesting_trees t@Node{ n_label = Essential } = [t]
interesting_trees t@Node{ n_label = Common } = interesting_trees (n_left t) ++ interesting_trees (n_right t)
interesting_trees _ = []

data Op a = Delete {-# UNPACK #-} !Word32 | Insert a | Match {-# UNPACK #-} !Word32 | Replace a | NoOp deriving Show

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
                       , out_seqdefs :: [ CodedSeqdef ] }
type Out a = StateT OutputState IO a 

newtype CodedSeqdef = CodedSeqdef B.ByteString

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
    st <- get
    put $! st { out_seqdefs = sd : out_seqdefs st }
    reduce_seqs xs
                
reduce_seq :: Tree -> [B.ByteString] -> Out CodedSeqdef
reduce_seq t bs = do
    lift $ putStrLn $ ".o0O( " ++ S.unpack (n_chrom t) ++ " )"
    p <- gets out_p
    cigars <- r $ consensus t bs
    p' <- gets out_p
    let genome Chimp = "pt2"
        genome Human = "hg18"
        genome _     = error "oops, picked wrong tree"

        mkhom c t' = defaultValue { Config.Homolog.genome = Just . Utf8 . B.pack . genome $ n_label t'
                                  , Config.Homolog.sequence = Just . Utf8 $ B.fromChunks [n_chrom t']
                                  , Config.Homolog.start = n_start t'
                                   , Config.Homolog.cigar = encodeCigar c }

        homologs = Q.fromList . zipWith mkhom cigars $ get_trees Human t ++ get_trees Chimp t

        co = defaultValue { Config.Contig.offset      = p+1
                          , Config.Contig.range_start = 0
                          , Config.Contig.range_end   = p'-p-1
                          , Config.Contig.homolog     = homologs }

        sd = defaultValue { Config.Sequence.contig = Q.singleton co
                          , Config.Sequence.name = Utf8 (B.fromChunks [n_chrom t]) }

        code = runPut $ messagePutM $ sd

    B.length code `seq` return $ CodedSeqdef code

  where
    r [] = return []
    r xs@( (_,ops) : _ ) = r' xs (replicate (length ops) Q.empty) 

    r' [              ] cs = return cs
    r' ( (x,ops) : xs ) cs = do putNybble x ; r' xs $! zipWith' app_op cs ops

        
app_op :: Cigar -> Op Word8 -> Cigar
app_op c o = case ( Q.viewr c, o ) of
    ( c' :> Match n,    Match m   ) -> (|>) c' $! Match $! n+m
    ( c' :> Delete n,   Delete m  ) -> (|>) c' $! Delete $! n+m
    ( c' :> Replace xs, Replace x ) -> (|>) c' $! Replace $! ((|>) xs $! x)
    ( c' :> Insert xs,  Insert x  ) -> (|>) c' $! Insert $! ((|>) xs $! x)
    ( _,                _         ) -> app2 o
  where 
    app2 (Replace x) = (|>) c $! Replace $! (Q.singleton $! x)
    app2 (Insert  x) = (|>) c $! Insert  $! (Q.singleton $! x)
    app2 (Match   n) = (|>) c $! Match n
    app2 (Delete  n) = (|>) c $! Delete n
    app2  NoOp       = c

encodeCigar :: Cigar -> Q.Seq Word32
encodeCigar = foldMap encodeOp
  where
    -- encoding CIGAR as sequence of ints: low two bits are operation (0 == match, 1 == delete, 2 == insert), rest is either a
    -- nucleotide (for insert) or a number (otherwise); replacement becomes deletion + insertion
    encodeOp (Match n)   = Q.singleton $ n `shiftL` 2 .|. 0
    encodeOp (Delete n)  = Q.singleton $ n `shiftL` 2 .|. 1
    encodeOp (Insert xs) = fmap (\x -> fromIntegral x `shiftL` 2 .|. 2) xs
    encodeOp (Replace xs) = encodeOp (Delete (fromIntegral (Q.length xs))) >< encodeOp (Insert xs)
    encodeOp NoOp = Q.empty

{-# INLINE zipWith' #-}
zipWith' :: ( a -> b -> c ) -> [a] -> [b] -> [c]
zipWith' f (x:xs) (y:ys) = ((:) $! f x y) $! zipWith' f xs ys
zipWith' _ _ _ = []

write_dna_file :: FilePath -> [ ( Tree, [B.ByteString] ) ] -> IO ()
write_dna_file path aseqs = 
    withBinaryFile path WriteMode $ \h -> do
        B.hPut h $ B.pack "DNA1\xD\xE\xA\xD\xB\xE\xE\xF"
        putStrLn "done writing header"

        Out { out_seqdefs = seqdefs, out_p = total } <- execStateT (reduce_seqs aseqs) $
                    Out { out_p = 24, out_seqdefs = [], out_byte = 0, out_hdl = h }
        putStrLn "done writing bases"

        index_start <- hTell h
        B.hPut h $ runPut $ messagePutM $ defaultValue {
                Config.Genome.name        = Just $ Utf8 $ B.pack "CA",
                Config.Genome.description = Just $ Utf8 $ B.pack "Human-Chimp common ancestor",
                Config.Genome.sequence    = foldr (\(CodedSeqdef sd) q -> q |> getFromBS messageGetM sd) Q.empty seqdefs,
                Config.Genome.total_size  = total }
        index_end <- hTell h
        putStrLn "done writing metainfo"
        
        hSeek h AbsoluteSeek 4
        B.hPut h $ runPut $ do putWord32le (fromIntegral index_start)
                               putWord32le (fromIntegral $ index_end-index_start)
        putStrLn "done writing pointers"


