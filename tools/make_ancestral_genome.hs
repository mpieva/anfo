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

import Control.Arrow ((&&&))
import Control.Monad.State
import Data.Attoparsec.Char8
import qualified Data.ByteString.Char8 as S
import qualified Data.ByteString.Lazy.Char8 as B
import qualified Data.ByteString.Lazy as BB
import Data.Char ( isSpace, toLower )
import Data.List ( transpose, foldl', groupBy, genericLength )
import qualified Data.Set as Set
import Data.Bits ( (.|.), shiftL )
import Data.Binary.Put
import qualified Data.Sequence as Q
import Data.Sequence ((|>),(<|),(><),ViewR(..),ViewL(..))
import Data.Word ( Word8 )
import System.IO
import Text.ProtocolBuffers

import qualified Config.Genome
import qualified Config.Sequence
import qualified Config.Contig
import qualified Config

import Debug.Trace

data Label = Boring | Chimp | Human | Common | Essential deriving (Show, Eq)

data Tree = Empty
          | Node { n_label    :: {-# UNPACK #-} !Label
                 , n_chrom    :: {-# UNPACK #-} !S.ByteString
                 , n_start    :: {-# UNPACK #-} !Int
                 , n_end      :: {-# UNPACK #-} !Int
                 , n_strand   :: {-# UNPACK #-} !Bool
                 , n_left     :: {-# UNPACK #-} !Tree
                 , n_right    :: {-# UNPACK #-} !Tree
                 , n_sequence :: {-# UNPACK #-} !Int
                 } deriving Show


parse_epo_entry :: B.ByteString -> (Tree, [B.ByteString], B.ByteString)
parse_epo_entry s = case parse epo_meta s of
    (rest, Right phylo) -> case split_sep rest of (bdef, more) -> phylo `seq` (phylo, B.lines bdef, more)
    (rest, Left err)    -> error $ err ++ "near " ++ show (B.take 30 rest)
  where
    epo_meta = do skipSpace ; many comment ; contig

    split_sep s = case B.span (/= '/') s of
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
make_anc_node l r genome = (new_label (n_label l) (n_label r), l, r)

make_leaf_node :: B.ByteString -> (Label, Tree, Tree)
make_leaf_node genome = (l (B.unpack genome), Empty, Empty)
  where l "Hsap" = Human ; l "Ptro" = Chimp ; l _ = Boring

find_seq :: [(B.ByteString, B.ByteString, B.ByteString, Int, Int, Bool)]
         -> (B.ByteString, B.ByteString, Int, Int, Bool)
         -> (B.ByteString -> (Label, Tree, Tree)) -> Int -> Tree
find_seq ((fam,spec,chr1,beg1,end1,str1):ss) d@(genome,chr,beg,end,str) make_children i
    | B.tail genome `B.isPrefixOf` spec && (B.head genome,chr1,beg1,end1,str1) == (B.head fam,chr,beg,end,str)
        = Node lbl (shelve chr) beg end str l r i
    | otherwise = find_seq ss d make_children (i+1)
  where (lbl, l, r) = make_children genome
find_seq [] (genome,chr,_,_,_) _ _
    = error $ "sequence not found: " ++ B.unpack genome ++ "_" ++ B.unpack chr

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

data Op a = Delete !Int | Insert a | Match !Int | Replace a | NoOp deriving Show

-- XXX: drops all meta data
type AncSeq = [( Word8, [Op Word8] )]

-- consensus sequence: common ancestor and all human or chimp seqs are included, structure is taken from ancestor
consensus :: Tree -> [B.ByteString] -> AncSeq
consensus = map . blockcons . getseqs
  where
    getseqs t = fromIntegral (n_sequence t) : get_seqs Human (n_left t) ++ get_seqs Human (n_right t)

    get_seqs _ Empty = []
    get_seqs s t@Node{ n_left = Empty, n_right = Empty, n_label = l } | l == s = [ fromIntegral $ n_sequence t ]
    get_seqs s t = get_seqs s (n_left t) ++ get_seqs s (n_right t)

    blockcons :: [Int64] -> B.ByteString -> ( Word8, [Op Word8] )
    blockcons cs@(c:_) ss = ( cigs `seq` k, cigs )
      where k | B.index ss c == '-' = 0 
              | otherwise           = foldl' (.|.) 0 . map (to_bits' . B.index ss) $ cs
            cigs = map' (mkcigar k . B.index ss) cs
            map' f [] = []
            map' f (x:xs) = ( (:) $! f x ) $! map' f xs

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
main = B.readFile "bar.emf" >>= write_dna_file "foo.dna" . parse_many

parse_many :: B.ByteString -> [ (S.ByteString, AncSeq) ]
parse_many = map mk_cons . drop_dups Set.empty . do_parse
  where 
    do_parse inp | B.all isSpace inp = []
                 | otherwise         = case parse_epo_entry inp of { (tree, bases, rest) ->
                                       [ (t, bases ) | t <- take 1 (interesting_trees tree) ] ++ do_parse rest }

    drop_dups seen [] = []
    drop_dups seen (t:ts) | n `Set.member` seen =     drop_dups               seen  ts
                          | otherwise           = t : drop_dups (Set.insert n seen) ts
        where n = n_chrom (fst t)

    mk_cons (tree, bases) = ( n_chrom tree, tree `consensus` bases )
        
    
{-
put_nybbles :: [Word8] -> Put
put_nybbles (x:y:zs) = putWord8 (x .|. (y `shiftL` 4)) >> put_nybbles zs
put_nybbles [x]      = putWord8 x
put_nybbles []       = return ()
-}

type Cigar = Q.Seq (Op (Q.Seq Word8))

type Out a = StateT OutputState PutM a 

data OutputState = Out { out_p :: !Word32, out_seqdefs :: !(Q.Seq Config.Sequence.Sequence), out_byte :: !Word8 }

putNybble :: Word8 -> Out ()
putNybble w = do
    s@Out { out_p = p, out_byte = v } <- get
    when (p `mod` 2 == 1) (lift $ putWord8 (w .|. (v `shiftL` 4)))
    put $! s { out_byte = w, out_p = succ p }

reduce_seqs :: [ ( S.ByteString, AncSeq ) ] -> Out ()
reduce_seqs [] = do
    putNybble 0
    p <- gets out_p 
    when (p `mod` 2 == 1) (putNybble 0)
    return ()
                         
reduce_seqs ((n,s):ss) = do
    putNybble 0
    _cigars <- reduce_seq s 
    -- something about the seqdef
    reduce_seqs ss

    -- ( 0 : ws, seqdefs {-cs-} ) -- , p'' )
  -- case (reduce_seq $! p+1) (
  -- where -- ( ws,  seqdefs {- , p'' -} ) = {- co `seq` seqdefs' `seq`-} reduce_seqs p {-'-} seqdefs ss

        -- ( ws', _, p' ) = 
        -- seqdefs' = seqdefs -- |> defaultValue { Config.Sequence.contig = Q.singleton co
                         --             , Config.Sequence.name = Utf8 (B.fromChunks [n]) }
        -- co = defaultValue { Config.Contig.offset      = p+1
                          -- , Config.Contig.range_start = 0
                          -- , Config.Contig.range_end   = p'-p-1 }
                
reduce_seq :: [ (Word8, [Op Word8]) ] -> Out [Cigar]
reduce_seq = r (repeat Q.empty) 
  where
    r cs [              ] = return cs
    r cs ( (x,ops) : xs ) = (r $! zipWith' app_op cs ops) xs 

app_op :: Cigar -> Op Word8 -> Cigar
app_op c o = case ( Q.viewr c, o ) of
    ( c' :> Match n,    Match m   ) -> (|>) c' $! Match (n+m)
    ( c' :> Delete n,   Delete m  ) -> (|>) c' $! Delete (n+m)
    ( c' :> Replace xs, Replace x ) -> (|>) c' $! Replace ((|>) xs $! x)
    ( c' :> Replace xs, _         ) -> (c' |> Delete (Q.length xs) |> Insert xs) `app2` o
    ( c' :> Insert xs,  Insert x  ) -> (|>) c' $! Insert ((|>) xs $! x)
    ( _,                _         ) -> app2 c o
  where 
    app2 c (Replace x) = c |> Replace (Q.singleton $! x)
    app2 c (Insert  x) = c |> Insert  (Q.singleton $! x)
    app2 c (Match   n) = c |> Match n
    app2 c (Delete  n) = c |> Delete n
    app2 c NoOp        = c

{-# INLINE zipWith' #-}
zipWith' :: ( a -> b -> c ) -> [a] -> [b] -> [c]
zipWith' f (x:xs) (y:ys) = ((:) $! f x y) (zipWith' f xs ys)
zipWith' f _ _ = []

write_dna_file :: FilePath -> [ ( S.ByteString, AncSeq ) ] -> IO ()
write_dna_file path aseqs = 
    withBinaryFile path WriteMode $ \h -> do
        B.hPut h $ B.pack "DNA1\xD\xE\xA\xD\xB\xE\xE\xF"
        putStrLn "done writing header"

        let ( Out { out_seqdefs = seqdefs, out_p = total }, bytes ) =
                runPutM $ execStateT (reduce_seqs aseqs) (Out { out_p = 24, out_seqdefs = Q.empty, out_byte = 0 })
        B.hPut h bytes
        putStrLn "done writing bases"

        index_start <- hTell h
        B.hPut h $ runPut $ messagePutM $ defaultValue {
                Config.Genome.name        = Just $ Utf8 $ B.pack "ca",
                Config.Genome.description = Just $ Utf8 $ B.pack "Human-Chimp common ancestor",
                Config.Genome.sequence    = seqdefs,
                Config.Genome.total_size  = total }
        index_end <- hTell h
        
        hSeek h AbsoluteSeek 4
        B.hPut h $ runPut $ do putWord32le (fromIntegral index_start)
                               putWord32le (fromIntegral $ index_end-index_start)

