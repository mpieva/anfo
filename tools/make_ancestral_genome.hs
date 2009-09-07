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

import Data.Attoparsec.Char8
import qualified Data.ByteString.Lazy.Char8 as B
import Data.Char (isSpace)
import Data.List (transpose)

data Label = Boring | Chimp | Human | Common | Essential deriving Show

data Tree = Empty
          | Node { n_label    :: Label
                 , n_chrom    :: B.ByteString
                 , n_start    :: Integer
                 , n_end      :: Integer
                 , n_strand   :: Bool
                 , n_left     :: Tree
                 , n_right    :: Tree
                 , n_sequence :: String } deriving Show


epo_file :: Parser [Tree]
epo_file = do
    skipSpace
    many comment
    contig `sepBy1` lexeme (string (B.pack "//"))


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
    bases <- many1 (many1 (satisfy (inClass "ACGT-")) <* char '\n')
    phylo (transpose bases)

seqline :: Parser (B.ByteString, B.ByteString, B.ByteString, Integer, Integer, Bool)
seqline = do
    try $ lexeme $ string (B.pack "SEQ")
    family <- takeTill (== '_')
    char '_'
    species <- word
    chrom <- word
    start <- lexeme (integer <?> "start coordinate")
    end <- lexeme (integer <?> "end coordinate")
    strand <- (>=0) <$> lexeme (integer <?> "strand indicator")
    lexeme $ skipMany (notChar '\n')
    return (family, species, chrom, start, end, strand)

treeline :: [(B.ByteString, B.ByteString, B.ByteString, Integer, Integer, Bool)]
         -> Parser ( [String] -> Parser Tree )
treeline seqs = try (lexeme $ string (B.pack "TREE")) *> node <* lexeme (char ';')
  where
    node :: Parser ( [String] -> Parser Tree )
    node =     do char '(' <?> "opening parenthesis"
                  l <- node
                  char ',' <?> "comma"
                  r <- node
                  char ')' <?> "closing parenthesis" 
                  def <- short_seq_def
                  return $ \bases -> do lt <- l bases ; rt <- r bases
                                        find_seq seqs def (make_anc_node lt rt) bases 
           <|> do def <- short_seq_def
                  return $ find_seq seqs def make_leaf_node

short_seq_def :: Parser (B.ByteString, B.ByteString, Integer, Integer, Bool)
short_seq_def = reorder <$> (takeTill (== '_') <* char '_' <?> "binary short name")
                        <*> (rest <* skipWhile (notInClass ",);"))
  where
    rest = try rest1 <|> (combine <$> (takeTill (== '_') <* char '_'<?> "sequence name") <*> rest)
    rest1 = (,,,) [] <$> ((integer <* char '_') <?> "start position")
                     <*> ((integer <* char '[') <?> "end position")
                     <*> ((const True <$> char '+' <|> const False <$> char '-') <?> "strand specifier")
    combine a (b,c,d,e) = (a:b,c,d,e)
    reorder a (bs,c,d,e) = (a, B.intercalate (B.pack "_") bs, c, d, e)
                
make_anc_node :: Tree -> Tree -> B.ByteString -> (Label, Tree, Tree)
make_anc_node l r genome = (new_label (n_label l) (n_label r), l, r)

make_leaf_node :: B.ByteString -> (Label, Tree, Tree)
make_leaf_node genome = (l (B.unpack genome), Empty, Empty)
  where l "Hsap" = Human ; l "Ptro" = Chimp ; l _ = Boring

find_seq :: [(B.ByteString, B.ByteString, B.ByteString, Integer, Integer, Bool)]
         -> (B.ByteString, B.ByteString, Integer, Integer, Bool)
         -> (B.ByteString -> (Label, Tree, Tree)) -> [String] -> Parser Tree
find_seq ((fam,spec,chr1,beg1,end1,str1):ss) d@(genome,chr,beg,end,str) make_children (bs:bases)
    | B.tail genome `B.isPrefixOf` spec && (B.head genome,chr1,beg1,end1,str1) == (B.head fam,chr,beg,end,str)
        = return $ Node lbl chr beg end str l r bs
    | otherwise = find_seq ss d make_children bases 
  where (lbl, l, r) = make_children genome
find_seq [] (genome,chr,_,_,_) _ [] 
    = fail $ "sequence not found: " ++ B.unpack genome ++ "_" ++ B.unpack chr
find_seq _ d _ _ = fail $ "SEQ lines inconsistent with TREE line" ++ show d

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

                       
 
main = do
    (rest, result) <- parse epo_file <$> B.readFile "example.emf"
    print rest
    case result of Right ts -> mapM_ print ts
                   Left err -> putStr err

