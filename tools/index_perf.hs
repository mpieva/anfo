{-

From calculating match probabilities on idealized sequences, we can
infer that...

- Megablast and Old Anfo are unable to achieve `five nines sensitivity'
  for queries shorter than ~96.

- Matching while allowing two mismatches is always cheaper and/or more
  sensitive than needing perfect matches.

- Below a query length of 23, we cannot realistically achieve `five
  nines'.

- For practical reasons, anything above ~24 match length is approximate
  anyway, so we'll make some educated guesses about the most sensible
  implementation.

- For long seeds, possibly with mismatches allowed, the index lookup
  becomes prohibitively expense, both in terms of time and memory.

- Longer seeds become impractical anyway, because of gaps.

- A rigorous analysis will take considerable effort, but that would be
  wasted, as it becomes increasingly clear that a seedless approach
  (based on a full text index) is superior in every way (except that it
  hasn't been implemented).


The following policy therefore appears as a sensible compromise, but I
won't go the trouble of trying to prove that.  We settle on indexing
words of length twelve every four positions of the genome, which yields
a manageable index and enough flexibility.  Then we require:

Length		Seed Length	    MM per word	    MM total
23		    12		        2		        2
35		    16		        2		        2
43		    20		        2		        2
59		    24		        1		        2
71		    28		        1		        2
87		    32		        1		        2
103		    36		        1		        2
119		    28		        1		        1
139		    32		        1		        1
163		    36		        1		        1

-}


import Data.Array
import Data.List

-- Try and estimate performance of indexing methods.  This code shows clearly what a tall order "five
-- nines" sensitivity is for a seed-based alignment.
--
-- Stolen from Daniel G. Brown: "A Survey of Seeding for Sequence Alignment"
--
--   The probability that an a-letter ungapped alignment, in which positions have probability p of being
--   matched and 1-p of being mismatched, includes k consecutive matching symbols can be computed in
--   O(ak) time.
--

prob_match :: Int -> Double -> [ (Int, Double, Double, Double) ]
prob_match k p = [ (i, memo0 ! (i,0), memo1 ! (i,0), memo2 ! (i,0) ) | i <- [0..50] ]
  where
    -- Let P0(i,j) be the probability that an alignment of length i that is forced to start with j
    -- matching symbols includes a length k with matching symbols.  We seek P0(a,0).
    pee0 i j | i < k     = 0
             | j == k    = 1
             | otherwise = p * memo0 ! (i,j+1) + (1-p) * memo0 ! (i-j-1,0)

    memo0 = listArray ( (0,0), (50,k) )
            [ pee0 i j | (i,j) <- range ( (0,0), (50,k) ) ]


    -- Let P1(i,j) be the probability that an alignment of length i that is forced to start with j
    -- matching symbols includes a length k with all but one matching symbols.  We seek P1(a,0).
    pee1 i j | i < k = 0
             | j == k = 1
             | otherwise =
        -- the symbol at (j+1) is a match with probability p, else we refer to an auxillary function
        p * memo1 ! (i,j+1) + (1-p) * memo'1 ! (i,j,j+1)

    memo1 = listArray ( (0,0), (50,k) )
            [ pee1 i j | (i,j) <- range ( (0,0), (50,k) ) ]
    
    -- Let P'1(i,j1,j2) be the probability that an alignment of length i that is forced to start with j1
    -- matching symbols and has no more than one mismatch in the first j2 symbols includes a length k
    -- with just one mismatching symbol.
    pee'1 i j1 j2 | i < k = 0        -- alignment too short -> 0
                  | j2 <= j1 = error "can't happen"
                  | j2 == k = 1      -- forced match in the beginning -> 1
                  | otherwise =
        -- the symbol at (j1+1) is definitely a mismatch. now the symbol at (j2+1) is a match with
        -- probability p
        p * memo'1 ! (i,j1,j2+1) + (1-p) * memo'1 ! (i-j1-1, j2-j1-1, j2-j1)

    memo'1 = listArray ( (0,0,0), (50,k,k) )
            [ pee'1 i j1 j2 | (i,j1,j2) <- range ( (0,0,0), (50,k,k) ) ]
    


    pee2 i j | i < k     = 0
             | j == k    = 1
             | otherwise = p * memo2 ! (i,j+1) + (1-p) * memo'2 ! (i,j,j+1)

    memo2 = listArray ( (0,0), (50,k) )
            [ pee2 i j | (i,j) <- range ( (0,0), (50,k) ) ]
    
    pee'2 i j1 j2 | i < k     = 0        -- alignment too short -> 0
                  | j2 == k   = 1      -- forced match in the beginning -> 1
                  | otherwise =
        p * memo'2 ! (i,j1,j2+1) + (1-p) * memo''2 ! (i,j1,j2,j2+1)

    memo'2 = listArray ( (0,0,0), (50,k,k) )
            [ pee'2 i j1 j2 | (i,j1,j2) <- range ( (0,0,0), (50,k,k) ) ]
    

    pee''2 i j1 j2 j3 | i < k     = 0
                      | j3 == k   = 1
                      | otherwise =
        p * memo''2 ! (i,j1,j2,j3+1) + (1-p) * memo''2 ! (i-j1-1,j2-j1-1,j3-j1-1,j3-j1)

    memo''2 = listArray ( (0,0,0,0), (50,k,k,k) )
              [ pee''2 i j1 j2 j3 | (i,j1,j2,j3) <- range ( (0,0,0,0), (50,k,k,k) ) ]


costed_match' k' p = map add_cost $ prob_match k' (p^4)
  where
    add_cost (i,p,q,r) = ( i', [(k,0,p,cp), (k,1,q,cq), (k,2,r,cr)] )
      where
        k = 4*k'
        i' = 4*i+3
        cp = fromIntegral ((i'-k+1) `max` 0) / 4^^k
        cq = cp * fromIntegral (1 + k*3)
        cr = cp * fromIntegral (1 + k*3 + k*(k-1)*9`div`2)

        
    
show_solution (i, ss) = shows i . (++) ":  " . foldr (.) id (map show_s ss)
  where
    show_s (k,m,p,c) = shows k . (:) '%' . shows m . (++) " (" . shows (c * 6e9) . (++) "), "

costed_match = map concat' $
    transpose [ costed_match' k 0.98 | k <- [3..12] ]
  where
    concat' xs = let (ks,vs) = unzip xs
                 in (the ks, find_best $ concat vs)
    the (x:xs) | all (x ==) xs = x

find_best = take 5 . sort_by_cost . filter precise_enough

precise_enough (_,_,p,_) = p >= 0.99999
sort_by_cost = sortBy (\x y -> cost x `compare` cost y)
  where cost (_,_,_,c) = c
