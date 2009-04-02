{- Partitions a list of chromosomes with lengths into bins of
 - approximately equal size.  Suggested use:
 -
 - twoBitInfo foo.2bit stdout | runhaskell bins.hs N        -}

import Data.List (sortBy)
import System.Environment (getArgs)

greedy_bin :: Integer -> [(a,Integer)] -> [([a],Integer)]
greedy_bin 1 xs = [(map fst xs, sum $ map snd xs)]
greedy_bin n xs = fill_one [] 0 xs []
  where
    total = sum (map snd xs)
    fill_one as cur ((b,s):bs) cs
        | cur + s < (total - cur) `div` (n-1) = fill_one (b:as) (cur+s) bs cs
        | otherwise                           = fill_one as cur bs ((b,s):cs)
    fill_one as cur [] cs = (as,cur) : greedy_bin (n-1) (reverse cs)

main = getArgs >>= \ [nbins] -> getContents >>=
       mapM_ (putStrLn . unwords . fst) . greedy_bin (read nbins) . 
       sortBy (\a b -> snd b `compare` snd a) .
       map (\[n,l] -> (n,read l)) . map words . lines
