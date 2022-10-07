module LinearEquationSolver where

import VectorsInternal
import Data.List (sortBy)
import Data.Function (on)
import Control.Monad (foldM)


-- Perfoms row operations on one row of a matrix 
oneRow :: Row -> Row -> Entry -> Row
oneRow row topRow multiplier = zipWith (-) row (map (*multiplier) topRow)

insertToRow :: Row -> Int -> Entry -> Row
insertToRow row i x = take i row ++ [x] ++ drop (i+1) row

{-
Together, innerLoop and outerLoop essentially do the following:
for i in range(n)
    for j in range(i+1,n)
        pivot = U[i,i]
        multiplier = U[j,i]/pivot
        L[j,i] = multiplier 
        U[j] = oneRow (upper!!j) (upper!!i) multiplier
-}

-- Simultaneously builds U & L for LU decomposition
-- Performs row reduction on column i of U, yielding
-- pivot at U[i,i] and zeros below
-- Then also fills column i of L with multiplier values
-- that make LU decomposition work.
-- @param weavedMatrices a list containing rows of L and U
--                   so weavedMatrices[j] = (L[j],U[j])
--                   where U[j] is row j of U.
-- @param i the column index to work on
-- pre: pivot is not zero
innerLoop :: Int -> [(Row, Row, Int)] -> [(Row, Row, Int)]
innerLoop i weavedMatrices = foldl cons [head weavedMatrices] (tail weavedMatrices)
    where cons acc (lowerE1, upperE1, trackingI) = acc ++ [(lowerRow, upperRow, trackingI)]
            where upperRow = oneRow upperE1 pivotRow multiplier
                  lowerRow = insertToRow lowerE1 i multiplier
                  pivotRow = (\(_,x,_)->x) (head acc)
                  multiplier = (upperE1!!i)/(pivotRow!!i)





-- Repeatedly applies innerLoop through all i values
-- as follows:
-- b = innerLoop 0 a
-- (head b) : innerLoop 1 (tail b)
-- etc ...
-- Performs partial pivoting for each pivot to ensure
-- non-zero pivots at all times (if possible).
-- pre: matrix is square
outerLoop :: [(Row,Row,Int)] -> Either String [(Row,Row,Int)]
outerLoop weavedMatrices = foldM cons weavedMatrices [0..length weavedMatrices-1]
    where cons acc i = let sanitised = sanitise (drop i acc) i
                       in case sanitised of
                            Nothing -> Left "no non-zero pivots"
                            Just xs -> Right (take i acc ++ innerLoop i xs)




-- Rearranges upper matrix part of weaved matrices to make
-- pivot non-zero, replacing with largest available pivot to
-- decrease round-off errors. Also rearranges tracking indices 
-- to ensure corresponding changes are encoded in permutation
-- matrix later. Returns `Nothing` if no non-zero pivots are
-- found, i.e. matrix is underdetermined or singular.
-- @param i the pivot index within a row
-- @param xs a weaved matrix (upper,lower,trackingIndex)
sanitise :: [(Row,Row,Int)] -> Int -> Maybe [(Row,Row,Int)]
sanitise xs i | abs(((!!i) . (\(_,x,_)->x) . head) xs) > 3e-16 = Just xs
              | abs(((!!i) . fst . head) sorted) <= 3e-16 = Nothing
              | otherwise = Just (zipWith (\l (u,i) -> (l,u,i)) lowerOnly sorted)
    where sorted = sortBy (flip compare `on` (\(xs,_)->xs!!i)) noLower -- decreasing sort
          noLower = map (\(_,upper,index)->(upper,index)) xs
          lowerOnly = map (\(lower,_,_)->lower) xs -- don't rearrange lower matrix 

-- Unweaves lower and upper rows and creates permutation
-- matrix P from trackingIndices such that PA=LU
unWeaver :: [(Row,Row,Int)] -> (Matrix, Matrix, Matrix)
unWeaver [] = ([],[],[])
unWeaver ((l,u,i):xs) = (l:ls, u:us, [if j==i then 1 else 0 | j<-[0..length l-1]]:ps)
    where (ls,us,ps) = unWeaver xs


-- Weaves an input matrix in with an identity matrix and row indices
-- @param n matrix size
weaver :: Int -> Matrix -> [(Row,Row,Int)]
weaver n = weaverInner 0 n
    where weaverInner i n [] = []
          weaverInner i n (xs:xss) = (idRow, xs, i) : weaverInner (i+1) n xss
            where idRow = replicate i 0 ++ [1] ++ replicate (n-i-1) 0


-- LUP decomposition with partial pivoting.
-- @param matrix an input matrix A
-- @returns (L,U,P) such that L is lower triangular,
--          U is upper triangular, P is a permutation
--          matrix and PA = LU.
-- pre: input matrix is square
decompLUP :: Matrix -> Either String (Matrix, Matrix, Matrix)
decompLUP matrix = let result = outerLoop . weaver (length matrix) $ matrix
                  in case result of
                    Left string -> Left string
                    Right weaved -> Right (unWeaver weaved)
    

-- Augments a matrix with a solution and a row index for easier folding
augmenter :: Matrix -> Vec -> [(Row,Entry,Int)]
augmenter matrix b = augmenterInner 0 matrix b
    where augmenterInner _     []       _    = []
          augmenterInner _     _        []   = []
          augmenterInner i (row:rows) (b:bs) = (row,b,i) : augmenterInner (i+1) rows bs

-- Computes undivided result of subtractions in a line of forward subsitution
-- @param line the matrix row we're working with
-- @param val the value in b this line is equal to
-- @param diag the index of the diagonal element in this line
-- @param x vector of solutions so far
solveLine :: Row -> Entry -> Int -> Vec -> Entry
solveLine line val diag x = foldl (-) val (take diag (zipWith (*) x line))

-- equivalent line solver for backward substitution, for upper triangular matrices
solveLineR :: Row -> Entry -> Int -> Vec -> Entry
solveLineR line val diag x = foldr (flip (-)) val (zipWith (*) x (drop (diag+1) line))

-- Solves Lx = b for x
-- pre: L is lower triangular && L,b have correct dimensions
forwardSub :: Matrix -> Vec -> Vec
forwardSub lower b = foldl cons [] (augmenter lower b)
    where cons acc (row,b,i) = acc ++ [unDividedResult/(row!!i)]
            where unDividedResult = solveLine row b i acc

-- Solves Ux = b for x
-- pre: U is upper triangular && U,b have correct dimensions
backwardSub :: Matrix -> Vec -> Vec
backwardSub upper b = foldr cons [] (augmenter upper b)
    where cons (row,b,i) acc = (unDividedResult/(row!!i)) : acc
            where unDividedResult = solveLineR row b i acc


-- Solves Ax = b for x if a solution exists.
-- Performs LU decomposition with partial pivoting 
-- then solves by forward and backward substitution
-- via the substitution Ux=y into LUx=b.
solve :: Matrix -> Vec -> Either String Vec
solve matrix b = let result = decompLUP matrix
                 in case result of 
                    Left _ -> Left "No solution"
                    Right (lower,upper,perm) -> Right (backwardSub upper y)
                            where y = forwardSub lower (multMV perm b)
