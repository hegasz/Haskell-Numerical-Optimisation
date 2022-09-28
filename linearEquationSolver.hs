import Vectors

-- solve Ax = b for x if any exists

type Entry = Double
type Row = [Entry]
type Matrix = [Row]

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
-- @param upperLower a list containing rows of L and U
--                   so upperLower[j] = (L[j],U[j])
--                   where U[j] is row j of U.
-- Performs row reduction on column i of U, yielding
-- pivot at U[i,i] and zeros below
-- Then also fills column i of L with multiplier values
-- that make LU decomposition work.
innerLoop :: Int -> [(Row, Row)] -> [(Row, Row)]
innerLoop i weavedMatrices = foldl cons [head weavedMatrices] (tail weavedMatrices)
    where cons acc (lowerE1, upperE1) = acc ++ [(lowerRow, upperRow)]
            where upperRow = oneRow upperE1 pivotRow multiplier
                  lowerRow = insertToRow lowerE1 i multiplier
                  pivotRow = snd (head acc)
                  multiplier = (upperE1!!i)/(pivotRow!!i)


-- Repeatedly applies innerLoop through all i values
-- as follows:
-- b = innerLoop 0 a
-- (head b) : innerLoop 1 (tail b)
-- etc ...
-- pre: matrix is square
outerLoop :: [(Row,Row)] -> [(Row,Row)]
outerLoop weavedMatrices = foldl cons weavedMatrices [0..length weavedMatrices-1]
    where cons acc i = take i acc ++ innerLoop i (drop i acc)


unWeaver :: [(Row,Row)] -> (Matrix, Matrix)
unWeaver [] = ([],[])
unWeaver ((l,u):xs) = (l:ls, u:us)
    where (ls,us) = unWeaver xs

-- Weaves a input matrix in with an identity matrix
weaver :: Int -> Int -> Matrix -> [(Row,Row)]
weaver i n [] = []
weaver i n (xs:xss) = (idRow, xs) : weaver (i+1) n xss
    where idRow = replicate i 0 ++ [1] ++ replicate (n-i-1) 0

-- @param matrix an input matrix A
-- @returns (L,U) such that L is lower triangular,
--          U is upper triangular, and A = LU.
decompLU :: Matrix -> (Matrix, Matrix)
decompLU matrix = unWeaver . outerLoop . weaver 0 (length matrix) $ matrix

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
-- Performs LU decomposition then solves by forward and
-- backward substitution via the substitution Ux=y into LUx=b
solve :: Matrix -> Vec -> Vec
solve matrix b = backwardSub upper y
    where (lower,upper) = decompLU matrix
          y = forwardSub lower b


a :: Matrix
a = [[2,4,5,6],
     [-1,2,8.5,1],
     [3,8,3,-3],
     [5,2,1.5,6.4]]

-- PARTIAL PIVOTING PLS
b :: Matrix
b = [[0,4,5,6],
     [-1,2,8.5,1],
     [3,8,3,-3],
     [5,2,1.5,6.4]]
