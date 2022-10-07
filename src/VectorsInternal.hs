module VectorsInternal
where

type Entry = Double
type Vec = [Entry]
type Row = [Entry]
type Matrix = [Row]

equalsMM :: Matrix -> Matrix -> Bool
equalsMM a b = all (== True) (zipWith (\x y -> if abs(x-y)<1e-15 then True else False) (concat a) (concat b))

equalsRR :: Row -> Row -> Bool
equalsRR r q = all (== True) (zipWith (\x y -> if abs(x-y)<1e-15 then True else False) r q)

equalsVV :: Vec -> Vec -> Bool
equalsVV x y = all (== True) (zipWith (\x y -> if abs(x-y)<1e-15 then True else False) x y)

vecToMatrix :: Vec -> Matrix
vecToMatrix = map (\x->[x])

matrixToVec :: Matrix -> Vec
matrixToVec = map head


-- Matrix transpose
transpose :: Matrix -> Matrix
transpose [] = []
transpose ([]:_) = []
transpose xss = heads : transpose tails
    where (heads, tails) = unzip [(x,xs) | x:xs <- xss]


-- Matrix-Matrix multiplication
-- pre: correct dimensions
multMM :: Matrix -> Matrix -> Matrix
multMM ass bss = map (\as -> map (dot as) tss) ass
    where tss = transpose bss


-- Matrix-Vector multiplication
-- pre: correct dimensions
multMV :: Matrix -> Vec -> Vec
multMV a x = matrixToVec $ multMM a (vecToMatrix x)


-- dot product
-- pre: correct dimensions
dot :: Vec -> Vec -> Double
dot xs ys | length xs == length ys = sum (zipWith (*) xs ys)
          | otherwise = error "Incompatible dimensions" 


-- subtraction for vectors
(-^) :: Vec -> Vec -> Vec
(-^) xs ys | length xs == length ys = zipWith (-) xs ys
           | otherwise = error "Incompatible dimensions" 

-- addition for vectors
(+^) :: Vec -> Vec -> Vec
(+^) xs ys | length xs == length ys = zipWith (+) xs ys
           | otherwise = error "Incompatible dimensions" 

-- scalar multiplication of vector
(*^) :: Double -> Vec -> Vec
(*^) a = map (*a)
