module Vectors(Vec, dot, (-^), (*^))
where

type Vec = [Double]

-- dot product
dot :: Vec -> Vec -> Maybe Double
dot xs ys | length xs == length ys = Just (sum (zipWith (*) xs ys))
          | otherwise = Nothing 


-- subtraction for vectors
(-^) :: Vec -> Vec -> Maybe Vec
(-^) xs ys | length xs == length ys = Just (zipWith (-) xs ys)
           | otherwise = Nothing 

-- scalar multiplication of vector
(*^) :: Double -> Vec -> Vec
(*^) a = map (*a)