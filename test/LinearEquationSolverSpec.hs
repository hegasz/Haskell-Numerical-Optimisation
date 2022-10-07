module LinearEquationSolverSpec where
import SpecHelper

import LinearEquationSolver
import VectorsInternal
import Data.Either
import Data.Maybe
import Test.QuickCheck


a :: [(Row, Row, Int)]
a = [([1,0,0,0],[2,4,5,6],0),
           ([0,1,0,0],[-1,2,8.5,1],1),
           ([0,0,1,0],[3,8,3,-3],2),
           ([0,0,0,1],[5,2,1.5,6.4],3)]

    
aZeroResult :: [(Row, Row, Int)]
aZeroResult = [([1,0,0,0],[2,4,5,6],0),
               ([-0.5,1,0,0],[0,4,11,4],1),
               ([1.5,0,1,0],[0,2,-4.5,-12],2),
               ([2.5,0,0,1],[0,-8,-11,-8.6],3)]


aOneResult :: [(Row, Row, Int)]
aOneResult = [([-0.5,1,0,0],[0,4,11,4],1),
              ([1.5,0.5,1,0],[0,0,-10,-14],2),
              ([2.5,-2,0,1],[0,0,11,-0.6],3)]

aReduced :: [(Row, Row, Int)]
aReduced = [([1,0,0,0],[2,4,5,6],0),
            ([-0.5,1,0,0],[0,4,11,4],1),
            ([1.5,0.5,1,0],[0,0,-10,-14],2),
            ([2.5,-2,-1.1,1],[0,0,0,-16],3)]

b :: [(Row, Row, Int)]
b = [([1,0,0],[0,3,5],0),
     ([0,1,0],[1,2,3],1),
     ([0,0,1],[2,16,20],2)]

bSanitised = [([1,0,0],[2,16,20],2),
              ([0,1,0],[1,2,3],1),
              ([0,0,1],[0,3,5],0)]


p, q, r :: Matrix
p = [[1,0,0],
     [0,1,0],
     [0,0,1]]
q = [[2,16,20],
     [1,2,3],
     [0,3,5]]
r = [[0,0,1],
     [0,1,0],
     [1,0,0]]


-- swaps zero pivot with LARGEST available pivot
bReduced :: [(Row, Row, Int)]
bReduced = [([1,0,0],[2,16,20],2),
     ([0.5,1,0],[0,-6,-7],1),
     ([0,-0.5,1],[0,0,1.5],0)]



weavedEqual :: [(Row, Row, Int)] -> [(Row, Row, Int)] -> Bool
weavedEqual xs ys = all (== True) (zipWith tupEq xs ys)
    where tupEq = (\(as,bs,i) (ps,qs,j) -> (equalsRR as ps) && (equalsRR bs qs) && i==j)



checkDecomp :: Matrix -> Bool
checkDecomp a = equalsMM (multMM perm a) (multMM lower upper)
    where (lower, upper, perm) = fromRight (p,p,p) (decompLUP a)

decompExample :: Matrix
decompExample = [[0,3,2,1],
                 [0,0.01,1,0],
                 [3.1,1,7,2],
                 [12532,10,1,2]]

spec :: Spec
spec = describe "All" $ do
        describe "OneRow" $ do
            it "works with zero multiplier" $ do
                oneRow [1,2,3] [4,5,6] 0  `shouldBe` [1,2,3]
            it "works on standard input"  $ do
                oneRow [3,5,7] [8,11,2] 1.5 `shouldBe` [-9,-11.5,4]

        describe "InnerLoop" $ do
            it "works on a random example" $ do
                weavedEqual (innerLoop 0 a) aZeroResult `shouldBe` True
            it "works a second time on a random example" $ do
                weavedEqual (innerLoop 1 (tail aZeroResult)) aOneResult `shouldBe` True
        
        describe "OuterLoop" $ do
            it "works on a standard example" $ do
                weavedEqual (fromRight [([],[],0)] $ outerLoop a) aReduced `shouldBe` True
            it "works on a zero pivot example" $ do
                weavedEqual (fromRight [([],[],0)] $ outerLoop b) bReduced

        describe "Sanitise" $ do
            it "works on a random example" $ do
                weavedEqual (fromMaybe [([],[],0)] $ sanitise b 0) bSanitised `shouldBe` True

        describe "UnWeaver" $ do
            it "works on a random example" $ do
                let (lower,upper,perm) = unWeaver bSanitised
                (equalsMM lower p) && (equalsMM upper q) && (equalsMM perm r) `shouldBe` True

        describe "DecompLUP" $ do
            it "works on a standard example" $ do
                checkDecomp q `shouldBe` True
            
            it "works with zero pivots" $ do
                checkDecomp decompExample `shouldBe` True
            
            it "throws an error on a singular matrix" $ do
                decompLUP [[1,3],[2,6]] `shouldBe` (Left "no non-zero pivots")

        describe "Solve"  $ do
            it "works on a standard example" $ do
                let (Right x) = solve [[1,1,1],[0,2,5],[2,5,-1]] [6,-4,27]
                equalsVV x [5,3,-2] `shouldBe` True

            it "works on another standard example" $ do
                let (Right x) = solve [[1,1,1,1],[2,3,0,-1],[-3,4,1,2],[1,2,-1,1]] [13,-1,10,1]
                equalsVV x [2,0,6,5] `shouldBe` True

            it "works if we rearrange previous example" $ do
                let (Right x) = solve [[0,3,2,-1],[1,1,1,1],[1,4,-3,2],[-1,2,1,1]] [-1,13,10,1]
                equalsVV x [6,0,2,5] `shouldBe` True

            it "works on a zero pivot example" $ do  
                let (Right x) = solve [[0,1,-4],[3,2,-2],[-1,2,-6]] [1,8,4]
                equalsVV x [0,5,1] `shouldBe` True

            it "doesn't work on a singular system" $ do
                solve [[0,0],[0,1]] [4,5] `shouldBe` (Left "No solution")
            


main :: IO ()
main = hspec spec