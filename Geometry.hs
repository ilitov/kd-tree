module Geometry
(
    Point(..)
,   PointCoordType
,   Point2D(..)
,   Point3D(..)
,   KDTree(..)
,   Region
,   generateListOfSortedPoints
,   makeRegion
,   median
,   euclideanDist2
,   regionContainsPoint
,   regionContainsRegion
,   regionIntersectRegion
,   addPointToRegion
,   calculateRegion
,   leftRegion
,   rightRegion
) where

import qualified Data.List as List

type PointCoordType = Float

class Point p where
    dim :: p -> Int
    coord :: p -> Int -> PointCoordType
    
    setCoord :: p -> Int -> PointCoordType -> p

    minPoint :: p -> p -> p
    maxPoint :: p -> p -> p

data Point2D = Point2D { x2d :: PointCoordType, y2d :: PointCoordType } deriving (Eq, Show, Read)

data Point3D = Point3D { x3d :: PointCoordType, y3d :: PointCoordType, z3d :: PointCoordType } deriving (Eq, Show, Read)

instance Point Point2D where
    dim _ = 2

    coord p 0 = x2d p
    coord p 1 = y2d p

    setCoord p 0 val = Point2D val (y2d p)
    setCoord p 1 val = Point2D (x2d p) val
    
    minPoint p1 p2 = Point2D (min (x2d p1) (x2d p2)) (min (y2d p1) (y2d p2))
    maxPoint p1 p2 = Point2D (max (x2d p1) (x2d p2)) (max (y2d p1) (y2d p2))

instance Point Point3D where
    dim _ = 3

    coord p 0 = x3d p
    coord p 1 = y3d p
    coord p 2 = z3d p

    setCoord p 0 val = Point3D val (y3d p) (z3d p)
    setCoord p 1 val = Point3D (x3d p) val (z3d p)
    setCoord p 2 val = Point3D (x3d p) (y3d p) val

    minPoint p1 p2 = Point3D (min (x3d p1) (x3d p2)) (min (y3d p1) (y3d p2)) (min (z3d p1) (z3d p2))
    maxPoint p1 p2 = Point3D (max (x3d p1) (x3d p2)) (max (y3d p1) (y3d p2)) (max (z3d p1) (z3d p2))

comparePointsByCoord :: (Point p) => Int -> p -> p -> Ordering
comparePointsByCoord c p1 p2 = compare (coord p1 c) (coord p2 c)

euclideanDist2 :: (Point p) => p -> p -> PointCoordType
euclideanDist2 p1 p2 = sum $ map (\c -> (coord p1 c - coord p2 c) ^ 2) [0..dim p1 - 1]

generateListOfSortedPoints :: (Point p) => [p] -> Int -> [[p]]
generateListOfSortedPoints [] _             = []
generateListOfSortedPoints points maxCoord  = map (\x -> List.sortBy (comparePointsByCoord x) points) [0..maxCoord - 1]

medianOfThree :: (Point p) => [p] -> Int -> Int -> p
medianOfThree points axis size =    if size < 2
                                    then head points
                                    else    let p1 = points !! 0
                                                p2 = points !! (size `div` 2)
                                                p3 = points !! (size - 1)
                                                sorted = List.sortBy (comparePointsByCoord axis) [p1, p2, p3] 
                                            in  sorted !! 1

middlePoint :: (Point p) => [p] -> Int -> Int -> p
middlePoint points axis k
        |   lengthSmaller == k  = pivot
        |   lengthSmaller < k   = middlePoint bigger axis (k - lengthSmaller)
        |   otherwise           = middlePoint (take (lengthSmaller - 1) smaller) axis k
            where   pivot = medianOfThree points axis size
                    smaller = filter (\p -> coord p axis <= coord pivot axis) points
                    bigger = filter (\p -> coord p axis > coord pivot axis) points
                    lengthSmaller = length smaller
                    size = length points

median :: (Point p) => [p] -> Int -> p
median points axis = middlePoint points axis k
                    where   k = size `div` 2
                            size = length points

data KDTree p = LeafKDTree { point :: [p] } |
                NodeKDTree {    axis :: Int
                            ,   pivot :: PointCoordType
                            ,   leftTree :: KDTree p
                            ,   rightTree :: KDTree p } deriving (Eq, Show, Read)

data Region p = Region {    minValue :: p
                        ,   maxValue :: p} deriving (Eq, Show)

makeRegion :: (Point p) => p -> p -> Region p
makeRegion p1 p2 = Region (minPoint p1 p2) (maxPoint p1 p2)

regionExtremeValues :: (Point p) => Region p -> Int -> (PointCoordType, PointCoordType)
regionExtremeValues region axis = (min val1 val2, max val1 val2)
                                    where   val1 = coord (minValue region) axis
                                            val2 = coord (maxValue region) axis

regionContainsPoint :: (Point p) => Region p -> p -> Bool
regionContainsPoint region point = 
    List.all contains [0..dim point - 1]
        where contains axis =   let (minVal, maxVal) = regionExtremeValues region axis
                                in  coord point axis >= minVal && coord point axis <= maxVal

regionIntersectRegion :: (Point p) => Region p -> Region p -> Bool
regionIntersectRegion region1 region2 =
    List.all intersect [0..dimension - 1]
        where   dimension = dim (minValue region1)
                intersect axis =    let (minVal1, maxVal1) = regionExtremeValues region1 axis
                                        (minVal2, maxVal2) = regionExtremeValues region2 axis
                                    in  minVal2 >= minVal1 && minVal2 <= maxVal1 ||
                                        maxVal2 >= minVal1 && maxVal2 <= maxVal1

regionContainsRegion :: (Point p) => Region p -> Region p -> Bool
regionContainsRegion region1 region2 = 
    List.all intersect [0..dimension - 1]
        where   dimension = dim (minValue region1)
                intersect axis =    let (minVal1, maxVal1) = regionExtremeValues region1 axis
                                        (minVal2, maxVal2) = regionExtremeValues region2 axis
                                    in  minVal2 >= minVal1 && minVal2 <= maxVal1 &&
                                        maxVal2 >= minVal1 && maxVal2 <= maxVal1

addPointToRegion :: (Point p) => Region p -> p -> Region p
addPointToRegion region point = Region (minPoint point (minValue region)) (maxPoint point (maxValue region))

mergeRegions :: (Point p, Eq p) => Region p -> Region p -> Region p
mergeRegions r1 r2 = Region (minPoint (minValue r1) (minValue r2)) (maxPoint (maxValue r1) (maxValue r2))

calculateRegion :: (Point p, Eq p) => KDTree p -> Maybe (Region p)
calculateRegion (LeafKDTree []) = Nothing
calculateRegion (LeafKDTree points) = Just (Region (getExtreme minPoint) (getExtreme maxPoint))
    where getExtreme f = foldr1 f points

calculateRegion (NodeKDTree _ _ leftTree rightTree) = 
    let resultLeft = calculateRegion leftTree
        resultRight = calculateRegion rightTree
    in  if resultLeft /= Nothing && resultRight /= Nothing
        then    let Just resLeft = resultLeft
                    Just resRight = resultRight
                in  Just $ mergeRegions resLeft resRight
        else    if resultLeft /= Nothing && resultRight == Nothing
                then resultLeft
        else    if resultLeft == Nothing && resultRight /= Nothing
                then resultRight
        else Nothing

leftRegion :: (Point p) => Region p -> Int -> PointCoordType -> Region p
leftRegion r axis val = Region (minValue r) (minPoint pivotAxis (maxValue r))
    where pivotAxis = setCoord (maxValue r) axis val

rightRegion :: (Point p) => Region p -> Int -> PointCoordType -> Region p
rightRegion r axis val = Region (maxPoint pivotAxis (minValue r)) (maxValue r)
    where pivotAxis = setCoord (minValue r) axis val
