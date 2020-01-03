import Geometry
import qualified Data.List as List
import qualified System.IO as IO

points = [Point2D 7 2, Point2D 5 4, Point2D 2 3, Point2D 4 7, Point2D 9 6, Point2D 8 1]

rp :: (Point p, Read p) => [Char] -> KDTree p
rp str = read str

maximumPointsInLeaf :: Int
maximumPointsInLeaf = 1

-- BUILD KD-TREE FROM LIST OF POINTS

makeKDTreeFromList :: (Point p, Eq p) => [p] -> KDTree p
makeKDTreeFromList listOfPoints = makeKDTreeFromListInner listOfPoints 0

makeKDTreeFromListInner :: (Point p, Eq p) => [p] -> Int -> KDTree p
makeKDTreeFromListInner [] _    = LeafKDTree []
makeKDTreeFromListInner points axisValue =
    if length points <= maximumPointsInLeaf
    then LeafKDTree points
    else NodeKDTree {   axis = currAxis
                    ,   pivot = midPointCoordValue
                    ,   leftTree = leftSubtree
                    ,   rightTree = rightSubtree }

     where  currAxis = mod axisValue (dim $ head points)
            midPoint = median points currAxis
            midPointCoordValue = coord midPoint currAxis
            smaller = filter (\p -> coord p currAxis <= midPointCoordValue) points
            bigger = filter (\p -> coord p currAxis > midPointCoordValue) points
            leftSubtree = makeKDTreeFromListInner smaller (axisValue + 1)
            rightSubtree = makeKDTreeFromListInner bigger (axisValue + 1)

-- BUILD KD-TREE FROM LIST OF POINTS - PRESORTED

makeKDTreeFromListPS :: (Point p) => [p] -> KDTree p
makeKDTreeFromListPS []         = LeafKDTree []
makeKDTreeFromListPS pts = makeKDTreeFromListPSInner 0 sorted
    where   sorted = generateListOfSortedPoints pts (dim $ head pts)

makeKDTreeFromListPSInner :: (Point p) => Int -> [[p]] -> KDTree p
makeKDTreeFromListPSInner axisValue sorted = 
    if ptsLength <= maximumPointsInLeaf
    then LeafKDTree pts
    else NodeKDTree {   axis = currAxis
                    ,   pivot = midPointCoordValue
                    ,   leftTree = leftSubtree
                    ,   rightTree = rightSubtree }
                    
    where   dimension = dim $ head $ head sorted
            currAxis = axisValue `mod` dimension
            pts = sorted !! currAxis
            ptsLength = length pts
            mid = ptsLength `div` 2 - 1
            midPoint = pts !! mid
            midPointCoordValue = coord midPoint currAxis

            filterSmaller pts = filter (\p -> coord p currAxis <= midPointCoordValue) pts
            filterBigger pts = filter (\p -> coord p currAxis > midPointCoordValue) pts

            smaller = map (\d -> filterSmaller . sorted !! d) [0..dimension - 1]
            bigger = map (\d -> filterBigger . sorted !! d) [0..dimension - 1]

            leftSubtree = makeKDTreeFromListPSInner (axisValue + 1) smaller
            rightSubtree = makeKDTreeFromListPSInner (axisValue + 1) bigger



-- SEARCH NEAREST NEIGHBOUR

nearestNeighbour :: (Point p, Eq p) => KDTree p -> (p -> p -> PointCoordType) -> p -> Maybe p
nearestNeighbour (LeafKDTree []) _ _         = Nothing
nearestNeighbour (LeafKDTree points) metric point = Just $ List.minimumBy (\x y -> metric point x `compare` metric point y) points
nearestNeighbour (NodeKDTree axis pivot leftTree rightTree) metric point =
    if pivot < currCoordPoint
    then nearestNeighbourInner rightTree leftTree
    else nearestNeighbourInner leftTree rightTree
     where  currCoordPoint = coord point axis

            betterCandidate :: (Point p, Eq p) => (p -> PointCoordType) -> p -> p -> p
            betterCandidate metric p1 p2 = List.minimumBy (\x y -> metric x `compare` metric y) [p1, p2]

            nearestNeighbourInner firstTree secondTree =
                let firstPartResult = nearestNeighbour firstTree metric point
                    secondPartResult =  if  firstPartResult == Nothing || 
                                            let Just firstRes = firstPartResult in 
                                                (pivot - currCoordPoint) ^ 2 < metric point firstRes
                                        then nearestNeighbour secondTree metric point
                                        else Nothing
                in  if firstPartResult /= Nothing && secondPartResult /= Nothing
                    then    let Just first = firstPartResult
                                Just second = secondPartResult
                            in Just $ betterCandidate (metric point) first second
                    else    if firstPartResult /= Nothing && secondPartResult == Nothing
                            then firstPartResult
                    else    if firstPartResult == Nothing && secondPartResult /= Nothing
                            then secondPartResult
                    else Nothing

--  SEARCH IN INTERVAL

-- Returns all points contained in that tree (the leaves of the tree)
reportSubtree :: (Point p) => KDTree p -> [p]
reportSubtree (LeafKDTree points) = points
reportSubtree (NodeKDTree _ _ leftTree rightTree) = reportSubtree leftTree ++ reportSubtree rightTree

intervalSearchInner :: (Point p) => KDTree p -> Region p -> Region p -> [p]
intervalSearchInner (LeafKDTree []) _ _ = []
intervalSearchInner (LeafKDTree points) _ region = filter (regionContainsPoint region) points
                             
intervalSearchInner (NodeKDTree axis pivot leftTree rightTree) treeRegion region =
    let leftR = leftRegion treeRegion axis pivot
        rightR = rightRegion treeRegion axis pivot
    in (findPoints leftR leftTree) ++ (findPoints rightR rightTree)
        where findPoints r tree =   if regionContainsRegion region r
                                    then reportSubtree tree
                                    else    if regionIntersectRegion treeRegion r
                                            then intervalSearchInner tree r region
                                    else []

intervalSearch :: (Point p, Eq p) => KDTree p -> Region p -> [p]
intervalSearch tree region =    if treeRegion == Nothing
                                then []
                                else intervalSearchInner tree r region
                                 where treeRegion = calculateRegion tree
                                       Just r = treeRegion

-- SAVE/READ TREE TO/FROM FILE

saveKDTree :: (Point p, Show p) => String -> KDTree p -> IO ()
saveKDTree filename tree = writeFile filename $ show tree

readKDTree :: (Point p, Read p) => String -> IO (KDTree p)
readKDTree filename = do
        contents <- readFile filename
        return (read contents)

-- OTHER

nearestNeighbourIO :: (Point p, Eq p) => IO (KDTree p) -> (p -> p -> PointCoordType) -> p -> IO (Maybe p)
nearestNeighbourIO tree metric p = do
                tmp <- tree
                return (nearestNeighbour tmp metric p)

intervalSearchIO :: (Point p, Eq p) => IO (KDTree p) -> Region p -> IO ([p])
intervalSearchIO tree region = do
                tmp <- tree
                return (intervalSearch tmp region)

