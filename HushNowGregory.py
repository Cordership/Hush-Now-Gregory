import numpy as np, copy, pygame, re, collections.abc

# The relative x and y coordinates for the moore neighborhood
mooreX = [-1, 0, 1, -1, 1, -1, 0, 1]
mooreY = [-1, -1, -1, 0, 0, 1, 1, 1]

# Transition table for 23/3/3
transitionTable=[
    [0,0,0,1,0,0,0,0,0],
    [2,2,1,1,2,2,2,2,2],
    [0,0,0,0,0,0,0,0,0]
]

# Relative coordinates of the cells of a parvoship to it's position.
parvoBodyCoords = [
    [  # phase 0
        [  # state 1
            [2, 0], [1, 1], [2, 1], [3, 1], [0, 2], [1, 2], [4, 2], [2, 3], [3, 3], [4, 3], [5, 3], [3, 4], [4, 4]
        ],
        [  # state 2
            [1, 3]
        ]
    ],
    [  # phase 1
        [
            [1, 0], [2, 0], [3, 0], [0, 1], [3, 1], [0, 2], [5, 2], [2, 3], [5, 3], [2, 4], [5, 4]
        ],
        [
            [1, 1], [2, 1], [1, 2], [4, 2], [3, 3], [4, 3], [3, 4], [4, 4]
        ]
    ],
    [  # phase 2
        [
            [2, 0], [1, 1], [2, 1], [3, 1], [0, 2], [3, 2], [4, 2], [1, 4], [5, 4], [6, 4]
        ],
        [
            [0, 3], [5, 3], [2, 4], [2, 5], [5, 5]
        ]
    ],
    [  # phase 3
        [
            [1, 0], [2, 0], [3, 0], [1, 1], [4, 1], [1, 2], [3, 2], [4, 2], [4, 3]
        ],
        [
            [2, 1], [3, 1], [0, 2], [1, 4], [5, 4], [6, 4]
        ]
    ]
]

# Arrays containing each phase of the parvoship
parvoBodyArray = [
    [  # phase 0
        [0, 0, 1, 0, 0],
        [0, 1, 1, 2, 0],
        [1, 1, 0, 1, 0],
        [0, 1, 0, 1, 1],
        [0, 0, 1, 1, 1],
        [0, 0, 0, 1, 0]
    ],
    [  # phase 1
        [0, 1, 1, 0, 0],
        [1, 2, 2, 0, 0],
        [1, 2, 0, 1, 1],
        [1, 1, 0, 2, 2],
        [0, 0, 2, 2, 2],
        [0, 0, 1, 1, 1]
    ],
    [  # phase 2
        [0, 0, 1, 2, 0, 0],
        [0, 1, 0, 0, 1, 0],
        [1, 1, 0, 0, 2, 2],
        [0, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 2, 1, 2],
        [0, 0, 0, 0, 1, 0]
    ],
    [  # phase 3
        [0, 0, 2, 0, 0],
        [1, 1, 1, 0, 2],
        [1, 2, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [0, 1, 1, 1, 0],
        [0, 0, 0, 0, 2],
        [0, 0, 0, 0, 2]
    ]
]

# The widths and heights of each phase of the parvoship
parvoWidth = [6, 6, 7, 7]
parvoHeight = [5, 5, 6, 5]

# Start (inclusive) and end (non-inclusive) points for each row of the spacing shadow cast by each phase of the parvoship
parvoInterspacing = [
    [  # same chirality
        [  # phase difference of 0
            [-5, 3], [-6, 4], [-7, 5], [-7, 6], [-7, 7], [-7, 8], [-7, 8], [-7, 8], [-6, 8], [-5, 8], [-4, 8], [-3, 7],
            [-2, 6]
        ],
        [  # phase difference of 1
            [-6, 4], [-7, 5], [-8, 5], [-8, 7], [-8, 7], [-8, 8], [-7, 8], [-6, 8], [-5, 8], [-5, 8], [-4, 8], [-3, 7]
        ],
        [  # phase difference of 2
            [-5, 3], [-6, 4], [-7, 5], [-7, 6], [-7, 7], [-7, 7], [-7, 8], [-6, 8], [-6, 8], [-5, 8], [-4, 8], [-3, 7],
            [-2, 6]
        ],
        [  # phase difference of 3
            [-6, 4], [-7, 5], [-7, 6], [-7, 6], [-7, 7], [-7, 8], [-7, 9], [-6, 9], [-6, 9], [-4, 9], [-4, 8], [-3, 7]
        ],
    ],
    [  # opposite chirality
        [
            [2, 10], [1, 11], [0, 12], [-1, 13], [-2, 14], [-2, 14], [-1, 14], [-2, 14], [-2, 14], [-1, 13], [0, 12],
            [1, 11], [2, 10]
        ],
        [
            [1, 11], [0, 12], [-1, 13], [-1, 14], [-2, 14], [-2, 14], [-2, 14], [-2, 14], [-1, 13], [0, 13], [0, 12],
            [1, 11]
        ],
        [
            [2, 10], [1, 11], [0, 12], [-1, 13], [-1, 13], [-1, 13], [-1, 13], [-1, 13], [-1, 13], [-1, 13], [0, 12],
            [1, 11], [2, 10]
        ],
        [
            [1, 11], [0, 12], [0, 13], [-1, 13], [-2, 14], [-2, 14], [-2, 14], [-2, 14], [-1, 14], [-1, 13], [0, 12],
            [1, 11]
        ],
    ]
]

# Left and right bounds for the parvoship spacing shadow
parvoInterspacingXBounds = [
    [  # same chirality
        [-7, 7],
        [-8, 7],
        [-7, 7],
        [-7, 8],
    ],
    [  # opposite chirality
        [-2, 13],
        [-2, 13],
        [-1, 12],
        [-2, 13],
    ]
]

# Start (inclusive) and end (non-inclusive) points for each row of the fatal interaction shadow cast by each living cell of the puff pattern
parvoCellspacing = [
    [  # phase 0
        [0, 5], [-1, 6], [-1, 7], [-1, 7], [0, 7], [0, 7], [0, 7]
    ],
    [  # phase 1
        [0, 5], [-1, 6], [-1, 6], [-1, 7], [-1, 7], [-1, 7], [0, 7], [1, 7]
    ],
    [  # phase 2
        [0, 5], [-1, 6], [-1, 6], [-1, 7], [-1, 7], [-1, 7], [0, 7], [1, 7]
    ],
    [  # phase 3
        [0, 5], [-1, 6], [-1, 7], [-1, 7], [-1, 7], [-1, 8], [0, 8], [0, 8]
    ]
]

# Coordinates of the search shadow cast by each living cell of the puff pattern (phase 3 doesn't have any safe interactions with the puff)
parvoSearchCoords = [
    [  # phase 0
        [-2, 1], [-2, 2], [7, 2], [-2, 3], [7, 4]
    ],
    [  # phase 1
        [-2, 1], [7, 2], [7, 4]
    ],
    [  # phase 2
        [8, 3], [7, 4], [8, 4], [-1, 5], [8, 5], [0, 6]
    ]
]

# Coordinates of the post-search spacing shadow cast by each living cell of the puff pattern (phase 3 doesn't have a post-search spacing shadow)
# These are all coordinates in the strict interaction shadow that are not in the fatal interaction shadow or the search shadow
parvoPostspacingCoords = [
    [  # phase 0
        [-1, 3], [-1, 4], [2, 6], [3, 6], [4, 6], [5, 6]
    ],
    [  # phase 1
        [-2, 0], [7, 1], [-2, 2], [-2, 3], [7, 3], [0, 5], [7, 5]
    ],
    [  # phase 2
        [7, 2], [7, 3], [7, 5], [7, 6]
    ]
]


# Parvoship class
class Parvo:
    phase = 0
    x = 0
    y = 0
    leftHanded = True

    def __init__(self, phase, x, y, leftHanded):
        self.phase = phase
        self.x = x
        self.y = y
        self.leftHanded = leftHanded

    def width(self):
        return parvoWidth[self.phase]

    def height(self):
        return parvoHeight[self.phase]

    # Y bound in negitive direction (inclusive)
    def top(self):
        return self.y

    # Y bound in positive direction (exclusive)
    def bottom(self):
        return self.y + parvoHeight[self.phase]

    # X bound in negitive direction (inclusive)
    def left(self):
        if self.leftHanded:
            return self.x
        return self.x - parvoWidth[self.phase] + 1

    # X bound in positive direction (exclusive)
    def right(self):
        if self.leftHanded:
            return self.x + parvoWidth[self.phase]
        return self.x + 1

    # Coordinates of the cell at the head of the parvoship
    def leadingEdge(self):
        return [self.x + 2 if self.leftHanded else self.x - 2, self.y]

    # Returns the state of the cell in the parvoship at (x,y) if (x,y) is in the bounds of the parvoship, otherwise returns false
    def cellAt(self, x, y):
        if self.left() <= x and x < self.right() and self.top() <= y and y < self.bottom():
            if self.leftHanded:
                return parvoBodyArray[self.phase][x - self.x][y - self.y]
            return parvoBodyArray[self.phase][self.x - x][y - self.y]
        return False

    # Advances the parvoship 1 generation
    def advance(self):
        self.phase = (self.phase + 1) % 4
        if self.phase % 2 == 0:
            self.y -= 1


# Class for simulating flotillas
class PatternGrid:
    array = np.array([])
    parvos = []
    xStart = 0
    yStart = 0

    # Array and Parvos must be duplicated before being inserted as arguments
    def __init__(self, array, parvos, x, y):
        self.array = array
        self.parvos = parvos
        # x and y coordinates of the top left corner of the grid
        self.xStart = x
        self.yStart = y
        # Draw all parvoships in parvos onto the array
        for parvo in parvos:
            if self.overlaps(parvo):
                writeXStart = max(self.xStart, parvo.left())
                writeXEnd = min(self.xEnd(), parvo.right())
                writeYStart = max(self.yStart, parvo.top())
                writeYEnd = min(self.yEnd(), parvo.bottom())
                for y in range(writeYStart, writeYEnd):
                    for x in range(writeXStart, writeXEnd):
                        state = parvo.cellAt(x, y)
                        if state:
                            self[x, y] = state

    # Get and Set functions that take into account xStart and yStart
    def __getitem__(self, key):
        return self.array[key[0] - self.xStart][key[1] - self.yStart]

    def __setitem__(self, key, value):
        self.array[key[0] - self.xStart][key[1] - self.yStart] = value

    # Checks if the overlapping portions of the two PatternGrids are the identical and that the non-overlapping portions are empty
    def __eq__(self, other):
        for y in range(self.yStart + 1, self.yEnd() - 1):
            for x in range(self.xStart + 1, self.xEnd() - 1):
                if x <= other.xStart or other.xEnd() - 1 <= x or y <= other.yStart or other.yEnd() - 1 <= y:
                    if self[x, y] != 0:
                        return False
                else:
                    if self[x, y] != other[x, y]:
                        return False
        for y in range(other.yStart + 1, other.yEnd() - 1):
            for x in range(other.xStart + 1, other.xEnd() - 1):
                if (x <= self.xStart or self.xEnd() - 1 <= x or y <= self.yStart or self.yEnd() - 1 <= y) and other[x, y] != 0:
                    return False
        return True

    # Checks if the two PatternGrids have the same coordinates, dimensions, and contents
    def strictlyEquals(self, other):
        if self.xStart != other.xStart or self.yStart != other.yStart:
            return False
        return np.array_equal(self.array, other.array)

    # Resizes grid to fit (x,y), filling with dead cells and ignoring parvoships
    def fit(self, x, y):
        if x < self.xStart:
            breadth = self.xStart - x
            edge = np.zeros([breadth, self.height()], dtype=np.int64)
            self.array = np.concatenate((edge, self.array), axis=0)
            self.xStart -= breadth
        elif self.xEnd() <= x:
            edge = np.zeros([x - self.xEnd() + 1, self.height()], dtype=np.int64)
            self.array = np.concatenate((self.array, edge), axis=0)
        if y < self.yStart:
            breadth = self.yStart - y
            edge = np.zeros([self.width(), breadth], dtype=np.int64)
            self.array = np.concatenate((edge, self.array), axis=1)
            self.yStart -= breadth
        elif self.yEnd() <= y:
            edge = np.zeros([self.width(), y - self.yEnd() + 1], dtype=np.int64)
            self.array = np.concatenate((self.array, edge), axis=1)

    def width(self):
        return self.array.shape[0]

    def height(self):
        return self.array.shape[1]

    # x and y bounds in the positive direction (exclusive)
    def xEnd(self):
        return self.xStart + self.array.shape[0]

    def yEnd(self):
        return self.yStart + self.array.shape[1]

    # Returns true if bounds of parvoship overlap bounds of PatternGrid
    def overlaps(self, parvo):
        return (parvo.left() < self.xEnd() or self.xStart < parvo.right()) and (
                parvo.top() < self.yEnd() or self.yStart < parvo.bottom())

    # Advances the PatternGrid 1 generation
    # Bounds of PatternGrid will always be two cells outside of the difference zone between the pattern and version of it containing only parvoships
    # This means that if a parvoship is destroyed, the bounds of the PatternGrid will follow the absence of a parvoship
    def advance(self):
        # newArray has edges cut off, seeing as the outermost edges of the PatternGrid can not be and are not simulated
        # This will not be a problem, seeing as the bounds of the grid expand so that all cells that could be affected by the difference zone are 1 cell inside the bounds
        # Thus, the edges will never need to be simulated, only referenced from the parvoship only version of the pattern to determine the state of the cells inside of the edges
        newArray = np.zeros([self.width() - 2, self.height() - 2], dtype=np.int64)
        for y in range(self.height()):
            for x in range(self.width()):
                if self.array[x][y]==1:
                    for n in range(8):
                        xN = x + mooreX[n] - 1
                        yN = y + mooreY[n] - 1
                        if 0<=xN and xN<self.width()-2 and 0<=yN and yN<self.height()-2:
                            newArray[xN][yN]+=1
        for y in range(self.height() - 2):
            for x in range(self.width() - 2):
                newArray[x][y]=transitionTable[self.array[x+1][y+1]][newArray[x][y]]
        # newArray now needs to be resized to fit the difference zone
        for parvo in self.parvos:
            parvo.advance()
        xMin = self.xEnd() - 1
        xMax = self.xStart
        yMin = self.yEnd() - 1
        yMax = self.yStart
        # Mark each cell that matches a non dead cell of one of the PatternGrid's parvoships by adding two to the state
        # If a non dead cell in one of the parvoships does not match the same cell in newArray, count it as part of the difference zone
        for parvo in self.parvos:
            if (parvo.left() < self.xEnd() - 1 or self.xStart + 1 < parvo.right()) and \
                    (parvo.top() < self.yEnd() - 1 or self.yStart + 1 < parvo.bottom()):
                for state in range(1, 3):
                    for coords in parvoBodyCoords[parvo.phase][state - 1]:
                        x = parvo.x + coords[0] if (parvo.leftHanded) else parvo.x - coords[0]
                        y = parvo.y + coords[1]
                        xArr = x - self.xStart - 1
                        yArr = y - self.yStart - 1
                        if 0 <= xArr and xArr < newArray.shape[0] and 0 <= yArr and yArr < newArray.shape[1]:
                            if newArray[xArr][yArr] == state:
                                newArray[xArr][yArr] += 2
                            else:
                                if x < xMin:
                                    xMin = x
                                if xMax < x:
                                    xMax = x
                                if y < yMin:
                                    yMin = y
                                if yMax < y:
                                    yMax = y
        # Scan for non dead cells in newArray that have not been marked; count them in the difference zone
        x = newArray.shape[0]
        scanning = True
        while scanning and xMax - self.xStart <= (x := x - 1):
            for y in range(newArray.shape[1]):
                if newArray[x][y] == 1 or newArray[x][y] == 2:
                    scanning = False
                    break
        xMax = self.xStart + x + 1
        # Calculate the width/height of the buffer that must be added on if some of the bounds of difference zone expanded or stayed the same
        edgeRightBreadth = xMax - self.xEnd() + 4
        # Crop newArray if the bounds of the difference zone shrunk by more than one cell, seeing as newArray has already been cropped by 1 cell compared to the current array
        if edgeRightBreadth < 0:
            newArray = newArray[:x + 3]
        # If the difference zone has no width, then the puff is dead.
        if x == -1:
            self.array = None
            self.xStart = 0
            self.yStart = 0
            return False
        # Same as above
        y = newArray.shape[1]
        scanning = True
        while scanning and yMax - self.yStart <= (y := y - 1):
            for x in range(newArray.shape[0]):
                if newArray[x][y] == 1 or newArray[x][y] == 2:
                    scanning = False
                    break
        yMax = self.yStart + y + 1
        edgeBottomBreadth = yMax - self.yEnd() + 4
        if edgeBottomBreadth < 0:
            newArray = newArray[:, :y + 3]
        # Same as above
        x = -1
        scanning = True
        while scanning and (x := x + 1) < xMin - self.xStart - 1:
            for y in range(newArray.shape[1]):
                if newArray[x][y] == 1 or newArray[x][y] == 2:
                    scanning = False
                    break
        xMin = self.xStart + x + 1
        edgeLeftBreadth = 2 - x
        if edgeLeftBreadth < 0:
            newArray = newArray[x - 2:]
        # Same as above
        y = -1
        scanning = True
        while scanning and (y := y + 1) < yMin - self.yStart - 1:
            for x in range(newArray.shape[0]):
                if newArray[x][y] == 1 or newArray[x][y] == 2:
                    scanning = False
                    break
        yMin = self.yStart + y + 1
        edgeTopBreadth = 2 - y
        if edgeTopBreadth < 0:
            newArray = newArray[:, y - 2:]
        # If bounds of difference zone have grown or stayed the same, add buffer referenced from parvoship only version of pattern
        if edgeLeftBreadth > 0:
            edge = PatternGrid(np.zeros((edgeLeftBreadth, newArray.shape[1]), dtype=np.int64),
                               self.parvos, xMin - 2, max(self.yStart + 1, yMin - 2))
            newArray = np.concatenate((edge.array, newArray), axis=0)
        if edgeTopBreadth > 0:
            edge = PatternGrid(np.zeros((newArray.shape[0], edgeTopBreadth), dtype=np.int64),
                               self.parvos, xMin - 2, yMin - 2)
            newArray = np.concatenate((edge.array, newArray), axis=1)
        if edgeRightBreadth > 0:
            edge = PatternGrid(np.zeros((edgeRightBreadth, newArray.shape[1]), dtype=np.int64),
                               self.parvos, self.xEnd() - 1, yMin - 2)
            newArray = np.concatenate((newArray, edge.array), axis=0)
        if edgeBottomBreadth > 0:
            edge = PatternGrid(np.zeros((newArray.shape[0], edgeBottomBreadth), dtype=np.int64),
                               self.parvos, xMin - 2, self.yEnd() - 1)
            newArray = np.concatenate((newArray, edge.array), axis=1)
        # Adjust start coordinates
        self.xStart = xMin - 2
        self.yStart = yMin - 2
        # Reassign array
        self.array = newArray
        # Unmark marked cells
        for parvo in self.parvos:
            if (parvo.left() < self.xEnd() - 1 or self.xStart + 1 < parvo.right()) and (
                    parvo.top() < self.yEnd() - 1 or self.yStart + 1 < parvo.bottom()):
                for state in range(2):
                    for coords in parvoBodyCoords[parvo.phase][state]:
                        x = parvo.x + coords[0] if (parvo.leftHanded) else parvo.x - coords[0]
                        y = parvo.y + coords[1]
                        if self.xStart <= x and x < self.xEnd() and self.yStart <= y and y < self.yEnd() and self[
                            x, y] > 2:
                            self[x, y] -= 2
        return True

    # Return copy of PatternGrid with parvoships removed
    def parvosRemoved(self):
        out = PatternGrid(copy.deepcopy(self.array), [], self.xStart, self.yStart)
        for parvo in self.parvos:
            if self.overlaps(parvo):
                writeXStart = max(self.xStart, parvo.left())
                writeXEnd = min(self.xEnd(), parvo.right())
                writeYStart = max(self.yStart, parvo.top())
                writeYEnd = min(self.yEnd(), parvo.bottom())
                for y in range(writeYStart, writeYEnd):
                    for x in range(writeXStart, writeXEnd):
                        if parvo.cellAt(x, y) == self[x, y]:
                            out[x, y] = 0
        return out


# Class that can be paired with PatternGrid to keep track of coordinates where parvoships cannot be added
# The first four bits of each cell correspond to the left-handed versions of phases 0-4; the next four correspond to the right-handed versions
# Bits set to 1 indicate that a parvoship of the corresponding phase and handedness cannot be placed at this coordinate at generation 0
class ShadowGrid:
    array = np.array([[0]])
    xStart = 0
    yStart = 0

    # Initialize the ShadowGrid as if the ash in startArray has already existed for several generations, given that the ash oscillates at the given period
    # This prevents parvoships from being added that would have interacted with the ash before the start of the search simulation
    # Does not insert the shadows of the parvoships into the ShadowGrid
    # startArray must be duplicated before being inserted as an argument
    def __init__(self, startArray, period):
        ashPhases = []
        ash = PatternGrid(startArray, [], -2, -2)
        ashPhases.append(copy.deepcopy(ash))
        for gen in range(period - 1):
            ash.advance()
            ashPhases.append(copy.deepcopy(ash))
        gen = -1
        # This basically runs until the bottom most cell of the newly added shadows is above the top most coordinate in parvoSearchCoords cast by the top most cell in the pattern at search generation 0
        while gen // 2 + ashPhases[gen % period].yEnd() >= -5:
            for yCell in range(ashPhases[gen % period].yStart, ashPhases[gen % period].yEnd()):
                for xCell in range(ashPhases[gen % period].xStart, ashPhases[gen % period].xEnd()):
                    if ashPhases[gen % period][xCell, yCell] == 1:
                        self.cellshadow(xCell, yCell, gen)
                        for leftHanded in [True, False]:
                            for phase in range(3):
                                phaseStart = (phase - gen) % 4
                                markShadow = (1 if leftHanded else 16) << phaseStart
                                yShiftBack = yCell + (gen + 1 - phase % 2) // 2
                                for coords in parvoSearchCoords[phase]:
                                    xShadow = xCell - coords[0] if leftHanded else xCell + coords[0]
                                    yShadow = yShiftBack - coords[1]
                                    self[xShadow, yShadow] = self[xShadow, yShadow] | markShadow
                        self.postshadow(xCell, yCell, gen)
            gen -= 1

    # Returns copy of ShadowGrid with the shadows of each parvoship in parvos added
    def parvoShadows(self, parvos):
        out = copy.deepcopy(self)
        for parvo in parvos:
            for dPhase in range(4):
                shift = parvo.phase & dPhase & 1
                phaseShadow = (parvo.phase + dPhase) % 4
                yShadowStart = parvo.y - [6, 5][dPhase % 2] - shift
                yShadowEnd = parvo.y + 7 - shift
                for mirrored in range(2):
                    out.fit(parvo.x + parvoInterspacingXBounds[mirrored][dPhase][0] if parvo.leftHanded else parvo.x - parvoInterspacingXBounds[mirrored][dPhase][0],
                            yShadowStart)
                    out.fit(parvo.x + parvoInterspacingXBounds[mirrored][dPhase][1] if parvo.leftHanded else parvo.x - parvoInterspacingXBounds[mirrored][dPhase][1],
                            yShadowEnd)
                    markShadow = [1, 16][(not parvo.leftHanded) ^ mirrored] << phaseShadow
                    row = 0
                    for yShadow in range(yShadowStart, yShadowEnd):
                        xShadowStart = parvoInterspacing[mirrored][dPhase][row][0]
                        xShadowEnd = parvoInterspacing[mirrored][dPhase][row][1]
                        for xShadow in range(parvo.x + xShadowStart if parvo.leftHanded else parvo.x - xShadowStart,
                                             parvo.x + xShadowEnd if parvo.leftHanded else parvo.x - xShadowEnd,
                                             [-1,1][parvo.leftHanded]):
                            out[xShadow, yShadow] = out[xShadow, yShadow] | markShadow
                        row += 1
        return out

    # Inserts the shadow cast onto generation 0 by a living cell at (x,y) that exists at the given gen
    def cellshadow(self, x, y, gen):
        self.fit(x - 8, y + gen // 2 - 6)
        self.fit(x + 8, y + (gen + 1) // 2 + 2)
        for leftHanded in [True, False]:
            for phase in range(4):
                phaseStart = (phase - gen) % 4
                markShadow = (1 if leftHanded else 16) << phaseStart
                yShiftBack = y + (gen + 1 - phase % 2) // 2
                row = 0
                for yShadow in range(yShiftBack + phase % 2 + 1, yShiftBack - parvoHeight[phase] - 1, -1):
                    for xShadow in range(
                            x - parvoCellspacing[phase][row][0] if leftHanded else x + parvoCellspacing[phase][row][0],
                            x - parvoCellspacing[phase][row][1] if leftHanded else x + parvoCellspacing[phase][row][1],
                            -1 if leftHanded else 1):
                        self[xShadow, yShadow] = self[xShadow, yShadow] | markShadow
                    row += 1

    # Inserts the post-search shadow cast onto generation 0 by a living cell at (x,y) that exists at the given gen
    def postshadow(self, x, y, gen):
        for leftHanded in [True, False]:
            for phase in range(3):
                phaseStart = (phase - gen) % 4
                markShadow = (1 if leftHanded else 16) << phaseStart
                yShiftBack = y + (gen + 1 - phase % 2) // 2
                for coords in parvoPostspacingCoords[phase]:
                    xShadow = x - coords[0] if leftHanded else x + coords[0]
                    yShadow = yShiftBack - coords[1]
                    self[xShadow, yShadow] = self[xShadow, yShadow] | markShadow

    # Get and Set functions that take into account xStart and yStart
    def __getitem__(self, key):
        return self.array[key[0] - self.xStart][key[1] - self.yStart]

    def __setitem__(self, key, value):
        self.array[key[0] - self.xStart][key[1] - self.yStart] = value

    # Resizes grid to fit (x,y), filling with zeros
    def fit(self, x, y):
        if x < self.xStart:
            breadth = self.xStart - x
            edge = np.zeros([breadth, self.height()], dtype=np.int64)
            self.array = np.concatenate((edge, self.array), axis=0)
            self.xStart -= breadth
        elif self.xEnd() <= x:
            edge = np.zeros([x - self.xEnd() + 1, self.height()], dtype=np.int64)
            self.array = np.concatenate((self.array, edge), axis=0)
        if y < self.yStart:
            breadth = self.yStart - y
            edge = np.zeros([self.width(), breadth], dtype=np.int64)
            self.array = np.concatenate((edge, self.array), axis=1)
            self.yStart -= breadth
        elif self.yEnd() <= y:
            edge = np.zeros([self.width(), y - self.yEnd() + 1], dtype=np.int64)
            self.array = np.concatenate((self.array, edge), axis=1)

    def width(self):
        return self.array.shape[0]

    def height(self):
        return self.array.shape[1]

    # x and y bounds in the positive direction (exclusive)
    def xEnd(self):
        return self.xStart + self.array.shape[0]

    def yEnd(self):
        return self.yStart + self.array.shape[1]


# Class for storing flotilla information
class Flotilla:
    startAsh = None
    # ShadowGrid without parvoship shadows
    startShadow = None
    parvos = []
    startGen = 0
    stopGen = 1
    parent = None

    # parvos must be duplicated before being inserted as an argument
    def __init__(self, startAsh, startShadow, parvos, startGen, stopGen, parent):
        self.startAsh = startAsh
        self.startShadow = startShadow
        self.parvos = parvos
        self.startGen = startGen
        self.stopGen = stopGen
        self.parent = parent


# Return Encoding of run
def encodeRun(match):
    match = match.group()
    return str(len(match)) + match[0]


# Return RLE of PatternGrid
def RLE(pattern):
    symbols = ['.', 'A', 'B']
    fullPattern = copy.deepcopy(pattern)
    for parvo in pattern.parvos:
        fullPattern.fit(parvo.left(), parvo.top())
        fullPattern.fit(parvo.right() - 1, parvo.bottom() - 1)
        for state in range(1, 3):
            for coords in parvoBodyCoords[parvo.phase][state - 1]:
                fullPattern[parvo.x + coords[0] if parvo.leftHanded else parvo.x - coords[0], parvo.y + coords[1]] = state
    array = fullPattern.array
    out = f"x = {array.shape[0]}, y = {array.shape[1]}, rule = 23/3/3\n"
    for y in range(array.shape[1]):
        for x in range(array.shape[0]):
            out += symbols[array[x][y]]
        out += "$"
    out = out + "!"
    out = re.sub(r"\.+\$", "$", out)
    out = re.sub(r"\$+!", "!", out)
    out = re.sub(r"(\.\.+|AA+|BB+|\$\$+)", encodeRun, out)
    return out


# Return RLEs of this flotilla and all of its ancestors
def reactionHistory(flotilla):
    if flotilla is None:
        return ""
    out = reactionHistory(flotilla.parent)
    return out + RLE(PatternGrid(copy.deepcopy(flotilla.startAsh), flotilla.parvos, -2, -2)) + "\n\n"


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Main function.
# Arguments:
#     startAshInput:
#         File containing RLE describing initial ash that the program will attempt to perturb to produce spaceships/specified stable output
#     spaceshipOutput:
#         File to output reactions that produce spaceships onto
# Keyword Arguments:
#     minDivGens:
#         Minimum number of consecutive generations a reaction must deviate based on its most recently added parvoship
#         for the addition of said parvoship to be considered a meaningful search branch
#         Defaults to 3
#     testRunLimit:
#         Number of generations after the projected first interaction between the puff and the most recently added parvoship that the program
#         will consider the reaction a spaceship producing reaction
#         Defaults to 300
#     resultAshInOut:
#         Tuple like object containing two filenames
#         The program will search for reactions that produce stable ash that matches the RLEs listed in first file
#         These reactions will be outputted in the second file
#         If one of the RLEs in the first file contains only dead cells, the program will search for a deleter flotilla
#         for the given start ash
# Example:
#     search("startAshInputRLE.txt","spaceshipOutputRLEs.txt",resultAshInOut=("resultAshInRLEs.txt","resultAshOutRLEs.txt"))
# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def search(startAshInput, spaceshipOutput, **kwargs):

    testRunLimit = kwargs.get('testRunLimit', 300)
    if type(testRunLimit) is not int:
        raise TypeError("testRunLimit must be int")

    if type(startAshInput) is not str:
        raise TypeError("startAshInput must be str")

    # Read startAshInput as RLE
    startAshRLE = re.sub(re.compile("^#.+", re.MULTILINE), "", re.sub("[^\S\n]", "", open(startAshInput, "r").read()))
    start = re.search(re.compile("[^\n].*", re.DOTALL), startAshRLE)
    if start is None:
        raise SyntaxError("file is empty", (startAshInput, None, None, None))
    header = re.match("x=[0-9]+,y=[0-9]+,rule=23/3/3\n", start.group(0))
    if header is None:
        raise SyntaxError("header is not properly formatted", (startAshInput, start.start() + 1, None, None))
    width = int(re.search("(?<=x=)[0-9]+(?=,)", header.group(0)).group(0))
    if width == 0:
        raise SyntaxError("width cannot be 0", (startAshInput, start.start() + 1, None, header.group(0)))
    height = int(re.search("(?<=y=)[0-9]+(?=,)", header.group(0)).group(0))
    if height == 0:
        raise SyntaxError("height cannot be 0", (startAshInput, start.start() + 1, None, header.group(0)))
    startPattern = np.zeros([width, height], dtype=np.int64)
    x = 0
    y = 0
    i = start.start() + header.end()
    j = startAshRLE.find("!", i)
    if j == -1:
        j = len(startAshRLE)
    line = start.start() + 2
    while i < j:
        if startAshRLE[i] == "\n":
            line += 1
            i += 1
            continue
        run = 1
        if startAshRLE[i].isdigit():
            runStr = ""
            k = i
            while startAshRLE[k].isdigit():
                runStr += startAshRLE[k]
                k += 1
                while k < len(startAshRLE) and startAshRLE[k] == "\n":
                    line += 1
                    k += 1
                if k == len(startAshRLE):
                    raise SyntaxError("file ends with dangling run", (startAshInput, line, None, None))
            run = int(runStr)
            i = k
        if startAshRLE[i] == "$":
            x = 0
            y += run
            if y > height:
                raise SyntaxError("RLE body overflows provided height", (startAshInput, line, None, None))
        else:
            if x + run > width:
                raise SyntaxError("RLE body overflows provided width", (startAshInput, line, None, None))
            if startAshRLE[i] == ".":
                x += run
            elif startAshRLE[i] == "A":
                for k in range(run):
                    startPattern[x][y] = 1
                    x += 1
            elif startAshRLE[i] == "B":
                for k in range(run):
                    startPattern[x][y] = 2
                    x += 1
            else:
                raise SyntaxError("file contains unrecognized symbol", (startAshInput, line, None, startAshRLE[i]))
        i += 1
    xMin = -1
    scanning = True
    while scanning and (xMin := xMin + 1) < startPattern.shape[0]:
        for y in range(startPattern.shape[1]):
            if startPattern[xMin][y] != 0:
                scanning = False
                break
    if xMin == startPattern.shape[0]:
        raise SyntaxError("RLE is empty", (startAshInput, start.start(), None, None))
    startPattern = startPattern[xMin:, :]
    yMin = -1
    scanning = True
    while scanning and (yMin := yMin + 1) < startPattern.shape[1]:
        for x in range(startPattern.shape[0]):
            if startPattern[x][yMin] != 0:
                scanning = False
                break
    startPattern = startPattern[:, yMin:]
    xMax = startPattern.shape[0]
    scanning = True
    while scanning and (xMax := xMax - 1) >= 0:
        for y in range(startPattern.shape[1]):
            if startPattern[xMax][y] != 0:
                scanning = False
                break
    startPattern = startPattern[:xMax + 1, :]
    yMax = startPattern.shape[1]
    scanning = True
    while scanning and (yMax := yMax - 1) >= 0:
        for x in range(startPattern.shape[0]):
            if startPattern[x][yMax] != 0:
                scanning = False
                break
    startPattern = startPattern[:, :yMax + 1]
    startPattern = np.concatenate((startPattern, np.zeros([2, startPattern.shape[1]], dtype=np.int64)), axis=0)
    startPattern = np.concatenate((startPattern, np.zeros([startPattern.shape[0], 2], dtype=np.int64)), axis=1)
    startPattern = np.concatenate((np.zeros([2, startPattern.shape[1]], dtype=np.int64), startPattern), axis=0)
    startPattern = np.concatenate((np.zeros([startPattern.shape[0], 2], dtype=np.int64), startPattern), axis=1)
    compareGrid = PatternGrid(copy.deepcopy(startPattern), [], -2, -2)
    periodGrid = copy.deepcopy(compareGrid)
    period = 1
    while True:
        if not periodGrid.advance():
            raise SyntaxError("RLE does not describe stable ash", (startAshInput, None, None, None))
        if periodGrid.strictlyEquals(compareGrid):
            break
        period += 1
        if period > testRunLimit:
            raise SyntaxError(f"pattern does not stabilize in {testRunLimit} generations", (startAshInput, None, None, None))
    startFlotilla = Flotilla(copy.deepcopy(startPattern), ShadowGrid(copy.deepcopy(startPattern), period), [], 0, period, None)

    if type(spaceshipOutput) is not str:
        raise TypeError("spaceshipOutput must be str")

    minDivGens = kwargs.get('minDivGens', 3)
    if type(minDivGens) is not int:
        raise TypeError("minDivGens must be int")


    # If the program is searching for a deleter reaction
    deleterSearch = False
    resultAshes = []
    resultAshPeriods = []
    resultAshInOut = kwargs.get('resultAshInOut', None)
    if resultAshInOut is not None:
        if not isinstance(resultAshInOut, collections.abc.Sequence):
            raise TypeError("resultAshInOut must be a sequence")
        if len(resultAshInOut) != 2:
            raise TypeError("resultAshInOut must contain two filenames")
        if type(resultAshInOut[0]) is not str:
            raise TypeError("input file (0) must be str")
        if type(resultAshInOut[1]) is not str:
            raise TypeError("output file (1) must be str")

        # Read resultAshInOut[0] as RLEs and trace ash periods
        start = re.sub(re.compile("^#.+", re.MULTILINE), "", re.sub("[^\S\n]", "", open(resultAshInOut[0], "r").read()))
        line = 1
        j = -1
        while True:
            match = re.search(re.compile("[^\n].*", re.DOTALL), start[j + 1:])
            if match is None:
                break
            start = match.group(0)
            line += match.start() - j - 1
            header = re.match("x=[0-9]+,y=[0-9]+,rule=23/3/3\n", start)
            if header is None:
                raise SyntaxError("header is not properly formatted", (resultAshInOut[0], line, None, None))
            width = int(re.search("(?<=x=)[0-9]+(?=,)", header.group(0)).group(0))
            if width == 0:
                raise SyntaxError("width cannot be 0", (resultAshInOut[0], line, None, header.group(0)))
            height = int(re.search("(?<=y=)[0-9]+(?=,)", header.group(0)).group(0))
            if height == 0:
                raise SyntaxError("height cannot be 0", (resultAshInOut[0], line, None, header.group(0)))
            newPattern = np.zeros([width, height], dtype=np.int64)
            x = 0
            y = 0
            i = header.end()
            j = start.find("!")
            if j == -1:
                j = len(start)
            line += 1
            while i < j:
                if start[i] == "\n":
                    line += 1
                    i += 1
                    continue
                run = 1
                if start[i].isdigit():
                    runStr = ""
                    k = i
                    while start[k].isdigit():
                        runStr += start[k]
                        k += 1
                        while k < len(start) and start[k] == "\n":
                            line += 1
                            k += 1
                        if k == len(start):
                            raise SyntaxError("file ends with dangling run", (resultAshInOut[0], line, None, None))
                    run = int(runStr)
                    i = k
                if start[i] == "$":
                    x = 0
                    y += run
                    if y > height:
                        raise SyntaxError("RLE body overflows provided height", (resultAshInOut[0], line, None, None))
                else:
                    if x + run > width:
                        raise SyntaxError("RLE body overflows provided width", (resultAshInOut[0], line, None, None))
                    if start[i] == ".":
                        x += run
                    elif start[i] == "A":
                        for k in range(run):
                            newPattern[x][y] = 1
                            x += 1
                    elif start[i] == "B":
                        for k in range(run):
                            newPattern[x][y] = 2
                            x += 1
                    else:
                        if re.match("x=[0-9]+,y=[0-9]+,rule=23/3/3\n", start[i:]) is not None:
                            start = start[i:]
                            j = -1
                            break
                        raise SyntaxError("file contains unrecognized symbol", (resultAshInOut[0], line, None, start[i]))
                i += 1
            xMin = -1
            scanning = True
            while scanning and (xMin := xMin + 1) < newPattern.shape[0]:
                for y in range(newPattern.shape[1]):
                    if newPattern[xMin][y] != 0:
                        scanning = False
                        break
            if xMin == newPattern.shape[0]:
                deleterSearch = True
                continue
            newPattern = newPattern[xMin:, :]
            yMin = -1
            scanning = True
            while scanning and (yMin := yMin + 1) < newPattern.shape[1]:
                for x in range(newPattern.shape[0]):
                    if newPattern[x][yMin] != 0:
                        scanning = False
                        break
            newPattern = newPattern[:, yMin:]
            xMax = newPattern.shape[0]
            scanning = True
            while scanning and (xMax := xMax - 1) >= 0:
                for y in range(newPattern.shape[1]):
                    if newPattern[xMax][y] != 0:
                        scanning = False
                        break
            newPattern = newPattern[:xMax + 1, :]
            yMax = newPattern.shape[1]
            scanning = True
            while scanning and (yMax := yMax - 1) >= 0:
                for x in range(newPattern.shape[0]):
                    if newPattern[x][yMax] != 0:
                        scanning = False
                        break
            newPattern = newPattern[:, :yMax + 1]
            newPattern = np.concatenate((newPattern, np.zeros([2, newPattern.shape[1]], dtype=np.int64)), axis=0)
            newPattern = np.concatenate((newPattern, np.zeros([newPattern.shape[0], 2], dtype=np.int64)), axis=1)
            newPattern = np.concatenate((np.zeros([2, newPattern.shape[1]], dtype=np.int64), newPattern), axis=0)
            newPattern = np.concatenate((np.zeros([newPattern.shape[0], 2], dtype=np.int64), newPattern), axis=1)
            resultAshes.append(newPattern)
            compareGrid = PatternGrid(copy.deepcopy(newPattern), [], -2, -2)
            periodGrid = copy.deepcopy(compareGrid)
            period = 1
            while True:
                if not periodGrid.advance():
                    raise SyntaxError("RLE does not describe stable ash", (resultAshInOut[0], line, None, None))
                if periodGrid.strictlyEquals(compareGrid):
                    break
                period += 1
                if period > testRunLimit:
                    raise SyntaxError(f"pattern does not stabilize in {testRunLimit} generations", (resultAshInOut[0], line, None, None))
            resultAshPeriods.append(period)
        if len(resultAshes) == 0 and not deleterSearch:
            raise SyntaxError("file is empty", (resultAshInOut[0], None, None, None))

    # Files for reactions that produce spaceships, and reactions that produce a stable ash
    # that is a desired result listed in the first file named in resultAshInOut, respectively
    sOFile = open(spaceshipOutput, "a")
    sOFile.truncate(0)
    sOFile.close()

    if resultAshInOut is not None:
        rAOFile = open(resultAshInOut[1], "a")
        rAOFile.truncate(0)
        rAOFile.close()

    # List of all stable ashes produced by flotillas
    ashes = [copy.deepcopy(startPattern)]
    # List of all flotillas in the previous search level where no parvoships are destroyed
    flotillas = [startFlotilla]
    # List of all flotillas in the previous search level where some parvoships are destroyed
    unstableFlotillas = []
    # Number of parvoships deep the program is into the search space
    parvoCount = 1

    # Counts for reactions added to sOFile and rAOFile
    sOCount = 0
    rAOCount = 0
    while True:
        sOFile = open(spaceshipOutput, "a")
        if resultAshInOut is not None:
            rAOFile = open(resultAshInOut[1], "a")

        newFlotillas = []
        newUnstableFlotillas = []
        banner = f"{'/' * 150}\n{parvoCount} parvo search space\n{'/' * 150}\n\n"
        sOFile.write(banner)
        if resultAshInOut is not None:
            rAOFile.write(banner)
        # Used to distinguish reactions in the out files
        sOID = 0
        rAOID = 0
        # Number of base flotillas searched in this search level
        flotillaCount = 0
        # Run through stable flotillas found in previous search level and attempt to add new parvoships
        for flotilla in flotillas:
            # SEARCH SIMULATION:
            searchGrid = PatternGrid(copy.deepcopy(flotilla.startAsh), copy.deepcopy(flotilla.parvos), -2, -2)
            shadows = flotilla.startShadow.parvoShadows(flotilla.parvos)
            # Advance searchGrid to first interaction between puff and the most recently added parvoship
            for gen in range(flotilla.startGen):
                puff = searchGrid.parvosRemoved()
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.cellshadow(xCell, yCell, gen)
                            for leftHanded in [True, False]:
                                for phase in range(3):
                                    phaseStart = (phase - gen) % 4
                                    markShadow = (1 if leftHanded else 16) << phaseStart
                                    yShiftBack = yCell + (gen + 1 - phase % 2) // 2
                                    for coords in parvoSearchCoords[phase]:
                                        xShadow = xCell - coords[0] if leftHanded else xCell + coords[0]
                                        yShadow = yShiftBack - coords[1]
                                        shadows[xShadow, yShadow] = shadows[xShadow, yShadow] | markShadow
                            shadows.postshadow(xCell, yCell, gen)
                searchGrid.advance()
            searchGen = flotilla.startGen
            # March through evolution of searchGrid
            while searchGen < flotilla.stopGen:
                puff = searchGrid.parvosRemoved()
                # Cast pre search shadow
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.cellshadow(xCell, yCell, searchGen)
                # Search for positions to insert parvoships and cast search shadow
                for phase in range(3):
                    for leftHanded in [True, False]:
                        # Cataloge all possible insertions of parvoships
                        phaseStart = (phase - searchGen) % 4
                        markShadow = (1 if leftHanded else 16) << phaseStart
                        testCoords = []
                        for yCell in range(puff.yStart, puff.yEnd()):
                            for xCell in range(puff.xStart, puff.xEnd()):
                                if puff[xCell, yCell] == 1:
                                    yShiftBack = yCell + (searchGen + 1 - phase % 2) // 2
                                    for coords in parvoSearchCoords[phase]:
                                        xShadow = xCell - coords[0] if leftHanded else xCell + coords[0]
                                        yShadow = yShiftBack - coords[1]
                                        if shadows[xShadow, yShadow] & markShadow == 0:
                                            testCoords.append([xShadow, yShadow])
                                            shadows[xShadow, yShadow] = shadows[xShadow, yShadow] | markShadow
                        # Test each insertion
                        for coords in testCoords:
                            # TEST SIMULATION:
                            # Test loop ignores all possible insertions that are above the most recently added parvoship,
                            # and prevents the same flotilla from being created in a different order
                            if len(flotilla.parvos) > 0 and \
                                    (coords[1] < flotilla.parvos[-1].y or
                                     coords[1] == flotilla.parvos[-1].y and coords[0] < flotilla.parvos[-1].x and searchGen == flotilla.startGen):
                                continue
                            testGrid = PatternGrid(copy.deepcopy(flotilla.startAsh), copy.deepcopy(flotilla.parvos) + [Parvo(phaseStart, coords[0], coords[1], leftHanded)], -2, -2)
                            testParvos = copy.deepcopy(testGrid.parvos)
                            # Track the bottom of the parvoship arrangement for four generations, to be referenced to track the bottom of the arrangement at any given generation
                            bottomsSequenceParvos = copy.deepcopy(testParvos)
                            i = len(bottomsSequenceParvos) - 2
                            endBottom = bottomsSequenceParvos[-1].bottom()
                            while i >= 0 and bottomsSequenceParvos[i].bottom() > endBottom - 2:
                                i -= 1
                            bottomsSequence = []
                            for slice in range(4):
                                bottom = endBottom
                                for iParvo in range(len(bottomsSequenceParvos) - 1, i, -1):
                                    tentativeBottom = bottomsSequenceParvos[iParvo].bottom()
                                    if tentativeBottom > bottom:
                                        bottom = tentativeBottom
                                bottomsSequence.append(bottom)
                                if slice != 3:
                                    for parvo in bottomsSequenceParvos:
                                        parvo.advance()
                            # Simulation of the flotilla being used as a base, to determine if testGrid diverges
                            immediateParent = PatternGrid(copy.deepcopy(flotilla.startAsh), copy.deepcopy(flotilla.parvos), -2, -2)
                            # If immediateParent still contains puff
                            parentAlive = True
                            # If testGrid diverges
                            diverges = False
                            # This flag tells the test loop to continue to the next possible insertion
                            verified = False
                            testGen = 0
                            # Tracks how many consecutive generations testGrid diverges from immediate parent
                            consecutiveDivergentGens = 0
                            # Advance testGrid to first interaction between puff and the most recently added parvoship
                            while testGen < searchGen:
                                # Track if testGrid diverges
                                if not diverges:
                                    testPuff = testGrid.parvosRemoved()
                                    parentPuff = immediateParent.parvosRemoved()
                                    if testPuff == parentPuff:
                                        consecutiveDivergentGens = 0
                                    else:
                                        consecutiveDivergentGens += 1
                                    if consecutiveDivergentGens == minDivGens:
                                        diverges = True
                                # Test for destroyed parvoships; if so continue to next possible insertion
                                for parvo in testGrid.parvos:
                                    head = parvo.leadingEdge()
                                    if testGrid.xStart < head[0] and head[0] < testGrid.xEnd() - 1 and testGrid.yStart < head[1] and head[1] < testGrid.yEnd() - 1 and testGrid[head] != 1:
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Advance testGrid; test if puff has died
                                if not testGrid.advance():
                                    # If searching for deleter reaction, add reaction to result ash output file
                                    if diverges and deleterSearch:
                                        rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, testGen + 1, flotilla.parent))}\n\n")
                                        rAOID += 1
                                        rAOCount += 1
                                    # Otherwise continue to next possible insertion
                                    verified = True
                                    break
                                # Advance immediateParent; if puff has died in immediateParent, cease to advance immediateParent and confirm that testGrid diverges
                                if parentAlive:
                                    parentAlive = immediateParent.advance()
                                    if not parentAlive:
                                        diverges = True
                                testGen += 1
                            if verified:
                                continue
                            # Record of testGrid puff after supposed first interaction with most recently added parvoship
                            puffRecords = []
                            mostRecentSeparation = -1
                            # Simulate testGrid; test for puff stabilization
                            while not verified:
                                # Track if testGrid diverges
                                if not diverges:
                                    testPuff = testGrid.parvosRemoved()
                                    parentPuff = immediateParent.parvosRemoved()
                                    if testPuff == parentPuff:
                                        consecutiveDivergentGens = 0
                                    else:
                                        consecutiveDivergentGens += 1
                                    if consecutiveDivergentGens == minDivGens:
                                        diverges = True
                                # Test for destroyed parvoships; if so add new flotilla to newUnstableFlotillas and continue to next possible insertion
                                for parvo in testGrid.parvos:
                                    head = parvo.leadingEdge()
                                    if testGrid.xStart < head[0] and head[0] < testGrid.xEnd() - 1 and testGrid.yStart < head[1] and head[1] < testGrid.yEnd() - 1 and testGrid[head] != 1:
                                        newUnstableFlotillas.append(Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, testGen, flotilla.parent))
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Test for puff stabilization
                                currentPuff = testGrid.parvosRemoved()
                                # iRecord scans puffRecords starting from the previous generation and going backwards
                                # to the last generation the bottom of the arrangement of the parvos did not intersect with testGrid
                                for iRecord in range(-1, mostRecentSeparation, -1):
                                    # Test if scan has found match to puff of current generation
                                    if currentPuff.strictlyEquals(puffRecords[iRecord]):
                                        # If testGrid diverges, add new flotilla and product ash to newFlotillas
                                        if diverges:
                                            # jRecord scans from the start of puffRecords onwards to find the earliest occurrence of the stable ash
                                            jRecord = 0
                                            while jRecord < len(puffRecords) and puffRecords[jRecord] != puffRecords[(jRecord - len(puffRecords) + 1) % iRecord - 1]:
                                                jRecord += 1
                                            # Add new flotilla to newFlotillas
                                            newFlotilla = Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, searchGen + jRecord, flotilla.parent)
                                            newFlotillas.append(newFlotilla)
                                            # If product ash matches result ash pattern: add reaction to result ash output file
                                            for iResult in range(len(resultAshes)):
                                                if resultAshPeriods[iResult] == -iRecord:
                                                    for kRecord in range(iRecord, 0):
                                                        if np.array_equal(puffRecords[kRecord].array, resultAshes[iResult]) or np.array_equal(puffRecords[kRecord].array, np.flipud(resultAshes[iResult])):
                                                            rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(newFlotilla)}\n\n")
                                                            rAOID += 1
                                                            rAOCount += 1
                                            # If product ash has already been found, do not add product ash to newFlotillas
                                            for ash in ashes:
                                                for kRecord in range(iRecord, 0):
                                                    if np.array_equal(ash, puffRecords[kRecord].array) or np.array_equal(np.flipud(ash), puffRecords[kRecord].array):
                                                        verified = True
                                                        break
                                                if verified:
                                                    break
                                            if verified:
                                                break
                                            # Add product ash to ashes, new flotilla to newFlotillas
                                            ashes.append(currentPuff.array)
                                            newFlotillas.append(Flotilla(currentPuff.array, ShadowGrid(copy.deepcopy(currentPuff.array), -iRecord), [], 0, -iRecord, newFlotilla))
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Test if bottom of parvoship arrangement has reintersected with testGrid
                                if bottomsSequence[testGen % 4] - testGen // 4 * 2 <= testGrid.yStart:
                                    mostRecentSeparation -= 1
                                else:
                                    mostRecentSeparation = -1
                                # If testGen has run for too long, assume flotilla has produced spaceship; add reaction to spaceship producing reactions file
                                if testGen > searchGen + testRunLimit:
                                    sOFile.write(f"Spaceship producing reaction #{parvoCount}-{sOID}:\n{reactionHistory(Flotilla(flotilla.startAsh, None, copy.deepcopy(testParvos), 0, 1, flotilla.parent))}\n\n")
                                    sOID += 1
                                    sOCount += 1
                                    break
                                puffRecords.append(currentPuff)
                                # If puff in testGrid dies and testGrid diverges, add new flotilla to newFlotillas
                                if not testGrid.advance():
                                    if diverges:
                                        newFlotilla = Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, testGen + 1, flotilla.parent)
                                        # If searching for deleter reaction, add reaction to result ash output file
                                        if deleterSearch:
                                            rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(newFlotilla)}\n\n")
                                            rAOID += 1
                                            rAOCount += 1
                                        # Do not add new flotilla to newFlotillas if searching for deleter reaction: will result in pileup of reactions in result ash output file
                                        else:
                                            newFlotillas.append(newFlotilla)
                                    break
                                # Advance immediateParent; if puff has died in immediateParent, cease to advance immediateParent and confirm that testGrid diverges
                                if parentAlive:
                                    parentAlive = immediateParent.advance()
                                    if not parentAlive:
                                        diverges = True
                                testGen += 1

                # Cast post search shadow
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.postshadow(xCell, yCell, searchGen)
                searchGen += 1
                searchGrid.advance()
            flotillaCount += 1
            print(f"{flotillaCount}/{len(flotillas)} stable base flotillas searched")

        # Number of base flotillas searched in this search level
        flotillaCount = 0
        # Run through unstable flotillas found in previous search level and attempt to add new parvoships to create a stable flotilla
        # Divergence is not tested here; if no parvoships are destroyed, by definition the new flotilla has diverged
        for flotilla in unstableFlotillas:
            # SEARCH SIMULATION:
            searchGrid = PatternGrid(copy.deepcopy(flotilla.startAsh), copy.deepcopy(flotilla.parvos), -2, -2)
            shadows = flotilla.startShadow.parvoShadows(flotilla.parvos)
            # Advance searchGrid to first interaction between puff and the most recently added parvoship
            for gen in range(flotilla.startGen):
                puff = searchGrid.parvosRemoved()
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.cellshadow(xCell, yCell, gen)
                            for leftHanded in [True, False]:
                                for phase in range(3):
                                    phaseStart = (phase - gen) % 4
                                    markShadow = (1 if leftHanded else 16) << phaseStart
                                    yShiftBack = yCell + (gen + 1 - phase % 2) // 2
                                    for coords in parvoSearchCoords[phase]:
                                        xShadow = xCell - coords[0] if leftHanded else xCell + coords[0]
                                        yShadow = yShiftBack - coords[1]
                                        shadows[xShadow, yShadow] = shadows[xShadow, yShadow] | markShadow
                            shadows.postshadow(xCell, yCell, gen)
                searchGrid.advance()
            searchGen = flotilla.startGen
            # March through evolution of searchGrid
            while searchGen < flotilla.stopGen:
                puff = searchGrid.parvosRemoved()
                # Cast pre search shadow
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.cellshadow(xCell, yCell, searchGen)
                # Search for positions to insert parvoships and cast search shadow
                for phase in range(3):
                    for leftHanded in [True, False]:
                        # Cataloge all possible insertions of parvoships
                        phaseStart = (phase - searchGen) % 4
                        markShadow = (1 if leftHanded else 16) << phaseStart
                        testCoords = []
                        for yCell in range(puff.yStart, puff.yEnd()):
                            for xCell in range(puff.xStart, puff.xEnd()):
                                if puff[xCell, yCell] == 1:
                                    yShiftBack = yCell + (searchGen + 1 - phase % 2) // 2
                                    for coords in parvoSearchCoords[phase]:
                                        xShadow = xCell - coords[0] if leftHanded else xCell + coords[0]
                                        yShadow = yShiftBack - coords[1]
                                        if shadows[xShadow, yShadow] & markShadow == 0:
                                            testCoords.append([xShadow, yShadow])
                                            shadows[xShadow, yShadow] = shadows[xShadow, yShadow] | markShadow
                        # Test each insertion
                        for coords in testCoords:
                            # TEST SIMULATION:
                            # Test loop ignores all possible insertions that are above the most recently added parvoship,
                            # and prevents the same flotilla from being created in a different order
                            if len(flotilla.parvos) > 0 and \
                                    (coords[1] < flotilla.parvos[-1].y or
                                     coords[1] == flotilla.parvos[-1].y and coords[0] < flotilla.parvos[-1].x and searchGen == flotilla.startGen):
                                continue
                            testGrid = PatternGrid(copy.deepcopy(flotilla.startAsh), copy.deepcopy(flotilla.parvos) + [Parvo(phaseStart, coords[0], coords[1], leftHanded)], -2, -2)
                            testParvos = copy.deepcopy(testGrid.parvos)
                            # Track the bottom of the parvoship arrangement for four generations, to be referenced to track the bottom of the arrangement at any given generation
                            bottomsSequenceParvos = copy.deepcopy(testParvos)
                            i = len(bottomsSequenceParvos) - 2
                            endBottom = bottomsSequenceParvos[-1].bottom()
                            while i >= 0 and bottomsSequenceParvos[i].bottom() > endBottom - 2:
                                i -= 1
                            bottomsSequence = []
                            for slice in range(4):
                                bottom = endBottom
                                for iParvo in range(len(bottomsSequenceParvos) - 1, i, -1):
                                    tentativeBottom = bottomsSequenceParvos[iParvo].bottom()
                                    if tentativeBottom < bottom:
                                        bottom = tentativeBottom
                                bottomsSequence.append(bottom)
                                if slice != 3:
                                    for parvo in bottomsSequenceParvos:
                                        parvo.advance()
                            # This flag tells the test loop to continue to the next possible insertion
                            verified = False
                            testGen = 0
                            # Advance testGrid to first interaction between puff and the most recently added parvoship
                            while testGen < searchGen:
                                for parvo in testGrid.parvos:
                                    head = parvo.leadingEdge()
                                    if testGrid.xStart < head[0] and head[0] < testGrid.xEnd() - 1 and testGrid.yStart < head[1] and head[1] < testGrid.yEnd() - 1 and testGrid[head] != 1:
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Advance testGrid; test if puff has died
                                if not testGrid.advance():
                                    # If searching for deleter reaction, add reaction to result ash output file
                                    if deleterSearch:
                                        rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, testGen + 1, flotilla.parent))}\n\n")
                                        rAOID += 1
                                        rAOCount += 1
                                    # Otherwise continue to next possible insertion
                                    verified = True
                                    break
                                testGen += 1
                            if verified:
                                continue
                            # Record of testGrid puff after supposed first interaction with most recently added parvoship
                            puffRecords = []
                            mostRecentSeparation = -1
                            # Simulate testGrid; test for puff stabilization
                            while not verified:
                                # Test for destroyed parvoships; if so continue to next possible insertion
                                for parvo in testGrid.parvos:
                                    head = parvo.leadingEdge()
                                    if testGrid.xStart < head[0] and head[0] < testGrid.xEnd() - 1 and testGrid.yStart < head[1] and head[1] < testGrid.yEnd() - 1 and testGrid[head] != 1:
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Test for puff stabilization
                                currentPuff = testGrid.parvosRemoved()
                                # iRecord scans puffRecords starting from the previous generation and going backwards
                                # to the last generation the bottom of the arrangement of the parvos did not intersect with testGrid
                                for iRecord in range(-1, mostRecentSeparation, -1):
                                    # Test if scan has found match to puff of current generation
                                    if currentPuff.strictlyEquals(puffRecords[iRecord]):
                                        # jRecord scans from the start of puffRecords onwards to find the earliest occurrence of the stable ash
                                        jRecord = 0
                                        while jRecord < len(puffRecords) and puffRecords[jRecord] != puffRecords[(jRecord - len(puffRecords) + 1) % iRecord - 1]:
                                            jRecord += 1
                                        # Add new flotilla to newFlotillas
                                        newFlotilla = Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, searchGen + jRecord, flotilla.parent)
                                        newFlotillas.append(newFlotilla)
                                        # If product ash matches result ash pattern: add reaction to result ash output file
                                        for iResult in range(len(resultAshes)):
                                            if resultAshPeriods[iResult] == -iRecord:
                                                for kRecord in range(iRecord, 0):
                                                    if np.array_equal(puffRecords[kRecord].array, resultAshes[iResult]) or np.array_equal(puffRecords[kRecord].array, np.flipud(resultAshes[iResult])):
                                                        rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(newFlotilla)}\n\n")
                                                        rAOID += 1
                                                        rAOCount += 1
                                        # If product ash has already been found, do not add product ash to newFlotillas
                                        for ash in ashes:
                                            for kRecord in range(iRecord, 0):
                                                if np.array_equal(ash, puffRecords[kRecord].array) or np.array_equal(np.flipud(ash), puffRecords[kRecord].array):
                                                    verified = True
                                                    break
                                            if verified:
                                                break
                                        if verified:
                                            break
                                        # Add product ash to ashes, new flotilla to newFlotillas
                                        ashes.append(currentPuff.array)
                                        newFlotillas.append(Flotilla(currentPuff.array, ShadowGrid(copy.deepcopy(currentPuff.array), -iRecord), [], 0, -iRecord, newFlotilla))
                                        verified = True
                                        break
                                if verified:
                                    break
                                # Test if bottom of parvoship arrangement has reintersected with testGrid
                                if bottomsSequence[testGen % 4] - testGen // 4 * 2 <= testGrid.yStart:
                                    mostRecentSeparation -= 1
                                else:
                                    mostRecentSeparation = -1
                                # If testGen has run for too long, assume flotilla has produced spaceship; add reaction to spaceship producing reactions file
                                if testGen > searchGen + testRunLimit:
                                    sOFile.write(f"Spaceship producing reaction #{parvoCount}-{sOID}:\n{reactionHistory(Flotilla(flotilla.startAsh, None, copy.deepcopy(testParvos), 0, 1, flotilla.parent))}\n\n")
                                    sOID += 1
                                    sOCount += 1
                                    break
                                puffRecords.append(currentPuff)
                                # If puff in testGrid dies, add new flotilla to newFlotillas
                                if not testGrid.advance():
                                    newFlotilla = Flotilla(flotilla.startAsh, flotilla.startShadow, copy.deepcopy(testParvos), searchGen, testGen + 1, flotilla.parent)
                                    # If searching for deleter reaction, add reaction to result ash output file
                                    if deleterSearch:
                                        rAOFile.write(f"Desired ash producing reaction #{parvoCount}-{rAOID}:\n{reactionHistory(newFlotilla)}\n\n")
                                        rAOID += 1
                                        rAOCount += 1
                                    # Do not add new flotilla to newFlotillas if searching for deleter reaction: will result in pileup of reactions in result ash output file
                                    else:
                                        newFlotillas.append(newFlotilla)
                                    break
                                testGen += 1

                # Cast post search shadow
                for yCell in range(puff.yStart, puff.yEnd()):
                    for xCell in range(puff.xStart, puff.xEnd()):
                        if puff[xCell, yCell] == 1:
                            shadows.postshadow(xCell, yCell, searchGen)
                searchGen += 1
                searchGrid.advance()
            flotillaCount += 1
            print(f"{flotillaCount}/{len(unstableFlotillas)} unstable base flotillas searched")
        flotillas = newFlotillas
        unstableFlotillas = newUnstableFlotillas
        # Query user for next steps regarding next level of search space
        sOFile.close()
        if resultAshInOut is not None:
            rAOFile.close()
        if len(flotillas)==0 and len(unstableFlotillas)==0:
            print(f"{parvoCount} parvo search space explored;\nAll possible branches explored; Search ended.\n{sOCount} total spaceship producing reactions found.\n{rAOCount} total desired ash producing reactions found.")
            break
        inp = ""
        while inp.lower() not in ["continue", "exit", "quit"]:
            inp = input(f"{parvoCount} parvo search space explored;\n{len(flotillas)} flotillas found in this space.\n{sOCount} total spaceship producing reactions found.\n{rAOCount} total desired ash producing reactions found.\nContinue to {parvoCount + 1} parvo search space? (continue/exit/montage): ").lower().strip()
            if inp == "montage":
                if len(flotillas) == 0:
                    print("No stable flotillas found.")
                    continue
                size = 256
                display = pygame.display.set_mode([size * 3, size * 3])
                disp = pygame.Surface((size, size))

                started = False
                iFlotillas = 0
                gen = 0
                grid = PatternGrid(copy.deepcopy(flotillas[iFlotillas].startAsh), copy.deepcopy(flotillas[iFlotillas].parvos), -2, -2)

                def flip():
                    disp.fill((0, 0, 0))
                    if grid.array is not None:
                        for y in range(grid.yStart, grid.yEnd()):
                            for x in range(grid.xStart, grid.xEnd()):
                                state = grid[x, y]
                                if state == 1:
                                    disp.set_at((x + 128, y + 128), (255, 0, 0))
                                if state == 2:
                                    disp.set_at((x + 128, y + 128), (255, 255, 0))
                    for parvo in grid.parvos:
                        if grid.array is None or parvo.left() < grid.xStart or grid.xEnd() < parvo.right() or parvo.top() < grid.yStart or grid.yEnd() < parvo.bottom():
                            for state in range(2):
                                for coords in parvoBodyCoords[parvo.phase][state]:
                                    x = parvo.x + coords[0] if (parvo.leftHanded) else parvo.x - coords[0]
                                    y = parvo.y + coords[1]
                                    if grid.array is None or x < grid.xStart or grid.xEnd() <= x or y < grid.yStart or grid.yEnd() <= y:
                                        disp.set_at((x + 128, y + 128), ((255, 0, 0), (255, 255, 0))[state])
                    display.blit(pygame.transform.scale(disp, (size * 3, size * 3)), (0, 0))

                clock = pygame.time.Clock()
                main = True
                flip()
                while main:
                    pygame.display.flip()
                    clock.tick(60)
                    if started and gen < flotillas[iFlotillas].stopGen:
                        grid.advance()
                        flip()
                        gen += 1
                    for event in pygame.event.get():
                        if event.type == pygame.QUIT:
                            pygame.quit()
                            main = False
                        if event.type == pygame.KEYDOWN:
                            if event.key == pygame.K_SPACE or event.key == pygame.K_RETURN:
                                if not started:
                                    started = True
                                elif gen == flotillas[iFlotillas].stopGen:
                                    gen = 0
                                    iFlotillas += 1
                                    if iFlotillas == len(flotillas):
                                        pygame.quit()
                                        main = False
                                        break
                                    else:
                                        grid = PatternGrid(copy.deepcopy(flotillas[iFlotillas].startAsh), copy.deepcopy(flotillas[iFlotillas].parvos), -2, -2)
                                        started = False
                                        flip()
                            if event.key == pygame.K_p:
                                print("Reaction RLE:")
                                print(RLE(PatternGrid(copy.deepcopy(flotillas[iFlotillas].startAsh), copy.deepcopy(flotillas[iFlotillas].parvos), -2, -2)))
        if inp in ["exit", "quit"]:
            break
        parvoCount += 1
#from HushNowGregory import search
if __name__ == '__main__':
    search("HNGFiles/startAshInputRLE.txt","HNGFiles/spaceshipOutputRLEs.txt",resultAshInOut=("HNGFiles/resultAshInRLEs.txt","HNGFiles/resultAshOutRLEs.txt"))
