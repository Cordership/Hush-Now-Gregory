# Hush-Now-Gregory
The program is set to run with default settings if run as \_\_main\_\_.
Doing so requires the folder HNGFiles to be in the path.
Otherwise, import the python file in a different script and run HushNowGregory.search(startAshInput, spaceshipOutput, \*\*kwargs).
Arguments go as follows:

Arguments:
    startAshInput:
        File containing RLE describing initial ash that the program will attempt to perturb to produce spaceships/specified stable output
    spaceshipOutput:
        File to output reactions that produce spaceships onto
Keyword Arguments:
    minDivGens:
        Minimum number of consecutive generations a reaction must deviate based on its most recently added parvoship
        for the addition of said parvoship to be considered a meaningful search branch
        Defaults to 3
    testRunLimit:
        Number of generations after the projected first interaction between the puff and the most recently added parvoship that the program
        will consider the reaction a spaceship producing reaction
        Defaults to 300
    resultAshInOut:
        Tuple like object containing two filenames
        The program will search for reactions that produce stable ash that matches the RLEs listed in first file
        These reactions will be outputted in the second file
        If one of the RLEs in the first file contains only dead cells, the program will search for a deleter flotilla
        for the given start ash
Example:
    search("startAshInputRLE.txt","spaceshipOutputRLEs.txt",resultAshInOut=("resultAshInRLEs.txt","resultAshOutRLEs.txt"))
