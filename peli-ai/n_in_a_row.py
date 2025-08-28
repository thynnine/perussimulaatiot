import numpy as np
import os
import copy
from numpy.random import default_rng

# A game of 4 (or n) in a row.
# In this game, you drop markers on a grid in turns 
# and try to get 4 (or n) of your markers in a row.


# Game parameters
nx = 7 # number of columns
ny = 5 # number of rows
nrow = 4 # number of consecutive tokens needed for a victory
grid = np.zeros([ny, nx], dtype=int) # the play area

n_players = 2 # number of players (typically 2, but why not try 3?)

# Symbols for different players
P1 = "X"
P2 = "O"
P3 = "#"
P4 = "@"
PSYM = [" ", P1, P2, P3, P4]



def new_game():
    """
    Creates an empty play area.
    
    returns: the empty play area (the initial game state)
    """
    return np.zeros([ny, nx], dtype=int)
    

def can_continue(grid):
    """
    Checks if there are available spaces in the play area.
    
    Does not check if the game is done due to someone winning.
    
    returns: True if there are empty spaces, False if not.
    """
    for i in range(nx):
        if grid[0, i] == 0:
            return True

    return False


def all_plays(grid, player=None):
    """
    Lists the indices of all columns that have empty slots.
    
    returns: list of possible play actions
    """

    empty_slots = []

    for i in range(nx):
        if grid[0, i] == 0:
            empty_slots.append( i )
                
    return empty_slots


def all_reasonable_plays(grid, player=None):
    """
    Same as all_plays().
    
    returns: list of possible play actions
    """
    
    return all_plays(grid)
    
    
def check_winner(grid):
    """
    Checks if someone has won the game.
    
    returns: the number of the winning player if there is one, 0 otherwise
    """
    for i in range(ny):
        for j in range(nx):
            row = count_consecutives(i,j, grid)
            
            if row >= nrow:
                return grid[i,j]
    return 0
   

def previous_player(player):
    """
    Number of the player whose turn comes before the given player.
    
    Note: Player numbering starts from 1 and ends at n_players.
    
    player: the player to check
    
    returns: the number of the previous player
    """
    return (player)%n_players + 1
    
def next_player(player):
    """
    Number of the player whose turn comes after the given player.
    
    Note: Player numbering starts from 1 and ends at n_players.
    
    player: the player to check
    
    returns: the number of the next player
    """
    return (player)%n_players + 1
    
    
def count_consecutives(i,j,grid):
    """
    Starting from the slot at grid[i,j], counts the highest number
    of consecutive squares with tokens of the same player.
    
    If the starting position has no token, returns 0.
    
    For ay given square, there are 8 possible directions for forming
    a line: up, down, left, rigth and the 4 diagonals. This function only checks
    half of them: down, right, down-right and down-left.
    
    i: vertical position (y coordinate) of the initial slot
    j: horizontal position (x coordinate) of the initial slot
    grid: the play area (current game state)
    
    returns: the highest number of similar consecutive tokens
    """
    
    # If there is no token at this position, return 0
    type = grid[i,j]
    if type == 0:
        return 0
    
    streak = True
    count = 1
    
    # count horizontally
    ii = i+1
    while streak:
        if ii < ny:
            if grid[ii,j] == type:
                count += 1
                ii += 1
            else:
                streak = False
        else:
            streak = False
            
    max_count = count   
    count = 1 
    streak = True
    
    # count vertically
    jj = j+1
    while streak:
        if jj < nx:
            if grid[i,jj] == type:
                count += 1
                jj += 1
            else:
                streak = False
        else:
            streak = False
         
    if count > max_count:
        max_count = count   
    count = 1 
    streak = True
    
    # count diagonally down-left
    ii = i+1
    jj = j+1
    while streak:
        if jj < nx and ii < ny:
            if grid[ii,jj] == type:
                count += 1
                ii += 1
                jj += 1
            else:
                streak = False
        else:
            streak = False
    
    if count > max_count:
        max_count = count   
    count = 1 
    streak = True
    
    # count diagonally down-right
    ii = i+1
    jj = j-1
    while streak:
        if jj >= 0 and ii < ny:
            if grid[ii,jj] == type:
                count += 1
                ii += 1
                jj -= 1
            else:
                streak = False
        else:
            streak = False
                
    if count > max_count:
        max_count = count   
        
    return max_count
    
    

def draw(grid):
    """
    Draw the play area with text graphics.
    
    grid: the play area (current game state)
    """

    # Clear the terminal.
    # The command depends on the operating system.
    os.system('cls' if os.name == 'nt' else 'clear')
    
    # Construct the header with column indices and
    # the horizontal lines.
    topheader = "    "
    vline = "    "
    for j in range(nx):
        if j < 10:
            topheader += " "+str(j)+"  "
        elif j < 100:
            topheader += " "+str(j)+" "
        vline += "---+"

    # print the header
    print(topheader)

    # loop over all horizontal rows
    for i in range(ny):

        # start each line with its index
        line = "    "
        if i < 10:
            line = "  "+str(i)+" "
        elif i < 100:
            line = " "+str(i)+" "
    
        # loop over all slots (columns) in this row
        for j in range(nx):
        
            # mark each slot as either empty or with a player token
            if grid[i,j] == 0: # empty slot
                line += "   |"
            else: # someone has put their token in this slot
                for type in range(1,n_players+1):
                    if grid[i,j] == type:
                        line += " "+PSYM[type]+" |"
        
        # print the row
        print(line[:-1])
        
        # after each row except the last one, print a separating line
        if i < ny-1:
            print(vline[:-1])
          

def make_move(play, player, grid):
    """
    Let a player have a turn.
    
    The function modifies the given game state (grid) to
    match the situation after the move has been carried out.
    
    play: the play action (move) to take
    player: number of the player who makes the move
    grid: the play area (current game state)
    """

    for i in range(ny):
        if grid[ny-i-1, play] == 0:
            grid[ny-i-1, play] = player
            return
                      

def can_take_slot(i,grid):
    """
    Checks if one can play in the given column.
    
    i: the column to check
    grid: the play area (current game state)
    
    returns: True if there are empty slots in column i, False otherwise
    """
    if grid[0, i] == 0:
        return True
    else:
        return False  

                
def ask_for_move(player,grid):
    """
    Asks a human player for a play action (move).
    
    This is done by visualizing the game and then
    asking for the index of the column where the player
    wants to play.
    
    player: number of the player whose turn it is
    grid: the play area (current game state)
    
    returns: the index of the column where the player decides to play
    """
    
    # visualize the situation
    draw(grid)   
     
    # notify whose turn it is
    print("You are "+PSYM[player])
    col = -1
    ok = False

    # keep asking for a decision until a valid decision is made    
    while not ok:
    
        # keep asking until the choice is a valid column index
        while col < 0 or col >= nx:
        
            try: # try to read input as an integer
                col = int(input("choose column: "))
                
            except: # if the input is invalid
                print("that's not an integer (enter -9 to quit)")
                col = -1
                
            if col == -9: # entering -9 quits the game
                quit()
        
        # check that the chosen column is not full
        if can_take_slot(col, grid):
            ok = True
        else:
            print("can't play there")
            col = -1
        
    return col


def copy_game(grid):
    """
    Returns a copy of the given game state.
    
    grid: the play area (current game state)
    
    returns: a copy of grid
    """
    return copy.deepcopy(grid)


def declare_winner(winner):
    """
    Celebrates the victory of a player or declares a draw.
    
    winner: number of the winning player or 0 for a draw
    """
    for i in range(1,len(PSYM)):
        if winner == i:
            print("")
            print(PSYM[i]+" won!")
            print("")

    if winner == 0:
        print("")
        print("The game ended in a draw.")
        print("")
        