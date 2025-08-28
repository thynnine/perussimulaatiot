import numpy as np
import os
import copy
from numpy.random import default_rng


# A game of tic-tac-toe.

# Game parameters
n = 4 # grid size
nrow = 3 # number of consecutive marks needed for a victory
grid = np.zeros([n,n], dtype=int) # the play area

n_players = 2 # number of players (will be set in game_ai)

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
    return np.zeros([n,n], dtype=int)
    

def can_continue(grid):
    """
    Checks if there are available spaces in the play area.
    
    Does not check if the game is done due to someone winning.
    
    returns: True if there are empty spaces, False if not.
    """
    for i in range(n):
        for j in range(n):
            if grid[i,j] == 0:
                return True

    return False


def all_plays(grid, player=None):
    """
    Lists the coordinates of all empty slots.
    
    returns: list of possible play actions
    """

    empty_slots = []

    for i in range(n):
        for j in range(n):
            if grid[i,j] == 0:
                empty_slots.append( [i,j] )
                
    return empty_slots


def all_reasonable_plays(grid, player=None):
    """
    Same as all_plays().
    
    returns: list of possible play actions
    """
    
    return all_plays(grid)
    

def all_reasonable_plays_alternative(grid, player=None):
    """
    Lists the coordinates of all empty slots that have an occupied neighbor.
    
    returns: list of possible play actions
    """
    
    good_slots = []

    for i in range(n):
        for j in range(n):
            if grid[i,j] == 0:
            
                good = False
                for di in [-1, 0, 1]:
                    ni = (i+di)%n

                    for dj in [-1, 0, 1]:
                        nj = (j+dj)%n
                        
                        if ni != 0 or nj != 0:
                            if grid[ni,nj] > 0:
                                good = True
                if good:
                    good_slots.append( [i,j] )
               
    if len(good_slots) == 0:
        return all_plays(grid)
    else:
        return good_slots


def check_winner(grid):
    """
    Checks if someone has won the game.
    
    returns: the number of the winning player if there is one, 0 otherwise
    """
    for i in range(n):
        for j in range(n):
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
    
    type = grid[i,j]
    if type == 0:
        return 0
    
    streak = True
    count = 1
    
    # count horizontally
    ii = i+1
    while streak:
        if ii < n:
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
        if jj < n:
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
        if jj < n and ii < n:
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
        if jj >= 0 and ii < n:
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
    for j in range(n):
        if j < 10:
            topheader += " "+str(j)+"  "
        elif j < 100:
            topheader += " "+str(j)+" "
        vline += "---+"

    # print the header
    print(topheader)

    # loop over all horizontal rows
    for i in range(n):

        # start each line with its index
        line = "    "
        if i < 10:
            line = "  "+str(i)+" "
        elif i < 100:
            line = " "+str(i)+" "
    
        # loop over all slots (columns) in this row
        for j in range(n):
        
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
        if i < n-1:
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
    
    row = play[0]
    col = play[1]
    grid[row, col] = player
              

def can_take_slot(i,j,grid):
    """
    Checks if one can play in the given slot.
    
    i: the row to check
    j: the column to check
    grid: the play area (current game state)
    
    returns: True if there are empty slots in column i, False otherwise
    """
    if grid[i,j] == 0:
        return True
    else:
        return False  

                
def ask_for_move(player,grid):
    """
    Asks a human player for a play action (move).
    
    This is done by visualizing the game and then
    asking for the index of the column and row where the player
    wants to play.
    
    player: number of the player whose turn it is
    grid: the play area (current game state)
    
    returns: the coordinates of the chosen slot as tuple (row, column)
    """
    draw(grid)   
     
    print("You are "+PSYM[player])
    row = -1
    col = -1
    ok = False
    
    while not ok:
        while col < 0 or col >= n:
            try:
                col = int(input("choose column: "))
            except:
                print("that's not an integer (enter -9 to quit)")
                col = -1
            if col == -9:
                quit()
        
        while row < 0 or row >= n:
            try:
                row = int(input("choose row: "))
            except:
                print("that's not an integer (enter -9 to quit)")
                col = -1  
            if row == -9:
                quit()              
        
        if can_take_slot(row, col, grid):
            ok = True
        else:
            print("can't play there")
            col = -1
            row = -1
        
    return row, col


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
            print(PSYM[i]+" won!")

    if winner == 0:
        print("The game ended in a draw.")
        
        