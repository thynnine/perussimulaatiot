import numpy as np
import os
from numpy.random import default_rng
#import tictac as game
import n_in_a_row as game


"""
A general AI for simple turn based games.

The game to be played must be imported as a module named game.
The module must implement these functions:

* new_game()
* can_continue(state)
* check_winner(state)
* previous_player(player)
* next_player(player)
* get_all_reasonable_plays(state)
* draw(state)
* make_move(play, player, state)
* copy_game(state)
* declare_winner(winner)

and define these variables
* n_players

"""

# Tags for marking human and cpu players
HUMAN = "bio"
CPU = "cpu"


class MC_node:
    """
    A node in the Monte Carlo tree.
    
    Each node records the win and play statistics and links to 
    the nodes that come before and after it in the tree.
    
    Args:
        player (int): player number for the player whose move lead to this node
        state: current state of the game
        layer (int): counts how deep in the tree this node is
        play: the play action ("move") that leads to this node
        explore (float): for large values the algorithm is more eager to try non-optimal plays
    """

    # Creates a node object
    def __init__(self, player, state, layer, play, explore=1.0):
        self.player = player
        self.wins = 0
        self.plays = 0
        self.layer = layer
        self.play = play
        self.state = state
        self.branches = []
        self.parent = None
        self.winner = game.check_winner(state)
        self.explore_factor = explore
        
    # String representation for the node object
    def __str__(self):
        info = "MC node:\n"
        info += " player: "+str(self.player)+"\n"        
        info += " layer: "+str(self.layer)+"\n"
        info += " play: "+str(self.play)+"\n"
        info += " state: \n"+str(self.state)+"\n"
        info += " wins / plays: "+str(self.wins)+" / "+str(self.plays)+"\n"
        info += " branches: "+str(len(self.branches))
        for b in self.branches:
            info += " "+str(b.play)
            
        return info+"\n"

                
    def has_branches(self):
        """
        Checks if more gameplay branches follow from this node.
        
        Returns: 
            bool: True if more nodes follow this node in the MC tree.
        """
        return len(self.branches) > 0
        
        
    def is_leaf(self):
        """
        Checks if this is a leaf node.
        
        A node is said to be a leaf node if
        
        (a) the game ends at this node (there are no following nodes) **or**
        (b) at least one of the following nodes has 0 recorded plays
        
        I.e., a leaf is a node at the outer edge of the tree.
        
        If this node is not terminal but the following nodes have not yet been
        created, the function calls :meth:`MC_node.create_branches` to create nodes
        for all possible plays.
        
        Returns: 
            bool: True if this is a leaf node.
        """
        
        # If there are no branches, this node may end the game.
        # Or maybe the branches just have not been created yet.
        if not self.has_branches():

            # If winner is 0, the game may have ended in a draw
            # or the following branches have not yet been created.
            if self.winner == 0: 
            
                self.create_branches()

                # If there are no following plays, so the game has ended in a draw.
                if len(self.branches) == 0:
                    return True

            else: # if winner is not 0, the game has been decided with someone winning.
                return True
    
        # If there are branches, the game is still ongoing.
        for b in self.branches:
            if b.plays == 0:
                return True

        return False

        
    def pick_unplayed_branch(self):
        """
        Randomly chooses one node out of the nodes that
        
        (a) follow this node in the MC tree **and**
        (b) have no recorded plays
        
        Returns: 
            MC_node: the chosen node
        """
        
        # first list all the unplayed branches
        unplayed = []
        for b in range(len(self.branches)):
            if self.branches[b].plays == 0:
                unplayed.append(b)
                
        # then randomly pick one
        chosen = random.choice(unplayed)
        
        return self.branches[chosen]
        
        
    def pick_played_branch(self):
        """
        Chooses the node that
        
        (a) follows this node in the MC tree **and**
        (b) has the highest play factor as given by :meth:`MC_node.bias_factor`.
        
        Returns: 
            MC_node: the chosen node
        """
                
        # calculate the bias factors for all branches
        # always save the best factor and branch
        best_branch = 0
        best_bias = self.branches[0].bias_factor()
    
        for b in range(1,len(self.branches)):
        
            bias = self.branches[b].bias_factor()
            
            # check if this is better than the previous best
            if bias > best_bias:
                best_bias = bias
                best_branch = b
                              
        return self.branches[best_branch]
        
        
    def pick_best_play(self):
        """
        Chooses the play action (i.e., move) recorded in the node that
        
        (a) follows this node **and**
        (b) has the most recorded plays
        
        Returns: 
            object: the most tried out play
        """
    
        # calculate the number of playouts for all branches
        # always save the highest play count and branch
        best_branch = 0
        most_plays = self.branches[0].plays
        
        for  b in range(1,len(self.branches)):
            plays = self.branches[b].plays
            if plays > most_plays:
                most_plays = plays
                best_branch = b
                
        # return the play action from the most tried branch
        return self.branches[best_branch].play
        
        
    def bias_factor(self):
        """
        Calculates the factor for choosing the next node for a node
        where all the following nodes have at least one recorded play.
        The higher this factor is, the more likely it is that this node is chosen.
        
        The factor is calculated as
        
        .. math ::
        
            \\frac{w}{n} + c \\sqrt{ \\frac{1}{n} \\ln(n_\\text{parent}) }
        
        where
        
            * :math:`w` = the recorded wins in all the playouts from this node
            * :math:`n` = the number of recorded playouts from this node
            * :math:`c` = the exploration factor
            * :math:`n_\\text{parent}` = the number of recorded playouts from the previous node
        
        The first term :math:`w / n` is the win ratio for the playouts from this node.
        If this node is likely to lead to victory, this ratio is close to 1
        and gets chosen often.
        
        The second term is large, if there are only a few plays from this node and
        many plays from the *other nodes* that branch off the same parent node (the competing
        play options). If this factor was not included, any node whose initial play lead
        to a loss would never be considered again by the algorithm. Since this may
        happen simply due to bad luck, it is important to give these nodes some additional
        plays.
        
        Returns: 
            float: the bias factor for this node
        """
    
        w = self.wins
        n = self.plays
        c = self.explore_factor
        N = self.parent.plays
    
        return w/n + c * np.sqrt( np.log(N) / n )
        
        
    def create_branches(self):
        """
        Creates new MC tree nodes corresponding to all the
        play options available at the game state represented by this node.
        The new nodes are created as :class:`MC_node` objects and
        saved in the list self.branches of this node.
        """
        
        # ask game for a list of all (reasonble) play actions (moves) in this situation
        the_plays = game.all_reasonable_plays(self.state, game.next_player(self.player))

        # loop over all the possible moves and create a new MC_node
        # representing the state of the game after that move
        for play in the_plays:
        
            # We must not change the current state of the game, so we
            # ask for a copy. Depending on the game, this may be just
            # a full copy of the current state, but it may also contain
            # some differences such as added randomization or obscured
            # information. (E.g. if there is a deck of cards, the AI
            # should not be allowed to know the order of the cards.)
            new_state = game.copy_game(self.state)

            # advance the game by taking the move
            game.make_move(play, game.next_player(self.player), new_state)
                        
            # create the new MC tree node
            new_branch = MC_node(player=game.next_player(self.player), 
                                 state=new_state, 
                                 layer=self.layer+1, 
                                 play=play)

            # record that the new node branches from this node
            new_branch.parent = self
            
            # add the new branch to the list of branches in this node
            self.branches.append(new_branch)
            
            
    def explore_branches(self):
        """
        Travels through the MC tree, starting from this node, looking for a leaf node.
        
        At each node, the algorithm checks if the current node is a leaf node
        using :meth:`MC_node.is_leaf`. If it is, the current node is returned.
        
        If the currently examined node is not a leaf, the algorithm calculates 
        bias factors for all the nodes following the current node and then picks 
        the node with the highest factor usin :meth:`MC_node.pick_played_branch`. 
        This node then becomes the current node and the process is repeated.
        
        Eventually, the process will end up at a leaf node
        At that point, the exploration stops and the found leaf is returned.
        
        Returns: 
            MC_node: the found leaf node
        """

        # starting point is this node
        active_node = self

        # if this is a leaf, there's nothing to do
        if active_node.is_leaf():
            return active_node
        
        # search as long as we haven't found a leaf
        while not active_node.is_leaf():

            # move down the tree
            active_node = active_node.pick_played_branch()
            
        return active_node


    def simulate(self):
        """
        Simulates the outcome of the game starting from the current node.
        
        The algorithm alternates between players making each player
        randomly choose one of their available play options until the game ends.
        
        Returns: 
            int: the player number of the winner or 0 if the game ended in a draw
        """
        
        # if the game has been decided, there's nothing to do
        if self.winner > 0:
            return self.winner
            
        # if the game is finished without a winner, it's a draw
        if not game.can_continue(self.state):
            return 0
    
        # the game is not finished, so we will simulate the rest of the game
        # by randomly choosing moves for all players until the game is finished
        finished = False
        
        # in order to not tamper with the game state, make a copy of it
        simulated_state = game.copy_game(self.state)
        
        player = game.next_player(self.player)
        winner = 0
        
        while not finished:
            
            # choose a random play
            the_plays = game.all_reasonable_plays(simulated_state, player)
            random_play = random.choice(the_plays)
            
            # advance the game
            game.make_move(random_play, player, simulated_state)
            
            # check if the game has now ended and return the winner if it has
            winner = game.check_winner(simulated_state)
            if winner > 0:
                return winner
            elif not game.can_continue(simulated_state):
                return 0
            else: # if the game is not finished, it's the next player's turn
                player = game.next_player(player)
                

    def backpropagate(self, winner):
        """
        Recursively records play statistics according to who won the simulated game.
        
        The number of plays from this node is increased by one regardless of
        the outcome of the game.
        
        If the player whose turn it is at this node won the game, the win
        count is also increased by one. If the game ended in a draw,
        the win count is increased by one half.
        
        Once the stats have been recorded for one node,
        the algorithm moves on to the parent node. The process is
        repeated until the root node (the node with no parent) is reached.
        
        Args:
            winner (int): number of the winning player or 0 if the result was a draw
        """
        
        # always record one playout
        self.plays += 1
        
        # Check if the player of this node won.
        # NOTE: self.player is the number of the player
        # whose action *lead* to this node. 
        # Therefore if this node leads to victory, it's
        # a good choice for the player whose decision made it happen
        # only if that player won.
        
        if winner == 0:
            # for a draw, record 0.5 wins
            self.wins += 0.5
            
        elif self.player == winner:
            # for a win, record 1.0 wins
            self.wins += 1
            
        # move down the tree to the previous node unless we are at the root
        if self.parent != None:
            self.parent.backpropagate(winner)
        
    

class AI:
    """
    A Monte Carlo tree search AI.
    
    The AI chooses its actions using a MC tree built out of
    MC_node objects.
    
    Args:
        player (int): the player number for this AI
        thinking_time (int): the number of iterations the AI is allowed when building the MC tree
    """
    

    def __init__(self, player, thinking_time = 5000):
        self.player = player       
        self.iterations = thinking_time


    def pick_move(self, state):
        """
        Given the state of the game, the AI chooses its next move.
        
        Args:
            state: the state of the game in the format defined by the game module
        
        Returns: 
            object: the chosen play action (move) in the format defined by the game module
        """

        # Create the starting node for the MC search tree.
        # Note that the player recorded at a node is the player whose move
        # *lead* to the node. For the current state of the game, that is the
        # previous player, not the AI player itself.
        root = MC_node(game.previous_player(self.player), state, layer=0, play=None)

        # empty line between the game and the AI progress bar
        print()
        
        # explore the MC tree as many times as allowed
        for i in range(self.iterations):
                        
            # progress bar
            self.print_progress(i)
                        
            # explore the existing tree until a leaf node is found
            leaf = root.explore_branches()
            
            # If the game has not finished at the leaf, pick one of the
            # possible moves for which no playouts have been recorded
            # and simulate the rest of the game afer that.
            # Once the game has ended, backpropagate the results to update
            # the play statistics of the leaf node and all the nodes leading
            # to it.
            if leaf.winner == 0 and len(leaf.branches) > 0:
                unplayed = leaf.pick_unplayed_branch()            
                result = unplayed.simulate()
                unplayed.backpropagate(result)

            else: # the game has already finished
                result = leaf.winner
                leaf.backpropagate(result)
             
        # after a set number of iterations, the algorithm must make a decision
                
        # complete progress bar
        self.print_progress(self.iterations)
            
        # The possible choices are the play actions corresponding to the
        # nodes in the first layer of the MC tree 
        # (i.e., the nodes directly following the root).
        # Choose the best one out of these options and return it.
        return root.pick_best_play()


    def print_progress(self, step):
        """
        Draws a progress bar.
        
        Args:
            step (int): the number of iterations taken so far
        """

        total = self.iterations
        message = "["
        total_bar_length = 60
        percentage = int(step / total * 100)
        bar_fill = int(step / total * total_bar_length)
        for i in range(total_bar_length):
            if i < bar_fill:
                message += "|"
            else:
                message += " "
    
        message += "] "+str(percentage)+" %"
        if step < total:
            print(message, end="\r")     
        else:
            print(message) 



def play_game(player_types, thinking_time = 5000):
    """
    Plays a full game.
    
    Args:
        player_types (list): a list of players that may be either HUMAN or AI
        thinking_time (int): the number of iterations allowed for the AI
    
    Returns 
        int: the player number of the winner or 0 for a draw
    """

    # Record the number of players in game
    game.n_players = len(player_types)

    # Start the game
    game_state = game.new_game()
    finished = False
        
    # create the AIs
    ais = [None]*game.n_players
    for p in range(game.n_players):
        if player_types[p] == CPU:
            ais[p] = AI(p+1, thinking_time)
        
    # visualize the starting state
    game.draw(game_state)
    
    # identify the starting player
    player_index = 0
    player_number = 1
    
    # let the game run
    while not finished:

        # ask for a move from either a human or the AI
        if player_types[player_index] == HUMAN:
            play = game.ask_for_move(player_number, game_state)
            
        elif player_types[player_index] == CPU:
            play = ais[player_index].pick_move(game_state)
            
        # play the chosen move
        game.make_move(play, player_number, game_state)
        
        # visualize the current state
        game.draw(game_state)

        # check if someone has won            
        win = game.check_winner(game_state)
            
        if win > 0: # we have a winner
            finished = True
        elif not game.can_continue(game_state): # it's a draw
            finished = True
        else: # still going, next player's turn
            player_number = game.next_player(player_number)
            player_index = player_number - 1     
        
               
    # visualize the end state
    game.draw(game_state) 
    
    # return the winner's number
    return win
        

if __name__ == "__main__":

    # initialize the rng
    random = default_rng()

    # play the game
    winner = play_game([CPU, CPU], 5000)

    # declare the winner
    game.declare_winner(winner)


