import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import copy


def print_progress(step, total):
    """
    Prints a progress bar.
    
    Args:
        step (int): progress counter
        total (int): counter at completion
    """

    message = "progress ["
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


class NeuralNetwork:
    """
    An artificial neural network.
    
    On creation, all weights are given random values between 0 and 1.
    
    Args:
        inputsize (int): number of input nodes
        layersizes (array): Lists the number of nodes in each hidden layer.
            For example, [5, 6] will result in two hidden layers where the
            first one has 5 and the second 6 nodes.
        outputsize (int): number of output nodes
        learning_rate (float): the learning rate, should be between 0 and 1
    """

    def __init__(self, inputsize, layersizes, outputsize, learning_rate=0.5):
        self.inputsize = inputsize
        self.outputsize = outputsize
        self.training_time = 0
        self.learning_rate = learning_rate
        
        self.weights = []
        
        self.weights.append( random.random( size=(inputsize,layersizes[0]) ).T )
        
        for i in range(1,len(layersizes)):
            self.weights.append( random.random( size=(layersizes[i-1], layersizes[i]) ).T )
            
        self.weights.append( random.random( size=(layersizes[-1], outputsize) ).T )

        self.signal_in = [0]*(len(self.weights)+1)
        
        
    def nudge(self,amount):
        """
        Randomly change weights.
        
        If learning gets stuck in a local optimum, one can try
        this to escape.
        
        Args:
            amount (float): the maximum change allowed in each weight.
        """
        for w in self.weights:
            w += 2*amount*random.random( size = w.shape )-amount
            
        
    def activation(self, signal):
        """
        The activation function.
        
        Neural networks can use different types of activation functions.
        This function implements the sigmoid function
        
        .. math::
            \\varphi(x) = \\frac{1}{1 + e^{-x}}.   
        
        Args:
            signal (array): input :math:`x` either as a float or an array of floats
            
        Returns:
            float or array: output :math:`\\varphi(x)`
        """
    
        return 1.0 / ( 1.0 + np.exp(-signal) )
        
        
    def activation_derivative(self, signal_out):
        """
        Derivative of the :meth:`activation` function.
        
        The derivative of the sigmoid, :math:`\\varphi(x) = \\frac{1}{1 + e^{-x}}` is
        
        .. math::
            \\varphi'(x) = -\\frac{e^{-x}}{(1 + e^{-x})^2}.
        
        However, since :math:`1 - \\varphi(x) = -\\frac{e^{-x}}{1 + e^{-x}}`,
        the derivative can also be written nicely in terms of the output value
        :math:`\\varphi` instead of the input :math:`x` as
        
        .. math::
            \\varphi'(x) = \\varphi(x) [1 - \\varphi(x)].
            
        Args:
            signal_out (array): sigmoid value :math:`\\varphi(x)` either as a float or an array of floats
        
        Returns:
            float or array: sigmoid derivative :math:`\\varphi'(x)`
        """

        return signal_out * (1.0 - signal_out)
        
    
    def feedforward(self, input):
        """
        Sends the signal through the network.
        
        In other words, produces output for the given input.
        
        The neurons in the input layer receive the given input :math:`x` as their
        activation signal. If the signal a neuron receives is strong enough, 
        the neuron activates and sends a new signal :math:`y` to the neurons in 
        the next layer.
        
        To simulate the strength of the connection between neurons, the 
        signal a neuron sends is multiplied by a coupling factor called a weight, :math:`w`.
        (If a weight is 0, there is no connection.)
        Neurons in layers other than the input layer receive signals
        from several neurons, and so for them the total activation signal
        is the sum of the weighted signals. If this sum of signals is strong enough,
        this neuron activates and sends a signal forward, etc.
        
        In this manner, the signal proceeds through the network.
        The signal sent by the final layer is the final output of the whole network.
        
        To be more precise, let us write the activation signal for 
        neuron :math:`i` in layer :math:`n` as :math:`x_i^n`.
        Activation of this neuron is represented by the :meth:`activation`
        function, which changes rapidly from 0 to 1 as the signal goes from
        negative to positive values. (So if :math:`x_i^n > 0`, the neuron activates.)
        The activation output of this neuron is therefore
        
        .. math ::
        
            y_i^n = \\varphi(x_i^n).
            
        The signal that is sent to neuron :math:`j` in layer :math:`n+1` is
        this output multiplied by the weight that connects the two
        neurons,
        
        .. math ::
        
            w_{i,j}^{n \\to n+1} y_i^n.
        
        The total activation signal for neuron :math:`j` is the sum of
        all signals it receives from layer :math:`n`,
        
        .. math ::
        
            x_j^{n+1} = \sum_{i} w_{i,j}^{n \\to n+1} y_i^n.
        
        This summation can be written efficiently with matrices.
        Define
        
            * input vector to layer :math:`n` as :math:`X^{n} = [x_0^{n}, x_1^{n}, \\ldots]^T`  
            * output vector from layer :math:`n` as :math:`Y^n = [y_0^n, y_1^n, \\ldots]^T`
            * weight matrix :math:`W^{n \\to n+1}` with elements :math:`w_{i,j}^{n \\to n+1}`.

        Then neutron activation in layer :math:`n` is calculated with
        
        .. math ::
            
            Y^n = \\varphi(X^n)
            
        and the activation signals for layer :math:`n+1` are obtained with
        
        .. math ::
        
            X^{n+1} = W^{n \\to n+1} Y^{n}.
        
        Args:
            input (array): input (for the input layer)
        
        Returns:
            array: output (from the output layer)
        """
        
        input = np.array([input]).T
        
        layer = 0
        signal = input
        self.signal_in[layer] = signal
        
        for w in self.weights:
        
            signal = self.activation( w @ signal )
            layer += 1
            self.signal_in[layer] = signal
        
        output = copy.copy(signal)
        output.shape = self.outputsize
        return output
        
        
    def backpropagate(self, target, output):
        """
        Compares the output to the target and adjusts weights to drive
        the output towards the target value.
        
        When this function is called, the weights of the network are 
        slightly adjusted so that the output of the network will
        resemble the given target somewhat better. When this function
        is repeatedly called with different learning samples, 
        the network gradually adjusts to reproduce the wanted results.
         
        Mathematically, backpropagation is a one-step gradient search for
        optimal weights :math:`w_{i,j}^{n \\to n+1}`.
        If :math:`E` is the error between the network output
        and the known result, the function calculates the derivatives
        :math:`\\frac{\\partial E}{\\partial w_{i,j}^{n \\to n+1}}` 
        and adjusts the weights by
        
        .. math ::
        
            \\Delta w_{i,j}^{n \\to n+1} = -\\eta \\frac{\\partial E}{\\partial w_{i,j}^{n \\to n+1}}.

        This means the weights are all adjusted in the direction
        that makes the error diminish. 
        
        Here :math:`\\eta` is the learning rate which controls how much
        the weights are adjusted. Typically, it should be between 0 and 1.
        
        Args:
            target (array): the known correct answer to some input
            output (array): the answer the network gives for the same input
        """
    
        # We use the following symbols:
        # x = input for a neuron
        # y = output from a neuron
        # w = network connection weight
        # t = output target
        # E = output error
        
        # Let's use sum of squares error E = sum (y-t)^2.
        # This has the derivative -dE/dy = 2 (t-y).
        # We save this as the vector "error".
        error = np.array( [target - output] ).T
        
        # number of weight matrices
        # this is same as the number of layers - 1
        n_weights = len(self.weights)
                
        # loop over all layers
        for i in range(n_weights):
        
            # the current weight layer
            # Note: we start from the output layer and go
            # towards the input layer.
            layer = n_weights - i - 1
            
            # For the output layer, the delta vector is defined as
            # delta = dE/dy dy/dx,
            # where dE/dy is stored in "error"
            # and dy/dx is given by the activation derivative.
            #
            # For other layers, the delta vector is
            # delta(n) = sum[ dE/dy(n+1) dy(n+1)/dx(n+1) dx(n+1)/dy(n) ] dy(n)/dx(n). 
            # Here the sum is over all neurons in the layer n+1.
            # But we have
            # dx(n+1)/dy(n) = w(n->n+1) and
            # dE/dy(n+1) dy(n+1)/dx(n+1) = delta(n+1),
            # and so
            # delta(n) = sum[ delta(n+1) w(n->n+1) ] dy(n)/dx(n)
            # The result of sum[ delta(n+1) w(n->n+1) ] should 
            # already be saved in "error" and
            # dy(n)/dx(n) is given by the activation derivative.
            #
            # Note that we calculate the derivative using the *output* at layer n,
            # y(n), which is the same as the *input* for layer n+1.
            #
            delta = error * self.activation_derivative(self.signal_in[layer+1])
            
            # Since we need sum[ delta(n+1) w(n->n+1) ]
            # to calculate the adjustments for the next layer n,
            # we pre-emptively save this sum in "error".
            #
            error = self.weights[layer].T @ delta
            
            # The weights are adjusted by
            #   -eta dE/dw(n-1->n) 
            # = -eta sum[ dE/dy(n+1) dy(n+1)/dx(n+1) dx(n+1)/dy(n) ] dy(n)/dx(n) dx(n)/dw(n-1->n)
            # = -eta delta(n) dx(n)/dw(n-1->n).
            #
            # But since x(n) = sum[ w(n-1->n) y(n-1) ], we have
            # dx(n)/dw(n-1->n) = y(n-1).
            # The correct adjustment is therefore
            #   -eta delta(n) y(n-1).
            #            
            self.weights[layer] += self.learning_rate * delta @ self.signal_in[layer].T
            
            
    def train(self, input, target):
        """Trains the network.
        
        The network takes the given input, calculates an output
        and compares the result to the given target output using
        :meth:`NeuralNetwork.backpropagate`.
        
        Calling this function several times with a large group of
        input - target pairs will make the network learn to reproduce
        the given target results.
        
        .. note ::
            This function is incomplete!
        
        Args:
            input (array): input to evaluate
            target (array): the correct answer
        """
    
        self.training_time += 1

        # todo
        
        
    def save_weights(self, filename="weights.txt"):
        """
        Print the current network weights in a file.
        
        Args:
            filename (str): name of the file to write
        """
    
        f = open(filename, "w")
        
        f.write(str(len(self.weights))+"\n")
        
        for w in self.weights:
        
            ni, nj = w.shape
            
            f.write(str(ni)+","+str(nj)+"\n")
            for i in range(ni):
                line = ""
                for j in range(nj):
                    line += str(w[i,j])+","
                    
                f.write(line[:-1]+"\n")
            
        f.close()

    def read_weights(self, filename="weights.txt"):
        """
        Reads network weights from a file.
        
        Args:
            filename (str): name of the file to read
        """
    
        f = open(filename)
        nw = int(f.readline())

        self.weights = [0]*nw
        
        for w in self.weights:
            shape = f.readline()
            parts = shape.split(",")
            ni = int(parts[0])
            nj = int(parts[1])
            
            w = np.array([ni,nj])
            for i in range(ni):
                line = f.readline()
                parts = line.split(",")
                for j in range(nj):
                    w[i,j] = float(parts[j])
                    
        f.close()
            
            
    def visualize(self):
        """
        Draws a visual representation of the network.
        
        Each node is represented as a circle and each layer as a row of circles.
        Input nodes are on the left, and output nodes are on the right.
        
        Weights between nodes are represented by arrows.
        Positive weights are red while negative ones are blue.
        The thicker the arrow, the larger the absolute value of the weight.
        """
        
        n_nodes = [0]
        node_max = 0
        w_max = 0

        i = 0
        for w in self.weights:
            n2, n1 = w.shape
            n_nodes[-1] = n1
            n_nodes.append(n2)
            if n1 > node_max:
                node_max = n1
            if n2 > node_max:
                node_max = n2
                
            for i in range(n2):
                for j in range(n1):
                    if w[i,j] > w_max:
                        w_max = w[i,j]
        
        n_layers = len(n_nodes)

        plt.clf()
        ax = plt.axes()
        ax.set_aspect('equal') 
        plt.xlim([0, 2*n_layers])
        plt.ylim([0, node_max+1])
        
        centers = np.zeros([n_layers, node_max, 2])
        
        for n in range(n_layers):
            m = n_nodes[n]
            for i in range(m):
                x = 2*n+1
                y = 0.5*(node_max-m+2)+i
                centers[ n, i, : ] = [x, y]
        
        
        for n in range(n_layers):
            for i in range(n_nodes[n]):
                x = centers[ n, i, 0 ]
                y = centers[ n, i, 1 ]
                
                if n < n_layers-1:
                    w = self.weights[n]
                    for j in range(n_nodes[n+1]):
                        dx = centers[ n+1, j, 0 ] - centers[ n, i, 0 ]
                        dy = centers[ n+1, j, 1 ] - centers[ n, i, 1 ]                        
                        
                        weight = w[j,i]
                        if weight > 0:
                            c = 'r'
                        else:
                            c = 'b'
                        rel = np.abs( weight/w_max )
                        a = min( np.abs(weight), 1 )
                        t = (0.8*rel+0.2)*0.05
                        
                        plt.arrow( x,y,dx,dy, color=c, 
                                    width = t,
                                    alpha = a,
                                    length_includes_head=True,
                                    head_length = 0.3 )
                        
                plt.gca().add_artist( plt.Circle( [x,y], 0.1, color='k' ) )
                    
                
        plt.show()

def pick_class(output):
    """
    Chooses the most likely class from the given output.
    
    Neural networks are often used to classify data.
    For instance, if we want to sort data instances in three classes,
    we can use a network with three outputs. Each output corresponds to
    a class and the output value (between 0 and 1) represents how
    likely the instance is from that class, according to the network. 
    If the output is [1,0,0], the instance is certainly from the 1st class.
    If the output is [0.1, 0.7, 0.1], the instance is likely from the 2nd class.
    
    This function looks at an output vector and gives the index of the
    class with the highest value.
    For [1,0,0], the function return 0.
    For [0.1, 0.7, 0.1], the function return 1.
    If there is a tie, the function returns the smallest of the tied indices.
    
    Args:
        output (array): neural network output
    
    Returns:
        int: index of the most likely class
    """

    pick = -1
    max = -1

    for i in range(len(output)):
        if output[i] > max:
            max = output[i]
            pick = i
            
    return pick
    

def check_performance(nn, inputs, targets, plot=False, printout=False, classify=False):
    """
    Checks how well a neural network has been trained.
    
    The inputs are given to the neural network and the results
    are compared to the target results (correct answers).
    
    The function may print or plot the results if required.
    It always returns the sum of squares error
    
    .. math::
    
        \\frac{1}{N} \\sum_{i=1}^{N} \\sum_{j=1}^{M} (y_{i,j} - t_{i,j} )^2
    
    where
     
    * :math:`N` is the amount of test data (number of inputs and targets),
    * :math:`M` is the length of the output vector,
    * :math:`y_{i,j}` is the jth component of the ith output and
    * :math:`t_{i,j}` is the jth component of the ith target
    
    Args:
        nn (NeuralNetwork): the network to evaluate
        inputs (list): inputs to test
        targets (list): target outputs to compare to
        plot (bool): If True, the results are visualized.
        printout (bool): If True, the results are printed on screen.
        classify (bool): If True, the network is used for classifying results using :meth:`pick_class`.

    Returns:
        float: the sum of squares error
    """

    outputs = []
    error_sq = 0.0

    j = 0
    for input in inputs:
        output = nn.feedforward(input)
        outputs.append(output)
        error_sq += np.sum( np.square(output[:]-targets[j,:]) )        
        
        if printout:
            if classify:
                predicted_class = pick_class(output)
                true_class = pick_class(targets[j,:])
                if predicted_class == true_class:
                    ok = "correct"
                else:
                    ok = "incorrect !!! "
                print("NN class / true class : ", 
                    predicted_class,
                    true_class,
                    ok )
            else:
                print("output , target , error: ", np.round( output, 2 ), 
                    np.round( targets[j,:], 2 ), np.round( output[:]-targets[j,:], 2) )
        
        j += 1
        
    outputs = np.array(outputs)
    error_sq /= targets.size
    
    if plot:      
        for axis in range(len(outputs[0,:])):
            plt.plot(targets[:, axis], targets[:, axis])
            plt.plot(targets[:, axis], outputs[:, axis], 'o')
            plt.xlabel("target") 
            plt.ylabel("output")
            plt.show() 
        
            if len(inputs[0,:]) == 1:
                plt.plot(inputs[:,0], outputs[:,axis], 'o', label='output')
                plt.plot(inputs[:,0], targets[:,axis], 'o', label='target')
                plt.xlabel("input")
                plt.ylabel("result")
                plt.legend()
                plt.show()

    return error_sq
    

def main(input_size, output_size, layers=[5], traincycles=5000, 
         trainfile="trainingdata.csv", testfile="testdata.csv",
         classify=False):
    """
    The main program.
    
    Creates a network, trains is using training data, 
    and tests the performance against separate test data.
    
    Args:
        input_size (int): number of input neurons
        output_size (int): number of output neurons
        layers (list): number of neurons in each hidden layer
        traincycles (int): how many times the training data is fed to the network
        trainfile (str): name of the file containing the training data
        testfile (str): name of the file containing the test data
        classify (bool): If True, the network is used for classifying results using :meth:`pick_class`.
    """

    # Read an shuffle training data.
    # Shuffling is done so that if the data is ordered, you
    # don't first train using only one type of data and then using only another type.
    # This could lead to bias towards the last type you use.
    trainingdata = np.genfromtxt(trainfile, delimiter=",")
    random.shuffle(trainingdata)

    # split the training data to inputs and target outputs
    inputs = trainingdata[:, 0:input_size]
    targets = trainingdata[:, input_size:input_size+output_size]
    
    # Create the ANN
    nn = NeuralNetwork(inputsize=input_size, layersizes=layers, outputsize=output_size)
    
    total_training_time = traincycles*len(trainingdata)
    errors = []
    lowest_error = check_performance(nn, inputs, targets)
    best_weights = nn.weights

    # start training
    for i in range(traincycles):
    
        # start with a fairly large learning rate but make it smaller as you progress
        nn.learning_rate = (1 - 0.9*i/traincycles)*0.5
        
        # for each cycle, have the ANN compare its output once 
        # to each datapoint in the training set
        for j in range(len(inputs)):
            input = inputs[j]
            target = targets[j]
            nn.train(input, target)
        
            print_progress(nn.training_time, total_training_time)
            
        # record how the squared error converges
        error_sq = check_performance(nn, inputs, targets)
        errors.append(error_sq)
        
        # if the current weights are the best yet, save them
        if error_sq < lowest_error:
            lowest_error = error_sq
            best_weights = nn.weights
        
        # The training may get stuck at a local minimum.
        # This will change the weights a little so that
        # the algorithm might find a better solution.
        if i%100 == 0:
            nn.nudge(0.5)
            
    # save the best set of weigths found during the training
    nn.weights = best_weights
    nn.save_weights()
    
    # draw the ANN
    nn.visualize()

    print("plotting error as function of training time")
    plt.plot([0]*len(errors))
    plt.plot(errors)
    plt.xlabel("training cycle")
    plt.ylabel("error $| Y - T |^2$")
    plt.show()

    # check how well the ANN handles the training data
    print("plotting performace against training data")    
    check_performance(nn, inputs, targets, plot=True, printout=True, classify=classify)
    print("")
    
    # read test data and check how well the ANN predicts it
    testdata = np.genfromtxt(testfile, delimiter=",")
    inputs = testdata[:, 0:input_size]
    targets = testdata[:, input_size:input_size+output_size]
    
    print("plotting performace against test data")  
    check_performance(nn, inputs, targets, plot=True, printout=True, classify=classify)
    
        
if __name__ == "__main__":
    
    random = default_rng()
    
    # These will affect performance.
    # You can try changing them.
    # DO NOT change input and output sizes.
    hidden_layers = [5]
    training_time = 2000
    
    # train the ANN to recognize flowers
    main(input_size=4, output_size=3, layers=hidden_layers, 
        traincycles=training_time,
        trainfile="iris-trainingdata.csv", testfile="iris-testdata.csv", classify=True)

    # alternatively, you can try to teach XOR or sin functions to the ANN
    #main(3,1, trainfile="xor-trainingdata.csv", testfile="xor-testdata.csv")
    #main(1,1, layers=[5,5], trainfile="sin-trainingdata.csv", testfile="sin-testdata.csv")
    

        