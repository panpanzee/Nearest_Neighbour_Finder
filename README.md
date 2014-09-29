Nearest_Neighbour_Finder
========================

An implementation of K nereast neighbor finder based on kd_tree.
The design ideas and implementatin details can be found in each of the header and source files.

Main prog usage:
  
  1. Check out the code to a local direcotry
  2. To build the main prog, run 
        $ g++ kd_tree.cpp nearest_neighbour.cpp -o NN_finder -O3 -lboost_thread-mt
  3. Read the help page of the main prog by calling
        $ ./NN_finder help
  4. Sample usage: (Noting that for large data set (i.e more than tens of million of people), the initial lauching time could be as much as tens of minutes. Please be patient as it's only to make the query quick. You can also select the fast lauch option. Please refer to ./NN_finder help or read the heading comments in the source file nearest_neighbour.cpp) 

        $ ./NN_finder 10
        Launching Nearest Neighbours Finder ...
        Launching finished.
        Please enter your age, latitude, and longtitude
        25
        51.534622
        -0.052503
        Your info: Anonymous
        Age: 25, Latitude: 51.5346, Longtitude: -0.052503.
        
        Your nearest neighbours are,
        Finlay Haynes
        Age: 25, Latitude: 51.5345, Longtitude: -0.0532591.
        
        Isabella Clark
        Age: 25, Latitude: 51.5351, Longtitude: -0.0516156.
        
        Toby Shah
        Age: 25, Latitude: 51.5345, Longtitude: -0.0510079.
        
        Daniel Wong
        Age: 25, Latitude: 51.5352, Longtitude: -0.0506018.
        
        Spencer Godfrey
        Age: 25, Latitude: 51.5359, Longtitude: -0.0513751.
        
        Alex Moore
        Age: 25, Latitude: 51.5361, Longtitude: -0.0531971.
        
        Tom Murphy
        Age: 25, Latitude: 51.5364, Longtitude: -0.0521158.
        
        Ellis George
        Age: 25, Latitude: 51.5331, Longtitude: -0.0502041.
        
        Isabella Moore
        Age: 25, Latitude: 51.5347, Longtitude: -0.0490735.
        
        Muhammad Kemp
        Age: 25, Latitude: 51.5355, Longtitude: -0.0560965.
        
        
        Entering Q for quit, others for continue.
        q

Test usage after checking out the code:
  1. To build the test prog, run 
     $ g++ kd_tree.cpp test.cpp -o NN_test -O3 -lboost_thread-mt
  2. Read the help page of the test prog by calling
        $ ./NN_test
  2. Sample usage:

        $ ./NN_test 3 1000000 20 10 10 1 2
        KDTree construction running time: 1.83658 secondes
        Brute force construction running time: 0.0777631 secondes
        KDTree searching running time: 0.17951 secondes
        Brute force searching running time: 0.182744 secondes
        Test point 1 passed
        
        KDTree searching running time: 0.185978 secondes
        Brute force searching running time: 0.173287 secondes
        Test point 2 passed
        
        KDTree searching running time: 0.180951 secondes
        Brute force searching running time: 0.174184 secondes
        Test point 3 passed
        
        KDTree searching running time: 0.177471 secondes
        Brute force searching running time: 0.183044 secondes
        Test point 4 passed
        
        KDTree searching running time: 0.182662 secondes
        Brute force searching running time: 0.18032 secondes
        Test point 5 passed
        
        KDTree searching running time: 0.174203 secondes
        Brute force searching running time: 0.177689 secondes
        Test point 6 passed
        
        KDTree searching running time: 0.190131 secondes
        Brute force searching running time: 0.188606 secondes
        Test point 7 passed
        
        KDTree searching running time: 0.17636 secondes
        Brute force searching running time: 0.176076 secondes
        Test point 8 passed
        
        KDTree searching running time: 0.172422 secondes
        Brute force searching running time: 0.177449 secondes
        Test point 9 passed
        
        KDTree searching running time: 0.192135 secondes
        Brute force searching running time: 0.182918 secondes
        Test point 10 passed
        
        KDTree searching running time: 0.184402 secondes
        Brute force searching running time: 0.181628 secondes
        Test point 11 passed
        
        KDTree searching running time: 0.182871 secondes
        Brute force searching running time: 0.181939 secondes
        Test point 12 passed
        
        KDTree searching running time: 0.183759 secondes
        Brute force searching running time: 0.185659 secondes
        Test point 13 passed
        
        KDTree searching running time: 0.190051 secondes
        Brute force searching running time: 0.173445 secondes
        Test point 14 passed
        
        KDTree searching running time: 0.176339 secondes
        Brute force searching running time: 0.192267 secondes
        Test point 15 passed
        
        KDTree searching running time: 0.182864 secondes
        Brute force searching running time: 0.185369 secondes
        Test point 16 passed
        
        KDTree searching running time: 0.18668 secondes
        Brute force searching running time: 0.185134 secondes
        Test point 17 passed
        
        KDTree searching running time: 0.179102 secondes
        Brute force searching running time: 0.180157 secondes
        Test point 18 passed
        
        KDTree searching running time: 0.19316 secondes
        Brute force searching running time: 0.20231 secondes
        Test point 19 passed
        
        KDTree searching running time: 0.220315 secondes
        Brute force searching running time: 0.190524 secondes
        Test point 20 passed
        
        All test passed.

