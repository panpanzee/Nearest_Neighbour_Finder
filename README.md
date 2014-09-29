Nearest_Neighbour_Finder
========================

An implementation of K nereast neighbor finder based on kd_tree.

Main prog usage:
  
  1. Check out the code to a local direcotry
  2. To build the main prog, run 
        $ g++ kd_tree.cpp nearest_neighbour.cpp -o NN_finder -O3 -lboost_thread-mt
  3. Read the help page of the main prog by calling
        $ ./NN_finder help
  4. Sample usage:
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
        $ ./NN_test 3 10000 20 10 10 1 2
