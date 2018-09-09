# SimpleFitter
This guy can provide simple sensitivity test for experiments such as Double Chooz, JUNO, T2K and DUNE. Not sure all files working currently due to missing input.

A way to compile: 
g++ -std=c++11 -o simple_Dc simple_DC.cc -lRooFit -lHtml -lMinuit -lRooFitCore `root-config --cflags --glibs`
