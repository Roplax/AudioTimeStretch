# AudioTimeStretch
AudioTimeStretch
The goal of this project is:
 1. Lightweight,audio  timestretching algorithms written in C.
 2. Input will be be wav file, output should be wav file.
 3. Example test input is: example_input1.wav, example_input2.wav
 4. First step is to simplify the existing codebase, last attempt at Wsola algorithm failed miserably, it is coded in jay_wsola.c
 5. Developing tests for submodules and documentation of key input/outputs for submodules is fairly important.


Breakingdown the Goal:
Refer to : https://www.audiolabs-erlangen.de/content/05-fau/professor/00-mueller/01-students/2011_DriedgerJonathan_TSM_MasterThesis.pdf
 1. Code: Algorithms 1 (page 20), 2(page 30)  and 3 (page 36)  
 2. Apply them on a simple sinusoidal  signal,  save the data produced to a text file
 3. Apply them to  sample wav files.
 
