This python script uses Ben's Burst detection algorithm on alex data in order to extract an alex figure and a corresponding fret proximity ratio. It is meant mainly to be a framework to modify, more than an end result, and you will notice a lot of commented out code as a result of that. If you have trouble running this, make sure you have all the requirements(listed below) installed, or ask Ben to help install the necessary tools.

To run this on a particular data file (for example '2011-09-05-run_022.timetag) the syntax is 
"BenBurstAlex 2011-09-05-run_022.timetag"

There are two parameters you can give it
1)Start Offset.
 This is the time after the delta channel goes high that the data is considered valid, in order to give the aotf some time to equilibrate. The defualt is 1e-6 seconds

2)Output.
 This is to be used if you want a file to be output.  By default the file will be shown on the screen which makes it more interactive, but if you are running this from a batch file, you may want to use this instead so you can save each file in the directory instead of doing it manually.

Program Flow:
   The code which is run directly is the last code, after "if __name__ == '__main__':".  A rough sketch follows.
   1) The parameters are loaded
   2) Ben's burst detection is run: "benSpans = ren_ben_burst_alg(f)"
   3) The output window is set up. Modify this to change the look of the figure
   4) The file is read in: "read_filtered_alex_data(..."
   5) The alex figure is plotted: "plot_alex_burst(..."
   6) The fret is plotted: "plot_fret_burst(..."

Requirements:
  1) Python
  2) Ben's suite of photon tools (in the repository)
  3) Haskel
  4) Ben's burst detection code (also in the repository)
