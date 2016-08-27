import sys

dir_out = sys.argv[1]
runID_i = sys.argv[2]

ff = open(dir_out+"/done_mm.txt", "a")
ff.write('%s\n' % (runID_i)) 
ff.close()
