import sys
import os.path
sys.path.append(
    os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

from square_marker_detect import *

def bench():
    cap = cv2.VideoCapture('/Users/patrickfuerst/Documents/Projects/Pupil-Laps/worldvideos/marker.mp4')
    status,img = cap.read()
    markers = []
    while status:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        markers = detect_markers_robust(img,5,markers,true_detect_every_frame=1)
        status,img = cap.read()



if __name__ == '__main__':
    import pstats, cProfile,subprocess,os
    cProfile.runctx("bench()",{},locals(),"world.pstats")
    loc = os.path.abspath(__file__).rsplit('pupil_src', 1)
    gprof2dot_loc = os.path.join(loc[0], 'pupil_src', 'shared_modules','gprof2dot.py')
    subprocess.call("python "+gprof2dot_loc+" -f pstats world.pstats | dot -Tpng -o world_cpu_time.png", shell=True)
    print "created  time graph for  process. Please check out the png next to this file"
    s = pstats.Stats("world.pstats")
    s.strip_dirs().sort_stats("time").print_stats()
