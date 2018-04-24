from sortedcontainers import SortedList
from utils import bound
from random import uniform
from structs import rounded, Point, Pixel, Segment, SweepLine as SL, mkpoint
import numpy as np

segments = []

def generate_segs(n, seed = None):
	if (seed is not None):
		np.random.seed(seed)

	rand = lambda : np.random.randint(0, 5 * bound)
	for i in range(n):
		q, w= rand(), rand()
		while q == w:
			w = rand()
		segments.append(Segment(Point(q, rand(), 5, homogeneous=True), Point(w, rand(), 5, homogeneous=True)))

generate_segs(9, 111111)

for q in segments:
	print(q.start, " ", q.end)

hot = []
current = SortedList()
line = SL(segments)
