from sortedcontainers import SortedListWithKey
from utils import bound, number_of_segments, seed
from random import uniform
from structs import rounded, Point, Pixel, Segment, SweepLine as SL
import numpy as np

segments = []

hot = []

def generate_segs(n, seed = None):
	if (seed is not None):
		np.random.seed(seed)

	rand = lambda : np.random.randint(0, 5 * bound)
	for i in range(n):
		q, w= rand(), rand()
		while q == w:
			w = rand()
		segments.append(Segment(Point(q, rand(), 5, homogeneous=True), Point(w, rand(), 5, homogeneous=True)))

generate_segs(number_of_segments, seed)

for q in segments:
	print(q.start, " ", q.end)

current = SortedListWithKey(key=lambda pix: pix.center.y)

line = SL(segments)
