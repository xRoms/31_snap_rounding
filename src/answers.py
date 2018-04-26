from sortedcontainers import SortedListWithKey
from structs import Segment, Pixel, SweepLine as SL, Point, normalize, get_pixel, add_pixel_to_seg
from common import hot, current, line, segments
from collections import deque



def bcheck(point):
	pixel = Pixel(point)
	return (pixel.is_on_top(point) and pixel.get_top_neighbour() in current) or (pixel.is_on_bottom(point) and pixel.get_bottom_neighbour() in current)


def segment_endpoint_answer(point, segment):
	pixel = get_pixel(point)
	if pixel not in current:
		# до этого в пикселе не было обнаружено критических точек
		heat_answer(pixel)
	# если событие -- начало отрезка
	if point == segment.start:
		add_pixel_to_seg(pixel, segment)
		pixel.segs.append(segment)
		# игнорируем ту часть отрезка, что лежит внутри пикселя
		for intersection in segment.intersections(pixel):
			if intersection is not None:
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break;
	# если событие -- конец отрезка
	else:
		line.remove(line.status, segment)
		line.remove(line.yasegments, segment)

def segseg_intersection_answer(point):
	pair = line.intersections_status[0]
	anotherpair = None
	if (len(line.intersections_segments) > 0):
		anotherpair = line.intersections_segments[0]
	line.sort_intersection(line.status, line.intersections_status)
	if (anotherpair is not None) and ((anotherpair[1] == pair[2] and anotherpair[2] == pair[1]) or (anotherpair[2] == pair[2] and anotherpair[1] == pair[1])):
		line.sort_intersection(line.yasegments, line.intersections_segments)
	new_point = normalize(point)
	pixel = get_pixel(new_point)
	if (pixel.is_on_top(point) and pixel.get_top_neighbour() in current):
		if (pair[1].isbound != 0):
			trueseg = pair[2]
		else:
			trueseg = pair[1]
		for intersection in trueseg.intersections(pixel.get_top_neighbour()):
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				segpix_intersection_answer(point, trueseg, pixel.get_top_neighbour())
		return
	if (pixel.is_on_bottom(point) and pixel.get_bottom_neighbour() in current):
		if (pair[1].isbound != 0):
			trueseg = pair[2]
		else:
			trueseg = pair[1]
		for intersection in trueseg.intersections(pixel.get_bottom_neighbour()):
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				segpix_intersection_answer(point, trueseg, pixel.get_bottom_neighbour())
		return
	if pixel not in current:
		heat_answer(pixel)

def segpix_intersection_answer(point, segment, pixel):
	if pixel.is_on_top(point):
		pixel.upper.append(segment)
	else:
		pixel.lower.append(segment)
	# добавляем отрезок в соответствующий список
	add_pixel_to_seg(pixel, segment)
	pixel.segs.append(segment)
	# игнорируем фрагмент отрезка, находящийся внутри пикселя
	line.remove(line.status, segment)
	lastinter = None
	for intersection in segment.intersections(pixel):
		if ((intersection is not None) and (intersection > point)) and ((lastinter is None) or (intersection > lastinter)):
			lastinter = intersection
	if (lastinter is not None):

		line.push(SL.Event(SL.Event.Type.SEG_REINSERT, lastinter, segment=segment, pixel=pixel))


def segment_reinsertion_answer(point, segment, pixel):
	line.insert(line.yasegments, segment, line.intersections_segments, shouldpush=False, msg="segment")
	line.insert(line.status, segment, line.intersections_status, shouldpush=True, msg="status")
	if pixel.is_on_top(point):
		pixel.upper.append(segment)
	else:
		pixel.lower.append(segment)
	neighbour = pixel.get_neighbour(point)
	if (neighbour is not None) and (neighbour in current):
		segpix_intersection_answer(point, segment, neighbour)

def pixel_end_answer(point, pixel):
	if pixel.is_on_top(point):
		line.remove(line.status, pixel.top())
	else:
		line.remove(line.status, pixel.bottom())
	hot.extend(current)
	current.clear()

def heat_answer(pixel):
	pixel.center = normalize(pixel.center)
	if (pixel in current):
		return
	l = []
	maybegood = []
	mypos = current.bisect(pixel)
	lower = None
	higher = None
	if not (0 == len(current)):
		if (mypos > 0 and mypos < len(current)):
			lower = current[mypos - 1]
			maybegood.extend(lower.upper)
			higher = current[mypos]
			maybegood.extend(higher.lower)
		if (mypos == 0):
			higher = current[mypos]
			maybegood.extend(higher.lower)
		else:
			if (mypos == len(current)):
				lower = current[len(current) - 1]
				maybegood.extend(lower.upper)
	
	current.add(pixel)
	
		
	li = 0
	if (lower is not None):

		li = max(li, line.bsearch(line.yasegments, lower.top()))
		while (li in range(1, len(line.yasegments)) and line.yasegments[li].atX(line.xpos) >= lower.top().atX(line.xpos)):
			li -= 1
	hi = len(line.yasegments)
	if (higher is not None):
		hi = min(hi, line.bsearch(line.yasegments, higher.bottom()))
		while (hi in range(0, len(line.yasegments)) and line.yasegments[hi].atX(line.xpos) <= higher.bottom().atX(line.xpos)):
			hi += 1
	extender = line.yasegments[li:hi]
	maybegood.extend(extender)

	for s in maybegood:
		if (s.isbound != 0):
			continue;
		for intersection in s.intersections(pixel) :
			if (intersection is not None):
				add_pixel_to_seg(pixel, s)
				pixel.segs.append(s)
				if (round(intersection.x / intersection.z, 6) <= round(line.xpos, 6)):
					l.append(s)
					break

	for segment in l:
		for intersection in segment.intersections(pixel):
			if intersection is not None and (intersection.x / intersection.z) > line.xpos:
				line.remove(line.status, segment)
				line.push(SL.Event(SL.Event.Type.SEG_REINSERT, intersection, segment=segment, pixel=pixel))
				break
	for segment in pixel.segs:
		inter = segment.intersects(pixel.bottom())
		if inter is not None:
			pixel.lower.append(segment)
		inter = segment.intersects(pixel.top())
		if inter is not None:
			pixel.upper.append(segment)
	top, bottom = pixel.top(), pixel.bottom()
	#top.start.coord = 
	line.insert(line.status, top, line.intersections_status, shouldpush=True, msg="status")
	line.push(SL.Event(SL.Event.Type.PIX_END, top.end, pixel=pixel))
	line.insert(line.status, bottom, line.intersections_status, shouldpush=True, msg="status")
	line.push(SL.Event(SL.Event.Type.PIX_END, bottom.end, pixel=pixel))
