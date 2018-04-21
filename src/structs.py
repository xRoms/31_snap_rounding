from __future__ import print_function
import sys
import numpy as np
from enum import Enum
from sortedcontainers import SortedList
import heapq
from numpy.linalg import det
from utils import bound
from cg import Point, turn
from cg.utils import gcd

from utils import eps, handlers

def smart_gcd(a, b):
	if a == 0:
		return b
	if b == 0:
		return a
	return gcd(a, b)

def normalize(smth):
	m = smart_gcd(smth.x, smart_gcd(smth.y, smth.z))
	q = 1
	if (smth.z <= 0):
		q = 1
	return Point(int(q * (smth.x // m)), int(q * (smth.y // m)), int(q * (smth.z // m)), homogeneous=True)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def vol(point, *hyperplane):
	return det(np.array([pt.coord for pt in hyperplane] + [point.coord]))

def gethomo(a):
	if (isinstance(a, int)):
		return 1
	c = 1
	while (((c * a)%1 != 0)):
		c *= 10
	return c

def mkpoint(a, b):
	c = max(gethomo(a), gethomo(b))
	return Point(np.array([int(a * c), int(b * c), int(c)]), homogeneous=True)

class Segment:
	"""
	Класс отрезка
	"""
	def __init__(self, a, b):
		"""
		Конструктор отрезка

		:param a, b: концы отрезка, экземпляры Point
		"""
		self.start = min(a, b)
		self.end = max(a, b)

	def __eq__(self, other):
		"""
		проверяет два отрезка на равенство

		:param other: экземпляр Segment
		:return: bool
		"""
		return self.start == other.start and self.end == other.end		

	def __lt__(self, other):
		"""
		для корректной работы status в SweepLine
		"""
		new_self, new_other = self.start.same_level(other.start)
		return new_self[1] < new_other[1]

	def __le__(self, other):
		return self < other or self == other

	def atX(self, x):
		intersection = self.intersects(Segment(mkpoint(x, 0), mkpoint(x, bound)))
		if intersection is None: 
			return self.start
		return intersection

	def intersects(self, other):
		"""
		пересечение отрезка с другим отрезком или с границами пикселя

		:param other: экземпляр Segment или экземпляр Pixel
		:return: Point точка пересечения или None, если пересечения нет
		"""
		if isinstance(other, Segment):  #vrode ok
			a, b = normalize(self.start), normalize(self.end)
			c, d = normalize(other.start), normalize(other.end)
			acd, bcd = turn(a, c, d), turn(b, c, d)
			cab, dab = turn(c, a, b), turn(d, a, b)

			do_intersect = False
			if acd == bcd == cab == dab == 0:
				do_intersect = (a <= c <= b or a <= d <= b or c <= a <= d or c <= b <= d)
			else:
				do_intersect = acd != bcd and cab != dab
			
			if do_intersect:
				cross = lambda a, b: np.cross([a.coord], [b.coord])
				prod = np.cross(cross(a, b), cross(c, d))[0]
				# point = Point(prod, homogeneous=True)
				# return Point(int(point.x / point.z), int(point.y / point.z))
				if (prod[0] == 0 and prod[1] == 0 and prod[2] == 0):
					return None
				return Point(prod, homogeneous = True)
			return None

		elif isinstance(other, Pixel):
			upper = self.intersects(other.top())
			right = self.intersects(other.right())
			lower = self.intersects(other.bottom())
			if upper is not None:
				return upper
			if right is not None:
				return right
			if lower is not None:
				return lower
			return None

		else:
			raise Exception('wrong argument type')

	def __str__(self):
		return '({}, {})'.format(self.start, self.end)

def rounddd(x):
	"""
	округляет точку до центра ближайшего пикселя

	:param x: экземпляр Point
	:return: Point, являющийся центром пикселя, которому принадлежит x
	"""
	zeps = gethomo(eps)
	m = zeps // gcd(zeps, x.z)
	new_eps = eps * x.z * m
	new_x = int(halfround(x.x * m, new_eps))
	new_y = int(halfround(x.y * m, new_eps))
	new_z = int(x.z * m)
	new_m = gcd(new_x, (gcd(new_y, new_z)))
	return Point(new_x // new_m, new_y // new_m, new_z // new_m, homogeneous= True)


def halfround(a, new_eps):
	print("hround", a, " ", new_eps)
	if (a % new_eps <= new_eps / 2): 
		return a - a % new_eps
	else:
		return a + new_eps - a % new_eps 


class Pixel:
	"""
	Класс пикселя
	"""
	

	def __init__(self, point):
		"""
		Конструктор пикселя

		:param point: экземпляр Point, точка внутри пикселя
		"""
		self.center = rounddd(point)
		self.upper = SortedList()
		self.lower = SortedList()
		self.segs = []

	def __eq__(self, other):
		"""
		проверяет два пикселя на равенство

		:param other: экземпляр Pixel
		:return: bool
		"""
		return self.center == other.center

	def __lt__(self, other):
		return self.center < other.center

	

	def top(self):
		"""
		возвращает отрезок, соответствующий верхней грани пикселя

		:return: Segment
		"""
		return Segment(self.nw, self.ne)
	
	def bottom(self):
		"""
		возвращает отрезок, соответствующий нижней грани пикселя

		:return: Segment
		"""
		return Segment(self.sw, self.se)

	def left(self):
		"""
		возвращает отрезок, соответствующий левой грани пикселя

		:return: Segment
		"""
		return Segment(self.nw, self.sw)

	def right(self):
		"""
		возвращает отрезок, соответствующий правой грани пикселя

		:return: Segment
		"""
		return Segment(self.ne, self.se)

	def is_on_top(self, point):
		"""
		проверяет, лежит ли точка на верхней грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		new_self, new_point = self.center.same_level(point)
		yself = new_self[1]
		ypoint = new_point[1]
		cureps = eps / 2 * new_self[2]

		return ypoint == yself + cureps

	def is_on_bottom(self, point):
		"""
		проверяет, лежит ли точка на нижней грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		new_self, new_point = self.center.same_level(point)
		yself = new_self[1]
		ypoint = new_point[1]
		cureps = eps / 2 * new_self[2]

		return ypoint == yself - cureps


	def is_on_left(self, point):
		"""
		проверяет, лежит ли точка на левой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		new_self, new_point = self.center.same_level(point)
		xself = new_self[0]
		xpoint = new_point[0]
		cureps = eps / 2 * new_self[2]

		return xpoint == xself - cureps

	def is_on_right(self, point):
		"""
		проверяет, лежит ли точка на правой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		new_self, new_point = self.center.same_level(point)
		xself = new_self[0]
		xpoint = new_point[0]
		cureps = eps / 2 * new_self[2]

		return xpoint == xself + cureps

	def get_top_neighbour(self, point):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Pixel(Point(int(new_self.x), int(new_self.y + 2 * new_eps), int(new_self.z * (zeps // m)), homogeneous = True))

	def get_bottom_neighbour(self, point):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Pixel(Point(int(new_self.x), int(new_self.y - 2 * new_eps), int(new_self.z * (zeps // m)), homogeneous = True))

	def get_left_neighbour(self, point):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Pixel(Point(int(new_self.x - 2 * new_eps), int(new_self.y), int(new_self.z * (zeps // m)), homogeneous = True))

	def get_right_neighbour(self, point):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Pixel(Point(int(new_self.x + 2 * new_eps), int(new_self.y), int(new_self.z * (zeps // m)), homogeneous = True))

	def get_neighbour(self, point):
		"""
		возвращает пиксель, смежный по той грани, которой принадлежит точка

		:param point: экземпляр Point
		:return: Pixel смежный пиксель
		"""
		eprint("start")
		eprint(self.center)
		eprint(point)
		eprint("end")
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		if self.is_on_top(point):
			return Pixel(Point(int(new_self.x), int(new_self.y + new_eps), int(new_self.z * (zeps // m)), homogeneous = True))
		if self.is_on_bottom(point):
			return Pixel(Point(int(new_self.x), int(new_self.y - new_eps), int(new_self.z * (zeps // m)), homogeneous = True))
		if self.is_on_left(point):
			return Pixel(Point(int(new_self.x - new_eps), int(new_self.y), int(new_self.z * (zeps // m)), homogeneous = True))
		if self.is_on_right(point):
			return Pixel(Point(int(new_self.x + new_eps), int(new_self.y), int(new_self.z * (zeps // m)), homogeneous = True))
		else:
			raise Exception('given point is not on pixel\'s side')

	@property
	def x(self):
		return self.center.x

	@property
	def y(self):
		return self.center.y

	@property
	def z(self):
		return self.center.z

	@property
	def nw(self):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Point(int(new_self.x - new_eps), int(new_self.y + new_eps), int(new_self.z * (zeps // m)), homogeneous = True)

	@property
	def ne(self):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Point(int(new_self.x + new_eps), int(new_self.y + new_eps), int(new_self.z * (zeps // m)), homogeneous = True)

	@property
	def sw(self):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Point(int(new_self.x - new_eps), int(new_self.y - new_eps), int(new_self.z * (zeps // m)), homogeneous = True)

	@property
	def se(self):
		cureps = eps / 2
		zeps = int(gethomo(cureps))
		trueps = int(cureps * zeps)
		m = gcd(self.z, zeps)
		new_eps = trueps * (self.z // m)
		new_self = self.center * (zeps // m)
		return Point(int(new_self.x + new_eps), int(new_self.y - new_eps), int(new_self.z * (zeps // m)), homogeneous = True)


class SweepLine:
	"""
	Класс заметающей прямой
	"""
	class Event:
		"""
		Класс события, обрабатываемого заметающей прямой
		"""
		class Type(Enum):
			"""
			Тип события
			"""
			SEG_END      = 0 # конец отрезка
			SEG_START    = 1 # начало отрезка
			SEG_SEG      = 2 # пересечение двух отрезков
			SEG_PIX      = 3 # пересечение отрезка с границей пикселя
			PIX_END      = 4 # конец пикселя
			SEG_REINSERT = 5 # повторное вхождение отрезка после выхода из пикселя


		def __init__(self, etype, point, segment=None, pixel=None):
			"""
			Конструктор события

			:param etype: SweepLine.Event.Type тип события
			:param point: Point точка, в которой происходит событие
			:param segment: Segment для событий типа SEG_START, SEG_END, SEG_PIX, SEG_REINSERT \
			                -- отрезок, участвующий в них
			:param pixel: Pixel для событий типа SEG_PIX, PIX_END -- пиксель, участвующий в них
			"""
			self.etype = etype
			self.point = normalize(point)
			print("new point ", self.point)
			self.segment = segment
			self.pixel = pixel

		def __eq__(self, other):
			"""
			проверяет два события на равенство

			:param other: экземпляр SweepLine.Event 
			:return: bool
			"""
			return self.point == other.point and self.etype == other.etype \
			and self.segment == other.segment and self.pixel == other.pixel

		def __lt__(self, other):
			"""
			лексикографическое "меньше" для событий

			:param other: экземпляр SweepLine.Event
			:return: bool
			"""
			new_self, new_other = self.point.same_level(other.point)
			if new_self[0] == new_other[0]:
				return self.etype.value < other.etype.value	
			return new_self[0] < new_other[0]

		def __str__(self):
			s = 'type: {}\npoint: {}'.format(self.etype.name, self.point)
			if self.segment is not None:
				s += '\nsegment: {}'.format(self.segment)
			if self.pixel is not None:
				s += '\npixel: {}'.format(self.pixel)
			return s

		def handle(self):
			"""
			вызывает функцию-обработчик для события данного типа
			"""
			handler = handlers[self.etype]
			if self.segment is not None and self.pixel is not None:
				handler(self.point, segment=self.segment, pixel=self.pixel)
			elif self.segment is not None:
				handler(self.point, segment=self.segment)
			elif self.pixel is not None:
				handler(self.point, pixel=self.pixel)
			else:
				handler(self.point)

		@property
		def x(self):
			return self.point.x

		@property
		def y(self):
			return self.point.y
		@property
		def z(self):
			return self.point.z


	def __init__(self, segments):
		"""
		Конструктор заметающей прямой

		:param segments: список отрезков Segment, для работы с которыми строится заметающая прямая
		"""
		self.xpos = 0
		self.status = [] # SortedList(key=lambda seg: seg.atX(self.xpos))
		self.events = []

		# инициализируем список событий началами и концами отрезков
		for s in segments:
			self.push(SweepLine.Event(SweepLine.Event.Type.SEG_START, s.start, segment=s))
			self.push(SweepLine.Event(SweepLine.Event.Type.SEG_END, s.end, segment=s))

	def insert(self, segment):
		"""
		добавляет отрезок в статус, а также вычисляет новые пересечения и добавляет их в очередь events

		:param segment: экземпляр Segment
		:param point: экземпляр Point
		"""
		print("-----", segment)
		self.status.append(segment)
		self.status.sort(key=lambda seg: seg.atX(self.xpos))
		i = self.status.index(segment)
		for seg in self.status:
			print(seg)
		print("--endstatus--")
		for eve in self.events:
			print(eve.etype, " ", eve.point)
		print("--endevents--")
		j = i
		while (j > 0):
			if (self.status[i].atX(self.xpos) == self.status[j - 1].atX(self.xpos)):
				j -= 1
				continue
			p = segment.intersects(self.status[j - 1])
			print("intersection of ", segment, " and ", self.status[j - 1], " is ", p)
			if p is not None and (p.x / p.z >= self.xpos):
				self.push(SweepLine.Event(SweepLine.Event.Type.SEG_SEG, p))
			break
		j = i
		while (j < len(self.status) - 1):
			if (self.status[i].atX(self.xpos) == self.status[j + 1].atX(self.xpos)):
				j += 1
				continue
			p = segment.intersects(self.status[j + 1])
			print("intersection of ", segment, " and ", self.status[j + 1], " is ", p)
			if p is not None and (p.x / p.z >= self.xpos):
				self.push(SweepLine.Event(SweepLine.Event.Type.SEG_SEG, p))
			break

	def remove(self, segment):
		"""
		удаляет отрезок из статуса

		:param segment: экземпляр Segment
		"""
		if segment in self.status:
			self.status.remove(segment)

	def push(self, event):
		"""
		добавляет событие в очередь events

		:param event: экземпляр SweepLine.Event
		"""
		heapq.heappush(self.events, event)

	def peek(self):
		"""
		возвращает самое раннее событие из events, не удаляя его из очереди

		:return: SweepLine.Event
		"""
		return self.events[0]

	def pop(self):
		"""
		возвращает самое раннее событие из events, удаляя его из очереди

		:return: SweepLine.Event
		"""
		return heapq.heappop(self.events)
