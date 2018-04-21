import numpy as np
from enum import Enum
from sortedcontainers import SortedList
import heapq
from numpy.linalg import det
from utils import bound
# from cg import Point, turn

from utils import eps, handlers

class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y
		self.coord = [x, y, 1]

	def __add__(self, other):
		return Point(self.x + other.x, self.y + other.y)

	def __sub__(self, other):
		return Point(self.x - other.x, self.y - other.y)

	def __neg__(self):
		return Point(-self.x, -self.y)

	def __eq__(self, other):
		return self.x == other.x and self.y == other.y

	def __lt__(self, other):
		if self.x == other.x:
			return self.y < other.y
		return self.x < other.x

	def __le__(self, other):
		return self < other or self == other

	def __str__(self):
		return '({}, {})'.format(self.x, self.y)


def vol(point, *hyperplane):
	return det(np.array([pt.coord for pt in hyperplane] + [point.coord]))

def turn(point, *hyperplane):
	return np.sign(vol(point, *hyperplane))


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
		return self.start.y < other.start.y

	def __le__(self, other):
		return self < other or self == other

	def atX(self, x):
		intersection = self.intersects(Segment(Point(x, 0), Point(x, bound)))
		if intersection is None: 
			return self.start
		return intersection

	def intersects(self, other):
		"""
		пересечение отрезка с другим отрезком или с границами пикселя

		:param other: экземпляр Segment или экземпляр Pixel
		:return: Point точка пересечения или None, если пересечения нет
		"""
		if isinstance(other, Segment):
			a, b = self.start, self.end
			c, d = other.start, other.end
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
				return Point(prod[0] / prod[2], prod[1] / prod[2])
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


class Pixel:
	"""
	Класс пикселя
	"""
	def __init__(self, point):
		"""
		Конструктор пикселя

		:param point: экземпляр Point, точка внутри пикселя
		"""
		self.center = Point(int(round(point.x)), int(round(point.y)))
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

	def round(x):
		"""
		округляет точку до центра ближайшего пикселя

		:param x: экземпляр Point
		:return: Point точка, являющаяся центром пикселя, которому принадлежит x
		"""
		if (x % eps <= eps / 2): 
			return x - x % eps
		else:
			return x + eps - x % eps

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
		return point.y == self.y + eps / 2

	def is_on_bottom(self, point):
		"""
		проверяет, лежит ли точка на нижней грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		return point.y == self.y - eps / 2

	def is_on_left(self, point):
		"""
		проверяет, лежит ли точка на левой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		return point.x == self.x - eps / 2

	def is_on_right(self, point):
		"""
		проверяет, лежит ли точка на правой грани пикселя

		:param point: экземпляр Point
		:return: bool
		"""
		return point.x == self.x + eps / 2

	def get_neighbour(self, point):
		"""
		возвращает пиксель, смежный по той грани, которой принадлежит точка

		:param point: экземпляр Point
		:return: Pixel смежный пиксель
		"""
		if self.is_on_top(point):
			return Pixel(Point(self.x, self.y + eps))
		if self.is_on_bottom(point):
			return Pixel(Point(self.x, self.y - eps))
		if self.is_on_left(point):
			return Pixel(Point(self.x - eps, self.y))
		if self.is_on_right(point):
			return Pixel(Point(self.x + eps, self.y))
		else:
			raise Exception('given point is not on pixel\'s side')

	@property
	def x(self):
		return self.center.x

	@property
	def y(self):
		return self.center.y

	@property
	def nw(self):
		return Point(self.x - eps / 2, self.y + eps / 2)

	@property
	def ne(self):
		return Point(self.x + eps / 2, self.y + eps / 2)

	@property
	def sw(self):
		return Point(self.x - eps / 2, self.y - eps / 2)

	@property
	def se(self):
		return Point(self.x + eps / 2, self.y - eps / 2)


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
			self.point = point
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
			if self.x == other.x:
				if self.y == other.y:
					return self.etype.value < other.etype.value
				return self.y < other.y	
			return self.x < other.x 

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
		print("-----")
		self.status.append(segment)
		self.status.sort(key=lambda seg: seg.atX(self.xpos))
		i = self.status.index(segment)
		for seg in self.status:
			print(seg)
		if i != 0:
			p = segment.intersects(self.status[i - 1])
			if p is not None:
				self.push(SweepLine.Event(SweepLine.Event.Type.SEG_SEG, p))
		if i != len(self.status) - 1:
			p = segment.intersects(self.status[i + 1])
			if p is not None:
				self.push(SweepLine.Event(SweepLine.Event.Type.SEG_SEG, p))

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
