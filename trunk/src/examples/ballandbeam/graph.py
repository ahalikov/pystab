# coding=utf-8

__author__="Artur"
__date__ ="$01.02.2010 14:35:33$"

from matplotlib.figure import Figure
from matplotlib.figure import SubplotParams
from numpy.oldnumeric.functions import arange

fig = Figure(figsize=(4,3), dpi=85,facecolor='white',edgecolor='lightblue',
    linewidth = 4.0, frameon = True,
    subplotpars=SubplotParams(left=0.1, bottom=0.1, right=0.9,top=0.9,wspace=0.1,hspace=0.1)
)

myplot = fig.add_subplot(2,2,3)
x1 = arange(0.0, 5.0, 0.1)
x2 = arange(0.0, 5.0, 0.02)

def f(t):
    s1 = sin(2*pi*t)
    e1 = exp(-t)
    return multiply(s1,e1)

line2 = self.subplot1.plot(x2, f(x2), color='blue')
line1 = self.subplot1.plot(x1, f(x1), 'ro')

#scrolled window

scrolledwindow1 = gtk.ScrolledWindow()
scrolledwindow1.show ()
vbox1.pack_start(self.scrolledwindow1, True, True, 0)
scrolledwindow1.set_border_width (8)
#
canvas = FigureCanvas(self.figure1) # «упаковать» диграмму внутрь gtk.DrawingArea
canvas.set_size_request(700, 500) # минимальнй размер области рисования
scrolledwindow1.add_with_viewport(self.canvas)
