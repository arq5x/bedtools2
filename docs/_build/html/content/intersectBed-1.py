a = """chr1 10      20\nchr1        30      40"""
b = """chr1 15      20"""

title = "bedtools intersect -a A.bed -b B.bed"
from matplotlib.pyplot import show
from pyplots.plotter import plot_a_b_tool
plot_a_b_tool(a, b, 'intersect', title, 'A.bed', 'B.bed')
show()