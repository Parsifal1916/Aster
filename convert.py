latex = r"""
\frac{G m_2}{r^2} *{n [-2 v_2^4+4 v_2^2(v_1 v_2)-2(v_1 v_2)^2..  +\frac{3}{2} v_1^2(n v_2)^2+\frac{9}{2} v_2^2(n v_2)^2-6(v_1 v_2)(n v_2)^2  -\frac{15}{8}(n v_2)^4+(\frac{G m_1}{r})(-\frac{15}{4} v_1^2+\frac{5}{4} v_2^2-\frac{5}{2} v_1 v_2. .+\frac{39}{2}(n v_1)^2-39(n v_1)(n v_2)+\frac{17}{2}(n v_2)^2)  +(\frac{G m_2}{r})(4 v_2^2-8 v_1 v_2+2(n v_1)^2.  ..-4(n v_1)(n v_2)-6(n v_2)^2)]  +(\underline{v}_1-v_2)[v_1^2(n v_2)+4 v_2^2(n v_1)-5 v_2^2(n v_2).  -4(v_1 v_2)(n v_1)+4(v_1 v_2)(n v_2)-6(n v_1)(n v_2)^2  +\frac{9}{2}(n v_2)^3+(\frac{G m_1}{r})(-\frac{63}{4} n v_1+\frac{55}{4} n v_2)  ..+(\frac{G m_2}{r})(-2 n v_1-2 n v_2)]\}  +\frac{G^3 m_2}{r^4} n[-\frac{57}{4} m_1^2-9 m_2^2-\frac{69}{2} m_1 m_2]
"""

from sympy import symbols
from sympy.parsing.latex import parse_latex
from sympy.codegen.ast import Assignment
from sympy.printing.c import ccode

v, m, g, h = symbols('v m g h')
expr = parse_latex(latex)
assignment = Assignment(symbols('E'), expr)
print(ccode(assignment))