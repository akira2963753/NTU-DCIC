"""Generate the Q3(b) block diagram (fig/q3b.drawio).

Export to PDF with the draw.io desktop app:
    draw.io.exe -x -f pdf --crop -o q3b.pdf q3b.drawio
"""
import sys
from xml.sax.saxutils import escape

FS = 17                                        # base font size (px)
FONT = f"fontFamily=Times New Roman;fontSize={FS};fontColor=#000000;"
SMALL = 'style="font-size:14px"'
INK = "#1a1a1a"
cells = []
_n = [1]


def _id():
    _n[0] += 1
    return f"c{_n[0]}"


def _cell(x, y, w, h, label, style):
    i = _id()
    value = escape(label, {'"': "&quot;"})
    cells.append(f'<mxCell id="{i}" value="{value}" style="{style}" vertex="1" parent="1">'
                 f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')
    return i


def box(cx, cy, w, h, label, extra=""):
    """Block centred at (cx, cy)."""
    st = f"rounded=0;whiteSpace=wrap;html=1;strokeWidth=1.5;strokeColor={INK};fillColor=#ffffff;" + FONT + extra
    return _cell(cx - w / 2, cy - h / 2, w, h, label, st)


def circle(cx, cy, label):
    st = f"ellipse;html=1;strokeWidth=1.5;strokeColor={INK};fillColor=#ffffff;" + FONT + "fontSize=20;"
    return _cell(cx - 14, cy - 14, 28, 28, label, st)


def dot(cx, cy):
    return _cell(cx - 3.5, cy - 3.5, 7, 7, "", f"ellipse;html=1;fillColor={INK};strokeColor={INK};")


def text(x, y, w, h, label, extra=""):
    st = "text;html=1;strokeColor=none;fillColor=none;whiteSpace=wrap;align=center;verticalAlign=middle;" + FONT + extra
    return _cell(x, y, w, h, label, st)


def group(x, y, w, h, label):
    st = ("rounded=1;arcSize=2;html=1;dashed=1;dashPattern=6 4;strokeWidth=1.3;strokeColor=#4a5a78;"
          "fillColor=#f3f5f9;verticalAlign=top;align=center;spacingTop=3;" + FONT + "fontStyle=1;")
    return _cell(x, y, w, h, label, st)


def wire(src, dst, pts=None, ex=None, en=None, arrow=True):
    i = _id()
    st = (f"edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;strokeWidth=1.5;strokeColor={INK};"
          + ("endArrow=block;endFill=1;endSize=7;" if arrow else "endArrow=none;"))
    if ex:
        st += f"exitX={ex[0]};exitY={ex[1]};exitDx=0;exitDy=0;"
    if en:
        st += f"entryX={en[0]};entryY={en[1]};entryDx=0;entryDy=0;"
    geo = '<mxGeometry relative="1" as="geometry">'
    if pts:
        geo += '<Array as="points">' + "".join(f'<mxPoint x="{a}" y="{b}"/>' for a, b in pts) + "</Array>"
    geo += "</mxGeometry>"
    cells.append(f'<mxCell id="{i}" value="" style="{st}" edge="1" parent="1" source="{src}" '
                 f'target="{dst}">{geo}</mxCell>')


def sub(main, small):
    return f'{main}<br><font {SMALL}>{small}</font>'


L, R, T, B = (0, 0.5), (1, 0.5), (0.5, 0), (0.5, 1)
X0, X1 = 20, 1100                              # left / right edge of the whole figure

# =================================================================== DCORR
YC, YE = 105, 255                              # numerator row, energy row
group(X0, 20, X1 - X0, 335, "DCORR: Delay Correlator and Window Energy &#183; S1&#8211;S3")

zin = text(26, YC - 24, 64, 48, "z<sub>m</sub><br>16b I/Q")
s1 = box(140, YC, 70, 36, "S1 reg")
wire(zin, s1, ex=R, en=L)
j1 = dot(195, YC)
wire(s1, j1, arrow=False)
j1b = dot(195, 170)
wire(j1, j1b, arrow=False)

zl = box(285, 170, 120, 46, sub("z<sup>&#8722;L</sup>", "128 &#215; 32b"))
cj = box(395, 170, 56, 36, "conj")
mu = circle(465, YC, "&#215;")
wire(j1b, zl, en=L)
wire(zl, cj, ex=R, en=L)
wire(cj, mu, pts=[(465, 170)], ex=R, en=B)
wire(j1, mu, en=L)

j2 = dot(545, YC)
wire(mu, j2, ex=R, arrow=False)
text(478, YC - 30, 66, 24, "p<sub>m</sub> 33b", "fontSize=15;")
zr = box(635, YC, 120, 46, sub("z<sup>&#8722;R</sup>", "32 &#215; 33b"))
wire(j2, zr, en=L)
sb = circle(760, YC, "+")
text(726, YC - 36, 20, 20, "&#8722;", "fontSize=20;fontStyle=1;")
wire(zr, sb, ex=R, en=L)
ad = circle(760, 170, "+")
wire(j2, ad, pts=[(545, 170)], en=L)
wire(ad, sb, ex=T, en=B)
j3 = dot(835, YC)
wire(sb, j3, ex=R, arrow=False)
z1 = box(835, 170, 52, 36, "z<sup>&#8722;1</sup>")
wire(j3, z1, en=T)
wire(z1, ad, ex=L, en=R)

pw = box(395, YE, 70, 36, "|z<sub>m</sub>|<sup>2</sup>")
wire(j1b, pw, pts=[(195, YE)], en=L)
j4 = dot(545, YE)
wire(pw, j4, ex=R, arrow=False)
text(462, YE - 30, 60, 24, "32b", "fontSize=15;")
zr2 = box(635, YE, 120, 46, sub("z<sup>&#8722;R</sup>", "32 &#215; 32b"))
wire(j4, zr2, en=L)
sb2 = circle(760, YE, "+")
text(726, YE - 36, 20, 20, "&#8722;", "fontSize=20;fontStyle=1;")
wire(zr2, sb2, ex=R, en=L)
ad2 = circle(760, 320, "+")
wire(j4, ad2, pts=[(545, 320)], en=L)
wire(ad2, sb2, ex=T, en=B)
j5 = dot(835, YE)
wire(sb2, j5, ex=R, arrow=False)
z12 = box(835, 320, 52, 36, "z<sup>&#8722;1</sup>")
wire(j5, z12, en=T)
wire(z12, ad2, ex=L, en=R)

text(862, YC - 30, 170, 24, "C(m): 38b signed", "fontSize=15;")
text(862, YE - 30, 170, 24, "E(m): 37b unsigned", "fontSize=15;")

# =================================================================== row 2
GY, GH = 400, 290
YA, YB, YM = 465, 625, 545                     # C path, E path, middle row

# ------------------------------------------------------------------ NORM_SQ
group(X0, GY, 480, GH, "NORM_SQ: Shift and Square &#183; S4&#8211;S5")
ab = box(100, YA, 120, 48, "|Re C|, |Im C|")
lod = box(225, YM, 140, 58, sub("leading-one detect", "&#8594; s, 4b"))
shA = box(355, YA, 64, 48, "&#8811; s")
shB = box(355, YB, 64, 48, "&#8811; s")
sqA = box(445, YA, 86, 56, sub("Re<sup>2</sup> + Im<sup>2</sup>", "num 47b"))
sqB = box(445, YB, 86, 56, sub("E<sup>2</sup>", "den 46b"))

jA = dot(225, YA)
wire(ab, jA, ex=R, arrow=False)
wire(jA, shA, en=L)
wire(jA, lod, en=T)
jB = dot(225, YB)
wire(jB, shB, en=L)
wire(jB, lod, en=B)
jS = dot(322, YM)
wire(lod, jS, ex=R, arrow=False)
wire(jS, shA, pts=[(322, 508), (355, 508)], en=B)
wire(jS, shB, pts=[(322, 582), (355, 582)], en=T)
wire(shA, sqA, ex=R, en=L)
wire(shB, sqB, ex=R, en=L)

# C(m): outer right channel -> lower channel -> down into |Re C|, |Im C|
wire(j3, ab, pts=[(1080, YC), (1080, 384), (60, 384)], en=(1 / 6, 0))
# E(m): inner right channel -> upper channel -> left of NORM_SQ -> lower path
wire(j5, jB, pts=[(1060, YE), (1060, 372), (8, 372), (8, YB)], arrow=False)

# ------------------------------------------------------------------ DIVIDER
DT, DH = 435, 220
group(525, GY, 190, GH, "DIVIDER &#183; S6&#8211;S21")
dv = box(620, DT + DH / 2, 166, DH,
         "16-stage pipelined<br>restoring divider<br><br>"
         "&#934; = &#8970;num &#183; 2<sup>14</sup>/den&#8971;<br><br>"
         f"<font {SMALL}>saturate if num &#8805; 4 den<br>&#934; = 0 if den = 0</font><br><br>output: 16b U2.14")
wire(sqA, dv, ex=R, en=(0, (YA - DT) / DH))
wire(sqB, dv, ex=R, en=(0, (YB - DT) / DH))

# ------------------------------------------------------------------ THRESH
group(740, GY, 130, GH, "THRESH &#183; S22")
eta = text(755, 440, 100, 48, sub("&#951; = 8520", "16b"))
th = box(805, YM, 96, 62, sub("&#934; &gt; &#951; ?", "&#8594; above"))
wire(eta, th, ex=B, en=T)
wire(dv, th, ex=R, en=L)

# ------------------------------------------------------------------ PEAK_SEARCH
group(895, GY, X1 - 895, GH, "PEAK_SEARCH")
ctl = box(990, 470, 160, 58, f"control: run, hold<br><font {SMALL}>emit at cnt = 32 or last</font>")
cmp_ = box(990, YM, 160, 44, "&#934; &gt; pk_phi ?")
reg = box(950, 625, 74, 50, "pk_phi<br>pk_idx")
cnt = box(1030, 625, 70, 50, "cnt 6b")
wire(th, cmp_, ex=R, en=L)
wire(cmp_, reg, ex=(0.25, 1), en=T)
wire(cmp_, cnt, ex=(0.75, 1), en=T)
text(954, 573, 46, 18, "update", "align=left;fontSize=13;")
text(1034, 573, 40, 18, "reset", "align=left;fontSize=13;")
wire(cnt, ctl, pts=[(1080, 625), (1080, 484)], ex=R, en=(1, 0.74))

out = text(900, 708, 200, 30, "det_valid, det_idx, det_phi", "fontSize=15;")
wire(reg, out, ex=B, en=(0.25, 0))
wire(ctl, out, pts=[(1092, 455), (1092, 700)], ex=(1, 0.24), en=(0.96, 0))

# =================================================================== note
text(X0, 748, X1 - X0, 26,
     "Every stage also carries the 10-bit sample index m, valid and last.  z<sub>m</sub> to &#934;: 22 cycles;  "
     "peak decision: W + 1 = 33 cycles after the peak.", "fontStyle=2;fontSize=15;")

xml = ('<mxfile host="drawio"><diagram id="q3b" name="Q3(b)"><mxGraphModel dx="1100" dy="800" grid="0" '
       'gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" '
       'pageWidth="1120" pageHeight="790" math="0" shadow="0"><root><mxCell id="0"/><mxCell id="1" parent="0"/>'
       + "".join(cells) + "</root></mxGraphModel></diagram></mxfile>")
with open(sys.argv[1] if len(sys.argv) > 1 else "q3b.drawio", "w", encoding="utf-8") as f:
    f.write(xml)
print("cells:", len(cells))
