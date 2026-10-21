"""Figure 1: cohort flow chart. Writes a draw.io file, which is exported with the draw.io CLI.

Counts: TriNetX Criteria Analysis exports taken on the same day as the primary analyses
  (<export root>/Figure1_AnalyzeCriteria/, one xlsx per cohort definition: S, G, F_S, F_G, M_S, M_G); the criteria are
  applied in sequence, so each exclusion count is the drop from the previous step. The final counts equal
  the before-matching cohorts of the primary analyses (checked below). Matched pairs: Main_Full,
  Female_Breast, Female_Endo_noHyst, Male exports. The number of healthcare organizations is the HCO count of
  the Network row of the Criteria Analysis export.
Usage: python3 -P figure1_flowchart.py <Figure1_AnalyzeCriteria dir> <export root> <out.drawio>
"""
import html
import re
import sys
import zipfile
from pathlib import Path

from lib_trinetx_export import before_n, export_dir, matched_n

STEPS = [  # (group tag in the criteria sheet, label of what is removed; codes are in Method S1)
    ("Group 2:", "No type 2 diabetes"),
    ("Group 4A", "Prior SGLT2i or GLP-1 RA use"),
    ("Group 3A", "{other} on the index date"),
    ("Group 9A", "Immunomodulator use"),
    ("Group 8A", "ESKD or dialysis"),
    ("Group 7A", "HIV"),
    ("Group 5A", "Transplant"),
    ("Group 1:", "Fewer than 2 visits"),
]
def read_criteria(path):
    z = zipfile.ZipFile(path)
    x = z.read("xl/sharedStrings.xml").decode()
    ss = [html.unescape("".join(re.findall(r"<t[^>]*>(.*?)</t>", si, flags=re.S))) for si in re.findall(r"<si>(.*?)</si>", x, flags=re.S)]
    sheet = z.read("xl/worksheets/sheet2.xml").decode()
    rows = []
    for r in re.findall(r"<row[^>]*>(.*?)</row>", sheet, flags=re.S):
        cells = {}
        for c in re.finditer(r'<c r="([A-Z]+)\d+"(?: [^>]*?t="(\w+)")?[^>]*?(?:/>|>(.*?)</c>)', r, flags=re.S):
            col, t, inner = c.groups()
            if inner:
                m = re.search(r"<v>(.*?)</v>", inner, flags=re.S)
                if m:
                    cells[col] = ss[int(m.group(1))] if t == "s" else m.group(1)
        if "A" in cells and "B" in cells:
            hcos = int(float(cells["D"])) if "D" in cells else None   # column D = HCOs
            rows.append((cells["A"], int(float(cells["B"])), hcos))
    return rows   # (criterion, patients, healthcare organizations)


def counts(path):
    rows = read_criteria(path)
    net, hcos = next((n, h) for a, n, h in rows if a == "Network")
    first = next(n for a, n, _ in rows if a.startswith("Group 6A"))
    seq = [first]
    for tag, _ in STEPS:
        seq.append(next(n for a, n, _ in rows if a.startswith(tag)))
    return net, hcos, seq


def esc(s):
    return html.escape(s, quote=True).replace("\n", "&lt;br&gt;")


class Diagram:
    def __init__(self):
        self.cells, self.n = [], 1

    def box(self, x, y, w, h, text, bold=False, stroke=1, font=11, align="center"):
        self.n += 1
        style = (f"rounded=0;whiteSpace=wrap;html=1;strokeWidth={stroke};fontFamily=Arial;fontSize={font};"
                 f"align={align};verticalAlign=middle;spacingLeft=4;spacingRight=4;fontStyle={1 if bold else 0};")
        self.cells.append(f'<mxCell id="b{self.n}" value="{esc(text)}" style="{style}" vertex="1" parent="1">'
                          f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')
        return f"b{self.n}", (x, y, w, h)

    def arrow(self, src, dst, exit_x=0.5, entry_x=0.5):
        self.n += 1
        style = (f"edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;endArrow=block;endFill=1;strokeWidth=1;"
                 f"exitX={exit_x};exitY=1;exitDx=0;exitDy=0;entryX={entry_x};entryY=0;entryDx=0;entryDy=0;")
        self.cells.append(f'<mxCell id="e{self.n}" style="{style}" edge="1" parent="1" source="{src}" target="{dst}">'
                          f'<mxGeometry relative="1" as="geometry"/></mxCell>')

    def xml(self):
        return ('<mxfile host="figure1_flowchart.py"><diagram name="Figure 1" id="fig1"><mxGraphModel grid="0" page="0" '
                'arrows="1" math="0" shadow="0"><root><mxCell id="0"/><mxCell id="1" parent="0"/>'
                + "".join(self.cells) + "</root></mxGraphModel></diagram></mxfile>")


def main(crit_dir, data_root, out):
    crit_dir = Path(crit_dir)
    # one workbook per cohort definition, named <prefix>_<S|G|F_S|F_G|M_S|M_G>_in ...xlsx
    f = {k: next(p for p in sorted(crit_dir.glob("*.xlsx")) if re.match(rf"[^_]+_{k}_in ", p.name))
         for k in ("S", "G", "F_S", "F_G", "M_S", "M_G")}
    net, hcos, seqS = counts(f["S"]); netG, _, seqG = counts(f["G"])
    assert net == netG
    fS, fG, mS, mG = (counts(f[k])[2][-1] for k in ("F_S", "F_G", "M_S", "M_G"))
    full, fem, endo, male = (export_dir(data_root, s) for s in ("Main_Full", "Female_Breast", "Female_Endo_noHyst", "Male"))
    assert (seqS[-1], seqG[-1]) == before_n(full), (seqS[-1], seqG[-1], before_n(full))
    assert (fS, fG) == before_n(fem) and (mS, mG) == before_n(male), ((fS, fG), before_n(fem), (mS, mG), before_n(male))
    mF, mW, mE, mM = matched_n(full)[0], matched_n(fem)[0], matched_n(endo)[0], matched_n(male)[0]
    assert matched_n(full) == (mF, mF)

    # Layout: short labels without codes, bus connectors, thick frames for the main nodes.
    # Width 760 units = 178 mm, so font 11 ≈ 7.3 pt.
    d = Diagram()
    W, colw, gap = 760, 330, 100
    xL, xR = 0, colw + gap
    y = 0
    top, _ = d.box((W - 430) / 2, y, 430, 46, f"TriNetX US Collaborative Network\n(N = {net:,}; {hcos} healthcare organizations)", stroke=2, font=12)
    y += 46 + 34
    cols = {}
    for key, x, seq, drug, other in (("S", xL, seqS, "SGLT2i", "GLP-1 RA"), ("G", xR, seqG, "GLP-1 RA", "SGLT2i")):
        yy = y
        h0, _ = d.box(x, yy, colw, 36, drug, stroke=2, font=12, bold=True)
        d.arrow(top, h0)
        yy += 36 + 26
        b1, _ = d.box(x, yy, colw, 56, f"New users, 2013–2024, without any prior neoplasm\n(n = {seq[0]:,})", font=11)
        d.arrow(h0, b1)
        yy += 56 + 26
        lines = ["Exclusions, applied in sequence:"]
        for (_, lab), a, b in zip(STEPS, seq[:-1], seq[1:]):
            lines.append(f"{lab.format(other=other)} (n = {a - b:,})")
        b2, _ = d.box(x, yy, colw, 136, "\n".join(lines), font=11)
        d.arrow(b1, b2)
        yy += 136 + 26
        b3, _ = d.box(x, yy, colw, 46, f"Eligible cohort\n(n = {seq[-1]:,})", stroke=2, font=12)
        d.arrow(b2, b3)
        cols[key] = b3
    y = yy + 46 + 34
    psm, _ = d.box((W - 430) / 2, y, 430, 40, "1:1 propensity score matching within each cohort", font=12)
    d.arrow(cols["S"], psm); d.arrow(cols["G"], psm)
    y += 40 + 34
    sw, sg = 178, 16
    texts = [
        f"Full cohort\n{mF:,} matched pairs\nOutcomes followed:\nkidney, bladder, colorectal,\nlung, liver, pancreatic,\nesophageal, thyroid and\noral cancers",
        f"Women\n{mW:,} matched pairs\nOutcome followed:\nbreast cancer",
        f"Women without\nprior hysterectomy\n{mE:,} matched pairs\nOutcome followed:\nendometrial cancer",
        f"Men\n{mM:,} matched pairs\nOutcome followed:\nprostate cancer",
    ]
    subs = []
    for i, t in enumerate(texts):
        b, _ = d.box(i * (sw + sg), y, sw, 120, t, font=11)
        subs.append(b)
        d.arrow(psm, b)
    y += 120 + 34
    end, _ = d.box((W - 430) / 2, y, 430, 40, "Follow-up from 1 to 10 years after the index date", stroke=2, font=12)
    for b in subs:
        d.arrow(b, end)
    Path(out).write_text(d.xml(), encoding="utf-8")
    print("network", net, "eligible", seqS[-1], seqG[-1], "women", fS, fG, "men", mS, mG, "matched", mF, mW, mE, mM)


if __name__ == "__main__":
    main(*sys.argv[1:4])
