"""
Uji cepat logika inti DrugFlow. Tidak butuh internet, tidak butuh RDKit.
Jalankan dengan:  python uji_cepat.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drugflow import prep, gridbox

# PDB tiruan: 3 atom protein + 2 air + ligan 4 atom
pdb = "\n".join([
"ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 20.00           N",
"ATOM      2  CA  MET A   1      11.000  10.500  10.200  1.00 20.00           C",
"ATOM      3  C   MET A   1      12.000  11.000  10.400  1.00 20.00           C",
"ATOM      4  CA  GLY B   2      30.000  30.000  30.000  1.00 20.00           C",
"HETATM    5  O   HOH A 101      50.000  50.000  50.000  1.00 30.00           O",
"HETATM    6  O   HOH A 102      51.000  51.000  51.000  1.00 30.00           O",
"HETATM    7  C1  N3  A 201      20.000  20.000  20.000  1.00 25.00           C",
"HETATM    8  C2  N3  A 201      22.000  21.000  20.500  1.00 25.00           C",
"HETATM    9  N1  N3  A 201      21.000  22.000  21.000  1.00 25.00           N",
"HETATM   10  O1  N3  A 201      20.500  20.500  23.000  1.00 25.00           O",
"HETATM   11 ZN    ZN A 301      60.000  60.000  60.000  1.00 30.00          ZN",
"END"])

b, c = prep.bersihkan_protein(pdb, "A")
print("bersih rantai A:", c)
assert c["atom_protein"]==3 and c["air_dibuang"]==2, c
b2, c2 = prep.bersihkan_protein(pdb, None)
print("bersih semua rantai:", c2); assert c2["atom_protein"]==4

nama, koord = gridbox.cari_ligan_asli(pdb)
print("ligan asli:", nama, len(koord)); assert nama=="N3" and len(koord)==4

k = gridbox.hitung_kotak(pdb)
print("grid box:", k)
assert k["center"]==[21.0,21.0,21.5], k["center"]
assert all(s>=20.0 for s in k["size"])

# tanpa ligan asli -> cadangan
tanpa = "\n".join(l for l in pdb.splitlines() if " N3  " not in l)
k2 = gridbox.hitung_kotak(tanpa)
print("cadangan:", k2["sumber"], k2["size"]); assert k2["size"]==[30.0,30.0,30.0]

# protein kosong harus memberi error yang jelas
try:
    prep.bersihkan_protein(pdb, "Z"); print("GAGAL: harusnya error")
except ValueError as e: print("error benar:", str(e)[:45])
print("\nSEMUA UJI LULUS")
