"""
TAHAP 4b — Menjalankan molecular docking dengan AutoDock Vina.

Program ini sengaja dibuat lentur, karena cara memasang Vina berbeda-beda
di tiap komputer. Urutan yang dicoba:

  Ligan    : Meeko  ->  Open Babel
  Protein  : Open Babel  ->  pengubah bawaan (mode darurat)
  Docking  : berkas vina.exe / vina  ->  pustaka Python 'vina'

Kalau semuanya gagal, pesan errornya menyebut perintah pemasangan yang perlu
dijalankan, bukan sekadar "gagal".
"""

import os
import re
import shutil
import subprocess

# Tipe atom AutoDock sederhana, dipakai hanya pada mode darurat.
_TIPE_DARURAT = {"C": "C", "N": "NA", "O": "OA", "S": "SA", "P": "P",
                 "F": "F", "CL": "Cl", "BR": "Br", "I": "I"}


# --------------------------------------------------------------------------
# Pemeriksaan perkakas
# --------------------------------------------------------------------------
def cari_vina(path_manual: str = "") -> str:
    """Cari berkas program Vina. Kembaliannya path, atau string kosong."""
    kandidat = [path_manual] if path_manual else []
    kandidat += ["vina", "vina.exe", "./vina", "./vina.exe",
                 os.path.join(os.getcwd(), "vina.exe"),
                 os.path.join(os.getcwd(), "vina")]
    for k in kandidat:
        if not k:
            continue
        if os.path.isfile(k) and os.access(k, os.X_OK):
            return os.path.abspath(k)
        ketemu = shutil.which(k)
        if ketemu:
            return ketemu
    return ""


def periksa_perkakas(path_vina: str = "") -> dict:
    """Laporkan perkakas apa saja yang tersedia di komputer ini."""
    def ada(nama):
        try:
            __import__(nama)
            return True
        except ImportError:
            return False

    return {
        "rdkit": ada("rdkit"),
        "meeko": ada("meeko"),
        "openbabel": bool(shutil.which("obabel")),
        "vina_program": cari_vina(path_vina),
        "vina_python": ada("vina"),
        "pdbfixer": ada("pdbfixer"),
    }


# --------------------------------------------------------------------------
# Konversi ke format .pdbqt
# --------------------------------------------------------------------------
def ligan_ke_pdbqt(mol, path_keluar: str) -> str:
    """Ubah molekul RDKit 3D menjadi berkas .pdbqt. Kembaliannya: cara yang dipakai."""
    # --- cara 1: Meeko (paling tepat, murni Python) ---
    try:
        from meeko import MoleculePreparation
        prep = MoleculePreparation()
        hasil = prep.prepare(mol)

        teks = None
        if isinstance(hasil, list) and hasil:
            try:
                from meeko import PDBQTWriterLegacy
                keluaran = PDBQTWriterLegacy.write_string(hasil[0])
                teks = keluaran[0] if isinstance(keluaran, tuple) else keluaran
            except ImportError:
                teks = hasil[0].write_pdbqt_string()
        else:
            teks = prep.write_pdbqt_string()      # Meeko versi lama

        if teks and "ATOM" in teks:
            with open(path_keluar, "w") as f:
                f.write(teks)
            return "Meeko"
    except ImportError:
        pass
    except Exception:
        pass

    # --- cara 2: Open Babel ---
    if shutil.which("obabel"):
        from rdkit import Chem
        sementara = path_keluar.replace(".pdbqt", ".sdf")
        Chem.MolToMolFile(mol, sementara)
        proses = subprocess.run(
            ["obabel", sementara, "-O", path_keluar, "-p", "7.4"],
            capture_output=True, text=True,
        )
        if os.path.exists(path_keluar) and os.path.getsize(path_keluar) > 0:
            return "Open Babel"
        raise RuntimeError(f"Open Babel gagal: {proses.stderr[:300]}")

    # --- cara 3: tulis sendiri dengan RDKit (ligan kaku) ---
    return _pdbqt_ligan_rdkit(mol, path_keluar)


def _pdbqt_ligan_rdkit(mol, path_keluar: str) -> str:
    """
    Cadangan terakhir: tulis .pdbqt ligan langsung dengan RDKit, tanpa Meeko
    maupun Open Babel.

    Keterbatasan: ligan diperlakukan KAKU (tidak ada torsi yang diputar), jadi
    Vina hanya menggeser dan memutar molekul utuh. Skornya cenderung sedikit
    lebih lemah daripada docking fleksibel. Pakai untuk demonstrasi; untuk
    angka yang dilaporkan, pasang Meeko.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.Mol(mol)
    try:
        AllChem.ComputeGasteigerCharges(mol)
    except Exception:
        pass

    def muatan(atom):
        try:
            q = float(atom.GetDoubleProp("_GasteigerCharge"))
            return 0.0 if (q != q) else q        # q != q menandai NaN
        except Exception:
            return 0.0

    konf = mol.GetConformer()
    baris, nomor = ["ROOT"], 0

    for atom in mol.GetAtoms():
        unsur = atom.GetSymbol()
        tetangga = atom.GetNeighbors()

        # Hidrogen nonpolar dilebur ke atom induknya, sesuai kelaziman AutoDock.
        if unsur == "H":
            if tetangga and tetangga[0].GetSymbol() == "C":
                continue
            tipe = "HD"
        elif unsur == "C":
            tipe = "A" if atom.GetIsAromatic() else "C"
        elif unsur == "N":
            tipe = "NA"
        elif unsur == "O":
            tipe = "OA"
        elif unsur == "S":
            tipe = "SA"
        elif unsur in ("F", "P", "I"):
            tipe = unsur
        elif unsur.upper() == "CL":
            tipe = "Cl"
        elif unsur.upper() == "BR":
            tipe = "Br"
        else:
            continue

        q = muatan(atom)
        # Muatan hidrogen nonpolar ditambahkan ke karbon pengikatnya.
        if unsur == "C":
            for t in tetangga:
                if t.GetSymbol() == "H":
                    q += muatan(t)

        p3 = konf.GetAtomPosition(atom.GetIdx())
        nomor += 1
        nama = (unsur + str(nomor))[:4]
        baris.append(
            f"ATOM  {nomor:5d} {nama:<4s} LIG A   1    "
            f"{p3.x:8.3f}{p3.y:8.3f}{p3.z:8.3f}  1.00  0.00    "
            f"{q:6.3f} {tipe:<2s}"
        )

    if nomor == 0:
        raise RuntimeError("Tidak ada atom ligan yang bisa ditulis ke .pdbqt.")

    baris += ["ENDROOT", "TORSDOF 0"]
    with open(path_keluar, "w") as f:
        f.write("\n".join(baris) + "\n")
    return "RDKit (ligan kaku — pasang Meeko untuk docking fleksibel)"


def protein_ke_pdbqt(path_pdb: str, path_keluar: str) -> str:
    """Ubah protein .pdb menjadi .pdbqt kaku. Kembaliannya: cara yang dipakai."""
    # --- cara 1: Open Babel (-xr = reseptor kaku) ---
    if shutil.which("obabel"):
        subprocess.run(
            ["obabel", path_pdb, "-O", path_keluar, "-xr"],
            capture_output=True, text=True,
        )
        if os.path.exists(path_keluar) and os.path.getsize(path_keluar) > 0:
            return "Open Babel"

    # --- cara 2: pengubah bawaan (mode darurat) ---
    return _pdbqt_darurat(path_pdb, path_keluar)


def _pdbqt_darurat(path_pdb: str, path_keluar: str) -> str:
    """
    Pengubah sederhana bawaan program.

    Catatan penting: fungsi skor Vina memakai TIPE ATOM, bukan muatan parsial,
    jadi berkas seperti ini masih bisa dipakai. Tetapi penentuan tipenya kasar
    (semua N dianggap akseptor), sehingga hasilnya kurang akurat dibanding
    Open Babel. Pakai hanya kalau terpaksa.
    """
    keluaran = []
    for baris in open(path_pdb):
        if not baris.startswith(("ATOM", "HETATM")):
            continue

        unsur = baris[76:78].strip().upper()
        if not unsur:
            unsur = baris[12:16].strip()[0].upper()
        if unsur == "H":
            continue                                   # hidrogen dilewati

        tipe = _TIPE_DARURAT.get(unsur)
        if tipe is None:
            continue
        keluaran.append(f"{baris[:66]:<66s}     0.000 {tipe:<2s}\n")

    if not keluaran:
        raise RuntimeError("Gagal mengubah protein ke .pdbqt: tidak ada atom yang terbaca.")

    with open(path_keluar, "w") as f:
        f.writelines(keluaran)
        f.write("TER\n")
    return "pengubah bawaan (mode darurat, akurasi lebih rendah)"


# --------------------------------------------------------------------------
# Menjalankan docking
# --------------------------------------------------------------------------
def jalankan_docking(reseptor_pdbqt, ligan_pdbqt, keluaran_pdbqt,
                     center, size, exhaustiveness=8, n_poses=5,
                     path_vina="") -> dict:
    """
    Jalankan docking. Kembaliannya:
    {"skor": [...], "mesin": "...", "log": "...", "berkas": "..."}
    """
    program = cari_vina(path_vina)

    # --- cara 1: berkas program Vina ---
    if program:
        perintah = [
            program,
            "--receptor", reseptor_pdbqt,
            "--ligand", ligan_pdbqt,
            "--out", keluaran_pdbqt,
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(size[0]),
            "--size_y", str(size[1]),
            "--size_z", str(size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses),
        ]
        proses = subprocess.run(perintah, capture_output=True, text=True, timeout=1800)
        if not os.path.exists(keluaran_pdbqt):
            raise RuntimeError(
                "Vina gagal dijalankan.\n\n"
                f"Pesan error:\n{(proses.stderr or proses.stdout)[:800]}"
            )
        return {
            "skor": baca_skor(keluaran_pdbqt),
            "mesin": f"program Vina ({os.path.basename(program)})",
            "log": proses.stdout,
            "berkas": keluaran_pdbqt,
        }

    # --- cara 2: pustaka Python 'vina' ---
    try:
        from vina import Vina
    except ImportError:
        raise RuntimeError(
            "AutoDock Vina belum tersedia.\n\n"
            "Pilih salah satu:\n"
            "  1. Unduh vina.exe dari halaman rilis AutoDock Vina di GitHub, "
            "lalu taruh di folder aplikasi ini (cara termudah di Windows).\n"
            "  2. pip install vina   (biasanya lancar di Linux dan macOS)"
        )

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(rigid_pdbqt_filename=reseptor_pdbqt)
    v.set_ligand_from_file(ligan_pdbqt)
    v.compute_vina_maps(center=list(map(float, center)), box_size=list(map(float, size)))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    v.write_poses(keluaran_pdbqt, n_poses=n_poses, overwrite=True)

    energi = v.energies(n_poses=n_poses)
    skor = [round(float(baris[0]), 2) for baris in energi]
    return {
        "skor": skor,
        "mesin": "pustaka Python vina",
        "log": "",
        "berkas": keluaran_pdbqt,
    }


def baca_skor(path_pdbqt: str) -> list:
    """Ambil daftar skor docking dari berkas keluaran Vina."""
    skor = []
    pola = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)")
    for baris in open(path_pdbqt):
        cocok = pola.search(baris)
        if cocok:
            skor.append(float(cocok.group(1)))
    return skor

