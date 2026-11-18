"""
TAHAP 4a — Menentukan grid box (kotak pencarian docking).

Cara paling andal: pakai posisi ligan asli yang sudah ada di berkas PDB.
Kalau protein tidak punya ligan asli, dipakai titik tengah seluruh protein
sebagai cadangan — hasilnya kurang tajam, tapi program tetap jalan.
"""

# Molekul hetero yang BUKAN ligan obat, jangan dijadikan acuan kotak.
BUKAN_LIGAN = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "DMS",
    "ACT", "MPD", "TRS", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE",
}


def _koordinat(baris: str):
    return (float(baris[30:38]), float(baris[38:46]), float(baris[46:54]))


def cari_ligan_asli(teks_pdb: str):
    """
    Cari ligan asli terbesar di berkas PDB.
    Kembaliannya: (nama_residu, daftar_koordinat) atau (None, []).
    """
    kumpulan = {}
    for baris in teks_pdb.splitlines():
        if not baris.startswith("HETATM"):
            continue
        resname = baris[17:20].strip()
        if resname in BUKAN_LIGAN:
            continue
        kunci = (resname, baris[21:22], baris[22:26].strip())
        try:
            kumpulan.setdefault(kunci, []).append(_koordinat(baris))
        except ValueError:
            continue

    if not kumpulan:
        return None, []

    # Ligan obat biasanya yang jumlah atomnya paling banyak.
    kunci, koord = max(kumpulan.items(), key=lambda x: len(x[1]))

    # Sebagian inhibitor (misalnya N3 pada 6LU7) tersimpan sebagai beberapa
    # potongan residu terpisah. Potongan lain yang letaknya berdekatan
    # digabungkan supaya grid box mencakup seluruh ligan, bukan sepotong saja.
    gabungan = list(koord)
    nama = [kunci[0]]
    pusat = _titik_tengah(koord)

    for k2, c2 in kumpulan.items():
        if k2 == kunci:
            continue
        if _jarak(pusat, _titik_tengah(c2)) <= 12.0:
            gabungan.extend(c2)
            if k2[0] not in nama:
                nama.append(k2[0])

    return "+".join(nama), gabungan


def _titik_tengah(koord):
    n = len(koord)
    return (sum(k[0] for k in koord) / n,
            sum(k[1] for k in koord) / n,
            sum(k[2] for k in koord) / n)


def _jarak(a, b):
    return sum((a[i] - b[i]) ** 2 for i in range(3)) ** 0.5


def hitung_kotak(teks_pdb: str, padding: float = 5.0, minimal: float = 20.0) -> dict:
    """
    Hitung titik pusat dan ukuran grid box, dalam satuan Angstrom.

    padding : ruang tambahan di sekeliling ligan asli
    minimal : ukuran sisi terkecil yang diizinkan
    """
    resname, koord = cari_ligan_asli(teks_pdb)
    sumber = f"ligan asli {resname}"

    if not koord:
        # Cadangan: pakai seluruh atom protein.
        koord = [
            _koordinat(b) for b in teks_pdb.splitlines() if b.startswith("ATOM")
        ]
        sumber = "titik tengah protein (ligan asli tidak ditemukan)"
        padding, minimal = 0.0, 25.0

    if not koord:
        raise ValueError("Tidak ada koordinat atom yang bisa dibaca dari berkas PDB.")

    xs, ys, zs = zip(*koord)

    def tengah(v):
        return round((min(v) + max(v)) / 2, 3)

    def sisi(v):
        return round(max(max(v) - min(v) + 2 * padding, minimal), 1)

    hasil = {
        "center": [tengah(xs), tengah(ys), tengah(zs)],
        "size": [sisi(xs), sisi(ys), sisi(zs)],
        "sumber": sumber,
        "jumlah_atom_acuan": len(koord),
    }

    # Grid box yang kelewat besar membuat docking lambat dan tidak fokus.
    if sumber.startswith("titik tengah"):
        hasil["size"] = [30.0, 30.0, 30.0]
    return hasil
