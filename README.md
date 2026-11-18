# 🧬 DrugFlow

Aplikasi web sederhana untuk alur penemuan obat berbasis fisika komputasi.
Masukkan ID protein dan ID senyawa, aplikasi mengurus sisanya sampai skor docking.

---

## Database yang dipakai (semuanya gratis, tanpa API key)

| Database | Dipakai untuk | Alamat |
|---|---|---|
| **RCSB Protein Data Bank** | Struktur 3D protein, berkas `.pdb` | `files.rcsb.org/download/{ID}.pdb` |
| **RCSB Data API** | Judul, resolusi, metode eksperimen | `data.rcsb.org/rest/v1/core/entry/{ID}` |
| **PubChem (NCBI) PUG REST** | Struktur senyawa dalam bentuk SMILES | `pubchem.ncbi.nlm.nih.gov/rest/pug` |

Tidak ada database yang perlu kamu unduh atau pasang sendiri — semuanya diambil
langsung dari internet saat aplikasi dijalankan. Jadi komputer harus online.

---

## Cara menjalankan

### Windows — cara tercepat
Klik dua kali **`jalankan.bat`**.
Berkas itu otomatis membuat virtual environment, memasang pustaka, dan membuka
aplikasi di browser.

### Windows / macOS / Linux — lewat terminal

```bash
cd drugflow_app

python -m venv .venv

# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

pip install -r requirements.txt

streamlit run app.py
```

Browser akan terbuka sendiri di `http://localhost:8501`.
Untuk menghentikannya, tekan `Ctrl + C` di terminal.

---

## Memasang AutoDock Vina (untuk tahap docking)

Tanpa Vina, aplikasi tetap berjalan: pengambilan data, pembersihan protein,
perhitungan Lipinski, grid box, dan tampilan 3D semuanya berfungsi. Hanya
tahap docking yang dilewati.

**Windows — paling mudah:** unduh `vina_1.2.x_win.exe` dari halaman *Releases*
AutoDock Vina di GitHub, ganti namanya jadi `vina.exe`, lalu taruh di folder
aplikasi ini. Aplikasi akan menemukannya sendiri.

**macOS / Linux:**
```bash
pip install vina
```

**Supaya hasilnya lebih akurat** (sangat disarankan), pasang juga Open Babel
dan PDBFixer lewat conda:
```bash
conda install -c conda-forge openbabel pdbfixer openmm
```
Tanpa Open Babel, aplikasi memakai pengubah `.pdbqt` bawaan yang lebih kasar —
tetap jalan, tapi penentuan tipe atomnya tidak setepat Open Babel. Status setiap
perkakas ditampilkan di panel kiri aplikasi.

---

## Contoh yang sudah terisi

| Kolom | Nilai | Keterangan |
|---|---|---|
| PDB ID | `6LU7` | Main protease SARS-CoV-2, punya ligan asli N3 |
| Rantai | `A` | |
| CID | `5280343` | Kuersetin |

Senyawa lain untuk dicoba: aspirin `2244`, kurkumin `969516`,
remdesivir `121304016`, nirmatrelvir `155903259`.

---

## Susunan berkas

```
drugflow_app/
├── app.py                 tampilan web Streamlit
├── drugflow/
│   ├── sources.py         tahap 1–2  ambil data PDB & PubChem
│   ├── prep.py            tahap 3    bersihkan protein, siapkan senyawa, Lipinski
│   ├── gridbox.py         tahap 4a   hitung grid box dari ligan asli
│   ├── docking.py         tahap 4b   jalankan AutoDock Vina
│   └── viz.py             tahap 6    tampilan 3D dan grafik skor
├── uji_cepat.py           uji mandiri, jalan tanpa internet
├── requirements.txt
├── jalankan.bat           peluncur Windows
└── jalankan.sh            peluncur macOS/Linux
```

Satu modul untuk satu tahap, sesuai rancangan modular di Borang 2 nomor 6.

---

## Kalau ada yang error

**`ModuleNotFoundError: No module named 'drugflow'`**
Kamu menjalankan dari folder yang salah. Pastikan `cd` dulu ke dalam
`drugflow_app`, lalu `streamlit run app.py`.

**`PDB ID tidak ditemukan`**
Kode PDB harus tepat 4 karakter. Cek di rcsb.org.

**Docking gagal dengan pesan dari Vina**
Paling sering karena grid box meleset dari sisi aktif. Protein tanpa ligan asli
otomatis memakai titik tengah protein, dan itu memang kurang tepat — pilih
struktur PDB yang sudah punya ligan asli.

**Aplikasi lambat saat docking**
Turunkan nilai *exhaustiveness* di panel kiri. Nilai 8 sudah cukup untuk
keperluan demonstrasi.

---

## Catatan kejujuran untuk laporan

Simulasi molecular dynamics (tahap 5 di borang) **tidak** disertakan dalam
aplikasi ini. Simulasi 20–100 ns butuh GPU dan waktu berjam-jam, tidak cocok
dijalankan dari aplikasi web yang ditunggu di depan layar. Kalau tahap itu perlu
ditampilkan, jalankan terpisah dengan OpenMM lalu sebutkan durasinya apa adanya.
