@echo off
echo === DrugFlow ===
if not exist .venv (
    echo Membuat virtual environment...
    python -m venv .venv
)
call .venv\Scripts\activate
echo Memasang pustaka...
pip install -q -r requirements.txt
echo Membuka aplikasi di browser...
streamlit run app.py
pause
